  subroutine InitiateModel(psi)
!--------------------------------------------------
! Initialize some variables to their required value
!--------------------------------------------------
  Use Sizes

  IMPLICIT NONE

  Interface
    subroutine MakeRandomFilter(InverseLengthscale)
    real, intent(IN) :: InverseLengthscale
    end subroutine MakeRandomFilter
  end Interface

  integer i, j,k,imax,jmax
  real  rho(nl), rho0
  real  pi, ga, ds
  real psi(nxx,nyy,nl)
  real umax,u,v,fnot,beta0
  real a,b

  write (*,*) ' Subroutine InitiateModel '
  
! Open a few output files for summary variables
  open(22, file='checks',status='unknown',form='formatted')

  pi    = 4.*atan(1.)     ! the constant pi

! Coriolis parameter on a beta plane, dependent on meridional (y) coordinate
  beta0 = 2.e-11
  fnot  = 7.28e-5
  do j = 1, nyy
    f(j) = (nyy/2.-j) * dx*beta0  + fnot
  enddo

! hulp variable for Laplace and Helmholtz solvers 1 over grid distance squared
  c1     = 1. / (dx*dx)

! Layer thicknesses in m
  H1   = 5000.
  H2   = 5000.
      
! density of the two layers, in kg/m^3
  rho0   = 1.05
  rho(1) = 1.01
  rho(2) = 1.05

! Parameters for stretching term in PV evolution equation, and for transition fields to modes and back
  ga = (rho(2)-rho(1)) * 9.81 / rho0
  F1 = fnot*fnot/(ga*H1)
  F2 = fnot*fnot/(ga*H2)
  lambda = (F1+F2)
   
! Write some numbers for verification
  write(*, *) 'f0, beta, lambda', fnot, beta0, lambda
  write(22, *)'f0, beta, lambda', fnot, beta0, lambda
       
! Set initial streamfunction field
  do i=1,nxx
    do j=1,nyy
      !psi(i,j,1) = - 5.e7*atan((j-nyy/2)*5.*4./real(nyy))&
      !             - 5.e6*sin(2*pi*(i-1)/real(nxx))*sin((j-1)*pi/real(nyy-1))*sin((j-1)*pi/real(nyy-1))
      psi(i,j,1) = - 25.e6*atan((j-nyy/2)*5.*4./real(nyy))&      !QG
                   - 25.e5*sin(2*pi*(i-1)/real(nxx))*sin((j-1)*pi/real(nyy-1))*sin((j-1)*pi/real(nyy-1))            !QG     
      psi(i,j,2) = 0.3*psi(i,j,1)                                                                                   !QG
    enddo
  enddo
! Ensure all sreamfunction values on northern boundary are the same, and same for southern boundary 
  do i=1,nxx
    do k=1,nl
      psi(i,1,k) = psi(1,1,k)
      psi(i,nyy,k) = psi(1,nyy,k)
    enddo
  enddo
! Ensure barotropic streamfunction is zero at southrn boundary, and baroclinic streamfunction averages to zero.
  a = (H1*psi(1,1,1)+H2*psi(1,1,2))/(H1+H2)
  b = a - sum(psi(:,:,1)-psi(:,:,2))/(nxx*nyy)
  do i=1,nxx
    do j=1,nyy
      psi(i,j,1) = psi(i,j,1)-a
      psi(i,j,2) = psi(i,j,2)-b
    enddo
  enddo 

! Store initial streamfuncrion field
  open(10,file='stream00',form = 'unformatted')
  write(10)psi
  close(10)

! A simple test to ensure the maximum velocity is realistic (between 10 and 100 m/s for atmosphere)
  umax=0.
  do k=1,nl
    do j=2,nyy-1
      do i=2,nxx-1
        v=(psi(i,j,k)-psi(i+1,j,k))/dx
        u=(psi(i,j,k)-psi(i,j+1,k))/dx
        u=sqrt(u*u+v*v)
        if(u.gt.umax)then
          umax=u
          imax=i
          jmax=j
        endif
      enddo
    enddo
  enddo
  write(*,*)'umax = ',umax,imax,jmax
  write(22,*)'umax = ',umax,imax,jmax

! Initialize Shapiro filter used to remove small-sclae noise
!  call shfact(8,sh)

! Initialize random field generator
! amplitudes of random fields
  qo(1) = sum(psi(:,:,1)*psi(:,:,1))/(nxx*nyy) - sum(psi(:,:,1))*sum(psi(:,:,1))/(nxx*nxx*nyy*nyy)
  qo(2) = sum(psi(:,:,2)*psi(:,:,2))/(nxx*nyy) - sum(psi(:,:,2))*sum(psi(:,:,2))/(nxx*nxx*nyy*nyy)
  qo(1) = 0.001*sqrt(qo(1))
  qo(2) = 0.001*sqrt(qo(2))
  alpha = 0.3
  call MakeRandomFilter(1./lengthscale)

  end subroutine InitiateModel

!--------------------------------------------------------------------------
!    time step routine
!--------------------------------------------------------------------------
  subroutine IntegrateModel (psi,nstep)

  Use Sizes

  IMPLICIT NONE
  
  integer nstep,ninitial
  integer i,j,k
  character*5 aaa

  real psi(nxx,nyy,nl)
  real p(nxx,nyy),r(nxx,nyy),q(nxx,nyy),pv00(nxx,nyy,nl),q0(nxx,nyy)
  real pv(nxx,nyy,nl),pv0(nxx,nyy),r1(nxx,nyy),r2(nxx,nyy),r3(nxx,nyy),r4(nxx,nyy)
  real psiRandom(nxx,nyy,nl)

  write(aaa,'(i5.5)')nstep
! calculate potential vorticity at next time level for all layers 
  do k = 1, nl
 
!   relative vorticity, store in q0 for use in diffusion term
    call laplacenew (psi(:,:,k), p)
    q0 = p

!   Add Stretching terms to pv
    if (k.eq.1) then
       r = - F1*(psi(:,:,1)-psi(:,:,2))
    else 
       r =   F2*(psi(:,:,1)-psi(:,:,2))
    endif
    p = p +  r

!   add Coriolis parameter to PV
    do i = 1, nxx
      p(i, :) = p(i, :) + f
    enddo
    if (nstep.eq.0) pv00(:,:,k) = p(:,:)

!   Jacobian j(p,psi), this is the advection of PV by the velocities from the streamfunction field using RK4
    pv0(:,:) = p(:,:)
    call xjac (p, psi(:,:,k), r)
    !pv(:,:,k) = p(:,:) + dt*(r(:,:) + diff * q(:,:))
    r1 = r 
    p = pv0 + 0.5*dt*(r + diff * q)
    call xjac(p,psi(:,:,k),r)
    r2 = r 
    p = pv0 + 0.5*dt*(r + diff * q)
    call xjac(p,psi(:,:,k),r)
    r3 = r
    p = pv0 + dt*(r + diff * q)
    call xjac(p,psi(:,:,k),r)
    r4 = r
    call laplacenew(q0,q)
    pv(:,:,k) = pv0 + dt*(r1+r2+r3+r4)/6.
   
!   Add diffusion term
    pv(:,:,k) = pv(:,:,k) + dt*diff*q

  enddo

! The helmholtz solver to obtain the stream function from new PV fields, include psi(..) in the call for boundary conditions
  if (nstep.eq.0) then
    call helmholtz(psi,pv00,0)
    open(10,file='stream1',form = 'unformatted')
    write(10)psi
    close(10)
  endif

! Store streamfunction fields every time step
  open(10,file='stream'//aaa,form = 'unformatted')
  write(10)psi
  close(10)

! Determine new streamfunction fields from new pv fields using Helmholtz solver
  call helmholtz (psi,pv,1)

! generate random fields and add to streamfunction fields
  call GenerateRandomField(psiRandom)
  psi = psi + psiRandom
! Store random field every time step
  open(10,file='random'//aaa,form = 'unformatted')
  write(10)psiRandom
  close(10)

! Shapiro filter to filter out small-scale noise in streamfunction fields, acts as small-scale 'diffusion term'
!  if (mod(nstep,1) .eq. 0) then
!    do k = 1, nl
!      call shfilt2 (8, psi(:,:,k), r)
!    enddo
!  endif

! Calculate the kinetic energy for verification
  if(mod(nstep,1) .eq. 0) then
    call energy (psi,  nstep)
  endif

  end subroutine IntegrateModel

!--------------------------------------------------------------------------
  subroutine energy (psi, nstep)
!-------------------------------------
!          Calculates kinetic energy
!-------------------------------------

  use Sizes

  IMPLICIT NONE

  integer i,j,k,nstep
  real psi(nxx,nyy,nl), xke(nl)

  xke = 0.
  do j = 2 , nyy-1
    do i = 2, nxx-1
      xke = xke + (psi(i+1,j,:) - psi(i-1,j,:))**2 + (psi(i,j+1,:) - psi(i,j-1,:))**2
    enddo
  enddo
  xke = xke/(nxx*nyy*dx*dx)

  write (22, *) " step ", nstep, " kinetic energy at top layer = ", xke(1)
  !if (xke(1) > 10000.) stop  ' xke(1) > 10000. This suggests an instability or a bug.'

  end subroutine energy
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine laplacenew(p,q)
!---------------------------------------
!         The calculation of q = Laplace(p) using 9 grid points for higher accuracy
!----------------------------------------

  use Sizes
  
  IMPLICIT NONE

  integer i,j
  real p(nxx,nyy),q(nxx,nyy)
  real a1,a2,a3

  a1 = -1./12.
  a2 = 4./3.
  a3 = -5.

! Model interior
  do j=3,nyy-2
    do i=3,nxx-2
      q(i,j) = c1*(a1*(p(i+2,j)+p(i-2,j)+p(i,j+2)+p(i,j-2)) + a2*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1))+a3*p(i,j)) 
    enddo
  enddo

! boundary and rows/columns next to boundaries
  do j=3,nyy-2
    i=2
    q(i,j) = c1*(a1*(p(i+2,j)+p(nxx,j)+p(i,j+2)+p(i,j-2)) + a2*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1))+a3*p(i,j)) 
    i=1
    q(i,j) = c1*(a1*(p(i+2,j)+p(nxx-1,j)+p(i,j+2)+p(i,j-2)) + a2*(p(i+1,j)+p(nxx,j)+p(i,j+1)+p(i,j-1))+a3*p(i,j)) 
    i=nxx-1
    q(i,j) = c1*(a1*(p(1,j)+p(i-2,j)+p(i,j+2)+p(i,j-2)) + a2*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1))+a3*p(i,j)) 
    i=nxx
    q(i,j) = c1*(a1*(p(2,j)+p(i-2,j)+p(i,j+2)+p(i,j-2)) + a2*(p(1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1))+a3*p(i,j)) 
  enddo
  do i = 2, nxx-1
    j=2
    q(i,j)= c1 * (p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1) - 4.*p(i,j))
    j=nyy-1
    q(i,j)= c1 * (p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1) - 4.*p(i,j))
  enddo

  i=1
  j=2
  q(i,j)= c1 * (p(i+1,j) + p(nxx,j) + p(i,j+1) + p(i,j-1) - 4.*p(i,j))
  j=nyy-1
  q(i,j)= c1 * (p(i+1,j) + p(nxx,j) + p(i,j+1) + p(i,j-1) - 4.*p(i,j))

  i=nxx
  j=2
  q(i,j)= c1 * (p(1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1) - 4.*p(i,j))
  j=nyy-1
  q(i,j)= c1 * (p(1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1) - 4.*p(i,j))
  
  do i=1,nxx
    q(i,1)   = 0.
    q(i,nyy) = 0.
    !q(i,1)   = q(i,2)
    !q(i,nyy) = q(i,nyy-1)
  enddo
 
  end subroutine laplacenew
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   subroutine for the calculation of the jacobian of 
!   p and q, result put in r. This is the advection term in the PV equation
!
  subroutine xjac(p,q,r)

  Use Sizes

  IMPLICIT NONE

  integer i,j
  real p(nxx,nyy),q(nxx,nyy),r(nxx,nyy)
  
  r = 0.
  do j=2,nyy-1
    do i=2,nxx-1
      r(i,j) =        (p(i,j-1)+p(i+1,j-1)-p(i+1,j+1)-p(i,j+1))*q(i+1,j)
      r(i,j) = r(i,j)-(p(i-1,j-1)+p(i,j-1)-p(i-1,j+1)-p(i,j+1))*q(i-1,j)
    enddo
  enddo
  do j=2,nyy-1
    do i=2,nxx-1
      r(i,j) = r(i,j) + (p(i+1,j)+p(i+1,j+1)-p(i-1,j)-p(i-1,j+1))*q(i,j+1)
      r(i,j) = r(i,j) - (p(i+1,j-1)+p(i+1,j)-p(i-1,j-1)-p(i-1,j))*q(i,j-1)
    enddo
  enddo
  do j=2,nyy-1
    do i=2,nxx-1
      r(i,j) = r(i,j) + (p(i+1,j)-p(i,j+1))*q(i+1,j+1)
      r(i,j) = r(i,j) - (p(i,j-1)-p(i-1,j))*q(i-1,j-1)
      r(i,j) = r(i,j) + (p(i,j+1)-p(i-1,j))*q(i-1,j+1)
      r(i,j) = c1*(r(i,j) - (p(i+1,j)-p(i,j-1))*q(i+1,j-1))/12.
    enddo
  enddo

  return
  end subroutine xjac
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine shfact(n,sh)
!-------------------------------------------
!         Initialization of Shapiro filter
!--------------------------------------------
       integer n  
       real sh(0:16) 
       real f2n,fj,f2nj,ff
 
       if (n.gt.8) then
          write(*,*)'error in "shfact"'
          write(*,*)'n is to large (>8): n=',n
          stop
       endif

       if ((n.eq.1).or.(n.eq.2).or.(n.eq.4).or.(n.eq.8)) then
       ff=2.0**(2*n)
       f2n=1.
       do j=1,2*n
          f2n=f2n*float(j)
       enddo

       fj=1.
       sh(0)=((-1.0)**(n-1))/ff
       do j=1,n

!        /* calculates $(2n-j)!$ for each $j$ */
         f2nj=1.
         do i=1,2*n-j
            f2nj=f2nj*float(i)
         enddo

!        /* calculates $j!$ */
         fj=fj*float(j)

!        /*  calculates the factors */
         sh(j)=(-1.0)**(n+j-1) * f2n/(ff*fj*f2nj) 
       enddo

       else
          write(*,*)'error in shfact.  n=',n
          stop
       endif

       return
       end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine shfilt2(ish,x,y)
!------------------------------------------------
!       Shapiro filter of order 8
!------------------------------------------------

      Use Sizes

      integer j,i,m,mm
      integer ish 

      real x(nxx,nyy) 
      real y(nxx,nyy)

      if ((2*ish+1.gt.nxx).or.(2*ish+1.gt.nyy)) then 
	write(*,*)'shfilt2:  the domain is to small for ish=',ish
          return
      endif
!
! filtering in the y direction.
!
!     /* setting the boundary points at the boundary $m=1$ and $m=n_y$ */
      do i=1,nxx
        do j=1,nyy
          y(i,j)=0.
        enddo
        y(i,1)=x(i,1)
        y(i,nyy)=x(i,nyy)
      enddo

!     /* all grid points close to the boundary |m=1| where |m<=ish| */
      do m=2,ish
        do j=0,ish-1
          mm=1-(m-ish+j)
          if (mm.ge.0) then
            do i=1,nxx
            y(i,m)=y(i,m)+sh(j)*(2.0*x(i,1)-x(i,1+mm)+x(i,m+ish-j))
            enddo
          else
            do i=1,nxx
              y(i,m)=y(i,m)+sh(j)*(x(i,m+ish-j)+x(i,m-ish+j))
            enddo
          endif
        enddo
        do i=1,nxx
          y(i,m)=(1.0+sh(ish))*x(i,m)+y(i,m)
        enddo
      enddo

!     /* all intermediate points |ish<m<=ny-ish| */
      do m=ish+1,nyy-ish
        do j=0,ish-1
          do i=1,nxx
            y(i,m)=y(i,m)+sh(j)*(x(i,m+ish-j)+x(i,m-ish+j))
          enddo
        enddo
        do i=1,nxx
          y(i,m)=(1.0+sh(ish))*x(i,m)+y(i,m)
        enddo
      enddo
         
!     /* all grid points close to the boundary |m=ny| where |m>ny-ish| */
      do m=nyy-ish+1,nyy-1
        do j=0,ish-1
          mm=(m+ish-j)-nyy
          if (mm.ge.0) then
            do i=1,nxx
          y(i,m)=y(i,m)+sh(j)*(2.0*x(i,nyy)-x(i,nyy-mm)+x(i,m-ish+j))
            enddo
          else
            do i=1,nxx
              y(i,m)=y(i,m)+sh(j)*(x(i,m+ish-j)+x(i,m-ish+j))
            enddo
          endif
        enddo
        do i=1,nxx
          y(i,m)=(1.0+sh(ish))*x(i,m)+y(i,m)
        enddo
      enddo
!
! filtering in the x direction.
!
      do j=1,nyy
        do i=1,nxx
          x(i,j)=0.
        enddo
        x(1,j)=y(1,j)
        x(nxx,j)=y(nxx,j)
      enddo

!     /* all grid points close to the boundary |i=1| where |i<=ish| */
      do i=2,ish
        do j=0,ish-1
          mm=1-(i-ish+j)
          if (mm.ge.0) then
            do jj=1,nyy
          x(i,jj)=x(i,jj)+sh(j)*(2.0*y(1,jj)-y(1+mm,jj)+y(i+ish-j,jj))
            enddo
          else
            do jj=1,nyy
              x(i,jj)=x(i,jj)+sh(j)*(y(i+ish-j,jj)+y(i-ish+j,jj))
            enddo
          endif
        enddo
        do jj=1,nyy
          x(i,jj)=(1.0+sh(ish))*y(i,jj)+x(i,jj)
        enddo
      enddo

!     /* all intermediate points |ish<i<=nxx-ish| */
      do i=ish+1,nxx-ish
        do j=0,ish-1
          do jj=1,nyy
            x(i,jj)=x(i,jj)+sh(j)*(y(i+ish-j,jj)+y(i-ish+j,jj))
          enddo
        enddo
        do jj=1,nyy
           x(i,jj)=(1.0+sh(ish))*y(i,jj)+x(i,jj)
        enddo
      enddo
         
!     /* all grid points close to the boundary |i=nx| where |i>nx-ish| */
      do i=nxx-ish+1,nxx-1
        do j=0,ish-1
          mm=(i+ish-j)-nxx
          if (mm.ge.0) then
            do jj=1,nyy
        x(i,jj)=x(i,jj)+sh(j)*(2.0*y(nxx,jj)-y(nxx-mm,jj)+y(i-ish+j,jj))
            enddo
          else
            do jj=1,nyy
              x(i,jj)=x(i,jj)+sh(j)*(y(i+ish-j,jj)+y(i-ish+j,jj))
            enddo
          endif
        enddo
        do jj=1,nyy
          x(i,jj)=(1.0+sh(ish))*y(i,jj)+x(i,jj)
        enddo
      enddo
      return
      end

