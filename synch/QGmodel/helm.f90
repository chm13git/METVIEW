!-----------------------------------------------------------
! The solution of the Helmholtz equation

  subroutine helmholtz (psi,pv,ninitial)

  use Sizes

  IMPLICIT NONE

  integer, parameter :: nwork = 6*nxx+5*nyy+15
  real pv(nxx,nyy,nl), psi(nxx,nyy,nl)

  integer i,j,ierror,ninitial
  integer bctype(4)
  real mode(nxx,nyy,nl),modeboundary(nxx,nyy,2)
  real bda(nyy),bdb(nxx),bdc(nyy),bdd(nxx)
  real pertrb
  real sl,denominator
  real q(nxx,nyy),q1(nxx+2,nyy+2)
  real n1,n2,n3,a1,a2,dn1,dn2,dn3,da1,da2,da3,num1,num2
  real work(nwork),gh(nxx+1,nyy+1)

! Determine modal pv fields and store in mode0
  mode = 0.
  mode(:,:,1) = (H1*pv(:,:,1) + H2*pv(:,:,2))/(H1+H2) 
  do i=1,nxx
    mode(i,:,1) = mode(i,:,1) - f(:)
  enddo
  mode(:,:,2) = pv(:,:,1) - pv(:,:,2)

!  ------------------
! Barotropic solve
!    Boundary conditions  3 is periodic and 1 is prescribed streamfunction
  bctype(1) = 3
  bctype(2) = 1
  bctype(3) = 3
  bctype(4) = 1
  do i=1,nxx
    bdb(i) = 0.
    bdd(i) = 0.
  enddo
  do j = 1,nyy    ! Note these are periodic boundaries so value is not used but needed in Helmholtz solver
    bda(j) = 0.
    bdc(j) = 0.
  enddo  

! Fields needed for Helmholtz solver
  gh = 0.
  do i=1,nxx-1
    do j=1,nyy-1
      gh(i+1,j+1) = 0.5*(mode(i,j,1)+mode(i+1,j+1,1))  ! Note, when iorder = 4 this needs to change!
    enddo
  enddo
  q1= 0
  do i=1,nxx
    do j=1,nyy
      q1(i,j) = mode(i,j,1)
    enddo
  enddo

! Helmholtz solver. Note that a small stretching term is needed for stability
  pertrb = 0.
  call HFFT2A(-1.e-25,nxx,nyy,dx,gh,nxx+1,bctype,bda,bdb,bdc,bdd,4,q1,nxx+2,work,nwork,ierror)
  if(ierror .ne. 0) write(*,*)'ierror is ',ierror

! Copy results to mode(i,j,1), which is now the barotropic streamfunction with zero at northern and southern boundaries
  do i=1,nxx
    do j=1,nyy
      mode(i,j,1) = q1(i,j)
    enddo
  enddo
 
! Analytic solution to Laplace(psi) = 0 and psi(south)=0 and psi(north) = 1
  do i=1,nxx
    do j=1,nyy
      modeboundary(i,j,1) = (j-1.)/real(nyy-1.)
    enddo
  enddo

! Determine magnitude of coefficient to form full barotropic solution that fulfils boundary conditions
  n1 = sum(mode(:,nyy,1)-mode(:,nyy-1,1))
  n2 = sum(modeboundary(:,nyy,1)-modeboundary(:,nyy-1,1))
  if (ninitial == 0) then
    cbt = (H1*psi(1,nyy,1)+H2*psi(1,nyy,2))/(H1+H2)
  else
    cbt = cbt - (n1 + cbt*n2)/sumpsiboundarybt + (sumpsibt+cbt*sumpsiboundarybt)/sumpsiboundarybt
  endif
  sumpsibt = n1
  sumpsiboundarybt = n2

! Full barotropic mode streamfunction
  mode(:,:,1) = mode(:,:,1) + cbt*modeboundary(:,:,1)
  
!  -------------------
! baroclinic solve
! Boundary condition, as for barotropic case
  bctype(1) = 3
  bctype(2) = 1
  bctype(3) = 3
  bctype(4) = 1
  do i=1,nxx
    bdb(i) = 0.
    bdd(i) = 0.
  enddo
  do j = 1,nyy
    bda(j) = 0.
    bdc(j) = 0.
  enddo  

! Fields needed by Helmholtz solver
  gh = 0.
  do i=1,nxx-1
    do j=1,nyy-1
      gh(i+1,j+1) = 0.5*(mode(i,j,2)+mode(i+1,j+1,2))
    enddo
  enddo
  q1= 0
  do i=1,nxx
    do j=1,nyy
      q1(i,j) = mode(i,j,2)
    enddo
  enddo

! Helmholtz solver
  call HFFT2A(-lambda,nxx,nyy,dx,gh,nxx+1,bctype,bda,bdb,bdc,bdd,4,q1,nxx+2,work,nwork,ierror)
  if(ierror .ne. 0) write(*,*)'ierror is ',ierror

! Copy results to mode(i,j,2), which is now the baroclinic streamfunction with zero at northern and southern boundaries
  do i=1,nxx
    do j=1,nyy
      mode(i,j,2) = q1(i,j)
    enddo
  enddo
  
! Analytic solution to Laplace(psi) - lambda psi = 0 with southern boundary zero and norther boundary value of 1
! Analytic solution to Laplace(psi) - lambda psi = 0 with southern boundary -1 and norther boundary value of 1
  sl = sqrt(lambda)*dx
  denominator = exp(sl*(nyy-1))-exp(-sl*(nyy-1))
  num1 = 1. + exp(-sl*(nyy-1))
  num2 = 1. + exp(sl*(nyy-1))
  do i=1,nxx
    do j=1,nyy
      modeboundary(i,j,1) = (exp(sl*(j-1))-exp(-sl*(j-1)))/denominator
      modeboundary(i,j,2) = (num1*exp(sl*(j-1))-num2*exp(-sl*(j-1)))/denominator
    enddo
  enddo
  
! Determine magnitude of coefficient to form full baroclinic solution that fulfils boundary conditions
  n1 = sum(mode(:,nyy,2)-mode(:,nyy-1,2))
  n2 = sum(modeboundary(:,nyy,1)-modeboundary(:,nyy-1,1))
  n3 = sum(modeboundary(:,nyy,2)-modeboundary(:,nyy-1,2))
  a1 = sum(mode(:,:,2))
  a2 = sum(modeboundary(:,:,1))
  !!write(*,*)'a1',a1,a2     !QG
  if (ninitial == 0) then
    cbc1 = - a1/a2
    cbc2 = -(psi(1,1,1)-psi(1,1,2))
  else
    dn1 = n1 - n1old
    dn2 = n2 - n2old
    dn3 = n3 - n3old
    da1 = a1 - a1old
    da2 = a2 - a2old
    cbc1 = - a1/a2
    cbc2 = (a2*n3-a2*dn3)*cbc2 - a2*dn1+a1*dn2+n2*da1-a1*n2*da2/a2
    cbc2 = cbc2 /(a2*n3)
  endif
  n1old = n1
  n2old = n2
  n3old = n3
  a1old = a1
  a2old = a2

! Full baroclinic mode streamfunction
  mode(:,:,2) = mode(:,:,2) + cbc1*modeboundary(:,:,1) + cbc2 * modeboundary(:,:,2)

! Determine layer fields from modal fields
  psi(:,:,1) = mode(:,:,1) + H2*mode(:,:,2)/(H1+H2)
  psi(:,:,2) = psi(:,:,1) - mode(:,:,2)

  end subroutine helmholtz

  
