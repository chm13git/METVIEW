! This subroutine runs the particle filter by first initialising the ensemble
! and then, in the main time loop, uses synchronization as proposal density
! It needs as input observations, a true initial state, an H operator, and error covariances
! for obs and model. 

Program SPF

!Use Variables
use Sizes         !QG

IMPLICIT NONE

integer i,j,mcn,t,k,l,iter,mm,ll,numiter,nbegin,tau,tt,nn,nycount,obsloc,ycount,ncolb
integer :: i2,jj,ii,kk,n,zz              !QG
integer time, status     !QG
!real*8 xx(nstate,nens,ntau),xxm(nstate,ntau),hxx,d(ntau*nobs,nens),dx(m,nens),dxx(m)
!real*8 yp(nobs),hx(ntau*nobs,nens),coefnew(nens,nens),xpold(m,nens),xpoldm,xinflation(nstate)
!real*8 sx(tstate),var(tstate),coef(nens,nens),gamma,hnew(nens,nens),varpost,varprior,xhulp2(m)
real :: psi(nxx,nyy,nl), psiRandom(nxx,nyy,nl)                            !QG
real :: psin(nxx,nyy,nl,nens)                      !QG 
real :: xpold(ntotal,nens)                              !QG
!real*8 dist,dist1,dist2,reldist,GC,xxmean,xzero(m),factor,dxfactor
real, allocatable, dimension(:,:)  :: dnew
real, allocatable, dimension(:,:)  :: hxnew
character*5 aa


! Initialize in and output files before and after IEWPF step.
open(20,file='statesin.txt')
open(30,file='statesout.txt')
open(40,file='truth.txt')     !FRP
open(50,file='states.txt')    !FRP

!nens = 10.            !QG
!m = 100.              !QG
!nxx = 10.             !QG 
!nyy = 10.             !QG
!nl = 1.               !QG

! Initialise ensemble by adding random perturbation to true state and add perturbations to the resulting state as ensemble members
!xzero = 0.                             ! This is the zero forcing term for model propagation
psi = 0.
!call randomv(randomx,m)
!xmp(:,1) = xtrue(:,1) + ivar * randomx
do zz = 1,nens
   call InitiateModel(psi)	
   write(*, *) 'Running psi...', zz, psi(1,3,1) 
   call GenerateRandomField(psiRandom) !QG 
   !write(*, *) 'Running psiRandom...', psiRandom(1,:,1)   !QG
   psin(:,:,:,zz) = psi + psiRandom    !QG
   write(*, *) 'Running psin...', zz, psin(1,3,1,1) !1st row/3rd column element for nl=1 (first non-zero random element)(nens=1) !QG
   !write(*, *) 'Running psin...', zz, psin(1,1,2,1) !First element for nl=2 (nens=1)    !QG
enddo   

write(*, *) 'Finished InitiateModel'

do time = time_to_start, time_to_stop
  write(*,*) ' time = ', time

  do zz = 1,nens
     call IntegrateModel (psin(:,:,:,zz),time)
  enddo

  ii=0
  do kk=1,nl
    do jj=1,nyy
      do i2=1,nxx
        ii = ii+1
        do n=1,nens
          xpold(ii,n) = psin(i2,jj,kk,n)
        enddo
      enddo
    enddo
  enddo
  !  call randomv(randomx,m)
  !  xp(:,k) = xmp(:,1) + ivar * randomx
  !enddo 
  !xpold = xp

  !write(40,*)xtrue(1,:)    !FRP
  !write(20,*)xpold(:,:)  !QG
  write(*, *) '131th element of vector xpold at time:', time, '=', xpold(131,1)     !131th element for nl=1 (the first non-zero element for the random field) (nens=1) 
  !write(*, *) 'Vector xpold:', xpold(4226,1)    !First element for nl=2 (nens=1)

enddo

print *, " Program is finished."


! Prepare synchronisation, generate ensemble over full synchronisation window
!xxm = 0.
!do k = 1,nens
!  tau = 1
!  xhulp = xp(:,k)
!  do t = 1,ntau*tsteps
!    call L96RK4(xhulp,xzero,0)
!    call randomv(randomx,m)
    !xhulp = xhulp+matmul(Qhalf,randomx)
    !if (mod(t,tsteps).eq.0) then
    !  xx(:,k,tau) = xhulp(:)
    !  xxm(:,tau) = xxm(:,tau)+xhulp(:)
    !  tau = tau + 1
    !endif
!  enddo
!enddo

! Main time loop
!nn = 1             ! Multiplication factor for synchronization term, increases with time in window
!w = 0              ! - log proposal weight
!ycount = 1
!do t = 2,tstate

!  write(50,*)xp(1,:)   !FRP

!  if(mod(t,tsteps).ne.0) then           ! at any time step exept assimilation time step

    ! propagate ensemble and calculate raw weights
!    xmp(:,t) = 0.
!    do k = 1,nens
!      dxx(:) = dx(:,k) ! * real(nn)       ! Synchronization term 
!      call L96RK4(xp(:,k),dxx(:),0)
!      call randomv(randomx,m)
!      xhulp = matmul(Qhalf,randomx)
!      xp(:,k) = xp(:,k)+ xhulp
!      xmp(:,t) = xmp(:,t) + xp(:,k)
!    enddo
!    xmp(:,t) = xmp(:,t)/nens
!    do k=1,nens
!      xinflation = (inflation-1.)*(xp(:,k)-xmp(:,t))
!      call randomv(randomx,m)
!      xhulp = matmul(Qhalf,randomx)
!      xp(:,k) = xp(:,k) + xinflation 

! Save RMSE and ensemble spread in file 10
!do t=1,tstate
!  write(10,*)t,sx(t),var(t)
!enddo
!write(*,*)sqrt(sum(sx*sx)/tstate), sqrt(sum(var*var)/tstate)
!!write(10,*)sqrt(sum(sx*sx)/tstate), sqrt(sum(var*var)/tstate)  !!FRP

end Program SPF 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

