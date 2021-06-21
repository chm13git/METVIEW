Module Sizes

    integer, parameter  ::  nxx = 65                ! zonal (x) dimension of the system
    integer, parameter  ::  nyy = 65                ! medridional (y) dimension of the system
    integer, parameter  ::  nl  = 2                  ! number of vertical layers (has to be 2 for this code)
    integer, parameter  ::  nx = 70               ! zonal (x) dimension of random field (has to be even) 
    integer, parameter  ::  ny = 70                ! medridional (y) dimension of random field (has to be even)
    integer, parameter  ::  nens = 10               ! QG - number of ensemble members
    integer, parameter  ::  ntotal = nxx*nyy*nl          ! QG - total number of grid points
    real, parameter     ::  dx = 1.e5                ! grid distance in m
    integer, parameter  ::  time_to_start = 0        ! start time of integration (should be zero for now)
    integer, parameter  ::  time_to_stop = 20      ! end time of integration in units of time steps
    real, parameter     ::  dt =  1800               ! time step in seconds
    real, parameter     ::  diff = 1.e5              ! diffusion coefficient in m^2/s in laplacian for vorticity. 
    real, parameter     ::  lengthscale = 3.       ! decorrelation lengthscale for random fields in grid points!
                                                     ! Note that a Shapiro filter is also active.

!   arrays used by many subroutines
    real :: f(nyy),H1,H2,F1,F2,lambda
    real :: sh(0:16)
    real :: c1
    real :: cbt,sumpsibt,sumpsiboundarybt
    real :: cbc1,cbc2,n1old,n2old,n3old,a1old,a2old
!    real :: psiRandom(nxx,nyy,nl)
    real :: xr(nx,ny),qo(nl),alpha

End Module Sizes



