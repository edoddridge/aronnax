!> @author
!> Ed Doddridge
!
!> Minimalist Isopycnal Model (MIM) with n layers
!!
!
!>
!>     @mainpage Documentation for MIM.f90
!>
!>     @section Overview
!>     This model is an isopycnal model on an Arakawa C-grid with n
!>     layers and arbitrary bathymetry.
!>
!>
!>
!>    @section Grid
!>
!>    /\ ------------
!>    |  |          |
!>    |  |          |
!>    dy U    H     |
!>    |  |          |
!>    |  |          |
!>    \/ Z----V------
!>        <---dx---->
!>
!>    H: tracer point - thickness, Bernoulli potential
!>    U: velocity point - u and v
!>    Z: vorticity point - zeta
!>

program MIM
  implicit none

  integer, parameter :: nx = 300 !< number of x grid points
  integer, parameter :: ny = 300 !< number of y grid points
  integer, parameter :: layers = 2 !< number of active layers in the model

  ! Layer thickness (h)
  double precision :: h(0:nx+1, 0:ny+1, layers)
  double precision :: dhdt(0:nx+1, 0:ny+1, layers)
  double precision :: dhdtold(0:nx+1, 0:ny+1, layers)
  double precision :: dhdtveryold(0:nx+1, 0:ny+1, layers)
  double precision :: hnew(0:nx+1, 0:ny+1, layers)
  ! for initialisation
  double precision :: hhalf(0:nx+1, 0:ny+1, layers)
  ! for saving average fields
  double precision :: hav(0:nx+1, 0:ny+1, layers)

  ! Velocity component u
  double precision :: u(0:nx+1, 0:ny+1, layers)
  double precision :: dudt(0:nx+1, 0:ny+1, layers)
  double precision :: dudtold(0:nx+1, 0:ny+1, layers)
  double precision :: dudtveryold(0:nx+1, 0:ny+1, layers)
  double precision :: unew(0:nx+1, 0:ny+1, layers)
  double precision :: dudt_bt(0:nx+1, 0:ny+1)
  ! barotropic velocity components (for pressure solver)
  double precision :: ub(0:nx+1, 0:ny+1)
  ! for initialisation
  double precision :: uhalf(0:nx+1, 0:ny+1, layers)
  ! for saving average fields
  double precision :: uav(0:nx+1, 0:ny+1, layers)

  ! Velocity component v
  double precision :: v(0:nx+1, 0:ny+1, layers)
  double precision :: dvdt(0:nx+1, 0:ny+1, layers)
  double precision :: dvdtold(0:nx+1, 0:ny+1, layers)
  double precision :: dvdtveryold(0:nx+1, 0:ny+1, layers)
  double precision :: vnew(0:nx+1, 0:ny+1, layers)
  double precision :: dvdt_bt(nx+1, ny)
  ! barotropic velocity components (for pressure solver)
  double precision :: vb(nx+1, ny)
  ! for initialisation
  double precision :: vhalf(0:nx+1, 0:ny+1, layers)
  ! for saving average fields
  double precision :: vav(0:nx+1, 0:ny+1, layers)

  ! Free surface (eta)
  double precision :: eta(0:nx+1, 0:ny+1)
  double precision :: etastar(0:nx+1, 0:ny+1)
  double precision :: etanew(0:nx+1, 0:ny+1)
  ! for saving average fields
  double precision :: etaav(0:nx+1, 0:ny+1)

  ! Bathymetry
  character(30) :: depthFile
  double precision :: depth(0:nx+1, 0:ny+1)
  double precision :: H0 ! default depth in no file specified
  ! Pressure solver variables
  double precision :: a(5, nx, ny)

  ! Bernoulli potential and relative vorticity
  double precision :: b(0:nx+1, 0:ny+1, layers)
  double precision :: zeta(0:nx+1, 0:ny+1, layers)


  ! Grid
  double precision :: dx, dy
  double precision :: wetmask(0:nx+1, 0:ny+1)
  double precision :: hfacW(0:nx+1, 0:ny+1)
  double precision :: hfacE(0:nx+1, 0:ny+1)
  double precision :: hfacN(0:nx+1, 0:ny+1)
  double precision :: hfacS(0:nx+1, 0:ny+1)
  ! Coriolis parameter at u and v grid-points respectively
  double precision :: fu(0:nx+1, 0:ny+1)
  double precision :: fv(0:nx+1, 0:ny+1)
  ! File names to read them from
  character(30) :: fUfile, fVfile
  character(30) :: wetMaskFile

  ! Numerics
  double precision :: pi, dt
  double precision :: au, ah(layers), ar, botDrag
  double precision :: slip
  double precision :: hmin
  integer nTimeSteps
  double precision :: dumpFreq, avFreq
  integer counter
  integer nwrite, avwrite
  double precision :: zeros(layers)
  integer maxits
  double precision :: eps, rjac
  double precision :: freesurfFac
  double precision :: h_norming(0:nx+1, 0:ny+1)
  double precision :: thickness_error

  ! Model
  double precision :: hmean(layers)
  ! Switch for using n + 1/2 layer physics, or using n layer physics
  logical :: RedGrav

  ! Loop variables
  integer :: i, j, k, n

  ! Character variable for numbering the outputs
  character(10) :: num


  ! Physics
  double precision :: g_vec(layers), rho0

  ! Wind
  double precision :: wind_x(0:nx+1, 0:ny+1)
  double precision :: wind_y(0:nx+1, 0:ny+1)
  double precision :: base_wind_x(0:nx+1, 0:ny+1)
  double precision :: base_wind_y(0:nx+1, 0:ny+1)
  logical :: UseSinusoidWind
  logical :: UseStochWind
  logical :: DumpWind
  double precision :: wind_alpha, wind_beta, wind_period, wind_t_offset
  integer :: n_stoch_wind
  double precision :: stoch_wind_mag
  character(30) :: wind_mag_time_series_file
  double precision, dimension(:), allocatable :: wind_mag_time_series


  ! Sponge
  double precision :: spongeHTimeScale(0:nx+1, 0:ny+1, layers)
  double precision :: spongeUTimeScale(0:nx+1, 0:ny+1, layers)
  double precision :: spongeVTimeScale(0:nx+1, 0:ny+1, layers)
  double precision :: spongeH(0:nx+1, 0:ny+1, layers)
  double precision :: spongeU(0:nx+1, 0:ny+1, layers)
  double precision :: spongeV(0:nx+1, 0:ny+1, layers)
  character(30) :: spongeHTimeScaleFile
  character(30) :: spongeUTimeScaleFile
  character(30) :: spongeVTimeScaleFile
  character(30) :: spongeHfile
  character(30) :: spongeUfile
  character(30) :: spongeVfile

  ! Input files
  character(30) :: initUfile, initVfile, initHfile, initEtaFile
  character(30) :: zonalWindFile, meridionalWindFile

  ! Set default values here

  ! TODO Possibly wait until the model is split into multiple files,
  ! then hide the long unsightly code there.

  namelist /NUMERICS/ au, ah, ar, botDrag, dt, slip, nTimeSteps, &
      dumpFreq, avFreq, hmin, maxits, freesurfFac, eps, &
      thickness_error

  namelist /MODEL/ hmean, depthFile, H0, RedGrav

  namelist /SPONGE/ spongeHTimeScaleFile, spongeUTimeScaleFile, &
      spongeVTimeScaleFile, spongeHfile, spongeUfile, spongeVfile

  namelist /PHYSICS/ g_vec, rho0

  namelist /GRID/ dx, dy, fUfile, fVfile, wetMaskFile

  namelist /INITIAL_CONDITONS/ initUfile, initVfile, initHfile, initEtaFile

  namelist /EXTERNAL_FORCING/ zonalWindFile, meridionalWindFile, &
      UseSinusoidWind, UseStochWind, wind_alpha, wind_beta, &
      wind_period, wind_t_offset, DumpWind, &
      wind_mag_time_series_file

  open(unit=8, file="parameters.in", status='OLD', recl=80)
  read(unit=8, nml=NUMERICS)
  read(unit=8, nml=MODEL)
  read(unit=8, nml=SPONGE)
  read(8, nml=PHYSICS)
  read(8, nml=GRID)
  read(8, nml=INITIAL_CONDITONS)
  read(8, nml=EXTERNAL_FORCING)
  close(unit=8)

  nwrite = int(dumpFreq/dt)
  avwrite = int(avFreq/dt)

  ! Pi, the constant
  pi = 3.1415926535897932384
  ! Zero vector - for internal use only
  zeros = 0d0

  ! Check that time-dependent forcing flags have been set correctly
  if (UseSinusoidWind .and. UseStochWind)  then
    ! Can't use both sinusiodal and stochastic wind variation.
    ! Write a file saying so
    open(unit=99, file='errors.txt', action="write", status="replace", &
        form="formatted")
    write(99, *) "Can't have both stochastic and sinusoidally varying &
        &wind forcings. Choose one."
    close(unit=99)

    ! Print it on the screen
    print *, "Can't have both stochastic and sinusoidally varying &
        &wind forcings. Choose one."
    ! Stop the program
    stop
  end if

  ! Read in arrays from the input files
  call read_input_fileU(initUfile, u, 0.d0, nx, ny, layers)
  call read_input_fileV(initVfile, v, 0.d0, nx, ny, layers)
  call read_input_fileH(initHfile, h, hmean, nx, ny, layers)

  if (.not. RedGrav) then
    call read_input_fileH_2D(depthFile, depth, H0, nx, ny)
    call read_input_fileH_2D(initEtaFile, eta, 0.d0, nx, ny)
    ! Check that depth is positive - it must be greater than zero
    if (minval(depth) .lt. 0) then
      print *, "depths must be positive - fix this and try again"
      stop
    end if
  end if

  ! TODO Should probably check that bathymetry file and layer
  ! thicknesses are consistent with each other.

  call read_input_fileU(fUfile, fu, 0.d0, nx, ny, 1)
  call read_input_fileV(fVfile, fv, 0.d0, nx, ny, 1)

  call read_input_fileU(zonalWindFile, base_wind_x, 0.d0, nx, ny, 1)
  call read_input_fileV(meridionalWindFile, base_wind_y, 0.d0, nx, ny, 1)
  call read_input_file_time_series(wind_mag_time_series_file, &
   wind_mag_time_series, 1.d0, nTimeSteps)

  wind_x = base_wind_x*wind_mag_time_series(1)
  wind_y = base_wind_y*wind_mag_time_series(1)

  call read_input_fileH(spongeHTimeScaleFile, spongeHTimeScale, &
      zeros, nx, ny, layers)
  call read_input_fileH(spongeHfile, spongeH, hmean, nx, ny, layers)
  call read_input_fileU(spongeUTimeScaleFile, spongeUTimeScale, &
      0.d0, nx, ny, layers)
  call read_input_fileU(spongeUfile, spongeU, 0.d0, nx, ny, layers)
  call read_input_fileV(spongeVTimeScaleFile, spongeVTimeScale, &
      0.d0, nx, ny, layers)
  call read_input_fileV(spongeVfile, spongeV, 0.d0, nx, ny, layers)
  call read_input_fileH_2D(wetMaskFile, wetmask, 1.d0, nx, ny)

  ! Initialise the average fields
  if (avwrite .ne. 0) then
    hav = 0.0
    uav = 0.0
    vav = 0.0
    if (.not. RedGrav) then
      etaav = 0.0
    end if
  end if

  ! For now enforce wetmask to have zeros around the edge.
  wetmask(0, :) = 0d0
  wetmask(nx+1, :) = 0d0
  wetmask(:, 0) = 0d0
  wetmask(:, ny+1) = 0d0
  ! TODO, this will change one day when the model can do periodic
  ! boundary conditions.

  call calc_boundary_masks(wetmask, hfacW, hfacE, hfacS, hfacN, nx, ny)

  ! Apply the boundary conditions
  ! - Enforce no normal flow boundary condition
  !   and no flow in dry cells.
  ! - no/free-slip is done inside the dudt and dvdt subroutines.
  ! - hfacW and hfacS are zero where the transition between
  !   wet and dry cells occurs.
  ! - wetmask is 1 in wet cells, and zero in dry cells.
  do k = 1, layers
    ! h(:, :, k) = h(:, :, k)*wetmask(:, :)
    u(:, :, k) = u(:, :, k) * hfacW * wetmask(:, :)
    v(:, :, k) = v(:, :, k) * hfacS * wetmask(:, :)
  end do

  ! ! If the winds are static, then set wind_ = base_wind_
  ! if (.not. UseSinusoidWind .and. .not. UseStochWind)  then
  !   wind_x = base_wind_x
  !   wind_y = base_wind_y
  ! end if

  ! ! Initialise random numbers for stochastic winds
  ! if (UseStochWind) then
  !   ! This ensures a different series of random numbers each time the
  !   ! model is run.
  !   call ranseed()
  !   ! Number of timesteps between updating the perturbed wind field.
  !   n_stoch_wind = int(wind_period/dt)
  ! end if

  if (.not. RedGrav) then
    ! Initialise arrays for pressure solver
    ! a = derivatives of the depth field
    do j = 1, ny
      do i = 1, nx
        a(1,i,j) = g_vec(1)*0.5*(depth(i+1,j)+depth(i,j))/dx**2
        a(2,i,j) = g_vec(1)*0.5*(depth(i,j+1)+depth(i,j))/dy**2
        a(3,i,j) = g_vec(1)*0.5*(depth(i,j)+depth(i-1,j))/dx**2
        a(4,i,j) = g_vec(1)*0.5*(depth(i,j)+depth(i,j-1))/dy**2
      end do
    end do
    do j = 1, ny
      a(1, nx, j) = 0.0
      a(3, 1, j) = 0.0
    end do
    do i = 1, nx
      a(2, i, ny) = 0.0
      a(4, i, 1) = 0.0
    end do
    do j = 1, ny
      do i = 1, nx
        a(5,i,j) = -a(1,i,j)-a(2,i,j)-a(3,i,j)-a(4,i,j)
      end do
    end do

    ! Calculate the spectral radius of the grid for use by the
    ! successive over-relaxation scheme
    rjac = (cos(pi/real(nx))*dy**2+cos(pi/real(ny))*dx**2) &
           /(dx**2+dy**2)
    ! If peridodic boundary conditions are ever implemented, then pi ->
    ! 2*pi in this calculation

    ! Check that the supplied free surface anomaly and layer
    ! thicknesses are consistent
    h_norming = (freesurfFac*eta+depth)/sum(h,3)
    do k = 1, layers
      h(:, :, k) = h(:, :, k) * h_norming
    end do

    if (maxval(abs(h_norming - 1d0)).gt. thickness_error) then
      print *, 'inconsistency between h and eta (in %):', &
          maxval(abs(h_norming - 1d0))*100d0
    end if
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Initialisation of the model STARTS HERE                              !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Do two initial time steps with Runge-Kutta second-order.
  ! These initialisation steps do NOT use or update the free surface.
  !
  ! -------------------------- negative 2 time step ---------------------------
  ! Code to work out dhdtveryold, dudtveryold and dvdtveryold

  ! Calculate baroclinic Bernoulli potential
  if (RedGrav) then
    call evaluate_b_RedGrav(b, h, u, v, nx, ny, layers, g_vec)
  else
    call evaluate_b_iso(b, h, u, v, nx, ny, layers, g_vec, depth)
  end if

  ! Calculate relative vorticity
  call evaluate_zeta(zeta, u, v, nx, ny, layers, dx, dy)

  ! Calculate dhdt, dudt, dvdt at current time step
  call evaluate_dhdt(dhdtveryold, h, u, v, ah, dx, dy, nx, ny, layers, &
      spongeHTimeScale, spongeH, wetmask, RedGrav)

  call evaluate_dudt(dudtveryold,  h, u, v, b, zeta, wind_x, fu, au, ar, &
      slip, dx, dy, hfacN, hfacS, nx, ny, layers, rho0, &
      spongeUTimeScale, spongeU, RedGrav, botDrag)

  call evaluate_dvdt(dvdtveryold,  h, u, v, b, zeta, wind_y, fv, &
      au, ar, slip, dx, dy, hfacW, hfacE, nx, ny, layers, rho0, &
      spongeVTimeScale, spongeV, RedGrav, botDrag)

  ! Calculate the values at half the time interval with Forward Euler
  hhalf = h+0.5d0*dt*dhdtveryold
  uhalf = u+0.5d0*dt*dudtveryold
  vhalf = v+0.5d0*dt*dvdtveryold

  ! Calculate Bernoulli potential
  if (RedGrav) then
    call evaluate_b_RedGrav(b, hhalf, uhalf, vhalf, nx, ny, layers, g_vec)
  else
    call evaluate_b_iso(b, hhalf, uhalf, vhalf, nx, ny, layers, g_vec, depth)
  end if

  ! Calculate relative vorticity
  call evaluate_zeta(zeta, uhalf, vhalf, nx, ny, layers, dx, dy)

  ! Now calculate d/dt of u, v, h and store as dhdtveryold, dudtveryold
  ! and dvdtveryold
  call evaluate_dhdt(dhdtveryold, hhalf, uhalf, vhalf, ah, dx, dy, &
      nx, ny, layers, spongeHTimeScale, spongeH, wetmask, RedGrav)

  call evaluate_dudt(dudtveryold, hhalf, uhalf, vhalf, b, zeta, &
      wind_x, fu, au, ar, slip, dx, dy, hfacN, hfacS, nx, ny, layers, rho0, &
      spongeUTimeScale, spongeU, RedGrav, botDrag)

  call evaluate_dvdt(dvdtveryold, hhalf, uhalf, vhalf, b, zeta, &
      wind_y, fv, au, ar, slip, dx, dy, hfacW, hfacE, nx, ny, layers, rho0, &
      spongeVTimeScale, spongeV, RedGrav, botDrag)
  ! These are the values to be stored in the 'veryold' variables ready
  ! to start the proper model run.

  ! Calculate h, u, v with these tendencies
  h = h + dt*dhdtveryold
  u = u + dt*dudtveryold
  v = v + dt*dvdtveryold

  ! Apply the boundary conditions
  ! - Enforce no normal flow boundary condition
  !   and no flow in dry cells.
  ! - no/free-slip is done inside the dudt and dvdt subroutines.
  ! - hfacW and hfacS are zero where the transition between
  !   wet and dry cells occurs.
  ! - wetmask is 1 in wet cells, and zero in dry cells.
  do k = 1, layers
    u(:, :, k) = u(:, :, k) * hfacW * wetmask(:, :)
    v(:, :, k) = v(:, :, k) * hfacS * wetmask(:, :)
  end do

  ! Wrap fields around for periodic simulations
  call wrap_fields_3D(u, nx, ny, layers)
  call wrap_fields_3D(v, nx, ny, layers)
  call wrap_fields_3D(h, nx, ny, layers)

  ! -------------------------- negative 1 time step ---------------------------
  ! Code to work out dhdtold, dudtold and dvdtold

  ! Calculate Bernoulli potential
  if (RedGrav) then
    call evaluate_b_RedGrav(b, h, u, v, nx, ny, layers, g_vec)
  else
    call evaluate_b_iso(b, h, u, v, nx, ny, layers, g_vec, depth)
  end if

  ! Calculate relative vorticity
  call evaluate_zeta(zeta, u, v, nx, ny, layers, dx, dy)

  ! Calculate dhdt, dudt, dvdt at current time step
  call evaluate_dhdt(dhdtold, h, u, v, ah, dx, dy, nx, ny, layers, &
      spongeHTimeScale, spongeH, wetmask, RedGrav)

  call evaluate_dudt(dudtold, h, u, v, b, zeta, wind_x, fu, &
      au, ar, slip, dx, dy, hfacN, hfacS, nx, ny, layers, rho0, &
      spongeUTimeScale, spongeU, RedGrav, botDrag)

  ! If the wind is changed to be meridional this code will need changing
  call evaluate_dvdt(dvdtold, h, u, v, b, zeta, wind_y, fv, &
      au, ar, slip, dx, dy, hfacW, hfacE, nx, ny, layers, rho0, &
      spongeVTimeScale, spongeV, RedGrav, botDrag)

  ! Calculate the values at half the time interval with Forward Euler
  hhalf = h+0.5d0*dt*dhdtold
  uhalf = u+0.5d0*dt*dudtold
  vhalf = v+0.5d0*dt*dvdtold

  ! Calculate Bernoulli potential
  if (RedGrav) then
    call evaluate_b_RedGrav(b, hhalf, uhalf, vhalf, nx, ny, layers, g_vec)
  else
    call evaluate_b_iso(b, hhalf, uhalf, vhalf, nx, ny, layers, g_vec, depth)
  end if

  ! Calculate relative vorticity
  call evaluate_zeta(zeta, uhalf, vhalf, nx, ny, layers, dx, dy)

  ! Now calculate d/dt of u, v, h and store as dhdtold, dudtold and dvdtold
  call evaluate_dhdt(dhdtold, hhalf, uhalf, vhalf, ah, dx, dy, &
      nx, ny, layers, spongeHTimeScale, spongeH, wetmask, RedGrav)
  call evaluate_dudt(dudtold, hhalf, uhalf, vhalf, b, zeta, wind_x, &
      fu, au, ar, slip, dx, dy, hfacN, hfacS, nx, ny, layers, rho0, &
      spongeUTimeScale, spongeU, RedGrav, botDrag)
  ! If the wind is changed to be meridional this code will need changing
  call evaluate_dvdt(dvdtold, hhalf, uhalf, vhalf, b, zeta, wind_y, &
      fv, au, ar, slip, dx, dy, hfacW, hfacE, nx, ny, layers, rho0, &
      spongeVTimeScale, spongeV, RedGrav, botDrag)
  ! These are the values to be stored in the 'old' variables ready to start
  ! the proper model run.

  ! Calculate h, u, v with these tendencies
  h = h + dt*dhdtold
  u = u + dt*dudtold
  v = v + dt*dvdtold

  ! Apply the boundary conditions
  ! - Enforce no normal flow boundary condition
  !   and no flow in dry cells.
  ! - no/free-slip is done inside the dudt and dvdt subroutines.
  ! - hfacW and hfacS are zero where the transition between
  !   wet and dry cells occurs.
  ! - wetmask is 1 in wet cells, and zero in dry cells.
  do k = 1, layers
    u(:, :, k) = u(:, :, k) * hfacW * wetmask(:, :)
    v(:, :, k) = v(:, :, k) * hfacS * wetmask(:, :)
  end do

  ! Wrap fields around for periodic simulations
  call wrap_fields_3D(u, nx, ny, layers)
  call wrap_fields_3D(v, nx, ny, layers)
  call wrap_fields_3D(h, nx, ny, layers)

  ! Now the model is ready to start.
  ! - We have h, u, v at the zeroth time step, and the tendencies at
  !   two older time steps.
  ! - The model then solves for the tendencies at the current step
  !   before solving for the fields at the next time step.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! MAIN LOOP OF THE MODEL STARTS HERE                                    !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do n = 1, nTimeSteps

    ! ! Time varying winds
    ! if (UseSinusoidWind .eqv. .true.) then
    !   wind_x = base_wind_x*(wind_alpha +  &
    !       wind_beta*sin(((2d0*pi*n*dt)/wind_period) - wind_t_offset))
    !   wind_y = base_wind_y*(wind_alpha +  &
    !       wind_beta*sin(((2d0*pi*n*dt)/wind_period) - wind_t_offset))
    ! else if (UseStochWind .eqv. .true.) then
    !   if (mod(n-1, n_stoch_wind).eq.0) then
    !     ! Gives a pseudorandom number in range 0 <= x < 1
    !     call random_number(stoch_wind_mag)
    !     ! Convert to -1 <= x < 1
    !     stoch_wind_mag = (stoch_wind_mag - 0.5d0)*2d0
    !     wind_x = base_wind_x*(wind_alpha +  &
    !         wind_beta*stoch_wind_mag)
    !     wind_y = base_wind_y*(wind_alpha +  &
    !         wind_beta*stoch_wind_mag)
    !   end if
    ! end if

    wind_x = base_wind_x*wind_mag_time_series(n)
    wind_y = base_wind_y*wind_mag_time_series(n)

    ! Calculate Bernoulli potential
    if (RedGrav) then
      call evaluate_b_RedGrav(b, h, u, v, nx, ny, layers, g_vec)
    else
      call evaluate_b_iso(b, h, u, v, nx, ny, layers, g_vec, depth)
    end if

    ! Calculate relative vorticity
    call evaluate_zeta(zeta, u, v, nx, ny, layers, dx, dy)

    ! Calculate dhdt, dudt, dvdt at current time step
    call evaluate_dhdt(dhdt, h, u, v, ah, dx, dy, nx, ny, layers, &
        spongeHTimeScale, spongeH, wetmask, RedGrav)

    call evaluate_dudt(dudt, h, u, v, b, zeta, wind_x, fu, au, ar, slip, &
        dx, dy, hfacN, hfacS, nx, ny, layers, rho0, &
        spongeUTimeScale, spongeU, RedGrav, botDrag)

    call evaluate_dvdt(dvdt, h, u, v, b, zeta, wind_y, fv, au, ar, slip, &
        dx, dy, hfacW, hfacE, nx, ny, layers, rho0, &
        spongeVTimeScale, spongeV, RedGrav, botDrag)

    ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time with
    ! the Adams-Bashforth third order linear multistep method

    unew = u+dt*(23d0*dudt - 16d0*dudtold + 5d0*dudtveryold)/12d0
    vnew = v+dt*(23d0*dvdt - 16d0*dvdtold + 5d0*dvdtveryold)/12d0
    hnew = h+dt*(23d0*dhdt - 16d0*dhdtold + 5d0*dhdtveryold)/12d0

    ! Apply the boundary conditions
    ! - Enforce no normal flow boundary condition
    !   and no flow in dry cells.
    ! - no/free-slip is done inside the dudt and dvdt subroutines.
    ! - hfacW and hfacS are zero where the transition between
    !   wet and dry cells occurs.
    ! - wetmask is 1 in wet cells, and zero in dry cells.
    do k = 1, layers
      unew(:, :, k) = unew(:, :, k) * hfacW * wetmask(:, :)
      vnew(:, :, k) = vnew(:, :, k) * hfacS * wetmask(:, :)
    end do

    ! Do the isopycnal layer physics
    if (.not. RedGrav) then
      ! Calculate the barotropic velocities
      call calc_baro_u(ub, unew, hnew, eta, freesurfFac, nx, ny, layers)
      call calc_baro_v(vb, vnew, hnew, eta, freesurfFac, nx, ny, layers)

      ! Calculate divergence of ub and vb, and solve for the pressure
      ! field that removes it
      call calc_eta_star(ub, vb, eta, etastar, freesurfFac, nx, ny, dx, dy, dt)
      ! print *, maxval(abs(etastar))

      ! Prevent barotropic signals from bouncing around outside the
      ! wet region of the model.
      ! etastar = etastar*wetmask

      call SOR_solver(a, etanew, etastar, freesurfFac, nx, ny, &
          dt, rjac, eps, maxits, n)
      ! print *, maxval(abs(etanew))

      etanew = etanew*wetmask
      call wrap_fields_2D(etanew, nx, ny)

      ! Now update the velocities using the barotropic tendency due to
      ! the pressure gradient.
      dudt_bt = 0d0
      dvdt_bt = 0d0

      do i = 1, nx
        do j = 0, ny
          do k = 1, layers
            dudt_bt(i,j) = -g_vec(1)*(etanew(i,j) - etanew(i-1,j))/(dx)
            unew(i,j,k) = unew(i,j,k) + dt*dudt_bt(i,j) ! 23d0*dt*dudt_bt(i,j)/12d0
          end do
        end do
      end do

      do i = 0, nx
        do j = 1, ny
          do k = 1, layers
            dvdt_bt(i,j) = -g_vec(1)*(etanew(i,j) - etanew(i,j-1))/(dy)
            vnew(i,j,k) = vnew(i,j,k) + dt*dvdt_bt(i,j)! 23d0*dt*dvdt_bt(i,j)/12d0
          end do
        end do
      end do

      ! We now have correct velocities at the next time step, but the
      ! layer thicknesses were updated using the velocities without
      ! the barotropic pressure contribution. Force consistency
      ! between layer thicknesses and ocean depth by scaling
      ! thicknesses to agree with free surface.

      h_norming = (freesurfFac*etanew + depth) / sum(hnew, 3)
      do k = 1, layers
        hnew(:, :, k) = hnew(:, :, k) * h_norming
      end do

      if (maxval(abs(h_norming - 1d0)) .gt. thickness_error) then
        print *, 'inconsistency between h and eta (in %).', &
            maxval(abs(h_norming - 1d0))*100
      end if

      ! Apply the boundary conditions
      ! - Enforce no normal flow boundary condition
      !   and no flow in dry cells.
      ! - no/free-slip is done inside the dudt and dvdt subroutines.
      ! - hfacW and hfacS are zero where the transition between
      !   wet and dry cells occurs.
      ! - wetmask is 1 in wet cells, and zero in dry cells.
      do k = 1, layers
        unew(:, :, k) = unew(:, :, k) * hfacW * wetmask(:, :)
        vnew(:, :, k) = vnew(:, :, k) * hfacS * wetmask(:, :)
      end do
    end if

    ! Code to stop layers from getting too thin
    counter = 0

    do k = 1, layers
      do j = 1, ny
        do i = 1, nx
          if (hnew(i, j, k) .lt. hmin) then
            hnew(i, j, k) = hmin
            counter = counter + 1
            if (counter .eq. 1) then
              ! Write a file saying that the layer thickness value
              ! dropped below hmin and this line has been used.
              open(unit=10, file='layer thickness dropped below hmin.txt', &
                  action="write", status="unknown", &
                  form="formatted", position="append")
              write(10, 1111) n
1111          format("layer thickness dropped below hmin at time step ", 1i10.10)
              close(unit=10)
            end if
          end if
        end do
      end do
    end do

    ! Wrap fields around for periodic simulations
    call wrap_fields_3D(unew, nx, ny, layers)
    call wrap_fields_3D(vnew, nx, ny, layers)
    call wrap_fields_3D(hnew, nx, ny, layers)
    call wrap_fields_2D(etanew, nx, ny)

    ! Accumulate average fields
    if (avwrite .ne. 0) then
      hav = hav + hnew
      uav = uav + unew
      vav = vav + vnew
      if (.not. RedGrav) then
        etaav = eta + etanew
      end if
    end if

    ! Shuffle arrays: old -> very old,  present -> old, new -> present
    ! Height and velocity fields
    h = hnew
    u = unew
    v = vnew
    if (.not. RedGrav) then
      eta = etanew
    end if

    ! Tendencies (old -> very old)
    dhdtveryold = dhdtold
    dudtveryold = dudtold
    dvdtveryold = dvdtold

    ! Tendencies (current -> old)
    dudtold = dudt
    dvdtold = dvdt
    dhdtold = dhdt

    ! Now have new fields in main arrays and old fields in very old arrays

    ! ------------------------ Code for generating output ---------------------
    ! Write data to file?
    if (mod(n-1, nwrite) .eq. 0) then

      !
      ! ----------------------- Snapshot Output --------------------
      write(num, '(i10.10)') n

      ! Output the data to a file
      open(unit=10, status='replace', file='output/snap.h.'//num, &
          form='unformatted')
      write(10) h(1:nx, 1:ny, :)
      close(10)
      open(unit=10, status='replace', file='output/snap.u.'//num, &
          form='unformatted')
      write(10) u(1:nx+1, 1:ny, :)
      close(10)
      open(unit=10, status='replace', file='output/snap.v.'//num, &
          form='unformatted')
      write(10) v(1:nx, 1:ny+1, :)
      close(10)
      if (.not. RedGrav) then
        open(unit=10, status='replace', file='output/snap.eta.'//num, &
            form='unformatted')
        write(10) eta(1:nx, 1:ny)
        close(10)
      end if

      if (DumpWind .eqv. .true.) then
        open(unit=10, status='replace', file='output/wind_x.'//num, &
            form='unformatted')
        write(10) wind_x(1:nx+1, 1:ny)
        close(10)
        open(unit=10, status='replace', file='output/wind_y.'//num, &
            form='unformatted')
        write(10) wind_y(1:nx, 1:ny+1)
        close(10)
      end if

      ! Check if there are NaNs in the data
      call break_if_NaN(h, nx, ny, layers, n)
      ! call break_if_NaN(u, nx, ny, layers, n)
      ! call break_if_NaN(v, nx, ny, layers, n)

    end if

    ! ----------------------- Average Output -----------------
    if (avwrite .eq. 0) then
      ! OK
    else if (mod(n-1, avwrite) .eq. 0) then

      hav = hav/real(avwrite)
      uav = uav/real(avwrite)
      vav = vav/real(avwrite)

      write(num, '(i10.10)') n

      ! output the data to a file

      open(unit=10, status='replace', file='output/av.h.'//num, &
          form='unformatted')
      write(10) hav(1:nx, 1:ny, :)
      close(10)
      open(unit=10, status='replace', file='output/av.u.'//num, &
          form='unformatted')
      write(10) uav(1:nx+1, 1:ny, :)
      close(10)
      open(unit=10, status='replace', file='output/av.v.'//num, &
          form='unformatted')
      write(10) vav(1:nx, 1:ny+1, :)
      close(10)
      if (.not. RedGrav) then
        open(unit=10, status='replace', file='output/av.eta.'//num, &
            form='unformatted')
        write(10) etaav(1:nx, 1:ny)
        close(10)
      end if

      ! Check if there are NaNs in the data
      call break_if_NaN(h, nx, ny, layers, n)
      ! call break_if_NaN(u, nx, ny, layers, n)
      ! call break_if_NaN(v, nx, ny, layers, n)

      ! Reset average quantities
      hav = 0.0
      uav = 0.0
      vav = 0.0
      if (.not. RedGrav) then
        etaav = 0.0
      end if
      ! h2av = 0.0

    end if

  end do

  open(unit=10, file='run_finished.txt', action="write", status="unknown", &
      form="formatted", position="append")
  write(10, 1112) n
1112 format( "run finished at time step ", 1i10.10)
  close(unit=10)

  print *, 'Execution ended normally'

  stop 0

end program MIM


! ------------------------- Beginning of the subroutines ----------------------
!
! -----------------------------------------------------------------------------

!> Evaluate the Bornoulli Potential for n-layer physics.
!! B is evaluated at the tracer point, for each grid box.
subroutine evaluate_b_iso(b, h, u, v, nx, ny, layers, g_vec, depth)
  implicit none

  ! Evaluate the baroclinic component of the Bernoulli Potential
  ! (u dot u + Montgomery potential) in the n-layer physics, at centre
  ! of grid box

  double precision, intent(out) :: b(0:nx+1, 0:ny+1, layers) !< Bernoulli Potential
  integer, intent(in) :: nx !< number of x grid points
  integer, intent(in) :: ny !< number of y grid points
  integer, intent(in) :: layers !< number of layers
  integer i, j, k
  double precision, intent(in) :: depth(0:nx+1, 0:ny+1) !< total depth of fluid
  double precision z(0:nx+1, 0:ny+1, layers)
  double precision, intent(in) :: h(0:nx+1, 0:ny+1, layers) !< layer thicknesses
  double precision, intent(in) :: u(0:nx+1, 0:ny+1, layers) !< zonal velocities
  double precision, intent(in) :: v(0:nx+1, 0:ny+1, layers) !< meridional velocities
  double precision, intent(in) :: g_vec(layers) !< reduced gravity at each interface
  double precision M(0:nx+1, 0:ny+1, layers)

  ! Calculate layer interface locations
  z = 0d0
  z(:, :, layers) = -depth

  do k = 1, layers-1
    z(:, :, layers - k) = z(:, :, layers-k+1) + h(:, :, layers-k+1)
  end do

  ! Calculate Montogmery potential
  ! The following loop is to get the baroclinic Montgomery potential
  ! in each layer
  M = 0d0
  do k = 2, layers
    M(:, :, k) = M(:, :, k-1) + g_vec(k) * z(:, :, k-1)
  end do

  b = 0d0
  ! No baroclinic pressure contribution to the first layer Bernoulli potential
  ! (the barotropic pressure contributes, but that's not done here).
  ! do j = 1, ny-1
  !     do i = 1, nx-1
  !         b(i,j,1) = (u(i,j,1)**2+u(i+1,j,1)**2+v(i,j,1)**2+v(i,j+1,1)**2)/4.0d0
  !     end do
  ! end do

  ! For the rest of the layers we get a baroclinic pressure contribution
  do k = 1, layers ! move through the different layers of the model
    do j = 1, ny ! move through longitude
      do i = 1, nx ! move through latitude
        b(i,j,k) = M(i,j,k) + (u(i,j,k)**2+u(i+1,j,k)**2+v(i,j,k)**2+v(i,j+1,k)**2)/4.0d0
        ! Add the (u^2 + v^2)/2 term to the Montgomery Potential
      end do
    end do
  end do

  return
end subroutine evaluate_b_iso

! -----------------------------------------------------------------------------

subroutine evaluate_b_RedGrav(b, h, u, v, nx, ny, layers, gr)
  implicit none

  ! Evaluate Bernoulli Potential at centre of grid box
  integer nx, ny, layers
  integer i, j, k, l, m
  double precision h(0:nx+1, 0:ny+1, layers)
  double precision u(0:nx+1, 0:ny+1, layers)
  double precision v(0:nx+1, 0:ny+1, layers)
  double precision b(0:nx+1, 0:ny+1, layers)
  double precision gr(layers)
  double precision h_temp, b_proto

  b = 0d0

  do k = 1, layers ! move through the different layers of the model
    do j = 1, ny ! move through longitude
      do i = 1, nx ! move through latitude
        ! The following loops are to get the pressure term in the
        ! Bernoulli Potential
        b_proto = 0d0
        do l = k, layers
          h_temp = 0d0
          do m = 1, l
            h_temp = h_temp + h(i, j, m) ! sum up the layer thicknesses
          end do
          ! Sum up the product of reduced gravity and summed layer
          ! thicknesses to form the pressure componenet of the
          ! Bernoulli Potential term
          b_proto = b_proto + gr(l)*h_temp
        end do
        ! Add the (u^2 + v^2)/2 term to the pressure componenet of the
        ! Bernoulli Potential
        b(i,j,k) = b_proto + (u(i,j,k)**2+u(i+1,j,k)**2+v(i,j,k)**2+v(i,j+1,k)**2)/4.0d0
      end do
    end do
  end do

  return
end subroutine evaluate_b_RedGrav

! -----------------------------------------------------------------------------
!> Evaluate relative vorticity at lower left grid boundary (du/dy
!! and dv/dx are at lower left corner as well)
subroutine evaluate_zeta(zeta, u, v, nx, ny, layers, dx, dy)
  implicit none

  integer nx, ny, layers
  integer i, j, k
  double precision u(0:nx+1, 0:ny+1, layers)
  double precision v(0:nx+1, 0:ny+1, layers)
  double precision zeta(0:nx+1, 0:ny+1, layers)
  double precision dx, dy

  zeta = 0d0

  do k = 1, layers
    do j = 1, ny+1
      do i = 1, nx+1
        zeta(i,j,k) = (v(i,j,k)-v(i-1,j,k))/dx-(u(i,j,k)-u(i,j-1,k))/dy
      end do
    end do
  end do

  return
end subroutine evaluate_zeta

! -----------------------------------------------------------------------------
!> Calculate the tendency of layer thickness for each of the active layers
!! dh/dt is in the centre of each grid point.
subroutine evaluate_dhdt(dhdt, h, u, v, ah, dx, dy, nx, ny, layers, &
    spongeTimeScale, spongeH, wetmask, RedGrav)
  implicit none

  ! dhdt is evaluated at the centre of the grid box
  integer nx, ny, layers
  integer i, j, k
  double precision dhdt(0:nx+1, 0:ny+1, layers)
  double precision h(0:nx+1, 0:ny+1, layers)
  double precision u(0:nx+1, 0:ny+1, layers)
  double precision v(0:nx+1, 0:ny+1, layers)
  ! Thickness tendency due to thickness diffusion (equivalent to Gent
  ! McWilliams in a z coordinate model)
  double precision dhdt_GM(0:nx+1, 0:ny+1, layers)
  double precision spongeTimeScale(0:nx+1, 0:ny+1, layers)
  double precision spongeH(0:nx+1, 0:ny+1, layers)
  double precision wetmask(0:nx+1, 0:ny+1)
  double precision dx, dy
  double precision ah(layers)
  logical :: RedGrav

  ! Calculate tendency due to thickness diffusion (equivalent
  ! to GM in z coordinate model with the same diffusivity).
  dhdt_GM = 0d0

  ! Loop through all layers except lowest and calculate
  ! thickness tendency due to diffusive mass fluxes
  do k = 1, layers-1
    do j = 1, ny
      do i = 1, nx
        dhdt_GM(i,j,k) = &
            ah(k)*(h(i+1,j,k)*wetmask(i+1,j)    &
              + (1d0 - wetmask(i+1,j))*h(i,j,k) & ! reflect around boundary
              + h(i-1,j,k)*wetmask(i-1,j)       &
              + (1d0 - wetmask(i-1,j))*h(i,j,k) & ! refelct around boundary
              - 2*h(i,j,k))/(dx*dx)             & ! x-component

            + ah(k)*(h(i,j+1,k)*wetmask(i,j+1) &
              + (1d0 - wetmask(i,j+1))*h(i,j,k) & ! reflect value around boundary
              + h(i,j-1,k)*wetmask(i,j-1)       &
              + (1d0 - wetmask(i,j-1))*h(i,j,k) & ! reflect value around boundary
              - 2*h(i,j,k))/(dy*dy) ! y-component horizontal diffusion
      end do
    end do
  end do

  ! Now do the lowest active layer, k = layers. If using reduced
  ! gravity physics, this is unconstrained and calculated as above. If
  ! using n-layer physics it is constrained to balance the layers
  ! above it.
  if (RedGrav) then
    do j = 1, ny
      do i = 1, nx
        dhdt_GM(i,j,layers) = &
            ah(layers)*(h(i+1,j,layers)*wetmask(i+1,j)   &
              + (1d0 - wetmask(i+1,j))*h(i,j,layers)     & ! boundary
              + h(i-1,j,layers)*wetmask(i-1,j)           &
              + (1d0 - wetmask(i-1,j))*h(i,j,layers)     & ! boundary
              - 2*h(i,j,layers))/(dx*dx)                 & ! x-component

            + ah(layers)*(h(i,j+1,layers)*wetmask(i,j+1) &
              + (1d0 - wetmask(i,j+1))*h(i,j,layers)     & ! reflect value around boundary
              + h(i,j-1,layers)*wetmask(i,j-1)           &
              + (1d0 - wetmask(i,j-1))*h(i,j,layers)     & ! reflect value around boundary
              - 2*h(i,j,layers))/(dy*dy) ! y-component horizontal diffusion
      end do
    end do
  else if (.not. RedGrav) then ! using n-layer physics
    ! Calculate bottom layer thickness tendency to balance layers above.
    ! In the flat-bottomed case this will give the same answer.
    dhdt_GM(:,:,layers) = -sum(dhdt_GM(:,:,:layers-1), 3)
  end if

  ! Now add this to the thickness tendency due to the flow field and
  ! sponge regions
  dhdt = 0d0

  do k = 1, layers
    do j = 1, ny
      do i = 1, nx
        dhdt(i,j,k) = &
            dhdt_GM(i,j,k) & ! horizontal thickness diffusion
            - ((h(i,j,k)+h(i+1,j,k))*u(i+1,j,k) &
               - (h(i-1,j,k)+h(i,j,k))*u(i,j,k))/(dx*2d0) & ! d(hu)/dx
            - ((h(i,j,k)+h(i,j+1,k))*v(i,j+1,k) &
              - (h(i,j-1,k)+h(i,j,k))*v(i,j,k))/(dy*2d0)  & ! d(hv)/dy
            + spongeTimeScale(i,j,k)*(spongeH(i,j,k)-h(i,j,k)) ! forced relaxtion in the sponge regions.
      end do
    end do
  end do

  ! Make sure the dynamics are only happening in the wet grid points.
  do k = 1, layers
    dhdt(:, :, k) = dhdt(:, :, k) * wetmask
  end do

  call wrap_fields_3D(dhdt, nx, ny, layers)

  return
end subroutine evaluate_dhdt

! -----------------------------------------------------------------------------
!> Calculate the tendency of zonal velocity for each of the active layers

subroutine evaluate_dudt(dudt, h, u, v, b, zeta, wind_x, fu, &
    au, ar, slip, dx, dy, hfacN, hfacS, nx, ny, layers, rho0, &
    spongeTimeScale, spongeU, RedGrav, botDrag)
  implicit none

  ! dudt(i, j) is evaluated at the centre of the left edge of the grid
  ! box, the same place as u(i, j).
  integer nx, ny, layers
  integer i, j, k
  double precision dudt(0:nx+1, 0:ny+1, layers)
  double precision h(0:nx+1, 0:ny+1, layers)
  double precision u(0:nx+1, 0:ny+1, layers)
  double precision v(0:nx+1, 0:ny+1, layers)
  double precision fu(0:nx+1, 0:ny+1)
  double precision zeta(0:nx+1, 0:ny+1, layers)
  double precision b(0:nx+1, 0:ny+1, layers)
  double precision spongeTimeScale(0:nx+1, 0:ny+1, layers)
  double precision spongeU(0:nx+1, 0:ny+1, layers)
  double precision wind_x(0:nx+1, 0:ny+1)
  double precision dx, dy
  double precision au, ar, rho0, slip, botDrag
  double precision hfacN(0:nx+1, 0:ny+1)
  double precision hfacS(0:nx+1, 0:ny+1)
  logical :: RedGrav

  dudt = 0d0

  do k = 1, layers
    do i = 1, nx
      do j = 1, ny
        dudt(i,j,k) = au*(u(i+1,j,k)+u(i-1,j,k)-2.0d0*u(i,j,k))/(dx*dx) & ! x-component
            + au*(u(i,j+1,k)+u(i,j-1,k)-2.0d0*u(i,j,k) &
              ! boundary conditions
              + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacN(i,j))*u(i,j,k) &
              + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacS(i,j))*u(i,j,k))/(dy*dy) & ! y-component
              ! Together make the horizontal diffusion term
            + 0.25d0*(fu(i,j)+0.5d0*(zeta(i,j,k)+zeta(i,j+1,k))) &
              *(v(i-1,j,k)+v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k)) & ! vorticity term
            - (b(i,j,k) - b(i-1,j,k))/dx & ! Bernoulli potential term
            + spongeTimeScale(i,j,k)*(spongeU(i,j,k)-u(i,j,k)) ! forced relaxtion in the sponge regions
        if (k .eq. 1) then ! only have wind forcing on the top layer
          ! This will need refining in the event of allowing outcropping.
          dudt(i,j,k) = dudt(i,j,k) + wind_x(i,j)/(rho0*h(i,j,k)) ! wind forcing
        end if
        if (layers .gt. 1) then ! only evaluate vertical momentum diffusivity if more than 1 layer
          if (k .eq. 1) then ! adapt vertical momentum diffusivity for 2+ layer model -> top layer
            dudt(i,j,k) = dudt(i,j,k) - 1.0d0*ar*(u(i,j,k) - 1.0d0*u(i,j,k+1))
          else if (k .eq. layers) then ! bottom layer
            dudt(i,j,k) = dudt(i,j,k) - 1.0d0*ar*(u(i,j,k) - 1.0d0*u(i,j,k-1))
            if (.not. RedGrav) then
              ! add bottom drag here in isopycnal version
              dudt(i,j,k) = dudt(i,j,k) - 1.0d0*botDrag*(u(i,j,k))
            end if
          else ! mid layer/s
            dudt(i,j,k) = dudt(i,j,k) - &
                1.0d0*ar*(2.0d0*u(i,j,k) - 1.0d0*u(i,j,k-1) - 1.0d0*u(i,j,k+1))
          end if
        end if
      end do
    end do
  end do

  call wrap_fields_3D(dudt, nx, ny, layers)

  return
end subroutine evaluate_dudt

! -----------------------------------------------------------------------------
!> Calculate the tendency of meridional velocity for each of the active layers

subroutine evaluate_dvdt(dvdt, h, u, v, b, zeta, wind_y, fv, &
    au, ar, slip, dx, dy, hfacW, hfacE, nx, ny, layers, rho0, &
    spongeTimeScale, spongeV, RedGrav, botDrag)
  implicit none

  ! dvdt(i, j) is evaluated at the centre of the bottom edge of the
  ! grid box, the same place as v(i, j)
  integer nx, ny, layers
  integer i, j, k
  double precision dvdt(0:nx+1, 0:ny+1, layers)
  double precision h(0:nx+1, 0:ny+1, layers)
  double precision u(0:nx+1, 0:ny+1, layers)
  double precision v(0:nx+1, 0:ny+1, layers)
  double precision fv(0:nx+1, 0:ny+1)
  double precision zeta(0:nx+1, 0:ny+1, layers)
  double precision b(0:nx+1, 0:ny+1, layers)
  double precision spongeTimeScale(0:nx+1, 0:ny+1, layers)
  double precision spongeV(0:nx+1, 0:ny+1, layers)
  double precision wind_y(0:nx+1, 0:ny+1)
  double precision dx, dy
  double precision au, ar, slip, botDrag
  double precision rho0
  double precision hfacW(0:nx+1, 0:ny+1)
  double precision hfacE(0:nx+1, 0:ny+1)
  logical :: RedGrav

  dvdt = 0d0

  do k = 1, layers
    do j = 1, ny
      do i = 1, nx
        dvdt(i,j,k) = &
            au*(v(i+1,j,k)+v(i-1,j,k)-2.0d0*v(i,j,k) &
              ! boundary conditions
              + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacW(i,j))*v(i,j,k) &
              + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacE(i,j))*v(i,j,k))/(dx*dx) & !x-component
            + au*(v(i,j+1,k) + v(i,j-1,k) - 2.0d0*v(i,j,k))/(dy*dy) & ! y-component.
            ! Together these make the horizontal diffusion term
            - 0.25d0*(fv(i,j)+0.5d0*(zeta(i,j,k)+zeta(i+1,j,k))) &
              *(u(i,j-1,k)+u(i,j,k)+u(i+1,j-1,k)+u(i+1,j,k)) & !vorticity term
            - (b(i,j,k)-b(i,j-1,k))/dy & ! Bernoulli Potential term
            + spongeTimeScale(i,j,k)*(spongeV(i,j,k)-v(i,j,k)) ! forced relaxtion to vsponge (in the sponge regions)
        if (k .eq. 1) then ! only have wind forcing on the top layer
          ! This will need refining in the event of allowing outcropping.
          dvdt(i,j,k) = dvdt(i,j,k) + wind_y(i,j)/(rho0*h(i,j,k)) ! wind forcing
        end if
        if (layers .gt. 1) then ! only evaluate vertical momentum diffusivity if more than 1 layer
          if (k .eq. 1) then ! adapt vertical momentum diffusivity for 2+ layer model -> top layer
            dvdt(i,j,k) = dvdt(i,j,k) - 1.0d0*ar*(v(i,j,k) - 1.0d0*v(i,j,k+1))
          else if (k .eq. layers) then ! bottom layer
            dvdt(i,j,k) = dvdt(i,j,k) - 1.0d0*ar*(v(i,j,k) - 1.0d0*v(i,j,k-1))
            if (.not. RedGrav) then
              ! add bottom drag here in isopycnal version
              dvdt(i,j,k) = dvdt(i,j,k) - 1.0d0*botDrag*(v(i,j,k))
            end if
          else ! mid layer/s
            dvdt(i,j,k) = dvdt(i,j,k) - &
                1.0d0*ar*(2.0d0*v(i,j,k) - 1.0d0*v(i,j,k-1) - 1.0d0*v(i,j,k+1))
          end if
        end if
      end do
    end do
  end do

  call wrap_fields_3D(dvdt, nx, ny, layers)

  return
end subroutine evaluate_dvdt

! -----------------------------------------------------------------------------
!> Calculate the barotropic u velocity
subroutine calc_baro_u(ub, u, h, eta, freesurfFac, nx, ny, layers)
  implicit none

  integer nx, ny, layers
  integer i, j, k
  double precision h(0:nx+1, 0:ny+1, layers)
  double precision h_temp(0:nx+1, 0:ny+1, layers)
  double precision eta(0:nx+1, 0:ny+1)
  double precision freesurfFac
  double precision u(0:nx+1, 0:ny+1, layers)
  double precision ub(0:nx+1, 0:ny+1)

  ub = 0d0

  h_temp = h
  ! add free surface elevation to the upper layer
  h_temp(:, :, 1) = h(:, :, 1) + eta*freesurfFac

  do i = 1, nx
    do j = 1, ny
      do k = 1, layers
        ub(i,j) = ub(i,j) + u(i,j,k)*(h_temp(i,j,k)+h_temp(i-1,j,k))/2d0
      end do
    end do
  end do

  call wrap_fields_2D(ub, nx, ny)

  return
end subroutine calc_baro_u

! -----------------------------------------------------------------------------

!> Calculate the barotropic v velocity
subroutine calc_baro_v(vb, v, h, eta, freesurfFac, nx, ny, layers)
  implicit none

  integer nx, ny, layers
  integer i, j, k
  double precision h(0:nx+1, 0:ny+1, layers)
  double precision h_temp(0:nx+1, 0:ny+1, layers)
  double precision eta(0:nx+1, 0:ny+1)
  double precision freesurfFac
  double precision v(0:nx+1, 0:ny+1, layers)
  double precision vb(0:nx+1, 0:ny+1)

  vb = 0d0

  h_temp = h
  ! add free surface elevation to the upper layer
  h_temp(:, :, 1) = h(:, :, 1) + eta*freesurfFac

  do i = 1, nx
    do j = 1, ny
      do k = 1, layers
        vb(i,j) = vb(i,j) + v(i,j,k)*(h_temp(i,j,k)+h_temp(i,j-1,k))/2d0
      end do
    end do
  end do

  call wrap_fields_2D(vb, nx, ny)

  return
end subroutine calc_baro_v

! -----------------------------------------------------------------------------

!> Calculate the free surface anomaly using the velocities
!! timestepped with the tendencies excluding the free surface
!! pressure gradient.
subroutine calc_eta_star(ub, vb, eta, etastar, freesurfFac, nx, ny, dx, dy, dt)
  implicit none

  integer nx, ny
  integer i, j
  double precision eta(0:nx+1, 0:ny+1)
  double precision etastar(0:nx+1, 0:ny+1)
  double precision ub(0:nx+1, 0:ny+1)
  double precision vb(0:nx+1, 0:ny+1)
  double precision freesurfFac, dx, dy, dt

  etastar = 0d0

  do i = 1, nx
    do j = 1, ny
      etastar(i,j) = freesurfFac*eta(i,j) - &
          dt*((ub(i+1,j) - ub(i,j))/dx + (vb(i,j+1) - vb(i,j))/dy)
    end do
  end do

  call wrap_fields_2D(etastar, nx, ny)

  return
end subroutine calc_eta_star

! -----------------------------------------------------------------------------
!> Use the successive over-relaxation algorithm to solve the backwards
!! Euler timestepping for the free surface anomaly, or for the surface
!! pressure required to keep the barotropic flow nondivergent.

subroutine SOR_solver(a, etanew, etastar, freesurfFac, nx, ny, dt, &
    rjac, eps, maxits, n)
  implicit none

  double precision a(5, nx, ny)
  double precision etanew(0:nx+1, 0:ny+1)
  double precision etastar(0:nx+1, 0:ny+1)
  double precision freesurfFac
  integer nx, ny, i, j, maxits, n, nit
  double precision dt
  double precision rjac, eps
  double precision rhs(nx, ny)
  double precision res(nx, ny)
  double precision norm, norm0
  double precision relax_param

  rhs = -etastar(1:nx,1:ny)/dt**2
  ! first guess for etanew
  etanew = etastar

  relax_param = 1.d0 ! successive over-relaxation parameter

  ! Calculate initial residual, so that we can stop the loop when the
  ! current residual = norm0*eps
  norm0 = 0.d0
  do i = 1, nx
    do j = 1, ny
      res(i,j) = &
          a(1,i,j)*etanew(i+1,j) &
          + a(2,i,j)*etanew(i,j+1) &
          + a(3,i,j)*etanew(i-1,j) &
          + a(4,i,j)*etanew(i,j-1) &
          + a(5,i,j)*etanew(i,j)   &
          - freesurfFac*etanew(i,j)/dt**2 &
          - rhs(i,j)
      norm0 = norm0 + abs(res(i,j))
      etanew(i,j) = etanew(i,j)-relax_param*res(i,j)/a(5,i,j)
    end do
  end do


  do nit = 1, maxits
    norm = 0.d0
    do i = 1, nx
      do j = 1, ny
        res(i,j) = &
            a(1,i,j)*etanew(i+1,j) &
            + a(2,i,j)*etanew(i,j+1) &
            + a(3,i,j)*etanew(i-1,j) &
            + a(4,i,j)*etanew(i,j-1) &
            + a(5,i,j)*etanew(i,j)   &
            - freesurfFac*etanew(i,j)/dt**2 &
            - rhs(i,j)
        norm = norm + abs(res(i,j))
        etanew(i,j) = etanew(i,j)-relax_param*res(i,j)/(a(5,i,j))
      end do
    end do
    if (nit.eq.1) then
      relax_param = 1.d0/(1.d0-0.5d0*rjac**2)
    else
      relax_param = 1.d0/(1.d0-0.25d0*rjac**2*relax_param)
    end if

    call wrap_fields_2D(etanew, nx, ny)

    if (nit.gt.1.and.norm.lt.eps*norm0) then

      return

    end if
  end do

  print *, 'warning: maximum iterations exceeded at time step ', n

  return
end subroutine SOR_solver

! -----------------------------------------------------------------------------
!> Check to see if there are any NaNs in the data field and stop the
!! calculation if any are found.
subroutine break_if_NaN(data, nx, ny, layers, n)
  implicit none

  ! To stop the program if it detects at NaN in the variable being checked

  integer nx, ny, layers, n, i, j, k
  double precision data(0:nx+1, 0:ny+1, layers)

  do k = 1, layers
    do j = 1, ny
      do i = 1, nx
        if (data(i,j,k) .ne. data(i,j,k)) then
          ! write a file saying so
          open(unit=10, file='NaN detected.txt', action="write", &
              status="replace", form="formatted")
          write(10, 1000) n
1000      format( "NaN detected at time step ", 1i10.10)
          close(unit=10)
          ! print it on the screen
          print *, 'NaN detected'
          ! Stop the code
          stop 'Nan detected'
        end if
      end do
    end do
  end do

  return
end subroutine break_if_NaN

!------------------------------------------------------------------------------
!> Define masks for boundary conditions in u and v.
!! This finds locations where neighbouring grid boxes are not the same
!! (i.e. one is land and one is ocean).
!! In the output,
!! 0 means barrier
!! 1 mean open

subroutine calc_boundary_masks(wetmask, hfacW, hfacE, hfacS, hfacN, nx, ny)
  implicit none

  integer nx !< number of grid points in x direction
  integer ny !< number of grid points in y direction
  double precision wetmask(0:nx+1, 0:ny+1)
  double precision temp(0:nx+1, 0:ny+1)
  double precision hfacW(0:nx+1, 0:ny+1)
  double precision hfacE(0:nx+1, 0:ny+1)
  double precision hfacN(0:nx+1, 0:ny+1)
  double precision hfacS(0:nx+1, 0:ny+1)
  integer i, j

  hfacW = 1d0

  temp = 0.0
  do j = 0, ny+1
    do i = 1, nx+1
      temp(i, j) = wetmask(i-1, j) - wetmask(i, j)
    end do
  end do

  do j = 0, ny+1
    do i = 1, nx+1
      if (temp(i, j) .ne. 0.0) then
        hfacW(i, j) = 0d0
      end if
    end do
  end do

  ! and now for all  western cells
  hfacW(0, :) = hfacW(nx, :)

  hfacE = 1d0

  temp = 0.0
  do j = 0, ny+1
    do i = 0, nx
      temp(i, j) = wetmask(i, j) - wetmask(i+1, j)
    end do
  end do

  do j = 0, ny+1
    do i = 0, nx
      if (temp(i, j) .ne. 0.0) then
        hfacE(i, j) = 0d0
      end if
    end do
  end do

  ! and now for all  eastern cells
  hfacE(nx+1, :) = hfacE(1, :)

  hfacS = 1

  temp = 0.0
  do j = 1, ny+1
    do i = 0, nx+1
      temp(i, j) = wetmask(i, j-1) - wetmask(i, j)
    end do
  end do

  do j = 1, ny+1
    do i = 0, nx+1
      if (temp(i, j) .ne. 0.0) then
        hfacS(i, j) = 0d0
      end if
    end do
  end do

  ! all southern cells
  hfacS(:, 0) = hfacS(:, ny)

  hfacN = 1
  temp = 0.0
  do j = 0, ny
    do i = 0, nx+1
      temp(i, j) = wetmask(i, j) - wetmask(i, j+1)
    end do
  end do

  do j = 0, ny
    do i = 0, nx+1
      if (temp(i, j) .ne. 0.0) then
        hfacN(i, j) = 0d0
      end if
    end do
  end do
  ! all northern cells
  hfacN(:, ny+1) = hfacN(:, 1)

  return
end subroutine calc_boundary_masks

! -----------------------------------------------------------------------------

subroutine read_input_fileH(name, array, default, nx, ny, layers)
  implicit none

  character(30) name
  integer nx, ny, layers, k
  double precision array(0:nx+1, 0:ny+1, layers), default(layers)
  double precision array_small(nx, ny, layers)

  if (name.ne.'') then
    open(unit=10, form='unformatted', file=name)
    read(10) array_small
    close(10)

    array(1:nx, 1:ny, :) = array_small
    ! wrap array around for periodicity
    array(0, :, :) = array(nx, :, :)
    array(nx+1, :, :) = array(1, :, :)
    array(:, 0, :) = array(:, ny, :)
    array(:, ny+1, :) = array(:, 1, :)
  else
    do k = 1, layers
      array(:, :, k) = default(k)
    end do
  end if

  return
end subroutine read_input_fileH

! -----------------------------------------------------------------------------

subroutine read_input_fileH_2D(name, array, default, nx, ny)
  implicit none

  character(30) name
  integer nx, ny
  double precision array(0:nx+1, 0:ny+1), default
  double precision array_small(nx, ny)

  if (name.ne.'') then
    open(unit=10, form='unformatted', file=name)
    read(10) array_small
    close(10)

    array(1:nx, 1:ny) = array_small
    ! wrap array around for periodicity
    array(0, :) = array(nx, :)
    array(nx+1, :) = array(1, :)
    array(:, 0) = array(:, ny)
    array(:, ny+1) = array(:, 1)
  else
    array = default
  end if

  return
end subroutine read_input_fileH_2D

! -----------------------------------------------------------------------------

subroutine read_input_fileU(name, array, default, nx, ny, layers)
  implicit none

  character(30) name
  integer nx, ny, layers
  double precision array(0:nx+1, 0:ny+1, layers), default
  double precision array_small(nx+1, ny, layers)

  if (name.ne.'') then
    open(unit=10, form='unformatted', file=name)
    read(10) array_small
    close(10)

    array(1:nx+1, 1:ny, :) = array_small
    ! wrap array around for periodicity
    array(0, :, :) = array(nx, :, :)
    array(nx+1, :, :) = array(1, :, :)
    array(:, 0, :) = array(:, ny, :)
    array(:, ny+1, :) = array(:, 1, :)
  else
    array = default
  end if

  return
end subroutine read_input_fileU

! -----------------------------------------------------------------------------

subroutine read_input_fileV(name, array, default, nx, ny, layers)
  implicit none

  character(30) name
  integer nx, ny, layers
  double precision array(0:nx+1, 0:ny+1, layers), default
  double precision array_small(nx, ny+1, layers)

  if (name.ne.'') then
    open(unit=10, form='unformatted', file=name)
    read(10) array_small
    close(10)

    array(1:nx, 1:ny+1, :) = array_small
    ! wrap array around for periodicity
    array(0, :, :) = array(nx, :, :)
    array(nx+1, :, :) = array(1, :, :)
    array(:, 0, :) = array(:, ny, :)
    array(:, ny+1, :) = array(:, 1, :)
  else
    array = default
  end if

  return
end subroutine read_input_fileV

! -----------------------------------------------------------------------------

subroutine read_input_file_time_series(name, array, default, nTimeSteps)
  implicit none

  character(30) name
  integer nTimeSteps
  double precision array(nTimeSteps), default

  if (name.ne.'') then
    open(unit=10, form='unformatted', file=name)
    read(10) array
    close(10)

  else
    array = default
  end if

  return
end subroutine read_input_file_time_series

!-----------------------------------------------------------------
!> Wrap 3D fields around for periodic boundary conditions

subroutine wrap_fields_3D(array, nx, ny, layers)
  implicit none

  double precision array(0:nx+1, 0:ny+1, layers)
  integer nx, ny, layers

  ! wrap array around for periodicity
  array(0, :, :) = array(nx, :, :)
  array(nx+1, :, :) = array(1, :, :)
  array(:, 0, :) = array(:, ny, :)
  array(:, ny+1, :) = array(:, 1, :)

  return
end subroutine wrap_fields_3D

!-----------------------------------------------------------------
!> Wrap 2D fields around for periodic boundary conditions

subroutine wrap_fields_2D(array, nx, ny)
  implicit none

  double precision array(0:nx+1, 0:ny+1)
  integer nx, ny

  ! wrap array around for periodicity
  array(0, :) = array(nx, :)
  array(nx+1, :) = array(1, :)
  array(:, 0) = array(:, ny)
  array(:, ny+1) = array(:, 1)

  return
end subroutine wrap_fields_2D

!-----------------------------------------------------------------
!> Robustly produce a different random seed for each run.
!! Code adapted from
!! http://web.ph.surrey.ac.uk/fortweb/glossary/random_seed.html

subroutine ranseed()
  implicit none

  ! ----- variables for portable seed setting -----
  integer :: i_seed
  integer, dimension(:), allocatable :: a_seed
  integer, dimension(1:8) :: dt_seed
  ! ----- end of variables for seed setting -----

  ! ----- Set up random seed portably -----
  call random_seed(size=i_seed)
  allocate(a_seed(1:i_seed))
  call random_seed(get=a_seed)
  call date_and_time(values=dt_seed)
  a_seed(i_seed) = dt_seed(8); a_seed(1) = dt_seed(8)*dt_seed(7)*dt_seed(6)
  call random_seed(put=a_seed)
  return
end subroutine ranseed
