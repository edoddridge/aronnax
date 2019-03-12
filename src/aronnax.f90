!> @author
!> Ed Doddridge
!
!> Aronnax, an idealized isopycnal model with n layers and variable bathymetry.
!!
!
!>
!>     @mainpage Documentation for aronnax.f90
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


program aronnax

  use model_main
  use declarations
  use mpi

  implicit none

  ! TODO Possibly wait until the model is split into multiple files,
  ! then hide the long unsightly code there.

  namelist /NUMERICS/ au, kh, kv, ar, botDrag, dt, slip, &
      niter0, nTimeSteps, hAdvecScheme, TS_algorithm, &
      dumpFreq, avFreq, checkpointFreq, diagFreq, hmin, maxits, &
      freesurfFac, eps, thickness_error, debug_level

  namelist /MODEL/ hmean, depthFile, H0, RedGrav

  namelist /PRESSURE_SOLVER/ nProcX, nProcY

  namelist /SPONGE/ spongeHTimeScaleFile, spongeUTimeScaleFile, &
      spongeVTimeScaleFile, spongeHfile, spongeUfile, spongeVfile

  namelist /PHYSICS/ g_vec, rho0

  namelist /GRID/ nx, ny, layers, OL, dx, dy, fUfile, fVfile, wetMaskFile

  namelist /INITIAL_CONDITIONS/ initUfile, initVfile, initHfile, initEtaFile

  namelist /EXTERNAL_FORCING/ zonalWindFile, meridionalWindFile, &
      RelativeWind, Cd, &
      DumpWind, wind_mag_time_series_file, wind_depth

  start_time = time()

  ! Set default values here

  ! io default frequencies
  dumpFreq = 1d9
  avFreq = 0d0
  checkpointFreq = 0d0
  diagFreq = 0d0

  debug_level = 0
  
  ! start from t = 0
  niter0 = 0

  ! use N/m^2
  RelativeWind = .FALSE.

  ! wind forcing only in the top layer
  wind_depth = 0d0

  ! use first-order centred differencing
  hAdvecScheme = 1

  ! use third-order AB time stepping
  TS_algorithm = 3

  ! No viscosity or diffusion
  au = 0d0
  ar = 0d0
  kh = 0d0
  kv = 0d0

  ! size of halo region
  OL = 1

  open(unit=8, file="parameters.in", status='OLD', recl=80)
  read(unit=8, nml=NUMERICS)
  read(unit=8, nml=MODEL)
  read(unit=8, nml=PRESSURE_SOLVER)
  read(unit=8, nml=SPONGE)
  read(unit=8, nml=PHYSICS)
  read(unit=8, nml=GRID)
  read(unit=8, nml=INITIAL_CONDITIONS)
  read(unit=8, nml=EXTERNAL_FORCING)
  close(unit=8)


  ! set timestepping order for linear multi-step methods
  ! based on TS_algorithm
  call set_AB_order(TS_algorithm, AB_order)


  ! optionally include the MPI code for parallel runs with external
  ! pressure solver
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  ! mpi_comm = MPI_COMM_WORLD

  if (num_procs .ne. nProcX * nProcY) then
    if (myid .eq. 0) then
       write(17, "(A)") "number of processors in run command must equal nProcX * nProcY - fix this and try again"
       write(17, "(A, I0)") 'num_procs = ', num_procs
       write(17, "(A, I0)") 'nProcX = ', nProcX
       write(17, "(A, I0)") 'nProcY = ', nProcY
    end if
    call clean_stop(0, .FALSE.)
  end if

  ! myid starts at zero, so index these variables from zero.
  ! i__(:,1) = indicies for x locations
  ! i__(:,2) = indicies for y locations
  allocate(ilower(0:num_procs-1, 2))
  allocate(iupper(0:num_procs-1, 2))

  allocate(xlower(0:num_procs-1))
  allocate(xupper(0:num_procs-1))
  allocate(ylower(0:num_procs-1))
  allocate(yupper(0:num_procs-1))

  ! calculate x and y limits for each tile
  do i = 0,nProcX-1
    xlower(i) = (i*nx/nProcX) + 1
    xupper(i) = (i+1)*nx/nProcX
  end do

  do i = 0,nProcY-1
    ylower(i) = (i*ny/nProcY) + 1
    yupper(i) = (i+1)*ny/nProcY
  end do

  ! Combine into ilower and iupper vectors for hypre
  do k = 0, num_procs-1
    ! i,j as horizontal idicies for tiles
    i = mod(k, nProcX)
    j = floor(real(k)/nProcX)

    ! x locations
    ilower(k,1) = xlower(i)
    iupper(k,1) = xupper(i)

    ! y ocations
    ilower(k,2) = ylower(j)
    iupper(k,2) = yupper(j)

  end do


#ifdef useExtSolver
  call create_Hypre_grid(MPI_COMM_WORLD, hypre_grid, ilower, iupper, &
          num_procs, myid, nx, ny, ierr)
#endif

  allocate(h(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(u(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(v(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(eta(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(depth(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))

  allocate(wetmask(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(hfacW(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(hfacE(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(hfacS(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(hfacN(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))

  allocate(fu(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(fv(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))

  allocate(zeros(layers))

  allocate(base_wind_x(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(base_wind_y(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(wind_mag_time_series(nTimeSteps))

  allocate(spongeHTimeScale(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(spongeUTimeScale(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(spongeVTimeScale(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(spongeH(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(spongeU(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(spongeV(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))

  ! Zero vector - for internal use only
  zeros = 0d0


  ! Read in arrays from the input files
  call read_input_fileU(initUfile, u, 0.d0, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileV(initVfile, v, 0.d0, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileH(initHfile, h, hmean, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)

  call read_input_fileU(fUfile, fu, 0.d0, nx, ny, 1, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileV(fVfile, fv, 0.d0, nx, ny, 1, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)

  call read_input_fileU(zonalWindFile, base_wind_x, 0.d0, nx, ny, 1, &
                        OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileV(meridionalWindFile, base_wind_y, 0.d0, nx, ny, &
                        1, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)

  call read_input_file_time_series(wind_mag_time_series_file, &
      wind_mag_time_series, 1d0, nTimeSteps)

  call read_input_fileH(spongeHTimeScaleFile, spongeHTimeScale, &
      zeros, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileH(spongeHfile, spongeH, hmean, nx, ny, layers, &
                        OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileU(spongeUTimeScaleFile, spongeUTimeScale, &
      0.d0, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileU(spongeUfile, spongeU, 0.d0, nx, ny, layers, &
                        OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileV(spongeVTimeScaleFile, spongeVTimeScale, &
      0.d0, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileV(spongeVfile, spongeV, 0.d0, nx, ny, layers, &
                        OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)

  call read_input_fileH_2D(wetMaskFile, wetmask, 1d0, nx, ny, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)

  call calc_boundary_masks(wetMaskFile, hfacW, hfacE, hfacS, hfacN, nx, ny, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2))


  if (.not. RedGrav) then
    call read_input_fileH_2D(depthFile, depth, H0, nx, ny, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
    call read_input_fileH_2D(initEtaFile, eta, 0.d0, nx, ny, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
    ! Check that depth is positive - it must be greater than zero
    if (minval(depth) .lt. 0) then
      write(17, "(A)") "Depths must be positive."
      call clean_stop(0, .FALSE.)
    end if
  end if


  call model_run(h, u, v, eta, depth, dx, dy, wetmask, &
      hfacW, hfacE, hfacS, hfacN, &
      fu, fv, &
      dt, au, ar, botDrag, kh, kv, slip, hmin, niter0, nTimeSteps, &
      dumpFreq, avFreq, checkpointFreq, diagFreq, &
      maxits, eps, freesurfFac, thickness_error, &
      debug_level, g_vec, rho0, &
      base_wind_x, base_wind_y, wind_mag_time_series, wind_depth, &
      spongeHTimeScale, spongeUTimeScale, spongeVTimeScale, &
      spongeH, spongeU, spongeV, &
      nx, ny, layers, OL, &
      ilower(myid,1), iupper(myid,1), &
      ilower(myid,2), iupper(myid,2), &
      RedGrav, hAdvecScheme, TS_algorithm, AB_order, &
      DumpWind, RelativeWind, Cd, start_time, &
      MPI_COMM_WORLD, myid, num_procs, ilower, iupper, &
      hypre_grid)

  ! Finalize MPI
  call clean_stop(nTimeSteps, .TRUE.)
end program aronnax
