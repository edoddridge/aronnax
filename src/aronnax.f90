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

  namelist /NUMERICS/ au, kh, kv, ar, bot_drag, dt, slip, &
      niter0, n_time_steps, h_advec_scheme, ts_algorithm, &
      dump_freq, av_freq, checkpoint_freq, diag_freq, hmin, maxits, &
      freesurf_fac, eps, thickness_error, debug_level

  namelist /MODEL/ hmean, depth_file, h0, red_grav, active_lower_layer

  namelist /PRESSURE_SOLVER/ n_proc_x, n_proc_y

  namelist /SPONGE/ sponge_h_time_scale_file, sponge_u_time_scale_file, &
      sponge_v_time_scale_file, sponge_h_file, sponge_u_file, sponge_v_file

  namelist /PHYSICS/ g_vec, rho0

  namelist /GRID/ nx, ny, layers, OL, dx, dy, f_u_file, f_v_file, wet_mask_file

  namelist /INITIAL_CONDITIONS/ init_u_file, init_v_file, init_h_file, init_eta_file

  namelist /EXTERNAL_FORCING/ zonal_wind_file, meridional_wind_file, &
      relative_wind, Cd, dump_wind, wind_mag_time_series_file, wind_depth, &
      wind_n_records, wind_period, wind_loop_fields, wind_interpolate

  start_time = time()

  ! Set default values here

  ! io default frequencies
  dump_freq = 1d9
  av_freq = 0d0
  checkpoint_freq = 0d0
  diag_freq = 0d0

  debug_level = 0
  
  ! start from t = 0
  niter0 = 0

  ! use N/m^2
  relative_wind = .FALSE.

  ! wind forcing only in the top layer
  wind_depth = 0d0

  ! use first-order centred differencing
  h_advec_scheme = 1

  ! use third-order AB time stepping
  ts_algorithm = 3

  ! No viscosity or diffusion
  au = 0d0
  ar = 0d0
  kh = 0d0
  kv = 0d0

  ! size of halo region
  OL = 3

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
  ! based on ts_algorithm
  call set_AB_order(ts_algorithm, AB_order)


  ! optionally include the MPI code for parallel runs with external
  ! pressure solver
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  ! mpi_comm = MPI_COMM_WORLD

  if (num_procs .ne. n_proc_x * n_proc_y) then
    if (myid .eq. 0) then
       write(17, "(A)") "number of processors in run command must equal n_proc_x * n_proc_y - fix this and try again"
       write(17, "(A, I0)") 'num_procs = ', num_procs
       write(17, "(A, I0)") 'n_proc_x = ', n_proc_x
       write(17, "(A, I0)") 'n_proc_y = ', n_proc_y
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
  do i = 0,n_proc_x-1
    xlower(i) = (i*nx/n_proc_x) + 1
    xupper(i) = (i+1)*nx/n_proc_x
  end do

  do i = 0,n_proc_y-1
    ylower(i) = (i*ny/n_proc_y) + 1
    yupper(i) = (i+1)*ny/n_proc_y
  end do

  ! Combine into ilower and iupper vectors for hypre
  do k = 0, num_procs-1
    ! i,j as horizontal idicies for tiles
    i = mod(k, n_proc_x)
    j = floor(real(k)/n_proc_x)

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
  allocate(hfac_w(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(hfac_e(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(hfac_s(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(hfac_n(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))

  allocate(fu(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))
  allocate(fv(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL))

  allocate(zeros(layers))

  allocate(base_wind_x(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, 0:wind_n_records-1))
  allocate(base_wind_y(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, 0:wind_n_records-1))
  allocate(wind_mag_time_series(n_time_steps))

  allocate(sponge_h_time_scale(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(sponge_u_time_scale(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(sponge_v_time_scale(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(sponge_h(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(sponge_u(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))
  allocate(sponge_v(ilower(myid,1)-OL:iupper(myid,1)+OL, &
             ilower(myid,2)-OL:iupper(myid,2)+OL, layers))

  ! Zero vector - for internal use only
  zeros = 0d0


  ! Read in arrays from the input files
  call read_input_fileU(init_u_file, u, 0.d0, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileV(init_v_file, v, 0.d0, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileH(init_h_file, h, hmean, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)

  call read_input_fileU(f_u_file, fu, 0.d0, nx, ny, 1, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileV(f_v_file, fv, 0.d0, nx, ny, 1, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)

  call read_forcing_fileU(zonal_wind_file, base_wind_x, 0.d0, &
                              wind_n_records, nx, ny, 1, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_forcing_fileV(meridional_wind_file, base_wind_y, 0.d0, &
                              wind_n_records, nx, ny, 1, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)

  call read_input_file_time_series(wind_mag_time_series_file, &
      wind_mag_time_series, 1d0, n_time_steps)

  call read_input_fileH(sponge_h_time_scale_file, sponge_h_time_scale, &
      zeros, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileH(sponge_h_file, sponge_h, hmean, nx, ny, layers, &
                        OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileU(sponge_u_time_scale_file, sponge_u_time_scale, &
      0.d0, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileU(sponge_u_file, sponge_u, 0.d0, nx, ny, layers, &
                        OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileV(sponge_v_time_scale_file, sponge_v_time_scale, &
      0.d0, nx, ny, layers, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
  call read_input_fileV(sponge_v_file, sponge_v, 0.d0, nx, ny, layers, &
                        OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)

  call read_input_fileH_2D(wet_mask_file, wetmask, 1d0, nx, ny, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)

  call calc_boundary_masks(wet_mask_file, hfac_w, hfac_e, hfac_s, hfac_n, nx, ny, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2))


  if (.not. red_grav) then
    call read_input_fileH_2D(depth_file, depth, h0, nx, ny, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
    call read_input_fileH_2D(init_eta_file, eta, 0.d0, nx, ny, OL, &
                              ilower(myid,1), iupper(myid,1), &
                              ilower(myid,2), iupper(myid,2), myid)
    ! Check that depth is positive - it must be greater than zero
    if (minval(depth) .lt. 0) then
      write(17, "(A)") "Depths must be positive."
      call clean_stop(0, .FALSE.)
    end if
  end if


  call model_run(h, u, v, eta, depth, dx, dy, wetmask, &
      hfac_w, hfac_e, hfac_s, hfac_n, &
      fu, fv, &
      dt, au, ar, bot_drag, kh, kv, slip, hmin, niter0, n_time_steps, &
      dump_freq, av_freq, checkpoint_freq, diag_freq, &
      maxits, eps, freesurf_fac, thickness_error, &
      debug_level, g_vec, rho0, &
      base_wind_x, base_wind_y, wind_mag_time_series, &
      wind_n_records, wind_period, wind_loop_fields, wind_interpolate, &
      wind_depth, &
      sponge_h_time_scale, sponge_u_time_scale, sponge_v_time_scale, &
      sponge_h, sponge_u, sponge_v, &
      nx, ny, layers, OL, &
      ilower(myid,1), iupper(myid,1), &
      ilower(myid,2), iupper(myid,2), &
      red_grav, h_advec_scheme, ts_algorithm, AB_order, &
      dump_wind, relative_wind, Cd, start_time, &
      MPI_COMM_WORLD, myid, num_procs, ilower, iupper, &
      hypre_grid)

  ! Finalize MPI
  call clean_stop(n_time_steps, .TRUE.)
end program aronnax
