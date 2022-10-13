module model_main
  use end_run
  use boundaries
  use barotropic_mode
  use state_deriv
  use time_stepping
  use enforce_thickness
  use exchange

  implicit none

  contains

  ! ------------------------------ Primary routine ----------------------------
  !> Run the model

  subroutine model_run(h, u, v, eta, depth, dx, dy, wetmask, &
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
      nx, ny, layers, OL, xlow, xhigh, ylow, yhigh, &
      red_grav, h_advec_scheme, ts_algorithm, AB_order, &
      dump_wind, relative_wind, Cd, start_time, &
      MPI_COMM_WORLD, myid, num_procs, ilower, iupper, &
      hypre_grid)
    implicit none

    ! Layer thickness (h)
    double precision, intent(inout) :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    ! Velocity component (u)
    double precision, intent(inout) :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    ! Velocity component (v)
    double precision, intent(inout) :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    ! Free surface (eta)
    double precision, intent(inout) :: eta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    ! Bathymetry
    double precision, intent(in) :: depth(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    ! Grid
    double precision, intent(in) :: dx, dy
    double precision, intent(in) :: wetmask(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfac_w(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfac_e(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfac_s(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfac_n(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    ! Coriolis parameter at u and v grid-points respectively
    double precision, intent(in) :: fu(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: fv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    ! Numerics
    double precision, intent(in) :: dt, au, ar, bot_drag
    double precision, intent(in) :: kh(layers), kv
    double precision, intent(in) :: slip, hmin
    integer,          intent(in) :: niter0, n_time_steps
    double precision, intent(in) :: dump_freq, av_freq, checkpoint_freq, diag_freq
    integer,          intent(in) :: maxits
    double precision, intent(in) :: eps, freesurf_fac, thickness_error
    integer,          intent(in) :: debug_level
    ! Physics
    double precision, intent(in) :: g_vec(layers)
    double precision, intent(in) :: rho0
    ! Wind
    double precision, intent(in) :: base_wind_x(xlow-OL:xhigh+OL, &
                                                ylow-OL:yhigh+OL, 0:wind_n_records-1)
    double precision, intent(in) :: base_wind_y(xlow-OL:xhigh+OL, &
                                                ylow-OL:yhigh+OL, 0:wind_n_records-1)
    double precision, intent(in) :: wind_mag_time_series(n_time_steps)
    integer         , intent(in) :: wind_n_records
    double precision, intent(in) :: wind_period
    logical         , intent(in) :: wind_loop_fields
    logical         , intent(in) :: wind_interpolate
    double precision, intent(in) :: wind_depth
    ! Sponge regions
    double precision, intent(in) :: sponge_h_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_u_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_v_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    ! Resolution
    integer,          intent(in) :: nx, ny, layers, OL
    integer,          intent(in) :: xlow, xhigh, ylow, yhigh
    ! Reduced gravity vs n-layer physics
    logical,          intent(in) :: red_grav
    integer,          intent(in) :: h_advec_scheme
    integer,          intent(in) :: ts_algorithm
    integer,          intent(in) :: AB_order
    ! Whether to write computed wind in the output
    logical,          intent(in) :: dump_wind
    logical,          intent(in) :: relative_wind
    double precision, intent(in) :: Cd
    integer*8,        intent(inout) :: start_time

    double precision :: dhdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision :: h_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    ! for saving average fields
    double precision :: hav(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)

    double precision :: dudt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision :: u_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    ! for saving average fields
    double precision :: uav(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)

    double precision :: dvdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision :: v_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    ! for saving average fields
    double precision :: vav(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)

    double precision :: etanew(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    ! for saving average fields
    double precision :: etaav(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)

    ! Pressure solver variables
    double precision :: a(5, xlow:xhigh, ylow:yhigh)

    ! Numerics
    double precision :: pi
    double precision :: rjac

    ! dumping output
    integer :: nwrite, avwrite, checkpointwrite, diagwrite

    ! External solver variables
    integer          :: offsets(2,5)
    integer          :: i, j ! loop variables
    double precision :: values((xhigh - xlow + 1)*(yhigh - ylow + 1))
    integer          :: indicies(2)
    integer*8        :: hypre_grid
    integer*8        :: stencil
    integer*8        :: hypre_A
    integer          :: ilower(0:num_procs-1,2), iupper(0:num_procs-1,2)
    integer          :: ierr
    integer          :: MPI_COMM_WORLD
    integer          :: myid
    integer          :: num_procs

    ! Time step loop variable
    integer :: n

    ! Wind
    double precision :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    integer          :: wind_n, wind_np1
    double precision :: wind_n_remainder

    ! Time
    integer*8 :: last_report_time, cur_time


    if (myid .eq. 0) then
      if (red_grav) then
        print "(A, I0, A, I0, A, I0, A, I0, A)", &
            "Running a reduced-gravity configuration of size ", &
            nx, "x", ny, "x", layers, " by ", n_time_steps, " time steps."
      else
        print "(A, I0, A, I0, A, I0, A, I0, A)", &
            "Running an n-layer configuration of size ", &
            nx, "x", ny, "x", layers, " by ", n_time_steps, " time steps."
      end if

    !   ! Show the domain decomposition
      do n = 0, num_procs-1
        print "(A, I0, A, I0, A, I0)", 'tile ', n, ': x limits = ', &
                  ilower(n,1), ", ", iupper(n,1)
        print "(A, I0, A, I0, A, I0)", 'tile ', n, ': y limits = ', &
                  ilower(n,2), ", ", iupper(n,2)
      end do
    end if

    last_report_time = start_time

    nwrite = int(dump_freq/dt)
    avwrite = int(av_freq/dt)
    checkpointwrite = int(checkpoint_freq/dt)
    diagwrite = int(diag_freq/dt)

    ! Pi, the constant
    pi = 3.1415926535897932384

    ! Initialize wind fields
    wind_x = base_wind_x(:,:,0)*wind_mag_time_series(1)
    wind_y = base_wind_y(:,:,0)*wind_mag_time_series(1)

    if (myid .eq. 0) then
      ! Initialise the diagnostic files
      call create_diag_file(layers, 'output/diagnostic.h.csv', 'h', niter0)
      call create_diag_file(layers, 'output/diagnostic.u.csv', 'u', niter0)
      call create_diag_file(layers, 'output/diagnostic.v.csv', 'v', niter0)
      if (.not. red_grav) then
        call create_diag_file(1, 'output/diagnostic.eta.csv', 'eta', niter0)
      end if
    end if

    ! Initialise the average fields
    if (avwrite .ne. 0) then
      hav = 0.0
      uav = 0.0
      vav = 0.0
      if (.not. red_grav) then
        etaav = 0.0
      end if
    end if

    ! initialise etanew
    etanew = 0d0

    call apply_boundary_conditions(u, hfac_w, wetmask, &
                          xlow, xhigh, ylow, yhigh, layers, OL)
    call apply_boundary_conditions(v, hfac_s, wetmask, &
                          xlow, xhigh, ylow, yhigh, layers, OL)


    if (.not. red_grav) then
      ! Initialise arrays for pressure solver
      ! a = derivatives of the depth field
        call calc_A_matrix(a, depth, g_vec(1), dx, dy, OL, &
            xlow, xhigh, ylow, yhigh, freesurf_fac, dt, &
            hfac_w, hfac_e, hfac_s, hfac_n)

#ifndef useExtSolver
      ! Calculate the spectral radius of the grid for use by the
      ! successive over-relaxation scheme
      rjac = (cos(pi/real(nx))*dy**2+cos(pi/real(ny))*dx**2) &
             /(dx**2+dy**2)
      ! If peridodic boundary conditions are ever implemented, then pi ->
      ! 2*pi in this calculation
#else
      ! use the external pressure solver
      call create_Hypre_A_matrix(MPI_COMM_WORLD, hypre_grid, hypre_A, &
            a, xlow, xhigh, ylow, yhigh, ierr)
#endif

      ! Check that the supplied free surface anomaly and layer
      ! thicknesses are consistent with the supplied depth field.
      ! If they are not, then scale the layer thicknesses to make
      ! them consistent.
      call enforce_depth_thickness_consistency(h, eta, depth, &
          freesurf_fac, thickness_error, &
          xlow, xhigh, ylow, yhigh, layers, OL)
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  Initialisation of the model STARTS HERE                            !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (niter0 .eq. 0) then

      n = 0
      call initialise_tendencies(dhdt, dudt, dvdt, h, u, v, depth, &
          dx, dy, dt, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, AB_order, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, debug_level)

    else if (niter0 .ne. 0) then
      n = niter0

      call load_checkpoint_files(dhdt, dudt, dvdt, h, u, v, eta, &
            red_grav, niter0, nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, AB_order)

    end if

    ! Now the model is ready to start.
    ! - We have h, u, v at the zeroth time step, and the tendencies at
    !   two older time steps.
    ! - The model then solves for the tendencies at the current step
    !   before solving for the fields at the next time step.

    cur_time = time()
    if (myid .eq. 0) then
      if (cur_time - start_time .eq. 1) then
        print "(A)", "Initialized in 1 second."
      else
        print "(A, I0, A)", "Initialized in " , cur_time - start_time, " seconds."
      end if
    end if

    last_report_time = cur_time

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! MAIN LOOP OF THE MODEL STARTS HERE                                  !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n = niter0+1, niter0+n_time_steps

      ! calculate wind_n
      if (wind_period .eq. 0) then
        wind_n = 0
        wind_np1 = 0
        wind_n_remainder = 0
      else
        wind_n = FLOOR(n*dt/wind_period)
        wind_np1 = wind_n + 1
        wind_n_remainder = (n*dt/wind_period) - wind_n
      end if

      ! prevent wind_n from being too big
      if (wind_loop_fields) then
        wind_n = MOD(wind_n, wind_n_records)
        wind_np1 = MOD(wind_np1, wind_n_records)
      else
        wind_n = MIN(wind_n, wind_n_records-1)
        wind_np1 = MIN(wind_np1, wind_n_records-1)
      end if

      if (wind_interpolate) then
        wind_x = ((1d0 - wind_n_remainder)*base_wind_x(:,:,wind_n) + &
                    wind_n_remainder*base_wind_x(:,:,wind_np1))* &
                    wind_mag_time_series(n-niter0)
        wind_y = ((1d0 - wind_n_remainder)*base_wind_y(:,:,wind_n) + &
                    wind_n_remainder*base_wind_y(:,:,wind_np1))* &
                    wind_mag_time_series(n-niter0)
      else
        wind_x = base_wind_x(:,:,wind_n)*wind_mag_time_series(n-niter0)
        wind_y = base_wind_y(:,:,wind_n)*wind_mag_time_series(n-niter0)
      end if

      call timestep(h_new, u_new, v_new, dhdt, dudt, dvdt, &
          h, u, v, depth, &
          dx, dy, dt, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, ts_algorithm, AB_order, &
          g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Apply the boundary conditions
      call apply_boundary_conditions(u_new, hfac_w, wetmask, &
                              xlow, xhigh, ylow, yhigh, layers, OL)
      call apply_boundary_conditions(v_new, hfac_s, wetmask, &
                              xlow, xhigh, ylow, yhigh, layers, OL)

      ! Do the isopycnal layer physics
      if (.not. red_grav) then
        call barotropic_correction(h_new, u_new, v_new, eta, etanew, depth, a, &
            dx, dy, wetmask, hfac_w, hfac_s, dt, &
            maxits, eps, rjac, freesurf_fac, thickness_error, &
            debug_level, g_vec, nx, ny, layers, OL, &
            xlow, xhigh, ylow, yhigh, n, &
            MPI_COMM_WORLD, myid, num_procs, ilower, iupper, &
            hypre_grid, hypre_A, ierr)

      end if

      ! Apply the boundary conditions
      call apply_boundary_conditions(u_new, hfac_w, wetmask, &
                              xlow, xhigh, ylow, yhigh, layers, OL)
      call apply_boundary_conditions(v_new, hfac_s, wetmask, &
                              xlow, xhigh, ylow, yhigh, layers, OL)

      ! update global and tile halos with new values
      call update_halos(h_new, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)
      call update_halos(u_new, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)
      call update_halos(v_new, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)
      ! eta_new halo is updated in barotropic_correction

      ! Accumulate average fields
      if (avwrite .ne. 0) then
        hav = hav + h_new
        uav = uav + u_new
        vav = vav + v_new
        if (.not. red_grav) then
          etaav = eta + etanew
        end if
      end if

      ! Shuffle arrays: old -> very old,  present -> old, new -> present
      ! Height and velocity fields
      h = h_new
      u = u_new
      v = v_new
      if (.not. red_grav) then
        eta = etanew
      end if

      call maybe_dump_output(h, hav, u, uav, v, vav, eta, etaav, &
          dudt, dvdt, dhdt, AB_order, &
          wind_x, wind_y, nx, ny, layers, ilower, iupper, &
          xlow, xhigh, ylow, yhigh, OL, num_procs, myid, &
          n, nwrite, avwrite, checkpointwrite, diagwrite, &
          red_grav, dump_wind, debug_level)


      cur_time = time()
      if (myid .eq. 0) then
        if (cur_time - last_report_time > 3) then
          ! Three seconds passed since last report
          last_report_time = cur_time
          print "(A, I0, A, I0, A)", "Completed time step ", &
              n, " at ", cur_time - start_time, " seconds."
        end if
      end if

    end do

    cur_time = time()
    if (myid .eq. 0) then
      print "(A, I0, A, I0, A)", "Run finished at time step ", &
          n, ", in ", cur_time - start_time, " seconds."
    end if

    ! save checkpoint at end of every simulation
    call write_checkpoint_output(h, nx, ny, layers, ilower, iupper, &
                      xlow, xhigh, ylow, yhigh, OL, 1, &
                      n, 'checkpoints/h.', num_procs, myid)
    call write_checkpoint_output(u, nx, ny, layers, ilower, iupper, &
                      xlow, xhigh, ylow, yhigh, OL, 1, &
                      n, 'checkpoints/u.', num_procs, myid)
    call write_checkpoint_output(v, nx, ny, layers, ilower, iupper, &
                      xlow, xhigh, ylow, yhigh, OL, 1, &
                      n, 'checkpoints/v.', num_procs, myid)

    call write_checkpoint_output(dhdt, nx, ny, layers, ilower, iupper, &
                      xlow, xhigh, ylow, yhigh, OL, AB_order, &
                      n, 'checkpoints/dhdt.', num_procs, myid)
    call write_checkpoint_output(dudt, nx, ny, layers, ilower, iupper, &
                      xlow, xhigh, ylow, yhigh, OL, AB_order, &
                      n, 'checkpoints/dudt.', num_procs, myid)
    call write_checkpoint_output(dvdt, nx, ny, layers, ilower, iupper, &
                      xlow, xhigh, ylow, yhigh, OL, AB_order, &
                      n, 'checkpoints/dvdt.', num_procs, myid)

    if (.not. red_grav) then
      call write_checkpoint_output(eta, nx, ny, 1, ilower, iupper, &
                      xlow, xhigh, ylow, yhigh, OL, 1, &
                      n, 'checkpoints/eta.', num_procs, myid)
    end if


    return
  end subroutine model_run

end module model_main
