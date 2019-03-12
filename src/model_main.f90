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
      hfacW, hfacE, hfacS, hfacN, &
      fu, fv, &
      dt, au, ar, botDrag, kh, kv, slip, hmin, niter0, nTimeSteps, &
      dumpFreq, avFreq, checkpointFreq, diagFreq, &
      maxits, eps, freesurfFac, thickness_error, &
      debug_level, g_vec, rho0, &
      base_wind_x, base_wind_y, wind_mag_time_series, wind_depth, &
      spongeHTimeScale, spongeUTimeScale, spongeVTimeScale, &
      spongeH, spongeU, spongeV, &
      nx, ny, layers, OL, xlow, xhigh, ylow, yhigh, &
      RedGrav, hAdvecScheme, TS_algorithm, AB_order, &
      DumpWind, RelativeWind, Cd, start_time, &
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
    double precision, intent(in) :: hfacW(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfacE(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfacS(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfacN(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    ! Coriolis parameter at u and v grid-points respectively
    double precision, intent(in) :: fu(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: fv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    ! Numerics
    double precision, intent(in) :: dt, au, ar, botDrag
    double precision, intent(in) :: kh(layers), kv
    double precision, intent(in) :: slip, hmin
    integer,          intent(in) :: niter0, nTimeSteps
    double precision, intent(in) :: dumpFreq, avFreq, checkpointFreq, diagFreq
    integer,          intent(in) :: maxits
    double precision, intent(in) :: eps, freesurfFac, thickness_error
    integer,          intent(in) :: debug_level
    ! Physics
    double precision, intent(in) :: g_vec(layers)
    double precision, intent(in) :: rho0
    ! Wind
    double precision, intent(in) :: base_wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: base_wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: wind_mag_time_series(nTimeSteps)
    double precision, intent(in) :: wind_depth
    ! Sponge regions
    double precision, intent(in) :: spongeHTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeUTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeVTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeH(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeU(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeV(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    ! Resolution
    integer,          intent(in) :: nx, ny, layers, OL
    integer,          intent(in) :: xlow, xhigh, ylow, yhigh
    ! Reduced gravity vs n-layer physics
    logical,          intent(in) :: RedGrav
    integer,          intent(in) :: hAdvecScheme
    integer,          intent(in) :: TS_algorithm
    integer,          intent(in) :: AB_order
    ! Whether to write computed wind in the output
    logical,          intent(in) :: DumpWind
    logical,          intent(in) :: RelativeWind
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

    ! Time
    integer*8 :: last_report_time, cur_time


    if (myid .eq. 0) then
      if (RedGrav) then
        print "(A, I0, A, I0, A, I0, A, I0, A)", &
            "Running a reduced-gravity configuration of size ", &
            nx, "x", ny, "x", layers, " by ", nTimeSteps, " time steps."
      else
        print "(A, I0, A, I0, A, I0, A, I0, A)", &
            "Running an n-layer configuration of size ", &
            nx, "x", ny, "x", layers, " by ", nTimeSteps, " time steps."
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

    nwrite = int(dumpFreq/dt)
    avwrite = int(avFreq/dt)
    checkpointwrite = int(checkpointFreq/dt)
    diagwrite = int(diagFreq/dt)

    ! Pi, the constant
    pi = 3.1415926535897932384

    ! Initialize wind fields
    wind_x = base_wind_x*wind_mag_time_series(1)
    wind_y = base_wind_y*wind_mag_time_series(1)

    if (myid .eq. 0) then
      ! Initialise the diagnostic files
      call create_diag_file(layers, 'output/diagnostic.h.csv', 'h', niter0)
      call create_diag_file(layers, 'output/diagnostic.u.csv', 'u', niter0)
      call create_diag_file(layers, 'output/diagnostic.v.csv', 'v', niter0)
      if (.not. RedGrav) then
        call create_diag_file(1, 'output/diagnostic.eta.csv', 'eta', niter0)
      end if
    end if

    ! Initialise the average fields
    if (avwrite .ne. 0) then
      hav = 0.0
      uav = 0.0
      vav = 0.0
      if (.not. RedGrav) then
        etaav = 0.0
      end if
    end if

    ! initialise etanew
    etanew = 0d0

    call apply_boundary_conditions(u, hfacW, wetmask, &
                          xlow, xhigh, ylow, yhigh, layers, OL)
    call apply_boundary_conditions(v, hfacS, wetmask, &
                          xlow, xhigh, ylow, yhigh, layers, OL)


    if (.not. RedGrav) then
      ! Initialise arrays for pressure solver
      ! a = derivatives of the depth field
        call calc_A_matrix(a, depth, g_vec(1), dx, dy, OL, &
            xlow, xhigh, ylow, yhigh, freesurfFac, dt, &
            hfacW, hfacE, hfacS, hfacN)

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
          freesurfFac, thickness_error, &
          xlow, xhigh, ylow, yhigh, layers, OL)
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  Initialisation of the model STARTS HERE                            !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (niter0 .eq. 0) then

      n = 0
      call initialise_tendencies(dhdt, dudt, dvdt, h, u, v, depth, &
          dx, dy, dt, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, AB_order, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, debug_level)

    else if (niter0 .ne. 0) then
      n = niter0

      call load_checkpoint_files(dhdt, dudt, dvdt, h, u, v, eta, &
            RedGrav, niter0, nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
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

    do n = niter0+1, niter0+nTimeSteps

      wind_x = base_wind_x*wind_mag_time_series(n-niter0)
      wind_y = base_wind_y*wind_mag_time_series(n-niter0)

      call timestep(h_new, u_new, v_new, dhdt, dudt, dvdt, &
          h, u, v, depth, &
          dx, dy, dt, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, TS_algorithm, AB_order, &
          g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Apply the boundary conditions
      call apply_boundary_conditions(u_new, hfacW, wetmask, &
                              xlow, xhigh, ylow, yhigh, layers, OL)
      call apply_boundary_conditions(v_new, hfacS, wetmask, &
                              xlow, xhigh, ylow, yhigh, layers, OL)

      ! Do the isopycnal layer physics
      if (.not. RedGrav) then
        call barotropic_correction(h_new, u_new, v_new, eta, etanew, depth, a, &
            dx, dy, wetmask, hfacW, hfacS, dt, &
            maxits, eps, rjac, freesurfFac, thickness_error, &
            debug_level, g_vec, nx, ny, layers, OL, &
            xlow, xhigh, ylow, yhigh, n, &
            MPI_COMM_WORLD, myid, num_procs, ilower, iupper, &
            hypre_grid, hypre_A, ierr)

      end if

      ! Apply the boundary conditions
      call apply_boundary_conditions(u_new, hfacW, wetmask, &
                              xlow, xhigh, ylow, yhigh, layers, OL)
      call apply_boundary_conditions(v_new, hfacS, wetmask, &
                              xlow, xhigh, ylow, yhigh, layers, OL)

      ! update global and tile halos with new values
      call update_halos(h_new, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)
      call update_halos(u_new, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)
      call update_halos(v_new, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)
      if (.not. RedGrav) then
        call update_halos(etanew, nx, ny, 1, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)
      end if   

      ! Accumulate average fields
      if (avwrite .ne. 0) then
        hav = hav + h_new
        uav = uav + u_new
        vav = vav + v_new
        if (.not. RedGrav) then
          etaav = eta + etanew
        end if
      end if

      ! Shuffle arrays: old -> very old,  present -> old, new -> present
      ! Height and velocity fields
      h = h_new
      u = u_new
      v = v_new
      if (.not. RedGrav) then
        eta = etanew
      end if

      call maybe_dump_output(h, hav, u, uav, v, vav, eta, etaav, &
          dudt, dvdt, dhdt, AB_order, &
          wind_x, wind_y, nx, ny, layers, ilower, iupper, &
          xlow, xhigh, ylow, yhigh, OL, num_procs, myid, &
          n, nwrite, avwrite, checkpointwrite, diagwrite, &
          RedGrav, DumpWind, debug_level)


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

    if (.not. RedGrav) then
      call write_checkpoint_output(eta, nx, ny, 1, ilower, iupper, &
                      xlow, xhigh, ylow, yhigh, OL, 1, &
                      n, 'checkpoints/eta.', num_procs, myid)
    end if
    

    return
  end subroutine model_run

end module model_main
