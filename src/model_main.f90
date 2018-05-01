module model_main
  use end_run
  use boundaries
  use barotropic_mode
  use state_deriv
  use time_stepping
  use enforce_thickness
  implicit none

  contains

  ! ------------------------------ Primary routine ----------------------------
  !> Run the model

  subroutine model_run(h, u, v, eta, depth, dx, dy, wetmask, fu, fv, &
      dt, au, ar, botDrag, kh, kv, slip, hmin, niter0, nTimeSteps, &
      dumpFreq, avFreq, checkpointFreq, diagFreq, &
      maxits, eps, freesurfFac, thickness_error, &
      debug_level, g_vec, rho0, &
      base_wind_x, base_wind_y, wind_mag_time_series, wind_depth, &
      spongeHTimeScale, spongeUTimeScale, spongeVTimeScale, &
      spongeH, spongeU, spongeV, &
      nx, ny, layers, RedGrav, hAdvecScheme, TS_algorithm, AB_order, &
      DumpWind, RelativeWind, Cd, &
      MPI_COMM_WORLD, myid, num_procs, ilower, iupper, &
      hypre_grid)
    implicit none

    ! Layer thickness (h)
    double precision, intent(inout) :: h(0:nx+1, 0:ny+1, layers)
    ! Velocity component (u)
    double precision, intent(inout) :: u(0:nx+1, 0:ny+1, layers)
    ! Velocity component (v)
    double precision, intent(inout) :: v(0:nx+1, 0:ny+1, layers)
    ! Free surface (eta)
    double precision, intent(inout) :: eta(0:nx+1, 0:ny+1)
    ! Bathymetry
    double precision, intent(in) :: depth(0:nx+1, 0:ny+1)
    ! Grid
    double precision, intent(in) :: dx, dy
    double precision, intent(in) :: wetmask(0:nx+1, 0:ny+1)
    ! Coriolis parameter at u and v grid-points respectively
    double precision, intent(in) :: fu(0:nx+1, 0:ny+1)
    double precision, intent(in) :: fv(0:nx+1, 0:ny+1)
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
    double precision, intent(in) :: base_wind_x(0:nx+1, 0:ny+1)
    double precision, intent(in) :: base_wind_y(0:nx+1, 0:ny+1)
    double precision, intent(in) :: wind_mag_time_series(nTimeSteps)
    double precision, intent(in) :: wind_depth
    ! Sponge regions
    double precision, intent(in) :: spongeHTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeUTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeVTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeH(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeU(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeV(0:nx+1, 0:ny+1, layers)
    ! Resolution
    integer,          intent(in) :: nx, ny, layers
    ! Reduced gravity vs n-layer physics
    logical,          intent(in) :: RedGrav
    integer,          intent(in) :: hAdvecScheme
    integer,          intent(in) :: TS_algorithm
    integer,          intent(in) :: AB_order
    ! Whether to write computed wind in the output
    logical,          intent(in) :: DumpWind
    logical,          intent(in) :: RelativeWind
    double precision,  intent(in) :: Cd

    double precision :: dhdt(0:nx+1, 0:ny+1, layers, AB_order)
    double precision :: h_new(0:nx+1, 0:ny+1, layers)
    ! for saving average fields
    double precision :: hav(0:nx+1, 0:ny+1, layers)

    double precision :: dudt(0:nx+1, 0:ny+1, layers, AB_order)
    double precision :: u_new(0:nx+1, 0:ny+1, layers)
    ! for saving average fields
    double precision :: uav(0:nx+1, 0:ny+1, layers)

    double precision :: dvdt(0:nx+1, 0:ny+1, layers, AB_order)
    double precision :: v_new(0:nx+1, 0:ny+1, layers)
    ! for saving average fields
    double precision :: vav(0:nx+1, 0:ny+1, layers)

    double precision :: etanew(0:nx+1, 0:ny+1)
    ! for saving average fields
    double precision :: etaav(0:nx+1, 0:ny+1)

    ! Pressure solver variables
    double precision :: a(5, nx, ny)

    ! Geometry
    double precision :: hfacW(0:nx+1, 0:ny+1)
    double precision :: hfacE(0:nx+1, 0:ny+1)
    double precision :: hfacN(0:nx+1, 0:ny+1)
    double precision :: hfacS(0:nx+1, 0:ny+1)

    ! Numerics
    double precision :: pi
    double precision :: rjac

    ! dumping output
    integer :: nwrite, avwrite, checkpointwrite, diagwrite

    ! External solver variables
    integer          :: offsets(2,5)
    integer          :: i, j ! loop variables
    double precision :: values(nx * ny)
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
    double precision :: wind_x(0:nx+1, 0:ny+1)
    double precision :: wind_y(0:nx+1, 0:ny+1)

    ! Time
    integer*8 :: start_time, last_report_time, cur_time


    start_time = time()
    if (RedGrav) then
      print "(A, I0, A, I0, A, I0, A, I0, A)", &
          "Running a reduced-gravity configuration of size ", &
          nx, "x", ny, "x", layers, " by ", nTimeSteps, " time steps."
    else
      print "(A, I0, A, I0, A, I0, A, I0, A)", &
          "Running an n-layer configuration of size ", &
          nx, "x", ny, "x", layers, " by ", nTimeSteps, " time steps."
    end if

    if (myid .eq. 0) then
      ! Show the domain decomposition
      print "(A)", "Domain decomposition:"
      print "(A, I0)", 'ilower (x) = ', ilower(:,1)
      print "(A, I0)", 'ilower (y) = ', ilower(:,2)
      print "(A, I0)", 'iupper (x) = ', iupper(:,1)
      print "(A, I0)", 'iupper (y) = ', iupper(:,2)
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

    ! Initialise the diagnostic files
    call create_diag_file(layers, 'output/diagnostic.h.csv', 'h', niter0)
    call create_diag_file(layers, 'output/diagnostic.u.csv', 'u', niter0)
    call create_diag_file(layers, 'output/diagnostic.v.csv', 'v', niter0)
    if (.not. RedGrav) then
      call create_diag_file(1, 'output/diagnostic.eta.csv', 'eta', niter0)
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

    call calc_boundary_masks(wetmask, hfacW, hfacE, hfacS, hfacN, nx, ny)

    call apply_boundary_conditions(u, hfacW, wetmask, nx, ny, layers)
    call apply_boundary_conditions(v, hfacS, wetmask, nx, ny, layers)


    if (.not. RedGrav) then
      ! Initialise arrays for pressure solver
      ! a = derivatives of the depth field
        call calc_A_matrix(a, depth, g_vec(1), dx, dy, nx, ny, freesurfFac, dt, &
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
            a, nx, ny, ierr)
#endif

      ! Check that the supplied free surface anomaly and layer
      ! thicknesses are consistent with the supplied depth field.
      ! If they are not, then scale the layer thicknesses to make
      ! them consistent.
      call enforce_depth_thickness_consistency(h, eta, depth, &
          freesurfFac, thickness_error, nx, ny, layers)
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
          nx, ny, layers, debug_level)

    else if (niter0 .ne. 0) then
      n = niter0

      call load_checkpoint_files(dhdt, dudt, dvdt, h, u, v, eta, &
            RedGrav, niter0, nx, ny, layers, AB_order)

    end if

    ! Now the model is ready to start.
    ! - We have h, u, v at the zeroth time step, and the tendencies at
    !   two older time steps.
    ! - The model then solves for the tendencies at the current step
    !   before solving for the fields at the next time step.

    cur_time = time()
    if (cur_time - start_time .eq. 1) then
      print "(A)", "Initialized in 1 second."
    else
      print "(A, I0, A)", "Initialized in " , cur_time - start_time, " seconds."
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
          nx, ny, layers, n, debug_level)

      ! Apply the boundary conditions
      call apply_boundary_conditions(u_new, hfacW, wetmask, nx, ny, layers)
      call apply_boundary_conditions(v_new, hfacS, wetmask, nx, ny, layers)

      ! Do the isopycnal layer physics
      if (.not. RedGrav) then
        call barotropic_correction(h_new, u_new, v_new, eta, etanew, depth, a, &
            dx, dy, wetmask, hfacW, hfacS, dt, &
            maxits, eps, rjac, freesurfFac, thickness_error, &
            debug_level, g_vec, nx, ny, layers, n, &
            MPI_COMM_WORLD, myid, num_procs, ilower, iupper, &
            hypre_grid, hypre_A, ierr)

      end if


      ! Stop layers from getting too thin
      ! call enforce_minimum_layer_thickness(h_new, hmin, nx, ny, layers, n)

      ! Wrap fields around for periodic simulations
      call wrap_fields_3D(u_new, nx, ny, layers)
      call wrap_fields_3D(v_new, nx, ny, layers)
      call wrap_fields_3D(h_new, nx, ny, layers)
      if (.not. RedGrav) then
        call wrap_fields_2D(etanew, nx, ny)
      end if    

      ! Apply the boundary conditions
      call apply_boundary_conditions(u_new, hfacW, wetmask, nx, ny, layers)
      call apply_boundary_conditions(v_new, hfacS, wetmask, nx, ny, layers)

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
          wind_x, wind_y, nx, ny, layers, &
          n, nwrite, avwrite, checkpointwrite, diagwrite, &
          RedGrav, DumpWind, debug_level)


      cur_time = time()
      if (cur_time - last_report_time > 3) then
        ! Three seconds passed since last report
        last_report_time = cur_time
        print "(A, I0, A, I0, A)", "Completed time step ", &
            n, " at ", cur_time - start_time, " seconds."
      end if

    end do

    cur_time = time()
    print "(A, I0, A, I0, A)", "Run finished at time step ", &
        n, ", in ", cur_time - start_time, " seconds."

    ! save checkpoint at end of every simulation
    call maybe_dump_output(h, hav, u, uav, v, vav, eta, etaav, &
        dudt, dvdt, dhdt, AB_order, &
        wind_x, wind_y, nx, ny, layers, &
        n, n, n, n-1, n, &
        RedGrav, DumpWind, 0)

    return
  end subroutine model_run

end module model_main
