module model_main
  use end_run
  use boundaries
  use barotropic_mode
  use state_deriv
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
      base_wind_x, base_wind_y, wind_mag_time_series, &
      spongeHTimeScale, spongeUTimeScale, spongeVTimeScale, &
      spongeH, spongeU, spongeV, &
      nx, ny, layers, RedGrav, hAdvecScheme, DumpWind, &
      RelativeWind, Cd, &
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
    ! Whether to write computed wind in the output
    logical,          intent(in) :: DumpWind
    logical,          intent(in) :: RelativeWind
    double precision,  intent(in) :: Cd

    double precision :: dhdt(0:nx+1, 0:ny+1, layers)
    double precision :: dhdtold(0:nx+1, 0:ny+1, layers)
    double precision :: dhdtveryold(0:nx+1, 0:ny+1, layers)
    double precision :: hnew(0:nx+1, 0:ny+1, layers)
    ! for initialisation
    double precision :: hhalf(0:nx+1, 0:ny+1, layers)
    ! for saving average fields
    double precision :: hav(0:nx+1, 0:ny+1, layers)

    double precision :: dudt(0:nx+1, 0:ny+1, layers)
    double precision :: dudtold(0:nx+1, 0:ny+1, layers)
    double precision :: dudtveryold(0:nx+1, 0:ny+1, layers)
    double precision :: unew(0:nx+1, 0:ny+1, layers)
    ! for initialisation
    double precision :: uhalf(0:nx+1, 0:ny+1, layers)
    ! for saving average fields
    double precision :: uav(0:nx+1, 0:ny+1, layers)

    double precision :: dvdt(0:nx+1, 0:ny+1, layers)
    double precision :: dvdtold(0:nx+1, 0:ny+1, layers)
    double precision :: dvdtveryold(0:nx+1, 0:ny+1, layers)
    double precision :: vnew(0:nx+1, 0:ny+1, layers)
    ! for initialisation
    double precision :: vhalf(0:nx+1, 0:ny+1, layers)
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

    ! dummy variable for loading checkpoints
    character(10)    :: num


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
      ! Do two initial time steps with Runge-Kutta second-order.
      ! These initialisation steps do NOT use or update the free surface.
      !
      ! ------------------------- negative 2 time step --------------------------
      ! Code to work out dhdtveryold, dudtveryold and dvdtveryold
      n = 0
      
      call state_derivative(dhdtveryold, dudtveryold, dvdtveryold, &
          h, u, v, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, n, debug_level)

      ! Calculate the values at half the time interval with Forward Euler
      hhalf = h+0.5d0*dt*dhdtveryold
      uhalf = u+0.5d0*dt*dudtveryold
      vhalf = v+0.5d0*dt*dvdtveryold

      call state_derivative(dhdtveryold, dudtveryold, dvdtveryold, &
          hhalf, uhalf, vhalf, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, n, debug_level)

      ! These are the values to be stored in the 'veryold' variables ready
      ! to start the proper model run.

      ! Calculate h, u, v with these tendencies
      h = h + dt*dhdtveryold
      u = u + dt*dudtveryold
      v = v + dt*dvdtveryold

      ! Apply the boundary conditions
      call apply_boundary_conditions(u, hfacW, wetmask, nx, ny, layers)
      call apply_boundary_conditions(v, hfacS, wetmask, nx, ny, layers)

      ! Wrap fields around for periodic simulations
      call wrap_fields_3D(u, nx, ny, layers)
      call wrap_fields_3D(v, nx, ny, layers)
      call wrap_fields_3D(h, nx, ny, layers)

      ! ------------------------- negative 1 time step --------------------------
      ! Code to work out dhdtold, dudtold and dvdtold

      call state_derivative(dhdtold, dudtold, dvdtold, &
          h, u, v, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, n, debug_level)

      ! Calculate the values at half the time interval with Forward Euler
      hhalf = h+0.5d0*dt*dhdtold
      uhalf = u+0.5d0*dt*dudtold
      vhalf = v+0.5d0*dt*dvdtold

      call state_derivative(dhdtold, dudtold, dvdtold, &
          hhalf, uhalf, vhalf, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, n, debug_level)

      ! These are the values to be stored in the 'old' variables ready to start
      ! the proper model run.

      ! Calculate h, u, v with these tendencies
      h = h + dt*dhdtold
      u = u + dt*dudtold
      v = v + dt*dvdtold

      ! Apply the boundary conditions
      call apply_boundary_conditions(u, hfacW, wetmask, nx, ny, layers)
      call apply_boundary_conditions(v, hfacS, wetmask, nx, ny, layers)

      ! Wrap fields around for periodic simulations
      call wrap_fields_3D(u, nx, ny, layers)
      call wrap_fields_3D(v, nx, ny, layers)
      call wrap_fields_3D(h, nx, ny, layers)

    else if (niter0 .ne. 0) then
      n = niter0

      ! load in the state and derivative arrays
      write(num, '(i10.10)') niter0

      open(unit=10, form='unformatted', file='checkpoints/h.'//num)
      read(10) h
      close(10)
      open(unit=10, form='unformatted', file='checkpoints/u.'//num)
      read(10) u
      close(10)
      open(unit=10, form='unformatted', file='checkpoints/v.'//num)
      read(10) v
      close(10)

      open(unit=10, form='unformatted', file='checkpoints/dhdt.'//num)
      read(10) dhdt
      close(10)
      open(unit=10, form='unformatted', file='checkpoints/dudt.'//num)
      read(10) dudt
      close(10)
      open(unit=10, form='unformatted', file='checkpoints/dvdt.'//num)
      read(10) dvdt
      close(10)

      open(unit=10, form='unformatted', file='checkpoints/dhdtold.'//num)
      read(10) dhdtold
      close(10)
      open(unit=10, form='unformatted', file='checkpoints/dudtold.'//num)
      read(10) dudtold
      close(10)
      open(unit=10, form='unformatted', file='checkpoints/dvdtold.'//num)
      read(10) dvdtold
      close(10)

      open(unit=10, form='unformatted', file='checkpoints/dhdtveryold.'//num)
      read(10) dhdtveryold
      close(10)
      open(unit=10, form='unformatted', file='checkpoints/dudtveryold.'//num)
      read(10) dudtveryold
      close(10)
      open(unit=10, form='unformatted', file='checkpoints/dvdtveryold.'//num)
      read(10) dvdtveryold
      close(10)

      if (.not. RedGrav) then
        open(unit=10, form='unformatted', file='checkpoints/eta.'//num)
        read(10) eta
        close(10)
      end if

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

      call state_derivative(dhdt, dudt, dvdt, h, u, v, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, n, debug_level)


      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time with
      ! the Adams-Bashforth third order linear multistep method

      unew = u + dt*(23d0*dudt - 16d0*dudtold + 5d0*dudtveryold)/12d0
      vnew = v + dt*(23d0*dvdt - 16d0*dvdtold + 5d0*dvdtveryold)/12d0
      hnew = h + dt*(23d0*dhdt - 16d0*dhdtold + 5d0*dhdtveryold)/12d0

      ! Apply the boundary conditions
      call apply_boundary_conditions(unew, hfacW, wetmask, nx, ny, layers)
      call apply_boundary_conditions(vnew, hfacS, wetmask, nx, ny, layers)

      ! Do the isopycnal layer physics
      if (.not. RedGrav) then
        call barotropic_correction(hnew, unew, vnew, eta, etanew, depth, a, &
            dx, dy, wetmask, hfacW, hfacS, dt, &
            maxits, eps, rjac, freesurfFac, thickness_error, &
            debug_level, g_vec, nx, ny, layers, n, &
            MPI_COMM_WORLD, myid, num_procs, ilower, iupper, &
            hypre_grid, hypre_A, ierr)

      end if


      ! Stop layers from getting too thin
      call enforce_minimum_layer_thickness(hnew, hmin, nx, ny, layers, n)

      ! Wrap fields around for periodic simulations
      call wrap_fields_3D(unew, nx, ny, layers)
      call wrap_fields_3D(vnew, nx, ny, layers)
      call wrap_fields_3D(hnew, nx, ny, layers)
      if (.not. RedGrav) then
        call wrap_fields_2D(etanew, nx, ny)
      end if    

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

      call maybe_dump_output(h, hav, u, uav, v, vav, eta, etaav, &
          dudt, dvdt, dhdt, &
          dudtold, dvdtold, dhdtold, &
          dudtveryold, dvdtveryold, dhdtveryold, &
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
        dudt, dvdt, dhdt, &
        dudtold, dvdtold, dhdtold, &
        dudtveryold, dvdtveryold, dhdtveryold, &
        wind_x, wind_y, nx, ny, layers, &
        n, n, n, n-1, n, &
        RedGrav, DumpWind, 0)

    return
  end subroutine model_run

end module model_main