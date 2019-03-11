module time_stepping
  use state_deriv
  use boundaries
  use adams_bashforth
  use exchange
  implicit none

  contains

  
  ! ---------------------------------------------------------------------------
  !> Parse the TS_algorithm parameter to determine whether multiple tendency
  !! arrays need to be stored
  !! 1 is Forward Euler (part of both AB and RK families)
  !! 2-9 are AB family algorithms
  !! 12-19 are RK family algorithms (no need to store multiple tendencies, 
  !! hence AB_order = 1)


  subroutine set_AB_order(TS_algorithm, AB_order)
    implicit none

    integer,          intent(in)  :: TS_algorithm
    integer,          intent(out) :: AB_order

    if (TS_algorithm .eq. 1) then
      ! Forward Euler
      AB_order = 1
    else if (TS_algorithm .eq. 2) then
      ! Second-order AB
      AB_order = 2
    else if (TS_algorithm .eq. 3) then
      ! Third-order AB
      AB_order = 3
    else if (TS_algorithm .eq. 4) then
      ! Fourth-order AB
      AB_order = 4
    else if (TS_algorithm .eq. 5) then
      ! Fifth-order AB
      AB_order = 5
    else if (TS_algorithm .eq. 12) then
      ! Second-order RK
      AB_order = 1
    else if (TS_algorithm .eq. 13) then
      ! Third-order RK
      AB_order = 1
    else if (TS_algorithm .eq. 14) then
      ! Fourth-order RK
      AB_order = 1
    else
      ! TS_algorithm not set correctly
      call clean_stop(0, .FALSE.)
    end if

  end subroutine set_AB_order

  
  ! ---------------------------------------------------------------------------
  !> Initialise tendencies for linear multi-step timestepping algorithms

  subroutine initialise_tendencies(dhdt, dudt, dvdt, h, u, v, depth, &
          dx, dy, dt, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, AB_order, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, debug_level)
    implicit none


    double precision, intent(out) :: dhdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision, intent(out) :: dudt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision, intent(out) :: dvdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision, intent(inout) :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(inout) :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(inout) :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: depth(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: dx, dy, dt
    double precision, intent(in) :: wetmask(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfacW(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfacE(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfacN(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfacS(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: fu(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: fv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: au, ar, botDrag
    double precision, intent(in) :: kh(layers), kv
    double precision, intent(in) :: hmin
    double precision, intent(in) :: slip
    logical,          intent(in) :: RedGrav
    integer,          intent(in) :: hAdvecScheme
    integer,          intent(in) :: AB_order
    double precision, intent(in) :: g_vec(layers)
    double precision, intent(in) :: rho0
    double precision, intent(in) :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: wind_depth
    logical,          intent(in) :: RelativeWind
    double precision, intent(in) :: Cd
    double precision, intent(in) :: spongeHTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeH(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeUTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeU(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeVTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeV(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    integer,          intent(in) :: nx, ny, layers
    integer,          intent(in) :: ilower(0:num_procs-1,2)
    integer,          intent(in) :: iupper(0:num_procs-1,2)
    integer,          intent(in) :: xlow, xhigh, ylow, yhigh, OL
    integer,          intent(in) :: num_procs, myid
    integer,          intent(in) :: debug_level



    double precision :: h_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: u_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: v_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    integer          :: t

    ! Do some initial time steps with Runge-Kutta second-order.
    ! These initialisation steps do NOT use or update the free surface.

    do t = 0, AB_order-2

      call RK2(h_new, u_new, v_new, &
          dhdt(:,:,:,AB_order-t), dudt(:,:,:,AB_order-t), dvdt(:,:,:,AB_order-t), &
          h, u, v, depth, &
          dx, dy, dt, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, 0, debug_level)

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
      ! Shuffle arrays: new -> present
      ! thickness and velocity fields
      h = h_new
      u = u_new
      v = v_new
    end do



  end subroutine initialise_tendencies


  ! ---------------------------------------------------------------------------
  !> Step the model forward one timestep. The timestepping scheme is set by the 
  !! TS_algorithm parameter

  subroutine timestep(h_new, u_new, v_new, dhdt, dudt, dvdt, &
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
    implicit none

    double precision, intent(out)   :: h_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out)   :: u_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out)   :: v_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(inout) :: dhdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision, intent(inout) :: dudt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision, intent(inout) :: dvdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision, intent(in)    :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: depth(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: dx, dy, dt
    double precision, intent(in)    :: wetmask(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: hfacW(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: hfacE(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: hfacN(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: hfacS(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: fu(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: fv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: au, ar, botDrag
    double precision, intent(in)    :: kh(layers), kv
    double precision, intent(in)    :: hmin
    double precision, intent(in)    :: slip
    logical,          intent(in)    :: RedGrav
    integer,          intent(in)    :: hAdvecScheme, TS_algorithm, AB_order
    double precision, intent(in)    :: g_vec(layers)
    double precision, intent(in)    :: rho0
    double precision, intent(in)    :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: wind_depth
    logical,          intent(in)    :: RelativeWind
    double precision, intent(in)    :: Cd
    double precision, intent(in)    :: spongeHTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: spongeH(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: spongeUTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: spongeU(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: spongeVTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: spongeV(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    integer,          intent(in)    :: nx, ny, layers
    integer,          intent(in)    :: ilower(0:num_procs-1,2)
    integer,          intent(in)    :: iupper(0:num_procs-1,2)
    integer,          intent(in)    :: xlow, xhigh, ylow, yhigh, OL
    integer,          intent(in)    :: num_procs, myid
    integer,          intent(in)    :: n
    integer,          intent(in)    :: debug_level



    if (TS_algorithm .eq. 1) then
      ! Forward Euler
      call state_derivative(dhdt(:,:,:,1), dudt(:,:,:,1), dvdt(:,:,:,1), h, u, v, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call ForwardEuler(h_new, dhdt, h, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call ForwardEuler(u_new, dudt, u, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call ForwardEuler(v_new, dvdt, v, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)


    else if (TS_algorithm .eq. 2) then
      ! Second-order AB
      call state_derivative(dhdt(:,:,:,1), dudt(:,:,:,1), dvdt(:,:,:,1), h, u, v, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call AB2(h_new, dhdt, h, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB2(u_new, dudt, u, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB2(v_new, dvdt, v, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)

    else if (TS_algorithm .eq. 3) then
      ! Third-order AB
      call state_derivative(dhdt(:,:,:,1), dudt(:,:,:,1), dvdt(:,:,:,1), h, u, v, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time with
      ! the Adams-Bashforth third order linear multistep method
      call AB3(h_new, dhdt, h, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB3(u_new, dudt, u, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB3(v_new, dvdt, v, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)

    else if (TS_algorithm .eq. 4) then
      ! Fourth-order AB

      call state_derivative(dhdt(:,:,:,1), dudt(:,:,:,1), dvdt(:,:,:,1), h, u, v, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call AB4(h_new, dhdt, h, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB4(u_new, dudt, u, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB4(v_new, dvdt, v, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)

    else if (TS_algorithm .eq. 5) then
      ! Fifth-order AB

      call state_derivative(dhdt(:,:,:,1), dudt(:,:,:,1), dvdt(:,:,:,1), h, u, v, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call AB5(h_new, dhdt, h, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB5(u_new, dudt, u, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB5(v_new, dvdt, v, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)

    else if (TS_algorithm .eq. 12) then
      ! Second-order RK
      call RK2(h_new, u_new, v_new, dhdt, dudt, dvdt, h, u, v, depth, &
          dx, dy, dt, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

    else if (TS_algorithm .eq. 13) then
      ! Third-order RK
      call clean_stop(0, .FALSE.)

    else if (TS_algorithm .eq. 14) then
      ! Fourth-order RK
      call clean_stop(0, .FALSE.)

    else
      ! TS_algorithm not set correctly
      call clean_stop(0, .FALSE.)
    end if

  end subroutine timestep

  ! ---------------------------------------------------------------------------
  !> A second-order Runge-Kutta algorithm
  !! This saves the tendency arrays so that it can be used to initialise
  !! the linear multi-step methods (Adams-Bashforth algorithms)

  subroutine RK2(h_new, u_new, v_new, dhdt, dudt, dvdt, h, u, v, depth, &
          dx, dy, dt, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)
    implicit none

    double precision, intent(out) :: h_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out) :: u_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out) :: v_new(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out) :: dhdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out) :: dudt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out) :: dvdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: depth(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: dx, dy, dt
    double precision, intent(in) :: wetmask(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfacW(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfacE(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfacN(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfacS(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: fu(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: fv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: au, ar, botDrag
    double precision, intent(in) :: kh(layers), kv
    double precision, intent(in) :: hmin
    double precision, intent(in) :: slip
    logical,          intent(in) :: RedGrav
    integer,          intent(in) :: hAdvecScheme
    double precision, intent(in) :: g_vec(layers)
    double precision, intent(in) :: rho0
    double precision, intent(in) :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: wind_depth
    logical,          intent(in) :: RelativeWind
    double precision, intent(in) :: Cd
    double precision, intent(in) :: spongeHTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeH(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeUTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeU(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeVTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: spongeV(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    integer,          intent(in) :: nx, ny, layers
    integer,          intent(in) :: ilower(0:num_procs-1,2)
    integer,          intent(in) :: iupper(0:num_procs-1,2)
    integer,          intent(in) :: xlow, xhigh, ylow, yhigh, OL
    integer,          intent(in) :: num_procs, myid
    integer,          intent(in) :: n
    integer,          intent(in) :: debug_level


    double precision :: hhalf(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: uhalf(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: vhalf(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)


      call state_derivative(dhdt, dudt, dvdt, &
          h, u, v, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Calculate the values at half the time interval with Forward Euler
      hhalf = h+0.5d0*dt*dhdt
      uhalf = u+0.5d0*dt*dudt
      vhalf = v+0.5d0*dt*dvdt

      call state_derivative(dhdt, dudt, dvdt, &
          hhalf, uhalf, vhalf, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Calculate h, u, v with these tendencies
      h_new = h + dt*dhdt
      u_new = u + dt*dudt
      v_new = v + dt*dvdt

  end subroutine RK2

end module time_stepping
