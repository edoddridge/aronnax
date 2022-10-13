module time_stepping
  use state_deriv
  use boundaries
  use adams_bashforth
  use exchange
  implicit none

  contains


  ! ---------------------------------------------------------------------------
  !> Parse the ts_algorithm parameter to determine whether multiple tendency
  !! arrays need to be stored
  !! 1 is Forward Euler (part of both AB and RK families)
  !! 2-9 are AB family algorithms
  !! 12-19 are RK family algorithms (no need to store multiple tendencies, 
  !! hence AB_order = 1)


  subroutine set_AB_order(ts_algorithm, AB_order)
    implicit none

    integer,          intent(in)  :: ts_algorithm
    integer,          intent(out) :: AB_order

    if (ts_algorithm .eq. 1) then
      ! Forward Euler
      AB_order = 1
    else if (ts_algorithm .eq. 2) then
      ! Second-order AB
      AB_order = 2
    else if (ts_algorithm .eq. 3) then
      ! Third-order AB
      AB_order = 3
    else if (ts_algorithm .eq. 4) then
      ! Fourth-order AB
      AB_order = 4
    else if (ts_algorithm .eq. 5) then
      ! Fifth-order AB
      AB_order = 5
    else if (ts_algorithm .eq. 12) then
      ! Second-order RK
      AB_order = 1
    else if (ts_algorithm .eq. 13) then
      ! Third-order RK
      AB_order = 1
    else if (ts_algorithm .eq. 14) then
      ! Fourth-order RK
      AB_order = 1
    else
      ! ts_algorithm not set correctly
      call clean_stop(0, .FALSE.)
    end if

  end subroutine set_AB_order


  ! ---------------------------------------------------------------------------
  !> Initialise tendencies for linear multi-step timestepping algorithms

  subroutine initialise_tendencies(dhdt, dudt, dvdt, h, u, v, depth, &
          dx, dy, dt, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, AB_order, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
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
    double precision, intent(in) :: hfac_w(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfac_e(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfac_n(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfac_s(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: fu(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: fv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: au, ar, bot_drag
    double precision, intent(in) :: kh(layers), kv
    double precision, intent(in) :: hmin
    double precision, intent(in) :: slip
    logical,          intent(in) :: red_grav
    integer,          intent(in) :: h_advec_scheme
    integer,          intent(in) :: AB_order
    double precision, intent(in) :: g_vec(layers)
    double precision, intent(in) :: rho0
    double precision, intent(in) :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: wind_depth
    logical,          intent(in) :: relative_wind
    double precision, intent(in) :: Cd
    double precision, intent(in) :: sponge_h_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_u_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_v_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
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
          dx, dy, dt, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, 0, debug_level)

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
      ! Shuffle arrays: new -> present
      ! thickness and velocity fields
      h = h_new
      u = u_new
      v = v_new
    end do



  end subroutine initialise_tendencies


  ! ---------------------------------------------------------------------------
  !> Step the model forward one timestep. The timestepping scheme is set by the 
  !! ts_algorithm parameter

  subroutine timestep(h_new, u_new, v_new, dhdt, dudt, dvdt, &
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
    double precision, intent(in)    :: hfac_w(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: hfac_e(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: hfac_n(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: hfac_s(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: fu(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: fv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: au, ar, bot_drag
    double precision, intent(in)    :: kh(layers), kv
    double precision, intent(in)    :: hmin
    double precision, intent(in)    :: slip
    logical,          intent(in)    :: red_grav
    integer,          intent(in)    :: h_advec_scheme, ts_algorithm, AB_order
    double precision, intent(in)    :: g_vec(layers)
    double precision, intent(in)    :: rho0
    double precision, intent(in)    :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: wind_depth
    logical,          intent(in)    :: relative_wind
    double precision, intent(in)    :: Cd
    double precision, intent(in)    :: sponge_h_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: sponge_h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: sponge_u_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: sponge_u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: sponge_v_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: sponge_v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    integer,          intent(in)    :: nx, ny, layers
    integer,          intent(in)    :: ilower(0:num_procs-1,2)
    integer,          intent(in)    :: iupper(0:num_procs-1,2)
    integer,          intent(in)    :: xlow, xhigh, ylow, yhigh, OL
    integer,          intent(in)    :: num_procs, myid
    integer,          intent(in)    :: n
    integer,          intent(in)    :: debug_level



    if (ts_algorithm .eq. 1) then
      ! Forward Euler
      call state_derivative(dhdt(:,:,:,1), dudt(:,:,:,1), dvdt(:,:,:,1), h, u, v, depth, &
          dx, dy, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call forward_euler(h_new, dhdt, h, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call forward_euler(u_new, dudt, u, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call forward_euler(v_new, dvdt, v, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)


    else if (ts_algorithm .eq. 2) then
      ! Second-order AB
      call state_derivative(dhdt(:,:,:,1), dudt(:,:,:,1), dvdt(:,:,:,1), h, u, v, depth, &
          dx, dy, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call AB2(h_new, dhdt, h, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB2(u_new, dudt, u, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB2(v_new, dvdt, v, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)

    else if (ts_algorithm .eq. 3) then
      ! Third-order AB
      call state_derivative(dhdt(:,:,:,1), dudt(:,:,:,1), dvdt(:,:,:,1), h, u, v, depth, &
          dx, dy, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time with
      ! the Adams-Bashforth third order linear multistep method
      call AB3(h_new, dhdt, h, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB3(u_new, dudt, u, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB3(v_new, dvdt, v, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)

    else if (ts_algorithm .eq. 4) then
      ! Fourth-order AB

      call state_derivative(dhdt(:,:,:,1), dudt(:,:,:,1), dvdt(:,:,:,1), h, u, v, depth, &
          dx, dy, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call AB4(h_new, dhdt, h, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB4(u_new, dudt, u, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB4(v_new, dvdt, v, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)

    else if (ts_algorithm .eq. 5) then
      ! Fifth-order AB

      call state_derivative(dhdt(:,:,:,1), dudt(:,:,:,1), dvdt(:,:,:,1), h, u, v, depth, &
          dx, dy, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call AB5(h_new, dhdt, h, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB5(u_new, dudt, u, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)
      call AB5(v_new, dvdt, v, dt, xlow, xhigh, ylow, yhigh, layers, OL, AB_order)

    else if (ts_algorithm .eq. 12) then
      ! Second-order RK
      call RK2(h_new, u_new, v_new, dhdt, dudt, dvdt, h, u, v, depth, &
          dx, dy, dt, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

    else if (ts_algorithm .eq. 13) then
      ! Third-order RK
      call clean_stop(0, .FALSE.)

    else if (ts_algorithm .eq. 14) then
      ! Fourth-order RK
      call clean_stop(0, .FALSE.)

    else
      ! ts_algorithm not set correctly
      call clean_stop(0, .FALSE.)
    end if

  end subroutine timestep

  ! ---------------------------------------------------------------------------
  !> A second-order Runge-Kutta algorithm
  !! This saves the tendency arrays so that it can be used to initialise
  !! the linear multi-step methods (Adams-Bashforth algorithms)

  subroutine RK2(h_new, u_new, v_new, dhdt, dudt, dvdt, h, u, v, depth, &
          dx, dy, dt, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
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
    double precision, intent(in) :: hfac_w(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfac_e(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfac_n(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: hfac_s(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: fu(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: fv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: au, ar, bot_drag
    double precision, intent(in) :: kh(layers), kv
    double precision, intent(in) :: hmin
    double precision, intent(in) :: slip
    logical,          intent(in) :: red_grav
    integer,          intent(in) :: h_advec_scheme
    double precision, intent(in) :: g_vec(layers)
    double precision, intent(in) :: rho0
    double precision, intent(in) :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: wind_depth
    logical,          intent(in) :: relative_wind
    double precision, intent(in) :: Cd
    double precision, intent(in) :: sponge_h_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_u_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_v_time_scale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: sponge_v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
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
          dx, dy, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Calculate the values at half the time interval with Forward Euler
      hhalf = h+0.5d0*dt*dhdt
      uhalf = u+0.5d0*dt*dudt
      vhalf = v+0.5d0*dt*dvdt

      call state_derivative(dhdt, dudt, dvdt, &
          hhalf, uhalf, vhalf, depth, &
          dx, dy, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
          au, ar, bot_drag, kh, kv, hmin, slip, &
          red_grav, h_advec_scheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, relative_wind, Cd, &
          sponge_h_time_scale, sponge_h, &
          sponge_u_time_scale, sponge_u, &
          sponge_v_time_scale, sponge_v, &
          nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
          OL, num_procs, myid, n, debug_level)

      ! Calculate h, u, v with these tendencies
      h_new = h + dt*dhdt
      u_new = u + dt*dudt
      v_new = v + dt*dvdt

  end subroutine RK2

end module time_stepping
