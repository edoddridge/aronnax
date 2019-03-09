module time_stepping
  use state_deriv
  use boundaries
  use adams_bashforth
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
          nx, ny, layers, OL, debug_level)
    implicit none

    double precision, intent(out) :: dhdt(0:nx+1, 0:ny+1, layers, AB_order)
    double precision, intent(out) :: dudt(0:nx+1, 0:ny+1, layers, AB_order)
    double precision, intent(out) :: dvdt(0:nx+1, 0:ny+1, layers, AB_order)
    double precision, intent(inout) :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(inout) :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(inout) :: v(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: depth(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: dx, dy, dt
    double precision, intent(in) :: wetmask(0:nx+1, 0:ny+1)
    double precision, intent(in) :: hfacW(0:nx+1, 0:ny+1)
    double precision, intent(in) :: hfacE(0:nx+1, 0:ny+1)
    double precision, intent(in) :: hfacN(0:nx+1, 0:ny+1)
    double precision, intent(in) :: hfacS(0:nx+1, 0:ny+1)
    double precision, intent(in) :: fu(0:nx+1, 0:ny+1)
    double precision, intent(in) :: fv(0:nx+1, 0:ny+1)
    double precision, intent(in) :: au, ar, botDrag
    double precision, intent(in) :: kh(layers), kv
    double precision, intent(in) :: hmin
    double precision, intent(in) :: slip
    logical,          intent(in) :: RedGrav
    integer,          intent(in) :: hAdvecScheme
    integer,          intent(in) :: AB_order
    double precision, intent(in) :: g_vec(layers)
    double precision, intent(in) :: rho0
    double precision, intent(in) :: wind_x(0:nx+1, 0:ny+1)
    double precision, intent(in) :: wind_y(0:nx+1, 0:ny+1)
    double precision, intent(in) :: wind_depth
    logical,          intent(in) :: RelativeWind
    double precision, intent(in) :: Cd
    double precision, intent(in) :: spongeHTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeH(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeUTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeU(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeVTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeV(0:nx+1, 0:ny+1, layers)
    integer, intent(in) :: nx, ny, layers, OL
    integer, intent(in) :: debug_level



    double precision :: h_new(0:nx+1, 0:ny+1, layers)
    double precision :: u_new(0:nx+1, 0:ny+1, layers)
    double precision :: v_new(0:nx+1, 0:ny+1, layers)
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
          nx, ny, layers, OL, 0, debug_level)

            ! Apply the boundary conditions
      call apply_boundary_conditions(u_new, hfacW, wetmask, nx, ny, layers, OL)
      call apply_boundary_conditions(v_new, hfacS, wetmask, nx, ny, layers, OL)

      ! Wrap fields around for periodic simulations
      call wrap_fields_3D(u_new, nx, ny, layers, OL)
      call wrap_fields_3D(v_new, nx, ny, layers, OL)
      call wrap_fields_3D(h_new, nx, ny, layers, OL)

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
      nx, ny, layers, OL, n, debug_level)
    implicit none


    double precision, intent(out)   :: h_new(0:nx+1, 0:ny+1, layers)
    double precision, intent(out)   :: u_new(0:nx+1, 0:ny+1, layers)
    double precision, intent(out)   :: v_new(0:nx+1, 0:ny+1, layers)
    double precision, intent(inout) :: dhdt(0:nx+1, 0:ny+1, layers, AB_order)
    double precision, intent(inout) :: dudt(0:nx+1, 0:ny+1, layers, AB_order)
    double precision, intent(inout) :: dvdt(0:nx+1, 0:ny+1, layers, AB_order)
    double precision, intent(in)    :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: v(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: depth(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: dx, dy, dt
    double precision, intent(in)    :: wetmask(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: hfacW(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: hfacE(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: hfacN(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: hfacS(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: fu(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: fv(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: au, ar, botDrag
    double precision, intent(in)    :: kh(layers), kv
    double precision, intent(in)    :: hmin
    double precision, intent(in)    :: slip
    logical,          intent(in)    :: RedGrav
    integer,          intent(in)    :: hAdvecScheme, TS_algorithm, AB_order
    double precision, intent(in)    :: g_vec(layers)
    double precision, intent(in)    :: rho0
    double precision, intent(in)    :: wind_x(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: wind_y(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: wind_depth
    logical,          intent(in)    :: RelativeWind
    double precision, intent(in)    :: Cd
    double precision, intent(in)    :: spongeHTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: spongeH(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: spongeUTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: spongeU(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: spongeVTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: spongeV(0:nx+1, 0:ny+1, layers)
    integer,          intent(in)    :: nx, ny, layers, OL
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
          nx, ny, layers, OL, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call ForwardEuler(h_new, dhdt, h, dt, nx, ny, layers, OL, AB_order)
      call ForwardEuler(u_new, dudt, u, dt, nx, ny, layers, OL, AB_order)
      call ForwardEuler(v_new, dvdt, v, dt, nx, ny, layers, OL, AB_order)


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
          nx, ny, layers, OL, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call AB2(h_new, dhdt, h, dt, nx, ny, layers, OL, AB_order)
      call AB2(u_new, dudt, u, dt, nx, ny, layers, OL, AB_order)
      call AB2(v_new, dvdt, v, dt, nx, ny, layers, OL, AB_order)

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
          nx, ny, layers, OL, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time with
      ! the Adams-Bashforth third order linear multistep method
      call AB3(h_new, dhdt, h, dt, nx, ny, layers, OL, AB_order)
      call AB3(u_new, dudt, u, dt, nx, ny, layers, OL, AB_order)
      call AB3(v_new, dvdt, v, dt, nx, ny, layers, OL, AB_order)

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
          nx, ny, layers, OL, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call AB4(h_new, dhdt, h, dt, nx, ny, layers, OL, AB_order)
      call AB4(u_new, dudt, u, dt, nx, ny, layers, OL, AB_order)
      call AB4(v_new, dvdt, v, dt, nx, ny, layers, OL, AB_order)

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
          nx, ny, layers, OL, n, debug_level)

      ! Use dh/dt, du/dt and dv/dt to step h, u and v forward in time
      call AB5(h_new, dhdt, h, dt, nx, ny, layers, OL, AB_order)
      call AB5(u_new, dudt, u, dt, nx, ny, layers, OL, AB_order)
      call AB5(v_new, dvdt, v, dt, nx, ny, layers, OL, AB_order)

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
          nx, ny, layers, OL, n, debug_level)

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
          nx, ny, layers, OL, n, debug_level)
    implicit none

    double precision, intent(out) :: h_new(0:nx+1, 0:ny+1, layers)
    double precision, intent(out) :: u_new(0:nx+1, 0:ny+1, layers)
    double precision, intent(out) :: v_new(0:nx+1, 0:ny+1, layers)
    double precision, intent(out) :: dhdt(0:nx+1, 0:ny+1, layers)
    double precision, intent(out) :: dudt(0:nx+1, 0:ny+1, layers)
    double precision, intent(out) :: dvdt(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: v(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: depth(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: dx, dy, dt
    double precision, intent(in) :: wetmask(0:nx+1, 0:ny+1)
    double precision, intent(in) :: hfacW(0:nx+1, 0:ny+1)
    double precision, intent(in) :: hfacE(0:nx+1, 0:ny+1)
    double precision, intent(in) :: hfacN(0:nx+1, 0:ny+1)
    double precision, intent(in) :: hfacS(0:nx+1, 0:ny+1)
    double precision, intent(in) :: fu(0:nx+1, 0:ny+1)
    double precision, intent(in) :: fv(0:nx+1, 0:ny+1)
    double precision, intent(in) :: au, ar, botDrag
    double precision, intent(in) :: kh(layers), kv
    double precision, intent(in) :: hmin
    double precision, intent(in) :: slip
    logical,          intent(in) :: RedGrav
    integer,          intent(in) :: hAdvecScheme
    double precision, intent(in) :: g_vec(layers)
    double precision, intent(in) :: rho0
    double precision, intent(in) :: wind_x(0:nx+1, 0:ny+1)
    double precision, intent(in) :: wind_y(0:nx+1, 0:ny+1)
    double precision, intent(in) :: wind_depth
    logical,          intent(in) :: RelativeWind
    double precision, intent(in) :: Cd
    double precision, intent(in) :: spongeHTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeH(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeUTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeU(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeVTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeV(0:nx+1, 0:ny+1, layers)
    integer, intent(in) :: nx, ny, layers, OL
    integer, intent(in) :: n
    integer, intent(in) :: debug_level


    double precision :: hhalf(0:nx+1, 0:ny+1, layers)
    double precision :: uhalf(0:nx+1, 0:ny+1, layers)
    double precision :: vhalf(0:nx+1, 0:ny+1, layers)


      call state_derivative(dhdt, dudt, dvdt, &
          h, u, v, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, hmin, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          wind_depth, RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, OL, n, debug_level)

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
          nx, ny, layers, OL, n, debug_level)

      ! Calculate h, u, v with these tendencies
      h_new = h + dt*dhdt
      u_new = u + dt*dudt
      v_new = v + dt*dvdt

  end subroutine RK2

end module time_stepping
