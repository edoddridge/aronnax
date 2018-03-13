module time_stepping
  use state_deriv
  use boundaries
  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> A second-order Runge-Kutta algorithm
  !! This saves the tendency arrays so that it can be used to initialise
  !! the linear multi-step methods (Adams-Bashforth algorithms)

  subroutine RK2(h_new, u_new, v_new, dhdt, dudt, dvdt, h, u, v, depth, &
          dx, dy, dt, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, n, debug_level)


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
    double precision, intent(in) :: slip
    logical,          intent(in) :: RedGrav
    integer,          intent(in) :: hAdvecScheme
    double precision, intent(in) :: g_vec(layers)
    double precision, intent(in) :: rho0
    double precision, intent(in) :: wind_x(0:nx+1, 0:ny+1)
    double precision, intent(in) :: wind_y(0:nx+1, 0:ny+1)
    logical,          intent(in) :: RelativeWind
    double precision, intent(in) :: Cd
    double precision, intent(in) :: spongeHTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeH(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeUTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeU(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeVTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: spongeV(0:nx+1, 0:ny+1, layers)
    integer, intent(in) :: nx, ny, layers
    integer, intent(in) :: n
    integer, intent(in) :: debug_level


    double precision :: hhalf(0:nx+1, 0:ny+1, layers)
    double precision :: uhalf(0:nx+1, 0:ny+1, layers)
    double precision :: vhalf(0:nx+1, 0:ny+1, layers)


      call state_derivative(dhdt, dudt, dvdt, &
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
      hhalf = h+0.5d0*dt*dhdt
      uhalf = u+0.5d0*dt*dudt
      vhalf = v+0.5d0*dt*dvdt

      call state_derivative(dhdt, dudt, dvdt, &
          hhalf, uhalf, vhalf, depth, &
          dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
          au, ar, botDrag, kh, kv, slip, &
          RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
          RelativeWind, Cd, &
          spongeHTimeScale, spongeH, &
          spongeUTimeScale, spongeU, &
          spongeVTimeScale, spongeV, &
          nx, ny, layers, n, debug_level)

      ! Calculate h, u, v with these tendencies
      h_new = h + dt*dhdt
      u_new = u + dt*dudt
      v_new = v + dt*dvdt

      ! Apply the boundary conditions
      call apply_boundary_conditions(u_new, hfacW, wetmask, nx, ny, layers)
      call apply_boundary_conditions(v_new, hfacS, wetmask, nx, ny, layers)

      ! Wrap fields around for periodic simulations
      call wrap_fields_3D(u_new, nx, ny, layers)
      call wrap_fields_3D(v_new, nx, ny, layers)
      call wrap_fields_3D(h_new, nx, ny, layers)

  end subroutine RK2

end module time_stepping