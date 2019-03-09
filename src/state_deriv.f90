module state_deriv
  use thickness
  use io
  use bernoulli
  use vorticity
  use momentum

  implicit none

  contains

  ! ----------------------------- Auxiliary routines --------------------------
  !> Compute the forward state derivative

  subroutine state_derivative(dhdt, dudt, dvdt, h, u, v, depth, &
      dx, dy, wetmask, hfacW, hfacE, hfacN, hfacS, fu, fv, &
      au, ar, botDrag, kh, kv, hmin, slip, &
      RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
      wind_depth, RelativeWind, Cd, &
      spongeHTimeScale, spongeH, &
      spongeUTimeScale, spongeU, &
      spongeVTimeScale, spongeV, &
      nx, ny, layers, OL, n, debug_level)
    implicit none

    double precision, intent(out) :: dhdt(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(out) :: dudt(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(out) :: dvdt(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: h(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: u(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: v(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: depth(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: dx, dy
    double precision, intent(in) :: wetmask(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(in) :: hfacW(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(in) :: hfacE(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(in) :: hfacN(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(in) :: hfacS(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(in) :: fu(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(in) :: fv(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(in) :: au, ar, botDrag
    double precision, intent(in) :: kh(layers), kv
    double precision, intent(in) :: hmin
    double precision, intent(in) :: slip
    logical,          intent(in) :: RedGrav
    integer,          intent(in) :: hAdvecScheme
    double precision, intent(in) :: g_vec(layers)
    double precision, intent(in) :: rho0
    double precision, intent(in) :: wind_x(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(in) :: wind_y(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(in) :: wind_depth
    logical,          intent(in) :: RelativeWind
    double precision, intent(in) :: Cd
    double precision, intent(in) :: spongeHTimeScale(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: spongeH(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: spongeUTimeScale(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: spongeU(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: spongeVTimeScale(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: spongeV(1-OL:nx+OL, 1-OL:ny+OL, layers)
    integer, intent(in) :: nx, ny, layers, OL
    integer, intent(in) :: n
    integer, intent(in) :: debug_level

    ! Bernoulli potential
    double precision :: b(1-OL:nx+OL, 1-OL:ny+OL, layers)
    ! Relative vorticity
    double precision :: zeta(1-OL:nx+OL, 1-OL:ny+OL, layers)

    ! Calculate Bernoulli potential
    if (RedGrav) then
      call evaluate_b_RedGrav(b, h, u, v, nx, ny, layers, OL, g_vec)
      if (debug_level .ge. 4) then
        call write_output_3d(b, nx, ny, layers, OL, 0, 0, &
          n, 'output/snap.BP.')
      end if
    else
      call evaluate_b_iso(b, h, u, v, nx, ny, layers, OL, g_vec, depth)
      if (debug_level .ge. 4) then
        call write_output_3d(b, nx, ny, layers, OL, 0, 0, &
          n, 'output/snap.BP.')
      end if
    end if

    ! Calculate relative vorticity
    call evaluate_zeta(zeta, u, v, nx, ny, layers, OL, dx, dy)
    if (debug_level .ge. 4) then
      call write_output_3d(zeta, nx, ny, layers, OL, 1, 1, &
        n, 'output/snap.zeta.')
    end if

    ! Calculate dhdt, dudt, dvdt at current time step
    call evaluate_dhdt(dhdt, h, u, v, kh, hmin, kv, dx, dy, nx, ny, &
      layers, OL, spongeHTimeScale, spongeH, wetmask, RedGrav, hAdvecScheme, n)

    call evaluate_dudt(dudt, h, u, v, b, zeta, wind_x, wind_y, wind_depth, &
        fu, au, ar, slip, dx, dy, hfacN, hfacS, nx, ny, layers, OL, rho0, &
        RelativeWind, Cd, spongeUTimeScale, spongeU, RedGrav, botDrag)

    call evaluate_dvdt(dvdt, h, u, v, b, zeta, wind_x, wind_y, wind_depth, &
        fv, au, ar, slip, dx, dy, hfacW, hfacE, nx, ny, layers, OL, rho0, &
        RelativeWind, Cd, spongeVTimeScale, spongeV, RedGrav, botDrag)

    return
  end subroutine state_derivative

end module state_deriv
