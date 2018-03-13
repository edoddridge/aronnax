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
      au, ar, botDrag, kh, kv, slip, &
      RedGrav, hAdvecScheme, g_vec, rho0, wind_x, wind_y, &
      RelativeWind, Cd, &
      spongeHTimeScale, spongeH, &
      spongeUTimeScale, spongeU, &
      spongeVTimeScale, spongeV, &
      nx, ny, layers, n, debug_level)
    implicit none

    double precision, intent(out) :: dhdt(0:nx+1, 0:ny+1, layers)
    double precision, intent(out) :: dudt(0:nx+1, 0:ny+1, layers)
    double precision, intent(out) :: dvdt(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: v(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: depth(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: dx, dy
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

    ! Bernoulli potential
    double precision :: b(0:nx+1, 0:ny+1, layers)
    ! Relative vorticity
    double precision :: zeta(0:nx+1, 0:ny+1, layers)

    ! Calculate Bernoulli potential
    if (RedGrav) then
      call evaluate_b_RedGrav(b, h, u, v, nx, ny, layers, g_vec)
      if (debug_level .ge. 4) then
        call write_output_3d(b, nx, ny, layers, 0, 0, &
          n, 'output/snap.BP.')
      end if
    else
      call evaluate_b_iso(b, h, u, v, nx, ny, layers, g_vec, depth)
      if (debug_level .ge. 4) then
        call write_output_3d(b, nx, ny, layers, 0, 0, &
          n, 'output/snap.BP.')
      end if
    end if

    ! Calculate relative vorticity
    call evaluate_zeta(zeta, u, v, nx, ny, layers, dx, dy)
    if (debug_level .ge. 4) then
      call write_output_3d(zeta, nx, ny, layers, 1, 1, &
        n, 'output/snap.zeta.')
    end if

    ! Calculate dhdt, dudt, dvdt at current time step
    call evaluate_dhdt(dhdt, h, u, v, kh, kv, dx, dy, nx, ny, layers, &
        spongeHTimeScale, spongeH, wetmask, RedGrav, hAdvecScheme, n)

    call evaluate_dudt(dudt, h, u, v, b, zeta, wind_x, wind_y, fu, au, ar, slip, &
        dx, dy, hfacN, hfacS, nx, ny, layers, rho0, RelativeWind, Cd, &
        spongeUTimeScale, spongeU, RedGrav, botDrag)

    call evaluate_dvdt(dvdt, h, u, v, b, zeta, wind_x, wind_y, fv, au, ar, slip, &
        dx, dy, hfacW, hfacE, nx, ny, layers, rho0, RelativeWind, Cd, &
        spongeVTimeScale, spongeV, RedGrav, botDrag)

    return
  end subroutine state_derivative

end module state_deriv
