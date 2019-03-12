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
      nx, ny, layers, ilower, iupper, xlow, xhigh, ylow, yhigh, &
      OL, num_procs, myid, n, debug_level)
    implicit none

    double precision, intent(out) :: dhdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out) :: dudt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out) :: dvdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: depth(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: dx, dy
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

    ! Bernoulli potential
    double precision :: b(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    ! Relative vorticity
    double precision :: zeta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)

    ! Calculate Bernoulli potential
    if (RedGrav) then
      call evaluate_b_RedGrav(b, h, u, v, g_vec, nx, ny, layers, xlow, xhigh, ylow, yhigh, OL, &
                                 ilower, iupper, num_procs, myid)
      if (debug_level .ge. 4) then
        call write_output_3d(b, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 0, &
          n, 'output/snap.BP.', num_procs, myid)
      end if
    else
      call evaluate_b_iso(b, h, u, v, g_vec, depth, &
                          nx, ny, layers, xlow, xhigh, ylow, yhigh, OL, &
                          ilower, iupper, num_procs, myid)
      if (debug_level .ge. 4) then
        call write_output_3d(b, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 0, &
          n, 'output/snap.BP.', num_procs, myid)
      end if
    end if

    ! Calculate relative vorticity
    call evaluate_zeta(zeta, u, v, dx, dy, nx, ny, layers, &
                        xlow, xhigh, ylow, yhigh, OL, &
                        ilower, iupper, num_procs, myid)
    if (debug_level .ge. 4) then
      call write_output_3d(zeta, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 1, 1, &
        n, 'output/snap.zeta.', num_procs, myid)
    end if

    ! Calculate dhdt, dudt, dvdt at current time step
    call evaluate_dhdt(dhdt, h, u, v, kh, hmin, kv, dx, dy, xlow, xhigh, ylow, yhigh, &
      nx, ny, layers, OL, spongeHTimeScale, spongeH, wetmask, RedGrav, hAdvecScheme, n, &
      ilower, iupper, num_procs, myid)

    call evaluate_dudt(dudt, h, u, v, b, zeta, wind_x, wind_y, wind_depth, &
        fu, au, ar, slip, dx, dy, hfacN, hfacS, xlow, xhigh, ylow, yhigh, &
        nx, ny, layers, OL, rho0, &
        RelativeWind, Cd, spongeUTimeScale, spongeU, RedGrav, botDrag, &
        ilower, iupper, num_procs, myid)

    call evaluate_dvdt(dvdt, h, u, v, b, zeta, wind_x, wind_y, wind_depth, &
        fv, au, ar, slip, dx, dy, hfacW, hfacE, xlow, xhigh, ylow, yhigh, &
        nx, ny, layers, OL, rho0, &
        RelativeWind, Cd, spongeVTimeScale, spongeV, RedGrav, botDrag, &
        ilower, iupper, num_procs, myid)

    return
  end subroutine state_derivative

end module state_deriv
