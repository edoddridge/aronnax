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
      dx, dy, wetmask, hfac_w, hfac_e, hfac_n, hfac_s, fu, fv, &
      au, ar, bot_drag, kh, kv, hmin, slip, &
      red_grav, h_advec_scheme, g_vec, rho0, wind_x, wind_y, &
      wind_depth, relative_wind, Cd, &
      sponge_h_time_scale, sponge_h, &
      sponge_u_time_scale, sponge_u, &
      sponge_v_time_scale, sponge_v, &
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

    ! Bernoulli potential
    double precision :: b(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    ! Relative vorticity
    double precision :: zeta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)

    ! Calculate Bernoulli potential
    if (red_grav) then
      call evaluate_b_red_grav(b, h, u, v, g_vec, nx, ny, layers, xlow, xhigh, ylow, yhigh, OL, &
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
      nx, ny, layers, OL, sponge_h_time_scale, sponge_h, wetmask, red_grav, h_advec_scheme, n, &
      ilower, iupper, num_procs, myid)

    call evaluate_dudt(dudt, h, u, v, b, zeta, wind_x, wind_y, wind_depth, &
        fu, au, ar, slip, dx, dy, hfac_n, hfac_s, xlow, xhigh, ylow, yhigh, &
        nx, ny, layers, OL, rho0, &
        relative_wind, Cd, sponge_u_time_scale, sponge_u, red_grav, bot_drag, &
        ilower, iupper, num_procs, myid)

    call evaluate_dvdt(dvdt, h, u, v, b, zeta, wind_x, wind_y, wind_depth, &
        fv, au, ar, slip, dx, dy, hfac_w, hfac_e, xlow, xhigh, ylow, yhigh, &
        nx, ny, layers, OL, rho0, &
        relative_wind, Cd, sponge_v_time_scale, sponge_v, red_grav, bot_drag, &
        ilower, iupper, num_procs, myid)

    return
  end subroutine state_derivative

end module state_deriv
