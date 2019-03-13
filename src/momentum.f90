module momentum
  use boundaries
  use end_run
  use exchange
  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Calculate the tendency of zonal velocity for each of the active layers

  subroutine evaluate_dudt(dudt, h, u, v, b, zeta, wind_x, wind_y, &
      wind_depth, fu, au, ar, slip, dx, dy, hfacN, hfacS, xlow, xhigh, ylow, yhigh, &
      nx, ny, layers, &
      OL, rho0, RelativeWind, Cd, spongeTimeScale, spongeU, RedGrav, botDrag, &
      ilower, iupper, num_procs, myid)
    implicit none

    ! dudt(i, j) is evaluated at the centre of the left edge of the grid
    ! box, the same place as u(i, j).
    double precision, intent(out) :: dudt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: b(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: zeta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: wind_depth
    double precision, intent(in)  :: fu(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: au, ar, slip, dx, dy
    double precision, intent(in)  :: hfacN(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: hfacS(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh
    integer,          intent(in)  :: nx, ny, layers, OL
    double precision, intent(in)  :: rho0
    logical,          intent(in)  :: RelativeWind
    double precision, intent(in)  :: Cd
    double precision, intent(in)  :: spongeTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: spongeU(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    logical,          intent(in)  :: RedGrav
    double precision, intent(in)  :: botDrag
    integer,          intent(in)  :: ilower(0:num_procs-1,2)
    integer,          intent(in)  :: iupper(0:num_procs-1,2)
    integer,          intent(in)  :: num_procs, myid

    integer          :: i, j, k
    double precision :: dudt_visc(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: dudt_vort(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: dudt_BP(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: dudt_sponge(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: dudt_wind(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: dudt_drag(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)


    dudt = 0d0
    dudt_visc = 0d0
    dudt_vort = 0d0
    dudt_BP = 0d0
    dudt_sponge = 0d0
    dudt_wind = 0d0
    dudt_drag = 0d0

    call evaluate_dudt_visc(dudt_visc, h, u, au, slip, dx, dy, hfacN, &
      hfacS, xlow, xhigh, ylow, yhigh, layers, OL)

    call evaluate_dudt_vort(dudt_vort, v, zeta, fu, xlow, xhigh, ylow, yhigh, layers, OL)

    call evaluate_dudt_BP(dudt_BP, b, dx, xlow, xhigh, ylow, yhigh, layers, OL)

    call evaluate_dudt_sponge(dudt_sponge, u, spongeTimeScale, spongeU, &
      xlow, xhigh, ylow, yhigh, layers, OL)

    call evaluate_dudt_wind(dudt_wind, h, u, v, wind_x, wind_y, &
      wind_depth, xlow, xhigh, ylow, yhigh, layers, OL, rho0, RelativeWind, Cd)

    call evaluate_dudt_drag(dudt_drag, u, ar, xlow, xhigh, ylow, yhigh, layers, OL, RedGrav, &
      botDrag)

    dudt = dudt + dudt_visc + dudt_vort + dudt_BP + dudt_sponge + &
            dudt_wind + dudt_drag

    call update_halos(dudt, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)

    return
  end subroutine evaluate_dudt


  ! ---------------------------------------------------------------------------
  ! Calculate the viscous tendency of zonal velocity for each of the
  !   active layers

  subroutine evaluate_dudt_visc(dudt_visc, h, u, au, slip, dx, dy, hfacN, &
      hfacS, xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    ! dudt_visc(i, j) is evaluated at the centre of the left edge of the grid
    ! box, the same place as u(i, j).
    double precision, intent(out) :: dudt_visc(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: au, slip, dx, dy
    double precision, intent(in)  :: hfacN(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: hfacS(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL

    integer i, j, k

    dudt_visc = 0d0

    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          dudt_visc(i,j,k) = au*(u(i+1,j,k)+u(i-1,j,k)-2.0d0*u(i,j,k))/(dx*dx) & ! x-component
              + au*(u(i,j+1,k)+u(i,j-1,k)-2.0d0*u(i,j,k) &
                ! boundary conditions
                + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacN(i,j))*u(i,j,k) &
                + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacS(i,j))*u(i,j,k))/(dy*dy) ! y-component
                ! Together make the horizontal diffusion term
        end do
      end do
    end do

    return
  end subroutine evaluate_dudt_visc



  ! ---------------------------------------------------------------------------
  ! Calculate the vorticity contribution to the tendency of zonal velocity for
  !   each of the active layers

  subroutine evaluate_dudt_vort(dudt_vort, v, zeta, fu, xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    ! dudt_vort(i, j) is evaluated at the centre of the left edge of the grid
    ! box, the same place as u(i, j).
    double precision, intent(out) :: dudt_vort(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: zeta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: fu(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL

    integer :: i, j, k

    dudt_vort = 0d0


    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          dudt_vort(i,j,k) = 0.25d0*(fu(i,j)+0.5d0*(zeta(i,j,k)+zeta(i,j+1,k))) &
                *(v(i-1,j,k)+v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k)) ! vorticity term
        end do
      end do
    end do

    return
  end subroutine evaluate_dudt_vort


  ! ---------------------------------------------------------------------------
  ! Calculate the Bernoulli potential contribution to the tendency of zonal
  !   velocity for each of the active layers

  subroutine evaluate_dudt_BP(dudt_BP, b, dx, xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    ! dudt_BP(i, j) is evaluated at the centre of the left edge of the grid
    ! box, the same place as u(i, j).
    double precision, intent(out) :: dudt_BP(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: b(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: dx
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL

    integer :: i, j, k

    dudt_BP = 0d0


    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          dudt_BP(i,j,k) = - (b(i,j,k) - b(i-1,j,k))/dx ! Bernoulli potential term
        end do
      end do
    end do

    return
  end subroutine evaluate_dudt_BP

! ---------------------------------------------------------------------------
  ! Calculate the sponge contribution to the tendency of zonal
  !   velocity for each of the active layers

  subroutine evaluate_dudt_sponge(dudt_sponge, u, spongeTimeScale, spongeU, &
      xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    ! dudt_sponge(i, j) is evaluated at the centre of the left edge of the grid
    ! box, the same place as u(i, j).
    double precision, intent(out) :: dudt_sponge(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: spongeTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: spongeU(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL

    integer :: i, j, k

    dudt_sponge = 0d0


    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          dudt_sponge(i,j,k) = + spongeTimeScale(i,j,k)*(spongeU(i,j,k)-u(i,j,k)) ! forced relaxtion in the sponge regions
        end do
      end do
    end do

    return
  end subroutine evaluate_dudt_sponge


! ---------------------------------------------------------------------------
  ! Calculate the wind contribution to the tendency of zonal
  !   velocity for each of the active layers

  subroutine evaluate_dudt_wind(dudt_wind, h, u, v, wind_x, wind_y, &
      wind_depth, xlow, xhigh, ylow, yhigh, layers, OL, rho0, RelativeWind, Cd)
    implicit none

    ! dudt(i, j) is evaluated at the centre of the left edge of the grid
    ! box, the same place as u(i, j).
    double precision, intent(out) :: dudt_wind(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: wind_depth
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL
    double precision, intent(in)  :: rho0
    logical,          intent(in)  :: RelativeWind
    double precision, intent(in)  :: Cd

    integer          :: i, j, k
    double precision :: z(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision :: forc_frac
    double precision :: recip_wind_depth

    dudt_wind = 0d0
    z = 0d0

    if (wind_depth .eq. 0d0) then
      ! momentum forcing acts only on the top layer (k=1), no matter how thin it gets
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          ! apply wind forcing
          if (RelativeWind) then
            dudt_wind(i,j,1) = (2d0*Cd* &
                 (wind_x(i,j) - u(i,j,1))* &
              sqrt((wind_x(i,j) - u(i,j,1))**2 + &
                   (wind_y(i,j) - v(i,j,1))**2))/((h(i,j,1) + h(i-1,j,1)))
          else
            dudt_wind(i,j,1) = 2d0*wind_x(i,j)/(rho0*(h(i,j,1) + h(i-1,j,1)))
          end if
        end do
      end do

    else if (wind_depth .gt. 0d0) then
      ! apply wind forcing to upper `wind_depth` m of the fluid
      recip_wind_depth = 1d0/wind_depth

      do k = 1, layers
        do j = ylow-OL+1, yhigh+OL-1
          do i = xlow-OL+1, xhigh+OL-1

            if (z(i,j) .le. wind_depth) then
              ! at least a portion of this layer is within wind_depth m
              ! of the surface
              if (z(i,j) + (h(i,j,k) + h(i-1,j,k))*0.5d0 .le. wind_depth) then
                ! all of this layer is within wind_depth m of the surface
                forc_frac = (h(i,j,k) + h(i-1,j,k))*0.5d0*recip_wind_depth
              else
                ! this means (z(i,j) + h(i,j,k) .gt. wind_depth)
                ! only a fraction of this layer is within wind_depth m of
                ! the surface
                forc_frac = (wind_depth - z(i,j))*recip_wind_depth
              end if
            else
              ! z(i,j) is greater than wind_depth. Therefore, this layer
              ! does not receive wind forcing.
              forc_frac = 0d0
            end if

            if (RelativeWind) then 
              dudt_wind(i,j,k) = forc_frac*(2d0*Cd* &
                   (wind_x(i,j) - u(i,j,k))* & 
                sqrt((wind_x(i,j) - u(i,j,k))**2 + &
                     (wind_y(i,j) - v(i,j,k))**2))/((h(i,j,k) + h(i-1,j,k)))
            else 
              dudt_wind(i,j,k) = forc_frac*2d0*wind_x(i,j)/(rho0*(h(i,j,k) &
                                  + h(i-1,j,k)))
            end if

            ! accumulate layer thickness for next pass through
            z(i,j) = z(i,j) + h(i,j,k)

          end do
        end do
      end do

    else
      ! wind_depth not set correctly, end program
      call clean_stop(0, .False.)
    end if

    return
  end subroutine evaluate_dudt_wind


! ---------------------------------------------------------------------------
  !> Calculate the tendency of zonal velocity for each of the active layers

  subroutine evaluate_dudt_drag(dudt_drag, u, ar, xlow, xhigh, ylow, yhigh, layers, OL, RedGrav, &
      botDrag)
    implicit none

    ! dudt_drag(i, j) is evaluated at the centre of the left edge of the grid
    ! box, the same place as u(i, j).
    double precision, intent(out) :: dudt_drag(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: ar
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL
    logical,          intent(in)  :: RedGrav
    double precision, intent(in)  :: botDrag

    integer  :: i, j, k

    dudt_drag = 0d0

    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          if (layers .gt. 1) then ! only evaluate vertical momentum diffusivity if more than 1 layer
            if (k .eq. 1) then ! adapt vertical momentum diffusivity for 2+ layer model -> top layer
              dudt_drag(i,j,k) = - 1.0d0*ar*(u(i,j,k) - 1.0d0*u(i,j,k+1))
            else if (k .eq. layers) then ! bottom layer
              dudt_drag(i,j,k) = - 1.0d0*ar*(u(i,j,k) - 1.0d0*u(i,j,k-1))
            else ! mid layer/s
              dudt_drag(i,j,k) = - &
                  1.0d0*ar*(2.0d0*u(i,j,k) - 1.0d0*u(i,j,k-1) - 1.0d0*u(i,j,k+1))
            end if
          end if
          if (k .eq. layers) then ! add bottom drag if not reduced gravity
            if (.not. RedGrav) then
                ! add bottom drag here in isopycnal version
                dudt_drag(i,j,k) = dudt_drag(i,j,k) - 1.0d0*botDrag*(u(i,j,k))
            end if
          end if
        end do
      end do
    end do

    return
  end subroutine evaluate_dudt_drag

  ! ---------------------------------------------------------------------------
  !> Calculate the tendency of meridional velocity for each of the
  !> active layers

  subroutine evaluate_dvdt(dvdt, h, u, v, b, zeta, wind_x, wind_y, &
      wind_depth, fv, au, ar, slip, dx, dy, hfacW, hfacE, xlow, xhigh, ylow, yhigh, &
      nx, ny, layers, &
      OL, rho0, RelativeWind, Cd, spongeTimeScale, spongeV, RedGrav, botDrag, &
      ilower, iupper, num_procs, myid)
    implicit none

    ! dvdt(i, j) is evaluated at the centre of the bottom edge of the
    ! grid box, the same place as v(i, j)
    double precision, intent(out) :: dvdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: b(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: zeta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: wind_depth
    double precision, intent(in)  :: fv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: au, ar, slip
    double precision, intent(in)  :: dx, dy
    double precision, intent(in)  :: hfacW(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: hfacE(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh
    integer,          intent(in)  :: nx, ny, layers, OL
    double precision, intent(in)  :: rho0
    logical,          intent(in)  :: RelativeWind
    double precision, intent(in)  :: Cd
    double precision, intent(in)  :: spongeTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: spongeV(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    logical,          intent(in)  :: RedGrav
    double precision, intent(in)  :: botDrag
    integer,          intent(in)  :: ilower(0:num_procs-1,2)
    integer,          intent(in)  :: iupper(0:num_procs-1,2)
    integer,          intent(in)  :: num_procs, myid

    integer          :: i, j, k
    double precision :: dvdt_visc(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: dvdt_vort(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: dvdt_BP(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: dvdt_sponge(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: dvdt_wind(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision :: dvdt_drag(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)


    dvdt = 0d0
    dvdt_visc = 0d0
    dvdt_vort = 0d0
    dvdt_BP = 0d0
    dvdt_sponge = 0d0
    dvdt_wind = 0d0
    dvdt_drag = 0d0

    call evaluate_dvdt_visc(dvdt_visc, v, au, slip, dx, dy, hfacW, hfacE, &
      xlow, xhigh, ylow, yhigh, layers, OL)

    call evaluate_dvdt_vort(dvdt_vort, u, zeta, fv, xlow, xhigh, ylow, yhigh, layers, OL)

    call evaluate_dvdt_BP(dvdt_BP, b, dy, xlow, xhigh, ylow, yhigh, layers, OL)

    call evaluate_dvdt_sponge(dvdt_sponge, v, xlow, xhigh, ylow, yhigh, layers, OL, &
      spongeTimeScale, spongeV)

    call evaluate_dvdt_wind(dvdt_wind, h, u, v, wind_x, wind_y, &
      wind_depth, xlow, xhigh, ylow, yhigh, layers, OL, rho0, RelativeWind, Cd)

    call evaluate_dvdt_drag(dvdt_drag, v, ar, xlow, xhigh, ylow, yhigh, layers, OL, RedGrav, &
      botDrag)

    dvdt = dvdt + dvdt_visc + dvdt_vort + dvdt_BP + dvdt_sponge + &
            dvdt_wind + dvdt_drag

    call update_halos(dvdt, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)

    return
  end subroutine evaluate_dvdt


  ! ---------------------------------------------------------------------------
  !> Calculate the viscous tendency of meridional velocity for each of the
  !> active layers

  subroutine evaluate_dvdt_visc(dvdt_visc, v, au, slip, dx, dy, hfacW, hfacE, &
      xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    ! dvdt_visc(i, j) is evaluated at the centre of the bottom edge of the
    ! grid box, the same place as v(i, j)
    double precision, intent(out) :: dvdt_visc(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: au, slip
    double precision, intent(in)  :: dx, dy
    double precision, intent(in)  :: hfacW(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: hfacE(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL

    integer :: i, j, k

    dvdt_visc = 0d0

    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          dvdt_visc(i,j,k) = &
              au*(v(i+1,j,k)+v(i-1,j,k)-2.0d0*v(i,j,k) &
                ! boundary conditions
                + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacW(i,j))*v(i,j,k) &
                + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacE(i,j))*v(i,j,k))/(dx*dx) & !x-component
              + au*(v(i,j+1,k) + v(i,j-1,k) - 2.0d0*v(i,j,k))/(dy*dy) ! y-component.
              ! Together these make the horizontal diffusion term
        end do
      end do
    end do

    return
  end subroutine evaluate_dvdt_visc

  ! ---------------------------------------------------------------------------
  ! Calculate the vorticity contribution to the tendency of meridional
  !   velocity for each of the active layers

  subroutine evaluate_dvdt_vort(dvdt_vort, u, zeta, fv, xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    ! dvdt_vort(i, j) is evaluated at the centre of the bottom edge of the
    ! grid box, the same place as v(i, j)
    double precision, intent(out) :: dvdt_vort(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: zeta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: fv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL

    integer :: i, j, k

    dvdt_vort = 0d0

    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          dvdt_vort(i,j,k) = &
              - 0.25d0*(fv(i,j)+0.5d0*(zeta(i,j,k)+zeta(i+1,j,k))) &
                *(u(i,j-1,k)+u(i,j,k)+u(i+1,j-1,k)+u(i+1,j,k)) !vorticity term
        end do
      end do
    end do


    return
  end subroutine evaluate_dvdt_vort


  ! ---------------------------------------------------------------------------
  ! Calculate the Bernoulli potential contribution to the tendency of
  !   meridional velocity for each of the active layers

  subroutine evaluate_dvdt_BP(dvdt_BP, b, dy, xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    ! dvdt_BP(i, j) is evaluated at the centre of the bottom edge of the
    ! grid box, the same place as v(i, j)
    double precision, intent(out) :: dvdt_BP(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: b(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: dy
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL

    integer :: i, j, k

    dvdt_BP = 0d0

    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          dvdt_BP(i,j,k) = - (b(i,j,k)-b(i,j-1,k))/dy ! Bernoulli Potential term
        end do
      end do
    end do

    return
  end subroutine evaluate_dvdt_BP


  ! ---------------------------------------------------------------------------
  ! Calculate the sponge contribution to the tendency of meridional velocity
  !   for each of the active layers

  subroutine evaluate_dvdt_sponge(dvdt_sponge, v, xlow, xhigh, ylow, yhigh, layers, OL, &
      spongeTimeScale, spongeV)
    implicit none

    ! dvdt_sponge(i, j) is evaluated at the centre of the bottom edge of the
    ! grid box, the same place as v(i, j)
    double precision, intent(out) :: dvdt_sponge(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL
    double precision, intent(in)  :: spongeTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: spongeV(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)

    integer :: i, j, k

    dvdt_sponge = 0d0

    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          dvdt_sponge(i,j,k) =  spongeTimeScale(i,j,k)*(spongeV(i,j,k)-v(i,j,k))
          ! forced relaxtion to vsponge (in the sponge regions)
        end do
      end do
    end do

    return
  end subroutine evaluate_dvdt_sponge


  ! ---------------------------------------------------------------------------
  ! Calculate the wind contribution to the tendency of meridional
  !   velocity for each of the active layers

  subroutine evaluate_dvdt_wind(dvdt_wind, h, u, v, wind_x, wind_y, &
      wind_depth, xlow, xhigh, ylow, yhigh, layers, OL, rho0, RelativeWind, Cd)
    implicit none

    ! dvdt_wind(i, j) is evaluated at the centre of the bottom edge of the grid
    ! box, the same place as v(i, j).
    double precision, intent(out) :: dvdt_wind(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: wind_depth
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL
    double precision, intent(in)  :: rho0
    logical,          intent(in)  :: RelativeWind
    double precision, intent(in)  :: Cd

    integer          :: i, j, k
    double precision :: z(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision :: forc_frac
    double precision :: recip_wind_depth

    dvdt_wind = 0d0
    z = 0d0

    if (wind_depth .eq. 0d0) then
      ! momentum forcing acts only on the top layer (k=1), no matter how thin it gets
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          ! apply wind forcing
          if (RelativeWind) then
            dvdt_wind(i,j,1) = (2d0*Cd* &
                 (wind_y(i,j) - v(i,j,1))* &
              sqrt((wind_y(i,j) - v(i,j,1))**2 + &
                   (wind_x(i,j) - u(i,j,1))**2))/((h(i,j,1) + h(i,j-1,1)))
          else
            dvdt_wind(i,j,1) = 2d0*wind_y(i,j)/(rho0*(h(i,j,1) + h(i,j-1,1)))
          end if
        end do
      end do

    else if (wind_depth .gt. 0d0) then
      ! apply wind forcing to upper `wind_depth` m of the fluid
      recip_wind_depth = 1d0/wind_depth

      do k = 1, layers
        do j = ylow-OL+1, yhigh+OL-1
          do i = xlow-OL+1, xhigh+OL-1
            if (z(i,j) .le. wind_depth) then
              ! at least a portion of this layer is within wind_depth m
              ! of the surface
              if (z(i,j) + (h(i,j,k) + h(i,j-1,k))*0.5d0 .le. wind_depth) then
                ! all of this layer is within wind_depth m of the surface
                forc_frac = (h(i,j,k) + h(i,j-1,k))*0.5d0*recip_wind_depth
              else
                ! this means (z(i,j) + h(i,j,k) .gt. wind_depth)
                ! only a fraction of this layer is within wind_depth m of
                ! the surface
                forc_frac = (wind_depth - z(i,j))*recip_wind_depth
              end if
            else
              ! z(i,j) is greater than wind_depth. Therefore, this layer
              ! does not receive wind forcing.
              forc_frac = 0d0
            end if

            if (RelativeWind) then 
              dvdt_wind(i,j,k) = forc_frac*(2d0*Cd* &
                   (wind_y(i,j) - v(i,j,k))* & 
                sqrt((wind_y(i,j) - v(i,j,k))**2 + &
                     (wind_x(i,j) - u(i,j,k))**2))/((h(i,j,k) + h(i,j-1,k)))
            else 
              dvdt_wind(i,j,k) = forc_frac*2d0*wind_y(i,j)/(rho0*(h(i,j,k) &
                                  + h(i,j-1,k)))
            end if

            ! accumulate layer thickness for next pass through
            z(i,j) = z(i,j) + h(i,j,k)

          end do
        end do
      end do

    else
      ! wind_depth not set correctly, end program
      call clean_stop(0, .False.)
    end if

    return
  end subroutine evaluate_dvdt_wind


  ! ---------------------------------------------------------------------------
  ! Calculate the linear drag contribution to the tendency of meridional
  !   velocity for each of the active layers

  subroutine evaluate_dvdt_drag(dvdt_drag, v, ar, xlow, xhigh, ylow, yhigh, layers, OL, RedGrav, &
      botDrag)
    implicit none

    ! dvdt_drag(i, j) is evaluated at the centre of the bottom edge of the
    ! grid box, the same place as v(i, j)
    double precision, intent(out) :: dvdt_drag(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: ar
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, layers, OL
    logical,          intent(in)  :: RedGrav
    double precision, intent(in)  :: botDrag

    integer :: i, j, k

    dvdt_drag = 0d0

    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          if (layers .gt. 1) then ! only evaluate vertical momentum diffusivity if more than 1 layer
            if (k .eq. 1) then ! adapt vertical momentum diffusivity for 2+ layer model -> top layer
              dvdt_drag(i,j,k) =  - 1.0d0*ar*(v(i,j,k) - 1.0d0*v(i,j,k+1))
            else if (k .eq. layers) then ! bottom layer
              dvdt_drag(i,j,k) =  - 1.0d0*ar*(v(i,j,k) - 1.0d0*v(i,j,k-1))
            else ! mid layer/s
              dvdt_drag(i,j,k) =  - &
                  1.0d0*ar*(2.0d0*v(i,j,k) - 1.0d0*v(i,j,k-1) - 1.0d0*v(i,j,k+1))
            end if
          end if
          if (k .eq. layers) then ! add bottom drag if not reduced gravity
            if (.not. RedGrav) then
              ! add bottom drag here in isopycnal version
              dvdt_drag(i,j,k) = dvdt_drag(i,j,k) - 1.0d0*botDrag*(v(i,j,k))
            end if
          end if
        end do
      end do
    end do

    return
  end subroutine evaluate_dvdt_drag

end module momentum
