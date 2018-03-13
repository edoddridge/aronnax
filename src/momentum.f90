module momentum
  use boundaries
  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Calculate the tendency of zonal velocity for each of the active layers

  subroutine evaluate_dudt(dudt, h, u, v, b, zeta, wind_x, wind_y, fu, &
      au, ar, slip, dx, dy, hfacN, hfacS, nx, ny, layers, rho0, & 
      RelativeWind, Cd, spongeTimeScale, spongeU, RedGrav, botDrag)
    implicit none

    ! dudt(i, j) is evaluated at the centre of the left edge of the grid
    ! box, the same place as u(i, j).
    double precision, intent(out) :: dudt(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: v(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: b(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: zeta(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: wind_x(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: wind_y(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: fu(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: au, ar, slip, dx, dy
    double precision, intent(in)  :: hfacN(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: hfacS(0:nx+1, 0:ny+1)
    integer, intent(in) :: nx, ny, layers
    double precision, intent(in)  :: rho0
    logical,          intent(in)  :: RelativeWind
    double precision, intent(in)  :: Cd
    double precision, intent(in)  :: spongeTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: spongeU(0:nx+1, 0:ny+1, layers)
    logical, intent(in) :: RedGrav
    double precision, intent(in)  :: botDrag

    integer i, j, k

    dudt = 0d0

    do k = 1, layers
      do j = 1, ny
        do i = 1, nx
          dudt(i,j,k) = au*(u(i+1,j,k)+u(i-1,j,k)-2.0d0*u(i,j,k))/(dx*dx) & ! x-component
              + au*(u(i,j+1,k)+u(i,j-1,k)-2.0d0*u(i,j,k) &
                ! boundary conditions
                + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacN(i,j))*u(i,j,k) &
                + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacS(i,j))*u(i,j,k))/(dy*dy) & ! y-component
                ! Together make the horizontal diffusion term
              + 0.25d0*(fu(i,j)+0.5d0*(zeta(i,j,k)+zeta(i,j+1,k))) &
                *(v(i-1,j,k)+v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k)) & ! vorticity term
              - (b(i,j,k) - b(i-1,j,k))/dx & ! Bernoulli potential term
              + spongeTimeScale(i,j,k)*(spongeU(i,j,k)-u(i,j,k)) ! forced relaxtion in the sponge regions
          if (k .eq. 1) then ! only have wind forcing on the top layer
            ! This will need refining in the event of allowing outcropping.
            ! apply wind forcing
            if (RelativeWind) then 
              dudt(i,j,k) = dudt(i,j,k) + (2d0*Cd* & 
                   (wind_x(i,j) - u(i,j,k))* & 
                sqrt((wind_x(i,j) - u(i,j,k))**2 + &
                     (wind_y(i,j) - v(i,j,k))**2))/((h(i,j,k) + h(i-1,j,k)))
            else 
              dudt(i,j,k) = dudt(i,j,k) + 2d0*wind_x(i,j)/(rho0*(h(i,j,k) + h(i-1,j,k))) 
            end if
          end if
          if (layers .gt. 1) then ! only evaluate vertical momentum diffusivity if more than 1 layer
            if (k .eq. 1) then ! adapt vertical momentum diffusivity for 2+ layer model -> top layer
              dudt(i,j,k) = dudt(i,j,k) - 1.0d0*ar*(u(i,j,k) - 1.0d0*u(i,j,k+1))
            else if (k .eq. layers) then ! bottom layer
              dudt(i,j,k) = dudt(i,j,k) - 1.0d0*ar*(u(i,j,k) - 1.0d0*u(i,j,k-1))
              if (.not. RedGrav) then
                ! add bottom drag here in isopycnal version
                dudt(i,j,k) = dudt(i,j,k) - 1.0d0*botDrag*(u(i,j,k))
              end if
            else ! mid layer/s
              dudt(i,j,k) = dudt(i,j,k) - &
                  1.0d0*ar*(2.0d0*u(i,j,k) - 1.0d0*u(i,j,k-1) - 1.0d0*u(i,j,k+1))
            end if
          end if
        end do
      end do
    end do

    call wrap_fields_3D(dudt, nx, ny, layers)

    return
  end subroutine evaluate_dudt

  ! ---------------------------------------------------------------------------
  !> Calculate the tendency of meridional velocity for each of the
  !> active layers

  subroutine evaluate_dvdt(dvdt, h, u, v, b, zeta, wind_x, wind_y, fv, &
      au, ar, slip, dx, dy, hfacW, hfacE, nx, ny, layers, rho0, &
      RelativeWind, Cd, spongeTimeScale, spongeV, RedGrav, botDrag)
    implicit none

    ! dvdt(i, j) is evaluated at the centre of the bottom edge of the
    ! grid box, the same place as v(i, j)
    double precision, intent(out) :: dvdt(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: v(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: b(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: zeta(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: wind_x(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: wind_y(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: fv(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: au, ar, slip
    double precision, intent(in)  :: dx, dy
    double precision, intent(in)  :: hfacW(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: hfacE(0:nx+1, 0:ny+1)
    integer, intent(in) :: nx, ny, layers
    double precision, intent(in)  :: rho0
    logical,          intent(in)  :: RelativeWind
    double precision, intent(in)  :: Cd
    double precision, intent(in)  :: spongeTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: spongeV(0:nx+1, 0:ny+1, layers)
    logical, intent(in) :: RedGrav
    double precision, intent(in)  :: botDrag

    integer i, j, k

    dvdt = 0d0

    do k = 1, layers
      do j = 1, ny
        do i = 1, nx
          dvdt(i,j,k) = &
              au*(v(i+1,j,k)+v(i-1,j,k)-2.0d0*v(i,j,k) &
                ! boundary conditions
                + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacW(i,j))*v(i,j,k) &
                + (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacE(i,j))*v(i,j,k))/(dx*dx) & !x-component
              + au*(v(i,j+1,k) + v(i,j-1,k) - 2.0d0*v(i,j,k))/(dy*dy) & ! y-component.
              ! Together these make the horizontal diffusion term
              - 0.25d0*(fv(i,j)+0.5d0*(zeta(i,j,k)+zeta(i+1,j,k))) &
                *(u(i,j-1,k)+u(i,j,k)+u(i+1,j-1,k)+u(i+1,j,k)) & !vorticity term
              - (b(i,j,k)-b(i,j-1,k))/dy & ! Bernoulli Potential term
              + spongeTimeScale(i,j,k)*(spongeV(i,j,k)-v(i,j,k)) ! forced relaxtion to vsponge (in the sponge regions)
          if (k .eq. 1) then ! only have wind forcing on the top layer
            ! This will need refining in the event of allowing outcropping.
            ! apply wind forcing
            if (RelativeWind) then 
              dvdt(i,j,k) = dvdt(i,j,k) + (2d0*Cd* & 
                   (wind_y(i,j) - v(i,j,k))* & 
                sqrt((wind_x(i,j) - u(i,j,k))**2 + &
                     (wind_y(i,j) - v(i,j,k))**2))/((h(i,j,k) + h(i,j-1,k)))
            else 
              dvdt(i,j,k) = dvdt(i,j,k) + 2d0*wind_y(i,j)/(rho0*(h(i,j,k) + h(i,j-1,k))) 
            end if
          end if
          if (layers .gt. 1) then ! only evaluate vertical momentum diffusivity if more than 1 layer
            if (k .eq. 1) then ! adapt vertical momentum diffusivity for 2+ layer model -> top layer
              dvdt(i,j,k) = dvdt(i,j,k) - 1.0d0*ar*(v(i,j,k) - 1.0d0*v(i,j,k+1))
            else if (k .eq. layers) then ! bottom layer
              dvdt(i,j,k) = dvdt(i,j,k) - 1.0d0*ar*(v(i,j,k) - 1.0d0*v(i,j,k-1))
              if (.not. RedGrav) then
                ! add bottom drag here in isopycnal version
                dvdt(i,j,k) = dvdt(i,j,k) - 1.0d0*botDrag*(v(i,j,k))
              end if
            else ! mid layer/s
              dvdt(i,j,k) = dvdt(i,j,k) - &
                  1.0d0*ar*(2.0d0*v(i,j,k) - 1.0d0*v(i,j,k-1) - 1.0d0*v(i,j,k+1))
            end if
          end if
        end do
      end do
    end do

    call wrap_fields_3D(dvdt, nx, ny, layers)

    return
  end subroutine evaluate_dvdt

end module momentum