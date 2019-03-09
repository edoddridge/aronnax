module advection_schemes

  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Use a first-order centered advection scheme to calculate the advctive
  !! thickness tendency

  subroutine h_advec_1_centered(dhdt_advec, h, u, v, dx, dy, nx, ny, &
                                layers, OL)
    implicit none

    ! dhdt is evaluated at the centre of the grid box
    double precision, intent(out) :: dhdt_advec(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: v(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: dx, dy
    integer, intent(in) :: nx, ny, layers, OL

    integer i, j, k

    dhdt_advec = 0d0

    do k = 1, layers
      do j = 1, ny
        do i = 1, nx
          dhdt_advec(i,j,k) = & 
              - ((h(i,j,k)+h(i+1,j,k))*u(i+1,j,k) &
               - (h(i-1,j,k)+h(i,j,k))*u(i,j,k))/(dx*2d0) & ! d(hu)/dx
              - ((h(i,j,k)+h(i,j+1,k))*v(i,j+1,k) &
               - (h(i,j-1,k)+h(i,j,k))*v(i,j,k))/(dy*2d0)   ! d(hv)/dy
        end do
      end do
    end do

    return
  end subroutine h_advec_1_centered

  ! ---------------------------------------------------------------------------
  !> Use a first-order upwind advection scheme to calculate the advctive
  !! thickness tendency

  subroutine h_advec_1_upwind(dhdt_advec, h, u, v, dx, dy, nx, ny, layers, OL)
    implicit none

    ! dhdt is evaluated at the centre of the grid box
    double precision, intent(out) :: dhdt_advec(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: v(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: dx, dy
    integer, intent(in) :: nx, ny, layers, OL

    integer i, j, k

    dhdt_advec = 0d0

    do k = 1, layers
      do j = 1, ny
        do i = 1, nx
          dhdt_advec(i,j,k) = &
              ! d(hu)/dx
              ! u(i+1,j,k) point
              - ( h(i+1,j,k)*(u(i+1,j,k) - abs(u(i+1,j,k)))/2d0 & ! u < 0
                + h(i,j,k)*(u(i+1,j,k)   + abs(u(i+1,j,k)))/2d0 & ! u > 0
                ! u(i,j,k) point
                - (h(i,j,k)*(u(i,j,k)   - abs(u(i,j,k)))/2d0 & ! u < 0
                  + h(i-1,j,k)*(u(i,j,k) + abs(u(i,j,k)))/2d0) & ! u > 0
                )/dx &

              ! d(hv)/dx
              ! v(i,j+1,k) point
              - ( h(i,j+1,k)*(v(i,j+1,k) - abs(v(i,j+1,k)))/2d0 & ! v < 0
                + h(i,j,k)*(v(i,j+1,k)   + abs(v(i,j+1,k)))/2d0 & ! v > 0
                ! v(i,j,k) point
                - ( h(i,j,k)*(v(i,j,k)   - abs(v(i,j,k)))/2d0 & ! v < 0
                  + h(i,j-1,k)*(v(i,j,k) + abs(v(i,j,k)))/2d0) & ! v > 0
                )/dy
        end do
      end do
    end do

    return
  end subroutine h_advec_1_upwind

end module advection_schemes
