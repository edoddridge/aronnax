module vorticity
  use boundaries


  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Evaluate relative vorticity at lower left grid boundary (du/dy
  !! and dv/dx are at lower left corner as well)
  subroutine evaluate_zeta(zeta, u, v, nx, ny, layers, dx, dy)
    implicit none

    double precision, intent(out) :: zeta(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: v(0:nx+1, 0:ny+1, layers)
    integer, intent(in) :: nx, ny, layers
    double precision, intent(in)  :: dx, dy

    integer i, j, k

    zeta = 0d0

    do k = 1, layers
      do j = 1, ny+1
        do i = 1, nx+1
          zeta(i,j,k) = (v(i,j,k)-v(i-1,j,k))/dx-(u(i,j,k)-u(i,j-1,k))/dy
        end do
      end do
    end do

    call wrap_fields_3D(zeta, nx, ny, layers)

    return
  end subroutine evaluate_zeta

end module vorticity
