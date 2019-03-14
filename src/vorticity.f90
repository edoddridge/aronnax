module vorticity
  use boundaries
  use exchange

  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Evaluate relative vorticity at lower left grid boundary (du/dy
  !! and dv/dx are at lower left corner as well)
  subroutine evaluate_zeta(zeta, u, v, dx, dy, nx, ny, layers, &
                        xlow, xhigh, ylow, yhigh, OL, &
                        ilower, iupper, num_procs, myid)
    implicit none

    double precision, intent(out) :: zeta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: dx, dy
    integer,          intent(in)  :: nx, ny, layers
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, OL
    integer,          intent(in)  :: ilower(0:num_procs-1,2)
    integer,          intent(in)  :: iupper(0:num_procs-1,2)
    integer,          intent(in)  :: num_procs, myid
    integer i, j, k

    zeta = 0d0

    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL
        do i = xlow-OL+1, xhigh+OL
          zeta(i,j,k) = (v(i,j,k)-v(i-1,j,k))/dx-(u(i,j,k)-u(i,j-1,k))/dy
        end do
      end do
    end do

    return
  end subroutine evaluate_zeta

end module vorticity
