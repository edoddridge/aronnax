module enforce_thickness

  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Check that the free surface anomaly and layer thicknesses are consistent with the depth field. If they're not, then scale the layer thicnkesses to make them fit.

  subroutine enforce_depth_thickness_consistency(h, eta, depth, &
      freesurfFac, thickness_error, nx, ny, layers)
    implicit none

    double precision, intent(inout) :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: eta(0:nx+1, 0:ny+1)
    double precision, intent(in) :: depth(0:nx+1, 0:ny+1)
    double precision, intent(in) :: freesurfFac, thickness_error
    integer, intent(in) :: nx, ny, layers

    integer k
    double precision h_norming(0:nx+1, 0:ny+1)

    h_norming = (freesurfFac*eta + depth) / sum(h,3)
    do k = 1, layers
      h(:, :, k) = h(:, :, k) * h_norming
    end do

    if (maxval(abs(h_norming - 1d0)) .gt. thickness_error) then
      write(17, "(A, F6.3, A)") 'Inconsistency between h and eta: ', &
          maxval(abs(h_norming - 1d0))*100d0, '%'
    end if

    return
  end subroutine enforce_depth_thickness_consistency

end module enforce_thickness
