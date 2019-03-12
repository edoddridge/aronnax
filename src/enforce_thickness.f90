module enforce_thickness

  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Check that the free surface anomaly and layer thicknesses are consistent with the depth field. If they're not, then scale the layer thicnkesses to make them fit.

  subroutine enforce_depth_thickness_consistency(h, eta, depth, &
      freesurfFac, thickness_error, xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    double precision, intent(inout) :: h(xlow-OL:xhigh+OL, &
                                          ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: eta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: depth(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: freesurfFac, thickness_error
    integer, intent(in) :: xlow, xhigh, ylow, yhigh, layers, OL

    integer k
    double precision h_norming(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)

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
