module enforce_thickness

  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Check that the free surface anomaly and layer thicknesses are consistent with the depth field. If they're not, then scale the layer thicnkesses to make them fit.

  subroutine enforce_depth_thickness_consistency(h, eta, depth, &
      freesurfFac, thickness_error, nx, ny, layers, OL)
    implicit none

    double precision, intent(inout) :: h(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: eta(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(in) :: depth(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(in) :: freesurfFac, thickness_error
    integer, intent(in) :: nx, ny, layers, OL

    integer k
    double precision h_norming(1-OL:nx+OL, 1-OL:ny+OL)

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

  ! ---------------------------------------------------------------------------
  !> Ensure that layer heights do not fall below the prescribed minimum

  subroutine enforce_minimum_layer_thickness(hnew, hmin, nx, ny, layers, OL, n)
    implicit none

    double precision, intent(inout) :: hnew(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: hmin
    integer, intent(in) :: nx, ny, layers, OL, n

    integer counter, i, j, k

    counter = 0

    do k = 1, layers
      do j = 1, ny
        do i = 1, nx
          if (hnew(i, j, k) .lt. hmin) then
            hnew(i, j, k) = hmin
            counter = counter + 1
            if (counter .eq. 1) then
              write(17, "(A, I0)") &
                  "Layer thickness dropped below hmin at time step ", n
            end if
          end if
        end do
      end do
    end do
    return
  end subroutine enforce_minimum_layer_thickness

end module enforce_thickness
