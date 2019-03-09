module boundaries

  implicit none

  contains

  !----------------------------------------------------------------------------
  !> Define masks for boundary conditions in u and v.
  !! This finds locations where neighbouring grid boxes are not the same
  !! (i.e. one is land and one is ocean).
  !! In the output,
  !! 0 means barrier
  !! 1 mean open

  subroutine calc_boundary_masks(wetmask, hfacW, hfacE, hfacS, hfacN, nx, ny, OL)
    implicit none

    double precision, intent(in)  :: wetmask(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(out) :: hfacW(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(out) :: hfacE(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(out) :: hfacN(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(out) :: hfacS(1-OL:nx+OL, 1-OL:ny+OL)
    integer, intent(in) :: nx !< number of grid points in x direction
    integer, intent(in) :: ny !< number of grid points in y direction
    integer, intent(in) :: OL !< size of halo region

    double precision temp(1-OL:nx+OL, 1-OL:ny+OL)
    integer i, j

    hfacW = 1d0

    temp = 0.0
    do j = 0, ny+1
      do i = 1, nx+1
        temp(i, j) = wetmask(i-1, j) - wetmask(i, j)
      end do
    end do

    do j = 0, ny+1
      do i = 1, nx+1
        if (temp(i, j) .ne. 0.0) then
          hfacW(i, j) = 0d0
        end if
      end do
    end do

    ! and now for all  western cells
    hfacW(0, :) = hfacW(nx, :)

    hfacE = 1d0

    temp = 0.0
    do j = 0, ny+1
      do i = 0, nx
        temp(i, j) = wetmask(i, j) - wetmask(i+1, j)
      end do
    end do

    do j = 0, ny+1
      do i = 0, nx
        if (temp(i, j) .ne. 0.0) then
          hfacE(i, j) = 0d0
        end if
      end do
    end do

    ! and now for all  eastern cells
    hfacE(nx+1, :) = hfacE(1, :)

    hfacS = 1

    temp = 0.0
    do j = 1, ny+1
      do i = 0, nx+1
        temp(i, j) = wetmask(i, j-1) - wetmask(i, j)
      end do
    end do

    do j = 1, ny+1
      do i = 0, nx+1
        if (temp(i, j) .ne. 0.0) then
          hfacS(i, j) = 0d0
        end if
      end do
    end do

    ! all southern cells
    hfacS(:, 0) = hfacS(:, ny)

    hfacN = 1
    temp = 0.0
    do j = 0, ny
      do i = 0, nx+1
        temp(i, j) = wetmask(i, j) - wetmask(i, j+1)
      end do
    end do

    do j = 0, ny
      do i = 0, nx+1
        if (temp(i, j) .ne. 0.0) then
          hfacN(i, j) = 0d0
        end if
      end do
    end do
    ! all northern cells
    hfacN(:, ny+1) = hfacN(:, 1)

    return
  end subroutine calc_boundary_masks

  ! ---------------------------------------------------------------------------
  !> Apply the boundary conditions

  subroutine apply_boundary_conditions(array, hfac, wetmask, nx, ny, layers, OL)
    implicit none

    double precision, intent(inout) :: array(1-OL:nx+OL,1-OL:ny+OL,layers)
    double precision, intent(in) :: hfac(1-OL:nx+OL,1-OL:ny+OL)
    double precision, intent(in) :: wetmask(1-OL:nx+OL,1-OL:ny+OL)
    integer, intent(in) :: nx, ny, layers, OL

    integer k

    ! - Enforce no normal flow boundary condition
    !   and no flow in dry cells.
    ! - no/free-slip is done inside the dudt and dvdt subroutines.
    ! - hfacW and hfacS are zero where the transition between
    !   wet and dry cells occurs.
    ! - wetmask is 1 in wet cells, and zero in dry cells.

    do k = 1, layers
      array(:, :, k) = array(:, :, k) * hfac * wetmask(:, :)
    end do

    return
  end subroutine apply_boundary_conditions

  !-----------------------------------------------------------------
  !> Wrap 3D fields around for periodic boundary conditions

  subroutine wrap_fields_3D(array, nx, ny, layers, OL)
    implicit none

    double precision, intent(inout) :: array(1-OL:nx+OL, 1-OL:ny+OL, layers)
    integer, intent(in) :: nx, ny, layers, OL

    ! wrap array around for periodicity
    array(0, :, :) = array(nx, :, :)
    array(nx+1, :, :) = array(1, :, :)
    array(:, 0, :) = array(:, ny, :)
    array(:, ny+1, :) = array(:, 1, :)

    return
  end subroutine wrap_fields_3D

  !-----------------------------------------------------------------
  !> Wrap 2D fields around for periodic boundary conditions

  subroutine wrap_fields_2D(array, nx, ny, OL)
    implicit none

    double precision, intent(inout) :: array(1-OL:nx+OL, 1-OL:ny+OL)
    integer, intent(in) :: nx, ny, OL

    ! wrap array around for periodicity
    array(0, :) = array(nx, :)
    array(nx+1, :) = array(1, :)
    array(:, 0) = array(:, ny)
    array(:, ny+1) = array(:, 1)

    return
  end subroutine wrap_fields_2D

end module boundaries
