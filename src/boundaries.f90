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

  subroutine calc_boundary_masks(wetMaskFile, hfacW, hfacE, hfacS, hfacN, nx, ny, OL, &
                                  xlow, xhigh, ylow, yhigh)
    implicit none

    character(60), intent(in)     :: wetMaskFile
    double precision, intent(out) :: hfacW(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(out) :: hfacE(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(out) :: hfacN(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(out) :: hfacS(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    integer, intent(in) :: nx !< number of grid points in x direction
    integer, intent(in) :: ny !< number of grid points in y direction
    integer, intent(in) :: OL !< size of halo region
    integer, intent(in) :: xlow, xhigh, ylow, yhigh ! tile limits

    double precision temp(1-OL:nx+OL, 1-OL:ny+OL)
    double precision wetmask_global(1-OL:nx+OL, 1-OL:ny+OL)
    double precision wetmask_global_small(nx, ny)
    double precision hfacW_global(1-OL:nx+OL, 1-OL:ny+OL)
    double precision hfacE_global(1-OL:nx+OL, 1-OL:ny+OL)
    double precision hfacN_global(1-OL:nx+OL, 1-OL:ny+OL)
    double precision hfacS_global(1-OL:nx+OL, 1-OL:ny+OL)
    integer i, j

    if (wetMaskFile.ne.'') then
      open(unit=10, form='unformatted', file=wetMaskFile)
      read(10) wetmask_global_small
      close(10)
      wetmask_global(1:nx, 1:ny) = wetmask_global_small
      call wrap_fields_2D(wetmask_global, nx, ny, OL)
    else
      wetmask_global = 1d0
    end if

    hfacW_global = 1d0

    temp = 0d0
    do j = 0, ny+1
      do i = 1, nx+1
        temp(i, j) = wetmask_global(i-1, j) - wetmask_global(i, j)
      end do
    end do

    do j = 0, ny+1
      do i = 1, nx+1
        if (temp(i, j) .ne. 0.0) then
          hfacW_global(i, j) = 0d0
        end if
      end do
    end do

    ! and now for all  western cells
    hfacW_global(0, :) = hfacW_global(nx, :)

    hfacE_global = 1d0

    temp = 0d0
    do j = 0, ny+1
      do i = 0, nx
        temp(i, j) = wetmask_global(i, j) - wetmask_global(i+1, j)
      end do
    end do

    do j = 0, ny+1
      do i = 0, nx
        if (temp(i, j) .ne. 0.0) then
          hfacE_global(i, j) = 0d0
        end if
      end do
    end do

    ! and now for all  eastern cells
    hfacE_global(nx+1, :) = hfacE_global(1, :)

    hfacS_global = 1d0

    temp = 0d0
    do j = 1, ny+1
      do i = 0, nx+1
        temp(i, j) = wetmask_global(i, j-1) - wetmask_global(i, j)
      end do
    end do

    do j = 1, ny+1
      do i = 0, nx+1
        if (temp(i, j) .ne. 0.0) then
          hfacS_global(i, j) = 0d0
        end if
      end do
    end do

    ! all southern cells
    hfacS_global(:, 0) = hfacS_global(:, ny)

    hfacN_global = 1d0
    temp = 0d0
    do j = 0, ny
      do i = 0, nx+1
        temp(i, j) = wetmask_global(i, j) - wetmask_global(i, j+1)
      end do
    end do

    do j = 0, ny
      do i = 0, nx+1
        if (temp(i, j) .ne. 0.0) then
          hfacN_global(i, j) = 0d0
        end if
      end do
    end do
    ! all northern cells
    hfacN_global(:, ny+1) = hfacN_global(:, 1)

    ! set tile sized arrays
    do j = ylow-OL, yhigh+OL
      do i = xlow-OL, xhigh+OL
        hfacW(i,j) = hfacW_global(i,j)
        hfacE(i,j) = hfacE_global(i,j)
        hfacS(i,j) = hfacS_global(i,j)
        hfacN(i,j) = hfacN_global(i,j)
      end do
    end do

    return
  end subroutine calc_boundary_masks

  ! ---------------------------------------------------------------------------
  !> Apply the boundary conditions

  subroutine apply_boundary_conditions(array, hfac, wetmask, &
                              xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    double precision, intent(inout) :: array(xlow-OL:xhigh+OL, &
                                              ylow-OL:yhigh+OL,layers)
    double precision, intent(in) :: hfac(xlow-OL:xhigh+OL, &
                                          ylow-OL:yhigh+OL)
    double precision, intent(in) :: wetmask(xlow-OL:xhigh+OL, &
                                            ylow-OL:yhigh+OL)
    integer, intent(in) :: xlow, xhigh, ylow, yhigh, layers, OL

    integer k

    ! - Enforce no normal flow boundary condition
    !   and no flow in dry cells.
    ! - no/free-slip is done inside the dudt and dvdt subroutines.
    ! - hfacW and hfacS are zero where the transition between
    !   wet and dry cells occurs.
    ! - wetmask is 1 in wet cells, and zero in dry cells.

    do k = 1, layers
      array(:, :, k) = array(:, :, k) * hfac(:,:) * wetmask(:, :)
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
