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

  subroutine calc_boundary_masks(wet_mask_file, hfac_w, hfac_e, hfac_s, hfac_n, nx, ny, OL, &
                                  xlow, xhigh, ylow, yhigh)
    implicit none

    character(60), intent(in)     :: wet_mask_file
    double precision, intent(out) :: hfac_w(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(out) :: hfac_e(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(out) :: hfac_n(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(out) :: hfac_s(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    integer, intent(in) :: nx !< number of grid points in x direction
    integer, intent(in) :: ny !< number of grid points in y direction
    integer, intent(in) :: OL !< size of halo region
    integer, intent(in) :: xlow, xhigh, ylow, yhigh ! tile limits

    double precision temp(1-OL:nx+OL, 1-OL:ny+OL)
    double precision wetmask_global(1-OL:nx+OL, 1-OL:ny+OL)
    double precision wetmask_global_small(nx, ny)
    double precision hfac_w_global(1-OL:nx+OL, 1-OL:ny+OL)
    double precision hfac_e_global(1-OL:nx+OL, 1-OL:ny+OL)
    double precision hfac_n_global(1-OL:nx+OL, 1-OL:ny+OL)
    double precision hfac_s_global(1-OL:nx+OL, 1-OL:ny+OL)
    integer i, j

    if (wet_mask_file.ne.'') then
      open(unit=10, form='unformatted', file=wet_mask_file)
      read(10) wetmask_global_small
      close(10)
      wetmask_global(1:nx, 1:ny) = wetmask_global_small
      call wrap_fields_2D(wetmask_global, nx, ny, OL)
    else
      wetmask_global = 1d0
    end if

    hfac_w_global = 1d0

    temp = 0d0
    do j = 0, ny+1
      do i = 1, nx+1
        temp(i, j) = wetmask_global(i-1, j) - wetmask_global(i, j)
      end do
    end do

    do j = 0, ny+1
      do i = 1, nx+1
        if (temp(i, j) .ne. 0.0) then
          hfac_w_global(i, j) = 0d0
        end if
      end do
    end do

    ! and now for all  western cells
    hfac_w_global(0, :) = hfac_w_global(nx, :)

    hfac_e_global = 1d0

    temp = 0d0
    do j = 0, ny+1
      do i = 0, nx
        temp(i, j) = wetmask_global(i, j) - wetmask_global(i+1, j)
      end do
    end do

    do j = 0, ny+1
      do i = 0, nx
        if (temp(i, j) .ne. 0.0) then
          hfac_e_global(i, j) = 0d0
        end if
      end do
    end do

    ! and now for all  eastern cells
    hfac_e_global(nx+1, :) = hfac_e_global(1, :)

    hfac_s_global = 1d0

    temp = 0d0
    do j = 1, ny+1
      do i = 0, nx+1
        temp(i, j) = wetmask_global(i, j-1) - wetmask_global(i, j)
      end do
    end do

    do j = 1, ny+1
      do i = 0, nx+1
        if (temp(i, j) .ne. 0.0) then
          hfac_s_global(i, j) = 0d0
        end if
      end do
    end do

    ! all southern cells
    hfac_s_global(:, 0) = hfac_s_global(:, ny)

    hfac_n_global = 1d0
    temp = 0d0
    do j = 0, ny
      do i = 0, nx+1
        temp(i, j) = wetmask_global(i, j) - wetmask_global(i, j+1)
      end do
    end do

    do j = 0, ny
      do i = 0, nx+1
        if (temp(i, j) .ne. 0.0) then
          hfac_n_global(i, j) = 0d0
        end if
      end do
    end do
    ! all northern cells
    hfac_n_global(:, ny+1) = hfac_n_global(:, 1)

    ! set tile sized arrays
    do j = ylow-OL, yhigh+OL
      do i = xlow-OL, xhigh+OL
        hfac_w(i,j) = hfac_w_global(i,j)
        hfac_e(i,j) = hfac_e_global(i,j)
        hfac_s(i,j) = hfac_s_global(i,j)
        hfac_n(i,j) = hfac_n_global(i,j)
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
    ! - hfac_w and hfac_s are zero where the transition between
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

    integer :: k

    do k = 1,layers
      ! wrap array around for periodicity
      array(1-OL:0, :, k) = array(nx-OL+1:nx, :, k)
      array(nx+1:nx+OL, :, k) = array(1:OL, :, k)
      array(:, 1-OL:0, k) = array(:, ny-OL+1:ny, k)
      array(:, ny+1:ny+OL, k) = array(:, 1:OL, k)
    end do

    return
  end subroutine wrap_fields_3D

  !-----------------------------------------------------------------
  !> Wrap 2D fields around for periodic boundary conditions

  subroutine wrap_fields_2D(array, nx, ny, OL)
    implicit none

    double precision, intent(inout) :: array(1-OL:nx+OL, 1-OL:ny+OL)
    integer, intent(in) :: nx, ny, OL

    ! wrap array around for periodicity
    array(1-OL:0, :) = array(nx-OL+1:nx, :)
    array(nx+1:nx+OL, :) = array(1:OL, :)
    array(:, 1-OL:0) = array(:, ny-OL+1:ny)
    array(:, ny+1:ny+OL) = array(:, 1:OL)

    return
  end subroutine wrap_fields_2D

end module boundaries
