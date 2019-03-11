module exchange
  use boundaries
  use mpi


  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Update global and tile halos
  !! Use MPI functions to move data between tiles.

  subroutine update_halos(array, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)
    implicit none

    double precision, intent(inout) :: array(xlow-OL:xhigh+OL, &
                                              ylow-OL:yhigh+OL, layers)
    integer,          intent(in) :: nx, ny, layers
    integer,          intent(in) :: ilower(0:num_procs-1,2)
    integer,          intent(in) :: iupper(0:num_procs-1,2)
    integer,          intent(in) :: xlow, xhigh, ylow, yhigh, OL
    integer,          intent(in) :: num_procs, myid

    integer          :: i, j, k, l
    integer          :: i_local, j_local, k_local
    double precision :: global_array(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision :: local_vector((xhigh-xlow+1)*(yhigh-ylow+1)*layers)
    double precision :: global_vector(nx*ny*layers)
    integer          :: ierr


    do k = 1, layers
      do j = ylow, yhigh
        do i = xlow, xhigh
          i_local = i - xlow
          j_local = j - ylow
          k_local = k - 1
     
     ! the 3D array is being laid out like
     ! [x1y1z1, x2y1z1, x3y1z1, x1y2z1, x2y2z1, x3y2z1, x1y3z1, x2y3z1, x3y3z1,
     !  x1y1z2, x2y1z2, x3y1z2, x1y2z2, x2y2z2, x3y2z2, x1y3z2, x2y3z2, x3y3z2]

          local_vector(i_local &
                     + j_local*(xhigh - xlow + 1) &
                     + k_local*(xhigh - xlow + 1)*(yhigh - ylow + 1) &
                     + 1) = array(i,j,k)
        end do
      end do
    end do


    ! bring all local_vectors together into global_vector
    call MPI_ALLGATHER(local_vector, (xhigh-xlow+1)*(yhigh-ylow+1)*layers, &
                      MPI_DOUBLE_PRECISION, &
                    global_vector, (xhigh-xlow+1)*(yhigh-ylow+1)*layers, &
                      MPI_DOUBLE_PRECISION, &
                      MPI_COMM_WORLD, ierr)

      ! lay global_vector out into global_array

    do l = 0, num_procs-1
      do k = 1, layers
        do j = ilower(l,2), iupper(l,2)
          do i = ilower(l,1), iupper(l,1)
            i_local = i - ilower(l,1)
            j_local = j - ilower(l,2)
            k_local = k - 1

            ! the 1D vector is laid out like
            ! [x1y1z1, x2y1z1, x3y1z1, x1y2z1, x2y2z1, x3y2z1, x1y3z1, x2y3z1, x3y3z1,
            !  x1y1z2, x2y1z2, x3y1z2, x1y2z2, x2y2z2, x3y2z2, x1y3z2, x2y3z2, x3y3z2,
            !  NEXT TILE
            !  x1y1z1, x2y1z1, x3y1z1, x1y2z1, x2y2z1, x3y2z1, x1y3z1, x2y3z1, x3y3z1,
            !  x1y1z2, x2y1z2, x3y1z2, x1y2z2, x2y2z2, x3y2z2, x1y3z2, x2y3z2, x3y3z2]
            ! where numbers refer to local tile indicies

            global_array(i,j,k) = global_vector(i_local &
                                + j_local*(xhigh - xlow + 1) &
                                + k_local*(xhigh - xlow + 1)*(yhigh - ylow + 1) &
                                + l*(xhigh - xlow + 1)*(yhigh - ylow + 1)*layers &
                                + 1)
          end do
        end do
      end do
    end do

    ! update global halos
    if (layers .eq. 1) then
      call wrap_fields_2D(global_array(:,:,1), nx, ny, OL)
    else
      call wrap_fields_3D(global_array, nx, ny, layers, OL)
    end if

    do k = 1, layers
      do j = ilower(myid,2) - OL, iupper(myid,2) + OL
        do i = ilower(myid,1) - OL, iupper(myid,1) + OL
          array(i,j,k) = global_array(i,j,k)
        end do
      end do
    end do

    return
  end subroutine update_halos

end module exchange
