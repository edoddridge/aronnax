module end_run
  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Check to see if there are any NaNs in the data field and stop the
  !! calculation if any are found.

  subroutine break_if_NaN(data, xlow, xhigh, ylow, yhigh, layers, OL, n)
    implicit none

    ! To stop the program if it detects a NaN in the variable being checked

    double precision, intent(in) :: data(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    integer, intent(in) :: xlow, xhigh, ylow, yhigh, layers, OL, n

    integer :: i, j, k

    do k = 1, layers
      do j = ylow, yhigh
        do i = xlow, xhigh
          if (data(i,j,k) .ne. data(i,j,k)) then
            write(17, "(A, I0)") "NaN detected at time step ", n
            call clean_stop(n, .FALSE.)
          end if
        end do
      end do
    end do

    return
  end subroutine break_if_NaN

  !-----------------------------------------------------------------
  !> finalise MPI and then stop the model

  subroutine clean_stop(n, happy)
    implicit none

    integer, intent(in) :: n
    logical, intent(in) :: happy

    integer :: ierr

    if (happy) then
      call MPI_Finalize(ierr)
      stop
    else
      print "(A, I0, A, I0, A)", "Unexpected termination at time step ", n
      call MPI_Finalize(ierr)
      stop 1
    end if

    return
  end subroutine clean_stop

end module end_run
