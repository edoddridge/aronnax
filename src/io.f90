module io
  use end_run
  use boundaries
  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Write output if it's time

  subroutine maybe_dump_output(h, hav, u, uav, v, vav, eta, etaav, &
          dudt, dvdt, dhdt, &
          dudtold, dvdtold, dhdtold, &
          dudtveryold, dvdtveryold, dhdtveryold, &
          wind_x, wind_y, nx, ny, layers, &
          n, nwrite, avwrite, checkpointwrite, diagwrite, &
          RedGrav, DumpWind, debug_level)
    implicit none

    double precision, intent(in)    :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(inout) :: hav(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(inout) :: uav(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: v(0:nx+1, 0:ny+1, layers)
    double precision, intent(inout) :: vav(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: eta(0:nx+1, 0:ny+1)
    double precision, intent(inout) :: etaav(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: dudt(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: dvdt(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: dhdt(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: dudtold(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: dvdtold(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: dhdtold(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: dudtveryold(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: dvdtveryold(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: dhdtveryold(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: wind_x(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: wind_y(0:nx+1, 0:ny+1)
    integer,          intent(in)    :: nx, ny, layers, n
    integer,          intent(in)    :: nwrite, avwrite, checkpointwrite, diagwrite
    logical,          intent(in)    :: RedGrav, DumpWind
    integer,          intent(in)    :: debug_level

    logical       :: dump_output


    ! Write snapshot to file?
    if (mod(n-1, nwrite) .eq. 0) then
      dump_output = .TRUE.
    else if (debug_level .ge. 4) then
      dump_output = .TRUE.
    else
      dump_output = .FALSE.
    end if

    if (dump_output) then 
      
      call write_output_3d(h, nx, ny, layers, 0, 0, &
      n, 'output/snap.h.')
      call write_output_3d(u, nx, ny, layers, 1, 0, &
      n, 'output/snap.u.')
      call write_output_3d(v, nx, ny, layers, 0, 1, &
      n, 'output/snap.v.')


      if (.not. RedGrav) then
        call write_output_2d(eta, nx, ny, 0, 0, &
          n, 'output/snap.eta.')
      end if

      if (DumpWind .eqv. .true.) then
        call write_output_2d(wind_x, nx, ny, 1, 0, &
          n, 'output/wind_x.')
        call write_output_2d(wind_y, nx, ny, 0, 1, &
          n, 'output/wind_y.')
      end if

      if (debug_level .ge. 1) then
        call write_output_3d(dhdt, nx, ny, layers, 0, 0, &
          n, 'output/debug.dhdt.')
        call write_output_3d(dudt, nx, ny, layers, 1, 0, &
          n, 'output/debug.dudt.')
        call write_output_3d(dvdt, nx, ny, layers, 0, 1, &
          n, 'output/debug.dvdt.')
      end if

      ! Check if there are NaNs in the data
      call break_if_NaN(h, nx, ny, layers, n)
      ! call break_if_NaN(u, nx, ny, layers, n)
      ! call break_if_NaN(v, nx, ny, layers, n)

    end if

    ! Write accumulated averages to file?
    if (avwrite .eq. 0) then
      ! OK
    else if (mod(n-1, avwrite) .eq. 0) then

      if (n .eq. 1) then
        ! pass, since dumping averages after first timestep isn't helpful
      else 
        hav = hav/real(avwrite)
        uav = uav/real(avwrite)
        vav = vav/real(avwrite)
        if (.not. RedGrav) then
          etaav = etaav/real(avwrite)
        end if

        call write_output_3d(hav, nx, ny, layers, 0, 0, &
        n, 'output/av.h.')
        call write_output_3d(uav, nx, ny, layers, 1, 0, &
        n, 'output/av.u.')
        call write_output_3d(vav, nx, ny, layers, 0, 1, &
        n, 'output/av.v.')


        if (.not. RedGrav) then
          call write_output_2d(etaav, nx, ny, 0, 0, &
            n, 'output/av.eta.')
        end if

        ! Check if there are NaNs in the data
        call break_if_NaN(h, nx, ny, layers, n)
        ! call break_if_NaN(u, nx, ny, layers, n)
        ! call break_if_NaN(v, nx, ny, layers, n)
      end if
      
      ! Reset average quantities
      hav = 0.0
      uav = 0.0
      vav = 0.0
      if (.not. RedGrav) then
        etaav = 0.0
      end if
      ! h2av = 0.0

    end if

    ! save a checkpoint?
    if (checkpointwrite .eq. 0) then
      ! not saving checkpoints, so move on
    else if (mod(n-1, checkpointwrite) .eq. 0) then
      call write_checkpoint_output(h, nx, ny, layers, &
      n, 'checkpoints/h.')
      call write_checkpoint_output(u, nx, ny, layers, &
      n, 'checkpoints/u.')
      call write_checkpoint_output(v, nx, ny, layers, &
      n, 'checkpoints/v.')

      call write_checkpoint_output(dhdt, nx, ny, layers, &
        n, 'checkpoints/dhdt.')
      call write_checkpoint_output(dudt, nx, ny, layers, &
        n, 'checkpoints/dudt.')
      call write_checkpoint_output(dvdt, nx, ny, layers, &
        n, 'checkpoints/dvdt.')

      call write_checkpoint_output(dhdtold, nx, ny, layers, &
        n, 'checkpoints/dhdtold.')
      call write_checkpoint_output(dudtold, nx, ny, layers, &
        n, 'checkpoints/dudtold.')
      call write_checkpoint_output(dvdtold, nx, ny, layers, &
        n, 'checkpoints/dvdtold.')

      call write_checkpoint_output(dhdtveryold, nx, ny, layers, &
        n, 'checkpoints/dhdtveryold.')
      call write_checkpoint_output(dudtveryold, nx, ny, layers, &
        n, 'checkpoints/dudtveryold.')
      call write_checkpoint_output(dvdtveryold, nx, ny, layers, &
        n, 'checkpoints/dvdtveryold.')

      if (.not. RedGrav) then
        call write_checkpoint_output(eta, nx, ny, 1, &
          n, 'checkpoints/eta.')
      end if

    end if

    if (diagwrite .eq. 0) then
      ! not saving diagnostics. Move one.
    else if (mod(n-1, diagwrite) .eq. 0) then
      call write_diag_output(h, nx, ny, layers, n, 'output/diagnostic.h.csv')
      call write_diag_output(u, nx, ny, layers, n, 'output/diagnostic.u.csv')
      call write_diag_output(v, nx, ny, layers, n, 'output/diagnostic.v.csv')
      if (.not. RedGrav) then
        call write_diag_output(eta, nx, ny, 1, n, 'output/diagnostic.eta.csv')
      end if
    end if

    return
  end subroutine maybe_dump_output

  ! ---------------------------------------------------------------------------

  subroutine read_input_fileH(name, array, default, nx, ny, layers)
    implicit none

    character(60), intent(in) :: name
    double precision, intent(out) :: array(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: default(layers)
    integer, intent(in) :: nx, ny, layers

    double precision array_small(nx, ny, layers)
    integer k



    if (name.ne.'') then
      open(unit=10, form='unformatted', file=name)
      read(10) array_small
      close(10)
      array(1:nx, 1:ny, :) = array_small
      call wrap_fields_3D(array, nx, ny, layers)
    else
      do k = 1, layers
        array(:, :, k) = default(k)
      end do
    end if

    return
  end subroutine read_input_fileH

  ! ---------------------------------------------------------------------------

  subroutine read_input_fileH_2D(name, array, default, nx, ny)
    implicit none

    character(60), intent(in) :: name
    double precision, intent(out) :: array(0:nx+1, 0:ny+1)
    double precision, intent(in) :: default
    integer, intent(in) :: nx, ny

    double precision array_small(nx, ny)

    if (name.ne.'') then
      open(unit=10, form='unformatted', file=name)
      read(10) array_small
      close(10)
      array(1:nx, 1:ny) = array_small
      call wrap_fields_2D(array, nx, ny)
    else
      array = default
    end if

    return
  end subroutine read_input_fileH_2D

  ! ---------------------------------------------------------------------------

  subroutine read_input_fileU(name, array, default, nx, ny, layers)
    implicit none

    character(60), intent(in) :: name
    double precision, intent(out) :: array(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: default
    integer, intent(in) :: nx, ny, layers

    double precision array_small(nx+1, ny, layers)

    if (name.ne.'') then
      open(unit=10, form='unformatted', file=name)
      read(10) array_small
      close(10)
      array(1:nx+1, 1:ny, :) = array_small
      call wrap_fields_3D(array, nx, ny, layers)
    else
      array = default
    end if

    return
  end subroutine read_input_fileU

  ! ---------------------------------------------------------------------------

  subroutine read_input_fileV(name, array, default, nx, ny, layers)
    implicit none

    character(60), intent(in) :: name
    double precision, intent(out) :: array(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: default
    integer, intent(in) :: nx, ny, layers

    double precision array_small(nx, ny+1, layers)

    if (name.ne.'') then
      open(unit=10, form='unformatted', file=name)
      read(10) array_small
      close(10)
      array(1:nx, 1:ny+1, :) = array_small
      call wrap_fields_3D(array, nx, ny, layers)
    else
      array = default
    end if

    return
  end subroutine read_input_fileV

  ! ---------------------------------------------------------------------------

  subroutine read_input_file_time_series(name, array, default, nTimeSteps)
    implicit none

    character(60), intent(in) :: name
    double precision, intent(out) :: array(nTimeSteps)
    double precision, intent(in) :: default
    integer, intent(in) :: nTimeSteps

    if (name.ne.'') then
      open(unit=10, form='unformatted', file=name)
      read(10) array
      close(10)
    else
      array = default
    end if

    return
  end subroutine read_input_file_time_series

  !-----------------------------------------------------------------
  !> Write snapshot output of 3d field

  subroutine write_output_3d(array, nx, ny, layers, xstep, ystep, &
      n, name)
    implicit none

    double precision, intent(in) :: array(0:nx+1, 0:ny+1, layers)
    integer,          intent(in) :: nx, ny, layers, xstep, ystep
    integer,          intent(in) :: n
    character(*),     intent(in) :: name

    character(10)  :: num

    write(num, '(i10.10)') n

    ! Output the data to a file
    open(unit=10, status='replace', file=name//num, &
        form='unformatted')
    write(10) array(1:nx+xstep, 1:ny+ystep, :)
    close(10)

    return
  end subroutine write_output_3d

  !-----------------------------------------------------------------
  !> Write snapshot output of 3d field

  subroutine write_checkpoint_output(array, nx, ny, layers, &
      n, name)
    implicit none

    double precision, intent(in) :: array(0:nx+1, 0:ny+1, layers)
    integer,          intent(in) :: nx, ny, layers
    integer,          intent(in) :: n
    character(*),     intent(in) :: name

    character(10)  :: num

    write(num, '(i10.10)') n

    ! Output the data to a file
    open(unit=10, status='replace', file=name//num, &
        form='unformatted')
    write(10) array
    close(10)

    return
  end subroutine write_checkpoint_output

  !-----------------------------------------------------------------
  !> Write snapshot output of 2d field

  subroutine write_output_2d(array, nx, ny, xstep, ystep, &
      n, name)
    implicit none

    double precision, intent(in) :: array(0:nx+1, 0:ny+1)
    integer,          intent(in) :: nx, ny, xstep, ystep
    integer,          intent(in) :: n
    character(*),     intent(in) :: name

    character(10)  :: num
    
    write(num, '(i10.10)') n

    ! Output the data to a file
    open(unit=10, status='replace', file=name//num, &
        form='unformatted')
    write(10) array(1:nx+xstep, 1:ny+ystep)
    close(10)

    return
  end subroutine write_output_2d


  !-----------------------------------------------------------------
  !> create a diagnostics file

  subroutine create_diag_file(layers, filename, arrayname, niter0)
    implicit none

    integer,          intent(in) :: layers
    character(*),     intent(in) :: filename
    character(*),     intent(in) :: arrayname
    integer,          intent(in) :: niter0

    integer        :: k
    logical        :: lex
    character(2)   :: layer_number
    character(17)   :: header((4*layers)+1)


    ! prepare header for file
    header(1) = 'timestep'
    do k = 1, layers
      write(layer_number, '(i2.2)') k
      header(2+(4*(k-1))) = 'mean'//layer_number
      header(3+(4*(k-1))) = 'max'//layer_number
      header(4+(4*(k-1))) = 'min'//layer_number
      header(5+(4*(k-1))) = 'std'//layer_number
    end do

    INQUIRE(file=filename, exist=lex)

    if (niter0 .eq. 0) then
      ! starting a nw run, intialise diagnostics files, but warn if 
      ! they were already there
      if (lex) then
        print "(A)", &
          "Starting a new run (niter0=0), but diagnostics file for "//arrayname//" already exists. Overwriting old file."
      else if (.not. lex) then
        print "(A)", &
          "Diagnostics file for "//arrayname//" does not exist. Creating it now."
      end if

      open(unit=10, status='replace', file=filename, &
        form='formatted')
      write (10,'(*(G0.4,:,","))') header
      close(10)

    else if (niter0 .ne. 0) then
      ! restarting from checkpoint, diagnostics file may or may not exist.
      if (lex) then
        print "(A)", &
          "Diagnostics file for "//arrayname//" already exists. Appending to it."
      else if (.not. lex) then
        print "(A)", &
          "Diagnostics file for "//arrayname//" does not exist. Creating it now."

        open(unit=10, status='new', file=filename, &
          form='formatted')
        write (10,'(*(G0.4,:,","))') header
        close(10)
      end if
    end if

    return
  end subroutine create_diag_file

  !-----------------------------------------------------------------
  !> Save diagnostistics of given fields

  subroutine write_diag_output(array, nx, ny, layers, &
      n, filename)
    implicit none

    double precision, intent(in) :: array(0:nx+1, 0:ny+1, layers)
    integer,          intent(in) :: nx, ny, layers
    integer,          intent(in) :: n
    character(*),     intent(in) :: filename

    double precision :: diag_out(4*layers)
    integer          :: k

    ! prepare data for file
    do k = 1, layers
      diag_out(1+(4*(k-1))) = sum(array(:,:,k))/dble(size(array(:,:,k))) ! mean
      diag_out(2+(4*(k-1))) = maxval(array(:,:,k))
      diag_out(3+(4*(k-1))) = minval(array(:,:,k))
      diag_out(4+(4*(k-1))) = sqrt( sum( (array(:,:,k) - diag_out(1+(4*(k-1))))**2)/ &
                            dble(size(array(:,:,k))))
    end do

    ! Output the data to a file
    open(unit=10, status='old', file=filename, &
        form='formatted', position='append')
    write (10,'(i10.10, ",", *(G22.15,:,","))') n, diag_out
    close(10)

    return
  end subroutine write_diag_output

end module io