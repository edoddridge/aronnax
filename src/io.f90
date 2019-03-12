module io
  use end_run
  use boundaries
  use exchange
  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Write output if it's time

  subroutine maybe_dump_output(h, hav, u, uav, v, vav, eta, etaav, &
          dudt, dvdt, dhdt, AB_order, &
          wind_x, wind_y, nx, ny, layers, ilower, iupper, &
          xlow, xhigh, ylow, yhigh, OL, num_procs, myid, &
          n, nwrite, avwrite, checkpointwrite, diagwrite, &
          RedGrav, DumpWind, debug_level)
    implicit none

    double precision, intent(in)    :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(inout) :: hav(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(inout) :: uav(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(inout) :: vav(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: eta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(inout) :: etaav(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: dudt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision, intent(in)    :: dvdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision, intent(in)    :: dhdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    integer,          intent(in)    :: AB_order
    double precision, intent(in)    :: wind_x(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: wind_y(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    integer,          intent(in)    :: nx, ny, layers
    integer,          intent(in)    :: ilower(0:num_procs-1,2)
    integer,          intent(in)    :: iupper(0:num_procs-1,2)
    integer,          intent(in)    :: xlow, xhigh, ylow, yhigh, OL
    integer,          intent(in)    :: num_procs, myid
    integer,          intent(in)    :: n
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
      
      call write_output_3d(h, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 0, &
                          n, 'output/snap.h.', num_procs, myid)
      call write_output_3d(u, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 1, 0, &
                          n, 'output/snap.u.', num_procs, myid)
      call write_output_3d(v, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 1, &
                          n, 'output/snap.v.', num_procs, myid)


      if (.not. RedGrav) then
        call write_output_3d(eta, nx, ny, 1, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 0, &
                          n, 'output/snap.eta.', num_procs, myid)
      end if

      if (DumpWind .eqv. .true.) then
        call write_output_3d(wind_x, nx, ny, 1, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 1, 0, &
                          n, 'output/wind_x.', num_procs, myid)
        call write_output_3d(wind_y, nx, ny, 1, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 1, &
                          n, 'output/wind_y.', num_procs, myid)
      end if

      if (debug_level .ge. 1) then
        call write_output_3d(dhdt(:,:,:,1), nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 0, &
          n, 'output/debug.dhdt.', num_procs, myid)
        call write_output_3d(dudt(:,:,:,1), nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 1, 0, &
          n, 'output/debug.dudt.', num_procs, myid)
        call write_output_3d(dvdt(:,:,:,1), nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 1, &
          n, 'output/debug.dvdt.', num_procs, myid)
      end if

      ! Check if there are NaNs in the data
      call break_if_NaN(h, xlow, xhigh, ylow, yhigh, layers, OL, n)
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

        call write_output_3d(hav, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 0, &
        n, 'output/av.h.', num_procs, myid)
        call write_output_3d(uav, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 1, 0, &
        n, 'output/av.u.', num_procs, myid)
        call write_output_3d(vav, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 1, &
        n, 'output/av.v.', num_procs, myid)


        if (.not. RedGrav) then
          call write_output_3d(etaav, nx, ny, 1, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 0, &
                          n, 'output/av.eta.', num_procs, myid)
        end if

        ! Check if there are NaNs in the data
        call break_if_NaN(h, xlow, xhigh, ylow, yhigh, layers, OL, n)
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
      call write_checkpoint_output(h, nx, ny, layers, ilower, iupper, &
                        xlow, xhigh, ylow, yhigh, OL, 1, &
                        n, 'checkpoints/h.', num_procs, myid)
      call write_checkpoint_output(u, nx, ny, layers, ilower, iupper, &
                        xlow, xhigh, ylow, yhigh, OL, 1, &
                        n, 'checkpoints/u.', num_procs, myid)
      call write_checkpoint_output(v, nx, ny, layers, ilower, iupper, &
                        xlow, xhigh, ylow, yhigh, OL, 1, &
                        n, 'checkpoints/v.', num_procs, myid)

      call write_checkpoint_output(dhdt, nx, ny, layers, ilower, iupper, &
                        xlow, xhigh, ylow, yhigh, OL, AB_order, &
                        n, 'checkpoints/dhdt.', num_procs, myid)
      call write_checkpoint_output(dudt, nx, ny, layers, ilower, iupper, &
                        xlow, xhigh, ylow, yhigh, OL, AB_order, &
                        n, 'checkpoints/dudt.', num_procs, myid)
      call write_checkpoint_output(dvdt, nx, ny, layers, ilower, iupper, &
                        xlow, xhigh, ylow, yhigh, OL, AB_order, &
                        n, 'checkpoints/dvdt.', num_procs, myid)

      if (.not. RedGrav) then
        call write_checkpoint_output(eta, nx, ny, 1, ilower, iupper, &
                        xlow, xhigh, ylow, yhigh, OL, 1, &
                        n, 'checkpoints/eta.', num_procs, myid)
      end if

    end if

    if (diagwrite .eq. 0) then
      ! not saving diagnostics. Move one.
    else if (mod(n-1, diagwrite) .eq. 0) then

      call write_diag_output(h, nx, ny, layers, ilower, iupper, &
                                xlow, xhigh, ylow, yhigh, OL, &
                                n, 'output/diagnostic.h.csv', num_procs, myid)
      call write_diag_output(u, nx, ny, layers, ilower, iupper, &
                                xlow, xhigh, ylow, yhigh, OL, &
                                n, 'output/diagnostic.u.csv', num_procs, myid)
      call write_diag_output(v, nx, ny, layers, ilower, iupper, &
                                xlow, xhigh, ylow, yhigh, OL, &
                                n, 'output/diagnostic.v.csv', num_procs, myid)
      if (.not. RedGrav) then
        call write_diag_output(eta, nx, ny, 1, ilower, iupper, &
                                xlow, xhigh, ylow, yhigh, OL, &
                                n, 'output/diagnostic.eta.csv', num_procs, myid)
      end if
    end if

    return
  end subroutine maybe_dump_output

  ! ---------------------------------------------------------------------------

  subroutine read_input_fileH(name, array, default, nx, ny, layers, OL, &
                                  xlow, xhigh, ylow, yhigh, myid)
    implicit none

    character(60), intent(in) :: name
     double precision, intent(out) :: array(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: default(layers)
    integer, intent(in) :: nx, ny, layers, OL
    integer, intent(in) :: xlow, xhigh, ylow, yhigh
    integer, intent(in) :: myid

    double precision global_array(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision array_small(nx, ny, layers)
    integer i, j, k

    if (name.ne.'') then
      open(unit=10, form='unformatted', file=name)
      read(10) array_small
      close(10)
      global_array(1:nx, 1:ny, :) = array_small
      call wrap_fields_3D(global_array, nx, ny, layers, OL)
    else
      do k = 1, layers
        global_array(:, :, k) = default(k)
      end do
    end if

    do i = xlow-OL, xhigh+OL
      do j = ylow-OL, yhigh+OL
        array(i,j,:) = global_array(i,j,:)
      end do
    end do

    return
  end subroutine read_input_fileH

  ! ---------------------------------------------------------------------------

  subroutine read_input_fileH_2D(name, array, default, nx, ny, OL, &
                                  xlow, xhigh, ylow, yhigh, myid)
    implicit none

    character(60), intent(in) :: name
    double precision, intent(out) :: array(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: default
    integer, intent(in) :: nx, ny, OL
    integer, intent(in) :: xlow, xhigh, ylow, yhigh
    integer, intent(in) :: myid

    double precision global_array(1-OL:nx+OL, 1-OL:ny+OL)
    double precision array_small(nx, ny)
    integer i, j


    if (name.ne.'') then
      open(unit=10, form='unformatted', file=name)
      read(10) array_small
      close(10)
      global_array(1:nx, 1:ny) = array_small
      call wrap_fields_2D(global_array, nx, ny, OL)
    else
      global_array = default
    end if

    do i = xlow-OL, xhigh+OL
      do j = ylow-OL, yhigh+OL
        array(i,j) = global_array(i,j)
      end do
    end do

    return
  end subroutine read_input_fileH_2D

  ! ---------------------------------------------------------------------------

  subroutine read_input_fileU(name, array, default, nx, ny, layers, OL, &
                              xlow, xhigh, ylow, yhigh, myid)
    implicit none

    character(60), intent(in) :: name
    double precision, intent(out) :: array(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: default
    integer, intent(in) :: nx, ny, layers, OL
    integer, intent(in) :: xlow, xhigh, ylow, yhigh
    integer, intent(in) :: myid

    double precision global_array(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision array_small(nx+1, ny, layers)
    integer i, j, k

    if (name.ne.'') then
      open(unit=10, form='unformatted', file=name)
      read(10) array_small
      close(10)
      global_array(1:nx+1, 1:ny, :) = array_small
      call wrap_fields_3D(global_array, nx, ny, layers, OL)
    else
      global_array = default
    end if

    do i = xlow-OL, xhigh+OL
      do j = ylow-OL, yhigh+OL
        array(i,j,:) = global_array(i,j,:)
      end do
    end do

    return
  end subroutine read_input_fileU

  ! ---------------------------------------------------------------------------

  subroutine read_input_fileV(name, array, default, nx, ny, layers, OL, &
                              xlow, xhigh, ylow, yhigh, myid)
    implicit none

    character(60), intent(in) :: name
    double precision, intent(out) :: array(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: default
    integer, intent(in) :: nx, ny, layers, OL
    integer, intent(in) :: xlow, xhigh, ylow, yhigh
    integer, intent(in) :: myid

    double precision global_array(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision array_small(nx, ny+1, layers)
    integer i, j, k

    if (name.ne.'') then
      open(unit=10, form='unformatted', file=name)
      read(10) array_small
      close(10)
      global_array(1:nx, 1:ny+1, :) = array_small
      call wrap_fields_3D(global_array, nx, ny, layers, OL)
    else
      global_array = default
    end if

    do i = xlow-OL, xhigh+OL
      do j = ylow-OL, yhigh+OL
        array(i,j,:) = global_array(i,j,:)
      end do
    end do

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

  subroutine write_output_3d(array, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, xstep, ystep, &
                          n, name, num_procs, myid)
    implicit none

    double precision, intent(in) :: array(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    integer,          intent(in) :: nx, ny, layers
    integer,          intent(in) :: ilower(0:num_procs-1,2)
    integer,          intent(in) :: iupper(0:num_procs-1,2)
    integer,          intent(in) :: xlow, xhigh, ylow, yhigh, OL, xstep, ystep
    integer,          intent(in) :: n
    character(*),     intent(in) :: name
    integer,          intent(in) :: num_procs, myid

    character(10)  :: num
    double precision :: global_array(1-OL:nx+OL, 1-OL:ny+OL, layers)

    write(num, '(i10.10)') n

    call collect_global_array(array, global_array, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)

    if (myid .eq. 0) then
      ! Output the data to a file
      open(unit=10, status='replace', file=name//num, &
          form='unformatted')
      write(10) global_array(1:nx+xstep, 1:ny+ystep, :)
      close(10)
    end if

    return
  end subroutine write_output_3d

  !-----------------------------------------------------------------
  !> Write checkpoint file

  subroutine write_checkpoint_output(array, nx, ny, layers, ilower, iupper, &
                        xlow, xhigh, ylow, yhigh, OL, AB_order, &
                        n, name, num_procs, myid)
    implicit none

    double precision, intent(in) :: array(1-OL:nx+OL, 1-OL:ny+OL, layers, AB_order)
    integer,          intent(in) :: nx, ny, layers
    integer,          intent(in) :: ilower(0:num_procs-1,2)
    integer,          intent(in) :: iupper(0:num_procs-1,2)
    integer,          intent(in) :: xlow, xhigh, ylow, yhigh, OL, AB_order
    integer,          intent(in) :: n
    character(*),     intent(in) :: name
    integer,          intent(in) :: num_procs, myid

    character(10)    :: num
    double precision :: global_array(1-OL:nx+OL, 1-OL:ny+OL, layers, AB_order)
    integer          :: k

    write(num, '(i10.10)') n

    do k = 1, AB_order
      call collect_global_array(array(:,:,:,k), global_array(:,:,:,k), &
                          nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)
    end do

    if (myid .eq. 0) then
      ! Output the data to a file
      open(unit=10, status='replace', file=name//num, &
          form='unformatted')
      write(10) global_array
      close(10)
    end if

    return
  end subroutine write_checkpoint_output

  ! ---------------------------------------------------------------------------
  !> Load in checkpoint files when restarting a simulation

  subroutine load_checkpoint_files(dhdt, dudt, dvdt, h, u, v, eta, &
            RedGrav, niter0, nx, ny, layers, ilower, iupper, &
            xlow, xhigh, ylow, yhigh, &
            OL, num_procs, myid, AB_order)
  implicit none

    double precision, intent(out) :: dhdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision, intent(out) :: dudt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision, intent(out) :: dvdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers, AB_order)
    double precision, intent(out) :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out) :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out) :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(out) :: eta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    logical,          intent(in)  :: RedGrav
    integer,          intent(in)  :: niter0
    integer,          intent(in)  :: nx, ny, layers
    integer,          intent(in)  :: ilower(0:num_procs-1,2)
    integer,          intent(in)  :: iupper(0:num_procs-1,2)
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, OL
    integer,          intent(in)  :: num_procs, myid, AB_order

    double precision :: h_global(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision :: u_global(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision :: v_global(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision :: dhdt_global(1-OL:nx+OL, 1-OL:ny+OL, layers, AB_order)
    double precision :: dudt_global(1-OL:nx+OL, 1-OL:ny+OL, layers, AB_order)
    double precision :: dvdt_global(1-OL:nx+OL, 1-OL:ny+OL, layers, AB_order)
    double precision :: eta_global(1-OL:nx+OL, 1-OL:ny+OL)

    integer :: i, j, k

    ! dummy variable for loading checkpoints
    character(10)    :: num

    ! load in the state and derivative arrays
    write(num, '(i10.10)') niter0

    open(unit=10, form='unformatted', file='checkpoints/h.'//num)
    read(10) h_global
    close(10)
    open(unit=10, form='unformatted', file='checkpoints/u.'//num)
    read(10) u_global
    close(10)
    open(unit=10, form='unformatted', file='checkpoints/v.'//num)
    read(10) v_global
    close(10)

    open(unit=10, form='unformatted', file='checkpoints/dhdt.'//num)
    read(10) dhdt_global
    close(10)
    open(unit=10, form='unformatted', file='checkpoints/dudt.'//num)
    read(10) dudt_global
    close(10)
    open(unit=10, form='unformatted', file='checkpoints/dvdt.'//num)
    read(10) dvdt_global
    close(10)

    if (.not. RedGrav) then
      open(unit=10, form='unformatted', file='checkpoints/eta.'//num)
      read(10) eta_global
      close(10)
    end if

    do k = 1, layers
      do j = ilower(myid,2) - OL, iupper(myid,2) + OL
        do i = ilower(myid,1) - OL, iupper(myid,1) + OL
          h(i,j,k) = h_global(i,j,k)
          u(i,j,k) = u_global(i,j,k)
          v(i,j,k) = v_global(i,j,k)
          dhdt(i,j,k,:) = dhdt_global(i,j,k,:)
          dudt(i,j,k,:) = dudt_global(i,j,k,:)
          dvdt(i,j,k,:) = dvdt_global(i,j,k,:)
          if (.not. RedGrav) then
            eta(i,j) = eta_global(i,j)
          end if
        end do
      end do
    end do

  end subroutine load_checkpoint_files


  !-----------------------------------------------------------------
  !> Write snapshot output of 2d field

  subroutine write_output_2d(array, nx, ny, OL, xstep, ystep, &
      n, name)
    implicit none

    double precision, intent(in) :: array(1-OL:nx+OL, 1-OL:ny+OL)
    integer,          intent(in) :: nx, ny, OL, xstep, ystep
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

  subroutine write_diag_output(array, nx, ny, layers, ilower, iupper, &
                                xlow, xhigh, ylow, yhigh, OL, &
                                n, filename, num_procs, myid)
    implicit none

    double precision, intent(in) :: array(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    integer,          intent(in) :: nx, ny, layers
    integer,          intent(in) :: ilower(0:num_procs-1,2)
    integer,          intent(in) :: iupper(0:num_procs-1,2)
    integer,          intent(in) :: xlow, xhigh, ylow, yhigh, OL
    integer,          intent(in) :: n
    character(*),     intent(in) :: filename
    integer,          intent(in) :: num_procs, myid

    double precision :: diag_out(4*layers)
    integer          :: k

    double precision :: global_array(1-OL:nx+OL, 1-OL:ny+OL, layers)

    call collect_global_array(array, global_array, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)

    if (myid .eq. 0) then
      ! prepare data for file
      do k = 1, layers
        diag_out(1+(4*(k-1))) = sum(global_array(:,:,k))/dble(size(global_array(:,:,k))) ! mean
        diag_out(2+(4*(k-1))) = maxval(global_array(:,:,k))
        diag_out(3+(4*(k-1))) = minval(global_array(:,:,k))
        diag_out(4+(4*(k-1))) = sqrt( sum( (global_array(:,:,k) - diag_out(1+(4*(k-1))))**2)/ &
                              dble(size(global_array(:,:,k))))
      end do

      ! Output the data to a file
      open(unit=10, status='old', file=filename, &
          form='formatted', position='append')
      write (10,'(i10.10, ",", *(G22.15,:,","))') n, diag_out
      close(10)

    end if

    return
  end subroutine write_diag_output

end module io
