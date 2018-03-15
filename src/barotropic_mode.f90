module barotropic_mode
  use io
  use vorticity
  use boundaries
  use enforce_thickness
  
  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Calculate the barotropic u velocity

  subroutine calc_baro_u(ub, u, h, eta, freesurfFac, nx, ny, layers)
    implicit none

    double precision, intent(out) :: ub(nx+1, ny)
    double precision, intent(in)  :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: eta(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: freesurfFac
    integer, intent(in) :: nx, ny, layers

    integer i, j, k
    double precision h_temp(0:nx+1, 0:ny+1, layers)

    ub = 0d0

    h_temp = h
    ! add free surface elevation to the upper layer
    h_temp(:, :, 1) = h(:, :, 1) + eta*freesurfFac

    do i = 1, nx+1
      do j = 1, ny
        do k = 1, layers
          ub(i,j) = ub(i,j) + u(i,j,k)*(h_temp(i,j,k)+h_temp(i-1,j,k))/2d0
        end do
      end do
    end do

    return
  end subroutine calc_baro_u

  ! ---------------------------------------------------------------------------
  !> Calculate the barotropic v velocity

  subroutine calc_baro_v(vb, v, h, eta, freesurfFac, nx, ny, layers)
    implicit none

    double precision, intent(out) :: vb(nx, ny+1)
    double precision, intent(in)  :: v(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: eta(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: freesurfFac
    integer, intent(in) :: nx, ny, layers

    integer i, j, k
    double precision h_temp(0:nx+1, 0:ny+1, layers)

    vb = 0d0

    h_temp = h
    ! add free surface elevation to the upper layer
    h_temp(:, :, 1) = h(:, :, 1) + eta*freesurfFac

    do i = 1, nx
      do j = 1, ny+1
        do k = 1, layers
          vb(i,j) = vb(i,j) + v(i,j,k)*(h_temp(i,j,k)+h_temp(i,j-1,k))/2d0
        end do
      end do
    end do

    return
  end subroutine calc_baro_v

  ! ---------------------------------------------------------------------------
  !> Calculate the free surface anomaly using the velocities
  !! timestepped with the tendencies excluding the free surface
  !! pressure gradient.

  subroutine calc_eta_star(ub, vb, eta, etastar, &
      freesurfFac, nx, ny, dx, dy, dt)
    implicit none

    double precision, intent(in)  :: ub(nx+1, ny)
    double precision, intent(in)  :: vb(nx, ny+1)
    double precision, intent(in)  :: eta(0:nx+1, 0:ny+1)
    double precision, intent(out) :: etastar(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: freesurfFac
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: dx, dy, dt

    integer i, j

    etastar = 0d0

    do i = 1, nx
      do j = 1, ny
        etastar(i,j) = freesurfFac*eta(i,j) - &
            dt*((ub(i+1,j) - ub(i,j))/dx + (vb(i,j+1) - vb(i,j))/dy)
      end do
    end do

    call wrap_fields_2D(etastar, nx, ny)

    return
  end subroutine calc_eta_star

  ! ---------------------------------------------------------------------------
  !> Use the successive over-relaxation algorithm to solve the backwards
  !! Euler timestepping for the free surface anomaly, or for the surface
  !! pressure required to keep the barotropic flow nondivergent.

  subroutine SOR_solver(a, etanew, etastar, nx, ny, dt, &
      rjac, eps, maxits, n)
    implicit none

    double precision, intent(in)  :: a(5, nx, ny)
    double precision, intent(out) :: etanew(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: etastar(0:nx+1, 0:ny+1)
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: dt
    double precision, intent(in) :: rjac, eps
    integer, intent(in) :: maxits, n

    integer i, j, nit
    double precision rhs(nx, ny)
    double precision res(nx, ny)
    double precision norm, norm0
    double precision relax_param

    rhs = -etastar(1:nx,1:ny)/dt**2
    ! first guess for etanew
    etanew = etastar

    relax_param = 1.d0 ! successive over-relaxation parameter

    ! Calculate initial residual, so that we can stop the loop when the
    ! current residual = norm0*eps
    norm0 = 0.d0
    do i = 1, nx
      do j = 1, ny
        res(i,j) = &
            a(1,i,j)*etanew(i+1,j) &
            + a(2,i,j)*etanew(i,j+1) &
            + a(3,i,j)*etanew(i-1,j) &
            + a(4,i,j)*etanew(i,j-1) &
            + a(5,i,j)*etanew(i,j)   &
            - rhs(i,j)
        norm0 = norm0 + abs(res(i,j))
        etanew(i,j) = etanew(i,j)-relax_param*res(i,j)/a(5,i,j)
      end do
    end do


    do nit = 1, maxits
      norm = 0.d0
      do i = 1, nx
        do j = 1, ny
          res(i,j) = &
              a(1,i,j)*etanew(i+1,j) &
              + a(2,i,j)*etanew(i,j+1) &
              + a(3,i,j)*etanew(i-1,j) &
              + a(4,i,j)*etanew(i,j-1) &
              + a(5,i,j)*etanew(i,j)   &
              - rhs(i,j)
          norm = norm + abs(res(i,j))
          etanew(i,j) = etanew(i,j)-relax_param*res(i,j)/(a(5,i,j))
        end do
      end do
      if (nit.eq.1) then
        relax_param = 1.d0/(1.d0-0.5d0*rjac**2)
      else
        relax_param = 1.d0/(1.d0-0.25d0*rjac**2*relax_param)
      end if

      call wrap_fields_2D(etanew, nx, ny)

      if (nit.gt.1.and.norm.lt.eps*norm0) then

        return

      end if
    end do

    write(17, "(A, I0)") 'Warning: maximum SOR iterations exceeded at time step ', n

    return
  end subroutine SOR_solver

  ! ---------------------------------------------------------------------------

  subroutine create_Hypre_grid(MPI_COMM_WORLD, hypre_grid, ilower, iupper, &
            num_procs, myid, nx, ny, ierr)
    implicit none

    integer,   intent(in)  :: MPI_COMM_WORLD
    integer*8, intent(out) :: hypre_grid
    integer,   intent(in)  :: ilower(0:num_procs-1,2)
    integer,   intent(in)  :: iupper(0:num_procs-1,2)
    integer,   intent(in)  :: num_procs
    integer,   intent(in)  :: myid
    integer,   intent(in)  :: nx
    integer,   intent(in)  :: ny
    integer,   intent(out)  :: ierr

#ifdef useExtSolver
    call Hypre_StructGridCreate(MPI_COMM_WORLD, 2, hypre_grid, ierr)

    !do i = 0, num_procs-1
    call HYPRE_StructGridSetExtents(hypre_grid, ilower(myid,:),iupper(myid,:), ierr)
    !end do

    call HYPRE_StructGridSetPeriodic(hypre_grid, [nx, ny], ierr)

    call HYPRE_StructGridAssemble(hypre_grid, ierr)
#endif

    return
  end subroutine create_Hypre_grid

  ! ---------------------------------------------------------------------------

  subroutine create_Hypre_A_matrix(MPI_COMM_WORLD, hypre_grid, hypre_A, &
            a, nx, ny, ierr)
    implicit none

    integer,          intent(in)  :: MPI_COMM_WORLD
    integer*8,        intent(in)  :: hypre_grid
    integer*8,        intent(out) :: hypre_A
    double precision, intent(in)  :: a(5, nx, ny)
    integer,          intent(in)  :: nx, ny
    integer,          intent(out) :: ierr

    ! Hypre stencil for creating the A matrix
    integer*8 :: stencil

    integer :: offsets(2,5)
    integer :: indicies(2)
    integer :: i, j

#ifdef useExtSolver

    ! Define the geometry of the stencil.  Each represents a relative
    ! offset (in the index space).
    offsets(1,1) =  0
    offsets(2,1) =  0
    offsets(1,2) = -1
    offsets(2,2) =  0
    offsets(1,3) =  1
    offsets(2,3) =  0
    offsets(1,4) =  0
    offsets(2,4) = -1
    offsets(1,5) =  0
    offsets(2,5) =  1


    call HYPRE_StructStencilCreate(2, 5, stencil, ierr)
    ! this gives a 2D, 5 point stencil centred around the grid point of interest.
    do i = 0, 4
      call HYPRE_StructStencilSetElement(stencil, i, offsets(:,i+1),ierr)
    end do

    call HYPRE_StructMatrixCreate(MPI_COMM_WORLD, hypre_grid, stencil, hypre_A, ierr)

    call HYPRE_StructMatrixInitialize(hypre_A, ierr)

    do i = 1, nx
      do j = 1, ny
        indicies(1) = i
        indicies(2) = j

        call HYPRE_StructMatrixSetValues(hypre_A, &
            indicies, 1, 0, &
            a(5,i,j), ierr)
        call HYPRE_StructMatrixSetValues(hypre_A, &
            indicies, 1, 1, &
            a(3,i,j), ierr)
        call HYPRE_StructMatrixSetValues(hypre_A, &
            indicies, 1, 2, &
            a(1,i,j), ierr)
        call HYPRE_StructMatrixSetValues(hypre_A, &
            indicies, 1, 3, &
            a(4,i,j), ierr)
        call HYPRE_StructMatrixSetValues(hypre_A, &
            indicies, 1, 4, &
            a(2,i,j), ierr)
      end do
    end do

    call HYPRE_StructMatrixAssemble(hypre_A, ierr)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

#endif

    return
  end subroutine create_Hypre_A_matrix

  ! ---------------------------------------------------------------------------

  subroutine Ext_solver(MPI_COMM_WORLD, hypre_A, hypre_grid, myid, num_procs, &
      ilower, iupper, etastar, &
      etanew, nx, ny, dt, maxits, eps, ierr)
    implicit none

    integer,          intent(in)  :: MPI_COMM_WORLD
    integer*8,        intent(in)  :: hypre_A
    integer*8,        intent(in)  :: hypre_grid
    integer,          intent(in)  :: myid
    integer,          intent(in)  :: num_procs
    integer,          intent(in)  :: ilower(0:num_procs-1,2)
    integer,          intent(in)  :: iupper(0:num_procs-1,2)
    double precision, intent(in)  :: etastar(0:nx+1, 0:ny+1)
    double precision, intent(out) :: etanew(0:nx+1, 0:ny+1)
    integer,          intent(in)  :: nx, ny
    double precision, intent(in)  :: dt
    integer,          intent(in)  :: maxits
    double precision, intent(in)  :: eps
    integer,          intent(out) :: ierr

    integer          :: i, j ! loop variables
    integer*8        :: hypre_b
    integer*8        :: hypre_x
    integer*8        :: hypre_solver
    integer*8        :: precond
    double precision, dimension(:),     allocatable :: values

    integer :: nx_tile, ny_tile

    nx_tile = iupper(myid,1)-ilower(myid,1) + 1
    ny_tile = iupper(myid,2)-ilower(myid,2) + 1

    allocate(values(nx_tile*ny_tile))
    ! just nx*ny for the tile this processor owns


    ! A currently unused variable that can be used to
    ! print information from the solver - see comments below.
  !  double precision :: hypre_out(2)


    ! wrap this code in preprocessing flags to allow the model to be compiled without the external library, if desired.
#ifdef useExtSolver
    ! Create the rhs vector, b
    call HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid, hypre_b, ierr)
    call HYPRE_StructVectorInitialize(hypre_b, ierr)

    ! set rhs values (vector b)
    do j = ilower(myid,2), iupper(myid,2) ! loop over every grid point
      do i = ilower(myid,1), iupper(myid,1)
    ! the 2D array is being laid out like
    ! [x1y1, x2y1, x3y1, x1y2, x2y2, x3y2, x1y3, x2y3, x3y3]
      values( ((j-1)*nx_tile + i) ) = -etastar(i,j)/dt**2
      end do
    end do

    call HYPRE_StructVectorSetBoxValues(hypre_b, &
      ilower(myid,:), iupper(myid,:), values, ierr)

    call HYPRE_StructVectorAssemble(hypre_b, ierr)

    ! now create the x vector
    call HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid, hypre_x, ierr)
    call HYPRE_StructVectorInitialize(hypre_x, ierr)

    call HYPRE_StructVectorSetBoxValues(hypre_x, &
      ilower(myid,:), iupper(myid,:), values, ierr)

    call HYPRE_StructVectorAssemble(hypre_x, ierr)

    ! now create the solver and solve the equation.
    ! Choose the solver
    call HYPRE_StructPCGCreate(MPI_COMM_WORLD, hypre_solver, ierr)

    ! Set some parameters
    call HYPRE_StructPCGSetMaxIter(hypre_solver, maxits, ierr)
    call HYPRE_StructPCGSetTol(hypre_solver, eps, ierr)
    ! other options not explained by user manual but present in examples
    ! call HYPRE_StructPCGSetMaxIter(hypre_solver, 50 );
    ! call HYPRE_StructPCGSetTol(hypre_solver, 1.0e-06 );
    call HYPRE_StructPCGSetTwoNorm(hypre_solver, 1 );
    call HYPRE_StructPCGSetRelChange(hypre_solver, 0 );
    call HYPRE_StructPCGSetPrintLevel(hypre_solver, 1 ); ! 2 will print each CG iteration
    call HYPRE_StructPCGSetLogging(hypre_solver, 1);

    ! use an algebraic multigrid preconditioner
    call HYPRE_BoomerAMGCreate(precond, ierr)
    ! values taken from hypre library example number 5
    ! print less solver info since a preconditioner
    call HYPRE_BoomerAMGSetPrintLevel(precond, 1, ierr);
    ! Falgout coarsening
    call HYPRE_BoomerAMGSetCoarsenType(precond, 6, ierr)
    ! old defaults
    call HYPRE_BoomerAMGSetOldDefault(precond, ierr)
    ! SYMMETRIC G-S/Jacobi hybrid relaxation
    call HYPRE_BoomerAMGSetRelaxType(precond, 6, ierr)
    ! Sweeeps on each level
    call HYPRE_BoomerAMGSetNumSweeps(precond, 1, ierr)
    ! conv. tolerance
    call HYPRE_BoomerAMGSetTol(precond, 0.0d0, ierr)
    ! do only one iteration!
    call HYPRE_BoomerAMGSetMaxIter(precond, 1, ierr)

    ! set amg as the pcg preconditioner
    call HYPRE_StructPCGSetPrecond(hypre_solver, 2, precond, ierr)


    ! now we set the system up and do the actual solve!
    call HYPRE_StructPCGSetup(hypre_solver, hypre_A, hypre_b, &
                              hypre_x, ierr)

    call HYPRE_ParCSRPCGSolve(hypre_solver, hypre_A, hypre_b, &
                              hypre_x, ierr)

    ! code for printing out results from the external solver
    ! Not being used, but left here since the manual isn't very helpful
    ! and this may be useful in the future.
    ! call HYPRE_ParCSRPCGGetNumIterations(hypre_solver, &
    !   hypre_out(1), ierr)
    ! print *, 'num iterations = ', hypre_out(1)

    ! call HYPRE_ParCSRPCGGetFinalRelative(hypre_solver, &
    !   hypre_out(2), ierr)
    ! print *, 'final residual norm = ', hypre_out(2)

    call HYPRE_StructVectorGetBoxValues(hypre_x, &
      ilower(myid,:), iupper(myid,:), values, ierr)

    do j = ilower(myid,2), iupper(myid,2) ! loop over every grid point
      do i = ilower(myid,1), iupper(myid,1)
      etanew(i,j) = values( ((j-1)*nx_tile + i) )
      end do
    end do

    ! debugging commands from hypre library - dump out a single
    ! copy of these two variables. Can be used to check that the
    ! values have been properly allocated.
    ! call HYPRE_StructVectorPrint(hypre_x, ierr)
    ! call HYPRE_StructMatrixPrint(hypre_A, ierr)

    call HYPRE_StructPCGDestroy(hypre_solver, ierr)
    call HYPRE_BoomerAMGDestroy(precond, ierr)
    call HYPRE_StructVectorDestroy(hypre_x, ierr)
    call HYPRE_StructVectorDestroy(hypre_b, ierr)

#endif

    return
  end subroutine Ext_solver

  ! ---------------------------------------------------------------------------

  !> Update velocities using the barotropic tendency due to the pressure
  !> gradient.

  subroutine update_velocities_for_barotropic_tendency(array, etanew, g_vec, &
      xstep, ystep, dspace, dt, nx, ny, layers)
    implicit none

    double precision, intent(inout) :: array(0:nx+1, 0:ny+1, layers)
    double precision, intent(in) :: etanew(0:nx+1, 0:ny+1)
    double precision, intent(in) :: g_vec(layers)
    integer, intent(in) :: xstep, ystep
    double precision, intent(in) :: dspace, dt
    integer, intent(in) :: nx, ny, layers

    integer i, j, k
    double precision baro_contrib

    ! TODO Assert that xstep and ystep are either 1, 0 or 0, 1.

    do i = xstep, nx
      do j = ystep, ny
        do k = 1, layers
          baro_contrib = &
              -g_vec(1)*(etanew(i,j) - etanew(i-xstep,j-ystep))/(dspace)
          array(i,j,k) = array(i,j,k) + dt*baro_contrib
        end do
      end do
    end do


    return
  end subroutine update_velocities_for_barotropic_tendency


  ! ---------------------------------------------------------------------------
  !> Compute derivatives of the depth field for the pressure solver

  subroutine calc_A_matrix(a, depth, g, dx, dy, nx, ny, freesurfFac, dt, &
            hfacW, hfacE, hfacS, hfacN)
    implicit none

    double precision, intent(out) :: a(5, nx, ny)
    double precision, intent(in)  :: depth(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: g, dx, dy
    integer, intent(in)           :: nx, ny
    double precision, intent(in)  :: freesurfFac
    double precision, intent(in)  :: dt
    double precision, intent(in)  :: hfacW(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: hfacE(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: hfacN(0:nx+1, 0:ny+1)
    double precision, intent(in)  :: hfacS(0:nx+1, 0:ny+1)

    integer i, j

    do j = 1, ny
      do i = 1, nx
        a(1,i,j) = g*0.5*(depth(i+1,j)+depth(i,j))*hfacE(i,j)/dx**2
        a(2,i,j) = g*0.5*(depth(i,j+1)+depth(i,j))*hfacN(i,j)/dy**2
        a(3,i,j) = g*0.5*(depth(i,j)+depth(i-1,j))*hfacW(i,j)/dx**2
        a(4,i,j) = g*0.5*(depth(i,j)+depth(i,j-1))*hfacS(i,j)/dy**2
      end do
    end do

    do j = 1, ny
      do i = 1, nx
        a(5,i,j) = -a(1,i,j)-a(2,i,j)-a(3,i,j)-a(4,i,j) - freesurfFac/dt**2
      end do
    end do

    return
  end subroutine calc_A_matrix

  ! ---------------------------------------------------------------------------
  !> Do the isopycnal layer physics

  subroutine barotropic_correction(hnew, unew, vnew, eta, etanew, depth, a, &
      dx, dy, wetmask, hfacW, hfacS, dt, &
      maxits, eps, rjac, freesurfFac, thickness_error, &
      debug_level, g_vec, nx, ny, layers, n, &
       MPI_COMM_WORLD, myid, num_procs, ilower, iupper, &
       hypre_grid, hypre_A, ierr)

    implicit none

    double precision, intent(inout) :: hnew(0:nx+1, 0:ny+1, layers)
    double precision, intent(inout) :: unew(0:nx+1, 0:ny+1, layers)
    double precision, intent(inout) :: vnew(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)    :: eta(0:nx+1, 0:ny+1)
    double precision, intent(out)   :: etanew(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: depth(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: a(5, nx, ny)
    double precision, intent(in)    :: dx, dy
    double precision, intent(in)    :: wetmask(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: hfacW(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: hfacS(0:nx+1, 0:ny+1)
    double precision, intent(in)    :: dt
    integer,          intent(in)    :: maxits
    double precision, intent(in)    :: eps, rjac, freesurfFac, thickness_error
    integer,          intent(in)    :: debug_level
    double precision, intent(in)    :: g_vec(layers)
    integer,          intent(in)    :: nx, ny, layers, n
    integer,          intent(in)    :: MPI_COMM_WORLD
    integer,          intent(in)    :: myid, num_procs
    integer,          intent(in)    :: ilower(0:num_procs-1,2)
    integer,          intent(in)    :: iupper(0:num_procs-1,2)
    integer*8,        intent(in)    :: hypre_grid
    integer*8,        intent(in)    :: hypre_A
    integer,          intent(out) :: ierr

    ! barotropic velocity components (for pressure solver)
    double precision :: ub(nx+1, ny)
    double precision :: vb(nx, ny+1)
    double precision :: etastar(0:nx+1, 0:ny+1)

    character(10)    :: num

    ! Calculate the barotropic velocities
    call calc_baro_u(ub, unew, hnew, eta, freesurfFac, nx, ny, layers)
    call calc_baro_v(vb, vnew, hnew, eta, freesurfFac, nx, ny, layers)
    
    if (debug_level .ge. 4) then

      write(num, '(i10.10)') n

      ! Output the data to a file
      open(unit=10, status='replace', file='output/'//'snap.ub.'//num, &
          form='unformatted')
      write(10) ub(1:nx+1, 1:ny)
      close(10)

      open(unit=10, status='replace', file='output/'//'snap.vb.'//num, &
          form='unformatted')
      write(10) vb(1:nx, 1:ny+1)
      close(10)
    end if


    ! Calculate divergence of ub and vb, and solve for the pressure
    ! field that removes it
    call calc_eta_star(ub, vb, eta, etastar, freesurfFac, nx, ny, dx, dy, dt)
    ! print *, maxval(abs(etastar))
    if (debug_level .ge. 4) then
      call write_output_2d(etastar, nx, ny, 0, 0, &
        n, 'output/snap.eta_star.')
    end if

    ! Prevent barotropic signals from bouncing around outside the
    ! wet region of the model.
    ! etastar = etastar*wetmask
#ifndef useExtSolver
    call SOR_solver(a, etanew, etastar, nx, ny, &
       dt, rjac, eps, maxits, n)
    ! print *, maxval(abs(etanew))
#endif

#ifdef useExtSolver
    call Ext_solver(MPI_COMM_WORLD, hypre_A, hypre_grid, myid, num_procs, &
      ilower, iupper, etastar, &
      etanew, nx, ny, dt, maxits, eps, ierr)
#endif

    if (debug_level .ge. 4) then
      call write_output_2d(etanew, nx, ny, 0, 0, &
        n, 'output/snap.eta_new.')
    end if

    etanew = etanew*wetmask

    call wrap_fields_2D(etanew, nx, ny)

    ! Now update the velocities using the barotropic tendency due to
    ! the pressure gradient.
    call update_velocities_for_barotropic_tendency(unew, etanew, g_vec, &
        1, 0, dx, dt, nx, ny, layers)
    call update_velocities_for_barotropic_tendency(vnew, etanew, g_vec, &
        0, 1, dy, dt, nx, ny, layers)

    ! We now have correct velocities at the next time step, but the
    ! layer thicknesses were updated using the velocities without
    ! the barotropic pressure contribution. Force consistency
    ! between layer thicknesses and ocean depth by scaling
    ! thicknesses to agree with free surface.
    call enforce_depth_thickness_consistency(hnew, etanew, depth, &
        freesurfFac, thickness_error, nx, ny, layers)

    ! Apply the boundary conditions
    call apply_boundary_conditions(unew, hfacW, wetmask, nx, ny, layers)
    call apply_boundary_conditions(vnew, hfacS, wetmask, nx, ny, layers)

    return
  end subroutine barotropic_correction

end module barotropic_mode
