module barotropic_mode
  use io
  use vorticity
  use boundaries
  use enforce_thickness

  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Calculate the barotropic u velocity

  subroutine calc_baro_u(ub, u, h, eta, freesurf_fac, xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    double precision, intent(out) :: ub(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: eta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: freesurf_fac
    integer, intent(in) :: xlow, xhigh, ylow, yhigh, layers, OL

    integer i, j, k
    double precision h_temp(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)

    ub = 0d0

    h_temp = h
    ! add free surface elevation to the upper layer
    h_temp(:, :, 1) = h(:, :, 1) + eta*freesurf_fac

    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL
          ub(i,j) = ub(i,j) + u(i,j,k)*(h_temp(i,j,k)+h_temp(i-1,j,k))/2d0
        end do
      end do
    end do

    return
  end subroutine calc_baro_u

  ! ---------------------------------------------------------------------------
  !> Calculate the barotropic v velocity

  subroutine calc_baro_v(vb, v, h, eta, freesurf_fac, xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    double precision, intent(out) :: vb(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: eta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: freesurf_fac
    integer, intent(in) :: xlow, xhigh, ylow, yhigh, layers, OL

    integer i, j, k
    double precision h_temp(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)

    vb = 0d0

    h_temp = h
    ! add free surface elevation to the upper layer
    h_temp(:, :, 1) = h(:, :, 1) + eta*freesurf_fac

    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL
        do i = xlow-OL+1, xhigh+OL-1
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
      freesurf_fac, xlow, xhigh, ylow, yhigh, OL, dx, dy, dt)
    implicit none

    double precision, intent(in)  :: ub(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: vb(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: eta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(out) :: etastar(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)  :: freesurf_fac
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh, OL
    double precision, intent(in)  :: dx, dy, dt

    integer i, j

    etastar = 0d0

    do j = ylow-OL+1, yhigh+OL-1
      do i = xlow-OL+1, xhigh+OL-1
        etastar(i,j) = freesurf_fac*eta(i,j) - &
            dt*((ub(i+1,j) - ub(i,j))/dx + (vb(i,j+1) - vb(i,j))/dy)
      end do
    end do

    ! call wrap_fields_2D(etastar, nx, ny, OL)

    return
  end subroutine calc_eta_star

  ! ---------------------------------------------------------------------------
  !> Use the successive over-relaxation algorithm to solve the backwards
  !! Euler timestepping for the free surface anomaly, or for the surface
  !! pressure required to keep the barotropic flow nondivergent.

  subroutine SOR_solver(a, etanew, etastar, nx, ny, OL, dt, &
      rjac, eps, maxits, n)
    implicit none

    double precision, intent(in)  :: a(5, nx, ny)
    double precision, intent(out) :: etanew(1-OL:nx+OL, 1-OL:ny+OL)
    double precision, intent(inout)  :: etastar(1-OL:nx+OL, 1-OL:ny+OL)
    integer, intent(in) :: nx, ny, OL
    double precision, intent(in) :: dt
    double precision, intent(in) :: rjac, eps
    integer, intent(in) :: maxits, n

    integer i, j, nit
    double precision rhs(nx, ny)
    double precision res(nx, ny)
    double precision norm, norm0
    double precision relax_param

    call wrap_fields_2D(etastar, nx, ny, OL)

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

      call wrap_fields_2D(etanew, nx, ny, OL)

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

    call HYPRE_StructGridSetExtents(hypre_grid, ilower(myid,:),iupper(myid,:), ierr)

    call HYPRE_StructGridSetPeriodic(hypre_grid, [nx, ny], ierr)

    call HYPRE_StructGridAssemble(hypre_grid, ierr)
#endif

    return
  end subroutine create_Hypre_grid

  ! ---------------------------------------------------------------------------

  subroutine create_Hypre_A_matrix(MPI_COMM_WORLD, hypre_grid, hypre_A, &
            a, xlow, xhigh, ylow, yhigh, ierr)
    implicit none

    integer,          intent(in)  :: MPI_COMM_WORLD
    integer*8,        intent(in)  :: hypre_grid
    integer*8,        intent(out) :: hypre_A
    double precision, intent(in)  :: a(5, xlow:xhigh, ylow:yhigh)
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh
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

    do i = xlow, xhigh
      do j = ylow, yhigh
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
      etanew, nx, ny, xlow, xhigh, ylow, yhigh, OL, dt, maxits, eps, ierr)
    implicit none

    integer,          intent(in)  :: MPI_COMM_WORLD
    integer*8,        intent(in)  :: hypre_A
    integer*8,        intent(in)  :: hypre_grid
    integer,          intent(in)  :: myid
    integer,          intent(in)  :: num_procs
    integer,          intent(in)  :: ilower(0:num_procs-1,2)
    integer,          intent(in)  :: iupper(0:num_procs-1,2)
    double precision, intent(in)  :: etastar(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(out) :: etanew(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    integer,          intent(in)  :: nx, ny, xlow, xhigh, ylow, yhigh, OL
    double precision, intent(in)  :: dt
    integer,          intent(in)  :: maxits
    double precision, intent(in)  :: eps
    integer,          intent(out) :: ierr

    integer          :: i, j ! loop variables
    integer*8        :: hypre_b
    integer*8        :: hypre_x
    integer*8        :: hypre_solver
    integer*8        :: precond
    double precision :: b_values((xhigh-xlow+1)*(yhigh-ylow+1))
    double precision :: init_values((xhigh-xlow+1)*(yhigh-ylow+1))
    double precision :: final_values((xhigh-xlow+1)*(yhigh-ylow+1))

    integer :: nx_tile, ny_tile
    integer :: i_local, j_local

        ! A currently unused variable that can be used to
    ! print information from the solver - see comments below.
   ! double precision :: hypre_out(2)

    nx_tile = xhigh-xlow + 1
    ny_tile = yhigh-ylow + 1

    ! wrap this code in preprocessing flags to allow the model to be compiled without the external library, if desired.
#ifdef useExtSolver
    ! Create the rhs vector, b
    call HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid, hypre_b, ierr)
    call HYPRE_StructVectorInitialize(hypre_b, ierr)

    ! set rhs values (vector b)
    do j = ylow, yhigh ! loop over every grid point
      do i = xlow, xhigh

          i_local = i - xlow
          j_local = j - ylow

     ! the 2D array is being laid out like
     ! [x1y1, x2y1, x3y1, x1y2, x2y2, x3y2, x1y3, x2y3, x3y3]

          b_values(i_local &
                     + j_local*(xhigh - xlow + 1) &
                     + 1) = -etastar(i,j)/dt**2

          init_values(i_local &
                     + j_local*(xhigh - xlow + 1) &
                     + 1) = etastar(i,j)

    ! the 2D array is being laid out like
    ! [x1y1, x2y1, x3y1, x1y2, x2y2, x3y2, x1y3, x2y3, x3y3]
      end do
    end do

    call HYPRE_StructVectorSetBoxValues(hypre_b, &
      ilower(myid,:), iupper(myid,:), b_values, ierr)

    call HYPRE_StructVectorAssemble(hypre_b, ierr)

    ! now create the x vector
    call HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid, hypre_x, ierr)
    call HYPRE_StructVectorInitialize(hypre_x, ierr)

    call HYPRE_StructVectorSetBoxValues(hypre_x, &
      ilower(myid,:), iupper(myid,:), init_values, ierr)

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
      ilower(myid,:), iupper(myid,:), final_values, ierr)

    do j = ylow, yhigh ! loop over every grid point
      do i = xlow, xhigh
        ! etanew(i,j) = final_values( ((j-1)*nx_tile + i) )
          i_local = i - xlow
          j_local = j - ylow

     ! the 2D array is being laid out like
     ! [x1y1, x2y1, x3y1, x1y2, x2y2, x3y2, x1y3, x2y3, x3y3]

          etanew(i,j) = final_values(i_local &
                                     + j_local*(xhigh - xlow + 1) &
                                     + 1)

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
      xstep, ystep, dspace, dt, xlow, xhigh, ylow, yhigh, layers, OL)
    implicit none

    double precision, intent(inout) :: array(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in) :: etanew(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in) :: g_vec(layers)
    integer, intent(in) :: xstep, ystep
    double precision, intent(in) :: dspace, dt
    integer, intent(in) :: xlow, xhigh, ylow, yhigh, layers, OL

    integer i, j, k
    double precision baro_contrib

    ! TODO Assert that xstep and ystep are either 1, 0 or 0, 1.

    do k = 1, layers
      do j = ylow-OL+ystep, yhigh+OL-1
        do i = xlow-OL+xstep, xhigh+OL-1

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

  subroutine calc_A_matrix(a, depth, g, dx, dy, OL, &
            xlow, xhigh, ylow, yhigh, freesurf_fac, dt, &
            hfac_w, hfac_e, hfac_s, hfac_n)
    implicit none

    double precision, intent(out) :: a(5, xlow:xhigh, ylow:yhigh)
    double precision, intent(in)  :: depth(xlow-OL:xhigh+OL, &
                                            ylow-OL:yhigh+OL)
    double precision, intent(in)  :: g, dx, dy
    integer,          intent(in)  :: OL, xlow, xhigh, ylow, yhigh
    double precision, intent(in)  :: freesurf_fac
    double precision, intent(in)  :: dt
    double precision, intent(in)  :: hfac_w(xlow-OL:xhigh+OL, &
                                            ylow-OL:yhigh+OL)
    double precision, intent(in)  :: hfac_e(xlow-OL:xhigh+OL, &
                                            ylow-OL:yhigh+OL)
    double precision, intent(in)  :: hfac_n(xlow-OL:xhigh+OL, &
                                            ylow-OL:yhigh+OL)
    double precision, intent(in)  :: hfac_s(xlow-OL:xhigh+OL, &
                                            ylow-OL:yhigh+OL)

    integer i, j

    do j = ylow, yhigh
      do i = xlow, xhigh
        a(1,i,j) = g*0.5*(depth(i+1,j)+depth(i,j))*hfac_e(i,j)/dx**2
        a(2,i,j) = g*0.5*(depth(i,j+1)+depth(i,j))*hfac_n(i,j)/dy**2
        a(3,i,j) = g*0.5*(depth(i,j)+depth(i-1,j))*hfac_w(i,j)/dx**2
        a(4,i,j) = g*0.5*(depth(i,j)+depth(i,j-1))*hfac_s(i,j)/dy**2
      end do
    end do

    do j = ylow, yhigh
      do i = xlow, xhigh
        a(5,i,j) = -a(1,i,j)-a(2,i,j)-a(3,i,j)-a(4,i,j) - freesurf_fac/dt**2
      end do
    end do

    return
  end subroutine calc_A_matrix

  ! ---------------------------------------------------------------------------
  !> Do the isopycnal layer physics

  subroutine barotropic_correction(hnew, unew, vnew, eta, etanew, depth, a, &
      dx, dy, wetmask, hfac_w, hfac_s, dt, &
      maxits, eps, rjac, freesurf_fac, thickness_error, hmin, &
      debug_level, g_vec, nx, ny, layers, OL, &
      xlow, xhigh, ylow, yhigh, &
      n, MPI_COMM_WORLD, myid, num_procs, ilower, iupper, &
      hypre_grid, hypre_A, ierr)

    implicit none

    double precision, intent(inout) :: hnew(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(inout) :: unew(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(inout) :: vnew(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)    :: eta(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(out)   :: etanew(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: depth(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: a(5, nx, ny)
    double precision, intent(in)    :: dx, dy
    double precision, intent(in)    :: wetmask(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: hfac_w(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: hfac_s(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision, intent(in)    :: dt
    integer,          intent(in)    :: maxits
    double precision, intent(in)    :: eps, rjac, freesurf_fac, thickness_error, hmin
    integer,          intent(in)    :: debug_level
    double precision, intent(in)    :: g_vec(layers)
    integer,          intent(in)    :: nx, ny, layers, OL
    integer,          intent(in)    :: xlow, xhigh, ylow, yhigh, n
    integer,          intent(in)    :: MPI_COMM_WORLD
    integer,          intent(in)    :: myid, num_procs
    integer,          intent(in)    :: ilower(0:num_procs-1,2)
    integer,          intent(in)    :: iupper(0:num_procs-1,2)
    integer*8,        intent(in)    :: hypre_grid
    integer*8,        intent(in)    :: hypre_A
    integer,          intent(out) :: ierr

    ! barotropic velocity components (for pressure solver)
    double precision :: ub(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision :: vb(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    double precision :: etastar(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)

    character(10)    :: num

    ! Calculate the barotropic velocities
    call calc_baro_u(ub, unew, hnew, eta, freesurf_fac, xlow, xhigh, ylow, yhigh, &
                      layers, OL)
    call calc_baro_v(vb, vnew, hnew, eta, freesurf_fac, xlow, xhigh, ylow, yhigh, &
                      layers, OL)

    if (debug_level .ge. 4) then
      ! Output the data to a file
      call write_output_3d(ub, nx, ny, 1, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 1, 0, &
                          n, 'output/snap.ub.', num_procs, myid)

      call write_output_3d(vb, nx, ny, 1, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 1, &
                          n, 'output/snap.vb.', num_procs, myid)
    end if

    ! Calculate divergence of ub and vb, and solve for the pressure
    ! field that removes it
    call calc_eta_star(ub, vb, eta, etastar, freesurf_fac, &
                    xlow, xhigh, ylow, yhigh, OL, dx, dy, dt)
    ! print *, maxval(abs(etastar))
    if (debug_level .ge. 4) then
      call write_output_3d(etastar, nx, ny, 1, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 0, &
                          n, 'output/snap.eta_star.', num_procs, myid)
    end if

    ! Prevent barotropic signals from bouncing around outside the
    ! wet region of the model.
    ! etastar = etastar*wetmask
#ifndef useExtSolver
    ! this call IS NOT PARALLEL SAFE
    ! It has been left with global indicies, since it only works on one core.
    if (num_procs .ne. 1) then
      print "(A)", "SOR solver can only be used on one core."
      call clean_stop(n, .FALSE.)
    end if
    call SOR_solver(a, etanew, etastar, nx, ny, OL, &
       dt, rjac, eps, maxits, n)
    ! print *, maxval(abs(etanew))
#endif

#ifdef useExtSolver
    call Ext_solver(MPI_COMM_WORLD, hypre_A, hypre_grid, myid, num_procs, &
      ilower, iupper, etastar, &
      etanew, nx, ny, xlow, xhigh, ylow, yhigh, OL, dt, maxits, eps, ierr)
#endif

    etanew = etanew*wetmask

    if (debug_level .ge. 4) then
      call write_output_3d(etanew, nx, ny, 1, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, 0, 0, &
                          n, 'output/snap.eta_new.', num_procs, myid)
    end if


    call update_halos(etanew, nx, ny, 1, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)

    ! Now update the velocities using the barotropic tendency due to
    ! the pressure gradient.
    call update_velocities_for_barotropic_tendency(unew, etanew, g_vec, &
        1, 0, dx, dt, xlow, xhigh, ylow, yhigh, layers, OL)
    call update_velocities_for_barotropic_tendency(vnew, etanew, g_vec, &
        0, 1, dy, dt, xlow, xhigh, ylow, yhigh, layers, OL)

    ! We now have correct velocities at the next time step, but the
    ! layer thicknesses were updated using the velocities without
    ! the barotropic pressure contribution. Force consistency
    ! between layer thicknesses and ocean depth by scaling
    ! thicknesses to agree with free surface.
    call enforce_depth_thickness_consistency(hnew, etanew, depth, &
        freesurf_fac, thickness_error, hmin, xlow, xhigh, ylow, yhigh, layers, OL)

    ! Apply the boundary conditions
    call apply_boundary_conditions(unew, hfac_w, wetmask, &
                                xlow, xhigh, ylow, yhigh, layers, OL)
    call apply_boundary_conditions(vnew, hfac_s, wetmask, &
                                xlow, xhigh, ylow, yhigh, layers, OL)

    return
  end subroutine barotropic_correction

end module barotropic_mode
