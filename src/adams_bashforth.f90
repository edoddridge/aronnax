module adams_bashforth
  implicit none

  contains

  
  ! ---------------------------------------------------------------------------
  !> A first-order Forward Euler algorithm
  !! This is not a good algorithm. Don't use it, except to show how
  !! how bad it is.

  subroutine ForwardEuler(X_new, dXdt, X, dt, nx, ny, layers, OL, AB_order)
    implicit none

    double precision, intent(out) :: X_new(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(inout) :: dXdt(1-OL:nx+OL, 1-OL:ny+OL, layers, AB_order)
    double precision, intent(in) :: X(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: dt
    integer,          intent(in) :: nx, ny, layers, OL, AB_order

    X_new = X + dt*dXdt(:,:,:,1)

  end subroutine ForwardEuler

! ---------------------------------------------------------------------------
  !> A second-order Adams-Bashforth algorithm
  subroutine AB2(X_new, dXdt, X, dt, nx, ny, layers, OL, AB_order)
    implicit none

    double precision, intent(out) :: X_new(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(inout) :: dXdt(1-OL:nx+OL, 1-OL:ny+OL, layers, AB_order)
    double precision, intent(in) :: X(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: dt
    integer,          intent(in) :: nx, ny, layers, OL, AB_order

    X_new = X + dt*(3d0*dXdt(:,:,:,1) - 1d0*dXdt(:,:,:,2))/2d0

    ! Cycle tendencies
    dXdt(:,:,:,2) = dXdt(:,:,:,1)
    
  end subroutine AB2

! ---------------------------------------------------------------------------
  !> A third-order Adams-Bashforth algorithm
  subroutine AB3(X_new, dXdt, X, dt, nx, ny, layers, OL, AB_order)
    implicit none

    double precision, intent(out) :: X_new(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(inout) :: dXdt(1-OL:nx+OL, 1-OL:ny+OL, layers, AB_order)
    double precision, intent(in) :: X(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: dt
    integer,          intent(in) :: nx, ny, layers, OL, AB_order

    X_new = X + dt*(23d0*dXdt(:,:,:,1) - 16d0*dXdt(:,:,:,2) + &
                    5d0*dXdt(:,:,:,3))/12d0

    ! Cycle tendencies
    dXdt(:,:,:,3) = dXdt(:,:,:,2)
    dXdt(:,:,:,2) = dXdt(:,:,:,1)
    
  end subroutine AB3

! ---------------------------------------------------------------------------
  !> A fourth-order Adams-Bashforth algorithm
  subroutine AB4(X_new, dXdt, X, dt, nx, ny, layers, OL, AB_order)
    implicit none

    double precision, intent(out) :: X_new(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(inout) :: dXdt(1-OL:nx+OL, 1-OL:ny+OL, layers, AB_order)
    double precision, intent(in) :: X(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: dt
    integer,          intent(in) :: nx, ny, layers, OL, AB_order

    X_new = X + dt*(55d0*dXdt(:,:,:,1) - 59d0*dXdt(:,:,:,2) + &
                    37d0*dXdt(:,:,:,3) - 9d0*dXdt(:,:,:,4))/24d0

    ! Cycle tendencies
    dXdt(:,:,:,4) = dXdt(:,:,:,3)
    dXdt(:,:,:,3) = dXdt(:,:,:,2)
    dXdt(:,:,:,2) = dXdt(:,:,:,1)
    
  end subroutine AB4


! ---------------------------------------------------------------------------
  !> A fifth-order Adams-Bashforth algorithm
  subroutine AB5(X_new, dXdt, X, dt, nx, ny, layers, OL, AB_order)
    implicit none

    double precision, intent(out) :: X_new(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(inout) :: dXdt(1-OL:nx+OL, 1-OL:ny+OL, layers, AB_order)
    double precision, intent(in) :: X(1-OL:nx+OL, 1-OL:ny+OL, layers)
    double precision, intent(in) :: dt
    integer,          intent(in) :: nx, ny, layers, OL, AB_order

    X_new = X + dt*(1901d0*dXdt(:,:,:,1) - 2774d0*dXdt(:,:,:,2) + &
                    2616d0*dXdt(:,:,:,3) - 1274d0*dXdt(:,:,:,4) + &
                    251d0*dXdt(:,:,:,5))/720d0

    ! Cycle tendencies
    dXdt(:,:,:,5) = dXdt(:,:,:,4)
    dXdt(:,:,:,4) = dXdt(:,:,:,3)
    dXdt(:,:,:,3) = dXdt(:,:,:,2)
    dXdt(:,:,:,2) = dXdt(:,:,:,1)
    
  end subroutine AB5

end module adams_bashforth
