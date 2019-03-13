module thickness
  use end_run
  use boundaries
  use advection_schemes
  use exchange
  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Calculate the tendency of layer thickness for each of the active layers
  !! dh/dt is in the centre of each grid point.
  subroutine evaluate_dhdt(dhdt, h, u, v, kh, hmin, kv, dx, dy, xlow, xhigh, ylow, yhigh, &
      nx, ny, layers, OL, spongeTimeScale, spongeH, wetmask, RedGrav, hAdvecScheme, n, &
      ilower, iupper, num_procs, myid)
    implicit none

    ! dhdt is evaluated at the centre of the grid box
    double precision, intent(out) :: dhdt(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: u(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: v(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: kh(layers), kv
    double precision, intent(in)  :: hmin
    double precision, intent(in)  :: dx, dy
    integer,          intent(in)  :: xlow, xhigh, ylow, yhigh
    integer,          intent(in)  :: nx, ny, layers, OL
    double precision, intent(in)  :: spongeTimeScale(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: spongeH(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: wetmask(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    logical,          intent(in)  :: RedGrav
    integer,          intent(in)  :: hAdvecScheme
    integer,          intent(in)  :: n
    integer,          intent(in)  :: ilower(0:num_procs-1,2)
    integer,          intent(in)  :: iupper(0:num_procs-1,2)
    integer,          intent(in)  :: num_procs, myid

    integer i, j, k
    ! Thickness tendency due to thickness diffusion (equivalent to Gent
    ! McWilliams in a z coordinate model)
    double precision dhdt_kh(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)

    ! Thickness tendency due to vertical diffusion of mass
    double precision dhdt_kv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)

    ! Thickness tendency due to thickness advection
    double precision dhdt_advec(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)

    ! Calculate tendency due to thickness diffusion (equivalent
    ! to GM in z coordinate model with the same diffusivity).
    dhdt_kh = 0d0

    call dhdt_hor_diff(dhdt_kh, h, kh, hmin, dx, dy, xlow, xhigh, ylow, yhigh, layers, &
      OL, wetmask, RedGrav)

    ! Calculate thickness tendency due to vertical diffusion
    dhdt_kv = 0d0
    call dhdt_vert_diff(dhdt_kv, h, kv, xlow, xhigh, ylow, yhigh, layers, OL, RedGrav)


    ! Calculate the thickness tendency due to the flow field
    dhdt_advec = 0d0

    if (hAdvecScheme .eq. 1) then
      ! first-order centered
      call h_advec_1_centered(dhdt_advec, h, u, v, dx, dy, xlow, xhigh, ylow, yhigh, layers, OL)
    else if (hAdvecScheme .eq. 2) then
      ! first-order upwind
      call h_advec_1_upwind(dhdt_advec, h, u, v, dx, dy, xlow, xhigh, ylow, yhigh, layers, OL)
    else
      call clean_stop(n, .FALSE.)
    end if

    ! Now add these together, along with the sponge contribution
    dhdt = 0d0

    do k = 1, layers
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          dhdt(i,j,k) = &
              dhdt_kh(i,j,k) & ! horizontal thickness diffusion
              + dhdt_kv(i,j,k) & ! vetical thickness diffusion 
              + dhdt_advec(i,j,k) & ! thickness advection
              + spongeTimeScale(i,j,k)*(spongeH(i,j,k)-h(i,j,k)) ! forced relaxtion in the sponge regions.
        end do
      end do
    end do

    ! Make sure the dynamics are only happening in the wet grid points.
    do k = 1, layers
      dhdt(:, :, k) = dhdt(:, :, k) * wetmask
    end do

    ! call wrap_fields_3D(dhdt, nx, ny, layers, OL)
    call update_halos(dhdt, nx, ny, layers, ilower, iupper, &
                          xlow, xhigh, ylow, yhigh, OL, num_procs, myid)
    return
  end subroutine evaluate_dhdt

  ! ---------------------------------------------------------------------------
  !> Calculate the tendency of layer thickness for each of the active layers
  !! due to horizontal thickness diffusion
  subroutine dhdt_hor_diff(dhdt_GM, h, kh, hmin, dx, dy, xlow, xhigh, ylow, yhigh, layers, &
      OL, wetmask, RedGrav)
    implicit none

    ! dhdt is evaluated at the centre of the grid box
    double precision, intent(out) :: dhdt_GM(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: kh(layers)
    double precision, intent(in)  :: hmin
    double precision, intent(in)  :: dx, dy
    integer, intent(in) :: xlow, xhigh, ylow, yhigh, layers, OL
    double precision, intent(in)  :: wetmask(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL)
    logical, intent(in) :: RedGrav

    integer i, j, k
    double precision :: kappa_local

    ! Calculate tendency due to thickness diffusion (equivalent
    ! to GM in z coordinate model with the same diffusivity).
    dhdt_GM = 0d0

    ! Loop through all layers except lowest and calculate
    ! thickness tendency due to horizontal diffusive mass fluxes
    do k = 1, layers-1
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          ! if h >> hmin, then kappa_local = kh(k)
          ! if h <~ hmin, then kappa_local is very big
          kappa_local = kh(k) + 100d0*exp(-h(i,j,k)/hmin)
          dhdt_GM(i,j,k) = &
              kappa_local*(h(i+1,j,k)*wetmask(i+1,j)    &
                + (1d0 - wetmask(i+1,j))*h(i,j,k) & ! reflect around boundary
                + h(i-1,j,k)*wetmask(i-1,j)       &
                + (1d0 - wetmask(i-1,j))*h(i,j,k) & ! refelct around boundary
                - 2*h(i,j,k))/(dx*dx)             & ! x-component

              + kappa_local*(h(i,j+1,k)*wetmask(i,j+1) &
                + (1d0 - wetmask(i,j+1))*h(i,j,k) & ! reflect value around boundary
                + h(i,j-1,k)*wetmask(i,j-1)       &
                + (1d0 - wetmask(i,j-1))*h(i,j,k) & ! reflect value around boundary
                - 2*h(i,j,k))/(dy*dy)               ! y-component horizontal diffusion
        end do
      end do
    end do


    ! Now do the lowest active layer, k = layers. If using reduced
    ! gravity physics, this is unconstrained and calculated as above. If
    ! using n-layer physics it is constrained to balance the layers
    ! above it.
    if (RedGrav) then
      do j = ylow-OL+1, yhigh+OL-1
        do i = xlow-OL+1, xhigh+OL-1
          ! if h >> hmin, then kappa_local = kh(k)
          ! if h <~ hmin, then kappa_local is very big
          kappa_local = kh(layers) + 100d0*exp(-h(i,j,layers)/hmin)
          dhdt_GM(i,j,layers) = &
              kappa_local*(h(i+1,j,layers)*wetmask(i+1,j)   &
                + (1d0 - wetmask(i+1,j))*h(i,j,layers)     & ! boundary
                + h(i-1,j,layers)*wetmask(i-1,j)           &
                + (1d0 - wetmask(i-1,j))*h(i,j,layers)     & ! boundary
                - 2*h(i,j,layers))/(dx*dx)                 & ! x-component

              + kappa_local*(h(i,j+1,layers)*wetmask(i,j+1) &
                + (1d0 - wetmask(i,j+1))*h(i,j,layers)     & ! reflect value around boundary
                + h(i,j-1,layers)*wetmask(i,j-1)           &
                + (1d0 - wetmask(i,j-1))*h(i,j,layers)     & ! reflect value around boundary
                - 2*h(i,j,layers))/(dy*dy) ! y-component horizontal diffusion
        end do
      end do
    else if (.not. RedGrav) then ! using n-layer physics
      ! Calculate bottom layer thickness tendency to balance layers above.
      ! In the flat-bottomed case this will give the same answer.
      dhdt_GM(:,:,layers) = -sum(dhdt_GM(:,:,:layers-1), 3)
    end if

    return
  end subroutine dhdt_hor_diff

  ! ---------------------------------------------------------------------------
  !> Calculate the tendency of layer thickness for each of the active layers
  !! due to vertical thickness diffusion
  subroutine dhdt_vert_diff(dhdt_kv, h, kv, xlow, xhigh, ylow, yhigh, layers, OL, RedGrav)
    implicit none

    ! dhdt is evaluated at the centre of the grid box
    double precision, intent(out) :: dhdt_kv(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: h(xlow-OL:xhigh+OL, ylow-OL:yhigh+OL, layers)
    double precision, intent(in)  :: kv
    integer, intent(in) :: xlow, xhigh, ylow, yhigh, layers, OL
    logical, intent(in) :: RedGrav

    integer i, j, k

    dhdt_kv = 0d0

    ! calculate vertical diffusive mass fluxes
    ! only evaluate vertical mass diff flux if more than 1 layer, or reduced gravity
    if (layers .eq. 1) then
      if (RedGrav) then
        do j = ylow, yhigh
          do i = xlow, xhigh
            dhdt_kv(i,j,1) = kv/h(i,j,1)
          end do
        end do
      end if
    else if (layers .gt. 1) then
      ! if more than one layer, need to have multiple fluxes
      do k = 1, layers
        do j = ylow, yhigh
          do i = xlow, xhigh
            if (k .eq. 1) then ! in top layer
              dhdt_kv(i,j,k) = kv/h(i,j,k) - kv/h(i,j,k+1)
            else if (k .eq. layers) then ! bottom layer
              dhdt_kv(i,j,k) = kv/h(i,j,k) - kv/h(i,j,k-1)
            else ! mid layer/s
              dhdt_kv(i,j,k) = 2d0*kv/h(i,j,k) -  &
                  kv/h(i,j,k-1) - kv/h(i,j,k+1)
            end if
          end do
        end do
      end do
    end if

    return
  end subroutine dhdt_vert_diff

end module thickness
