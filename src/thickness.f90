module thickness
  use end_run
  use boundaries
  use advection_schemes
  implicit none

  contains

  ! ---------------------------------------------------------------------------
  !> Calculate the tendency of layer thickness for each of the active layers
  !! dh/dt is in the centre of each grid point.
  subroutine evaluate_dhdt(dhdt, h, u, v, kh, kv, dx, dy, nx, ny, layers, &
      spongeTimeScale, spongeH, wetmask, RedGrav, hAdvecScheme, n)
    implicit none

    ! dhdt is evaluated at the centre of the grid box
    double precision, intent(out) :: dhdt(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: u(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: v(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: kh(layers), kv
    double precision, intent(in)  :: dx, dy
    integer, intent(in) :: nx, ny, layers
    double precision, intent(in)  :: spongeTimeScale(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: spongeH(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: wetmask(0:nx+1, 0:ny+1)
    logical, intent(in) :: RedGrav
    integer, intent(in) :: hAdvecScheme
    integer, intent(in) :: n

    integer i, j, k
    ! Thickness tendency due to thickness diffusion (equivalent to Gent
    ! McWilliams in a z coordinate model)
    double precision dhdt_kh(0:nx+1, 0:ny+1, layers)

    ! Thickness tendency due to vertical diffusion of mass
    double precision dhdt_kv(0:nx+1, 0:ny+1, layers)

    ! Thickness tendency due to thickness advection
    double precision dhdt_advec(0:nx+1, 0:ny+1, layers)

    ! Calculate tendency due to thickness diffusion (equivalent
    ! to GM in z coordinate model with the same diffusivity).
    dhdt_kh = 0d0

    call dhdt_hor_diff(dhdt_kh, h, kh, dx, dy, nx, ny, layers, &
      wetmask, RedGrav)

    ! Calculate thickness tendency due to vertical diffusion
    dhdt_kv = 0d0
    call dhdt_vert_diff(dhdt_kv, h, kv, nx, ny, layers, RedGrav)


    ! Calculate the thickness tendency due to the flow field
    dhdt_advec = 0d0

    if (hAdvecScheme .eq. 1) then
      ! first-order centered
      call h_advec_1_centered(dhdt_advec, h, u, v, dx, dy, nx, ny, layers)
    else if (hAdvecScheme .eq. 2) then
      ! first-order upwind
      call h_advec_1_upwind(dhdt_advec, h, u, v, dx, dy, nx, ny, layers)
    else
      call clean_stop(n, .FALSE.)
    end if

    ! Now add these together, along with the sponge contribution
    dhdt = 0d0

    do k = 1, layers
      do j = 1, ny
        do i = 1, nx
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

    call wrap_fields_3D(dhdt, nx, ny, layers)

    return
  end subroutine evaluate_dhdt

  ! ---------------------------------------------------------------------------
  !> Calculate the tendency of layer thickness for each of the active layers
  !! due to horizontal thickness diffusion
  subroutine dhdt_hor_diff(dhdt_GM, h, kh, dx, dy, nx, ny, layers, &
      wetmask, RedGrav)
    implicit none

    ! dhdt is evaluated at the centre of the grid box
    double precision, intent(out) :: dhdt_GM(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: kh(layers)
    double precision, intent(in)  :: dx, dy
    integer, intent(in) :: nx, ny, layers
    double precision, intent(in)  :: wetmask(0:nx+1, 0:ny+1)
    logical, intent(in) :: RedGrav

    integer i, j, k

    ! Calculate tendency due to thickness diffusion (equivalent
    ! to GM in z coordinate model with the same diffusivity).
    dhdt_GM = 0d0

    ! Loop through all layers except lowest and calculate
    ! thickness tendency due to horizontal diffusive mass fluxes
    do k = 1, layers-1
      do j = 1, ny
        do i = 1, nx
          dhdt_GM(i,j,k) = &
              kh(k)*(h(i+1,j,k)*wetmask(i+1,j)    &
                + (1d0 - wetmask(i+1,j))*h(i,j,k) & ! reflect around boundary
                + h(i-1,j,k)*wetmask(i-1,j)       &
                + (1d0 - wetmask(i-1,j))*h(i,j,k) & ! refelct around boundary
                - 2*h(i,j,k))/(dx*dx)             & ! x-component

              + kh(k)*(h(i,j+1,k)*wetmask(i,j+1) &
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
      do j = 1, ny
        do i = 1, nx
          dhdt_GM(i,j,layers) = &
              kh(layers)*(h(i+1,j,layers)*wetmask(i+1,j)   &
                + (1d0 - wetmask(i+1,j))*h(i,j,layers)     & ! boundary
                + h(i-1,j,layers)*wetmask(i-1,j)           &
                + (1d0 - wetmask(i-1,j))*h(i,j,layers)     & ! boundary
                - 2*h(i,j,layers))/(dx*dx)                 & ! x-component

              + kh(layers)*(h(i,j+1,layers)*wetmask(i,j+1) &
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
  subroutine dhdt_vert_diff(dhdt_kv, h, kv, nx, ny, layers, RedGrav)
    implicit none

    ! dhdt is evaluated at the centre of the grid box
    double precision, intent(out) :: dhdt_kv(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: h(0:nx+1, 0:ny+1, layers)
    double precision, intent(in)  :: kv
    integer, intent(in) :: nx, ny, layers
    logical, intent(in) :: RedGrav

    integer i, j, k

    dhdt_kv = 0d0

    ! calculate vertical diffusive mass fluxes
    ! only evaluate vertical mass diff flux if more than 1 layer, or reduced gravity
    if (layers .eq. 1) then
      if (RedGrav) then
        do j = 1, ny
          do i = 1, nx
            dhdt_kv(i,j,1) = kv/h(i,j,1)
          end do
        end do
      end if
    else if (layers .gt. 1) then
      ! if more than one layer, need to have multiple fluxes
      do k = 1, layers
        do j = 1, ny
          do i = 1, nx
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
