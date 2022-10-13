module declarations

  implicit none

  integer, parameter :: layerwise_input_length = 10000
  ! Resolution
  integer :: nx ! number of x grid points
  integer :: ny ! number of y grid points
  integer :: layers ! number of active layers in the model
  integer :: OL ! size of halo region
  ! Layer thickness (h)
  double precision, dimension(:,:,:), allocatable :: h
  ! Velocity component (u)
  double precision, dimension(:,:,:), allocatable :: u
  ! Velocity component (v)
  double precision, dimension(:,:,:), allocatable :: v
  ! Free surface (eta)
  double precision, dimension(:,:),   allocatable :: eta
  ! Bathymetry
  character(60) :: depth_file
  double precision, dimension(:,:),   allocatable :: depth
  double precision :: h0 ! default depth in no file specified
  ! Grid
  double precision :: dx, dy
  double precision, dimension(:,:),   allocatable :: wetmask
  double precision, dimension(:,:),   allocatable :: hfac_w
  double precision, dimension(:,:),   allocatable :: hfac_e
  double precision, dimension(:,:),   allocatable :: hfac_s
  double precision, dimension(:,:),   allocatable :: hfac_n

  ! Coriolis parameter at u and v grid-points respectively
  double precision, dimension(:,:),   allocatable :: fu
  double precision, dimension(:,:),   allocatable :: fv
  ! File names to read them from
  character(60) :: f_u_file, f_v_file
  character(60) :: wet_mask_file
  ! Numerics
  double precision :: dt ! delta_t in seconds
  double precision :: au ! viscosity
  double precision :: ar ! linear drag between layers
  double precision :: bot_drag ! linear bottom drag
  double precision :: kh(layerwise_input_length) ! horizontal thickness diffusivity
  double precision :: kv ! vertical thickness diffusivity
  double precision :: slip ! tangential momentum boundary condition
  double precision :: hmin ! minimum layer thickness
  integer          :: niter0 ! timestep to start from
  integer          :: n_time_steps ! timesteps to simulate
  double precision :: dump_freq ! time period between snapshot outputs
  double precision :: av_freq ! time period between averaged outputs
  double precision :: checkpoint_freq ! time period between checkpoint outputs
  double precision :: diag_freq ! time period between diagnostic outputs
  double precision, dimension(:), allocatable :: zeros
  integer          :: maxits ! maximum iterations for pressure solver
  double precision :: eps ! tolerance for pressure solver
  double precision :: freesurf_fac ! controls rigid lid or free surface
  double precision :: thickness_error ! max error between layer thickness
  ! and depth before a warning is printed
  integer          :: debug_level ! how much output should there be?
  integer          :: h_advec_scheme ! selects thickness advection scheme
  integer          :: ts_algorithm ! selects timestepping algorithm
  integer          :: AB_order ! used to construct tendency arrays

  ! Model
  ! shortcut for flat initial conditions
  double precision :: hmean(layerwise_input_length)
  ! Switch for using n + 1/2 layer physics, or using n layer physics
  logical :: red_grav
  ! Physics
  double precision :: g_vec(layerwise_input_length) ! gravity between layers
  double precision :: rho0 ! background density
  ! Wind
  double precision, dimension(:,:,:),   allocatable :: base_wind_x
  double precision, dimension(:,:,:),   allocatable :: base_wind_y
  logical          :: dump_wind
  character(60)    :: wind_mag_time_series_file
  double precision :: wind_depth ! depth over which the wind forcing is spread
  double precision, dimension(:),     allocatable :: wind_mag_time_series
  integer          :: wind_n_records
  double precision :: wind_period
  logical          :: wind_loop_fields
  logical          :: wind_interpolate

  ! Sponge regions
  double precision, dimension(:,:,:), allocatable :: sponge_h_time_scale
  double precision, dimension(:,:,:), allocatable :: sponge_u_time_scale
  double precision, dimension(:,:,:), allocatable :: sponge_v_time_scale
  double precision, dimension(:,:,:), allocatable :: sponge_h
  double precision, dimension(:,:,:), allocatable :: sponge_u
  double precision, dimension(:,:,:), allocatable :: sponge_v
  character(60) :: sponge_h_time_scale_file
  character(60) :: sponge_u_time_scale_file
  character(60) :: sponge_v_time_scale_file
  character(60) :: sponge_h_file
  character(60) :: sponge_u_file
  character(60) :: sponge_v_file
  ! Main input files
  character(60) :: init_u_file, init_v_file, init_h_file, init_eta_file
  character(60) :: zonal_wind_file, meridional_wind_file
  logical       :: relative_wind
  double precision :: Cd

  integer*8 :: start_time

  ! External pressure solver variables
  integer :: n_proc_x, n_proc_y


  integer :: ierr
  integer :: num_procs, myid

  integer, dimension(:),   allocatable :: xlower, xupper
  integer, dimension(:),   allocatable :: ylower, yupper
  integer, dimension(:,:), allocatable :: ilower, iupper
  integer, dimension(:,:), allocatable :: jlower, jupper
  integer*8 :: hypre_grid
  integer   :: i, j, k
  integer   :: offsets(2,5)

  contains

end module declarations
