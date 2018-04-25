module declarations

  implicit none

  integer, parameter :: layerwise_input_length = 10000
  ! Resolution
  integer :: nx ! number of x grid points
  integer :: ny ! number of y grid points
  integer :: layers ! number of active layers in the model
  ! Layer thickness (h)
  double precision, dimension(:,:,:), allocatable :: h
  ! Velocity component (u)
  double precision, dimension(:,:,:), allocatable :: u
  ! Velocity component (v)
  double precision, dimension(:,:,:), allocatable :: v
  ! Free surface (eta)
  double precision, dimension(:,:),   allocatable :: eta
  ! Bathymetry
  character(60) :: depthFile
  double precision, dimension(:,:),   allocatable :: depth
  double precision :: H0 ! default depth in no file specified
  ! Grid
  double precision :: dx, dy
  double precision, dimension(:,:),   allocatable :: wetmask
  ! Coriolis parameter at u and v grid-points respectively
  double precision, dimension(:,:),   allocatable :: fu
  double precision, dimension(:,:),   allocatable :: fv
  ! File names to read them from
  character(60) :: fUfile, fVfile
  character(60) :: wetMaskFile
  ! Numerics
  double precision :: dt ! delta_t in seconds
  double precision :: au ! viscosity
  double precision :: ar ! linear drag between layers
  double precision :: botDrag ! linear bottom drag
  double precision :: kh(layerwise_input_length) ! horizontal thickness diffusivity
  double precision :: kv ! vertical thickness diffusivity
  double precision :: slip ! tangential momentum boundary condition
  double precision :: hmin ! minimum layer thickness
  integer          :: niter0 ! timestep to start from 
  integer          :: nTimeSteps ! timesteps to simulate
  double precision :: dumpFreq ! time period between snapshot outputs
  double precision :: avFreq ! time period between averaged outputs
  double precision :: checkpointFreq ! time period between checkpoint outputs
  double precision :: diagFreq ! time period between diagnostic outputs
  double precision, dimension(:), allocatable :: zeros
  integer          :: maxits ! maximum iterations for pressure solver
  double precision :: eps ! tolerance for pressure solver
  double precision :: freesurfFac ! controls rigid lid or free surface
  double precision :: thickness_error ! max error between layer thickness
  ! and depth before a warning is printed
  integer          :: debug_level ! how much output should there be?
  integer          :: hAdvecScheme ! selects thickness advection scheme
  integer          :: TS_algorithm ! selects timestepping algorithm
  integer          :: AB_order ! used to construct tendency arrays

  ! Model
  ! shortcut for flat initial conditions
  double precision :: hmean(layerwise_input_length) 
  ! Switch for using n + 1/2 layer physics, or using n layer physics
  logical :: RedGrav
  ! Physics
  double precision :: g_vec(layerwise_input_length) ! gravity between layers
  double precision :: rho0 ! background density
  ! Wind
  double precision, dimension(:,:),   allocatable :: base_wind_x
  double precision, dimension(:,:),   allocatable :: base_wind_y
  logical       :: DumpWind
  character(60) :: wind_mag_time_series_file
  double precision, dimension(:),     allocatable :: wind_mag_time_series
  ! Sponge regions
  double precision, dimension(:,:,:), allocatable :: spongeHTimeScale
  double precision, dimension(:,:,:), allocatable :: spongeUTimeScale
  double precision, dimension(:,:,:), allocatable :: spongeVTimeScale
  double precision, dimension(:,:,:), allocatable :: spongeH
  double precision, dimension(:,:,:), allocatable :: spongeU
  double precision, dimension(:,:,:), allocatable :: spongeV
  character(60) :: spongeHTimeScaleFile
  character(60) :: spongeUTimeScaleFile
  character(60) :: spongeVTimeScaleFile
  character(60) :: spongeHfile
  character(60) :: spongeUfile
  character(60) :: spongeVfile
  ! Main input files
  character(60) :: initUfile, initVfile, initHfile, initEtaFile
  character(60) :: zonalWindFile, meridionalWindFile
  logical       :: RelativeWind
  double precision :: Cd

  ! External pressure solver variables
  integer :: nProcX, nProcY


  integer :: ierr
  integer :: num_procs, myid

  integer, dimension(:,:), allocatable :: ilower, iupper
  integer, dimension(:,:), allocatable :: jlower, jupper
  integer*8 :: hypre_grid
  integer   :: i, j
  integer   :: offsets(2,5)

  contains

end module declarations
