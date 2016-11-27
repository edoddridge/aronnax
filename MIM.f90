!> @author
!> Ed Doddridge
!
!> Minimalist Isopycnal Model (MIM) with n layers
!!
!
!>       
!>     @mainpage Documentation for MIM.f90
!>
!>     @section Overview
!>     This model is an isopycnal model on an Arakawa C-grid with n
!>     layers and a rigid lid.
!>
!>
!>            
!>    @section Grid
!>
!>    /\ ------------
!>    |  |          |
!>    |  |          |
!>    dy U    H     |
!>    |  |          |
!>    |  |          |
!>    \/ Z----V------
!>        <---dx---->
!>
!>    H: tracer point - thickness, Bernoulli potential
!>    U: velocity point - u and v
!>    Z: vorticity point - zeta
!>

program MIM
!     
!
    implicit none

!> number of x and y grid points
    integer nx,ny !< number of x and y grid points
    integer layers !< number of active layers in the model
    parameter(nx=300,ny=300)

!   number of active layers
    parameter(layers = 2)

!     layer thickness (h), velocity components (u,v), free surface (eta)
    double precision h(0:nx,0:ny,layers)
    double precision u(nx,0:ny,layers),v(0:nx,ny,layers)
    double precision eta(0:nx,0:ny)
    double precision dhdt(0:nx,0:ny,layers)
    double precision dudt(nx,0:ny,layers),dvdt(0:nx,ny,layers) 
    double precision dudt_bt(nx,0:ny),dvdt_bt(0:nx,ny) 
    double precision dhdtold(0:nx,0:ny,layers), dudtold(nx,0:ny,layers),dvdtold(0:nx,ny,layers)         
    double precision dhdtveryold(0:nx,0:ny,layers), dudtveryold(nx,0:ny,layers),dvdtveryold(0:nx,ny,layers) 
    double precision hnew(0:nx,0:ny,layers)
    double precision unew(nx,0:ny,layers),vnew(0:nx,ny,layers)
    double precision etastar(0:nx,0:ny) 
    double precision etanew(0:nx,0:ny)
! barotropic velocity components (for pressure solver)
    double precision ub(nx,0:ny),vb(0:nx,ny)
!    Arrays for initialisation
    double precision hhalf(0:nx, 0:ny,layers), uhalf(nx,0:ny,layers), vhalf(0:nx,ny,layers)
!    for saving average fields
    double precision hav(0:nx,0:ny,layers)
    double precision uav(nx,0:ny,layers), vav(0:nx,ny,layers)
    double precision etaav(0:nx,0:ny)
!   bathymetry
    character(30) depthFile
    double precision depth(0:nx,0:ny)
    double precision H0 ! default depth in no file specified
!   Pressure solver variables
    double precision a(5,nx-1,ny-1)
    double precision phi(0:nx,0:ny), phiold(0:nx,0:ny)
!
!     Bernoulli potential and relative vorticity 
    double precision b(nx,ny,layers),zeta(nx,ny,layers)

!    Grid
    double precision dx, dy
    double precision x_u(nx), y_u(0:ny)
    double precision x_v(0:nx), y_v(ny)
    double precision wetmask(0:nx,0:ny)
    double precision hfacW(nx,0:ny)
    double precision hfacE(nx,0:ny)
    double precision hfacN(0:nx,ny)
    double precision hfacS(0:nx,ny)
!    Coriolis parameter at u and v grid-points respectively
    double precision fu(nx,0:ny),fv(0:nx,ny)
    character(30) fUfile, fVfile
    character(30) wetMaskFile

!   Numerics
    double precision pi,au,ah(layers),ar,dt
    double precision slip
    double precision hmin
    integer nTimeSteps
    double precision dumpFreq, avFreq
    integer counter
    integer nwrite, avwrite
    double precision zeros(layers)
    integer maxits
    double precision eps, rjac
    double precision freesurfFac
    double precision h_norming(0:nx,0:ny)

!   model
    double precision hmean(layers)
    logical :: RedGrav ! logical switch for using n + 1/2 layer physics, or using n layer physics

!   loop variables
    integer i,j,k,n

!   character variable for numbering the outputs
    character(10) num

!   Physics
    double precision g_vec(layers), rho0


!     Wind
    double precision wind_x(nx, 0:ny),wind_y(0:nx, ny)
    double precision base_wind_x(nx, 0:ny),base_wind_y(0:nx, ny)
    logical :: UseSinusoidWind
    logical :: UseStochWind
    logical :: DumpWind
    double precision :: wind_alpha, wind_beta, wind_period, wind_t_offset
    integer :: n_stoch_wind
    double precision :: stoch_wind_mag


!   Sponge
    double precision spongeHTimeScale(0:nx,0:ny,layers)
    double precision spongeUTimeScale(nx,0:ny,layers)
    double precision spongeVTimeScale(0:nx,ny,layers)
    double precision spongeH(0:nx,0:ny,layers)
    double precision spongeU(nx,0:ny,layers)
    double precision spongeV(0:nx,ny,layers)
    character(30) spongeHTimeScaleFile
    character(30) spongeUTimeScaleFile
    character(30) spongeVTimeScaleFile
    character(30) spongeHfile
    character(30) spongeUfile
    character(30) spongeVfile

!   input files
    character(30) initUfile,initVfile,initHfile,initEtaFile
    character(30) zonalWindFile,meridionalWindFile

!   Set default values here
!   Possibly wait until the model is split into multiple files, then hide the long unsightly code there.

    NAMELIST /NUMERICS/ au,ah,ar,dt,slip,nTimeSteps,dumpFreq, & 
                        avFreq,hmin, maxits, freesurfFac, eps

    NAMELIST /MODEL/ hmean, depthFile, H0, RedGrav

    NAMELIST /SPONGE/ spongeHTimeScaleFile,spongeUTimeScaleFile,spongeVTimeScaleFile, &
    spongeHfile,spongeUfile,spongeVfile
    
    NAMELIST /PHYSICS/ g_vec, rho0
    
    NAMELIST /GRID/ dx, dy, fUfile, fVfile, wetMaskFile
    
    NAMELIST /INITIAL_CONDITONS/ initUfile,initVfile,initHfile,initEtaFile

    NAMELIST /EXTERNAL_FORCING/ zonalWindFile,meridionalWindFile, &
        UseSinusoidWind, UseStochWind, wind_alpha, wind_beta, &
        wind_period, wind_t_offset, DumpWind


    open(unit=8,file="parameters.in", status='OLD', recl=80)
    read(unit=8,nml=NUMERICS)
    read(unit=8,nml=MODEL)
    read(unit=8,nml=SPONGE)
    read(8,nml=PHYSICS)
    read(8,nml=GRID)
    read(8,nml=INITIAL_CONDITONS)
    read(8,nml=EXTERNAL_FORCING)
    close (unit = 8)


    nwrite = int(dumpFreq/dt)
    avwrite = int(avFreq/dt)

!    Pi, the constant
    pi=3.14159
!   Zero vector - for internal use only
    zeros = 0d0

! check that time dependent forcing flags have been set correctly
    if (UseSinusoidWind .and. UseStochWind)  then
        
        ! write a file saying so
        OPEN(UNIT=99, FILE='errors.txt', ACTION="write", STATUS="replace", &
            FORM="formatted")
        write(99,*) "Can't have both stochastic and sinusoidally varying wind forcings. Choose one."
        CLOSE(UNIT=99)

        ! print it on the screen
        print *, "Can't have both stochastic and sinusoidally varying wind forcings. Choose one."
        ! Stop the code
        STOP
    endif

!   Read in arrays from the input files
    call read_input_fileU(initUfile,u,0.d0,nx,ny,layers)
    call read_input_fileV(initVfile,v,0.d0,nx,ny,layers)
    call read_input_fileH(initHfile,h,hmean,nx,ny,layers)

    if (.not. RedGrav) then
        call read_input_fileH_2D(depthFile,depth,H0,nx,ny)
        call read_input_fileH_2D(initEtaFile,eta,0.d0,nx,ny)
    endif
    ! Should probably check that bathymetry file and layer thicknesses are consistent with each other.

    call read_input_fileU(fUfile,fu,0.d0,nx,ny,1)
    call read_input_fileV(fVfile,fv,0.d0,nx,ny,1)

    call read_input_fileU(zonalWindFile,base_wind_x,0.d0,nx,ny,1)
    call read_input_fileV(meridionalWindFile,base_wind_y,0.d0,nx,ny,1)

    call read_input_fileH(spongeHTimeScaleFile,spongeHTimeScale, zeros,nx,ny,layers)
    call read_input_fileH(spongeHfile,spongeH,hmean,nx,ny,layers)
    call read_input_fileU(spongeUTimeScaleFile,spongeUTimeScale,0.d0,nx,ny,layers)
    call read_input_fileU(spongeUfile,spongeU,0.d0,nx,ny,layers)
    call read_input_fileV(spongeVTimeScaleFile,spongeVTimeScale,0.d0,nx,ny,layers)
    call read_input_fileV(spongeVfile,spongeV,0.d0,nx,ny,layers)
    call read_input_fileH_2D(wetMaskFile,wetmask,0.d0,nx,ny)


! Initialise the average fields
    if (avwrite .ne. 0) then
        hav = 0.0
        uav = 0.0
        vav = 0.0
        if (.not. RedGrav) then
            etaav = 0.0
        endif
    endif

!   For now the model can't do periodic boundaries - would be a good addition.
    if (wetMaskFile.eq.'') then
    ! 1 means ocean, 0 means land
        wetmask = 1d0 
        wetmask(0,:) = 0d0
        wetmask(:,0) = 0d0
        wetmask(nx,:) = 0d0
        wetmask(:,ny) = 0d0
    endif

    call calc_boundary_masks(wetmask,hfacW,hfacE,hfacS,hfacN,nx,ny)
     


!   If the winds are static, then set wind_ = base_wind_
    if (.not. UseSinusoidWind .and. .not. UseStochWind)  then
        wind_x = base_wind_x
        wind_y = base_wind_y
    endif


! initialise random numbers for stochastic winds
    if (UseStochWind) then
    ! this ensures a different series of random 
    ! numbers each time the model is run
        call ranseed()
    ! how many timesteps between updating the perturbed wind field.
        n_stoch_wind = int(wind_period/dt)
    endif

    if (.not. RedGrav) then
        !   initialise arrays for pressure solver 
        ! a =  derivatives of the depth field
        do j=1,ny-1
            do i=1,nx-1
                a(1,i,j)=g_vec(1)*0.5*(depth(i+1,j)+depth(i,j))/dx**2
                a(2,i,j)=g_vec(1)*0.5*(depth(i,j+1)+depth(i,j))/dy**2
                a(3,i,j)=g_vec(1)*0.5*(depth(i,j)+depth(i-1,j))/dx**2
                a(4,i,j)=g_vec(1)*0.5*(depth(i,j)+depth(i,j-1))/dy**2
            end do 
        end do
        do j=1,ny-1
            a(1,nx-1,j)=0.0
            a(3,1,j)=0.0
        end do
        do i=1,nx-1
            a(2,i,ny-1)=0.0
            a(4,i,1)=0.0
        end do
        do j=1,ny-1
            do i=1,nx-1
                a(5,i,j)=-a(1,i,j)-a(2,i,j)-a(3,i,j)-a(4,i,j)
            end do 
        end do

        ! calculate the spectral radius of the grid for use by the successive over relaxation scheme
        rjac=(cos(pi/real(nx))*dy**2+cos(pi/real(ny))*dx**2) & 
                /(dx**2+dy**2)
        ! if peridodic boundary conditions are ever implemented, then pi -> 2*pi in this calculation 

        ! Model currently only works with rigid lid (freesurfFac = 0)
        if (freesurfFac .ne. 0) then
            freesurfFac = 0d0
            print *, 'free surface not working yet. freesurfFac set to 0.'
        end if
        
        ! check that the supplied free surface anomaly and layer thicknesses are consistent
        h_norming = (freesurfFac*eta+depth)/sum(h,3)
        do k = 1,layers
            h(:,:,k) = h(:,:,k)*h_norming
        end do

        if (maxval(abs(h_norming - 1d0)).gt. 1e-2) then
            print *, 'inconsistency between h and eta (in %):', maxval(abs(h_norming - 1d0))*100d0
        end if
    endif

!    Do two initial time steps with Runge-Kutta second-order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!    Initialisation of the model STARTS HERE            !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! these initialisation steps do NOT use or update the free surface.
!
!--------------------------- negative 2 time step -----------------------------
!    Code to work out dhdtveryold, dudtveryold and dvdtveryold



!    calculate baroclinic Bernoulli potential
    if (RedGrav) then
        call evaluate_b_RedGrav(b,h,u,v,nx,ny,layers,g_vec)
    else
        call evaluate_b_iso(b,h,u,v,nx,ny,layers,g_vec)
    endif

!    calculate relative vorticity
    call evaluate_zeta(zeta,u,v,nx,ny,layers,dx,dy)

!     Calculate dhdt, dudt, dvdt at current time step
    call evaluate_dhdt(dhdtveryold, h,u,v,ah,dx,dy,nx,ny,layers,&
        spongeHTimeScale,spongeH,wetmask)

    call evaluate_dudt(dudtveryold, h,u,v,b,zeta,wind_x,fu, au,ar,&
        slip,dx,dy,hfacN,hfacS,nx,ny,layers,rho0,&
        spongeUTimeScale,spongeU)

    call evaluate_dvdt(dvdtveryold, h,u,v,b,zeta,wind_y,fv, &
        au,ar,slip,dx,dy,hfacW,hfacE,nx,ny,layers,rho0,spongeVTimeScale,spongeV)

!    Calculate the values at half the time interval with Forward Euler
    hhalf = h+0.5d0*dt*dhdtveryold
    uhalf = u+0.5d0*dt*dudtveryold
    vhalf = v+0.5d0*dt*dvdtveryold


!    calculate Bernoulli potential
    if (RedGrav) then
        call evaluate_b_RedGrav(b,hhalf,uhalf,vhalf,nx,ny,layers,g_vec)
    else
        call evaluate_b_iso(b,hhalf,uhalf,vhalf,nx,ny,layers,g_vec)
    endif

!    calculate relative vorticity
    call evaluate_zeta(zeta,uhalf,vhalf,nx,ny,layers,dx,dy)

!    Now calculate d/dt of u,v,h and store as dhdtveryold, dudtveryold and dvdtveryold
    call evaluate_dhdt(dhdtveryold, hhalf,uhalf,vhalf,ah,dx,dy,nx,ny,layers, spongeHTimeScale,spongeH,wetmask)

    call evaluate_dudt(dudtveryold, hhalf,uhalf,vhalf,b,zeta, &
        wind_x,fu, au,ar,slip,dx,dy,hfacN,hfacS,nx,ny,layers,rho0, &
        spongeUTimeScale,spongeU)

    call evaluate_dvdt(dvdtveryold, hhalf,uhalf,vhalf,b,zeta, &
        wind_y,fv, au,ar,slip ,dx,dy,hfacW,hfacE,nx,ny,layers,rho0,& 
        spongeVTimeScale,spongeV)
!    These are the values to be stored in the 'veryold' variables ready to start 
!    the proper model run.

!    Calculate h,u,v with these tendencies
    h = h + dt*dhdtveryold
    u = u + dt*dudtveryold
    v = v + dt*dvdtveryold

!    Apply the boundary conditions
!    Enforce no normal flow boundary condition
!    and no flow in dry cells.
!    no/free-slip is done inside dudt and dvdt subroutines.
!    hfacW and hfacS are zero where the transition between
!    wet and dry cells occurs. wetmask is 1 in wet cells,
!    and zero in dry cells.
    do k = 1,layers
        u(:,:,k) = u(:,:,k)*hfacW*wetmask(1:nx,:)
        v(:,:,k) = v(:,:,k)*hfacS*wetmask(:,1:ny)
    enddo



!--------------------------- negative 1 time step -----------------------------
!    Code to work out dhdtold, dudtold and dvdtold
!    calculate Bernoulli potential
    if (RedGrav) then
        call evaluate_b_RedGrav(b,h,u,v,nx,ny,layers,g_vec)
    else
        call evaluate_b_iso(b,h,u,v,nx,ny,layers,g_vec)
    endif
!
!    calculate relative vorticity!
    call evaluate_zeta(zeta,u,v,nx,ny,layers,dx,dy)

!     Calculate dhdt, dudt, dvdt at current time step
    call evaluate_dhdt(dhdtold, h,u,v,ah,dx,dy,nx,ny,layers, spongeHTimeScale,spongeH,wetmask)

    call evaluate_dudt(dudtold, h,u,v,b,zeta,wind_x,fu, au,ar,slip,dx,dy,hfacN,hfacS,nx,ny,layers,rho0,spongeUTimeScale,spongeU)

    !if the wind is changed to be meridional this code will need changing

    call evaluate_dvdt(dvdtold, h,u,v,b,zeta,wind_y,fv, au,ar,slip, &
        dx,dy,hfacW,hfacE,nx,ny,layers,rho0,spongeVTimeScale,spongeV)

!    Calculate the values at half the time interval with Forward Euler
    hhalf = h+0.5d0*dt*dhdtold
    uhalf = u+0.5d0*dt*dudtold
    vhalf = v+0.5d0*dt*dvdtold


!    calculate Bernoulli potential
    if (RedGrav) then
        call evaluate_b_RedGrav(b,hhalf,uhalf,vhalf,nx,ny,layers,g_vec)
    else
        call evaluate_b_iso(b,hhalf,uhalf,vhalf,nx,ny,layers,g_vec)
    endif

!    calculate relative vorticity
    call evaluate_zeta(zeta,uhalf,vhalf,nx,ny,layers,dx,dy)

!    Now calculate d/dt of u,v,h and store as dhdtold, dudtold and dvdtold
    call evaluate_dhdt(dhdtold, hhalf,uhalf,vhalf,ah,dx,dy,nx,ny,layers, spongeHTimeScale,spongeH,wetmask)
    call evaluate_dudt(dudtold, hhalf,uhalf,vhalf,b,zeta,wind_x,&
        fu, au,ar,slip,dx,dy,hfacN,hfacS,nx,ny,layers,rho0,&
        spongeUTimeScale,spongeU)
    !if the wind is changed to be meridional this code will need changing
    call evaluate_dvdt(dvdtold, hhalf,uhalf,vhalf,b,zeta,wind_y,&
        fv, au,ar,slip, dx,dy,hfacW,hfacE,nx,ny,layers,rho0,& 
        spongeVTimeScale,spongeV)
!    These are the values to be stored in the 'old' variables ready to start 
!    the proper model run.

!    Calculate h,u,v with these tendencies
    h = h + dt*dhdtold
    u = u + dt*dudtold
    v = v + dt*dvdtold

!    Apply the boundary conditions
!    Enforce no normal flow boundary condition
!    and no flow in dry cells.
!    no/free-slip is done inside dudt and dvdt subroutines.
!    hfacW and hfacS are zero where the transition between
!    wet and dry cells occurs. wetmask is 1 in wet cells,
!    and zero in dry cells.
    do k = 1,layers
        u(:,:,k) = u(:,:,k)*hfacW*wetmask(1:nx,:)
        v(:,:,k) = v(:,:,k)*hfacS*wetmask(:,1:ny)
    enddo


!    Now the model is ready to start.
!    We have h, u, v at the zeroth time step, and the tendencies at two older time steps.
!    The model then solves for the
!    tendencies at the current step before solving for the fields at the 
!    next time step.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!        MAIN LOOP OF THE MODEL STARTS HERE        !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do 9999 n=1,nTimeSteps
      
!   time varying winds
    if (UseSinusoidWind .eqv. .TRUE.) then

        wind_x = base_wind_x*(wind_alpha +  &
            wind_beta*sin(((2d0*pi*n*dt)/wind_period) - wind_t_offset))

        wind_y = base_wind_y*(wind_alpha +  &
            wind_beta*sin(((2d0*pi*n*dt)/wind_period) - wind_t_offset))

    else if (UseStochWind .eqv. .TRUE.) then
        
        if (mod(n-1,n_stoch_wind).eq.0) then 

!           gives a pseudorandom number in range 0 <= x < 1
            call random_number(stoch_wind_mag)
!           convert to -1 <= x < 1
            stoch_wind_mag = (stoch_wind_mag - 0.5d0)*2d0

            wind_x = base_wind_x*(wind_alpha +  &
                wind_beta*stoch_wind_mag)

            wind_y = base_wind_y*(wind_alpha +  &
                wind_beta*stoch_wind_mag)
        endif
    endif



!     calculate Bernoulli potential
    if (RedGrav) then
        call evaluate_b_RedGrav(b,h,u,v,nx,ny,layers,g_vec)
    else
        call evaluate_b_iso(b,h,u,v,nx,ny,layers,g_vec)
    endif

!    calculate relative vorticity
    call evaluate_zeta(zeta,u,v,nx,ny,layers,dx,dy)
!
!     Calculate dhdt, dudt, dvdt at current time step
    call evaluate_dhdt(dhdt, h,u,v,ah,dx,dy,nx,ny,layers, &
        spongeHTimeScale,spongeH,wetmask)

    call evaluate_dudt(dudt, h,u,v,b,zeta,wind_x,fu, au,ar,slip,&
        dx,dy,hfacN,hfacS,nx,ny,layers,rho0,&
        spongeUTimeScale,spongeU)

    call evaluate_dvdt(dvdt, h,u,v,b,zeta,wind_y,fv, au,ar,slip, &
        dx,dy,hfacW,hfacE,nx,ny,layers,rho0,spongeVTimeScale,spongeV)

!    Use dh/dt, du/dt and dv/dt to step h, u and v forward in time with
!    Adams-Bashforth third order linear multistep method


    unew=u+dt*(23d0*dudt - 16d0*dudtold + 5d0*dudtveryold)/12d0

    vnew=v+dt*(23d0*dvdt - 16d0*dvdtold + 5d0*dvdtveryold)/12d0

    hnew=h+dt*(23d0*dhdt - 16d0*dhdtold + 5d0*dhdtveryold)/12d0

!    Apply the boundary conditions
!    Enforce no normal flow boundary condition
!    and no flow in dry cells.
!    no/free-slip is done inside dudt and dvdt subroutines.
!    hfacW and hfacS are zero where the transition between
!    wet and dry cells occurs. wetmask is 1 in wet cells,
!    and zero in dry cells.
    do k = 1,layers
        unew(:,:,k) = unew(:,:,k)*hfacW*wetmask(1:nx,:)
        vnew(:,:,k) = vnew(:,:,k)*hfacS*wetmask(:,1:ny)
    enddo


!   Do the isopycnal layer physics
    if (.not. RedGrav) then
        ! calculate the barotropic velocities
        call calc_baro_u(ub,unew,hnew,eta,freesurfFac,nx,ny,layers)
        call calc_baro_v(vb,vnew,hnew,eta,freesurfFac,nx,ny,layers)

        ! calculate divergence of ub and vb, and solve for the pressure field that removes it
        call calc_eta_star(ub,vb,eta,etastar,freesurfFac,nx,ny,dx,dy,dt)
        !print *, maxval(abs(etastar))

        ! prevent barotropic signals from bouncing around outside the wet region of the model.
        etastar = etastar*wetmask

        call SOR_solver(a,etanew,etastar,freesurfFac,nx,ny, &
            dt,rjac,eps,maxits,n)
        !print *, maxval(abs(etanew))

        ! now update the velocities using the barotropic tendency due to the pressure gradient
        dudt_bt = 0d0
        dvdt_bt = 0d0

        do i = 2,nx-1
            do j = 1,ny-1
                dudt_bt(i,j) = -g_vec(1)*(etanew(i,j) - etanew(i-1,j))/(dx)
                unew(i,j,:) = unew(i,j,:) + 23d0*dt*dudt_bt(i,j)/12d0
            end do
        end do

        do i = 1,nx-1
            do j = 2,ny-1
                dvdt_bt(i,j) = -g_vec(1)*(etanew(i,j) - etanew(i,j-1))/(dy)
                vnew(i,j,:) = vnew(i,j,:) + 23d0*dt*dvdt_bt(i,j)/12d0
            end do
        end do

        ! We now have correct velocities at the next time level, but the layer thicknesses were updated with the old velocities. force consistency by scaling thicknesses to agree with free surface.

        h_norming = (freesurfFac*etanew+depth)/sum(hnew,3)
        do k = 1,layers
            hnew(:,:,k) = hnew(:,:,k)*h_norming
        end do

        if (maxval(abs(h_norming - 1d0)).gt. 1e-2) then
            print *, 'substantial inconsistency between h and eta (in %).', maxval(abs(h_norming - 1d0))*100
        end if

    !    Enforce no normal flow boundary condition
    !    and no flow in dry cells.
    !    no/free-slip is done inside dudt and dvdt subroutines.
    !    hfacW and hfacS are zero where the transition between
    !    wet and dry cells occurs. wetmask is 1 in wet cells,
    !    and zero in dry cells.
        do k = 1,layers
            unew(:,:,k) = unew(:,:,k)*hfacW*wetmask(1:nx,:)
            vnew(:,:,k) = vnew(:,:,k)*hfacS*wetmask(:,1:ny)
        enddo
    endif
    

    
! code to stop layers getting too thin
    counter = 0
    
    do k = 1,layers
      do j=1,ny-1
        do i=1,nx-1    
          if (hnew(i,j,k) .lt. hmin) then
        hnew(i,j,k) = hmin
        counter = counter + 1
        if (counter .eq. 1) then
          ! write a file saying that the layer thickness value dropped below hmin and this line has been used
          OPEN(UNIT=10, FILE='layer thickness dropped below hmin.txt', ACTION="write", STATUS="unknown", &
          FORM="formatted", POSITION = "append")
          write(10,1111) n
1111              format( "layer thickness dropped below hmin at time step ", 1i10.10)  
          CLOSE(UNIT=10)
        endif
          endif
        end do
      end do
    end do

!    Accumulate average fields
    if (avwrite .ne. 0) then
        hav = hav + hnew
        uav = uav + unew
        vav = vav + vnew
        if (.not. RedGrav) then
            etaav = eta + etanew
        endif
    end if

!     shuffle arrays: old -> very old,  present -> old, new -> present
!    height and velocity fields
    h = hnew
    u = unew
    v = vnew
    if (.not. RedGrav) then
        eta = etanew
    endif

!    tendencies (old -> very old)
    dhdtveryold = dhdtold
    dudtveryold = dudtold
    dvdtveryold = dvdtold

!    tendencies (current -> old)
    do k = 1,layers
        dudtold(:,:,k) = dudt(:,:,k) + dudt_bt
        dvdtold(:,:,k) = dvdt(:,:,k) + dvdt_bt
    end do
    dhdtold = dhdt

!     now have new fields in main arrays and old fields in old arrays

!--------------------------- Code for generating output -------------------------------------
!     write data to file? 
      if (mod(n-1,nwrite).eq.0) then 

        !     
        !----------------------- Snapshot Output -----------------
        write(num,'(i10.10)') n

        !output the data to a file

        open(unit = 10, status='replace',file='output/snap.h.'//num,form='unformatted')  
        write(10) h
        close(10) 
        open(unit = 10, status='replace',file='output/snap.u.'//num,form='unformatted')  
        write(10) u
        close(10) 
        open(unit = 10, status='replace',file='output/snap.v.'//num,form='unformatted')  
        write(10) v
        close(10) 
        if (.not. RedGrav) then
            open(unit = 10, status='replace',file='output/snap.eta.'//num,form='unformatted')  
            write(10) eta
            close(10) 
        endif

        if (DumpWind .eqv. .TRUE.) then
            open(unit = 10, status='replace',file='output/wind_x.'//num,form='unformatted')  
            write(10) wind_x
            close(10) 
            open(unit = 10, status='replace',file='output/wind_y.'//num,form='unformatted')  
            write(10) wind_y
            close(10)
        endif 

        ! Check if there are NaNs in the data
        call break_if_NaN(h,nx,ny,layers,n)
        !     call break_if_NaN(u,nx,ny,layers,n)
        !     call break_if_NaN(v,nx,ny,layers,n)     
    
    endif

        !----------------------- Average Output -----------------
      if (avwrite .eq. 0) then 
        go to 120
      else if (mod(n-1,avwrite).eq.0) then 

        hav=hav/real(avwrite)
        uav=uav/real(avwrite)
        vav=vav/real(avwrite)

        write(num,'(i10.10)') n

        !output the data to a file

        open(unit = 10, status='replace',file='output/av.h.'//num,form='unformatted')  
        write(10) hav
        close(10) 
        open(unit = 10, status='replace',file='output/av.u.'//num,form='unformatted')  
        write(10) uav
        close(10) 
        open(unit = 10, status='replace',file='output/av.v.'//num,form='unformatted')  
        write(10) vav
        close(10) 
        if (.not. RedGrav) then
            open(unit = 10, status='replace',file='output/av.eta.'//num,form='unformatted')  
            write(10) etaav
            close(10) 
        endif

        ! Check if there are NaNs in the data
        call break_if_NaN(h,nx,ny,layers,n)
        !     call break_if_NaN(u,nx,ny,layers,n)
        !     call break_if_NaN(v,nx,ny,layers,n)     

        !    Reset average quantities
        hav = 0.0
        uav = 0.0
        vav = 0.0
        if (.not. RedGrav) then
            etaav = 0.0
        endif
        !       h2av=0.0
    
120 endif

9999  continue

    OPEN(UNIT=10, FILE='run_finished.txt', ACTION="write", STATUS="unknown", &
    FORM="formatted", POSITION = "append")
    write(10,1112) n
    1112      format( "run finished at time step ", 1i10.10)  
    CLOSE(UNIT=10)

    stop

    END PROGRAM MIM



! ----------------------------- Beginning of the subroutines --------------------
!
!-----------------------------------------------------------------------------------
!> Evaluate the Bornoulli Potential for each grid box.
!! B is evaluated at the tracer point, centre of the grid, for each grid box.
!!
!! 
!!!
    subroutine evaluate_b_iso(b,h,u,v,nx,ny,layers,g_vec)
!    Evaluate baroclinic component of the Bernoulli Potential in the n-layer physics, at centre of grid box
    integer nx,ny,layers
    integer i,j,k
    double precision h(0:nx,0:ny,layers),u(nx,0:ny,layers),v(0:nx,ny,layers)
    double precision b(nx,ny,layers), g_vec(layers)
    double precision b_proto

    b = 0d0
!   no pressure contribution to the first layer Bernoulli potential 
!   (the barotropic pressure contributes, but that's not done here).
    do j = 1,ny-1
        do i = 1,nx-1
            b(i,j,1) = (u(i,j,1)**2+u(i+1,j,1)**2+v(i,j,1)**2+v(i,j+1,1)**2)/4.0d0
        end do
    end do

!   for the rest of the layers we get a baroclinic pressure contribution
    do k = 2,layers !move through the different layers of the model
      do j=1,ny-1 !move through longitude
        do i=1,nx-1 ! move through latitude
          ! The following loop is to get the baroclinic pressure term in the Bernoulli Potential
          b_proto = 0d0

          do l = 2,k
            b_proto = b_proto + g_vec(l)*h(i,j,l) !sum up the product of reduced gravity and layer thicknesses to form the baroclinic pressure componenet of the Bernoulli Potential term.
          end do
          b(i,j,k)= b_proto + (u(i,j,k)**2+u(i+1,j,k)**2+v(i,j,k)**2+v(i,j+1,k)**2)/4.0d0
          ! Add the (u^2 + v^2)/2 term to the pressure componenet of the Bernoulli Potential
        end do 
      end do 
    end do

    return
    end subroutine
!---------------------------------------------------------------------------------------

    subroutine evaluate_b_RedGrav(b,h,u,v,nx,ny,layers,gr)
!    Evaluate Bernoulli Potential at centre of grid box
    integer nx,ny,layers
    integer i,j,k
    double precision h(0:nx,0:ny,layers),u(nx,0:ny,layers),v(0:nx,ny,layers)
    double precision b(nx,ny,layers), gr(layers)
    double precision h_temp, b_proto

    b = 0d0

    do k = 1,layers !move through the different layers of the model
      do j=1,ny-1 !move through longitude
        do i=1,nx-1 ! move through latitude
          ! The following loops are to get the pressure term in the Bernoulli Potential
          b_proto = 0d0
          do l = k,layers
            h_temp = 0d0
            do m = 1,l
              h_temp = h_temp + h(i,j,m) !sum up the layer thicknesses
            end do
            b_proto = b_proto + gr(l)*h_temp !sum up the product of reduced gravity and summed layer thicknesses to form the pressure componenet of the Bernoulli Potential term
          end do
          b(i,j,k)= b_proto + (u(i,j,k)**2+u(i+1,j,k)**2+v(i,j,k)**2+v(i,j+1,k)**2)/4.0d0
          ! Add the (u^2 + v^2)/2 term to the pressure componenet of the Bernoulli Potential
        end do 
      end do 
    end do


    return
    end subroutine

!---------------------------------------------------------------------------------------
!>    Evaluate relative vorticity at lower left grid boundary (du/dy 
!!    and dv/dx are at lower left corner as well)
    subroutine evaluate_zeta(zeta,u,v,nx,ny,layers,dx,dy)
    integer nx,ny,layers
    integer i,j,k
    double precision h(0:nx,0:ny,layers),u(nx,0:ny,layers),v(0:nx,ny,layers)
    double precision zeta(nx,ny,layers)
    double precision dx, dy

    zeta = 0d0
    
    do k = 1,layers
      do j=1,ny
        do i=1,nx
          zeta(i,j,k)=(v(i,j,k)-v(i-1,j,k))/dx-(u(i,j,k)-u(i,j-1,k))/dy
        end do 
      end do 
    end do

    return
    end subroutine
!------------------------------------------------------------------------------------------
!> Calculate the tendency of layer thickness for each of the active layers
!! dh/dt is in the centre of each grid point.
    subroutine evaluate_dhdt(dhdt, h,u,v,ah,dx,dy,nx,ny,layers, spongeTimeScale,spongeH,wetmask)
!    dhdt is evaluated at the centre of the grid box
    integer nx,ny,layers
    integer i,j,k
    double precision dhdt(0:nx,0:ny,layers), h(0:nx,0:ny,layers)
    double precision u(nx,0:ny,layers),v(0:nx,ny,layers)
    double precision spongeTimeScale(0:nx,0:ny,layers)
    double precision spongeH(0:nx,0:ny,layers)
    double precision wetmask(0:nx,0:ny)
    double precision dx, dy
    double precision ah(layers)

    dhdt = 0d0

    do k = 1,layers
      do j=1,ny-1
        do i=1,nx-1
          dhdt(i,j,k)=ah(k)*(h(i+1,j,k)*wetmask(i+1,j) + &
                            (1d0 - wetmask(i+1,j))*h(i,j,k) & ! boundary
                        + h(i-1,j,k)*wetmask(i-1,j) + &
                            (1d0 - wetmask(i-1,j))*h(i,j,k) & ! boundary
                        - 2*h(i,j,k))/ &
        (dx*dx) + & ! x-componenet 
          ah(k)*(h(i,j+1,k)*wetmask(i,j+1) + &
                (1d0 - wetmask(i,j+1))*h(i,j,k) & ! reflect value around boundary
            + h(i,j-1,k)*wetmask(i,j-1) + &
                (1d0 - wetmask(i,j-1))*h(i,j,k) & ! reflect value around boundary
            - 2*h(i,j,k))/ &
        (dy*dy) - & !y-component horizontal diffusion
        ((h(i,j,k)+h(i+1,j,k))*u(i+1,j,k) - (h(i-1,j,k)+h(i,j,k))*u(i,j,k))/(dx*2d0) - & !d(hu)/dx
        ((h(i,j,k)+h(i,j+1,k))*v(i,j+1,k) - (h(i,j-1,k)+h(i,j,k))*v(i,j,k))/(dy*2d0) + & !d(hv)/dy
        spongeTimeScale(i,j,k)*(spongeH(i,j,k)-h(i,j,k)) ! forced relaxtion in the sponge regions.
        end do 
      end do 
    end do

    do k = 1,layers
        dhdt(:,:,k) = dhdt(:,:,k)*wetmask
    end do

    return
    end subroutine
!--------------------------------------------------------------------------------------------
!> Calculate the tendency of zonal velocity for each of the active layers

    subroutine evaluate_dudt(dudt, h,u,v,b,zeta,wind_x,fu, au,ar,slip,dx,dy,hfacN,hfacS,nx,ny,layers,rho0,spongeTimeScale,spongeU)
!    dudt(i,j) is evaluated at the centre of the left edge of the grid box,
!     the same place as u(i,j)
    integer nx,ny,layers
    integer i,j,k
    double precision dudt(nx,0:ny,layers), h(0:nx,0:ny,layers)
    double precision u(nx,0:ny,layers),v(0:nx,ny,layers)
    double precision fu(nx,0:ny), zeta(nx,ny,layers)
    double precision b(nx,ny,layers)
    double precision spongeTimeScale(nx,0:ny,layers)
    double precision spongeU(nx,0:ny,layers)
    double precision wind_x(nx, 0:ny)
    double precision dx, dy
    double precision au, ar, rho0, slip
    double precision hfacN(0:nx,ny)
    double precision hfacS(0:nx,ny)

    dudt = 0d0

    do k = 1,layers
      do i=2,nx-1  
        do j=1,ny-1
            dudt(i,j,k)= au*(u(i+1,j,k)+u(i-1,j,k)-2.0d0*u(i,j,k))/ & 
            (dx*dx)  + &   !    x-component
            au*(u(i,j+1,k)+u(i,j-1,k)-2.0d0*u(i,j,k)+ &
            ! boundary conditions
            (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacN(i,j))*u(i,j,k) + &
            (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacS(i,j))*u(i,j,k))/ &
            (dy*dy) + & !y-component
    !    together make the horizontal diffusion term
            0.25d0*(fu(i,j)+0.5d0*(zeta(i,j,k)+zeta(i,j+1,k)))* &
            (v(i-1,j,k)+v(i,j,k)+v(i-1,j+1,k)+v(i,j+1,k)) - & !vorticity term
            (b(i,j,k) - b(i-1,j,k))/dx + & !Bernoulli Potential term
            spongeTimeScale(i,j,k)*(spongeU(i,j,k)-u(i,j,k)) ! forced relaxtion in the sponge regions)
            if (k .eq. 1) then ! only have wind forcing on the top layer
                !This will need refining in the event of allowing outcropping.
                dudt(i,j,k) = dudt(i,j,k) + wind_x(i,j)/(rho0*h(i,j,k)) !wind forcing
            endif
            if (layers .gt. 1) then ! only evaluate vertical momentum diffusivity if more than 1 layer
                if (k .eq. 1) then ! adapt vertical momentum diffusivity for 2+ layer model -> top layer
                    dudt(i,j,k) = dudt(i,j,k) - 1.0d0*ar*(u(i,j,k) - 1.0d0*u(i,j,k+1))
                else if (k .eq. layers) then ! bottom layer
                    dudt(i,j,k) = dudt(i,j,k) - 1.0d0*ar*(u(i,j,k) - 1.0d0*u(i,j,k-1))
                else ! mid layer/s
                    dudt(i,j,k) = dudt(i,j,k) - 1.0d0*ar*(2.0d0*u(i,j,k) - 1.0d0*u(i,j,k-1) - 1.0d0*u(i,j,k+1))
                endif
            endif
        end do
      end do
    end do
    return
    end subroutine
!------------------------------------------------------------------------------------
!> Calculate the tendency of meridional velocity for each of the active layers

    subroutine evaluate_dvdt(dvdt, h,u,v,b,zeta,wind_y,fv, au,ar,slip,dx,dy,hfacW,hfacE,nx,ny,layers,rho0,spongeTimeScale,spongeV)
!    dvdt(i,j) is evaluated at the centre of the bottom edge of the grid box,
!     the same place as v(i,j)
    integer nx,ny,layers
    integer i,j,k
    double precision dvdt(0:nx,ny,layers), h(0:nx,0:ny,layers),u(nx,0:ny,layers),v(0:nx,ny,layers)
    double precision fv(0:nx,ny), zeta(nx,ny,layers), b(nx,ny,layers)
    double precision spongeTimeScale(0:nx,ny,layers)
    double precision spongeV(0:nx,ny,layers)
    double precision wind_y(0:nx, ny)
    double precision dx, dy
    double precision au, ar, slip
    double precision rho0
    double precision hfacW(nx,0:ny)
    double precision hfacE(nx,0:ny)

    dvdt = 0d0
    
    do k = 1,layers
      do j=2,ny-1
        do i=1,nx-1
            dvdt(i,j,k)=au*(v(i+1,j,k)+v(i-1,j,k)-2.0d0*v(i,j,k) + & 
            ! boundary conditions
            (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacW(i,j))*v(i,j,k) + &
            (1.0d0 - 2.0d0*slip)*(1.0d0 - hfacE(i,j))*v(i,j,k))/ &
            (dx*dx) + & !x-component
              au*(v(i,j+1,k) + v(i,j-1,k) - 2.0d0*v(i,j,k))/ &
            (dy*dy) - & 
    !    y-component. Together these make the horizontal diffusion term
            0.25d0*(fv(i,j)+0.5d0*(zeta(i,j,k)+zeta(i+1,j,k)))* &
            (u(i,j-1,k)+u(i,j,k)+u(i+1,j-1,k)+u(i+1,j,k)) - & !vorticity term
            (b(i,j,k)-b(i,j-1,k))/dy + & ! Bernoulli Potential term
            spongeTimeScale(i,j,k)*(spongeV(i,j,k)-v(i,j,k)) ! forced relaxtion to vsponge (in the sponge regions)
            if (k .eq. 1) then ! only have wind forcing on the top layer
                !This will need refining in the event of allowing outcropping.
                dvdt(i,j,k) = dvdt(i,j,k) + wind_y(i,j)/(rho0*h(i,j,k)) !wind forcing
            endif
            if (layers .gt. 1) then ! only evaluate vertical momentum diffusivity if more than 1 layer
                if (k .eq. 1) then ! adapt vertical momentum diffusivity for 2+ layer model -> top layer
                    dvdt(i,j,k) = dvdt(i,j,k) - 1.0d0*ar*(v(i,j,k) - 1.0d0*v(i,j,k+1))
                else if (k .eq. layers) then ! bottom layer
                    dvdt(i,j,k) = dvdt(i,j,k) - 1.0d0*ar*(v(i,j,k) - 1.0d0*v(i,j,k-1))
                else ! mid layer/s
                    dvdt(i,j,k) = dvdt(i,j,k) - 1.0d0*ar*(2.0d0*v(i,j,k) - 1.0d0*v(i,j,k-1) - 1.0d0*v(i,j,k+1))
                endif
            endif
        end do 
      end do 
    end do

    return
    end subroutine
!--------------------------------------------------------------------------------------
!> Calculate the barotropic u velocity
    subroutine calc_baro_u(ub,u,h,eta,freesurfFac,nx,ny,layers)

    integer nx,ny,layers
    integer i,j,k
    double precision h(0:nx,0:ny,layers)
    double precision h_temp(0:nx,0:ny,layers)
    double precision eta(0:nx,0:ny)
    double precision freesurfFac
    double precision u(nx,0:ny,layers),ub(nx,0:ny)
    double precision proto_ub

    ub = 0d0

    ! add free surface elevation to the upper layer
    h_temp(:,:,1) = h(:,:,1) + eta*freesurfFac
    
    do i = 2,nx-1
        do j = 1,ny-1
            proto_ub = 0d0
            do k = 1,layers
                proto_ub = proto_ub + u(i,j,k)*(h_temp(i,j,k)+h_temp(i-1,j,k))/2d0
            end do
            ub(i,j) = proto_ub
        end do
    end do

    return
    end subroutine


! --------------------------------------------------------------

!> Calculate the barotropic v velocity
    subroutine calc_baro_v(vb,v,h,eta,freesurfFac,nx,ny,layers)

    integer nx,ny,layers
    integer i,j,k
    double precision h(0:nx,0:ny,layers)
    double precision h_temp(0:nx,0:ny,layers)
    double precision eta(0:nx,0:ny)
    double precision freesurfFac
    double precision v(0:nx,ny,layers),vb(0:nx,ny)
    double precision proto_vb

    vb = 0d0

    ! add free surface elevation to the upper layer
    h_temp(:,:,1) = h(:,:,1) + eta*freesurfFac
    
    do i = 1,nx-1
        do j = 2,ny-1
            proto_vb = 0d0
            do k = 1,layers
                proto_vb = proto_vb + v(i,j,k)*(h_temp(i,j,k)+h_temp(i,j-1,k))/2d0
            end do
            vb(i,j) = proto_vb
        end do
    end do

    return
    end subroutine


! --------------------------------------------------------------

!> calculate the free surface anomaly using the velocities 
! timestepped with the tendencies excluding the free surface
! pressure gradient.
    subroutine calc_eta_star(ub,vb,eta,etastar,freesurfFac,nx,ny,dx,dy,dt)
    integer nx,ny,layers
    integer i,j,k
    double precision eta(0:nx,0:ny)
    double precision etastar(0:nx,0:ny)
    double precision ub(nx,0:ny), vb(0:nx,ny)
    double precision freesurfFac, dx,dy,dt

    etastar = 0d0

    do i = 1,nx-1
        do j = 1,ny-1
            etastar(i,j) = freesurfFac*eta(i,j) - &
                    dt*((ub(i+1,j) - ub(i,j))/dx + (vb(i,j+1) - vb(i,j))/dy)
        end do
    end do

    return
    end subroutine

! --------------------------------------------------------------

!> Use successive over relaxation algorithm to solve the backwards Euler timestepping for the free surface anomaly, or for the surface pressure required to keep the barotropic flow nondivergent.

subroutine SOR_solver(a,etanew,etastar,freesurfFac,nx,ny,dt,rjac,eps,maxits,n)
    double precision a(5,nx-1,ny-1)
    double precision etanew(0:nx,0:ny), etastar(0:nx,0:ny)
    double precision freesurfFac
    integer nx,ny, i,j, maxits, n
    double precision dt
    double precision rjac, eps
    double precision rhs(nx-1,ny-1), res(nx-1,ny-1)
    double precision norm, norm0, norm_old
    double precision relax_param
    character(30) Format

    rhs = -etastar(1:nx-1,1:ny-1)/dt**2
    ! first guess for etanew
    etanew = etastar


    relax_param=1.d0 ! successive over relaxation parameter

! calculate initial residual, so that we can stop the loop when the current residual = norm0*eps
    norm0 = 0.d0
    do i = 1,nx-1
        do j = 1,ny-1
          res(i,j) = a(1,i,j)*etanew(i+1,j) &
                   + a(2,i,j)*etanew(i,j+1) &
                   + a(3,i,j)*etanew(i-1,j) &
                   + a(4,i,j)*etanew(i,j-1) &
                   + a(5,i,j)*etanew(i,j)   &
                   - freesurfFac*etanew(i,j)/dt**2 & 
                   - rhs(i,j)
          norm0 = norm0 + abs(res(i,j))
          etanew(i,j) = etanew(i,j)-relax_param*res(i,j)/a(5,i,j)
        end do
    end do

    !norm_old = norm0


    do nit = 1,maxits
        !norm_old = norm
        norm=0.d0
        do i=1,nx-1
            do j=1,ny-1
              res(i,j) = a(1,i,j)*etanew(i+1,j) &
                       + a(2,i,j)*etanew(i,j+1) &
                       + a(3,i,j)*etanew(i-1,j) &
                       + a(4,i,j)*etanew(i,j-1) &
                       + a(5,i,j)*etanew(i,j)   &
                       - freesurfFac*etanew(i,j)/dt**2 & 
                       - rhs(i,j)
              norm = norm + abs(res(i,j))
              etanew(i,j)=etanew(i,j)-relax_param*res(i,j)/(a(5,i,j))
            end do
        end do
        if (nit.eq.1) then 
            relax_param=1.d0/(1.d0-0.5d0*rjac**2)
        else
            relax_param=1.d0/(1.d0-0.25d0*rjac**2*relax_param)
        endif
      
        if (nit.gt.1.and.norm.lt.eps*norm0) then
            ! print 10014,  eps,  nit
            ! 10014   format('Res less than eps of original. eps, iteration: ', F5.4, I7)
            return
        ! elseif (nit.gt.1.and.abs(norm - norm_old).lt.norm_old*eps) then
        !     print *, 'change in residual less than eps of previous residual. iteration:', eps, nit
        !     return
        ! This was not a good way of exiting the solver - it would occasionally leave after 2 or 3 iterations.
        endif
    end do

    print *,'warning: maximum iterations exceeded at time step ', n

    return
    end subroutine


! ------------------------------------------------
!> Export data to file
    subroutine write_data(data,nx,ny,filename)

!    A subroutine to output data as a 2D array for plotting and things

    implicit none

    integer nx,ny
    double precision data(0:nx,0:ny)
    character (len=*) filename
    integer i,j

!    to show which field is being output
!    print *,filename 

    OPEN(UNIT=13, FILE=filename, ACTION="write", STATUS="replace", &
        FORM="formatted")

    do j=1,ny-1
        write(13,1000) (data(i,j),i=1,nx-1)
1000         format( 400f32.16)  
    end do  
    CLOSE(UNIT=13)
    return
    end
!------------------------------------------------------------------------------------
!> Check to see if there are any NaNs in the data field and stop the calculation
!! if any are found.
    subroutine break_if_NaN(data,nx,ny,layers,n)

!     To stop the program if it detects at NaN in the variable being checked

    integer nx, ny,layers,n
    double precision data(0:nx, 0:ny,layers)


      do k=1,layers
      do j=1,ny-1
      do i=1,nx-1
        if (data(i,j,k).ne.data(i,j,k))  then
        ! write a file saying so
        OPEN(UNIT=10, FILE='NaN detected.txt', ACTION="write", STATUS="replace", &
        FORM="formatted")
        write(10,1000) n
1000         format( "NaN detected at time step ", 1i10.10)  
        CLOSE(UNIT=10)
        ! print it on the screen
        print *, 'NaN detected' 
        ! Stop the code
        STOP
        endif
      end do
      end do
      end do
    
    return
    end

!-----------------------------------------------------------------------------------
!> Define masks for boundary conditions in u and v
!! this finds locations where neighbouring grid boxes are not the same (i.e. one is land and one is ocean). 
!! in the output
!! 0 means barrier
!! 1 mean open

    subroutine calc_boundary_masks(wetmask,hfacW,hfacE,hfacS,hfacN,nx,ny)

    implicit none
    integer nx, ny
    double precision wetmask(0:nx,0:ny)
    double precision temp(0:nx,0:ny)
    double precision hfacW(nx,0:ny)
    double precision hfacE(nx,0:ny)
    double precision hfacN(0:nx,ny)
    double precision hfacS(0:nx,ny)
    integer i,j


    hfacW = 1d0
! and now for all  western cells
    hfacW(1,:) = 0d0
    
    temp = 0.0
    do j = 0,ny
      do i = 1,nx
        temp(i,j) = wetmask(i-1,j)- wetmask(i,j)
      enddo
    enddo
        
    do j = 0,ny
      do i = 1,nx
        if (temp(i,j) .ne. 0.0) then
          hfacW(i,j) = 0d0
        endif
      enddo
    enddo

    hfacE = 1
! and now for all  eastern cells
    hfacE(nx,:) = 0
    
    temp = 0.0
    do j = 0,ny
      do i = 1,nx-1
        temp(i,j) = wetmask(i,j)- wetmask(i+1,j)
      enddo
    enddo
        
    do j = 0,ny
      do i = 1,nx
        if (temp(i,j) .ne. 0.0) then
          hfacE(i,j) = 0d0
        endif
      enddo
    enddo


    hfacS = 1
!   all southern cells
    hfacS(:,1) = 0
    temp = 0.0
    do j = 1,ny
      do i = 0,nx
        temp(i,j) = wetmask(i,j-1)- wetmask(i,j)
      enddo
    enddo
    
    do j = 1,ny
      do i = 0,nx
        if (temp(i,j) .ne. 0.0) then
          hfacS(i,j) = 0d0
        endif
      enddo
    enddo

    hfacN = 1
!   all northern cells
    hfacN(:,ny) = 0
    temp = 0.0
    do j = 1,ny-1
      do i = 0,nx
        temp(i,j) = wetmask(i,j)- wetmask(i,j+1)
      enddo
    enddo
    
    do j = 1,ny
      do i = 0,nx
        if (temp(i,j) .ne. 0.0) then
          hfacN(i,j) = 0d0
        endif
      enddo
    enddo

    return
    end

!-----------------------------------------------------------------


    subroutine read_input_fileH(name,array,default,nx,ny,layers)

    implicit none
    character(30) name
    integer nx, ny, layers, k
    double precision array(0:nx,0:ny,layers), default(layers)

    if (name.ne.'') then
        open(unit = 10, form='unformatted', file=name)  
        read(10) array
        close(10) 
    else
      do k = 1,layers
        array(:,:,k) = default(k)
      enddo
    endif

    return
    end


    subroutine read_input_fileH_2D(name,array,default,nx,ny)

    implicit none
    character(30) name
    integer nx, ny
    double precision array(0:nx,0:ny), default

    if (name.ne.'') then
        open(unit = 10, form='unformatted', file=name)  
        read(10) array
        close(10) 
    else
        array = default
    endif

    return
    end



    subroutine read_input_fileU(name,array,default,nx,ny,layers)

    implicit none
    character(30) name
    integer nx, ny,layers
    double precision array(nx,0:ny,layers), default

    if (name.ne.'') then
        open(unit = 10, form='unformatted', file=name)  
        read(10) array
        close(10) 
    else
        array = default
    endif

    return
    end


    subroutine read_input_fileV(name,array,default,nx,ny,layers)

    implicit none
    character(30) name
    integer nx, ny, layers
    double precision array(0:nx,ny,layers), default

    if (name.ne.'') then
        open(unit = 10, form='unformatted', file=name)  
        read(10) array
        close(10) 
    else
        array = default
    endif

    return
    end


!-----------------------------------------------------------------
!> Robustly produce a different random seed for each run.
!! Code adapted from http://web.ph.surrey.ac.uk/fortweb/glossary/random_seed.html


    subroutine ranseed()

    implicit none

    ! ----- variables for portable seed setting -----
    INTEGER :: i_seed
    INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
    INTEGER, DIMENSION(1:8) :: dt_seed
    ! ----- end of variables for seed setting -----

    ! ----- Set up random seed portably -----
    CALL RANDOM_SEED(size=i_seed)
    ALLOCATE(a_seed(1:i_seed))
    CALL RANDOM_SEED(get=a_seed)
    CALL DATE_AND_TIME(values=dt_seed)
    a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
    call random_seed(put=a_seed)
    return
    end

