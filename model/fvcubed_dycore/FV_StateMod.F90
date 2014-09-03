#include "MAPL_Generic.h"

module FV_StateMod
!BOP
!
! !MODULE: FV_StateMod --- GEOS5/CAM Cubed-Sphere fvcore state variables/init/run/finalize
!
! !USES:
#if defined( MAPL_MODE )
   use ESMF_Mod            ! ESMF base class
   use MAPL_Mod            ! MAPL base class
#endif

!   use MAPL_ConstantsMod

   use constants_mod,      only: grav, pi, kappa, cp_air

   use fms_mod, only: fms_init, set_domain, nullify_domain, file_exist
   use mpp_domains_mod, only: mpp_get_data_domain, mpp_get_compute_domain, &
                              mpp_update_domains, DGRID_NE, mpp_get_boundary
   use mpp_parameter_mod, only: AGRID_PARAM=>AGRID, CORNER
   use fv_timing_mod,    only: timing_on, timing_off, timing_init, timing_prt
   use mpp_mod, only: mpp_pe, mpp_npes, mpp_get_current_pelist
   use fv_mp_mod, only: gid, masterproc, tile, mp_start, domain_decomp, domain, npes_x, npes_y

   use fv_grid_utils_mod,  only: grid_utils_init, grid_utils_end, inner_prod, &
                          mid_pt_sphere, cubed_to_latlon, &
                          vlon, vlat, es, ew, edge_vect_s,edge_vect_n,edge_vect_w,edge_vect_e, &
                          f0, fC, ptop
   use fv_grid_tools_mod,  only: init_grid, cosa, sina, rarea, area, area_c, globalarea, &
                          dx, dy, rdxa, rdya, rotate_winds, atob_s, grid_type, &
                          get_unit_vector
   use fv_arrays_mod, only: fv_atmos_type
   use fv_restart_mod, only: fv_restart, fv_restart_init
   use fv_control_mod, only: fv_init, uniform_ppm, p_ref, n_sponge, &
                          DDDMP,   &     ! Del-2 Divergence Damping Coeficient
                          DDDM4,   &     ! Del-4 Divergence Damping Coeficient
                          D_EXT,   &     ! External divergence damping 
                          HORD_MT, &     ! Horizontal advection scheme (Momentum transport)
                          HORD_VT, &     ! Horizontal advection scheme (Vorticity transport)
                          HORD_TM, &     ! Horizontal advection scheme (Heat transport)
                          HORD_TR, &     ! Horizontal advection scheme (Tracer transport)
                          KORD_MT, &     ! monotonicity order for vertical remapping (Momentum)
                          KORD_TM, &     ! monotonicity order for vertical remapping (Heat)
                          KORD_TR        ! monotonicity order for vertical remapping (Tracers)

   use test_cases_mod, only: test_case
   use fv_dynamics_mod, only: fv_dynamics
   use sw_core_mod, only: d2a2c_vect
   use external_ic_mod, only: cubed_a2d

   use fv_diagnostics_mod, only: prt_maxmin

implicit none
private

#include "mpif.h"

  integer,  parameter :: r8           = 8
  integer,  parameter :: r4           = 4

  real(r8), save :: elapsed_time = 0

! !PUBLIC DATA MEMBERS:

  logical :: COLDSTART = .false.
  logical :: SW_DYNAMICS = .false.
  logical :: PURE_ADVECTION = .false.
  integer :: NXQ = 0

  public NXQ
  public COLDSTART, SW_DYNAMICS, PURE_ADVECTION
  public FV_InitState, FV_Run, FV_Add_Incs, FV_Finalize
  public FV_To_State, State_To_FV
  public T_TRACERS, T_FVDYCORE_VARS, T_FVDYCORE_GRID, T_FVDYCORE_STATE
  public fv_getTopography
  public fv_getMassFluxes
  public fv_getPK
  public fv_getOmega
  public fv_getAgridWinds
  public fv_getVorticity
  public a2d3d
  public latlon_to_cubed_winds

  integer,  parameter :: ntiles_per_pe = 1
  type(fv_atmos_type), save :: FV_Atm(ntiles_per_pe)

  type T_TRACERS
       logical                                   :: is_r4
       real(r8), dimension(:,:,:  ), pointer     :: content
       real(r4), dimension(:,:,:  ), pointer     :: content_r4
  end type T_TRACERS

! T_FVDYCORE_VARS contains the prognostic variables for FVdycore
  type T_FVDYCORE_VARS
       real(r8), dimension(:,:,:  ), pointer     :: U      ! U winds (D-grid)
       real(r8), dimension(:,:,:  ), pointer     :: V      ! V winds (D-grid)
       real(r8), dimension(:,:,:  ), pointer     :: PT     ! scaled virtual pot. temp.
       real(r8), dimension(:,:,:  ), pointer     :: PE     ! Pressure at layer edges
       real(r8), dimension(:,:,:  ), pointer     :: PKZ    ! P^kappa mean
       type(T_TRACERS), dimension(:), pointer    :: tracer ! Tracers
  end type T_FVDYCORE_VARS

! T_FVDYCORE_GRID contains information about the horizontal and vertical
! discretization, unlike in ARIES where these data are split into HORZ_GRID
! and VERT_GRID.  The reason for this: currently all of this information is
! initialized in one call to FVCAM dynamics_init.

  type T_FVDYCORE_GRID
!
#if defined( MAPL_MODE )
    type (MAPL_MetaComp),   pointer :: FVgenstate
    type (ESMF_Grid)                :: GRID           ! The 'horizontal' grid (2D decomp only)
#endif

    integer                         :: NG               ! Ghosting
!
    integer                         :: IS               ! Start X-index (exclusive, unghosted)
    integer                         :: IE               ! End X-index (exclusive, unghosted)
    integer                         :: JS               ! Start Y-index (exclusive, unghosted)
    integer                         :: JE               ! End Y-index (exclusive, unghosted)
!
    integer                         :: ISD              ! Start X-index (exclusive, ghosted)
    integer                         :: IED              ! End X-index (exclusive, ghosted)
    integer                         :: JSD              ! Start Y-index (exclusive, ghosted)
    integer                         :: JED              ! End Y-index (exclusive, ghosted)
!
    integer                         :: NPX             ! Full X- dim
    integer                         :: NPY             ! Full Y- dim
    integer                         :: NPZ             ! Numer of levels
    integer                         :: NPZ_P1          ! NPZ+1 (?)

    integer                         :: NTILES          ! How many log-rectangular tiles does my grid Have
                                                       ! lat-lon      = 1
                                                       ! cubed-sphere = 6

    real(r8), allocatable           :: AREA(:,:)       ! local cell area
    real(r8)                        :: GLOBALAREA      ! global area

    integer                         :: KS              ! Number of true pressure levels (out of NPZ+1)
    real(r8)                        :: PTOP            ! pressure at top (ak(1))
    real(r8)                        :: PINT            ! initial pressure (ak(npz+1))
    real(r8), dimension(:), pointer :: AK              ! Sigma mapping
    real(r8), dimension(:), pointer :: BK              ! Sigma mapping
    integer                         :: N_SPONGE        ! Number of sponge layers at top-of-atmosphere
    real(r8)                        :: f_coriolis_angle
!
! Tracers
!
    integer                         :: NQ              ! Number of advected tracers
    integer                         :: NTOTQ           ! Total number of tracers (NQ <= NC)
  end type T_FVDYCORE_GRID

! Constants used by fvcore
  type T_FVDYCORE_CONSTANTS
    real(r8)                             :: pi
    real(r8)                             :: omega    ! angular velocity of earth's rotation
    real(r8)                             :: cp       ! heat capacity of air at constant pressure
    real(r8)                             :: ae       ! radius of the earth (m)
    real(r8)                             :: rair     ! Gas constant of the air
    real(r8)                             :: cappa    ! Cappa?
    real(r8)                             :: zvir     ! RWV/RAIR-1
  end type T_FVDYCORE_CONSTANTS

  integer, parameter :: NUM_FVDYCORE_ALARMS        = 3
  integer, parameter :: NUM_TIMES      = 8
  integer, parameter :: TIME_TO_RUN  = 1
  integer, parameter :: CHECK_MAXMIN = 2

  type T_FVDYCORE_STATE
!!!    private
    type (T_FVDYCORE_VARS)               :: VARS
    type (T_FVDYCORE_GRID )              :: GRID
    type (T_FVDYCORE_CONSTANTS)          :: CONSTANTS
#if defined( MAPL_MODE )
    type (ESMF_Clock), pointer           :: CLOCK
    type (ESMF_Alarm)                    :: ALARMS(NUM_FVDYCORE_ALARMS)
#endif
    integer(kind=8)                      :: RUN_TIMES(4,NUM_TIMES)
    logical                              :: DOTIME, DODYN
    real(r8)                             :: DT          ! Large time step
    real(r8)                             :: CHECK_DT    ! Time step to check maxmin
    logical                              :: CONSV       ! dycore conserves tot. en.
    integer                              :: NSPLIT
    integer                              :: NUM_CALLS
  end type T_FVDYCORE_STATE

  real(kind=4), pointer             :: phis(:,:)
  logical :: phis_initialized = .false.

!
! !DESCRIPTION:
!
!      This module provides variables which are specific to the Lin-Rood
!      dynamical core.  Most of them were previously SAVE variables in
!      different routines and were set with an "if (first)" statement.
!
!      \begin{tabular}{|l|l|} \hline \hline
!        lr\_init    &  Initialize the Lin-Rood variables  \\ \hline
!        lr\_clean   &  Deallocate all internal data structures \\ \hline
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   2007.07.17   Putman     Created from lat-lon core
!
!EOP
!-----------------------------------------------------------------------
   real(r8), parameter ::  D0_0                    =   0.0_r8
   real(r8), parameter ::  D0_5                    =   0.5_r8
   real(r8), parameter ::  D1_0                    =   1.0_r8
   real(r8), parameter ::  D2_0                    =   2.0_r8
   real(r8), parameter ::  D4_0                    =   4.0_r8
   real(r8), parameter ::  D180_0                  = 180.0_r8
   real(r8), parameter ::  ratmax                  =  0.81_r8

contains

 subroutine FV_InitState (STATE, CLOCK, INTERNAL, IMPORT, GC, LAYOUT_FILE)

  use fv_control_mod, only : npx,npy,npz, ntiles, ncnst, nwat
  use fv_control_mod, only : hord_mt, hord_vt, kord_mt, hord_tm, hord_dp, hord_ze, kord_tm, hord_tr, kord_tr
  use fv_control_mod, only : n_split, m_split, q_split, master
  use fv_control_mod, only : dddmp, dddm4, d_ext
  use fv_control_mod, only : k_top, m_riem, n_sponge, p_ref
  use fv_control_mod, only : uniform_ppm, remap_t,  z_tracer, fv_debug
  use fv_control_mod, only : external_ic, ncep_ic, res_latlon_dynamics, res_latlon_tracers, fv_land
  use fv_control_mod, only : consv_te, fv_sg_adj, tau, tau_h2o, rf_center
  use fv_control_mod, only : fv_init, fv_end
  use fv_control_mod, only : adiabatic, nf_omega, moist_phys
  use fv_control_mod, only : hydrostatic, phys_hydrostatic,  hybrid_z, quick_p_c, quick_p_d, m_grad_p, a2b_ord, c2l_ord
  use fv_control_mod, only : nt_prog, nt_phys

  type (T_FVDYCORE_STATE),pointer              :: STATE

  type (ESMF_Clock), target,     intent(INOUT) :: CLOCK
  type (ESMF_GridComp)         , intent(IN   ) :: GC
  type (ESMF_State)            , intent(INOUT) :: INTERNAL
  type (ESMF_State)            , intent(INOUT) :: IMPORT
  character(LEN=*)             , intent(IN   ) :: LAYOUT_FILE

! Local variables

! Pointers to geography info in the MAPL MetaComp

  real,                 pointer :: LATS (:,:)
  real,                 pointer :: LONS (:,:)

  type (ESMF_Time)     :: currentTime
  type (ESMF_TimeInterval)     :: Time2Run
  type (ESMF_TimeInterval)     :: CheckMaxMin
  type (ESMF_Config)           :: cf
  type (ESMF_VM)               :: VM
  type (T_FVDYCORE_GRID) , pointer :: GRID
  integer              :: rc
  integer              :: status
  integer              :: len
  real(r8) :: REAL_PACK(6)
  real(r8) :: DT

  integer :: nstep, nymd, nhms
  integer :: yr, mm, dd, h, m, s
  integer :: INT_PACK(6)

  integer   :: is ,ie , js ,je    !  Local dims
  integer   :: isc,iec, jsc,jec   !  Local dims
  integer   :: isd,ied, jsd,jed   !  Local dims
  integer   :: ks                 !  True # press. levs
  integer   :: k                  !  Vertical loop index
  integer   :: ntotq
  integer   :: ng

  integer   :: i,j,n

  type (ESMF_Time) :: fv_time
  integer :: days, seconds

  character(len=ESMF_MAXSTR)       :: IAm='FV:FV_InitState'

  real(r8), pointer                   :: AK(:), BK(:)
  real(r8), dimension(:,:,:), pointer :: U, V, PT, PE, PKZ
  type (MAPL_MetaComp),       pointer :: mapl

  real(r8), ALLOCATABLE :: UA(:,:,:)
  real(r8), ALLOCATABLE :: VA(:,:,:)
  real(r8), ALLOCATABLE :: UD(:,:,:)
  real(r8), ALLOCATABLE :: VD(:,:,:)

  integer, allocatable :: pelist(:) 
  integer :: commID       

  real(r8):: p2(2), p3(2), p4(2)
  real(r8):: e1(3), e2(3), ex(3), ey(3)
  real(r8):: utmp, vtmp


! BEGIN

! Retrieve the pointer to the state
! ---------------------------------

  call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
  VERIFY_(STATUS)

! Save the mapl state for FVperf_module
! -------------------------------------

  STATE%GRID%FVgenstate => MAPL

! READ LAYOUT FILE
!
! Get the layout and store directly in the GRID data structure
!
  GRID => STATE%GRID     ! For convenience

  cf = ESMF_ConfigCreate(rc)
  call ESMF_ConfigLoadFile( cf, LAYOUT_FILE, rc = rc )
!
  call ESMF_ConfigGetAttribute   ( cf, ntiles, label = 'ntiles:', default=6, rc = rc )

  call ESMF_ConfigGetAttribute   ( cf, npx, label = 'npx:', default=32, rc = rc )
  call ESMF_ConfigGetAttribute   ( cf, npy, label=  'npy:', default=32, rc = rc )
  call ESMF_ConfigGetAttribute   ( cf, npz, label = 'npz:', default=72, rc = rc )
  if (npz == 1) SW_DYNAMICS = .true.
 
  call ESMF_ConfigGetAttribute   ( cf, npes_x, label = 'npes_x:', default=0, rc = rc )
  call ESMF_ConfigGetAttribute   ( cf, npes_y, label=  'npes_y:', default=0, rc = rc )

! Get other scalars
! -----------------

  call ESMF_ConfigGetAttribute( cf, dt, label='dt:', default=1800.0_r8, rc = rc )

  call ESMF_ConfigGetAttribute( cf, hord_mt, label='hord_mt:', default= 9, rc = rc )
  call ESMF_ConfigGetAttribute( cf, hord_vt, label='hord_vt:', default= 9, rc = rc )
  call ESMF_ConfigGetAttribute( cf, hord_tm, label='hord_tm:', default= 9, rc = rc )
  call ESMF_ConfigGetAttribute( cf, hord_tr, label='hord_tr:', default=12, rc = rc )
  call ESMF_ConfigGetAttribute( cf, kord_mt, label='kord_mt:', default= 8, rc = rc )
  call ESMF_ConfigGetAttribute( cf, kord_tm, label='kord_tm:', default=-8, rc = rc )
  call ESMF_ConfigGetAttribute( cf, kord_tr, label='kord_tr:', default= 8, rc = rc )

  call ESMF_ConfigGetAttribute( cf,   dddmp, label='dddmp:', default= 0.005_r8, rc = rc )
  call ESMF_ConfigGetAttribute( cf,   dddm4, label='dddm4:', default= 0.05_r8 , rc = rc )
  call ESMF_ConfigGetAttribute( cf,   d_ext, label='d_ext:', default= 0.02_r8 , rc = rc )
  call ESMF_ConfigGetAttribute( cf, a2b_ord, label='a2b_ord:', default= 4, rc = rc )

  if (SW_DYNAMICS) then
     call ESMF_ConfigGetAttribute( cf, n_sponge, label='n_sponge:', default=0, rc = rc )
  else
     call ESMF_ConfigGetAttribute( cf, n_sponge, label='n_sponge:', default=1, rc = rc )
  endif

  call ESMF_ConfigGetAttribute( cf, n_split, label='nsplit:', default=0, rc = rc )
  call ESMF_ConfigGetAttribute( cf, ntotq,   label='ntotq:', default=1, rc = rc )   ! default: 1 tracer

  call ESMF_ConfigGetAttribute( cf, nwat     , label='nwat:', default=0     , rc = rc )
  call ESMF_ConfigGetAttribute( cf, grid_type, label='grid_type:', default=0     , rc = rc )
  call ESMF_ConfigGetAttribute( cf, consv_te , label='consv_te:', default=0.0_r8, rc = rc )
  call ESMF_ConfigGetAttribute( cf, remap_t , label='remap_t:', default=.true., rc = rc )
  call ESMF_ConfigGetAttribute( cf, m_grad_p , label='m_grad_p:', default=1     , rc = rc )

  call ESMF_ConfigGetAttribute( cf, test_case, label='case:', default= 13, rc = rc )
  if ( (SW_DYNAMICS) .and. (test_case>10) ) then
     call WRITE_PARALLEL(test_case ,format='("Shallow Water Test Case TOO Large : ",(   I3))')
     ASSERT_(test_case < 10)
  endif
!
  GRID%NG     = 3 ; ng = 3
  GRID%NPX    = NPX
  GRID%NPY    = NPY
  GRID%NPZ    = NPZ
  GRID%NPZ_P1 = NPZ+1
  GRID%NTILES = 6
  GRID%N_SPONGE  = N_SPONGE
  GRID%NTOTQ  = MAX(1,NTOTQ)
  GRID%NQ     = MAX(1,NTOTQ)
  ncnst       = MAX(1,NTOTQ)
!
! FV likes npx;npy in terms of cell vertices
!
  npx=npx+1 ; npy=npy+1

!
! All done with configuration
!
  call ESMF_ConfigDestroy( cf, rc = rc )

!
! Calculate N_SPLIT if it was specified as 0
  if ( N_SPLIT == 0 ) n_split= INIT_NSPLIT(dt,npx-1,npy-1)

    call MAPL_TimerOn(MAPL,"--FMS_INIT")

! Start up FMS/MPP
    call fms_init(MPI_COMM_WORLD)
! Start up FV                   
    call fv_init(FV_Atm, DT)
    call MAPL_TimerOff(MAPL,"--FMS_INIT")

  ASSERT_(DT > 0.0                           )

  call WRITE_PARALLEL("Dynamics PE Layout ")
  call WRITE_PARALLEL(NPES_X    ,format='("NPES_X  : ",(   I3))')
  call WRITE_PARALLEL(NPES_Y    ,format='("NPES_Y  : ",(   I3))')

  call WRITE_PARALLEL(HORD_MT,format='(  "HORD_MT: ",(I3))')
  call WRITE_PARALLEL(HORD_VT,format='(  "HORD_VT: ",(I3))')
  call WRITE_PARALLEL(HORD_TM,format='(  "HORD_TM: ",(I3))')
  call WRITE_PARALLEL(HORD_TR,format='(  "HORD_TR: ",(I3))')
  call WRITE_PARALLEL(KORD_MT,format='(  "KORD_MT: ",(I3))')
  call WRITE_PARALLEL(KORD_TM,format='(  "KORD_TM: ",(I3))')
  call WRITE_PARALLEL(KORD_TR,format='(  "KORD_TR: ",(I3))')

  call WRITE_PARALLEL(DDDMP  ,format='(  "  DDDMP: ",(f6.4))')
  call WRITE_PARALLEL(DDDM4  ,format='(  "  DDDM4: ",(f6.4))')
  call WRITE_PARALLEL(D_EXT  ,format='(  "  D_EXT: ",(f6.4))')

  STATE%DOTIME= .TRUE.
  STATE%CHECK_DT  = 21600.   ! Check max and min of arrays every 6 hours.
  STATE%DT        = DT
  STATE%NSPLIT    = N_SPLIT

  call WRITE_PARALLEL(' ')
  call WRITE_PARALLEL(STATE%DT,format='("Dynamics time step : ",(F10.4))')
  call WRITE_PARALLEL(' ')

! 
! Get the main GRIDXY grid from the application (no longer set in this module)
! 
  call ESMF_GridCompGet(gc, grid=GRID%GRID, vm=vm, rc=STATUS)

! Get size, grid, and coordinate specifications

  call WRITE_PARALLEL((/npx,npy,npz/)       , &
    format='("Resolution of dynamics restart     =",3I5)'  )

  ks = FV_Atm(1)%ks ! ALT: this was the value when we read "old" style FV_internal restart
                    !      if needed, we could compute, ks by count(BK==0.0)
                    !      then FV will try to run slightly more efficient code
                    !      So far, GEOS-5 has used ks = 0
  ASSERT_(ks <= NPZ+1)
  call WRITE_PARALLEL(ks                          , &
     format='("Number of true pressure levels =", I5)'   )

! Get pointers to internal state vars
  call MAPL_GetPointer(internal, ak, "AK",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, bk, "BK",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, u, "U",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, v, "V",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, pt, "PT",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, pe, "PE",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, pkz, "PKZ",rc=status)
  VERIFY_(STATUS)

!
! Define STATE%VARS pointing to cubed-sphere fvdycore
!
! STATE%VARS%U => U     
! STATE%VARS%V => V     
! STATE%VARS%PT => PT   
! STATE%VARS%PE => PE   
! STATE%VARS%PKZ => PKZ 

  call CREATE_VARS ( FV_Atm(1)%isc, FV_Atm(1)%iec, FV_Atm(1)%jsc, FV_Atm(1)%jec,     &
                     1, npz, npz+1,          &
                     U, V, PT, PE, PKZ,      &
                     STATE%VARS  )

  GRID%IS     = FV_Atm(1)%isc
  GRID%IE     = FV_Atm(1)%iec
  GRID%JS     = FV_Atm(1)%jsc
  GRID%JE     = FV_Atm(1)%jec
  GRID%ISD    = FV_Atm(1)%isd
  GRID%IED    = FV_Atm(1)%ied
  GRID%JSD    = FV_Atm(1)%jsd
  GRID%JED    = FV_Atm(1)%jed

! Local Copy of dimensions

  IS     = FV_Atm(1)%isc
  IE     = FV_Atm(1)%iec
  JS     = FV_Atm(1)%jsc
  JE     = FV_Atm(1)%jec
  ISC    = FV_Atm(1)%isc         
  IEC    = FV_Atm(1)%iec 
  JSC    = FV_Atm(1)%jsc 
  JEC    = FV_Atm(1)%jec
  ISD    = FV_Atm(1)%isd
  IED    = FV_Atm(1)%ied
  JSD    = FV_Atm(1)%jsd
  JED    = FV_Atm(1)%jed

  allocate( GRID%AREA(IS:IE,JS:JE) )
  GRID%AREA = area(IS:IE,JS:JE)
  GRID%GLOBALAREA = globalarea

   if (GRID%f_coriolis_angle == -999) then
     fC(:,:) = 0.0
     f0(:,:) = 0.0
   else
     do j=jsd,jed+1
        do i=isd,ied+1
           fC(i,j) = 2.*MAPL_OMEGA*( -COS(FV_Atm(1)%grid(i,j,1))*COS(FV_Atm(1)%grid(i,j,2))*SIN(GRID%f_coriolis_angle) + &
                                               SIN(FV_Atm(1)%grid(i,j,2))*COS(GRID%f_coriolis_angle) )
        enddo
     enddo
     do j=jsd,jed
        do i=isd,ied
           f0(i,j) = 2.*MAPL_OMEGA*( -COS(FV_Atm(1)%agrid(i,j,1))*COS(FV_Atm(1)%agrid(i,j,2))*SIN(GRID%f_coriolis_angle) + &
                                               SIN(FV_Atm(1)%agrid(i,j,2))*COS(GRID%f_coriolis_angle) )
        enddo
     enddo
   endif

! Check coordinate information from MAPL_MetaComp
!--------------------------------------------
    call MAPL_Get(MAPL,                &
       LATS          = LATS,           & ! These are in radians
       LONS          = LONS,           & ! These are in radians
       INTERNAL_ESMF_STATE=INTERNAL,   &
                             RC=STATUS )
    VERIFY_(STATUS)

  STATE%CLOCK => CLOCK
  call ESMF_TimeIntervalSet(Time2Run, &
                            S=nint(STATE%DT), rc=status)
  VERIFY_(status)

  STATE%ALARMS(TIME_TO_RUN) = ESMF_AlarmCreate("Time2Run", clock, &
                              ringInterval=Time2Run, &
                              Enabled=.TRUE., rc=status) ; VERIFY_(status)
  call ESMF_AlarmEnable(STATE%ALARMS(TIME_TO_RUN), rc=status); VERIFY_(status)
  call ESMF_AlarmRingerOn(STATE%ALARMS(TIME_TO_RUN), rc=status); VERIFY_(status)

  call WRITE_PARALLEL(' ')
  call WRITE_PARALLEL(STATE%DT, &
    format='("INITIALIZED ALARM: DYN_TIME_TO_RUN EVERY ",F9.1," secs.")')

  call ESMF_TimeIntervalSet(CheckMaxMin, S=nint(STATE%CHECK_DT), rc=status)
  STATE%ALARMS(CHECK_MAXMIN) = ESMF_AlarmCreate("CheckMaxMin", clock, &
                               RingInterval=CheckMaxMin, &
                               Enabled=.TRUE., rc=status); VERIFY_(status)

  call WRITE_PARALLEL(' ')
  call WRITE_PARALLEL(STATE%CHECK_DT, &
    format='("INITIALIZED ALARM: CHECK MAX AND MIN EVERY ",F9.1," secs.")')


!  Clear wall clock time clocks and global budgets

  STATE%RUN_TIMES = 0
  STATE%NUM_CALLS = 0

  call ESMF_ClockGet( CLOCK, currTime=fv_time, rc=STATUS )
  VERIFY_(STATUS)
  call ESMF_TimeGet( fv_time, dayOfYear=days, s=seconds, rc=STATUS )
  VERIFY_(STATUS)

  if (COLDSTART) then
     VERIFY_(STATUS)
     ALLOCATE( UA(isd:ied,jsd:jed,1:npz) )
     ALLOCATE( VA(isd:ied,jsd:jed,1:npz) )
     ALLOCATE( UD(isd:ied,jsd:jed+1,1:npz) )
     ALLOCATE( VD(isd:ied+1,jsd:jed,1:npz) )
     UA(isc:iec,jsc:jec,:) = STATE%VARS%U(isc:iec,jsc:jec,:)
     VA(isc:iec,jsc:jec,:) = STATE%VARS%V(isc:iec,jsc:jec,:)
     call a2d3d( STATE%grid, UA, VA, UD, VD )

#if defined(ALT_CASE5_WINDS)
! Compute u-wind
  do j=js,je
     do i=is,ie
        call mid_pt_sphere(FV_Atm(1)%grid(i,j,1:2), FV_Atm(1)%grid(i+1,j,1:2), p2)
        utmp = 20.0*COS(FV_Atm(1)%grid(i,j,2))
        vtmp = 0.0
        p3(1) = FV_Atm(1)%grid(i,j,  1)
        p3(2) = FV_Atm(1)%grid(i,j,  2)
        p4(1) = FV_Atm(1)%grid(i+1,j,1)
        p4(2) = FV_Atm(1)%grid(i+1,j,2)
        call get_unit_vector(p3, p4, e1)
        call get_latlon_vector(p2, ex, ey)  ! note: p2 shifted
        STATE%VARS%U(i,j,:) = utmp*inner_prod(e1,ex) + vtmp*inner_prod(e1,ey)
      enddo
  enddo
! Compute v-wind
  do j=js,je
     do i=is,ie
        call mid_pt_sphere(FV_Atm(1)%grid(i,j,1:2), FV_Atm(1)%grid(i,j+1,1:2), p2)
        utmp = 20.0*COS(p2(2))
        vtmp = 0.0
        p3(1) = FV_Atm(1)%grid(i,j,  1)
        p3(2) = FV_Atm(1)%grid(i,j,  2)
        p4(1) = FV_Atm(1)%grid(i,j+1,1)
        p4(2) = FV_Atm(1)%grid(i,j+1,2)
        call get_unit_vector(p3, p4, e2)
        call get_latlon_vector(p2, ex, ey)  ! note: p2 shifted
        STATE%VARS%V(i,j,:) = utmp*inner_prod(e2,ex) + vtmp*inner_prod(e2,ey)
      enddo
  enddo
#endif


     STATE%VARS%U(isc:iec,jsc:jec,:) = UD(isc:iec,jsc:jec,:)
     STATE%VARS%V(isc:iec,jsc:jec,:) = VD(isc:iec,jsc:jec,:)
     DEALLOCATE ( UA )
     DEALLOCATE ( VA )
     DEALLOCATE ( UD )
     DEALLOCATE ( VD )
  endif


     call MAPL_GetPointer ( import, phis, 'PHIS', RC=STATUS )
     VERIFY_(STATUS)
    !FV_Atm(1)%phis(isc:iec,jsc:jec) = phis
    !call mpp_update_domains(FV_Atm(1)%phis, domain)

     FV_Atm(1)%ak = ak
     FV_Atm(1)%bk = bk
     ptop = FV_Atm(1)%ak(1)
     call State_To_FV( STATE )

#if defined(DEBUG)
  if (gid==0) print*,''
  if (gid==0) print*,'-------------- Initialize 1 ----------------------------'
  if (gid==0) print*,''
  call prt_maxmin('PS', FV_Atm(1)%ps, isc, iec, jsc, jec, ng,   1, 1.d-2, gid==0)
  call prt_maxmin('U ', FV_Atm(1)%u , isc, iec, jsc, jec+1, ng, npz, 1.d00, gid==0)
  call prt_maxmin('V ', FV_Atm(1)%v , isc, iec+1, jsc, jec, ng, npz, 1.d00, gid==0)
  call prt_maxmin('TA', FV_Atm(1)%pt, isc, iec, jsc, jec, ng, npz, 1.d00, gid==0)
  call prt_maxmin('Q1', FV_Atm(1)%q(isd,jsd,1,1), isc, iec, jsc, jec, ng, npz, 1.d0, gid==0)
  call prt_maxmin('PZ', FV_Atm(1)%pkz,isc, iec, jsc, jec, 0, npz, 1.d00, gid==0)
  if (gid==0) print*,'-------------- Initialize 1 ----------------------------'
#endif

  if ( (npz==1) .and. (COLDSTART) ) then
     call fv_restart(domain, FV_Atm, state%dt, seconds, days, .true., grid_type)
     call FV_To_State ( STATE )
    !phis = FV_Atm(1)%phis(isc:iec,jsc:jec)
  else
     call fv_restart(domain, FV_Atm, state%dt, seconds, days, .false., grid_type)
  endif

#if defined(DEBUG)
  if (gid==0) print*,''
  if (gid==0) print*,'-------------- Initialize 2 ----------------------------'
  if (gid==0) print*,''
  call prt_maxmin('PS', FV_Atm(1)%ps, isc, iec, jsc, jec, ng,   1, 1.d-2, gid==0)
  call prt_maxmin('U ', FV_Atm(1)%u , isc, iec, jsc, jec+1, ng, npz, 1.d00, gid==0)
  call prt_maxmin('V ', FV_Atm(1)%v , isc, iec+1, jsc, jec, ng, npz, 1.d00, gid==0)
  call prt_maxmin('TA', FV_Atm(1)%pt, isc, iec, jsc, jec, ng, npz, 1.d00, gid==0)
  call prt_maxmin('Q1', FV_Atm(1)%q(isd,jsd,1,1), isc, iec, jsc, jec, ng, npz, 1.d0, gid==0)
  call prt_maxmin('PZ', FV_Atm(1)%pkz,isc, iec, jsc, jec, 0, npz, 1.d00, gid==0)
  if (gid==0) print*,'-------------- Initialize 2 ----------------------------'
#endif

!
! Write the vertical coordinate to STDOUT
!
  if( gid.eq.0 ) then
        print *
        write(6,100)
100     format(2x,' k ','      A(k)    ',2x,' B(k)   ',2x,'  Pref    ',2x,'  DelP',/, &
               1x,'----',3x,'----------',2x,'--------',2x,'----------',2x,'---------' )
           k=0
        write(6,101) k+1,ak(k)*0.01, bk(k), ak(k)*0.01 + 1000.0*bk(k)
        do k=1,ubound(ak,1)
        write(6,102) k+1,ak(k)*0.01, bk(k), ak(k)*0.01 + 1000.0*bk(k), &
                     (ak(k)-ak(k-1))*0.01 + 1000.0*(bk(k)-bk(k-1))
        enddo

        print *
101     format(2x,i3,2x,f10.6,2x,f8.4,2x,f10.4)
102     format(2x,i3,2x,f10.6,2x,f8.4,2x,f10.4,3x,f8.4)
  endif

  return

contains

!-----------------------------------------------------------------------
! BOP
! !IROUTINE:  init_nsplit --- find proper value for nsplit if not specified
!
! !INTERFACE:
  integer function INIT_NSPLIT(dtime,npx,npy)
!
! !USES:
    use fv_control_mod, only : a2b_ord

    implicit none

! !INPUT PARAMETERS:
    real (r8), intent(in) :: dtime      !  time step
    integer, intent(in)   :: npx,npy    !  Global horizontal resolution

! !DESCRIPTION:
! 
!    If nsplit=0 (module variable) then determine a good value
!    for ns (used in fvdycore) based on resolution and the large-time-step
!    (dtime). The user may have to set this manually if instability occurs.
! 
! !REVISION HISTORY:
!   00.10.19   Lin     Creation
!   01.03.26   Sawyer  ProTeX documentation
!   01.06.10   Sawyer  Modified for dynamics_init framework
!   03.12.04   Sawyer  Moved here from dynamics_vars.  Now a function
!   07.16.07   Putman  Modified for cubed-sphere
!
! EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
    real (r8)   dimx
    real (r8)   dim0                      ! base dimension
    real (r8)   dt0                       ! base time step              
    real (r8)   ns0                       ! base nsplit for base dimension
    integer     ns                        ! final value to be returned
                     
    parameter ( dim0 = 180.  )
    parameter ( dt0  = 1800. )
    
    ns0  = 5.
    dimx = 4.0*npx
    if ( npx >= 120 ) ns0 = 6
    ns = nint ( ns0*abs(dtime)*dimx/(dt0*dim0) + 0.49 )
    ns = max ( 1, ns )

    init_nsplit = ns

    call WRITE_PARALLEL ( ns ,format='("Dynamics NSPLIT: ",(i2))' )

    return
  end function INIT_NSPLIT
!---------------------------------------------------------------------

end subroutine FV_InitState

subroutine FV_Run (STATE, CLOCK)

  type (T_FVDYCORE_STATE),pointer              :: STATE

  type (ESMF_Clock), target,     intent(IN   ) :: CLOCK

! Local variables
  integer              :: rc
  integer              :: status
  character(len=ESMF_MAXSTR)       :: IAm='FV:FV_Run'

  type (ESMF_Time) :: fv_time
  integer  :: days, seconds
  real(r8) :: time_total
  real(r8) :: zvir

  integer :: i,j,n,k
  integer :: isc,iec,jsc,jec,ngc
  integer :: isd,ied,jsd,jed
  integer :: npz

! Splitting for Pure Advection
  integer  :: myLOOPS, myNSPLIT
  real(r8) :: myDT

! Begin

  call ESMF_ClockGet( CLOCK, currTime=fv_time, rc=STATUS ) 
  VERIFY_(STATUS)
  call ESMF_TimeGet( fv_time, dayOfYear=days, s=seconds, rc=STATUS )
  VERIFY_(STATUS)

  time_total = days*86400. + seconds

  isc = FV_Atm(1)%isc
  iec = FV_Atm(1)%iec
  jsc = FV_Atm(1)%jsc
  jec = FV_Atm(1)%jec
  isd = FV_Atm(1)%isd
  ied = FV_Atm(1)%ied
  jsd = FV_Atm(1)%jsd
  jed = FV_Atm(1)%jed
  npz = FV_Atm(1)%npz
  ngc = FV_Atm(1)%ng

  zvir=1.0

  ! Be sure we have the correct PHIS
   if (.not. phis_initialized) then
    !call fv_getTopography(phis)
    !print*, phis(:,jsc)
     FV_Atm(1)%phis(isc:iec,jsc:jec) = real(phis,kind=r8)
     call mpp_update_domains(FV_Atm(1)%phis, domain)
     phis_initialized = .true.
   endif

  ! Pull Tracers
   if (STATE%GRID%NQ > 0) then
    FV_Atm(1)%ncnst = STATE%GRID%NQ
    deallocate( FV_Atm(1)%q )
    allocate  ( FV_Atm(1)%q(isd:ied  ,jsd:jed  ,npz, FV_Atm(1)%ncnst) )
    do n=1,FV_Atm(1)%ncnst
       if (state%vars%tracer(n)%is_r4) then
          FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,n) = state%vars%tracer(n)%content_r4(:,:,:)
       else
          FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,n) = state%vars%tracer(n)%content(:,:,:)
       endif
    enddo
   endif

   if ( (PURE_ADVECTION) .and. (FV_Atm(1)%n_split > 0) ) then
      myDT = state%dt/FV_Atm(1)%n_split
      myNSPLIT = 1
      myLOOPS = FV_Atm(1)%n_split
   else
      myDT = state%dt
      myNSPLIT = FV_Atm(1)%n_split
      myLOOPS = 1
   endif 

   do n=1,myLOOPS

    elapsed_time = elapsed_time + myDT

! Update FV with Internal State
    call State_To_FV( STATE )

    call set_domain(FV_Atm(1)%domain)  ! needed for diagnostic output done in fv_dynamics
    call fv_dynamics(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,   &
                     myDT, FV_Atm(1)%consv_te, FV_Atm(1)%fill, FV_Atm(1)%reproduce_sum, kappa,   &
                     cp_air, zvir, FV_Atm(1)%ks, FV_Atm(1)%ncnst, myNSPLIT, FV_Atm(1)%q_split, &
                     FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%w, FV_Atm(1)%delz,                             &
                     FV_Atm(1)%hydrostatic, FV_Atm(1)%pt, FV_Atm(1)%delp, FV_Atm(1)%q, FV_Atm(1)%ps,       &
                     FV_Atm(1)%pe, FV_Atm(1)%pk, FV_Atm(1)%peln, FV_Atm(1)%pkz,                         &
                     FV_Atm(1)%phis, FV_Atm(1)%omga, FV_Atm(1)%ua, FV_Atm(1)%va, FV_Atm(1)%uc, FV_Atm(1)%vc,  &
                     FV_Atm(1)%ak, FV_Atm(1)%bk, FV_Atm(1)%mfx, FV_Atm(1)%mfy, FV_Atm(1)%cx, FV_Atm(1)%cy,    &
                     FV_Atm(1)%ze0, FV_Atm(1)%hybrid_z, time_total, elapsed_time, PURE_ADVECTION)
    call nullify_domain()

! Copy FV to internal State
    if (.not. PURE_ADVECTION) call FV_To_State ( STATE )
   enddo

  ! Push Tracers
   if (STATE%GRID%NQ > 0) then
    do n=1,FV_Atm(1)%ncnst
       if (state%vars%tracer(n)%is_r4) then
          state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,n)
       else
          state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,n)
       endif
    enddo
   endif

#if defined(DEBUG) || defined(DEBUG1)
    if (mod(time_total,1.0)==0) then
    if (gid==0) print*,''
    if (gid==0) print*,'-------------- After Dynamics Execution ----------------------------'
    if (gid==0) print*,''
    call prt_maxmin('PS', FV_Atm(1)%ps, isc, iec, jsc, jec, ngc,   1, 1.d-2, gid==0)
!   call prt_maxmin('U ', FV_Atm(1)%u , isc, iec, jsc, jec+1, ngc, npz, 1.d00, gid==0)
!   call prt_maxmin('V ', FV_Atm(1)%v , isc, iec+1, jsc, jec, ngc, npz, 1.d00, gid==0)
!   call prt_maxmin('UA', FV_Atm(1)%ua, isc, iec, jsc, jec, ngc, npz, 1.d00, gid==0)
!   call prt_maxmin('VA', FV_Atm(1)%va, isc, iec, jsc, jec, ngc, npz, 1.d00, gid==0)
!   call prt_maxmin('TA', FV_Atm(1)%pt, isc, iec, jsc, jec, ngc, npz, 1.d00, gid==0)
!   do n=1,FV_Atm(1)%ncnst
!      call prt_maxmin('Q ', FV_Atm(1)%q(isc-ngc,jsc-ngc,1,n), isc, iec, jsc, jec, ngc, npz, 1.d0, gid==0)
!   enddo
    if (gid==0) print*,''
    if (gid==0) print*,'--------------------------------------------------------------------'
    if (gid==0) print*,''
    endif
#endif

end subroutine FV_Run

subroutine FV_Add_Incs ( STATE,IMPORT,DT,REPLACE, RC )

!
! !INPUT PARAMETERS:

   type(T_FVDYCORE_STATE),      pointer   :: STATE
   type(ESMF_State),       intent(INOUT)  :: IMPORT   !  A-grid XY
   real(r8),               intent(IN   )  :: DT
   logical, optional,      intent(IN   )  :: REPLACE
   integer, optional,      intent(  OUT)  :: RC

!
! !DESCRIPTION:  This routine adds the tendencies to the state,
!                weighted appropriately by the time step.  Temperature
!                tendencies are pressure weighted (ie., DELP*DT/Dt).
!                All tendencies are on the A-grid, and have an XY decomposition.
!

    integer               :: I1, IN, J1, JN, K, im, jm, km
    integer               :: KL, KU
    integer               :: I,J
    real(r8)              :: FAC, DTFAC, SUMOUT
    real(r8), allocatable ::    dum(:,:,:)
    real(r8), allocatable ::    pke(:,:,:),  dpinv(:,:,:)
    real(r8), allocatable :: tend_ua(:,:,:), tend_va(:,:,:)
    real(r8), allocatable :: tend_ud(:,:,:), tend_vd(:,:,:)
    real(kind=4), pointer :: tend(:,:,:)

    character(len=ESMF_MAXSTR)         :: IAm="FV_ADD_INCS"
    integer                            :: status
    type (T_FVDYCORE_GRID),    pointer :: GRID

    grid  => state%grid   ! direct handle to grid

    IF ( present(REPLACE) ) THEN
        FAC = 0.0     ! Replace the state with the tendencies
      DTFAC = 1.0
    ELSE
        FAC = 1.0D0   ! Add the tendencies to the state
      DTFAC = DT
    ENDIF

    i1 = state%grid%is
    in = state%grid%ie
    j1 = state%grid%js
    jn = state%grid%je
    im = state%grid%npx
    jm = state%grid%npy
    km = state%grid%npz

! **********************************************************************
! ****           Compute Pressure Thickness Inverse                 ****
! **********************************************************************

    ALLOCATE( dpinv(i1:in,j1:jn,km) )
    do k=1,km
       do j=j1,jn
          do i=i1,in
             if (( STATE%VARS%PE(i,j,k+1)-STATE%VARS%PE(i,j,k) ) /= 0.0) &
                 dpinv(i,j,k) = 1.0/( STATE%VARS%PE(i,j,k+1)-STATE%VARS%PE(i,j,k) )
          enddo
       enddo
    enddo

! *********************************************************************************
! ****                      Wind Tendencies                                    ****
! ****         Note: State Variables are on the D-Grid cubed-sphere oriented,  ****
! ****        while IMPORT Tendencies are on the A-Grid latlon oriented        ****
! *********************************************************************************

    ALLOCATE( tend_ua(state%grid%isd:state%grid%ied  ,state%grid%jsd:state%grid%jed  ,km) )
    ALLOCATE( tend_va(state%grid%isd:state%grid%ied  ,state%grid%jsd:state%grid%jed  ,km) )
    ALLOCATE( tend_ud(state%grid%isd:state%grid%ied  ,state%grid%jsd:state%grid%jed+1,km) )
    ALLOCATE( tend_vd(state%grid%isd:state%grid%ied+1,state%grid%jsd:state%grid%jed  ,km) )

    call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DUDT',RC=STATUS )
    VERIFY_(STATUS)
    tend_ua(i1:in,j1:jn,1:km) = tend

    call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DVDT',RC=STATUS )
    VERIFY_(STATUS)
    tend_va(i1:in,j1:jn,1:km) = tend

! Put the wind tendencies on the D-grid
! -------------------------------------
    call a2d3d( grid, tend_ua, tend_va, tend_ud, tend_vd )

#if defined(DEBUG)
    call prt_maxmin('TEND_UA', tend_ua, i1, in, j1, jn, FV_Atm(1)%ng, km, 1.d00, gid==0)
    call prt_maxmin('TEND_VA', tend_va, i1, in, j1, jn, FV_Atm(1)%ng, km, 1.d00, gid==0)
#endif

! Add the wind tendencies to the control variables
! ------------------------------------------------
    STATE%VARS%U(:,:,1:km) = FAC*STATE%VARS%U(:,:,1:km) + DTFAC*TEND_UD(i1:in,j1:jn,1:km)
    STATE%VARS%V(:,:,1:km) = FAC*STATE%VARS%V(:,:,1:km) + DTFAC*TEND_VD(i1:in,j1:jn,1:km)

    DEALLOCATE( tend_ua )
    DEALLOCATE( tend_va )
    DEALLOCATE( tend_ud )
    DEALLOCATE( tend_vd )

! **********************************************************************
! ****                     Pressure Tendency                        ****
! **********************************************************************

    call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DPEDT',RC=STATUS )
    VERIFY_(STATUS)

    KL = lbound( tend,3 )
    KU = ubound( tend,3 )
    ALLOCATE( dum(i1:in,j1:jn,KL:KU) )

    DUM = DTFAC*TEND

    STATE%VARS%PE = FAC*STATE%VARS%PE + DUM
    DEALLOCATE (DUM)

! **********************************************************************
! ****            Update Diagnostic State Variable PKZ              ****
! **********************************************************************

    ALLOCATE( pke(i1:in,j1:jn,1:km+1) )

    pke = STATE%VARS%PE**kappa
    do k=1,km
       do j=j1,jn
          do i=i1,in
             if (( STATE%VARS%PE(i,j,k+1)-STATE%VARS%PE(i,j,k) ) /= 0.0) &
                 STATE%VARS%PKZ(i,j,k) = ( pke(i,j,k+1)-pke(i,j,k) ) &
                          / ( kappa*log( STATE%VARS%PE(i,j,k+1)/STATE%VARS%PE(i,j,k) ) )
          enddo
       enddo
    enddo
    DEALLOCATE( PKE   )

! *********************************************************************
! ****                  Dry Temperature Tendency                   ****
! ****                  ------------------------                   ****
! ****  Note: State Variable is Potential Temperature T/P**kappa   ****
! ****             while IMPORT Coupling is (Delta_P)*DTDt         ****
! *********************************************************************

    call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DTDT',RC=STATUS )
    VERIFY_(STATUS)

    KL = lbound( tend,3 )
    KU = ubound( tend,3 )
    ALLOCATE( dum(i1:in,j1:jn,KL:KU) )

    DUM = DTFAC*TEND*DPINV/STATE%VARS%PKZ

    STATE%VARS%PT = FAC*STATE%VARS%PT + DUM
    DEALLOCATE (DUM)

    DEALLOCATE( DPINV )

   return

 end subroutine FV_Add_Incs

 subroutine FV_Finalize (STATE)

  type (T_FVDYCORE_STATE),pointer              :: STATE

  integer isc, iec, jsc, jec
  integer isd, ied, jsd, jed
  integer npz, ng

  isc = FV_Atm(1)%isc
  iec = FV_Atm(1)%iec
  jsc = FV_Atm(1)%jsc
  jec = FV_Atm(1)%jec
  isd = FV_Atm(1)%isd
  ied = FV_Atm(1)%ied
  jsd = FV_Atm(1)%jsd
  jed = FV_Atm(1)%jed
  npz = FV_Atm(1)%npz
  ng  = FV_Atm(1)%ng

#if defined(DEBUG)
  call State_To_FV( STATE )
  if (gid==0) print*,''
  if (gid==0) print*,'-------------- Finalize ----------------------------'
  if (gid==0) print*,''
  call prt_maxmin('PS', FV_Atm(1)%ps, isc, iec, jsc, jec, ng,   1, 1.d-2, gid==0)
  call prt_maxmin('U ', FV_Atm(1)%u , isc, iec, jsc, jec+1, ng, npz, 1.d00, gid==0)
  call prt_maxmin('V ', FV_Atm(1)%v , isc, iec+1, jsc, jec, ng, npz, 1.d00, gid==0)
  call prt_maxmin('TA', FV_Atm(1)%pt, isc, iec, jsc, jec, ng, npz, 1.d00, gid==0)
  call prt_maxmin('Q1', FV_Atm(1)%q(isd,jsd,1,1), isc, iec, jsc, jec, ng, npz, 1.d0, gid==0)
  if (gid==0) print*,'-------------- Finalize ----------------------------'
#endif

      call timing_off('TOTAL')
      call timing_prt( mpp_pe() )

      deallocate (    FV_Atm(1)%u )
      deallocate (    FV_Atm(1)%v )
      deallocate (   FV_Atm(1)%pt )
      deallocate ( FV_Atm(1)%delp )
      deallocate (    FV_Atm(1)%q )
      deallocate (   FV_Atm(1)%ps )
      deallocate (   FV_Atm(1)%pe )
      deallocate (   FV_Atm(1)%pk )
      deallocate ( FV_Atm(1)%peln )
      deallocate (  FV_Atm(1)%pkz )
      deallocate ( FV_Atm(1)%phis )
      deallocate ( FV_Atm(1)%omga )
      deallocate (   FV_Atm(1)%ua )
      deallocate (   FV_Atm(1)%va )
      deallocate (   FV_Atm(1)%uc )
      deallocate (   FV_Atm(1)%vc )
      deallocate ( FV_Atm(1)%mfx )
      deallocate ( FV_Atm(1)%mfy )
      deallocate (  FV_Atm(1)%cx )
      deallocate (  FV_Atm(1)%cy )
      deallocate (  FV_Atm(1)%ak )
      deallocate (  FV_Atm(1)%bk )

      deallocate ( FV_Atm(1)%u_srf )
      deallocate ( FV_Atm(1)%v_srf )
#ifdef FV_LAND
      deallocate ( FV_Atm(1)%sgh )
#endif
      deallocate ( FV_Atm(1)%oro )

! Non-hydrostatic:
      deallocate ( FV_Atm(1)%w )
      deallocate ( FV_Atm(1)%delz  )
      deallocate ( FV_Atm(1)%ze0   )

#if defined( MAPL_MODE )
!    call ESMF_GridDestroy  (STATE%GRID%GRID)
#endif

 end subroutine FV_Finalize

subroutine State_To_FV ( STATE )

!
! !INPUT PARAMETERS:

   type(T_FVDYCORE_STATE),      pointer   :: STATE

    integer               :: ISC,IEC, JSC,JEC, IM,JM,KM, NG
    integer               :: I,J,K,N
    real(r8)              :: akap

    real(r8) :: wbuffer(state%grid%npy+2,state%grid%npz)
    real(r8) :: sbuffer(state%grid%npx+2,state%grid%npz)

    character(len=ESMF_MAXSTR)         :: IAm="State_To_FV"

    ISC = state%grid%is
    IEC = state%grid%ie
    JSC = state%grid%js
    JEC = state%grid%je
    KM  = state%grid%npz
    NG  = state%grid%ng

    akap  = kappa
    if (SW_DYNAMICS) akap  = 1.

#if defined(DEBUG)
  if (gid==0) print*,''
  if (gid==0) print*,'-------------- State To FV 1 ----------------------------'
  if (gid==0) print*,''
  call prt_maxmin('PS', FV_Atm(1)%ps, isc, iec, jsc, jec, ng,   1, 1.d-2, gid==0)
  call prt_maxmin('U ', FV_Atm(1)%u , isc, iec, jsc, jec+1, ng, km, 1.d00, gid==0)
  call prt_maxmin('V ', FV_Atm(1)%v , isc, iec+1, jsc, jec, ng, km, 1.d00, gid==0)
  call prt_maxmin('TA', FV_Atm(1)%pt, isc, iec, jsc, jec, ng, km, 1.d00, gid==0)
  call prt_maxmin('PZ', FV_Atm(1)%pkz,isc, iec, jsc, jec, 0, km, 1.d00, gid==0)
  if (gid==0) print*,'-------------- State To FV 1 ----------------------------'
#endif

    FV_Atm(1)%u(isc:iec,jsc:jec,:) = STATE%VARS%U
    FV_Atm(1)%v(isc:iec,jsc:jec,:) = STATE%VARS%V
    call mpp_get_boundary(FV_Atm(1)%u, FV_Atm(1)%v, domain, &
                          wbuffery=wbuffer, ebuffery=FV_Atm(1)%v(iec+1,jsc:jec,1:km), &
                          sbufferx=sbuffer, nbufferx=FV_Atm(1)%u(isc:iec,jec+1,1:km), &
                          gridtype=DGRID_NE )

   if (SW_DYNAMICS) then
      do k=1,km+1
        do j=jsc,jec
          do i=isc,iec
            FV_Atm(1)%pe(i,k,j)   = STATE%VARS%PE(i,j,k)
          enddo
        enddo
      enddo
   else
      do k=1,km+1
        do j=jsc,jec
          do i=isc,iec
            FV_Atm(1)%pe(i,k,j)   = STATE%VARS%PE(i,j,k)
            FV_Atm(1)%pk(i,j,k)   = STATE%VARS%PE(i,j,k) ** akap
            FV_Atm(1)%peln(i,k,j) = log(STATE%VARS%PE(i,j,k))
          enddo
        enddo
      enddo
   endif

    do k=1,km
      do j=jsc,jec
        do i=isc,iec
          FV_Atm(1)%delp(i,j,k) = FV_Atm(1)%pe(i,k+1,j) - FV_Atm(1)%pe(i,k,j) 
        enddo
      enddo
    enddo

    if (.not. SW_DYNAMICS) then
       FV_Atm(1)%ps(isc:iec,jsc:jec) = STATE%VARS%PE(:,:,km+1)
       FV_Atm(1)%pkz(isc:iec,jsc:jec,:) = STATE%VARS%PKZ
       FV_Atm(1)%pt( isc:iec,jsc:jec,:) = STATE%VARS%PT*STATE%VARS%PKZ
    endif

#if defined(DEBUG)
  if (gid==0) print*,''
  if (gid==0) print*,'-------------- State To FV 2 ----------------------------'
  if (gid==0) print*,''
  call prt_maxmin('PS', FV_Atm(1)%ps, isc, iec, jsc, jec, ng,   1, 1.d-2, gid==0)
  call prt_maxmin('U ', FV_Atm(1)%u , isc, iec, jsc, jec+1, ng, km, 1.d00, gid==0)
  call prt_maxmin('V ', FV_Atm(1)%v , isc, iec+1, jsc, jec, ng, km, 1.d00, gid==0)
  call prt_maxmin('TA', FV_Atm(1)%pt, isc, iec, jsc, jec, ng, km, 1.d00, gid==0)
  call prt_maxmin('PZ', FV_Atm(1)%pkz,isc, iec, jsc, jec, 0, km, 1.d00, gid==0)
  if (gid==0) print*,'-------------- State To FV 2 ----------------------------'
#endif

   return

end subroutine State_To_FV

subroutine FV_To_State ( STATE )

!
! !INPUT PARAMETERS:

   type(T_FVDYCORE_STATE),      pointer   :: STATE

    integer               :: ISC,IEC, JSC,JEC, IM,JM,KM
    integer               :: I,J,K,N

    character(len=ESMF_MAXSTR)         :: IAm="FV_To_State"

    ISC = state%grid%is
    IEC = state%grid%ie
    JSC = state%grid%js
    JEC = state%grid%je
    KM  = state%grid%npz

! Copy updated FV data to internal state
    STATE%VARS%U(:,:,:) = FV_Atm(1)%u(isc:iec,jsc:jec,:)
    STATE%VARS%V(:,:,:) = FV_Atm(1)%v(isc:iec,jsc:jec,:)
    if (SW_DYNAMICS) then
       STATE%VARS%PE(:,:,1) = FV_Atm(1)%phis(isc:iec,jsc:jec)
       STATE%VARS%PE(:,:,2) = FV_Atm(1)%phis(isc:iec,jsc:jec) + FV_Atm(1)%delp(isc:iec,jsc:jec,1)
    else
       do j=jsc,jec
          do i=isc,iec
             STATE%VARS%PE(i,j,:) = FV_Atm(1)%pe(i,:,j)
          enddo
       enddo
       STATE%VARS%PKZ = FV_Atm(1)%pkz(isc:iec,jsc:jec,:)
       STATE%VARS%PT  = FV_Atm(1)%pt(isc:iec,jsc:jec,:)
       STATE%VARS%PT  = STATE%VARS%PT/STATE%VARS%PKZ
    endif
   return

end subroutine FV_To_State


subroutine a2d3d(grid, ua, va, ud, vd)

! Move A-Grid winds/tendencies oriented on lat/lon to the D-grid cubed-sphere orientation

      type (T_FVDYCORE_GRID), intent(in) :: grid
! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: ua(grid%isd:grid%ied  ,grid%jsd:grid%jed  ,grid%npz) ! U-Wind
      real(r8), intent(inout) :: va(grid%isd:grid%ied  ,grid%jsd:grid%jed  ,grid%npz) ! V-Wind
      real(r8), intent(inout) :: ud(grid%isd:grid%ied  ,grid%jsd:grid%jed+1,grid%npz) ! U-Wind
      real(r8), intent(inout) :: vd(grid%isd:grid%ied+1,grid%jsd:grid%jed  ,grid%npz) ! V-Wind
! !Local Variables
      integer :: isd,ied, jsd,jed
      integer :: is ,ie , js ,je 
      integer :: npx, npy, npz
      integer :: i,j,k, im2,jm2
      real(r8) :: p1(2), p2(2), p3(2), p4(2)

      real(r8) :: v3(grid%is-1:grid%ie+1,grid%js-1:grid%je+1,3)
      real(r8) :: ue(grid%is-1:grid%ie+1,grid%js  :grid%je+1,3)    ! 3D winds at edges
      real(r8) :: ve(grid%is  :grid%ie+1,grid%js-1:grid%je+1,3)    ! 3D winds at edges
      real(r8), dimension(grid%is:grid%ie):: ut1, ut2, ut3
      real(r8), dimension(grid%js:grid%je):: vt1, vt2, vt3

      npx = grid%npx
      npy = grid%npy
      npz = grid%npz
      isd = grid%isd
      ied = grid%ied
      jsd = grid%jsd
      jed = grid%jed
      is  = grid%is
      ie  = grid%ie
      js  = grid%js
      je  = grid%je

#if !defined(OLD_A2D)

      call cubed_a2d(npx, npy, npz, ua, va, ud, vd )

#else
      im2 = (npx-1)/2
      jm2 = (npy-1)/2

    call mpp_update_domains(ua, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
    call mpp_update_domains(va, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)

    do k=1, npz
! Compute 3D wind tendency on A grid
       do j=js-1,je+1
          do i=is-1,ie+1
             v3(i,j,1) = ua(i,j,k)*vlon(i,j,1) + va(i,j,k)*vlat(i,j,1)
             v3(i,j,2) = ua(i,j,k)*vlon(i,j,2) + va(i,j,k)*vlat(i,j,2)
             v3(i,j,3) = ua(i,j,k)*vlon(i,j,3) + va(i,j,k)*vlat(i,j,3)
          enddo
       enddo

! A --> D
! Interpolate to cell edges
       do j=js,je+1
          do i=is-1,ie+1
             ue(i,j,1) = v3(i,j-1,1) + v3(i,j,1)
             ue(i,j,2) = v3(i,j-1,2) + v3(i,j,2)
             ue(i,j,3) = v3(i,j-1,3) + v3(i,j,3)
          enddo
       enddo

       do j=js-1,je+1
          do i=is,ie+1
             ve(i,j,1) = v3(i-1,j,1) + v3(i,j,1)
             ve(i,j,2) = v3(i-1,j,2) + v3(i,j,2)
             ve(i,j,3) = v3(i-1,j,3) + v3(i,j,3)
          enddo
       enddo

! --- E_W edges (for v-wind):
     if ( is==1 ) then
       i = 1
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_w(j)*ve(i,j-1,1)+(1.-edge_vect_w(j))*ve(i,j,1)
             vt2(j) = edge_vect_w(j)*ve(i,j-1,2)+(1.-edge_vect_w(j))*ve(i,j,2)
             vt3(j) = edge_vect_w(j)*ve(i,j-1,3)+(1.-edge_vect_w(j))*ve(i,j,3)
        else
             vt1(j) = edge_vect_w(j)*ve(i,j+1,1)+(1.-edge_vect_w(j))*ve(i,j,1)
             vt2(j) = edge_vect_w(j)*ve(i,j+1,2)+(1.-edge_vect_w(j))*ve(i,j,2)
             vt3(j) = edge_vect_w(j)*ve(i,j+1,3)+(1.-edge_vect_w(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
     if ( (ie+1)==npx ) then
       i = npx
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_e(j)*ve(i,j-1,1)+(1.-edge_vect_e(j))*ve(i,j,1)
             vt2(j) = edge_vect_e(j)*ve(i,j-1,2)+(1.-edge_vect_e(j))*ve(i,j,2)
             vt3(j) = edge_vect_e(j)*ve(i,j-1,3)+(1.-edge_vect_e(j))*ve(i,j,3)
        else
             vt1(j) = edge_vect_e(j)*ve(i,j+1,1)+(1.-edge_vect_e(j))*ve(i,j,1)
             vt2(j) = edge_vect_e(j)*ve(i,j+1,2)+(1.-edge_vect_e(j))*ve(i,j,2)
             vt3(j) = edge_vect_e(j)*ve(i,j+1,3)+(1.-edge_vect_e(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
! N-S edges (for u-wind):
     if ( js==1 ) then
       j = 1
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_s(i)*ue(i-1,j,1)+(1.-edge_vect_s(i))*ue(i,j,1)
             ut2(i) = edge_vect_s(i)*ue(i-1,j,2)+(1.-edge_vect_s(i))*ue(i,j,2)
             ut3(i) = edge_vect_s(i)*ue(i-1,j,3)+(1.-edge_vect_s(i))*ue(i,j,3)
        else
             ut1(i) = edge_vect_s(i)*ue(i+1,j,1)+(1.-edge_vect_s(i))*ue(i,j,1)
             ut2(i) = edge_vect_s(i)*ue(i+1,j,2)+(1.-edge_vect_s(i))*ue(i,j,2)
             ut3(i) = edge_vect_s(i)*ue(i+1,j,3)+(1.-edge_vect_s(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif
     if ( (je+1)==npy ) then
       j = npy
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_n(i)*ue(i-1,j,1)+(1.-edge_vect_n(i))*ue(i,j,1)
             ut2(i) = edge_vect_n(i)*ue(i-1,j,2)+(1.-edge_vect_n(i))*ue(i,j,2)
             ut3(i) = edge_vect_n(i)*ue(i-1,j,3)+(1.-edge_vect_n(i))*ue(i,j,3)
        else
             ut1(i) = edge_vect_n(i)*ue(i+1,j,1)+(1.-edge_vect_n(i))*ue(i,j,1)
             ut2(i) = edge_vect_n(i)*ue(i+1,j,2)+(1.-edge_vect_n(i))*ue(i,j,2)
             ut3(i) = edge_vect_n(i)*ue(i+1,j,3)+(1.-edge_vect_n(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif

! Update:
       do j=js,je+1
          do i=is,ie
             ud(i,j,k) = 0.5*( ue(i,j,1)*es(1,i,j,1) +  &
                               ue(i,j,2)*es(2,i,j,1) +  &
                               ue(i,j,3)*es(3,i,j,1) )
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             vd(i,j,k) = 0.5*( ve(i,j,1)*ew(1,i,j,2) +  &
                               ve(i,j,2)*ew(2,i,j,2) +  &
                               ve(i,j,3)*ew(3,i,j,2) )
          enddo
       enddo

    enddo         ! k-loop

#endif

end subroutine a2d3d

subroutine fv_getTopography(phis)
  use fv_control_mod,     only: npx,npy
  use fv_grid_tools_mod,  only: grid, agrid, dxc, dyc
  use fv_surf_map_mod,    only: surfdrv
  real(kind=4), intent(INOUT) ::   phis(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec)
  real(r8)                    :: r8phis(FV_Atm(1)%isd:FV_Atm(1)%ied,FV_Atm(1)%jsd:FV_Atm(1)%jed)

  call surfdrv(npx, npy, grid, agrid, area, dx, dy, dxc, dyc,  &
               r8phis, gid==masterproc)
  call mpp_update_domains( r8phis, domain )

  phis = r8phis(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec)

end subroutine fv_getTopography

subroutine fv_getMassFluxes(mfx, mfy, mfz, dt)
  real(r8), intent(OUT) :: mfx(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)
  real(r8), intent(OUT) :: mfy(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)
  real(r8), intent(OUT) :: mfz(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz+1)
  real(r8), intent(IN)  :: dt
  integer isc,iec,jsc,jec,npz,i,j,k

  real(r8) :: tmp_mfx(FV_Atm(1)%isd:FV_Atm(1)%ied,FV_Atm(1)%jsd:FV_Atm(1)%jed,1:FV_Atm(1)%npz)
  real(r8) :: tmp_mfy(FV_Atm(1)%isd:FV_Atm(1)%ied,FV_Atm(1)%jsd:FV_Atm(1)%jed,1:FV_Atm(1)%npz)
  real(r8) :: conv(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)
  real(r8) :: pit(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec)

  real(r8) :: fac

  isc=FV_Atm(1)%isc ; iec=FV_Atm(1)%iec
  jsc=FV_Atm(1)%jsc ; jec=FV_Atm(1)%jec
  npz=FV_Atm(1)%npz

  fac = 1.0/(dt*grav)

! call cubed_to_latlon(FV_Atm(1)%mfx, FV_Atm(1)%mfy, &
!                            tmp_mfx,       tmp_mfy, &
!                      dx, dy, rdxa, rdya, npz)
  do k=1,npz
     mfx(:,:,k) = FV_Atm(1)%mfx(isc:iec,jsc:jec,k) * fac  ! Pa m^2 / s
     mfy(:,:,k) = FV_Atm(1)%mfy(isc:iec,jsc:jec,k) * fac  ! Pa m^2 / s
  enddo 
!
! Compute the vertical mass flux
!
!   Compute Convergence of the horizontal Mass flux
    do k=1,npz
       do j=jsc,jec
          do i=isc,iec
             conv(i,j,k) = ( FV_Atm(1)%mfx(i,j,k) - FV_Atm(1)%mfx(i+1,j,k) +  &
                             FV_Atm(1)%mfy(i,j,k) - FV_Atm(1)%mfy(i,j+1,k) ) * fac
          enddo
       enddo
    enddo
!   Surface pressure tendency
    pit(:,:) = 0.0
    do k=1,npz
       do j=jsc,jec
          do i=isc,iec
             pit(i,j) = pit(i,j) + conv(i,j,k)
          enddo
       enddo
    enddo
!   Sum over levels
    do k=2,npz
       do j=jsc,jec
          do i=isc,iec
             conv(i,j,k) = conv(i,j,k) + conv(i,j,k-1)
          enddo
       enddo
    enddo
    mfz(:,:,:) = 0.0
    do k=2,npz
       do j=jsc,jec
          do i=isc,iec
             mfz(i,j,k) = ( conv(i,j,k-1)  - FV_Atm(1)%bk(k)*pit(i,j) )/(grav*area(i,j))  ! Kg/m^2/s
          enddo
       enddo
    enddo

return
end subroutine fv_getMassFluxes

subroutine fv_getOmega(omga)
  real(r8), intent(OUT) :: omga(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)
  integer isc,iec,jsc,jec

  isc=FV_Atm(1)%isc ; iec=FV_Atm(1)%iec
  jsc=FV_Atm(1)%jsc ; jec=FV_Atm(1)%jec

  omga(:,:,:) = FV_Atm(1)%omga(isc:iec,jsc:jec,:)

return
end subroutine fv_getOmega

subroutine fv_getPK(pkxyz)
  real(r8), intent(OUT) :: pkxyz(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz+1)
  integer isc,iec,jsc,jec
  isc=FV_Atm(1)%isc ; iec=FV_Atm(1)%iec
  jsc=FV_Atm(1)%jsc ; jec=FV_Atm(1)%jec
  pkxyz(:,:,:) = FV_Atm(1)%pk(isc:iec,jsc:jec,:)
  return
end subroutine fv_getPK

subroutine fv_getVorticity(u, v, vort)
  real(r8), intent(IN)  ::  u(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)
  real(r8), intent(IN)  ::  v(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)
  real(r8), intent(OUT) ::  vort(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)

  real(r8) :: wbuffer(FV_Atm(1)%npy+2,FV_Atm(1)%npz)
  real(r8) :: sbuffer(FV_Atm(1)%npx+2,FV_Atm(1)%npz)

  integer isc,iec,jsc,jec
  integer npz
  integer i,j,k

  isc=FV_Atm(1)%isc ; iec=FV_Atm(1)%iec
  jsc=FV_Atm(1)%jsc ; jec=FV_Atm(1)%jec
  npz = FV_Atm(1)%npz

  FV_Atm(1)%u(isc:iec,jsc:jec,:) = u
  FV_Atm(1)%v(isc:iec,jsc:jec,:) = v
  call mpp_get_boundary(FV_Atm(1)%u, FV_Atm(1)%v, domain, &
                        wbuffery=wbuffer, ebuffery=FV_Atm(1)%v(iec+1,jsc:jec,1:npz), &
                        sbufferx=sbuffer, nbufferx=FV_Atm(1)%u(isc:iec,jec+1,1:npz), &
                        gridtype=DGRID_NE )
! Calc Vorticity
! Convert winds to circulation elements:
    do k=1,npz
       do j=jsc,jec
          do i=isc,iec
             vort(i,j,k) = rarea(i,j)*(FV_Atm(1)%u(i,j,k)*dx(i,j)-FV_Atm(1)%u(i,j+1,k)*dx(i,j+1) - &
                                       FV_Atm(1)%v(i,j,k)*dy(i,j)+FV_Atm(1)%v(i+1,j,k)*dy(i+1,j))
          enddo
       enddo
    enddo
end subroutine fv_getVorticity

subroutine fv_getAgridWinds(u, v, ua, va, rotate)
  real(r8), intent(IN)  ::  u(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)
  real(r8), intent(IN)  ::  v(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)
  real(r8), intent(OUT) :: ua(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)
  real(r8), intent(OUT) :: va(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)
  logical, optional, intent(IN) :: rotate

  real(r8) :: wbuffer(FV_Atm(1)%npy+2,FV_Atm(1)%npz)
  real(r8) :: sbuffer(FV_Atm(1)%npx+2,FV_Atm(1)%npz)

  integer isc,iec,jsc,jec
  integer npz
  integer i,j,k

  real(r8) :: p1(2), p2(2), p3(2), p4(2) 

  real(r8) :: ut(FV_Atm(1)%isd:FV_Atm(1)%ied, FV_Atm(1)%jsd:FV_Atm(1)%jed)
  real(r8) :: vt(FV_Atm(1)%isd:FV_Atm(1)%ied, FV_Atm(1)%jsd:FV_Atm(1)%jed)

  isc=FV_Atm(1)%isc ; iec=FV_Atm(1)%iec
  jsc=FV_Atm(1)%jsc ; jec=FV_Atm(1)%jec
  npz = FV_Atm(1)%npz

  FV_Atm(1)%u(isc:iec,jsc:jec,:) = u
  FV_Atm(1)%v(isc:iec,jsc:jec,:) = v
  call mpp_get_boundary(FV_Atm(1)%u, FV_Atm(1)%v, domain, &
                        wbuffery=wbuffer, ebuffery=FV_Atm(1)%v(iec+1,jsc:jec,1:npz), &
                        sbufferx=sbuffer, nbufferx=FV_Atm(1)%u(isc:iec,jec+1,1:npz), &
                        gridtype=DGRID_NE )
  call mpp_update_domains(FV_Atm(1)%u, FV_Atm(1)%v, domain, gridtype=DGRID_NE)

  do k=1,npz
     call d2a2c_vect(FV_Atm(1)%u(:,:,k),  FV_Atm(1)%v(:,:,k), &
                     FV_Atm(1)%ua(:,:,k), FV_Atm(1)%va(:,:,k), &
                     FV_Atm(1)%uc(:,:,k), FV_Atm(1)%vc(:,:,k), ut, vt, .false.)
  enddo
  if (present(rotate)) then
     if (rotate) call cubed_to_latlon(FV_Atm(1)%u  , FV_Atm(1)%v  , &
                                      FV_Atm(1)%ua , FV_Atm(1)%va , &
                                      dx, dy, rdxa, rdya, npz)
  endif
  ua(:,:,:) = FV_Atm(1)%ua(isc:iec,jsc:jec,:) 
  va(:,:,:) = FV_Atm(1)%va(isc:iec,jsc:jec,:)

  return
end subroutine fv_getAgridWinds

subroutine latlon_to_cubed_winds(ua, va)
  real(r8), intent(INOUT)  ::  ua(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)
  real(r8), intent(INOUT)  ::  va(FV_Atm(1)%isc:FV_Atm(1)%iec,FV_Atm(1)%jsc:FV_Atm(1)%jec,1:FV_Atm(1)%npz)

  real(r8) :: p1(2), p2(2), p3(2), p4(2)

  integer :: i,j,k
  integer :: IS,IE,JS,JE,KM
 
  IS = FV_Atm(1)%isc
  IE = FV_Atm(1)%iec
  JS = FV_Atm(1)%jsc
  JE = FV_Atm(1)%jec
  KM = FV_Atm(1)%npz

  do k=1,km
    do j=js,je
      do i=is,ie
         call mid_pt_sphere(FV_Atm(1)%grid(i,j,1:2), FV_Atm(1)%grid(i,j+1,1:2), p1)
         call mid_pt_sphere(FV_Atm(1)%grid(i,j,1:2), FV_Atm(1)%grid(i+1,j,1:2), p2)
         call mid_pt_sphere(FV_Atm(1)%grid(i+1,j,1:2), FV_Atm(1)%grid(i+1,j+1,1:2), p3)
         call mid_pt_sphere(FV_Atm(1)%grid(i,j+1,1:2), FV_Atm(1)%grid(i+1,j+1,1:2), p4)
         call rotate_winds(ua(i,j,k), va(i,j,k), p1,p2,p3,p4, FV_Atm(1)%agrid(i,j,1:2), 2, 1)
      enddo
    enddo
  enddo

end subroutine latlon_to_cubed_winds

  subroutine CREATE_VARS (I1, IN, J1, JN, K1, KN, KP, &
       U, V, PT, PE, PKZ, VARS )

    integer, intent(IN   ) :: I1, IN, J1, JN, K1, KN, KP
    real(r8), target ::   U(I1:IN,J1:JN,K1:KN  )
    real(r8), target ::   V(I1:IN,J1:JN,K1:KN  )
    real(r8), target ::  PT(I1:IN,J1:JN,K1:KN  )
    real(r8), target ::  PE(I1:IN,J1:JN,K1:KP  )
    real(r8), target :: PKZ(I1:IN,J1:JN,K1:KN  )

    type (T_FVDYCORE_VARS), intent(INOUT) :: VARS

    VARS%U => U
    VARS%V => V
    VARS%PT => PT
    VARS%PE => PE
    VARS%PKZ => PKZ

    return
  end subroutine CREATE_VARS

      subroutine get_latlon_vector (pp, elon, elat)
      real(r8), intent(IN)  :: pp(2)
      real(r8), intent(OUT) :: elon(3), elat(3)

         elon(1) = -SIN(pp(1))
         elon(2) =  COS(pp(1))
         elon(3) =  0.0
         elat(1) = -SIN(pp(2))*COS(pp(1))
         elat(2) = -SIN(pp(2))*SIN(pp(1))
#ifdef RIGHT_HAND
         elat(3) =  COS(pp(2))
#else
! Left-hand system needed to be consistent with rest of the codes
         elat(3) = -COS(pp(2))
#endif

      end subroutine get_latlon_vector

end module FV_StateMod

