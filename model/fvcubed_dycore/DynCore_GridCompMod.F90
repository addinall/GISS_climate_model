#include "MAPL_Generic.h"

#define USE_FVcubed
#define DynCore_GridCompMod FVcubed_dycore_GridCompMod

!-----------------------------------------------------------------------
!              ESMA - Earth System Modeling Applications
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: DynCore_GridCompMod --- Dynamical Core Grid Component
!
! !INTERFACE:

   Module DynCore_GridCompMod

! !USES:

   use ESMF_Mod            ! ESMF base class
   use MAPL_Mod            ! GEOS base class

#if defined(USE_FV)
   use dynamics_vars, only : DynTracers      => T_TRACERS,        &
                             DynVars         => T_FVDYCORE_VARS,  &
                             DynGrid         => T_FVDYCORE_GRID,  &
                             DynState        => T_FVDYCORE_STATE
   use FV_StateMod, only :   DynInit         => FV_InitState,       &
                             DynFinalize     => FV_Finalize,        &
                             getAgridWinds   => stateGetAgridWinds, &
                             getOmega        => stateGetOmega,      &
                             getPK           => stateGetPK,         &
                             getVorticity    => stateGetVorticity,  &
                             DYN_COLDSTART   => COLDSTART,          &
                             SW_DYNAMICS, PURE_ADVECTION, NXQ
#endif

#if defined(USE_FVcubed)
! FV Specific Module
   use FV_StateMod, only : DynTracers      => T_TRACERS,        &
                           DynVars         => T_FVDYCORE_VARS,  &
                           DynGrid         => T_FVDYCORE_GRID,  &
                           DynState        => T_FVDYCORE_STATE, &
                           DynInit         => FV_InitState,     &
                           DynRun          => FV_Run,           &
                           DynFinalize     => FV_Finalize,      &
                           getAgridWinds   => fv_getAgridWinds, &
                           getMassFluxes   => fv_getMassFluxes, &
                           getOmega        => fv_getOmega,      &
                           getPK           => fv_getPK,         &
                           getVorticity    => fv_getVorticity,  &
                           Agrid_To_Native => a2d3d,            &
                           DYN_COLDSTART   => COLDSTART,        &
                           SW_DYNAMICS, PURE_ADVECTION, NXQ
#endif

#if defined(USE_GISS)
   use GISS_StateMod, only : DynTracers      => T_TRACERS,          &
                             DynVars         => T_DYCORE_VARS,      &
                             DynGrid         => T_DYCORE_GRID,      &
                             DynState        => T_DYCORE_STATE,     &
                             DynInit         => INITIALIZE,         &
                             DynRun          => RUN,                &
                             DynFinalize     => FINALIZE,           &
                             getAgridWinds   => stateGetAgridWinds, &
                             getOmega        => stateGetOmega,      &
                             getPK           => stateGetPK,         &
                             getVorticity    => stateGetVorticity,  &
                             DYN_COLDSTART   => COLDSTART,          &
                             SW_DYNAMICS, PURE_ADVECTION, NXQ
#endif

! !PUBLIC MEMBER FUNCTIONS:

  implicit none
  private

  public  SetServices      ! Register component methods

! !DESCRIPTION: This module implements the Dynamical Core as
!               an ESMF gridded component.
!
! \paragraph{Overview}
!
!   This module contains an ESMF wrapper for a generic
!   Dynamical Core.
!
! \paragraph{Internal State}
!
!  FVdycore maintains an internal state consisting of the
!  following fields:  control variables
!
!   \begin{itemize}
!     \item {\tt U}:    U winds on the native grid  (m/s)
!     \item {\tt V}:    V winds on the native grid (m/s)
!     \item {\tt PT}:   Scaled Virtual Potential Temperature (T_v/PKZ)
!     \item {\tt PE}:   Edge pressures
!     \item {\tt Q}:    Tracers
!     \item {\tt PKZ}:  Consistent mean for p^kappa
!   \end{itemize}
!
!  as well as a GRID (to be mentioned later) 
!  and same additional run-specific variables 
!
! Note: {\tt PT} is not updated if the flag {\tt CONVT} is true.
!
! The internal state is updated each time FVdycore is called.
!
! \paragraph{Import State}
!
! The import state consists of the tendencies of the 
! control variables plus the surface geopotential heights:
!
!   \begin{itemize}
!     \item {\tt DUDT}:    U wind tendency on a A-grid (m/s)
!     \item {\tt DVDT}:    V wind tendency on a A-grid (m/s)
!     \item {\tt DTDT}:    Delta-pressure-weighted temperature tendency
!     \item {\tt DPEDT}:   Edge pressure tendency
!     \item {\tt PHIS}:    Surface Geopotential Heights
!   \end{itemize}
!
! These are by definition on an A-grid and have an XY
! domain decomposition.
!
! \paragraph{Export State}
!
!   The export state can provide the following variables:
!
!   \begin{itemize}
!     \item {\tt U}:         U winds on a A-grid (m/s)
!     \item {\tt V}:         V winds on a A-grid (m/s)
!     \item {\tt U_CGRID}:   U winds on a C-grid (m/s)
!     \item {\tt V_CGRID}:   V winds on a C-grid (m/s)
!     \item {\tt U_DGRID}:   U winds on a D-grid (m/s)
!     \item {\tt V_DGRID}:   V winds on a D-grid (m/s)
!     \item {\tt T}:         Temperature (K)
!     \item {\tt Q}:         Tracers
!     \item {\tt TH}:        Potential Temperature (K)
!     \item {\tt ZL}:        Mid-Layer Heights (m)
!     \item {\tt ZLE}:       Edge Heights (m)
!     \item {\tt PLE}:       Edge pressures (Pa)
!     \item {\tt PLK}:       P^Kappa at Mid-Layers
!     \item {\tt OMEGA}:     Vertical pressure velocity (pa/s)
!     \item {\tt MFX_UR}:    Mass-Weighted U-Wind on C-Grid (Pa m^2/s)
!     \item {\tt MFY_UR}:    Mass-Weighted V-wind on C-Grid (Pa m^2/s)
!     \item {\tt MFX}:       Remapped Mass-Weighted U-Wind on C-Grid (Pa m^2/s)
!     \item {\tt MFY}:       Remapped Mass-Weighted V-wind on C-Grid (Pa m^2/s)
!     \item {\tt MFZ}:       Remapped Vertical mass flux (kg/(m^2*s))
!     \item {\tt MFX_A}:     Remapped Mass-Weighted U-Wind on A-Grid (Pa m/s)
!     \item {\tt MFY_A}:     Remapped Mass-Weighted V-wind on A-Grid (Pa m/s)
!     \item {\tt PV}:        Ertel's Potential Vorticity (m^2 / kg*s)
!     \item {\tt DUDT}:      U-wind Tendency (m/s/s)
!     \item {\tt DVDT}:      V-wind Tendency (m/s/s)
!     \item {\tt DTDT}:      Mass-Weighted Temperature Tendency (Pa K/s)
!   \end{itemize}
!
!   All variables are on an A-grid with points at the poles, and have an XY decomposition.
!
! \paragraph{Grids and Decompositions}
!
!   The current version supports only a 1D latitude-based
!   decomposition of the domain (with OMP task-parallelism
!   in the vertical, resulting in reasonable scalability 
!   on large PE configurations).  In the near future it will 
!   support a 2D domain decomposition, in which import and
!   export state are decomposed in longitude and latitude,
!   while the internal state (for the most part) is 
!   decomposed in latitude and level.  When needed, 
!   the data is redistributed (``transposed'') internally.
!
!   There are two fundamental ESMF grids in use;
!   \begin{itemize}
!     \item {GRIDXY}: longitude-latitude ESMF grid (public)
!     \item {GRIDYZ}: A latitude-level cross-sectional
!                     decomposition (private to this module) 
!   \end{itemize}
!
!   PILGRIM will be used for communication until ESMF has 
!   sufficient functionality and performance to take over 
!   the task.  The use of pilgrim requires a call to 
!   {\tt INIT\_SPMD} to set SPMD parameters, decompositions,
!   etc.
!
! \paragraph{Required Files}
!
!  The following files are needed for a standard restart run:
!
!  \begin{itemize}
!    \item Layout file
!      \begin{itemize}
!        \item {\tt nprxy_x, nprxy_y, npryz_y, npryz_z}:
!          process dimensions in XY and YZ.
!        \item {\tt imxy, jmxy, jmyz, kmyz}: distributions for XY and YZ
!        \item {\tt iord, jord}: the order of the lon. and lat. algorithms
!        \item {\tt dtime}:  The large (advection) time step
!        \item {\tt nsplit}: the ratio between the large and small time step
!          (possibly zero for automatic determination),
!      \end{itemize}
!    \item Restart file
!      \begin{itemize}
!        \item date in standard format yy, mm, dd, hh, mm, ss
!        \item dimensions im, jm, km, nq
!        \item control variables {\tt U, V, PT, PE, Q}
!      \end{itemize}
!    \item Topography file
!
!  \end{itemize}
!
! \paragraph{Future Additions}
!
!  \begin{itemize}
!    \item  Conservation of energy (CONSV  == .TRUE. )
!    \item  2D decomposition (requires transposes in the coupler)
!    \item  Use r8 instead of r4 (currently supported in StopGap)
!  \end{itemize}
!
! !REVISION HISTORY:
!
! 11Jul2003  Sawyer    From Trayanov/da Silva EVAC 
! 23Jul2003  Sawyer    First informal tiptoe-through
! 29Jul2003  Sawyer    Modifications based on comments from 23Jul2003
! 28Aug2003  Sawyer    First check-in; Internal state to D-grid
! 15Sep2003  Sawyer    Extensive bug fixes, revisions
! 24Sep2003  Sawyer    Modified names; corrected weighting of T, Q
! 22Oct2003  Sawyer    pmgrid removed (data now in spmd\_dyn)
! 25Nov2003  Sawyer    Optimization for 1D decomposition
! 03Dec2003  Sawyer    Switched over to specified decompositions
! 04Dec2003  Sawyer    Moved T_FVDYCORE_GRID to dynamics_vars
! 21Jan2004  Takacs    Modified Import/Export, Added Generic State, Added TOPO utility
! 20Sep2004  Sawyer    Revised cd_core, trac2d interfaces, refactoring
! 06Oct2004  Sawyer    More refactoring, removed spmd_dyn
! 17Feb2005  Sawyer    Added Ertel's potential vorticity to diagnostics
! 20Mar2005  Sawyer    Tracers are now pointers into import state
! 12Apr2005  Sawyer    Extensive changes to minimize tracer memory
! 18May2005  Sawyer    Put FVdycore_wrapper in separate file; CAM/GEOS5 merge
! 16Nov2005  Takacs    Added option for DCADJ, Merge with Daedalus_p5
! 18Jan2006  Putman    Added mass fluxes to export state
!
!EOP
!----------------------------------------------------------------------
!BOC

  integer,  parameter :: r8           = 8
  integer,  parameter :: r4           = 4

  real(r8), parameter :: RADIUS       = MAPL_RADIUS
  real(r8), parameter :: CP           = MAPL_CP
  real(r8), parameter :: PI           = MAPL_PI_R8
  real(r8), parameter :: OMEGA        = MAPL_OMEGA
  real(r8), parameter :: KAPPA        = MAPL_KAPPA
  real(r8), parameter :: P00          = MAPL_P00
  real(r8), parameter :: GRAV         = MAPL_GRAV
  real(r8), parameter :: RGAS         = MAPL_RGAS
  real(r8), parameter :: RVAP         = MAPL_RVAP
  real(r8), parameter :: EPS          = RVAP/RGAS-1.0

  integer,  parameter :: TIME_TO_RUN  = 1
  integer,  parameter :: CHECK_MAXMIN = 2

  integer :: I, J, K  !  Default declaration for loops.

! Tracer I/O History stuff
! -------------------------------------
    integer, parameter         :: nlevs=5
    integer, parameter         :: ntracers=6                       
    integer                    :: nlev, ntracer                    
    integer                    :: plevs(nlevs)          
    character(len=ESMF_MAXSTR) :: myTracer
    data plevs /850,700,600,500,300/


! Wrapper for extracting internal state
! -------------------------------------

  type DYN_wrap
     type (DynState), pointer :: DYN_STATE
  end type DYN_wrap

  interface addTracer
     module procedure addTracer_r4
     module procedure addTracer_r8
  end interface

!#define DEBUG
#if defined(DEBUG)         
  interface Write_Profile
     module procedure Write_Profile_R4
     module procedure Write_Profile_R8
  end interface
#endif

contains

!----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices

! !DESCRIPTION:  SetServices registers Initialize, Run, and Finalize
!   methods for FV. Two stages of the FV run method are registered. The
!   first one does the dynamics calculations, and the second adds 
!   increments from external sources that appear in the Import state.
!   SetServices also creates a private internal state in which FV
!   keeps invariant or auxilliary state variables, as well as pointers to
!   the true state variables. The MAPL internal state contains the
!   true state variables and is managed by MAPL.
!
! !INTERFACE:

   Subroutine SetServices ( gc, rc )

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: gc     ! gridded component
   integer,             intent(out)   :: rc     ! return code
!         
!EOP         
!----------------------------------------------------------------------
  
   type (DynState), pointer :: dyn_internal_state 
    type (DYN_wrap)                  :: wrap

    integer                          :: status
    character(len=ESMF_MAXSTR)       :: IAm
    character(len=ESMF_MAXSTR)       :: COMP_NAME

    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "SetServices"
    
! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( dyn_internal_state, stat=status )
    VERIFY_(STATUS)
    wrap%dyn_state => dyn_internal_state
 
! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC,'DYNstate',wrap,status )
    VERIFY_(STATUS)

!BOP

! !IROUTINE: State Descriptions

! !DESCRIPTION: The component uses all three states (Import, Export
!  and Internal), in addition to a Private (non-ESMF) Internal state. All
!  three are managed by MAPL. 
!
!  The Private Internal state contains invariant
!  quantities defined by an FV specific routine, as well as pointers 
!  to the true state variables, kept in the MAPL Internal state. 
!  The MAPL Internal is kept at FV's real*8 precision.
!
!  The Import State conatins tendencies to be added in the second 
!  run stage, the geopotential at the lower boundary, and a bundle
!  of Friendly tracers to be advected. The Import and Export states
!  are both at the default precision.
!

! !IMPORT STATE:

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_tendency',                    &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_tendency',                   &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'delta-p_weighted_temperature_tendency',     &
         UNITS      = 'Pa K s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQVANA',                                    &
         LONG_NAME  = 'specific_humidity_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_tendency',                    &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'PHIS',                                      &
         LONG_NAME  = 'surface_geopotential_height',               &
         UNITS      = 'm+2 sec-2',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec( gc,                              &
        SHORT_NAME = 'QTR',                                        &
        LONG_NAME  = 'advected_quantities',                        &
        UNITS      = 'unknown',                                    &
        VLOCATION  = MAPL_VLocationCenter,               &
        DATATYPE   = MAPL_BundleItem,               &
        RC=STATUS  )
    VERIFY_(STATUS)

! !EXPORT STATE:

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'TAVE',                                                               &
         LONG_NAME  = 'vertically_averaged_dry_temperature',                                &
         UNITS      = 'K',                                                                  &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'UAVE',                                                               &
         LONG_NAME  = 'vertically_averaged_zonal_wind',                                     &
         UNITS      = 'm sec-1',                                                            &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'KEPHY',                                                              &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_physics',       &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
       SHORT_NAME   = 'PEPHY',                                      &
       LONG_NAME    = 'total_potential_energy_tendency_due_to_physics',          &
       UNITS        = 'W m-2',                                      &
       DIMS         = MAPL_DimsHorzOnly,                            &
       VLOCATION    = MAPL_VLocationNone,                 RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                   &
       SHORT_NAME   = 'TEPHY',                                      &
       LONG_NAME    = 'mountain_work_tendency_due_to_physics',      &
       UNITS        = 'W m-2',                                      &
       DIMS         = MAPL_DimsHorzOnly,                            &
       VLOCATION    = MAPL_VLocationNone,                 RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                    &
       SHORT_NAME         = 'KEANA',                                   &
       LONG_NAME          = 'total_kinetic_energy_tendency_due_to_analysis',&
       UNITS              = 'W m-2',                                      &
       DIMS               = MAPL_DimsHorzOnly,                            &
       VLOCATION          = MAPL_VLocationNone,                 RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                    &
       SHORT_NAME         = 'PEANA',                                   &
       LONG_NAME          = 'total_potential_energy_tendency_due_to_analysis',&
       UNITS              = 'W m-2',                                      &
       DIMS               = MAPL_DimsHorzOnly,                            &
       VLOCATION          = MAPL_VLocationNone,                 RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                    &
       SHORT_NAME         = 'TEANA',                                   &
       LONG_NAME          = 'mountain_work_tendency_due_to_analysis',&
       UNITS              = 'W m-2',                                      &
       DIMS               = MAPL_DimsHorzOnly,                            &
       VLOCATION          = MAPL_VLocationNone,                 RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'KEDYN',                                                              &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_across_dynamics',      &
         UNITS      = 'Pa m s-1',                                                           &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                       &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'PEDYN',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_across_dynamics',    &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'TEDYN',                                                              &
         LONG_NAME  = 'mountain_work_tendency_across_dynamics',                             &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'KECDCOR',                                                            &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_across_cdcore',        &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'PECDCOR',                                                            &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_across_cdcore',      &
         UNITS      = 'Pa m s-1',                                                           &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                       &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'TECDCOR',                                                            &
         LONG_NAME  = 'mountain_work_tendency_across_cdcore',                               &
         UNITS      = 'Pa m s-1',                                                           &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                       &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'KETEMAP',                                                            &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_across_temap',         &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'PETEMAP',                                                            &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_across_temap',       &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'TETEMAP',                                                            &
         LONG_NAME  = 'mountain_work_tendency_across_temap',                                &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'KEGEN',                                                              &
         LONG_NAME  = 'vertically_integrated_generation_of_kinetic_energy',                 &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                       &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'DQVDTINT',                                                           &
         LONG_NAME  = 'vertically_integrated_water_vapor_tendency_across_dynamics',         &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'DMDTDYN',                                                            &
         LONG_NAME  = 'vertically_integrated_mass_tendency_across_dynamics',                &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'DQLDTINT',                                                           &
         LONG_NAME  = 'vertically_integrated_liquid_water_tendency_across_dynamics',        &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'DQIDTINT',                                                           &
         LONG_NAME  = 'vertically_integrated_ice_water_tendency_across_dynamics',           &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                      &
         SHORT_NAME = 'DOXDTINT',                                                           &
         LONG_NAME  = 'vertically_integrated_ozone_tendency_across_dynamics',               &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                               &
         SHORT_NAME = 'CONVKE',                                                      &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_convergence',            &
         UNITS      = 'W m-2',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                             &
         VLOCATION  = MAPL_VLocationCenter,                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                               &
         SHORT_NAME = 'CONVCPT',                                                     &
         LONG_NAME  = 'vertically_integrated_enthalpy_convergence',                  &
         UNITS      = 'W m-2',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                             &
         VLOCATION  = MAPL_VLocationCenter,                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                               &
         SHORT_NAME = 'CONVPHI',                                                     &
         LONG_NAME  = 'vertically_integrated_geopotential_convergence',              &
         UNITS      = 'W m-2',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                             &
         VLOCATION  = MAPL_VLocationCenter,                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VORT',                                      &
         LONG_NAME  = 'vorticity',                                 &
         UNITS      = 's-1',                                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,            &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'T',                                         &
         LONG_NAME  = 'air_temperature',                           &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,            &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'PL',                                        &
         LONG_NAME  = 'mid_level_pressure',                        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'ZLE',                                       &
         LONG_NAME  = 'edge_heights',                              &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'ZL',                                        &
         LONG_NAME  = 'mid_layer_heights',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'S',                                         &
         LONG_NAME  = 'mid_layer_dry_static_energy',               &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,            &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'PLE',                                       &
         LONG_NAME  = 'edge_pressure',                             &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'TH',                                        &
         LONG_NAME  = 'potential_temperature',                     &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'PLK',                                       &
         LONG_NAME  = 'mid-layer_p^kappa',                         &
         UNITS      = 'Pa^kappa',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'OMEGA',                                     &
         LONG_NAME  = 'vertical_pressure_velocity',                &
         UNITS      = 'Pa sec-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                               &
         SHORT_NAME = 'MFX_UR',                                      &
         LONG_NAME  = 'pressure_weighted_eastward_wind_unremapped',  &
         UNITS      = 'Pa m s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,                          &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                               &
         SHORT_NAME = 'MFY_UR',                                      &
         LONG_NAME  = 'pressure_weighted_northward_wind_unremapped', &
         UNITS      = 'Pa m s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,                          &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'MFX',                                       &
         LONG_NAME  = 'pressure_weighted_eastward_wind',           &
         UNITS      = 'Pa m+2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'MFY',                                       &
         LONG_NAME  = 'pressure_weighted_northward_wind',          &
         UNITS      = 'Pa m+2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'MFZ',                                       &
         LONG_NAME  = 'vertical_mass_flux',                        &
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'MFX_A',                                     &
         LONG_NAME  = 'zonal_mass_flux',                           &
         UNITS      = 'Pa m s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,            &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'MFY_A',                                     &
         LONG_NAME  = 'meridional_mass_flux',                      &
         UNITS      = 'Pa m s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,            &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'PV',                                        &
         LONG_NAME  = 'ertels_isentropic_potential_vorticity',                &
         UNITS      = 'm+2 kg-1 sec-1',                            &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'EPV',                                       &
         LONG_NAME  = 'ertels_potential_vorticity',                &
         UNITS      = 'K m+2 kg-1 sec-1',                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Q',                                         &
         LONG_NAME  = 'specific_humidity',                         &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    do ntracer=1,ntracers
        write(myTracer, "('TRACER_',i1.1)") ntracer
        call MAPL_AddExportSpec ( gc,                             &
            SHORT_NAME = TRIM(myTracer),                              &
            LONG_NAME  = TRIM(myTracer),                              &
            UNITS      = '1',                                         &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
        VERIFY_(STATUS)
    enddo

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DUDTDYN',                                   &
         LONG_NAME  = 'tendency_of_eastward_wind_due_to_dynamics', &
         UNITS      = 'm sec-2',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,            &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DVDTDYN',                                   &
         LONG_NAME  = 'tendency_of_northward_wind_due_to_dynamics',&
         UNITS      = 'm sec-2',                                 &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,            &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                               &
         SHORT_NAME = 'DTDTDYN',                                     &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_dynamics', &
         UNITS      = 'K sec-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'DQVDTDYN',                                      &
         LONG_NAME  = 'tendency_of_specific_humidity_due_to_dynamics', &
         UNITS      = 'kg/kg/sec',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'DQIDTDYN',                                      &
         LONG_NAME  = 'tendency_of_ice_water_due_to_dynamics',         &
         UNITS      = 'kg/kg/sec',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'DQLDTDYN',                                      &
         LONG_NAME  = 'tendency_of_liquid_water_due_to_dynamics',      &
         UNITS      = 'kg/kg/sec',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'DOXDTDYN',                                      &
         LONG_NAME  = 'tendency_of_ozone_due_to_dynamics',             &
         UNITS      = 'kg/kg/sec',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'PREF',                                      &
         LONG_NAME  = 'reference_air_pressure',                    &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'PS',                                  &
       LONG_NAME          = 'surface_pressure',                    &
       UNITS              = 'Pa',                                  &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'TA',                                  &
       LONG_NAME          = 'surface_air_temperature',             &
       UNITS              = 'K',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'QA',                                  &
       LONG_NAME          = 'surface_specific_humidity',           &
       UNITS              = '1',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,        &
       RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'US',                                  &
       LONG_NAME          = 'surface_eastward_wind',               &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,        &
       RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'VS',                                  &
       LONG_NAME          = 'surface_northward_wind',              &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'SPEED',                               &
       LONG_NAME          = 'surface_wind_speed',                  &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'PHIS',                                &
       LONG_NAME          = 'surface_geopotential_height',         &
       UNITS              = 'm',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'DZ',                                  &
       LONG_NAME          = 'surface_layer_height',                &
       UNITS              = 'm',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'SLP',                                 &
       LONG_NAME          = 'sea_level_pressure',                  &
       UNITS              = 'Pa',                                  &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,        &
       RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'TROPP',                               &
       LONG_NAME          = 'tropopause_pressure',                 &
       UNITS              = 'Pa',                                  &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,        &
       RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'TROPT',                               &
       LONG_NAME          = 'tropopause_temperature',              &
       UNITS              = 'K',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
       SHORT_NAME         = 'TROPQ',                               &
       LONG_NAME          = 'tropopause_specific_humidity',        &
       UNITS              = 'kg/kg',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DELP',                                      &
         LONG_NAME  = 'pressure_thickness',                        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U_CGRID',                                   &
         LONG_NAME  = 'eastward_wind_on_C-Grid',                   &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V_CGRID',                                   &
         LONG_NAME  = 'northward_wind_on_C-Grid',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U_DGRID',                                   &
         LONG_NAME  = 'eastward_wind_on_native_D-Grid',            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V_DGRID',                                   &
         LONG_NAME  = 'northward_wind_on_native_D-Grid',           &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'TV',                                        &
         LONG_NAME  = 'air_virtual_temperature',                   &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'THV',                                       &
         LONG_NAME  = 'scaled_virtual_potential_temperature',      &
         UNITS      = 'K/Pa**kappa',                               &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,            &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DDELPDTDYN',                                     &
         LONG_NAME  = 'tendency_of_pressure_thickness_due_to_dynamics', &
         UNITS      = 'Pa sec-1',                                       &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,       &
         RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'UKE',                                       &
         LONG_NAME  = 'eastward_flux_of_atmospheric_kinetic_energy',&
         UNITS      = 'J m-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VKE',                                       &
         LONG_NAME  = 'northward_flux_of_atmospheric_kinetic_energy',          &
         UNITS      = 'J m-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'UCPT',                                      &
         LONG_NAME  = 'eastward_flux_of_atmospheric_enthalpy',     &
         UNITS      = 'J m-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VCPT',                                      &
         LONG_NAME  = 'northward_flux_of_atmospheric_enthalpy',    &
         UNITS      = 'J m-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'UPHI',                                      &
         LONG_NAME  = 'eastward_flux_of_atmospheric_potential_energy',&
         UNITS      = 'J m-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VPHI',                                      &
         LONG_NAME  = 'northward_flux_of_atmospheric_potential_energy',&
         UNITS      = 'J m-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'UQV',                                       &
         LONG_NAME  = 'eastward_flux_of_atmospheric_water_vapor',  &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VQV',                                       &
         LONG_NAME  = 'northward_flux_of_atmospheric_water_vapor', &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'UQL',                                       &
         LONG_NAME  = 'eastward_flux_of_atmospheric_liquid_water', &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VQL',                                       &
         LONG_NAME  = 'northward_flux_of_atmospheric_liquid_water',&
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'UQI',                                       &
         LONG_NAME  = 'eastward_flux_of_atmospheric_ice',          &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VQI',                                       &
         LONG_NAME  = 'northward_flux_of_atmospheric_ice',         &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DKE',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_kinetic_energy_content_due_to_dynamics',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DCPT',                                      &
         LONG_NAME  = 'tendency_of_atmosphere_dry_energy_content_due_to_dynamics',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DPET',                                      &
         LONG_NAME  = 'tendency_of_atmosphere_topographic_potential_energy_due_to_dynamics',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'WRKT',                                      &
         LONG_NAME  = 'work_done_by_atmosphere_at_top',            &
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DQV',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_water_vapor_content_due_to_dynamics',&
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DQL',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_liquid_water_content_due_to_dynamics',&
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'DQI',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_ice_content_due_to_dynamics',&
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'CNV',                                       &
         LONG_NAME  = 'generation_of_atmosphere_kinetic_energy_content',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     do ntracer=1,ntracers
        do nlev=1,nlevs
           write(myTracer, "('TRACER_',i1.1,'_',i3.3)") ntracer, plevs(nlev)
           call MAPL_AddExportSpec ( gc,                             &     
                SHORT_NAME = TRIM(myTracer),                              &
                LONG_NAME  = TRIM(myTracer),                             &
                UNITS      = '1',                                         &
                DIMS       = MAPL_DimsHorzOnly,                           &
                VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
           VERIFY_(STATUS)
        enddo
     enddo         

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VORT850',                                   &
         LONG_NAME  = 'vorticity_at_850_hPa',                      &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VORT700',                                   &
         LONG_NAME  = 'vorticity_at_700_hPa',                      &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'VORT200',                                   &
         LONG_NAME  = 'vorticity_at_200_hPa',                      &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U850',                                      &
         LONG_NAME  = 'eastward_wind_at_850_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &     
         SHORT_NAME = 'U700',                                      &
         LONG_NAME  = 'eastward_wind_at_700_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)         

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U500',                                      &
         LONG_NAME  = 'eastward_wind_at_500_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U250',                                      &
         LONG_NAME  = 'eastward_wind_at_250_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U200',                                      &
         LONG_NAME  = 'eastward_wind_at_200_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V850',                                      &
         LONG_NAME  = 'northward_wind_at_850_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V700',                                      &
         LONG_NAME  = 'northward_wind_at_700_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V500',                                      &
         LONG_NAME  = 'northward_wind_at_500_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V250',                                      &
         LONG_NAME  = 'northward_wind_at_250_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V200',                                      &
         LONG_NAME  = 'northward_wind_at_200_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'T850',                                      &
         LONG_NAME  = 'air_temperature_at_850_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'T700',                                      &
         LONG_NAME  = 'air_temperature_at_700_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'T500',                                      &
         LONG_NAME  = 'air_temperature_at_500_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'T300',                                      &
         LONG_NAME  = 'air_temperature_at_300_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'T250',                                      &
         LONG_NAME  = 'air_temperature_at_250_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Q850',                                      &
         LONG_NAME  = 'specific_humidity_at_850_hPa',              &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Q500',                                      &
         LONG_NAME  = 'specific_humidity_at_500_hPa',              &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Q250',                                      &
         LONG_NAME  = 'specific_humidity_at_250_hPa',              &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Z700',                                      &
         LONG_NAME  = 'geopotential_height_at_700_hPa',            &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Z500',                                      &
         LONG_NAME  = 'geopotential_height_at_500_hPa',            &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'Z300',                                      &
         LONG_NAME  = 'geopotential_height_at_300_hPa',            &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'H850',                                      &
         LONG_NAME  = 'height_at_850_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'H700',                                      &
         LONG_NAME  = 'height_at_700_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'H500',                                      &
         LONG_NAME  = 'height_at_500_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'H300',                                      &
         LONG_NAME  = 'height_at_300_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'H250',                                      &
         LONG_NAME  = 'height_at_250_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'OMEGA850',                                  &
         LONG_NAME  = 'omega_at_850_hPa',                          &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'OMEGA500',                                  &
         LONG_NAME  = 'omega_at_500_hPa',                          &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'U50M',                                      &
         LONG_NAME  = 'eastward_wind_at_50_meters',                &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                             &
         SHORT_NAME = 'V50M',                                      &
         LONG_NAME  = 'northward_wind_at_50_meters',               &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,           &
         RC=STATUS  )
     VERIFY_(STATUS)

! !INTERNAL STATE:

!ALT: technically the first 2 records of "old" style FV restart have 
!     6 ints: YYYY MM DD H M S
!     5 ints: I,J,K, KS (num true pressure levels), NQ (num tracers) headers

    call MAPL_AddInternalSpec ( gc,                           &
         SHORT_NAME = 'AK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_a',                  &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,            &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                           &
         SHORT_NAME = 'BK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_b',                  &
         UNITS      = '1',                                       &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,            &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                           &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                           &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                           &
         SHORT_NAME = 'PT',                                        &
         LONG_NAME  = 'scaled_potential_temperature',              &
         UNITS      = 'K Pa$^{-\kappa}$',                          &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,     &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                           &
         SHORT_NAME = 'PE',                                        &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                           &
         SHORT_NAME = 'PKZ',                                       &
         LONG_NAME  = 'pressure_to_kappa',                         &
         UNITS      = 'Pa$^\kappa$',                                 &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )


!EOP

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="INITIALIZE"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN"         ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN2"        ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-FV_INIT"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--FMS_INIT"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-DYN_CORE"     ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="FINALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( gc, ESMF_SETINIT,  Initialize, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_SETRUN,   Run, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_SETRUN,   RunAddIncs, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_SETFINAL, Finalize, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, "ESMF_ReadRestart", Coldstart, rc=status)
    VERIFY_(STATUS)
 
! Generic SetServices
!--------------------

    call MAPL_GenericSetServices( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine Initialize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: gc       ! composite gridded component 
  type(ESMF_State),    intent(inout) :: import   ! import state
  type(ESMF_State),    intent(inout) :: export   ! export state
  type(ESMF_Clock),    intent(inout) :: clock    ! the clock
  
  integer, intent(out), OPTIONAL     :: rc       ! Error code:
                                                 ! = 0 all is well
                                                 ! otherwise, error
  integer        :: I
  type (ESMF_Grid)                   :: grid
  type (ESMF_Config)                 :: cf
  type (ESMF_Config), pointer        :: config

  type (DYN_wrap)                    :: wrap
  type (DynState),  pointer  :: STATE

  type (MAPL_MetaComp),      pointer :: mapl 

  character (len=ESMF_MAXSTR)        :: restart_file,     &
                                        layout_file


  type (ESMF_Field)                  :: field
  type (ESMF_Array)                  :: array
  type (ESMF_VM)                     :: VM
  real, pointer                      :: pref(:)
  real(r8), pointer                  :: ak(:), bk(:)

  
  integer                            :: status
  character(len=ESMF_MAXSTR)         :: IAm
  character(len=ESMF_MAXSTR)         :: COMP_NAME

  type (ESMF_State)                  :: INTERNAL
  
! Begin
!------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Call Generic Initialize
!------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

! Start the timers
!-----------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"INITIALIZE")

! Get the private internal state
!-------------------------------

    call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, status)
    VERIFY_(STATUS)
    state => wrap%dyn_state

! Get file names from the configuration
!--------------------------------------

!BOR
! !RESOURCE_ITEM: none :: name of layout file
    call MAPL_GetResource ( MAPL, layout_file, 'LAYOUT:', default='fvcore_layout.rc', rc=status )
!EOR
    VERIFY_(STATUS)

! Set Private Internal State from Restart File
! --------------------------------------------

    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"-FV_INIT")
    call DynInit ( STATE, CLOCK, INTERNAL, IMPORT, GC, LAYOUT_FILE)
    call MAPL_TimerOff(MAPL,"-FV_INIT")

! Create PREF EXPORT Coupling (Needs to be done only once per run)
! ----------------------------------------------------------------

    call MAPL_GetPointer(EXPORT,PREF,'PREF',ALLOC=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, AK, 'AK', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, BK, 'BK', RC=STATUS)
    VERIFY_(STATUS)

    PREF = AK + BK * P00

    call MAPL_TimerOff(MAPL,"INITIALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!BOP

! !IROUTINE: Run

! !DESCRIPTION: This is the first Run stage of FV. It is the container
!    for the dycore calculations. Subroutines from the core are 
!    invoked to do most of the work. A second run method, descibed below,
!    adds the import tendencies from external sources to the FV 
!    variables. 
!
!    In addition to computing and adding all dynamical contributions
!    to the FV variables (i.e., winds, pressures, and temperatures),
!    this method advects an arbitrary number of  tracers. These appear
!    in a ``Friendly'' bundle in the IMPORT state and are updated with
!    the advective tendency.
!
!
! !INTERFACE:

subroutine Run(gc, import, export, clock, rc)

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: gc
  type (ESMF_State),   intent(inout) :: import
  type (ESMF_State),   intent(inout) :: export
  type (ESMF_Clock),   intent(in)    :: clock
  integer, intent(out), optional     :: rc 

!EOP

! !Local Variables:
  
    integer                                          :: status
    type (ESMF_Bundle)                               :: bundle
    type (ESMF_Field)                                :: field
    type (ESMF_Array)                                :: array
    type (ESMF_Config)                               :: cf

    type (MAPL_MetaComp), pointer :: genstate 

    type (DYN_wrap) :: wrap
    type (DynState), pointer :: STATE
    type (DynGrid),  pointer :: GRID
    type (DynVars),  pointer :: VARS
    
    integer  :: J1, JN, K1, KN, NQ
    integer  :: IM, JM, KM
    integer  :: NKE, NCPT, NPHI, NP
    integer  :: ifirstxy, ilastxy, jfirstxy, jlastxy
    integer  :: I, J, K, n
    logical, parameter :: convt = .false. ! Until this is run with full physics
    logical  :: is_ringing

    real(r8),     pointer :: phisxy(:,:)
    real(kind=4), pointer ::   phis(:,:)

    real(r8), allocatable ::   pkxy(:,:,:) ! pe**kappa
    real(r8), allocatable ::     pl(:,:,:) ! mid-level pressure
    real(r8), allocatable :: tempxy(:,:,:) ! mid-level temperature
    real(r8), allocatable ::     ua(:,:,:) ! temporary array
    real(r8), allocatable ::     va(:,:,:) ! temporary array
    real(r8), allocatable ::     qv(:,:,:) ! temporary array
    real(r8), allocatable ::    zle(:,:,:) ! temporary array
    real(r8), allocatable ::    pke(:,:,:) ! temporary array
    real(r8), allocatable ::   delp(:,:,:) ! temporary array
    real(r8), allocatable ::   dudt(:,:,:) ! temporary array
    real(r8), allocatable ::   dvdt(:,:,:) ! temporary array
    real(r8), allocatable ::   dtdt(:,:,:) ! temporary array
    real(r8), allocatable ::   dqdt(:,:,:) ! temporary array
    real(r8), allocatable ::  ddpdt(:,:,:) ! temporary array
    real(r8), allocatable ::   dmdt(:,:)   ! temporary array
    real(r8), allocatable :: tmp2d (:,:)   ! temporary array
    real(r8), allocatable :: tmp3d (:,:,:) ! temporary array

    real(r8), allocatable, target :: ke    (:,:,:) ! Kinetic    Energy
    real(r8), allocatable, target :: cpt   (:,:,:) ! Internal   Energy
    real(r8), allocatable, target :: phi   (:,:,:) ! Potential  Energy
    real(r8), allocatable, target :: pkinv (:,:,:) ! 1.0 / Mid-Level P**kappa
    real(r8), allocatable :: qke   (:,:)   ! Vertically Integrated Kinetic   Energy Tracer
    real(r8), allocatable :: qcpt  (:,:)   ! Vertically Integrated Internal  Energy Tracer
    real(r8), allocatable :: qphi  (:,:)   ! Vertically Integrated Potential Energy Tracer

    real(r8), allocatable :: phi00 (:,:)   ! Vertically Integrated phi
    real(r8), allocatable :: penrg (:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrg (:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrg (:,:)   ! PHIS*(Psurf-Ptop)
    real(r8), allocatable :: penrg0(:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrg0(:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrg0(:,:)   ! PHIS*(Psurf-Ptop)
    real(r8), allocatable :: penrga(:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrga(:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrga(:,:)   ! PHIS*(Psurf-Ptop)
    real(r8), allocatable :: penrgb(:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrgb(:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrgb(:,:)   ! PHIS*(Psurf-Ptop)

    real(r8), allocatable :: pelnxz(:,:,:) ! log pressure (pe) at layer edges
    real(r8), allocatable :: omgxyz(:,:,:) ! vertical pressure velocity (pa/sec)
    real(r8), allocatable :: epvxyz(:,:,:) ! ertel's potential vorticity
    real(r8), allocatable :: cxxyz(:,:,:)  ! Accumulated zonal winds
    real(r8), allocatable :: cyxyz(:,:,:)  ! Accumulated meridional winds
    real(r8), allocatable :: mfxxyz_ur(:,:,:) ! zonal mass flux
    real(r8), allocatable :: mfyxyz_ur(:,:,:) ! meridional mass flux
    real(r8), allocatable :: mfxxyz(:,:,:) ! zonal mass flux
    real(r8), allocatable :: mfyxyz(:,:,:) ! meridional mass flux
    real(r8), allocatable :: mfzxyz(:,:,:) ! vertical mass flux
    real(r8), allocatable :: mfxxyz_a(:,:,:) ! zonal mass flux A-Grid
    real(r8), allocatable :: mfyxyz_a(:,:,:) ! meridional mass flux A-Grid
    real(r8)              :: dt            ! Dynamics time step
    real(r8)              :: alpha
    real(r8)              :: kinetic       ! local kinetic   energy
    real(r8)              :: potential     ! local potential energy

    real(kind=4), pointer :: dqldt (:,:,:)
    real(kind=4), pointer :: dqidt (:,:,:)
    real(kind=4), pointer :: doxdt (:,:,:)
    real(kind=4), pointer :: qtend (:,:,:)
    real(kind=4), pointer :: vort  (:,:,:)
    real(kind=4), pointer :: temp3d(:,:,:)
    real(kind=4), pointer :: temp2d(:,:)
    real(kind=4), pointer :: tempu (:,:)
    real(kind=4), pointer :: tempv (:,:)

    character(len=ESMF_MAXSTR), ALLOCATABLE       :: NAMES (:)
    character(len=ESMF_MAXSTR), ALLOCATABLE, save :: NAMES0(:)
    character(len=ESMF_MAXSTR) :: IAm
    character(len=ESMF_MAXSTR) :: COMP_NAME
    character(len=ESMF_MAXSTR) :: STRING

    type(DynTracers)            :: qqq       ! Specific Humidity

  Iam = "Run"
  call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
  VERIFY_(STATUS)
  Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the generic state
! -----------------------------------------

  call MAPL_GetObjectFromGC (GC, GENSTATE,  RC=STATUS )
  VERIFY_(STATUS)

  call MAPL_TimerOn(GENSTATE,"TOTAL")
  call MAPL_TimerOn(GENSTATE,"RUN")

! Retrieve the pointer to the internal state
! ------------------------------------------


  call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, status)
  VERIFY_(STATUS)
  state => wrap%dyn_state

  vars  => state%vars   ! direct handle to control variables
  grid  => state%grid   ! direct handle to grid
  dt    =  state%dt     ! dynamics time step (large)

  ifirstxy = grid%is
  ilastxy  = grid%ie
  jfirstxy = grid%js
  jlastxy  = grid%je

  im       = grid%npx
  jm       = grid%npy
  km       = grid%npz

  is_ringing = ESMF_AlarmIsRinging( STATE%ALARMS(TIME_TO_RUN),rc=status); VERIFY_(status) 
  if (.not. is_ringing) return


! Allocate Arrays
! ---------------
      ALLOCATE(   delp(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dudt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dvdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dtdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dqdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  ddpdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE( tempxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     pl(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     ua(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     va(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     qv(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )

      ALLOCATE(     ke(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(    cpt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(    phi(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  pkinv(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )

      ALLOCATE(   qke (ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(   qcpt(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(   qphi(ifirstxy:ilastxy,jfirstxy:jlastxy)    )

      ALLOCATE(   dmdt(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(  phi00(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(  kenrg(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(  penrg(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(  tenrg(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( kenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( penrg0(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( tenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy)    )

      ALLOCATE( kenrga(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( penrga(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( tenrga(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( kenrgb(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( penrgb(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( tenrgb(ifirstxy:ilastxy,jfirstxy:jlastxy)    )

      ALLOCATE(  tmp3d(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  tmp2d(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( phisxy(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(   pkxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE(    zle(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE(    pke(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE( pelnxz(ifirstxy:ilastxy,km+1,jfirstxy:jlastxy) )
      ALLOCATE( omgxyz(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( epvxyz(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE(  cxxyz(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE(  cyxyz(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfxxyz_ur(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfyxyz_ur(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfxxyz(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfyxyz(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfzxyz(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE( mfxxyz_a(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfyxyz_a(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )


! Surface Geopotential from IMPORT state
!---------------------------------------

      call MAPL_GetPointer ( IMPORT, PHIS, 'PHIS', RC=STATUS )
      VERIFY_(STATUS)

      phisxy = real(phis,kind=r8)


! Energetics before analysis increments of u, v, t, p, AND Q.
!------------------------------------------------------------

      call getAgridWinds(vars%u, vars%v, ua, va, rotate=.true.)

      delp   = vars%pe(:,:,2:)  -vars%pe(:,:,:km)   ! Pressure Thickness

      call MAPL_GetPointer ( IMPORT, qtend, 'DQVANA', RC=STATUS )
      VERIFY_(STATUS)
      call PULL_Q ( STATE, IMPORT, qqq, NXQ, rc )
      tempxy = vars%pt
      if( associated(qtend) ) then
         if (associated(qqq%content_r4)) then
            if ( size(qqq%content_r4) == size(tempxy) ) qqq%content_r4 = max(qqq%content_r4,0._r4)
            if ( size(qqq%content_r4) == size(tempxy) ) tempxy = tempxy * (1.0+eps*(qqq%content_r4-qtend))
         elseif (associated(qqq%content)) then
            if ( size(qqq%content) == size(tempxy) ) qqq%content    = max(qqq%content   ,0._r8)
            if ( size(qqq%content) == size(tempxy) ) tempxy = tempxy * (1.0+eps*(qqq%content   -qtend))
         endif
      endif

      call Energetics (ua,va,tempxy,vars%pe,delp,vars%pkz,phisxy,kenrg,penrg,tenrg,ke,cpt,phi,phi00)

! Add Diabatic Forcing from Analysis to State Variables
! -----------------------------------------------------

      call ADD_INCS ( STATE,IMPORT,DT )

! Clear mass fluxes
!------------------

      mfxxyz_ur  (:,:,:) = 0.
      mfyxyz_ur  (:,:,:) = 0.
      mfxxyz  (:,:,:) = 0.
      mfyxyz  (:,:,:) = 0.
      mfzxyz  (:,:,:) = 0.
      mfxxyz_a(:,:,:) = 0.
      mfyxyz_a(:,:,:) = 0.


! Report advected friendlies
!---------------------------

      call ESMF_StateGetBundle ( IMPORT, 'QTR' , BUNDLE,   RC=STATUS )
      VERIFY_(STATUS)
      call ESMF_BundleGet ( BUNDLE, fieldCount=NQ, RC=STATUS )
      VERIFY_(STATUS)

    if (NQ > 0) then
      allocate( NAMES(NQ),STAT=STATUS )
      VERIFY_(STATUS)
      call ESMF_BundleGetFieldNames ( BUNDLE, nameList=NAMES, rc=STATUS )
      VERIFY_(STATUS)

      if( .not.allocated( names0 ) ) then
           allocate( NAMES0(NQ),STAT=STATUS )
           VERIFY_(STATUS)
           write(STRING,'(A,I3,A)') "Advecting the following ", nq, " tracers in DYN:"
              call WRITE_PARALLEL( trim(STRING)   )
           do k=1,nq
              call WRITE_PARALLEL( trim(NAMES(K)) )
           end do
           NAMES0 = NAMES
      endif

      if( size(names0).ne.size(names) ) then
           deallocate( NAMES0 )
             allocate( NAMES0(NQ),STAT=STATUS )
           VERIFY_(STATUS)
           write(STRING,'(A,I3,A)') "Advecting the following ", nq, " tracers in DYN:"
              call WRITE_PARALLEL( trim(STRING)   )
           do k=1,nq
              call WRITE_PARALLEL( trim(NAMES(K)) )
           end do
           NAMES0 = NAMES
      endif
    endif

! Compute Dry Temperature
!------------------------

      tempxy = vars%pt * vars%pkz

! Get A-grid winds
! ----------------

      call getAgridWinds(vars%u, vars%v, ua, va, rotate=.true.)

! Load Tracer Friendlies
! ----------------------

      delp   = vars%pe(:,:,2:)  -vars%pe(:,:,:km)   ! Pressure Thickness
      dmdt   = vars%pe(:,:,km+1)-vars%pe(:,:,1)

! Get specific humidity
!----------------------

      if ( (qqq%is_r4) .and. associated(qqq%content_r4) ) then
         if (size(     qv)==size(qqq%content_r4)) then
            qv = qqq%content_r4
            ASSERT_(all(qv >= 0.0))
            vars%pt = vars%pt * (1.0+eps*qqq%content_r4)
            dqdt    = qqq%content_r4 * delp
         endif
      elseif (associated(qqq%content)) then
         if (size(     qv)==size(qqq%content)) then
            qv = qqq%content
            ASSERT_(all(qv >= 0.0))
            vars%pt = vars%pt * (1.0+eps*qqq%content   )
            dqdt    = qqq%content    * delp
         endif
      endif

! Initialize Diagnostic Dynamics Tendencies
! -----------------------------------------

      ddpdt  =          delp       ! Pressure Thickness Tendency
      dudt   =     ua * delp       ! U-Wind on A-Grid   Tendency
      dvdt   =     va * delp       ! V-Wind on A-Grid   Tendency
      dtdt   = tempxy * delp       ! Dry Temperature    Tendency

! Initialize Tracer Dynamics Tendencies
! -------------------------------------

      call MAPL_GetPointer( export,dqldt,'DQLDTDYN', rc=status )
      VERIFY_(STATUS)
      call MAPL_GetPointer( export,dqidt,'DQIDTDYN', rc=status )
      VERIFY_(STATUS)
      call MAPL_GetPointer( export,doxdt,'DOXDTDYN', rc=status )
      VERIFY_(STATUS)

     if (allocated(names)) then
      if( associated(dqldt) ) then
          dqldt = 0.0
          do k = 1,size(names)
             if( trim(names(k)).eq.'QLCN' .or. &
                 trim(names(k)).eq.'QLLS' ) then
                 if( state%vars%tracer(k)%is_r4 ) then
                     if (size(dqldt)==size(state%vars%tracer(k)%content_r4)) &
                              dqldt = dqldt + state%vars%tracer(k)%content_r4
                 else
                     if (size(dqldt)==size(state%vars%tracer(k)%content)) &
                              dqldt = dqldt + state%vars%tracer(k)%content
                 endif
             endif
          enddo
          dqldt = dqldt * delp
      endif

      if( associated(dqidt) ) then
          dqidt = 0.0
          do k = 1,size(names)
             if( trim(names(k)).eq.'QICN' .or. &
                 trim(names(k)).eq.'QILS' ) then    
                 if( state%vars%tracer(k)%is_r4 ) then 
                     if (size(dqidt)==size(state%vars%tracer(k)%content_r4)) &
                              dqidt = dqidt + state%vars%tracer(k)%content_r4
                 else
                     if (size(dqidt)==size(state%vars%tracer(k)%content)) &
                              dqidt = dqidt + state%vars%tracer(k)%content
                 endif
             endif
          enddo
          dqidt = dqidt * delp
      endif

      if( associated(doxdt) ) then
          doxdt = 0.0
          do k = 1,size(names)
             if( trim(names(k)).eq.'OX' ) then
                 if( state%vars%tracer(k)%is_r4 ) then 
                     if (size(doxdt)==size(state%vars%tracer(k)%content_r4)) &
                              doxdt = doxdt + state%vars%tracer(k)%content_r4
                 else
                     if (size(doxdt)==size(state%vars%tracer(k)%content)) &
                              doxdt = doxdt + state%vars%tracer(k)%content
                 endif
             endif
          enddo
          doxdt = doxdt * delp
      endif
    endif

! Compute Energetics Before Dycore
! --------------------------------

      call Energetics (ua,va,vars%pt,vars%pe,delp,vars%pkz,phisxy, kenrg0,penrg0,tenrg0,ke,cpt,phi,phi00)

    kenrg = (kenrg0-kenrg)/DT
    penrg = (penrg0-penrg)/DT
    tenrg = (tenrg0-tenrg)/DT

    call FILLOUT2 (export, 'KEANA', kenrg, rc=status); VERIFY_(STATUS)
    call FILLOUT2 (export, 'PEANA', penrg, rc=status); VERIFY_(STATUS)
    call FILLOUT2 (export, 'TEANA', tenrg, rc=status); VERIFY_(STATUS)

! Add Passive Tracers for KE, CPT, and PHI
! ----------------------------------------
    nq = STATE%GRID%NQ
    if (NXQ == 4) then
      NKE  = nq-3
      NCPT = nq-2
      NPHI = nq-1
      NP   = nq

      state%vars%tracer(NKE )%content => KE
      state%vars%tracer(NCPT)%content => CPT
      state%vars%tracer(NPHI)%content => PHI

      state%vars%tracer(NKE )%is_r4 = .false.
      state%vars%tracer(NCPT)%is_r4 = .false.
      state%vars%tracer(NPHI)%is_r4 = .false.

      pkinv = 1.0 / vars%pkz

      state%vars%tracer(NP)%content => pkinv
      state%vars%tracer(NP)%is_r4   = .false.
    else
      NKE  = -1
      NCPT = -1
      NPHI = -1
      NP   = -1
    endif
 
      call MAPL_TimerOn(GENSTATE,"-DYN_CORE")
#if defined(USE_FV)
      call FVdycore_wrapper( phisxy, tempxy, qqq, STATE,    &
           pkxy, pelnxz, omgxyz, epvxyz,  &
           cxxyz, cyxyz,                  &
           mfxxyz_ur, mfyxyz_ur,          &
           mfxxyz, mfyxyz, mfzxyz, convt, &
           kenrga, penrga, tenrga,        &
           kenrgb, penrgb, tenrgb,        &
           PI, CP, KAPPA, OMEGA,          &
           RADIUS, GRAV, RGAS, EPS,       &
           NKE, NCPT, NPHI, NP            )
#else
      call DynRun(STATE, CLOCK)

#if defined(DEBUG)         
  call Write_Profile(grid, vars%u, 'U-after-DynRun')
  call Write_Profile(grid, vars%v, 'V-after-DynRun')
#endif

      call getMassFluxes(mfxxyz, mfyxyz, mfzxyz, dt)
#endif
      call MAPL_TimerOff(GENSTATE,"-DYN_CORE")

      call getPK ( pkxy )

    if (SW_DYNAMICS) then

      call MAPL_GetPointer(export,temp2d,'PHIS', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = phisxy

      call MAPL_GetPointer(export,temp2d,'PS',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d =  vars%pe(:,:,km+1)/GRAV

      call getAgridWinds(vars%u, vars%v, ua, va, rotate=.true.)
      call FILLOUT3 (export, 'U'      , ua      , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V'      , va      , rc=status); VERIFY_(STATUS)

      call FILLOUT3 (export, 'U_DGRID', vars%u  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V_DGRID', vars%v  , rc=status); VERIFY_(STATUS)

      call MAPL_GetPointer(export, vort, 'VORT', rc=status)
      VERIFY_(STATUS)
      if(associated(vort)) then
          call getVorticity(vars%u, vars%v, tmp3d)
          vort=tmp3d
      endif
 
    else               ! .not. SW_DYNAMICS
! Load Local Variable with Vapor Specific Humidity
! ------------------------------------------------


    if (STATE%GRID%NQ > 0) then
      if ( qqq%is_r4 ) then
         if (size(qv)==size(qqq%content_r4)) qv = qqq%content_r4
      else
         if (size(qv)==size(qqq%content)) qv = qqq%content
      endif
    else
      qv = 0.0
    endif

! Compute Dry Theta and T with Unified Poles
! ------------------------------------------

      vars%pt = vars%pt / (1.0+eps*qv )
      tempxy  = vars%pt * vars%pkz

! Compute Mid-Layer Pressure and Pressure Thickness
! -------------------------------------------------

      delp = ( vars%pe(:,:,2:) - vars%pe(:,:,:km) )
      pl   = ( vars%pe(:,:,2:) + vars%pe(:,:,:km) ) * 0.5

! Compute Tropopause Pressure, Temperature, and Moisture
! ------------------------------------------------------

      alpha = 0.03  ! Empirical Value for Nominal Conditions

      call tropopause (alpha,ilastxy-ifirstxy+1,jlastxy-jfirstxy+1,km, &
                       vars%pe,pl,tempxy,qv,       &
                       zle(:,:,1),zle(:,:,2),zle(:,:,3) )

      call FILLOUT2 (export, 'TROPP',zle(:,:,1), rc=status); VERIFY_(STATUS)
      call FILLOUT2 (export, 'TROPT',zle(:,:,2), rc=status); VERIFY_(STATUS)
      call FILLOUT2 (export, 'TROPQ',zle(:,:,3), rc=status); VERIFY_(STATUS)

! Compute A-Grid Winds
! --------------------

      call getAgridWinds(vars%u, vars%v, ua, va, rotate=.true.)

#if defined(LATLON_CORE)
! Compute A-Grid Mass Fluxes
! --------------------------
      call c2a3d( grid, mfxxyz, mfyxyz, mfxxyz_a, mfyxyz_a )
#endif

! Compute Diagnostic Dynamics Tendencies
!  (Note: initial values of d(m,u,v,T,q)/dt are progs m,u,v,T,q)
! --------------------------------------------------------------

      dmdt = ( vars%pe(:,:,km+1)-vars%pe(:,:,1) - dmdt )/(grav*dt)

      dudt = (    ua*delp  -dudt )/(delp*dt)
      dvdt = (    va*delp  -dvdt )/(delp*dt)
      dtdt = (tempxy*delp - dtdt )/(delp*dt)
      dqdt = (    qv*delp - dqdt )/(delp*dt)

      ddpdt = ( delp - ddpdt )/dt ! Pressure Thickness Tendency


      call FILLOUT3 (export, 'DELP'      ,delp , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'DUDTDYN'   ,dudt , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'DVDTDYN'   ,dvdt , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'DTDTDYN'   ,dtdt , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'DQVDTDYN'  ,dqdt , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'DDELPDTDYN',ddpdt, rc=status); VERIFY_(STATUS)
 
      call FILLOUT3 (export, 'U_CGRID'   ,cxxyz, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V_CGRID'   ,cyxyz, rc=status); VERIFY_(STATUS)

      call FILLOUT3 (export, 'MFX_UR' , mfxxyz_ur  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFY_UR' , mfyxyz_ur  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFX'    , mfxxyz  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFY'    , mfyxyz  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFZ'    , mfzxyz  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFX_A'  , mfxxyz_a, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFY_A'  , mfyxyz_a, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'U'      , ua      , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V'      , va      , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'T'      , tempxy  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'Q'      , qv      , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'PL'     , pl      , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'PLE'    , vars%pe , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'PLK'    , vars%pkz, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'U_DGRID', vars%u  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V_DGRID', vars%v  , rc=status); VERIFY_(STATUS)

      call MAPL_GetPointer(export, vort, 'VORT', rc=status)
      VERIFY_(STATUS)
      if(associated(vort)) then
          call getVorticity(vars%u, vars%v, tmp3d)
          vort=tmp3d
      endif

      call MAPL_GetPointer(export, temp3D, 'EPV', rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = epvxyz*(p00**kappa)

      do ntracer=1,ntracers
         write(myTracer, "('TRACER_',i1.1)") ntracer
         call MAPL_GetPointer(export, temp3D, TRIM(myTracer), rc=status)
         VERIFY_(STATUS)
         if((associated(temp3d)) .and. (NQ>=ntracer+1)) then
            if (state%vars%tracer(ntracer+1)%is_r4) then
               temp3d = state%vars%tracer(ntracer+1)%content_r4
            else
               temp3d = state%vars%tracer(ntracer+1)%content
            endif
         endif
      enddo
 
      call MAPL_GetPointer(export, temp3D, 'PV', rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = epvxyz/vars%pt
 
      call MAPL_GetPointer(export, temp3D, 'S', rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = tempxy*cp

      call MAPL_GetPointer(export, temp3d, 'TH',rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = vars%pt*(p00**kappa)

      call MAPL_GetPointer(export, temp2d, 'DMDTDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = dmdt

! Compute Tracer Dynamics Tendencies
! ----------------------------------
    
    if (allocated(names)) then
      if( associated(dqldt) ) then
          dqldt = -dqldt
          do k = 1,size(names)
             if( trim(names(k)).eq.'QLCN' .or. &
                 trim(names(k)).eq.'QLLS' ) then
                 if( state%vars%tracer(k)%is_r4 ) then 
                     if (size(dqldt)==size(state%vars%tracer(k)%content_r4)) &
                              dqldt = dqldt + state%vars%tracer(k)%content_r4 * delp
                 else
                     if (size(dqldt)==size(state%vars%tracer(k)%content)) &
                              dqldt = dqldt + state%vars%tracer(k)%content    * delp
                 endif
             endif
          enddo
          dqldt = dqldt / (delp*dt)
      endif

      if( associated(dqidt) ) then
          dqidt = -dqidt
          do k = 1,size(names)
             if( trim(names(k)).eq.'QICN' .or. &
                 trim(names(k)).eq.'QILS' ) then
                 if( state%vars%tracer(k)%is_r4 ) then 
                     if (size(dqidt)==size(state%vars%tracer(k)%content_r4)) &
                              dqidt = dqidt + state%vars%tracer(k)%content_r4 * delp
                 else
                     if (size(dqidt)==size(state%vars%tracer(k)%content)) &
                              dqidt = dqidt + state%vars%tracer(k)%content    * delp
                 endif
             endif
          enddo
          dqidt = dqidt / (delp*dt)
      endif

      if( associated(doxdt) ) then
          doxdt = -doxdt
          do k = 1,size(names)
             if( trim(names(k)).eq.'OX' ) then
                 if( state%vars%tracer(k)%is_r4 ) then
                     if (size(doxdt)==size(state%vars%tracer(k)%content_r4)) &
                              doxdt = doxdt + state%vars%tracer(k)%content_r4 * delp
                 else
                     if (size(doxdt)==size(state%vars%tracer(k)%content)) &
                              doxdt = doxdt + state%vars%tracer(k)%content    * delp
                 endif
             endif
          enddo
          doxdt = doxdt / (delp*dt)
      endif

! Compute Vertically Integrated Tracer Dynamics Tendencies
! --------------------------------------------------------

      call MAPL_GetPointer ( export, temp2D, 'DQVDTINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          temp2d = 0.0
          do k=1,km
             temp2d = temp2d + dqdt(:,:,k)*delp(:,:,k)
          enddo
             temp2d = temp2d * (1.0/grav)
      endif

      call MAPL_GetPointer ( export, temp2D, 'DQLDTINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D)    ) then
          if( associated(dqldt) ) then
          temp2d = 0.0
          do k=1,km
             temp2d = temp2d + dqldt(:,:,k)*delp(:,:,k)
          enddo
             temp2d = temp2d * (1.0/grav)
          else
          temp2d = MAPL_UNDEF
          endif
      endif

      call MAPL_GetPointer ( export, temp2D, 'DQIDTINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D)    ) then
          if( associated(dqidt) ) then
          temp2d = 0.0
          do k=1,km
             temp2d = temp2d + dqidt(:,:,k)*delp(:,:,k)
          enddo
             temp2d = temp2d * (1.0/grav)
          else
          temp2d = MAPL_UNDEF
          endif
      endif

      call MAPL_GetPointer ( export, temp2D, 'DOXDTINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D)    ) then
          if( associated(doxdt) ) then
          temp2d = 0.0
          do k=1,km
             temp2d = temp2d + doxdt(:,:,k)*delp(:,:,k)
          enddo
             temp2d = temp2d * (1.0/grav)
          else
          temp2d = MAPL_UNDEF
          endif
      endif
    endif

! Fill Surface and Near-Surface Variables
! ---------------------------------------

      call MAPL_GetPointer(export,temp2d,'PS',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d =  vars%pe(:,:,km+1)    

      call MAPL_GetPointer(export,temp2d,'US',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d =       ua(:,:,km)      

      call MAPL_GetPointer(export,temp2d,'VS'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = va(:,:,km)      

      call MAPL_GetPointer(export,temp2d,'TA'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d =   tempxy(:,:,km)

      call MAPL_GetPointer(export,temp2d,'QA'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d =       qv(:,:,km)

      call MAPL_GetPointer(export,temp2d,'SPEED',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = sqrt( ua(:,:,km)**2 + va(:,:,km)**2 )


! Virtual temperature
! -------------------

      tempxy =  tempxy*(1.0+eps*qv)

      call MAPL_GetPointer(export,temp3D,'TV'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp3D)) temp3D = tempxy


! Funny Tracers
! -------------

    if (STATE%GRID%NQ > 0) then
      call MAPL_GetPointer(export,temp2d,'CONVKE',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         qke = 0.0
         do k=1,km
            qke = qke + ke(:,:,k)*delp(:,:,k)
         enddo
         temp2d = (qke/grav-kenrg0)/dt
      end if

      call MAPL_GetPointer(export,temp2d,'CONVCPT',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         qcpt = 0.0
         do k=1,km
            qcpt = qcpt + cpt(:,:,k)*delp(:,:,k)
         enddo
         temp2d = (qcpt/grav-penrg0)/dt
      end if

      call MAPL_GetPointer(export,temp2d,'CONVPHI',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         qphi = 0.0
         do k=1,km
            qphi = qphi + phi(:,:,k)*delp(:,:,k)
         enddo
!        temp2d = ( qphi/grav-(tenrg0+(rgas/cp)*penrg0) )/dt
         temp2d = ( qphi/grav-phi00 )/dt
      end if
     endif

! Fluxes
!-------

      call MAPL_GetPointer(export,temp2d,'UCPT',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
         do k=1,km
            temp2d = temp2d + ua(:,:,k)*tempxy(:,:,k)*delp(:,:,k)
         enddo
         temp2d = temp2d*(cp/grav)
      end if

      call MAPL_GetPointer(export,temp2d,'VCPT',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
         do k=1,km
            temp2d = temp2d + va(:,:,k)*tempxy(:,:,k)*delp(:,:,k)
         enddo
         temp2d = temp2d*(cp/grav)
      end if

      call MAPL_GetPointer(export,tempu,'UKE',rc=status); VERIFY_(STATUS)
      call MAPL_GetPointer(export,tempv,'VKE',rc=status); VERIFY_(STATUS)

      if(associated(tempu) .or. associated(tempv)) then
         ke = 0.5*(ua**2 + va**2)
      end if

      if(associated(tempu)) then
         tempu = 0.0
         do k=1,km
            tempu = tempu + ua(:,:,k)*ke(:,:,k)*delp(:,:,k)
         enddo
         tempu = tempu / grav
      end if

      if(associated(tempv)) then
         tempv = 0.0
         do k=1,km
            tempv = tempv + va(:,:,k)*ke(:,:,k)*delp(:,:,k)
         enddo
         tempv = tempv / grav
      end if

      call MAPL_GetPointer(export,temp2d,'UQV',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
         do k=1,km
            temp2d = temp2d + ua(:,:,k)*QV(:,:,k)*delp(:,:,k)
         enddo
         temp2d = temp2d / grav
      end if

      call MAPL_GetPointer(export,temp2d,'VQV',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
         do k=1,km
            temp2d = temp2d + va(:,:,k)*QV(:,:,k)*delp(:,:,k)
         enddo
         temp2d = temp2d / grav
      end if

    if (allocated(names)) then
      call MAPL_GetPointer(export,temp2d,'UQL',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QLCN' .or. &
                 trim(names(n)).eq.'QLLS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then 
                     if (size(ua)==size(state%vars%tracer(n)%content_r4)) &
                              temp2d = temp2d + ua(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                     if (size(ua)==size(state%vars%tracer(n)%content)) &
                              temp2d = temp2d + ua(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo

             endif
          enddo
         temp2d = temp2d / grav
      end if

      call MAPL_GetPointer(export,temp2d,'VQL',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QLCN' .or. &
                 trim(names(n)).eq.'QLLS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then 
                     if (size(ua)==size(state%vars%tracer(n)%content_r4)) &
                              temp2d = temp2d + va(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                     if (size(ua)==size(state%vars%tracer(n)%content)) &
                              temp2d = temp2d + va(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo

             endif
          enddo
         temp2d = temp2d / grav
      end if

      call MAPL_GetPointer(export,temp2d,'UQI',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QICN' .or. &
                 trim(names(n)).eq.'QILS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then 
                     if (size(ua)==size(state%vars%tracer(n)%content_r4)) &
                              temp2d = temp2d + ua(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                     if (size(ua)==size(state%vars%tracer(n)%content)) &
                              temp2d = temp2d + ua(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo

             endif
          enddo
         temp2d = temp2d / grav
      end if

      call MAPL_GetPointer(export,temp2d,'VQI',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QICN' .or. &
                 trim(names(n)).eq.'QILS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then
                     if (size(va)==size(state%vars%tracer(n)%content_r4)) &
                              temp2d = temp2d + va(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                     if (size(va)==size(state%vars%tracer(n)%content)) &
                              temp2d = temp2d + va(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo

             endif
          enddo
         temp2d = temp2d / grav
      end if
    endif

! Compute Energetics After Dycore
! -------------------------------

      tempxy = vars%pt*(1.0+eps*qv)  ! Convert TH to THV

      call Energetics (ua,va,tempxy,vars%pe,delp,vars%pkz,phisxy,kenrg,penrg,tenrg,ke,cpt,phi,phi00)

      call MAPL_GetPointer(export,temp2d,'KEDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = (kenrg-kenrg0)/DT

      call MAPL_GetPointer(export,temp2d,'PEDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = (penrg-penrg0)/DT

      call MAPL_GetPointer(export,temp2d,'TEDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = (tenrg-tenrg0)/DT

       kenrg = (kenrga-kenrg0)/DT
       penrg = (penrga-penrg0)/DT
       tenrg = (tenrga-tenrg0)/DT

      call FILLOUT2 (export,'KECDCOR', kenrg)
      call FILLOUT2 (export,'PECDCOR', penrg)
      call FILLOUT2 (export,'TECDCOR', tenrg)

       kenrg = (kenrgb-kenrga)/DT
       penrg = (penrgb-penrga)/DT
       tenrg = (tenrgb-tenrga)/DT

      call FILLOUT2 (export,'KETEMAP', kenrg)
      call FILLOUT2 (export,'PETEMAP', penrg)
      call FILLOUT2 (export,'TETEMAP', tenrg)


! Height related diagnostics
! --------------------------
      zle(:,:,km+1) = phisxy(:,:)
      do k=km,1,-1
        zle(:,:,k) = zle(:,:,k+1) + cp*tempxy(:,:,k)*( pkxy(:,:,k+1)-pkxy(:,:,k) )
      enddo
      zle = zle/grav

      call MAPL_GetPointer(export,temp3d,'THV',rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = tempxy

      call MAPL_GetPointer(export,temp3d,'ZLE',rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = zle

      call MAPL_GetPointer(export,temp2d,'PHIS', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = zle(:,:,km+1)*grav

      call MAPL_GetPointer(export,temp2d,'DZ', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = 0.5*( zle(:,:,km)-zle(:,:,km+1) )

      call MAPL_GetPointer(export,temp3d,'ZL' ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = 0.5*( zle(:,:,:km)+zle(:,:,2:) )

      call MAPL_GetPointer(export,temp3d,'S'  ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = temp3d + grav*(0.5*( zle(:,:,:km)+zle(:,:,2:) ))

! Compute Kinetic Energy Generation and Omega
! -------------------------------------------

#if defined(USE_FV)
      call getOmega ( omgxyz, delp, zle )
#else
      call getOmega ( omgxyz )
#endif

      call MAPL_GetPointer(export,temp2d,'KEGEN', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         qke = 0.0
      do k=1,km
         qke = qke - omgxyz(:,:,k)*delp(:,:,k)
      enddo
         temp2d = qke / grav
      end if

      call FILLOUT3 (export, 'OMEGA'  , omgxyz     , rc=status)
      VERIFY_(STATUS)

      zle = log(vars%pe)

      call MAPL_GetPointer(export,temp2d,'OMEGA850', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(temp2d,omgxyz,zle,log(85000.)  , status)
         VERIFY_(STATUS)
      end if

      call MAPL_GetPointer(export,temp2d,'OMEGA500', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(temp2d,omgxyz,zle,log(50000.)  , status)
         VERIFY_(STATUS)
      end if

     end if   ! SW_DYNAMICS
      
  
 
! De-Allocate Arrays
! ------------------

      DEALLOCATE( KENRG  )
      DEALLOCATE( PENRG  )
      DEALLOCATE( TENRG  )
      DEALLOCATE( KENRG0 )
      DEALLOCATE( PENRG0 )
      DEALLOCATE( TENRG0 )

      DEALLOCATE( KENRGA )
      DEALLOCATE( PENRGA )
      DEALLOCATE( TENRGA )
      DEALLOCATE( KENRGB )
      DEALLOCATE( PENRGB )
      DEALLOCATE( TENRGB )

      DEALLOCATE( qke  , ke  )
      DEALLOCATE( qcpt , cpt )
      DEALLOCATE( qphi , phi )
      DEALLOCATE( pkinv )

      DEALLOCATE( pke    )
      DEALLOCATE( ZLE    )
      DEALLOCATE( PKXY   )
      DEALLOCATE( tmp3d  )
      DEALLOCATE( tmp2d  )
      DEALLOCATE( pelnxz )
      DEALLOCATE( omgxyz )
      DEALLOCATE( epvxyz )
      DEALLOCATE(  cxxyz )
      DEALLOCATE(  cyxyz )
      DEALLOCATE( mfxxyz_ur )
      DEALLOCATE( mfyxyz_ur )
      DEALLOCATE( mfxxyz )
      DEALLOCATE( mfyxyz )
      DEALLOCATE( mfzxyz )
      DEALLOCATE( mfxxyz_a )
      DEALLOCATE( mfyxyz_a )
      DEALLOCATE( tempxy )
      DEALLOCATE( pl     )
      DEALLOCATE( va     )
      DEALLOCATE( ua     )
      DEALLOCATE( qv     )
      DEALLOCATE( delp   )
      DEALLOCATE( dmdt   )
      DEALLOCATE( dudt   )
      DEALLOCATE( dvdt   )
      DEALLOCATE( dtdt   )
      DEALLOCATE( dqdt   )
      DEALLOCATE( ddpdt  )
      DEALLOCATE( phisxy )

      if (allocated(names)) DEALLOCATE( names  )

      DEALLOCATE( phi00  )
  
      call freeTracers(state)

  call MAPL_TimerOff(GENSTATE,"RUN")
  call MAPL_TimerOff(GENSTATE,"TOTAL")

  RETURN_(ESMF_SUCCESS)

end subroutine RUN

!-----------------------------------------------------------------------

  subroutine PULL_Q(STATE, IMPORT, QQQ, iNXQ, RC)

    type (DynState)        :: STATE
    type (ESMF_State)              :: IMPORT
    type (DynTracers)               :: QQQ       ! Specific Humidity
    integer,           intent(IN)  :: iNXQ
    integer, optional, intent(OUT) :: RC

    integer                          :: STATUS
    character(len=ESMF_MAXSTR)       :: IAm="Pull_Q"
    character(len=ESMF_MAXSTR)       :: FIELDNAME
    type (ESMF_Bundle)               :: BUNDLE
    type (ESMF_Field)                :: field
    type (ESMF_Array)                :: array
    type (ESMF_DataKind)             :: kind
    real(r4),              pointer   :: ptr_r4(:,:,:), humidity(:,:,:)
    real(r8),              pointer   :: ptr_r8(:,:,:)
    integer                          :: I,K,N,NQ
    logical                          :: EMPTY
    integer                          :: i1,in,j1,jn,im,jm,km
    real(r8)                         :: sumout

    i1 = state%grid%is
    in = state%grid%ie
    j1 = state%grid%js
    jn = state%grid%je
    im = state%grid%npx
    jm = state%grid%npy
    km = state%grid%npz

    call ESMF_StateGetBundle(IMPORT, 'QTR' , BUNDLE,   RC=STATUS)
    VERIFY_(STATUS)

! Count the friendlies
!---------------------

    call ESMF_BundleGet(BUNDLE, fieldCount=NQ, RC=STATUS)
    VERIFY_(STATUS)

               NQ = NQ + iNXQ
    STATE%GRID%NQ = NQ       ! GRID%NQ is now the "official" NQ

!
! Tracer pointer array
!
    ALLOCATE(STATE%VARS%tracer(nq), STAT=STATUS)
    VERIFY_(STATUS)

    DO n = 1, NQ-iNXQ
       call ESMF_BundleGetField(bundle, fieldIndex=n, field=field, rc=status)
       VERIFY_(STATUS)
       call ESMF_FieldGet(field, name=fieldname, rc=status)
       VERIFY_(STATUS)
       call ESMF_FieldGetArray(field, array, rc=status)
       VERIFY_(STATUS)

       call ESMF_ArrayGet(array,kind=kind,rc=status)
       VERIFY_(STATUS)

       STATE%VARS%TRACER(N)%IS_R4  = (kind == ESMF_R4)   ! Is real*4?

       if ( STATE%VARS%TRACER(N)%IS_R4 ) then

          call ESMF_ArrayGetData(array, ptr_r4, ESMF_DATA_REF, rc=status)
          VERIFY_(STATUS)
          state%vars%tracer(n)%content_r4 => MAPL_RemapBounds(PTR_R4, i1,in,j1,jn, &
                                                              1, km)
          if (fieldname == "Q") then
             qqq%is_r4 = .true.
             qqq%content_r4 => state%vars%tracer(n)%content_r4
          end if

       else

          if (fieldname == "Q") then
             print*, "Q is assumed to real4, I don't know why"
          end if

       endif
       END DO

  end subroutine PULL_Q

!-----------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!BOP

! !IROUTINE: RunAddIncs

! !DESCRIPTION: This is the second registered stage of FV. 
!    It calls an Fv supplied routine to add external contributions 
!    to FV's state variables. It does not touch the Friendly tracers.
!    It also computes additional diagnostics and updates the 
!    FV internal state to reflect the added tendencies.
!
!
! !INTERFACE:

  subroutine RunAddIncs(gc, import, export, clock, rc)

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc
    type (ESMF_State),   intent(inout) :: import
    type (ESMF_State),   intent(inout) :: export
    type (ESMF_Clock),   intent(in)    :: clock
    integer, intent(out), optional     :: rc 

!EOP

! !Local Variables:
  
    integer                                          :: status
    character(len=ESMF_MAXSTR) :: IAm

    type (MAPL_MetaComp), pointer :: genstate 

    type (DYN_wrap) :: wrap
    type (DynState), pointer :: STATE
    type (DynGrid),  pointer :: GRID
    type (DynVars),  pointer :: VARS
    type (DynTracers)                 :: qqq     ! Specific Humidity
    
    real(r8), allocatable :: penrg (:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrg (:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrg (:,:)   ! PHIS*(Psurf-Ptop)
    real(r8), allocatable :: penrg0(:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrg0(:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrg0(:,:)   ! PHIS*(Psurf-Ptop)

    real(r8),     pointer :: phisxy(:,:)
    real(r4),     pointer ::   phis(:,:)
    real(r8), allocatable ::    slp(:,:)
    real(r8), allocatable ::  tmp2d(:,:)
    real(r8), allocatable ::  tmp3d(:,:,:)
    real(r8), allocatable ::    pke(:,:,:)
    real(r8), allocatable ::     pl(:,:,:)
    real(r8), allocatable ::     ua(:,:,:)
    real(r8), allocatable ::     va(:,:,:)
    real(r8), allocatable ::     qv(:,:,:)
    real(r8), allocatable ::     dp(:,:,:)
    real(r8), allocatable ::    thv(:,:,:)
    real(r8), allocatable ::    zle(:,:,:)
    real(r8), allocatable :: tempxy(:,:,:)

    real(r8)              :: dt
    real(r8)              :: delp          ! delta pressure thickness
    real(r8)              :: kinetic       ! local kinetic   energy
    real(r8)              :: potential     ! local potential energy


    real(r4), pointer     :: QOLD(:,:,:)
    real(r4), pointer     :: temp3d(:,:,:)
    real(r4), pointer     :: temp2d(:,:  )

    integer ifirstxy, ilastxy
    integer jfirstxy, jlastxy
    integer im,jm,km, inxq
    integer i,j,k

    character(len=ESMF_MAXSTR) :: COMP_NAME

    Iam = "RunAddIncs"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the generic state
! -----------------------------------------

    call MAPL_GetObjectFromGC (GC, GENSTATE,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(GENSTATE,"TOTAL")
    call MAPL_TimerOn(GENSTATE,"RUN2")

! Retrieve the pointer to the internal state
! ------------------------------------------

    call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, status)
    VERIFY_(STATUS)
    state => wrap%dyn_state

    vars  => state%vars   ! direct handle to control variables
    grid  => state%grid   ! direct handle to grid
    dt    =  state%dt     ! dynamics time step (large)

    ifirstxy = grid%is
    ilastxy  = grid%ie
    jfirstxy = grid%js
    jlastxy  = grid%je

    im  = grid%npx
    jm  = grid%npy
    km  = grid%npz
    inxq = 0

    if (.not. SW_DYNAMICS) then

    ALLOCATE(  kenrg(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE(  penrg(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE(  tenrg(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE( kenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE( penrg0(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE( tenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy) )

    ALLOCATE(  tmp3d(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
    ALLOCATE(  tmp2d(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE( phisxy(ifirstxy:ilastxy,jfirstxy:jlastxy) )

    ALLOCATE(     ua(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE(     va(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE(     qv(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE(     pl(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE(     dp(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE(    thv(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE( tempxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )

    ALLOCATE(    pke(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
    ALLOCATE(    zle(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )


    call MAPL_GetPointer ( IMPORT, PHIS, 'PHIS', RC=STATUS )
    VERIFY_(STATUS)

    phisxy = real(phis,kind=r8)

! Compute Pressure Thickness
! --------------------------

    dp = ( vars%pe(:,:,2:) - vars%pe (:,:,:km) )

! Get A-grid winds
! ----------------
  call getAgridWinds(vars%u, vars%v, ua, va, rotate=.true.)

! Load Specific Humidity
! ----------------------


    call MAPL_GetPointer(export,QOLD,'Q',  rc=status)

    call PULL_Q ( STATE, IMPORT, qqq, inxq, rc )
    if (STATE%GRID%NQ > 0) then
      if ( (qqq%is_r4) .and. (associated(qqq%content_r4)) ) then
       if (size(qv)==size(qqq%content_r4)) qv = qqq%content_r4
      elseif (associated(qqq%content)) then
       if (size(qv)==size(qqq%content)) qv = qqq%content
      endif
    else
      qv = 0.0
    endif

! Compute Energetics Before Diabatic Forcing
! ------------------------------------------
    if (associated(QOLD)) then
       thv = vars%pt*(1.0+eps*QOLD)
    else
       thv = vars%pt
    endif

    call Energetics (ua,va,thv,vars%pe,dp,vars%pkz,phisxy,kenrg0,penrg0,tenrg0,tempxy,tempxy,tempxy,tmp2d)


! Add Diabatic Forcing to State Variables
! ---------------------------------------
    call ADD_INCS ( STATE,IMPORT,DT )

! Update Mid-Layer Pressure and Pressure Thickness
! ------------------------------------------------

    dp = ( vars%pe(:,:,2:) - vars%pe (:,:,:km) )
    pl = ( vars%pe(:,:,2:) + vars%pe (:,:,:km) )*0.5

! Create A-Grid Winds
! -------------------
    call getAgridWinds(vars%u, vars%v, ua, va, rotate=.true.)

! Compute Energetics After Diabatic Forcing
! -----------------------------------------

    thv = vars%pt*(1.0+eps*qv)

    call Energetics (ua,va,thv,vars%pe,dp,vars%pkz,phisxy,kenrg,penrg,tenrg,tempxy,tempxy,tempxy,tmp2d)

    kenrg = (kenrg-kenrg0)/DT
    penrg = (penrg-penrg0)/DT
    tenrg = (tenrg-tenrg0)/DT

    call FILLOUT2 (export, 'KEPHY', kenrg, rc=status); VERIFY_(STATUS)
    call FILLOUT2 (export, 'PEPHY', penrg, rc=status); VERIFY_(STATUS)
    call FILLOUT2 (export, 'TEPHY', tenrg, rc=status); VERIFY_(STATUS)

! Fill Diagnostics
! ----------------

    tempxy = vars%pt * vars%pkz   ! Dry Temperature

    call FILLOUT3 (export, 'DELP'   , dp      , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'U'      , ua      , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'V'      , va      , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'T'      , tempxy  , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'Q'      , qv      , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PL'     , pl      , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PLE'    , vars%pe , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PLK'    , vars%pkz, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'THV'    , thv     , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'U_DGRID', vars%u  , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'V_DGRID', vars%v  , rc=status); VERIFY_(STATUS)

    call MAPL_GetPointer(export,temp3d,'TH',rc=status)
    VERIFY_(STATUS)
    if(associated(temp3d)) temp3d = vars%pt*(p00**kappa)

      do ntracer=1,ntracers
         write(myTracer, "('TRACER_',i1.1)") ntracer
         call MAPL_GetPointer(export, temp3D, TRIM(myTracer), rc=status)
         VERIFY_(STATUS)
         if((associated(temp3d)) .and. (STATE%GRID%NQ>=ntracer+1)) then
            if (state%vars%tracer(ntracer+1)%is_r4) then
               temp3d = state%vars%tracer(ntracer+1)%content_r4
            else
               temp3d = state%vars%tracer(ntracer+1)%content
            endif
         endif
      enddo

! Compute Edge Heights
! --------------------

    pke           = vars%pe**kappa

    zle(:,:,km+1) = phisxy(:,:)/grav
    do k=km,1,-1
       zle(:,:,k) = zle(:,:,k+1) + (cp/grav)*thv(:,:,k)*( pke(:,:,k+1)-pke(:,:,k) )
    enddo
    call FILLOUT3 (export, 'ZLE', zle, rc=status); VERIFY_(STATUS)

! Compute Mid-Layer Heights
! -------------------------

    call MAPL_GetPointer(export,temp3d,'ZL',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp3d)) temp3d =  0.5*( zle(:,:,2:) + zle(:,:,:km) )

! Fill Single Level Variables
! ---------------------------

    pke = log(vars%pe)

! Write Pressure Level Tracer Data
! --------------------------------
    do ntracer=1,ntracers
       do nlev=1,nlevs
          write(myTracer, "('TRACER_',i1.1,'_',i3.3)") ntracer, plevs(nlev)
          call MAPL_GetPointer(export, temp2d, TRIM(myTracer), rc=status)
          VERIFY_(STATUS)
          if((associated(temp2d)) .and. (STATE%GRID%NQ>=ntracer+1)) then
             if (state%vars%tracer(ntracer+1)%is_r4) then
               tmp3d = state%vars%tracer(ntracer+1)%content_r4(:,:,:)
             else
               tmp3d = state%vars%tracer(ntracer+1)%content(:,:,:)
             endif
             call VertInterp(temp2d,tmp3d,pke,log(plevs(nlev)*100.0)  ,  status)
             VERIFY_(STATUS)
          endif
       enddo
    enddo

    call MAPL_GetPointer(export,temp2d,'VORT200',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call getVorticity(vars%u, vars%v, tmp3d)
       call VertInterp(temp2d,tmp3d,pke,log(20000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'VORT700',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call getVorticity(vars%u, vars%v, tmp3d)
       call VertInterp(temp2d,tmp3d,pke,log(70000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'VORT850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call getVorticity(vars%u, vars%v, tmp3d)
       call VertInterp(temp2d,tmp3d,pke,log(85000.)  ,  status)
       VERIFY_(STATUS)
    end if
 
    call MAPL_GetPointer(export,temp2d,'U200',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,ua,pke,log(20000.)  ,  status)
       VERIFY_(STATUS)
    end if
 
    call MAPL_GetPointer(export,temp2d,'U250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,ua,pke,log(25000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'U500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,ua,pke,log(50000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'U700',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,ua,pke,log(70000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'U850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,ua,pke,log(85000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V200',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,va,pke,log(20000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,va,pke,log(25000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,va,pke,log(50000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V700',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,va,pke,log(70000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,va,pke,log(85000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,tempxy,pke,log(25000.)  , status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T300',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,tempxy,pke,log(30000.)  , status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,tempxy,pke,log(50000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T700',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,tempxy,pke,log(70000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,tempxy,pke,log(85000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Q250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,qv,pke,log(25000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Q500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,qv,pke,log(50000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Q850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,qv,pke,log(85000.)  ,  status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Z700',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle*grav,pke,log(70000.)  , status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Z500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle*grav,pke,log(50000.)  , status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Z300',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle*grav,pke,log(30000.)  , status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle,pke,log(25000.)  , status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H300',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle,pke,log(30000.)  , status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle,pke,log(50000.)  , status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H700',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle,pke,log(70000.)  , status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle,pke,log(85000.)  , status)
       VERIFY_(STATUS)
    end if

! Compute Heights Above Surface
! -----------------------------
    do k=1,km+1
    zle(:,:,k) = zle(:,:,k) - zle(:,:,km+1)
    enddo
    
    call MAPL_GetPointer(export,temp2d,'U50M',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,ua,-zle,-50., status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V50M',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,va,-zle,-50., status)
       VERIFY_(STATUS)
    end if

! Compute Surface Pressure
! ------------------------

    call MAPL_GetPointer(export,temp2d,'PS',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) temp2d=vars%pe(:,:,km+1)

! Compute Vertically Averaged T,U
! -------------------------------
    call MAPL_GetPointer(export,temp2d,'TAVE',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       temp2d = 0.0
       do k=1,km
       temp2d = temp2d + tempxy(:,:,k)*dp(:,:,k)
       enddo
       temp2d = temp2d / (vars%pe(:,:,km+1)-vars%pe(:,:,1))
    endif

    call MAPL_GetPointer(export,temp2d,'UAVE',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       temp2d = 0.0
       do k=1,km
       temp2d = temp2d + ua(:,:,k)*dp(:,:,k)
       enddo
       temp2d = temp2d / (vars%pe(:,:,km+1)-vars%pe(:,:,1))
    endif

! Convert T to Tv
! ---------------

    tempxy = tempxy*(1.0+eps*qv)

    call MAPL_GetPointer(export,temp3d,'TV',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp3d)) temp3d=tempxy

! Compute Sea-Level Pressure
! --------------------------
    
    call MAPL_GetPointer(export,temp2d,'SLP',rc=status)
    VERIFY_(STATUS)

    if(associated(temp2d)) then
       ALLOCATE( slp(ifirstxy:ilastxy,jfirstxy:jlastxy) )
       do j=jfirstxy,jlastxy
          do i=ifirstxy,ilastxy
             call get_slp ( km,vars%pe (i,j,  km+1),phisxy(i,j),  slp(i,j), &
                               vars%pe (i,j,1:km+1),                        &
                               vars%pkz(i,j,1:km  ),                        &
                                 tempxy(i,j,1:km  ) )
          enddo
       enddo
       temp2d =   slp
       DEALLOCATE(slp)
    end if

! Deallocate Memory
! -----------------

    DEALLOCATE(  kenrg )
    DEALLOCATE(  penrg )
    DEALLOCATE(  tenrg )
    DEALLOCATE( kenrg0 )
    DEALLOCATE( penrg0 )
    DEALLOCATE( tenrg0 )
    DEALLOCATE(  tmp2d )
    DEALLOCATE(  tmp3d )

    DEALLOCATE( phisxy )

    DEALLOCATE(     ua )
    DEALLOCATE(     va )
    DEALLOCATE(     qv )
    DEALLOCATE(     pl )
    DEALLOCATE(     dp )
    DEALLOCATE( tempxy )

    DEALLOCATE(    thv )
    DEALLOCATE(    pke )
    DEALLOCATE(    zle )

    end if ! .not. SW_DYNAMICS

    call MAPL_TimerOff(GENSTATE,"RUN2")
    call MAPL_TimerOff(GENSTATE,"TOTAL")

    RETURN_(ESMF_SUCCESS)
end subroutine RunAddIncs

!-----------------------------------------------------------------------
  subroutine ADD_INCS ( STATE,IMPORT,DT,QOLD,QNEW,RC )

!
! !INPUT PARAMETERS:

   type(DynState), intent(INOUT)  :: STATE
   type(ESMF_State),       intent(INOUT)  :: IMPORT
   real(r8),               intent(IN   )  :: DT
   real(r4), optional,     intent(IN   )  :: QOLD (:,:,:)
   real(r8), optional,     intent(IN   )  :: QNEW (:,:,:)
   integer,  optional,     intent(OUT  )  :: RC

!
! !DESCRIPTION:  This routine adds the tendencies to the state,
!                weighted appropriately by the time step.  Temperature
!                tendencies are pressure weighted (ie., DELP*DT/Dt).
!                All tendencies are on the A-grid, and have an XY decomposition.
!

    integer                          :: status

    integer               :: I1, IN, J1, JN, K, im, jm, km
    integer               :: KL, KU
    real(r8)              :: SUMOUT
    real(r8), allocatable ::    dum(:,:,:)
    real(r8), allocatable ::    pke(:,:,:),  dpinv(:,:,:)
    real(r8), allocatable :: tend_ua(:,:,:), tend_va(:,:,:)
    real(r8), allocatable :: tend_un(:,:,:), tend_vn(:,:,:)
    real(kind=4), pointer :: tend(:,:,:)

    character(len=ESMF_MAXSTR)         :: IAm="ADD_INCS"


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
       dpinv(:,:,k) = 1.0/( state%vars%pe(:,:,k+1)-state%vars%pe(:,:,k) )
    enddo

!#define ADIABATIC
#if !defined(ADIABATIC)
! **********************************************************************
! ****                      Wind Tendencies                         ****
! ****         Note: State Variables are on the D-Grid,             ****
! ****        while IMPORT Tendencies are on the A-Grid             ****
! **********************************************************************

    ALLOCATE( tend_ua(state%grid%isd:state%grid%ied  ,state%grid%jsd:state%grid%jed  ,km) )
    ALLOCATE( tend_va(state%grid%isd:state%grid%ied  ,state%grid%jsd:state%grid%jed  ,km) )
    ALLOCATE( tend_un(state%grid%isd:state%grid%ied  ,state%grid%jsd:state%grid%jed+1,km) )
    ALLOCATE( tend_vn(state%grid%isd:state%grid%ied+1,state%grid%jsd:state%grid%jed  ,km) )

    call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DUDT',RC=STATUS )
    VERIFY_(STATUS)
    tend_ua(i1:in,j1:jn,1:km) = tend

    call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DVDT',RC=STATUS )
    VERIFY_(STATUS)
    tend_va(i1:in,j1:jn,1:km) = tend

! Put the wind tendencies on the Native Dynamics grid
! ---------------------------------------------------
    call Agrid_To_Native( state%grid, tend_ua, tend_va, tend_un, tend_vn )

! Add the wind tendencies to the control variables
! ------------------------------------------------
    STATE%VARS%U = STATE%VARS%U + DT*TEND_UN(i1:in,j1:jn,1:km)
    STATE%VARS%V = STATE%VARS%V + DT*TEND_VN(i1:in,j1:jn,1:km)

    DEALLOCATE( tend_ua )
    DEALLOCATE( tend_va )
    DEALLOCATE( tend_un )
    DEALLOCATE( tend_vn )

! **********************************************************************
! ****                     Pressure Tendency                        ****
! **********************************************************************

    call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DPEDT',RC=STATUS )
    VERIFY_(STATUS)

    KL = lbound( tend,3 )
    KU = ubound( tend,3 )
    allocate( dum(i1:in,j1:jn,KL:KU) )

    DUM = DT*TEND
    STATE%VARS%PE = STATE%VARS%PE + DUM
    DEALLOCATE (DUM)

! **********************************************************************
! ****            Update Diagnostic State Variable PKZ              ****
! **********************************************************************
  
    ALLOCATE( pke(i1:in,j1:jn,km+1) )

    pke = STATE%VARS%PE**kappa
    do k=1,km
    STATE%VARS%PKZ(:,:,k) = ( pke(:,:,k+1)-pke(:,:,k) ) &
                          / ( kappa*log( STATE%VARS%PE(:,:,k+1)/STATE%VARS%PE(:,:,k) ) )
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
    allocate( dum(i1:in,j1:jn,KL:KU) )

    DUM = DT*TEND*DPINV/STATE%VARS%PKZ

    STATE%VARS%PT =  STATE%VARS%PT              + DUM
!   STATE%VARS%PT = (STATE%VARS%PT*(1+EPS*QOLD) + DUM)/(1+EPS*QNEW)
    DEALLOCATE (DUM)
#endif

    DEALLOCATE( DPINV )

   return

 end subroutine ADD_INCS


   subroutine FILLOUT3(export, name, V, RC)
     type (ESMF_State),  intent(inout) :: export
     character(len=*),   intent(IN   ) :: name
     real(r8),           intent(IN   ) :: V(:,:,:)
     integer, optional,  intent(  out) :: rc

     real(r4), pointer          :: CPL(:,:,:)
     integer                    :: status
     character(len=ESMF_MAXSTR) :: IAm="Fillout3"

     call MAPL_GetPointer(export, cpl, name, RC=STATUS)
     VERIFY_(STATUS)
     if(associated(cpl)) cpl=v

   end subroutine FILLOUT3

!-----------------------------------------------------------------------

   subroutine FILLOUT2(export, name, V, rc)
     type (ESMF_State),  intent(inout) :: export
     character(len=*),   intent(IN   ) :: name
     real(r8),           intent(IN   ) :: V(:,:)
     integer, optional,  intent(  out) :: rc

     real(kind=4), pointer      :: CPL(:,:)
     integer                    :: status
     character(len=ESMF_MAXSTR) :: IAm="Fillout2"

     call MAPL_GetPointer(export, cpl, name, RC=STATUS)
     VERIFY_(STATUS)
     if(associated(cpl)) cpl=v

     return
   end subroutine FILLOUT2

!-----------------------------------------------------------------------

  subroutine Energetics (ua,va,thv,ple,delp,pk,phis,ke,pe,te,tke,cpt,phi,phii)

  real(8)  cpt(:,:,:)
  real(8)  phi(:,:,:)
  real(8)  tke(:,:,:)
  real(8)   ua(:,:,:)
  real(8)   va(:,:,:)
  real(8)  thv(:,:,:)
  real(8)  ple(:,:,:)
  real(8) delp(:,:,:)
  real(8)   pk(:,:,:)
  real(8)   ke(:,:)
  real(8)   pe(:,:)
  real(8)   te(:,:)
  real(8) phis(:,:)
  real(8) phii(:,:)

  real(8) kinetic, potential
  integer i,ifirst,ilast
  integer j,jfirst,jlast
  integer k,km

  real(8), allocatable ::   pke(:,:,:)
  real(8), allocatable :: gztop(:,:)

  ifirst = lbound( ua,1 ) 
  ilast  = ubound( ua,1 ) 
  jfirst = lbound( ua,2 ) 
  jlast  = ubound( ua,2 ) 
  km     = ubound( ua,3 ) 

  allocate( gztop( ifirst:ilast, jfirst:jlast          ) )
  allocate( pke  ( ifirst:ilast, jfirst:jlast , 1:km+1 ) )

! Compute Model Top Height
! ------------------------
      pke   = ple**kappa
      gztop = phis
      phi(:,:,km) = gztop
      do k=km,1,-1
               gztop(:,:)     = gztop(:,:) + cp*thv(:,:,k)*( pke(:,:,k+1)-pke(:,:,k) )
                 phi(:,:,k)   = ( phi(:,:,k)+gztop )*0.5
      if(k.ne.1) phi(:,:,k-1) = gztop
      enddo

! Compute Energetics:  Cp*Tv + K + PHI
! ------------------------------------
       ke = 0.0
       pe = 0.0
     phii = 0.0
  do k=1,km
  do j=jfirst,jlast
  do i=ifirst,ilast
       kinetic   = 0.5*( ua(i,j,k)**2 + va(i,j,k)**2 )
       potential = cp*thv(i,j,k)*pk(i,j,k)
      tke(i,j,k) = kinetic
      cpt(i,j,k) = potential
       ke(i,j)   =   ke(i,j) +    kinetic*delp(i,j,k)
       pe(i,j)   =   pe(i,j) +  potential*delp(i,j,k)
     phii(i,j)   = phii(i,j) + phi(i,j,k)*delp(i,j,k)
  enddo
  enddo
  enddo
       ke(:,:) =    ke(:,:)/grav
       pe(:,:) =    pe(:,:)/grav
     phii(:,:) =  phii(:,:)/grav
       te(:,:) = (phis(:,:)*ple(:,:,km+1)-gztop(:,:)*ple(:,:,1))/grav

  deallocate ( gztop )
  deallocate ( pke   )

  return
  end subroutine Energetics


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Finalize

! !DESCRIPTION: Writes restarts and cleans-up through MAPL\_GenericFinalize and
!   deallocates memory from the Private Internal state. 
!
! !INTERFACE:

subroutine Finalize(gc, import, export, clock, rc)

! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: gc
    type (ESMF_State),    intent(inout) :: import
    type (ESMF_State),    intent(inout) :: export
    type (ESMF_Clock),    intent(inout) :: clock
    integer, optional,    intent(  out) :: rc
 
!EOP

! Local variables
    type (DYN_wrap) :: wrap
    type (DynState), pointer  :: STATE
    character (len=ESMF_MAXSTR)       :: restart_file
 
    character(len=ESMF_MAXSTR)        :: IAm
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    integer                           :: status

    type (MAPL_MetaComp),     pointer :: MAPL 
    type (ESMF_Config)                :: cf


! BEGIN

    Iam = "Finalize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, config=cf, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"FINALIZE")

! Retrieve the pointer to the state
!----------------------------------

    call ESMF_UserCompGetInternalState(gc, 'DYNstate', wrap, status)
    VERIFY_(STATUS)
  
    state => wrap%dyn_state
 
    call DynFinalize( STATE )
 
! Call Generic Finalize
!----------------------

    call MAPL_TimerOff(MAPL,"FINALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL")

    call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

  RETURN_(ESMF_SUCCESS)
 
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine PRINT_TIMES(TIMES,DAYS)
  integer(kind=8), intent(INOUT) :: TIMES(:,:)
  real(r8),        intent(IN   ) :: DAYS
  TIMES = 0
 
  return
 end subroutine PRINT_TIMES
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine FINALIZE

      subroutine get_slp ( km,ps,phis,slp,pe,pk,tv)
      implicit   none
      integer  km
      real(r8)   pk(km)    ! layer-mean P**kappa
      real(r8)   tv(km)    ! layer-mean virtual Temperature
      real(r8)   pe(km+1)  ! press at layer edges (Pa)
      real(r8)   ps        ! surface pressure (Pa)
      real(r8) phis        ! surface geopotential
      real(r8)  slp        ! sea-level pressure (hPa)
      real(r8)  tstar                 ! extrapolated temperature (K)
      real(r8) p_bot
      real(r8) tref                   ! Reference virtual temperature (K)
      real(r8) pref                   ! Reference pressure level (Pa)
      real(r8) pkref                  ! Reference pressure level (Pa) ** kappa
      real(r8) dp1, dp2

      real(r8), parameter :: gamma    = 6.5e-3
      real(r8), parameter :: p_offset = 15000.
      real(r8), parameter :: gg       = gamma/MAPL_GRAV

      real(r8), parameter :: factor   = MAPL_grav / ( MAPL_Rgas * gamma ) 
      real(r8), parameter :: yfactor  = MAPL_Rgas * gg

      integer k_bot, k, k1, k2

      p_bot = ps - p_offset
      k_bot = -1

      do k = km, 2, -1
         if ( pe(k+1) .lt. p_bot ) then
              k_bot = k
              exit
         endif
      enddo

      k1    = k_bot - 1
      k2    = k_bot
      dp1   = pe(k_bot)   - pe(k_bot-1)
      dp2   = pe(k_bot+1) - pe(k_bot)
      pkref = ( pk(k1)*dp1 + pk(k2)*dp2 ) / (dp1+dp2)
       tref = ( tv(k1)*dp1 + tv(k2)*dp2 ) / (dp1+dp2)
       pref = 0.5 * ( pe(k_bot+1) + pe(k_bot-1) )
      tstar = tref*( ps/pref )**yfactor

      slp   = ps*( 1.0+gg*phis/tstar )**factor

      return
  end subroutine get_slp

      subroutine tropopause (alpha,im,jm,lm,ple,pl,      &
                             tmpu,sphu,TropP,TropT,TropQ )

!********************************************************************
!!                     Subroutine tropopause
!********************************************************************
!
!!IROUTINE:  tropopause
!
!!DESCRIPTION:
!
!     This routine finds the tropopause pressure and temperature from
!     model temperature profiles.  It finds these values for a
!     given 3-dimensional grid.  The algorithm is based on the
!     temperature profile only.  It is similar to visually inspecting
!     a Skew-T Log P diagram for the left-most kink in the temperature
!     profile.  In a standard Skew-T Log P diagram the isotherms
!     intersect the Log P scale at a 45 degree angle.  To more
!     effectively isolate the tropopause,  this angle can be adjusted.
!     That is how this program works.  It adjusts the angle at which
!     the isotherms intersect the Log P axis via the coeffiecient
!     alpha. A It simply looks for the furthest point to the left on
!     the profile.  The routine defines the leftward position of each
!     temperature point as xfact:
!
!         xfact = (alpha * temperature) - log (pressure)
!
!     The tropopause for a given profile is where xfact is a minimum.
!     Uncertainty can occur when the kink is not very distinct.  This
!     is where the selection of alpha becomes important.  Optimal
!     values for alpha appear to be between .02 and .04.  For lower
!     values of alpha, the tropopause selection will favor lower
!     levels (higher P).  For higher values of alpha, the selection
!     will favor higher levels (lower P).  A value of alpha = .03
!     appears to be optimal in generating tropopause values within
!     the range indicated in the Handbook of Geophysics and the Space
!     Environment, AFGL, chapt 14, 1985.
!
!!INPUT PARAMETERS:
!
!     alpha   = see discussion above ((log mb)/deg K)
!     im      = number of longitude grid points
!     jm      = number of latitude grid points
!     lm      = number of model levels
!     ple     = pressures at model edges      (Pa)
!     pl      = pressures at model mid-layers (Pa)
!     TmpU    = 3-d array of gridded temperature on model level (K)
!     SphU    = 3-d array of gridded spec.hum.   on model level (kg/kg)
!
!!OUTPUT PARAMETERS:
!
!     TropP = array of gridded tropopause pressures (Pa)
!     TropT = array of gridded tropopause temperatures (K)
!     TropQ = array of gridded tropopause specific humidity (kg/kg)
!
!!REVISION HISTORY:
!
!     Created 25 Jun 97 by Jim Stobie
!
!********************************************************************

      implicit none

! Passed Variables

      integer im
      integer jm
      integer lm

      real(r8) alpha
      real(r8) pl   (im,jm,lm)    ! Pressure at model mid-layers
      real(r8) ple  (im,jm,lm+1)  ! Pressure at model edges
      real(r8) TmpU (im,jm,lm)
      real(r8) Sphu (im,jm,lm)
      real(r8) TropP(im,jm)
      real(r8) TropT(im,jm)
      real(r8) TropQ(im,jm)

! Local Variables

      integer i,j,k          !loop variables
      integer kend           !end index for model level search
      integer kstart         !start index for model level search
      integer ktrop          !index for tropopause level

      real(r8) phigh         !highest pressure for search
      real(r8) plow          !lowest pressure for search
      real(r8) undef         !value for undefined variables
      real(r8) xfacmn        !minimum x-factor, see prologue
      real(r8) xfact(lm)     !x-factor, see prologue

      undef = MAPL_UNDEF

!----------------------------------------------------------------
! Set vertical limits on search.  Tropopause search will be
! limited to the range between plow and phigh (mb).
! According to Handbook of Geophysics and the Space Environment,
! AFGL, 1985, pg 14-6, the lowest tropopause values are near 8 km
! in the polar winter (approx 350 mb) and the highest near 18 km
! in the tropics (approx 80 mb).
!----------------------------------------------------------------

        plow  =  4000.0
        phigh = 55000.0

!--------------------
! Loop over lat/lon Grid
!--------------------

        do j = 1, jm
        do i = 1, im

!-------------------------------------------------------
! Find pressure range for search.  Search will begin
! at first model level edge above phigh, and end
! at the first model level edge below plow.
!-------------------------------------------------------

               kend = lm
            do while (ple(i,j,kend).GE.phigh)
               kend = kend-1
            enddo

               kstart = 1
            do while (ple(i,j,kstart).le.plow)
               kstart = kstart+1
            enddo

!-----------------------------------------------------
! Calculate pressure of the model layer midpoints.
! Then calculate xfact for these points.  See prologue
! for description of xfact.
!-----------------------------------------------------

            do k = kstart, kend
              xfact(k) = alpha * TmpU(i,j,k) - log10(pl(i,j,k)*0.01)
            end do

!-----------------------------------------------
! Tropopause is level for which xfact is minimum
!-----------------------------------------------

            xfacmn = 100000.
            ktrop  = 0

            do k = kstart, kend
               if (xfact(k).LT.xfacmn) then
                   xfacmn = xfact(k)
                   ktrop  = k
               endif
            enddo

!-------------------------------------------------------
! If the minimum value of xfact is at the upper or lower
! boundary of the search, then the tropopause has not
! been sucessfully isolated.  In this case a warning
! message is printed and the grid point value filled
! with the undefined value.
!-------------------------------------------------------

            if( ktrop.EQ.kstart .or. &
                ktrop.EQ.kend   .or. &
                ktrop.EQ.0    ) then

              tropp(i,j) = undef
              tropt(i,j) = undef
              TropQ(i,j) = undef

            else

!------------------------------------------------------
! If the tropopause has been successfully isolated
!     store tropopause pressure    in TropP
! and store tropopause temperature in TropT.
!------------------------------------------------------

              tropp(i,j) =   pl(i,j,ktrop)
              tropt(i,j) = TmpU(i,j,ktrop)
              tropq(i,j) = SphU(i,j,ktrop)

            end if

        end do
        end do

      return

      end subroutine tropopause

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine VertInterp(v2,v3,ple,pp,rc)

    real(r4), intent(OUT) :: v2(:,:)
    real(r8), intent(IN ) :: v3(:,:,:)
    real(r8), intent(IN ) :: ple(:,:,:)
    real    , intent(IN ) :: pp
    integer, optional, intent(OUT) :: rc

    real, dimension(size(v2,1),size(v2,2)) :: al,PT,PB
    integer km
    logical edge

    integer        :: status
    character*(10) :: Iam='VertInterp'

    km   = size(ple,3)-1
    edge = size(v3,3)==km+1

    ASSERT_(edge .or. size(v3,3)==km)

    v2   = MAPL_UNDEF

    if(EDGE) then
       pb   = ple(:,:,km+1)
       do k=km,1,-1
          pt = ple(:,:,k)
          if(all(pb<pp)) exit
          where(pp>pt .and. pp<=pb)
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k)*al + v3(:,:,k+1)*(1.0-al)
          end where
          pb = pt
       end do
    else
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
       do k=km,2,-1
          pt = 0.5*(ple(:,:,k-1)+ple(:,:,k))
          if(all(pb<pp)) exit
          where( (pp>pt.and.pp<=pb) )
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k-1)*al + v3(:,:,k)*(1.0-al)
          end where
          pb = pt
       end do
       pt = 0.5*(ple(:,:,km)+ple(:,:,km-1))
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
          where( (pp>pb.and.pp<=ple(:,:,km+1)) )
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,km-1)*al + v3(:,:,km)*(1.0-al)
          end where
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine VertInterp


!BOP

! !IROUTINE: Coldstart

! !DESCRIPTION:
!   Routine to coldstart from an isothermal state of rest.
!   The temperature can be specified in the config, otherwise
!   it is 300K. The surface pressure is assumed to be 1000 hPa.
!
! !INTERFACE:

subroutine Coldstart(gc, import, export, clock, rc)

    USE jw, only : temperature, u_wind, v_wind, surface_geopotential
    USE jw, only : tracer_q, tracer_q1_q2, tracer_q3
    USE testcases_3_4_5_6, only : advection, Rossby_Haurwitz, mountain_Rossby, gravity_wave

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc
    type(ESMF_State),    intent(inout) :: import
    type(ESMF_State),    intent(inout) :: export
    type (ESMF_Clock),   intent(in)    :: clock
    integer, intent(out), optional     :: rc
 
!EOP

    character(len=ESMF_MAXSTR)        :: IAm
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    integer                           :: status

    type (MAPL_MetaComp),     pointer :: MAPL 
    type (ESMF_State)                 :: INTERNAL

    real(r8), pointer                 :: AK(:), BK(:)
    real(r8), pointer                 :: U      (:,:,:)
    real(r8), pointer                 :: V      (:,:,:)
    real(r8), pointer                 :: PT     (:,:,:)
    real(r8), pointer                 :: PE     (:,:,:)
    real(r8), pointer                 :: PKL    (:,:,:)
    real(kind=4), pointer             :: phis   (:,:)
    real, pointer                     :: LONS   (:,:)
    real, pointer                     :: LATS   (:,:)
    real                              :: T0
    integer                           :: L
    type(ESMF_Config)                 :: CF
    integer                           :: i,j,k
    integer                           :: IS,IE, JS,JE, KS,KE, KM

    integer                     :: case_id
    integer                     :: case_rotation
    character(len=ESMF_MAXSTR)  :: case_tracers

    integer  :: p_ks
    real(r8) :: dummyr8_1, dummyr8_2, dummyr8_3, dummyr8_4, dummyr8_5, dummyr8_6
    real(r8) :: dz, ztop, height, pressure
    real(r8) :: LONc,LATc
    real(r8) :: ptop, eta, eta_top, rot_ang
    real(r8), allocatable :: PS(:,:)
    logical :: perturb
    logical :: ak_is_missing = .false.
    logical :: bk_is_missing = .false.

    type (DYN_wrap) :: wrap
    type (DynState), pointer :: STATE
    type (DynGrid),  pointer :: GRID

! Tracer Stuff
    real(r4), pointer                :: TRACER(:,:,:)
    real(r8), allocatable            :: Q5(:,:,:)
    real(r8), allocatable            :: Q6(:,:,:)
    type (ESMF_Grid)                 :: esmfGRID 
    integer :: grid_size(3)
    type (ESMF_Bundle)               :: QTR_BUNDLE
    character(len=ESMF_MAXSTR)       :: FIELDNAME

! Begin

    Iam = "Coldstart"
    call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_UserCompGetInternalState(GC, 'DYNstate', wrap, status)
    VERIFY_(STATUS)
    state => wrap%dyn_state
    grid  => state%grid   ! direct handle to grid

    call MAPL_TimerOn(MAPL,"TOTAL")
!BOR    
! !RESOURCE_ITEM: K :: Value of isothermal temperature on coldstart
    call MAPL_GetResource ( MAPL, T0, 'T0:', default=273., RC=STATUS )
    VERIFY_(STATUS)
!EOR
    call MAPL_Get ( MAPL,                &
           INTERNAL_ESMF_STATE=INTERNAL, &
           lats = LATS,                  &
           lons = LONS,                  &
                               RC=STATUS )
    VERIFY_(STATUS)


! A-Grid U Wind
        call MAPL_GetPointer(Internal,U,'U'  ,rc=STATUS)
        VERIFY_(STATUS)
! A-Grid V Wind
        call MAPL_GetPointer(Internal,V,'V'  ,rc=STATUS)
! Surface Geopotential
        call MAPL_GetPointer ( IMPORT, phis, 'PHIS', RC=STATUS )
        VERIFY_(STATUS)
! Potential-Temperature
        call MAPL_GetPointer(Internal,PT,'PT',rc=STATUS)
        VERIFY_(STATUS)
! Edge Pressures
        call MAPL_GetPointer(Internal,PE  ,'PE',rc=STATUS)
        VERIFY_(STATUS)
! Presssure ^ kappa at mid-layers
        call MAPL_GetPointer(Internal,PKL ,'PKZ',rc=STATUS)
        VERIFY_(STATUS)
! AK and BK for vertical coordinate
        call MAPL_GetPointer(Internal,ak  ,'AK' ,rc=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(Internal,bk  ,'BK' ,rc=STATUS)
        VERIFY_(STATUS)


    U = 0.0

    IS = lbound(U,1)
    IE = ubound(U,1)
    JS = lbound(U,2)
    JE = ubound(U,2)
    KS = lbound(U,3)
    KE = ubound(U,3)
    KM = KE-KS+1

    ALLOCATE( PS(IS:IE,JS:JE) )

  if (KM<=2) then   ! Shallow Water

    ! Init State in DynInit based on case number from imput resource file

  else              ! 3-D Baroclinic

    U(IS:IE,JS:JE,KE) = .001*abs(lats(:,:))
    V = 0.0

    call ESMF_ConfigFindLabel( cf, 'AK:', rc = status )
    if (STATUS == 0) then
       do L = 0, SIZE(AK)-1
          call ESMF_ConfigNextLine  ( CF, rc=STATUS )
          call ESMF_ConfigGetAttribute( cf, AK(L), rc = status )
          VERIFY_(STATUS)
       enddo
    else
       ak_is_missing = .true.
    endif

    call ESMF_ConfigFindLabel( cf, 'BK:', rc = status )
    if (STATUS == 0) then
       do L = 0, SIZE(bk)-1
          call ESMF_ConfigNextLine  ( CF, rc=STATUS )
          call ESMF_ConfigGetAttribute( cf, BK(L), rc = status )
          VERIFY_(STATUS)
       enddo
    else
       bk_is_missing = .true.
    endif

    if (ak_is_missing .or. bk_is_missing) call set_eta(km, ptop, AK, BK, p_ks)

    ASSERT_(ANY(AK /= 0.0) .or. ANY(BK /= 0.0))
    do L=lbound(PE,3),ubound(PE,3)
       PE(:,:,L) = AK(L) + BK(L)*MAPL_P00
    enddo

    PKL = 0.5*(PE(:,:,lbound(PE,3)  :ubound(PE,3)-1) + &
               PE(:,:,lbound(PE,3)+1:ubound(PE,3)  ) )
    PKL = PKL**MAPL_KAPPA

    PT = T0/PKL


! 3D Baroclinic Test Cases

    call ESMF_ConfigGetAttribute( cf, case_id      , label='CASE_ID:'      , default=1 , rc = rc )
    call ESMF_ConfigGetAttribute( cf, case_rotation, label='CASE_ROTATION:', default=0 , rc = rc )
    call ESMF_ConfigGetAttribute( cf, case_tracers , label='CASE_TRACERS:' , default='', rc = rc )

! Parse case_rotation
    if (case_rotation == -1) rot_ang =  0
    if (case_rotation ==  0) rot_ang =  0
    if (case_rotation ==  1) rot_ang = 15
    if (case_rotation ==  2) rot_ang = 30
    if (case_rotation ==  3) rot_ang = 45
    if (case_rotation ==  4) rot_ang = 60
    if (case_rotation ==  5) rot_ang = 75
    if (case_rotation ==  6) rot_ang = 90
    if (case_rotation == -1) then
       grid%f_coriolis_angle = -999
    else
       grid%f_coriolis_angle = rot_ang*PI/180.0
    endif

    if (case_id == 1) then ! Steady State

      perturb = .false.
      do k=1,KM
         eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               U(i,j,k) = u_wind(LONc,LATc,eta,perturb,rot_ang)
               V(i,j,k) = v_wind(LONc,LATc,eta,perturb,rot_ang)
     if (k==1) phis(i,j) = surface_geopotential(LONc,LATc,rot_ang)
               PT(i,j,k) = temperature(LONc,LATc,eta,rot_ang)
            enddo
         enddo
      enddo
      PT = PT/PKL

    elseif (case_id == 2) then ! Baroclinic Wave

      perturb = .true.
      do k=1,KM
         eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               U(i,j,k) = u_wind(LONc,LATc,eta,perturb,rot_ang)
               V(i,j,k) = v_wind(LONc,LATc,eta,perturb,rot_ang)
     if (k==1) phis(i,j) = surface_geopotential(LONc,LATc,rot_ang)
               PT(i,j,k) = temperature(LONc,LATc,eta,rot_ang)
            enddo
         enddo
      enddo
      PT = PT/PKL

    elseif (case_id == 3) then ! Advection

      PURE_ADVECTION = .true.

      allocate( Q5(IS:IE, JS:JE, 1:KM), STAT=STATUS)
      VERIFY_(STATUS)
      allocate( Q6(IS:IE, JS:JE, 1:KM), STAT=STATUS)
      VERIFY_(STATUS)

      ztop = 12000.0
      dz   = ztop/KM
      do k=1,KM
         height = (ztop - 0.5*dz) - (k-1)*dz  ! Layer middle height
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               call  advection('56', LONc, LATc, height, rot_ang,  &
                        U(i,j,k), V(i,j,k), PT(i,j,k), dummyr8_1, &
                        PS(i,j), Q5(i,j,k), Q6(i,j,k))
               phis(i,j) = dummyr8_1
            enddo
         enddo
      enddo
      do L=lbound(PE,3),ubound(PE,3)
         PE(:,:,L) = AK(L) + BK(L)*PS(:,:)
      enddo


      do k=1,KM
         do j=JS,JE
            do i=IS,IE
               PKL(i,j,k) = ( (PE(i,j,k)**kappa) - (PE(i,j,k-1)**kappa) ) /  &
                            ( kappa*(log(PE(i,j,k))-log(PE(i,j,k-1))) )

            enddo
         enddo
      enddo

!!!      PKL = 0.5*(PE(:,:,lbound(PE,3)  :ubound(PE,3)-1) + &
!!!                 PE(:,:,lbound(PE,3)+1:ubound(PE,3)  ) )
!!!      PKL = PKL**MAPL_KAPPA
      PT = PT/PKL

    elseif (case_id == 4) then ! 3D Rossby-Haurwitz

      do j=JS,JE
         do i=IS,IE
            LONc = LONS(i,j)
            LATc = LATS(i,j)
            pressure = 500._r8
            call Rossby_Haurwitz(LONc,LATc, pressure, U(i,j,1), V(i,j,1), PT(i,j,1), dummyr8_1, PS(i,j))
            phis(i,j) = dummyr8_1
         enddo
      enddo
      do L=lbound(PE,3),ubound(PE,3)
         PE(:,:,L) = AK(L) + BK(L)*PS(:,:)
      enddo
      do k=1,KM
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               pressure = 0.5*(PE(i,j,k-1)+PE(i,j,k))
               call Rossby_Haurwitz(LONc,LATc, pressure, U(i,j,k), V(i,j,k), PT(i,j,k), dummyr8_1, PS(i,j))
               phis(i,j) = dummyr8_1
            enddo
         enddo
      enddo

      do k=1,KM
         do j=JS,JE
            do i=IS,IE
               PKL(i,j,k) = ( (PE(i,j,k)**kappa) - (PE(i,j,k-1)**kappa) ) /  &
                            ( kappa*(log(PE(i,j,k))-log(PE(i,j,k-1))) )

            enddo
         enddo
      enddo

!!!      PKL = 0.5*(PE(:,:,lbound(PE,3)  :ubound(PE,3)-1) + &
!!!                 PE(:,:,lbound(PE,3)+1:ubound(PE,3)  ) )
!!!      PKL = PKL**MAPL_KAPPA
      PT = PT/PKL

    elseif (case_id == 5) then ! Mountain-Induced Rossby Wave

      do k=1,KM
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               pressure = 0.5*(PE(i,j,k-1)+PE(i,j,k))
               call mountain_Rossby(case_rotation,LONc,LATc, pressure, U(i,j,k), V(i,j,k), PT(i,j,k), dummyr8_1, PS(i,j))
               phis(i,j) = dummyr8_1
            enddo
         enddo
      enddo
      do L=lbound(PE,3),ubound(PE,3)
         PE(:,:,L) = AK(L) + BK(L)*PS(:,:)
      enddo

      do k=1,KM
         do j=JS,JE
            do i=IS,IE
               PKL(i,j,k) = ( (PE(i,j,k)**kappa) - (PE(i,j,k-1)**kappa) ) /  &
                            ( kappa*(log(PE(i,j,k))-log(PE(i,j,k-1))) )

            enddo
         enddo
      enddo

!!!      PKL = 0.5*(PE(:,:,lbound(PE,3)  :ubound(PE,3)-1) + &
!!!                 PE(:,:,lbound(PE,3)+1:ubound(PE,3)  ) )
!!!      PKL = PKL**MAPL_KAPPA
      PT = PT/PKL

    elseif (case_id == 6) then ! Gravity Waves

   ! case_rotation index has different meaning for this test
      if (case_rotation < 3) then
         grid%f_coriolis_angle = -999
      else
         grid%f_coriolis_angle = 0.0
      endif
   ! Get ICs
      ztop = 10000.d0
      dz   = ztop/KM
      do k=1,KM
         height = (ztop - 0.5d0*dz) - (k-1)*dz  ! Layer middle height
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               call gravity_wave(case_rotation, LONc,LATc, height, U(i,j,k), V(i,j,k), PT(i,j,k), dummyr8_1, PS(i,j))
               phis(i,j) = dummyr8_1
            enddo
         enddo
      enddo
   ! Reconstruct Edge Pressures and AK BK arrays for rotation=0, otherwise use values from set_eta which are OK
      if (case_rotation == 0) then
      PTOP = 27381.905d0
      do k=lbound(PE,3),ubound(PE,3)
         height = ztop - k*dz  ! Layer edge height
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               call gravity_wave(case_rotation, LONc,LATc, height, dummyr8_1, dummyr8_2, dummyr8_3, dummyr8_4, dummyr8_5, pressure=PE(i,j,k))
               eta     = PE(i,j,k)/PS(i,j)
               eta_top = PTOP/PS(i,j)
               BK(k) = (eta - eta_top)/(1.d0 - eta_top)
               AK(k) = 100000.d0 * (eta - BK(k))
            enddo
         enddo
      enddo
      endif
    ! Update PE, PKL and PT
      do L=lbound(PE,3),ubound(PE,3)
         PE(:,:,L) = AK(L) + BK(L)*PS(:,:)
      enddo

      do k=1,KM
         do j=JS,JE
            do i=IS,IE
               PKL(i,j,k) = ( (PE(i,j,k)**kappa) - (PE(i,j,k-1)**kappa) ) /  &
                            ( kappa*(log(PE(i,j,k))-log(PE(i,j,k-1))) )

            enddo
         enddo
      enddo

!!!      PKL = 0.5d0*(PE(:,:,lbound(PE,3)  :ubound(PE,3)-1) + &
!!!                   PE(:,:,lbound(PE,3)+1:ubound(PE,3)  ) )
!!!      PKL = PKL**MAPL_KAPPA
      PT = PT/PKL

    endif ! case_id

!--------------------
! Parse Tracers
!--------------------


      call ESMF_StateGetBundle(IMPORT, 'QTR' , QTR_BUNDLE,   RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_GridCompGet(gc, grid=esmfGRID, rc=STATUS)
      VERIFY_(STATUS)

      allocate( TRACER(IS:IE, JS:JE, 1:KM), STAT=STATUS)
      VERIFY_(STATUS)

      TRACER(:,:,:)  = 0.0_r4
      FIELDNAME = 'Q'
      call addTracer(STATE, QTR_BUNDLE, TRACER, esmfGRID, FIELDNAME)

    if (case_tracers(1:1) /= '0') then
!-----------------------------------------------------------------------
!     tracer q1
!-----------------------------------------------------------------------
      if (case_tracers(1:1) == '1' .or. case_tracers(2:2) == '1' .or. case_tracers(3:3) == '1' .or. &
          case_tracers(4:4) == '1' .or. case_tracers(5:5) == '1') then
      TRACER(:,:,:) = 0.0_r8
      do k=1,KM
         eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               dummyr8_1 = tracer_q1_q2(LONc,LATc,eta,rot_ang,0.6_r8)
               TRACER(i,j,k) = dummyr8_1
            enddo
         enddo
      enddo
      FIELDNAME = 'Q1'
      call addTracer(STATE, QTR_BUNDLE, TRACER, esmfGRID, FIELDNAME)
      endif

!-----------------------------------------------------------------------
!     tracer q2
!-----------------------------------------------------------------------
      if (case_tracers(1:1) == '2' .or. case_tracers(2:2) == '2' .or. case_tracers(3:3) == '2' .or. &
          case_tracers(4:4) == '2' .or. case_tracers(5:5) == '2') then
      do k=1,KM
         eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               dummyr8_1 = tracer_q1_q2(LONc,LATc,eta,rot_ang,1.0_r8)
               TRACER(i,j,k) = dummyr8_1
            enddo
         enddo
      enddo
      FIELDNAME = 'Q2'
      call addTracer(STATE, QTR_BUNDLE, TRACER, esmfGRID, FIELDNAME)
      endif

!-----------------------------------------------------------------------
!     tracer q3
!-----------------------------------------------------------------------
      if (case_tracers(1:1) == '3' .or. case_tracers(2:2) == '3' .or. case_tracers(3:3) == '3' .or. &
          case_tracers(4:4) == '3' .or. case_tracers(5:5) == '3') then
      do k=1,KM
         eta = 0.5*( (ak(k-1)+ak(k))/1.e5 + bk(k-1)+bk(k) )
         do j=JS,JE
            do i=IS,IE
               LONc = LONS(i,j)
               LATc = LATS(i,j)
               dummyr8_1 = tracer_q3(LONc,LATc,eta,rot_ang)
               TRACER(i,j,k) = dummyr8_1
            enddo
         enddo
      enddo
      FIELDNAME = 'Q3'
      call addTracer(STATE, QTR_BUNDLE, TRACER, esmfGRID, FIELDNAME)
      endif

!-----------------------------------------------------------------------
!     tracer q4
!-----------------------------------------------------------------------
      if (case_tracers(1:1) == '4' .or. case_tracers(2:2) == '4' .or. case_tracers(3:3) == '4' .or. &
          case_tracers(4:4) == '4' .or. case_tracers(5:5) == '4') then
      TRACER(:,:,:)  = 1.0_r4
      FIELDNAME = 'Q4'
      call addTracer(STATE, QTR_BUNDLE, TRACER, esmfGRID, FIELDNAME)
      VERIFY_(STATUS)
      endif

!-----------------------------------------------------------------------
!     tracer q5
!-----------------------------------------------------------------------
      if (allocated(Q5)) then
      TRACER(:,:,:)  = Q5(:,:,:)
      FIELDNAME = 'Q5'
      call addTracer(STATE, QTR_BUNDLE, TRACER, esmfGRID, FIELDNAME)
      VERIFY_(STATUS)
      deallocate( Q5, STAT=STATUS)
      VERIFY_(STATUS)
      endif

!-----------------------------------------------------------------------
!     tracer q6
!-----------------------------------------------------------------------
      if (allocated(Q6)) then
      TRACER(:,:,:)  = Q6(:,:,:)
      FIELDNAME = 'Q6'
      call addTracer(STATE, QTR_BUNDLE, TRACER, esmfGRID, FIELDNAME)
      VERIFY_(STATUS)
      deallocate( Q6, STAT=STATUS)
      VERIFY_(STATUS)
      endif

    endif

      deallocate( TRACER, STAT=STATUS)
      VERIFY_(STATUS)

    endif

    DEALLOCATE( PS )

    call MAPL_TimerOff(MAPL,"TOTAL")

    DYN_COLDSTART=.true.

    RETURN_(ESMF_SUCCESS)
  end subroutine COLDSTART

 subroutine set_eta(km, ptop, ak, bk, ks)

      integer,  intent(in   )::  km          ! vertical dimension
      real*8,   intent(  out):: ptop         ! model top (Pa)
      integer,  intent(  out):: ks           ! number of pure p layers
      real*8,   intent(inout):: ak(km+1)
      real*8,   intent(inout):: bk(km+1)

! local
      real*8 a20_01(21),b20_01(21)      ! NCAR Colloquium 20-levels N=0.01
      real*8 a20_0178(21),b20_0178(21)  ! NCAR Colloquium 20-levels N=0.0178
      real*8 a26(27),b26(27)            ! NCAR Colloquium 26-levels

      real*8 :: p0=1000.E2
      real*8 :: pc=200.E2
      real*8 pt, pint, lnpe, dlnp
      real*8 press(km+1)
      integer  k

      data a20_01 / 0.27381905404907E+05,  0.26590539035976E+05,  0.25752394878279E+05,  0.24865429808716E+05, &
                 0.23927536347865E+05,  0.22936541085572E+05,  0.21890203071294E+05,  0.20786212168493E+05, &
                 0.19622187372385E+05,  0.18395675090318E+05,  0.17104147384052E+05,  0.15745000173179E+05, &
                 0.14315551398919E+05,  0.12813039147516E+05,  0.11234619732416E+05,  0.95773657344247E+04, &
                 0.78382639990006E+04,  0.60142135898353E+04,  0.41020236978492E+04,  0.20984115047143E+04, &
                 0.00000000000000E+00 /

      data b20_01 / 0.00000000000000E+00,  0.28901070149364E-01,  0.59510487036309E-01,  0.91902866472543E-01, &
                 0.12615517459290E+00,  0.16234678535331E+00,  0.20055953931639E+00,  0.24087780374962E+00, &
                 0.28338853406205E+00,  0.32818133660555E+00,  0.37534853286773E+00,  0.42498522508382E+00, &
                 0.47718936329560E+00,  0.53206181388604E+00,  0.58970642961892E+00,  0.65023012121324E+00, &
                 0.71374293048299E+00,  0.78035810507338E+00,  0.85019217482527E+00,  0.92336502980036E+00, &
                 0.10000000000000E+01 /

      data a20_0178 / 0.32021324453921E+05,  0.31137565415634E+05,  0.30202026400316E+05,  0.29211673587770E+05, &
                      0.28163295404433E+05,  0.27053492108706E+05,  0.25878664766072E+05,  0.24635003578258E+05, &
                      0.23318475528610E+05,  0.21924811303582E+05,  0.20449491447964E+05,  0.18887731708932E+05, &
                      0.17234467521390E+05,  0.15484337584307E+05,  0.13631666474783E+05,  0.11670446243450E+05, &
                      0.95943169315531E+04,  0.73965459465018E+04,  0.50700062290314E+04,  0.26071531411601E+04, &
                      0.00000000000000E+00 /

      data b20_0178 / 0.00000000000000E+00,  0.27599078219223E-01,  0.56815203138214E-01,  0.87743118501982E-01, & 
                      0.12048311914891E+00,  0.15514137625266E+00,  0.19183028162025E+00,  0.23066881216269E+00, &
                      0.27178291572025E+00,  0.31530591949337E+00,  0.36137896240390E+00,  0.41015145278854E+00, &
                      0.46178155290889E+00,  0.51643669184922E+00,  0.57429410846515E+00,  0.63554142614418E+00, &
                      0.70037726124166E+00,  0.76901186716541E+00,  0.84166781619770E+00,  0.91858072126555E+00, &
                      0.10000000000000E+01 /


      data a26 /  219.4067,   489.5209,   988.2418,  1805.2010,  2983.7240,  4462.3340,   &
                 6160.5870,  7851.2430,  7731.2710,  7590.1310,  7424.0860,   &
                 7228.7440,  6998.9330,  6728.5740,  6410.5090,  6036.3220,   &
                 5596.1110,  5078.2250,  4468.9600,  3752.1910,  2908.9490,   &
                  2084.739,   1334.443,    708.499,   252.1360,  0.0, 0.0     /

      data b26 / 0.0, 0.0, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,&
                 0.0000000, 0.01505309, 0.03276228, 0.05359622, 0.07810627,      &
                 0.1069411, 0.1408637, 0.1807720, 0.2277220, 0.2829562,       &
                 0.3479364, 0.4243822, 0.5143168, 0.6201202, 0.7235355,       &
                 0.8176768, 0.8962153, 0.9534761, 0.9851122, 1.0000000        /

      SELECT CASE(km)
  
      CASE(20)

          do k=1,km+1
            ak(k) = a20_0178(k)
            bk(k) = b20_0178(k)
          enddo
! Search KS
            ks = 0
          do k=1,km
             if(bk(k) > 0) then
                ks = k-1
                goto 120
             endif
          enddo
120   continue

      CASE(26)

          do k=1,km+1
            ak(k) = a26(k)
            bk(k) = b26(k)
          enddo
! Search KS
            ks = 0
          do k=1,km
             if(bk(k) > 0) then
                ks = k-1
                goto 126
             endif
          enddo
126   continue
 
      CASE(40)
!--------------------------------------------------
! Pure sigma-coordinate with uniform spacing in "z"
!--------------------------------------------------
!
         ptop = 27381.905404907        ! model top pressure (pascal)
         press(1) = ptop
         press(km+1) = p0
         dlnp = (log(p0) - log(ptop)) / real(km)

            lnpe = log(press(km+1))
         do k=km,2,-1
            lnpe = lnpe - dlnp
            press(k) = exp(lnpe)
         enddo

! Search KS
            ks = 0
         do k=1,km
            if(press(k) >= pc) then
               ks = k-1
               goto 140
            endif
         enddo
140   continue

         if(ks /= 0) then
            do k=1,ks
               ak(k) = press(k)
               bk(k) = 0.
            enddo
          endif

             pint = press(ks+1)
          do k=ks+1,km
             ak(k) =  pint*(press(km)-press(k))/(press(km)-pint)
             bk(k) = (press(k) - ak(k)) / press(km+1)
          enddo
             ak(km+1) = 0.
             bk(km+1) = 1.

      CASE(60)
!--------------------------------------------------
! Pure sigma-coordinate with uniform spacing in "z"
!--------------------------------------------------
!
         ptop = 25499.234876157        ! model top pressure (pascal)
         press(1) = ptop
         press(km+1) = p0
         dlnp = (log(p0) - log(ptop)) / real(km)

            lnpe = log(press(km+1))
         do k=km,2,-1
            lnpe = lnpe - dlnp
            press(k) = exp(lnpe)
         enddo

! Search KS
            ks = 0
         do k=1,km
            if(press(k) >= pc) then
               ks = k-1
               goto 160
            endif
         enddo
160   continue

         if(ks /= 0) then
            do k=1,ks
               ak(k) = press(k)
               bk(k) = 0.
            enddo
          endif

             pint = press(ks+1)
          do k=ks+1,km
             ak(k) =  pint*(press(km)-press(k))/(press(km)-pint)
             bk(k) = (press(k) - ak(k)) / press(km+1)
          enddo
             ak(km+1) = 0.
             bk(km+1) = 1.

     CASE DEFAULT

        print*, 'Bad KM in DynCore_GridCompMod:set_eta', km

     END SELECT


 end subroutine set_eta

subroutine addTracer_r8(state, bundle, var, grid, fieldname)
  type (DynState), pointer         :: STATE
  type (ESMF_Bundle)               :: BUNDLE
  real(r8), pointer                :: var(:,:,:)
  type (ESMF_Grid)                 :: GRID
  character(len=ESMF_MAXSTR)       :: FIELDNAME

  integer :: isd,ied,jsd,jed,npz,nq, rc,status
  type(DynTracers), pointer        :: t(:)

  character(len=ESMF_MAXSTR)       :: IAm='FV:addTracer_r8'

  type (ESMF_Field)                :: field
  type (ESMF_Array)                :: array
  type (ESMF_FieldDataMap) :: DATAMAP

! Count the friendlies
!---------------------
      call ESMF_BundleGet(BUNDLE, fieldCount=NQ, RC=STATUS)
      VERIFY_(STATUS)

      NQ = NQ + 1

      ARRAY = ESMF_ArrayCreate(var, ESMF_DATA_COPY, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldDataMapSetDefault(DATAMAP, 3, rc=status)
      VERIFY_(STATUS)
      field = ESMF_FieldCreate( GRID, ARRAY, &
                                horzRelloc=ESMF_CELL_CENTER, vertRelloc=ESMF_CELL_CELL, DATAMAP=DATAMAP, &
                                name=fieldname, rc=STATUS)
      VERIFY_(STATUS)
      call ESMF_BundleAddField ( bundle, field, rc=STATUS )
      VERIFY_(STATUS)

      if (NQ == 1) then
         ALLOCATE(STATE%VARS%tracer(nq), STAT=STATUS)
         VERIFY_(STATUS)
         call ESMF_FieldGetDataPointer( field, state%vars%tracer(nq)%content, rc=STATUS )
         VERIFY_(STATUS)
         state%vars%tracer(nq  )%is_r4 = .false.
      else
         allocate(t(nq))
         t(1:nq-1) = state%vars%tracer
         deallocate(state%vars%tracer)
         state%vars%tracer => t
         call ESMF_FieldGetDataPointer( field, state%vars%tracer(nq)%content, rc=STATUS )
         state%vars%tracer(nq  )%is_r4 = .false.
      endif

      STATE%GRID%NQ = NQ

  return
end subroutine addTracer_r8

subroutine addTracer_r4(state, bundle, var, grid, fieldname)
  type (DynState), pointer         :: STATE
  type (ESMF_Bundle)               :: BUNDLE
  real(r4), pointer                :: var(:,:,:)
  type (ESMF_Grid)                 :: GRID
  character(len=ESMF_MAXSTR)       :: FIELDNAME

  integer :: isd,ied,jsd,jed,npz,nq, rc,status
  type(DynTracers), pointer        :: t(:)

  character(len=ESMF_MAXSTR)       :: IAm='FV:addTracer_r4'

  type (ESMF_Field)                :: field
  type (ESMF_Array)                :: array
  type (ESMF_FieldDataMap) :: DATAMAP

! Count the friendlies
!---------------------
      call ESMF_BundleGet(BUNDLE, fieldCount=NQ, RC=STATUS)
      VERIFY_(STATUS)

      NQ = NQ + 1

      ARRAY = ESMF_ArrayCreate(var, ESMF_DATA_COPY, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldDataMapSetDefault(DATAMAP, 3, rc=status)
      VERIFY_(STATUS)

      field = ESMF_FieldCreate( GRID, ARRAY, &
                                horzRelloc=ESMF_CELL_CENTER, vertRelloc=ESMF_CELL_CELL, DATAMAP=DATAMAP, &
                                name=fieldname, rc=STATUS)
      VERIFY_(STATUS)
      call ESMF_BundleAddField ( bundle, field, rc=STATUS )
      VERIFY_(STATUS)

      if (NQ == 1) then
         ALLOCATE(STATE%VARS%tracer(nq), STAT=STATUS)
         VERIFY_(STATUS)
         call ESMF_FieldGetDataPointer( field, state%vars%tracer(nq)%content_r4, rc=STATUS )
         VERIFY_(STATUS)
         state%vars%tracer(nq  )%is_r4 = .true.
      else
         allocate(t(nq))
         t(1:nq-1) = state%vars%tracer
         deallocate(state%vars%tracer)
         state%vars%tracer => t
         call ESMF_FieldGetDataPointer( field, state%vars%tracer(nq)%content_r4, rc=STATUS )
         state%vars%tracer(nq  )%is_r4 = .true.
      endif

      STATE%GRID%NQ = NQ

  return
end subroutine addTracer_r4

subroutine freeTracers(state)
  type (DynState) :: STATE

  if (associated(STATE%VARS%tracer)) then
     DEALLOCATE( STATE%VARS%tracer)   ! Comment out to output tracer to checkpoint file
  end if

  return
end subroutine freeTracers

#if defined(DEBUG)         
  Subroutine Write_Profile_R8(grid, arr, name)
    type (DynGrid),   intent(IN) :: grid
    real(r8),         intent(IN) :: arr(grid%is:grid%ie,grid%js:grid%je,1:grid%npz)
    character(len=*), intent(IN) :: name

    integer  :: istrt,iend, jstrt,jend, kstrt,kend
    integer  :: im, jm, km, k
    real(r8) :: arr_global(grid%npx,grid%ntiles*grid%npy,grid%npz)
    real(r8) :: rng(3,grid%npz)
    real(r8) :: GSUM

    type (ESMF_Grid)            :: esmfGRID
    real(kind=ESMF_KIND_R8)     :: locArr(grid%is:grid%ie,grid%js:grid%je)
    real(kind=ESMF_KIND_R8)     :: glbArr(grid%npx,grid%ntiles*grid%npy)

    istrt = grid%is
    iend  = grid%ie
    jstrt = grid%js
    jend  = grid%je
    kstrt = 1
    kend  = grid%npz
    im    = grid%npx
    jm    = grid%npy*grid%ntiles
    km    = grid%npz

    call write_parallel('GlobalSUm')
    do k=kstrt,kend
       locArr(:,:) = arr(:,:,k)
       call ArrayGather(locArr, glbArr, grid%grid)
       arr_global(:,:,k) = glbArr
    enddo

    IF (MAPL_AM_I_ROOT()) Then
       rng(1,:) = MINVAL(MINVAL(arr_global,DIM=1),DIM=1)
       rng(2,:) = MAXVAL(MAXVAL(arr_global,DIM=1),DIM=1)
       rng(3,:) = SUM(SUM(arr_global,DIM=1),DIM=1)/(IM*JM)
       GSUM     = SUM(SUM(SUM(arr_global,DIM=1),DIM=1),DIM=1)

       print*,'***********'
       print*,'stats for ',trim(name)

       Do k = 1, km
          Write(*,'(a,i4.0,3(f21.9,1x))')'k:',k,rng(:,k)
       End Do
       Write(*,"('GlobalSum: ',f21.9)") GSUM
       print*,'***********'
       print*,' '
    End IF

  End Subroutine Write_Profile_R8

  Subroutine Write_Profile_R4(grid, arr, name, delp)
    type (DynGrid),   intent(IN) :: grid
    real(r4),         intent(IN) :: arr(grid%is:grid%ie,grid%js:grid%je,1:grid%npz)
    character(len=*), intent(IN) :: name
    real(r8), optional, intent(IN) :: delp(grid%is:grid%ie,grid%js:grid%je,1:grid%npz)

    integer  :: istrt,iend, jstrt,jend, kstrt,kend
    integer  :: im, jm, km, k
    real(r4) :: arr_global(grid%npx,grid%ntiles*grid%npy,grid%npz)
    real(r4) :: rng(3,grid%npz)
    real(r8) :: gsum_p
    real(r4) :: GSUM
    
    type (ESMF_Grid)            :: esmfGRID
    real(kind=ESMF_KIND_R8)     :: locArr(grid%is:grid%ie,grid%js:grid%je)
    real(kind=ESMF_KIND_R8)     :: glbArr(grid%npx,grid%ntiles*grid%npy)
      
    istrt = grid%is
    iend  = grid%ie
    jstrt = grid%js
    jend  = grid%je 
    kstrt = 1
    kend  = grid%npz
    im    = grid%npx
    jm    = grid%npy*grid%ntiles      
    km    = grid%npz
    
    do k=kstrt,kend
       locArr(:,:) = arr(:,:,k)       
       call ArrayGather(locArr, glbArr, grid%grid)
       arr_global(:,:,k) = glbArr
    enddo
    IF (MAPL_AM_I_ROOT()) Then
       rng(1,:) = MINVAL(MINVAL(arr_global,DIM=1),DIM=1)
       rng(2,:) = MAXVAL(MAXVAL(arr_global,DIM=1),DIM=1)
       rng(3,:) = SUM(SUM(arr_global,DIM=1),DIM=1)/(IM*JM)
       print*,'***********'
       print*,'stats for ',trim(name)
       Do k = 1, km
          Write(*,'(a,i4.0,3(f21.9,1x))')'k:',k,rng(:,k)
       End Do
       print*,'***********'
       print*,' '
    End IF

    if (present(delp)) then
    gsum_p = 0
    do k=kstrt,kend
       locArr(:,:) = arr(:,:,k)*grid%area(:,:)*delp(:,:,k)
       call ArrayGather(locArr, glbArr, grid%grid)
       arr_global(:,:,k) = glbArr
       locArr(:,:) = delp(:,:,k)
       call ArrayGather(locArr, glbArr, grid%grid)
       gsum_p = gsum_p + SUM(SUM(glbArr,DIM=1),DIM=1)
    enddo
    IF (MAPL_AM_I_ROOT()) Then
       GSUM     = SUM(SUM(SUM(arr_global,DIM=1),DIM=1),DIM=1)
       print*,'***********'
       Write(*,"('GlobalSum: ',e21.9)") GSUM/(grid%globalarea*gsum_p)
       print*,'***********'
       print*,' '
    End IF
    endif

  End Subroutine Write_Profile_R4
#endif

end module DynCore_GridCompMod

