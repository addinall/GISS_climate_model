E4TampF40.R GISS Model E Run with MATRIX Aerosols

E4F40: modelE as frozen (or not yet) in July 2009
modelE4 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean data: prescribed, 1876-1885 climatology
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define TRAC_ADV_CPU
#define USE_ENT                  ! include dynamic vegetation model
#define TRACERS_ON               ! include tracers code
#define TRACERS_WATER            ! wet deposition and water tracer
!#define TRACERS_DUST             ! include dust tracers
!#define TRACERS_DUST_Silt4       ! include 4th silt size class of dust
#define TRACERS_DRYDEP           ! default dry deposition
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
#define NO_HDIURN                ! exclude hdiurn diagnostics
!#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
!#define SHINDELL_STRAT_CHEM         ! turns on stratospheric chemistry
#define RAD_O3_GCM_HRES     ! Use GCM horiz resl to input rad code clim Ozone
!  OFF #define AUXILIARY_OX_RADF ! radf diags for climatology or tracer Ozone
!#define TRACERS_TERP                ! include terpenes in gas-phase chemistry
!#define BIOGENIC_EMISSIONS       ! turns on interactive isoprene emissions
!#define INITIAL_GHG_SETUP        ! only for setup hour to get ghg IC file
!#define TRACERS_AEROSOLS_Koch    ! Dorothy Koch's tracers (aerosols, etc)
!#define TRACERS_AEROSOLS_SOA     ! Secondary Organic Aerosols
!  OFF #define SOA_DIAGS                ! Additional diagnostics for SOA
!#define TRACERS_NITRATE
!#define TRACERS_HETCHEM
#define BC_ALB                      !optional tracer BC affects snow albedo
#define TRACERS_AMP
#define TRACERS_AMP_M1
#define CLD_AER_CDNC                ! aerosol - cloud
#define BLK_2MOM              ! aerosol - cloud
#define NEW_IO
!  OFF #define WATER_MISC_GRND_CH4_SRC ! adds lake, ocean, misc. ground sources for CH4
!  OFF #define CALCULATE_FLAMMABILITY  ! activated code to determine flammability of surface veg
!  OFF #define CALCULATE_LIGHTNING ! turn on Colin Price lightning when TRACERS_SPECIAL_Shindell off
!  OFF #define SHINDELL_STRAT_EXTRA     ! non-chemistry stratospheric tracers
!  OFF #define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
!  OFF #define NUDGE_ON                 ! nudge the meteorology
!  OFF #define GFED_3D_BIOMASS          ! turns on IIASA AR4 GFED biomass burning
!  OFF #define HTAP_LIKE_DIAGS    ! adds many diags, changes OH diag, adds Air tracer
!  OFF #define ACCMIP_LIKE_DIAGS  ! adds many diags as defined by ACCMIP project
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
RES_stratF40                        ! horiz/vert resolution, 2x2.5, top at 0.1mb, 40 layers
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform

    ! lat-lon grid specific source codes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT_netcdf                ! diagn/post-processing output
IO_DRV                               ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV                             ! advection of T
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

QUS3D                               ! advection of Q and tracers
TRDUST_COM TRDUST TRDUST_DRV        ! dust tracer specific code

#include "tracer_shared_source_files" 
TRDIAG 

!#include "tracer_shindell_source_files"
!#include "tracer_aerosols_source_files"
TRACERS_AEROSOLS_Koch_e4          
CLD_AEROSOLS_Menon_MBLK_MAT_E29q BLK_DRV ! aerosol-cloud interactions
!                                     Aerosol Micro Physics
TRAMP_drv        |-extend_source  |  
TRAMP_actv       |-extend_source  |  
TRAMP_diam       |-extend_source  | 
TRAMP_nomicrophysics |-extend_source  |  
TRAMP_subs       |-extend_source  |  
TRAMP_coag       |-extend_source  |  
TRAMP_depv       |-extend_source  | 
TRAMP_param_GISS |-extend_source  |    
TRAMP_config  
TRAMP_dicrete    |-extend_source  |    
TRAMP_init       |-extend_source  |  
TRAMP_quad       |-extend_source  |  
TRAMP_matrix     |-extend_source  |        
TRAMP_setup      |-extend_source  |  
TRAMP_npf        |-extend_source  |  
TRAMP_rad        |-extend_source  |
! When using ISORROPIA Thermodynamics
! AMP_thermo_isorr.f
! AMP_isocom.f          
! AMP_isrpia.ext        
! AMP_isofwd          
! AMP_isorev
! When using EQSAM Thermodynamics
TRAMP_thermo_eqsam |-extend_source  |  
TRAMP_eqsam_v03d
! ----------------------------------

#include "modelE4_source_files"
RAD_native_O3                       ! for reading ozone to rad code at native GCM horiz res.
lightning                           ! Colin Price lightning model
! flammability_drv flammability       ! Olga's fire model

#include "static_ocn_source_files"

Components:
shared ESMF_Interface solvers giss_LSM
Ent
dd2d

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB    
OPTS_giss_LSM = USE_ENT=YES           

Data input files:

!  different GIC   #include "IC_144x90_input_files"
    ! start from observed conditions AIC(,OIC), model ground data GIC   ISTART=2
AIC=AIC.RES_F40.D771201      ! observed init cond (atm. only)
GIC=GIC.144X90.DEC01.1.ext.nc   ! initial ground conditions

#include "static_ocn_2000_144x90_input_files"
VEG_DENSE=gsin/veg_dense_2x2.5 ! vegetation density for flammability calculations
RVR=RD_modelE_Fa.RVR.bin          ! river direction file

#include "land144x90_input_files"
#include "rad_input_files"
#include "TAero2008_input_files"
#include "O3_2010_144x90_input_files"
#include "chemistry_144x90_input_files"
#include "dust_tracer_144x90_input_files"
#include "dry_depos_144x90_input_files"
#include "aeros_AMPconstSRC_input_files"

SO2_AIRCRAFT=NOy_sources/aircraft_4x5_1940-2000 ! zero in 1940 and before.
OFFLINE_HNO3.nc=HNO3_dummy_2000_GISS2x2.nc

MSU_wts=MSU.RSS.weights.data      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
E4TampF40 (E4F40 with MATRIX aerosols)


&&PARAMETERS
#include "static_ocn_params"
#include "sdragF40_params"
#include "gwdragF40_params"
CDEF=1.5  ! overwrites value above

! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.74 ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.65 ! below 850mb and MC regions; tune this last  to get rad.balance

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=0.      ! activates strat.H2O generated by CH4
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

!#include "atmCompos_1850_params"
#include "atmCompos_2000_params"

!!!!!!!!!!!!!!!!!!!!!!!
! Please note that making o3_yr non-zero tells the model
! to override the transient chemistry tracer emissions'
! use of model year and use abs(o3_yr) instead!
!!!!!!!!!!!!!!!!!!!!!!!
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols

#include "aerosol_params"
rad_interact_aer=1 ! 1: couples aerosols to radiation, 0: use climatology
diag_rad=1         ! 1: additional radiation diagnostics
OCB_om2oc=1.4      ! biomass burning organic matter to organic carbon ratio (default is 1.4)
BBinc=1.4          ! enhancement factor for carbonaceous aerosols (1.4 for AR5 emissions, 1.0 elsewhere)
!--- number of biomass burning sources (per tracer)
NH3_nBBsources=2
SO2_nBBsources=2
M_BC1_BC_nBBsources=2
M_OCC_OC_nBBsources=2

imAER=5         !3 historic; 1 AEROCOM ; 0,2 for standard or sector inputs (not working)
imPI=0          !for pre-industrial aerosols (natural-only) use imPI=1, imAER=5, aer_int_yr=1850
aer_int_yr=2000    !used for imAER=3, select desired year (1890 to 2000) or 0 to use JYEAR
ad_interact_aer=1  ! 1=couples aerosols to radiation, 0=use climatology

OFFLINE_DMS_SS=0   ! 1= read offline DMS and dust from aerocom file

!------------------  AMP parameters
AMP_DIAG_FC=  2    ! 2=nmode radiation calls  ||  1=one radiation call
AMP_RAD_KEY = 2    ! 1=Volume Mixing || 2=Core-Shell || 3=Maxwell Garnett

#include "dust_params"
imDust=0           ! 0: PDF emission scheme, 1: AEROCOM
COUPLED_CHEM=0     ! to couple chemistry and aerosols


DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=1         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)
#include "diag_params"
KCOPY=1          ! saving acc + rsf

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960
&&END_PARAMETERS

 &INPUTZ
 YEARI=1999,MONTHI=6,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=2010,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1999,MONTHE=6,DATEE=1,HOURE=1,
 &END
