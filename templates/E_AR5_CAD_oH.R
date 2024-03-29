E_AR5_CAD_oH.R GISS Model E  1850 ocn/atm/tracers  tnl   10/15/2010

E_AR5_CAD_oH: E_AR5_CAD  coupled to 1x1deg 26-layer Hybrid-Isopycnal Coordinate Ocean Model (HYCOM)

E_AR5_CAD_oH : sibling E_AR5_CAD_oR
modelE4 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean: coupled to 1x1deg 26-layer HYCOM
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define USE_ENT
#define NEW_IO
#define CHECK_OCEAN                  ! needed to compile aux/file CMPE002
#define TRAC_ADV_CPU
! #define TRACERS_GASEXCH_Natassa    ! special tracers to be passed to ocean
! #define TRACERS_HYCOM_Ventilation
#define ATM2x2h                      ! 2x2.5 40 layer atm
! #define ATM4x5                     ! 4x5 20 layer atm
#define HYCOM1deg                    ! 1deg 26 layer hycom (387x360x26)
! #define HYCOM2deg                  ! 2deg 26 layer hycom (195x180x26)
#define TRACERS_WATER            ! wet deposition and water tracer
#define TRACERS_DUST             ! include dust tracers
#define TRACERS_DUST_Silt4       ! include 4th silt size class of dust
#define TRACERS_DRYDEP           ! default dry deposition
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
#define NO_HDIURN                ! exclude hdiurn diagnostics
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
#define SHINDELL_STRAT_CHEM         ! turns on stratospheric chemistry
#define RAD_O3_GCM_HRES     ! Use GCM horiz resl to input rad code clim Ozone
!  OFF #define AUXILIARY_OX_RADF ! radf diags for climatology or tracer Ozone
#define TRACERS_TERP                ! include terpenes in gas-phase chemistry
#define BIOGENIC_EMISSIONS       ! turns on interactive isoprene emissions
#define TRACERS_AEROSOLS_Koch    ! Dorothy Koch's tracers (aerosols, etc)
#define TRACERS_AEROSOLS_SOA     ! Secondary Organic Aerosols
!  OFF #define SOA_DIAGS                ! Additional diagnostics for SOA
#define TRACERS_NITRATE
#define TRACERS_HETCHEM
#define BC_ALB                      !optional tracer BC affects snow albedo
!  OFF #define CLD_AER_CDNC              !aerosol-cloud interactions
!  OFF #define BLK_2MOM                  !aerosol-cloud interactions
!  OFF #define WATER_MISC_GRND_CH4_SRC ! adds lake, ocean, misc. ground sources for CH4
!  OFF #define CALCULATE_FLAMMABILITY  ! activated code to determine flammability of surface veg
!  OFF #define DYNAMIC_BIOMASS_BURNING  ! alter biomas burning my flammability
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
IO_DRV TRDIAG                       ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV                             ! advection of T
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

QUS3D                               ! advection of Q and tracers
TRDUST_COM TRDUST TRDUST_DRV        ! dust tracer specific code
#include "tracer_shared_source_files"
#include "tracer_shindell_source_files"
#include "tracer_aerosols_source_files"
! CLD_AEROSOLS_Menon_MBLK_MAT BLK_DRV ! aerosol-cloud interactions

#include "modelE4_source_files"
RAD_native_O3                       ! for reading ozone to rad code at native GCM horiz res.
lightning                           ! Colin Price lightning model
! flammability_drv flammability       ! Olga's fire model

#include "hycom_source_files"

Components:
#include "E4_components_nc"    /* without "Ent" */
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB    /* needed for "Ent" only */
OPTS_giss_LSM = USE_ENT=YES           /* needed for "Ent" only */
  OPTS_dd2d = NC_IO=PNETCDF           /* an OPTION for new i/o */

Data input files:
#include "IC_144x90_input_files_AR5"
#include "hycom_387x360_input_files"
RVR=RD_modelE_Fa.RVR_1deghycom_may10.bin ! river direction file
VEG_DENSE=gsin/veg_dense_2x2.5 ! vegetation density for flammability calculations

#include "land144x90_input_files"
#include "rad_input_files"
#include "TAero2008_input_files"
#include "O3_2010_144x90_input_files"
!#include "O3_2005_input_files"

#include "chemistry_input_files"
#include "chemistry_144x90_input_files"

#include "dust_tracer_144x90_input_files"
#include "dry_depos_144x90_input_files"

#include "chem_emiss_144x90_input_files"

#include "aeros_input_files"

MSU_wts=MSU.RSS.weights.data      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
E_AR5_CAD_oH (E_AR5_NINT + chem.,aer.,dust tracers + Ocean Hycom)


&&PARAMETERS
#include "dynamic_ocn_params"
#include "sdragF40_params"
#include "gwdragF40_params"

! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.56  ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.00  ! below 850mb and MC regions; tune this last  to get rad.balance
WMUI_multiplier = 2.

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=0.      ! activates strat.H2O generated by CH4
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

initial_GHG_setup = 1 ! Set to 0 after initial setup.

#include "atmCompos_1850_params"
!!!!!!!!!!!!!!!!!!!!!!!
! Please note that making o3_yr non-zero tells the model
! to override the transient chemistry tracer emissions'
! use of model year and use abs(o3_yr) instead!
!!!!!!!!!!!!!!!!!!!!!!!
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols
#include "aerosol_params"
imAER=5         !3 historic; 1 AEROCOM ; 0,2 for standard or sector inputs (not working)
imPI=0          !for pre-industrial aerosols (natural-only) use imPI=1, imAER=5, aer_int_yr=1850
aer_int_yr=1850    !used for imAER=3, select desired year (1890 to 2000) or 0 to use JYEAR
#include "dust_params"
#include "chemistry_params"

DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)
nradfrc=0        ! no repeated radiation calculations for inst. rad. forcing
#include "diag_params"

Nssw=48          ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960

itest=-1         ! default is -1
jtest=-1         ! default is -1
iocnmx=2         ! default is 2
brntop=50.       ! default is 50.
brnbot=200.      ! default is 200.
diapyn=2.e-7     ! default is 2.e-7
diapyc=.2e-4     ! default is .2e-4
jerlv0=1         ! default is 1
&&END_PARAMETERS

 &INPUTZ
 YEARI=1899,MONTHI=12,DATEI=01,HOURI=00, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1900,MONTHE=12,DATEE=02,HOURE=00, KDIAG=13*0,
 ISTART=2,IRANDI=0, YEARE=1899,MONTHE=12,DATEE=02,HOURE=00,IWRITE=1,JWRITE=1,
 &END
