E4TcadF40oQ32.R GISS Model E  1850 ocn/atm/tracers   M. Kelley  07/2010

E4TcadF40oQ32: E4TcadF40 coupled to 1x1.25deg 32-layer GISS ocean model
E4TcadF40: E4F40 + Dust, Chemistry and aerosol tracers

E4F40: modelE as frozen (or not yet) in July 2009
modelE4 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define CHECK_OCEAN              ! needed to compile aux/file CMPE002
#define TRAC_ADV_CPU
#define USE_ENT                  ! include dynamic vegetation model
!  OFF #define NEW_IO
#define TRACERS_ON               ! include tracers code
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
ORES_1Qx1_L32 OFFT288E              ! ocean horiz res 1x1.25deg, 32 layers
OSTRAITS_1QX1_COM                   ! ocean straits settings for this res.

    ! lat-lon grid specific source codes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT_netcdf                ! diagn/post-processing output
IORSF                               ! i/o
! or IO_DRV TRDIAG                       ! new i/o

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

#include "dynamic_ocn_source_files"
OCN_Int_LATLON                      ! atm-ocn regrid routines

Components:
#include "E4_components"    /* without "Ent" */
Ent
! dd2d  /* needed for new i/o only */

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB    /* needed for "Ent" only */
OPTS_giss_LSM = USE_ENT=YES           /* needed for "Ent" only */
! OPTS_dd2d = NC_IO=PNETCDF           /* an OPTION for new i/o */

Data input files:
#include "IC_144x90_input_files"
#include "dynamic_ocn_288x180_input_files"
TOPO=Z2HX2fromZ1QX1N            ! atm-grid surface fractions and topography
RVR=RD_Fb.RVR.bin               ! river direction file

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
E4TcadF40oQ32 (E4F40 + chem,aero,dust tracers + Q32 Ocean R)


&&PARAMETERS
#include "dynamic_ocn_params"
#include "sdragF40_params"
#include "gwdragF40_params"

! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.74 ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.65 ! below 850mb and MC regions; tune this last  to get rad.balance

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

DTO=112.5        ! ocean dynamics timestep
DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)
#include "diag_params"

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960
&&END_PARAMETERS

 &INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1949,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
 &END
