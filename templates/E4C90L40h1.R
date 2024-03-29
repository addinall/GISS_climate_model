E4C90L40h1.R GISS Model E  coupled version        M. Kelley   05/2010

E4C90L40h1 = E4C90L40 coupled to coupled to 1-deg 26-layer HYCOM
E4C90L40 = E4F40 on cubed sphere grid 6x90x90 grid boxes
E4F40 = modelE as frozen in April 2010:
Cubed Sphere grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean: coupled to 1-deg 26-layer Hybrid-Isopycnal Coordinate Ocean Model (HYCOM)
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics: finite Volume ; physics 30 min.; radiation 2.5 hrs
filters: none

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define USE_ENT
#define HYCOM1deg                    ! 26 layer 1deg hycom (387x360)
#define NEW_IO
#define CALC_GWDRAG
#define SET_SOILCARBON_GLOBAL_TO_ZERO
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
RES_CS90L40                         ! C90 horiz resolution, top at 0.1mb, 40 layers

     ! Codes used by the cubed-atmosphere configuration (FV dynamics)
#include "cubed_sphere_source_files"

STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)
#include "modelE4_source_files"
#include "hycom_source_files"
cplercs                             ! coupler between cubed atmosphere and hycom

Components:
#include "E4_components"    /* without "Ent" */
Ent
dd2d

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB    /* needed for "Ent" only */
OPTS_giss_LSM = USE_ENT=YES           /* needed for "Ent" only */
OPTS_dd2d = NC_IO=PNETCDF
FVCUBED = YES

Data input files:
#include "IC_CS90_input_files"
#include "hycom_387x360_input_files"
TOPO=Z_C90N.1deghycom                 ! atm topog and surf type fractions
aoremap=remap_c90_hyc387x360.nc      ! cubed-sphere <-> hycom remapping info
taui2o=taua2o_a2x1_o1degr.8bin       ! coupler weights for vector: ice B -> ocean A
tauo2i=ssto2a_a2x1_o1degr.8bin       ! coupler weights for vector: ocean A -> ice B

RVR=RDdistocean_C90.1deghycom         ! river direction file

#include "landCS90_input_files"
#include "rad_input_files"
#include "TAero2008_input_files"
#include "O3_2005_input_files"

MSU_wts=MSU.RSS.weights.data      ! MSU-diag
REG=REG.txt                       ! special regions-diag

Label and Namelist:  (next 2 lines)
E4C90L40h1 (E4C90L40, 1850 atm.;  1-deg HYCOM ocean)


&&PARAMETERS
#include "dynamic_ocn_params"

#include "sdragCS90_params"
#include "gwdragCS90_params"

! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.60 ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.10 ! below 850mb and MC regions; tune this last  to get rad.balance

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=1.      ! activates strat.H2O generated by CH4
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

#include "atmCompos_1850_params"
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols

DTsrc=1800.      ! cannot be changed after a run has been started
DT=1800.         ! for FV dynamics, set same as DTsrc

NIsurf=1         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)

#include "diag_params"

itest=-1         ! default is -1
jtest=-1         ! default is -1
iocnmx=2         ! default is 2
brntop=50.       ! default is 50.
brnbot=200.      ! default is 200.
diapyn=3.e-7     ! default is 3.e-7
diapyc=.5e-4     ! default is .5e-4
jerlv0=1         ! default is 1

Nssw=48          ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960
&&END_PARAMETERS

 &INPUTZ
 YEARI=1899,MONTHI=12,DATEI=01,HOURI=00, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1900,MONTHE=12,DATEE=02,HOURE=00, KDIAG=13*0,
 ISTART=2,IRANDI=0, YEARE=1899,MONTHE=12,DATEE=2,HOURE=0,IWRITE=1,JWRITE=1,
 &END
