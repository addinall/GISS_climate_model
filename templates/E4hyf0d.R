E4hyf0d.R GISS Model E  1850 ocn/atm                 ssun 06/01/2009

E4hyf0d: replace this section by a description of what distinguishes this run ?
       Use as many lines as you need. Look carefully at all the possible    ?
       choices, particularly the lines containing '?'.
       The final rundeck should contain no '?'
       Check and modify the rest of the description below:                  ?
modelE1 (3.0) 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)     ?
atmospheric composition from year 1979
ocean data: prescribed, 1975-1984 climatology
uses turbulence scheme (no dry conv), simple strat.drag (no grav.wave drag) ?
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs ?
filters: U,V in E-W direction (after every dynamics time step)              ?
         sea level pressure (after every physics time step)                 ?

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define USE_ENT
#define HYCOM_RESOLUTION_1deg
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_stratF40                        ! horiz/vert resolution, 2x2.5, top at 0.1mb, 40 layers
MODEL_COM GEOM_B IORSF              ! model variables and geometry
MODELE                              ! Main and model overhead
ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV !GHY                 ! land surface and soils
ENT_DRV ENT_COM
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB_E1                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO                    ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_F                          ! diagnostics (resolution dependent)
FFT144 OFFT288E
POUT                                ! post-processing output
SparseCommunicator_mod              ! sparse gather/scatter modules

hycom_arrays|-r8| hycom_dim|-r8| kprf_arrays|-r8|
kprf_arrays_loc_renamer|-r8| hycom_atm|-r8|
hycom_arrays_glob|-r8| hycom_arrays_glob_renamer|-r8|
hycom_scalars|-r8| hycom_dim_glob|-r8|
hycom |-r8| OCEAN_hycom|-r8|        ! ocean model - driver
advfct|-r8|                         ! advection
archyb|-r8|                         ! continuity eqn.
barotp|-r8|                         ! barotropic eqn.
bigrid|-r8|                         ! basin grid
blkprf|-r8|                         ! block data
cnuity|-r8|                         ! continuity eqn.
convec|-r8|                         ! convection
cpler |-r8|                         ! coupler
diapfx|-r8|                         ! diapycnal diffusion
dpthuv|-r8| dpudpv|-r8|             ! off-center depth
eice  |-r8|                         ! ice forming
geopar|-r8|                         ! geography related parameters
hybgn1|-r8|                         ! grid generator
inicon|-r8| inigis|-r8| inikpp|-r8| ! initial conditions
matinv|-r8| mxkprf|-r8| mxlayr|-r8| ! mixing scheme
momtum|-r8|                         ! momemtum Eqn.
prtetc|-r8|                         ! print routines, etc.
reflux|-r8|                         ! flux conversion
sigetc|-r8|                         ! eqn.of state, etc.
thermf|-r8|                         ! thermal forcing
trcadv|-r8|                         ! tracer advection
tsadvc|-r8| advem|-r8|              ! advecting t/s

Components:
Ent shared ESMF_Interface solvers giss_LSM

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB
OPTS_giss_LSM = USE_ENT=YES

Data input files:
AIC=AIC.RES_F40.D771201  ! observed init cond (atm. only) ISTART=2
GIC=GIC.144X90.DEC01.1.ext   ! initial ground conditions      ISTART=2
CDN=CD144X90.ext
VEG=V144X90_no_crops.ext
CROPS=CROPS2007_144X90N_nocasp
SOIL=S144X900098M.ext
TOPO=Z144X90N.1deghycom_1
REG=REG2X2.5                      ! special regions-diag
RVR=RD_modelE_Fa.RVR_1deghycom_1.bin
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=LWTables33k.1a              ! rad.tables and history files
RADN4=LWTables33k.1b              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2
! updated aerosols need MADAER=3
TAero_SUL=SUL_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_SSA=SSA_Koch2008_kg_m2_72x46x20h
TAero_NIT=NIT_Bauer2008_kg_m2_72x46x20_1890-2000h
TAero_OCA=OCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCA=BCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCB=BCB_Koch2008_kg_m2_72x46x20_1890-2000h
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust_mass_CakmurMillerJGR06_72x46x20x7x12
RADN7=STRATAER.VOL.1850-1999.Apr02_hdr
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr       ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
! ozone files (minimum 1, maximum 9 files + 1 trend file)
O3file_01=mar2004_o3_shindelltrop_72x46x49x12_1850
O3file_02=mar2004_o3_shindelltrop_72x46x49x12_1890
O3file_03=mar2004_o3_shindelltrop_72x46x49x12_1910
O3file_04=mar2004_o3_shindelltrop_72x46x49x12_1930
O3file_05=mar2004_o3_shindelltrop_72x46x49x12_1950
O3file_06=mar2004_o3_shindelltrop_72x46x49x12_1960
O3file_07=mar2004_o3_shindelltrop_72x46x49x12_1970
O3file_08=mar2005_o3_shindelltrop_72x46x49x12_1980
O3file_09=mar2005_o3_shindelltrop_72x46x49x12_1990
O3trend=mar2005_o3timetrend_46x49x2412_1850_2050
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
BC_dep=BC.Dry+Wet.depositions.ann
TOP_INDEX=top_index_144x90_a.ij.ext
ZVAR=ZVAR2X25A             ! topographic variation for gwdrag
MSU_wts=MSU.RSS.weights.data
GLMELT=GLMELT_144X90_gas.OCN   ! glacial melt distribution
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_nmaps_2x2.5bin.dat

latlonij=latlon387x360.4bin          ! lat & lon at each i,j
hycomtopo=depth387x360.4bin_1        ! topography used in ocean model: with baltic sea
temp_ini=temp387x360x26jan_hv_z1.txt ! 3-d temperature as initial condition
salt_ini=salt387x360x26jan_hv_z1.txt ! 3-d salinity as initial condition
pout_ini=pout387x360x26jan_hv_z1.txt ! 3-d layer pressure as initial condition
ibasin=ibasin387x360.txt_1  ! basin mask
flxa2o=flxa2o387x360.8bin_1 ! coupler weights for flux from atm to ocean
taua2o=taua2o387x360.8bin_1 ! coupler weights for vector from atm to ocean
ssto2a=ssto2a387x360.8bin_1 ! coupler weights for sst from ocean to atm
e_o2a=e_o2a387x360.8bin_1   ! coupler weights for eastward vel from ocean to atm
n_o2a=n_o2a387x360.8bin_1   ! coupler weights for northward vel from ocean to atm
cososino=cososino387x360.8bin           ! cos/sin of i,j axis angle on ocean grid
kpar=seawifs_kpar_387x360.tbin          ! monthly/annual seawifs_kpar data


Label and Namelist:
E4hyf0d (ModelE 2x2.5, 40 lyrs + 1deg hycom + pkpp_min(3/N,.3);[30:200])

DTFIX=180.
&&PARAMETERS
! parameters set for prescribed ocean runs:
KOCEAN=1        ! ocn is prognostic
variable_lk=1
init_flake=1

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:24 25 26 27 28 29 30 31 32 33   34 35 36 37 38 39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000045 !the default is 15d-6
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=10.     ! Shear drag coefficient
CMTN=0.2       ! default is 0.5
CDEF=1.5       ! deformation drag coefficient
XCDNST=400.,10000.   ! strat. gw drag parameters
QGWMTN=1 ! mountain waves ON
QGWDEF=1 ! deformation waves ON
QGWSHR=0 ! shear drag OFF
QGWCNV=0 ! convective drag OFF

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)


U00a=0.73    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.55    ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.57      ! tune this first to get: glob. ann. mean plan.alb=30%   (U00ice up=>albedo down)
U00wtrX=1.46    ! this to get: glob. ann. mean net heat at surf. = 0   (U00wtrX+.01=>NetHtSrf+.7)

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2
madaer=3    ! updated aerosols
aer_rad_forc=0
cloud_rad_forc=1

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1850  ! if -1, crops in VEG-file is used
s0_yr=1850
s0_day=182
ghg_yr=1850
ghg_day=182
volc_yr=-1
volc_day=182
aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1850
dalbsnX=.024
o3_yr=-1850

variable_orb_par=0
orb_par_year_bp=100 !   (which is 1850 AD)

! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

DTsrc=1800.     ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=225.
NIsurf=1        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=480
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags 0hrly
KCOPY=2         ! saving acc + rsf
isccp_diags=0   ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
Nssw=2   ! until diurnal diags are fixed, Nssw has to be even
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

itest=-1        ! default is -1
jtest=-1        ! default is -1
iocnmx=2        ! default is 0
brntop=30.      ! default is 0.
brnbot=200.     ! default is 300.
diapyn=3.e-7    ! default is 3.e-7
diapyc=.3e-4    ! default is 1.e-4
jerlv0=1
&&END_PARAMETERS

 &INPUTZ
   YEARI=1800,MONTHI=1,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default) or earlier
   YEARE=1801,MONTHE=1,DATEE=1,HOURE=0, KDIAG=12*0,9,
   ISTART=2,IRANDI=0, YEARE=1800,MONTHE=1,DATEE=1,HOURE=1,
 &END
