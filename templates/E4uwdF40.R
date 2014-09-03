E4uwdF40.R GISS Model E  1850 ocn/atm                     tzhou 03/03/2010
E4uwdF40 : modelE with alternative gravity wave drag: Unresolved Wave Drag
modelE1 (3.0) 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean data: prescribed, 1876-1885 climatology
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define USE_ENT
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_F40                             ! horiz/vert resolution, 2x2.5, top at 0.1mb, 40 layers
MODEL_COM GEOM_B IORSF              ! model variables and geometry
MODELE                              ! Main and model overhead
ALLOC_DRV                           ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
! UNRDRAG_COM UNRDRAG UNRDRAG_DRV     ! unresolved wave drag (alternative gravity wave drag)
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV ! + component giss_LSM: land surface and soils
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + component Ent: new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB_E1                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
OCEAN OCNML                         ! ocean modules
!SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO                    ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_F                          ! diagnostics (resolution dependent)
FFT144                              ! utilities
POUT                                ! post-processing output

Components:
Ent shared ESMF_Interface solvers giss_LSM

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB
OPTS_giss_LSM = USE_ENT=YES

Data input files:
    ! resolution dependent files
    ! start up from the restart file of an earlier run ...
! AIC=1....rsfE... ! initial conditions, no GIC needed, use ISTART=8
    ! ... or from observed conditions AIC and model ground data GIC
AIC=AIC.RES_F40.D771201  ! observed init cond (atm. only) ISTART=2
GIC=GIC.144X90.DEC01.1   ! initial ground conditions      ISTART=2
OSST=OST_144x90.1876-1885avg.HadISST1.1    ! prescr. climatological ocean
SICE=SICE_144x90.1876-1885avg.HadISST1.1   ! prescr. climatological sea ice
! For q-flux ocean, replace all the above by the next 2 lines, set KOCEAN=1, ISTART=8
!! AIC=1JAN1961.rsfE4F40.MXL65m        ! end of preliminary run with KOCEAN=0
!! OHT=OTSPEC.E4F40.MXL65m.1956-1960   ! ocean horizontal heat transports
OCNML=Z1O.B144x90                    ! mixed layer depth (not used if KOCEAN=0)
CDN=CD144X90.ext
VEG=V144X90_no_crops.ext
CROPS=CROPS2007_144X90N_nocasp
SOIL=S144X900098M.ext
TOPO=Z144X90N_nocasp              ! bdy.cond
REG=REG2X2.5                      ! special regions-diag
RVR=RD_modelE_Fa.RVR.bin          ! river direction file
TOP_INDEX=top_index_144x90_a.ij.ext
Z4var=Z4var144x89                 ! topographic variances multiplied by four
GLMELT=GLMELT_144X90_gas.OCN   ! glacial melt distribution
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_nmaps_2x2.5bin.dat
    ! resolution independent files
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
GHG=GHG.Mar2009.txt ! use GHG.Jul2009.txt for runs that start before 1850
dH2O=dH2O_by_CH4_monthly
BC_dep=BC.Dry+Wet.depositions.ann
MSU_wts=MSU.RSS.weights.data

Label and Namelist:  (next 2 lines)
E4uwdF40 (E4F40 with alternative gravity wave drag: unresolved wave drag)

&&PARAMETERS
! parameters set for choice of ocean model:
KOCEAN=0        ! ocean is prescribed
!! KOCEAN=1        ! ocean is computed
Kvflxo=0        ! usually set to 1 only during a prescr.ocn run by editing "I"
!  Kvflxo=1     ! saves VFLXO files to prepare for q-flux runs (mkOTSPEC)
ocn_cycl=1      ! =0 if ocean varies from year to year

variable_lk=1   ! variable lakes

USE_UNR_DRAG=1      !if 1 => SDRAG is turned off and unresolved drag is applied.
                    !if 0 => SDRAG is intact and alternative gwd is not employed.

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.73 ! affects clouds above 850mb w/o MC region;  tune this first to get about 30% high cloud
U00b=1.68 ! affects clouds below 850mb and MC regions; tune this last  to get rad.balance
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
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
albsn_yr=1850
dalbsnX=.024
o3_yr=-1850
CO2X=1.

variable_orb_par=0
orb_par_year_bp=100  !  BP i.e. 1950-orb_par_year_bp AD = 1850 AD

! parameters that control the Shapiro filter
DT_XUfilter=180. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=180. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

DTsrc=1800.     ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=180.
NIsurf=1        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=480
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags 0hrly
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
Nssw=2   ! until diurnal diags are fixed, Nssw has to be even
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs
&&END_PARAMETERS

 &INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default) or earlier
 YEARE=1971,MONTHE=1,DATEE=2,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
 &END
