SCMSGPCONT.R GISS Model E      awolf 01/2013   

SCM: RUN GISS Model E as SCM  using one latitude band                   
scm run using SGP Continuous Forcing Data from Jan 1999 - Dec 2001  
modelE1 (3.0) 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)  
atmospheric composition from year 1979
ocean data: prescribed, 1975-1984 climatology
uses turbulence scheme (no dry conv), simple strat.drag (no grav.wave drag) ?
time steps: physics 30 min.; radiation 30 min.
filters: U,V in E-W direction (after every dynamics time step)              ?
         sea level pressure (after every physics time step)                 ?

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define SCM                          ! run as Single Column Model
#define USE_ENT
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
RES_F40                             ! horiz/vert resolution, 2x2.5, top at 0.1mb, 40 layers
DIAG_RES_F                          ! diagnostics 
FFT144                              ! Fast Fourier Transform 

     ! lat-lon grid specific codes
GEOM_B                              ! model geometry 
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                       ! diagn/post-processing output
IORSF                               ! model variables and geometry

     ! GISS dynamics
ATMDYN_SCM MOMEN2ND                 ! replace atmospheric dynamics wth SCM routines
QUS_DRV TQUS_DRV                    ! advection of Q/tracers

MODEL_COM                           ! model variables
MODELE                              ! Main and model overhead
ALLOC_DRV                           ! domain decomposition, allocate global distributed arrays
ATMDYN_SCM_COM                      ! atmospheric dynamics 
ATMDYN_SCM_EXT ATM_UTILS            ! utilities for some atmospheric quantities
SCM_COM SCMDATA_SGP                 ! routines for reading and processing SCM forcings and IC's
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules 
SURFACE  FLUXES                     ! surface calculation and fluxes
GHY_COM GHY_DRV ! + giss_LSM        ! land surface and soils + snow model     
VEG_DRV                             ! vegetation 
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV ENT_COM ! + Ent             ! new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB                               ! turbulence in whole atmosphere 
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DUM                          ! ice dynamics modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO                    ! radiation and albedo
DIAG_COM DIAG DEFACC                ! diagnostics (diag, diag_prt dummies in scm_diag) 
SCM_DIAG_COM SCM_DIAG               ! SCM diagnostics
OCEAN OCNML                         ! ocean modules

Components:
shared ESMF_Interface solvers giss_LSM
dd2d
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB
OPTS_giss_LSM = USE_ENT=YES

Data input files:
AIC=AIC.RES_F40.D771201  ! observed init cond (atm. only) ISTART=2
GIC=GIC.144X90.DEC01.1   ! initial ground conditions      ISTART=2
OSST=OST_144x90.B.1975-1984avg.Hadl1 ! prescr. climatological ocean (1 yr data)
SICE=SICE_144x90.B.1975-1984avg.Hadl1 ! prescr. climatological sea ice
!! q-flux ocean: use the next line instead,       set KOCEAN=1
!! OHT=OTSPEC.E4F40.MXL65m.1956-1960            ! ocean horizontal heat transports
!! OCNML is not used if KOCEAN=0, but needed in and to prepare for q-flux model
OCNML=Z1O.B144x90                               ! mixed layer depth
!! TOPO=Z144X90N_nocasp ! bdy.cond
TOPO=Z2HX2fromZ1QX1N

!!RVR=RD_modelE_F.RVR.bin      ! river direction file
RVR=RD_Fb.RVR.bin          ! river direction file

CDN=CD144X90.ext
VEG=V144X90_no_crops.ext
!!CROPS=CROPS2007_144X90N_nocasp
CROPS=CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp
SOIL=S144X900098M.ext
TOP_INDEX=top_index_144x90_a.ij.ext
ZVAR=ZVAR2X25A             ! topographic variation for gwdrag
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_nmaps_2x2.5bin.dat
GLMELT=GLMELT_144X90_gas.OCN   ! glacial melt distribution
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
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust_mass_CakmurMillerJGR06_72x46x20x7x12
RADN7=STRATAER.VOL.SATO.1850-1999.Apr02_hdr
RADN8=cloud.epsilon4.72x46
RADN9=solar.DBbglean.ann850-2000.uvflux_hdr       ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
! ozone files (minimum 1, maximum 9 files + 1 trend file)
GHG=GHG.Mar2009.txt ! use GHG.Jul2009.txt for runs that start before 1850
dH2O=dH2O_by_CH4_monthly
BC_dep=BC.Dry+Wet.depositions.ann
! updated aerosols need MADAER=3
!!TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
!!TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
!!TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
!!TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
TAero_SUL=SUL_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_SSA=SSA_Koch2008_kg_m2_72x46x20h
TAero_NIT=NIT_Bauer2008_kg_m2_72x46x20_1890-2000h
TAero_OCA=OCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCA=BCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCB=BCB_Koch2008_kg_m2_72x46x20_1890-2000h
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

MSU_wts=MSU.RSS.weights.data
REG=REG2X2.5          ! special regions-diag

SCMSRF=SGP.surface.0507.dat 
SCMLAY=SGP.layer.0507.dat
SCMcape=psadilookup.dat

Label and Namelist:
SCMSGPCONT (ModelE 2x2.5, 40 lyrs, 1979 atm/ocn;
up to 60 (or 52) columns here to describe your run)?<--col 53  to  72-->to 80-->
DTFIX=180.
&&PARAMETERS
! Target Coordinates for SCM
I_TARG=34        !Southern Great Plains 
J_TARG=64

! parameters set for prescribed ocean runs:
KOCEAN=0        ! ocn is prescribed
Kvflxo=0        ! use =1 to save VFLXO daily ONLY to prepare for q-flux runs
ocn_cycl=1      ! =0 if ocean varies from year to year

variable_lk=0   ! let lakes grow or shrink in horizontal extent

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38   39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.


xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

U00a=.60    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
wmui_multiplier=2.0
entrainment_cont1=.4

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1979  ! if -1, crops in VEG-file is used
s0_yr=1979
s0_day=182
ghg_yr=1979
ghg_day=182
volc_yr=1979
volc_day=182
aero_yr=1979
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1979
dalbsnX=.024
o3_yr=-1979

variable_orb_par=0
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols



DTsrc=1800.     ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=1        ! increase as layer 1 gets thinner
NRAD=1          ! Radiation called every dynamics time step
! parameters that affect at most diagn. output:
aer_rad_forc=0    ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1  ! use=0 to save cpu time (diagnostics doubls radiation calls 
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
Ndisk=4800 
&&END_PARAMETERS

 &INPUTZ
   YEARI=2005,MONTHI=7,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default) or earlier
   YEARE=2005,MONTHE=7,DATEE=31,HOURE=23,     KDIAG=12*0,9,
   ISTART=2,IRANDI=0, YEARE=2005,MONTHE=7,DATEE=3,HOURE=23
 &END
