E1F40o32.R GISS Model E  coupled version          larissa   04/06/2009

E1F40o32: 2x2.5x40 layers modelE version, 1850 atm.; 32 layers in the ocean
          NOTE: new ocean initial condition OIC=OIC.WOA98.2HX2.L32.D1201

E1F40o13: replace this section by a description of what distinguishes this run ?
          Use as many lines as you need. Look carefully at all the possible    ?
          choices, particularly the lines containing '?'. In some cases, you   ?
          will have to pick the appropriate choice to make this rundeck work   ?
          The final rundeck should contain no '?'
          Check and modify the rest of the description below:                  ?
modelE1 (3.0) 4x5 hor. grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)       ?
atmospheric composition from year 1880 ? 1979                               ?
ocean: coupled to GISS ocean model (Russell - Schmidt)                      ?
uses turbulence scheme (no dry conv), simple strat.drag (no grav.wave drag) ?
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs  ?
filters: U,V in E-W direction (after every dynamics time step)              ?
         sea level pressure (after every physics time step)                 ?

Run Options
STACKSIZE=524288

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define CHECK_OCEAN                 ! needed to compile aux/file CMPE002
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_F40                             ! horiz/vert resolution, 2x2.5deg, 40 layers -> .1mb
RES_2Hx2_L32                        ! ocean horiz res 2x2.5deg, 32 vert layers
MODEL_COM GEOM_B IORSF              ! model variables and geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
                                    ! parameter database
              ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
ODIAG_COM OCEAN_COM OSTRAITS_F_COM OGEOM  ! dynamic ocean modules
OCNDYN OCNDYN2                            ! dynamic ocean routines
OCN_Interp OCN_Int_LATLON                 ! dynamic ocean routines
OSTRAITS OCNGM OCNKPP                     ! dynamic ocean routines
OCEANR_DIM AFLUXES OFLUXES
ODIAG_PRT                              ! ocean diagnostic print out
OCNFUNTAB                              ! ocean function look up table
SNOW_DRV SNOW                          ! snow model
RAD_COM RAD_DRV RADIATION              ! radiation modules
RAD_UTILS ALBEDO                       ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT          ! diagnostics
DIAG_ZONAL GCDIAGb                     ! grid-dependent code for lat-circle diags
DIAG_RES_F                             ! diagnostics (resolution dependent)
      FFT144 OFFT144E                  ! utilities
POUT                                   ! post-processing output
SparseCommunicator_mod                 ! sparse gather/scatter module

Components:
ESMF_Interface shared

Data input files:
AIC=AIC.RES_F40.D771201         ! observed init cond (atm. only) ISTART=2
GIC=GIC.144X90.DEC01.1.ext      ! initial ground conditions      ISTART=2
OIC=OIC.E2HX2.L32.D1201         ! Levitus ocean intial conditions
TOPO=Z2HX2fromZ1QX1N            ! surface fractions and topography
TOPO_OC=Z2HX2fromZ1QX1N         ! ocean fraction and topography
OFTAB=OFTABLE_NEW               ! ocean function table
AVR=OPF.E2HX2.L32               ! ocean filter
KBASIN=KB144X90.modelE          ! ocean basin designations
CDN=CD144X90.ext                ! neutral drag coefficient
VEG=V144X90_no_crops.ext        ! vegatation file
CROPS=CROPS2007_144X90N_nocasp  ! crops
SOIL=S144X900098M.ext           ! soil properties
REG=REG2X2.5                    ! special regions-diag
RVR=RD_modelE_F.RVR.bin         ! river direction file
RADN1=sgpgxg.table8             ! rad.tables and history files
RADN2=LWTables33k.1a            ! rad.tables and history files
RADN4=LWTables33k.1b            ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
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
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr     ! need KSOLAR=2
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
MSU_wts=MSU.RSS.weights.data
GLMELT=GLMELT_144X90_gas.OCN   ! glacial melt distribution

Label and Namelist:
E1F40o32 (32 ocean layers; 1850 atm.,the current modelE version)

DTFIX=180
&&PARAMETERS
! parameters set for coupled ocean runs:
KOCEAN=1        ! ocn is prognostic
variable_lk=1
init_flake=1

! parameters usually not changed when switching to coupled ocean:

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38  39  40
vsdragl=0.021,0.041,0.077,0.125,0.22,0.275,0.276,0.447,0.96,0.92,  0.91,1.22,1.53,0.3,0.6,0.83, 1.
ANG_SDRAG=1         ! conserve ang. mom.

OBottom_drag=1      !  Drags at the ocean bottom (NO drags -> OBottom_drag=0)
OCoastal_drag=1     !  Drags at the ocean coasts (NO drags -> OCoastal_drag=0)

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)


U00a=0.72   ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.40   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.60       ! tune this first to get: glob. ann. mean plan.alb=30%   (U00ice up=>albedo down)
U00wtrX=1.47     ! this to get: glob. ann. mean net heat at surf. = 0   (U00wtrX+.01=>NetHtSrf+.7)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2
madaer=3    ! updated aerosols

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
dalbsnX=.015
o3_yr=-1850

! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

! parameters that may have to be changed in emergencies:
DTsrc=1800.
DT=225.
DTO=225.
NIsurf=2        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=480       ! use =48 except on halem
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time if isccp-diags are not essential
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
&&END_PARAMETERS

 &INPUTZ
   YEARI=1900,MONTHI=12,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
   YEARE=1951,MONTHE=12,DATEE=1,HOURE=0, KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1900,MONTHE=12,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,
 &END
