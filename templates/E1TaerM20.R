E1TaerM20.R GISS Model E  2004 modelE             Dorothy Koch  11/08

E1TaerM20: Sample rundeck with interactive mass aerosols
  including sulfate, nitrate, sea salt, black and organic carbon
 Set imAER=3, aer_int_yr = 1890 to 2000 for historic 
 Set imAER=5, aer_int_yr = 1850 or 2000 for AR5
 Set imAER=1 for AEROCOM emissions (but this hasn't been tested lately)
 Set emission input files to match (see below)
 Using _E1 'slush' model since it is tested
modelE 4x5 hor. grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850 (or 1979)      (look below for "_yr")  ?
ocean data: prescribed, 1876-1885 (or 1975-1984) climatology  (see OSST/SICE) ?
uses turbulence scheme, simple strat.drag (not grav.wave drag)
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs
filters:    U,V in E-W direction (after every dynamics time step)
            sea level pressure (after every physics time step)

Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define TRACERS_WATER               ! tracers can interact with water
#define TRACERS_DRYDEP              ! include tracer dry deposition
#define TRACERS_AEROSOLS_Koch       ! Dorothy Koch's tracers (aerosols, etc)
#define BC_ALB                      !optional tracer BC affects snow albedo
#define TRACERS_NITRATE
!#define BIN_TRACERS
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M20AT DIAG_RES_M FFT72          ! horiz/vert resolution, 4x5deg, 20 layers -> .1mb
MODEL_COM GEOM_B IORSF              ! model variables and geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
ATM_UTILS
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
TRACER_COM TRACERS_DRV              ! common and driver for tracers
TRDRYDEP                            ! tracer dry deposition from Harvard CTM
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
CLOUDS2_E1 CLOUDS2_DRV CLOUDS_COM   ! clouds modules
!CLOUDS2 CLOUDS2_DRV CLOUDS_COM   ! clouds modules
! ---TRACER SPECIFIC CODES----------
TRACERS_AEROSOLS_Koch_e4            ! BC/OC/sulfate/seasalt
TRACER_NITRATE                      ! Nitrate aerosol
TRAMP_eqsam_v03d                    ! EQSAM module for inorganic aerosol thermodynamic equilibrium
! ---------------------------------
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL_E1                ! atmospheric pbl
!PBL_COM PBL_DRV PBL                ! atmospheric pbl
ATURB_E1                               ! turbulence in whole atmosphere
!ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV_E1 RADIATION_E1           ! radiation modules
!RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO                    ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
                                    ! utilities
POUT                                ! post-processing output

Components:
ESMF_Interface shared

Data input files:
    ! start up from restart file of earlier run
! AIC=1DECxxxx.rsfEyyyy           ! initial conditions (atm./ground), no GIC, ISTART=8
    ! or start up from observed conditions
AIC=AIC.RES_M20A.D771201          ! initial conditions (atm.)      needs GIC, ISTART=2
GIC=GIC.E046D3M20A.1DEC1955.ext   ! initial conditions (ground)
    ! ocean data for "prescribed ocean" runs : climatological ocean
OSST=OST4X5.B.1876-85avg.Hadl1.1  ! prescr. climatological ocean (1 yr of data)
SICE=SICE4X5.B.1876-85avg.Hadl1.1 ! prescr. climatological sea ice
!? for 1979 OSST=OST4X5.B.1975-84avg.Hadl1.1
!? for 1979 SICE=SICE4X5.B.1975-84avg.Hadl1.1
OCNML=Z1O.B4X5.cor                ! mixed layer depth (needed for post processing)
!                                             (end of section 1 of data input files)
    ! resolution dependent files
TOPO=Z72X46N.cor4_nocasp SOIL=S4X50093.ext ! soil/topography bdy.conds
! VEG=V72X46.1.cor2   ! or:       ! vegetation fractions  (sum=1), need crops_yr=-1
VEG=V72X46.1.cor2_no_crops.ext CROPS=CROPS2007_72X46N.cor4_nocasp  ! veg. fractions, crops history
CDN=CD4X500S.ext                  ! surf.drag coefficient
REG=REG4X5                        ! special regions-diag
RVR=RD_modelE_M.RVR.bin               ! river direction file
TOP_INDEX=top_index_72x46_a.ij.ext  ! only used if #define DO_TOPMODEL_RUNOFF
!                                             (end of section 2 of data input files)
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=radfil33k      ! 8/2003 version
RADN3=miescatpar.abcdv2
!RADN2=LWTables33k.1a              ! rad.tables and history files
!RADN4=LWTables33k.1b              ! rad.tables and history files
RADN5=H2Ocont_Ma_2000             ! H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
!  MADAER=1 (default) needs:
TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
! MADAER=3 needs: (temporary version)
!TAero_SUL=SUL_Koch2008_kg_m2_72x46x20_1890-2000h
!TAero_SSA=SSA_Koch2008_kg_m2_72x46x20h
!TAero_NIT=NIT_Bauer2008_kg_m2_72x46x20_1890-2000h
!TAero_OCA=OCA_Koch2008_kg_m2_72x46x20_1890-2000h
!TAero_BCA=BCA_Koch2008_kg_m2_72x46x20_1890-2000h
!TAero_BCB=BCB_Koch2008_kg_m2_72x46x20_1890-2000h
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust_mass_CakmurMillerJGR06_72x46x20x7x12
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr      ! need KSOLAR=2
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
MSU_wts=MSU.RSS.weights.data
GLMELT=GLMELT_4X5.OCN   ! glacial melt distribution
!------- Needed for dry deposition ---------
ATMFJ=chem_files/jv_atms.dat
DRYCOEFF=chem_files/drydep.coef
VEGTYPE=chem_files/vegtype.global
VEGTYPEBIN=chem_files/vegtype.global4X5.bin
OLSON=chem_files/drydep.table
LAI01=chem_files/lai01.global
LAI02=chem_files/lai02.global
LAI03=chem_files/lai03.global
LAI04=chem_files/lai04.global
LAI05=chem_files/lai05.global
LAI06=chem_files/lai06.global
LAI07=chem_files/lai07.global
LAI08=chem_files/lai08.global
LAI09=chem_files/lai09.global
LAI10=chem_files/lai10.global
LAI11=chem_files/lai11.global
LAI12=chem_files/lai12.global
Ox_IC=gsin/Ox_init_cond_M23_4x5 !see README in /usr/people/cmrun/gsin
! ------- Dorothy's inputs needed for imAER= 1,3,5 -----------
SO2_VOLCANO=SO2_volc_conti2000.AEROCOM_FEB12
! ------- Dorothy's inputs needed for imAER= 1 or 3 -----------
!SO2_BIOMASS=SO2_bio2000.AEROCOM_DEC03
! ------- Dorothy's inputs needed for imAER= 3 -----------
!BC_SHIPS=Shipping_BC_2000_4x5.bin
!POM_SHIPS=Shipping_POM_2000_4x5.bin
!SO2_SHIPS=Shipping_SO2_2000_4x5.bin
! --------Dorothy's inputs needed for imAER = 3 or 5 --------------
DMS_SEA=DMS_Night_4x5
AER_CHEM=Sulf_chem_drewE_20
AER_OH_STRAT=Strat_OH_drewE_20
! The following line is only needed for imAER=1
! Note that if gas-phase chemistry is enabled, both
! TERPENE and Terpenes_01 are used, which differ significantly
!TERPENE=terp_Guenther_4x5
Terpenes_01=ORCHIDEE_Terpenes_1990_4x5_h
! ------Dorothy's inputs needed for imAER = 3 -old historic emissions------------
!SO2_EM_1=SO2_EDGAR_Feb09_4x5_h_1890-2000 
!OC_EM_1=OC_Bond_Feb09_4x5_h_1850-2000 !BC/OC Bond  
!BC_EM_1=BC_Bond_Feb09_4x5_h_1850-2000
!BC_BIOMASS=BC_GFED_97-06_4x5
!OC_BIOMASS=OC_GFED_97-06_4x5
!SO2_AIRCRAFT=NOy_sources/aircraft_4x5_1940-2000 ! zero in 1940 and before.
!--------Dorothy's inputs for imAER=5 ---AR5 emissions-----------------
!   for now all years in this block must be set to either 1850 or 2000
BC_EM_1=BCII_awb_AR5_2000_4x5_h
BC_EM_2=BCII_dom_AR5_2000_4x5_h
BC_EM_3=BCII_ene_AR5_2000_4x5_h
BC_EM_4=BCII_ind_AR5_2000_4x5_h
BC_EM_5=BCII_tra_AR5_2000_4x5_h
BC_EM_6=BCII_wst_AR5_2000_4x5_h
BC_EM_7=BCII_shp_AR5_2000_4x5_h
OC_EM_1=OCII_awb_AR5_2000_4x5_h
OC_EM_2=OCII_dom_AR5_2000_4x5_h
OC_EM_3=OCII_ene_AR5_2000_4x5_h
OC_EM_4=OCII_ind_AR5_2000_4x5_h
OC_EM_5=OCII_tra_AR5_2000_4x5_h
OC_EM_6=OCII_wst_AR5_2000_4x5_h
OC_EM_7=OCII_shp_AR5_2000_4x5_h
SO2_EM_1=SO2_shp_AR5_2000_4x5_h
SO2_EM_2=SO2_awb_AR5_2000_4x5_h
SO2_EM_3=SO2_dom_AR5_2000_4x5_h
SO2_EM_4=SO2_ind_AR5_2000_4x5_h
SO2_EM_5=SO2_tra_AR5_2000_4x5_h
SO2_EM_6=SO2_wst_AR5_2000_4x5_h
SO2_EM_7=SO2_ene_AR5_2000_4x5_h
BCB_EM_1=BCB_forestfire_AR5_2000_4x5_h
BCB_EM_2=BCB_grassfire_AR5_2000_4x5_h
OCB_EM_1=OCB_forestfire_AR5_2000_4x5_h
OCB_EM_2=OCB_grassfire_AR5_2000_4x5_h
SO2B_EM_1=SO2_forestfire_AR5_2000_4x5_h
SO2B_EM_2=SO2_grassfire_AR5_2000_4x5_h
! -------Dorothy's inputs for imAER=1 AeroCom (not tested lately)---------
!SO2_VOLCANO_EXP=SO2_volc_expl1750.AEROCOM
!TERPENE=SOA_2000.AEROCOM_DEC03
!BC_BIOMASS=BC_fire2000.AEROCOM_DEC03
!OC_BIOMASS=OC_fire2000.AEROCOM_DEC03
!BC_BIOFUEL=BC_bf2000.AEROCOM_DEC03
!BC_FOSSIL_FUEL=BC_ff2000.AEROCOM_DEC03
!OC_BIOFUEL=OC_bf2000.AEROCOM_DEC03
!OC_FOSSIL_FUEL=OC_ff2000.AEROCOM_DEC03
!SO2_IND=SO2_ind2000.AEROCOM_DEC03
!SALT1=SALT_bin1_2000_new.nc
!SALT2=SALT_bin2_2000_new.nc
!-----------------------------------------------
!  AEROSOL INPUT NITRATE oxidants, fields use for imAER=3,5
O3_FIELD=Ox_3D_field_bell
OFFLINE_HNO3.nc=HNO3_E70_GISS4x5.nc
OFFLINE_SEAS.nc=SEASALT_EK1su_GISS4x5.nc
!------- for imAER=3 ---------
!NH3_CON_01=NH3hCON_IIASA_CLE_Apr09_4x5_h_1890-2030
!NH3_CYC_01=NH3hCYC_IIASA_CLE_Apr09_4x5_h_1890-2030
!------- for imAER=5 -----set date to 1850 or 2000-----
NH3_CON_01=NH3hCON_OCEANflux_Jan10_4x5_h
NH3_CON_02=NH3_forestfire_AR5_2000_4x5_h
NH3_CON_03=NH3_grassfire_AR5_2000_4x5_h
NH3_CON_04=NH3_agr_AR5_2000_4x5_h
NH3_CON_05=NH3_awb_AR5_2000_4x5_h
NH3_CON_06=NH3_dom_AR5_2000_4x5_h
NH3_CON_07=NH3_ind_AR5_2000_4x5_h
NH3_CON_08=NH3_ene_AR5_2000_4x5_h
NH3_CON_09=NH3_tra_AR5_2000_4x5_h


Label and Namelist:
E1TaerM20 (ModelE1 4x5, 20 lyrs, with interactive aerosol tracers) 
DTFIX=300

&&PARAMETERS
! parameters set for prescribed ocean runs:
KOCEAN=0 ! 0 or 1 , use =0 if ocn is prescribed, use =1 if ocn is predicted
Kvflxo=0 ! use 1 ONLY to save VFLXO daily to prepare for q-flux run ?
ocn_cycl=1  ! ? use =0 if prescribed ocean varies from year to year

variable_lk=1 ! let lakes grow or shrink in horizontal extent
wsn_max=2.   ! restrict snow depth to 2 m-h2o (if 0. snow depth is NOT restricted)

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

! tuning param.: this setting works for 1850; use U00wtrX=1.28 for 1979

!     if CLOUDS2_E1 is replaced by CLOUDS2, use: 
! U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
! U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
!     instead of:
U00ice=.62      ! U00ice+.01 =>dBal=1.5,dPl.alb=-.9%   goals:Bal=0,plan.alb=30%
U00wtrX=1.29    ! U00wtrX+.01=>dBal=0.7,dPl.alb=-.25%  Bal=glb.ann NetHt at z0
! HRMAX=500.    ! not needed unless do_blU00=1, HRMAX up => nethtz0 down (alb up)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=-1 no land ice fixup, 1 Lacis' scheme)
KSOLAR=2

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr= -1  ! if -1, crops in VEG-file is used   ! =1979 , also change OSST,SICE
s0_yr=1990                                         ! =1979 , also change OSST,SICE
s0_day=182
ghg_yr=1990                                        ! =1979 , also change OSST,SICE
ghg_day=182
volc_yr=-1  ! 1850-1999 mean strat.aeros           ! =1979 , also change OSST,SICE
volc_day=182
aero_yr=-1850                                       ! =-1979 , also change OSST,SICE
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1850                                      ! =1979 , also change OSST,SICE
dalbsnX=.024
o3_yr=-1850                                        ! =-1979 , also change OSST,SICE

! parameters that control the Shapiro filter
DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

DTsrc=1800.     ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=450.
NIsurf=1        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=480
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
KCOPY=2         ! saving acc + rsf  ? =3 to also save "oda"-files
isccp_diags=1   ! use =0 to save cpu time, but you lose some key diagnostics
cloud_rad_forc=0 ! use =1 to activate this diagnostic (doubles radiation calls !)
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
nssw=2          ! until diurnal diagn. are fixed, nssw should be even
!MADAER=3
  
!--------- Aerosol parameters----------------
imAER=5         !3 historic; 1 AEROCOM ; 0,2 for standard or sector inputs (not working)
aer_int_yr=2000    !used for imAER=3,5 select desired year (1890 to 2000) or 0 to use JYEAR
rad_interact_aer=1 ! 1=couples aerosols to radiation, 0=use climatology
rad_forc_lev=1     ! 0 for TOA, 1 for tropopause for rad forcing diags.
prather_limits=1   ! to avoid some negative tracers in sub-gridscale

aircrafts_Tyr1=1940 ! for aircraft emissions. For non-transient emissions,
aircrafts_Tyr2=2000 !          set these two equal or omit them.
!----------------------------------------------------------------------

&&END_PARAMETERS

 &INPUTZ
   QCHECK=.false.
   kdiag = 0,0,0,0,0,0,0,0,0,0,0,0,0,
   YEARI=1998,MONTHI=12,DATEI=1,HOURI=0, IYEAR1=1998 ! or earlier
   YEARE=2001,MONTHE=1,DATEE=1,HOURE=0,     
   ISTART=2,IRANDI=0, YEARE=1998,MONTHE=12,DATEE=1,HOURE=1,
 &END
