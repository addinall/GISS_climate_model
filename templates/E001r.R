E001r.R GISS Model E  strat.H2O added inst.frc          rar 2/01/02

E001r: modelE1 (3.0) radiation only, instantaneous forcing

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M12                             ! horiz/vert resolution, 4x5deg, 12 layers -> 10mb
MODEL_COM GEOM_B IORSF              ! model variables and geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
                                    ! parameter database
              ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS2 CLOUDS2_DRV CLOUDS_COM        ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices: ATURB or DRYCNV
! ATURB                             ! turbulence in whole atmosphere
DRYCNV                              ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN  ! or: ICEDYN_DUM ! dynamic sea ice modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO                    ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_M                          ! diagnostics (resolution dependent)
      FFT72                         ! utilities
POUT                                ! post-processing output

Components:
ESMF_Interface shared

Data input files:
AIC=1JAN1951.rsfE001.frc ! just 3 label records + radia-record (use fcop1)
OSST=OST4X5.B.1975-84avg.Hadl1.1  ! prescr. climatological ocean (1 yr of data)
SICE=SICE4X5.B.1975-84avg.Hadl1.1 ! prescr. climatological sea ice
CDN=CD4X500S    ! surf.drag coefficient
! VEG=V72X46.1.cor2   ! or:       ! vegetation fractions  (sum=1), need crops_yr=-1
VEG=V72X46.1.cor2_no_crops CROPS=CROPS2007_72X46N.cor4_nocasp  ! veg. fractions, crops history
SOIL=S4X50093 TOPO=Z72X46N.cor4_nocasp   ! soil/topography bdy.conds
REG=REG4X5                        ! special regions-diag
RVR=RD_modelE_M.RVR.bin                   ! river direction file
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=LWTables33k.1a              ! rad.tables and history files
RADN4=LWTables33k.1b              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_Ma_2008
RADN3=miescatpar.abcdv2
! RADNA,RADNB are no longer used
TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust_mass_CakmurMillerJGR06_72x46x20x7x12
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr      ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
! new ozone files (minimum 1, maximum 9 files)
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
TOP_INDEX=top_index_72x46_a.ij.ext
GLMELT=GLMELT_4X5.OCN   ! glacial melt distribution
RADJAN=E001/RADJAN1950    ! replace E001 by the control run
RADFEB=E001/RADFEB1950    ! assuming that the saved data are
RADMAR=E001/RADMAR1950    ! still in the original directory
RADAPR=E001/RADAPR1950
RADMAY=E001/RADMAY1950
RADJUN=E001/RADJUN1950
RADJUL=E001/RADJUL1950
RADAUG=E001/RADAUG1950
RADSEP=E001/RADSEP1950
RADOCT=E001/RADOCT1950
RADNOV=E001/RADNOV1950
RADDEC=E001/RADDEC1950

Label and Namelist:
E001r (ModelE1 (3.0) inst.forcing run - control)


&&PARAMETERS
KOCEAN=0        ! ocn is prescribed
Kvflxo=0        ! don't touch this line
ocn_cycl=1      ! =0 if ocean varies from year to year

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1979 ! if -1, crops in VEG-file is used
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
dalbsnX=.015
o3_yr=-1979

! parameters that control the Shapiro filter
DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

         ! most other parameters are irrelevant
NSUBDD=0            ! don't touch this line
KCOPY=2             ! saving acc + rsf
Kradia=1            ! use Kradia=2 for adj. forcing run
&&END_PARAMETERS

 &INPUTZ
   YEARE=1951,MONTHE=12,DATEE=1,HOURE=0,  ! assumed start: 12/1/1950
   ISTART=8, YEARE=1950,MONTHE=12,DATEE=1,HOURE=6,IWRITE=1,JWRITE=1,ITWRITE=23,
 &END

It's important that DTsrc, NIrad, ItimeI are left as in the control
run so the input data are available when they are needed ! That's
why in this case ISTART=8 (start) looks more like ISTART=9 (restart).
The only parameters/input files that matter are the ones that affect
the radiation. It is their change whose instantaneous/adjusted
forcing is computed.
