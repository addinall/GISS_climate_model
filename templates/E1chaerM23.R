E1chaerM23.R GISS Model E                                gas 06/00

modelE1 (3.0) with:
- Drew Shindell tropospheric and stratospheric chemistry (based on Greg's E001TdsSM23.R)
- Dorothy Koch and Kostas Tsigaridis aerosols (based on Dorothy's E1TaerM20.R)

 Set imAER=3, aer_int_yr = 1890 to 2000 for historic 
 Set imAER=5, aer_int_yr = 1850 or 2000 for AR5
 Set imAER=1 for AEROCOM emissions (but this hasn't been tested lately)
 Set emission input files to match (see below)
 Using _E1 'slush' model since it is tested
modelE 4x5 hor. grid with 23 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850 (or 1979)      (look below for "_yr")  ?
ocean data: prescribed, 1876-1885 (or 1975-1984) climatology  (see OSST/SICE) ?
uses turbulence scheme, simple strat.drag (not grav.wave drag)
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs
filters:    U,V in E-W direction (after every dynamics time step)
            sea level pressure (after every physics time step)

************************************************************
Please note this rundeck is provided to show how to run with
chemistry and aerosol tracers on. You should check on the GCM tuning and
boundary conditions before using in production simulations.
************************************************************

You might need to increase your stacksize to get it to run in
parallel. Try 'unlimit stacksize' or increase to the proper amount.

Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define TRACERS_WATER               ! tracers can interact with water
#define TRACERS_DRYDEP              ! include tracer dry deposition
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
#define TRACERS_TERP                ! include terpenes in gas-phase chemistry
#define SHINDELL_STRAT_CHEM         ! turns on stratospheric chemistry
#define WATER_MISC_GRND_CH4_SRC ! adds lake, ocean, misc. ground sources for CH4
!  OFF #define CALCULATE_FLAMMABILITY  ! activated code to determine flammability of surface veg
!  OFF #define CALCULATE_LIGHTNING ! turn on Colin Price lightning when TRACERS_SPECIAL_Shindell off
!  OFF #define SHINDELL_STRAT_EXTRA     ! non-chemistry stratospheric tracers
!  OFF #define INITIAL_GHG_SETUP        ! only for setup hour to get ghg IC file
#define TRACERS_AEROSOLS_Koch    ! Dorothy Koch's tracers (aerosols, etc)
#define TRACERS_NITRATE          ! Nitrate aerosol
#define TRACERS_DUST             ! Dust aerosol
#define TRACERS_AEROSOLS_SOA     ! Secondary Organic Aerosols
#define BC_ALB                   ! Optional tracer BC affects snow albedo
!#define CLD_AER_CDNC             ! Aerosol-cloud interactions
!#define BLK_2MOM                 ! Aerosol-cloud interactions
!  OFF #define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
!  OFF #define NUDGE_ON                 ! nudge the meteorology
!  OFF #define GFED_3D_BIOMASS          ! turns on IIASA AR4 GFED biomass burning
!  OFF #define BIOGENIC_EMISSIONS       ! turns on interactive isoprene emissions
!  OFF #define SULF_ONLY_AEROSOLS        ! when using Koch aerosols,omit BC,OC,SS
!  OFF #define HTAP_LIKE_DIAGS    ! adds many diags, changes OH diag, adds Air tracer
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M23                             ! horiz/vert resolution
MODEL_COM GEOM_B IORSF              ! model variables and geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
STRATDYN STRAT_DIAG                 ! strospheric dynamics (incl. gw drag)
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
TRACER_COM TRACERS_DRV              ! common and driver for tracers
TRDRYDEP                            ! tracer dry deposition from Harvard CTM
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
! ---TRACER SPECIFIC CODES----------
TRACERS_SPECIAL_Shindell            ! routines specific to drew's 15-tracers
TRCHEM_Shindell_COM                 ! Drew Shindell's tracers common
TRCHEM_calc                         ! chemical reaction calculations
TRCHEM_init                         ! chemistry initialization, I/O
TRCHEM_family                       ! tracer family chemistry
!TRCHEM_fastj ! troposphere-only version of tracer chem photlysis code/rad transf
TRCHEM_fastj2                       ! used for trop+strat chem version
TRCHEM_master                       ! trop chem "driver"/strat prescrioption
TRACERS_AEROSOLS_Koch_e4
TRACER_NITRATE                      ! Nitrate aerosol
TRAMP_eqsam_v03d                    ! EQSAM module for inorganic aerosol thermodynamic equilibrium
TRDUST_COM TRDUST TRDUST_DRV        ! DUST
TRACERS_AEROSOLS_SOA                ! Secondary Organic Aerosols
! COSMO_SOURCES
! BIOGENIC_EMISSIONS                  ! old N.Unger interactive isoprene emissions
! ----------------------------------
CLOUDS2_E1 CLOUDS2_DRV CLOUDS_COM   ! clouds modules
!CLD_AEROSOLS_Menon_BLK_MAT BLK_DRV  ! Aerosol-cloud interactions
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils
VEG_DRV                             ! both offline and online vegetation need this
VEG_COM VEGETATION                  ! offline vegetation
PBL_COM PBL_DRV PBL_E1              ! atmospheric pbl
lightning                           ! Colin Price lightning model
!! flammability_drv flammability       ! Olga's fire model
! pick exactly one of the next 2 choices: ATURB or DRYCNV
ATURB_E1                            ! turbulence in whole atmosphere
! DRYCNV                            ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN ICEDYN_DRV                   ! dynamic ice modules
SparseCommunicator_mod              ! sparse gather/scatter module
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV_E1 RADIATION_E1        ! radiation modules
RAD_UTILS ALBEDO                    ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_M                          ! diagnostics (resolution dependent)
      FFT72                         ! utilities
! NUDGE                             ! S. Bauer's code to nudge meteorology
POUT_netcdf                         ! post-processing output

Components:
ESMF_Interface shared

Data input files:
    ! start up from restart file of earlier run
! AIC=1JAN2000.rsfE1chaerM23_terpenes_nitrate_dust        ! initial conditions (atm./ground), no GIC, ISTART=8
    ! or start up from observed conditions
AIC=AIC.RES_M23.D771201          ! initial conditions (atm.)      needs GIC, ISTART=2
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
ZVAR=ZVAR4X5         ! topographic variation for gwdrag
!                                             (end of section 2 of data input files)
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=radfil33k      ! 8/2003 version
RADN3=miescatpar.abcdv2
MSU_wts=MSU.RSS.weights.data
GLMELT=GLMELT_4X5.OCN   ! glacial melt distribution
TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust_mass_CakmurMillerJGR06_72x46x20x7x12
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr    ! need KSOLAR=2
! RADN9=Solar_spectrum.1500-2004_gsf ! need KSOLAR=9
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
delta_O3=do3_shindell_72x46x49x12_2037-2107_E009TdsSxHp-SnHpM23
GHG=GHG.Mar2004.txt
! e.g. GHGic=GHG_IC_1999
dH2O=dH2O_by_CH4_monthly
TOP_INDEX=top_index_72x46_a.ij.ext ! only used if #define DO_TOPMODEL_RUNOFF
BC_dep=BC.Dry+Wet.depositions.ann
!-----------------------------------------------
! choose these for trop-only chem model
!-----------------------------------------------
! Note in the paramters, also set fix_CH4_chemistry=0 please
!  MOLEC=chem_files/ds3ch4_moleculesE_terp_soa
!  JPLRX=chem_files/jpl00_T15trop_SEP08_terp
!  JPLPH=chem_files/ds_photlist_trop_15
!  RATJ=chem_files/ratj.giss_15
!  SPECFJ=chem_files/jv_spec00_15.dat
!-----------------------------------------------
! choose these for strat+trop chem model
!-----------------------------------------------
MOLEC=chem_files/ds4_moleculesE_terp_soa
JPLRX=chem_files/jpl00_T25_SEP08_fastterp
JPLPH=chem_files/ds4_photlist_T25
RATJ=chem_files/ratj.giss_25
SPECFJ=chem_files/jv_spec_AV.dat
N2O_IC=gsin/N2O_IC_M23_4x5_6.17_conc
CFC_IC=gsin/CFC_IC_M23_4x5_6.17_conc
CH4_IC=gsin/CH4_IC_M23_4x5_6.17_conc
!-----------------------------------------------
!------- Needed for dry deposition ---------
VEG_DENSE=gsin/veg_dense_4x5
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
Ox_IC=gsin/Ox_init_cond_M23_4x5_conc !see README in /usr/people/cmrun/gsin
CO_IC=gsin/CO_init_cond_M23_conc
Ox_ref=gsin/O3ref_O3JDAY_1850_182.dat
! fltran file used if rad_FL.ne.0:
! FLTRAN=chem_files/Solar_spectrum.1500-2004_fastj2 ! KSOLAR=9
! FLTRAN=chem_files/solar.lean02.ann.uvflux_hdr_fastj2  ! KSOLAR=2
!----------Default emissions case (mostly AR5 late apr 09)------------
CO_01=AR5_emis/M/NOV09/2000/CO_agr_AR5_2000_4x5_h
CO_02=AR5_emis/M/NOV09/2000/CO_awb_AR5_2000_4x5_h
CO_03=AR5_emis/M/NOV09/2000/CO_dom_AR5_2000_4x5_h
CO_04=AR5_emis/M/NOV09/2000/CO_ene_AR5_2000_4x5_h
CO_05=AR5_emis/M/NOV09/2000/CO_ind_AR5_2000_4x5_h
CO_06=AR5_emis/M/NOV09/2000/m_CO_shp_AR5_2000_4x5_h
CO_07=AR5_emis/M/NOV09/2000/CO_slv_AR5_2000_4x5_h
CO_08=AR5_emis/M/NOV09/2000/CO_tra_AR5_2000_4x5_h
CO_09=AR5_emis/M/NOV09/2000/CO_wst_AR5_2000_4x5_h
CO_10=AR5_emis/M/NOV09/2000/CO_forestfire_AR5_2000_4x5_h
CO_11=AR5_emis/M/NOV09/2000/CO_grassfire_AR5_2000_4x5_h
! ========= please remember that Alkenes =================
! ========= and Paraffin emissions files =================
! ========= must now be in Kmole units,  =================
! ========= not Kg units ...             =================
Alkenes_01=AR5_emis/M/NOV09/2000/Alkenes_agr_AR5_2000_4x5_h
Alkenes_02=AR5_emis/M/NOV09/2000/Alkenes_awb_AR5_2000_4x5_h
Alkenes_03=AR5_emis/M/NOV09/2000/Alkenes_dom_AR5_2000_4x5_h
Alkenes_04=AR5_emis/M/NOV09/2000/Alkenes_ene_AR5_2000_4x5_h
Alkenes_05=AR5_emis/M/NOV09/2000/Alkenes_ind_AR5_2000_4x5_h
Alkenes_06=AR5_emis/M/NOV09/2000/m_Alkenes_shp_AR5_2000_4x5_h
Alkenes_07=AR5_emis/M/NOV09/2000/Alkenes_slv_AR5_2000_4x5_h
Alkenes_08=AR5_emis/M/NOV09/2000/Alkenes_tra_AR5_2000_4x5_h
Alkenes_09=AR5_emis/M/NOV09/2000/Alkenes_wst_AR5_2000_4x5_h 
Alkenes_10=gsin/Alkenes_GEIA_vegetation_head_1
Alkenes_11=AR5_emis/M/NOV09/2000/Alkenes_forestfire_AR5_2000_4x5_h
Alkenes_12=AR5_emis/M/NOV09/2000/Alkenes_grassfire_AR5_2000_4x5_h
Paraffin_01=AR5_emis/M/NOV09/2000/Paraffin_agr_AR5_2000_4x5_h
Paraffin_02=AR5_emis/M/NOV09/2000/Paraffin_awb_AR5_2000_4x5_h
Paraffin_03=AR5_emis/M/NOV09/2000/Paraffin_dom_AR5_2000_4x5_h
Paraffin_04=AR5_emis/M/NOV09/2000/Paraffin_ene_AR5_2000_4x5_h
Paraffin_05=AR5_emis/M/NOV09/2000/Paraffin_ind_AR5_2000_4x5_h
Paraffin_06=AR5_emis/M/NOV09/2000/m_Paraffin_shp_AR5_2000_4x5_h
Paraffin_07=AR5_emis/M/NOV09/2000/Paraffin_slv_AR5_2000_4x5_h
Paraffin_08=AR5_emis/M/NOV09/2000/Paraffin_tra_AR5_2000_4x5_h
Paraffin_09=AR5_emis/M/NOV09/2000/Paraffin_wst_AR5_2000_4x5_h
Paraffin_10=gsin/Paraffin_GEIA_vegetation_head_1
Paraffin_11=AR5_emis/M/NOV09/2000/Paraffin_forestfire_AR5_2000_4x5_h
Paraffin_12=AR5_emis/M/NOV09/2000/Paraffin_grassfire_AR5_2000_4x5_h
NOx_01=AR5_emis/M/NOV09/2000/NOx_agr_AR5_2000_4x5_h
NOx_02=AR5_emis/M/NOV09/2000/NOx_awb_AR5_2000_4x5_h
NOx_03=AR5_emis/M/NOV09/2000/NOx_dom_AR5_2000_4x5_h
NOx_04=AR5_emis/M/NOV09/2000/NOx_ene_AR5_2000_4x5_h
NOx_05=AR5_emis/M/NOV09/2000/NOx_ind_AR5_2000_4x5_h
NOx_06=AR5_emis/M/NOV09/2000/m_NOx_shp_AR5_2000_4x5_h
NOx_07=AR5_emis/M/NOV09/2000/NOx_tra_AR5_2000_4x5_h
NOx_08=AR5_emis/M/NOV09/2000/NOx_wst_AR5_2000_4x5_h
NOx_09=NOy_sources/NOx_GEIA_soil_half_head
NOx_10=AR5_emis/M/NOV09/2000/NOx_forestfire_AR5_2000_4x5_h
NOx_11=AR5_emis/M/NOV09/2000/NOx_grassfire_AR5_2000_4x5_h
NOx_AIRC=AR5_emis/M/NOV09/2000/NOx_air_AR5_2000_4x5
CH4_01=AR5_emis/M/NOV09/2000/CH4_agr_AR5_2000_4x5_h
CH4_02=AR5_emis/M/NOV09/2000/CH4_awb_AR5_2000_4x5_h
CH4_03=AR5_emis/M/NOV09/2000/CH4_dom_AR5_2000_4x5_h
CH4_04=AR5_emis/M/NOV09/2000/CH4_ene_AR5_2000_4x5_h
CH4_05=AR5_emis/M/NOV09/2000/CH4_ind_AR5_2000_4x5_h
CH4_06=AR5_emis/M/NOV09/2000/m_CH4_shp_AR5_2000_4x5_h
CH4_07=AR5_emis/M/NOV09/2000/CH4_tra_AR5_2000_4x5_h
CH4_08=AR5_emis/M/NOV09/2000/CH4_wst_AR5_2000_4x5_h
CH4_09=methane/gcm_data/CH4_GEIA_Soil_Absorption_header
CH4_10=methane/gcm_data/CH4_GEIA_Termites_header
CH4_11=AR5_emis/M/2000/CH4_forestfire_AR5_2000_4x5_h
CH4_12=AR5_emis/M/2000/CH4_grassfire_AR5_2000_4x5_h
CH4_13=methane/gcm_data/CH4_GEIA_Wetlands_and_Tundra_header
! if BIOGENIC_EMISSIONS or PS_BVOC are defined, and only one
! Terpenes file is available, the Isoprene file is needed too
Isoprene_01=ORCHIDEE_Isoprene_1990_4x5_h
Terpenes_01=ORCHIDEE_Terpenes_1990_4x5_h
Terpenes_02=ORCHIDEE_ORVOC_1990_4x5_h
SULFATE_SA=NOy_sinks/sulfate_fakeM23_M_SA
DMS_FIELD=dms_conc
SO2_FIELD=so2_conc
! ----- for interactive wetlands -----
PREC_NCEP=gsin/mean_prec_ncep_4x5
TEMP_NCEP=gsin/mean_temp_ncep_4x5
BETA_NCEP=gsin/beta_p_ch4_4x5
ALPHA_NCEP=gsin/alpha_t_ch4_4x5
!-----------------------------------------------
! ------- Dorothy's inputs needed for imAER= 1 or 3 -----------
SO2_VOLCANO=SO2_volc_conti2000.AEROCOM_FEB12
! ------- Dorothy's inputs needed for imAER= 1 or 3 -----------
!SO2_BIOMASS=SO2_bio2000.AEROCOM_DEC03
! ------- Dorothy's inputs needed for imAER= 3 -----------
!BC_SHIPS=Shipping_BC_2000_4x5.bin
!POM_SHIPS=Shipping_POM_2000_4x5.bin
!SO2_SHIPS=Shipping_SO2_2000_4x5.bin
! --------Dorothy's inputs needed for imAER = 3 or 5 --------------
DMS_SEA=DMS_Night_4x5
! AER_CHEM: inputs for sulfur and H2O2 chemistry, saved from the
! Shindell chemistry runs: OH, HO2, photolysis rate of H2O2 and NO3
! Only needed if COUPLED_CHEM = 0
!AER_CHEM=Sulf_chem_drewE_20      ! 20 layer version
AER_CHEM=trace_gas_3D_fields_E   ! 23 layer version
AER_OH_STRAT=Strat_OH_drewE_20
! The following line is only needed for imAER=1
! Note that if gas-phase chemistry is enabled, both
! TERPENE and Terpenes_01 are used, which differ significantly
!TERPENE=terp_Guenther_4x5
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

! files for dust tracers
VTRSH=fake_144x90_dust_wind_speed_thresholds ! modelIIprime varying thresholds (obsolete)
FRCLAY=fake_144x90_dust_FrClay ! modelIIprime soil clay fractin (obsolete)
FRSILT=fake_144x90_dust_FrSilt ! modelIIprime model silt fraction (obsolete)
DRYHR=fake_144x90_dust_DryHrPminusE ! emission only if E exceeds P (obsolete)
ERS=ERS1_1993_MONTHLY.72x46.threshold-13 ! ERS data
DSRC=Ginoux2001_source_VegMask_72x46     ! preferred dust sources
                                         ! optimized use: prefDustSources=0
                                         ! (only choice so far)
LKTAB=log_dust_emission_60ms-1 ! look up table for emission calculations
LKTAB1=table_wspdf             ! look up table for wind speed probabilities


Label and Namelist:
E1chaerM23 (sample rundeck with Shindell chemistry tracers)
R=00BG/B
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

H2ObyCH4=0.     ! (=1) activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 no land icea fixup, A.Lacis orig. 6-band alb)
KSOLAR=2 ! 9

!--- define emission sectors above files belong to ---
NOx_AIRC_sect='AIR' ! special 3D source case
CH4_13_sect='WETL'
!      (careful; they're allowed to overlap):
!       ---------define-REGIONS------------
!        global S.Asia E.Asia Europe N.Amer
REG_S=    -90.,    5.,   15.,   25.,   15.
REG_N=     90.,   35.,   50.,   65.,   55.
REG_W=   -180.,   50.,   95.,  -10., -125.
REG_E=    180.,   95.,  160.,   50.,  -60.
!       ---define-regions-names/order------
REGIONS_ARE='global S_Asia E_Asia Europe N_America'
!-fit-here--|                                                              |---
!       ---define-factors-by-sector--------
!        global S.Asia E.Asia Europe N.Amer
SECT_01= 1.000, 1.000, 1.000, 1.000, 1.000 ! WETL
SECT_02= 1.000, 1.000, 1.000, 1.000, 1.000 ! AIR 
!       ---define-sectors-names/order------
SECTORS_ARE='WETL AIR'
!-fit-here--|
!-----
aircraft_Tyr1=0 ! for non-transient emissions,
aircraft_Tyr2=0 ! set these two equal or omit them.
biomass_Tyr1= 0 ! for non-transient emissions,
biomass_Tyr2= 0 ! set these two equal or omit them.

! factor to tune the base isoprene emissions globally,
! only when #defined BIOGENIC_EMISSIONS, otherwise use
! a sector/region method above...
base_isopreneX=1.d0

! Colin Price lightning model needs resolution-dependant tuning:
tune_lt_land=4.774d0 ! =2.2d0*2.17d0 for 4x5 model setting
tune_lt_sea=8.463d0  ! =3.9d0*2.17d0 for 4x5 model setting

! ---- for interactive wetlands -----
nn_or_zon=1     ! int dist method 1=zonal avg, 0=nearest neighbor
int_wet_dist=1  ! turn on(1)/off(0) interacive SPATIAL wetlands
ice_age=0.      ! if not 0 no wetl emis for lats poleward of +/- this in deg
ns_wet=13       ! index of CH4 source that is the wetlands (dumb, I know)
exclude_us_eu=0 ! to exclude (=1) the U.S. and E.U. from inter wetl dist
topo_lim=205.d0 ! upper limit of topographic variation for new wetlands 
sat_lim=-9.d0   ! lower limit on surf air temp for new wetlants
gw_ulim=100.d0  ! upper limit on ground wetness for new wetlands
gw_llim=18.d0   ! lower limit on ground wetness for new wetlands
SW_lim=27.d0    ! lower limit on SW downward flux for new wetlands
! -----------------------------------
X_SDRAG=.0005,.00005  ! used for lin. sdrag above P_SDRAG mb
C_SDRAG=0.     ! no constant sdrag
P_SDRAG=.01     ! lin. sdrag above p_sdrag mb (top layer for M23) except near poles
PP_SDRAG=4.6   ! lin. sdrag above  pp_sdrag mb near poles (top 5 layers for M23)
ANG_SDRAG=1    ! if =1: sdrag conserves ang mom.
WMAX=300.     ! maximum wind velocity in sdrag; default=200 when GW drag not used. For ISTART=2, use WMAX=300.
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000030 !the default is 15d-6
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=1.      ! Shear drag coefficient
CMTN=0.25      ! default is 0.5
CDEF=1.5       ! deformation drag coefficient
LMCM=16              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters

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

DTsrc=1800.        ! half-hour physics time step (default: DTsrc=3600.). Cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=450.            ! dynamics time step
NIsurf=1           ! number of surface time steps. Increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=480       ! e.g. every 10 days
! pick sub-daily frequency diags, next line ::
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
SUBDD1=' '
LmaxSUBDD=18         ! only need to save up to level LmaxSUBDD
KCOPY=2              ! saving acc + rsf
isccp_diags=1        ! use =0 to save cpu time, but you lose some key diagnostics
cloud_rad_forc=0   ! use =1 to activate this diagnostic (doubles radiation calls !)
nda5d=13           ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13           ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
nssw=2          ! until diurnal diagn. are fixed, nssw should be even
  
!--------- Aerosol parameters----------------
imAER=5         !3 historic; 1 AEROCOM ; 0,2 for standard or sector inputs (not working)
aer_int_yr=2000    !used for imAER=3,5 select desired year (1890 to 2000) or 0 to use JYEAR
rad_interact_aer=1 ! 1=couples aerosols to radiation, 0=use climatology
                   ! (either case does the rad-forcing calculation)
rad_forc_lev=1     ! 0 for TOA, 1 for tropopause for rad forcing diags.
                   ! use LTROPO(I,J) level for rad forcing diags.
diag_rad=0         ! 1: additional radiation diagnostics
diag_wetdep=0      ! 1: additional wet deposition diagnostics

!--------- dust aerosol parameters----------------
imDust=0              ! 0: PDF emission scheme, 1: AEROCOM
adiurn_dust=0         ! 1: daily dust diagnostics at certain grid points
prefDustSources=0     ! 0: Ginoux 2001 w/ vegetation mask
                      ! 1-4: No files and optimization available yet
                      ! >4: Free choice of emis. parameters
!fracClayPDFscheme=1. ! Frac. clay emis, only effective for prefDustSources > 4
!fracSiltPDFscheme=1. ! Frac. silt emis, only effective for prefDustSources > 4
                      ! set internally for 0-4

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
COUPLED_CHEM=1     ! to couple chemistry and aerosols
use_sol_Ox_cycle=0 ! (=1) apply ozone changes in radiation, based on solar cycle
clim_interact_chem=1! 1=use calculated Ox/CH4 in radiation, 0=use climatology
                   ! (either case does the rad-forcing calculation)
                   ! 0 also turns off chem(H2O)-->Q( ) feedback
rad_forc_lev=1     ! use LTROPO(I,J) level for rad forcing diags.
use_rad_n2o=0      ! use the radiation code's N2O 
use_rad_cfc=0      ! use rad code cfc11+cfc12, adjusted 
use_rad_ch4=0      ! use rad code CH4, shut off sfc sources 
rad_FL=0           ! use rad code insolation getting fastj2 photon flux
prather_limits=1   ! to avoid some negative tracers in sub-gridscale
which_trop=1       ! choose tropopause for chemistry purposes:
                   ! 0=LTROPO(I,J), 1=LS1-1
fix_CH4_chemistry=-1   ! for setting fixed methane value for chemistry:
pfix_CH4_S=1.750d-6    ! Southern Hemisphere (fix_CH4_chemistry=1)
pfix_CH4_N=1.855d-6    ! Northern Hemisphere (fix_CH4_chemistry=1)
ch4_init_sh=1.750      ! init cond, S.Hemi (fix_CH4_chemistry=0 ?)
ch4_init_nh=1.855      ! init cond, N.Hemi (fix_CH4_chemistry=0 ?)
scale_ch4_IC_file=1.d0 ! multiplicative factor on CH4 IC file (fix_CH4_chemistry=-1)

aircrafts_Tyr1=1940 ! for aircraft emissions. For non-transient emissions,
aircrafts_Tyr2=2000 !          set these two equal or omit them.

! To run a preindustrial case, alter the sulfate surface area, SST, & seaice
! files and the radiation years in this rundeck. Also, fix the ch4 chemistry
! to the appropriate value (0.73ppmv in both hemispheres?) Then, set PI_run=1
! and choose the various ratios for altering initial conditions, stratospheric
! overwriting. Probably best to alter sources via the sectors set in the
! rundeck instead, use these for initial conditions/overwriting:
PI_run        =    0    ! 0 =no, 1=yes for running pre-industrial cases
PIratio_N     = 0.667d0 ! {NOx, HNO3, N2O5, HO2NO2}
PIratio_CO_T  = 0.333d0 ! CO in troposphere
PIratio_CO_S  = 0.500d0 ! CO in stratosphere
PIratio_other = 0.500d0 ! {PAN,Isoprene,AlkyNit,Alkenes,Paraffin}
PIratio_N2O   = 0.896d0 ! {N2O IC's and L=1 overwriting}
PIratio_CFC   = 0.000d0 ! {CFC IC's and L=1 overwriting}
!----------------------------------------------------------------------
!----------------------------------------------------------------------
&&END_PARAMETERS

 &INPUTZ
   QCHECK=.false.
   kdiag = 0,0,0,0,0,0,0,0,0,0,0,0,0,
   YEARI=1999,MONTHI=12,DATEI=1,HOURI=0,
   YEARE=2000,MONTHE=1,DATEE=1,HOURE=0,
   ISTART=2,IRANDI=0, YEARE=1999, MONTHE=12,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,
 &END

