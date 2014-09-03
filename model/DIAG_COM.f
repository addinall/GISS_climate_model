#include "rundeck_opts.h"
      MODULE DIAG_COM
!@sum  DIAG_COM Diagnostic model variables
!@auth Original Development Team
!@ver  2010/11/12
      use resolution, only : ls1
      USE MODEL_COM, only : im,jm,lm,ntype,kep,istrat,lm_req
      use diag_zonal, only : jm_budg,imlonh,jmlat,xwon
      use socpbl, only : npbl=>n
#ifdef NEW_IO
      use cdl_mod
#endif
      IMPLICIT NONE
      SAVE
      private

      public LM_REQ,im,jm,lm,imlonh,ntype,istrat,kep,jm_budg,jmlat,xwon

C**** Accumulating_period information
      INTEGER, DIMENSION(12), public :: MONACC  !@var MONACC(1)=#Januaries, etc
      CHARACTER*12, public :: ACC_PERIOD='PARTIAL'    !@var string MONyyr1-yyr2
!@var AMON0,JMON0,JDATE0,JYEAR0,JHOUR0,Itime0  beg.of acc-period

!!  WARNING: if new diagnostics are added, change io_diags/reset_DIAG !!
C**** ACCUMULATING DIAGNOSTIC ARRAYS

!@var LAT_BUDG latitudes of budget grid
!@var DXYP_BUDG area array of budget grid
      REAL*8, DIMENSION(JM_BUDG), public :: LAT_BUDG,DXYP_BUDG

!@var local area array of budget grid (locally owned by 1 processor but has global size)
      REAL*8, DIMENSION(JM_BUDG), public :: DXYP_BUDG_LOC
!@var area weight for diagnostics on budget grid
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: WTBUDG,WTBUDG2

!@param KAJ number of accumulated zonal budget diagnostics
      INTEGER, PARAMETER, public :: KAJ=85
#ifdef HEALY_LM_DIAGS
     &                                 + 3
      REAL*8, public :: CROPS_DIAG(IM,JM)
#endif

!@var AJ zonal budget diagnostics for each surface type
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AJ,AJ_loc
     &     ,AJ_out

!@var SQRTM moved from DIAG5A where it was a saved local array to this
!@var place so its size could be allocated dynamically and still have
!@var it preserved from call to call of DIAG5A
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: SQRTM
!@param NREG number of regions for budget diagnostics
      INTEGER, PARAMETER, public :: NREG=24
!@var AREG regional budget diagnostics
      REAL*8, DIMENSION(NREG,KAJ), public :: AREG,AREG_loc,
     &     AREG_out

!@var TITREG,NAMREG title and names of regions for AREG diagnostics
      CHARACTER*4, public :: TITREG*80,NAMREG(2,23)
!@var JREG lat/lon array defining regions for AREG diagnostics
      INTEGER, ALLOCATABLE, DIMENSION(:,:), public :: JREG
cmax      INTEGER, DIMENSION(IM,JM), public :: JREG
!@var SAREA_REG areas of the special regions
      REAL*8, DIMENSION(NREG), public :: SAREA_REG

!@param KAJL number of AJL diagnostics
      INTEGER, PARAMETER, public :: KAJL=77
!@var AJL latitude/height diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AJL,AJL_loc

!@param KASJL number of ASJL diagnostics
      INTEGER, PARAMETER, public :: KASJL=5
!@var ASJL latitude/height supplementary diagnostics (merge with AJL?)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: ASJL,ASJL_loc

!@param KAIJ number of AIJ diagnostics
      INTEGER, PARAMETER, public :: KAIJ=415
#ifdef ACCMIP_LIKE_DIAGS
     &                                   + 8
#endif
#ifdef HEALY_LM_DIAGS
     &                                   + 1
#endif
C**** RDF 
#ifdef TES_SIM
     &                                  +  26
#endif
!@var AIJ latitude/longitude diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AIJ,AIJ_loc

!@param KAIJL number of AIJL accumulations
      INTEGER, PARAMETER, public :: KAIJL=20
#ifdef CLD_AER_CDNC
     &                                  + 16
#endif
#ifdef HTAP_LIKE_DIAGS
     &                                  +  1
#endif
C**** RDF
#ifdef TES_SIM
     &                                  +  20
#endif
!@var IJL_xxx,IJK_xxx AIJL diagnostic indices
!@+   IJL/IJK refer to model versus constant-pressure levels
      INTEGER, public ::
     &     IJL_DP,IJK_DP,IJL_U,IJL_V,IJK_TX,IJK_Q,
     &     IJL_W,IJK_RH,IJL_RC,IJL_MC
      INTEGER, public :: IJL_CF, IJL_MCamFX, IJL_cldwtr,IJL_cldice
     &    ,IJL_LLH,IJL_MCTLH,IJL_MCDLH,IJL_MCSLH
     &    ,IJL_REWM,IJL_REWS,IJL_CDWM,IJL_CDWS,IJL_CWWM,IJL_CWWS
     &    ,IJL_REIM,IJL_REIS,IJL_CDIM,IJL_CDIS,IJL_CWIM,IJL_CWIS
     &    ,IJL_CFWM,IJL_CFIM,IJL_CFWS,IJL_CFIS
     &    ,IJL_TEMPL,IJL_GRIDH,IJL_HUSL,IJL_ZL
C**** RDF
#ifdef TES_SIM	 
C    for categorical AK selection
     &    ,IJL_DP_QF,IJL_H2O_QF,IJL_HDO_QF
     &    ,IJL_DP_TES,IJL_H2O_TES,IJL_HDO_TES
     &    ,IJL_H2O_TES_OLD,IJL_HDO_TES_OLD
     &    ,IJL_H2O_PRIOR
     &    ,IJL_TES_SEN,IJL_TES_SEN_HT
C    for AK selection from retrievals
     &    ,IJL_DP_QF_R,IJL_H2O_QF_R,IJL_HDO_QF_R
     &    ,IJL_DP_TES_R,IJL_H2O_TES_R,IJL_HDO_TES_R
     &    ,IJL_H2O_TES_OLD_R,IJL_HDO_TES_OLD_R
     &    ,IJL_H2O_PRIOR_R
     &    ,IJL_TES_SEN_R,IJL_TES_SEN_HT_R
#endif

!@var AIJL 3D accumulations for longitude/latitude/level diagnostics
      REAL*8, DIMENSION(:,:,:,:), allocatable, public :: AIJL,AIJL_loc

C NEHIST=(TROPO/L STRAT/M STRAT/U STRAT)X(ZKE/EKE/SEKE/ZPE/EPE)X(SH/NH)
!@param NED number of different energy history diagnostics
!@param NEHIST,HIST_DAYS number of energy history columns,rows (max)
      INTEGER, PARAMETER, public :: NED=10
      INTEGER, PARAMETER, public :: NEHIST=NED*(2+ISTRAT)
      INTEGER, PARAMETER, public :: HIST_DAYS=100
!@var ENERGY energy diagnostics
      REAL*8, DIMENSION(NEHIST,HIST_DAYS), public :: ENERGY

!@var NPTS number of points at which standard conserv. diags are called
      INTEGER, PARAMETER, public :: NPTS = 11
!@param NQUANT Number of conserved quantities in conservation diags
      INTEGER, PARAMETER, public :: NQUANT=24
!@param KCON number of conservation diagnostics
      INTEGER, PARAMETER, public :: KCON=170
!@var CONSRV conservation diagnostics
      REAL*8, DIMENSION(JM_BUDG,KCON), public :: CONSRV_loc,CONSRV
!@var SCALE_CON scales for conservation diagnostics
      REAL*8, DIMENSION(KCON), public :: SCALE_CON
!@var TITLE_CON titles for conservation diagnostics
      CHARACTER*32, DIMENSION(KCON), public :: TITLE_CON
!@var NSUM_CON indices for summation of conservation diagnostics
!@var IA_CON IDACC numbers for conservation diagnostics
      INTEGER, DIMENSION(KCON), public :: NSUM_CON, IA_CON
!@var NOFM indices for CONSRV array
      INTEGER, DIMENSION(NPTS+1,NQUANT), public :: NOFM
!@var icon_xx indexes for conservation quantities
      INTEGER, public ::
     &     icon_AM,icon_KE,icon_MS,icon_TPE,icon_WM,icon_LKM,
     *     icon_LKE,icon_EWM,icon_WTG,icon_HTG,icon_OCE,icon_OMSI,
     *     icon_OHSI,icon_OSSI,icon_LMSI,icon_LHSI,icon_MLI,icon_HLI,
     *     icon_MICB,icon_HICB
!@var KCMX actual number of conservation diagnostics
      INTEGER, public :: KCMX = 25 ! take up first 25 indexes for special cases
!@var CONPT0 default titles for each point where conserv diags. are done
      CHARACTER*10, DIMENSION(NPTS), public :: CONPT0 = (/
     *     "DYNAMICS  ","CONDENSATN","RADIATION ","PRECIPITAT",
     *     "LAND SURFC","SURFACE   ","FILTER    ","OCEAN     ",
     *     "DAILY     ","SRF OCN FL","OCN DYNAM "/)

!@param KSPECA,NSPHER number of spectral diagnostics, and harmonics used
      INTEGER, PARAMETER, public :: KSPECA=20
      INTEGER, PARAMETER, public :: NSPHER=4*(2+ISTRAT)
!@var SPECA spectral diagnostics
      REAL*8, DIMENSION((IMLONH+1),KSPECA,NSPHER), public :: SPECA
!@var KLAYER index for dividing up atmosphere into layers for spec.anal.
      INTEGER, DIMENSION(LM), public :: KLAYER
!@param PSPEC pressure levels at which layers are seperated and defined
C**** 1000 - 150: troposphere           150 - 10 : low strat.
C****   10 - 1: mid strat               1 and up : upp strat.
      REAL*8, DIMENSION(4), PARAMETER, public ::
     &     PSPEC = (/ 150., 10., 1., 0. /)
!@var LSTR level of interface between low and mid strat. (approx 10 mb)
      INTEGER, public :: LSTR = LM   ! defaults to model top.

!@param KTPE number of spectral diagnostics for pot. enthalpy
      INTEGER, PARAMETER, public :: KTPE=8
      integer, parameter, public :: NHEMI=2
!@var ATPE pot. enthalpy spectral diagnostics
      REAL*8, DIMENSION(KTPE,NHEMI), public :: ATPE

!@param HR_IN_DAY hours in day
      INTEGER, PARAMETER, public :: HR_IN_DAY=24
!@param lmax_dd2 most upper layer for which multilayer diurnal diagnostics
!+               is written, currently: up to first constant pressure layer
      integer, parameter, public :: lmax_dd2=ls1
!@param NDIUVAR number of diurnal diagnostics
#ifdef TRACERS_AMP
      INTEGER, PARAMETER, public :: NDIUVAR=73+16+16+100+40+40+40+40
#else
#ifdef TRACERS_DUST
      INTEGER, PARAMETER, public :: NDIUVAR=73+14*lmax_dd2+6*npbl
     &     +4*(npbl-1)
#else
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      INTEGER, PARAMETER, public :: NDIUVAR=62
#else
      INTEGER, PARAMETER, public :: NDIUVAR=56
#endif
#endif
#endif
!@param NDIUPT number of points where diurnal diagnostics are kept
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      INTEGER, PARAMETER, public :: NDIUPT=34
#else
      INTEGER, PARAMETER, public :: NDIUPT=4
#endif
!@dbparam adiurn_dust  flag to switch on/off intra daily diagnostics for dust
!@+                    default=0 (off)
      INTEGER, public :: adiurn_dust=0
!@dbparam save3dAOD flag to calculate ttasv_save in rad code even when
!@+ adiurn_dust is off.
      INTEGER, public :: save3dAOD=0
!@dbparam IJDD,NAMDD (i,j)-coord.,names of boxes w/diurnal cycle diag
!@+       defaults set in DIAG_RES (depends on resolution)
      INTEGER, DIMENSION(2,NDIUPT), public :: IJDD
      CHARACTER*4, DIMENSION(NDIUPT), public :: NAMDD
!@dbparam LLDD (lon,lat)-coords (deg) of boxes w/diurnal cycle diag
!@+       defaults set in init_DIAG
      REAL*8, DIMENSION(2,NDIUPT), public :: LLDD

!@var ADIURN diurnal diagnostics (24 hour cycles at selected points)
      REAL*8, DIMENSION(NDIUVAR,NDIUPT,HR_IN_DAY), public :: ADIURN
     &     ,ADIURN_loc
!@param HR_IN_MONTH hours in month
      INTEGER, PARAMETER, public :: HR_IN_MONTH=HR_IN_DAY*31
#ifndef NO_HDIURN
!@var HDIURN hourly diagnostics (hourly value at selected points)
!@+     Same quantities as ADIURN but not averaged over the month
      REAL*8, DIMENSION(NDIUVAR,NDIUPT,HR_IN_MONTH), public :: HDIURN
     &     ,HDIURN_loc
#endif
!@param KAGC number of latitude-height General Circulation diags
!@param KAGCX number of accumulated+derived GC diagnostics
      INTEGER, PARAMETER, public :: KAGC=82+KEP, KAGCX=KAGC+100
!@var AGC latitude-height General Circulation diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AGC,AGC_loc
     &     ,AGC_out

!@param KAIJK number of lat/lon constant pressure diagnostics
      INTEGER, PARAMETER, public :: KAIJK=15
!@var AIJK lat/lon constant pressure diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:), public :: AIJK,AIJK_loc

!@param NWAV_DAG number of components in spectral diagnostics
      INTEGER, PARAMETER, public :: NWAV_DAG=min(9,imlonh)
!@param Max12HR_sequ,Min12HR_sequ lengths of time series for wave powers
      INTEGER, PARAMETER, public :: Max12HR_sequ=2*31, Min12HR_sequ=2*28
!@param RE_AND_IM complex components of wave power diagnostics
      INTEGER, PARAMETER, public :: RE_AND_IM=2
!@param KWP number of wave power diagnostics
      INTEGER, PARAMETER, public :: KWP=12
!@var WAVE frequency diagnostics (wave power)
      REAL*8, DIMENSION(RE_AND_IM,Max12HR_sequ,NWAV_DAG,KWP), public ::
     &     WAVE

C**** parameters and variables for ISCCP diags
!@param ntau,npress number of ISCCP optical depth,pressure categories
      integer, parameter, public :: ntau=7,npres=7
!@param nisccp number of ISCCP histogram regions
      integer, parameter, public :: nisccp = 5
!@var isccp_press pressure mid points for isccp histogram
      INTEGER, PARAMETER, public ::
     &     isccp_press(npres) = (/ 90, 245, 375, 500,
     *     630, 740, 900 /)
!@var isccp_tau lower bound of optical depth for each isccp tau category
      REAL*8, PARAMETER, public ::
     &     isccp_tau(ntau) = (/ 0d0,.1d0,1.3d0,3.6d0,
     *     9.4d0,23d0,60d0 /)
!@var isccp_taum mid point of optical depth for each isccp tau category
      REAL*8, PARAMETER, public :: isccp_taum(ntau-1) = (/ .05d0,0.7d0,2
     *     .95d0,6.5d0,16.2d0,41.5d0 /)
!@var isccp_late edge latitudes for each isccp lat category (region)
!@var isccp_lat midpoint latitudes for each isccp lat category (region)
! calculation of midpoints is hard-coded until all fortran compilers
! allow array arithmetic in PARAMETERS
      REAL*8, PARAMETER, public ::
     &  isccp_late(nisccp+1)=(/-60,-30,-15,15,30,60/)
     & ,isccp_lat(nisccp)=(/-45.,-22.5,0.,22.5,45./)
!@var AISCCP accumlated array of ISCCP histogram
      real*8, public, dimension(ntau,npres,nisccp) :: AISCCP,AISCCP_loc
!@var WISCCP denominator array for ISCCP histograms
      real*8, public, dimension(nisccp) :: WISCCP

!@param KGZ number of pressure levels for some diags
      INTEGER, PARAMETER, public :: KGZ = 13
!@param kgz_max is the actual number of geopotential heights saved
      INTEGER, public :: kgz_max
!@param PMB pressure levels for geopotential heights (extends to strat)
!@param GHT ~mean geopotential heights at PMB level (extends to strat)
!@param PMNAME strings describing PMB pressure levels
      REAL*8, DIMENSION(KGZ), PARAMETER, public ::
     &     PMB=(/1000d0,850d0,700d0,500d0,300d0,100d0,30d0,10d0,
     *     3.4d0,.7d0,.16d0,.07d0,.03d0/),
     *     GHT=(/0.,1500.,3000.,5600.,9500.,16400.,24000.,30000.,
     *     40000.,50000.,61000.,67000.,72000./)
      CHARACTER*4, DIMENSION(KGZ), PARAMETER, public :: PMNAME=(/
     *     "1000","850 ","700 ","500 ","300 ","100 ","30  ","10  ",
     *     "3.4 ","0.7 ",".16 ",".07 ",".03 " /)
#ifdef TRACERS_SPECIAL_Shindell
!@var O_inst saved instantaneous Ox tracer (at PMB lvls)
!@var X_inst saved instantaneous NOx tracer (at PMB lvls)
!@var N_inst saved instantaneous NO2 non-tracer (at PMB lvls)
!@var M_inst saved instantaneous CO tracer (at PMB lvls)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public ::
     & O_inst,N_inst,X_inst,M_inst
#endif
#ifdef TES_LIKE_DIAGS
!@param kgz_max_more is the actual number of TES pressure levels saved
      INTEGER, public :: kgz_max_more
!@param KGZmore more than KGZ number of pressure levels for some diags
      INTEGER, PARAMETER, public :: KGZmore = 32
!@param PMBmore like PMB but more pressure levels for some diags
      REAL*8, DIMENSION(KGZmore), PARAMETER, public ::
     &     PMBmore=(/
     & 1000.0000, 825.40198, 681.29102, 562.34198, 464.16000, 383.11700,
     & 316.22699, 261.01599, 215.44400, 177.82899, 146.77901, 121.15200,
     & 100.00000, 82.540604, 68.129501, 56.233898, 46.415798, 38.311901,
     & 31.622900, 26.101700, 21.544300, 17.782801, 14.678000, 12.115300,
     & 10.000000, 8.2540197, 5.1089802, 3.1622701, 2.1544299, 1.3335201,
     &  0.681292, 0.2154430/)
!@param PMNAMEmore strings describing PMBmore pressure levels
      CHARACTER*4, DIMENSION(KGZmore),PARAMETER,public::PMNAMEmore=(/
     & "1000", "825 ", "681 ", "562 ", "464 ", "383 ",
     & "316 ", "261 ", "215 ", "178 ", "147 ", "121 ",
     & "100 ", "82.5", "68.1", "56.2", "46.4", "38.3",
     & "31.6", "26.1", "21.5", "17.8", "14.7", "12.1",
     & "10.0", "8.25", "5.11", "3.16", "2.15", "1.33",
     & "0.68", "0.22" /)
!@var Q_more saved instantaneous specific hum (at PMBmore lvls)
!@var T_more saved instantaneous temperature(at PMBmore lvls)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: Q_more,T_more
#ifdef TRACERS_SPECIAL_Shindell
!@var O_more saved instantaneous Ox tracer (at PMBmore lvls)
!@var X_more saved instantaneous NOx tracer (at PMBmore lvls)
!@var N_more saved instantaneous NO2 non-tracer (at PMBmore lvls)
!@var M_more saved instantaneous CO tracer (at PMBmore lvls)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public ::
     & O_more,N_more,X_more,M_more
#endif
#endif /* TES_LIKE_DIAGS */
!@var save{H,M,L}CLDI,saveCTPI,saveTAUI,save{S,T}CLDI: SUBDDiag
!@+  instantaneous save arrays for ISCCP cloud variables
!@var saveMCCLDTP instnt.SUBDD moist convective cloud top pressure
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: saveHCLDI,
     & saveMCLDI,saveLCLDI,saveCTPI,saveTAUI,saveSCLDI,saveTCLDI,
     & saveMCCLDTP
C**** RDF: vars for ISCCP subgrid tau and ctp
#ifdef TES_SIM
!@var saveBOXTAU,saveBOXCTP ISCCP subgrid tau and ctp
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: saveBOXTAU,
     & saveBOXCTP
#endif
C**** Instantaneous constant pressure level fields
!@var Z_inst saved instantaneous height field (at PMB levels)
!@var RH_inst saved instantaneous relative hum (at PMB levels)
!@var T_inst saved instantaneous temperature(at PMB levels)
!@var P_acc accumulated precip (special for SUBDD)
!@var PM_acc accumulated moist convective precip (special for SUBDD)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public ::
     &     Z_inst,RH_inst,T_inst
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: P_acc,PM_acc
!@var R_acc accumulated total runoff (special for SUBDD)
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: R_acc

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:), public :: AFLX_ST

!@param KTSF number of freezing temperature diagnostics
      integer, parameter, public :: ktsf=4
!@var TSFREZ freezing temperature diagnostics
C****   1  FIRST DAY OF GROWING SEASON (JULIAN DAY)
C****   2  LAST DAY OF GROWING SEASON (JULIAN DAY)
C****   3  LAST DAY OF ICE-FREE LAKE (JULIAN DAY)
C****   4  LAST DAY OF ICED-UP LAKE  (JULIAN DAY)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: TSFREZ,TSFREZ_loc

!@param KTD number of diurnal range diagnostics (incl. temperature)
      INTEGER, PARAMETER, public :: KTD=12
!@var TDIURN diurnal range temperature diagnostics
C****   1  MIN TG1 OVER EARTH FOR CURRENT DAY (C)
C****   2  MAX TG1 OVER EARTH FOR CURRENT DAY (C)
C****   3  MIN TS OVER EARTH FOR CURRENT DAY (K)
C****   4  MAX TS OVER EARTH FOR CURRENT DAY (K)
C****   5  SUM OF COMPOSITE TS OVER TIME FOR CURRENT DAY (C)
C****   6  MAX COMPOSITE TS FOR CURRENT DAY (K)
C****   7  MAX TG1 OVER OCEAN ICE FOR CURRENT DAY (C)
C****   8  MAX TG1 OVER LAND ICE FOR CURRENT DAY (C)
C****   9  MIN COMPOSITE TS FOR CURRENT DAY (K)
C****  10  MIN RH SURF FOR CURRENT DAY (%)
C****  11  MAX RH SURF FOR CURRENT DAY (%)
C****  12  MAX SURF WIND SPEED FOR CURRENT DAY (M/S)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: TDIURN
      REAL*8,allocatable, public :: TDIURN_glob(:,:,:)

!@nlparam KDIAG array of flags to control diagnostics printout
      INTEGER, DIMENSION(13), public :: KDIAG

!@param NKEYNR number of key number diagnostics
      INTEGER, PARAMETER, public :: NKEYNR=43
!@param NKEYMO number of months key diagnostics are saved
      INTEGER, PARAMETER, public :: NKEYMO=50
!@var KEYNR time-series of key numbers
      INTEGER, DIMENSION(NKEYNR,NKEYMO), public :: KEYNR = 0
!@var KEYCT next index in KEYNR to be used (1->nkeymo)
      INTEGER, public :: KEYCT = 1

!@nlparam IWRITE,JWRITE,ITWRITE control rad.debug output (i,j,amount)
      INTEGER, public :: IWRITE = 0, JWRITE = 0, ITWRITE = 23
!@nlparam QDIAG TRUE for outputting binary diagnostics
      LOGICAL, public :: QDIAG = .FALSE.
!@nlparam QDIAG_ratios TRUE for forming ratios if title="q1 x q2"
      LOGICAL, public :: QDIAG_ratios = .TRUE.

!@var OA generic diagnostic array for ocean heat transport calculations
C****
C****       DATA SAVED IN ORDER TO CALCULATE OCEAN TRANSPORTS
C****
C****       1  ACE1I+SNOWOI  (INSTANTANEOUS AT NOON GMT)
C****       2  MSI2   (INSTANTANEOUS AT NOON GMT)
C****       3  HSIT   (INSTANTANEOUS AT NOON GMT)
C****       4  ENRGP  (INTEGRATED OVER THE DAY)
C****       5  SRHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       6  TRHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       7  SHDT   (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       8  EVHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       9  TRHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      10  SHDT   (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      11  EVHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      12  SRHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****
C**** Extra array needed for dealing with advected ice
C****      13  HCHSI  (HORIZ CONV SEA ICE ENRG, INTEGRATED OVER THE DAY)
C****
!@param KOA number of diagnostics needed for ocean heat transp. calcs
      INTEGER, PARAMETER, public :: KOA = 13
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: OA
      ! REAL*8 :: OA_glob    (IM, JM, KOA)
      REAL*8,allocatable, public :: OA_glob(:,:,:)

C****
C**** Information about acc-arrays:
C****      names, indices, units, idacc-numbers, etc.

!@var iparm/dparm int/double global parameters written to acc-file
      integer, parameter, public :: niparm_max=100
      character(len=20), dimension(niparm_max), public :: iparm_name
      integer, dimension(niparm_max), public :: iparm
      integer, public :: niparm=0
      integer, parameter, public :: ndparm_max=100
      character(len=20), dimension(ndparm_max), public :: dparm_name
      REAL*8, dimension(ndparm_max), public :: dparm
      integer, public :: ndparm=0

      integer, parameter, public ::
     &     sname_strlen=30,units_strlen=30,lname_strlen=80

!@var J_xxx zonal J diagnostic names
      INTEGER, public ::
     &     J_SRINCP0, J_SRNFP0, J_SRNFP1, J_SRABS, J_SRINCG,
     *     J_SRNFG, J_TRNFP0, J_TRNFP1, J_TRHDT, J_RNFP0, J_RNFP1,
     *     J_RHDT, J_SHDT, J_EVHDT, J_HZ1, J_TG2, J_TG1, J_EVAP,
     *     J_PRCP, J_TX, J_TX1, J_TSRF, J_DTSGST, J_DTDGTR, J_RICST,
     *     J_RICTR, J_ROSST, J_ROSTR, J_RSI, J_TYPE, J_RSNOW,
     *     J_OHT, J_DTDJS, J_DTDJT, J_LSTR, J_LTRO, J_EPRCP,
     *     J_RUN, J_ERUN, J_HZ0, J_H2OCH4, J_LWCORR,
     *     J_RVRD,J_ERVR,J_IMELT, J_HMELT, J_SMELT,J_IMPLM, J_IMPLH,
     *     J_WTR1,J_ACE1, J_WTR2,J_ACE2, J_SNOW, J_BRTEMP, J_HZ2,
     *     J_PCLDSS,J_PCLDMC, J_PCLD,J_CTOPP, J_PRCPSS, J_PRCPMC, J_QP,
     *     J_GAM,J_GAMM, J_GAMC,J_TRINCG, J_FTHERM, J_HSURF, J_HATM,
     *     J_PLAVIS,J_PLANIR,J_ALBVIS, J_ALBNIR, J_SRRVIS, J_SRRNIR,
     *     J_SRAVIS,J_SRANIR,J_CLDDEP, J_CLRTOA, J_CLRTRP, J_TOTTRP,
     *     J_ALBP0,J_ALBG,J_IRGW,J_IRGWE
#ifdef HEALY_LM_DIAGS
     *     ,J_VTAU,J_CROPS,J_GHG
#endif
!@var NAME_J,UNITS_J Names/Units of zonal J diagnostics
      character(len=sname_strlen), dimension(kaj), public :: name_j
      character(len=units_strlen), dimension(kaj), public :: units_j
!@var LNAME_J Long names of zonal J diagnostics
      character(len=lname_strlen), dimension(kaj), public :: lname_j
!@var STITLE_J short titles for print out for zonal J diagnostics
      character(len=16), dimension(kaj), public :: stitle_j
!@var SCALE_J scale for zonal J diagnostics
      real*8, dimension(kaj), public :: scale_j
!@var IA_J IDACC indexes for zonal J diagnostics
      integer, dimension(kaj), public :: ia_j
!@var FMT_J Format strings for zonal J diagnostics
      character(len=30), dimension(kaj), public :: fmt_j,fmt_reg
!@var iden_j denominators for zonal J diagnostics
      integer, dimension(kaj), public :: iden_j,iden_reg
!@var HEMIS_J hemispheric/global averages of AJ
      real*8, dimension(:,:,:), allocatable, public :: hemis_j

      character(len=sname_strlen), dimension(kaj), public :: name_reg

!@var IJ_xxx AIJ diagnostic names
      INTEGER, public ::
     &     IJ_RSOI, IJ_RSNW, IJ_SNOW, IJ_SHDT, IJ_PREC, IJ_EVAP,
     *     IJ_SSAT, IJ_BETA,  IJ_SLP1,  IJ_P4UV, IJ_PRES, IJ_PHI1K,
     *     IJ_PHI850, IJ_PHI700, IJ_PHI500, IJ_PHI300, IJ_PHI100,
     *     IJ_PHI30, IJ_PHI10, IJ_PHI3p4, IJ_PHI0p7, IJ_PHI0p16,
     *     IJ_PHI0p07, IJ_PHI0p03, IJ_T850, IJ_T500, IJ_T300, IJ_Q850,
     *     IJ_Q500, IJ_Q300, IJ_PMCCLD, IJ_CLDTPPR, IJ_CLDCV, IJ_DSEV,
     *     IJ_CLDTPT, IJ_CLDCV1, IJ_CLDT1T,IJ_CLDT1P, IJ_T700, IJ_Q700,
     *     ij_wtrcld,ij_icecld,ij_optdw,ij_optdi,IJ_PBLHT, IJ_RH700,
     *     IJ_RH1, IJ_RH850, IJ_RH500, IJ_RH300, IJ_SWCRF, IJ_LWCRF,
     *     IJ_SWCRF2, IJ_LWCRF2, IJ_T100,IJ_Q100,IJ_RH100,
     *     IJ_TRNFP0, IJ_SRTR, IJ_NETH, IJ_SRNFP0, IJ_SRINCP0, IJ_SRNFG,
     *     IJ_SRINCG, IJ_TG1, IJ_RSIT, IJ_TDSL, IJ_TDCOMP, IJ_DTDP,
     *     IJ_RUNE, IJ_TS1, IJ_RUNLI, IJ_WS, IJ_TS, IJ_US, IJ_VS,
     *     IJ_SLP, IJ_UJET, IJ_VJET, IJ_PCLDL, IJ_PCLDM, IJ_PCLDH,
     *     IJ_BTMPW, IJ_SRREF, IJ_SRVIS, IJ_TOC2, IJ_TAUS, IJ_TAUUS,
     *     IJ_TAUVS, IJ_GWTR, IJ_QS, IJ_STRNGTS, IJ_ARUNU, IJ_DTGDTS,
     *     IJ_PUQ, IJ_PVQ, IJ_TGO, IJ_MSI, IJ_TGO2, IJ_EVAPO, ij_RHs,
     *     IJ_EVAPI, IJ_EVAPLI,IJ_EVAPE, IJ_F0OC,IJ_F0OI,IJ_F0LI,IJ_F0E,
     *     IJ_F1LI, IJ_SNWF, IJ_TSLI, IJ_SHDTLI, IJ_EVHDT,
     *     IJ_TRHDT, IJ_TMAXE, IJ_TMAXC, IJ_TMINC, IJ_TMNMX, IJ_PEVAP,
     *     IJ_WMSUM, IJ_PSCLD, IJ_PDCLD, IJ_DCNVFRQ, IJ_SCNVFRQ,
     *     IJ_EMTMOM, IJ_SMTMOM, IJ_FMU, IJ_FMV, IJ_SSTABX,
     *     IJ_FGZU, IJ_FGZV, IJ_ERVR, IJ_MRVR, IJ_SSS, IJ_PRECMC,
     *     IJ_LKON, IJ_LKOFF, IJ_LKICE, IJ_PTROP, IJ_TTROP, IJ_TSI,
     *     IJ_SSI1,IJ_SSI2,IJ_SMFX,     ! IJ_MSU2,IJ_MSU2R,
     *     IJ_MLTP,IJ_FRMP, IJ_P850, IJ_CLR_SRINCG,
     *     IJ_GPP, IJ_IPP, IJ_RAUTO, IJ_CLAB, IJ_DLEAF, IJ_LAI, !VEG DIAGNOSTICS
     *     IJ_SOILRESP, IJ_SOILCPOOLSUM, !additional veg diags (soil bgc)
     *     IJ_GICE, IJ_GWTR1, IJ_ZSNOW, IJ_AFLMLT, IJ_AERUNS, IJ_AERUNU,
     *     IJ_HTSOIL, IJ_HTSNOW, IJ_AINTRCP, IJ_MCCLDTP, IJ_MCCLDBS,
     *     IJ_SRNTP,IJ_TRNTP,IJ_CLR_SRNTP,IJ_CLR_TRNTP, IJ_TRSDN,
     *     IJ_TRSUP, IJ_CLR_SRNFG,IJ_CLR_TRDNG,IJ_CLR_SRUPTOA,
     *     IJ_CLR_TRUPTOA, IJ_CLDW, IJ_CLDI, IJ_QM, IJ_SSH, IJ_FWOC,
     *     IJ_FWIO, IJ_HTIO, IJ_STIO, IJ_DSKIN, IJ_MCCVTP, IJ_MCCVBS,
     *     IJ_SWDCLS,IJ_SWNCLS,IJ_LWDCLS,IJ_SWNCLT,IJ_LWNCLT,
     *     IJ_P1000,IJ_P925,IJ_P700,IJ_P600,IJ_P500, IJ_LI, IJ_LK,
     &     IJ_FVEG,IJ_GUSTI, IJ_MCCON, IJ_SRVDIR, IJ_SRVISSURF
     *     ,IJ_WISUM, IJ_SLPQ, IJ_PRESQ, IJ_RNFP1
     *     ,ij_dzwm,ij_dzim,ij_dzws,ij_dzis
     *     ,ij_3dnwm,ij_3dnim,ij_3dnws,ij_3dnis
     *     ,ij_3drwm,ij_3drim,ij_3drws,ij_3dris
     *     ,ij_3dlwm,ij_3dlim,ij_3dlws,ij_3dlis
     *     ,ij_ssprec,ij_mcprec,IJ_WMCLWP,IJ_WMCTWP
     &     ,ij_wdry,ij_wtke,ij_wmoist,ij_wsgcm,ij_wspdf
     &     ,ij_flam,ij_CtoG,ij_flash,ij_chl,ij_swaerrf,ij_lwaerrf
     *     ,ij_swaersrf,ij_lwaersrf,ij_swaerabs,ij_lwaerabs,ij_swaerrfnt
     *     ,ij_lwaerrfnt,ij_swaersrfnt,ij_lwaersrfnt,ij_swaerabsnt
     *     ,ij_lwaerabsnt,ij_evapsn,ij_irrW,ij_irrE,ij_irrW_tot
     *     ,ij_mwl,ij_gml,ij_mwlir,ij_gmlir,ij_irrgw,ij_irrgwE
     *     ,ij_kw, ij_alpha, ij_gasx, ij_rvrflo
     *     ,ij_sisnd,ij_tsice,ij_sisnwf,ij_sigrfr
     *     ,ij_sigrcg,ij_sigrlt,ij_sntosi,ij_sitopmlt,ij_sibotmlt
     *     ,ij_sihc,ij_siswd,ij_siswu,ij_silwd,ij_silwu,ij_sish
     *     ,ij_impmli,ij_imphli,ij_eicb,ij_micb, IJ_ERVRO, IJ_MRVRO
     *     ,IJ_IMPMGR,IJ_IMPHGR,IJ_IMPMKI,IJ_IMPHKI
     *     ,IJ_MLKtoGR,IJ_HLKtoGR, IJ_MSNFLOOD,IJ_HSNFLOOD
     *     ,ij_precli,ij_precsi,ij_precoo,ij_precgr
#ifdef HEALY_LM_DIAGS
     &     ,IJ_CROPS
#endif
C**** RDF: new aij indices for sunlight,subgrid ISCCP counts,precipitable water
#ifdef TES_SIM
     &     ,IJ_TES_SUNL
     &     ,IJ_TES_TTL,IJ_TES_TTL_R
     &     ,IJ_TES_GOOD,IJ_TES_GOOD_R
     &     ,IJ_TES_TTL_CLDY
     &         ,IJ_TES_TTL_CLDY_GOOD,IJ_TES_CLDY_GOOD_R
     &     ,IJ_TES_TGRND,IJ_TES_TGRND_GOOD,IJ_TES_TGRND_GOOD_R
     &     ,IJ_TES_PWL,IJ_TES_PWL_GOOD,IJ_TES_PWL_GOOD_R
     &     ,IJ_TES_PWH,IJ_TES_PWH_GOOD,IJ_TES_PWH_GOOD_R
     &     ,IJ_TES_CTP,IJ_TES_CTP_GOOD,IJ_TES_CTP_GOOD_R
     &     ,IJ_TES_COD,IJ_TES_COD_GOOD,IJ_TES_COD_GOOD_R
     &     ,IJ_TES_CODCLDY,IJ_TES_CODCLDY_GOOD,IJ_TES_CODCLDY_GOOD_R

#endif
!@var IJ_Gxx names for old AIJG arrays
      INTEGER, public ::
     *     ij_gbsw, ij_gbsbet, ij_gbetpen, ij_gvsw, ij_gbvswt,
     *     ij_gconatm, ij_gconcan, ij_gevppen, ij_gbst, ij_gbsevp,
     *     ij_gdcevp, ij_gwcevp, ij_gvst, ij_gwtbl, ij_gvswet,
     *     ij_gbetat, ij_gbssnd, ij_gvssnd
!@var IJ_GWx names for gravity wave diagnostics
      INTEGER, public ::
     &     IJ_GW1,IJ_GW2,IJ_GW3,IJ_GW4,IJ_GW5,IJ_GW6,IJ_GW7,IJ_GW8
     *     ,IJ_GW9
!@var IJ_[MHS][UV]SI indices for sea ice mass/heat/salt transport diags
      INTEGER, public ::
     &     IJ_MUSI,IJ_MVSI,IJ_HUSI,IJ_HVSI,IJ_SUSI,IJ_SVSI
!@var IJ_xxxI names for ISCCP diagnostics
      INTEGER, public ::
     &     IJ_CTPI,IJ_TAUI,IJ_LCLDI,IJ_MCLDI,IJ_HCLDI,IJ_TCLDI,IJ_SCLDI
#ifdef ACCMIP_LIKE_DIAGS
!@var IJ_fcghg GHG forcing diagnostics (2=LW,SW, 4=CH4,N2O,CFC11,CFC12)
      INTEGER, public, dimension(2,4) :: IJ_fcghg
#endif
c weighting fractions
      INTEGER, public :: IJ_PSOIL,IJ_CLRSKY,IJ_POCEAN,IJ_POPOCN,IJ_VSFR
     *     ,IJ_BSFR,IJ_POPWAT,IJ_PWATER
c derived/composite diagnostics
      INTEGER, public ::
     *  ij_topo, ij_jet, ij_wsmn, ij_jetdir, ij_wsdir, ij_grow,
     *  ij_netrdp, ij_albp, ij_albg, ij_albv, ij_ntdsese, ij_ntdsete,
     *  ij_fland, ij_dzt1, ij_albgv, ij_msu2,ij_msu3,ij_msu4,
     *  ij_Tatm, ij_RTSE, ij_HWV, ij_PVS


!@param LEGEND "contour levels" for ij-maps
      CHARACTER(LEN=40), DIMENSION(25), PARAMETER, public :: LEGEND=(/ !
     1  '0=0,1=5...9=45,A=50...K=100             ', ! ir_pct    fac=.2
     2  '0=0...9=90,A=100...I=180...R=270        ', ! ir_angl       .1
     3  '1=.5...9=4.5,A=5...Z=17.5,+=MORE        ', ! ir_0_18        2
     4  '1=.1...9=.9,A=1...Z=3.5,+=MORE          ', ! ir_0_4        10
     5  '1=2...9=18,A=20...Z=70,+=MORE           ', ! ir_0_71       .5
     6  '1=50...9=450,A=500...Z=1750,+=MORE      ', ! ir_0_1775     .02
     7  '1=100...9=900,A=1000...Z=3500,+=MORE    ', ! ir_0_3550     .01
     8  '1=20...9=180,A=200...Z=700,+=MORE       ', ! ir_0_710      .05
     9  'A=1...Z=26,3=30...9=90,+=100-150,*=MORE ', ! ir_0_26_150    1
     O  '0=0,A=.1...Z=2.6,3=3...9=9,+=10-15      ', ! ir_0_3_15     10
     1  '-=LESS,Z=-78...0=0...9=27,+=MORE        ', ! ir_m80_28     .33
     2  '-=LESS,Z=-260...0=0...9=90,+=MORE       ', ! ir_m265_95    .1
     3  '-=LESS,Z=-520...0=0...9=180,+=MORE      ', ! ir_m530_190   .05
     4  '-=LESS,Z=-1300...0=0...9=450,+=MORE     ', ! ir_m1325_475  .02
     5  '-=LESS,Z=-2600...0=0...9=900,+=MORE     ', ! ir_m2650_950  .01
     6  '-=LESS,Z=-3900...0=0...9=1350,+=MORE    ', ! ir_m3975_1425 .007
     7  '-=LESS,Z=-5200...0=0...9=1800,+=MORE    ', ! ir_m5300_1900 .005
     8  '-=LESS,9=-.9...0=0,A=.1...Z=2.6,+=MORE  ', ! ir_m1_3       10
     9  '-=LESS,9=-45...0=0,A=5...I=45...+=MORE  ', ! ir_m45_130    .2
     O  '-=LESS,9=-90...0=0,A=10...Z=260,+=MORE  ', ! ir_m95_265    .1
     1  '-=LESS,9=-180...A=20...Z=520,+=MORE     ', ! ir_m190_530   .05
     2  '-=LESS,9=-9...0=0,A=1...Z=26,+=MORE     ', ! ir_m9_26       1
     3  '-=LESS,9=-36...0=0,A=4...Z=104,+=MORE   ', ! ir_m38_106    .25
     4  '1=5...9=45,A=50...Z=175,+=MORE          ', ! ir_0_180      .2
     5  '9=-512...1=-2,0=0,A=2,B=4,C=8...+=MORE  '/)! ir_log2       1.
!@var ir_xxxx names for indices to LEGEND indicating the (rounded) range
      integer, parameter, public ::
     &     ir_pct=1, ir_angl=2, ir_0_18=3, ir_0_4=4,
     * ir_0_71=5, ir_0_1775=6, ir_0_3550=7, ir_0_710=8, ir_0_26_150=9,
     * ir_0_3_15=10, ir_m80_28=11, ir_m265_95=12, ir_m530_190=13,
     * ir_m1325_475=14, ir_m2650_950=15, ir_m3975_1425=16,
     * ir_m5300_1900=17, ir_m1_3=18, ir_m45_130=19, ir_m95_265=20,
     * ir_m190_530=21, ir_m9_26=22, ir_m38_106=23, ir_0_180=24,
     * ir_log2=25
!@var fac_legnd = 1/(range_of_1_colorbox)
      real*8, dimension(25), public :: fac_legnd=(/
     1      1d0/5,  1d0/10,    2.d0,   10.d0,   1d0/2,
     6     1d0/50, 1d0/100,  1d0/20,    1.d0,   10.d0,
     1      1d0/3,  1d0/10,  1d0/20,  1d0/50, 1d0/100,
     6    1d0/150, 1d0/200,   10.d0,   1d0/5,  1d0/10,
     1     1d0/20,    1.d0,   1d0/4,   1d0/5,     1d0  /)

!@param CBAR "color bars" for ij-maps
      CHARACTER(LEN=38), PARAMETER, DIMENSION(5), public :: CBAR=(/
     &     ' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+',  ! ib_pos
     &     ' 0123456789ABCDEFGHIJKX               ',  ! ib_pct
     &     '-9876543210ABCDEFGHIJKLMNOPQRSTUVWXYZ+',  ! ib_npp,ib_ntr
     &     ' 0ABCDEFGHIJKLMNOPQRSTUVWXYZ3456789+* ',  ! ib_hyb
     &     '-ZYXWVUTSRQPONMLKJIHGFEDCBA0123456789+'/) ! ib_nnp
!@var ib_xxx indices for color bars
      integer, parameter, public ::
     &     ib_pos=1,ib_pct=2,ib_npp=3,ib_hyb=4,ib_nnp=5
     &     ,ib_ntr=6

!@dbparam isccp_diags: if 1 accumulate ISCCP cloud data (default 0)
      INTEGER, public :: isccp_diags = 0

!@dbparam lh_diags: if 1 accumulate 3D latent heating profiles (default 0)
      INTEGER, public :: lh_diags = 0

!@var SCALE_IJ scaling for weighted AIJ diagnostics
      REAL*8, DIMENSION(KAIJ), public :: SCALE_IJ
!@var NAME_IJ,UNITS_IJ Names/Units of lat/lon IJ diagnostics
      character(len=sname_strlen), dimension(kaij), public :: name_ij
      character(len=units_strlen), dimension(kaij), public :: units_ij
!@var LNAME_IJ Long names of lat/lon IJ diagnostics
      character(len=lname_strlen), dimension(kaij), public :: lname_ij
!@var HEMIS_IJ hemispheric/global averages of AIJ
      real*8, dimension(:,:,:), allocatable, public :: hemis_ij
!@var nwts_ij = number of weight-ij-arrays used in IJ-diagnostics
      integer, parameter, public :: nwts_ij = 8
!@var wt_ij various weight-arrays use in ij-diagnostics
      real*8, dimension(IM,JM,NWTS_IJ), public :: wt_ij
!@var IW_xxx index for weight-array
      integer, parameter, public :: iw_all=1 , iw_ocn=2 , iw_lake=3,
     *   iw_lice=4 , iw_soil=5 , iw_bare=6 , iw_veg=7, iw_land=8
!@var IR_IJ range indices for IJ diagnostics
      integer, dimension(kaij), public :: ir_ij
!@var IA_IJ IDACC indexes for lat/lon IJ diagnostics
      integer, dimension(kaij), public :: ia_ij
!@var [ij]grid_ij 1=primary grid  2=secondary grid
      integer, dimension(kaij), public :: igrid_ij,jgrid_ij
!@var denom_ij index of AIJ element to use as time/area weight
      integer, dimension(kaij), public :: denom_ij

!@var JL_xxx, JK_xxx names for AJL indices
!@+   JL/JK refer to model versus constant-pressure levels
      INTEGER, public ::
     &     jl_mcmflx,jl_srhr,jl_trcr,jl_sshr,jl_trbhr,jl_mchr
     &     ,jl_dtdyn,jl_totcld,jl_mcdflx,jl_sscld,jl_mccld
     &     ,jl_rhe,jl_damdc,jl_dammc,jl_mchphas,jl_mcdtotw
     &     ,jl_mcldht,jl_trbke,jl_trbdlht,jl_mcheat,jl_mcdry
     &     ,jl_mcdeep,jl_mcshlw,jl_cldmc,jl_cldss,jl_csizmc,jl_csizss
     &     ,jl_wcld,jl_icld,jl_wcod,jl_icod,jl_wcsiz,jl_icsiz
     &     ,jl_cnumwm,jl_cnumim,jl_cnumws,jl_cnumis
     &     ,jl_dpa,jl_dpasrc,jl_dwasrc,jl_wcldwt,jl_icldwt
     &     ,jl_rad_cool
     &     ,jl_epacwt,jl_wpacwt
     &     ,jl_uepac,jl_vepac,jl_wepac,jl_uwpac,jl_vwpac,jl_wwpac
     &     ,jl_dudtsdrg,jl_dtdtsdrg
      INTEGER, public ::
     &      jl_dudfmdrg,jl_dumtndrg,jl_dushrdrg,jl_dumcdrgm10
     &     ,jl_dumcdrgp10,jl_dumcdrgm40,jl_dumcdrgp40,jl_dumcdrgm20
     &     ,jl_dumcdrgp20,jl_sdifcoef,jl_dudtsdif,jl_dudtvdif
     &     ,jl_gwfirst,jl_mcdrgpm10,jl_mcdrgpm40,jl_mcdrgpm20,jl_sumdrg

      INTEGER, public ::
     &     JK_hght, JK_dpwt, JK_tx, JK_q, JK_cldh2o ,JK_rh
     &     ,JK_cldwtr, JK_cldice


!@var JGRID_U, JGRID_KE latitudes at which U-wind and KE diags are defined
!@+   (1 for primary latitudes, 2 for secondary latitudes)
!@+   In future, can be other than the values in GEOM if requested.
      integer, public :: jgrid_u, jgrid_ke
      integer, parameter, public ::
     &     igridc=0,igride=1,jgridc=0,jgride=2,kgridc=0,kgride=4,
     &     ijkgridc=igridc+jgridc+kgridc

!@var SNAME_JL Names of lat-sigma JL diagnostics
      character(len=sname_strlen), dimension(kajl), public :: sname_jl
!@var LNAME_JL,UNITS_JL Descriptions/Units of JL diagnostics
      character(len=lname_strlen), dimension(kajl), public :: lname_jl
      character(len=units_strlen), dimension(kajl), public :: units_jl
!@var SCALE_JL printout scaling factors for JL diagnostics
      REAL*8, dimension(kajl), public :: scale_jl
!@var IA_JL,JGRID_JL,LGRID_JL idacc-numbers,gridtypes for JL diagnostics
      integer, dimension(kajl), public :: ia_jl,jgrid_jl,lgrid_jl
!@var POW_JL printed output scaled by 10**(-pow_jl)
      integer, dimension(kajl), public :: pow_jl
!@var DENOM_JL index of AJL element to use as weight
      integer, dimension(kajl), public :: denom_jl
!@var HEMIS_JL hemispheric/global averages of AJL
!@var VMEAN_JL vertical sums of AJL
      real*8, dimension(:,:,:), allocatable, public :: hemis_jl,vmean_jl
!@param [CTR,EDG]_[ML,CP] tags for center,edge model-layer,const-pres
!@+   vertical grids
      integer, parameter, public :: ctr_ml=1,edg_ml=2,ctr_cp=3,edg_cp=4

!@var NAME_SJL Names of radiative-layer-only SJL diagnostics
      character(len=sname_strlen), dimension(kasjl), public :: name_sjl
!@var LNAME_SJL,UNITS_SJL Descriptions/Units of SJL diagnostics
      character(len=lname_strlen), dimension(kasjl), public :: lname_sjl
      character(len=units_strlen), dimension(kasjl), public :: units_sjl
!@var SCALE_SJL printout scaling factors for SJL diagnostics
      REAL*8, dimension(kasjl), public :: scale_sjl
!@var IA_SJL idacc-numbers for SJL diagnostics
      integer, dimension(kasjl), public :: ia_sjl

!@var SNAME_GC Names of lat-pressure GC diagnostics
      character(len=sname_strlen), dimension(kagcx), public :: sname_gc
!@var LNAME_GC,UNITS_GC Descriptions/Units of GC diagnostics
      character(len=lname_strlen), dimension(kagcx), public :: lname_gc
      character(len=units_strlen), dimension(kagcx), public :: units_gc
!@var SCALE_GC printout scaling factors for GC diagnostics
      REAL*8, dimension(kagcx), public :: scale_gc
!@var IA_GC,JGRID_GC,LGRID_GC idacc-numbers,gridtypes for GC diagnostics
      integer, dimension(kagcx), public :: ia_gc,jgrid_gc,lgrid_gc
!@var DENOM_GC index of AGC element to use as weight
      integer, dimension(kagcx), public :: denom_gc
!@var POW_GC printed output scaled by 10**(-pow_gc)
      integer, dimension(kagcx), public :: pow_gc
!@var HEMIS_GC hemispheric/global averages of AGC
!@var VMEAN_GC vertical sums of AGC
      real*8, dimension(:,:,:), allocatable, public :: hemis_gc,vmean_gc
!@var lat_gc latitudes of the primary grid for GC diagnostics
      real*8, dimension(jmlat), public :: lat_gc,lat_gc2

!@var SCALE_IJK scaling for weighted AIJK diagnostics
      REAL*8, DIMENSION(Kaijk), public :: SCALE_IJK
!@var OFF_IJK offset for weighted AIJK diagnostics
      REAL*8, DIMENSION(Kaijk), public :: OFF_IJK
!@var NAME_IJK Names of lon-lat-pressure IJK diagnostics
      character(len=sname_strlen), dimension(kaijk), public :: name_ijk
!@var LNAME_IJK,UNITS_IJK Descriptions/Units of IJK diagnostics
      character(len=lname_strlen), dimension(kaijk), public ::
     &     lname_ijk
      character(len=units_strlen), dimension(kaijk), public ::
     &     units_ijk
!@var jgrid_ijk 1=primary grid  2=secondary grid
      integer, dimension(Kaijk), public :: jgrid_ijk,denom_ijk,ia_ijk

!@var SCALE_IJL scale factor for AIJL diagnostics
      REAL*8, DIMENSION(KAIJL), public :: SCALE_IJL
!@var IA_IJL,DENOM_IJL  idacc-numbers,weights for AIJL diagnostics
      INTEGER, DIMENSION(KAIJL), public :: IA_IJL,DENOM_IJL,LGRID_IJL
     &                                                     ,JGRID_IJL
!@var NAME_IJL Names of lon-lat-level IJL diagnostics
      character(len=sname_strlen), dimension(kaijl), public :: name_ijl
!@var LNAME_IJL,UNITS_IJL Descriptions/Units of IJL diagnostics
      character(len=lname_strlen), dimension(kaijl), public ::
     &     lname_ijl
      character(len=units_strlen), dimension(kaijl), public ::
     &     units_ijl

      character(len=sname_strlen), dimension(kwp), public :: name_wave
      character(len=units_strlen), dimension(kwp), public :: units_wave
      character(len=lname_strlen), dimension(kwp), public :: lname_wave

      character(len=sname_strlen), dimension(kcon), public ::
     &     name_consrv = 'unused'
      character(len=units_strlen), dimension(kcon), public ::
     &     units_consrv
      character(len=lname_strlen), dimension(kcon), public ::
     &     lname_consrv
!@var HEMIS_CONSRV hemispheric/global averages of CONSRV
      real*8, dimension(:,:), allocatable, public :: hemis_consrv

!@var J50N,J70N,J5NUV,J5SUV,J5S,J5N special latitudes for AIL diags
      INTEGER, PARAMETER, public :: J50N  = (50.+90.)*(JM-1)/180.+1.5
      INTEGER, PARAMETER, public :: J70N  = (70.+90.)*(JM-1)/180.+1.5
      INTEGER, PARAMETER, public :: J5NUV = (90.+5.)*(JM-1.)/180.+2.
      INTEGER, PARAMETER, public :: J5SUV = (90.-5.)*(JM-1.)/180.+2.
      INTEGER, PARAMETER, public :: J5N   = (90.+5.)*(JM-1.)/180.+1.5
      INTEGER, PARAMETER, public :: J5S   = (90.-5.)*(JM-1.)/180.+1.5

      character(len=sname_strlen), dimension(ndiuvar), public :: name_dd
      character(len=units_strlen), dimension(ndiuvar), public ::
     &     units_dd
      character(len=lname_strlen), dimension(ndiuvar), public ::
     &     lname_dd
      real*8, dimension(ndiuvar), public :: scale_dd
      integer, dimension(ndiuvar), public :: denom_dd

!@var IDD_xxx names for diurnal diagnostics
      INTEGER, public ::
c     standard set of names
     &     IDD_ISW, IDD_PALB, IDD_GALB, IDD_ABSA, IDD_ECND,
     *     IDD_SPR, IDD_PT5, IDD_TS, IDD_TG1, IDD_Q5, IDD_QS,
     *     IDD_QG, IDD_SWG, IDD_LWG, IDD_SH, IDD_LH, IDD_HZ0, IDD_UG,
     *     IDD_VG, IDD_WG, IDD_US, IDD_VS, IDD_WS, IDD_CIA, IDD_RIS,
     *     IDD_RIG, IDD_CM, IDD_CH, IDD_CQ, IDD_EDS, IDD_DBL, IDD_DCF,
     *     IDD_LDC, IDD_PR, IDD_EV, IDD_DMC, IDD_SMC, IDD_CL7, IDD_W,
     *     IDD_CCV, IDD_SSP, IDD_MCP ! 56
c     names for one layer dust diagnostics
     &     ,idd_wtke,idd_wd,idd_wm,idd_wsgcm,idd_wspdf,idd_wtrsh
     &     ,idd_emis,idd_emis2,idd_ws2,idd_ustar,idd_us3,idd_stress
     &     ,idd_lmon,idd_rifl,idd_wet,idd_grav,idd_turb ! +17
c     names for llmax_dd2 layers dust diagnostics
     &     ,idd_u1,idd_v1,idd_uv1,idd_t1,idd_qq1,idd_p1,idd_w1,idd_phi1
     &     ,idd_sr1,idd_tr1,idd_load1,idd_conc1,idd_tau1,idd_tau_cs1 ! +14*ls1
c     names for npbl layers dust diagnostics
     &     ,idd_zpbl1,idd_uabl1,idd_vabl1,idd_uvabl1,idd_tabl1
     &     ,idd_qabl1           ! +6*npbl
c     names for npbl-1 layers dust diagnostics
     &     ,idd_zhat1,idd_e1,idd_km1,idd_ri1 ! +4*(npbl-1)
c    hourly AMP diagnostics
     *     ,idd_diam, idd_aot, idd_lwp, idd_ccn, idd_cdnc
     *     ,idd_mass, idd_numb, idd_so2, idd_lwc, idd_ncL

!@var tf_xxx tsfrez diagnostic names
      INTEGER, public :: tf_day1,tf_last,tf_lkon,tf_lkoff
      character(len=sname_strlen), dimension(ktsf), public :: name_tsf
      character(len=units_strlen), dimension(ktsf), public :: units_tsf
      character(len=lname_strlen), dimension(ktsf), public :: lname_tsf

      character(len=8), dimension(ntype), public :: stype_names=
     &     (/ 'OCEAN   ','OCEANICE','EARTH   ',
     &        'LANDICE ','LAKE    ','LAKEICE ' /)

!@param NTYPE_OUT number of output budgets pages
      INTEGER, PARAMETER :: NTYPE_OUT=NTYPE+3  ! to include comp/regio
C**** Expanded version of surfaces (including composites)
!@var TERRAIN name of surface type
      CHARACTER*16, DIMENSION(NTYPE_OUT), PARAMETER :: TERRAIN = (/
     *     '    (GLOBAL)','(OPEN OCEAN)',' (OCEAN ICE)','     (OCEAN)',
     *     '      (LAND)','  (LAND ICE)',' (OPEN LAKE)','  (LAKE ICE)',
     *     '     (LAKES)'/)
C**** weighting functions for surface types
      REAL*8, DIMENSION(NTYPE_OUT,NTYPE), PARAMETER ::
     *     WTJ_COMP=RESHAPE(          ! separate types + composites
     *     (/1.,1.,0.,1.,0.,0.,0.,0.,0., 1.,0.,1.,1.,0.,0.,0.,0.,0.,
     *       1.,0.,0.,0.,1.,0.,0.,0.,0., 1.,0.,0.,0.,0.,1.,0.,0.,0.,
     *       1.,0.,0.,0.,0.,0.,1.,0.,1., 1.,0.,0.,0.,0.,0.,0.,1.,1./),
     *     (/NTYPE_OUT,NTYPE/) )
      public :: ntype_out,terrain,wtj_comp

c idacc-indices of various processes
      integer, parameter, public ::
     &     ia_src=1, ia_rad=2, ia_srf=3, ia_dga=4, ia_d4a=5, ia_d5f=6,
     *     ia_d5d=7, ia_d5s=8, ia_12hr=9, ia_filt=10, ia_rad_frc=11,
     *     ia_inst=12

!@var PLE,PLM, PLE_DN ref pressures at upper, middle and lower edge
      REAL*8, DIMENSION(LM), public :: PLE
      REAL*8, DIMENSION(LM), public :: PLE_DN
      REAL*8, DIMENSION(LM+LM_REQ), public :: PLM
!@var P1000K scaling to change reference pressure from 1mb to 1000mb
      REAL*8, public :: P1000K
CXXXX inci,incj NOT GRID-INDPENDENT
!@var inci,incj print increments for i and j, so maps/tables fit on page
      integer, parameter, public ::
     &     inci=(im+35)/36,incj=(JM+23)/24, jmby2=jm/2
!@var linect = current line on page of print out
      integer, public :: linect

!@var LMOMAX max no. of layers in any ocean
      INTEGER, PARAMETER, public :: LMOMAX=50
!@var ZOC, ZOC1 ocean depths for diagnostics (m) (ONLY FOR DEEP OCEAN)
      REAL*8, public :: ZOC(LMOMAX) = 0. , ZOC1(LMOMAX+1) = 0.

!@param L_ROSSBY_NUMBER length scale for budget-page Rossby number
      real*8, parameter, public :: l_rossby_number=1d6 ! 1000 km


#ifdef NEW_IO
      type(cdl_type), public :: cdl_latbudg,cdl_heights

!@var CDL_J consolidated metadata for AJ output fields in CDL notation
!@+   CDL_REG                         AREG
      type(cdl_type), public :: cdl_j,cdl_reg
!@var CDL_IJ consolidated metadata for AIJ output fields in CDL notation
      type(cdl_type), public ::
     &     cdl_ij_template,cdl_ij_latlon_template,
     &     cdl_ij,cdl_ij_latlon
!@var CDL_JL consolidated metadata for AJL output fields in CDL notation
      type(cdl_type), public :: cdl_jl,cdl_jl_template
!@var CDL_GC consolidated metadata for AGC output fields in CDL notation
      type(cdl_type), public :: cdl_gc
!@var CDL_IJL consolidated metadata for AIJL output fields in CDL notation
      type(cdl_type), public ::
     &     cdl_ijl_template,cdl_ijl_latlon_template,
     &     cdl_ijl,cdl_ijl_latlon
!@var CDL_IJK consolidated metadata for AIJK output fields in CDL notation
!@+   (latlon-only for the moment)
      type(cdl_type), public :: cdl_ijk
!@var CDL_CONSRV consolidated metadata for CONSRV output fields in
!@+   CDL notation
      type(cdl_type), public :: cdl_consrv
!@var CDL_DD consolidated metadata for ADIURN output fields in CDL notation
      type(cdl_type), public :: cdl_dd,cdl_hd


c declarations that facilitate switching between restart and acc
c instances of arrays
      target :: aj,aj_out,areg,areg_out,agc,agc_out
      REAL*8, dimension(:,:,:), public, pointer ::
     &     AJ_ioptr
      REAL*8, dimension(:,:), public, pointer ::
     &     AREG_ioptr
      REAL*8, dimension(:,:,:), public, pointer ::
     &     AGC_ioptr
#endif

      END MODULE DIAG_COM

      SUBROUTINE ALLOC_DIAG_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID,GET,AM_I_ROOT
      USE RESOLUTION, ONLY : IM,LM
      USE MODEL_COM, ONLY : NTYPE,lm_req
      USE DIAG_COM, ONLY : KAJ,KCON,KAJL,KASJL,KAIJ,KAGC,KAIJK,
     &                   KGZ,KOA,KTSF,nwts_ij,KTD,NREG,KAIJL,JM_BUDG
      USE DIAG_COM, ONLY : SQRTM,AJ_loc,JREG,AJL_loc,ASJL_loc
     *     ,AIJ_loc,AGC_loc,AIJK_loc,AIJL_loc,AFLX_ST
     *     ,Z_inst,RH_inst,T_inst,TDIURN,TSFREZ_loc,OA,P_acc,PM_acc
     *     ,saveHCLDI,saveMCLDI,saveLCLDI,saveCTPI,saveTAUI,saveSCLDI
     *     ,saveTCLDI,saveMCCLDTP,R_acc
      USE DIAG_COM, ONLY : JMLAT,AJ,AJL,ASJL,AGC,AJ_OUT,ntype_out
     &     ,AGC_out
      USE DIAG_COM, ONLY : hemis_j,hemis_jl,vmean_jl,hemis_consrv
     &     ,hemis_gc,vmean_gc,hemis_ij
#ifdef TRACERS_SPECIAL_Shindell
      USE DIAG_COM, ONLY : o_inst,n_inst,m_inst,x_inst
#endif
#ifdef TES_LIKE_DIAGS
      USE DIAG_COM, ONLY : T_more,Q_more,KGZmore
#ifdef TRACERS_SPECIAL_Shindell
     &     ,o_more,n_more,m_more,x_more
#endif
#endif
C**** RDF 
#ifdef TES_SIM
      USE DIAG_COM, ONLY : saveBOXTAU,saveBOXCTP
      USE CLOUDS_COM, only: ncol
#endif
      use diag_zonal, only : get_alloc_bounds

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid
      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER
      LOGICAL, SAVE :: init = .false.
      integer :: j_0budg,j_1budg,j_0jk,j_1jk

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

      CALL GET( grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H  )
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      call get_alloc_bounds(grid,
     &     j_strt_budg=j_0budg,j_stop_budg=j_1budg,
     &     j_strt_jk=j_0jk,j_stop_jk=j_1jk)

      ALLOCATE(  JREG(I_0H:I_1H, J_0H:J_1H),
     &         SQRTM(I_0H:I_1H, J_0H:J_1H),
     &         STAT = IER)

      ALLOCATE(
     &         AJ_loc(J_0BUDG:J_1BUDG, KAJ, NTYPE),
     &         AJL_loc(J_0BUDG:J_1BUDG, LM, KAJL),
     &         ASJL_loc(J_0BUDG:J_1BUDG,LM_REQ,KASJL),
     &         AGC_loc(J_0JK:J_1JK,LM,KAGC),
     &         AIJ_loc(I_0H:I_1H,J_0H:J_1H,KAIJ),
     &         Z_inst(KGZ,I_0H:I_1H,J_0H:J_1H),
     &         RH_inst(KGZ,I_0H:I_1H,J_0H:J_1H),
     &         T_inst(KGZ,I_0H:I_1H,J_0H:J_1H),
     &         TSFREZ_loc(I_0H:I_1H,J_0H:J_1H,KTSF),
     &         P_acc(I_0H:I_1H,J_0H:J_1H),
     &         PM_acc(I_0H:I_1H,J_0H:J_1H),
     &         R_acc(I_0H:I_1H,J_0H:J_1H),
     &         TDIURN(I_0H:I_1H,J_0H:J_1H,KTD),
     &         OA(I_0H:I_1H,J_0H:J_1H,KOA),
     &         saveHCLDI(I_0H:I_1H,J_0H:J_1H),
     &         saveMCLDI(I_0H:I_1H,J_0H:J_1H),
     &         saveLCLDI(I_0H:I_1H,J_0H:J_1H),
     &         saveCTPI(I_0H:I_1H,J_0H:J_1H),
     &         saveTAUI(I_0H:I_1H,J_0H:J_1H),
     &         saveSCLDI(I_0H:I_1H,J_0H:J_1H),
     &         saveTCLDI(I_0H:I_1H,J_0H:J_1H),
     &         saveMCCLDTP(I_0H:I_1H,J_0H:J_1H),
C**** RDF
#ifdef TES_SIM
     &         saveBOXTAU(I_0H:I_1H,J_0H:J_1H,ncol),
     &         saveBOXCTP(I_0H:I_1H,J_0H:J_1H,ncol),
#endif 
#ifdef TRACERS_SPECIAL_Shindell
     &         O_inst(KGZ,I_0H:I_1H,J_0H:J_1H),
     &         M_inst(KGZ,I_0H:I_1H,J_0H:J_1H),
     &         N_inst(KGZ,I_0H:I_1H,J_0H:J_1H),
     &         X_inst(KGZ,I_0H:I_1H,J_0H:J_1H),
#endif
#ifdef TES_LIKE_DIAGS
     &         Q_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
     &         T_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
#ifdef TRACERS_SPECIAL_Shindell
     &         O_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
     &         N_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
     &         M_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
     &         X_more(KGZmore,I_0H:I_1H,J_0H:J_1H),
#endif
#endif
     &         STAT = IER)

      P_acc=0.d0; PM_acc=0.d0
      R_acc=0.d0

      ALLOCATE( AIJK_loc(I_0H:I_1H,J_0H:J_1H,LM,KAIJK),
     &         AFLX_ST(LM+LM_REQ+1,I_0H:I_1H,J_0H:J_1H,5),
     &         STAT = IER)


      ALLOCATE( AIJL_loc(I_0H:I_1H,J_0H:J_1H,LM,KAIJL))

c allocate master copies of budget- and JK-arrays on root
      if(am_i_root()) then
        ALLOCATE(AJ(JM_BUDG, KAJ, NTYPE),
     &           AJL(JM_BUDG, LM, KAJL),
     &           ASJL(JM_BUDG,LM_REQ,KASJL),
     &           AGC(JMLAT,LM,KAGC),
     &           AGC_out(JMLAT,LM,KAGC),
     &           STAT = IER)
        allocate(aj_out(jm_budg,kaj,ntype_out))
        allocate(hemis_j(3,kaj,ntype_out))
        allocate(hemis_jl(3,lm,kajl))
        allocate(vmean_jl(jm_budg+3,1,kajl))
        allocate(hemis_consrv(3,kcon))
        allocate(hemis_gc(3,lm,kagc))
        allocate(vmean_gc(jmlat+3,1,kagc))
        allocate(hemis_ij(1,3,kaij))
      else
        ALLOCATE(AJ(1,1,1),
     &           AJL(1,1,1),
     &           ASJL(1,1,1),
     &           AGC(1,1,1),
     &           AGC_out(1,1,1),
     &        STAT = IER)
        allocate(aj_out(1,1,1))
        allocate(hemis_j(1,1,1))
        allocate(hemis_jl(1,1,1))
        allocate(vmean_jl(1,1,1))
        allocate(hemis_consrv(1,1))
        allocate(hemis_gc(1,1,1))
        allocate(vmean_gc(1,1,1))
        allocate(hemis_ij(1,1,1))

      endif

      RETURN
      END SUBROUTINE ALLOC_DIAG_COM

      SUBROUTINE ALLOC_ijdiag_glob
!@sum  To allocate large global arrays only when needed
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE RESOLUTION, ONLY : IM,JM,LM
      USE DOMAIN_DECOMP_ATM, Only : AM_I_ROOT
      USE DIAG_COM, ONLY : KAIJ,KAIJK,KOA,KTSF,KTD,KAIJL
      USE DIAG_COM, ONLY : AIJ,AIJK,AIJL,TSFREZ,TDIURN_GLOB,OA_GLOB
      IMPLICIT NONE
      INTEGER :: IER

      if(AM_I_ROOT()) then
         ALLOCATE(AIJ(IM,JM,KAIJ),
     &        TSFREZ(IM,JM,KTSF),
     &        AIJK(IM,JM,LM,KAIJK),
     &        AIJL(IM,JM,LM,KAIJL),
     &        TDIURN_glob(IM, JM, KTD),
     &        OA_glob(IM, JM, KOA))
      else
         ALLOCATE(AIJ(1,1,1),
     &        TSFREZ(1,1,1),
     &        AIJK(1,1,1,1),
     &        AIJL(1,1,1,1),
     &        TDIURN_glob(1,1,1),
     &        OA_glob(1,1,1))
      end if


      RETURN
      END SUBROUTINE ALLOC_ijdiag_glob

      SUBROUTINE DEALLOC_ijdiag_glob
!@sum  To deallocate large global arrays not currently needed
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE DOMAIN_DECOMP_ATM, Only : AM_I_ROOT
      USE DIAG_COM, ONLY : AIJ,AIJK,AIJL,TSFREZ,TDIURN_GLOB,OA_GLOB

      IMPLICIT NONE
      DEALLOCATE(AIJ,AIJK,AIJL,TSFREZ,TDIURN_glob,OA_glob)

      RETURN
      END SUBROUTINE DEALLOC_ijdiag_glob

      SUBROUTINE io_diags(kunit,it,iaction,ioerr)
!@sum  io_diag reads and writes diagnostics to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,ioread_single,irerun
     *    ,iowrite,iowrite_mon,iowrite_single,lhead, idacc,nsampl
     *    ,Kradia
      USE DIAG_COM
      USE DOMAIN_DECOMP_1D, Only : grid, GET, PACK_DATA, UNPACK_DATA
      USE DOMAIN_DECOMP_1D, Only : PACK_COLUMN, UNPACK_COLUMN
      USE DOMAIN_DECOMP_1D, Only : AM_I_ROOT
      USE DOMAIN_DECOMP_1D, Only : ESMF_BCAST
      IMPLICIT NONE

!@param KACC total number of diagnostic elements
      INTEGER, PARAMETER :: KACC= JM_BUDG*KAJ*NTYPE + NREG*KAJ
     *     + JM_BUDG*LM*KAJL + JM_BUDG*LM_REQ*KASJL + IM*JM*KAIJ +
     *     IM*JM*LM*KAIJL + NEHIST*HIST_DAYS + JM_BUDG*KCON +
     *     (IMLONH+1)*KSPECA*NSPHER + KTPE*NHEMI +
     *     HR_IN_DAY*NDIUVAR*NDIUPT +
     *     RE_AND_IM*Max12HR_sequ*NWAV_DAG*KWP + JM*LM*KAGC +
     *     IM*JM*LM*KAIJK+ntau*npres*nisccp
#ifndef NO_HDIURN
     *     + NDIUVAR*NDIUPT*HR_IN_MONTH
#endif
!@var AJ4,...,AFLX4 real*4 dummy arrays needed for postprocessing only
      ! REAL*4 AJ4(JM_BUDG,KAJ,NTYPE),AREG4(NREG,KAJ)
      ! REAL*4 AJL4(JM_BUDG,LM,KAJL),ASJL4(JM_BUDG,LM_REQ,KASJL),AIJ4(IM,JM,KAIJ)
      ! REAL*4 AIJL4(IM,JM,LM,KAIJL),ENERGY4(NEHIST,HIST_DAYS)
      ! REAL*4 CONSRV4(JM_BUDG,KCON),SPECA4(IMLONH+1,KSPECA,NSPHER)
      ! REAL*4 ATPE4(KTPE,NHEMI),ADIURN4(NDIUVAR,NDIUPT,HR_IN_DAY)
      ! REAL*4 WAVE4(RE_AND_IM,Max12HR_sequ,NWAV_DAG,KWP)
      ! REAL*4 AGC4(JM,LM,KAGC),AIJK4(IM,JM,LM,KAIJK)
      ! REAL*4 AISCCP4(ntau,npres,nisccp)
      ! REAL*4 TSFREZ4(IM,JM,KTSF),AFLX4(LM+LM_REQ+1,IM,JM,5)
      REAL*4,allocatable, dimension(:,:,:) :: AJ4,AJL4,ASJL4
     *            ,AIJ4,AGC4, SPECA4, ADIURN4, AISCCP4, TSFREZ4
      REAL*4,allocatable, dimension(:,:) :: AREG4,ENERGY4,CONSRV4,ATPE4
      REAL*4,allocatable, dimension(:,:,:,:) :: WAVE4, AIJK4,AIJL4,AFLX4
#ifndef NO_HDIURN
      ! REAL*4 HDIURN4(NDIUVAR,NDIUPT,HR_IN_MONTH)
      REAL*4,allocatable ::  HDIURN4(:,:,:)
#endif
      integer monac1(12),i_ida,i_xtra,it_check
!@var Kcomb counts acc-files as they are added up
      INTEGER, SAVE :: Kcomb=0

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER idac1(12)
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "DIAG01"
!@var header_subdd string label for subdaily accumulation arrays
      character(len=80) :: header_subdd
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
      ! REAL*8 :: AFLX_ST_GLOB(LM+LM_REQ+1,IM,JM,5)
      REAL*8,allocatable :: AFLX_ST_GLOB(:,:,:,:)
      real(kind=8),allocatable,dimension(:,:) :: P_acc_glob,PM_acc_glob
      real(kind=8),allocatable,dimension(:,:) :: R_acc_glob

      INTEGER :: J_0, J_1

      CALL GET( grid, J_STRT=J_0, J_STOP=J_1  )

      if(kradia.gt.0) then
        write (MODULE_HEADER(LHEAD+1:80),'(a6,i8,a20,i3,a7)')
     *   '#acc(=',idacc(2),') R8:SU.SD.TU.TD.dT(',lm+lm_req+1,',ijM,5)'

        IF (AM_I_ROOT()) then
           allocate (AFLX_ST_glob(LM+LM_REQ+1,IM,JM,5))
        else
           allocate (AFLX_ST_glob(LM+LM_REQ+1,1,1,5))
        end if

        SELECT CASE (IACTION)
        CASE (IOWRITE)            ! output to standard restart file
          CALL PACK_COLUMN(grid, AFLX_ST, AFLX_ST_glob)
          IF (AM_I_ROOT()) THEN
            WRITE (kunit,err=10) MODULE_HEADER,idacc(2),AFLX_ST_glob,it
          END IF
        CASE (IOWRITE_SINGLE)     ! output in single precision
          MODULE_HEADER(LHEAD+18:LHEAD+18) = '4'
          MODULE_HEADER(LHEAD+44:80) = ',monacc(12)'
          CALL PACK_COLUMN(grid, AFLX_ST, AFLX_ST_glob)
          IF (AM_I_ROOT()) THEN
            WRITE (kunit,err=10) MODULE_HEADER,idacc(2),
     *           REAL(AFLX_ST_glob,KIND=4), monacc,it
          END IF
        CASE (IOWRITE_MON)        ! output to end-of-month restart file
          MODULE_HEADER(LHEAD+1:80) = 'itime '
          IF (AM_I_ROOT()) THEN
            WRITE (kunit,err=10) MODULE_HEADER,it
          END IF
        CASE (ioread)           ! input from restart file
          if (AM_I_ROOT()) THEN
           READ (kunit,err=10) HEADER,idacc(2),AFLX_ST_glob,it
           IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
           END IF
          endif
          CALL UNPACK_COLUMN(grid, AFLX_ST_glob, AFLX_ST)
          CALL ESMF_BCAST(grid, idacc)
          CALL ESMF_BCAST(grid, it   )

        CASE (IOREAD_SINGLE)      !
!ESMF-- Allow all processes to read to avoid scattering monac1 and idac1.
          CALL PACK_COLUMN(grid, AFLX_ST, AFLX_ST_glob)
          if (AM_I_ROOT()) then
             allocate (AFLX4(LM+LM_REQ+1,IM,JM,5))
             READ (kunit,err=10) HEADER,idac1(2),AFLX4,monac1
             AFLX_ST_glob = AFLX_ST_glob + AFLX4
             deallocate (AFLX4)
          end if
          CALL UNPACK_COLUMN(grid, AFLX_ST_glob, AFLX_ST)
          CALL ESMF_BCAST(grid,idac1)
          CALL ESMF_BCAST(grid,monac1)
          IDACC(2) = IDACC(2) + IDAC1(2)
          monacc = monacc + monac1
        END SELECT

        deallocate (AFLX_ST_glob)

        return
      end if

C**** The regular model (Kradia le 0)
      write (MODULE_HEADER(LHEAD+1:LHEAD+15),'(a10,i4,a1)')
     *   'I/R8 keys(',1+NKEYNR*NKEYMO,')'             ! keyct,keynr(:,:)
      i_ida = Lhead + 10+4+1 + 10+2+1 + 1
      write (MODULE_HEADER(LHEAD+10+4+1+1:i_ida-1),'(a10,i2,a1)')
     *   ',TSFR(IJM,',KTSF,')'
      write (MODULE_HEADER(i_ida:i_ida+9),'(a7,i2,a1)')
     *   ',idacc(',nsampl,')'
      write (MODULE_HEADER(i_ida+9+1:i_ida+9 + 5+8+1),'(a5,i8,a1)')
     *   ',acc(',kacc,')'
      i_xtra = i_ida+9 + 5+8+1 + 1

      call alloc_ijdiag_glob

      allocate(P_acc_glob(im,jm),PM_acc_glob(im,jm))
      allocate(R_acc_glob(im,jm))

      SELECT CASE (IACTION)
      CASE (IOWRITE)            ! output to standard restart file
        write (MODULE_HEADER(i_xtra:80),             '(a7,i2,a)')
     *   ',x(IJM,',KTD+KOA,')'  ! make sure that i_xtra+7+2 < 80

        Call Gather_Diagnostics()

        If (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,keyct,KEYNR,TSFREZ,
     *     idacc, AJ,AREG,AJL,ASJL,AIJ,
     *     AIJL, ENERGY,CONSRV,
     *     SPECA,ATPE,ADIURN,WAVE,AGC,AIJK,AISCCP,
#ifndef NO_HDIURN
     *     HDIURN,
#endif
     *     TDIURN_glob,OA_glob,it
        END IF

c**** write accumulation arrays for subdaily diagnostics
        call pack_data(grid,P_acc,P_acc_glob)
        call pack_data(grid,PM_acc,PM_acc_glob)
        call pack_data(grid,R_acc,R_acc_glob)
        header_subdd='accumulation variables for subdaily diagnostics'
        if (am_i_root()) write(kunit,err=10) header_subdd,P_acc_glob
     &       ,PM_acc_glob, R_acc_glob

      CASE (IOWRITE_SINGLE)     ! output in single precision
        MODULE_HEADER(LHEAD+1:LHEAD+4) = 'I/R4'
        MODULE_HEADER(i_xtra:80) = ',monacc(12)'

        Call Gather_Diagnostics()

        If (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,
     *     keyct,KEYNR,REAL(TSFREZ,KIND=4),   idacc,
     *     REAL(AJ,KIND=4),REAL(AREG,KIND=4),
     *     REAL(AJL,KIND=4),REAL(ASJL,KIND=4),
     *     REAL(AIJ,KIND=4),REAL(AIJL,KIND=4),
     *     REAL(ENERGY,KIND=4), REAL(CONSRV,KIND=4),
     *     REAL(SPECA,KIND=4),REAL(ATPE,KIND=4),REAL(ADIURN,KIND=4),
     *     REAL(WAVE,KIND=4),REAL(AGC,KIND=4),
     *     REAL(AIJK,KIND=4),REAL(AISCCP,KIND=4),
#ifndef NO_HDIURN
     *     REAL(HDIURN,KIND=4),
#endif
     *     monacc,it
        END IF
      CASE (IOWRITE_MON)        ! output to end-of-month restart file
        MODULE_HEADER(i_ida:80) = ',it '

        CALL PACK_DATA(grid, TSFREZ_loc, TSFREZ)
        If (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,keyct,KEYNR,TSFREZ,it
        END IF
      CASE (ioread)           ! input from restart file
        if (AM_I_ROOT()) Then
          READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZ,
     *       idacc, AJ,AREG,AJL,ASJL,AIJ,AIJL,
     *       ENERGY,CONSRV,SPECA,ATPE,ADIURN,WAVE,AGC,AIJK,AISCCP,
#ifndef NO_HDIURN
     *       HDIURN,
#endif
     *       TDIURN_glob,OA_glob,it
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
        END IF

        Call BCAST_Scalars()
        Call Scatter_Diagnostics()

c**** read accumulation arrays for subdaily diagnostics
        if (am_i_root()) read(kunit,err=10) header_subdd,P_acc_glob
     &       ,PM_acc_glob,R_acc_glob
        call unpack_data(grid,P_acc_glob,P_acc)
        call unpack_data(grid,PM_acc_glob,PM_acc)
        call unpack_data(grid,R_acc_glob,R_acc)

      CASE (IOREAD_SINGLE)      !
        call Gather_Diagnostics()  ! to keep global arrays in up-to-date
        If (AM_I_ROOT()) Then
          call alloc_diag_r4
          READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZ4,
     *         idac1, AJ4,AREG4,AJL4,ASJL4,AIJ4,AIJL4,ENERGY4
     *         ,CONSRV4,SPECA4,ATPE4,ADIURN4,WAVE4,AGC4,AIJK4,AISCCP4,
#ifndef NO_HDIURN
     *       HDIURN4,
#endif
     *       monac1,it_check
          if(it.ne.it_check) then
            PRINT*,"io_diags: compare aj,aj4, ... dimensions"
            GO TO 10            ! or should that be just a warning ??
          end if

        ! copy or add in to full precision global variables
        ! First "non-distributed" arrays
          ENERGY=ENERGY+ENERGY4
          SPECA=SPECA+SPECA4 ; ATPE=ATPE+ATPE4 ; ADIURN=ADIURN+ADIURN4
          WAVE=WAVE+WAVE4
          AISCCP=AISCCP+AISCCP4
#ifndef NO_HDIURN
          HDIURN=HDIURN+HDIURN4
#endif
          AREG  = AREG   + AREG4
          CONSRV= CONSRV + CONSRV4
          ! Now for the global versions of the distributed arrays
          TSFREZ= TSFREZ4       ! not accumulated
          AJ    = AJ     + AJ4
          AJL   = AJL    + AJL4
          ASJL  = ASJL   + ASJL4
          AIJ   = AIJ    + AIJ4
          AGC   = AGC    + AGC4
          AIJK  = AIJK   + AIJK4
          AIJL  = AIJL   + AIJL4

          IDACC = IDACC + IDAC1
          call dealloc_diag_r4
        End If

        ! Send to other processors - update distributed arrays
        Call BCAST_Scalars()

        Call Scatter_Diagnostics()

!@var idacc(5) is the length of a time series (daily energy history).
!****   If combining acc-files, rather than concatenating these series,
!****   we average their beginnings (up to the length of the shortest)
        Kcomb = Kcomb + 1          ! reverse addition, take min instead
        if (Kcomb.gt.1) IDACC(5) = MIN(IDACC(5)-IDAC1(5),IDAC1(5))
        monacc = monacc + monac1
      CASE (irerun)      ! only keynr,tsfrez needed at beg of acc-period
        If (AM_I_ROOT()) Then
          READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZ  ! 'it' not read
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        End If
        CALL ESMF_BCAST(grid, keyct )
        CALL ESMF_BCAST(grid, KEYNR )
        CALL UNPACK_DATA(grid,  TSFREZ, TSFREZ_loc)
      END SELECT

      call dealloc_ijdiag_glob

      deallocate(P_acc_glob,PM_acc_glob)
      deallocate(R_acc_glob)
      RETURN

 10   IOERR=1
      call dealloc_ijdiag_glob
      deallocate(P_acc_glob,PM_acc_glob)
      deallocate(R_acc_glob)
      RETURN

      Contains

      Subroutine BCAST_Scalars()
        CALL ESMF_BCAST(grid, keyct )
        CALL ESMF_BCAST(grid, KEYNR )
        CALL ESMF_BCAST(grid, idacc )
        CALL ESMF_BCAST(grid, ENERGY)
        CALL ESMF_BCAST(grid, SPECA )
        CALL ESMF_BCAST(grid, ATPE  )
c        CALL ESMF_BCAST(grid, ADIURN)
        CALL ESMF_BCAST(grid, WAVE  )
        CALL ESMF_BCAST(grid, AISCCP)
#ifndef NO_HDIURN
c        CALL ESMF_BCAST(grid, HDIURN)
#endif
        CALL ESMF_BCAST(grid, it    )
      End Subroutine BCAST_Scalars

      Subroutine Scatter_Diagnostics()
        call scatter_zonal_diags()
        CALL UNPACK_DATA(grid,  TSFREZ, TSFREZ_loc)
        CALL UNPACK_DATA(grid,  AIJ,    AIJ_loc)
        CALL UNPACK_DATA(grid,  AIJK,   AIJK_loc)
        CALL UNPACK_DATA(grid,  AIJL,   AIJL_loc)
        CALL UNPACK_DATA(grid,  TDIURN_glob, TDIURN)
        CALL UNPACK_DATA(grid,  OA_glob,     OA)
      End Subroutine Scatter_Diagnostics

      Subroutine alloc_diag_r4
        allocate (AJ4(JM_BUDG,KAJ,NTYPE),AREG4(NREG,KAJ))
        allocate (AJL4(JM_BUDG,LM,KAJL),ASJL4(JM_BUDG,LM_REQ,KASJL))
        allocate (AIJ4(IM,JM,KAIJ))
        allocate (AIJL4(IM,JM,LM,KAIJL),ENERGY4(NEHIST,HIST_DAYS))
        allocate (CONSRV4(JM_BUDG,KCON),SPECA4(IMLONH+1,KSPECA,NSPHER))
        allocate (ATPE4(KTPE,NHEMI),ADIURN4(NDIUVAR,NDIUPT,HR_IN_DAY))
        allocate (WAVE4(RE_AND_IM,Max12HR_sequ,NWAV_DAG,KWP))
        allocate (AGC4(JM,LM,KAGC),AIJK4(IM,JM,LM,KAIJK))
        allocate (AISCCP4(ntau,npres,nisccp))
        allocate (TSFREZ4(IM,JM,KTSF))
#ifndef NO_HDIURN
        allocate (HDIURN4(NDIUVAR,NDIUPT,HR_IN_MONTH))
#endif
      End Subroutine alloc_diag_r4

      Subroutine dealloc_diag_r4
        deallocate (AJ4,AREG4,  AJL4,ASJL4,AIJ4,AIJL4)
        deallocate (ENERGY4,CONSRV4,SPECA4,ATPE4,ADIURN4,WAVE4)
        deallocate (AGC4,AIJK4, AISCCP4,TSFREZ4)
#ifndef NO_HDIURN
        deallocate (HDIURN4)
#endif
      End Subroutine dealloc_diag_r4


      END SUBROUTINE io_diags

      Subroutine Gather_Diagnostics()
      use domain_decomp_1d, only : grid,pack_data
      use diag_com
      implicit none
      call gather_zonal_diags
      call collect_scalars
      CALL PACK_DATA(grid,  TSFREZ_loc, TSFREZ)
      CALL PACK_DATA(grid,  AIJ_loc,    AIJ)
      CALL PACK_DATA(grid,  AIJK_loc,   AIJK)
      CALL PACK_DATA(grid,  AIJL_loc,   AIJL)
      CALL PACK_DATA(grid,  TDIURN, TDIURN_glob)
      CALL PACK_DATA(grid,  OA, OA_glob)
      return
      End Subroutine Gather_Diagnostics

      Subroutine Collect_Scalars()
      use precision_mod
      use domain_decomp_atm, only : grid,sumxpe,am_i_root
      use diag_com
      implicit none
      integer :: k
      real*8, dimension(:,:), allocatable :: consrv_sv
      CALL SUMXPE(AREG_loc, AREG, increment=.true.)
      AREG_loc(:,:)=0.
      CALL SUMXPE(ADIURN_loc, ADIURN, increment=.true.)
      ADIURN_loc=0
#ifndef NO_HDIURN
      CALL SUMXPE(HDIURN_loc, HDIURN, increment=.true.)
      HDIURN_loc=0
#endif
      CALL SUMXPE(AISCCP_loc, AISCCP, increment=.true.)
      AISCCP_loc=0

      if(am_i_root()) then
c for reproducibility on different numbers of processors
        call reduce_precision(areg,1d-9)
        call reduce_precision(aisccp,1d-9)
      endif

c CONSRV is a mixture of accumulations and instantaneous values, hence the
c more complicated logic
      if(am_i_root()) then
        allocate(consrv_sv(jm_budg,kcon)); consrv_sv = consrv ! pde hack
        do k=1,kcon
          if(ia_con(k).eq.ia_inst) consrv(:,k)=0. ! collect instantaneous values
        enddo
      endif
      CALL SUMXPE(CONSRV_loc, CONSRV, increment=.true.)
      do k=1,kcon
        if(ia_con(k).ne.ia_inst) consrv_loc(:,k)=0. ! keep instantaneous values
      enddo
      if(am_i_root()) then ! pde hack
        do k=1,kcon
          if(ia_con(k).eq.ia_inst .and. all(consrv(:,k)==0.))
     &         consrv(:,k) = consrv_sv(:,k)
        enddo
        deallocate(consrv_sv)
      endif

      return
      End Subroutine Collect_Scalars

      Subroutine Gather_zonal_diags()
      use domain_decomp_atm, only : grid
      use diag_com
      use diag_zonal, only : pack_lc
      implicit none
      call pack_lc(grid, AJ_loc,     AJ)
      call pack_lc(grid, AJL_loc,    AJL)
      call pack_lc(grid, ASJL_loc,   ASJL)
      call pack_lc(grid, AGC_loc,    AGC)
      return
      End Subroutine Gather_zonal_diags

      Subroutine Scatter_zonal_diags()
      use domain_decomp_atm, only : grid
      use diag_com
      use diag_zonal, only : unpack_lc
      implicit none
      call unpack_lc  (grid, aj,     aj_loc)
      call unpack_lc  (grid, ajl,    ajl_loc)
      call unpack_lc  (grid, asjl,   asjl_loc)
      call unpack_lc  (grid, agc,    agc_loc)
      return
      End Subroutine Scatter_zonal_diags

      SUBROUTINE aPERIOD (JMON1,JYR1,months,years,moff,  aDATE,LDATE)
!@sum  aPERIOD finds a 7 or 12-character name for an accumulation period
!@+   if the earliest month is NOT the beginning of the 2-6 month period
!@+   the name will reflect that fact ONLY for 2 or 3-month periods
!@auth Reto A. Ruedy
!@ver  1.0
      USE MODEL_COM, only : AMONTH
      implicit none
!@var JMON1,JYR1 month,year of beginning of period 1
      INTEGER JMON1,JYR1
!@var JMONM,JMONL middle,last month of period
      INTEGER JMONM,JMONL
!@var months,years length of 1 period,number of periods
      INTEGER months,years
!@var moff = # of months from beginning of period to JMON1 if months<12
      integer moff
!@var yr1,yr2 (end)year of 1st and last period
      INTEGER yr1,yr2
!@var aDATE date string: MONyyr1(-yyr2)
      character*12 aDATE
!@var LDATE length of date string (7 or 12)
      INTEGER LDATE

      LDATE = 7                  ! if years=1
      if(years.gt.1) LDATE = 12

      aDATE(1:12)=' '
      aDATE(1:3)=AMONTH(JMON1)        ! letters 1-3 of month IF months=1
      yr1=JYR1
      JMONL=JMON1+months-1
      if(JMONL.GT.12) then
         yr1=yr1+1
         JMONL=JMONL-12
      end if
      if (moff.gt.0.and.months.le.3) then  ! earliest month is NOT month
        JMONL = 1 + mod(10+jmon1,12)       ! 1 of the 2-3 month period
        yr1=JYR1
        if (jmon1.gt.1) yr1=yr1+1
      end if
      yr2=yr1+years-1
      write(aDATE(4:7),'(i4.4)') yr1
      if(years.gt.1) write(aDATE(8:12),'(a1,i4.4)') '-',yr2

      if(months.gt.12) aDATE(1:1)='x'                ! should not happen
      if(months.le.1 .or. months.gt.12) return

!**** 1<months<13: adjust characters 1-3 of aDATE (=beg) if necessary:
!**** beg=F?L where F/L=letter 1 of First/Last month for 2-11 mo.periods
!****    =F+L                                        for 2 month periods
!****    =FML where M=letter 1 of Middle month       for 3 month periods
!****    =FnL where n=length of period if n>3         4-11 month periods
      aDATE(3:3)=AMONTH(JMONL)(1:1)            ! we know: months>1
      IF (months.eq.2) then
        aDATE(2:2)='+'
        return
      end if
      if (months.eq.3) then
        JMONM = JMONL-1
        if (moff.eq.1) jmonm = jmon1+1
        if (jmonm.gt.12) jmonm = jmonm-12
        if (jmonm.le.0 ) jmonm = jmonm+12
        aDATE(2:2)=AMONTH(JMONM)(1:1)
        return
      end if
      if (moff.gt.0) then  ! can't tell non-consec. from consec. periods
        jmon1 = jmon1-moff
        if (jmon1.le.0) jmon1 = jmon1+12
        JMONL=JMON1+months-1
        if (jmonl.gt.12) jmonl = jmonl-12
        aDATE(1:1)=AMONTH(JMON1)(1:1)
        aDATE(3:3)=AMONTH(JMONL)(1:1)
      end if
      IF (months.ge.4.and.months.le.9) write (aDATE(2:2),'(I1)') months
      IF (months.eq.10) aDATE(2:2)='X'         ! roman 10
      IF (months.eq.11) aDATE(2:2)='B'         ! hex   11
      IF (months.eq.6) THEN                    !    exceptions:
         IF (JMON1.eq. 5) aDATE(1:3)='NHW'     ! NH warm season May-Oct
         IF (JMON1.eq.11) aDATE(1:3)='NHC'     ! NH cold season Nov-Apr
      END IF
      IF (months.eq.7) THEN                    !    to avoid ambiguity:
         IF (JMON1.eq. 1) aDATE(1:3)='J7L'     ! Jan-Jul J7J->J7L
         IF (JMON1.eq. 7) aDATE(1:3)='L7J'     ! Jul-Jan J7J->L7J
      END IF
      IF (months.eq.12) THEN
C****    beg=ANn where the period ends with month n if n<10 (except 4)
         aDATE(1:3)='ANN'                      ! regular annual mean
         IF (JMONL.le. 9) WRITE(aDATE(3:3),'(I1)') JMONL
         IF (JMONL.eq. 4) aDATE(1:3)='W+C'     ! NH warm+cold seasons
         IF (JMONL.eq.10) aDATE(1:3)='C+W'     ! NH cold+warm seasons
         IF (JMONL.eq.11) aDATE(1:3)='ANM'     ! meteor. annual mean
      END IF
      return
      end SUBROUTINE aPERIOD

C**** Routines associated with the budget grid

      SUBROUTINE SET_J_BUDG
!@sum set_j_budg definition for grid points map to budget-grid zonal means
!@auth Gavin Schmidt
      USE GEOM, only : j_budg,j_0b,j_1b,lat2d_dg
      USE DIAG_COM, only : jm_budg
      USE DOMAIN_DECOMP_ATM, only :GRID,GET
      IMPLICIT NONE
!@var I,J are atm grid point values for the accumulation
      INTEGER :: I,J,J_0,J_1,I_0,I_1,J_0H,J_1H,I_0H,I_1H
      INTEGER :: IER

C**** define atmospheric grid
      CALL GET(grid, J_STRT=J_0,J_STOP=J_1, J_STRT_HALO=J_0H,
     *     J_STOP_HALO=J_1H  )
      I_0H = grid%I_STRT_HALO ; I_1H = grid%I_STOP_HALO
      I_0 = grid%I_STRT ; I_1 = grid%I_STOP

      ALLOCATE( J_BUDG(I_0H:I_1H, J_0H:J_1H), STAT = IER)

C**** Define mapping from actual lon/lat point to budget grid
C**** this should be valid for all grids (lat/lon, cubed sphere,...)
      DO J=J_0H,J_1H
        DO I=I_0,I_1
           J_BUDG(I,J)=NINT(1+(lat2d_dg(I,J)+90)*(JM_BUDG-1)/180.)
        END DO
      END DO

C**** define limits on budget indices for each processor
      j_0b=MINVAL( J_BUDG(I_0:I_1,J_0:J_1) )
      j_1b=MAXVAL( J_BUDG(I_0:I_1,J_0:J_1) )

      RETURN
      END SUBROUTINE SET_J_BUDG

      subroutine set_wtbudg()
!@sum Precomputes area weights for zonal means on budget grid
!auth Denis Gueyffier
      USE GEOM, only : j_budg, axyp, imaxj
      USE MODEL_COM, only : fim
      USE DIAG_COM, only : jm_budg,wtbudg,wtbudg2,lat_budg,dxyp_budg
      USE DOMAIN_DECOMP_ATM, only :GRID,GET
      IMPLICIT NONE
      INTEGER :: I,J,J_0,J_1,I_0,I_1
      INTEGER :: IER
      real*8 :: dlat_budg,fjeq_budg

      CALL GET(grid, J_STRT=J_0,J_STOP=J_1)
      I_0 = grid%I_STRT ; I_1 = grid%I_STOP

      ALLOCATE( wtbudg(I_0:I_1, J_0:J_1), STAT = IER)  !deallocated near very end, stays in memory all the time
      ALLOCATE( wtbudg2(I_0:I_1, J_0:J_1), STAT = IER)

C**** get the nominal latitudes of the budget grid
      dlat_budg = 180./REAL(JM_budg)  ! for full polar box
      if(jm_budg.eq.46) dlat_budg=4. ! 1/2 box at pole for 4 deg res.
      if(jm_budg.eq.24) dlat_budg=8. ! 1/4 box at pole for 8 deg res.
      lat_budg(1) = -90.; lat_budg(jm_budg) = +90.;
      fjeq_budg = .5*(1+jm_budg)
      do j=2,jm_budg-1
        lat_budg(j) = dlat_budg*(j-fjeq_budg)
      enddo

      call set_budg_area()

C**** Compute area weights of zonal bands, including adjustment for polar
C**** box in lat/lon case.
      do J=J_0,J_1
         do I=I_0,I_1
            wtbudg2(I,J) = axyp(I,J)/dxyp_budg(J_BUDG(I,J))
#ifdef CUBED_SPHERE
            wtbudg (I,J) = wtbudg2(I,J)
#else
            wtbudg (I,J) = wtbudg2(I,J)*fim/real(imaxj(j),kind=8)
#endif
         enddo
      enddo

      RETURN
      END SUBROUTINE set_wtbudg

      subroutine set_budg_area()
!@sum  pre-computes area of budget-grid band (zig-zag band on the cubed sphere)
!@+    accross processors. Several processors can contribute to same band
!@auth Denis Gueyffier
      use GEOM, only: J_BUDG,axyp
      use DIAG_COM, only : dxyp_budg,dxyp_budg_loc,JM_BUDG
      USE DOMAIN_DECOMP_ATM, only :grid,GET,sumxpe,esmf_bcast
      IMPLICIT NONE
      INTEGER :: I,J,J_0,J_1,I_0,I_1
      logical :: increment

      CALL GET(grid, J_STRT=J_0,J_STOP=J_1)
      I_0 = grid%I_STRT ; I_1 = grid%I_STOP

      dxyp_budg_loc(:)=0.0

      do J=J_0,J_1
         do I=I_0,I_1
            dxyp_budg_loc(J_BUDG(I,J))=dxyp_budg_loc(J_BUDG(I,J))
     &       +axyp(I,J)
         enddo
      enddo

      increment = .false.
      call SUMXPE(dxyp_budg_loc,dxyp_budg,increment)   !summing global area
      call esmf_bcast(grid,dxyp_budg)

      return
      end subroutine set_budg_area

      SUBROUTINE INC_AJ(I,J,ITYPE,J_DIAG,ACC)
!@sum inc_aj grid dependent incrementer for zonal mean budget diags
!@auth Gavin Schmidt
      USE DIAG_COM, only : aj=>aj_loc, wtbudg
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J are atm grid point values for the accumulation
!@var ITYPE is the surface type
      INTEGER, INTENT(IN) :: I,J,ITYPE
!@var J_DIAG is placing of diagnostic element
      INTEGER, INTENT(IN) :: J_DIAG
!@var ACC value of the diagnostic to be accumulated
      REAL*8, INTENT(IN) :: ACC

C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reproducible for MPI).
      !   wtbudg area-weight =1 on lat-lon, <1 on cubed sphere
      AJ(J_BUDG(I,J),J_DIAG,ITYPE) = AJ(J_BUDG(I,J),J_DIAG,ITYPE)
     &     +wtbudg(I,J)*ACC

      RETURN
      END SUBROUTINE INC_AJ

      SUBROUTINE INC_AREG(I,J,JR,J_DIAG,ACC)
!@sum inc_areg incrementer for regional budget diags
!@auth Gavin Schmidt
      USE DIAG_COM, only : areg=>areg_loc,wtbudg
      USE GEOM, only : j_budg,axyp
      IMPLICIT NONE
!@var I,J are atm grid point values for the accumulation
!@var JR is the region
      INTEGER, INTENT(IN) :: I,J,JR
!@var J_DIAG identifier of diagnostic element
      INTEGER, INTENT(IN) :: J_DIAG
!@var ACC value of the diagnostic to be accumulated
      REAL*8, INTENT(IN) :: ACC

      AREG(JR,J_DIAG) = AREG(JR,J_DIAG)+ACC*AXYP(I,J)

      RETURN
      END SUBROUTINE INC_AREG

      SUBROUTINE INC_AJL(I,J,L,JL_INDEX,ACC)
!@sum inc_ajl adds ACC located at atmospheric gridpoint I,J,L
!@+   to the latitude-height zonal sum AJL(J,L,JL_INDEX).
!@auth M. Kelley
      USE DIAG_COM, only : ajl=>ajl_loc,wtbudg
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: JL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC

      AJL(J_BUDG(I,J),L,JL_INDEX) = AJL(J_BUDG(I,J),L,JL_INDEX)
     &     +wtbudg(I,J)*ACC

      RETURN
      END SUBROUTINE INC_AJL

      SUBROUTINE INC_AJL2(I,J,L,JL_INDEX,ACC)
c temporary variant of inc_ajl without any weighting
      USE DIAG_COM, only : ajl=>ajl_loc
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: JL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC

      AJL(J_BUDG(I,J),L,JL_INDEX) = AJL(J_BUDG(I,J),L,JL_INDEX)
     &     +ACC

      RETURN
      END SUBROUTINE INC_AJL2

      SUBROUTINE INC_ASJL(I,J,L,SJL_INDEX,ACC)
!@sum inc_asjl adds ACC located at atmospheric gridpoint I,J,L
!@+   to the latitude-height zonal sum ASJL(J,L,JL_INDEX).
!@+   This is a trivial version for the latlon grid.
!@auth M. Kelley
      USE DIAG_COM, only : asjl=>asjl_loc,wtbudg
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L
!@var SJL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: SJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC

      ASJL(J_BUDG(I,J),L,SJL_INDEX) = ASJL(J_BUDG(I,J),L,SJL_INDEX)
     *     +wtbudg(I,J)*ACC

      RETURN
      END SUBROUTINE INC_ASJL

#ifdef NEW_IO
      subroutine def_rsf_acc(fid,r4_on_disk)
!@sum  def_rsf_acc defines accumulation array structure in restart/acc files
!@auth M. Kelley
!@ver  beta
      use model_com, only : idacc
      use diag_com, only : monacc,
     &     aj=>aj_ioptr,areg=>areg_ioptr,agc=>agc_ioptr,
     &     aij=>aij_loc,aijl=>aijl_loc,aijk=>aijk_loc, ! dist
     &     oa,tdiurn,                                  ! dist
     &     ajl,asjl,consrv,
     &     speca,atpe,adiurn,energy,wave,aisccp
#ifndef NO_HDIURN
      use diag_com, only :  hdiurn
#endif
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_attr
      implicit none
      integer fid   !@var fid file id
      logical :: r4_on_disk  !@var r4_on_disk if true, real*8 stored as real*4
      real*8 :: r8dum
      character(len=20) :: tmpstr

      call defvar(grid,fid,idacc,'monacc(twelve)')
      call defvar(grid,fid,idacc,'idacc(nsampl)')
      call defvar(grid,fid,tdiurn,'tdiurn(dist_im,dist_jm,ktd)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,oa,'oa(dist_im,dist_jm,koa)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,aij,'aij(dist_im,dist_jm,kaij)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,aijl,'aijl(dist_im,dist_jm,lm,kaijl)',
     &     r4_on_disk=r4_on_disk)
#ifndef CUBED_SPHERE
      call defvar(grid,fid,aijk,'aijk(dist_im,dist_jm,lm,kaijk)',
     &     r4_on_disk=r4_on_disk)
#endif

#ifdef CACHED_SUBDD
      if(.not.r4_on_disk) call def_rsf_subdd_acc(fid,r4_on_disk)
#endif

      call defvar(grid,fid,aj,'aj(jm_budg,kaj,ntype)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,ajl,'ajl(jm_budg,lm,kajl)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,asjl,'asjl(jm_budg,lm_req,kasjl)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,agc,'agc(jmlat,lm,kagc)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,consrv,'consrv(jm_budg,kcon)',
     &     r4_on_disk=r4_on_disk)

      call defvar(grid,fid,areg,'areg(nreg,kaj)',r4_on_disk=r4_on_disk)

      call defvar(grid,fid,energy,'energy(nehist,hist_days)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,speca,
     &     'speca(imlonh_plus_1,kspeca,nspher)',r4_on_disk=r4_on_disk)
      call defvar(grid,fid,atpe,'atpe(ktpe,nhemi)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,wave,
     &     'wave(re_and_im,max12hr_sequ,nwav_dag,kwp)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,aisccp,'aisccp(ntau,npres,nisccp)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,adiurn,
     &     'adiurn(ndiuvar,ndiupt,hr_in_day)',r4_on_disk=r4_on_disk)
#ifndef NO_HDIURN
      call defvar(grid,fid,hdiurn,
     &     'hdiurn(ndiuvar,ndiupt,hr_in_month)',r4_on_disk=r4_on_disk)
#endif

      if(.not.r4_on_disk) then          ! reproducibility info
        tmpstr = 'is_npes_reproducible' ! for checkpoint files
        call defvar(grid,fid,r8dum,tmpstr)
        call write_attr(grid,fid,tmpstr,'aisccp' ,'no')
        call write_attr(grid,fid,tmpstr,'areg'   ,'no')
        call write_attr(grid,fid,tmpstr,'speca'  ,'no')
        call write_attr(grid,fid,tmpstr,'energy' ,'no')
c       if(using Ocean R not on atm grid)
c       call write_attr(grid,fid,tmpstr,'consrv' ,'no')
#ifdef CUBED_SPHERE
        call write_attr(grid,fid,tmpstr,'aj'     ,'no')
        call write_attr(grid,fid,tmpstr,'ajl'    ,'no')
        call write_attr(grid,fid,tmpstr,'asjl'   ,'no')
        call write_attr(grid,fid,tmpstr,'agc'    ,'no')
        call write_attr(grid,fid,tmpstr,'consrv' ,'no')
#ifdef TRACERS_ON /* move these declarations to the tracer code? */
        call write_attr(grid,fid,tmpstr,'tajln'  ,'no')
        call write_attr(grid,fid,tmpstr,'tajls'  ,'no')
        call write_attr(grid,fid,tmpstr,'tconsrv','no')
#endif
#endif
      endif

      return
      end subroutine def_rsf_acc

      subroutine new_io_acc(fid,iaction)
!@sum  new_io_acc read/write accumulation arrays from/to restart/acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : im,jm,lm,ioread,iowrite,iowrite_single,
     &     idacc
c i/o pointers point to:
c    primary instances of arrays when writing restart files
c    extended/rescaled instances of arrays when writing acc files
      use diag_com, only : monacc,kaijl,
     &     aj=>aj_ioptr,areg=>areg_ioptr,agc=>agc_ioptr,
     &     aij=>aij_loc,aijl=>aijl_loc,aijk=>aijk_loc, ! dist
     &     oa,tdiurn,                                  ! dist
     &     ajl,asjl,consrv,
     &     speca,atpe,adiurn,energy,wave,aisccp
#ifndef NO_HDIURN
      use diag_com, only :  hdiurn
#endif
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : hasNorthPole, hasSouthPole
      use pario, only : write_dist_data,read_dist_data
     &     ,write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: k,l
      select case (iaction)
      case (iowrite,iowrite_single) ! output to restart or acc file
        if(iaction.eq.iowrite) then ! already done for iowrite_single
          call gather_zonal_diags
          call collect_scalars
        else  ! pole fills needed for acc-files
          do k=1,kaijl
            do l=1,lm
              if(hasSouthPole(grid)) then
                aijl(2:im, 1,l,k) = aijl(1, 1,l,k)
              endif
              if(hasNorthPole(grid)) then
                aijl(2:im,jm,l,k) = aijl(1,jm,l,k)
              endif
            enddo
          enddo
        endif
        call write_data(grid,fid,'monacc',monacc)
        call write_data(grid,fid,'idacc',idacc)
        call write_data(grid,fid,'energy',energy)
        call write_data(grid,fid,'speca',speca)
        call write_data(grid,fid,'atpe',atpe)
        call write_data(grid,fid,'wave',wave)
        call write_data(grid,fid,'aisccp',aisccp)
        call write_data(grid,fid,'adiurn',adiurn)
#ifndef NO_HDIURN
        call write_data(grid,fid,'hdiurn',hdiurn)
#endif
        call write_dist_data(grid,fid,'tdiurn',tdiurn)
        call write_dist_data(grid,fid,'oa',oa)
        call write_dist_data(grid,fid,'aij',aij)
        call write_dist_data(grid,fid,'aijl',aijl)
#ifndef CUBED_SPHERE
        call write_dist_data(grid,fid,'aijk',aijk)
#endif

        call write_data(grid,fid,'aj',aj)
        call write_data(grid,fid,'ajl',ajl)
        call write_data(grid,fid,'asjl',asjl)
        call write_data(grid,fid,'agc',agc)
        call write_data(grid,fid,'consrv',consrv)

        call write_data(grid,fid,'areg',areg)

      case (ioread)            ! input from restart or acc file
c for which scalars is bcast_all=.true. necessary?
        call read_data(grid,fid,'monacc',monacc,bcast_all=.true.)
        call read_data(grid,fid,'idacc',idacc,bcast_all=.true.)
        call read_data(grid,fid,'energy',energy,bcast_all=.true.)
        call read_data(grid,fid,'speca',speca,bcast_all=.true.)
        call read_data(grid,fid,'atpe',atpe,bcast_all=.true.)
        call read_data(grid,fid,'wave',wave,bcast_all=.true.)
        call read_data(grid,fid,'aisccp',aisccp,bcast_all=.true.)
        call read_data(grid,fid,'adiurn',adiurn,bcast_all=.true.)
#ifndef NO_HDIURN
        call read_data(grid,fid,'hdiurn',hdiurn,bcast_all=.true.)
#endif
        call read_dist_data(grid,fid,'tdiurn',tdiurn)
        call read_dist_data(grid,fid,'oa',oa)
        call read_dist_data(grid,fid,'aij',aij)
        call read_dist_data(grid,fid,'aijl',aijl)
#ifndef CUBED_SPHERE
        call read_dist_data(grid,fid,'aijk',aijk)
#endif

        call read_data(grid,fid,'aj',aj)
        call read_data(grid,fid,'ajl',ajl)
        call read_data(grid,fid,'asjl',asjl)
        call read_data(grid,fid,'agc',agc)
        call read_data(grid,fid,'consrv',consrv)
        call scatter_zonal_diags

        call read_data(grid,fid,'areg',areg)

      end select

#ifdef CACHED_SUBDD
      if(iaction.eq.iowrite) then
        call write_subdd_accdata(fid,iaction)
      endif
#endif

      return
      end subroutine new_io_acc

      subroutine def_rsf_longacc(fid,r4_on_disk)
!@sum  def_rsf_longacc defines accumulation array structure in restart/acc files
!@auth M. Kelley
!@ver  beta
      use diag_com, only : tsfrez=>tsfrez_loc,keyct,keynr
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid           !@var fid file id
      logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4
      call defvar(grid,fid,keyct,'keyct')
      call defvar(grid,fid,keynr,'keynr(nkeynr,nkeymo)')
      call defvar(grid,fid,tsfrez,'tsfrez(dist_im,dist_jm,ktsf)',
     &     r4_on_disk=r4_on_disk)
      return
      end subroutine def_rsf_longacc

      subroutine new_io_longacc(fid,iaction)
!@sum  new_io_longacc read/write accumulation arrays from/to restart+acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use diag_com, only : tsfrez=>tsfrez_loc,keyct,keynr
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
     &     ,write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart or acc file
        call write_data(grid,fid,'keyct',keyct)
        call write_data(grid,fid,'keynr',keynr)
        call write_dist_data(grid,fid,'tsfrez',tsfrez)
      case (ioread)            ! input from restart or acc file
c for which scalars is bcast_all=.true. necessary?
        call read_data(grid,fid,'keyct',keyct,bcast_all=.true.)
        call read_data(grid,fid,'keynr',keynr,bcast_all=.true.)
        call read_dist_data(grid,fid,'tsfrez',tsfrez)
      end select
      return
      end subroutine new_io_longacc

c def_rsf_subdd
      subroutine def_rsf_subdd(fid)
!@sum def_rsf_subdd defines write/read structure of accumulation arrays
!@+                 for subdaily diagnostics to/from restart files
!@auth Jan Perlwitz

      use domain_decomp_atm, only: grid
      use pario, only: defvar
      use diag_com, only: P_acc,PM_acc,R_acc

      implicit none

      integer,intent(in) :: fid

      call defvar(grid,fid,P_acc,'P_acc(dist_im,dist_jm)')
      call defvar(grid,fid,PM_acc,'PM_acc(dist_im,dist_jm)')
      call defvar(grid,fid,R_acc,'R_acc(dist_im,dist_jm)')

      return
      end subroutine def_rsf_subdd

c new_io_subdd
      subroutine new_io_subdd(fid,iaction)
!@sum new_io_subdd write/read of accumulation arrays for subdaily diagnostics
!+                 to/from restart files
!@auth Jan Perlwitz

      use domain_decomp_atm, only: grid
      use pario, only: read_dist_data,write_dist_data
      use model_com, only: ioread,iowrite
      use diag_com, only: P_acc,PM_acc,R_acc

      implicit none

      integer,intent(in) :: iaction,fid

      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'P_acc',P_acc)
        call write_dist_data(grid,fid,'PM_acc',PM_acc)
        call write_dist_data(grid,fid,'R_acc',R_acc)

      case (ioread)             ! input from restart file
        call read_dist_data(grid,fid,'P_acc',P_acc)
        call read_dist_data(grid,fid,'PM_acc',PM_acc)
        call read_dist_data(grid,fid,'R_acc',R_acc)

      end select

      return
      end subroutine new_io_subdd

      subroutine def_meta_atmacc(fid)
!@sum  def_meta_atmacc defines metadata in atm acc files
!@auth M. Kelley
!@ver  beta
      use model_com, only : im
      use diag_com, only : kagc,
     &     ia_j,ia_jl,ia_ij,ia_ijl,ia_con,ia_gc,ia_ijk,
     &     name_j,name_reg,sname_jl,name_ij,name_ijl,name_dd,
     &     name_consrv,sname_gc,name_ijk,
     &     cdl_j,cdl_reg,cdl_jl,
     &     cdl_ij,cdl_ijl,cdl_ij_latlon,cdl_ijl_latlon,
     &     cdl_dd,cdl_hd,cdl_consrv,cdl_gc,cdl_ijk,
     &     hemis_j,hemis_jl,vmean_jl,hemis_consrv,hemis_gc,vmean_gc,
     &     hemis_ij,
     &     scale_j,scale_jl,scale_ij,scale_ijl,scale_dd,scale_con,
     &     scale_gc,scale_ijk,
     &     iden_j,iden_reg,denom_jl,denom_ijl,denom_ij,denom_dd,
     &     denom_gc,denom_ijk,
     &     lm,
     &     isccp_press,isccp_tau,isccp_late,wisccp
      use geom, only : axyp
#ifdef CUBED_SPHERE
      use geom, only : lon2d_dg,lat2d_dg,lonbds,latbds
#endif
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_attr
      use cdl_mod, only : defvar_cdl
      implicit none
      integer :: fid         !@var fid file id
      integer :: int_dummy

#ifdef CUBED_SPHERE
      call defvar(grid,fid,lon2d_dg,'lon(dist_im,dist_jm)')
      call defvar(grid,fid,lat2d_dg,'lat(dist_im,dist_jm)')
      call defvar(grid,fid,lonbds,'lonbds(four,dist_im,dist_jm)')
      call defvar(grid,fid,latbds,'latbds(four,dist_im,dist_jm)')
#endif

      call defvar(grid,fid,axyp,'axyp(dist_im,dist_jm)')

      call write_attr(grid,fid,'aj','reduction','sum')
      call write_attr(grid,fid,'aj','split_dim',2)
      call defvar(grid,fid,hemis_j,'hemis_aj(shnhgm,kaj,ntype)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_aj','reduction','sum')
      call defvar(grid,fid,ia_j,'ia_aj(kaj)')
      call defvar(grid,fid,scale_j,'scale_aj(kaj)')
      call defvar(grid,fid,iden_j,'denom_aj(kaj)')
      call defvar(grid,fid,name_j,'sname_aj(sname_strlen,kaj)')
      call defvar_cdl(grid,fid,cdl_j,'cdl_aj(cdl_strlen,kcdl_aj)')

      call write_attr(grid,fid,'areg','reduction','sum')
      call write_attr(grid,fid,'areg','split_dim',2)
      call defvar(grid,fid,ia_j,'ia_areg(kaj)')
      call defvar(grid,fid,scale_j,'scale_areg(kaj)')
      call defvar(grid,fid,iden_reg,'denom_areg(kaj)')
      call defvar(grid,fid,name_reg,'sname_areg(sname_strlen,kaj)')
      call defvar_cdl(grid,fid,cdl_reg,
     &     'cdl_areg(cdl_strlen,kcdl_areg)')

      call write_attr(grid,fid,'consrv','reduction','sum')
      call write_attr(grid,fid,'consrv','split_dim',2)
      call defvar(grid,fid,hemis_consrv,'hemis_consrv(shnhgm,kcon)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_consrv','reduction','sum')
      call defvar(grid,fid,ia_con,'ia_consrv(kcon)')
      call defvar(grid,fid,scale_con,'scale_consrv(kcon)')
      call defvar(grid,fid,name_consrv,
     &     'sname_consrv(sname_strlen,kcon)')
      call defvar_cdl(grid,fid,cdl_consrv,
     &     'cdl_consrv(cdl_strlen,kcdl_consrv)')

      call write_attr(grid,fid,'ajl','reduction','sum')
      call write_attr(grid,fid,'ajl','split_dim',3)
      call defvar(grid,fid,hemis_jl,'hemis_ajl(shnhgm,lm,kajl)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_ajl','reduction','sum')
      call defvar(grid,fid,vmean_jl,'vmean_ajl(jm_budg_plus3,one,kajl)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'vmean_ajl','reduction','sum')
      call defvar(grid,fid,ia_jl,'ia_ajl(kajl)')
      call defvar(grid,fid,scale_jl,'scale_ajl(kajl)')
      call defvar(grid,fid,denom_jl,'denom_ajl(kajl)')
      call defvar(grid,fid,sname_jl,'sname_ajl(sname_strlen,kajl)')
      call defvar_cdl(grid,fid,cdl_jl,'cdl_ajl(cdl_strlen,kcdl_ajl)')

      call write_attr(grid,fid,'agc','reduction','sum')
      call write_attr(grid,fid,'agc','split_dim',3)
      call defvar(grid,fid,hemis_gc,'hemis_agc(shnhgm,lm,kagc)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_agc','reduction','sum')
      call defvar(grid,fid,vmean_gc,'vmean_agc(jmlat_plus3,one,kagc)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'vmean_agc','reduction','sum')
      call defvar(grid,fid,ia_gc(1:kagc),'ia_agc(kagc)')
      call defvar(grid,fid,scale_gc(1:kagc),'scale_agc(kagc)')
      call defvar(grid,fid,denom_gc(1:kagc),'denom_agc(kagc)')
      call defvar(grid,fid,sname_gc(1:kagc),
     &     'sname_agc(sname_strlen,kagc)')
      call defvar_cdl(grid,fid,cdl_gc,
     &     'cdl_agc(cdl_strlen,kcdl_agc)')

      call write_attr(grid,fid,'aij','reduction','sum')
      call write_attr(grid,fid,'aij','split_dim',3)
      call defvar(grid,fid,ia_ij,'ia_aij(kaij)')
      call defvar(grid,fid,scale_ij,'scale_aij(kaij)')
      call defvar(grid,fid,denom_ij,'denom_aij(kaij)')
      call defvar(grid,fid,name_ij,'sname_aij(sname_strlen,kaij)')
      call defvar_cdl(grid,fid,cdl_ij,'cdl_aij(cdl_strlen,kcdl_aij)')
      call defvar(grid,fid,hemis_ij,'hemis_aij(one,shnhgm,kaij)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_aij','reduction','sum')
#ifdef CUBED_SPHERE
      call defvar_cdl(grid,fid,cdl_ij_latlon,
     &     'cdl_aij_latlon(cdl_strlen,kcdl_aij_latlon)')
#endif

      call write_attr(grid,fid,'aijl','reduction','sum')
      call write_attr(grid,fid,'aijl','split_dim',4)
      call defvar(grid,fid,ia_ijl,'ia_aijl(kaijl)')
      call defvar(grid,fid,scale_ijl,'scale_aijl(kaijl)')
      call defvar(grid,fid,denom_ijl,'denom_aijl(kaijl)')
      call defvar(grid,fid,name_ijl,'sname_aijl(sname_strlen,kaijl)')
      call defvar_cdl(grid,fid,cdl_ijl,
     &     'cdl_aijl(cdl_strlen,kcdl_aijl)')
#ifdef CUBED_SPHERE
      call defvar_cdl(grid,fid,cdl_ijl_latlon,
     &     'cdl_aijl_latlon(cdl_strlen,kcdl_aijl_latlon)')
#endif

#ifndef CUBED_SPHERE
      call write_attr(grid,fid,'aijk','reduction','sum')
      call write_attr(grid,fid,'aijk','split_dim',4)
      call defvar(grid,fid,ia_ijk,'ia_aijk(kaijk)')
      call defvar(grid,fid,scale_ijk,'scale_aijk(kaijk)')
      call defvar(grid,fid,denom_ijk,'denom_aijk(kaijk)')
      call defvar(grid,fid,name_ijk,'sname_aijk(sname_strlen,kaijk)')
      call defvar_cdl(grid,fid,cdl_ijk,
     &     'cdl_aijk(cdl_strlen,kcdl_aijk)')
#endif

      call write_attr(grid,fid,'adiurn','reduction','sum')
      call write_attr(grid,fid,'adiurn','split_dim',1)
      call defvar(grid,fid,int_dummy,'ntime_adiurn')
      call write_attr(grid,fid,'ntime_adiurn','reduction','sum')
      call defvar(grid,fid,denom_dd,'denom_adiurn(ndiuvar)')
      call defvar(grid,fid,scale_dd,'scale_adiurn(ndiuvar)')
      call defvar(grid,fid,name_dd,'sname_adiurn(sname_strlen,ndiuvar)')
      call defvar_cdl(grid,fid,cdl_dd,
     &     'cdl_adiurn(cdl_strlen,kcdl_adiurn)')
#ifndef NO_HDIURN
      call write_attr(grid,fid,'hdiurn','split_dim',1)
      call defvar(grid,fid,int_dummy,'ntime_hdiurn')
      call defvar(grid,fid,denom_dd,'denom_hdiurn(ndiuvar)')
      call defvar(grid,fid,scale_dd,'scale_hdiurn(ndiuvar)')
      call defvar(grid,fid,name_dd,'sname_hdiurn(sname_strlen,ndiuvar)')
      call defvar_cdl(grid,fid,cdl_hd,
     &     'cdl_hdiurn(cdl_strlen,kcdl_hdiurn)')
#endif

      call write_attr(grid,fid,'aisccp','reduction','sum')
      call defvar(grid,fid,wisccp,'wisccp(nisccp)')
      call write_attr(grid,fid,'wisccp','reduction','sum')
      call defvar(grid,fid,isccp_press,'isccp_press(npres)')
      call defvar(grid,fid,isccp_tau,'isccp_tau(ntau)')
      call defvar(grid,fid,isccp_late,'isccp_late(nisccp_plus_1)')

      return
      end subroutine def_meta_atmacc

      subroutine write_meta_atmacc(fid)
!@sum  write_meta_atmacc write atm accumulation metadata to file
!@auth M. Kelley
      use model_com, only : im,nday,idacc
      use diag_com, only : kagc,
     &     ia_j,ia_jl,ia_ij,ia_ijl,ia_con,ia_gc,ia_ijk,
     &     name_j,name_reg,sname_jl,name_ij,name_ijl,name_dd,
     &     name_consrv,sname_gc,name_ijk,
     &     cdl_j,cdl_reg,cdl_jl,
     &     cdl_ij,cdl_ijl,cdl_ij_latlon,cdl_ijl_latlon,
     &     cdl_dd,cdl_hd,cdl_consrv,cdl_gc,cdl_ijk,
     &     hemis_j,hemis_jl,vmean_jl,hemis_consrv,hemis_gc,vmean_gc,
     &     hemis_ij,
     &     scale_j,scale_jl,scale_ij,scale_ijl,scale_dd,scale_con,
     &     scale_gc,scale_ijk,
     &     iden_j,iden_reg,denom_jl,denom_ij,denom_ijl,denom_dd,
     &     denom_gc,denom_ijk,
     &     lm,ia_12hr,
     &     isccp_press,isccp_tau,isccp_late,wisccp
      use geom, only : axyp
#ifdef CUBED_SPHERE
      use geom, only : lon2d_dg,lat2d_dg,lonbds,latbds
#endif
      use domain_decomp_atm, only : grid
      use pario, only : write_data,write_dist_data
      use cdl_mod, only : write_cdl
      implicit none
      integer fid   !@var fid unit number of read/write
      integer :: i,n,ntime_dd,ntime_hd

#ifdef CUBED_SPHERE
      call write_dist_data(grid,fid,'lon',lon2d_dg)
      call write_dist_data(grid,fid,'lat',lat2d_dg)
      call write_dist_data(grid,fid,'lonbds',lonbds,jdim=3)
      call write_dist_data(grid,fid,'latbds',latbds,jdim=3)
#endif

      call write_dist_data(grid,fid,'axyp',axyp)

      call write_data(grid,fid,'hemis_aj',hemis_j)
      call write_data(grid,fid,'ia_aj',ia_j)
      call write_data(grid,fid,'scale_aj',scale_j)
      call write_data(grid,fid,'denom_aj',iden_j)
      call write_data(grid,fid,'sname_aj',name_j)
      call write_cdl(grid,fid,'cdl_aj',cdl_j)

      call write_data(grid,fid,'ia_areg',ia_j)
      call write_data(grid,fid,'scale_areg',scale_j)
      call write_data(grid,fid,'denom_areg',iden_reg)
      call write_data(grid,fid,'sname_areg',name_reg)
      call write_cdl(grid,fid,'cdl_areg',cdl_reg)

      call write_data(grid,fid,'hemis_consrv',hemis_consrv)
      call write_data(grid,fid,'ia_consrv',ia_con)
      call write_data(grid,fid,'scale_consrv',scale_con)
      call write_data(grid,fid,'sname_consrv',name_consrv)
      call write_cdl(grid,fid,'cdl_consrv',cdl_consrv)

      call write_data(grid,fid,'hemis_ajl',hemis_jl)
      call write_data(grid,fid,'vmean_ajl',vmean_jl)
      call write_data(grid,fid,'ia_ajl',ia_jl)
      call write_data(grid,fid,'scale_ajl',scale_jl)
      call write_data(grid,fid,'denom_ajl',denom_jl)
      call write_data(grid,fid,'sname_ajl',sname_jl)
      call write_cdl(grid,fid,'cdl_ajl',cdl_jl)

      call write_data(grid,fid,'hemis_agc',hemis_gc)
      call write_data(grid,fid,'vmean_agc',vmean_gc)
      call write_data(grid,fid,'ia_agc',ia_gc(1:kagc))
      call write_data(grid,fid,'scale_agc',scale_gc(1:kagc))
      call write_data(grid,fid,'denom_agc',denom_gc(1:kagc))
      call write_data(grid,fid,'sname_agc',sname_gc(1:kagc))
      call write_cdl(grid,fid,'cdl_agc',cdl_gc)

      call write_data(grid,fid,'hemis_aij',hemis_ij)
      call write_data(grid,fid,'ia_aij',ia_ij)
      call write_data(grid,fid,'scale_aij',scale_ij)
      call write_data(grid,fid,'denom_aij',denom_ij)
      call write_data(grid,fid,'sname_aij',name_ij)
      call write_cdl(grid,fid,'cdl_aij',cdl_ij)
#ifdef CUBED_SPHERE
      call write_cdl(grid,fid,'cdl_aij_latlon',cdl_ij_latlon)
#endif

      call write_data(grid,fid,'ia_aijl',ia_ijl)
      call write_data(grid,fid,'scale_aijl',scale_ijl)
      call write_data(grid,fid,'denom_aijl',denom_ijl)
      call write_data(grid,fid,'sname_aijl',name_ijl)
      call write_cdl(grid,fid,'cdl_aijl',cdl_ijl)
#ifdef CUBED_SPHERE
      call write_cdl(grid,fid,'cdl_aijl_latlon',cdl_ijl_latlon)
#endif

#ifndef CUBED_SPHERE
      call write_data(grid,fid,'ia_aijk',ia_ijk)
      call write_data(grid,fid,'scale_aijk',scale_ijk)
      call write_data(grid,fid,'denom_aijk',denom_ijk)
      call write_data(grid,fid,'sname_aijk',name_ijk)
      call write_cdl(grid,fid,'cdl_aijk',cdl_ijk)
#endif

      ntime_dd = (idacc(ia_12hr)/2)*(nday/24)
      call write_data(grid,fid,'ntime_adiurn',ntime_dd)
      call write_data(grid,fid,'scale_adiurn',scale_dd)
      call write_data(grid,fid,'denom_adiurn',denom_dd)
      call write_data(grid,fid,'sname_adiurn',name_dd)
      call write_cdl(grid,fid,'cdl_adiurn',cdl_dd)

#ifndef NO_HDIURN
      ntime_hd = nday/24
      call write_data(grid,fid,'ntime_hdiurn',ntime_hd)
      call write_data(grid,fid,'scale_hdiurn',scale_dd)
      call write_data(grid,fid,'denom_hdiurn',denom_dd)
      call write_data(grid,fid,'sname_hdiurn',name_dd)
      call write_cdl(grid,fid,'cdl_hdiurn',cdl_hd)
#endif

      call write_data(grid,fid,'isccp_press',isccp_press)
      call write_data(grid,fid,'isccp_tau',isccp_tau)
      call write_data(grid,fid,'isccp_late',isccp_late)
      call write_data(grid,fid,'wisccp',wisccp)

      return
      end subroutine write_meta_atmacc

      subroutine set_ioptrs_atmacc_default
c point i/o pointers for diagnostic accumlations to the
c instances of the arrays used during normal operation.
      use diag_com
      implicit none
      aj_ioptr     => aj
      areg_ioptr   => areg
      agc_ioptr    => agc
      return
      end subroutine set_ioptrs_atmacc_default

      subroutine set_ioptrs_atmacc_extended
c point i/o pointers for diagnostic accumlations to the
c instances of the arrays containing derived outputs
      use diag_com
      implicit none
      aj_ioptr     => aj_out
      areg_ioptr   => areg_out
      agc_ioptr    => agc_out
      return
      end subroutine set_ioptrs_atmacc_extended

#endif /* NEW_IO */
