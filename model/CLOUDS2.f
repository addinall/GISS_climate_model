#include "rundeck_opts.h"

#define ALT_CDNC_INPUTS

      MODULE CLOUDS

!@sum  CLOUDS column physics of moist conv. and large-scale condensation
!@auth M.S.Yao/A. Del Genio (modifications by Gavin Schmidt)
!@cont MSTCNV,LSCOND,ANVIL_OPTICAL_THICKNESS,MC_CLOUD_FRACTION,
!@+    CONVECTIVE_MICROPHYSICS,MC_PRECIP_PHASE,MASS_FLUX,PRECIP_MP
      USE CONSTANT, only : rgas,grav,lhe,lhs,lhm,sha,bysha,pi,by6
     *     ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd,rhow,twopi
      USE MODEL_COM, only : lm,dtsrc,itime,coupled_chem
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
      USE CONSTANT, only : kapa,mair,gasc
      USE MODEL_COM, only : ptop,psf,ls1,sig,sige
#endif
#ifdef SCM
      USE MODEL_COM, only : I_TARG,J_TARG,NSTEPSCM
      USE SCMCOM, only: SG_P,SCM_SAVE_T,SCM_SAVE_Q,SCM_DEL_T,
     &                  SCM_DEL_Q,SCM_ATURB_FLAG,iu_scm_prt,NRINIT,
     &                  SG_T,SG_Q 
      USE SCMDIAG, only : WCUSCM,WCUALL,WCUDEEP,PRCCDEEP,NPRCCDEEP,
     &           MPLUMESCM,MPLUMEALL,MPLUMEDEEP,ENTSCM,ENTALL,
     &           ENTDEEP,DETRAINDEEP,TPALL,PRCCGRP,PRCCICE,MCCOND,
     &           PRESAV,LHPSAV,PREMC,LHPMC,CUMFLX,DWNFLX,SCM_LWP_MC,
     &           SCM_LWP_SS,SCM_IWP_MC,SCM_IWP_SS,SCM_WM_MC,PLUME_MAX,
     &           PLUME_MIN,SCM_GZ,SCM_H,SCM_HSAT,ARM_H,ARM_HSAT
#endif

      USE CLOUDS_COM, only : ncol

#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
      USE mo_bulk2m_driver_gcm, ONLY: execute_bulk2m_driver
#endif
      USE QUSDEF, only : nmom,xymoms,zmoms,zdir
#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm,trname,t_qlimit,ntm_soa
#ifdef TRACERS_WATER
     &     ,nGAS, nPART, nWATER, tr_wd_TYPE, tr_RKD, tr_DHD,
     *     tr_evap_fact
     &     ,gases_list,gases_count
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
     &     ,aqchem_list,aqchem_count
#endif
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,Ntm_dust
#endif
#endif
#endif
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
#ifdef TRACERS_AMP
      USE CLOUDS_COM, only: NACTC,NAERC
      USE AERO_CONFIG, only: NMODES
#endif
#endif
      IMPLICIT NONE
      SAVE
C**** parameters and constants
      REAL*8, PARAMETER :: TI=233.16d0   !@param TI pure ice limit
      REAL*8, PARAMETER :: CLDMIN=.10d0 !@param CLDMIN min MC/LSC region
!@param WMU critical cloud water content for rapid conversion (g m**-3)
      REAL*8, PARAMETER :: WMU=.25
      REAL*8, PARAMETER :: WMUL=.5       !@param WMUL WMU over land
C     REAL*8, PARAMETER :: WMUI=.1d0     !@param WMUI WMU for ice clouds
      REAL*8 WMUI                          !@param WMUI WMU for ice clouds
      REAL*8 WMUSI        !@param WMUSI WMU for liquid clouds over sea-ice
      REAL*8, PARAMETER :: BRCLD=.2d0    !@param BRCLD for cal. BYBR
      REAL*8, PARAMETER :: FDDET=.25d0 !@param FDDET remainder of downdraft
      REAL*8, PARAMETER :: DTMIN1=1.d0 !@param DTMIN1 min DT to stop downdraft drop
      REAL*8, PARAMETER :: SLHE=LHE*BYSHA
      REAL*8, PARAMETER :: SLHS=LHS*BYSHA
!@param CCMUL multiplier for convective cloud cover
!@param CCMUL1 multiplier for deep anvil cloud cover
!@param CCMUL2 multiplier for shallow anvil cloud cover
!@param COETAU multiplier for convective cloud optical thickness
      REAL*8, PARAMETER :: CCMUL=2.,CCMUL1=5.,CCMUL2=3.,COETAU=.08d0

      REAL*8 :: RTEMP,CMX,RCLDX,WMUIX,CONTCE1,CONTCE2,TNX,QNX
      REAL*8 :: BYBR,BYDTsrc,XMASS,PLAND
!@var BYBR factor for converting cloud particle radius to effect. radius
!@var XMASS dummy variable
!@var PLAND land fraction

C**** Set-able variables
!@dbparam LMCM max level for originating MC plumes
      INTEGER :: LMCM = -1 ! defaults to LS1-1 if not set in rundeck
!@dbparam ISC integer to turn on computation of stratocumulus clouds
      INTEGER :: ISC = 0  ! set ISC=1 to compute stratocumulus clouds
C     REAL*8 :: U00MAX = .99d0      ! maximum U00 for water clouds
C****
C**** WARNING: U00wtrX AND U00ice ARE NO LONGER USED BY THE GCM. USE U00a and U00b INSTEAD
C****
!@dbparam U00wtrX multiplies U00ice for critical humidity for water clds
      REAL*8 :: U00wtrX = 1.0d0     ! default, needed for AR4 runs
!@dbparam U00ice critical humidity for ice cloud condensation
      REAL*8 :: U00ice = .7d0       ! default, needed for AR4 runs
!@dbparam U00a tuning knob for U00 above 850 mb without moist convection
!@dbparam U00b tuning knob for U00 below 850 mb and in convective regions
      REAL*8 :: U00a = 0.55d0       ! default
      REAL*8 :: U00b = 1.00d0       ! default
!@dbparam funio_denominator funio denominator
      REAL*8 :: funio_denominator=22.d0  ! default
!@dbparam autoconv_multiplier autoconversion rate multiplier
      REAL*8 :: autoconv_multiplier=1.d0 ! default
!@dbparam radius_multiplier cloud particle radius multiplier
      REAL*8 :: radius_multiplier=1.d0   ! default
!@dbparam wmui_multiplier critical ice cloud water multiplier
      REAL*8 :: wmui_multiplier=1.d0     ! default
!@dbparam entrainment_cont1 constant for entrainment rate, plume 1
      REAL*8 :: entrainment_cont1=.3d0   ! default
!@dbparam entrainment_cont2 constant for entrainment rate, plume 2
      REAL*8 :: entrainment_cont2=.6d0   ! default
!@dbparam HRMAX maximum distance an air parcel rises from surface
      REAL*8 :: HRMAX = 1000.d0     ! default (m)
!@dbparam RIMAX maximum ice cloud size
!@dbparam RWMAX maximum water cloud size
      REAL*8 :: RIMAX = 100.d0, RWMAX = 20.d0      ! microns
!@dbparam RWCldOX multiplies part.size of water clouds over ocean
      REAL*8 :: RWCldOX=1.d0
!@dbparam RICldX multiplies part.size of ice clouds at 1000mb
!@+       RICldX changes linearly to 1 as p->0mb
      REAL*8 :: RICldX=1.d0 , xRICld
!@dbparam do_blU00 =1 if boundary layer U00 is treated differently
      INTEGER :: do_blU00=0     ! default is to disable this

#ifdef TRACERS_ON
!@var ntx,NTIX: Number and Indices of active tracers used in convection
      integer, dimension(ntm) :: ntix
      integer ntx
#endif
C**** ISCCP diag related variables
!!! @parameter ncol used to be set here  20 for gcm runs and 100 for scm runs
!!!     now moved to CLOUDS_COM.f  for portability
!@var tautab look-up table to convert count value to optical thickness
!@var invtau look-up table to convert optical thickness to count value
      real*8 :: tautab(0:255)
      integer :: invtau(-20:45000)
#ifdef SCM
      real*8 scmqsat,armqsat
#endif

C**** input variables
      LOGICAL DEBUG
!@var RA ratio of primary grid box to secondary gridbox
      REAL*8, DIMENSION(:), ALLOCATABLE :: RA !(KMAX)
!@var UM,VM,UM1,VM1,U_0,V_0 velocity related variables(UM,VM)=(U,V)*AIRM
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: UM,VM,UM1,VM1 !(KMAX,LM)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: U_0,V_0       !(KMAX,LM)

!@var Miscellaneous vertical arrays set in driver
!@var PLE pressure at layer edge
!@var LHP array of precip phase ! may differ from LHX
      REAL*8, DIMENSION(LM+1) :: PLE,LHP
      REAL*8, DIMENSION(LM) :: PL,PLK,AIRM,BYAM,ETAL,TL,QL,TH,RH,WMX
     *     ,VSUBL,MCFLX,DGDSM,DPHASE,DTOTW,DQCOND,DGDQM,AQ,DPDT,RH1
     *     ,FSSL,VLAT,DDMFLX,WTURB,TVL,W2L,GZL
     *     ,SAVWL,SAVWL1,SAVE1L,SAVE2L,DPHASHLW,DPHADEEP,DGSHLW,DGDEEP
     *     ,QDNL,TDNL,U00L
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
     *     ,WMXICE
#endif
!@var PL layer pressure (mb)
!@var PLK PL**KAPA
!@var AIRM the layer's pressure depth (mb)
!@var BYAM 1./AIRM
!@var ETAL fractional entrainment rate
!@var TL, QL temperature, specific humidity of the layer
!@var TH potential temperature (K)
!@var RH relative humidity
!@var RH1 relative humidity to compare with the threshold humidity
!@var WMX cloud water mixing ratio (kg/kg)
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
!@var WMXICE ice water mixing ratio (kg/kg)
#endif
!@var VSUBL downward vertical velocity due to cumulus subsidence (cm/s)
!@var MCFLX, DGDSM, DPHASE, DQCOND, DGDQM dummy variables
!@var DDMFLX accumulated downdraft mass flux (mb)
!@var AQ time change rate of specific humidity (s**-1)
!@var DPDT time change rate of pressure (mb/s)
!@var FSSL grid fraction for large-scale clouds
!@var VLAT dummy variable
      REAL*8, DIMENSION(LM+1) :: PRECNVL
!@var WTURB turbulent vertical velocity (m)
!@var PRECNVL convective precip entering the layer top
C**** new arrays must be set to model arrays in driver (before MSTCNV)
      REAL*8, DIMENSION(LM) :: SDL,WML
!@var SDL vertical velocity in sigma coordinate
!@var WML cloud water mixing ratio (kg/kg)
C**** new arrays must be set to model arrays in driver (after MSTCNV)
      REAL*8, DIMENSION(LM) :: TAUMCL,SVLATL,CLDMCL,SVLHXL,SVWMXL,SVLAT1
!@var TAUMCL convective cloud optical thickness
!@var SVLATL saved LHX for convective cloud
!@var CLDMCL convective cloud cover
!@var SVLHXL saved LHX for large-scale cloud
!@var SVWMXL saved detrained convective cloud water
      REAL*8, DIMENSION(LM) :: CSIZEL
!@var CSIZEL cloud particle radius (micron)
#ifdef CLD_AER_CDNC
      REAL*8, DIMENSION(LM) :: ACDNWM,ACDNIM
!@var ACDNWM,ACDNIM -CDNC - warm and cold moist cnv clouds (cm^-3)
      REAL*8, DIMENSION(LM) :: ACDNWS,ACDNIS
!@var ACDNWS,ACDNIS -CDNC - warm and cold large scale clouds (cm^-3)
      REAL*8, DIMENSION(LM) :: AREWS,AREIS,AREWM,AREIM  ! for diag
!@var AREWS and AREWM are moist cnv, and large scale Reff arrays (um)
      REAL*8, DIMENSION(LM) :: ALWWS,ALWIS,ALWWM,ALWIM  ! for diag
      REAL*8, DIMENSION(LM) :: CDN3DL,CRE3DL
!@var ALWWM and ALWIM  etc are liquid water contents
!@var SMLWP is LWP
      REAL*8 SMLWP
!@var SME is the TKE in 1 D from e(l) = egcm(l,i,j)  (m^2/s^2)
      REAL*8, DIMENSION(LM)::SME
      INTEGER NLSW,NLSI,NMCW,NMCI
#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
!@var CTTEM,CD3DL,CL3DL,CI3DL are cld temp, cld thickness,cld water
      REAL*8, DIMENSION(LM) ::CTEML,CD3DL,CL3DL,CI3DL
#endif
C**** new arrays must be set to model arrays in driver (before LSCOND)
      REAL*8, DIMENSION(LM) :: TTOLDL,CLDSAVL,CLDSV1
!@var TTOLDL previous potential temperature
!@var CLDSAVL saved large-scale cloud cover
#ifdef CLD_AER_CDNC
      REAL*8, DIMENSION(LM)::OLDCDL,OLDCDI
!@var OLDCDL is saved CDNC
!@var OLDCDI is saved ice crystal numbe
#endif
C**** new arrays must be set to model arrays in driver (after LSCOND)
      REAL*8, DIMENSION(LM) :: SSHR,DCTEI,TAUSSL,CLDSSL
!@var SSHR,DCTEI height diagnostics of dry and latent heating by MC
!@var TAUSSL large-scale cloud optical thickness
!@var CLDSSL large-scale cloud cover

!@var SM,QM Vertical profiles of (T/p**kappa)*AIRM, q*AIRM
      REAL*8, DIMENSION(LM) :: SM,QM
      REAL*8, DIMENSION(NMOM,LM) :: SMOM,QMOM,SMOMMC,QMOMMC,
     *  SMOMLS,QMOMLS

#ifdef TRACERS_ON
!@var TM Vertical profiles of tracers
      REAL*8, DIMENSION(LM,NTM) :: TM
      REAL*8, DIMENSION(nmom,lm,ntm) :: TMOM
!@var TRDNL tracer concentration in lowest downdraft (kg/kg)
      REAL*8, DIMENSION(NTM,LM) :: TRDNL
      COMMON/CLD_TRCCOM/TM,TMOM,TRDNL
#ifdef TRACERS_WATER
!@var TRWML Vertical profile of liquid water tracers (kg)
!@var TRSVWML New liquid water tracers from m.c. (kg)
      REAL*8, DIMENSION(NTM,LM) :: TRWML, TRSVWML
!@var TRPRSS super-saturated tracer precip (kg)
!@var TRPRMC moist convective tracer precip (kg)
      REAL*8, DIMENSION(NTM)    :: TRPRSS,TRPRMC
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
c for diagnostics
      REAL*8, DIMENSION(NTM,LM) :: DT_SULF_MC,DT_SULF_SS
#endif
#ifdef TRDIAG_WETDEPO
!@dbparam diag_wetdep switches on/off special diags for wet deposition
      INTEGER :: diag_wetdep=0 ! =off (default) (on: 1)
!@var trcond_mc saves tracer condensation in MC clouds [kg]
!@var trdvap_mc saves tracers evaporated in downdraft of MC clouds [kg]
!@var trflcw_mc saves tracers condensed in cloud water of MC clouds [kg]
!@var trprcp_mc saves tracer precipitated from MC clouds [kg]
!@var trnvap_mc saves reevaporated tracer of MC clouds precip [kg]
!@var trwash_mc saves tracers washed out by collision for MC clouds [kg]
      REAL*8,DIMENSION(Lm,Ntm) :: trcond_mc,trdvap_mc,trflcw_mc,
     &     trprcp_mc,trnvap_mc,trwash_mc
!@var trwash_ls saves tracers washed out by collision for LS clouds [kg]
!@var trprcp_ls saves tracer precipitation from LS clouds [kg]
!@var trclwc_ls saves tracers condensed in cloud water of LS clouds [kg]
!@var trevap_ls saves reevaporated tracers of LS cloud precip [kg]
!@var trclwe_ls saves tracers evaporated from cloud water of LS clouds [kg]
!@var trcond_ls saves tracer condensation in LS clouds [kg]
      REAL*8,DIMENSION(Lm,Ntm) :: trwash_ls,trevap_ls,trclwc_ls,
     &     trprcp_ls,trclwe_ls,trcond_ls
#endif
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
!@var tm_dust vertical profile of dust/mineral tracers [kg]
      REAL*8,DIMENSION(Lm,Ntm_dust) :: tm_dust
!@var tmom_dust vertical profiles of dust/mineral tracer moments [kg]
      REAL*8,DIMENSION(nmom,Lm,Ntm_dust) :: tmom_dust
!@var trprc_dust dust/mineral tracer precip [kg]
      REAL*8,DIMENSION(Lm,Ntm_dust) :: trprc_dust
#endif
#endif
#ifdef TRACERS_WATER
      COMMON/CLD_WTRTRCCOM/TRWML, TRSVWML,TRPRSS,TRPRMC
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
     *     ,DT_SULF_MC,DT_SULF_SS
#endif
#ifdef TRDIAG_WETDEPO
     &     ,trcond_mc,trdvap_mc,trflcw_mc,trprcp_mc,trnvap_mc,trwash_mc
     &     ,trwash_ls,trevap_ls,trclwc_ls,trprcp_ls,trclwe_ls,trcond_ls
#endif
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      COMMON/CLD_PRECDUST/ tm_dust,tmom_dust,trprc_dust
#endif
#endif
#ifdef TRACERS_WATER
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
#endif
#endif
#endif

!@var KMAX index for surrounding velocity
!@var LP50 50mb level
      INTEGER ::  KMAX,LP50
!@var PEARTH fraction of land in grid box
!@var TS average surface temperture (C)
!@var RIS, RI1, RI2 Richardson numbers
      REAL*8 :: PEARTH,TS,QS,US,VS,RIS,RI1,RI2,DXYPIJ,ROICE
!@var DCL max level of planetary boundary layer
      INTEGER :: DCL

C**** output variables
      REAL*8 :: PRCPMC,PRCPSS,HCNDSS,WMSUM
!@var PRCPMC precip due to moist convection
!@var PRCPSS precip due to large-scale condensation
!@var HCNDSS heating due to large-scale condensation
!@var WMSUM cloud liquid water path
#ifdef CLD_AER_CDNC
      REAL*8 :: WMCLWP,WMCTWP
!@var WMCLWP , WMCTWP moist convective LWP and total water path
#endif
      REAL*8 :: CLDSLWIJ,CLDDEPIJ
!@var CLDSLWIJ shallow convective cloud cover
!@var CLDDEPIJ deep convective cloud cover
      INTEGER :: LMCMAX,LMCMIN
!@var LMCMAX upper-most convective layer
!@var LMCMIN lowerest convective layer
!@var AIRXL is convective mass flux (mb)
      REAL*8 AIRXL,PRHEAT
!@var RNDSSL stored random number sequences
      REAL*8  RNDSSL(3,LM)
!@var prebar1 copy of variable prebar
      REAL*8 prebar1(Lm+1)
#ifdef CFMIP
!@var CFMIP_CCL Mixing ratio of convective cloud [kg/kg]
      REAL*8 cfmip_ccl(LM)
!@var CFMIP_REFF Effective radius of large-scale cloud [um]
      REAL*8 cfmip_reff(LM)
!@var CFMIP_CCP Mixing ratio of convective precipitation [kg/kg]
      REAL*8 cfmip_ccp(LM+1)
#endif

      COMMON/CLDPRV/PLE,PL,PLK,AIRM,BYAM,ETAL
     *  ,TL,QL,TH,RH,WMX,VSUBL,MCFLX,SSHR,DGDSM,DPHASE,LHP
     *  ,DPHASHLW,DPHADEEP,DGSHLW,DGDEEP,SVLAT1
     *  ,DTOTW,DQCOND,DCTEI,DGDQM,DXYPIJ,DDMFLX,PLAND
     *  ,AQ,DPDT,PRECNVL,SDL,WML,SVLATL,SVLHXL,SVWMXL,CSIZEL,RH1
     *  ,TTOLDL,CLDSAVL,TAUMCL,CLDMCL,TAUSSL,CLDSSL,RNDSSL
     *  ,SM,QM,SMOM,QMOM,PEARTH,TS,QS,US,VS,RIS,RI1,RI2, AIRXL
     *  ,SMOMMC,QMOMMC,SMOMLS,QMOMLS,CLDSV1,PRHEAT,TDNL,QDNL,U00L
     *  ,PRCPMC,PRCPSS,HCNDSS,WMSUM,CLDSLWIJ,CLDDEPIJ,VLAT
#ifdef CLD_AER_CDNC
     *  ,ACDNWM,ACDNIM,ACDNWS,ACDNIS
     *  ,AREWM,AREIM,AREWS,AREIS,ALWIM,ALWWM
     *  ,OLDCDL,OLDCDI
     *  ,SME
     *  ,SMLWP,CDN3DL,CRE3DL
     *  ,WMCLWP,WMCTWP
#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
     *  ,CTEML,CD3DL,CL3DL,CI3DL
#endif
     *  ,TNX,QNX,RTEMP,CMX,RCLDX,WMUIX,CONTCE1,CONTCE2
     *  ,FSSL,WTURB,TVL,W2L,GZL
     *  ,SAVWL,SAVWL1,SAVE1L,SAVE2L
#ifdef CLD_AER_CDNC
     *  ,NLSW,NLSI,NMCW,NMCI
#endif
     *  ,prebar1,LMCMAX,LMCMIN,KMAX,DCL,DEBUG  ! int/logic last (alignment)

#ifdef TRACERS_ON
! The following tracer arrays are workspace for MSTCNV.  They are
! declared as permanent arrays here to avoid the expense of initializing
! temporary LM,NTM arrays to zero each time MSTCNV is called.  After
! completion of MC calculations, MSTCNV resets these arrays to zero in
! the layers in which they were used.
!@var DTM,DTMR: Vertical profiles of Tracers changes
      REAL*8, DIMENSION(LM,NTM)      :: DTM=0, DTMR=0, TMDNL=0
      REAL*8, DIMENSION(NMOM,LM,NTM) :: DTMOM=0, DTMOMR=0, TMOMDNL=0
!@var TPOLD saved plume temperature after condensation for tracers
!@+   (this is slightly different from TPSAV)
      REAL*8, DIMENSION(LM)       :: TPOLD=0
#ifdef TRACERS_WATER
!@var TRCOND tracer mass in condensate
!@var TRCONDV tracer mass in lofted condensate
      REAL*8, DIMENSION(NTM,LM)   :: TRCOND=0,TRCONDV=0
#endif
#endif

      CONTAINS

      SUBROUTINE MSTCNV(IERR,LERR,i_debug,j_debug)

!@sum  MSTCNV moist convective processes (precip, convective clouds,...)
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
!@calls adv1d,QSAT,DQSATDT,THBAR, MASS_FLUX, CONVECTIVE_MICROPHYSICS, PRECIP_MP,
!@+     MC_PRECIP_PHASE, MC_CLOUD_FRACTION, ANVIL_OPTICAL_THICKNESS

C**** FREE PARAMETERS THAT USERS MIGHT WANT TO VARY INCLUDE
C****    (1) ADJUSTMENT TIME FOR STABILIZATION OF CLOUD BASE BY CUMULUS MASS
C****        FLUX (TADJ): DEFAULT VALUE = 1.0 (HOUR)
C****    (2) SCALING FACTOR FOR ENTRAINMENT STRENGTH (CONTCE): DEFAULT VALUES
C****        CONTCE1 = 0.3, CONTCE2 = 0.6
C****    (3) SCALING FACTOR FOR EFFECT OF PRESSURE GRADIENT ON CONVECTIVE
C****        HORIZONTAL MOMENTUM TRANSPORT (PGRAD): DEFAULT VALUE = 0.7
C****    (4) FRACTION OF UPDRAFT MASS ASSUMED FOR INITIAL DOWNDRAFT MASS (ETADN):
C****        DEFAULT VALUE = 1./3.
C****    (5) FRACTION OF PRECIPITATING CONDENSATE AVAILABLE FOR RE-EVAPORATION
C****        RATHER THAN INCORPORATION INTO DOWNDRAFT (FDDRT): DEFAULT VALUE = 1.
C****
C**** NOTE: THESE PARAMETERS ARE INTENDED FOR RESEARCH PURPOSES ONLY AND
C**** SHOULD NOT BE USED TO ADJUST THE MODEL TO RADIATION BALANCE
C****

      IMPLICIT NONE
C
C****
C              *******************************
C              ***   DECLARATION SECTION   ***
C              *******************************
C                        FUNCTIONS
C                        PARAMETERS
C                        VARIABLES
C
C              *******************************
C
C                  ***               ***
C                  ***   FUNCTIONS   ***
C                  ***               ***
C
      REAL*8 ::   DQSATDT, PRECIP_MP, QSAT, THBAR
C
!@var DQSATDT     dQSAT/dT
!@var PRECIP_MP   mass density of precipitating condensate (kg/m^3)
!@var QSAT        saturation specific humidity
!@var THBAR       virtual temperature at layer edge
C
C              *******************************
C
C                  ***                ***
C                  ***   PARAMETERS   ***
C                  ***                ***
C
      REAL*8,   PARAMETER :: AIRM0=100.d0
      REAL*8,   PARAMETER :: CK1 = 1.,  CN0=8.d6,  CN0I=8.d6,  CN0G=8.d6
      INTEGER,  PARAMETER :: ITMAX=50
      REAL*8,   PARAMETER :: FITMAX=1d0/ITMAX
      REAL*8,   PARAMETER :: PN=1.d0,   RHOG=400., RHOIP=100.
C
!@param AIRM0 air mass used to compute convective cloud cover
!@param CK1 a tunning const.
!@param CN0, CN0I, CN0G intercepts of Marshall-Palmer particle size dist.
!@param ITMAX max iteration indices
!@param FITMAX set to 1/ITMAX
!@param PN tuning exponential for computing WV
!@param RHOG,RHOIP density of graupel and ice particles
C
#ifdef CLD_AER_CDNC
      INTEGER, PARAMETER :: SNTM=29 !17+ntm_soa/2  !for tracers for CDNC
#endif
C
C              *******************************
C
C                  ***               ***
C                  ***   VARIABLES   ***
C                  ***               ***
C
C              Sections:
C                        ARRAYS   - DECLARATIONS
C                                 - DEFINITIONS
C                        SCALARS  - DECLARATIONS
C                                 - DEFINITIONS
C                        VARIABLES FOR CONDITIONAL COMPILIATIONS
C
C
C        *********   ARRAY DECLARATIONS   **********
C
      REAL*8, DIMENSION(0:LM) :: CM
C
      REAL*8, DIMENSION(LM) ::
     *   BUOY,
     *   CCM,  CDHEAT, CMNEG, COND, CONDP, CONDP1, CONDV, CONDGP,
     *   CONDIP,
     *   DQM,  DSM,    DQMR,  DSMR, DET,   DM,     DMR,   DDM,    DDR,
ccc  *   ENT,  F,      HEAT1, LHP,  ML,    QM1,    QMT,   QMDNL,  QMOLD,
     *   ENT,  F,      HEAT1,       ML,    QM1,    QMT,   QMDNL,  QMOLD,
     *   SM1,  SMT,    SMDNL, SMOLD,TPSAV, TAUMC1, WCU,   WCU2
C
      REAL*8, DIMENSION(KMAX) ::
     *        SUMU,SUMV,SUMU1,SUMV1,UMP,VMP,UMDN,VMDN
C
      REAL*8, DIMENSION(KMAX,LM) :: DUM,DVM,UMDNL,VMDNL
C
      REAL*8, DIMENSION(NMOM) ::
     *        SMOMP,QMOMP, SMOMPMAX,QMOMPMAX, SMOMDN,QMOMDN
C
      REAL*8, DIMENSION(NMOM,LM) ::
     *     DQMOM, DQMOMR, DSMOM, DSMOMR, FMOM, QMOMDNL, QMOMOLD,
     *                                         SMOMDNL, SMOMOLD
C
C        *********   ARRAY DEFINITIONS   **********
C
!        (0:LM)
!@var CM            air mass of subsidence
!        (LM)
!@var BUOY
!@var CCM           convective plume mass (mb)
!@var CDHEAT        heating due to condensation
!@var CMNEG
!@var COND,CONDP, CONDP1, CONDV, CONDGP,CONDIP condensate mass density (kg/m^3)
!@var DQM,DSM,DQMR,DSMR Vertical profiles of T/Q and changes
!@var DDM           downdraft mass (mb)
!@var DM, DMR       change in air mass
!@var ENT
!@var F
!@var HEAT1         heating needed for phase change
!@var ML            layer air mass (mb)
!@var SM1, QM1, SMT, QMT dummy variables
!@var SMOLD,QMOLD   profiles prior to any moist convection
!@var SMDNL
!@var TPSAV         array to save plume temperature (set once) (K)
!@var TAUMC1
!@var WCU,  WCU2
!        (KMAX)
!@var SUMU,SUMV,SUMU1,SUMV1
!@var UMP, VMP      momentum carried by convective plumes
!@var UMDN,VMDN     dummy variables
!        (KMAX,LM)
!@var DUM, DVM      changes of UM,VM
!@var UMDNL,VMDN
!        (NMOM)
!@var SMOMP,SMOMPMAX,SMOMDN
!@var QMOMP,QMOMPMAX,QMOMDN
!@var
!        (NMOM,LM)
!@var DQMOM, DQMOMR, DSMOM, DSMOMR
!@var FMOM
!@var QMOMDNL, QMOMOLD, SMOMDNL, SMOMOLD
!
C
C          *********   SCALAR DECLARATIONS   *********
C
      LOGICAL  BELOW_CLOUD, MC1
      INTEGER, INTENT(IN)  :: I_DEBUG, J_DEBUG
      INTEGER, INTENT(OUT) :: IERR,LERR
      INTEGER  IC, ITER, ITYPE, IERRT, K, KSUB,
     *     L, LLMIN, LFRZ,  LERRT, LDRAFT, LMIN, LMAX, LDMIN, LM1,
     *     MCCONT, MAXLVL, MINLVL, N, NPPL, NSUB

      REAL*8
     *     ALPHA,ALPHAU, BETA,BETAU,BYKSUB,BYPBLM,
     *     CDHM,CDHSUM,CLDM,CLDREF,CDHSUM1,CONTCE,CDHDRT,CONDMU,
C
     *     DCW,DCG,DCI,DWCU,DMSE,DFP,DQSUM,DMMIX,DQ,DMSE1,DQSUM1,
     *     DDRAFT,DELTA,DQEVP,DDRUP,DDROLD,
     *     EPLUME,ETADN,ETAL1,EVPSUM,EDRAFT,
C
     *     FLAMW,FLAMG,FLAMI,FG,FI,FMC1,FCLW,
     *     FMP0,FPLUME,FMP2,FRAT1,FRAT2,FCTYPE,FQCOND1,FPOLD,
     *     FDDRT,FCDH,FCDH1,FCLD,FCLOUD,FDDL,FDDP,FENTR,FENTRA,
     *     FEVAP,FLEFT,FQCOND,FQCONDV,FQEVP,FPRCP,FSEVP,FSSUM,
     *     FCONV_tmp,FSUB_tmp,FSSL_tmp,
C
     *     GAMA, HDEP,HPBL, LHX,LHX1,
     *     MPLUME,MPLUM1,MCLOUD,MPMAX,MPOLD,
     *     MNdO,MNdL,MNdI,MCDNCW,MCDNCI,      ! Menon
     *     PBLM,PRCP,PGRAD,
C
     *     QENV,QMO1,QMO2,QDN,QUP,QEDGE,QMN1,QMN2,QMP,QMDN,QMIX,
     *     QMPMAX,QMPT,QSATC,QSATMP, RCLD,RCLDE,RHO,
C
     *     SENV,SMO1,SMO2,SDN,SUP,SEDGE,SVDN,SVUP,SVEDG,SMN1,SMN2,
     *     SMP,SLH,SMDN,SMIX,SMPMAX,SMPT,SUMAJ,SVMIX,SVM1,SUMDP,
C
     *     TADJ,TEMWM,TEM,TIG,TNX1,TP,TVP,TOLD,TOLD1,TTURB,TRATIO,
     *     UMTEMP,VMTEMP,VT,
     *     WMDN,WMUP,WMEDG,WMIX,W2TEM,WTEM,WCONST,WCUFRZ,WORK,WMAX,WV
!
C
C            *********   SCALAR DEFINITIONS   *********
C
!
!@var     Logical Variables
!@var BELOW_CLOUD      is the current level below cloud?
!@var MC1              true for the first convective event
!
!@var     Integer Variables
!@var IC               integer for cloud types
!@var ITER             number for iteration
!@var ITYPE            convective cloud types
!@var IERR,LERR        error reports from advection
!@var IERRT,LERRT      error reports from advection
!@var K,L,N            loop variables
!@var KSUB, BYKSUB     number, 1/number of subsidence iterations
!@var LDRAFT           the layer at which the downdraft orginates
!@var LMIN, LMAX       the base, top layers of a convective event
!@var LDMIN            the lowest layer to which the downdraft descends
!@var LFRZ             freezing level
!@var LM1
!@var MCCONT           integer to count convective events
!@var MAXLVL, MINLVL   the lowest, the highest layer of convective events
!@var NPPL             iteration counter for computing MC area partition
!@var NSUB             LMAX - LMIN + 1
!
!@var     REAL*8  Dummy Variables
!@var       ALPHA,ALPHAU, BETA,BETAU,BYPBLM, CDHDRT,CDHM,CDHSUM,CLDREF,
!@var       DQSUM,DQ,DDRUP,DDROLD, EVPSUM,
!@var       FCDH,FCDH1,FCLD,FCLOUD,FDDL,FDDP,FMP2,FRAT1,FRAT2,  GAMA,
!@var       QMN1,QMN2,QMDN,QMIX,QMPT,QENV,QMO1,QMO2,QDN,QUP,QEDGE,QMPT,QNX,
!@var       SMN1,SMN2,SMDN,SMIX,SMPT,SENV,SMO1,SMO2,SDN,SUP,SEDGE,SMPT,
!@var       SUMAJ,SUMDP,SVDN,SVUP,SVEDG,
!@var       TEMWM,TEM,TIG,TNX,TTURB,TRATIO,  WTEM,WCONST,WMDN,WMUP,WMEDG
!
!@var     REAL*8 Variables
!@var CONDMU        convective condensate in Kg/m^3
!@var CLDM          subsidence due to convection (mb)
!@var CONTCE        scaling factor for entrainment strength
!@var DCW,DCG,DCI   critical cloud particle sizes for onset of precip
!@var DMSE, DMSE1   difference in moist static energy
!@var DFP           an iterative increment
!@var DMMIX
!@var DDRAFT        downdraft mass (mb)
!@var DELTA         fraction of plume that stays in the layer
!@var DQEVP         amount of condensate that evaporates in downdrafts
!@var EPLUME        mass of entrained air (mb)
!@var ETADN         initial downdraft mass / updraft mass
!@var ETAL1         fractional entrainment rate
!@var EDRAFT        mass of air entrained into downdrafts (mb)
!@var FLAMW,FLAMG,FLAMI Marshall-Palmer lambda for water, graupel and ice
!@var FG, FI        fractions of glaciated mass in graupel and ice
!@var FMC1          fraction of grdibox occupied by moist convection + subsidence
!@var FCLW          fraction of condensate in plume that remains as CLW
!@var FMP0          less entraining convective mass (mb)
!@var FPLUME        convective plume mass / layer mass
!@var FCTYPE        fraction of total convective mass for each convective cloud type
!@var FDDRT         fraction of precipitating condensate avaliable for re-evaporation
!@var FENTR         entrained mass / convective mass
!@var FENTRA        entrained mass / layer mass
!@var FEVAP         fraction of layer mass available for precip evaporation
!@var FLEFT         fraction of plume after removing downdraft mass
!@var FQCOND        fraction of water vapor that condenses in plume
!@var FQCONDV       fraction of condensate that is lofted
!@var FQEVP         fraction of water vapor that evaporates in downdraft
!@var FPRCP         fraction of evaporated precipitation
!@var FSEVP         fraction of energy lost to evaporate
!@var FSSUM         fraction of energy lost to evaporate
!@var HDEP
!@var HPBL, PBLM    PBL height (m) and air mass in PBL (mb)
!@var LHX           latent heat of evaporation or sublimation (J/Kg)
!@var MPLUME,MPLUM1 mass of convective plume (mb)
!@var MCLOUD        air mass available for re-evaporation of precip (mb)
!@var MPMAX         mass of convective plume at the detrainment level (mb)
!@var MPOLD
!@var MNdO,MNdL,MNdI,MCDNCW,MCDNCI      Menon Stuff
!@var QMPMAX,SMPMAX values of QMP, SMP in detrained air
!@var QSATC         saturation vapor mixing ratio
!@var QSATMP        plume's saturation vapor mixing ratio
!@var RCLD,RCLDE    cloud particle radius, effective radius (microns)
!@var RHO           air density
!@var SMP, QMP      plume's SM, QM
!@var SLH           LHX/SHA
!@var TADJ          adjustment time for stablization of cloud base by cumulus mass flux
!@var TP            plume's temperature (K)
!@var TVP
!@var TOLD,TOLD1    old temperatures
!@var VT            precip terminal velocity (m/s)
!@var WMIX
!@var W2TEM
!@var WORK          work done on convective plume
!@var WMAX specified maximum convective updraft speed (m/s)
!@var WV convective updraft speed (m/s)
!
C          *******    VARIABLES DEFINED FOR     *******
C          *******  CONDITIONAL COMPILIATIONS   *******
C
#ifdef TRACERS_ON
!@var TMOLD: old TM (tracer mass)
      REAL*8, DIMENSION(LM,NTM)      :: TMOLD, TM1
      REAL*8, DIMENSION(NMOM,LM,NTM) :: TMOMOLD
      REAL*8, DIMENSION(NTM) :: TMP, TMPMAX, TENV, TMDN, TM_dum, DTR
      REAL*8, DIMENSION(NMOM,NTM) :: TMOMP, TMOMPMAX, TMOMDN
#ifdef TRACERS_WATER
!@var TRPCRP tracer mass in precip
      REAL*8, DIMENSION(NTM)      :: TRPRCP
!@var FQCONDT fraction of tracer that condenses
!@var FQEVPT  fraction of tracer that evaporates (in downdrafts)
!@var FPRCPT fraction of tracer that evaporates (in net re-evaporation)
!@var FWASHT  fraction of tracer scavenged by below-cloud precipitation
      REAL*8 :: FQCONDT(NTM), FWASHT(NTM), FPRCPT(NTM), FQEVPT(NTM)
!@var WMXTR available water mixing ratio for tracer condensation ( )?
!@var b_beta_DT precipitating gridbox fraction from lowest precipitating
!@+   layer. The name was chosen to correspond to Koch et al. p. 23,802.
!@var precip_mm precipitation (mm) from the grid box above for washout
      REAL*8 WMXTR, b_beta_DT, precip_mm
c for tracers in general, added by Koch
      REAL*8, DIMENSION(NTM) :: THLAW,THWASH,TR_LEF,TMFAC,TR_LEFT
      REAL*8 CLDSAVT
      INTEGER :: IGAS
!@var TR_LEF limits precurser dissolution following sulfate formation
!@var THLAW Henry's Law determination of amount of tracer dissolution
!@var TMFAC used to adjust tracer moments
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      REAL*8 TMP_SUL(LM,NTM)
c for sulfur chemistry
!@var WA_VOL Cloud water volume (L). Used by GET_SULFATE.
      REAL*8 WA_VOL
      REAL*8, DIMENSION(NTM) ::SULFOUT,SULFIN,SULFINC
      INTEGER :: IAQCH
#endif
      REAL*8 HEFF
#endif
#endif
!
#ifdef CLD_AER_CDNC
       REAL*8, DIMENSION(LM) ::  CONDPC
#endif
!
#ifdef CLD_AER_CDNC
       REAL*8
     *      MCDNO1,MCDNL1,CDNCB,fcnv,ATEMP,VVEL
     *     ,DSGL(LM,SNTM),DSS(SNTM)
     *     ,Repsis,Repsi,Rbeta,RCLD_C
     *     ,AIRM_CDNC,TM_CDNC(NTM)
     *     ,MNdO_max(LM),MNdL_max(LM)
     *     ,MNdO_min(LM),MNdL_min(LM)
#ifdef TRACERS_AMP
      real*8                    :: ncaero (nmodes)
      integer                   ::nm
#endif
#endif
#ifdef CLD_AER_CDNC
!@var MCDNCW,MCDNCI cloud droplet # for warm,cold moist conv clouds (cm^-3)
#endif
#ifdef CLD_AER_CDNC
C     REAL*8 RHO   ! air density
!CN0 is the No parameter in the Marshall-Palmer distribution
#endif
C
C          *******                                          *******
C          *******         END DECLARATION SECTION          *******
C          *******                                          *******
C
!
C****
C**** MOIST CONVECTION
C****
C**** CONVECTION USES A MASS FLUX APPROACH WITH THE MASS PARTITIONED
C**** INTO TWO PLUMES WITH DIFFERENT ENTRAINMENT RATES.  CONVECTION
C**** IS INITIATED AT THE LOWEST MOIST CONVECTIVELY UNSTABLE LEVEL
C**** AND TERMINATES AT THE LEVEL AT WHICH THE UPDRAFT SPEED IS NO
C**** LONGER POSITIVE.  FURTHER CONVECTION MAY THEN ARISE FROM
C**** SUCCESSIVELY HIGHER BASE LEVELS.
C****
C**** THE PARAMETERIZATION CONTAINS 7 PRIMARY PHYSICS LOOPS:
C**** (1) MAJOR OUTER LOOP OVER SUCCESSIVE POTENTIAL CLOUD BASE LAYERS IN
C****     WHICH THE CUMULUS MASS FLUX AND ITS PARTITIONING ARE CALCULATED.
C**** (2) MAJOR INNER LOOP OVER THE TWO CLOUD TYPES CONTAINING MOST OF THE
C****     OTHER PHYSICS LOOPS.
C**** (3) PARCEL ASCENT LOOP IN WHICH CONVECTION IS TRIGGERED, WATER IS
C****     CONDENSED, AIR IS ENTRAINED AND DETRAINED, DOWNDRAFT INITIATION
C****     LEVELS ARE DEFINED, UPDRAFT SPEED IS CALCULATED, CONDENSATE IS
C****     PARTITIONED INTO CLOUD AND PRECIPITATION, AND ASCENT IS TERMINATED.
C**** (4) DOWNDRAFT DESCENT LOOP.
C**** (5) COMPENSATING ENVIRONMENTAL SUBSIDENCE LOOP.
C**** (6) PRECIPITATION AND RE-EVAPORATION LOOP, INCLUDING CONVECTIVE CLOUD
C****     FRACTION CALCULATION.
C**** (7) CONVECTIVE CLOUD OPTICAL THICKNESS LOOP, THE ONLY ONE OUTSIDE
C****     THE MAJOR OUTER LOOP.
C****
      ierr=0 ; lerr=0
      LMCMIN=0
      LMCMAX=0
      MCCONT=0
      FMC1=0.
      FSSL=1.
      RCLDX=radius_multiplier
C**** SET USER-CONTROLLED FREE PARAMETERS
      TADJ = 1.d0
C     PGRAD = 0.7d0
      PGRAD = 0.7                  ! to maintain bit compatibility
      CONTCE1=entrainment_cont1
      CONTCE2=entrainment_cont2
      FDDRT = 1.d0
C**** initiallise arrays of computed output
      TAUMCL=0
      SVWMXL=0
      SVLATL=0
      SVLAT1=0
      VSUBL=0
      PRECNVL=0
      CLDMCL=0
      CLDSLWIJ=0
      CLDDEPIJ=0
      PRCPMC=0.
      TPSAV=0
      CSIZEL=RWCLDOX*10.*(1.-PEARTH)+10.*PEARTH ! droplet rad in stem
      VLAT=LHE
      SAVWL=0.
      SAVWL1=0.
      LHP=0
#ifdef SCM
      if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
        LHPMC=0. ; PREMC = 0.
        WCUDEEP=0. ; MPLUMEDEEP=0. ; ENTDEEP=0.
        DETRAINDEEP = 0. ; PRCCDEEP = 0. ; NPRCCDEEP = 0.
        TPALL = 0. ; MCCOND = 0. ; PRCCGRP = 0. ; PRCCICE = 0.
        PLUME_MAX = 0. ; PLUME_MIN = 0.
      endif
#endif
      SAVE1L=0.
      SAVE2L=0.
#ifdef TRACERS_WATER
      trsvwml = 0.
      TRPRCP = 0.
      TRPRMC = 0.
#ifdef TRDIAG_WETDEPO
      IF (diag_wetdep == 1) THEN
c**** initialize diagnostic arrays
        trcond_mc=0.D0
        trdvap_mc=0.D0
        trflcw_mc=0.D0
        trprcp_mc=0.D0
        trnvap_mc=0.D0
        trwash_mc=0.D0
      END IF
#endif
#endif
#ifdef CFMIP
        cfmip_ccp(:) = 0.
        cfmip_ccl(:) = 0.
#endif
C**** zero out diagnostics
         MCFLX =0.
         DGDSM=0.
         DGDEEP=0.
         DGSHLW=0.
         DPHASE=0.
         DPHADEEP=0.
         DPHASHLW=0.
         DTOTW=0.
         DQCOND=0.
         DGDQM=0.
         DDMFLX=0.
         TDNL=0.
         QDNL=0.
#ifdef SCM
         if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
             CUMFLX=0.
             DWNFLX=0.
         endif
#endif
C**** save initial values (which will be updated after subsid)
      SM1=SM
      QM1=QM
#ifdef TRACERS_ON
      TM1(:,1:NTX) = TM(:,1:NTX)
      TRDNL = 0.
#ifdef TRACERS_WATER
      CLDSAVT=0.
#endif
#endif
C**** SAVE ORIG PROFILES
      SMOLD(:) = SM(:)
      SMOMOLD(:,:) = SMOM(:,:)
      QMOLD(:) = QM(:)
      QMOMOLD(:,:) = QMOM(:,:)
#ifdef TRACERS_ON
      TMOLD(:,1:NTX) = TM(:,1:NTX)
      TMOMOLD(:,:,1:NTX) = TMOM(:,:,1:NTX)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      DT_SULF_MC(1:NTM,:)=0.
#endif
#ifdef TRACERS_WATER
! TR_LEF is an input to get_cond_factor not currently used for MC clouds
      TR_LEF(:)=1.D0
      thlaw(:) = 0.  ! nonzero only for gas tracers
      thwash(:) = 0. ! nonzero only for gas tracers
      tmfac(:) = 0.  ! nonzero only for gas tracers
      fwasht(:) = 0. ! nonzero only for aerosols
      fqcondt(:) = 0.
#endif
#endif
C**** CALULATE PBL HEIGHT AND MASS
      PBLM=0.
      HPBL=0.
      DO L=1,DCL
        IF(L.LT.DCL) THEN
          PBLM=PBLM+AIRM(L)
          HPBL=HPBL+AIRM(L)*TL(L)*RGAS/(GRAV*PL(L))
        ELSE
          PBLM=PBLM+.5d0*AIRM(L)
          HPBL=HPBL+.5d0*AIRM(L)*TL(L)*RGAS/(GRAV*PL(L))
        END IF
      END DO
      BYPBLM=1.d0/PBLM
C**** CALCULATE THRESHOLD RH FOR LARGE-SCALE CLOUD FORMATION IN PBL
C**** BASED ON SIEBESMA ET AL. (2003, JAS)
      DO L=1,LM
        U00L(L)=0.
        IF(PL(L).GE.850.d0) THEN
          LHX=LHE
C         IF(TL(L).LT.TF) LHX=LHS ! use 10C
          U00L(L)=1.d0-2.*(U00b*.001*.050*(HPBL/500.)*
C    *         (.001*SQRT(DXYPIJ))**.33)/QSAT(TL(L),LHX,PL(L))
     *         (.001*SQRT(DXYPIJ))**.33)/QSAT(283.16d0,LHX,PL(L))
        END IF
      END DO
C**** CALCULATE DEL WCU TO TRAVEL HALF LAYER THICKNESS IN ONE
C**** TIMESTEP; USED LATER TO DETERMINE FRACTION OF CLOUD WATER
C**** THAT DETRAINS AT CURRENT LAYER
      DWCU=0.
      DO L=1,LMCM
        DWCU=DWCU+AIRM(L)*TL(L)*RGAS/(GRAV*PL(L))
      END DO
      DWCU=0.5*DWCU*BYDTsrc/real(LMCM)

#ifdef CLD_AER_CDNC
      MNdO_max(:)=teeny
      MNdL_max(:)=teeny
      MNdO_min(:)=teeny
      MNdL_min(:)=teeny
#endif

c      save Moist static energy profile
#ifdef SCM
      if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
         do L=1,LM
            SCM_H(L) = sha*TL(L) + SCM_GZ(L) + LHE*QL(L)
            SCM_H(L) = SCM_H(L)/1000.
            ARM_H(L) = sha*SG_T(L) + SCM_GZ(L) + LHE*SG_Q(L)
            ARM_H(L) = ARM_H(L)/1000.
            scmqsat = QSAT(TL(L),LHE,PL(L))
            armqsat = QSAT(SG_T(L),LHE,PL(L))
            SCM_HSAT(L) = sha*TL(L) + SCM_GZ(L)
     &                  + LHE*scmqsat
            SCM_HSAT(L) = SCM_HSAT(L)/1000.
            ARM_HSAT(L) = sha*SG_T(L) + SCM_GZ(L) + LHE*armqsat
            ARM_HSAT(L) = ARM_HSAT(L)/1000.
c           write(iu_scm_prt,7878) NSTEPSCM,L,sha,TL(L),SCM_GZ(L),
c    &        QL(L),scmqsat,SCM_H(L),SCM_HSAT(L),ARM_H(L),ARM_HSAT(L)
c7878       format(1x,'nstep L cp T GZ Q Qsat H Hsat ARM H Hsat',2(i5),
c    &             f8.2,f8.2,f10.2,f8.5,f8.5,4(f8.2))
         enddo
      endif
#endif

C****
C**** BEGIN OUTER LOOP (1) OVER BASE LAYERS
C****

      CLOUD_BASE: DO LMIN=1,LMCM-1
      MAXLVL=0
      MINLVL=LM

c****
C**** THE MASS FLUX CLOSURE CALCULATES THE MASS REQUIRED TO RESTORE THE
c**** CLOUD BASE LEVEL TO NEUTRAL BUOYANCY AFTER SUBSIDENCE (YAO AND
c**** DEL GENIO 1989, J. CLIM.).  THE RESULTING CUMULUS MASS FLUX IS
c**** ASSUMED TO OCCUR OVER A CONVECTIVE ADJUSTMENT TIME OF 1 HOUR
c****

C**** COMPUTE THE CONVECTIVE MASS OF THE LESS ENTRAINING PART BASED ON THE
c**** LARGE-SCALE VERTICAL VELOCITY AT CLOUD BASE
      FMP0=-10.*CK1*SDL(LMIN+1)*BYGRAV*XMASS
      IF(FMP0.LE.0.) FMP0=0.

C**** CREATE A PLUME IN THE BOTTOM LAYER
C****
C**** DETERMINE IF CLOUD BASE LEVEL IS UNSTABLE
C**** ITERATION TO FIND FPLUME WHICH RESTORES CLOUD BASE TO NEUTRAL STATE
C****
      SMO1=SM(LMIN)
      QMO1=QM(LMIN)
      SMO2=SM(LMIN+1)
      QMO2=QM(LMIN+1)
      SDN=SMO1*BYAM(LMIN)
      SUP=SMO2*BYAM(LMIN+1)
      SEDGE=THBAR(SUP,SDN)
      QDN=QMO1*BYAM(LMIN)
      QUP=QMO2*BYAM(LMIN+1)
      WMDN=WML(LMIN)
      WMUP=WML(LMIN+1)
      SVDN=SDN*(1.+DELTX*QDN-WMDN)
      SVUP=SUP*(1.+DELTX*QUP-WMUP)
      QEDGE=.5*(QUP+QDN)
      WMEDG=.5*(WMUP+WMDN)
      SVEDG=SEDGE*(1.+DELTX*QEDGE-WMEDG)
      LHX=LHE
      SLH=LHX*BYSHA
      DMSE=(SVUP-SVEDG)*PLK(LMIN+1)+(SVEDG-SVDN)*PLK(LMIN)+
     *  SLHE*(QSAT(SUP*PLK(LMIN+1),LHX,PL(LMIN+1))-QDN)
      IF(DMSE.GT.-1d-10) CYCLE  ! try next level

C**** MASS_FLUX PERFORMS THE ITERATIONS
C**** COMPUTES FPLUME, FMP2
      CALL MASS_FLUX (FPLUME, FMP2, DQSUM, LMIN,
     *  LHX, QMO1, QMO2, SLH, SMO1, SMO2, WMDN, WMUP, WMEDG)

      IF(FPLUME.LE..001) CYCLE ! try next level

C****
C**** BEGIN LOOP (2) THROUGH CLOUD TYPES OF DIFFERENT ENTRAINMENT RATE
C****

      ITYPE=2           ! always 2 types of clouds: less and more entraining
      FCTYPE=1.

c**** STABILIZATION IS ASSUMED TO OCCUROVER 1 HOUR, SO ONLY APPLY A FRACTION
c**** OF THE REQUIRED MASS FLUX IN ONE PHYSICS TIMESTEP
      FMP2=FMP2*min(1d0,DTsrc/(TADJ*3600.d0))
      WMAX=50.
C****
      CLOUD_TYPES:  DO IC=1,ITYPE
C****
C**** Initialise plume characteristics
      MC1=.FALSE.    ! flag for first convection event
      LHX=LHE
      MPLUME=MIN(AIRM(LMIN),AIRM(LMIN+1))
      IF(MPLUME.GT.FMP2) MPLUME=FMP2

      IF(ITYPE.EQ.2) THEN     ! cal. MPLUME for 1st plume and 2nd plume
        FCTYPE=1.
        IF(MPLUME.GT.FMP0) FCTYPE=FMP0/MPLUME
        IF(IC.EQ.2) FCTYPE=1.-FCTYPE
        IF(FCTYPE.LT.0.001) CYCLE CLOUD_TYPES
      END IF
      MPLUM1=MPLUME

C****
C**** Calculate gridbox area partition between (moist convection + subsidence)
C**** and (large-scale cloud + clear sky).  A first calculation of the subsidence
C**** area is performed.  If the convection rises 2 layers or more (MCCONT > 1)
C**** subsidence is assumed to occur over a greater area relative to the updraft
C****
C**** To test code with only one estimate of area partition
C**** remove NPPL Loop, or set NPPL=1,1
C****
      AREA_PARTITION: DO NPPL=1,2

        IF(NPPL.EQ.2)   THEN                  ! Perform 2nd pass if MC
          IF(.NOT.MC1 .OR. MCCONT.LT.2)  THEN ! went more than 2 layers
            CYCLE AREA_PARTITION              ! else - skip 2nd pass
          ELSE
#ifdef TRACERS_ON
          call reset_tracer_work_arrays(lmin,lmax)
#endif
          END IF
        END IF

      MPLUME=MPLUM1*FCTYPE

C**** Calculate fraction for plume either for the first time, or if
C**** we are redoing the calculation for deeper convection
      IF (MCCONT.eq.0  .OR. MC1) THEN

        IF (MCCONT.eq.0) THEN   ! first time
          FSUB_tmp=1.d0+(AIRM(LMIN+1)-100.d0)/200.d0
        ELSE                    ! redoing plume calc
          FSUB_tmp=1.d0+(PL(LMIN)-PL(LMAX)-100.d0)/200.d0
        END IF

        FCONV_tmp=MIN(MPLUM1*BYAM(LMIN+1),1d0)
        IF(FSUB_tmp.GT.1.d0/(FCONV_tmp+1.d-20)-1.d0)
     *       FSUB_tmp=1.d0/(FCONV_tmp+1.d-20)-1.d0
        FSUB_tmp=MAX(1.d0,MIN(FSUB_tmp,5.d0))

        FSSL_tmp=1.d0-(1.d0+FSUB_tmp)*FCONV_tmp
        FSSL_tmp=MAX(CLDMIN,MIN(FSSL_tmp,1d0-CLDMIN))

        FMC1=(1.d0-FSSL_tmp)+teeny
      END IF

C**** guard against possibility of too big a plume
      IF (MC1 .or. MCCONT.GT.0) THEN
        MPLUME=MIN(0.95d0*AIRM(LMIN)*FMC1,MPLUME)
      END IF

      DO L=1,LM
        COND(L)=0.    ;    CDHEAT(L)=0. !;  VLAT(L)=LHE
        CONDP(L)=0.   ;  CONDP1(L)=0. ;    CONDGP(L)=0.
        CONDIP(L)=0.  ;  CONDV(L)=0.  ;    HEAT1(L)=0.
        DM(L)=0.      ;  DMR(L)=0.    ;    DDR(L)=0.
        CCM(L)=0.     ;  DDM(L)=0.    ;    TAUMC1(L)=0.
        ENT(L)=0.     ;  DET(L)=0.    ;    BUOY(L)=0.
        WCU(L)=0.     ; SMDNL(L)=0.   ;  QMDNL(L)=0.
#ifdef SCM
        if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
          WCUALL(L,IC,LMIN) = 0.
          MPLUMEALL(L,IC,LMIN) = 0.
          ENTALL(L,IC,LMIN) = 0.
        endif
#endif
      END DO
      SMOMDNL(:,:)=0.   ;  QMOMDNL(:,:)=0.
      UMDN(1:KMAX)=0.   ;  VMDN(1:KMAX)=0.
      UMDNL(1:KMAX,:)=0.  ;  VMDNL(1:KMAX,:)=0.
      DUM(1:KMAX,:)=0.    ;  DVM(1:KMAX,:)=0.
      DSM(:) = 0. ; DSMOM(:,:) = 0. ; DSMR(:) = 0. ; DSMOMR(:,:) = 0.
      DQM(:) = 0. ; DQMOM(:,:) = 0. ; DQMR(:) = 0. ; DQMOMR(:,:) = 0.
c#ifdef TRACERS_ON
c (re)zeroing now done post-calculation with calls to
c reset_tracer_work_arrays
c      DTM(:,1:NTX) = 0.   ;   DTMOM(:,:,1:NTX) = 0.
c      DTMR(:,1:NTX) = 0.  ;  DTMOMR(:,:,1:NTX) = 0.
c      TMDNL(:,1:NTX) = 0.  ;  TMOMDNL(:,:,1:NTX) = 0.
c      TPOLD = 0.
c#endif
c#ifdef TRACERS_WATER
c      TRCOND = 0.  ; TRCONDV = 0.
c#endif

C**** adjust MPLUME to take account of restricted area of subsidence
C**** (i.e. MPLUME is now a greater fraction of the relevant airmass.
      MPLUME=MIN( MPLUME/FMC1,
     *   AIRM(LMIN)*0.95d0*QM(LMIN)/(QMOLD(LMIN) + teeny) )
      IF(MPLUME.LE..001*AIRM(LMIN)) CYCLE CLOUD_TYPES
      FPLUME=MPLUME*BYAM(LMIN)
      SMP  =  SMOLD(LMIN)*FPLUME
      SMOMP(xymoms)=SMOMOLD(xymoms,LMIN)*FPLUME
      QMP  =  QMOLD(LMIN)*FPLUME
      QMOMP(xymoms)=QMOMOLD(xymoms,LMIN)*FPLUME
      IF (TPSAV(LMIN).eq.0) TPSAV(LMIN)=SMP*PLK(LMIN)/MPLUME
      DMR(LMIN)=-MPLUME
        DSMR(LMIN)=-SMP
      DSMOMR(xymoms,LMIN)=-SMOMP(xymoms)
      DSMOMR(zmoms,LMIN)=-SMOMOLD(zmoms,LMIN)*FPLUME
        DQMR(LMIN)=-QMP
      DQMOMR(xymoms,LMIN)=-QMOMP(xymoms)
      DQMOMR(zmoms,LMIN)=-QMOMOLD(zmoms,LMIN)*FPLUME
#ifdef TRACERS_ON
C**** This is a fix to prevent very occasional plumes that take out
C**** too much tracer mass. This can impact tracers with very sharp
C**** vertical gradients
      do n=1,ntx
        TMP(n) = TMOLD(LMIN,n)*FPLUME
        if(t_qlimit(n)) TMP(n) = MIN(TMP(n),0.95d0*TM(LMIN,n))
      enddo
      TMOMP(xymoms,1:NTX)=TMOMOLD(xymoms,LMIN,1:NTX)*FPLUME
        DTMR(LMIN,1:NTX)=-TMP(1:NTX)
      DTMOMR(xymoms,LMIN,1:NTX)=-TMOMP(xymoms,1:NTX)
      DTMOMR( zmoms,LMIN,1:NTX)=-TMOMOLD(zmoms,LMIN,1:NTX)*FPLUME
      TPOLD(LMIN)=TPSAV(LMIN)  ! initial plume temperature
#endif
      DO K=1,KMAX
         UMP(K)=UM(K,LMIN)*FPLUME
         DUM(K,LMIN)=-UMP(K)
         VMP(K)=VM(K,LMIN)*FPLUME
         DVM(K,LMIN)=-VMP(K)
      END DO
C****
C**** RAISE THE PLUME TO THE TOP OF CONVECTION AND CALCULATE
C**** ENTRAINMENT, CONDENSATION, AND SECONDARY MIXING
C****
      CDHSUM=0.
      CDHSUM1=0.
      CDHDRT=0.
      ETADN=0.
      LDRAFT=LM
      EVPSUM=0.
      DDRAFT=0.
      LFRZ=0
      LMAX=LMIN
      CONTCE=CONTCE1
      IF(IC.EQ.2) CONTCE=CONTCE2
C**** INITIAL CUMULUS UPDRAFT SPEED DETERMINED FROM SQRT(TKE), BUT
C**** THE LARGER OF 0.5 M/S OR THE TURBULENT VERTICAL VELOCITY FOR
C**** THE MORE ENTRAINING PLUME; THE LARGER OF 0.5 M/S OR TWICE THE
C**** TURBULENT VERTICAL VELOCITY FOR THE LESS ENTRAINING PLUME
      WCU(LMIN)=MAX(.5D0,WTURB(LMIN))
      IF(IC.EQ.1) WCU(LMIN)=MAX(.5D0,2.D0*WTURB(LMIN))
      WCU2(LMIN)=WCU(LMIN)*WCU(LMIN)

C****
C**** BEGIN LOOP (3) OVER POSSIBLE CLOUD TOP LEVELS
C****

      CLOUD_TOP:  DO L=LMIN+1,LM

C****
C**** TRIGGERING CONDITIONS FOR MOIST CONVECTION
C****

C**** (1) TEST WHETHER MASS OF AIR REQUIRED TO STABILIZE CLOUD BASE
C**** LARGE ENOUGH TO WARRANT PERFORMING CALCULATIONS
      IF(MPLUME.LE..001*AIRM(L)) EXIT CLOUD_TOP

C**** (2) TEST WHETHER VIRTUAL MOIST STATIC ENERGY OF PARCEL LIFTED
C**** TO NEXT LEVEL EXCEEDS THAT OF ENVIRONMENT AT THAT LEVEL, I.E.,
C**** WHETHER PARCEL IS BUOYANT
      SDN=SMP/MPLUME
      SUP=SM1(L)*BYAM(L)
      QDN=QMP/MPLUME
      QUP=QM1(L)*BYAM(L)
      WMDN=0.
      WMUP=WML(L)
      SVDN=SDN*(1.+DELTX*QDN-WMDN)
      SVUP=SUP*(1.+DELTX*QUP-WMUP)
      IF(PLK(L-1)*(SVUP-SVDN)+SLHE*(QUP-QDN).GE.0.) EXIT CLOUD_TOP
C**** Only set TPSAV if it hasn't already been set
      IF (TPSAV(L).eq.0) TPSAV(L)=SMP*PLK(L)/MPLUME
      TP=TPSAV(L)
      IF(TPSAV(L-1).GE.TF.AND.TPSAV(L).LT.TF) LFRZ=L-1

C**** (3)TEST TO SEE WHETHER LIFTED PARCEL IS ABOVE LIFTING CONDENSATION
C**** W.R.T. LIQUID WATER, OR W.R.T. ICE FOR HOMOGENEOUS NUCLEATION (T<-40)
      LHX=LHE
      IF(TP.LT.TI) LHX=LHS
      QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
      IF(QMP.LT.QSATMP) EXIT CLOUD_TOP    !  no plume

      IF(TP.LT.TF.AND.LHX.EQ.LHE) THEN
        LHX=LHS
        QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
      END IF

C**** DEFINE DUMMY LATENT HEAT VARIABLE TO AVOID PHASE DISCREPANCY BETWEEN PLUMES
      IF (VLAT(L).EQ.LHS) LHX=LHS
      VLAT(L)=LHX
      SLH=LHX*BYSHA

C**** THRESHOLD RH FOR PBL STRATIFORM CLOUDS
      IF(TL(L).GE.TF .AND. U00L(L).NE.U00a) U00L(L)=
     *    1.d0-2.*(U00b*2.d-4*MPLUME/QSATMP)

      MCCONT=MCCONT+1
      IF(MCCONT.EQ.1) MC1=.TRUE.
C****
C**** IF PLUME MASS IS TOO LARGE FOR UPPER LAYER, LEAVE PART BEHIND IN LOWER LAYER
C**** AND ADJUST TEMPERATURE, HUMIDITY, AND MOMENTUM THERE
C****
      IF(MPLUME.GT..95*AIRM(L)) THEN
        DELTA=(MPLUME-.95*AIRM(L))/MPLUME
        DM(L-1)=DM(L-1)+DELTA*MPLUME
        MPLUME=.95*AIRM(L)

        DSM(L-1)=  DSM(L-1)+DELTA*SMP
        SMP = SMP  *(1.-DELTA)
        DSMOM(xymoms,L-1)=DSMOM(xymoms,L-1)+DELTA*SMOMP(xymoms)
        SMOMP(xymoms) = SMOMP(xymoms)*(1.-DELTA)

        DQM(L-1)=  DQM(L-1)+DELTA*QMP
        QMP = QMP  *(1.-DELTA)
        DQMOM(xymoms,L-1)=DQMOM(xymoms,L-1)+DELTA*QMOMP(xymoms)
        QMOMP(xymoms) = QMOMP(xymoms)*(1.-DELTA)

        DO K=1,KMAX
          DUM(K,L-1)=DUM(K,L-1)+UMP(K)*DELTA
          DVM(K,L-1)=DVM(K,L-1)+VMP(K)*DELTA
          UMP(K)=UMP(K)-UMP(K)*DELTA
          VMP(K)=VMP(K)-VMP(K)*DELTA
        END DO

#ifdef TRACERS_ON
        DTM(L-1,1:NTX) = DTM(L-1,1:NTX)+DELTA*TMP(1:NTX)
        DTMOM(xymoms,L-1,1:NTX)=DTMOM(xymoms,L-1,1:NTX)+DELTA
     *       *TMOMP(xymoms,1:NTX)
        TMP(1:NTX) = TMP(1:NTX)*(1.-DELTA)
        TMOMP(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)*(1.-DELTA)
#endif
      END IF

C**** WORK DONE BY CONVECTION IN UPPER LAYER REMOVES ENERGY FROM THE PLUME
      WORK=MPLUME*(SUP-SDN)*(PLK(L-1)-PLK(L))/PLK(L-1)
      DSM(L-1)=DSM(L-1)-WORK
      CCM(L-1)=MPLUME

C**** CONDENSE VAPOR IN THE PLUME AND ADD LATENT HEAT
      call get_dq_cond(smp,qmp,plk(l),mplume,lhx,pl(l),dqsum,fqcond)

      IF(DQSUM.GT.0. .AND. QMP.gt.teeny) THEN
        QMOMP(xymoms) =  QMOMP(xymoms)*(1.-FQCOND)
        SMP=SMP+SLH*DQSUM/PLK(L)
        QMP=QMP-DQSUM
      END IF

#ifdef TRACERS_ON
C**** save plume temperature after possible condensation
      TPOLD(L)=SMP*PLK(L)/MPLUME
#endif

C**** TOTAL CONDENSATE IN LAYER = NEWLY CONDENSED + ADVECTED FROM LOWER LAYER;
C**** THE LATTER IS CALCULATED IN THE MICROPHYSICS SECTION
      COND(L)=DQSUM
      CDHEAT(L)=SLH*COND(L)          ! calculate CDHEAT before add CONDV
      CDHSUM=CDHSUM+CDHEAT(L)
      COND(L)=COND(L)+CONDV(L-1)     ! add in the vertical transported COND

      IF (VLAT(L-1).ne.VLAT(L)) THEN
        SMP=SMP-(VLAT(L-1)-VLAT(L))*CONDV(L-1)*BYSHA/PLK(L) ! phase change of uplifted condensate
        CDHEAT(L)=CDHEAT(L)-(VLAT(L-1)-VLAT(L))*CONDV(L-1)*BYSHA
        CDHSUM=CDHSUM-(VLAT(L-1)-VLAT(L))*CONDV(L-1)*BYSHA
      END IF

C**** CALCULATE SLOPE PARAMETER FOR MARSHALL-PALMER SIZE DISTRIBUTION FOR
C**** LIQUID, GRAUPEL, ICE (DEL GENIO ET AL. 2005, J. CLIM.)
      CONDMU=100.*COND(L)*PL(L)/(CCM(L-1)*TL(L)*RGAS)
      FLAMW=(1000.d0*PI*CN0/(CONDMU+teeny))**.25
      FLAMG=(400.d0*PI*CN0G/(CONDMU+teeny))**.25
      FLAMI=(100.d0*PI*CN0I/(CONDMU+teeny))**.25

#if (defined CLD_AER_CDNC) && (defined TRACERS_AEROSOLS_Koch)
!@auth Menon  saving aerosols mass for CDNC prediction
      DO N=1,SNTM
       DSS(N)=1.d-10
       DSGL(L,N)=1.d-10
      ENDDO
#endif

#ifdef CLD_AER_CDNC
#ifdef ALT_CDNC_INPUTS
c aerosols in the updraft
      tm_cdnc(:) = tmp(:)
      airm_cdnc = mplume
#else
c ambient aerosols at this level
      tm_cdnc(:) = tm(l,:)
      airm_cdnc = airm(l)
#endif
#endif
C**** Here we change convective precip due to aerosols
#if (defined CLD_AER_CDNC) && (defined TRACERS_AEROSOLS_Koch)
      DO N=1,NTX
       select case (trname(ntix(n)))
        case('SO4')
          DSGL(L,1)=tm_cdnc(n)     !n=19
          DSS(1) = DSGL(L,1)
        case('seasalt1')
          DSGL(L,2)=tm_cdnc(n)     !n=21
          DSS(2) = DSGL(L,2)
        case('seasalt2')
          DSGL(L,3)=tm_cdnc(n)     !n=22
          DSS(3) = DSGL(L,3)
        case('OCIA')
          DSGL(L,4)=tm_cdnc(n)     !n=27
          DSS(4) = DSGL(L,4)
        case('OCB')
          DSGL(L,5)=tm_cdnc(n)     !n=28
          DSS(5) = DSGL(L,5)
        case('BCIA')
          DSGL(L,6)=tm_cdnc(n)     !n=24
          DSS(6) = DSGL(L,6)
        case('BCB')
          DSGL(L,7)=tm_cdnc(n)     !n=25
          DSS(7) = DSGL(L,7)
        case('OCII')
          DSGL(L,8)=tm_cdnc(n)     !n=26
          DSS(8) = DSGL(L,8)
        case('BCII')
          DSGL(L,9)=tm_cdnc(n)     !n=23
          DSS(9) = DSGL(L,9)
#ifdef TRACERS_DUST
        case('Clay')
          DSGL(L,10)=tm_cdnc(n)    !n=23
          DSS(10) = DSGL(L,10)
        case('Silt1')
          DSGL(L,11)=tm_cdnc(n)    !n=23
          DSS(11) = DSGL(L,11)
        case('Silt2')
          DSGL(L,12)=tm_cdnc(n)    !n=23
          DSS(12) = DSGL(L,12)
        case('Silt3')
          DSGL(L,13)=tm_cdnc(n)    !n=23
          DSS(13) = DSGL(L,13)
#endif
#ifdef TRACERS_NITRATE
        case('NO3p')
          DSGL(L,14)=tm_cdnc(n)    !n=23
          DSS(14) = DSGL(L,14)
#endif
#ifdef TRACERS_HETCHEM
C**** Here are dust particles coated with sulfate
        case('SO4_d1')
          DSGL(L,15)=tm_cdnc(n)    !n=20
          DSS(15) = DSGL(L,15)
        case('SO4_d2')
          DSGL(L,16)=tm_cdnc(n)    !n=21
          DSS(16) = DSGL(L,16)
        case('SO4_d3')
          DSGL(L,17)=tm_cdnc(n)    !n=22
          DSS(17) = DSGL(L,17)
#endif
#ifdef TRACERS_AEROSOLS_SOA
        case('isopp1a')
          DSGL(L,18)=tm_cdnc(n)
          DSS(18) = DSGL(L,18)
        case('isopp2a')
          DSGL(L,19)=tm_cdnc(n)
          DSS(19) = DSGL(L,19)
#ifdef TRACERS_TERP
        case('apinp1a')
          DSGL(L,20)=tm_cdnc(n)
          DSS(20) = DSGL(L,20)
        case('apinp2a')
          DSGL(L,21)=tm_cdnc(n)
          DSS(21) = DSGL(L,21)
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
        end select
      END DO      !end of n loop for tracers
#endif  /* (TRACERS_AEROSOLS_Koch) and (CLD_AER_CDNC) */
C** Use MATRIX AMP_actv to decide what the aerosol number conc. is
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
#ifdef TRACERS_AMP
        do nm=1,nmodes
         ncaero(nm)=naerc(l,nm)*1.d-6
c         if(naerc(l,nm).gt.1.d-30) write(6,*)"mat",ncaero(nm),nm
        enddo
       CALL GET_CC_CDNC_MX(L,nmodes,ncaero,MCDNL1,MCDNO1)
#else
C** This is for the old mass to number calculations nc. is
       CALL GET_CC_CDNC(L,AIRM_CDNC,DXYPIJ,PL(L),TL(L),DSS,
     &      MCDNL1,MCDNO1)

#endif  /* (TRACERS_AMP) */
      MNdO=MCDNO1
      MNdL=MCDNL1
      MNdO_max(L)=max(MNdO_max(L),MCDNO1)
      MNdL_max(L)=max(MNdL_max(L),MCDNL1)
      if(MNdO_min(L)==teeny) then
        MNdO_min(L)=MCDNO1
        MNdL_min(L)=MCDNL1
      else
        MNdO_min(L)=min(MNdO_min(L),MCDNO1)
        MNdL_min(L)=min(MNdL_min(L),MCDNL1)
      endif
      MNdI = 0.06417127d0
      MCDNCW=MNdO*(1.-PEARTH)+MNdL*PEARTH
      MCDNCI=MNdI
c      if(MCDNCW.gt.0.)write(6,*)"Mst CDNC,a",MCDNCW,MNdO,MNdL,L
      if (MCDNCW.le.20.d0) MCDNCW=20.d0     !set min CDNC, sensitivity test
      if (MCDNCW.ge.2000.d0) MCDNCW=2000.d0     !set max CDNC, sensitivity test
C** Using the Liu and Daum paramet, Nature, 2002, Oct 10, Vol 419
C** for spectral dispersion effects on droplet size distribution
C** central value of 0.003 for alfa:Rotstayn&Daum, J.Clim,2003,16,21, Nov 2003.
      Repsi=1.d0 - 0.7d0*exp(-0.003d0*MCDNCW)
      Repsis=Repsi*Repsi
      Rbeta=(((1.d0+2.d0*Repsis)**0.667d0))/((1.d0+Repsis)**by3)
      RCLD_C=14.d0/Rbeta       !set Reff to threshold size =14 um (Rosenfeld)
#endif /* (CLD_AER_CDNC) */
      TAUMC1(L)=TAUMC1(L)+COND(L)*FMC1

#ifdef TRACERS_WATER
C**** CONDENSING TRACERS
      WMXTR=DQSUM*BYAM(L)

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      WA_VOL=COND(L)*1.d2*BYGRAV*DXYPIJ

      IF (FPLUME.GT.teeny) then
        TMP_SUL(L,aqchem_list)=TMP(aqchem_list)/FPLUME
      else
        TMP_SUL(L,aqchem_list)=0.
      END IF

      CALL GET_SULFATE(L,TPOLD(L),FPLUME,WA_VOL,WMXTR,SULFIN,
     *     SULFINC,SULFOUT,TR_LEFT,TMP_SUL,TRCOND(1,L),
     *     AIRM,LHX,DT_SULF_MC(1,L),CLDSAVT)

      do iaqch=1,aqchem_count
        n = aqchem_list(iaqch)
c first apply chemistry
c removal of precursers
        TMP(N)=TMP(N)*(1.+SULFIN(N))
        TMOMP(xymoms,N)= TMOMP(xymoms,N)*(1.+SULFIN(N))
c formation of sulfate
        TRCOND(N,L) = TRCOND(N,L)+SULFOUT(N)
      enddo

#endif

      TM_dum(:) = TMP(:)
c     CALL GET_COND_FACTOR_array(
c    &      NTX,WMXTR,TPOLD(L),TPOLD(L-1),LHX,FPLUME
c    &     ,FQCOND,FQCONDT,.true.,TRCOND(:,L),TM_dum,THLAW,TR_LEF,PL(L)
c    &     ,ntix,CLDSAVT)
      CALL GET_COND_FACTOR_array(
     &      NTX,WMXTR,TPOLD(L),TPOLD(L-1),LHX,FPLUME
     &     ,FQCOND,FQCONDT,.true.,TRCOND(:,L),TM_dum,THLAW,TR_LEF,PL(L)
     &     ,ntix,FPLUME)
      dtr(1:ntx) = fqcondt(1:ntx)*tmp(1:ntx)
#ifdef TRDIAG_WETDEPO
      IF (diag_wetdep == 1) trcond_mc(l,1:ntx)=trcond_mc(l,1:ntx)
     &     +dtr(1:ntx)
#endif
      DO N=1,NTX
        TRCOND(N,L) = DTR(N) + TRCOND(N,L) + TRCONDV(N,L-1)
        TMP(N)         = TMP(N)         *(1.-FQCONDT(N))
        TMOMP(xymoms,N)= TMOMP(xymoms,N)*(1.-FQCONDT(N))
      END DO
#endif

C****
C**** ENTRAINMENT
C****

C**** ENTRAINMENT USES THE PARAMETERIZATION OF GREGORY (2001, QJRMS)
C**** WITH CONSTANT C = 0.3 (0.6) FOR THE LESS (MORE) ENTRAINING PLUME
C****
      TVP=(SMP/MPLUME)*PLK(L)*(1.+DELTX*QMP/MPLUME)
      BUOY(L)=(TVP-TVL(L))/TVL(L)-COND(L)/MPLUME
      ENT(L)=.16667D0*CONTCE*GRAV*BUOY(L)/(WCU(L-1)*WCU(L-1)+teeny)

C**** WHEN PARCEL LOSES BUOYANCY TERMINATE ENTRAINMENT AND DETRAIN AN
C**** AMOUNT OF AIR EQUAL TO THAT WHICH WOULD HAVE BEEN ENTRAINED HAD
C**** THE BUOYANCY BEEN POSITIVE
      IF (ENT(L).LT.0.D0) THEN
        DET(L)=-ENT(L)
        ENT(L)=0.D0
      END IF

      IF(ENT(L).GT.0.D0) THEN    ! non-zero entrainment
      FENTR=1000.D0*ENT(L)*GZL(L)*FPLUME
      IF(FENTR+FPLUME.GT.1.) THEN
        FENTR=1.-FPLUME
        ENT(L)=0.001d0*FENTR/(GZL(L)*FPLUME)
      END IF
      IF(FENTR.GE.teeny) THEN    !  Big Enough, Proceed
      MPOLD=MPLUME
      FPOLD=FPLUME
      ETAL1=FENTR/(FPLUME+teeny)
      EPLUME=MPLUME*ETAL1
C**** Reduce EPLUME so that mass flux is less than mass in box
      IF (EPLUME.GT.AIRM(L)*0.975d0-MPLUME) THEN
        EPLUME=AIRM(L)*0.975d0-MPLUME
      END IF
      MPLUME=MPLUME+EPLUME
      ETAL1=EPLUME/MPOLD
      FENTR=ETAL1*FPOLD
      ENT(L)=0.001d0*FENTR/(GZL(L)*FPOLD)
      FPLUME=FPLUME+FENTR      ! to increase mass flux, remove this formula
c      FPLUME = MPLUME*BYAM(L) ! and use this instead
      FENTRA = EPLUME*BYAM(L)
      DSMR(L)=DSMR(L)-EPLUME*SUP        ! = DSM(L)-SM(L)*FENTRA
      DSMOMR(:,L)=DSMOMR(:,L)-SMOM(:,L)*FENTRA
      DQMR(L)=DQMR(L)-EPLUME*QUP        ! = DQM(L)-QM(L)*FENTRA
      DQMOMR(:,L)=DQMOMR(:,L)-QMOM(:,L)*FENTRA
      DMR(L)=DMR(L)-EPLUME
      SMP=SMP+EPLUME*SUP
      SMOMP(xymoms)= SMOMP(xymoms)+ SMOM(xymoms,L)*FENTRA
      QMP=QMP+EPLUME*QUP
      QMOMP(xymoms)= QMOMP(xymoms)+ QMOM(xymoms,L)*FENTRA
#ifdef TRACERS_ON
      DTMR(L,1:NTX) = DTMR(L,1:NTX)-TM(L,1:NTX)*FENTRA
      DTMOMR(:,L,1:NTX) = DTMOMR(:,L,1:NTX)-TMOM(:,L,1:NTX)*FENTRA
      TMP(1:NTX) = TMP(1:NTX)+TM(L,1:NTX)*FENTRA
      TMOMP(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)+TMOM(xymoms,L,1:NTX)
     *     *FENTRA
#endif

C****
C**** CONVECTIVE MOMENTUM TRANSPORT IS BASED ON GREGORY ET AL. (1997, QJRMS).
C**** THE MOMENTUM OF THE RISING PARCEL IS CHANGED BY ENTRAINMENT AND BY THE
C**** EFFECT OF THE CONVECTIVE SCALE HORIZONTAL PRESSURE GRADIENT, THE
C**** LATTER PARAMETERIZED AS PROPORTIONAL TO THE CUMULUS MASS FLUX TIMES THE
C**** VERTICAL SHEAR OF THE HORIZONTAL WIND
c****
      DO K=1,KMAX
         UMTEMP=PGRAD*MPLUME**2*(U_0(K,L+1)-U_0(K,L))/(PL(L)-PL(L+1))
         VMTEMP=PGRAD*MPLUME**2*(V_0(K,L+1)-V_0(K,L))/(PL(L)-PL(L+1))
         UMP(K)=UMP(K)+U_0(K,L)*EPLUME+UMTEMP
         DUM(K,L)=DUM(K,L)-U_0(K,L)*EPLUME-UMTEMP
         VMP(K)=VMP(K)+V_0(K,L)*EPLUME+VMTEMP
         DVM(K,L)=DVM(K,L)-V_0(K,L)*EPLUME-VMTEMP
      END DO

      END IF         ! big enough check
      END IF         ! non-zero entrainment check

C****
C**** DETRAINMENT ONLY OCCURS WHEN PARCEL LOSES BUOYANCY; MASS
C**** CALCULATED ABOVE
C****

C**** CHANGE LAYER PROPERTIES DUE TO DETRAINMENT
      IF(DET(L).GT.0.D0) THEN
        DELTA=1000.D0*DET(L)*GZL(L)
        IF(DELTA.GT..95D0) THEN
          DELTA=.95d0
          DET(L)=.001d0*DELTA/GZL(L)
        END IF
        DM(L)=DM(L)+DELTA*MPLUME
        MPLUME=MPLUME*(1.D0-DELTA)
        DSM(L)=  DSM(L)+DELTA*SMP
        SMP = SMP  *(1.-DELTA)
        DSMOM(xymoms,L)=DSMOM(xymoms,L)+DELTA*SMOMP(xymoms)
        SMOMP(xymoms) = SMOMP(xymoms)*(1.-DELTA)
        DQM(L)=  DQM(L)+DELTA*QMP
        QMP = QMP  *(1.-DELTA)
        DQMOM(xymoms,L)=DQMOM(xymoms,L)+DELTA*QMOMP(xymoms)
        QMOMP(xymoms) = QMOMP(xymoms)*(1.-DELTA)
        DO K=1,KMAX
          UMTEMP=PGRAD*MPLUME**2*(U_0(K,L+1)-U_0(K,L))/(PL(L)-PL(L+1))
          VMTEMP=PGRAD*MPLUME**2*(V_0(K,L+1)-V_0(K,L))/(PL(L)-PL(L+1))
          DUM(K,L)=DUM(K,L)+UMP(K)*DELTA-UMTEMP
          DVM(K,L)=DVM(K,L)+VMP(K)*DELTA-VMTEMP
          UMP(K)=UMP(K)-UMP(K)*DELTA+UMTEMP
          VMP(K)=VMP(K)-VMP(K)*DELTA+VMTEMP
        END DO

#ifdef TRACERS_ON
        DTM(L,1:NTX) = DTM(L,1:NTX)+DELTA*TMP(1:NTX)
        DTMOM(xymoms,L,1:NTX)=DTMOM(xymoms,L,1:NTX)+DELTA
     *        *TMOMP(xymoms,1:NTX)
        TMP(1:NTX) = TMP(1:NTX)*(1.-DELTA)
        TMOMP(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)*(1.-DELTA)
#endif
      END IF

C****
C**** CHECK THE POSSIBILITY OF DOWNDRAFT.  A DOWNDRAFT FORMS WHEN AN
c**** EQUAL MIXTURE OF CLOUDY AND CLEAR AIR IS NEGATIVELY BUOYANT
c****
C**** Define downdraft properties
C****
C**** DOWNDRAFT BEGINS 2 LEVELS ABOVE CLOUD BASE
C****
      IF(L-LMIN.GT.1) THEN

      SMIX=.5*(SUP+SMP/MPLUME)
      QMIX=.5*(QUP+QMP/MPLUME)
C     WMIX=.5*(WMUP+WMDN)
C     SVMIX=SMIX*(1.+DELTX*QMIX-WMIX)
      SVMIX=SMIX                    ! *(1.+DELTX*QMIX)
      SVUP=SUP                      ! *(1.+DELTX*QUP)
      DMMIX=(SVUP-SVMIX)*PLK(L)
C    *  +SLHE*(QSAT(SUP*PLK(L),LHX,PL(L))-QMIX)
      IF(DMMIX.LT.1d-10) CDHDRT=CDHDRT+CDHEAT(L)

C**** NO DOWNDRAFT IF BUOYANT
      IF(DMMIX.GE.1d-10) THEN       ! the mixture is negatively buoyant

      LDRAFT=L                      ! the highest downdraft level
C**** INITIATE DOWNDRAFT WITH SPECIFIED FRACTION OF UPDRAFT MASS
      ETADN=1.d0/3.d0    ! .20d0 ! reduce ETADN to improve computational stability

C**** To test with code with no downdrafts, set etadn=0. here
C**** etadn=0.  ! test

      FLEFT=1.-.5*ETADN
      DDRAFT=ETADN*MPLUME
      DDR(L)=DDRAFT
      CDHSUM1=CDHSUM1+CDHDRT*.5*ETADN      ! calculate before CDHDRT
      CDHDRT=CDHDRT-CDHDRT*.5*ETADN+CDHEAT(L)    ! SLH*COND(L)
      FDDP = .5*DDRAFT/MPLUME
      FDDL = .5*DDRAFT*BYAM(L)
      MPLUME=FLEFT*MPLUME
      SMDNL(L)=DDRAFT*SMIX
      SMOMDNL(xymoms,L)=SMOM(xymoms,L)*FDDL +  SMOMP(xymoms)*FDDP
      SMP=FLEFT*SMP
      SMOMP(xymoms)= SMOMP(xymoms)*FLEFT
      QMDNL(L)=DDRAFT*QMIX
      QMOMDNL(xymoms,L)=QMOM(xymoms,L)*FDDL +  QMOMP(xymoms)*FDDP
      QMP=FLEFT*QMP
      QMOMP(xymoms)= QMOMP(xymoms)*FLEFT
      DMR(L) = DMR(L)-.5*DDRAFT
      DSMR(L)=DSMR(L)-.5*DDRAFT*SUP        ! = DSM(L)-SM(L)*FDDL
      DSMOMR(:,L)=DSMOMR(:,L) - SMOM(:,L)*FDDL
      DQMR(L)=DQMR(L)-.5*DDRAFT*QUP        ! = DQM(L)-QM(L)*FDDL
      DQMOMR(:,L)=DQMOMR(:,L) - QMOM(:,L)*FDDL
#ifdef TRACERS_ON
      Tmdnl(l,1:NTX) = tm(l,1:NTX)*fddl+Tmp(1:NTX)*fddp
      tmomdnl(xymoms,l,1:NTX) = tmom(xymoms,l,1:NTX)*fddl+
     *     tmomp(xymoms,  1:NTX)*fddp
      dtmr    (l,1:NTX) = dtmr    (l,1:NTX)-fddl *tm    (l,1:NTX)
      dtmomr(:,l,1:NTX) = dtmomr(:,l,1:NTX)-fddl *tmom(:,l,1:NTX)
      Tmp         (1:NTX) = Tmp         (1:NTX)*fleft
      tmomp(xymoms,1:NTX) = tmomp(xymoms,1:NTX)*fleft
#endif
      DO K=1,KMAX
         UMDNL(K,L)=.5*(ETADN*UMP(K)+DDRAFT*U_0(K,L))
         UMP(K)=UMP(K)*FLEFT
         DUM(K,L)=DUM(K,L)-.5*DDRAFT*U_0(K,L)
         VMDNL(K,L)=.5*(ETADN*VMP(K)+DDRAFT*V_0(K,L))
         VMP(K)=VMP(K)*FLEFT
         DVM(K,L)=DVM(K,L)-.5*DDRAFT*V_0(K,L)
      END DO

      END IF    !   Buoyancy Test
      END IF    !   More Than 1 Layer Check
C****
C**** COMPUTE CUMULUS UPDRAFT SPEED BASED ON GREGORY (2001, QJRMS)
C****

      W2TEM=.16667D0*GRAV*BUOY(L)-WCU(L-1)*WCU(L-1)*
     *      (.66667D0*DET(L)+ENT(L))
      HDEP=AIRM(L)*TL(L)*RGAS/(GRAV*PL(L))
      WCU2(L)=WCU2(L-1)+2.*HDEP*W2TEM    ! use WCU2
C     WCU(L)=WCU(L-1)+HDEP*W2TEM
      WCU(L)=0.
      IF (WCU2(L).GT.0.D0) WCU(L)=SQRT(WCU2(L))
      IF (WCU(L).GE.0.D0) WCU(L)=MIN(50.D0,WCU(L))
      IF (WCU(L).LT.0.D0) WCU(L)=MAX(-50.D0,WCU(L))
#ifdef SCM
      if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
c     save cumulus updraft speed and plume temperature
        WCUALL(L,IC,LMIN) = WCU(L)
        MPLUMEALL(L,IC,LMIN) = CCM(L-1)
        ENTALL(L,IC,LMIN) = 1000.D0*ENT(L)
        TPALL(L,IC,LMIN) = TP   ! Save plume temp ????? THIS IS NOT CORRECT
c     write(iu_scm_prt,885) LMIN,ic,L,WCU(L), WCUALL(L,IC,LMIN)
c885  format(1x,'mstcnv  lmin ic l wcu wcuall ',
c    *      i5,i5,i5,2(f12.6))
      endif
#endif
c**** UPDATE ALL QUANTITIES CARRIED BY THE PLUME
c      SVLATL(L)=VLAT(L)
      SMPMAX=SMP
      SMOMPMAX(xymoms) =  SMOMP(xymoms)
      QMPMAX=QMP
      QMOMPMAX(xymoms) =  QMOMP(xymoms)
#ifdef TRACERS_ON
C**** Tracers at top of plume
      TMPMAX(1:NTX) = TMP(1:NTX)
      TMOMPMAX(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)
#endif
      MPMAX=MPLUME
      LMAX = LMAX + 1
      IF(WCU2(L).LT.0.D0) EXIT CLOUD_TOP  !  use WCU(L)*WCU(L)

C**** CUMULUS MICROPHYSICS
      WCUFRZ=0 ; if(LFRZ>0) WCUFRZ=WCU(LFRZ)
      CALL CONVECTIVE_MICROPHYSICS(PL(L),WCU(L),DWCU,LFRZ,WCUFRZ,TP,
     *  TI,FITMAX,PLAND,CN0,CN0I,CN0G,FLAMW,FLAMG,FLAMI,RHOIP,RHOG,
     *  ITMAX,TL(LMIN),TL(LMIN+1),WMAX,
#ifdef CLD_AER_CDNC
     *  TL(L),RCLD_C,MCDNCW,CONDMU,CONDPC(L),
#endif
     *  CONDP(L),CONDP1(L),CONDIP(L),CONDGP(L))

C**** convert condp to the same units as cond
      CONDP(L)=.01d0*CONDP(L)*CCM(L-1)*TL(L)*RGAS/PL(L)
      CONDP1(L)=.01d0*CONDP1(L)*CCM(L-1)*TL(L)*RGAS/PL(L)

C**** To test code with no precipitation from deep convection
C**** set CONDP(L)=0. here.
C**** To test code with no vertical advection of condensate
C**** set CONDP1(L)=COND(L) here.

c#ifdef SCM
c     if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
c         write(iu_scm_prt,339)  LMIN,IC,L,COND(L),CONDP1(L),CONDP(L),
c    *               CONDGP(L),CONDIP(L)
c339      format(1x,'LMIN IC L COND p1 p gp ip ',3(i5),5(f10.6))
c     endif
c#endif

C**** CALCULATE VERTICALLY TRANSPORTED PART OF CONDENSATE AND REMOVE
C**** FROM REMAINDER OF CONDENSATE (PRECIPITATING PART REMOVED LATER)
      IF(CONDP1(L).GT.COND(L)) CONDP1(L)=COND(L)
      IF(CONDP(L).GT.CONDP1(L)) CONDP(L)=CONDP1(L)
      CONDV(L)=COND(L)-CONDP1(L)       ! part of COND transported up
      COND(L)=COND(L)-CONDV(L)         ! CONDP1(L)

#ifdef TRACERS_WATER
      FQCONDV=CONDV(L)/((COND(L)+CONDV(L))+teeny)
      TRCONDV(:,L)=FQCONDV*TRCOND(:,L)
      TRCOND (:,L)=TRCOND(:,L)-TRCONDV(:,L)
#endif
      END DO  CLOUD_TOP     !   End Loop 3  (L)

      IF(LMIN.EQ.LMAX) THEN
#ifdef TRACERS_ON
        call reset_tracer_work_arrays(lmin,lmin)
#endif
        CYCLE CLOUD_TYPES
      ENDIF

      END DO AREA_PARTITION    ! end area partition loop (NPPL)

C**** UPDATE CHANGES CARRIED BY THE PLUME IN THE TOP CLOUD LAYER
      TAUMCL(LMIN:LMAX)=TAUMCL(LMIN:LMAX)+TAUMC1(LMIN:LMAX)

#ifdef SCM
      if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
        IF(PLE(LMIN)-PLE(LMAX+1).GE.450.) THEN
c       write(iu_scm_prt,887) LMIN,LMAX,IC
c887    format(1x,'mstcnv -- deep lmin lmax ic ',i5,i5,i5)
          DO L=LMIN,LMAX
            WCUDEEP(L,IC) = WCU(L)
            MPLUMEDEEP(L,IC) = CCM(L-1)
            ENTDEEP(L,IC) = 1000.D0*ENT(L)
            IF(IC.EQ.1) THEN
              SAVWL(L)=WCU(L)
            ELSE
              SAVWL1(L)=WCU(L)
            END IF
c       write(iu_scm_prt,888) ic,L,WCUDEEP(L,ic)
c888    format(1x,'mstcnv--  deep  ic l wcudeep ',
c    *                     i5,i5,f10.4)
          END DO
        END IF
      endif
#endif

      IF(PL(LMIN).LT.850.d0) THEN
        LHX1=LHE
        IF(TL(LMIN).LT.TF) LHX1=LHS
        U00L(LMIN)=1.d0-2.*(U00b*2.d-4/QSAT(TL(LMIN),LHX1,PL(LMIN)))
      ELSE
        DO L=1,LMIN
          LHX1=LHE
          IF(TL(L).LT.TF) LHX1=LHS
          U00L(L)=1.d0-2.*(U00b*2.d-4/QSAT(TL(L),LHX1,PL(L)))
        END DO
      END IF
      IF(TPSAV(LMAX).GE.TF) LFRZ=LMAX
      U00L(LMAX)=U00a

C**** Add plume characteristics at LMAX
      DM(LMAX)=DM(LMAX)+MPMAX
      DSM(LMAX)=DSM(LMAX)+SMPMAX
      DSMOM(xymoms,LMAX)=DSMOM(xymoms,LMAX) + SMOMPMAX(xymoms)
      DQM(LMAX)=DQM(LMAX)+QMPMAX
      DQMOM(xymoms,LMAX)=DQMOM(xymoms,LMAX) + QMOMPMAX(xymoms)
#ifdef TRACERS_ON
#ifdef TRACERS_NH3_UPDRAFT
       DO N=1,NTX
        select case (trname(ntix(n)))
        case('NH3')
      DTM(LMIN,N) = DTM(LMIN,N) + TMPMAX(N)
      DTMOM(xymoms,LMIN,N) = DTMOM(xymoms,LMIN,N) + TMOMPMAX(xymoms,N)

        case default
      DTM(LMAX,N) = DTM(LMAX,N) + TMPMAX(N)
      DTMOM(xymoms,LMAX,N) = DTMOM(xymoms,LMAX,N) + TMOMPMAX(xymoms,N)
        end select
       end do
#else
      DTM(LMAX,1:NTX) = DTM(LMAX,1:NTX) + TMPMAX(1:NTX)
      DTMOM(xymoms,LMAX,1:NTX) =
     *   DTMOM(xymoms,LMAX,1:NTX) + TMOMPMAX(xymoms,1:NTX)
#endif
#endif
      CCM(LMAX)=0.
      DO K=1,KMAX
         DUM(K,LMAX)=DUM(K,LMAX)+UMP(K)
         DVM(K,LMAX)=DVM(K,LMAX)+VMP(K)
      END DO
      CDHM=0.
      IF(MINLVL.GT.LMIN) MINLVL=LMIN
      IF(MAXLVL.LT.LMAX) MAXLVL=LMAX
      IF(LMCMIN.EQ.0) LMCMIN=LMIN
      IF(LMCMAX.LT.MAXLVL) LMCMAX=MAXLVL
#ifdef SCM
      if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
          PLUME_MAX(IC,LMIN) = SG_P(LMAX)
          PLUME_MIN(IC,LMIN) = SG_P(LMIN)
c         write(iu_scm_prt,899) LMIN,IC,LMAX,MINLVL,MAXLVL,
c    &         PLUME_MIN(IC,LMIN),PLUME_MAX(IC,LMIN)
c899      format(1x,'plume  lmin ic lmax ',i5,i5,i5,
c    &       ' min max ',i5,i5,' PLUME min max ',2(f10.2))
      endif
#endif
C****
C**** DOWNDRAFT DESCENT AND TRANSPORT LOOP (4)
C****
      LDMIN=LDRAFT-1
      LLMIN=LDMIN
      EDRAFT=0.

      IF(ETADN.GT.1d-10) THEN ! downdraft possible

C**** DOWNDRAFT MASS, TEMPERATURE, HUMIDITY, MOMENTUM
      DDRAFT=DDR(LDRAFT)
      DDROLD=DDRAFT
      SMDN=SMDNL(LDRAFT)
      QMDN=QMDNL(LDRAFT)
      SMOMDN(xymoms)=SMOMDNL(xymoms,LDRAFT)
      QMOMDN(xymoms)=QMOMDNL(xymoms,LDRAFT)
      DO K=1,KMAX
        UMDN(K)=UMDNL(K,LDRAFT)
        VMDN(K)=VMDNL(K,LDRAFT)
      END DO
#ifdef TRACERS_ON
      TMDN(:)=TMDNL(LDRAFT,:)
      TMOMDN(xymoms,:)=TMOMDNL(xymoms,LDRAFT,:)
#endif

C**** LOOP FROM TOP DOWN OVER POSSIBLE DOWNDRAFTS
C****
      DOWNDRAFT: DO L=LDRAFT,1,-1
      LHX=VLAT(L)               ! LHX consistency
      SLH=LHX*BYSHA
      TNX1=SMDN*PLK(L)/DDRAFT   ! save for tracers

      call get_dq_evap(smdn,qmdn,plk(l),ddraft,lhx,pl(l),cond(l)
     *     ,dqsum,fqcond1)

C**** EVAPORATE CONVECTIVE CONDENSATE IN DOWNDRAFT AND UPDATE DOWNDRAFT
C**** TEMPERATURE AND HUMIDITY; CURRENTLY ALL CONDENSATE IS ALLOWED TO
C**** EVAPORATE IF POSSIBLE (FDDRT = 1)
      DQEVP=FDDRT*COND(L)       ! limit evap from condensate to fddrt of amount
      IF(DQEVP.GT.DQSUM) DQEVP=DQSUM           ! limit evaporation
      IF(DQEVP.GT.SMDN*PLK(L)/SLH) DQEVP=SMDN*PLK(L)/SLH
      IF (L.LT.LMIN) DQEVP=0.

      FSEVP = 0
      IF (PLK(L)*SMDN.gt.teeny) FSEVP = SLH*DQEVP/(PLK(L)*SMDN)
      SMDN=SMDN-SLH*DQEVP/PLK(L)
      SMOMDN(xymoms)=SMOMDN(xymoms)*(1.-FSEVP)

      FQEVP = 0
      IF (COND(L).gt.0.) FQEVP = DQEVP/COND(L)
      QMDN=QMDN+DQEVP

c**** REMOVE EVAPORATED WATER FROM CONVECTIVE CONDENSATE AMOUNT
      COND(L)=COND(L)-DQEVP
      TAUMCL(L)=TAUMCL(L)-DQEVP*FMC1
      CDHEAT(L)=CDHEAT(L)-DQEVP*SLH
      EVPSUM=EVPSUM+DQEVP*SLH

#ifdef TRACERS_WATER
C**** RE-EVAPORATION OF TRACERS IN DOWNDRAFTS
C**** (If 100% evaporation, allow all tracers to evaporate completely.)
      IF(FQEVP.eq.1.) THEN  ! total evaporation
          TMDN(1:NTX)     = TMDN(1:NTX) + TRCOND(1:NTX,L)
#ifdef TRDIAG_WETDEPO
          IF (diag_wetdep == 1)
     &         trdvap_mc(l,1:ntx)=trdvap_mc(l,1:ntx)+trcond(1:ntx,l)
#endif
          TRCOND(1:NTX,L) = 0.d0
      ELSE            ! otherwise, tracers evaporate dependent on type of tracer
        CALL GET_EVAP_FACTOR_array(
     &         NTX,TNX1,LHX,.FALSE.,1d0,FQEVP,FQEVPT,ntix)
        dtr(1:ntx) = fqevpt(1:ntx)*trcond(1:ntx,l)
#ifdef TRDIAG_WETDEPO
        IF (diag_wetdep == 1) trdvap_mc(l,1:ntx)=trdvap_mc(l,1:ntx)
     &       +dtr(1:ntx)
#endif
        TMDN(1:NTX)     = TMDN(1:NTX)     + DTR(1:NTX)
        TRCOND(1:NTX,L) = TRCOND(1:NTX,L) - DTR(1:NTX)
      END IF
#endif
C**** ENTRAINMENT INTO DOWNDRAFTS
      IF(L.LT.LDRAFT.AND.L.GT.1) THEN
        DDRUP=DDRAFT
        DDRAFT=DDRAFT+DDROLD*ETAL(L)   ! add in entrainment
        IF(DDRUP.GT.DDRAFT) DDRAFT=DDRUP
        SMIX=SMDN/(DDRUP+teeny)
        QMIX=QMDN/(DDRUP+teeny)
        SVMIX=SMIX*PLK(L-1)              ! *(1.+DELTX*QMIX)
        SVM1=SM1(L-1)*BYAM(L-1)*PLK(L-1) ! *(1.+DELTX*QM1(L-1)*BYAM(L-1))

        IF ((SVMIX-SVM1).GE.DTMIN1) THEN
          DDRAFT=FDDET*DDRUP             ! detrain downdraft if buoyant
        END IF

C**** LIMIT SIZE OF DOWNDRAFT IF NEEDED
        IF(DDRAFT.GT..95d0*(AIRM(L-1)+DMR(L-1)))
     *    DDRAFT=.95d0*(AIRM(L-1)+DMR(L-1))
        EDRAFT=DDRAFT-DDRUP

C**** ENTRAIN INTO DOWNDRAFT, UPDATE TEMPERATURE AND HUMIDITY
        IF (EDRAFT.gt.0) THEN  ! usual case, entrainment into downdraft
          FENTRA=EDRAFT*BYAM(L)
          SENV=SM(L)*BYAM(L)
          QENV=QM(L)*BYAM(L)
          SMDN=SMDN+EDRAFT*SENV
          QMDN=QMDN+EDRAFT*QENV
          SMOMDN(xymoms)= SMOMDN(xymoms)+ SMOM(xymoms,L)*FENTRA
          QMOMDN(xymoms)= QMOMDN(xymoms)+ QMOM(xymoms,L)*FENTRA
          DSMR(L)=DSMR(L)-EDRAFT*SENV
          DSMOMR(:,L)=DSMOMR(:,L)-SMOM(:,L)*FENTRA
          DQMR(L)=DQMR(L)-EDRAFT*QENV
          DQMOMR(:,L)=DQMOMR(:,L)-QMOM(:,L)*FENTRA
          DMR(L)=DMR(L)-EDRAFT
C**** EFFECT OF ENTRAINMENT AND CONVECTIVE PRESSURE GRADIENT ON
C**** HORIZONTAL MOMENTUM OF DOWNDRAFT AIR
          DO K=1,KMAX        ! add in momentum entrainment
            UMTEMP=PGRAD*DDRAFT**2*(U_0(K,L+1)-U_0(K,L))/(PL(L)-PL(L+1))
            VMTEMP=PGRAD*DDRAFT**2*(V_0(K,L+1)-V_0(K,L))/(PL(L)-PL(L+1))
            UMDN(K)=UMDN(K)+FENTRA*UM(K,L)-UMTEMP
            VMDN(K)=VMDN(K)+FENTRA*VM(K,L)-VMTEMP
            DUM(K,L)=DUM(K,L)-FENTRA*UM(K,L)+UMTEMP
            DVM(K,L)=DVM(K,L)-FENTRA*VM(K,L)+VMTEMP
          END DO
#ifdef TRACERS_ON
          Tenv(1:NTX)=tm(l,1:NTX)/airm(l)
          TMDN(1:NTX)=TMDN(1:NTX)+EDRAFT*Tenv(1:NTX)
          TMOMDN(xymoms,1:NTX)= TMOMDN(xymoms,1:NTX)+ TMOM(xymoms,L
     *         ,1:NTX)*FENTRA
          DTMR(L,1:NTX)=DTMR(L,1:NTX)-EDRAFT*TENV(1:NTX)
          DTMOMR(:,L,1:NTX)=DTMOMR(:,L,1:NTX)-TMOM(:,L,1:NTX)*FENTRA
#endif
        ELSE  ! occasionally detrain into environment if ddraft too big
          FENTRA=EDRAFT/(DDRUP+teeny)  ! < 0
          DSM(L)=DSM(L)-FENTRA*SMDN
          DSMOM(xymoms,L)=DSMOM(xymoms,L)-SMOMDN(xymoms)*FENTRA
          DQM(L)=DQM(L)-FENTRA*QMDN
          DQMOM(xymoms,L)=DQMOM(xymoms,L)-QMOMDN(xymoms)*FENTRA
          SMDN=SMDN*(1+FENTRA)
          QMDN=QMDN*(1+FENTRA)
          SMOMDN(xymoms)= SMOMDN(xymoms)*(1+FENTRA)
          QMOMDN(xymoms)= QMOMDN(xymoms)*(1+FENTRA)
          DM(L)=DM(L)-EDRAFT       ! EDRAFT is negative here
          DO K=1,KMAX          ! add in momentum detrainment
            UMTEMP=PGRAD*DDRAFT**2*(U_0(K,L+1)-U_0(K,L))/(PL(L)-PL(L+1))
            VMTEMP=PGRAD*DDRAFT**2*(V_0(K,L+1)-V_0(K,L))/(PL(L)-PL(L+1))
            UMDN(K)=UMDN(K)*(1.+FENTRA)-UMTEMP
            VMDN(K)=VMDN(K)*(1.+FENTRA)-VMTEMP
            DUM(K,L)=DUM(K,L)-FENTRA*UMDN(K)+UMTEMP
            DVM(K,L)=DVM(K,L)-FENTRA*VMDN(K)+VMTEMP
          END DO
#ifdef TRACERS_ON
          DTM(L,1:NTX)=DTM(L,1:NTX)-FENTRA*TMDN(1:NTX)
          DTMOM(xymoms,L,1:NTX)=DTMOM(xymoms,L,1:NTX)-
     *         TMOMDN(xymoms,1:NTX)*FENTRA
          TMDN(1:NTX)=TMDN(1:NTX)*(1.+FENTRA)
          TMOMDN(xymoms,1:NTX)= TMOMDN(xymoms,1:NTX)*(1.+FENTRA)
#endif
        END IF
      END IF

      LDMIN=L
      LLMIN=LDMIN        ! save LDMIN for diagnostics
C**** ALLOW DOWNDRAFT TO DESCEND BELOW CLOUD BASE IF IT IS NEGATIVELY BUOYANT
      IF (L.GT.1) THEN
        SMIX=SMDN/(DDRAFT+teeny)
        QMIX=QMDN/(DDRAFT+teeny)
        SVMIX=SMIX*PLK(L-1)              ! *(1.+DELTX*QMIX)
        SVM1=SM1(L-1)*BYAM(L-1)*PLK(L-1) ! *(1.+DELTX*QM1(L-1)*BYAM(L-1))
        IF (L.LE.LMIN.AND.SVMIX.GE.SVM1) EXIT
        DDM(L-1)=DDRAFT
        DDROLD=DDRAFT
        DDRAFT=DDRAFT+DDR(L-1)    ! add in downdraft one layer below
        SMDN=SMDN+SMDNL(L-1)
        QMDN=QMDN+QMDNL(L-1)
        SMOMDN(xymoms)=SMOMDN(xymoms)+SMOMDNL(xymoms,L-1)
        QMOMDN(xymoms)=QMOMDN(xymoms)+QMOMDNL(xymoms,L-1)
        DO K=1,KMAX
        UMDN(K)=UMDN(K)+UMDNL(K,L-1)
        VMDN(K)=VMDN(K)+VMDNL(K,L-1)
        END DO
#ifdef TRACERS_ON
        DO N=1,NTX
          TMDN(N)=TMDN(N)+TMDNL(L-1,N)
          TMOMDN(xymoms,N)=TMOMDN(xymoms,N)+TMOMDNL(xymoms,L-1,N)
        END DO
#endif
      END IF

      END DO DOWNDRAFT
C**** end of loop over downdrafts

C**** UPDATE PROPERTIES OF LOWEST LAYER TO WHICH DOWNDRAFT PENETRATES
      DSM(LDMIN)=DSM(LDMIN)+SMDN
      DSMOM(xymoms,LDMIN)=DSMOM(xymoms,LDMIN) + SMOMDN(xymoms)
      DQM(LDMIN)=DQM(LDMIN)+QMDN

      DQMOM(xymoms,LDMIN)=DQMOM(xymoms,LDMIN) + QMOMDN(xymoms)
      TDNL(LDMIN)=SMDN*PLK(LDMIN)/(DDRAFT+teeny)
      QDNL(LDMIN)=QMDN/(DDRAFT+teeny)
#ifdef TRACERS_ON
      DTM(LDMIN,1:NTX) = DTM(LDMIN,1:NTX) + TMDN(1:NTX)
      DTMOM(xymoms,LDMIN,1:NTX) = DTMOM(xymoms,LDMIN,1:NTX) +
     *     TMOMDN(xymoms,1:NTX)
      TRDNL(1:NTX,LDMIN)=TMDN(1:NTX)/(DDRAFT+teeny)
#endif
      DO K=1,KMAX
        DUM(K,LDMIN)=DUM(K,LDMIN)+UMDN(K)
        DVM(K,LDMIN)=DVM(K,LDMIN)+VMDN(K)
      END DO
      DM(LDMIN)=DM(LDMIN)+DDRAFT
      END IF

C****
C**** SUBSIDENCE AND MIXING LOOP (5)
C****
C**** Calculate vertical mass fluxes (Note CM for subsidence is defined
C**** in opposite sense than normal (positive is down))
      IF(LDMIN.GT.LMIN) LDMIN=LMIN    ! some loops require LMIN to LMAX
      DO L=0,LDMIN-1
        CM(L) = 0.
      END DO
      DO L=LDMIN,LMAX
        CM(L) = CM(L-1) - DM(L) - DMR(L)
        SMT(L)=SM(L)    ! Save profiles for diagnostics
        QMT(L)=QM(L)
      END DO
      DO L=LMAX,LM
        CM(L) = 0.
      END DO
C**** simple upwind scheme for momentum
      DO K=1,KMAX
        SUMU(K)=sum(UM(K,LDMIN:LMAX))
        SUMV(K)=sum(VM(K,LDMIN:LMAX))
      END DO
      SUMDP=sum(AIRM(LDMIN:LMAX))
      ALPHA=0.
      DO L=LDMIN,LMAX
        CLDM=CCM(L)
        IF(L.LT.LDRAFT.AND.L.GE.LLMIN.AND.ETADN.GT.1d-10)
     *       CLDM=CCM(L)-DDM(L)
        IF(MC1) VSUBL(L)=100.*CLDM*RGAS*TL(L)/(PL(L)*GRAV*DTsrc)
        BETA=CLDM*BYAM(L+1)
        IF(CLDM.LT.0.) BETA=CLDM*BYAM(L)
        BETAU=BETA
        ALPHAU=ALPHA
        IF(BETA.LT.0.) BETAU=0.
        IF(ALPHA.LT.0.) ALPHAU=0.
        DO K=1,KMAX
          UM(K,L)=
     *         UM(K,L)+RA(K)*(-ALPHAU*UM(K,L)+BETAU*UM(K,L+1)+DUM(K,L))
          VM(K,L)=
     *         VM(K,L)+RA(K)*(-ALPHAU*VM(K,L)+BETAU*VM(K,L+1)+DVM(K,L))
        END DO
        ALPHA=BETA
      END DO
      DO K=1,KMAX
        SUMU1(K)=sum(UM(K,LDMIN:LMAX))
        SUMV1(K)=sum(VM(K,LDMIN:LMAX))
      END DO
      DO K=1,KMAX                          ! momentum adjustment
        UM(K,LDMIN:LMAX)=UM(K,LDMIN:LMAX)-(SUMU1(K)-SUMU(K))*
     *    AIRM(LDMIN:LMAX)/SUMDP
        VM(K,LDMIN:LMAX)=VM(K,LDMIN:LMAX)-(SUMV1(K)-SUMV(K))*
     *    AIRM(LDMIN:LMAX)/SUMDP
      END DO

C****

c Determine the number of subsidence sub-timesteps such that
c courant numbers in the QUS do not exceed 1
      ksub = 1
      do l=ldmin,lmax-1
        if(    +cm(l) > airm(l+1)+dmr(l+1)) then
          ksub = max(ksub, 1+int((+cm(l)-dmr(l+1))/airm(l+1)) )
        elseif(-cm(l) > airm(l  )+dmr(l  )) then
          ksub = max(ksub, 1+int((-cm(l)-dmr(l  ))/airm(l  )) )
        endif
      enddo
      ksub = min(ksub,2) ! max 2 iterations allowed currently
c      ksub = 2 ! non-interactive default

      byksub = 1d0/ksub
      nsub = lmax-ldmin+1
      cmneg=0.          ! initialization
      cmneg(ldmin:lmax-1) = -cm(ldmin:lmax-1)*byksub
      cmneg(lmax) = 0.  ! avoid roundoff error (esp. for qlimit)

C**** Subsidence uses Quadratic Upstream Scheme for QM and SM
      DO ITER=1,KSUB ! subsidence sub-timesteps to improve stability
      ML(LDMIN:LMAX) = AIRM(LDMIN:LMAX) +   DMR(LDMIN:LMAX)*BYKSUB
      SM(LDMIN:LMAX) =  SM(LDMIN:LMAX)  +  DSMR(LDMIN:LMAX)*BYKSUB
      SMOM(:,LDMIN:LMAX)=SMOM(:,LDMIN:LMAX)+DSMOMR(:,LDMIN:LMAX)*BYKSUB
      call adv1d(sm(ldmin),smom(1,ldmin), f(ldmin),fmom(1,ldmin),
     &     ml(ldmin),cmneg(ldmin), nsub,.false.,1, zdir,ierrt,lerrt)
      SM(LDMIN:LMAX) =   SM(LDMIN:LMAX) +   DSM(LDMIN:LMAX)*BYKSUB
      SMOM(:,LDMIN:LMAX)=SMOM(:,LDMIN:LMAX)+DSMOM(:,LDMIN:LMAX)*BYKSUB
      ierr=max(ierrt,ierr) ; lerr=max(lerrt+ldmin-1,lerr)
C****
      ML(LDMIN:LMAX) = AIRM(LDMIN:LMAX) +   DMR(LDMIN:LMAX)*BYKSUB
      QM(LDMIN:LMAX) =   QM(LDMIN:LMAX) +  DQMR(LDMIN:LMAX)*BYKSUB
      QMOM(:,LDMIN:LMAX)=QMOM(:,LDMIN:LMAX)+DQMOMR(:,LDMIN:LMAX)*BYKSUB
      call adv1d(qm(ldmin),qmom(1,ldmin), f(ldmin),fmom(1,ldmin),
     &     ml(ldmin),cmneg(ldmin), nsub,.true.,1, zdir,ierrt,lerrt)
      QM(LDMIN:LMAX) =   QM(LDMIN:LMAX) +   DQM(LDMIN:LMAX)*BYKSUB
      QMOM(:,LDMIN:LMAX)=QMOM(:,LDMIN:LMAX)+DQMOM(:,LDMIN:LMAX)*BYKSUB
      ierr=max(ierrt,ierr) ; lerr=max(lerrt+ldmin-1,lerr)
#ifdef TRACERS_ON
C**** Subsidence of tracers by Quadratic Upstream Scheme
      DO N=1,NTX
      ML(LDMIN:LMAX) =  AIRM(LDMIN:LMAX) +    DMR(LDMIN:LMAX)*BYKSUB
      TM(LDMIN:LMAX,N) =  TM(LDMIN:LMAX,N) + DTMR(LDMIN:LMAX,N)*BYKSUB
      TMOM(:,LDMIN:LMAX,N) = TMOM(:,LDMIN:LMAX,N)+DTMOMR(:,LDMIN:LMAX,N)
     &                      *BYKSUB
      call adv1d(tm(ldmin,n),tmom(1,ldmin,n), f(ldmin),fmom(1,ldmin),
     &     ml(ldmin),cmneg(ldmin), nsub,t_qlimit(n),1, zdir,ierrt,lerrt)
      TM(LDMIN:LMAX,N) = TM(LDMIN:LMAX,N) +   DTM(LDMIN:LMAX,N)*BYKSUB
      TMOM(:,LDMIN:LMAX,N) = TMOM(:,LDMIN:LMAX,N) +DTMOM(:,LDMIN:LMAX,N)
     &                      *BYKSUB
      ierr=max(ierrt,ierr) ; lerr=max(lerrt+ldmin-1,lerr)
      END DO
#endif
      END DO        ! end of sub-timesteps for subsidence
C**** Check for v. rare negative humidity error condition
      DO L=LDMIN,LMAX
        IF(QM(L).LT.0.d0) then
          WRITE(6,*) ' Q neg: it,i,j,l,q,cm',itime,i_debug,j_debug,l
     $          ,qm(l),cmneg(l)
C**** reduce subsidence post hoc.
          LM1=max(1,L-1)
          IF (QM(LM1)+QM(L).lt.0) then
            write(6,*) "Q neg cannot be fixed!",L,QM(LM1:L)
          ELSE
            QM(L-1)=QM(L-1)+QM(L)
            QM(L)=0.
#ifdef TRACERS_WATER
C**** corresponding water tracer adjustment
            DO N=1,NTX
              IF (tr_wd_type(n) .eq. nWater) then
                TM(L-1,N)=TM(L-1,N)+TM(L,N)
                TM(L,N)=0.
              END IF
            END DO
#endif
          END IF
        END IF
      END DO
#ifdef TRACERS_ON
C**** check for independent tracer errors
      DO N=1,NTX
        IF (.not.t_qlimit(n)) cycle
        DO L=LDMIN,LMAX
          IF (TM(L,N).lt.0.) then
            WRITE(6,*) trname(n),' neg: it,i,j,l,tr,cm',itime,i_debug
     $             ,j_debug,l,tm(l,n),cmneg(l)
C**** reduce subsidence post hoc.
            LM1=max(1,L-1)
            IF (TM(LM1,N)+TM(L,N).lt.0) THEN
              write(6,*) trname(n)," neg cannot be fixed!",L,TM(LM1:L,N)
            ELSE
              TM(L-1,N)=TM(L-1,N)+TM(L,N)
              TM(L,N)=0.
            END IF
          END IF
        END DO
      END DO
#endif
C**** diagnostics
      DO L=LDMIN,LMAX
        FCDH=0.
        IF(L.EQ.LMAX) FCDH=CDHSUM-CDHSUM1+CDHM
        FCDH1=0.
        IF(L.EQ.LLMIN) FCDH1=CDHSUM1-EVPSUM
        MCFLX(L)=MCFLX(L)+CCM(L)*FMC1
        DGDSM(L)=DGDSM(L)+(PLK(L)*(SM(L)-SMT(L))-FCDH-FCDH1)*FMC1
        IF(PLE(LMAX+1).GT.700.d0) DGSHLW(L)=DGSHLW(L)+
     *    (PLK(L)*(SM(L)-SMT(L))-FCDH-FCDH1)*FMC1
        IF(PLE(LMIN)-PLE(LMAX+1).GE.450.d0) DGDEEP(L)=DGDEEP(L)+
     *    (PLK(L)*(SM(L)-SMT(L))-FCDH-FCDH1)*FMC1
        DTOTW(L)=DTOTW(L)+SLHE*(QM(L)-QMT(L)+COND(L))*FMC1
        DGDQM(L)=DGDQM(L)+SLHE*(QM(L)-QMT(L))*FMC1
        DDMFLX(L)=DDMFLX(L)+DDM(L)*FMC1
#ifdef SCM
        if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
            CUMFLX(L,IC,LMIN) = 100.*MCFLX(L)*bygrav/dtsrc
            DWNFLX(L,IC,LMIN) = 100.*DDMFLX(L)*bygrav/dtsrc
c           write(iu_scm_prt,*) 'LMIN IC L CUMFLX DWNFLX ',
c    &            LMIN,IC,L,CUMFLX(L,IC,LMIN),DWNFLX(L,IC,LMIN)
        endif
#endif
      END DO
C**** save new 'environment' profile for static stability calc.
      DO L=1,LM
        SM1(L)=SM(L)
        QM1(L)=QM(L)
      END DO
#ifdef TRACERS_ON
      TM1(:,1:NTX) = TM(:,1:NTX)
#endif

C****
C**** Partition condensate into precipitation and cloud water
C****
      COND(LMAX)=COND(LMAX)+CONDV(LMAX)
#ifdef TRACERS_WATER
      TRCOND(:,LMAX)=TRCOND(:,LMAX)+TRCONDV(:,LMAX)
#endif
      IF(PLE(LMIN)-PLE(LMAX+1).GE.450.) THEN
        DO L=LMAX,LMIN,-1
          IF(COND(L).LT.CONDP(L)) CONDP(L)=COND(L)
          FCLW=0.
          IF (COND(L).GT.0) FCLW=(COND(L)-CONDP(L))/COND(L)
#ifdef SCM
          if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
cccc  in g/m3
c            PRCCDEEP(L,IC,LMIN) = FMC1*1.d5*CONDP(L) *         ! save COND
c    *             BYAM(L)*PL(L)/(RGAS*TL(L))
c            NPRCCDEEP(L,IC,LMIN) = FMC1*1.d5*(COND(l)-CONDP(L)) *
c    *             BYAM(L)*PL(L)/(RGAS*TL(L))
cccc  in kg/kg
             PRCCDEEP(L,IC,LMIN) = FMC1*CONDP(L)*BYAM(L)       ! save COND
             NPRCCDEEP(L,IC,LMIN) = FMC1*(COND(l)-CONDP(L))*BYAM(L)
             PRCCGRP(L,IC,LMIN) = FMC1*CONDGP(L)*BYAM(L)
             PRCCICE(L,IC,LMIN) = FMC1*CONDIP(L)*BYAM(L)
             DETRAINDEEP(L,IC,LMIN) = FCLW*COND(L)*BYAM(L)*FMC1
c            write(iu_scm_prt,889) lmin,ic,FMC1,L,
c    *              PRCCDEEP(L,IC,LMIN)*1000.,
c    *              NPRCCDEEP(L,IC,LMIN)*1000.,TPALL(L,IC,LMIN),
c    *              PRCCGRP(L,IC,LMIN)*1000.,
c    *              PRCCICE(L,IC,LMIN)*1000.
c889         format(1x,
c    *        'mc--dp lmin ic fmc1 l prcc nprcc tp grp ice',
c    *              2(i4),f8.3,i4,f9.5,f9.5,f8.2,f9.5,f9.5)
          endif
#endif
C**** check in case saved condensate is different phase
          if (SVLATL(L).gt.0 .and. SVLATL(L).ne.VLAT(L)) HEAT1(L)=
     *      HEAT1(L)+(SVLATL(L)-VLAT(L))*svwmxl(l)*airm(l)*BYSHA/FMC1
          SVLATL(L)=VLAT(L)         ! moved from above
          SVWMXL(L)=SVWMXL(L)+FCLW*COND(L)*BYAM(L)*FMC1
          COND(L)=CONDP(L)

#ifdef TRACERS_WATER
C**** Apportion cloud tracers and condensation
C**** Note that TRSVWML is in mass units unlike SVWMX
          TRSVWML(1:NTX,L) = TRSVWML(1:NTX,L) + FCLW*TRCOND(1:NTX,L)
     *         *FMC1
#ifdef TRDIAG_WETDEPO
          IF (diag_wetdep == 1)
     &        trflcw_mc(l,1:ntx)=trflcw_mc(l,1:ntx)+fclw*trcond(1:ntx,l)
#endif
          TRCOND(1:NTX,L) = (1.-FCLW)*TRCOND(1:NTX,L)
#endif
        END DO
      END IF

C****
C**** RE-EVAPORATION AND PRECIPITATION LOOP (6); INCLUDES CALCULATION OF
C**** CONVECTIVE CLOUD FRACTION
C****

#ifdef SCM
      if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
          do L=1,LM
             MCCOND(L,IC,LMIN) = COND(L)*FMC1*BYAM(L)
          enddo
      endif
#endif
      PRCP=COND(LMAX)
      PRHEAT=CDHEAT(LMAX)
#ifdef CFMIP
      DO L=1,LM
        cfmip_ccl(L) = COND(L)*FMC1*BYAM(L)
      ENDDO
#endif
C**** check whether environment is the same phase as cond
      TOLD=SMOLD(LMAX)*PLK(LMAX)*BYAM(LMAX)
      lhp(lmax)=lhe
      if (TOLD.le.TF) lhp(lmax)=lhs

C**** adjust phase of precip to that of environment
      if ((TOLD.gt.TF .and. vlat(Lmax).eq.lhs) .or. (TOLD.LE.TF .and.
     *     vlat(lmax).eq.lhe)) then
        FSSUM = 0
        IF (ABS(PLK(LMAX)*SM(LMAX)).gt.teeny .and. ((lhp(lmax)-vlat(lmax
     *       ))*PRCP*BYSHA).lt.0) FSSUM = -(lhp(lmax)-vlat(lmax))*PRCP
     *       *BYSHA/(PLK(LMAX)*SM(LMAX))
c        if (debug) print*,"cnv0",lmax,(lhp(lmax)-vlat(lmax))*PRCP*BYSHA
c     $       /PLK(LMAX),lhp(lmax),vlat(lmax)
        SM(LMAX)=SM(LMAX)+(lhp(lmax)-vlat(lmax))*PRCP*BYSHA/PLK(LMAX)
        SMOM(:,LMAX) =  SMOM(:,LMAX)*(1.-FSSUM)
      end if

C**** add in heat1 from lmax
c      if (heat1(lmax).ne.0) print*,"cnvA",i_debug,j_debug,lmax
c     $     ,heat1(lmax),vlat(lmax),lhp(lmax),prcp

c      FSSUM = 0
c      IF (ABS(PLK(LMAX)*SM(LMAX)).gt.teeny .and. HEAT1(lmax).gt.0) FSSUM
c     *     = -heat1(lmax)/(PLK(LMAX)*SM(LMAX))
c      SM(LMAX)=SM(LMAX)-heat1(lmax)/PLK(LMAX)
c      SMOM(:,LMAX) =  SMOM(:,LMAX)*(1.-FSSUM)


#ifdef TRACERS_WATER
C**** Tracer precipitation
C Note that all of the tracers that condensed do not precipitate here,
C since a fraction (FCLW) of TRCOND was removed above.
      TRPRCP(1:NTX) = TRCOND(1:NTX,LMAX)
#ifdef TRDIAG_WETDEPO
      IF (diag_wetdep == 1)
     &     trprcp_mc(lmax,1:ntx)=trprcp_mc(lmax,1:ntx)
     &     +trcond(1:ntx,lmax)
#endif
#endif
         DPHASE(LMAX)=DPHASE(LMAX)+(CDHSUM-CDHSUM1+
     *                CDHM)*FMC1
         IF(PLE(LMAX+1).GT.700.d0) DPHASHLW(LMAX)=DPHASHLW(LMAX)+
     *     (CDHSUM-CDHSUM1+CDHM)*FMC1
         IF(PLE(LMIN)-PLE(LMAX+1).GE.450.d0) DPHADEEP(LMAX)=
     *     DPHADEEP(LMAX)+(CDHSUM-CDHSUM1+CDHM)*FMC1

C**** Loop down from top of plume
      EVAP_PRECIP: DO L=LMAX-1,1,-1
#ifdef SCM
C     if (SCM_ATURB_FLAG.eq.0) then  ! commented out by yao
c         for SCM running with Dry Convection instead of ATURB WTURB
c         not filled -- set TRATIO=1
C         TRATIO = 1.d0
C     else
c         for SCM running with ATURB
C         TTURB=HPBL/WTURB(L)
C         TRATIO=TTURB/DTsrc
C     endif
#else
C     TTURB=HPBL/WTURB(L)
C     TRATIO=TTURB/DTsrc
#endif

C**** CALCULATE CONVECTIVE CLOUD FRACTION
      CALL MC_CLOUD_FRACTION(L,TL(L),PL(L),CCM(L),WCU(L),PLE(LMIN),
     *  PLE(L+2),PLE(LMAX+1),CCM(LMIN),WCU(LMIN),LMIN,LMAX,
     *  CCMUL,CCMUL2,FCLOUD)

C**** PRECIPITATION IS ALLOWED TO EVAPORATE FULLY INTO A FRACTION OF
C**** THE GRIDBOX HALF AS LARGE AS THE FRACTION OF GRIDBOX MASS THAT
C**** CONVECTS
      FEVAP=.5*CCM(L)*BYAM(L+1)
      IF(L.LT.LMIN) FEVAP=.5*CCM(LMIN)*BYAM(LMIN+1)
      IF(FEVAP.GT..5) FEVAP=.5
      CLDMCL(L+1)=MIN(CLDMCL(L+1)+FCLOUD*FMC1,FMC1)
      CLDREF=CLDMCL(L+1)
      IF(PLE(LMAX+1).GT.700..AND.CLDREF.GT.CLDSLWIJ)
     *  CLDSLWIJ=CLDREF
      IF(PLE(LMIN)-PLE(LMAX+1).GE.450..AND.CLDREF.GT.CLDDEPIJ)
     *  CLDDEPIJ=CLDREF
      TOLD=SMOLD(L)*PLK(L)*BYAM(L)
      TOLD1=SMOLD(L+1)*PLK(L+1)*BYAM(L+1)

C**** FORWARD STEP COMPUTES HUMIDITY CHANGE BY RECONDENSATION
C**** Q = Q + F(TOLD,PRHEAT,QOLD+EVAP)
      PRECNVL(L+1)=PRECNVL(L+1)+PRCP*BYGRAV

C**** CALCULATE CONVECTIVE PRECIP PHASE
      CALL MC_PRECIP_PHASE(L,AIRM(L),FEVAP,
     *  LHP(L+1),PRCP,TOLD,TOLD1,VLAT(L),COND(L),LMIN,
     *  LHP(L),MCLOUD,HEAT1(L))

C**** set phase of precip based on local environment temperature
      LHX=LHP(L)
#ifdef SCM
c     save MC precip phase
      if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
          LHPMC(L) = LHP(L)
c         write(iu_scm_prt,*) 'L LHPMC ',L,LHPMC(L)
      endif
#endif

      DQSUM=0.
      FPRCP=0.
      SLH=LHX*BYSHA

      IF (PRCP.GT.0.) THEN

      if (mcloud.gt.0) call get_dq_evap(smold(l),qmold(l),plk(l),airm(l)
     *     ,lhx,pl(l),prcp*AIRM(L)/MCLOUD,dqsum,fprcp)
      dqsum=dqsum*MCLOUD*BYAM(L)

      PRCP=PRCP-DQSUM
      QM(L)=QM(L)+DQSUM

      END IF

C**** UPDATE TEMPERATURE DUE TO NET REEVAPORATION IN CLOUDS
      FSSUM = 0
      IF (ABS(PLK(L)*SM(L)).gt.teeny .and. (SLH*DQSUM+HEAT1(L)).gt.0)
     *     FSSUM = (SLH*DQSUM+HEAT1(L))/(PLK(L)*SM(L))
c      if (debug) print*,"cnv4",l,SLH*DQSUM,HEAT1(L)
      SM(L)=SM(L)-(SLH*DQSUM+HEAT1(L))/PLK(L)
      SMOM(:,L) =  SMOM(:,L)*(1.-FSSUM)

         FCDH1=0.
         IF(L.EQ.LLMIN) FCDH1=CDHSUM1-EVPSUM
         DPHASE(L)=DPHASE(L)-(SLH*DQSUM-FCDH1+HEAT1(L))*FMC1
         IF(PLE(LMAX+1).GT.700.d0) DPHASHLW(L)=DPHASHLW(L)-
     *     (SLH*DQSUM-FCDH1+HEAT1(L))*FMC1
         IF(PLE(LMIN)-PLE(LMAX+1).GE.450.d0) DPHADEEP(L)=DPHADEEP(L)-
     *     (SLH*DQSUM-FCDH1+HEAT1(L))*FMC1
         DQCOND(L)=DQCOND(L)-SLH*DQSUM*FMC1

#ifdef TRACERS_WATER
      IF (PRCP+DQSUM.GT.0.) THEN
C**** Tracer net re-evaporation
C**** (If 100% evaporation, allow all tracers to evaporate completely.)
      BELOW_CLOUD = L.lt.LMIN
      IF(FPRCP.eq.1.) THEN      !total evaporation
#ifdef TRDIAG_WETDEPO
        IF (diag_wetdep == 1) trnvap_mc(l,1:ntx) = trnvap_mc(l,1:ntx)
     &       +trprcp(1:ntx)
#endif
        DO N=1,NTX
          TM(L,N)   = TM(L,N)  + TRPRCP(N)
c         if (debug .and.n.eq.1) print*,"cld2",L,TM(L,N),TRPRCP(N),2
c     *         *FEVAP
          TRPRCP(N) = 0.d0
        END DO
      ELSE ! otherwise, tracers evaporate dependent on type of tracer
C**** estimate effective humidity
        if (below_cloud) then
          TNX1=(SM(L)*PLK(L)-SLH*DQSUM*(1./(2.*MCLOUD)-1.))*BYAM(L)
          HEFF=MIN(1d0,(QM(L)+DQSUM*(1./(2.*MCLOUD)-1.))*BYAM(L)
     *         /QSAT(TNX1,LHX,PL(L)))
        else
          heff=1.
        end if
        CALL GET_EVAP_FACTOR_array(
     &         NTX,TOLD,LHX,BELOW_CLOUD,HEFF,FPRCP,FPRCPT,ntix)
        dtr(1:ntx) = fprcpt(1:ntx)*trprcp(1:ntx)
#ifdef TRDIAG_WETDEPO
        IF (diag_wetdep == 1) trnvap_mc(l,1:ntx)=trnvap_mc(l,1:ntx)
     &       +dtr(1:ntx)
#endif
        TM(L,1:NTX) = TM(L,1:NTX)     + DTR(1:NTX)
c          if (debug .and.n.eq.1) print*,"cld3",L,TM(L,N),FPRCP
c     *         ,FPRCPT(N),TRPRCP(N)
        TRPRCP(1:NTX) = TRPRCP(1:NTX) - DTR(1:NTX)
      END IF
#ifndef NO_WASHOUT_IN_CLOUDS
      IF (.NOT. below_cloud .AND. prcp > teeny) THEN
c**** Washout of tracers in cloud
        wmxtr=prcp*byam(l)
        precip_mm=prcp*100.*bygrav
        b_beta_DT=fplume
        TM_dum(:) = TM(L,:)
        CALL GET_WASH_factor_array(
     &       ntx,b_beta_dt,precip_mm,fwasht,told,lhx,
     &   wmxtr,fplume,tm_dum,trprcp,thwash,pl(l),ntix,.true.)
        dtr(1:ntx) = fwasht(1:ntx)*tm_dum(1:ntx)
#ifdef TRDIAG_WETDEPO
        IF (diag_wetdep == 1) trwash_mc(l,1:ntx)=trwash_mc(l,1:ntx)
     &       +dtr(1:ntx)+thwash(1:ntx)
#endif
        do igas=1,gases_count
          n = gases_list(igas)
          IF (tm(l,n) > teeny) THEN
            tmfac(n)=thwash(n)/tm(l,n)
          ELSE
            tmfac(n)=0.
          END IF
        enddo
        DO n=1,ntx
          trprcp(n)=dtr(n)+trprcp(n)+thwash(n)
          tm(l,n)=tm(l,n)*(1.-fwasht(n))-thwash(n)
          tmom(xymoms,l,n)=tmom(xymoms,l,n)*(1.-fwasht(n)-tmfac(n))
        END DO
      ELSE IF (below_cloud .AND. prcp > teeny) THEN
#else
      IF (below_cloud .AND. prcp > teeny) THEN
#endif
C**** WASHOUT of TRACERS BELOW CLOUD
        WMXTR = PRCP*BYAM(L)
        precip_mm = PRCP*100.*bygrav
        b_beta_DT = FPLUME
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
        WA_VOL= precip_mm*DXYPIJ

        CALL GET_SULFATE(L,TOLD,FPLUME,WA_VOL,WMXTR,SULFIN,
     *       SULFINC,SULFOUT,TR_LEFT,TM,TRPRCP,AIRM,LHX,
     *       DT_SULF_MC(1,L),CLDSAVT)

        do iaqch=1,aqchem_count
          n = aqchem_list(iaqch)
          TRPRCP(N)=TRPRCP(N)*(1.+SULFINC(N))
          TM(L,N)=TM(L,N)*(1.+SULFIN(N))
          TMOM(xymoms,L,N)=TMOM(xymoms,L,N) *(1.+SULFIN(N))
          TRCOND(N,L) = TRCOND(N,L)+SULFOUT(N)
        enddo
#endif

cdmk Here I took out GET_COND, since we are below cloud.
cdmk GET_WASH now has gas dissolution, extra arguments
        TM_dum(:) = TM(L,:)
        CALL GET_WASH_FACTOR_array(NTX,b_beta_DT,precip_mm,FWASHT,
     *  TOLD,LHX,WMXTR,FPLUME,TM_dum,TRPRCP,THWASH,pl(l),ntix,.true.)
        dtr(1:ntx) = fwasht(1:ntx)*tm_dum(1:ntx)
#ifdef TRDIAG_WETDEPO
        IF (diag_wetdep == 1) trwash_mc(l,1:ntx)=trwash_mc(l,1:ntx)
     &       +dtr(1:ntx)+thwash(1:ntx)
#endif
        do igas=1,gases_count
          n = gases_list(igas)
          IF (TM(L,N).GT.teeny) THEN
            TMFAC(N)=THWASH(N)/TM(L,N)
          ELSE
            TMFAC(N)=0.
          END IF
        enddo
        DO N=1,NTX
          TRPRCP(N) = DTR(N)+TRPRCP(N)+THWASH(N)
          TM(L,N)=TM(L,N)*(1.-FWASHT(N))-THWASH(N)
c          if (debug .and.n.eq.1) print*,"cld4",L,TM(L,N),FWASHT(N)
c     *         ,THWASH(N)
          TMOM(xymoms,L,N)=TMOM(xymoms,L,N) *
     &                (1.-FWASHT(N)-TMFAC(N))
        END DO
      END IF

      END IF
#endif

C**** ADD PRECIPITATION AND LATENT HEAT BELOW
      PRHEAT=CDHEAT(L)+SLH*PRCP
c      if (debug) print*,"cnv5",l,prcp,cond(l)
      PRCP=PRCP+COND(L)
#ifdef SCM
c     save precip by layer for MC (kg/kg)
      if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
          PREMC(L) = PRCP/AIRM(L)
          if (PREMC(L).lt.0.0) premc(l) = 0.0d0
c         write(iu_scm_prt,*) 'L PREMC ',L,PREMC(L)*1000.,PRCP,
c    &        LHPMC(L)
      endif
#endif
#ifdef TRACERS_WATER
      TRPRCP(1:NTX) = TRPRCP(1:NTX) + TRCOND(1:NTX,L)
#ifdef TRDIAG_WETDEPO
      IF (diag_wetdep == 1)
     &     trprcp_mc(l,1:ntx)=trprcp_mc(l,1:ntx)+trcond(1:ntx,l)
#endif
#ifdef TRACERS_SPECIAL_O18
C**** Isotopic equilibration of liquid precip with water vapour
      IF (LHX.eq.LHE .and. PRCP.gt.0) THEN
        DO N=1,NTX
          CALL ISOEQUIL(NTIX(N),TOLD,.TRUE.,QM(L),PRCP,TM(L,N),TRPRCP(N)
     *         ,0.5d0)
        END DO
      END IF
#endif
#endif
C**** end of loop down from top of plume
      END DO EVAP_PRECIP
C****
      IF(PRCP.GT.0.) THEN
        IF(PLE(LMIN)-PLE(LMAX+1).LT.450.) THEN
          CLDMCL(1)=MIN(CLDMCL(1),FMC1)
        ELSE
          RHO=PL(1)/(RGAS*TL(1))
          CLDMCL(1)=MIN(CLDMCL(1)+FMC1*CCM(LMIN)/(RHO*GRAV*WCU(LMIN)*
     *              DTsrc+teeny),FMC1)
        END IF
      END IF
      PRCPMC=PRCPMC+PRCP*FMC1
#ifdef CFMIP
      cfmip_ccp(LMIN) = FMC1*PRCP*BYAM(LMIN)
#endif
#ifdef TRACERS_WATER
      TRPRMC(1:NTX) = TRPRMC(1:NTX) + TRPRCP(1:NTX)*FMC1
#endif
      IF(LMCMIN.GT.LDMIN) LMCMIN=LDMIN
C****
C**** END OF INNER LOOP (2) OVER CLOUD TYPES
C****
      MC1=.FALSE.

#ifdef TRACERS_ON
      call reset_tracer_work_arrays(ldmin,lmax)
#endif

      END DO CLOUD_TYPES

C****
C**** END OF OUTER LOOP (1) OVER CLOUD BASE
C****
      END DO CLOUD_BASE

      IF(LMCMIN.GT.0) THEN

C**** set fssl array
        DO L=1,LM
          FSSL(L)=1.-FMC1
        END DO
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
        IF (diag_wetdep == 1) THEN
          do n=1,ntx
            do l=1,lmcmax
              trcond_mc(l,n)=trcond_mc(l,n)*fmc1
              trdvap_mc(l,n)=trdvap_mc(l,n)*fmc1
              trflcw_mc(l,n)=trflcw_mc(l,n)*fmc1
              trprcp_mc(l,n)=trprcp_mc(l,n)*fmc1
              trnvap_mc(l,n)=trnvap_mc(l,n)*fmc1
              trwash_mc(l,n)=trwash_mc(l,n)*fmc1
            enddo
          enddo
        END IF
#endif

C**** ADJUSTMENT TO CONSERVE CP*T DURING SUBSIDENCE
        SUMAJ=0.
        SUMDP=0.
        DO L=LMCMIN,LMCMAX
          SUMDP=SUMDP+AIRM(L)*FMC1
          SUMAJ=SUMAJ+DGDSM(L)
        END DO
        DO L=LMCMIN,LMCMAX
          DGDSM(L)=DGDSM(L)-SUMAJ*AIRM(L)*FMC1/SUMDP
          SM(L)=SM(L)-SUMAJ*AIRM(L)/(SUMDP*PLK(L))
c          if (debug) print*,"cnv6",l,SUMAJ*AIRM(L)/SUMDP
      END DO

C**** LOAD MASS EXCHANGE ARRAY FOR GWDRAG
#ifdef CUBED_SPHERE
      ! maximum at any level
        AIRXL = MAXVAL(MCFLX(LMCMIN:LMCMAX))
#else
      ! sum over levels
        AIRXL = 0.
        DO L=LMCMIN,LMCMAX
          AIRXL = AIRXL+MCFLX(L)
        END DO
#endif
      END IF

C****
C**** CALCULATE CONVECTIVE CLOUD OPTICAL THICKNESS
C****
      WCONST=WMU*(1.-PEARTH)+WMUL*PEARTH
      WMSUM=0.
#ifdef SCM
      if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
          SCM_LWP_MC = 0.d0
          SCM_IWP_MC = 0.d0
          SCM_WM_MC = 0.d0
      endif
#endif
#ifdef CLD_AER_CDNC
      WMCLWP=0.  ; WMCTWP=0. ; ACDNWM=0. ; ACDNIM=0.
      AREWM=0.   ; AREIM=0.  ; ALWWM=0.  ; ALWIM=0.
      NMCW=0     ; NMCI=0
#endif
      OPTICAL_THICKNESS: DO L=1,LMCMAX
        TL(L)=(SM(L)*BYAM(L))*PLK(L)
        TEMWM=(TAUMCL(L)-SVWMXL(L)*AIRM(L))*1.d2*BYGRAV
        IF(TL(L).GE.TF) WMSUM=WMSUM+TEMWM ! pick up water path
#ifdef SCM
        if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
            SCM_WM_MC(L) = TAUMCL(L)*BYAM(L)-SVWMXL(L)
            if (TL(L).GE.TF) SCM_LWP_MC = SCM_LWP_MC + TEMWM
            if (TL(L).LT.TF) SCM_IWP_MC = SCM_IWP_MC + TEMWM
        endif
#endif



#ifdef CLD_AER_CDNC
        WMCTWP=WMCTWP+TEMWM
        IF(TL(L).GE.TF) WMCLWP=WMCLWP+TEMWM
#endif
C**** DEFAULT OPTICAL THICKNESS = 8 PER 100 MB CLOUD DEPTH, BUT 2 PER
C**** 100 MB INSTEAD FOR DETRAINMENT LEVEL OF SHALLOW/MIDLEVEL
C**** CONVECTION.  ALSO, AN OPTICAL THICKNESS OF 2 PER 100 MB IS
C**** ASSUMED BELOW CLOUD BASE FOR RAIN FALLING FROM DEEP CONVECTION
        IF(CLDMCL(L).GT.0.) THEN
          TAUMCL(L)=AIRM(L)*COETAU
          IF(L.EQ.LMCMAX .AND. PLE(LMCMIN)-PLE(LMCMAX+1).LT.450.)
     *         TAUMCL(L)=AIRM(L)*.02d0
          IF(L.LE.LMCMIN .AND. PLE(LMCMIN)-PLE(LMCMAX+1).GE.450.)
     *         TAUMCL(L)=AIRM(L)*.02d0
        END IF
        SVLAT1(L) = SVLATL(L)   ! used in large-scale clouds
        IF(SVLATL(L).EQ.0.) THEN
          SVLATL(L)=LHE
          IF ( (TPSAV(L).gt.0. .and. TPSAV(L).LT.TF) .or.
     *         (TPSAV(L).eq.0. .and. TL(L).lt.TF) ) SVLATL(L)=LHS
        END IF

C**** FOR DEEP CONVECTIVE ANVILS USE DETRAINED CLOUD WATER AND
C**** EFFECTIVE RADIUS BASED ON FIXED CLOUD DROPLET NUMBER
C**** CONCENTRATION TO CALCULATE OPTICAL THICKNESS (ANALOGOUS TO
C**** CALCULATION FOR STRATIFORM CLOUDS)
        IF(SVWMXL(L).GT.0.) THEN
          FCLD=CLDMCL(L)+1.E-20
          TEM=1.d5*SVWMXL(L)*AIRM(L)*BYGRAV
          WTEM=1.d5*SVWMXL(L)*PL(L)/(FCLD*TL(L)*RGAS)
          IF(SVLATL(L).EQ.LHE.AND.SVWMXL(L)/FCLD.GE.WCONST*1.d-3)
     *         WTEM=1d2*WCONST*PL(L)/(TL(L)*RGAS)
          IF(WTEM.LT.1.d-10) WTEM=1.d-10
C**   Set CDNC for moist conv. clds (const at present)
          MNdO = 59.68d0/(RWCLDOX**3)
          MNdL = 174.d0
          MNdI = 0.06417127d0
#ifdef CLD_AER_CDNC
#ifdef ALT_CDNC_INPUTS
c (min+max)/2 of values from all updrafts at this level.
c could average over all updrafts with appropriate weighting instead
          MNdO=.5*(MNdO_max(L)+MNdO_min(L))
          MNdL=.5*(MNdL_max(L)+MNdL_min(L))
#else
c values from the last updraft computation at its detrainment level
          MNdO=MCDNO1
          MNdL=MCDNL1
#endif
#endif
          MCDNCW=MNdO*(1.-PEARTH)+MNdL*PEARTH
c          if(MCDNCW.gt.0.) write(6,*)"CDNC MC cld",MNdO,MNdL,l
          MCDNCI=MNdI

C**** COMPUTE ANVIL OPTICAL THICKNESS AND DROPLET RADIUS
          CALL ANVIL_OPTICAL_THICKNESS(SVLATL(L),RCLDX,MCDNCW,MCDNCI,
     *      RIMAX,BYBR,FCLD,TEM,WTEM,RCLD,TAUMCL(L))

          RCLDE=RCLD/BYBR       !  effective droplet radius in anvil
#ifdef CLD_AER_CDNC
C** Using the Liu and Daum paramet, Nature, 2002, Oct 10, Vol 419
C** for spectral dispersion effects on droplet size distribution
C** central value of 0.003 for alfa:Rotstayn&Daum, J.Clim,2003,16,21, Nov 2003.
          Repsi=1.d0 - 0.7d0*exp(-0.003d0*MCDNCW)
          Repsis=Repsi*Repsi
          Rbeta=(((1.d0+2.d0*Repsis)**0.667d0))/((1.d0+Repsis)**by3)
          RCLDE=RCLD*Rbeta
          TAUMCL(L)=1.5*TEM/(FCLD*RCLDE+1.E-20) ! evaluate TAUMCL again
          IF(TAUMCL(L).GT.100.) TAUMCL(L)=100.
c     write(6,*)"RCLD",RCLDE,RCLD,Rbeta,WTEM,L,MCDNCW
#endif
          IF (CLDMCL(L).gt.0) CSIZEL(L)=RCLDE  !  effective droplet radius in anvil
#ifdef CLD_AER_CDNC
          if (FCLD.gt.1.d-5.and.SVLATL(L).eq.LHE) then
            ACDNWM(L)= MCDNCW
            AREWM(L) = RCLDE
            ALWWM(L)= WTEM      ! cld water density in g m-3
            NMCW  = NMCW+1
c           if(MCDNCW.gt.0.) write(6,*)"MST CDNC b",ACDNWM(L),MCDNCW,L
          elseif(FCLD.gt.1.d-5.and.SVLATL(L).eq.LHS) then
            ACDNIM(L)= MCDNCI
            AREIM(L) = RCLDE
            ALWIM(L) = WTEM
            NMCI  = NMCI+1
          END IF
#endif
        END IF
        IF(TAUMCL(L).LT.0..and.CLDMCL(L).le.0.) TAUMCL(L)=0.
      END DO OPTICAL_THICKNESS

      IF(LMCMAX.LE.1) THEN
        DO L=1,LM
          IF(PL(L).LT.850.d0) U00L(L)=0.
        END DO
      END IF

      RETURN

#ifdef TRACERS_ON

      CONTAINS

      subroutine reset_tracer_work_arrays(l1,l2)
      integer :: l1,l2
      integer :: l,n
      tpold(1:l2) = 0.
      do n=1,ntx
        do l=l1,l2
          dtm(l,n) = 0.
          dtmom(:,l,n) = 0.
          dtmr(l,n) = 0.
          dtmomr(:,l,n) = 0.
          tmdnl(l,n) = 0.
          tmomdnl(:,l,n) = 0.
        enddo
      enddo
#ifdef TRACERS_WATER
      do l=1,l2
        trcond(1:ntx,l) = 0.
        trcondv(1:ntx,l) = 0.
      enddo
#endif
      end subroutine reset_tracer_work_arrays
#endif

      END SUBROUTINE MSTCNV

c ****************************************************************************************

      SUBROUTINE LSCOND(IERR,WMERR,LERR,i_debug,j_debug)
!@sum  LSCOND column physics of large scale condensation
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
!@calls CTMIX,QSAT,DQSATDT,THBAR
      IMPLICIT NONE

!@var IERR,WMERR,LERR error reporting
      INTEGER, INTENT(OUT) :: IERR,LERR
      REAL*8, INTENT(OUT) :: WMERR
      INTEGER, INTENT(IN) :: i_debug,j_debug
      logical :: debug_out
      REAL*8 LHX

C**** functions
      REAL*8 :: QSAT, DQSATDT,ERMAX
!@var QSAT saturation humidity
!@var DQSATDT dQSAT/dT

!@param CM00 upper limit for autoconversion rate
!@param AIRM0 scaling factor for computing rain evaporation
!@param HEFOLD e-folding length for computing HDEP
!@param COEFM coefficient for ratio of cloud water amount and WCONST
!@param COEFT coefficient used in computing PRATM
!@param COESIG coefficient for equ. 23 of Del Genio et al. (1996)
!@param COEEC coefficient for computing cloud evaporation
!@param ERP exponential power for computing ER
C**** Adjust COEFT and COEFM to change proportion of super-cooled rain
C**** to snow. Increasing COEFT reduces temperature range of super
C**** -cooled rain, increasing COEFM enhances probability of snow.
      REAL*8, PARAMETER :: AIRM0=100.d0, GbyAIRM0=GRAV/AIRM0
      REAL*8 CM00,TEM1,RCLDE1
      REAL*8, PARAMETER :: HEFOLD=500.,COEFM=10.,COEFT=2.5
      REAL*8, PARAMETER :: COESIG=1d-3,COEEC=1000.
      INTEGER, PARAMETER :: ERP=2
      REAL*8, DIMENSION(KMAX) :: UMO1,UMO2,UMN1,UMN2 !@var dummy variables
      REAL*8, DIMENSION(KMAX) :: VMO1,VMO2,VMN1,VMN2 !@var dummy variables
!@var Miscellaneous vertical arrays
      REAL*8, DIMENSION(LM) ::
     *     QSATL,RHF,ATH,SQ,ER,QHEAT,WMPR,
     *     CLEARA,PREP,RH00,EC,WMXM
!@var QSATL saturation water vapor mixing ratio
!@var RHF environmental relative humidity
!@var ATH change in potential temperature
!@var SQ ERMAX dummy variables
!@var ER evaporation of precip
!@var QHEAT change of latent heat
!@var CLEARA fraction of clear region
!@var PREP precip conversion rate
!@var RH00 threshold relative humidity
!@var EC cloud evaporation rate
!@var WMXM cloud water mass (mb)
      REAL*8, DIMENSION(LM+1) :: PREBAR,PREICE
!@var PREBAR,PREICE precip entering layer top for total, snow

#ifdef TRACERS_WATER
!@var TRPRBAR tracer precip entering layer top for total (kg)
      REAL*8, DIMENSION(NTM,LM+1) :: TRPRBAR
!@var DTER change of tracer by evaporation (kg)
!@var FWTOQ fraction of CLW that goes to water vapour
!@var FPR fraction of CLW that precipitates
!@var FER fraction of precipitate that evaporates
      REAL*8 TWMTMP,FWTOQ,FPR,FER
!@var DTPRT tracer-specific change of tracer by precip (kg)
!@var DTERT tracer-specific change of tracer by evaporation (kg)
!@var DTWRT tracer-specific change of tracer by washout (kg)
!@var DTQWT tracer-specific change of tracer by condensation (kg)
!@var FWTOQT tracer-specific fraction of tracer in CLW that evaporates
!@var FQTOWT tracer-specific fraction of gas tracer condensing in CLW
!@var FPRT tracer-specific fraction of tracer in CLW that precipitates
!@var FERT tracer-specific fraction of tracer in precipitate evaporating
!@var FQCONDT fraction of tracer that condenses
      REAL*8, DIMENSION(NTM) ::
     &     FQTOWT,FQCONDT,FERT,FWTOQT,DTERT,DTPRT,DTQWT
      REAL*8 :: DTWRT,FPRT,PRLIQ
!@var BELOW_CLOUD logical- is the current level below cloud?
!@var CLOUD_YET logical- in L loop, has there been any cloud so far?
      LOGICAL BELOW_CLOUD,CLOUD_YET
!@var FWASHT  fraction of tracer scavenged by below-cloud precipitation
      REAL*8 :: FWASHT(NTM),TM_dum(NTM), DTR(NTM)
!@var WMXTR available water mixing ratio for tracer condensation ( )?
!@var b_beta_DT precipitating gridbox fraction from lowest precipitating
!@+   layer. The name was chosen to correspond to Koch et al. p. 23,802.
!@var precip_mm precipitation (mm) from the grid box above for washout
      REAL*8 WMXTR, b_beta_DT, precip_mm
c for tracers in general, added by Koch
      real*8, dimension(ntm) :: THLAW,TR_LEF,TR_LEFT
      REAL*8 THWASH(NTM),TMFAC(NTM),TMFAC2(NTM),CLDSAVT
      INTEGER :: IGAS
!@var TR_LEF limits precurser dissolution following sulfate formation
!@var THLAW Henry's Law determination of amount of tracer dissolution
!@var THWASH Henry's Law for below cloud dissolution
!@var TMFAC,TMFAC2 used to adjust tracer moments
!@var CLDSAVT is present cloud fraction, saved for tracer use
!@var cldprec cloud fraction at lowest precipitating level
      REAL*8 :: cldprec
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
c for sulfur chemistry
!@var WA_VOL Cloud water volume (L). Used by GET_SULFATE.
      REAL*8 WA_VOL
      REAL*8, DIMENSION(NTM) ::SULFOUT,SULFIN,SULFINC
      INTEGER :: IAQCH
#endif
#endif

      REAL*8 AIRMR,BETA,BMAX
     *     ,CBF,CBFC0,CK,CKIJ,CK1,CK2,CKM,CKR,CM,CM0,CM1,DFX,DQ,DQSDT
     *     ,DQSUM,DQUP,DRHDT,DSE,DSEC,DSEDIF,DWDT,DWDT1,DWM,ECRATE,EXPST
     *     ,FCLD,FMASS,FMIX,FPLUME,FPMAX,FQTOW,FRAT,FUNI,dqsum1,fqcond1
     *     ,FUNIL,FUNIO,HCHANG,HDEP,HPHASE,OLDLAT,OLDLHX,PFR,PMI,PML
     *     ,HPBL,PRATIO,QCONV,QHEATC,QLT1,QLT2,QMN1,QMN2,QMO1,QMO2,QNEW
     *     ,QNEWU,QOLD,QOLDU,QSATC,QSATE,RANDNO,RCLDE,RHI,RHN,RHO,RHT1
     *     ,RHW,SEDGE,SIGK,SLH,SMN1,SMN2,SMO1,SMO2,TEM,TEMP,TEVAP,THT1
     *     ,THT2,TLT1,TNEW,TNEWU,TOLD,TOLDU,TOLDUP,VDEF,WCONST,WMN1,WMN2
     *     ,WMNEW,WMO1,WMO2,WMT1,WMT2,WMX1,WTEM,VVEL,RCLD,FCOND
     *     ,PRATM,SMN12,SMO12,QF
       real*8 SNdO,SNdL,SNdI,SCDNCW,SCDNCI
#ifdef CLD_AER_CDNC
!@auth Menon  - storing var for cloud droplet number
       integer, PARAMETER :: SNTM=29 !17+ntm_soa/2
       real*8 Repsis,Repsi,Rbeta,CDNL1,QAUT,DSU(SNTM),QCRIT
     * ,CDNL0,NEWCDN,OLDCDN,SNd
       real*8 dynvis(LM),DSGL(LM,SNTM),DSS(SNTM),r6,r6c
       real*8 NEWCLD,SAVCLD
       real*8, dimension(lm) :: vvel_sv,CLDSAV0
       real*8, dimension(sntm,lm) :: dsu_sv
#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
       real*8 DPP,TEMPR,RHODK,PPRES,PRS        ! for 3 hrly diag
       real*8 D3DL(LM),CWCON(LM)               ! for 3 hrly diag
#endif
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
      integer,PARAMETER         :: mkx=1   ! lm
      real*8,PARAMETER          :: mw0 = 2.094395148947515E-15
      real*8,PARAMETER          :: mi0 = 2.094395148947515E-15
      logical, PARAMETER        :: lSCM=.false.
      logical, PARAMETER        :: wSCM=.false.
      REAL(8), PARAMETER        :: tiny = 1.0D-30
      REAL*8,dimension(mkx)     :: tk0,qk0,pk0,w0,v0,r0,rablk
      REAL*8,dimension(mkx)     :: tk0new,qk0new
      REAL*8,dimension(mkx)     :: ndrop,mdrop,ncrys,mcrys
      REAL*8,dimension(mkx)     :: ndrop_old,mdrop_old
      REAL*8,dimension(mkx)     :: ndrop_new,mdrop_new
      REAL*8,dimension(mkx)     :: ndrop_blk,mdrop_blk
      REAL*8,dimension(mkx)     :: ndrop_res,mdrop_res
      REAL*8,dimension(mkx)     :: ncrys_old,mcrys_old
      REAL*8,dimension(mkx)     :: ncrys_new,mcrys_new
      REAL*8,dimension(mkx)     :: ncrys_blk,mcrys_blk
      REAL*8,dimension(mkx)     :: ncrys_res,mcrys_res
      REAL*8,dimension(mkx)     :: npccn,nprc,nnuccc,nnucci
      REAL*8,dimension(mkx)     :: mpccn,mprc,mnuccc,mnucci
      REAL*8,dimension(mkx)     :: nnuccd,nnucmd,nnucmt
      REAL*8,dimension(mkx)     :: mnuccd,mnucmd,mnucmt
      REAL*8,dimension(mkx)     :: nc_tnd,qc_tnd,ni_tnd,qi_tnd
      REAL*8,dimension(mkx)     :: nc_tot,qc_tot,ni_tot,qi_tot

      REAL*8                    :: DTB2M,QAUT_B2M
      real*8                    :: NEWCDNC,OLDCDNC
      real*8,dimension(LM)      :: DCLD
      logical                   :: ldummy=.false.
      character*8               :: sname='lscond: '
      integer                   :: nm,iuo=801
#ifdef TRACERS_AMP
      real*8                    :: naero (mkx,nmodes)
c     real*8,dimension(lm,nmodes)   :: nactc
#endif
#endif

!@var BETA,BMAX,CBFC0,CKIJ,CK1,CK2,PRATM dummy variabls
!@var SMN12,SMO12 dummy variables
!@var AIRMR
!@var CBF enhancing factor for precip conversion
!@var CK ratio of cloud top jumps in moist static energy and total water
!@var CKM CTEI threshold in MacVean and Mason theory
!@var CKR CTEI threshold in Randall theory
!@var CM conversion rate for large cloud water content
!@var CM0,CM1 limiting autoconversion rate for large cloud water content
!@var DFX iteration increment
!@var DQ condensed water vapor
!@var DQSDT derivative of saturation vapor pressure w.r.t. temperature
!@var DQSUM,DQUP dummy variables
!@var DRHDT time change of relative humidity
!@var DSE moist static energy jump at cloud top
!@var DSEC critical DSE for CTEI to operate
!@var DSFDIF DSE-DSEC
!@var DWDT,DWDT1 time change rates of cloud water
!@var ECRATE cloud droplet evaporation rate
!@var EXPST exponential term in determining the fraction in CTEI
!@var FCLD cloud fraction
!@var FMIX,FRAT fraction of mixing in CTEI
!@var FMASS mass of mixing in CTEI
!@var FPLUME fraction of mixing in CTEI
!@var FPMAX max fraction of mixing in CTEI
!@var FQTOW fraction of water vapour that goes to CLW
!@var FUNI the probablity for ice cloud to form
!@var FUNIL,FUNIO FUNI over land, ocean
!@var HCHANG,HPHASE latent heats for changing phase
!@var HDEP,HPBL layer depth (m)  (note change of unit km-->m)
!@var OLDLAT,OLDLHX previous LHX
!@var PFR PROBABLITY OF GLACIATION OF SUPER-COOLED WATER
!@var PMI icy precip entering the layer top
!@var PML layer's cloud water devided by GRAV
!@var PRATIO PMI/PML
!@var QCONV convergence of latent heat
!@var QHEATC,QLT1,QLT2,QMN1,QMN2,QMO1,QMO2 dummy variables
!@var QNEW,QNEWU updated specific humidity
!@var QOLD,QOLDU previous specific humidity
!@var QSATC saturation vapor mixing ratio
!@var QSATE saturation vapor mixing ratio w.r.t. water
!@var RANDNO random number
!@var RCLD,RCLDE cloud particle's radius, effective radius
!@var RHI relative humidity w.r.t. ice
!@var RHN dummy variable
!@var RHO air density
!@var RHT1,RHW,SIGK dummy variables
!@var SEDGE potential temperature at layer edge
!@var SMN1,SMN2,SMO1,SMO2,TEM,TEMP,TEVAP,THT1,THT2,TLT1 dummy variables
!@var TNEW,TNEWU updated tempertures
!@var TOLD,TOLDU,TOLDUP previous temperatures
!@var VDEF = VVEL - VSUB
!@var WCONST,WMN1,WMN2,WMO1,WMO2,WMT1,WMT2,WMX1 dummy variables
!@var WMNEW updated cloud water mixing ratio
!@var WTEM cloud water density (g m**-3)
!@var VVEL vertical velocity (cm/s)
!@var FCOND QF dummy variables
      INTEGER LN,ITER !@var LN,ITER loop variables
      LOGICAL BANDF  !@var BANDF true if Bergero-Findeisen proc. occurs
      LOGICAL FORM_CLOUDS !@var FORM_CLOUDS true if clouds are formed

      INTEGER K,L,N  !@var K,L,N loop variables

      REAL*8 THBAR !@var THBAR potential temperature at layer edge

C****
C**** LARGE-SCALE CLOUDS AND PRECIPITATION
C**** THE LIQUID WATER CONTENT IS PREDICTED
C****
      IERR=0
C****
      debug_out=.false.
      PRCPSS=0.
      HCNDSS=0.
      CKIJ=1.
      RCLDX=radius_multiplier
      RTEMP=funio_denominator
      CMX=autoconv_multiplier
      WMUIX=wmui_multiplier
C**** initialise vertical arrays
      ER=0.
      EC=0.
      PREP=0.
      PREBAR=0.
      LHP=0.
#ifdef SCM
      if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
          LHPSAV=0.
          PRESAV=0.
      endif
#endif
      QHEAT=0.
      CLDSSL=0
      TAUSSL=0
      WMPR=0.
#ifdef TRACERS_WATER
      TRPRSS = 0.
      TRPRBAR = 0.
      BELOW_CLOUD=.false.
      CLOUD_YET=.false.
      CLDSAVT = 0.
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      DT_SULF_SS(1:NTM,:)=0.
#endif
      TR_LEF(:)=1. ! currently, only aqchem_list elements vary with layer
      thlaw(:) = 0.  ! nonzero only for gas tracers
      thwash(:) = 0. ! nonzero only for gas tracers
      tmfac(:) = 0.  ! nonzero only for gas tracers
      tmfac2(:) = 0. ! nonzero only for gas tracers
      fwasht(:) = 0. ! nonzero only for aerosols
      fqcondt(:) = 0.
      fqtowt(:) = 0.
#endif
#ifdef CLD_AER_CDNC
       CDN3DL=0.
       CRE3DL=0.
       SMLWP=0.
c      AERTAU=0.
       DSGL(:,1:SNTM)=0.
#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
       CTEML=0.
       CD3DL=0.
       CL3DL=0.
       CI3DL=0.
#endif
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
       WMXICE(:)=0.
c      print *,sname,'WMX, WMXICE = ', WMX, WMXICE
#endif
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
      IF (diag_wetdep == 1) THEN
C**** initialize diagnostic arrays
        trwash_ls=0.D0
        trprcp_ls=0.D0
        trclwc_ls=0.D0
        trclwe_ls=0.D0
        trcond_ls=0.D0
      END IF
#endif
      DO L=1,LP50
        CLEARA(L)=1.-CLDSAVL(L)
        IF(WMX(L).LE.0.) CLEARA(L)=1.
#ifdef CLD_AER_CDNC
        CLDSAV0(L) = 1.-CLEARA(L)
#endif
      END DO
      DQUP=0.
      TOLDUP=TL(LP50)
      PREICE(LP50+1)=0.
      WCONST=WMU*(1.-PEARTH)+WMUL*PEARTH
         SSHR=0.
         DCTEI=0.
#ifdef CFMIP
      cfmip_reff(L) = 0.
#endif
C****
C**** MAIN L LOOP FOR LARGE-SCALE CONDENSATION, PRECIPITATION AND CLOUDS
C****
      DO L=LP50,1,-1
      TOLD=TL(L)
      QOLD=QL(L)
      OLDLHX=SVLHXL(L)
      OLDLAT=SVLATL(L)
C**** COMPUTE VERTICAL VELOCITY IN CM/S
      TEMP=100.*RGAS*TL(L)/(PL(L)*GRAV)
      IF(L.EQ.1)  THEN
         VVEL=-SDL(L+1)*TEMP
      ELSE IF(L.EQ.LM)  THEN
         VVEL=-SDL(L)*TEMP
      ELSE
         VVEL=-.5*(SDL(L)+SDL(L+1))*TEMP
      END IF
      VDEF=VVEL-VSUBL(L)
#ifdef CLD_AER_CDNC
      vvel_sv(l) = vvel ! save for opt. depth calc.
#endif

      FCLD=(1.-CLEARA(L))*FSSL(L)+teeny

C**** COMPUTE THE PROBABILITY OF ICE FORMATION, FUNI, AND
C**** THE PROBABLITY OF GLACIATION OF SUPER-COOLED WATER, PFR
C**** DETERMINE THE PHASE MOISTURE CONDENSES TO
C**** DETERMINE THE POSSIBILITY OF B-F PROCESS
      BANDF=.FALSE.
      LHX=LHE
      CBF=1. + EXP(-((TL(L)-258.16d0)/10.)**2)

      IF (TL(L).LE.238.16) THEN     ! below -35C: force ice
        LHX=LHS
      ELSEIF (TL(L).GE.TF) THEN ! above freezing: force water
        LHX=LHE
      ELSE                      ! in between: compute probability
        IF(TL(L).GT.269.16) THEN ! OC/SI/LI clouds: water above -4
          FUNIO=0.
        ELSE
          FUNIO=1.-EXP(-((TL(L)-269.16d0)/RTEMP)**4)
        END IF
        IF(TL(L).GT.263.16) THEN ! land clouds water: above -10
          FUNIL=0.
        ELSE
          FUNIL=1.-EXP(-((TL(L)-263.16d0)/RTEMP)**4)
        END IF
        FUNI=FUNIO*(1.-PEARTH)+FUNIL*PEARTH
        RANDNO=RNDSSL(1,L)       !  RANDNO=RANDU(XY)
        IF(RANDNO.LT.FUNI) LHX=LHS

        IF (OLDLHX.EQ.LHS.AND.TL(L).LT.TF) LHX=LHS   ! keep old phase
        IF (OLDLHX.EQ.LHE.AND.TL(L).GT.269.16d0) LHX=LHE ! keep old phase
C**** special case 1) if ice previously then stay as ice (if T<Tf)
C       IF((OLDLHX.EQ.LHS.OR.OLDLAT.EQ.LHS).AND.TL(L).LT.TF) THEN
        IF(OLDLAT.EQ.LHS.AND.TL(L).LT.TF.and.SVLAT1(L).gt.0.) THEN
          IF(LHX.EQ.LHE) BANDF=.TRUE.
          LHX=LHS
        END IF
c        if (debug) print*,"ls0",l,oldlhx,oldlat,lhx,lhp(l)

        IF (L.LT.LP50) THEN
C**** Decide whether precip initiates B-F process
          PML=WMX(L)*AIRM(L)*BYGRAV
          PMI=PREICE(L+1)*DTsrc
          RANDNO=RNDSSL(2,L)     !  RANDNO=RANDU(XY)
C**** Calculate probability of ice precip seeding a water cloud
          IF (LHX.EQ.LHE.AND.PMI.gt.0) THEN
            PRATIO=MIN(PMI/(PML+1.E-20),10d0)
            CM00=3.d-5           ! reduced by a factor of 3
            IF(ROICE.GT..1d0) CM00=3.d-4
            CM0=CM00
            IF(VDEF.GT.0.) CM0=CM00*10.**(-0.2*VDEF)
            CBFC0=.5*CM0*CBF*DTsrc
            PFR=(1.-EXP(-(PRATIO*PRATIO)))*(1.-EXP(-(CBFC0*CBFC0)))
            IF(PFR.GT.RANDNO) THEN
              BANDF=.TRUE.
              LHX=LHS
            END IF
          END IF
C**** If liquid rain falls into an ice cloud, B-F must occur
          IF (LHP(L+1).EQ.LHE .AND. LHX.EQ.LHS .AND. PML.GT.0.)
     *         BANDF=.TRUE.
        END IF
      END IF
      IF(LHX.EQ.LHS .AND. (OLDLHX.EQ.LHE.OR.OLDLAT.EQ.LHE)) BANDF=.TRUE.

C**** COMPUTE THE LIMITING AUTOCONVERSION RATE FOR CLOUD WATER CONTENT
      CM00=1.d-4
      IF(LHX.EQ.LHS.AND.SVWMXL(L).LE.0d0) CM00=1.d-3
      IF(LHX.EQ.LHE) THEN                 ! reduced by a factor of 3
        CM00=3.d-5
        IF(ROICE.GT..1d0) CM00=3.d-4
      END IF
      CM0=CM00
      IF(VDEF.GT.0.) CM0=CM00*10.**(-0.2*VDEF)

C**** COMPUTE RELATIVE HUMIDITY
      QSATL(L)=QSAT(TL(L),LHX,PL(L))
      RH1(L)=QL(L)/QSATL(L)
      IF(LHX.EQ.LHS.AND.WMX(L).LE.0d0) THEN          ! Karcher and Lohmann formula
        QSATE=QSAT(TL(L),LHE,PL(L))
        RHW=(2.583d0-TL(L)/207.83)*(QSAT(TL(L),LHS,PL(L))/QSATE)
        IF(TL(L).LT.238.16) RH1(L)=QL(L)/(QSATE*RHW)
      END IF
C**** PHASE CHANGE OF CLOUD WATER CONTENT
      HCHANG=0.
      IF(OLDLHX.EQ.LHE.AND.LHX.EQ.LHS) HCHANG= WML(L)*LHM
      IF(OLDLHX.EQ.LHS.AND.LHX.EQ.LHE) HCHANG=-WML(L)*LHM
      IF(OLDLAT.EQ.LHE.AND.LHX.EQ.LHS) HCHANG=HCHANG+SVWMXL(L)*LHM
      IF(OLDLAT.EQ.LHS.AND.LHX.EQ.LHE) HCHANG=HCHANG-SVWMXL(L)*LHM
c        if (debug) print*,"ls1",l,hchang

      SVLHXL(L)=LHX
      TL(L)=TL(L)+HCHANG/(SHA*FSSL(L)+teeny)
      TH(L)=TL(L)/PLK(L)
#ifdef SCM
c     preserving T for difference from before updating with ARM data
      if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
          if (NRINIT.ne.0) then
             TH(L) = SCM_SAVE_T(L)*PLK(L)+SCM_DEL_T(L)+
     &            HCHANG/(SHA*FSSL(L)+teeny)
             TH(L) = TH(L)/PLK(L)
          endif
      endif
#endif

      ATH(L)=(TH(L)-TTOLDL(L))*BYDTsrc
C**** COMPUTE RH IN THE CLOUD-FREE AREA, RHF
      RHI=QL(L)/QSAT(TL(L),LHS,PL(L))
    ! this formulation is used for consistency with current practice
      RH00(L)=U00a
      IF(PL(L).LT.850.d0) THEN

        RH00(L) = RH00(L)/(RH00(L) + (1.-RH00(L))*AIRM(L)/35.)

        IF(VDEF.GT..2d0.AND.LMCMAX.LE.1) RH00(L)=
     *    RH00(L)*MIN(SQRT(.2d0/VDEF),.5d0) ! dependece on vertical velocity
      END IF
      IF(U00L(L).GT.RH00(L)) RH00(L)=U00L(L)
C**** Option to treat boundary layer differently
      IF (do_blU00.eq.1) then
        IF (L.LE.DCL) THEN      ! boundary layer clouds
C**** calculate total pbl depth
          HPBL=0.
          DO LN=1,DCL
            HPBL=HPBL+AIRM(LN)*TL(LN)*RGAS/(GRAV*PL(LN))
          END DO
C**** Scale HPBL by HRMAX to provide tuning control for PBL clouds
          HDEP = MIN(HPBL,HRMAX*(1.-EXP(-HPBL/HEFOLD)))
C**** Special conditions for boundary layer contained wholly in layer 1
          IF (DCL.LE.1) THEN
            IF (RIS.GT.1.) HDEP=10d0
            IF (RIS.LE.1..AND.RI1.GT.1.) HDEP=50d0
            IF (RIS.LE.1..AND.RI1.LE.1..AND.RI2.GT.1.) HDEP=100d0
          END IF
C**** Estimate critical rel. hum. based on parcel lifting argument
          RH00(L)=1.-GAMD*LHE*HDEP/(RVAP*TS*TS)
          IF(RH00(L).LT.0.) RH00(L)=0.
        END IF
      END IF
C****
      IF(RH00(L).LT.0.) RH00(L)=0.
      IF(RH00(L).GT.1.) RH00(L)=1.
      RHF(L)=RH00(L)+(1.-CLEARA(L))*(1.-RH00(L))
C**** Set precip phase to be the same as the cloud, unless precip above
C**** is ice and temperatures after ice melt would still be below TFrez
      LHP(L)=LHX
      IF (LHP(L+1).eq.LHS .and.
     *     TL(L).lt.TF+DTsrc*LHM*PREICE(L+1)*GRAV*BYAM(L)*BYSHA)
     *     LHP(L)=LHP(L+1)
#ifdef SCM
c     save precip phase LS
      if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
          LHPSAV(L) = LHP(L)
c         write(iu_scm_prt,*) 'L LHPSAV ',L,LHPSAV(L)
      endif
#endif

#if (defined CLD_AER_CDNC) && (defined TRACERS_AEROSOLS_Koch)
!@auth Menon  saving aerosols mass for CDNC prediction
      DO N=1,SNTM
        DSS(N)=1.d-10
        DSGL(L,N)=1.d-10
      END DO
      DO N=1,NTX
        select case (trname(ntix(n)))
        case('SO4')
          DSGL(L,1)=tm(l,n)     !n=4
          DSS(1) = DSGL(L,1)
        case('seasalt1')
          DSGL(L,2)=tm(l,n)     !n=6
          DSS(2) = DSGL(L,2)
        case('seasalt2')
          DSGL(L,3)=tm(l,n)     !n=7
          DSS(3) = DSGL(L,3)
        case('OCIA')
          DSGL(L,4)=tm(l,n)     !n=12
          DSS(4) = DSGL(L,4)
        case('OCB')
          DSGL(L,5)=tm(l,n)     !n=13
          DSS(5) = DSGL(L,5)
        case('BCIA')
          DSGL(L,6)=tm(l,n)     !n=9
          DSS(6) = DSGL(L,6)
        case('BCB')
          DSGL(L,7)=tm(l,n)     !n=10
          DSS(7) = DSGL(L,7)
        case('OCII')
          DSGL(L,8)=tm(l,n)     !n=11
          DSS(8) = DSGL(L,8)
        case('BCII')
          DSGL(L,9)=tm(l,n)     !n=8
          DSS(9) = DSGL(L,9)
#ifdef TRACERS_DUST
        case('Clay')
          DSGL(L,10)=tm(l,n)    !n=23
          DSS(10) = DSGL(L,10)
        case('Silt1')
          DSGL(L,11)=tm(l,n)    !n=23
          DSS(11) = DSGL(L,11)
        case('Silt2')
          DSGL(L,12)=tm(l,n)    !n=23
          DSS(12) = DSGL(L,12)
        case('Silt3')
          DSGL(L,13)=tm(l,n)    !n=23
          DSS(13) = DSGL(L,13)
#endif
#ifdef TRACERS_NITRATE
        case('NO3p')
          DSGL(L,14)=tm(l,n)    !n=23
          DSS(14) = DSGL(L,14)
#endif
#ifdef TRACERS_HETCHEM
C**** Here are dust particles coated with sulfate
        case('SO4_d1')
          DSGL(L,15)=tm(l,n)    !n=20
          DSS(15) = DSGL(L,15)
        case('SO4_d2')
          DSGL(L,16)=tm(l,n)    !n=21
          DSS(16) = DSGL(L,16)
        case('SO4_d3')
          DSGL(L,17)=tm(l,n)    !n=22
          DSS(17) = DSGL(L,17)
#endif
#ifdef TRACERS_AEROSOLS_SOA
        case('isopp1a')
          DSGL(L,18)=tm(l,n)
          DSS(18) = DSGL(L,18)
        case('isopp2a')
          DSGL(L,19)=tm(l,n)
          DSS(19) = DSGL(L,19)
#ifdef TRACERS_TERP
        case('apinp1a')
          DSGL(L,20)=tm(l,n)
          DSS(20) = DSGL(L,20)
        case('apinp2a')
          DSGL(L,21)=tm(l,n)
          DSS(21) = DSGL(L,21)
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
        end select
      END DO      !end of n loop for tracers
#endif   /* tracerpart and cld-aer part */
C***Setting constant values of CDNC over land and ocean to get RCLD=f(CDNC,LWC)
      SNdO = 59.68d0/(RWCLDOX**3)
      SNdL = 174.d0
      SNdI = 0.06417127d0
      SCDNCW=SNdO*(1.-PEARTH)+SNdL*PEARTH
      SCDNCI=SNdI
      WMUI=WMUIX*.001         ! .0001
      WMUSI=0.1
#if (defined CLD_AER_CDNC) && (defined TRACERS_AEROSOLS_Koch)
      CALL GET_CDNC(L,LHX,WCONST,WMUI,AIRM(L),WMX(L),DXYPIJ,
     *FCLD,CLEARA(L),CLDSAVL(L),DSS,PL(L),TL(L),
     *OLDCDL(L),VVEL,SME(L),DSU,CDNL0,CDNL1)
      DSU_SV(:,L) = DSU(:) ! save for opt. depth calc.
c     write(6,*)"Where is",DSU(L),l
      SNd=CDNL1
cC** Pass old and new cloud droplet number
      NEWCDN=SNd
      OLDCDN=CDNL0
c     if(SNd.gt.20.) write(6,*)"SM 11 CDNC",NEWCDN,OLDCDN,L
#endif
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
C Microphysical time step
       dtB2M=DTsrc
C Set all tendencies to zero
       ldummy=execute_bulk2m_driver('all','2zero',mkx)
C Set thermodynamics
       tk0=TL(L)                 ! Temperature, [K]
       qk0=QL(L)                 ! Water vapor mixing ratio, [kq/kq]
       pk0=PL(L)                 ! Pressure, [hPa]
       w0=VVEL*1.d-02            ! Large-scale velocity, [m/s]
       v0=WTURB(L) !; v0=w3d(k,i,j) ! Sub-grid velocity, [m/s]
       r0=0.0  !RDTNDL(L)        ! tendency due to radiation, [K/s], not needed
c      print*,"rad tendency",w0,v0,r0,L
       ldummy=execute_bulk2m_driver('all'
     *           ,tk0,qk0,pk0,w0,v0,r0)
c Set microphysics
       IF(LHX.EQ.LHE)  THEN
          mdrop=WMX(L)            ! drop content, [kg water/kg air]
          ndrop =OLDCDL(L)*1.d6   !convert from cm-3 to m-3
          if (WMX(L).eq.0.) ndrop=0.d0
          ncrys=0.d0;mcrys=0.0d0
        ELSE
          WMXICE(L) = WMX(L)
          mcrys=WMXICE(L)         ! crys content, [kg water/kg air]
          ncrys=OLDCDI(L)*1.d6     ! convert cm-3 to m-3; set at 0.1 l-1 = 1.d-4 cm-3
          if (WMX(L).eq.0.) ncrys=0.d0
          ndrop=0.0d0;mdrop=0.0d0
        ENDIF
c
        ndrop_old=ndrop;mdrop_old=mdrop;ncrys_old=ncrys;mcrys_old=mcrys
        ndrop_new=0.0d0;mdrop_new=0.0d0;ncrys_new=0.0d0;mcrys_new=0.0d0
        nc_tnd=0.0d0;qc_tnd=0.0d0;ni_tnd=0.0d0;qi_tnd=0.0d0
        nc_tot=0.0d0;qc_tot=0.0d0;ni_tot=0.0d0;qi_tot=0.0d0
c
C** Convert from l-1 to cm-3====>  1 l^-1 = 10^-3 cm^-3 = 10^3 m^-3
c      ldummy=execute_bulk2m_driver('all'
c    *           ,ndrop,mdrop,ncrys,mcrys,'end')
#ifdef TRACERS_AMP
        do nm=1,nmodes
          naero(mkx,nm)=nactc(l,nm)
c       if(l.eq.1) then
c       if(nactc(l,nm).gt.1.)write(6,*)"Callmatrix",naero(mkx,nm)*1.e-6
c        endif
        enddo
       ldummy=execute_bulk2m_driver('all'
     *           ,ndrop,mdrop,ncrys,mcrys,naero,nmodes,'end')
#endif
#ifdef TRACERS_AEROSOLS_Koch
       ldummy=execute_bulk2m_driver('all'
     *           ,ndrop,mdrop,ncrys,mcrys,'end')
#endif
c Make calls to calculate hydrometeors' growth rates due to microphysical
c processes
c content      :       [kg water/kg air/s]
c concentration:       [No/kg air/s]
c Activation of cloud droplets: prescribed AP spectrum
C*** For the originial HM scheme with fixed distributions for amm. sulfate
c       ldummy=execute_bulk2m_driver('hugh','drop_nucl',dtB2M,mkx)

C*** Use this if using the Lohmann or Gultepe scheme  for mass to number
#ifdef TRACERS_AEROSOLS_Koch
        OLDCDNC=OLDCDN*1.d6  !convert from cm-3 to m-3
        NEWCDNC=NEWCDN*1.d6  !convert from cm-3 to m-3
        ldummy=execute_bulk2m_driver('gult','drop_nucl',dtB2M,mkx,
     &OLDCDNC,NEWCDNC)
#endif
#ifdef TRACERS_AMP
C*** Using the AMP_actv interface from MATRIX
        ldummy=execute_bulk2m_driver('matr','drop_nucl',dtB2M,mkx)
#endif
c Droplets' autoconversion: Beheng (concentration and content)
        ldummy=execute_bulk2m_driver('hugh','drop_auto',dtB2M,mkx)
c Freezing of cloud droplets (contact and immersion)
        ldummy=execute_bulk2m_driver('hugh','drop_frzn',dtB2M,mkx)
c Crystal nucleation: Ncrys=anuc(k)*wef**bnuc(k)
        ldummy=execute_bulk2m_driver('hugh','crys_nucl',dtB2M,mkx)
c Numerous processes of water-water,water-ice, ice-ice interaction,
c condensation/evaporation/deposition/sublimation, ice multiplication
c and sedimention are ready to be called. There are only a few examples:
c        ldummy=execute_bulk2m_driver('hugh','drop_rain',dtB2M,mkx)
c        ldummy=execute_bulk2m_driver('hugh','drop_snow',dtB2M,mkx)
c        ldummy=execute_bulk2m_driver('hugh','crys_auto',dtB2M,mkx)
c        ldummy=execute_bulk2m_driver('hugh','crys_snow',dtB2M,mkx)
c        ldummy=execute_bulk2m_driver('hugh','crys_cond',dtB2M,mkx)
c        ldummy=execute_bulk2m_driver('hugh','snow_melt',dtB2M,mkx)
c
c In this chain of events the very last call, which applies saturation
c adjustment to keep environment about water saturation, is supposed to be
c        ldummy=execute_bulk2m_driver('hugh','drop_cond',dtB2M,mkx)
c
c To ensure calculated growth rates don't lead to negative contents:
c Previous call ('drop_cond') HAS to be used.
c        ldummy=execute_bulk2m_driver('all','make_balance',mkx)
c Otherwise, the caller is responsible to handle the problem.
c
c To get particular growth rate:
c        rablk = execute_bulk2m_driver('get','rate','mprc')
c return value "rablk" is real*8 array whose dimension is equal to "mkx"
c
c To get tendencies of temperature, water vapor mixing ratio or
c tendencies of concentration/content of particular hydrometeor:
c        rablk = execute_bulk2m_driver('get','tnd','qc_tnd')
c return value "rablk" is real*8 array whose dimension is equal to "mkx"
c
c To get parameters of hydrometeors size distributions:
c        rablk = execute_bulk2m_driver('get','val','ec')
c return value "rablk" is real*8 array whose dimension is equal to "mkx"
c
c To update concentration and contents due to uncommented processes
c If "make_balance" was used:
c        ndrop=ndrop+dtB2M*execute_bulk2m_driver('get','tnd','nc_tnd')
c        mdrop=mdrop+dtB2M*execute_bulk2m_driver('get','tnd','qc_tnd')
c        ncrys=ncrys+dtB2M*execute_bulk2m_driver('get','tnd','ni_tnd')
c        mcrys=mcrys+dtB2M*execute_bulk2m_driver('get','tnd','qi_tnd')

c
c This is our case
c Without "make_balance":
c Droplet concentration
c       if(l.eq.1) then
c        npccn=execute_bulk2m_driver('get','npccn')
c        nprc=execute_bulk2m_driver('get','nprc')
c        nnuccc=execute_bulk2m_driver('get','nnuccc')
c        nnucci=execute_bulk2m_driver('get','nnucci')
c         if(npccn(l).gt.1.)write(6,*)"check BLK",ndrop*1.e-6,
c    *npccn*1.e-6,nprc*1.e-6,nnuccc*1.e-6,nnucci*1.e-6
c       endif

        ndrop=ndrop+(
c       npccn              ! change n droplets activation
     * +execute_bulk2m_driver('get','npccn')
c       nprc               ! change n autoconversion of droplets:
     * -execute_bulk2m_driver('get','nprc')
c       nnuccc             ! change n due to con droplets freez
     * -execute_bulk2m_driver('get','nnuccc')
c       nnucci             ! change n due to imm droplets freez
     * -execute_bulk2m_driver('get','nnucci')
c
     *              )*dtB2M
c
c Droplet content
        mdrop=mdrop+(
c       mpccn              ! change q droplets activation
     * +execute_bulk2m_driver('get','mpccn')
c       mprc               ! change q autoconversion of droplets:
     * -execute_bulk2m_driver('get','mprc')
c       mnuccc             ! change q due to con droplets freez
     * -execute_bulk2m_driver('get','mnuccc')
c       mnucci             ! change q due to imm droplets freez
     * -execute_bulk2m_driver('get','mnucci')
c
     *              )*dtB2M
c
c Crystal concentration
        ncrys=ncrys+(
c       nnuccc             ! change n due to contact droplets freez
     * +execute_bulk2m_driver('get','nnuccc')
c       nnucci             ! change n due to immersion droplets freez
     * +execute_bulk2m_driver('get','nnucci')
c       nnuccd             ! change n freezing aerosol (prim ice nuc)
     * +execute_bulk2m_driver('get','nnuccd')
     *              )*dtB2M
c      nnucmd              ! change n cond freezing Meyer's (prim ice nuc)
     * +execute_bulk2m_driver('get','nnucmd')
c      nnucmt              ! change n cont freezing Meyer's (prim ice nuc)
     * +execute_bulk2m_driver('get','nnucmt')
c
c Crystal content
        mcrys=mcrys+(
c       mnuccc             ! change q due to con droplets freez
     * +execute_bulk2m_driver('get','mnuccc')
c       mnucci             ! change q due to imm droplets freez
     * +execute_bulk2m_driver('get','mnucci')
c       mnuccd             ! change q freezing aerosol (prim ice nuc)
     * +execute_bulk2m_driver('get','mnuccd')
     *              )*dtB2M
c      mnucmd              ! change q cond freezing Meyer's (prim ice nuc)
     * +execute_bulk2m_driver('get','mnucmd')
c      mnucmt              ! change q cont freezing Meyer's (prim ice nuc)
     * +execute_bulk2m_driver('get','mnucmt')
c
        if_balance:
     *  if( (ndrop(mkx) .lt. 0) .or. (mdrop(mkx) .lt. 0)  .or.
     *      (ncrys(mkx) .lt. 0) .or. (mcrys(mkx) .lt. 0)) then
c
        if(lSCM) then
         write(6,*)"stop BLK: ndrop_old,mdrop_old,ncrys_old,mcrys_old"
     *,l,ndrop_old*1.e-6,mdrop_old*1.e+3,ncrys_old*1.e-3,mcrys_old*1.e+3
c
         write(6,*)"stop BLK: ndrop,mdrop,ncrys,mcrys"
     *,l,ndrop*1.e-6,mdrop*1.e+3,ncrys*1.e-3,mcrys*1.e+3
         endif
c
c No/m^3
c
         npccn  =             ! change n droplets activation
     * +execute_bulk2m_driver('get','npccn')*dtB2M
         nprc   =             ! change n autoconversion of droplets:
     * -execute_bulk2m_driver('get','nprc')*dtB2M
         nnuccc =             ! change n due to con droplets freez
     * -execute_bulk2m_driver('get','nnuccc')*dtB2M
         nnucci =             ! change n due to imm droplets freez
     * -execute_bulk2m_driver('get','nnucci')*dtB2M

         nc_tot = npccn + nprc + nnuccc + nnucci
c
c No/cc
c
        if(lSCM) then
         write(6,*)"stop BLK: ndrop_old,nc_tot,ndrop"
     *,l,ndrop_old*1.e-6,nc_tot*1.e-6,ndrop*1.e-6
         write(6,*)"stop BLK: npccn,nprc,nnuccc,nnucci"
     *,l,npccn*1.e-6,nprc*1.e-6,nnuccc*1.e-6,nnucci*1.e-6
         endif
c
c kg/kg
c
         mpccn  =              ! change q droplets activation
     *    execute_bulk2m_driver('get','mpccn')*dtB2M
         mprc   =              ! change q autoconversion of droplets:
     *   -execute_bulk2m_driver('get','mprc')*dtB2M
         mnuccc =              ! change q due to con droplets freez
     *   -execute_bulk2m_driver('get','mnuccc')*dtB2M
         mnucci =              ! change q due to imm droplets freez
     *   -execute_bulk2m_driver('get','mnucci')*dtB2M

         qc_tot = mpccn + mprc + mnuccc + mnucci
c
c g/kg
c
        if(lSCM) then
         write(6,*)"stop BLK: mdrop_old,qc_tot,mdrop"
     *,l,mdrop_old*1.e+3,qc_tot*1.e+3,mdrop*1.e+3
         write(6,*)"stop BLK: mpccn,mprc,mnuccc,mnucci"
     *,l,mpccn*1.e+3,mprc*1.e+3,mnuccc*1.e+3,mnucci*1.e+3
         endif
c
c No/m^3
c
         nnuccc  =             ! change n due to contact droplets freez
     *   +execute_bulk2m_driver('get','nnuccc')*dtB2M
         nnucci  =             ! change n due to immersion droplets freez
     *   +execute_bulk2m_driver('get','nnucci')*dtB2M
         nnuccd  =             ! change n freezing aerosol (prim ice nuc)
     *   +execute_bulk2m_driver('get','nnuccd') ! *dtB2M
         nnucmd =              ! change n cond freezing Meyer's (prim ice nuc)
     *   +execute_bulk2m_driver('get','nnucmd') ! *dtB2M
         nnucmt =              ! change n cont freezing Meyer's (prim ice nuc)
     *   +execute_bulk2m_driver('get','nnucmt') ! *dtB2M

         ni_tot = nnuccc + nnucci + nnuccd + nnucmd + nnucmt
c
c No/l
c
        if(lSCM) then
         write(6,*)"stop BLK: ncrys_old,ni_tot,ncrys"
     *,l,ncrys_old*1.e-3,ni_tot*1.e-3,ncrys*1.e-3
         write(6,*)"stop BLK: nnuccc,nnucci,nnuccd,nnucmd,nnucmt"
     *,l,nnuccc*1.e-3,nnucci*1.e-3,nnuccd*1.e-3,nnucmd*1.e-3
     *,nnucmt*1.e-3
         endif
c
c kg/kg
c
         mnuccc  =             ! change q due to contact droplets freez
     *   +execute_bulk2m_driver('get','mnuccc')*dtB2M
         mnucci  =             ! change q due to immersion droplets freez
     *   +execute_bulk2m_driver('get','mnucci')*dtB2M
         mnuccd  =             ! change q freezing aerosol (prim ice nuc)
     *   +execute_bulk2m_driver('get','mnuccd') ! *dtB2M
         mnucmd =              ! change q cond freezing Meyer's (prim ice nuc)
     *   +execute_bulk2m_driver('get','mnucmd') ! *dtB2M
         mnucmt =              ! change q cont freezing Meyer's (prim ice nuc)
     *   +execute_bulk2m_driver('get','mnucmt') ! *dtB2M

         qi_tot = mnuccc + mnucci + mnuccd + mnucmd + mnucmt
c
c g/m^3
c
        if(lSCM) then
         write(6,*)"stop BLK: mcrys_old,qi_tot,mcrys"
     *,l,mcrys_old*1.e+3,qi_tot*1.e+3,mcrys*1.e+3
         write(6,*)"stop BLK: mnuccc,mnucci,mnuccd,mnucmd,mnucmt"
     *,l,mnuccc*1.e+3,mnucci*1.e+3,mnuccd*1.e+3,mnucmd*1.e+3
     *,mnucmt*1.e+3
         endif
c
c balanced tendecies:
c
         ldummy=execute_bulk2m_driver('all','make_balance',mkx)
c
         nc_tnd=dtB2M*execute_bulk2m_driver('get','tnd','nc_tnd')
        if(lSCM) then
         write(6,*)"stop BLK:00: nc_tnd",l,nc_tnd*1.e-6
         endif
c
         npccn  =             ! change n droplets activation
     * +execute_bulk2m_driver('get','npccn')*dtB2M
         nprc   =             ! change n autoconversion of droplets:
     * -execute_bulk2m_driver('get','nprc')*dtB2M
         nnuccc =             ! change n due to con droplets freez
     * -execute_bulk2m_driver('get','nnuccc')*dtB2M
         nnucci =             ! change n due to imm droplets freez
     * -execute_bulk2m_driver('get','nnucci')*dtB2M

         nc_tnd = npccn + nprc + nnuccc + nnucci
        if(lSCM) then
         write(6,*)"stop BLK:01: nc_tnd",l,nc_tnd*1.e-6
        endif
c
         qc_tnd=dtB2M*execute_bulk2m_driver('get','tnd','qc_tnd')
        if(lSCM) then
         write(6,*)"stop BLK:00: qc_tnd",l,qc_tnd*1.e+3
        endif
c
         mpccn  =              ! change q droplets activation
     *    execute_bulk2m_driver('get','mpccn')*dtB2M
         mprc   =              ! change q autoconversion of droplets:
     *   -execute_bulk2m_driver('get','mprc')*dtB2M
         mnuccc =              ! change q due to con droplets freez
     *   -execute_bulk2m_driver('get','mnuccc')*dtB2M
         mnucci =              ! change q due to imm droplets freez
     *   -execute_bulk2m_driver('get','mnucci')*dtB2M

         qc_tnd = mpccn + mprc + mnuccc + mnucci
        if(lSCM) then
         write(6,*)"stop BLK:01: qc_tnd",l,qc_tnd*1.e+3
        endif
c
         ni_tnd=dtB2M*execute_bulk2m_driver('get','tnd','ni_tnd')
        if(lSCM) then
         write(6,*)"stop BLK:01: ni_tnd",l,ni_tnd*1.e-3
        endif
c
         qi_tnd=dtB2M*execute_bulk2m_driver('get','tnd','qi_tnd')
        if(lSCM) then
         write(6,*)"stop BLK:01: qi_tnd",l,qi_tnd*1.e+3
        endif
c
         ndrop_new = ndrop_old + nc_tnd
         mdrop_new = mdrop_old + qc_tnd
         ncrys_new = ncrys_old + ni_tnd
         mcrys_new = mcrys_old + qi_tnd
c
         ndrop_blk = execute_bulk2m_driver('get','val','nc')
         mdrop_blk = execute_bulk2m_driver('get','val','qc')
         ncrys_blk = execute_bulk2m_driver('get','val','ni')
         mcrys_blk = execute_bulk2m_driver('get','val','qi')
c
         ndrop_res = ndrop_blk + nc_tnd
         mdrop_res = mdrop_blk + qc_tnd
         ncrys_res = ncrys_blk + ni_tnd
         mcrys_res = mcrys_blk + qi_tnd
c
        if(wSCM) then
         write(6,*)
     *   "stop BLK: ndrop_old,nc_tnd,ndrop_new"
     *  ,l,ndrop_old*1.e-6,nc_tnd*1.e-6,ndrop_new*1.e-6
c
         write(6,*)
     *   "stop BLK: mdrop_old,qc_tnd,mdrop_new"
     *  ,l,mdrop_old*1.e+3,qc_tnd*1.e3,mdrop_new*1.e+3
c
         write(6,*)
     *   "stop BLK: ncrys_old,ni_tnd,ncrys_new"
     *  ,l,ncrys_old*1.e-3,ni_tnd*1.e-3,ncrys_new*1.e-3
c
         write(6,*)
     *   "stop BLK: mcrys_old,qi_tnd,mcrys_new"
     *  ,l,mcrys_old*1.e+3,qi_tnd*1.e3,mcrys_new*1.e+3
c
         write(6,*)
     *   "stop BLK: ndrop_blk,nc_tnd,ndrop_res"
     *  ,l,ndrop_blk*1.e-6,nc_tnd*1.e-6,ndrop_res*1.e-6
c
         write(6,*)
     *   "stop BLK: mdrop_blk,qc_tnd,mdrop_res"
     *  ,l,mdrop_blk*1.e+3,qc_tnd*1.e3,mdrop_res*1.e+3
c
         write(6,*)
     *   "stop BLK: ndrop_old,ndrop,ndrop_new,nc_tot,nc_tnd"
     *  ,l,ndrop_old*1.e-6,ndrop_new*1.e-6,ndrop*1.e-6
     *  ,nc_tot*1.e-6,nc_tnd*1.e-6
c
         write(6,*)
     *   "stop BLK: mdrop_old,mdrop,mdrop_new,qc_tot,qc_tnd"
     *  ,l,mdrop_old*1.e+3,mdrop_new*1.e+3,mdrop*1.e+3
     *  ,qc_tot*1.e+3,qc_tnd*1.e+3
        endif
c
c output for standalone internal variables
c
         if(lSCM) then
         write(iuo,*) 'l dtB2M'
         write(iuo,*) l
         write(iuo,*) dtB2M

         write(iuo,*) 'wmx wmxice tl ql pl svlhxl lhx '
         write(iuo,*) wmx
         write(iuo,*) wmxice
         write(iuo,*) tl
         write(iuo,*) ql
         write(iuo,*) pl
         write(iuo,*) svlhxl
         write(iuo,*) lhx

         write(iuo,*) 'ndrop_old mdrop_old ncrys_old mcrys_old'
         write(iuo,*) ndrop_old
         write(iuo,*) mdrop_old
         write(iuo,*) ncrys_old
         write(iuo,*) mcrys_old

         write(iuo,*) 'nc_tnd qc_tnd ni_tnd qi_tnd'
         write(iuo,*) nc_tnd
         write(iuo,*) qc_tnd
         write(iuo,*) ni_tnd
         write(iuo,*) qi_tnd

         write(iuo,*) 'ndrop_new mdrop_new ncrys_new mcrys_new'
         write(iuo,*) ndrop_new
         write(iuo,*) mdrop_new
         write(iuo,*) ncrys_new
         write(iuo,*) mcrys_new

         write(iuo,*) 'ndrop mdrop ncrys mcrys'
         write(iuo,*) ndrop
         write(iuo,*) mdrop
         write(iuo,*) ncrys
         write(iuo,*) mcrys

         write(iuo,*) 'ndrop_res mdrop_res ncrys_res mcrys_res'
         write(iuo,*) ndrop_res
         write(iuo,*) mdrop_res
         write(iuo,*) ncrys_res
         write(iuo,*) mcrys_res

         write(iuo,*) 'npccn nprc nnuccc nnucci nc_tot'
         write(iuo,*) npccn
         write(iuo,*) nprc
         write(iuo,*) nnuccc
         write(iuo,*) nnucci
         write(iuo,*) nc_tot

         write(iuo,*) 'mpccn mprc mnuccc mnucci qc_tot'
         write(iuo,*) mpccn
         write(iuo,*) mprc
         write(iuo,*) mnuccc
         write(iuo,*) mnucci
         write(iuo,*) qc_tot

         write(iuo,*) 'nnuccd nnucmd nnucmt ni_tot'
         write(iuo,*) nnuccd
         write(iuo,*) nnucmd
         write(iuo,*) nnucmt
         write(iuo,*) ni_tot

         write(iuo,*) 'mnuccd mnucmd mnucmt qi_tot'
         write(iuo,*) mnuccd
         write(iuo,*) mnucmd
         write(iuo,*) mnucmt
         write(iuo,*) qi_tot
c
         endif
c
        if( (ndrop_res(mkx) .lt. 0) .or. (mdrop_res(mkx) .lt. 0)  .or.
     *      (ncrys_res(mkx) .lt. 0) .or. (mcrys_res(mkx) .lt. 0)) then
c        call stop_model("BLK2MOM: Negative conc/cont...", 255)
c         write(6,*)"We reached -ve con.",ndrop_res(mkx),mdrop_res(mkx),
c     * ncrys_res(mkx), mcrys_res(mkx),l
         ndrop_res(mkx)=20.*1.d06
         mdrop_res(mkx)=1*1.d-06
         ncrys_res(mkx)=1*1.d-06
         mcrys_res(mkx)=1*1.d02
        else
         ndrop=ndrop_res;mdrop=mdrop_res;ncrys=ncrys_res;mcrys=mcrys_res
        endif
c
        endif if_balance
c
c To calculate "new" temperature and vapor mixing ratio:
        ldummy=execute_bulk2m_driver('tkqv','tk_qv',mkx)
        tk0new=tk0+dtB2M*execute_bulk2m_driver('get','tnd','tk_tnd')
        qk0new=qk0+dtB2M*execute_bulk2m_driver('get','tnd','qv_tnd')
c
c At this point you have 2 phases separately.
c Almost all processes are switched off, but you can calculate also
c accreation of droplets/crystal  by rain/snow, for example, and use
c rain/snow as diagnostic variables. But you need one additional
c long-storage array to keep ice crystal hydrometeor content as a minimum
c
!     IF(LHX.EQ.LHE)  THEN
!        WMX(L)=mdrop(mkx)
!      ELSE
!        WMX(L)=mcrys(mkx)
!      ENDIF
c
c GCM logics...........  SNd0, SNdL [No/cc; ]SNdI Units also in /cc
c
c      SNdI=ncrys(mkx)*1.0d-6          ! ncrys, [No/m^3]
       SNdI = 0.06417127d0
c      if(SNdI.gt.0.) write(6,*)"ICE CRY",SNdI, SNdI/dtB2M
       if(SNdI.gt.1.d0) SNdI=1.d0      !try to limit to 1000 /l
       SNd=ndrop(mkx)*1.d-6                 ! ndrop, [No/m^3]
c      if(SNd.gt.20.) write(6,*)"SM 12 CDNC",SNd   ,l
C**** Old treatment to get CDNC for cloud changes within time steps
       DCLD(L) = FCLD-CLDSAVL(L) ! cloud fraction change
C** If previous time step is clear sky
       if(CLDSAVL(L).eq.0.) then
         SNd=SNd
       elseif (DCLD(L).le.0.d0) then
         SNd=OLDCDL(L)
       elseif(DCLD(L).gt.0.d0) then
         SNd=( (OLDCDL(L)*CLDSAVL(L)) + (SNd*DCLD(L)) )/FCLD
       endif
C* If using an alternate definition for QAUT
       rablk=execute_bulk2m_driver('get','mprc')
       QAUT_B2M=rablk(mkx)
#endif  /* (cld-aer and BLK_2MOM) */
#ifdef CLD_AER_CDNC
      SCDNCW=SNd      ! we have already passed the grid box value
      SCDNCI=SNdI
      if (SCDNCI.le.0.0d0) SCDNCI=teeny         !set min ice crystal, do we need this, please check
      if (SCDNCW.le.20.d0) SCDNCW=20.d0         !set min CDNC, sensitivity test
c     if(SCDNCW.gt.2000.) write(6,*)"PROBLEM",SCDNCW,L
      if (SCDNCW.ge.1400.d0) SCDNCW=1400.d0     !set max CDNC, sensitivity test
c     if (SNd.gt.20.) write(6,*)"CDNC LSS",SCDNCW,SNd,L
#endif
C**** COMPUTE THE AUTOCONVERSION RATE OF CLOUD WATER TO PRECIPITATION
      IF(WMX(L).GT.0.) THEN
        RHO=1d5*PL(L)/(RGAS*TL(L))
        TEM=RHO*WMX(L)/(WCONST*FCLD+teeny)
        IF(LHX.EQ.LHS ) TEM=RHO*WMX(L)/(WMUI*FCLD+teeny)
C       IF(LHX.EQ.LHE.AND.ROICE.GT..1d0) TEM=RHO*WMX(L)
C    *    /(WMUSI*FCLD+teeny)
        TEM=TEM*TEM
        IF(TEM.GT.10.) TEM=10.
        IF(VDEF.GT.0..AND.RHO*WMX(L).GE.10.0d0) CM0=CM00
        CM1=CM0
        IF(BANDF) CM1=CM0*CBF
        IF(LHX.EQ.LHS) CM1=CM0
        CM=CM1*(1.-1./EXP(TEM*TEM))+100.*(PREBAR(L+1)+
     *       PRECNVL(L+1)*BYDTsrc)
CC#ifdef CLD_AER_CDNC
C** Choice of 2 different routines to get the autoconversion rate
CC#ifdef BLK_2MOM
C*** using an alternate QAUT definition based on Beheng (1994)
CC        CM=QAUT_B2M/(WMX(L)+1.d-20)+1.d0*100.d0*(PREBAR(L+1)+
CC     *     PRECNVL(L+1)*BYDTsrc)
c     if (QAUT_B2M.lt.0.) write(6,*)"QAUT BLK_2M",QAUT_B2M,CM,WMX(L),L
c       if(L.eq.1) write(6,*)"4th check BLK_2M",CM,QAUT_B2M,WMX(L)
CC#else
C** Use Qaut definition based on Rotstayn and Liu (2005, GRL)
CC         WTEM=1d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+teeny)
CC          IF(LHX.EQ.LHE)  THEN
CC            RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCW))**BY3
CC          ELSE
CC            RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCI))**BY3
C    *         *(1.+pl(l)*xRICld)
CC          END IF
CC       CALL GET_QAUT(L,PL(L),TL(L),FCLD,WMX(L),SCDNCW,RCLD,RHOW,
CC     *r6,r6c,QCRIT,QAUT)
C** Can also use other Qaut definitions if BLK_2MOM is not defined by switching to this call
!     CALL GET_QAUT(L,TL(L),FCLD,WMX(L),SCDNCW,RHO,QCRIT,QAUT)
!      CALL GET_QAUT(L,FCLD,WMX(L),SCDNCW,RHO,QAUT)
C*** If 6th moment of DSD is greater than critical radius r6c start QAUT
!     if (r6.gt.r6c) then
CC      if ((WMX(L)/(FCLD+teeny)).GT.QCRIT) then
CC        CM=QAUT/(WMX(L)+1.d-20)+1.d0*100.d0*(PREBAR(L+1)+
CC     *     PRECNVL(L+1)*BYDTsrc)
CC      else
CC        CM=0.d0
CC      endif
C** end routine for QAUT as a function of N,LWC
CC#endif
CC#endif
        CM=CM*CMX
        IF(CM.GT.BYDTsrc) CM=BYDTsrc
        PREP(L)=WMX(L)*CM
#ifdef SCM
       if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
           PRESAV(L)=PREP(L)*DTsrc
       endif
#endif
        IF(TL(L).LT.TF.AND.LHX.EQ.LHE) THEN ! check snowing pdf
          PRATM=1d5*COEFM*WMX(L)*PL(L)/(WCONST*FCLD*TL(L)*RGAS+teeny)
          PRATM=MIN(PRATM,1d0)*(1.-EXP(MAX(-1d2,(TL(L)-TF)/COEFT)))
          IF(PRATM.GT.RNDSSL(3,L)) LHP(L)=LHS
        END IF
      ELSE
        CM=0.
      END IF
C**** DECIDE WHETHER TO FORM CLOUDS
C**** FORM CLOUDS ONLY IF RH GT RH00
      IF (RH1(L).LT.RH00(L)) THEN
        FORM_CLOUDS=.FALSE.
      ELSE   ! COMPUTE THE CONVERGENCE OF AVAILABLE LATENT HEAT
        SQ(L)=LHX*QSATL(L)*DQSATDT(TL(L),LHX)*BYSHA
        TEM=-LHX*DPDT(L)/PL(L)
        QCONV=LHX*AQ(L)-RH(L)*SQ(L)*SHA*PLK(L)*ATH(L)
     *       -TEM*QSATL(L)*RH(L)
        FORM_CLOUDS= (QCONV.GT.0. .OR. WMX(L).GT.0.)
      END IF
C****
      ERMAX=LHX*PREBAR(L+1)*GRAV*BYAM(L)
      IF (FORM_CLOUDS) THEN
C**** COMPUTE EVAPORATION OF RAIN WATER, ER
        RHN=MIN(RH(L),RHF(L))
        IF(WMX(L).GT.0.)  THEN
          ER(L)=(1.-RHN)**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
        ELSE                    !  WMX(l).le.0.
          IF(PREICE(L+1).GT.0..AND.TL(L).LT.TF)  THEN
            ER(L)=(1.-RHI)**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
          ELSE
            ER(L)=(1.-RH(L))**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
          END IF
        END IF
        ER(L)=MAX(0d0,MIN(ER(L),ERMAX))
C**** COMPUTATION OF CLOUD WATER EVAPORATION
        IF (CLEARA(L).GT.0.) THEN
          WTEM=1d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+teeny)
          IF(LHX.EQ.LHE.AND.WMX(L)/FCLD.GE.WCONST*1d-3)
     *         WTEM=1d2*WCONST*PL(L)/(TL(L)*RGAS)
          IF(WTEM.LT.1d-10) WTEM=1d-10
          IF(LHX.EQ.LHE)  THEN
!           RCLD=1d-6*(RWCLDOX*10.*(1.-PEARTH)+7.*PEARTH)*(WTEM*4.)**BY3
            RCLD=RCLDX*1d-6*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCW))**BY3
          ELSE
!           RCLD=25.d-6*(WTEM/4.2d-3)**BY3 * (1.+pl(l)*xRICld)
            RCLD=RCLDX*100.d-6*(WTEM/(2.d0*BY3*TWOPI*SCDNCI))**BY3
            RCLD=MIN(RCLD,RIMAX)
          END IF
          CK1=1000.*LHX*LHX/(2.4d-2*RVAP*TL(L)*TL(L))
          CK2=1000.*RGAS*TL(L)/(2.4d-3*QSATL(L)*PL(L))
          TEVAP=COEEC*(CK1+CK2)*RCLD*RCLD
          WMX1=WMX(L)-PREP(L)*DTsrc
          WMPR(L)=PREP(L)*DTsrc             ! precip water
          ECRATE=(1.-RHF(L))/(TEVAP*FCLD+teeny)
          IF(ECRATE.GT.BYDTsrc) ECRATE=BYDTsrc
          EC(L)=WMX1*ECRATE*LHX
        END IF
C**** COMPUTE NET LATENT HEATING DUE TO STRATIFORM CLOUD PHASE CHANGE,
C**** QHEAT, AND NEW CLOUD WATER CONTENT, WMNEW
        DRHDT=2.*CLEARA(L)*CLEARA(L)*(1.-RH00(L))*(QCONV+ER(L))/LHX/
     *       (WMX(L)/(FCLD+teeny)+2.*CLEARA(L)*QSATL(L)*(1.-RH00(L))
     *       +teeny)
        IF(ER(L).EQ.0.AND.WMX(L).LE.0.) DRHDT=0.
        QHEAT(L)=FSSL(L)*(QCONV-LHX*DRHDT*QSATL(L))/(1.+RH(L)*SQ(L))
c        if (debug) print*,"ls2",l,qheat(l)
        DWDT=QHEAT(L)/LHX-PREP(L)+CLEARA(L)*FSSL(L)*ER(L)/LHX
        WMNEW =WMX(L)+DWDT*DTsrc
        IF(WMNEW.LT.0.) THEN
          WMNEW=0.
          QHEAT(L)=(-WMX(L)*BYDTsrc+PREP(L))*LHX-CLEARA(L)*FSSL(L)*ER(L)
c        if (debug) print*,"ls3",l,qheat(l)
        END IF
      ELSE
C**** UNFAVORABLE CONDITIONS FOR CLOUDS TO EXIT, PRECIP OUT CLOUD WATER
        QHEAT(L)=0.
        IF (WMX(L).GT.0.) THEN

          call get_dq_evap(tl(l)*rh00(l)/plk(l),ql(l)*rh00(l),plk(l)
     *         ,rh00(l),lhx,pl(l),wmx(l)/(fssl(l)*rh00(l)),dqsum,fqcond1
     *         )
          DWDT=DQSUM*RH00(L)*FSSL(L)

C**** DWDT is amount of water going to vapour, store LH (sets QNEW below)
          QHEAT(L)=-DWDT*LHX*BYDTsrc
c        if (debug) print*,"ls4",l,qheat(l)
          PREP(L)=MAX(0d0,(WMX(L)-DWDT)*BYDTsrc) ! precip out cloud water
          WMPR(L)=PREP(L)*DTsrc ! precip water (for opt. depth calculation)
#ifdef SCM
            if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
               PRESAV(L)=PREP(L)*DTsrc
            endif
#endif
        END IF
        ER(L)=(1.-RH(L))**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
        IF(PREICE(L+1).GT.0..AND.TL(L).LT.TF)
     *       ER(L)=(1.-RHI)**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
        ER(L)=MAX(0d0,MIN(ER(L),ERMAX))
        QHEAT(L)=QHEAT(L)-CLEARA(L)*FSSL(L)*ER(L)
c        if (debug) print*,"ls5",l,qheat(l),lhx
        WMNEW=0.
      END IF

C**** PHASE CHANGE OF PRECIPITATION, FROM ICE TO WATER
C**** This occurs if 0 C isotherm is crossed, or ice is falling into
C**** a super-cooled water cloud (that has not had B-F occur).
C**** Note: on rare occasions we have ice clouds even if T>0
C**** In such a case, no energy of phase change is needed.
      HPHASE=0.
      IF (LHP(L+1).EQ.LHS.AND.LHP(L).EQ.LHE.AND.PREICE(L+1).GT.0) THEN
        HPHASE=HPHASE+LHM*PREICE(L+1)*GRAV*BYAM(L)
        PREICE(L+1)=0.
c        if (debug) print*,"ls6",l,hphase,lhp(l),lhp(l+1),lhe
      END IF
C**** PHASE CHANGE OF PRECIP, FROM WATER TO ICE
      IF (LHP(L+1).EQ.LHE.AND.LHP(L).EQ.LHS.AND.PREBAR(L+1).GT.0) THEN
         HPHASE=HPHASE-LHM*PREBAR(L+1)*GRAV*BYAM(L)
c         if (debug) print*,"ls7",l,hphase,lhp(l),lhp(l+1),lhe
      END IF

C**** Make sure energy is conserved for transfers between P and CLW
      IF (LHP(L).NE.LHX) THEN
         HPHASE=HPHASE+(ER(L)*CLEARA(L)*FSSL(L)/LHX-PREP(L))*LHM
c         if (debug) print*,"ls8",l,hphase
      END IF
C**** COMPUTE THE PRECIP AMOUNT ENTERING THE LAYER TOP
      IF (ER(L).eq.ERMAX) THEN ! to avoid round off problem
        PREBAR(L)=PREBAR(L+1)*(1.-CLEARA(L)*FSSL(L))+
     *            AIRM(L)*PREP(L)*BYGRAV
      ELSE
        PREBAR(L)=MAX(0d0,PREBAR(L+1)+
     *       AIRM(L)*(PREP(L)-ER(L)*CLEARA(L)*FSSL(L)/LHX)*BYGRAV)
      END IF
C**** UPDATE NEW TEMPERATURE AND SPECIFIC HUMIDITY
      QNEW =QL(L)-DTsrc*QHEAT(L)/(LHX*FSSL(L)+teeny)
      IF(QNEW.LT.0.) THEN
        QNEW=0.
        QHEAT(L)=QL(L)*LHX*BYDTsrc*FSSL(L)
c        if (debug) print*,"ls9",l,qheat(l)
        DWDT1=QHEAT(L)/LHX-PREP(L)+CLEARA(L)*FSSL(L)*ER(L)/LHX
        WMNEW=WMX(L)+DWDT1*DTsrc
C**** IF WMNEW .LT. 0., THE COMPUTATION IS UNSTABLE
        IF(WMNEW.LT.0.) THEN
          IERR=1
          LERR=L
          WMERR=WMNEW
          WMNEW=0.
        END IF
      END IF
C**** Only Calculate fractional changes of Q to W
#ifdef TRACERS_WATER
      FPR=0.
      IF (WMX(L).gt.0.) FPR=PREP(L)*DTsrc/WMX(L)              ! CLW->P
      FPR=MIN(1d0,FPR)
      FER=0.
      IF (PREBAR(L+1).gt.0.) FER=CLEARA(L)*FSSL(L)*ER(L)*AIRM(L)/
     *     (GRAV*LHX*PREBAR(L+1))                             ! P->Q
      FER=MIN(1d0,FER)
      FWTOQ=0.                                                ! CLW->Q
#endif
      FQTOW=0.                                                ! Q->CLW
      IF (FSSL(L).gt.0) THEN
      IF (QHEAT(L)+CLEARA(L)*FSSL(L)*ER(L).gt.0) THEN
        IF (LHX*QL(L)+DTsrc*CLEARA(L)*ER(L).gt.0.) FQTOW=(QHEAT(L
     *      )+CLEARA(L)*FSSL(L)*ER(L))*DTsrc/((LHX*QL(L)+DTsrc*CLEARA(L)
     *       *ER(L))*FSSL(L))
#ifdef TRACERS_WATER
      ELSE
        IF (WMX(L)-PREP(L)*DTsrc.gt.0.) FWTOQ=-(QHEAT(L)
     *      +CLEARA(L)*FSSL(L)*ER(L))*DTsrc/(LHX*(WMX(L)-PREP(L)*DTsrc))
        FWTOQ=MIN(1d0,FWTOQ)
#endif
      END IF
      END IF
      QL(L)=QNEW
C**** adjust gradients down if Q decreases
      QMOM(:,L)= QMOM(:,L)*(1.-FQTOW)
      WMX(L)=WMNEW
!     if(abs(DTsrc*(QHEAT(L)-HPHASE)).gt.100.*SHA*FSSL(L)) then ! warn
!       write(0,*) 'it,i,j,l,tlold,dtl,qht,hph',itime,i_debug,j_debug,
!    *    L,TL(L),DTsrc*(QHEAT(L)-HPHASE)/(SHA*FSSL(L)+teeny),
!    *    QHEAT(L),HPHASE
!       debug_out=.true.
!     end if
      TL(L)=TL(L)+DTsrc*(QHEAT(L)-HPHASE)/(SHA*FSSL(L)+teeny)
c        if (debug) print*,"lsA",l,qheat(l),hphase,lhx

      TH(L)=TL(L)/PLK(L)
      QSATC=QSAT(TL(L),LHX,PL(L))
      RH(L)=QL(L)/QSATC
#ifdef TRACERS_WATER
C**** update tracers from cloud formation (in- and below-cloud
C****    precipitation, evaporation, condensation, and washout)
c CLDSAVT is current FCLD
        IF(RH(L).LE.1.) THEN
          IF (RH00(L).lt.1.) then
            CLDSAVT=1.-DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
          ELSE
            CLDSAVT=0.
          END IF
        END IF
        IF(CLDSAVT.LT.0.) CLDSAVT=0.
        IF(RH(L).GT.1.) CLDSAVT=1.
        IF (CLDSAVT.GT.1.) CLDSAVT=1.
        IF (WMX(L).LE.0.) CLDSAVT=0.
        CLDSAVT=CLDSAVT*FSSL(L)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      WA_VOL=0.
      IF (WMNEW.GT.teeny) THEN
        WA_VOL=WMNEW*AIRM(L)*1.D2*BYGRAV*DXYPIJ
      END IF
      WMXTR = WMX(L)
      IF (BELOW_CLOUD.and.WMX(L).LT.teeny) THEN
        precip_mm = PREBAR(L+1)*100.*DTsrc
        if (precip_mm.lt.0.) precip_mm=0.
        WMXTR = PREBAR(L+1)*grav*BYAM(L)*dtsrc
        if (wmxtr.lt.0.) wmxtr=0.
        WA_VOL=precip_mm*DXYPIJ
      END IF

      CALL GET_SULFATE(L,TL(L),FCLD,WA_VOL
     * ,WMXTR,SULFIN,SULFINC,SULFOUT,TR_LEFT,TM,TRWML(1,l),AIRM,LHX
     *  ,DT_SULF_SS(1,L),CLDSAVT)

      do iaqch=1,aqchem_count
        n = aqchem_list(iaqch)
        TR_LEF(n)=TR_LEFT(N)
        TRWML(N,L)=TRWML(N,L)*(1.+SULFINC(N))
        TM(L,N)=TM(L,N)*(1.+SULFIN(N))
        TMOM(:,L,N)  = TMOM(:,L,N)*(1. +SULFIN(N))
        if (WMX(L).LT.teeny.and.BELOW_CLOUD) then
          TRPRBAR(N,L+1)=TRPRBAR(N,L+1)+SULFOUT(N)
        else
          TRWML(N,L) = TRWML(N,L)+SULFOUT(N)
        endif
      enddo

#endif

c precip. tracer evap
      FPRT = FPR
      DTPRT(1:NTX) = FPRT  *TRWML(1:NTX,L)

      if(fer.ne.0.) then
        CALL GET_EVAP_FACTOR_array(
     &       NTX,TL(L),LHX,.FALSE.,1d0,FER,FERT,ntix)
        DTERT(1:NTX) = FERT(1:NTX)  *TRPRBAR(1:NTX,L+1)
      else
        FERT(1:NTX) = 0.
        DTERT(1:NTX) = 0.
      endif

      if(fwtoq.ne.0.) then
        CALL GET_EVAP_FACTOR_array(
     &       NTX,TL(L),LHX,.FALSE.,1d0,FWTOQ,FWTOQT,ntix)
        FPRT=FPR
        DO N=1,NTX
          DTQWT(N) = -FWTOQT(N)*TRWML(N,L)*(1.-FPRT)
        ENDDO
      else
        FWTOQT(1:NTX) = 0.
        DTQWT(1:NTX) = 0.
      endif

      TM_dum(:) = TM(L,:)
      IF(BELOW_CLOUD.and.WMX(L).lt.teeny) THEN
        FQTOWT(:)=0.
        THLAW(gases_list)=0.
        precip_mm = PREBAR(L+1)*100.*dtsrc
        WMXTR = PREBAR(L+1)*grav*BYAM(L)*dtsrc
        if (precip_mm.lt.0.) precip_mm=0.
        if (wmxtr.lt.0.) wmxtr=0.
        CALL GET_WASH_FACTOR_array(NTX,b_beta_DT,precip_mm,FWASHT,
     &       tl(l),LHX,WMXTR,cldprec,TM_dum,TRPRBAR(:,l),
     &       THWASH,pl(l),ntix,.true.) !washout
      ELSE
c         b_beta_DT is needed at the lowest precipitating level,
c         so saving it here for below cloud case:
        b_beta_DT = cldsavt  !*CM*dtsrc
c         saves cloud fraction at lowest precipitating level for washout
        cldprec=cldsavt
cdmk added arguments above; THLAW added below (no way to factor this)
        WMXTR = WMX(L)
        CALL GET_COND_FACTOR_array(
     &        NTX,WMXTR,TL(L),TL(L),LHX,FCLD,FQTOW
     &       ,FQTOWT,.false.,TRWML(:,L),TM_dum,THLAW,TR_LEF,PL(L)
     &       ,ntix,CLDSAVT)
        DO N=1,NTX
          DTQWT(N) = DTQWT(N)
     &         +FQTOWT(N)*TR_LEF(N)*(TM_dum(N)+DTERT(N))
        ENDDO
#ifdef NO_WASHOUT_IN_CLOUDS
        THWASH(:)=0.
        FWASHT(:)=0.
#endif
      ENDIF

#ifndef NO_WASHOUT_IN_CLOUDS
c**** washout in clouds
c apply certain removal processes before calculating washout in clouds
      tm_dum(1:ntx)=tm_dum(1:ntx)-dtqwt(1:ntx)-thlaw(1:ntx)
      IF(.NOT.(BELOW_CLOUD.and.WMX(L).lt.teeny)) THEN
        precip_mm = PREBAR(L+1)*100.*dtsrc
        WMXTR = PREBAR(L+1)*grav*BYAM(L)*dtsrc
        if (precip_mm.lt.0.) precip_mm=0.
        if (wmxtr.lt.0.) wmxtr=0.
        CALL GET_WASH_FACTOR_array(NTX,b_beta_DT,precip_mm,FWASHT,
     &       tl(l),LHX,WMXTR,cldprec,TM_dum,TRPRBAR(:,l),
     &       THWASH,pl(l),ntix,.true.) !washout
      END IF
#endif

#ifdef TRDIAG_WETDEPO
      IF (diag_wetdep == 1) THEN
        FPRT=FPR
        DO N=1,NTX
          trevap_ls(l,n)=dtert(n)
          trwash_ls(l,n)=fwasht(n)*tm_dum(n)+thwash(n)
          trclwc_ls(l,n)=fqtowt(n)*tr_lef(n)*(tm(l,n)+dtert(n))+thlaw(n)
          trprcp_ls(l,n)=dtprt(n)
          trclwe_ls(l,n)=fwtoqt(n)*trwml(n,l)*(1.-fprt)
        ENDDO
      ENDIF
#endif

      do igas=1,gases_count
        n = gases_list(igas)
        IF (TM(L,N).GT.teeny) THEN
          TMFAC(N)=THLAW(N)/TM(L,N)
        ELSE
          TMFAC(N)=0.
        END IF
#ifndef NO_WASHOUT_IN_CLOUDS
        IF (tm_dum(n).GT.teeny) THEN
          TMFAC2(N)=THWASH(N)/tm_dum(n) ! use tm instead?
        ELSE
          TMFAC2(N)=0.
        END IF
#endif
      enddo

      DO N=1,NTX

        dtwrt=fwasht(n)*tm_dum(n)

c ---------------------- apply fluxes ------------------------
        FPRT=FPR
        TRWML(N,L) = TRWML(N,L)*(1.-FPRT)  + DTQWT(N)+THLAW(N)

        TM(L,N) = MAX(0., TM(L,N)
     &       + DTERT(N) - DTWRT - DTQWT(N) - THLAW(N) - THWASH(N) )

        TRPRBAR(N,L)=TRPRBAR(N,L+1)*(1.-FERT(N))
     &       +DTPRT(N)+DTWRT+THWASH(N)

        TMOM(:,L,N)  = TMOM(:,L,N)*
     &       (1. - FQTOWT(N) - FWASHT(N) - TMFAC(N) - TMFAC2(N))

#ifdef TRACERS_SPECIAL_O18
C**** Isotopic equilibration of the CLW and water vapour
        IF (LHX.eq.LHE .and. WMX(L).gt.0) THEN  ! only if liquid
          CALL ISOEQUIL(NTIX(N),TL(L),.TRUE.,QL(L)*FSSL(L),WMX(L),
     *         TM(L,N),TRWML(N,L),1d0)
        END IF
C**** Isotopic equilibration of Precip (if liquid) and water vapour
C**** Note that precip is either all water or all ice
        IF (LHP(L).eq.LHE .AND. PREBAR(L).gt.0 .AND. QL(L).gt.0) THEN
          PRLIQ=PREBAR(L)*DTSrc*BYAM(L)*GRAV
          CALL ISOEQUIL(NTIX(N),TL(L),.TRUE.,QL(L)*FSSL(L),PRLIQ,
     *         TM(L,N),TRPRBAR(N,L),1d0)
        END IF
#endif
      END DO
      IF (PREBAR(L).eq.0) TRPRBAR(1:ntx,L)=0. ! remove round off error
      IF (WMX(L).eq.0)    TRWML(1:ntx,L)=0.   ! remove round off error
#endif
C**** CONDENSE MORE MOISTURE IF RELATIVE HUMIDITY .GT. 1
      RH1(L)=QL(L)/QSATC
      IF(LHX.EQ.LHS) THEN
        IF(RH(L).LE.1.) THEN
          IF (RH00(L).lt.1.) then
            CLEARA(L)=DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
          ELSE
            CLEARA(L)=1.
          END IF
        END IF
        IF(CLEARA(L).GT.1.) CLEARA(L)=1.
        IF(RH(L).GT.1.) CLEARA(L)=0.
        IF(WMX(L).LE.0.) CLEARA(L)=1.
        QF=(QL(L)-QSATC*(1.-CLEARA(L)))/(CLEARA(L)+teeny)
        IF(QF.LT.0.) WRITE(6,*) 'L CA QF Q QSA=',L,CLEARA(L),QF,QL(L),
     *               QSATC
        QSATE=QSAT(TL(L),LHE,PL(L))
        RHW=(2.583d0-TL(L)/207.83)*(QSAT(TL(L),LHS,PL(L))/QSATE)
        IF(TL(L).LT.238.16.AND.WMX(L).LE.0d0) RH1(L)=QF/(QSATE*RHW)
      END IF
      IF(RH1(L).GT.1.) THEN    ! RH was used in old versions
      SLH=LHX*BYSHA

      call get_dq_cond(tl(l),ql(l),1d0,1d0,lhx,pl(l),dqsum,fcond)

      IF(DQSUM.GT.0.) THEN
c        if (debug) print*,"lsB",l,slh*dqsum

      TL(L)=TL(L)+SLH*DQSUM
      QL(L)=QL(L)-DQSUM
      WMX(L)=WMX(L)+DQSUM*FSSL(L)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      WA_VOL=0.
      IF (WMX(L).GT.teeny) THEN
        WA_VOL=WMX(L)*AIRM(L)*1.D2*BYGRAV*DXYPIJ
      END IF
#endif
C**** adjust gradients down if Q decreases
      QMOM(:,L)= QMOM(:,L)*(1.-FCOND)
#ifdef TRACERS_WATER
C**** CONDENSING MORE TRACERS
      WMXTR = WMX(L)
        IF(RH(L).LE.1.) THEN
          IF (RH00(L).lt.1.) then
            CLDSAVT=1.-DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
          ELSE
            CLDSAVT=0.
          END IF
        END IF
        IF(CLDSAVT.LT.0.) CLDSAVT=0.
        IF(RH(L).GT.1.) CLDSAVT=1.
        IF (CLDSAVT.GT.1.) CLDSAVT=1.
        IF (WMX(L).LE.0.) CLDSAVT=0.
        CLDSAVT=CLDSAVT*FSSL(L)
cdmks  I took out some code above this that was for below cloud
c   processes - this should be all in-cloud
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)

      CALL GET_SULFATE(L,TL(L),FCLD,WA_VOL,WMXTR,SULFIN,
     *     SULFINC,SULFOUT,TR_LEFT,TM,TRWML(1,L),AIRM,LHX,
     *     DT_SULF_SS(1,L),CLDSAVT)

      do iaqch=1,aqchem_count
        n = aqchem_list(iaqch)
        TRWML(N,L)=TRWML(N,L)*(1.+SULFINC(N))
        TM(L,N)=TM(L,N)*(1.+SULFIN(N))
        TMOM(:,L,N) =TMOM(:,L,N)*(1.+SULFIN(N))
        TRWML(N,L) = TRWML(N,L)+SULFOUT(N)
        TR_LEF(N)=TR_LEFT(N)
      enddo

#endif
c below TR_LEFT(N) limits the amount of available tracer in gridbox
cdmkf and below, extra arguments for GET_COND, addition of THLAW
      TM_dum(:) = TM(L,:)
      CALL GET_COND_FACTOR_array(NTX,WMXTR,TL(L),TL(L),LHX,FCLD,FCOND
     &     ,FQCONDT,.false.,TRWML(:,L),TM_dum,THLAW,TR_LEF,pl(l)
     &     ,ntix,CLDSAVT)
      dtr(1:ntx) = fqcondt(1:ntx)*tm_dum(1:ntx)
#ifdef TRDIAG_WETDEPO
      IF (diag_wetdep == 1) trcond_ls(l,1:ntx) =
     &     dtr(1:ntx)+thlaw(1:ntx)
#endif
      do igas=1,gases_count
        n = gases_list(igas)
        IF (TM_dum(N).GT.teeny) THEN
          TMFAC(N)=THLAW(N)/TM_dum(N)
        ELSE
          TMFAC(N)=0.
        END IF
      enddo

      DO N=1,NTX
        TRWML(N,L)  =TRWML(N,L)+ DTR(N)+THLAW(N)
        TM(L,N)     =TM(L,N)    *(1.-FQCONDT(N))   -THLAW(N)
        TMOM(:,L,N) =TMOM(:,L,N)*(1.-FQCONDT(N) - TMFAC(N))
      END DO
#endif
      END IF
      RH(L)=QL(L)/QSAT(TL(L),LHX,PL(L))
      TH(L)=TL(L)/PLK(L)
!     if (debug_out) write(0,*) 'after condensation: l,tlnew,',l,tl(l)
      END IF

      IF(RH(L).LE.1.) THEN
        IF (RH00(L).lt.1.) THEN
          CLEARA(L)=DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
        ELSE
          CLEARA(L)=1.
        END IF
      END IF
      IF(CLEARA(L).GT.1.) CLEARA(L)=1.
      IF(RH(L).GT.1.) CLEARA(L)=0.
      IF(WMX(L).LE.0.) CLEARA(L)=1.
      IF(CLEARA(L).LT.0.) CLEARA(L)=0.
      RHF(L)=RH00(L)+(1.-CLEARA(L))*(1.-RH00(L))
      IF(RH(L).LE.RHF(L).AND.RH(L).LT..999999.AND.WMX(L).gt.0.) THEN
C**** PRECIP OUT CLOUD WATER IF RH LESS THAN THE RH OF THE ENVIRONMENT
#ifdef TRACERS_WATER
        TRPRBAR(1:NTX,L) = TRPRBAR(1:NTX,L) + TRWML(1:NTX,L)
#ifdef TRDIAG_WETDEPO
        IF (diag_wetdep == 1)
     &       trprcp_ls(l,1:ntx)=trprcp_ls(l,1:ntx)+trwml(1:ntx,l)
#endif
        TRWML(1:NTX,L) = 0.
#endif
        PREBAR(L)=PREBAR(L)+WMX(L)*AIRM(L)*BYGRAV*BYDTsrc
        IF(LHP(L).EQ.LHS .AND. LHX.EQ.LHE) THEN
          HCHANG=WMX(L)*LHM
c          if (debug) print*,"lsC",l,hchang
          TL(L)=TL(L)+HCHANG/(SHA*FSSL(L)+teeny)
!         if(debug_out) write(0,*) 'after rain out: l,tlnew',l,tl(l)
          TH(L)=TL(L)/PLK(L)
        END IF
        WMX(L)=0.
      END IF
      prebar1(l)=prebar(l)
#ifdef SCM
c     save precip LS (kg/kg)
      if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
          PRESAV(L) = PREBAR(L)*DTSRC*GRAV/AIRM(L)
          if (PREBAR(L).le.0) LHPSAV(L) = 0.
c         write(iu_scm_prt,*) 'L PRESAV PREBAR LHPSAV ',
c    &           L,PRESAV(l)*1000.,PREBAR(L),LHPSAV(L)
      endif
#endif
C**** set phase of condensation for next box down
      PREICE(L)=0.
      IF (PREBAR(L).gt.0 .AND. LHP(L).EQ.LHS) PREICE(L)=PREBAR(L)
      IF (PREBAR(L).le.0) LHP(L)=0.
C**** COMPUTE THE LARGE-SCALE CLOUD COVER
      IF(RH(L).LE.1.) THEN
        IF (RH00(L).lt.1.) THEN
          CLEARA(L)=DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
        ELSE
          CLEARA(L)=1.
        END IF
      END IF
      IF(CLEARA(L).GT.1.) CLEARA(L)=1.
      IF(RH(L).GT.1.) CLEARA(L)=0.
      IF(WMX(L).LE.teeny) WMX(L)=0.
      IF(WMX(L).LE.0.) CLEARA(L)=1.
      IF(CLEARA(L).LT.0.) CLEARA(L)=0.
      CLDSSL(L)=FSSL(L)*(1.-CLEARA(L))
      CLDSAVL(L)=1.-CLEARA(L)
#ifdef TRACERS_WATER
      IF(CLDSSL(L).gt.0.) CLOUD_YET=.true.
      IF(CLOUD_YET.and.CLDSSL(L).eq.0.) BELOW_CLOUD=.true.
#endif
      TOLDUP=TOLD
C**** ACCUMULATE SOME DIAGNOSTICS
         HCNDSS=HCNDSS+FSSL(L)*(TL(L)-TOLD)*AIRM(L)
         SSHR(L)=SSHR(L)+FSSL(L)*(TL(L)-TOLD)*AIRM(L)
      END DO  ! end of loop over L

      PRCPSS=MAX(0d0,PREBAR(1)*GRAV*DTsrc) ! fix small round off err
#ifdef TRACERS_WATER
      do n=1,ntx
        TRPRSS(n)=TRPRBAR(n,1)
        if(t_qlimit(n)) TRPRSS(n)=MAX(0d0,TRPRSS(n))
      enddo
#endif
C****
C**** CLOUD-TOP ENTRAINMENT INSTABILITY
C****
      DO L=LP50-1,1,-1
        SM(L)=TH(L)*AIRM(L)
        QM(L)=QL(L)*AIRM(L)
        WMXM(L)=WMX(L)*AIRM(L)
        SM(L+1)=TH(L+1)*AIRM(L+1)
        QM(L+1)=QL(L+1)*AIRM(L+1)
        WMXM(L+1)=WMX(L+1)*AIRM(L+1)
        IF(WMX(L+1).GT.teeny) CYCLE
        TOLD=TL(L)
        TOLDU=TL(L+1)
        QOLD=QL(L)
        QOLDU=QL(L+1)
        FCLD=(1.-CLEARA(L))*FSSL(L)+teeny
        IF(CLEARA(L).EQ.1. .OR. (CLEARA(L).LT.1..AND.CLEARA(L+1).LT.1.))
     *       CYCLE
C       IF(WMX(L).EQ.0. .OR. (WMX(L).GT.0..AND.WMX(L+1).GT.0.)) CYCLE
        SEDGE=THBAR(TH(L+1),TH(L))
        DSE=(TH(L+1)-SEDGE)*PLK(L+1)+(SEDGE-TH(L))*PLK(L)+
     *       SLHE*(QL(L+1)-QL(L))
        DWM=QL(L+1)-QL(L)+(WMX(L+1)-WMX(L))/FCLD
        DQSDT=DQSATDT(TL(L),LHE)*QL(L)/(RH(L)+1d-30)
        BETA=(1.+BYMRAT*TL(L)*DQSDT)/(1.+SLHE*DQSDT)
        CKM=(1.+SLHE*DQSDT)*(1.+(1.-DELTX)*TL(L)/SLHE)/
     *       (2.+(1.+BYMRAT*TL(L)/SLHE)*SLHE*DQSDT)
        CKR=TL(L)/(BETA*SLHE)
        CK=DSE/(SLHE*DWM+teeny)
        SIGK=0.
        IF(CKR.GT.CKM) CYCLE
        IF(CK.GT.CKR) SIGK=COESIG*((CK-CKR)/((CKM-CKR)+teeny))**5
        EXPST=EXP(-SIGK*DTsrc)
        IF(L.LE.1) CKIJ=EXPST
        DSEC=DWM*TL(L)/BETA
        IF(CK.LT.CKR) CYCLE
        FPMAX=MIN(1d0,1.-EXPST)       ! full CTE strength
        IF(FPMAX.LE.0.) CYCLE
        IF(DSE.GE.DSEC) CYCLE
C**** MIXING TO REMOVE CLOUD-TOP ENTRAINMENT INSTABILITY
c        if (debug) print*,"lse",l,wmxm(l:l+1),svlhxl(l:l+1)

        AIRMR=(AIRM(L+1)+AIRM(L))*BYAM(L+1)*BYAM(L)
        SMO1=SM(L)
        QMO1=QM(L)
        WMO1=WMXM(L)
        SMO2=SM(L+1)
        QMO2=QM(L+1)
        WMO2=WMXM(L+1)
        SMO12=SMO1*PLK(L)+SMO2*PLK(L+1)
        DO K=1,KMAX
          UMO1(K)=UM(K,L)
          VMO1(K)=VM(K,L)
          UMO2(K)=UM(K,L+1)
          VMO2(K)=VM(K,L+1)
        END DO
        FPLUME=FPMAX*FSSL(L)
        DFX=FPLUME ! not FPMAX
        DO ITER=1,9
          DFX=DFX*0.5
          FMIX=FPLUME*FCLD
          FMASS=FMIX*AIRM(L)
          FMASS=MIN(FMASS,(AIRM(L+1)*AIRM(L))/(AIRM(L+1)+AIRM(L)))
          FMIX=FMASS*BYAM(L)
          FRAT=FMASS*BYAM(L+1)
          SMN1=SMO1*(1.-FMIX)+FRAT*SMO2
          QMN1=QMO1*(1.-FMIX)+FRAT*QMO2
          WMN1=WMO1*(1.-FMIX)+FRAT*WMO2
          SMN2=SMO2*(1.-FRAT)+FMIX*SMO1
          QMN2=QMO2*(1.-FRAT)+FMIX*QMO1
          WMN2=WMO2*(1.-FRAT)+FMIX*WMO1
          THT1=SMN1*BYAM(L)/PLK(L)
          QLT1=QMN1*BYAM(L)
          TLT1=THT1*PLK(L)
          LHX=SVLHXL(L)
          RHT1=QLT1/(QSAT(TLT1,LHX,PL(L)))
          WMT1=WMN1*BYAM(L)
          THT2=SMN2*BYAM(L+1)/PLK(L+1)
          QLT2=QMN2*BYAM(L+1)
          WMT2=WMN2*BYAM(L+1)
          SEDGE=THBAR(THT2,THT1)
          DSE=(THT2-SEDGE)*PLK(L+1)+(SEDGE-THT1)*PLK(L)+SLHE*(QLT2-QLT1)
          DWM=QLT2-QLT1+(WMT2-WMT1)/FCLD
          DQSDT=DQSATDT(TLT1,LHE)*QLT1/(RHT1+1d-30)
          BETA=(1.+BYMRAT*TLT1*DQSDT)/(1.+SLHE*DQSDT)
          CKM=(1.+SLHE*DQSDT)*(1.+(1.-DELTX)*TLT1/SLHE)/
     *         (2.+(1.+BYMRAT*TLT1/SLHE)*SLHE*DQSDT)
          DSEC=DWM*TLT1/BETA
          DSEDIF=DSE-DSEC
c          if (debug) print*,"lsf",iter,dsedif,fplume,dfx,smn1,smn2
          IF(DSEDIF.GT.1d-3) FPLUME=FPLUME-DFX
          IF(DSEDIF.LT.-1d-3) FPLUME=FPLUME+DFX
          IF(ABS(DSEDIF).LE.1d-3.OR.FPLUME.GT.FPMAX*FSSL(L)) EXIT
        END DO
c        IF (FPLUME.GT.FPMAX) print*,"lsH",i_debug,j_debug,fplume,fpmax
c     $       ,fssl(l),dfx,iter

C**** UPDATE TEMPERATURE, SPECIFIC HUMIDITY AND MOMENTUM DUE TO CTEI
        SMN12=SMN1*PLK(L)+SMN2*PLK(L+1)
        SMN1=SMN1-(SMN12-SMO12)*AIRM(L)/((AIRM(L)+AIRM(L+1))*PLK(L))
        SMN2=SMN2-(SMN12-SMO12)*AIRM(L+1)/((AIRM(L)+AIRM(L+1))*PLK(L+1))
        TH(L)=SMN1*BYAM(L)
        TL(L)=TH(L)*PLK(L)
!       if(debug_out) write(0,*) 'after CTEI: l,tlnew',l,tl(l)
        QL(L)=QMN1*BYAM(L)
        LHX=SVLHXL(L)
        RH(L)=QL(L)/QSAT(TL(L),LHX,PL(L))
        WMX(L)=WMN1*BYAM(L)
        TH(L+1)=SMN2*BYAM(L+1)
        QL(L+1)=QMN2*BYAM(L+1)
        WMX(L+1)=WMN2*BYAM(L+1)
        CALL CTMIX (SM(L),SMOM(1,L),FMASS*AIRMR,FMIX,FRAT)
        CALL CTMIX (QM(L),QMOM(1,L),FMASS*AIRMR,FMIX,FRAT)
C****
#ifdef TRACERS_ON
        DO N=1,NTX
          CALL CTMIX (TM(L,N),TMOM(1,L,N),FMASS*AIRMR,FMIX,FRAT)
#ifdef TRACERS_WATER
C**** mix cloud liquid water tracers as well
          TWMTMP      = TRWML(N,L  )*(1.-FMIX)+FRAT*TRWML(N,L+1)
          TRWML(N,L+1)= TRWML(N,L+1)*(1.-FRAT)+FMIX*TRWML(N,L  )
          TRWML(N,L)  = TWMTMP
#endif
        END DO
#endif
        DO K=1,KMAX
          UMN1(K)=(UMO1(K)*(1.-FMIX)+FRAT*UMO2(K))
          VMN1(K)=(VMO1(K)*(1.-FMIX)+FRAT*VMO2(K))
          UMN2(K)=(UMO2(K)*(1.-FRAT)+FMIX*UMO1(K))
          VMN2(K)=(VMO2(K)*(1.-FRAT)+FMIX*VMO1(K))
          UM(K,L)=UM(K,L)+(UMN1(K)-UMO1(K))*RA(K)
          VM(K,L)=VM(K,L)+(VMN1(K)-VMO1(K))*RA(K)
          UM(K,L+1)=UM(K,L+1)+(UMN2(K)-UMO2(K))*RA(K)
          VM(K,L+1)=VM(K,L+1)+(VMN2(K)-VMO2(K))*RA(K)
        END DO
C**** RE-EVAPORATION OF CLW IN THE UPPER LAYER
        QL(L+1)=QL(L+1)+WMX(L+1)/(FSSL(L)+teeny)
c        if (wmxm(l+1).gt.0. .and. svlhxl(l+1).gt.0 .and. svlhxl(l+1).ne
c     $       .svlhxl(l)) print*,"lsT",i_debug,j_debug,wmxm(l)*svlhxl(l)
c     $       +wmxm(l+1)*svlhxl(l+1)-wmx(l)*airm(l)*svlhxl(l)-wmx(l+1)
c     $       *airm(l+1)*svlhxl(l+1),wmxm(l:l+1),svlhxl(l:l+1)
c**** assumes that wmx(l+1) is same phase as wmx(l)?
c**** energy fix?
        TH(L+1)=TH(L+1)-((SVLHXL(L+1)-SVLHXL(L))*BYSHA)*WMXM(L+1)/AIRM(L
     $       +1)/(PLK(L+1)*FSSL(L)+teeny)

        TH(L+1)=TH(L+1)-(LHX*BYSHA)*WMX(L+1)/(PLK(L+1)*FSSL(L)+teeny)
c        if (debug) print*,"lsD",l+1,LHX*WMX(L+1)
        TL(L+1)=TH(L+1)*PLK(L+1)
!       if(debug_out) write(0,*) 'after re-evap: l,tlnew',l+1,tl(l+1)
        RH(L+1)=QL(L+1)/QSAT(TL(L+1),LHX,PL(L+1))
        WMX(L+1)=0.
#ifdef TRACERS_WATER
        TM(L+1,1:NTX)=TM(L+1,1:NTX)+TRWML(1:NTX,L+1)
#ifdef TRDIAG_WETDEPO
        IF (diag_wetdep == 1)
     &       trclwe_ls(l+1,1:ntx)=trclwe_ls(l+1,1:ntx)+trwml(1:ntx,l+1)
#endif
        TRWML(1:NTX,L+1)=0.
#endif
        IF(RH(L).LE.1.) THEN
          IF (RH00(L).lt.1.) THEN
            CLEARA(L)=DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
          ELSE
            CLEARA(L)=1.
          END IF
        END IF
        IF(CLEARA(L).GT.1.) CLEARA(L)=1.
        IF(RH(L).GT.1.) CLEARA(L)=0.
        CLDSSL(L)=FSSL(L)*(1.-CLEARA(L))
        CLDSAVL(L)=1.-CLEARA(L)
        TNEW=TL(L)
        TNEWU=TL(L+1)
        QNEW=QL(L)
        QNEWU=QL(L+1)
        HCNDSS=HCNDSS+FSSL(L)*(TNEW-TOLD)*AIRM(L)+
     *         FSSL(L+1)*(TNEWU-TOLDU)*AIRM(L+1)
        SSHR(L)=SSHR(L)+FSSL(L)*(TNEW-TOLD)*AIRM(L)
        SSHR(L+1)=SSHR(L+1)+FSSL(L+1)*(TNEWU-TOLDU)*AIRM(L+1)
       DCTEI(L)=DCTEI(L)+FSSL(L)*(QNEW-QOLD)*AIRM(L)*LHX*BYSHA
       DCTEI(L+1)=DCTEI(L+1)+FSSL(L+1)*(QNEWU-QOLDU)*AIRM(L+1)*LHX*BYSHA
      END DO

C**** COMPUTE CLOUD PARTICLE SIZE AND OPTICAL THICKNESS
      WMSUM=0.
#ifdef SCM
      if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
          SCM_LWP_SS = 0.d0
          SCM_IWP_SS = 0.d0
      endif
#endif
#ifdef CLD_AER_CDNC
      ACDNWS=0.
      ACDNIS=0.
      AREWS=0.
      AREIS=0.
      ALWWS=0.
      ALWIS=0.
      NLSW = 0
      NLSI = 0
#endif
      DO L=1,LP50
        FCLD=CLDSSL(L)+teeny
        WTEM=1.d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+teeny)
        LHX=SVLHXL(L)
        IF(WTEM.LT.1d-10) WTEM=1.d-10
c        if (debug) print*,"lsY",l,wtem,wmx(l),fcld,tl(l)
C***Setting constant values of CDNC over land and ocean to get RCLD=f(CDNC,LWC)
      SNdO = 59.68d0/(RWCLDOX**3)
      SNdL = 174.d0
      SNdI = 0.06417127d0
      SCDNCW=SNdO*(1.-PEARTH)+SNdL*PEARTH
      SCDNCI=SNdI
#if (defined CLD_AER_CDNC ) || (defined BLK_2MOM)
#ifdef ALT_CDNC_INPUTS
       VVEL = VVEL_sv(l)    ! retrieve value at this level
       DSU(:) = DSU_SV(:,L) ! retrieve
       NEWCLD = 1.-CLEARA(L)
       SAVCLD = CLDSAV0(L) ! from prev. timestep
#else
c These choices always produce DCLD <= 0.
       NEWCLD = CLDSSL(L)  ! = (1-CLEARA)*FSSL (updated)
       SAVCLD = CLDSAVL(L) ! = (1-CLEARA)      (already updated)
#endif
!@auth Menon for CDNC prediction
#ifdef TRACERS_AEROSOLS_Koch
      CALL GET_CDNC_UPD(L,LHX,WCONST,WMUI,WMX(L),FCLD,NEWCLD,
     *SAVCLD,VVEL,SME(L),DSU,OLDCDL(L),
     *CDNL0,CDNL1)
      OLDCDL(L) = CDNL1
      SNd=CDNL1
C** Pass old and new cloud droplet number
      NEWCDN=SNd
      OLDCDN=CDNL0
c     if (L.eq.1)write(6,*)"BLK_2M NUPD",NEWCDN,OLDCDN
#endif
#ifdef TRACERS_AMP
       OLDCDL(L)=SNd
	   OLDCDI(L)=SNdi
#endif
#endif
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
c Update thermo if environment was changed
       tk0=TL(L)                 ! Temperature, [K]
       qk0=QL(L)                 ! Water vapor mixing ratio, [kq/kq]
       pk0=PL(L)                 ! Pressure, [hPa]
       w0=VVEL*1.d-02           ! Large-scale velocity, [m/s]
       v0=WTURB(L) !; v0=w3d(k,i,j) ! Sub-grid velocity, [m/s]
       r0= 0.0  !RDTNDL(L)               ! T tendency due to radiation, [K/s]
       ldummy=execute_bulk2m_driver('all'
     *           ,tk0,qk0,pk0,w0,v0,r0)
c Update micro if contents were changed
c        mdrop=WMX(L)            ! drop content, [kg water/kg air]
c        ndrop=mdrop/mw0         ! drop concent, [No/m3]
c        mcrys=WMXICE(L)         ! crys content, [kg water/kg air]
c        ncrys=mcrys/mi0         ! crys concent, [No/m3]
      IF(LHX.EQ.LHE)  THEN
         mdrop =WMX(L)
         ndrop= OLDCDL(L)*1.d6  !mdrop/mw0         ! drop concent, [No/m3]
         if(WMX(L).eq.0.) ndrop=0.0
       ELSE
         mcrys =WMX(L)
         WMXICE(L) = WMX(L)
         ncrys= OLDCDI(L)*1.d6  !mcrys/mi0         ! crystal concent, [No/m3]
         if(WMX(L).eq.0.) ncrys=0.0
       ENDIF
c      if(L.eq.1)write(6,*)"5th check BLK_2M",
c    *WMX(L),OLDCDL(L),OLDCDI(L)
c
       ldummy=execute_bulk2m_driver('all'
     *           ,ndrop,mdrop,ncrys,mcrys,'end')
c Get new drop & crys concentration
c     ldummy=execute_bulk2m_driver('surabi','GET_CDNC_UPD',dtB2M,mkx)
#ifdef TRACERS_AEROSOLS_Koch
C*** Call Lohmann's or Gultepe's scheme for CDNC
        OLDCDNC=OLDCDN*1.d6  !convert from cm-3 to m-3
        NEWCDNC=NEWCDN*1.d6  !convert from cm-3 to m-3
	    ldummy=execute_bulk2m_driver('gult','drop_nucl',dtB2M,mkx,
     &OLDCDNC,NEWCDNC)
#endif
#ifdef TRACERS_AMP
C*** Using the AMP_actv interface from MATRIX
        ldummy=execute_bulk2m_driver('matr','drop_nucl',dtB2M,mkx)
#endif
      rablk=execute_bulk2m_driver('get','value','nc') + (
     *  +execute_bulk2m_driver('get','npccn')
     *                                       )*dtB2M
      SNd=rablk(mkx)*1.0d-6            ! ndrop, [No/m^3], SNdL, [No/cc]
c     if(SNd.gt.20.) write(6,*)"Finally out",SNd, l
C**** Old treatment to get CDNC for cloud changes within time steps
       DCLD(L) = NEWCLD-SAVCLD ! cloud fraction change
C** If previous time step is clear sky
       if(SAVCLD.eq.0.) then
         SNd=SNd
       elseif (DCLD(L).le.0.d0) then
         SNd=OLDCDL(L)
       elseif(DCLD(L).gt.0.d0) then
         SNd=( (OLDCDL(L)*SAVCLD) + (SNd*DCLD(L)) )/NEWCLD
        endif
      rablk=execute_bulk2m_driver('get','value','ni') + (
c       nnuccc             ! change n due to contact droplets freez
     * +execute_bulk2m_driver('get','nnuccc')
c       nnucci             ! change n due to immersion droplets freez
     * +execute_bulk2m_driver('get','nnucci')
c       nnuccd             ! change n freezing aerosol (prim ice nuc)
     * +execute_bulk2m_driver('get','nnuccd')
     *                                       )*dtB2M
c      nnucmd              ! change n cond freezing Meyer's (prim ice nuc)
     * +execute_bulk2m_driver('get','nnucmd')
c      nnucmt              ! change n cont freezing Meyer's (prim ice nuc)
     * +execute_bulk2m_driver('get','nnucmt')

c      SNdI=rablk(mkx)*1.0d-6             ! from ncrys [No/m^3] to SNdI in [No/cc]
      SNdI = 0.06417127d0
       if(SNdI.gt.1.d0) SNdI=1.d0      !try to limit to 1000 /l
      OLDCDL(L) = SNd
      OLDCDI(L) = SNdI
#ifdef TRACERS_AMP
       nactc(l,1:nmodes) =  naero(mkx,1:nmodes)
c      do nm=1,nmodes
c        if(nactc(l,nm).gt.0.)
c    *   write(6,*)"NMOD1",nactc(l,nm),l,nm
c       enddo
#endif
c      if(L.eq.1) write(6,*)"6_LO check BLK_2M",SNd,SNdI
c To get effective radii in micron
      rablk=execute_bulk2m_driver('get','value','ec')  ! [micron]
#endif
#ifdef CLD_AER_CDNC
      SCDNCW=SNd
      SNdI = 0.06417127d0
      SCDNCI=SNdI
      If (SCDNCW.le.20.d0) SCDNCW=20.d0   !set min CDNC sensitivity test
c     If (SCDNCI.le.0.06d0) SCDNCI=0.06417127d0   !set min ice crystal
      If (SCDNCI.le.0.0d0) SCDNCI=teeny           !set min ice crystal
      if(SCDNCW.gt.1400.d0) SCDNCw=1400.d0
c     if (SCDNCW.gt.20.) write(6,*) "SCND CDNC",SCDNCW,OLDCDL(l),l
#endif
c      if (debug) print*,"ls0",l,rclde
        IF(LHX.EQ.LHE) THEN

!         RCLD=(RWCLDOX*10.*(1.-PEARTH)+7.0*PEARTH)*(WTEM*4.)**BY3
          RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCW))**BY3
          QHEATC=(QHEAT(L)+FSSL(L)*CLEARA(L)*(EC(L)+ER(L)))/LHX
          IF(RCLD.GT.RWMAX.AND.PREP(L).GT.QHEATC) RCLD=RWMAX
          RCLDE=RCLD/BYBR
c          if (debug) print*,"ls1",l,rclde,rcld,rwmax,wtem,rcldx,
c     *         scdncw
#ifdef CLD_AER_CDNC
C** Using the Liu and Daum paramet
C** for spectral dispersion effects on droplet size distribution
      Repsi=1.d0 - 0.7d0*exp(-0.003d0*SCDNCW)
      Repsis=Repsi*Repsi
      Rbeta=(((1.d0+2.d0*Repsis)**0.667d0))/((1.d0+Repsis)**0.333d0)
!     write(6,*)"RCLD",Rbeta,RCLD,SCDNCW,Repsis
      RCLDE=RCLD*Rbeta
c      if (debug) print*,"ls2",l,rclde
!@auth Menon    end of addition  comment out the RCLDE definition below
#endif
#ifdef BLK_2MOM
c        if(l.eq.1) write(6,*)"7th check BLK_2M",RCLDE,RCLD
c    *   ,SCDNCW,WTEM
!        rablk=execute_bulk2m_driver('get','value','ec')  ! [micron]
!        RCLDE=rablk(mkx)
c        if(l.eq.1) write(6,*)"8th check BLK_2M",RCLDE
#endif
        ELSE
!         RCLD=25.0*(WTEM/4.2d-3)**BY3 * (1.+pl(l)*xRICld)
          RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCI))**BY3
          RCLD=MIN(RCLD,RIMAX)
          RCLDE=RCLD/BYBR
c      if (debug) print*,"ls3",l,rclde

#ifdef BLK_2MOM
c       if(L.eq.1)  write(6,*)"9th check BLK_2M",RCLDE
!        rablk=execute_bulk2m_driver('get','value','ei')  ! [micron]
!        RCLDE=rablk(mkx)
c        if(l.eq.1) write(6,*)"10th check BLK_2M",RCLDE
#endif
        END IF
        RCLDE1=5.*RCLDE          ! for precip optical thickness
        if (cldssl(l).gt.0) CSIZEL(L)=RCLDE
c      if (debug) print*,"lsX",l,csizel(l)
#ifdef CFMIP
        !@auth Mike Bauer
        ! Fix for cases where FCLD is essentially zero for the column and
        ! RCLD is very large (~0.385E+11 microns) this might have been
        ! reduced to RWMAX by IF(RCLD.GT.RWMAX.AND.PREP(L).GT.QHEATC)
        ! except in this case PREP(L) & QHEATC are essentially zero as well.
        IF(FCLD.LE.teeny) THEN
            cfmip_reff(L) = 0.
        ELSE
            cfmip_reff(L) = RCLDE
        ENDIF
#endif
#ifdef CLD_AER_CDNC  /* save for diag purposesi */
        IF (FCLD.gt.1.d-5.and.LHX.eq.LHE) then
            ACDNWS(L)= SCDNCW
            AREWS(L) = RCLDE
            ALWWS(L) = WTEM
            CDN3DL(L)=SCDNCW
            CRE3DL(L)=RCLDE
            NLSW  = NLSW + 1
c      if(ACDNWS(L).gt.20.d0) write(6,*)"INWCLD",ACDNWS(L),
c    * SCDNCW,NLSW,AREWS(L),RCLDE,LHX
        elseif(FCLD.gt.1.d-5.and.LHX.eq.LHS) then
            ACDNIS(L)= SCDNCI
            AREIS(L) = RCLDE
            ALWIS(L) = WTEM
            CDN3DL(L)=SCDNCI
            CRE3DL(L)=RCLDE
            NLSI  = NLSI + 1
c      if(ACDNIS(L).gt.0.d0)    write(6,*)"INICLD",ACDNIS(L),
c    * SCDNCI,NLSI,AREIS(L),RCLDE,LHX
        END IF
#endif
        TEM=AIRM(L)*WMX(L)*1.d2*BYGRAV
        TAUSSL(L)=1.5d3*TEM/(FCLD*RCLDE+teeny)
        TEM1=AIRM(L)*WMPR(L)*1.d2*BYGRAV      ! precip contribution
        TAUSSL(L)=TAUSSL(L)+1.5d3*TEM1/(FCLD*RCLDE1+teeny)
        IF(FCLD.LE.teeny) TAUSSL(L)=0.
        IF(TAUSSL(L).GT.100.) TAUSSL(L)=100.
        IF(LHX.EQ.LHE) WMSUM=WMSUM+TEM      ! pick up water path
#ifdef SCM
        if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
            if (LHX.eq.LHE) SCM_LWP_SS = SCM_LWP_SS + TEM
            if (LHX.eq.LHS) SCM_IWP_SS = SCM_IWP_SS + TEM
        endif
#endif
#ifdef CLD_AER_CDNC
        SMLWP=WMSUM
#endif
      END DO

C**** CALCULATE OPTICAL THICKNESS
      DO L=1,LP50
        CLDSV1(L)=CLDSSL(L)
        IF(WMX(L).LE.0.) SVLHXL(L)=0.
        IF(TAUMCL(L).EQ.0..OR.CKIJ.NE.1.) THEN
          BMAX=1.-EXP(-(CLDSV1(L)/.3d0))
          IF(CLDSV1(L).GE..95d0) BMAX=CLDSV1(L)
          IF(L.EQ.1.OR.L.LE.DCL) THEN
            CLDSSL(L)=MIN(CLDSSL(L)+(BMAX-CLDSSL(L))*CKIJ,FSSL(L))
            TAUSSL(L)=TAUSSL(L)*CLDSV1(L)/(CLDSSL(L)+teeny)
          END IF
          IF(TAUSSL(L).LE.0.) CLDSSL(L)=0.
          IF(L.GT.DCL .AND. TAUMCL(L).LE.0.) THEN
            CLDSSL(L)=MIN(CLDSSL(L)**(2.*BY3),FSSL(L))
            TAUSSL(L)=TAUSSL(L)*CLDSV1(L)**BY3
          END IF
        END IF
        IF(TAUSSL(L).LT.0.) THEN
          WRITE(6,*) 'negative TAUSS CLDSS WMX I J L=',TAUSSL(L),
     *      CLDSSL(L),WMX(L),i_debug,j_debug,L
          TAUSSL(L)=0.
          CLDSSL(L)=0.
          WMX(L)=0.
        END IF
      END DO


#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
!Save variables for 3 hrly diagnostics
c     AAA=1
c     DO L=LP50,1,-1
c       if (CLDSSL(L).gt.0.d0.and.CLDMCL(L).gt.0.d0) then
c        if (TAUSSL(L).gt.0.d0.and.TAUMCL(L).gt.0.d0) then
c         if (CLDSSL(L).LE.randu(xx))GO TO 7
c          AERTAUSS=TAUSSL(L)
c          AERTAU(L)=AERTAUSS
c   7     if (CLDMCL(L).LE.randu(xx))GO TO 8
c          AERTAUMC=TAUMCL(L)
c          IF(AERTAUMC.GT.AERTAU(L)) AERTAU(L)=AERTAUMC
c   8      SUMTAU=SUMTAU+AERTAU(L)
C**
C**   DETECT AT WHAT LAYER SUMTAU (COLUMN OPTICAL THICKNESS)
C**   IS ABOVE THOLD (THRESHOLD SET AT .1 TAU)
C**   THIS LAYER BECOMES CLOUD TOP LAYER (ILTOP)
C**
c         THOLD=(.1)
c         IF(SUMTAU.GT.THOLD)THEN
c         IFLAG=(AAA)
c           IF(IFLAG.EQ.1)THEN
c             ILTOP=L
c             AAA=0
c           ENDIF
c         ENDIF
c        ENDIF
c       ENDIF
c     ENDDO

      DO L=1,LP50
      PRS = (PL(1)-PTOP)/SIG(1)

       IF (L.GE.ls1) THEN
         PPRES = (SIG(L)*(PSF-PTOP)+PTOP)         !in hPa
         DPP= (SIGE(L+1)-SIGE(L))*(PSF-PTOP)      !in hPa
         TEMPR=(TL(L)/PLK(L))*(SIG(L)*(PSF-PTOP)+PTOP)**KAPA
       ELSE
         PPRES= (SIG(L)*PRS+PTOP)
         DPP= (SIGE(L+1)-SIGE(L))*PRS
         TEMPR=(TL(L)/PLK(L))*(SIG(L)*PRS+PTOP)**KAPA
       ENDIF

       CTEML(L)=TEMPR                                        ! Cloud temperature(K)
       D3DL(L)=DPP/PPRES*TEMPR/GRAV*(gasc*1.d03)/mair        ! For Cloud thickness (m)
       IF(CLDSSL(L).GT.0.d0) CD3DL(L)=-1.d0*D3DL(L)*CLDSAVL(L)/CLDSSL(L)
       RHODK=100.d0*PPRES/(gasc*1.d03)/TEMPR*mair
c       IF (SVLHXL(L).EQ.LHE) CL3DL(L) = WMX(L)*RHODK*CD3DL(L) ! cld water kg m-2
c       IF (SVLHXL(L).EQ.LHS) CI3DL(L) = WMX(L)*RHODK*CD3DL(L) ! ice water kg m-2
       IF (SVLHXL(L).EQ.LHE) CL3DL(L) = WMX(L)*RHODK          ! cld water kg m-3
       IF (SVLHXL(L).EQ.LHS) CI3DL(L) = WMX(L)*RHODK          ! ice water kg m-3
c      write(6,*)"CT",L,WMX(L),CD3DL(l),CL3DL(L),CI3DL(L)

c      CTEML(L)=TL(L)                            ! Cloud temperature(K)
c      D3DL(L)=AIRM(L)*TL(L)*bygrav*rgas/PL(L)   ! For Cloud thickness (m)
c      IF(CLDSSL(L).GT.0.d0) CD3DL(L)= 1.d0*D3DL(L)*CLDSAVL(L)/CLDSSL(L)
c      RHO=1d2*PL(L)/(RGAS*TL(L))
c      IF (SVLHXL(L).EQ.LHE) CL3DL(L) = WMX(L)*RHO*CD3DL(L) ! cld water kg m-2
c      IF (SVLHXL(L).EQ.LHS) CI3DL(L) = WMX(L)*RHO*CD3DL(L) ! ice water kg m-2
c      write(6,*)"CT",L,WMX(L),CD3DL(l),CL3DL(L),CI3DL(L)
      END DO
#endif

      RETURN
      END SUBROUTINE LSCOND

c**************************************************************************************

      SUBROUTINE ANVIL_OPTICAL_THICKNESS(SVLATL,                         ! input
     *  RCLDX,MCDNCW,MCDNCI,RIMAX,BYBR,FCLD,TEM,WTEM,                    ! input
     *  RCLD,TAUMC)                                                      ! output

!@sum  MC_OPTICAL_THICKNESS calculates anvil cloud optical thickness and droplet radius
!@auth M.S. Yao

      USE CONSTANT, only :rgas,grav,lhe,lhs,lhm,sha,bysha,pi,by6
     *     ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd,rhow,twopi
      IMPLICIT NONE

!@var SVLATL latent heat
      REAL*8, INTENT(IN) :: SVLATL
!@var RCLDX RCLD multiplier
      REAL*8, INTENT(IN) :: RCLDX
!@var MCDNCW,MCDNCI,RIMAX cloud droplet for liquid and ice, max ice RCLD
      REAL*8, INTENT(IN) :: MCDNCW,MCDNCI,RIMAX
!@var BYBR factor for converting cloud particle radius to effective radius
      REAL*8, INTENT(IN) :: BYBR
!@var RCLD,TAUMC cloud droplet radius and optical thickness
      REAL*8, INTENT(OUT) :: RCLD
      REAL*8, INTENT(INOUT) :: TAUMC
!@var RCLDE effective cloud droplet radius
!@var FCLD,TEM,WTEM cloud cover and dummy variables
      REAL*8 FCLD,TEM,WTEM,RCLDE

C**** CALCULATE ANVIL CLOUD OPTICAL THICKNESS AND DROPLET RADIUS
          IF(SVLATL.EQ.LHE)  THEN
            RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*MCDNCW))**BY3
          ELSE
            RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*MCDNCI))**BY3
            RCLD=MIN(RCLD,RIMAX)
          END IF
          RCLDE=RCLD/BYBR       !  effective droplet radius in anvil
          TAUMC=1.5*TEM/(FCLD*RCLDE+1.E-20)
          IF(TAUMC.GT.100.) TAUMC=100.

      RETURN
      END SUBROUTINE ANVIL_OPTICAL_THICKNESS

c*****************************************************************************************

      SUBROUTINE MC_CLOUD_FRACTION(L,TL,PL,CCM,WCU,PLEMIN,        ! input
     *  PLEL2,PLEMAX,CCMMIN,WCUMIN,LMIN,LMAX,CCMUL,CCMUL2,        ! input
     *  FCLOUD)                                                   ! output

!@sum  MC_CLOUD_FRACTION calculates convective cloud fraction
!@auth M.S. Yao

      USE CONSTANT, only :rgas,grav,lhe,lhs,lhm,sha,bysha,pi,by6
     *     ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd,rhow,twopi
      USE MODEL_COM, only : dtsrc
      IMPLICIT NONE

!@var L index for layer
      INTEGER, INTENT(IN) :: L
!@var TL,PL temperature and pressure
      REAL*8, INTENT(IN) :: TL,PL
!@var CCM,WCU convective mass and cumulus updraft speed
      REAL*8, INTENT(IN) :: CCM,WCU
!@var PLEMIN,PLEL2,PLEMAX pressure at LMIN, L+2 and LMAX+1
      REAL*8, INTENT(IN) :: PLEMIN,PLEL2,PLEMAX
!@var CCMMIN,WCUMIN CCM and WCU at LMIN
      REAL*8, INTENT(IN) :: CCMMIN,WCUMIN
!@var LMIN,LMAX convective base and maximum level
      INTEGER, INTENT(IN) :: LMIN,LMAX
!@var CCMUL,CCMUL2 cloud cover multipliers
      REAL*8, INTENT(IN) :: CCMUL,CCMUL2
!@var FCLOUD cloud fraction
      REAL*8, INTENT(OUT) :: FCLOUD
!@var RHO air density
      REAL*8 RHO

C**** CONVECTIVE CLOUD FRACTION IS CALCULATED BASED ON RATIO OF MASS
C**** FLUX TO UPDRAFT SPEED; AREAL FRACTION BELOW ANVIL CLOUD BASE IS
C**** ASSUMED TO BE TWICE THE AREA OF ACTIVE UPDRAFT CORES
      RHO=PL/(RGAS*TL)
      FCLOUD=CCMUL*CCM/(RHO*GRAV*WCU*DTsrc+teeny)

C**** DEEP ANVIL CLOUD FRACTION IS ASSUMED TO BE 5 TIMES THE NOMINAL
C**** CONVECTIVE CLOUD FRACTION; WILL BE OVERRIDDEN IN LSCOND IF RH
C**** IS ABOVE THRESHOLD FOR STRATIFORM CLOUD FORMATION
      IF(PLEMIN-PLEL2.GE.450.) FCLOUD=5.d0*FCLOUD

C**** "CLOUDY" AREA BELOW CLOUD BASE INCLUDED TO REPRESENT VIRGA
      IF(L.LT.LMIN) THEN
        FCLOUD=CCMUL*CCMMIN/(RHO*GRAV*WCUMIN*DTsrc+teeny)
      END IF

C**** CLOUD FRACTION AT DETRAINMENT LEVEL OF SHALLOW/MIDLEVEL
C**** CONVECTION IS ASSUMED TO BE 3 TIMES NOMINAL CLOUD FRACTION
      IF(PLEMIN-PLEMAX.LT.450.) THEN
        IF(L.EQ.LMAX-1) THEN
          FCLOUD=CCMUL2*CCM/(RHO*GRAV*WCU*DTsrc+teeny)
        END IF

C**** NO VIRGA ASSUMED BELOW SHALLOW/MIDLEVEL CONVECTION
        IF(L.LT.LMIN) FCLOUD=0.
      END IF
      IF(FCLOUD.GT.1.) FCLOUD=1.

      RETURN
      END SUBROUTINE MC_CLOUD_FRACTION

c***************************************************************************************

      SUBROUTINE CONVECTIVE_MICROPHYSICS(PL,WCU,DWCU,LFRZ,WCUFRZ,TP,      ! input
     *  TI,FITMAX,PLAND,CN0,CN0I,CN0G,FLAMW,FLAMG,FLAMI,RHOIP,            ! input
     *  RHOG,ITMAX,TLMIN,TLMIN1,WMAX,                                     ! input
#ifdef CLD_AER_CDNC
     *  TL,RCLD_C,MCDNCW,CONDMU,CONDPC,                                   ! input
#endif
     *  CONDP,CONDP1,CONDIP,CONDGP)                                       ! output

!@sum  CONVECTIVE_MICROPHYSICS calculates convective microphysics and precip condensate
!@auth M.S. Yao

      USE CONSTANT, only :rgas,grav,lhe,lhs,lhm,sha,bysha,pi,by6
     *     ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd,rhow,twopi
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: PL,WCU
!@var PL,WCU pressure and cumulus updraft speed
      INTEGER, INTENT(IN) :: LFRZ,ITMAX
!@var LFRZ,ITMAX freezing level and maximum iterations
      REAL*8, INTENT(IN) :: DWCU,WCUFRZ
!@var DWCU,WCUFRZ increment for WCU and WCU at freezing level
      REAL*8, INTENT(IN) :: TP,TI
!@var TP,TI convective plume termperature and a constant 233.16
      REAL*8, INTENT(IN) :: TLMIN,TLMIN1
!@var TLMIN,TLMIN1 temperature at LMIN and LMIN+1
      REAL*8, INTENT(IN) :: FITMAX
!@var FITMAX 1.d0/ITMAX
      REAL*8, INTENT(IN) :: PLAND,WMAX
!@var PLAND,WMAX land fraction and a constant 50.
      REAL*8, INTENT(IN) :: CN0,CN0I,CN0G
!@var CN0,CN0I,CN0G intercepts of Marshall-Palmer particle size dist.
      REAL*8, INTENT(IN) :: FLAMW,FLAMG,FLAMI
!@var FLAMW,FLAMG,FLAMI Marshall-Palmer lambda for water, graupel and ice
      REAL*8, INTENT(IN) :: RHOIP,RHOG
!@var RHOG,RHOIP density of graupel and ice particles
      REAL*8, INTENT(OUT) :: CONDP,CONDP1
!@var CONDP,CONDP1 precipitating part of convective condensate
      REAL*8 WV,VT
!@var WV,VT convective updraft speed and precip terminal velocity
      REAL*8 FG,FI
!@var FG,FI fraction for graupel and ice
      REAL*8 DCI,DCW,DCG
!@var DCI,DCW,DCG critical cloud particle sizes for onset of precip
      REAL*8 TIG,DDCW
!@var TIG,DDCW dummy variables
      REAL*8 CONDIP,CONDGP
!@var CONDIP,CONDGP precipitating part of convective condensate
      INTEGER ITER
!@var ITER index for iteration
#ifdef CLD_AER_CDNC
      REAL*8 TL,RCLD_C,MCDNCW,CONDMU,CONDPC,RHO
!@var TL temperature
!@var RHO air density
!@var CONDMU,CONDPC convective condensate
!@var RCLD_C,MCDNCW
#endif
C****
C**** CUMULUS MICROPHYSICS IS BASED ON DEL GENIO ET AL. (2005, J. CLIM.).  THE
C**** PARAMETERIZATION ASSUMES A MARSHALL-PALMER PARTICLE SIZE DISTRIBUTION
C**** AND EMPIRICAL SIZE-FALLSPEED RELATIONS FOR LIQUID, GRAUPEL AND FLUFFY ICE.
C**** COMPARISON OF THE CUMULUS UPDRAFT SPEED TO THE FALLSPEEDS OF PARTICLES
C**** OF DIFFERENT SIZES IS USED TO PARTITION THE PARTICLE SIZE DISTRIBUTION
C**** INTO PARTS THAT PRECIPITATE, DETRAIN AT THE CURRENT LEVEL, OR ARE
C**** CARRIED UPWARD TO THE NEXT LEVEL (PUBLISHED VERSION ONLY HAD A 2-WAY
C**** SEPARATION BETWEEN PRECIPITATING AND DETRAINING PARTS).
C****

C**** PARTICLES WITH FALLSPEEDS WITHIN +/- DWCU OF THE UPDRAFT SPEED ARE
C**** ASSUMED TO DETRAIN AT THE CURRENT LEVEL; FIRST FIND SMALL PARTICLE END OF
C**** THIS RANGE (LOWER CRITICAL SIZE) FOR GRAUPEL AND FLUFFY ICE

      WV=WCU-DWCU          ! apply the calculated cumulus updraft speed
      IF(WV.LT.0.d0) WV=0.d0

C**** LOWER CRITICAL SIZE FOR GRAUPEL (EQ. 4B OF DEL GENIO ET AL. (2005))
      DCG=((WV/19.3)*(PL/1000.)**.4)**2.7
      DCG=MIN(DCG,1.D-2)

C**** LOWER CRITICAL SIZE FOR FLUFFY ICE (MODIFIED FROM EQ. 4C OF DEL GENIO ET AL.
C**** (2005) TO REPRESENT UNRIMED RADIATING ASSEMBLAGES RATHER THAN DENDRITES)
      DCI=((WV/11.72)*(PL/1000.)**.4)**2.439
      DCI=MIN(DCI,1.D-2)

C**** ABOVE MELTING LEVEL CONDENSATE IS ASSUMED TO BE FROZEN, TRANSITIONING FROM
C**** PARTLY GRAUPEL TO PURE FLUFFY ICE OVER A RANGE OF TEMPERATURES THAT INCREASES
C**** AS MELTING LEVEL UPDRAFT SPEED STRENGTHENS
      IF(LFRZ.EQ.0) THEN
        TIG=TF           ! initialization
      ELSE
        TIG=TF-4.*WCUFRZ
      END IF
      IF(TIG.LT.TI-10.d0) TIG=TI-10.d0

C**** LOWER CRITICAL SIZE FOR LIQUID WATER (EQ. 4A OF DEL GENIO ET AL. (2005))
      IF (TP.GE.TF) THEN ! water phase
        DDCW=6d-3*FITMAX
        IF(PLAND.LT..5) DDCW=1.5d-3*FITMAX
        DCW=0.
        DO ITER=1,ITMAX-1
          VT=(-.267d0+DCW*(5.15D3-DCW*(1.0225D6-7.55D7*DCW)))*
     *       (1000./PL)**.4d0
          IF(VT.GE.0..AND.ABS(VT-WV).LT..3) EXIT
          IF(VT.GT.WMAX) EXIT
          DCW=DCW+DDCW
        END DO

C**** PART OF LIQUID CONDENSATE NOT ADVECTED UP (EQ. 5, DEL GENIO ET AL. (2005))
        CONDP1 = PRECIP_MP (RHOW,FLAMW,DCW,CN0)

C**** UPPER CRITICAL SIZE FOR LIQUID WATER
        WV=WCU+DWCU          ! apply the calculated cumulus updraft speed
        DCW=0.
        DO ITER=1,ITMAX-1
          VT=(-.267d0+DCW*(5.15D3-DCW*(1.0225D6-7.55D7*DCW)))*
     *       (1000./PL)**.4d0
          IF(VT.GE.0..AND.ABS(VT-WV).LT..3) EXIT
          IF(VT.GT.WMAX) EXIT
          DCW=DCW+DDCW
        END DO

C**** PRECIPITATING PART OF LIQUID CONDENSATE (RAIN)
        CONDP = PRECIP_MP (RHOW,FLAMW,DCW,CN0)

#ifdef CLD_AER_CDNC
        RHO=1d2*PL/(RGAS*TL) !air density in kg/m3
C** Here we calculate condensate converted to precip
C** only if drop size (Reff) exceeds 14 um
C** Conver Rvol to m and CDNC to m from um and cm, respectively
C** Precip condensate is simply the LWC
        CONDPC=4.d0*by3*PI*RHOW*((RCLD_C*1.d-6)**3)*1.d6*MCDNCW/RHO
        IF(CONDMU.gt.CONDPC)  then
          CONDP=CONDMU-CONDPC !
        ELSE
          CONDP=0.d0
        END IF
c     if (CONDP(L).lt.0.d0)
c    *write(6,*)"Mup",CONDP(L),CONDPC(L),CONDMU,DCW,MCDNCW,RCLD_C,L
c      write(6,*)"Mup",CONDP(L),DCW,FLAMW,CONDMU,L
#endif

C**** ABOVE MIXED PHASE REGION, PART OF FLUFFY ICE NOT ADVECTED UP
      ELSE IF (TP.LE.TIG) THEN ! pure ice phase, use TIG
        CONDP1 = PRECIP_MP (RHOIP,FLAMI,DCI,CN0I)

C**** NOW FIND LARGE PARTICLE END OF RANGE FOR DETRAINMENT AT CURRENT
C**** LEVEL FOR EACH FROZEN PARTICLE TYPE (UPPER CRITICAL SIZE)
        WV=WCU+DWCU

C**** UPPER CRITICAL SIZE FOR FLUFFY ICE
        DCI=((WV/11.72)*(PL/1000.)**.4)**2.439
        DCI=MIN(DCI,1.D-2)

C**** PRECIPITATING PART OF FLUFFY ICE (SNOW)
        CONDP = PRECIP_MP (RHOIP,FLAMI,DCI,CN0I)

C**** WITHIN MIXED PHASE REGION, CALCULATE FRACTION OF FROZEN CONDENSATE
C**** THAT IS GRAUPEL VS. FLUFFY ICE
      ELSE ! mixed phase
        FG=(TP-TIG)/((TF-TIG)+teeny)
        IF(FG.GT.1d0) FG=1d0
        IF(FG.LT.0d0) FG=0d0
        IF(TLMIN.LE.TF.OR.TLMIN1.LE.TF) FG=0.d0
        FI=1.-FG

C**** PART OF FROZEN CONDENSATE NOT ADVECTED UP IN MIXED PHASE REGION
        CONDIP = PRECIP_MP (RHOIP,FLAMI,DCI,CN0I)
        CONDGP = PRECIP_MP (RHOG, FLAMG,DCG,CN0G)
        CONDP1=FG*CONDGP+FI*CONDIP
        WV=WCU+DWCU

C**** UPPER CRITICAL SIZE FOR FLUFFY ICE
        DCI=((WV/11.72)*(PL/1000.)**.4)**2.439
        DCI=MIN(DCI,1.D-2)

C**** UPPER CRITICAL SIZE FOR GRAUPEL
        DCG=((WV/19.3)*(PL/1000.)**.4)**2.7
        DCG=MIN(DCG,1.D-2)

C**** PRECIPITATING PART OF FROZEN CONDENSATE (HAIL + SNOW) IN MIXED PHASE REGION
        CONDIP = PRECIP_MP (RHOIP,FLAMI,DCI,CN0I)
        CONDGP = PRECIP_MP (RHOG, FLAMG,DCG,CN0G)
        CONDP=FG*CONDGP+FI*CONDIP
      END IF

      RETURN
      END SUBROUTINE CONVECTIVE_MICROPHYSICS

c**************************************************************************************

      SUBROUTINE MC_PRECIP_PHASE(L,AIRM,FEVAP,LHP1,                     ! input
     *  PRCP,TOLD,TOLD1,VLAT,COND,LMIN,                                 ! input
     *  LHP,MCLOUD,HEAT1)                                               ! output

!@sum  MC_PRECIP_PHASE calculates the convective precip phase
!@auth M.S. Yao

      USE CONSTANT, only :rgas,grav,lhe,lhs,lhm,sha,bysha,pi,by6
     *     ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd,rhow,twopi
      IMPLICIT NONE

!@var L,LMIN index for layer and the convective base layer
      INTEGER, INTENT(IN) :: L,LMIN
!@var AIRM air mass
      REAL*8, INTENT(IN) :: AIRM
!@var FEVAP,LHP1 fraction for evaporation, precip phase at L+1
      REAL*8, INTENT(IN) :: FEVAP,LHP1
!@var TOLD,TOLD1 old temperature at L and L+1
      REAL*8, INTENT(IN) :: TOLD,TOLD1
!@var VLAT,COND save LHX temporarily, condensate
      REAL*8, INTENT(IN) :: VLAT,COND
!@var LHP precip phase
      REAL*8, INTENT(OUT) :: LHP
!@var MCLOUD air mass for re-evaporation of precip
      REAL*8, INTENT(OUT) :: MCLOUD
!@var HEAT1 heating for phase change
      REAL*8, INTENT(OUT) :: HEAT1
!@var PRCP precipitation
      REAL*8, INTENT(INOUT) :: PRCP

C**** DETERMINE THE AIR MASS FOR RAIN EVAPORATION TO TAKE PLACE
      MCLOUD=0.
      IF(L.LE.LMIN) MCLOUD=2.*FEVAP*AIRM
      IF(MCLOUD.GT.AIRM) MCLOUD=AIRM
      lhp=lhp1

C**** phase of precip is defined by original environment (i.e. TOLD)
C**** phase of condensate is defined by VLAT (based on plume T).

C**** decide whether to melt frozen precip
      IF (lhp.eq.lhs .and. TOLD.GT.TF.AND.TOLD1.LE.TF) then
c         if (debug) print*,"cnv1",l,lhp,told,told1,LHM*PRCP*BYSHA
        HEAT1=HEAT1+LHM*PRCP*BYSHA
        lhp=lhe
      end if
C**** and deal with possible inversions and re-freezing of rain
      IF (lhp.eq.lhe .and. TOLD.LE.TF.AND.TOLD1.GT.TF) then
c         if (debug) print*,"cnv2",l,lhp,told,told1,-LHM*PRCP*BYSHA
        HEAT1=HEAT1-LHM*PRCP*BYSHA
        lhp=lhs
      end if
C**** check for possible inconsistency with cond at this level
      IF (LHP.NE.VLAT.AND. COND.GT.0) then ! convert to cond to precip phase
        HEAT1=HEAT1+(VLAT-LHP)*COND*BYSHA
c        if (debug) print*,"cnv3",l,lhp,vlat,cond,(VLAT-LHP
c     $       )*COND*BYSHA
      END IF

      RETURN
      END SUBROUTINE MC_PRECIP_PHASE

c***************************************************************************************

      FUNCTION  PRECIP_MP (RHO, FLAM, DC, CN)
C
c     Calculates the mass density of the precipitating condensate part of the
C     Marshall-Palmer distribution for specified phase (liquid,graupel,ice)
C
C     INPUTS:
C        RHO  -  Density pf Condenssed Phase
C        FLAM -  Marshall-Palmer Slope Parameter
C        DC   -  Critical Size for Precipitation
C        CN   -  Marshall-Palmer Intercept
C
      USE CONSTANT, only : PI, BY6
      IMPLICIT NONE
C
      REAL*8   PRECIP_MP, RHO, FLAM, DC, CN
C
      PRECIP_MP = RHO*(PI*BY6)*CN*EXP(-FLAM*DC)*
     *            (DC*DC*DC/FLAM+3.*DC*DC/(FLAM*FLAM)+
     *            6.*DC/(FLAM*FLAM*FLAM)+6./FLAM**4.)
C
      RETURN
      END FUNCTION PRECIP_MP

C***************************************************************************************

      SUBROUTINE MASS_FLUX ( FPLUME, FMP2, DQSUM, LMIN,
     *           LHX, QMO1, QMO2, SLH, SMO1, SMO2, WMDN, WMUP, WMEDG )
C
      USE CONSTANT, only : deltx
C         deltx = coeff. of humidity in virtual temperature defn. (0.6078)
C
      IMPLICIT NONE
C
C ****************************************************************************
C * MASS_FLUX  - PERFORMS ITERATION TO FIND FPLUME MASS WHICH RESTORES CLOUD *
C *              BASE TO NEUTRAL STATE                                       *
C ****************************************************************************
C
C **************
C **  INPUTS  **
C **************
C        REAL*8
C             AIRM                        the layer's pressure depth (mb)
C             BYAM                        1./AIRM
C             PL                          layer pressure (mb)
C             PLK                         PL**KAPA
C             SM,QM                       Vertical profiles of T/Q
         INTEGER
     *        LMIN                     !  base layer of a convective event
         REAL*8
     *        LHX,                     !  latent heat of evaporation (J/Kg)
     *        QMO1,QMO2, SMO1,SMO2,    !  LMIN & LMIN+1 elements of QM & SM
     *        SLH,                     !  LHX/SHA
     *        WMDN,WMUP,WMEDG          !  LMIN & LMIN+1 elements of WML & middle
C
C ***************
C **  OUTPUTS  **
C ***************
         REAL*8
     *        FPLUME,                  !  fraction of convective plume
     *        FMP2,                    !  FPLUME*AIRM(LMIN)
     *        DQSUM                    !  SUM OF DQs DURING ITERATION
C
C **************
C **  LOCAL   **
C **************
         INTEGER ITER        !   number for iteration
         REAL*8
     *        DMSE1,         !   difference in moist static energy
     *        DFP,           !   an iterative increment
     *        FEVAP,         !   fraction of air available for precip evaporation
     *        MCLOUD,        !   air mass available for re-evaporation of precip
     *        SMP, QMP,      !   plume's SM, QM
     *        TP,            !   plume's temperature (K)
     *        QSATC,         !   saturation vapor mixing ratio
C             QNX,TNX        !   Dummy Variables in /CLDPRV/
     *        FRAT1,FRAT2,QMN1,QMN2,SMN1,SMN2,   !   Dummy Variables
     *        GAMA,DQ,QDN,QUP,SDN,SUP,SEDGE,     !   Dummy Variables
     *        SVDN,SVUP,SVEDG,QEDGE              !   Dummy Variables
C
C *************
C   Functions
C     QSATMP  plume's saturation vapor mixing ratio (function)
C     QSAT    saturation specific humidity
C     DQSATDT dQSAT/dT
C     THBAR   virtual temperature at layer edge
      REAL*8 QSATMP, QSAT, DQSATDT, THBAR
C
C
      FPLUME=.25
      DFP = .25
      DO ITER=1,8
        DFP=DFP*0.5
        FMP2=FPLUME*AIRM(LMIN)
        FRAT1=FMP2*BYAM(LMIN+1)
        FRAT2=FMP2*BYAM(LMIN+2)
        SMN1=SMO1*(1.-FPLUME)+FRAT1*SMO2
        QMN1=QMO1*(1.-FPLUME)+FRAT1*QMO2
        SMN2=SMO2*(1.-FRAT1)+FRAT2*SM(LMIN+2)
        QMN2=QMO2*(1.-FRAT1)+FRAT2*QM(LMIN+2)
        SMP=SMO1*FPLUME
        QMP=QMO1*FPLUME
        TP=SMO1*PLK(LMIN+1)*BYAM(LMIN)
        QSATMP=FMP2*QSAT(TP,LHX,PL(LMIN+1))
        GAMA=SLH*QSATMP*DQSATDT(TP,LHX)/FMP2
        DQSUM=(QMP-QSATMP)/(1.+GAMA)
        IF(DQSUM.GT.0.)  THEN
           FEVAP=.5*FPLUME
           MCLOUD=FEVAP*AIRM(LMIN+1)
           TNX=SMO2*PLK(LMIN+1)*BYAM(LMIN+1)
           QNX=QMO2*BYAM(LMIN+1)
           QSATC=QSAT(TNX,LHX,PL(LMIN+1))
           DQ=MCLOUD*(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
           IF(DQ.GT.DQSUM) DQ=DQSUM
           SMN2=SMN2-SLH*DQ/PLK(LMIN+1)
           QMN2=QMN2+DQ
           DQSUM=DQSUM-DQ
           IF(DQSUM.GT.0.)  THEN
              MCLOUD=FEVAP*AIRM(LMIN)
              TNX=SMO1*PLK(LMIN)*BYAM(LMIN)
              QNX=QMO1*BYAM(LMIN)
              QSATC=QSAT(TNX,LHX,PL(LMIN))
              DQ=MCLOUD*(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
              IF(DQ.GT.DQSUM) DQ=DQSUM
              SMN1=SMN1-SLH*DQ/PLK(LMIN)
              QMN1=QMN1+DQ
           END IF
        END IF
        SDN=SMN1*BYAM(LMIN)
        SUP=SMN2*BYAM(LMIN+1)
        SEDGE=THBAR(SUP,SDN)
        QDN=QMN1*BYAM(LMIN)
        QUP=QMN2*BYAM(LMIN+1)
        SVDN=SDN*(1.+DELTX*QDN-WMDN)
        SVUP=SUP*(1.+DELTX*QUP-WMUP)
        QEDGE=.5*(QUP+QDN)
        SVEDG=SEDGE*(1.+DELTX*QEDGE-WMEDG)
        DMSE1=(SVUP-SVEDG)*PLK(LMIN+1)+(SVEDG-SVDN)*PLK(LMIN)+
     *       SLHE*(QSAT(SUP*PLK(LMIN+1),LHX,PL(LMIN+1))-QDN)
        IF (ABS(DMSE1).LE.1.d-3) EXIT ! finish iteration
        IF(DMSE1.GT. 1.d-3) FPLUME=FPLUME-DFP
        IF(DMSE1.LT.-1.d-3) FPLUME=FPLUME+DFP
      END DO
C
      RETURN
C
      END SUBROUTINE MASS_FLUX

C****************************************************************************************

      END MODULE CLOUDS

C----------

      SUBROUTINE ISCCP_CLOUD_TYPES(sunlit,pfull
     *     ,phalf,qv,cc,conv,dtau_s,dtau_c,skt,at,dem_s,dem_c,itrop
     *     ,fq_isccp,meanptop,meantaucld,boxtau,boxptop,nbox,jerr)
!@sum  ISCCP_CLOUD_TYPES calculate isccp cloud diagnostics in a column
!@auth Gavin Schmidt
!@ver  2.0 (from isccp version 3.5)
!@     GISS modifications:
!@       1) used modelE random number generator
!@       2) real-->real*8, e-2 -> d-2 etc.
!@       3) ncol, nlev taken from modules, add jerr flag
!@       4) emsfc_lw,top_height,overlap,npoints fixed
!@       5) sub in for constants
!@       6) use pre-calculated tropopause
!@       7) tautab/invtau from module
!@       8) removed boxtau,boxptop from output
!@       9) added back nbox for backwards compatibility
!@@@@@ added boxtau,boxptop back to output      Audrey Wolf
! *****************************COPYRIGHT*******************************
! (c) COPYRIGHT Steve Klein and Mark Webb 2004, All Rights Reserved.
! Steve Klein klein21@mail.llnl.gov
! Mark Webb mark.webb@metoffice.gov.uk markwebb@mail.com
! ISCCP SIMULATOR icarus-scops version 3.5
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
! See also http://www.gnu.org/copyleft/lesser.txt
!
! The Met Office hereby disclaims all copyright interest in
! the ISCCP Simulator ( a library which converts model clouds into
! direct equivalents to the satellite observations of the ISCCP)
! written by Steve Klein and Mark Webb
!
! Catherine Senior, August 2004
! Manager Climate Sensitivy Group
! Met Office   Hadley Centre for Climate Prediction and Research
! FitzRoy Road  Exeter  EX1 3PB  United Kingdom
! *****************************COPYRIGHT*******************************
      USE CONSTANT, only : wtmair=>mair,Navo=>avog,bygrav,bymrat
      USE RANDOM, only : randu
      USE MODEL_COM, only : nlev=>lm,qcheck
      USE CLOUDS, only : tautab,invtau
      USE CLOUDS_COM, only : ncol
      implicit none
!@var  emsfc_lw    longwave emissivity of surface at 10.5 microns
      REAL*8, PARAMETER :: emsfc_lw=0.99d0
!@var pc1bylam Planck constant c1 by wavelength (10.5 microns)
      REAL*8, PARAMETER :: pc1bylam = 1.439d0/10.5d-4
!@var boxarea fractional area of each sub-grid scale box
      REAL*8  :: boxarea
!@var t0 ave temp  (K)
      REAL*8, PARAMETER :: t0 = 296.
!@var t0bypstd ave temp by sea level press
      REAL*8, PARAMETER :: t0bypstd = t0/1.013250d6
      real*8, parameter :: bywc = 1./2.56d0 , byic= 1./2.13d0
      real*8, parameter :: isccp_taumin = 0.3d0  ! default is 0.3, used to be 0.1

!  number of model points in the horizontal (FIXED FOR ModelE CODE)
      INTEGER, PARAMETER :: npoints=1

!     -----
!     Input
!     -----

c      INTEGER nlev          !  number of model levels in column
c      INTEGER ncol          !  number of subcolumns

      INTEGER sunlit(npoints) !  1 for day points, 0 for night time

c      INTEGER seed(npoints)
      !  seed values for marsaglia  random number generator
      !  It is recommended that the seed is set
      !  to a different value for each model
      !  gridbox it is called on, as it is
      !  possible that the choice of the same
      !  seed value every time may introduce some
      !  statistical bias in the results, particularly
      !  for low values of NCOL.

      REAL*8 pfull(npoints,nlev)
                       !  pressure of full model levels (Pascals)
                  !  pfull(npoints,1) is top level of model
                  !  pfull(npoints,nlev) is bot of model

      REAL*8 phalf(npoints,nlev+1)
                  !  pressure of half model levels (Pascals)
                  !  phalf(npoints,1) is top of model
                  !  phalf(npoints,nlev+1) is the surface pressure

      REAL*8 qv(npoints,nlev)
                  !  water vapor specific humidity (kg vapor/ kg air)
                  !         on full model levels

      REAL*8 cc(npoints,nlev)
                  !  input cloud cover in each model level (fraction)
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by clouds

      REAL*8 conv(npoints,nlev)
                  !  input convective cloud cover in each model
                  !   level (fraction)
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by convective clouds

      REAL*8 dtau_s(npoints,nlev)
                  !  mean 0.67 micron optical depth of stratiform
                !  clouds in each model level
                  !  NOTE:  this the cloud optical depth of only the
                  !  cloudy part of the grid box, it is not weighted
                  !  with the 0 cloud optical depth of the clear
                  !         part of the grid box

      REAL*8 dtau_c(npoints,nlev)
                  !  mean 0.67 micron optical depth of convective
                !  clouds in each
                  !  model level.  Same note applies as in dtau_s.

c      INTEGER overlap           !  overlap type
                                !  1=max
                                !  2=rand
                                !  3=max/rand

c      INTEGER top_height        !  1 = adjust top height using both a computed
                                !  infrared brightness temperature and the visible
                                !  optical depth to adjust cloud top pressure. Note
                                !  that this calculation is most appropriate to compare
                                !  to ISCCP data during sunlit hours.
                                !  2 = do not adjust top height, that is cloud top
                                !  pressure is the actual cloud top pressure
                                !  in the model
                                !  3 = adjust top height using only the computed
                                !  infrared brightness temperature. Note that this
                                !  calculation is most appropriate to compare to ISCCP
                                !  IR only algortihm (i.e. you can compare to nighttime
                                !  ISCCP data with this option)
      INTEGER, PARAMETER :: top_height=1, overlap=3


c      REAL*8 tautab(0:255)      !  ISCCP table for converting count value to
                                !  optical thickness

c      INTEGER invtau(-20:45000) !  ISCCP table for converting optical thickness
                                !  to count value
!
!     The following input variables are used only if top_height = 1 or top_height = 3
!
      REAL*8 skt(npoints)       !  skin Temperature (K)
c      REAL*8 emsfc_lw           !  10.5 micron emissivity of surface (fraction)

      REAL*8 at(npoints,nlev)   !  temperature in each model level (K)
      REAL*8 dem_s(npoints,nlev) !  10.5 micron longwave emissivity of stratiform
                                !  clouds in each
                                !  model level.  Same note applies as in dtau_s.
      REAL*8 dem_c(npoints,nlev) !  10.5 micron longwave emissivity of convective
                                !  clouds in each
                                !  model level.  Same note applies as in dtau_s.
      INTEGER itrop(npoints)    ! model level containing tropopause (WMO defn)

!     ------
!     Output
!     ------

      REAL*8 fq_isccp(npoints,7,7) !  the fraction of the model grid box covered by
                                !  each of the 49 ISCCP D level cloud types

      REAL*8 totalcldarea(npoints) !  the fraction of model grid box columns
                                !  with cloud somewhere in them.  This should
                                !  equal the sum over all entries of fq_isccp
      INTEGER nbox(npoints) !  the number of model grid box columns
                                !  with cloud somewhere in them.

! The following three means are averages over the cloudy areas only.  If no
! clouds are in grid box all three quantities should equal zero.

      REAL*8 meanptop(npoints)  !  mean cloud top pressure (mb) - linear averaging
                                !  in cloud top pressure.

      REAL*8 meantaucld(npoints) !  mean optical thickness
                                !  linear averaging in albedo performed.

      REAL*8 boxtau(npoints,ncol) !  optical thickness in each column

      REAL*8 boxptop(npoints,ncol) !  cloud top pressure (mb) in each column

      INTEGER JERR  ! error flag

!
!     ------
!     Working variables added when program updated to mimic Mark Webb's PV-Wave code
!     ------

      REAL*8 frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
                              ! Equivalent of BOX in original version, but
                              ! indexed by column then row, rather than
                              ! by row then column

      REAL*8 tca(npoints,0:nlev) ! total cloud cover in each model level (fraction)
                                ! with extra layer of zeroes on top
                                ! in this version this just contains the values input
                                ! from cc but with an extra level
      REAL*8 cca(npoints,nlev)  ! convective cloud cover in each model level (fraction)
                                ! from conv

      REAL*8 threshold(npoints,ncol) ! pointer to position in gridbox
      REAL*8 maxocc(npoints,ncol) ! Flag for max overlapped conv cld
      REAL*8 maxosc(npoints,ncol) ! Flag for max overlapped strat cld

      REAL*8 boxpos(npoints,ncol) ! ordered pointer to position in gridbox

      REAL*8 threshold_min(npoints,ncol) ! minimum value to define range in with new threshold
                                ! is chosen

      REAL*8 dem(npoints,ncol),bb(npoints) !  working variables for 10.5 micron longwave
                                !  emissivity in part of
                                !  gridbox under consideration

      REAL*8 ran(npoints)       ! vector of random numbers
      REAL*8 ptrop(npoints)
      REAL*8 attrop(npoints)
c      REAL*8 attropmin (npoints)
      REAL*8 atmax(npoints)
      REAL*8 atmin(npoints)
      REAL*8 btcmin(npoints)
      REAL*8 transmax(npoints)

      INTEGER i,j,ilev,ibox
      INTEGER ipres(npoints)
      INTEGER itau(npoints),ilev2
      INTEGER acc(nlev,ncol)
      INTEGER match(npoints,nlev-1)
      INTEGER nmatch(npoints)
      INTEGER levmatch(npoints,ncol)

      !variables needed for water vapor continuum absorption
      real*8 fluxtop_clrsky(npoints),trans_layers_above_clrsky(npoints)
      real*8 taumin(npoints)
      real*8 dem_wv(npoints,nlev)
      real*8 press(npoints), dpress(npoints), atmden(npoints)
      real*8 rvh20(npoints), wk(npoints), rhoave(npoints)
      real*8 rh20s(npoints), rfrgn(npoints)
      real*8 tmpexp(npoints),tauwv(npoints)

      character*1 cchar(6),cchar_realtops(6)
      integer icycle
      REAL*8 tau(npoints,ncol)
      LOGICAL box_cloudy(npoints,ncol)
      REAL*8 tb(npoints,ncol)
      REAL*8 ptop(npoints,ncol)
      REAL*8 emcld(npoints,ncol)
      REAL*8 fluxtop(npoints,ncol)
      REAL*8 trans_layers_above(npoints,ncol)
      real*8 fluxtopinit(npoints),tauir(npoints)
      real*8 meanalbedocld(npoints)
      REAL*8 albedocld(npoints,ncol)
c      integer debug       ! set to non-zero value to print out inputs
c                    ! with step debug
c      integer debugcol    ! set to non-zero value to print out column
c                    ! decomposition with step debugcol
      integer rangevec(npoints),rangeerror

      real*8 tauchk,xx

c      character*10 ftn09

      DATA cchar / ' ','-','1','+','I','+'/
      DATA cchar_realtops / ' ',' ','1','1','I','I'/

      JERR=0
      boxarea = 1d0/ncol
c      INTEGER irand,i2_16,huge32,overflow_32  ! variables for RNG
c      PARAMETER(huge32=2147483647)
c      i2_16=65536

      tauchk = -1.*log(0.9999999)

c      ncolprint=0
c
c      if ( debug.ne.0 ) then
c          j=1
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write(6,'(a10)') 'debug='
c          write(6,'(8I10)') debug
c          write(6,'(a10)') 'debugcol='
c          write(6,'(8I10)') debugcol
c          write(6,'(a10)') 'npoints='
c          write(6,'(8I10)') npoints
c          write(6,'(a10)') 'nlev='
c          write(6,'(8I10)') nlev
c          write(6,'(a10)') 'ncol='
c          write(6,'(8I10)') ncol
c          write(6,'(a10)') 'top_height='
c          write(6,'(8I10)') top_height
c          write(6,'(a10)') 'overlap='
c          write(6,'(8I10)') overlap
c          write(6,'(a10)') 'emsfc_lw='
c          write(6,'(8f10.2)') emsfc_lw
c          write(6,'(a10)') 'tautab='
c          write(6,'(8f10.2)') tautab
c          write(6,'(a10)') 'invtau(1:100)='
c          write(6,'(8i10)') (invtau(i),i=1,100)
c        do j=1,npoints,debug
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write(6,'(a10)') 'sunlit='
c          write(6,'(8I10)') sunlit(j)
c          write(6,'(a10)') 'seed='
c          write(6,'(8I10)') seed(j)
c          write(6,'(a10)') 'pfull='
c          write(6,'(8f10.2)') (pfull(j,i),i=1,nlev)
c          write(6,'(a10)') 'phalf='
c          write(6,'(8f10.2)') (phalf(j,i),i=1,nlev+1)
c          write(6,'(a10)') 'qv='
c          write(6,'(8f10.3)') (qv(j,i),i=1,nlev)
c          write(6,'(a10)') 'cc='
c          write(6,'(8f10.3)') (cc(j,i),i=1,nlev)
c          write(6,'(a10)') 'conv='
c          write(6,'(8f10.2)') (conv(j,i),i=1,nlev)
c          write(6,'(a10)') 'dtau_s='
c          write(6,'(8g12.5)') (dtau_s(j,i),i=1,nlev)
c          write(6,'(a10)') 'dtau_c='
c          write(6,'(8f10.2)') (dtau_c(j,i),i=1,nlev)
c          write(6,'(a10)') 'skt='
c          write(6,'(8f10.2)') skt(j)
c          write(6,'(a10)') 'at='
c          write(6,'(8f10.2)') (at(j,i),i=1,nlev)
c          write(6,'(a10)') 'dem_s='
c          write(6,'(8f10.3)') (dem_s(j,i),i=1,nlev)
c          write(6,'(a10)') 'dem_c='
c          write(6,'(8f10.3)') (dem_c(j,i),i=1,nlev)
c        enddo
c      endif

!     ---------------------------------------------------!

!     assign 2d tca array using 1d input array cc

      do j=1,npoints
        tca(j,0)=0
      enddo

      do ilev=1,nlev
        do j=1,npoints
          tca(j,ilev)=cc(j,ilev)
        enddo
      enddo

!     assign 2d cca array using 1d input array conv

      do ilev=1,nlev
          do j=1,npoints
           cca(j,ilev)=conv(j,ilev)
          enddo
      enddo

c      if (ncolprint.ne.0) then
c      do j=1,npoints,1000
c        write(6,'(a10)') 'j='
c        write(6,'(8I10)') j
c        write (6,'(a)') 'seed:'
c        write (6,'(I3.2)') seed(j)
c
c        write (6,'(a)') 'tca_pp_rev:'
c        write (6,'(8f5.2)')
c     &   ((tca(j,ilev)),
c     &      ilev=1,nlev)
c
c        write (6,'(a)') 'cca_pp_rev:'
c        write (6,'(8f5.2)')
c     &   ((cca(j,ilev),ibox=1,ncolprint),ilev=1,nlev)
c      enddo
c      endif

      if (top_height .eq. 1 .or. top_height .eq. 3) then

      do j=1,npoints ! use pre-computed tropopause level
        ptrop(j) = pfull(j,itrop(j))
        attrop(j) = at(j,itrop(j))
        atmin(j) = 400.
        atmax(j) = 0.
      enddo

      do 12 ilev=1,nlev
        do j=1,npoints
          if (at(j,ilev) .gt. atmax(j)) atmax(j)=at(j,ilev)
          if (at(j,ilev) .lt. atmin(j)) atmin(j)=at(j,ilev)
        enddo
 12   continue

      end if

!     -----------------------------------------------------!

!     ---------------------------------------------------!

      do ilev=1,nlev
        do j=1,npoints

          rangevec(j)=0

          if (cc(j,ilev) .lt. 0. .or. cc(j,ilev) .gt. 1.) then
!           error = cloud fraction less than zero
!           error = cloud fraction greater than 1
            rangevec(j)=rangevec(j)+1
          endif

          if (conv(j,ilev) .lt. 0. .or. conv(j,ilev) .gt. 1.) then
!           ' error = convective cloud fraction less than zero'
!           ' error = convective cloud fraction greater than 1'
            rangevec(j)=rangevec(j)+2
          endif

          if (dtau_s(j,ilev) .lt. 0.) then
!           ' error = stratiform cloud opt. depth less than zero'
            rangevec(j)=rangevec(j)+4
          endif

          if (dtau_c(j,ilev) .lt. 0.) then
!           ' error = convective cloud opt. depth less than zero'
            rangevec(j)=rangevec(j)+8
          endif

          if (dem_s(j,ilev) .lt. 0. .or. dem_s(j,ilev) .gt. 1.) then
!             ' error = stratiform cloud emissivity less than zero'
!             ' error = stratiform cloud emissivity greater than 1'
            rangevec(j)=rangevec(j)+16
          endif

          if (dem_c(j,ilev) .lt. 0. .or. dem_c(j,ilev) .gt. 1.) then
!             ' error = convective cloud emissivity less than zero'
!             ' error = convective cloud emissivity greater than 1'
              rangevec(j)=rangevec(j)+32
          endif
        enddo

        rangeerror=0
        do j=1,npoints
            rangeerror=rangeerror+rangevec(j)
        enddo

        if (rangeerror.ne.0) then
              write (6,*) 'Input variable out of range'
              write (6,*) 'rangevec:'
              write (6,*) rangevec
              call sys_flush(6)
              JERR=1 ; return  !STOP
        endif
      enddo

      do ibox=1,ncol
        do j=1,npoints
        boxpos(j,ibox)=(ibox-.5)*boxarea
        enddo
      enddo

!     ---------------------------------------------------!
!     Initialise working variables
!     ---------------------------------------------------!

!     Initialised frac_out to zero

      do ilev=1,nlev
        do ibox=1,ncol
          do j=1,npoints
          frac_out(j,ibox,ilev)=0.0
          enddo
        enddo
      enddo

c      if (ncolprint.ne.0) then
c        write (6,'(a)') 'frac_out_pp_rev:'
c          do j=1,npoints,1000
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write (6,'(8f5.2)')
c     &     ((frac_out(j,ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)
c
c          enddo
c        write (6,'(a)') 'ncol:'
c        write (6,'(I3)') ncol
c      endif
c      if (ncolprint.ne.0) then
c        write (6,'(a)') 'last_frac_pp:'
c          do j=1,npoints,1000
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write (6,'(8f5.2)') (tca(j,0))
c          enddo
c      endif

!     ---------------------------------------------------!
!     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
!     frac_out is the array that contains the information
!     where 0 is no cloud, 1 is a stratiform cloud and 2 is a
!     convective cloud

      !loop over vertical levels
      DO 200 ilev = 1,nlev

!     Initialise threshold

        IF (ilev.eq.1) then
          ! If max overlap
          IF (overlap.eq.1) then
            ! select pixels spread evenly
            ! across the gridbox
              DO ibox=1,ncol
                do j=1,npoints
                  threshold(j,ibox)=boxpos(j,ibox)
                enddo
              enddo
          ELSE
              DO ibox=1,ncol
c                include 'congvec.f'
          do j=1,npoints
            ran(j)=randu(xx)
          end do
                ! select random pixels from the non-convective
                ! part the gridbox ( some will be converted into
                ! convective pixels below )
                do j=1,npoints
                  threshold(j,ibox)=
     &            cca(j,ilev)+(1-cca(j,ilev))*ran(j)
                enddo
              enddo
            END IF
c            IF (ncolprint.ne.0) then
c              write (6,'(a)') 'threshold_nsf2:'
c                do j=1,npoints,1000
c                write(6,'(a10)') 'j='
c                write(6,'(8I10)') j
c                write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
c                enddo
c            END IF
        END IF

c        IF (ncolprint.ne.0) then
c            write (6,'(a)') 'ilev:'
c            write (6,'(I2)') ilev
c        END IF

        DO ibox=1,ncol

          ! All versions
          do j=1,npoints
            if (boxpos(j,ibox).le.cca(j,ilev)) then
              maxocc(j,ibox) = 1
            else
              maxocc(j,ibox) = 0
            end if
          enddo

          ! Max overlap
          if (overlap.eq.1) then
            do j=1,npoints
              threshold_min(j,ibox)=cca(j,ilev)
              maxosc(j,ibox)=1
            enddo
          endif

          ! Random overlap
          if (overlap.eq.2) then
            do j=1,npoints
              threshold_min(j,ibox)=cca(j,ilev)
              maxosc(j,ibox)=0
            enddo
          endif

          ! Max/Random overlap
          if (overlap.eq.3) then
            do j=1,npoints
              threshold_min(j,ibox)=max(cca(j,ilev),
     &          min(tca(j,ilev-1),tca(j,ilev)))
              if (threshold(j,ibox)
     &          .lt.min(tca(j,ilev-1),tca(j,ilev))
     &          .and.(threshold(j,ibox).gt.cca(j,ilev))) then
                   maxosc(j,ibox)= 1
              else
                   maxosc(j,ibox)= 0
              end if
            enddo
          endif

          ! Reset threshold

c          include 'congvec.f'
          do j=1,npoints
            ran(j)=randu(xx)
          end do

          do j=1,npoints
            threshold(j,ibox)=
              !if max overlapped conv cloud
     &        maxocc(j,ibox) * (
     &            boxpos(j,ibox)
     &        ) +
              !else
     &        (1-maxocc(j,ibox)) * (
                  !if max overlapped strat cloud
     &            (maxosc(j,ibox)) * (
                      !threshold=boxpos
     &                threshold(j,ibox)
     &            ) +
                  !else
     &            (1-maxosc(j,ibox)) * (
                      !threshold_min=random[thrmin,1]
     &                threshold_min(j,ibox)+
     &                  (1-threshold_min(j,ibox))*ran(j)
     &           )
     &        )
          enddo

        END DO ! ibox

!          Fill frac_out with 1's where tca is greater than the threshold

           DO ibox=1,ncol
             do j=1,npoints
               if (tca(j,ilev).gt.threshold(j,ibox)) then
               frac_out(j,ibox,ilev)=1
               else
               frac_out(j,ibox,ilev)=0
               end if
             enddo
           END DO

!         Code to partition boxes into startiform and convective parts
!         goes here

           DO ibox=1,ncol
             do j=1,npoints
                if (threshold(j,ibox).le.cca(j,ilev)) then
                    ! = 2 IF threshold le cca(j)
                    frac_out(j,ibox,ilev) = 2
                else
                    ! = the same IF NOT threshold le cca(j)
                    frac_out(j,ibox,ilev) = frac_out(j,ibox,ilev)
                end if
             enddo
           END DO

!         Set last_frac to tca at this level, so as to be tca
!         from last level next time round

c          if (ncolprint.ne.0) then
c
c            do j=1,npoints ,1000
c            write(6,'(a10)') 'j='
c            write(6,'(8I10)') j
c            write (6,'(a)') 'last_frac:'
c            write (6,'(8f5.2)') (tca(j,ilev-1))
c
c            write (6,'(a)') 'cca:'
c            write (6,'(8f5.2)') (cca(j,ilev),ibox=1,ncolprint)
c
c            write (6,'(a)') 'max_overlap_cc:'
c            write (6,'(8f5.2)') (maxocc(j,ibox),ibox=1,ncolprint)
c
c            write (6,'(a)') 'max_overlap_sc:'
c            write (6,'(8f5.2)') (maxosc(j,ibox),ibox=1,ncolprint)
c
c            write (6,'(a)') 'threshold_min_nsf2:'
c            write (6,'(8f5.2)') (threshold_min(j,ibox),ibox=1,ncolprint)
c
c            write (6,'(a)') 'threshold_nsf2:'
c            write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
c
c            write (6,'(a)') 'frac_out_pp_rev:'
c            write (6,'(8f5.2)')
c     &       ((frac_out(j,ibox,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
c          enddo
c          endif

200   CONTINUE    !loop over nlev

!
!     ---------------------------------------------------!


!
!     ---------------------------------------------------!
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau

      !initialize tau and albedocld to zero
      do 15 ibox=1,ncol
        do j=1,npoints
            tau(j,ibox)=0.
          albedocld(j,ibox)=0.
          boxtau(j,ibox)=0.
          boxptop(j,ibox)=0.
          box_cloudy(j,ibox)=.false.
        enddo
15    continue

      !compute total cloud optical depth for each column
      do ilev=1,nlev
            !increment tau for each of the boxes
            do ibox=1,ncol
              do j=1,npoints
                 if (frac_out(j,ibox,ilev).eq.1) then
                        tau(j,ibox)=tau(j,ibox)
     &                     + dtau_s(j,ilev)
                 endif
                 if (frac_out(j,ibox,ilev).eq.2) then
                        tau(j,ibox)=tau(j,ibox)
     &                     + dtau_c(j,ilev)
                 end if
              enddo
            enddo ! ibox
      enddo ! ilev
c          if (ncolprint.ne.0) then
c
c              do j=1,npoints ,1000
c                write(6,'(a10)') 'j='
c                write(6,'(8I10)') j
c                write(6,'(i2,1X,8(f7.2,1X))')
c     &          ilev,
c     &          (tau(j,ibox),ibox=1,ncolprint)
c              enddo
c          endif
!
!     ---------------------------------------------------!



!
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERATURES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done if top_height = 1 or 3
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     trans_layers_above is the total transmissivity in the layers
!             above the current layer
!     fluxtop_clrsky(j) and trans_layers_above_clrsky(j) are the clear
!             sky versions of these quantities.

      if (top_height .eq. 1 .or. top_height .eq. 3) then


        !----------------------------------------------------------------------
        !
        !             DO CLEAR SKY RADIANCE CALCULATION FIRST
        !
        !compute water vapor continuum emissivity
        !this treatment follows Schwarkzopf and Ramasamy
        !JGR 1999,vol 104, pages 9467-9499.
        !the emissivity is calculated at a wavenumber of 955 cm-1,
        !or 10.47 microns
c        wtmair = 28.9644
c        wtmh20 = 18.01534
c        Navo = 6.023E+23
c        grav = 9.806650E+02
c        pstd = 1.013250E+06
c        t0 = 296.
c        if (ncolprint .ne. 0)
c     &         write(6,*)  'ilev   pw (kg/m2)   tauwv(j)      dem_wv'
        do 125 ilev=1,nlev
          do j=1,npoints
               !press and dpress are dyne/cm2 = Pascals *10
               press(j) = pfull(j,ilev)*10.
               dpress(j) = (phalf(j,ilev+1)-phalf(j,ilev))*10
               !atmden = g/cm2 = kg/m2 / 10
! minor GISS changes correct for unit difference in bygrav, bymrat,t0bypstd
               atmden(j) = dpress(j)*bygrav*0.01d0
               rvh20(j) = qv(j,ilev)*bymrat    !wtmair/wtmh20
               wk(j) = rvh20(j)*Navo*atmden(j)/wtmair
c               rhoave(j) = (press(j)/pstd)*(t0/at(j,ilev))
               rhoave(j) = (press(j)/at(j,ilev))*t0bypstd
               rh20s(j) = rvh20(j)*rhoave(j)
               rfrgn(j) = rhoave(j)-rh20s(j)
               tmpexp(j) = exp(-0.02d0*(at(j,ilev)-t0))
               tauwv(j) = wk(j)*1.d-20*(
     &           (0.0224697d0*rh20s(j)*tmpexp(j)) +
     &                (3.41817d-7*rfrgn(j)) )*0.98d0
               dem_wv(j,ilev) = 1. - exp( -1. * tauwv(j))
          enddo
c               if (ncolprint .ne. 0) then
c               do j=1,npoints ,1000
c               write(6,'(a10)') 'j='
c               write(6,'(8I10)') j
c               write(6,'(i2,1X,3(f8.3,3X))') ilev,
c     &           qv(j,ilev)*(phalf(j,ilev+1)-phalf(j,ilev))*bygrav,
c     &           tauwv(j),dem_wv(j,ilev)
c               enddo
c             endif
125     continue

        !initialize variables
        do j=1,npoints
          fluxtop_clrsky(j) = 0.
          trans_layers_above_clrsky(j)=1.
        enddo

        do ilev=1,nlev
          do j=1,npoints

            ! Black body emission at temperature of the layer
            ! 1307.27 = Planck constant c1 by wavelength 10.5
              bb(j)=1 / ( exp(pc1bylam/at(j,ilev)) - 1. )

              ! increase TOA flux by flux emitted from layer
              ! times total transmittance in layers above

                fluxtop_clrsky(j) = fluxtop_clrsky(j)
     &            + dem_wv(j,ilev)*bb(j)*trans_layers_above_clrsky(j)

                ! update trans_layers_above with transmissivity
              ! from this layer for next time around loop

                trans_layers_above_clrsky(j)=
     &            trans_layers_above_clrsky(j)*(1.-dem_wv(j,ilev))


          enddo
c            if (ncolprint.ne.0) then
c             do j=1,npoints ,1000
c              write(6,'(a10)') 'j='
c              write(6,'(8I10)') j
c              write (6,'(a)') 'ilev:'
c              write (6,'(I2)') ilev
c
c              write (6,'(a)')
c     &        'emiss_layer,100.*bb(j),100.*f,total_trans:'
c              write (6,'(4(f7.2,1X))') dem_wv(j,ilev),100.*bb(j),
c     &             100.*fluxtop_clrsky(j),trans_layers_above_clrsky(j)
c             enddo
c            endif

        enddo   !loop over level

        do j=1,npoints
          !add in surface emission
          bb(j)=1/( exp(pc1bylam/skt(j)) - 1. )
          !bb(j)=5.67d-8*skt(j)**4

          fluxtop_clrsky(j) = fluxtop_clrsky(j) + emsfc_lw * bb(j)
     &     * trans_layers_above_clrsky(j)
        enddo

c        if (ncolprint.ne.0) then
c        do j=1,npoints ,1000
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write (6,'(a)') 'id:'
c          write (6,'(a)') 'surface'
c
c          write (6,'(a)') 'emsfc,100.*bb(j),100.*f,total_trans:'
c          write (6,'(4(f7.2,1X))') emsfc_lw,100.*bb(j),
c     &      100.*fluxtop_clrsky(j),
c     &       trans_layers_above_clrsky(j)
c        enddo
c      endif


        !
        !           END OF CLEAR SKY CALCULATION
        !
        !----------------------------------------------------------------



c        if (ncolprint.ne.0) then
c
c        do j=1,npoints ,1000
c            write(6,'(a10)') 'j='
c            write(6,'(8I10)') j
c            write (6,'(a)') 'ts:'
c            write (6,'(8f7.2)') (skt(j),ibox=1,ncolprint)
c
c            write (6,'(a)') 'ta_rev:'
c            write (6,'(8f7.2)')
c     &       ((at(j,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
c
c        enddo
c        endif
        !loop over columns
        do ibox=1,ncol
          do j=1,npoints
            fluxtop(j,ibox)=0.
            trans_layers_above(j,ibox)=1.
          enddo
        enddo

        do ilev=1,nlev
              do j=1,npoints
                ! Black body emission at temperature of the layer

              bb(j)=1 / ( exp(pc1bylam/at(j,ilev)) - 1. )
              !bb(j)= 5.67d-8*at(j,ilev)**4
              enddo

            do ibox=1,ncol
              do j=1,npoints

              ! emissivity for point in this layer
                if (frac_out(j,ibox,ilev).eq.1) then
                dem(j,ibox)= 1. -
     &             ( (1. - dem_wv(j,ilev)) * (1. -  dem_s(j,ilev)) )
                else if (frac_out(j,ibox,ilev).eq.2) then
                dem(j,ibox)= 1. -
     &             ( (1. - dem_wv(j,ilev)) * (1. -  dem_c(j,ilev)) )
                else
                dem(j,ibox)=  dem_wv(j,ilev)
                end if


                ! increase TOA flux by flux emitted from layer
              ! times total transmittance in layers above

                fluxtop(j,ibox) = fluxtop(j,ibox)
     &            + dem(j,ibox) * bb(j)
     &            * trans_layers_above(j,ibox)

                ! update trans_layers_above with transmissivity
              ! from this layer for next time around loop

                trans_layers_above(j,ibox)=
     &            trans_layers_above(j,ibox)*(1.-dem(j,ibox))

              enddo ! j
            enddo ! ibox

c            if (ncolprint.ne.0) then
c              do j=1,npoints,1000
c              write (6,'(a)') 'ilev:'
c              write (6,'(I2)') ilev
c
c              write(6,'(a10)') 'j='
c              write(6,'(8I10)') j
c              write (6,'(a)') 'emiss_layer:'
c              write (6,'(8f7.2)') (dem(j,ibox),ibox=1,ncolprint)
c
c              write (6,'(a)') '100.*bb(j):'
c              write (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
c
c              write (6,'(a)') '100.*f:'
c              write (6,'(8f7.2)')
c     &         (100.*fluxtop(j,ibox),ibox=1,ncolprint)
c
c              write (6,'(a)') 'total_trans:'
c              write (6,'(8f7.2)')
c     &          (trans_layers_above(j,ibox),ibox=1,ncolprint)
c            enddo
c          endif

        enddo ! ilev


          do j=1,npoints
            !add in surface emission
            bb(j)=1/( exp(pc1bylam/skt(j)) - 1. )
            !bb(j)=5.67d-8*skt(j)**4
          end do

        do ibox=1,ncol
          do j=1,npoints

            !add in surface emission

            fluxtop(j,ibox) = fluxtop(j,ibox)
     &         + emsfc_lw * bb(j)
     &         * trans_layers_above(j,ibox)

          end do
        end do

c        if (ncolprint.ne.0) then
c
c          do j=1,npoints ,1000
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write (6,'(a)') 'id:'
c          write (6,'(a)') 'surface'
c
c          write (6,'(a)') 'emiss_layer:'
c          write (6,'(8f7.2)') (dem(1,ibox),ibox=1,ncolprint)
c
c          write (6,'(a)') '100.*bb(j):'
c          write (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
c
c          write (6,'(a)') '100.*f:'
c          write (6,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)
c          end do
c      endif

        !now that you have the top of atmosphere radiance account
        !for ISCCP procedures to determine cloud top temperature

        !account for partially transmitting cloud recompute flux
        !ISCCP would see assuming a single layer cloud
        !note choice here of 2.13, as it is primarily ice
        !clouds which have partial emissivity and need the
        !adjustment performed in this section
        !
      !If it turns out that the cloud brightness temperature
      !is greater than 260K, then the liquid cloud conversion
        !factor of 2.56 is used.
      !
        !Note that this is discussed on pages 85-87 of
        !the ISCCP D level documentation (Rossow et al. 1996)

          do j=1,npoints
            !compute minimum brightness temperature and optical depth
            btcmin(j) = 1. /  ( exp(pc1bylam/(attrop(j)-5.)) - 1. )
          enddo
        do ibox=1,ncol
          do j=1,npoints
            transmax(j) = (fluxtop(j,ibox)-btcmin(j))
     &                /(fluxtop_clrsky(j)-btcmin(j))
          !note that the initial setting of tauir(j) is needed so that
          !tauir(j) has a realistic value should the next if block be
          !bypassed
            tauir(j) = tau(j,ibox) * byic
            taumin(j) = -1. * log(max(min(transmax(j),0.9999999d0),0
     *           .001d0))

          enddo

          if (top_height .eq. 1) then
            do j=1,npoints
              if (transmax(j) .gt. 0.001 .and.
     &          transmax(j) .le. 0.9999999) then
                fluxtopinit(j) = fluxtop(j,ibox)
              tauir(j) = tau(j,ibox) *byic
              endif
            enddo
            do icycle=1,2
              do j=1,npoints
                if (tau(j,ibox) .gt. (tauchk            )) then
                if (transmax(j) .gt. 0.001 .and.
     &            transmax(j) .le. 0.9999999) then
                  emcld(j,ibox) = 1. - exp(-1. * tauir(j)  )
                  fluxtop(j,ibox) = fluxtopinit(j) -
     &              ((1.-emcld(j,ibox))*fluxtop_clrsky(j))
                  fluxtop(j,ibox)=max(1.d-06,
     &              (fluxtop(j,ibox)/emcld(j,ibox)))
                  tb(j,ibox)= pc1bylam
     &              / (log(1. + (1./fluxtop(j,ibox))))
                  if (tb(j,ibox) .gt. 260.) then
                  tauir(j) = tau(j,ibox) * bywc
                  end if
                end if
                end if
              enddo
            enddo

          endif

          do j=1,npoints
            if (tau(j,ibox) .gt. (tauchk            )) then
                !cloudy box
                tb(j,ibox)= pc1bylam/ (log(1. + (1./fluxtop(j,ibox))))
                if (top_height.eq.1.and.tauir(j).lt.taumin(j)) then
                         tb(j,ibox) = attrop(j) - 5.
                   tau(j,ibox) = 2.13d0*taumin(j)
                end if
            else
                !clear sky brightness temperature
                tb(j,ibox) = pc1bylam/(log(1.+(1./fluxtop_clrsky(j))))
            end if
          enddo ! j
        enddo ! ibox

c        if (ncolprint.ne.0) then
c
c          do j=1,npoints,1000
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c
c          write (6,'(a)') 'attrop:'
c          write (6,'(8f7.2)') (attrop(j))
c
c          write (6,'(a)') 'btcmin:'
c          write (6,'(8f7.2)') (btcmin(j))
c
c          write (6,'(a)') 'fluxtop_clrsky*100:'
c          write (6,'(8f7.2)')
c     &      (100.*fluxtop_clrsky(j))
c
c          write (6,'(a)') '100.*f_adj:'
c          write (6,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)
c
c          write (6,'(a)') 'transmax:'
c          write (6,'(8f7.2)') (transmax(ibox),ibox=1,ncolprint)
c
c          write (6,'(a)') 'tau:'
c          write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
c
c          write (6,'(a)') 'emcld:'
c          write (6,'(8f7.2)') (emcld(j,ibox),ibox=1,ncolprint)
c
c          write (6,'(a)') 'total_trans:'
c          write (6,'(8f7.2)')
c     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)
c
c          write (6,'(a)') 'total_emiss:'
c          write (6,'(8f7.2)')
c     &        (1.0-trans_layers_above(j,ibox),ibox=1,ncolprint)
c
c          write (6,'(a)') 'total_trans:'
c          write (6,'(8f7.2)')
c     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)
c
c          write (6,'(a)') 'ppout:'
c          write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
c          enddo ! j
c      endif

      end if

!     ---------------------------------------------------!

!
!     ---------------------------------------------------!
!     DETERMINE CLOUD TOP PRESSURE
!
!     again the 2 methods differ according to whether
!     or not you use the physical cloud top pressure (top_height = 2)
!     or the radiatively determined cloud top pressure (top_height = 1 or 3)
!

      !compute cloud top pressure
      do 30 ibox=1,ncol
        !segregate according to optical thickness
        if (top_height .eq. 1 .or. top_height .eq. 3) then
          !find level whose temperature
          !most closely matches brightness temperature
          do j=1,npoints
            nmatch(j)=0
          enddo
          do 29 ilev=1,nlev-1
            !cdir nodep
            do j=1,npoints
              if ((at(j,ilev)   .ge. tb(j,ibox) .and.
     &          at(j,ilev+1) .lt. tb(j,ibox)) .or.
     &          (at(j,ilev) .le. tb(j,ibox) .and.
     &          at(j,ilev+1) .gt. tb(j,ibox))) then

                nmatch(j)=nmatch(j)+1
                if(abs(at(j,ilev)-tb(j,ibox)) .lt.
     &            abs(at(j,ilev+1)-tb(j,ibox))) then
                  match(j,nmatch(j))=ilev
                else
                  match(j,nmatch(j))=ilev+1
                end if
              end if
            enddo
29        continue

          do j=1,npoints
            if (nmatch(j) .ge. 1) then
              ptop(j,ibox)=pfull(j,match(j,nmatch(j)))
              levmatch(j,ibox)=match(j,nmatch(j))
            else
              if (tb(j,ibox) .lt. atmin(j)) then
                ptop(j,ibox)=ptrop(j)
                levmatch(j,ibox)=itrop(j)
              end if
              if (tb(j,ibox) .gt. atmax(j)) then
                ptop(j,ibox)=pfull(j,nlev)
                levmatch(j,ibox)=nlev
              end if
            end if
          enddo ! j

        else ! if (top_height .eq. 1 .or. top_height .eq. 3)

          do j=1,npoints
            ptop(j,ibox)=0.
          enddo
          do ilev=1,nlev
            do j=1,npoints
              if ((ptop(j,ibox) .eq. 0. )
     &           .and.(frac_out(j,ibox,ilev) .ne. 0)) then
                ptop(j,ibox)=pfull(j,ilev)
              levmatch(j,ibox)=ilev
              end if
            end do
          end do
        end if

        do j=1,npoints
          if (tau(j,ibox) .le. (tauchk            )) then
            ptop(j,ibox)=0.
            levmatch(j,ibox)=0
          endif
        enddo

30    continue

!
!
!     ---------------------------------------------------!


!
!     ---------------------------------------------------!
!     DETERMINE ISCCP CLOUD TYPE FREQUENCIES
!
!     Now that ptop and tau have been determined,
!     determine amount of each of the 49 ISCCP cloud
!     types
!
!     Also compute grid box mean cloud top pressure and
!     optical thickness.  The mean cloud top pressure and
!     optical thickness are averages over the cloudy
!     area only. The mean cloud top pressure is a linear
!     average of the cloud top pressures.  The mean cloud
!     optical thickness is computed by converting optical
!     thickness to an albedo, averaging in albedo units,
!     then converting the average albedo back to a mean
!     optical thickness.
!

      !compute isccp frequencies

      !reset frequencies
      do 38 ilev=1,7
      do 38 ilev2=1,7
        do j=1,npoints !
             fq_isccp(j,ilev,ilev2)=0.
        enddo
38    continue

      !reset variables need for averaging cloud properties
      do j=1,npoints
        totalcldarea(j) = 0.
        meanalbedocld(j) = 0.
        meanptop(j) = 0.
        meantaucld(j) = 0.
      enddo ! j

c      boxarea = 1./real(ncol,kind=8)

      do 39 ibox=1,ncol
        do j=1,npoints

c          if (tau(j,ibox) .gt. (tauchk            )
c     &      .and. ptop(j,ibox) .gt. 0.) then
c              box_cloudy(j,ibox)=.true.
c          endif

          box_cloudy(j,ibox)= (tau(j,ibox) .gt. tauchk .and. ptop(j,ibox
     *         ) .gt. 0.)

          if (box_cloudy(j,ibox)) then

              ! totalcldarea always diagnosed day or night
              totalcldarea(j) = totalcldarea(j) + boxarea

              if (sunlit(j).eq.1) then

                ! tau diagnostics only with sunlight

                boxtau(j,ibox) = tau(j,ibox)

                !convert optical thickness to albedo
                albedocld(j,ibox)
     &            =real(invtau(min(nint(100.*tau(j,ibox)),45000)))

                !contribute to averaging
              meanalbedocld(j) = meanalbedocld(j)
     &            +albedocld(j,ibox)*boxarea

            endif

          endif

          if (sunlit(j).eq.1 .or. top_height .eq. 3) then

            !convert ptop to millibars
            ptop(j,ibox)=ptop(j,ibox) / 100.

            !save for output cloud top pressure and optical thickness
            boxptop(j,ibox) = ptop(j,ibox)

            if (box_cloudy(j,ibox)) then

              meanptop(j) = meanptop(j) + ptop(j,ibox)*boxarea

              !reset itau(j), ipres(j)
              itau(j) = 0
              ipres(j) = 0

              !determine optical depth category
              if (tau(j,ibox) .lt. isccp_taumin) then
                  itau(j)=1
              else if (tau(j,ibox) .ge. isccp_taumin
     &          .and. tau(j,ibox) .lt. 1.3) then
                itau(j)=2
              else if (tau(j,ibox) .ge. 1.3
     &          .and. tau(j,ibox) .lt. 3.6) then
                itau(j)=3
              else if (tau(j,ibox) .ge. 3.6
     &          .and. tau(j,ibox) .lt. 9.4) then
                  itau(j)=4
              else if (tau(j,ibox) .ge. 9.4
     &          .and. tau(j,ibox) .lt. 23.) then
                  itau(j)=5
              else if (tau(j,ibox) .ge. 23.
     &          .and. tau(j,ibox) .lt. 60.) then
                  itau(j)=6
              else if (tau(j,ibox) .ge. 60.) then
                  itau(j)=7
              end if

              !determine cloud top pressure category
              if (    ptop(j,ibox) .gt. 0.
     &          .and.ptop(j,ibox) .lt. 180.) then
                  ipres(j)=1
              else if(ptop(j,ibox) .ge. 180.
     &          .and.ptop(j,ibox) .lt. 310.) then
                  ipres(j)=2
              else if(ptop(j,ibox) .ge. 310.
     &          .and.ptop(j,ibox) .lt. 440.) then
                  ipres(j)=3
              else if(ptop(j,ibox) .ge. 440.
     &          .and.ptop(j,ibox) .lt. 560.) then
                  ipres(j)=4
              else if(ptop(j,ibox) .ge. 560.
     &          .and.ptop(j,ibox) .lt. 680.) then
                  ipres(j)=5
              else if(ptop(j,ibox) .ge. 680.
     &          .and.ptop(j,ibox) .lt. 800.) then
                  ipres(j)=6
              else if(ptop(j,ibox) .ge. 800.) then
                  ipres(j)=7
              end if

              !update frequencies
              if(ipres(j) .gt. 0.and.itau(j) .gt. 0)
     *             fq_isccp(j,itau(j),ipres(j))=
     &             fq_isccp(j,itau(j),ipres(j))+ boxarea
            end if

          end if

        enddo ! j
39    continue

      !compute mean cloud properties
      do j=1,npoints
        nbox(j)=totalcldarea(j)*ncol    ! for backwards compatibility
        if (totalcldarea(j) .gt. 0.) then
         meanptop(j) = meanptop(j) / totalcldarea(j)
          if (sunlit(j).eq.1) then
            meanalbedocld(j) = meanalbedocld(j) / totalcldarea(j)
          meantaucld(j) = tautab(min(255,max(1,nint(meanalbedocld(j)))))
          end if
        end if
      enddo ! j

!
!     ---------------------------------------------------!

!     ---------------------------------------------------!
!     OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
!
c      if (debugcol.ne.0) then
c
c         do j=1,npoints,debugcol
c
c            !produce character output
c            do ilev=1,nlev
c              do ibox=1,ncol
c                   acc(ilev,ibox)=0
c              enddo
c            enddo
c
c            do ilev=1,nlev
c              do ibox=1,ncol
c                   acc(ilev,ibox)=frac_out(j,ibox,ilev)*2
c                   if (levmatch(j,ibox) .eq. ilev)
c     &                 acc(ilev,ibox)=acc(ilev,ibox)+1
c              enddo
c            enddo
c
c             !print test
c
c          write(ftn09,11) j
c11        format('ftn09.',i4.4)
c          open(9, FILE=ftn09, FORM='FORMATTED')
c
c             write(9,'(a1)') ' '
c             write(9,'(10i5)')
c     &                  (ilev,ilev=5,nlev,5)
c             write(9,'(a1)') ' '
c
c             do ibox=1,ncol
c               write(9,'(40(a1),1x,40(a1))')
c     &           (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev)
c     &           ,(cchar(acc(ilev,ibox)+1),ilev=1,nlev)
c             end do
c             close(9)
c
c             if (ncolprint.ne.0) then
c               write(6,'(a1)') ' '
c                    write(6,'(a2,1X,5(a7,1X),a50)')
c     &                  'ilev',
c     &                  'pfull','at',
c     &                  'cc*100','dem_s','dtau_s',
c     &                  'cchar'

!               do 4012 ilev=1,nlev
!                    write(6,'(60i2)') (box(i,ilev),i=1,ncolprint)
!                   write(6,'(i2,1X,5(f7.2,1X),50(a1))')
!     &                  ilev,
!     &                  pfull(j,ilev)/100.,at(j,ilev),
!     &                  cc(j,ilev)*100.0,dem_s(j,ilev),dtau_s(j,ilev)
!     &                  ,(cchar(acc(ilev,ibox)+1),ibox=1,ncolprint)
!4012           continue
c               write (6,'(a)') 'skt(j):'
c               write (6,'(8f7.2)') skt(j)
c
c               write (6,'(8I7)') (ibox,ibox=1,ncolprint)
c
c               write (6,'(a)') 'tau:'
c               write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
c
c               write (6,'(a)') 'tb:'
c               write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
c
c               write (6,'(a)') 'ptop:'
c               write (6,'(8f7.2)') (ptop(j,ibox),ibox=1,ncolprint)
c             endif
c
c        enddo
c
c      end if

      return
      end

      subroutine get_dq_cond(sm,qm,plk,mass,lhx,pl,   ! input
     *                       dqsum,fcond)         ! output
!@sum get_dq calculation of condensation from vapour
      USE CONSTANT, only : bysha
      IMPLICIT NONE
!@var SM,QM heat and water content
      REAL*8, INTENT(IN) :: sm,qm
!@var MASS air mass
!@var PLK P**K for mid point
!@var LHX relevant latent heat
!@var PL mid point pressure
      REAL*8, INTENT(IN) :: plk,mass,lhx,pl
!@var DQSUM amount of water change (+ve is condensation)
!@var FCOND fractional amount of vapour that condenses
      REAL*8, INTENT(OUT) :: dqsum,fcond
      REAL*8 TP,QST,QSAT,DQSATDT,DQ,QMT,SLH
      INTEGER N

      DQSUM=0.
      FCOND=0.
      IF (QM.gt.0) THEN
        SLH=LHX*BYSHA
        QMT=QM
        TP=SM*PLK/MASS
        DO N=1,3
          QST=QSAT(TP,LHX,PL)
          DQ=(QMT-MASS*QST)/(1.+SLH*QST*DQSATDT(TP,LHX))
          TP=TP+SLH*DQ/MASS
          QMT=QMT-DQ
          DQSUM=DQSUM+DQ
        END DO
C**** condensing (DQSUM > 0)
        DQSUM=MAX(0d0,MIN(DQSUM,QM))
        FCOND=DQSUM/QM
      END IF

      return
      end subroutine get_dq_cond

      subroutine get_dq_evap(sm,qm,plk,mass,lhx,pl,cond,   ! input
     *                       dqsum,fevp) ! output
!@sum get_dq calculation of evaporation from condensate
      USE CONSTANT, only : bysha
      IMPLICIT NONE
!@var SM,QM heat and water content
      REAL*8, INTENT(IN) :: sm,qm
!@var MASS air mass
!@var PLK P**K for mid point
!@var LHX relevant latent heat
!@var PL mid point pressure
!@var COND amount of condensate
      REAL*8, INTENT(IN) :: plk,mass,lhx,pl,cond
!@var DQSUM amount of water change (+ve is evaporation)
!@var FEVP fractional amount of condensate that evaporates
      REAL*8, INTENT(OUT) :: dqsum,fevp
      REAL*8 TP,QST,QSAT,DQSATDT,DQ,QMT,SLH
      INTEGER N

      DQSUM=0.
      FEVP=0.
      IF (COND.gt.0) THEN
        SLH=LHX*BYSHA
        QMT=QM
        TP=SM*PLK/MASS
        DO N=1,3
          QST=QSAT(TP,LHX,PL)
          DQ=(QMT-MASS*QST)/(1.+SLH*QST*DQSATDT(TP,LHX))
          TP=TP+SLH*DQ/MASS
          QMT=QMT-DQ
          DQSUM=DQSUM-DQ
        END DO
C**** evaporating (DQ < 0, DQSUM > 0)
        DQSUM=MAX(0d0,MIN(DQSUM,COND))
        FEVP=DQSUM/COND
      END IF

      return
      end subroutine get_dq_evap
