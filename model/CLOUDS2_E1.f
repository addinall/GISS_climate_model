#include "rundeck_opts.h"

      MODULE CLOUDS
!@sum  CLOUDS column physics of moist conv. and large-scale condensation
!@auth M.S.Yao/A. Del Genio (modifications by Gavin Schmidt)
!@cont MSTCNV,LSCOND
      USE CONSTANT, only : rgas,grav,lhe,lhs,lhm,sha,bysha,pi,by6
     *     ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd,rhow
     *  ,twopi
#ifdef CLD_AER_CDNC
     *,kapa,mair,gasc
#endif
      USE MODEL_COM, only : lm,dtsrc,itime,coupled_chem
#ifdef CLD_AER_CDNC
     * ,ptop,psf,ls1,sig,sige
#endif
#ifdef BLK_2MOM
      USE mo_bulk2m_driver_gcm, ONLY: execute_bulk2m_driver
#endif
      USE QUSDEF, only : nmom,xymoms,zmoms,zdir
#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm,trname,t_qlimit,ntm_soa
#ifdef TRACERS_WATER
     &     ,nGAS, nPART, nWATER, tr_wd_TYPE,
     *     tr_evap_fact
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,Ntm_dust
#endif
#endif
#ifdef RUNTIME_NTM
      USE TRACER_COM, only:
     &   ntix=>ntix_cld,tm=>tm_col,tmom=>tmom_col
     &   ,trdnl=>trdn_col,trwml=>trwm_col,trsvwml=>trsvwm_col
     &   ,trprss=>trprss_cld,trprmc=>trprmc_cld
#endif
#endif
#ifdef BLK_2MOM
#ifdef TRACERS_AMP
      USE AMP_AEROSOL, only: NACTC,NAERC
      USE AERO_CONFIG, only: NMODES
#endif
#endif

CCC   USE RANDOM
      IMPLICIT NONE
      SAVE
C**** parameters and constants
      REAL*8, PARAMETER :: TI=233.16d0   !@param TI pure ice limit
      REAL*8, PARAMETER :: CLDMIN=.25d0 !@param CLDMIN min MC/LSC region
!@param WMU critical cloud water content for rapid conversion (g m**-3)
      REAL*8, PARAMETER :: WMU=.25
      REAL*8, PARAMETER :: WMUL=.5       !@param WMUL WMU over land
      REAL*8, PARAMETER :: WMUI=.1d0     !@param WMUI WMU for ice clouds
      REAL*8, PARAMETER :: BRCLD=.2d0    !@param BRCLD for cal. BYBR
      REAL*8, PARAMETER :: SLHE=LHE*BYSHA
      REAL*8, PARAMETER :: SLHS=LHS*BYSHA
!@param CCMUL multiplier for convective cloud cover
!@param CCMUL1 multiplier for deep anvil cloud cover
!@param CCMUL2 multiplier for shallow anvil cloud cover
!@param COETAU multiplier for convective cloud optical thickness
      REAL*8, PARAMETER :: CCMUL=1.,CCMUL1=5.,CCMUL2=3.,COETAU=.08d0

      REAL*8 :: BYBR,BYDTsrc,XMASS,PLAND
      REAL*8 :: RTEMP,CMX,RCLDX,CONTCE1,CONTCE2
!@var BYBR factor for converting cloud particle radius to effect. radius
!@var XMASS dummy variable
!@var PLAND land fraction

C**** Set-able variables
!@dbparam LMCM max level for originating MC plumes
      INTEGER :: LMCM = -1 ! defaults to LS1-1 if not set in rundeck
!@dbparam ISC integer to turn on computation of stratocumulus clouds
      INTEGER :: ISC = 0  ! set ISC=1 to compute stratocumulus clouds
!@dbparam U00wtrX multiplies U00ice for critical humidity for water clds
      REAL*8 :: U00wtrX = 1.0d0     ! default
!@dbparam U00ice critical humidity for ice cloud condensation
      REAL*8 :: U00ice = .7d0       ! default
!@dbparam U00a tuning knob for U00 above 850 mb without moist convection
!@dbparam U00b tuning knob for U00 below 850 mb and in convective regions
      REAL*8 :: U00a = 0.55d0       ! default
      REAL*8 :: U00b = 1.00d0       ! default
!@dbparam funio_denominator funio denominator
      REAL*8 :: funio_denominator=15.d0  ! default
!@dbparam autoconv_multiplier autoconversion rate multiplier
      REAL*8 :: autoconv_multiplier=1.d0 ! default
!@dbparam radius_multiplier cloud particle radius multiplier
      REAL*8 :: radius_multiplier=1.d0   ! default
!@dbparam wmui_multiplier critical ice cloud water multiplier
      REAL*8 :: wmui_multiplier=1.d0     ! default
!@dbparam entrainment_cont1 constant for entrainment rate, plume 1
      REAL*8 :: entrainment_cont1=.2d0   ! default
!@dbparam entrainment_cont2 constant for entrainment rate, plume 2
      REAL*8 :: entrainment_cont2=.2d0   ! default
!@dbparam HRMAX maximum distance an air parcel rises from surface
      REAL*8 :: HRMAX = 1000.d0     ! default (m)
!@dbparam RWCldOX multiplies part.size of water clouds over ocean
      REAL*8 :: RWCldOX=1.d0
!@dbparam RICldX multiplies part.size of ice clouds at 1000mb
!@+       RICldX changes linearly to 1 as p->0mb
      REAL*8 :: RICldX=1.d0 , xRICld
!@dbparam do_blU00 =1 if boundary layer U00 is treated differently
      INTEGER :: do_blU00=0     ! default is to disable this

#ifdef TRACERS_ON
!@var ntx,NTIX: Number and Indices of active tracers used in convection
#ifndef RUNTIME_NTM
      integer, dimension(ntm) :: ntix
#endif
      integer ntx
#endif
C**** ISCCP diag related variables
      INTEGER,PARAMETER :: ncol =20    !@var ncol number of subcolumns
!@var tautab look-up table to convert count value to optical thickness
!@var invtau look-up table to convert optical thickness to count value
      real*8 :: tautab(0:255)
      integer :: invtau(-20:45000)

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
     *     ,FSSL,FSUB,FCONV,FMCL,VLAT,DDMFLX,WTURB,TVL,W2L,GZL
     *     ,SAVWL,SAVWL1,SAVE1L,SAVE2L,DPHASHLW,DPHADEEP,DGSHLW,DGDEEP
     *     ,QDNL,TDNL
#ifdef BLK_2MOM
     *     ,WMICE
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
#ifdef BLK_2MOM
!@var WMICE ice water mixing ratio (kg/kg)
#endif
!@var VSUBL downward vertical velocity due to cumulus subsidence (cm/s)
!@var MCFLX, DGDSM, DPHASE, DQCOND, DGDQM dummy variables
!@var DDMFLX accumulated downdraft mass flux (mb)
!@var AQ time change rate of specific humidity (s**-1)
!@var DPDT time change rate of pressure (mb/s)
!@var FSSL grid fraction for large-scale clouds
!@var FCONV convective fraction
!@var FSUB subsiding fraction
!@var FMCL grid fraction for moist convection
!@var VLAT dummy variable
      REAL*8, DIMENSION(LM+1) :: PRECNVL
!@var PRECNVL convective precip entering the layer top
C**** new arrays must be set to model arrays in driver (before MSTCNV)
      REAL*8, DIMENSION(LM) :: SDL,WML
!@var SDL vertical velocity in sigma coordinate
!@var WML cloud water mixing ratio (kg/kg)
C**** new arrays must be set to model arrays in driver (after MSTCNV)
      REAL*8, DIMENSION(LM) :: TAUMCL,SVLATL,CLDMCL,SVLHXL,SVWMXL
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
!@var ALWWM and ALWIM  etc are liquid water contents
!@var CTTEM,CD3DL,CL3DL,CI3DL,SMLWP are cld temp, cld thickness,cld water,LWP
      REAL*8, DIMENSION(LM) ::CTEML,CD3DL,CL3DL,CI3DL,CDN3DL,CRE3DL
      REAL*8 SMLWP
!@var SME is the TKE in 1 D from e(l) = egcm(l,i,j)  (m^2/s^2)
      REAL*8, DIMENSION(LM)::SME
      INTEGER NLSW,NLSI,NMCW,NMCI
#endif
C**** new arrays must be set to model arrays in driver (before LSCOND)
      REAL*8, DIMENSION(LM) :: TTOLDL,CLDSAVL,CLDSV1
!@var TTOLDL previous potential temperature
!@var CLDSAVL saved large-scale cloud cover
#ifdef CLD_AER_CDNC
      REAL*8, DIMENSION(LM)::OLDCDL,OLDCDI
!@var OLDCDL is saved CDNC
!@var OLDCDI is saved ice crystal number
#endif
C**** new arrays must be set to model arrays in driver (after LSCOND)
      REAL*8, DIMENSION(LM) :: SSHR,DCTEI,TAUSSL,CLDSSL
!@var SSHR,DCTEI height diagnostics of dry and latent heating by MC
!@var TAUSSL large-scale cloud optical thickness
!@var CLDSSL large-scale cloud cover

!@var SM,QM Vertical profiles of T/Q
      REAL*8, DIMENSION(LM) :: SM,QM
      REAL*8, DIMENSION(NMOM,LM) :: SMOM,QMOM,SMOMMC,QMOMMC,
     *  SMOMLS,QMOMLS

#ifdef TRACERS_ON
#ifndef RUNTIME_NTM
!@var TM Vertical profiles of tracers
      REAL*8, DIMENSION(LM,NTM) :: TM
      REAL*8, DIMENSION(nmom,lm,ntm) :: TMOM
!@var TRDNL tracer concentration in lowest downdraft (kg/kg)
      REAL*8, DIMENSION(NTM,LM) :: TRDNL
      COMMON/CLD_TRCCOM/TM,TMOM,TRDNL
!$OMP  THREADPRIVATE (/CLD_TRCCOM/)
#endif /*RUNTIME_NTM*/
#ifdef TRACERS_WATER
!@var TRWML Vertical profile of liquid water tracers (kg)
!@var TRSVWML New liquid water tracers from m.c. (kg)
#ifndef RUNTIME_NTM
      REAL*8, DIMENSION(NTM,LM) :: TRWML, TRSVWML
!@var TRPRSS super-saturated tracer precip (kg)
!@var TRPRMC moist convective tracer precip (kg)
      REAL*8, DIMENSION(NTM)    :: TRPRSS,TRPRMC
#endif /*RUNTIME_NTM*/
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
#ifndef RUNTIME_NTM
      COMMON/CLD_WTRTRCCOM/TRWML, TRSVWML,TRPRSS,TRPRMC
#endif
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
!$OMP  THREADPRIVATE (/CLD_WTRTRCCOM/)
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
!$OMP  THREADPRIVATE (/CLD_PRECDUST/)
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
CCOMP  does not work yet:
CCOMP  THREADPRIVATE (RA,UM,VM,U_0,V_0,PLE,PL,PLK,AIRM,BYAM,ETAL
CCOMP*  ,TL,QL,TH,RH,WMX,VSUBL,MCFLX,SSHR,DGDSM,DPHASE
CCOMP*  ,DTOTW,DQCOND,DCTEI,DGDQM,dxypij,DDMFLX
CCOMP*  ,AQ,DPDT,PRECNVL,SDL,WML,SVLATL,SVLHXL,SVWMXL,CSIZEL,RH1
CCOMP*  ,TTOLDL,CLDSAVL,TAUMCL,CLDMCL,TAUSSL,CLDSSL,RNDSSL
CCOMP*  ,SM,QM,SMOM,QMOM,PEARTH,TS,QS,US,VS,DCL,RIS,RI1,RI2, AIRXL
CCOMP*  ,PRCPMC,PRCPSS,HCNDSS,WMSUM,CLDSLWIJ,CLDDEPIJ,LMCMAX
CCOMP*  ,LMCMIN,KMAX,DEBUG)
CCOMP* RA,UM,VM,UM1,VM1,U_0,V_0 are no longer part of this COMMON block!
      COMMON/CLDPRV/PLE,PL,PLK,AIRM,BYAM,ETAL
     *  ,TL,QL,TH,RH,WMX,VSUBL,MCFLX,SSHR,DGDSM,DPHASE,LHP
     *  ,DPHASHLW,DPHADEEP,DGSHLW,DGDEEP
     *  ,DTOTW,DQCOND,DCTEI,DGDQM,DXYPIJ,DDMFLX,PLAND
     *  ,AQ,DPDT,PRECNVL,SDL,WML,SVLATL,SVLHXL,SVWMXL,CSIZEL,RH1
     *  ,TTOLDL,CLDSAVL,TAUMCL,CLDMCL,TAUSSL,CLDSSL,RNDSSL
     *  ,SM,QM,SMOM,QMOM,PEARTH,TS,QS,US,VS,RIS,RI1,RI2, AIRXL
     *  ,SMOMMC,QMOMMC,SMOMLS,QMOMLS,CLDSV1,PRHEAT,QDNL,TDNL
     *  ,PRCPMC,PRCPSS,HCNDSS,WMSUM,CLDSLWIJ,CLDDEPIJ,VLAT
#ifdef CLD_AER_CDNC
     *  ,ACDNWM,ACDNIM,ACDNWS,ACDNIS
     *  ,AREWM,AREIM,AREWS,AREIS,ALWIM,ALWWM
     *  ,OLDCDL,OLDCDI
     *  ,SME
     *  ,CTEML,CD3DL,CL3DL,CI3DL,SMLWP,CDN3DL,CRE3DL
     *  ,WMCLWP,WMCTWP
#endif
     *  ,FSUB,FCONV,FSSL,FMCL
#ifdef CLD_AER_CDNC
     *  ,NLSW,NLSI,NMCW,NMCI
#endif
     *  ,prebar1,LMCMAX,LMCMIN,KMAX,DCL,DEBUG  ! int/logic last (alignment)
!$OMP  THREADPRIVATE (/CLDPRV/)

      CONTAINS

      SUBROUTINE MSTCNV(IERR,LERR,i_debug,j_debug)
!@sum  MSTCNV moist convective processes (precip, convective clouds,...)
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
!@calls adv1d,QSAT,DQSATDT,THBAR
      IMPLICIT NONE
      REAL*8 LHX,MPLUME,MPLUM1,MCLOUD,MPMAX,SENV,QENV
!@var LHX latent heat of evaporation (J/Kg)
!@var MPLUME,MPLUM1 mass of convective plume (mb)
!@var MCLOUD air mass available for re-evaporation of precip
!@var MPMAX convective plume at the detrainment level
!@var SENV,QENV dummy variables
!@var FMC1 grid fraction for moist convection
      REAL*8 FMC1

C**** functions
      REAL*8 :: QSAT, DQSATDT
!@var QSAT saturation humidity
!@var DQSATDT dQSAT/dT

      REAL*8, DIMENSION(0:LM) :: CM     !@var CM air mass of subsidence
      REAL*8, DIMENSION(KMAX) :: UMP,VMP,UMDN,VMDN
!@var UMP, VMP momentum carried by convective plumes
!@var UMDN,VMDN dummy variables
!@var DQM,DSM,DQMR,DSMR Vertical profiles of T/Q and changes
      REAL*8, DIMENSION(LM) ::
     * SMOLD,QMOLD, DQM,DSM,DQMR,DSMR
!@var SMOLD,QMOLD profiles prior to any moist convection
      REAL*8, DIMENSION(LM) :: F,CMNEG
      REAL*8, DIMENSION(NMOM,LM) :: FMOM
      REAL*8, DIMENSION(NMOM,LM) :: SMOMOLD,QMOMOLD
      REAL*8, DIMENSION(NMOM,LM) :: DSMOM,DQMOM,DSMOMR,DQMOMR
      REAL*8, DIMENSION(NMOM) ::
     &     SMOMP,QMOMP, SMOMPMAX,QMOMPMAX, SMOMDN,QMOMDN

#ifdef TRACERS_ON
!@var TMOLD: old TM (tracer mass)
!@var DTM,DTMR: Vertical profiles of Tracers mass and changes
      REAL*8, DIMENSION(LM,NTM)      :: TMOLD,  DTM,  DTMR,TM1
      REAL*8, DIMENSION(NMOM,LM,NTM) :: TMOMOLD,DTMOM,DTMOMR
      REAL*8, DIMENSION(NTM)      :: TMP,  TMPMAX,  TMDN, TENV
      REAL*8, DIMENSION(NMOM,NTM) :: TMOMP,TMOMPMAX,TMOMDN
!@var TPOLD saved plume temperature after condensation for tracers
!@+   (this is slightly different from TPSAV)
      REAL*8, DIMENSION(LM)       :: TPOLD
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM)      :: TRPRCP
      REAL*8, DIMENSION(NTM,LM)   :: TRCOND
!@var FQCONDT fraction of tracer that condenses
!@var FQEVPT  fraction of tracer that evaporates (in downdrafts)
!@var FPRCPT fraction of tracer that evaporates (in net re-evaporation)
!@var FWASHT  fraction of tracer scavenged by below-cloud precipitation
      REAL*8 :: FQCONDT, FWASHT, FPRCPT, FQEVPT
!@var WMXTR available water mixing ratio for tracer condensation ( )?
!@var b_beta_DT precipitating gridbox fraction from lowest precipitating
!@+   layer. The name was chosen to correspond to Koch et al. p. 23,802.
!@var precip_mm precipitation (mm) from the grid box above for washout
      REAL*8 WMXTR, b_beta_DT, precip_mm
c for tracers in general, added by Koch
      REAL*8 THLAW,TR_LEF,TMFAC,TR_LEFT(ntm),CLDSAVT
!@var TR_LEF limits precurser dissolution following sulfate formation
!@var THLAW Henry's Law determination of amount of tracer dissolution
!@var TMFAC used to adjust tracer moments
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      REAL*8 TMP_SUL(LM,NTM)
c for sulfur chemistry
!@var WA_VOL Cloud water volume (L). Used by GET_SULFATE.
      REAL*8 WA_VOL
      REAL*8, DIMENSION(NTM) ::SULFOUT,SULFIN,SULFINC
#endif
      REAL*8 DTSUM,HEFF
#endif
#endif

      REAL*8, DIMENSION(LM) ::
     *     DM,COND,CDHEAT,CCM,SM1,QM1,DMR,ML,SMT,QMT,TPSAV,DDM,CONDP
#ifdef CLD_AER_CDNC
     *     ,CONDPC
#endif
      REAL*8 :: CONDGP,CONDIP
!@var DM change in air mass
!@var COND,CONDGP,CONDIP condensate
!@var CDHEAT heating due to condensation
!@var CCM convective plume mass
!@var SM1, QM1 dummy variables
!@var DMR change of air mass
!@var ML layer air mass
!@var SMT, QMT dummy variables
!@var TPSAV array to save plume temperature
!@var DDM downdraft mass
!@param FCLW fraction of condensate in plume that remains as CLW
      REAL*8 :: FCLW

!@var IERRT,LERRT error reports from advection
      INTEGER :: IERRT,LERRT
      INTEGER, INTENT(IN) :: i_debug,j_debug
#ifdef CLD_AER_CDNC
      INTEGER, PARAMETER :: SNTM=17+ntm_soa/2  !for tracers for CDNC
#endif
      INTEGER LDRAFT,LMAX,LMIN,MCCONT,MAXLVL
     *     ,MINLVL,ITER,IC,LFRZ,NSUB,LDMIN
!@var LDRAFT the layer the downdraft orginates
!@var LEVAP the layer evaporation of precip starts
!@var LMAX, LMIN the lowest, the highest layer of a convective event
!@var MCCONT integer to count convective events
!@var MAXLVL, MINLVL the lowest, the highest layer of convective events
!@var ITER number for iteration
!@var IC integer for cloud types
!@var LFRZ freezing level
!@var nsub = LMAX - LMIN + 1
!@var LDMIN the lowest layer the downdraft drops
      REAL*8 TERM1,FMP0,SMO1
     *     ,QMO1,SMO2,QMO2,SDN,QDN,SUP,QUP,SEDGE,QEDGE,WMDN,WMUP,SVDN
     *     ,SVUP,WMEDG,SVEDG,DMSE,FPLUME,DFP,FMP2,FRAT1,FRAT2,SMN1
     *     ,QMN1,SMN2,QMN2,SMP,QMP,TP,GAMA,DQSUM,TNX,TNX1
     *     ,DQ,DMSE1,FCTYPE,BETAU,ALPHAU
     *     ,CDHDRT,DDRAFT,DELTA
     *     ,ALPHA,BETA,CDHM,CDHSUM,CLDM,CLDREF,CONSUM,DQEVP
     *     ,DQRAT,EPLUME,ETADN,ETAL1,EVPSUM,FCDH
     *     ,FCDH1,FCLD,FCLOUD,FDDL,FDDP,FENTR,FENTRA,FEVAP,FLEFT
     *     ,FQCOND,FQEVP,FPRCP,FSEVP,FSSUM,HEAT1
     *     ,PRCP
     *     ,QMDN,QMIX,QMPMAX,QMPT,QNX,QSATC,QSATMP
     *     ,RCLD,RCLDE,SLH,SMDN,SMIX,SMPMAX,SMPT,SUMAJ
     *     ,SUMDP,DDRUP,EDRAFT
     *     ,TOLD,TOLD1,TEMWM,TEM,WTEM,WCONST,WORK
     *     ,FCONV_tmp,FSUB_tmp,FSSL_tmp
     *     , MNdO,MNdL,MNdI,MCDNCW,MCDNCI   !Menon
#ifdef CLD_AER_CDNC
     *     ,MCDNO1,MCDNL1,CDNCB,fcnv,ATEMP,VVEL
     *     ,DSGL(LM,SNTM),DSS(SNTM)
     *     ,Repsis,Repsi,Rbeta,RCLD_C
#ifdef BLK_2MOM
#ifdef TRACERS_AMP
      real*8                    :: ncaero (nmodes)
C** To get the model to work with right dependencies uncomment the following declaration
C** and for nactc dimension declaration
CC* Comment the stmt in MODULE CLOUDS beginning part  USE AMP_AEROSOL, only: NACTC,NAERC
C** Once you save the right dependency (cp amp_aerosol.mod amp_aerosol.modsave)
C** (this needs to be done only once)
C** recomment the declaration below and for nactc and uncomment the USE AMP_AEROSOL stmt
C** and after gmake clean vclean
C** you will need to cp amp_aerosol.modsave as amp_aerosol.mod
C** This fix is till MATRIX dependencies solved
c     real*8,dimension(lm,nmodes)   :: naerc
      integer                   ::nm
#endif
#endif
#endif
!@var TERM1 contribution to non-entraining convective cloud
!@var FMP0 non-entraining convective mass
!@var SMO1,QMO1,SMO2,QMO2,SDN,QDN,SUP,QUP,SEDGE,QEDGE dummy variables
!@var WMDN,WMUP,SVDN,SVUP,WMEDG,SVEDG,DDRUP dummy variables
!@var DMSE difference in moist static energy
!@var FPLUME fraction of convective plume
!@var DFP an iterative increment
!@var FMP2,FRAT1,FRAT2,SMN1,QMN1,SMN2,QMN2 dummy variables
!@var SMP, QMP plume's SM, QM
!@var TP plume's temperature (K)
!@var GAMA,DQSUM,TNX,DQ,CONSUM,BETAU,ALPHAU dummy variables
!@var DMSE1 difference in moist static energy
!@var FCTYPE fraction for convective cloud types
!@var CDHDRT,ALPHA,BETA,CDHM,CDHSUM,CLDREF dummy variables
!@var DDRAFT downdraft mass
!@var DELTA fraction of plume stays in the layer
!@var CLDM subsidence due to convection
!@var DQEVP amount of condensate that evaporates in downdrafts
!@var DQRAT fraction for the condensate to evaporate
!@var EPLUME air mass of entrainment
!@var ETADN fraction of the downdraft
!@var ETAL1 fractional entrainment rate
!@var EVPSUM,FCDH,FCDH1,FCLD,FCLOUD,FDDL,FDDP dummy variables
!@var FENTR fraction of entrainment
!@var FENTRA fraction of entrainment
!@var FEVAP fraction of air available for precip evaporation
!@var FLEFT fraction of plume after removing downdraft mass
!@var FQCOND fraction of water vapour that condenses in plume
!@var FQEVP fraction of water vapour that evaporates in downdraft
!@var FPRCP fraction of evaporated precipitation
!@var FSEVP fraction of energy lost to evaporate
!@var FSSUM fraction of energy lost to evaporate
!@var HEAT1 heating due to phase change
!@var PRHEAT energy of condensate
!@var PRCP precipipation
!@var QMDN,QMIX,SMDN,SMIX dummy variables
!@var QMPMAX,SMPMAX detrainment of QMP, SMP
!@var QMPT,SMPT dummy variables
!@var QNX,SUMAJ,SUMDP dummy variables
!@var QSATC saturation vapor mixing ratio
!@var QSATMP plume's saturation vapor mixing ratio
!@var RCLD,RCLDE cloud particle's radius, effective radius
#ifdef CLD_AER_CDNC
!@var MCDNCW,MCDNCI cloud droplet # for warm,cold moist conv clouds (cm^-3)
#endif
!@var SLH LHX/SHA
!@var EDRAFT entrainment into downdrafts
!@var TOLD,TOLD1 old temperatures
!@var TEMWM,TEM,WTEM,WCONST dummy variables
!@var WORK work done on convective plume

      LOGICAL MC1  !@var MC1 true for the first convective event

      REAL*8,  PARAMETER :: CK1 = 1.       !@param CK1 a tunning const.
!@param RHOG,RHOIP density of graupel and ice particles
!@param ITMAX max iteration index
!@param CN0 constant use in computing FLAMW, etc
!@param PN tuning exponential for computing WV
      REAL*8,  PARAMETER :: CN0=8.d6,PN=1.d0
#ifdef CLD_AER_CDNC
      REAL*8 RHO   ! air density
!CN0 is the No parameter in the Marshall-Palmer distribution
#endif
      REAL*8,  PARAMETER :: RHOG=400.,RHOIP=100.
      INTEGER,  PARAMETER :: ITMAX=50
!@var FLAMW,FLAMG,FLAMI lamda for water, graupel and ice, respectively
!@var WMAX specified maximum convective vertical velocity
!@var WV convetive vertical velocity
!@var VT precip terminal velocity
!@var DCW,DCG,DCI critical cloud particle sizes
!@var FG, FI fraction for graupel and ice
!@var FITMAX set to ITMAX
!@var CONDMU convective condensate in Kg/m^3
      REAL*8 :: FLAMW,FLAMG,FLAMI,WMAX,WV,DCW,DCG,DCI,FG,FI,FITMAX,DDCW
     *     ,VT,CONDMU
      INTEGER K,L,N  !@var K,L,N loop variables
      INTEGER ITYPE  !@var convective cloud types
!@var IERR,LERR error reports from advection
      INTEGER, INTENT(OUT) :: IERR,LERR
!@var DUM, DVM changes of UM,VM
      REAL*8, DIMENSION(KMAX,LM) :: DUM,DVM

      REAL*8 THBAR  !@var THBAR virtual temperature at layer edge
!@var BELOW_CLOUD logical- is the current level below cloud?
      LOGICAL BELOW_CLOUD
C****
C**** MOIST CONVECTION
C****
C**** CONVECTION OCCURS AT THE LOWEST MOIST CONVECTIVELY UNSTABLE
C**** LEVEL AND CONTINUES UNTIL A STABLE LAYER PAIR IS REACHED.  RE-
C**** EVAPORATION AND PRECIPITATION ARE COMPUTED AT THE END OF THIS
C**** CYCLE.  THE WHOLE PROCESS MAY BE REPEATED FROM A NEW LOWEST
C**** UNSTABLE LEVEL.
C****
      ierr=0 ; lerr=0
      LMCMIN=0
      LMCMAX=0
      MCCONT=0
      FMC1=0.
      FSSL=1.
      FITMAX=ITMAX
      RCLDX=radius_multiplier
! allow variable entrainment: scale by ENTCON=0.2 (consistency w/ devel model)
      CONTCE1=entrainment_cont1*5d0
C**** initiallise arrays of computed ouput
      TAUMCL=0
      SVWMXL=0
      SVLATL=0
      VSUBL=0
      PRECNVL=0
      CLDMCL=0
      CLDSLWIJ=0
      CLDDEPIJ=0
      PRCPMC=0.
      TPSAV=0
      CSIZEL=RWCLDOX*10.*(1.-PEARTH)+10.*PEARTH ! droplet rad in stem
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
C**** save initial values (which will be updated after subsid)
      SM1=SM
      QM1=QM
      FSUB=0.
      FCONV=0.
      FMCL=0.
      VLAT=LHE
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
#endif
C**** OUTER MC LOOP OVER BASE LAYER
      DO 600 LMIN=1,LMCM-1
      MAXLVL=0
      MINLVL=LM
C****
C**** COMPUTE THE CONVECTIVE MASS OF THE NON-ENTRAINING PART
C****
      TERM1=-10.*CK1*SDL(LMIN+1)*BYGRAV
      FMP0=TERM1*XMASS
      IF(FMP0.LE.0.) FMP0=0.
C**** CREATE A PLUME IN THE BOTTOM LAYER
C****
C**** ITERATION TO FIND FPLUME WHICH RESTORES THE ATM TO NEUTRAL STATE
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
      IF(DMSE.GT.-1d-10) GO TO 600
C****
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
      IF(DQSUM.LE.0.) GO TO 205
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
      IF(DQSUM.LE.0.) GO TO 205
      MCLOUD=FEVAP*AIRM(LMIN)
      TNX=SMO1*PLK(LMIN)*BYAM(LMIN)
      QNX=QMO1*BYAM(LMIN)
      QSATC=QSAT(TNX,LHX,PL(LMIN))
      DQ=MCLOUD*(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX)) ! new func'n
      IF(DQ.GT.DQSUM) DQ=DQSUM
      SMN1=SMN1-SLH*DQ/PLK(LMIN)
      QMN1=QMN1+DQ
  205 SDN=SMN1*BYAM(LMIN)
      SUP=SMN2*BYAM(LMIN+1)
      SEDGE=THBAR(SUP,SDN)
      QDN=QMN1*BYAM(LMIN)
      QUP=QMN2*BYAM(LMIN+1)
      SVDN=SDN*(1.+DELTX*QDN-WMDN)
      SVUP=SUP*(1.+DELTX*QUP-WMUP)
      QEDGE=.5*(QUP+QDN)
      SVEDG=SEDGE*(1.+DELTX*QEDGE-WMEDG)
      DMSE1=(SVUP-SVEDG)*PLK(LMIN+1)+(SVEDG-SVDN)*PLK(LMIN)+
     *  SLHE*(QSAT(SUP*PLK(LMIN+1),LHX,PL(LMIN+1))-QDN)
      IF (ABS(DMSE1).LE.1.d-3) GO TO 411
      IF(DMSE1.GT.1.d-3) FPLUME=FPLUME-DFP
      IF(DMSE1.LT.-1.d-3) FPLUME=FPLUME+DFP
      END DO
  411 IF(FPLUME.LE..001) GO TO 600
C****
C**** ITERATION THROUGH CLOUD TYPES
C****
      ITYPE=2                        ! always 2 types of clouds:
C     IF(LMIN.LE.2) ITYPE=2          ! entraining & non-entraining
      FCTYPE=1.

      DO 570 IC=1,ITYPE
C**** INITIALLISE VARIABLES USED FOR EACH TYPE
      DO L=1,LM
        COND(L)=0.
        CDHEAT(L)=0.
        CONDP(L)=0.
        DM(L)=0.
        DMR(L)=0.
        CCM(L)=0.
        DDM(L)=0.
      END DO
      DUM(1:KMAX,:)=0.
      DVM(1:KMAX,:)=0.
      DSM(:) = 0.
      DSMOM(:,:) = 0.
      DSMR(:) = 0.
      DSMOMR(:,:) = 0.
      DQM(:) = 0.
      DQMOM(:,:) = 0.
      DQMR(:) = 0.
      DQMOMR(:,:) = 0.
#ifdef TRACERS_ON
      DTM(:,1:NTX) = 0.
      DTMOM(:,:,1:NTX) = 0.
      DTMR(:,1:NTX) = 0.
      DTMOMR(:,:,1:NTX) = 0.
      TPOLD = 0.
#endif
#ifdef TRACERS_WATER
      TRCOND = 0.
#endif
      MC1=.FALSE.
      IF (IC.EQ.1) THEN
        WMAX=2.
        IF(PLAND.GE..5) WMAX=5.
      ELSE
        WMAX=1.
        IF(PLAND.GE..5) WMAX=2.5
      ENDIF
      LHX=LHE
      MPLUME=MIN(AIRM(LMIN),AIRM(LMIN+1))
      IF(MPLUME.GT.FMP2) MPLUME=FMP2
      IF(ITYPE.EQ.2) THEN
      FCTYPE=1.
      IF(MPLUME.GT.FMP0) FCTYPE=FMP0/MPLUME
      IF(IC.EQ.2) FCTYPE=1.-FCTYPE
      IF(FCTYPE.LT.0.001) GO TO 570
      END IF
      MPLUM1=MPLUME
      MPLUME=MPLUME*FCTYPE

C**** calculate subsiding fraction here. Then we can use FMC1 from the
C**** beginning. The analagous arrays are only set if plume is actually
C**** moved up.
      IF (MCCONT.le.0) THEN
         FCONV_tmp=MPLUM1/AIRM(LMIN+1)
         IF(FCONV_tmp.GT.1.d0) FCONV_tmp=1.d0
         FSUB_tmp=1.d0+(AIRM(LMIN+1)-100.d0)/200.d0
         IF(FSUB_tmp.GT.1.d0/(FCONV_tmp+1.d-20)-1.d0)
     *        FSUB_tmp=1.d0/(FCONV_tmp+1.d-20)-1.d0
         IF(FSUB_tmp.LT.1.d0) FSUB_tmp=1.d0
         IF(FSUB_tmp.GT.5.d0) FSUB_tmp=5.d0
         FSSL_tmp=1.d0-(1.d0+FSUB_tmp)*FCONV_tmp
         IF(FSSL_tmp.LT.CLDMIN) FSSL_tmp=CLDMIN
         IF(FSSL_tmp.GT.1.d0-CLDMIN) FSSL_tmp=1.d0-CLDMIN
         FMC1=1.d0-FSSL_tmp+teeny
      ELSE
C**** guard against possibility of too big a plume
        MPLUME=MIN(0.95d0*AIRM(LMIN)*FMC1,MPLUME)
      END IF

C**** adjust MPLUME to take account of restricted area of subsidence
C**** (i.e. MPLUME is now a greater fraction of the relevant airmass.
      MPLUME=MIN( MPLUME/FMC1,
     *   AIRM(LMIN)*0.95d0*QM(LMIN)/(QMOLD(LMIN) + teeny) )
      IF(MPLUME.LE..001*AIRM(LMIN)) GO TO 570
      FPLUME=MPLUME*BYAM(LMIN)
      SMP  =  SMOLD(LMIN)*FPLUME
      SMOMP(xymoms)=SMOMOLD(xymoms,LMIN)*FPLUME
      QMP  =  QMOLD(LMIN)*FPLUME
      QMOMP(xymoms)=QMOMOLD(xymoms,LMIN)*FPLUME
      TPSAV(LMIN)=SMP*PLK(LMIN)/MPLUME
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
      ENDDO
C****
C**** RAISE THE PLUME TO THE TOP OF CONVECTION AND CALCULATE
C**** ENTRAINMENT, CONDENSATION, AND SECONDARY MIXING
C****
      CDHSUM=0.
      CDHDRT=0.
      ETADN=0.
      LDRAFT=LM
      EVPSUM=0.
      DDRAFT=0.
      LFRZ=0
      LMAX=LMIN
 220  L=LMAX+1
C**** TEST FOR SUFFICIENT AIR, MOIST STATIC STABILITY AND ENERGY
C     IF(L.GT.LMIN+1.AND.SDL(L).GT.0.) GO TO 340
      IF(MPLUME.LE..001*AIRM(L)) GO TO 340
      SDN=SMP/MPLUME
      SUP=SM1(L)*BYAM(L)
      QDN=QMP/MPLUME
      QUP=QM1(L)*BYAM(L)
      WMDN=0.
      WMUP=WML(L)
      SVDN=SDN*(1.+DELTX*QDN-WMDN)
      SVUP=SUP*(1.+DELTX*QUP-WMUP)
      IF(LMAX.GT.LMIN) THEN
      SEDGE=THBAR(SUP,SDN)
      QEDGE=.5*(QUP+QDN)
      WMEDG=.5*(WMUP+WMDN)
      SVEDG=SEDGE*(1.+DELTX*QEDGE-WMEDG)
      LHX=LHE
      DMSE=(SVUP-SVEDG)*PLK(L)+(SVEDG-SVDN)*PLK(L-1)+
     *  SLHE*(QSAT(SUP*PLK(L),LHX,PL(L))-QDN)
      IF(DMSE.GT.-1d-10) GO TO 340
      END IF
      IF(PLK(L-1)*(SVUP-SVDN)+SLHE*(QUP-QDN).GE.0.) GO TO 340
C**** TEST FOR CONDENSATION ALSO DETERMINES IF PLUME REACHES UPPER LAYER
      TP=SMP*PLK(L)/MPLUME
      TPSAV(L)=TP
      IF(TPSAV(L-1).GE.TF.AND.TPSAV(L).LT.TF) LFRZ=L-1
      IF(TP.LT.TI) LHX=LHS
      QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
      IF(QMP.LT.QSATMP) GO TO 340
      IF(TP.GE.TF.OR.LHX.EQ.LHS) GO TO 290
      LHX=LHS
      QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
  290 CONTINUE
C**** DEFINE VLAT TO AVOID PHASE DISCREPANCY BETWEEN TWO PLUMES
      IF (VLAT(L).EQ.LHS) LHX=LHS
      VLAT(L)=LHX
      SLH=LHX*BYSHA
      MCCONT=MCCONT+1
      IF(MCCONT.EQ.1) MC1=.TRUE.
C     IF(MC1.AND.L.EQ.LMIN+1) THEN
C        FCONV(L)=FCONV_tmp   ! these are set here but do not make
C        FSSL(L)=FSSL_tmp     ! much sense at the moment...
C        FSUB(L)=FSUB_tmp
C     ENDIF
C****
C**** DEPOSIT PART OF THE PLUME IN LOWER LAYER
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
        ENDDO

#ifdef TRACERS_ON
        DTM(L-1,1:NTX) = DTM(L-1,1:NTX)+DELTA*TMP(1:NTX)
        DTMOM(xymoms,L-1,1:NTX)=DTMOM(xymoms,L-1,1:NTX)+DELTA
     *       *TMOMP(xymoms,1:NTX)
        TMP(1:NTX) = TMP(1:NTX)*(1.-DELTA)
        TMOMP(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)*(1.-DELTA)
#endif
      END IF
C****
C**** CONVECTION IN UPPER LAYER   (WORK DONE COOLS THE PLUME)
C****
      WORK=MPLUME*(SUP-SDN)*(PLK(L-1)-PLK(L))/PLK(L-1)
C     SMP=SMP-WORK
      DSM(L-1)=DSM(L-1)-WORK
      CCM(L-1)=MPLUME
      WV=WMAX
      IF(PL(L).GT.700..AND.PLE(1).GT.700.)
     *  WV=WMAX*(PLE(1)-PL(L))/(PLE(1)-700.)
      IF(PL(L).LT.400.) WV=WMAX*((PL(L)-PLE(LM+1))/(400.-PLE(LM+1)))**PN
      DCG=((WV/19.3)*(PL(L)/1000.)**.4)**2.7
      DCG=MIN(DCG,1.D-2)
      DCI=((WV/1.139)*(PL(L)/1000.)**.4)**9.09
      DCI=MIN(DCI,1.D-2)
C****
C**** ENTRAINMENT
C****
      IF(IC.EQ.2.OR.(IC.EQ.1.AND.PL(L).GE.800.)) THEN
      FENTR=CONTCE1*ETAL(L)*FPLUME   ! optional scaling
      IF(FENTR+FPLUME.GT.1.) FENTR=1.-FPLUME
      IF(FENTR.LT.teeny) GO TO 293
      ETAL1=FENTR/(FPLUME+teeny)
      FPLUME=FPLUME+FENTR
      EPLUME=MPLUME*ETAL1
C**** Reduce EPLUME so that mass flux is less than mass in box
      IF (EPLUME.GT.AIRM(L)*0.975d0-MPLUME) THEN
        EPLUME=AIRM(L)*0.975d0-MPLUME
      END IF
      MPLUME=MPLUME+EPLUME
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
      DO K=1,KMAX
         UMP(K)=UMP(K)+U_0(K,L)*EPLUME
         DUM(K,L)=DUM(K,L)-U_0(K,L)*EPLUME
         VMP(K)=VMP(K)+V_0(K,L)*EPLUME
         DVM(K,L)=DVM(K,L)-V_0(K,L)*EPLUME
      ENDDO
  293 CONTINUE
      END IF
C****
C**** CHECK THE DOWNDRAFT POSSIBILITY
C****
      IF(L-LMIN.LE.1) GO TO 291
      IF(ETADN.GT.1d-10) GO TO 291
      SMIX=.5*(SUP+SMP/MPLUME)
      QMIX=.5*(QUP+QMP/MPLUME)
C     WMIX=.5*(WMUP+WMDN)
C     SVMIX=SMIX*(1.+DELTX*QMIX-WMIX)
C     DMMIX=(SVUP-SVMIX)*PLK(L)+
C    *  SLHE*(QSAT(SUP*PLK(L),LHX,PL(L))-QMIX)
C     IF(DMMIX.LT.1d-10) GO TO 291
      IF(SMIX.GE.SUP) GO TO 291
      IF(PL(L).GT.700.) GO TO 291
      LDRAFT=L
      ETADN=BY3
      FLEFT=1.-.5*ETADN
      DDRAFT=ETADN*MPLUME
      FDDP = .5*DDRAFT/MPLUME
      FDDL = .5*DDRAFT*BYAM(L)
      MPLUME=FLEFT*MPLUME
      SMDN=DDRAFT*SMIX         ! = SM(L)*FDDL +  SMP*FDDP
      SMOMDN(xymoms)=SMOM(xymoms,L)*FDDL +  SMOMP(xymoms)*FDDP
      SMP=FLEFT*SMP
      SMOMP(xymoms)= SMOMP(xymoms)*FLEFT
      QMDN=DDRAFT*QMIX         ! = QM(L)*FDDL +  QMP*FDDP
      QMOMDN(xymoms)=QMOM(xymoms,L)*FDDL +  QMOMP(xymoms)*FDDP
      QMP=FLEFT*QMP
      QMOMP(xymoms)= QMOMP(xymoms)*FLEFT
      DMR(L) = DMR(L)-.5*DDRAFT
      DSMR(L)=DSMR(L)-.5*DDRAFT*SUP        ! = DSM(L)-SM(L)*FDDL
      DSMOMR(:,L)=DSMOMR(:,L) - SMOM(:,L)*FDDL
      DQMR(L)=DQMR(L)-.5*DDRAFT*QUP        ! = DQM(L)-QM(L)*FDDL
      DQMOMR(:,L)=DQMOMR(:,L) - QMOM(:,L)*FDDL
#ifdef TRACERS_ON
      Tmdn(1:NTX) = tm(l,1:NTX)*fddl+Tmp(1:NTX)*fddp
      tmomdn(xymoms,1:NTX) = tmom(xymoms,l,1:NTX)*fddl+
     *     tmomp(xymoms,  1:NTX)*fddp
      dtmr    (l,1:NTX) = dtmr    (l,1:NTX)-fddl *tm    (l,1:NTX)
      dtmomr(:,l,1:NTX) = dtmomr(:,l,1:NTX)-fddl *tmom(:,l,1:NTX)
      Tmp         (1:NTX) = Tmp         (1:NTX)*fleft
      tmomp(xymoms,1:NTX) = tmomp(xymoms,1:NTX)*fleft
#endif
      DO K=1,KMAX
         UMDN(K)=.5*(ETADN*UMP(K)+DDRAFT*U_0(K,L))
         UMP(K)=UMP(K)*FLEFT
         DUM(K,L)=DUM(K,L)-.5*DDRAFT*U_0(K,L)
         VMDN(K)=.5*(ETADN*VMP(K)+DDRAFT*V_0(K,L))
         VMP(K)=VMP(K)*FLEFT
         DVM(K,L)=DVM(K,L)-.5*DDRAFT*V_0(K,L)
      ENDDO
C****
C**** CONDENSE VAPOR IN THE PLUME AND ADD LATENT HEAT
C****
  291 DQSUM=0.
      SMPT=SMP
      QMPT=QMP
      DO 292 N=1,3
      TP=SMP*PLK(L)/MPLUME
      QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
      GAMA=SLH*QSATMP*DQSATDT(TP,LHX)/MPLUME
      DQ=(QMP-QSATMP)/(1.+GAMA)
      SMP=SMP+SLH*DQ/PLK(L)
      QMP=QMP-DQ
  292 DQSUM=DQSUM+DQ
#ifdef TRACERS_ON
C**** save plume temperature after possible condensation
      TPOLD(L)=SMP*PLK(L)/MPLUME
#endif
      FQCOND = 0
      IF(DQSUM.GE.0.) THEN
        IF (QMPT.gt.teeny) FQCOND = DQSUM/QMPT
        QMOMP(xymoms) =  QMOMP(xymoms)*(1.-FQCOND)
      ELSE  ! no change
        DQSUM=0.
        SMP=SMPT
        QMP=QMPT
      END IF
      COND(L)=DQSUM
      CONDMU=100.*COND(L)*PL(L)/(CCM(L-1)*TL(L)*RGAS)
      FLAMW=(1000.d0*PI*CN0/(CONDMU+teeny))**.25
      FLAMG=(400.d0*PI*CN0/(CONDMU+teeny))**.25
      FLAMI=(100.d0*PI*CN0/(CONDMU+teeny))**.25
#ifdef CLD_AER_CDNC
!@auth Menon  saving aerosols mass for CDNC prediction
      DO N=1,SNTM
        DSS(N)=1.d-10
        DSGL(L,N)=1.d-10
      ENDDO
C** Here we change convective precip due to aerosols
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
       DO N=1,NTX
        select case (trname(ntix(n)))
         case('SO4')
         DSGL(L,1)=tm(l,n)  !n=19
         DSS(1) = DSGL(L,1)
         case('seasalt1')
         DSGL(L,2)=tm(l,n)  !n=21
         DSS(2) = DSGL(L,2)
         case('seasalt2')
         DSGL(L,3)=tm(l,n)  !n=22
         DSS(3) = DSGL(L,3)
         case('OCIA')
         DSGL(L,4)=tm(l,n)    !n=27
         DSS(4) = DSGL(L,4)
         case('OCB')
         DSGL(L,5)=tm(l,n)  !n=28
         DSS(5) = DSGL(L,5)
         case('BCIA')
         DSGL(L,6)=tm(l,n)  !n=24
         DSS(6) = DSGL(L,6)
         case('BCB')
         DSGL(L,7)=tm(l,n)  !n=25
         DSS(7) = DSGL(L,7)
         case('OCII')
         DSGL(L,8)=tm(l,n)  !n=26
         DSS(8) = DSGL(L,8)
         case('BCII')
         DSGL(L,9)=tm(l,n)  !n=23
         DSS(9) = DSGL(L,9)
#ifdef TRACERS_DUST
         case('Clay')
         DSGL(L,10)=tm(l,n)  !n=23
         DSS(10) = max(1.d-10,DSGL(L,10))
         case('Silt1')
         DSGL(L,11)=tm(l,n)  !n=23
         DSS(11) = max(1.d-10,DSGL(L,11))
         case('Silt2')
         DSGL(L,12)=tm(l,n)  !n=23
         DSS(12) = max(1.d-10,DSGL(L,12))
         case('Silt3')
         DSGL(L,13)=tm(l,n)  !n=23
         DSS(13) = max(1.d-10,DSGL(L,13))
#endif
#ifdef TRACERS_NITRATE
         case('NO3p')
         DSGL(L,14)=tm(l,n)  !n=23
         DSS(14) = DSGL(L,14)
#endif
#ifdef TRACERS_HETCHEM
c!*** Here are dust particles coated with sulfate
         case('SO4_d1')
         DSGL(L,15)=tm(l,n)  !n=20
         DSS(15) = max(1.d-10,DSGL(L,15))
         case('SO4_d2')
         DSGL(L,16)=tm(l,n)  !n=21
         DSS(16) = max(1.d-10,DSGL(L,16))
         case('SO4_d3')
         DSGL(L,17)=tm(l,n)  !n=22
         DSS(17) = max(1.d-10,DSGL(L,17))
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
#endif
C** Use MATRIX AMP_actv to decide what the aerosol number conc. is
#ifdef BLK_2MOM
#ifdef TRACERS_AMP
        do nm=1,nmodes
         ncaero(nm)=naerc(l,nm)*1.d-6
c         if(naerc(l,nm).gt.1.d-30) write(6,*)"mat",ncaero(nm),nm
        enddo
       CALL GET_CC_CDNC_MX(L,nmodes,ncaero,MCDNL1,MCDNO1)
#else
C** This is for the old mass to number calculations nc. is
       CALL GET_CC_CDNC(L,AIRM(L),DXYPIJ,PL(L),TL(L),DSS,MCDNL1,MCDNO1)

#endif
#endif
       MNdO=MCDNO1
       MNdL=MCDNL1
c      if(MNdL.gt.10.) write(6,*)"Here is CDNC in E1",MCDNO1,MCDNL1,L
c    *,ncaero(1),ncaero(2),ncaero(3)
      MNdI = 0.06417127d0
      MCDNCW=MNdO*(1.-PEARTH)+MNdL*PEARTH
      MCDNCI=MNdI
      if (MCDNCW.le.20.d0) MCDNCW=20.d0     !set min CDNC, sensitivity test
      if (MCDNCW.ge.2000.d0) MCDNCW=2000.d0     !set max CDNC, sensitivity test
c     write(6,*)"CCONV",MCDNCW,MNdO,MNdL,L,LMIN
C** Using the Liu and Daum paramet, Nature, 2002, Oct 10, Vol 419
C** for spectral dispersion effects on droplet size distribution
C** central value of 0.003 for alfa:Rotstayn&Daum, J.Clim,2003,16,21, Nov 2003.
      Repsi=1.d0 - 0.7d0*exp(-0.003d0*MCDNCW)
      Repsis=Repsi*Repsi
      Rbeta=(((1.d0+2.d0*Repsis)**0.667d0))/((1.d0+Repsis)**by3)
c     RCLDE =RCLD*Rbeta
      RCLD_C=14.d0/Rbeta       !set Reff to threshold size =14 um (Rosenfeld)
#endif

      IF (TP.GE.TF) THEN ! water phase
        DDCW=6d-3/FITMAX
        IF(PLAND.LT..5) DDCW=1.5d-3/FITMAX
        DCW=0.
        DO ITER=1,ITMAX-1
          VT=(-.267d0+DCW*(5.15D3-DCW*(1.0225D6-7.55D7*DCW)))*
     *       (1000./PL(L))**.4d0
          IF(VT.GE.0..AND.ABS(VT-WV).LT..3) EXIT
          IF(VT.GT.WMAX) EXIT
          DCW=DCW+DDCW
        END DO
        CONDP(L)=RHOW*(PI*by6)*CN0*EXP(-FLAMW*DCW)*
     *     (DCW*DCW*DCW/FLAMW+3.*DCW*DCW/(FLAMW*FLAMW)+
     *     6.*DCW/(FLAMW*FLAMW*FLAMW)+6./FLAMW**4)
#ifdef CLD_AER_CDNC
        RHO=1d2*PL(L)/(RGAS*TL(L)) !air density in kg/m3
C** Here we calculate condensate converted to precip
C** only if drop size (Reff) exceeds 14 um
C** Conver Rvol to m and CDNC to m from um and cm, respectively
C** Precip condensate is simply the LWC
        CONDPC(L)=4.d0*by3*PI*RHOW*((RCLD_C*1.d-6)**3)*1.d6*MCDNCW/RHO
        IF(CONDMU.gt.CONDPC(L))  then
          CONDP(L)=CONDMU-CONDPC(L) !
        ELSE
          CONDP(L)=0.d0
        ENDIF
c     if (CONDP(L).lt.0.d0)
c    *write(6,*)"Mup",CONDP(L),CONDPC(L),CONDMU,DCW,MCDNCW,RCLD_C,L
#endif
      ELSE IF (TP.LE.TI) THEN ! pure ice phase
        CONDP(L)=RHOIP*(PI*by6)*CN0*EXP(-FLAMI*DCI)*
     *    (DCI*DCI*DCI/FLAMI+3.*DCI*DCI/(FLAMI*FLAMI)+
     *    6.*DCI/(FLAMI*FLAMI*FLAMI)+6./FLAMI**4)
      ELSE ! mixed phase
        FG=(TP-TF+40.)*0.025d0
        FI=1.-FG
        CONDIP=RHOIP*(PI*by6)*CN0*EXP(-FLAMI*DCI)*
     *    (DCI*DCI*DCI/FLAMI+3.*DCI*DCI/(FLAMI*FLAMI)+
     *    6.*DCI/(FLAMI*FLAMI*FLAMI)+6./FLAMI**4)
        CONDGP=RHOG*(PI*by6)*CN0*EXP(-FLAMG*DCG)*
     *    (DCG*DCG*DCG/FLAMG+3.*DCG*DCG/(FLAMG*FLAMG)+
     *    6.*DCG/(FLAMG*FLAMG*FLAMG)+6./FLAMG**4)
        CONDP(L)=FG*CONDGP+FI*CONDIP
      ENDIF
c convert condp to the same units as cond
      CONDP(L)=.01d0*CONDP(L)*CCM(L-1)*TL(L)*RGAS/PL(L)
#ifdef TRACERS_WATER
C**** CONDENSING TRACERS
      WMXTR=DQSUM*BYAM(L)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      WA_VOL=COND(L)*1.d2*BYGRAV*DXYPIJ
      DO N=1,NTX
      select case (trname(ntix(n)))
      case('SO2','SO4','H2O2_s','H2O2')
      if (trname(ntix(n)).eq."H2O2" .and. coupled_chem.eq.0) goto 400
      if (trname(ntix(n)).eq."H2O2_s" .and. coupled_chem.eq.1) goto 400

        IF (FPLUME.GT.teeny) then
          TMP_SUL(L,N)=TMP(N)/FPLUME
        else
          TMP_SUL(L,N)=0.
        ENDIF

 400   CONTINUE

      end select

      END DO
      CALL GET_SULFATE(L,TPOLD(L),FPLUME,WA_VOL,WMXTR,SULFIN,
     *     SULFINC,SULFOUT,TR_LEFT,TMP_SUL,TRCOND(1,L),
     *     AIRM,LHX,
     *     DT_SULF_MC(1,L),CLDSAVT)
#endif
      DO N=1,NTX
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      select case (trname(ntix(n)))
      case('SO2','SO4','H2O2_s','H2O2')
      if (trname(ntix(n)).eq."H2O2" .and. coupled_chem.eq.0) goto 401
      if (trname(ntix(n)).eq."H2O2_s" .and. coupled_chem.eq.1) goto 401

c first apply chemistry
c removal of precursers
        TMP(N)=TMP(N)*(1.+SULFIN(N))
        TMOMP(xymoms,N)= TMOMP(xymoms,N)*(1.+SULFIN(N))
c formation of sulfate
        TRCOND(N,L) = TRCOND(N,L)+SULFOUT(N)

 401    CONTINUE

      end select

#endif
        TR_LEF=1.D0
!       CALL GET_COND_FACTOR(L,N,WMXTR,TPOLD(L),TPOLD(L-1),LHX,FPLUME
!     *       ,FQCOND,FQCONDT,.true.,TRCOND,TM,THLAW,TR_LEF,PL(L),ntix
!     *       ,CLDSAVT)
       CALL GET_COND_FACTOR(L,N,WMXTR,TPOLD(L),TPOLD(L-1),LHX,0.
     *       ,FQCOND,FQCONDT,.true.,TRCOND,TM,THLAW,TR_LEF,PL(L),ntix
     *       ,FPLUME)
        TRCOND(N,L) = FQCONDT * TMP(N) + TRCOND(N,L)
#ifdef TRDIAG_WETDEPO
        IF (diag_wetdep == 1)
     &     trcond_mc(l,n)=trcond_mc(l,n)+fqcondt*tmp(n)
#endif
        TMP(N)         = TMP(N)         *(1.-FQCONDT)
        TMOMP(xymoms,N)= TMOMP(xymoms,N)*(1.-FQCONDT)
      END DO
#endif
      TAUMCL(L)=TAUMCL(L)+DQSUM*FMC1
      CDHEAT(L)=SLH*COND(L)
      CDHSUM=CDHSUM+CDHEAT(L)
      IF(ETADN.GT.1d-10) CDHDRT=CDHDRT+SLH*COND(L)
C****
C**** UPDATE ALL QUANTITIES CARRIED BY THE PLUME
C****
C     MCCONT=MCCONT+1
C     IF(MCCONT.EQ.1) MC1=.TRUE.
C     IF(MC1.AND.PLE(LMIN)-PLE(L+2).GE.450.) SVLATL(L)=LHX
      SVLATL(L)=VLAT(L)
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
      IF (LMAX.LT.LM) GO TO 220   ! CHECK FOR NEXT POSSIBLE LMAX
C**** UPDATE CHANGES CARRIED BY THE PLUME IN THE TOP CLOUD LAYER
  340 IF(LMIN.EQ.LMAX) GO TO 600
      IF(TPSAV(LMAX).GE.TF) LFRZ=LMAX
      DM(LMAX)=DM(LMAX)+MPMAX
      DSM(LMAX)=DSM(LMAX)+SMPMAX
      DSMOM(xymoms,LMAX)=DSMOM(xymoms,LMAX) + SMOMPMAX(xymoms)
      DQM(LMAX)=DQM(LMAX)+QMPMAX
      DQMOM(xymoms,LMAX)=DQMOM(xymoms,LMAX) + QMOMPMAX(xymoms)
#ifdef TRACERS_ON
C**** Add plume tracers at LMAX
      DTM(LMAX,1:NTX) = DTM(LMAX,1:NTX) + TMPMAX(1:NTX)
      DTMOM(xymoms,LMAX,1:NTX) =
     *   DTMOM(xymoms,LMAX,1:NTX) + TMOMPMAX(xymoms,1:NTX)
#endif
      CCM(LMAX)=0.
      DO K=1,KMAX
         DUM(K,LMAX)=DUM(K,LMAX)+UMP(K)
         DVM(K,LMAX)=DVM(K,LMAX)+VMP(K)
      ENDDO
      CDHM=0.
      IF(MINLVL.GT.LMIN) MINLVL=LMIN
      IF(MAXLVL.LT.LMAX) MAXLVL=LMAX
      IF(LMCMIN.EQ.0) LMCMIN=LMIN
      IF(LMCMAX.LT.MAXLVL) LMCMAX=MAXLVL
C****
C**** PROCESS OF DOWNDRAFTING
C****
      LDMIN=LMIN
      EDRAFT=0.
      IF(ETADN.GT.1d-10) THEN
      CONSUM=0.
      DO 347 L=LMIN,LDRAFT
  347 CONSUM=CONSUM+COND(L)
      TNX=SMDN*PLK(LMIN)/DDRAFT
      QNX=QMDN/DDRAFT
      LHX=LHE
      IF(TPSAV(LMIN).LT.TF) LHX=LHS
      SLH=LHX*BYSHA
      QSATC=QSAT(TNX,LHX,PL(LMIN))
      DQ=(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
      DQRAT=DQ*DDRAFT/(CONSUM+teeny)
      DO 346 L=LDRAFT,1,-1
      LHX=LHE
      IF (L.GE.LMIN) THEN  ! avoids reference to uninitiallised value
        IF(TPSAV(L).LT.TF) LHX=LHS
      END IF
      SLH=LHX*BYSHA
      DQEVP=DQRAT*COND(L)
      IF(DQEVP.GT.COND(L)) DQEVP=COND(L)
      IF (L.LT.LMIN) DQEVP=0.
      FSEVP = 0
      IF (ABS(PLK(L)*SMDN).gt.teeny) FSEVP = SLH*DQEVP/(PLK(L)*SMDN)
      SMDN=SMDN-SLH*DQEVP/PLK(L)
      SMOMDN(xymoms)=SMOMDN(xymoms)*(1.-FSEVP)
      FQEVP = 0
      IF (COND(L).gt.0.) FQEVP = DQEVP/COND(L)
      QMDN=QMDN+DQEVP
      COND(L)=COND(L)-DQEVP
      TAUMCL(L)=TAUMCL(L)-DQEVP*FMC1
      CDHEAT(L)=CDHEAT(L)-DQEVP*SLH
      EVPSUM=EVPSUM+DQEVP*SLH
#ifdef TRACERS_WATER
C**** RE-EVAPORATION OF TRACERS IN DOWNDRAFTS
C**** (If 100% evaporation, allow all tracers to evaporate completely.)
      DO N=1,NTX
        IF(FQEVP.eq.1.) THEN                 ! total evaporation
          TMDN(N)     = TMDN(N) + TRCOND(N,L)
#ifdef TRDIAG_WETDEPO
          IF (diag_wetdep == 1)
     &         trdvap_mc(l,n)=trdvap_mc(l,n)+trcond(n,l)
#endif
          TRCOND(N,L) = 0.d0
        ELSE ! otherwise, tracers evaporate dependent on type of tracer
          CALL GET_EVAP_FACTOR(N,TNX,LHX,.FALSE.,1d0,FQEVP,FQEVPT,ntix)
          TMDN(N)     = TMDN(N)     + FQEVPT * TRCOND(N,L)
          TRCOND(N,L) = TRCOND(N,L) - FQEVPT * TRCOND(N,L)
#ifdef TRDIAG_WETDEPO
          IF (diag_wetdep == 1)
     &         trdvap_mc(l,n)=trdvap_mc(l,n)+fqevpt*trcond(n,l)
#endif
        END IF
      END DO
#endif
C**** ENTRAINMENT OF DOWNDRAFTS
      IF(L.LT.LDRAFT.AND.L.GT.1) THEN
        DDRUP=DDRAFT
        DDRAFT=DDRAFT*(1.+ETAL(L))
        IF(DDRUP.GT.DDRAFT) DDRAFT=DDRUP
        IF(DDRAFT.GT..95d0*(AIRM(L-1)+DMR(L-1)))
     *    DDRAFT=.95d0*(AIRM(L-1)+DMR(L-1))
        EDRAFT=DDRAFT-DDRUP
        IF (EDRAFT.gt.0) THEN  ! usual case, entrainment into downdraft
          FENTRA=EDRAFT*BYAM(L)
          SENV=SM(L)/AIRM(L)
          QENV=QM(L)/AIRM(L)
          SMDN=SMDN+EDRAFT*SENV
          QMDN=QMDN+EDRAFT*QENV
          SMOMDN(xymoms)= SMOMDN(xymoms)+ SMOM(xymoms,L)*FENTRA
          QMOMDN(xymoms)= QMOMDN(xymoms)+ QMOM(xymoms,L)*FENTRA
          DSMR(L)=DSMR(L)-EDRAFT*SENV
          DSMOMR(:,L)=DSMOMR(:,L)-SMOM(:,L)*FENTRA
          DQMR(L)=DQMR(L)-EDRAFT*QENV
          DQMOMR(:,L)=DQMOMR(:,L)-QMOM(:,L)*FENTRA
          DMR(L)=DMR(L)-EDRAFT
#ifdef TRACERS_ON
          Tenv(1:NTX)=tm(l,1:NTX)/airm(l)
          TMDN(1:NTX)=TMDN(1:NTX)+EDRAFT*Tenv(1:NTX)
          TMOMDN(xymoms,1:NTX)= TMOMDN(xymoms,1:NTX)+ TMOM(xymoms,L
     *         ,1:NTX)*FENTRA
          DTMR(L,1:NTX)=DTMR(L,1:NTX)-EDRAFT*TENV(1:NTX)
          DTMOMR(:,L,1:NTX)=DTMOMR(:,L,1:NTX)-TMOM(:,L,1:NTX)*FENTRA
#endif
        ELSE  ! occasionally detrain into environment if ddraft too big
          FENTRA=EDRAFT/DDRUP  ! < 0
          DSM(L)=DSM(L)-FENTRA*SMDN
          DSMOM(xymoms,L)=DSMOM(xymoms,L)-SMOMDN(xymoms)*FENTRA
          DQM(L)=DQM(L)-FENTRA*QMDN
          DQMOM(xymoms,L)=DQMOM(xymoms,L)-QMOMDN(xymoms)*FENTRA
          SMDN=SMDN*(1+FENTRA)
          QMDN=QMDN*(1+FENTRA)
          SMOMDN(xymoms)= SMOMDN(xymoms)*(1+FENTRA)
          QMOMDN(xymoms)= QMOMDN(xymoms)*(1+FENTRA)
          DM(L)=DM(L)-EDRAFT
#ifdef TRACERS_ON
          DTM(L,1:NTX)=DTM(L,1:NTX)-FENTRA*TMDN(1:NTX)
          DTMOM(xymoms,L,1:NTX)=DTMOM(xymoms,L,1:NTX)-
     *         TMOMDN(xymoms,1:NTX)*FENTRA
          TMDN(1:NTX)=TMDN(1:NTX)*(1.+FENTRA)
          TMOMDN(xymoms,1:NTX)= TMOMDN(xymoms,1:NTX)*(1.+FENTRA)
#endif
        END IF
      ENDIF
      IF(L.GT.1) DDM(L-1)=DDRAFT
      LDMIN=L
C**** ALLOW FOR DOWNDRAFT TO DROP BELOW LMIN, IF IT'S NEGATIVE BUOYANT
      IF (L.LE.LMIN.AND.L.GT.1) THEN
        SMIX=SMDN/DDRAFT
        IF (SMIX.GE.(SM1(L-1)/AIRM(L-1))) GO TO 345
      ENDIF
  346 CONTINUE
  345 DSM(LDMIN)=DSM(LDMIN)+SMDN
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
      ENDDO
      DM(LDMIN)=DM(LDMIN)+DDRAFT
      END IF
C****
C**** SUBSIDENCE AND MIXING
C****
C**** Calculate vertical mass fluxes (Note CM for subsidence is defined
C**** in opposite sense than normal (positive is down))
      DO L=0,LDMIN-1
        CM(L) = 0.
      END DO
      DO L=LDMIN,LMAX
        CM(L) = CM(L-1) - DM(L) - DMR(L)
        SMT(L)=SM(L)    ! Save profiles for diagnostics
        QMT(L)=QM(L)
      END DO
      DO L=LMAX+1,LM
        CM(L) = 0.
      END DO
C**** simple upwind scheme for momentum
      ALPHA=0.
      DO 380 L=LDMIN,LMAX
      CLDM=CCM(L)
      IF(L.LT.LDRAFT.AND.ETADN.GT.1d-10) CLDM=CCM(L)-DDM(L)
      IF(MC1) VSUBL(L)=100.*CLDM*RGAS*TL(L)/(PL(L)*GRAV*DTsrc)
      BETA=CLDM*BYAM(L+1)
      IF(CLDM.LT.0.) BETA=CLDM*BYAM(L)
      BETAU=BETA
      ALPHAU=ALPHA
      IF(BETA.LT.0.) BETAU=0.
      IF(ALPHA.LT.0.) ALPHAU=0.
      DO K=1,KMAX
       UM(K,L)=
     *    UM(K,L)+RA(K)*(-ALPHAU*UM(K,L)+BETAU*UM(K,L+1)+DUM(K,L))
       VM(K,L)=
     *    VM(K,L)+RA(K)*(-ALPHAU*VM(K,L)+BETAU*VM(K,L+1)+DVM(K,L))
      ENDDO
  380 ALPHA=BETA
C**** Subsidence uses Quadratic Upstream Scheme for QM and SM
      ML(LDMIN:LMAX) = AIRM(LDMIN:LMAX) +   DMR(LDMIN:LMAX)
      SM(LDMIN:LMAX) =   SM(LDMIN:LMAX) +  DSMR(LDMIN:LMAX)
      SMOM(:,LDMIN:LMAX) =  SMOM(:,LDMIN:LMAX) + DSMOMR(:,LDMIN:LMAX)
C****
      nsub = lmax-ldmin+1
      cmneg(ldmin:lmax)=-cm(ldmin:lmax)
      cmneg(lmax)=0. ! avoid roundoff error (esp. for qlimit)
      call adv1d(sm(ldmin),smom(1,ldmin), f(ldmin),fmom(1,ldmin),
     &     ml(ldmin),cmneg(ldmin), nsub,.false.,1, zdir,ierrt,lerrt)
      SM(LDMIN:LMAX) =   SM(LDMIN:LMAX) +   DSM(LDMIN:LMAX)
      SMOM(:,LDMIN:LMAX) =  SMOM(:,LDMIN:LMAX) +  DSMOM(:,LDMIN:LMAX)
      ierr=max(ierrt,ierr) ; lerr=max(lerrt+ldmin-1,lerr)
C****
      ML(LDMIN:LMAX) = AIRM(LDMIN:LMAX) +   DMR(LDMIN:LMAX)
      QM(LDMIN:LMAX) =   QM(LDMIN:LMAX) +  DQMR(LDMIN:LMAX)
      QMOM(:,LDMIN:LMAX) =  QMOM(:,LDMIN:LMAX) + DQMOMR(:,LDMIN:LMAX)
      call adv1d(qm(ldmin),qmom(1,ldmin), f(ldmin),fmom(1,ldmin),
     &     ml(ldmin),cmneg(ldmin), nsub,.true.,1, zdir,ierrt,lerrt)
      QM(LDMIN:LMAX) =   QM(LDMIN:LMAX) +   DQM(LDMIN:LMAX)
      QMOM(:,LDMIN:LMAX) =  QMOM(:,LDMIN:LMAX) +  DQMOM(:,LDMIN:LMAX)
      ierr=max(ierrt,ierr) ; lerr=max(lerrt+ldmin-1,lerr)
C**** diagnostics
      DO L=LDMIN,LMAX
        FCDH=0.
        IF(L.EQ.LMAX) FCDH=CDHSUM-(CDHSUM-CDHDRT)*.5*ETADN+CDHM
        FCDH1=0.
        IF(L.EQ.LDMIN) FCDH1=(CDHSUM-CDHDRT)*.5*ETADN-EVPSUM
        MCFLX(L)=MCFLX(L)+CCM(L)*FMC1
        DGDSM(L)=DGDSM(L)+(PLK(L)*(SM(L)-SMT(L))-FCDH-FCDH1)*FMC1
        IF(PLE(LMAX+1).GT.700.d0) DGSHLW(L)=DGSHLW(L)+
     *    (PLK(L)*(SM(L)-SMT(L))-FCDH-FCDH1)*FMC1
        IF(PLE(LMIN)-PLE(LMAX+1).GE.450.d0) DGDEEP(L)=DGDEEP(L)+
     *    (PLK(L)*(SM(L)-SMT(L))-FCDH-FCDH1)*FMC1
        DTOTW(L)=DTOTW(L)+SLHE*(QM(L)-QMT(L)+COND(L))*FMC1
        DGDQM(L)=DGDQM(L)+SLHE*(QM(L)-QMT(L))*FMC1
        DDMFLX(L)=DDMFLX(L)+DDM(L)*FMC1
        IF(QM(L).LT.0.d0) WRITE(6,*) ' Q neg: it,i,j,l,q',
     *   itime,i_debug,j_debug,l,qm(l)
      END DO
#ifdef TRACERS_ON
C**** Subsidence of tracers by Quadratic Upstream Scheme
      DO N=1,NTX
        if (debug.and.n.eq.1) print*,"cld0",i_debug,ldmin,lmax
     *       ,DTMR(LDMIN:LMAX,N),DTM(LDMIN:LMAX,N)
        if (debug.and.n.eq.1) print*,"cld1",TM(LDMIN:LMAX,N)
      ML(LDMIN:LMAX) =  AIRM(LDMIN:LMAX) +    DMR(LDMIN:LMAX)
      TM(LDMIN:LMAX,N) =  TM(LDMIN:LMAX,N) + DTMR(LDMIN:LMAX,N)
      TMOM(:,LDMIN:LMAX,N) = TMOM(:,LDMIN:LMAX,N)+DTMOMR(:,LDMIN:LMAX,N)
      call adv1d(tm(ldmin,n),tmom(1,ldmin,n), f(ldmin),fmom(1,ldmin),
     &     ml(ldmin),cmneg(ldmin), nsub,t_qlimit(n),1, zdir,ierrt,lerrt)
      TM(LDMIN:LMAX,N) = TM(LDMIN:LMAX,N) +   DTM(LDMIN:LMAX,N)
        if (debug .and.n.eq.1) print*,"cld2",TM(LDMIN:LMAX,N)
      TMOM(:,LDMIN:LMAX,N) = TMOM(:,LDMIN:LMAX,N) +DTMOM(:,LDMIN:LMAX,N)
      ierr=max(ierrt,ierr) ; lerr=max(lerrt+ldmin-1,lerr)
      END DO
#endif
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
      IF(PLE(LMIN)-PLE(LMAX+1).GE.450.) THEN
        DO L=LMAX,LMIN,-1
          IF(COND(L).LT.CONDP(L)) CONDP(L)=COND(L)
          FCLW=0.
          IF (COND(L).GT.0) FCLW=(COND(L)-CONDP(L))/COND(L)
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
C**** REEVAPORATION AND PRECIPITATION
C****
      PRCP=COND(LMAX)
      PRHEAT=CDHEAT(LMAX)
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
         DPHASE(LMAX)=DPHASE(LMAX)+(CDHSUM-(CDHSUM-CDHDRT)*.5*ETADN+
     *                CDHM)*FMC1
         IF(PLE(LMAX+1).GT.700.d0) DPHASHLW(LMAX)=DPHASHLW(LMAX)+
     *     (CDHSUM-(CDHSUM-CDHDRT)*.5*ETADN+CDHM)*FMC1
         IF(PLE(LMIN)-PLE(LMAX+1).GE.450.d0) DPHADEEP(LMAX)=
     *     DPHADEEP(LMAX)+(CDHSUM-(CDHSUM-CDHDRT)*.5*ETADN+CDHM)*FMC1
      DO 540 L=LMAX-1,1,-1
      IF(PRCP.LE.0.) GO TO 530
      FCLOUD=CCMUL*CCM(L)*BYAM(L+1)
      IF(PLE(LMIN)-PLE(L+2).GE.450.) FCLOUD=CCMUL1*CCM(L)*BYAM(L+1)
      IF(L.LT.LMIN) FCLOUD=CCMUL*CCM(LMIN)*BYAM(LMIN+1)
      IF(PLE(LMIN)-PLE(LMAX+1).LT.450.) THEN
        IF(L.EQ.LMAX-1) FCLOUD=CCMUL2*CCM(L)*BYAM(L+1)
        IF(L.LT.LMIN) FCLOUD=0.
      ENDIF
      IF(FCLOUD.GT.1.) FCLOUD=1.
      FEVAP=.5*CCM(L)*BYAM(L+1)
      IF(L.LT.LMIN) FEVAP=.5*CCM(LMIN)*BYAM(LMIN+1)
      IF(FEVAP.GT..5) FEVAP=.5
      CLDMCL(L+1)=CLDMCL(L+1)+FCLOUD*FMC1
      CLDREF=CLDMCL(L+1)
      IF(PLE(LMAX+1).GT.700..AND.CLDREF.GT.CLDSLWIJ)
     *  CLDSLWIJ=CLDREF
      IF(PLE(LMIN)-PLE(LMAX+1).GE.450..AND.CLDREF.GT.CLDDEPIJ)
     *  CLDDEPIJ=CLDREF
C**** FORWARD STEP COMPUTES HUMIDITY CHANGE BY RECONDENSATION
C**** Q = Q + F(TOLD,PRHEAT,QOLD+EVAP)
      PRECNVL(L+1)=PRECNVL(L+1)+PRCP*BYGRAV
      MCLOUD=0.
      IF(L.LE.LMIN) MCLOUD=2.*FEVAP*AIRM(L)
      TOLD=SMOLD(L)*PLK(L)*BYAM(L)
      TOLD1=SMOLD(L+1)*PLK(L+1)*BYAM(L+1)
C**** decide whether to melt snow/ice
      HEAT1=0.
      IF(L.EQ.LFRZ.OR.(L.LE.LMIN.AND.TOLD.GE.TF.AND.TOLD1.LT.TF.AND.L.GT
     *     .LFRZ)) HEAT1=LHM*PRCP*BYSHA
C**** and deal with possible inversions and re-freezing of rain
      IF (TOLD.LT.TF.AND.TOLD1.GT.TF) HEAT1=-LHM*PRCP*BYSHA
      TNX=TOLD
      QNX=QMOLD(L)*BYAM(L)
      LHX=LHE
      IF(TNX.LT.TF) LHX=LHS
      IF(L.GT.LMIN) THEN    ! needed to avoid unintiallised value
        IF (TPSAV(L).GE.TF) LHX=LHE
      END IF
      SLH=LHX*BYSHA
      DQSUM=0.
      DO 510 N=1,3
      QSATC=QSAT(TNX,LHX,PL(L))
      DQ=(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
      TNX=TNX-SLH*DQ
      QNX=QNX+DQ
  510 DQSUM=DQSUM+DQ*MCLOUD
      IF(DQSUM.LT.0.) DQSUM=0.
      IF(DQSUM.GT.PRCP) DQSUM=PRCP
      FPRCP=DQSUM/PRCP
      PRCP=PRCP-DQSUM
C**** UPDATE TEMPERATURE AND HUMIDITY DUE TO NET REEVAPORATION IN CLOUDS
      FSSUM = 0
      IF (ABS(PLK(L)*SM(L)).gt.teeny) FSSUM = (SLH*DQSUM+HEAT1)/
     *     (PLK(L)*SM(L))
      SM(L)=SM(L)-(SLH*DQSUM+HEAT1)/PLK(L)
      SMOM(:,L) =  SMOM(:,L)*(1.-FSSUM)
      QM(L)=QM(L)+DQSUM
#ifdef TRACERS_WATER
C**** Tracer net re-evaporation
C**** (If 100% evaporation, allow all tracers to evaporate completely.)
      BELOW_CLOUD = L.lt.LMIN
      IF(FPRCP.eq.1.) THEN      !total evaporation
        DO N=1,NTX
          TM(L,N)   = TM(L,N)  + TRPRCP(N)
          if (debug .and.n.eq.1) print*,"cld2",L,TM(L,N),TRPRCP(N),2
     *         *FEVAP
#ifdef TRDIAG_WETDEPO
          IF (diag_wetdep == 1)
     &         trnvap_mc(l,n)=trnvap_mc(l,n)+trprcp(n)
#endif
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
        DO N=1,NTX
          CALL GET_EVAP_FACTOR(N,TNX,LHX,BELOW_CLOUD,HEFF,FPRCP,FPRCPT
     *         ,ntix)
          TM(L,N) = TM(L,N)     + FPRCPT*TRPRCP(N)
          if (debug .and.n.eq.1) print*,"cld3",L,TM(L,N),FPRCP
     *         ,FPRCPT,TRPRCP(N)
          TRPRCP(N) = TRPRCP(N) - FPRCPT*TRPRCP(N)
#ifdef TRDIAG_WETDEPO
          IF (diag_wetdep == 1)
     &         trnvap_mc(l,n)=trnvap_mc(l,n)+fprcpt*trprcp(n)
#endif
        END DO
      END IF
#ifndef NO_WASHOUT_IN_CLOUDS
      IF (.NOT. below_cloud .AND. prcp > teeny) THEN
c**** Washout of tracers in cloud
        wmxtr=prcp*byam(l)
        precip_mm=prcp*100.*bygrav
        b_beta_DT=fplume
        DO N=1,NTX
          CALL get_wash_factor(n,b_beta_dt,precip_mm,fwasht,tnx,lhx,
     &         wmxtr,fplume,l,tm,trprcp,thlaw,pl(l),ntix)
          trprcp(n)=fwasht*tm(l,n)+trprcp(n)+thlaw
#ifdef TRDIAG_WETDEPO
          IF (diag_wetdep == 1)
     &         trwash_mc(l,n)=trwash_mc(l,n)+fwasht*tm(l,n)+thlaw
#endif
          IF (tm(l,n) > teeny) THEN
            tmfac=thlaw/tm(l,n)
          ELSE
            tmfac=0.
          ENDIF
          tm(l,n)=tm(l,n)*(1.-fwasht)-thlaw
          tmom(xymoms,l,n)=tmom(xymoms,l,n)*(1.-fwasht-tmfac)
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

        CALL GET_SULFATE(L,TNX,FPLUME,WA_VOL,WMXTR,SULFIN,
     *       SULFINC,SULFOUT,TR_LEFT,TM,TRPRCP,AIRM,LHX,
     *       DT_SULF_MC(1,L),CLDSAVT)

#endif
        DO N=1,NTX
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      select case (trname(ntix(n)))
      case('SO2','SO4','H2O2_s','H2O2')
      if (trname(ntix(n)).eq."H2O2" .and. coupled_chem.eq.0) goto 402
      if (trname(ntix(n)).eq."H2O2_s" .and. coupled_chem.eq.1) goto 402

          TRPRCP(N)=TRPRCP(N)*(1.+SULFINC(N))
          TM(L,N)=TM(L,N)*(1.+SULFIN(N))
          TMOM(xymoms,L,N)=TMOM(xymoms,L,N) *(1.+SULFIN(N))
          TRCOND(N,L) = TRCOND(N,L)+SULFOUT(N)

 402      CONTINUE

      end select

#endif
cdmk Here I took out GET_COND, since we are below cloud.
cdmk GET_WASH now has gas dissolution, extra arguments
          CALL GET_WASH_FACTOR(N,b_beta_DT,precip_mm,FWASHT
     *         ,TNX,LHX,WMXTR,FPLUME,L,TM,TRPRCP,THLAW,pl(l),ntix)
          TRPRCP(N) = FWASHT*TM(L,N)+TRPRCP(N)+THLAW
#ifdef TRDIAG_WETDEPO
          IF (diag_wetdep == 1)
     &         trwash_mc(l,n)=trwash_mc(l,n)+fwasht*tm(l,n)+thlaw
#endif
          IF (TM(L,N).GT.teeny) THEN
            TMFAC=THLAW/TM(L,N)
          ELSE
            TMFAC=0.
          ENDIF
          TM(L,N)=TM(L,N)*(1.-FWASHT)-THLAW
          if (debug .and.n.eq.1) print*,"cld4",L,TM(L,N),FWASHT
     *         ,THLAW
          TMOM(xymoms,L,N)=TMOM(xymoms,L,N) *
     &                (1.-FWASHT-TMFAC)
        END DO
      END IF
#endif
         FCDH1=0.
         IF(L.EQ.LDMIN) FCDH1=(CDHSUM-CDHDRT)*.5*ETADN-EVPSUM
         DPHASE(L)=DPHASE(L)-(SLH*DQSUM-FCDH1+HEAT1)*FMC1
         IF(PLE(LMAX+1).GT.700.d0) DPHASHLW(L)=DPHASHLW(L)-
     *     (SLH*DQSUM-FCDH1+HEAT1)*FMC1
         IF(PLE(LMIN)-PLE(LMAX+1).GE.450.d0) DPHADEEP(L)=DPHADEEP(L)-
     *     (SLH*DQSUM-FCDH1+HEAT1)*FMC1
         DQCOND(L)=DQCOND(L)-SLH*DQSUM*FMC1
C**** ADD PRECIPITATION AND LATENT HEAT BELOW
  530 PRHEAT=CDHEAT(L)+SLH*PRCP
      PRCP=PRCP+COND(L)
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
          CALL ISOEQUIL(NTIX(N),TNX,.TRUE.,QM(L),PRCP,TM(L,N),TRPRCP(N)
     *         ,0.5d0)
        END DO
      END IF
#endif
#endif
  540 CONTINUE
C****
      IF(PRCP.GT.0.) CLDMCL(1)=CLDMCL(1)+CCM(LMIN)*BYAM(LMIN+1)*FMC1
      PRCPMC=PRCPMC+PRCP*FMC1
#ifdef TRACERS_WATER
      TRPRMC(1:NTX) = TRPRMC(1:NTX) + TRPRCP(1:NTX)*FMC1
#endif
      IF(LMCMIN.GT.LDMIN) LMCMIN=LDMIN
C****
C**** END OF LOOP OVER CLOUD TYPES
C****
      MC1=.FALSE.                        !!!
  570 CONTINUE
C****
C**** END OF OUTER LOOP OVER CLOUD BASE
C****
  600 CONTINUE
      IF(LMCMIN.GT.0) THEN
C**** set fssl array
        DO L=1,LM
          FSSL(L)=1.-FMC1
C         FMCL(L)=FMC1                 ! may be generalized
        END DO
C**** ADJUSTMENT TO CONSERVE CP*T
        SUMAJ=0.
        SUMDP=0.
        DO L=LMCMIN,LMCMAX
          SUMDP=SUMDP+AIRM(L)*FMC1
          SUMAJ=SUMAJ+DGDSM(L)
        END DO
        DO L=LMCMIN,LMCMAX
          DGDSM(L)=DGDSM(L)-SUMAJ*AIRM(L)*FMC1/SUMDP
          SM(L)=SM(L)-SUMAJ*AIRM(L)/(SUMDP*PLK(L))
        END DO
C**** LOAD MASS EXCHANGE ARRAY FOR GWDRAG
        AIRXL = 0.
        DO L=LMCMIN,LMCMAX
          AIRXL = AIRXL+MCFLX(L)
        END DO
      END IF
C**** CALCULATE OPTICAL THICKNESS
      WCONST=WMU*(1.-PEARTH)+WMUL*PEARTH
      WMSUM=0.
#ifdef CLD_AER_CDNC
      WMCLWP=0.
      WMCTWP=0.
      ACDNWM=0.
      ACDNIM=0.
      AREWM=0.
      AREIM=0.
      ALWWM=0.
      ALWIM=0.
      NMCW=0
      NMCI=0
#endif
      DO L=1,LMCMAX
         TL(L)=(SM(L)*BYAM(L))*PLK(L)
         TEMWM=(TAUMCL(L)-SVWMXL(L)*AIRM(L))*1.d2*BYGRAV
         IF(TL(L).GE.TF) WMSUM=WMSUM+TEMWM
#ifdef CLD_AER_CDNC
         WMCTWP=WMCTWP+TEMWM
         IF(TL(L).GE.TF) WMCLWP=WMCLWP+TEMWM
#endif
         IF(CLDMCL(L).GT.0.) THEN
               TAUMCL(L)=AIRM(L)*COETAU
            IF(L.EQ.LMCMAX .AND. PLE(LMCMIN)-PLE(LMCMAX+1).LT.450.)
     *         TAUMCL(L)=AIRM(L)*.02d0
            IF(L.LE.LMCMIN .AND. PLE(LMCMIN)-PLE(LMCMAX+1).GE.450.)
     *         TAUMCL(L)=AIRM(L)*.02d0
         END IF
         IF(SVLATL(L).EQ.0.) THEN
            SVLATL(L)=LHE
            IF ( (TPSAV(L).gt.0. .and. TPSAV(L).LT.TF) .or.
     *           (TPSAV(L).eq.0. .and. TL(L).lt.TF) ) SVLATL(L)=LHS
         ENDIF
         IF(SVWMXL(L).GT.0.) THEN
            FCLD=CLDMCL(L)+1.E-20
            TEM=1.d5*SVWMXL(L)*AIRM(L)*BYGRAV
            WTEM=1.d5*SVWMXL(L)*PL(L)/(FCLD*TL(L)*RGAS)
            IF(SVLATL(L).EQ.LHE.AND.SVWMXL(L)/FCLD.GE.WCONST*1.d-3)
     *           WTEM=1d2*WCONST*PL(L)/(TL(L)*RGAS)
            IF(SVLATL(L).EQ.LHS.AND.SVWMXL(L)/FCLD.GE.WMUI*1.d-3)
     *           WTEM=1d2*WMUI*PL(L)/(TL(L)*RGAS)
            IF(WTEM.LT.1.d-10) WTEM=1.d-10
C**   Set CDNC for moist conv. clds (const at present)
              MNdO = 59.68d0/(RWCLDOX**3)
              MNdL = 174.d0
              MNdI = 0.06417127d0
#ifdef CLD_AER_CDNC
             MNdO=MCDNO1
             MNdL=MCDNL1
#endif
              MCDNCW=MNdO*(1.-PEARTH)+MNdL*PEARTH
              MCDNCI=MNdI
            IF(SVLATL(L).EQ.LHE)  THEN
!              RCLD=(RWCLDOX*10.*(1.-PEARTH)+7.0*PEARTH)*(WTEM*4.)**BY3
               RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*MCDNCW))**BY3
             ELSE
!              RCLD=25.0*(WTEM/4.2d-3)**BY3 * (1.+pl(l)*xRICld)
               RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*MCDNCI))**BY3
     *              *(1.+pl(l)*xRICld)
            END IF
            RCLDE=RCLD/BYBR   !  effective droplet radius in anvil
#ifdef CLD_AER_CDNC
C** Using the Liu and Daum paramet, Nature, 2002, Oct 10, Vol 419
C** for spectral dispersion effects on droplet size distribution
C** central value of 0.003 for alfa:Rotstayn&Daum, J.Clim,2003,16,21, Nov 2003.
      Repsi=1.d0 - 0.7d0*exp(-0.003d0*MCDNCW)
      Repsis=Repsi*Repsi
      Rbeta=(((1.d0+2.d0*Repsis)**0.667d0))/((1.d0+Repsis)**by3)
      RCLDE=RCLD*Rbeta
c     write(6,*)"RCLD",RCLDE,RCLD,Rbeta,WTEM,L,MCDNCW
#endif

            CSIZEL(L)=RCLDE   !  effective droplet radius in anvil
#ifdef CLD_AER_CDNC
            if (FCLD.gt.1.d-5.and.SVLATL(L).eq.LHE) then
              ACDNWM(L)= MCDNCW
              AREWM(L) = RCLDE
              ALWWM(L)= WTEM   ! cld water density in g m-3
              NMCW  = NMCW+1
            elseif(FCLD.gt.1.d-5.and.SVLATL(L).eq.LHS) then
              ACDNIM(L)= MCDNCI
              AREIM(L) = RCLDE
              ALWIM(L) = WTEM
              NMCI  = NMCI+1
            ENDIF
#endif
            TAUMCL(L)=1.5*TEM/(FCLD*RCLDE+1.E-20)
            IF(TAUMCL(L).GT.100.) TAUMCL(L)=100.
         END IF
         IF(TAUMCL(L).LT.0..and.CLDMCL(L).LE.0.) TAUMCL(L)=0.
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
         IF (diag_wetdep == 1) THEN
           trcond_mc(l,1:ntx)=trcond_mc(l,1:ntx)*fmc1
           trdvap_mc(l,1:ntx)=trdvap_mc(l,1:ntx)*fmc1
           trflcw_mc(l,1:ntx)=trflcw_mc(l,1:ntx)*fmc1
           trprcp_mc(l,1:ntx)=trprcp_mc(l,1:ntx)*fmc1
           trnvap_mc(l,1:ntx)=trnvap_mc(l,1:ntx)*fmc1
           trwash_mc(l,1:ntx)=trwash_mc(l,1:ntx)*fmc1
         END IF
#endif
      END DO

      RETURN
      END SUBROUTINE MSTCNV

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
      REAL*8, PARAMETER :: CM00=1.d-4, AIRM0=100.d0, GbyAIRM0=GRAV/AIRM0
      REAL*8, PARAMETER :: HEFOLD=500.,COEFM=10.,COEFT=2.5
      REAL*8, PARAMETER :: COESIG=1d-3,COEEC=1000.
      INTEGER, PARAMETER :: ERP=2
      REAL*8, DIMENSION(KMAX) :: UMO1,UMO2,UMN1,UMN2 !@var dummy variables
      REAL*8, DIMENSION(KMAX) :: VMO1,VMO2,VMN1,VMN2 !@var dummy variables
!@var Miscellaneous vertical arrays
      REAL*8, DIMENSION(LM) ::
     *     QSATL,RHF,ATH,SQ,ER,QHEAT,
     *     CAREA,PREP,RH00,EC,WMXM
!@var QSATL saturation water vapor mixing ratio
!@var RHF environmental relative humidity
!@var ATH change in potential temperature
!@var SQ ERMAX dummy variables
!@var ER evaporation of precip
!@var QHEAT change of latent heat
!@var CAREA fraction of clear region
!@var PREP precip conversion rate
!@var RH00 threshold relative humidity
!@var EC cloud evaporation rate
!@var WMXM cloud water mass (mb)
      REAL*8, DIMENSION(LM+1) :: PREBAR,PREICE
!@var PREBAR,PREICE precip entering layer top for total, snow

#ifdef TRACERS_WATER
!@var TRPRBAR tracer precip entering layer top for total (kg)
      REAL*8, DIMENSION(NTM,LM+1) :: TRPRBAR
!@var DTPR change of tracer by precip (kg)
!@var DTER change of tracer by evaporation (kg)
!@var DTQW change of tracer by condensation (kg)
!@var FWTOQ fraction of CLW that goes to water vapour
!@var FPR fraction of CLW that precipitates
!@var FER fraction of precipitate that evaporates
      REAL*8 DTPR,DTER,DTQW,TWMTMP,DTSUM,FWTOQ,FPR,FER
!@var DTPRT tracer-specific change of tracer by precip (kg)
!@var DTERT tracer-specific change of tracer by evaporation (kg)
!@var DTWRT tracer-specific change of tracer by washout (kg)
!@var DTQWT tracer-specific change of tracer by condensation (kg)
!@var FWTOQT tracer-specific fraction of tracer in CLW that evaporates
!@var FQTOWT tracer-specific fraction of gas tracer condensing in CLW
!@var FPRT tracer-specific fraction of tracer in CLW that precipitates
!@var FERT tracer-specific fraction of tracer in precipitate evaporating
      REAL*8 ::DTWRT,DTPRT,DTERT,DTQWT,FWTOQT,FQTOWT,FPRT,FERT,PRLIQ
!@var BELOW_CLOUD logical- is the current level below cloud?
!@var CLOUD_YET logical- in L loop, has there been any cloud so far?
      LOGICAL BELOW_CLOUD,CLOUD_YET
!@var FQCONDT fraction of tracer that condenses
!@var FQEVPT  fraction of tracer that evaporates (in downdrafts)
!@var FPRCPT fraction of tracer that evaporates (in net re-evaporation)
!@var FWASHT  fraction of tracer scavenged by below-cloud precipitation
      REAL*8 :: FQCONDT, FWASHT, FPRCPT, FQEVPT
!@var WMXTR available water mixing ratio for tracer condensation ( )?
!@var b_beta_DT precipitating gridbox fraction from lowest precipitating
!@+   layer. The name was chosen to correspond to Koch et al. p. 23,802.
!@var precip_mm precipitation (mm) from the grid box above for washout
      REAL*8 WMXTR, b_beta_DT, precip_mm
c for tracers in general, added by Koch
      REAL*8 THLAW,THWASH,TR_LEF,TMFAC,TMFAC2,TR_LEFT(ntm),CLDSAVT
!@var TR_LEF limits precurser dissolution following sulfate formation
!@var THLAW Henry's Law determination of amount of tracer dissolution
!@var THWASH Henry's Law for below cloud dissolution
!@var TMFAC,TMFAC2 used to adjust tracer moments
!@var CLDSAVT is present cloud fraction, saved for tracer use
!@var cldprec cloud fraction at lowest precipitating level
      REAL*8 :: cldprec
#ifndef NO_WASHOUT_IN_CLOUDS
!@var tm_temp temporary array for tracer mass after applying other removal
!@+           and re-evaporation processes to calculate washout in clouds
      REAL(8),DIMENSION(Lm,Ntm) :: tm_temp
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
c for sulfur chemistry
!@var WA_VOL Cloud water volume (L). Used by GET_SULFATE.
      REAL*8 WA_VOL
      REAL*8, DIMENSION(NTM) ::SULFOUT,SULFIN,SULFINC
#endif
#endif

      REAL*8 AIRMR,BETA,BMAX
     *     ,CBF,CBFC0,CK,CKIJ,CK1,CK2,CKM,CKR,CM,CM0,CM1,DFX,DQ,DQSDT
     *     ,DQSUM,DQUP,DRHDT,DSE,DSEC,DSEDIF,DWDT,DWDT1,DWM,ECRATE,EXPST
     *     ,FCLD,FMASS,FMIX,FPLUME,FPMAX,FQTOW,FRAT,FUNI
     *     ,FUNIL,FUNIO,HCHANG,HDEP,HPHASE,OLDLAT,OLDLHX,PFR,PMI,PML
     *     ,HPBL,PRATIO,QCONV,QHEATC,QLT1,QLT2,QMN1,QMN2,QMO1,QMO2,QNEW
     *     ,QNEWU,QOLD,QOLDU,QSATC,QSATE,RANDNO,RCLDE,RHI,RHN,RHO,RHT1
     *     ,RHW,SEDGE,SIGK,SLH,SMN1,SMN2,SMO1,SMO2,TEM,TEMP,TEVAP,THT1
     *     ,THT2,TLT1,TNEW,TNEWU,TOLD,TOLDU,TOLDUP,VDEF,WCONST,WMN1,WMN2
     *     ,WMNEW,WMO1,WMO2,WMT1,WMT2,WMX1,WTEM,VVEL,RCLD,FCOND
     *     ,PRATM,SMN12,SMO12
       real*8 SNdO,SNdL,SNdI,SCDNCW,SCDNCI
#ifdef CLD_AER_CDNC
!@auth Menon  - storing var for cloud droplet number
       integer, PARAMETER :: SNTM=17+ntm_soa/2
       real*8 Repsis,Repsi,Rbeta,CDNL1,QAUT,DSU(SNTM),QCRIT
     * ,CDNL0,NEWCDN,OLDCDN,SNd
       real*8 dynvis(LM),DSGL(LM,SNTM),DSS(SNTM),r6,r6c
       real*8 DPP,TEMPR,RHODK,PPRES,PRS        ! for 3 hrly diag
       real*8 D3DL(LM),CWCON(LM)               ! for 3 hrly diag
#endif
#ifdef BLK_2MOM
      integer,PARAMETER         :: mkx=1   ! lm
      real*8,PARAMETER          :: mw0 = 2.094395148947515E-15
      real*8,PARAMETER          :: mi0 = 2.094395148947515E-15
      logical, PARAMETER        :: lSCM=.false.
      logical, PARAMETER        :: wSCM=.false.
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
!@var FCOND dummy variables
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
C**** initialise vertical arrays
      ER=0.
      EC=0.
      PREP=0.
      PREBAR=0.
      LHP=0.
      QHEAT=0.
      CLDSSL=0
      TAUSSL=0
#ifdef TRACERS_WATER
      TRPRSS = 0.
      TRPRBAR = 0.
      BELOW_CLOUD=.false.
      CLOUD_YET=.false.
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      DT_SULF_SS(1:NTM,:)=0.
#endif
#endif
#ifdef CLD_AER_CDNC
       CTEML=0.
       CD3DL=0.
       CL3DL=0.
       CI3DL=0.
       CDN3DL=0.
       CRE3DL=0.
       SMLWP=0.
       DSGL(:,1:SNTM)=0.
#endif
#ifdef BLK_2MOM
       WMICE(:)=0.
c      print *,sname,'WMX, WMICE = ', WMX, WMICE
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
        CAREA(L)=1.-CLDSAVL(L)
C       IF(RH(L).LE.1.) CAREA(L)=DSQRT((1.-RH(L))/(1.-RH00(L)+teeny))
C       IF(RH(L).LE.1.) CAREA(L)=DSQRT((1.-RH(L))/(1.-U00ice+teeny))
C       IF(CAREA(L).GT.1.) CAREA(L)=1.
C       IF(RH(L).GT.1.) CAREA(L)=0.
        IF(WMX(L).LE.0.) CAREA(L)=1.
      END DO
      DQUP=0.
      TOLDUP=TL(LP50)
      PREICE(LP50+1)=0.
      WCONST=WMU*(1.-PEARTH)+WMUL*PEARTH
         SSHR=0.
         DCTEI=0.
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
C**** COMPUTE THE LIMITING AUTOCONVERSION RATE FOR CLOUD WATER CONTENT
      CM0=CM00
      VDEF=VVEL-VSUBL(L)
      IF(VDEF.GT.0.) CM0=CM00*10.**(-VDEF)
c      FCLD=1.-CAREA(L)+1.E-20            !!!
      FCLD=(1.-CAREA(L))*FSSL(L)+teeny
C**** COMPUTE THE PROBABILITY OF ICE FORMATION, FUNI, AND
C**** THE PROBABLITY OF GLACIATION OF SUPER-COOLED WATER, PFR
C**** DETERMINE THE PHASE MOISTURE CONDENSES TO
C**** DETERMINE THE POSSIBILITY OF B-F PROCESS
      BANDF=.FALSE.
      LHX=LHE
      CBF=1. + EXP(-((TL(L)-258.16d0)/10.)**2)

      IF (TL(L).LE.TI) THEN     ! below -40: force ice
        LHX=LHS
      ELSEIF (TL(L).GE.TF) THEN ! above freezing: force water
        LHX=LHE
      ELSE                      ! in between: compute probability
        IF(TL(L).GT.269.16) THEN ! OC/SI/LI clouds: water above -4
          FUNIO=0.
        ELSE
C         FUNIO=1.-EXP(-((TL(L)-269.16d0)/15.)**2)
          FUNIO=1.-EXP(-((TL(L)-269.16d0)/RTEMP)**2)
        END IF
        IF(TL(L).GT.263.16) THEN ! land clouds water: above -10
          FUNIL=0.
        ELSE
          FUNIL=1.-EXP(-((TL(L)-263.16d0)/15.)**2)
        END IF
        FUNI=FUNIO*(1.-PEARTH)+FUNIL*PEARTH
        RANDNO=RNDSSL(1,L)       !  RANDNO=RANDU(XY)
        IF(RANDNO.LT.FUNI) LHX=LHS

C**** special case 1) if ice previously then stay as ice (if T<Tf)
        IF((OLDLHX.EQ.LHS.OR.OLDLAT.EQ.LHS).AND.TL(L).LT.TF) THEN
          IF(LHX.EQ.LHE) BANDF=.TRUE.
          LHX=LHS
        ENDIF

        IF (L.LT.LP50) THEN
C**** Decide whether precip initiates B-F process
          PML=WMX(L)*AIRM(L)*BYGRAV
          PMI=PREICE(L+1)*DTsrc
          RANDNO=RNDSSL(2,L)     !  RANDNO=RANDU(XY)
C**** Calculate probability of ice precip seeding a water cloud
          IF (LHX.EQ.LHE.AND.PMI.gt.0) THEN
            PRATIO=MIN(PMI/(PML+1.E-20),10d0)
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

C**** COMPUTE RELATIVE HUMIDITY
      QSATL(L)=QSAT(TL(L),LHX,PL(L))
      RH1(L)=QL(L)/QSATL(L)
      IF(LHX.EQ.LHS) THEN
        QSATE=QSAT(TL(L),LHE,PL(L))
        RHW=.00536d0*TL(L)-.276d0
        RH1(L)=QL(L)/QSATE
        IF(TL(L).LT.238.16) RH1(L)=QL(L)/(QSATE*RHW)
      END IF
C**** PHASE CHANGE OF CLOUD WATER CONTENT
      HCHANG=0.
      IF(OLDLHX.EQ.LHE.AND.LHX.EQ.LHS) HCHANG= WML(L)*LHM
      IF(OLDLHX.EQ.LHS.AND.LHX.EQ.LHE) HCHANG=-WML(L)*LHM
      IF(OLDLAT.EQ.LHE.AND.LHX.EQ.LHS) HCHANG=HCHANG+SVWMXL(L)*LHM
      IF(OLDLAT.EQ.LHS.AND.LHX.EQ.LHE) HCHANG=HCHANG-SVWMXL(L)*LHM
      SVLHXL(L)=LHX
      TL(L)=TL(L)+HCHANG/(SHA*FSSL(L)+teeny)
      TH(L)=TL(L)/PLK(L)
      ATH(L)=(TH(L)-TTOLDL(L))*BYDTsrc
C**** COMPUTE RH IN THE CLOUD-FREE AREA, RHF
      RHI=QL(L)/QSAT(TL(L),LHS,PL(L))
    ! this formulation is used for consistency with current practice
      RH00(L)=U00wtrX*U00ice
      IF(LHX.EQ.LHS) RH00(L)=U00ice
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
      IF(RH00(L).GT.1.) RH00(L)=1.
      RHF(L)=RH00(L)+(1.-CAREA(L))*(1.-RH00(L))
C**** Set precip phase to be the same as the cloud, unless precip above
C**** is ice and temperatures after ice melt would still be below TFrez
      LHP(L)=LHX
      IF (LHP(L+1).eq.LHS .and.
     *     TL(L).lt.TF+DTsrc*LHM*PREICE(L+1)*GRAV*BYAM(L)*BYSHA/(FSSL(L)
     *     +teeny)) LHP(L)=LHP(L+1)
#ifdef CLD_AER_CDNC
!@auth Menon  saving aerosols mass for CDNC prediction
      DO N=1,SNTM
        DSS(N)=1.d-10
        DSGL(L,N)=1.d-10
      ENDDO
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      DO N=1,NTX
       select case (trname(ntix(n)))
        case('SO4')
        DSGL(L,1)=tm(l,n)  !n=19
        DSS(1) = DSGL(L,1)
        case('seasalt1')
        DSGL(L,2)=tm(l,n)  !n=21
        DSS(2) = DSGL(L,2)
        case('seasalt2')
        DSGL(L,3)=tm(l,n)  !n=22
        DSS(3) = DSGL(L,3)
        case('OCIA')
        DSGL(L,4)=tm(l,n)    !n=27
        DSS(4) = DSGL(L,4)
        case('OCB')
        DSGL(L,5)=tm(l,n)  !n=28
        DSS(5) = DSGL(L,5)
        case('BCIA')
        DSGL(L,6)=tm(l,n)  !n=24
        DSS(6) = DSGL(L,6)
        case('BCB')
        DSGL(L,7)=tm(l,n)  !n=25
        DSS(7) = DSGL(L,7)
        case('OCII')
        DSGL(L,8)=tm(l,n)  !n=26
        DSS(8) = DSGL(L,8)
        case('BCII')
        DSGL(L,9)=tm(l,n)  !n=23
        DSS(9) = DSGL(L,9)
#ifdef TRACERS_DUST
        case('Clay')
        DSGL(L,10)=tm(l,n)  !n=23
        DSS(10) = max(1.d-10,DSGL(L,10))
        case('Silt1')
        DSGL(L,11)=tm(l,n)  !n=23
        DSS(11) = max(1.d-10,DSGL(L,11))
        case('Silt2')
        DSGL(L,12)=tm(l,n)  !n=23
        DSS(12) = max(1.d-10,DSGL(L,12))
        case('Silt3')
        DSGL(L,13)=tm(l,n)  !n=23
        DSS(13) = max(1.d-10,DSGL(L,13))
#endif
#ifdef TRACERS_NITRATE
        case('NO3p')
        DSGL(L,14)=tm(l,n)  !n=23
        DSS(14) = DSGL(L,14)
#endif
#ifdef TRACERS_HETCHEM
C*** Here are dust particles coated with sulfate
       case('SO4_d1')
       DSGL(L,15)=tm(l,n)  !n=20
       DSS(15) = max(1.d-10,DSGL(L,15))
       case('SO4_d2')
       DSGL(L,16)=tm(l,n)  !n=21
       DSS(16) = max(1.d-10,DSGL(L,16))
       case('SO4_d3')
       DSGL(L,17)=tm(l,n)  !n=22
       DSS(17) = max(1.d-10,DSGL(L,17))
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
#endif
#endif
C***Setting constant values of CDNC over land and ocean to get RCLD=f(CDNC,LWC)
      SNdO = 59.68d0/(RWCLDOX**3)
      SNdL = 174.d0
      SNdI = 0.06417127d0
      SCDNCW=SNdO*(1.-PEARTH)+SNdL*PEARTH
      SCDNCI=SNdI
#ifdef CLD_AER_CDNC
      CALL GET_CDNC(L,LHX,WCONST,WMUI,AIRM(L),WMX(L),DXYPIJ,
     *FCLD,CAREA(L),CLDSAVL(L),DSS,PL(L),TL(L),
     *OLDCDL(L),VVEL,SME(L),DSU,CDNL0,CDNL1)
      SNd=CDNL1
cC** Pass old and new cloud droplet number
      NEWCDN=SNd
      OLDCDN=CDNL0
cc    if(L.eq.1) write(6,*)"SM CDNC",NEWCDN,OLDCDN,SNdL,CDNL0
#endif
#ifdef BLK_2MOM
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
          WMICE(L) = WMX(L)
          mcrys=WMICE(L)         ! crys content, [kg water/kg air]
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
c        if(l.eq.1) then
c       if(nactc(l,nm).gt.1.)write(6,*)"Call matrix",naero(mkx,nm)*1.e-6
c     *,nm
c        endif
        enddo
       ldummy=execute_bulk2m_driver('all'
     *           ,ndrop,mdrop,ncrys,mcrys,naero,nmodes,'end')
#else
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
        OLDCDNC=OLDCDN*1.d6  !convert from cm-3 to m-3
        NEWCDNC=NEWCDN*1.d6  !convert from cm-3 to m-3
#ifdef TRACERS_AMP
!#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
C*** Using the AMP_actv interface from MATRIX
        ldummy=execute_bulk2m_driver('matr','drop_nucl',dtB2M,mkx)
#else
        ldummy=execute_bulk2m_driver('gult','drop_nucl',dtB2M,mkx,
     &OLDCDNC,NEWCDNC)
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

         write(iuo,*) 'wmx wmice tl ql pl svlhxl lhx '
         write(iuo,*) wmx
         write(iuo,*) wmice
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
         write(6,*)"We reached -ve con.",ndrop_res(mkx),mdrop_res(mkx),
     * ncrys_res(mkx), mcrys_res(mkx),l
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
       SNdI=ncrys(mkx)*1.0d-6          ! ncrys, [No/m^3]
c      if(SNdI.gt.0.) write(6,*)"ICE CRY",SNdI, SNdI/dtB2M
       if(SNDI.gt.1.d0) SNdI=1.d0      !try to limit to 1000 /l
       SNd=ndrop(mkx)*1.d-6                 ! ndrop, [No/m^3]
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
#endif
#ifdef CLD_AER_CDNC
      SCDNCW=SNd      ! we have already passed the grid box value
      SCDNCI=SNdI
c     if (SCDNCI.le.0.06d0) SCDNCI=0.06417127d0 !set min ice crystal, do we need this, please check
      if (SCDNCW.le.20.d0) SCDNCW=20.d0         !set min CDNC, sensitivity test
      if (SCDNCW.ge.1400.d0) SCDNCW=1400.d0     !set max CDNC, sensitivity test
c     write(6,*)"CDNC LSS",SCDNCW,SNdO,SNdL,L
#endif
C**** COMPUTE THE AUTOCONVERSION RATE OF CLOUD WATER TO PRECIPITATION
      IF(WMX(L).GT.0.) THEN
        RHO=1d5*PL(L)/(RGAS*TL(L))
        TEM=RHO*WMX(L)/(WCONST*FCLD+teeny)
        IF(LHX.EQ.LHS ) TEM=RHO*WMX(L)/(WMUI*FCLD+teeny)
        TEM=TEM*TEM
        IF(TEM.GT.10.) TEM=10.
        CM1=CM0
        IF(BANDF) CM1=CM0*CBF
        IF(LHX.EQ.LHS) CM1=CM0
        CM=CM1*(1.-1./EXP(TEM*TEM))+100.*(PREBAR(L+1)+
     *       PRECNVL(L+1)*BYDTsrc)
#ifdef CLD_AER_CDNC
C** Choice of 2 different routines to get the autoconversion rate
#ifdef BLK_2MOM
C*** using an alternate QAUT definition based on Beheng (1994)
        CM=QAUT_B2M/(WMX(L)+1.d-20)+1.d0*100.d0*(PREBAR(L+1)+
     *     PRECNVL(L+1)*BYDTsrc)
c     if (QAUT_B2M.lt.0.) write(6,*)"QAUT BLK_2M",QAUT_B2M,CM,WMX(L),L
c       if(L.eq.1) write(6,*)"4th check BLK_2M",CM,QAUT_B2M,WMX(L)
#else
C** Use Qaut definition based on Rotstayn and Liu (2005, GRL)
         WTEM=1d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+teeny)
          IF(LHX.EQ.LHE)  THEN
            RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCW))**BY3
          ELSE
            RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCI))**BY3
     *         *(1.+pl(l)*xRICld)
          END IF
       CALL GET_QAUT(L,PL(L),TL(L),FCLD,WMX(L),SCDNCW,RCLD,RHOW,
     *r6,r6c,QCRIT,QAUT)
C** Can also use other Qaut definitions if BLK_2MOM is not defined by switching to this call
!     CALL GET_QAUT(L,TL(L),FCLD,WMX(L),SCDNCW,RHO,QCRIT,QAUT)
!      CALL GET_QAUT(L,FCLD,WMX(L),SCDNCW,RHO,QAUT)
C*** If 6th moment of DSD is greater than critical radius r6c start QAUT
!     if (r6.gt.r6c) then
      if ((WMX(L)/(FCLD+teeny)).GT.QCRIT) then
        CM=QAUT/(WMX(L)+1.d-20)+1.d0*100.d0*(PREBAR(L+1)+
     *     PRECNVL(L+1)*BYDTsrc)
      else
        CM=0.d0
      endif
C** end routine for QAUT as a function of N,LWC
#endif
#endif
        CM=CM*CMX
        IF(CM.GT.BYDTsrc) CM=BYDTsrc
        PREP(L)=WMX(L)*CM
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
        IF (CAREA(L).GT.0.) THEN
          WTEM=1d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+teeny)
          IF(LHX.EQ.LHE.AND.WMX(L)/FCLD.GE.WCONST*1d-3)
     *         WTEM=1d2*WCONST*PL(L)/(TL(L)*RGAS)
          IF(LHX.EQ.LHS.AND.WMX(L)/FCLD.GE.WMUI*1d-3)
     *         WTEM=1d2*WMUI*PL(L)/(TL(L)*RGAS)
          IF(WTEM.LT.1d-10) WTEM=1d-10
          IF(LHX.EQ.LHE)  THEN
!           RCLD=1d-6*(RWCLDOX*10.*(1.-PEARTH)+7.*PEARTH)*(WTEM*4.)**BY3
            RCLD=RCLDX*1d-6*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCW))**BY3
          ELSE
!           RCLD=25.d-6*(WTEM/4.2d-3)**BY3 * (1.+pl(l)*xRICld)
            RCLD=RCLDX*100.d-6*(WTEM/(2.d0*BY3*TWOPI*SCDNCI))**BY3
     *         *(1.+pl(l)*xRICld)
          END IF
          CK1=1000.*LHX*LHX/(2.4d-2*RVAP*TL(L)*TL(L))
          CK2=1000.*RGAS*TL(L)/(2.4d-3*QSATL(L)*PL(L))
          TEVAP=COEEC*(CK1+CK2)*RCLD*RCLD
          WMX1=WMX(L)-PREP(L)*DTsrc
          ECRATE=(1.-RHF(L))/(TEVAP*FCLD+teeny)
          IF(ECRATE.GT.BYDTsrc) ECRATE=BYDTsrc
          EC(L)=WMX1*ECRATE*LHX
        END IF
C**** COMPUTE NET LATENT HEATING DUE TO STRATIFORM CLOUD PHASE CHANGE,
C**** QHEAT, AND NEW CLOUD WATER CONTENT, WMNEW
        DRHDT=2.*CAREA(L)*CAREA(L)*(1.-RH00(L))*(QCONV+ER(L))/LHX/
     *       (WMX(L)/(FCLD+teeny)+2.*CAREA(L)*QSATL(L)*(1.-RH00(L))
     *       +teeny)
        IF(ER(L).EQ.0.AND.WMX(L).LE.0.) DRHDT=0.
        QHEAT(L)=FSSL(L)*(QCONV-LHX*DRHDT*QSATL(L))/(1.+RH(L)*SQ(L))
        DWDT=QHEAT(L)/LHX-PREP(L)+CAREA(L)*FSSL(L)*ER(L)/LHX
        WMNEW =WMX(L)+DWDT*DTsrc
        IF(WMNEW.LT.0.) THEN
          WMNEW=0.
          QHEAT(L)=(-WMX(L)*BYDTsrc+PREP(L))*LHX-CAREA(L)*FSSL(L)*ER(L)
        END IF
      ELSE
C**** UNFAVORABLE CONDITIONS FOR CLOUDS TO EXIT, PRECIP OUT CLOUD WATER
        IF(WMX(L).GT.0.) PREP(L)=WMX(L)*BYDTsrc
        ER(L)=(1.-RH(L))**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
        IF(PREICE(L+1).GT.0..AND.TL(L).LT.TF)
     *       ER(L)=(1.-RHI)**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
        ER(L)=MAX(0d0,MIN(ER(L),ERMAX))
        QHEAT(L)=-CAREA(L)*FSSL(L)*ER(L)
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
      ENDIF
C**** PHASE CHANGE OF PRECIP, FROM WATER TO ICE
      IF (LHP(L+1).EQ.LHE.AND.LHP(L).EQ.LHS.AND.PREBAR(L+1).GT.0)
     *     HPHASE=HPHASE-LHM*PREBAR(L+1)*GRAV*BYAM(L)
C**** Make sure energy is conserved for transfers between P and CLW
      IF (LHP(L).NE.LHX)
     *     HPHASE=HPHASE+(ER(L)*CAREA(L)*FSSL(L)/LHX-PREP(L))*LHM
C**** COMPUTE THE PRECIP AMOUNT ENTERING THE LAYER TOP
      IF (ER(L).eq.ERMAX) THEN ! to avoid round off problem
        PREBAR(L)=PREBAR(L+1)*(1.-CAREA(L)*FSSL(L))+
     *            AIRM(L)*PREP(L)*BYGRAV
      ELSE
        PREBAR(L)=MAX(0d0,PREBAR(L+1)+
     *       AIRM(L)*(PREP(L)-ER(L)*CAREA(L)*FSSL(L)/LHX)*BYGRAV)
      END IF
C**** UPDATE NEW TEMPERATURE AND SPECIFIC HUMIDITY
      QNEW =QL(L)-DTsrc*QHEAT(L)/(LHX*FSSL(L)+teeny)
      IF(QNEW.LT.0.) THEN
        QNEW=0.
        QHEAT(L)=QL(L)*LHX*BYDTsrc*FSSL(L)
        DWDT1=QHEAT(L)/LHX-PREP(L)+CAREA(L)*FSSL(L)*ER(L)/LHX
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
      IF (PREBAR(L+1).gt.0.) FER=CAREA(L)*FSSL(L)*ER(L)*AIRM(L)/
     *     (GRAV*LHX*PREBAR(L+1))                             ! P->Q
      FER=MIN(1d0,FER)
      FWTOQ=0.                                                ! CLW->Q
#endif
      FQTOW=0.                                                ! Q->CLW
      IF (FSSL(L).gt.0) THEN
      IF (QHEAT(L)+CAREA(L)*FSSL(L)*ER(L).gt.0) THEN
        IF (LHX*QL(L)+DTsrc*CAREA(L)*ER(L).gt.0.) FQTOW=(QHEAT(L
     *       )+CAREA(L)*FSSL(L)*ER(L))*DTsrc/((LHX*QL(L)+DTsrc*CAREA(L)
     *       *ER(L))*FSSL(L))
#ifdef TRACERS_WATER
      ELSE
        IF (WMX(L)-PREP(L)*DTsrc.gt.0.) FWTOQ=-(QHEAT(L)
     *       +CAREA(L)*FSSL(L)*ER(L))*DTsrc/(LHX*(WMX(L)-PREP(L)*DTsrc))
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
      TH(L)=TL(L)/PLK(L)
      TNEW=TL(L)
      QSATC=QSAT(TL(L),LHX,PL(L))
      RH(L)=QL(L)/QSATC
#ifdef TRACERS_WATER
C**** update tracers from cloud formation (in- and below-cloud
C****    precipitation, evaporation, condensation, and washout)
c CLDSAVT is current FCLD
        IF(RH(L).LE.1.) CLDSAVT=1.-DSQRT((1.-RH(L))/(1.-RH00(L)+teeny))
        IF(CLDSAVT.LT.0.) CLDSAVT=0.
        IF(RH(L).GT.1.) CLDSAVT=1.
        IF (CLDSAVT.GT.1.) CLDSAVT=1.
        IF (WMX(L).LE.0.) CLDSAVT=0.
        CLDSAVT=CLDSAVT*FSSL(L)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      WA_VOL=0.
      IF (WMNEW.GT.teeny) THEN
        WA_VOL=WMNEW*AIRM(L)*1.D2*BYGRAV*DXYPIJ
      ENDIF
      WMXTR = WMX(L)
      IF (BELOW_CLOUD.and.WMX(L).LT.teeny) THEN
        precip_mm = PREBAR(L+1)*100.*DTsrc
        if (precip_mm.lt.0.) precip_mm=0.
        WMXTR = PREBAR(L+1)*grav*BYAM(L)*dtsrc
        if (wmxtr.lt.0.) wmxtr=0.
        WA_VOL=precip_mm*DXYPIJ
      ENDIF
      CALL GET_SULFATE(L,TL(L),FCLD,WA_VOL
     * ,WMXTR,SULFIN,SULFINC,SULFOUT,TR_LEFT,TM,TRWML(1,l),AIRM,LHX
     *  ,DT_SULF_SS(1,L),CLDSAVT)

      DO N=1,NTX
      select case (trname(ntix(n)))
      case('SO2','SO4','H2O2_s','H2O2')
      if (trname(ntix(n)).eq."H2O2" .and. coupled_chem.eq.0) goto 403
      if (trname(ntix(n)).eq."H2O2_s" .and. coupled_chem.eq.1) goto 403

        TRWML(N,L)=TRWML(N,L)*(1.+SULFINC(N))
        TM(L,N)=TM(L,N)*(1.+SULFIN(N))
        TMOM(:,L,N)  = TMOM(:,L,N)*(1. +SULFIN(N))
        if (WMX(L).LT.teeny.and.BELOW_CLOUD) then
          TRPRBAR(N,L+1)=TRPRBAR(N,L+1)+SULFOUT(N)
        else
          TRWML(N,L) = TRWML(N,L)+SULFOUT(N)
        endif

 403    CONTINUE

      end select

      END DO
#endif
      DO N=1,NTX
c ---------------------- initialize fractions ------------------------
        FPRT  =0.
        FERT  =0.
        FWASHT=0.
        FQTOWT=0.
        FWTOQT=0.
        THLAW=0.
        THWASH=0.
c ----------------------- calculate fractions --------------------------
c precip. tracer evap
        CALL GET_EVAP_FACTOR(N,TL(L),LHX,.FALSE.,1d0,FER,FERT,ntix)
        CALL GET_EVAP_FACTOR(N,TL(L),LHX,.FALSE.,1d0,FWTOQ,FWTOQT,ntix)
        TR_LEF=1.D0
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      select case (trname(ntix(n)))
      case('SO2','SO4','H2O2_s','H2O2')
      if (trname(ntix(n)).eq."H2O2" .and. coupled_chem.eq.0) goto 404
      if (trname(ntix(n)).eq."H2O2_s" .and. coupled_chem.eq.1) goto 404

        TR_LEF=TR_LEFT(N)
 404    CONTINUE

      end select

#endif
        IF(BELOW_CLOUD.and.WMX(L).lt.teeny) THEN
          precip_mm = PREBAR(L+1)*100.*dtsrc
          WMXTR = PREBAR(L+1)*grav*BYAM(L)*dtsrc
          if (precip_mm.lt.0.) precip_mm=0.
          if (wmxtr.lt.0.) wmxtr=0.
cdmk change GET_WASH below - extra arguments
          CALL GET_WASH_FACTOR(N,b_beta_DT,precip_mm,FWASHT
     *     ,tl(l),LHX,WMXTR,cldprec,L,TM,TRPRBAR(1,l),THWASH,pl(l),ntix) !washout
        ELSE
          WMXTR = WMX(L)
c         b_beta_DT is needed at the lowest precipitating level,
c         so saving it here for below cloud case:
          b_beta_DT = cldsavt*CM*dtsrc
          CALL GET_COND_FACTOR(L,N,WMXTR,TL(L),TL(L),LHX,FCLD,FQTOW
     *         ,FQTOWT,.false.,TRWML,TM,THLAW,TR_LEF,PL(L),ntix,CLDSAVT)
cdmk added arguments above; THLAW added below (no way to factor this)
        END IF
        IF (TM(L,N).GT.teeny) THEN
          TMFAC=THLAW/TM(L,N)
        ELSE
          TMFAC=0.
        ENDIF
!gone:  CALL GET_PREC_FACTOR(N,BELOW_CLOUD,CM,CLDSAVT,FPR,FPRT,ntix) !precip CLM
        FPRT=FPR
c ---------------------- calculate fluxes ------------------------
        DTERT = FERT  *TRPRBAR(N,L+1)
        DTPRT = FPRT  *TRWML(N,L)
        DTQWT =
     &  FQTOWT*TR_LEF*(TM(L,N)+DTERT)-FWTOQT*TRWML(N,L)*(1.-FPRT)
#ifndef NO_WASHOUT_IN_CLOUDS
c**** washout in clouds
        tm_temp(l,n)=tm(l,n)-dtqwt-thlaw
        IF(.NOT.(BELOW_CLOUD.and.WMX(L).lt.teeny)) THEN
          precip_mm = prebar(l+1)*100.*dtsrc
          wmxtr=prebar(l+1)*grav*byam(l)*dtsrc
          IF (precip_mm < 0.) precip_mm=0.
          IF (wmxtr < 0.) wmxtr=0.
          CALL get_wash_factor(n,b_beta_dt,precip_mm,fwasht,tl(l),lhx,
     &         wmxtr,cldprec,l,tm_temp,trprbar(1,l),thwash,pl(l),ntix) !washout
c         saves cloud fraction at lowest precipitating level for washout
          cldprec=cldsavt
        END IF
        dtwrt=fwasht*tm_temp(l,n)
        IF (tm_temp(l,n).GT.teeny) THEN
          TMFAC2=THWASH/tm_temp(l,n)
        ELSE
          TMFAC2=0.
        ENDIF
#else
        dtwrt=fwasht*tm(l,n)
        IF (TM(L,N).GT.teeny) THEN
          TMFAC2=THWASH/TM(L,N)
        ELSE
          TMFAC2=0.
        ENDIF
#endif
c ---------------------- apply fluxes ------------------------
        TRWML(N,L) = TRWML(N,L)*(1.-FPRT)  + DTQWT+THLAW
#ifdef TRDIAG_WETDEPO
        IF (diag_wetdep == 1) THEN
          trevap_ls(l,n)=dtert
          trwash_ls(l,n)=dtwrt+thwash
          trclwc_ls(l,n)=fqtowt*tr_lef*(tm(l,n)+dtert)+thlaw
          trprcp_ls(l,n)=dtprt
          trclwe_ls(l,n)=fwtoqt*trwml(n,l)*(1.-fprt)
        END IF
#endif
        TM(L,N)    = TM(L,N)  + DTERT - DTWRT - DTQWT - THLAW - THWASH
        TRPRBAR(N,L)=TRPRBAR(N,L+1)*(1.-FERT) + DTPRT+DTWRT+THWASH
        IF (PREBAR(L).eq.0) TRPRBAR(N,L)=0.  ! remove round off error
        IF (WMX(L).eq.0) TRWML(N,L)=0.       ! remove round off error
        TMOM(:,L,N)  = TMOM(:,L,N)*(1. - FQTOWT - FWASHT
     *    - TMFAC - TMFAC2)
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
#endif
C**** CONDENSE MORE MOISTURE IF RELATIVE HUMIDITY .GT. 1
      IF(RH(L).GT.1.) THEN
      SLH=LHX*BYSHA
      DQSUM=0.
      DO N=1,3
        IF(N.NE.1) QSATC=QSAT(TL(L),LHX,PL(L))
        DQ=(QL(L)-QSATC)/(1.+SLH*QSATC*DQSATDT(TL(L),LHX))
        TL(L)=TL(L)+SLH*DQ
        QL(L)=QL(L)-DQ
        DQSUM=DQSUM+DQ
      END DO
      IF(DQSUM.GT.0.) THEN
      WMX(L)=WMX(L)+DQSUM*FSSL(L)
      FCOND=DQSUM/QNEW
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      WA_VOL=0.
      IF (WMX(L).GT.teeny) THEN
        WA_VOL=WMX(L)*AIRM(L)*1.D2*BYGRAV*DXYPIJ
      ENDIF
#endif
C**** adjust gradients down if Q decreases
      QMOM(:,L)= QMOM(:,L)*(1.-FCOND)
#ifdef TRACERS_WATER
C**** CONDENSING MORE TRACERS
      WMXTR = WMX(L)
        IF(RH(L).LE.1.) CLDSAVT=1.-DSQRT((1.-RH(L))/(1.-RH00(L)+teeny))
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
#endif
      DO N=1,NTX
        TR_LEF=1.
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      select case (trname(ntix(n)))
      case('SO2','SO4','H2O2_s','H2O2')
      if (trname(ntix(n)).eq."H2O2" .and. coupled_chem.eq.0) goto 405
      if (trname(ntix(n)).eq."H2O2_s" .and. coupled_chem.eq.1) goto 405

        TRWML(N,L)=TRWML(N,L)*(1.+SULFINC(N))
        TM(L,N)=TM(L,N)*(1.+SULFIN(N))
        TMOM(:,L,N) =TMOM(:,L,N)*(1.+SULFIN(N))
        TRWML(N,L) = TRWML(N,L)+SULFOUT(N)
        TR_LEF=TR_LEFT(N)

 405    CONTINUE

      end select

#endif
c below TR_LEFT(N) limits the amount of available tracer in gridbox
cdmkf and below, extra arguments for GET_COND, addition of THLAW
        CALL GET_COND_FACTOR(L,N,WMXTR,TL(L),TL(L),LHX,FCLD,FCOND
     *       ,FQCONDT,.false.,TRWML,TM,THLAW,TR_LEF,pl(l),ntix,CLDSAVT)
        IF (TM(L,N).GT.teeny) THEN
          TMFAC=THLAW/TM(L,N)
        ELSE
          TMFAC=0.
        ENDIF
        TRWML(N,L)  =TRWML(N,L)+ FQCONDT*TM(L,N)+THLAW
        TM(L,N)     =TM(L,N)    *(1.-FQCONDT)   -THLAW
        TMOM(:,L,N) =TMOM(:,L,N)*(1.-FQCONDT - TMFAC)
#ifdef TRDIAG_WETDEPO
        IF (diag_wetdep == 1) trcond_ls(l,n)=fqcondt*tm(l,n)+thlaw
#endif
      END DO
#endif
      ELSE
      TL(L)=TNEW
      QL(L)=QNEW
      END IF
      RH(L)=QL(L)/QSAT(TL(L),LHX,PL(L))
      TH(L)=TL(L)/PLK(L)
      TNEW=TL(L)
!     if (debug_out) write(0,*) 'after condensation: l,tlnew,',l,tl(l)
      END IF
      IF(RH(L).LE.1.) CAREA(L)=DSQRT((1.-RH(L))/(1.-RH00(L)+teeny))
      IF(CAREA(L).GT.1.) CAREA(L)=1.
      IF(RH(L).GT.1.) CAREA(L)=0.
      IF(WMX(L).LE.0.) CAREA(L)=1.
      IF(CAREA(L).LT.0.) CAREA(L)=0.
      RHF(L)=RH00(L)+(1.-CAREA(L))*(1.-RH00(L))
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
          TL(L)=TL(L)+HCHANG/(SHA*FSSL(L)+teeny)
!         if(debug_out) write(0,*) 'after rain out: l,tlnew',l,tl(l)
          TH(L)=TL(L)/PLK(L)
        END IF
        WMX(L)=0.
      END IF
      prebar1(l)=prebar(l)
C**** set phase of condensation for next box down
      PREICE(L)=0.
      IF (PREBAR(L).gt.0 .AND. LHP(L).EQ.LHS) PREICE(L)=PREBAR(L)
      IF (PREBAR(L).le.0) LHP(L)=0.
C**** COMPUTE THE LARGE-SCALE CLOUD COVER
      IF(RH(L).LE.1.) CAREA(L)=DSQRT((1.-RH(L))/(1.-RH00(L)+teeny))
      IF(CAREA(L).GT.1.) CAREA(L)=1.
      IF(RH(L).GT.1.) CAREA(L)=0.
      IF(WMX(L).LE.0.) CAREA(L)=1.
      IF(CAREA(L).LT.0.) CAREA(L)=0.
      CLDSSL(L)=FSSL(L)*(1.-CAREA(L))
      CLDSAVL(L)=1.-CAREA(L)
#ifdef TRACERS_WATER
      IF(CLDSSL(L).gt.0.) CLOUD_YET=.true.
      IF(CLOUD_YET.and.CLDSSL(L).eq.0.) BELOW_CLOUD=.true.
#endif
      TOLDUP=TOLD
C**** ACCUMULATE SOME DIAGNOSTICS
         HCNDSS=HCNDSS+FSSL(L)*(TNEW-TOLD)*AIRM(L)
         SSHR(L)=SSHR(L)+FSSL(L)*(TNEW-TOLD)*AIRM(L)
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
        TOLD=TL(L)
        TOLDU=TL(L+1)
        QOLD=QL(L)
        QOLDU=QL(L+1)
        FCLD=(1.-CAREA(L))*FSSL(L)+teeny
        IF(CAREA(L).EQ.1. .OR. (CAREA(L).LT.1..AND.CAREA(L+1).LT.1.))
     *       CYCLE
        SEDGE=THBAR(TH(L+1),TH(L))
        DSE=(TH(L+1)-SEDGE)*PLK(L+1)+(SEDGE-TH(L))*PLK(L)+
     *       SLHE*(QL(L+1)-QL(L))
        DWM=QL(L+1)-QL(L)+(WMX(L+1)-WMX(L))/FCLD
        DQSDT=DQSATDT(TL(L),LHE)*QL(L)/(RH(L)+1d-30)
        BETA=(1.+BYMRAT*TL(L)*DQSDT)/(1.+SLHE*DQSDT)
        CKM=(1.+SLHE*DQSDT)*(1.+(1.-DELTX)*TL(L)/SLHE)/
     *       (2.+(1.+BYMRAT*TL(L)/SLHE)*SLHE*DQSDT)
        CKR=TL(L)/(BETA*SLHE)
        CK=DSE/(SLHE*DWM)
        SIGK=0.
        IF(CKR.GT.CKM) CYCLE
        IF(CK.GT.CKR) SIGK=COESIG*((CK-CKR)/(CKM-CKR+teeny))**5
        EXPST=EXP(-SIGK*DTsrc)
        IF(L.LE.1) CKIJ=EXPST
        DSEC=DWM*TL(L)/BETA
        IF(CK.LT.CKR) CYCLE
        FPMAX=MIN(1d0,1.-EXPST)
        IF(FPMAX.LE.0.) CYCLE
        IF(DSE.GE.DSEC) CYCLE
C**** MIXING TO REMOVE CLOUD-TOP ENTRAINMENT INSTABILITY
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
        ENDDO
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
          IF(DSEDIF.GT.1d-3) FPLUME=FPLUME-DFX
          IF(DSEDIF.LT.-1d-3) FPLUME=FPLUME+DFX
          IF(ABS(DSEDIF).LE.1d-3.OR.FPLUME.GT.FPMAX) EXIT
        END DO
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
        TH(L+1)=TH(L+1)-(LHX*BYSHA)*WMX(L+1)/(PLK(L+1)*FSSL(L)+teeny)
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
        IF(RH(L).LE.1.) CAREA(L)=DSQRT((1.-RH(L))/(1.-RH00(L)+teeny))
        IF(CAREA(L).GT.1.) CAREA(L)=1.
        IF(RH(L).GT.1.) CAREA(L)=0.
        CLDSSL(L)=FSSL(L)*(1.-CAREA(L))
        CLDSAVL(L)=1.-CAREA(L)
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
        IF(LHX.EQ.LHS.AND.WMX(L)/FCLD.GE.WMUI*1d-3)
     *       WTEM=1d5*WMUI*1.d-3*PL(L)/(TL(L)*RGAS)
        IF(WTEM.LT.1d-10) WTEM=1.d-10
C***Setting constant values of CDNC over land and ocean to get RCLD=f(CDNC,LWC)
      SNdO = 59.68d0/(RWCLDOX**3)
      SNdL = 174.d0
      SNdI = 0.06417127d0
      SCDNCW=SNdO*(1.-PEARTH)+SNdL*PEARTH
      SCDNCI=SNdI
#ifdef CLD_AER_CDNC
!@auth Menon for CDNC prediction
      CALL GET_CDNC_UPD(L,LHX,WCONST,WMUI,WMX(L),FCLD,CLDSSL(L),
     *CLDSAVL(L),VVEL,SME(L),DSU,OLDCDL(L),
     *CDNL0,CDNL1)
      OLDCDL(L) = CDNL1
      SNd=CDNL1
C** Pass old and new cloud droplet number
      NEWCDN=SNd
      OLDCDN=CDNL0
c     if (L.eq.1)write(6,*)"BLK_2M NUPD",NEWCDN,OLDCDN
#endif
#ifdef BLK_2MOM
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
c        mcrys=WMICE(L)         ! crys content, [kg water/kg air]
c        ncrys=mcrys/mi0         ! crys concent, [No/m3]
      IF(LHX.EQ.LHE)  THEN
         mdrop =WMX(L)
         ndrop= OLDCDL(L)*1.d6  !mdrop/mw0         ! drop concent, [No/m3]
         if(WMX(L).eq.0.) ndrop=0.0
       ELSE
         mcrys =WMX(L)
         WMICE(L) = WMX(L)
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
C*** Call Lohmann's or Gultepe's scheme for CDNC
        OLDCDNC=OLDCDN*1.d6  !convert from cm-3 to m-3
        NEWCDNC=NEWCDN*1.d6  !convert from cm-3 to m-3
C***
#ifdef TRACERS_AMP
C*** Using the AMP_actv interface from MATRIX
        ldummy=execute_bulk2m_driver('matr','drop_nucl',dtB2M,mkx)
#else
        ldummy=execute_bulk2m_driver('gult','drop_nucl',dtB2M,mkx,
     &OLDCDNC,NEWCDNC)
#endif

      rablk=execute_bulk2m_driver('get','value','nc') + (
     *  +execute_bulk2m_driver('get','npccn')
     *                                       )*dtB2M
      SNd=rablk(mkx)*1.0d-6            ! ndrop, [No/m^3], SNdL, [No/cc]
C**** Old treatment to get CDNC for cloud changes within time steps
       DCLD(L) = CLDSSL(L)-CLDSAVL(L) ! cloud fraction change
C** If previous time step is clear sky
       if(CLDSAVL(L).eq.0.) then
         SNd=SNd
       elseif (DCLD(L).le.0.d0) then
         SNd=OLDCDL(L)
       elseif(DCLD(L).gt.0.d0) then
         SNd=( (OLDCDL(L)*CLDSAVL(L)) + (SNd*DCLD(L)) )/CLDSSL(L)
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

       SNdI=rablk(mkx)*1.0d-6             ! from ncrys [No/m^3] to SNdI in [No/cc]
       if(SNDI.gt.1.d0) SNdI=1.d0      !try to limit to 1000 /l
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
      SCDNCI=SNdI
      If (SCDNCW.le.20.d0) SCDNCW=20.d0   !set min CDNC sensitivity test
c     If (SCDNCI.le.0.06d0) SCDNCI=0.06417127d0   !set min ice crystal
      if(SCDNCW.gt.1400.d0) SCDNCw=1400.d0
c     write(6,*) "SCND CDNC",SCDNCW,OLDCDL(l),OLDCDO(l),l
#endif

        IF(LHX.EQ.LHE) THEN

!         RCLD=(RWCLDOX*10.*(1.-PEARTH)+7.0*PEARTH)*(WTEM*4.)**BY3
          RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCW))**BY3
          QHEATC=(QHEAT(L)+FSSL(L)*CAREA(L)*(EC(L)+ER(L)))/LHX
          IF(RCLD.GT.20..AND.PREP(L).GT.QHEATC) RCLD=20.
          RCLDE=RCLD/BYBR
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
     *         *(1.+pl(l)*xRICld)
          RCLDE=RCLD/BYBR
#ifdef BLK_2MOM
c       if(L.eq.1)  write(6,*)"9th check BLK_2M",RCLDE
!        rablk=execute_bulk2m_driver('get','value','ei')  ! [micron]
!        RCLDE=rablk(mkx)
c        if(l.eq.1) write(6,*)"10th check BLK_2M",RCLDE
#endif
        ENDIF
#ifdef CLD_AER_CDNC
C** Using the Liu and Daum paramet
C** for spectral dispersion effects on droplet size distribution
      Repsi=1.d0 - 0.7d0*exp(-0.003d0*SCDNCW)
      Repsis=Repsi*Repsi
      Rbeta=(((1.d0+2.d0*Repsis)**0.667d0))/((1.d0+Repsis)**0.333d0)
!     write(6,*)"RCLD",Rbeta,RCLD,SCDNCW,Repsis
      RCLDE=RCLD*Rbeta
!@auth Menon    end of addition  comment out the RCLDE definition below
#endif
        CSIZEL(L)=RCLDE
#ifdef CLD_AER_CDNC
C**** save for diag purposes
        IF (FCLD.gt.1.d-5.and.LHX.eq.LHE) then
            ACDNWS(L)= SCDNCW   ! cloud droplets in cm-3
            AREWS(L) = RCLDE
            ALWWS(L) = WTEM
            CDN3DL(L)=SCDNCW
            CRE3DL(L)=RCLDE
            NLSW  = NLSW + 1
c      if(NLSW.ge.1.and.l.eq.1) write(6,*)"INCLD",ACDNWS(L),
c    * SCDNCW,NLSW
        elseif(FCLD.gt.1.d-5.and.LHX.eq.LHS) then
            ACDNIS(L)= SCDNCI     ! ice crysal in cm-3, *1.d-3 for l-1
            AREIS(L) = RCLDE
            ALWIS(L) = WTEM
            CDN3DL(L)=SCDNCI
            CRE3DL(L)=RCLDE
            NLSI  = NLSI + 1
        ENDIF
#endif
        TEM=AIRM(L)*WMX(L)*1.d2*BYGRAV
        TAUSSL(L)=1.5d3*TEM/(FCLD*RCLDE+teeny)
        IF(TAUSSL(L).GT.100.) TAUSSL(L)=100.
        IF(LHX.EQ.LHE) WMSUM=WMSUM+TEM
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
            CLDSSL(L)=CLDSSL(L)+(BMAX-CLDSSL(L))*CKIJ
            TAUSSL(L)=TAUSSL(L)*CLDSV1(L)/(CLDSSL(L)+teeny)
          ENDIF
          IF(TAUSSL(L).LE.0.) CLDSSL(L)=0.
          IF(L.GT.DCL .AND. TAUMCL(L).LE.0.) THEN
            CLDSSL(L)=CLDSSL(L)**(2.*BY3)
            TAUSSL(L)=TAUSSL(L)*CLDSV1(L)**BY3
          END IF
        END IF
      END DO

#ifdef CLD_AER_CDNC
!Save variables for 3 hrly diagnostics

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
       IF (SVLHXL(L).EQ.LHE) CL3DL(L) = WMX(L)*RHODK*CD3DL(L) ! cld water kg m-2
       IF (SVLHXL(L).EQ.LHS) CI3DL(L) = WMX(L)*RHODK*CD3DL(L) ! ice water kg m-2
c      write(6,*)"CT",L,WMX(L),CD3DL(l),CL3DL(L),CI3DL(L)

      END DO

#endif

      RETURN
      END SUBROUTINE LSCOND

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
!$Id: CLOUDS2_E1.f,v 1.41 2010/11/04 23:01:33 cdrar Exp $
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
      USE CLOUDS, only : ncol,tautab,invtau
      implicit none
!@var  emsfc_lw    longwave emissivity of surface at 10.5 microns
      REAL*8, PARAMETER :: emsfc_lw=0.99d0
!@var pc1bylam Planck constant c1 by wavelength (10.5 microns)
      REAL*8, PARAMETER :: pc1bylam = 1.439d0/10.5d-4
!@var boxarea fractional area of each sub-grid scale box
      REAL*8, PARAMETER :: boxarea = 1d0/ncol
!@var t0 ave temp  (K)
      REAL*8, PARAMETER :: t0 = 296.
!@var t0bypstd ave temp by sea level press
      REAL*8, PARAMETER :: t0bypstd = t0/1.013250d6
      real*8, parameter :: bywc = 1./2.56d0 , byic= 1./2.13d0
      real*8, parameter :: isccp_taumin = 0.3d0  ! used to be 0.1

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
            ENDIF
c            IF (ncolprint.ne.0) then
c              write (6,'(a)') 'threshold_nsf2:'
c                do j=1,npoints,1000
c                write(6,'(a10)') 'j='
c                write(6,'(8I10)') j
c                write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
c                enddo
c            ENDIF
        ENDIF

c        IF (ncolprint.ne.0) then
c            write (6,'(a)') 'ilev:'
c            write (6,'(I2)') ilev
c        ENDIF

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

        ENDDO ! ibox

!          Fill frac_out with 1's where tca is greater than the threshold

           DO ibox=1,ncol
             do j=1,npoints
               if (tca(j,ilev).gt.threshold(j,ibox)) then
               frac_out(j,ibox,ilev)=1
               else
               frac_out(j,ibox,ilev)=0
               end if
             enddo
           ENDDO

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
           ENDDO

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


        !---------------------------------------------------------------
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
        !---------------------------------------------------------------



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

