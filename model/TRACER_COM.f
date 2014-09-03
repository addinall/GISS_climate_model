#include "rundeck_opts.h"
      

!@sum  TRACER_COM: Exists alone to minimize the number of dependencies
!@+    This version for simple trace gases/chemistry/isotopes
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      MODULE TRACER_COM
!@sum  TRACER_COM tracer variables
!@auth Jean Lerner
!@ver  1.0
C
      USE QUSDEF, only: nmom
      USE MODEL_COM, only: im,jm,lm
C
      IMPLICIT NONE
      SAVE

C**** Each tracer has a variable name and a unique index
!@param NTM number of tracers
!@var TRNAME: Name for each tracer >>> MUST BE LEFT-JUSTIFIED <<<

!@var ntm_O18: Number of TRACERS_SPECIAL_O18 tracers.
#ifdef TRACERS_SPECIAL_O18
      integer, parameter :: ntm_o18=2
#else
      integer, parameter :: ntm_o18=0
#endif  /* TRACERS_SPECIAL_O18 */
!@var ntm_gasexch: Number of TRACERS_GASEXCH_ocean tracers.
#if defined(TRACERS_GASEXCH_ocean) || defined(TRACERS_GASEXCH_land)
      integer, parameter :: ntm_gasexch=1
#else
      integer, parameter :: ntm_gasexch=0
#endif  /* TRACERS_GASEXCH_ocean */
!@var ntm_lerner: Number of TRACERS_SPECIAL_Lerner tracers.
#ifdef TRACERS_SPECIAL_Lerner
      integer, parameter :: ntm_lerner=9
#else
      integer, parameter :: ntm_lerner=0
#endif  /* TRACERS_SPECIAL_Lerner */
!@var ntm_water: Number of TRACERS_WATER tracers.
#if (defined TRACERS_WATER) 
      integer, parameter :: ntm_water=1
#else
      integer, parameter :: ntm_water=0
#endif  /* TRACERS_WATER */
!@var ntm_shindell_trop: Number of TRACERS_SPECIAL_Shindell tracers.
#ifdef TRACERS_SPECIAL_Shindell
      integer, parameter :: ntm_shindell_trop=15
#else
      integer, parameter :: ntm_shindell_trop=0
#endif  /* TRACERS_SPECIAL_Shindell */
!@var ntm_terp: Number of TRACERS_TERP tracers.
#ifdef TRACERS_TERP
      integer, parameter :: ntm_terp=1
#else
      integer, parameter :: ntm_terp=0
#endif  /* TRACERS_TERP */
!@var ntm_shindell_strat: Number of SHINDELL_STRAT_CHEM tracers.
#ifdef SHINDELL_STRAT_CHEM
      integer, parameter :: ntm_shindell_strat=10
#else
      integer, parameter :: ntm_shindell_strat=0
#endif  /* SHINDELL_STRAT_CHEM */
!@var ntm_shindell_extra: Number of SHINDELL_STRAT_EXTRA tracers.
#ifdef SHINDELL_STRAT_EXTRA
#ifdef ACCMIP_LIKE_DIAGS
      integer, parameter :: ntm_shindell_extra=3
#else
      integer, parameter :: ntm_shindell_extra=1
#endif /* ACCMIP_LIKE_DIAGS*/
#else
      integer, parameter :: ntm_shindell_extra=0
#endif  /* SHINDELL_STRAT_EXTRA */
!@var ntm_koch: Number of TRACERS_AEROSOLS_Koch tracers.
#ifdef TRACERS_AEROSOLS_Koch
#ifdef SULF_ONLY_AEROSOLS
      integer, parameter :: ntm_koch=5
#else
      integer, parameter :: ntm_koch=13
#endif  /* SULF_ONLY_AEROSOLS */
#else
      integer, parameter :: ntm_koch=0
#endif  /* TRACERS_AEROSOLS_Koch */

!@var ntm_dust: Number of dust aerosol tracers.
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
#if (defined TRACERS_DUST) || (defined TRACERS_AMP)
#ifdef TRACERS_DUST_Silt4
      integer, parameter :: ntm_dust=5
#else
      integer, parameter :: ntm_dust=4
#endif  /* TRACERS_DUST_Silt4 */
#else /* TRACERS_MINERALS || TRACERS_QUARZHEM */
!@var ntm_minerals: Number of TRACERS_MINERALS tracers.
#ifdef TRACERS_MINERALS
      integer, parameter :: ntm_minerals = 20
#else
      integer, parameter :: ntm_minerals = 0
#endif  /* TRACERS_MINERALS */
!@var ntm_quarzhem: Number of TRACERS_QUARZHEM tracers.
#ifdef TRACERS_QUARZHEM
      integer, parameter :: ntm_quarzhem = 3
#else
      integer, parameter :: ntm_quarzhem = 0
#endif  /* TRACERS_QUARZHEM */
      integer, parameter :: ntm_dust = ntm_minerals + ntm_quarzhem
#endif /* TRACERS_MINERALS || TRACERS_QUARZHEM */
#else
      integer, parameter :: ntm_dust=0
#endif

!@var ntm_het: Number of TRACERS_HETCHEM tracers.
#ifdef TRACERS_HETCHEM
#ifdef TRACERS_NITRATE
      integer, parameter :: ntm_het=6
#else
      integer, parameter :: ntm_het=3
#endif  /* TRACERS_NITRATE */
#else
      integer, parameter :: ntm_het=0
#endif  /* TRACERS_HETCHEM */
!@var ntm_nitrate: Number of TRACERS_NITRATE tracers.
#ifdef TRACERS_NITRATE
      integer, parameter :: ntm_nitrate=3
#else
      integer, parameter :: ntm_nitrate=0
#endif  /* TRACERS_NITRATE */
!@param ntm_soa: Number of SOA tracers.
#ifdef TRACERS_AEROSOLS_SOA
!@+            Index g means gas phase, index p means particulate
!@+            (aerosol) phase. g should ALWAYS be EXACTLY BEFORE
!@+            the corresponding p and all SOA species should be
!@+            one after the other. Do not forget to declare them
!@+            to the chemistry species as well.
#ifdef TRACERS_TERP
      integer, parameter :: ntm_soa=8
#else
      integer, parameter :: ntm_soa=4
#endif  /* TRACERS_TERP */
!@param nsoa the total number of aerosol-phase SOA related species (=ntm_soa/2)
      integer, parameter :: nsoa=ntm_soa/2
#else
      integer, parameter :: ntm_soa=0
#endif  /* TRACERS_AEROSOLS_SOA */
!@var ntm_nitrate: Number of TRACERS_COSMO tracers.
!@var ntm_nitrate: Number of TRACERS_RADON tracers.
#ifdef TRACERS_COSMO
#ifdef TRACERS_RADON
      integer, parameter :: ntm_cosmo=4
#else
      integer, parameter :: ntm_cosmo=2
#endif  /* TRACERS_RADON */
#else
      integer, parameter :: ntm_cosmo=0
#endif  /* TRACERS_COSMO */
!@var ntm_om_sp: Number of TRACERS_OM_SP tracers.
#ifdef TRACERS_OM_SP
      integer, parameter :: ntm_om_sp=7
#else
      integer, parameter :: ntm_om_sp=0
#endif  /* TRACERS_OM_SP */
!@var ntm_ocean: Number of TRACERS_OCEAN tracers.
#ifdef TRACERS_OCEAN
      integer, parameter :: ntm_ocean=0
#else
      integer, parameter :: ntm_ocean=0
#endif  /* TRACERS_OCEAN */
!@var ntm_air: Number of TRACERS_AIR tracers.
#if defined TRACERS_AIR || defined HTAP_LIKE_DIAGS
      integer, parameter :: ntm_air=1
#else
      integer, parameter :: ntm_air=0
#endif  /* TRACERS_AIR */
#ifdef TRACERS_AMP
#ifdef TRACERS_AMP_M1
      integer, parameter :: ntmAMP=51
      integer, parameter :: ntm_amp=5+ntmAMP
#endif  /* TRACERS_AMP_M1 */
#ifdef TRACERS_AMP_M2
      integer, parameter :: ntmAMP=51
      integer, parameter :: ntm_amp=5+ntmAMP
#endif  /* TRACERS_AMP_M2 */
#ifdef TRACERS_AMP_M3
      integer, parameter :: ntmAMP=41
      integer, parameter :: ntm_amp=5+ntmAMP
#endif  /* TRACERS_AMP_M3 */
#ifdef TRACERS_AMP_M4
      integer, parameter :: ntmAMP=34
      integer, parameter :: ntm_amp=5+ntmAMP
#endif  /* TRACERS_AMP_M4 */
#ifdef TRACERS_AMP_M5
      integer, parameter :: ntmAMP=45
      integer, parameter :: ntm_amp=5+ntmAMP
#endif  /* TRACERS_AMP_M5 */
#ifdef TRACERS_AMP_M6
      integer, parameter :: ntmAMP=45
      integer, parameter :: ntm_amp=5+ntmAMP
#endif  /* TRACERS_AMP_M6 */
#ifdef TRACERS_AMP_M7
      integer, parameter :: ntmAMP=35
      integer, parameter :: ntm_amp=5+ntmAMP
#endif  /* TRACERS_AMP_M7 */
#ifdef TRACERS_AMP_M8
      integer, parameter :: ntmAMP=28
      integer, parameter :: ntm_amp=5+ntmAMP
#endif  /* TRACERS_AMP_M8 */
#else
      integer, parameter :: ntm_amp=0
#endif  /* TRACERS_AMP */

!@param ntm_chem number of drew-only tracers
      integer, parameter :: ntm_chem=ntm_shindell_trop+
     *                               ntm_terp+
     *                               ntm_shindell_strat+
     *                               ntm_soa
#ifdef TRACERS_AMP 
! This is kept seperate, as ntm_dust needs to be set 
c          (in order to calculate dust emissions), but not added to ntm.
      integer, parameter :: ntm=ntm_amp
#else
!@param ntm number of tracers
      integer, parameter :: ntm=ntm_O18+ntm_gasexch+ntm_lerner+
     *                          ntm_water+ntm_koch+ntm_dust+ntm_het+
     *                          ntm_nitrate+ntm_cosmo+ntm_om_sp+
     *                          ntm_ocean+ntm_air+ntm_chem+
     *                          ntm_shindell_extra
#endif
C**** Each tracer has a variable name and a unique index
C**** The chemistry species need to be declared first, until the
C**** do igas=1,ntm_chem instances get corrected.
!@var trname_p1: dummy matrix for the definitions only
!@var trname: Name for each tracer >>> MUST BE LEFT-JUSTIFIED <<<
      character*8, parameter :: trname_p1(ntm+1)=(/
#ifdef TRACERS_SPECIAL_Shindell
     *    'Ox      ','NOx     ',
#ifdef SHINDELL_STRAT_CHEM
     *    'ClOx    ','BrOx    ',
#endif  /* SHINDELL_STRAT_CHEM */
     *                          'N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
#ifdef TRACERS_TERP
     *    'Terpenes',
#endif  /* TRACERS_TERP */
#ifdef TRACERS_AEROSOLS_SOA
     *    'isopp1g ','isopp1a ','isopp2g ','isopp2a ',
#ifdef TRACERS_TERP
     *    'apinp1g ','apinp1a ','apinp2g ','apinp2a ',
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef SHINDELL_STRAT_CHEM
     *                          'HCl     ','HOCl    ','ClONO2  ',
     *    'HBr     ','HOBr    ','BrONO2  ','N2O     ','CFC     ',
#ifdef SHINDELL_STRAT_EXTRA
#ifdef ACCMIP_LIKE_DIAGS
     *    'codirect','stratOx ','GLT     ',
#else
     *    'GLT     ', ! used to also do Be7, Be10
#endif /* ACCMIP_LIKE_DIAGS */
#else
!kt     *    'Water   ',
#endif  /* SHINDELL_STRAT_EXTRA */
#endif  /* SHINDELL_STRAT_CHEM */
#endif  /* TRACERS_SPECIAL_Shindell */
#ifndef TRACERS_AMP
#ifdef TRACERS_WATER
     *    'Water   ',
#endif  /* TRACERS_WATER */
#endif
#ifdef TRACERS_SPECIAL_O18
     *     'H2O18   ','HDO     ',   !'H2O17   ',
#endif  /* TRACERS_SPECIAL_O18 */
#ifdef TRACERS_GASEXCH_ocean_CFC
     *     'CFCn    ',
#endif  /* TRACERS_GASEXCH_ocean_CFC */
#if defined(TRACERS_GASEXCH_ocean_CO2) || defined(TRACERS_GASEXCH_land_CO2)
     *     'CO2n    ',
#endif  /* TRACERS_GASEXCH_ocean_CO2 */
#ifdef TRACERS_SPECIAL_Lerner
     *     'SF6     ','Rn222   ','CO2     ','N2O     ',
     *     'CFC11   ','14CO2   ','CH4     ','O3      ','SF6_c   ',
#endif  /* TRACERS_SPECIAL_LERNER */
#ifdef TRACERS_AEROSOLS_Koch
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
#ifndef SULF_ONLY_AEROSOLS
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ',
#endif  /* SULF_ONLY_AEROSOLS */
#endif  /* TRACERS_AEROSOLS_Koch */
#ifdef TRACERS_DUST
     *    'Clay    ','Silt1   ','Silt2   ','Silt3   ',
#ifdef TRACERS_DUST_Silt4
     *    'Silt4   ',
#endif  /* TRACERS_DUST_Silt4 */
#endif  /* TRACERS_DUST */
#ifdef TRACERS_NITRATE
     *    'NH3     ','NH4     ','NO3p    ',
#endif  /* TRACERS_NITRATE */
#ifdef TRACERS_HETCHEM
     *    'SO4_d1  ','SO4_d2  ','SO4_d3  ',
#ifdef TRACERS_NITRATE
     *    'N_d1    ','N_d2    ','N_d3    ',
#endif  /* TRACERS_NITRATE */
#endif  /* TRACERS_HETCHEM */
#ifdef TRACERS_COSMO
#ifdef TRACERS_RADON
     *    'Pb210   ',
#endif  /* TRACERS_RADON */
     *    'Be7     ','Be10    ',
#ifdef TRACERS_RADON
     *               'Rn222   ',
#endif  /* TRACERS_RADON */
#endif  /* TRACERS_COSMO */
#ifdef TRACERS_OM_SP
     *    'OCI1    ','OCI2    ','OCI3    ','OCA1    ','OCA2    ',
     *    'OCA3    ','OCA4    ',
#endif  /* TRACERS_OM_SP */
#ifdef TRACERS_MINERALS
     *     'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     *     'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     *     'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     *     'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
#endif  /* TRACERS_MINERALS */
#ifdef TRACERS_QUARZHEM
     *     'Sil1QuHe','Sil2QuHe','Sil3QuHe',
#endif  /* TRACERS_QUARZHEM */
#ifdef TRACERS_OCEAN
!#ifdef TRACERS_AGE_OCEAN
!     *     'Age     ',
!#endif  /* TRACERS_AGE_OCEAN */
#endif  /* TRACERS_OCEAN */
#if defined TRACERS_AIR || defined HTAP_LIKE_DIAGS
     *     'Air     ',
#endif  /* TRACERS_AIR */
#ifdef TRACERS_AMP
! The order of the ntmAMP aerosols matters!!!
     *    'M_NO3   ','M_NH4   ','M_H2O   ',
#if defined TRACERS_AMP_M1 || defined TRACERS_AMP_M2 || defined TRACERS_AMP_M3 \
 || defined TRACERS_AMP_M5 || defined TRACERS_AMP_M6 || defined TRACERS_AMP_M7
     *    'M_AKK_SU','N_AKK_1 ',                                 !AKK
#endif
     *    'M_ACC_SU','N_ACC_1 ',                                 !ACC
     *    'M_DD1_SU','M_DD1_DU','N_DD1_1 ',                      !DD1  
     *    'M_DS1_SU','M_DS1_DU','N_DS1_1 ',                      !DS1
#if defined TRACERS_AMP_M1 || defined TRACERS_AMP_M2 || defined TRACERS_AMP_M3 \
 || defined TRACERS_AMP_M4
     *    'M_DD2_SU','M_DD2_DU','N_DD2_1 ',                      !DD2
     *    'M_DS2_SU','M_DS2_DU','N_DS2_1 ',                      !DS2
#endif
#if defined TRACERS_AMP_M1 || defined TRACERS_AMP_M2 || defined TRACERS_AMP_M3 \
 || defined TRACERS_AMP_M5 || defined TRACERS_AMP_M6 || defined TRACERS_AMP_M7
     *    'M_SSA_SU','M_SSA_SS',                                 !SSA
     *    'M_SSC_SS',                                            !SSC  
#endif
#if defined TRACERS_AMP_M4 || defined TRACERS_AMP_M8
     *    'M_SSS_SU','M_SSS_SS',                                 !SSS
#endif
     *    'M_OCC_SU','M_OCC_OC','N_OCC_1 ',                      !OCC
     *    'M_BC1_SU','M_BC1_BC','N_BC1_1 ',                      !BC1
     *    'M_BC2_SU','M_BC2_BC','N_BC2_1 ',                      !BC2
#if defined TRACERS_AMP_M1 || defined TRACERS_AMP_M2 || defined TRACERS_AMP_M3 \
 || defined TRACERS_AMP_M5 || defined TRACERS_AMP_M6 || defined TRACERS_AMP_M7
#endif
#if defined TRACERS_AMP_M1 || defined TRACERS_AMP_M5
     *    'M_BC3_SU','M_BC3_BC','N_BC3_1 ',                      !BC3
#endif
#if defined TRACERS_AMP_M2 || defined TRACERS_AMP_M6
     *    'M_OCS_SU','M_OCS_OC','N_OCS_1 ',                      !OCS
#endif
#if defined TRACERS_AMP_M1 || defined TRACERS_AMP_M2 || defined TRACERS_AMP_M6
     *    'M_DBC_SU','M_DBC_BC','M_DBC_DU','N_DBC_1 ',           !DBC
#endif
#if defined TRACERS_AMP_M1 || defined TRACERS_AMP_M2 || defined TRACERS_AMP_M3 \
 || defined TRACERS_AMP_M6 || defined TRACERS_AMP_M7
     *    'M_BOC_SU','M_BOC_BC','M_BOC_OC','N_BOC_1 ',           !BOC
#endif
#if defined TRACERS_AMP_M1 || defined TRACERS_AMP_M2 || defined TRACERS_AMP_M5 \
 || defined TRACERS_AMP_M6
     *    'M_BCS_SU','M_BCS_BC','N_BCS_1 ',                      !BCS
#endif
     *    'M_MXX_SU',                                            !MXX
     *    'M_MXX_BC','M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ',
     *    'H2SO4   ','DMS     ','SO2     ','H2O2_s  ','NH3     ',
#endif  /* TRACERS_AMP */
     *     'Dummyspc'/) ! This line should always be last!

      character*8, parameter :: trname(ntm)=(/trname_p1(1:ntm)/)

#ifdef TRACERS_AMP
#ifdef TRACERS_AMP_M1
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,5 ,5,  !DS1,DD2
     *    5 ,6 ,6 ,6 ,7,  !DD2,DS2,SSA
     *    7 ,8,           !SSA,SSC
     *    9 ,9 ,9 ,10,10, !OCC,BC1
     *    10,11,11,11,12, !BC1,BC2,BC3
     *    12,12,13,13,13, !BC3,DBC
     *    13,14,14,14,14, !DBC,BOC
     *    15,15,15,16,16, !BCS,MXX
     *    16,16,16,16/)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,0 ,0,  !DS1,DD2
     *    5 ,0 ,0 ,6 ,0,  !DD2,DS2,SSA
     *    0 ,0,           !SSA,SSC
     *    0 ,0 ,9 ,0 ,0 , !OCC,BC1
     *    10,0 ,0 ,11,0 , !BC1,BC2,BC3
     *    0 ,12,0 ,0 ,0 , !BC3,DBC
     *    13,0 ,0 ,0 ,14, !DBC,BOC
     *    0 ,0 ,15,0 ,0 , !BCS,MXX
     *    0 ,0 , 0,16/)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15,16,17,18,19,20,
     *    21,      24,   26,27,28,29,30,        
     *    31,32,33,34,35,36,37,38,39,40,
     *    41,42,43,44,45,46,47,48,49,50,
     *    51,52,53,54   /)

      integer, parameter :: AMP_trm_nm1(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !AKK
     *    6 ,6 ,8 ,8 ,8,  !ACC,DD1
     *    11,11,11,14,14, !DS1,DD2
     *    14,17,17,17,20, !DD2,DS2,SSA
     *    20,22,          !SSA,SSC
     *    23,23,23,26,26, !OCC,BC1
     *    26,29,29,29,32, !BC1,BC2,BC3
     *    32,32,35,35,35, !BC3,DBC
     *    35,39,39,39,39, !DBC,BOC
     *    43,43,43,46,46, !BCS,MXX
     *    46,46,46,46/)
      integer, parameter :: AMP_trm_nm2(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !AKK
     *    6 ,6 ,9 ,9 ,9,  !ACC,DD1
     *    12,12,12,15,15, !DS1,DD2
     *    15,18,18,18,21, !DD2,DS2,SSA
     *    21,22,          !SSA,SSC
     *    24,24,24,27,27, !OCC,BC1
     *    27,30,30,30,33, !BC1,BC2,BC3
     *    33,33,37,37,37, !BC3,DBC
     *    37,41,41,41,41, !DBC,BOC
     *    44,44,44,50,50, !BCS,MXX
     *    50,50,50,50/)
#endif
#ifdef TRACERS_AMP_M2
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,5 ,5,  !DS1,DD2
     *    5 ,6 ,6 ,6 ,7,  !DD2,DS2,SSA
     *    7 ,8,           !SSA,SSC
     *    9 ,9 ,9 ,10,10, !OCC,BC1
     *    10,11,11,11,12, !BC1,BC2,OSC
     *    12,12,13,13,13, !BC3,DBC
     *    13,14,14,14,14, !DBC,BOC
     *    15,15,15,16,16, !BCS,MXX
     *    16,16,16,16/)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,0 ,0,  !DS1,DD2
     *    5 ,0 ,0 ,6 ,0,  !DD2,DS2,SSA
     *    0 ,0,           !SSA,SSC
     *    0 ,0 ,9 ,0 , 0, !OCC,BC1
     *    10,0 ,0 ,11,0 , !BC1,BC2,OSC
     *    0 ,12,0 ,0 ,0 , !DBC
     *    13,0 ,0 ,0 ,14, !DBC,BOC
     *    0 ,0 ,15,0 ,0 , !BCS,MXX
     *    0 ,0 ,0 ,16/)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15,16,17,18,19,20,
     *    21,      24,  26,27,28,29,30,        
     *    31,32,33,34,35,36,37,38,39,40,
     *    41,42,43,44,45,46,47,48,49,50,
     *    51,52,53,54   /)

      integer, parameter :: AMP_trm_nm1(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !AKK
     *    6 ,6 ,8 ,8 ,8,  !ACC,DD1
     *    11,11,11,14,14, !DS1,DD2
     *    14,17,17,17,20, !DD2,DS2,SSA
     *    20,22,          !SSA,SSC
     *    23,23,23,26,26, !OCC,BC1
     *    26,29,29,29,32, !BC1,BC2
     *    32,32,35,35,35, !OSC,DBC
     *    35,39,39,39,39, !DBC,BOC
     *    43,43,43,46,46, !BCS,MXX
     *    46,46,46,46/)
      integer, parameter :: AMP_trm_nm2(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !AKK
     *    6 ,6 ,9 ,9 ,9,  !ACC,DD1
     *    12,12,12,15,15, !DS1,DD2
     *    15,18,18,18,21, !DD2,DS2,SSA
     *    21,22,          !SSA,SSC
     *    24,24,24,27,27, !OCC,BC1
     *    27,30,30,30,33, !BC1,BC2
     *    33,33,37,37,37, !OCS,DBC
     *    37,41,41,41,41, !DBC,BOC
     *    44,44,44,50,50, !BCS,MXX
     *    50,50,50,50/)
#endif
#ifdef TRACERS_AMP_M3
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,5 ,5,  !DS1,DD2
     *    5 ,6 ,6 ,6 ,7,  !DD2,DS2,SSA
     *    7 ,8,           !SSA,SSC
     *    9 ,9 ,9 ,10,10, !OCC,BC1
     *    10,11,11,11,12, !BC1,BC2,BOc
     *    12,12,12,13,13, !BOC,MXX
     *    13,13,13,13/)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,0 ,0,  !DS1,DD2
     *    5 ,0 ,0 ,6 ,0,  !DD2,DS2,SSA
     *    0 ,0,           !SSA,SSC
     *    0 ,0 ,9 ,0 ,0 , !OCC,BC1
     *    10,0 ,0 ,11,0 , !BC1,BC2,BOc
     *    0 ,0 ,12,0 ,0 , !BOC,MXX
     *    0 ,0 ,0 ,13/)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15      ,18   ,20,
     *    21,22,23,24,25,26,27,28,29,30,        
     *    31,32,33,34,35,36,37,38,39,40,
     *    41,42,43,44  /)

      integer, parameter :: AMP_trm_nm1(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !AKK
     *    6 ,6 ,8 ,8 ,8,  !ACC,DD1
     *    11,11,11,14,14, !DS1,DD2
     *    14,17,17,17,20, !DD2,DS2,SSA
     *    20,22,          !SSA,SSC
     *    23,23,23,26,26, !OCC,BC1
     *    26,29,29,29,32, !BC1,BC2,BOC
     *    32,32,32,36,36, !BOC,MXX
     *    36,36,36,36/)

      integer, parameter :: AMP_trm_nm2(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !AKK
     *    6 ,6 ,9 ,9 ,9,  !ACC,DD1
     *    12,12,12,15,15, !DS1,DD2
     *    15,18,18,18,21, !DD2,DS2,SSA
     *    21,22,          !SSA,SSC
     *    24,24,24,27,27, !OCC,BC1
     *    27,30,30,30,34, !BC1,BC2,BOC
     *    34,34,34,40,40, !BOC,MXX
     *    40,40,40,40/)
 
#endif
#ifdef TRACERS_AMP_M4
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !ACC
     *    2 ,2 ,2 ,3 ,3,  !DD1,DS1
     *    3 ,4 ,4 ,4 ,5,  !DS1,DD2,DS2
     *    5 ,5 ,6 ,6,     !DS2,SSS
     *    7 ,7 ,7 ,8 ,8,  !OCC,BC1
     *    8 ,9 ,9 ,9 ,10, !BC1,BC2,MXX
     *    10,10,10,10,10  !MXX
     *    /)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !ACC
     *    0 ,0 ,2 ,0 ,0,  !DD1,DS1
     *    3 ,0 ,0 ,4 ,0,  !DS1,DD2,DS2
     *    0 ,5 ,0 ,0,     !DS2,SSS
     *    0 ,0 ,6 ,0 ,0,  !OCC,BC1
     *    7 ,0 ,0 ,8 ,0,  !BC1,BC2,MXX
     *    0, 0, 0 ,0 ,9  !MXX
     *    /)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15,16,17,18,19,
     *    21,22,23,24,25,26,27,28,29,30,        
     *    31,32,33,34,35   /)
      integer, parameter :: AMP_trm_nm1(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !ACC
     *    6 ,6 ,6 ,9 ,9,  !DD1
     *    9 ,12,12,12,15, !DS1,DD2
     *    15,15,18,18,    !DD2,DS2,SSA
     *    20,20,20,23,23, !SSA,SSC
     *    23,26,26,26,29, !OCC,BC1
     *    29,29,29,29,29/)

      integer, parameter :: AMP_trm_nm2(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !ACC
     *    7 ,7 ,7 ,10,10,  !DD1
     *    10,13,13,13,16, !DS1,DD2
     *    16,16,19,19,    !DD2,DS2,SSA
     *    21,21,21,24,24, !SSA,SSC
     *    24,27,27,27,33, !OCC,BC1
     *    33,33,33,33,33/)
 
#endif
#ifdef TRACERS_AMP_M5
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,       !DS1
     *    5 ,             !SSA
     *    5 ,6 ,          !SSA,SSC
     *    7 ,7 ,7, 8, 8,  !OCC,BC1
     *    8, 9, 9, 9,10,  !BC1,BC2,BC3
     *    10,10,11,11,11, !BC3,DBC
     *    11,12,12,12,12, !DBC,BOC
     *    13,13,13,14,14, !BCS,MXX
     *    14,14,14,14/)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,       !DS1
     *    0 ,             !SSA
     *    0 ,0 ,          !SSA,SSC
     *    0 ,0 ,7, 0, 0,  !OCC,BC1
     *    8, 0, 0, 9, 0,  !BC1,BC2,BC3
     *    0, 10,0 ,0 ,0, !BC3,DBC
     *    11,0 ,0 ,0 ,12, !DBC,BOC
     *    0, 0, 13,0 ,0, !BCS,MXX
     *    0,0,0,14/)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15,      18,   20,
     *    21,22,23,24,25,26,27,28,29,30,        
     *    31,32,33,34,35,36,37,38,39,40,
     *    41,42,43,44,45,46,47,48   /)
#endif
#ifdef TRACERS_AMP_M6
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,       !DS1
     *    5 ,             !SSA
     *    5 ,6 ,          !SSA,SSC
     *    7 ,7 ,7, 8, 8,  !OCC,BC1
     *    8, 9, 9, 9,10,  !BC1,BC2,OCS
     *    10,10,11,11,11, !OCS,DBC
     *    11,12,12,12,12, !DBC,BOC
     *    13,13,13,14,14, !BCS,MXX
     *    14,14,14,14/)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,       !DS1
     *    0 ,             !SSA
     *    0 ,0 ,          !SSA,SSC
     *    0 ,0 ,7, 0, 0,  !OCC,BC1
     *    8, 0, 0, 9, 0,  !BC1,BC2,OCS
     *    0 ,10,0 ,0 ,0 , !OCS,DBC
     *    11,0 ,0 ,0 ,12, !DBC,BOC
     *    0 ,0 ,13,0 ,0 , !BCS,MXX
     *    0 ,0 ,0 ,14/)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15,      18,   20,
     *    21,22,23,24,25,26,27,28,29,30,        
     *    31,32,33,34,35,36,37,38,39,40,
     *    41,42,43,44,45,46,47,48   /)
#endif
#ifdef TRACERS_AMP_M7
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,       !DS1
     *    5 ,             !SSA
     *    5 ,6 ,          !SSA,SSC
     *    7 ,7 ,7, 8, 8,  !OCC,BC1
     *    8, 9, 9, 9,     !BC1,BC2
     *    10,10,10,10,    !OCS,DBC
     *    11,11,          !MXX
     *    11,11,11,11/)   !MXX
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,       !DS1
     *    0 ,             !SSA
     *    0 ,0 ,          !SSA,SSC
     *    0 ,0 ,7, 0, 0,  !OCC,BC1
     *    8, 0, 0, 9,     !BC1,BC2
     *    0, 0, 0,10,    !OCS,DBC
     *    0, 0,          !MXX
     *    0, 0, 0,11/)   !MXX
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15,      18,   20,
     *    21,22,23,24,25,26,27,28,29,30,        
     *    31,32,33,34,35,36,37,38/)
#endif
#ifdef TRACERS_AMP_M8
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,       !
     *    1 ,1 ,2 ,2 ,2,  !ACC,DD1
     *    3 ,3 ,3 ,       !DS1
     *    4 ,             !SSS
     *    4 ,5 ,          !SSS
     *    5 ,5 ,5, 6, 6,  !OCC,BC1
     *    6, 7, 7, 7, 8,  !BC1,BC2
     *    8,8,8,8,8,/)    !MXX
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,       !
     *    0 ,1 ,0 ,0 ,2,  !ACC,DD1
     *    0 ,0 ,3 ,       !DS1
     *    0 ,             !SSS
     *    0 ,0 ,          !SSS
     *    0 ,0 ,5, 0, 0,  !OCC,BC1
     *    6, 0, 0, 7, 0,  !BC1,BC2
     *    0,0,0,0,8,/)    !MXX
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,        
     *    6 ,7 ,8 ,9 ,10,        
     *    11,12,13,15,16,        
     *    17,18,19,20,21,
     *    22,23,24,25,26,        
     *    27,28,29   /)
#endif
#endif   /* TRACERS_AMP */

!@var N_XXX: variable names of indices for tracers (init = 0)
      integer ::
     *                 n_SF6_c=0,                                        
     *     n_Air=0,    n_SF6=0,   n_Rn222=0, n_CO2=0,      n_N2O=0,
     *     n_CFC11=0,  n_14CO2=0, n_CH4=0,   n_O3=0,       n_water=0,
     *     n_H2O18=0,  n_HDO=0,   n_HTO=0,   n_Ox=0,       n_NOx=0,
     *     n_N2O5=0,   n_HNO3=0,  n_H2O2=0,  n_CH3OOH=0,   n_HCHO=0,
     *     n_HO2NO2=0, n_CO=0,    n_PAN=0,   n_H2O17=0,
     *     n_Isoprene=0, n_AlkylNit=0, n_Alkenes=0, n_Paraffin=0,
     *     n_stratOx=0, n_Terpenes=0,n_codirect=0,
     *     n_isopp1g=0,n_isopp1a=0,n_isopp2g=0,n_isopp2a=0,
     *     n_apinp1g=0,n_apinp1a=0,n_apinp2g=0,n_apinp2a=0,
     *     n_DMS=0,    n_MSA=0,   n_SO2=0,   n_SO4=0,    n_H2O2_s=0,
     *     n_ClOx=0,   n_BrOx=0,  n_HCl=0,   n_HOCl=0,   n_ClONO2=0,
     *     n_HBr=0,    n_HOBr=0,  n_BrONO2=0,n_CFC=0,    n_GLT=0,
     *     n_Pb210 = 0,n_Be7=0,   n_Be10=0,
     .     n_CFCn=0,   n_CO2n=0,  n_Age=0,
     *     n_seasalt1=0,  n_seasalt2=0, n_SO4_d1=0,  n_SO4_d2=0,
     *     n_SO4_d3=0,n_N_d1=0,  n_N_d2=0,  n_N_d3=0,
     &     n_NH3=0,   n_NH4=0,   n_NO3p=0,
     *     n_BCII=0,  n_BCIA=0,  n_BCB=0,
     *     n_OCII=0,  n_OCIA=0,  n_OCB=0,
     *     n_OCI1=0,  n_OCI2=0,  n_OCI3=0,
     *     n_OCA1=0,  n_OCA2=0,  n_OCA3=0, n_OCA4,
     &     n_clay=0,   n_silt1=0, n_silt2=0, n_silt3=0, n_silt4=0,
     &     n_clayilli=0,n_claykaol=0,n_claysmec=0,n_claycalc=0,
     &     n_clayquar=0,n_sil1quar=0,n_sil1feld=0,n_sil1calc=0,
     &     n_sil1hema=0,n_sil1gyps=0,n_sil2quar=0,n_sil2feld=0,
     &     n_sil2calc=0,n_sil2hema=0,n_sil2gyps=0,n_sil3quar=0,
     &     n_sil3feld=0,n_sil3calc=0,n_sil3hema=0,n_sil3gyps=0,
     &     n_sil1quhe=0,n_sil2quhe=0,n_sil3quhe=0,
     *     n_M_NO3=0,   n_M_NH4=0,   n_M_H2O=0,   n_M_AKK_SU=0,
     *     n_N_AKK_1=0, n_M_ACC_SU=0,n_N_ACC_1=0, n_M_DD1_SU=0,
     *     n_M_DD1_DU=0,n_N_DD1_1=0, n_M_DS1_SU=0,n_M_DS1_DU=0,
     *     n_N_DS1_1 =0,n_M_DD2_SU=0,n_M_DD2_DU=0,n_N_DD2_1 =0,
     *     n_M_DS2_SU=0,n_M_DS2_DU=0,n_N_DS2_1 =0,n_M_SSA_SU=0,
     *     n_M_SSA_SS=0,n_M_SSC_SS=0,
     *     n_M_OCC_SU=0,n_M_OCC_OC=0,n_N_OCC_1 =0,
     *     n_M_BC1_SU=0,n_M_BC1_BC=0,n_N_BC1_1 =0,n_M_BC2_SU=0,
     *     n_M_BC2_BC=0,n_N_BC2_1 =0,n_M_BC3_SU=0,n_M_BC3_BC=0,
     *     n_N_BC3_1 =0,n_M_DBC_SU=0,n_M_DBC_BC=0,n_M_DBC_DU=0,
     *     n_N_DBC_1 =0,n_M_BOC_SU=0,n_M_BOC_BC=0,n_M_BOC_OC=0,
     *     n_N_BOC_1=0, n_M_BCS_SU=0,n_M_BCS_BC=0,n_N_BCS_1 =0,
     *     n_M_MXX_SU=0,n_M_MXX_BC=0,n_M_MXX_OC=0,n_M_MXX_DU=0,
     *     n_M_MXX_SS=0,n_N_MXX_1 =0,n_M_OCS_SU=0,n_M_OCS_OC=0,
     *     n_N_OCS_1=0,n_M_SSS_SS=0,n_M_SSS_SU=0,
     *     n_H2SO4=0, n_N_SSA_1=0, n_N_SSC_1=0

C**** standard tracer and tracer moment arrays

!@var TRM: Tracer array (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: trm

!@var TRMOM: Second order moments for tracers (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: trmom

!@var TRDN1: lowest level downdraft tracer concentration (kg/kg)
       REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: trdn1

#ifdef TRACERS_WATER
!@var TRWM tracer in cloud liquid water amount (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: trwm
#endif

!@var daily_z: altitude of model levels (m), updated once per day
!@+   and used for vertical distribution of 3D emissions
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: daily_z

C**** The following general arrays are set in tracer_IC

!@var T_QLIMIT: if t_qlimit=.true. tracer is maintained as positive
      logical, dimension(ntm) :: t_qlimit
!@var trdecay radioactive decay constant (1/s) (=0 for stable tracers)
      real*8, dimension(ntm) :: trdecay
!@dbparam ITIME_TR0: start time for each tracer (hours)
      integer, dimension(ntm) :: itime_tr0
#ifdef TRACERS_ON
!@var NTM_POWER: Power of 10 associated with each tracer (for printing)
      integer, dimension(ntm) :: ntm_power
!@var TR_MM: molecular mass of each tracer (g/mole)
      real*8, dimension(ntm) :: tr_mm
!@var mass2vol: mass to volume ratio = mair/tr_mm
      real*8, dimension(ntm) :: mass2vol
!@var vol2mass: volume to mass ratio = tr_mm/mair
      real*8, dimension(ntm) :: vol2mass
!@var needtrs: true if surface tracer value from PBL is required
      logical, dimension(ntm) :: needtrs
#ifdef TRACERS_DRYDEP
!@var dodrydep: true if tracer should undergo dry deposition
      logical, dimension(ntm) :: dodrydep
!@var F0 reactivity factor for oxidation of biological substances
      real*8, dimension(ntm)  :: F0
!@var HSTAR Henry's Law const for tracer dry deposition. mole/(L atm)
!@+   Same as the tr_RKD wet dep variable, except for units.
!@+   If F0 & HSTAR both 0, & tracer not particulate, then no drydep.
      real*8, dimension(ntm)  :: HSTAR
#endif
#endif
!@var do_fire: true if tracer should have emissions via flammability
      logical, dimension(ntm) :: do_fire
!@var nBBsources number of sources attributed to biomass burning
      integer, dimension(ntm) :: nBBsources
!@var emisPerFireByVegType tracer emissions per fire count as a
!@+ function of 12 standard GISS (VDATA) vegetation types, or Ent
!@+ remapped to them. Thus the "12" hardcode.
      real*8, dimension(ntm,12) :: emisPerFireByVegType=0.d0
!@var trpdens tracer particle density (kg/m^3)
!@+               (=0 for non-particle tracers)
      real*8, dimension(ntm) :: trpdens
!@var trradius tracer effective radius (m) (=0 for non particle tracers)
      real*8, dimension(ntm) :: trradius

!@var tr_wd_TYPE: tracer wet dep type (gas, particle, water)
      integer, dimension(ntm) :: tr_wd_TYPE
!@var tr_RKD: Henry's Law coefficient (in mole/Joule please !)
      real*8, dimension(ntm) :: tr_RKD
!@var tr_DHD: coefficient of temperature-dependence term of Henry's
!@+   Law coefficient (in Joule/mole please !)
      real*8, dimension(ntm) :: tr_DHD
!@var fq_aer fraction of aerosol that condenses
      real*8 fq_aer(ntm)
!@var rc_washt aerosol washout rate
      REAL*8 :: rc_washt(ntm)=1.D-1

!@var isdust: index array for testing if tracer is a dust type
      integer, dimension(ntm) :: isdust

C**** parameters for tr_wd_TYPE
!@param nGAS   index for wetdep tracer type = gas
!@param nPART  index for wetdep tracer type = particle/aerosol
!@param nWATER index for wetdep tracer type = water
      integer, parameter :: nGAS=1, nPART=2, nWATER=3

#ifdef TRACERS_WATER
!@param nWD_TYPES number of tracer types for wetdep purposes
      integer, parameter :: nWD_TYPES=3 !(gas,particle,water)
!@param tr_evap_fact fraction of re-evaporation of raindrops by tracer type
C note, tr_evap_fact is not dimensioned as NTM:
      REAL*8, parameter, dimension(nWD_TYPES) :: tr_evap_fact=
     *     (/1.d0, 0.5d0,  1.d0/)
!@var tr_H2ObyCH4 conc. of tracer in water from methane oxidation
      real*8, dimension(ntm) :: tr_H2ObyCH4
!@var dowetdep true if tracer has some form of wet deposition
      logical, dimension(ntm) :: dowetdep
#endif

#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
!@var TRW0 default tracer concentration in water (kg/kg)
      real*8, dimension(ntm) :: trw0
!@var NTROCN scaling power factor for ocean/ice tracer concentrations
      integer, dimension(ntm) :: ntrocn
!@var conc_from_fw true if ocean conc is defined using fresh water
      logical, dimension(ntm) :: conc_from_fw = .true.
#endif

#ifdef TRACERS_OCEAN
!@var TRGLAC tracer ratio in glacial runoff to ocean (kg/kg)
      real*8, dimension(ntm) :: trglac
#endif

C**** Diagnostic indices and meta-data

!@var ntsurfsrcmax maximum number of surface 2D sources/sinks
      integer, parameter :: ntsurfsrcmax=16
!@var ntsurfsrc no. of non-interactive surface sources for each tracer
      integer, dimension(ntm) :: ntsurfsrc
!@var ntisurfsrc no. of interactive surface sources for each tracer
      integer, dimension(ntm) :: ntisurfsrc
!@var nt3Dsrcmax maximum number of 3D tracer sources/sinks
      integer, parameter :: nt3Dsrcmax=7
!@var sfc_src array holds tracer sources that go into trsource( )
!@+ maybe wasteful of memory, but works for now...
      real*8, allocatable, dimension(:,:,:,:) :: sfc_src

!@var MTRACE: timing index for tracers
!@var MCHEM: timing index for chemistry (if needed)
      integer mtrace, mchem
!@var MTRADV: timing index for tracer advection (if requested)
#ifdef TRAC_ADV_CPU
      integer mtradv
#endif
!@dbparam no_emis_over_ice, switch whether (>0) or not (<=0) to
!@+ set surface tracer emissions to zero over >90% ice boxes
      integer :: no_emis_over_ice=0

C**** EVERYTHING BELOW HERE IS TRACER SPECIFIC. PLEASE THINK 
C**** ABOUT MOVING IT ELSEWHERE

!@var 3D on-line radical array for interactive aerosol and gas
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oh_live
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: no3_live

#ifdef TRACERS_SPECIAL_O18
C**** Water isotope specific parameters

#ifdef TRACERS_OCEAN
      integer :: water_tracer_ic=1 ! Read water tracers ic from H2O18ic (=1) or set all to SMOW (=0)
#endif

!@dbparam supsatfac factor controlling super saturation for isotopes
      real*8 :: supsatfac = 2d-3
!@var iso_index indexing taking actual tracer number to isotope
!@+   fractionation number (1=water,2=h2o18,3=hdo,4=hto,5=h2o17)
      integer :: iso_index(ntm)
#endif

#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_AMP)
C**** Chemistry specific 3-D arrays

!@var RSULF1, RSULF2, RSULF3, RSULF4: rate coefficients
c for gas phase sulfur chemistry used by aerosol and chemistry models
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)::rsulf1,rsulf2,rsulf3,rsulf4
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP)
C**** Aerosol specific switches and arrays

!@dbparam imAER 0 determines emission choice: 0 present-day, 1 AEROCOM,
c  2 Bond/Streets past, present, future using sector inputs; 3 historic
      integer :: imAER = 1
!!@dbparam OFFLINE_DMS_SS is 0 for standard case, 1 for offline dms, seasalt emission
      integer :: OFFLINE_DMS_SS = 0
!!@dbparam OFFLINE_SS is 0 for standard case, 1 for offline seasalt emission
      integer :: OFFLINE_SS = 0
!@dbparam imPI is 0 for industrial simulations, 1 for pre-industrial
      integer :: imPI = 0
!@dbparam aer_int_yr indicates year of emission
      integer :: aer_int_yr = 0
!@var SNFST0,TNFST0 are instantaneous SW, LW aerosol forcings for AEROCOM
c      real*8 SNFST0(2,NTM,IM,JM),TNFST0(2,NTM,IM,JM)
#endif

C**** tracer specific switches

!@dbparam rnsrc is a switch for Radon sources 
!@+       0=standard, 1=Conen&Robertson, 2=modified Schery&Wasiolek
      integer :: rnsrc = 0

C**** arrays that could be general, but are only used by chemistry

!@var ssname holds source name, read from file header, e.g. to be
!@+ placed into lname and sname arrays.
!@var freq frequency (annual? monthly?) read from emis file header
!@var res horiz. resolution (M? F?) read from emis file header
!@var nameT tracer name read from emis file header (should match trname)
!@var ty_start starting year/decade for a transient emissions file
!@var ty_end     ending year/decade for a transient emissions file
!! dbparam trans_emis_overr_yr year for overriding tracer transient emis
!! dbparam trans_emis_overr_day day for overriding tracer transient emis
!@var trans_emis_overr_yr year for overriding tracer transient emis
!@var trans_emis_overr_day day for overriding tracer transient emis
      character*30, dimension(ntm,ntsurfsrcmax) :: ssname ! some maybe
      character*10, dimension(ntm,ntsurfsrcmax) :: nameT  ! need not be
      character*1, dimension(ntm,ntsurfsrcmax) :: freq,res ! arrays
      character*9, dimension(ntm,ntsurfsrcmax) :: Tyears   ! here...
      integer, dimension(ntm,ntsurfsrcmax) :: ty_start,ty_end,delTyr
      integer :: trans_emis_overr_yr=0, trans_emis_overr_day=0
! ---- section for altering tracers sources by sector/region ----
!@param n_max_sect maximum number of sectors for emissions altering
!@param n_max_reg  maximum number of regions for emissions altering
      integer, parameter :: n_max_sect=10, n_max_reg=10
!@var num_tr_sectors number of sectors for a particular tracer and source
      integer, dimension(ntm,ntsurfsrcmax) :: num_tr_sectors
!@var num_tr_sectors3D number of sectors for a tracer's 3D source
      integer, dimension(ntm,nt3Dsrcmax) :: num_tr_sectors3D
!@var num_regions the number of source-altering regions from rundeck
!@var num_sectors the number of source-altering sectors from rundeck
      integer :: num_regions, num_sectors
!@var alter_sources true if any source altering factors are on
      logical :: alter_sources
!@var reg_N the north edge of rectangular regions for emissions altering
!@var reg_S the south edge of rectangular regions for emissions altering
!@var reg_E the east  edge of rectangular regions for emissions altering
!@var reg_W the west  edge of rectangular regions for emissions altering
      real*8, dimension(n_max_reg) :: reg_N,reg_S,reg_E,reg_W
!@var tr_sect_index array hold the sector index for given tracer/source
      integer, dimension(ntm,ntsurfsrcmax,n_max_sect) :: tr_sect_index
!@var tr_sect_index3D holds 3d source sector index for given tracer/source
      integer, dimension(ntm,nt3Dsrcmax,n_max_sect) :: tr_sect_index3D
!@var tr_sect_name array hold the sector name for given tracer/source
      character*10,dimension(ntm,ntsurfsrcmax,n_max_sect):: tr_sect_name
!@var tr_sect_name3D holds 3d source sector name for given tracer/source
      character*10,dimension(ntm,nt3Dsrcmax,n_max_sect):: tr_sect_name3D
!@var sect_name array hold the sector names (all)
      character*10,dimension(n_max_sect):: sect_name
!@var ef_fact the actual factors that alter sources by region/sector
!@var ef_fact3d factors used to alter 3D sources (these are more
!@+ hard-coded for now...)
      real*8, dimension(n_max_sect,n_max_reg) :: ef_fact,ef_fact3D
! variables for outputting a map of the regions:
      real*8, allocatable, dimension(:,:) :: ef_REG_IJ
! --- end of source-altering section ----------------------------
!@param nChemistry index for tracer chemistry 3D source
!@param nOverwrite index for tracer overwrite 3D source
!@param nOther index for tracer misc. 3D source
!@param nAircraft index for tracer aircraft 3D source
!@param nBiomass index for tracer biomass burning 3D source
!@param nVolcanic index for tracer volcano 3D source
!@param nChemloss index for tracer chemistry 3D loss
! Must be a better way to do this, but for now, it is better than
! hardcoding indicies with a number like "3":
      INTEGER, PARAMETER :: nChemistry = 1, nOverwrite = 2,
     &     nOther = 3, nAircraft = 4, nBiomass = 5,
     &     nVolcanic = 6, nChemloss = 7

#ifdef TRACERS_GASEXCH_ocean
#ifdef TRACERS_GASEXCH_ocean_CFC
!@var ocmip_cfc: CFC-11 emissions estimated from OCMIP surf.conc.
      !60years (1939--1998) OCMIP surfc. concentr. converted to
      !global averaged emission rates
      !each value corresponds to the annual value
      REAL*8, DIMENSION(67,NTM) :: ocmip_cfc
#endif
#endif

#if (defined TRACERS_HETCHEM) || (defined TRACERS_NITRATE)
      integer, parameter :: rhet=3
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: rxts,rxts1,rxts2,rxts3
     *                                         ,rxts4
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: krate
#endif

!@var xyz_count,xyz_list count/list of tracers in category xyz.
!@+   A tracer can belong to more than one list.
c note: not applying CPP when declaring counts/lists.
      integer ::  ! counts of tracers meeting the following criteria:
     &      active_count  ! itime_tr0 <= itime
     &     ,gases_count   ! tr_wd_type == nGas
     &     ,aero_count    ! tr_wd_type == nPart
     &     ,water_count   ! tr_wd_type == nWater
     &     ,hlawt_count   ! tr_wd_type == nGas and tr_DHD != 0
     &     ,aqchem_count  ! participates in cloud aqueous chemistry
      integer, dimension(:), allocatable ::
     &     active_list,gases_list,aero_list,water_list,
     &     hlawt_list,aqchem_list

      contains

      subroutine remake_tracer_lists()
!@sum regenerates the counts and lists of tracers in various categories
      use model_com, only : itime,coupled_chem
      implicit none
      integer :: n,nactive
      integer, dimension(1000) ::
     &     tmplist_active,tmplist_gases,tmplist_aero,tmplist_water,
     &     tmplist_hlawt,tmplist_aqchem
      active_count = 0
      gases_count = 0
      aero_count = 0
      water_count = 0
      hlawt_count = 0
      aqchem_count = 0
      do n=1,ntm

        if(itime.lt.itime_tr0(n)) cycle

        active_count = active_count + 1
        tmplist_active(active_count) = n

        select case(tr_wd_type(n))
        case(nGAS)
          gases_count = gases_count + 1
          tmplist_gases(gases_count) = active_count
          if(tr_DHD(n).ne.0.) then
            hlawt_count = hlawt_count + 1
            tmplist_hlawt(hlawt_count) = active_count
          endif
        case(nPART)
          aero_count = aero_count + 1
          tmplist_aero(aero_count) = active_count
        case(nWATER)
          water_count = water_count + 1
          tmplist_water(water_count) = active_count
        end select
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
        select case (trname(n))
        case('SO2','SO4','H2O2_s','H2O2')
          if ( .not.(
     *         (trname(n).eq."H2O2" .and. coupled_chem.eq.0).or.
     *         (trname(n).eq."H2O2_s" .and. coupled_chem.eq.1)) )
     *         then
            aqchem_count = aqchem_count + 1
            tmplist_aqchem(aqchem_count) = active_count
          end if
        end select
#endif
      enddo
      if(allocated(active_list)) deallocate(active_list)
      allocate(active_list(active_count))
      active_list = tmplist_active(1:active_count)

      if(allocated(gases_list)) deallocate(gases_list)
      allocate(gases_list(gases_count))
      gases_list = tmplist_gases(1:gases_count)

      if(allocated(aero_list )) deallocate(aero_list )
      allocate(aero_list(aero_count))
      aero_list  = tmplist_aero(1:aero_count)

      if(allocated(water_list)) deallocate(water_list)
      allocate(water_list(water_count))
      water_list = tmplist_water(1:water_count)

      if(allocated(hlawt_list)) deallocate(hlawt_list)
      allocate(hlawt_list(hlawt_count))
      hlawt_list = tmplist_hlawt(1:hlawt_count)

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      if(allocated(aqchem_list)) deallocate(aqchem_list)
      allocate(aqchem_list(aqchem_count))
      aqchem_list = tmplist_aqchem(1:aqchem_count)
#endif

      return
      end subroutine remake_tracer_lists

      END MODULE TRACER_COM

      SUBROUTINE ALLOC_TRACER_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE TRACER_COM
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, GET
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H, I_1H, I_0H

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO

      ALLOCATE(   ef_REG_IJ(I_0H:I_1H,J_0H:J_1H) )
      ALLOCATE(     oh_live(I_0H:I_1H,J_0H:J_1H,LM),
     *             no3_live(I_0H:I_1H,J_0H:J_1H,LM),
     *                  trm(I_0H:I_1H,J_0H:J_1H,LM,NTM),
     *                trmom(NMOM,I_0H:I_1H,J_0H:J_1H,LM,NTM),
     *                trdn1(NTM,I_0H:I_1H,J_0H:J_1H),
     *              sfc_src(I_0H:I_1H,J_0H:J_1H,ntm,ntsurfsrcmax))

      ALLOCATE(  daily_z(I_0H:I_1H,J_0H:J_1H,LM) )
      daily_z = 0.

#ifdef TRACERS_WATER
      ALLOCATE(        trwm(I_0H:I_1H,J_0H:J_1H,LM,NTM) )
#endif
#ifdef TRACERS_HETCHEM
      ALLOCATE( rxts(I_0H:I_1H,J_0H:J_1H,LM),
     *          rxts1(I_0H:I_1H,J_0H:J_1H,LM),
     *          rxts2(I_0H:I_1H,J_0H:J_1H,LM),
     *          rxts3(I_0H:I_1H,J_0H:J_1H,LM),
     *          rxts4(I_0H:I_1H,J_0H:J_1H,LM),
     *          krate(I_0H:I_1H,J_0H:J_1H,LM,8,rhet))
#endif
#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_AMP)
      ALLOCATE(  rsulf1(I_0H:I_1H,J_0H:J_1H,LM),
     *           rsulf2(I_0H:I_1H,J_0H:J_1H,LM),
     *           rsulf3(I_0H:I_1H,J_0H:J_1H,LM),
     *           rsulf4(I_0H:I_1H,J_0H:J_1H,LM) )
#endif

      END SUBROUTINE ALLOC_TRACER_COM

