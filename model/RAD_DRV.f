#include "rundeck_opts.h"
#ifdef SKIP_TRACERS_RAD
#undef TRACERS_ON
#endif
!@sum RAD_DRV contains drivers for the radiation related routines
!@ver  2009/05/11
!@cont init_RAD, RADIA
C**** semi-random cloud overlap (computed opt.d+diagn)
C**** to be used with R99E or later radiation  routines.  carbon/2
C****

      SUBROUTINE CALC_ZENITH_ANGLE
!@sum calculate zenith angle for current time step
!@auth Gavin Schmidt (from RADIA)
      USE CONSTANT, only : twopi,sday
      USE MODEL_COM, only : itime,nday,dtsrc
      USE RAD_COM, only : cosz1
      USE RAD_COSZ0, only : coszt
      IMPLICIT NONE
      INTEGER JTIME
      REAL*8 ROT1,ROT2

      JTIME=MOD(ITIME,NDAY)
      ROT1=(TWOPI*JTIME)/NDAY
      ROT2=ROT1+TWOPI*DTsrc/SDAY
      CALL COSZT (ROT1,ROT2,COSZ1)

      END SUBROUTINE CALC_ZENITH_ANGLE

      SUBROUTINE init_RAD(istart)
!@sum  init_RAD initialises radiation code
!@auth Original Development Team
!@ver  1.0
!@calls RADPAR:RCOMP1, ORBPAR
      USE FILEMANAGER
      USE PARAM
      USE CONSTANT, only : grav,bysha,twopi
      USE MODEL_COM, only : jm,lm,dtsrc,nrad
     *     ,kradia,lm_req,pednl00,jyear,iyear1,itime
      USE DOMAIN_DECOMP_ATM, only : grid, get, write_parallel, am_i_root
#ifndef CUBED_SPHERE
      USE GEOM, only : lat_dg,imaxj
#endif
      USE GEOM, only : lonlat_to_ij
      USE RADPAR, only : !rcomp1,writer,writet       ! routines
     &      PTLISO ,KTREND ,LMR=>NL, PLB, LS1_loc
     *     ,KCLDEM,KSIALB,KSOLAR, SHL, snoage_fac_max, KZSNOW
     *     ,KYEARS,KJDAYS,MADLUV, KYEARG,KJDAYG,MADGHG
     *     ,KYEARO,KJDAYO,MADO3M, KYEARA,KJDAYA,MADAER , O3YR_max
     *     ,KYEARD,KJDAYD,MADDST, KYEARV,KJDAYV,MADVOL
     *     ,KYEARE,KJDAYE,MADEPS, KYEARR,KJDAYR
!g95     *     ,FSXAER,FTXAER    ! scaling (on/off) for default aerosols
     *     ,ITR,NTRACE        ! turning on options for extra aerosols
     *     ,FS8OPX,FT8OPX, TRRDRY,KRHTRA,TRADEN,REFDRY
     *     ,rcomp1, writer, writet
     *     ,FSTASC,FTTASC
#ifdef ALTER_RADF_BY_LAT
     *     ,FS8OPX_orig,FT8OPX_orig
#endif
      USE RAD_COM, only : kliq,s0x, co2x,n2ox,ch4x,cfc11x,cfc12x,xGHGx
     *     ,s0_yr,s0_day,ghg_yr,ghg_day,volc_yr,volc_day,aero_yr,O3_yr
     *     ,H2ObyCH4,dH2O,h2ostratx,O3x,RHfix,CLDx,ref_mult
     *     ,obliq,eccn,omegt,obliq_def,eccn_def,omegt_def
     *     ,CC_cdncx,OD_cdncx,cdncl,pcdnc,vcdnc
     *     ,cloud_rad_forc,aer_rad_forc   ,ijradfrc ! only for kradia<0
     *     ,PLB0,shl0  ! saved to avoid OMP-copyin of input arrays
     *     ,albsn_yr,dALBsnX,depoBC,depoBC_1990, nradfrc
     *     ,rad_interact_aer,clim_interact_chem,rad_forc_lev,ntrix,wttr
     *     ,nrad_clay,variable_orb_par,orb_par_year_bp,orb_par
#ifdef TRACERS_ON
     &     ,nTracerRadiaActive,tracerRadiaActiveFlag
#endif
#ifdef TRACERS_SPECIAL_Shindell
     *     ,maxNtraceFastj
#endif
#ifdef ALTER_RADF_BY_LAT
     *     ,FULGAS_lat,FS8OPX_lat,FT8OPX_lat
#endif
#ifdef CHL_from_SeaWIFs
     *     ,iu_CHL
#endif
#if (defined OBIO_RAD_coupling) || (defined CHL_from_SeaWIFs)
     *     ,wfac
#endif
      use RAD_COSZ0, only : cosz_init
      USE CLOUDS_COM, only : llow
      USE DIAG_COM, only : iwrite,jwrite,itwrite,save3dAOD
#ifdef TRACERS_ON
      USE TRACER_COM
#endif
#ifdef TRACERS_AMP
      USE AERO_CONFIG, only: nmodes
#endif
      use IndirectAerParam_mod, only: dCDNC_est, aermix
      IMPLICIT NONE

      INTEGER L,LR,n1,istart,n,nn,iu2 ! LONR,LATR
      REAL*8 PLBx(LM+1),pyear,lonlat(2)
!@var NRFUN indices of unit numbers for radiation routines
      INTEGER NRFUN(14),IU,i,j,ic,jc,ij(2)
!@dbparam Ikliq 0,1,-1 initialize kliq as dry,equil,current model state
      INTEGER :: Ikliq = -1  !  get kliq-array from restart file
!@var RUNSTR names of files for radiation routines
      CHARACTER*5 :: RUNSTR(14) = (/"RADN1","RADN2","RADN3",
     *     "RADN4","RADN5","RADN6","RADN7","RADN8",
     *     "RADN9","RADNA","RADNB","RADNC","RADND",
     *     "RADNE"/)
!@var QBIN true if files for radiation input files are binary
      LOGICAL :: QBIN(14)=(/.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.
     *     ,.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE./)

#if (defined OBIO_RAD_coupling) || (defined CHL_from_SeaWIFs)
      integer, parameter :: nlt=33
      real*8 :: aw(nlt), bw(nlt), saw, sbw
      real*8 :: b0, b1, b2, b3, a0, a1, a2, a3, t, tlog, fac, rlam
      integer :: nl,ic , iu_bio, lambda, lam(nlt)
      character title*50
      data a0,a1,a2,a3 /0.9976d0, 0.2194d0,  5.554d-2,  6.7d-3 /
      data b0,b1,b2,b3 /5.026d0, -0.01138d0, 9.552d-6, -2.698d-9/
#endif

      character(len=300) :: out_line
      character*6 :: skip
      character*12 :: fn='IiiiJjjj.rad'

C**** sync radiation parameters from input
      call sync_param( "S0X", S0X )
      call sync_param( "CO2X", CO2X )
      call sync_param( "N2OX", N2OX )
      call sync_param( "CH4X", CH4X )
      call sync_param( "CFC11X", CFC11X )
      call sync_param( "CFC12X", CFC12X )
      call sync_param( "XGHGX", XGHGX )
      call sync_param( "H2OstratX", H2OstratX )
      call sync_param( "O3X", O3X )
      call sync_param( "CLDX", CLDX )
      call sync_param( "H2ObyCH4", H2ObyCH4 )
      call sync_param( "S0_yr", S0_yr )
      call sync_param( "ghg_yr", ghg_yr )
      call sync_param( "ghg_day", ghg_day )
      call sync_param( "S0_day", S0_day )
      call sync_param( "volc_yr", volc_yr )
      call sync_param( "volc_day", volc_day )
      call sync_param( "aero_yr", aero_yr )
      call sync_param( "madaer", madaer )
      call sync_param( "dALBsnX", dALBsnX )
      call sync_param( "albsn_yr", albsn_yr )
      call sync_param( "aermix", aermix , 13 )
      call sync_param( "REFdry", REFdry , 8 )
      call sync_param( "FS8OPX", FS8OPX , 8 )
      call sync_param( "FT8OPX", FT8OPX , 8 )
      call sync_param( "RHfix", RHfix )
      call sync_param( "CC_cdncx", CC_cdncx )
      call sync_param( "OD_cdncx", OD_cdncx )
      call sync_param( "O3_yr", O3_yr )
      call sync_param( "PTLISO", PTLISO )
      call sync_param( "O3YR_max", O3YR_max )
      call sync_param( "KSOLAR", KSOLAR )
      call sync_param( "KSIALB", KSIALB )
      call sync_param( "KZSNOW", KZSNOW )
      call sync_param( "snoage_fac_max", snoage_fac_max )
      call sync_param( "nradfrc", nradfrc )
      if(snoage_fac_max.lt.0. .or. snoage_fac_max.gt.1.) then
        write(out_line,*) 'set 0<snoage_fac_max<1, not',snoage_fac_max
        call write_parallel(trim(out_line),unit=6)
        call stop_model('init_RAD: snoage_fac_max out of range',255)
      end if
      call sync_param( "rad_interact_aer", rad_interact_aer )
      call sync_param( "clim_interact_chem", clim_interact_chem )
      call sync_param( "rad_forc_lev", rad_forc_lev )
      call sync_param( "cloud_rad_forc", cloud_rad_forc )
      call sync_param( "aer_rad_forc", aer_rad_forc )
      call sync_param( "ref_mult", ref_mult )
      call sync_param( "save3dAOD", save3dAOD)
      REFdry = REFdry*ref_mult
      if (istart < 9) then
        call sync_param( "Ikliq", Ikliq )
        if(Ikliq.eq.1) kliq=1  ! hysteresis init: equilibrium
        if(Ikliq.eq.0) kliq=0  ! hysteresis init: dry
      end if
      if (istart.le.0) return

C**** Set orbital parameters appropriately
      select case (variable_orb_par)
      case(1) ! use parameters for model_year-orb_par_year_bp
        pyear = JYEAR-orb_par_year_bp ! bp=before present model year
        call orbpar(pyear,eccn, obliq, omegt)
        if (am_i_root()) then
          write(6,*) 'Variable orbital parameters, updated each year'
          write(6,*) 'Current orbital parameters from year',pyear
          write(6,*) '  Eccentricity:',eccn
          write(6,*) '  Obliquity (degs):',obliq
          write(6,*) '  Precession (degs from ve):',omegt
        end if
      case(0)  ! orbital parameters fixed from year orb_par_year_bp
        pyear=1950.-orb_par_year_bp ! here "present" means "1950"
        call orbpar(pyear, eccn, obliq, omegt)
        if (am_i_root()) then
          write(6,*) 'Fixed orbital parameters from year',pyear,' CE'
          write(6,*) '  Eccentricity:',eccn
          write(6,*) '  Obliquity (degs):',obliq
          write(6,*) '  Precession (degs from ve):',omegt
        end if
      case(-1) ! orbital parameters fixed, directly set
        eccn= orb_par(1) ; obliq=orb_par(2) ; omegt=orb_par(3)
        if (am_i_root()) then
          write(6,*) 'Orbital Parameters Specified:'
          write(6,*) '  Eccentricity:',eccn
          write(6,*) '  Obliquity (degs):',obliq
          write(6,*) '  Precession (degs from ve):',omegt
        end if
      case default  ! set from defaults (defined in CONSTANT module)
        omegt=omegt_def
        obliq=obliq_def
        eccn=eccn_def
      end select

      call cosz_init

C****
C**** SET THE CONTROL PARAMETERS FOR THE RADIATION (need mean pressures)
C****
      LMR=LM+LM_REQ
      PLB(1:LMR+1)=PEDNL00(1:LMR+1)
      DO L=1,LM
        PLBx(L)=PLB(L)           ! needed for CH4 prod. H2O
      END DO
      PLBx(LM+1)=0.
      DO LR=LM+1,LMR
        PLB0(LR-LM) = PLB(LR+1)
      END DO
      cdncl=0 ; call reterp(vcdnc,pcdnc,7, cdncl,plb,llow+2)

      KTREND=1   !  GHgas trends are determined by input file
!note KTREND=0 is a possible but virtually obsolete option
C****
C             Model Add-on Data of Extended Climatology Enable Parameter
C     MADO3M  = -1   Reads                      Ozone data the GCM way
C     MADAER  =  1   Reads   Tropospheric Aerosol climatology 1850-2050
C     MADAER  =  3   uses Koch,Bauer 2008 aerosol climatology 1890-2000
C     MADDST  =  1   Reads   Dust-windblown mineral climatology   RFILE6
C     MADVOL  =  1   Reads   Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS  =  1   Reads   Epsilon cloud heterogeniety data     RFILE8
C     MADLUV  =  1   Reads   Lean's SolarUV 1882-1998 variability RFILE9
C**** Radiative forcings are either constant = obs.value at given yr/day
C****    or time dependent (year=0); if day=0 an annual cycle is used
C****                                         even if the year is fixed
      KYEARS=s0_yr   ; KJDAYS=s0_day ;  MADLUV=1   ! solar 'constant'
      KYEARG=ghg_yr  ; KJDAYG=ghg_day              ! well-mixed GHGases
#ifndef ALTER_RADF_BY_LAT
      if(ghg_yr.gt.0)  MADGHG=0                    ! skip GHG-updating
#endif
      KYEARO=O3_yr   ; KJDAYO=0 ;       MADO3M=-1  ! ozone (ann.cycle)
      if(KYEARO.gt.0) KYEARO=-KYEARO              ! use ONLY KYEARO-data
      KYEARA=Aero_yr ; KJDAYA=0 ! MADAER=1 or 3, trop.aeros (ann.cycle)
      if(KYEARA.gt.0) KYEARA=-KYEARA              ! use ONLY KYEARA-data
      KYEARV=Volc_yr ; KJDAYV=Volc_day; MADVOL=1   ! Volc. Aerosols
!***  KYEARV=0 : use current year
!***  KYEARV<0 : use long term mean stratospheric aerosols (use -1)
!     Hack: KYEARV= -2000 and -2010 were used for 2 specific runs that
!           ended in 2100 and repeated some 20th century volcanos
!***  KYEARV=-2000: use volcanos from 100 yrs ago after 2000
!***  KYEARV=-2010: repeat 2nd half, then first half of 20th century
      if(KYEARV.le.-2000) KYEARV=0   ! use current year (before 2000)
C**** NO time history (yet), except for ann.cycle, for forcings below;
C****  if KJDAY?=day0 (1->365), data from that day are used all year
      KYEARD=0       ; KJDAYD=0 ;       MADDST=1   ! Desert dust
      KYEARE=0       ; KJDAYE=0 ;       MADEPS=1   !cloud Epsln - KCLDEP
      KYEARR=0       ; KJDAYR=0           ! surf.reflectance (ann.cycle)
      KCLDEM=1                  ! 0:old 1:new LW cloud scattering scheme

C**** Aerosols:
C**** Currently there are five different default aerosol controls
C****   1:total 2:background+tracer 3:Climatology 4:dust 5:volcanic
C**** By adjusting FSXAER,FTXAER you can remove the default
C**** aerosols and replace them with your version if required
C**** (through TRACER in RADIA).
C**** FSXAER is for the shortwave,    FTXAER is for the longwave effects
caer  FSXAER = (/ 1.,1.,1.,1.,1. /) ; FTXAER = (/ 1.,1.,1.,1.,1. /)

C**** climatology aerosols are grouped into 6 types from 13 sources:
C****  Pre-Industrial+Natural 1850 Level  Industrial Process  BioMBurn
C****  ---------------------------------  ------------------  --------
C****   1    2    3    4    5    6    7    8    9   10   11   12   13
C****  SNP  SBP  SSP  ANP  ONP  OBP  BBP  SUI  ANI  OCI  BCI  OCB  BCB
C**** using the following default scaling/tuning factors  AERMIX(1-13)
C****  1.0, 1.0, .26, 1.0, 2.5, 2.5, 1.9, 1.0, 1.0, 2.5, 1.9, 2.5, 1.9
C**** The 8 groups are (adding dust and volcanic aerosols as 7. and 8.)
C**** 1. Sulfates (industr and natural), 2. Sea Salt, 3. Nitrates
C**** 4. Organic Carbons, 5. industr Black Carbons(BC), 6. Biomass BC
C**** 7. Dust aerosols, 8. Volcanic aerosols
C**** use FS8OPX and FT8OPX to enhance the optical effect; defaults:
caer  FS8OPX = (/1., 1., 1., 1., 2., 2.,    1.   ,   1./)     solar
caer  FT8OPX = (/1., 1., 1., 1., 1., 1.,    1.3d0,   1./)     thermal
!!!!! Note: FS|T8OPX(7-8) makes FS|TXAER(4-5) redundant.
C**** Particle sizes of the first 4 groups have RelHum dependence

C**** To add up to 8 further aerosols:
C****  1) set NTRACE to the number of extra aerosol fields
C****  2) ITR defines which set of Mie parameters get used, choose
C****     from the following:
C****     1 SO4,  2 seasalt, 3 nitrate, 4 OCX organic carbons
C****     5 BCI,  6 BCB,     7 dust,    8 H2SO4 volc
C****  2b) set up the indexing array NTRIX to map the RADIATION tracers
C****      to the main model tracers
C****  2c) set up the weighting array WTTR to weight main model tracers
C****
C****  3) Use FSTOPX/FTTOPX(1:NTRACE) to scale them in RADIA
C****  4) Set TRRDRY to dry radius
C****  5) Set KRHTRA=1 if aerosol has RH dependence, 0 if not
C**** Note: whereas FSXAER/FTXAER are global (shared), FSTOPX/FTTOPX
C****       have to be reset for each grid box to allow for the way it
C****       is used in RADIA (TRACERS_AEROSOLS_Koch)
caer   NTRACE = 0
caer   ITR = (/ 0,0,0,0, 0,0,0,0 /)
caer   TRRDRY=(/ .1d0, .1d0, .1d0, .1d0, .1d0, .1d0, .1d0, .1d0/)
caer   KRHTRA=(/1,1,1,1,1,1,1,1/)

#if (defined TRACERS_AMP) || (defined TRACERS_AMP_M1)
      if (rad_interact_aer > 0) then
C                  SO4    SEA    NO3    OCX    BCI    BCB    DST   VOL
        FS8OPX = (/0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0 , 1d0/)
        FT8OPX = (/0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0,  1d0/)
      end if
      NTRACE=nmodes
      NTRIX(1:NMODES)=
     *     (/ n_N_AKK_1 ,n_N_ACC_1 ,n_N_DD1_1 ,n_N_DS1_1 ,n_N_DD2_1,
     *        n_N_DS2_1, n_N_SSA_1, n_N_SSC_1, n_N_OCC_1, n_N_BC1_1,
     *        n_N_BC2_1 ,n_N_BC3_1,
     *        n_N_DBC_1, n_N_BOC_1, n_N_BCS_1, n_N_MXX_1/)
#endif
#ifdef TRACERS_OM_SP
      if (rad_interact_aer > 0) then  ! if BC's sol.effect are doubled:
        FS8OPX = (/1d0, 1d0, 1d0, 0d0, 2d0, 2d0,  1d0 , 1d0/)
        FT8OPX = (/1d0, 1d0, 1d0, 0d0, 1d0, 1d0, 1.3d0, 1d0/)
      end if
      NTRACE=1
      TRRDRY(1:NTRACE)=(/ .3d0/)
      NTRIX(1:NTRACE)=
     *     (/n_OCA4/)
      WTTR(1:NTRACE) = 1d0
#endif
#ifdef TRACERS_AEROSOLS_Koch
      if (rad_interact_aer > 0) then  ! if BC's sol.effect are doubled:
#ifdef SULF_ONLY_AEROSOLS
        NTRACE=1
        FS8OPX(1:NTRACE) = (/0d0/)
        FT8OPX(1:NTRACE) = (/0d0/)
#else
c       FS8OPX = (/0d0, 0d0, 1d0, 0d0, 2d0, 2d0,  1d0 , 1d0/)
        FS8OPX = (/0d0, 0d0, 1d0, 0d0, 0d0, 0d0,  1d0 , 1d0/)
        FT8OPX = (/0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1.3d0, 1d0/)
#endif
      end if
#ifndef TRACERS_NITRATE
#ifdef SULF_ONLY_AEROSOLS
      NTRACE=1
      TRRDRY(1:NTRACE)=(/.15d0/)
      ITR(1:NTRACE) = (/1/)
      KRHTRA(1:NTRACE)=(/1/)
      NTRIX(1:NTRACE)=(/n_SO4/)
#else /* SULF aerosol allowed */
#ifdef TRACERS_AEROSOLS_SOA
      NTRACE=8
      TRRDRY(1:NTRACE)=
     &(/.15d0, .44d0, 1.7d0, .2d0, .2d0, .2d0, .08d0, .08d0/)
c augment BC by 50%
      FSTASC(1:NTRACE)=
     &(/1.d0,   1.d0, 1.d0 , 1.d0, 1.d0, 1.d0, 1.5d0, 1.5d0/)
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2, 4,4,4, 5,6/)
      KRHTRA(1:NTRACE)=(/1,1,1, 1,1,1, 0,0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=(
     &/n_sO4,n_seasalt1,n_seasalt2,n_OCIA,n_OCB,n_isopp1a,n_BCIA,n_BCB/)
#else /* OFF: TRACERS_AEROSOLS_SOA */
      NTRACE=7
      TRRDRY(1:NTRACE)=(/.15d0, .44d0, 1.7d0, .2d0, .2d0, .08d0, .08d0/)
c augment BC by 50%
      FSTASC(1:NTRACE)=(/1.d0,   1.d0, 1.d0 , 1.d0, 1.d0, 1.5d0, 1.5d0/)
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2,4,4, 5,6/)
      KRHTRA(1:NTRACE)=(/1,1,1,1,1, 0,0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=
     *(/ n_sO4, n_seasalt1, n_seasalt2, n_OCIA, n_OCB, n_BCIA, n_BCB/)
#endif /* TRACERS_AEROSOLS_SOA */
#endif /* NOT defied: SULF_ONLY_AEROSOLS */
C**** define weighting (only used for clays so far)
      WTTR(1:NTRACE) = 1d0

#else /* TRACERS_NITRATE ON */

#ifdef SULF_ONLY_AEROSOLS
      call stop_model('SULF_ONLY_AEROSOLS and TRACERS_NITRATE on',255)
#else /* SULF aerosol allowed */
#ifdef TRACERS_AEROSOLS_SOA
      NTRACE=9
      if (rad_interact_aer > 0) then ! turn off default nitrate
        FS8OPX(3) = 0. ; FT8OPX(3) = 0.
      end if
      TRRDRY(1:NTRACE)=
     * (/.15d0,.44d0, 1.7d0, .2d0, .2d0, .2d0, .08d0, .08d0,0.15d0/)
c augment BC by 50% (solar)
#ifdef TRACERS_NH3_UPDRAFT
c No scaling of Nitrate
      FSTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.d0, 1.5d0, 1.5d0,1.0d0/)
      FTTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.d0, 1.0d0, 1.0d0,1.0d0/)
#else
c scale down NO3p (solar and thermal) until tracer burden is reasonable
      FSTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.d0, 1.5d0, 1.5d0,0.2d0/)
      FTTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.d0, 1.0d0, 1.0d0,0.2d0/)
#endif
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2, 4,4,4, 5,6,3/)
      KRHTRA(1:NTRACE)=(/1,1,1, 1,1,1, 0,0,1/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=
     * (/ n_sO4, n_seasalt1, n_seasalt2, n_OCIA, n_OCB, n_isopp1a,
     *   n_BCIA, n_BCB, n_NO3p/)
#else /* OFF: TRACERS_AEROSOLS_SOA */
      NTRACE=8
      if (rad_interact_aer > 0) then ! turn off default nitrate
        FS8OPX(3) = 0. ; FT8OPX(3) = 0.
      end if
      TRRDRY(1:NTRACE)=
     * (/.15d0,.44d0, 1.7d0, .2d0, .2d0, .08d0, .08d0,0.15d0/)
c augment BC by 50% (solar)
#ifdef TRACERS_NH3_UPDRAFT
c No scaling of Nitrate
      FSTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.5d0, 1.5d0,1.0d0/)
      FTTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.0d0, 1.0d0,1.0d0/)
#else
c scale down NO3p (solar and thermal) until tracer burden is reasonable
      FSTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.5d0, 1.5d0,0.2d0/)
      FTTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.0d0, 1.0d0,0.2d0/)
#endif
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2,4,4, 5,6,3/)
      KRHTRA(1:NTRACE)=(/1,1,1,1,1, 0,0,1/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=
     * (/ n_sO4, n_seasalt1, n_seasalt2, n_OCIA, n_OCB, n_BCIA, n_BCB,
     *   n_NO3p/)
#endif /* TRACERS_AEROSOLS_SOA */
C**** define weighting (only used for clays so far)
      WTTR(1:NTRACE) = 1d0
#endif /* OFF: SULF_ONLY_AEROSOLS */
#endif /* TRACERS_NITRATE ON */
#endif /* TRACERS_AEROSOLS_Koch ON */

#ifdef TRACERS_DUST
C**** add dust optionally to radiatively active aerosol tracers
C**** should also work if other aerosols are not used
      if (rad_interact_aer > 0) then ! turn off default dust
        FS8OPX(7) = 0. ; FT8OPX(7) = 0.
      end if
      n1=NTRACE+1
      nrad_clay=n1
      NTRACE=NTRACE+ntm_dust+3  ! add dust tracers
c tracer 7 is dust
      ITR(n1:NTRACE) = 7
      KRHTRA(n1:NTRACE)= 0.  ! no deliq for dust

      SELECT CASE (ntm_dust)
      CASE (4)
C**** effective radii for dust
        TRRDRY(n1:NTRACE)=(/0.132D0,0.23D0,0.416D0,0.766D0,1.386D0,
     &       2.773D0,5.545D0/)
C**** Particle density of dust
        TRADEN(n1:NTRACE)=(/2.5D0,2.5D0,2.5D0,2.5D0,2.65D0,2.65D0,
     &       2.65D0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics. Adjust if number of dust tracers changes.
        NTRIX(n1:NTRACE)=(/n_clay,n_clay,n_clay,n_clay,n_silt1,n_silt2,
     &       n_silt3/)
C**** define weighting for different clays
        WTTR(n1:NTRACE)=(/0.009D0,0.081D0,0.234D0,0.676D0,1D0,1D0,1D0/)
      CASE (5)
        TRRDRY(n1:NTRACE)=(/0.132D0,0.23D0,0.416D0,0.766D0,1.386D0,
     &       2.773D0,5.545D0,11.090D0/)
        TRADEN(n1:NTRACE)=(/2.5D0,2.5D0,2.5D0,2.5D0,2.65D0,2.65D0,
     &       2.65D0,2.65D0/)
        NTRIX(n1:NTRACE)=(/n_clay,n_clay,n_clay,n_clay,n_silt1,n_silt2,
     &       n_silt3,n_silt4/)
        WTTR(n1:NTRACE)=(/0.009D0,0.081D0,0.234D0,0.676D0,1D0,1D0,1D0,
     &       1D0/)
      END SELECT
#else
#ifdef TRACERS_MINERALS
C**** add minerals optionally to radiatively active aerosol tracers
C**** so far all minerals have the properties of far traveled Saharan
C**** dust - to be changed soon
      if (rad_interact_aer > 0) then ! turn off default dust
        FS8OPX(7) = 0. ; FT8OPX(7) = 0.
      end if
      n1=NTRACE+1
      NTRACE = NTRACE + ntm_minerals ! add mineral tracers
c tracer 7 is dust
      ITR(n1:NTRACE) = (/7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
     &     7,7,7,7,7,7,7,7,7,7,7,7,7,7/)
      KRHTRA(n1:NTRACE)= 0.  ! no deliq for minerals
C**** effective radii for minerals
      TRRDRY(n1:NTRACE)=(/0.132D0,0.23D0,0.416D0,0.766D0,0.132D0,0.23D0,
     &     0.416D0,0.766D0,0.132D0,0.23D0,0.416D0,0.766D0,0.132D0,
     &     0.23D0,0.416D0,0.766D0,0.132D0,0.23D0,0.416D0,0.766D0,
     &     1.386D0,1.386D0,1.386D0,1.386D0,1.386D0,2.773D0,2.773D0,
     &     2.773D0,2.773D0,2.773D0,5.545D0,5.545D0,5.545D0,5.545D0,
     &     5.545D0/)
C**** Particle density of dust
      TRADEN(n1:NTRACE)=(/2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,
     &     2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,
     &     2.5D0,2.5D0,2.5D0,2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,
     &     2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,
     &     2.65D0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics. Adjust if number of dust tracers changes.
      NTRIX(n1:NTRACE)=(/n_clayilli,n_clayilli,n_clayilli,n_clayilli,
     &     n_claykaol,n_claykaol,n_claykaol,n_claykaol,n_claysmec,
     &     n_claysmec,n_claysmec,n_claysmec,n_claycalc,n_claycalc,
     &     n_claycalc,n_claycalc,n_clayquar,n_clayquar,n_clayquar,
     &     n_clayquar,n_sil1quar,n_sil1feld,n_sil1calc,n_sil1hema,
     &     n_sil1gyps,n_sil2quar,n_sil2feld,n_sil2calc,n_sil2hema,
     &     n_sil2gyps,n_sil3quar,n_sil3feld,n_sil3calc,n_sil3hema,
     &     n_sil3gyps/)
C**** define weighting for different clays
      WTTR(n1:NTRACE)=(/0.009D0,0.081D0,0.234D0,0.676D0,0.009D0,0.081D0,
     &     0.234D0,0.676D0,0.009D0,0.081D0,0.234D0,0.676D0,0.009D0,
     &     0.081D0,0.234D0,0.676D0,0.009D0,0.081D0,0.234D0,0.676D0,1D0,
     &     1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0/)
#endif
#ifdef TRACERS_QUARZHEM
C**** add quartz/hematite to radiatively active aerosol tracers
C**** so far all minerals have the properties of far traveled Saharan
C**** dust - to be changed soon
      if (rad_interact_aer > 0) then ! turn off default dust
        FS8OPX(7) = 0. ; FT8OPX(7) = 0.
      end if
      n1=NTRACE+1
      NTRACE = NTRACE + ntm_quarzhem ! add quartz/hematite aggregate tracers
c tracer 7 is dust
      ITR(n1:NTRACE) = (/7,7,7/)
      KRHTRA(n1:NTRACE)= 0.  ! no deliq for quartz/hematite aggregates
C**** effective radii for quartz/hematite
      TRRDRY(n1:NTRACE)=(/1.386D0,2.773D0,5.545D0/)
C**** Particle density of quartz/hematite
      TRADEN(n1:NTRACE)=(/2.65D0,2.65D0,2.65D0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics. Adjust if number of dust tracers changes.
      NTRIX(n1:NTRACE)=(/n_sil1quhe,n_sil2quhe,n_sil3quhe/)
C**** define weighting
      WTTR(n1:NTRACE)=(/1.D0,1.D0,1.D0/)
#endif
#endif

#ifdef TRACERS_SPECIAL_Shindell
      if(NTRACE > maxNtraceFastj)
     &call stop_model("NTRACE > maxNtraceFastj in init_Rad",13)
#endif

      if (ktrend.ne.0) then
C****   Read in time history of well-mixed greenhouse gases
        call openunit('GHG',iu,.false.,.true.)
        call ghghst(iu)
        call closeunit(iu)
        if (H2ObyCH4.ne.0..and.Kradia.le.0) then
C****     Read in dH2O: H2O prod.rate in kg/m^2 per day and ppm_CH4
          call openunit('dH2O',iu,.false.,.true.)
#if defined(CUBED_SPHERE)
          call read_qma(iu,plbx)
#else
          call getqma(iu,lat_dg,plbx,dh2o,lm,jm)
#endif
          call closeunit(iu)
        end if
      end if
      call updBCd(1990) ; depoBC_1990 = depoBC
C**** set up unit numbers for 14 more radiation input files
      DO IU=1,14
        IF (IU==12.OR.IU==13) CYCLE                    ! not used in GCM
        IF (IU==10.OR.IU==11) CYCLE                   ! obsolete O3 data
#ifdef USE_RADIATION_E1
        IF (IU==4 .OR.IU==5 ) CYCLE
#endif
        call openunit(RUNSTR(IU),NRFUN(IU),QBIN(IU),.true.)
      END DO

#ifdef CHL_from_SeaWIFs
C**** open chlorophyll data
      call openunit("CHL_DATA",iu_CHL,.true.,.true.)
#endif

      LS1_loc=1  ! default
C***********************************************************************
C     Main Radiative Initializations
C     ------------------------------------------------------------------
      CALL RCOMP1 (NRFUN)
      if (am_i_root()) CALL WRITER(6,0)  ! print rad. control parameters
C***********************************************************************
      DO IU=1,14
        IF (IU==12.OR.IU==13) CYCLE                    ! not used in GCM
        IF (IU==10.OR.IU==11) CYCLE                   ! obsolete O3 data
#ifdef USE_RADIATION_E1
        IF (IU==4 .OR.IU==5 ) CYCLE
#endif
        call closeunit(NRFUN(IU))
      END DO
C**** Save initial (currently permanent and global) Q in rad.layers
      do LR=1,LM_REQ
        shl0(LR) = shl(LM+LR)
      end do
      write(out_line,*) 'spec.hum in rad.equ.layers:',shl0
      call write_parallel(trim(out_line),unit=6)

#ifdef ALTER_RADF_BY_LAT
C**** Save initial rad forcing alterations:
      FS8OPX_orig(:)=FS8OPX(:); FT8OPX_orig(:)=FT8OPX(:) ! aerosols

C**** Read in the factors used for alterations:
      call openunit('ALT_GHG_LAT',iu2,.false.,.true.)
      read(iu2,*) ! skip first line
      do n=1,46
        read(iu2,'(a6,13D8.3)') skip,(FULGAS_lat(nn,n),nn=1,13)
      enddo
      call closeunit(iu2)
      call openunit('ALT_AER_LAT',iu2,.false.,.true.)
      read(iu2,*) ! skip first line
      do n=1,46
        read(iu2,'(a6,8D8.3)') skip,(FS8OPX_lat(nn,n),nn=1,8)
      enddo
      read(iu2,*) ! skip first line
      do n=1,46
        read(iu2,'(a6,8D8.3)') skip,(FT8OPX_lat(nn,n),nn=1,8)
      enddo
      call closeunit(iu2)
#endif

#if (defined OBIO_RAD_coupling) || (defined CHL_from_SeaWIFs)
      call openunit('cfle1',iu_bio,.false.,.true.)
      do ic = 1,6
        read(iu_bio,'(a50)')title
      enddo
      do nl = 1,nlt
        read(iu_bio,20) lambda,saw,sbw
        lam(nl) = lambda
        aw(nl) = saw
        bw(nl) = sbw
        if (lam(nl) .lt. 900) then
          t = exp(-(aw(nl)+0.5*bw(nl)))
          tlog = dlog(1.0D-36+t)
          fac = a0 + a1*tlog + a2*tlog*tlog + a3*tlog*tlog*tlog
          wfac(nl) = max(0d0,min(fac,1d0))
        else
          rlam = float(lam(nl))
          fac = b0 + b1*rlam + b2*rlam*rlam + b3*rlam*rlam*rlam
          wfac(nl) = max(fac,0d0)
        endif
      enddo
      print*,'RAD_DRV, wfac initializ= ', wfac
      call closeunit(iu_bio)
 20   format(i5,f15.4,f10.4)
#endif

#ifdef TRACERS_ON
c**** set tracerRadiaActiveFlag for radiatively active tracer
      do n=1,ntrace
        if (ntrix(n) > 0) tracerRadiaActiveFlag(ntrix(n))=.true.
      end do
      nTracerRadiaActive=count(tracerRadiaActiveFlag)
#endif

#ifndef CUBED_SPHERE
      if(kradia>0) then
!**** kradia>0 works on lat-lon grid only !!! open input files
        do j=grid%J_strt,grid%J_stop
        do i=grid%I_strt,imaxJ(J)
          write(fn,'(a1,i3.3,a1,i3.3,a4)') 'I',i,'J',j,'.rad'
          iu=1000*i+j
          open(iu,file=fn,form='unformatted')
        end do
        end do
      end if
#endif

      if (kradia .ge. 0) return
C**** pick selected representatives for 36x24 grid boxes if kradia<0
      allocate
     *   (ijradfrc(grid%I_strt:grid%I_stop,grid%J_strt:grid%J_stop))
!***  defaults (non-selected points)
      ijradfrc = 0
      do jc=1,24  ! 1,JcM
       do ic=1,36  ! 1,IcM
         if (ic > 1 .and. (jc==1 .or. jc==24)) cycle
         if (jc==1) then
           lonlat(1) = -179.9 ; lonlat(2) = -90.
         else if (jc==24) then
           lonlat(1) = -179.9 ; lonlat(2) = 90.
         else
           lonlat(1) = -176. + (ic-1)*10 ! lon
           lonlat(2) =  -92. + (jc-1)*8  ! lat
         end if
         call lonlat_to_ij (lonlat,ij)
         i = ij(1) ; j = ij(2)
         if (j < grid%J_strt .or. j > grid%J_stop) cycle
         if (i < grid%I_strt .or. i > grid%I_stop) cycle
         ijradfrc(i,j) = ic*1000 + jc
         write(fn,'(a1,i3.3,a1,i3.3,a4)') 'I',ic,'J',jc,'.rad'
         open(ijradfrc(i,j),file=fn,form='unformatted')
       end do
      end do

      if (istart < 10) return
!**** position the radfrc output files - kradia < 0 only
      do j=grid%J_strt,grid%J_stop
      do i=grid%I_strt,grid%I_stop
        if (ijradfrc(i,j)==0) cycle
        call io_POS_many(ijradfrc(i,j),Itime-1,Nrad)
      end do
      end do

      RETURN
      END SUBROUTINE init_RAD

      subroutine io_POS_many (iunit,it,itdif)
!@sum   io_POS  positions a seq. output file for the next write operat'n
!@auth  Original Development Team
!@ver   1.0

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT !@var IUNIT  file unit number
      INTEGER, INTENT(IN) :: IT,ITdif !@var IT,ITdif current time,dt

      INTEGER :: IT1,IT2   !@var  time_tags at start,end of each record
      INTEGER :: N              !@var  N      loop variable

      read (iunit,end=10,err=50) it1
      if(it1 .le. it) go to 30
   10 continue ! starting new file
      rewind iunit
      return

   20 read (iunit,end=35,err=50) it1
   30 continue
      if (it1 .le. it) go to 20
      it1=it1-itdif
   35 backspace iunit
      if (it1+itdif .le. it) go to 40
      return
   40 write (6,*) "file ",IUNIT," too short, it1/it=",it1,it
      call stop_model('io_POS: file too short',255)
   50 write (6,*) "Read error on: ",IUNIT,", it1/it=",it1,it
      call stop_model('io_POS: read error',255)
      END subroutine io_POS_many

      subroutine SETATM             ! dummy routine in gcm
      end subroutine SETATM

      subroutine GETVEG(LONR,LATR)  ! dummy routine in gcm
      integer LONR,LATR
      end subroutine GETVEG

      SUBROUTINE daily_RAD(end_of_day)
!@sum  daily_RAD sets radiation parameters that change every day
!@auth G. Schmidt
!@calls RADPAR:RCOMPT
      USE CONSTANT, only : by12
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP_ATM, only : am_I_root,GRID,GET,REWIND_PARALLEL
     *     ,READT_PARALLEL
      USE MODEL_COM, only : jday,jyear,im,jm,focean,jmon,JDendOfM
     *     ,JDmidOfM, jdate
      USE GEOM, only : imaxj
      USE RADPAR, only : FULGAS,JYEARR=>JYEAR,JDAYR=>JDAY
     *     ,xref,KYEARV
#ifdef ALTER_RADF_BY_LAT
     *     ,FULGAS_orig
#endif
      USE RADPAR, only : rcompt,writet
      USE RAD_COM, only : co2x,n2ox,ch4x,cfc11x,cfc12x,xGHGx,h2ostratx
     *     ,o3x,o3_yr,ghg_yr,co2ppm,Volc_yr,albsn_yr
#ifdef CHL_from_SeaWIFs
     *     ,iu_CHL,achl,echl1,echl0,bchl,cchl
      USE FLUXES, only : chl
#endif
      use DIAG_COM, only : iwrite,jwrite,itwrite
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: end_of_day
!@var TEMP_LOCAL stores ACHL+ECHL1 to avoid the use of common block
      REAL*8 :: TEMP_LOCAL(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                     GRID%J_STRT_HALO:GRID%J_STOP_HALO,2)
!@var IMON0 current month for CHL climatology reading
      INTEGER, SAVE :: IMON0 = 0
      INTEGER :: LSTMON,I,J
      REAL*8 TIME

      INTEGER :: J_0,J_1, I_0,I_1
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Update time dependent radiative parameters each day
!     Get black carbon deposition data for the appropriate year
!     (does nothing except at a restart or the beginning of a new year)
      if (albsn_yr.eq.0) then
        call updBCd (JYEAR)
      else
        call updBCd (albsn_yr)
      end if

!     Hack: 2 specific volc. eruption scenarios for 2000-2100 period
      if(volc_yr.eq.-2010) then              ! repeat some old volcanos
         KYEARV=JYEAR
         if(JYEAR.GT.2010) KYEARV=JYEAR-100  ! go back 100 years
      end if
      if(volc_yr.eq.-2000) then
         KYEARV=JYEAR
         if(JYEAR.GT.2000) KYEARV=JYEAR-50   ! go back 50 years til 2050
         if(JYEAR.GT.2050) KYEARV=JYEAR-150  ! then go back 150 years
      end if

      JDAYR=JDAY
      JYEARR=JYEAR
      CALL RCOMPT
!     FULGAS(2:) is set only in the first call to RCOMPT unless ghg_yr=0
!     Optional scaling of the observed value only in case it was (re)set
      if(ghg_yr.eq.0 .or. .not. end_of_day) then
         FULGAS(2)=FULGAS(2)*CO2X
         FULGAS(6)=FULGAS(6)*N2OX
         FULGAS(7)=FULGAS(7)*CH4X
         FULGAS(8)=FULGAS(8)*CFC11X
         FULGAS(9)=FULGAS(9)*CFC12X
         FULGAS(11)=FULGAS(11)*XGHGX
      end if
      IF(.not. end_of_day .and. H2OstratX.GE.0.)
     *     FULGAS(1)=FULGAS(1)*H2OstratX
      IF(.not. end_of_day .or. O3_yr==0.) FULGAS(3)=FULGAS(3)*O3X

C**** write trend table for forcing 'itwrite' for years iwrite->jwrite
C**** itwrite: 1-2=GHG 3=So 4-5=O3 6-9=aerosols: Trop,DesDust,Volc,Total
      if (am_i_root() .and.
     *  jwrite.gt.1500) call writet (6,itwrite,iwrite,jwrite,1,0)

#ifdef ALTER_RADF_BY_LAT
C**** Save initial rad forcing alterations:
      FULGAS_orig(:)=FULGAS(:) ! GHGs
#endif

C**** Define CO2 (ppm) for rest of model
      co2ppm = FULGAS(2)*XREF(1)

#ifdef CHL_from_SeaWIFs
C**** Read in Seawifs files here
      IF (JMON.NE.IMON0) THEN
      IF (IMON0==0) THEN
C**** READ IN LAST MONTH'S END-OF-MONTH DATA
        LSTMON=JMON-1
        if (lstmon.eq.0) lstmon = 12
        CALL READT_PARALLEL
     *       (grid,iu_CHL,NAMEUNIT(iu_CHL),TEMP_LOCAL,LSTMON)
        ECHL0 = TEMP_LOCAL(:,:,2)
      ELSE
C**** COPY END-OF-OLD-MONTH DATA TO START-OF-NEW-MONTH DATA
        ECHL0=ECHL1
      END IF
C**** READ IN CURRENT MONTHS DATA: MEAN AND END-OF-MONTH
      IMON0=JMON
      if (jmon.eq.1) CALL REWIND_PARALLEL( iu_CHL )
      CALL READT_PARALLEL
     *     (grid,iu_CHL,NAMEUNIT(iu_CHL),TEMP_LOCAL,1)
      ACHL  = TEMP_LOCAL(:,:,1)
      ECHL1 = TEMP_LOCAL(:,:,2)

C**** FIND INTERPOLATION COEFFICIENTS (LINEAR/QUADRATIC FIT)
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          BCHL(I,J)=ECHL1(I,J)-ECHL0(I,J)
          CCHL(I,J)=3.*(ECHL1(I,J)+ECHL0(I,J))-6.*ACHL(I,J)
        END DO
      END DO
      END IF
C**** Calculate CHL for current day
      TIME=(JDATE-.5)/(JDendOFM(JMON)-JDendOFM(JMON-1))-.5 ! -.5<TIME<.5
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          IF (FOCEAN(I,J).gt.0) THEN
C**** CHL always uses quadratic fit
            CHL(I,J)=ACHL(I,J)+BCHL(I,J)*TIME
     *           +CCHL(I,J)*(TIME**2-BY12)
            if (CHL(I,J).lt. 0) CHL(I,J)=0. ! just in case
          END IF
        END DO
      END DO
C**** REPLICATE VALUES AT POLE
      IF(HAVE_NORTH_POLE) then
       if (FOCEAN(1,JM).gt.0) CHL(2:IM,JM)=CHL(1,JM)
      ENDIF
      IF(HAVE_SOUTH_POLE) then
       if (FOCEAN(1, 1).gt.0) CHL(2:IM, 1)=CHL(1, 1)
      ENDIF
#endif

      RETURN
      END SUBROUTINE daily_RAD

      SUBROUTINE RADIA
!@sum  RADIA adds the radiation heating to the temperatures
!@auth Original Development Team
!@ver  1.0
!@calls tropwmo,coszs,coszt, RADPAR:rcompx ! writer,writet
      USE CONSTANT, only : sday,lhe,lhs,twopi,tf,stbo,rhow,mair,grav
     *     ,bysha,pi,radian
      USE MODEL_COM
      USE GEOM, only : imaxj, axyp, areag, byaxyp
     &     ,lat2d,lon2d
c      USE ATMDYN, only : CALC_AMPK
      USE RADPAR
     &  , only :  ! routines
     &           lx  ! for threadprivate copyin common block
     &          ,tauwc0,tauic0 ! set in radpar block data
     &          ,writer,rcompx,updghg, table
C     INPUT DATA         ! not (i,j) dependent
     X          ,S00WM2,RATLS0,S0,JYEARR=>JYEAR,JDAYR=>JDAY,FULGAS
     &          ,use_tracer_chem,FS8OPX,FT8OPX,use_o3_ref,KYEARG,KJDAYG
#ifdef ALTER_RADF_BY_LAT
     &          ,FS8OPX_orig,FT8OPX_orig,FULGAS_orig
#endif
C     INPUT DATA  (i,j) dependent
     &             ,JLAT,ILON, L1,LMR=>NL, PLB ,TLB,TLM ,SHL,RHL
     &             ,ltopcl,TAUWC ,TAUIC ,SIZEWC ,SIZEIC, kdeliq
     &             ,POCEAN,PEARTH,POICE,PLICE,PLAKE,COSZ,PVT
     &             ,TGO,TGE,TGOI,TGLI,TSL,WMAG,WEARTH
     &             ,AGESN,SNOWE,SNOWOI,SNOWLI,dALBsn, ZSNWOI,ZOICE
     &             ,zmp,fmp,flags,LS1_loc,snow_frac,zlake
     *             ,TRACER,NTRACE,FSTOPX,FTTOPX,chem_IN,O3natL,O3natLref
     *             ,FTAUC,LOC_CHL,FSTASC,FTTASC
#ifdef HEALY_LM_DIAGS
     *             ,VTAULAT 
#endif

C     OUTPUT DATA
     &          ,TRDFLB ,TRNFLB ,TRUFLB, TRFCRL ,chem_out
     &          ,SRDFLB ,SRNFLB ,SRUFLB, SRFHRL
     &          ,PLAVIS ,PLANIR ,ALBVIS ,ALBNIR ,FSRNFG
     &          ,SRRVIS ,SRRNIR ,SRAVIS ,SRANIR ,SRXVIS ,SRDVIS
     &          ,BTEMPW ,TTAUSV ,SRAEXT ,SRASCT ,SRAGCB
     &          ,SRDEXT ,SRDSCT ,SRDGCB ,SRVEXT ,SRVSCT ,SRVGCB
     &          ,aesqex,aesqsc,aesqcb
     &          ,SRXNIR,SRDNIR
      USE RAD_COM, only : rqt,srhr,trhr,fsf,cosz1,s0x,rsdist,nradfrc
     *     ,plb0,shl0,tchg,alb,fsrdir,srvissurf,srdn,cfrac,rcld
     *     ,chem_tracer_save,rad_interact_aer,kliq,RHfix,CLDx
     *     ,ghg_yr,CO2X,N2OX,CH4X,CFC11X,CFC12X,XGHGX,rad_forc_lev,ntrix
     *     ,wttr,cloud_rad_forc,CC_cdncx,OD_cdncx,cdncl,nrad_clay
     *     ,dALBsnX,depoBC,depoBC_1990,rad_to_chem,trsurf
     *     ,FSRDIF,DIRNIR,DIFNIR,aer_rad_forc,clim_interact_chem
     *     ,TAUSUMW,TAUSUMI, ijradfrc ! needed only if kradia<0
#ifdef ALTER_RADF_BY_LAT
     *     ,FULGAS_lat,FS8OPX_lat,FT8OPX_lat
#endif
#ifdef TRACERS_DUST
     &     ,srnflb_save,trnflb_save
#endif
#ifdef TRACERS_ON
     &     ,ttausv_sum,ttausv_sum_cs,ttausv_count,ttausv_save
     &     ,ttausv_cs_save,aerAbs6SaveInst
#endif
#ifdef TRACERS_SPECIAL_Shindell
     &     ,ttausv_ntrace
#endif
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     &     ,stratO3_tracer_save
#endif
      USE RANDOM
      USE CLOUDS_COM, only : tauss,taumc,svlhx,rhsav,svlat,cldsav,
     *     cldmc,cldss,csizmc,csizss,llow,lmid,lhi,fss
      USE PBLCOM, only : wsavg,tsavg
#ifdef SCM
      USE SCMCOM, only : SCM_SURF_ALBEDO_FLAG,iu_scm_prt
      USE SCMDIAG, only : SRDFLBTOP,SRNFLBTOP,SRUFLBTOP,TRUFLBTOP,
     *                    SRDFLBBOT,SRNFLBBOT,SRUFLBBOT,TRUFLBBOT,
     *                    TRDFLBBOT,TRDFLBTOP,SRFHRLCOL,TRFCRLCOL,
     *                    CSSRNTOP,CSTRUTOP,CSSRNBOT,CSTRNBOT,
     *                    CSSRDBOT,TRNFLBBOT,dTradlw,dTradsw,
     *                    scm_FDIF_vis,scm_FDIR_vis,scm_FDIF_nir,
     *                    scm_FDIR_nir
      USE RADPAR, only : KEEPAL,SRBALB
#endif
      USE DIAG_COM, only : ia_rad,jreg,aij=>aij_loc,aijl=>aijl_loc
     *     ,adiurn=>adiurn_loc,ndiuvar,ia_rad_frc,
#ifndef NO_HDIURN
     *     hdiurn=>hdiurn_loc,
#endif
     *     iwrite,jwrite,itwrite,ndiupt,j_pcldss,j_pcldmc,ij_pmccld,
     *     j_clddep,j_pcld,ij_cldcv,ij_pcldl,ij_pcldm,ij_pcldh,
     *     ij_cldtppr,j_srincp0,j_srnfp0,j_srnfp1,j_srincg,
     *     j_srnfg,j_brtemp,j_trincg,j_hsurf,j_hatm,ij_trnfp0,
     *     j_plavis, j_planir, j_albvis, j_albnir,
     *     j_srrvis, j_srrnir, j_sravis, j_sranir,
     *     ij_srnfp0,ij_srincp0,ij_srnfg,ij_srincg,ij_btmpw,ij_srref
     *     ,ij_srvis,ij_rnfp1,j_clrtoa,j_clrtrp,j_tottrp,ijl_rc
     *     ,ijdd,idd_cl7,idd_ccv,idd_isw,idd_palb,idd_galb, idd_aot
     *     ,idd_absa,jl_srhr,jl_trcr,jl_totcld,jl_sscld,jl_mccld
     *     ,ij_frmp,jl_wcld,jl_icld,jl_wcod,jl_icod,jl_wcsiz,jl_icsiz
     *     ,jl_wcldwt,jl_icldwt
     *     ,ij_clr_srincg,ij_CLDTPT,ij_cldt1t,ij_cldt1p,ij_cldcv1
     *     ,ij_wtrcld,ij_icecld,ij_optdw,ij_optdi,ij_swcrf,ij_lwcrf
     *     ,AFLX_ST, hr_in_day,hr_in_month,ij_srntp,ij_trntp
     *     ,ij_clr_srntp,ij_clr_trntp,ij_clr_srnfg,ij_clr_trdng
     *     ,ij_clr_sruptoa,ij_clr_truptoa,ijl_cf
     *     ,ij_swdcls,ij_swncls,ij_lwdcls,ij_swnclt,ij_lwnclt, NREG
     *     ,adiurn_dust,j_trnfp0,j_trnfp1,ij_srvdir, ij_srvissurf
     *     ,ij_chl, ij_swaerrf, ij_lwaerrf,ij_swaersrf,ij_lwaersrf
     *     ,ij_swaerrfnt,ij_lwaerrfnt,ij_swaersrfnt,ij_lwaersrfnt
     *     ,ij_swcrf2,ij_lwcrf2, ij_siswd, ij_siswu, save3dAOD
#ifdef ACCMIP_LIKE_DIAGS
     *     ,ij_fcghg ! array
#endif
#ifdef HEALY_LM_DIAGS
     *     ,j_vtau,j_ghg
#endif

      USE DYNAMICS, only : pk,pedn,pmid,pdsig,ltropo,am,byam
      USE SEAICE, only : rhos,ace1i,rhoi
      USE SEAICE_COM, only : rsi,snowi,pond_melt,msi,flag_dsws
      USE GHY_COM, only : snowe_com=>snowe,snoage,wearth_com=>wearth
     *     ,aiearth,fr_snow_rad_ij,fearth
#ifdef USE_ENT
      use ent_com, only : entcells
      use ent_mod, only : ent_get_exports
     &                    ,N_COVERTYPES !YKIM-temp hack
      use ent_drv, only : map_ent2giss  !YKIM-temp hack
#else
      USE VEG_COM, only : vdata
#endif
      USE LANDICE_COM, only : snowli_com=>snowli
      USE LAKES_COM, only : flake,mwl
      USE FLUXES, only : nstype,gtempr,bare_soil_wetness
#if (defined CHL_from_SeaWIFs) || (defined TRACERS_OceanBiology)
     .                  ,chl
#endif
      USE DOMAIN_DECOMP_ATM, ONLY: grid,GET, write_parallel
      USE DOMAIN_DECOMP_ATM, ONLY: GLOBALSUM
      USE RAD_COSZ0, only : COSZT,COSZS

#ifdef TRACERS_ON
      USE TRACER_COM, only: NTM,n_Ox,trm,trname,n_OCB,n_BCII,n_BCIA
     *     ,n_OCIA,N_OCII,n_so4_d2,n_so4_d3,trpdens,n_SO4,n_stratOx
     *     ,n_OCI1,n_OCI2,n_OCI3,n_OCA1,n_OCA2,n_OCA3,n_OCA4,n_N_AKK_1
#ifdef TRACERS_NITRATE
     *     ,tr_mm,n_NH4,n_NO3p
#endif
#ifdef TRACERS_AEROSOLS_SOA
     *     ,n_isopp1a,n_isopp2a
#ifdef TRACERS_TERP
     *     ,n_apinp1a,n_apinp2a
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_Koch
c    *     ,SNFST0,TNFST0
#endif  /* TRACERS_AEROSOLS_Koch */
      USE TRDIAG_COM, only: taijs=>taijs_loc,taijls=>taijls_loc,ijts_fc
     *     ,ijts_tau,ijts_tausub,ijts_fcsub,ijlt_3dtau,ijts_sqex
     *     ,ijts_sqexsub,ijts_sqsc,ijts_sqscsub,ijts_sqcb,ijts_sqcbsub
     *     ,diag_rad
#ifdef AUXILIARY_OX_RADF
     *     ,ijts_auxfc
#endif /* AUXILIARY_OX_RADF */
#ifdef BC_ALB
     *     ,ijts_alb
#endif  /* BC_ALB */
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: Lmax_rad_O3,Lmax_rad_CH4
#endif /* TRACERS_SPECIAL_Shindell */
#endif /* TRACERS_ON */
#ifdef TRACERS_AMP
      USE AERO_CONFIG, only: nmodes
      USE AMP_AEROSOL, only: AMP_DIAG_FC
#endif
#ifdef RAD_O3_GCM_HRES
      use RAD_native_O3, only: O3JDAY_native,O3JREF_native
#endif
      use IndirectAerParam_mod, only: dCDNC_est
      USE TimerPackage_mod, only: startTimer => start, stopTimer => stop
      USE PARAM, only : get_param, is_set_param
#ifdef CACHED_SUBDD
      use subdd_mod, only : sched_rad, subdd_groups, subdd_type
     &     ,subdd_ngroups,inc_subdd,find_groups, lmaxsubdd
#endif
      IMPLICIT NONE
C
C     INPUT DATA   partly (i,j) dependent, partly global
      REAL*8 U0GAS,taulim, xdalbs,sumda,tauda,fsnow,fsoot
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     sumda_psum,tauda_psum
      COMMON/RADPAR_hybrid/U0GAS(LX,13)
!$OMP  THREADPRIVATE(/RADPAR_hybrid/)

      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     COSZ2,COSZA,TRINCG,BTMPW
      REAL*8, DIMENSION(4,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                    grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     SNFS,TNFS
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     SNFSCRF,TNFSCRF,SNFSCRF2,TNFSCRF2
      REAL*8, DIMENSION(18,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     SNFSAERRF,TNFSAERRF

#ifdef CACHED_SUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups)
      type(subdd_type), pointer :: subdd
!@var OLR_ij, OLR_cnt --  Net thermal radiation at TOA (W/m^2) for SUBDD
      REAL*8  :: OLR_cnt
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     OLR_ij,SWU_ij,SDDARR
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &     SDDARR3D
#endif
#ifdef ACCMIP_LIKE_DIAGS
!@var snfs_ghg,tnfs_ghg like SNFS/TNFS but with reference GHG for
!@+   radiative forcing calculations. TOA only.
!@+   index 1=CH4, 2=N2O, 3=CFC11, 4=CFC12.
      real*8,dimension(4,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                   grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &                   snfs_ghg,tnfs_ghg
      real*8,dimension(4) :: sv_fulgas_ref,sv_fulgas_now
      integer :: nf,GFrefY,GFrefD,GFnowY,GFnowD
!@var nfghg fulgas( ) index of radf diag ghgs:
      integer, dimension(4) :: nfghg=(/7,6,8,9/)
#endif
#ifdef HEALY_LM_DIAGS
C  GHG Effective forcing relative to 1850
       real*8 :: ghg_totforc, CO2I=285.2,N2OI=.2754,CH4I=.791  ! 1850  GHG's
       real*8 ::              CO2R=337.9,N2OR=.3012,CH4R=1.547 ! RAD's 1979 Reference values
       real*8 ::              FCO2,FN2O,FCH4                   ! Current Model GHG
       real*8 :: Fe !! Function
#endif
#ifdef TRACERS_ON
!@var StauL sum of ttausv over L
      real*8 :: StauL
!@var SNFST,TNFST like SNFS/TNFS but with/without specific tracers for
!@+   radiative forcing calculations
      REAL*8,DIMENSION(2,NTRACE,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                          grid%J_STRT_HALO:grid%J_STOP_HALO)::
     *     SNFST,TNFST
!@var SNFST_o3ref,TNFST_o3ref like snfst,tnfst for special case ozone for
!@+   which ntrace fields are not defined. Indicies are :
!@+   1=LTROPO,reference, 2=TOA,reference; not saving surface forcing.
!@+   3=LTROPO or LS1-1,auxiliary, 4=TOA,auxiliary; 5=LS1-1,reference
      REAL*8,DIMENSION(5,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                   grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     SNFST_o3ref,TNFST_o3ref
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     &    ,snfst_stratOx,tnfst_stratOx
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS */
#ifdef BC_ALB
      REAL*8,DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                 grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     ALBNBC,NFSNBC
#endif /* BC_ALB */
#endif /* TRACERS_ON */
      REAL*8, DIMENSION(LM_REQ,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                         grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     TRHRS,SRHRS
      REAL*8, DIMENSION(LM+LM_REQ) :: TRHRA,SRHRA ! kradia>1 only
      REAL*8, DIMENSION(LM+LM_REQ) :: aml,pdsigl,pmidl,pcoll ! kradia>0
      REAL*8, DIMENSION(LM) :: TOTCLD,dcc_cdncl,dod_cdncl
      INTEGER I,J,L,K,KR,LR,JR,IH,IHM,INCH,JK,IT,iy,iend,N,onoff_aer
     *     ,onoff_chem,LFRC,JTIME,n1,tmpS(8),tmpT(8),moddrf
      REAL*8 ROT1,ROT2,PLAND,CSS,CMC,DEPTH,QSS,TAUSSL,RANDSS
     *     ,TAUMCL,ELHX,CLDCV,X,OPNSKY,CSZ2,tauup,taudn,ptype4(4)
     *     ,taucl,wtlin,MSTRAT,STRATQ,STRJ,MSTJ,optdw,optdi,rsign_aer
     *     ,rsign_chem,tauex5,tauex6,tausct,taugcb,dcdnc
      REAL*8 RANDXX ! temporary
      REAL*8 QSAT
#ifdef BC_ALB
      REAL*8 dALBsn1
#endif
      LOGICAL NO_CLOUD_ABOVE, set_clayilli,set_claykaol,set_claysmec,
     &     set_claycalc,set_clayquar
C
      REAL*8  RDSS(LM,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,RDMC(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &           grid%J_STRT_HALO:grid%J_STOP_HALO)

      REAL*8 :: TMP(NDIUVAR)
      INTEGER, PARAMETER :: NLOC_DIU_VAR = 8
      INTEGER :: idx(NLOC_DIU_VAR)

      INTEGER, PARAMETER :: NLOC_DIU_VARb = 3
      INTEGER :: idxb(NLOC_DIU_VARb)

      integer :: aj_alb_inds(8)
      real*8, dimension(lm_req) :: bydpreq

c     INTEGER ICKERR,JCKERR,KCKERR
      INTEGER :: J_0, J_1, I_0, I_1
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      character(len=300) :: out_line
#ifdef SCM
      integer iscm
#endif

      integer :: nij_before_j0,nij_after_j1,nij_after_i1
      integer :: initial_GHG_setup

#ifdef USE_ENT
      real*8 :: PVT0(N_COVERTYPES)
#endif
#ifdef TRACERS_NITRATE
      real*8 :: nh4_on_no3
#endif

C
C****
      call startTimer('RADIA()')

      idx = (/ (IDD_CL7+i-1,i=1,7), IDD_CCV /)
      idxb = (/ IDD_PALB, IDD_GALB, IDD_ABSA /)

      Call GET(grid, HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP


C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** GTEMPR RADIATIVE TEMPERATURE ARRAY OVER ALL SURFACE TYPES (K)
C****   RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****
C**** VDATA  1-11 RATIOS FOR THE 11 VEGETATION TYPES (1)
C****

C**** limit optical cloud depth from below: taulim
      taulim=min(tauwc0,tauic0) ! currently both .001
      tauwc0 = taulim ; tauic0 = taulim
C**** Calculate mean cosine of zenith angle for the current physics step
      JTIME=MOD(ITIME,NDAY)
      ROT1=(TWOPI*JTIME)/NDAY
c      ROT2=ROT1+TWOPI*DTsrc/SDAY
c      CALL COSZT (ROT1,ROT2,COSZ1)
      CALL CALC_ZENITH_ANGLE   ! moved to main loop

      IF (MODRD.NE.0) GO TO 900
      IDACC(ia_rad)=IDACC(ia_rad)+1
      moddrf=1    ! skip rad.forcing diags if nradfrc.le.0
      if (nradfrc>0) moddrf=mod(itime-itimei,nrad*nradfrc)
C****
      if (moddrf==0) IDACC(ia_rad_frc)=IDACC(ia_rad_frc)+1
C**** Interface with radiation routines, done only every NRAD time steps
C****
C**** Calculate mean cosine of zenith angle for the full radiation step
      ROT2=ROT1+TWOPI*NRAD*DTsrc/SDAY
      CALL COSZS (ROT1,ROT2,COSZ2,COSZA)
      JDAYR=JDAY
      JYEARR=JYEAR

      S0=S0X*S00WM2*RATLS0/RSDIST

      if (kradia .le. 0) then ! full model
c**** find scaling factors for surface albedo reduction
#ifdef SCM
      xdalbs = 0.d0
#else
      IF (HAVE_SOUTH_POLE) THEN
         sumda_psum(:,1)=axyp(1,1)
         tauda_psum(:,1)=axyp(1,1)*depobc_1990(1,1)
      End If
      do j=J_0S,J_1S
      do i=I_0,I_1
c ilon, jlat are indices w.r.t 72x46 grid
c      JLAT=INT(1.+(J-1.)*0.25*DLAT_DG+.5)   ! slightly more general
c      ILON=INT(.5+(I-.5)*72./IM+.5)
        ilon = 1 + int( 72d0*lon2d(i,j)/twopi )
        jlat = 1 + int( 45d0*(lat2d(i,j)+92d0*radian)/pi )
        fsnow = flice(i,j) + rsi(i,j)*(1-fland(i,j))
        if(SNOWE_COM(I,J).gt.0.) fsnow = fsnow+fearth(i,j)
        sumda_psum(i,j) = axyp(i,j)*fsnow
        tauda_psum(i,j) = axyp(i,j)*fsnow*depobc_1990(ilon,jlat)
      end do
      end do
      IF (HAVE_NORTH_POLE) THEN
         sumda_psum(:,JM)=axyp(1,jm)*rsi(1,jm)
         tauda_psum(:,JM)=axyp(1,jm)*rsi(1,jm)*depobc_1990(1,46)
      END IF
      CALL GLOBALSUM(grid, sumda_psum,sumda,all=.true.)
      CALL GLOBALSUM(grid, tauda_psum,tauda,all=.true.)

      xdalbs=-dalbsnX*sumda/tauda ; fsoot = sumda/tauda
      IF(QCHECK) write(6,*) 'coeff. for snow alb reduction',xdalbs
#endif

      IF (QCHECK) THEN
C****   Calculate mean strat water conc
        STRATQ=0.
        MSTRAT=0.
        DO J=J_0,J_1
          STRJ=0.
          MSTJ=0.
          DO I=I_0,IMAXJ(J)
            DO L=LTROPO(I,J)+1,LM
              STRJ=STRJ+Q(I,J,L)*AM(L,I,J)*AXYP(I,J)
              MSTJ=MSTJ+AM(L,I,J)*AXYP(I,J)
            END DO
          END DO
          IF (J.eq.1 .or. J.eq.JM) THEN
            STRJ=STRJ*FIM
            MSTJ=MSTJ*FIM
          END IF
          STRATQ=STRATQ+STRJ
          MSTRAT=MSTRAT+MSTJ
        END DO
        PRINT*,"Strat water vapour (ppmv), mass (mb)",1d6*STRATQ*mair
     *       /(18.*MSTRAT),PMTOP+1d-2*GRAV*MSTRAT/AREAG
      END IF

C**** Get the random numbers outside openMP parallel regions
C**** but keep MC calculation separate from SS clouds
C**** To get parallel consistency also with mpi, force each process
C**** to generate random numbers for all latitudes (using BURN_RANDOM)

C**** MC clouds are considered as a block for each I,J grid point

      CALL BURN_RANDOM(nij_before_j0(J_0))

      DO J=J_0,J_1                    ! complete overlap
      CALL BURN_RANDOM((I_0-1))
      DO I=I_0,IMAXJ(J)
        RDMC(I,J) = RANDU(X)
      END DO
      CALL BURN_RANDOM(nij_after_i1(I_1))
      END DO

      CALL BURN_RANDOM((nij_after_j1(J_1)))

C**** SS clouds are considered as a block for each continuous cloud
      CALL BURN_RANDOM(nij_before_j0(j_0)*LM)

      DO J=J_0,J_1                    ! semi-random overlap
      CALL BURN_RANDOM((I_0-1)*LM)
      DO I=I_0,IMAXJ(J)
        NO_CLOUD_ABOVE = .TRUE.
        DO L=LM,1,-1
          IF(TAUSS(L,I,J).le.taulim) CLDSS(L,I,J)=0.
          IF(TAUMC(L,I,J).le.taulim) CLDMC(L,I,J)=0.
          RANDXX = RANDU(X)
          IF(CLDSS(L,I,J).GT.0.) THEN
            IF (NO_CLOUD_ABOVE) THEN
              RANDSS = RANDXX
              NO_CLOUD_ABOVE = .FALSE.
            END IF
          ELSE
            RANDSS = 1.
            NO_CLOUD_ABOVE = .TRUE.
          END IF
          RDSS(L,I,J) = RANDSS
        END DO
      END DO
      CALL BURN_RANDOM(nij_after_i1(I_1)*LM)
      END DO

      CALL BURN_RANDOM(nij_after_j1(j_1)*LM)

      end if                    ! kradia le 0

#ifdef TRACERS_ON
      ttausv_count=ttausv_count+1.d0
#endif
#ifdef ACCMIP_LIKE_DIAGS
! because of additional updghg calls, these factors won't apply:
      if(CO2X.ne.1.)  call stop_model('CO2x.ne.1 accmip diags',255)
      if(N2OX.ne.1.)  call stop_model('N2Ox.ne.1 accmip diags',255)
      if(CH4X.ne.1.)  call stop_model('CH4x.ne.1 accmip diags',255)
      if(CFC11X.ne.1.)call stop_model('CFC11x.ne.1 accmip diags',255)
      if(CFC12X.ne.1.)call stop_model('CFC12x.ne.1 accmip diags',255)
      if(XGHGX.ne.1.) call stop_model('XGHGx.ne.1 accmip diags',255)
      GFrefY=1850; GFrefD=182     ! ghg forcing refrnce year, day
      GFnowY=JyearR; GFnowD=JdayR ! ghg current desired year, day
      if(KJDAYG > 0) GFnowD=KJDAYG ! unless presribed in deck
      if(KYEARG > 0) GFnowY=KYEARG !           "
      call updghg(GFrefY,GFrefD)
      sv_fulgas_ref(1:4)=fulgas(nfghg(1:4))
      call updghg(GFnowY,GFnowD)
      sv_fulgas_now(1:4)=fulgas(nfghg(1:4))
#endif
#ifdef HEALY_LM_DIAGS
      FCO2=FULGAS(2)*CO2R
      FN2O=FULGAS(6)*N2OR
      FCH4=FULGAS(7)*CH4R
c      
c      write(6,*) 'RJH: GHG: CONC=',
c     * FCO2,FN2O,FCH4
      ghg_totforc=5.35d0*log(FCO2/CO2I)
     *  +0.036d0*(sqrt(FCH4)-sqrt(CH4I))
     *  -(Fe(FCH4,N2OI)-Fe(CH4I,N2OI))
     *  + 0.12d0*(sqrt(FN2O)-sqrt(N2OI))
     *  -(Fe(CH4I,FN2O)-Fe(CH4I,N2OI))
c      write(6,*) 'RJH: GHG: FORC=',ghg_totforc
#endif

      aj_alb_inds = (/ J_PLAVIS, J_PLANIR, J_ALBVIS, J_ALBNIR,
     &                 J_SRRVIS, J_SRRNIR, J_SRAVIS, J_SRANIR /)

C****
C**** MAIN J LOOP
C****
      DO J=J_0,J_1

c     ICKERR=0
c     JCKERR=0
c     KCKERR=0

!$OMP  PARALLEL DEFAULT(NONE)
!$OMP*   PRIVATE(CSS,CMC,CLDCV, DEPTH,OPTDW,OPTDI, ELHX,
!$OMP*   tmp,tmpS,tmpT,StauL,dALBsn1,
!$OMP*   I,INCH,IH,IHM,IT,JR, K,KR, L,LR,LFRC, N, onoff_aer,
!$OMP*   onoff_chem, OPNSKY, CSZ2, PLAND,ptype4,tauex5,tauex6,tausct,
!$OMP*   taugcb,set_clayilli,set_claykaol,set_claysmec,set_claycalc,
!$OMP*   set_clayquar,dcc_cdncl,dod_cdncl,dCDNC,n1,pvt0,
#ifdef ALTER_RADF_BY_LAT
!$OMP*   fulgas,FS8OPX,FT8OPX,
#endif
!$OMP*   QSS, TOTCLD,TAUSSL,TAUMCL,tauup,taudn,taucl,wtlin)
!$OMP*   COPYIN(/RADPAR_hybrid/)
!$OMP*   SHARED(ITWRITE, J,I_0,IMAXJ,lon2d,lat2d,jreg,fland,rsi,
!$OMP*    flake,flice,fearth,gtempr,itime,ltropo,t,pk,cldx,
!$OMP*    cc_cdncx,cdncl,od_cdncx,q,rhsav,fss,cldsav,cldss,rdss,
!$OMP*    tauss,kradia,cldmc,rdmc,pdsig,taumc,pmid,fsoot,
!$OMP*   AIJL,CSIZMC,svlat,csizss,svlhx,cfrac,ftype,
!$OMP*   AIJ,llow,lmid,lhi,
!$OMP*   nrad,jtime,nday,
!$OMP*   idx,pedn,rhfix,ntrace,ntrix,n_ocb, trm,n_ocii,n_ocia,
!$OMP*   idd_cl7,idd_ccv,adiurn,byaxyp,n_oci1,n_oca1,n_oci2,
!$OMP*   n_oca3,n_oca4,n_bcii,n_bcia,wttr,rqt,plb0,shl0,tchg,aflx_st,
!$OMP*   cosza, tsavg, snowi, snowli_com,snowe_com,fr_snow_rad_ij,
!$OMP*   snoage,xdalbs,depobc,msi,n_oca2,n_oci3,flag_dsws,
!$OMP*   pond_melt,mwl,axyp,bare_soil_wetness,
!$OMP*   entcells,wsavg,rad_forc_lev,rad_interact_aer,
!$OMP*   clim_interact_chem,chem_tracer_save,lmax_rad_o3,lmax_rad_ch4,
!$OMP*   kliq,snfst,tnfst,snfst_ozone,tnfst_ozone,nfsnbc,albnbc,
!$OMP*   cloud_rad_forc,snfscrf,tnfscrf,
#ifdef ALTER_RADF_BY_LAT
!$OMP*   fulgas_orig,FS8OPX_orig,FT8OPX_orig,
#endif
!$OMP*   cosz2,aer_rad_forc,
!$OMP*   FS8OPX,FT8OPX,SNFSAERRF,TNFSAERRF,nrad_clay,diag_rad,
!$OMP*   ttausv_cs_save, iwrite, jwrite, rad_to_chem, srhra,trhra,
!$OMP*   DTsrc,byam,byaml00,fsf,srhr,trhr,trsurf,srhrs,trhrs,snfs,tnfs,
!$OMP*   trincg,btmpw,alb,aj_alb_inds,srnflb_save,trnflb_save,srdn,
!$OMP*   FSRDIR,SRVISSURF,FSRDIF,DIRNIR,DIFNIR,taijls,taijs,
!$OMP*   ttausv_sum,ttausv_sum_cs,adiurn_dust,ttausv_save,
#ifdef TRACERS_SPECIAL_Shindell
!$OMP*   ttausv_ntrace,
#endif
!$OMP*   rcld,IJ_PMCCLD,ij_cldcv,ij_pcldl,ij_pcldh,jl_sscld,jl_mccld,
!$OMP*   jl_totcld,IJL_CF,jl_wcld,jl_wcldwt,jl_icld,ij_swdcls,ij_frmp,
!$OMP*   jl_icldwt,jl_wcod,jl_icod,jl_wcsiz,jl_icsiz,ijdd,ij_pcldm,
!$OMP*   j_pcldss,j_pcldmc,j_clddep,j_pcld,ij_swnclt,ij_lwnclt,
!$OMP*   ijts_tausub,ijlt_3dtau,ijts_sqexsub,ij_swncls,ij_lwdcls,
!$OMP*   ijts_sqscsub,ijts_sqcbsub,ijts_tau,ijts_sqex,ijts_sqsc,
!$OMP*   IJ_CLR_SRINCG,ijts_sqcb,ij_optdw,ij_wtrcld,ij_optdi,ij_icecld,
!$OMP*   IJ_CLR_SRNFG,IJ_CLR_TRDNG,IJ_CLR_SRUPTOA,IJ_CLR_TRUPTOA,
!$OMP*   IJ_CLR_SRNTP,IJ_CLR_TRNTP,IJ_SRNTP,IJ_TRNTP,J_CLRTOA,J_CLRTRP,
!$OMP*   IJ_CLDTPT,IJ_CLDCV1,IJ_CLDT1T,IJ_CLDT1P,IJ_CLDTPPR,J_TOTTRP)
!$OMP    DO SCHEDULE(dynamic,8)
!!OMP*   REDUCTION(+:ICKERR,JCKERR,KCKERR)

C****
C**** MAIN I LOOP
C****
      DO I=I_0,IMAXJ(J)
C**** Radiation input files use a 72x46 grid independent of IM and JM
C**** (ilon,jlat) is the 4x5 box containing the center of box (i,j)
c      JLAT=INT(1.+(J-1.)*45./(JM-1.)+.5)  !  lat_index w.r.to 72x46 grid
c      JLAT=INT(1.+(J-1.)*0.25*DLAT_DG+.5) ! slightly more general
c      ILON=INT(.5+(I-.5)*72./IM+.5)  ! lon_index w.r.to 72x46 grid
      ilon = 1 + int( 72d0*lon2d(i,j)/twopi )
      jlat = 1 + int( 45d0*(lat2d(i,j)+92d0*radian)/pi )
#ifdef ALTER_RADF_BY_LAT
      FULGAS(:)=FULGAS_orig(:)*FULGAS_lat(:,JLAT)
      FS8OPX(:)=FS8OPX_orig(:)*FS8OPX_lat(:,JLAT)
      FT8OPX(:)=FT8OPX_orig(:)*FT8OPX_lat(:,JLAT)
#endif
      use_o3_ref = 0
#ifdef RAD_O3_GCM_HRES
        O3natL(:)=O3JDAY_native(:,I,J)
        O3natLref(:)=O3JREF_native(:,I,J)
#endif

      L1 = 1                         ! lowest layer above ground
      LMR=LM+LM_REQ                  ! radiation allows var. # of layers
      JR=JREG(I,J)

      if (kradia.gt.0) then    ! read in all rad. input data (frc.runs)
        iend = 0
    5   it = itime-1           ! make sure, at least 1 record is read
        do while (mod(itime-it,NDAY*JDPERY).ne.0)
          read(i*1000+j,end=10,err=10) it
C****   input data:          WARNINGS
C****        1 - any changes here also go in later (look for 'iu_rad')
C****        2 - keep "dimrad_sv" up-to-date:         dimrad_sv={
     *     ,Tlm,Tsl,Tgo,Tgoi,Tgli,Tge                  ! LMR+5+
     *     ,shl,wearth,P(i,j),tauwc,tauic,sizeic        ! LMR+1+1+3*LMR+
     *     ,poice,zoice,fsoot ! xdalbs=-dalbsnX*fsoot  ! 1+1+1+
     *     ,wmag,snowoi,snowli,snowe                   ! 1+1+1+1+
     *     ,AgeSn,fmp,flags,ls1_loc                    ! 3+1+.5+.5+
     *     ,snow_frac,zlake,plake                      ! 2+1+1+
     *     ,pocean,plice,pearth,pvt   ! for convenience  1+1+1+12+ 
C****   output data: really needed only if kradia>1
     *     ,srhra,trhra                                ! 2*LMR }
C****   total: dimrad_sv= (7*(LM+LM_REQ) + 38)
     *     ,iy
          if (qcheck) then
            write(out_line,*) 'reading RADfile at Itime',Itime,it,iy
            call write_parallel(trim(out_line),unit=6)
          endif
        end do
        if (it.ne.iy) then
          write(out_line,*) 'RAD input file bad:',itime,it,iy,iend
          call write_parallel(trim(out_line),unit=6)
          call stop_model('RADIA: input file bad',255)
        end if
        go to 20
   10   rewind (i*1000+j)
        iend = iend + 1
        if(iend == 1) go to 5
        call stop_model('RADIA: input file too short',255)
   20   continue
C****   kradia=1: instantaneous forcing - LS1_loc is not used
C****   kradia>1: adjusted forcing, i.e. T adjusts in L=LS1_loc->LM+3
        if(kradia>1) LS1_loc=LS1_loc+3-kradia     ! favorite:kradia=3
        if(kradia>3) LS1_loc=1                    ! favorite:kradia=3
C****   Find arrays derived from P(i,j) : plb,pmidl,aml
        call calc_vert_amp(p(i,j),lmr,pcoll,aml,pdsigl,plb,pmidl)
!****   find more derived fields and do some optional scalings
        RHL(:) = shl(:)/QSAT(TLm(:),LHE,PMIDL(:))
        if(RHfix.ge.0.) RHL(:)=RHfix
        tauwc = cldx*tauwc
        tauic = cldx*tauic
        sizewc = sizeic
        if (kradia.gt.1) then
          do l=1,lm+lm_req
            tlm(l) = tlm(l) + Tchg(l,i,j)
            AFLX_ST(L,I,J,5)=AFLX_ST(L,I,J,5)+Tchg(L,I,J)
          end do
        end if
        COSZ=COSZA(I,J)
        zsnwoi=snowoi/rhos
        dALBsn = -dalbsnX*fsoot*depobc(ilon,jlat)
        zmp=min(0.8d0*fmp,0.9d0*zoice)
        FSTOPX(:)=1. ; FTTOPX(:)=1. ; FTAUC=1. ! deflt (aeros/clouds on)
        kdeliq=0   ! initialize mainly for l>lm
      else ! full model only (kradia not >0)

C**** DETERMINE FRACTIONS FOR SURFACE TYPES AND COLUMN PRESSURE
      PLAND=FLAND(I,J)
      POICE=RSI(I,J)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      PLAKE=FLAKE(I,J)
      PLICE=FLICE(I,J)
      PEARTH=FEARTH(I,J)
      ptype4(1) = pocean ! open ocean and open lake
      ptype4(2) = poice  ! ocean/lake ice
      ptype4(3) = plice  ! glacial ice
      ptype4(4) = pearth ! non glacial ice covered soil

C**** CHECK SURFACE TEMPERATURES
      DO IT=1,4
        IF(ptype4(IT) > 0.) then
          IF(GTEMPR(IT,I,J).LT.124..OR.GTEMPR(IT,I,J).GT.370.) then
            WRITE(6,*) 'In Radia: Time,I,J,IT,TG1',ITime,I,J,IT
     *         ,GTEMPR(IT,I,J)
CCC         STOP 'In Radia: Grnd Temp out of range'
c           ICKERR=ICKERR+1
          END IF
        END IF
      END DO

#if (defined CHL_from_SeaWIFs) || (defined TRACERS_OceanBiology)
C**** Set Chlorophyll concentration
      if (POCEAN.gt.0) then
          LOC_CHL = chl(I,J)
          AIJ(I,J,IJ_CHL)=AIJ(I,J,IJ_CHL)+CHL(I,J)*FOCEAN(I,J)
!         write(*,'(a,3i5,e12.4)')'RAD_DRV:',
!    .    itime,i,j,chl(i,j)
      endif
#endif

      LS1_loc=LTROPO(I,J)+1  ! define stratosphere for radiation
      kdeliq=0   ! initialize mainly for L>LM
C****
C**** DETERMINE CLOUDS (AND THEIR OPTICAL DEPTHS) SEEN BY RADIATION
C****
      CSS=0. ; CMC=0. ; CLDCV=0. ; DEPTH=0. ; OPTDW=0. ; OPTDI=0.
      call dCDNC_EST(ilon,jlat,pland, dCDNC, table)
      dCC_CDNCL = CC_cdncx*dCDNC*CDNCL
      dOD_CDNCL = OD_cdncx*dCDNC*CDNCL
      DO L=1,LM
        if(q(i,j,l)<0) then
           WRITE(6,*)'In Radia: Time,I,J,L,Q<0',ITime,I,J,L,Q,'->0'
           Q(I,J,L)=0.
        end if
        QSS=Q(I,J,L)/(RHSAV(L,I,J)+1.D-20)
        shl(L)=QSS
        IF(FSS(L,I,J)*CLDSAV(L,I,J).LT.1.)
     *       shl(L)=(Q(I,J,L)-QSS*FSS(L,I,J)*CLDSAV(L,I,J))/
     /              (1.-FSS(L,I,J)*CLDSAV(L,I,J))
        TLm(L)=T(I,J,L)*PK(L,I,J)
        TAUSSL=0.
        TAUMCL=0.
        TAUWC(L)=0.
        TAUIC(L)=0.
        SIZEWC(L)=0.
        SIZEIC(L)=0.
        TOTCLD(L)=0.
C**** Determine large scale and moist convective cloud cover for radia
        IF (CLDSS(L,I,J)*(1.+dcc_cdncl(l)).GT.RDSS(L,I,J)) THEN
          TAUSSL=TAUSS(L,I,J)*(1.+dod_cdncl(l))
          shl(L)=QSS
          CSS=1.
          call inc_ajl(i,j,l,jl_sscld,css)
        END IF
        IF (CLDMC(L,I,J).GT.RDMC(I,J)) THEN
          CMC=1.
          call inc_ajl(i,j,l,jl_mccld,cmc)
          DEPTH=DEPTH+PDSIG(L,I,J)
          IF(TAUMC(L,I,J).GT.TAUSSL) THEN
            TAUMCL=TAUMC(L,I,J)
            ELHX=LHE
            IF(TLm(L).LE.TF) ELHX=LHS
            shl(L)=QSAT(TLm(L),ELHX,PMID(L,I,J))
          END IF
        END IF
        IF(TAUSSL+TAUMCL.GT.0.) THEN
             CLDCV=1.
          TOTCLD(L)=1.
          call inc_ajl(i,j,l,jl_totcld,1d0)
C**** save 3D cloud fraction as seen by radiation
          if(cldx>0) AIJL(I,J,L,IJL_CF)=AIJL(I,J,L,IJL_CF)+1.
          IF(TAUMCL.GT.TAUSSL) THEN
            SIZEWC(L)=CSIZMC(L,I,J)
            SIZEIC(L)=CSIZMC(L,I,J)
            IF(SVLAT(L,I,J).EQ.LHE) THEN
              TAUWC(L)=cldx*TAUMCL
              OPTDW=OPTDW+TAUWC(L)
              call inc_ajl(i,j,l,jl_wcld,1d0)
              call inc_ajl(i,j,l,jl_wcldwt,pdsig(l,i,j))
            ELSE
              TAUIC(L)=cldx*TAUMCL
              OPTDI=OPTDI+TAUIC(L)
              call inc_ajl(i,j,l,jl_icld,1d0)
              call inc_ajl(i,j,l,jl_icldwt,pdsig(l,i,j))
            END IF
          ELSE
            SIZEWC(L)=CSIZSS(L,I,J)
            SIZEIC(L)=CSIZSS(L,I,J)
            IF(SVLHX(L,I,J).EQ.LHE) THEN
              TAUWC(L)=cldx*TAUSSL
              OPTDW=OPTDW+TAUWC(L)
              call inc_ajl(i,j,l,jl_wcld,1d0)
              call inc_ajl(i,j,l,jl_wcldwt,pdsig(l,i,j))
            ELSE
              TAUIC(L)=cldx*TAUSSL
              OPTDI=OPTDI+TAUIC(L)
              call inc_ajl(i,j,l,jl_icld,1d0)
              call inc_ajl(i,j,l,jl_icldwt,pdsig(l,i,j))
            END IF
          END IF
          call inc_ajl(i,j,l,jl_wcod,tauwc(l))
          call inc_ajl(i,j,l,jl_icod,tauic(l))
          call inc_ajl(i,j,l,jl_wcsiz,sizewc(l)*tauwc(l))
          call inc_ajl(i,j,l,jl_icsiz,sizeic(l)*tauic(l))
        END IF
C**** save some radiation/cloud fields for wider use
        RCLD(L,I,J)=TAUWC(L)+TAUIC(L)
      END DO
      CFRAC(I,J) = CLDCV    ! cloud fraction consistent with radiation
C**** effective cloud cover diagnostics
         OPNSKY=1.-CLDCV
         DO IT=1,NTYPE
           call inc_aj(i,j,it,J_PCLDSS,CSS  *FTYPE(IT,I,J))
           call inc_aj(i,j,it,J_PCLDMC,CMC  *FTYPE(IT,I,J))
           call inc_aj(i,j,it,J_CLDDEP,DEPTH*FTYPE(IT,I,J))
           call inc_aj(i,j,it,J_PCLD  ,CLDCV*FTYPE(IT,I,J))
         END DO
         call inc_areg(i,j,jr,J_PCLDSS,CSS  )
         call inc_areg(i,j,jr,J_PCLDMC,CMC  )
         call inc_areg(i,j,jr,J_CLDDEP,DEPTH)
         call inc_areg(i,j,jr,J_PCLD  ,CLDCV)
         AIJ(I,J,IJ_PMCCLD)=AIJ(I,J,IJ_PMCCLD)+CMC
         AIJ(I,J,IJ_CLDCV) =AIJ(I,J,IJ_CLDCV) +CLDCV
         DO L=1,LLOW
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDL)=AIJ(I,J,IJ_PCLDL)+1.
           exit
         end do
         DO L=LLOW+1,LMID
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDM)=AIJ(I,J,IJ_PCLDM)+1.
           exit
         end do
         DO L=LMID+1,LHI
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDH)=AIJ(I,J,IJ_PCLDH)+1.
           exit
         end do

         TAUSUMW(I,J) = OPTDW
         TAUSUMI(I,J) = OPTDI
         if(optdw.gt.0.) then
            AIJ(I,J,IJ_optdw)=AIJ(I,J,IJ_optdw)+optdw
            AIJ(I,J,IJ_wtrcld)=AIJ(I,J,IJ_wtrcld)+1.
         end if
         if(optdi.gt.0.) then
            AIJ(I,J,IJ_optdi)=AIJ(I,J,IJ_optdi)+optdi
            AIJ(I,J,IJ_icecld)=AIJ(I,J,IJ_icecld)+1.
         end if

         DO KR=1,NDIUPT
           IF (I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
C**** Warning: this replication may give inaccurate results for hours
C****          1->(NRAD-1)*DTsrc (ADIURN) or skip them (HDIURN)
             TMP(IDD_CL7:IDD_CL7+6)=TOTCLD(1:7)
             TMP(IDD_CCV)=CLDCV
             DO INCH=1,NRAD
               IHM=1+(JTIME+INCH-1)*HR_IN_DAY/NDAY
               IH=IHM
               IF(IH.GT.HR_IN_DAY) IH = IH - HR_IN_DAY
               ADIURN(IDX(:),KR,IH)=ADIURN(IDX(:),KR,IH)+TMP(IDX(:))
#ifndef NO_HDIURN
               IHM = IHM+(JDATE-1)*HR_IN_DAY
               IF(IHM.LE.HR_IN_MONTH) THEN
                 HDIURN(IDX(:),KR,IHM)=HDIURN(IDX(:),KR,IHM)+TMP(IDX(:))
               ENDIF
#endif
             ENDDO
           END IF
         END DO
C****
C**** SET UP VERTICAL ARRAYS OMITTING THE I AND J INDICES
C****
C**** EVEN PRESSURES
      PLB(LM+1)=PEDN(LM+1,I,J)
      DO L=1,LM
        PLB(L)=PEDN(L,I,J)
C**** TEMPERATURES
C---- TLm(L)=T(I,J,L)*PK(L,I,J)     ! already defined
        IF(TLm(L).LT.124..OR.TLm(L).GT.370.) THEN
          WRITE(6,*) 'In Radia: Time,I,J,L,TL',ITime,I,J,L,TLm(L)
          WRITE(6,*) 'GTEMPR:',GTEMPR(:,I,J)
CCC       STOP 'In Radia: Temperature out of range'
c         ICKERR=ICKERR+1
        END IF
C**** MOISTURE VARIABLES
C---- shl(L)=Q(I,J,L)        ! already defined and reset to 0 if <0
c       if(shl(l).lt.0.) then
c         WRITE(0,*)'In Radia: Time,I,J,L,QL<0',ITime,I,J,L,shl(L),'->0'
c         KCKERR=KCKERR+1
c         shl(l)=0.
c       end if
        RHL(L) = shl(L)/QSAT(TLm(L),LHE,PMID(L,I,J))
        if(RHfix.ge.0.) RHL(L)=RHfix
C**** Extra aerosol data
C**** For up to NTRACE aerosols, define the aerosol amount to
C**** be used (kg/m^2)
C**** Only define TRACER is individual tracer is actually defined.
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM) ||\
    (defined TRACERS_OM_SP)
C**** loop over tracers that are passed to radiation.
C**** Some special cases for black carbon, organic carbon, SOAs where
C**** more than one tracer is lumped together for radiation purposes
      do n=1,NTRACE
        if (NTRIX(n).gt.0) then
          select case (trname(NTRIX(n)))
          case ("OCIA")
           TRACER(L,n)=(trm(i,j,l,n_OCII)+trm(i,j,l,n_OCIA))*BYAXYP(I,J)
          case ("OCB") 
           TRACER(L,n)=trm(i,j,l,n_OCB)*BYAXYP(I,J)
           ! The only reason this is a special case is, following
           ! how OCIA is done, it has no wttr factor like default does.
#ifdef TRACERS_AEROSOLS_SOA
          case ("isopp1a")
           TRACER(L,n)=( trm(i,j,l,n_isopp1a)+trm(i,j,l,n_isopp2a)
#ifdef TRACERS_TERP
     &                  +trm(i,j,l,n_apinp1a)+trm(i,j,l,n_apinp2a)
#endif /* TRACERS_TERP */
     &                 )*BYAXYP(I,J)
#endif /* TRACERS_AEROSOLS_SOA */
          case ("OCA4")
           TRACER(L,n)=(trm(i,j,l,n_OCI1)+trm(i,j,l,n_OCA1)+
     *          trm(i,j,l,n_OCI2)+trm(i,j,l,n_OCA2)+
     *          trm(i,j,l,n_OCI3)+trm(i,j,l,n_OCA3)+
     *          trm(i,j,l,n_OCA4))*BYAXYP(I,J)
          case ("BCIA")
           TRACER(L,n)=(trm(i,j,l,n_BCII)+trm(i,j,l,n_BCIA))*BYAXYP(I,J)
          case default
#ifdef TRACERS_NITRATE
! assume full neutralization of NO3p, if NH4 suffice
           select case (trname(NTRIX(n)))
           case ("NO3p")
            if (trm(i,j,l,NTRIX(n)) > 0.d0) then
              nh4_on_no3=min(trm(i,j,l,n_NO3p)*(tr_mm(n_NO3p)+
     *            tr_mm(n_NH4))/tr_mm(n_NO3p)-trm(i,j,l,n_NO3p),
     *                       trm(i,j,l,n_NH4))
              wttr(n)=(nh4_on_no3+trm(i,j,l,NTRIX(n)))/
     *                trm(i,j,l,NTRIX(n))
            endif
           case ("SO4")
            if (trm(i,j,l,NTRIX(n)) > 0.d0) then
              nh4_on_no3=min(trm(i,j,l,n_NO3p)*(tr_mm(n_NO3p)+
     *            tr_mm(n_NH4))/tr_mm(n_NO3p)-trm(i,j,l,n_NO3p),
     *                       trm(i,j,l,n_NH4))
              wttr(n)=(trm(i,j,l,n_NH4)-nh4_on_no3+trm(i,j,l,NTRIX(n)))/
     *                trm(i,j,l,NTRIX(n))
            endif
           end select
#endif
           TRACER(L,n)=wttr(n)*trm(i,j,l,NTRIX(n))*BYAXYP(I,J)
          end select
        end if
      end do
#endif /* TRACERS_AEROSOLS_Koch/DUST/MINERALS/QUARZHEM/OM_SP */

#ifdef TRACERS_AMP
      CALL SETAMP_LEV(i,j,l)
#endif

      END DO
C**** Radiative Equilibrium Layer data
      DO K=1,LM_REQ
        IF(RQT(K,I,J).LT.124..OR.RQT(K,I,J).GT.370.) THEN
        WRITE(6,*) 'In RADIA: Time,I,J,L,TL',ITime,I,J,LM+K,RQT(K,I,J)
CCC     STOP 'In Radia: RQT out of range'
c       JCKERR=JCKERR+1
        END IF
        TLm(LM+K)=RQT(K,I,J)
        PLB(LM+k+1) = PLB0(k)
        shl(LM+k)    = shl0(k)
        RHL(LM+k) = shl(LM+k)/QSAT(TLm(LM+k),LHE,
     *                                .5d0*(PLB(LM+k)+PLB(LM+k+1)) )
        tauwc(LM+k) = 0.
        tauic(LM+k) = 0.
        sizewc(LM+k)= 0.
        sizeic(LM+k)= 0.
#ifdef TRACERS_ON
C**** set radiative equilibirum extra tracer amount to zero
        IF (NTRACE.gt.0) TRACER(LM+k,1:NTRACE)=0.
#endif
      END DO
C**** Zenith angle and GROUND/SURFACE parameters
      COSZ=COSZA(I,J)
      TGO =GTEMPR(1,I,J)
      TGOI=GTEMPR(2,I,J)
      TGLI=GTEMPR(3,I,J)
      TGE =GTEMPR(4,I,J)
      TSL=TSAVG(I,J)
      SNOWOI=SNOWI(I,J)
      SNOWLI=SNOWLI_COM(I,J)
      SNOWE=SNOWE_COM(I,J)                    ! snow depth (kg/m**2)
      snow_frac(:) = fr_snow_rad_ij(:,i,j)    ! snow cover (1)
      AGESN(1)=SNOAGE(3,I,J)    ! land         ! ? why are these numbers
      AGESN(2)=SNOAGE(1,I,J)    ! ocean ice        so confusing ?
      AGESN(3)=SNOAGE(2,I,J)    ! land ice
c      print*,"snowage",i,j,SNOAGE(1,I,J)
C**** set up parameters for new sea ice and snow albedo
      zsnwoi=snowoi/rhos
      dALBsn = xdalbs*depobc(ilon,jlat)
c to use on-line tracer albedo impact, set dALBsnX=0. in rundeck
#if (defined TRACERS_AEROSOLS_Koch) && (defined BC_ALB)
      call GET_BC_DALBEDO(i,j,dALBsn1)
      dALBsn=dALBsn1
#endif
#if (defined TRACERS_AMP) && (defined BC_ALB)
      call GET_BC_DALBEDO(i,j,dALBsn1)
      dALBsn=dALBsn1
#endif
      if (poice.gt.0.) then
        zoice=(ace1i+msi(i,j))/rhoi
        flags=flag_dsws(i,j)
        fmp=min(1.6d0*sqrt(pond_melt(i,j)/rhow),1d0)
           AIJ(I,J,IJ_FRMP) = AIJ(I,J,IJ_FRMP) + fmp*POICE
        zmp=min(0.8d0*fmp,0.9d0*zoice)
      else
        zoice=0. ; flags=.FALSE. ; fmp=0. ; zmp=0.
      endif
C**** set up new lake depth parameter to incr. albedo for shallow lakes
      zlake=0.
      if (plake.gt.0) then
        zlake = MWL(I,J)/(RHOW*PLAKE*AXYP(I,J))
      end if
C****
        !WEARTH=(WEARTH_COM(I,J)+AIEARTH(I,J))/(WFCS(I,J)+1.D-20)
        WEARTH=bare_soil_wetness(i,j)
        if (wearth.gt.1.) wearth=1.
#ifdef USE_ENT
      if ( fearth(i,j) > 0.d0 ) then
        call ent_get_exports( entcells(i,j),
     &       vegetation_fractions=PVT0 )
        call map_ent2giss(PVT0,PVT) !temp hack: ent pfts->giss veg
      else
        PVT(:) = 0.d0  ! actually PVT is not supposed to be used in this case
      endif
#else
      DO K=1,12
        PVT(K)=VDATA(I,J,K)
      END DO
#endif
      WMAG=WSAVG(I,J)
C****
C**** Radiative interaction and forcing diagnostics:
C**** If no radiatively active tracers are defined, nothing changes.
C**** Currently this works for aerosols and ozone but should be extended
C**** to cope with all trace gases.
C****
      FTAUC=1. ! deflt (clouds on)
      use_tracer_chem(:) = 0 ! by default use climatological ozone/ch4
C**** Set level for inst. rad. forc. calcs for aerosols/trace gases
C**** This is set from the rundeck.
      LFRC=LM+LM_REQ+1          ! TOA
      if (rad_forc_lev.gt.0) LFRC=LTROPO(I,J) ! TROPOPAUSE
#ifdef ACCMIP_LIKE_DIAGS
      if(rad_forc_lev.gt.0)call stop_model
     &('ACCMIP_LIKE_DIAGS desires TOA RADF diags',255)
#endif
C**** The calculation of the forcing is slightly different.
C**** depending on whether full radiative interaction is turned on
C**** or not.
      onoff_aer=0; onoff_chem=0
      if (rad_interact_aer > 0) onoff_aer=1
      if (clim_interact_chem > 0) onoff_chem=1
      FSTOPX(:)=onoff_aer ; FTTOPX(:)=onoff_aer
      use_o3_ref=0

#ifdef TRACERS_SPECIAL_Shindell
C**** Ozone and Methane: 
      CHEM_IN(1:2,1:LM)=chem_tracer_save(1:2,1:LM,I,J)
      if (clim_interact_chem > 0) then
        use_tracer_chem(1)=Lmax_rad_O3  ! O3
        use_tracer_chem(2)=Lmax_rad_CH4 ! CH4
      endif
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      if(clim_interact_chem<=0)
     &call stop_model("stratOx RADF on, clim_interact_chem<=0",255)
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_SPECIAL_Shindell */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM) ||\
    (defined TRACERS_OM_SP)
        
C**** Aerosols incl. Dust:        set up for radiative forcing diagnostics
      if (NTRACE>0 .and. moddrf==0) then
        set_clayilli=.FALSE.
        set_claykaol=.FALSE.
        set_claysmec=.FALSE.
        set_claycalc=.FALSE.
        set_clayquar=.FALSE.
        do n=1,NTRACE
          IF (ntrix(n) > 0) THEN
            IF (trname(NTRIX(n)).eq."seasalt2") CYCLE ! not for seasalt2
            IF (trname(ntrix(n)) == 'ClayIlli' .AND. set_clayilli) cycle
            IF (trname(ntrix(n)) == 'ClayKaol' .AND. set_claykaol) cycle
            IF (trname(ntrix(n)) == 'ClaySmec' .AND. set_claysmec) cycle
            IF (trname(ntrix(n)) == 'ClayCalc' .AND. set_claycalc) cycle
            IF (trname(ntrix(n)) == 'ClayQuar' .AND. set_clayquar) cycle
            FSTOPX(n)=1-onoff_aer ; FTTOPX(n)=1-onoff_aer ! turn on/off tracer
C**** Warning: small bit of hardcoding assumes that seasalt2 immediately
C****          succeeds seasalt1 in NTRACE array
            IF (trname(NTRIX(n)).eq."seasalt1") THEN          !add seasalt2
              FSTOPX(n+1)=1-onoff_aer;FTTOPX(n+1)=1-onoff_aer !to seasalt1
            END IF
C**** Do radiation calculations for all clay classes at once
C**** Assumes that 4 clay tracers are adjacent in NTRACE array
            SELECT CASE (trname(ntrix(n)))
            CASE ('ClayIlli')
              fstopx(n+1:n+3)=1-onoff_aer; fttopx(n+1:n+3)=1-onoff_aer
              set_clayilli=.true.
            CASE ('ClayKaol')
              fstopx(n+1:n+3)=1-onoff_aer; fttopx(n+1:n+3)=1-onoff_aer
              set_claykaol=.true.
            CASE ('ClaySmec')
              fstopx(n+1:n+3)=1-onoff_aer; fttopx(n+1:n+3)=1-onoff_aer
              set_claysmec=.true.
            CASE ('ClayCalc')
              fstopx(n+1:n+3)=1-onoff_aer; fttopx(n+1:n+3)=1-onoff_aer
              set_claycalc=.true.
            CASE ('ClayQuar')
              fstopx(n+1:n+3)=1-onoff_aer; fttopx(n+1:n+3)=1-onoff_aer
              set_clayquar=.true.
            END SELECT
            kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
            CALL RCOMPX  ! tr.aero.Koch/dust/miner./quarz/om_sp
            SNFST(1,n,I,J)=SRNFLB(1) ! surface forcing
            TNFST(1,n,I,J)=TRNFLB(1)
            SNFST(2,n,I,J)=SRNFLB(LFRC)
            TNFST(2,n,I,J)=TRNFLB(LFRC)
            FSTOPX(n)=onoff_aer ; FTTOPX(n)=onoff_aer   ! back to default
            IF (trname(NTRIX(n)).eq."seasalt1") THEN    ! also for seasalt2
              FSTOPX(n+1)=onoff_aer ; FTTOPX(n+1)=onoff_aer
            END IF
            SELECT CASE (trname(ntrix(n)))           ! also for clays
            CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc',
     &           'ClayQuar')
              fstopx(n+1:n+3)=onoff_aer ; fttopx(n+1:n+3)=onoff_aer
            END SELECT
          END IF
        end do
      end if
#endif /* TRACERS_AEROSOLS_Koch/DUST/MINERALS/QUARZHEM/OM_SP */

      if (moddrf==0) then
#ifdef TRACERS_SPECIAL_Shindell
C**** Ozone:
! ozone rad forcing diags now use a constant reference year
! for this first call. And no tracer values...
        use_o3_ref=1 ; use_tracer_chem(1)=0
        kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
        CALL RCOMPX        ! tr_Shindell Ox tracer
        SNFST_o3ref(1,I,J)=SRNFLB(LTROPO(I,J)) ! meteo.tropopause
        TNFST_o3ref(1,I,J)=TRNFLB(LTROPO(I,J))
        SNFST_o3ref(2,I,J)=SRNFLB(LM+LM_REQ+1) ! T.O.A.
        TNFST_o3ref(2,I,J)=TRNFLB(LM+LM_REQ+1)
        SNFST_o3ref(5,I,J)=SRNFLB(LS1-1) ! fixed tropopause
        TNFST_o3ref(5,I,J)=TRNFLB(LS1-1)
#ifdef AUXILIARY_OX_RADF
! if needed, also save the auxiliary ozone field (i.e. climatology
! if tracer is used in final call, tracers if climatology is used.)
#ifdef AUX_OX_RADF_TROP
        ! forces use of tracer from L=1,LS1-1 and reference above that:
        use_o3_ref=1 ; use_tracer_chem(1)=LS1-1 
#else
        ! use tracer or climatology, whichever won't be used in final call:
        use_o3_ref=0 ; use_tracer_chem(1)=(1-onoff_chem)*Lmax_rad_O3
#endif
        kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
        CALL RCOMPX        ! tr_Shindell Ox tracer
#ifdef AUX_OX_RADF_TROP
        SNFST_o3ref(3,I,J)=SRNFLB(LS1-1)       ! fixed tropopause
        TNFST_o3ref(3,I,J)=TRNFLB(LS1-1)
#else
        SNFST_o3ref(3,I,J)=SRNFLB(LTROPO(I,J)) ! meteor.tropopause
        TNFST_o3ref(3,I,J)=TRNFLB(LTROPO(I,J))
#endif
        SNFST_o3ref(4,I,J)=SRNFLB(LM+LM_REQ+1) ! T.O.A.
        TNFST_o3ref(4,I,J)=TRNFLB(LM+LM_REQ+1)
#endif
! ... after which it can use either climatological or tracer O3:
        use_o3_ref=0 ; use_tracer_chem(1)=onoff_chem*Lmax_rad_O3
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
! Optional intermediate call with stratOx tracer:
!NEED chem_IN(1,1:LM)=stratO3_tracer_save(1:LM,I,J)
!NEED kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
!NEED CALL RCOMPX        ! stratOx diag tracer
        SNFST_stratOx(1,I,J)=SRNFLB(LTROPO(I,J)) ! tropopause
        TNFST_stratOx(1,I,J)=TRNFLB(LTROPO(I,J))
        SNFST_stratOx(2,I,J)=SRNFLB(LM+LM_REQ+1) ! T.O.A.
        TNFST_stratOx(2,I,J)=TRNFLB(LM+LM_REQ+1)
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS */
        chem_IN(1,1:LM)=chem_tracer_save(1,1:LM,I,J)  ! Ozone
        chem_IN(2,1:LM)=chem_tracer_save(2,1:LM,I,J)  ! Methane
#ifdef ACCMIP_LIKE_DIAGS
! TOA GHG rad forcing: nf=1,4 are CH4, N2O, CFC11, and CFC12:
! Initial calls are reference year/day:
      do nf=1,4
        if(nf==1)then ! CH4 reference call must not use tracer
          use_tracer_chem(2)=0
        else ! N2O and CFC call's CH4 should match final call
          use_tracer_chem(2)=onoff_chem*Lmax_rad_CH4
        endif
        fulgas(nfghg(nf))=sv_fulgas_ref(nf)
        kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
        CALL RCOMPX
        SNFS_ghg(nf,I,J)=SRNFLB(LM+LM_REQ+1)
        TNFS_ghg(nf,I,J)=TRNFLB(LM+LM_REQ+1)
        fulgas(nfghg(nf))=sv_fulgas_now(nf)
      enddo
#endif /* ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_SPECIAL_Shindell */
      end if ! moddrf=0
#ifdef TRACERS_SPECIAL_Shindell
! final (main) RCOMPX call can use tracer methane (or not):
      use_tracer_chem(2)=onoff_chem*Lmax_rad_CH4
      if (is_set_param('initial_GHG_setup')) then
        call get_param('initial_GHG_setup', initial_GHG_setup)
        if (initial_GHG_setup == 1 .and. itime == itimeI) then
          use_tracer_chem(2)=0  ! special case; model outputs climatology
        end if
      end if
#endif /* TRACERS_SPECIAL_Shindell */


      if (moddrf==0) then
#ifdef BC_ALB
        dalbsn=0.d0
        CALL RCOMPX
        NFSNBC(I,J)=SRNFLB(LM+LM_REQ+1)
c       NFSNBC(I,J)=SRNFLB(LFRC)
        ALBNBC(I,J)=SRNFLB(1)/(SRDFLB(1)+1.D-20)
c set for BC-albedo effect
        dALBsn=dALBsn1
#endif
#ifdef TRACERS_AMP
        IF (AMP_DIAG_FC == 2) THEN
          Do n = 1,nmodes
            FSTOPX(n) = 1-onoff_aer !turns off online tracer
            FTTOPX(n) = 1-onoff_aer !
            if (n.eq.1) FSTOPX(:) = 1-onoff_aer
            if (n.eq.1) FTTOPX(:) = 1-onoff_aer
            CALL RCOMPX
            SNFST(1,n,I,J)=SRNFLB(1) ! surface forcing
            TNFST(1,n,I,J)=TRNFLB(1)
            SNFST(2,n,I,J)=SRNFLB(LFRC) ! Tropopause forcing
            TNFST(2,n,I,J)=TRNFLB(LFRC)
            FSTOPX(:) = onoff_aer !turns on online tracer
            FTTOPX(:) = onoff_aer !
          ENDDO
        ELSE
           n = 1
          FSTOPX(:) = 1-onoff_aer !turns off online tracer
          FTTOPX(:) = 1-onoff_aer !
          CALL RCOMPX
          SNFST(1,n,I,J)=SRNFLB(1) ! surface forcing
          TNFST(1,n,I,J)=TRNFLB(1)
          SNFST(2,n,I,J)=SRNFLB(LFRC) ! Tropopause forcing
          TNFST(2,n,I,J)=TRNFLB(LFRC)
          FSTOPX(:) = onoff_aer !turns on online tracer
          FTTOPX(:) = onoff_aer !
        ENDIF
#endif
C**** Optional calculation of CRF using a clear sky calc.
        if (cloud_rad_forc.gt.0) then
          FTAUC=0.   ! turn off cloud tau (tauic +tauwc)
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
#ifdef SCM
          KEEPAL=0
          if (I.eq.I_TARG.and.J.eq.J_TARG) then
c             if SCM using prescribed surface albedo
              if (SCM_SURF_ALBEDO_FLAG.eq.1) KEEPAL = 1
          endif
#endif
          CALL RCOMPX          ! cloud_rad_forc>0 : clr sky
      
          SNFSCRF(I,J)=SRNFLB(LM+LM_REQ+1)   ! always TOA
          TNFSCRF(I,J)=TRNFLB(LM+LM_REQ+1)   ! always TOA
#ifdef SCM
          KEEPAL = 0
          if (I.eq.I_TARG.and.J.eq.J_TARG) then
              CSSRNTOP = SRNFLB(LM+LM_REQ+1)*COSZ2(I,J)
              CSTRUTOP = TRUFLB(LM+LM_REQ+1)
              CSSRNBOT = SRNFLB(1)*COSZ2(I,J)
              CSTRNBOT = TRNFLB(1)
              CSSRDBOT = SRDFLB(1)*COSZ2(I,J)
c             write(iu_scm_prt,8887) (SRBALB(iscm),iscm=1,6)
c8887         format(//1x,'albedo ',6(f8.4))
c             write(iu_scm_prt,8888) CSSRNTOP,CSTRUTOP,CSSRNBOT,
c    &              CSTRNBOT,CSSRDBOT,COSZ2(I,J)
c8888         format(1x,'clr sky swntop lwutop swnbot lwnbot swdbot',
c    &            5(f10.4),' cosz2 ',f10.4)
          endif
#endif
C         BEGIN AMIP
          AIJ(I,J,IJ_SWDCLS)=AIJ(I,J,IJ_SWDCLS)+SRDFLB(1)*COSZ2(I,J)
          AIJ(I,J,IJ_SWNCLS)=AIJ(I,J,IJ_SWNCLS)+SRNFLB(1)*COSZ2(I,J)
          AIJ(I,J,IJ_LWDCLS)=AIJ(I,J,IJ_LWDCLS)+TRDFLB(1)
          AIJ(I,J,IJ_SWNCLT)=AIJ(I,J,IJ_SWNCLT)+SRNFLB(LM+LM_REQ+1)
     *     *COSZ2(I,J)
          AIJ(I,J,IJ_LWNCLT)=AIJ(I,J,IJ_LWNCLT)+TRNFLB(LM+LM_REQ+1)
C       END AMIP
        end if
        FTAUC=1.     ! default: turn on cloud tau


C**** 2nd Optional calculation of CRF using a clear sky calc. without aerosols and Ox
        if (cloud_rad_forc.gt.0) then
          FTAUC=0.   ! turn off cloud tau (tauic +tauwc)
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
c Including turn off of aerosols and Ox during crf calc.+++++++++++++++++++ 
#ifdef TRACERS_SPECIAL_Shindell
       use_o3_ref=1 ; use_tracer_chem(1)=0  !turns off ozone
#endif
       FSTOPX(:) = 1-onoff_aer !turns off aerosol tracer
       FTTOPX(:) = 1-onoff_aer !
        CALL RCOMPX          ! cloud_rad_forc>0 : clr sky
       FSTOPX(:) = onoff_aer !turns on aerosol tracer
       FTTOPX(:) = onoff_aer !
#ifdef TRACERS_SPECIAL_Shindell
       use_o3_ref=0 ; use_tracer_chem(1)=onoff_chem*Lmax_rad_O3 ! turns on ozone tracers
#endif
          SNFSCRF2(I,J)=SRNFLB(LM+LM_REQ+1)   ! always TOA
          TNFSCRF2(I,J)=TRNFLB(LM+LM_REQ+1)   ! always TOA

c          AIJ(I,J,IJ_SWDCLS2)=AIJ(I,J,IJ_SWDCLS2)+SRDFLB(1)*COSZ2(I,J)
c          AIJ(I,J,IJ_SWNCLS2)=AIJ(I,J,IJ_SWNCLS2)+SRNFLB(1)*COSZ2(I,J)
c          AIJ(I,J,IJ_LWDCLS2)=AIJ(I,J,IJ_LWDCLS2)+TRDFLB(1)
c          AIJ(I,J,IJ_SWNCLT2)=AIJ(I,J,IJ_SWNCLT2)+SRNFLB(LM+LM_REQ+1)
c     *     *COSZ2(I,J)
c          AIJ(I,J,IJ_LWNCLT2)=AIJ(I,J,IJ_LWNCLT2)+TRNFLB(LM+LM_REQ+1)
        end if
        FTAUC=1.     ! default: turn on cloud tau


C**** Optional calculation of the impact of default aerosols
        if (aer_rad_forc.gt.0) then
C**** first, separate aerosols
          DO N=1,8
          tmpS(N)=FS8OPX(N)   ; tmpT(N)=FT8OPX(N)
          FS8OPX(N)=0.     ; FT8OPX(N)=0.
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
          CALL RCOMPX           ! aer_rad_forc>0 : no aerosol N
          SNFSAERRF(N,I,J)=SRNFLB(LM+LM_REQ+1) ! TOA
          TNFSAERRF(N,I,J)=TRNFLB(LM+LM_REQ+1) ! TOA
          SNFSAERRF(N+8,I,J)=SRNFLB(1) ! SURF
          TNFSAERRF(N+8,I,J)=TRNFLB(1) ! SURF
          FS8OPX(N)=tmpS(N)   ; FT8OPX(N)=tmpT(N)
          END DO
C**** second, net aerosols
          tmpS(:)=FS8OPX(:)   ; tmpT(:)=FT8OPX(:)
          FS8OPX(:)=0.     ; FT8OPX(:)=0.
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
          CALL RCOMPX             ! aer_rad_forc>0 : no aerosols
          SNFSAERRF(17,I,J)=SRNFLB(LM+LM_REQ+1) ! TOA
          TNFSAERRF(17,I,J)=TRNFLB(LM+LM_REQ+1) ! TOA
          SNFSAERRF(18,I,J)=SRNFLB(1) ! SURF
          TNFSAERRF(18,I,J)=TRNFLB(1) ! SURF
          FS8OPX(:)=tmpS(:)   ; FT8OPX(:)=tmpT(:)
        end if
      end if  ! moddrf=0
C**** End of initial computations for optional forcing diagnostics
      end if  ! kradia not >0, i.e. full model only
C**** Localize fields that are modified by RCOMPX
      kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)

#ifdef SCM
      KEEPAL=0
      if (I.eq.I_TARG.and.J.eq.J_TARG) then
c         if SCM using prescribed surface albedo
          if (SCM_SURF_ALBEDO_FLAG.eq.1) KEEPAL = 1
      endif
#endif
C*****************************************************
C     Main RADIATIVE computations, SOLAR and THERM(A)L

      CALL RCOMPX
C*****************************************************
#ifdef SCM
      KEEPAL = 0
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM) ||\
    (defined TRACERS_OM_SP)    || (defined TRACERS_AMP)

#ifdef TRACERS_AMP
      NTRACE = nmodes
#endif

C**** Save optical depth diags
      do n=1,NTRACE
        IF (ntrix(n) > 0) THEN
          SELECT CASE (trname(ntrix(n)))
          CASE ('Clay')
            n1=n-nrad_clay+1
            IF (diag_rad /= 1) THEN
              IF (ijts_tausub(1,ntrix(n),n1) > 0)
     &             taijs(i,j,ijts_tausub(1,ntrix(n),n1))
     &             =taijs(i,j,ijts_tausub(1,ntrix(n),n1))
     &             +SUM(ttausv(1:Lm,n))
              IF (ijts_tausub(2,ntrix(n),n1) > 0)
     &             taijs(i,j,ijts_tausub(2,ntrix(n),n1))
     &             =taijs(i,j,ijts_tausub(2,ntrix(n),n1))
     &             +SUM(ttausv(1:Lm,n))*OPNSKY
            END IF
            if (ijlt_3Dtau(NTRIX(n)).gt.0)
     *           taijls(i,j,1:lm,ijlt_3Dtau(NTRIX(n)))
     *           =taijls(i,j,1:lm,ijlt_3Dtau(NTRIX(n)))+TTAUSV(1:lm,n)
            IF (diag_rad == 1) THEN
              DO kr=1,6
                IF (ijts_sqexsub(1,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqexsub(1,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqexsub(1,kr,ntrix(n),n1))
     &               +SUM(aesqex(1:Lm,kr,n))
                IF (ijts_sqexsub(2,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqexsub(2,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqexsub(2,kr,ntrix(n),n1))
     &               +SUM(aesqex(1:Lm,kr,n))*OPNSKY
                IF (ijts_sqscsub(1,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqscsub(1,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqscsub(1,kr,ntrix(n),n1))
     &               +SUM(aesqsc(1:Lm,kr,n))
                IF (ijts_sqscsub(2,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqscsub(2,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqscsub(2,kr,ntrix(n),n1))
     &               +SUM(aesqsc(1:Lm,kr,n))*OPNSKY
                IF (ijts_sqcbsub(1,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqcbsub(1,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqcbsub(1,kr,ntrix(n),n1))
     &               +SUM(aesqcb(1:Lm,kr,n))
     &               /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)
                IF (ijts_sqcbsub(2,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqcbsub(2,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqcbsub(2,kr,ntrix(n),n1))
     &               +SUM(aesqcb(1:Lm,kr,n))
     &               /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)*OPNSKY
              END DO
            END IF
          CASE DEFAULT

            IF (diag_rad /= 1) THEN
              if (ijts_tau(1,NTRIX(n)).gt.0)
     &             taijs(i,j,ijts_tau(1,NTRIX(n)))
     &             =taijs(i,j,ijts_tau(1,NTRIX(n)))+SUM(TTAUSV(1:lm,n))
              if (ijts_tau(2,NTRIX(n)).gt.0)
     &             taijs(i,j,ijts_tau(2,NTRIX(n)))
     &             =taijs(i,j,ijts_tau(2,NTRIX(n)))
     &             +SUM(TTAUSV(1:lm,n))*OPNSKY
            END IF
            if (ijlt_3Dtau(NTRIX(n)).gt.0)
     &           taijls(i,j,1:lm,ijlt_3Dtau(NTRIX(n)))
     &           =taijls(i,j,1:lm,ijlt_3Dtau(NTRIX(n)))+TTAUSV(1:lm,n)
            IF (diag_rad == 1) THEN
              DO kr=1,6
c                 print*,'SUSA  diag',SUM(aesqex(1:Lm,kr,n))
                IF (ijts_sqex(1,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqex(1,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqex(1,kr,ntrix(n)))
     &               +SUM(aesqex(1:Lm,kr,n))
                IF (ijts_sqex(2,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqex(2,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqex(2,kr,ntrix(n)))
     &               +SUM(aesqex(1:Lm,kr,n))*OPNSKY
                IF (ijts_sqsc(1,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqsc(1,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqsc(1,kr,ntrix(n)))
     &               +SUM(aesqsc(1:Lm,kr,n))
                IF (ijts_sqsc(2,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqsc(2,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqsc(2,kr,ntrix(n)))
     &               +SUM(aesqsc(1:Lm,kr,n))*OPNSKY
                IF (ijts_sqcb(1,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqcb(1,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqcb(1,kr,ntrix(n)))
     &               +SUM(aesqcb(1:Lm,kr,n))
     &               /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)
                IF (ijts_sqcb(2,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqcb(2,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqcb(2,kr,ntrix(n)))
     &               +SUM(aesqcb(1:Lm,kr,n))
     &               /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)*OPNSKY
              END DO
            END IF
          END SELECT
        END IF
      end do

      DO KR=1,NDIUPT
        IF (I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
          TMP(idd_aot) =SUM(aesqex(1:Lm,6,1:NTRACE)) !*OPNSKY
          DO INCH=1,NRAD
            IHM=1+(JTIME+INCH-1)*HR_IN_DAY/NDAY
            IH=IHM
            IF(IH.GT.HR_IN_DAY) IH = IH - HR_IN_DAY
            ADIURN(idd_aot,KR,IH)=ADIURN(idd_aot,KR,IH)+TMP(idd_aot)
#ifndef NO_HDIURN
            IHM = IHM+(JDATE-1)*HR_IN_DAY
            IF(IHM.LE.HR_IN_MONTH) THEN
              HDIURN(idd_aot,KR,IHM)=HDIURN(idd_aot,KR,IHM)+TMP(idd_aot)
            ENDIF
#endif
          ENDDO
        ENDIF
      ENDDO

#endif

#ifdef TRACERS_ON
! this to accumulate daily SUM of optical thickness for
! each active tracer. Will become average in DIAG.f.
! Also saving the aerosol absorption (band 6) 3D, summed over species. 
      aerAbs6SaveInst(i,j,:)=0.d0
      do n=1,NTRACE
        if(ntrix(n) > 0) then
          StauL=sum(ttausv(1:LM,n))
          ttausv_sum(i,j,ntrix(n))=ttausv_sum(i,j,ntrix(n))+StauL
          ttausv_sum_cs(i,j,ntrix(n))=ttausv_sum_cs(i,j,ntrix(n))
     &         +StauL*OPNSKY
          aerAbs6SaveInst(i,j,1:lm)=aerAbs6SaveInst(i,j,1:lm) +
     &    (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))
        endif
      enddo
      IF (adiurn_dust == 1 .or. save3dAOD == 1) THEN
        ttausv_save(i,j,:,:)=0.D0
        DO n=1,NTRACE
          IF (ntrix(n) > 0) THEN
            do k=1,LM
              ttausv_save(i,j,ntrix(n),k)=ttausv_save(i,j,ntrix(n),k)
     &             +ttausv(k,n)
              ttausv_cs_save(i,j,ntrix(n),k)
     &             =ttausv_cs_save(i,j,ntrix(n),k)+ttausv(k,n)*OPNSKY
            end do
          END IF
        END DO
      END IF
#ifdef TRACERS_SPECIAL_Shindell
      do n=1,NTRACE
        if(ntrix(n) > 0) then
          do k=1,LM
            ttausv_ntrace(i,j,n,k)=ttausv(k,n)
          end do
        end if
      end do
#endif
#endif /* TRACERS_ON */

      IF (I.EQ.IWRITE .and. J.EQ.JWRITE) CALL WRITER(6,ITWRITE)
      CSZ2=COSZ2(I,J)
      do L=1,LM
        rad_to_chem(:,L,i,j)=chem_out(L,:)
        do k=1,4
          kliq(L,k,i,j)=kdeliq(L,k) ! save updated flags
        end do
      end do
      if (kradia.gt.0) then  ! rad. forc. model; acc diagn
        do L=1,LM+LM_REQ+1
          AFLX_ST(L,I,J,1)=AFLX_ST(L,I,J,1)+SRUFLB(L)*CSZ2
          AFLX_ST(L,I,J,2)=AFLX_ST(L,I,J,2)+SRDFLB(L)*CSZ2
          AFLX_ST(L,I,J,3)=AFLX_ST(L,I,J,3)+TRUFLB(L)
          AFLX_ST(L,I,J,4)=AFLX_ST(L,I,J,4)+TRDFLB(L)
        end do
        if(kradia.eq.1) then
          tauex6=0. ; tauex5=0. ; tausct=0. ; taugcb=0.
          do L=1,LM
            AFLX_ST(L,I,J,5)=AFLX_ST(L,I,J,5)+1.d2*RHL(L)
            tauex6=tauex6+SRAEXT(L,6)+SRDEXT(L,6)+SRVEXT(L,6)
            tauex5=tauex5+SRAEXT(L,5)+SRDEXT(L,5)+SRVEXT(L,5)
            tausct=tausct+SRASCT(L,6)+SRDSCT(L,6)+SRVSCT(L,6)
            taugcb=taugcb+SRASCT(L,6)*SRAGCB(L,6)+
     +        SRDSCT(L,6)*SRDGCB(L,6)+SRVSCT(L,6)*SRVGCB(L,6)
          end do
          AFLX_ST(LM+1,I,J,5)=AFLX_ST(LM+1,I,J,5)+tauex5
          AFLX_ST(LM+2,I,J,5)=AFLX_ST(LM+2,I,J,5)+tauex6
          AFLX_ST(LM+3,I,J,5)=AFLX_ST(LM+3,I,J,5)+tausct
          AFLX_ST(LM+4,I,J,5)=AFLX_ST(LM+4,I,J,5)+taugcb
          cycle
        end if
        do l=LS1_loc,lmr
          tchg(l,i,j) = tchg(l,i,j) + ( srfhrl(l)*csz2-srhra(l) +
     +      (-trfcrl(l)-trhra(l)) )*nrad*DTsrc*bysha/aml(l)
        end do
        cycle
      else if (kradia.lt.0) then ! save i/o data for frc.runs
        do L=1,LM+LM_REQ                 ! output data (for adj frc)
          SRHRA(L)=SRFHRL(L)*CSZ2
          TRHRA(L)=-TRFCRL(L)
        end do
C****   Save all input data to disk if kradia<0
        if(ijradfrc(i,j)>0) write(ijradfrc(i,j)) itime
     &     ,Tlm,tsl,tgo,tgoi,tgli,tge                  ! LMR+5+
     &     ,shl,wearth,P(i,j),tauwc,tauic,sizeic       ! LMR+1+1+3*LMR+
     &     ,poice,zoice,fsoot ! xdalbs=-dalbsnX*fsoot  ! 1+1+1+
     &     ,wmag,snowoi,snowli,snowe                   ! 1+1+1+1+
     &     ,AgeSn,fmp,flags,ls1_loc                    ! 3+1+.5+.5+
     &     ,snow_frac,zlake,plake                      ! 2+1+1+
     *     ,pocean,plice,pearth,pvt   ! for convenience  1+1+1+12+
C**** output data: really needed only if kradia>1
     &     ,srhra,trhra                                ! 2*LMR
     &     ,itime
      end if ! kradia not 0
C****
C**** Save relevant output in model arrays
C****
C**** (some generalisation and coherence needed in the rad surf type calc)
      FSF(1,I,J)=FSRNFG(1)   !  ocean
      FSF(2,I,J)=FSRNFG(3)   !  ocean ice
      FSF(3,I,J)=FSRNFG(4)   !  land ice
      FSF(4,I,J)=FSRNFG(2)   !  soil
      SRHR(0,I,J)=SRNFLB(1)
      TRHR(0,I,J)=STBO*(POCEAN*GTEMPR(1,I,J)**4+POICE*GTEMPR(2,I,J)**4+
     +  PLICE*GTEMPR(3,I,J)**4+PEARTH*GTEMPR(4,I,J)**4)-TRNFLB(1)
      TRSURF(1,I,J) = STBO*GTEMPR(1,I,J)**4  !  ocean
      TRSURF(2,I,J) = STBO*GTEMPR(2,I,J)**4  !  ocean ice
      TRSURF(3,I,J) = STBO*GTEMPR(3,I,J)**4  !  land ice
      TRSURF(4,I,J) = STBO*GTEMPR(4,I,J)**4  !  soil
      DO L=1,LM
        SRHR(L,I,J)=SRFHRL(L)
        TRHR(L,I,J)=-TRFCRL(L)
      END DO
      DO LR=1,LM_REQ
        SRHRS(LR,I,J)= SRFHRL(LM+LR)
        TRHRS(LR,I,J)=-TRFCRL(LM+LR)
      END DO
C**** Save fluxes at four levels surface, P0, P1, LTROPO
      SNFS(1,I,J)=SRNFLB(1)     ! Surface
      TNFS(1,I,J)=TRNFLB(1)
      SNFS(2,I,J)=SRNFLB(LM+1)  ! P1
      TNFS(2,I,J)=TRNFLB(LM+1)
      SNFS(3,I,J)=SRNFLB(LM+LM_REQ+1) ! P0 = TOA
      TNFS(3,I,J)=TRNFLB(LM+LM_REQ+1)
      SNFS(4,I,J)=SRNFLB(LTROPO(I,J)) ! LTROPO
      TNFS(4,I,J)=TRNFLB(LTROPO(I,J))

#ifdef SCM
      if (I.eq.I_TARG .and. J.eq.J_TARG) then
          do L=1,LM
             SRFHRLCOL(L) = SRFHRL(L) * COSZ1(I,J)
             TRFCRLCOL(L) = TRFCRL(L)
          enddo
          SRNFLBBOT = SRNFLB(1) * COSZ1(I,J)            ! Surface
          SRNFLBTOP = SRNFLB(LM+LM_REQ+1) * COSZ1(I,J)  ! P0 = TOA
          SRDFLBBOT = SRDFLB(1) * COSZ1(I,J)
          SRDFLBTOP = SRDFLB(LM+LM_REQ+1) * COSZ1(I,J)
          SRUFLBBOT = SRUFLB(1) * COSZ1(I,J)
          SRUFLBTOP = SRUFLB(LM+LM_REQ+1) * COSZ1(I,J)
          TRUFLBTOP = TRUFLB(LM+LM_REQ+1)
          TRDFLBTOP = TRDFLB(LM+LM_REQ+1)
          TRUFLBBOT = TRUFLB(1)
          TRDFLBBOT = TRDFLB(1)
      endif
#endif

C****
      TRINCG(I,J)=TRDFLB(1)
      BTMPW(I,J)=BTEMPW-TF
      ALB(I,J,1)=SRNFLB(1)/(SRDFLB(1)+1.D-20)
      ALB(I,J,2)=PLAVIS
      ALB(I,J,3)=PLANIR
      ALB(I,J,4)=ALBVIS
      ALB(I,J,5)=ALBNIR
      ALB(I,J,6)=SRRVIS
      ALB(I,J,7)=SRRNIR
      ALB(I,J,8)=SRAVIS
      ALB(I,J,9)=SRANIR

#ifdef TRACERS_DUST
      IF (adiurn_dust == 1) THEN
        srnflb_save(i,j,1:lm)=srnflb(1:lm)
        trnflb_save(i,j,1:lm)=trnflb(1:lm)
      END IF
#endif
#ifdef CACHED_SUBDD
      SWU_ij(I,J)=SRUFLB(1)*CSZ2
#endif

      SRDN(I,J) = SRDFLB(1)     ! save total solar flux at surface
C**** SALB(I,J)=ALB(I,J,1)      ! save surface albedo (pointer)
      FSRDIR(I,J)=SRXVIS        ! direct visible solar at surface **coefficient
      SRVISSURF(I,J)=SRDVIS     ! total visible solar at surface
      FSRDIF(I,J)=SRDVIS*(1-SRXVIS) ! diffuse visible solar at surface

      DIRNIR(I,J)=SRXNIR*SRDNIR     ! direct beam nir solar at surface
      DIFNIR(I,J)=SRDNIR*(1-SRXNIR) ! diffuse     nir solar at surface

cdiag write(*,'(a,2i5,6e12.4)')'RAD_DRV: ',
cdiag.    I,J,FSRDIR(I,J),SRVISSURF(I,J),FSRDIF(I,J),
cdiag.        DIRNIR(I,J),SRDNIR,DIFNIR(I,J)

#ifdef SCM   
      if (I.eq.I_TARG.and.J.eq.J_TARG) then
c         do L=1,LM
c            write(iu_scm_prt,884) L,TOTCLD(L)
c884         format(1x,'in RAD_DRV  L totcld ',i5,f8.2)
c         enddo
c         write(iu_scm_prt,885)  J,I,SRDVIS,SRXVIS,SRDNIR,SRXNIR
c885      format(1x,'in RAD_DRV  J I SRDVIS SRXVIS SRDNIR SRXNIR ',
c    &      i5,i5,2(f10.4,f10.6))
          scm_FDIR_vis = SRXVIS*SRDVIS
          scm_FDIF_vis = (1.-SRXVIS)*SRDVIS
          scm_FDIR_nir = SRXNIR*SRDNIR
          scm_FDIF_nir = (1.-SRXNIR)*SRDNIR
c         write(iu_scm_prt,886) scm_FDIR_vis,scm_FDIF_vis,SCM_FDIR_nir,
c    &                          scm_FDIF_nir
c886      format(1x,'RAD_DRV  scm vis FDIR FDIF ',2(f10.4),
c    &             '     nir FDIR FDIF ',2(f10.4))
      endif
#endif

C**** Save clear sky/tropopause diagnostics here
      AIJ(I,J,IJ_CLR_SRINCG)=AIJ(I,J,IJ_CLR_SRINCG)+OPNSKY*
     *     SRDFLB(1)*CSZ2
      AIJ(I,J,IJ_CLR_SRNFG)=AIJ(I,J,IJ_CLR_SRNFG)+OPNSKY*
     *     SRNFLB(1)*CSZ2
      AIJ(I,J,IJ_CLR_TRDNG)=AIJ(I,J,IJ_CLR_TRDNG)+OPNSKY*TRHR(0,I,J)
      AIJ(I,J,IJ_CLR_SRUPTOA)=AIJ(I,J,IJ_CLR_SRUPTOA)+OPNSKY*
     *     SRUFLB(LM+LM_REQ+1)*CSZ2
      AIJ(I,J,IJ_CLR_TRUPTOA)=AIJ(I,J,IJ_CLR_TRUPTOA)+OPNSKY*
     *     TRUFLB(LM+LM_REQ+1)
      AIJ(I,J,IJ_CLR_SRNTP)=AIJ(I,J,IJ_CLR_SRNTP)+OPNSKY*
     *     SRNFLB(LTROPO(I,J))*CSZ2
      AIJ(I,J,IJ_CLR_TRNTP)=AIJ(I,J,IJ_CLR_TRNTP)+OPNSKY*
     *     TRNFLB(LTROPO(I,J))
      AIJ(I,J,IJ_SRNTP)=AIJ(I,J,IJ_SRNTP)+SRNFLB(LTROPO(I,J))*CSZ2
      AIJ(I,J,IJ_TRNTP)=AIJ(I,J,IJ_TRNTP)+TRNFLB(LTROPO(I,J))
      AIJ(I,J,IJ_SISWD)=AIJ(I,J,IJ_SISWD)+POICE*SRDFLB(1)*CSZ2 
      AIJ(I,J,IJ_SISWU)=AIJ(I,J,IJ_SISWU)+
     *     POICE*(SRDFLB(1)-FSRNFG(3))*CSZ2

      DO IT=1,NTYPE
         call inc_aj(i,j,it,J_CLRTOA,OPNSKY*(SRNFLB(LM+LM_REQ+1)
     *        *CSZ2-TRNFLB(LM+LM_REQ+1))*FTYPE(IT,I,J))
         call inc_aj(i,j,it,J_CLRTRP,OPNSKY*(SRNFLB(LTROPO(I,J))
     *        *CSZ2-TRNFLB(LTROPO(I,J)))*FTYPE(IT,I,J))
         call inc_aj(i,j,it,J_TOTTRP,(SRNFLB(LTROPO(I,J))
     *        *CSZ2-TRNFLB(LTROPO(I,J)))*FTYPE(IT,I,J))
      END DO
      call inc_areg(i,j,jr,J_CLRTOA,OPNSKY*(SRNFLB(LM+LM_REQ+1)
     *     *CSZ2-TRNFLB(LM+LM_REQ+1)))
      call inc_areg(i,j,jr,J_CLRTRP,OPNSKY*(SRNFLB(LTROPO(I,J))
     *     *CSZ2-TRNFLB(LTROPO(I,J))))
      call inc_areg(i,j,jr,J_TOTTRP,(SRNFLB(LTROPO(I,J))
     *     *CSZ2-TRNFLB(LTROPO(I,J))))
C**** Save cloud top diagnostics here
      if (CLDCV.le.0.) go to 590
      AIJ(I,J,IJ_CLDTPPR)=AIJ(I,J,IJ_CLDTPPR)+plb(ltopcl+1)
      AIJ(I,J,IJ_CLDTPT)=AIJ(I,J,IJ_CLDTPT)+(tlb(ltopcl+1) - tf)
C**** Save cloud tau=1 related diagnostics here (opt.depth=1 level)
      tauup=0.
      DO L=LM,1,-1
         taucl=tauwc(l)+tauic(l)
         taudn=tauup+taucl
         if (taudn.gt.1.) then
            aij(i,j,ij_cldcv1)=aij(i,j,ij_cldcv1)+1.
            wtlin=(1.-tauup)/taucl
            aij(i,j,ij_cldt1t)=aij(i,j,ij_cldt1t)+( tlb(l+1)-tf +
     +           (tlb(l)-tlb(l+1))*wtlin )
            aij(i,j,ij_cldt1p)=aij(i,j,ij_cldt1p)+( plb(l+1)+
     +           (plb(l)-plb(l+1))*wtlin )
            go to 590
         end if
      end do
 590  continue

      END DO
C****
C**** END OF MAIN LOOP FOR I INDEX
C****
!$OMP  END DO
!$OMP  END PARALLEL

      END DO
C****
C**** END OF MAIN LOOP FOR J INDEX
C****

      if(kradia.gt.0) then
         call stopTimer('RADIA()')
         return
      end if
C**** Stop if temperatures were out of range
C**** Now only warning messages are printed for T,Q errors
c     IF(ICKERR.GT.0)
c     &     call stop_model('In Radia: Temperature out of range',11)
c     IF(JCKERR.GT.0)  call stop_model('In Radia: RQT out of range',11)
c     IF(KCKERR.GT.0)  call stop_model('In Radia: Q<0',255)
C****
C**** ACCUMULATE THE RADIATION DIAGNOSTICS
C****
      bydpreq(:) = 1d0/(req_fac_d(:)*pmtop)
      DO 780 J=J_0,J_1
         DO 770 I=I_0,IMAXJ(J)
            do l=1,lm
              call inc_ajl(i,j,l,jl_srhr,SRHR(L,I,J)*COSZ2(I,J))
              call inc_ajl(i,j,l,jl_trcr,TRHR(L,I,J))
            enddo
            CSZ2=COSZ2(I,J)
            JR=JREG(I,J)
            DO LR=1,LM_REQ
              call inc_asjl(i,j,lr,3,bydpreq(lr)*SRHRS(LR,I,J)*CSZ2)
              call inc_asjl(i,j,lr,4,bydpreq(lr)*TRHRS(LR,I,J))
            END DO
            DO KR=1,NDIUPT
            IF (I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
              TMP(IDD_PALB)=(1.-SNFS(3,I,J)/S0)
              TMP(IDD_GALB)=(1.-ALB(I,J,1))
              TMP(IDD_ABSA)=(SNFS(3,I,J)-SRHR(0,I,J))*CSZ2
              DO INCH=1,NRAD
                IHM=1+(JTIME+INCH-1)*HR_IN_DAY/NDAY
                IH=IHM
                IF(IH.GT.HR_IN_DAY) IH = IH - HR_IN_DAY
                ADIURN(IDXB(:),KR,IH)=ADIURN(IDXB(:),KR,IH)+TMP(IDXB(:))
#ifndef NO_HDIURN
                IHM = IHM+(JDATE-1)*HR_IN_DAY
                IF(IHM.LE.HR_IN_MONTH) THEN
                  HDIURN(IDXB(:),KR,IHM)=HDIURN(IDXB(:),KR,IHM)+
     &                 TMP(IDXB(:))
                ENDIF
#endif
              ENDDO
            END IF
            END DO

      DO IT=1,NTYPE
         call inc_aj(I,J,IT,J_SRINCP0,(S0*CSZ2)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_SRNFP0 ,(SNFS(3,I,J)*CSZ2)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_SRINCG ,(SRHR(0,I,J)*CSZ2/(ALB(I,J,1)+1.D
     *        -20))*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_BRTEMP ,BTMPW(I,J) *FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_TRINCG ,TRINCG(I,J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_HSURF  ,-(TNFS(3,I,J)-TNFS(1,I,J))
     *        *FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_TRNFP0 ,-TNFS(3,I,J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_TRNFP1 ,-TNFS(2,I,J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_SRNFP1 ,SNFS(2,I,J)*CSZ2*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_HATM   ,-(TNFS(2,I,J)-TNFS(1,I,J))
     *        *FTYPE(IT,I,J))
#ifdef HEALY_LM_DIAGS
         call inc_aj(I,J,IT,j_vtau ,10.*-20*VTAULAT(J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,j_ghg ,10.*ghg_totforc*FTYPE(IT,I,J))
#endif

      END DO
C**** Note: confusing because the types for radiation are a subset
         call inc_aj(I,J,ITOCEAN,J_SRNFG,(FSF(1,I,J)*CSZ2)*FOCEAN(I,J)
     *        *(1.-RSI(I,J)))
         call inc_aj(I,J,ITLAKE ,J_SRNFG,(FSF(1,I,J)*CSZ2)* FLAKE(I,J)
     *        *(1.-RSI(I,J)))
         call inc_aj(I,J,ITEARTH,J_SRNFG,(FSF(4,I,J)*CSZ2)*FEARTH(I,J))
         call inc_aj(I,J,ITLANDI,J_SRNFG,(FSF(3,I,J)*CSZ2)* FLICE(I,J))
         call inc_aj(I,J,ITOICE ,J_SRNFG,(FSF(2,I,J)*CSZ2)*FOCEAN(I,J)
     *        *RSI(I,J))
         call inc_aj(I,J,ITLKICE,J_SRNFG,(FSF(2,I,J)*CSZ2)* FLAKE(I,J)
     *        *RSI(I,J))
C****
         call inc_areg(I,J,JR,J_SRINCP0,(S0*CSZ2))
         call inc_areg(I,J,JR,J_SRNFP0 ,(SNFS(3,I,J)*CSZ2))
         call inc_areg(I,J,JR,J_SRNFP1 ,(SNFS(2,I,J)*CSZ2))
         call inc_areg(I,J,JR,J_SRINCG ,(SRHR(0,I,J)*CSZ2/
     *        (ALB(I,J,1)+1.D-20)))
         call inc_areg(I,J,JR,J_HATM,-(TNFS(2,I,J)-TNFS(1,I,J)))
         call inc_areg(I,J,JR,J_SRNFG,(SRHR(0,I,J)*CSZ2))
         call inc_areg(I,J,JR,J_HSURF,-(TNFS(3,I,J)-TNFS(1,I,J)))
         call inc_areg(I,J,JR,J_BRTEMP, BTMPW(I,J))
         call inc_areg(I,J,JR,J_TRINCG,TRINCG(I,J))
         call inc_areg(I,J,JR,J_TRNFP0,-TNFS(3,I,J))
         call inc_areg(I,J,JR,J_TRNFP1,-TNFS(2,I,J))
         DO K=2,9
           JK=AJ_ALB_INDS(K-1) ! accumulate 8 radiation diags.
           DO IT=1,NTYPE
             call inc_aj(I,J,IT,JK,(S0*CSZ2)*ALB(I,J,K)*FTYPE(IT,I,J))
           END DO
           call inc_areg(I,J,JR,JK,(S0*CSZ2)*ALB(I,J,K))
         END DO
         AIJ(I,J,IJ_SRINCG) =AIJ(I,J,IJ_SRINCG) +(SRHR(0,I,J)*CSZ2/
     /        (ALB(I,J,1)+1.D-20))
         AIJ(I,J,IJ_SRNFG) =AIJ(I,J,IJ_SRNFG) +(SRHR(0,I,J)*CSZ2)
         AIJ(I,J,IJ_BTMPW) =AIJ(I,J,IJ_BTMPW) +BTMPW(I,J)
         AIJ(I,J,IJ_SRREF) =AIJ(I,J,IJ_SRREF) +S0*CSZ2*ALB(I,J,2)
         AIJ(I,J,IJ_SRVIS) =AIJ(I,J,IJ_SRVIS) +S0*CSZ2*ALB(I,J,4)
         AIJ(I,J,IJ_TRNFP0)=AIJ(I,J,IJ_TRNFP0)-TNFS(3,I,J)
         AIJ(I,J,IJ_SRNFP0)=AIJ(I,J,IJ_SRNFP0)+(SNFS(3,I,J)*CSZ2)
         AIJ(I,J,IJ_RNFP1) =AIJ(I,J,IJ_RNFP1) +(SNFS(2,I,J)*CSZ2
     *                                         -TNFS(2,I,J))
         AIJ(I,J,ij_srvdir)=AIJ(I,J,ij_srvdir)
     &        + FSRDIR(I,J)*SRVISSURF(I,J)
         AIJ(I,J,IJ_SRVISSURF)=AIJ(I,J,IJ_SRVISSURF)+SRVISSURF(I,J)
#ifdef CACHED_SUBDD
         OLR_ij(I,J) = -TNFS(3,I,J)
#endif
C**** CRF diags if required
         if (moddrf .ne. 0) go to 770
         if (cloud_rad_forc > 0) then
c    CRF diagnostics 
           AIJ(I,J,IJ_SWCRF)=AIJ(I,J,IJ_SWCRF)+
     +          (SNFS(3,I,J)-SNFSCRF(I,J))*CSZ2
           AIJ(I,J,IJ_LWCRF)=AIJ(I,J,IJ_LWCRF)-
     -          (TNFS(3,I,J)-TNFSCRF(I,J))
c    CRF diagnostics without aerosols and Ox
           AIJ(I,J,IJ_SWCRF2)=AIJ(I,J,IJ_SWCRF2)+
     +          (SNFS(3,I,J)-SNFSCRF2(I,J))*CSZ2
           AIJ(I,J,IJ_LWCRF2)=AIJ(I,J,IJ_LWCRF2)-
     -          (TNFS(3,I,J)-TNFSCRF2(I,J))
         end if

C**** AERRF diags if required
         if (aer_rad_forc > 0) then
           do N=1,8
             AIJ(I,J,IJ_SWAERRF+N-1)=AIJ(I,J,IJ_SWAERRF+N-1)+
     *            (SNFS(3,I,J)-SNFSAERRF(N,I,J))*CSZ2
             AIJ(I,J,IJ_LWAERRF+N-1)=AIJ(I,J,IJ_LWAERRF+N-1)-
     *            (TNFS(3,I,J)-TNFSAERRF(N,I,J))
             AIJ(I,J,IJ_SWAERSRF+N-1)=AIJ(I,J,IJ_SWAERSRF+N-1)+
     *            (SNFS(1,I,J)-SNFSAERRF(N+8,I,J))*CSZ2
             AIJ(I,J,IJ_LWAERSRF+N-1)=AIJ(I,J,IJ_LWAERSRF+N-1)-
     *            (TNFS(1,I,J)-TNFSAERRF(N+8,I,J))
           end do
           AIJ(I,J,IJ_SWAERRFNT)=AIJ(I,J,IJ_SWAERRFNT)+
     *          (SNFS(3,I,J)-SNFSAERRF(17,I,J))*CSZ2
           AIJ(I,J,IJ_LWAERRFNT)=AIJ(I,J,IJ_LWAERRFNT)-
     *          (TNFS(3,I,J)-TNFSAERRF(17,I,J))
           AIJ(I,J,IJ_SWAERSRFNT)=AIJ(I,J,IJ_SWAERSRFNT)+
     *          (SNFS(1,I,J)-SNFSAERRF(18,I,J))*CSZ2
           AIJ(I,J,IJ_LWAERSRFNT)=AIJ(I,J,IJ_LWAERSRFNT)-
     *          (TNFS(1,I,J)-TNFSAERRF(18,I,J))
         end if

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP)
C**** Generic diagnostics for radiative forcing calculations
C**** Depending on whether tracers radiative interaction is turned on,
C**** diagnostic sign changes (for aerosols)
         rsign_aer=1. ; rsign_chem=-1.
         if (rad_interact_aer > 0) rsign_aer=-1.
C**** define SNFS/TNFS level (TOA/TROPO) for calculating forcing
         LFRC=3                 ! TOA
         if (rad_forc_lev.gt.0) LFRC=4 ! TROPOPAUSE
#ifdef  TRACERS_AMP
         IF (AMP_DIAG_FC == 2) THEN
            NTRACE = nmodes
         ELSE
            NTRACE = 1
            NTRIX(1) = 1
         ENDIF
#endif /* TRACERS_AMP */
         if (ntrace > 0) then
#ifdef BC_ALB
      if (ijts_alb(1).gt.0)
     * TAIJS(I,J,ijts_alb(1))=TAIJS(I,J,ijts_alb(1))
     *   + 100.d0*(ALBNBC(I,J)-ALB(I,J,1))
      if (ijts_alb(2).gt.0)
     & taijs(i,j,ijts_alb(2))
     &     =taijs(i,j,ijts_alb(2))
     &         +(SNFS(3,I,J)-NFSNBC(I,J))*CSZ2
#endif /* BC_ALB */
#ifdef TRACERS_AEROSOLS_Koch
c          snfst0(:,:,i,j)=0.D0
c          tnfst0(:,:,i,j)=0.D0
#endif /* TRACERS_AEROSOLS_Koch */
c     ..........
c     accumulation of forcings for tracers for which ntrace fields are
c     defined
c     ..........
           set_clayilli=.FALSE.
           set_claykaol=.FALSE.
           set_claysmec=.FALSE.
           set_claycalc=.FALSE.
           set_clayquar=.FALSE.
           do n=1,ntrace
             IF (ntrix(n) > 0) THEN
               SELECT CASE (trname(ntrix(n)))
               CASE ('Clay')
                 n1=n-nrad_clay+1
c shortwave forcing (TOA or TROPO) of Clay sub size classes
                 if (ijts_fcsub(1,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(1,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(1,ntrix(n),n1))
     &                +rsign_aer*(snfst(2,n,i,j)-snfs(lfrc,i,j))*csz2
c longwave forcing  (TOA or TROPO) of Clay size sub classes
                 if (ijts_fcsub(2,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(2,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(2,ntrix(n),n1))
     &                -rsign_aer*(tnfst(2,n,i,j)-tnfs(lfrc,i,j))
c shortwave forcing (TOA or TROPO) clear sky of Clay sub size classes
                 if (ijts_fcsub(5,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(5,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(5,ntrix(n),n1))
     &                +rsign_aer*(snfst(2,n,i,j)-snfs(lfrc,i,j))*csz2
     &                *(1.D0-cfrac(i,j))
c longwave forcing  (TOA or TROPO) clear sky of Clay sub size classes
                 if (ijts_fcsub(6,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(6,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(6,ntrix(n),n1))
     &                -rsign_aer*(tnfst(2,n,i,j)-tnfs(lfrc,i,j))
     &                *(1.D0-cfrac(i,j))
c shortwave forcing at surface (if required) of Clay sub size classes
                 if (ijts_fcsub(3,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(3,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(3,ntrix(n),n1))
     &                +rsign_aer*(snfst(1,n,i,j)-snfs(1,i,j))*csz2
c longwave forcing at surface (if required) of Clay sub size classes
                 if (ijts_fcsub(4,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(4,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(4,ntrix(n),n1))
     &                -rsign_aer*(tnfst(1,n,i,j)-tnfs(1,i,j))
c shortwave forcing at surface clear sky (if required) of Clay sub size classes
                 if (ijts_fcsub(7,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(7,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(7,ntrix(n),n1))
     &                +rsign_aer*(snfst(1,n,i,j)-snfs(1,i,j))*csz2
     &                *(1.D0-cfrac(i,j))
c longwave forcing at surface clear sky (if required) of Clay sub size classes
                 if (ijts_fcsub(8,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(8,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(8,ntrix(n),n1))
     &                -rsign_aer*(tnfst(1,n,i,j)-tnfs(1,i,j))
     &                *(1.D0-cfrac(i,j))
               CASE DEFAULT
                 SELECT CASE (trname(ntrix(n)))
                 CASE ('seasalt2')
                   CYCLE
                 CASE ('ClayIlli')
                   IF (set_clayilli) CYCLE
                   set_clayilli=.TRUE.
                 CASE ('ClayKaol')
                   IF (set_claykaol) CYCLE
                   set_claykaol=.TRUE.
                 CASE ('ClaySmec')
                   IF (set_claysmec) CYCLE
                   set_claysmec=.TRUE.
                 CASE ('ClayCalc')
                   IF (set_claycalc) CYCLE
                   set_claycalc=.TRUE.
                 CASE ('ClayQuar')
                   IF (set_clayquar) CYCLE
                   set_clayquar=.TRUE.
                 END SELECT
c shortwave forcing (TOA or TROPO)
                 if (ijts_fc(1,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(1,ntrix(n)))
     &                =taijs(i,j,ijts_fc(1,ntrix(n)))
     &                +rsign_aer*(SNFST(2,N,I,J)-SNFS(LFRC,I,J))*CSZ2
c longwave forcing  (TOA or TROPO)
                 if (ijts_fc(2,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(2,ntrix(n)))
     &                =taijs(i,j,ijts_fc(2,ntrix(n)))
     &                -rsign_aer*(TNFST(2,N,I,J)-TNFS(LFRC,I,J))
c shortwave forcing (TOA or TROPO) clear sky
                 if (ijts_fc(5,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(5,ntrix(n)))
     &                =taijs(i,j,ijts_fc(5,ntrix(n)))
     &                +rsign_aer*(SNFST(2,N,I,J)-SNFS(LFRC,I,J))*CSZ2
     &                *(1.d0-CFRAC(I,J))
c longwave forcing  (TOA or TROPO) clear sky
                 if (ijts_fc(6,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(6,ntrix(n)))
     &                =taijs(i,j,ijts_fc(6,ntrix(n)))
     &                -rsign_aer*(TNFST(2,N,I,J)-TNFS(LFRC,I,J))
     &                *(1.d0-CFRAC(I,J))
c shortwave forcing at surface (if required)
                 if (ijts_fc(3,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(3,ntrix(n)))
     &                =taijs(i,j,ijts_fc(3,ntrix(n)))
     &                +rsign_aer*(SNFST(1,N,I,J)-SNFS(1,I,J))*CSZ2
c longwave forcing at surface (if required)
                 if (ijts_fc(4,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(4,ntrix(n)))
     &                =taijs(i,j,ijts_fc(4,ntrix(n)))
     &                -rsign_aer*(TNFST(1,N,I,J)-TNFS(1,I,J))
c shortwave forcing at surface clear sky (if required)
                 if (ijts_fc(7,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(7,ntrix(n)))
     &                =taijs(i,j,ijts_fc(7,ntrix(n)))
     &                +rsign_aer*(SNFST(1,N,I,J)-SNFS(1,I,J))*CSZ2
     &                *(1.d0-CFRAC(I,J))
c longwave forcing at surface clear sky (if required)
                 if (ijts_fc(8,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(8,ntrix(n)))
     &                =taijs(i,j,ijts_fc(8,ntrix(n)))
     &                -rsign_aer*(TNFST(1,N,I,J)-TNFS(1,I,J))
     &                *(1.d0-CFRAC(I,J))
               END SELECT
#ifdef  TRACERS_AMP
         IF (AMP_DIAG_FC == 2) THEN
         ELSE
         NTRIX(1)=  n_N_AKK_1
         ENDIF
#endif /* TRACERS_AMP */
#ifdef TRACERS_AEROSOLS_Koch
c              SNFST0(1,ntrix(n),I,J)=SNFST0(1,ntrix(n),I,J)
c    &              +rsign_aer*(SNFST(2,n,I,J)-SNFS(LFRC,I,J))*CSZ2
c              SNFST0(2,ntrix(n),I,J)=SNFST0(2,ntrix(n),I,J)
c    &              +rsign_aer*(SNFST(1,n,I,J)-SNFS(1,I,J))*CSZ2
c              TNFST0(1,ntrix(n),I,J)=TNFST0(1,ntrix(n),I,J)
c    &              -rsign_aer*(TNFST(2,n,I,J)-TNFS(LFRC,I,J))
c              TNFST0(2,ntrix(n),I,J)=TNFST0(2,ntrix(n),I,J)
c    &              -rsign_aer*(TNFST(1,n,I,J)-TNFS(1,I,J))
#endif /* TRACERS_AEROSOLS_Koch */
             END IF   ! ntrix(n)>0
           end do     ! n=1,ntrace
         end if       ! ntrace>0

c ..........
c accumulation of forcings for special case ozone (ntrace fields
c not defined) Warning: indicies used differently, since we don't
c need CS or Surface, but are doing both TOA and Ltropo:
c ..........
         if(n_Ox > 0) then ! ------ main Ox tracer -------
c shortwave forcing at tropopause
           if (ijts_fc(1,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(1,n_Ox))=taijs(i,j,ijts_fc(1,n_Ox))
     &     +rsign_chem*(SNFST_o3ref(1,I,J)-SNFS(4,I,J))*CSZ2
c longwave forcing at tropopause
           if (ijts_fc(2,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(2,n_Ox))=taijs(i,j,ijts_fc(2,n_Ox))
     &     -rsign_chem*(TNFST_o3ref(1,I,J)-TNFS(4,I,J))
c shortwave forcing at TOA
           if (ijts_fc(3,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(3,n_Ox))=taijs(i,j,ijts_fc(3,n_Ox))
     &     +rsign_chem*(SNFST_o3ref(2,I,J)-SNFS(3,I,J))*CSZ2
c longwave forcing at TOA
           if (ijts_fc(4,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(4,n_Ox))=taijs(i,j,ijts_fc(4,n_Ox))
     &     -rsign_chem*(TNFST_o3ref(2,I,J)-TNFS(3,I,J))
         endif
#ifdef AUXILIARY_OX_RADF
c shortwave forcing at tropopause
         if (ijts_auxfc(1)>0)
     &   taijs(i,j,ijts_auxfc(1))=taijs(i,j,ijts_auxfc(1))
#ifdef AUX_OX_RADF_TROP
     &   +rsign_chem*(SNFST_o3ref(5,I,J)-SNFST_o3ref(3,I,J))*CSZ2
#else
     &   +rsign_chem*(SNFST_o3ref(1,I,J)-SNFST_o3ref(3,I,J))*CSZ2
#endif
c longwave forcing at tropopause
         if (ijts_auxfc(2)>0)
     &   taijs(i,j,ijts_auxfc(2))=taijs(i,j,ijts_auxfc(2))
#ifdef AUX_OX_RADF_TROP
     &   -rsign_chem*(TNFST_o3ref(5,I,J)-TNFST_o3ref(3,I,J))
#else
     &   -rsign_chem*(TNFST_o3ref(1,I,J)-TNFST_o3ref(3,I,J))
#endif
c shortwave forcing at TOA
         if (ijts_auxfc(3)>0)
     &   taijs(i,j,ijts_auxfc(3))=taijs(i,j,ijts_auxfc(3))
     &   +rsign_chem*(SNFST_o3ref(2,I,J)-SNFST_o3ref(4,I,J))*CSZ2
c longwave forcing at TOA
         if (ijts_auxfc(4)>0)
     &   taijs(i,j,ijts_auxfc(4))=taijs(i,j,ijts_auxfc(4))
     &   -rsign_chem*(TNFST_o3ref(2,I,J)-TNFST_o3ref(4,I,J))
#endif /* AUXILIARY_OX_RADF */
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
                      ! ------ diag stratOx tracer -------
! note for now for this diag, there is a failsafe that stops model
! if clim_interact_chem .le. 0 when the below would be wrong:
c shortwave forcing at tropopause
         if (ijts_fc(1,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(1,n_stratOx))=taijs(i,j,ijts_fc(1,n_stratOx))
     &   +rsign_chem*(SNFST_o3ref(1,I,J)-SNFST_stratOx(1,I,J))*CSZ2
c longwave forcing at tropopause
         if (ijts_fc(2,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(2,n_stratOx))=taijs(i,j,ijts_fc(2,n_stratOx))
     &   -rsign_chem*(TNFST_o3ref(1,I,J)-TNFST_stratOx(1,I,J))
c shortwave forcing at TOA
         if (ijts_fc(3,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(3,n_stratOx))=taijs(i,j,ijts_fc(3,n_stratOx))
     &   +rsign_chem*(SNFST_o3ref(2,I,J)-SNFST_stratOx(2,I,J))*CSZ2
c longwave forcing at TOA
         if (ijts_fc(4,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(4,n_stratOx))=taijs(i,j,ijts_fc(4,n_stratOx))
     &   -rsign_chem*(TNFST_o3ref(2,I,J)-TNFST_stratOx(2,I,J))
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS*/
#endif /* any of various tracer groups defined */

#ifdef ACCMIP_LIKE_DIAGS
         do nf=1,4 ! CH4, N2O, CFC11, and CFC12:
c shortwave GHG forcing at TOA
           if(ij_fcghg(1,nf).gt.0)aij(i,j,ij_fcghg(1,nf))=
     &     aij(i,j,ij_fcghg(1,nf))+(SNFS(3,I,J)-SNFS_ghg(nf,I,J))
     &     *CSZ2
c longwave GHG forcing at TOA
           if(ij_fcghg(2,nf).gt.0)aij(i,j,ij_fcghg(2,nf))=
     &     aij(i,j,ij_fcghg(2,nf))+(TNFS_ghg(nf,I,J)-TNFS(3,I,J))
         enddo
#endif /* ACCMIP_LIKE_DIAGS */

  770    CONTINUE
  780    CONTINUE

      DO J=J_0,J_1
      DO I=I_0,I_1
      DO L=1,LM
        AIJL(i,j,l,IJL_RC)=AIJL(i,j,l,IJL_RC)+
     &       (SRHR(L,I,J)*COSZ2(I,J)+TRHR(L,I,J))
      END DO
      END DO
      END DO

#ifdef CACHED_SUBDD
      do k=1,subdd_ngroups
        subdd => subdd_groups(k)
        subdd%nacc(subdd%subdd_period,sched_rad) =
     &  subdd%nacc(subdd%subdd_period,sched_rad) + 1
      enddo
C****
C**** Collect some high-frequency outputs
C****
      call find_groups('aijh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
      case ('olrrad')
        call inc_subdd(subdd,k,OLR_ij)
        OLR_ij(:,:) = 0.
c
      case ('swu')
        call inc_subdd(subdd,k,SWU_ij)
        SWU_ij(:,:) = 0.
c
      case ('swd')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j)=srdn(i,j)*cosz1(i,j)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
c
      case ('swdf')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j)=fsrdif(i,j)+difnir(i,j)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
c
      case ('tcld')
        call inc_subdd(subdd,k,cfrac)
      end select
      enddo
      enddo
#endif  /* CACHED_SUBDD */

C****
C**** Update radiative equilibrium temperatures
C****
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          DO LR=1,LM_REQ
            RQT(LR,I,J)=RQT(LR,I,J)+(SRHRS(LR,I,J)*COSZ2(I,J)
     &           +TRHRS(LR,I,J))*NRAD*DTsrc*bysha*byaml00(lr+lm)
          END DO
        END DO
      END DO

#ifdef CACHED_SUBDD
C**** ACCUMULATE SUBDAILY ON MODEL LEVELS.
C**** for interpolating to pressure levels in RADIA, more
C**** consideration is needed, speak with M.Kelly
      call find_groups('aijlh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
         case ('swhr')
            do j=j_0,j_1; do i=i_0,imaxj(j); do l=1,lmaxsubdd
              sddarr3d(i,j,l) = SRHR(L,I,J)*bysha*byam(L,I,J)*COSZ2(I,J)
            enddo;                enddo;    enddo
            call inc_subdd(subdd,k,sddarr3d)
         case ('lwhr')
            do j=j_0,j_1; do i=i_0,imaxj(j); do l=1,lmaxsubdd
              sddarr3d(i,j,l) = TRHR(L,I,J)*bysha*byam(L,I,J)
            enddo;           enddo;         enddo
            call inc_subdd(subdd,k,sddarr3d)
         end select
      enddo
      enddo
#endif
C****
C**** Update other temperatures every physics time step
C****
  900 DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          DO L=1,LM
            T(I,J,L)=T(I,J,L)+(SRHR(L,I,J)*COSZ1(I,J)+TRHR(L,I,J))*
     *           DTsrc*bysha*byam(l,i,j)/PK(L,I,J)
#ifdef SCM
            if (I.eq.I_TARG.and.J.eq.J_TARG) then
              dTradlw(L)=TRHR(L,I,J)*DTsrc*bysha*byam(l,i,j)/PK(L,I,J)
              dTradsw(L)=SRHR(L,I,J)*COSZ1(I,J)
     &                  *DTsrc*bysha*byam(L,I,J)/PK(L,I,J)
            endif
#endif
          END DO
          AIJ(I,J,IJ_SRINCP0)=AIJ(I,J,IJ_SRINCP0)+(S0*COSZ1(I,J))
        END DO
      END DO

C**** daily diagnostics
      IH=1+JHOUR
      IHM = IH+(JDATE-1)*24
      DO KR=1,NDIUPT
        I = IJDD(1,KR)
        J = IJDD(2,KR)
        IF ((J >= J_0) .AND. (J <= J_1) .AND.
     &      (I >= I_0) .AND. (I <= I_1)) THEN
          ADIURN(IDD_ISW,KR,IH)=ADIURN(IDD_ISW,KR,IH)+S0*COSZ1(I,J)
#ifndef NO_HDIURN
          HDIURN(IDD_ISW,KR,IHM)=HDIURN(IDD_ISW,KR,IHM)+S0*COSZ1(I,J)
#endif
        ENDIF
      ENDDO

      call stopTimer('RADIA()')
      RETURN
      END SUBROUTINE RADIA

      SUBROUTINE RESET_SURF_FLUXES(I,J,ITYPE_OLD,ITYPE_NEW,FTYPE_ORIG,
     *     FTYPE_NOW)
!@sum set incident solar and upward thermal fluxes appropriately
!@+   as fractions change to conserve energy, prevent restart problems
!@auth Gavin Schmidt
      use rad_com, only : fsf,trsurf
      implicit none
!@var itype_old, itype_new indices for the old type turning to new type
      integer, intent(in) :: i,j,itype_old,itype_new
!@var ftype_orig, ftype_now original and current fracs of the 'new' type
      real*8, intent(in) :: ftype_orig, ftype_now
      real*8 :: delf ! change in fraction from old to new

      if (( (ITYPE_OLD==1 .and. ITYPE_NEW==2)
     *  .or.(ITYPE_OLD==2 .and. ITYPE_NEW==1))
     *  .and. (FTYPE_NOW .le. 0. .or. FTYPE_NOW .gt. 1.)) then
        write (6,*) ' RESET_SURF_FLUXES:
     *    I, J, ITYPE_OLD, ITYPE_NEW, FTYPE_ORIG, FTYPE_NOW = ',
     *    I, J, ITYPE_OLD, ITYPE_NEW, FTYPE_ORIG, FTYPE_NOW
        call stop_model('RESET_SURF_FLUXES: INCORRECT RESET',255)
      end if

      delf = FTYPE_NOW-FTYPE_ORIG
C**** Constrain fsf_1*ftype_1+fsf_2*ftype_2 to be constant
      FSF(ITYPE_NEW,I,J)=(FSF(ITYPE_NEW,I,J)*FTYPE_ORIG+
     *     FSF(ITYPE_OLD,I,J)*DELF)/FTYPE_NOW

C**** Same for upward thermal
      TRSURF(ITYPE_NEW,I,J)=(TRSURF(ITYPE_NEW,I,J)*FTYPE_ORIG+
     *     TRSURF(ITYPE_OLD,I,J)*DELF)/FTYPE_NOW

      RETURN
      END SUBROUTINE RESET_SURF_FLUXES

      SUBROUTINE GHGHST(iu)
!@sum  reads history for nghg well-mixed greenhouse gases
!@auth R. Ruedy
!@ver  1.0

      use domain_decomp_atm, only : write_parallel
      USE RADPAR, only : nghg,ghgyr1,ghgyr2,ghgam
      USE RAD_COM, only : ghg_yr
      IMPLICIT NONE
      INTEGER :: iu,n,k,nhead=4,iyr
      CHARACTER*80 title
      character(len=300) :: out_line

      write(out_line,*)  ! print header lines and first data line
      call write_parallel(trim(out_line),unit=6)
      do n=1,nhead+1
        read(iu,'(a)') title
        write(out_line,'(1x,a80)') title
        call write_parallel(trim(out_line),unit=6)
      end do
      if(title(1:2).eq.'--') then                 ! older format
        read(iu,'(a)') title
        write(out_line,'(1x,a80)') title
        call write_parallel(trim(out_line),unit=6)
        nhead=5
      end if

!**** find range of table: ghgyr1 - ghgyr2
      read(title,*) ghgyr1
      do ; read(iu,'(a)',end=20) title ; end do
   20 read(title,*) ghgyr2
      rewind iu  !   position to data lines
      do n=1,nhead ; read(iu,'(a)') ; end do

      allocate (ghgam(nghg,ghgyr2-ghgyr1+1))
      do n=1,ghgyr2-ghgyr1+1
        read(iu,*) iyr,(ghgam(k,n),k=1,nghg)
        do k=1,nghg ! replace -999. by reasonable numbers
          if(ghgam(k,n).lt.0.) ghgam(k,n)=ghgam(k,n-1)
        end do
        if(ghg_yr>0 .and. abs(ghg_yr-iyr).le.1) then
          write(out_line,'(i5,6f10.4)') iyr,(ghgam(k,n),k=1,nghg)
          call write_parallel(trim(out_line),unit=6)
        endif
      end do
      write(out_line,*) 'read GHG table for years',ghgyr1,' - ',ghgyr2
      call write_parallel(trim(out_line),unit=6)
      return
      end SUBROUTINE GHGHST

#if defined(CUBED_SPHERE)
      subroutine read_qma (iu,plb)
!@sum  reads H2O production rates induced by CH4 (Tim Hall)
!@auth R. Ruedy
!@ver  1.0
      use domain_decomp_atm, only : write_parallel
      use rad_com, only : dH2O,jma=>jm_dh2o,lat_dh2o
      use model_com, only : lm
      use constant, only : radian
      implicit none
      integer, parameter:: lma=24
      integer m,iu,j,l,ll,ldn(lm),lup(lm)
      real*8 :: plb(lm+1)
      real*4 pb(0:lma+1),h2o(jma,0:lma),z(lma),dz(0:lma)
      character*100 title
      real*4 pdn,pup,dh,fracl
      character(len=300) :: out_line

C**** read headers/latitudes
      read(iu,'(a)') title
      write(out_line,'(''0'',a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
      write(out_line,'(1x,a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
c      write(6,'(1x,a100)') title
      read(title(10:100),*) (lat_dh2o(j),j=1,jma)
      lat_dh2o(:) = lat_dh2o(:)*radian

C**** read heights z(km) and data (kg/km^3/year)
      do m=1,12
        read(iu,'(a)') title
        write(out_line,'(1x,a100)') title
        call write_parallel(trim(out_line),unit=6)
        do l=lma,1,-1
          read(iu,'(a)') title
c          write(6,'(1x,a100)') title
          read(title,*) z(l),(H2O(j,l),j=1,jma)
        end do
        do j=1,jma
          h2o(j,0) = 0.
        end do

C**** Find edge heights and pressures
        dz(0) = 0.
        dz(1) = z(2)-z(1)
        do l=2,lma-1
           dz(l)=.5*(z(l+1)-z(l-1))
        end do
        dz(lma) = z(lma)-z(lma-1)

        pb(0) = plb(1)
        do l=1,lma
           Pb(l)=1000.*10.**(-(z(l)-.5*dz(l))/16.)
        end do
C**** extend both systems vertically to p=0
        pb(lma+1)=0.
        plb(lm+1)=0.

C**** Interpolate vertical resolution to model layers
        ldn(:) = 0
        do l=1,lm
          do while (pb(ldn(l)+1).ge.plb(l) .and. ldn(l).lt.lma)
            ldn(l)=ldn(l)+1
          end do
          lup(l)=ldn(l)
          do while (pb(lup(l)+1).gt.plb(l+1) .and. lup(l).lt.lma)
            lup(l)=lup(l)+1
          end do
        end do

C**** Interpolate (extrapolate) vertically
        do j=1,jma
          do l=1,lm
            dh = 0.
            pdn = plb(l)
            if (lup(l).gt.0) then
              do ll=ldn(l),lup(l)
                pup = max(REAL(pb(ll+1),KIND=8),plb(l+1))
                fracl= (pdn-pup)/(pb(ll)-pb(ll+1))
                dh = dh+h2o(j,ll)*fracl*dz(ll)
                pdn = pup
              end do
            end if
            dh2o(j,l,m) = 1.d-6*dh/1.74d0/365. !->(kg/m^2/ppm_CH4/day)
          end do
        end do
      end do
      return
      end subroutine read_qma

      subroutine lat_interp_qma (rlat,lev,mon,dh2o_interp)
!@sum  interpolate CH4->H2O production rates in latitude
!@auth R. Ruedy
!@ver  1.0
      use rad_com, only : jma=>jm_dh2o,xlat=>lat_dh2o,dh2o
      implicit none
      real*8 :: rlat ! input latitude (radians)
      integer :: lev,mon ! input level, month
      real*8 :: dh2o_interp  ! output
      real*8 w1,w2
      integer :: j1,j2

C**** Interpolate (extrapolate) horizontally
      j2 = 2+(jma-1)*(rlat-xlat(1))/(xlat(jma)-xlat(1)) ! first guess
      j2 = min(max(2,j2),jma)
      j1 = j2-1
      if(rlat.gt.xlat(j2)) then ! j guess was too low
         do while (j2.lt.jma .and. rlat.gt.xlat(j2))
          j2 = j2+1
         end do
         j1 = j2-1
      elseif(rlat.lt.xlat(j1)) then ! j guess was too high
         do while (j1.gt.1 .and. rlat.lt.xlat(j1))
          j1 = j1-1
         end do
         j2 = j1+1
      endif
C**** coeff. for latitudinal linear inter/extrapolation
      w1 = (xlat(j2)-rlat)/(xlat(j2)-xlat(j1))
C**** for extrapolations, only use half the slope
      if(w1.gt.1.) w1=.5+.5*w1
      if(w1.lt.0.) w1=.5*w1
      w2 = 1.-w1
      dh2o_interp = w1*dh2o(j1,lev,mon)+w2*dh2o(j2,lev,mon)
      return
      end subroutine lat_interp_qma

#endif /* CUBED_SPHERE */

      subroutine getqma (iu,dglat,plb,dh2o,lm,jm)
!@sum  reads H2O production rates induced by CH4 (Tim Hall)
!@auth R. Ruedy
!@ver  1.0
      use domain_decomp_atm, only : grid,get,write_parallel
      implicit none
      integer, parameter:: jma=18,lma=24
      integer m,iu,jm,lm,j,j1,j2,l,ll,ldn(lm),lup(lm)
      real*8 PLB(lm+1),dH2O(grid%j_strt_halo:grid%j_stop_halo,lm,12)
     &     ,dglat(jm)
      real*4 pb(0:lma+1),h2o(jma,0:lma),xlat(jma),z(lma),dz(0:lma)
      character*100 title
      real*4 pdn,pup,w1,w2,dh,fracl
      integer :: j_0,j_1
      character(len=300) :: out_line
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

C**** read headers/latitudes
      read(iu,'(a)') title
      write(out_line,'(''0'',a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
      write(out_line,'(1x,a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
c      write(6,'(1x,a100)') title
      read(title(10:100),*) (xlat(j),j=1,jma)

C**** read heights z(km) and data (kg/km^3/year)
      do m=1,12
        read(iu,'(a)') title
        write(out_line,'(1x,a100)') title
        call write_parallel(trim(out_line),unit=6)
        do l=lma,1,-1
          read(iu,'(a)') title
c          write(6,'(1x,a100)') title
          read(title,*) z(l),(H2O(j,l),j=1,jma)
        end do
        do j=1,jma
          h2o(j,0) = 0.
        end do

C**** Find edge heights and pressures
        dz(0) = 0.
        dz(1) = z(2)-z(1)
        do l=2,lma-1
           dz(l)=.5*(z(l+1)-z(l-1))
        end do
        dz(lma) = z(lma)-z(lma-1)

        pb(0) = plb(1)
        do l=1,lma
           Pb(l)=1000.*10.**(-(z(l)-.5*dz(l))/16.)
        end do
C**** extend both systems vertically to p=0
        pb(lma+1)=0.
        plb(lm+1)=0.

C**** Interpolate vertical resolution to model layers
        ldn(:) = 0
        do l=1,lm
          do while (pb(ldn(l)+1).ge.plb(l) .and. ldn(l).lt.lma)
            ldn(l)=ldn(l)+1
          end do
          lup(l)=ldn(l)
          do while (pb(lup(l)+1).gt.plb(l+1) .and. lup(l).lt.lma)
            lup(l)=lup(l)+1
          end do
        end do

C**** Interpolate (extrapolate) horizontally and vertically
        j2=2
        do j=j_0,j_1
C**** coeff. for latitudinal linear inter/extrapolation
          do while (j2.lt.jma .and. dglat(j).gt.xlat(j2))
            j2 = j2+1
          end do
          j1 = j2-1
          w1 = (xlat(j2)-dglat(j))/(xlat(j2)-xlat(j1))
C**** for extrapolations, only use half the slope
          if(w1.gt.1.) w1=.5+.5*w1
          if(w1.lt.0.) w1=.5*w1
          w2 = 1.-w1
          do l=1,lm
            dh = 0.
            pdn = plb(l)
            if (lup(l).gt.0) then
              do ll=ldn(l),lup(l)
                pup = max(REAL(pb(ll+1),KIND=8),plb(l+1))
                fracl= (pdn-pup)/(pb(ll)-pb(ll+1))
                dh = dh+(w1*h2o(j1,ll)+w2*h2o(j2,ll))*fracl*dz(ll)
                pdn = pup
              end do
            end if
            dh2o(j,l,m) = 1.d-6*dh/1.74d0/365. !->(kg/m^2/ppm_CH4/day)
          end do
        end do
      end do
      return
      end subroutine getqma

      Subroutine ORBPAR (YEAR, ECCEN,OBLIQ,OMEGVP)
C****
!@sum ORBPAR calculates the three orbital parameters as a function of
!@+   YEAR.  The source of these calculations is: Andre L. Berger,
!@+   1978, "Long-Term Variations of Daily Insolation and Quaternary
!@+   Climatic Changes", JAS, v.35, p.2362.  Also useful is: Andre L.
!@+   Berger, May 1978, "A Simple Algorithm to Compute Long Term
!@+   Variations of Daily Insolation", published by Institut
!@+   D'Astronomie de Geophysique, Universite Catholique de Louvain,
!@+   Louvain-la Neuve, No. 18.
!@+
!@+   Tables and equations refer to the first reference (JAS).  The
!@+   corresponding table or equation in the second reference is
!@+   enclosed in parentheses.  The coefficients used in this
!@+   subroutine are slightly more precise than those used in either
!@+   of the references.  The generated orbital parameters are precise
!@+   within plus or minus 1000000 years from present.
C****
!@auth Gary L. Russell (with extra terms from D. Thresher)
C****
      USE CONSTANT, only : twopi,PI180=>radian
      IMPLICIT NONE
C**** Input:
!@var YEAR   = years C.E. are positive, B.C.E are -ve (i.e 4BCE = -3)
      REAL*8, INTENT(IN) :: YEAR
C**** Output:
!@var ECCEN  = eccentricity of orbital ellipse
!@var OBLIQ  = latitude of Tropic of Cancer in degrees
!@var OMEGVP = longitude of perihelion =
!@+          = spatial angle from vernal equinox to perihelion
!@+            in degrees with sun as angle vertex
      REAL*8, INTENT(OUT) :: ECCEN,OBLIQ,OMEGVP
C**** Table 1 (2).  Obliquity relative to mean ecliptic of date: OBLIQD
      REAL*8, PARAMETER, DIMENSION(3,47) :: TABL1 = RESHAPE( (/
     1            -2462.2214466d0, 31.609974d0, 251.9025d0,
     2             -857.3232075d0, 32.620504d0, 280.8325d0,
     3             -629.3231835d0, 24.172203d0, 128.3057d0,
     4             -414.2804924d0, 31.983787d0, 292.7252d0,
     5             -311.7632587d0, 44.828336d0,  15.3747d0,
     6              308.9408604d0, 30.973257d0, 263.7951d0,
     7             -162.5533601d0, 43.668246d0, 308.4258d0,
     8             -116.1077911d0, 32.246691d0, 240.0099d0,
     9              101.1189923d0, 30.599444d0, 222.9725d0,
     O              -67.6856209d0, 42.681324d0, 268.7809d0,
     1               24.9079067d0, 43.836462d0, 316.7998d0,
     2               22.5811241d0, 47.439436d0, 319.6024d0,
     3              -21.1648355d0, 63.219948d0, 143.8050d0,
     4              -15.6549876d0, 64.230478d0, 172.7351d0,
     5               15.3936813d0,  1.010530d0,  28.9300d0,
     6               14.6660938d0,  7.437771d0, 123.5968d0,
     7              -11.7273029d0, 55.782177d0,  20.2082d0,
     8               10.2742696d0,   .373813d0,  40.8226d0,
     9                6.4914588d0, 13.218362d0, 123.4722d0,
     O                5.8539148d0, 62.583231d0, 155.6977d0,
     1               -5.4872205d0, 63.593761d0, 184.6277d0,
     2               -5.4290191d0, 76.438310d0, 267.2772d0,
     3                5.1609570d0, 45.815258d0,  55.0196d0,
     4                5.0786314d0,  8.448301d0, 152.5268d0,
     5               -4.0735782d0, 56.792707d0,  49.1382d0,
     6                3.7227167d0, 49.747842d0, 204.6609d0,
     7                3.3971932d0, 12.058272d0,  56.5233d0,
     8               -2.8347004d0, 75.278220d0, 200.3284d0,
     9               -2.6550721d0, 65.241008d0, 201.6651d0,
     O               -2.5717867d0, 64.604291d0, 213.5577d0,
     1               -2.4712188d0,  1.647247d0,  17.0374d0,
     2                2.4625410d0,  7.811584d0, 164.4194d0,
     3                2.2464112d0, 12.207832d0,  94.5422d0,
     4               -2.0755511d0, 63.856665d0, 131.9124d0,
     5               -1.9713669d0, 56.155990d0,  61.0309d0,
     6               -1.8813061d0, 77.448840d0, 296.2073d0,
     7               -1.8468785d0,  6.801054d0, 135.4894d0,
     8                1.8186742d0, 62.209418d0, 114.8750d0,
     9                1.7601888d0, 20.656133d0, 247.0691d0,
     O               -1.5428851d0, 48.344406d0, 256.6114d0,
     1                1.4738838d0, 55.145460d0,  32.1008d0,
     2               -1.4593669d0, 69.000539d0, 143.6804d0,
     3                1.4192259d0, 11.071350d0,  16.8784d0,
     4               -1.1818980d0, 74.291298d0, 160.6835d0,
     5                1.1756474d0, 11.047742d0,  27.5932d0,
     6               -1.1316126d0,  0.636717d0, 348.1074d0,
     7                1.0896928d0, 12.844549d0,  82.6496d0/),(/3,47/) )
C**** Table 2 (4).  Precessional parameter: ECCEN sin(omega) (unused)
      REAL*8, PARAMETER, DIMENSION(3,46) :: TABL2 = RESHAPE(  (/
     1     .0186080D0,  54.646484D0,   32.012589D0,
     2     .0162752D0,  57.785370D0,  197.181274D0,
     3    -.0130066D0,  68.296539D0,  311.699463D0,
     4     .0098883D0,  67.659821D0,  323.592041D0,
     5    -.0033670D0,  67.286011D0,  282.769531D0,
     6     .0033308D0,  55.638351D0,   90.587509D0,
     7    -.0023540D0,  68.670349D0,  352.522217D0,
     8     .0014002D0,  76.656036D0,  131.835892D0,
     9     .0010070D0,  56.798447D0,  157.536392D0,
     O     .0008570D0,  66.649292D0,  294.662109D0,
     1     .0006499D0,  53.504456D0,  118.253082D0,
     2     .0005990D0,  67.023102D0,  335.484863D0,
     3     .0003780D0,  68.933258D0,  299.806885D0,
     4    -.0003370D0,  56.630219D0,  149.162415D0,
     5     .0003334D0,  86.256454D0,  283.915039D0,
     6     .0003334D0,  23.036499D0,  320.110107D0,
     7     .0002916D0,  89.395340D0,   89.083817D0,
     8     .0002916D0,  26.175385D0,  125.278732D0,
     9     .0002760D0,  69.307068D0,  340.629639D0,
     O    -.0002330D0,  99.906509D0,  203.602081D0,
     1    -.0002330D0,  36.686569D0,  239.796982D0,
     2     .0001820D0,  67.864838D0,  155.484787D0,
     3     .0001772D0,  99.269791D0,  215.494690D0,
     4     .0001772D0,  36.049850D0,  251.689606D0,
     5    -.0001740D0,  56.625275D0,  130.232391D0,
     6    -.0001240D0,  68.856720D0,  214.059708D0,
     7     .0001153D0,  87.266983D0,  312.845215D0,
     8     .0001153D0,  22.025970D0,  291.179932D0,
     9     .0001008D0,  90.405869D0,  118.013870D0,
     O     .0001008D0,  25.164856D0,   96.348694D0,
     1     .0000912D0,  78.818680D0,  160.318298D0,
     2     .0000912D0,  30.474274D0,   83.706894D0,
     3    -.0000806D0, 100.917038D0,  232.532120D0,
     4    -.0000806D0,  35.676025D0,  210.866943D0,
     5     .0000798D0,  81.957565D0,  325.487061D0,
     6     .0000798D0,  33.613159D0,  248.875565D0,
     7    -.0000638D0,  92.468735D0,   80.005234D0,
     8    -.0000638D0,  44.124329D0,    3.393823D0,
     9     .0000612D0, 100.280319D0,  244.424728D0,
     O     .0000612D0,  35.039322D0,  222.759552D0,
     1    -.0000603D0,  98.895981D0,  174.672028D0,
     2    -.0000603D0,  35.676025D0,  210.866943D0,
     3     .0000597D0,  87.248322D0,  342.489990D0,
     4     .0000597D0,  24.028381D0,   18.684967D0,
     5     .0000559D0,  86.630264D0,  324.737793D0,
     6     .0000559D0,  22.662689D0,  279.287354D0/), (/3,46/) )
C**** Table 3 (5).  Eccentricity: ECCEN (unused)
      REAL*8, PARAMETER, DIMENSION(3,42) :: TABL3 = RESHAPE( (/
     1     .01102940D0,   3.138886D0,  165.168686D0,
     2    -.00873296D0,  13.650058D0,  279.687012D0,
     3    -.00749255D0,  10.511172D0,  114.518250D0,
     4     .00672394D0,  13.013341D0,  291.579590D0,
     5     .00581229D0,   9.874455D0,  126.410858D0,
     6    -.00470066D0,   0.636717D0,  348.107422D0,
     7    -.00254464D0,  12.639528D0,  250.756897D0,
     8     .00231485D0,   0.991874D0,   58.574905D0,
     9    -.00221955D0,   9.500642D0,   85.588211D0,
     O     .00201868D0,   2.147012D0,  106.593765D0,
     1    -.00172371D0,   0.373813D0,   40.822647D0,
     2    -.00166112D0,  12.658154D0,  221.112030D0,
     3     .00145096D0,   1.010530D0,   28.930038D0,
     4     .00131342D0,  12.021467D0,  233.004639D0,
     5     .00101442D0,   0.373813D0,   40.822647D0,
     6    -.00088343D0,  14.023871D0,  320.509521D0,
     7    -.00083395D0,   6.277772D0,  330.337402D0,
     8     .00079475D0,   6.277772D0,  330.337402D0,
     9     .00067546D0,  27.300110D0,  199.373871D0,
     O    -.00066447D0,  10.884985D0,  155.340912D0,
     1     .00062591D0,  21.022339D0,  229.036499D0,
     2     .00059751D0,  22.009552D0,   99.823303D0,
     3    -.00053262D0,  27.300110D0,  199.373871D0,
     4    -.00052983D0,   5.641055D0,  342.229980D0,
     5    -.00052983D0,   6.914489D0,  318.444824D0,
     6     .00052836D0,  12.002811D0,  262.649414D0,
     7     .00051457D0,  16.788940D0,   84.855621D0,
     8    -.00050748D0,  11.647654D0,  192.181992D0,
     9    -.00049048D0,  24.535049D0,   75.027847D0,
     O     .00048888D0,  18.870667D0,  294.654541D0,
     1     .00046278D0,  26.026688D0,  223.159103D0,
     2     .00046212D0,   8.863925D0,   97.480820D0,
     3     .00046046D0,  17.162750D0,  125.678268D0,
     4     .00042941D0,   2.151964D0,  125.523788D0,
     5     .00042342D0,  37.174576D0,  325.784668D0,
     6     .00041713D0,  19.748917D0,  252.821732D0,
     7    -.00040745D0,  21.022339D0,  229.036499D0,
     8    -.00040569D0,   3.512699D0,  205.991333D0,
     9    -.00040569D0,   1.765073D0,  124.346024D0,
     O    -.00040385D0,  29.802292D0,   16.435165D0,
     1     .00040274D0,   7.746099D0,  350.172119D0,
     2     .00040068D0,   1.142024D0,  273.759521D0/), (/3,42/) )
C**** Table 4 (1).  Fundamental elements of the ecliptic: ECCEN sin(pi)
      REAL*8, PARAMETER, DIMENSION(3,19) :: TABL4 = RESHAPE( (/
     1     .01860798D0,   4.207205D0,   28.620089D0,
     2     .01627522D0,   7.346091D0,  193.788772D0,
     3    -.01300660D0,  17.857263D0,  308.307024D0,
     4     .00988829D0,  17.220546D0,  320.199637D0,
     5    -.00336700D0,  16.846733D0,  279.376984D0,
     6     .00333077D0,   5.199079D0,   87.195000D0,
     7    -.00235400D0,  18.231076D0,  349.129677D0,
     8     .00140015D0,  26.216758D0,  128.443387D0,
     9     .00100700D0,   6.359169D0,  154.143880D0,
     O     .00085700D0,  16.210016D0,  291.269597D0,
     1     .00064990D0,   3.065181D0,  114.860583D0,
     2     .00059900D0,  16.583829D0,  332.092251D0,
     3     .00037800D0,  18.493980D0,  296.414411D0,
     4    -.00033700D0,   6.190953D0,  145.769910D0,
     5     .00027600D0,  18.867793D0,  337.237063D0,
     6     .00018200D0,  17.425567D0,  152.092288D0,
     7    -.00017400D0,   6.186001D0,  126.839891D0,
     8    -.00012400D0,  18.417441D0,  210.667199D0,
     9     .00001250D0,   0.667863D0,   72.108838D0/), (/3,19/) )
C**** Table 5 (3).  General precession in longitude: psi
      REAL*8, PARAMETER, DIMENSION(3,78) :: TABL5 = RESHAPE( (/
     1    7391.0225890d0,  31.609974d0,   251.9025d0,
     2    2555.1526947d0,  32.620504d0,   280.8325d0,
     3    2022.7629188d0,  24.172203d0,   128.3057d0,
     4   -1973.6517951d0,   0.636717d0,   348.1074d0,
     5    1240.2321818d0,  31.983787d0,   292.7252d0,
     6     953.8679112d0,   3.138886d0,   165.1686d0,
     7    -931.7537108d0,  30.973257d0,   263.7951d0,
     8     872.3795383d0,  44.828336d0,    15.3747d0,
     9     606.3544732d0,   0.991874d0,    58.5749d0,
     O    -496.0274038d0,   0.373813d0,    40.8226d0,
     1     456.9608039d0,  43.668246d0,   308.4258d0,
     2     346.9462320d0,  32.246691d0,   240.0099d0,
     3    -305.8412902d0,  30.599444d0,   222.9725d0,
     4     249.6173246d0,   2.147012d0,   106.5937d0,
     5    -199.1027200d0,  10.511172d0,   114.5182d0,
     6     191.0560889d0,  42.681324d0,   268.7809d0,
     7    -175.2936572d0,  13.650058d0,   279.6869d0,
     8     165.9068833d0,   0.986922d0,    39.6448d0,
     9     161.1285917d0,   9.874455d0,   126.4108d0,
     O     139.7878093d0,  13.013341d0,   291.5795d0,
     1    -133.5228399d0,   0.262904d0,   307.2848d0,
     2     117.0673811d0,   0.004952d0,    18.9300d0,
     3     104.6907281d0,   1.142024d0,   273.7596d0,
     4      95.3227476d0,  63.219948d0,   143.8050d0,
     5      86.7824524d0,   0.205021d0,   191.8927d0,
     6      86.0857729d0,   2.151964d0,   125.5237d0,
     7      70.5893698d0,  64.230478d0,   172.7351d0,
     8     -69.9719343d0,  43.836462d0,   316.7998d0,
     9     -62.5817473d0,  47.439436d0,   319.6024d0,
     O      61.5450059d0,   1.384343d0,    69.7526d0,
     1     -57.9364011d0,   7.437771d0,   123.5968d0,
     2      57.1899832d0,  18.829299d0,   217.6432d0,
     3     -57.0236109d0,   9.500642d0,    85.5882d0,
     4     -54.2119253d0,   0.431696d0,   156.2147d0,
     5      53.2834147d0,   1.160090d0,    66.9489d0,
     6      52.1223575d0,  55.782177d0,    20.2082d0,
     7     -49.0059908d0,  12.639528d0,   250.7568d0,
     8     -48.3118757d0,   1.155138d0,    48.0188d0,
     9     -45.4191685d0,   0.168216d0,     8.3739d0,
     O     -42.2357920d0,   1.647247d0,    17.0374d0,
     1     -34.7971099d0,  10.884985d0,   155.3409d0,
     2      34.4623613d0,   5.610937d0,    94.1709d0,
     3     -33.8356643d0,  12.658184d0,   221.1120d0,
     4      33.6689362d0,   1.010530d0,    28.9300d0,
     5     -31.2521586d0,   1.983748d0,   117.1498d0,
     6     -30.8798701d0,  14.023871d0,   320.5095d0,
     7      28.4640769d0,   0.560178d0,   262.3602d0,
     8     -27.1960802d0,   1.273434d0,   336.2148d0,
     9      27.0860736d0,  12.021467d0,   233.0046d0,
     O     -26.3437456d0,  62.583231d0,   155.6977d0,
     1      24.7253740d0,  63.593761d0,   184.6277d0,
     2      24.6732126d0,  76.438310d0,   267.2772d0,
     3      24.4272733d0,   4.280910d0,    78.9281d0,
     4      24.0127327d0,  13.218362d0,   123.4722d0,
     5      21.7150294d0,  17.818769d0,   188.7132d0,
     6     -21.5375347d0,   8.359495d0,   180.1364d0,
     7      18.1148363d0,  56.792707d0,    49.1382d0,
     8     -16.9603104d0,   8.448301d0,   152.5268d0,
     9     -16.1765215d0,   1.978796d0,    98.2198d0,
     O      15.5567653d0,   8.863925d0,    97.4808d0,
     1      15.4846529d0,   0.186365d0,   221.5376d0,
     2      15.2150632d0,   8.996212d0,   168.2438d0,
     3      14.5047426d0,   6.771027d0,   161.1199d0,
     4     -14.3873316d0,  45.815258d0,    55.0196d0,
     5      13.1351419d0,  12.002811d0,   262.6495d0,
     6      12.8776311d0,  75.278220d0,   200.3284d0,
     7      11.9867234d0,  65.241008d0,   201.6651d0,
     8      11.9385578d0,  18.870667d0,   294.6547d0,
     9      11.7030822d0,  22.009553d0,    99.8233d0,
     O      11.6018181d0,  64.604291d0,   213.5577d0,
     1     -11.2617293d0,  11.498094d0,   154.1631d0,
     2     -10.4664199d0,   0.578834d0,   232.7153d0,
     3      10.4333970d0,   9.237738d0,   138.3034d0,
     4     -10.2377466d0,  49.747842d0,   204.6609d0,
     5      10.1934446d0,   2.147012d0,   106.5938d0,
     6     -10.1280191d0,   1.196895d0,   250.4676d0,
     7      10.0289441d0,   2.133898d0,   332.3345d0,
     8     -10.0034259d0,   0.173168d0,    27.3039d0/), (/3,78/) )
C****
      REAL*8 :: YM1950,SUMC,ARG,ESINPI,ECOSPI,PIE,PSI,FSINFD
      INTEGER :: I
C****
      YM1950 = YEAR-1950.
C****
C**** Obliquity from Table 1 (2):
C****   OBLIQ# = 23.320556 (degrees)             Equation 5.5 (15)
C****   OBLIQ  = OBLIQ# + sum[A cos(ft+delta)]   Equation 1 (5)
C****
      SUMC = 0.
      DO I=1,47
        ARG  = PI180*(YM1950*TABL1(2,I)/3600.+TABL1(3,I))
        SUMC = SUMC + TABL1(1,I)*COS(ARG)
      END DO
      OBLIQ = 23.320556D0 + SUMC/3600.
!      OBLIQ  = OBLIQ*PI180 ! not needed for output in degrees
C****
C**** Eccentricity from Table 4 (1):
C****   ECCEN sin(pi) = sum[M sin(gt+beta)]           Equation 4 (1)
C****   ECCEN cos(pi) = sum[M cos(gt+beta)]           Equation 4 (1)
C****   ECCEN = ECCEN sqrt[sin(pi)^2 + cos(pi)^2]
C****
      ESINPI = 0.
      ECOSPI = 0.
      DO I=1,19
        ARG    = PI180*(YM1950*TABL4(2,I)/3600.+TABL4(3,I))
        ESINPI = ESINPI + TABL4(1,I)*SIN(ARG)
        ECOSPI = ECOSPI + TABL4(1,I)*COS(ARG)
      END DO
      ECCEN  = SQRT(ESINPI*ESINPI+ECOSPI*ECOSPI)
C****
C**** Perihelion from Equation 4,6,7 (9) and Table 4,5 (1,3):
C****   PSI# = 50.439273 (seconds of degree)         Equation 7.5 (16)
C****   ZETA =  3.392506 (degrees)                   Equation 7.5 (17)
C****   PSI = PSI# t + ZETA + sum[F sin(ft+delta)]   Equation 7 (9)
C****   PIE = atan[ECCEN sin(pi) / ECCEN cos(pi)]
C****   OMEGVP = PIE + PSI + 3.14159                 Equation 6 (4.5)
C****
      PIE = ATAN2(ESINPI,ECOSPI)
      FSINFD = 0.
      DO I=1,78
        ARG    = PI180*(YM1950*TABL5(2,I)/3600.+TABL5(3,I))
        FSINFD = FSINFD + TABL5(1,I)*SIN(ARG)
      END DO
      PSI    = PI180*(3.392506D0+(YM1950*50.439273D0+FSINFD)/3600.)
      OMEGVP = MOD(PIE+PSI+.5*TWOPI,TWOPI)
      IF(OMEGVP.lt.0.)  OMEGVP = OMEGVP + TWOPI
      OMEGVP = OMEGVP/PI180  ! for output in degrees
C****
      RETURN
      END SUBROUTINE ORBPAR

      SUBROUTINE ORBIT (DOBLIQ,ECCEN,DOMEGVP,VEDAY,EDPY, DAY,
     *                  SDIST,SIND,COSD,SUNLON,SUNLAT,EQTIME)
C****
C**** ORBIT receives orbital parameters and time of year, and returns
C**** distance from Sun, declination angle, and Sun's overhead position.
C**** Reference for following caculations is:  V.M.Blanco and
C**** S.W.McCuskey, 1961, "Basic Physics of the Solar System", pages
C**** 135 - 151.  Existence of Moon and heavenly bodies other than
C**** Earth and Sun are ignored.  Earth is assumed to be spherical.
C****
C**** Program author: Gary L. Russell 2004/11/16
C**** Angles, longitude and latitude are measured in radians.
C****
C**** Input: ECCEN  = eccentricity of the orbital ellipse
C****        OBLIQ  = latitude of Tropic of Cancer
C****        OMEGVP = longitude of perihelion (sometimes Pi is added) =
C****               = spatial angle from vernal equinox to perihelion
C****                 with Sun as angle vertex
C****        DAY    = days measured since 2000 January 1, hour 0
C****
C****        EDPY  = Earth days per year
C****                tropical year = 365.2425 (Gregorgian Calendar)
C****                tropical year = 365      (Generic Year)
C****        VEDAY = Vernal equinox
C****                79.0 (Generic year Mar 21 hour 0)
C****                79.5 (Generic year Mar 21 hour 12 - PMIP standard)
C****                79.3125d0 for days from 2000 January 1, hour 0 till vernal
C****                     equinox of year 2000 = 31 + 29 + 19 + 7.5/24
C****
C**** Intermediate quantities:
C****    BSEMI = semi minor axis in units of semi major axis
C****   PERIHE = perihelion in days since 2000 January 1, hour 0
C****            in its annual revolution about Sun
C****       TA = true anomaly = spatial angle from perihelion to
C****            current location with Sun as angle vertex
C****       EA = eccentric anomaly = spatial angle measured along
C****            eccentric circle (that circumscribes Earth's orbit)
C****            from perihelion to point above (or below) Earth's
C****            absisca (where absisca is directed from center of
C****            eccentric circle to perihelion)
C****       MA = mean anomaly = temporal angle from perihelion to
C****            current time in units of 2*Pi per tropical year
C****   TAofVE = TA(VE) = true anomaly of vernal equinox = - OMEGVP
C****   EAofVE = EA(VE) = eccentric anomaly of vernal equinox
C****   MAofVE = MA(VE) = mean anomaly of vernal equinox
C****   SLNORO = longitude of Sun in Earth's nonrotating reference frame
C****   VEQLON = longitude of Greenwich Meridion in Earth's nonrotating
C****            reference frame at vernal equinox
C****   ROTATE = change in longitude in Earth's nonrotating reference
C****            frame from point's location on vernal equinox to its
C****            current location where point is fixed on rotating Earth
C****   SLMEAN = longitude of fictitious mean Sun in Earth's rotating
C****            reference frame (normal longitude and latitude)
C****
C**** Output: SIND = sine of declination angle = sin(SUNLAT)
C****         COSD = cosine of the declination angle = cos(SUNLAT)
C****       SUNDIS = distance to Sun in units of semi major axis
C****       SUNLON = longitude of point on Earth directly beneath Sun
C****       SUNLAT = latitude of point on Earth directly beneath Sun
C****       EQTIME = Equation of Time =
C****              = longitude of fictitious mean Sun minus SUNLON
C****
C**** From the above reference:
C**** (4-54): [1 - ECCEN*cos(EA)]*[1 + ECCEN*cos(TA)] = (1 - ECCEN^2)
C**** (4-55): tan(TA/2) = sqrt[(1+ECCEN)/(1-ECCEN)]*tan(EA/2)
C**** Yield:  tan(EA) = sin(TA)*sqrt(1-ECCEN^2) / [cos(TA) + ECCEN]
C****    or:  tan(TA) = sin(EA)*sqrt(1-ECCEN^2) / [cos(EA) - ECCEN]
C****
      USE CONSTANT, only : twopi,pi,radian
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: DOBLIQ,ECCEN,DOMEGVP,DAY,VEDAY,EDPY
      REAL*8, INTENT(OUT) :: SIND,COSD,SDIST,SUNLON,SUNLAT,EQTIME

      REAL*8 MA,OMEGVP,OBLIQ,EA,DEA,BSEMI
     *     ,TAofVE,EAofVE,MAofVE,SUNDIS,TA,SUNX,SUNY,SLNORO
     *     ,VEQLON,ROTATE,SLMEAN
c      REAL*8, PARAMETER :: EDAYzY=365.2425d0, VE2000=79.3125d0
c      REAL*8, PARAMETER :: EDAYzY=365d0, VE2000=79d0  ! original parameters
      REAL*8  EDAYzY,VE2000
C****
      VE2000=VEDAY
      EDAYzY=EDPY
      OMEGVP=DOMEGVP*radian
      OBLIQ=DOBLIQ*radian
C**** Determine EAofVE from geometry: tan(EA) = b*sin(TA) / [e+cos(TA)]
C**** Determine MAofVE from Kepler's equation: MA = EA - e*sin(EA)
C**** Determine MA knowing time from vernal equinox to current day
C****
      BSEMI  = SQRT (1 - ECCEN*ECCEN)
      TAofVE = - OMEGVP
      EAofVE = ATAN2 (BSEMI*SIN(TAofVE), ECCEN+COS(TAofVE))
      MAofVE = EAofVE - ECCEN*SIN(EAofVE)
C     PERIHE = VE2000 - MAofVE*EDAYzY/TWOPI
      MA     = MODULO (TWOPI*(DAY-VE2000)/EDAYzY + MAofVE, TWOPI)
C****
C**** Numerically invert Kepler's equation: MA = EA - e*sin(EA)
C****
      EA  = MA + ECCEN*(SIN(MA) + ECCEN*SIN(2*MA)/2)
   10 dEA = (MA - EA + ECCEN*SIN(EA)) / (1 - ECCEN*COS(EA))
      EA  = EA + dEA
      IF(ABS(dEA).gt.1d-10)  GO TO 10
C****
C**** Calculate distance to Sun and true anomaly
C****
      SUNDIS = 1 - ECCEN*COS(EA)
      TA     = ATAN2 (BSEMI*SIN(EA), COS(EA)-ECCEN)
      SDIST  = SUNDIS*SUNDIS   ! added for compatiblity
C****
C**** Change reference frame to be nonrotating reference frame, angles
C**** fixed according to stars, with Earth at center and positive x
C**** axis be ray from Earth to Sun were Earth at vernal equinox, and
C**** x-y plane be Earth's equatorial plane.  Distance from current Sun
C**** to this x axis is SUNDIS sin(TA-TAofVE).  At vernal equinox, Sun
C**** is located at (SUNDIS,0,0).  At other times, Sun is located at:
C****
C**** SUN = (SUNDIS cos(TA-TAofVE),
C****        SUNDIS sin(TA-TAofVE) cos(OBLIQ),
C****        SUNDIS sin(TA-TAofVE) sin(OBLIQ))
C****
      SIND   = SIN(TA-TAofVE) * SIN(OBLIQ)
      COSD   = SQRT (1 - SIND*SIND)
      SUNX   = COS(TA-TAofVE)
      SUNY   = SIN(TA-TAofVE) * COS(OBLIQ)
      SLNORO = ATAN2 (SUNY,SUNX)
C****
C**** Determine Sun location in Earth's rotating reference frame
C**** (normal longitude and latitude)
C****
      VEQLON = TWOPI*VE2000 - PI + MAofVE - TAofVE  !  modulo 2*Pi
      ROTATE = TWOPI*(DAY-VE2000)*(EDAYzY+1)/EDAYzY
      SUNLON = MODULO (SLNORO-ROTATE-VEQLON, TWOPI)
      IF(SUNLON.gt.PI)  SUNLON = SUNLON - TWOPI
      SUNLAT = ASIN (SIN(TA-TAofVE)*SIN(OBLIQ))
C****
C**** Determine longitude of fictitious mean Sun
C**** Calculate Equation of Time
C****
      SLMEAN = PI - TWOPI*(DAY-FLOOR(DAY))
      EQTIME = MODULO (SLMEAN-SUNLON, TWOPI)
      IF(EQTIME.gt.PI)  EQTIME = EQTIME - TWOPI
C****
      RETURN
      END SUBROUTINE ORBIT

      subroutine updBCd (year)
!@sum  reads appropriate Black Carbon deposition data if necessary
!@auth R. Ruedy
!@ver  1.0
      USE FILEMANAGER
      USE RAD_COM, only : depoBC
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      implicit none
      integer, intent(in)   :: year

      integer,parameter :: imr=72,jmr=46
      real*8  BCdep1(imr,jmr),BCdep2(imr,jmr),wt   ! to limit i/o

      integer :: iu,year1,year2,year0,yearL=-2,year_old=-1
      save       iu,year1,year2,year0,yearL,   year_old,BCdep1,BCdep2

      character*80 title
      real*4 BCdep4(imr,jmr)

C**** check whether update is needed
      if (year.eq.year_old) return
      if (year_old.eq.yearL.and.year.gt.yearL) return
      if (year_old.eq.year0.and.year.lt.year0) return

      call openunit('BC_dep',iu,.true.,.true.)

      if (year_old.lt.0) then
C****   read whole input file and find range: year0->yearL
   10   read(iu,end=20) title
        read(title,*) yearL
        go to 10
      end if

   20 rewind (iu)
      read(iu) title,BCdep4
      read(title,*) year0
      BCdep1=BCdep4 ; BCdep2=BCdep4 ; year2=year0 ; year1=year0
      if (year.le.year1)              year2=year+1

      do while (year2.lt.year .and. year2.ne.yearL)
         year1 = year2 ; BCdep1 = BCdep2
         read (iu) title,BCdep4
         read(title,*) year2
         BCdep2 = BCdep4
      end do

      if(year.le.year1) then
        wt = 0.
      else if (year.ge.yearL) then
        wt = 1.
      else
        wt = (year-year1)/(real(year2-year1,kind=8))
      end if

      if (AM_I_ROOT())  write(6,*)
     &     'Using BCdep data from year',year1+wt*(year2-year1)
      call closeunit(iu)

C**** Set the Black Carbon deposition array
      depoBC(:,:) = BCdep1(:,:) + wt*(BCdep2(:,:)-BCdep1(:,:))

      year_old = year

      return
      end subroutine updBCd
#ifdef HEALY_LM_DIAGS
      real*8 function Fe(M,N)
      real*8 M,N

      Fe=0.47d0*log(1.+2.01d-5*(M*N)**(0.75)
     .   +5.31d-15*M*(M*N)**(1.52)) 
      return
      end
#endif

