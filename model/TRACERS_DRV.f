#include "rundeck_opts.h"
!@sum  TRACERS_DRV: tracer-dependent routines for air/water mass
!@+    and ocean tracers
!@+    Routines included:
!@+      Those that MUST EXIST for all tracers:
!@+        Diagnostic specs: init_tracer
!@+        Tracer initialisation + sources: tracer_ic, set_tracer_source
!@+        Entry points: daily_tracer
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      SUBROUTINE init_tracer
!@sum init_tracer initializes trace gas attributes and diagnostics
!@auth J. Lerner
!@calls sync_param, SET_TCON, RDLAND, RDDRYCF
      USE DOMAIN_DECOMP_ATM, only:GRID,GET,AM_I_ROOT,
     &     write_parallel,readt8_parallel
      USE CONSTANT, only: mair,mwat,sday
#ifdef TRACERS_AEROSOLS_SOA
     &                   ,gasc
#endif  /* TRACERS_AEROSOLS_SOA */
      USE MODEL_COM, only: dtsrc,byim,lm,jm,itime,pmidl00,nisurf
      USE GEOM, only: axyp,byaxyp
      USE DYNAMICS, only: am,byam  ! Air mass of each box (kg/m^2)
      USE TRACER_COM
#ifdef TRACERS_ON
      USE TRDIAG_COM
#endif
      USE PARAM
#ifdef TRACERS_SPECIAL_Lerner
      USE TRACERS_MPchem_COM, only: n_MPtable,tcscale
!@dbparam dsol describes portion of solar cycle being modeled for linoz
!@+      +1.0 = solar max, 0.0 = neutral, -1.0 = solar min
      USE LINOZ_CHEM_COM, only: dsol
#endif
#ifdef TRACERS_WATER
      USE LANDICE_COM, only : trli0    ! should these be in tracer_com?
      USE SEAICE_COM, only : trsi0
#ifdef TRDIAG_WETDEPO
      USE CLOUDS, ONLY : diag_wetdep
#endif
#endif /* TRACERS_WATER */
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      use tracers_dust,only : imDust,prefDustSources,fracClayPDFscheme
     &   ,fracSiltPDFscheme
      use trdust_drv, only : init_soildust
#endif
#ifdef TRACERS_QUARZHEM
     &     ,DensityHematite, DensityQuartz, FreeFe, frHemaInQuarAggr,
     &     pureByTotalHematite
#endif
#ifdef TRACERS_SPECIAL_Shindell
      use tracer_sources, only: aircraft_Tyr1,aircraft_Tyr2
#ifdef GFED_3D_BIOMASS
     & ,biomass_Tyr1,biomass_Tyr2
#endif
      USE TRCHEM_Shindell_COM,only:LCOalt,PCOalt,
     &     CH4altINT,CH4altINX,LCH4alt,PCH4alt,checktracer_on,
     &     CH4altX,CH4altT,ch4_init_sh,ch4_init_nh,scale_ch4_IC_file,
     &     OxICIN,OxIC,OxICINL,OxICL,
     &     fix_CH4_chemistry,which_trop,PI_run,PIratio_N,PIratio_CO_T,
     &     PIratio_CO_S,PIratio_other,allowSomeChemReinit,
     &     CH4ICIN,CH4ICX,CH4ICINL,CH4ICL,rad_FL,use_rad_ch4,
     &     COICIN,COIC,COICINL,COICL,Lmax_rad_O3,Lmax_rad_CH4
#ifdef SHINDELL_STRAT_CHEM
     &     ,BrOxaltIN,ClOxaltIN,ClONO2altIN,HClaltIN,BrOxalt,
     &     ClOxalt,ClONO2alt,HClalt,N2OICIN,N2OICX,N2OICINL,N2OICL,
     &     CFCICIN,CFCIC,CFCICINL,CFCICL,PIratio_N2O,PIratio_CFC,
     &     use_rad_n2o,use_rad_cfc,cfc_rad95,PltOx,Tpsc_offset_N,
     &     Tpsc_offset_S
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      USE TRACER_SOURCES, only:int_wet_dist,topo_lim,sat_lim,gw_ulim,
     &gw_llim,sw_lim,exclude_us_eu,nn_or_zon,ice_age,nday_ch4,max_days,
     &ns_wet,nra_ch4
#endif
#ifdef BIOGENIC_EMISSIONS
      use biogenic_emis, only: base_isopreneX
#endif
#endif /* TRACERS_SPECIAL_Shindell */
#ifdef TRACERS_AEROSOLS_SOA
      USE TRACERS_SOA, only: n_soa_i,n_soa_e,soa_init
#endif  /* TRACERS_AEROSOLS_SOA */
#if (defined TRACERS_COSMO)
      USE COSMO_SOURCES, only: be7_src_param
#endif
#if (defined TRACERS_AMP)
      USE AERO_PARAM, only: DG_DD1, DG_DD2, DG_AKK,
     & DG_DS1, DG_DS2, DG_SSA, DG_SSC, DG_ACC,
     & DG_SSS, DG_OCC, DG_BC1, DG_BC2, DG_BC3,
     & DG_DBC, DG_BOC, DG_BCS, DG_OCS, DG_MXX,

     & SOLU_DD1, SOLU_DD2, SOLU_AKK, SOLU_ACC,
     & SOLU_DS1, SOLU_DS2, SOLU_SSA, SOLU_SSC,
     & SOLU_SSS, SOLU_OCC, SOLU_BC1, SOLU_BC2, SOLU_BC3,
     & SOLU_DBC, SOLU_BOC, SOLU_BCS, SOLU_OCS, SOLU_MXX

      USE AERO_ACTV, only: DENS_SULF, DENS_DUST,
     &          DENS_SEAS, DENS_BCAR, DENS_OCAR
      USE AERO_CONFIG, only: nbins
      USE AMP_AEROSOL, only: AMP_DIAG_FC, AMP_RAD_KEY
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      USE obio_forc, only : atmCO2
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      USE AEROSOL_SOURCES, only: tune_ss1, tune_ss2, om2oc, BBinc
#endif

C**** RDF: TES variables to retrieve from rundeck
#ifdef TES_SIM
      USE TES_SIMULATOR,only: TESObsDiagnosticsOn,TESObsMinDegFr,
     &    TESObsDataDir,TESCatDiagnosticsOn,TESCatMinQuality,
     &    TESCatUseColloc,TESCatVariableHDOAKOnly
#endif

      USE FILEMANAGER, only: openunit,closeunit,nameunit
      implicit none
      integer :: l,k,n,kr,m,ns
#ifdef TRACERS_SPECIAL_O18
      real*8 fracls
#endif
#if (defined TRACERS_WATER) || (defined TRACERS_DRYDEP)
!@param convert_HSTAR converts from mole/Joule to mole/(L*atm)
      real*8, parameter :: convert_HSTAR = 1.01325d2
#endif
#ifdef TRACERS_SPECIAL_Shindell
!@var PRES local nominal pressure for vertical interpolations
!@var iu_data unit number
!@var title header read in from file
      REAL*8, DIMENSION(LM) :: PRES
      integer iu_data,i,j,nq
      character*80 title
      character(len=300) :: out_line
#ifdef SHINDELL_STRAT_CHEM
      real*8, dimension(6) :: temp_ghg
      integer :: temp_year
#endif /* SHINDELL_STRAT_CHEM */
#endif /* TRACERS_SPECIAL_Shindell */

#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_GASEXCH_ocean_CFC)
      integer i, iu_data
#endif
#if (!defined(TRACERS_GASEXCH_ocean_CO2)) && defined(TRACERS_GASEXCH_land_CO2)
      real*8 :: atmCO2 = 280.d0
#endif

      INTEGER J_0, J_1, I_0, I_1
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Set defaults for tracer attributes (all dimensioned ntm)

      itime_tr0=itime
      t_qlimit = .true.
      trdecay = 0.
      do_fire = .false.
#ifdef TRACERS_ON
      needtrs = .false.
      ntsurfsrc = 0
      trpdens = 0.
      trradius = 0.
#ifdef TRACERS_SPECIAL_Lerner
      n_MPtable = 0
      tcscale = 0.
#endif
#endif /* TRACERS_ON */
      tr_wd_TYPE = nGas       !other options are nPART or nWATER
      tr_RKD = 0.
      tr_DHD = 0.
      fq_aer = 0.
      isdust = 0
#ifdef TRACERS_WATER
      trli0 = 0.
      trsi0 = 0.
      tr_H2ObyCH4 = 0.
      dowetdep = .false.
#endif
#ifdef TRACERS_SPECIAL_O18
      iso_index = 1             ! act like water by default
#endif
#ifdef TRACERS_DRYDEP
      dodrydep = .false.
      F0 = 0.
      HSTAR = 0.   ! tr_RKD * convert_HSTAR
#endif
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
      trw0 = 0.
      ntrocn = 0
#endif
#ifdef TRACERS_OCEAN
      trglac = 0.
#endif

C**** Synchronise tracer related paramters from rundeck

C**** Get itime_tr0 from rundeck if it exists
      call sync_param("itime_tr0",itime_tr0,ntm)
#ifdef TRACERS_WATER
C**** Decide on water tracer conc. units from rundeck if it exists
      call sync_param("to_per_mil",to_per_mil,ntm)
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      call sync_param("tune_ss1",tune_ss1)
      call sync_param("tune_ss2",tune_ss2)
      call sync_param("BBinc",BBinc)
#endif
#ifdef TRACERS_SPECIAL_O18
C**** set super saturation parameter for isotopes if needed
      call sync_param("supsatfac",supsatfac)
#if (defined TRACERS_OCEAN)
      call sync_param ( "water_tracer_ic",water_tracer_ic  )
#endif

#endif

C**** RDF: Retrieve rundeck parameters for TES
#ifdef TES_SIM
      call sync_param("TESObsDiagnosticsOn",TESObsDiagnosticsOn)
      call sync_param("TESObsMinDegFr",TESObsMinDegFr)
      call sync_param("TESObsDataDir",TESObsDataDir)
      call sync_param("TESCatDiagnosticsOn",TESCatDiagnosticsOn)
      call sync_param("TESCatMinQuality",TESCatMinQuality)
      call sync_param("TESCatUseColloc",TESCatUseColloc)
      call sync_param("TESCatVariableHDOAKOnly",TESCatVariableHDOAKOnly)
C      call sync_param("ISCCPncol",ncol)

#endif

#ifdef TRACERS_ON
      CALL sync_param("diag_rad",diag_rad)
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
      CALL sync_param("diag_wetdep",diag_wetdep)
#endif
!not params call sync_param("trans_emis_overr_day",trans_emis_overr_day)
!not params call sync_param("trans_emis_overr_yr", trans_emis_overr_yr )
#endif /* TRACERS_ON */
#ifdef TRACERS_SPECIAL_Shindell
      call sync_param("allowSomeChemReinit",allowSomeChemReinit)
      call sync_param("which_trop",which_trop)
      call sync_param("PI_run",PI_run)
      call sync_param("PIratio_N",PIratio_N)
      call sync_param("PIratio_CO_T",PIratio_CO_T)
      call sync_param("PIratio_CO_S",PIratio_CO_S)
      call sync_param("PIratio_other",PIratio_other)
      call sync_param("rad_FL",rad_fl)
      call sync_param("checktracer_on",checktracer_on)
      call sync_param("use_rad_ch4",use_rad_ch4)
      call sync_param("Lmax_rad_O3",Lmax_rad_O3)
      call sync_param("Lmax_rad_CH4",Lmax_rad_CH4)
      call sync_param("aircraft_Tyr1",aircraft_Tyr1)
      call sync_param("aircraft_Tyr2",aircraft_Tyr2)
#ifdef GFED_3D_BIOMASS
      call sync_param("biomass_Tyr1",biomass_Tyr1)
      call sync_param("biomass_Tyr2",biomass_Tyr2)
#endif
#ifdef SHINDELL_STRAT_CHEM
      call sync_param("use_rad_n2o",use_rad_n2o)
      call sync_param("use_rad_cfc",use_rad_cfc)
      call sync_param("PIratio_N2O",PIratio_N2O)
      call sync_param("PIratio_CFC",PIratio_CFC)
      call sync_param("PltOx",PltOx)
      call sync_param("Tpsc_offset_N",Tpsc_offset_N)
      call sync_param("Tpsc_offset_S",Tpsc_offset_S)
#endif
#ifdef BIOGENIC_EMISSIONS
      call sync_param("base_isopreneX",base_isopreneX)
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      call sync_param("ice_age",ice_age)
      call sync_param("ns_wet",ns_wet)
      call sync_param("int_wet_dist",int_wet_dist)
      call sync_param("topo_lim",topo_lim)
      call sync_param("sat_lim",sat_lim)
      call sync_param("gw_ulim",gw_ulim)
      call sync_param("gw_llim",gw_llim)
      call sync_param("sw_lim",sw_lim)
      call sync_param("exclude_us_eu",exclude_us_eu)
      call sync_param("nn_or_zon",nn_or_zon)
      do n=1,nra_ch4
       if(nday_ch4(n) > max_days .or. nday_ch4(n) < 1)
     & call stop_model('nday_ch4 out of range',255)
      end do
#endif
      PRES(1:LM)=PMIDL00(1:LM)
#endif /* TRACERS_SPECIAL_Shindell */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP)
C**** DMS, seasalt from offline fields
      call sync_param("OFFLINE_DMS_SS",OFFLINE_DMS_SS)
C**** seasalt from offline fields
      call sync_param("OFFLINE_SS",OFFLINE_SS)
C**** decide on emissions
      call sync_param("imAER",imAER)
C**** decide if preindustrial emissions
      call sync_param("imPI",imPI)
C**** determine year of emissions
      call sync_param("aer_int_yr",aer_int_yr)
#endif
#if (defined TRACERS_AMP)
C**** Decide on how many times Radiation is called for aerosols once or nmode, default one call
      call sync_param("AMP_DIAG_FC",AMP_DIAG_FC)
C**** Decide Radiative Mixing Rules - Volume - Core Shell - Maxwell Garnett, default Volume
      call sync_param("AMP_RAD_KEY",AMP_RAD_KEY)
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
C**** decide on AEROCOM or interactive emissions
      CALL sync_param('imDUST',imDUST)
      call sync_param('prefDustSources', prefDustSources)
      call sync_param('fracClayPDFscheme', fracClayPDFscheme)
      call sync_param('fracSiltPDFscheme', fracSiltPDFscheme)
#endif
#ifdef TRACERS_QUARZHEM
      CALL sync_param('FreeFe',FreeFe)
      call sync_param( 'frHemaInQuarAggr', frHemaInQuarAggr )
      call sync_param( 'pureByTotalHematite', pureByTotalHematite )
#endif

#if (defined TRACERS_COSMO)
C**** get rundeck parameter for cosmogenic source factor
      call sync_param("be7_src_param", be7_src_param)
#endif
      call sync_param("no_emis_over_ice",no_emis_over_ice)

! call routine to read/set up regions and sectors for emissions:
      call setup_emis_sectors_regions

! initialize 3D source factors:
      ef_fact3d(:,:)=1.d0

! The following section will check for rundeck file of
! the form: trname_01, trname_02... and thereby define
! the ntsurfsrc(n). If those files exist it reads an
! 80 char header to get information including the
! source name (ssame-->{sname,lname,etc.}. ntsurfsrc(n)
! get set to zero if those files aren't found:
! (I can enclose this in an ifdef if it causes problems
! for people). num_srf_sources routine also assigns
! sources to sectors, if desired:
      nBBsources(:)=0
      do n=1,ntm
! general case:
#ifdef TRACERS_AMP
        call num_srf_sources(n,ntsurfsrc(n),.false.)
#else
        call num_srf_sources(n,ntsurfsrc(n),.true.)
#endif

#ifdef DYNAMIC_BIOMASS_BURNING
! allow some tracers to have biomass burning based on fire model:
        select case (trname(n))
        case('NOx','CO','Alkenes','Paraffin','BCB','OCB','NH3','SO2'
#ifdef TRACERS_SPECIAL_Shindell
     &  ,'CH4' ! in here to avoid potential Lerner tracers conflict
#endif
     &  )
          do_fire(n)=.true.
        end select
#endif /* DYNAMIC_BIOMASS_BURNING */

! allow some tracers to have biomass burning sources that mix over
! PBL layers (these become 3D sources no longer within ntsurfsrc(n)):
        select case (trname(n))
        case ('Alkenes', 'CO', 'NOx', 'Paraffin',
#ifdef TRACERS_SPECIAL_Shindell
     &  'CH4', ! in here to avoid potential Lerner tracers conflict
#endif
     &  'NH3', 'SO2', 'BCB', 'OCB', ! do not include sulfate here
     &  'M_BC1_BC', 'M_OCC_OC', 'M_BOC_BC', 'M_BOC_OC')
          call sync_param(trim(trname(n))//"_nBBsources",
     &    nBBsources(n))
          if(nBBsources(n)>0)then
            if(do_fire(n))then
              if(am_i_root())write(6,*)
     &        'nBBsource>0 for ',trim(trname(n)),' do_fire=t'
              call stop_model('nBBsource do_fire conflict',13)
            else
              ntsurfsrc(n)=ntsurfsrc(n)-nBBsources(n)
            end if
          end if
        end select
        if(do_fire(n) .and.  (ntsurfsrc(n)+1 > ntsurfsrcmax))then
          write(6,*)trname(n),'ntsurfsrc+1 > max of ',ntsurfsrcmax
          call stop_model('do_fire+ntsurfsrc too large',13)
        end if
        if(ntsurfsrc(n)+nBBsources(n) > ntsurfsrcmax)then
          write(6,*)trname(n),'ntsurfsrc+nBBsources > max of ',
     &    ntsurfsrcmax
          call stop_model('ntsurfsrc+nBBsources too large',13)
        end if

! other special cases:
#ifndef TRACERS_AEROSOLS_SOA
#ifdef TRACERS_AMP
        select case (trname(n))
        case ('M_OCC_OC') ! this handles OCT_src (terpene source)
          ntsurfsrc(n)=ntsurfsrc(n)+1
          ssname(n,ntsurfsrc(n))="Terpene_source"
        end select
#endif
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_SPECIAL_Shindell
        if(trname(n)=='CH4' .and. use_rad_ch4/=0) ntsurfsrc(n)=0
#endif
#if (defined TRACERS_SPECIAL_Shindell) && (defined GFED_3D_BIOMASS)
        call check_GFED_sectors(n) ! possible special 3D source sector
#endif
      enddo

C**** Define individual tracer characteristics
      do n=1,ntm

      select case (trname(n))

#ifdef TRACERS_ON
      case ('Air')
      n_Air = n
          ntm_power(n) = -2
          tr_mm(n) = mair

      case ('CO2n')
      n_CO2n = n
          ntm_power(n) = -6
          tr_mm(n) = 44.d0     !grams
          t_qlimit(n) = .false.
          ntsurfsrc(n) = 1
          needtrs(n)=.true.

      case ('CFCn')
      n_CFCn = n
          ntm_power(n) = -12
          !!tr_mm(n) = 136.0d0    !NCAR value
          tr_mm(n) = 137.37d0     !note units are in gr
          ntsurfsrc(n) = 1
          needtrs(n)=.true.

      case ('SF6')
      n_SF6 = n
          ntm_power(n) = -14
          tr_mm(n) = 146.01d0
          ntsurfsrc(n) = 1

      case ('Rn222')
      n_Rn222 = n
          ntm_power(n) = -21
          tr_mm(n) = 222.d0
          trdecay(n) = 2.1d-6
          ntsurfsrc(n) = 1

      case ('CO2')
      n_CO2 = n
          ntm_power(n) = -6
          tr_mm(n) = 44.d0
          t_qlimit(n) = .false.
          ntsurfsrc(n) = 6

      case ('N2O')
      n_N2O = n
          ntm_power(n) = -9
          tr_mm(n) = 44.d0
#ifdef TRACERS_SPECIAL_Lerner
          ntsurfsrc(n) = 1
          n_MPtable(n) = 1
          tcscale(n_MPtable(n)) = 1.
#endif
#ifdef SHINDELL_STRAT_CHEM
          call openunit('N2O_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),N2OICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           N2OICINL(:)=N2OICIN(i,j,:) ! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,N2OICINL,LM,PRES,N2OICL,.true.)
           N2OICX(I,J,:)=N2OICL(:)*am(:,i,j)*axyp(i,j)
          end do     ; end do
#endif

      case ('CFC11')   !!! should start April 1
      n_CFC11 = n
          ntm_power(n) = -12
          tr_mm(n) = 137.4d0
          ntsurfsrc(n) = 1
#ifdef TRACERS_SPECIAL_Lerner
          n_MPtable(n) = 2
          tcscale(n_MPtable(n)) = 1.
#endif

      case ('14CO2')   !!! should start 10/16
      n_14CO2 = n
          ntm_power(n) = -18
          tr_mm(n) = 46.d0
          ntsurfsrc(n) = 1

      case ('CH4')
      n_CH4 = n
          tr_mm(n) = 16.d0
#ifdef TRACERS_SPECIAL_Lerner
          ntsurfsrc(n) = 14
          ntm_power(n) = -9
          n_MPtable(n) = 3
          tcscale(n_MPtable(n)) = 1.
#endif
#ifdef TRACERS_SPECIAL_Shindell
C**** determine initial CH4 distribution if set from rundeck
C**** This is only effective with a complete restart.
          call sync_param("ch4_init_sh",ch4_init_sh)
          call sync_param("ch4_init_nh",ch4_init_nh)
          call sync_param("fix_CH4_chemistry",fix_CH4_chemistry)
          call sync_param("scale_ch4_IC_file",scale_ch4_IC_file)
          ntm_power(n) = -8
C         Interpolate CH4 altitude-dependence to model resolution:
          CALL LOGPINT(LCH4alt,PCH4alt,CH4altINT,LM,PRES,CH4altT,.true.)
          CALL LOGPINT(LCH4alt,PCH4alt,CH4altINX,LM,PRES,CH4altX,.true.)
          if(fix_CH4_chemistry.eq.-1)then
            call openunit('CH4_IC',iu_data,.true.,.true.)
            CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),
     &           CH4ICIN,0)
            call closeunit(iu_data)
            do j=J_0,J_1  ; do i=I_0,I_1
             CH4ICINL(:)=CH4ICIN(I,J,:)! now in PPPM
             CALL LOGPINT(LCOalt,PCOalt,CH4ICINL,LM,PRES,CH4ICL,.true.)
             CH4ICX(I,J,:)=CH4ICL(:)*scale_ch4_IC_file*am(:,i,j)*
     *            axyp(i,j)
            end do     ; end do
          end if
#ifdef DYNAMIC_BIOMASS_BURNING
          ! 12 below are the 12 VDATA veg types or Ent remapped to them,
          ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
          emisPerFireByVegType(n,1:12)=(/0.0000000d+00,1.0864168d-06,
     &    6.3624935d-07, 4.7021388d-07, 1.0293500d-06, 1.7132404d-06,
     &    1.4364367d-06, 3.0849296d-06, 0.0000000d+00, 0.0000000d+00,
     &    0.0000000d+00, 0.0000000d+00/)
#endif /* DYNAMIC_BIOMASS_BURNING */
#endif /* TRACERS_SPECIAL_Shindell */

      case ('O3')
      n_O3 = n
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
          ntsurfsrc(n) = 1
#ifdef TRACERS_SPECIAL_Lerner
C**** Get solar variability coefficient from namelist if it exits
          dsol = 0.
          call sync_param("dsol",dsol)
#endif

      case ('SF6_c')
      n_SF6_c = n
          ntm_power(n) = -14
          tr_mm(n) = 146.01d0
          ntsurfsrc(n) = 1
#endif /* TRACERS_ON */
      case ('Water')
      n_Water = n
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
          trw0(n) = 1.
          ntrocn(n)= 0
#endif
#ifdef TRACERS_WATER
          ntm_power(n) = -4
          tr_mm(n) = mwat
          needtrs(n) = .true.
          tr_wd_TYPE(n) = nWater
          trli0(n) = 1.
          trsi0(n) = 1.
          tr_H2ObyCH4(n) = 1.
#endif
#ifdef TRACERS_OCEAN
          trglac(n) = 1.
#endif
#ifdef TRACERS_SPECIAL_O18
          iso_index(n) = 1  ! indexing for isotopic fractionation calcs
#endif

#ifdef TRACERS_ON
#ifdef TRACERS_SPECIAL_O18
      case ('H2O18')
      n_H2O18 = n
          ntm_power(n) = -7
          tr_mm(n) = 20.
          needtrs(n) = .true.
          tr_wd_TYPE(n) = nWater
          iso_index(n) = 2  ! indexing for isotopic fractionation calcs
          trw0(n) = 2.228d-3   ! SMOW mass ratio of water molecules
          trli0(n) = 0.980d0*trw0(n)  ! d=-20
          trsi0(n) = fracls(n)*trw0(n)
          tr_H2ObyCH4(n) = trw0(n)*1.023d0 ! d=+23 (ie. no frac from O2)
          ntrocn(n) = -3
#ifdef TRACERS_OCEAN
          trglac(n) = trw0(n)*0.98d0   ! d=-20
#endif

      case ('HDO')
      n_HDO = n
          ntm_power(n) = -8
          tr_mm(n) = 19.
          needtrs(n) = .true.
          tr_wd_TYPE(n) = nWater
          iso_index(n) = 3  ! indexing for isotopic fractionation calcs
          trw0(n) = 3.29d-4    ! SMOW mass ratio of water molecules
          trli0(n) = 0.830d0*trw0(n)  ! d=-170
          trsi0(n) = fracls(n)*trw0(n)
          tr_H2ObyCH4(n) = trw0(n)*0.93d0  ! d=-70
          ntrocn(n) = -4
#ifdef TRACERS_OCEAN
          trglac(n) = trw0(n)*0.84d0   ! d=-160
#endif
C**** RDF
#ifdef TES_SIM
          call ReadTESAvgKrnlFile
#endif

      case ('HTO')
      n_HTO = n
          ntm_power(n) = -18
          tr_mm(n) = 20.
          needtrs(n) = .true.
          tr_wd_TYPE(n) = nWater
          iso_index(n) = 4  ! indexing for isotopic fractionation calcs
          trw0(n) = 0. !2.22d-18   ! SMOW mass ratio of water molecules
          trli0(n) = 0.
          trsi0(n) = 0.
          tr_H2ObyCH4(n) = 0.
          trdecay(n) = 1.77d-9      ! =5.59d-2 /yr
          ntrocn(n) = -18
#ifdef TRACERS_OCEAN
          trglac(n) = 0.
#endif

      case ('H2O17')
      n_H2O17 = n
          ntm_power(n) = -7
          tr_mm(n) = 19.
          needtrs(n) = .true.
          tr_wd_TYPE(n) = nWater
          iso_index(n) = 5  ! indexing for isotopic fractionation calcs
          trw0(n) = 4.020d-5   ! SMOW mass ratio of water molecules
          trli0(n) = 0.98937d0*trw0(n)  ! d=-10.63 D17O=0
          trsi0(n) = fracls(n)*trw0(n)
          tr_H2ObyCH4(n) = trw0(n)*1.011596d0 ! d=+11.596 (some frac from O2)
          ntrocn(n) = -3
#ifdef TRACERS_OCEAN
          trglac(n) = trw0(n)*0.98937d0   ! d=-10.63 D17O=
#endif
#endif  /* TRACERS_SPECIAL_O18 */

#ifdef TRACERS_SPECIAL_Shindell
      case ('Ox')
      n_Ox = n
          call openunit('Ox_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),OxICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           OxICINL(:)=OxICIN(I,J,:)! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,OxICINL,LM,PRES,OxICL,.true.)
           OxIC(I,J,:)=OxICL(:)*am(:,i,j)*axyp(i,j)
          end do     ; end do
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
#ifdef TRACERS_DRYDEP
          F0(n) = 1.4d0
          HSTAR(n) = 1.d-2
#endif
      case ('NOx')
      n_NOx = n
          ntm_power(n) = -11
          tr_mm(n) = 14.01d0
#ifdef TRACERS_DRYDEP
          F0(n) = 1.d-1
          HSTAR(n) = 1.d-2
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
          ! 12 below are the 12 VDATA veg types or Ent remapped to them,
          ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
          emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 1.1378230d-07,
     &    3.2166037d-07, 1.5559274d-07, 4.1611088d-07, 5.7316458d-07,
     &    2.1700112d-07, 3.0054335d-07, 0.0000000d+00, 0.0000000d+00,
     &    0.0000000d+00, 0.0000000d+00/)
#endif
          call check_aircraft_sectors ! special 3D source case

      case ('N2O5')
      n_N2O5 = n
          ntm_power(n) = -12
          tr_mm(n) = 108.02d0

#ifdef SHINDELL_STRAT_CHEM
      case ('ClOx')
      n_ClOx = n
          ntm_power(n) = -11
          tr_mm(n) = 51.5d0
C         Interpolate ClOx altitude-dependence to model resolution:
          CALL LOGPINT(LCOalt,PCOalt,ClOxaltIN,LM,PRES,ClOxalt,.true.)

      case ('BrOx')
      n_BrOx = n
          ntm_power(n) = -14
          tr_mm(n) = 95.9d0
C         Interpolate BrOx altitude-dependence to model resolution:
          CALL LOGPINT(LCOalt,PCOalt,BrOxaltIN,LM,PRES,BrOxalt,.true.)

      case ('HCl')
      n_HCl = n
          ntm_power(n) = -10
          tr_mm(n) = 36.5d0
C         Interpolate HCl altitude-dependence to model resolution:
          CALL LOGPINT(LCOalt,PCOalt,HClaltIN,LM,PRES,HClalt,.true.)

      case ('ClONO2')
      n_ClONO2 = n
          ntm_power(n) = -11
          tr_mm(n) = 97.5d0
C         Interpolate ClONO2 altitude-dependence to model resolution:
          CALL
     &    LOGPINT(LCOalt,PCOalt,ClONO2altIN,LM,PRES,ClONO2alt,.true.)
#endif

      case ('HOCl')
      n_HOCl = n
          ntm_power(n) = -12
          tr_mm(n) = 52.5d0

      case ('HBr')
      n_HBr = n
          ntm_power(n) = -14
          tr_mm(n) = 80.9d0

      case ('HOBr')
      n_HOBr = n
          ntm_power(n) = -14
          tr_mm(n) = 96.9d0

      case ('BrONO2')
      n_BrONO2 = n
          ntm_power(n) = -14
          tr_mm(n) = 141.9d0

      case ('CFC')
      n_CFC = n
          ntm_power(n) = -12
          tr_mm(n) = 137.4d0 !CFC11
#ifdef SHINDELL_STRAT_CHEM
          if(AM_I_ROOT( ))then
C          check on GHG file's 1995 value for CFCs:
           call openunit('GHG',iu_data,.false.,.true.)
           do i=1,5; read(iu_data,'(a80)') title; enddo
           temp_year=0
           do while(temp_year <= 1995)
             read(iu_data,*,end=101) temp_year,(temp_ghg(j),j=1,6)
             if(temp_year==1995)then
               temp_ghg(1)=cfc_rad95*0.95d0
               temp_ghg(2)=cfc_rad95*1.05d0
               temp_ghg(3)=(temp_ghg(4)+temp_ghg(5))*1.d-9
               if(temp_ghg(3) < temp_ghg(1) .or.
     &         temp_ghg(3) > temp_ghg(2))then
                 call stop_model('please check on cfc_rad95 2',255)
               endif
             endif
           enddo
 101       continue
           if(temp_year<1995)
     &     call stop_model('please check on cfc_rad95 1',255)
           call closeunit(iu_data)
          endif
C          read the CFC initial conditions:
          call openunit('CFC_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),CFCICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           CFCICINL(:)=CFCICIN(I,J,:)! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,CFCICINL,LM,PRES,CFCICL,.true.)
           CFCIC(I,J,:)=CFCICL(:)*am(:,i,j)*axyp(i,j)
          end do     ; end do
#endif

      case ('HNO3')
      n_HNO3 = n
          ntm_power(n) = -11
          tr_mm(n) = 63.018d0
          tr_RKD(n) = 2.073d3 ! in mole/J = 2.1d5 mole/(L atm)
#ifdef TRACERS_DRYDEP
          HSTAR(n)=1.d14
#endif

      case ('H2O2')
      n_H2O2 = n
          ntm_power(n) = -11
          tr_mm(n) = 34.016d0
          tr_RKD(n) = 9.869d2    ! in mole/J = 1.d5 mole/(L atm)
          tr_DHD(n) = -5.52288d4 ! in J/mole = -13.2 kcal/mole.
#ifdef TRACERS_DRYDEP
          HSTAR(n)=tr_RKD(n)*convert_HSTAR
          F0(n) = 1.d0
#endif

#ifdef SHINDELL_STRAT_EXTRA
      case ('GLT') ! generic linear tracer
      n_GLT = n
          ntm_power(n) = -11
          tr_mm(n) = mair

#ifdef ACCMIP_LIKE_DIAGS
      case ('stratOx')
      n_stratOx = n
          ! assumes initial Ox conditions read in for Ox tracer
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
#ifdef TRACERS_DRYDEP
          F0(n) = 1.0d0
          HSTAR(n) = 1.d-2
#endif
      case('codirect')
      n_codirect = n
          ntm_power(n) = -8
          tr_mm(n) = 28.01d0
          trdecay(n) = 2.31482d-7 ! 1/(50 days)
          ! not a radiactive decay, but functionally identical

#endif /* ACCMIP_LIKE_DIAGS */
#endif /* SHINDELL_STRAT_EXTRA */

      case ('CH3OOH')
      n_CH3OOH = n
          ntm_power(n) = -11
          tr_mm(n) = 48.042d0
#ifdef TRACERS_DRYDEP
          HSTAR(n) = 3.d2
#endif

      case ('HCHO')
      n_HCHO = n
          ntm_power(n) = -11
          tr_mm(n) = 30.026d0
          tr_RKD(n) = 6.218d1 ! mole/J = 6.3d3 mole/(L atm)
#ifdef TRACERS_DRYDEP
          HSTAR(n)=6.d3
#endif

      case ('HO2NO2')
      n_HO2NO2 = n
          ntm_power(n) = -12
          tr_mm(n) = 79.018d0

      case ('CO')
      n_CO = n
          call openunit('CO_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),COICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           COICINL(:)=COICIN(I,J,:)! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,COICINL,LM,PRES,COICL,.true.)
           COIC(I,J,:)=COICL(:)*am(:,i,j)*axyp(i,j)
          end do     ; end do
          ntm_power(n) = -8
          tr_mm(n) = 28.01d0
#ifdef DYNAMIC_BIOMASS_BURNING
          ! 12 below are the 12 VDATA veg types or Ent remapped to them,
          ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
          emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 7.0401156d-06,
     &    1.8708386d-05, 1.0678024d-05, 2.6742857d-05, 4.0226296d-05,
     &    2.3661527d-05, 4.4639346d-05, 0.0000000d+00, 0.0000000d+00,
     &    0.0000000d+00, 0.0000000d+00/)
#endif

      case ('PAN')
      n_PAN = n
          ntm_power(n) = -11
          tr_mm(n) = 121.054d0   ! assuming CH3COOONO2 = PAN
#ifdef TRACERS_DRYDEP
          HSTAR(n) = 3.6d0
#endif

      case ('Isoprene')
      n_Isoprene = n
          ntm_power(n) = -11
          tr_mm(n) = 60.05d0 ! i.e. 5 carbons
#ifdef TRACERS_DRYDEP
          HSTAR(n) = 1.3d-2
#endif

      case ('AlkylNit')
      n_AlkylNit = n
          ntm_power(n) = -11
          tr_mm(n) = mair !unknown molecular weight, so use air and make
                          ! note in the diagnostics write-out...

      case ('Alkenes')
      n_Alkenes = n
          ntm_power(n) = -10
          tr_mm(n) = 1.0d0 ! So, careful: source files now in Kmole/m2/s or
                           ! equivalently, kg/m2/s for species with tr_mm=1
#ifdef DYNAMIC_BIOMASS_BURNING
          ! 12 below are the 12 VDATA veg types or Ent remapped to them,
          ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
          emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 6.1516259d-09,
     &    1.1544214d-08, 6.9501711d-09, 1.7481154d-08, 2.5840087d-08,
     &    1.5709551d-08, 3.5913079d-08, 0.0000000d+00, 0.0000000d+00,
     &    0.0000000d+00, 0.0000000d+00/)
#endif

      case ('Paraffin')
      n_Paraffin = n
          ntm_power(n) = -10
          tr_mm(n) = 1.0d0 ! So, careful: source files now in Kmole/m2/s or
                           ! equivalently, kg/m2/s for species with tr_mm=1
#ifdef DYNAMIC_BIOMASS_BURNING
          ! 12 below are the 12 VDATA veg types or Ent remapped to them,
          ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
          emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 1.5258348d-09,
     &    5.6236904d-09, 3.1752858d-09, 1.0662656d-08, 1.5271524d-08,
     &    8.0735774d-09, 2.6055675d-08, 0.0000000d+00, 0.0000000d+00,
     &    0.0000000d+00, 0.0000000d+00/)
#endif

#endif  /* TRACERS_SPECIAL_Shindell */

#ifdef TRACERS_TERP
      case ('Terpenes')
      n_Terpenes = n
          ntm_power(n) = -11
          tr_mm(n) = 120.10d0 ! i.e. 10 carbons
#ifdef TRACERS_DRYDEP
          HSTAR(n) = 1.3d-2
#endif
#endif  /* TRACERS_TERP */

#ifdef TRACERS_AEROSOLS_SOA
      case ('isopp1g')
          n_isopp1g = n
          n_soa_i = n_isopp1g        !the first from the soa species
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          tr_RKD(n) = 1.d4 / convert_HSTAR !Henry; from mole/(L atm) to mole/J
          tr_DHD(n) = -12.d0 * gasc        !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
          tr_wd_TYPE(n) = nGAS
#ifdef TRACERS_DRYDEP
          HSTAR(n) = tr_RKD(n) * convert_HSTAR
#endif

      case ('isopp1a')
          n_isopp1a = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          trpdens(n) = 1.5d3         !kg/m3
          trradius(n) = 3.d-7        !m
          fq_aer(n) = 0.8d0          !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART

      case ('isopp2g')
          n_isopp2g = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          tr_RKD(n) = 1.d4 / convert_HSTAR !Henry; from mole/(L atm) to mole/J
          tr_DHD(n) = -12.d0 * gasc        !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
          tr_wd_TYPE(n) = nGAS
#ifdef TRACERS_DRYDEP
          HSTAR(n) = tr_RKD(n) * convert_HSTAR
#endif

      case ('isopp2a')
          n_isopp2a = n
#ifndef TRACERS_TERP
          n_soa_e = n_isopp2a        !the last from the soa species
#endif  /* TRACERS_TERP */
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          trpdens(n) = 1.5d3         !kg/m3
          trradius(n) = 3.d-7        !m
          fq_aer(n) = 0.8d0          !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART

#ifdef TRACERS_TERP
      case ('apinp1g')
          n_apinp1g = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          tr_RKD(n) = 1.d4 / convert_HSTAR !Henry; from mole/(L atm) to mole/J
          tr_DHD(n) = -12.d0 * gasc        !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
          tr_wd_TYPE(n) = nGAS
#ifdef TRACERS_DRYDEP
          HSTAR(n) = tr_RKD(n) * convert_HSTAR
#endif

      case ('apinp1a')
          n_apinp1a = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          trpdens(n) = 1.5d3         !kg/m3
          trradius(n) = 3.d-7        !m
          fq_aer(n) = 0.8d0          !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART

      case ('apinp2g')
          n_apinp2g = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          tr_RKD(n) = 1.d4 / convert_HSTAR !Henry; from mole/(L atm) to mole/J
          tr_DHD(n) = -12.d0 * gasc        !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
          tr_wd_TYPE(n) = nGAS
#ifdef TRACERS_DRYDEP
          HSTAR(n) = tr_RKD(n) * convert_HSTAR
#endif

      case ('apinp2a')
          n_apinp2a = n
          n_soa_e = n_apinp2a        !the last from the soa species
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          trpdens(n) = 1.5d3         !kg/m3
          trradius(n) = 3.d-7        !m
          fq_aer(n) = 0.8d0          !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */

      case ('DMS')
      n_DMS = n
          ntm_power(n) = -12
! the ocean source of DMS is actually interactive and therefore should
! not count for ntsurfsrc....
          ntsurfsrc(n) = 0   ! ocean DMS concentration
          ntisurfsrc(n)=1.
          tr_mm(n) = 62.
          needtrs(n)=.true.

      case ('MSA')
      n_MSA = n
          ntm_power(n) = -13
          tr_mm(n) = 96.  !(H2O2 34;SO2 64)
          trpdens(n)=1.7d3   !kg/m3 this is sulfate value
          trradius(n)=5.d-7 !m (SO4 3;BC 1;OC 3)
          fq_aer(n)=1.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART

      case ('SO2')
      n_SO2 = n
          ntm_power(n) = -11
          tr_mm(n) = 64.
          tr_RKD(n) =0.0118d0 !mole/J or  1.2  M/atm
          tr_DHD(n) =-2.62d4! in J/mole= -6.27 kcal/mol
          tr_wd_TYPE(n) = nGAS
#ifdef TRACERS_DRYDEP
c         HSTAR(n)=tr_RKD(n)*convert_HSTAR
          HSTAR(N)=1.D5
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
          ! 12 below are the 12 VDATA veg types or Ent remapped to them,
          ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
          emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 2.6305068d-07,
     &    1.3513656d-07, 1.0093965d-07, 1.6911058d-07, 3.2019645d-07,
     &    3.1232341d-07, 4.1607765d-07, 0.0000000d+00, 0.0000000d+00,
     &    0.0000000d+00, 0.0000000d+00/)
#endif

      case ('SO4')
      n_SO4 = n
          ntm_power(n) = -11
          tr_mm(n) = 96.
          trpdens(n)=1.7d3   !kg/m3 this is sulfate value
          trradius(n)=3.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART

#ifdef TRACERS_HETCHEM
      case ('SO4_d1')
      n_SO4_d1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.   !!!! Sulfat
          trpdens(n)=2.5d3   !kg/m3 this is clay density
          trradius(n)=0.75D-06 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('SO4_d2')
      n_SO4_d2 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)=2.65d3   !kg/m3 this is silt1 value
          trradius(n)=2.2D-06 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('SO4_d3')
      n_SO4_d3 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)=2.65d3   !this is silt2 value
          trradius(n)=4.4D-06 !m this is silt2 value
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('N_d1')
      n_N_d1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 62.   ! NO3
          trpdens(n)=2.5d3   !kg/m3 this is clay density
          trradius(n)=0.75D-06 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('N_d2')
      n_N_d2 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 62.
          trpdens(n)=2.65d3   !kg/m3 this is silt1 value
          trradius(n)=2.2D-06 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('N_d3')
      n_N_d3 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 62.
          trpdens(n)=2.65d3   !this is silt2 value
          trradius(n)=4.4D-06 !m this is silt2 value
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
#endif

      case ('BCII')  !Insoluble industrial BC
      n_BCII = n
          ntm_power(n) = -12
          ntsurfsrc(n) = 2  !Industrial, ships
          tr_mm(n) = 12.
          trpdens(n)=1.3d3   !kg/m3
          trradius(n)=1.d-7 !m
          fq_aer(n)=0.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('BCIA') !Soluble (aged) industrial BC
      n_BCIA = n
          ntm_power(n) = -12
          ntsurfsrc(n) = 0
          tr_mm(n) = 12.
          trpdens(n)=1.3d3   !kg/m3
          trradius(n)=1.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('BCB') !Biomass BC
      n_BCB = n
          ntm_power(n) = -12
          tr_mm(n) = 12.
          trpdens(n)=1.3d3   !kg/m3
          trradius(n)=1.d-7 !m
          fq_aer(n)=0.6   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
#ifdef DYNAMIC_BIOMASS_BURNING
          ! 12 below are the 12 VDATA veg types or Ent remapped to them,
          ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
          emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 8.2731990d-09,
     &    1.7817767d-07, 8.5456378d-08, 2.5662467d-07, 3.3909114d-07,
     &    1.1377826d-07, 2.9145593d-07, 0.0000000d+00, 0.0000000d+00,
     &    0.0000000d+00, 0.0000000d+00/)
#endif
      case ('OCII') !Insoluble industrial organic mass
      n_OCII = n
          ntm_power(n) = -11
#ifdef TRACERS_AEROSOLS_SOA
          ntsurfsrc(n) = 2 !industrial, ships
#else
          ntsurfsrc(n) = 3 !Terpene, industrial, ships
#endif  /* TRACERS_AEROSOLS_SOA */
          tr_mm(n) = 15.6d0
          trpdens(n)=1.5d3   !kg/m3
          trradius(n)=3.d-7 !m
          fq_aer(n)=0.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCIA') !Aged industrial organic mass
      n_OCIA = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          trpdens(n)=1.5d3   !kg/m3
          trradius(n)=3.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCI1') !Insoluble organic mass, 1st type: WIOC-SA
      n_OCI1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0 !3D emissions
          tr_mm(n) = 15.6d0
          trpdens(n)=1.5d3   !kg/m3
          trradius(n)=3.d-7 !m
          fq_aer(n)=0.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCA1') !Aged organic mass, 1st type
      n_OCA1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          trpdens(n)=1.5d3   !kg/m3
          trradius(n)=3.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCI2') !Insoluble organic mass, 2nd type: WIOC-NA
      n_OCI2 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0 !industrial emissions
          tr_mm(n) = 15.6d0
          trpdens(n)=1.5d3   !kg/m3
          trradius(n)=3.d-7 !m
          fq_aer(n)=0.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCA2') !Aged organic mass, 2nd type
      n_OCA2 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          trpdens(n)=1.5d3   !kg/m3
          trradius(n)=3.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCI3') !Insoluble organic mass, 3rd type: WSOC-LS
      n_OCI3 = n
          ntm_power(n) = -11
#ifdef TRACERS_AEROSOLS_SOA
          ntsurfsrc(n) = 0
#else
          ntsurfsrc(n) = 1 ! terpene emissions
#endif  /* TRACERS_AEROSOLS_SOA */
          tr_mm(n) = 15.6d0
          trpdens(n)=1.5d3   !kg/m3
          trradius(n)=3.d-7 !m
          fq_aer(n)=0.2   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCA3') !Aged organic mass, 3rd type
      n_OCA3 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          trpdens(n)=1.5d3   !kg/m3
          trradius(n)=3.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCA4') !soluble organic mass, 4th type: WSOC-MS
      n_OCA4 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6d0
          trpdens(n)=1.5d3   !kg/m3
          trradius(n)=3.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCB') !Biomass organic mass
      n_OCB = n
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
          call sync_param("OCB_om2oc",om2oc(n_OCB))
#endif
          ntm_power(n) = -11
          tr_mm(n) = 15.6d0
          trpdens(n)=1.5d3   !kg/m3
          trradius(n)=3.d-7 !m
          fq_aer(n)=0.8   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
#ifdef DYNAMIC_BIOMASS_BURNING
          ! 12 below are the 12 VDATA veg types or Ent remapped to them,
          ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
          emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 6.4230818d-07,
     &    1.6844633d-06, 8.7537586d-07, 2.5902559d-06, 4.3200689d-06,
     &    2.6824284d-06, 3.9549395d-06, 0.0000000d+00, 0.0000000d+00,
     &    0.0000000d+00, 0.0000000d+00/)
#endif
      case ('Be7')
      n_Be7 = n
CCC#ifdef SHINDELL_STRAT_EXTRA
CCC          ntm_power(n) = -21        ! power of ten for tracer
CCC#else
          ntm_power(n) = -23        ! power of ten for tracer
CCC#endif
          tr_mm(n) = 7.d0
          trdecay(n) = 1.51d-7
          trpdens(n) = 1.7d3    !kg/m3 this is SO4 value
          trradius(n) = 1.d-7  !appropriate for stratosphere
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART ! same as SO4

      case ('Be10')
      n_Be10 = n
CCC#ifdef SHINDELL_STRAT_EXTRA
CCC          ntm_power(n) = -21
CCC#else
          ntm_power(n) = -23
CCC#endif
          tr_mm(n) = 10.d0
          trpdens(n) = 1.7d3   !kg/m3 this is SO4 value
          trradius(n) = 1.d-7  !appropriate for stratosphere
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART ! same as SO4

      case ('Pb210')
          n_Pb210 = n
          ntm_power(n) = -23
          tr_mm(n) = 210.d0
          trdecay(n) = 9.85d-10
          trpdens(n) = 1.7d3    !kg/m3 this is SO4 value
          trradius(n) = 3.d-7  !again S04 value
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART ! same as SO4

      case ('H2O2_s')
      n_H2O2_s = n
          ntm_power(n) = -10
          tr_mm(n) = 34.016
          tr_RKD(n) = 986.9
          tr_DHD(n) = -5.52288d4 ! in J/mole = -13.2 kcal/mole.
          tr_wd_TYPE(n) = nGAS
#ifdef TRACERS_DRYDEP
          HSTAR(n)=tr_RKD(n)*convert_HSTAR
          F0(n) = 1.d0
#endif
      case ('seasalt1')
      n_seasalt1 = n
          ntsurfsrc(n) = 0   ! ocean bubbles
          ntisurfsrc(n)=1.
          ntm_power(n) = -10
          tr_mm(n) = 75. !Na x 3.256
          trpdens(n)=2.2d3  !kg/m3 This is for non-hydrated
          trradius(n)=4.4d-7 ! This is non-hydrated
          fq_aer(n)=1.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('seasalt2')
      n_seasalt2 = n
          ntsurfsrc(n) = 0   ! ocean bubbles
          ntisurfsrc(n)=1.
          ntm_power(n) = -9
          tr_mm(n) = 75. !Na x 3.256
          trpdens(n)=2.2d3  !kg/m3 This is for non-hydrated
          trradius(n)=5.0d-6                 ! This is non-hydrated
#ifdef TRACERS_AEROSOLS_Koch
          if (imAER.ne.1) trradius(n)=1.7d-6 ! This is non-hydrated
#endif
          fq_aer(n)=1.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('NH3')
      n_NH3 = n
          ntm_power(n) = -10
          tr_mm(n) = 17.
          tr_RKD(n) = 0.7303   !tr_RKD=74 M/atm
          tr_DHD(n) = -2.84d4  !tr_DHD=-6.80 kcal/mole
          tr_wd_TYPE(n) = nGAS
#ifdef TRACERS_DRYDEP
          HSTAR(n)=tr_RKD(n)*convert_HSTAR
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
          ! 12 below are the 12 VDATA veg types or Ent remapped to them,
          ! from Olga Pechony's AR5_EPFC_factors_corrected_NH3.xlsx file.
          emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 1.2274993d-06,
     &    3.8813269d-07, 3.5230462d-07, 4.0781484d-07, 8.8901584d-07,
     &    1.1341459d-06, 1.4117913d-06, 0.0000000d+00, 0.0000000d+00,
     &    0.0000000d+00, 0.0000000d+00/)
#endif
      case ('NH4')
      n_NH4 = n
          ntsurfsrc(n) = 0
          ntm_power(n) = -10
          tr_mm(n) = 18.
          trpdens(n)=1.7d3
          trradius(n)=3.d-7
          fq_aer(n)=1.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('NO3p')
      n_NO3p = n
          ntsurfsrc(n) = 0
          ntm_power(n) = -11
          tr_mm(n) = 62.
          trpdens(n)=1.7d3
          trradius(n)=3.d-7
          fq_aer(n)=1.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART

#ifdef TRACERS_DUST
      CASE('Clay')
      n_clay=n
          ntm_power(n)=-9
          trpdens(n)=2.5d3
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Silt1')
      n_silt1=n
          ntm_power(n)=-9
          trpdens(n)=2.65d3
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Silt2')
      n_silt2=n
          ntm_power(n)=-9
          trpdens(n)=2.65d3
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Silt3')
      n_silt3=n
          ntm_power(n)=-9
          trpdens(n)=2.65d3
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Silt4')
      n_silt4=n
          ntm_power(n)=-9
          trpdens(n)=2.65d3
#ifdef TRACERS_DRYDEP
          trradius(n)=11.77D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

#else
#ifdef TRACERS_MINERALS
      CASE('ClayIlli')          ! http://webmineral.com/data/Illite.shtml
      n_clayilli=n
          ntm_power(n)=-9
          trpdens(n) = 2.795d3  ! measured; http://www.mindat.org/min-2011.html
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('ClayKaol')       ! http://www.webmineral.com/data/Kaolinite.shtml
      n_claykaol=n
          ntm_power(n)=-9
          trpdens(n) = 2.63D3  ! calculated; http://www.mindat.org/min-2011.html
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('ClaySmec')       ! http://www.webmineral.com/data/Rectorite.shtml
      n_claysmec=n
          ntm_power(n)=-9
          trpdens(n)=2.35D3     ! for Montmorillonite
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('ClayCalc')       ! http://www.webmineral.com/data/Calcite.shtml
      n_claycalc=n
          ntm_power(n)=-9
          trpdens(n) = 2.71d3   ! measured; http://www.mindat.org/min-859.html
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('ClayQuar')       ! http://www.webmineral.com/data/Quartz.shtml
      n_clayquar=n
          ntm_power(n)=-9
          trpdens(n) = DensityQuartz
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil1Quar')       ! http://www.webmineral.com/data/Quartz.shtml
      n_sil1quar=n
          ntm_power(n)=-9
          trpdens(n) = DensityQuartz
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil1Feld')       ! http://webmineral.com/data/Plagioclase.shtml
      n_sil1feld=n
          ntm_power(n)=-9
          trpdens(n) = 2.68d3     ! average Plagioclase
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil1Calc')       ! http://www.webmineral.com/data/Calcite.shtml
      n_sil1calc=n
          ntm_power(n)=-9
          trpdens(n) = 2.71d3   ! measured; http://www.mindat.org/min-859.html
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil1Hema')       ! http://www.webmineral.com/data/Hematite.shtml
      n_sil1hema=n
          ntm_power(n)=-9
          trpdens(n) = DensityHematite
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil1Gyps')       ! http://www.webmineral.com/data/Gypsum.shtml
      n_sil1gyps=n
          ntm_power(n)=-9
          trpdens(n) = 2.312d3  ! measured; http://www.mindat.org/min-1784.html
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil2Quar')       ! http://www.webmineral.com/data/Quartz.shtml
      n_sil2quar=n
          ntm_power(n)=-9
          trpdens(n) = DensityQuartz
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil2Feld')       ! http://webmineral.com/data/Plagioclase.shtml
      n_sil2feld=n
          ntm_power(n)=-9
          trpdens(n) = 2.68d3     ! average Plagioclase
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil2Calc')       ! http://www.webmineral.com/data/Calcite.shtml
      n_sil2calc=n
          ntm_power(n)=-9
          trpdens(n) = 2.71d3   ! measured; http://www.mindat.org/min-859.html
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil2Hema')       ! http://www.webmineral.com/data/Hematite.shtml
      n_sil2hema=n
          ntm_power(n)=-9
          trpdens(n) = DensityHematite
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil2Gyps')       ! http://www.webmineral.com/data/Gypsum.shtml
      n_sil2gyps=n
          ntm_power(n)=-9
          trpdens(n) = 2.312d3  ! measured; http://www.mindat.org/min-1784.html
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil3Quar')       ! http://www.webmineral.com/data/Quartz.shtml
      n_sil3quar=n
          ntm_power(n)=-9
          trpdens(n) = DensityQuartz
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil3Feld')       ! http://webmineral.com/data/Plagioclase.shtml
      n_sil3feld=n
          ntm_power(n)=-9
          trpdens(n) = 2.68d3     ! average Plagioclase
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil3Calc')       ! http://www.webmineral.com/data/Calcite.shtml
      n_sil3calc=n
          ntm_power(n)=-9
          trpdens(n) = 2.71d3   ! measured; http://www.mindat.org/min-859.html
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil3Hema')       ! http://www.webmineral.com/data/Hematite.shtml
      n_sil3hema=n
          ntm_power(n)=-9
          trpdens(n) = DensityHematite
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil3Gyps')       ! http://www.webmineral.com/data/Gypsum.shtml
      n_sil3gyps=n
          ntm_power(n)=-9
          trpdens(n) = 2.312d3  ! measured; http://www.mindat.org/min-1784.html
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

#endif  /* TRACERS_MINERALS */
#ifdef TRACERS_QUARZHEM
      CASE('Sil1QuHe')
      n_sil1quhe=n
          ntm_power(n)=-9
          trpdens(n) = (1.d0 - frHemaInQuarAggr) * DensityQuartz +
     &         frHemaInQuarAggr * DensityHematite
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil2QuHe')
      n_sil2quhe=n
          ntm_power(n)=-9
          trpdens(n) = (1.d0 - frHemaInQuarAggr) * DensityQuartz +
     &         frHemaInQuarAggr * DensityHematite
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

      CASE('Sil3QuHe')
      n_sil3quhe=n
          ntm_power(n)=-9
          trpdens(n) = (1.d0 - frHemaInQuarAggr) * DensityQuartz +
     &         frHemaInQuarAggr * DensityHematite
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=5.D-1
          rc_washt(n)=5.D-1
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
          isdust(n) = 1

#endif  /* TRACERS_QUARZHEM */
#endif  /* TRACERS_DUST */
#ifdef TRACERS_AMP
C**** Tracers for Scheme AMP: Aerosol Microphysics (Mechanism M1 - M8)
      case ('H2SO4')
      n_H2SO4 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 98.
          trpdens(n)= DENS_SULF
          trradius(n)= DG_ACC * .5d-6
          fq_aer(n) = SOLU_ACC
          tr_wd_TYPE(n) = nPART
      case ('M_NO3')
      n_M_NO3 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 62.
          trpdens(n)=1.7d3
          trradius(n)=3.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('M_NH4')
      n_M_NH4 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 18.
          trpdens(n)=1.7d3
          trradius(n)=3.d-7
          fq_aer(n)=1.
          tr_wd_TYPE(n) = nPART
      case ('M_H2O')
      n_M_H2O = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = mwat
          trpdens(n)=1.d3
          trradius(n)=3.d-7
          fq_aer(n)=1.
          tr_wd_TYPE(n) = nPART !nWater
      case ('M_AKK_SU')
      n_M_AKK_SU = n
          ntm_power(n) = -11
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_AKK * .5d-6
          fq_aer(n)= SOLU_AKK
          tr_wd_TYPE(n) = nPART
       case ('N_AKK_1')
      n_N_AKK_1 = n
          ntm_power(n) = -11
          tr_mm(n) = 1.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_AKK * .5d-6
          fq_aer(n) = SOLU_AKK
          tr_wd_TYPE(n) = nPART
       case ('M_ACC_SU')
      n_M_ACC_SU = n
          ntm_power(n) = -11
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)= DG_ACC * .5d-6
          fq_aer(n) = SOLU_ACC
          tr_wd_TYPE(n) = nPART
       case ('N_ACC_1')
      n_N_ACC_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_ACC * .5d-6
          fq_aer(n) = SOLU_ACC
          tr_wd_TYPE(n) = nPART
       case ('M_DD1_SU')
      n_M_DD1_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_DD1 * .5d-6
          fq_aer(n) = SOLU_DD1
          tr_wd_TYPE(n) = nPART
       case ('M_DD1_DU')
      n_M_DD1_DU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_DUST
          trradius(n)=DG_DD1 * .5d-6
          fq_aer(n)= SOLU_DD1
          tr_wd_TYPE(n) = nPART
       case ('N_DD1_1')
      n_N_DD1_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_DUST
          trradius(n)=DG_DD1 * .5d-6
          fq_aer(n)= SOLU_DD1
          tr_wd_TYPE(n) = nPART
       case ('M_DS1_SU')
      n_M_DS1_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_DS1 * .5d-6
          fq_aer(n) = SOLU_DS1
          tr_wd_TYPE(n) = nPART
       case ('M_DS1_DU')
      n_M_DS1_DU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_DUST
          trradius(n)=DG_DS1 * .5d-6
          fq_aer(n) = SOLU_DS1
          tr_wd_TYPE(n) = nPART
       case ('N_DS1_1')
      n_ N_DS1_1= n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_DUST
          trradius(n)=DG_DS1 * .5d-6
          fq_aer(n) = SOLU_DS1
          tr_wd_TYPE(n) = nPART
       case ('M_DD2_SU')
      n_M_DD2_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_DD2 * .5d-6
          fq_aer(n) = SOLU_DD2
          tr_wd_TYPE(n) = nPART
       case ('M_DD2_DU')
      n_M_DD2_DU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_DUST
          trradius(n)=DG_DD2 * .5d-6
          fq_aer(n)= SOLU_DD2
          tr_wd_TYPE(n) = nPART
       case ('N_DD2_1')
      n_N_DD2_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_DUST
          trradius(n)=DG_DD2 * .5d-6
          fq_aer(n)= SOLU_DD2
          tr_wd_TYPE(n) = nPART
        case ('M_DS2_SU')
      n_M_DS2_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_DS2 * .5d-6
          fq_aer(n) = SOLU_DS2
          tr_wd_TYPE(n) = nPART
      case ('M_DS2_DU')
      n_M_DS2_DU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_DUST
          trradius(n)=DG_DS2 * .5d-6
          fq_aer(n) = SOLU_DS2
          tr_wd_TYPE(n) = nPART
      case ('N_DS2_1')
      n_N_DS2_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_DUST
          trradius(n)=DG_DS2 * .5d-6
          fq_aer(n) = SOLU_DS2
          tr_wd_TYPE(n) = nPART
      case ('M_SSA_SU')
      n_M_SSA_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_SSA * .5d-6
          fq_aer(n) = SOLU_SSA
          tr_wd_TYPE(n) = nPART
      case ('M_SSA_SS')
      n_M_SSA_SS = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 75.
          trpdens(n)= DENS_SEAS
          trradius(n)=DG_SSA * .5d-6
          fq_aer(n) = SOLU_SSA
          tr_wd_TYPE(n) = nPART
      case ('M_SSC_SS')
      n_M_SSC_SS = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 75.
          trpdens(n)= DENS_SEAS
          trradius(n)=DG_SSC * .5d-6
          fq_aer(n) = SOLU_SSC
          tr_wd_TYPE(n) = nPART
      case ('M_SSS_SS')
      n_M_SSS_SS = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 75.
          trpdens(n)= DENS_SEAS
          trradius(n)=DG_SSS * .5d-6
          fq_aer(n) = SOLU_SSS
          tr_wd_TYPE(n) = nPART
      case ('M_SSS_SU')
      n_M_SSS_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 75.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_SSS * .5d-6
          fq_aer(n) = SOLU_SSS
          tr_wd_TYPE(n) = nPART
      case ('M_OCC_SU')
      n_M_OCC_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_OCC * .5d-6
          fq_aer(n) = SOLU_OCC
          tr_wd_TYPE(n) = nPART
      case ('M_OCC_OC')
      n_M_OCC_OC = n
          ntm_power(n) = -11
          tr_mm(n) = 15.6
          trpdens(n)= DENS_OCAR
          trradius(n)=DG_OCC * .5d-6
          fq_aer(n)= SOLU_OCC
          tr_wd_TYPE(n) = nPART
      case ('N_OCC_1')
      n_N_OCC_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_OCAR
          trradius(n)=DG_OCC * .5d-6
          fq_aer(n)= SOLU_OCC
          tr_wd_TYPE(n) = nPART
      case ('M_BC1_SU')
      n_M_BC1_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_BC1 * .5d-6
          fq_aer(n) = SOLU_BC1
          tr_wd_TYPE(n) = nPART
      case ('M_BC1_BC')
      n_M_BC1_BC = n
          ntm_power(n) = -11
          tr_mm(n) = 12.
          trpdens(n)= DENS_BCAR
          trradius(n)=DG_BC1 * .5d-6
          fq_aer(n)= SOLU_BC1
          tr_wd_TYPE(n) = nPART
      case ('N_BC1_1')
      n_N_BC1_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_BCAR
          trradius(n)=DG_BC1 * .5d-6
          fq_aer(n)= SOLU_BC1
          tr_wd_TYPE(n) = nPART
      case ('M_BC2_SU')
      n_M_BC2_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_BC2 * .5d-6
          fq_aer(n) = SOLU_BC2
          tr_wd_TYPE(n) = nPART
      case ('M_BC2_BC')
      n_M_BC2_BC = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 12.
          trpdens(n)= DENS_BCAR
          trradius(n)=DG_BC2 * .5d-6
          fq_aer(n) = SOLU_BC2
          tr_wd_TYPE(n) = nPART
      case ('N_BC2_1')
      n_N_BC2_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_BCAR
          trradius(n)=DG_BC2 * .5d-6
          fq_aer(n) = SOLU_BC2
          tr_wd_TYPE(n) = nPART
      case ('M_BC3_SU')
      n_M_BC3_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_BC3 * .5d-6
          fq_aer(n) = SOLU_BC3
          tr_wd_TYPE(n) = nPART
      case ('M_BC3_BC')
      n_M_BC3_BC = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 12.
          trpdens(n)= DENS_BCAR
          trradius(n)=DG_BC3 * .5d-6
          fq_aer(n) = SOLU_BC3
          tr_wd_TYPE(n) = nPART
      case ('N_BC3_1')
      n_N_BC3_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_BCAR
          trradius(n)= DG_BC3 * .5d-6
          fq_aer(n) = SOLU_BC3
          tr_wd_TYPE(n) = nPART
      case ('M_DBC_SU')
      n_M_DBC_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_DBC * .5d-6
          fq_aer(n) = SOLU_DBC
          tr_wd_TYPE(n) = nPART
      case ('M_DBC_BC')
      n_M_DBC_BC = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_BCAR
          trradius(n)=DG_DBC * .5d-6
          fq_aer(n)= SOLU_DBC
          tr_wd_TYPE(n) = nPART
      case ('M_DBC_DU')
      n_M_DBC_DU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_DUST
          trradius(n)=DG_DBC * .5d-6
          fq_aer(n)= SOLU_DBC
          tr_wd_TYPE(n) = nPART
      case ('N_DBC_1')
      n_N_DBC_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_DUST
          trradius(n)=DG_DBC * .5d-6
          fq_aer(n)= SOLU_DBC
          tr_wd_TYPE(n) = nPART
      case ('M_BOC_SU')
      n_M_BOC_SU = n
          ntm_power(n) = -11
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_BOC * .5d-6
          fq_aer(n) = SOLU_BOC
          tr_wd_TYPE(n) = nPART
      case ('M_BOC_BC')
      n_M_BOC_BC = n
          ntm_power(n) = -11
          tr_mm(n) = 12.
          trpdens(n)= DENS_BCAR
          trradius(n)=DG_BOC * .5d-6
          fq_aer(n)= SOLU_BOC
          tr_wd_TYPE(n) = nPART
      case ('M_BOC_OC')
      n_M_BOC_OC = n
          ntm_power(n) = -11
          tr_mm(n) = 15.6
          trpdens(n)= DENS_OCAR
          trradius(n)=DG_BOC * .5d-6
          fq_aer(n)= SOLU_BOC
          tr_wd_TYPE(n) = nPART
      case ('N_BOC_1')
      n_N_BOC_1 = n
          ntm_power(n) = -11
          tr_mm(n) = 1.
          trpdens(n)= DENS_BCAR
          trradius(n)=DG_BOC * .5d-6
          fq_aer(n)= SOLU_BOC
          tr_wd_TYPE(n) = nPART
      case ('M_BCS_SU')
      n_M_BCS_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_BCS * .5d-6
          fq_aer(n) = SOLU_BCS
          tr_wd_TYPE(n) = nPART
      case ('M_BCS_BC')
      n_M_BCS_BC = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 12.
          trpdens(n)= DENS_BCAR
          trradius(n)=DG_BCS * .5d-6
          fq_aer(n) = SOLU_BCS
          tr_wd_TYPE(n) = nPART
      case ('N_BCS_1')
      n_N_BCS_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_BCAR
          trradius(n)=DG_BCS * .5d-6
          fq_aer(n) = SOLU_BCS
          tr_wd_TYPE(n) = nPART
      case ('M_MXX_SU')
      n_M_MXX_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_MXX * .5d-6
          fq_aer(n) = SOLU_MXX
          tr_wd_TYPE(n) = nPART
      case ('M_MXX_BC')
      n_M_MXX_BC = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 12.
          trpdens(n)= DENS_BCAR
          trradius(n)=DG_MXX * .5d-6
          fq_aer(n) = SOLU_MXX
          tr_wd_TYPE(n) = nPART
      case ('M_MXX_OC')
      n_M_MXX_OC = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6
          trpdens(n)= DENS_OCAR
          trradius(n)=DG_MXX * .5d-6
          fq_aer(n) = SOLU_MXX
          tr_wd_TYPE(n) = nPART
      case ('M_MXX_DU')
      n_M_MXX_DU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_DUST
          trradius(n)=DG_MXX * .5d-6
          fq_aer(n) = SOLU_MXX
          tr_wd_TYPE(n) = nPART
      case ('M_MXX_SS')
      n_M_MXX_SS = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 75.
          trpdens(n)= DENS_SEAS
          trradius(n)=DG_MXX * .5d-6
          fq_aer(n) = SOLU_MXX
          tr_wd_TYPE(n) = nPART
      case ('N_MXX_1')
      n_N_MXX_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_SEAS
          trradius(n)=DG_MXX * .5d-6
          fq_aer(n) = SOLU_MXX
          tr_wd_TYPE(n) = nPART
      case ('M_OCS_SU')
      n_M_OCS_SU = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)= DENS_SULF
          trradius(n)=DG_OCS * .5d-6
          fq_aer(n) = SOLU_OCS
          tr_wd_TYPE(n) = nPART
      case ('M_OCS_OC')
      n_M_OCS_OC = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 15.6
          trpdens(n)= DENS_OCAR
          trradius(n)=DG_OCS * .5d-6
          fq_aer(n) = SOLU_OCS
          tr_wd_TYPE(n) = nPART
      case ('N_OCS_1')
      n_N_OCS_1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 1.
          trpdens(n)= DENS_OCAR
          trradius(n)=DG_OCS * .5d-6
          fq_aer(n) = SOLU_OCS
          tr_wd_TYPE(n) = nPART
#endif /* TRACERS_AMP */
#endif /* TRACERS_ON */
      end select

#ifdef TRACERS_WATER
C**** Tracers that are soluble or are scavenged or are water => wet dep
      if (tr_wd_TYPE(n).eq.nWater.or.tr_wd_TYPE(n) .EQ. nPART .or.
     *  tr_RKD(n).gt.0) then
        dowetdep(n)=.true.
      end if
#endif
#ifdef TRACERS_DRYDEP
C**** If tracers are particles or have non-zero HSTAR or F0 do dry dep:
C**** Any tracers that dry deposits needs the surface concentration:
      if(HSTAR(n).GT.0..OR.F0(n).GT.0..OR.tr_wd_TYPE(n).eq.nPART) then
        dodrydep(n)=.true.
        needtrs(n)=.true.
#ifdef TRACERS_WATER
        if (tr_wd_TYPE(n).eq.nWATER) call stop_model
     &       ('A water tracer should not undergo dry deposition.',255)
#endif
      end if
#endif /* TRACERS_DRYDEP */

#ifdef TRACERS_ON
C**** Define the conversion from mass to volume units here
      mass2vol(n) = mair/tr_mm(n)
      vol2mass(n) = tr_mm(n)/mair
      to_conc(n)  = 0
#ifdef TRACERS_SPECIAL_Shindell
C**** Aerosol tracer output should be mass mixing ratio
      select case (tr_wd_TYPE(n))
        case (nGAS)
          to_volume_MixRat(n) = 1    !gas output to volume mixing ratio
        case (nPART)
          to_volume_MixRat(n) = 0    !aerosol output to mass mixing ratio
        case (nWATER)
          to_volume_MixRat(n) = 0    !water output to mass mixing ratio
        case default
          to_volume_MixRat(n) = 0    !default output to mass mixing ratio
      end select
#endif
#if defined(TRACERS_GASEXCH_ocean_CO2) || defined(TRACERS_GASEXCH_land_CO2)
          to_volume_MixRat(n) = 1    !gas output to volume mixing ratio
#endif
#endif /* TRACERS_ON */

      end do
#ifdef TRACERS_ON
#ifdef TRACERS_AEROSOLS_SOA
      call soa_init
#endif  /* TRACERS_AEROSOLS_SOA */
C**** Get to_volume_MixRat from rundecks if it exists
      call sync_param("to_volume_MixRat",to_volume_MixRat,ntm)
C**** Get to_conc from rundecks if it exists
      call sync_param("to_conc",to_conc,ntm)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
c**** soil dust aerosol initializations
      call init_soildust
#endif

C**** DIAGNOSTIC DEFINTIONS

C**** Set some diags that are the same regardless
      call set_generic_tracer_diags

C**** Zonal mean/height diags
      call init_jls_diag

C**** lat/lon tracer sources, sinks and specials
      call init_ijts_diag

C**** lat/lon/height tracer specials
      call init_ijlts_diag

C**** Initialize conservation diagnostics
      call init_tracer_cons_diag

C**** Miscellaneous initialisations

#ifdef TRACERS_DRYDEP
C Read landuse parameters and coefficients for tracer dry deposition:
      CALL RDLAND
      CALL RDDRYCF
#endif
#ifdef BIOGENIC_EMISSIONS
      CALL RDISOPCF
      CALL RDISOBASE
#endif
#ifdef TRACERS_SPECIAL_Shindell
      call cheminit ! **** Initialize the chemistry ****
      call special_layers_init
#endif
#ifdef TRACERS_COSMO
      do n=1,ntm
        if (trname(n) .eq. "Be7" .OR. trname(n) .eq. "Be10") then
          call init_cosmo
          exit
        end if
      end do
#endif
#endif /* TRACERS_ON */

#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_GASEXCH_ocean_CFC)
      !read in OCMIP based CFC-11 global emissions
      !=sum(dC/dt) for each hemisphere
      !these are *annual global averages* and need to be
      !converted to our timestep value
      print*, 'opening file=OCMIP_cfc.dat'
      call openunit('OCMIP_cfc',iu_data,.false.,.true.)
      do n=1,ntm
        do i=1,67
          read(iu_data,'(5x,e12.4)')ocmip_cfc(i,n)
        enddo
      enddo
      call closeunit(iu_data)
#endif

#if defined(TRACERS_GASEXCH_ocean_CO2) || defined(TRACERS_GASEXCH_land_CO2)
      call sync_param("atmCO2",atmCO2)
#endif

      return
      end subroutine init_tracer

      subroutine init_tracer_cons_diag
!@sum init_tracer_cons_diag Initialize tracer conservation diagnostics
!@auth Gavin Schmidt
      USE TRACER_COM
#ifdef TRACERS_ON
      USE TRDIAG_COM
#endif /* TRACERS_ON */
      implicit none
      character*20 sum_unit(ntm),inst_unit(ntm)   ! for conservation
      character*50 :: unit_string
#ifdef TRACERS_ON
      logical :: qcon(KTCON-1), qsum(KTCON-1), T=.TRUE. , F=.FALSE.
      logical :: Qf
      integer n,k,g,kk

C**** To add a new conservation diagnostic:
C****       Set up a QCON, and call SET_TCON to allocate array numbers,
C****       set up scales, titles, etc.
C**** QCON denotes when the conservation diags should be accumulated
C**** QSUM says whether that diag is to be used in summation (if the
C****      routine DIAGTCB is used, this must be false).
C**** 1:NPTS+1 ==> INST,  DYN,   COND,   RAD,   PREC,   LAND,  SURF,
C****            FILTER,STRDG/OCEAN, DAILY, OCEAN1, OCEAN2,
C**** First 12 are standard for all tracers and GCM
C**** Later indices are configurable - you provide title and itcon
C**** index (which is used wherever you want to check point)
C**** For example, separate Moist convection/Large scale condensation
c      itcon_mc(n)=xx
c      qcon(itcon_mc(n))=.true.  ; conpts(yy) = 'MOIST CONV'
c      qsum(itcon_mc(n)) = .false.
c      itcon_ss(n)=xx
c      qcon(itcon_ss(n))=.true.  ; conpts(yy) = 'LS COND'
c      qsum(itcon_ss(n)) = .false.

#ifdef CUBED_SPHERE
      Qf = .false.  ! no SLP filter
#else
      Qf = .true.   ! SLP filter on
#endif

      QCON=(/ t,                                           !instant.
     *        T,  T,  F,  F,  T,  T, Qf,  T,  F,  F,  F,   !2-12 (npts)
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F,       !13-22
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F        !21-ktcon-1
#ifdef TRACERS_AMP
     *      , F,  F,  F
#endif
     *      /)
      QSUM=(/ f,                                           !instant.
     *        T,  T,  F,  F,  T,  T, Qf,  T,  F,  F,  F,   !2-12 (npts)
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F,       !13-22
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F        !21-ktcon-1
#ifdef TRACERS_AMP
     *      , F,  F,  F
#endif
     *      /)
      do n=1,ntm
        kt_power_inst(n)   = ntm_power(n)+2
        kt_power_change(n) = ntm_power(n)-4
      end do

C**** set some defaults
      itcon_mc(:)=0
      itcon_AMP(:,:)=0
      itcon_AMPe(:)=0
      itcon_AMPm(:,:)=0
      itcon_ss(:)=0
      itcon_surf(:,:)=0
      itcon_3Dsrc(:,:)=0
      itcon_decay(:)=0
      itcon_wt(:)=0
#ifdef TRACERS_DRYDEP
      itcon_dd(:,:)=0
#endif

      k = 0
      do n=1,ntm
        select case (trname(n))

        case ('Air','CFCn', 'SF6', 'SF6_c')
               ! nothing to do: use defaults

        case ('CO2n')
          qcon(10) = .true.
          qsum(10) = .true.

        case ('Rn222')
          itcon_decay(n) = 13
          qcon(itcon_decay(n)) = .true.; conpts(1) = 'DECAY'
          qsum(itcon_decay(n)) = .true.

        case ('CO2')
          itcon_surf(1,N) = 13
          qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'FossilFuel'
          itcon_surf(2,N) = 14
          qcon(itcon_surf(2,N)) = .true.; conpts(2) = 'Fertilization'
          itcon_surf(3,N) = 15
          qcon(itcon_surf(3,N)) = .true.; conpts(3) = 'Forest Regrowth'
          itcon_surf(4,N) = 16
          qcon(itcon_surf(4,N)) = .true.; conpts(4) = 'Land Use'
          itcon_surf(5,N) = 17
          qcon(itcon_surf(5,N)) = .true.; conpts(5) = 'Ecosystem Exch'
          itcon_surf(6,N) = 18
          qcon(itcon_surf(6,N)) = .true.; conpts(6) = 'Ocean Exch'

        case ('N2O')   ! two versions dependent on configuration
#ifdef TRACERS_SPECIAL_Lerner
          itcon_surf(1,N) = 13
          qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Reset in L1'
          itcon_3Dsrc(1,N) = 14
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(2) = 'Strat. Chem'
          qsum(itcon_3Dsrc(1,N)) = .true.
#endif
#ifdef TRACERS_SPECIAL_Shindell
          kt_power_change(n) = -14

          g=13; itcon_3Dsrc(nChemistry,N) = g
          qcon(itcon_3Dsrc(nChemistry,N))=.true.; conpts(g-12)=
     *         'Chemistry'
          qsum(itcon_3Dsrc(nChemistry,N))=.true.
          g=g+1; itcon_3Dsrc(nOverwrite,N) = g
          qcon(itcon_3Dsrc(nOverwrite,N))=.true.; conpts(g-12)
     *         ='Overwrite'
          qsum(itcon_3Dsrc(nOverwrite,N)) = .true.
          do kk=1,ntsurfsrc(n)
            g=g+1; itcon_surf(kk,N) = g
            qcon(itcon_surf(kk,N))=.true.; conpts(g-12)=trim(ssname(N,kk
     *           ))
          enddo
#endif

        case ('CFC11')
          itcon_surf(1,N) = 13
          qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'L1 Source'
          itcon_3Dsrc(1,N) = 14
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(2) = 'Strat. Chem'
          qsum(itcon_3Dsrc(1,N)) = .true.

        case ('14CO2')
          itcon_surf(1,N) = 13
          qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Bombs and drift'

        case ('CH4')            ! two versions
#ifdef TRACERS_SPECIAL_Shindell
          kt_power_change(n) = -13

          g=13; itcon_3Dsrc(nChemistry,N) = g
          qcon(itcon_3Dsrc(nChemistry,N)) =.true.
          conpts(g-12)='Chemistry'
          qsum(itcon_3Dsrc(nChemistry,N)) = .true.
          g=g+1; itcon_3Dsrc(nOverwrite,N) = g
          qcon(itcon_3Dsrc(nOverwrite,N))=.true.
          conpts(g-12)='Overwrite'
          qsum(itcon_3Dsrc(nOverwrite,N)) = .true.
          do kk=1,ntsurfsrc(n)
            g=g+1; itcon_surf(kk,N) = g
            qcon(itcon_surf(kk,N))=.true.
            conpts(g-12)=trim(ssname(N,kk))
          enddo
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            g=g+1; itcon_dd(n,1)= g
            qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'DRY DEP'
            qsum(itcon_dd(n,1)) = .false.
          end if
#endif
           g=g+1; itcon_3Dsrc(nBiomass,N) = g
           qcon(g) = .true.; conpts(g-12) = 'Biomass src'
           qsum(g) = .true.
#else  /* not TRACERS_SPECIAL_Shindell */
          itcon_surf(1,N) = 13
          qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Animal source'
          itcon_surf(2,N) = 14
          qcon(itcon_surf(2,N)) = .true.; conpts(2) = 'Coal Mine source'
          itcon_surf(3,N) = 15
          qcon(itcon_surf(3,N)) = .true.; conpts(3) = 'Gas Leak source'
          itcon_surf(4,N) = 16
          qcon(itcon_surf(4,N)) = .true.; conpts(4) = 'Gas Vent source'
          itcon_surf(5,N) = 17
          qcon(itcon_surf(5,N)) = .true.; conpts(5) = 'City Dump source'
          itcon_surf(6,N) = 18
          qcon(itcon_surf(6,N)) = .true.; conpts(6) = 'Soil sink'
          itcon_surf(7,N) = 19
          qcon(itcon_surf(7,N)) = .true.; conpts(7) = 'Termite Source'
          itcon_surf(8,N) = 20
          qcon(itcon_surf(8,N)) = .true.; conpts(8) = 'Coal Combustion'
          itcon_surf(9,N) = 21
          qcon(itcon_surf(9,N)) = .true.; conpts(9) = 'Ocean source'
          itcon_surf(10,N) = 22
          qcon(itcon_surf(10,N)) = .true.; conpts(10) = 'Lake source'
          itcon_surf(11,N) = 23
          qcon(itcon_surf(11,N)) = .true.; conpts(11)
     *         ='Misc. Ground source'
          itcon_surf(12,N) = 24
          qcon(itcon_surf(12,N)) = .true.; conpts(12) =
     *         'Biomass Burning'
          itcon_surf(13,N) = 25
          qcon(itcon_surf(13,N)) = .true.; conpts(13) = 'Rice source'
          itcon_surf(14,N) = 26
          qcon(itcon_surf(14,N)) = .true.; conpts(14) =
     *         'Wetlands+Tundra'
          itcon_3Dsrc(1,N) = 27
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(15) = 'Tropos. Chem'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_3Dsrc(2,N) = 28
          qcon(itcon_3Dsrc(2,N)) = .true.; conpts(16) = 'Stratos. Chem'
          qsum(itcon_3Dsrc(2,N)) = .true.
#endif /* TRACERS_SPECIAL_Shindell */

        case ('O3')
          itcon_surf(1,N) = 13
          qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Deposition'
          itcon_3Dsrc(1,N) = 14
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(2) = 'Stratos. Chem'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_3Dsrc(2,N) = 15
          qcon(itcon_3Dsrc(2,N)) = .true.; conpts(3) = 'Trop.Chem Prod'
          qsum(itcon_3Dsrc(2,N)) = .true.
          itcon_3Dsrc(3,N) = 16
          qcon(itcon_3Dsrc(3,N)) = .true.; conpts(4) = 'Trop.Chem Loss'
          qsum(itcon_3Dsrc(3,N)) = .true.

        case ('Ox','N2O5','HNO3','H2O2','CH3OOH','HCHO','HO2NO2','PAN'
     *       ,'AlkylNit','ClOx','BrOx','HCl','HOCl','ClONO2','HBr'
     *       ,'HOBr','BrONO2','CFC','NOx','CO','Isoprene','Alkenes'
     *       ,'Paraffin','stratOx','Terpenes') ! N2O done above
          select case (trname(n))
            case ('N2O5','CH3OOH','HCHO','HO2NO2','PAN','AlkylNit','CFC'
     *           ,'ClOx','BrOx','HCl','HOCl','ClONO2','HBr','HOBr'
     *           ,'BrONO2','NOx')
              kt_power_change(n) = -14
            case ('HNO3','H2O2','CO','Isoprene','Alkenes','Paraffin'
     *           ,'Terpenes')
              kt_power_change(n) = -13
            case default
              kt_power_change(n) = -12
          end select

          g=13; itcon_3Dsrc(nChemistry,N) = g
          qcon(itcon_3Dsrc(nChemistry,N)) = .true.
          conpts(g-12) = 'Chemistry'
          qsum(itcon_3Dsrc(nChemistry,N)) = .true.
          g=g+1; itcon_3Dsrc(nOverwrite,N) = g
          qcon(itcon_3Dsrc(nOverwrite,N)) = .true.
          conpts(g-12)='Overwrite'
          qsum(itcon_3Dsrc(nOverwrite,N)) = .true.
          select case(trname(n))
            case ('NOx')
              g=g+1; itcon_3Dsrc(nOther,N) = g
              qcon(itcon_3Dsrc(nOther,N)) = .true.
              conpts(g-12) = 'Lightning'
              qsum(itcon_3Dsrc(nOther,N)) = .true.
              g=g+1; itcon_3Dsrc(nAircraft,N) = g
              qcon(itcon_3Dsrc(nAircraft,N)) = .true.
              conpts(g-12) = 'Aircraft'
              qsum(itcon_3Dsrc(nAircraft,N)) = .true.
          end select
          select case(trname(n))
            case('NOx','CO','Alkenes','Paraffin')
              g=g+1; itcon_3Dsrc(nBiomass,N) = g
              qcon(g) = .true.; conpts(g-12) = 'Biomass src'
              qsum(g) = .true.
          end select
#ifdef TRACERS_NITRATE
          select case (trname(n))
            case ('HNO3')
              g=g+1; itcon_3Dsrc(3,N) = g
              qcon(g) = .true.; conpts(g-12)='Nitrate Chemistry'
              qsum(g) = .true.
          end select
#endif
          do kk=1,ntsurfsrc(n)
            g=g+1; itcon_surf(kk,N) = g
            qcon(itcon_surf(kk,N))=.true.
            conpts(g-12)=trim(ssname(N,kk))
          end do
#ifdef TRACERS_WATER
          if(dowetdep(n)) then
            g=g+1; itcon_mc(n) = g
            qcon(itcon_mc(n)) = .true.  ; conpts(g-12) = 'MOIST CONV'
            g=g+1; itcon_ss(n) = g
            qcon(itcon_ss(n)) = .true.  ; conpts(g-12) = 'LS COND'
          end if
#endif
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            g=g+1; itcon_dd(n,1)=g
            qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'DRY DEP'
          end if
#endif

        case ('codirect')
          kt_power_change(n) = -13
          g=13; itcon_decay(n) = g
          qcon(itcon_decay(n)) = .true.
          conpts(g-12) = 'DECAY'
          qsum(itcon_decay(n)) = .true.
          do kk=1,ntsurfsrc(n)
            g=g+1; itcon_surf(kk,n) = g
            qcon(itcon_surf(kk,n))=.true.
            conpts(g-12)=trim(ssname(n,kk))
          end do

#ifdef TRACERS_AEROSOLS_SOA
        case ('isopp1g','isopp1a','isopp2g','isopp2a',
     &        'apinp1g','apinp1a','apinp2g','apinp2a')
          g=13; itcon_3Dsrc(nChemistry,N) = g
          qcon(itcon_3Dsrc(nChemistry,N)) = .true.
          conpts(g-12) = 'Chemistry'
          qsum(itcon_3Dsrc(nChemistry,N)) = .true.
          g=g+1; itcon_mc(n) = g
          qcon(itcon_mc(n)) = .true.  ; conpts(g-12) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          g=g+1; itcon_ss(n) = g
          qcon(itcon_ss(n)) = .true.  ; conpts(g-12) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            g=g+1; itcon_dd(n,1) = g
            qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            g=g+1; itcon_dd(n,2) = g
            qcon(itcon_dd(n,2)) = .true. ; conpts(g-12) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif
#endif  /* TRACERS_AEROSOLS_SOA */

        case ('GLT')
          kt_power_change(n) = -17

          g=13 ; itcon_3Dsrc(1,N) = g
          qcon(itcon_3Dsrc(1,N)) = .true.
          conpts(g-12) = 'L1 overwriting'
          qsum(itcon_3Dsrc(1,N)) = .true.

        case ('Water','H2O18', 'HDO', 'H2O17' )
          itcon_mc(n) = 13
          qcon(itcon_mc(n)) = .true.  ; conpts(1) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) = 14
          qcon(itcon_ss(n)) = .true.  ; conpts(2) = 'LS COND'
          qsum(itcon_ss(n)) = .false.

        case ('HTO')
          itcon_mc(n) = 13
          qcon(itcon_mc(n)) = .true.  ; conpts(1) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) = 14
          qcon(itcon_ss(n)) = .true.  ; conpts(2) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
          itcon_decay(n) = 15
          qcon(itcon_decay(n)) = .true.; conpts(3) = 'DECAY'
          qsum(itcon_decay(n)) = .true.

        case ('DMS')
          itcon_surf(1,N) = 13
          qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Ocean Src'
          qsum(itcon_surf(1,N))=.false.
          itcon_3Dsrc(1,N) = 14
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(2) = 'Chem'
          qsum(itcon_3Dsrc(1,N))= .true.

        case ('MSA')
          itcon_3Dsrc(1,N) = 13
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Chem'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_mc(n) = 14
          qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) = 15
          qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=16
            qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=17
            qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif

        case ('SO2')
          g=13; itcon_3Dsrc(nVolcanic,N) = g
          qcon(g) = .true.; conpts(g-12) = 'Volcanic src'
          qsum(g) = .true.
          g=g+1; itcon_3Dsrc(nAircraft,N) = g
          qcon(g) = .true.; conpts(g-12) = 'Aircraft src'
          qsum(g)=.true.
          g=g+1; itcon_3Dsrc(nBiomass,N) = g
          qcon(g) = .true.; conpts(g-12) = 'Biomass src'
          qsum(g)=.true.
          g=g+1; itcon_3Dsrc(nChemistry,N) = g
          qcon(g) = .true.; conpts(g-12) = 'Chem src'
          qsum(g) = .true.
          g=g+1; itcon_3Dsrc(nChemloss,N) = g
          qcon(g) = .true.; conpts(g-12) = 'Chem sink'
          qsum(g) = .true.
          do kk=1,ntsurfsrc(n)
            g=g+1; itcon_surf(kk,N) = g
            qcon(itcon_surf(kk,N))=.true.
            conpts(g-12)=trim(ssname(N,kk))
            qsum(g)=.false.
          end do
          g=g+1; itcon_mc(n) = g
          qcon(g) = .true.; conpts(g-12) = 'MOIST CONV'
          qsum(g) = .false.
          g=g+1; itcon_ss(n) = g
          qcon(g) = .true.; conpts(g-12) = 'LS COND'
          qsum(g) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            g=g+1; itcon_dd(n,1)= g
            qcon(g) = .true.; conpts(g-12) = 'TURB DEP'
            qsum(g) = .false.
          end if
#endif

        case ('SO4')
          g=13; itcon_3Dsrc(nChemistry,N) = g
          qcon(g) = .true.; conpts(g-12) = 'Gas phase src'
          qsum(g) = .true.
          g=g+1; itcon_3Dsrc(nVolcanic,N) = g
          qcon(g) = .true.; conpts(g-12) = 'Volcanic src'
          qsum(g) = .true.
          g=g+1; itcon_3Dsrc(nBiomass,N) = g
          qcon(g) = .true.; conpts(g-12) = 'Biomass src'
          qsum(g) = .true.
          do kk=1,ntsurfsrc(n_SO2)
            g=g+1; itcon_surf(kk,N) = g
            qcon(itcon_surf(kk,N))=.true.
            conpts(g-12)=trim(ssname(n_SO2,kk))
            qsum(g)=.false.
          end do
          g=g+1; itcon_mc(n) = g
          qcon(g) = .true.; conpts(g-12) = 'MOIST CONV'
          qsum(g) = .false.
          g=g+1; itcon_ss(n) = g
          qcon(g) = .true.; conpts(g-12) = 'LS COND'
          qsum(g) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            g=g+1; itcon_dd(n,1)= g
            qcon(g) = .true.; conpts(g-12) = 'TURB DEP'
            qsum(g) = .false.
            g=g+1; itcon_dd(n,2)= g
            qcon(g) = .true.; conpts(g-12) = 'GRAV SET'
            qsum(g) = .false.
          end if
#endif

        case ('BCII')
          itcon_3Dsrc(1,N) = 13
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Aging loss'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_3Dsrc(2,N) = 14
          qcon(itcon_3Dsrc(2,N)) = .true.; conpts(2) =
     *         'Aircraft source'
          qsum(itcon_3Dsrc(2,N)) = .true.
          itcon_surf(1,N) = 15
          qcon(itcon_surf(1,N)) = .true.; conpts(3) = 'Industrial src'
          qsum(itcon_surf(1,N)) = .false.
          itcon_surf(2,N) = 16
          qcon(itcon_surf(2,N)) = .true.; conpts(4) = 'ship src'
          qsum(itcon_surf(2,N)) = .false.
          itcon_mc(n) = 17
          qcon(itcon_mc(n)) = .true.  ; conpts(5) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) = 18
          qcon(itcon_ss(n)) = .true.  ; conpts(6) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=19
            qcon(itcon_dd(n,1)) = .true. ; conpts(7) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=20
            qcon(itcon_dd(n,2)) = .true. ; conpts(8) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif

        case ('BCIA')
          itcon_3Dsrc(1,N) = 13
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Aging source'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_3Dsrc(2,N) = 14
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(2) =
     *         'Aircraft source'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_mc(n) = 15
          qcon(itcon_mc(n)) = .true.  ; conpts(3) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) = 16
          qcon(itcon_ss(n)) = .true.  ; conpts(4) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=17
            qcon(itcon_dd(n,1)) = .true. ; conpts(5) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=18
            qcon(itcon_dd(n,2)) = .true. ; conpts(6) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif

        case ('OCIA')
          itcon_3Dsrc(1,N) = 13
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Aging source'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_mc(n) = 14
          qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) = 15
          qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=16
            qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=17
            qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif

        case ('BCB', 'OCB')
          g=13; itcon_3Dsrc(nBiomass,N) = g
          qcon(g) = .true.; conpts(g-12) = 'Biomass src'
          qsum(g) = .true.
          g=g+1; itcon_mc(n) = g
          qcon(g) = .true.; conpts(g-12) = 'MOIST CONV'
          qsum(g) = .false.
          g=g+1; itcon_ss(n) = g
          qcon(g) = .true.; conpts(g-12) = 'LS COND'
          qsum(g) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            g=g+1; itcon_dd(n,1)=g
            qcon(g) = .true.; conpts(g-12) = 'TURB DEP'
            qsum(g) = .false.
            g=g+1; itcon_dd(n,2)=g
            qcon(g) = .true.; conpts(g-12) = 'GRAV SET'
            qsum(g) = .false.
          end if
#endif

        case ('OCII')
          g=13; itcon_3Dsrc(1,N) = g
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(g-12) = 'Aging loss'
          qsum(itcon_3Dsrc(1,N)) = .true.
          g=g+1; itcon_surf(1,N) = g
          qcon(itcon_surf(1,N)) = .true.; conpts(g-12)='Industrial src'
          qsum(itcon_surf(1,N)) = .false.
#ifndef TRACERS_AEROSOLS_SOA
          g=g+1; itcon_surf(3,N) = g
          qcon(itcon_surf(3,N)) = .true.; conpts(g-12) = 'Terpene src'
          qsum(itcon_surf(3,N)) = .false.
#endif  /* TRACERS_AEROSOLS_SOA */
          g=g+1; itcon_surf(2,N) = g
          qcon(itcon_surf(2,N)) = .true.; conpts(g-12) = 'ship src'
          qsum(itcon_surf(2,N)) = .false.
          g=g+1; itcon_mc(n) = g
          qcon(itcon_mc(n)) = .true.  ; conpts(g-12) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          g=g+1; itcon_ss(n) = g
          qcon(itcon_ss(n)) = .true.  ; conpts(g-12) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            g=g+1; itcon_dd(n,1)= g
            qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            g=g+1; itcon_dd(n,2)= g
            qcon(itcon_dd(n,2)) = .true. ; conpts(g-12) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif

        case ('OCI3')
          g=13; itcon_3Dsrc(1,N) = g
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(g-12) = 'Emission src'
          qsum(itcon_3Dsrc(1,N)) = .true.
          g=g+1; itcon_3Dsrc(2,N) = g
          qcon(itcon_3Dsrc(2,N)) = .true.; conpts(g-12) = 'Aging loss'
          qsum(itcon_3Dsrc(2,N)) = .true.
#ifndef TRACERS_AEROSOLS_SOA
          g=g+1; itcon_surf(1,N) = g
          qcon(itcon_surf(1,N)) = .true.; conpts(g-12) = 'Terpene src'
          qsum(itcon_surf(1,N)) = .false.
#endif  /* TRACERS_AEROSOLS_SOA */
          g=g+1; itcon_mc(n) = g
          qcon(itcon_mc(n)) = .true.  ; conpts(g-12) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          g=g+1; itcon_ss(n) = g
          qcon(itcon_ss(n)) = .true.  ; conpts(g-12) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            g=g+1; itcon_dd(n,1)= g
            qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            g=g+1; itcon_dd(n,2)= g
            qcon(itcon_dd(n,2)) = .true. ; conpts(g-12) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif

        case ('OCI1','OCI2')
          itcon_3Dsrc(1,N) = 13
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Emission src'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_3Dsrc(2,N) = 14
          qcon(itcon_3Dsrc(2,N)) = .true.; conpts(2) = 'Aging loss'
          qsum(itcon_3Dsrc(2,N)) = .true.
          itcon_mc(n) = 15
          qcon(itcon_mc(n)) = .true.  ; conpts(3) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) = 16
          qcon(itcon_ss(n)) = .true.  ; conpts(4) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=17
            qcon(itcon_dd(n,1)) = .true. ; conpts(5) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=18
            qcon(itcon_dd(n,2)) = .true. ; conpts(6) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif

        case ('OCA4')
          itcon_3Dsrc(1,N) = 13
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Emission src'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_mc(n) = 14
          qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) = 15
          qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=16
            qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=17
            qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif

        case ('OCA1','OCA2','OCA3')
          itcon_3Dsrc(1,N) = 13
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Decay src'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_mc(n) = 14
          qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) = 15
          qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=16
            qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=17
            qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif

        case ('SO4_d1', 'SO4_d2','SO4_d3','N_d1','N_d2','N_d3')
          itcon_3Dsrc(1,N) = 13
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) =
     *         'Gas phase change'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_mc(n) =14
          qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) =15
          qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=16
            qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=17
            qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif

        case ('NH3','H2SO4')
          g=13; itcon_3Dsrc(nChemistry,N) = g
          qcon(g) = .true.; conpts(g-12) = 'Gas phase change'
          qsum(g) = .true.
          select case (trname(n))
          case ('NH3')
            g=g+1; itcon_3Dsrc(nBiomass,N) = g
            qcon(g) = .true.; conpts(g-12) = 'Biomass src'
            qsum(g) = .true.
            do kk=1,ntsurfsrc(n)
              g=g+1; itcon_surf(kk,N) = g
              qcon(itcon_surf(kk,N))=.true.
              conpts(g-12)=trim(ssname(N,kk))
              qsum(g)=.false.
            end do
          end select
          g=g+1; itcon_mc(n) =g
          qcon(g) = .true.  ; conpts(g-12) = 'MOIST CONV'
          qsum(g) = .false.
          g=g+1; itcon_ss(n) =g
          qcon(g) = .true.  ; conpts(g-12) = 'LS COND'
          qsum(g) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            g=g+1; itcon_dd(n,1)=g
            qcon(g) = .true. ; conpts(g-12) = 'TURB DEP'
            qsum(g) = .false.
            g=g+1; itcon_dd(n,2)=g
            qcon(g) = .true. ; conpts(g-12) = 'GRAV SET'
            qsum(g) = .false.
          end if
#endif

        case ('NH4', 'NO3p')
          itcon_3Dsrc(1,N) = 13
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) =
     *         'Gas phase change'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_mc(n) =14
          qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) =15
          qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=16
            qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=17
            qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif

        case ('Be7', 'Be10')
          itcon_3Dsrc(1,N) =13
          qcon(itcon_3Dsrc(1,N)) = .true.  ; conpts(1) = 'COSMO SRC'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_mc(n) =14
          qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) =15
          qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=16
            qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=17
            qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif
          if (trname(n).eq."Be7") then
            itcon_decay(n) = 18
            qcon(itcon_decay(n)) = .true.; conpts(6) = 'DECAY'
            qsum(itcon_decay(n)) = .true.
          end if

        case ('Pb210')
          itcon_3Dsrc(1,N) =13
          qcon(itcon_3Dsrc(1,N)) = .true.  ; conpts(1) = 'RADIO SRC'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_mc(n) =14
          qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) =15
          qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=16
            qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=17
            qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif
          itcon_decay(n) = 18
          qcon(itcon_decay(n)) = .true.; conpts(6) = 'DECAY'
          qsum(itcon_decay(n)) = .true.

        case ('H2O2_s')
          itcon_3Dsrc(1,N) = 13
          qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Gas phase src'
          qsum(itcon_3Dsrc(1,N)) = .true.
          itcon_3Dsrc(2,N) = 14
          qcon(itcon_3Dsrc(2,N)) = .true.; conpts(2) =
     *         'Gas phase sink'
          qsum(itcon_3Dsrc(2,N)) = .true.
          itcon_mc(n) =15
          qcon(itcon_mc(n)) = .true.  ; conpts(3) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) =16
          qcon(itcon_ss(n)) = .true.  ; conpts(4) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=17
            qcon(itcon_dd(n,1)) = .true. ; conpts(5) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
          end if
#endif

        case ('seasalt1','seasalt2', 'Clay','Silt1','Silt2','Silt3'
     *       ,'Silt4','ClayIlli','ClayKaol','ClaySmec','ClayCalc'
     *       ,'ClayQuar','Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema'
     *       ,'Sil1Gyps','Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema'
     *       ,'Sil2Gyps','Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema'
     *       ,'Sil3Gyps','Sil1QuHe','Sil2QuHe','Sil3QuHe')
          itcon_mc(n) =13
          qcon(itcon_mc(n)) = .true.  ; conpts(1) = 'MOIST CONV'
          qsum(itcon_mc(n)) = .false.
          itcon_ss(n) =14
          qcon(itcon_ss(n)) = .true.  ; conpts(2) = 'LS COND'
          qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            itcon_dd(n,1)=15
            qcon(itcon_dd(n,1)) = .true. ; conpts(3) = 'TURB DEP'
            qsum(itcon_dd(n,1)) = .false.
            itcon_dd(n,2)=16
            qcon(itcon_dd(n,2)) = .true. ; conpts(4) = 'GRAV SET'
            qsum(itcon_dd(n,2)) = .false.
          end if
#endif
#ifndef TRACERS_WATER
          itcon_wt(n)=17
          qcon(itcon_wt(n)) = .true. ; conpts(5) = 'WET DEP'
          qsum(itcon_wt(n)) = .false.
#endif  /* not TRACERS_WATER */

c- Species including AMP  emissions - 2D sources and 3D sources
        case('M_AKK_SU','M_ACC_SU','M_OCC_OC','M_BC1_BC',
     *       'M_SSA_SS','M_SSC_SS','M_SSS_SS','M_DD1_DU','M_DD2_DU',
     *       'M_BOC_BC','M_BOC_OC',
     *       'M_NO3   ','M_NH4   ','M_H2O   ','M_DD1_SU',
     *       'M_DS1_SU','M_DS1_DU','M_DD2_SU',
     *       'M_DS2_SU','M_DS2_DU','M_SSA_SU',
     *       'M_OCC_SU','M_BC1_SU',
     *       'M_BC2_SU','M_BC2_BC','M_BC3_SU',
     *       'M_BC3_BC','M_DBC_SU','M_DBC_BC','M_DBC_DU',
     *       'M_BOC_SU',
     *       'M_BCS_SU','M_BCS_BC','M_MXX_SU','M_MXX_BC',
     *       'M_MXX_OC','M_MXX_DU','M_MXX_SS','M_OCS_SU',
     *       'M_OCS_OC','M_SSS_SU')
          g=13; itcon_3Dsrc(nChemistry,N) = g
          qcon(g) = .true.;conpts(g-12)='Gas phase change'
          qsum(g) = .true.
          g=g+1; itcon_mc(n) = g
          qcon(g) = .true.  ; conpts(g-12) = 'MOIST CONV'
          qsum(g) = .false.
          g=g+1; itcon_ss(n) = g
          qcon(g) = .true.  ; conpts(g-12) = 'LS COND'
          qsum(g) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            g=g+1; itcon_dd(n,1)= g
            qcon(g) = .true. ; conpts(g-12) = 'TURB DEP'
            qsum(g) = .false.
            g=g+1; itcon_dd(n,2)= g
            qcon(g) = .true. ; conpts(g-12) = 'GRAV SET'
            qsum(g) = .false.
          end if
#endif
          g=g+1; itcon_3Dsrc(nChemistry,n) = g
          qcon(g) = .true.; conpts(g-12) = 'AMP source'
          qsum(g) = .true.
          select case (trname(n))
            case ('M_SSA_SS','M_SSC_SS','M_SSS_SS','M_DD1_DU'
     *           ,'M_DD2_DU')
              g=g+1; itcon_surf(1,n) = g
              qcon(g) = .true.; conpts(g-12) = 'Emission 2D AMP'
              qsum(g) = .true.
            case ('M_AKK_SU','M_ACC_SU',
     &            'M_BC1_BC','M_OCC_OC','M_BOC_BC','M_BOC_OC')
              select case (trname(n))
              case ('M_AKK_SU','M_ACC_SU')
                do kk=1,ntsurfsrc(n_SO2)
                  g=g+1; itcon_surf(kk,N) = g
                  qcon(itcon_surf(kk,N))=.true.
                  conpts(g-12)=trim(ssname(n_SO2,kk))//' 2D AMP'
                  qsum(g)=.true.
                end do
                g=g+1; itcon_3Dsrc(nVolcanic,n) = g
                qcon(g) = .true.; conpts(g-12) = 'Volcano 3D AMP'
                qsum(g) = .true.
              case ('M_BC1_BC','M_OCC_OC','M_BOC_BC','M_BOC_OC')
                do kk=1,ntsurfsrc(n)
                  g=g+1; itcon_surf(kk,N) = g
                  qcon(itcon_surf(kk,N))=.true.
                  conpts(g-12)=trim(ssname(n,kk))//' Emission 2D AMP'
                  qsum(g)=.true.
                end do
              end select
              g=g+1; itcon_3Dsrc(nBiomass,n) = g
              qcon(g) = .true.; conpts(g-12) = 'Biomass 3D AMP'
              qsum(g) = .true.
          end select
c Processes AMP Budget
          g=g+1; itcon_AMP(1,n)= g
          qcon(g) = .true.; conpts(g-12)='P1 Nucleation'
          qsum(g) = .true.
          g=g+1; itcon_AMP(2,n)= g
          qcon(g) = .true.; conpts(g-12)='P2 Coagulation'
          qsum(g) = .true.
          g=g+1; itcon_AMP(3,n)= g
          qcon(g) = .true.;conpts(g-12)='P3 Condensation'
          qsum(g) = .true.
          g=g+1; itcon_AMP(4,n)= g
          qcon(g) = .true.;conpts(g-12)='P4 Incloud'
          qsum(g) = .true.
          g=g+1; itcon_AMP(5,n)= g
          qcon(g) = .true.;conpts(g-12)='P5 Intermode Loss'
          qsum(g) = .true.
          g=g+1; itcon_AMP(6,n)= g
          qcon(g) = .true.;conpts(g-12)='P6 Mode Transf'
          qsum(g) = .true.
          g=g+1; itcon_AMP(7,n)= g
          qcon(g) = .true.;conpts(g-12)='P7 AMP Budget'
          qsum(g) = .true.

        case ('N_AKK_1 ','N_ACC_1 ','N_DD1_1 ','N_DS1_1 ','N_DD2_1 '
     *       ,'N_DS2_1 ','N_OCC_1 ','N_BC1_1 ','N_BC2_1 ','N_BC3_1 '
     *       ,'N_DBC_1 ','N_BOC_1 ','N_BCS_1 ','N_MXX_1 ','N_OCS_1 ')

          kt_power_change(n) = 5
          kt_power_inst(n) = 3

          g=13; itcon_3Dsrc(nChemistry,N) = g
          qcon(g) = .true.; conpts(g-12) = 'Gas phase change'
          qsum(g) = .true.
          g=g+1; itcon_mc(n) = g
          qcon(g) = .true.  ; conpts(g-12) = 'MOIST CONV'
          qsum(g) = .false.
          g=g+1; itcon_ss(n) = g
          qcon(g) = .true.  ; conpts(g-12) = 'LS COND'
          qsum(g) = .false.
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            g=g+1; itcon_dd(n,1)= g
            qcon(g) = .true. ; conpts(g-12) = 'TURB DEP'
            qsum(g) = .false.
            g=g+1; itcon_dd(n,2)= g
            qcon(g) = .true. ; conpts(g-12) = 'GRAV SET'
            qsum(g) = .false.
          end if
#endif
          g=g+1; itcon_AMPm(2,n)= g
          qcon(g) = .true. ; conpts(g-12) = 'Mode AktivPart'
          qsum(g) = .false.
c     Processes AMP Budget
          g=g+1; itcon_AMP(1,n)= g
          qcon(g) = .true. ; conpts(g-12) = 'P1 Nucleation'
          qsum(g) = .true.
          g=g+1; itcon_AMP(2,n)= g
          qcon(g) = .true. ; conpts(g-12) = 'P2 Coagulation'
          qsum(g) = .true.
          g=g+1; itcon_AMP(3,n)= g
          qcon(g) = .true.;conpts(g-12) ='P3 NOTHING'
          qsum(g) = .true.
          g=g+1; itcon_AMP(4,n)= g
          qcon(g) = .true.;conpts(g-12)='P4 Intermode Coag'
          qsum(g) = .true.
          g=g+1; itcon_AMP(5,n)= g
          qcon(g) = .true.;conpts(g-12)='P5 Intramode Tr'
          qsum(g) = .true.
          g=g+1; itcon_AMP(6,n)= g
          qcon(g) = .true.;conpts(g-12)='P6 Mode Transf'
          qsum(g) = .true.
          g=g+1; itcon_AMP(7,n)= g
          qcon(g) = .true. ; conpts(g-12) = 'P7 AMP Budget'
          qsum(g) = .true.

        end select

        scale_inst(n)   = 10d0**(-kt_power_inst(n))
        scale_change(n) = 10d0**(-kt_power_change(n))
        inst_unit(n) = unit_string(kt_power_inst(n),  'kg/m^2)')
        sum_unit(n)  = unit_string(kt_power_change(n),'kg/m^2 s)')

        CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *       sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
        qcon(13:) = .false.     ! reset to defaults for next tracer
        qsum(13:) = .false.     ! reset to defaults for next tracer
        qcon(10)  = .false.     ! reset to defaults for next tracer
        qsum(10)  = .false.     ! reset to defaults for next tracer

        natmtrcons=N
      end do
#endif /* TRACERS_ON */

      return
      end subroutine init_tracer_cons_diag

      subroutine init_jls_diag
!@sum init_jls_diag Initialise zonal mean/height tracer diags
!@auth Gavin Schmidt
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      USE CONSTANT, only: sday
      USE MODEL_COM, only: dtsrc,byim
      USE TRACER_COM
      USE DIAG_COM
#ifdef TRACERS_ON
      USE TRDIAG_COM
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      use tracers_dust, only: nDustEmjl, nDustEm2jl, nDustEv1jl,
     &   nDustEv2jl, nDustWthjl, imDust
#endif
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
      USE CLOUDS, ONLY : diag_wetdep
#endif
#endif /* TRACERS_ON */
      implicit none
      integer k,n,kk,ltop
      character*50 :: unit_string

C**** Please note that short names for diags i.e. sname_jls are used
C**** in special ways and MUST NOT contain spaces, commas or % signs.
C**** Underscores and minus signs are allowed.

C**** Define a max layer for some optionally trop/strat tracers
      LTOP = LM

#ifdef TRACERS_ON
C**** Tracer sources and sinks
C**** Defaults for jls (sources, sinks, etc.)
C**** These need to be 'hand coded' depending on circumstances
      do k=1,ktajls             ! max number of sources and sinks
        jgrid_jls(k) = 1
        jwt_jls(k) = 1
        ia_jls(k) = ia_src
        scale_jls(k) = 1./DTsrc
      end do
      jls_grav=0
#ifdef TRACERS_WATER
C**** set defaults for some precip/wet-dep related diags
      jls_prec(:,:)=0
#endif

      k = 0
      do n=1,ntm
      select case (trname(n))

      case ('SF6','SF6_c','CFCn')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Layer_1_source_of_'//trname(n)
        lname_jls(k) = trim(trname(n))//' CFC-GRID SOURCE, LAYER 1'
        jls_ltop(k) = 1
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('CO2n')
        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Ocean_Gas_Exchange_'//trname(n)
        lname_jls(k) = trim(trname(n))//' Ocean/Atmos. Gas Exchange'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('Rn222')
        k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'LOSS OF RADON-222 BY DECAY'
        jls_ltop(k) = lm
        jls_power(k) = -26
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k)=3

        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Ground_Source_of_'//trname(n)
        lname_jls(k) = 'RADON-222 SOURCE, LAYER 1'
        jls_ltop(k) = 1
        jls_power(k) = -10
        units_jls(k) = unit_string(jls_power(k),'kg/s')

! keep AIJ and AJL CO2 sources in same order !!
      case ('CO2')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Fossil_fuel_source_'//trname(n)
        lname_jls(k) = 'CO2 Fossil fuel source (Marland)'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'fertilization_sink_'//trname(n)
        lname_jls(k) = 'CO2 fertilization sink (Friedlingstein)'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'Northern_forest_regrowth_'//trname(n)
        lname_jls(k) = 'CO2 Northern forest regrowth sink'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'Land_Use_Modification_'//trname(n)
        lname_jls(k) = 'CO2 from Land use modification (Houton)'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'Ecosystem_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ecosystem exchange (Matthews)'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(6,n) = k
        sname_jls(k) = 'Ocean_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ocean exchange'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('N2O')
#ifdef TRACERS_SPECIAL_Shindell
        do kk=1,ntsurfsrc(n)
          k = k + 1
          jls_source(kk,n) = k
          sname_jls(k) = trim(trname(n))//'_'//trim(ssname(n,kk))
          lname_jls(k) = trname(n)//' source from '//trim(ssname(n,kk))
          jls_ltop(k) = 1
          jls_power(k) = -1
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        end do
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF '//trname(n)//' BY CHEMISTRY'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(nOverwrite,n) = k
        sname_jls(k) = 'overwrite_source_of'//trname(n)
        lname_jls(k) =
     &  'CHANGE OF '//trname(n)//' BY OVERWRITE'
        jls_ltop(k) = 1 ! really L=1 overwrite only
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
#ifdef TRACERS_SPECIAL_Lerner
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF N20 BY RESETTING TO 462.2d-9, L1'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Stratos_chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF N2O BY CHEMISTRY IN STRATOS'
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif

      case ('CFC11')   !!! should start April 1
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF CFC-11 BY SOURCE, L1'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Stratos_chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CFC-11 BY CHEMISTRY IN STRATOS'
        jls_ltop(k) = lm
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('14CO2')   !!! should start 10/16
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF 14CO2 by SINK, L1'
        jls_ltop(k) = 1
        jls_power(k) = -4
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('CH4')
#ifdef TRACERS_SPECIAL_Shindell
        do kk=1,ntsurfsrc(n)
          k = k + 1
          jls_source(kk,n) = k
          sname_jls(k) = trim(trname(n))//'_'//trim(ssname(n,kk))
          lname_jls(k) = trname(n)//' source from '//trim(ssname(n,kk))
          jls_ltop(k) = 1
          jls_power(k) = -1
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        end do
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY CHEMISTRY'
        jls_ltop(k) = LTOP
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(nOverwrite,n) = k
        sname_jls(k) = 'overwrite_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY OVERWRITE'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(nBiomass,n) = k
        sname_jls(k) = 'Biomass_src_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' biomass source'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#else
        k = k + 1
        jls_source(6,n) = k
        sname_jls(k) = 'Soil_sink_'//trname(n)
        lname_jls(k) = 'CH4 sink due to soil absorption'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(7,n) = k
        sname_jls(k) = 'Termite_source_'//trname(n)
        lname_jls(k) = 'CH4 Termite source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(9,n) = k
        sname_jls(k) = 'Ocean_source_'//trname(n)
        lname_jls(k) = 'CH4 Ocean source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(10,n) = k
        sname_jls(k) = 'Fresh_Water_lake_source_'//trname(n)
        lname_jls(k) = 'CH4 Fresh Water lake source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(11,n) = k
        sname_jls(k) = 'Misc_Ground_source_'//trname(n)
        lname_jls(k) = 'CH4 Misc_Ground source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(14,n) = k
        sname_jls(k) = 'Wetlands+Tundra_source_'//trname(n)
        lname_jls(k) = 'CH4 Wetlands+Tundra source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Animal_source_of'//trname(n)
        lname_jls(k) = 'CH4 Animal source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'Coal_Mine_source_'//trname(n)
        lname_jls(k) = 'CH4 Coal Mine source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'Gas_Leak_source_'//trname(n)
        lname_jls(k) = 'CH4 Gas Leak source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'Gas_Venting_source_'//trname(n)
        lname_jls(k) = 'CH4 Gas Venting source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'Municipal_solid_waste_source_'//trname(n)
        lname_jls(k) = 'CH4 Municipal solid waste source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(8,n) = k
        sname_jls(k) = 'Coal_combustion_source_'//trname(n)
        lname_jls(k) = 'CH4 Coal combustion source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(12,n) = k
        sname_jls(k) = 'Biomass_burning_source_'//trname(n)
        lname_jls(k) = 'CH4 Biomass burning source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(13,n) = k
        sname_jls(k) = 'Rice_Cultivation_source_'//trname(n)
        lname_jls(k) = 'CH4 Rice Cultivation source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Tropos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY CHEMISTRY IN TROPOSPHERE'
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Stratos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY CHEMISTRY IN STRATOS'
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif

      case ('O3')
       k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Deposition_L1'//trname(n)
        lname_jls(k) = 'Change of O3 by Deposition in Layer 1'
        jls_ltop(k) = 1
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Strat_Chem_change_'//trname(n)
        lname_jls(k) = 'Change of O3 by Chemistry in Stratos'
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Trop_Chem_Prod_change_'//trname(n)
        lname_jls(k) = 'Change of O3 by Chem Prod. in Troposphere'
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_3Dsource(3,n) = k
        sname_jls(k) = 'Trop_Chem_Loss_change_'//trname(n)
        lname_jls(k) = 'Change of O3 by Chem Loss in Troposphere'
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

#ifdef TRACERS_WATER
C**** generic ones for many water tracers
      case ('Water', 'H2O18', 'HDO', 'HTO', 'H2O17' )
       k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Evap_'//trname(n)
        lname_jls(k) = 'EVAPORATION OF '//trname(n)
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_isrc(2,n) = k
        sname_jls(k) = 'Ocn_Evap_'//trname(n)
        lname_jls(k) = 'OCEAN EVAP OF '//trname(n)
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_prec(1,n)=k
        sname_jls(k) = 'Precip_'//trname(n)
        lname_jls(k) = 'PRECIPITATION OF '//trname(n)
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_prec(2,n)=k
        sname_jls(k) = 'Ocn_Precip_'//trname(n)
        lname_jls(k) = 'OCEAN PRECIP OF '//trname(n)
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm/day')

C**** special one unique to HTO
      if (trname(n).eq."HTO") then
       k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'LOSS OF '//TRIM(trname(n))//' BY DECAY'
        jls_ltop(k) = lm
        jls_power(k) = ntm_power(n)+8
        scale_jls(k) = 1./DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg/s')
      end if
#ifdef FAA_emission 
      if (trname(n).eq."Water") then
          k = k + 1
          jls_3Dsource(1,n) = k
          sname_jls(k) = 'aircraft_source_of'//trname(n)
          lname_jls(k) = 'CHANGE OF '//trname(n)//' BY AIRCRAFT'
          jls_ltop(k) = LM
          jls_power(k) = -2
          units_jls(k) = unit_string(jls_power(k),'kg/s')
       endif
#endif
#endif

!#ifdef TRACERS_NITRATE
!       case ('HNO3')
!        k = k + 1
!        jls_3Dsource(nChemistry,n) = k
!        sname_jls(k) = 'chemistry_nitrat_of'//trname(n)
!        lname_jls(k) = 'CHANGE OF HNO3 BY NITRAT C'
!        jls_ltop(k) = LTOP
!        jls_power(k) = 0
!        units_jls(k) = unit_string(jls_power(k),'kg/s')
!#endif

#ifdef SHINDELL_STRAT_EXTRA
      case ('GLT')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'L1_overwrite_soure_'//trname(n)
        lname_jls(k) = trname(n)//'L1 overwrite source'
        jls_ltop(k) = 1
        jls_power(k) = -5
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif

      case ('codirect')
        do kk=1,ntsurfsrc(n)
          k = k + 1
          jls_source(kk,n) = k
          sname_jls(k) = trim(trname(n))//'_'//trim(ssname(n,kk))
          lname_jls(k) = trim(trname(n))//' source from '//
     &                   trim(ssname(n,kk))
          jls_ltop(k) = 1
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        end do
        k = k + 1
        jls_decay(n) = k   ! decay loss
        sname_jls(k) = 'Decay_of_'//trim(trname(n))
        lname_jls(k) = 'LOSS OF '//trim(trname(n))//' BY DECAY'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('HCl','HOCl','ClONO2','HBr','HOBr','BrONO2','CFC',
     &      'BrOx','ClOx','Alkenes','Paraffin','Isoprene','CO',
     &      'N2O5','HNO3','H2O2','CH3OOH','HCHO','HO2NO2','PAN',
     &      'AlkylNit','Ox','NOx','stratOx','Terpenes')
        do kk=1,ntsurfsrc(n)
          k = k + 1
          jls_source(kk,n) = k
          sname_jls(k) = trim(trname(n))//'_'//trim(ssname(n,kk))
          lname_jls(k) = trim(trname(n))//' source from '//
     &                   trim(ssname(n,kk))
          jls_ltop(k) = 1
          select case(trname(n))
          case ('Paraffin','Isoprene','Terpenes')
            jls_power(k) = 0
          case default
            jls_power(k) = -2
          end select
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        end do
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF '//trname(n)//' BY CHEMISTRY'
        jls_ltop(k) = LM
        select case(trname(n))
        case ('Ox','stratOx')
          jls_power(k) = 1
        case default
          jls_power(k) = -1
        end select
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        select case(trname(n))
        case ('Alkenes','Paraffin','Isoprene','CO','N2O5','HNO3',
     &  'H2O2','CH3OOH','HCHO','HO2NO2','PAN','AlkylNit','Ox',
     &  'Terpenes','NOx','stratOx','BrOx','ClOx')
          k = k + 1
          jls_3Dsource(nOverwrite,n) = k
          sname_jls(k) = 'overwrite_source_of'//trname(n)
          lname_jls(k) =
     &    'CHANGE OF '//trname(n)//' BY OVERWRITE'
          jls_ltop(k) = LM
          jls_power(k) = -1
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        case ('CFC')  ! L=1 overwrite only.
          k = k + 1
          jls_3Dsource(nOverwrite,n) = k
          sname_jls(k) = 'overwrite_source_of'//trname(n)
          lname_jls(k) =
     &    'CHANGE OF '//trname(n)//' BY OVERWRITE'
          jls_ltop(k) = 1 ! L=1 overwrite only
          jls_power(k) = -1
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        end select
        select case(trname(n))
        case('NOx')
          k = k + 1
          jls_3Dsource(nOther,n) = k
          sname_jls(k) = 'lightning_source_of'//trname(n)
          lname_jls(k) = 'CHANGE OF '//trname(n)//' BY LIGHTNING'
          jls_ltop(k) = LM
          jls_power(k) = -2
          units_jls(k) = unit_string(jls_power(k),'kg/s')
          k = k + 1
          jls_3Dsource(nAircraft,n) = k
          sname_jls(k) = 'aircraft_source_of'//trname(n)
          lname_jls(k) = 'CHANGE OF '//trname(n)//' BY AIRCRAFT'
          jls_ltop(k) = LM
          jls_power(k) = -2
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        end select
#ifdef FAA_emission 
        select case(trname(n))
        case('CO','Alkenes','Paraffin')
          k = k + 1
          jls_3Dsource(nAircraft,n) = k
          sname_jls(k) = 'aircraft_source_of'//trname(n)
          lname_jls(k) = 'CHANGE OF '//trname(n)//' BY AIRCRAFT'
          jls_ltop(k) = LM
          jls_power(k) = -2
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        end select
#endif
        select case(trname(n))
        case('NOx','CO','Alkenes','Paraffin')
          k = k + 1
          jls_3Dsource(nBiomass,n) = k
          sname_jls(k) = 'Biomass_src_of_'//trim(trname(n))
          lname_jls(k) = trim(trname(n))//' biomass source'
          jls_ltop(k) = LM
          jls_power(k) = -2
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        end select

#ifdef TRACERS_AEROSOLS_SOA
      case ('isopp1g','isopp2g','apinp1g','apinp2g')
c put in chemical production
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of_'//trname(n)
        lname_jls(k) = 'CHANGE OF '//trname(n)//' BY CHEMISTRY'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('isopp1a','isopp2a','apinp1a','apinp2a')
c put in chemical production
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of_'//trname(n)
        lname_jls(k) = 'CHANGE OF '//trname(n)//' BY CHEMISTRY'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of SOA
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//trname(n)
        lname_jls(k) = 'Gravitational Settling of '//trname(n)
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif  /* TRACERS_AEROSOLS_SOA*/

      case ('DMS')
        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Ocean_source_of'//trname(n)
        lname_jls(k) = 'DMS ocean source'
        jls_ltop(k) = 1
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
C
        k = k + 1
        jls_isrc(2,n) = k
        sname_jls(k) = 'TKE_Contribution'//trname(n)
        lname_jls(k) = 'SGSWSP TKE'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_isrc(3,n) = k
        sname_jls(k) = 'Wet_Conv_Contr'//trname(n)
        lname_jls(k) = 'SGSWSP Wet Conv'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_isrc(4,n) = k
        sname_jls(k) = 'Dry_Conv_Contr'//trname(n)
        lname_jls(k) = 'SGSWSP Dry Conv'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_isrc(5,n) = k
        sname_jls(k) = 'SGSWSP-old'//trname(n)
        lname_jls(k) = 'DMS SGSWP-old/old'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Chemical_sink_of'//trname(n)
        lname_jls(k) = 'DMS chemical loss'
        jls_ltop(k) =LM
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg/s')

       case ('MSA')
c put in chemical production of MSA
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'Chemical production of MSA'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of MSA
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of MSA'
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

       case ('NH3')
c industrial source
        do kk=1,ntsurfsrc(n)
          k = k + 1
          jls_source(kk,n) = k
          sname_jls(k) = trim(trname(n))//'_src_'//trim(ssname(n,kk))
          lname_jls(k) = trim(trname(n))//' source from '//
     &                   trim(ssname(n,kk))
          jls_ltop(k) = 1
          jls_power(k) =0
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        enddo
c biomass burning source
        k = k + 1
        jls_3Dsource(nBiomass,n) = k
        sname_jls(k) = 'Biomass_src_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' biomass source'
        jls_ltop(k) = LM
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg/s')

       case ('SO2')
c industrial source
        do kk=1,ntsurfsrc(n)
          k = k + 1
          jls_source(kk,n) = k
          sname_jls(k) = trim(trname(n))//'_src_'//trim(ssname(n,kk))
          lname_jls(k) = trim(trname(n))//' source from '//
     &                   trim(ssname(n,kk))
          jls_ltop(k) = 1
          jls_power(k) =0
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        enddo
c volcanic production of SO2
        k = k + 1
        jls_3Dsource(nVolcanic,n) = k
        sname_jls(k) = 'volcanic_source_of_'//trname(n)
        lname_jls(k) = 'production of SO2 from volcanos'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c aircraft production of SO2
        k = k + 1
        jls_3Dsource(nAircraft,n) = k
        sname_jls(k) = 'aircraft_source_of_'//trname(n)
        lname_jls(k) = 'production of SO2 from aircraft'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c biomass burning source
        k = k + 1
        jls_3Dsource(nBiomass,n) = k
        sname_jls(k) = 'Biomass_src_of_'//trname(n)
        lname_jls(k) = 'SO2 biomass source'
        jls_ltop(k) = LM
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c put in chemical production of SO2
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'dms_source_of_'//trname(n)
        lname_jls(k) = 'production of SO2 from DMS'
        jls_ltop(k) = LM
        jls_power(k) =  1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c put in chemical sink of SO2
        k = k + 1
        jls_3Dsource(nChemloss,n) = k
        sname_jls(k) = 'chem_sink_of_'//trname(n)
        lname_jls(k) = 'chemical sink of SO2'
        jls_ltop(k) = LM
        jls_power(k) =  1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
c convective chem cloud phase sink of SO2
        k = k + 1
        jls_incloud(1,n) = k
        sname_jls(k) = 'mc_cloud_chem_sink_of'//trname(n)
        lname_jls(k) = 'SO2 used in convective cloud chem'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c stratiform chem cloud phase sink of SO2
        k = k + 1
        jls_incloud(2,n) = k
        sname_jls(k) = 'ss_cloud_chem_sink_of'//trname(n)
        lname_jls(k) = 'SO2 used in stratiform cloud chem'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
        case ('SO4')
c gas phase source of SO4
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'gas_phase_source_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c volcanic source of SO4
        k = k + 1
        jls_3Dsource(nVolcanic,n) = k
        sname_jls(k) = 'volcanic_source_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' volcanic source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifdef FAA_emission
c aircraft production of SO4
        k = k + 1
        jls_3Dsource(nAircraft,n) = k
        sname_jls(k) = 'aircraft_source_of_'//trname(n)
        lname_jls(k) = 'production of SO4 from aircraft'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
c biomass source of SO4
        k = k + 1
        jls_3Dsource(nBiomass,n) = k
        sname_jls(k) = 'biomass_source_of_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' biomass source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#if (defined TRACERS_AEROSOLS_Koch)
c convective cloud phase source of SO4
        k = k + 1
        jls_incloud(1,n) = k
        sname_jls(k) = 'mc_cloud_source_of'//trname(n)
        lname_jls(k) = 'SO4 made in convective clouds'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c stratiform cloud phase source of SO4
        k = k + 1
        jls_incloud(2,n) = k
        sname_jls(k) = 'ss_cloud_source_of'//trname(n)
        lname_jls(k) = 'SO4 made in stratiform clouds'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
c industrial source
        do kk=1,ntsurfsrc(n_SO2)
          k = k + 1
          jls_source(kk,n) = k
          sname_jls(k) = trim(trname(n))//'_src_'//
     &                   trim(ssname(n_SO2,kk))
          lname_jls(k) = trim(trname(n))//' source from '//
     &                   trim(ssname(n_SO2,kk))
          jls_ltop(k) = 1
          jls_power(k) =0
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        enddo
c gravitational settling of SO4
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//trim(trname(n))
        lname_jls(k) = 'Gravitational Settling of '//trim(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

        case ('SO4_d1', 'SO4_d2', 'SO4_d3')
c gas phase source
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'gas_phase_source_of'//trname(n)
        lname_jls(k) = trim(trname(n))//' gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of '//trname(n)
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

        case ('Be7')
c cosmogenic source from file
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Cosmogenic_src_of'//trname(n)
        lname_jls(k) = 'Be7 cosmogenic src'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3
c radioactive decay
        k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'Loss of Be7 by decay'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3
c gravitational settling
        k = k + 1
        jls_grav(n) = k   ! special array grav. settling sinks
        sname_jls(k) = 'Grav_Settle_of_'//trname(n)
        lname_jls(k) = 'Loss of Be7 by grav settling'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3

        case ('Be10')
c cosmogenic source from file/same as Be7
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Cosmogenic_src_of'//trname(n)
        lname_jls(k) = 'Be10 cosmogenic src'
        jls_ltop(k) = lm
        jls_power(k) = -28  !may need changing around
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3
c gravitational settling
        k = k + 1
        jls_grav(n) = k   ! special array grav. settling sinks
        sname_jls(k) = 'Grav_Settle_of_'//trname(n)
        lname_jls(k) = 'Loss of Be10 by grav settling'
        jls_ltop(k) = lm
        jls_power(k) = -28  !may need changing around
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3

        case ('Pb210')
c source of Pb210 from Rn222 decay
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Radioactive_src_of'//trname(n)
        lname_jls(k) = 'Pb210 radioactive src'
        jls_ltop(k) = lm
        jls_power(k) =-26   ! -10  !may need to be changed
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3
c radioactive decay
        k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'Loss of Pb210 by decay'
        jls_ltop(k) = lm
        jls_power(k) =-26   ! -10  !may need to be changed
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3
c gravitational settling
        k = k + 1
        jls_grav(n) = k   ! special array grav. settling sinks
        sname_jls(k) = 'Grav_Settle_of_'//trname(n)
        lname_jls(k) = 'Loss of Pb210 by grav settling'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3

        case ('H2O2_s')
c gas phase source and sink of H2O2
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'gas_phase_source_of'//trname(n)
        lname_jls(k) = 'H2O2 gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'gas_phase_sink_of'//trname(n)
        lname_jls(k) = 'H2O2 gas phase sink'
        jls_ltop(k) = LM
        jls_power(k) = 2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
c convective chem cloud phase sink of H2O2
        k = k + 1
        jls_incloud(1,n) = k
        sname_jls(k) = 'mc_cloud_chem_sink_of'//trname(n)
        lname_jls(k) = 'H2O2 used in convective cloud chem'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c stratiform chem cloud phase sink of H2O2
        k = k + 1
        jls_incloud(2,n) = k
        sname_jls(k) = 'ss_cloud_chem_sink_of'//trname(n)
        lname_jls(k) = 'H2O2 used in stratiform cloud chem'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c photolysis rate
        k = k + 1
        jls_phot = k
        sname_jls(k) = 'photolysis_rate_of_H2O2'//trname(n)
        lname_jls(k) = 'photolysis rate of H2O2'
        jls_ltop(k) =LM
        jls_power(k) =-9
        units_jls(k) = unit_string(jls_power(k),'/s')
#endif
      case ('BCII')
        k = k + 1
        jls_3Dsource(1,n) = k   ! defined but not output
        sname_jls(k) = 'unused'   ! 'Aging_sink_of'//trname(n)
        lname_jls(k) = 'unused'   ! 'BCII aging sink'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of BCII
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of BCII'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Indust_src_'//trname(n)
        lname_jls(k) = 'BCII Industrial source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'Ship_src_'//trname(n)
        lname_jls(k) = 'BCII ship source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('BCIA')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Aging_source_of'//trname(n)
        lname_jls(k) = 'BCIA aging source'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Aircraft_source_of'//trname(n)
        lname_jls(k) = 'BCIA aircraft source'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of BCIA
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of BCIA'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('BCB', 'OCB')
c gravitational settling of BCB and OCB
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//trim(trname(n))
        lname_jls(k) = 'Gravitational Settling of '//trim(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c emission
        k = k + 1
        jls_3Dsource(nBiomass,n) = k
        sname_jls(k) = 'biomass_src_'//trim(trname(n))
        lname_jls(k) = trim(trname(n))//' Biomass source'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('OCII')
c gravitational settling of OCII
       k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//TRIM(trname(n))
        lname_jls(k) = 'Grav Settling of '//TRIM(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k   ! defined but not output
        sname_jls(k) = 'unused'    ! 'Aging_sink_'//trname(n)
        lname_jls(k) = 'unused'    ! 'OCIA aging sink'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'indust_src_'//TRIM(trname(n))
        lname_jls(k) = TRIM(trname(n))//' Industrial source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'Ship_src_'//TRIM(trname(n))
        lname_jls(k) = 'OCII ship source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifndef TRACERS_AEROSOLS_SOA
        k = k + 1
        jls_source(3,n) = k ! Terpene source should always be last
        sname_jls(k) = 'terpene_src_'//TRIM(trname(n))
        lname_jls(k) = 'OCII Terpene source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif  /* TRACERS_AEROSOLS_SOA */
      case ('OCI3')
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//TRIM(trname(n))
        lname_jls(k) = 'Grav Settling of '//TRIM(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k   ! 3D source
        sname_jls(k) = 'emission_of_'//TRIM(trname(n))
        lname_jls(k) = 'Emission of'//TRIM(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(2,n) = k   ! defined but not output
        sname_jls(k) = 'unused'    ! 'Aging_sink_'//trname(n)
        lname_jls(k) = 'unused'    ! 'OCIA aging sink'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifndef TRACERS_AEROSOLS_SOA
        k = k + 1
        jls_source(1,n) = k ! Terpene source should always be last
        sname_jls(k) = 'terpene_src_'//TRIM(trname(n))
        lname_jls(k) = 'OCI3 Terpene source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif  /* TRACERS_AEROSOLS_SOA */
      case ('OCI1','OCI2')
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//TRIM(trname(n))
        lname_jls(k) = 'Grav Settling of '//TRIM(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k   ! 3D source
        sname_jls(k) = 'emission_of_'//TRIM(trname(n))
        lname_jls(k) = 'Emission of'//TRIM(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(2,n) = k   ! defined but not output
        sname_jls(k) = 'unused'    ! 'Aging_sink_'//trname(n)
        lname_jls(k) = 'unused'    ! 'OCIA aging sink'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       case ('OCA4')
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//TRIM(trname(n))
        lname_jls(k) = 'Grav Settling of '//TRIM(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k   ! 3D source
        sname_jls(k) = 'emission_of_'//TRIM(trname(n))
        lname_jls(k) = 'Emission of'//TRIM(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('OCIA','OCA1','OCA2','OCA3')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Aging_source_'//TRIM(trname(n))
        lname_jls(k) = TRIM(trname(n))//' aging source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifdef FAA_emission
        select case(trname(n))
        case ('OCIA')
          k = k + 1
          jls_3Dsource(2,n) = k
          sname_jls(k) = 'Aircraft_source_of'//trname(n)
          lname_jls(k) = 'OCIA aircraft source'
          jls_ltop(k) = LM
          jls_power(k) = -1
          units_jls(k) = unit_string(jls_power(k),'kg/s')
        end select
#endif
c gravitational settling of OCIA,OCA1,OCA2,OCA3
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of_'//TRIM(trname(n))
        lname_jls(k) ='Grav Settling of '//TRIM(trname(n))
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('seasalt1')
        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Ocean_source_of'//trname(n)
        lname_jls(k) = 'seasalt1 ocean source'
        jls_ltop(k) = 1
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

c gravitational settling of ss1
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of seasalt1'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('seasalt2')
        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Ocean_source_of'//trname(n)
        lname_jls(k) = 'seasalt2 ocean source'
        jls_ltop(k) = 1
        jls_power(k) =1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

c gravitational settling of ss2
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of seasalt2'
        jls_ltop(k) = LM
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg/s')

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
        CASE('Clay','Silt1','Silt2','Silt3','Silt4',
     &       'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &       'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &       'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &       'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &       'Sil1QuHe','Sil2QuHe','Sil3QuHe')

        k=k+1
          jls_isrc(nDustEmjl,n)=k
          lname_jls(k)='Emission of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_emission'
          jls_ltop(k)=1
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg/s')
        IF (imDust == 0) THEN
          k=k+1
          jls_isrc(nDustEm2jl,n)=k
          lname_jls(k)='Cubic emission of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_emission2'
          jls_ltop(k)=1
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg/s')
        END IF
#ifndef TRACERS_DRYDEP
        k=k+1
          jls_isrc(nDustTurbjl,n)=k
          lname_jls(k)='Turbulent deposition of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_turb_depo'
          jls_ltop(k)=1
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg/s')
#endif
        k=k+1
          jls_grav(n)=k
          lname_jls(k)='Gain by gravitational settling of '
     &         //TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_grav_sett'
          jls_ltop(k)=Lm
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg/s')
#ifndef TRACERS_WATER
        k=k+1
          jls_wet(n)=k
          lname_jls(k)='Loss by wet deposition of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_wet_depo'
          jls_ltop(k)=Lm
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg/s')
#endif
#endif /* TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM */

C**** Here are some more examples of generalised diag. configuration
c      n = n_dust
c        k = k + 1
c        jls_grav(n) = k   ! special array grav. settling sinks
c        sname_jls(k) = 'Grav_Settle_of_'//trname(n)
c        lname_jls(k) = 'LOSS OF DUST BY SETTLING'
c        jls_ltop(k) = lm
c        jls_power(k) = -11
c        units_jls(k) = unit_string(jls_power(k),'kg/s')

      end select

#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
c**** additional wet deposition diagnostics
      IF (diag_wetdep == 1) THEN
        k=k+1
        jls_trdpmc(1,n)=k
        lname_jls(k)='MC Condensation of '//TRIM(trname(n))
        sname_jls(k)=TRIM(trname(n))//'_cond_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
        k=k+1
        jls_trdpmc(2,n)=k
        lname_jls(k)='Evaporated '//TRIM(trname(n))//' in MC Downdrafts'
        sname_jls(k)=TRIM(trname(n))//'_downeva_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
        k=k+1
        jls_trdpmc(3,n)=k
        lname_jls(k)='Condensed '//TRIM(trname(n))//' in MC CLW'
        sname_jls(k)=TRIM(trname(n))//'_conclw_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
        k=k+1
        jls_trdpmc(4,n)=k
        lname_jls(k)='Precipitated '//TRIM(trname(n))//' by MC'
        sname_jls(k)=TRIM(trname(n))//'_precip_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
        k=k+1
        jls_trdpmc(5,n)=k
        lname_jls(k)='Reevaporated '//TRIM(trname(n))//' from MC Precip'
        sname_jls(k)=TRIM(trname(n))//'_reevap_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
        k=k+1
        jls_trdpmc(6,n)=k
        lname_jls(k)='MC Washout of '//TRIM(trname(n))
        sname_jls(k)=TRIM(trname(n))//'_washout_mc'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
        k=k+1
        jls_trdpls(1,n)=k
        lname_jls(k)='LS Washout of '//TRIM(trname(n))
        sname_jls(k)=TRIM(trname(n))//'_washout_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
        k=k+1
        jls_trdpls(2,n)=k
        lname_jls(k)='Precipitated '//TRIM(trname(n))//' by LS'
        sname_jls(k)=TRIM(trname(n))//'_precip_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
        k=k+1
        jls_trdpls(3,n)=k
        lname_jls(k)='Condensed '//TRIM(trname(n))// ' in LS CLW'
        sname_jls(k)=TRIM(trname(n))//'_conclw_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
        k=k+1
        jls_trdpls(4,n)=k
        lname_jls(k)='Reevaporated '//TRIM(trname(n))//' from LS Precip'
        sname_jls(k)=TRIM(trname(n))//'_reevap_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
        k=k+1
        jls_trdpls(5,n)=k
        lname_jls(k)='Evaporated '//TRIM(trname(n))//' from LS CLW'
        sname_jls(k)=TRIM(trname(n))//'_clwevap_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
        k=k+1
        jls_trdpls(6,n)=k
        lname_jls(k)='LS Condensation of '//TRIM(trname(n))
        sname_jls(k)=TRIM(trname(n))//'_cond_ls'
        jls_ltop(k)=Lm
        jls_power(k)=1
        units_jls(k)=unit_string(jls_power(k),'kg/s')
      END IF
#endif

c
C**** Checks
      if (ntsurfsrc(n).gt.ntsurfsrcmax) then
!       write(6,*) ' ntsurfsrc too large for ',trname(n)
        if (am_i_root())
     &      write(6,*) ' Increase ntsurfsrcmax to at least',ntsurfsrc(n)
        call stop_model(
     &       ' Ntsurfsrc too large.  Increase ntsurfsrcmax',255)
      end if

      end do

C**** Additional Special JL diagnostics
C**** (not necessary associated with a particular tracer)
#ifdef TRACERS_SPECIAL_Shindell
#ifdef SHINDELL_STRAT_CHEM
        k = k + 1
        jls_ClOcon=k
        sname_jls(k) = 'ClO_conc'
        lname_jls(k) = 'ClO concentration'
        jls_ltop(k)  = LTOP
        jls_power(k) = -11
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'V/V air')
        k = k + 1
        jls_H2Ocon=k
        sname_jls(k) = 'H2O_conc'
        lname_jls(k) = 'H2O concentration'
        jls_ltop(k)  = LTOP
        jls_power(k) = -7
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'V/V air')
        k = k + 1
        jls_H2Ochem=k
        sname_jls(k) = 'H2O_chem'
        lname_jls(k) = 'H2O change due to chemistry'
        jls_ltop(k)  = LTOP
        jls_power(k) = -4
        scale_jls(k) = 1./DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
        k = k + 1
        jls_Oxp=k
        sname_jls(k) = 'Ox_chem_prod'
        lname_jls(k) = 'Ox production due to chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 2
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_Oxd=k
        sname_jls(k) = 'Ox_chem_dest'
        lname_jls(k) = 'Ox destruction due to chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 2
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_OxpT=k
        sname_jls(k) = 'trop_Ox_chem_prod'
        lname_jls(k) = 'Troposphere Ox prod by chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 2
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_OxdT=k
        sname_jls(k) = 'trop_Ox_chem_dest'
        lname_jls(k) = 'Troposphere Ox dest by chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 2
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_COp=k
        sname_jls(k) = 'CO_chem_prod'
        lname_jls(k) = 'CO production due to chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 1
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_COd=k
        sname_jls(k) = 'CO_chem_dest'
        lname_jls(k) = 'CO destruction due to chemistry'
        jls_ltop(k)  = LM
        jls_power(k) = 1
        scale_jls(k) = 1.d0/DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_OHcon=k
        sname_jls(k) = 'OH_conc'
        lname_jls(k) = 'OH concentration'
        jls_ltop(k)  = LTOP
        jls_power(k) = 5
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'molecules/cm3')
c
        k = k + 1
        jls_H2Omr=k
        sname_jls(k) = 'H2O_mr'
        lname_jls(k) = 'H2O mixing ratio (weighted by daylight)'
        jls_ltop(k)  = LTOP
        jls_power(k) = -4
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'parts/vol')
c
        k = k + 1
        jls_day=k
        sname_jls(k) = 'daylight'   ! not output
        lname_jls(k) = 'Daylight weighting'
        jls_ltop(k)  = 1
        jls_power(k) = 0
        scale_jls(k) = 100.
        units_jls(k) = unit_string(jls_power(k),'%')
c
        k = k + 1
        jls_N2O5sulf=k
        sname_jls(k) = 'N2O5_sulf'
        lname_jls(k) = 'N2O5 sulfate sink'
        jls_ltop(k)  = LTOP
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')

#endif  /* TRACERS_SPECIAL_Shindell */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
c Oxidants
#ifndef TRACERS_SPECIAL_Shindell
        k = k + 1
        jls_OHconk = k
        sname_jls(k) = 'OH_conc'
        lname_jls(k) = 'OH Concentration'
        jls_ltop(k) = LM
        jls_power(k) =5
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'molec/cm3')
#endif
        k = k + 1
        jls_HO2con = k
        sname_jls(k) = 'HO2_conc'
        lname_jls(k) = 'HO2 Concentration'
        jls_ltop(k) =LM
        jls_power(k) =7
        scale_jls(k) =1.
        units_jls(k) = unit_string(jls_power(k),'molec/cm3')

        k = k + 1
        jls_NO3 = k
        sname_jls(k) = 'NO3_conc'
        lname_jls(k) = 'NO3 Concentration'
        jls_ltop(k) =LM
        jls_power(k) =5
        scale_jls(k) =1.
        units_jls(k) = unit_string(jls_power(k),'molec/cm3')
#endif  /* TRACERS_AEROSOLS_Koch || TRACERS_AMP */

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      k = k + 1
      jls_spec(nDustEv1jl)=k
      lname_jls(k)='No. dust events'
      sname_jls(k)='no_dust_ev1'
      jls_ltop(k)=1
      scale_jls(k)=Sday/Dtsrc
      units_jls(k)='1/d'
      k = k + 1
      jls_spec(nDustEv2jl)=k
      lname_jls(k)='No. dust events above threshold wind'
      sname_jls(k)='no_dust_ev2'
      jls_ltop(k)=1
      scale_jls(k)=Sday/Dtsrc
      units_jls(k)='1/d'
      k = k + 1
      jls_spec(nDustWthjl)=k
      lname_jls(k)='Threshold velocity for dust emission'
      sname_jls(k)='wtrsh'
      jls_ltop(k)=1
      scale_jls(k)=1.
      units_jls(k)='m/s'
#endif

      if (k.gt. ktajls) then
        if (AM_I_ROOT()) write (6,*)
     &   'tjl_defs: Increase ktajls=',ktajls,' to at least ',k
         print *,'k should be ',k
        call stop_model('ktajls too small',255)
      end if
#endif /* TRACERS_ON */

      return
      end subroutine init_jls_diag

      subroutine init_ijts_diag
!@sum init_ijts_diag Initialise lat/lon tracer diags
!@auth Gavin Schmidt
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      USE CONSTANT, only: sday
      USE MODEL_COM, only: dtsrc
      USE TRACER_COM
      USE DIAG_COM
#ifdef TRACERS_ON
      USE TRDIAG_COM
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      use tracers_dust, only: nDustEmij, nDustEm2ij, nDustEv1ij
     &   ,nDustEv2ij, nDustWthij, imDust
#endif
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
      USE CLOUDS, ONLY : diag_wetdep
#endif
#endif /* TRACERS_ON */
#ifdef TRACERS_AMP
      USE AMP_AEROSOL, only: AMP_DIAG_FC
#endif
      implicit none
      integer k,n,n1,kr,ktaijs_out
      character*50 :: unit_string
      CHARACTER*17 :: cform

#ifdef TRACERS_ON
C**** Defaults for ijts (sources, sinks, etc.)
      ijts_fc(:,:)=0
      ijts_3Dsource(:,:)=0
      ijts_aq(:)=0
      ijts_isrc(:,:)=0
      ijts_gasex(:,:)=0
      ijts_HasArea(:) = .true. ! default applies to >50% of cases ???
      denom_ijts(:) = 0
#ifdef TRACERS_AMP
      ijts_AMPe(:)=0
      ijts_AMPp(:,:)=0
      ijts_AMPpdf(:,:)=0
#endif
C**** This needs to be 'hand coded' depending on circumstances
      k = 0

      do n=1,ntm
      select case (trname(n))

      case ('CFCn')
      k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CFCn Ocean source'
        sname_ijts(k) = 'CFCn_Ocean_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      k = k+1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Layer 1 SOURCE'
        sname_ijts(k) = trim(trname(n))//'_CFC-GRID_SOURCE_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      k = k+1  ! Gas Exchange Coefficient (piston velocity) (open ocean only)
        ijts_gasex(1,n)  = k
        ia_ijts(k) = ia_srf
        sname_ijts(k) = 'Piston_Veloc_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Piston Velocity'
        ijtc_power(n) = -5
        units_ijts(k) = unit_string(ijtc_power(n),'m/s')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.

      k = k+1  ! Gas Exchange Solubility coefficient
        ijts_gasex(2,n)  = k
        ia_ijts(k) = ia_srf
        sname_ijts(k) = 'Solubility_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Solubility'
        ijtc_power(n) = -5
        units_ijts(k) = unit_string(ijtc_power(n),'mol/m3/uatm')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.

      k = k+1  ! Gas exchange 
        ijts_gasex(3,n)  = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Gas_Exchange_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Gas Exchange'
        ijtc_power(n) = 0
        units_ijts(k) = unit_string(ijtc_power(n),'molCFC/m2/yr')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.

      case ('CO2n')
      k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'CO2_O_GASX'
        lname_ijts(k) = 'AO GASEX CO2'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg,CO2/m2/s')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      k = k+1  ! Gas Exchange Coefficient (piston velocity) (open ocean only)
        ijts_gasex(1,n)  = k
        ia_ijts(k) = ia_srf
        sname_ijts(k) = 'Piston_Veloc_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Piston Velocity'
        ijtc_power(n) = -5
        units_ijts(k) = unit_string(ijtc_power(n),'m/s')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.

      k = k+1  ! Gas Exchange Solubility coefficient
        ijts_gasex(2,n)  = k
        ia_ijts(k) = ia_srf
        sname_ijts(k) = 'Solubility_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Solubility'
        ijtc_power(n) = -5
        units_ijts(k) = unit_string(ijtc_power(n),'mol/m3/uatm')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.
      
      k = k+1  ! Gas exchange 
        ijts_gasex(3,n)  = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Gas_Exchange_'//trim(trname(n))
        dname_ijts(k) = 'ocnfr'
        lname_ijts(k) = trim(trname(n))//' Gas Exchange'
        ijtc_power(n) = 0
        units_ijts(k) = unit_string(ijtc_power(n),'molCO2/m2/yr')
        scale_ijts(k) = 10.**(-ijtc_power(n))
        ijts_HasArea(k) = .false.

      case ('SF6','SF6_c')
      k = k+1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Layer 1 SOURCE'
        sname_ijts(k) = trim(trname(n))//'_CFC-GRID_SOURCE_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('Rn222')
      k = k+1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Rn222 L 1 SOURCE'
        sname_ijts(k) = 'Rn222_SOURCE_Layer_1'
        ijts_power(k) = -21
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('CO2')
! keep AIJ and AJL CO2 sources in same order !!
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Fossil_fuel_source_'//trname(n)
        lname_ijts(k) = 'CO2 Fossil fuel src'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(2,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'fertilization_sink_'//trname(n)
        lname_ijts(k) = 'CO2 fertilization'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(3,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Northern_forest_regrowth_'//trname(n)
        lname_ijts(k) = 'CO2 North forest regrowth'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(4,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Land_Use_Modification_'//trname(n)
        lname_ijts(k) = 'CO2 from Land use mods'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(5,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Ecosystem_exchange_'//trname(n)
        lname_ijts(k) = 'CO2 Ecosystem exch'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(6,n) = k
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Ocean_exchange_'//trname(n)
        lname_ijts(k) = 'CO2 Ocean exchange'
        ijts_power(k) = -11
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('N2O')
#ifdef TRACERS_SPECIAL_Shindell
        do kr=1,ntsurfsrc(n)
          k = k+1
          ijts_source(kr,n) = k
          ia_ijts(k) = ia_src
          sname_ijts(k) = trim(trname(n))//'_'//trim(ssname(n,kr))
          lname_ijts(k) = trname(n)//' source from '//trim(ssname(n,kr))
          ijts_power(k) = -14
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        end do
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' Chemistry'
        sname_ijts(k) = trim(trname(n))//'_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(nOverwrite,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' Overwrite'
        sname_ijts(k) = trim(trname(n))//'_overw'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
#ifdef TRACERS_SPECIAL_Lerner
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'N2O CHANGE IN L 1'
        sname_ijts(k) = 'N2O_CHANGE_IN_L_1'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' Chemistry'
        sname_ijts(k) = trim(trname(n))//'_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif

      case ('CFC11')
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CFC_11 L 1 SOURCE'
        sname_ijts(k) = 'CFC_11_SOURCE_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CFC_11 Stratospheric Chemistry Sink'
        sname_ijts(k) = 'CFC_11_strat_sink'
        ijts_power(k) = -18
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('14CO2')
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = '14CO2 L 1 Sink'
        sname_ijts(k) = '14CO2_L1_Sink'
        ijts_power(k) = -21
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('codirect')
        do kr=1,ntsurfsrc(n)
          k = k+1
          ijts_source(kr,n) = k
          ia_ijts(k) = ia_src
          sname_ijts(k) = trim(trname(n))//'_'//trim(ssname(n,kr))
          lname_ijts(k) = trname(n)//' source from '//trim(ssname(n,kr))
          ijts_power(k) = -14
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        end do

      case ('NOx','CO','Isoprene','Alkenes','Paraffin',
     &'ClOx','BrOx','HCl','HOCl','ClONO2','HBr','HOBr','BrONO2',
     &'CFC','H2O2','CH3OOH','Ox','N2O5','HNO3','HCHO','Terpenes',
     &'HO2NO2','PAN','AlkylNit','stratOx')
        do kr=1,ntsurfsrc(n)
          k = k+1
          ijts_source(kr,n) = k
          ia_ijts(k) = ia_src
          sname_ijts(k) = trim(trname(n))//'_'//trim(ssname(n,kr))
          lname_ijts(k) = trim(trname(n))//' source from '//
     &                    trim(ssname(n,kr))
          ijts_power(k) = -14
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        end do
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' Chemistry'
        sname_ijts(k) = trim(trname(n))//'_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        select case(trname(n))
        case('NOx','CO','Isoprene','Alkenes','Paraffin',
     &  'CFC','H2O2','CH3OOH','Ox','N2O5','HNO3','HCHO',
     &  'Terpenes','HO2NO2','PAN','AlkylNit','stratOx')
          k = k + 1
          ijts_3Dsource(nOverwrite,n) = k
          ia_ijts(k) = ia_src
          lname_ijts(k) = trname(n)//' Overwrite'
          sname_ijts(k) = trim(trname(n))//'_overw'
          ijts_power(k) = -12
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        end select
        select case(trname(n))
        case('NOx')
          k = k + 1
          ijts_3Dsource(nOther,n) = k
          ia_ijts(k) = ia_src
          lname_ijts(k) = trname(n)//' Lightning Source'
          sname_ijts(k) = trim(trname(n))//'_lightning'
          ijts_power(k) = -12
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
          k = k + 1
          ijts_3Dsource(nAircraft,n) = k
          ia_ijts(k) = ia_src
          lname_ijts(k) = trname(n)//' Aircraft Source'
          sname_ijts(k) = trim(trname(n))//'_aircraft'
          ijts_power(k) = -12
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef FAA_emission
       case('CO','Alkenes','Paraffin')
          k = k + 1
          ijts_3Dsource(nAircraft,n) = k
          ia_ijts(k) = ia_src
          lname_ijts(k) = trname(n)//' Aircraft Source'
          sname_ijts(k) = trim(trname(n))//'_aircraft'
          ijts_power(k) = -12
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
        case('Ox','stratOx')
          k = k + 1
          ijts_fc(1,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trname(n)//' tropopause SW rad forc'
          sname_ijts(k) = 'swf_tp_'//trim(trname(n))
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
          k = k + 1
          ijts_fc(2,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trname(n)//' tropopause LW rad forc'
          sname_ijts(k) = 'lwf_tp_'//trim(trname(n))
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
          k = k + 1
          ijts_fc(3,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trname(n)//' TOA SW rad forc'
          sname_ijts(k) = 'swf_toa_'//trim(trname(n))
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
          k = k + 1
          ijts_fc(4,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trname(n)//' TOA LW rad forc'
          sname_ijts(k) = 'lwf_toa_'//trim(trname(n))
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
#ifdef AUX_OX_RADF_TROP
#ifndef AUXILIARY_OX_RADF
          call stop_model('AUX_OX_RADF_TROP needs AUXILIARY_OX_RADF',13)
#endif
#endif
#ifdef AUXILIARY_OX_RADF
          if(trname(n)=='Ox')then
            k = k + 1
            ijts_auxfc(1) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//' AUX tropp SW rad forc'
            sname_ijts(k) = 'swfauxtp_'//trim(trname(n))
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
            k = k + 1
            ijts_auxfc(2) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//' AUX tropp LW rad forc'
            sname_ijts(k) = 'lwfauxtp_'//trim(trname(n))
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
            k = k + 1
            ijts_auxfc(3) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//' AUX TOA SW rad forc'
            sname_ijts(k) = 'swfauxtoa_'//trim(trname(n))
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
            k = k + 1
            ijts_auxfc(4) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//' AUX TOA LW rad forc'
            sname_ijts(k) = 'lwfauxtoa_'//trim(trname(n))
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
          endif
#endif /* AUXILIARY_OX_RADF */
#ifdef ACCMIP_LIKE_DIAGS
          if(trname(n)=='Ox' .and. dodrydep(n))then
            k = k+1
            ijts_Sdrydep = k
            ia_ijts(k) = ia_src
            lname_ijts(k) = trim(trname(n))//' stomatal drydep flux'
            sname_ijts(k) = 'stomatal_'//trim(trname(n))
            ijts_power(k) = ntm_power(n)-4
            units_ijts(k) = unit_string(ijts_power(k),'kg/m^2/s')
            scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
            ijts_HasArea(k) = .false.
          end if
#endif /* ACCMIP_LIKE_DIAGS */
        end select
        select case(trname(n))
        case('NOx','CO','Alkenes','Paraffin')
          k = k + 1
          ijts_3Dsource(nBiomass,n) = k
          ia_ijts(k) = ia_src
          lname_ijts(k) = trim(trname(n))//' Biomass source'
          sname_ijts(k) = trim(trname(n))//'_Biomass_source'
          ijts_power(k) = -12
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        end select

#ifdef TRACERS_AEROSOLS_SOA
      case ('isopp1g','isopp1a','isopp2g','isopp2a',
     &      'apinp1g','apinp1a','apinp2g','apinp2a')
c chemical production
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Chemistry'
        sname_ijts(k) = trim(trname(n))//'_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        select case(trname(n))
        case('isopp1a')
          ! In the radiation code the RCOMPX call for isopp1a
          ! currently contains isopp1a+isopp2a and if TRACERS_TERP
          ! also apinp1a+apinp2a, so using SOA in stead of trname:
           
          call set_diag_rad(n,k)

c SOA shortwave radiative forcing
          k = k + 1
          ijts_fc(1,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = 'SOA SW radiative forcing'
          sname_ijts(k) = 'swf_SOA'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c SOA longwave radiative forcing
          k = k + 1
          ijts_fc(2,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = 'SOA LW radiative forcing'
          sname_ijts(k) = 'lwf_SOA'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c SOA shortwave surface radiative forcing
          k = k + 1
          ijts_fc(3,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = 'SOA SW surface rad forcing'
          sname_ijts(k) = 'swf_surf_SOA'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c SOA longwave surface radiative forcing
          k = k + 1
          ijts_fc(4,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = 'SOA LW surface rad forcing'
          sname_ijts(k) = 'lwf_surf_SOA'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c SOA clear sky shortwave radiative forcing
          k = k + 1
          ijts_fc(5,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = 'SOA clr sky SW rad forcing'
          sname_ijts(k) = 'swf_CS_SOA'
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c SOA clear sky longwave radiative forcing
          k = k + 1
          ijts_fc(6,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = 'SOA clr sky LW rad forcing'
          sname_ijts(k) = 'lwf_CS_SOA'
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c SOA clear sky shortwave surface radiative forcing
          k = k + 1
          ijts_fc(7,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = 'SOA clr sky SW surface rad forcing'
          sname_ijts(k) = 'swf_CS_surf_SOA'
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c SOA clear sky longwave surface radiative forcing
          k = k + 1
          ijts_fc(8,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = 'SOA clr sky LW surface rad forcing'
          sname_ijts(k) = 'lwf_CS_surf_SOA'
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
        end select ! isopp1a representing SOA as a group
#endif  /* TRACERS_AEROSOLS_SOA*/

      case ('CH4')
#ifdef TRACERS_SPECIAL_Shindell
        do kr=1,ntsurfsrc(n)
          k = k+1
          ijts_source(kr,n) = k
          ia_ijts(k) = ia_src
          sname_ijts(k) = trim(trname(n))//'_'//trim(ssname(n,kr))
          lname_ijts(k) = trname(n)//' source from '//trim(ssname(n,kr))
          ijts_power(k) = -14
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        end do
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' Chemistry'
        sname_ijts(k) = trim(trname(n))//'_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(nOverwrite,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' Overwrite'
        sname_ijts(k) = trim(trname(n))//'_overw'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(nBiomass,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Biomass source'
        sname_ijts(k) = trim(trname(n))//'_Biomass_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#else
      k = k + 1
        ijts_source(6,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 sink due to soil absorp.'
        sname_ijts(k) = 'CH4_soil_sink.'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(7,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Termite source'
        sname_ijts(k) = 'CH4_Termite_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(9,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Ocean source'
        sname_ijts(k) = 'CH4_Ocean_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(10,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Fresh Water lake source'
        sname_ijts(k) = 'CH4_lake_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(11,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Misc. Ground source'
        sname_ijts(k) = 'CH4_Misc._Ground_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(14,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Wetlands+Tundra source'
        sname_ijts(k) = 'CH4_Wetlands+Tundra_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Animal source'
        sname_ijts(k) = 'CH4_Animal_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Coal Mine source'
        sname_ijts(k) = 'CH4_Coal_Mine_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(3,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Gas Leak source'
        sname_ijts(k) = 'CH4_Gas_Leak_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(4,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Gas Venting source'
        sname_ijts(k) = 'CH4_Gas_Venting_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(5,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Municipal solid waste src'
        sname_ijts(k) = 'CH4_MSW_src'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(8,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Coal combustion source'
        sname_ijts(k) = 'CH4_Coal_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(12,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Biomass burning source'
        sname_ijts(k) = 'CH4_Biomass_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(13,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Rice cultivation source'
        sname_ijts(k) = 'CH4_Rice_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Tropospheric Chemistry'
        sname_ijts(k) = 'CH4_trop_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Stratospheric Chemistry'
        sname_ijts(k) = 'CH4_strat_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif

      case ('O3')
      k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'O3 deposition, layer 1'
        sname_ijts(k) = 'O3_deposition_L1'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'O3 Stratospheric Chem.'
        sname_ijts(k) = 'O3_strat_chem'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'O3 Tropo. Chem. Production'
        sname_ijts(k) = 'O3_trop_chem_prod'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(3,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'O3 Tropo. Chemistry Loss'
        sname_ijts(k) = 'O3_trop_chem_loss'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

#ifdef TRACERS_WATER
      case ('Water', 'H2O18', 'H2O17', 'HDO', 'HTO' )
          ! nothing I can think of....

#ifdef FAA_emission
      if (trname(n).eq."Water") then
          k = k + 1
          ijts_3Dsource(1,n) = k
          ia_ijts(k) = ia_src
          lname_ijts(k) = trname(n)//' Aircraft Source'
          sname_ijts(k) = trim(trname(n))//'_aircraft'
          ijts_power(k) = -12
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
         endif
#endif
#endif

#ifdef SHINDELL_STRAT_EXTRA
      case ('GLT')
      k = k+1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//'L1 overwrite source'
        sname_ijts(k) = trim(trname(n))//'L1_overwrite'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif

      case ('BCII')
        k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'BC Industrial source'
        sname_ijts(k) = 'BC_Industrial_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_source(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'BC ship source'
        sname_ijts(k) = 'BC_Ship_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(1,n) = k   ! defined but not output
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'unused'    ! 'BCII Aging sink'
        sname_ijts(k) = 'unused'    ! 'BCII_Aging_sink'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('BCIA')
        k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'BC Aging source'
        sname_ijts(k) = 'BC_Aging_Source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        k = k + 1
        ijts_3Dsource(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'BCIA Aircraft source'
        sname_ijts(k) = 'BCIA_Aircraft_src'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

#ifdef TRACERS_AEROSOLS_Koch

        call set_diag_rad(n,k)

c BCI shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'BCI SW radiative forcing'
        sname_ijts(k) = 'swf_BCI'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCI longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'BCI LW radiative forcing'
        sname_ijts(k) = 'lwf_BCI'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCI shortwave surface radiative forcing
        k = k + 1
        ijts_fc(3,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'BCI SW surface rad forcing'
        sname_ijts(k) = 'swf_surf_BCI'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCI longwave surface radiative forcing
        k = k + 1
        ijts_fc(4,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'BCI LW surface rad forcing'
        sname_ijts(k) = 'lwf_surf_BCI'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCI clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'BCI clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS_BCI'
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCI clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'BCI clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS_BCI'
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCI clear sky shortwave surface radiative forcing
        k = k + 1
        ijts_fc(7,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'BCI clr sky SW surface rad forcing'
        sname_ijts(k) = 'swf_CS_surf_BCI'
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCI clear sky longwave surface radiative forcing
        k = k + 1
        ijts_fc(8,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'BCI clr sky LW surface rad forcing'
        sname_ijts(k) = 'lwf_CS_surf_BCI'
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
#endif

      case ('BCB')
        k = k + 1
        ijts_3Dsource(nBiomass,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Biomass source'
        sname_ijts(k) = trim(trname(n))//'_Biomass_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_AEROSOLS_Koch

        call set_diag_rad(n,k)

c BCB shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' SW radiative forcing'
        sname_ijts(k) = 'swf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCB longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' LW radiative forcing'
        sname_ijts(k) = 'lwf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCB shortwave surface radiative forcing
        k = k + 1
        ijts_fc(3,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' SW surface rad forcing'
        sname_ijts(k) = 'swf_surf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCB longwave surface radiative forcing
        k = k + 1
        ijts_fc(4,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' LW surface rad forcing'
        sname_ijts(k) = 'lwf_surf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCB clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCB clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCB clear sky shortwave surface radiative forcing
        k = k + 1
        ijts_fc(7,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n)) //
     &       ' clr sky SW surface rad forcing'
        sname_ijts(k) = 'swf_CS_surf_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c BCB clear sky longwave surface radiative forcing
        k = k + 1
        ijts_fc(8,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n)) //
     &       ' clr sky LW surface rad forcing'
        sname_ijts(k) = 'lwf_CS_surf_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
#endif

      case ('OCII')
        k = k + 1
        ijts_source(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'OC Industrial source'
        sname_ijts(k) = 'OC_Industrial_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_source(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'OC ship source'
        sname_ijts(k) = 'OC_Ship_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifndef TRACERS_AEROSOLS_SOA
        k = k + 1
        ijts_source(3,n) = k ! Terpene source should always be last
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'OC Terpene source'
        sname_ijts(k) = 'OC_Terpene_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif  /* TRACERS_AEROSOLS_SOA */
        k = k + 1
        ijts_3Dsource(1,n) = k   ! defined but not output
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'unused'       ! 'OCII Aging sink'
        sname_ijts(k) = 'unused'       ! 'OCII_Aging_Sink'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('OCI3')
#ifndef TRACERS_AEROSOLS_SOA
        k = k + 1
        ijts_source(1,n) = k ! Terpene source should always be last
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'OC Terpene source'
        sname_ijts(k) = 'OC_Terpene_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif  /* TRACERS_AEROSOLS_SOA */
        k = k + 1
        ijts_3Dsource(1,n) = k   ! Emissions
        ia_ijts(k) = ia_src
        lname_ijts(k) =  TRIM(trname(n))//' Emissions'
        sname_ijts(k) = TRIM(trname(n))//'_emissions'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(2,n) = k   ! defined but not output
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'unused'       ! 'OCII Aging sink'
        sname_ijts(k) = 'unused'       ! 'OCII_Aging_Sink'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      case ('OCI1','OCI2')
        k = k + 1
        ijts_3Dsource(1,n) = k   ! Emissions
        ia_ijts(k) = ia_src
        lname_ijts(k) =  TRIM(trname(n))//' Emissions'
        sname_ijts(k) = TRIM(trname(n))//'_emissions'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(2,n) = k   ! defined but not output
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'unused'       ! 'Aging sink'
        sname_ijts(k) = 'unused'       ! 'OCII_Aging_Sink'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      case ('OCA4')
        k = k + 1
        ijts_3Dsource(1,n) = k   ! Emissions
        ia_ijts(k) = ia_src
        lname_ijts(k) =  TRIM(trname(n))//' Emissions'
        sname_ijts(k) = TRIM(trname(n))//'_emissions'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_OM_SP

        call set_diag_rad(n,k)

c OC shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' SW rad forcing'
        sname_ijts(k) = 'swf_'//TRIM(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OC longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' LW rad forcing'
        sname_ijts(k) = 'lwf_'//TRIM(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OC shortwave surface radiative forcing
        k = k + 1
        ijts_fc(3,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' SW surf rad forcing'
        sname_ijts(k) = 'swf_surf_'//TRIM(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OC longwave surface radiative forcing
        k = k + 1
        ijts_fc(4,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' LW surf rad forcing'
        sname_ijts(k) = 'lwf_surf_'//TRIM(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OC clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS_OC'
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OC clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS_'//TRIM(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OC clear sky shortwave surface radiative forcing
        k = k + 1
        ijts_fc(7,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' clr sky SW surf rad forcing'
        sname_ijts(k) = 'swf_CS_surf_'//TRIM(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OC clear sky longwave surface radiative forcing
        k = k + 1
        ijts_fc(8,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' clr sky LW surf rad forcing'
        sname_ijts(k) = 'lwf_CS_surf_'//TRIM(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
#endif
      case ('OCIA','OCA1','OCA2','OCA3')
        k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = TRIM(trname(n))//' Aging source'
        sname_ijts(k) = TRIM(trname(n))//'_Aging_Source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef FAA_emission
        select case (trname(n))
        case ('OCIA')
          k = k + 1
          ijts_3Dsource(2,n) = k
          ia_ijts(k) = ia_src
          lname_ijts(k) = 'OCIA Aircraft source'
          sname_ijts(k) = 'OCIA_Aircraft_src'
          ijts_power(k) = -12
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        end select
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP)

        call set_diag_rad(n,k)

c OCI shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' SW rad forcing'
        sname_ijts(k) = 'swf_'//TRIM(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCI longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' LW rad forcing'
        sname_ijts(k) = 'lwf_'//TRIM(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCI shortwave surface radiative forcing
        k = k + 1
        ijts_fc(3,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' SW surf rad forcing'
        sname_ijts(k) = 'swf_surf_'//TRIM(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCI longwave surface radiative forcing
        k = k + 1
        ijts_fc(4,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' LW surf rad forcing'
        sname_ijts(k) = 'lwf_surf_'//TRIM(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCI clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS_'//TRIM(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCI clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS_'//TRIM(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCI clear sky shortwave surface radiative forcing
        k = k + 1
        ijts_fc(7,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' clr sky SW surf rad forcing'
        sname_ijts(k) = 'swf_CS_surf_'//TRIM(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCI clear sky longwave surface radiative forcing
        k = k + 1
        ijts_fc(8,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' clr sky LW surf rad forcing'
        sname_ijts(k) = 'lwf_CS_surf_'//TRIM(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
#endif
      case ('OCB')
        k = k + 1
        ijts_3Dsource(nBiomass,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = trim(trname(n))//' Biomass source'
        sname_ijts(k) = trim(trname(n))//'_Biomass_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

#ifdef TRACERS_AEROSOLS_Koch

        call set_diag_rad(n,k)

c OCB shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' SW radiative forcing'
        sname_ijts(k) = 'swf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCB longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' LW radiative forcing'
        sname_ijts(k) = 'lwf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCB shortwave surface radiative forcing
        k = k + 1
        ijts_fc(3,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' SW surface rad forcing'
        sname_ijts(k) = 'swf_surf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCB longwave surface radiative forcing
        k = k + 1
        ijts_fc(4,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' LW surface rad forcing'
        sname_ijts(k) = 'lwf_surf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCB clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCB clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n))//' clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCB clear sky shortwave surface radiative forcing
        k = k + 1
        ijts_fc(7,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n)) //
     &       ' clr sky SW surface rad forcing'
        sname_ijts(k) = 'swf_CS_surf_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c OCB clear sky longwave surface radiative forcing
        k = k + 1
        ijts_fc(8,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = trim(trname(n)) //
     &       ' clr sky LW surface rad forcing'
        sname_ijts(k) = 'lwf_CS_surf_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
#endif


      case ('DMS')
        k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'DMS Ocean source'
        sname_ijts(k) = 'DMS_Ocean_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'DMS Chem sink'
        sname_ijts(k) = 'DMS_Chem_sink'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('MSA')
c put in chemical production of MSA
        k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'MSA Chemical source'
        sname_ijts(k) = 'MSA_Chemical_source'
        ijts_power(k) = -17
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('SO2')
c production of SO2 from volcanic emissions
        k = k + 1
        ijts_3Dsource(nVolcanic,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO2 source from volcanos'
        sname_ijts(k) = 'SO2_source_from_volcanos'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c production of SO2 from aircraft
        k = k + 1
        ijts_3Dsource(nAircraft,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO2 source from aircraft'
        sname_ijts(k) = 'SO2_source_from_aircraft'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c emissions of biomass SO2
        k = k + 1
        ijts_3Dsource(nBiomass,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Biomass SO2 source'
        sname_ijts(k) = 'SO2_source_from_biomass'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c put in production of SO2 from DMS
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO2 source from DMS'
        sname_ijts(k) = 'SO2_source_from_DMS'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c put in chemical loss of SO2
        k = k + 1
        ijts_3Dsource(nChemloss,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO2 Chemical sink'
        sname_ijts(k) = 'SO2_chem_sink'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c emissions of industrial SO2
        do kr=1,ntsurfsrc(n)
          k = k + 1
          ijts_source(kr,n) = k
          ia_ijts(k) = ia_src
          sname_ijts(k) = trim(trname(n))//'_src_'//trim(ssname(n,kr))
          lname_ijts(k) = trim(trname(n))//' source from '//
     &                   trim(ssname(n,kr))
          ijts_power(k) = -15
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        enddo

        case ('SO4')
c put in production of SO4 from gas phase
        k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4 gas phase source'
        sname_ijts(k) = 'SO4_gas_phase_source'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(nVolcanic,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4 volcanic source'
        sname_ijts(k) = 'SO4_volcanic_source'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef FAA_emission
c production of SO4 from aircraft
        k = k + 1
        ijts_3Dsource(nAircraft,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4 source from aircraft'
        sname_ijts(k) = 'SO4_source_from_aircraft'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
        k = k + 1
        ijts_3Dsource(nBiomass,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4 biomass source'
        sname_ijts(k) = 'SO4_biomass_source'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c SO4 from industrial emissions
        do kr=1,ntsurfsrc(n_SO2)
          k = k + 1
          ijts_source(kr,n) = k
          ia_ijts(k) = ia_src
          sname_ijts(k) = trim(trname(n))//'_src_'//
     &                    trim(ssname(n_SO2,kr))
          lname_ijts(k) = trim(trname(n))//' source from '//
     &                    trim(ssname(n_SO2,kr))
          ijts_power(k) = -15
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        enddo
#ifdef TRACERS_AEROSOLS_Koch
c put in source of SO4 from aqueous chem
        k = k + 1
        ijts_aq(n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4 aqueous chem source'
        sname_ijts(k) = 'SO4_aq_chem_source'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        call set_diag_rad(n,k)

c SO4 shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SO4 SW radiative forcing'
        sname_ijts(k) = 'swf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SO4 longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SO4 LW radiative forcing'
        sname_ijts(k) = 'lwf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SO4 shortwave surface radiative forcing
        k = k + 1
        ijts_fc(3,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SO4 SW surface rad forcing'
        sname_ijts(k) = 'swf_surf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SO4 longwave surface radiative forcing
        k = k + 1
        ijts_fc(4,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SO4 LW surface rad forcing'
        sname_ijts(k) = 'lwf_surf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SO4 clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SO4 clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SO4 clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SO4 clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SO4 clear sky shortwave surface radiative forcing
        k = k + 1
        ijts_fc(7,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SO4 clr sky SW surface rad forcing'
        sname_ijts(k) = 'swf_CS_surf_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SO4 clear sky longwave surface radiative forcing
        k = k + 1
        ijts_fc(8,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SO4 clr sky LW surface rad forcing'
        sname_ijts(k) = 'lwf_CS_surf_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
#endif

c#ifdef TRACERS_NITRATE
      case ('NH3')
c emissions of biomass NH3
        k = k + 1
        ijts_3Dsource(nBiomass,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Biomass NH3 source'
        sname_ijts(k) = 'NH3_source_from_biomass'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c emissions of industrial NH3
        do kr=1,ntsurfsrc(n)
          k = k + 1
          ijts_source(kr,n) = k
          ia_ijts(k) = ia_src
          sname_ijts(k) = trim(trname(n))//'_src_'//trim(ssname(n,kr))
          lname_ijts(k) = trim(trname(n))//' source from '//
     &                   trim(ssname(n,kr))
          ijts_power(k) = -15
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        enddo

      case ('NO3p')

        call set_diag_rad(n,k)

c NO3 shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'NO3 SW radiative forcing'
        sname_ijts(k) = 'swf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c NO3 longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'NO3 LW radiative forcing'
        sname_ijts(k) = 'lwf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c NO3 shortwave surface radiative forcing
        k = k + 1
        ijts_fc(3,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'NO3 SW surface rad forcing'
        sname_ijts(k) = 'swf_surf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c NO3 longwave surface radiative forcing
        k = k + 1
        ijts_fc(4,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'NO3 LW surface rad forcing'
        sname_ijts(k) = 'lwf_surf_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c NO3 clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'NO3 clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c NO3 clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'NO3 clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c NO3 clear sky shortwave surface radiative forcing
        k = k + 1
        ijts_fc(7,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'NO3 clr sky SW surface rad forcing'
        sname_ijts(k) = 'swf_CS_surf_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c NO3 clear sky longwave surface radiative forcing
        k = k + 1
        ijts_fc(8,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'NO3 clr sky LW surface rad forcing'
        sname_ijts(k) = 'lwf_CS_surf_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c#endif
#ifdef TRACERS_AMP
        case ('M_NO3   ','M_NH4   ','M_H2O   ','M_AKK_SU','N_AKK_1 ',!AKK
     *    'M_ACC_SU','N_ACC_1 ','M_DD1_SU','M_DD1_DU','N_DD1_1 ',!ACC,DD1
     *    'M_DS1_SU','M_DS1_DU','N_DS1_1 ','M_DD2_SU','M_DD2_DU',!DS1,DD2
     *    'N_DD2_1 ','M_DS2_SU','M_DS2_DU','N_DS2_1 ','M_SSA_SU',!DD2,DS2,SSA
     *    'M_SSA_SS','M_SSC_SS'                                 ,!SSA,SSC
     *    'M_OCC_SU','M_OCC_OC','N_OCC_1 ','M_BC1_SU','M_BC1_BC',!OCC,BC1
     *    'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ','M_BC3_SU',!BC1,BC2,BC3
     *    'M_BC3_BC','N_BC3_1 ','M_DBC_SU','M_DBC_BC','M_DBC_DU',!BC3,DBC
     *    'N_DBC_1 ','M_BOC_SU','M_BOC_BC','M_BOC_OC','N_BOC_1 ',!DBC,BOC
     *    'M_BCS_SU','M_BCS_BC','N_BCS_1 ','M_MXX_SU','M_MXX_BC',!BCS,MXX
     *    'M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ','M_OCS_SU',
     *    'M_OCS_OC','N_OCS_1 ','M_SSS_SS','M_SSS_SU')
       k = k + 1
         ijts_3Dsource(nChemistry,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'AMP_src_'//trim(trname(n))
         sname_ijts(k) = 'AMP_src_'//trim(trname(n))
         ijts_power(k) = -15
         units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
         scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
       k = k + 1
         ijts_AMPp(1,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P1_Nucl_'//trim(trname(n))
         sname_ijts(k) = 'P1_Nucl_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(2,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P2_Coag_'//trim(trname(n))
         sname_ijts(k) = 'P2_Coag_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(3,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P3_Cond_'//trim(trname(n))
         sname_ijts(k) = 'P3_Cond_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(4,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P4_Incld_NIMC_'//trim(trname(n))
         sname_ijts(k) = 'P4_Incld_NIMC_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(5,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P5_IMLoss_NIAC_'//trim(trname(n))
         sname_ijts(k) = 'P5_IMLoss_NIAC_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(6,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P6_Mode_Trans_'//trim(trname(n))
         sname_ijts(k) = 'P6_Mode_Trans_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
       k = k + 1
         ijts_AMPp(7,n)=k
         ia_ijts(k) = ia_src
         lname_ijts(k) = 'P7_Total_Change_'//trim(trname(n))
         sname_ijts(k) = 'P7_Total_Change_'//trim(trname(n))
         ijts_power(k) = -11
         units_ijts(k) = unit_string(ijts_power(k),' ')
         scale_ijts(k) = 10.**(-ijts_power(k))
#endif
#ifdef TRACERS_HETCHEM
      case ('SO4_d1')
c chemical production of SO4 from SO2 on dust
        k = k + 1
        ijts_source(1,n) = k  ! 3dsource?
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d1 Chemical source'
        sname_ijts(k) = 'SO4d1_Chemical_source'
        ijts_power(k) = -10
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      case ('SO4_d2')
c chemical production of SO4 from SO2 on dust
        k = k + 1
        ijts_source(1,n) = k  ! 3dsource?
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d2 Chemical source'
        sname_ijts(k) = 'SO4d2_Chemical_source'
        ijts_power(k) = -10
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      case ('SO4_d3')
c chemical production of SO4 from SO2 on dust
        k = k + 1
        ijts_source(1,n) = k  ! 3dsource?
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d3 Chemical source'
        sname_ijts(k) = 'SO4d3_Chemical_source'
        ijts_power(k) = -10
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
      case ('H2O2_s')
c put in production of H2O2 from gas phase
        k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'H2O2 gas phase source'
        sname_ijts(k) = 'H2O2_gas_phase_source'
        ijts_power(k) = -10
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c put in production of H2O2 from gas phase
        k = k + 1
        ijts_3Dsource(2,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'H2O2 gas phase sink'
        sname_ijts(k) = 'H2O2_gas_phase_sink'
        ijts_power(k) = -10
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('Be7','Be10')
c cosmogenic source from file
        k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Cosmogenic source of '//trname(n)
        sname_ijts(k) = trim(trname(n))//'_cosmo_src'
        ijts_power(k) = -25
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('Pb210')
c source of Pb210 from Rn222 decay
        k = k + 1
        ijts_3Dsource(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Radioactive source of '//trname(n)
        sname_ijts(k) = 'Pb210_radio_src'
        ijts_power(k) = -24
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('seasalt1')
        k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'seasalt1 Ocean source'
        sname_ijts(k) = 'seasalt1_Ocean_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_AEROSOLS_Koch
        call set_diag_rad(n,k)

c SS shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SS SW radiative forcing'
        sname_ijts(k) = 'swf_SS'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SS longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SS LW radiative forcing'
        sname_ijts(k) = 'lwf_SS'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SS shortwave surface radiative forcing
        k = k + 1
        ijts_fc(3,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SS SW surface rad forcing'
        sname_ijts(k) = 'swf_surf_SS'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SS longwave surface radiative forcing
        k = k + 1
        ijts_fc(4,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SS LW surface rad forcing'
        sname_ijts(k) = 'lwf_surf_SS'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SS clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SS clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS_SS'
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SS clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SS clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS_SS'
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SS clear sky shortwave surface radiative forcing
        k = k + 1
        ijts_fc(7,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SS clr sky SW surface rad forcing'
        sname_ijts(k) = 'swf_CS_surf_SS'
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c SS clear sky longwave surface radiative forcing
        k = k + 1
        ijts_fc(8,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'SS clr sky LW surface rad forcing'
        sname_ijts(k) = 'lwf_CS_surf_SS'
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
#endif
       case ('seasalt2')
        k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'seasalt2 Ocean source'
        sname_ijts(k) = 'seasalt2_Ocean_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_AEROSOLS_Koch

       call set_diag_rad(n,k)

#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      CASE('Clay','Silt1','Silt2','Silt3','Silt4',
     &   'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &   'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &   'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &   'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &   'Sil1QuHe','Sil2QuHe','Sil3QuHe')
        k=k+1
        ijts_isrc(nDustEmij,n)=k
        lname_ijts(k)='Emission of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_emission'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        IF (imDust == 0) THEN
        k=k+1
        ijts_isrc(nDustEm2ij,n)=k
        lname_ijts(k)='Cubic emission of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_emission2'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        END IF
#ifndef TRACERS_DRYDEP
      k=k+1
        ijts_isrc(nDustTurbij,n)=k
        lname_ijts(k)='Turbulent Deposition of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_turb_depo'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
#ifndef TRACERS_WATER
      k=k+1
        ijts_wet(n)=k
        lname_ijts(k)='Wet deposition of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_wet_depo'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
        SELECT CASE (trname(n))
        CASE ('Clay')

C???? can this be replaced with calls to set_diag_rad?
          IF (diag_rad /= 1) THEN
c dust optical thickness of four clay sub size classes
            do kr=1,4
              k = k + 1
              ijts_tausub(1,n,kr) = k
              ia_ijts(k) = ia_rad
              lname_ijts(k) = trim(trname(n))//char(48+kr)/
     *             /' optical thickness'
              sname_ijts(k) = 'tau_'//trim(trname(n))//char(48+kr)
              ijts_power(k) = -2
              units_ijts(k) = unit_string(ijts_power(k),' ')
              scale_ijts(k) = 10.**(-ijts_power(k))
              ijts_HasArea(k) = .false.
            end do
c dust clear sky optical thickness of four clay sub size classes
            do kr=1,4
              k = k + 1
              ijts_tausub(2,n,kr) = k
              ia_ijts(k) = ia_rad
              lname_ijts(k) = trim(trname(n))//char(48+kr)/
     *             /' CS optical thickness'
              sname_ijts(k) = 'tau_CS_'//trim(trname(n))//char(48+kr)
              dname_ijts(k) = 'clrsky'
              ijts_power(k) = -2
              units_ijts(k) = unit_string(ijts_power(k),' ')
              scale_ijts(k) = 10.**(-ijts_power(k))
              ijts_HasArea(k) = .false.
            end do
          ELSE
            DO kr=1,6
              DO n1=1,4
c extinction optical thickness in six solar bands for four clay sub classes
                k=k+1
                ijts_sqexsub(1,kr,n,n1)=k
                ia_ijts(k)=ia_rad
                WRITE(cform,'(A2,I1,A11)') '(A',LEN_TRIM(trname(n)),
     &               ',I1,A26,I1)'
                WRITE(lname_ijts(k),cform) TRIM(trname(n)),n1,
     &               ' SW total extinction band ',kr
                WRITE(cform,'(A11,I1,A4)') '(A8,I1,A1,A',
     &               LEN_TRIM(trname(n)),',I1)'
                WRITE(sname_ijts(k),cform) 'ext_band',kr,'_',
     &               TRIM(trname(n)),n1
                ijts_power(k) = -4
                units_ijts(k) = unit_string(ijts_power(k),' ')
                scale_ijts(k) = 10.**(-ijts_power(k))
                ijts_HasArea(k) = .false.
c clear sky extinction optical thickness in six solar bands for four clay
c sub classes
                k=k+1
                ijts_sqexsub(2,kr,n,n1)=k
                ia_ijts(k)=ia_rad
                WRITE(cform,'(A2,I1,A11)') '(A',LEN_TRIM(trname(n)),
     &               ',I1,A29,I1)'
                WRITE(lname_ijts(k),cform) TRIM(trname(n)),n1,
     &               ' CS SW total extinction band ',kr
                WRITE(cform,'(A12,I1,A4)') '(A11,I1,A1,A',
     &               LEN_TRIM(trname(n)),',I1)'
                WRITE(sname_ijts(k),cform) 'ext_CS_band',kr,'_',
     &               TRIM(trname(n)),n1
                dname_ijts(k) = 'clrsky'
                ijts_power(k) = -4
                units_ijts(k) = unit_string(ijts_power(k),' ')
                scale_ijts(k) = 10.**(-ijts_power(k))
                ijts_HasArea(k) = .false.
c scattering optical thickness in six solar bands for four clay sub classes
                k=k+1
                ijts_sqscsub(1,kr,n,n1)=k
                ia_ijts(k)=ia_rad
                WRITE(cform,'(A2,I1,A11)') '(A',LEN_TRIM(trname(n)),
     &               ',I1,A28,I1)'
                WRITE(lname_ijts(k),cform) TRIM(trname(n)),n1,
     &               ' SW scatter extinction band ',kr
                WRITE(cform,'(A11,I1,A4)') '(A8,I1,A1,A',
     &               LEN_TRIM(trname(n)),',I1)'
                WRITE(sname_ijts(k),cform) 'sct_band',kr,'_',
     &               TRIM(trname(n)),n1
                ijts_power(k) = -4
                units_ijts(k) = unit_string(ijts_power(k),' ')
                scale_ijts(k) = 10.**(-ijts_power(k))
                ijts_HasArea(k) = .false.
c clear sky scattering optical thickness in six solar bands for four clay
c sub classes
                k=k+1
                ijts_sqscsub(2,kr,n,n1)=k
                ia_ijts(k)=ia_rad
                WRITE(cform,'(A2,I1,A11)') '(A',LEN_TRIM(trname(n)),
     &               ',I1,A31,I1)'
                WRITE(lname_ijts(k),cform) TRIM(trname(n)),n1,
     &               ' CS SW scatter extinction band ',kr
                WRITE(cform,'(A12,I1,A4)') '(A11,I1,A1,A',
     &               LEN_TRIM(trname(n)),',I1)'
                WRITE(sname_ijts(k),cform) 'sct_CS_band',kr,'_',
     &               TRIM(trname(n)),n1
                dname_ijts(k) = 'clrsky'
                ijts_power(k) = -4
                units_ijts(k) = unit_string(ijts_power(k),' ')
                scale_ijts(k) = 10.**(-ijts_power(k))
                ijts_HasArea(k) = .false.
c scattering asymmetry factor in six solar bands for four clay sub classes
                k=k+1
                ijts_sqcbsub(1,kr,n,n1)=k
                ia_ijts(k)=ia_rad
                WRITE(cform,'(A2,I1,A11)') '(A',LEN_TRIM(trname(n)),
     &               ',I1,A26,I1)'
                WRITE(lname_ijts(k),cform) TRIM(trname(n)),n1,
     &               ' SW asymmetry factor band ',kr
                WRITE(cform,'(A11,I1,A4)') '(A8,I1,A1,A',
     &               LEN_TRIM(trname(n)),',I1)'
                WRITE(sname_ijts(k),cform) 'asf_band',kr,'_',
     &               TRIM(trname(n)),n1
                ijts_power(k) = -2
                units_ijts(k) = unit_string(ijts_power(k),' ')
                scale_ijts(k) = 10.**(-ijts_power(k))
                ijts_HasArea(k) = .false.
c clear sky scattering asymmetry factor in six solar bands for four clay
c sub classes
                k=k+1
                ijts_sqcbsub(2,kr,n,n1)=k
                ia_ijts(k)=ia_rad
                WRITE(cform,'(A2,I1,A11)') '(A',LEN_TRIM(trname(n)),
     &               ',I1,A29,I1)'
                WRITE(lname_ijts(k),cform) TRIM(trname(n)),n1,
     &               ' CS SW asymmetry factor band ',kr
                WRITE(cform,'(A12,I1,A4)') '(A11,I1,A1,A',
     &               LEN_TRIM(trname(n)),',I1)'
                WRITE(sname_ijts(k),cform) 'asf_CS_band',kr,'_',
     &               TRIM(trname(n)),n1
                dname_ijts(k) = 'clrsky'
                ijts_power(k) = -2
                units_ijts(k) = unit_string(ijts_power(k),' ')
                scale_ijts(k) = 10.**(-ijts_power(k))
                ijts_HasArea(k) = .false.
              END DO
            END DO
          END IF
c dust shortwave radiative forcing of four clay sub size classes
          do kr=1,4
            k = k + 1
            ijts_fcsub(1,n,kr) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//char(48+kr)/
     *           /' SW radiative forcing'
            sname_ijts(k) = 'swf_'//trim(trname(n))//char(48+kr)
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
          end do
c dust longwave radiative forcing of four clay sub size classes
          do kr=1,4
            k = k + 1
            ijts_fcsub(2,n,kr) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//char(48+kr)/
     *           /' LW radiative forcing'
            sname_ijts(k) = 'lwf_'//trim(trname(n))//char(48+kr)
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
          end do
c dust shortwave radiative forcing at surface of four clay sub size classes
          do kr=1,4
            k = k + 1
            ijts_fcsub(3,n,kr) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//char(48+kr)/
     *           /' SW Surf radiative forcing'
            sname_ijts(k) = 'swf_surf_'//trim(trname(n))//char(48+kr)
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
          end do
c dust longwave radiative forcing at surface of four sub size classes
          do kr=1,4
            k = k + 1
            ijts_fcsub(4,n,kr) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//char(48+kr)/
     *           /' LW Surf radiative forcing'
            sname_ijts(k) = 'lwf_surf_'//trim(trname(n))//char(48+kr)
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
          end do
c dust clear sky shortwave radiative forcing of four clay sub size classes
          do kr=1,4
            k = k + 1
            ijts_fcsub(5,n,kr) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//char(48+kr)/
     *           /' clr sky SW radiative forcing'
            sname_ijts(k) = 'swf_CS_'//trim(trname(n))//char(48+kr)
            dname_ijts(k) = 'clrsky'
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
          end do
c dust clear sky longwave radiative forcing of four clay sub size classes
          do kr=1,4
            k = k + 1
            ijts_fcsub(6,n,kr) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//char(48+kr)/
     *           /' clr sky LW radiative forcing'
            sname_ijts(k) = 'lwf_CS_'//trim(trname(n))//char(48+kr)
            dname_ijts(k) = 'clrsky'
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
          end do
c dust clear sky shortwave radiative forcing at surface of four clay sub size classes
          do kr=1,4
            k = k + 1
            ijts_fcsub(7,n,kr) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//char(48+kr)/
     *           /' clr sky SW Surf radiative forcing'
            sname_ijts(k) = 'swf_CS_surf_'//trim(trname(n))//char(48+kr)
            dname_ijts(k) = 'clrsky'
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
          end do
c dust clear sky longwave radiative forcing at surface of four sub size classes
          do kr=1,4
            k = k + 1
            ijts_fcsub(8,n,kr) = k
            ia_ijts(k) = ia_rad_frc
            lname_ijts(k) = trim(trname(n))//char(48+kr)/
     *           /' clr sky LW Surf radiative forcing'
            sname_ijts(k) = 'lwf_CS_surf_'//trim(trname(n))//char(48+kr)
            dname_ijts(k) = 'clrsky'
            ijts_power(k) = -2
            units_ijts(k) = unit_string(ijts_power(k),'W/m2')
            scale_ijts(k) = 10.**(-ijts_power(k))
            ijts_HasArea(k) = .false.
          end do
        CASE('Silt1','Silt2','Silt3','Silt4',
     &     'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &     'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &     'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &     'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &     'Sil1QuHe','Sil2QuHe','Sil3QuHe')

          call set_diag_rad(n,k)

c dust shortwave radiative forcing
          k = k + 1
          ijts_fc(1,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trim(trname(n))//' SW radiative forcing'
          sname_ijts(k) = 'swf_'//trim(trname(n))
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c dust longwave radiative forcing
          k = k + 1
          ijts_fc(2,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trim(trname(n))//' LW radiative forcing'
          sname_ijts(k) = 'lwf_'//trim(trname(n))
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c dust shortwave radiative forcing at surface
          k = k + 1
          ijts_fc(3,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trim(trname(n))//' SW Surf radiative forcing'
          sname_ijts(k) = 'swf_surf_'//trim(trname(n))
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c dust longwave radiative forcing at surface
          k = k + 1
          ijts_fc(4,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trim(trname(n))//' LW Surf radiative forcing'
          sname_ijts(k) = 'lwf_surf_'//trim(trname(n))
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c dust clear sky shortwave radiative forcing
          k = k + 1
          ijts_fc(5,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trim(trname(n)) //
     &         ' clr sky SW radiative forcing'
          sname_ijts(k) = 'swf_CS_'//trim(trname(n))
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c dust clear sky longwave radiative forcing
          k = k + 1
          ijts_fc(6,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trim(trname(n)) //
     &         ' clr sky LW radiative forcing'
          sname_ijts(k) = 'lwf_CS_'//trim(trname(n))
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c dust clear sky shortwave radiative forcing at surface
          k = k + 1
          ijts_fc(7,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trim(trname(n))
     &         //' clr sky SW surf radiative forcing'
          sname_ijts(k) = 'swf_CS_surf_'//trim(trname(n))
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c dust clear sky longwave radiative forcing at surface
          k = k + 1
          ijts_fc(8,n) = k
          ia_ijts(k) = ia_rad_frc
          lname_ijts(k) = trim(trname(n))
     &         //' clr sky LW surf radiative forcing'
          sname_ijts(k) = 'lwf_CS_surf_'//trim(trname(n))
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
        END SELECT
#endif  /* TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM */

      end select

#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
c**** additional wet deposition diagnostics
      IF (diag_wetdep == 1) THEN
        k=k+1
        ijts_trdpmc(1,n)=k
        lname_ijts(k)='MC Condensation of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_cond_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpmc(2,n)=k
        lname_ijts(k)='Evaporated '//TRIM(trname(n))
     &       //' in MC Downdrafts'
        sname_ijts(k)=TRIM(trname(n))//'_downeva_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpmc(3,n)=k
        lname_ijts(k)='Condensed '//TRIM(trname(n))//' in MC CLW'
        sname_ijts(k)=TRIM(trname(n))//'_conclw_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpmc(4,n)=k
        lname_ijts(k)='Precipitated '//TRIM(trname(n))//' by MC'
        sname_ijts(k)=TRIM(trname(n))//'_precip_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpmc(5,n)=k
        lname_ijts(k)='Reevaporated '//TRIM(trname(n))
     &       //' from MC Precip'
        sname_ijts(k)=TRIM(trname(n))//'_reevap_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpmc(6,n)=k
        lname_ijts(k)='MC Washout of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_washout_mc'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(1,n)=k
        lname_ijts(k)='LS Washout of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_washout_ls'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(2,n)=k
        lname_ijts(k)='Precipitated '//TRIM(trname(n))//' by LS'
        sname_ijts(k)=TRIM(trname(n))//'_precip_ls'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(3,n)=k
        lname_ijts(k)='Condensed '//TRIM(trname(n))// ' in LS CLW'
        sname_ijts(k)=TRIM(trname(n))//'_conclw_ls'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(4,n)=k
        lname_ijts(k)='Reevaporated '//TRIM(trname(n))
     &       //' from LS Precip'
        sname_ijts(k)=TRIM(trname(n))//'_reevap_ls'

        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(5,n)=k
        lname_ijts(k)='Evaporated '//TRIM(trname(n))//' from LS CLW'
        sname_ijts(k)=TRIM(trname(n))//'_clwevap_ls'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k=k+1
        ijts_trdpls(6,n)=k
        lname_ijts(k)='LS Condensation of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_cond_ls'
        ia_ijts(k)=ia_src
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      END IF
#endif
      end do

C**** Additional Special IJ diagnostics
C**** (not necessary associated with a particular tracer)
#ifdef BC_ALB
c BC impact on grain size
c         k = k + 1
c         ijts_alb(2,n) = k
c         ia_ijts(k) = ia_rad????
c         lname_ijts(k) = 'BC impact on grain size'
c         sname_ijts(k) = 'grain_BC'
c         ijts_power(k) = -9
c         units_ijts(k) = unit_string(ijts_power(k),' ')
c         scale_ijts(k) = 10.**(-ijts_power(k))
c BC impact on albedo
        k = k + 1
        ijts_alb(1) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'BC impact on albedo (%)'
        sname_ijts(k) = 'alb_BC'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))

c SW forcing from albedo change
        k = k + 1
        ijts_alb(2) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'BCalb SW radiative forcing'
        sname_ijts(k) = 'swf_BCALB'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.

#endif
#ifdef TRACERS_AEROSOLS_Koch
      IF (diag_rad.eq.1) THEN
        k = k + 1
          ijs_ai = k           ! unused ?????
          ia_ijts(k) = ia_rad
          lname_ijts(k) = 'Aerosol Index'
          sname_ijts(k) = 'ain_CSN'
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
      ENDIF
#endif
#ifdef TRACERS_SPECIAL_Shindell
#ifdef BIOGENIC_EMISSIONS
      k = k+1
        ijs_isoprene=k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Interactive isoprene source'
        sname_ijts(k) = 'Int_isop'
        ijts_power(k) = -10
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
      k = k + 1
        ijs_NO2_1030=k
        ia_ijts(k) = ia_src
        write(lname_ijts(k),'(a18)')'NO2 10:30 trop col'
        sname_ijts(k) = 'NO2_1030'
        dname_ijts(k) = 'NO2_1030c'
        ijts_power(k) = 15
        units_ijts(k) = unit_string(ijts_power(k),'molecules/cm2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
      k = k + 1
        ijs_NO2_1030c=k
        ia_ijts(k) = ia_src ! overridden in TRACER_PRT...
        write(lname_ijts(k),'(a24)')'count NO2 10:30 trop col'
        write(sname_ijts(k),'(a9)')'NO2_1030c'
        ijts_power(k) = 0
        units_ijts(k) = unit_string(ijts_power(k),'number of accum')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
      k = k + 1
        ijs_NO2_1330=k
        ia_ijts(k) = ia_src
        write(lname_ijts(k),'(a18)')'NO2 13:30 trop col'
        sname_ijts(k) = 'NO2_1330'
        dname_ijts(k) = 'NO2_1330c'
        ijts_power(k) = 15
        units_ijts(k) = unit_string(ijts_power(k),'molecules/cm2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
      k = k + 1
        ijs_NO2_1330c=k
        ia_ijts(k) = ia_src ! overridden in TRACER_PRT...
        write(lname_ijts(k),'(a24)')'count NO2 13:30 trop col'
        write(sname_ijts(k),'(a9)')'NO2_1330c'
        ijts_power(k) = 0
        units_ijts(k) = unit_string(ijts_power(k),'number of accum')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
#endif  /* TRACERS_SPECIAL_Shindell */
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      k = k + 1
      ijts_spec(nDustEv1ij)=k
      lname_ijts(k)='No. dust events'
      sname_ijts(k)='no_dust_ev1'
      ia_ijts(k)=ia_src
      scale_ijts(k)=Sday/Dtsrc
      units_ijts(k)='1/d'
      ijts_HasArea(k) = .false.
      k = k + 1
      ijts_spec(nDustEv2ij)=k
      lname_ijts(k)='No. dust events above threshold wind'
      sname_ijts(k)='no_dust_ev2'
      ia_ijts(k)=ia_src
      scale_ijts(k)=Sday/Dtsrc
      units_ijts(k)='1/d'
      ijts_HasArea(k) = .false.
      k = k + 1
      ijts_spec(nDustWthij)=k
      lname_ijts(k)='Threshold velocity for dust emission'
      sname_ijts(k)='wtrsh'
      ia_ijts(k)=ia_src
      scale_ijts(k)=1.
      units_ijts(k)='m/s'
      ijts_HasArea(k) = .false.
#endif

#ifdef TRACERS_AMP
      do n=1,ntm
      select case(trname(n))
      CASE('M_AKK_SU','M_ACC_SU',
     &     'M_BC1_BC','M_OCC_OC','M_BOC_BC','M_BOC_OC')
        select case(trname(n))
        case('M_AKK_SU','M_ACC_SU')
        k = k + 1
          ijts_3Dsource(nVolcanic,n)=k
          ia_ijts(k) = ia_src
          lname_ijts(k) = 'Emission volcano '//trim(trname(n))
          sname_ijts(k) = 'Emission_volcano_'//trim(trname(n))
          ijts_power(k) = -15
          units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
          scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c Surface industrial emissions
        do kr=1,ntsurfsrc(n_SO2)
          k = k + 1
            ijts_source(kr,n) = k
            ia_ijts(k) = ia_src
            sname_ijts(k) = trim(trname(n))//'_src_'//
     &                      trim(ssname(n_SO2,kr))
            lname_ijts(k) = trim(trname(n))//' source from '//
     &                      trim(ssname(n_SO2,kr))
            ijts_power(k) = -15
            units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
            scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        enddo
        case('M_BC1_BC','M_OCC_OC')
c Surface industrial emissions
       do kr=1,ntsurfsrc(n)
        k = k + 1
        ijts_source(kr,n) = k  
        ia_ijts(k) = ia_src
            ia_ijts(k) = ia_src
            sname_ijts(k) = trim(trname(n))//'_src_'//
     &                      trim(ssname(n,kr))
            lname_ijts(k) = trim(trname(n))//' source from '//
     &                      trim(ssname(n,kr))
        ijts_power(k) = -15.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        end do
        end select
      k = k + 1
        ijts_3Dsource(nBiomass,n)=k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Emission biomass '//trim(trname(n))
        sname_ijts(k) = 'Emission_biomass_'//trim(trname(n))
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c- interactive sources diagnostic
      CASE('M_DD1_DU','M_SSA_SS','M_SSC_SS','M_DD2_DU','M_SSS_SS')
      k = k + 1
        ijts_isrc(1,n) = k
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Emission_'//trim(trname(n))
        sname_ijts(k) = 'Emission_'//trim(trname(n))
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      CASE('N_AKK_1 ','N_ACC_1 ','N_DD1_1 ','N_DS1_1 ','N_DD2_1 ',
     *     'N_DS2_1 ','N_OCC_1 ','N_BC1_1 ',
     *     'N_BC2_1 ','N_BC3_1 ','N_DBC_1 ','N_BOC_1 ','N_BCS_1 ',
     *     'N_MXX_1 ','N_OCS_1 ')
      IF ( AMP_DIAG_FC == 2 ) THEN
cc shortwave radiative forcing
      k = k + 1
        ijts_fc(1,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' SW rad forcing'
        sname_ijts(k) = 'swf_'//TRIM(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c longwave radiative forcing
      k = k + 1
        ijts_fc(2,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' LW rad forcing'
        sname_ijts(k) = 'lwf_'//TRIM(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c shortwave surface radiative forcing
      k = k + 1
        ijts_fc(3,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' SW surf forc'
        sname_ijts(k) = 'swf_surf_'//TRIM(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c longwave surface radiative forcing
      k = k + 1
        ijts_fc(4,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' LW surf forc'
        sname_ijts(k) = 'lwf_surf_'//TRIM(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c clear sky shortwave radiative forcing
      k = k + 1
        ijts_fc(5,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' SW cs forc'
        sname_ijts(k) = 'swf_CS_'//TRIM(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c clear sky longwave radiative forcing
      k = k + 1
        ijts_fc(6,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = TRIM(trname(n))//' LW CS forc'
        sname_ijts(k) = 'lwf_CS_'//TRIM(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
      ENDIF

c Special Radiation Diagnostic
      call set_diag_rad(n,k)

      end select
      end do

c - Tracer independent Diagnostic
      IF ( AMP_DIAG_FC == 1 ) THEN
        n=1    !  really? why use ijts_fc then?
cc shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'AMP SW radiative forcing'
        sname_ijts(k) = 'swf_AMP'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'AMP LW radiative forcing'
        sname_ijts(k) = 'lwf_AMP'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c shortwave surface radiative forcing
        k = k + 1
        ijts_fc(3,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'AMP SW surface rad forcing'
        sname_ijts(k) = 'swf_surf_AMP'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c longwave surface radiative forcing
        k = k + 1
        ijts_fc(4,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'AMP LW surface rad forcing'
        sname_ijts(k) = 'lwf_surf_AMP'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'AMP clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS_AMP'
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ia_ijts(k) = ia_rad_frc
        lname_ijts(k) = 'AMP clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS_AMP'
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
      ENDIF
c end special radiation diagnostic

c - Tracer independent Diagnostic (stays here if 2D, moves to ijlt if 3D)
c      do L=1,1    !LTOP
c      do m=1,NBINS
c        k = k + 1
c         ijts_AMPpdf(l,m)=k
c         ia_ijts(k) = ia_src
c         write(lname_ijts(k),'(a15,i2.2,i2.2)') 'NUMB_PDF BIN L=',L,M
c         write(sname_ijts(k),'(a9,i2.2,i2.2)') 'N_PDF_BIN',L,M
c         ijts_power(k) = -2
c         units_ijts(k) = unit_string(ijts_power(k),'Numb.')
c         scale_ijts(k) = 10.**(-ijts_power(k))
c      end do
c      end do
#endif

c
c Append some denominator fields if necessary
c
      if(any(dname_ijts(1:k).eq.'clrsky')) then
        k = k + 1
        ijts_clrsky = k
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'CLEAR SKY FRACTION'
        sname_ijts(k) = 'clrsky'
        units_ijts(k) = '%'
        scale_ijts(k) = 100.
        ijts_HasArea(k) = .false.
      endif

      if(any(dname_ijts(1:k).eq.'ocnfr')) then
        k = k + 1
        ijts_pocean = k
        lname_ijts(k) = 'OCEAN FRACTION'
        units_ijts(k) = '%'
        sname_ijts(k) = 'ocnfr'
        ia_ijts(k) = ia_src     ! ia_ij(ij_pocean) is not initialized yet :(
        scale_ijts(k) = 100.
        ijts_HasArea(k) = .false.
      endif

      ktaijs_out = k
      if (ktaijs_out .gt. ktaijs) then
        if (AM_I_ROOT())
     *       write (6,*)'ijt_defs: Increase ktaijs=',ktaijs
     *       ,' to at least ',ktaijs_out
        call stop_model('ktaijs too small',255)
      end if

c find indices of denominators
      call FindStrings(dname_ijts,sname_ijts,denom_ijts,ktaijs_out)
c      do k=1,ktaijs_out
c        if(len_trim(dname_ijts(k)).gt.0) then
c          do kk=ktaijs_out,1,-1
c            if(trim(sname_ijts(kk)).eq.trim(dname_ijts(k))) then
c              denom_ijts(k) = kk
c              exit
c            endif
c          enddo
c          if(denom_ijts(k).eq.0) then
c            if(am_i_root()) then
c              write(6,*)
c     &             'init_ijts_diag: denominator '//trim(dname_ijts(k))
c              write(6,*) 'not found for field '//trim(sname_ijts(k))
c            endif
c            call stop_model('init_ijts_diag: denominator not found',255)
c          endif
c        endif
c      enddo

#endif /* TRACERS_ON */

      return
      end subroutine init_ijts_diag

#ifdef TRACERS_ON
      subroutine FindStrings(StringsToFind,ListOfStrings,Indices,n)
!@sum FindStrings finds the positions of a list of strings in a 2nd list.
!     Needs optimization.
      use diag_com, only : sname_strlen
      implicit none
      integer :: n
      character(len=sname_strlen), dimension(n) ::
     &     StringsToFind,ListOfStrings
      integer, dimension(n) :: Indices
      integer :: k,kk
      logical :: found
      do k=1,n
        if(len_trim(StringsToFind(k)).gt.0) then
          found = .false.
          do kk=1,n
            if(trim(ListOfStrings(kk)).eq.trim(StringsToFind(k))) then
              Indices(k) = kk
              found = .true.
              exit
            endif
          enddo
          if(.not.found) then
            write(6,*) 'FindStrings: string '//
     &           trim(StringsToFind(k))//' not found'
            call stop_model('FindStrings: string not found',255)
          endif
        endif
      enddo
      end subroutine FindStrings
#endif

      subroutine set_diag_rad(n,k)
!@sum set_diag_rad sets special rad diags for aerosols
!@auth Dorothy Koch
      USE TRACER_COM
#ifdef TRACERS_ON
      USE TRDIAG_COM
#endif /* TRACERS_ON */
      USE DIAG_COM
      implicit none
      integer, intent(inout) :: k
      integer, intent(in) :: n
      integer kr
      character*50 :: unit_string
      character*17 :: cform

#ifdef TRACERS_ON
      IF (diag_rad /= 1) THEN
c optical thickness
        k = k + 1
        ijts_tau(1,n) = k
        ia_ijts(k) = ia_rad
        lname_ijts(k) = trim(trname(n))//' optical thickness'
        sname_ijts(k) = 'tau_'//trim(trname(n))
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
c clear sky optical thickness
        k = k + 1
        ijts_tau(2,n) = k
        ia_ijts(k) = ia_rad
        lname_ijts(k) = trim(trname(n))//' clr sky optical thickness'
        sname_ijts(k) = 'tau_CS_'//trim(trname(n))
        dname_ijts(k) = 'clrsky'
        ijts_power(k) = -2
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
        ijts_HasArea(k) = .false.
      ELSE
        DO kr=1,6
c extinction optical thickness in six solar bands
          k=k+1
          ijts_sqex(1,kr,n)=k
          ia_ijts(k)=ia_rad
          WRITE(cform,'(A2,I1,A8)') '(A',LEN_TRIM(trname(n)),
     &         ',A26,I1)'
          WRITE(lname_ijts(k),cform) TRIM(trname(n)),
     &         ' SW total extinction band ',kr
          WRITE(cform,'(A11,I1,A1)') '(A8,I1,A1,A',
     &         LEN_TRIM(trname(n)),')'
          WRITE(sname_ijts(k),cform) 'ext_band',kr,'_',TRIM(trname(n))
          ijts_power(k) = -4
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c clear sky extinction optical thickness in six solar bands
          k=k+1
          ijts_sqex(2,kr,n)=k
          ia_ijts(k)=ia_rad
          WRITE(cform,'(A2,I1,A8)') '(A',LEN_TRIM(trname(n)),
     &         ',A29,I1)'
          WRITE(lname_ijts(k),cform) TRIM(trname(n)),
     &         ' CS SW total extinction band ',kr
          WRITE(cform,'(A12,I1,A1)') '(A11,I1,A1,A',
     &         LEN_TRIM(trname(n)),')'
          WRITE(sname_ijts(k),cform) 'ext_CS_band',kr,'_',
     &         TRIM(trname(n))
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -4
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c scattering optical thickness in six solar bands
          k=k+1
          ijts_sqsc(1,kr,n)=k
          ia_ijts(k)=ia_rad
          WRITE(cform,'(A2,I1,A8)') '(A',LEN_TRIM(trname(n)),
     &         ',A28,I1)'
          WRITE(lname_ijts(k),cform) TRIM(trname(n)),
     &         ' SW scatter extinction band ',kr
          WRITE(cform,'(A11,I1,A1)') '(A8,I1,A1,A',
     &         LEN_TRIM(trname(n)),')'
          WRITE(sname_ijts(k),cform) 'sct_band',kr,'_',TRIM(trname(n))
          ijts_power(k) = -4
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c clear sky scattering optical thickness in six solar bands
          k=k+1
          ijts_sqsc(2,kr,n)=k
          ia_ijts(k)=ia_rad
          WRITE(cform,'(A2,I1,A8)') '(A',LEN_TRIM(trname(n)),
     &         ',A31,I1)'
          WRITE(lname_ijts(k),cform) TRIM(trname(n)),
     &         ' CS SW scatter extinction band ',kr
          WRITE(cform,'(A12,I1,A1)') '(A11,I1,A1,A',
     &         LEN_TRIM(trname(n)),')'
          WRITE(sname_ijts(k),cform) 'sct_CS_band',kr,'_',
     &         TRIM(trname(n))
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -4
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c scattering asymmetry factor in six solar bands
          k=k+1
          ijts_sqcb(1,kr,n)=k

          ia_ijts(k)=ia_rad
          WRITE(cform,'(A2,I1,A8)') '(A',LEN_TRIM(trname(n)),
     &         ',A26,I1)'
          WRITE(lname_ijts(k),cform) TRIM(trname(n)),
     &         ' SW asymmetry factor band ',kr
          WRITE(cform,'(A11,I1,A1)') '(A8,I1,A1,A',
     &         LEN_TRIM(trname(n)),')'
          WRITE(sname_ijts(k),cform) 'asf_band',kr,'_',TRIM(trname(n))
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
c clear sky scattering asymmetry factor in six solar bands
          k=k+1
          ijts_sqcb(2,kr,n)=k
          ia_ijts(k)=ia_rad
          WRITE(cform,'(A2,I1,A8)') '(A',LEN_TRIM(trname(n)),
     &         ',A29,I1)'
          WRITE(lname_ijts(k),cform) TRIM(trname(n)),
     &         ' CS SW asymmetry factor band ',kr
          WRITE(cform,'(A12,I1,A1)') '(A11,I1,A1,A',
     &         LEN_TRIM(trname(n)),')'
          WRITE(sname_ijts(k),cform) 'asf_CS_band',kr,'_',
     &         TRIM(trname(n))
          dname_ijts(k) = 'clrsky'
          ijts_power(k) = -2
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          ijts_HasArea(k) = .false.
        END DO
      END IF
#endif

      return
      end subroutine set_diag_rad

      subroutine init_ijlts_diag
!@sum init_ijlts_diag Initialise lat/lon/height tracer diags
!@auth Gavin Schmidt
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      USE TRACER_COM
#ifdef TRACERS_ON
      USE TRDIAG_COM
#endif /* TRACERS_ON */
#ifdef HTAP_LIKE_DIAGS
      USE MODEL_COM, only: dtsrc
#endif
#ifdef FAA_emission
      USE MODEL_COM, only: dtsrc
#endif
      USE DIAG_COM
#ifdef SOA_DIAGS
      use tracers_soa, only: issoa
#endif  /* SOA_DIAGS */
      implicit none
      integer k,n,i
      character*50 :: unit_string

#ifdef TRACERS_ON
      ir_ijlt = ir_log2  ! default
      ia_ijlt = ia_src   ! default
#ifdef TRACERS_AMP
      ijlt_AMPext(:)=0
      ijlt_AMPm(:,:)=0
#endif

      k=0
C**** use this routine to set 3D tracer-related diagnostics.

C**** some tracer specific 3D arrays
      do n=1,ntm
        select case(trname(n))

#ifdef TRACERS_DUST
      CASE('Clay','Silt1','Silt2','Silt3','Silt4')
        k = k + 1
         ijlt_3Dtau(n)=k
         ia_ijlt(k) = ia_rad
         lname_ijlt(k) = trim(trname(n))//' tau'
         sname_ijlt(k) = 'tau_3D_'//trname(n)
         ijlt_power(k) = -2
         units_ijlt(k) = unit_string(ijlt_power(k),' ')
         scale_ijlt(k) = 10.**(-ijlt_power(k))
#endif

#ifdef TRACERS_AMP
c- 3D diagnostic per mode
      CASE('N_AKK_1 ','N_ACC_1 ','N_DD1_1 ','N_DS1_1 ','N_DD2_1 ',
     *     'N_DS2_1 ','N_OCC_1 ','N_BC1_1 ',
     *     'N_BC2_1 ','N_BC3_1 ','N_DBC_1 ','N_BOC_1 ','N_BCS_1 ',
     *     'N_MXX_1 ','N_OCS_1 ')
        k = k + 1
         ijlt_AMPm(1,n)=k
         lname_ijlt(k) = TRIM(trname(n))//' DIAM'
         sname_ijlt(k) = 'DIAM_'//TRIM(trname(n))
         ijlt_power(k) = -2.
         units_ijlt(k) = unit_string(ijlt_power(k),'m')
         scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
         ijlt_AMPm(2,n)=k
         lname_ijlt(k) = TRIM(trname(n))//' ACTI'
         sname_ijlt(k) = 'ACTI3D_'//TRIM(trname(n))
         ijlt_power(k) = -2.
         units_ijlt(k) = unit_string(ijlt_power(k),'Numb.')
         scale_ijlt(k) = 10.**(-ijlt_power(k))
#endif
      end select
      end do

C**** 3D tracer-related arrays but not attached to any one tracer

#ifdef TRACERS_SPECIAL_Shindell
#ifdef HTAP_LIKE_DIAGS
      k = k + 1
        ijlt_OH=k
        lname_ijlt(k) = 'OH mixing ratio'
        sname_ijlt(k) = 'OH_vmr'
        ijlt_power(k) = -10
        units_ijlt(k) = unit_string(ijlt_power(k),'V/V air')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
#else
      k = k + 1
        ijlt_OH=k
        lname_ijlt(k) = 'OH concentration'
        sname_ijlt(k) = 'OH_con'
        ijlt_power(k) = 5
        units_ijlt(k) = unit_string(ijlt_power(k),'molecules/cm3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
#endif
      k = k + 1
        ijlt_NO3=k
        lname_ijlt(k) = 'NO3 concentration'
        sname_ijlt(k) = 'NO3_con'
        ijlt_power(k) = 5
        units_ijlt(k) = unit_string(ijlt_power(k),'molecules/cm3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_HO2=k
        lname_ijlt(k) = 'HO2 concentration'
        sname_ijlt(k) = 'HO2_con'
        ijlt_power(k) = 7
        units_ijlt(k) = unit_string(ijlt_power(k),'molecules/cm3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_JO1D=k
        lname_ijlt(k) = 'Ox to O1D photolysis rate'
        sname_ijlt(k) = 'JO1D'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_JNO2=k
        lname_ijlt(k) = 'NO2 photolysis rate'
        sname_ijlt(k) = 'JNO2'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_JH2O2=k
        lname_ijlt(k) = 'H2O2 photolysis rate'
        sname_ijlt(k) = 'JH2O2'
        ijlt_power(k) = 2
        units_ijlt(k) = unit_string(ijlt_power(k),'s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
! #ifdef FAA_emission
!       k = k + 1
!         ijlt_FAANOX=k
!         lname_ijlt(k) = 'NOx_FAA_emission'
!         sname_ijlt(k) = 'NOX_FAA'
!        ijlt_power(k) = -12
!        units_ijlt(k) = unit_string(ijlt_power(k),'kg/s*m^2')
!        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
! #endif
#ifdef HTAP_LIKE_DIAGS
      k = k + 1
        ijlt_COp=k
        lname_ijlt(k) = 'CO production rate'
        sname_ijlt(k) = 'COprod'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
      k = k + 1
        ijlt_COd=k
        lname_ijlt(k) = 'CO destruction rate'
        sname_ijlt(k) = 'COdest'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
      k = k + 1
        ijlt_Oxp=k
        lname_ijlt(k) = 'Ox production rate'
        sname_ijlt(k) = 'Oxprod'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
      k = k + 1
        ijlt_Oxd=k
        lname_ijlt(k) = 'Ox destruction rate'
        sname_ijlt(k) = 'Oxdest'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
      k = k + 1
        ijlt_CH4d=k
        lname_ijlt(k) = 'CH4 destruction rate'
        sname_ijlt(k) = 'CH4dest'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
#endif /* HTAP_LIKE_DIAGS */
#ifdef ACCMIP_LIKE_DIAGS
      k = k + 1
        ijlt_OxpHO2=k
        lname_ijlt(k) = 'Ox prod rate via HO2+NO'
        sname_ijlt(k) = 'OxpHO2'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_OxpCH3O2=k
        lname_ijlt(k) = 'Ox prod rate via CH3O2+NO'
        sname_ijlt(k) = 'OxpCH3O2'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_OxpRO2=k
        lname_ijlt(k) = 'Ox prod rate via RO2+NO'
        sname_ijlt(k) = 'OxpRO2'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_OxlOH=k
        lname_ijlt(k) = 'Ox loss rate via OH'
        sname_ijlt(k) = 'OxlOH'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_OxlHO2=k
        lname_ijlt(k) = 'Ox loss rate via HO2'
        sname_ijlt(k) = 'OxlHO2'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_OxlALK=k
        lname_ijlt(k) = 'Ox loss rate via Alkenes'
        sname_ijlt(k) = 'OxlALK'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_pO1D=k
        lname_ijlt(k) = 'O1D production from ozone'
        sname_ijlt(k) = 'pO1d'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_pOH=k
        lname_ijlt(k) = 'OH production from O1D+H2O'
        sname_ijlt(k) = 'pOH'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'mole m-3 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_NOxLgt=k
        lname_ijlt(k) = 'NOx production from Lightning'
        sname_ijlt(k) = 'NOx_Lightn'
        ijlt_power(k) = -15
        units_ijlt(k) = unit_string(ijlt_power(k),'kg(N) m-2 s-1')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_NOvmr=k
        lname_ijlt(k) = 'NO mixing ratio'
        sname_ijlt(k) = 'NO_vmr'
        ijlt_power(k) = -10 ! to match NOx
        units_ijlt(k) = unit_string(ijlt_power(k),'V/V air')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_NO2vmr=k
        lname_ijlt(k) = 'NO2 mixing ratio'
        sname_ijlt(k) = 'NO2_vmr'
        ijlt_power(k) = -10 ! to match NOx
        units_ijlt(k) = unit_string(ijlt_power(k),'V/V air')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
#ifdef TRACERS_AEROSOLS_Koch
      k = k + 1
        ijlt_prodSO4aq=k
        lname_ijlt(k) = 'SO4 aqueous chem source 3D'
        sname_ijlt(k) = 'SO4aqSrc3D'
        ijlt_power(k) = -15 ! to match ijts 2D
        units_ijlt(k) = unit_string(ijlt_power(k),'kg/s*m^2')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
      k = k + 1
        ijlt_prodSO4gs=k
        lname_ijlt(k) = 'SO4 gas phase source 3D'
        sname_ijlt(k) = 'SO4gasSrc3D'
        ijlt_power(k) = -15 ! to match ijts 2D
        units_ijlt(k) = unit_string(ijlt_power(k),'kg/s*m^2')
        scale_ijlt(k) = 10.**(-ijlt_power(k))/DTsrc
#endif /* TRACERS_AEROSOLS_Koch */
#endif /* ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_SPECIAL_Shindell */

#ifdef SOA_DIAGS
      k = k + 1
        ijlt_soa_changeL_isoprene=k
        lname_ijlt(k) = 'changeL of isoprene'
        sname_ijlt(k) = 'SOA_changeL_isoprene'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_changeL_terpenes=k
        lname_ijlt(k) = 'changeL of terpenes'
        sname_ijlt(k) = 'SOA_changeL_terpenes'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_voc2nox=k
        lname_ijlt(k) = 'VOC/NOx ratio'
        sname_ijlt(k) = 'SOA_voc2nox'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ppbC/ppb')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_pcp=k
        lname_ijlt(k) = 'Total non-volatile SOA-absorbing mass'
        sname_ijlt(k) = 'SOA_pcp'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_aerotot=k
        lname_ijlt(k) = 'PCP plus SOA'
        sname_ijlt(k) = 'SOA_aerotot'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3 per MW')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_aerotot_gas=k
        lname_ijlt(k) = 'Gas-phase semivolatile potential SOA'
        sname_ijlt(k) = 'SOA_aerotot_gas'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3 per MW')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_xmf_isop=k
        lname_ijlt(k) = 'Molar fraction of isoprene SOA'
        sname_ijlt(k) = 'SOA_xmf_isop'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'fraction')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_xmf_apin=k
        lname_ijlt(k) = 'Molar fraction of a-pinene SOA'
        sname_ijlt(k) = 'SOA_xmf_apin'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'fraction')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_zcoef_isop=k
        lname_ijlt(k) = 'Activity coefficient for isoprene SOA'
        sname_ijlt(k) = 'SOA_zcoef_isop'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'dimensionless')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_zcoef_apin=k
        lname_ijlt(k) = 'Activity coefficient for a-pinene SOA'
        sname_ijlt(k) = 'SOA_zcoef_apin'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'dimensionless')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_meanmw=k
        lname_ijlt(k) = 'Mean organic aerosol molecular weight'
        sname_ijlt(k) = 'SOA_meanmw'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'g/mol')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_iternum=k
        lname_ijlt(k) = 'Total iterations for SOA calculations'
        sname_ijlt(k) = 'SOA_iternum'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'count')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_soa_m0=k
        lname_ijlt(k) = 'Final M0 value'
        sname_ijlt(k) = 'SOA_M0'
        ijlt_power(k) = 0
        units_ijlt(k) = unit_string(ijlt_power(k),'ug m-3')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      do i=1,nsoa
        k = k + 1
          ijlt_soa_y0_ug_g(i)=k
          lname_ijlt(k) = 'y0_ug of '//trim(trname(issoa(i)-1))
          sname_ijlt(k) = 'SOA_y0_ug_'//trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_y0_ug_a(i)=k
          lname_ijlt(k) = 'y0_ug of '//trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_y0_ug_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_y_ug_g(i)=k
          lname_ijlt(k) = 'y_ug of '//trim(trname(issoa(i)-1))
          sname_ijlt(k) = 'SOA_y_ug_'//trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_y_ug_a(i)=k
          lname_ijlt(k) = 'y_ug of '//trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_y_ug_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_changeL_g_before(i)=k
          lname_ijlt(k) = 'changeL of '//trim(trname(issoa(i)-1))//
     &                    ' before SOA'
          sname_ijlt(k) = 'SOA_changeL_before_'//
     &                    trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_changeL_a_before(i)=k
          lname_ijlt(k) = 'changeL of '//trim(trname(issoa(i)))//
     &                    ' before SOA'
          sname_ijlt(k) = 'SOA_changeL_before_'//
     &                    trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_changeL_g_after(i)=k
          lname_ijlt(k) = 'changeL of '//trim(trname(issoa(i)-1))//
     &                    ' after SOA'
          sname_ijlt(k) = 'SOA_changeL_after_'//
     &                    trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_changeL_a_after(i)=k
          lname_ijlt(k) = 'changeL of '//trim(trname(issoa(i)))//
     &                    ' after SOA'
          sname_ijlt(k) = 'SOA_changeL_after_'//
     &                    trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_apartmass(i)=k
          lname_ijlt(k) = 'Effective apartmass of '//
     &                     trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_apartmass_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'dimensionless')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_kpart(i)=k
          lname_ijlt(k) = 'Partitioning coefficient of '//
     &                    trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_kpart_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'m3/ug')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_kp(i)=k
          lname_ijlt(k) = 'Final partitioning coefficient of '//
     &                    trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_kp_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'m3/ug')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_soamass(i)=k
          lname_ijlt(k) = 'Potential SOA mass of '//
     &                    trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_soamass_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_partfact(i)=k
          lname_ijlt(k) = 'Final partfact value of '//
     &                    trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_partfact_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'dimensionless')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_evap(i)=k
          lname_ijlt(k) = 'Evaporation of '//
     &                    trim(trname(issoa(i)))
          sname_ijlt(k) = 'SOA_evap_'//trim(trname(issoa(i)))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_cond(i)=k
          lname_ijlt(k) = 'Condensation of pre-existing '//
     &                    trim(trname(issoa(i)-1))
          sname_ijlt(k) = 'SOA_cond_'//trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
        k = k + 1
          ijlt_soa_chem(i)=k
          lname_ijlt(k) = 'Condensation of same-step produced '//
     &                    trim(trname(issoa(i)-1))
          sname_ijlt(k) = 'SOA_chem_'//trim(trname(issoa(i)-1))
          ijlt_power(k) = 0
          units_ijlt(k) = unit_string(ijlt_power(k),'ug/m3')
          scale_ijlt(k) = 10.**(-ijlt_power(k))
      enddo
#endif  /* SOA_DIAGS */

#ifdef TRACERS_AMP
      k = k + 1
        ijlt_AMPext(1)=k
        lname_ijlt(k) = 'N_SSA ACTI'
        sname_ijlt(k) = 'ACTI3D_N_SSA_1'
        ijlt_power(k) = -2
        units_ijlt(k) = unit_string(ijlt_power(k),'Numb.')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_AMPext(2)=k
        lname_ijlt(k) = 'N_SSC ACTI'
        sname_ijlt(k) = 'ACTI3D_N_SSC_1'
        ijlt_power(k) = -2
        units_ijlt(k) = unit_string(ijlt_power(k),'Numb.')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_AMPext(3)=k
        lname_ijlt(k) = 'N_SSA DIAM'
        sname_ijlt(k) = 'DIAM_N_SSA_1'
        ijlt_power(k) = -2
        units_ijlt(k) = unit_string(ijlt_power(k),'m')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_AMPext(4)=k
        lname_ijlt(k) = 'N_SSC DIAM'
        sname_ijlt(k) = 'DIAM_N_SSC_1'
        ijlt_power(k) = -2
        units_ijlt(k) = unit_string(ijlt_power(k),'m')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_AMPext(5)=k
        lname_ijlt(k) = 'N_SSA_1'
        sname_ijlt(k) = 'N_SSA_1'
        ijlt_power(k) = -10
        units_ijlt(k) = unit_string(ijlt_power(k),'Numb.')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
      k = k + 1
        ijlt_AMPext(6)=k
        lname_ijlt(k)= 'N_SSC_1'
        sname_ijlt(k)= 'N_SSC_1'
        ijlt_power(k) = -10
        units_ijlt(k) = unit_string(ijlt_power(k),'Numb.')
        scale_ijlt(k) = 10.**(-ijlt_power(k))
#endif

      if (k .gt. ktaijl) then
       if (AM_I_ROOT())
     *       write (6,*)'ijlt_defs: Increase ktaijl=',ktaijl
     *       ,' to at least ',k
        call stop_model('ktaijl too small',255)
      end if
#endif /* TRACERS_ON */

      return
      end subroutine init_ijlts_diag

      SUBROUTINE tracer_IC
!@sum tracer_IC initializes tracers when they are first switched on
!@auth Jean Lerner
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT,readt_parallel,
     &     readt8_column, skip_parallel
      USE PARAM, only : get_param, is_set_param
#ifdef TRACERS_ON
      USE CONSTANT, only: mair,rhow,sday,grav,tf,avog,rgas
      USE resolution,ONLY : Im,Jm,Lm,Ls1
      USE MODEL_COM, only: itime,jday,dtsrc,ptop,q,wm,flice,jyear,
     & PMIDL00, itimeI
#ifdef TRACERS_WATER
     &     ,focean
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID,GET,write_parallel
      USE SOMTQ_COM, only : qmom,mz,mzz
      USE TRACER_COM, only: ntm,trm,trmom,itime_tr0,trname,needtrs,
     *     tr_mm,rnsrc,vol2mass
#if (defined TRACERS_AEROSOLS_Koch)||(defined TRACERS_OM_SP)||\
    (defined TRACERS_AMP)
     *     ,imAER,n_SO2,imPI,aer_int_yr,OFFLINE_DMS_SS,OFFLINE_SS
#endif
#ifdef TRACERS_WATER
     *     ,trwm,trw0,tr_wd_TYPE,nWATER,n_HDO,n_H2O18
      USE LANDICE, only : ace1li,ace2li
      USE LANDICE_COM, only : trli0,trsnowli,trlndi,snowli
      USE SEAICE, only : xsi,ace1i
      USE SEAICE_COM, only : rsi,msi,snowi,trsi,trsi0,ssi
      USE LAKES_COM, only : trlake,mwl,mldlk,flake
      USE GHY_COM, only : tr_w_ij,tr_wsn_ij,w_ij
     &     ,wsn_ij,nsn_ij,fr_snow_ij,fearth
      USE FLUXES, only : gtracer
#endif
      USE GEOM, only: axyp,byaxyp,lat2d_dg,lonlat_to_ij
      USE DYNAMICS, only: am,byam  ! Air mass of each box (kg/m^2)
      USE PBLCOM, only: npbl,trabl,qabl,tsavg
#ifdef TRACERS_SPECIAL_Lerner
      USE LINOZ_CHEM_COM, only: tlt0m,tltzm, tltzzm
      USE PRATHER_CHEM_COM, only: nstrtc
#endif
      USE FILEMANAGER, only: openunit,closeunit,nameunit
#ifdef TRACERS_SPECIAL_Shindell
      USE RAD_COM, only : chem_tracer_save,rad_to_file,ghg_yr
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     &  ,stratO3_tracer_save
#endif
      USE TRCHEM_Shindell_COM,only:O3MULT,MDOFM,ch4icx,
     &  OxIC,COIC,byO3MULT,PI_run,fix_CH4_chemistry,
     &  PIratio_N,PIratio_CO_T,PIratio_CO_S,PIratio_other
     &  ,use_rad_n2o,use_rad_cfc,use_rad_ch4
#ifdef SHINDELL_STRAT_CHEM
     &  ,ClOxalt,BrOxalt,ClONO2alt,HClalt,N2OICX,CFCIC
     &  ,PIratio_N2O,PIratio_CFC
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      USE TRACER_SOURCES, only:first_mod,first_ncep,avg_model,avg_ncep,
     & PRS_ch4,sum_ncep
#endif
#ifdef SHINDELL_STRAT_EXTRA
      USE TRACER_SOURCES, only:GLTic
#endif
#endif /* TRACERS_SPECIAL_Shindell */
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP)
      USE AEROSOL_SOURCES, only: DMSinput,BCI_src,OCI_src,
     * BCB_src,OCB_src,om2oc,
#ifndef TRACERS_AEROSOLS_SOA
     * OCT_src,n_OCII,
#endif  /* TRACERS_AEROSOLS_SOA */
     * DMS_AER,SS1_AER,SS2_AER,
     * SO2_src_3D,SO2_biosrc_3D,SO2_src,bci_src_3D,
     * lmAER,craft,NH3_src_con,NH3_src_cyc
#endif
#ifdef TRACERS_RADON
       USE AEROSOL_SOURCES, only: rn_src
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      USE tracers_dust,ONLY : hbaij,ricntd
      use trdust_drv, only: tracer_ic_soildust
#endif
#ifdef TRACERS_AMP
      USE AMP_AEROSOL
#endif
#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_GASEXCH_ocean_CO2)
      USE FLUXES, only : gtracer
#ifdef OBIO_ON_GARYocean
      Use AFLUXES, Only: aTRAC
#else
      USE TRACER_GASEXCH_COM, only : atrac=>atrac_loc
#endif
      USE MODEL_COM, only : nstep=>itime
#ifdef constCO2
      USE obio_forc, only : atmCO2
#else
      USE RADPAR, only : xnow
#endif
#endif
#if (!defined(TRACERS_GASEXCH_ocean_CO2)) && defined(TRACERS_GASEXCH_land_CO2)
      USE FLUXES, only : gtracer
      USE RADPAR, only : xnow
#endif
      IMPLICIT NONE
      real*8,parameter :: d18oT_slope=0.45,tracerT0=25
      INTEGER i,n,l,j,iu_data,ipbl,it,lr,m,ls,lt
      CHARACTER*80 title
      CHARACTER*300 out_line
      REAL*8 CFC11ic,conv
      REAL*8 :: trinit =1., tmominit=0.
      real*8 tracerTs

      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) ::
     *                                CO2ic,ic14CO2
      REAL*4, DIMENSION(jm,lm)    ::  N2Oic   !each proc. reads global array
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) ::
     *                                                      CH4ic
#ifdef TRACERS_SPECIAL_Lerner
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: icCFC
      REAL*8 stratm,xlat,pdn,pup
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      REAL*4, DIMENSION(GRID%I_STRT:GRID%I_STOP,
     &                  GRID%J_STRT:GRID%J_STOP,366) ::
     &     DMS_AER_nohalo, SS1_AER_nohalo, SS2_AER_nohalo
#endif
#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_GASEXCH_ocean_CFC)
      REAL*8 :: dummy
#endif

!@param bymair 1/molecular wt. of air = 1/mair
!@param byjm 1./JM
      REAL*8, PARAMETER :: bymair = 1.d0/mair, byjm =1.d0/JM
#ifdef TRACERS_SPECIAL_Shindell
      character*4 ghg_name
      character*80 ghg_file
      real*8, dimension(LM,GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                     GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: ghg_in
!@var imonth dummy index for choosing the right month
!@var ICfactor varying factor for altering initial conditions
      INTEGER imonth, J2
      REAL*8 ICfactor
!@var PRES local nominal pressure for vertical interpolations
      REAL*8, DIMENSION(LM) :: PRES
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      include 'netcdf.inc'
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP)
      integer start(3),count(3),status,ncidu,id1
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     * bci_src1,bci_src2,bci_src3,bci_src4
     * ,oci_src1,oci_src2,oci_src3,oci_src4
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,8) ::
     *  OCI_src5
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) ::
     *  volc_exp
      INTEGER mon_unit, mont,ii,jj,ir,mm,iuc,mmm,ll,iudms
      INTEGER iuc2,lmax
      real*8 carbstuff,ccnv,carb(8)
#endif
#if defined (TRACERS_AEROSOLS_Koch) || defined (TRACERS_AMP)
#ifdef CUBED_SPHERE /* 1-deg volc. emiss */
      real*8 :: volc_lons(360),volc_lats(180),
     &     volc_pup(360,180),volc_emiss(360,180)
#else /* volc. emiss on model grid, 1 extra lat at SP */
      real*8 :: volc_lons(Im),volc_lats(Jm+1),
     &     volc_pup(Im,Jm+1),volc_emiss(Im,Jm+1)
#endif
      real*8 :: x1d(lm),amref(lm),pednref(lm+1),amsum
      real*8, allocatable, dimension(:,:) :: psref
      integer :: iu_ps,file_id,vid,ilon,jlat,volc_ij(2)
#endif

      INTEGER J_0, J_1, I_0, I_1
      INTEGER J_0H, J_1H
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      integer :: initial_GHG_setup

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1,
     *               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     *               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

#ifdef SHINDELL_STRAT_CHEM
      PRES(1:LM)=PMIDL00(1:LM)
#endif
      do n=1,ntm
      if (itime.eq.itime_tr0(n)) then

C**** set some defaults for air mass tracers
      trm(:,J_0:J_1,:,n) = 0.
      trmom(:,:,J_0:J_1,:,n) = 0.

#ifdef TRACERS_WATER
C**** set some defaults for water tracers
      trwm(:,J_0:J_1,:,n)=0. ! cloud liquid water
      trlake(n,:,:,J_0:J_1)=0.
      trsi(n,:,:,J_0:J_1)=0.
      trlndi(n,:,J_0:J_1)=0.
      trsnowli(n,:,J_0:J_1)=0.
      tr_w_ij(n,:,:,:,J_0:J_1)=0.
      tr_wsn_ij(n,:,:,:,J_0:J_1)=0.
#endif
      select case (trname(n))

        case default
          write(6,*) 'In TRACER_IC:',trname(n),' does not exist '
          call stop_model("TRACER_IC",255)

        case ('Air')
          do l=1,lm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*axyp(:,j)
          end do; enddo

        case ('SF6','SF6_c')
          ! defaults ok

        case ('Be7', 'Be10', 'Pb210', 'Rn222')
          ! defaults ok

        case ('CO2')
          call openunit('CO2_IC',iu_data,.true.,.true.)
          CALL READT_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),CO2IC,0)
          call closeunit(iu_data)
          do l=1,lm         !ppmv==>ppmm
          do j=J_0,J_1
            trm(:,j,l,n) = co2ic(:,j,l)*am(l,:,j)*axyp(:,j)*1.54d-6
          enddo; enddo

        case ('N2O')
#ifdef TRACERS_SPECIAL_Lerner
          call openunit('N2O_IC',iu_data,.true.,.true.)
C**** ESMF: Each processor reads the global array: N2Oic
          read (iu_data) title,N2Oic     ! unit is PPMM/(M*AXYP)
          call closeunit(iu_data)
          if (AM_I_ROOT()) write(6,*) title,' read from N2O_IC'
          do l=1,lm         !ppmv==>ppmm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*axyp(:,j)*N2Oic(j,l)
          enddo; enddo
#endif
#ifdef SHINDELL_STRAT_CHEM
         if(use_rad_n2o <= 0)then
           select case(PI_run)
           case(1)     ; ICfactor=PIratio_N2O
           case default; ICfactor=1.d0
           end select
           do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
             trm(i,j,l,n) = N2OICX(i,j,l)*ICfactor
           end do   ; end do   ; end do
         else
           if (is_set_param('initial_GHG_setup')) then
             call get_param('initial_GHG_setup', initial_GHG_setup)
             if (initial_GHG_setup == 1 .and. itime == itimeI) then
               select case(PI_run)
               case(1)     ; ICfactor=PIratio_N2O
               case default; ICfactor=1.d0
               end select
               do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
                 trm(i,j,l,n) = N2OICX(i,j,l)*ICfactor
               end do   ; end do   ; end do
             else
               if(ghg_yr/=0)then; write(ghg_name,'(I4)') ghg_yr
               else; write(ghg_name,'(I4)') jyear; endif
               ghg_file='GHG_IC_'//ghg_name
               call openunit(ghg_file,iu_data,.true.,.true.)
               do m=1,3
                 CALL READT8_COLUMN
     &           (grid,iu_data,NAMEUNIT(iu_data),GHG_IN,0)
                 rad_to_file(m,:,I_0:I_1,J_0:J_1)=
     &           ghg_in(:,I_0:I_1,J_0:J_1)
               enddo
               call closeunit(iu_data)
               do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
                 trm(I,J,L,n) = rad_to_file(3,l,i,j)
               end do   ; end do   ; end do
             end if
           endif
         end if
#endif

        case ('CFC11')   !!! should start April 1
          CFC11ic = 268.D-12*136.5/29.029    !268 PPTV
          do l=1,lm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*axyp(:,j)*CFC11ic
          enddo; enddo
#ifdef TRACERS_SPECIAL_Lerner
C****
C**** Read in first layer distribution; This is used up to about 100 mb
C****
      call openunit('CFCic_Lerner',iu_data,.true.,.true.)
      CALL READT_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),icCFC,0)
      call closeunit(iu_data)
C**** Fill in the tracer; above 100 mb interpolate linearly with P to 0 at top
      stratm = 101.9368
      DO J=J_0,J_1
      DO I=I_0,I_1
        trm(i,j,:,n) = 0.
        trm(i,j,:,n) = 0.
        PUP = STRATM*GRAV
        DO LS=LM,1,-1
          PDN = PUP + AM(ls,I,J)*GRAV
          IF(PDN.GT.10000.d0)  GO TO 450
          trm(I,J,LS,N) =
     *      AM(ls,I,J)*AXYP(I,J)*icCFC(i,j)*.5*(PUP+PDN)/10000.d0
          PUP = PDN
        enddo
  450   CONTINUE
        trm(I,J,LS,N) = AM(ls,I,J)*AXYP(I,J)*icCFC(i,j)*
     *    (1.-.5*(10000.-PUP)*(10000.-PUP)/(10000.*(PDN-PUP)))
        DO LT=1,LS-1
          trm(I,J,LT,N) = AM(lt,I,J)*AXYP(I,J)*icCFC(i,j)
        enddo
      enddo; enddo
#endif


        case ('14CO2')   !!! this tracer is supposed to start 10/16
#ifdef TRACERS_SPECIAL_Lerner
          call get_14CO2_IC(ic14CO2)
          do l=1,lm         !ppmv==>ppmm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*axyp(:,j)*ic14CO2(:,j,l)*1.d-18
          enddo; enddo
#endif

        case ('CH4')
#ifdef TRACERS_SPECIAL_Shindell
         if(use_rad_ch4 <= 0)then
          select case (fix_CH4_chemistry)
          case default
            call get_CH4_IC(0) ! defines trm(:,:,:,n_CH4) within
          case(-1) ! ICs from file...
            call get_CH4_IC(0) ! defines trm(:,:,:,n_CH4) within
            do l=ls1,lm; do j=J_0,J_1; do i=I_0,I_1
              trm(I,J,L,n) = CH4ICX(I,J,L)
            end do   ; end do   ; end do
          end select
#ifdef INTERACTIVE_WETLANDS_CH4
          first_mod(:,:,:)=1
          first_ncep(:)=1
          avg_model(:,:,:)=0.d0
          avg_ncep(:,:,:)=0.d0
          PRS_ch4(:,:,:)=0.d0
          sum_ncep(:,:,:)=0.d0
#endif
         else
           if (is_set_param('initial_GHG_setup')) then
             call get_param('initial_GHG_setup', initial_GHG_setup)
             if (initial_GHG_setup == 1 .and. itime == itimeI) then
               select case (fix_CH4_chemistry)
               case default
                 call get_CH4_IC(0) ! defines trm(:,:,:,n_CH4) within
               case(-1)         ! ICs from file...
                 call get_CH4_IC(0) ! defines trm(:,:,:,n_CH4) within
                 do l=ls1,lm; do j=J_0,J_1; do i=I_0,I_1
                   trm(I,J,L,n) = CH4ICX(I,J,L)
                 end do   ; end do   ; end do
               end select
             else
               if(ghg_yr/=0)then; write(ghg_name,'(I4)') ghg_yr
               else; write(ghg_name,'(I4)') jyear; endif
               ghg_file='GHG_IC_'//ghg_name
               call openunit(ghg_file,iu_data,.true.,.true.)
               do m=1,4
                 CALL READT8_COLUMN(grid,iu_data,NAMEUNIT(iu_data),
     &                GHG_IN,0)
                 rad_to_file(m,:,I_0:I_1,J_0:J_1)=
     &                ghg_in(:,I_0:I_1,J_0:J_1)
               enddo
               call closeunit(iu_data)
               do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
                 trm(I,J,L,n) = rad_to_file(4,l,i,j)
               end do   ; end do   ; end do
             end if
           end if
         end if
         do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
           chem_tracer_save(2,L,I,J)=trm(I,J,L,n)
     &          *byaxyp(i,j)*avog/(tr_mm(n)*2.69e20) ! to atm*cm
         end do   ; end do   ; end do
#endif /* TRACERS_SPECIAL_Shindell */
#ifdef TRACERS_SPECIAL_Lerner
          call get_wofsy_gas_IC(trname(n),CH4ic)
          do l=1,lm         !ppbv==>ppbm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*axyp(:,j)*CH4ic(j,l)*0.552d-9
          enddo; enddo
#endif

        case ('O3')
          do l=1,lm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*axyp(:,j)*20.d-9*vol2mass(n)
          enddo; enddo
#ifdef TRACERS_SPECIAL_Lerner
          do l=lm,lm+1-nstrtc,-1
          lr = lm+1-l
            do j=J_0,J_1
            if (tlt0m(j,lr,5) /= 0.) then
            trm(:,j,l,n) =
     *          tlt0m(j,lr,1)*am(l,:,j)*axyp(:,j)*vol2mass(n)
            trmom(mz,:,j,l,n)  =
     *          tltzm(j,lr,1)*am(l,:,j)*axyp(:,j)*vol2mass(n)
            trmom(mzz,:,j,l,n)  =
     *         tltzzm(j,lr,1)*am(l,:,j)*axyp(:,j)*vol2mass(n)
            end if
            end do
          end do
#endif

#ifdef TRACERS_WATER
      case ('Water', 'H2O18', 'HDO', 'HTO', 'H2O17')

C**** initial atmospheric conc. needs to be defined for each tracer
        select case (trname(n))
        case ('Water')
          trinit=1.
C**** for gradients defined on air mass
          tmominit = 1.
C**** for gradients defined on water mass (should be an option?)
c     tmominit = 0.
        case ('H2O18')        ! d18O=-80
          trinit=0.92d0*trw0(n)
          tmominit = trinit
        case ('H2O17')        ! d18O=-43.15  (D17O=0)
          trinit=.95685d0*trw0(n)
          tmominit = trinit
        case ('HDO')   ! dD=-630
          trinit=0.37d0*trw0(n)
          tmominit = trinit
        case ('HTO')
          trinit=0.
          tmominit = trinit
        end select

        do l=1,lm
        do j=J_0,J_1
          do i=I_0,I_1
            trm(i,j,l,n) =  q(i,j,l)*am(l,i,j)*axyp(i,j)*trinit
            trwm(i,j,l,n)= wm(i,j,l)*am(l,i,j)*axyp(i,j)*trinit
            trmom(:,i,j,l,n) = qmom(:,i,j,l)*am(l,i,j)*axyp(i,j)
     *           *tmominit
          end do
        end do
        end do
        if (HAVE_SOUTH_POLE) then
           do i=2,im
              trm(i,1,:,n) =  trm(1,1,:,n) !poles
              trwm(i, 1,:,n)= trwm(1, 1,:,n) !poles
              trmom(:,i, 1,:,n)=0.
           enddo
        endif
        if (HAVE_NORTH_POLE) then
           do i=2,im
              trm(i,jm,:,n) = trm(1,jm,:,n) !poles
              trwm(i,jm,:,n)= trwm(1,jm,:,n) !poles
              trmom(:,i,jm,:,n)=0.
          enddo
        endif
        if (trname(n).eq."HTO") then ! initialise bomb source
          do l=ls1-1,ls1+1      ! strat. source lat 44 N - 56 N
          do j=J_0,J_1
          do i=I_0,I_1
            if(nint(lat2d_dg(i,j)).ge.44.and.nint(lat2d_dg(i,j)).le.56)
     *           trm(i,j,l,n)= q(i,j,l)*am(l,i,j)*axyp(i,j)*1d10*1d-18
          end do
          end do
          end do
        end if

        do j=J_0,J_1
          do i=I_0,I_1
            tracerTs=trw0(n)
#ifdef TRACERS_SPECIAL_O18
c Define a simple d18O based on Tsurf for GIC, put dD on meteoric water line
            if(trname(n).eq."H2O18") tracerTs=TRW0(n_H2O18)*(1.+1d-3*
     *           ((tsavg(i,j)-(tf+tracerT0))*d18oT_slope))
            if(trname(n).eq."HDO") tracerTs=TRW0(n_HDO)*(1.+(1d-3*
     *           (((tsavg(i,j)-(tf+tracerT0))*d18oT_slope)*8+1d1)))
#endif
C**** lakes
            if (flake(i,j).gt.0) then
              trlake(n,1,i,j)=tracerTs*mldlk(i,j)*rhow*flake(i,j)
     *             *axyp(i,j)
              if (mwl(i,j)-mldlk(i,j)*rhow*flake(i,j)*axyp(i,j).gt.1d-10
     *             *mwl(i,j)) then
                trlake(n,2,i,j)=tracerTs*mwl(i,j)-trlake(n,1,i,j)
              else
                trlake(n,2,i,j)=0.
              end if
              gtracer(n,1,i,j)=trw0(n)
            else !if (focean(i,j).eq.0) then
              trlake(n,1,i,j)=trw0(n)*mwl(i,j)
              trlake(n,2,i,j)=0.
c            else
c              trlake(n,1:2,i,j)=0.
            end if
c**** ice
            if (msi(i,j).gt.0) then
              trsi(n,1,i,j)=trsi0(n)*
     *             (xsi(1)*(snowi(i,j)+ace1i)-ssi(1,i,j))
              trsi(n,2,i,j)=trsi0(n)*
     *             (xsi(2)*(snowi(i,j)+ace1i)-ssi(2,i,j))
              trsi(n,3,i,j)=trsi0(n)*(xsi(3)*msi(i,j)-ssi(3,i,j))
              trsi(n,4,i,j)=trsi0(n)*(xsi(4)*msi(i,j)-ssi(4,i,j))
              gtracer(n,2,i,j)=trsi0(n)
            end if
c**** landice
            if (flice(i,j).gt.0) then
              trlndi(n,i,j)=trli0(n)*(ace1li+ace2li)
              trsnowli(n,i,j)=trli0(n)*snowli(i,j)
              gtracer(n,3,i,j)=trli0(n)
            else
              trlndi(n,i,j)=0.
              trsnowli(n,i,j)=0.
              gtracer(n,3,i,j)=0.
            end if
c**** earth
            !!!if (fearth(i,j).gt.0) then
            if (focean(i,j) < 1.d0) then
              conv=rhow         ! convert from m to kg/m^2
              tr_w_ij  (n,:,:,i,j)=tracerTs*w_ij (:,:,i,j)*conv
              tr_wsn_ij(n,1:nsn_ij(1,i,j),1,i,j)=
     &             tracerTs*wsn_ij(1:nsn_ij(1,i,j),1,i,j)
     &             *fr_snow_ij(1,i,j)*conv
              tr_wsn_ij(n,1:nsn_ij(2,i,j),2,i,j)=
     &             tracerTs*wsn_ij(1:nsn_ij(2,i,j),2,i,j)
     &             *fr_snow_ij(2,i,j)*conv
              !trsnowbv(n,2,i,j)=trw0(n)*snowbv(2,i,j)*conv
              gtracer (n,4,i,j)=trw0(n)
            else
              tr_w_ij  (n,:,:,i,j)=0.
              tr_wsn_ij(n,:,:,i,j)=0.
              !trsnowbv(n,1,i,j)=0.
              !trsnowbv(n,2,i,j)=0.
              gtracer(n,4,i,j)=0.
            end if
          end do
          end do
#ifdef TRACERS_SPECIAL_O18
          if (AM_I_ROOT()) then
            if(trname(n).eq."H2O18") write(6,'(A52,f6.2,A15,f8.4,A18)')
     *            "Initialized trlake tr_w_ij tr_wsn_ij using Tsurf at"
     *           ,tracerT0,"degC, 0 permil",d18oT_slope
     *           ,"permil d18O/degC"
          endif
#endif

#endif /* TRACERS_WATER */

#ifdef TRACERS_SPECIAL_Shindell
        case ('Ox')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(I,J,L,n) = OxIC(I,J,L)
            chem_tracer_save(1,L,I,J)=OxIC(I,J,L)*byO3MULT*byaxyp(i,j)
          end do   ; end do   ; end do
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        case ('stratOx')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(I,J,L,n) = OxIC(I,J,L)
            stratO3_tracer_save(L,I,J)=OxIC(I,J,L)*byO3MULT*byaxyp(i,j)
          end do   ; end do   ; end do
#endif
#endif /* TRACERS_SPECIAL_Shindell */

        case ('NOx')
#ifdef TRACERS_SPECIAL_Shindell
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*1.d-11*ICfactor
#ifdef SHINDELL_STRAT_CHEM
            if(PRES(L).lt.10.)trm(i,j,l,n)=trm(i,j,l,n)*3.d2
#endif
          end do; end do; end do
#endif /* TRACERS_SPECIAL_Shindell */

#if (defined TRACERS_SPECIAL_Shindell) && (defined SHINDELL_STRAT_CHEM)
        case ('ClOx')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      am(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-11*ClOxalt(l)
          end do; end do; end do

        case ('BrOx')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      am(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-11*BrOxalt(l)
          end do; end do; end do

        case ('HCl')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      am(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-11*HClalt(l)
          end do; end do; end do

        case ('ClONO2')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      am(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-11*ClONO2alt(l)
          end do; end do; end do
#endif
        case ('N2O5')
#ifdef TRACERS_SPECIAL_Shindell
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*1.d-12*ICfactor
          end do; end do; end do
#endif
        case ('HNO3')
#ifdef TRACERS_SPECIAL_Shindell
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*1.d-10*ICfactor
#ifdef SHINDELL_STRAT_CHEM
            if(PRES(L).lt.50.and.PRES(L).gt.10.)
     &      trm(i,j,l,n)=trm(i,j,l,n)*1.d2
#endif
          end do; end do; end do
#endif /* TRACERS_SPECIAL_Shindell */
        case ('H2O2')
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*5.d-10
          end do; end do; end do

#ifdef SHINDELL_STRAT_EXTRA
        case ('GLT')
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = GLTic*vol2mass(n)*am(l,i,j)*axyp(i,j)
          end do; end do; end do
#endif

        case ('CH3OOH', 'HCHO')
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*1.d-11
          end do; end do; end do

#ifdef TRACERS_SPECIAL_Shindell
        case ('HO2NO2')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=i_0,i_1
            trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*1.d-12*ICfactor
          end do; end do; end do

        case ('CO')
          do l=1,lm
            select case(PI_run)
            case(1) ! ise scaling
              if(L.le.LS1-1) then
                ICfactor=PIratio_CO_T ! troposphere
              else
                ICfactor=PIratio_CO_S ! stratosphere
              end if
            case default; ICfactor=1.d0
            end select
            do j=J_0,J_1; do i=I_0,I_1
              trm(I,J,L,n) = COIC(I,J,L)*ICfactor
            end do   ; end do
          end do

        case ('codirect')
          ! supposed to start from zero conc for the HTAP TP anyway...
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) = 0.d0
          end do; end do; end do

        case ('PAN')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      am(l,i,j)*axyp(i,j)*vol2mass(n)*4.d-11*ICfactor
          end do; end do; end do

        case ('Isoprene')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      am(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-11*ICfactor
          end do; end do; end do

        case ('AlkylNit')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      am(l,i,j)*axyp(i,j)*vol2mass(n)*2.d-10*ICfactor
          end do; end do; end do

        case('Alkenes')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      am(l,i,j)*axyp(i,j)*vol2mass(n)*4.d-10*ICfactor
          end do; end do; end do

        case('Paraffin')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      am(l,i,j)*axyp(i,j)*vol2mass(n)*5.d-10*ICfactor
          end do; end do; end do

#ifdef TRACERS_TERP
        case ('Terpenes')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      am(l,i,j)*axyp(i,j)*vol2mass(n)*0.d0*ICfactor
          end do; end do; end do
#endif  /* TRACERS_TERP */
#endif /* TRACERS_SPECIAL_Shindell */

#ifdef TRACERS_AEROSOLS_SOA
        case('isopp1g','isopp1a','isopp2g','isopp2a',
     &       'apinp1g','apinp1a','apinp2g','apinp2a')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) =
     &      am(l,i,j)*axyp(i,j)*vol2mass(n)*0.d0*ICfactor
          end do; end do; end do
#endif  /* TRACERS_AEROSOLS_SOA */

#if defined(TRACERS_GASEXCH_ocean_CO2) || defined(TRACERS_GASEXCH_land_CO2)
        case ('CO2n')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
             !units: [am]=kg_air/m2, [axyp]=m2, [tr_mm]=kg_CO2,
             !       [bymair]=1/kg_air, [atmCO2]=ppmv=10^(-6)kg_CO2/kg_air
             !       [vol2mass]=(gr,CO2/moleCO2)/(gr,air/mole air)
#ifdef constCO2
             trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*vol2mass(n)
     .                    * atmCO2*1.d-6
             !!!!trmom(:,i,j,l,n)=0.d0
             gtracer(n,1,i,j) = vol2mass(n)
     .                    * atmCO2 * 1.d-6      !initialize gtracer
#else
             trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*vol2mass(n)
     .                    * xnow(1) * 1.d-6
             gtracer(n,1,i,j) = vol2mass(n)
     .                    * xnow(1) * 1.d-6      !initialize gtracer
#endif
#if defined(TRACERS_GASEXCH_ocean_CO2)
             atrac(i,j,n) = gtracer(n,1,i,j)
#endif
          end do; end do; end do
#endif

#ifdef TRACERS_GASEXCH_ocean_CFC
        case ('CFCn')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            if(l.ge.LS1) then
              trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*vol2mass(n)*2.d-13
            else
              trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-13
            end if
          end do; end do; end do
#endif

        case ('CFC')
#ifdef SHINDELL_STRAT_CHEM
         if(use_rad_cfc.le.0)then
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_CFC
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(I,J,L,n) = CFCIC(I,J,L)*ICfactor
          end do   ; end do   ; end do
         else
           if (is_set_param('initial_GHG_setup')) then
             call get_param('initial_GHG_setup', initial_GHG_setup)
             if (initial_GHG_setup == 1 .and. itime == itimeI) then
               select case(PI_run)
             case(1)     ; ICfactor=PIratio_CFC
               case default; ICfactor=1.d0
             end select
             do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
               trm(I,J,L,n) = CFCIC(I,J,L)*ICfactor
             end do   ; end do   ; end do
           else
             if(ghg_yr/=0)then; write(ghg_name,'(I4)') ghg_yr
             else; write(ghg_name,'(I4)') jyear; endif
             ghg_file='GHG_IC_'//ghg_name
             call openunit(ghg_file,iu_data,.true.,.true.)
             do m=1,5
               CALL READT8_COLUMN(grid,iu_data,NAMEUNIT(iu_data),GHG_IN,
     &              0)
               rad_to_file(m,:,I_0:I_1,J_0:J_1)=
     &              ghg_in(:,I_0:I_1,J_0:J_1)
             enddo
             call closeunit(iu_data)
             do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
               trm(I,J,L,n) = rad_to_file(5,l,i,j)
             end do   ; end do   ; end do
           end if
         endif
       end if
#endif /* SHINDELL_STRAT_CHEM */

        case ('BrONO2','HBr','HOBr')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            if(l.ge.LS1) then
              trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*vol2mass(n)*2.d-13
            else
              trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-13
            end if
          end do; end do; end do

        case ('HOCl')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            if(l.ge.LS1) then
              trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*vol2mass(n)*5.d-11
            else
              trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*vol2mass(n)*1.d-11
            end if
          end do; end do; end do

        case('DMS')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*vol2mass(n)*5.d-13
          end do; end do; end do

        case('MSA', 'SO2', 'SO4', 'SO4_d1', 'SO4_d2', 'SO4_d3',
     *         'N_d1','N_d2','N_d3','NH3','NH4','NO3p',
     *         'BCII', 'BCIA', 'BCB', 'OCII', 'OCIA', 'OCB', 'H2O2_s',
     *         'OCI1', 'OCI2', 'OCI3', 'OCA1','OCA2', 'OCA3', 'OCA4',
     *         'seasalt1', 'seasalt2',
     *         'M_NO3   ','M_NH4   ','M_H2O   ','N_AKK_1 ',
     *         'N_ACC_1 ','M_DD1_SU','N_DD1_1 ',
     *         'M_DS1_SU','M_DS1_DU','N_DS1_1 ','M_DD2_SU','M_DD2_DU',
     *         'N_DD2_1 ','M_DS2_SU','M_DS2_DU','N_DS2_1 ','M_SSA_SU',
     *         'M_OCC_SU','N_OCC_1 ','M_BC1_SU',
     *         'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ','M_BC3_SU',
     *         'M_BC3_BC','N_BC3_1 ','M_DBC_SU','M_DBC_BC','M_DBC_DU',
     *         'N_DBC_1 ','M_BOC_SU','M_BOC_BC','M_BOC_OC','N_BOC_1 ',
     *         'M_BCS_SU','M_BCS_BC','N_BCS_1 ','M_MXX_SU','M_MXX_BC',
     *         'M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ','M_OCS_SU',
     *         'M_OCS_OC','N_OCS_1 ','H2SO4',
     *         'M_AKK_SU','M_ACC_SU','M_DD1_DU',
     *         'M_SSA_SS','M_SSC_SS','M_BC1_BC','M_OCC_OC',
     *         'M_SSS_SS','M_SSS_SU')
          do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
            trm(i,j,l,n) = am(l,i,j)*axyp(i,j)*vol2mass(n)*5.d-14
          end do; end do; end do

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
        CASE('Clay','Silt1','Silt2','Silt3','Silt4',
     &       'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &       'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &       'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &       'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &       'Sil1QuHe','Sil2QuHe','Sil3QuHe')
          ! defaults ok
          hbaij=0D0
          ricntd=0D0
#endif

      end select

C**** Initialise pbl profile if necessary
      if (needtrs(n)) then
        do it=1,4
        do j=J_0,J_1
        do ipbl=1,npbl
#ifdef TRACERS_WATER
          if(tr_wd_TYPE(n).eq.nWATER)THEN
            trabl(ipbl,n,it,:,j) = trinit*qabl(ipbl,it,:,j)
          ELSE
#endif
            trabl(ipbl,n,it,:,j) = trm(:,j,1,n)*byam(1,:,j)*byaxyp(:,j)


#ifdef TRACERS_WATER
          END IF
#endif
        end do
        end do
        end do
      end if

      write(out_line,*) ' Tracer ',trname(n),' initialized at itime='
     *     ,itime
      call write_parallel(trim(out_line))

      end if
      end do
#endif /* TRACERS_ON */
#ifdef TRACERS_OCEAN
C**** Initialise ocean tracers if necessary
      call tracer_ic_ocean
#endif
C****
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
c read in DMS source
      if (OFFLINE_DMS_SS.ne.1) then !initialize interactive DMS (non-AeroCom)
      call openunit('DMS_SEA',iudms,.true.,.true.)
        DMSinput(:,:,:)= 0.d0
       do mm=1,12
       call readt_parallel(grid,iudms,nameunit(iudms),
     *                     DMSinput(:,:,mm),0)
       end do
       call closeunit(iudms)
      else  ! AEROCOM DMS
c these netcdf reads are still latlon-specific.
c will call read_dist_data for cubed sphere compatibility
        status=NF_OPEN('DMS_SEA',NCNOWRIT,ncidu)
        status=NF_INQ_VARID(ncidu,'dms',id1)
        start(1)=i_0
        start(2)=j_0
        start(3)=1
        count(1)=1+(i_1-i_0)
        count(2)=1+(j_1-j_0)
        count(3)=366
        status=NF_GET_VARA_REAL(ncidu,id1,start,count,DMS_AER_nohalo)
        status=NF_CLOSE(ncidu)
        DMS_AER(I_0:I_1,J_0:J_1,:) = DMS_AER_nohalo(I_0:I_1,J_0:J_1,:)
      endif
 901  FORMAT(3X,3(I4),E11.3)
c read in AEROCOM seasalt
      if (OFFLINE_DMS_SS.eq.1.or.OFFLINE_SS.eq.1) then
        status=NF_OPEN('SALT1',NCNOWRIT,ncidu)
        status=NF_INQ_VARID(ncidu,'salt',id1)
        start(1)=i_0
        start(2)=j_0
        start(3)=1
        count(1)=1+(i_1-i_0)
        count(2)=1+(j_1-j_0)
        count(3)=366
        status=NF_GET_VARA_REAL(ncidu,id1,start,count,SS1_AER_nohalo)
        status=NF_CLOSE(ncidu)
        SS1_AER(I_0:I_1,J_0:J_1,:) = SS1_AER_nohalo(I_0:I_1,J_0:J_1,:)

        status=NF_OPEN('SALT2',NCNOWRIT,ncidu)
        status=NF_INQ_VARID(ncidu,'salt',id1)
        start(1)=i_0
        start(2)=j_0
        start(3)=1
        count(1)=1+(i_1-i_0)
        count(2)=1+(j_1-j_0)
        count(3)=366
        status=NF_GET_VARA_REAL(ncidu,id1,start,count,SS2_AER_nohalo)
        status=NF_CLOSE(ncidu)
        SS2_AER(I_0:I_1,J_0:J_1,:) = SS2_AER_nohalo(I_0:I_1,J_0:J_1,:)
      endif

c read in SO2 emissions
c Industrial
      if (imPI.eq.0) then
        if (imAER.eq.0) then
C    Initialize:
          so2_src(:,:,:)= 0.d0
          call openunit('SO2_IND',iuc,.false.,.true.)
          do
            READ(iuc,*) ii,jj,carbstuff
            if (ii.eq.0) exit
            if (jj<j_0 .or. jj>j_1) cycle
            so2_src(ii,jj,1)=carbstuff/(sday*30.4)/12.d0
          end do
          call closeunit(iuc)
        endif
#ifndef FAA_emission
C Put aircraft for so2 and BC
        so2_src_3d(:,:,:,2)= 0.d0
        bci_src_3d(:,:,:)=0.d0
        craft(:,:,:) = 0d0
        if (imAER.ne.1.and.imAER.ne.3.and.imAER.ne.5) then

        if ( AM_I_ROOT() ) then
          write(6,*) 'Reading of the AIRCRAFT file now uses ',
     &         'READT_PARALLEL which expects 80-character titles. ',
     &         'Please expand the 56-character AIRCRAFT titles ',
     &         'and then delete this warning/stop statement.'
        endif
        call stop_model('AIRCRAFT file',255)

        call openunit('AIRCRAFT',iuc2,.true.,.true.)
        if (lm.eq.23) then
          lmax = ls1+1
        else
          lmax = lm
        endif
        DO L=1,LMAX
          CALL READT_PARALLEL(grid,iuc2,NAMEUNIT(iuc2),craft(:,:,l),0)
        END DO
        call closeunit(iuc2)

c craft is Kg fuel/day. Convert to Kg SO2/s. 2.3 factor
c adjusts 2015 source to 1990 source.
c 4.d0d-4 converts to kg S
c (for BC multiply this by 0.1)
        so2_src_3d(:,j_0:j_1,:,2)=craft(:,j_0:j_1,:)*4.0d-4
     *    *tr_mm(n_SO2)/32.d0/2.3d0/sday
        bci_src_3D(:,j_0:j_1,:)=so2_src_3D(:,j_0:j_1,:,2)*0.1d0/2.d0
c divide by 2 because BC = 0.1 x S, not SO2
        endif
#endif
      endif ! imPI eq 0
c volcano - continuous
C    Initialize:
      so2_src_3D(:,:,:,1)= 0.d0
c read lat-lon netcdf file and convert lat,lon,pres to i,j,l.
c NOTE: the input file specifies integrals over its gridboxes.
      ALLOCATE(  psref(grid%i_strt_halo:grid%i_stop_halo,
     &                 grid%j_strt_halo:grid%j_stop_halo) )
      call openunit('PSREF',iu_ps,.true.,.true.)
      CALL READT_PARALLEL(grid,iu_ps,'PSREF',psref,1)
      call closeunit(iu_ps)
      status = nf_open('SO2_VOLCANO',nf_nowrite,file_id)
      status = nf_inq_varid(file_id,'lon',vid)
      status = nf_get_var_double(file_id,vid,volc_lons)
      status = nf_inq_varid(file_id,'lat',vid)
      status = nf_get_var_double(file_id,vid,volc_lats)
      status = nf_inq_varid(file_id,'Pres_CONTmax',vid)
      status = nf_get_var_double(file_id,vid,volc_pup)
      status = nf_inq_varid(file_id,'VOLC_CONT',vid)
      status = nf_get_var_double(file_id,vid,volc_emiss)
      status = nf_close(file_id)
      do jlat=1,ubound(volc_lats,1)
        do ilon=1,ubound(volc_lons,1)
          if(volc_emiss(ilon,jlat) <= 0.) cycle
          call lonlat_to_ij(
     &         (/volc_lons(ilon),volc_lats(jlat)/),volc_ij)
          ii = volc_ij(1); jj = volc_ij(2)
          if(jj<j_0 .or. jj>j_1) cycle
          if(ii<i_0 .or. ii>i_1) cycle
          call CALC_VERT_AMP(psref(ii,jj)-ptop,lm,
     &         x1d,amref,x1d,pednref,x1d)
          lmax = 1
          do while(pednref(lmax) > volc_pup(ilon,jlat))
            lmax = lmax + 1
          enddo
          amsum = sum(amref(1:lmax))
          do ll=1,lmax ! add source between surf and max height
            so2_src_3d(ii,jj,ll,1) = so2_src_3d(ii,jj,ll,1)
     &          +(amref(ll)/amsum)*
     &           volc_emiss(ilon,jlat)/(sday*30.4d0)/12.d0
          enddo
        enddo
      enddo
      deallocate(psref)
c read ASCII file specifying the i,j,l of point sources
c      call openunit('SO2_VOLCANO',iuc,.false.,.true.)
c      do
c      read(iuc,*) ii,jj,ll,carbstuff
c      if (ii.eq.0) exit
c      if (jj<j_0 .or. jj>j_1) cycle
c      so2_src_3d(ii,jj,ll,1)=carbstuff/(sday*30.4d0)/12.d0
c      end do
c      call closeunit(iuc)

c Biomass for AeroCom
      if (imAER.eq.1) then
        so2_biosrc_3D(:,:,:,:)= 0.d0
        call openunit('SO2_BIOMASS',iuc,.false.)
        do
          read(iuc,*) ii,jj,mmm,ll,carbstuff
          if (ii.eq.0) exit
          if (jj<j_0 .or. jj>j_1) cycle
          if (imPI.eq.1) carbstuff=carbstuff*0.5d0
          SO2_biosrc_3D(ii,jj,ll,mmm)=carbstuff/(sday*30.4d0)
        end do
        call closeunit(iuc)

c volcano - explosive, for AeroCom
        volc_exp(:,:,:)= 0.d0
        call openunit('SO2_VOLCANO_EXP',iuc,.false.,.true.)
        do
          read(iuc,*) ii,jj,ll,carbstuff
          if (ii.eq.0) exit
          if (jj<j_0 .or. jj>j_1) cycle
          volc_exp(ii,jj,ll)=carbstuff/(sday*30.4d0)/12.d0
        end do
        call closeunit(iuc)
        so2_src_3d(:,j_0:j_1,:,1)=so2_src_3d(:,j_0:j_1,:,1)+
     *                            volc_exp(:,j_0:j_1,:)
      endif
! ----------------------------------------------------
c Biomass BC/OC for AEROCOM, with their own 3D distribution
c  (otherwise it is done in get_hist_BM, across boundary layer)
      if (imAER.eq.1) then ! AEROCOM
        BCB_src(:,:,:,:)= 0.d0
        OCB_src(:,:,:,:)=0.d0
        call openunit('BC_BIOMASS',iuc,.false.,.true.)
        do
          read(iuc,*) ii,jj,mmm,ll,carbstuff
          if (ii.eq.0) exit
          if (jj<j_0 .or. jj>j_1) cycle
          if (imPI.eq.1) carbstuff=carbstuff*0.5d0
          BCB_src(ii,jj,ll,mmm)=carbstuff
        end do
        call closeunit(iuc)
        call openunit('OC_BIOMASS',iuc,.false.,.true.)
        do
          read(iuc,*) ii,jj,mmm,ll,carbstuff
          if (ii.eq.0.) exit
          if (jj<j_0 .or. jj>j_1) cycle
          if (imPI.eq.1) carbstuff=carbstuff*0.5d0
          OCB_src(ii,jj,ll,mmm)=carbstuff
        end do
        call closeunit(iuc)
        ccnv=1.d0/(sday*30.4)
        BCB_src(:,j_0:j_1,1:7,:)=BCB_src(:,j_0:j_1,1:7,:)*ccnv
        OCB_src(:,j_0:j_1,1:7,:)=OCB_src(:,j_0:j_1,1:7,:)*ccnv  !*1.3d0
      endif
#endif
! ---------------------------------------------------
#ifndef TRACERS_AEROSOLS_SOA
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP)
c Terpenes
      if (imAER.ne.1) then
        OCT_src(:,:,:)=0.d0
        call openunit('Terpenes_01',iuc,.true.,.true.)
        call skip_parallel(iuc)
        do mm=1,12
          call readt_parallel(grid,iuc,nameunit(iuc),OCT_src(:,:,mm),0)
        end do
        call closeunit(iuc)
c units are mg Terpene/m2/month
        do i=I_0,I_1; do j=J_0,J_1; do mm=1,12
! 10% of terpenes end up being SOA

          OCT_src(i,j,mm)=OCT_src(i,j,mm)*axyp(i,j)*0.1d0
#ifdef TRACERS_AMP
     +                    *om2oc(n_M_OCC_OC)
#else
     +                    *om2oc(n_OCII)
#endif
        end do; end do; end do
      else ! AEROCOM
c This assumes 10% emission yield (Chin, Penner)
c 1.3 converts OC to OM
        OCT_src(:,:,:)=0.d0
        call openunit('TERPENE',mon_unit,.false.,.true.)
        do
          read(mon_unit,*) ii,jj,mm,carbstuff
          if (ii.eq.0) exit
          if (jj<j_0 .or. jj>j_1) cycle
          carbstuff=carbstuff/(sday*30.4d0)
          OCT_src(ii,jj,mm)=carbstuff   !*1.3d0
        end do
        call closeunit(mon_unit)
      endif
#endif
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_OM_SP
       OCI_src(:,:,1:8)= 0.d0
       call openunit('OC_FOSSIL_FUEL',iuc,.false.,.true.)
       do
       read(iuc,*) ii,jj,(carb(ir),ir=1,8)
       if (ii.eq.0) exit
       if (jj<j_0 .or. jj>j_1) cycle
       DO mm=1,8
       OCI_src5(ii,jj,mm)=carb(mm)
       END DO
       call closeunit(iuc)
       OCI_src(:,j_0:j_1,1:4)=OCI_src5(:,j_0:j_1,1:4)
     * *1.3d0/(sday*365.d0)
#endif
#if (defined TRACERS_NITRATE) || (defined TRACERS_AMP)
      if (imAER.ne.3.and.imAER.ne.5) then  !else read in daily
c read in NH3 emissions
c Industrial
        NH3_src_con(:,:)=0.d0
        NH3_src_cyc(:,:)=0.d0

        call openunit('NH3SOURCE_CYC',iuc,.false.)
        do
          read(iuc,*) ii,jj,carbstuff
          if (ii.eq.0.) exit
          if (jj<j_0 .or. jj>j_1) cycle
! conversion [kg N/gb/a] -> [kg NH3 /gb/s]
          NH3_src_cyc(ii,jj)=carbstuff*1.2142/(sday*30.4*12.)!/axyp(ii,jj)
        end do
        call closeunit(iuc)

        call openunit('NH3SOURCE_CON',iuc,.false.)
        do
          read(iuc,*) ii,jj,carbstuff
          if (ii.eq.0.) exit
          if (jj<j_0 .or. jj>j_1) cycle
          NH3_src_con(ii,jj)=carbstuff*1.2142/(sday*30.4*12.)!/axyp(ii,jj)
        end do
        call closeunit(iuc)
      endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
c **** reads in files for dust/mineral tracers
      call tracer_ic_soildust
#endif

      end subroutine tracer_IC


      subroutine daily_tracer(end_of_day)
!@sum daily_tracer is called once a day for tracers
!@+   SUBROUTINE tracer_IC is called from daily_tracer to allow for
!@+     tracers that 'turn on' on different dates.
!@auth Jean Lerner
C**** Note this routine must always exist (but can be a dummy routine)
      USE MODEL_COM, only:jmon,jday,itime,coupled_chem,fearth0,focean
     $     ,flake0,jyear,p,t,lm
      USE SOMTQ_COM, only : tmom,mz
      USE DOMAIN_DECOMP_ATM, only : grid, get, write_parallel, am_i_root
      USE RAD_COM, only: o3_yr
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only : variable_phi
#endif
      USE TRACER_COM, only: daily_z
      USE CONSTANT, only: grav
      USE TRACER_COM, only: ntm,trname,itime_tr0,nOther,nAircraft,
     & n_CH4,n_Isoprene,n_codirect,sfc_src,ntsurfsrc,ssname,do_fire,
     & trans_emis_overr_yr,trans_emis_overr_day,nBBsources
#ifdef TRACERS_SPECIAL_Shindell
     & ,ntm_chem
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
     * ,aer_int_yr,imAER,imPI,n_NH3,n_SO2,n_SO4,n_BCB,n_OCB
     * ,n_M_ACC_SU,n_M_AKK_SU,n_M_BC1_BC,n_M_OCC_OC,n_M_BOC_BC
     * ,n_M_BOC_OC
      USE AEROSOL_SOURCES, only: SO2_biosrc_3D
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      USE resolution,ONLY : lm
      USE GEOM, only: axyp
      USE DYNAMICS, only: am  ! Air mass of each box (kg/m^2)
      USE TRACER_COM, only: trm,vol2mass,trmom,n_CO2n
      USE DOMAIN_DECOMP_ATM, only: GRID,GLOBALSUM
      USE MODEL_COM, only : nday,nstep=>itime,psf
#ifdef constCO2
      USE obio_forc, only : atmCO2
#else
      USE RADPAR, only : xnow
#endif
#endif
#ifdef TRACERS_SPECIAL_Lerner
      USE TRACERS_MPchem_COM, only: n_MPtable,tcscale,STRATCHEM_SETUP
      USE LINOZ_CHEM_COM, only: LINOZ_SETUP
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE FLUXES, only: tr3Dsource
      USE TRCHEM_Shindell_COM,only: PI_run, use_rad_ch4, rad_FL,
     & dms_offline,so2_offline,sulfate,fix_CH4_chemistry
#endif
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only : variable_phi
#endif

      IMPLICIT NONE
      INTEGER n,last_month,kk,nread,xday,xyear,ns
      LOGICAL, INTENT(IN) :: end_of_day
#ifdef TRACERS_GASEXCH_ocean_CO2
      integer i,j,l
      real*8 :: sarea,
     &          trm_vert(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                   GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     &          trm_glbavg,factor,atm_glbavg
#endif
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
     &     :: daily_gz
      data last_month/-1/
      INTEGER J_0, J_1, I_0, I_1
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      if(end_of_day) then
        call COMPUTE_GZ(p,t,tmom(mz,:,:,:),daily_z)
        daily_z = daily_z/grav
      endif
      daily_gz = grav*daily_z

#ifdef TRACERS_SPECIAL_Lerner
      if (.not. end_of_day) then
C**** Initialize tables for linoz
      do n=1,ntm
        if (trname(n).eq."O3" .and. itime.ge.itime_tr0(n)) then
          call linoz_setup(n)
          exit
        end if
      end do

C**** Initialize tables for Prather StratChem tracers
      do n=1,ntm
        if (trname(n).eq."N2O" .or. trname(n).eq."CH4" .or.
     *      trname(n).eq."CFC11")
     *    call stratchem_setup(n_MPtable(n),trname(n))
      end do
      end if  ! not end of day

C**** Prather StratChem tracers and linoz tables change each month
      IF (JMON.NE.last_month) THEN
        do n=1,ntm
          if ((trname(n).eq."N2O" .or. trname(n).eq."CH4" .or.
     *         trname(n).eq."CFC11") .and. itime.ge.itime_tr0(n)) then
            CALL STRTL  ! one call does all based on n_MPtable_max
            exit
          end if
        end do
        do n=1,ntm
          if (trname(n).eq."O3" .and. itime.ge.itime_tr0(n)) then
            CALL linoz_STRATL
            exit
          end if
        end do
        last_month = JMON
      END IF

C**** Tracer specific call for CO2
      do n=1,ntm
        if (trname(n).eq."CO2") then
          call read_CO2_sources(n)
          exit
        end if
      end do

C**** Tracer specific call for CH4
      do n=1,ntm
        if (trname(n).eq."CH4") then
          call read_CH4_sources(n)
          exit
        end if
      end do
#endif

#ifdef TRACERS_COSMO
      if (variable_phi .eq. 0) then
         call read_Be_source_noAlpha
         print*, "called old version of Be source"
      end if

      if (variable_phi .eq. 1) then
         call read_Be_source
         print*, "called new version of Be source"
      end if

      if (variable_phi .eq. 2) then
         if ((jday .eq. 1) .or. (.not. end_of_day)) then
            call update_annual_phi
            print*, "called update_annual_phi"
         end if
      end if

      if (variable_phi .eq. 3) then
         call update_daily_phi
         print*, "called update_daily_phi"
      end if
#endif


#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) \
 || (defined TRACERS_AMP)
C**** Allow overriding of transient emissions date:
! for now, tying this to O3_yr becasue Gavin
! didn't want a new parameter, also not allowing
! day overriding yet, because of that.
      trans_emis_overr_yr=ABS(o3_yr)
      if(trans_emis_overr_yr > 0)then
        xyear=trans_emis_overr_yr
      else
        xyear=jyear
      endif
!!    if(trans_emis_overr_day > 0)then
!!      xday=trans_emis_overr_day
!!    else
        xday=jday
!!    endif
#ifdef TRACERS_SPECIAL_Shindell
C**** Next line for fastj photon fluxes to vary with time:
      if(rad_FL.gt.0) call READ_FL(end_of_day)
C**** Daily tracer-specific calls to read 2D and 3D sources:
      if (COUPLED_CHEM.ne.1) then
        call read_aero(dms_offline,'DMS_FIELD') !not applied directly to tracer
        call read_aero(so2_offline,'SO2_FIELD') !not applied directly to tracer
      endif
#ifdef CUBED_SPHERE
      call get_aircraft_NOx(xyear,xday,daily_gz,.true.)
#endif
#endif /* TRACERS_SPECIAL_Shindell */

      do n=1,ntm
        if(trname(n)=='CH4')then ! ---------- methane --------------
#ifdef TRACERS_SPECIAL_Shindell
         nread=ntsurfsrc(n)+nBBsources(n)
         if(nread>0) call read_sfc_sources(n,nread,xyear,xday,.true.)
#ifdef WATER_MISC_GRND_CH4_SRC
         do ns=1,ntsurfsrc(n) 
           if(ssname(n,ns)=='gsfMGOLjal')sfc_src(I_0:I_1,J_0:J_1,n,ns)=
     &     1.698d-12*fearth0(I_0:I_1,J_0:J_1) + ! 5.3558e-5 Jean
     &     5.495d-11*flake0(I_0:I_1,J_0:J_1)  + ! 17.330e-4 Jean
     &     1.141d-12*focean(I_0:I_1,J_0:J_1)    ! 3.5997e-5 Jean
         end do
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
         if(nread>0) call read_ncep_for_wetlands(end_of_day)
#endif
#ifdef GFED_3D_BIOMASS
         if(fix_CH4_chemistry/=1.or.ntsurfsrc(n)/=0)
     &   call get_GFED_biomass_burning(n,xyear,xday)
#endif
#endif /* TRACERS_SPECIAL_Shindell */
        else !-------------------------------------- general ---------
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
          if ( ! this if statement is needed, since some tracers have ntsurfsrc>0 but no trname_XX files. It should dissappear one day.
#ifdef TRACERS_SPECIAL_Shindell
     &        n<=ntm_chem .or.
#endif
     &        n==n_SO2
#ifdef TRACERS_AEROSOLS_Koch
     &        .or. n==n_BCB .or. n==n_OCB
#endif
#if (defined TRACERS_NITRATE) || (defined TRACERS_AMP)
     &        .or. n==n_NH3
#endif
#ifdef TRACERS_AMP
     &        .or. n==n_M_BC1_BC .or. n==n_M_OCC_OC
     &        .or. n==n_M_BOC_BC .or. n==n_M_BOC_OC
#endif
     &        ) then
#ifdef TRACERS_SPECIAL_Shindell
            if (n>ntm_chem) then
#endif
              if(aer_int_yr > 0) then
                xyear=aer_int_yr
              else
                xyear=jyear
              endif
#ifdef TRACERS_SPECIAL_Shindell
            endif
#endif
#endif /* (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) */

            nread=ntsurfsrc(n) ! default
            select case (trname(n))
            case ('Alkenes', 'CO', 'NOx', 'Paraffin', ! CH4 done above
     &      'NH3', 'SO2', 'BCB', 'OCB', ! do not include sulfate here
     &      'M_BC1_BC', 'M_OCC_OC', 'M_BOC_BC', 'M_BOC_OC')
              nread=nread+nBBsources(n)
            end select

#ifndef TRACERS_AEROSOLS_SOA
#ifdef TRACERS_AMP
            select case (trname(n))
            case ('M_OCC_OC')
              nread=nread-1
            end select
#endif
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AMP
            if(nread>0)call read_sfc_sources(n,nread,xyear,xday,.false.)
#ifndef TRACERS_AEROSOLS_SOA
            select case (trname(n))
            case ('M_OCC_OC')
              sfc_src(:,J_0:J_1,n,ntsurfsrc(n):
     &                            ntsurfsrc(n)+nBBsources(n))=
     &          sfc_src(:,J_0:J_1,n,ntsurfsrc(n)-1:
     &                              ntsurfsrc(n)+nBBsources(n)-1)
              sfc_src(:,J_0:J_1,n,ntsurfsrc(n))=0.d0 ! this will become terpene sources
            end select
#endif  /* TRACERS_AEROSOLS_SOA */
#else
            if(nread>0)call read_sfc_sources(n,nread,xyear,xday,.true.)
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
            if (trname(n).eq.'SO2') then ! set this AFTER reading
#ifdef TRACERS_AEROSOLS_Koch
              ntsurfsrc(n_SO4)=ntsurfsrc(n)
#endif
#ifdef TRACERS_AMP
              ntsurfsrc(n_M_ACC_SU)=ntsurfsrc(n)
#ifndef TRACERS_AMP_M4
              ntsurfsrc(n_M_AKK_SU)=ntsurfsrc(n)
#endif
#endif
            endif
#endif
#ifdef TRACERS_SPECIAL_Shindell
            select case (trname(n))
            case ('NOx')
!           (lightning and aircraft called from tracer_3Dsource)
            case ('N2O5')
              tr3Dsource(I_0:I_1,J_0:J_1,:,:,n) = 0.
              if (COUPLED_CHEM.ne.1)
     &        call read_aero(sulfate,'SULFATE_SA') !not applied directly
            end select
#ifdef GFED_3D_BIOMASS
            select case (trname(n))
            case ('NOx','CO','Alkenes','Paraffin')
              call get_GFED_biomass_burning(n,xyear,xday)
            end select
#endif
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
          endif ! n=n_...
#endif
        endif !------------------------------------------------------
      end do ! ntm

#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      call read_sfc_sources(n_codirect,ntsurfsrc(n_codirect),xyear,
     & xday,.true.)
#endif

#endif /* TRACERS_SPECIAL_Shindell || TRACERS_AEROSOLS_Koch || TRACERS_AMP) */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
        if (imAER.eq.3.or.imAER.eq.5) then
          do n=1,ntm
            select case (trname(n))
            case ('SO2')
              if (imPI.eq.0) then    ! these are all pollution
c               if (imAER.eq.3) call get_ships(end_of_day,jday)  !reads in SO2, BC and POM sources
c               if (imAER.eq.3) call get_aircraft_SO2   ! this does SO2 and BCIA
              endif
            case ('BCII')  !this does BC and OC
              if (imPI.eq.0) call get_BCOC(end_of_day,jday)
            case ('BCB', 'M_BOC_BC')  ! this does BCB and OCB and SO2
              if (imAER==3) call get_hist_BMB(end_of_day,jday)
            case ('NH3')
              if (imAER==3) call read_hist_NH3(end_of_day,jday)
            end select
          end do
        endif
c      if (COUPLED_CHEM.ne.1) call get_O3_offline
#endif
#ifdef TRACERS_OM_SP
        call get_hist_BMB(end_of_day)
#endif

C US FAA aviation emissions (must have chemistry and aerosols switched on) nu/mk
#ifdef FAA_emission 
         call read_daily_faa
#endif

C****
C**** Initialize tracers here to allow for tracers that 'turn on'
C**** at the start of any day

      call tracer_IC


!gas exchange CO2 case reset trm here
!for the constCO2 case just reset to atmCO2 which is defined in the rundeck
!for the variable case (presently default) reset to the value scaled by
!the xnow value.
#ifdef TRACERS_GASEXCH_ocean_CO2

      if (end_of_day) then ! only at end of day
         n=n_CO2n

         !area weighted tracer global average
         do j=J_0,J_1 ; do i=I_0,I_1
             trm_vert(i,j) = sum(trm(i,j,1:lm,n))
         enddo; enddo

         CALL GLOBALSUM(grid,axyp,    sarea,     all=.true.)
         CALL GLOBALSUM(grid,trm_vert,trm_glbavg,all=.true.)

         !total atm mass
         atm_glbavg = PSF*sarea*100.d0/grav

#ifdef constCO2
         !current concentration to new concentration
         factor = atmCO2*atm_glbavg/trm_glbavg *vol2mass(n)*1.d-6

         if(AM_I_ROOT( ))then
         write(*,'(a,i5,8e12.4)')
     .           "TRACER_DRV, factor", nstep,factor,
     .            atm_glbavg,vol2mass(n),
     .            atmCO2,sarea,
     .            PSF,grav,
     .            trm_glbavg/(atm_glbavg*vol2mass(n)*1.d-6)
         endif
#else
         !current concentration to new concentration
         factor = xnow(1)*atm_glbavg/trm_glbavg *vol2mass(n)*1.d-6

         if(AM_I_ROOT( ))then
         write(*,'(a,i5,8e12.4)')
     .           "TRACER_DRV, factor", nstep,factor,
     .            atm_glbavg,vol2mass(n),
     .            xnow(1),sarea,
     .            PSF,grav,
     .            trm_glbavg/(atm_glbavg*vol2mass(n)*1.d-6)
         endif
#endif
         do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
             trm(i,j,l,n) = factor*trm(i,j,l,n)
         enddo; enddo; enddo

         if (factor .lt. 1.d0) then ! adjust moments
           do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
             trmom(:,i,j,l,n)=factor*trmom(:,i,j,l,n)
             !?? do we need trmom=0 for the atmco2 case?
           enddo; end do; end do
         end if

      endif  ! end_of_day
#endif

      return
      end subroutine daily_tracer

#ifdef TRACERS_ON
      SUBROUTINE set_tracer_2Dsource
!@sum tracer_source calculates non-interactive sources for tracers
!@auth Jean Lerner/Gavin Schmidt
      USE MODEL_COM, only: itime,JDperY,fland,psf,pmtop,jmpery
     *  ,dtsrc,jmon,nday,flice,focean,im,jm
      USE DOMAIN_DECOMP_ATM, only : GRID, GET, GLOBALSUM,AM_I_ROOT
     *   ,globalmax

      USE GEOM, only: axyp,areag,lat2d_dg,lon2d_dg,imaxj,lat2d
      USE QUSDEF
      USE DYNAMICS, only: am  ! Air mass of each box (kg/m^2)
      USE TRACER_COM
      USE FLUXES, only: trsource
      USE SEAICE_COM, only: rsi
      USE GHY_COM, only : fearth
      USE CONSTANT, only: tf,sday,hrday,bygrav,mair,pi,teeny
      USE PBLCOM, only: tsavg
#if (defined INTERACTIVE_WETLANDS_CH4) && (defined TRACERS_SPECIAL_Shindell)
      USE TRACER_SOURCES, only: ns_wet,add_wet_src
#endif
#ifdef TRACERS_SPECIAL_Lerner
      USE CO2_SOURCES, only: co2_src
      USE CH4_SOURCES, only: ch4_src
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP)
      USE AEROSOL_SOURCES, only: so2_src,BCI_src,OCI_src,BBinc,
#ifndef TRACERS_AEROSOLS_SOA
     * OCT_src,
#endif  /* TRACERS_AEROSOLS_SOA */
     * NH3_src_con,NH3_src_cyc,
     * BC_ship,POM_ship,SO2_ship
#endif
#ifdef TRACERS_RADON
      USE AEROSOL_SOURCES, only: rn_src
#endif
#if (defined TRACERS_NITRATE) || (defined TRACERS_AMP) || \
    (defined TRACERS_SPECIAL_Shindell)
      USE apply3d, only : apply_tracer_3Dsource
      USE RAD_COM,  only : cosz1,cosz_day
#endif
#ifdef TRACERS_AMP
      USE AERO_SETUP, only : RECIP_PART_MASS
      USE TRDIAG_COM, only : taijs=>taijs_loc,ijts_AMPe
#endif
      USE LAKES_COM, only : flake
      implicit none
      integer :: i,j,ns,ns_isop,l,ky,n,nsect,kreg
      REAL*8 :: source,sarea,steppy,base,steppd,x,airm,anngas,
     *  tmon,bydt,tnew,scca(im),fice
      REAL*8 :: sarea_prt(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                    GRID%J_STRT_HALO:GRID%J_STOP_HALO)
#ifdef TRACERS_SPECIAL_Shindell
c      real*8 :: factj(GRID%J_STRT_HALO:GRID%J_STOP_HALO)
c      real*8 :: nlight, max_COSZ1, fact0
#endif
      real*8 :: lon_w,lon_e,lat_s,lat_n
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      integer :: src_index
      real*8 :: src_fact
#endif

#ifdef TRACERS_TERP
!@param orvoc_fact Fraction of ORVOC added to Terpenes, for SOA production (Griffin et al., 1999)
      real*8, parameter :: orvoc_fact=0.32d0
      real*8 :: max_isop_flux
#endif  /* TRACERS_TERP */
#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_GASEXCH_ocean_CFC)
      integer :: i_ocmip,imax
      real*8  :: factor
      real*8  :: trsource_prt(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                        GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      real*8, dimension(ntm) :: trsource_glbavg
#endif
      INTEGER I_0, I_1, J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1, I_STRT=I_0, I_STOP=I_1)

      bydt = 1./DTsrc
C**** All sources are saved as kg/s
      do n=1,ntm
      if (itime.lt.itime_tr0(n)) cycle
      select case (trname(n))

      case default
!     write(6,*) ' Sources for ',trname(n),' are not in this routine'
C****
C**** Surface Sources of SF6 and CFCn (Same grid as CFC11)
C****
      case ('SF6','CFC11','CFCn','SF6_c')
        trsource(:,:,:,n)=0
C**** SF6 source increases each year by .3pptv/year
C**** SF6_c source is constant, same as first year SF6, but always
C**** CFCn source increases each year so that the glbavg is from obs
C**** CFC source is the same each year
C**** Distribute source over ice-free land
        steppy = 1./(sday*JDperY)
        if (trname(n).eq.'SF6' .or. trname(n).eq.'CFCn' .or.
     *      trname(n).eq.'SF6_c') then
C         Make sure index KY=1 in year that tracer turns on
          ky = 1 + (itime-itime_tr0(n))/(nday*JDperY)
          if (trname(n).eq.'SF6_c') ky = 1
          base = (0.3d-12)*vol2mass(n) !pptm
          x = base*ky
          airm = (psf-pmtop)*100.*bygrav*AREAG !(kg/m**2 X m**2 = kg)
          anngas = x*airm
        else if (trname(n).eq.'CFC11') then
          anngas = 310.d6
        endif

c Could the masks for latlon rectangles be precomputed at
c initialization and their areas saved? (Or do whenever
c fearth changes.)

C**** Source over United States and Canada
        source = .37d0*anngas*steppy
        lon_e =  -70.d0
        lon_w = -125.d0
        lat_n =   50.d0
        lat_s =   30.d0
        call get_latlon_mask(lon_w,lon_e,lat_s,lat_n,sarea_prt)
        do j=j_0,j_1; do i=i_0,i_1
            sarea_prt(i,j) = sarea_prt(i,j)*axyp(i,j)*fearth(i,j)
        enddo; enddo
        call globalsum(grid, sarea_prt, sarea, all=.true.)
        do j=j_0,j_1; do i=i_0,i_1
            trsource(i,j,1,n) = trsource(i,j,1,n) +
     &         source*sarea_prt(i,j)/sarea
        enddo; enddo
C**** Source over Europe and Russia
        source = .37d0*anngas*steppy
        lon_e =  45.d0
        lon_w = -10.d0
        lat_n =  65.d0
        lat_s =  36.1d0 ! 0.1 deg offset avoids overlap with Middle East
        call get_latlon_mask(lon_w,lon_e,lat_s,lat_n,sarea_prt)
        do j=j_0,j_1; do i=i_0,i_1
            sarea_prt(i,j) = sarea_prt(i,j)*axyp(i,j)*fearth(i,j)
        enddo; enddo
        call globalsum(grid, sarea_prt, sarea, all=.true.)
        do j=j_0,j_1; do i=i_0,i_1
            trsource(i,j,1,n) = trsource(i,j,1,n) +
     &         source*sarea_prt(i,j)/sarea
        enddo; enddo
C**** Source over Far East
        source = .13d0*anngas*steppy
        lon_e = 150.d0
        lon_w = 120.d0
        lat_n =  45.d0
        lat_s =  20.d0
        call get_latlon_mask(lon_w,lon_e,lat_s,lat_n,sarea_prt)
        do j=j_0,j_1; do i=i_0,i_1
            sarea_prt(i,j) = sarea_prt(i,j)*axyp(i,j)*fearth(i,j)
        enddo; enddo
        call globalsum(grid, sarea_prt, sarea, all=.true.)
        do j=j_0,j_1; do i=i_0,i_1
            trsource(i,j,1,n) = trsource(i,j,1,n) +
     &         source*sarea_prt(i,j)/sarea
        enddo; enddo
C**** Source over Middle East
        source = .05d0*anngas*steppy
        lon_e = 75.d0
        lon_w = 30.d0
        lat_n = 35.9d0 ! 0.1 deg offset avoids overlap with Europe
        lat_s = 15.d0
        call get_latlon_mask(lon_w,lon_e,lat_s,lat_n,sarea_prt)
        do j=j_0,j_1; do i=i_0,i_1
            sarea_prt(i,j) = sarea_prt(i,j)*axyp(i,j)*fearth(i,j)
        enddo; enddo
        call globalsum(grid, sarea_prt, sarea, all=.true.)
        do j=j_0,j_1; do i=i_0,i_1
            trsource(i,j,1,n) = trsource(i,j,1,n) +
     &         source*sarea_prt(i,j)/sarea
        enddo; enddo
C**** Source over South America
        source = .04d0*anngas*steppy
        lon_e = -40.d0
        lon_w = -50.d0
        lat_n = -22.5d0
        lat_s = -23.5d0
        call get_latlon_mask(lon_w,lon_e,lat_s,lat_n,sarea_prt)
        do j=j_0,j_1; do i=i_0,i_1
            sarea_prt(i,j) = sarea_prt(i,j)*axyp(i,j)*fearth(i,j)
        enddo; enddo
        call globalsum(grid, sarea_prt, sarea, all=.true.)
        do j=j_0,j_1; do i=i_0,i_1
            trsource(i,j,1,n) = trsource(i,j,1,n) +
     &         source*sarea_prt(i,j)/sarea
        enddo; enddo
C**** Source over South Africa
        source = .02d0*anngas*steppy
        lat_n = -24.d0
        lat_s = -28.d0
        lon_e =  30.d0
        lon_w =  25.d0
        call get_latlon_mask(lon_w,lon_e,lat_s,lat_n,sarea_prt)
        do j=j_0,j_1; do i=i_0,i_1
            sarea_prt(i,j) = sarea_prt(i,j)*axyp(i,j)*fearth(i,j)
        enddo; enddo
        call globalsum(grid, sarea_prt, sarea, all=.true.)
        do j=j_0,j_1; do i=i_0,i_1
            trsource(i,j,1,n) = trsource(i,j,1,n) +
     &         source*sarea_prt(i,j)/sarea
        enddo; enddo
C**** Source over Australia and New Zealand
        source = .02d0*anngas*steppy
        lat_n = -33.5d0
        lat_s = -34.5d0
        lon_e = 150.5d0
        lon_w = 149.5d0
        call get_latlon_mask(lon_w,lon_e,lat_s,lat_n,sarea_prt)
        do j=j_0,j_1; do i=i_0,i_1
            sarea_prt(i,j) = sarea_prt(i,j)*axyp(i,j)*fearth(i,j)
        enddo; enddo
        call globalsum(grid, sarea_prt, sarea, all=.true.)
        do j=j_0,j_1; do i=i_0,i_1
            trsource(i,j,1,n) = trsource(i,j,1,n) +
     &         source*sarea_prt(i,j)/sarea
        enddo; enddo

#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_GASEXCH_ocean_CFC)
        !print out global average for each time step before weighing
        !in the OCMIP values
        sarea  = 0.
        trsource_glbavg(n)=0.
        sarea_prt(:,:)  = 0.
        trsource_prt(:,:) = 0.
        do j=J_0,J_1
         imax=72
         if (j .eq. 1 .or. j .eq. 46) imax=1
          do i=I_0,I_1
           factor = axyp(i,j)*fearth(i,j)
           sarea_prt(i,j)= FACTOR
           trsource_prt(i,j) = trsource(i,j,1,n)*FACTOR
          enddo
        enddo

        CALL GLOBALSUM(grid, sarea_prt,    sarea, all=.true.)
        CALL GLOBALSUM(grid,trsource_prt,trsource_glbavg(n),all=.true.)

        trsource_glbavg(n)=trsource_glbavg(n)/sarea

        !weight trsource by ocmip_cfc global average
        !number of steps/year=JDperY*sday/dtsrc=365*86400/1800=17520
        i_ocmip=(itime-itime_tr0(n))/JDperY/int(sday/dtsrc)+1
        if (mod(itime,JDperY*int(sday/dtsrc)) .eq. 0.)
     .     write(6,'(a,2i5)'),
     .             'TRACERS_DRV, new year: itime, i_ocmip=',
     .             itime,i_ocmip
        do j=J_0,J_1 ! TNL
          do i=1,72

cdiag     write(6,'(a,2i5,2e12.4,i5,4e12.4)')'TRACERS_DRV '
cdiag.         ,i,j,trsource(i,j,1,n),ocmip_cfc(i_ocmip,n),
cdiag.          JDperY,hrday,dtsrc,3600,trsource_glbavg(n)

             trsource(i,j,1,n) = trsource(i,j,1,n)
     .        * (ocmip_cfc(i_ocmip,n)/(JDperY*sday/dtsrc))
     .        / trsource_glbavg(n)
          enddo
        enddo

        write(6,'(a,3i5,3e12.4)')'TRACERS_DRV, CFC source at pt(15,33)'
     .               ,15,33,i_ocmip,trsource(15,33,1,n)
     .               ,ocmip_cfc(i_ocmip,n),trsource_glbavg(n)

        !recompute global average after weighting in OCMIP
        sarea  = 0.
        trsource_glbavg(n)=0.
        sarea_prt(:,:)  = 0.
        trsource_prt(:,:) = 0.
        do j=J_0,J_1
         imax=72
         if (j .eq. 1 .or. j .eq. 46) imax=1
          do i=I_0,I_1
           factor = axyp(i,j)*fearth(i,j)
           sarea_prt(i,j)= FACTOR
           trsource_prt(i,j) = trsource(i,j,1,n)*FACTOR
          enddo
        enddo

        CALL GLOBALSUM(grid, sarea_prt,    sarea, all=.true.)
        CALL GLOBALSUM(grid,trsource_prt,trsource_glbavg(n),all=.true.)

        trsource_glbavg(n)=trsource_glbavg(n)/sarea

#endif

C****
C**** Surface Sources for Radon-222
C****
      case ('Rn222')
        trsource(:,J_0:J_1,:,n)=0
C**** ground source
        steppd = 1./sday
        do j=J_0,J_1
          do i=I_0,I_1
          if (rnsrc.eq.0) then !standard source
C**** source from ice-free land
            if(tsavg(i,j).lt.tf) then !composite surface air temperature
              trsource(i,j,1,n) = 1.0d-16*steppd*axyp(i,j)*fearth(i,j)
            else  ! 1 atom/cm^2/s
              trsource(i,j,1,n) = 3.2d-16*steppd*axyp(i,j)*fearth(i,j)
            end if
         else if (rnsrc.eq.1) then !Conen and Robertson
              trsource(i,j,1,n) = 3.2d-16*steppd*axyp(i,j)*fearth(i,j)
c add code to implement Conen and Robertson - linear decrease in Rn222
c   emission from 1 at 30N to 0.2 at 70N and 0.2 north of 70N
           if (nint(lat2d_dg(i,j)).gt.30 .and.
     &         nint(lat2d_dg(i,j)).lt.70) then
             trsource(i,j,1,n)=trsource(i,j,1,n)*
     &            (1.d0-(lat2d_dg(i,j)-30.d0)/40.d0*0.8d0)
           else if (nint(lat2d_dg(i,j)).ge.70) then
             trsource(i,j,1,n)=0.2*trsource(i,j,1,n)
           endif
          else if (rnsrc.eq.2) then !Schery and Wasiolek
#ifdef TRACERS_RADON
c Schery source
          trsource(i,j,1,n)=rn_src(i,j,jmon)
#endif
          endif
          if (rnsrc.le.1) then
C**** source from ice-free ocean
            trsource(i,j,1,n) =trsource(i,j,1,n)+ 1.6d-18*steppd*axyp(i
     *           ,j)*(1.-fland(i,j))*(1.-rsi(i,j))
          endif
          enddo                 !i
        enddo                   !j

#ifdef TRACERS_SPECIAL_Lerner
C****
C**** Sources and sinks for CO2 (kg/s)
C****
      case ('CO2')
        do ns=1,ntsurfsrc(n)
          do j=J_0,J_1
            trsource(:,j,ns,n) = co2_src(:,j,ns)*axyp(:,j)
          end do
        end do

C****
C**** Sources and sinks for CH4 (kg/s)
C****
      case ('CH4')
        do ns=1,ntsurfsrc(n)
          do j=J_0,J_1
            trsource(:,j,ns,n) = ch4_src(:,j,ns)*axyp(:,j)
          end do
        end do
C****
C**** Sources and sinks for N2O:
C**** First layer is set to a constant 462.2 ppbm. (300 PPB V)
C****
      case ('N2O')
      do j=J_0,J_1
        trsource(:,j,1,n) = (am(1,:,j)*axyp(:,j)*462.2d-9
     *   -trm(:,j,1,n))*bydt
      end do
C****
C**** Linoz Deposition from layer 1
C****
      case ('O3')
      call linoz_depo(1,n)
#endif
C****
C**** Sources and sinks for 14CO2
C**** NOTE: This tracer is supposed to start on 10/16
C**** Decay is a function of the number of months since itime_tr0
C**** The tracer is reset to specific values in layer 1 only if
C****   this results in a sink
C****
      case ('14CO2')
      tmon = (itime-itime_tr0(n))*jmpery/(nday*jdpery)
      trsource(:,J_0:J_1,1,n) = 0.
      do j=J_0,J_1
      do i=I_0,I_1
         if (lat2d(i,j).lt.0.) then
               tnew = am(1,i,j)*axyp(i,j)*(4.82d-18*46./mair)*
     *          (44.5 + tmon*(1.02535d0 - tmon*
     *                  (2.13565d-2 - tmon*8.61853d-5)))
               if (tnew.lt.trm(i,j,1,n))
     *             trsource(i,j,1,n) = (tnew-trm(i,j,1,n))*bydt
         else
               tnew = am(1,i,j)*axyp(i,j)*(4.82d-18*46./mair)*
     *          (73.0 - tmon*(0.27823d0 + tmon*
     *                  (3.45648d-3 - tmon*4.21159d-5)))
               if (tnew.lt.trm(i,j,1,n))
     *             trsource(i,j,1,n) = (tnew-trm(i,j,1,n))*bydt
         endif
      end do
      end do

C****
C**** No non-interactive surface sources of Water
C****
      case ('Water')
        trsource(:,J_0:J_1,:,n)=0.d0

#ifdef TRACERS_SPECIAL_Shindell
      case ('Ox','NOx','ClOx','BrOx','N2O5','HNO3','H2O2','CH3OOH',
     &      'HCHO','HO2NO2','CO','PAN','AlkylNit','Alkenes','Paraffin',
     &      'HCl','HOCl','ClONO2','HBr','HOBr','BrONO2','N2O','CFC',
     &      'stratOx','codirect')
#ifdef DYNAMIC_BIOMASS_BURNING
        if(do_fire(n))call dynamic_biomass_burning(n,ntsurfsrc(n)+1)
#endif
        do ns=1,ntsurfsrc(n); do j=J_0,J_1
          trsource(I_0:I_1,j,ns,n)=
     &    sfc_src(I_0:I_1,j,n,ns)*axyp(I_0:I_1,j)
        end do ; end do
#ifdef TRACERS_TERP
      case ('Terpenes')
        do ns=1,ntsurfsrc(n)
          if(ns==1) then
            do j=J_0,J_1
              trsource(I_0:I_1,j,ns,n)=
     &        sfc_src(I_0:I_1,j,n,ns)*axyp(I_0:I_1,j)
            end do
! If no orvoc file provided, scale up the terpenes one instead.
! 0.4371 is the ratio of orvoc/isoprene emissions in the Lathiere et al. (2005) results
            if(ntsurfsrc(n)==1) then ! no orvoc file exists
              call globalmax(grid,
     &             maxval(sfc_src(I_0:I_1,J_0:J_1,n_Isoprene,
     &                            1:ntsurfsrc(n_Isoprene))),
     &             max_isop_flux)
              if (max_isop_flux <= 0.d0) call stop_model(
     &          'Offline isoprene sources are needed', 255)
              do ns_isop=1,ntsurfsrc(n_Isoprene) ! use all Isoprene sources for orvoc scaling
                do j=J_0,J_1
                  trsource(I_0:I_1,j,ns,n)=trsource(I_0:I_1,j,ns,n)+
     &            orvoc_fact*0.4371*axyp(I_0:I_1,j)*
     &            sfc_src(I_0:I_1,j,n_Isoprene,ns_isop)
                end do
              end do
            end if
          else ! use the orvoc file
            do j=J_0,J_1
              trsource(I_0:I_1,j,ns,n)=orvoc_fact*
     &        sfc_src(I_0:I_1,j,n,ns)*axyp(I_0:I_1,j)
            end do
          endif
        end do
#endif  /* TRACERS_TERP */
      case ('CH4')
#ifdef DYNAMIC_BIOMASS_BURNING
        if(do_fire(n))call dynamic_biomass_burning(n,ntsurfsrc(n)+1)
#endif
        do ns=1,ntsurfsrc(n); do j=J_0,J_1
          trsource(I_0:I_1,j,ns,n)=
     &    sfc_src(I_0:I_1,j,n,ns)*axyp(I_0:I_1,j)
        end do ; end do
#ifdef INTERACTIVE_WETLANDS_CH4
        if(ntsurfsrc(n) > 0) then
          call alter_wetlands_source(n,ns_wet)
          do j=J_0,J_1
            trsource(I_0:I_1,j,ns_wet,n)=trsource(I_0:I_1,j,ns_wet,n)+
     &      add_wet_src(I_0:I_1,j)*axyp(I_0:I_1,j)
          enddo
        endif
#endif
#if !defined(PS_BVOC) && !defined(BIOGENIC_EMISSIONS)
      case ('Isoprene')
! Isoprene sources to be emitted only during sunlight, and
! weighted by cos of solar zenith angle:
        do ns=1,ntsurfsrc(n); do j=J_0,J_1; do i=I_0,I_1
          if(COSZ1(i,j)>0.)then
            trsource(i,j,ns,n)=(COSZ1(i,j)/(COSZ_day(i,j)+teeny))*
     &      sfc_src(i,j,n,ns)*axyp(i,j)
          else
            trsource(i,j,ns,n)=0.d0
          endif
        end do ; end do; enddo
#endif
#endif /* TRACERS_SPECIAL_Shindell */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      case ('SO2', 'SO4', 'M_ACC_SU', 'M_AKK_SU',
     &      'BCB', 'OCB',
     &      'M_BC1_BC', 'M_OCC_OC', 'M_BOC_BC', 'M_BOC_OC')
        src_fact=1.d0               ! factor to multiply emissions with
        src_index=n                 ! index to be used for emissions
        select case (trname(n))
        case ('SO2')
          src_fact=0.975d0 ! the rest goes to sulfate (SO4 or M_ACC_SU)
        case ('SO4')
          src_fact=0.0375d0 ! (1.-SO2 fraction)*tr_mm(n_SO4)/tr_mm(n_SO4)
          src_index=n_SO2
        case ('M_ACC_SU')
          src_fact=0.0375d0
#ifndef TRACERS_AMP_M4
     &            *0.99d0 ! the rest goes to M_AKK_SU
#endif
          src_index=n_SO2
#ifndef TRACERS_AMP_M4
        case ('M_AKK_SU')
          src_fact=0.0375d0
     &            *0.01d0
          src_index=n_SO2
#endif
        case ('OCB', 'M_OCC_OC', 'M_BOC_OC')
          src_fact=BBinc
        end select
#ifdef DYNAMIC_BIOMASS_BURNING
        if(do_fire(n))call dynamic_biomass_burning(n,ntsurfsrc(n)+1) 
#endif
#ifndef TRACERS_AEROSOLS_SOA
#ifdef TRACERS_AMP
        select case (trname(n))
        case ('M_OCC_OC')
          sfc_src(:,J_0:J_1,src_index,ntsurfsrc(n))=
     &      OCT_src(:,J_0:J_1,jmon)/axyp(:,J_0:J_1)/src_fact
        end select
#endif
#endif  /* TRACERS_AEROSOLS_SOA */
        do ns=1,ntsurfsrc(src_index)
          trsource(:,J_0:J_1,ns,n)=
     &      sfc_src(:,J_0:J_1,src_index,ns)
     &      *axyp(:,J_0:J_1)*src_fact
        enddo ! ns
#endif /* (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) */
#if (defined TRACERS_NITRATE) || (defined TRACERS_AMP)
      case ('NH3')
#ifdef DYNAMIC_BIOMASS_BURNING
        if(do_fire(n))call dynamic_biomass_burning(n,ntsurfsrc(n)+1) 
#endif
        if (imAER.ne.5) then
          do j=J_0,J_1
            scca(:) = 0.d0
            do i = I_0,I_1
              if (cosz1(i,j) > 0.) scca(i) = cosz1(i,j) * 4.
            enddo
            trsource(:,j,1,n) = NH3_src_con(:,j) +
     &      NH3_src_cyc(:,j)*SCCA(:)
          enddo
        else
          do ns=1,ntsurfsrc(n)
            trsource(:,J_0:J_1,ns,n)=sfc_src(:,J_0:J_1,n,ns)
     &        *axyp(:,J_0:J_1)
          enddo
        endif
#endif /* TRACERS_NITRATE || TRACERS_AMP */
#ifdef TRACERS_AEROSOLS_Koch
c! TRACERS_AMP
      case ('BCII')
         do j=J_0,J_1
            trsource(:,j,1,n) = BCI_src(:,j)
            trsource(:,j,2,n) = BC_ship(:,j,jmon)
         end do
      case ('OCII')
         do j=J_0,J_1
            trsource(:,j,1,n) = OCI_src(:,j,1)
            trsource(:,j,2,n) = POM_ship(:,j,jmon)
#ifndef TRACERS_AEROSOLS_SOA
            trsource(:,j,3,n) = OCT_src(:,j,jmon)
#endif  /* TRACERS_AEROSOLS_SOA */
         end do
#endif  /* TRACERS_AEROSOLS_Koch */

#ifdef TRACERS_OM_SP
       case ('OCI3')
#ifndef TRACERS_AEROSOLS_SOA
         do j=J_0,J_1
            trsource(:,j,1,n) = OCT_src(:,j,jmon)
         end do
#endif  /* TRACERS_AEROSOLS_SOA */
#endif
      end select

! please keep at end of tracer loop :
      if(alter_sources)then               ! if altering requested
        do ns=1,ntsurfsrc(n)              ! loop over source
          do nsect=1,num_tr_sectors(n,ns) ! and sectors for that source
            do j=J_0,J_1                  ! and latitudes
              do i=I_0,imaxj(j)           ! and longitudes
                do kreg=1,num_regions     ! and defined regions
          if(lat2d_dg(i,j) >= reg_S(kreg) .and. lat2d_dg(i,j)! check if
     &    <= reg_N(kreg) .and. lon2d_dg(i,j) >= reg_W(kreg)  ! in region
     &    .and. lon2d_dg(i,j) < reg_E(kreg) ) then
            if(ef_fact(tr_sect_index(n,ns,nsect),kreg) > -1.e20)
     &      trsource(i,j,ns,n)=trsource(i,j,ns,n)*
     &      ef_FACT(tr_sect_index(n,ns,nsect),kreg)
          endif
                enddo
              enddo
            enddo
          enddo
        enddo
      endif

! optionally set sources to zero over (>90%) ice:
      if(no_emis_over_ice > 0)then
        do j=J_0,J_1
          do i=I_0,imaxj(j)
            fice=flice(i,j)+rsi(i,j)*(focean(i,j)+flake(i,j))
            if(fice > 0.9d0) trsource(i,j,:,:)=0.d0
          enddo
        enddo
      endif

      end do ! n - main tracer loop


      END SUBROUTINE set_tracer_2Dsource

      subroutine get_latlon_mask(lon_w,lon_e,lat_s,lat_n,latlon_mask)
!@sum Set mask array to 1 for all cells overlapping a lat-lon rectangle
!@auth Kelley
      use domain_decomp_atm, only : get,grid
#ifdef CUBED_SPHERE
      use geom, only : lon2d_dg,lat2d_dg
#else
      use geom, only : lon_to_i,lat_to_j
#endif
      implicit none
      real*8 :: lon_w,lon_e,lat_s,lat_n
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     latlon_mask
      integer :: i,j, i_0,i_1,j_0,j_1
      integer :: ie,iw,js,jn

      CALL GET(grid, I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)

      latlon_mask(:,:) = 0d0

#ifdef CUBED_SPHERE
      do j=j_0,j_1
      do i=i_0,i_1
        if(lon2d_dg(i,j) >= lon_w .and.
     &     lon2d_dg(i,j) <= lon_e .and.
     &     lat2d_dg(i,j) >= lat_s .and.
     &     lat2d_dg(i,j) <= lat_n) then
          latlon_mask(i,j) = 1d0
        endif
      enddo
      enddo
c The above approach to defining masks does not work for latlon
c rectangles whose lon and/or lat lengths are less than the
c gridsize.  Solution: if a rectangle is narrow in either
c dimension, trace a discretized path from its SW to its NE
c corner and call lonlat_to_ij for each latlon point along the
c path, setting the mask to 1 for the returned i,j that are in
c the local domain.  Or, trace all 4 edges of the rectangle.
c      if(narrow rectangle) then
c        lon = lon_w
c        lat = lat_s
c        do point=1,npoints_traverse
c          lon = lon + (lon_e-lon_w)/npoints_traverse
c          lat = lat + (lat_n-lat_s)/npoints_traverse
c          call lonlat_to_ij((/lon,lat/),ij)
c          if(i,j in local domain) then
c            latlon_mask(i,j) = 1d0
c          endif
c        enddo
c      endif
#else
c latlon grid
      ie = lon_to_i(lon_e)
      iw = lon_to_i(lon_w)
      jn = lat_to_j(lat_n)
      js = lat_to_j(lat_s)
      do j=max(js,j_0),min(jn,j_1)
        latlon_mask(iw:ie,j) = 1d0
      enddo
#endif
      return
      end subroutine get_latlon_mask

      SUBROUTINE tracer_3Dsource
!@sum tracer_3Dsource calculates interactive sources for tracers
!@+   Please note that if the generic routine 'apply_tracer_3Dsource'
!@+   is used, all diagnostics and moments are updated automatically.
!@auth Jean Lerner/Greg Faluvegi
!@calls DIAGTCA, masterchem, apply_tracer_3Dsource
      USE DOMAIN_DECOMP_ATM, only : GRID, GET, write_parallel,AM_I_ROOT
      USE TRACER_COM
      USE CONSTANT, only : mair, avog,bygrav
      USE FLUXES, only: tr3Dsource
      USE MODEL_COM, only: itime,jmon, dtsrc,jday,jyear,itimeI,
     *  jhour,zatmo
      USE DYNAMICS, only: am,byam ! Air mass of each box (kg/m^2)
      use dynamics, only: phi
      USE apply3d, only : apply_tracer_3Dsource
      USE GEOM, only : byaxyp,axyp
      USE RAD_COM, only: o3_yr
      USE PARAM, only : get_param, is_set_param
      use trdiag_com, only : trcsurf,trcSurfByVol,taijls=>taijls_loc,
     & ijlt_prodSO4gs
#ifdef FAA_emission
     *,ijlt_FAANOX
#endif
CCC#if (defined TRACERS_COSMO) || (defined SHINDELL_STRAT_EXTRA)
#if (defined TRACERS_COSMO)
      USE COSMO_SOURCES, only: be7_src_3d, be10_src_3d, be7_src_param
#endif
#ifdef TRACERS_AEROSOLS_SOA
      USE TRACERS_SOA, only: n_soa_i,n_soa_e
#endif  /* TRACERS_AEROSOLS_SOA */
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP)
      USE AEROSOL_SOURCES, only: so2_src_3d,BCI_src_3d,BCB_src,
     *     OCB_src,SO2_biosrc_3D,lmAER,OCBt_src,BCBt_src,OCI_src
     *     ,so2t_src,BBinc,om2oc

c Laki emissions
c     USE DYNAMICS, only: LTROPO
c     USE CONSTANT, only: sday
c     USE MODEL_COM, only: jday,jyear
c     USE LAKI_SOURCE, only: LAKI_MON,LAKI_DAY,LAKI_AMT_T,LAKI_AMT_S
#endif
      USE PBLCOM, only: dclev
#ifdef TRACERS_AMP
      USE AERO_SETUP, only : RECIP_PART_MASS
      USE TRDIAG_COM, only : itcon_AMP, itcon_AMPe,itcon_AMPm
      USE TRDIAG_COM, only : taijs=>taijs_loc,ijts_AMPe
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: fix_CH4_chemistry,sOx_acc,sNOx_acc,
     & sCO_acc,l1Ox_acc,l1NO2_acc,mNO2
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
      USE TRDIAG_COM, only: sPM2p5_acc,sPM10_acc,l1PM2p5_acc,l1PM10_acc,
     &                      csPM2p5_acc,csPM10_acc
#endif

#ifdef TRACERS_SPECIAL_Shindell
      use RAD_COM, only: rad_to_chem
      use TRCHEM_Shindell_COM, only: fact_cfc, 
     &     use_rad_n2o, use_rad_ch4, use_rad_cfc
#endif
#ifdef FAA_emission
       USE FAA_COM
#endif

      implicit none
      INTEGER n,ns,najl,i,j,l,blay,xyear,xday   ; real*8 now
      INTEGER J_0, J_1, I_0, I_1
      integer :: src_index,bb_i,bb_e
      integer :: initial_ghg_setup
      real*8 :: src_fact,bb_fact
!@var blsrc (m2/s) tr3Dsource (kg/s) in boundary layer,
!@+                per unit of air mass (kg/m2)
      real*8 :: blsrc
      real*8 :: byavog
#ifdef CUBED_SPHERE
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
     &     :: dummy3d
#endif
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP


C**** All sources are saved as kg/s
      do n=1,ntm
      if (itime.lt.itime_tr0(n)) cycle

      select case (trname(n))

      case default
#ifdef TRACERS_SPECIAL_Lerner
C****
      case ('CH4')
      tr3Dsource(:,J_0:J_1,:,:,n) = 0.
      call Trop_chem_CH4(1,n)
      call apply_tracer_3Dsource(1,n)
      call Strat_chem_Prather(2,n)
      call apply_tracer_3Dsource(2,n,.false.)
C****
      case ('O3')
      tr3Dsource(:,J_0:J_1,:,:,n) = 0.
      call Trop_chem_O3(2,3,n)
        call apply_tracer_3Dsource(2,n,.false.)
        call apply_tracer_3Dsource(3,n,.false.)
      call Strat_chem_O3(1,n)
        call apply_tracer_3Dsource(1,n,.false.)
C****
      case ('N2O')
      tr3Dsource(:,J_0:J_1,:,:,n) = 0.
      call Strat_chem_Prather(1,n)
      call apply_tracer_3Dsource(1,n,.FALSE.)
C****
      case ('CFC11')
      tr3Dsource(:,J_0:J_1,:,:,n) = 0.
      call Strat_chem_Prather(1,n)
      call apply_tracer_3Dsource(1,n,.FALSE.)
C****
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      case ('BCIA')
C**** aircraft source for fresh industrial BC
#ifndef FAA_emission
        if (imAER.ne.3.and.imAER.ne.5) tr3Dsource(:,J_0:J_1,:,2,n)
     *    = BCI_src_3d(:,J_0:J_1,:)
        if (imAER.eq.0.or.imAER.eq.2.or.imAER.eq.3.or.imAER.eq.5)
     *    call apply_tracer_3Dsource(2,n) ! aircraft
#endif
#endif /* TRACERS_AEROSOLS_Koch || TRACERS_AMP */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_SPECIAL_Shindell)
      case ('Alkenes', 'CO', 'NOx', 'Paraffin','CH4',
     &      'NH3', 'SO2', 'SO4', 'BCB', 'OCB', 
     &      'M_ACC_SU', 'M_AKK_SU',
     &      'M_BC1_BC', 'M_OCC_OC', 'M_BOC_BC', 'M_BOC_OC')
          src_fact=1.d0 ! factor to multiply emissions with
          bb_fact=1.d0 ! factor to multiply biomass_burning emissions with
          src_index=n   ! index to be used for emissions
          select case (trname(n))
          case ('SO2')
            src_fact=0.975d0 ! the rest goes to sulfate (SO4 or M_ACC_SU)
          case ('SO4')
            src_fact=0.0375d0 ! (1.-SO2 fraction)*tr_mm(n_SO4)/tr_mm(n_SO4)
            src_index=n_SO2
          case ('M_ACC_SU')
            src_fact=0.0375d0
#ifndef TRACERS_AMP_M4
     &              *0.99d0 ! the rest goes to M_AKK_SU
#endif
            src_index=n_SO2
#ifndef TRACERS_AMP_M4
          case ('M_AKK_SU')
            src_fact=0.0375d0
     &              *0.01d0
            src_index=n_SO2
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
          case ('BCB', 'M_BC1_BC', 'M_BOC_BC')
            if(.not.do_fire(n))bb_fact=BBinc
          case ('OCB', 'M_OCC_OC', 'M_BOC_OC')
            if(.not.do_fire(n))then
              src_fact=om2oc(n)
              bb_fact=BBinc
            end if
#endif
          end select

#ifndef FAA_emission
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
C**** 3D aircraft source - only SO2, no sulfate
        select case (trname(n))
        case ('SO2')
          if (imAER.ne.3.and.imAER.ne.5)
     &    tr3Dsource(:,J_0:J_1,:,nAircraft,n) =
     &      so2_src_3d(:,J_0:J_1,:,2)
          if (imAER.eq.0.or.imAER.eq.2.or.imAER.eq.3.or.imAER.eq.5)
     &     call apply_tracer_3Dsource(nAircraft,n)
        end select
#endif
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
C**** 3D volcanic source
        select case (trname(n))
        case ('SO2', 'SO4', 'M_ACC_SU', 'M_AKK_SU')
          tr3Dsource(:,J_0:J_1,:,nVolcanic,n) =
     &      so2_src_3d(:,J_0:J_1,:,1)*src_fact
          call apply_tracer_3Dsource(nVolcanic,n)
        end select
#endif

C**** 3D biomass source
        tr3Dsource(:,J_0:J_1,:,nBiomass,n) = 0.
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
        if ((imAER.ne.1) .or. (n<=ntm_chem)) then
#endif
          if(do_fire(src_index) .or. nBBsources(src_index) > 0) then
            bb_i=ntsurfsrc(src_index)+1 ! index of first BB source
            if(do_fire(src_index))then
              bb_e=bb_i ! index last BB source
            else
              bb_e=ntsurfsrc(src_index)+nBBsources(src_index) ! index last BB source
            end if
            do j=J_0,J_1; do i=I_0,I_1
              blay=int(dclev(i,j)+0.5d0)
              blsrc = axyp(i,j)*src_fact*bb_fact*
     &         sum(sfc_src(i,j,src_index,bb_i:bb_e))/sum(am(1:blay,i,j))
              do l=1,blay
                tr3Dsource(i,j,l,nBiomass,n) = blsrc*am(l,i,j)
              end do
            end do; end do
          end if 
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
        else
          select case(trname(n))
!kt - not ready - SO2t_src has one extra dimension
!          case ('SO4')
!            tr3Dsource(:,J_0:J_1,1:lmAER,nBiomass,n) =
!     &      SO2t_src(:,J_0:J_1,:,jmon)*0.0375d0
!          case ('M_ACC_SU')
!            tr3Dsource(:,J_0:J_1,1:lmAER,nBiomass,n) =
!#ifndef TRACERS_AMP_M4
!     &      0.99d0*
!#endif
!     &      SO2t_src(:,J_0:J_1,:,jmon)*0.0375d0
!          case ('M_AKK_SU')
!            tr3Dsource(:,J_0:J_1,1:lmAER,nBiomass,n) =
!     &      0.01d0*
!     &      SO2t_src(:,J_0:J_1,:,jmon)*0.0375d0
!kt
          case ('BCB', 'M_BC1_BC', 'M_BOC_BC')
            tr3Dsource(:,J_0:J_1,1:lmAER,nBiomass,n) =
     &        BCB_src(:,J_0:J_1,:,jmon)
#ifndef FAA_emission
            select case (trname(n))
            case ('M_BC1_BC')
              tr3Dsource(:,J_0:J_1,:,nAircraft,n) =
     &          BCI_src_3d(:,J_0:J_1,:)
            end select
#endif
          case ('OCB', 'M_OCC_OC', 'M_BOC_OC')
            tr3Dsource(:,J_0:J_1,1:lmAER,nBiomass,n) =
     &        OCB_src(:,J_0:J_1,:,jmon)
          end select
        endif ! imAER
#ifndef FAA_emission
        if(n>ntm_chem)call apply_tracer_3Dsource(nAircraft,n)
#endif
#endif
#ifndef GFED_3D_BIOMASS
        call apply_tracer_3Dsource(nBiomass,n)
#else
        if(n>ntm_chem)call apply_tracer_3Dsource(nBiomass,n)
#endif /* note that is a NOT defined */
#endif /* TRACERS_AEROSOLS_Koch || TRACERS_AMP || TRACERS_SPECIAL_Shindell */

#ifdef TRACERS_OM_SP
      case ('OCI1')
c**** biomass + industrial for OC species
      tr3Dsource(:,J_0:J_1,:,1,n) = 0.
      do j=J_0,J_1; do i=I_0,I_1
      blay=int(dclev(i,j)+0.5)
      do l=1,blay
      tr3Dsource(i,j,l,1,n) = OCBt_src(i,j)/real(blay)*OCI_src(i,j,5)
      end do
      tr3Dsource(i,j,1,1,n) = tr3Dsource(i,j,1,1,n)+OCI_src(i,j,1)
      end do; end do
      call apply_tracer_3Dsource(1,n) ! biomass+industrial
      case ('OCI2')
c**** biomass + industrial for OC species
      tr3Dsource(:,J_0:J_1,:,1,n) = 0.
      do j=J_0,J_1; do i=I_0,I_1
      blay=int(dclev(i,j)+0.5)
      do l=1,blay
      tr3Dsource(i,j,l,1,n) = OCBt_src(i,j)/real(blay)*OCI_src(i,j,6)
      end do
      tr3Dsource(i,j,1,1,n) = tr3Dsource(i,j,1,1,n)+OCI_src(i,j,2)
      end do; end do
      call apply_tracer_3Dsource(1,n) ! biomass+industrial
      case ('OCI3')
c**** biomass + industrial for OC species
      tr3Dsource(:,J_0:J_1,:,1,n) = 0.
      do j=J_0,J_1; do i=I_0,I_1
      blay=int(dclev(i,j)+0.5)
      do l=1,blay
      tr3Dsource(i,j,l,1,n) = OCBt_src(i,j)/real(blay)*OCI_src(i,j,7)
      end do
      tr3Dsource(i,j,1,1,n) = tr3Dsource(i,j,1,1,n)+OCI_src(i,j,3)
      end do; end do
      call apply_tracer_3Dsource(1,n) ! biomass+industrial
      case ('OCA4')
c**** biomass + industrial for OC species
      tr3Dsource(:,J_0:J_1,:,1,n) = 0.
      do j=J_0,J_1; do i=I_0,I_1
      blay=int(dclev(i,j)+0.5)
      do l=1,blay
      tr3Dsource(i,j,l,1,n) = OCBt_src(i,j)/real(blay)*OCI_src(i,j,8)
      end do
      tr3Dsource(i,j,1,1,n) = tr3Dsource(i,j,1,1,n)+OCI_src(i,j,4)
      end do; end do
      call apply_tracer_3Dsource(1,n) ! biomass+industrial
#endif
CCC#if (defined TRACERS_COSMO) || (defined SHINDELL_STRAT_EXTRA)
#if (defined TRACERS_COSMO)
C****
      case ('Be7')
c cosmogenic src
        do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
          tr3Dsource(i,j,l,1,n) = am(l,i,j)*be7_src_3d(i,j,l)
        end do; end do; end do

        call apply_tracer_3Dsource(1,n)
C****
      case ('Be10')
c cosmogenic src
        do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
           tr3Dsource(i,j,l,1,n) = am(l,i,j)*be10_src_3d(i,j,l)
        end do; end do; end do

        call apply_tracer_3Dsource(1,n)
C****
#endif
      case('Pb210')
        call apply_tracer_3Dsource(1,n) !radioactive decay of Rn222

      end select

      end do

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AMP)
c Calculation of gas phase reaction rates for sulfur chemistry
      CALL GET_SULF_GAS_RATES
#endif

#ifdef TRACERS_AEROSOLS_Koch
       call aerosol_gas_chem
       call apply_tracer_3Dsource(1,n_DMS)  ! DMS chem sink
       call apply_tracer_3Dsource(1,n_MSA)  ! MSA chem source
       call apply_tracer_3Dsource(nChemistry,n_SO2)  ! SO2 chem source
       call apply_tracer_3Dsource(nChemloss,n_SO2)  ! SO2 chem sink
#ifdef ACCMIP_LIKE_DIAGS
       do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
         taijls(i,j,l,ijlt_prodSO4gs)=taijls(i,j,l,ijlt_prodSO4gs)+
     &   tr3Dsource(i,j,l,nChemistry,n_SO4)*dtsrc*byaxyp(i,j)
       end do; end do; end do
#endif
       call apply_tracer_3Dsource(nChemistry,n_SO4)  ! SO4 chem source
       call apply_tracer_3Dsource(1,n_H2O2_s) ! H2O2 chem source
       call apply_tracer_3Dsource(2,n_H2O2_s) ! H2O2 chem sink
       call apply_tracer_3Dsource(1,n_BCII) ! BCII aging sink
       call apply_tracer_3Dsource(1,n_BCIA) ! BCIA aging source
       call apply_tracer_3Dsource(1,n_OCII) ! OCII aging sink
       call apply_tracer_3Dsource(1,n_OCIA) ! OCIA aging source

#ifdef TRACERS_HETCHEM
       call apply_tracer_3Dsource(1,n_SO4_d1) ! SO4 chem prod on dust
       call apply_tracer_3Dsource(1,n_SO4_d2) ! SO4 chem prod on dust
       call apply_tracer_3Dsource(1,n_SO4_d3) ! SO4 chem prod on dust
#endif
#endif
#ifdef TRACERS_OM_SP
       call aerosol_gas_chem
       call apply_tracer_3Dsource(2,n_OCI1) ! OCI1 aging sink
       call apply_tracer_3Dsource(1,n_OCA1) ! OCA1 aging source
       call apply_tracer_3Dsource(2,n_OCI2) ! OCI2 aging sink
       call apply_tracer_3Dsource(1,n_OCA2) ! OCA2 aging source
       call apply_tracer_3Dsource(2,n_OCI3) ! OCI3 aging sink
       call apply_tracer_3Dsource(1,n_OCA3) ! OCA3 aging source
#endif
C nu/mk FAA hourly resolution aviation emmissions
C Must have Shindell chemistry and Koch aerosols switched on
#ifdef FAA_emission
C water (g -> kg/s)
      tr3Dsource(I_0:I_1,J_0:J_1,:,1,n_Water)  = 0.
      tr3Dsource(I_0:I_1,J_0:J_1,:,1,n_Water)  =
     * emiss_5d(I_0:I_1,J_0:J_1,:,jhour+1,1)*1231.0*
     * (1.0d-3/3600.0)
      call apply_tracer_3Dsource(1,n_Water)
C Units NOx (g) -> kgN /s 
C NOx
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_NOx)  = 0.
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_NOx)  =
     * emiss_5d(I_0:I_1,J_0:J_1,:,jhour+1,4)*(14.0/46.0)
     * *(1.0d-3/3600.0)
      call apply_tracer_3Dsource(nAircraft,n_NOx)
!C ADD IJLT DIAG FOR FAA NOX -/divide by grid box area (to give kg/m2/s)
!       do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
!       taijls(i,j,l,ijlt_FAANOX)=taijls(i,j,l,ijlt_FAANOX)+
!     * tr3Dsource(i,j,l,nAircraft,n_NOx)*dtsrc*byaxyp(i,j)
!        end do; end do; end do
C SO2 from fuel burn, fuel sulfur content and % emission as S(IV)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_SO2)  = 0.
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_SO2)  =
     * emiss_5d(I_0:I_1,J_0:J_1,:,jhour+1,1)
     * *((600.d0/1000.d0)*((100.0d0-2.0d0)/100)*(64.0/32.0))
     * *(1.0d-3/3600.0)
      call apply_tracer_3Dsource(nAircraft,n_SO2)
C BCIA
      tr3Dsource(I_0:I_1,J_0:J_1,:,2,n_BCIA)  = 0.
      do i=I_0,I_1; do j=J_0,J_1; do l=1,lm 
      if ((daily_z(i,j,l)-(zatmo(i,j)*bygrav)) <= 914.4) then
      tr3Dsource(i,j,l,2,n_BCIA)  =
     * emiss_5d(i,j,l,jhour+1,5)      
     * *(1.0d-3/3600.0)
      else
      tr3Dsource(i,j,l,2,n_BCIA)  =
     * emiss_5d(i,j,l,jhour+1,1)*0.03      
     * *(1.0d-3/3600.0)
      endif
      end do; end do; end do
      call apply_tracer_3Dsource(2,n_BCIA)
C CO
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_CO)  = 0.
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_CO)  =
     * emiss_5d(I_0:I_1,J_0:J_1,:,jhour+1,2)
     * *(1.0d-3/3600.0)
      call apply_tracer_3Dsource(nAircraft,n_CO)
C OCIA 
      tr3Dsource(I_0:I_1,J_0:J_1,:,2,n_OCIA)  = 0.
      do i=I_0,I_1; do j=J_0,J_1; do l=1,lm 
      if ((daily_z(i,j,l)-(zatmo(i,j)*bygrav)) <= 914.4) then
      tr3Dsource(i,j,l,2,n_OCIA)  =
     * emiss_5d(i,j,l,jhour+1,6)      
     * *1.0d-3/3600.0
      else
      tr3Dsource(i,j,l,2,n_OCIA)  =
     * emiss_5d(i,j,l,jhour+1,1)*0.03d0      
     * *(1.0d-3/3600.0d0)
      endif
      end do; end do; end do
      call apply_tracer_3Dsource(2,n_OCIA)
C Sulfate
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_SO4)  = 0.
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_SO4)  =
     * emiss_5d(I_0:I_1,J_0:J_1,:,jhour+1,1)
     * *((600.d0/1000.d0)*((2.0d0)/100)*(96.0d0/32.0d0))
     * *(1.0d-3/3600.0d0)
      call apply_tracer_3Dsource(nAircraft,n_SO4)
C Alkenes (as TOG in CH4-equivalents)
C CAREFUL - GF's change so tr_mm=1g/mol for Alkenes and Paraffin
C (1/16=0.0625)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_Alkenes)  = 0.
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_Alkenes)  =
     * emiss_5d(I_0:I_1,J_0:J_1,:,jhour+1,3)*1.16d0
     * *(1.0d-3/3600.0d0)*0.381*(0.0625)
      call apply_tracer_3Dsource(nAircraft,n_Alkenes)
C Paraffin (as TOG in CH4-equivalents)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_Paraffin)  = 0.
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_Paraffin)  =
     * emiss_5d(I_0:I_1,J_0:J_1,:,jhour+1,3)*1.16d0
     * *(1.0d-3/3600.0d0)*0.203*(0.0625)
      call apply_tracer_3Dsource(nAircraft,n_Paraffin)
#endif
#ifdef TRACERS_SPECIAL_Shindell
C Apply non-chemistry 3D sources, so they can be "seen" by chemistry:
C (Note: using this method, tracer moments are changed just like they
C are done for chemistry.  It might be better to do it like surface
C sources are done? -- GSF 11/26/02)
c
      CALL TIMER (NOW,MTRACE)
C**** Allow overriding of transient emissions date:
! for now, tying this to o3_yr becasue Gavin
! didn't want a new parameter, also not allowing
! day overriding yet, because of that.
      trans_emis_overr_yr=ABS(o3_yr)
      if(trans_emis_overr_yr > 0)then
        xyear=trans_emis_overr_yr
      else
        xyear=jyear
      endif
!!    if(trans_emis_overr_day > 0)then
!!      xday=trans_emis_overr_day
!!    else
        xday=jday
!!    endif
#ifdef SHINDELL_STRAT_EXTRA
      tr3Dsource(I_0:I_1,J_0:J_1,:,1,n_GLT) = 0.d0
      call overwrite_GLT
      call apply_tracer_3Dsource(1,n_GLT)
#endif
#ifndef FAA_emission
      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,n_NOx)  = 0.
#ifdef CUBED_SPHERE
      call get_aircraft_NOx(xyear,xday,dummy3d,.false.)
#else
      call get_aircraft_NOx(xyear,xday,phi,.true.) ! read from disk
#endif
      call apply_tracer_3Dsource(nAircraft,n_NOx)
#endif
      tr3Dsource(I_0:I_1,J_0:J_1,:,nOther,n_NOx) = 0.d0
      call get_lightning_NOx
      call apply_tracer_3Dsource(nOther,n_NOx)
#ifdef GFED_3D_BIOMASS
      tr3Dsource(I_0:I_1,J_0:J_1,:,nBiomass,n_CH4) = 0.d0
      if(fix_CH4_chemistry/=1.or.ntsurfsrc(n_CH4)/=0)then
        call dist_GFED_biomass_burning(n_CH4)
        call apply_tracer_3Dsource(nBiomass,n_CH4)
      endif
      tr3Dsource(I_0:I_1,J_0:J_1,:,nBiomass,n_NOx) = 0.d0
      call dist_GFED_biomass_burning(n_NOx)
      call apply_tracer_3Dsource(nBiomass,n_NOx)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nBiomass,n_CO) = 0.d0
      call dist_GFED_biomass_burning(n_CO)
      call apply_tracer_3Dsource(nBiomass,n_CO)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nBiomass,n_Alkenes) = 0.d0
      call dist_GFED_biomass_burning(n_Alkenes)
      call apply_tracer_3Dsource(nBiomass,n_Alkenes)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nBiomass,n_Paraffin) = 0.d0
      call dist_GFED_biomass_burning(n_Paraffin)
      call apply_tracer_3Dsource(nBiomass,n_Paraffin)
#endif /*GFED_3D_BIOMASS*/

C**** Make sure that these 3D sources for all chem tracers start at 0.:
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,1:ntm_chem)  = 0.d0
      tr3Dsource(I_0:I_1,J_0:J_1,:,nOverwrite,1:ntm_chem) = 0.d0
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_stratOx)  = 0.d0
      tr3Dsource(I_0:I_1,J_0:J_1,:,nOverwrite,n_stratOx) = 0.d0
#endif
#if (defined TRACERS_HETCHEM) && (defined TRACERS_NITRATE)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_N_d1)  = 0.d0
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_N_d2)  = 0.d0
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_N_d3)  = 0.d0
#endif
#ifdef TRACERS_AEROSOLS_SOA
      tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_soa_i:n_soa_e)  = 0.d0
#endif  /* TRACERS_AEROSOLS_SOA */


      if (is_set_param('initial_ghg_setup')) then
        call get_param('initial_GHG_setup', initial_GHG_setup)
        if (initial_GHG_setup == 1 .and. itime == itimeI) then
          byavog = 1.d+0/avog

          if (use_rad_n2o > 0) call applyRadChem(3, n_N2O, 1.d+0)
          if (use_rad_ch4 > 0) call applyRadChem(4, n_CH4, 1.d+0)
          if (use_rad_cfc > 0) call applyRadChem(5, n_CFC, fact_CFC)

        end if
      end if

C**** Call the model CHEMISTRY and OVERWRITEs:
      call masterchem ! does chemistry and over-writing.
                      ! tr3Dsource defined within, for both processes

C**** Apply chemistry and overwrite changes:
      do n=1,ntm_chem
        call apply_tracer_3Dsource(nChemistry,n)
        call apply_tracer_3Dsource(nOverwrite,n)
      end do
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      call apply_tracer_3Dsource(nChemistry,n_stratOx)
      call apply_tracer_3Dsource(nOverwrite,n_stratOx)
#endif
#if (defined TRACERS_HETCHEM) && (defined TRACERS_NITRATE)
       call apply_tracer_3Dsource(nChemistry,n_N_d1) ! NO3 chem prod on dust
       call apply_tracer_3Dsource(nChemistry,n_N_d2) ! NO3 chem prod on dust
       call apply_tracer_3Dsource(nChemistry,n_N_d3) ! NO3 chem prod on dust
#endif
      CALL TIMER (NOW,MCHEM)
#endif /* TRACERS_SPECIAL_Shindell */
#ifdef TRACERS_NITRATE
#ifdef TRACERS_SPECIAL_Shindell
       tr3Dsource(I_0:I_1,J_0:J_1,:,3,n_HNO3) = 0.d0
#endif
       tr3Dsource(I_0:I_1,J_0:J_1,:,1,n_NO3p) = 0.d0
       tr3Dsource(I_0:I_1,J_0:J_1,:,1,n_NH4)  = 0.d0
       tr3Dsource(I_0:I_1,J_0:J_1,:,nChemistry,n_NH3)  = 0.d0

       call EQSAM_DRV

#ifdef TRACERS_SPECIAL_Shindell
       call apply_tracer_3Dsource(3,n_HNO3) ! NO3 chem prod
#endif
       call apply_tracer_3Dsource(1,n_NO3p) ! NO3 chem prod
       call apply_tracer_3Dsource(1,n_NH4)  ! NO3 chem prod
       call apply_tracer_3Dsource(nChemistry,n_NH3)  ! NH3

#endif /* TRACERS_NITRATE */
#ifdef TRACERS_AMP

       call aerosol_gas_chem

       call apply_tracer_3Dsource(2,n_H2SO4) ! H2SO4 chem prod
       call apply_tracer_3Dsource(1,n_DMS)  ! DMS chem sink
       call apply_tracer_3Dsource(nChemistry,n_SO2)  ! SO2 chem source
       call apply_tracer_3Dsource(nChemloss,n_SO2)  ! SO2 chem sink
       call apply_tracer_3Dsource(1,n_H2O2_s)! H2O2 chem source
       call apply_tracer_3Dsource(2,n_H2O2_s)! H2O2 chem sink
      DO n=1,ntmAMP
        tr3Dsource(:,J_0:J_1,:,1,n)  = 0.d0! Aerosol Mirophysics
      ENDDO
        tr3Dsource(:,J_0:J_1,:,1,n_H2SO4)  = 0.d0! Aerosol Mirophysics
        tr3Dsource(:,J_0:J_1,:,nChemistry,n_NH3)  = 0.d0! Aerosol Mirophysics
#ifdef  TRACERS_SPECIAL_Shindell
        tr3Dsource(:,J_0:J_1,:,3,n_HNO3)  = 0.d0! Aerosol Mirophysics
#endif
        call MATRIX_DRV

      DO n=1,ntmAMP
       call apply_tracer_3Dsource(nChemistry,n) ! Aerosol Mirophysics !kt is the index correct?
      ENDDO

       call apply_tracer_3Dsource(nChemistry,n_NH3)  ! NH3
       call apply_tracer_3Dsource(1,n_H2SO4) ! H2SO4 chem prod
#ifdef  TRACERS_SPECIAL_Shindell
       call apply_tracer_3Dsource(3,n_HNO3) ! H2SO4 chem prod
#endif
#endif /* TRACERS_AMP */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_SOA)
! This section is to accumulate/aggregate certain tracers' SURFACE and
! L=1 values into particulate matter PM2.5 and PM10 for use in the sub-
! daily diags. Saved in ppmm or kg/m3. Also save Ox and NO2 in ppmv:
      do n=1,ntm
        select case (trname(n))
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_AEROSOLS_SOA)
! 100% of these: <-----------------------------------------------------
        case('BCII','BCIA','BCB','OCII','OCIA','OCB','SO4','NO3p',
#ifdef TRACERS_AEROSOLS_SOA
     &       'isopp1a', 'isopp2a', 'apinp1a', 'apinp2a',
#endif  /* TRACERS_AEROSOLS_SOA */
     &       'Clay','seasalt1','N_d1','SO4_d1')
          sPM2p5_acc(:,:)=sPM2p5_acc(:,:)  + 1.d6*trcsurf(:,:,n)
          sPM10_acc(:,:)=sPM10_acc(:,:)    + 1.d6*trcsurf(:,:,n)
          csPM2p5_acc(:,:)=csPM2p5_acc(:,:)  + trcSurfByVol(:,:,n)
          csPM10_acc(:,:)=csPM10_acc(:,:)    + trcSurfByVol(:,:,n)
          L1PM2p5_acc(:,:)=L1PM2p5_acc(:,:)+
     &                 trm(:,:,1,n)*byam(1,:,:)*byaxyp(:,:)*1.d6
          L1PM10_acc(:,:)=L1PM10_acc(:,:)  +
     &                 trm(:,:,1,n)*byam(1,:,:)*byaxyp(:,:)*1.d6
! then, conditional or partial of these: <=============================
        case('Silt1','N_d2','SO4_d2')
          sPM2p5_acc(:,:)=sPM2p5_acc(:,:)  + 0.322d0*1.d6*trcsurf(:,:,n)
          sPM10_acc(:,:)=sPM10_acc(:,:)    + 1.d6*trcsurf(:,:,n)
          csPM2p5_acc(:,:)=csPM2p5_acc(:,:)+ 0.322d0*trcSurfByVol(:,:,n)
          csPM10_acc(:,:)=csPM10_acc(:,:)  + trcSurfByVol(:,:,n)
          L1PM2p5_acc(:,:)=L1PM2p5_acc(:,:)+ 0.322d0*
     &                         trm(:,:,1,n)*byam(1,:,:)*byaxyp(:,:)*1.d6
          L1PM10_acc(:,:)=L1PM10_acc(:,:)  +
     &                         trm(:,:,1,n)*byam(1,:,:)*byaxyp(:,:)*1.d6
        case('Silt2','N_d3','SO4_d3')
          sPM10_acc(:,:)=sPM10_acc(:,:)    + 1.d6*trcsurf(:,:,n)
          csPM10_acc(:,:)=csPM10_acc(:,:)  + trcSurfByVol(:,:,n)
          L1PM10_acc(:,:)=L1PM10_acc(:,:)  +
     &                 trm(:,:,1,n)*byam(1,:,:)*byaxyp(:,:)*1.d6
        case('Silt3')
          sPM10_acc(:,:)=sPM10_acc(:,:)    + 0.322d0*1.d6*trcsurf(:,:,n)
          csPM10_acc(:,:)=csPM10_acc(:,:)  + 0.322d0*trcSurfByVol(:,:,n)
          L1PM10_acc(:,:)=L1PM10_acc(:,:)  + 0.322d0*
     &                         trm(:,:,1,n)*byam(1,:,:)*byaxyp(:,:)*1.d6
        case('seasalt2')
          sPM2p5_acc(:,:)=sPM2p5_acc(:,:)  + 0.500d0*1.d6*trcsurf(:,:,n)
          sPM10_acc(:,:)=sPM10_acc(:,:)    +         1.d6*trcsurf(:,:,n)
          csPM2p5_acc(:,:)=csPM2p5_acc(:,:)+ 0.500d0*trcSurfByVol(:,:,n)
          csPM10_acc(:,:)=csPM10_acc(:,:)  +         trcSurfByVol(:,:,n)
          L1PM2p5_acc(:,:)=L1PM2p5_acc(:,:)+ 0.500d0*
     &                        trm(:,:,1,n)*byam(1,:,:)*byaxyp(:,:)*1.d6
          L1PM10_acc(:,:)=L1PM10_acc(:,:)  +
     &                        trm(:,:,1,n)*byam(1,:,:)*byaxyp(:,:)*1.d6
#endif
#ifdef TRACERS_SPECIAL_Shindell
        case('Ox')
          sOx_acc(:,:)=sOx_acc(:,:)+1.d9*trcsurf(:,:,n)*mass2vol(n)
          l1Ox_acc(:,:)=l1Ox_acc(:,:)+trm(:,:,1,n)*byam(1,:,:)*
     &    byaxyp(:,:)*1.d9*mass2vol(n)
        case('NOx')
          sNOx_acc(:,:)=sNOx_acc(:,:)+1.d9*trcsurf(:,:,n)*mass2vol(n)
          l1NO2_acc(:,:)=l1NO2_acc(:,:)+mNO2(:,:,1)*1.d9  ! note: NO2 not NOx
        case('CO') 
          sCO_acc(:,:)=sCO_acc(:,:)+1.d9*trcsurf(:,:,n)*mass2vol(n)
#endif
        end select
      enddo
#endif

#ifdef CACHED_SUBDD
      ! Accumulate the tracer-related subdaily diagnostics
      ! (seems like a reasonable place to put this, as tracer_3Dsource
      ! is called each time step, and here chemistry has been done, etc.,
      ! but we could move it):
      call accumCachedTracerSUBDDs
#endif

      return

#ifdef TRACERS_SPECIAL_Shindell
      contains

      subroutine applyRadChem(index, n, factor)
      integer, intent(in) :: index
      integer, intent(in) :: n
      real*8, intent(in) :: factor
      
      integer :: L
      do L = 1, LM
        tr3Dsource(I_0:I_1,J_0:J_1,L,nOverwrite,n) = 
     &       (rad_to_chem(index,L,I_0:I_1,J_0:J_1)*2.69e20*byavog*
     &       axyp(I_0:I_1,J_0:J_1)*tr_mm(n) * factor - 
     &       trm(I_0:I_1,J_0:J_1,L,n)) / dtsrc
      end do
      call apply_tracer_3Dsource(nOverwrite,n)
      tr3Dsource(I_0:I_1,J_0:J_1,:,nOverwrite,n) = 0.d0

      end subroutine applyRadChem
#endif

      END SUBROUTINE tracer_3Dsource
#endif /* TRACERS_ON */

#ifdef TRACERS_WATER
C---SUBROUTINES FOR TRACER WET DEPOSITION-------------------------------

      SUBROUTINE GET_COND_FACTOR(L,N,WMXTR,TEMP,TEMP0,LHX,FCLOUD,FQ0,fq,
     *  TR_CONV,TRWML,TM,THLAW,TR_LEF,pl,ntix,CLDSAVT)
!@sum  GET_COND_FACTOR calculation of condensate fraction for tracers
!@+    within or below convective or large-scale clouds. Gas
!@+    condensation uses Henry's Law if not freezing.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 CLOUDCHCC and CLOUDCHEM subroutines)
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only: BYGASC, MAIR,teeny,LHE,tf,by3
      USE TRACER_COM, only: tr_RKD,tr_DHD,nWATER,nGAS,nPART,tr_wd_TYPE
     *     ,trname,ntm,lm,t_qlimit,fq_aer,trpdens
#ifdef TRACERS_SPECIAL_O18
     &     ,supsatfac
#endif
#ifdef TRACERS_HETCHEM
     *     ,trm ,n_SO4_d1, n_SO4_d2, n_SO4_d3,n_SO4
     *     ,n_N_d1,n_N_d2,n_N_d3,n_NO3p
      USE MODEL_COM, only  : dtsrc
#endif
      IMPLICIT NONE
C**** Local parameters and variables and arguments:
!@param BY298K unknown meaning for now (assumed= 1./298K)
!@var Ppas pressure at current altitude (in Pascal=kg/s2/m)
!@var TFAC exponential coeffiecient of tracer condensation temperature
!@+   dependence (mole/joule)
!@var FCLOUD fraction of cloud available for tracer condensation
!@var SSFAC dummy variable (assumed units= kg water?)
!@var FQ            fraction of tracer that goes into condensate
!@var FQ0 default fraction of water tracer that goes into condensate
!@var L index for altitude loop
!@var N index for tracer number loop
!@var WMXTR mixing ratio of water available for tracer condensation?
!@var SUPSAT super-saturation ratio for cloud droplets
!@var LHX latent heat flag for whether condensation is to ice or water
!@var RKD dummy variable (= tr_RKD*EXP[ ])
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac, RKD,CLDINC
#ifdef TRACERS_SPECIAL_O18
      real*8 tdegc,alph,fracvs,fracvl,kin_cond_ice,fqi,gint
      integer i
!@param nstep no. of steps for integration of Rayleigh condensation
      integer, parameter :: nstep=6   !8
!@param wgt weightings for simpson's rule integration
      real*8, parameter, dimension(nstep+1) ::
     *     wgt = (/ by3, 4*by3, 2*by3, 4*by3, 2*by3, 4*by3, by3 /)
c     *     wgt = (/ by3, 4*by3, 2*by3, 4*by3, 2*by3, 4*by3, 2*by3, 4*by3
c     *     , by3 /)
#endif
      REAL*8,  INTENT(IN) :: fq0, FCLOUD, WMXTR, TEMP, TEMP0,LHX, TR_LEF
     *     , pl,CLDSAVT
      REAL*8,  INTENT(IN), DIMENSION(ntm,lm) :: trwml
      REAL*8,  INTENT(IN), DIMENSION(lm,ntm) :: TM
      REAL*8,  INTENT(OUT):: fq,thlaw
      INTEGER, INTENT(IN) :: L, N, ntix(ntm)
      LOGICAL TR_CONV
      REAL*8 :: SUPSAT
c
C**** CALCULATE the fraction of tracer mass that becomes condensate:
c
      thlaw=0.
      SELECT CASE(tr_wd_TYPE(NTIX(N)))
        CASE(nGAS)                            ! gas tracer
          fq = 0.D0                           ! frozen and default case
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_SPECIAL_Shindell) ||\
    (defined TRACERS_AMP)
          IF(LHX.eq.LHE) THEN                 ! if not frozen then:
            Ppas = PL*1.D2                 ! pressure to pascals
            tfac = (1.D0/TEMP - BY298K)*BYGASC
            IF(tr_DHD(NTIX(N)).ne.0.D0) THEN
              RKD=tr_RKD(NTIX(N))*DEXP(-tr_DHD(NTIX(N))*tfac)
            ELSE
              RKD=tr_RKD(NTIX(N))
            END IF
c           clwc=WMXTR*MAIR*1.D-6*Ppas*BYGASC/(TEMP*FCLOUD)
c           ssfac=RKD*GASC*TEMP*clwc   ! Henry's Law
            ssfac=RKD*WMXTR*MAIR*1.D-6*Ppas/(CLDSAVT+teeny)
            if (.not.tr_conv) then  !stratiform
              thlaw=(ssfac*tr_lef*tm(l,NTIX(N))
     *        -TRWML(NTIX(N),L))/(1.D0+ssfac)
              if (thlaw.lt.0.) thlaw=0.d0
              if (thlaw.gt.tm(l,NTIX(N))) thlaw=tm(l,NTIX(N))
            else  !if convection
              fq=ssfac / (1.D0 + ssfac)
              if (fq.ge.1.) fq=1.d0
              thlaw=0.
            endif
            if (FCLOUD.LT.1.D-16) fq=0.d0
            if (FCLOUD.LT.1.D-16) thlaw=0.d0
            if (fq0.eq.0.) fq=0.d0
#ifdef TRACERS_SPECIAL_Shindell
            if(t_qlimit(NTIX(N)).and.fq.gt.1.)fq=1.!no negative tracers
#endif
          END IF
#endif
        CASE(nWATER)                          ! water tracer
#ifdef TRACERS_SPECIAL_O18
          if (fq0.gt.0. .and. fq0.lt.1.) then
C**** If process occurs at constant temperature, calculate condensate
C**** in equilibrium with source vapour. Otherwise, use mid-point
C**** temperature and estimate instantaneous fractionation. This gives
C**** a very good estimate to complete integral
C****
            if (abs(temp-temp0).gt.1d-14) then  ! use instantaneous frac
              tdegc=0.5*(temp0 + temp) -tf
C**** Calculate alpha (fractionation coefficient)
                if (LHX.eq.LHE) then ! cond to water
                  alph=1./fracvl(tdegc,ntix(n))
                else            ! cond to ice
                  alph=1./fracvs(tdegc,ntix(n))
C**** kinetic fractionation can occur as a function of supersaturation
C**** this is a parameterisation from Georg Hoffmann
                  supsat=1d0-supsatfac*tdegc
                  if (supsat .gt. 1.) alph=kin_cond_ice(alph,supsat
     *                 ,ntix(n))
                end if
                fq = 1.- (1.-fq0)**alph
            else
C**** assume condensate in equilibrium with vapour at temp
              tdegc=temp -tf
              if (LHX.eq.LHE) then ! cond to water
                alph=1./fracvl(tdegc,ntix(n))
              else              ! cond to ice
                alph=1./fracvs(tdegc,ntix(n))
C**** kinetic fractionation can occur as a function of supersaturation
C**** this is a parameterisation from Georg Hoffmann
                supsat=1d0-supsatfac*tdegc
                if (supsat .gt. 1.) alph=kin_cond_ice(alph,supsat
     *               ,ntix(n))
              end if
              fq = alph * fq0/(1.+(alph-1.)*fq0)
            end if
          else
            fq = fq0
          end if
#else
          fq = fq0
#endif
        CASE(nPART)                           ! particulate tracer
          fq = 0.D0                           ! defaults to zero.
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_COSMO) ||\
    (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP) || (defined TRACERS_RADON)
c only dissolve if the cloud has grown
#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_DUST) &&\
    (defined TRACERS_HETCHEM)
      select case(trname(ntix(n)))
      case('Clay')
         if ( ( TM(l,ntix(n_SO4_d1)) /trpdens(n_SO4)) >
     *      (( TM(l,ntix(n))  /trpdens(n)) * 0.03 ) ) then
            fq_aer(NTIX(N))  = 1.
         else
            fq_aer(NTIX(N))  = 0.
         endif
#ifdef TRACERS_NITRATE
         if ( ( TM(l,ntix(n_N_d1)) /trpdens(n_NO3p)) >
     *      (( TM(l,ntix(n))  /trpdens(n)) * 0.03 ) ) then
            fq_aer(NTIX(N))  = 1.
         endif
#endif

      case('Silt1')
        if ( ( TM(l,ntix(n_SO4_d2)) /trpdens(n_SO4)) >
     *      (( TM(l,ntix(n))  /trpdens(n)) * 0.03 ) ) then
            fq_aer(NTIX(N))  = 1.
         else
            fq_aer(NTIX(N))  = 0.
        endif
#ifdef TRACERS_NITRATE
         if ( ( TM(l,ntix(n_N_d2)) /trpdens(n_NO3p)) >
     *      (( TM(l,ntix(n))  /trpdens(n)) * 0.03 ) ) then
            fq_aer(NTIX(N))  = 1.
         endif
#endif

      case('Silt2')
        if ( ( TM(l,ntix(n_SO4_d3)) /trpdens(n_SO4)) >
     *      (( TM(l,ntix(n))  /trpdens(n)) * 0.03 ) ) then
            fq_aer(NTIX(N))  = 1.
         else
            fq_aer(NTIX(N))  = 0.
        endif
#ifdef TRACERS_NITRATE
         if ( ( TM(l,ntix(n_N_d3)) /trpdens(n_NO3p)) >
     *      (( TM(l,ntix(n))  /trpdens(n)) * 0.03 ) ) then
            fq_aer(NTIX(N))  = 1.
         endif
#endif
      end select
#endif
           CLDINC=CLDSAVT-FCLOUD
          if (fq0.gt.0.and.CLDINC.gt.0.) then
          if(LHX.EQ.LHE) then !liquid cloud
            fq = fq_aer(NTIX(N))*CLDINC
           else ! ice cloud - small dissolution
            fq = fq_aer(NTIX(N))*CLDINC*0.12d0
           endif
          endif
c complete dissolution in convective clouds
c with double dissolution if partially soluble
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
          IF (fq_aer(ntix(n)) > 0. .AND. tr_conv) THEN
#else
          if (TR_CONV) then
#endif
           if (LHX.EQ.LHE) then !liquid cloud
c
               fq=fq_aer(ntix(n))
c              fq=(1.d0+fq_aer(ntix(n)))/2.d0
c              fq=(1.d0+3.d0*fq_aer(ntix(n)))/4.d0
           else
               fq=fq_aer(ntix(n))*0.12d0
c              fq=(1.d0+fq_aer(ntix(n)))/2.d0*0.05d0
c              fq=(1.d0+3.d0*fq_aer(ntix(n)))/4.d0*0.05d0
           endif
          endif
          if (FCLOUD.LT.1.D-16) fq=0.d0
          if (fq.ge.1.d0) fq=0.9999
          if (fq0.eq.0.) fq=0.d0
#endif

        CASE DEFAULT                                ! error
          call stop_model(
     &    'tr_wd_TYPE(NTIX(N)) out of range in GET_COND_FACTOR',255)
      END SELECT
c
      RETURN
      END SUBROUTINE GET_COND_FACTOR

      SUBROUTINE GET_COND_FACTOR_array(
     &     NTX,WMXTR,TEMP,TEMP0,LHX,FCLOUD,
     &     FQ0,fq,TR_CONV,TRWML,TM,THLAW,TR_LEF,pl,ntix,CLDSAVT)
!@sum  GET_COND_FACTOR calculation of condensate fraction for tracers
!@+    within or below convective or large-scale clouds. Gas
!@+    condensation uses Henry's Law if not freezing.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 CLOUDCHCC and CLOUDCHEM subroutines)
! NOTE: THLAW is only computed for the tracers in gases_list!
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only: BYGASC, MAIR,teeny,LHE,tf,by3

      USE TRACER_COM, only :
     &     gases_count,aero_count,water_count,hlawt_count,
! NB: these lists are often used for implicit loops
     &     gases_list,aero_list,water_list,hlawt_list

      USE TRACER_COM, only: tr_RKD,tr_DHD,nWATER,nGAS,nPART,tr_wd_TYPE
     *     ,trname,ntm,t_qlimit,fq_aer,trpdens,n_SO2,n_H2O2,n_H2O2_s
#ifdef TRACERS_SPECIAL_O18
     &     ,supsatfac
#endif
#ifdef TRACERS_HETCHEM
     *     ,trm ,n_SO4_d1, n_SO4_d2, n_SO4_d3,n_SO4
     *     ,n_N_d1,n_N_d2,n_N_d3,n_NO3p, n_Clay,n_Silt1,n_Silt2
      USE MODEL_COM, only  : dtsrc
#endif
      IMPLICIT NONE
C**** Local parameters and variables and arguments:
!@param BY298K unknown meaning for now (assumed= 1./298K)
!@var Ppas pressure at current altitude (in Pascal=kg/s2/m)
!@var TFAC exponential coeffiecient of tracer condensation temperature
!@+   dependence (mole/joule)
!@var FCLOUD fraction of cloud available for tracer condensation
!@var SSFAC dummy variable (assumed units= kg water?)
!@var FQ            fraction of tracer that goes into condensate
!@var FQ0 default fraction of water tracer that goes into condensate
!@var L index for altitude loop
!@var N index for tracer number loop
!@var WMXTR mixing ratio of water available for tracer condensation?
!@var SUPSAT super-saturation ratio for cloud droplets
!@var LHX latent heat flag for whether condensation is to ice or water
!@var RKD dummy variable (= tr_RKD*EXP[ ])
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, RKD,CLDINC,trlef
#ifdef TRACERS_SPECIAL_O18
      real*8 tdegc,alph,fracvs,fracvl,kin_cond_ice
#endif
      REAL*8,  INTENT(IN) :: fq0, FCLOUD, WMXTR, TEMP, TEMP0,LHX
     &     , TR_LEF(NTM), pl,CLDSAVT
      REAL*8,  INTENT(IN), DIMENSION(ntm) :: trwml
      REAL*8,  INTENT(IN), DIMENSION(ntm) :: TM
      REAL*8,  INTENT(OUT):: fq(ntm),thlaw(ntm)
      INTEGER, INTENT(IN) :: NTX, ntix(ntm)
      LOGICAL TR_CONV
      REAL*8 :: FQ0FAC,SUPSAT,SSFAC(NTM),SSFAC0
      INTEGER :: N,IGAS,IAERO,IWAT

c      thlaw(:) = 0. ! default

#if (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_SPECIAL_Shindell) ||\
    (defined TRACERS_AMP)
c
c gases
c
c     cldinc=max(0.,cldsavt-fcloud)
      if(lhx.eq.lhe .and. fcloud.ge.1d-16) then
        Ppas = PL*1.D2          ! pressure to pascals
        tfac = (1.D0/TEMP - BY298K)*BYGASC
        ssfac0 = WMXTR*MAIR*1.D-6*Ppas/(CLDSAVT+teeny)
        ssfac(gases_list) = ssfac0*tr_RKD(gases_list)
        do igas=1,hlawt_count
          n = hlawt_list(igas)
          ssfac(n) = ssfac(n)*exp(-tr_DHD(n)*tfac)
        enddo
        if(tr_conv) then ! convective cloud
          fq0fac = 1.
          if (fq0.eq.0.) fq0fac=0.d0
          do igas=1,gases_count
            n = gases_list(igas)
c not sure why the min was necessary here
c            fq(n) = min(1d0, fq0fac*ssfac(n) / (1d0 + ssfac(n)))
            fq(n) = fq0fac*ssfac(n) / (1d0 + ssfac(n))
            thlaw(n) = 0.
          enddo
        else             ! stratiform cloud
          do igas=1,gases_count
            n = gases_list(igas)
            fq(n) = 0.
c limit gas dissolution to incremental cloud change after cloud forms
c   only apply to non-aqueous sulfur species since this is already
c   done in GET_SULFATE
c but H2O2 should be limited if not coupled with sulfate, haven't done this 
c           if (n.ne.n_h2O2.and.n.ne.n_so2.and.n.ne.n_h2O2_s) then
c           if (FCLOUD.ne.0.) tr_lef(n)=cldinc
c           endif
            trlef=min(tr_lef(n),cldsavt)
            thlaw(n) = min(tm(n),max(0d0,
     &           (ssfac(n)*trlef*tm(n)-TRWML(n))
     &           /(1.D0+ssfac(n)) ))
          enddo
        endif
      else
        fq(gases_list) = 0.
        thlaw(gases_list) = 0.
      endif
#endif /* dissolved gases */

c
c loop over water species
c
      do iwat=1,water_count
        n = water_list(iwat)
        fq(n) = fq0
#ifdef TRACERS_SPECIAL_O18
          if (fq0.gt.0. .and. fq0.lt.1.) then
C**** If process occurs at constant temperature, calculate condensate
C**** in equilibrium with source vapour. Otherwise, use mid-point
C**** temperature and estimate instantaneous fractionation. This gives
C**** a very good estimate to complete integral
C****
            if (abs(temp-temp0).gt.1d-14) then  ! use instantaneous frac
              tdegc=0.5*(temp0 + temp) -tf
C**** Calculate alpha (fractionation coefficient)
                if (LHX.eq.LHE) then ! cond to water
                  alph=1./fracvl(tdegc,ntix(n))
                else            ! cond to ice
                  alph=1./fracvs(tdegc,ntix(n))
C**** kinetic fractionation can occur as a function of supersaturation
C**** this is a parameterisation from Georg Hoffmann
                  supsat=1d0-supsatfac*tdegc
                  if (supsat .gt. 1.) alph=kin_cond_ice(alph,supsat
     *                 ,ntix(n))
                end if
                fq(n) = 1.- (1.-fq0)**alph
            else
C**** assume condensate in equilibrium with vapour at temp
              tdegc=temp -tf
              if (LHX.eq.LHE) then ! cond to water
                alph=1./fracvl(tdegc,ntix(n))
              else              ! cond to ice
                alph=1./fracvs(tdegc,ntix(n))
C**** kinetic fractionation can occur as a function of supersaturation
C**** this is a parameterisation from Georg Hoffmann
                supsat=1d0-supsatfac*tdegc
                if (supsat .gt. 1.) alph=kin_cond_ice(alph,supsat
     *               ,ntix(n))
              end if
              fq(n) = alph * fq0/(1.+(alph-1.)*fq0)
            end if
          else
            fq(n) = fq0
          end if
#endif
      enddo ! end loop over water species

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_COSMO) ||\
    (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP) || (defined TRACERS_RADON)
c
c aerosols
c

c     if (FCLOUD.lt.1.D-16 .or. fq0.eq.0.) then
      if (CLDSAVT.lt.1.D-16 .or. fq0.eq.0.) then
        fq(aero_list) = 0.      ! defaults to zero
      else

#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_DUST) &&\
    (defined TRACERS_HETCHEM)

      n = n_Clay
      if ( ( TM(ntix(n_SO4_d1)) /trpdens(n_SO4)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        fq_aer(NTIX(N))  = 1.
      else
        fq_aer(NTIX(N))  = 0.
      endif
#ifdef TRACERS_NITRATE
      if ( ( TM(ntix(n_N_d1)) /trpdens(n_NO3p)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        fq_aer(NTIX(N))  = 1.
      endif
#endif

      n = n_Silt1
      if ( ( TM(ntix(n_SO4_d2)) /trpdens(n_SO4)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        fq_aer(NTIX(N))  = 1.
      else
        fq_aer(NTIX(N))  = 0.
      endif
#ifdef TRACERS_NITRATE
      if ( ( TM(ntix(n_N_d2)) /trpdens(n_NO3p)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        fq_aer(NTIX(N))  = 1.
      endif
#endif

      n = n_Silt2
      if ( ( TM(ntix(n_SO4_d3)) /trpdens(n_SO4)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        fq_aer(NTIX(N))  = 1.
      else
        fq_aer(NTIX(N))  = 0.
      endif
#ifdef TRACERS_NITRATE
      if ( ( TM(ntix(n_N_d3)) /trpdens(n_NO3p)) >
     *     (( TM(ntix(n))  /trpdens(n)) * 0.03 ) ) then
        fq_aer(NTIX(N))  = 1.
      endif
#endif

#endif

      cldinc=cldsavt-fcloud
      if(tr_conv) then          ! convective cloud
c complete dissolution in convective clouds
c with double dissolution if partially soluble
        if(lhx.eq.lhe) then
          fq(aero_list) = fq_aer(aero_list)
        else
          fq(aero_list) = fq_aer(aero_list)*0.12d0
        endif
c this shouldn't work because cldinc should be fcld for
c    when cloud first forms
      elseif(fq0.gt.0 .and. cldinc.gt.0.) then ! growing stratiform cloud
        if(lhx.eq.lhe) then
          fq(aero_list) = fq_aer(aero_list)*cldinc
        else
          fq(aero_list) = fq_aer(aero_list)*cldinc*0.12d0
        endif
      else
        fq(aero_list) = 0.
      endif
      where(fq(aero_list).ge.1.d0) fq(aero_list)=0.9999
c
c use this code in place of the above if the commented-out formulas
c for (dust?) fq are reinstated
c

c      do iaero=1,aero_count ! loop over aerosols
c        n = aero_list(iaero)
cc complete dissolution in convective clouds
cc with double dissolution if partially soluble
c          if (TR_CONV) then ! convective cloud
c            if (LHX.EQ.LHE) then !liquid cloud
ccdust #if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
ccdust     (defined TRACERS_QUARZHEM)
ccdust           IF (fq_aer(ntix(n)) > 0.)
ccdust #endif
c              fq(n)=fq_aer(ntix(n))
ccdust?              fq(n)=(1.d0+fq_aer(ntix(n)))/2.d0
ccdust?              fq(n)=(1.d0+3.d0*fq_aer(ntix(n)))/4.d0
c            else
ccdust #if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
ccdust     (defined TRACERS_QUARZHEM)
ccdust           IF (fq_aer(ntix(n)) > 0.)
ccdust #endif
c              fq(n)=fq_aer(ntix(n))*0.12d0
ccdust?              fq(n)=(1.d0+fq_aer(ntix(n)))/2.d0*0.05d0
ccdust?              fq(n)=(1.d0+3.d0*fq_aer(ntix(n)))/4.d0*0.05d0
c
c            endif
c          elseif (fq0.gt.0.and.CLDINC.gt.0.) then ! stratiform cloud.
cc only dissolve if the cloud has grown
c            if(LHX.EQ.LHE) then !liquid cloud
c              fq(n) = fq_aer(ntix(n))*CLDINC
c            else                ! ice cloud - small dissolution
c              fq(n) = fq_aer(ntix(n))*CLDINC*0.12d0
c            endif
c          endif
c          if (fq(n).ge.1.d0) fq(n)=0.9999
c      enddo ! end loop over aerosols

      endif ! fcloud>0 and fq0.ne.0

#endif /* aerosols */

      RETURN
      END SUBROUTINE GET_COND_FACTOR_array


      SUBROUTINE GET_WASH_FACTOR(N,b_beta_DT,PREC,fq
     * ,TEMP,LHX,WMXTR,FCLOUD,L,TM,TRPR,THLAW,pl,ntix)
!@sum  GET_WASH_FACTOR calculation of the fraction of tracer
!@+    scavanged by precipitation below convective clouds ("washout").
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 CWASH and WASH_EVAP routines)
c
C**** GLOBAL parameters and variables:
      USE TRACER_COM, only: nWATER, nGAS, nPART, tr_wd_TYPE,
     * tr_RKD,tr_DHD,LM,NTM,rc_washt
#ifdef TRACERS_AEROSOLS_Koch
     * ,trname,n_seasalt1,n_seasalt2
c     USE PBLCOM, only: wsavg
#endif
c      USE CLOUDS, only: NTIX,PL
      USE CONSTANT, only: BYGASC,LHE,MAIR,teeny
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var FQ fraction of tracer scavenged by below-cloud precipitation
!@param rc_wash aerosol washout rate constant (mm-1)
!@var PREC precipitation amount from layer above for washout (mm)
!@var b_beta_DT precipitating grid box fraction from lowest
!@+   percipitating layer.
!@+   The name was chosen to correspond to Koch et al. p. 23,802.
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: N,L,ntix(ntm)
      REAL*8, INTENT(OUT):: FQ,THLAW
      REAL*8, INTENT(IN) :: PREC,b_beta_DT,TEMP,LHX,WMXTR,FCLOUD,
     *  TM(LM,NTM),pl
      REAL*8, PARAMETER :: rc_wash = 1.D-1, BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac, RKD, TRPR(ntm)
C
      thlaw=0.
      SELECT CASE(tr_wd_TYPE(NTIX(N)))
        CASE(nGAS)                            ! gas
          fq = 0.D0                           ! frozen and default case
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_SPECIAL_Shindell) ||\
    (defined TRACERS_AMP)
          IF(LHX.EQ.LHE) THEN                 ! if not frozen then:
            Ppas = PL*1.D2                 ! pressure to pascals
            tfac = (1.D0/TEMP - BY298K)*BYGASC
            ssfac=WMXTR*MAIR*1.D-6*Ppas/(FCLOUD+teeny)
            IF(tr_DHD(NTIX(N)).ne.0.D0) THEN
              ssfac=(ssfac*tr_RKD(NTIX(N)))*DEXP(-tr_DHD(NTIX(N))*tfac)
            ELSE
              ssfac=(ssfac*tr_RKD(NTIX(N)))
            END IF
c            ssfac=RKD*WMXTR*MAIR*1.D-6*Ppas/(FCLOUD+teeny)
            thlaw=(ssfac*tm(l,NTIX(N))-TRPR(NTIX(N)))
     *            /(1.D0+ssfac)
            if (thlaw.lt.0.) thlaw=0.d0
            if (thlaw.gt.tm(l,NTIX(N))) thlaw=tm(l,NTIX(N))
            if (FCLOUD.lt.1.D-16) fq=0.d0
            if (FCLOUD.LT.1.D-16) thlaw=0.d0
          ENDIF
#endif
        CASE(nWATER)                          ! water/original method
          fq = 0.D0
        CASE(nPART)                           ! aerosols
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_COSMO) ||\
    (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_OM_SP) ||\
    (defined SHINDELL_STRAT_EXTRA) || (defined TRACERS_AMP) ||\
    (defined TRACERS_RADON)
          fq = -b_beta_DT*(EXP(-PREC*rc_washt(ntix(n)))-1.D0)
          if (FCLOUD.lt.1.D-16) fq=0.d0
          if (fq.lt.0.) fq=0.d0
c         if (wsavg(i,j).gt.10.and.PREC.gt.0.and.l.eq.1) then
c         select case (trname(n))
c          case('seasalt1')
c          fq=0.
c          case('seasalt2')
c          fq=0.
c         end select
c         endif
#endif
        CASE DEFAULT                          ! error
          call stop_model(
     &    'tr_wd_TYPE(NTIX(N)) out of range in WASHOUT_TRACER',255)
      END SELECT
c
      RETURN
      END SUBROUTINE GET_WASH_FACTOR

      SUBROUTINE GET_WASH_FACTOR_array(NTX,b_beta_DT,PREC,fq
     * ,TEMP,LHX,WMXTR,FCLOUD,TM,TRPR,THLAW,pl,ntix,BELOW_CLOUD)
!@sum  GET_WASH_FACTOR calculation of the fraction of tracer
!@+    scavanged by precipitation below convective clouds ("washout").
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 CWASH and WASH_EVAP routines)
! NOTE: THLAW is only computed for the tracers in gases_list!
! NOTE: FQ is only computed for the tracers in aero_list!
c
C**** GLOBAL parameters and variables:
      USE TRACER_COM, only: nWATER, nGAS, nPART, tr_wd_TYPE,
     * tr_RKD,tr_DHD,NTM,rc_washt
#ifdef TRACERS_AEROSOLS_Koch
     * ,trname,n_seasalt1,n_seasalt2
c     USE PBLCOM, only: wsavg
#endif

      USE TRACER_COM, only :
     &     gases_count,aero_count,water_count,hlawt_count,
! NB: these lists are often used for implicit loops
     &     gases_list,aero_list,water_list,hlawt_list

c      USE CLOUDS, only: NTIX,PL
      USE CONSTANT, only: BYGASC,LHE,MAIR,teeny
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var FQ fraction of tracer scavenged by below-cloud precipitation
!@param rc_wash aerosol washout rate constant (mm-1)
!@var PREC precipitation amount from layer above for washout (mm)
!@var b_beta_DT precipitating grid box fraction from lowest
!@+   percipitating layer.
!@+   The name was chosen to correspond to Koch et al. p. 23,802.
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: NTX,ntix(ntm)
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: THLAW
      REAL*8, INTENT(INOUT), DIMENSION(NTM) :: FQ
      REAL*8, INTENT(IN) :: PREC,b_beta_DT,TEMP,LHX,WMXTR,FCLOUD,
     *  TM(NTM),pl, TRPR(ntm)
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac0, ssfac(ntm), bb_tmp
      INTEGER :: N,IGAS,IAERO
      LOGICAL BELOW_CLOUD
C

c      thlaw(:)=0.

c
c gases
c
c      fq(gases_list) = 0.D0
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_SPECIAL_Shindell) ||\
    (defined TRACERS_AMP)
      if(      LHX.EQ.LHE ! if not frozen
     &   .AND. FCLOUD.GE.1D-16 .AND. WMXTR.GT.0. AND . BELOW_CLOUD) THEN 
        bb_tmp = max(b_beta_DT,0.) ! necessary check?
        Ppas = PL*1.D2          ! pressure to pascals
        tfac = (1.D0/TEMP - BY298K)*BYGASC
        ssfac0 = WMXTR*MAIR*1.D-6*Ppas/(FCLOUD+teeny)
        ssfac(gases_list) = ssfac0*tr_RKD(gases_list)
        do igas=1,hlawt_count
          n = hlawt_list(igas)
          ssfac(n) = ssfac(n)*exp(-tr_DHD(n)*tfac)
        enddo
        do igas=1,gases_count
          n = gases_list(igas)
          thlaw(n) = min( tm(n),max( 0d0,(FCLOUD*
     &               ssfac(n)*tm(n)-TRPR(n))/(1.D0+ssfac(n)) ))
        enddo
      else
        thlaw(gases_list) = 0.
      endif
#endif

c
c water species
c
c      fq(water_list) = 0d0

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_COSMO) ||\
    (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_OM_SP) ||\
    (defined SHINDELL_STRAT_EXTRA) || (defined TRACERS_AMP) ||\
    (defined TRACERS_RADON)
c
c aerosols
c
      if(FCLOUD.GE.1.D-16 .and. prec.gt.0.) then
        bb_tmp = max(b_beta_DT,0.) ! necessary check?
        do iaero=1,aero_count
          n = aero_list(iaero)
          fq(n) = bb_tmp*(1d0-exp(-prec*rc_washt(n)))
c         if (wsavg(i,j).gt.10.and.PREC.gt.0.and.l.eq.1) then
c         select case (trname(n))
c          case('seasalt1')
c          fq(n)=0.
c          case('seasalt2')
c          fq(n)=0.
c         end select
c         endif
        enddo
      else
        fq(aero_list) = 0.
      endif
#endif

      RETURN
      END SUBROUTINE GET_WASH_FACTOR_array

      SUBROUTINE GET_EVAP_FACTOR(N,TEMP,LHX,QBELOW,HEFF,FQ0,fq,ntix)
!@sum  GET_EVAP_FACTOR calculation of the evaporation fraction
!@+    for tracers.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 EVAPD and WASH_EVAP routines)
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only : tf,lhe
      USE TRACER_COM, only: ntm,tr_evap_fact, tr_wd_TYPE,nwater,trname
c      USE CLOUDS, only: NTIX
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var FQ            fraction of tracer evaporated
!@var FQ0 [default] fraction of tracer evaporated
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: N,ntix(ntm)
      REAL*8,  INTENT(OUT):: FQ
      REAL*8,  INTENT(IN) :: FQ0,TEMP,LHX
!@var QBELOW true if evap is occuring below cloud
      LOGICAL, INTENT(IN) :: QBELOW
!@var HEFF effective relative humidity for evap occuring below cloud
      REAL*8, INTENT(IN) :: HEFF
#ifdef TRACERS_SPECIAL_O18
      real*8 tdegc,alph,fracvl,fracvs,kin_evap_prec
#endif
c
      select case (tr_wd_TYPE(NTIX(N)))
      case default
        fq=FQ0*tr_evap_fact(tr_wd_TYPE(NTIX(N)))
        if(FQ0.ge.1.) fq=1.D0 ! total evaporation
c
      case (nWater)
#ifdef TRACERS_SPECIAL_O18
          tdegc=temp-tf
          if (lhx.eq.lhe) then
            alph=fracvl(tdegc,ntix(n))
C**** below clouds kinetic effects with evap into unsaturated air
            if (QBELOW.and.heff.lt.1.) alph=kin_evap_prec(alph,heff
     *           ,ntix(n))
          else
C**** no fractionation for ice evap
            alph=1.
          end if
          if (fq0.ne.1.) then
             fq = 1. - (1.-fq0)**alph
          else
            fq = fq0
          end if
#else
          fq=FQ0*tr_evap_fact(tr_wd_TYPE(NTIX(N)))
          if(FQ0.ge.1.) fq=1.D0 ! total evaporation
#endif
      end select
      RETURN
      END SUBROUTINE GET_EVAP_FACTOR

      SUBROUTINE GET_EVAP_FACTOR_array(
     &     NTX,TEMP,LHX,QBELOW,HEFF,FQ0,fq,ntix)
!@sum  GET_EVAP_FACTOR calculation of the evaporation fraction
!@+    for tracers.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 EVAPD and WASH_EVAP routines)
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only : tf,lhe
      USE TRACER_COM, only: ntm,tr_evap_fact, tr_wd_TYPE,nwater,trname
     &     ,water_count,water_list
c      USE CLOUDS, only: NTIX
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var FQ            fraction of tracer evaporated
!@var FQ0 [default] fraction of tracer evaporated
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: NTX,ntix(ntm)
      REAL*8,  INTENT(OUT):: FQ(NTM)
      REAL*8,  INTENT(IN) :: FQ0,TEMP,LHX
!@var QBELOW true if evap is occuring below cloud
      LOGICAL, INTENT(IN) :: QBELOW
!@var HEFF effective relative humidity for evap occuring below cloud
      REAL*8, INTENT(IN) :: HEFF
#ifdef TRACERS_SPECIAL_O18
      real*8 tdegc,alph,fracvl,fracvs,kin_evap_prec
      integer :: iwat,n
#endif
c

      if(fq0.ge.1d0) then
        fq(1:ntx) = 1d0 ! total evaporation
      else
        fq(1:ntx) = fq0*tr_evap_fact(tr_wd_type(ntix(1:ntx)))
      endif

#ifdef TRACERS_SPECIAL_O18
c overwrite fq for water isotopes
      tdegc=temp-tf
      do iwat=1,water_count
        n = water_list(iwat)
        if (lhx.eq.lhe) then
          alph=fracvl(tdegc,ntix(n))
C**** below clouds kinetic effects with evap into unsaturated air
          if (QBELOW.and.heff.lt.1.)
     &         alph=kin_evap_prec(alph,heff,ntix(n))
        else
C**** no fractionation for ice evap
          alph=1.
        end if
        if (fq0.ne.1.) then
          fq(n) = 1. - (1.-fq0)**alph
        else
          fq(n) = fq0
        end if
      enddo
#endif

      RETURN
      END SUBROUTINE GET_EVAP_FACTOR_array

#endif

#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_AMP)
      SUBROUTINE GET_SULF_GAS_RATES
!@sum  GET_SULF_GAS_RATES calculation of rate coefficients for
!@+    gas phase sulfur oxidation chemistry
!@auth Bell
!@ver  1.0
      USE MODEL_COM, only: im,jm,lm,t,ls1
      USE DOMAIN_DECOMP_ATM, only : GRID, GET, write_parallel
      USE DYNAMICS, only: pmid,am,pk,LTROPO
      USE GEOM, only: axyp,imaxj
      USE TRACER_COM, only: rsulf1,rsulf2,rsulf3,rsulf4
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: which_trop
#endif
      INTEGER J_0, J_1, I_0, I_1
      real*8 ppres,te,tt,mm,dmm,rk4,ek4,f
CGreg: certain things now done outside the loops for speed:
      real*8, parameter ::  a= 73.41463d20, ! 6.02d20/.082d0
     *     aa=1.d-20,
     *     b= 0.357d-22,         ! 1.7d-22*0.21d0*1.d-20/aa
     *     c= 1.155d-11,         ! 5.5d-20*0.21d0*1.d-11/aa
     *     d= 4.0d-11            ! 4.0d-20*1.d-11/aa

#ifdef TRACERS_SPECIAL_Shindell
!@var maxl chosen tropopause 0=LTROPO(I,J), 1=LS1-1
      integer maxl
#endif

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C Initialise
      rsulf1(:,j_0:j_1,:)=0.d0
      rsulf2(:,j_0:j_1,:)=0.d0
      rsulf3(:,j_0:j_1,:)=0.d0
      rsulf4(:,j_0:j_1,:)=0.d0

C Reactions
C***1.DMS + OH -> 0.75SO2 + 0.25MSA
C***2.DMS + OH -> SO2
C***3.DMS + NO3 -> HNO3 + SO2
C***4.SO2 + OH -> SO4 + HO2


      do l=1,LM
      do j=J_0,J_1
      do i=I_0,imaxj(j)
c
      maxl = ltropo(i,j)
#ifdef TRACERS_SPECIAL_Shindell
      if(which_trop.eq.1)maxl=ls1-1
#endif
      if(l.le.maxl) then

C Calculate effective temperature

        ppres=pmid(l,i,j)*9.869d-4 !in atm
        te=pk(l,i,j)*t(i,j,l)
        mm=am(l,i,j)*axyp(i,j)
        tt = 1.d0/te

c DMM is number density of air in molecules/cm3

        dmm=ppres*tt*a
        rsulf1(i,j,l) =
     & b*dmm*exp(7810.d0*tt)*aa/(1.d0+c*exp(7460.d0*tt)*dmm*aa)

        rsulf2(i,j,l) = 9.6d-12*exp(-234.d0*tt)

        rsulf3(i,j,l) = 1.9d-13*exp(520.d0*tt)

        rk4 = aa*((tt*300.d0)**(3.3d0))*dmm*d
        f=log10(0.5d12*rk4)
        ek4 = 1.d0/(1.d0 + (f*f))

        rsulf4(i,j,l) = (rk4/(1.d0 + 0.5d12*rk4  ))*(0.45d0**ek4)

      endif

      end do
      end do
      end do

      END SUBROUTINE GET_SULF_GAS_RATES
#endif
