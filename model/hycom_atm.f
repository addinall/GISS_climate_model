#include "rundeck_opts.h"

      module hycom_atm
!@sum module for atmospheric variables to be passed to/from hylom.
!@+   hycom will see them as global arrays
      USE HYCOM_DIM, only : aI_0, aI_1, aI_0H, aI_1H
      USE HYCOM_DIM, only : aJ_0, aJ_1, aJ_0H, aJ_1H
      use FLUXES, only: NSTYPE

      USE FLUXES, only : PREC_loc => PREC
      USE FLUXES, only : EVAPOR_loc => EVAPOR
      USE FLUXES, only : FLOWO_loc => FLOWO
      USE FLUXES, only : GMELT_loc => GMELT
      USE FLUXES, only : MELTI_loc => MELTI
      USE FLUXES, only : RUNOSI_loc => RUNOSI
      USE FLUXES, only : RUNPSI_loc => RUNPSI
      USE FLUXES, only : E0_loc => E0
      USE FLUXES, only : EPREC_loc => EPREC
      USE FLUXES, only : EFLOWO_loc => EFLOWO
      USE FLUXES, only : EGMELT_loc => EGMELT
      USE FLUXES, only : EMELTI_loc => EMELTI
      USE FLUXES, only : ERUNOSI_loc => ERUNOSI
      USE FLUXES, only : ERUNPSI_loc => ERUNPSI
      USE FLUXES, only : SRUNOSI_loc => SRUNOSI
      USE FLUXES, only : SRUNPSI_loc => SRUNPSI
      USE FLUXES, only : SMELTI_loc => SMELTI
      USE FLUXES, only : DMUA_loc => DMUA
      USE FLUXES, only : DMVA_loc => DMVA
      USE ICEDYN, only : grid_icdyn
#ifndef CUBED_SPHERE
      USE ICEDYN_COM, only : pack_i2a ! to send dmui,dmvi to the atm grid
      USE ICEDYN_COM, only : pack_a2i
#endif
      USE ICEDYN_COM, only : DMUI_loc => DMUI
      USE ICEDYN_COM, only : DMVI_loc => DMVI
      USE FLUXES, only : SOLAR_loc => SOLAR
      USE FLUXES, only : SSS_loc => SSS
      USE FLUXES, only : UOSURF_loc => UOSURF
      USE FLUXES, only : VOSURF_loc => VOSURF
      USE ICEDYN_COM, only : UOSURF_4DYNSI_loc => UOSURF
      USE ICEDYN_COM, only : VOSURF_4DYNSI_loc => VOSURF
      USE FLUXES, only : OGEOZA_loc => OGEOZA
      USE FLUXES, only : GTEMP_loc => GTEMP
      USE FLUXES, only : GTEMPR_loc => GTEMPR
      USE FLUXES, only : DMSI_loc => DMSI
      USE FLUXES, only : DHSI_loc => DHSI
      USE FLUXES, only : DSSI_loc => DSSI

      USE SEAICE_COM, only : RSI_loc => RSI ! seems to be used for diags only?
      USE MODEL_COM, only : FOCEAN_loc => FOCEAN

#ifdef TRACERS_GASEXCH_ocean
      USE TRACER_COM, only : ntm
      USE FLUXES, only : GTRACER_loc => GTRACER, TRGASEX_loc => TRGASEX
#endif
#ifdef TRACERS_OceanBiology
      USE RAD_COM, only : COSZ1_loc => COSZ1
      USE PBLCOM, only : wsavg_loc => wsavg
#endif
#ifdef OBIO_RAD_coupling
      USE RAD_COM, only: FSRDIR_loc => FSRDIR
     .                  ,FSRDIF_loc => FSRDIF
     .                  ,DIRNIR_loc => DIRNIR
     .                  ,DIFNIR_loc => DIFNIR
     .                  ,SRVISSURF_loc => SRVISSURF
#endif
#ifdef pCO2_ONLINE
      USE PBL_DRV, only : t_pbl_args
#endif


      implicit none
      private

      public PREC_loc
      public EVAPOR_loc
      public FLOWO_loc
      public GMELT_loc
      public MELTI_loc
      public RUNOSI_loc
      public RUNPSI_loc
      public E0_loc
      public EPREC_loc
      public EFLOWO_loc
      public EGMELT_loc
      public EMELTI_loc
      public ERUNOSI_loc
      public ERUNPSI_loc
      public SRUNOSI_loc
      public SRUNPSI_loc
      public SMELTI_loc
      public DMUA_loc
      public grid_icdyn
#ifndef CUBED_SPHERE
      public pack_i2a
      public pack_a2i
#endif
      public DMUI_loc
      public DMVA_loc
      public DMVI_loc
      public SOLAR_loc
      public SSS_loc
      public UOSURF_loc
      public VOSURF_loc
      public UOSURF_4DYNSI_loc
      public VOSURF_4DYNSI_loc
      public OGEOZA_loc
      public GTEMP_loc
      public GTEMPR_loc
      public DMSI_loc
      public DHSI_loc
      public DSSI_loc

      public RSI_loc
      public FOCEAN_loc

      public alloc_hycom_atm

      public ataux_loc,atauy_loc,aflxa2o_loc
     .     ,aemnp_loc,aice_loc,asalt_loc
     .     ,austar_loc,aswflx_loc
     .     ,admui_loc,admvi_loc
     .     ,asst_loc,atempr_loc

#ifdef TRACERS_GASEXCH_ocean
      public atracflx_loc
      public GTRACER_loc
      public TRGASEX_loc
#endif
#ifdef TRACERS_OceanBiology
      public asolz_loc
      public awind_loc
      public wsavg_loc
      public COSZ1_loc
#endif
#ifdef OBIO_RAD_coupling
      public avisdir_loc
      public avisdif_loc
      public anirdir_loc
      public anirdif_loc
      public FSRDIR_loc
      public FSRDIF_loc
      public DIRNIR_loc
      public DIFNIR_loc
      public SRVISSURF_loc
#endif

      ! accumulators for output on atmospheric grid
      ! (shouldn't these actually be accumulated on ocean grid
      !  and then interpolated to atmospheric grid?)
      real*8, allocatable, dimension(:,:) :: ataux_loc,atauy_loc
     .     ,aflxa2o_loc,aemnp_loc,aice_loc,asalt_loc
     .     ,austar_loc,aswflx_loc
     .     ,admui_loc,admvi_loc ! == dmui_loc,dmvi_loc on atm. domain
     .     ,asst_loc,atempr_loc

#ifdef TRACERS_GASEXCH_ocean
      real, ALLOCATABLE, DIMENSION(:,:,:) :: atracflx_loc
#endif
#ifdef TRACERS_OceanBiology
      real, ALLOCATABLE, DIMENSION(:,:)    :: asolz_loc
!wind speed from modelE (see hycom2.f)
      real, ALLOCATABLE, DIMENSION(:,:)    :: awind_loc
#endif
#ifdef OBIO_RAD_coupling
      real*8, ALLOCATABLE, DIMENSION(:,:) :: avisdir_loc,avisdif_loc
     .     ,anirdir_loc,anirdif_loc
#endif


      contains

      subroutine alloc_hycom_atm

      ALLOCATE(
     &     ataux_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     atauy_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aflxa2o_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aemnp_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aice_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     asalt_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     austar_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aswflx_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     asst_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     atempr_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     admui_loc(aI_0H:aI_1H,aJ_0H:aJ_1H), ! temporary
     &     admvi_loc(aI_0H:aI_1H,aJ_0H:aJ_1H)  ! temporary
     &     )

#ifdef TRACERS_GASEXCH_ocean
      ALLOCATE( atracflx_loc(aI_0H:aI_1H,aJ_0H:aJ_1H,ntm) )
#endif
#ifdef TRACERS_OceanBiology
      ALLOCATE(asolz_loc(aI_0H:aI_1H,aJ_0H:aJ_1H))
      ALLOCATE(awind_loc(aI_0H:aI_1H,aJ_0H:aJ_1H))
#endif
#ifdef OBIO_RAD_coupling
      ALLOCATE(avisdir_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &         avisdif_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &         anirdir_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &         anirdif_loc(aI_0H:aI_1H,aJ_0H:aJ_1H) )
#endif

      end subroutine alloc_hycom_atm

cddd      subroutine alloc_locals
cddd
cddd      ALLOCATE( PREC_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( EVAPOR_loc( I_0H:I_1H , J_0H:J_1H , NSTYPE ) )
cddd      ALLOCATE( FLOWO_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( GMELT_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( MELTI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( RUNOSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( RUNPSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( E0_loc( I_0H:I_1H , J_0H:J_1H , NSTYPE ) )
cddd      ALLOCATE( EPREC_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( EFLOWO_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( EGMELT_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( EMELTI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( ERUNOSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( ERUNPSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( SRUNOSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( SRUNPSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( SMELTI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( DMUA_loc( I_0H:I_1H , J_0H:J_1H , NSTYPE ) )
cddd      ALLOCATE( DMUI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( DMVA_loc( I_0H:I_1H , J_0H:J_1H, NSTYPE ) )
cddd      ALLOCATE( DMVI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( SOLAR_loc(  3  , I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( SSS_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( UOSURF_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( VOSURF_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( OGEOZA_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( GTEMP_loc( 2 , NSTYPE, I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( GTEMPR_loc( NSTYPE, I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( DMSI_loc(  2  , I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( DHSI_loc(  2  , I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( DSSI_loc(  2  , I_0H:I_1H , J_0H:J_1H ) )
cddd
cddd
cddd      end subroutine alloc_locals

      end module hycom_atm
