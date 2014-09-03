#include "rundeck_opts.h"

      MODULE FLUXES
!@sum  FLUXES contains the fluxes between various atm-grid components
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
      USE DOMAIN_DECOMP_ATM, ONLY : grid

#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm
#ifndef SKIP_TRACER_SRCS
     *     ,ntsurfsrcmax,nt3Dsrcmax
#endif
#endif
      IMPLICIT NONE

!@var RUNOSI run off from sea/lake ice after surface (kg/m^2)
!@var ERUNOSI energy of run off from sea/lake ice after surface (J/m^2)
!@var SRUNOSI salt in run off from sea/lake ice after surface (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RUNOSI, ERUNOSI, SRUNOSI
!@var RUNPSI run off from sea/lake ice after precip (kg/m^2)
!@var ERUNPSI energy of run off from sea/lake ice after precip (J/m^2)
!@var SRUNPSI salt in run off from sea/lake ice after precip (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RUNPSI, SRUNPSI, ERUNPSI
!@var RUNOE run off from earth (kg/m^2)
!@var ERUNOE energy of run off from earth (J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RUNOE, ERUNOE
C**** DMSI,DHSI,DSSI are fluxes for ice formation within water column
!@var DMSI mass flux of sea ice 1) open water and 2) under ice (kg/m^2)
!@var DHSI energy flux of sea ice 1) open water and 2) under ice (J/m^2)
!@var DSSI salt flux in sea ice 1) open water and 2) under ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DMSI, DHSI, DSSI
!@var fmsi_io,fhsi_io,fssi_io basal ice-ocean fluxes (kg or J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: fmsi_io,fhsi_io,fssi_io
!@var RUNOLI run off from land ice (kg/m^2) (Energy always=0)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RUNOLI

C**** surface energy fluxes defined over type
!@param NSTYPE number of surface types for radiation purposes
      INTEGER, PARAMETER :: NSTYPE=4
!@var E0 net energy flux at surface for each type (J/m^2)
!@var E1 net energy flux at layer 1 for each type (J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: E0,E1
!@var EVAPOR evaporation over each type (kg/m^2) 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: EVAPOR
!@var SOLAR absorbed solar radiation (J/m^2)
!@+   SOLAR(1)  absorbed by open water
!@+   SOLAR(2)  absorbed by ice
!@+   SOLAR(3)  absorbed by water under the ice
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SOLAR

C**** Momemtum stresses are calculated as if they were over whole box
!@var DMUA,DMVA momentum flux from atmosphere for each type (kg/m s) 
!@+   On atmospheric A grid (tracer point)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DMUA,DMVA
!@var UI2rho Ustar*2*rho ice-ocean friction velocity on atmospheric grid
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: UI2rho
!@var OGEOZA ocean surface height geopotential (m^2/s^2) (on ATM grid)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: OGEOZA

C**** currently saved - should be replaced by fluxed quantities
!@var DTH1,DQ1 heat/water flux from atmos. summed over type (C, kg/kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DTH1,DQ1

!@var uflux1 surface turbulent u-flux (=-<uw>) 
!@var vflux1 surface turbulent v-flux (=-<vw>)
!@var tflux1 surface turbulent t-flux (=-<tw>)
!@var qflux1 surface turbulent q-flux (=-<qw>)
      real*8, allocatable, dimension(:,:) :: 
     &        uflux1,vflux1,tflux1,qflux1

C**** The E/FLOWO, E/S/MELTI, E/GMELT arrays are used to flux quantities 
C**** to the ocean that are not tied to the open water/ice covered 
C**** fractions. This is done separately for river flow, complete
C**** sea ice melt and iceberg/glacial melt.
!@var FLOWO,EFLOWO mass, energy from rivers into ocean (kg, J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FLOWO,EFLOWO
!@var MELTI,EMELTI,SMELTI mass,energy,salt from simelt into ocean (kg,J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: MELTI,EMELTI,SMELTI
!@var GMELT,EGMELT mass,energy from glacial melt into ocean (kg,J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: GMELT,EGMELT

!@var PREC precipitation (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PREC
!@var EPREC energy of preciptiation (J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: EPREC
!@var PRECSS precipitation from super-saturation (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PRECSS

#ifdef IRRIGATION_ON
!@var Actual irrigation rate (& energy) [m/s] [W/m2]
      REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:) :: irrig_water_act
      REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:) :: irrig_energy_act
#ifdef TRACERS_WATER
!@var Actual irrigation tracer rate [kg/s] 
      REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: irrig_tracer_act
#endif
#endif

!@var GTEMP ground temperature (upper two levels) over surface type (C)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: GTEMP
!@var GTEMPR radiative ground temperature over surface type (K)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GTEMPR
!@var SSS sea surface salinity on atmospheric grid (ppt)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SSS
!@var MLHC ocean mixed layer heat capacity (J/m^2 C) 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: MLHC
!@var UOSURF, VOSURF ocean surface velocity (Atm A grid) (m/s)
!@+   components defined along true N/S and E/W directions
!@+   At the NP, U points from 90E to 90W, V from IDL to GM
!@+   At the SP, U points from 9OW to 90E, V from GM to IDL 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: UOSURF,VOSURF
!@var UISURF, VISURF dynamic ice surface velocity (Atm A grid) (m/s)
!@+   directions as for UOSURF/VOSURF
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: UISURF,VISURF
!@var APRESS total atmos + sea ice pressure (at base of sea ice) (Pa)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: APRESS
!@var FWSIM fresh water sea ice mass (kg/m^2) (used for qflux model)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FWSIM
!@var MSICNV fresh water sea ice mass convergence after advsi (kg/m^2)
!@+   (used for qflux model)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: MSICNV

C**** fluxes associated with variable lake fractions
!@var DMWLDF  water deficit over land surface (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DMWLDF
!@var DGML energy associated with DMWLDF (J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DGML

!@var bare_soil_wetness bare_soil_wetness (1)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: bare_soil_wetness

#if (defined CHL_from_SeaWIFs) || (defined TRACERS_OceanBiology)
C**** array of Chlorophyll data for use in ocean albedo calculation
!@var CHL Chlorophyll concentration data (mgr/m**3)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: CHL
#endif

#ifdef TRACERS_ON
!@var TRSOURCE non-interactive surface sources/sinks for tracers (kg/s)
#ifndef SKIP_TRACER_SRCS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: trsource
#endif
!@var TRSRFFLX interactive surface sources/sinks for tracers (kg/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: trsrfflx
!@var TRFLUX1 total surface flux for each tracer (kg/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: trflux1
!@var GTRACER ground concentration of tracer on atmospheric grid (kg/kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:):: GTRACER
!@var TR3DSOURCE 3D sources/sinks for tracers (kg/s)
#ifndef SKIP_TRACER_SRCS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:):: tr3Dsource
#endif

#ifdef TRACERS_GASEXCH_ocean
!@var TRGASEX  tracer gas exchange over each type (mol,CO2/m^2/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRGASEX
#endif

#ifdef TRACERS_WATER
!@var TRPREC tracers in precip (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: TRPREC
!@var TREVAPOR tracer evaporation over each type (kg/m^2) 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TREVAPOR
!@var TRUNPSI tracer in run off from sea/lake ice after precip (kg/m^2)
!@var TRUNOSI tracer in run off from sea/lake ice after surface (kg/m^2)
!@var TRUNOE tracer runoff from earth (kg/m^2)
!@var TRUNOLI tracer runoff from land ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: 
     &        TRUNPSI, TRUNOSI, TRUNOE, TRUNOLI

!@var TRFLOWO tracer in river runoff into ocean (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRFLOWO
!@var TRMELTI tracer from simelt into ocean (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRMELTI

C**** fluxes associated with variable lake fractions
!@var DTRL tracers associate with DMWLDF (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DTRL

#ifdef TRACERS_OCEAN
!@var TRGMELT tracer from glacial melt into ocean (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRGMELT
#endif
!@var ftrsi_io ice-ocean tracer fluxes under ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ftrsi_io
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
!@var trprec_dust dust/mineral tracers in precip [kg]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:):: trprec_dust
#endif

#ifdef TRACERS_DRYDEP
!@var TRDRYDEP tracer dry deposition by type (kg/m^2) (positive down)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRDRYDEP 
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
!@var pprec precipitation at previous time step [kg/m^2]
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: pprec
!@var pevap evaporation at previous time step [kg/m^2]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: pevap
!@var dust_flux_glob global array of dust emission flux [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: dust_flux_glob
!@var dust_flux2_glob global array of cubic dust emission flux (for diags only)
!@+   [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: dust_flux2_glob
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) 
#ifdef TRACERS_DRYDEP
!@var depo_turb_glob global array of flux due to dry turb. dep. of tracers
!@+   [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: depo_turb_glob
!@var depo_grav_glob global array of flux due to gravit. settling of tracers
!@+   [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: depo_grav_glob
#endif
#endif

#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
!@var DTRSI tracer flux in sea ice under ice and on open water (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: DTRSI
#endif
#endif

      END MODULE FLUXES

      SUBROUTINE ALLOC_FLUXES(grd_dum)
!@sum   Initializes FLUXES''s arrays
!@auth  Rosalinda de Fainchtein
!@ver  1.0
      USE CONSTANT, only : tf
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
      USE FLUXES
#ifdef TRACERS_ON
      USE tracer_com,ONLY : Ntm
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
     &     ,Ntm_dust
#endif
#endif
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO

      !I-J arrays
      ALLOCATE( RUNOSI  ( I_0H:I_1H , J_0H:J_1H ), 
     &          ERUNOSI ( I_0H:I_1H , J_0H:J_1H ), 
     &          SRUNOSI ( I_0H:I_1H , J_0H:J_1H ),
     &          RUNPSI  ( I_0H:I_1H , J_0H:J_1H ), 
     &          SRUNPSI ( I_0H:I_1H , J_0H:J_1H ),
     &          ERUNPSI ( I_0H:I_1H , J_0H:J_1H ),
     &          RUNOE   ( I_0H:I_1H , J_0H:J_1H ), 
     &          ERUNOE  ( I_0H:I_1H , J_0H:J_1H ),
     &          fmsi_io ( I_0H:I_1H , J_0H:J_1H ),
     &          fhsi_io ( I_0H:I_1H , J_0H:J_1H ),
     &          fssi_io ( I_0H:I_1H , J_0H:J_1H ),
     &          RUNOLI  ( I_0H:I_1H , J_0H:J_1H ),
     &          UI2rho  ( I_0H:I_1H , J_0H:J_1H ),
     &          OGEOZA  ( I_0H:I_1H , J_0H:J_1H ),
     &          DTH1    ( I_0H:I_1H , J_0H:J_1H ),
     &          DQ1     ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT=IER )
      ALLOCATE( uflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          vflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          tflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          qflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          FLOWO   ( I_0H:I_1H , J_0H:J_1H ),
     &          EFLOWO  ( I_0H:I_1H , J_0H:J_1H ),
     &          MELTI   ( I_0H:I_1H , J_0H:J_1H ),
     &          EMELTI  ( I_0H:I_1H , J_0H:J_1H ),
     &          SMELTI  ( I_0H:I_1H , J_0H:J_1H ),
     &          GMELT   ( I_0H:I_1H , J_0H:J_1H ),
     &          EGMELT  ( I_0H:I_1H , J_0H:J_1H ),
     &          PREC    ( I_0H:I_1H , J_0H:J_1H ),
     &          EPREC   ( I_0H:I_1H , J_0H:J_1H ),
     &          PRECSS  ( I_0H:I_1H , J_0H:J_1H ),
     &          SSS     ( I_0H:I_1H , J_0H:J_1H ),
     &          MLHC    ( I_0H:I_1H , J_0H:J_1H ),
     &          UOSURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          VOSURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          UISURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          VISURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          APRESS  ( I_0H:I_1H , J_0H:J_1H ),
     &          FWSIM   ( I_0H:I_1H , J_0H:J_1H ),
     &          MSICNV  ( I_0H:I_1H , J_0H:J_1H ),
     &          DMWLDF  ( I_0H:I_1H , J_0H:J_1H ),
     &          DGML    ( I_0H:I_1H , J_0H:J_1H ),
     &          bare_soil_wetness( I_0H:I_1H , J_0H:J_1H ),
#ifdef IRRIGATION_ON
     &          irrig_water_act ( I_0H:I_1H , J_0H:J_1H ),
     &          irrig_energy_act( I_0H:I_1H , J_0H:J_1H ),
#ifdef TRACERS_WATER
     &          irrig_tracer_act(NTM, I_0H:I_1H , J_0H:J_1H),
#endif
#endif
#if (defined CHL_from_SeaWIFs) || (defined TRACERS_OceanBiology)
     &          CHL     ( I_0H:I_1H , J_0H:J_1H ),
#endif
     &   STAT=IER)


       !n-I-J arrays
       ALLOCATE( DMSI    (  2  , I_0H:I_1H , J_0H:J_1H ), 
     &           DHSI    (  2  , I_0H:I_1H , J_0H:J_1H ), 
     &           DSSI    (  2  , I_0H:I_1H , J_0H:J_1H ),
     &           SOLAR   (  3  , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)


      !I-J-: arrays
      ALLOCATE( E0      ( I_0H:I_1H , J_0H:J_1H , NSTYPE ),
     &          E1      ( I_0H:I_1H , J_0H:J_1H , NSTYPE ),
     &          EVAPOR  ( I_0H:I_1H , J_0H:J_1H , NSTYPE ),
     &          DMUA    ( I_0H:I_1H , J_0H:J_1H , NSTYPE ),
     &          DMVA    ( I_0H:I_1H , J_0H:J_1H , NSTYPE ),
     &   STAT = IER)


       !:,:,I,J array
       ALLOCATE( GTEMP( 2 , NSTYPE, I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
       ALLOCATE( GTEMPR( NSTYPE, I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
       GTEMP=0.    ! initialize at 0 C
       GTEMPR=TF   ! initialize at 273 K

!TRACERS_ON**********************

#ifdef TRACERS_ON
      !(I,J,:,:)  array
#ifndef SKIP_TRACER_SRCS
      ALLOCATE(trsource (I_0H:I_1H,J_0H:J_1H,ntsurfsrcmax,NTM)
     &  ,STAT = IER)
      trsource = 0.
#endif

      !(I,J,:) arrays
      ALLOCATE( trsrfflx( I_0H:I_1H , J_0H:J_1H , NTM    ),
     &          trflux1 ( I_0H:I_1H , J_0H:J_1H , NTM    ),
     &   STAT = IER)

      !I-J-L-:-: array
#ifndef SKIP_TRACER_SRCS
      ALLOCATE( tr3Dsource(I_0H:I_1H,J_0H:J_1H,LM,nt3Dsrcmax,NTM)
     &  ,STAT = IER)
#endif

      !:-:-I-J arrays
      ALLOCATE( GTRACER ( NTM , NSTYPE , I_0H:I_1H , J_0H:J_1H ),
     &  STAT = IER)
      GTRACER=0.

#ifdef TRACERS_GASEXCH_ocean

      ALLOCATE( TRGASEX( NTM , NSTYPE , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
      TRGASEX=0.     !initialize to zero

#endif

#ifdef TRACERS_WATER
                                                    !(:)-(:)-I-J arrays
      ALLOCATE( TREVAPOR( NTM , NSTYPE , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)


       !:-I-J arrays
       ALLOCATE( TRPREC  ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           TRUNPSI ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           TRUNOSI ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           TRUNOE  ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           TRUNOLI ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           TRFLOWO ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           TRMELTI ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           ftrsi_io( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           DTRL    ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
#ifdef TRACERS_OCEAN
       ALLOCATE( TRGMELT ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
#endif
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      ALLOCATE(trprec_dust(Ntm_dust,I_0H:I_1H ,J_0H:J_1H),STAT=ier)
#endif

#ifdef TRACERS_DRYDEP
       ALLOCATE(TRDRYDEP( NTM , NSTYPE , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
         TRDRYDEP = 0.   !Initialize to 0.
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      ALLOCATE(pprec(I_0H:I_1H,J_0H:J_1H),STAT = IER)
      ALLOCATE(pevap(I_0H:I_1H,J_0H:J_1H,NSTYPE),STAT = IER)
      ALLOCATE(dust_flux_glob(I_0H:I_1H,J_0H:J_1H,Ntm_dust),STAT = IER)
      ALLOCATE(dust_flux2_glob(I_0H:I_1H,J_0H:J_1H,Ntm_dust),STAT = IER)
#ifdef TRACERS_DRYDEP
      ALLOCATE(depo_turb_glob(I_0H:I_1H,J_0H:J_1H,Nstype,Ntm)
     &     ,STAT = IER)
      ALLOCATE(depo_grav_glob(I_0H:I_1H,J_0H:J_1H,Nstype,Ntm)
     &     ,STAT = IER)
#endif
#endif

#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      ALLOCATE( DTRSI( NTM ,    2   , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
#endif
#endif

      END SUBROUTINE ALLOC_FLUXES
