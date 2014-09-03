#include "rundeck_opts.h"

      MODULE OFLUXES
!@sum  OFLUXES contains the fluxes between various components
!       that are going to be used on ocean grid
!@auth Larissa Nazarenko
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
      USE OCEANRES,  only : imo,jmo
!      USE DOMAIN_DECOMP_1D, ONLY : grid
      USE OCEANR_DIM, only : grid=>ogrid

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: ntm
#ifdef TRACERS_GASEXCH_ocean
      USE TRACER_COM, only: ntm_gasexch
#endif
#endif

      IMPLICIT NONE

!@var SOLAR absorbed solar radiation (J/m^2)
!@+   SOLAR(1)  absorbed by open water
!@+   SOLAR(2)  absorbed by ice
!@+   SOLAR(3)  absorbed by water under the ice
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oSOLAR
!@param NSTYPE number of surface types for radiation purposes
      INTEGER, PARAMETER :: NSTYPE=4
!@var E0 net energy flux at surface for each type (J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oE0
!@var EVAPOR evaporation over each type (kg/m^2) 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oEVAPOR
C**** DMSI,DHSI,DSSI are fluxes for ice formation within water column
!@var DMSI mass flux of sea ice 1) open water and 2) under ice (kg/m^2)
!@var DHSI energy flux of sea ice 1) open water and 2) under ice (J/m^2)
!@var DSSI salt flux in sea ice 1) open water and 2) under ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oDMSI, oDHSI, oDSSI
!@var RUNOSI run off from sea/lake ice after surface (kg/m^2)
!@var ERUNOSI energy of run off from sea/lake ice after surface (J/m^2)
!@var SRUNOSI salt in run off from sea/lake ice after surface (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oRUNOSI, oERUNOSI, oSRUNOSI
!@var FLOWO,EFLOWO mass, energy from rivers into ocean (kg, J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oFLOWO, oEFLOWO
!@var APRESS total atmos + sea ice pressure (at base of sea ice) (Pa)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oAPRESS
!@var MELTI,EMELTI,SMELTI mass,energy,salt from simelt into ocean (kg,J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oMELTI, oEMELTI, oSMELTI

!@var DMUA,DMVA momentum flux from atmosphere for each type (kg/m s) 
!@+   On OCN A grid (tracer point)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oDMUA, oDMVA
!@var DMUI,DMVI momentum flux from sea ice to ocean (kg/m s)
!@+   On OCN C grid 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oDMUI, oDMVI

!@var GMELT,EGMELT mass,energy from glacial melt into ocean (kg,J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oGMELT, oEGMELT

!@var RUNPSI run off from sea/lake ice after precip (kg/m^2)
!@var ERUNPSI energy of run off from sea/lake ice after precip (J/m^2)
!@var SRUNPSI salt in run off from sea/lake ice after precip (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oRUNPSI, oERUNPSI, oSRUNPSI   
!@var PREC precipitation (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oPREC
!@var EPREC energy of preciptiation (J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oEPREC

!@var RSI fraction of open water area covered in ice
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oRSI

#ifdef TRACERS_OCEAN

#ifdef TRACERS_GASEXCH_ocean
!@var TRGASEX  tracer gas exchange over each type (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: oTRGASEX
#endif

#ifdef TRACERS_WATER
!@var TRFLOWO tracer in river runoff into ocean (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oTRFLOWO
!@var TREVAPOR tracer evaporation over each type (kg/m^2) 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: oTREVAPOR
!@var TRUNPSI tracer in run off from sea/lake ice after precip (kg/m^2)
!@var TRUNOSI tracer in run off from sea/lake ice after surface (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: oTRUNOSI, oTRUNPSI
!@var TRMELTI tracer from simelt into ocean (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oTRMELTI
!@var TRPREC tracers in precip (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: oTRPREC

!@var TRGMELT tracer from glacial melt into ocean (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oTRGMELT
#endif
!@var DTRSI tracer flux in sea ice under ice and on open water (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: oDTRSI

#ifdef TRACERS_DRYDEP
!@var TRDRYDEP tracer dry deposition by type (kg/m^2) (positive down)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: oTRDRYDEP 
#endif
#endif

      END MODULE OFLUXES

      SUBROUTINE ALLOC_OFLUXES(grd_dum)
!@sum   Initializes FLUXES''s arrays
!@auth  Larissa Nazarenko
!@ver  1.0
      USE CONSTANT, only : tf
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      USE OFLUXES
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO

      !I-J arrays
      ALLOCATE( oRUNOSI  ( I_0H:I_1H , J_0H:J_1H ), 
     &          oERUNOSI ( I_0H:I_1H , J_0H:J_1H ), 
     &          oSRUNOSI ( I_0H:I_1H , J_0H:J_1H ),
     &          oRUNPSI  ( I_0H:I_1H , J_0H:J_1H ), 
     &          oSRUNPSI ( I_0H:I_1H , J_0H:J_1H ),
     &          oERUNPSI ( I_0H:I_1H , J_0H:J_1H ),
     &          oDMUI    ( I_0H:I_1H , J_0H:J_1H ),
     &          oDMVI    ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT=IER )
      ALLOCATE( oFLOWO   ( I_0H:I_1H , J_0H:J_1H ),
     &          oEFLOWO  ( I_0H:I_1H , J_0H:J_1H ),
     &          oMELTI   ( I_0H:I_1H , J_0H:J_1H ),
     &          oEMELTI  ( I_0H:I_1H , J_0H:J_1H ),
     &          oSMELTI  ( I_0H:I_1H , J_0H:J_1H ),
     &          oGMELT   ( I_0H:I_1H , J_0H:J_1H ),
     &          oEGMELT  ( I_0H:I_1H , J_0H:J_1H ),
     &          oPREC    ( I_0H:I_1H , J_0H:J_1H ),
     &          oEPREC   ( I_0H:I_1H , J_0H:J_1H ),
     &          oAPRESS  ( I_0H:I_1H , J_0H:J_1H ),
     &          oRSI     ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT=IER)


       !n-I-J arrays
       ALLOCATE( oDMSI    (  2  , I_0H:I_1H , J_0H:J_1H ), 
     &           oDHSI    (  2  , I_0H:I_1H , J_0H:J_1H ), 
     &           oDSSI    (  2  , I_0H:I_1H , J_0H:J_1H ),
     &           oSOLAR   (  3  , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)


      !I-J-: arrays
      ALLOCATE( oE0      ( I_0H:I_1H , J_0H:J_1H , 1 ),
     &          oEVAPOR  ( I_0H:I_1H , J_0H:J_1H , 1 ),
     &          oDMUA    ( I_0H:I_1H , J_0H:J_1H , 1 ),
     &          oDMVA    ( I_0H:I_1H , J_0H:J_1H , 1 ),
     &   STAT = IER)


!TRACERS **********************
#ifdef TRACERS_OCEAN

      !:-:-I-J arrays
#ifdef TRACERS_GASEXCH_ocean

      ALLOCATE( oTRGASEX(ntm_gasexch, 1 , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
      oTRGASEX=0.     !initialize to zero

#endif

#ifdef TRACERS_WATER
      !(:)-(:)-I-J arrays
      ALLOCATE( oTREVAPOR( NTM , 1, I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)


       !:-I-J arrays
       ALLOCATE( oTRPREC  ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           oTRUNPSI ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           oTRUNOSI ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           oTRFLOWO ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           oTRMELTI ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)

       ALLOCATE( oTRGMELT ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
#endif

#ifdef TRACERS_DRYDEP
       ALLOCATE(oTRDRYDEP( NTM , 1 , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
         oTRDRYDEP = 0.   !Initialize to 0.
#endif

      ALLOCATE( oDTRSI( NTM ,    2   , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
#endif


      END SUBROUTINE ALLOC_OFLUXES

      module hntrp_mod
!@sum hntrp_mod contains domain-decomposed procedures for
!@+   conservative regridding of lat-lon fields.
!@+   Pole rotations are not supported.
!@auth G. Russell   original code
!@+    M. Kelley    F90+parallel repackaging 08/2009

c import domain decomposition structures/procedures for parallelization
      Use domain_decomp_1d, only :
     &     dist_grid,
     &     band_pack_type,init_band_pack_type,
     &     band_pack,band_pack_column

      implicit none

c common /hntrcb/ was converted to a derived type to
c facilitate keeping several instances of regrid
c and communication info.
      type hntrp_type
        Real*8 :: SINA(0:361),SINB(0:361),
     *     FMIN(720),FMAX(720),GMIN(361),GMAX(361), DATMIS
        Integer :: IMIN(720),IMAX(720),JMIN(361),JMAX(361),
     *       IMA,JMA, IMB,JMB,J1B,JNB
        Integer :: J1B_halo,JNB_halo
        type(band_pack_type) :: bpack ! communication info
      end type hntrp_type

      contains

      Subroutine Init_Hntrp_Type(htype,
     *     grid_A,OFFIA,DLATA,
     *     grid_B,OFFIB,DLATB,
     *     DATMIS,
! set these optional args to JMA-1 or JMB-1 for secondary lats
     *     JMA_4interp,JMB_4interp)
C****
C**** HNTR80 fills in the common block HNTRCB with coordinate
C**** parameters that will be used by subsequent calls to HNTR8.
C**** The 5 Real input values are expected to be Real*8.
C****
C**** Input: IMA = number of east-west cell in global grid A
C****        JMA = number of north-south cell in global grid A
C****        J1A = first latitude cell of grid A for present processor
C****        JNA = last latitude cell of grid A for present processor
C****      OFFIA = number of cells of grid A in east-west direction
C****              from IDL (180) to western edge of cell IA=1
C****      DLATA = minutes of latitude for non-polar cells on grid A
C****        IMB = number of cells in east-west direction of grid B
C****        JMB = number of cells in north-south direction of grid B
C****        J1B = first latitude cell of grid B for present processor
C****        JNB = last latitude cell of grid B for present processor
C****      OFFIB = number of cells of grid B in east-west direction
C****              from IDL (180) to western edge of cell IB=1
C****      DLATB = minutes of latitude for non-polar cells on grid B
C****     DATMIS = missing data value inserted in output array B when
C****              cell (IB,JB) has integrated value 0 of WTA
C****
C**** Output: common block /HNTRCB/
C**** SINA(JA) = sine of latitude of northern edge of cell JA on grid A
C**** SINB(JB) = sine of latitude of northern edge of cell JB on grid B
C**** FMIN(IB) = fraction of cell IMIN(IB) on grid A west of cell IB
C**** FMAX(IB) = fraction of cell IMAX(IB) on grid A east of cell IB
C**** GMIN(JB) = fraction of cell JMIN(JB) on grid A south of cell JB
C**** GMAX(JB) = fraction of cell JMAX(JB) on grid A north of cell JB
C**** IMIN(IB) = western most cell of grid A that intersects cell IB
C**** IMAX(IB) = eastern most cell of grid A that intersects cell IB
C**** JMIN(JB) = southern most cell of grid A that intersects cell JB
C**** JMAX(JB) = northern most cell of grid A that intersects cell JB
C****
      Implicit None
      Real*8,Parameter :: TWOPI = 6.283185307179586477d0
      Real*8 OFFIA,DLATA, OFFIB,DLATB, DATMIS
      type(hntrp_type) :: htype
      type(dist_grid) :: grid_A, grid_B
      integer, optional, intent(in) :: JMA_4interp,JMB_4interp
C**** Local vars
      Real*8 :: DIA,DIB,RIA,RIB,RJA,RJB,FJEQA,FJEQB
      Integer :: IMA,JMA,J1A,JNA, IMB,JMB,J1B,JNB
      Integer :: IA,IB,JA,JB,IBp1
      Integer :: jmin_pack,jmax_pack
C****
      IMA = grid_A%IM_WORLD
      JMA = grid_A%JM_WORLD
      if(present(JMA_4interp)) JMA=JMA_4interp
      J1A = 1
      JNA = JMA
      IMB = grid_B%IM_WORLD
      JMB = grid_B%JM_WORLD
      if(present(JMB_4interp)) JMB=JMB_4interp
      J1B = grid_B%J_STRT
      JNB = grid_B%J_STOP
      htype%IMA = IMA
      htype%JMA = JMA
      htype%IMB = IMB
      htype%JMB = JMB
      htype%J1B = J1B
      htype%JNB = JNB
      htype%J1B_halo = grid_B%J_STRT_halo
      htype%JNB_halo = grid_B%J_STOP_halo
      htype%DATMIS = DATMIS
      If (IMA<1 .or. IMA>720 .or. JMA<1 .or. JMA>361 .or.
     *    IMB<1 .or. IMB>720 .or. JMB<1 .or. JMB>361)  GoTo 800
C****
C**** Partitions in east-west (I) direction
C**** Domain, around the globe, is scaled to fit from 0 to IMA*IMB
C****
      DIA = IMB  !  width of single A grid cell in scaled domain
      DIB = IMA  !  width of single B grid cell in scaled domain
      IA  = 1
      RIA = (IA+OFFIA - IMA)*IMB  !  scaled longitude of eastern edge
      IB  = IMB
      Do 150 IBp1=1,IMB
      RIB = (IBp1-1+OFFIB)*IMA    !  scaled longitude of eastern edge
  110 If (RIA-RIB)  120,130,140
  120 IA  = IA  + 1
      RIA = RIA + DIA
      GoTo 110
C**** Eastern edges of cells IA of grid A and IB of grid B coincide
  130 htype%IMAX(IB) = IA
      htype%FMAX(IB) = 0
      IA  = IA  + 1
      RIA = RIA + DIA
      htype%IMIN(IBp1) = IA
      htype%FMIN(IBp1) = 0
      GoTo 150
C**** Cell IA of grid A contains western edge of cell IB of grid B
  140 htype%IMAX(IB) = IA
      htype%FMAX(IB) = (RIA-RIB)/DIA
      htype%IMIN(IBp1) = IA
      htype%FMIN(IBp1) = 1-htype%FMAX(IB)
  150 IB = IBp1
      htype%IMAX(IMB) = htype%IMAX(IMB) + IMA
C       WRITE (0,915) 'IMIN=',htype%IMIN(1:IMB)
C       WRITE (0,915) 'IMAX=',htype%IMAX(1:IMB)
C       WRITE (0,916) 'FMIN=',htype%FMIN(1:IMB)
C       WRITE (0,916) 'FMAX=',htype%FMAX(1:IMB)
C****
C**** Partitions in the north-south (J) direction
C**** Domain is measured in minutes (1/60-th of a degree)
C****
c init JMIN/JMAX if JMA_4interp or JMB_4interp not the full global JM
      htype%JMIN(:) = 2*JMA; htype%JMAX(:) = -1
      FJEQA = .5*(1+JMA)
      Do 210 JA=1,JMA-1 !J1A-1,JNA
      RJA = (JA+.5-FJEQA)*DLATA  !  latitude in minutes of northern edge
  210 htype%SINA(JA) = Sin (RJA*TWOPI/(360*60))
      htype%SINA(0)  = -1
      htype%SINA(JMA) = 1
C****
      FJEQB = .5*(1+JMB)
      Do 220 JB=1,JMB-1 !J1B-1,JNB
      RJB = (JB+.5-FJEQB)*DLATB  !  latitude in minutes of northern edge
  220 htype%SINB(JB) = Sin (RJB*TWOPI/(360*60))
      htype%SINB(0)  = -1
      htype%SINB(JMB) = 1
C****
C**** Calculate JMIN,GMIN,JMAX,GMAX
C****
      htype%JMIN(1) = 1
      htype%GMIN(1) = 0
      JA = 1
      Do 350 JB=1,JMB-1
  310 If (htype%SINA(JA)-htype%SINB(JB))  320,330,340
  320 JA = JA + 1
      GoTo 310
C**** Northern edges of cells JA of grid A and JB of grid B coincide
  330 htype%JMAX(JB) = JA
      htype%GMAX(JB) = 0
      JA = JA + 1
      htype%JMIN(JB+1) = JA
      htype%GMIN(JB+1) = 0
      GoTo 350
C**** Cell JA of grid A contains northern edge of cell JB of grid B
  340 htype%JMAX(JB) = JA
      htype%GMAX(JB) = htype%SINA(JA) - htype%SINB(JB)
      htype%JMIN(JB+1) = JA
      htype%GMIN(JB+1) = htype%SINB(JB) - htype%SINA(JA-1)
  350 Continue
      htype%JMAX(JMB) = JMA
      htype%GMAX(JMB) = 0

c      If (htype%SINB(J1B-1) < htype%SINA(J1A-1))  GoTo 830
c      If (htype%SINB(JNB)   > htype%SINA(JNA)  )  GoTo 840
c      JA = J1A-1
c      JB = J1B-1
c  300 If (htype%SINA(JA)-htype%SINB(JB))  310,320,330
c  310 JA = JA+1
c      GoTo 300
cC**** Northern edges of cells JA of grid A and JB of grid B coincide
c  320 JA = JA+1
c      JB = JB+1
c      htype%JMIN(JB) = JA
c      htype%GMIN(JB) = 0
c      GoTo 400
cC**** Cell JA of grid A contains northern edge of cell JB of grid B
c  330 JB = JB+1
c      htype%JMIN(JB) = JA
c      htype%GMIN(JB) = htype%SINB(JB-1) - htype%SINA(JA-1)
cC****
cC**** Calculate JMIN(J1B+1:JNB), GMIN(J1B+1:JNB),
cC****       and JMAX(J1B:JNB-1), GMAX(J1B:JNB-1)
cC****
c  400 If (htype%SINA(JA)-htype%SINB(JB))  410,420,430
c  410 JA = JA+1
c      GoTo 400
cC**** Northern edges of cells JA of grid A and JB of grid B coincide
c  420 htype%JMAX(JB) = JA
c      htype%GMAX(JB) = 0
c      If (JB == JNB)  GoTo 500
c      JA = JA+1
c      JB = JB+1
c      htype%JMIN(JB) = JA
c      htype%GMIN(JB) = 0
c      GoTo 400
cC**** Cell JA of grid A contains northern edge of cell JB of grid B
c  430 htype%JMAX(JB) = JA
c      htype%GMAX(JB) = htype%SINA(JA) - htype%SINB(JB)
c      If (JB == JNB)  GoTo 500
c      JB = JB+1
c      htype%JMIN(JB) = JA
c      htype%GMIN(JB) = htype%SINB(JB-1) - htype%SINA(JA-1)
c      GoTo 400
c  500 Continue
C 500   WRITE (0,951) 'JMIN=',htype%JMIN(1:JMB)
C       WRITE (0,951) 'JMAX=',htype%JMAX(1:JMB)
C       WRITE (0,952) 'GMIN=',htype%GMIN(1:JMB)
C       WRITE (0,952) 'GMAX=',htype%GMAX(1:JMB)
C 951   FORMAT (/ 1X,A5 / (20I6))
C 952   FORMAT (/ 1X,A5 / (20(1X,F5.4)))

      if(grid_B%have_domain) then
        jmin_pack = minval(htype%jmin(j1b:jnb))
        jmax_pack = maxval(htype%jmax(j1b:jnb))
      else ! request nonexistent data which will not be sent
        jmin_pack = grid_B%JM_WORLD+1
        jmax_pack = grid_B%JM_WORLD+1
      endif
      call init_band_pack_type(grid_A, grid_B, jmin_pack,jmax_pack,
     &     htype%bpack)

      RETURN
C****
C**** Invalid parameters or dimensions out of range
C****
  800 Write (0,980) IMA,JMA,J1A,JNA,OFFIA,DLATA,
     *              IMB,JMB,J1B,JNB,OFFIB,DLATB, DATMIS
      Stop 800
  830 Write (0,*)
     *'Southern edge of grid B is south of southern edge of grid A.'
      GoTo 800
  840 Write (0,*)
     *'Northern edge of grid B is north of northern edge of grid A.'
      GoTo 800
C****
  980 Format (/ ' Arguments received by HNTRP0 in order:'/
     *   2I12,' = IMA,JMA = array dimensions for A grid'/
     *   2I12,' = J1A,JNA = local latitude band for A grid'/
     *  E24.8,' = OFFIA   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATA   = minutes of latitude for interior grid cell'/
     *   2I12,' = IMB,JMB = array dimensions for B grid'/
     *   2I12,' = J1B,JNB = local latitude band for B grid'/
     *  E24.8,' = OFFIB   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATB   = minute of latitude for interior grid cell'/
     *  E24.8,' = DATMIS  = missing data value to be put in B array',
     *                    ' when integrated WTA = 0'/
     *  ' These arguments are invalid or out of range.' /
     *  ' IMA and IMB may not exceed 720.' /
     *  ' JMA and JMB may not exceed 360.')
      End Subroutine Init_Hntrp_Type

      Subroutine HNTR8_band (WTA,A,htype,B)
C****
C**** HNTR8 performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.  B grid values that cannot be calculated because the
C**** covering A grid boxes have WTA = 0, are set to the value DATMIS.
C**** The area weighted integral of the quantity is conserved.
C**** The 3 Real input values are expected to be Real*8.
C****
C**** Input: WTA = weighting array for values on the A grid
C****          A = per unit area or per unit mass quantity
C**** Output:  B = horizontally interpolated quantity on B grid
C****
C**** The algorithm boils down to:
C****   dAREA(J) = SIN(J) - SIN(J-1)
C****   B(I,J) = Sum[dAREA(J)*WTA(I,J)*A(I,J)] / Sum[dAREA(J)*WTA(I,J)]
C****
      Implicit None
      type(hntrp_type) :: htype
      Real*8, dimension(htype%IMA,
     &     htype%bpack%jband_strt:htype%bpack%jband_stop) :: WTA,A
      Real*8, dimension(htype%IMB,htype%J1B_halo:htype%JNB_halo) :: B
C**** Local vars
      Integer :: IA,JA,IAREV,IB,JB,IAMIN,IAMAX,JAMIN,JAMAX
      Real*8 :: WEIGHT,VALUE,F,G
      if(htype%J1B > htype%JNB) return
C****
C**** Interpolate from grid A to grid B
C****
      Do JB=htype%J1B,htype%JNB
        JAMIN = htype%JMIN(JB)
        JAMAX = htype%JMAX(JB)
        Do IB=1,htype%IMB
          WEIGHT= 0
          VALUE = 0
          IAMIN = htype%IMIN(IB)
          IAMAX = htype%IMAX(IB)
          Do JA=JAMIN,JAMAX
            G = htype%SINA(JA)-htype%SINA(JA-1)
            If (JA==JAMIN)  G = G - htype%GMIN(JB)
            If (JA==JAMAX)  G = G - htype%GMAX(JB)
            Do IAREV=IAMIN,IAMAX
              IA  = 1 + Mod(IAREV-1,htype%IMA)
              F   = 1
              If (IAREV==IAMIN)  F = F - htype%FMIN(IB)
              If (IAREV==IAMAX)  F = F - htype%FMAX(IB)
              WEIGHT = WEIGHT + F*G*WTA(IA,JA)
              VALUE  = VALUE  + F*G*WTA(IA,JA)*A(IA,JA)
            End Do
          End Do
          B(IB,JB) = htype%DATMIS
          If (WEIGHT /= 0)  B(IB,JB) = VALUE/WEIGHT
        End Do
      End Do
      Return
      End Subroutine HNTR8_band

      Subroutine HNTR8_band_lij (A,htype,B)
C****
C**** HNTR8 performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.
C**** The area weighted integral of the quantity is conserved.
C****
C**** Input:   A = per unit area or per unit mass quantity
C**** Output:  B = horizontally interpolated quantity on B grid
C****
C**** The algorithm boils down to:
C****   dAREA(J) = SIN(J) - SIN(J-1)
C****   B(I,J) = Sum[dAREA(J)*A(I,J)] / Sum[dAREA(J)]
C****
      Implicit None
      type(hntrp_type) :: htype
      Real*8, dimension(:,:,htype%bpack%jband_strt:) :: A
      Real*8, dimension(:,:,htype%J1B_halo:) :: B
C**** Local vars
      Integer :: IA,JA,IAREV,IB,JB,IAMIN,IAMAX,JAMIN,JAMAX
      Real*8 :: WEIGHT,F,G
      if(htype%J1B > htype%JNB) return
C****
C**** Interpolate from grid A to grid B
C****
      Do JB=htype%J1B,htype%JNB
        JAMIN = htype%JMIN(JB)
        JAMAX = htype%JMAX(JB)
        Do IB=1,htype%IMB
          WEIGHT= 0
          B(:,IB,JB) = 0
          IAMIN = htype%IMIN(IB)
          IAMAX = htype%IMAX(IB)
          Do JA=JAMIN,JAMAX
            G = htype%SINA(JA)-htype%SINA(JA-1)
            If (JA==JAMIN)  G = G - htype%GMIN(JB)
            If (JA==JAMAX)  G = G - htype%GMAX(JB)
            Do IAREV=IAMIN,IAMAX
              IA  = 1 + Mod(IAREV-1,htype%IMA)
              F   = 1
              If (IAREV==IAMIN)  F = F - htype%FMIN(IB)
              If (IAREV==IAMAX)  F = F - htype%FMAX(IB)
              WEIGHT = WEIGHT + F*G
              B(:,IB,JB) = B(:,IB,JB) + F*G*A(:,IA,JA)
            End Do
          End Do
          If(WEIGHT /= 0) B(:,IB,JB) = B(:,IB,JB)/WEIGHT
        End Do
      End Do
      Return
      End Subroutine HNTR8_band_lij

      Subroutine HNTR8P_band (WTA,A,htype,B)
C****
C**** HNTR8P is similar to HNTR8 but polar values are replaced by
C**** their longitudinal mean.
C**** The 3 Real input values are expected to be Real*8.
C****
      Implicit None
      type(hntrp_type) :: htype
      Real*8, dimension(htype%IMA,
     &     htype%bpack%jband_strt:htype%bpack%jband_stop) :: WTA,A
      Real*8, dimension(htype%IMB,htype%J1B_halo:htype%JNB_halo) :: B
C**** Local vars
      Integer :: I, IMB,JMB,J1B,JNB
      Real*8 :: WEIGHT,VALUE,BMEAN,DATMIS
C****
      if(htype%J1B > htype%JNB) return
      Call HNTR8_band (WTA,A,htype,B)
C****
C**** Replace individual values at a pole by longitudinal mean
C****
      IMB = htype%IMB
      JMB = htype%JMB
      J1B = htype%J1B
      JNB = htype%JNB
      DATMIS = htype%DATMIS
C**** South pole
      If (J1B == 1) then
        WEIGHT = 0
        VALUE  = 0
        Do I=1,IMB
          If (DATMIS /= 0. .and. B(I,1) == DATMIS) cycle
          WEIGHT = WEIGHT + 1
          VALUE  = VALUE  + B(I,1)
        End Do
        BMEAN = DATMIS
        If (WEIGHT /= 0)  BMEAN = VALUE/WEIGHT
        B(1:IMB,1) = BMEAN
      End If
C**** North pole
      If (JNB == JMB) then
        WEIGHT = 0
        VALUE  = 0
        Do I=1,IMB
          If (DATMIS /= 0. .and. B(I,JMB) == DATMIS) cycle
          WEIGHT = WEIGHT + 1
          VALUE  = VALUE  + B(I,JMB)
        End Do
        BMEAN  = DATMIS
        If (WEIGHT.ne.0)  BMEAN = VALUE/WEIGHT
        B(1:IMB,JMB) = BMEAN
      End If
      Return
      End Subroutine HNTR8P_band

      end module hntrp_mod
