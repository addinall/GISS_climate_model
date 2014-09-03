      MODULE DYNAMICS
!@sum  DYNAMICS contains all the pressure and momentum related variables
!@auth Original development team
!@ver  1.0
      USE DOMAIN_DECOMP_ATM, ONLY : grid
      USE RESOLUTION , ONLY : LM
      IMPLICIT NONE
      SAVE
C**** Some helpful arrays (arrays should be L first)
!@var  PLIJ  Surface pressure: P(I,J) or PSF-PTOP (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PLIJ
!@var  PDSIG  Surface pressure * DSIG(L) (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PDSIG
!@var  AM  Air mass of each box (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: AM     ! PLIJ*DSIG(L)*100/grav
!@var  BYAM  1/Air mass (m^2/kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: BYAM
!@var  PMID  Pressure at mid point of box (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PMID    ! SIG(L)*PLIJ+PTOP
!@var  PK   PMID**KAPA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PK
!@var  PEUP  Pressure at lower edge of box (incl. surface) (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PEDN  ! SIGE(L)*PLIJ+PTOP
!@var  PEK  PEUP**KAPA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PEK
!@var  SQRTP  square root of P (used in diagnostics)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SQRTP
!@var  PTROPO  Pressure at mid point of tropopause level (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PTROPO
!@var  LTROPO  Tropopause layer
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LTROPO

C**** module should own dynam variables used by other routines
!@var PTOLD pressure at beginning of dynamic time step (for clouds)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)    :: PTOLD
!@var SD_CLOUDS vert. integrated horizontal convergence (for clouds)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SD_CLOUDS
!@var GZ geopotential height (for Clouds and Diagnostics)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GZ
!@var DPDX_BY_RHO,DPDY_BY_RHO (pressure gradients)/density at L=1
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDX_BY_RHO,DPDY_BY_RHO
!@var DPDX_BY_RHO_0,DPDY_BY_RHO_0 surface (pressure gradients)/density
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDX_BY_RHO_0
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDY_BY_RHO_0

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PU,PV
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: CONV
cgsfc      REAL*8, DIMENSION(IM,JM,LM-1) :: SD
      REAL*8, ALLOCATABLE :: SD(:,:,:)
!@var PIT  pressure tendency (mb m^2/s)
      REAL*8, ALLOCATABLE :: PIT(:,:)
cgsfc      EQUIVALENCE (SD(1,1,1),CONV(1,1,2))

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PHI,SPA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DUT,DVT
!!@var xAVRX scheme-depend. coefficient for AVRX: 1,byrt2 (2nd,4th order)
!      REAL*8 xAVRX

!@var PUA,PVA,SDA,PS save PU,PV,SD,P for hourly tracer advection
!@var MB Air mass array for tracers (before advection)
!@var MA Air mass array for tracers (updated during advection)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PUA,PVA,SDA,MB,MA
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PS

!@var WCP vertical mass flux in a constant-pressure vertical
!@+   coordinate whose pressure levels are the global means of
!@+   each layer.
!@var WCPsig: WCP interpolated to terrain-following coordinate
!@+   surfaces (sigma levels)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WCP,WCPsig

!@var DKE change in KE due to dissipation (SURF/DC/MC) (m^2/s^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DKE
!@var KEA KE on the A grid (m^2/s^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: KEA ! ke on A grid
!@var UALIJ,VALIJ U,V on the A grid (m/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: UALIJ,VALIJ
!@var WSAVE vertical velocity (m/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WSAVE
!@var SMASS local but "SAVE"d array ADVECV in MOMEN2ND made global
!@    here since its use does not go beyond ATMDYN that calls ADVECV
      REAL*8, ALLOCATABLE:: SMASS(:)

      INTEGER, DIMENSION(2) :: t_dyn_a, t_dyn_b, t_dyn_c, t_dyn_d
      INTEGER, DIMENSION(2) :: t_aflux, t_advecm, t_advecv
      INTEGER, DIMENSION(2) :: t_PGF, t_filter, t_fltruv, t_calc_pijl
      INTEGER, DIMENSION(2) :: t_sdrag, t_calc_trop

      REAL*8, PARAMETER :: COS_LIMIT = 0.15d0

      END MODULE DYNAMICS


      SUBROUTINE ALLOC_DYNAMICS(grid)
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
      USE RESOLUTION , ONLY : LM
      USE DYNAMICS, ONLY : PLIJ,PDSIG,AM,BYAM,PMID,PK,PEDN,PEK,
     $                     SD_CLOUDS,GZ,PU,PV,CONV,PHI,SPA,DUT,
     $                     DVT,PUA,PVA,SDA,MB,MA,DKE,WSAVE,SD,PIT,
     $                     SQRTP,PTROPO,LTROPO,PTOLD,DPDX_BY_RHO,
     $                     DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0,
     $                     PS,SMASS,KEA,UALIJ,VALIJ,WCP,WCPsig
      USE DYNAMICS, only : t_dyn_a, t_dyn_b, t_dyn_c, t_dyn_d
      USE DYNAMICS, only : t_aflux, t_advecm, t_advecv

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO
      
      ! K-I-J arrays
      ALLOCATE ( PLIJ(LM,I_0H:I_1H,J_0H:J_1H), 
     $           PDSIG(LM,I_0H:I_1H,J_0H:J_1H),
     $             AM(LM,I_0H:I_1H,J_0H:J_1H),  
     $           BYAM(LM,I_0H:I_1H,J_0H:J_1H),
     $           PMID(LM,I_0H:I_1H,J_0H:J_1H),    
     $           PK(LM,I_0H:I_1H,J_0H:J_1H),
     $           PEDN(LM+1,I_0H:I_1H,J_0H:J_1H), 
     $            PEK(LM+1,I_0H:I_1H,J_0H:J_1H),
     $   STAT = IER)
      ! I-J-K arrays
      ALLOCATE( SD_CLOUDS(I_0H:I_1H,J_0H:J_1H,LM),  
     $                 GZ(I_0H:I_1H,J_0H:J_1H,LM), 
     $                 PU(I_0H:I_1H,J_0H:J_1H,LM),  
     $                 PV(I_0H:I_1H,J_0H:J_1H,LM), 
     $               CONV(I_0H:I_1H,J_0H:J_1H,LM), 
     $               PHI(I_0H:I_1H,J_0H:J_1H,LM), 
     $                SPA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DUT(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DVT(I_0H:I_1H,J_0H:J_1H,LM), 
     $                PUA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                PVA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                SDA(I_0H:I_1H,J_0H:J_1H,LM),  
     $                MB(I_0H:I_1H,J_0H:J_1H,LM), 
     $                MA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DKE(I_0H:I_1H,J_0H:J_1H,LM), 
     $                KEA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                UALIJ(LM,I_0H:I_1H,J_0H:J_1H), 
     $                VALIJ(LM,I_0H:I_1H,J_0H:J_1H), 
     $              WSAVE(I_0H:I_1H,J_0H:J_1H,LM-1), 
     $              PIT(I_0H:I_1H,J_0H:J_1H),
     $              SD(I_0H:I_1H,J_0H:J_1H,LM-1),
     $   STAT = IER)

      !hack to remove NaNs from "SD" - possible mistake somewhere
      CONV(:,:,:) = 0.d0

      ! I-J arrays
      ALLOCATE(  SQRTP(I_0H:I_1H,J_0H:J_1H), 
     $          PTROPO(I_0H:I_1H,J_0H:J_1H),
     $          LTROPO(I_0H:I_1H,J_0H:J_1H),  
     $          PTOLD(I_0H:I_1H,J_0H:J_1H),
     $          DPDX_BY_RHO(I_0H:I_1H,J_0H:J_1H), 
     $          DPDY_BY_RHO(I_0H:I_1H,J_0H:J_1H),
     $          DPDX_BY_RHO_0(I_0H:I_1H,J_0H:J_1H), 
     $          DPDY_BY_RHO_0(I_0H:I_1H,J_0H:J_1H),
     $          PS(I_0H:I_1H,J_0H:J_1H),
     $          WCP(I_0H:I_1H,J_0H:J_1H,LM),
     $          WCPsig(I_0H:I_1H,J_0H:J_1H,LM),
     $   STAT = IER)

      ! J arrays
      ALLOCATE(  SMASS(J_0H:J_1H), 
     $   STAT = IER)

! correct or wrong, but being static all arrays were initialized
! to zero by default. They have to be initialized to something now 
! to avoid floating point exceptions...
      DPDX_BY_RHO(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDY_BY_RHO(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDX_BY_RHO_0(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDY_BY_RHO_0(I_0H:I_1H,J_0H:J_1H) = 0.d0

      PU(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0
      PV(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0

      SD_CLOUDS(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0
      PIT(:,:) = 0.d0

      END SUBROUTINE ALLOC_DYNAMICS


