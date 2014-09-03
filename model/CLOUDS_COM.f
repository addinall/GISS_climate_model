#include "rundeck_opts.h"
      MODULE CLOUDS_COM
!@sum  CLOUDS_COM model variables for moist convction and
!@+    large-scale condensation
!@+    cloud droplet number added to list of saves for rsf files
!@auth M.S.Yao/T. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
      USE MODEL_COM, only : IM,JM,LM
      IMPLICIT NONE
      SAVE
!@var TTOLD,QTOLD previous potential temperature, humidity
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TTOLD,QTOLD
!@var SVLHX,SVLAT previous latent heat of evaporation
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SVLHX,SVLAT
!@var RHSAV previous relative humidity
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: RHSAV
C**** some arrays here for compatility with new clouds
!@var CLDSAV, CLDSAV1 previous cloud cover area (percent)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CLDSAV,CLDSAV1
!@var ULS,VLS,UMC,VMC velocity work arrays
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ULS,VLS,UMC,VMC
!@var TLS,QLS,TMC,QMC temperature and humidity work arrays
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TLS,QLS,TMC,QMC
!@var FSS grid fraction for large-scale clouds
!@+   initialised as 1. for compatibility with previous clouds
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: FSS
#ifdef CLD_AER_CDNC
!@var OLDNL old CDNC,OLDNI old ice crystal
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: OLDNL,OLDNI
!@var N, Re, LWP for 3 hrly diag save
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CDN3D,CRE3D
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: CLWP
#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
!@var LWC,Cld depth, cld tem for 3 hrly diag save
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CL3D,CI3D,CD3D,CTEM
#endif

C**** variables saved for radiation calculations
!@var TAUSS optical depth from super-saturated clouds
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAUSS
!@var TAUMC optical depth from moist-convective clouds
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAUMC
!@var CLDSS super-saturated cloud cover area (percent)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CLDSS
!@var CLDMC moist convective cloud cover area (percent)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CLDMC
!@var CSIZMC,CSIZSS mc,ss effective cloud droplet radius (microns)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CSIZMC,CSIZSS

C**** variables saved for surface wind spectrum calculations
!@var DDM1 downdraft mass flux / rho at lowest level (m/s)
!@var DDML lowest level of downdraft 
!@var DDMS downdraft mass flux at level 1 (kg/s/m**2)
!@var TDN1 downdraft temperature (K)
!@var QDN1 downdraft humidity (kg/kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DDM1,DDMS,TDN1,QDN1
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: DDML

C**** variables used (and saved) for gravity wave drag calculations
!@var AIRX, AIRMX*AREA convective mass flux (kg/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: AIRX
!@var LMC max layer of mc convective mass flux.
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: LMC

!@var LLOW,LMID,LHI max levels for low, mid and high clouds
      INTEGER LLOW,LMID,LHI

!@var ISCCP_REG2D distributed version of ISCCP_REG
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISCCP_REG2D

!@var UKM,VKM arrays for vertical momentum mixing
      REAL*8, dimension(:,:,:,:), allocatable :: UKM,VKM
#ifdef BLK_2MOM
#ifdef TRACERS_AMP
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: NACTC      ! = 1.0D-30  ![#/m3](l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: NAERC      ! = 1.0D-30  ![#/m3](l,nmodes)
#endif
#endif

C**** ISCCP diagnostics related parameter
#ifdef SCM
      INTEGER, parameter, public :: ncol = 100    !@var ncol number of subcolumns
#else
      INTEGER, parameter, public :: ncol = 20    !@var ncol number of subcolumns
#endif
#ifdef CFMIP
!@dbparam Nsubdd_for_cfmip: sub-daily diag freq for CFMIP/COSP
      INTEGER :: Nsubdd_for_cfmip
#endif
      END MODULE CLOUDS_COM

      SUBROUTINE ALLOC_CLOUDS_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
      USE MODEL_COM, ONLY : IM,LM
#ifdef BLK_2MOM
#ifdef TRACERS_AMP
      USE AERO_CONFIG, ONLY: NMODES
#endif
#endif
      USE CLOUDS_COM, ONLY : TTOLD,QTOLD,SVLHX,SVLAT,RHSAV,CLDSAV,
     *                       CLDSAV1,FSS,
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
     *                       CL3D,CI3D,CD3D,CTEM,
#endif
#ifdef CLD_AER_CDNC
     *                       OLDNL,OLDNI,
     *                       CDN3D,CRE3D,CLWP,
#endif
     *                       TAUSS,TAUMC, CLDSS,CLDMC,CSIZMC,CSIZSS,
     *                       ULS,VLS,UMC,VMC,TLS,QLS,
     *                       TMC,QMC,DDM1,AIRX,LMC,DDMS,TDN1,QDN1,DDML
#ifdef BLK_2MOM
#ifdef TRACERS_AMP
     *                       ,NACTC,NAERC
#endif
#endif
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE(     TTOLD(LM,I_0H:I_1H,J_0H:J_1H),
     *              QTOLD(LM,I_0H:I_1H,J_0H:J_1H),
     *              SVLHX(LM,I_0H:I_1H,J_0H:J_1H),
     *              SVLAT(LM,I_0H:I_1H,J_0H:J_1H),
     *              RHSAV(LM,I_0H:I_1H,J_0H:J_1H),
     *             CLDSAV(LM,I_0H:I_1H,J_0H:J_1H),
     *            CLDSAV1(LM,I_0H:I_1H,J_0H:J_1H),
     *                FSS(LM,I_0H:I_1H,J_0H:J_1H),
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
     *             CTEM(LM,I_0H:I_1H,J_0H:J_1H),
     *             CD3D(LM,I_0H:I_1H,J_0H:J_1H),
     *             CL3D(LM,I_0H:I_1H,J_0H:J_1H),
     *             CI3D(LM,I_0H:I_1H,J_0H:J_1H),
#endif
#ifdef CLD_AER_CDNC
     *             OLDNL(LM,I_0H:I_1H,J_0H:J_1H),
     *             OLDNI(LM,I_0H:I_1H,J_0H:J_1H),
     *             CDN3D(LM,I_0H:I_1H,J_0H:J_1H),
     *             CRE3D(LM,I_0H:I_1H,J_0H:J_1H),
     *             CLWP(I_0H:I_1H,J_0H:J_1H),
#endif
     *              TAUSS(LM,I_0H:I_1H,J_0H:J_1H),
     *              TAUMC(LM,I_0H:I_1H,J_0H:J_1H),
     *              CLDSS(LM,I_0H:I_1H,J_0H:J_1H),
     *              CLDMC(LM,I_0H:I_1H,J_0H:J_1H),
     *             CSIZMC(LM,I_0H:I_1H,J_0H:J_1H),
     *             CSIZSS(LM,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(     ULS(I_0H:I_1H,J_0H:J_1H,LM),
     *              VLS(I_0H:I_1H,J_0H:J_1H,LM),
     *              UMC(I_0H:I_1H,J_0H:J_1H,LM),
     *              VMC(I_0H:I_1H,J_0H:J_1H,LM),
     *              TLS(I_0H:I_1H,J_0H:J_1H,LM),
     *              QLS(I_0H:I_1H,J_0H:J_1H,LM),
     *              TMC(I_0H:I_1H,J_0H:J_1H,LM),
     *              QMC(I_0H:I_1H,J_0H:J_1H,LM),
     *         STAT=IER)

!@var FSS initialized to 1.
      FSS = 1.
#ifdef CLD_AER_CDNC
!@var OLDNL is initialized to 10.0 cm-3
!@var OLDNI is initialised to 0.1 l^-1 or 10^-4 cm-3
      OLDNL = 10.
      OLDNI = 1.d-4
#endif

      ALLOCATE(     DDM1(I_0H:I_1H,J_0H:J_1H),
     *              AIRX(I_0H:I_1H,J_0H:J_1H),
     *              DDMS(I_0H:I_1H,J_0H:J_1H),
     *              TDN1(I_0H:I_1H,J_0H:J_1H),
     *              QDN1(I_0H:I_1H,J_0H:J_1H),
     *              DDML(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(     LMC(2,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

C**** Initialise some output used in dynamics
      LMC(:,:,J_0H:J_1H)=0
      AIRX(:,J_0H:J_1H)=0.

#ifdef BLK_2MOM
#ifdef TRACERS_AMP
      allocate(  NACTC(LM,nmodes) )
      NACTC      = 1.0D-30
      allocate(  NAERC(LM,nmodes) )
      NAERC      = 1.0D-30
#endif
#endif
      END SUBROUTINE ALLOC_CLOUDS_COM

      SUBROUTINE io_clouds(kunit,iaction,ioerr)
!@sum  io_clouds reads and writes cloud arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT, GRID
      USE DOMAIN_DECOMP_1D, only : PACK_COLUMN, UNPACK_COLUMN
      USE CLOUDS_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "CLD01"
      REAL*8, ALLOCATABLE,  DIMENSION(:,:,:) :: TTOLD_glob,QTOLD_glob
     &                              ,SVLHX_glob,RHSAV_glob,CLDSAV_glob
#ifdef CLD_AER_CDNC
     &     ,OLDNL_glob,OLDNI_glob
#endif

      call allocate_me

      write(MODULE_HEADER(lhead+1:80),'(a)')
     &'R8 dim(im,jm,lm):potT,Hum,LatHeat,RHum,CldCv,NO,NL,SM (all old)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output to standard restart file
        CALL PACK_COLUMN(grid, TTOLD,  TTOLD_glob)
        CALL PACK_COLUMN(grid, QTOLD,  QTOLD_glob)
        CALL PACK_COLUMN(grid, SVLHX,  SVLHX_glob)
        CALL PACK_COLUMN(grid, RHSAV,  RHSAV_glob)
        CALL PACK_COLUMN(grid, CLDSAV, CLDSAV_glob)
#ifdef CLD_AER_CDNC
        CALL PACK_COLUMN(grid, OLDNL, OLDNL_glob)
        CALL PACK_COLUMN(grid, OLDNI, OLDNI_glob)
#endif
        IF (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,
     *       TTOLD_glob,QTOLD_glob,SVLHX_glob,RHSAV_glob,CLDSAV_glob
#ifdef CLD_AER_CDNC
     *       ,OLDNL_glob,OLDNI_glob
#endif
        END IF

      CASE (IOREAD:)            ! input from restart file
        if (AM_I_ROOT()) THEN
           READ (kunit,err=10) HEADER,
     *          TTOLD_glob,QTOLD_glob,SVLHX_glob,RHSAV_glob,CLDSAV_glob
#ifdef CLD_AER_CDNC
     *          ,OLDNL_glob,OLDNI_glob
#endif
          IF (HEADER(1:15).NE.MODULE_HEADER(1:15)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
       END IF
C***ESMF: Unpack global arrays into distributed local arrays.
        CALL UNPACK_COLUMN(grid, TTOLD_glob , TTOLD)
        CALL UNPACK_COLUMN(grid, QTOLD_glob , QTOLD)
        CALL UNPACK_COLUMN(grid, SVLHX_glob , SVLHX)
        CALL UNPACK_COLUMN(grid, RHSAV_glob , RHSAV)
        CALL UNPACK_COLUMN(grid, CLDSAV_glob, CLDSAV)
#ifdef CLD_AER_CDNC
        CALL UNPACK_COLUMN(grid, OLDNL_glob , OLDNL)
        CALL UNPACK_COLUMN(grid, OLDNI_glob , OLDNI)
#endif
      END SELECT

      call deallocate_me
      RETURN
 10   IOERR=1
      call deallocate_me
      RETURN

      contains
      subroutine allocate_me
      integer :: img, jmg, lmg
      if (AM_I_ROOT()) then
         img = IM
         jmg = JM
         lmg = LM
      else
         img = 1
         jmg = 1
         lmg = 1
      end if         
      ALLOCATE( TTOLD_glob(lmg,img,jmg),
     &     QTOLD_glob(lmg,img,jmg),
     &     SVLHX_glob(lmg,img,jmg),
     &     RHSAV_glob(lmg,img,jmg),
     &     CLDSAV_glob(lmg,img,jmg)  
#ifdef CLD_AER_CDNC
     &     ,OLDNL_glob(lmg,img,jmg)  
     &     ,OLDNI_glob(lmg,img,jmg)  
#endif
     & )
      end subroutine allocate_me
      subroutine deallocate_me
      DEALLOCATE( TTOLD_glob,
     &     QTOLD_glob,
     &     SVLHX_glob,
     &     RHSAV_glob,
     &     CLDSAV_glob 
#ifdef CLD_AER_CDNC
     &     ,OLDNL_glob 
     &     ,OLDNI_glob 
#endif
     & )
      end subroutine deallocate_me
      END SUBROUTINE io_clouds

#ifdef NEW_IO
      subroutine def_rsf_clouds(fid)
!@sum  def_rsf_clouds defines cloud array structure in restart files
!@auth M. Kelley
!@ver  beta
      use clouds_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      character(len=20) :: lijstr
      lijstr='(lm,dist_im,dist_jm)'
      call defvar(grid,fid,ttold,'ttold'//lijstr)
      call defvar(grid,fid,qtold,'qtold'//lijstr)
      call defvar(grid,fid,svlhx,'svlhx'//lijstr)
      call defvar(grid,fid,rhsav,'rhsav'//lijstr)
      call defvar(grid,fid,cldsav,'cldsav'//lijstr)
#ifdef CLD_AER_CDNC
      call defvar(grid,fid,oldnl,'oldnl'//lijstr)
      call defvar(grid,fid,oldni,'oldni'//lijstr)
#endif
      call defvar(grid,fid,airx,'airx(dist_im,dist_jm)')
      call defvar(grid,fid,lmc,'lmc(two,dist_im,dist_jm)')
      return
      end subroutine def_rsf_clouds

      subroutine new_io_clouds(fid,iaction)
!@sum  new_io_clouds read/write cloud arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use clouds_com
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to restart file
        call write_dist_data(grid, fid, 'ttold', ttold, jdim=3)
        call write_dist_data(grid, fid, 'qtold', qtold, jdim=3)
        call write_dist_data(grid, fid, 'svlhx', svlhx, jdim=3)
        call write_dist_data(grid, fid, 'rhsav', rhsav, jdim=3)
        call write_dist_data(grid, fid, 'cldsav', cldsav, jdim=3)
#ifdef CLD_AER_CDNC
        call write_dist_data(grid, fid, 'oldnl', oldnl, jdim=3)
        call write_dist_data(grid, fid, 'oldni', oldni, jdim=3)
#endif
        call write_dist_data(grid, fid, 'airx', airx)
        call write_dist_data(grid, fid, 'lmc', lmc, jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'ttold', ttold, jdim=3)
        call read_dist_data(grid, fid, 'qtold', qtold, jdim=3)
        call read_dist_data(grid, fid, 'svlhx', svlhx, jdim=3)
        call read_dist_data(grid, fid, 'rhsav', rhsav, jdim=3)
        call read_dist_data(grid, fid, 'cldsav', cldsav, jdim=3)
#ifdef CLD_AER_CDNC
        call read_dist_data(grid, fid, 'oldnl', oldnl, jdim=3)
        call read_dist_data(grid, fid, 'oldni', oldni, jdim=3)
#endif
        call read_dist_data(grid, fid, 'airx', airx)
        call read_dist_data(grid, fid, 'lmc', lmc, jdim=3)
      end select
      return
      end subroutine new_io_clouds
#endif /* NEW_IO */
