#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif
      MODULE LAKES_COM
!@sum  LAKES_COM model variables for Lake/Rivers module
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : IM,JM,ioread,iowrite,lhead,irerun,irsfic
     *     ,irsficno
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE
      SAVE
!@var MWL mass of lake water (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: MWL
!@var GML total enthalpy of lake (J)
      REAL*8, ALLOCATABLE,  DIMENSION(:,:) :: GML
!@var TLAKE temperature of lake (C)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: TLAKE
!@var MLDLK mixed layer depth in lake (m)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: MLDLK
!@var FLAKE variable lake fraction (1)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: FLAKE
!@var TANLK tan(alpha) = slope for conical lake (1)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: TANLK
!@var SVFLAKE previous lake fraction (1)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: SVFLAKE

#ifdef TRACERS_WATER
!@var TRLAKE tracer amount in each lake level (kg)      
Crgr      REAL*8,  ALLOCATABLE, DIMENSION(NTM,2,:,:) :: TRLAKE
      REAL*8,  ALLOCATABLE, DIMENSION(:,:,:,:) :: TRLAKE
#endif

!@param NRVRMX Max No. of named rivers
      INTEGER, PARAMETER :: NRVRMX = 42
!@var NRVR actual No. of named rivers
      INTEGER :: NRVR
!@var IRVRMTH,JRVRMTH indexes for named river mouths
      INTEGER, DIMENSION(NRVRMX) :: IRVRMTH,JRVRMTH
!@var NAMERVR Named rivers
      CHARACTER*8, DIMENSION(NRVRMX) :: NAMERVR
!@var RVROUT Discharges from named rivers
      REAL*8, DIMENSION(NRVRMX) :: RVROUT

      END MODULE LAKES_COM


       SUBROUTINE ALLOC_LAKES_COM (GRID)
C23456789012345678901234567890123456789012345678901234567890123456789012
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth Raul Garza-Robles
!@ver  1.0
      USE DOMAIN_DECOMP_ATM, only: DIST_GRID, GET
      USE MODEL_COM, only : IM, JM
      USE LAKES_COM, ONLY: MWL, GML, TLAKE, MLDLK, FLAKE, TANLK, SVFLAKE
#ifdef TRACERS_WATER
      USE TRACER_COM, only : NTM
      USE LAKES_COM, ONLY:  TRLAKE
#endif
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid
      INTEGER :: I_0H,I_1H, J_0H,J_1H
      INTEGER IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE ( MWL(I_0H:I_1H,J_0H:J_1H),
     *           GML(I_0H:I_1H,J_0H:J_1H),
     *           TLAKE(I_0H:I_1H,J_0H:J_1H),
     *           MLDLK(I_0H:I_1H,J_0H:J_1H),
     *           FLAKE(I_0H:I_1H,J_0H:J_1H),
     *           TANLK(I_0H:I_1H,J_0H:J_1H),
     *           SVFLAKE(I_0H:I_1H,J_0H:J_1H),
     *           STAT=IER
     *            )

#ifdef TRACERS_WATER
!@var TRLAKE tracer amount in each lake level (kg)
      ALLOCATE( TRLAKE(NTM,2,I_0H:I_1H,J_0H:J_1H)
     * , STAT=IER)
#endif
      RETURN
      END SUBROUTINE ALLOC_LAKES_COM

      SUBROUTINE io_lakes(kunit,iaction,ioerr)
!@sum  io_lakes reads and writes lake arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT, grid
      USE DOMAIN_DECOMP_1D, only : PACK_DATA  , PACK_BLOCK
      USE DOMAIN_DECOMP_1D, only : UNPACK_DATA, UNPACK_BLOCK,
     *     BACKSPACE_PARALLEL
      USE LAKES_COM
      IMPLICIT NONE
      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "LAKE02"
      REAL*8, DIMENSION(IM,JM):: MLDLK_glob,MWL_glob,TLAKE_glob,GML_glob
     *     ,FLAKE_glob
#ifdef TRACERS_WATER
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRLAK01"
      REAL*8 :: TRLAKE_GLOB(NTM,2,IM,JM)
      IF (AM_I_ROOT())
     *   write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a7,i3,a)')'R8 dim(',NTM,',2,im,jm):TRLAKE'
#endif

      MODULE_HEADER(lhead+1:80) = 'R8 dim(im,jm):MixLD,MWtr,Tlk,Enth,Fl'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_DATA(grid, MLDLK, MLDLK_GLOB)
        CALL PACK_DATA(grid, MWL  ,   MWL_GLOB)
        CALL PACK_DATA(grid, TLAKE, TLAKE_GLOB)
        CALL PACK_DATA(grid, GML  ,   GML_GLOB)
        CALL PACK_DATA(grid, FLAKE, FLAKE_GLOB)
#ifdef TRACERS_WATER
        CALL PACK_BLOCK(grid, TRLAKE, TRLAKE_glob)
#endif
        IF (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,MLDLK_glob,MWL_glob,
     &                         TLAKE_glob,GML_glob,FLAKE_glob
#ifdef TRACERS_WATER
          WRITE (kunit,err=10) TRMODULE_HEADER,TRLAKE_glob
#endif
        END IF
      CASE (IOREAD:)            ! input from restart file
        if ( AM_I_ROOT() ) then
          READ (kunit,err=10) HEADER
          CALL BACKSPACE_PARALLEL(kunit)
          if (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            READ (kunit,err=10) HEADER,MLDLK_glob,MWL_glob,TLAKE_glob,
     &           GML_glob       ! no FLAKE
            FLAKE_glob = 0.d0
          else
            READ (kunit,err=10) HEADER,MLDLK_glob,MWL_glob,TLAKE_glob,
     &           GML_glob,FLAKE_glob

c          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
c            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
c            GO TO 10
          END IF
        end if
        CALL UNPACK_DATA(grid,  MLDLK_GLOB, MLDLK)
        CALL UNPACK_DATA(grid,    MWL_GLOB, MWL  )
        CALL UNPACK_DATA(grid,  TLAKE_GLOB, TLAKE)
        CALL UNPACK_DATA(grid,    GML_GLOB, GML  )
        CALL UNPACK_DATA(grid,  FLAKE_GLOB, FLAKE)

#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRERUN,IOREAD,IRSFIC,IRSFICNO)    ! reruns/restarts
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) TRHEADER,TRLAKE_glob
            IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRHEADER
     *             ,TRMODULE_HEADER
              GO TO 10
            END IF
          end if
          CALL UNPACK_BLOCK(grid, TRLAKE_GLOB, TRLAKE)
        END SELECT
#endif
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_lakes

#ifdef NEW_IO
      subroutine def_rsf_lakes(fid)
!@sum  def_rsf_lakes defines lake array structure in restart files
!@auth M. Kelley
!@ver  beta
      use lakes_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use conserv_diags
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,mldlk,'mldlk(dist_im,dist_jm)')
      call defvar(grid,fid,mwl,'mwl(dist_im,dist_jm)')
      call defvar(grid,fid,tlake,'tlake(dist_im,dist_jm)')
      call defvar(grid,fid,gml,'gml(dist_im,dist_jm)')
      call defvar(grid,fid,flake,'flake(dist_im,dist_jm)')
#ifdef TRACERS_WATER
      call defvar(grid,fid,trlake,'trlake(ntm,d2,dist_im,dist_jm)')
#endif
      call declare_conserv_diags( grid, fid, 'wliql(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'eliql(dist_im,dist_jm)' )
      return
      end subroutine def_rsf_lakes

      subroutine new_io_lakes(fid,iaction)
!@sum  new_io_lakes read/write lake arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use lakes_com
      use conserv_diags
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      external conserv_LKM, conserv_LKE
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 'mldlk', mldlk)
        call write_dist_data(grid, fid, 'mwl',   mwl)
        call write_dist_data(grid, fid, 'tlake', tlake)
        call write_dist_data(grid, fid, 'gml',   gml)
        call write_dist_data(grid, fid, 'flake', flake)
#ifdef TRACERS_WATER
        call write_dist_data(grid, fid, 'trlake', trlake, jdim=4)
#endif
        call dump_conserv_diags( grid, fid, 'wliql', conserv_LKM )
        call dump_conserv_diags( grid, fid, 'eliql', conserv_LKE )
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'mldlk', mldlk)
        call read_dist_data(grid, fid, 'mwl',   mwl)
        call read_dist_data(grid, fid, 'tlake', tlake)
        call read_dist_data(grid, fid, 'gml',   gml)
        call read_dist_data(grid, fid, 'flake', flake)
#ifdef TRACERS_WATER
        call read_dist_data(grid, fid, 'trlake', trlake, jdim=4)
#endif
      end select
      return
      end subroutine new_io_lakes

      subroutine def_meta_rvracc(fid)
!@sum  def_meta_rvracc defines river metadata in acc files
!@auth M. Kelley
!@ver  beta
      use lakes_com, only : nrvr,rvrout,namervr
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_attr
      implicit none
      integer :: fid         !@var fid file id
      call defvar(grid,fid,rvrout(1:nrvr),'rvr(nrvr)')
      call write_attr(grid,fid,'rvr','reduction','sum')
      call defvar(grid,fid,namervr(1:nrvr),'namervr(rvr_strlen,nrvr)')
      return
      end subroutine def_meta_rvracc

      subroutine write_meta_rvracc(fid)
!@sum  write_meta_rvracc write river accumulation metadata to file
!@auth M. Kelley
      use lakes_com, only : nrvr,rvrout,namervr
      use domain_decomp_atm, only : grid
      use pario, only : write_data
      implicit none
      integer fid   !@var fid unit number of read/write
      call write_data(grid,fid,'rvr',rvrout(1:nrvr))
      call write_data(grid,fid,'namervr',namervr(1:nrvr))
      return
      end subroutine write_meta_rvracc

#endif /* NEW_IO */
