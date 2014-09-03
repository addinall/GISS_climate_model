#include "rundeck_opts.h"
!@sum OCNML contains routines used for Qflux mixed layer,no deep diff.
!@auth G. Schmidt
!@ver  1.0

      SUBROUTINE CHECKO(SUBR)
!@sum  CHECKO Checks whether Ocean are reasonable
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm,focean
      USE STATIC_OCEAN, only : tocean
      USE DOMAIN_DECOMP_1D, only : GRID
      USE DOMAIN_DECOMP_1D, only : GET
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
      LOGICAL QCHECKO
      INTEGER I,J
      integer :: J_0, J_1, J_0H, J_1H, I_0, I_1, I_0H, I_1H, njpol
C****
C**** Extrack useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,  J_STRT_HALO = J_0H,
     *     J_STOP_HALO = J_1H) 
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      njpol = grid%J_STRT_SKP-grid%J_STRT

C**** Check for NaN/INF in ocean data
      CALL CHECK3C(TOCEAN(:,I_0:I_1,J_0:J_1),3,I_0,I_1,J_0,J_1,NJPOL,
     &     SUBR,'toc')

      QCHECKO = .FALSE.
C**** Check for reasonable values for ocean variables
      DO J=J_0, J_1
        DO I=I_0, I_1
          IF (FOCEAN(I,J).gt.0 .and. (TOCEAN(1,I,J).lt.-2. .or. 
     *         TOCEAN(1,I,J).gt.50.)) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,TOCEAN=',I,J,TOCEAN(1:3,I,J)
            QCHECKO = .TRUE.
          END IF
       END DO
      END DO
      IF (QCHECKO)
     *     call stop_model("CHECKO: Ocean variables out of bounds",255)

      END SUBROUTINE CHECKO

      SUBROUTINE io_ocean(kunit,iaction,ioerr)
!@sum  io_ocean reads and writes ocean arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE STATIC_OCEAN
      USE DOMAIN_DECOMP_1D, only : grid, GET, AM_I_ROOT
      USE DOMAIN_DECOMP_1D, only : PACK_COLUMN, PACK_DATA
      USE DOMAIN_DECOMP_1D, only : UNPACK_COLUMN, UNPACK_DATA
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCN01"
      REAL*8, DIMENSION(3,IM,JM) :: TOCEAN_glob
      REAL*8, DIMENSION(IM,JM)   :: Z1O_glob
      INTEGER :: J_0,J_1

      MODULE_HEADER(lhead+1:80) = 'R8 Tocn(3,im,jm),MixLD(im,jm)'

      CALL GET(grid, J_STRT=J_0,J_STOP=J_1)

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_COLUMN(grid, TOCEAN, TOCEAN_glob)
        CALL PACK_DATA(grid, Z1O, Z1O_glob)
        IF (AM_I_ROOT())
     &    WRITE (kunit,err=10) MODULE_HEADER,TOCEAN_glob,Z1O_glob
      CASE (IOREAD:)            ! input from restart file
        If (AM_I_ROOT()) Then
          READ (kunit,err=10) HEADER,TOCEAN_glob,Z1O_glob
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        END IF

        CALL UNPACK_COLUMN(grid, TOCEAN_glob, TOCEAN)
        CALL UNPACK_DATA(grid, Z1O_glob, Z1O)
c***        TOCEAN(:,:,J_0:J_1) = TOCEAN_glob(:,:,J_0:J_1)
c***        Z1O(:,J_0:J_1) = Z1O_glob(:,J_0:J_1)
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_ocean

#ifdef NEW_IO
      subroutine def_rsf_ocean(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      use static_ocean
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,tocean,'tocean(d3,dist_im,dist_jm)')
      call defvar(grid,fid,z1o,'z1o(dist_im,dist_jm)')
      return
      end subroutine def_rsf_ocean

      subroutine new_io_ocean(fid,iaction)
!@sum  new_io_ocean read/write ocean arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use static_ocean
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 'tocean', tocean, jdim=3)
        call write_dist_data(grid, fid, 'z1o', z1o)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'tocean', tocean, jdim=3)
        call read_dist_data(grid, fid, 'z1o', z1o)
      end select
      return
      end subroutine new_io_ocean

c
c the following are dummy routines
c
      subroutine def_rsf_ocdiag(fid)
      implicit none
      integer fid
      return
      end subroutine def_rsf_ocdiag
      subroutine new_io_ocdiag(fid,iaction)
      implicit none
      integer fid,iaction
      return
      end subroutine new_io_ocdiag
      subroutine set_ioptrs_ocnacc_default
      return
      end subroutine set_ioptrs_ocnacc_default
      subroutine set_ioptrs_ocnacc_extended
      return
      end subroutine set_ioptrs_ocnacc_extended
      subroutine def_meta_ocdiag(fid)
      implicit none
      integer fid
      return
      end subroutine def_meta_ocdiag
      subroutine write_meta_ocdiag(fid)
      implicit none
      integer fid
      return
      end subroutine write_meta_ocdiag
#endif /* NEW_IO */

      SUBROUTINE conserv_OCE(OCEANE)
!@sum  conserv_OCE calculates ocean energy for Qflux ocean
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : shw,rhows
      USE MODEL_COM, only : im,jm,fim,focean
      USE GEOM, only : imaxj
      USE STATIC_OCEAN, only : tocean,z1o,z12o
      USE DOMAIN_DECOMP_1D, only : GRID
      USE DOMAIN_DECOMP_1D, only : GET
      IMPLICIT NONE
!@var OCEANE ocean energy (J/M^2)
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: OCEANE
      INTEGER I,J
      integer :: J_0, J_1, I_0, I_1
      logical :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &          HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &          HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0, J_1
      DO I=I_0,IMAXJ(J)
        IF (FOCEAN(I,J).gt.0) THEN
          OCEANE(I,J)=(TOCEAN(1,I,J)*Z1O(I,J)
     *         +TOCEAN(2,I,J)*(Z12O(I,J)-Z1O(I,J)))*SHW*RHOWS
        ELSE
          OCEANE(I,J)=0
        END IF
      ENDDO
      ENDDO
      IF (HAVE_SOUTH_POLE) OCEANE(2:im,1) =OCEANE(1,1)
      IF (HAVE_NORTH_POLE) OCEANE(2:im,JM)=OCEANE(1,JM)
C****
      END SUBROUTINE conserv_OCE

      SUBROUTINE DUMMY_OCN
!@sum  DUMMY necessary entry points for non-dynamic/non-deep oceans
!@auth Gavin Schmidt
!@ver  1.0
      ENTRY ODIFS
      ENTRY io_ocdiag
      ENTRY reset_ODIAG
      ENTRY gather_odiags ()
      ENTRY diag_OCEAN
      ENTRY diag_OCEAN_prep

      ENTRY init_ODEEP
      ENTRY alloc_ODEEP

      RETURN
      END SUBROUTINE DUMMY_OCN
