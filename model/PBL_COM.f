#include "rundeck_opts.h"

      MODULE PBLCOM
!@sum  PBLCOM contains the arrays used by the Boundary Layer code
!@auth Greg Hartke/Ye Cheng
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm
#endif
      USE SOCPBL, only : npbl=>n
      IMPLICIT NONE
      SAVE
!@var uabl boundary layer profile for zonal wind
!@var vabl boundary layer profile for meridional wind
!@var tabl boundary layer profile for temperature
!@var qabl boundary layer profile for humidity
      real*8, allocatable, dimension(:,:,:,:) :: uabl,vabl,tabl,qabl
!@var eabl boundary layer profile for turbulent KE (calc. on sec. grid)
      real*8, allocatable, dimension(:,:,:,:) :: eabl
#ifdef TRACERS_ON
!@var trabl boundary layer profile for tracers
      real*8, allocatable, dimension(:,:,:,:,:) :: trabl
#endif

!@var cmgs drag coefficient (dimensionless surface momentum flux)
!@var chgs Stanton number   (dimensionless surface heat flux)
!@var cqgs Dalton number    (dimensionless surface moisture flux)
      real*8, allocatable, dimension(:,:,:) :: cmgs,chgs,cqgs

!@var ipbl flag for whether pbl properties were found at last timestep
      integer, allocatable, dimension(:,:,:) :: ipbl

!@var ROUGHL log10(zgs/roughness length), prescribed with zgs=30 m.
      REAL*8, allocatable, dimension(:,:) :: roughl
!@var WSAVG     COMPOSITE SURFACE WIND MAGNITUDE (M/S)
!@var TSAVG     COMPOSITE SURFACE AIR TEMPERATURE (K)
!@var QSAVG     COMPOSITE SURFACE AIR SPECIFIC HUMIDITY (1)
!@var DCLEV     LAYER TO WHICH DRY CONVECTION MIXES (1)
!@var USAVG     COMPOSITE SURFACE U WIND
!@var VSAVG     COMPOSITE SURFACE V WIND
!@var TAUAVG    COMPOSITE SURFACE MOMENTUM TRANSFER (TAU)
!@var w2_l1     COMPOSITE vertical component of t.k.e. at gcm layer 1
!@var USTAR_pbl friction velocity (sqrt of srfc mom flux) (m/s)
      REAL*8, allocatable, dimension(:,:) ::
     &     wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,tgvavg,qgavg
     &    ,w2_l1,gustiwind
      REAL*8, allocatable, dimension(:,:,:) :: ustar_pbl

!@var [tuv]1_after_aturb first-layer temp/winds after ATURB completes
!@+   (used to compute tendencies seen by the PBL code)
      REAL*8, allocatable, dimension(:,:) ::
     &     t1_after_aturb,u1_after_aturb,v1_after_aturb

!@var egcm  3-d turbulent kinetic energy in the whole atmosphere
!@var w2gcm vertical component of egcm
!@var t2gcm 3-d turbulent temperature variance in the whole atmosphere
      real*8, allocatable, dimension(:,:,:) :: egcm,w2gcm,t2gcm

!@var uflux surface turbulent u-flux (=-<uw>)
!@var vflux surface turbulent v-flux (=-<vw>)
!@var tflux surface turbulent t-flux (=-<tw>)
!@var qflux surface turbulent q-flux (=-<qw>)
      real*8, allocatable, dimension(:,:) :: uflux,vflux,tflux,qflux

      END MODULE PBLCOM

      SUBROUTINE io_pbl(kunit,iaction,ioerr)
!@sum  io_pbl reads and writes model variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,irsfic,irerun,iowrite,irsficno,lhead
      USE PBLCOM
      USE DOMAIN_DECOMP_1D, only : grid, GET, AM_I_ROOT
      USE DOMAIN_DECOMP_1D, only : pack_column, pack_data
      USE DOMAIN_DECOMP_1D, only : unpack_column, unpack_data
      USE DOMAIN_DECOMP_1D, only : pack_block , unpack_block
#ifdef TRACERS_ON
      use tracer_com, only : trname
#endif
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
      real*8, dimension(:,:,:,:), allocatable :: ! (npbl,4,IM,JM)
     &     uabl_glob,vabl_glob,tabl_glob,qabl_glob,eabl_glob
      REAL*8, DIMENSION(:,:,:), allocatable ::   ! (4,IM,JM)
     &     cmgs_glob, chgs_glob, cqgs_glob
      INTEGER, DIMENSION(:,:,:), allocatable ::  ! (4,IM,JM)
     &     ipbl_glob
      INTEGER :: J_0, J_1, j_0h, j_1h, n, i_0h, i_1h
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "PBL01"
#ifdef TRACERS_ON
!@var TR_HEADER Character string label for tracer record
      CHARACTER*80 :: TR_HEADER, TR_MODULE_HEADER = "TRPBL01"
      REAL*8, DIMENSION(:,:,:,:), allocatable :: trabl_glob,trabl_loc
#endif
      integer :: img, jmg

#ifdef PBL_USES_GCM_TENDENCIES
      call stop_model('io_pbl: need to save [tuv]1_after_aturb',255)
#endif

      write (MODULE_HEADER(lhead+1:80),'(a7,i2,a)') 'R8 dim(',npbl,
     *  ',4,ijm):Ut,Vt,Tt,Qt,Et dim(4,ijm,3):Cmhq, I:Ipb(4,ijm)'

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &     J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      if(am_i_root()) then
         img = IM
         jmg = JM
      else
         img = 1
         jmg = 1
      end if
      allocate(uabl_glob(npbl,4,img,jmg))
      allocate(vabl_glob(npbl,4,img,jmg))
      allocate(tabl_glob(npbl,4,img,jmg))
      allocate(qabl_glob(npbl,4,img,jmg))
      allocate(eabl_glob(npbl,4,img,jmg))
      allocate(cmgs_glob(4,img,jmg))
      allocate(chgs_glob(4,img,jmg))
      allocate(cqgs_glob(4,img,jmg))
      allocate(ipbl_glob(4,img,jmg))
#ifdef TRACERS_ON
      allocate(trabl_glob(npbl,4,im,jm))
#endif
#ifdef TRACERS_ON
      allocate(trabl_loc(npbl,4,i_0h:i_1h,j_0h:j_1h))
#endif

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_BLOCK(grid, uabl, uabl_glob)
        CALL PACK_BLOCK(grid, vabl, vabl_glob)
        CALL PACK_BLOCK(grid, tabl, tabl_glob)
        CALL PACK_BLOCK(grid, qabl, qabl_glob)
        CALL PACK_BLOCK(grid, eabl, eabl_glob)

        CALL PACK_COLUMN(grid, cmgs, cmgs_glob)
        CALL PACK_COLUMN(grid, chgs, chgs_glob)
        CALL PACK_COLUMN(grid, cqgs, cqgs_glob)
        CALL PACK_COLUMN(grid, ipbl, ipbl_glob)

        IF (AM_I_ROOT()) THEN
          WRITE (KUNIT,ERR=10) MODULE_HEADER,UABL_GLOB,VABL_GLOB
     *       ,TABL_GLOB,QABL_GLOB,EABL_GLOB,CMGS_GLOB
     *       ,CHGS_GLOB,CQGS_GLOB,IPBL_GLOB
        END IF
#ifdef TRACERS_ON
        do n=1,ntm
          trabl_loc(:,:,:,:) = trabl(:,n,:,:,:)
          CALL PACK_BLOCK(grid, trabl_loc, trabl_glob)
          write (TR_MODULE_HEADER(lhead+1:80),'(a8,a8,i2,a)')
     &         trname(n),
     &         ' R8 dim(',npbl,',4,ijm):TRt'
          IF (AM_I_ROOT()) write(kunit,err=10)
     &         TR_MODULE_HEADER,TRABL_GLOB
        enddo
#endif

      CASE (IOREAD:)            ! input from restart file or restart
        if ( AM_I_ROOT() ) then
          READ (KUNIT,ERR=10) HEADER,UABL_glob,VABL_glob,TABL_glob,
     &         QABL_glob,EABL_glob,CMGS_glob,CHGS_glob,CQGS_glob,
     &         IPBL_glob
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        end if

        call UNPACK_BLOCK(grid, uabl_glob, uabl)
        call UNPACK_BLOCK(grid, vabl_glob, vabl)
        call UNPACK_BLOCK(grid, tabl_glob, tabl)
        call UNPACK_BLOCK(grid, qabl_glob, qabl)
        call UNPACK_BLOCK(grid, eabl_glob, eabl)

        call UNPACK_COLUMN(grid, cmgs_glob, cmgs)
        call UNPACK_COLUMN(grid, chgs_glob, chgs)
        call UNPACK_COLUMN(grid, cqgs_glob, cqgs)
        call UNPACK_COLUMN(grid, ipbl_glob, ipbl)
#ifdef TRACERS_ON
        SELECT CASE (IACTION)
        CASE (IOREAD,IRERUN,IRSFIC,IRSFICNO)    ! restarts
        do n=1,ntm
          IF (AM_I_ROOT()) then
            read(kunit,err=10) TR_HEADER,TRABL_GLOB
            IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in tracer module version ",TR_HEADER
     *             ,TR_MODULE_HEADER
              GO TO 10
            END IF
          ENDIF
          CALL UNPACK_BLOCK(grid, trabl_glob, trabl_loc)
          trabl(:,n,:,:,:) = trabl_loc(:,:,:,:)
        enddo
        END SELECT
#endif
      END SELECT
      call freemem
      RETURN
 10   IOERR=1
      call freemem
      RETURN
      contains
      subroutine freemem
      deallocate(uabl_glob)
      deallocate(vabl_glob)
      deallocate(tabl_glob)
      deallocate(qabl_glob)
      deallocate(eabl_glob)
      deallocate(cmgs_glob)
      deallocate(chgs_glob)
      deallocate(cqgs_glob)
      deallocate(ipbl_glob)
#ifdef TRACERS_ON
        deallocate(trabl_glob)
#endif
#ifdef TRACERS_ON
      deallocate(trabl_loc)
#endif
      end subroutine freemem
      END SUBROUTINE io_pbl

      SUBROUTINE io_bldat(kunit,iaction,ioerr)
!@sum  io_bldat reads and writes boundary layer data to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE DOMAIN_DECOMP_1D, only : GET, grid, AM_I_ROOT
      USE DOMAIN_DECOMP_1D, only : UNPACK_DATA, UNPACK_COLUMN
      USE DOMAIN_DECOMP_1D, only : PACK_DATA, PACK_COLUMN
      USE PBLCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "BLD02"
      INTEGER :: J_0, J_1
      REAL*8, DIMENSION(IM,JM) :: wsavg_glob,tsavg_glob,
     *       qsavg_glob,dclev_glob,usavg_glob,
     *       vsavg_glob,tauavg_glob,qgavg_glob,tgvavg_glob
      REAL*8, DIMENSION(LM,IM,JM) :: egcm_glob, w2gcm_glob
      REAL*8, DIMENSION(4,IM,JM) :: ustar_pbl_glob

      MODULE_HEADER(lhead+1:80) = 'R8 dim(ijm):ws,ts,qs,'//
     *  'LvlDC,us,vs,tau,u*(4,.),ke;w2(lijm),tgv,qg'

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_DATA(grid, wsavg, wsavg_glob)
        CALL PACK_DATA(grid, tsavg, tsavg_glob)
        CALL PACK_DATA(grid, qsavg, qsavg_glob)
        CALL PACK_DATA(grid, dclev, dclev_glob)
        CALL PACK_DATA(grid, usavg, usavg_glob)
        CALL PACK_DATA(grid, vsavg, vsavg_glob)
        CALL PACK_DATA(grid, tauavg, tauavg_glob)
        CALL PACK_DATA(grid, tgvavg, tgvavg_glob)
        CALL PACK_DATA(grid, qgavg, qgavg_glob)

        CALL PACK_COLUMN(grid, ustar_pbl, ustar_pbl_glob)
        CALL PACK_COLUMN(grid, egcm, egcm_glob)
        CALL PACK_COLUMN(grid, w2gcm, w2gcm_glob)

        IF (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,wsavg_glob,tsavg_glob
     *         ,qsavg_glob,dclev_glob,usavg_glob,vsavg_glob,tauavg_glob
     *         ,ustar_pbl_glob,egcm_glob,w2gcm_glob,tgvavg_glob
     *         ,qgavg_glob
        END IF

      CASE (IOREAD:)            ! input from restart file
        if ( AM_I_ROOT() ) then
          READ (kunit,err=10) HEADER,wsavg_glob,tsavg_glob,
     *       qsavg_glob,dclev_glob,usavg_glob,
     *       vsavg_glob,tauavg_glob,ustar_pbl_glob,egcm_glob,
     *       w2gcm_glob,tgvavg_glob,qgavg_glob
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        end if

        call UNPACK_DATA(grid, wsavg_glob, wsavg)
        call UNPACK_DATA(grid, tsavg_glob, tsavg)
        call UNPACK_DATA(grid, qsavg_glob, qsavg)
        call UNPACK_DATA(grid, dclev_glob, dclev)
        call UNPACK_DATA(grid, usavg_glob, usavg)
        call UNPACK_DATA(grid, vsavg_glob, vsavg)
        call UNPACK_DATA(grid, tauavg_glob, tauavg)
        call UNPACK_DATA(grid, tgvavg_glob, tgvavg)
        call UNPACK_DATA(grid, qgavg_glob, qgavg)

        call UNPACK_COLUMN(grid, ustar_pbl_glob, ustar_pbl)
        call UNPACK_COLUMN(grid, egcm_glob, egcm)
        call UNPACK_COLUMN(grid, w2gcm_glob, w2gcm)

      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_bldat

#ifdef NEW_IO
      subroutine def_rsf_pbl(fid)
!@sum  def_rsf_pbl defines pbl array structures in restart files
!@auth M. Kelley
!@ver  beta
      use pblcom
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      character(len=29) :: dimstr
      dimstr='(npbl,nstype,dist_im,dist_jm)'
      call defvar(grid,fid,uabl,'uabl'//dimstr)
      call defvar(grid,fid,vabl,'vabl'//dimstr)
      call defvar(grid,fid,tabl,'tabl'//dimstr)
      call defvar(grid,fid,qabl,'qabl'//dimstr)
      call defvar(grid,fid,eabl,'eabl'//dimstr)
      dimstr='(nstype,dist_im,dist_jm)'
      call defvar(grid,fid,cmgs,'cmgs'//dimstr)
      call defvar(grid,fid,chgs,'chgs'//dimstr)
      call defvar(grid,fid,cqgs,'cqgs'//dimstr)
      call defvar(grid,fid,ipbl,'ipbl'//dimstr)
#ifdef TRACERS_ON
      call defvar(grid,fid,trabl,
     &     'trabl(npbl,ntm,nstype,dist_im,dist_jm)')
#endif
      dimstr='(dist_im,dist_jm)'
      call defvar(grid,fid,t1_after_aturb,'t1_after_aturb'//dimstr)
      call defvar(grid,fid,u1_after_aturb,'u1_after_aturb'//dimstr)
      call defvar(grid,fid,v1_after_aturb,'v1_after_aturb'//dimstr)
      return
      end subroutine def_rsf_pbl

      subroutine new_io_pbl(fid,iaction)
!@sum  new_io_pbl read/write pbl arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use pblcom
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 'uabl', uabl, jdim=4)
        call write_dist_data(grid, fid, 'vabl', vabl, jdim=4)
        call write_dist_data(grid, fid, 'tabl', tabl, jdim=4)
        call write_dist_data(grid, fid, 'qabl', qabl, jdim=4)
        call write_dist_data(grid, fid, 'eabl', eabl, jdim=4)
        call write_dist_data(grid, fid, 'cmgs', cmgs, jdim=3)
        call write_dist_data(grid, fid, 'chgs', chgs, jdim=3)
        call write_dist_data(grid, fid, 'cqgs', cqgs, jdim=3)
        call write_dist_data(grid, fid, 'ipbl', ipbl, jdim=3)
#ifdef TRACERS_ON
        call write_dist_data(grid, fid, 'trabl', trabl, jdim=5)
#endif
        call write_dist_data(grid, fid, 't1_after_aturb',t1_after_aturb)
        call write_dist_data(grid, fid, 'u1_after_aturb',u1_after_aturb)
        call write_dist_data(grid, fid, 'v1_after_aturb',v1_after_aturb)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'uabl', uabl, jdim=4)
        call read_dist_data(grid, fid, 'vabl', vabl, jdim=4)
        call read_dist_data(grid, fid, 'tabl', tabl, jdim=4)
        call read_dist_data(grid, fid, 'qabl', qabl, jdim=4)
        call read_dist_data(grid, fid, 'eabl', eabl, jdim=4)
        call read_dist_data(grid, fid, 'cmgs', cmgs, jdim=3)
        call read_dist_data(grid, fid, 'chgs', chgs, jdim=3)
        call read_dist_data(grid, fid, 'cqgs', cqgs, jdim=3)
        call read_dist_data(grid, fid, 'ipbl', ipbl, jdim=3)
#ifdef TRACERS_ON
        call read_dist_data(grid, fid, 'trabl', trabl, jdim=5)
#endif
        call read_dist_data(grid, fid, 't1_after_aturb',t1_after_aturb)
        call read_dist_data(grid, fid, 'u1_after_aturb',u1_after_aturb)
        call read_dist_data(grid, fid, 'v1_after_aturb',v1_after_aturb)
      end select
      return
      end subroutine new_io_pbl

      subroutine def_rsf_bldat(fid)
!@sum  def_rsf_bldat defines bldat array structure in restart files
!@auth M. Kelley
!@ver  beta
      use pblcom
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,wsavg,'wsavg(dist_im,dist_jm)')
      call defvar(grid,fid,tsavg,'tsavg(dist_im,dist_jm)')
      call defvar(grid,fid,qsavg,'qsavg(dist_im,dist_jm)')
      call defvar(grid,fid,dclev,'dclev(dist_im,dist_jm)')
      call defvar(grid,fid,usavg,'usavg(dist_im,dist_jm)')
      call defvar(grid,fid,vsavg,'vsavg(dist_im,dist_jm)')
      call defvar(grid,fid,tauavg,'tauavg(dist_im,dist_jm)')
      call defvar(grid,fid,tgvavg,'tgvavg(dist_im,dist_jm)')
      call defvar(grid,fid,qgavg,'qgavg(dist_im,dist_jm)')
      call defvar(grid,fid,ustar_pbl,
     &     'ustar_pbl(nstype,dist_im,dist_jm)')
      call defvar(grid,fid,egcm,'egcm(lm,dist_im,dist_jm)')
      call defvar(grid,fid,w2gcm,'w2gcm(lm,dist_im,dist_jm)')
      return
      end subroutine def_rsf_bldat

      subroutine new_io_bldat(fid,iaction)
!@sum  new_io_bldat read/write bldat arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use pblcom
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'wsavg',wsavg)
        call write_dist_data(grid,fid,'tsavg',tsavg)
        call write_dist_data(grid,fid,'qsavg',qsavg)
        call write_dist_data(grid,fid,'dclev',dclev)
        call write_dist_data(grid,fid,'usavg',usavg)
        call write_dist_data(grid,fid,'vsavg',vsavg)
        call write_dist_data(grid,fid,'tauavg',tauavg)
        call write_dist_data(grid,fid,'tgvavg',tgvavg)
        call write_dist_data(grid,fid,'qgavg',qgavg)
        call write_dist_data(grid,fid,'ustar_pbl',ustar_pbl,jdim=3)
        call write_dist_data(grid,fid,'egcm',egcm, jdim=3)
        call write_dist_data(grid,fid,'w2gcm',w2gcm, jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'wsavg',wsavg)
        call read_dist_data(grid,fid,'tsavg',tsavg)
        call read_dist_data(grid,fid,'qsavg',qsavg)
        call read_dist_data(grid,fid,'dclev',dclev)
        call read_dist_data(grid,fid,'usavg',usavg)
        call read_dist_data(grid,fid,'vsavg',vsavg)
        call read_dist_data(grid,fid,'tauavg',tauavg)
        call read_dist_data(grid,fid,'tgvavg',tgvavg)
        call read_dist_data(grid,fid,'qgavg',qgavg)
        call read_dist_data(grid,fid,'ustar_pbl',ustar_pbl,jdim=3)
        call read_dist_data(grid,fid,'egcm',egcm, jdim=3)
        call read_dist_data(grid,fid,'w2gcm',w2gcm, jdim=3)
      end select
      return
      end subroutine new_io_bldat
#endif /* NEW_IO */

      SUBROUTINE ALLOC_PBL_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE PBLCOM
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, GET
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      ALLOCATE(    uabl(npbl,4,    I_0H:I_1H,J_0H:J_1H),
     *             vabl(npbl,4,    I_0H:I_1H,J_0H:J_1H),
     *             tabl(npbl,4,    I_0H:I_1H,J_0H:J_1H),
     *             qabl(npbl,4,    I_0H:I_1H,J_0H:J_1H),
     *             eabl(npbl,4,    I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)
      qabl=0.  ! initialise to make life easier

#ifdef TRACERS_ON
      ALLOCATE(    trabl(npbl,ntm,4,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)
#endif

      ALLOCATE(    cmgs(4,I_0H:I_1H,J_0H:J_1H),
     *             chgs(4,I_0H:I_1H,J_0H:J_1H),
     *             cqgs(4,I_0H:I_1H,J_0H:J_1H),
     *             ipbl(4,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)
      ipbl(:,:,J_0H:J_1H) = 0

      ALLOCATE(    roughl(I_0H:I_1H,J_0H:J_1H),
     *              wsavg(I_0H:I_1H,J_0H:J_1H),
     *              tsavg(I_0H:I_1H,J_0H:J_1H),
     *              qsavg(I_0H:I_1H,J_0H:J_1H),
     *              dclev(I_0H:I_1H,J_0H:J_1H),
     *              usavg(I_0H:I_1H,J_0H:J_1H),
     *              vsavg(I_0H:I_1H,J_0H:J_1H),
     *             tauavg(I_0H:I_1H,J_0H:J_1H),
     *             tgvavg(I_0H:I_1H,J_0H:J_1H),
     *              qgavg(I_0H:I_1H,J_0H:J_1H),
     *              w2_l1(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)
      w2_l1(:,J_0H:J_1H) =0.

      ALLOCATE(    ustar_pbl(4,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(    egcm(lm,I_0H:I_1H,J_0H:J_1H),
     *            w2gcm(lm,I_0H:I_1H,J_0H:J_1H),
     *            t2gcm(lm,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(   uflux(I_0H:I_1H,J_0H:J_1H),
     *            vflux(I_0H:I_1H,J_0H:J_1H),
     *            tflux(I_0H:I_1H,J_0H:J_1H),
     *            qflux(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)
      uflux(:,J_0H:J_1H) = 0. ! define defaults
      vflux(:,J_0H:J_1H) = 0. 
      tflux(:,J_0H:J_1H) = 0. 
      qflux(:,J_0H:J_1H) = 0. 

      ALLOCATE( gustiwind(I_0H:I_1H,J_0H:J_1H),STAT=IER)
      gustiwind(:,J_0H:J_1H) = 0.

      ALLOCATE(t1_after_aturb(I_0H:I_1H,J_0H:J_1H),
     &         u1_after_aturb(I_0H:I_1H,J_0H:J_1H),
     &         v1_after_aturb(I_0H:I_1H,J_0H:J_1H))
      t1_after_aturb(:,J_0H:J_1H) = 0.
      u1_after_aturb(:,J_0H:J_1H) = 0.
      v1_after_aturb(:,J_0H:J_1H) = 0.

      END SUBROUTINE ALLOC_PBL_COM

