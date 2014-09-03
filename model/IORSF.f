#include "rundeck_opts.h"

      SUBROUTINE io_rsf(filenm,it,iaction,ioerr)
!@sum   io_rsf controls the reading and writing of the restart files
!@auth  Gavin Schmidt
!@ver   1.0
!@calls io_model,io_ocean,io_lakes,io_seaice,io_earth,io_soils,io_snow
!@+     io_landice,io_bldat,io_pbl,io_clouds,io_somtq,io_rad,io_diags
!@+     io_ocdiag,io_icedyn,io_icdiag
      USE FILEMANAGER, only : openunit,closeunit
      USE DOMAIN_DECOMP_1D, only : am_i_root !REWIND_PARALLEL
      USE MODEL_COM, only : ioread_single,iowrite_single,Kradia
     *                     ,ioread,ioread_nodiag,iowrite

      IMPLICIT NONE
!@var filenm name of file to be read or written
      character(len=*) :: filenm
!@var iaction flag for reading or writing rsf file
      INTEGER, INTENT(IN) :: iaction
!@var it hour of model run
      INTEGER, INTENT(INOUT) :: it
!@var IOERR (1,0,-1) if there (is, is maybe, is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var IT1 hour for correct reading check
!@var ITM maximum hour for post-processing
!@var kunit Fortran unit number of file i/o
      INTEGER IT1,itm,iact,kunit
      logical skip_diag

      iact=iaction   ; skip_diag=.false.
      if(iaction==ioread_nodiag) then
         iact=ioread ; skip_diag=.true.
      end if

      ioerr=-1

      if(iaction.le.iowrite) then
c open output files with status = 'UNKNOWN'
        if(am_i_root()) ! only root PE can write
     &       call openunit(trim(filenm),kunit,.true.,.false.)
      elseif(iaction.ge.ioread) then
c open input files with status = 'OLD'
c all PEs need to read the label records, so no am_i_root check
        call openunit(trim(filenm),kunit,.true.,.true.)
      endif
c      if (kunit.gt.0) CALL REWIND_PARALLEL( kunit )

C**** For all iaction < 0  ==> WRITE, For all iaction > 0  ==> READ
C**** Particular values may produce variations in indiv. i/o routines

C**** Calls to individual i/o routines
      call io_label  (kunit,it,itm,iact,ioerr)
      it1=it
      if (Kradia>9) go to 10 ! used for testing daily stuff
      if (Kradia.gt.0) then
        if (iaction.ne.ioread_single .and.
     *   iaction.ne.iowrite_single) call io_rad (kunit,iact,ioerr)
        call io_diags  (kunit,it,iact,ioerr)
        if(it1.ne.it .or. ioerr.eq.1)
     &       call stop_model('restart problem',255)
        go to 10
      end if
      if(iaction.ne.ioread_single.and.iaction.ne.iowrite_single) then
        call io_model  (kunit,iact,ioerr)
        call io_strat  (kunit,iact,ioerr)
        call io_ocean  (kunit,iact,ioerr)
        call io_lakes  (kunit,iact,ioerr)
        call io_seaice (kunit,iact,ioerr)
        call io_earth  (kunit,iact,ioerr)
        call io_soils  (kunit,iact,ioerr)
        call io_vegetation  (kunit,iact,ioerr)
#ifdef USE_ENT
        !!! actually not sure if this call is needed
        !!! (seems like it is duplicated in io_vegetation...)
        call io_veg_related  (kunit,iact,ioerr)
        !call io_ent    (kunit,iaction,ioerr)
#endif
        call io_snow   (kunit,iact,ioerr)
        call io_landice(kunit,iact,ioerr)
        call io_bldat  (kunit,iact,ioerr)
        call io_pbl    (kunit,iact,ioerr)
        call io_clouds (kunit,iact,ioerr)
        call io_somtq  (kunit,iact,ioerr)
        call io_rad    (kunit,iact,ioerr)
        call io_icedyn (kunit,iact,ioerr)
#ifdef CALCULATE_FLAMMABILITY
        call io_flammability(kunit,iact,ioerr)
#endif
#ifdef TRACERS_ON
        call io_tracer (kunit,iact,ioerr)
#endif
      end if
      if(skip_diag) go to 10
      call io_diags  (kunit,it,iaction,ioerr)
      call io_ocdiag (kunit,it,iaction,ioerr)
      call io_icdiag (kunit,it,iaction,ioerr)
#ifdef TRACERS_ON
      call io_trdiag (kunit,it,iaction,ioerr)
#endif
      if (it1.ne.it) THEN
        WRITE(6,*) "TIMES DO NOT MATCH READING IN RSF FILE",it,it1
        ioerr=1
      END IF
      if (ioerr.eq.1) WRITE(6,*) "I/O ERROR IN RESTART FILE: KUNIT="
     *     ,kunit

C**** return maximum time
   10 it=itm

      if(am_i_root() .or. iaction.ge.ioread) call closeunit(kunit)

      RETURN
      END SUBROUTINE io_rsf

      subroutine read_ground_ic
!@sum   read_ground_ic read initial conditions file for
!@+     sea ice, land surface, land ice.  Transplanted from INPUT.
!@auth  M. Kelley
!@ver   1.0
!@calls io_seaice,io_earth,io_soils,io_landice
      use model_com, only : ioreadnt
      use filemanager, only : openunit,closeunit
      use domain_decomp_1d, only : am_i_root
      implicit none
      integer :: iu_GIC,ioerr
      call openunit("GIC",iu_GIC,.true.,.true.)
      ioerr=-1
      read(iu_GIC) ! ignore first line (ocean ic done in init_OCEAN)
      call io_seaice (iu_GIC,ioreadnt,ioerr)
      call io_earth  (iu_GIC,ioreadnt,ioerr)
      call io_soils  (iu_GIC,ioreadnt,ioerr)
      call io_landice(iu_GIC,ioreadnt,ioerr)
      if (ioerr.eq.1) then
        if (AM_I_ROOT())
     *       WRITE(6,*) "I/O ERROR IN GIC FILE: KUNIT=",iu_GIC
        call stop_model("INPUT: GIC READ IN ERROR",255)
      end if
      call closeunit (iu_GIC)
      return
      end subroutine read_ground_ic

      subroutine find_later_rsf(kdisk)
!@sum set kdisk such that Itime in rsf_file_name(kdisk) is
!@+   the larger of the Itimes in rsf_file_name(1:2).
!@+   Transplanted from INPUT.
      USE FILEMANAGER, only : openunit,closeunit
      use model_com, only : rsf_file_name
      implicit none
      integer, intent(out) :: kdisk
      integer :: Itime1,Itime2,kunit
      Itime1=-1
      call openunit(rsf_file_name(1),kunit,.true.,.true.)
      READ (kunit,ERR=410) Itime1
 410  continue
      call closeunit(kunit)
      Itime2=-1
      call openunit(rsf_file_name(2),kunit,.true.,.true.)
      READ (kunit,ERR=420) Itime2
 420  continue
      call closeunit(kunit)
      if (Itime1+Itime2.LE.-2) then
        call stop_model(
     &       'FIND_LATER_RSF: ERRORS ON BOTH RESTART FILES',255)
      endif
      KDISK=1
      IF (Itime2.GT.Itime1) KDISK=2
      return
      end subroutine find_later_rsf
