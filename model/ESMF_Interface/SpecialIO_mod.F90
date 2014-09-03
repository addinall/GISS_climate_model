#include "rundeck_opts.h"
#ifdef MPI_DEFS_HACK
#include "mpi_defs.h"
#endif

#if ( defined USE_ESMF )  || ( defined USE_MPP )
#define USE_MPI
#endif

module SpecialIO_mod
  use dist_grid_mod
  use GatherScatter_mod
#ifdef USE_ESMF
  use ESMF_Mod
#endif
  use MpiSupport_mod

  implicit none
  private
  public :: BACKSPACE_PARALLEL
  public :: REWIND_PARALLEL
  public :: SKIP_PARALLEL
  public :: DREAD_PARALLEL,DREAD8_PARALLEL,DWRITE8_PARALLEL
  public :: MREAD_PARALLEL
  public :: READT_PARALLEL,WRITET_PARALLEL,READT_PARALLEL_COLUMN
  public :: READT8_PARALLEL,READT8_COLUMN,WRITET8_COLUMN
  public :: READ_PARALLEL
  public :: WRITE_PARALLEL
  public :: WRITEI_PARALLEL
  public :: WRITEI8_PARALLEL
  
  interface DREAD_PARALLEL
    module procedure DREAD_PARALLEL_2D
    module procedure DREAD_PARALLEL_3D
  end interface
  
  interface DREAD8_PARALLEL
    !        MODULE PROCEDURE DREAD_PARALLEL_2D
    module procedure DREAD8_PARALLEL_3D
  end interface
  
  interface DWRITE8_PARALLEL
    module procedure DWRITE8_PARALLEL_3D
  end interface
  
  interface MREAD_PARALLEL
    module procedure MREAD_PARALLEL_2D
    module procedure MREAD_PARALLEL_3D
  end interface
  
  interface READT_PARALLEL
    module procedure READT_PARALLEL_2D
    module procedure READT_PARALLEL_3D
  end interface
  
  interface WRITET_PARALLEL
    module procedure WRITET_PARALLEL_2D
    !        MODULE PROCEDURE READT_PARALLEL_3D
  end interface
  
  interface READT_PARALLEL_COLUMN
    module procedure READT_PARALLEL_COLUMN_3D
  end interface
  
  interface READT8_PARALLEL
    !        MODULE PROCEDURE READT8_PARALLEL_2D
    module procedure READT8_PARALLEL_3D
  end interface
  
  interface READT8_COLUMN
    module procedure READT8_COLUMN_3D
  end interface
  
  interface WRITET8_COLUMN
    module procedure WRITET8_COLUMN_3D
  end interface
  
  interface READ_PARALLEL
    module procedure READ_PARALLEL_INTEGER_0
    module procedure READ_PARALLEL_INTEGER_1
    module procedure READ_PARALLEL_REAL8_1
  end interface
  
  interface WRITE_PARALLEL
    module procedure WRITE_PARALLEL_INTEGER_0
    module procedure WRITE_PARALLEL_INTEGER_1
    module procedure WRITE_PARALLEL_REAL8_0
    module procedure WRITE_PARALLEL_REAL8_1
    module procedure WRITE_PARALLEL_STRING_0
    module procedure WRITE_PARALLEL_STRING_1
    
  end interface
  
  interface WRITEI_PARALLEL
    module procedure WRITEI_PARALLEL_2D
  end interface

  interface WRITEI8_PARALLEL
    module procedure WRITEI8_PARALLEL_3D
  end interface
  
#ifdef USE_MPI
  include 'mpif.h'
#endif
contains
  
  subroutine BACKSPACE_PARALLEL(IUNIT)
    integer, intent(IN) :: IUNIT

    if (AM_I_ROOT()) backspace IUNIT

  end subroutine BACKSPACE_PARALLEL
  
  subroutine REWIND_PARALLEL(IUNIT)
    integer, intent(IN) :: IUNIT

    if (AM_I_ROOT()) rewind IUNIT
    
  end subroutine REWIND_PARALLEL
  
  subroutine SKIP_PARALLEL(IUNIT)
    integer, intent(IN) :: IUNIT
    
    if (AM_I_ROOT()) read(IUNIT)
    
  end subroutine SKIP_PARALLEL
  
  subroutine DREAD_PARALLEL_2D (grd_dum,IUNIT,NAME,AVAR, &
       &     recs_to_skip)
!@sum DREAD_PARALLEL  Parallel version of UTILDBL.f:DREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
    type (DIST_GRID),  intent(IN) :: grd_dum
    integer,      intent(IN)  :: IUNIT     !@var  IUNIT file unit number
    character*(*), intent(IN)  :: NAME      !@var  NAME  name of file being read
    real*8,       intent(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:) !@var  AOUT real*8 array
    integer, intent(IN), optional :: recs_to_skip
    real*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD)  !@var  AIN  real*4 array
    real*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var  AOUT real*8 array
    
    integer :: N                        !@var  N loop variable
    integer :: IM,JM
    integer :: rc
    ! now local
    
    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD
    !      write(*,*) "DD dread parallel IM,JM",IM,JM
    if (AM_I_ROOT()) then
      if(present(recs_to_skip)) then
        do n=1,recs_to_skip
          read (IUNIT,IOSTAT=rc)
        enddo
      endif
      read (IUNIT,IOSTAT=rc) AIN
      !****  convert from real*4 to real*8
      AOUT=AIN
    endif
    
    call scatterReal8(grd_dum, AOUT, AVAR, shape(AVAR), 2)
    call checkReadStatus(rc, name, 'DREAD_PARALLEL_2D')

  END SUBROUTINE DREAD_PARALLEL_2D

  SUBROUTINE DREAD_PARALLEL_3D (grd_dum,IUNIT,NAME,AVAR, &
       &     recs_to_skip)
!@sum DREAD_PARALLEL  Parallel version of UTILDBL.f:DREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT       !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME        !@var  NAME  name of file being read
    REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:,:) !@var  AOUT real*8 array
    INTEGER, INTENT(IN), OPTIONAL :: recs_to_skip
    REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3))  !@var  AIN  real*4 array
    REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var  AOUT real*8 array
    INTEGER :: N                        !@var  N loop variable
    INTEGER :: IM,JM,NM
    INTEGER :: rc
    ! now local

    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD
    NM   = size(AVAR,3)
    !      write(*,*) "DD dread parallel IM,JM,NM",IM,JM,NM

    If (AM_I_ROOT()) then
      if(present(recs_to_skip)) then
        do n=1,recs_to_skip
          READ (IUNIT,IOSTAT=rc)
        enddo
      endif
      READ (IUNIT,IOSTAT=rc) AIN
      !****  convert from real*4 to real*8
      AOUT=AIN
    EndIf

    Call scatterReal8(grd_dum, AOUT, AVAR, shape(AVAR), 2)
    call checkReadStatus(rc, name, 'DREAD_PARALLEL_3D')

  END SUBROUTINE DREAD_PARALLEL_3D

  SUBROUTINE DREAD8_PARALLEL_3D (grd_dum,IUNIT,NAME,AVAR, &
       &     recs_to_skip)
!@sum DREAD_PARALLEL  read an array real*8 avar(im,jm,:)
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT       !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME        !@var  NAME  name of file being read
    REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:,:) !@var AVAR real*8 array
    INTEGER, INTENT(IN), OPTIONAL :: recs_to_skip
    REAL*8 :: AGLOB(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var AGLOB global array
    INTEGER :: N,rc

    If (AM_I_ROOT()) then
      if(present(recs_to_skip)) then
        do n=1,recs_to_skip
          READ (IUNIT,IOSTAT=rc)
        enddo
      endif
      READ (IUNIT,IOSTAT=rc) AGLOB
    EndIf

    Call scatterReal8(grd_dum, AGLOB, AVAR, shape(AVAR), 2)
    call checkReadStatus(rc, name, 'DREAD8_PARALLEL_3D')

  END SUBROUTINE DREAD8_PARALLEL_3D

  SUBROUTINE MREAD_PARALLEL_2D (grd_dum,IUNIT,NAME,M,AVAR)
!@sum MREAD_PARALLEL  Parallel version of UTILDBL.f:MREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT     !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME      !@var  NAME  name of file being read
    INTEGER,      INTENT(OUT) :: M         !@var  M      initial integer
    REAL*8,      INTENT(OUT)  :: AVAR(:,grd_dum%J_STRT_HALO:) !@var  AOUT real*8 array
    REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD)  !@var  AIN  real*4 array
    REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var  AOUT real*8 array
    INTEGER :: N                        !@var  N loop variable
    INTEGER :: IM,JM
    INTEGER :: rc
    ! now local
    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD

    If (AM_I_ROOT()) then
      READ (IUNIT,IOSTAT=rc) M, AIN
      !****  convert from real*4 to real*8
      AOUT=AIN
    EndIf

    Call scatterReal8(grd_dum, AOUT, AVAR, shape(AVAR), 2)
#ifdef USE_MPI
    CALL ESMF_BCAST(grd_dum, M   )
#endif

    call checkReadStatus(rc, name, 'MREAD_PARALLEL_2D')

  END SUBROUTINE MREAD_PARALLEL_2D

  SUBROUTINE MREAD_PARALLEL_3D (grd_dum,IUNIT,NAME,M,AVAR)
!@sum MREAD_PARALLEL  Parallel version of UTILDBL.f:MREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT       !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME        !@var  NAME  name of file being read
    INTEGER,      INTENT(OUT) :: M           !@var  M      initial integer
    REAL*8,      INTENT(OUT)  :: AVAR(:,grd_dum%J_STRT_HALO:,:) !@var  AOUT real*8 array
    REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3))  !@var  AIN  real*4 array
    REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var  AOUT real*8 array
    INTEGER :: N                        !@var  N loop variable
    INTEGER :: IM,JM,NM
    INTEGER :: rc 
    ! now local

    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD
    NM   = size(AVAR,3)

    If (AM_I_ROOT()) then
      READ (IUNIT,IOSTAT=rc) M, AIN
      !****  convert from real*4 to real*8
      AOUT=AIN
    EndIf

    Call scatterReal8(grd_dum, AOUT, AVAR, shape(AVAR), 2)
#ifdef USE_MPI
    CALL ESMF_BCAST(grd_dum, M   )
#endif
    call checkReadStatus(rc, name, 'MREAD_PARALLEL_3D')

  END SUBROUTINE MREAD_PARALLEL_3D

  SUBROUTINE READT_PARALLEL_2D (grd_dum,IUNIT,NAME,AVAR,IPOS)
!@sum READT_PARALLEL  Parallel version of UTILDBL.f:READT for (im,jm) arrays
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of file being read
    REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:)  !@var  AOUT real*8 array
    INTEGER,      INTENT(IN)  :: IPOS       !@var  IPOS  no. of recs. to advance
    REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD)  !@var  AIN  real*4 array
    REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var  AOUT real*8 array
    INTEGER :: N                        !@var  N loop variable
    INTEGER :: IM,JM
    CHARACTER*80 :: TITLE               !@var  TITLE title of file record
    INTEGER :: rc
    ! now local

    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD

    If (AM_I_ROOT()) then
      DO N=1,IPOS-1
        READ (IUNIT,IOSTAT=rc)
      END DO
      READ (IUNIT, IOSTAT=rc) TITLE, AIN
      !****  convert from real*4 to real*8
      AOUT=AIN
    EndIf

    call scatterReal8(grd_dum, AOUT, AVAR, shape(AVAR), 2)
    call checkReadStatus(rc, name, 'READT_PARALLEL_2D', title)

  END SUBROUTINE READT_PARALLEL_2D


  SUBROUTINE READT_PARALLEL_3D (grd_dum,IUNIT,NAME,AVAR,IPOS)
!@sum READT_PARALLEL  Parallel version of UTILDBL.f:READT for (im,jm) arrays
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT        !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME         !@var  NAME  name of file being read
    REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:,:)  !@var  AOUT real*8 array
    INTEGER,      INTENT(IN)  :: IPOS         !@var  IPOS  no. of recs. to advance
    REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3))  !@var  AIN  real*4 array
    REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var  AOUT real*8 array
    INTEGER :: N                        !@var  N loop variable
    INTEGER :: IM,JM,NM
    CHARACTER*80 :: TITLE               !@var  TITLE title of file record
    INTEGER :: rc
    ! now local
    IM   = grd_dum%IM_WORLD
    JM   = grd_dum%JM_WORLD
    NM   = size(AVAR,3)

    If (AM_I_ROOT()) then
      DO N=1,IPOS-1
        READ (IUNIT,IOSTAT=rc)
      END DO
      READ (IUNIT, IOSTAT=rc) TITLE, AIN
      !****  convert from real*4 to real*8
      AOUT=AIN
    EndIf

    Call scatterReal8(grd_dum, AOUT, AVAR, shape(AVAR), 2)
    call checkReadStatus(rc, name, 'READT_PARALLEL_3D', title)

  END SUBROUTINE READT_PARALLEL_3D

  subroutine READT_PARALLEL_COLUMN_3D (grd_dum,IUNIT,NAME,AVAR,IPOS)
!@sum READT_PARALLEL_COLUMN read in real*4 (:,im,jm) arrays
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT        !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME         !@var  NAME  name of file being read
    REAL*8,       INTENT(OUT) :: AVAR(:,:,grd_dum%J_STRT_HALO:)  !@var  AVAR real*8 array
    INTEGER,      INTENT(IN)  :: IPOS         !@var  IPOS  no. of recs. to advance
    REAL*4 :: AGLOB4(size(AVAR,1),grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var AGLOB global array
    REAL*8 :: AGLOB (size(AVAR,1),grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var AGLOB global array
    INTEGER :: N                        !@var  N loop variable
    CHARACTER*80 :: TITLE               !@var  TITLE title of file record
    INTEGER :: rc

    If (AM_I_ROOT()) then
      DO N=1,IPOS-1
        READ (IUNIT,IOSTAT=rc)
      END DO
      READ (IUNIT, IOSTAT=rc) TITLE, AGLOB4
      !****  convert from real*4 to real*8
      AGLOB = AGLOB4
    EndIf

    Call scatterReal8(grd_dum, AGLOB, AVAR, shape(AVAR), 3)
    call checkReadStatus(rc, name, 'READT_PARALLEL_COLUMN_3D', title)

  end subroutine READT_PARALLEL_COLUMN_3D

  subroutine READT8_PARALLEL_3D (grd_dum,IUNIT,NAME,AVAR,IPOS)
!@sum READT8_PARALLEL read in real*8 (im,jm,:) arrays
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT        !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME         !@var  NAME  name of file being read
    REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:,:)  !@var  AOUT real*8 array
    INTEGER,      INTENT(IN)  :: IPOS         !@var  IPOS  no. of recs. to advance
    REAL*8 :: AGLOB(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var AGLOB global array
    INTEGER :: N                        !@var  N loop variable
    CHARACTER*80 :: TITLE               !@var  TITLE title of file record
    INTEGER :: rc

    If (AM_I_ROOT()) then
      DO N=1,IPOS-1
        READ (IUNIT,IOSTAT=rc)
      END DO
      READ (IUNIT, IOSTAT=rc) TITLE, AGLOB
    EndIf

    Call scatterReal8(grd_dum, AGLOB, AVAR, shape(AVAR), 2)
    call checkReadStatus(rc, name, 'READT8_PARALLEL_3D', title)

  END SUBROUTINE READT8_PARALLEL_3D

  subroutine READT8_COLUMN_3D (grd_dum,IUNIT,NAME,AVAR,IPOS)
!@sum READT8_COLUMN read in real*8 (:,im,jm) arrays
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT        !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME         !@var  NAME  name of file being read
    REAL*8,       INTENT(OUT) :: AVAR(:,:,grd_dum%J_STRT_HALO:)  !@var  AVAR real*8 array
    INTEGER,      INTENT(IN)  :: IPOS         !@var  IPOS  no. of recs. to advance
    REAL*8 :: AGLOB(size(AVAR,1),grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var AGLOB global array
    INTEGER :: N                        !@var  N loop variable
    CHARACTER*80 :: TITLE               !@var  TITLE title of file record
    INTEGER :: rc

    If (AM_I_ROOT()) then
      DO N=1,IPOS-1
        READ (IUNIT,IOSTAT=rc)
      END DO
      READ (IUNIT, IOSTAT=rc) TITLE, AGLOB
    EndIf

    Call scatterReal8(grd_dum, AGLOB, AVAR, shape(AVAR), 3)
    call checkReadStatus(rc, name, 'READT8_COLUMN_3D', title)

  end subroutine READT8_COLUMN_3D

  SUBROUTINE WRITET8_COLUMN_3D (grd_dum,IUNIT,NAME,buf,title)
!@sum WRITET8_COLUMN write title*80, real*8 buf(:,im,jm)
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of file being read
    REAL*8,       INTENT(IN) :: buf(:,:,grd_dum%J_STRT_HALO:)  !@var  buf real*8 array
    CHARACTER*80,  INTENT(IN)  :: title
    REAL*8 :: buf_glob(size(buf,1),grd_dum%IM_WORLD,grd_dum%JM_WORLD)
    INTEGER :: rc

    Call gatherReal8(grd_dum, buf, buf_glob, shape(buf), 3, .false.)
    if (am_i_root()) write (IUNIT, IOSTAT=rc) title, buf_glob
    call checkWriteStatus(rc, name, 'WRITET8_COLUMN_3D', title)

  END SUBROUTINE WRITET8_COLUMN_3D

  SUBROUTINE WRITEI_PARALLEL_2D (grd_dum,IUNIT,NAME,buf,it)
!@sum WRITEI_PARALLEL  Parallel version of UTILDBL.f:WRITEI for (im,jm) arrays
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of record being read
    REAL*4,       INTENT(IN) :: buf(:,grd_dum%J_STRT_HALO:)  !@var  buf real*8 array
    INTEGER,      INTENT(IN)  :: it       !@var  it iteration
    REAL*8 :: buf8(grd_dum%IM_WORLD, &
         &     grd_dum%J_STRT_HALO:grd_dum%J_STOP_HALO)
!@var buf_glob real*4 array
    REAL*4 :: buf_glob(grd_dum%IM_WORLD,grd_dum%JM_WORLD)
    REAL*8 :: buf_glob8(grd_dum%IM_WORLD,grd_dum%JM_WORLD)
    INTEGER :: rc

! Use real*8 interface with copy to/from real*4
    buf8 = buf
    Call gatherReal8(grd_dum, buf8, buf_glob8, shape(buf), 2, .false.)
    if (am_i_root()) buf_glob = buf_glob8

    if (am_i_root()) write (iunit, iostat=rc) it, buf_glob, it
    call checkWriteStatus(rc, name, 'WRITEI_PARALLEL_2D')

  END SUBROUTINE WRITEI_PARALLEL_2D

  SUBROUTINE WRITET_PARALLEL_2D (grd_dum,IUNIT,NAME,buf,title)
!@sum WRITET_PARALLEL  write character*80 title, real(buf(1:im,1:jm),kind=4)
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of file being written
    REAL*8,       INTENT(IN) :: buf(:,grd_dum%J_STRT_HALO:)  !@var  buf real*8 array
    CHARACTER*80,  INTENT(IN)  :: title
    REAL*8 :: buf_glob(grd_dum%IM_WORLD,grd_dum%JM_WORLD)
    INTEGER :: rc

    Call gatherReal8(grd_dum, buf, buf_glob, shape(buf), 2, .false.)
    if (am_i_root()) write (iunit, iostat=rc) title, real(buf_glob,kind=4)
    call checkWriteStatus(rc, name, 'WRITET_PARALLEL_2D')

  END SUBROUTINE WRITET_PARALLEL_2D

  SUBROUTINE WRITEI8_PARALLEL_3D (grd_dum,IUNIT,NAME,buf,it)
!@sum WRITEI8_PARALLEL  Parallel version of UTILDBL.f:WRITEI8 for (im,jm,:) arrays
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of file being written
    REAL*8,       INTENT(IN) :: buf(:,grd_dum%J_STRT_HALO:,:)  !@var  buf real*8 array
    INTEGER,      INTENT(IN)  :: it       !@var  it iteration
    REAL*8 :: buf_glob(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(buf,3))
    INTEGER :: rc 

    Call gatherReal8(grd_dum, buf, buf_glob, shape(buf), 2, .false.)

    if (am_i_root()) write (iunit, iostat=rc) it, buf_glob, it
    call checkWriteStatus(rc, name, 'WRITEI8_PARALLEL_3D')

  END SUBROUTINE WRITEI8_PARALLEL_3D

  SUBROUTINE DWRITE8_PARALLEL_3D (grd_dum,IUNIT,NAME,buf)
!@sum DWRITE8_PARALLEL  Writes real*8 (im,jm,:) arrays as real*8 on disk
!@auth NCCS-ESMF Development Team
    TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
    INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
    CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of file being written
    REAL*8,       INTENT(IN) :: buf(:,grd_dum%J_STRT_HALO:,:)  !@var  buf real*8 array
    REAL*8, dimension (:,:,:), allocatable :: buf_glob ! global array written to disk
    INTEGER :: rc

    if (am_i_root()) then
      allocate(buf_glob(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(buf,3)))
    else
      allocate(buf_glob(1,1,1))
    end if

    if (am_i_root()) write (iunit, IOSTAT=rc) buf_glob
    call checkWriteStatus(rc, name, 'DWRITE8_PARALLEL_3D')
    deallocate(buf_glob)

  END SUBROUTINE DWRITE8_PARALLEL_3D
  
  subroutine READ_PARALLEL_INTEGER_0 ( DATA, UNIT, FORMAT)
    ! Wrapper for ESMF_READ_PARALLEL_INTEGER_0

    integer, intent(out  )                       :: DATA
    integer,            intent(in   ),  optional :: UNIT
    character(len=*),   intent(in   ),  optional :: FORMAT

    ! now local
    type (ESMF_DELayout)                         :: layout

    layout=ESMF_LAYOUT
    call ESMF_READ_PARALLEL_INTEGER_0( layout, DATA, UNIT, FORMAT)

  END SUBROUTINE READ_PARALLEL_INTEGER_0

  !---------------------------
  subroutine ESMF_READ_PARALLEL_INTEGER_0(layout,DATA,UNIT,FORMAT)

    type (ESMF_DELayout)                         :: layout
    integer, intent(out  )                       :: DATA
    integer,            intent(in   ),  optional :: UNIT
    character(len=*),   intent(in   ),  optional :: FORMAT

    character(len=ESMF_MAXSTR) :: FORMATTED
    character(LEN=ESMF_MAXSTR) :: FILENAME
    logical                :: IS_NAMED
    integer                :: IOSTAT
    integer                :: ierr

    if (amRootESMF(layout)) then
      if (present(UNIT)) then
        inquire(unit=UNIT, formatted=FORMATTED)
        if   (FORMATTED == "YES") then
          if (present(FORMAT)) then
            read(UNIT, FORMAT, IOSTAT=IOSTAT) DATA
          else
            read(UNIT, *, IOSTAT=IOSTAT) DATA
          end if
        elseif(FORMATTED == "NO") then
          read(UNIT,       IOSTAT=IOSTAT) DATA
        end if

        if (IOSTAT < 0) then
          inquire(unit=UNIT, NAMED=IS_NAMED, NAME=FILENAME)
          if (IS_NAMED) then
            print *, "Premature end of file ",FILENAME
            goto 9999
          else
            print *, "Premature end of file UNNAMED"
            goto 9999
          end if
        end if
      else
        if (present(FORMAT)) then; read(*, FORMAT ) DATA
        else; read          *, DATA
        end if
      end if
    end if

#ifdef USE_MPI
    CALL MPI_BCAST(data, 1, MPI_INTEGER, ROOT_PROCESS, &
         &       MPI_COMM_WORLD, ierr)
#endif

    return
9999 continue
    call stop_model('ESMF_READ_PARALLEL: READ ERROR',255)
  END SUBROUTINE ESMF_READ_PARALLEL_INTEGER_0

  !---------------------------

  subroutine READ_PARALLEL_INTEGER_1 ( DATA, UNIT, FORMAT)
    ! Wrapper for ESMF_READ_PARALLEL_INTEGER_1

    integer, intent(out  )                       :: DATA(:)
    integer,            intent(in   ),  optional :: UNIT
    character(len=*),   intent(in   ),  optional :: FORMAT

    ! now local
    type (ESMF_DELayout)                         :: layout

    layout = ESMF_LAYOUT
    call ESMF_READ_PARALLEL_INTEGER_1 ( layout, DATA, UNIT, FORMAT)

  END SUBROUTINE READ_PARALLEL_INTEGER_1

  !---------------------------
  subroutine ESMF_READ_PARALLEL_INTEGER_1(layout,DATA,UNIT,FORMAT)

    type (ESMF_DELayout)                         :: layout
    integer, intent(out  )                       :: DATA(:)
    integer,            intent(in   ),  optional :: UNIT
    character(len=*),   intent(in   ),  optional :: FORMAT

    character(len=ESMF_MAXSTR) :: FORMATTED
    character(LEN=ESMF_MAXSTR) :: FILENAME
    logical                :: IS_NAMED
    integer                :: IOSTAT
    integer                :: ierr

    if (amRootESMF(layout)) then
      if (present(UNIT)) then
        inquire(unit=UNIT, formatted=FORMATTED)
        if   (FORMATTED == "YES") then
          if (present(FORMAT)) then
            read(UNIT, FORMAT, IOSTAT=IOSTAT) DATA
          else
            read(UNIT, *, IOSTAT=IOSTAT) DATA
          end if
        elseif(FORMATTED == "NO") then
          read(UNIT,       IOSTAT=IOSTAT) DATA
        end if

        if (IOSTAT < 0) then
          inquire(unit=UNIT, NAMED=IS_NAMED, NAME=FILENAME)
          if (IS_NAMED) then
            print *, "Premature end of file ",FILENAME
            goto 9999
          else
            print *, "Premature end of file UNNAMED"
            goto 9999
          end if
        end if
      else
        if (present(FORMAT)) then; read(*, FORMAT ) DATA
        else; read          *, DATA
        end if
      end if
    end if

#ifdef USE_MPI
    CALL MPI_BCAST(data, size(data), MPI_INTEGER, ROOT_PROCESS, &
         &       MPI_COMM_WORLD, ierr)
#endif

    return
9999 continue
    call stop_model('ESMF_READ_PARALLEL: READ ERROR',255)
  END SUBROUTINE ESMF_READ_PARALLEL_INTEGER_1

  !---------------------------

  subroutine READ_PARALLEL_REAL8_1 ( DATA, UNIT, FORMAT)
    ! Wrapper for READ_PARALLEL_REAL8_1

    real(kind=8),       intent(out  )            :: DATA(:)
    integer,            intent(in   ),  optional :: UNIT
    character(len=*),   intent(in   ),  optional :: FORMAT

    ! now local
    type (ESMF_DELayout)                         :: layout

    layout = ESMF_LAYOUT
    call ESMF_READ_PARALLEL_REAL8_1 ( layout, DATA, UNIT, FORMAT)

  END SUBROUTINE READ_PARALLEL_REAL8_1

  !---------------------------

  subroutine ESMF_READ_PARALLEL_REAL8_1 (layout,DATA,UNIT,FORMAT)

    type (ESMF_DELayout)                         :: layout
    real(kind=8),       intent(out  )            :: DATA(:)
    integer,            intent(in   ),  optional :: UNIT
    character(len=*),   intent(in   ),  optional :: FORMAT

    character(len=ESMF_MAXSTR) :: FORMATTED
    character(LEN=ESMF_MAXSTR) :: FILENAME
    logical                :: IS_NAMED
    integer                :: IOSTAT
    integer                :: ierr

    if (amRootESMF(layout)) then
      if (present(UNIT)) then
        inquire(unit=UNIT, formatted=FORMATTED)
        if   (FORMATTED == "YES") then
          if (present(FORMAT)) then
            read(UNIT, FORMAT, IOSTAT=IOSTAT) DATA
          else
            read(UNIT, *, IOSTAT=IOSTAT) DATA
          end if
        elseif(FORMATTED == "NO") then
          read(UNIT,       IOSTAT=IOSTAT) DATA
        end if

        if (IOSTAT < 0) then
          inquire(unit=UNIT, NAMED=IS_NAMED, NAME=FILENAME)
          if (IS_NAMED) then
            print *, "Premature end of file ",FILENAME
            goto 9999
          else
            print *, "Premature end of file UNNAMED"
            goto 9999
          end if
        end if
      else
        if (present(FORMAT)) then; read(*, FORMAT ) DATA
        else; read          *, DATA
        end if
      end if
    end if

#ifdef USE_MPI
    CALL MPI_BCAST(data, Size(data), MPI_DOUBLE_PRECISION, ROOT_PROCESS, &
         &       MPI_COMM_WORLD, ierr)
#endif

    return
9999 continue
    call stop_model('ESMF_READ_PARALLEL: READ ERROR',255)
  END SUBROUTINE ESMF_READ_PARALLEL_REAL8_1

  subroutine WRITE_PARALLEL_INTEGER_0 ( data, UNIT, format, CRIT)

    INTEGER, intent(in   )            :: data
    integer,            intent(in   ),  optional :: UNIT
    character(len=*),   intent(in   ),  optional :: format
    logical,            intent(in   ), optional  :: CRIT

    character(len=ESMF_MAXSTR) :: FORMATTED

    logical :: crit_,crit2_ ! local

    crit_ = .false.
    crit2_= .false.
    if (present(CRIT)) then
      crit_ = crit
    else
      crit2_ = AM_I_ROOT()
    endif

    if (crit_ .or. crit2_) then
      if (present(UNIT)) then
        inquire(unit=UNIT, formatted=FORMATTED)
        if   (FORMATTED == "YES") then
          if (present(format)) then;    write(UNIT, format) data
          else;    write(UNIT, *     ) data
          end if
        elseif(FORMATTED == "NO") then; write(UNIT        ) data
        end if
      else
        if (present(format)) then; write(*, format) data
        else; write(*,      *) data
        end if
      end if
    end if

    return

  end subroutine WRITE_PARALLEL_INTEGER_0

  !---------------------------
  subroutine WRITE_PARALLEL_INTEGER_1 ( data, UNIT, format, CRIT)

    INTEGER, intent(in   )            :: data (:)
    integer,            intent(in   ),  optional :: UNIT
    character(len=*),   intent(in   ),  optional :: format
    logical,            intent(in   ),  optional :: CRIT

    character(len=ESMF_MAXSTR) :: FORMATTED
    logical :: crit_,crit2_ ! local

    crit_ = .false.
    crit2_= .false.
    if (present(CRIT)) then
      crit_ = crit
    else
      crit2_ = AM_I_ROOT()
    endif

    if (crit_ .or. crit2_) then
      if (present(UNIT)) then
        inquire(unit=UNIT, formatted=FORMATTED)
        if   (FORMATTED == "YES") then
          if (present(format)) then;    write(UNIT, format) data
          else;    write(UNIT, *     ) data
          end if
        elseif(FORMATTED == "NO") then; write(UNIT        ) data
        end if
      else
        if (present(format)) then; write(*, format) data
        else; write(*,      *) data
        end if
      end if
    end if

    return

  end subroutine WRITE_PARALLEL_INTEGER_1

  !---------------------------

  subroutine WRITE_PARALLEL_REAL8_0 ( data, UNIT, format, CRIT)

    REAL (KIND=8), intent(in   )            :: data
    integer,            intent(in   ),  optional :: UNIT
    character(len=*),   intent(in   ),  optional :: format
    logical,            intent(in   ),  optional :: CRIT

    character(len=ESMF_MAXSTR) :: FORMATTED

    logical :: crit_,crit2_ ! local

    crit_ = .false.
    crit2_= .false.
    if (present(CRIT)) then
      crit_ = crit
    else
      crit2_ = AM_I_ROOT()
    endif

    if (crit_ .or. crit2_) then
      if (present(UNIT)) then
        inquire(unit=UNIT, formatted=FORMATTED)
        if   (FORMATTED == "YES") then
          if (present(format)) then;    write(UNIT, format) data
          else;    write(UNIT, *     ) data
          end if
        elseif(FORMATTED == "NO") then; write(UNIT        ) data
        end if
      else
        if (present(format)) then; write(*, format) data
        else; write(*,      *) data
        end if
      end if
    end if

    return

  end subroutine WRITE_PARALLEL_REAL8_0

  !---------------------------
  subroutine WRITE_PARALLEL_REAL8_1 ( data, UNIT, format, CRIT)

    REAL (KIND=8), intent(in   )            :: data(:)
    integer,            intent(in   ),  optional :: UNIT
    character(len=*),   intent(in   ),  optional :: format
    logical,            intent(in   ), optional  :: CRIT

    character(len=ESMF_MAXSTR) :: FORMATTED

    logical :: crit_,crit2_ ! local

    crit_ = .false.
    crit2_= .false.
    if (present(CRIT)) then
      crit_ = crit
    else
      crit2_ = AM_I_ROOT()
    endif

    if (crit_ .or. crit2_) then
      if (present(UNIT)) then
        inquire(unit=UNIT, formatted=FORMATTED)
        if   (FORMATTED == "YES") then
          if (present(format)) then;    write(UNIT, format) data
          else;    write(UNIT, *     ) data
          end if
        elseif(FORMATTED == "NO") then; write(UNIT        ) data
        end if
      else
        if (present(format)) then; write(*, format) data
        else; write(*,      *) data
        end if
      end if
    end if

    return

  end subroutine WRITE_PARALLEL_REAL8_1

  !---------------------------

  subroutine WRITE_PARALLEL_STRING_0 (STRING, FORMAT, UNIT, CRIT)
    character(LEN=*), intent(in   )            :: STRING
    character(LEN=*), intent(in   ), optional  :: FORMAT
    integer,          intent(in   ),  optional :: UNIT
    logical,          intent(in   ), optional  :: CRIT
    logical :: crit_,crit2_ ! local

    crit_ = .false.
    crit2_= .false.
    if (present(CRIT)) then
      crit_ = crit
    else
      crit2_ = AM_I_ROOT()
    endif

    if (crit_ .or. crit2_) then
      if (present(UNIT)) then
        if (present(FORMAT)) then
          write(UNIT, FMT=FORMAT) STRING
        else
          write(UNIT, FMT='(A)') STRING
        endif
      else
        if (present(FORMAT)) then
          write(*, FMT=FORMAT) STRING
        else
          write(*, FMT='(A)') STRING
        endif
      end if
    end if

  end subroutine WRITE_PARALLEL_STRING_0

  !---------------------------

  subroutine WRITE_PARALLEL_STRING_1 (STRING, FORMAT, UNIT, CRIT)
    character(LEN=*), intent(in   )            :: STRING(:)
    character(LEN=*), intent(in   ), optional  :: FORMAT
    integer,          intent(in   ), optional :: UNIT
    logical,          intent(in   ), optional  :: CRIT
    logical :: crit_,crit2_ ! local

    crit_ = .false.
    crit2_= .false.
    if (present(CRIT)) then
      crit_ = crit
    else
      crit2_ = AM_I_ROOT()
    endif

    if (crit_ .or. crit2_) then
      if (present(UNIT)) then
        if (present(FORMAT)) then
          write(UNIT, FMT=FORMAT) STRING
        else
          write(UNIT, FMT='(A)') STRING
        endif
      else
        if (present(FORMAT)) then
          write(*, FMT=FORMAT) STRING
        else
          write(*, FMT='(A)') STRING
        endif
      end if
    end if

  end subroutine WRITE_PARALLEL_STRING_1

  !---------------------------

  subroutine checkReadStatus(rc, fileName, procedureName, title)
    integer, intent(in) :: rc
    character(len=*), intent(in) :: fileName
    character(len=*), intent(in) :: procedureName
    character(len=*), optional, intent(in) :: title

    if (am_i_root()) then
      if (rc==0) then
        write(6,'(a,a)',advance='no') "Read from file ",trim(fileName)
        if (present(title)) then
          write(6,*)': ', trim(title)
        else
          write(6,*)
        end if
      else
        write(6,'(a,a)',advance='no') 'READ ERROR ON FILE ',trim(fileName)
        if (present(title)) write(6,'(a,a)',advance='no') ': ',  trim(title)
        write(6,*)': IOSTAT=',rc
        call stop_model(trim(procedureName)//': READ ERROR',255)
      endif
    end if

  end subroutine checkReadStatus

  subroutine checkWriteStatus(rc, fileName, procedureName, title)
    integer, intent(in) :: rc
    character(len=*), intent(in) :: fileName
    character(len=*), intent(in) :: procedureName
    character(len=*), optional, intent(in) :: title

    if (am_i_root()) then
      if (rc==0) then
        write(6,'(a,a)',advance='no') "Wrote to file ",trim(fileName)
        if (present(title)) then
          write(6,*)': ', trim(title)
        else
          write(6,*)
        end if
      else
        write(6,'(a,a)',advance='no') 'WRITE ERROR ON FILE ',trim(fileName)
        if (present(title)) write(6,'(a,a)',advance='no')': ', trim(title)
        write(6,*)': IOSTAT=',rc
        call stop_model(trim(procedureName)//': WRITE ERROR',255)
      endif
    end if

  end subroutine checkWriteStatus

end module SpecialIO_mod
