
!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function AppGridCreateF(cf, vm, rc) result(grid)
#include "MAPL_Generic.h"

  use ESMF_Mod
  use MAPL_BaseMod
  use MAPL_GenericMod
  use MAPL_ConstantsMod, only : pi=> MAPL_PI
  use fv_grid_utils_mod, only: gnomonic_grids, cell_center2
  use fv_grid_tools_mod, only: mirror_grid
  implicit none

! !ARGUMENTS:
    type(ESMF_Config), intent(INOUT) :: cf
    type (ESMF_VM),    intent(IN   ) :: VM
    integer, optional, intent(OUT)   :: rc
    type (ESMF_Grid)                 :: grid

! ErrLog variables
!-----------------

 integer                      :: STATUS
 character(len=ESMF_MAXSTR), parameter :: Iam="AppGridCreateF"

! Local variables
!-----------------

#ifdef EIGHT_BYTE
  integer, parameter:: f_p = selected_real_kind(15)   ! same as 12 on Altix
#else
! Higher precisions for grid geometrical factors:
  integer, parameter:: f_p = selected_real_kind(20)
#endif

  integer, parameter            :: grid_type = 0
  type(ESMF_DELayout)           :: layout
  integer                       :: DECOUNT(2)
  integer                       :: GLOBAL_DE_START(2)
  integer                       :: MYID
  integer                       :: NPES, NPES_X, NPES_Y
  integer                       :: NPX, NPY
  integer                       :: NX, NY
  integer                       :: isg, ieg
  integer                       :: jsg, jeg
  integer                       :: is, ie
  integer                       :: js, je
  integer                       :: myTile
  integer                       :: npts
  integer                       :: ntiles=6
  integer                       :: ndims=2
  integer                       :: I, J, N
  integer                       :: IG, JG

  real(ESMF_KIND_R8), allocatable :: xs(:,:)
  real(ESMF_KIND_R8), allocatable :: ys(:,:)
  real(ESMF_KIND_R8), allocatable :: grid_global(:,:,:,:)
  type (ESMF_Array),  pointer     :: coords(:)
  real(ESMF_KIND_R8), pointer     :: lons(:,:)
  real(ESMF_KIND_R8), pointer     :: lats(:,:)
  type (ESMF_Array),  target      :: tarray(2)
  real(ESMF_KIND_R8)              :: alocs(2)

  integer                         :: IM_WORLD
  integer                         :: JM_WORLD
  integer                         :: LM
  integer                         :: L
  integer, allocatable            :: IMS(:), JMS(:)
  character(len=ESMF_MAXSTR)      :: gridname, FMT, FMTIM, FMTJM
  real(ESMF_KIND_R8)              :: minCoord(3)
  real(ESMF_KIND_R8)              :: deltaX, deltaY, deltaZ

! grid create

    call ESMF_ConfigGetAttribute( cf, IM_WORLD, label ='IM:', rc = status )
    VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute( cf, JM_WORLD, label ='JM:', rc = status )
    VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute( cf, LM,       label ='LM:', rc = status )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( cf, NX,       label ='NX:', rc = status )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( cf, NY,       label ='NY:', rc = status )
    VERIFY_(STATUS)

! Give the IMS and JMS the MAPL default distribution
! --------------------------------------------------

       allocate( IMS(0:NX-1) )
       allocate( JMS(0:NY-1) )

       call DecomposeDim ( IM_WORLD  , IMS             , NX   )
       call DecomposeDim ( JM_WORLD/6, JMS(0:NY/6 -1)  , NY/6 )
       do n=2,6
          JMS((n-1)*NY/6 : n*NY/6 -1) = JMS(0:NY/6 -1)
       enddo

! We should have a much simpler create with the next ESMF grid design
!--------------------------------------------------------------------

       layout = ESMF_DELayoutCreate(vm, deCountList=(/NX, NY/), rc=status)
       VERIFY_(STATUS)

       deltaX      = 2.0*PI/IM_WORLD
       deltaY = PI/(JM_WORLD  )
       minCoord(1) = 0.0
       minCoord(2) = -PI/2

       deltaZ = 1.0D0
       minCoord(3) = deltaZ/2

       call GET_INT_FORMAT(IM_WORLD, FMTIM)
       call GET_INT_FORMAT(JM_WORLD, FMTJM)
       FMT = '(A,' // trim(FMTIM) //',A,' // trim(FMTJM) // ',A)'
       write(gridname,trim(FMT)) 'PE',IM_WORLD,'x',JM_WORLD,'-CF'

       grid = ESMF_GridCreateHorzLatLonUni(         &
            counts = (/IM_WORLD, JM_WORLD/),        &
            minGlobalCoordPerDim=minCoord(1:2),     &
            deltaPerDim=(/deltaX, deltaY /),        &
            horzStagger=ESMF_Grid_Horz_Stagger_A,   &
            periodic=(/ESMF_FALSE, ESMF_FALSE/),     &
            name=gridname,                          &
            rc=status)
       VERIFY_(STATUS)

       if (LM>1) then
         call ESMF_GridAddVertHeight(grid,            &
              delta=(/(deltaZ, L=1,LM) /),            &
              vertStagger=ESMF_GRID_VERT_STAGGER_CENTER, &
              rc=status)
         VERIFY_(STATUS)
       else
         call ESMF_GridAddVertHeight(grid,            &
              delta=(/(deltaZ, L=1,2) /),            &
              vertStagger=ESMF_GRID_VERT_STAGGER_CENTER, &
              rc=status)
         VERIFY_(STATUS)
       endif

       call ESMF_GridDistribute(grid,               &
            deLayout=layout,                        &
            countsPerDEDim1=ims,                    &
            countsPerDEDim2=jms,                    &
            rc=status)
       VERIFY_(STATUS)

! -------------


 call ESMF_GridGet(GRID, deLayout=layout, RC=STATUS)
 VERIFY_(STATUS)

 call ESMF_DELayoutGet(layout, deCountPerDim=DECOUNT, localDE=MYID, RC=STATUS)
 VERIFY_(STATUS)
 NPES_X = DECOUNT(1)
 NPES_Y = DECOUNT(2)
 NPES = NPES_X+NPES_Y

 call ESMF_GridGetDELocalInfo(GRID, horzRelLoc=ESMF_CELL_CENTER, &
      localCellCountPerDim=DECOUNT, globalStartPerDim=GLOBAL_DE_START,RC=STATUS)
 VERIFY_(STATUS)
 NX = DECOUNT(1)
 NY = DECOUNT(2)

 npx = IM_WORLD
 npy = JM_WORLD
 isg = GLOBAL_DE_START(1) + 1
 ieg = GLOBAL_DE_START(1) + NX
 jsg = GLOBAL_DE_START(2) + 1
 jeg = GLOBAL_DE_START(2) + NY
 myTile = jsg/(npy/ntiles)

 is = isg
 ie = ieg
 js = jsg - myTile*(npy/ntiles)
 je = jeg - myTile*(npy/ntiles)

 npts = (npy/ntiles)
 if (npts /= npx) then
    print*, 'Error npts /= npx', npts, npx
    STATUS=1
 endif
 VERIFY_(STATUS)

!print*, MYID, myTile, is, ie, js, je

 allocate( xs(npts+1,npts+1) )
 allocate( ys(npts+1,npts+1) )
 allocate( grid_global(npts+1,npts+1,ndims,ntiles) )

  call gnomonic_grids(grid_type, npts, xs, ys)
  do j=1,npts+1
     do i=1,npts+1
        grid_global(i,j,1,1) = xs(i,j)
        grid_global(i,j,2,1) = ys(i,j)
     enddo
  enddo
! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi]
  call mirror_grid(grid_global, 0, npts+1, npts+1, 2, 6)

  do n=1,ntiles
     do j=1,npts+1
        do i=1,npts+1
!---------------------------------
! Shift the corner away from Japan
!---------------------------------
#ifdef SHIFT_WEST
! This will result in the corner close to east coast of China
           grid_global(i,j,1,n) = grid_global(i,j,1,n) - pi/18.
#endif
           if ( grid_global(i,j,1,n) < 0. )              &
                grid_global(i,j,1,n) = grid_global(i,j,1,n) + 2.*pi
           if (ABS(grid_global(i,j,1,1)) < 1.e-10) grid_global(i,j,1,1) = 0.0
           if (ABS(grid_global(i,j,2,1)) < 1.e-10) grid_global(i,j,2,1) = 0.0
        enddo
     enddo
  enddo
!---------------------------------
! Clean Up Corners
!---------------------------------
  grid_global(  1,1:npts+1,:,2)=grid_global(npts+1,1:npts+1,:,1)
  grid_global(  1,1:npts+1,:,3)=grid_global(npts+1:1:-1,npts+1,:,1)
  grid_global(1:npts+1,npts+1,:,5)=grid_global(1,npts+1:1:-1,:,1)
  grid_global(1:npts+1,npts+1,:,6)=grid_global(1:npts+1,1,:,1)
  grid_global(1:npts+1,  1,:,3)=grid_global(1:npts+1,npts+1,:,2)
  grid_global(1:npts+1,  1,:,4)=grid_global(npts+1,npts+1:1:-1,:,2)
  grid_global(npts+1,1:npts+1,:,6)=grid_global(npts+1:1:-1,1,:,2)
  grid_global(  1,1:npts+1,:,4)=grid_global(npts+1,1:npts+1,:,3)
  grid_global(  1,1:npts+1,:,5)=grid_global(npts+1:1:-1,npts+1,:,3)
  grid_global(npts+1,1:npts+1,:,3)=grid_global(1,1:npts+1,:,4)
  grid_global(1:npts+1,  1,:,5)=grid_global(1:npts+1,npts+1,:,4)
  grid_global(1:npts+1,  1,:,6)=grid_global(npts+1,npts+1:1:-1,:,4)
  grid_global(  1,1:npts+1,:,6)=grid_global(npts+1,1:npts+1,:,5)

! Fill lat/lons at cell center locations
  coords=>tarray
  call ESMF_GridGetCoord(GRID, horzRelloc=ESMF_CELL_CENTER, centerCoord=coords, rc=status)
  VERIFY_(STATUS)
! Cell Center Locations
 ! Longitudes
  call ESMF_ArrayGetData(coords(1), lons, RC=status)
  VERIFY_(STATUS)
 ! Latitudes
  call ESMF_ArrayGetData(coords(2), lats, RC=status)
  VERIFY_(STATUS)
  do jg=jsg,jeg
     do ig=isg,ieg
        i=ig
        j=jg-myTile*npts
        call cell_center2(grid_global(i,j,  1:2,myTile+1), grid_global(i+1,j,  1:2,myTile+1),   &
                          grid_global(i,j+1,1:2,myTile+1), grid_global(i+1,j+1,1:2,myTile+1),   &
                          alocs)
        i=ig-isg+1
        j=jg-jsg+1
        lons(i,j) = alocs(1)
        lats(i,j) = alocs(2)
     enddo
  enddo

  deallocate( xs )
  deallocate( ys )
  deallocate( grid_global )
  deallocate( IMS )
  deallocate( JMS )

    RETURN_(STATUS)
  end function AppGridCreateF

subroutine AppGridCreate (META, GRID, RC)

#include "MAPL_Generic.h"

  use ESMF_Mod
  use MAPL_BaseMod
  use MAPL_GenericMod
  use MAPL_ConstantsMod, only : pi=> MAPL_PI
  use fv_grid_utils_mod, only: gnomonic_grids, cell_center2
  use fv_grid_tools_mod, only: mirror_grid
  implicit none

! !ARGUMENTS:

 type(MAPL_MetaComp), intent(INOUT) :: META
 type (ESMF_Grid),    intent(  OUT) :: grid
 integer, optional,   intent(  OUT) :: rc

! ErrLog variables
!-----------------

 integer                      :: STATUS
 character(len=ESMF_MAXSTR), parameter :: Iam="AppGridCreate"

! Local variables
!-----------------

#ifdef EIGHT_BYTE
  integer, parameter:: f_p = selected_real_kind(15)   ! same as 12 on Altix
#else
! Higher precisions for grid geometrical factors:
  integer, parameter:: f_p = selected_real_kind(20)
#endif

  integer, parameter            :: grid_type = 0
  type(ESMF_DELayout)           :: layout
  integer                       :: DECOUNT(2)
  integer                       :: GLOBAL_DE_START(2)
  integer                       :: MYID
  integer                       :: NPES, NPES_X, NPES_Y
  integer                       :: NPX, NPY
  integer                       :: NX, NY 
  integer                       :: isg, ieg
  integer                       :: jsg, jeg
  integer                       :: is, ie
  integer                       :: js, je
  integer                       :: myTile
  integer                       :: npts
  integer                       :: ntiles=6
  integer                       :: ndims=2
  integer                       :: I, J, N
  integer                       :: IG, JG

  real(ESMF_KIND_R8), allocatable :: xs(:,:)
  real(ESMF_KIND_R8), allocatable :: ys(:,:)
  real(ESMF_KIND_R8), allocatable :: grid_global(:,:,:,:)
  type (ESMF_Array),  pointer     :: coords(:)
  real(ESMF_KIND_R8), pointer     :: lons(:,:)
  real(ESMF_KIND_R8), pointer     :: lats(:,:)
  type (ESMF_Array),  target      :: tarray(2)
  real(ESMF_KIND_R8)              :: alocs(2)

  integer                         :: IM_WORLD
  integer                         :: JM_WORLD
  integer                         :: LM
  integer                         :: L
  integer, allocatable            :: IMS(:), JMS(:)
  character(len=ESMF_MAXSTR)      :: gridname
  real(ESMF_KIND_R8)              :: minCoord(3)
  real(ESMF_KIND_R8)              :: deltaX, deltaY, deltaZ
  type (ESMF_VM)                  :: VM

! ------------


! We need the VM to create the grid ??
!------------------------------------

       call ESMF_VMGetCurrent(vm, rc=STATUS)
       VERIFY_(STATUS)

! Get Decomposition from CF
!--------------------------

! !RESOURCE_ITEM: none :: Processing elements in 1st dimension
       call MAPL_GetResource( META, NX,       label ='NX:', default=1, rc = status )
       VERIFY_(STATUS)
! !RESOURCE_ITEM: none :: Processing elements in 2nd dimension
       call MAPL_GetResource( META, NY,       label ='NY:', default=1, rc = status )
       VERIFY_(STATUS)

! Get World problem size from CF
!-------------------------------

! !RESOURCE_ITEM: none :: Grid size in 1st dimension
       call MAPL_GetResource( META, IM_WORLD, 'IM:',            rc = status )
       VERIFY_(STATUS)
! !RESOURCE_ITEM: none :: Grid size in 2nd dimension
       call MAPL_GetResource( META, JM_WORLD, 'JM:',            rc = status )
       VERIFY_(STATUS)
       JM_WORLD = 6*JM_WORLD
! !RESOURCE_ITEM: none :: Grid size in 3rd dimension
       call MAPL_GetResource( META, LM,       'LM:', default=1, rc = status )
       VERIFY_(STATUS)

! The grid's name is optional
!----------------------------

! !RESOURCE_ITEM: none :: Optional grid name
       call MAPL_GetResource( META, GRIDNAME, 'GRIDNAME:', default='APPGRID', rc = status )
       VERIFY_(STATUS)

! Give the IMS and JMS the MAPL default distribution
! --------------------------------------------------

       allocate( IMS(0:NX-1) )
       allocate( JMS(0:NY-1) )

       call DecomposeDim ( IM_WORLD  , IMS             , NX   )
       call DecomposeDim ( JM_WORLD/6, JMS(0:NY/6 -1)  , NY/6 )
       do n=2,6
          JMS((n-1)*NY/6 : n*NY/6 -1) = JMS(0:NY/6 -1)
       enddo

! Override them with alues in CF, if any
!---------------------------------------

!! !RESOURCE_ITEM: none :: gridpoints in each PE along 1st dimension
!       call MAPL_GetResource( META, IMS, 'IMS:', default=IMS, rc = status )
!       VERIFY_(STATUS)
!! !RESOURCE_ITEM: none :: gridpoints in each PE along 2nd dimension
!       call MAPL_GetResource( META, JMS, 'JMS:', default=JMS, rc = status )
!       VERIFY_(STATUS)

! We should have a much simpler create with the next ESMF grid design
!--------------------------------------------------------------------

       layout = ESMF_DELayoutCreate(vm, deCountList=(/NX, NY/), rc=status)
       VERIFY_(STATUS)

       deltaX      = 2.0*PI/IM_WORLD
       deltaY = PI/(JM_WORLD  )
       minCoord(1) = 0.0
       minCoord(2) = -PI/2

       deltaZ = 1.0D0
       minCoord(3) = deltaZ/2

       write(gridname,201) 'PE',IM_WORLD,'x',JM_WORLD,'-CF'
 201   format(A,i4.4,A,i4.4,A)
!      print*, gridname

       grid = ESMF_GridCreateHorzLatLonUni(         &
            counts = (/IM_WORLD, JM_WORLD/),        &
            minGlobalCoordPerDim=minCoord(1:2),     &
            deltaPerDim=(/deltaX, deltaY /),        &
            horzStagger=ESMF_Grid_Horz_Stagger_A,   &
            periodic=(/ESMF_FALSE, ESMF_FALSE/),     &
            name=gridname,                          &
            rc=status)
       VERIFY_(STATUS)

       if (LM>1) then
         call ESMF_GridAddVertHeight(grid,            &
              delta=(/(deltaZ, L=1,LM) /),            &
              vertStagger=ESMF_GRID_VERT_STAGGER_CENTER, &
              rc=status)
         VERIFY_(STATUS)
       else
         call ESMF_GridAddVertHeight(grid,            &
              delta=(/(deltaZ, L=1,2) /),            &
              vertStagger=ESMF_GRID_VERT_STAGGER_CENTER, &
              rc=status)
         VERIFY_(STATUS)
       endif

       call ESMF_GridDistribute(grid,               &
            deLayout=layout,                        &
            countsPerDEDim1=ims,                    &
            countsPerDEDim2=jms,                    &
            rc=status)
       VERIFY_(STATUS)

! -------------


 call ESMF_GridGet(GRID, deLayout=layout, RC=STATUS)
 VERIFY_(STATUS)

 call ESMF_DELayoutGet(layout, deCountPerDim=DECOUNT, localDE=MYID, RC=STATUS)
 VERIFY_(STATUS)
 NPES_X = DECOUNT(1)
 NPES_Y = DECOUNT(2)
 NPES = NPES_X+NPES_Y

 call ESMF_GridGetDELocalInfo(GRID, horzRelLoc=ESMF_CELL_CENTER, &
      localCellCountPerDim=DECOUNT, globalStartPerDim=GLOBAL_DE_START,RC=STATUS)
 VERIFY_(STATUS)
 NX = DECOUNT(1)
 NY = DECOUNT(2)

 npx = IM_WORLD
 npy = JM_WORLD
 isg = GLOBAL_DE_START(1) + 1
 ieg = GLOBAL_DE_START(1) + NX
 jsg = GLOBAL_DE_START(2) + 1
 jeg = GLOBAL_DE_START(2) + NY
 myTile = jsg/(npy/ntiles) 

 is = isg 
 ie = ieg
 js = jsg - myTile*(npy/ntiles)
 je = jeg - myTile*(npy/ntiles)

 npts = (npy/ntiles)
 if (npts /= npx) then 
    print*, 'Error npts /= npx', npts, npx
    STATUS=1
 endif
 VERIFY_(STATUS)

!print*, MYID, myTile, ie-is+1, je-js+1

 allocate( xs(npts+1,npts+1) )
 allocate( ys(npts+1,npts+1) )
 allocate( grid_global(npts+1,npts+1,ndims,ntiles) )

  call gnomonic_grids(grid_type, npts, xs, ys)
  do j=1,npts+1
     do i=1,npts+1
        grid_global(i,j,1,1) = xs(i,j)
        grid_global(i,j,2,1) = ys(i,j)
     enddo
  enddo
! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi] 
  call mirror_grid(grid_global, 0, npts+1, npts+1, 2, 6)

  do n=1,ntiles
     do j=1,npts+1
        do i=1,npts+1
!---------------------------------
! Shift the corner away from Japan
!---------------------------------
#ifdef SHIFT_WEST
! This will result in the corner close to east coast of China
           grid_global(i,j,1,n) = grid_global(i,j,1,n) - pi/18.
#endif
           if ( grid_global(i,j,1,n) < 0. )              &
                grid_global(i,j,1,n) = grid_global(i,j,1,n) + 2.*pi
           if (ABS(grid_global(i,j,1,1)) < 1.e-10) grid_global(i,j,1,1) = 0.0
           if (ABS(grid_global(i,j,2,1)) < 1.e-10) grid_global(i,j,2,1) = 0.0
        enddo
     enddo
  enddo
!---------------------------------
! Clean Up Corners
!---------------------------------
  grid_global(  1,1:npts+1,:,2)=grid_global(npts+1,1:npts+1,:,1)
  grid_global(  1,1:npts+1,:,3)=grid_global(npts+1:1:-1,npts+1,:,1)
  grid_global(1:npts+1,npts+1,:,5)=grid_global(1,npts+1:1:-1,:,1)
  grid_global(1:npts+1,npts+1,:,6)=grid_global(1:npts+1,1,:,1)
  grid_global(1:npts+1,  1,:,3)=grid_global(1:npts+1,npts+1,:,2)
  grid_global(1:npts+1,  1,:,4)=grid_global(npts+1,npts+1:1:-1,:,2)
  grid_global(npts+1,1:npts+1,:,6)=grid_global(npts+1:1:-1,1,:,2)
  grid_global(  1,1:npts+1,:,4)=grid_global(npts+1,1:npts+1,:,3)
  grid_global(  1,1:npts+1,:,5)=grid_global(npts+1:1:-1,npts+1,:,3)
  grid_global(npts+1,1:npts+1,:,3)=grid_global(1,1:npts+1,:,4)
  grid_global(1:npts+1,  1,:,5)=grid_global(1:npts+1,npts+1,:,4)
  grid_global(1:npts+1,  1,:,6)=grid_global(npts+1,npts+1:1:-1,:,4)
  grid_global(  1,1:npts+1,:,6)=grid_global(npts+1,1:npts+1,:,5)

! Fill lat/lons at cell center locations
  coords=>tarray
  call ESMF_GridGetCoord(GRID, horzRelloc=ESMF_CELL_CENTER, centerCoord=coords, rc=status)
  VERIFY_(STATUS)
! Cell Center Locations
 ! Longitudes
  call ESMF_ArrayGetData(coords(1), lons, RC=status)
  VERIFY_(STATUS)
 ! Latitudes
  call ESMF_ArrayGetData(coords(2), lats, RC=status)
  VERIFY_(STATUS)
  do jg=jsg,jeg
     do ig=isg,ieg
        i=ig
        j=jg-myTile*npts
        call cell_center2(grid_global(i,j,  1:2,myTile+1), grid_global(i+1,j,  1:2,myTile+1),   &
                          grid_global(i,j+1,1:2,myTile+1), grid_global(i+1,j+1,1:2,myTile+1),   &
                          alocs)
        i=ig-isg+1
        j=jg-jsg+1
        lons(i,j) = alocs(1)
        lats(i,j) = alocs(2)
     enddo
  enddo

  deallocate( xs )
  deallocate( ys )
  deallocate( grid_global )
  deallocate( IMS )
  deallocate( JMS )

 end subroutine AppGridCreate

  subroutine GET_INT_FORMAT(N, FMT)
    integer          :: N
    character(len=*) :: FMT

    IF(N < 10) THEN
       FMT = 'I1'
    ELSE IF (N< 100) THEN
       FMT = 'I2'
    ELSE IF (N< 1000) THEN
       FMT = 'I3'
    ELSE IF (N< 10000) THEN
       FMT = 'I4'
    else
       FMT = 'I5'
    end IF
  end subroutine GET_INT_FORMAT

  subroutine DecomposeDim ( dim_world,dim,NDEs )

!
! From FMS/MPP
!

      implicit   none
      integer    dim_world, NDEs
      integer    dim(0:NDEs-1)

      integer :: is,ie,isg,ieg
      integer :: ndiv,ndivs,imax,ndmax,ndmirror,n
      integer :: ibegin(0:NDEs-1)
      integer :: iend(0:NDEs-1)

      logical :: symmetrize
      logical :: even, odd
      even(n) = (mod(n,2).EQ.0)
      odd (n) = (mod(n,2).EQ.1)

       isg = 1
       ieg = dim_world
       ndivs = NDEs

       is = isg
       n = 0
       do ndiv=0,ndivs-1
             !modified for mirror-symmetry
             !original line
             !                 ie = is + CEILING( float(ieg-is+1)/(ndivs-ndiv) ) - 1

             !problem of dividing nx points into n domains maintaining symmetry
             !i.e nx=18 n=4 4554 and 5445 are solutions but 4455 is not.
             !this will always work for nx even n even or odd
             !this will always work for nx odd, n odd
             !this will never  work for nx odd, n even: for this case we supersede the mirror calculation
             !                 symmetrize = .NOT. ( mod(ndivs,2).EQ.0 .AND. mod(ieg-isg+1,2).EQ.1 )
             !nx even n odd fails if n>nx/2
             symmetrize = ( even(ndivs) .AND. even(ieg-isg+1) ) .OR. &
                  (  odd(ndivs) .AND.  odd(ieg-isg+1) ) .OR. &
                  (  odd(ndivs) .AND. even(ieg-isg+1) .AND. ndivs.LT.(ieg-isg+1)/2 )

             !mirror domains are stored in the list and retrieved if required.
             if( ndiv.EQ.0 )then
                !initialize max points and max domains
                imax = ieg
                ndmax = ndivs
             end if
             !do bottom half of decomposition, going over the midpoint for odd ndivs
             if( ndiv.LT.(ndivs-1)/2+1 )then
                !domain is sized by dividing remaining points by remaining domains
                ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
                ndmirror = (ndivs-1) - ndiv !mirror domain
                if( ndmirror.GT.ndiv .AND. symmetrize )then !only for domains over the midpoint
                   !mirror extents, the max(,) is to eliminate overlaps
                   ibegin(ndmirror) = max( isg+ieg-ie, ie+1 )
                   iend(ndmirror)   = max( isg+ieg-is, ie+1 )
                   imax = ibegin(ndmirror) - 1
                   ndmax = ndmax - 1
                end if
             else
                if( symmetrize )then
                   !do top half of decomposition by retrieving saved values
                   is = ibegin(ndiv)
                   ie = iend(ndiv)
                else
                   ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
                end if
             end if
          dim(ndiv) = ie-is+1
          is = ie + 1
       end do

  end subroutine DecomposeDim


