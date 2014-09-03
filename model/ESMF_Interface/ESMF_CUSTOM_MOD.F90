!***********************************************************
!* WARNING                                                 *
!*                                                         *
!* This file has been automatically generated via a rather *
!* perverse set of preprocessors.  Do _NOT_ directly edit  *
!* this file.  Instead, on should obtain the $?.nw (noweb) *
!* and make changes there.                                 *
!*                                                         *
!* For more information, contact:                          *
!*        Tom Clune <Thomas.L.Clune@nasa.gov>              *
!*                                                         *
!***********************************************************
#include "assert.h"
Module ESMF_CUSTOM_MOD
  Use ESMF_MOD, Only: tGrid => ESMF_Grid
  Use ESMF_MOD, Only: Field => ESMF_Field
  Use ESMF_MOD, Only: ESMF_Grid
  Use ESMF_MOD, Only: ESMF_Field
  Use ESMF_MOD, Only: KIND_R4 => ESMF_KIND_R4
  Use ESMF_MOD, Only: KIND_R8 => ESMF_KIND_R8
  Use ESMF_MOD, Only: ESMF_SUCCESS

!!$$ USE ESMF_MOD, Only: Field_HaloStore => ESMF_FieldHaloStore
!!$$ USE ESMF_MOD, Only: Field_Halo => ESMF_FieldHalo

!AOO  Use FieldComm, Only: Field_Halo
!AOO  Use FieldComm, Only: Field_Reduce
!!$$ Use FieldComm, Only: Field_Redist
!AOO  Use FieldComm, Only: Field_Gather
!!$$ Use FieldComm, Only: Field_Scatter

!AOO  USE FieldComm, Only: FIELD_SUM, FIELD_MAX, FIELD_MIN

  Use ESMF_MOD, Only: ESMF_CELL_NFACE
  Use ESMF_MOD, Only: ESMF_CELL_SFACE
  Use ESMF_MOD, Only: ESMF_CELL_CENTER
  Use ESMF_MOD, Only: ESMF_CELL_SWCORNER

  Use ESMF_MOD_private, Only: Field_GetDataPointer
  Use ESMF_MOD_private, Only: Field_SetDataPointer
!!$  Use ESMF_MOD_private, Only: NORTH, SOUTH


#ifdef USE_ESMF
  Use ESMF_MOD, Only: VM => ESMF_VM
  Use ESMF_MOD, Only: ESMF_DELayout
#else
#endif

  Implicit None
  Private

  Public :: tGrid
  Public :: Field

  Public :: Initialize_App
  Public :: Finalize_App

  Public :: Grid_GetExtents

  Public :: Field_Allocate
  Public :: Field_Deallocate
  Public :: Field_Clone
  Public :: Field_GetDataPointer
  Public :: Field_SetDataPointer

!AOO  Public :: Field_Halo
!AOO  Public :: Field_Reduce
  !!$$Public :: Field_Redist
!AOO  Public :: Field_Gather
  !!$$Public :: Field_Scatter

  Public :: ESMF_CELL_NFACE
  Public :: ESMF_CELL_SFACE
  Public :: ESMF_CELL_CENTER
  Public :: ESMF_CELL_SWCORNER


  Public :: KIND_R4
  Public :: KIND_R8
!AOO  Public :: FIELD_SUM, FIELD_MAX, FIELD_MIN
!!$  Public :: NORTH, SOUTH
  Public :: ESMF_SUCCESS

  ! To avoid the need to modify modelE interfaces for passing
  ! the ESMF grid throughout the source code, this one object is
  ! mode publicly available.
  Public :: modelE_grid
#ifdef USE_ESMF
  Public :: modelE_vm
  Public :: modelE_default_layout
#endif

  Logical :: initQ = .false.
  Type (ESMF_Grid), Target, Save :: modelE_grid
#ifdef USE_ESMF
  Type (VM), Target, Save :: modelE_vm
  Type (ESMF_DELayout), Target, Save :: modelE_default_layout
#endif


  Interface Field_Allocate
     Module Procedure Field_Allocate_Real4_3D
     Module Procedure Field_Allocate_Real8_1D
     Module Procedure Field_Allocate_Real8_2D
     Module Procedure Field_Allocate_Real8_3D
     Module Procedure Field_Allocate_Real8_4D
     Module Procedure Field_Allocate_Real8_5D
     Module Procedure Field_Allocate_Integer_1D
     Module Procedure Field_Allocate_Integer_2D
     Module Procedure Field_Allocate_Integer_3D
  End Interface


  Interface Field_Deallocate
     Module Procedure Field_Deallocate_Real4_3D
     Module Procedure Field_Deallocate_Real8_1D
     Module Procedure Field_Deallocate_Real8_2D
     Module Procedure Field_Deallocate_Real8_3D
     Module Procedure Field_Deallocate_Real8_4D
     Module Procedure Field_Deallocate_Real8_5D
     Module Procedure Field_Deallocate_Integer_1D
     Module Procedure Field_Deallocate_Integer_2D
     Module Procedure Field_Deallocate_Integer_3D
  End Interface

  Integer, Parameter :: N_DIMENSIONS = 3

Contains

  Subroutine Initialize_App(im, jm, lm, rc)
    Use ESMF_MOD, Only: ESMF_KIND_R8

#ifdef USE_ESMF
    Use ESMF_MOD, Only: ESMF_GRID_HORZ_STAGGER_A
    USE ESMF_MOD, Only: ESMF_Initialize
    USE ESMF_MOD, Only: ESMF_GridCreateHorzLatLonUni
    USE ESMF_MOD, Only: ESMF_SUCCESS
    USE ESMF_MOD, Only: ESMF_VMGet
    USE ESMF_MOD, Only: ESMF_DELayoutCreate
    USE ESMF_MOD, Only: ESMF_GridDistribute
    USE ESMF_MOD, Only: ESMF_GridGetDELocalAI
    USE ESMF_MOD, Only: ESMF_TRUE, ESMF_FALSE
    USE ESMF_MOD, Only: ESMF_GridAddVertHeight
#else
    USE ESMF_MOD_private, Only: Grid_Create
#endif
    USE ESMF_MOD, Only: ESMF_AxisIndex
    USE ESMF_MOD, Only: ESMF_CELL_CENTER

    Implicit None
    Integer, Intent(In) :: im ! num longitudes
    Integer, Intent(In) :: jm ! num latitudes (incl. poles)
    Integer, Intent(In) :: lm ! num levels
    Integer, Intent(Out), Optional :: rc ! return code
    Type (ESMF_AxisIndex) :: AI(N_DIMENSIONS)


    Real(Kind=ESMF_KIND_R8), Parameter :: lat_min = -90
    Real(Kind=ESMF_KIND_R8), Parameter :: lat_max = +90
    Real(Kind=ESMF_KIND_R8), Parameter :: lon_min = 0
    Real(Kind=ESMF_KIND_R8), Parameter :: lon_max = 360

#ifdef USE_ESMF
    Integer :: my_pet, n_pet, pet
    Integer :: L
    Type (ESMF_DELayout) :: modelE_DElayout
    REAL*8 :: deltaZ
#endif


#ifdef USE_ESMF
    Call ESMF_Initialize(vm=modelE_vm, rc=rc)
#endif

#ifdef USE_ESMF
    modelE_grid = ESMF_GridCreateHorzLatLonUni(counts = (/ IM, JM /), &
         & minGlobalCoordPerDim = (/ lon_min, lat_min /), &
         & maxGlobalCoordPerDim = (/ lon_max, lat_max /), &
         & horzStagger=ESMF_GRID_HORZ_STAGGER_A, &
         & name = "modelE default ESMF grid", &
         & periodic = (/ ESMF_TRUE, ESMF_FALSE /), rc=rc)

    deltaZ = 1.0d0
    call ESMF_GridAddVertHeight(modelE_grid, &
         delta=(/(deltaZ, L=1,LM) /),        &
         rc=rc)
!!$$    modelE_grid = ESMF_GridCreateHorzXYUni(counts = (/ IM, JM /), &
!!$$         & minGlobalCoordPerDim = (/ lon_min, lat_min /), &
!!$$         & maxGlobalCoordPerDim = (/ lon_max, lat_max /), &
!!$$         & horzStagger=ESMF_GRID_HORZ_STAGGER_B_SW, &
!!$$         & name = "modelE default ESMF grid", rc=rc)

    Call ESMF_VMGet(modelE_vm, localPET=my_pet, petCount=n_pet, rc=rc)
    modelE_DElayout = ESMF_DELayoutCreate(modelE_vm, deCountList = (/ 1, n_PET /))
    call ESMF_GridDistribute(grid = modelE_grid, delayout=modelE_DElayout, rc=rc)

    !!$$ Call ESMF_GridGetDELocalAI(modelE_grid, AIPerDim=AI, horzrelloc=ESMF_CELL_CENTER, rc=rc)

#else
    modelE_grid = Grid_Create(im,jm,lm)
#endif


    rc = ESMF_SUCCESS
  End Subroutine Initialize_App

  Subroutine Finalize_App()
#ifdef USE_ESMF
    USE ESMF_MOD, Only: ESMF_Finalize
    Implicit None
    Integer :: rc

    call ESMF_Finalize(rc=rc)
#endif

  End Subroutine Finalize_App


  Subroutine Grid_GetExtents(aGrid, default, stagger, skip_poles, halo, Have_South_pole, Have_North_Pole, J_equator)
    Use ESMF_MOD, Only: ESMF_AxisIndex
    Use ESMF_MOD, Only: ESMF_Grid
#ifdef USE_ESMF
    USE ESMF_MOD, Only: ESMF_GridGet
    USE ESMF_MOD, Only: ESMF_GridGetDELocalAI
    USE ESMF_MOD, Only: ESMF_RelLoc
    USE ESMF_MOD, Only: ESMF_CELL_CENTER
    USE ESMF_MOD, Only: ESMF_CELL_CELL
#else
    USE ESMF_MOD_private, Only : Grid_Get
#endif
    Implicit None
    ! Input
    Type (ESMF_Grid), Optional, Target, Intent(In) :: aGrid
    ! Output
    Type (ESMF_AxisIndex), Optional, Intent(Out) :: default(N_DIMENSIONS)
    Type (ESMF_AxisIndex), Optional, Intent(Out) :: stagger(N_DIMENSIONS)
    Type (ESMF_AxisIndex), Optional, Intent(Out) :: skip_poles(N_DIMENSIONS)
    Type (ESMF_AxisIndex), Optional, Intent(Out) :: halo(N_DIMENSIONS)
    Logical, Optional, Intent(Out) :: Have_North_Pole
    Logical, Optional, Intent(Out) :: Have_South_Pole
    Integer, Optional, Intent(Out) :: J_equator

    Type (ESMF_Grid), Pointer :: pGrid
    Integer :: grid_size(3) ! global
    Integer :: rc
#ifdef USE_ESMF
    Type (ESMF_RelLoc) :: horzRelLoc = ESMF_CELL_CENTER
#else
    Integer :: horzRelLoc = ESMF_CELL_CENTER
#endif
    Type (ESMF_AxisIndex) :: AI_local(3)

    pGrid => modelE_grid ! default
    If (Present(aGrid)) pGrid => aGrid

#ifdef USE_ESMF
    Call ESMF_GridGet(pGrid, horzRelLoc, globalCellCountPerDim = grid_size)
    Call ESMF_GridGetDELocalAI(pGrid, AIPerDim=AI_local, horzrelloc=horzRelLoc, &
    & vertrelloc=ESMF_CELL_CELL,rc=rc)
#else
    Call Grid_Get(pGrid, AI_Local)
    grid_size = AI_Local%max
#endif

    If (Present(Have_South_Pole)) Have_South_Pole = (AI_local(2)%min == 1)
    If (Present(Have_North_Pole)) Have_North_Pole = (AI_local(2)%max == grid_size(2))

    If (Present(default)) default = AI_local
    If (Present(stagger)) Then
       stagger = AI_Local
       stagger(2)%min = MAX(2, stagger(2)%min)
    End If
    If (Present(skip_poles)) Then
       skip_poles = AI_local
       skip_poles(2)%min = MAX(AI_local(2)%min, 2)
       skip_poles(2)%max = MIN(AI_local(2)%max, grid_size(2)-1)
    End If

    If (Present(halo)) Then
#ifdef USE_ESMF
       halo(1) = AI_local(1)
  !!$$ halo(3) = AI_local(3)

       Call ESMF_GridGetDELocalAI(pGrid, AIPerDim=AI_local,  horzrelloc=horzRelLoc, vertrelloc=ESMF_CELL_CELL, total=.true., rc=rc)
       halo(2)%min = AI_local(2)%min - 1
       halo(2)%max = AI_local(2)%max - 1
#else
       halo = AI_local ! halo = usual domain
#endif
    End If

    If (Present(J_equator)) J_equator = grid_size(2)/2

  End Subroutine Grid_GetExtents
  
  Function Field_Allocate_Real4_3D ( & 
    &     aPtr, aGrid, shape, halo_width, grid_rank_index, horzRelLoc, name, rc) Result(aField) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_MAXSTR 
#ifdef USE_ESMF 
    Use ESMF_MOD, Only: ESMF_FieldDataMap 
    Use ESMF_MOD, Only: ESMF_RelLoc 
    Use ESMF_MOD, Only: ESMF_FieldCreate 
    Use ESMF_MOD, Only: ESMF_FieldDataMapSet 
    Use ESMF_MOD, Only: ESMF_Array 
    Use ESMF_MOD, Only: ESMF_ArrayCreate 
    Use ESMF_MOD, Only: ESMF_CELL_CENTER 
    Use ESMF_MOD, Only: ESMF_CELL_CELL
    Use ESMF_MOD, Only: ESMF_DATA_REF 
    Use ESMF_MOD, Only: ESMF_GridGetDELocalAI 
    Use ESMF_MOD, Only: ESMF_ArrayGet 
#else 
    USE ESMF_MOD, Only: ESMF_FieldSetDataPointer 
    Use ESMF_MOD_private, Only: Field_SetDistRank 
#endif 
    Implicit None 
    Type (ESMF_Field) :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer           :: aPtr (:,:,:) 
    Type (ESMF_Grid), Intent(In)          :: aGrid 
    Integer, Intent(In), Optional    :: shape(:) 
    Integer, Intent(In), Optional    :: halo_width 
    Integer, Intent(In), Optional    :: grid_rank_index(:) 
#ifdef USE_ESMF 
    Type (ESMF_RelLoc), Intent(in), Optional    :: horzRelLoc 
#else 
    Integer, Intent(In), Optional :: horzRelLoc 
#endif 
    Character(Len=*), Intent(In), Optional :: name 
    Integer, Intent(Out), Optional   :: rc 
 
    Type (ESMF_AxisIndex) :: AI(N_DIMENSIONS) 
#ifdef USE_ESMF 
    Type (ESMF_FieldDataMap) :: map 
    Type (ESMF_Array) :: tmpArray 
    Type (ESMF_RelLoc) :: hzRelLoc 
#endif 
    Integer            :: hwidth 
    Character(Len=ESMF_MAXSTR) :: name_ 
 
    Integer, Dimension(3) :: l, u ! lower and upper bounds 
    Integer, Dimension(3) :: dataIndices 
    Integer :: j_min, j_max 
    Integer :: stat 
    Integer :: idx, gidx 
    Integer :: j_index, i_index 
 
    !------------------------------------------- 
 
#ifdef USE_ESMF 
      hzRelLoc = ESMF_CELL_CENTER 
      If (Present(horzRelLoc)) hzRelLoc = horzRelLoc 
#endif 
      hwidth   = 1 
      name_    = 'unnamed' 
 
      If (Present(halo_width)) Then 
         ASSERT(halo_width >= 0, 'illegal value for hwidth') 
         hwidth   = halo_width 
      End If 
 
      If (Present(name)) name_ = Trim(name) 
 
    l = 1      ! default lower bound 
    u = shape  ! default upper bound 
 
#ifdef USE_ESMF 
    Call ESMF_GridGetDELocalAI(aGrid, AIPerDim=AI, horzrelloc=hzRelLoc, rc=rc) 
 
    ! Extract bounds of computational domain 
    ! Adjust lower/upper bounds with halo size 
    Do idx = 1, 3 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = AI(gidx)%min - hwidth ! Would like hwidth to be independent of axis, but ESMF.1.0.7 
          u(idx) = AI(gidx)%max + hwidth ! does not support this. 
       End If 
    End Do 
#else 
    Do idx = 1, 3 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = l(idx)-hwidth 
          u(idx) = u(idx)+hwidth 
       End If 
    End Do 
#endif 
 
 
    Allocate(aPtr(l(1):u(1),l(2):u(2),l(3):u(3)), STAT=stat) 
    ASSERT_ALWAYS(stat == 0, 'Allocate Failure') 
 
#ifdef USE_ESMF 
    Call ESMF_FieldDataMapSet(map, 3, grid_rank_index, shape, horzRelloc=hzRelLoc, rc=rc) 
    tmpArray = ESMF_ArrayCreate(aPtr, ESMF_DATA_REF, haloWidth=hwidth, rc=rc) 
    Call ESMF_ArrayGet(tmparray, haloWidth=hwidth) 
    aField   = ESMF_FieldCreate(aGrid,  tmpArray, halowidth=hwidth, horzRelloc=hzRelLoc, name=name_, datamap=map, rc=rc) 
#else 
 
    Call ESMF_FieldSetDataPointer(aField,aPtr) 
    Do idx = 1, 3 
       gidx =grid_rank_index(idx) 
       If (gidx > 0) Exit 
    End Do 
 
    Call Field_SetDistRank(aField,idx) 
#endif 
 
 
    rc = ESMF_SUCCESS 
 
  End Function Field_Allocate_Real4_3D 

  
  Function Field_Allocate_Real8_1D ( & 
    &     aPtr, aGrid, shape, halo_width, grid_rank_index, horzRelLoc, name, rc) Result(aField) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_MAXSTR 
#ifdef USE_ESMF 
    Use ESMF_MOD, Only: ESMF_FieldDataMap 
    Use ESMF_MOD, Only: ESMF_RelLoc 
    Use ESMF_MOD, Only: ESMF_FieldCreate 
    Use ESMF_MOD, Only: ESMF_FieldDataMapSet 
    Use ESMF_MOD, Only: ESMF_Array 
    Use ESMF_MOD, Only: ESMF_ArrayCreate 
    Use ESMF_MOD, Only: ESMF_CELL_CENTER 
    Use ESMF_MOD, Only: ESMF_DATA_REF 
    Use ESMF_MOD, Only: ESMF_GridGetDELocalAI 
    Use ESMF_MOD, Only: ESMF_ArrayGet 
#else 
    USE ESMF_MOD, Only: ESMF_FieldSetDataPointer 
    Use ESMF_MOD_private, Only: Field_SetDistRank 
#endif 
    Implicit None 
    Type (ESMF_Field) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer           :: aPtr (:) 
    Type (ESMF_Grid), Intent(In)          :: aGrid 
    Integer, Intent(In), Optional    :: shape(:) 
    Integer, Intent(In), Optional    :: halo_width 
    Integer, Intent(In), Optional    :: grid_rank_index(:) 
#ifdef USE_ESMF 
    Type (ESMF_RelLoc), Intent(in), Optional    :: horzRelLoc 
#else 
    Integer, Intent(In), Optional :: horzRelLoc 
#endif 
    Character(Len=*), Intent(In), Optional :: name 
    Integer, Intent(Out), Optional   :: rc 
 
    Type (ESMF_AxisIndex) :: AI(N_DIMENSIONS) 
#ifdef USE_ESMF 
    Type (ESMF_FieldDataMap) :: map 
    Type (ESMF_Array) :: tmpArray 
    Type (ESMF_RelLoc) :: hzRelLoc 
#endif 
    Integer            :: hwidth 
    Character(Len=ESMF_MAXSTR) :: name_ 
 
    Integer, Dimension(1) :: l, u ! lower and upper bounds 
    Integer, Dimension(1) :: dataIndices 
    Integer :: j_min, j_max 
    Integer :: stat 
    Integer :: idx, gidx 
    Integer :: j_index, i_index 
 
    !------------------------------------------- 
 
#ifdef USE_ESMF 
      hzRelLoc = ESMF_CELL_CENTER 
      If (Present(horzRelLoc)) hzRelLoc = horzRelLoc 
#endif 
      hwidth   = 1 
      name_    = 'unnamed' 
 
      If (Present(halo_width)) Then 
         ASSERT(halo_width >= 0, 'illegal value for hwidth') 
         hwidth   = halo_width 
      End If 
 
      If (Present(name)) name_ = Trim(name) 
 
    l = 1      ! default lower bound 
    u = shape  ! default upper bound 
 
#ifdef USE_ESMF 
    Call ESMF_GridGetDELocalAI(aGrid, AIPerDim=AI, horzrelloc=hzRelLoc, rc=rc) 
 
    ! Extract bounds of computational domain 
    ! Adjust lower/upper bounds with halo size 
    Do idx = 1, 1 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = AI(gidx)%min - hwidth ! Would like hwidth to be independent of axis, but ESMF.1.0.7 
          u(idx) = AI(gidx)%max + hwidth ! does not support this. 
       End If 
    End Do 
#else 
    Do idx = 1, 1 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = l(idx)-hwidth 
          u(idx) = u(idx)+hwidth 
       End If 
    End Do 
#endif 
 
 
    Allocate(aPtr(l(1):u(1)), STAT=stat) 
    ASSERT_ALWAYS(stat == 0, 'Allocate Failure') 
 
#ifdef USE_ESMF 
    Call ESMF_FieldDataMapSet(map, 1, grid_rank_index, shape, horzRelloc=hzRelLoc, rc=rc) 
    tmpArray = ESMF_ArrayCreate(aPtr, ESMF_DATA_REF, haloWidth=hwidth, rc=rc) 
    Call ESMF_ArrayGet(tmparray, haloWidth=hwidth) 
    aField   = ESMF_FieldCreate(aGrid,  tmpArray, halowidth=hwidth, horzRelloc=hzRelLoc, name=name_, datamap=map, rc=rc) 
#else 
 
    Call ESMF_FieldSetDataPointer(aField,aPtr) 
    Do idx = 1, 1 
       gidx =grid_rank_index(idx) 
       If (gidx > 0) Exit 
    End Do 
 
    Call Field_SetDistRank(aField,idx) 
#endif 
 
 
    rc = ESMF_SUCCESS 
 
  End Function Field_Allocate_Real8_1D 

  
  Function Field_Allocate_Real8_2D ( & 
    &     aPtr, aGrid, shape, halo_width, grid_rank_index, horzRelLoc, name, rc) Result(aField) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_MAXSTR 
#ifdef USE_ESMF 
    Use ESMF_MOD, Only: ESMF_FieldDataMap 
    Use ESMF_MOD, Only: ESMF_RelLoc 
    Use ESMF_MOD, Only: ESMF_FieldCreate 
    Use ESMF_MOD, Only: ESMF_FieldDataMapSet 
    Use ESMF_MOD, Only: ESMF_Array 
    Use ESMF_MOD, Only: ESMF_ArrayCreate 
    Use ESMF_MOD, Only: ESMF_CELL_CENTER 
    Use ESMF_MOD, Only: ESMF_DATA_REF 
    Use ESMF_MOD, Only: ESMF_GridGetDELocalAI 
    Use ESMF_MOD, Only: ESMF_ArrayGet 
#else 
    USE ESMF_MOD, Only: ESMF_FieldSetDataPointer 
    Use ESMF_MOD_private, Only: Field_SetDistRank 
#endif 
    Implicit None 
    Type (ESMF_Field) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer           :: aPtr (:,:) 
    Type (ESMF_Grid), Intent(In)          :: aGrid 
    Integer, Intent(In), Optional    :: shape(:) 
    Integer, Intent(In), Optional    :: halo_width 
    Integer, Intent(In), Optional    :: grid_rank_index(:) 
#ifdef USE_ESMF 
    Type (ESMF_RelLoc), Intent(in), Optional    :: horzRelLoc 
#else 
    Integer, Intent(In), Optional :: horzRelLoc 
#endif 
    Character(Len=*), Intent(In), Optional :: name 
    Integer, Intent(Out), Optional   :: rc 
 
    Type (ESMF_AxisIndex) :: AI(N_DIMENSIONS) 
#ifdef USE_ESMF 
    Type (ESMF_FieldDataMap) :: map 
    Type (ESMF_Array) :: tmpArray 
    Type (ESMF_RelLoc) :: hzRelLoc 
#endif 
    Integer            :: hwidth 
    Character(Len=ESMF_MAXSTR) :: name_ 
 
    Integer, Dimension(2) :: l, u ! lower and upper bounds 
    Integer, Dimension(2) :: dataIndices 
    Integer :: j_min, j_max 
    Integer :: stat 
    Integer :: idx, gidx 
    Integer :: j_index, i_index 
 
    !------------------------------------------- 
 
#ifdef USE_ESMF 
      hzRelLoc = ESMF_CELL_CENTER 
      If (Present(horzRelLoc)) hzRelLoc = horzRelLoc 
#endif 
      hwidth   = 1 
      name_    = 'unnamed' 
 
      If (Present(halo_width)) Then 
         ASSERT(halo_width >= 0, 'illegal value for hwidth') 
         hwidth   = halo_width 
      End If 
 
      If (Present(name)) name_ = Trim(name) 
 
    l = 1      ! default lower bound 
    u = shape  ! default upper bound 
 
#ifdef USE_ESMF 
    Call ESMF_GridGetDELocalAI(aGrid, AIPerDim=AI, horzrelloc=hzRelLoc, rc=rc) 
 
    ! Extract bounds of computational domain 
    ! Adjust lower/upper bounds with halo size 
    Do idx = 1, 2 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = AI(gidx)%min - hwidth ! Would like hwidth to be independent of axis, but ESMF.1.0.7 
          u(idx) = AI(gidx)%max + hwidth ! does not support this. 
       End If 
    End Do 
#else 
    Do idx = 1, 2 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = l(idx)-hwidth 
          u(idx) = u(idx)+hwidth 
       End If 
    End Do 
#endif 
 
 
    Allocate(aPtr(l(1):u(1),l(2):u(2)), STAT=stat) 
    ASSERT_ALWAYS(stat == 0, 'Allocate Failure') 
 
#ifdef USE_ESMF 
    Call ESMF_FieldDataMapSet(map, 2, grid_rank_index, shape, horzRelloc=hzRelLoc, rc=rc) 
    tmpArray = ESMF_ArrayCreate(aPtr, ESMF_DATA_REF, haloWidth=hwidth, rc=rc) 
    Call ESMF_ArrayGet(tmparray, haloWidth=hwidth) 
    aField   = ESMF_FieldCreate(aGrid,  tmpArray, halowidth=hwidth, horzRelloc=hzRelLoc, name=name_, datamap=map, rc=rc) 
#else 
 
    Call ESMF_FieldSetDataPointer(aField,aPtr) 
    Do idx = 1, 2 
       gidx =grid_rank_index(idx) 
       If (gidx > 0) Exit 
    End Do 
 
    Call Field_SetDistRank(aField,idx) 
#endif 
 
 
    rc = ESMF_SUCCESS 
 
  End Function Field_Allocate_Real8_2D 

  
  Function Field_Allocate_Real8_3D ( & 
    &     aPtr, aGrid, shape, halo_width, grid_rank_index, horzRelLoc, name, rc) Result(aField) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_MAXSTR 
#ifdef USE_ESMF 
    Use ESMF_MOD, Only: ESMF_FieldDataMap 
    Use ESMF_MOD, Only: ESMF_RelLoc 
    Use ESMF_MOD, Only: ESMF_FieldCreate 
    Use ESMF_MOD, Only: ESMF_FieldDataMapSet 
    Use ESMF_MOD, Only: ESMF_Array 
    Use ESMF_MOD, Only: ESMF_ArrayCreate 
    Use ESMF_MOD, Only: ESMF_CELL_CENTER 
    Use ESMF_MOD, Only: ESMF_DATA_REF 
    Use ESMF_MOD, Only: ESMF_GridGetDELocalAI 
    Use ESMF_MOD, Only: ESMF_ArrayGet 
#else 
    USE ESMF_MOD, Only: ESMF_FieldSetDataPointer 
    Use ESMF_MOD_private, Only: Field_SetDistRank 
#endif 
    Implicit None 
    Type (ESMF_Field) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer           :: aPtr (:,:,:) 
    Type (ESMF_Grid), Intent(In)          :: aGrid 
    Integer, Intent(In), Optional    :: shape(:) 
    Integer, Intent(In), Optional    :: halo_width 
    Integer, Intent(In), Optional    :: grid_rank_index(:) 
#ifdef USE_ESMF 
    Type (ESMF_RelLoc), Intent(in), Optional    :: horzRelLoc 
#else 
    Integer, Intent(In), Optional :: horzRelLoc 
#endif 
    Character(Len=*), Intent(In), Optional :: name 
    Integer, Intent(Out), Optional   :: rc 
 
    Type (ESMF_AxisIndex) :: AI(N_DIMENSIONS) 
#ifdef USE_ESMF 
    Type (ESMF_FieldDataMap) :: map 
    Type (ESMF_Array) :: tmpArray 
    Type (ESMF_RelLoc) :: hzRelLoc 
#endif 
    Integer            :: hwidth 
    Character(Len=ESMF_MAXSTR) :: name_ 
 
    Integer, Dimension(3) :: l, u ! lower and upper bounds 
    Integer, Dimension(3) :: dataIndices 
    Integer :: j_min, j_max 
    Integer :: stat 
    Integer :: idx, gidx 
    Integer :: j_index, i_index 
 
    !------------------------------------------- 
 
#ifdef USE_ESMF 
      hzRelLoc = ESMF_CELL_CENTER 
      If (Present(horzRelLoc)) hzRelLoc = horzRelLoc 
#endif 
      hwidth   = 1 
      name_    = 'unnamed' 
 
      If (Present(halo_width)) Then 
         ASSERT(halo_width >= 0, 'illegal value for hwidth') 
         hwidth   = halo_width 
      End If 
 
      If (Present(name)) name_ = Trim(name) 
 
    l = 1      ! default lower bound 
    u = shape  ! default upper bound 
 
#ifdef USE_ESMF 
    Call ESMF_GridGetDELocalAI(aGrid, AIPerDim=AI, horzrelloc=hzRelLoc, rc=rc) 
 
    ! Extract bounds of computational domain 
    ! Adjust lower/upper bounds with halo size 
    Do idx = 1, 3 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = AI(gidx)%min - hwidth ! Would like hwidth to be independent of axis, but ESMF.1.0.7 
          u(idx) = AI(gidx)%max + hwidth ! does not support this. 
       End If 
    End Do 
#else 
    Do idx = 1, 3 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = l(idx)-hwidth 
          u(idx) = u(idx)+hwidth 
       End If 
    End Do 
#endif 
 
 
    Allocate(aPtr(l(1):u(1),l(2):u(2),l(3):u(3)), STAT=stat) 
    ASSERT_ALWAYS(stat == 0, 'Allocate Failure') 
 
#ifdef USE_ESMF 
    Call ESMF_FieldDataMapSet(map, 3, grid_rank_index, shape, horzRelloc=hzRelLoc, rc=rc) 
    tmpArray = ESMF_ArrayCreate(aPtr, ESMF_DATA_REF, haloWidth=hwidth, rc=rc) 
    Call ESMF_ArrayGet(tmparray, haloWidth=hwidth) 
    aField   = ESMF_FieldCreate(aGrid,  tmpArray, halowidth=hwidth, horzRelloc=hzRelLoc, name=name_, datamap=map, rc=rc) 
#else 
 
    Call ESMF_FieldSetDataPointer(aField,aPtr) 
    Do idx = 1, 3 
       gidx =grid_rank_index(idx) 
       If (gidx > 0) Exit 
    End Do 
 
    Call Field_SetDistRank(aField,idx) 
#endif 
 
 
    rc = ESMF_SUCCESS 
 
  End Function Field_Allocate_Real8_3D 

  
  Function Field_Allocate_Real8_4D ( & 
    &     aPtr, aGrid, shape, halo_width, grid_rank_index, horzRelLoc, name, rc) Result(aField) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_MAXSTR 
#ifdef USE_ESMF 
    Use ESMF_MOD, Only: ESMF_FieldDataMap 
    Use ESMF_MOD, Only: ESMF_RelLoc 
    Use ESMF_MOD, Only: ESMF_FieldCreate 
    Use ESMF_MOD, Only: ESMF_FieldDataMapSet 
    Use ESMF_MOD, Only: ESMF_Array 
    Use ESMF_MOD, Only: ESMF_ArrayCreate 
    Use ESMF_MOD, Only: ESMF_CELL_CENTER 
    Use ESMF_MOD, Only: ESMF_DATA_REF 
    Use ESMF_MOD, Only: ESMF_GridGetDELocalAI 
    Use ESMF_MOD, Only: ESMF_ArrayGet 
#else 
    USE ESMF_MOD, Only: ESMF_FieldSetDataPointer 
    Use ESMF_MOD_private, Only: Field_SetDistRank 
#endif 
    Implicit None 
    Type (ESMF_Field) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer           :: aPtr (:,:,:,:) 
    Type (ESMF_Grid), Intent(In)          :: aGrid 
    Integer, Intent(In), Optional    :: shape(:) 
    Integer, Intent(In), Optional    :: halo_width 
    Integer, Intent(In), Optional    :: grid_rank_index(:) 
#ifdef USE_ESMF 
    Type (ESMF_RelLoc), Intent(in), Optional    :: horzRelLoc 
#else 
    Integer, Intent(In), Optional :: horzRelLoc 
#endif 
    Character(Len=*), Intent(In), Optional :: name 
    Integer, Intent(Out), Optional   :: rc 
 
    Type (ESMF_AxisIndex) :: AI(N_DIMENSIONS) 
#ifdef USE_ESMF 
    Type (ESMF_FieldDataMap) :: map 
    Type (ESMF_Array) :: tmpArray 
    Type (ESMF_RelLoc) :: hzRelLoc 
#endif 
    Integer            :: hwidth 
    Character(Len=ESMF_MAXSTR) :: name_ 
 
    Integer, Dimension(4) :: l, u ! lower and upper bounds 
    Integer, Dimension(4) :: dataIndices 
    Integer :: j_min, j_max 
    Integer :: stat 
    Integer :: idx, gidx 
    Integer :: j_index, i_index 
 
    !------------------------------------------- 
 
#ifdef USE_ESMF 
      hzRelLoc = ESMF_CELL_CENTER 
      If (Present(horzRelLoc)) hzRelLoc = horzRelLoc 
#endif 
      hwidth   = 1 
      name_    = 'unnamed' 
 
      If (Present(halo_width)) Then 
         ASSERT(halo_width >= 0, 'illegal value for hwidth') 
         hwidth   = halo_width 
      End If 
 
      If (Present(name)) name_ = Trim(name) 
 
    l = 1      ! default lower bound 
    u = shape  ! default upper bound 
 
#ifdef USE_ESMF 
    Call ESMF_GridGetDELocalAI(aGrid, AIPerDim=AI, horzrelloc=hzRelLoc, rc=rc) 
 
    ! Extract bounds of computational domain 
    ! Adjust lower/upper bounds with halo size 
    Do idx = 1, 4 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = AI(gidx)%min - hwidth ! Would like hwidth to be independent of axis, but ESMF.1.0.7 
          u(idx) = AI(gidx)%max + hwidth ! does not support this. 
       End If 
    End Do 
#else 
    Do idx = 1, 4 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = l(idx)-hwidth 
          u(idx) = u(idx)+hwidth 
       End If 
    End Do 
#endif 
 
 
    Allocate(aPtr(l(1):u(1),l(2):u(2),l(3):u(3),l(4):u(4)), STAT=stat) 
    ASSERT_ALWAYS(stat == 0, 'Allocate Failure') 
 
#ifdef USE_ESMF 
    Call ESMF_FieldDataMapSet(map, 4, grid_rank_index, shape, horzRelloc=hzRelLoc, rc=rc) 
    tmpArray = ESMF_ArrayCreate(aPtr, ESMF_DATA_REF, haloWidth=hwidth, rc=rc) 
    Call ESMF_ArrayGet(tmparray, haloWidth=hwidth) 
    aField   = ESMF_FieldCreate(aGrid,  tmpArray, halowidth=hwidth, horzRelloc=hzRelLoc, name=name_, datamap=map, rc=rc) 
#else 
 
    Call ESMF_FieldSetDataPointer(aField,aPtr) 
    Do idx = 1, 4 
       gidx =grid_rank_index(idx) 
       If (gidx > 0) Exit 
    End Do 
 
    Call Field_SetDistRank(aField,idx) 
#endif 
 
 
    rc = ESMF_SUCCESS 
 
  End Function Field_Allocate_Real8_4D 

  
  Function Field_Allocate_Real8_5D ( & 
    &     aPtr, aGrid, shape, halo_width, grid_rank_index, horzRelLoc, name, rc) Result(aField) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_MAXSTR 
#ifdef USE_ESMF 
    Use ESMF_MOD, Only: ESMF_FieldDataMap 
    Use ESMF_MOD, Only: ESMF_RelLoc 
    Use ESMF_MOD, Only: ESMF_FieldCreate 
    Use ESMF_MOD, Only: ESMF_FieldDataMapSet 
    Use ESMF_MOD, Only: ESMF_Array 
    Use ESMF_MOD, Only: ESMF_ArrayCreate 
    Use ESMF_MOD, Only: ESMF_CELL_CENTER 
    Use ESMF_MOD, Only: ESMF_DATA_REF 
    Use ESMF_MOD, Only: ESMF_GridGetDELocalAI 
    Use ESMF_MOD, Only: ESMF_ArrayGet 
#else 
    USE ESMF_MOD, Only: ESMF_FieldSetDataPointer 
    Use ESMF_MOD_private, Only: Field_SetDistRank 
#endif 
    Implicit None 
    Type (ESMF_Field) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer           :: aPtr (:,:,:,:,:) 
    Type (ESMF_Grid), Intent(In)          :: aGrid 
    Integer, Intent(In), Optional    :: shape(:) 
    Integer, Intent(In), Optional    :: halo_width 
    Integer, Intent(In), Optional    :: grid_rank_index(:) 
#ifdef USE_ESMF 
    Type (ESMF_RelLoc), Intent(in), Optional    :: horzRelLoc 
#else 
    Integer, Intent(In), Optional :: horzRelLoc 
#endif 
    Character(Len=*), Intent(In), Optional :: name 
    Integer, Intent(Out), Optional   :: rc 
 
    Type (ESMF_AxisIndex) :: AI(N_DIMENSIONS) 
#ifdef USE_ESMF 
    Type (ESMF_FieldDataMap) :: map 
    Type (ESMF_Array) :: tmpArray 
    Type (ESMF_RelLoc) :: hzRelLoc 
#endif 
    Integer            :: hwidth 
    Character(Len=ESMF_MAXSTR) :: name_ 
 
    Integer, Dimension(5) :: l, u ! lower and upper bounds 
    Integer, Dimension(5) :: dataIndices 
    Integer :: j_min, j_max 
    Integer :: stat 
    Integer :: idx, gidx 
    Integer :: j_index, i_index 
 
    !------------------------------------------- 
 
#ifdef USE_ESMF 
      hzRelLoc = ESMF_CELL_CENTER 
      If (Present(horzRelLoc)) hzRelLoc = horzRelLoc 
#endif 
      hwidth   = 1 
      name_    = 'unnamed' 
 
      If (Present(halo_width)) Then 
         ASSERT(halo_width >= 0, 'illegal value for hwidth') 
         hwidth   = halo_width 
      End If 
 
      If (Present(name)) name_ = Trim(name) 
 
    l = 1      ! default lower bound 
    u = shape  ! default upper bound 
 
#ifdef USE_ESMF 
    Call ESMF_GridGetDELocalAI(aGrid, AIPerDim=AI, horzrelloc=hzRelLoc, rc=rc) 
 
    ! Extract bounds of computational domain 
    ! Adjust lower/upper bounds with halo size 
    Do idx = 1, 5 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = AI(gidx)%min - hwidth ! Would like hwidth to be independent of axis, but ESMF.1.0.7 
          u(idx) = AI(gidx)%max + hwidth ! does not support this. 
       End If 
    End Do 
#else 
    Do idx = 1, 5 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = l(idx)-hwidth 
          u(idx) = u(idx)+hwidth 
       End If 
    End Do 
#endif 
 
 
    Allocate(aPtr(l(1):u(1),l(2):u(2),l(3):u(3),l(4):u(4),l(5):u(5)), STAT=stat) 
    ASSERT_ALWAYS(stat == 0, 'Allocate Failure') 
 
#ifdef USE_ESMF 
    Call ESMF_FieldDataMapSet(map, 5, grid_rank_index, shape, horzRelloc=hzRelLoc, rc=rc) 
    tmpArray = ESMF_ArrayCreate(aPtr, ESMF_DATA_REF, haloWidth=hwidth, rc=rc) 
    Call ESMF_ArrayGet(tmparray, haloWidth=hwidth) 
    aField   = ESMF_FieldCreate(aGrid,  tmpArray, halowidth=hwidth, horzRelloc=hzRelLoc, name=name_, datamap=map, rc=rc) 
#else 
 
    Call ESMF_FieldSetDataPointer(aField,aPtr) 
    Do idx = 1, 5 
       gidx =grid_rank_index(idx) 
       If (gidx > 0) Exit 
    End Do 
 
    Call Field_SetDistRank(aField,idx) 
#endif 
 
 
    rc = ESMF_SUCCESS 
 
  End Function Field_Allocate_Real8_5D 

  
  Function Field_Allocate_Integer_1D ( & 
    &     aPtr, aGrid, shape, halo_width, grid_rank_index, horzRelLoc, name, rc) Result(aField) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_MAXSTR 
#ifdef USE_ESMF 
    Use ESMF_MOD, Only: ESMF_FieldDataMap 
    Use ESMF_MOD, Only: ESMF_RelLoc 
    Use ESMF_MOD, Only: ESMF_FieldCreate 
    Use ESMF_MOD, Only: ESMF_FieldDataMapSet 
    Use ESMF_MOD, Only: ESMF_Array 
    Use ESMF_MOD, Only: ESMF_ArrayCreate 
    Use ESMF_MOD, Only: ESMF_CELL_CENTER 
    Use ESMF_MOD, Only: ESMF_DATA_REF 
    Use ESMF_MOD, Only: ESMF_GridGetDELocalAI 
    Use ESMF_MOD, Only: ESMF_ArrayGet 
#else 
    USE ESMF_MOD, Only: ESMF_FieldSetDataPointer 
    Use ESMF_MOD_private, Only: Field_SetDistRank 
#endif 
    Implicit None 
    Type (ESMF_Field) :: aField 
    Integer , Pointer           :: aPtr (:) 
    Type (ESMF_Grid), Intent(In)          :: aGrid 
    Integer, Intent(In), Optional    :: shape(:) 
    Integer, Intent(In), Optional    :: halo_width 
    Integer, Intent(In), Optional    :: grid_rank_index(:) 
#ifdef USE_ESMF 
    Type (ESMF_RelLoc), Intent(in), Optional    :: horzRelLoc 
#else 
    Integer, Intent(In), Optional :: horzRelLoc 
#endif 
    Character(Len=*), Intent(In), Optional :: name 
    Integer, Intent(Out), Optional   :: rc 
 
    Type (ESMF_AxisIndex) :: AI(N_DIMENSIONS) 
#ifdef USE_ESMF 
    Type (ESMF_FieldDataMap) :: map 
    Type (ESMF_Array) :: tmpArray 
    Type (ESMF_RelLoc) :: hzRelLoc 
#endif 
    Integer            :: hwidth 
    Character(Len=ESMF_MAXSTR) :: name_ 
 
    Integer, Dimension(1) :: l, u ! lower and upper bounds 
    Integer, Dimension(1) :: dataIndices 
    Integer :: j_min, j_max 
    Integer :: stat 
    Integer :: idx, gidx 
    Integer :: j_index, i_index 
 
    !------------------------------------------- 
 
#ifdef USE_ESMF 
      hzRelLoc = ESMF_CELL_CENTER 
      If (Present(horzRelLoc)) hzRelLoc = horzRelLoc 
#endif 
      hwidth   = 1 
      name_    = 'unnamed' 
 
      If (Present(halo_width)) Then 
         ASSERT(halo_width >= 0, 'illegal value for hwidth') 
         hwidth   = halo_width 
      End If 
 
      If (Present(name)) name_ = Trim(name) 
 
    l = 1      ! default lower bound 
    u = shape  ! default upper bound 
 
#ifdef USE_ESMF 
    Call ESMF_GridGetDELocalAI(aGrid, AIPerDim=AI, horzrelloc=hzRelLoc, rc=rc) 
 
    ! Extract bounds of computational domain 
    ! Adjust lower/upper bounds with halo size 
    Do idx = 1, 1 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = AI(gidx)%min - hwidth ! Would like hwidth to be independent of axis, but ESMF.1.0.7 
          u(idx) = AI(gidx)%max + hwidth ! does not support this. 
       End If 
    End Do 
#else 
    Do idx = 1, 1 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = l(idx)-hwidth 
          u(idx) = u(idx)+hwidth 
       End If 
    End Do 
#endif 
 
 
    Allocate(aPtr(l(1):u(1)), STAT=stat) 
    ASSERT_ALWAYS(stat == 0, 'Allocate Failure') 
 
#ifdef USE_ESMF 
    Call ESMF_FieldDataMapSet(map, 1, grid_rank_index, shape, horzRelloc=hzRelLoc, rc=rc) 
    tmpArray = ESMF_ArrayCreate(aPtr, ESMF_DATA_REF, haloWidth=hwidth, rc=rc) 
    Call ESMF_ArrayGet(tmparray, haloWidth=hwidth) 
    aField   = ESMF_FieldCreate(aGrid,  tmpArray, halowidth=hwidth, horzRelloc=hzRelLoc, name=name_, datamap=map, rc=rc) 
#else 
 
    Call ESMF_FieldSetDataPointer(aField,aPtr) 
    Do idx = 1, 1 
       gidx =grid_rank_index(idx) 
       If (gidx > 0) Exit 
    End Do 
 
    Call Field_SetDistRank(aField,idx) 
#endif 
 
 
    rc = ESMF_SUCCESS 
 
  End Function Field_Allocate_Integer_1D 

  
  Function Field_Allocate_Integer_2D ( & 
    &     aPtr, aGrid, shape, halo_width, grid_rank_index, horzRelLoc, name, rc) Result(aField) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_MAXSTR 
#ifdef USE_ESMF 
    Use ESMF_MOD, Only: ESMF_FieldDataMap 
    Use ESMF_MOD, Only: ESMF_RelLoc 
    Use ESMF_MOD, Only: ESMF_FieldCreate 
    Use ESMF_MOD, Only: ESMF_FieldDataMapSet 
    Use ESMF_MOD, Only: ESMF_Array 
    Use ESMF_MOD, Only: ESMF_ArrayCreate 
    Use ESMF_MOD, Only: ESMF_CELL_CENTER 
    Use ESMF_MOD, Only: ESMF_DATA_REF 
    Use ESMF_MOD, Only: ESMF_GridGetDELocalAI 
    Use ESMF_MOD, Only: ESMF_ArrayGet 
#else 
    USE ESMF_MOD, Only: ESMF_FieldSetDataPointer 
    Use ESMF_MOD_private, Only: Field_SetDistRank 
#endif 
    Implicit None 
    Type (ESMF_Field) :: aField 
    Integer , Pointer           :: aPtr (:,:) 
    Type (ESMF_Grid), Intent(In)          :: aGrid 
    Integer, Intent(In), Optional    :: shape(:) 
    Integer, Intent(In), Optional    :: halo_width 
    Integer, Intent(In), Optional    :: grid_rank_index(:) 
#ifdef USE_ESMF 
    Type (ESMF_RelLoc), Intent(in), Optional    :: horzRelLoc 
#else 
    Integer, Intent(In), Optional :: horzRelLoc 
#endif 
    Character(Len=*), Intent(In), Optional :: name 
    Integer, Intent(Out), Optional   :: rc 
 
    Type (ESMF_AxisIndex) :: AI(N_DIMENSIONS) 
#ifdef USE_ESMF 
    Type (ESMF_FieldDataMap) :: map 
    Type (ESMF_Array) :: tmpArray 
    Type (ESMF_RelLoc) :: hzRelLoc 
#endif 
    Integer            :: hwidth 
    Character(Len=ESMF_MAXSTR) :: name_ 
 
    Integer, Dimension(2) :: l, u ! lower and upper bounds 
    Integer, Dimension(2) :: dataIndices 
    Integer :: j_min, j_max 
    Integer :: stat 
    Integer :: idx, gidx 
    Integer :: j_index, i_index 
 
    !------------------------------------------- 
 
#ifdef USE_ESMF 
      hzRelLoc = ESMF_CELL_CENTER 
      If (Present(horzRelLoc)) hzRelLoc = horzRelLoc 
#endif 
      hwidth   = 1 
      name_    = 'unnamed' 
 
      If (Present(halo_width)) Then 
         ASSERT(halo_width >= 0, 'illegal value for hwidth') 
         hwidth   = halo_width 
      End If 
 
      If (Present(name)) name_ = Trim(name) 
 
    l = 1      ! default lower bound 
    u = shape  ! default upper bound 
 
#ifdef USE_ESMF 
    Call ESMF_GridGetDELocalAI(aGrid, AIPerDim=AI, horzrelloc=hzRelLoc, rc=rc) 
 
    ! Extract bounds of computational domain 
    ! Adjust lower/upper bounds with halo size 
    Do idx = 1, 2 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = AI(gidx)%min - hwidth ! Would like hwidth to be independent of axis, but ESMF.1.0.7 
          u(idx) = AI(gidx)%max + hwidth ! does not support this. 
       End If 
    End Do 
#else 
    Do idx = 1, 2 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = l(idx)-hwidth 
          u(idx) = u(idx)+hwidth 
       End If 
    End Do 
#endif 
 
 
    Allocate(aPtr(l(1):u(1),l(2):u(2)), STAT=stat) 
    ASSERT_ALWAYS(stat == 0, 'Allocate Failure') 
 
#ifdef USE_ESMF 
    Call ESMF_FieldDataMapSet(map, 2, grid_rank_index, shape, horzRelloc=hzRelLoc, rc=rc) 
    tmpArray = ESMF_ArrayCreate(aPtr, ESMF_DATA_REF, haloWidth=hwidth, rc=rc) 
    Call ESMF_ArrayGet(tmparray, haloWidth=hwidth) 
    aField   = ESMF_FieldCreate(aGrid,  tmpArray, halowidth=hwidth, horzRelloc=hzRelLoc, name=name_, datamap=map, rc=rc) 
#else 
 
    Call ESMF_FieldSetDataPointer(aField,aPtr) 
    Do idx = 1, 2 
       gidx =grid_rank_index(idx) 
       If (gidx > 0) Exit 
    End Do 
 
    Call Field_SetDistRank(aField,idx) 
#endif 
 
 
    rc = ESMF_SUCCESS 
 
  End Function Field_Allocate_Integer_2D 

  
  Function Field_Allocate_Integer_3D ( & 
    &     aPtr, aGrid, shape, halo_width, grid_rank_index, horzRelLoc, name, rc) Result(aField) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_MAXSTR 
#ifdef USE_ESMF 
    Use ESMF_MOD, Only: ESMF_FieldDataMap 
    Use ESMF_MOD, Only: ESMF_RelLoc 
    Use ESMF_MOD, Only: ESMF_FieldCreate 
    Use ESMF_MOD, Only: ESMF_FieldDataMapSet 
    Use ESMF_MOD, Only: ESMF_Array 
    Use ESMF_MOD, Only: ESMF_ArrayCreate 
    Use ESMF_MOD, Only: ESMF_CELL_CENTER 
    Use ESMF_MOD, Only: ESMF_DATA_REF 
    Use ESMF_MOD, Only: ESMF_GridGetDELocalAI 
    Use ESMF_MOD, Only: ESMF_ArrayGet 
#else 
    USE ESMF_MOD, Only: ESMF_FieldSetDataPointer 
    Use ESMF_MOD_private, Only: Field_SetDistRank 
#endif 
    Implicit None 
    Type (ESMF_Field) :: aField 
    Integer , Pointer           :: aPtr (:,:,:) 
    Type (ESMF_Grid), Intent(In)          :: aGrid 
    Integer, Intent(In), Optional    :: shape(:) 
    Integer, Intent(In), Optional    :: halo_width 
    Integer, Intent(In), Optional    :: grid_rank_index(:) 
#ifdef USE_ESMF 
    Type (ESMF_RelLoc), Intent(in), Optional    :: horzRelLoc 
#else 
    Integer, Intent(In), Optional :: horzRelLoc 
#endif 
    Character(Len=*), Intent(In), Optional :: name 
    Integer, Intent(Out), Optional   :: rc 
 
    Type (ESMF_AxisIndex) :: AI(N_DIMENSIONS) 
#ifdef USE_ESMF 
    Type (ESMF_FieldDataMap) :: map 
    Type (ESMF_Array) :: tmpArray 
    Type (ESMF_RelLoc) :: hzRelLoc 
#endif 
    Integer            :: hwidth 
    Character(Len=ESMF_MAXSTR) :: name_ 
 
    Integer, Dimension(3) :: l, u ! lower and upper bounds 
    Integer, Dimension(3) :: dataIndices 
    Integer :: j_min, j_max 
    Integer :: stat 
    Integer :: idx, gidx 
    Integer :: j_index, i_index 
 
    !------------------------------------------- 
 
#ifdef USE_ESMF 
      hzRelLoc = ESMF_CELL_CENTER 
      If (Present(horzRelLoc)) hzRelLoc = horzRelLoc 
#endif 
      hwidth   = 1 
      name_    = 'unnamed' 
 
      If (Present(halo_width)) Then 
         ASSERT(halo_width >= 0, 'illegal value for hwidth') 
         hwidth   = halo_width 
      End If 
 
      If (Present(name)) name_ = Trim(name) 
 
    l = 1      ! default lower bound 
    u = shape  ! default upper bound 
 
#ifdef USE_ESMF 
    Call ESMF_GridGetDELocalAI(aGrid, AIPerDim=AI, horzrelloc=hzRelLoc, rc=rc) 
 
    ! Extract bounds of computational domain 
    ! Adjust lower/upper bounds with halo size 
    Do idx = 1, 3 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = AI(gidx)%min - hwidth ! Would like hwidth to be independent of axis, but ESMF.1.0.7 
          u(idx) = AI(gidx)%max + hwidth ! does not support this. 
       End If 
    End Do 
#else 
    Do idx = 1, 3 
       gidx = grid_rank_index(idx) 
       If (gidx > 0) Then  ! halo for grid 
          l(idx) = l(idx)-hwidth 
          u(idx) = u(idx)+hwidth 
       End If 
    End Do 
#endif 
 
 
    Allocate(aPtr(l(1):u(1),l(2):u(2),l(3):u(3)), STAT=stat) 
    ASSERT_ALWAYS(stat == 0, 'Allocate Failure') 
 
#ifdef USE_ESMF 
    Call ESMF_FieldDataMapSet(map, 3, grid_rank_index, shape, horzRelloc=hzRelLoc, rc=rc) 
    tmpArray = ESMF_ArrayCreate(aPtr, ESMF_DATA_REF, haloWidth=hwidth, rc=rc) 
    Call ESMF_ArrayGet(tmparray, haloWidth=hwidth) 
    aField   = ESMF_FieldCreate(aGrid,  tmpArray, halowidth=hwidth, horzRelloc=hzRelLoc, name=name_, datamap=map, rc=rc) 
#else 
 
    Call ESMF_FieldSetDataPointer(aField,aPtr) 
    Do idx = 1, 3 
       gidx =grid_rank_index(idx) 
       If (gidx > 0) Exit 
    End Do 
 
    Call Field_SetDistRank(aField,idx) 
#endif 
 
 
    rc = ESMF_SUCCESS 
 
  End Function Field_Allocate_Integer_3D 

  
  Subroutine Field_Deallocate_Real4_3D ( aField, aPtr) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_FieldDestroy 
    Implicit None 
    Type (ESMF_Field), Intent(Out)        :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer :: aPtr (:,:,:) 
    !------------------------------------------- 
    Integer :: rc 
 
    Call ESMF_FieldDestroy(aField, rc=rc) 
    Nullify(aPtr) 
 
  End Subroutine Field_Deallocate_Real4_3D 

  
  Subroutine Field_Deallocate_Real8_1D ( aField, aPtr) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_FieldDestroy 
    Implicit None 
    Type (ESMF_Field), Intent(Out)        :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer :: aPtr (:) 
    !------------------------------------------- 
    Integer :: rc 
 
    Call ESMF_FieldDestroy(aField, rc=rc) 
    Nullify(aPtr) 
 
  End Subroutine Field_Deallocate_Real8_1D 

  
  Subroutine Field_Deallocate_Real8_2D ( aField, aPtr) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_FieldDestroy 
    Implicit None 
    Type (ESMF_Field), Intent(Out)        :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer :: aPtr (:,:) 
    !------------------------------------------- 
    Integer :: rc 
 
    Call ESMF_FieldDestroy(aField, rc=rc) 
    Nullify(aPtr) 
 
  End Subroutine Field_Deallocate_Real8_2D 

  
  Subroutine Field_Deallocate_Real8_3D ( aField, aPtr) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_FieldDestroy 
    Implicit None 
    Type (ESMF_Field), Intent(Out)        :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer :: aPtr (:,:,:) 
    !------------------------------------------- 
    Integer :: rc 
 
    Call ESMF_FieldDestroy(aField, rc=rc) 
    Nullify(aPtr) 
 
  End Subroutine Field_Deallocate_Real8_3D 

  
  Subroutine Field_Deallocate_Real8_4D ( aField, aPtr) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_FieldDestroy 
    Implicit None 
    Type (ESMF_Field), Intent(Out)        :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer :: aPtr (:,:,:,:) 
    !------------------------------------------- 
    Integer :: rc 
 
    Call ESMF_FieldDestroy(aField, rc=rc) 
    Nullify(aPtr) 
 
  End Subroutine Field_Deallocate_Real8_4D 

  
  Subroutine Field_Deallocate_Real8_5D ( aField, aPtr) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_FieldDestroy 
    Implicit None 
    Type (ESMF_Field), Intent(Out)        :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer :: aPtr (:,:,:,:,:) 
    !------------------------------------------- 
    Integer :: rc 
 
    Call ESMF_FieldDestroy(aField, rc=rc) 
    Nullify(aPtr) 
 
  End Subroutine Field_Deallocate_Real8_5D 

  
  Subroutine Field_Deallocate_Integer_1D ( aField, aPtr) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_FieldDestroy 
    Implicit None 
    Type (ESMF_Field), Intent(Out)        :: aField 
    Integer , Pointer :: aPtr (:) 
    !------------------------------------------- 
    Integer :: rc 
 
    Call ESMF_FieldDestroy(aField, rc=rc) 
    Nullify(aPtr) 
 
  End Subroutine Field_Deallocate_Integer_1D 

  
  Subroutine Field_Deallocate_Integer_2D ( aField, aPtr) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_FieldDestroy 
    Implicit None 
    Type (ESMF_Field), Intent(Out)        :: aField 
    Integer , Pointer :: aPtr (:,:) 
    !------------------------------------------- 
    Integer :: rc 
 
    Call ESMF_FieldDestroy(aField, rc=rc) 
    Nullify(aPtr) 
 
  End Subroutine Field_Deallocate_Integer_2D 

  
  Subroutine Field_Deallocate_Integer_3D ( aField, aPtr) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_FieldDestroy 
    Implicit None 
    Type (ESMF_Field), Intent(Out)        :: aField 
    Integer , Pointer :: aPtr (:,:,:) 
    !------------------------------------------- 
    Integer :: rc 
 
    Call ESMF_FieldDestroy(aField, rc=rc) 
    Nullify(aPtr) 
 
  End Subroutine Field_Deallocate_Integer_3D 

  Function Field_Clone(aField, new_name) Result(aClone)
  Use ESMF_MOD, Only: ESMF_Field
  Use ESMF_MOD, Only: ESMF_Array
  Use ESMF_MOD, Only: ESMF_Grid
  Use ESMF_MOD, Only: ESMF_MAXSTR
  Use ESMF_MOD, Only: ESMF_DataType
  Use ESMF_MOD, Only: ESMF_DataKind
#ifdef USE_ESMF
  Use ESMF_MOD, Only: ESMF_FieldDataMap
  Use ESMF_MOD, Only: ESMF_RelLoc
  Use ESMF_MOD, Only: ESMF_FieldGet
  Use ESMF_MOD, Only: ESMF_ArrayGet
  Use ESMF_MOD, Only: ESMF_ArrayCreate
  Use ESMF_MOD, Only: ESMF_FieldCreate
  Use ESMF_MOD, Only: ESMF_ArrayDestroy
#else
  USE ESMF_MOD_private, Only : Clone => Field_Clone
#endif
  Implicit None
  Type (ESMF_Field) :: aField
  Character(Len=*), Optional, Intent(In) :: new_name
  Type (ESMF_Field) :: aClone

  Integer :: rank
  Type (ESMF_DataType) :: type
  Type (ESMF_DataKind) :: kind

#ifdef USE_ESMF
  Integer :: haloWidth
  Character(Len=ESMF_MAXSTR) :: orig_name, name_
  Integer :: rc
  Integer, Parameter :: MAX_RANK = 7
  Integer, Dimension(MAX_RANK) :: lbounds, ubounds, counts
  Type (ESMF_Grid) :: grid
  Type (ESMF_Array) :: array, new_array
  Type (ESMF_FieldDataMap) :: datamap
  Type (ESMF_RelLoc) :: horzrelloc, vertrelloc


  Call ESMF_FieldGet(aField, grid, array, datamap=datamap, &
    & horzrelloc=horzrelloc, vertrelloc=vertrelloc, name=orig_name, rc=rc)

  name_ = 'clone of '//orig_name
  If (Present(new_name)) name_ = new_name

  halowidth=1
  Call ESMF_ArrayGet(array, rank, type, kind, counts=counts, haloWidth=halowidth, &
   & lbounds=lbounds, ubounds=ubounds, rc=rc)
  print*,'haloWidth = ',halowidth,' should be 1'
  halowidth = 1
  new_array = ESMF_ArrayCreate(rank, type, kind, counts=counts, haloWidth=haloWidth, &
     & lbounds=lbounds, ubounds=ubounds, rc=rc)
  aClone = ESMF_FieldCreate(grid, new_array, horzrelloc=horzrelloc, &
    & name=name_, haloWidth=haloWidth, datamap=datamap, rc=rc)
  !!$$Call ESMF_ArrayDestroy(new_array, rc)
#else
  aClone = Clone(aField)
#endif
  End Function Field_Clone


End Module ESMF_CUSTOM_MOD
