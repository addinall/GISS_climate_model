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
Module ESMF_MOD_private
  Implicit None
  Private

  ! Derived Types
  Public :: Grid
  Public :: Array
  Public :: Field
  Public :: AxisIndex
  Public :: RouteHandle

  ! Methods
  Public :: Field_SetDataPointer
  Public :: Field_GetDataPointer
  Public :: Field_GetArray
  Public :: Field_SetDistRank, Field_GetDistRank
  Public :: Field_GetRank, Field_GetShape
  Public :: Field_Clone
  Public :: Field_Destroy
  Public :: Grid_Create
  Public :: Grid_Destroy
  Public :: Grid_Get

  ! Constants
  Public :: ESMF_KIND_R4
  Public :: ESMF_KIND_R8
  Integer, Parameter :: ESMF_KIND_R4 = Selected_Real_Kind(P=6)
  Integer, Parameter :: ESMF_KIND_R8 = Selected_Real_Kind(P=14)

  ! Constants
  Public :: ESMF_FAILURE
  Public :: ESMF_SUCCESS
  Integer, Parameter :: ESMF_FAILURE = -1
  Integer, Parameter :: ESMF_SUCCESS = 0

  Public :: TYPE_REAL
  Public :: TYPE_INTEGER
  Integer, Parameter :: TYPE_REAL = 1
  Integer, Parameter :: TYPE_INTEGER = 2

  Public :: ESMF_CELL_NFACE
  Public :: ESMF_CELL_SFACE
  Public :: ESMF_CELL_CENTER
  Public :: ESMF_CELL_SWCORNER

  Integer, Parameter :: ESMF_CELL_NFACE=1
  Integer, Parameter :: ESMF_CELL_SFACE=2
  Integer, Parameter :: ESMF_CELL_CENTER=3
  Integer, Parameter :: ESMF_CELL_SWCORNER=4

  Public :: MAXSTR
  Integer, Parameter :: MAXSTR = 100
  Integer, Parameter :: N_DIMENSIONS = 3 ! dimensions of grid

  Public :: ESMF_DATATYPE, ESMF_DATAKIND
  Public :: ESMF_DATA_REAL, ESMF_DATA_INTEGER
  Public :: ESMF_R8, ESMF_R4
  Public :: ESMF_I8, ESMF_I4

  Public :: Operator(==)
  interface operator (==)
     module procedure ESMF_type_eq
     module procedure ESMF_kind_eq
  End interface

  Type ESMF_DATATYPE
     Integer :: type
  End Type ESMF_DATATYPE

  Type ESMF_DATAKIND
     Integer :: kind
  End Type ESMF_DATAKIND

  Type (ESMF_DATATYPE), Parameter :: ESMF_DATA_REAL=ESMF_DATATYPE(1)
  Type (ESMF_DATATYPE), Parameter :: ESMF_DATA_INTEGER=ESMF_DATATYPE(2)
  Type (ESMF_DATAKIND), Parameter :: ESMF_R8=ESMF_DATAKIND(1)
  Type (ESMF_DATAKIND), Parameter :: ESMF_R4=ESMF_DATAKIND(2)
  Type (ESMF_DATAKIND), Parameter :: ESMF_I8=ESMF_DATAKIND(3)
  Type (ESMF_DATAKIND), Parameter :: ESMF_I4=ESMF_DATAKIND(4)

  type AxisIndex
     sequence
     integer :: min
     integer :: max
     integer :: stride
  end type AxisIndex

  Type Grid
     Private
     Integer :: IM
     Integer :: JM
     Integer :: LM
     Type (AxisIndex) :: global(N_DIMENSIONS)
  End Type Grid

!!$  Type HaloDirection
!!$     Private
!!$     Integer :: direction ! NORTH, SOUTH, EAST, WEST
!!$  End Type HaloDirection
!!$  Integer, Public, Parameter :: NORTH = 2**0
!!$  Integer, Public, Parameter :: SOUTH = 2**1
  Integer, Parameter :: MAX_RANK=5
  Type Field
     Private

     Logical :: active = .false.
     Integer :: esmf_kind = -1
     Integer :: esmf_type = -2
     Integer :: rank, dist_rank
     Integer :: data_shape(MAX_RANK)
     Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_1D (:) => Null()
     Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_2D (:,:) => Null()
     Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_3D (:,:,:) => Null()
     Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_4D (:,:,:,:) => Null()
     Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_5D (:,:,:,:,:) => Null()
     Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_1D (:) => Null()
     Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_2D (:,:) => Null()
     Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_3D (:,:,:) => Null()
     Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_4D (:,:,:,:) => Null()
     Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_5D (:,:,:,:,:) => Null()
     Integer , Pointer :: ptr_Integer_1D (:) => Null()
     Integer , Pointer :: ptr_Integer_2D (:,:) => Null()
     Integer , Pointer :: ptr_Integer_3D (:,:,:) => Null()
     Integer , Pointer :: ptr_Integer_4D (:,:,:,:) => Null()
     Integer , Pointer :: ptr_Integer_5D (:,:,:,:,:) => Null()
  End Type Field

  Type Array
     private
     Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_1D (:) => Null()
     Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_2D (:,:) => Null()
     Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_3D (:,:,:) => Null()
     Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_4D (:,:,:,:) => Null()
     Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_5D (:,:,:,:,:) => Null()
     Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_1D (:) => Null()
     Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_2D (:,:) => Null()
     Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_3D (:,:,:) => Null()
     Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_4D (:,:,:,:) => Null()
     Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_5D (:,:,:,:,:) => Null()
     Integer , Pointer :: ptr_Integer_1D (:) => Null()
     Integer , Pointer :: ptr_Integer_2D (:,:) => Null()
     Integer , Pointer :: ptr_Integer_3D (:,:,:) => Null()
     Integer , Pointer :: ptr_Integer_4D (:,:,:,:) => Null()
     Integer , Pointer :: ptr_Integer_5D (:,:,:,:,:) => Null()
  End Type Array

  Type RouteHandle
     Private
     Integer :: cannot_have_empty_types
  End Type RouteHandle



  Interface Field_SetDataPointer
     Module Procedure Field_SetDataPointer_Real4_1D
     Module Procedure Field_SetDataPointer_Real4_2D
     Module Procedure Field_SetDataPointer_Real4_3D
     Module Procedure Field_SetDataPointer_Real4_4D
     Module Procedure Field_SetDataPointer_Real4_5D
     Module Procedure Field_SetDataPointer_Real8_1D
     Module Procedure Field_SetDataPointer_Real8_2D
     Module Procedure Field_SetDataPointer_Real8_3D
     Module Procedure Field_SetDataPointer_Real8_4D
     Module Procedure Field_SetDataPointer_Real8_5D
     Module Procedure Field_SetDataPointer_Integer_1D
     Module Procedure Field_SetDataPointer_Integer_2D
     Module Procedure Field_SetDataPointer_Integer_3D
     Module Procedure Field_SetDataPointer_Integer_4D
     Module Procedure Field_SetDataPointer_Integer_5D
  End Interface



  Interface Field_GetDataPointer
     Module Procedure Field_GetDataPointer_Real4_1D
     Module Procedure Field_GetDataPointer_Real4_2D
     Module Procedure Field_GetDataPointer_Real4_3D
     Module Procedure Field_GetDataPointer_Real4_4D
     Module Procedure Field_GetDataPointer_Real4_5D
     Module Procedure Field_GetDataPointer_Real8_1D
     Module Procedure Field_GetDataPointer_Real8_2D
     Module Procedure Field_GetDataPointer_Real8_3D
     Module Procedure Field_GetDataPointer_Real8_4D
     Module Procedure Field_GetDataPointer_Real8_5D
     Module Procedure Field_GetDataPointer_Integer_1D
     Module Procedure Field_GetDataPointer_Integer_2D
     Module Procedure Field_GetDataPointer_Integer_3D
     Module Procedure Field_GetDataPointer_Integer_4D
     Module Procedure Field_GetDataPointer_Integer_5D
  End Interface

Contains

  ! This routine is a simplified rendition of DataMap. It
  ! is used simply to specify which rank of a field is decomposed.
  Subroutine Field_SetDistRank(aField, dist_rank)
    Type (Field), Intent(InOut) :: aField
    Integer, Intent(In) :: dist_rank

    aField%dist_rank = dist_rank
  End Subroutine Field_SetDistRank

  ! This routine is a simplified rendition of DataMap. It
  ! is used simply to determine which rank of a field is decomposed.
  Function Field_GetDistRank(aField) Result(dist_rank)
    Type (Field), Intent(In) :: aField
    Integer :: dist_rank

    dist_rank = aField%dist_rank
  End Function Field_GetDistRank

  Function Field_GetRank(aField) Result(rank)
    Type (Field), Intent(In) :: aField
    Integer :: rank

    rank = aField%rank
  End Function Field_GetRank

  Subroutine Field_GetShape(aField, data_shape)
    Type (Field), Intent(In) :: aField
    Integer :: data_shape(aField%rank)

    data_shape = aField%data_shape(1:aField%rank)

  End Subroutine Field_GetShape
  
  Subroutine Field_SetDataPointer_Real4_1D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer       :: dataPointer (:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 1 
    aField% ptr_Real4_1D => dataPointer 
    aField%data_shape(1:1) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Real4_1D 

  
  Subroutine Field_SetDataPointer_Real4_2D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer       :: dataPointer (:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 2 
    aField% ptr_Real4_2D => dataPointer 
    aField%data_shape(1:2) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Real4_2D 

  
  Subroutine Field_SetDataPointer_Real4_3D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer       :: dataPointer (:,:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 3 
    aField% ptr_Real4_3D => dataPointer 
    aField%data_shape(1:3) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Real4_3D 

  
  Subroutine Field_SetDataPointer_Real4_4D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer       :: dataPointer (:,:,:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 4 
    aField% ptr_Real4_4D => dataPointer 
    aField%data_shape(1:4) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Real4_4D 

  
  Subroutine Field_SetDataPointer_Real4_5D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer       :: dataPointer (:,:,:,:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 5 
    aField% ptr_Real4_5D => dataPointer 
    aField%data_shape(1:5) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Real4_5D 

  
  Subroutine Field_SetDataPointer_Real8_1D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer       :: dataPointer (:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 1 
    aField% ptr_Real8_1D => dataPointer 
    aField%data_shape(1:1) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Real8_1D 

  
  Subroutine Field_SetDataPointer_Real8_2D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer       :: dataPointer (:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 2 
    aField% ptr_Real8_2D => dataPointer 
    aField%data_shape(1:2) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Real8_2D 

  
  Subroutine Field_SetDataPointer_Real8_3D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer       :: dataPointer (:,:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 3 
    aField% ptr_Real8_3D => dataPointer 
    aField%data_shape(1:3) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Real8_3D 

  
  Subroutine Field_SetDataPointer_Real8_4D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer       :: dataPointer (:,:,:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 4 
    aField% ptr_Real8_4D => dataPointer 
    aField%data_shape(1:4) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Real8_4D 

  
  Subroutine Field_SetDataPointer_Real8_5D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer       :: dataPointer (:,:,:,:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 5 
    aField% ptr_Real8_5D => dataPointer 
    aField%data_shape(1:5) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Real8_5D 

  
  Subroutine Field_SetDataPointer_Integer_1D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Integer , Pointer       :: dataPointer (:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 1 
    aField% ptr_Integer_1D => dataPointer 
    aField%data_shape(1:1) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Integer_1D 

  
  Subroutine Field_SetDataPointer_Integer_2D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Integer , Pointer       :: dataPointer (:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 2 
    aField% ptr_Integer_2D => dataPointer 
    aField%data_shape(1:2) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Integer_2D 

  
  Subroutine Field_SetDataPointer_Integer_3D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Integer , Pointer       :: dataPointer (:,:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 3 
    aField% ptr_Integer_3D => dataPointer 
    aField%data_shape(1:3) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Integer_3D 

  
  Subroutine Field_SetDataPointer_Integer_4D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Integer , Pointer       :: dataPointer (:,:,:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 4 
    aField% ptr_Integer_4D => dataPointer 
    aField%data_shape(1:4) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Integer_4D 

  
  Subroutine Field_SetDataPointer_Integer_5D (aField, dataPointer, haloWidth, rc) 
    Implicit None 
    Type (Field), Intent(InOut) :: aField 
    Integer , Pointer       :: dataPointer (:,:,:,:,:) 
    Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
    Integer, Optional, Intent(out) :: rc 
     
    aField%rank = 5 
    aField% ptr_Integer_5D => dataPointer 
    aField%data_shape(1:5) = shape(dataPointer) 
 
    afield%active = .true. 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_SetDataPointer_Integer_5D 

  
  Subroutine Field_GetDataPointer_Real4_1D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer           :: ptr (:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Real4_1D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Real4_1D 

  
  Subroutine Field_GetDataPointer_Real4_2D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer           :: ptr (:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Real4_2D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Real4_2D 

  
  Subroutine Field_GetDataPointer_Real4_3D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer           :: ptr (:,:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Real4_3D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Real4_3D 

  
  Subroutine Field_GetDataPointer_Real4_4D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer           :: ptr (:,:,:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Real4_4D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Real4_4D 

  
  Subroutine Field_GetDataPointer_Real4_5D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Real (Kind=ESMF_KIND_R4), Pointer           :: ptr (:,:,:,:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Real4_5D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Real4_5D 

  
  Subroutine Field_GetDataPointer_Real8_1D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer           :: ptr (:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Real8_1D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Real8_1D 

  
  Subroutine Field_GetDataPointer_Real8_2D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer           :: ptr (:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Real8_2D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Real8_2D 

  
  Subroutine Field_GetDataPointer_Real8_3D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer           :: ptr (:,:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Real8_3D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Real8_3D 

  
  Subroutine Field_GetDataPointer_Real8_4D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer           :: ptr (:,:,:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Real8_4D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Real8_4D 

  
  Subroutine Field_GetDataPointer_Real8_5D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Real (Kind=ESMF_KIND_R8), Pointer           :: ptr (:,:,:,:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Real8_5D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Real8_5D 

  
  Subroutine Field_GetDataPointer_Integer_1D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Integer , Pointer           :: ptr (:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Integer_1D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Integer_1D 

  
  Subroutine Field_GetDataPointer_Integer_2D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Integer , Pointer           :: ptr (:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Integer_2D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Integer_2D 

  
  Subroutine Field_GetDataPointer_Integer_3D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Integer , Pointer           :: ptr (:,:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Integer_3D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Integer_3D 

  
  Subroutine Field_GetDataPointer_Integer_4D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Integer , Pointer           :: ptr (:,:,:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Integer_4D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Integer_4D 

  
  Subroutine Field_GetDataPointer_Integer_5D (aField, ptr, rc) 
    Implicit None 
    Type (Field), Intent(In) :: aField 
    Integer , Pointer           :: ptr (:,:,:,:,:) 
    integer, intent(out), optional :: rc  
     
    ptr => aField% ptr_Integer_5D 
 
    If (Present(rc)) rc = ESMF_SUCCESS 
     
  End Subroutine Field_GetDataPointer_Integer_5D 

  Function Grid_Create(im, jm, lm) Result(aGrid)
    Implicit None
    Integer :: im, jm, lm
    Type (Grid) :: aGrid

    aGrid%im = im
    aGrid%jm = jm
    aGrid%lm = lm

  End Function Grid_Create

  Subroutine Grid_Destroy(agrid, rc)
    Type (Grid) :: agrid
    Integer, Optional, Intent(Out) :: rc

    If (present(rc)) rc = Esmf_Success
  End Subroutine Grid_Destroy

  Subroutine Grid_Get(aGrid, ai)
    Implicit None
    Type (Grid) :: aGrid
    Type (AxisIndex) :: ai(N_DIMENSIONS)

  !!$$ ai%min = (/ 1, 1, 1 /)
  !!$$ ai%max = (/ aGrid%im, aGrid%jm, aGrid%lm /)
    ai%min = (/ 1, 1, 1 /)
    ai%max = (/ aGrid%im, aGrid%jm, aGrid%lm /)

  End Subroutine Grid_Get

  !< Method: [[Field_HaloStore]] >>
  !< Method: [[Field_Halo]] >>

  Subroutine Field_GetArray(aField, anArray, rc)
    Type (Field) :: aField
    Type (Array) :: anArray
    Integer, Optional, Intent(Out) :: rc

    ! Not implemented
    If (Present(rc)) rc = ESMF_FAILURE

  End Subroutine Field_GetArray

  Function Field_Clone(f1) Result(f2)
    Type (Field), Intent(In), Target :: f1
    Type (Field) :: f2
    Integer, Pointer :: s(:)

    f2%active = f1%active
    f2%esmf_kind = f1%esmf_kind
    f2%esmf_type = f1%esmf_type
    f2%rank = f1%rank
    f2%dist_rank = f1%dist_rank
    f2%data_shape = f1%data_shape

    s => f1%data_shape(:f1%rank)
     
    If (Associated(f1%ptr_Real4_1D)) Then 
      Allocate(f2%ptr_Real4_1D(s(1))) 
    End If 

     
    If (Associated(f1%ptr_Real4_2D)) Then 
      Allocate(f2%ptr_Real4_2D(s(1),s(2))) 
    End If 

     
    If (Associated(f1%ptr_Real4_3D)) Then 
      Allocate(f2%ptr_Real4_3D(s(1),s(2),s(3))) 
    End If 

     
    If (Associated(f1%ptr_Real4_4D)) Then 
      Allocate(f2%ptr_Real4_4D(s(1),s(2),s(3),s(4))) 
    End If 

     
    If (Associated(f1%ptr_Real4_5D)) Then 
      Allocate(f2%ptr_Real4_5D(s(1),s(2),s(3),s(4),s(5))) 
    End If 

     
    If (Associated(f1%ptr_Real8_1D)) Then 
      Allocate(f2%ptr_Real8_1D(s(1))) 
    End If 

     
    If (Associated(f1%ptr_Real8_2D)) Then 
      Allocate(f2%ptr_Real8_2D(s(1),s(2))) 
    End If 

     
    If (Associated(f1%ptr_Real8_3D)) Then 
      Allocate(f2%ptr_Real8_3D(s(1),s(2),s(3))) 
    End If 

     
    If (Associated(f1%ptr_Real8_4D)) Then 
      Allocate(f2%ptr_Real8_4D(s(1),s(2),s(3),s(4))) 
    End If 

     
    If (Associated(f1%ptr_Real8_5D)) Then 
      Allocate(f2%ptr_Real8_5D(s(1),s(2),s(3),s(4),s(5))) 
    End If 

     
    If (Associated(f1%ptr_Integer_1D)) Then 
      Allocate(f2%ptr_Integer_1D(s(1))) 
    End If 

     
    If (Associated(f1%ptr_Integer_2D)) Then 
      Allocate(f2%ptr_Integer_2D(s(1),s(2))) 
    End If 

     
    If (Associated(f1%ptr_Integer_3D)) Then 
      Allocate(f2%ptr_Integer_3D(s(1),s(2),s(3))) 
    End If 

     
    If (Associated(f1%ptr_Integer_4D)) Then 
      Allocate(f2%ptr_Integer_4D(s(1),s(2),s(3),s(4))) 
    End If 

     
    If (Associated(f1%ptr_Integer_5D)) Then 
      Allocate(f2%ptr_Integer_5D(s(1),s(2),s(3),s(4),s(5))) 
    End If 

  End Function Field_Clone

  Subroutine Field_Destroy(f, rc)
    Type (Field), Intent(InOut), Target :: f
    Integer, Optional, Intent(Out) :: rc

    f%active = .false.
    f%rank = 0
    f%dist_rank = 0
    f%data_shape = 0

    If (Present(rc)) rc = ESMF_SUCCESS
     
    If (Associated(f%ptr_Real4_1D)) Then 
       Deallocate(f%ptr_Real4_1D) 
    End If 

     
    If (Associated(f%ptr_Real4_2D)) Then 
       Deallocate(f%ptr_Real4_2D) 
    End If 

     
    If (Associated(f%ptr_Real4_3D)) Then 
       Deallocate(f%ptr_Real4_3D) 
    End If 

     
    If (Associated(f%ptr_Real4_4D)) Then 
       Deallocate(f%ptr_Real4_4D) 
    End If 

     
    If (Associated(f%ptr_Real4_5D)) Then 
       Deallocate(f%ptr_Real4_5D) 
    End If 

     
    If (Associated(f%ptr_Real8_1D)) Then 
       Deallocate(f%ptr_Real8_1D) 
    End If 

     
    If (Associated(f%ptr_Real8_2D)) Then 
       Deallocate(f%ptr_Real8_2D) 
    End If 

     
    If (Associated(f%ptr_Real8_3D)) Then 
       Deallocate(f%ptr_Real8_3D) 
    End If 

     
    If (Associated(f%ptr_Real8_4D)) Then 
       Deallocate(f%ptr_Real8_4D) 
    End If 

     
    If (Associated(f%ptr_Real8_5D)) Then 
       Deallocate(f%ptr_Real8_5D) 
    End If 

     
    If (Associated(f%ptr_Integer_1D)) Then 
       Deallocate(f%ptr_Integer_1D) 
    End If 

     
    If (Associated(f%ptr_Integer_2D)) Then 
       Deallocate(f%ptr_Integer_2D) 
    End If 

     
    If (Associated(f%ptr_Integer_3D)) Then 
       Deallocate(f%ptr_Integer_3D) 
    End If 

     
    If (Associated(f%ptr_Integer_4D)) Then 
       Deallocate(f%ptr_Integer_4D) 
    End If 

     
    If (Associated(f%ptr_Integer_5D)) Then 
       Deallocate(f%ptr_Integer_5D) 
    End If 

   End Subroutine Field_Destroy

   Logical Function ESMF_type_eq(t1, t2)
     Type (ESMF_DATATYPE), Intent(In) :: t1, t2
     ESMF_type_eq = (t1%type == t2%type)
   End Function ESMF_type_eq

   Logical Function ESMF_kind_eq(t1, t2)
     Type (ESMF_DATAKIND), Intent(In) :: t1, t2
     ESMF_kind_eq = (t1%kind == t2%kind)
   End Function ESMF_kind_eq

End Module ESMF_MOD_private

#ifndef USE_ESMF
Module ESMF_MOD
  Use ESMF_MOD_private, Only: ESMF_Grid => Grid
  Use ESMF_MOD_private, Only: ESMF_Field => Field
  Use ESMF_MOD_private, Only: ESMF_Array => Array
  Use ESMF_MOD_private, Only: ESMF_AxisIndex => AxisIndex
  Use ESMF_MOD_private, Only: ESMF_KIND_R4
  Use ESMF_MOD_private, Only: ESMF_KIND_R8

  Use ESMF_MOD_private, Only: ESMF_FAILURE, ESMF_SUCCESS
  Use ESMF_MOD_private, Only: ESMF_CELL_NFACE, ESMF_CELL_SFACE
  Use ESMF_MOD_private, Only: ESMF_CELL_CENTER, ESMF_CELL_SWCORNER

  Use ESMF_MOD_private, Only: ESMF_MAXSTR => MAXSTR

  Use ESMF_MOD_private, Only: ESMF_FieldGetDataPointer => Field_GetDataPointer
  Use ESMF_MOD_private, Only: ESMF_FieldGetArray => Field_GetArray
  Use ESMF_MOD_private, Only: ESMF_FieldDestroy => Field_Destroy

  Use ESMF_MOD_private, Only: ESMF_DATATYPE, ESMF_DATAKIND
  Use ESMF_MOD_private, Only: ESMF_DATA_INTEGER, ESMF_DATA_REAL
  Use ESMF_MOD_private, Only: ESMF_R8, ESMF_R4
  Use ESMF_MOD_private, Only: ESMF_I8, ESMF_I4

  Use ESMF_Mod_private, Only: Operator(==)

  Use ESMF_Mod_private, Only: ESMF_GridCreate => Grid_Create
  Use ESMF_Mod_private, Only: ESMF_GridDestroy => Grid_Destroy
  Use ESMF_MOD_private, Only: ESMF_FieldSetDataPointer => Field_SetDataPointer

  Implicit None
  Private

  Public :: ESMF_Grid
  Public :: ESMF_Array
  Public :: ESMF_Field
  Public :: ESMF_AxisIndex

  Public :: ESMF_KIND_R4
  Public :: ESMF_KIND_R8

  Public :: ESMF_FieldGetDataPointer
  Public :: ESMF_FieldSetDataPointer
  Public :: ESMF_FieldGetArray
  Public :: ESMF_FieldDestroy

  Public :: ESMF_FAILURE, ESMF_SUCCESS

  Public :: ESMF_CELL_NFACE
  Public :: ESMF_CELL_SFACE
  Public :: ESMF_CELL_CENTER
  Public :: ESMF_CELL_SWCORNER

  Public :: ESMF_MAXSTR

  Public :: ESMF_DATATYPE, ESMF_DATAKIND
  Public :: ESMF_DATA_INTEGER, ESMF_DATA_REAL
  Public :: ESMF_R8, ESMF_R4
  Public :: ESMF_I8, ESMF_I4

  Public :: Operator(==)

  Public :: ESMF_Initialize, ESMF_Finalize
  Public :: ESMF_GridCreate, ESMF_GridDestroy

Contains

  Subroutine ESMF_Initialize(rc)
    Integer, Optional, Intent(Out) :: rc
    If (Present(rc)) rc = ESMF_SUCCESS
  End Subroutine ESMF_Initialize

  Subroutine ESMF_Finalize(rc)
    Integer, Optional, Intent(Out) :: rc
    If (Present(rc)) rc = ESMF_SUCCESS
  End Subroutine ESMF_Finalize

End Module ESMF_MOD
#endif
