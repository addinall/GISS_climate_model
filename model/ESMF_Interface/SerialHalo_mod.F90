module SerialHalo_mod
  implicit none
  private

  public :: Halo_type
  public :: update
  
  public :: NORTH, SOUTH, EAST, WEST
  type Halo_type
    integer :: emptyPlaceholder
  end type Halo_type

  integer, parameter :: NORTH = 2^0
  integer, parameter :: SOUTH = 2^1
  integer, parameter :: EAST = 2^2
  integer, parameter :: WEST = 2^3

contains

  function newHalo(...) result(this)
  end function newHalo

  subroutine update(this, distributedArray, directions)
    type (Halo_type) :: this
    real*8, intent(in) :: distributedArray(:)
    integer, intent(in) :: directions


  end subroutine update

end module SerialHalo_mod
