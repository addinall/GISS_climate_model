module Tracers_mod
  use Dictionary_mod, only: Dictionary_type
  implicit none
  private

  public :: Tracer_type       ! derived type
  public :: TracerBundle_type ! derived type
  public :: Tracer            ! constructor
  public :: TracerBundle      ! constructor
  public :: addTracer
  public :: readFromText      
  public :: writeFormatted
  public :: readUnformatted
  public :: writeUnformatted
  public :: getTracer
  public :: getProperty
  public :: getProperties
  public :: setProperty
  public :: getIndex
  public :: getCount
  public :: makeSubset
  public :: hasProperty
  public :: checkMandatory
  public :: operator(==)
  public :: clean

  public :: NOT_FOUND

  type Tracer_type
    private
    type (Dictionary_type) :: properties
  end type Tracer_type

  type TracerBundle_type
    private
    type (Tracer_type), pointer :: tracers(:) => null()
  end type TracerBundle_type

  interface clean
    module procedure cleanTracer
    module procedure cleanBundle
  end interface

  integer, parameter :: NOT_FOUND = -1

  interface getProperty
    module procedure getProperty_single
    module procedure getProperty_multi
  end interface

  interface getProperties
    module procedure getProperties_single
    module procedure getProperties_multi
  end interface

  interface hasProperty
    module procedure hasProperty_single
    module procedure hasProperty_multi
  end interface

  interface setProperty
    module procedure setProperty_integer, setProperty_integerArr
    module procedure setProperty_real64,  setProperty_real64Arr
    module procedure setProperty_logical, setProperty_logicalArr
    module procedure setProperty_string,  setProperty_stringArr

    module procedure setProperty_multi_integer, setProperty_multi_integerArr
    module procedure setProperty_multi_real64,  setProperty_multi_real64Arr
    module procedure setProperty_multi_logical, setProperty_multi_logicalArr
    module procedure setProperty_multi_string,  setProperty_multi_stringArr
  end interface

  interface getCount
    module procedure getCount_
    module procedure getCount_property
  end interface

  interface operator(==)
    module procedure equals
  end interface

  interface writeFormatted
    module procedure writeFormatted_bundle
  end interface

  interface writeUnformatted
    module procedure writeUnformatted_bundle
    module procedure writeUnformatted_tracer
  end interface

  interface readUnformatted
    module procedure readUnformatted_bundle
    module procedure readUnformatted_tracer
  end interface

  integer, parameter :: LEN_HEADER = 80
  integer, parameter :: VERSION = 1
  character(len=*), parameter :: DESCRIPTION = 'TracerBundle'

contains

  function Tracer()
!@sum Construct empty tracer    
    use Dictionary_mod, only: Dictionary
    type (Tracer_type) :: Tracer
    Tracer%properties = Dictionary()
  end function Tracer

  function TracerBundle()
!@sum Construct empty tracer    
    use Dictionary_mod, only: Dictionary
    type (TracerBundle_type) :: TracerBundle
    allocate(TracerBundle%tracers(0))
  end function TracerBundle

  function readFromText(unit, defaultValues) result(bundle)
!@sum Populate a TracerBundle from a unit attached to a formatted
!@+ file.  Optionally apply default values.
    use Parser_mod
    use Dictionary_mod
    integer, intent(in) :: unit
    type (Dictionary_type), optional, intent(in) :: defaultValues
    type (TracerBundle_type) :: bundle
    type (Tracer_type) :: tracer
    type (Parser_type) :: parser

    integer :: status
  
    allocate(bundle%tracers(0))

    do
      tracer = readOneTracer(unit, status)
      if (status /= 0) exit
      if (present(defaultValues)) call merge(tracer%properties, defaultValues)
      call addTracer(bundle, tracer)
      call clean(tracer)
    end do

  end function readFromText

  subroutine addTracer(this, tracer)
!@sum Insert one tracer into a bundle.  Optionally
!@+ apply default values to the tracer dictionary.
    use Dictionary_mod
    type (TracerBundle_type), intent(inout) :: this
    type (Tracer_type), intent(in) :: tracer

    type (Tracer_type), allocatable :: oldList(:)
    integer :: count

    count = getNumTracers(this)
    allocate(oldList(count))
    oldList = this%tracers ! shallow copy ok
    deallocate(this%tracers)
    allocate(this%tracers(count+1))
    this%tracers(1:count) = oldList
    deallocate(oldList)

    this%tracers(count+1)%properties = Dictionary(tracer%properties)

  end subroutine addTracer

  integer function getNumTracers(this)
!@sum Returns number of tracers in a bundle.
    type (TracerBundle_type), intent(in) :: this
    getNumTracers = size(this%tracers)
  end function getNumTracers

  subroutine writeUnformatted_bundle(this, unit)
!@sum Write a bundle to a unit attached to an unformatted sequential file.
    type (TracerBundle_type), intent(in) :: this
    integer, intent(in) :: unit

    integer :: i, n

    character(len=LEN_HEADER) :: header

    write(header, '(a," - version:",1x,i5.0)') DESCRIPTION, VERSION
    write(unit) header

    n = getNumTracers(this)
    write(unit) n
    do i = 1, n
      call writeUnformatted(this%tracers(i), unit)
    end do

  end subroutine writeUnformatted_bundle

  subroutine writeUnformatted_tracer(this, unit)
!@sum Write a tracer to a unit attached to an unformatted sequential file.
    use Dictionary_mod, only: writeUnformatted
    type (Tracer_type), intent(in) :: this
    integer, intent(in) :: unit

    call writeUnformatted(this%properties, unit)
  end subroutine writeUnformatted_tracer

  subroutine readUnformatted_bundle(this, unit)
!@sum Read a bundle to a unit attached to an unformatted sequential file.
    type (TracerBundle_type), intent(out) :: this
    integer, intent(in) :: unit

    integer :: i, n
    integer :: oldVersion
    character(len=len(DESCRIPTION)) :: tag
    character(len=LEN_HEADER) :: header

    read(unit) header
    read(header, '(a," - version:",i5.0)') tag, oldVersion

    if (tag /= DESCRIPTION) then
      call throwException(DESCRIPTION // '::readUnformatted() - incorrect header.', 14)
    end if

    if (oldVersion /= VERSION) then
      call throwException(DESCRIPTION // '::readUnformatted() - unsupported format.', 14)
    end if
      
    read(unit) n
    allocate(this%tracers(n))
    do i = 1, n
      call readUnformatted(this%tracers(i), unit)
    end do

  end subroutine readUnformatted_bundle

  subroutine readUnformatted_tracer(this, unit)
!@sum Read a bundle to a unit attached to an unformatted sequential file.
    use Dictionary_mod, only: readUnformatted
    type (Tracer_type), intent(out) :: this
    integer, intent(in) :: unit
    call readUnformatted(this%properties, unit)
  end subroutine readUnformatted_tracer

  subroutine writeFormatted_bundle(this, unit)
    use Parser_mod
    type (TracerBundle_type), intent(in) :: this
    integer, intent(in) :: unit

    type (Parser_type) :: parser
    type (Dictionary_type), pointer :: properties
    integer :: i

    call setBeginData(parser, '{')
    call setEndData(parser, '}')
    call setTokenSeparators(parser, '=,')

    do i = 1, getCount(this)
      properties => getProperties(this%tracers(i))
      call writeFormatted(parser, unit, properties)
    end do
    
  end subroutine writeFormatted_bundle

  function readOneTracer(unit, status) result(tracer)
    use Parser_mod, only: Parser_type
    use Parser_mod, only: setBeginData, setEndData
    use Parser_mod, only: setCommentCharacters, setTokenSeparators
    use Parser_mod, only: parse

    integer, intent(in) :: unit
    integer, intent(out) :: status
    type (Tracer_type) :: tracer

    type (Parser_type) :: parser

    call setBeginData(parser, '{')
    call setEndData(parser, '}')
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')

    tracer%properties = parse(parser, unit, status)
    if (status /= 0) return

  end function readOneTracer

  integer function getIndex(this, name) result(index)
    use Dictionary_mod, only: lookup
    use KeyValuePair_mod, only: MAX_LEN_KEY
    use GenericType_mod
    type (TracerBundle_type), intent(in) :: this
    character(len=*), intent(in) :: name

    integer :: i
    character(len=MAX_LEN_KEY) :: aName

    index = NOT_FOUND
    do i = 1, getCount(this)
      if (all(name == lookup(this%tracers(i)%properties, 'name'))) then
        index = i
        exit
      end if
    end do

  end function getIndex

  logical function hasProperty_single(this, property) result(has)
    use Dictionary_mod, only: hasKey
    type (Tracer_type), intent(in) :: this
    character(len=*), intent(in) :: property
    has = hasKey(this%properties, property)
  end function hasProperty_single

  function hasProperty_multi(this, property) result(has)
    use Dictionary_mod, only: hasKey
    type (TracerBundle_type), intent(in) :: this
    character(len=*), intent(in) :: property
    logical :: has(size(this%tracers))

    integer :: i

    has = (/ (hasProperty(this%tracers(i), property), i = 1, getCount(this)) /)
  end function hasProperty_multi

  subroutine checkMandatory(this, properties)
    type (TracerBundle_type), intent(in) :: this
    character(len=*), intent(in) :: properties(:)

    integer :: i, j
    do i = 1, getCount(this)
      do j = 1, size(properties)
        if (hasProperty(this%tracers(i), properties(j))) cycle
        call throwException("Tracers_mod::checkMandatory() - " // &
             &  trim(getName(this%tracers(i))) // " missing property '" // &
             &  trim(properties(j)) // "'.", 14)
        return
      end do
    end do
  end subroutine checkMandatory

  function getName(this) result(name)
    use GenericType_mod
    use KeyValuePair_mod, only: MAX_LEN_KEY
    type (Tracer_type), intent(in) :: this
    character(len=MAX_LEN_KEY) :: name
    name = getProperty(this, 'name')
  end function getName

  function getTracer(this, species) result(tracer)
    type (TracerBundle_type), target, intent(in) :: this
    character(len=*), intent(in) :: species
    type (Tracer_type), pointer :: tracer


    integer :: index

    index = getIndex(this, species)
    if (index /= NOT_FOUND) then
      tracer => this%tracers(index)
    else
      tracer => null()
    end if
    
  end function getTracer

  function getProperties_single(this) result(properties)
    type (Tracer_type), target, intent(in) :: this
    type (Dictionary_type), pointer :: properties
    properties => this%properties
  end function getProperties_single

  function getProperties_multi(this, species) result(properties)
    type (TracerBundle_type), target, intent(in) :: this
    character(len=*), intent(in) :: species
    type (Dictionary_type), pointer :: properties

    integer :: index

    index = getIndex(this, species)
    if (index /= NOT_FOUND) then
      properties => this%tracers(index)%properties
    else
      nullify(properties)
    end if

  end function getProperties_multi

  function getProperty_single(this, property) result (propertyValues)
    use GenericType_mod
    use Dictionary_mod, only: lookup
    type (Tracer_type), intent(in) :: this
    character(len=*), intent(in) :: property
    type (GenericType_type), pointer :: propertyValues(:)

    propertyValues => lookup(this%properties, property)
    
  end function getProperty_single

  function getProperty_multi(this, species, property) result (propertyValues)
    use GenericType_mod
    use Dictionary_mod, only: lookup
    type (TracerBundle_type), intent(in) :: this
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: property
    type (GenericType_type), pointer :: propertyValues(:)

    integer :: index

    index = getIndex(this, species)
    if (index /= NOT_FOUND) then
      propertyValues => lookup(this%tracers(index)%properties, property)
    else
      allocate(propertyValues(0))
    end if
    
  end function getProperty_multi

  subroutine setProperty_integer(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    integer, intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_integer

  subroutine setProperty_integerArr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    integer, intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_integerArr

  subroutine setProperty_real64(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    real*8, intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_real64

  subroutine setProperty_real64Arr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    real*8, intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_real64Arr

  subroutine setProperty_logical(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    logical, intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_logical

  subroutine setProperty_logicalArr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    logical, intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_logicalArr

  subroutine setProperty_string(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    character(len=*), intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_string

  subroutine setProperty_stringArr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    character(len=*), intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_stringArr


  ! Select a tracer from list
  subroutine setProperty_multi_integer(this, name, property, value)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    integer, intent(in) :: value
    
    integer :: index

    index = getIndex(this, name)
    call insert(this%tracers(index)%properties, property, value)
  end subroutine setProperty_multi_integer

  subroutine setProperty_multi_integerArr(this, name, property, values)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    integer, intent(in) :: values(:)
    
    integer :: index

    index = getIndex(this, name)
    call insert(this%tracers(index)%properties, property, values)
  end subroutine setProperty_multi_integerArr

  subroutine setProperty_multi_real64(this, name, property, value)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    real*8, intent(in) :: value
    
    integer :: index

    index = getIndex(this, name)
    call insert(this%tracers(index)%properties, property, value)
  end subroutine setProperty_multi_real64

  subroutine setProperty_multi_real64Arr(this, name, property, values)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    real*8, intent(in) :: values(:)
    
    integer :: index

    index = getIndex(this, name)
    call insert(this%tracers(index)%properties, property, values)
  end subroutine setProperty_multi_real64Arr

  subroutine setProperty_multi_logical(this, name, property, value)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    logical, intent(in) :: value
    
    integer :: index

    index = getIndex(this, name)
    call insert(this%tracers(index)%properties, property, value)
  end subroutine setProperty_multi_logical

  subroutine setProperty_multi_logicalArr(this, name, property, values)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    logical, intent(in) :: values(:)
    
    integer :: index

    index = getIndex(this, name)
    call insert(this%tracers(index)%properties, property, values)
  end subroutine setProperty_multi_logicalArr

  subroutine setProperty_multi_string(this, name, property, value)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    character(len=*), intent(in) :: value
    
    integer :: index

    index = getIndex(this, name)
    call insert(this%tracers(index)%properties, property, value)
  end subroutine setProperty_multi_string

  subroutine setProperty_multi_stringArr(this, name, property, values)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    character(len=*), intent(in) :: values(:)
    
    integer :: index

    index = getIndex(this, name)
    call insert(this%tracers(index)%properties, property, values)
  end subroutine setProperty_multi_stringArr

  integer function getCount_(this) result(number)
    type (TracerBundle_type), intent(in) :: this

    number = size(this%tracers)
    
  end function getCount_

  integer function getCount_property(this, withProperty) result(number)
    type (TracerBundle_type), intent(in) :: this
    character(len=*), intent(in) :: withProperty

    number = count(hasProperty(this, withProperty))
    
  end function getCount_property

  subroutine makeSubset(this, withProperty, subset)
    type (TracerBundle_type), target, intent(in) :: this
    character(len=*), intent(in) :: withProperty
    type (TracerBundle_type) :: subset

    integer :: i, j

    allocate(subset%tracers(getCount(this, withProperty)))

    j = 1
    do i = 1, getCount(this)
      if (hasProperty(this%tracers(i), withProperty)) then
        subset%tracers(j)%properties = this%tracers(i)%properties
        j = j + 1
      end if
    end do

  end subroutine makeSubset

  logical function equals(bundleA, bundleB) result(isEqual)
    use GenericType_mod
    use Dictionary_mod, only: operator(==)
    use Parser_mod, only: MAX_LEN_TOKEN
    type (TracerBundle_type), target, intent(in) :: bundleA
    type (TracerBundle_type), target, intent(in) :: bundleB

    integer :: i
    character(len=MAX_LEN_TOKEN) :: name

    isEqual = .true.
    do i = 1, getCount(bundleA)
      name = getProperty(bundleA%tracers(1), 'name')
      if (.not. (getProperties(bundleA, name) == getProperties(bundleB, name))) then
        isEqual = .false.
        exit
      end if
    end do

  end function equals

  subroutine cleanTracer(this)
    use Dictionary_mod, only: clean
    type (Tracer_type), intent(inout) :: this
    call clean(this%properties)
  end subroutine cleanTracer

  subroutine cleanBundle(this)
    use Dictionary_mod, only: clean
    type (TracerBundle_type), intent(inout) :: this
    integer :: i
    do i = 1, getCount(this)
      call clean(this%tracers(i)%properties)
    end do
    deallocate(this%tracers)
  end subroutine cleanBundle

end module Tracers_mod
