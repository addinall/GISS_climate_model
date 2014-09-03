module testTracers_mod
  use pFUnit
  use Tracers_mod
  use Dictionary_mod, only: Dictionary_type
  implicit none
  private

  public :: testReadFromText
  public :: testGetIndex
  public :: testHasPropertySingle
  public :: testHasPropertyMulti
  public :: testSetProperty
  public :: testSetPropertyMulti
  public :: testGetPropertySingle
  public :: testGetPropertyMulti
  public :: testCountHasProperty
  public :: testGetSubset
  public :: testCheckMandatory
  public :: testDefaultValues
  public :: testWriteFormatted
  public :: testWriteUnformatted

  public :: setUp
  public :: tearDown
  public :: fixture


  type fixture
    type (TracerBundle_type) :: bundle
    type (Dictionary_type) :: defaultValues
    integer :: unit
    integer :: numTracers
  end type fixture

contains

  subroutine setUp(this)
    use Dictionary_mod, only: Dictionary
    use Dictionary_mod, only: insert
    use FileManager
    type (fixture) :: this
    integer :: unit

    call openUnit('testTracers.txt', unit, qold=.false., qbin=.false.)
    this%unit = unit

    write(unit,*)'{'
    write(unit,*)' name = speciesA'
    write(unit,*)' molecularMass = 1.2 ! g'
    write(unit,*)' optionA = True'
    write(unit,*)' optionC = False'
    write(unit,*)'}'
    write(unit,*)' '
    write(unit,*)'{'
    write(unit,*)' name = speciesB'
    write(unit,*)' molecularMass = 2.3 ! g'
    write(unit,*)' optionB = False'
    write(unit,*)'}'

    rewind(unit)

    this%numTracers = 2

    this%defaultValues = Dictionary()
    call insert(this%defaultValues, 'optionC', .true.)
    call insert(this%defaultValues, 'optionD',  1)

  end subroutine setUp

  subroutine tearDown(this)
    use Dictionary_mod, only: clean
    type (fixture) :: this

    close(this%unit, status='delete')
    this%unit = -1
    call clean(this%defaultValues)
    call clean(this%bundle)
  end subroutine tearDown

  subroutine testReadFromText(this)
    type (fixture) :: this

    integer :: unit

    unit = this%unit
    this%bundle = readFromText(unit, this%defaultValues)

    call assertEqual(this%numTracers, getCount(this%bundle), &
         & 'Incorrect number of tracers found.')

  end subroutine testReadFromText

  ! support routine
  subroutine readTracers(this)
    type (fixture) :: this
    this%bundle = readFromText(this%unit, this%defaultValues)
  end subroutine readTracers

  subroutine testGetIndex(this)
    type (fixture) :: this

    call readTracers(this)

    call assertEqual(1, getIndex(this%bundle, 'speciesA'))
    call assertEqual(2, getIndex(this%bundle, 'speciesB'))
    call assertEqual(NOT_FOUND, getIndex(this%bundle, 'speciesZ'))

  end subroutine testGetIndex

  subroutine testGetPropertySingle(this)
    use GenericType_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found
    type (Tracer_type), pointer :: tracer

    call readTracers(this)

    expected = 1.2d+0
    tracer => getTracer(this%bundle, 'speciesA')
    found = getProperty(tracer, 'molecularMass')
    call assertEqual(expected, found)

  end subroutine testGetPropertySingle

  subroutine testGetPropertyMulti(this)
    use GenericType_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found

    call readTracers(this)

    expected = 1.2d+0
    found = getProperty(this%bundle, 'speciesA', 'molecularMass')
    call assertEqual(expected, found)

  end subroutine testGetPropertyMulti

  subroutine testSetProperty(this)
    use GenericType_mod
    use Dictionary_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found
    type (Dictionary_type), pointer :: properties
    type (Tracer_type), pointer :: tracer

    call readTracers(this)
    tracer => getTracer(this%bundle, 'speciesA')

    call setProperty(tracer, 'otherProperty', expected)
    found = getProperty(tracer, 'otherProperty')
    call assertEqual(expected, found)

  end subroutine testSetProperty

  subroutine testSetPropertyMulti(this)
    use GenericType_mod
    use Dictionary_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found
    type (Dictionary_type), pointer :: properties
    type (Tracer_type), pointer :: tracer

    call readTracers(this)

    expected = 2.1d+0
    tracer => getTracer(this%bundle, 'speciesB')
    call setProperty(tracer, 'otherProperty', expected)

    found = getProperty(this%bundle, 'speciesB', 'otherProperty')
    call assertEqual(expected, found)

  end subroutine testSetPropertyMulti

  subroutine testHasPropertySingle(this)
    use GenericType_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found
    type (Tracer_type), pointer :: tracer

    call readTracers(this)
    tracer => getTracer(this%bundle,'speciesA')

    call assertTrue(hasProperty(tracer, 'optionA'))
    call assertFalse(hasProperty(tracer, 'optionB'))

  end subroutine testHasPropertySingle

  subroutine testHasPropertyMulti(this)
    use GenericType_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found

    call readTracers(this)

    call assertTrue(all(hasProperty(this%bundle, 'molecularMass')))
    call assertTrue(all([.true.,.false.] .eqv. hasProperty(this%bundle, 'optionA')))

  end subroutine testHasPropertyMulti

  subroutine testCountHasProperty(this)
    use GenericType_mod
    type (fixture) :: this

    type (Tracer_type), pointer :: subset(:)

    call readTracers(this)
    call assertEqual(2, getCount(this%bundle, 'name'))
    call assertEqual(1, getCount(this%bundle, 'optionA'))

  end subroutine testCountHasProperty

  subroutine testGetSubset(this)
    use GenericType_mod
    use KeyValuePair_mod
    type (fixture) :: this

    type (TracerBundle_type) :: subset
    character(len=MAX_LEN_KEY) :: name

    call readTracers(this)

    call makeSubset(this%bundle, withProperty='optionA', subset=subset)
    call assertEqual(1, getCount(subset))
    call assertTrue(associated(getTracer(subset, 'speciesA')))

  end subroutine testGetSubset

  subroutine testCheckMandatory(this)
    use KeyValuePair_mod
    type (fixture) :: this
    character(len=MAX_LEN_KEY) :: mandatoryProperties(2)

    call readTracers(this)

    mandatoryProperties(1) = 'name'
    mandatoryProperties(2) = 'molecularMass'

    call checkMandatory(this%bundle, mandatoryProperties)
    if (catch(preserve=.true.)) return

    mandatoryProperties(1) = 'optionA'
    call checkMandatory(this%bundle, mandatoryProperties)
    call assertFailedAssert( &
         & "Tracers_mod::checkMandatory() - speciesB missing property 'optionA'.", &
         & "Failed to detect missing mandatory option.")

  end subroutine testCheckMandatory

  subroutine testDefaultValues(this)
    use GenericType_mod
    type (fixture) :: this
    logical :: expectedOptionC, foundOptionC
    integer :: expectedOptionD, foundOptionD
    
    call readTracers(this)

    foundOptionC = getProperty(this%bundle, 'speciesA', 'optionC')
    expectedOptionC = .false.
    call assertTrue(expectedOptionC .eqv. foundOptionC)

    foundOptionC = getProperty(this%bundle, 'speciesB', 'optionC')
    expectedOptionC = .true.
    call assertTrue(expectedOptionC .eqv. foundOptionC)

    foundOptionD = getProperty(this%bundle, 'speciesA', 'optionD')
    expectedOptionD = 1
    call assertTrue(expectedOptionD == foundOptionD)
    
  end subroutine testDefaultValues

  subroutine testWriteFormatted(this)
    use FileManager
    type (fixture) :: this
    type (TracerBundle_type) :: bundle
    integer :: unit

    call readTracers(this)

    call openUnit('testTracersOut.txt', unit, qold=.false., qbin=.false.)
    call writeFormatted(this%bundle, unit)
    rewind(unit)

    bundle = readFromText(unit)
    call assertTrue(bundle == this%bundle)

    close(unit, status ='delete')
    call clean(bundle)
    
  end subroutine testWriteFormatted

  subroutine testWriteUnformatted(this)
    use FileManager
    type (fixture) :: this
    type (TracerBundle_type) :: bundle
    integer :: unit

    call readTracers(this)

    call openUnit('testTracersOut.bin', unit, qold=.false., qbin=.true.)
    call writeUnformatted(this%bundle, unit)
    rewind(unit)

    call readUnformatted(bundle, unit)
    call assertTrue(bundle == this%bundle)

    close(unit, status ='delete')
    call clean(bundle)
    
  end subroutine testWriteUnformatted

end module testTracers_mod
