 module testDictionary_mod
  use pfunit
  use Dictionary_mod
  use KeyValuePair_mod
  use GenericType_mod
  implicit none
  private

  ! tests of legacy functionality
  public :: testSetInteger
  public :: testSetReal
  public :: testSetString

  public :: testSetIntegerList
  public :: testSetRealList
  public :: testSetStringList

  public :: testQuery
  public :: testReadWrite
  public :: testReadWriteUnformatted

  ! tests of new implementation
  public :: testKeyNotFound
  public :: testGetValueA
  public :: testGetValueB
  public :: testGetValueReal64

  public :: testGetNumEntries
  public :: testHasKey
  public :: testGetKeys
  
  public :: testFailDuplicateKey

  public :: testEqualsA
  public :: testEqualsB
  public :: testEqualsC
  public :: testEqualsD
  public :: testEqualsE
  public :: testEqualsF
  public :: testEqualsG
  public :: testEqualsH

  public :: testMergeInteger
  public :: testMergeDictionary

contains

  subroutine testSetInteger()
    integer :: expected
    integer :: found

    call reset()
    expected = 1
    call set_param('testParam_Int', expected)
    call get_param('testParam_Int', found)
    call assertEqual(expected, found)

  end subroutine testSetInteger

  subroutine testSetReal()
    real(kind=r64) :: expected
    real(kind=r64) :: found

    call reset()
    expected = 1.234
    call set_param('testParam_Real', expected)
    call get_param('testParam_Real', found)
    call assertEqual(expected, found)
    
  end subroutine testSetReal

  subroutine testSetString()
    character(len=20) :: expected
    character(len=20) :: found

    call reset()
    expected = 'testValue'
    call set_param('testParam_String', expected)
    call get_param('testParam_String', found)
    call assertEqual(expected, found)
    
  end subroutine testSetString

  subroutine testSetIntegerList()
    integer, parameter :: NUM_ELEMENTS = 2
    integer :: expected(NUM_ELEMENTS)
    integer :: found(NUM_ELEMENTS)

    call reset()
    expected = [1,2]
    call set_param('testParam_IntList', expected, NUM_ELEMENTS)
    call get_param('testParam_IntList', found, NUM_ELEMENTS)
    call assertEqual(expected, found)

  end subroutine testSetIntegerList

  subroutine testSetRealList()
    integer, parameter :: NUM_ELEMENTS = 2
    real(kind=r64) :: expected(NUM_ELEMENTS)
    real(kind=r64) :: found(NUM_ELEMENTS)

    expected = [1.1, 2.2]
    call set_param('testParam_RealList', expected, NUM_ELEMENTS)
    call get_param('testParam_RealList', found, NUM_ELEMENTS)
    call assertEqual(expected, found)

  end subroutine testSetRealList

  subroutine testSetStringList()
    integer, parameter :: NUM_ELEMENTS = 2
    character*(128) :: expected(NUM_ELEMENTS)
    character*(128) :: found(NUM_ELEMENTS)
    integer :: i

    call reset()
    expected = ['hello ','goodbye']
    call set_param('testParam_StringList', expected, NUM_ELEMENTS)
    call get_param('testParam_StringList', found, NUM_ELEMENTS)
    do i = 1, NUM_ELEMENTS
      call assertEqual(expected(i), found(i))
    end  do

  end subroutine testSetStringList

  subroutine testQuery()
    integer, parameter :: NUM_ELEMENTS = 2
    integer :: values(NUM_ELEMENTS)

    character(len=20) :: expectedName, foundName
    character(len=1) :: ptype
    integer :: dim

    call reset()
    values = [1,2]
    expectedName = 'testquery'

    call set_param(expectedName, values, NUM_ELEMENTS)
    call query_param(1, foundName, dim, ptype)

    call assertEqual(expectedName, foundName, 'check order that tests are run')
    call assertEqual(NUM_ELEMENTS, dim)
    call assertEqual('i', ptype)

  end subroutine testQuery

  subroutine testReadWrite()
    use FileManager
    integer :: unit
    character(len=20) :: foundName
    integer :: dim
    character(len=1) :: pType

    call reset()
    call set_param('testReadWriteA', 2)
    call set_param('testReadWriteB', [3,4], 2)

    call openUnit('dictionary.dat',unit,qold=.false.,qbin=.true.)
    call write_param(unit)
    close(unit)

    call reset()
    ! start with one parameter just to make sure things were cleared
    ! before reading
    call set_param('testReadWriteC', 'aString')

    call openUnit('dictionary.dat',unit,qold=.true.,qbin=.true.)
    call read_param(unit, ovrwrt = .false.)
    call query_param(3, foundName, dim, pType)
    close(unit)

    call assertEqual('testreadwriteb', foundname)
    call assertEqual(2, dim)
    call assertEqual('i', pType)
    call deleteFile('dictionary.dat')

  end subroutine testReadWrite

  subroutine testReadWriteUnformatted()
    use FileManager
    integer :: unit
    type (Dictionary_type) :: dictionaryA
    type (Dictionary_type) :: dictionaryB

    dictionaryA = Dictionary()
    call insert(dictionaryA, 'key1', 1)
    call insert(dictionaryA, 'key2', [.true., .false.])
    call insert(dictionaryA, 'key3', [1.2d+0,2.d+0])
    call insert(dictionaryA, 'key4', 'string')
    
    call openUnit('dictionary.dat',unit,qold=.false.,qbin=.true.)
    call writeUnformatted(dictionaryA, unit)

    rewind(unit)
    call readUnformatted(dictionaryB, unit)

    call assertTrue(dictionaryA == dictionaryB)
    close(unit, status='delete')
    
  end subroutine testReadWriteUnformatted

  subroutine testKeyNotFound()
    type (Dictionary_type) :: aDictionary
    integer :: i

    aDictionary = Dictionary()
    i = lookup(aDictionary, 'alpha')
    call assertFailedAssert('GenericType_mod: nonconforming shapes.')
    call assertFailedAssert('Key not found: <alpha>.', &
         & 'Failed to detect missing key.')

  end subroutine testKeyNotFound

  subroutine testGetValueA()
    type (Dictionary_type) :: aDictionary
    integer :: expected, found

    expected = 1
    aDictionary = Dictionary()
    call insert(aDictionary, 'key', expected)
    found = lookup(aDictionary, 'key')
    call assertEqual(expected, found)

    call clean(aDictionary)

  end subroutine testGetValueA

  subroutine testGetValueB()
    type (Dictionary_type) :: aDictionary
    integer :: expectedA, foundA
    integer :: expectedB, foundB

    expectedA = 1
    expectedB = 2
    aDictionary = Dictionary()
    call insert(aDictionary, 'keyA', expectedA)
    call insert(aDictionary, 'keyB', expectedB)

    foundB = lookup(aDictionary, 'keyB')
    call assertEqual(expectedB, foundB)

    foundA = lookup(aDictionary, 'keyA')
    call assertEqual(expectedA, foundA)

    call clean(aDictionary)

  end subroutine testGetValueB

  subroutine testGetValueReal64()
    type (Dictionary_type) :: aDictionary
    real(kind=r64) :: expected, found

    expected = 1
    aDictionary = Dictionary()
    call insert(aDictionary, 'key', expected)
    found = lookup(aDictionary, 'key')
    call assertEqual(expected, found)

    call clean(aDictionary)

  end subroutine testGetValueReal64

  subroutine testGetNumEntries()
    type (Dictionary_type) :: aDictionary
    real(kind=r64) :: expected, found

    expected = 1
    aDictionary = Dictionary()

    call assertEqual(0, getNumEntries(aDictionary))

    call insert(aDictionary, 'key1', expected)
    call assertEqual(1, getNumEntries(aDictionary))

    call insert(aDictionary, 'key2', expected)
    call assertEqual(2, getNumEntries(aDictionary))

  end subroutine testGetNumEntries

  subroutine testHasKey()
    type (Dictionary_type) :: aDictionary

    aDictionary = Dictionary()
    call assertFalse(hasKey(aDictionary, 'a'))

    call insert(aDictionary, 'a', 1)
    call assertTrue(hasKey(aDictionary, 'a'))
    call assertFalse(hasKey(aDictionary, 'b'))

    call insert(aDictionary, 'b', 1)
    call assertTrue(hasKey(aDictionary, 'a'))
    call assertTrue(hasKey(aDictionary, 'b'))

    call clean(aDictionary)

  end subroutine testHasKey

  subroutine testGetKeys()
    type (Dictionary_type) :: aDictionary
    character(len=MAX_LEN_KEY), pointer :: keys(:)

    aDictionary = Dictionary()
    call insert(aDictionary, 'a', 1)
    call insert(aDictionary, 'b', 1)

    keys => getKeys(aDictionary)
    call assertEqual(2, size(keys))
    call assertTrue(any(keys == 'b'))
    call assertTrue(any(keys == 'a'))

    call clean(aDictionary)
  end subroutine testGetKeys

  subroutine testFailDuplicateKey()
    type (Dictionary_type) :: aDictionary
    aDictionary = Dictionary()
    call insert(aDictionary, 'a', 1)
    call insert(aDictionary, 'a', 2)
    call assertFailedAssert('Dictionary: duplicate key - <a>.', &
         & 'Failed to diagnose duplicate key.')
    call clean(aDictionary)
  end subroutine testFailDuplicateKey

  subroutine testEqualsA()
    type (Dictionary_type) :: a
    type (Dictionary_type) :: b
    
    a = Dictionary()
    b = Dictionary()
    call assertTrue(a == b)

  end subroutine testEqualsA

  subroutine testEqualsB()
    type (Dictionary_type) :: a
    type (Dictionary_type) :: b
    
    a = Dictionary()
    b = Dictionary()
    call insert(a, 'key1', 1)
    call assertFalse(a == b)

  end subroutine testEqualsB

  subroutine testEqualsC()
    type (Dictionary_type) :: a
    type (Dictionary_type) :: b
    
    a = Dictionary()
    b = Dictionary()
    call insert(a, 'key1', 1)
    call insert(b, 'key1', 1)
    call assertTrue(a == b)

  end subroutine testEqualsC

  subroutine testEqualsD()
    type (Dictionary_type) :: a
    type (Dictionary_type) :: b
    
    a = Dictionary()
    b = Dictionary()
    call insert(a, 'key1', 1)
    call insert(b, 'key2', 1)
    call assertFalse(a == b)

  end subroutine testEqualsD

  subroutine testEqualsE()
    type (Dictionary_type) :: a
    type (Dictionary_type) :: b
    
    a = Dictionary()
    b = Dictionary()
    call insert(a, 'key1', 1)
    call insert(b, 'key1', 2)
    call assertFalse(a == b)

  end subroutine testEqualsE

  subroutine testEqualsF()
    type (Dictionary_type) :: a
    type (Dictionary_type) :: b
    
    a = Dictionary()
    b = Dictionary()
    call insert(a, 'key1', [1,2])
    call insert(b, 'key1', [2])
    call assertFalse(a == b)

  end subroutine testEqualsF

  subroutine testEqualsG()
    type (Dictionary_type) :: a
    type (Dictionary_type) :: b
    
    a = Dictionary()
    b = Dictionary()
    call insert(a, 'key1', [1,2])
    call insert(b, 'key1', [1,2])
    call assertTrue(a == b)

  end subroutine testEqualsG

  subroutine testEqualsH()
    type (Dictionary_type) :: a
    type (Dictionary_type) :: b
    
    a = Dictionary()
    b = Dictionary()
    call insert(a, 'key1', [1,2])
    call insert(a, 'key2', .true.)
    call insert(a, 'key3', 'hello')

    call insert(b, 'key3', 'hello')
    call insert(b, 'key2', .true.)
    call insert(b, 'key1', [1,2])

    call assertTrue(a == b)

    call clean(a)
    call clean(b)

  end subroutine testEqualsH

  subroutine testMergeInteger()
    type (Dictionary_type) :: aDictionary
    integer :: probe
    integer :: two

    aDictionary = Dictionary()
    call insert(aDictionary, 'a', 1)
    
    ! 
    probe = 2
    call merge(aDictionary, 'a', probe)
    call assertEqual(1, probe)
    call assertTrue(all(1 == lookup(aDictionary, 'a')))

    two = 2
    call merge(aDictionary, 'b', two)
    call assertEqual(2, two)
    call assertTrue(all(2 == lookup(aDictionary, 'b')))

    call clean(aDictionary)
  end subroutine testMergeInteger

  subroutine testMergeDictionary()
    use GenericType_mod
    type (Dictionary_type) :: aDictionary
    type (Dictionary_type) :: otherDictionary

    integer :: expectedA
    real(kind=r64) :: expectedB
    logical :: expectedC

    aDictionary = Dictionary()
    otherDictionary = Dictionary()

    call insert(aDictionary, 'a', 1)
    call insert(aDictionary, 'b', 2.34d+0)

    call insert(otherDictionary, 'a', 1) ! same
    call insert(otherDictionary, 'b', 1.23d+0) ! different
    call insert(otherDictionary, 'c', .true.) ! new

    call merge(aDictionary, otherDictionary)

    expectedA = 1
    expectedB = 2.34d+0
    expectedC = .true.

    call assertTrue(all(expectedA == lookup(aDictionary, 'a')))
    call assertTrue(all(expectedB == lookup(aDictionary, 'b')))
    call assertTrue(all(expectedC == lookup(aDictionary, 'c')))

    call clean(otherDictionary)
    call clean(aDictionary)
  end subroutine testMergeDictionary

end module testDictionary_mod
