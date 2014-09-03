module testKeyValuePair_mod
  use pFUnit
  use GenericType_mod
  use KeyValuePair_mod
  implicit none
  private

  public :: testGetKey
  public :: testGetNumValues
  public :: testGetValue
  public :: testEquals
  public :: testReadWriteUnformatted

contains

  subroutine testGetKey()
    type (KeyValuePair_type) :: pair

    pair = KeyValuePair('aKey', GenericType(1))
    call assertEqual('aKey', getKey(pair)) 

  end subroutine testGetKey

  subroutine testGetNumValues()
    type (KeyValuePair_type) :: pair

    pair = KeyValuePair('aKey', GenericType(1))
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType(1.234d+0))
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType(.true.))
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType('hello'))
    call assertEqual(1, getNumValues(pair))

    ! lists
    pair = KeyValuePair('aKey', GenericType([1,2,3]))
    call assertEqual(3, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType([1.234d+0]))
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType([.true.,.false.,.true.]))
    call assertEqual(3, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType(['hello', 'goodbye']))
    call assertEqual(2, getNumValues(pair))

  end subroutine testGetNumValues

  subroutine testGetValue()
    integer :: expected(2)

    call assertTrue(1 == getValue(KeyValuePair('aKey', GenericType(1))))

    expected = [1,2]
    call assertTrue(all(expected == getValues(KeyValuePair('aKey', GenericType(expected)))))

  end subroutine testGetValue

  subroutine testEquals()
    integer, parameter :: NUM_CASES = 6
    type (KeyValuePair_type) :: pairs(6)
    integer :: i, j
    logical :: flag

    pairs(1) = KeyValuePair('a', GenericType(1)) 
    pairs(2) = KeyValuePair('b', GenericType(1))  ! different key
    pairs(3) = KeyValuePair('b', GenericType(2))  ! different value
    pairs(4) = KeyValuePair('a', GenericType([1,2])) ! different length
    pairs(5) = KeyValuePair('a', GenericType([2,1])) ! different order
    pairs(6) = KeyValuePair('a', GenericType([.true.,.false.])) ! different type

    do i = 1, NUM_CASES
      do j = 1, NUM_CASES
        flag = (pairs(i) == pairs(j))
        if (i == j) then
          call assertTrue(flag)
        else
          call assertFalse(flag)
        end if
      end do
    end do

  end subroutine testEquals

  subroutine testReadWriteUnformatted()
    use GenericType_mod
    use FileManager
    integer :: unit

    type (KeyValuePair_type) :: expectedA
    type (KeyValuePair_type) :: expectedB
    type (KeyValuePair_type) :: expectedC

    type (KeyValuePair_type) :: foundA
    type (KeyValuePair_type) :: foundB
    type (KeyValuePair_type) :: foundC

    call openUnit('testKeyValuePair.txt', unit, qold=.false., qbin=.true.)

    expectedA = KeyValuePair('one', GenericType([1,2,3]))
    expectedB = KeyValuePair('two', GenericType(.true.))
    expectedC = KeyValuePair('three', GenericType([1.23d+0, 2.34d+0]))

    call writeUnformatted(expectedA, unit)
    call writeUnformatted(expectedB, unit)
    call writeUnformatted(expectedC, unit)

    rewind(unit)

    call readUnformatted(foundA, unit)
    call readUnformatted(foundB, unit)
    call readUnformatted(foundC, unit)

    call assertTrue(expectedA == foundA, 'case A')
    call assertTrue(expectedB == foundB, 'case B')
    call assertTrue(expectedC == foundC, 'case C')
    close(unit, status='delete')

  end subroutine testReadWriteUnformatted

end module testKeyValuePair_mod
