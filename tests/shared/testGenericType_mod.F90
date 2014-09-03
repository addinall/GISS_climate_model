module testGenericType_mod
!@sum Tests GenericType implementation.
!@auth T. Clune
  use pFUnit
  use GenericType_mod
  implicit none

contains

  subroutine testGetInteger()
!@sum Test that GenericType can both be assigned to an integer,
!@+ and also can be compared against an integer.
    integer :: found

    found = GenericType(1)
    call assertEqual(1, found)

    call assertTrue(1 == GenericType(1))
    call assertTrue(2 == GenericType(2))
    call assertFalse(1 == GenericType(2))
  end subroutine testGetInteger

  subroutine testGetIntegerArray()
!@sum Test that GenericType can both be assigned to an array of integers,
!@+ and also can be compared against an array of integers.
    integer, parameter :: NUM_VALUES=2
    integer :: expected(NUM_VALUES), found(NUM_VALUES)

    expected = [1,2]
    found = GenericType(expected)
    call assertEqual(expected, found)
    call assertTrue(all(expected == GenericType(expected)))

  end subroutine testGetIntegerArray

  subroutine testGetReal64()
!@sum Test that GenericType can both be assigned to a double,
!@+ and also can be compared against a double.
    real(kind=r64) :: expected, found

    expected = 1.234d+0
    found = GenericType(expected)
    call assertEqual(expected, found)
    call assertTrue(expected == GenericType(expected))
    call assertFalse(expected-1 == GenericType(expected))

  end subroutine testGetReal64

  subroutine testToStringInteger()
!@sum Test that GenericType can be converted to a string.
    call assertEqual('1', toString(GenericType(1)))
    call assertEqual('2', toString(GenericType(2)))
  end subroutine testToStringInteger

  subroutine testToStringReal64()
!@sum Test that GenericType can be converted to a string.
!@+ Need to ensure that value -> string -> value reproduces
!@+ exactly.
    character(len=MAX_LEN_VALUE) :: string
    real(kind=r64) :: expectedValue, foundValue

    expectedValue = 321.1234567890123456789d+0
    string = toString(GenericType(expectedValue))
    read(string,'(g)') foundValue
    call assertEqual(expectedValue, foundValue)

  end subroutine testToStringReal64

  subroutine testToStringLogical()
!@sum Test that GenericType can be converted to a string.
    call assertEqual('.true.', toString(GenericType(.true.)))
    call assertEqual('.false.', toString(GenericType(.false.)))
  end subroutine testToStringLogical

  subroutine testFromStringInteger()
!@sum Test that GenericType can be constructed from a string.
    call assertTrue(1 == GenericType('1', INTEGER_TYPE))
    call assertTrue(1234 == GenericType(' 1234', INTEGER_TYPE))
  end subroutine testFromStringInteger
    
  subroutine testFromStringReal64()
!@sum Test that GenericType can be constructed from a string.
    real(kind=r64) :: expected
    type (GenericType_type) :: generic

    expected = 1.23d+0
    generic = GenericType(expected)
    call assertTrue(expected == GenericType(toString(generic), REAL64_TYPE))

    expected = -123.789d+0
    generic = GenericType(expected)
    call assertTrue(expected == GenericType(toString(generic), REAL64_TYPE))

  end subroutine testFromStringReal64

  subroutine testFromStringLogical()
!@sum Test that GenericType can be constructed from a string.
    call assertTrue(.false. == GenericType('f', LOGICAL_TYPE))
    call assertTrue(.false. == GenericType('false', LOGICAL_TYPE))
    call assertTrue(.false. == GenericType('.False.', LOGICAL_TYPE))

    call assertTrue(.true. == GenericType('T', LOGICAL_TYPE))
    call assertTrue(.true. == GenericType('true', LOGICAL_TYPE))
    call assertTrue(.true. == GenericType('.True.', LOGICAL_TYPE))
  end subroutine testFromStringLogical
    
  subroutine testFailFromStringInteger()
!@sum Test that GenericType() returs appropriate errors when args are
!@+ inconsistent.
    integer :: i

    i = GenericType('1.234', INTEGER_TYPE)
    call assertFailedAssert('GenericType::GenericType() - cannot convert string "1.234" to integer.', &
         & 'Failed to detect conversion error.')

    i = GenericType('1.234', -1)
    call assertFailedAssert('GenericType::GenericType() - no such type.')

  end subroutine testFailFromStringInteger
    
  subroutine testFailFromStringLogical()
!@sum Test that GenericType() returs appropriate errors when args are
!@+ inconsistent.
    logical :: flag

    flag = GenericType('fa', LOGICAL_TYPE)
    call assertFailedAssert('GenericType::GenericType() - cannot convert string "fa" to logical.', &
         & 'Failed to detect conversion error.')

  end subroutine testFailFromStringLogical
    
end module testGenericType_mod
