module testStringUtilities_mod
!@sum Tests procedures in StringUtilities_mod.
!@auth T. Clune
  use pFUnit
  use StringUtilities_mod
  implicit none

  public :: testToLowerCase

contains

  subroutine testToLowerCase()
!@sum Test that letters are switched to lower case, but that other symbols are unchanged.
    call assertEqual('a', toLowerCase('a'), 'trivial case')
    call assertEqual('a', toLowerCase('A'), 'simple case')
    call assertEqual('abcd', toLowerCase('AbCd'),'multiple letters')
    call assertEqual('abcdef123abcdef{!}%#', toLowerCase('AbCdEf123aBcDeF{!}%#'), &
         & 'other symbols')
  end subroutine testToLowerCase

end module testStringUtilities_mod

