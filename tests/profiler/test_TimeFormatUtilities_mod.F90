module test_TimeFormatUtilities_mod
   use pFUnit
   use TimeFormatUtilities_mod

   implicit none
   private

   public :: test_formatHMS
   public :: test_formatSeconds

contains

   subroutine test_formatHMS()
      real (kind=r64) :: time
      character(len=80) :: expected

      time = 0.
      expected = '  0:00:00.00'
      call assertEqual(expected, format('HMS',time))

      time = 0.12
      expected = '  0:00:00.12'
      call assertEqual(expected, format('HMS',time))

      time = time + 3
      expected = '  0:00:03.12'
      call assertEqual(expected, format('HMS',time))

      time = time + 4 * 60
      expected = '  0:04:03.12'
      call assertEqual(expected, format('HMS',time))

      time = time + 123 * 3600
      expected = '123:04:03.12'
      call assertEqual(expected, format('HMS',time))

      ! rounding should produce correct results
      time = 9.999999
      expected = '  0:00:10.00'
      call assertEqual(expected, format('HMS',time))

      time = 3599.999999
      expected = '  1:00:00.00'
      call assertEqual(expected, format('HMS',time))

   end subroutine test_formatHMS

   subroutine test_formatSeconds()
      real (kind=r64) :: time
      character(len=80) :: expected

      time = 0.12
      expected = '    0.12'
      call assertEqual(expected, format('seconds', time))

      time = time + 23
      expected = '   23.12'
      call assertEqual(expected, format('seconds', time))

      time = time + 34 * 60
      expected = ' 2063.12'
      call assertEqual(expected, format('seconds', time))

   end subroutine test_formatSeconds

end module test_TimeFormatUtilities_mod
