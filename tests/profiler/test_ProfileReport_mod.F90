module test_ProfileReport_mod
   use pFUnit
   use ProfileReport_mod
   use ReportColumn_mod
   use TimerList_mod
   implicit none
   private

   public :: test_oneColumn
   public :: test_nameTrip
   public :: test_percent
   public :: test_minPerDay
   public :: test_generateReport

contains

   subroutine test_oneColumn()
      type (ProfileReport_type) :: report
      type (TimerList_type) :: timers

      report = newReport()
      call addColumn(report, newColumn(NAME_COLUMN, fieldWidth=10))

      call reset(timers)
      call addTimer(timers,'A')

      call assertEqual('   Name', getHeader(report))
      call assertEqual('', getUnits(report))
      call assertEqual('A',    getLine(report, timers, 1))

      call reset(timers)
   end subroutine test_oneColumn

   subroutine test_nameTrip()
      use Timer_mod
      type (ProfileReport_type) :: report
      type (TimerList_type) :: timers

      report = newReport()
      call addColumn(report, newColumn(NAME_COLUMN, fieldWidth=10))
      call addColumn(report, newColumn(TRIP_COUNTS_COLUMN, fieldWidth=3))

      call reset(timers)
      call addTimer(timers,'A')
      call start(timers, 'A')
      call stop(timers, 'A')
      call start(timers, 'A')
      call stop(timers, 'A')
      !                       1234567890123456
      call assertEqual('   Name    Tri ', getHeader(report))
      call assertEqual('               ', getUnits(report))
      call assertEqual('A            2 ', getLine(report, timers, 1))

      call reset(timers)
   end subroutine test_nameTrip

   subroutine test_percent()
      use Timer_mod
      type (ProfileReport_type) :: report
      type (TimerList_type) :: timers
      type (ReportColumn_type) :: column
      real(kind=r64) :: totalTime

      report = newReport()
      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth=7)
      call reset(timers)
      call addTimer(timers,'A')
      call addTimer(timers,'B')
      call start(timers, 'A', 0._r64)
      call stop(timers, 'A',  1._r64)
      call start(timers, 'B', 1._r64)
      call stop(timers, 'B',  4._r64)

      totalTime = 4._r64
      call setScale(column, 100/totalTime, 'Percent')
      call addColumn(report, column)


      !                       1234567890123
      call assertEqual('Inclusi', getHeader(report))
      call assertEqual('Percent', getUnits(report))
      call assertEqual('  25.00', getLine(report, timers, 1))
      call assertEqual('  75.00', getLine(report, timers, 2))

      call reset(timers)
   end subroutine test_percent

   subroutine test_minPerDay()
      use Timer_mod
      type (ProfileReport_type) :: report
      type (TimerList_type) :: timers
      type (ReportColumn_type) :: column
      real(kind=r64) :: simulatedTime
      integer, parameter :: MINUTES_PER_DAY = 60*24
      integer, parameter :: SECONDS_PER_DAY = 60 * MINUTES_PER_DAY

      report = newReport()
      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth=7)
      call reset(timers)
      call addTimer(timers,'A')
      call addTimer(timers,'B')
      call start(timers, 'A', 0._r64)
      call stop(timers, 'A',  60._r64)
      call start(timers, 'B', 60._r64)
      call stop(timers, 'B', 180._r64)

      simulatedTime = SECONDS_PER_DAY ! day
      call setScale(column, MINUTES_PER_DAY/simulatedTime, 'min/day')
      call addColumn(report, column)

      !                       1234567890123
      call assertEqual('Inclusi', getHeader(report))
      call assertEqual('min/day', getUnits(report))
      call assertEqual('   1.00', getLine(report, timers, 1))
      call assertEqual('   2.00', getLine(report, timers, 2))

      call reset(timers)
   end subroutine test_minPerDay

   subroutine test_generateReport()
      use Timer_mod
      type (ProfileReport_type) :: report
      type (TimerList_type) :: timers
      type (ReportColumn_type) :: column
      real(kind=r64) :: totalTime
      character(len=MAX_RECORD_LENGTH), pointer :: lines(:)

      report = newReport()
      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth=8)
      call reset(timers)
      call addTimer(timers,'A')
      call addTimer(timers,'B')
      call start(timers, 'A', 0._r64)
      call stop(timers, 'A',  1._r64)
      call start(timers, 'B', 1._r64)
      call stop(timers, 'B',  4._r64)

      totalTime = 4._r64
      call setScale(column, 100/totalTime, ' %')
      call addColumn(report, newColumn(NAME_COLUMN, 10))
      call addColumn(report, column)
      call addColumn(report, newColumn(TRIP_COUNTS_COLUMN, 4))

      lines => generateReport(report, timers)

      !                       12345678901234567890123456
      call assertEqual('-------------------------', trim(lines(1)))
      call assertEqual('   Name    Inclusiv Trip ', trim(lines(2)))
      call assertEqual('               %         ', trim(lines(3)))
      call assertEqual('-------------------------', trim(lines(4)))
      call assertEqual('A             25.00    1 ', trim(lines(5)))
      call assertEqual('B             75.00    1 ', trim(lines(6)))
      call assertEqual('-------------------------', trim(lines(7)))
      deallocate(lines)
      call reset(timers)
   end subroutine test_generateReport

end module test_ProfileReport_mod
