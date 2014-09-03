module test_TimerList_mod
   use TimerList_mod
   use pFUnit
   implicit none
   private

   public :: fixture, setUp, tearDown

   public :: test_initialize
   public :: test_addTimer
   public :: test_getTimer
   public :: test_startTimerNotFound
   public :: test_stopTimerNotFound
   public :: test_detectUnbalanced

   type fixture
      private
      integer :: placeholder
   end type fixture

contains

   subroutine setUp(this)
      type (fixture) :: this
      call initialize(mockTime = 0._r64)
   end subroutine setUp

   subroutine tearDown(this)
      use Timer_mod
      use TimeFormatUtilities_mod
      type (fixture) :: this
      type (Timer_type), pointer :: pTimer
      call finalize(mockTime = 10*3600._r64)
      pTimer => getTimer('main')

      call reset()
   end subroutine tearDown

   subroutine test_initialize(this)
      type (fixture) :: this
      call assertEqual(1,getNumTimers())
      call assertTrue(isTimer('main'))
   end subroutine test_initialize

   subroutine test_addTimer(this)
      type (fixture) :: this
      call assertFalse(isTimer('timer A'),'timer A should not be found')
      call addTimer('timer A')
      call assertTrue(isTimer('timer A'),'timer A not found')
      call assertEqual(2, getNumTimers())

      call assertFalse(isTimer('timer B'),'timer B should not be found')
      call addTimer('timer B')
      call assertTrue(isTimer('timer A'),'timer A not found')
      call assertTrue(isTimer('timer B'),'timer B not found')
      call assertEqual(3, getNumTimers())

   end subroutine test_addTimer

   subroutine test_getTimer(this)
      use Timer_mod
      type (fixture) :: this
      type (Timer_type), pointer :: timer_A

      call addTimer('timer A')
      call addTimer('timer B')

      timer_A => getTimer('timer A')

      call assertEqual(0, getNumTrips(timer_A))
      call start('timer A')
      call stop('timer A')
      call assertEqual(1, getNumTrips(timer_A))
      call start('timer A')
      call stop('timer A')
      call assertEqual(2, getNumTrips(timer_A))

   end subroutine test_getTimer

   ! For convenience create a timer if it does not exist
   ! Could abort instead, but that would be annoying.
   subroutine test_startTimerNotFound(this)
      type (fixture) :: this
      character(len=*), parameter :: name = 'undeclared Timer'
      call start(name)
      call assertTrue(isTimer(name))
      call stop(name)
   end subroutine test_startTimerNotFound

   ! Cannot stop a timer that does not exist.  Even creating it
   ! makes no sense at this point.
   subroutine test_stopTimerNotFound(this)
      type (fixture) :: this
      character(len=*), parameter :: name = 'undeclared Timer'
#ifdef USE_PFUNIT
      call stop(name)
      call assertTrue(catch('Timer <' // trim(name) // '> has not been declared prior to use.'))
#endif
   end subroutine test_stopTimerNotFound

   subroutine test_detectUnbalanced(this)
      type (fixture) :: this
#ifdef USE_PFUNIT
      call addTimer('timer')
      call start('timer')
      call finalize()
      call assertTrue(catch('Unbalanced start/stop for timer <timer>.'))
      call reset('timer')
#endif
   end subroutine test_detectUnbalanced

end module test_TimerList_mod
