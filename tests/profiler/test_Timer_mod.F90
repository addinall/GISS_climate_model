module test_Timer_mod
   use pFUnit, only: assertEqual, assertTrue, assertFalse, raiseException
   use Timer_mod
   implicit none
   private

   public :: setUp
   public :: tearDown
   public :: Fixture

   type Fixture
      private
      type (Timer_type) :: timer
      type (Timer_type) :: otherTimer
   end type Fixture

   public :: test_getNumTrips0
   public :: test_getNumTrips1
   public :: test_getNumTrips2
   public :: test_getNumTripsAfterReset

   public :: test_getTimeNoTrips

   public :: test_addTime
   public :: test_addTime2
   public :: test_getTimeAfterReset

   public :: test_reallyMeasuresTime

   public :: test_isActiveUnused
   public :: test_isActiveAfterStart
   public :: test_isActiveAfterStop
   public :: test_isActiveAfterReset

   public :: test_getMaximumTimeNoTrips
   public :: test_getMaximumTime
   public :: test_getMinimumTimeNoTrips
   public :: test_getMinimumTime
   public :: test_getAverageTripTimeNoTrips
   public :: test_getAverageTripTime
   public :: test_getExclusiveTime

   real(kind=r64), parameter :: EPS = EPSilon(1.)

contains

   subroutine setUp(this)
      type (Fixture) :: this
      call reset(this%timer)
      call reset(this%otherTimer)
   end subroutine setUp

   subroutine tearDown(this)
      type (Fixture) :: this
      call reset(this%timer)
      call reset(this%otherTimer)
   end subroutine tearDown

   subroutine test_getNumTrips0(this)
      type (Fixture) :: this
      call assertEqual(0, getNumTrips(this%timer))
   end subroutine test_getNumTrips0

   subroutine test_getNumTrips1(this)
      type (Fixture) :: this
      call addTrip(this%Timer)
      call assertEqual(1, getNumTrips(this%timer))
   end subroutine test_getNumTrips1

   subroutine test_getNumTrips2(this)
      type (Fixture) :: this
      call addTrip(this%Timer)
      call addTrip(this%Timer)
      call assertEqual(2, getNumTrips(this%timer))
   end subroutine test_getNumTrips2

   subroutine test_getNumTripsAfterReset(this)
      type (Fixture) :: this
      call start(this%timer)
      call stop(this%timer)
      call reset(this%timer)
      call assertEqual(0, getNumTrips(this%timer))
   end subroutine test_getNumTripsAfterReset

   subroutine test_getTimeNoTrips(this)
      type (Fixture) :: this
      call assertEqual(0._r64, getInclusiveTime(this%timer))
      call assertEqual(0._r64, getExclusiveTime(this%timer))
   end subroutine test_getTimeNoTrips

   subroutine test_addTime(this)
      type (Fixture) :: this
      real(kind=r64) :: dt1

      dt1 = 1.234_r64
      call start(this%timer, time = 0._r64)
      call stop(this%timer,  time = dt1)
      call assertEqual(dt1, getInclusiveTime(this%timer)) ! no tolerance
      call assertEqual(dt1, getExclusiveTime(this%timer)) ! no tolerance

   end subroutine test_addTime

   subroutine test_addTime2(this)
      type (Fixture) :: this
      real(kind=r64) :: dt1, dt2
      real(kind=r64) :: timeSecondTrip

      dt1 = 1.234_r64
      dt2 = 2.345_r64
      timeSecondTrip = 10.

      call start(this%timer, time = 0._r64)
      call stop(this%timer, time = dt1)
      call start(this%timer, time = timeSecondTrip)
      call stop(this%timer, time = timeSecondTrip+dt2)
      call assertEqual(dt1+dt2, getInclusiveTime(this%timer), EPS, 'inclusive')
      call assertEqual(dt1+dt2, getExclusiveTime(this%timer), EPS, 'exclusive')

   end subroutine test_addTime2

   subroutine test_getTimeAfterReset(this)
      type (Fixture) :: this
      real(kind=r64) :: dt1

      dt1 = 1.234_r64

      call start(this%timer, time = 0._r64)
      call stop(this%timer, time = dt1)
      call reset(this%timer)
      call assertEqual(0._r64, getInclusiveTime(this%timer))
      call assertEqual(0._r64, getExclusiveTime(this%timer))

   end subroutine test_getTimeAfterReset

   ! Use trusted system_clock(this) timer to ensure that a minimum time
   ! has elapsed between start/stop.
   subroutine test_reallyMeasuresTime(this)
      type (Fixture) :: this

      real (kind=r64), parameter :: SHORT_TIME = 0.1    ! seconds - a longish test
      real (kind=r64), parameter :: ACCURACY   = 0.0001 ! seconds
      real (kind=r64) :: time

      call start(this%timer)
      call pause(SHORT_TIME)
      call stop(this%timer)
      call assertTrue(getInclusiveTime(this%timer) >= SHORT_TIME)
      call assertEqual(getInclusiveTime(this%timer), getExclusiveTime(this%timer))

   contains

      subroutine pause(delay)
         real(kind=r64) :: delay
         integer :: count0, count1, countRate
         integer :: counter
         integer :: i
         integer, parameter :: MAX_COUNTER=HUGE(1)


         call system_clock(count0, countRate)
         counter = 0
         do
            counter = counter + 1
            if (counter == MAX_COUNTER) then ! possibly counter has cycled
               call raiseException('failed to complete timing loop')
               exit
            end if
            call system_clock(count1)
            time = real(count1 - count0)/countRate
            if (time > delay) exit
         end do
         return
      end subroutine pause

   end subroutine test_reallyMeasuresTime

   subroutine test_isActiveUnused(this)
      type (Fixture) :: this
      call assertFalse(isActive(this%timer))
   end subroutine test_isActiveUnused

   subroutine test_isActiveAfterStart(this)
      type (Fixture) :: this

      call start(this%timer)
      call assertTrue(isActive(this%timer))

   end subroutine test_isActiveAfterStart

   subroutine test_isActiveAfterStop(this)
      type (Fixture) :: this

      call start(this%timer)
      call stop(this%timer)
      call assertFalse(isActive(this%timer))

   end subroutine test_isActiveAfterStop

   subroutine test_isActiveAfterReset(this)
      type (Fixture) :: this

      call start(this%timer)
      call reset(this%timer)
      call assertFalse(isActive(this%timer))

   end subroutine test_isActiveAfterReset

   subroutine test_getMaximumTimeNoTrips(this)
      type (Fixture) :: this
      call assertEqual(0._r64, getMaximumTime(this%timer))
   end subroutine test_getMaximumTimeNoTrips

   subroutine test_getMaximumTime(this)
      type (Fixture) :: this
      real(kind=r64) :: dt1, dt2

      dt1 = 0.5
      dt2 = 1.5

      call start(this%timer, time = 0._r64)
      call stop(this%timer, time = dt1)
      call start(this%timer, time = 0._r64)
      call stop(this%timer, time = dt2)

      call assertEqual(max(dt1,dt2), getMaximumTime(this%timer))

   end subroutine test_getMaximumTime

   subroutine test_getMinimumTimeNoTrips(this)
      type (Fixture) :: this
      call assertEqual(0._r64, getMinimumTime(this%timer))
   end subroutine test_getMinimumTimeNoTrips

   subroutine test_getMinimumTime(this)
      type (Fixture) :: this
      real(kind=r64) :: dt1, dt2

      dt1 = 0.5
      dt2 = 1.5

      call start(this%timer, time = 0._r64)
      call stop(this%timer, time = dt1)
      call start(this%timer, time = 0._r64)
      call stop(this%timer, time = dt2)

      call assertEqual(min(dt1,dt2), getMinimumTime(this%timer))

   end subroutine test_getMinimumTime

   subroutine test_getAverageTripTimeNoTrips(this)
      type (Fixture) :: this
      call assertEqual(0., getAverageTripTime(this%timer), EPS)
   end subroutine test_getAverageTripTimeNoTrips

   subroutine test_getAverageTripTime(this)
      type (Fixture) :: this
      real(kind=r64) :: dt1, dt2

      dt1 = 0.5
      dt2 = 1.5

      call start(this%timer, time = 0._r64)
      call stop(this%timer, time = dt1)
      call start(this%timer, time = 0._r64)
      call stop(this%timer, time = dt2)

      call assertEqual((dt1+dt2)/2, getAverageTripTime(this%timer), EPS)

   end subroutine test_getAverageTripTime

   subroutine test_getExclusiveTime(this)
      type (Fixture) :: this
      real(kind=r64) :: dt1, dt2, dt3

      dt1 = 1.0
      dt2 = 3.0
      dt3 = 4.0

      ! add time to other timer _while_ 1st timer is still active
      call start(this%timer,       time = 0._r64)
      call start(this%otherTimer,  time = dt1)
      call stop(this%otherTimer,   time = dt2)
      call stop(this%timer,        time = dt3)

      call assertEqual(dt2-dt1, getExclusiveTime(this%otherTimer), EPS) 
      call assertEqual(dt3-(dt2-dt1), getExclusiveTime(this%timer), EPS) 

   end subroutine test_getExclusiveTime

end module test_Timer_mod
