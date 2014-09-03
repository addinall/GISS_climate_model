! *** MPI Test Cases ***
! Do _NOT_ delete the line above

module test_TimerParallel_mod
   use pFUnit, only: assertEqual
   use pFUnit, only: TestInfo_type, mpiCommunicator, numProcesses, processRank
   use Timer_mod
   implicit none
   private

   public :: test_gather
   public :: nproc_test_gather

   integer, parameter :: nproc_test_gather(2) = (/ 2, 3 /)

contains

   subroutine test_gather(info)
      type (TestInfo_type), intent(in) :: info

      type (Timer_type) :: localTimer
      type (Timer_type) :: globalTimer
      character(len=200) :: report, expected

      integer :: rank, npes, ier
      integer :: i
      real (kind=r64) :: dt
      integer, parameter :: root = 0

      integer :: numTripsExpected
      real (kind=r64) :: exclusiveTimeExpected
      real (kind=r64) :: inclusiveTimeExpected
      real (kind=r64) :: minimumTimeExpected
      real (kind=r64) :: maximumTimeExpected
      real (kind=r64) :: time

      rank = processRank(info)
      npes = numProcesses(info)

      dt = 1 + rank
      time = 0

      do i = 0, rank
         call start(localTimer, time)
         time = time + dt
         call stop(localTimer, time)
      end do

      globalTimer = gather(localTimer, mpiCommunicator(info))

      numTripsExpected = sumN(1, npes)
      call assertEqual(numTripsExpected, getNumTrips(globalTimer),'numtrips')

      inclusiveTimeExpected = sumNsquared(1, npes) / npes
      exclusiveTimeExpected = sumNsquared(1, npes) / npes

      call assertEqual(inclusiveTimeExpected, getInclusiveTime(globalTimer),0.01_r64,'inclusive')
      call assertEqual(exclusiveTimeExpected, getExclusiveTime(globalTimer),0.01_r64,'exclusive')

      minimumTimeExpected = 1
      maximumTimeExpected = npes**2

      call assertEqual(minimumTimeExpected, getMinimumTime(globalTimer),0.01_r64,'min')
      call assertEqual(maximumTimeExpected, getMaximumTime(globalTimer),0.01_r64,'max')

      call assertEqual(0, getMinProcess(globalTimer), 'minProcess')
      call assertEqual(npes-1, getMaxProcess(globalTimer), 'maxProcess')

      call reset(globalTimer)
      call reset(localTimer)

   contains

      real function sumN(i0, i1) result (total)
         integer, intent(in) :: i0, i1
         total = (i1 - i0 + 1) * (i0 + i1) / 2
      end function sumN

      real function sumNsquared(i0, i1) result (total)
         integer, intent(in) :: i0, i1
         integer :: i
         total = sum ( (/ (i**2, i=i0,i1) /) )
      end function sumNsquared

   end subroutine test_gather
end module test_TimerParallel_mod
