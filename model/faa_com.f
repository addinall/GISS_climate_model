#include "rundeck_opts.h"

      MODULE faa_com
!@sum variables for high resolution aviation emissions
      use resolution, only : lm 
      integer, parameter :: nqty=6
      integer, parameter :: nhours=24

C**************  Latitude-Dependant (allocatable) *******************
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: emiss_5d



c      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
c     &                  grid%j_strt_halo:grid%j_stop_halo,lm,24,nqty) ::
c     &     emiss_5d  


          END MODULE faa_com

      subroutine alloc_faa_com(grid)
      use domain_decomp_atm, only : dist_grid, get
      use model_com, only     : im,lm
      use faa_com, only: emiss_5d

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H
      logical :: init = .false.
      integer, parameter :: nqty=6
      integer, parameter :: nhours=24


      if(init)return
      init=.true.
    
      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO

      allocate(    emiss_5d(I_0H:I_1H,J_0H:J_1H,LM,nhours,nqty)   )

      return
      end subroutine alloc_faa_com
