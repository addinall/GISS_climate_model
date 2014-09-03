      module growthallometry
!@sum Routines to calculate uptake and allocation of nutrients within an
!@sum individual plant, growth, and allometry.

      use ent_types
      use ent_const

      implicit none


      contains
      !*********************************************************************
      
      subroutine uptake_N(dtsec, pp)
      real*8 :: dtsec           !dt in seconds
      type(patch),pointer :: pp

      !*****************************************************************
      !*              FILL IN CODE                                     *
      !*****************************************************************

      
      end subroutine uptake_N
      !*********************************************************************

      subroutine init_rootdistr(fracroot, pft)
      !use ent_GISSveg, only : GISS_calc_fracroot
      real*8 :: fracroot(N_DEPTH)
      integer :: pft

      !GISS TEMPORARY - prescribed roots should be passed in
!      call GISS_calc_fracroot(fracroot, pft)

      ! at least set it to zero, since it is called in zero_cohort...
      fracroot(:) = 0.d0

      !* Prognostic roots calculated here *!.

      end subroutine init_rootdistr
      !*********************************************************************
      end module growthallometry
