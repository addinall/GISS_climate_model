#include "rundeck_opts.h"

module DOMAIN_DECOMP_1D
  use dist_grid_mod
  use Halo_mod
  use SpecialIO_mod
  use GatherScatter_mod
  use GlobalSum_mod
  implicit none
  public

end module DOMAIN_DECOMP_1D

#ifdef MPI_DEFS_HACK
#include "mpi_defs.h"
#endif

#ifndef CUBED_SPHEREa
#define DOMAIN_DECOMP_ATM_IS_1D
#endif

#ifdef NEW_IO
#define USE_DD2D_UTILS
#endif

#ifdef DOMAIN_DECOMP_ATM_IS_1D
! If the atmosphere has a 1D domain decomposition, pass along the contents
! of DOMAIN_DECOMP_1D
      MODULE DOMAIN_DECOMP_ATM
        use dist_grid_mod
        use Halo_mod
        use SpecialIO_mod
        use GatherScatter_mod
        use GlobalSum_mod
      implicit none
      public
      END MODULE DOMAIN_DECOMP_ATM
#endif
