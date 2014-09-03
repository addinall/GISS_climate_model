#include "rundeck_opts.h"

      MODULE ICEDYN_COM
      integer :: imic=0,kticij=0
      real*8,dimension(:),allocatable :: rsix,rsiy,usi,vsi
      real*8,dimension(:),allocatable :: USIDT,VSIDT,RSISAVE
      real*8,dimension(:),allocatable :: icij
      contains
      subroutine alloc_icedyn_com(grid)
      use DOMAIN_DECOMP_ATM, only : DIST_GRID
      TYPE (DIST_GRID), INTENT(IN) :: grid
      return
      end subroutine alloc_icedyn_com
      END MODULE ICEDYN_COM

      MODULE ICEDYN
      use ICEDYN_COM, only : usi, vsi
      integer :: imicdyn=0
      contains
      subroutine alloc_icedyn(grid)
      use DOMAIN_DECOMP_ATM, only : DIST_GRID
      TYPE (DIST_GRID), INTENT(IN) :: grid
      return
      end subroutine alloc_icedyn
      END MODULE ICEDYN

      SUBROUTINE ICEDYN_DUM
!@sum ICEDYN_DUM dummy routines to replace ice dynamics
      !ENTRY alloc_icedyn
      ENTRY alloc_icedyn_com
      ENTRY gather_icdiags
      ENTRY io_icedyn
      ENTRY io_icdiag
      ENTRY reset_icdiag
      ENTRY ADVSI
      ENTRY init_icedyn
      ENTRY diag_ICEDYN
#ifdef NEW_IO
      entry def_rsf_icedyn
      entry new_io_icedyn
      entry def_rsf_icdiag
      entry new_io_icdiag
      entry set_ioptrs_iceacc_default
      entry def_meta_icdiag
      entry write_meta_icdiag
#endif
      RETURN
      END SUBROUTINE ICEDYN_DUM

      SUBROUTINE alloc_icedyn
      use ICEDYN_COM
      implicit none
      allocate( rsix(2), rsiy(2), usi(2), vsi(2))
      allocate( USIDT(2), VSIDT(2), RSISAVE(2) )
      allocate( icij(2) )
      rsix=0; rsiy=0; usi=0; vsi=0
      USIDT=0; VSIDT=0; RSISAVE=0
      icij=0
      END SUBROUTINE alloc_icedyn


      SUBROUTINE DYNSI
!@sum DYNSI simple coding to estimate ice-ocean friction velocity
!@auth Gavin Schmidt
      USE CONSTANT, only : rhows
      USE MODEL_COM, only : im,jm,kocean,focean,dtsrc
      USE DOMAIN_DECOMP_ATM, only : grid, get
      USE GEOM, only : imaxj
      USE SEAICE, only : oi_ustar0
      USE SEAICE_COM, only : rsi
      USE FLUXES, only : UI2rho,dmua,dmva

      IMPLICIT NONE
      INTEGER I,J
      REAL*8 ustar1
      INTEGER :: I_0,I_1, J_0,J_1


      IF (KOCEAN.eq.1) THEN

        CALL GET (grid, J_STRT=J_0,   J_STOP=J_1 )
        I_0 = GRID%I_STRT
        I_1 = GRID%I_STOP
        DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
c          UI2rho(I,J) = rhows*(oi_ustar0)**2  ! default
C**** with wind stress dependence
          if (rsi(i,j)*focean(i,j).gt.0) then
            ustar1= SQRT(SQRT(DMUA(I,J,2)**2+DMVA(I,J,2)**2)/(RSI(i,j)
     *           *focean(i,j)*DTSRC*RHOWS))
            UI2rho(I,J)=rhows*(oi_ustar0*max(1d0,1d3*ustar1))**2
          end if
        END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE DYNSI
