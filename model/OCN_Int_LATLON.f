#include "rundeck_opts.h"

      MODULE INT_AG2OG_MOD

!@sum INT_AG2OG_MOD contains subroutines for conversion 2D, 3D, etc. 
!!    arrays from atm. to the ocean grid 
!@auth Larissa Nazarenko
!@ver  1.0
      use hntrp_mod
      IMPLICIT NONE
      SAVE
      PRIVATE
      PUBLIC INT_AG2OG

      Interface INT_AG2OG
c _par are parallelized versions
      Module Procedure INT_AG2OG_2Da_par
      Module Procedure INT_AG2OG_2Db ! need _par for tracers
      Module Procedure INT_AG2OG_3Da ! need _par for tracers
      Module Procedure INT_AG2OG_3Db_par
      Module Procedure INT_AG2OG_3Dc_par
      Module Procedure INT_AG2OG_4Da ! need _par for tracers
      Module Procedure INT_AG2OG_4Db ! need _par for tracers
      Module Procedure INT_AG2OG_Vector1_par
      Module Procedure INT_AG2OG_Vector2
      End Interface

      type(hntrp_type) :: hntrp_a2o ! instance for atm A -> ocn A
      logical :: hntrp_a2o_needs_init = .true.

      contains

      SUBROUTINE INT_AG2OG_2Da(aA,oA,aWEIGHT_loc)

!@sum INT_AG2OG_2D is for conversion 2D arrays from atm. to the ocean grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_DATA, UNPACK_DATA
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER

      REAL*8, INTENT(IN)  :: 
     *        aWEIGHT_loc(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(IN)  :: aA(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(OUT) :: oA(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, ALLOCATABLE :: aA_glob(:,:), oA_glob(:,:)
     *                     , aFtemp(:,:),oFtemp(:,:), aWEIGHT(:,:)

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oA(1:aIM,1:aJM) = aA(1:aIM,1:aJM)
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aA_glob(aIM,aJM), STAT = IER)
        ALLOCATE(oA_glob(oIM,oJM), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(aWEIGHT(aIM,aJM), STAT = IER)
      end if

C***  Gather 2D array on atmospheric grid into the global array

      CALL PACK_DATA (aGRID, aA, aA_glob)

      CALL PACK_DATA (aGRID, aWEIGHT_loc, aWEIGHT)

C***  Interpolate aA_glob from atmospheric grid to ocean grid 

      if(AM_I_ROOT()) then

        call HNTR80 (aIM,aJM,0.d0,aDLATM, oIM,oJM,0.d0,oDLATM, 0.d0)

        aFtemp(:,:) = aA_glob(:,:)
        aFtemp(2:aIM,aJM) = aFtemp(1,aJM)
        call HNTR8P (aWEIGHT, aFtemp, oFtemp)
        oA_glob(:,:) = oFtemp(:,:)
      end if

C***  Scatter global array oA_glob to the ocean grid

      CALL UNPACK_DATA (ogrid, oA_glob, oA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob, aFtemp,oFtemp, aWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_AG2OG_2Da

      SUBROUTINE INT_AG2OG_2Db(aA1,aA2,oA,aWEIGHT_loc) 

!@sum INT_AG2OG_2D is for conversion 2D arrays from atm. to the ocean grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      USE GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_DATA, UNPACK_DATA
      USE OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER

      REAL*8, INTENT(IN)  :: 
     *        aWEIGHT_loc(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(IN) :: aA1(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(IN) :: aA2(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(OUT) :: oA(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, ALLOCATABLE :: aA1_glob(:,:),aA2_glob(:,:), oA_glob(:,:)
     *                     , aFtemp(:,:),oFtemp(:,:), aWEIGHT(:,:)

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oA(1:aIM,1:aJM) = aA1(1:aIM,1:aJM)*aA2(1:aIM,1:aJM)
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aA1_glob(aIM,aJM), STAT = IER)
        ALLOCATE(aA2_glob(aIM,aJM), STAT = IER)
        ALLOCATE(oA_glob(oIM,oJM), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(aWEIGHT(aIM,aJM), STAT = IER)
      end if

C***  Gather 2D array on atmospheric grid into the global array

      CALL PACK_DATA (aGRID, aA1, aA1_glob)
      CALL PACK_DATA (aGRID, aA2, aA2_glob)

      CALL PACK_DATA (aGRID, aWEIGHT_loc, aWEIGHT)

C***  Interpolate aA_glob from atmospheric grid to ocean grid 

      if(AM_I_ROOT()) then

        call HNTR80 (aIM,aJM,0.d0,aDLATM, oIM,oJM,0.d0,oDLATM, 0.d0)

        aFtemp(:,:) = 0.
        aFtemp(:,:) = aA1_glob(:,:)*aA2_glob(:,:)
        aFtemp(2:aIM,aJM) = aFtemp(1,aJM)
        call HNTR8P (aWEIGHT, aFtemp, oFtemp)
        oA_glob(:,:) = oFtemp(:,:)
      end if
C***  Scatter global array oA_glob to the ocean grid

      CALL UNPACK_DATA (ogrid, oA_glob, oA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA1_glob,aA2_glob, oA_glob, aFtemp,oFtemp, aWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_AG2OG_2Db

      SUBROUTINE INT_AG2OG_3Da(aA,oA,aWEIGHT_loc,NT)

!@sum INT_AG2OG_3D is for conversion 3D arrays from atm. to the ocean grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT, PACK_DATA
     *                           , PACK_COLUMN, UNPACK_COLUMN
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER, NT,N

      REAL*8, INTENT(IN)  :: 
     *        aWEIGHT_loc(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(IN)  :: 
     *  aA(NT,aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(OUT) :: 
     *  oA(NT,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, ALLOCATABLE :: aA_glob(:,:,:), oA_glob(:,:,:)
     *                     , aFtemp(:,:),oFtemp(:,:), aWEIGHT(:,:)

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oA(1:aIM,1:aJM,:) = aA(1:aIM,1:aJM,:)
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aA_glob(NT,aIM,aJM), STAT = IER)
        ALLOCATE(oA_glob(NT,oIM,oJM), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(aWEIGHT(aIM,aJM), STAT = IER)
      end if

C***  Gather 3D array on atmospheric grid into the global array

      CALL PACK_COLUMN (aGRID, aA, aA_glob)

      CALL PACK_DATA (aGRID, aWEIGHT_loc, aWEIGHT)

C***  Interpolate aA_glob from atmospheric grid to ocean grid 

      if(AM_I_ROOT()) then

        call HNTR80 (aIM,aJM,0.d0,aDLATM, oIM,oJM,0.d0,oDLATM, 0.d0)

        DO N=1,NT
          aFtemp(:,:) = aA_glob(N,:,:)
          aFtemp(2:aIM,aJM) = aFtemp(1,aJM)
          call HNTR8P (aWEIGHT, aFtemp, oFtemp)
          oA_glob(N,:,:) = oFtemp(:,:)
        END DO
      end if

C***  Scatter global array oA_glob to the ocean grid

      CALL UNPACK_COLUMN (ogrid, oA_glob, oA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob, aFtemp,oFtemp, aWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_AG2OG_3Da

      SUBROUTINE INT_AG2OG_3Db(aA,oA,aWEIGHT_loc, aN,oN)

!@sum INT_AG2OG_3D is for conversion 3D arrays from atm. to the ocean grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_DATA, UNPACK_DATA
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER, aN,oN

      REAL*8, INTENT(IN)  :: 
     *        aWEIGHT_loc(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(IN)  :: 
     *  aA(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN) 
      REAL*8, INTENT(OUT) :: 
     *  oA(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN)

      REAL*8, ALLOCATABLE :: aA_glob(:,:,:), oA_glob(:,:,:)
     *                     , aFtemp(:,:),oFtemp(:,:), aWEIGHT(:,:)

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oA(1:aIM,1:aJM,oN) = aA(1:aIM,1:aJM,oN)
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aA_glob(aIM,aJM,aN), STAT = IER)
        ALLOCATE(oA_glob(oIM,oJM,oN), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(aWEIGHT(aIM,aJM), STAT = IER)
      end if

C***  Gather 3D array on atmospheric grid into the global array

      CALL PACK_DATA (aGRID, aA, aA_glob)

      CALL PACK_DATA (aGRID, aWEIGHT_loc, aWEIGHT)

C***  Interpolate aA_glob from atmospheric grid to ocean grid 

      if(AM_I_ROOT()) then

        call HNTR80 (aIM,aJM,0.d0,aDLATM, oIM,oJM,0.d0,oDLATM, 0.d0)

        aFtemp(:,:) = 0.
        aFtemp(:,:) = aA_glob(:,:,oN)
        aFtemp(2:aIM,aJM) = aFtemp(1,aJM)
        call HNTR8P (aWEIGHT, aFtemp, oFtemp)
        oA_glob(:,:,oN) = oFtemp(:,:)
      end if

C***  Scatter global array oA_glob to the ocean grid

      CALL UNPACK_DATA (ogrid, oA_glob, oA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob, aFtemp,oFtemp, aWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_AG2OG_3Db

      SUBROUTINE INT_AG2OG_3Dc(aA,oA,aWEIGHT_loc, aN,oN,oNin)

!@sum INT_AG2OG_3D is for conversion 3D arrays from atm. to the ocean grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT, PACK_DATA
     *                           , PACK_COLUMN, UNPACK_COLUMN
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER, aN,oN, oNin 

      REAL*8, INTENT(IN)  :: 
     *        aWEIGHT_loc(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(IN)  :: 
     *  aA(aN,aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(OUT) :: 
     *  oA(oN,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, ALLOCATABLE :: aA_glob(:,:,:), oA_glob(:,:,:)
     *                     , aFtemp(:,:),oFtemp(:,:), aWEIGHT(:,:)

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oA(oNin,1:aIM,1:aJM) = aA(oNin,1:aIM,1:aJM)
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aA_glob(aN,aIM,aJM), STAT = IER)
        ALLOCATE(oA_glob(oN,oIM,oJM), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(aWEIGHT(aIM,aJM), STAT = IER)
      end if

C***  Gather 3D array on atmospheric grid into the global array

      CALL PACK_COLUMN (aGRID, aA, aA_glob)
      CALL PACK_COLUMN (oGRID, oA, oA_glob)

      CALL PACK_DATA (aGRID, aWEIGHT_loc, aWEIGHT)

C***  Interpolate aA_glob from atmospheric grid to ocean grid 

      if(AM_I_ROOT()) then

        call HNTR80 (aIM,aJM,0.d0,aDLATM, oIM,oJM,0.d0,oDLATM, 0.d0)

        aFtemp(:,:) = 0.
        aFtemp(:,:) = aA_glob(oNin,:,:)
        aFtemp(2:aIM,aJM) = aFtemp(1,aJM)
        call HNTR8P (aWEIGHT, aFtemp, oFtemp)
        oA_glob(oNin,:,:) = oFtemp(:,:)
      end if

C***  Scatter global array oA_glob to the ocean grid

      CALL UNPACK_COLUMN (ogrid, oA_glob, oA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob, aFtemp,oFtemp, aWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_AG2OG_3Dc

      SUBROUTINE INT_AG2OG_4Da(aA,oA,aWEIGHT_loc, NT, aN,oN)

!@sum INT_AG2OG_4D is for conversion 4D arrays from atm. to the ocean grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT, PACK_DATA
     *                           , PACK_BLOCK, UNPACK_BLOCK
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER, NT,N, aN,oN 

      REAL*8, INTENT(IN)  :: 
     *        aWEIGHT_loc(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(IN)  :: 
     *  aA(NT,aN,aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(OUT) :: 
     *  oA(NT,oN,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, ALLOCATABLE :: aA_glob(:,:,:,:), oA_glob(:,:,:,:)
     *                     , aFtemp(:,:),oFtemp(:,:), aWEIGHT(:,:)

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oA(:,oN,1:aIM,1:aJM) = aA(:,oN,1:aIM,1:aJM)
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aA_glob(NT,aN,aIM,aJM), STAT = IER)
        ALLOCATE(oA_glob(NT,oN,oIM,oJM), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(aWEIGHT(aIM,aJM), STAT = IER)
      end if

C***  Gather 4D array on atmospheric grid into the global array

      CALL PACK_BLOCK (aGRID, aA, aA_glob)

      CALL PACK_DATA (aGRID, aWEIGHT_loc, aWEIGHT)

C***  Interpolate aA_glob from atmospheric grid to ocean grid 

      if(AM_I_ROOT()) then

        call HNTR80 (aIM,aJM,0.d0,aDLATM, oIM,oJM,0.d0,oDLATM, 0.d0)

        DO N=1,NT
          aFtemp(:,:) = 0.
          aFtemp(:,:) = aA_glob(N,oN,:,:)
          aFtemp(2:aIM,aJM) = aFtemp(1,aJM)
          call HNTR8P (aWEIGHT, aFtemp, oFtemp)
          oA_glob(N,oN,:,:) = oFtemp(:,:)
        END DO
      end if

C***  Scatter global array oA_glob to the ocean grid

      CALL UNPACK_BLOCK (ogrid, oA_glob, oA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob, aFtemp,oFtemp, aWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_AG2OG_4Da

      SUBROUTINE INT_AG2OG_4Db(aA,oA,oB,aWEIGHT_loc, NT, aN,oN)

!@sum INT_AG2OG_4D is for conversion 4D arrays from atm. to the ocean grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT, PACK_DATA
     *                           , PACK_BLOCK, UNPACK_BLOCK
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER, NT,N, aN,oN 

      REAL*8, INTENT(IN)  :: 
     *        aWEIGHT_loc(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(IN)  :: 
     *  aA(NT,aN,aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(OUT) :: 
     *  oA(NT,oN,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, DIMENSION(oIM,oJM,NT) :: oB

      REAL*8, ALLOCATABLE :: aA_glob(:,:,:,:), oA_glob(:,:,:,:)
     *                     , aFtemp(:,:),oFtemp(:,:), aWEIGHT(:,:)

      if(AM_I_ROOT()) then
        ALLOCATE(aA_glob(NT,aN,aIM,aJM), STAT = IER)
        ALLOCATE(oA_glob(NT,oN,oIM,oJM), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(aWEIGHT(aIM,aJM), STAT = IER)
      end if

C***  Gather 4D array on atmospheric grid into the global array

      CALL PACK_BLOCK (aGRID, aA, aA_glob)

      CALL PACK_DATA (aGRID, aWEIGHT_loc, aWEIGHT)

C***  Interpolate aA_glob from atmospheric grid to ocean grid 

      if(AM_I_ROOT()) then

        call HNTR80 (aIM,aJM,0.d0,aDLATM, oIM,oJM,0.d0,oDLATM, 0.d0)

        DO N=1,NT
          aFtemp(:,:) = 0.
          aFtemp(:,:) = aA_glob(N,oN,:,:)
          aFtemp(2:aIM,aJM) = aFtemp(1,aJM)
          call HNTR8P (aWEIGHT, aFtemp, oFtemp)
          oA_glob(N,oN,:,:) = oFtemp(:,:)
          oB(:,:,N) = oA_glob(N,oN,:,:) 
        END DO
      end if

C***  Scatter global array oA_glob to the ocean grid

      CALL UNPACK_BLOCK (ogrid, oA_glob, oA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob, aFtemp,oFtemp, aWEIGHT)
      end if

      RETURN
      END SUBROUTINE INT_AG2OG_4Db

      SUBROUTINE INT_AG2OG_Vector1(aU,aV,oU,oV,aWEIGHT_loc,aFOCEAN_loc
     &     ,aN,oN) 

!@sum INT_AG2OG_Vector1 is for conversion vector from atm. A grid to ocean A grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM,oSINI=>SINIC,oCOSI=>COSIC
      USE GEOM,  only : aDLATM=>DLATM,aIMAXJ=>IMAXJ,aSINI=>SINIP,
     &                  aCOSI=>COSIP

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_DATA, UNPACK_DATA
      USE OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER, I,J, aN,oN

      REAL*8, INTENT(IN)  :: 
     *        aFOCEAN_loc(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(IN)  :: 
     *        aWEIGHT_loc(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(IN)  ::  
     *  aU(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN) 
     * ,aV(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN) 
      REAL*8, INTENT(OUT) ::
     *  oU(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN)
     * ,oV(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN)

      REAL*8  aUsp, aVsp, aUnp, aVnp
      REAL*8  oUsp, oVsp, oUnp, oVnp

      REAL*8, ALLOCATABLE :: aU_glob(:,:,:),aV_glob(:,:,:)
     *                      ,oU_glob(:,:,:),oV_glob(:,:,:)
     *                     , aFtemp(:,:),oFtemp(:,:)
     *                     , aWEIGHT(:,:), aFOCEAN(:,:)

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oU(1:aIM,1:aJM,oN) = aU(1:aIM,1:aJM,oN)
        oV(1:aIM,1:aJM,oN) = aV(1:aIM,1:aJM,oN)
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aU_glob(aIM,aJM,aN), STAT = IER)
        ALLOCATE(aV_glob(aIM,aJM,aN), STAT = IER)
        ALLOCATE(oU_glob(oIM,oJM,oN), STAT = IER)
        ALLOCATE(oV_glob(oIM,oJM,oN), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(aWEIGHT(aIM,aJM), STAT = IER)
        ALLOCATE(aFOCEAN(aIM,aJM), STAT = IER)
      end if

C***  Gather 2D array on atmospheric grid into the global array

      CALL PACK_DATA (aGRID, aU, aU_glob)
      CALL PACK_DATA (aGRID, aV, aV_glob)

      CALL PACK_DATA (aGRID, aWEIGHT_loc, aWEIGHT)

      CALL PACK_DATA (aGRID, aFOCEAN_loc, aFOCEAN)

C***  Interpolate aA_glob from atmospheric grid to ocean grid 

      if(AM_I_ROOT()) then

        aUsp = aU_glob(1,1,oN)
        aVsp = aV_glob(1,1,oN)
        aUnp = aU_glob(1,aJM,oN)
        aVnp = aV_glob(1,aJM,oN)

        aFtemp(:,1  ) = aUsp*aCOSI(:) - aVsp*aSINI(:)
        aFtemp(:,aJM) = aUnp*aCOSI(:) + aVnp*aSINI(:)

        DO J=2,aJM-1
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aU_glob(I,J,oN)                !  kg/m s
          END IF
        END DO
        END DO
        call HNTR80 (aIM,aJM,0.d0,aDLATM, oIM,oJM,0.d0,oDLATM, 0.d0)
        call HNTR8  (aWEIGHT, aFtemp, oU_glob)          !!  A-grid => A-grid

        oUsp =  SUM(oU_glob(:,  1,oN)*oCOSI(:))*2/oIM
        oVsp = -SUM(oU_glob(:,  1,oN)*oSINI(:))*2/oIM
        oUnp =  SUM(oU_glob(:,oJM,oN)*oCOSI(:))*2/oIM
        oVnp =  SUM(oU_glob(:,oJM,oN)*oSINI(:))*2/oIM

        oU_glob(1,  1,oN) = oUsp
        oU_glob(1,oJM,oN) = oUnp

        aFtemp(:,  1) = aVsp*aCOSI(:) + aUsp*aSINI(:)
        aFtemp(:,aJM) = aVnp*aCOSI(:) - aUnp*aSINI(:)

        DO J=2,aJM-1
          DO I=1,aIMAXJ(J)
            IF (aFOCEAN(I,J).gt.0.) THEN
              aFtemp(I,J) = aV_glob(I,J,oN)              !  kg/m s
            END IF
          END DO
        END DO
        call HNTR80 (aIM,aJM,0.d0,aDLATM, oIM,oJM,0.d0,oDLATM, 0.d0)
        call HNTR8  (aWEIGHT, aFtemp, oV_glob)          !!  A-grid => A-grid

        oV_glob(1,  1,oN) = oVsp
        oV_glob(1,oJM,oN) = oVnp
      end if

C***  Scatter global array oA_glob to the ocean grid

      CALL UNPACK_DATA (ogrid, oU_glob, oU)
      CALL UNPACK_DATA (ogrid, oV_glob, oV)

      if(AM_I_ROOT()) then
        DEALLOCATE(aU_glob,aV_glob, oU_glob,oV_glob
     *           , aFtemp,oFtemp, aWEIGHT, aFOCEAN)
      end if

      end if

      RETURN
      END SUBROUTINE INT_AG2OG_Vector1

      SUBROUTINE INT_AG2OG_Vector2(aU,aV, oU,oV, aWEIGHT_loc
     *                           , IVSPO,IVNPO)

!@sum INT_AG2OG_Vector2 is for conversion vector from atm. (U,V) grid to ocean (U,V) grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      USE GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_DATA, UNPACK_DATA
      USE OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER, I,J, IVSPO,IVNPO

      REAL*8, INTENT(IN)  :: 
     *        aWEIGHT_loc(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(IN)  ::  
     *  aU(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
     * ,aV(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(OUT) ::
     *  oU(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
     * ,oV(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, DIMENSION(aIM) :: aSINI,aCOSI
      REAL*8, DIMENSION(oIM) :: oSINI,oCOSI

      REAL*8  aUsp, aVsp, aUnp, aVnp
      REAL*8  oUsp, oVsp, oUnp, oVnp

      REAL*8, ALLOCATABLE :: aU_glob(:,:),aV_glob(:,:)
     *                      ,oU_glob(:,:),oV_glob(:,:)
     *                     , aFtemp(:,:),oFtemp(:,:), aWEIGHT(:,:)

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oU(1:aIM,1:aJM) = aU(1:aIM,1:aJM)
        oV(1:aIM,1:aJM) = aV(1:aIM,1:aJM)
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aU_glob(aIM,aJM), STAT = IER)
        ALLOCATE(aV_glob(aIM,aJM), STAT = IER)
        ALLOCATE(oU_glob(oIM,oJM), STAT = IER)
        ALLOCATE(oV_glob(oIM,oJM), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(aWEIGHT(aIM,aJM), STAT = IER)
      end if

C***  Gather 2D array on atmospheric grid into the global array

      CALL PACK_DATA (aGRID, aU, aU_glob)
      CALL PACK_DATA (aGRID, aV, aV_glob)

      CALL PACK_DATA (aGRID, aWEIGHT_loc, aWEIGHT)

C***  Interpolate aA_glob from atmospheric grid to ocean grid 

      if(AM_I_ROOT()) then

        aUsp = 0.0 
        aVsp = 0.0  
        aUnp = 0.0  
        aVnp = 0.0  

        aFtemp(:,:)   = aU_glob(:,:)
        aFtemp(:,aJM) = 0.0
        call HNTR80 (aIM,aJM,.5d0,aDLATM, oIM,oJM,.5d0,oDLATM, 0.d0)
        call HNTR8  (aWEIGHT, aFtemp, oU_glob)     !!  U-grid => U-grid


        oUsp = 0.0  
        oVsp = 0.0  
        oUnp = 0.0   
        oVnp = 0.0   

        oU_glob(  oIM,  1) = oUsp
        oU_glob(IVSPO,  1) = oVsp
        oU_glob(  oIM,oJM) = oUnp
        oU_glob(IVNPO,oJM) = oVnp

        call HNTR80 (aIM,aJM-1,0.d0,aDLATM, oIM,oJM-1,0.d0,oDLATM, 0.d0)
        call HNTR8  (aWEIGHT, aV_glob, oV_glob)      !!  V-grid => V-grid
        oV_glob(:,oJM) = 0.0                       !!  should not be used
      end if

C***  Scatter global array oA_glob to the ocean grid

      CALL UNPACK_DATA (ogrid, oU_glob, oU)
      CALL UNPACK_DATA (ogrid, oV_glob, oV)

      if(AM_I_ROOT()) then
        DEALLOCATE(aU_glob,aV_glob, oU_glob,oV_glob
     *           , aFtemp,oFtemp, aWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_AG2OG_Vector2

      SUBROUTINE INT_AG2OG_2Da_par(aA,oA,aWEIGHT,copypole)

!@sum INT_AG2OG_2Da_par parallel version of INT_AG2OG_2Da
!@auth M. Kelley

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm
      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM
      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, 
     &     hasNorthPole, hasSouthPole
      Use OCEANR_DIM,       only : oGRID
      IMPLICIT NONE

      real*8, dimension(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     intent(in)  :: aWEIGHT, aA
      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(out)  :: oA
      integer :: jmin,jmax
      real*8, dimension(:,:), allocatable :: aA_band,aWEIGHT_band
      logical, intent(in), optional :: copypole

      logical :: copypole_

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oA(1:aIM,1:aJM) = aA(1:aIM,1:aJM)
        return
      endif

      if(hntrp_a2o_needs_init) then
        call Init_Hntrp_Type(hntrp_a2o,
     &     aGRID, 0.d0,aDLATM,
     &     oGRID, 0.d0,oDLATM,
     &     0.d0)
        hntrp_a2o_needs_init = .false.
      endif

C***  Gather the requisite atm latitude bands
      jmin = hntrp_a2o%bpack%jband_strt
      jmax = hntrp_a2o%bpack%jband_stop
      ALLOCATE(aA_band(aIM,jmin:jmax),aWEIGHT_band(aIM,jmin:jmax))
      CALL BAND_PACK (hntrp_a2o%bpack, aA, aA_band)
      CALL BAND_PACK (hntrp_a2o%bpack, aWEIGHT, aWEIGHT_band)

C***  Interpolate aA from atmospheric grid to ocean grid 
      copypole_ = .true.
      if(present(copypole)) copypole_ = copypole
      if(copypole_ .and. hasNorthPole(aGRID)) then
        aA_band(2:aIM,aJM) = aA_band(1,aJM)
      endif
      if(copypole_) then
        call HNTR8P_band(aWEIGHT_band, aA_band, hntrp_a2o, oA)
      else
        call HNTR8_band (aWEIGHT_band, aA_band, hntrp_a2o, oA)
      endif
      DEALLOCATE(aA_band,aWEIGHT_band)

      RETURN
      END SUBROUTINE INT_AG2OG_2Da_par

      SUBROUTINE INT_AG2OG_3Db_par(aA,oA,aWEIGHT, aN,oN)
!@sum INT_AG2OG_3Db_par parallel version of INT_AG2OG_3Db
!@auth M. Kelley
      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID
      Use OCEANR_DIM,       only : oGRID
      IMPLICIT NONE
      integer, intent(in) :: aN,oN
      real*8, dimension(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN),
     &     intent(in)  :: aA
      real*8, dimension(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     intent(in)  :: aWEIGHT
      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN),
     &     intent(out)  :: oA
      call INT_AG2OG(aA(:,:,oN),oA(:,:,oN),aWEIGHT)
      RETURN
      END SUBROUTINE INT_AG2OG_3Db_par

      SUBROUTINE INT_AG2OG_3Dc_par(aA,oA,aWEIGHT, aN,oN,oNin)
!@sum INT_AG2OG_3Dc_par parallel version of INT_AG2OG_3Dc
!@auth M. Kelley
      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID
      Use OCEANR_DIM,       only : oGRID
      IMPLICIT NONE
      integer, intent(in) :: aN,oN,oNin
      real*8, dimension(aN,aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     intent(in)  :: aA
      real*8, dimension(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     intent(in)  :: aWEIGHT
      real*8, dimension(oN,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(out)  :: oA
      real*8, dimension(:,:), allocatable :: aA2d,oA2d

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oA(oNin,1:aIM,1:aJM) = aA(oNin,1:aIM,1:aJM)
      else
        allocate(aA2d(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
        allocate(oA2d(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
        aA2d = aA(oNin,:,:)
        call INT_AG2OG(aA2d,oA2d,aWEIGHT)
        oA(oNin,:,:) = oA2d
        deallocate(aA2d,oA2d)
      endif
      RETURN
      END SUBROUTINE INT_AG2OG_3Dc_par

      SUBROUTINE INT_AG2OG_Vector1_par(aU,aV,oU,oV,aWEIGHT,
     &     aFOCEAN, aN,oN)
!@sum INT_AG2OG_Vector1_par parallel version of INT_AG2OG_Vector1
!@auth M. Kelley
      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM,oSINI=>SINIC,oCOSI=>COSIC
      USE GEOM,  only : aDLATM=>DLATM,aIMAXJ=>IMAXJ,aSINI=>SINIP,
     &                  aCOSI=>COSIP

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID,
     &     hasNorthPole, hasSouthPole
      USE OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      integer, intent(in) :: aN,oN
      real*8, dimension(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     intent(in) :: aFOCEAN, aWEIGHT
      real*8, dimension(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN),
     &     intent(in) :: aU, aV
      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN),
     &     intent(out) :: oU, oV

C***  Local vars
      integer :: i,j
      real*8 :: aUnp,aVnp, oUnp,oVnp
      real*8, dimension(:,:), allocatable :: aUtmp,aVtmp

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        oU(1:aIM,1:aJM,oN) = aU(1:aIM,1:aJM,oN)
        oV(1:aIM,1:aJM,oN) = aV(1:aIM,1:aJM,oN)
        return
      endif


      allocate(aUtmp(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(aVtmp(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))

      do j=max(2,aGRID%J_STRT),min(aGRID%J_STOP,aJM-1)
        do i=1,aIMAXJ(J)
          if(aFOCEAN(i,j).gt.0.) then
            aUtmp(i,j) = aU(i,j,oN)
            aVtmp(i,j) = aV(i,j,oN)
          else
            aUtmp(i,j) = 0.
            aVtmp(i,j) = 0.
          endif
        enddo
      enddo
      if(hasSouthPole(aGRID)) then
        aUtmp(:,1) = 0.
        aVtmp(:,1) = 0.
      endif
      if(hasNorthPole(aGRID)) then
        aUnp = aU(1,aJM,oN)
        aVnp = aV(1,aJM,oN)
        aUtmp(:,aJM) = aUnp*aCOSI(:) + aVnp*aSINI(:)
        aVtmp(:,aJM) = aVnp*aCOSI(:) - aUnp*aSINI(:)
      endif
c      call INT_AG2OG(aUtmp,oU(:,:,oN),aWEIGHT,copypole=.false.)
c      call INT_AG2OG(aVtmp,oV(:,:,oN),aWEIGHT,copypole=.false.)
      call INT_AG2OG_2Da_par(aUtmp,oU(:,:,oN),aWEIGHT,copypole=.false.)
      call INT_AG2OG_2Da_par(aVtmp,oV(:,:,oN),aWEIGHT,copypole=.false.)
      if(hasNorthPole(oGRID)) then
        oUnp = SUM(oU(:,oJM,oN)*oCOSI(:))*2/oIM
        oVnp = SUM(oU(:,oJM,oN)*oSINI(:))*2/oIM
        oU(1,oJM,oN) = oUnp
        oV(1,oJM,oN) = oVnp
      endif

      DEALLOCATE(aUtmp,aVtmp)

      RETURN
      END SUBROUTINE INT_AG2OG_Vector1_par

      END MODULE INT_AG2OG_MOD

      MODULE INT_OG2AG_MOD

!@sum INT_OG2AG_MOD contains subroutines for conversion 2D, 3D, etc. 
!!    arrays from ocean to the atm. grid 
!@auth Larissa Nazarenko
!@ver  1.0
      use hntrp_mod
      IMPLICIT NONE
      SAVE
      PRIVATE
      PUBLIC INT_OG2AG

      Interface INT_OG2AG
c _par are parallelized versions
      Module Procedure INT_OG2AG_2Da_par
      Module Procedure INT_OG2AG_3Da_par
      Module Procedure INT_OG2AG_3Db_par
      Module Procedure INT_OG2AG_4Da ! need _par for tracers
c      Module Procedure INT_OG2AG_Vector1

      End Interface

      type(hntrp_type) :: hntrp_o2a ! instance for ocn A -> atm A
      logical :: hntrp_o2a_needs_init = .true.

      contains

      SUBROUTINE INT_OG2AG_2Da(oA,aA,oWEIGHT_loc, CopyPole, AvgPole)

!@sum INT_OG2AG_3D is for conversion 3D arrays from ocean to the atm. grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM, aIMAXJ=>IMAXJ

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_DATA, UNPACK_DATA
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: CopyPole
      LOGICAL, INTENT(IN), OPTIONAL :: AvgPole

      INTEGER :: IER

      REAL*8, INTENT(IN)  :: 
     *        oWEIGHT_loc(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) 
      REAL*8  :: 
     *  aA(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8 :: 
     *  oA(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, ALLOCATABLE :: aA_glob(:,:), oA_glob(:,:)
     *                     , aFtemp(:,:),oFtemp(:,:), oWEIGHT(:,:)

      LOGICAL :: AvgPole_

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   

        aA(1:aIM,1:aJM) = oA(1:aIM,1:aJM)
      else

      AvgPole_ = .true.
      if(present(AvgPole)) AvgPole_ = AvgPole

      if(AM_I_ROOT()) then
        ALLOCATE(aA_glob(aIM,aJM), STAT = IER)
        ALLOCATE(oA_glob(oIM,oJM), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(oWEIGHT(oIM,oJM), STAT = IER)
      end if

C***  Gather 3D array on ocean grid into the global array

      CALL PACK_DATA (oGRID, oA, oA_glob)

      CALL PACK_DATA (oGRID, oWEIGHT_loc, oWEIGHT)

C***  Interpolate aA_glob from ocean grid to atmospheric grid 

      if(AM_I_ROOT()) then

        call HNTR80 (oIM,oJM,0.d0,oDLATM, aIM,aJM,0.d0,aDLATM, 0.d0)

        IF (CopyPole) oWEIGHT(2:oIM,oJM) = oWEIGHT(1,oJM)
        oFtemp(:,:) = oA_glob(:,:)
        if(AvgPole_) then
          oFtemp(2:oIM,oJM) = oFtemp(1,oJM)
          call HNTR8P (oWEIGHT, oFtemp, aFtemp)
        else
          call HNTR8 (oWEIGHT, oFtemp, aFtemp)
        endif
        aA_glob(:,:) = aFtemp(:,:)  
      end if

C***  Scatter global array aA_glob to the atmospheric grid

      CALL UNPACK_DATA (agrid, aA_glob, aA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob, aFtemp,oFtemp, oWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_OG2AG_2Da

      SUBROUTINE INT_OG2AG_3Da(oA,aA,oWEIGHT_loc,NT)

!@sum INT_OG2AG_3D is for conversion 3D arrays from ocean to the atm. grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM, aIMAXJ=>IMAXJ

      USE MODEL_COM, only : aFOCEAN_loc=>FOCEAN

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_COLUMN, UNPACK_COLUMN
     *                           , PACK_DATA
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER, NT, I,J,N
      INTEGER :: aJ_0,aJ_1, aI_0,aI_1

      REAL*8, INTENT(IN)  :: 
     *        oWEIGHT_loc(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) 
      REAL*8  :: 
     *  aA(NT,aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8 :: 
     *  oA(NT,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, ALLOCATABLE :: aA_glob(:,:,:), oA_glob(:,:,:)
     *                     , aFOCEAN(:,:) 
     *                     , aFtemp(:,:),oFtemp(:,:), oWEIGHT(:,:)

      aJ_0 = aGRID%j_STRT
      aJ_1 = aGRID%j_STOP
      aI_0 = aGRID%I_STRT
      aI_1 = aGRID%I_STOP

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   

        DO N=1,NT
          DO J=aJ_0,aJ_1
            DO I=aI_0,aIMAXJ(J)
              IF (aFOCEAN_loc(I,J).gt.0.) THEN
                aA(N,I,J) = oA(N,I,J)
              END IF
            END DO
          END DO
        END DO
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aFOCEAN(aIM,aJM), STAT = IER)
        ALLOCATE(aA_glob(NT,aIM,aJM), STAT = IER)
        ALLOCATE(oA_glob(NT,oIM,oJM), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(oWEIGHT(oIM,oJM), STAT = IER)
      end if

C***  Gather 3D array on both atmospheric and ocean grids into the global array

      CALL PACK_COLUMN (aGRID, aA, aA_glob)
      CALL PACK_COLUMN (oGRID, oA, oA_glob)

      CALL PACK_DATA  (aGRID, aFOCEAN_loc, aFOCEAN)

      CALL PACK_DATA (oGRID, oWEIGHT_loc, oWEIGHT)

C***  Interpolate aA_glob from ocean grid to atmospheric grid 

      if(AM_I_ROOT()) then

        call HNTR80 (oIM,oJM,0.d0,oDLATM, aIM,aJM,0.d0,aDLATM, 0.d0)

        DO N=1,NT
          oFtemp(:,:) = oA_glob(N,:,:)
          oFtemp(2:oIM,oJM) = oFtemp(1,oJM)
          call HNTR8P (oWEIGHT, oFtemp, aFtemp)
          DO J=1,aJM
            DO I=1,aIMAXJ(J)
              IF (aFOCEAN(I,J).gt.0.) THEN
                aA_glob(N,I,J) = aFtemp(I,J)  
              END IF
            END DO
          END DO
        END DO
      end if

C***  Scatter global array aA_glob to the atmospheric grid

      CALL UNPACK_COLUMN (agrid, aA_glob, aA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob, aFOCEAN, aFtemp,oFtemp, oWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_OG2AG_3Da

      SUBROUTINE INT_OG2AG_3Db(oA,aA,oWEIGHT_loc, oN,aN,CopyPole)

!@sum INT_OG2AG_3D is for conversion 3D arrays from ocean to the atm. grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM, aIMAXJ=>IMAXJ

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_DATA, UNPACK_DATA
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: CopyPole

      INTEGER :: IER, N, oN,aN

      REAL*8, INTENT(IN)  :: 
     *        oWEIGHT_loc(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) 
      REAL*8  :: 
     *  aA(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN) 
      REAL*8 :: 
     *  oA(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN)

      REAL*8, ALLOCATABLE :: aA_glob(:,:,:), oA_glob(:,:,:)
     *                     , aFtemp(:,:),oFtemp(:,:), oWEIGHT(:,:)

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   

        DO N=1,aN
          aA(1:aIM,1:aJM,N) = oA(1:aIM,1:aJM,N)
        END DO
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aA_glob(aIM,aJM,aN), STAT = IER)
        ALLOCATE(oA_glob(oIM,oJM,oN), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(oWEIGHT(oIM,oJM), STAT = IER)
      end if

C***  Gather 3D array on ocean grid into the global array

      CALL PACK_DATA (oGRID, oA, oA_glob)

      CALL PACK_DATA (oGRID, oWEIGHT_loc, oWEIGHT)

C***  Interpolate aA_glob from ocean grid to atmospheric grid 

      if(AM_I_ROOT()) then

        call HNTR80 (oIM,oJM,0.d0,oDLATM, aIM,aJM,0.d0,aDLATM, 0.d0)

        IF (CopyPole) oWEIGHT(2:oIM,oJM) = oWEIGHT(1,oJM)
        DO N=1,aN
          oFtemp(:,:) = oA_glob(:,:,N)
          oFtemp(2:oIM,oJM) = oFtemp(1,oJM)
          call HNTR8P (oWEIGHT, oFtemp, aFtemp)
          aA_glob(:,:,N) = aFtemp(:,:)  
        END DO
      end if

C***  Scatter global array aA_glob to the atmospheric grid

      CALL UNPACK_DATA (agrid, aA_glob, aA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob, aFtemp,oFtemp, oWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_OG2AG_3Db

      SUBROUTINE INT_OG2AG_4Da(oA,aA,oWEIGHT_loc,NT,NTM)

!@sum INT_OG2AG_4D is for conversion 4D arrays from ocean to the atm. grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM, aIMAXJ=>IMAXJ

      USE MODEL_COM, only : aFOCEAN_loc=>FOCEAN

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_BLOCK, UNPACK_BLOCK
     *                           , PACK_DATA
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER, NTM,NT,N1,N2, I,J,N
      INTEGER :: aJ_0,aJ_1, aI_0,aI_1

      REAL*8, INTENT(IN)  :: 
     *        oWEIGHT_loc(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) 
      REAL*8  :: 
     *  aA(NTM,NT,aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8 :: 
     *  oA(NTM,NT,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, ALLOCATABLE :: aA_glob(:,:,:,:), oA_glob(:,:,:,:)
     *                     , aFOCEAN(:,:) 
     *                     , aFtemp(:,:),oFtemp(:,:), oWEIGHT(:,:)

      aJ_0 = aGRID%j_STRT
      aJ_1 = aGRID%j_STOP
      aI_0 = aGRID%I_STRT
      aI_1 = aGRID%I_STOP

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   

        DO N1=1,NTM
          DO N2=1,NT
          DO J=aJ_0,aJ_1
            DO I=aI_0,aIMAXJ(J)
              IF (aFOCEAN_loc(I,J).gt.0.) THEN
                aA(N1,N2,I,J) = oA(N1,N2,I,J)
              END IF
            END DO
          END DO
          END DO
        END DO
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aFOCEAN(aIM,aJM), STAT = IER)
        ALLOCATE(aA_glob(NTM,NT,aIM,aJM), STAT = IER)
        ALLOCATE(oA_glob(NTM,NT,oIM,oJM), STAT = IER)
        ALLOCATE(aFtemp(aIM,aJM), STAT = IER)
        ALLOCATE(oFtemp(oIM,oJM), STAT = IER)
        ALLOCATE(oWEIGHT(oIM,oJM), STAT = IER)
      end if

C***  Gather 3D array on ocean grid into the global array

      CALL PACK_BLOCK (aGRID, aA, aA_glob)
      CALL PACK_BLOCK (oGRID, oA, oA_glob)

      CALL PACK_DATA  (agrid, aFOCEAN_loc, aFOCEAN)

      CALL PACK_DATA (oGRID, oWEIGHT_loc, oWEIGHT)

C***  Interpolate aA_glob from ocean grid to atmospheric grid 

      if(AM_I_ROOT()) then

        call HNTR80 (oIM,oJM,0.d0,oDLATM, aIM,aJM,0.d0,aDLATM, 0.d0)

        DO N1=1,NTM
        DO N2=1,NT
          oFtemp(:,:) = oA_glob(N1,N2,:,:)
          oFtemp(2:oIM,oJM) = oFtemp(1,oJM)
          call HNTR8P (oWEIGHT, oFtemp, aFtemp)
          DO J=1,aJM
            DO I=1,aIMAXJ(J)
              IF (aFOCEAN(I,J).gt.0.) THEN
                aA_glob(N1,N2,I,J) = aFtemp(I,J)  
              END IF
            END DO
          END DO
        END DO
        END DO
      end if

C***  Scatter global array aA_glob to the atm. grid

      CALL UNPACK_BLOCK (agrid, aA_glob, aA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob, aFOCEAN, aFtemp,oFtemp, oWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_OG2AG_4Da

      SUBROUTINE INT_OG2AG_Vector1(oUO1,oVO1,aUO1,aVO1, oWEIGHT_loc 
     *                           , IVSPO,IVNPO) 

!@sum INT_OG2AG_3D is for conversion 3D arrays from ocean to the atm. grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM, oCOSU=>COSU,oSINU=>SINU
      Use GEOM,  only : aDLATM=>DLATM, aIMAXJ=>IMAXJ
     *                , aCOSI=>COSIP,aSINI=>SINIP

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_DATA, UNPACK_DATA
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER, IVSPO,IVNPO

      REAL*8  aUsp, aVsp, aUnp, aVnp
      REAL*8  oUsp, oVsp, oUnp, oVnp

      REAL*8, INTENT(IN)  :: 
     *        oWEIGHT_loc(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) 
      REAL*8  :: 
     *  aUO1(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8  :: 
     *  aVO1(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8 :: 
     *  oUO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8 :: 
     *  oVO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, ALLOCATABLE :: aUO1_glob(:,:),aVO1_glob(:,:)
     *                     , oUO1_glob(:,:),oVO1_glob(:,:)
     *                     , oWEIGHT(:,:)

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   

        aUO1(1:aIM,1:aJM) = oUO1(1:aIM,1:aJM)
        aVO1(1:aIM,1:aJM) = oVO1(1:aIM,1:aJM)
      else

      if(AM_I_ROOT()) then
        ALLOCATE(aUO1_glob(aIM,aJM), STAT = IER)
        ALLOCATE(aVO1_glob(aIM,aJM), STAT = IER)
        ALLOCATE(oUO1_glob(oIM,oJM), STAT = IER)
        ALLOCATE(oVO1_glob(oIM,oJM), STAT = IER)
        ALLOCATE(oWEIGHT(oIM,oJM), STAT = IER)
      end if

C***  Gather 3D array on ocean grid into the global array

      CALL PACK_DATA (oGRID, oUO1, oUO1_glob)
      CALL PACK_DATA (oGRID, oVO1, oVO1_glob)

      CALL PACK_DATA (oGRID, oWEIGHT_loc, oWEIGHT)

C***  Interpolate aA_glob from ocean grid to atmospheric grid 

      if(AM_I_ROOT()) then

!!!  U velocity for the 1st ocean layer.

        oVsp = oUO1_glob(IVSPO,  1)
        oVnp = oUO1_glob(IVNPO,oJM)

        oUO1_glob(:,  1) = oUO1_glob(oIM,  1)*oCOSU(:) - oVsp*oSINU(:)
        oUO1_glob(:,oJM) = oUO1_glob(oIM,oJM)*oCOSU(:) + oVnp*oSINU(:)

        call HNTR80 (oIM,oJM,0.5d0,oDLATM, aIM,aJM,0.d0,aDLATM, 0.d0)
        call HNTR8  (oWEIGHT, oUO1_glob, aUO1_glob)    !!  U-grid => A-grid

        aUsp = SUM(aUO1_glob(:,  1)*aCOSI(:))*2/aIM
        aVsp = SUM(aUO1_glob(:,  1)*aSINI(:))*2/aIM
        aUnp = SUM(aUO1_glob(:,aJM)*aCOSI(:))*2/aIM
        aVnp = SUM(aUO1_glob(:,aJM)*aSINI(:))*2/aIM

        aUO1_glob(1,  1) = aUsp
        aUO1_glob(1,aJM) = aUnp

!!!  V velocity for the 1st ocean layer

        call HNTR80 (oIM,oJM-1,0.d0,oDLATM, aIM,aJM,0.d0,aDLATM, 0.d0)
        call HNTR8  (oWEIGHT, oVO1_glob, aVO1_glob)      !!  V-grid => A-grid

        aVO1_glob(1,  1) = aVsp
        aVO1_glob(1,aJM) = aVnp
      end if

C***  Scatter global arrays to the atmospheric grid

      CALL UNPACK_DATA (agrid, aUO1_glob, aUO1)
      CALL UNPACK_DATA (agrid, aVO1_glob, aVO1)

      if(AM_I_ROOT()) then
        DEALLOCATE(aUO1_glob,aVO1_glob, oUO1_glob,oVO1_glob, oWEIGHT)
      end if

      end if

      RETURN
      END SUBROUTINE INT_OG2AG_Vector1

      SUBROUTINE INT_OG2AG_2Da_par(oA,aA,oWEIGHT, CopyPole, AvgPole)
!@sum INT_OG2AG_2Da_par parallel version of INT_OG2AG_2Da
!@auth M. Kelley

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID,
     &     hasNorthPole, hasSouthPole
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: CopyPole
      LOGICAL, INTENT(IN), OPTIONAL :: AvgPole

      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(in)  :: oWEIGHT, oA
      real*8, dimension(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     intent(out)  :: aA

      integer :: jmin,jmax
      real*8, dimension(:,:), allocatable :: oA_band,oWEIGHT_band
      logical :: AvgPole_

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        aA(1:aIM,1:aJM) = oA(1:aIM,1:aJM)
        return
      endif

      if(hntrp_o2a_needs_init) then
        call Init_Hntrp_Type(hntrp_o2a,
     &     oGRID, 0.d0,oDLATM,
     &     aGRID, 0.d0,aDLATM,
     &     0.d0)
        hntrp_o2a_needs_init = .false.
      endif

      AvgPole_ = .true.
      if(present(AvgPole)) AvgPole_ = AvgPole

      jmin = hntrp_o2a%bpack%jband_strt
      jmax = hntrp_o2a%bpack%jband_stop
      ALLOCATE(oA_band(oIM,jmin:jmax),oWEIGHT_band(oIM,jmin:jmax))

C***  Gather the requisite ocean latitude bands
      CALL BAND_PACK (hntrp_o2a%bpack, oA, oA_band)
      CALL BAND_PACK (hntrp_o2a%bpack, oWEIGHT, oWEIGHT_band)

C***  Interpolate oA_band from ocean grid to atmospheric grid 

      if(hasNorthPole(oGRID)) then
        if(CopyPole) oWEIGHT_band(2:oIM,oJM) = oWEIGHT_band(1,oJM)
        if(AvgPole_) oA_band(2:oIM,oJM) = oA_band(1,oJM)
      endif
      if(AvgPole_) then
        call HNTR8P_band(oWEIGHT_band, oA_band, hntrp_o2a, aA)
      else
        call HNTR8_band(oWEIGHT_band, oA_band, hntrp_o2a, aA)
      endif
      DEALLOCATE(oA_band, oWEIGHT_band)

      RETURN
      END SUBROUTINE INT_OG2AG_2Da_par

      SUBROUTINE INT_OG2AG_3Da_par(oA,aA,oWEIGHT,NT)
!@sum INT_OG2AG_3Da_par parallel version of INT_OG2AG_3Da
!@auth M. Kelley
      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM, aIMAXJ=>IMAXJ

      USE MODEL_COM, only : aFOCEAN_loc=>FOCEAN

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID,
     &     hasNorthPole, hasSouthPole
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NT

      real*8, dimension(NT,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(in)  :: oA
      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(in)  :: oWEIGHT
      real*8, dimension(NT,aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
     &     ,intent(out) :: aA

      real*8, dimension(:,:,:), allocatable :: oA_band
      real*8, dimension(:,:), allocatable :: oWEIGHT_band,oA2d,aA2d
      integer :: i,j,n, jmin,jmax

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        DO J=aGRID%j_STRT,aGRID%j_STOP
          DO I=1,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              DO N=1,NT
                aA(N,I,J) = oA(N,I,J)
              END DO
            END IF
          END DO
        END DO
        return
      endif

      if(hntrp_o2a_needs_init) then
        call Init_Hntrp_Type(hntrp_o2a,
     &     oGRID, 0.d0,oDLATM,
     &     aGRID, 0.d0,aDLATM,
     &     0.d0)
        hntrp_o2a_needs_init = .false.
      endif

      allocate(aA2d(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))

      jmin = hntrp_o2a%bpack%jband_strt
      jmax = hntrp_o2a%bpack%jband_stop
      ALLOCATE(oA_band(NT,oIM,jmin:jmax),oWEIGHT_band(oIM,jmin:jmax))
      allocate(oA2d(oIM,jmin:jmax))

C***  Gather the requisite ocean latitude bands
      CALL BAND_PACK_COLUMN (hntrp_o2a%bpack, oA, oA_band)
      CALL BAND_PACK (hntrp_o2a%bpack, oWEIGHT, oWEIGHT_band)

C***  Interpolate oA from ocean grid to atmospheric grid 
      do n=1,NT
        oA2d = oA_band(n,:,:)
        if(hasNorthPole(oGRID)) then
          oA2D(2:oIM,oJM) = oA2D(1,oJM)
        endif
        call HNTR8P_band(oWEIGHT_band, oA2D, hntrp_o2a, aA2d)
        do j=aGRID%j_STRT,aGRID%j_STOP
          do i=1,aIMAXJ(j)
            if (afocean_loc(i,j).gt.0.) then
              aA(n,i,j) = aA2d(i,j)
            endif
          enddo
        enddo
      enddo

      DEALLOCATE(oA_band, oWEIGHT_band, oA2d,aA2d)

      RETURN
      END SUBROUTINE INT_OG2AG_3Da_par

      SUBROUTINE INT_OG2AG_3Db_par(oA,aA,oWEIGHT,  oN,aN, CopyPole)
!@sum INT_OG2AG_3Db_par parallel version of INT_OG2AG_3Db
!@auth M. Kelley
      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID,
     &     hasNorthPole, hasSouthPole
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: oN,aN
      LOGICAL, INTENT(IN) :: CopyPole

      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN),
     &     intent(in)  :: oA
      real*8, dimension(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     intent(in)  :: oWEIGHT
      real*8, dimension(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN),
     &     intent(out)  :: aA

      real*8, dimension(:,:,:), allocatable :: oA_band
      real*8, dimension(:,:), allocatable :: oWEIGHT_band
      integer :: n, jmin,jmax

      if (oIM .eq. aIM .and. oJM .eq. aJM) then   
        aA(1:aIM,1:aJM,1:aN) = oA(1:aIM,1:aJM,1:aN)
        return
      endif

      if(hntrp_o2a_needs_init) then
        call Init_Hntrp_Type(hntrp_o2a,
     &     oGRID, 0.d0,oDLATM,
     &     aGRID, 0.d0,aDLATM,
     &     0.d0)
        hntrp_o2a_needs_init = .false.
      endif

      jmin = hntrp_o2a%bpack%jband_strt
      jmax = hntrp_o2a%bpack%jband_stop
      ALLOCATE(oA_band(oIM,jmin:jmax,aN),oWEIGHT_band(oIM,jmin:jmax))

C***  Gather the requisite ocean latitude bands
      CALL BAND_PACK (hntrp_o2a%bpack, oA(:,:,1:aN), oA_band)
      CALL BAND_PACK (hntrp_o2a%bpack, oWEIGHT, oWEIGHT_band)
      if(CopyPole .and. hasNorthPole(oGRID)) then
        oWEIGHT_band(2:oIM,oJM) = oWEIGHT_band(1,oJM)
      endif

C***  Interpolate oA from ocean grid to atmospheric grid 
      do n=1,aN
        if(hasNorthPole(oGRID)) then
          oA_band(2:oIM,oJM,n) = oA_band(1,oJM,n)
        endif
        call HNTR8P_band(oWEIGHT_band,oA_band(:,:,n),hntrp_o2a,
     &       aA(:,:,n))
      enddo

      DEALLOCATE(oA_band, oWEIGHT_band)

      RETURN
      END SUBROUTINE INT_OG2AG_3Db_par

      END MODULE INT_OG2AG_MOD

      Subroutine HNTR80 (IMA,JMA,OFFIA,DLATA,
     *                   IMB,JMB,OFFIB,DLATB, DATMIS)
C****
C**** HNTR80 fills in the common block HNTRCB with coordinate
C**** parameters that will be used by subsequent calls to HNTR8.
C**** The 5 Real input values are expected to be Real*8.
C****
C**** Input: IMA = number of cells in east-west direction of grid A
C****        JMA = number of cells in north-south direction of grid A
C****      OFFIA = number of cells of grid A in east-west direction
C****              from IDL (180) to western edge of cell IA=1
C****      DLATA = minutes of latitude for non-polar cells on grid A
C****        IMB = number of cells in east-west direction of grid B
C****        JMB = number of cells in north-south direction of grid B
C****      OFFIB = number of cells of grid B in east-west direction
C****              from IDL (180) to western edge of cell IB=1
C****      DLATB = minutes of latitude for non-polar cells on grid B
C****     DATMIS = missing data value inserted in output array B when
C****              cell (IB,JB) has integrated value 0 of WTA
C****
C**** Output: common block /HNTRCB/
C**** SINA(JA) = sine of latitude of northern edge of cell JA on grid A
C**** SINB(JB) = sine of latitude of northern edge of cell JB on grid B
C**** FMIN(IB) = fraction of cell IMIN(IB) on grid A west of cell IB
C**** FMAX(IB) = fraction of cell IMAX(IB) on grid A east of cell IB
C**** GMIN(JB) = fraction of cell JMIN(JB) on grid A south of cell JB
C**** GMAX(JB) = fraction of cell JMAX(JB) on grid A north of cell JB
C**** IMIN(IB) = western most cell of grid A that intersects cell IB
C**** IMAX(IB) = eastern most cell of grid A that intersects cell IB
C**** JMIN(JB) = southern most cell of grid A that intersects cell JB
C**** JMAX(JB) = northern most cell of grid A that intersects cell JB
C****
      Implicit Real*8 (A-H,O-Z)
      Parameter (TWOPI=6.283185307179586477d0)
      Real*8 OFFIA,DLATA, OFFIB,DLATB, DATMIS,DATMCB
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMCB, INA,JNA, INB,JNB
C****
      INA = IMA  ;  JNA = JMA
      INB = IMB  ;  JNB = JMB
      DATMCB = DATMIS
      If (IMA<1 .or. IMA>10800 .or. JMA<1 .or. JMA>5401 .or.
     *    IMB<1 .or. IMB>10800 .or. JMB<1 .or. JMB>5401)  GoTo 400
C****
C**** Partitions in east-west (I) direction
C**** Domain, around the globe, is scaled to fit from 0 to IMA*IMB
C****
      DIA = IMB  !  width of single A grid cell in scaled domain
      DIB = IMA  !  width of single B grid cell in scaled domain
      IA  = 1
      RIA = (IA+OFFIA - IMA)*IMB  !  scaled longitude of eastern edge
      IB  = IMB
      Do 150 IBp1=1,IMB
      RIB = (IBp1-1+OFFIB)*IMA    !  scaled longitude of eastern edge
  110 If (RIA-RIB)  120,130,140
  120 IA  = IA  + 1
      RIA = RIA + DIA
      GoTo 110
C**** Eastern edges of cells IA of grid A and IB of grid B coincide
  130 IMAX(IB) = IA
      FMAX(IB) = 0
      IA  = IA  + 1
      RIA = RIA + DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 0
      GoTo 150
C**** Cell IA of grid A contains western edge of cell IB of grid B
  140 IMAX(IB) = IA
      FMAX(IB) = (RIA-RIB)/DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 1-FMAX(IB)
  150 IB = IBp1
      IMAX(IMB) = IMAX(IMB) + IMA
C       WRITE (0,915) 'IMIN=',IMIN(1:IMB)
C       WRITE (0,915) 'IMAX=',IMAX(1:IMB)
C       WRITE (0,916) 'FMIN=',FMIN(1:IMB)
C       WRITE (0,916) 'FMAX=',FMAX(1:IMB)
C****
C**** Partitions in the north-south (J) direction
C**** Domain is measured in minutes (1/60-th of a degree)
C****
      FJEQA = .5*(1+JMA)
      Do 210 JA=1,JMA-1
      RJA = (JA+.5-FJEQA)*DLATA  !  latitude in minutes of northern edge
  210 SINA(JA) = Sin (RJA*TWOPI/(360*60))
      SINA(0)  = -1
      SINA(JMA)=  1
C****
      FJEQB = .5*(1+JMB)
      Do 220 JB=1,JMB-1
      RJB = (JB+.5-FJEQB)*DLATB  !  latitude in minutes of northern edge
  220 SINB(JB) = Sin (RJB*TWOPI/(360*60))
      SINB(0)  = -1
      SINB(JMB)=  1
C****
      JMIN(1) = 1
      GMIN(1) = 0
      JA = 1
      Do 350 JB=1,JMB-1
  310 If (SINA(JA)-SINB(JB))  320,330,340
  320 JA = JA + 1
      GoTo 310
C**** Northern edges of cells JA of grid A and JB of grid B coincide
  330 JMAX(JB) = JA
      GMAX(JB) = 0
      JA = JA + 1
      JMIN(JB+1) = JA
      GMIN(JB+1) = 0
      GoTo 350
C**** Cell JA of grid A contains northern edge of cell JB of grid B
  340 JMAX(JB) = JA
      GMAX(JB) = SINA(JA) - SINB(JB)
      JMIN(JB+1) = JA
      GMIN(JB+1) = SINB(JB) - SINA(JA-1)
  350 Continue
      JMAX(JMB) = JMA
      GMAX(JMB) = 0
C       WRITE (0,915) 'JMIN=',JMIN(1:JMB)
C       WRITE (0,915) 'JMAX=',JMAX(1:JMB)
C       WRITE (0,916) 'GMIN=',GMIN(1:JMB)
C       WRITE (0,916) 'GMAX=',GMAX(1:JMB)
      Return
C****
C**** Invalid parameters or dimensions out of range
C****
  400 Write (0,940) IMA,JMA,OFFIA,DLATA, IMB,JMB,OFFIB,DLATB, DATMIS
      Stop 400
C****
C 915 Format (/ 1X,A5 / (20I6))
C 916 Format (/ 1X,A5 / (20F6.2))
  940 Format ('0Arguments received by HNTRP0 in order:'/
     *   2I12,' = IMA,JMA = array dimensions for A grid'/
     *  E24.8,' = OFFIA   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATA   = minutes of latitude for interior grid cell'/
     *   2I12,' = IMB,JMB = array dimensions for B grid'/
     *  E24.8,' = OFFIB   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATB   = minute of latitude for interior grid cell'/
     *  E24.8,' = DATMIS  = missing data value to be put in B array',
     *                    ' when integrated WTA = 0'/
     *  '0These arguments are invalid or out of range.')
      End Subroutine HNTR80

      Subroutine HNTR8 (WTA,A,B)
C****
C**** HNTR8 performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.  B grid values that cannot be calculated because the
C**** covering A grid boxes have WTA = 0, are set to the value DATMIS.
C**** The area weighted integral of the quantity is conserved.
C**** The 3 Real input values are expected to be Real*8.
C****
C**** Input: WTA = weighting array for values on the A grid
C****          A = per unit area or per unit mass quantity
C**** Output:  B = horizontally interpolated quantity on B grid
C****
      Implicit Real*8 (A-H,O-Z)
      Real*8 WTA(*), A(*), B(*), DATMIS
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMIS, IMA,JMA, IMB,JMB
C****
C**** Interpolate the A grid onto the B grid
C****
      Do 20 JB=1,JMB
      JAMIN = JMIN(JB)
      JAMAX = JMAX(JB)
      Do 20 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      WEIGHT= 0
      VALUE = 0
      IAMIN = IMIN(IB)
      IAMAX = IMAX(IB)
      Do 10 JA=JAMIN,JAMAX
      G = SINA(JA)-SINA(JA-1)
      If (JA==JAMIN)  G = G - GMIN(JB)
      If (JA==JAMAX)  G = G - GMAX(JB)
      Do 10 IAREV=IAMIN,IAMAX
      IA  = 1 + Mod(IAREV-1,IMA)
      IJA = IA + IMA*(JA-1)
      F   = 1
      If (IAREV==IAMIN)  F = F - FMIN(IB)
      If (IAREV==IAMAX)  F = F - FMAX(IB)
      WEIGHT = WEIGHT + F*G*WTA(IJA)
   10 VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
      B(IJB) = DATMIS
      If (WEIGHT.ne.0)  B(IJB) = VALUE/WEIGHT
   20 Continue
      Return
      End Subroutine HNTR8

      Subroutine HNTR8P (WTA,A,B)
C****
C**** HNTR8P is similar to HNTR8 but polar values are replaced by
C**** their longitudinal mean.
C**** The 3 Real input values are expected to be Real*8.
C****
      Implicit Real*8 (A-H,O-Z)
      Real*8 WTA(*), A(*), B(*), DATMIS
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMIS, IMA,JMA, IMB,JMB
C****
      Call HNTR8 (WTA,A,B)
C****
C**** Replace individual values near the poles by longitudinal mean
C****
      Do 40 JB=1,JMB,JMB-1
      BMEAN  = DATMIS
      WEIGHT = 0
      VALUE  = 0
      Do 10 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      If (B(IJB) == DATMIS)  GoTo 10
      WEIGHT = WEIGHT + 1
      VALUE  = VALUE  + B(IJB)
   10 Continue
      If (WEIGHT.ne.0)  BMEAN = VALUE/WEIGHT
   20 Do 30 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
   30 B(IJB) = BMEAN
   40 Continue
      Return
      End Subroutine HNTR8P

