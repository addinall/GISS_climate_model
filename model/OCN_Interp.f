#include "rundeck_opts.h"

!#define AG2OG_PRECIP_BUNDLE
!#define OG2AG_TOC2SST_BUNDLE
!#define AG2OG_OCEANS_BUNDLE 
!#define OG2AG_OCEANS_BUNDLE

!#define BUNDLE_INTERP

#ifdef BUNDLE_INTERP

#ifdef CUBED_SPHERE
      subroutine bundle_interpolation(lstr, remap)
!@auth Igor Aleinov, Denis Gueyffier, Maxwell Kelley 
      USE cs2ll_utils, only : xgridremap_type,xgridremap_lij
      USE array_bundle, only : lookup_str, ab_bundle, ab_unbundle
      implicit none
      type (lookup_str) :: lstr
      type (xgridremap_type) :: remap
c
      real*8, pointer :: buf_s(:,:,:), buf_d(:,:,:)

      call ab_bundle( lstr, buf_s, buf_d )

      call xgridremap_lij( remap, buf_s, buf_d)

      call ab_unbundle( lstr, buf_s, buf_d )

      end subroutine bundle_interpolation
#else
      subroutine bundle_interpolation(lstr, htype)
!@auth Igor Aleinov, Denis Gueyffier, Maxwell Kelley 
      USE hntrp_mod, only : hntrp_type,band_pack_column,hntr8_band_lij
      USE array_bundle, only : lookup_str,get_bounds,ab_copy,
     &     ab_bundle,ab_unbundle
      implicit none
      type (lookup_str) :: lstr
      type (hntrp_type) :: htype
c
      real*8, pointer :: buf_s(:,:,:), buf_d(:,:,:), buf_band(:,:,:)

      integer :: si_0, si_1, sj_0, sj_1,
     &           di_0, di_1, dj_0, dj_1, lm, jmin,jmax

      call get_bounds(lstr, si_0, si_1, sj_0, sj_1,
     &                      di_0, di_1, dj_0, dj_1)

      if (dj_0 .eq. sj_0 .and. dj_1 .eq. sj_1) then   
         call ab_copy(lstr)
         return
      endif

      call ab_bundle( lstr, buf_s, buf_d )

      lm = size(buf_s,1)
      jmin = htype%bpack%jband_strt
      jmax = htype%bpack%jband_stop
      allocate(buf_band(lm,htype%ima,jmin:jmax))
      call band_pack_column(htype%bpack, buf_s, buf_band)
      call hntr8_band_lij(buf_band, htype, buf_d)
      deallocate(buf_band)

      call ab_unbundle( lstr, buf_s, buf_d )

      end subroutine bundle_interpolation
#endif

      subroutine lon_avg(arr,IM)
!@sum longitudinal average
      implicit none
      real*8 :: arr(1:IM),asum
      integer :: IM

      asum=sum(arr)
      arr(:)=asum/real(IM)

      end subroutine lon_avg

      subroutine copy_pole(arr,IM)
      implicit none
      real*8 :: arr(1:IM)
      integer :: IM
      
      arr(:)=arr(1)

      end subroutine copy_pole

#endif /* BUNDLE_INTERP */

#ifdef AG2OG_PRECIP_BUNDLE
      SUBROUTINE AG2OG_precip
!@sum  INT_AG2OG_precip is for interpolation of precipitation
!!       arrays from atmospheric grid to the ocean grid
!@auth Larissa Nazarenko, Denis Gueyffier
!@ver  1.0

      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE OCEAN,      only : oIM=>im, oJM=>jm

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE OCEANR_DIM, only : ogrid

      USE OCEAN, only : oDXYPO=>DXYPO,oDLATM=>DLATM, OXYP
      Use GEOM,  only : AXYP

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : NTM
#ifdef TRACERS_WATER
      USE FLUXES,  only : aTRPREC=>TRPREC, aTRUNPSI=>TRUNPSI
      USE OFLUXES, only : oTRPREC, oTRUNPSI
#endif
#endif
      USE SEAICE_COM, only : aRSI=>RSI
      USE FLUXES, only : aPREC=>PREC, aEPREC=>EPREC
     *     , aRUNPSI=>RUNPSI, aSRUNPSI=>SRUNPSI, aERUNPSI=>ERUNPSI

      USE OFLUXES, only : oRSI, oPREC, oEPREC
     *     , oRUNPSI, oSRUNPSI, oERUNPSI
      use ocean, only : remap_a2o
      USE array_bundle
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole

      IMPLICIT NONE
      INTEGER N
      REAL*8, allocatable :: aWEIGHT(:,:),oWEIGHT(:,:)
      REAL*8, allocatable :: aWEIGHT1(:,:),oWEIGHT1(:,:)
      REAL*8, allocatable :: atmp(:,:,:)
      type (lookup_str) :: lstr
      integer :: i,j,l

      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      !-- need to allocate separate array for different weights
      allocate(aweight1(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight1(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      !-- initializing "lstr"
      call ab_init( lstr, aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO,
     &     oGRID%I_STRT_HALO, oGRID%I_STOP_HALO,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO )

      aWEIGHT(:,:) = 1.-aRSI(:,:) !!  open ocean fraction
      call ab_add(lstr, aWEIGHT, oWEIGHT, shape(aWEIGHT), 'ij' )

      if (hasNorthPole(agrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         call copy_pole(aPREC(:,aJM),aIM)
         call copy_pole(aEPREC(:,aJM),aIM)
         call copy_pole(aRUNPSI(:,aJM),aIM)
         call copy_pole(aSRUNPSI(:,aJM),aIM)
         call copy_pole(aERUNPSI(:,aJM),aIM)
         call copy_pole(aRSI(:,aJM),aIM)
      endif

      call ab_add(lstr, aPREC, oPREC, shape(aPREC), 'ij',
     &     aWEIGHT, oWEIGHT )

      call ab_add(lstr, aEPREC, oEPREC, shape(aEPREC), 'ij',
     &     aWEIGHT, oWEIGHT)
#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)
      !-- allocate tmp array for trprec
      allocate(atmp(ntm,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))

      do N=1,NTM
        atmp(N,:,:) = aTRPREC(N,:,:)/AXYP(:,:)
        if (hasNorthPole(agrid)
     &       .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &       call copy_pole(atmp(N,:,aJM),aIM)
      enddo
     
      call ab_add(lstr, atmp, oTRPREC, shape(atmp), 'lij', 
     &     aWEIGHT, oWEIGHT )
#endif

      aWEIGHT1(:,:) = aRSI(:,:)   
!     oWEIGHT1(:,:) = 1.-oWEIGHT(:,:)            ! using the property REGRID(1-RSI)=1-REGRID(RSI)  
      call ab_add(lstr, aWEIGHT1, oWEIGHT1, shape(aWEIGHT1), 'ij')        ! this line should be removed when previous line is uncommented 

      call ab_add(lstr, aRUNPSI, oRUNPSI, shape(aRUNPSI), 'ij',
     &     aWEIGHT1, oWEIGHT1)

      call ab_add(lstr, aSRUNPSI, oSRUNPSI, shape(aSRUNPSI), 
     &     'ij', aWEIGHT1, oWEIGHT1)

      call ab_add(lstr, aERUNPSI, oERUNPSI, shape(aERUNPSI), 
     &     'ij', aWEIGHT1, oWEIGHT1)

#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)
      do N=1,NTM
         if (hasNorthPole(agrid)
     &        .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &        call copy_pole(aTRUNPSI(N,:,aJM),aIM)
      enddo
      call ab_add( lstr, aTRUNPSI, oTRUNPSI, shape(aTRUNPSI), 'lij',
     &     aWEIGHT1, oWEIGHT1)
#endif

      call ab_add( lstr, aRSI, oRSI, shape(aRSI), 'ij')

c*   actual interpolation here
      call bundle_interpolation(lstr,remap_A2O)

c*   polar values are replaced by their longitudinal mean
      if (hasNorthPole(ogrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         call lon_avg(oPREC(:,oJM), oIM)
         call lon_avg(oEPREC(:,oJM), oIM)
         call lon_avg(oRUNPSI(:,oJM), oIM)
         call lon_avg(oSRUNPSI(:,oJM), oIM)
         call lon_avg(oERUNPSI(:,oJM), oIM)
         call lon_avg(oRSI(:,oJM), oIM)
      endif
      if (hasSouthPole(ogrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM)) then
         call lon_avg(oPREC(:,1), oIM)
         call lon_avg(oEPREC(:,1), oIM)
         call lon_avg(oRUNPSI(:,1), oIM)
         call lon_avg(oSRUNPSI(:,1), oIM)
         call lon_avg(oERUNPSI(:,1), oIM)
         call lon_avg(oRSI(:,1), oIM)
      endif

      oEPREC(:,:) = oEPREC(:,:)*OXYP(:,:)
      oSRUNPSI(:,:) = oSRUNPSI(:,:)*OXYP(:,:)
      oERUNPSI(:,:) = oERUNPSI(:,:)*OXYP(:,:)

#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)
      DO N=1,NTM
      if (hasNorthPole(ogrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         call lon_avg(oTRPREC(N,:,oJM), oIM)
         call lon_avg(oTRUNPSI(N,:,oJM), oIM)
      endif
      if (hasSouthPole(ogrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM)) then
         call lon_avg(oTRPREC(N,:,1), oIM)
         call lon_avg(oTRUNPSI(N,:,1), oIM)
      endif
        oTRPREC(N,:,:) = oTRPREC(N,:,:)*OXYP(:,:)
        oTRUNPSI(N,:,:) = oTRUNPSI(N,:,:)*OXYP(:,:)
      END DO
      deallocate(atmp)
#endif

      deallocate(aweight,oweight)
      deallocate(aweight1,oweight1)

      RETURN
      END SUBROUTINE AG2OG_precip

#else
      SUBROUTINE AG2OG_precip

!@sum  INT_AG2OG_precip is for interpolation of precipitation
!!       arrays from atmospheric grid to the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE OCEAN,      only : oIM=>im, oJM=>jm

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid

      USE OCEAN, only : oDXYPO=>DXYPO,oDLATM=>DLATM, OXYP
      Use GEOM,  only : AXYP

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : NTM
#ifdef TRACERS_WATER
      USE FLUXES,  only : aTRPREC=>TRPREC, aTRUNPSI=>TRUNPSI
      USE OFLUXES, only : oTRPREC, oTRUNPSI
#endif
#endif
      USE SEAICE_COM, only : aRSI=>RSI
      USE FLUXES, only : aPREC=>PREC, aEPREC=>EPREC
     *     , aRUNPSI=>RUNPSI, aSRUNPSI=>SRUNPSI, aERUNPSI=>ERUNPSI

      USE OFLUXES, only : oRSI, oPREC, oEPREC
     *     , oRUNPSI, oSRUNPSI, oERUNPSI

      USE INT_AG2OG_MOD, only : INT_AG2OG

      IMPLICIT NONE

      INTEGER N

      REAL*8, allocatable :: aWEIGHT(:,:),atmp(:,:,:)


      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))

      aWEIGHT(:,:) = 1.- aRSI(:,:)  !!  open ocean fraction
      CALL INT_AG2OG(aPREC,oPREC, aWEIGHT)

      CALL INT_AG2OG(aEPREC,oEPREC, aWEIGHT)
      oEPREC(:,:) = oEPREC(:,:)*OXYP(:,:)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aRUNPSI,oRUNPSI, aWEIGHT)

      CALL INT_AG2OG(aSRUNPSI,oSRUNPSI, aWEIGHT)
      oSRUNPSI(:,:) = oSRUNPSI(:,:)*OXYP(:,:)

      CALL INT_AG2OG(aERUNPSI,oERUNPSI, aWEIGHT)
      oERUNPSI(:,:) = oERUNPSI(:,:)*OXYP(:,:)

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aRSI,oRSI, aWEIGHT)

#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)

      aWEIGHT(:,:) = 1.- aRSI(:,:)
      !-- allocate tmp array for trprec
      allocate(atmp(ntm,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))

      DO N=1,NTM
        atmp(N,:,:) = aTRPREC(N,:,:)/AXYP(:,:)
      END DO
      CALL INT_AG2OG(atmp,oTRPREC, aWEIGHT, NTM)
      DO N=1,NTM
        oTRPREC(N,:,:) = oTRPREC(N,:,:)*OXYP(:,:)
      END DO

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aTRUNPSI,oTRUNPSI, aWEIGHT, NTM)
      DO N=1,NTM
        oTRUNPSI(N,:,:) = oTRUNPSI(N,:,:)*OXYP(:,:)
      END DO
 
      deallocate(atmp)
#endif

      deallocate(aweight)

      RETURN
      END SUBROUTINE AG2OG_precip
#endif /* ag2og_precip_bundle */

#ifdef OG2AG_TOC2SST_BUNDLE
      SUBROUTINE OG2AG_TOC2SST
!@sum  OG2AG_oceans: ocean arrays used in the subr. TOC2SST are gathered
!       on the ocean grid, interpolated to the atm. grid, and scattered
!!      on the atm. grid
!@auth Larissa Nazarenko, Denis Gueyffier
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM, oJM=>JM, oLM=>LMO
     *                , oFOCEAN_loc=>FOCEAN_loc
     *                , oDXYPO=>DXYPO, OXYP, oIMAXJ=>IMAXJ
     *                , oCOSI=>COSIC,oSINI=>SINIC
     *                , IVSPO=>IVSP,IVNPO=>IVNP
     *                , sinpo, sinvo
#ifndef CUBED_SPHERE
      Use GEOM,  only : aCOSI=>COSIP,aSINI=>SINIP
#endif
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE,SOUTH,
     &     hasNorthPole, hasSouthPole
      USE OCEANR_DIM, only : ogrid

      USE OCEAN, only : MO, UO,VO, G0M
     *     , S0M, OGEOZ,OGEOZ_SV
#ifdef TRACERS_OCEAN
     *     , TRMO
#endif
      USE AFLUXES, only : aMO, aUO1,aVO1, aG0
     *     , aS0, aOGEOZ,aOGEOZ_SV
#ifdef TRACERS_OCEAN
     *     , aTRAC
#endif
      USE OFLUXES, only : oRSI

#ifdef TRACERS_OCEAN
      USE TRACER_COM, only: NTM
      USE OCN_TRACER_COM, only: conc_from_fw
#endif
#ifdef TRACERS_OceanBiology
!only for TRACERS_OceanBiology and not for seawifs
!/bc we interpolate an internal field
      USE obio_com, only: tot_chlo
      USE FLUXES, only : CHL
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      USE MODEL_COM, only: nstep=>itime
      USE obio_com, only: pCO2
      USE TRACER_COM, only : vol2mass,ntm_gasexch
      USE FLUXES, only : gtracer
#endif
      use ocean, only : remap_O2A
      USE array_bundle
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole

      IMPLICIT NONE
      INTEGER N
      INTEGER IER, I,J,K,L, NT
      INTEGER oJ_0,oJ_1, oI_0,oI_1, oJ_0S,oJ_1S
      REAL*8 :: UNP,VNP,AWT1,AWT2
      REAL*8, ALLOCATABLE :: oG0(:,:,:), oS0(:,:,:)
     *                     , oUO1(:,:), oVO1(:,:), oTRAC(:,:,:)
     *                     , oTOT_CHLO_loc(:,:),opCO2_loc(:,:)
     *                     , oMOtmp(:,:,:),OGEOZtmp(:,:)
     *                     , OGEOZ_SVtmp(:,:)
      REAL*8, allocatable :: aWEIGHT(:,:),oWEIGHT(:,:)
      REAL*8, allocatable :: aWEIGHT1(:,:),oWEIGHT1(:,:)
      REAL*8, allocatable :: aWEIGHT2(:,:),oWEIGHT2(:,:)
      REAL*8, allocatable :: aWEIGHT3(:,:),oWEIGHT3(:,:)
      type (lookup_str) :: lstr

      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &                ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      !-- need to allocate separate array for different weights
      allocate(aweight1(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &                 ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aweight2(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &                 ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight2(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aweight3(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &                 ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight3(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))


      !-- initializing "lstr"
      call ab_init( lstr, 1, oIM,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO,
     &     aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO)

      oI_0 = oGRID%I_STRT
      oI_1 = oGRID%I_STOP
      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP
      oJ_0S = oGRID%j_STRT_SKP
      oJ_1S = oGRID%j_STOP_SKP

      ALLOCATE
     *  (oG0(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2), STAT = IER)
      ALLOCATE
     *  (oS0(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2), STAT = IER)
      ALLOCATE
     *  (oUO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), STAT = IER)
      ALLOCATE
     *  (oVO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), STAT = IER)
#ifdef TRACERS_OCEAN
      ALLOCATE
     *  (oTRAC(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,NTM), STAT = IER)
#endif
      ALLOCATE
     *  (oTOT_CHLO_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
     * , STAT = IER)
      ALLOCATE
     *  (opCO2_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) ,STAT = IER)


      allocate(oMOtmp(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2),
     &     OGEOZtmp(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     OGEOZ_SVtmp(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      call ab_add( lstr, oWEIGHT, aWEIGHT, shape(oWEIGHT),'ij')

      do L=1,2
      oMOtmp(:,:,L) = MO(:,:,L)
      if (hasNorthPole(oGrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(oMOtmp(:,oJM,L),oIM)
      enddo
      call ab_add( lstr, oMOtmp, aMO, shape(oMOtmp), 
     &     'ijk', oWEIGHT, aWEIGHT) 

      OGEOZtmp=OGEOZ
      if (hasNorthPole(oGrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(OGEOZtmp(:,oJM),oIM)
      call ab_add( lstr, OGEOZtmp, aOGEOZ, shape(OGEOZtmp),'ij',
     &     oWEIGHT, aWEIGHT)

      OGEOZ_SVtmp=OGEOZ_SV
      if (hasNorthPole(oGrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(OGEOZ_SVtmp(:,oJM),oIM)
      call ab_add( lstr, OGEOZ_SVtmp, aOGEOZ_SV, shape(OGEOZ_SVtmp),
     &     'ij',oWEIGHT, aWEIGHT)

      oWEIGHT1(:,:) = MO(:,:,1)*oFOCEAN_loc(:,:)
      if (hasNorthPole(oGrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(oWEIGHT1(:,oJM),oIM)
      call ab_add( lstr, oWEIGHT1, aWEIGHT1, shape(oWEIGHT1),'ij')

      oG0(:,:,:) = 0.d0
      oS0(:,:,:) = 0.d0
      DO L = 1,2
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oG0(I,J,L) = G0M(I,J,L)/(MO(I,J,L)*OXYP(I,J))
              oS0(I,J,L) = S0M(I,J,L)/(MO(I,J,L)*OXYP(I,J))
            END IF
          END DO
        END DO
      END DO
      if (hasNorthPole(oGrid)
     &        .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
        do L=1,2
          call copy_pole(oG0(:,oJM,L),oIM)
          call copy_pole(oS0(:,oJM,L),oIM)
        enddo
      endif
      call ab_add( lstr, oG0, aG0, shape(oG0), 'ijk', 
     &     oWEIGHT1, aWEIGHT1) 
      call ab_add( lstr, oS0, aS0, shape(oS0), 'ijk',
     &     oWEIGHT1, aWEIGHT1) 

c Discontinued method for ocean C-grid -> atm A-grid:
c use a variant of INT_OG2AG aware of C-grid staggering.
c May be reinstated later when cubed-sphere INT_OG2AG
c has this capability.
c      oUO1(:,:) = UO(:,:,1)
c      oVO1(:,:) = VO(:,:,1)
c      oWEIGHT(:,:) = 1.d0
c      CALL INT_OG2AG(oUO1,oVO1,aUO1,aVO1, oWEIGHT
c     *             , IVSPO,IVNPO)
c
c ocean C-grid -> atm A-grid method requiring fewer INT_OG2AG variants:
c ocean C -> ocean A followed by ocean A -> atm A via INT_OG2AG
c
      ocean_processors_only: if(oGRID%have_domain) then
        call halo_update(ogrid,vo(:,:,1),from=south)
      endif ocean_processors_only
      do j=oJ_0S,oJ_1S
c area weights that would have been used by HNTRP for ocean C -> ocean A
        awt1 = (sinpo(j)-sinvo(j-1))/(sinvo(j)-sinvo(j-1))
        awt2 = 1.-awt1
        i=1
          oUO1(i,j) = .5*(UO(i,j,1)+UO(oIM,j,1))
          oVO1(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        do i=2,oIM
          oUO1(i,j) = .5*(UO(i,j,1)+UO(i-1,j,1))
          oVO1(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        enddo
      enddo
      if(hasSouthPole(oGRID)) then
        oUO1(:,1) = 0.; oVO1(:,1) = 0.
      endif
      if(hasNorthPole(oGRID)) then ! NP U,V from prognostic polar U,V
        oUO1(:,oJM) = UO(oIM,oJM,1)*oCOSI(:) + UO(IVNPO,oJM,1)*oSINI(:)
! oVO1 currently has no effect when atm is lat-lon
        oVO1(:,oJM) = UO(IVNPO,oJM,1)*oCOSI(:) - UO(oIM,oJM,1)*oSINI(:)
      endif

      call ab_add(lstr, oUO1, aUO1, shape(oUO1),'ij')
      call ab_add(lstr, oVO1, aVO1, shape(oVO1),'ij')

#ifdef TRACERS_OCEAN
C**** surface tracer concentration
      DO NT = 1,NTM
        if (conc_from_fw(nt)) then  ! define conc from fresh water
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTRAC(I,J,NT)=TRMO(I,J,1,NT)/(MO(I,J,1)*OXYP(I,J)
     *             -S0M(I,J,1))
            ELSE
              oTRAC(I,J,NT)=0.
            END IF
          END DO
        END DO
        else  ! define conc from total sea water mass
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTRAC(I,J,NT)=TRMO(I,J,1,NT)/(MO(I,J,1)*OXYP(I,J))
            ELSE
              oTRAC(I,J,NT)=0.
            END IF
          END DO
        END DO
        end if
      if (hasNorthPole(oGrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(oTRAC(:,oJM,NT),oIM)
      END DO

      call ab_add( lstr, oTRAC, aTRAC, shape(oTRAC), 'ijk', 
     &     oWEIGHT1, aWEIGHT1) 

#ifdef TRACERS_OceanBiology
!total ocean chlorophyll. Units are kg,chlorophyll/m3 of seawater
!tot_chlo is defined over all ocean points. Here only use open water
!chorophyll, because that is what is seen by radiation
      oWEIGHT2(:,:) = oFOCEAN_loc(:,:) * (1.d0 - oRSI(:,:))
      call ab_add( lstr, oWEIGHT2, aWEIGHT2, shape(oWEIGHT2),'ij')
      DO J=oJ_0,oJ_1
        DO I=oI_0,oIMAXJ(J)
          IF (oFOCEAN_loc(I,J).gt.0.) THEN
            oTOT_CHLO_loc(I,J) = tot_chlo(I,J)
          ELSE
            oTOT_CHLO_loc(I,J)=0.
          END IF
        END DO
      END DO

      if (hasNorthPole(oGrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(oTOT_CHLO_loc(:,oJM),oIM)

      call ab_add( lstr, oTOT_CHLO_loc, CHL, 
     &     shape(oTOT_CHLO_loc), 'ij', oWEIGHT2, aWEIGHT2) 
#endif

#ifdef TRACERS_GASEXCH_ocean_CO2
!partial CO2 pressure in seawater. Units are uatm.
!defined only over open ocean cells, because this is what is
!involved in gas exchage with the atmosphere.
      oWEIGHT3(:,:) = oFOCEAN_loc(:,:)*(1.d0-oRSI(:,:))
      call ab_add( lstr, oWEIGHT3, aWEIGHT3, shape(oWEIGHT3),'ij')
      DO J=oJ_0,oJ_1
        DO I=oI_0,oIMAXJ(J)
          IF (oFOCEAN_loc(I,J).gt.0.) THEN
            !pco2 is in uatm, convert to kg,CO2/kg,air
            opCO2_loc(I,J) = pCO2(I,J)
     .        * vol2mass(ntm_gasexch)* 1.d-6 ! ppmv (uatm) -> kg,CO2/kg,air
          ELSE
            opCO2_loc(I,J)=0.
          END IF
        END DO
      END DO

      if (hasNorthPole(oGrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(opCO2_loc(:,oJM),oIM)

      DO NT = 1,ntm_gasexch
         call ab_add( lstr, opCO2_loc, aTRAC(:,:,NT), 
     &        shape(opCO2_loc), 'ij', oWEIGHT2, aWEIGHT2) 
      END DO
#endif
#endif


c*    actual interpolation here
      call bundle_interpolation(lstr,remap_O2A)


#ifndef CUBED_SPHERE
      if(hasNorthPole(agrid)) then ! latlon atm needs single polar vector
        UNP = SUM(aUO1(:,aJM)*aCOSI(:))*2/aIM
        VNP = SUM(aUO1(:,aJM)*aSINI(:))*2/aIM
        aUO1(1,aJM) = UNP
        aVO1(1,aJM) = VNP
      endif
#endif

      if(hasNorthPole(agrid) .and.
     &     (aIM .ne. oIM .or. aJM .ne. oJM)) then
         call lon_avg(aOGEOZ(:,aJM),aIM)
         call lon_avg(aOGEOZ_sv(:,aJM),aIM)
        do l=1,2
         call lon_avg(aMO(:,aJM,l),aIM)
         call lon_avg(aG0(:,aJM,l),aIM)
         call lon_avg(aS0(:,aJM,l),aIM)
        enddo
      endif

      if(hasSouthPole(agrid) .and.
     &     (aIM .ne. oIM .or. aJM .ne. oJM)) then
         call lon_avg(aOGEOZ(:,1),aIM)
         call lon_avg(aOGEOZ_sv(:,1),aIM)
        do l=1,2
         call lon_avg(aMO(:,1,l),aIM)
         call lon_avg(aG0(:,1,l),aIM)
         call lon_avg(aS0(:,1,l),aIM)
        enddo
      endif


#ifdef TRACERS_OCEAN

#ifdef TRACERS_GASEXCH_ocean_CO2
      do l=1,ntm_gasexch
         if(hasNorthPole(agrid) .and.
     &        (aIM .ne. oIM .or. aJM .ne. oJM))
     &        call lon_avg(aTRAC(:,aJM,l),aIM)
         if(hasSouthPole(agrid) .and.
     &        (aIM .ne. oIM .or. aJM .ne. oJM)) 
     &        call lon_avg(aTRAC(:,1,l),aIM)
         
         aTRAC(:,:,l) = aTRAC(:,:,l)*vol2mass(l)*1.d-6 ! ppmv (uatm) -> kg,CO2/kg,air
         if (nstep.eq.0) aTRAC(:,:,l) = gtracer(l,1,:,:)
      enddo
#endif
#ifdef TRACERS_OceanBiology
      if(hasNorthPole(agrid) .and.
     &     (aIM .ne. oIM .or. aJM .ne. oJM))
     &     call lon_avg(CHL(:,aJM),aIM)
      if(hasSouthPole(agrid) .and.
     &     (aIM .ne. oIM .or. aJM .ne. oJM))
     &     call lon_avg(CHL(:,1),aIM)
#endif

#ifdef TRACERS_GASEXCH_ocean_CO2
      DO l = 1,ntm_gasexch
         if(hasNorthPole(agrid) .and.
     &        (aIM .ne. oIM .or. aJM .ne. oJM))
     &        call lon_avg(aTRAC(:,aJM,l),aIM)
         if(hasSouthPole(agrid) .and.
     &        (aIM .ne. oIM .or. aJM .ne. oJM)) 
     &        call lon_avg(aTRAC(:,1,l),aIM)
      END DO
#endif

      DEALLOCATE(oTRAC)
#endif


      DEALLOCATE(oG0, oS0, oUO1,oVO1, oTOT_CHLO_loc, opCO2_loc)

      deallocate(oweight,aweight,oweight1,aweight1,
     &     oweight2,aweight2,oweight3,aweight3,
     &     oMOtmp,OGEOZtmp,OGEOZ_SVtmp)

      RETURN
      END SUBROUTINE OG2AG_TOC2SST

#else
      SUBROUTINE OG2AG_TOC2SST
!@sum  OG2AG_oceans: ocean arrays used in the subr. TOC2SST are gathered
!       on the ocean grid, interpolated to the atm. grid, and scattered
!!      on the atm. grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM, oJM=>JM, oLM=>LMO
     *                , oFOCEAN_loc=>FOCEAN_loc
     *                , oDXYPO=>DXYPO, OXYP, oIMAXJ=>IMAXJ
     *                , oCOSI=>COSIC,oSINI=>SINIC
     *                , IVSPO=>IVSP,IVNPO=>IVNP
     *                , sinpo, sinvo
#ifndef CUBED_SPHERE
      Use GEOM,  only : aCOSI=>COSIP,aSINI=>SINIP
#endif
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE,SOUTH,
     &     hasNorthPole, hasSouthPole
      USE OCEANR_DIM, only : ogrid

      USE OCEAN, only : MO, UO,VO, G0M
     *     , S0M, OGEOZ,OGEOZ_SV
#ifdef TRACERS_OCEAN
     *     , TRMO
#endif
      USE AFLUXES, only : aMO, aUO1,aVO1, aG0
     *     , aS0, aOGEOZ,aOGEOZ_SV
#ifdef TRACERS_OCEAN
     *     , aTRAC
#endif
      USE OFLUXES, only : oRSI

#ifdef TRACERS_OCEAN
      USE TRACER_COM, only: NTM
      USE OCN_TRACER_COM, only: conc_from_fw
#endif
#ifdef TRACERS_OceanBiology
!only for TRACERS_OceanBiology and not for seawifs
!/bc we interpolate an internal field
      USE obio_com, only: tot_chlo
      USE FLUXES, only : CHL
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      USE MODEL_COM, only: nstep=>itime
      USE obio_com, only: pCO2
      USE TRACER_COM, only : vol2mass,ntm_gasexch
      USE FLUXES, only : gtracer
#endif

      USE INT_OG2AG_MOD, only : INT_OG2AG

      IMPLICIT NONE

      INTEGER IER, I,J, L, NT
      INTEGER oJ_0,oJ_1, oI_0,oI_1, oJ_0S,oJ_1S
      REAL*8, allocatable :: oWEIGHT(:,:)
      REAL*8 :: UNP,VNP,AWT1,AWT2
      REAL*8, ALLOCATABLE :: oG0(:,:,:), oS0(:,:,:)
     *                     , oUO1(:,:), oVO1(:,:), oTRAC(:,:,:)
     *                     , oTOT_CHLO_loc(:,:),opCO2_loc(:,:)
     *                     ,atest(:,:)

      oI_0 = oGRID%I_STRT
      oI_1 = oGRID%I_STOP
      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP
      oJ_0S = oGRID%j_STRT_SKP
      oJ_1S = oGRID%j_STOP_SKP

      ALLOCATE
     *  (oG0(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2), STAT = IER)
      ALLOCATE
     *  (oS0(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2), STAT = IER)
      ALLOCATE
     *  (oUO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), STAT = IER)
      ALLOCATE
     *  (oVO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), STAT = IER)
#ifdef TRACERS_OCEAN
      ALLOCATE
     *  (oTRAC(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,NTM), STAT = IER)
#endif
      ALLOCATE
     *  (oTOT_CHLO_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
     * , STAT = IER)
      ALLOCATE
     *  (opCO2_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) ,STAT = IER)

      allocate(oWEIGHT(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(MO,aMO,oWEIGHT,oLM,2,.FALSE.)

      oG0(:,:,:) = 0.d0
      oWEIGHT(:,:) = MO(:,:,1)*oFOCEAN_loc(:,:)
      DO L = 1,2
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oG0(I,J,L) = G0M(I,J,L)/(MO(I,J,L)*OXYP(I,J))
            END IF
          END DO
        END DO
      END DO
      CALL INT_OG2AG(oG0,aG0, oWEIGHT, 2,2,.TRUE.)

      oS0(:,:,:) = 0.d0
      DO L = 1,2
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oS0(I,J,L) = S0M(I,J,L)/(MO(I,J,L)*OXYP(I,J))
            END IF
          END DO
        END DO
      END DO
      CALL INT_OG2AG(oS0,aS0, oWEIGHT, 2,2,.TRUE.)

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(OGEOZ,aOGEOZ, oWEIGHT, .FALSE.)

      CALL INT_OG2AG(OGEOZ_SV,aOGEOZ_SV, oWEIGHT, .FALSE.)

c Discontinued method for ocean C-grid -> atm A-grid:
c use a variant of INT_OG2AG aware of C-grid staggering.
c May be reinstated later when cubed-sphere INT_OG2AG
c has this capability.
c      oUO1(:,:) = UO(:,:,1)
c      oVO1(:,:) = VO(:,:,1)
c      oWEIGHT(:,:) = 1.d0
c      CALL INT_OG2AG(oUO1,oVO1,aUO1,aVO1, oWEIGHT
c     *             , IVSPO,IVNPO)
c
c ocean C-grid -> atm A-grid method requiring fewer INT_OG2AG variants:
c ocean C -> ocean A followed by ocean A -> atm A via INT_OG2AG
c
      call halo_update(ogrid,vo(:,:,1),from=south)
      do j=oJ_0S,oJ_1S
c area weights that would have been used by HNTRP for ocean C -> ocean A
        awt1 = (sinpo(j)-sinvo(j-1))/(sinvo(j)-sinvo(j-1))
        awt2 = 1.-awt1
        i=1
          oUO1(i,j) = .5*(UO(i,j,1)+UO(oIM,j,1))
          oVO1(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        do i=2,oIM
          oUO1(i,j) = .5*(UO(i,j,1)+UO(i-1,j,1))
          oVO1(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        enddo
      enddo
      if(hasSouthPole(oGRID)) then
        oUO1(:,1) = 0.; oVO1(:,1) = 0.
      endif
      if(hasNorthPole(oGRID)) then ! NP U,V from prognostic polar U,V
        oUO1(:,oJM) = UO(oIM,oJM,1)*oCOSI(:) + UO(IVNPO,oJM,1)*oSINI(:)
! oVO1 currently has no effect when atm is lat-lon
        oVO1(:,oJM) = UO(IVNPO,oJM,1)*oCOSI(:) - UO(oIM,oJM,1)*oSINI(:)
      endif
      oWEIGHT(:,:) = 1.d0
      CALL INT_OG2AG(oUO1, aUO1, oWEIGHT, .FALSE., AvgPole=.FALSE.)
      CALL INT_OG2AG(oVO1, aVO1, oWEIGHT, .FALSE., AvgPole=.FALSE.)
#ifndef CUBED_SPHERE
      if(hasNorthPole(aGRID)) then ! latlon atm needs single polar vector
        UNP = SUM(aUO1(:,aJM)*aCOSI(:))*2/aIM
        VNP = SUM(aUO1(:,aJM)*aSINI(:))*2/aIM
        aUO1(1,aJM) = UNP
        aVO1(1,aJM) = VNP
      endif
#endif

#ifdef TRACERS_OCEAN
C**** surface tracer concentration
      oWEIGHT(:,:) = MO(:,:,1)*oFOCEAN_loc(:,:)
      DO NT = 1,NTM
        if (conc_from_fw(nt)) then  ! define conc from fresh water
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTRAC(I,J,NT)=TRMO(I,J,1,NT)/(MO(I,J,1)*OXYP(I,J)
     *             -S0M(I,J,1))
            ELSE
              oTRAC(I,J,NT)=0.
            END IF
          END DO
        END DO
        else  ! define conc from total sea water mass
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTRAC(I,J,NT)=TRMO(I,J,1,NT)/(MO(I,J,1)*OXYP(I,J))
            ELSE
              oTRAC(I,J,NT)=0.
            END IF
          END DO
        END DO
        end if
      END DO
      CALL INT_OG2AG(oTRAC,aTRAC, oWEIGHT, NTM,NTM,.TRUE.)

#ifdef TRACERS_OceanBiology
!total ocean chlorophyll. Units are kg,chlorophyll/m3 of seawater
!tot_chlo is defined over all ocean points. Here only use open water
!chorophyll, because that is what is seen by radiation
      oWEIGHT(:,:) = oFOCEAN_loc(:,:) * (1.d0 - oRSI(:,:))
      DO J=oJ_0,oJ_1
        DO I=oI_0,oIMAXJ(J)
          IF (oFOCEAN_loc(I,J).gt.0.) THEN
            oTOT_CHLO_loc(I,J) = tot_chlo(I,J)
          ELSE
            oTOT_CHLO_loc(I,J)=0.
          END IF
        END DO
      END DO
      CALL INT_OG2AG(oTOT_CHLO_loc,CHL, oWEIGHT, .FALSE.)
#endif

#ifdef TRACERS_GASEXCH_ocean_CO2
!partial CO2 pressure in seawater. Units are uatm.
!defined only over open ocean cells, because this is what is
!involved in gas exchage with the atmosphere.
      oWEIGHT(:,:) = oFOCEAN_loc(:,:)*(1.d0-oRSI(:,:))
      DO J=oJ_0,oJ_1
        DO I=oI_0,oIMAXJ(J)
          IF (oFOCEAN_loc(I,J).gt.0.) THEN
            !pco2 is in uatm, convert to kg,CO2/kg,air
            opCO2_loc(I,J) = pCO2(I,J)
     .          * vol2mass(ntm_gasexch)* 1.d-6 ! ppmv (uatm) -> kg,CO2/kg,air
          ELSE
            opCO2_loc(I,J)=0.
          END IF
        END DO
      END DO
      DO NT = 1,ntm_gasexch
         CALL INT_OG2AG(opCO2_loc,aTRAC(:,:,NT), oWEIGHT, .FALSE.)

         !gtracer is first set in TRACER_DRV, then atrac is interpolated
         !here from pco2 in the ocean and later OCNDYN sets gtracer=atrac
         !Therefore in timesetep 0 (cold start) pco2 has not yet been defined
         !and atrac has to be hard coded in here, so that we do not have
         !urealistic tracer flux at the air-sea interface during step0.

         if (nstep.eq.0) aTRAC(:,:,NT) = gtracer(NT,1,:,:)

      END DO
#endif
#endif

      DEALLOCATE(oG0, oS0, oUO1,oVO1, oTOT_CHLO_loc, opCO2_loc)
#ifdef TRACERS_OCEAN
      DEALLOCATE(oTRAC)
#endif

      deallocate(oweight)

      RETURN
      END SUBROUTINE OG2AG_TOC2SST

#endif /* OG2AG_TOC2SST_BUNDLE */

#ifdef AG2OG_OCEANS_BUNDLE

      SUBROUTINE AG2OG_oceans
!@sum  AG2OG_oceans: all atmospheric arrays used in the subr. OCEANS are gathered
!       on the atmospheric grid, interpolated to the ocean grid, and scattered
!!      on the ocean grid
!@auth Larissa Nazarenko, Denis Gueyffier
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM, oJM=>JM
     &                , oDXYPO=>DXYPO, OXYP
     &                , IVSPO=>IVSP,IVNPO=>IVNP
     &                , oSINI=>SINIC, oCOSI=>COSIC 
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
#if defined (TRACERS_OceanBiology) && !defined (TRACERS_GASEXCH_ocean)
      USE OCN_TRACER_COM, only: NTM
#else
      USE TRACER_COM, only: NTM
#endif
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole
      USE OCEANR_DIM, only : ogrid

      USE SEAICE_COM, only : aRSI=>RSI
      USE FLUXES, only : aSOLAR=>SOLAR, aE0=>E0, aEVAPOR=>EVAPOR
     *     , aRUNOSI=>RUNOSI,aERUNOSI=>ERUNOSI,aSRUNOSI=>SRUNOSI
     *     , aFLOWO=>FLOWO,aEFLOWO=>EFLOWO, aAPRESS=>APRESS
     *     , aMELTI=>MELTI,aEMELTI=>EMELTI,aSMELTI=>SMELTI
     *     , aDMUA=>DMUA, aDMVA=>DMVA
     *     , aGMELT=>GMELT, aEGMELT=>EGMELT
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , aTRFLOWO=>TRFLOWO, aTREVAPOR=>TREVAPOR
     *     , aTRUNOSI=>TRUNOSI, aTRMELTI=>TRMELTI
     *     , aTRGMELT=>TRGMELT
#ifdef TRACERS_DRYDEP
     *     , aTRDRYDEP=>TRDRYDEP
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE FLUXES, only : aTRGASEX=>TRGASEX
      USE DOMAIN_DECOMP_1D, only : pack_data
#endif
#ifdef OBIO_RAD_coupling
      USE RAD_COM, only : avisdir    => FSRDIR
     *                   ,asrvissurf => SRVISSURF
     *                   ,avisdif    => FSRDIF
     *                   ,anirdir    => DIRNIR
     *                   ,anirdif    => DIFNIR
#endif
#ifdef TRACERS_OceanBiology
      USE RAD_COM, only : aCOSZ1=>COSZ1
      USE PBLCOM, only : awind=>wsavg
#endif

      USE OFLUXES, only : oRSI, oSOLAR, oE0, oEVAPOR
     *     , oRUNOSI, oERUNOSI, oSRUNOSI
     *     , oFLOWO, oEFLOWO, oAPRESS
     *     , oMELTI, oEMELTI, oSMELTI
     *     , oDMUA,oDMVA
     *     , oGMELT, oEGMELT
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , oTRFLOWO, oTREVAPOR
     *     , oTRUNOSI, oTRMELTI
     *     , oTRGMELT
#ifdef TRACERS_DRYDEP
     *     , oTRDRYDEP
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE MODEL_COM, only: nstep=>itime
      USE OFLUXES, only : oTRGASEX
#endif
#ifdef OBIO_RAD_coupling
      USE obio_forc, only : ovisdir,ovisdif,onirdir,onirdif
#endif
#ifdef TRACERS_OceanBiology
      USE obio_forc, only : osolz, owind
#endif

      Use GEOM,  only : AXYP,aIMAXJ=>IMAXJ
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole
#ifndef CUBED_SPHERE
       Use GEOM,  only : aCOSI=>COSIP,aSINI=>SINIP
#endif
      USE MODEL_COM, only : aFOCEAN_loc=>FOCEAN

      use ocean, only : remap_A2O
      USE array_bundle

      IMPLICIT NONE
      REAL*8, allocatable :: aWEIGHT(:,:),oWEIGHT(:,:)
      REAL*8, allocatable :: aWEIGHT1(:,:),oWEIGHT1(:,:)
      REAL*8, allocatable :: aWEIGHT2(:,:),oWEIGHT2(:,:)
      REAL*8, allocatable :: aWEIGHT3(:,:),oWEIGHT3(:,:)
      REAL*8, allocatable :: aWEIGHT4(:,:),oWEIGHT4(:,:)
      REAL*8, allocatable :: aWEIGHT5(:,:),oWEIGHT5(:,:)
      REAL*8, allocatable :: aWEIGHT6(:,:),oWEIGHT6(:,:)
      REAL*8, allocatable :: aE0tmp(:,:),oE0tmp(:,:),
     &     aEVAPORtmp(:,:),oEVAPORtmp(:,:),
     &     aSOLAR1tmp(:,:),oSOLAR1tmp(:,:),
     &     aSOLAR3tmp(:,:),oSOLAR3tmp(:,:),
     &     aDMUA1tmp(:,:),oDMUA1tmp(:,:),
     &     aDMVA1tmp(:,:),oDMVA1tmp(:,:),
     &     atmp01(:,:),atmp02(:,:),
     &     atmp03(:,:),atmp04(:,:),
     &     atmp05(:,:),atmp06(:,:),
     &     atmp07(:,:),
     &     atmp08(:,:,:),
     &     atmp09(:,:,:),atmp10(:,:,:),
     &     atmp(:,:,:),otmp(:,:,:),
     &     atmp2(:,:,:),otmp2(:,:,:),
     &     atmp3(:,:,:),otmp3(:,:,:),
     &     atmp4(:,:),otmp4(:,:)
      REAL*8 :: aDMUAnp,aDMVAnp,aDMUAsp,aDMVAsp
      INTEGER, PARAMETER :: NSTYPE=4
      INTEGER I,J,N
      INTEGER aJ_0,aJ_1, aI_0,aI_1
      INTEGER oJ_0,oJ_1
      REAL*8,
     *     DIMENSION(aGRID%I_STRT:aGRID%I_STOP,
     *     aGRID%J_STRT:aGRID%J_STOP)::
     *     aFact
      REAL*8 :: oDMUA1sp,oDMVA1sp,oDMUA1np,oDMVA1np
      type (lookup_str) :: lstr


      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aweight1(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight1(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aweight2(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight2(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      allocate(aweight3(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight3(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aweight4(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight4(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aweight5(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight5(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aweight6(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight6(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      !-- initializing "lstr"
      call ab_init( lstr, aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO,
     &     oGRID%I_STRT_HALO, oGRID%I_STOP_HALO,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO )

      aJ_0 = aGRID%j_STRT
      aJ_1 = aGRID%j_STOP
      aI_0 = aGRID%I_STRT
      aI_1 = aGRID%I_STOP
      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP

      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(aRSI(:,aJM),aIM)
      call ab_add( lstr, aRSI, oRSI, shape(aRSI),'ij')

      allocate(atmp01(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &     ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     atmp02(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &     ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     atmp03(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &     ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     atmp04(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &     ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     atmp05(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &     ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     atmp06(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &     ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     atmp07(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &     ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)  )

      atmp01=0.d0
      atmp02=0.d0
      atmp03=0.d0
      atmp04=0.d0
      atmp05=0.d0
      atmp06=0.d0
      atmp07=0.d0

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
            atmp01(I,J) = aFLOWO(I,J)*aFact(I,J)
            atmp02(I,J) = aEFLOWO(I,J)*aFact(I,J)
            atmp03(I,J) = aMELTI(I,J)*aFact(I,J)
            atmp04(I,J) = aEMELTI(I,J)*aFact(I,J)
            atmp05(I,J) = aSMELTI(I,J)*aFact(I,J)
            atmp06(I,J) = aGMELT(I,J)*aFact(I,J)
          END IF
        END DO
      END DO

      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
        call copy_pole(atmp01(:,aJM),aIM)
        call copy_pole(atmp02(:,aJM),aIM)
        call copy_pole(atmp03(:,aJM),aIM)
        call copy_pole(atmp04(:,aJM),aIM)
        call copy_pole(atmp05(:,aJM),aIM)
        call copy_pole(atmp06(:,aJM),aIM)
      endif

      call ab_add( lstr, atmp01, oFLOWO, shape(atmp01), 'ij')
      call ab_add( lstr, atmp02, oEFLOWO, shape(atmp02), 'ij')
      call ab_add( lstr, atmp03, oMELTI, shape(atmp03), 'ij')
      call ab_add( lstr, atmp04, oEMELTI, shape(atmp04), 'ij')
      call ab_add( lstr, atmp05, oSMELTI, shape(atmp05), 'ij')
      call ab_add(lstr, atmp06, oGMELT, shape(atmp06), 'ij')

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            atmp07(I,J) = aEGMELT(I,J)*aFact(I,J)
            aFact(I,J) = 1.d0/AXYP(I,J) ! not used?
          END IF
        END DO
      END DO
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(atmp07(:,aJM),aIM)
      call ab_add( lstr, atmp07, oEGMELT, shape(atmp07), 'ij')

      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(aAPRESS(:,aJM),aIM)
      call ab_add( lstr, aAPRESS, oAPRESS, shape(aAPRESS), 'ij')

      aWEIGHT(:,:) = aRSI(:,:)
      call ab_add( lstr, aWEIGHT, oWEIGHT, shape(aWEIGHT), 'ij')

      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
        call copy_pole(aRUNOSI(:,aJM),aIM)
        call copy_pole(aERUNOSI(:,aJM),aIM)
        call copy_pole(aSRUNOSI(:,aJM),aIM)
      endif
      call ab_add( lstr, aRUNOSI, oRUNOSI, shape(aRUNOSI),
     &     'ij', aWEIGHT, oWEIGHT)
      call ab_add( lstr, aERUNOSI, oERUNOSI, shape(aERUNOSI), 
     &     'ij', aWEIGHT, oWEIGHT)
      call ab_add( lstr, aSRUNOSI, oSRUNOSI, shape(aSRUNOSI), 
     &     'ij', aWEIGHT, oWEIGHT)

      aWEIGHT1(:,:) = 1.d0 - aRSI(:,:)
      call ab_add( lstr, aWEIGHT1, oWEIGHT1, shape(aWEIGHT1), 'ij')

!     allocate temporary arrays 
      allocate(aE0tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oE0tmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aEVAPORtmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oEVAPORtmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aSOLAR1tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oSOLAR1tmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aSOLAR3tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oSOLAR3tmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aDMUA1tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oDMUA1tmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aDMVA1tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oDMVA1tmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

!     copy E0(:,:,1) in temporary 2d array on atm grid, 
!     and return temporary 2d array on ocean grid
      aE0tmp(:,:)=aE0(:,:,1)
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(aE0tmp(:,aJM),aIM)
      call ab_add(lstr, aE0tmp, oE0tmp, shape(aE0tmp), 'ij',
     &     aWEIGHT1, oWEIGHT1)

!     copy EVAPOR(:,:,1) in temporary 2d array on atm grid, 
!     and return temporary 2d array on ocean grid
      aEVAPORtmp(:,:)=aEVAPOR(:,:,1)
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(aEVAPORtmp(:,aJM),aIM)
      call ab_add(lstr, aEVAPORtmp, oEVAPORtmp, shape(aEVAPORtmp),
     &     'ij', aWEIGHT1, oWEIGHT1)

!     copy SOLAR(1,:,:) in temporary 2d array on atm grid, 
!     and return temporary 2d array on ocean grid
      aSOLAR1tmp(:,:)=aSOLAR(1,:,:)
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(aSOLAR1tmp(:,aJM),aIM)
      call ab_add(lstr, aSOLAR1tmp, oSOLAR1tmp, shape(aSOLAR1tmp),
     &     'ij', aWEIGHT1, oWEIGHT1)

      aWEIGHT2(:,:) = aRSI(:,:)
      call ab_add( lstr, aWEIGHT2, oWEIGHT2, shape(aWEIGHT2), 'ij')

!     copy SOLAR(1,:,:) in temporary 2d array on atm grid, 
!     and return temporary 2d array on ocean grid
      aSOLAR3tmp(:,:)=aSOLAR(3,:,:)
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(aSOLAR3tmp(:,aJM),aIM)
      call ab_add( lstr, aSOLAR3tmp, oSOLAR3tmp, shape(aSOLAR3tmp), 
     &     'ij', aWEIGHT2, oWEIGHT2)

#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      allocate(atmp08(NTM,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &     ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     atmp09(NTM,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &     ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     atmp10(NTM,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &     ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)      )
      atmp08=0.d0
      atmp09=0.d0
      atmp10=0.d0

      aWEIGHT(:,:) = 1.d0
      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
              atmp08(N,I,J) = aTRFLOWO(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(atmp08(N,:,aJM),aIM)
      END DO
      call ab_add( lstr, atmp08,oTRFLOWO,shape(atmp08),'lij')

      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              atmp09(N,I,J) = aTRMELTI(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(atmp09(N,:,aJM),aIM)
      END DO
      call ab_add( lstr, atmp09,oTRMELTI,shape(atmp09),'lij')


      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
              atmp10(N,I,J) = aTRGMELT(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
        if ( (hasNorthPole(agrid))
     &       .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &       call copy_pole(atmp10(N,:,aJM),aIM)
      END DO
      call ab_add( lstr, atmp10,oTRGMELT,shape(atmp10),'lij')

      aWEIGHT3(:,:) = aRSI(:,:)
      call ab_add(lstr, aWEIGHT3, oWEIGHT3, shape(aWEIGHT3), 'ij')

      DO N=1,NTM
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(aTRUNOSI(N,:,aJM),aIM)
      END DO
      call ab_add( lstr, aTRUNOSI,oTRUNOSI,shape(aTRUNOSI),'lij',
     &     aWEIGHT3,oWEIGHT3)

      allocate(otmp(NTM,oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )
      allocate(atmp(NTM,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )

      aWEIGHT4(:,:) = 1.d0 - aRSI(:,:)
      call ab_add(lstr, aWEIGHT4, oWEIGHT4, shape(aWEIGHT4), 'ij')
      do N=1,NTM
         atmp(N,:,:)=aTREVAPOR(N,1,:,:)
         if ( (hasNorthPole(agrid))
     &        .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &        call copy_pole(atmp(N,:,aJM),aIM)
      enddo
      call ab_add( lstr, atmp,otmp,shape(atmp),'lij',
     &     aWEIGHT4,oWEIGHT4)

#ifdef TRACERS_DRYDEP
      allocate(otmp2(NTM,oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )
      allocate(atmp2(NTM,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )

      do N=1,NTM
      atmp2(N,:,:)=aTRDRYDEP(N,1,:,:)
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(atmp2(N,:,aJM),aIM)
      enddo
      call ab_add( lstr,atmp2,otmp2,shape(atmp2),'lij')
#endif

#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      allocate(otmp3(NTM,oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )
      allocate(atmp3(NTM,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )

      do N=1,ntm_gasexch
       atmp3(N,:,:)=aTRGASEX(N,1,:,:)
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(atmp3(N,:,aJM),aIM)
      enddo
      call ab_add( lstr,atmp3,otmp3,shape(atmp3),'lij')
#endif

#ifdef OBIO_RAD_coupling
      allocate(otmp4(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )
      allocate(atmp4(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )

      aWEIGHT5(:,:) = aFOCEAN_loc(:,:)
      call ab_add(lstr, aWEIGHT5, oWEIGHT5, shape(aWEIGHT5), 'ij')

      atmp4(:,:)=aVISDIR(:,:)*aSRVISSURF(:,:)
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(atmp4(:,aJM),aIM)

      call ab_add( lstr,atmp4,oVISDIR,shape(atmp4),'ij',
     &     aWEIGHT5,oWEIGHT5)

      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(aVISDIF(:,aJM),aIM)

      call ab_add( lstr,aVISDIF,oVISDIF,shape(aVISDIF),'ij',
     &     aWEIGHT5,oWEIGHT5)

      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(aNIRDIR(:,aJM),aIM)

      call ab_add( lstr,aNIRDIR,oNIRDIR,shape(aNIRDIR),'ij',
     &     aWEIGHT5,oWEIGHT5)

      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) )
     &     call copy_pole(aNIRDIF(:,aJM),aIM)

      call ab_add( lstr,aNIRDIF,oNIRDIF,shape(aNIRDIF),'ij',
     &     aWEIGHT5,oWEIGHT5)
#endif

#ifdef TRACERS_OceanBiology
      aWEIGHT6(:,:) = aFOCEAN_loc(:,:)

      call ab_add(lstr, aWEIGHT6, oWEIGHT6, shape(aWEIGHT6), 'ij')

      call ab_add( lstr,aCOSZ1,oSOLZ,shape(aCOSZ1),'ij',
     &     aWEIGHT6,oWEIGHT6)

      call ab_add( lstr,aWIND,oWIND,shape(aWIND),'ij',
     &     aWEIGHT6,oWEIGHT6)
#endif

#ifndef CUBED_SPHERE
      if (oIM .eq. aIM .and. oJM .eq. aJM) then
         aDMUA1tmp(:,:) = aDMUA(:,:,1)
         aDMVA1tmp(:,:) = aDMVA(:,:,1)
      else
         
         if (hasNorthPole(agrid)) then
            aDMUAnp = aDMUA(1,aJM,1)
            aDMVAnp = aDMVA(1,aJM,1)
            aDMUA1tmp(:,aJM) = aDMUAnp*aCOSI(:) + aDMVAnp*aSINI(:)
            aDMVA1tmp(:,aJM) = aDMVAnp*aCOSI(:) - aDMUAnp*aSINI(:)
         endif
         if (hasSouthPole(agrid)) then
            aDMUAsp = aDMUA(1,1,1)
            aDMVAsp = aDMVA(1,1,1)
            aDMUA1tmp(:,1) = aDMUAsp*aCOSI(:) - aDMVAsp*aSINI(:)
            aDMVA1tmp(:,1) = aDMVAsp*aCOSI(:) + aDMUAsp*aSINI(:)
         endif

         do j=max(2,aGRID%J_STRT_HALO),min(aJM-1,aGRID%J_STOP_HALO) ! exclude poles
            do i=1,aIMAXJ(J)
               if (aFOCEAN_loc(i,j).gt.0.) then
                  aDMUA1tmp(i,j) = aDMUA(i,j,1) 
               endif
            enddo
         enddo
         
         do j=max(2,aGRID%J_STRT_HALO),min(aJM-1,aGRID%J_STOP_HALO) ! exclude poles
            do i=1,aIMAXJ(J)
               if (aFOCEAN_loc(i,j).gt.0.) then
                  aDMVA1tmp(i,j) = aDMVA(i,j,1) 
               endif
            enddo
         enddo
         
      endif
#else
      do j=aGRID%J_STRT_HALO,aGRID%J_STOP_HALO
         do i=aGRID%I_STRT_HALO,aGRID%I_STOP_HALO
            if (aFOCEAN_loc(i,j).gt.0.) then
               aDMUA1tmp(i,j) = aDMUA(i,j,1) 
            endif
         enddo
      enddo 
      
      do j=aGRID%J_STRT_HALO,aGRID%J_STOP_HALO
         do i=aGRID%I_STRT_HALO,aGRID%I_STOP_HALO
            if (aFOCEAN_loc(i,j).gt.0.) then
               aDMVA1tmp(i,j) = aDMVA(i,j,1) 
            endif
         enddo
      enddo
#endif /* not CUBED_SPHERE */

      call ab_add(lstr, aDMUA1tmp, oDMUA1tmp, shape(aDMUA1tmp), 'ij')
      call ab_add(lstr, aDMVA1tmp, oDMVA1tmp, shape(aDMVA1tmp), 'ij')


c*   actual interpolation here
      call bundle_interpolation(lstr,remap_A2O)
c*


c*   polar values are replaced by their longitudinal mean
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         call lon_avg( oRSI(:,oJM), oIM)
         call lon_avg( oFLOWO(:,oJM), oIM)
         call lon_avg( oEFLOWO(:,oJM), oIM)
         call lon_avg( oMELTI(:,oJM), oIM)
         call lon_avg( oEMELTI(:,oJM), oIM)
         call lon_avg( oSMELTI(:,oJM), oIM)
         call lon_avg( oGMELT(:,oJM), oIM)
         call lon_avg( oEGMELT(:,oJM), oIM)
         call lon_avg( oAPRESS(:,oJM), oIM)
         call lon_avg( oRUNOSI(:,oJM), oIM)
         call lon_avg( oERUNOSI(:,oJM), oIM)
         call lon_avg( oSRUNOSI(:,oJM), oIM)
         call lon_avg( oE0tmp(:,oJM), oIM)
         call lon_avg( oEVAPORtmp(:,oJM), oIM)
         call lon_avg( oSOLAR1tmp(:,oJM), oIM)
         call lon_avg( oSOLAR3tmp(:,oJM), oIM)
      endif

      if ( (hasSouthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         call lon_avg( oRSI(:,1), oIM)
         call lon_avg( oFLOWO(:,1), oIM)
         call lon_avg( oEFLOWO(:,1), oIM)
         call lon_avg( oMELTI(:,1), oIM)
         call lon_avg( oEMELTI(:,1), oIM)
         call lon_avg( oSMELTI(:,1), oIM)
         call lon_avg( oGMELT(:,1), oIM)
         call lon_avg( oEGMELT(:,1), oIM)
         call lon_avg( oAPRESS(:,1), oIM)
         call lon_avg( oRUNOSI(:,1), oIM)
         call lon_avg( oERUNOSI(:,1), oIM)
         call lon_avg( oSRUNOSI(:,1), oIM)
         call lon_avg( oE0tmp(:,1), oIM)
         call lon_avg( oEVAPORtmp(:,1), oIM)
         call lon_avg( oSOLAR1tmp(:,1), oIM)
         call lon_avg( oSOLAR3tmp(:,1), oIM)
      endif
c*

      oEGMELT(:,:) = oEGMELT(:,:)*OXYP(:,:)
      oGMELT(:,:) =  oGMELT(:,:)*OXYP(:,:)

      oE0(:,:,1)=oE0tmp(:,:)  

      oEVAPOR(:,:,1)=oEVAPORtmp(:,:)  

      oSOLAR(1,:,:)=oSOLAR1tmp(:,:)  

      oSOLAR(3,:,:)=oSOLAR3tmp(:,:) 

#ifdef CUBED_SPHERE
      if (hasSouthPole(ogrid)) then
         oDMUA1sp =  SUM(oDMUA1tmp(:,  1)*oCOSI(:))*2/oIM
         oDMVA1sp = -SUM(oDMUA1tmp(:,  1)*oSINI(:))*2/oIM
         oDMUA1tmp(1,  1) = oDMUA1sp
         oDMVA1tmp(1,  1) = oDMVA1sp
      endif
      if (hasNorthPole(oGrid)) then
         oDMUA1np =  SUM(oDMUA1tmp(:,oJM)*oCOSI(:))*2/oIM
         oDMVA1np =  SUM(oDMUA1tmp(:,oJM)*oSINI(:))*2/oIM
         oDMUA1tmp(1,oJM) = oDMUA1np
         oDMVA1tmp(1,oJM) = oDMVA1np
      endif
#else
      if (oIM .ne. aIM .or. oJM .ne. aJM) then
         
         if (hasSouthPole(ogrid)) then
            oDMUA1sp =  SUM(oDMUA1tmp(:,  1)*oCOSI(:))*2/oIM
            oDMVA1sp = -SUM(oDMUA1tmp(:,  1)*oSINI(:))*2/oIM
            oDMUA1tmp(1,  1) = oDMUA1sp
            oDMVA1tmp(1,  1) = oDMVA1sp
         endif
         if (hasNorthPole(oGrid) then
            oDMUA1np =  SUM(oDMUA1tmp(:,oJM)*oCOSI(:))*2/oIM
            oDMVA1np =  SUM(oDMUA1tmp(:,oJM)*oSINI(:))*2/oIM
            oDMUA1tmp(1,oJM) = oDMUA1np
            oDMVA1tmp(1,oJM) = oDMVA1np
         endif

      endif
#endif

#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER

      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         do N=1,NTM
            call lon_avg( oTRFLOWO(N,:,oJM), oIM)
            call lon_avg( oTRMELTI(N,:,oJM), oIM)
            call lon_avg( oTRGMELT(N,:,oJM), oIM)
            call lon_avg( oTRUNOSI(N,:,oJM), oIM)
            call lon_avg( otmp(N,:,oJM), oIM)
         enddo
      endif

      if ( (hasSouthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         do N=1,NTM
            call lon_avg( oTRFLOWO(N,:,1), oIM)
            call lon_avg( oTRMELTI(N,:,1), oIM)
            call lon_avg( oTRGMELT(N,:,1), oIM)
            call lon_avg( oTRUNOSI(N,:,1), oIM)
            call lon_avg( otmp(N,:,1), oIM)
         enddo
      endif

      DO N=1,NTM
        oTRGMELT(N,:,:) = oTRGMELT(N,:,:)*OXYP(:,:)
      END DO

      do N=1,NTM
         oTREVAPOR(N,1,:,:)=otmp(N,:,:)
      enddo
      deallocate(atmp,otmp,atmp08,atmp09,atmp10)

#ifdef TRACERS_DRYDEP
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         do N=1,NTM
            call lon_avg( otmp2(N,:,oJM), oIM)
         enddo
      endif
      if ( (hasSouthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         do N=1,NTM
            call lon_avg( otmp2(N,:,1), oIM)
         enddo
      endif

      do N=1,NTM
         oTRDRYDEP(N,1,:,:)=otmp2(N,:,:)
      enddo
      deallocate(atmp2,otmp2)
#endif
#ifdef TRACERS_GASEXCH_ocean
      if ( (hasNorthPole(agrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         do N=1,ntm_gasexch
            call lon_avg( otmp3(N,:,oJM), oIM)
         enddo
      endif
      if ( (hasSouthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         do N=1,ntm_gasexch
            call lon_avg( otmp3(N,:,1), oIM)
         enddo
      endif

      do N=1,ntm_gasexch
         oTRGASEX(N,1,:,:)=otmp3(N,:,:)
      enddo
      deallocate(atmp3,otmp3)
#endif
#ifdef OBIO_RAD_coupling
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         call lon_avg( oVISDIR(:,oJM), oIM)
         call lon_avg( oVISDIF(:,oJM), oIM)
         call lon_avg( oNIRDIR(:,oJM), oIM)
         call lon_avg( oNIRDIF(:,oJM), oIM)
      endif
      if ( (hasSouthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         call lon_avg( oVISDIR(:,1), oIM)
         call lon_avg( oVISDIF(:,1), oIM)
         call lon_avg( oNIRDIR(:,1), oIM)
         call lon_avg( oNIRDIF(:,1), oIM)
      endif
      deallocate(atmp4,otmp4)
#endif
#ifdef TRACERS_OceanBiology
      if ( (hasNorthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         call lon_avg( oSOLZ(:,oJM), oIM)
         call lon_avg( oWIND(:,oJM), oIM)
      endif
      if ( (hasSouthPole(agrid))
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         call lon_avg( oSOLZ(:,1), oIM)
         call lon_avg( oWIND(:,1), oIM)
      endif
#endif
#endif
#endif

      oDMUA(:,:,1)=oDMUA1tmp(:,:)
      oDMVA(:,:,1)=oDMVA1tmp(:,:)
      
      deallocate(aweight,oweight,
     &     aweight1,oweight1,
     &     aweight2,oweight2,
     &     aweight3,oweight3,
     &     aweight4,oweight4,
     &     aweight5,oweight5,
     &     aweight6,oweight6)


      deallocate(atmp01,atmp02,
     &     atmp03,atmp04,
     &     atmp05,atmp06,
     &     atmp07,
     &     aE0tmp, oE0tmp,
     &     aEVAPORtmp,oEVAPORtmp,
     &     aSOLAR1tmp,oSOLAR1tmp,
     &     aSOLAR3tmp,oSOLAR3tmp,
     &     aDMUA1tmp,oDMUA1tmp,
     &     aDMVA1tmp,oDMVA1tmp)


      END SUBROUTINE AG2OG_oceans
#else

      SUBROUTINE AG2OG_oceans
!@sum  AG2OG_oceans: all atmospheric arrays used in the subr. OCEANS are gathered
!       on the atmospheric grid, interpolated to the ocean grid, and scattered
!!      on the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM, oJM=>JM
     *                , oDXYPO=>DXYPO, OXYP
     *                , IVSPO=>IVSP,IVNPO=>IVNP

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
#if defined (TRACERS_OceanBiology) && !defined (TRACERS_GASEXCH_ocean)
      USE OCN_TRACER_COM, only: NTM
#else
      USE TRACER_COM, only: NTM
#endif
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE OCEANR_DIM, only : ogrid

      USE SEAICE_COM, only : aRSI=>RSI
      USE FLUXES, only : aSOLAR=>SOLAR, aE0=>E0, aEVAPOR=>EVAPOR
     *     , aRUNOSI=>RUNOSI,aERUNOSI=>ERUNOSI,aSRUNOSI=>SRUNOSI
     *     , aFLOWO=>FLOWO,aEFLOWO=>EFLOWO, aAPRESS=>APRESS
     *     , aMELTI=>MELTI,aEMELTI=>EMELTI,aSMELTI=>SMELTI
     *     , aDMUA=>DMUA, aDMVA=>DMVA
     *     , aGMELT=>GMELT, aEGMELT=>EGMELT
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , aTRFLOWO=>TRFLOWO, aTREVAPOR=>TREVAPOR
     *     , aTRUNOSI=>TRUNOSI, aTRMELTI=>TRMELTI
     *     , aTRGMELT=>TRGMELT
#ifdef TRACERS_DRYDEP
     *     , aTRDRYDEP=>TRDRYDEP
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE FLUXES, only : aTRGASEX=>TRGASEX
      USE TRACER_GASEXCH_COM, only: tracflx
#endif
#ifdef OBIO_RAD_coupling
      USE RAD_COM, only : avisdir    => FSRDIR
     *                   ,asrvissurf => SRVISSURF
     *                   ,avisdif    => FSRDIF
     *                   ,anirdir    => DIRNIR
     *                   ,anirdif    => DIFNIR
#endif
#ifdef TRACERS_OceanBiology
      USE RAD_COM, only : aCOSZ1=>COSZ1
      USE PBLCOM, only : awind=>wsavg
#endif

      USE OFLUXES, only : oRSI, oSOLAR, oE0, oEVAPOR
     *     , oRUNOSI, oERUNOSI, oSRUNOSI
     *     , oFLOWO, oEFLOWO, oAPRESS
     *     , oMELTI, oEMELTI, oSMELTI
     *     , oDMUA,oDMVA
     *     , oGMELT, oEGMELT
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , oTRFLOWO, oTREVAPOR
     *     , oTRUNOSI, oTRMELTI
     *     , oTRGMELT
#ifdef TRACERS_DRYDEP
     *     , oTRDRYDEP
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE MODEL_COM, only: nstep=>itime
      USE OFLUXES, only : oTRGASEX
      USE TRACER_COM, only : ntm_gasexch
#endif
#ifdef OBIO_RAD_coupling
      USE obio_forc, only : ovisdir,ovisdif,onirdir,onirdif
#endif
#ifdef TRACERS_OceanBiology
      USE obio_forc, only : osolz,owind
#endif

      Use GEOM,  only : AXYP,aIMAXJ=>IMAXJ

      USE MODEL_COM, only : aFOCEAN_loc=>FOCEAN

      USE INT_AG2OG_MOD, only : INT_AG2OG

      IMPLICIT NONE

      INTEGER, PARAMETER :: NSTYPE=4
      INTEGER I,J, N
      INTEGER aJ_0,aJ_1, aI_0,aI_1
      INTEGER oJ_0,oJ_1,oI_0,oI_1

      REAL*8,
     * DIMENSION(aGRID%I_STRT:aGRID%I_STOP,aGRID%J_STRT:aGRID%J_STOP)::
     * aFact

      REAL*8, allocatable :: aweight(:,:),atmp(:,:),atmp2(:,:,:)
      real*8, dimension(oIM,oJM) :: otest_glob
      real*4, dimension(oIM,oJM) :: tr4

      allocate (aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )

      aJ_0 = aGRID%j_STRT
      aJ_1 = aGRID%j_STOP
      aI_0 = aGRID%I_STRT
      aI_1 = aGRID%I_STOP

      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP
      oI_0 = oGRID%i_STRT
      oI_1 = oGRID%i_STOP

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aRSI,oRSI, aWEIGHT)
      
      allocate (atmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )
      atmp=0.

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
            atmp(I,J) = aFLOWO(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(atmp,oFLOWO, aWEIGHT)
      deallocate(atmp)

      allocate (atmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )
      atmp=0.

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            atmp(I,J) = aEFLOWO(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(atmp,oEFLOWO, aWEIGHT)
      deallocate(atmp)

      allocate (atmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )
      atmp=0.

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            atmp(I,J) = aMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(atmp,oMELTI, aWEIGHT)
      deallocate(atmp)

      allocate (atmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )
      atmp=0.

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            atmp(I,J) = aEMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(atmp,oEMELTI, aWEIGHT)
      deallocate(atmp)

      allocate (atmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )
      atmp=0.

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            atmp(I,J) = aSMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(atmp,oSMELTI, aWEIGHT)
      deallocate(atmp)

      allocate (atmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )
      atmp=0.

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            atmp(I,J) = aGMELT(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(atmp,oGMELT, aWEIGHT)
      oGMELT(:,:) = oGMELT(:,:)*OXYP(:,:)
      deallocate(atmp)

      allocate (atmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )
      atmp=0.

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            atmp(I,J) = aEGMELT(I,J)*aFact(I,J)
            aFact(I,J) = 1.d0/AXYP(I,J) ! not used?
          END IF
        END DO
      END DO
      CALL INT_AG2OG(atmp,oEGMELT, aWEIGHT)
      oEGMELT(:,:) = oEGMELT(:,:)*OXYP(:,:)
      deallocate(atmp)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aRUNOSI,oRUNOSI, aWEIGHT)

      CALL INT_AG2OG(aERUNOSI,oERUNOSI, aWEIGHT)

      CALL INT_AG2OG(aSRUNOSI,oSRUNOSI, aWEIGHT)

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aAPRESS,oAPRESS, aWEIGHT)

      aWEIGHT(:,:) = 1.d0 - aRSI(:,:)
      CALL INT_AG2OG(aE0,oE0, aWEIGHT, NSTYPE,1)

      CALL INT_AG2OG(aEVAPOR,oEVAPOR, aWEIGHT, NSTYPE,1)

      CALL INT_AG2OG(aSOLAR,oSOLAR, aWEIGHT, 3,3,1)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aSOLAR,oSOLAR, aWEIGHT, 3,3,3)

#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      aWEIGHT(:,:) = 1.d0

      allocate (atmp2(NTM,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )
      atmp2=0.

      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
              atmp2(N,I,J) = aTRFLOWO(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(atmp2,oTRFLOWO, aWEIGHT, NTM)
      deallocate(atmp2)

      allocate (atmp2(NTM,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )
      atmp2=0.

      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              atmp2(N,I,J) = aTRMELTI(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(atmp2,oTRMELTI, aWEIGHT, NTM)
      deallocate(atmp2)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aTRUNOSI,oTRUNOSI, aWEIGHT, NTM)

      aWEIGHT(:,:) = 1.d0
      allocate (atmp2(NTM,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )
      atmp2=0.

      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
              atmp2(N,I,J) = aTRGMELT(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(atmp2,oTRGMELT, aWEIGHT, NTM)
      DO N=1,NTM
        oTRGMELT(N,:,:) = oTRGMELT(N,:,:)*OXYP(:,:)
      END DO

      aWEIGHT(:,:) = 1.d0 - aRSI(:,:)
      CALL INT_AG2OG(aTREVAPOR,oTREVAPOR, aWEIGHT, NTM, NSTYPE,1)

#ifdef TRACERS_DRYDEP
      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aTRDRYDEP,oTRDRYDEP, aWEIGHT, NTM, NSTYPE,1)
#endif
#endif
#endif

#ifdef TRACERS_GASEXCH_ocean
      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aTRGASEX,oTRGASEX, aWEIGHT,
     .               ntm_gasexch, NSTYPE,1)

      DO N=1,ntm_gasexch
        DO J=oJ_0,oJ_1
          DO I=oI_0,oI_1
             tracflx(I,J,N) = oTRGASEX(N,1,I,J) 
          ENDDO
        ENDDO
      ENDDO
#endif

#ifdef OBIO_RAD_coupling
      aWEIGHT(:,:) = aFOCEAN_loc(:,:)
      CALL INT_AG2OG(aVISDIR,aSRVISSURF,oVISDIR, aWEIGHT)

      CALL INT_AG2OG(aVISDIF,oVISDIF, aWEIGHT)

      CALL INT_AG2OG(aNIRDIR,oNIRDIR, aWEIGHT)

      CALL INT_AG2OG(aNIRDIF,oNIRDIF, aWEIGHT)
#endif

#ifdef TRACERS_OceanBiology
      aWEIGHT(:,:) = aFOCEAN_loc(:,:)
      CALL INT_AG2OG(aCOSZ1,oSOLZ, aWEIGHT)

      CALL INT_AG2OG(aWIND,oWIND, aWEIGHT)
#endif
      aWEIGHT(:,:) = 1.d0
      
      CALL INT_AG2OG(aDMUA,aDMVA,oDMUA,oDMVA, aWEIGHT,aFOCEAN_loc
     *              ,NSTYPE,1)

      deallocate(aweight)

      END SUBROUTINE AG2OG_oceans
#endif /*ag2og_oceans_bundle */


#ifdef OG2AG_OCEANS_BUNDLE      
      SUBROUTINE OG2AG_oceans
!@sum  OG2AG_oceans: ocean arrays for sea ice formation calculated in the
!       subr. OCEANS are gathered on the ocean grid, interpolated to the
!!      atmospheric grid, and scattered on the atmospheric ocean grid
!@auth Larissa Nazarenko, Denis Gueyffier
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM,oJM=>JM
     *                , oFOCEAN_loc=>FOCEAN_loc

#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm_atm=>ntm
#endif
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: ntm
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole
      USE OCEANR_DIM, only : ogrid

      USE MODEL_COM, ONLY : aFOCEAN_loc=>FOCEAN

      USE FLUXES, only : aDMSI=>DMSI,aDHSI=>DHSI,aDSSI=>DSSI
#ifdef TRACERS_ON
     *     , aDTRSI=>DTRSI
      Use GEOM,  only : IMAXJ
#endif

      USE OFLUXES, only : oDMSI, oDHSI, oDSSI
#ifdef TRACERS_OCEAN
     *     , oDTRSI
#endif

      use ocean, only : remap_O2A
      USE array_bundle
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole

      IMPLICIT NONE

      REAL*8, allocatable :: aWEIGHT(:,:),oWEIGHT(:,:)
      REAL*8, allocatable :: oDMSItmp(:,:,:),aDMSItmp(:,:,:)
      REAL*8, allocatable :: oDHSItmp(:,:,:),aDHSItmp(:,:,:)
      REAL*8, allocatable :: oDSSItmp(:,:,:),aDSSItmp(:,:,:)
      REAL*8, allocatable :: otmp(:,:,:,:),atmp(:,:,:,:)
      integer :: i,j,l,N1,N2
      type (lookup_str) :: lstr


      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      allocate(oDMSItmp(2,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(oDHSItmp(2,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(oDSSItmp(2,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      allocate(aDMSItmp(2,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(aDHSItmp(2,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(aDSSItmp(2,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))


      !-- initializing "lstr"
      call ab_init( lstr, 1, oIM,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO,
     &     aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO)


      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      call ab_add(lstr, oWEIGHT, aWEIGHT, shape(oWEIGHT), 'ij')

      oDMSItmp=oDMSI
      oDHSItmp=oDHSI
      oDSSItmp=oDSSI
      if (hasNorthPole(oGrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
        do L=1,2
           call copy_pole(oDMSItmp(L,:,oJM),oIM)
           call copy_pole(oDHSItmp(L,:,oJM),oIM)
           call copy_pole(oDSSItmp(L,:,oJM),oIM)
        enddo
      endif
      call ab_add(lstr, oDMSItmp, aDMSItmp, shape(oDMSItmp),
     &     'lij', oWEIGHT, aWEIGHT)
      call ab_add(lstr, oDHSItmp, aDHSItmp, shape(oDHSItmp),
     &     'lij', oWEIGHT, aWEIGHT) 
      call ab_add(lstr, oDSSItmp, aDSSItmp, shape(oDSSItmp),
     &     'lij', oWEIGHT, aWEIGHT) 

#if (defined TRACERS_OCEAN) && (defined TRACERS_ON)

      IF (NTM == NTM_ATM) THEN
         allocate(otmp(NTM,2,ogrid%I_STRT_HALO:ogrid%I_STOP_HALO,
     &        ogrid%J_STRT_HALO:ogrid%J_STOP_HALO))
         allocate(atmp(NTM,2,agrid%I_STRT_HALO:agrid%I_STOP_HALO,
     &        agrid%J_STRT_HALO:agrid%J_STOP_HALO))
         if (oIM .eq. aIM .and. oJM .eq. aJM) then  
            DO N1=1,NTM
               DO N2=1,2
                  DO J=ogrid%J_STRT,oGRID%J_STOP
                     DO I=ogrid%I_STRT,IMAXJ(J)
                        IF (aFOCEAN_loc(I,J).gt.0.) THEN
                           otmp(N1,N2,I,J) = oDTRSI(N1,N2,I,J)
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         else
            otmp=oDTRSI
            do N1=1,NTM
               do N2=1,2
                  call copy_pole(otmp(N1,N2,:,oJM),oIM)
               enddo
            enddo
         endif

         do N1=1,NTM
            do N2=1,2
               call ab_add(lstr,otmp(N1,N2,:,:),atmp(N1,N2,:,:),
     &              shape(otmp(N1,N2,:,:)),'ij',oWEIGHT,aWEIGHT)
            enddo
         enddo
         
      ELSE
C**** if number of ocean and atm. tracers are not the same
C**** do something in here

      END IF
#endif

c*   actual interpolation here
      call bundle_interpolation(lstr,remap_O2A)
c*

      if(hasNorthPole(agrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         do l=1,2
            call lon_avg(aDMSItmp(l,:,aJM),aIM)
            call lon_avg(aDHSItmp(l,:,aJM),aIM)
            call lon_avg(aDSSItmp(l,:,aJM),aIM)
         enddo
      endif

      if(hasSouthPole(agrid)
     &     .and. (aIM .ne. oIM .or. aJM .ne. oJM) ) then
         do l=1,2
            call lon_avg(aDMSItmp(l,:,1),aIM)
            call lon_avg(aDHSItmp(l,:,1),aIM)
            call lon_avg(aDSSItmp(l,:,1),aIM)
         enddo
      endif


      do j=agrid%j_strt,agrid%j_stop
      do i=agrid%i_strt,agrid%i_stop
        if(aFOCEAN_loc(i,j) > 0.) then
           aDMSI(:,i,j) = aDMSItmp(:,i,j)
           aDHSI(:,i,j) = aDHSItmp(:,i,j)
           aDSSI(:,i,j) = aDSSItmp(:,i,j)
        endif
      enddo
      enddo

#if (defined TRACERS_OCEAN) && (defined TRACERS_ON)
      aDTRSI=atmp
      deallocate(otmp,atmp)
#endif

      deallocate(oWEIGHT, aWEIGHT)
      deallocate(oDMSItmp,aDMSItmp)
      deallocate(oDHSItmp,aDHSItmp)
      deallocate(oDSSItmp,aDSSItmp)

      RETURN
      END SUBROUTINE OG2AG_oceans

#else
      SUBROUTINE OG2AG_oceans
!@sum  OG2AG_oceans: ocean arrays for sea ice formation calculated in the
!       subr. OCEANS are gathered on the ocean grid, interpolated to the
!!      atmospheric grid, and scattered on the atmospheric ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM,oJM=>JM
     *                , oFOCEAN_loc=>FOCEAN_loc

#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm_atm=>ntm
#endif
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: ntm
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE OCEANR_DIM, only : ogrid

      USE MODEL_COM, ONLY : aFOCEAN_loc=>FOCEAN

      USE FLUXES, only : aDMSI=>DMSI,aDHSI=>DHSI,aDSSI=>DSSI
#ifdef TRACERS_ON
     *     , aDTRSI=>DTRSI
#endif

      USE OFLUXES, only : oDMSI, oDHSI, oDSSI
#ifdef TRACERS_OCEAN
     *     , oDTRSI
#endif

      USE INT_OG2AG_MOD, only : INT_OG2AG

      IMPLICIT NONE

      REAL*8,
     * DIMENSION(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)::
     * oWEIGHT

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(oDMSI,aDMSI, oWEIGHT, 2)

      CALL INT_OG2AG(oDHSI,aDHSI, oWEIGHT, 2)

      CALL INT_OG2AG(oDSSI,aDSSI, oWEIGHT, 2)

#if (defined TRACERS_OCEAN) && (defined TRACERS_ON)

      IF (NTM == NTM_ATM) THEN

        oWEIGHT(:,:) = oFOCEAN_loc(:,:)
        CALL INT_OG2AG(oDTRSI,aDTRSI, oWEIGHT, 2, NTM)

      ELSE
C**** if number of ocean and atm. tracers are not the same
C**** do something in here

      END IF
#endif
      RETURN
      END SUBROUTINE OG2AG_oceans

#endif /* OG2AG_OCEANS_BUNDLE */

      MODULE IGOG_regrid_info
!@sum IGOG_info saves instances of hntrp_type for use in
!@+   ice <-> ocean regrids
      use hntrp_mod, only : hntrp_type
      implicit none
      save

      type(hntrp_type) :: hntrp_i2o_u ! ice u C -> ocn u C
      type(hntrp_type) :: hntrp_i2o_v ! ice v C -> ocn v C
      logical :: hntrp_i2o_uv_need_init = .true.

      type(hntrp_type) :: hntrp_o2i_u ! ocn u C -> ice u A
      type(hntrp_type) :: hntrp_o2i_v ! ocn v C -> ice v A
      logical :: hntrp_o2i_uv_need_init = .true.

      END MODULE IGOG_regrid_info


      SUBROUTINE IG2OG_oceans
!@sum IG2OG_oceans interpolates relevant DYNSI outputs to the
!@+   ocean grid
!@auth M. Kelley
      use hntrp_mod
      use IGOG_regrid_info
      use icedyn_com, only : iDMUI=>DMUI, iDMVI=>DMVI
      use ofluxes,    only : oDMUI,oDMVI
      USE ICEDYN, only : iIM=>imicdyn,iJM=>jmicdyn
      USE OCEAN,  only : oIM=>im,oJM=>jm
      USE ICEDYN, only : iDLATM=>DLATM
      USE OCEAN,  only : oDLATM=>DLATM
      USE ICEDYN,     only : iGRID=>grid_icdyn
      USE OCEANR_DIM, only : oGRID
      USE DOMAIN_DECOMP_1D, only : BAND_PACK,
     &      hasNorthPole, hasSouthPole
      implicit none
      real*8, dimension(:,:), allocatable ::
     &     ones_band,idmui_band,idmvi_band
      integer :: jmin,jmax

c assumption: every dynsi PE is an ocean PE
      if(.not. oGRID%have_domain) return

c      if(iIM.eq.oIM .and. iJM.eq.oJM) then
c        oDMUI(:,:) = iDMUI(:,:)
c        oDMVI(:,:) = iDMVI(:,:)
c        return
c      endif

      if(hntrp_i2o_uv_need_init) then
        call Init_Hntrp_Type(hntrp_i2o_u,
     &     iGRID, .5d0,iDLATM,
     &     oGRID, .5d0,oDLATM,
     &     0.d0)
        call Init_Hntrp_Type(hntrp_i2o_v,
     &     iGRID, 0.d0,iDLATM,
     &     oGRID, 0.d0,oDLATM,
     &     0.d0,
     &     JMA_4interp=iJM-1,JMB_4interp=oJM-1) ! for secondary lats
        hntrp_i2o_uv_need_init = .false.
      endif

c regrid DMUI from ice C to ocn C
      jmin = hntrp_i2o_u%bpack%jband_strt
      jmax = hntrp_i2o_u%bpack%jband_stop
      ALLOCATE(iDMUI_band(iIM,jmin:jmax),ones_band(iIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_i2o_u%bpack, iDMUI, iDMUI_band)
      if(jmax == iJM) then
        iDMUI_band(:,iJM) = 0.
      endif
      if(iIM.eq.oIM .and. iJM.eq.oJM) then ! no interp necessary
        oDMUI(:,jmin:jmax) = iDMUI_band(:,jmin:jmax)
      else
        call HNTR8_band (ones_band, iDMUI_band, hntrp_i2o_u, oDMUI)
      endif
      deallocate(iDMUI_band,ones_band)

c regrid DMVI from ice C to ocn C
      jmin = hntrp_i2o_v%bpack%jband_strt
      jmax = hntrp_i2o_v%bpack%jband_stop
      ALLOCATE(iDMVI_band(iIM,jmin:jmax),ones_band(iIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_i2o_v%bpack, iDMVI, iDMVI_band)
      if(jmax == iJM) then
        iDMVI_band(:,iJM) = 0.
      endif
      if(iIM.eq.oIM .and. iJM.eq.oJM) then ! no interp necessary
        oDMVI(:,jmin:jmax) = iDMVI_band(:,jmin:jmax)
      else
        call HNTR8_band (ones_band, iDMVI_band, hntrp_i2o_v, oDMVI)
      endif
      deallocate(iDMVI_band,ones_band)

      if(hasNorthPole(oGRID)) then ! INT_AG2OG_Vector2 set them to zero
        oDMUI(:,oJM) = 0.0
        oDMVI(:,oJM) = 0.0
      endif

      return
      END SUBROUTINE IG2OG_oceans

      SUBROUTINE OG2IG_uvsurf
!@sum OG2IG_uvsurf interpolates ocean surface velocity to the
!@+   DYNSI A-grid
!@auth M. Kelley
      use hntrp_mod
      use IGOG_regrid_info
      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE ICEDYN, only     : iIM=>imicdyn,iJM=>jmicdyn
      USE OCEAN,  only     : oIM=>im,oJM=>jm
      USE FLUXES, only : atm_uosurf=>uosurf,atm_vosurf=>vosurf
      USE ICEDYN_COM, only : uosurf,vosurf
      USE ICEDYN, only : iDLATM=>DLATM
      USE OCEAN,  only : oDLATM=>DLATM
      USE ICEDYN,     only : iGRID=>grid_icdyn
      USE OCEANR_DIM, only : oGRID
      USE DOMAIN_DECOMP_1D, only : BAND_PACK, 
     &     hasNorthPole, hasSouthPole
      USE OCEAN, only : UO,VO
#ifndef CUBED_SPHERE
      USE ICEDYN_COM, only : pack_a2i
#endif
      implicit none
      real*8, dimension(:,:), allocatable ::
     &     ones_band,ocnu_band,ocnv_band
      integer :: jmin,jmax

c assumption: every dynsi PE is an ocean PE
      if(.not. oGRID%have_domain) return

#ifndef CUBED_SPHERE
c DYNSI grid == ATM grid, so simply replicate ATM copy
c (Note this requires that TOC2SST has been called first)
      call band_pack(pack_a2i, atm_uosurf, uosurf)
      call band_pack(pack_a2i, atm_vosurf, vosurf)
      return
#endif

c call HNTRP
      if(hntrp_o2i_uv_need_init) then
        call Init_Hntrp_Type(hntrp_o2i_u,
     &       oGRID, .5d0,oDLATM,
     &       iGRID, 0.d0,iDLATM,
     &       0.d0)
        call Init_Hntrp_Type(hntrp_o2i_v,
     &       oGRID, 0.d0,oDLATM,
     &       iGRID, 0.d0,iDLATM,
     &       0.d0,
     &       JMA_4interp=oJM-1) ! secondary ocean lats
        hntrp_o2i_uv_need_init = .false.
      endif
c regrid uosurf
      jmin = hntrp_o2i_u%bpack%jband_strt
      jmax = hntrp_o2i_u%bpack%jband_stop
      ALLOCATE(ocnu_band(oIM,jmin:jmax),ones_band(oIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_o2i_u%bpack, uo(:,:,1), ocnu_band)
      call HNTR8_band (ones_band, ocnu_band, hntrp_o2i_u, uosurf)
      deallocate(ocnu_band,ones_band)
c regrid vosurf
      jmin = hntrp_o2i_v%bpack%jband_strt
      jmax = hntrp_o2i_v%bpack%jband_stop
      ALLOCATE(ocnv_band(oIM,jmin:jmax),ones_band(oIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_o2i_v%bpack, vo(:,:,1), ocnv_band)
      call HNTR8_band (ones_band, ocnv_band, hntrp_o2i_v, vosurf)
      deallocate(ocnv_band,ones_band)
      if(hasNorthPole(igrid)) then
c Eventually, uosurf/vosurf will be interpolated directly to
c the ice B grid, which has no polar point.  For now, setting
c A-grid values at the polar point to zero
        uosurf(:,iJM) = 0.
        vosurf(:,iJM) = 0.
      endif
      return
      END SUBROUTINE OG2IG_uvsurf
