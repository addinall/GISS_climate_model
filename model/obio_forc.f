#include "rundeck_opts.h"

      MODULE obio_forc

#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only : kdm=>lmo 
#else
      USE hycom_dim_glob
#endif
      USE obio_dim


      implicit none


      integer, ALLOCATABLE, DIMENSION(:,:) :: ihra            !counter for daylight hours

!@var COSZ1 Mean Solar Zenith angle for curr. physics(not rad) time step
!@var WSAVG     COMPOSITE SURFACE WIND MAGNITUDE (M/S)
      real, ALLOCATABLE, DIMENSION(:,:)    :: osolz
      real, ALLOCATABLE, DIMENSION(:,:)    :: owind           !wind speed in ocean grid (see hycom2.f)

      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tirrq3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: avgq            !mean daily irradiance in quanta
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: atmFe
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: atmFe_glob      !surface iron deposition
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: alk             !alkalinity in 'umol/kg'
#ifdef TRACERS_Alkalinity
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: alk_glob        !alkalinity in 'umol/kg'
#endif

#ifdef OBIO_RAD_coupling
      real*8, ALLOCATABLE, DIMENSION(:,:)    :: ovisdir,ovisdif
     .                                         ,onirdir,onirdif
#endif
#ifndef OBIO_RAD_coupling
      real, ALLOCATABLE, DIMENSION(:,:,:,:,:):: Eda,Esa       !direct,diffuse downwelling irradiance
#endif


      real solz               !mean cosine solar zenith angle
      real sunz               !solar zenith angle
      common /brod1/ solz,sunz
!$OMP THREADPRIVATE(/brod1/)

#ifdef OBIO_RAD_coupling 
      real eda_frac,esa_frac
      common /frac_oasim/eda_frac(nlt),esa_frac(nlt)

      real ovisdir_ij,ovisdif_ij,onirdir_ij,onirdif_ij
      common /rada2o_ij/ ovisdir_ij,ovisdif_ij,onirdir_ij,onirdif_ij
!$OMP THREADPRIVATE(/rada2o_ij/)
#else
      real Eda2,Esa2
      common /beda2/ Eda2(nlt,nhn),Esa2(nlt,nhn)
!$OMP THREADPRIVATE(/beda2/)
#endif

      real Ed,Es 
      common /beds/  Ed(nlt),Es(nlt)
!$OMP THREADPRIVATE(/beds/)

      real wind               !surface wind from atmos
      common /bwind/ wind
!$OMP THREADPRIVATE(/bwind/)

      real tirrq                   !total mean irradiance in quanta
      common /blte/ tirrq(kdm)
!$OMP THREADPRIVATE(/blte/)

      real, parameter ::  tirrq_critical=10.      !in quanta threshold at compensation depth

      real rmud                    !downwelling irradiance average cosine
      common /bmud /rmud 
!$OMP THREADPRIVATE(/bmud/)

      real atmCO2

      real rhosrf       !surface air density which comes from PBL.f

      contains

!------------------------------------------------------------------------------
      subroutine alloc_obio_forc

      USE obio_dim
#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
      USE OCEANRES, only : idm=>imo,jdm=>jmo,kdm=>lmo
#else
      USE hycom_dim_glob
      USE hycom_dim, only : ogrid,i_0h,i_1h,j_0h,j_1h
#endif


      implicit none

#ifdef OBIO_ON_GARYocean
      INTEGER :: j_0h,j_1h,i_0h,i_1h

      I_0H = ogrid%I_STRT_HALO
      I_1H = ogrid%I_STOP_HALO
      J_0H = ogrid%J_STRT_HALO
      J_1H = ogrid%J_STOP_HALO
#endif

      ALLOCATE(osolz(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(owind(i_0h:i_1h,j_0h:j_1h))

      ALLOCATE(tirrq3d(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(   ihra(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(   avgq(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(    alk(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(atmFe(i_0h:i_1h,j_0h:j_1h,12),atmFe_glob(idm,jdm,12))
#ifdef TRACERS_Alkalinity
      ALLOCATE(alk_glob(idm,jdm,kdm))
#endif


#ifdef OBIO_RAD_coupling
      ALLOCATE(ovisdir(i_0h:i_1h,j_0h:j_1h)
     .        ,ovisdif(i_0h:i_1h,j_0h:j_1h)
     .        ,onirdir(i_0h:i_1h,j_0h:j_1h)
     .        ,onirdif(i_0h:i_1h,j_0h:j_1h))
#endif

      end subroutine alloc_obio_forc

      END MODULE obio_forc
