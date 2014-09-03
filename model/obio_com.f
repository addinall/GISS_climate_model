#include "rundeck_opts.h"
      MODULE obio_com
!@sum  obio_com contains the parameters, arrays and definitions
!@+    necessary for the OceanBiology routines
!@auth NR
!@ver  1.0e-11

      USE obio_dim

#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only : kdm=>lmo
#else
      USE hycom_dim_glob
      USE hycom_scalars, only: baclin
#endif

      implicit none

c --- dobio       activate Watson Gregg's ocean biology code
      logical dobio
      data dobio/.true./
c

      real, ALLOCATABLE, DIMENSION(:,:)    :: tzoo2d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tfac3d,wshc3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: Fescav3d
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: rmuplsr3d,rikd3d
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: acdom3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: gcmax         !cocco max growth rate
      real, ALLOCATABLE, DIMENSION(:,:)    :: pCO2          !partial pressure of CO2
      real, ALLOCATABLE, DIMENSION(:,:)    :: pCO2_glob
      real, ALLOCATABLE, DIMENSION(:,:)    :: cexpij        !detritus term (Pg,C/yr)
      real, ALLOCATABLE, DIMENSION(:,:)    :: cexp_glob
      real, ALLOCATABLE, DIMENSION(:,:)    :: caexpij       !CaCO3 export term (Pg,C/yr)
      real, ALLOCATABLE, DIMENSION(:,:)    :: caexp_glob
      real, ALLOCATABLE, DIMENSION(:,:)    :: pp2tot_day    !net pp total per day
      real, ALLOCATABLE, DIMENSION(:,:)    :: pp2tot_day_glob !net pp total per day
      real, ALLOCATABLE, DIMENSION(:,:)    :: tot_chlo      !tot chlorophyl at surf. layer
      real, ALLOCATABLE, DIMENSION(:,:)    :: tot_chlo_glob !tot chlorophyl at surf. layer
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: rhs_obio      !rhs matrix
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: chng_by       !integr tendency for total C
#ifndef OBIO_ON_GARYocean   /* NOT for Russell ocean */
      real, ALLOCATABLE, DIMENSION(:,:,:,:) :: tracav, tracav_loc
      real, ALLOCATABLE, DIMENSION(:,:,:)   :: plevav, plevav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: pCO2av, pCO2av_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2tot_dayav
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2tot_dayav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: cexpav, cexpav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: caexpav, caexpav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: ao_co2fluxav,ao_co2fluxav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: ao_co2flux_loc  !ao CO2 on the ocean grid ***NOT for GASEXCH runs****
      real, ALLOCATABLE, DIMENSION(:,:) :: ao_co2flux_glob
#endif
#ifdef OBIO_ON_GARYocean
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: tracer_loc    !only for gary ocean
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: tracer        !only for gary ocean

      integer nstep0


      !test point
!!    integer, parameter :: itest=16, jtest=45    !equatorial Pacific                  2deg ocean
!!    integer, parameter :: itest=32, jtest=20    !southern ocean; Pacific          
      integer, parameter :: itest=1,  jtest=90     !equator Pacific
#else
!     integer, parameter :: itest=(220,320) equator Atlant; (245,275) 0.6S;274.5E Nino3
!     integer, parameter :: itest=316, jtest=258    !257.5E;-50.7S
      integer, parameter :: itest=243, jtest=1      !equator,dateline
      real :: diag_counter
#endif


      integer, parameter :: EUZ_DEFINED=1


      real, parameter :: rlamz=1.0            !Ivlev constant
      real, parameter :: greff=0.25           !grazing efficiency     
      real, parameter :: drate=0.05/24.0      !phytoplankton death rate/hr
      real, parameter :: dratez1=0.1/24.0     !zooplankton death rate/hr
      real, parameter :: dratez2=0.5/24.0     !zooplankton death rate/hr
      real, parameter :: regen=0.25           !regeneration fraction

      real ::  obio_deltath,obio_deltat       !time steps in hours because all rates are in hrs
      real ::  co2mon(26,12)        !26 years 1979-2004, 12 months

      integer npst,npnd   !starting and ending array index for PAR
      data npst,npnd /3,17/

      real WtoQ           !Watts/m2 to quanta/m2/s conversion
      common /bwq/ WtoQ(nlt)

! reduced rank arrays for obio_model calculations
      integer ihra_ij
      common /reducarr2/ihra_ij
!$OMP THREADPRIVATE(/reducarr2/)

      real temp1d,dp1d,obio_P,det,car,avgq1d,gcmax1d
     .    ,saln1d,p1d,alk1d,cexp,caexp,flimit
      common /reducarr1/temp1d(kdm),dp1d(kdm),obio_P(kdm,ntyp+n_inert)
     .                 ,det(kdm,ndet),car(kdm,ncar),avgq1d(kdm)
     .                 ,gcmax1d(kdm),saln1d(kdm),p1d(kdm+1)
     .                 ,alk1d(kdm),flimit(kdm,nchl,5)
!$OMP THREADPRIVATE(/reducarr1/)

      real atmFe_ij,covice_ij
      common /reducarr3/atmFe_ij,covice_ij
!$OMP THREADPRIVATE(/reducarr3/)


      integer inwst,inwnd,jnwst,jnwnd     !starting and ending indices 
                                          !for daylight 
                                          !in i and j directions
      common /bnwi/  inwst,inwnd,jnwst,jnwnd

      real    acdom
      common /bcdom/ acdom(kdm,nlt)       !absorptio coefficient of CDOM
!$OMP THREADPRIVATE(/bcdom/)


      real P_tend                             !bio tendency (dP/dt)
      common /bkpt/ P_tend(kdm,ntyp+n_inert)
!$OMP THREADPRIVATE(/bkpt/)

#ifdef TRACERS_Alkalinity
      real A_tend                             !alk tendency (dA/dt)
      common /akpt/ A_tend(kdm)
!$OMP THREADPRIVATE(/akpt/)
#endif

      real rmuplsr                            !growth+resp rate
      common /bgro/ rmuplsr(kdm,nchl)   
!$OMP THREADPRIVATE(/bgro/)

      real D_tend                             !detrtial tendency
      common /bdtend/ D_tend(kdm,ndet)     
!$OMP THREADPRIVATE(/bdtend/)

      real obio_ws                            !phyto sinking rate
      common /bkws/ obio_ws(kdm+1,nchl)   
!$OMP THREADPRIVATE(/bkws/)

      real tfac                               !phyto T-dependence
      common /btfac/ tfac(kdm)          
!$OMP THREADPRIVATE(/btfac/)

      real pnoice                    !pct ice-free
      common /bpnoice/ pnoice(kdm)  
!$OMP THREADPRIVATE(/bpnoice/)

      real wsdet                              !detrital sinking rate
      common /bkwsdet2/ wsdet(kdm+1,ndet) 
!$OMP THREADPRIVATE(/bkwsdet2/)

      real rikd                               !photoadaption state
      common /bikd/ rikd(kdm,nchl)      
!$OMP THREADPRIVATE(/bikd/)

      real tzoo                               !herbivore T-dependence
      common /bzoo/ tzoo                
!$OMP THREADPRIVATE(/bzoo/)

      real Fescav
      common /bscav/ Fescav(kdm)           !iron scavenging rate
!$OMP THREADPRIVATE(/bscav/)


C if NCHL_DEFINED > 3
      real wshc                            !cocco sinking rate
      common /bwshc/ wshc(kdm)             
!$OMP THREADPRIVATE(/bwshc/)
C endif


      real :: C_tend                        !carbon tendency
      common /bctend/ C_tend(kdm,ncar)      
!$OMP THREADPRIVATE(/bctend/)

      real :: pCO2_ij                 !partial pressure of CO2
      common /bco2ij/ pCO2_ij               
!$OMP THREADPRIVATE(/bco2ij/)

      real :: gro                           !realized growth rate
      common /bgro2D/ gro(kdm,nchl)         
!$OMP THREADPRIVATE(/bgro2D/)

      integer :: day_of_month, hour_of_day

      real :: rhs
      common /brhs/ rhs(kdm,ntrac,16)      !secord arg-refers to tracer definition 
                                        !we are not using n_inert (always ntrac-1)
                                        !second argument refers to process that
                                        !contributes to tendency 
!$OMP THREADPRIVATE(/brhs/)

      real :: pp2_1d          
      common /bpp2/ pp2_1d(kdm,nchl)           !net primary production
!$OMP THREADPRIVATE(/bpp2/)

#ifndef TRACERS_GASEXCH_ocean_CO2
#ifdef TRACERS_OceanBiology
      real*8 :: ao_co2flux
#endif
#endif
      integer kzc
      real*8 :: carb_old,iron_old    !prev timesetep total carbon inventory

#ifdef restoreIRON
!per year change dI/I, + for sink, - for source
      !!!real*8 :: Iron_BC = 0.002
      real*8 :: Iron_BC = -0.005
#endif
      contains

!------------------------------------------------------------------------------
      subroutine alloc_obio_com

      USE obio_dim

#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, only : get
      USE OCEANRES, only :idm=>imo,jdm=>jmo,kdm=>lmo
#else
      USE hycom_dim_glob 
      USE hycom_dim, only : i_0h,i_1h,j_0h,j_1h
#endif

      implicit none

#ifdef OBIO_ON_GARYocean
c**** Extract domain decomposition info
      INTEGER :: j_0h,j_1h,i_0h,i_1h

      I_0H = ogrid%I_STRT_HALO
      I_1H = ogrid%I_STOP_HALO
      J_0H = ogrid%J_STRT_HALO
      J_1H = ogrid%J_STOP_HALO


      ALLOCATE(tracer_loc(i_0h:i_1h,j_0h:j_1h,kdm,ntrac))
      ALLOCATE(tracer(idm,jdm,kdm,ntrac))
#endif

      ALLOCATE(tzoo2d(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(wshc3d(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(Fescav3d(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(rmuplsr3d(i_0h:i_1h,j_0h:j_1h,kdm,nchl),
     &            rikd3d(i_0h:i_1h,j_0h:j_1h,kdm,nchl))
      ALLOCATE(acdom3d(i_0h:i_1h,j_0h:j_1h,kdm,nlt))
      ALLOCATE(tfac3d(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(gcmax(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(pCO2(i_0h:i_1h,j_0h:j_1h))           
      ALLOCATE(pCO2_glob(idm,jdm))           
      ALLOCATE(pp2tot_day(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(pp2tot_day_glob(idm,jdm))
      ALLOCATE(tot_chlo(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(tot_chlo_glob(idm,jdm))
      ALLOCATE(rhs_obio(i_0h:i_1h,j_0h:j_1h,ntrac,16))
      ALLOCATE(chng_by(i_0h:i_1h,j_0h:j_1h,14))
#ifndef OBIO_ON_GARYocean   /* NOT for Russell ocean */
      ALLOCATE(ao_co2flux_loc(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(ao_co2flux_glob(idm,jdm))
      ALLOCATE(tracav(idm,jdm,kdm,ntrac))
      ALLOCATE(tracav_loc(i_0h:i_1h,j_0h:j_1h,kdm,ntrac))
      ALLOCATE(plevav(idm,jdm,kdm))
      ALLOCATE(plevav_loc(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(pCO2av(idm,jdm))
      ALLOCATE(pp2tot_dayav(idm,jdm))
      ALLOCATE(cexpij(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(cexp_glob(idm,jdm))           
      ALLOCATE(caexpij(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(caexp_glob(idm,jdm))           
      ALLOCATE(cexpav(idm,jdm))
      ALLOCATE(caexpav(idm,jdm))
      ALLOCATE(pCO2av_loc(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(pp2tot_dayav_loc(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(cexpav_loc(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(caexpav_loc(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(ao_co2fluxav(idm,jdm))
      ALLOCATE(ao_co2fluxav_loc(i_0h:i_1h,j_0h:j_1h))
#endif

      end subroutine alloc_obio_com

!------------------------------------------------------------------------------
      subroutine gather_pCO2

#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
#else
      USE HYCOM_DIM, only : ogrid
#endif
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA
 
      call pack_data( ogrid, pCO2, pCO2_glob )
#ifndef OBIO_ON_GARYocean
      call pack_data( ogrid, pp2tot_day, pp2tot_day_glob )
      call pack_data( ogrid, cexpij, cexp_glob )
      call pack_data( ogrid, caexpij, caexp_glob )
#endif

      end subroutine gather_pCO2

      subroutine gather_chl

#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
#else
      USE HYCOM_DIM, only : ogrid
#endif
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA

      call pack_data( ogrid, tot_chlo, tot_chlo_glob )

      end subroutine gather_chl

      END MODULE obio_com
