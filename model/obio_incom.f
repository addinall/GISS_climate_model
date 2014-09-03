#include "rundeck_opts.h"

      MODULE obio_incom

! parameters and arrays neccessary for obio_init and obio_bioinit

      USE obio_dim

      implicit none

      integer, ALLOCATABLE, DIMENSION(:,:)   :: ir
      real,    ALLOCATABLE, DIMENSION(:,:,:) :: Fer,dicmod,dic
#ifndef pCO2_ONLINE
      real,    ALLOCATABLE, DIMENSION(:,:,:,:):: pco2tab
#endif

      real :: rmumax                !max phyto growth rate at 20oC, d/
      real :: rik                   !light saturation parameter umol quanta/m2/s
      common /bphy/ rmumax(nchl),rik(3,nchl)
      real :: obio_wsd,obio_wsh     !phyto sinking rate m/d, m/h
      common /bphy2/ obio_wsd(nchl),obio_wsh(nchl)
!$OMP THREADPRIVATE(/bphy2/)

      real :: rkn,rks     !half-saturation constants for nitrogen, silica (uM)
      real :: rkf         !half-saturation constant for iron (nM)
      common /brk/ rkn(nchl),rks(nchl),rkf(nchl)

      real rad, pi2
      common /bcnst/ rad,pi2

      real :: pHsfc,pHmin,pHmax   !pH at surface,minimun for iteration,
                                  !max for iteration
      common /cpH/pHsfc,pHmin,pHmax

      real Pdeep,detdeep,cardeep
      common /bpdp/ Pdeep(ntyp)            !deep BC
      common /bddp/ detdeep(ndet)          !detrital deep BC
      common /bcdp/ cardeep(ncar)          !carbon deep BC

      real :: cchl        !C:chl ratio g:g for 3 light states
      common /bcchl/ cchl(3)

      real :: cnratio     !C:N ratio
      real :: csratio     !C:Si ratio
      real :: cfratio     !C:Fe ratio
      common /beratio/ cnratio,csratio,cfratio

      real :: mgchltouMC

      real :: wsdeth      !sinking rate of detritus (m/d)
      real :: remin       !detrital remineralization rate /d
      common /bkwsdet/ wsdeth(ndet),remin(ndet)

      real :: Fescavrate  !scavenging rate for dissolved iron
      common /bfescav/ Fescavrate(2)

C if CARBON == 1
      real, parameter :: excp=0.05              !excretion of DOC by phyto %growth
      real, parameter :: resp=0.05              !respiration of DIC by phyto %growth
      real, parameter :: excz=0.05/24.0         !excretion of DOC by zoopl/hr
      real, parameter :: resz=0.05/24.0         !respiration of DIC by zoopl/hr
      real, parameter :: phygross=1.0-(excp+resp) !factor to derive gross PP!

!change: March 15, 2010
!     real, parameter :: rlamdoc=0.017/24.0    !N-dependent DOC remin/hr
      real, parameter :: rlamdoc=0.005/24.0    !N-dependent DOC remin/hr

      real, parameter :: rkdoc1=0.3*10.0       !N-dep half sat uM(PO4) modified
                                               !for nitrate by mult*10.0,
                                               !based on Conkright et al. 1994
      real, parameter :: rkdoc2=15.0           !DOC-dep half-sat uM(C)
      real, parameter :: rlampoc=0.05/24.0     !detrital breakdown/hr
      real, parameter :: uMtomgm3=12.0         !conversion uM to mg/m3 C
      real, parameter :: Pzo=1.0*uMtomgm3/50.0 !zoopl half-sat for
                                               !DOC excretion mg/m3(chl,assuming
                                               !C:chl ratio of 50))
      real, parameter :: awan=0.337/(3.6E+5)   !piston vel coeff., from
                                               !Wanninkof 1992, but adjusted
                                               !by OCMIP, and converted from
                                               !cm/hr to m/s
      real, parameter :: stdslp=1013.25        !standard sea level pressure in mb

      real, parameter :: Rm=1.20/24.0          !max zoopl. growth rate/hr
                                               !increase to account for excretion
                                               !and respiration
C if CARBON /=1    parameter(Rm=1.0/24.0)      !max zoopl. growth rate/hr


      integer lam               !wavelength in nm
      common /blam/ lam(nlt)  

      real facirr               !array of factors to compute mean irradiance w/in water column
      common /bfac/ facirr(nh,nch,5,ncd)  

      real aw,bw                !absorption,scattering coefficients of water
      common /bwat/ aw(nlt),bw(nlt)

      real ac,bc                !absorption and scattering coefficients of chlorophyll
      common /bopt/ ac(nchl,nlt),bc(nchl,nlt)

!     real, parameter :: solFe=0.02    !solubility of iron: this is the default
!     real, parameter :: solFe=0.05    !solubility of iron
      real solFe                       !now defined in rundeck

c     parameter(bn=0.5,bs=0.5)        !N/chl and Si/chl ratios
      
      real bn,bf,cchlratio

!save from moved ifst parts
       integer, parameter :: it0inc=1,nt0=80/it0inc,isalinc=1,nsal=20
       integer, parameter :: idicinc=2,ndic=(650+idicinc)/idicinc
       integer, parameter :: itainc=2,nta=(500+itainc)/itainc


#ifndef pCO2_ONLINE
       !real pco2tab
       !common /bpco2tab/pco2tab(nt0,nsal,ndic,nta)
#endif

      integer nl450
      common/exifst1/nl450

      real excdom,bbw,Dmax,rd,ru,rmus,rmuu,rn,roair
      common/exifst2/excdom(nlt),bbw,Dmax,rd,ru,rmus,rmuu
     .             ,rn,roair

#ifndef OBIO_RAD_coupling
      !if obio-rad-coupling is defined then this part is done 
      !inside RAD_COM.f and RAD_DRV.f
      real wfac(nlt)
      common/exifst3/wfac
#endif

!define compensation depth
      real, parameter ::  zc = 75. ! in meters (from OCMIP)

#ifdef TRACERS_Alkalinity
      real, parameter ::  rain_ratio=0.07
      real, parameter ::  npratio   =16     ! N:P Redfield ratio
      real, parameter ::  cpratio   =117    ! C:P Redfield ratio,
                                            ! these values are from OCMIP 
                                            ! protocol (Najjar and Orr, 1998)
                   ! Yamanaka and Tajika(1996) suggest R=0.08, rCP=106

! sigma  fraction of net downward flux of organic matter across the
!        compensation depth that is in dissolved form, ie form of
!        export production that is in dissolved form
      real, parameter ::  sigma_Ca  = 0.67 ! (Yamanaka and Tajika, 1997)

      real, parameter ::  d_Ca = 3500.    ! in meters (Yamanaka and Tajika, 1996)
      real, parameter ::  kappa_Ca = 2.    ! (1/0.5years)^-1   !OCMIP
#endif

      contains

      subroutine alloc_obio_incom

#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only :idm=>imo,jdm=>jmo,kdm=>lmo
#else
      USE hycom_dim_glob
#endif

      implicit none

!*******NEED TO DEALLOCATE LATER: only the first two lines
      ALLOCATE(ir(idm,jdm))
      ALLOCATE(Fer(idm,jdm,kdm),dicmod(idm,jdm,kdm),dic(idm,jdm,kdm))


      end subroutine alloc_obio_incom


      END MODULE obio_incom
