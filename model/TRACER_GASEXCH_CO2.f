#include "rundeck_opts.h"
!@sum  TRACERS_GASEXCH_COM: module for ocean-atmosphere gas exchange
!@+    special case for CO2
!@auth Natassa Romanou
!@ver  1.0

      MODULE TRACER_GASEXCH_COM

      USE TRACER_COM, only : ntm=>ntm_gasexch    !tracers in air-sea gas exch

      implicit none

#include "dimension2.h"

      private

      public alloc_gasexch_com

      public
     .  tracflx !  tracer flux at air-sea intfc

      public
     .  tracflx1d

      public atrac_loc

      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: atrac_loc

      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: tracflx !  tracer flux at air-sea intfc
  
      real*8 tracflx1d(ntm)
      common /gasexch3/tracflx1d
!$OMP THREADPRIVATE(/gasexch3/)


      contains

!------------------------------------------------------------------------------
      subroutine alloc_gasexch_com

      USE TRACER_COM, only : ntm=>ntm_gasexch    !tracers in air-sea gas exch
      use domain_decomp_atm, only : agrid=>grid
#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
#else
      USE HYCOM_DIM, only : ogrid
#endif
      integer i_0h,i_1h,j_0h,j_1h

      i_0h=ogrid%I_STRT_HALO
      i_1h=ogrid%I_STOP_HALO
      j_0h=ogrid%J_STRT_HALO
      j_1h=ogrid%J_STOP_HALO

      ALLOCATE(tracflx(i_0h:i_1h,j_0h:j_1h,ntm))

      allocate(atrac_loc(agrid%i_strt_halo:agrid%i_stop_halo,
     &                   agrid%j_strt_halo:agrid%j_stop_halo,
     &                   ntm ))

      end subroutine alloc_gasexch_com

      END MODULE TRACER_GASEXCH_COM

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
         !ONLY FOR HYCOM
#ifndef OBIO_ON_GARYocean
      subroutine init_gasexch_co2
      USE hycom_atm, only : gtracer_loc,focean_loc
      USE TRACER_COM, only : ntm_gasexch  !tracers involved in air-sea gas exch
      USE TRACER_GASEXCH_COM, only : atrac_loc
      use domain_decomp_atm, only : agrid=>grid
      implicit none
      integer nt,i,j

      do j=agrid%j_strt,agrid%j_stop
      do i=agrid%i_strt,agrid%i_stop
      if (focean_loc(i,j).gt.0.) then
          do nt=1,ntm_gasexch
            GTRACER_loc(nt,1,i,j)=atrac_loc(i,j,nt)
          enddo
      endif
      enddo
      enddo

      end subroutine init_gasexch_co2
#endif

c ---------------------------------------------------------------------

      SUBROUTINE TRACERS_GASEXCH_ocean_CO2_PBL(tg1,ws,
     . alati,psurf,itr,trconstflx,byrho,Kw_gas,alpha_gas,
     . beta_gas,trsf,trcnst,ilong,jlat)
  

      USE CONSTANT, only:    rhows,mair
      USE TRACER_COM, only : trname,tr_mm,vol2mass,ntm_gasexch
      USE obio_incom, only : awan
#ifdef OBIO_ON_GARYocean
      USE MODEL_COM,  only : nstep=>itime
#else
      USE HYCOM_SCALARS, only : nstep
#endif
      
      implicit none

      integer :: ilong,jlat,itr
      real*8  :: tg1,ws,psurf,trconstflx,byrho,trsf,trcnst
      real*8  :: alati,Kw_gas,alpha_gas,beta_gas
      real*8  :: Sc_gas
      real*8, external :: sc_co2,sol_co2

!@var  psurf surface pressure
!@var  alati SSS at i,j
!@var  Kw_gas  gas exchange transfer velocity at i,j only over ocean
!@var  alpha_gas  solubility of gas
!@var  beta_gas  conversion term  that includes solubility
!@var  trconstflx  constant component of surface tracer flux
!      byrho   1/surface air density

!routine for PBL calculations of CO2 gas exchange

      !---------------------------------------------------------------
      !TRANSFER VELOCITY
      !---------------------------------------------------------------
      !Schmidt number for gas
       Sc_gas=sc_co2(tg1)

      !wind speed ws: magn. of surf. wind modified by buoyancy flux (m/s)
      !compute transfer velocity Kw only over ocean
      if (Sc_gas .le. 0.) then
        write(*,'(a,2i4,a,2f9.3)')
     .          'warning: Sc_gas negtv, at ',ilong,jlat,
     .          ', Sc_gas,temp_c=',Sc_gas,tg1
         Kw_gas=1.e-10
      else
         Kw_gas=(Sc_gas/660.d0)**(-0.5d0) * ws * ws * awan !units of m/s
      endif

      !---------------------------------------------------------------
      !gas SOLUBILITY
      !---------------------------------------------------------------
      !alpha --solubility of CO2 in seawater
      alpha_gas=sol_co2(tg1,alati)    !mol/m^3/picoatm
      alpha_gas = alpha_gas * 1.d-6   !mol,CO2/m3/uatm
      alpha_gas = alpha_gas * 1024.5  !Watson Gregg
      !---------------------------------------------------------------
      !psurf is in mb. 
      beta_gas = alpha_gas * psurf/1013.25      !stdslp and psurf in mb, no need to change units

      !trsf is really sfac = Kw_gas * beta_gas
      !the term 1.0d6 / vol2mass(ntm_gasexch) is needed to convert uatm -> kg,co2/kg,air 
      !in the denominator of alpha
      trsf = Kw_gas * beta_gas * 1.0d6 / vol2mass(ntm_gasexch)

      !trconstflx comes in from SURFACE.f and has units kg,co2/kg,air/m2
      !therefore trcnst needs to be multiplied by byrho before it is sent to  PBL.f
      trcnst = Kw_gas * alpha_gas * trconstflx * byrho     
     .                * 1.0d6 / vol2mass(ntm_gasexch)    

        if (ilong.eq.1. .and. jlat.eq.45) then
        write(*,'(a,3i7,11e12.4)')'PBL, TRACER_GASEXCH_CO2:',
!       write(*,'(a,3i7,11e12.4)')'44444444444444444444444',
     .   nstep,ilong,jlat,tg1,(Sc_gas/660.d0)**(-0.5d0),ws*ws,
     .   Kw_gas,alpha_gas,beta_gas,trsf,trcnst,trconstflx,byrho,rhows
        endif


      RETURN
      END SUBROUTINE TRACERS_GASEXCH_ocean_CO2_PBL

c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c used with TRACERS_GASEXCH_ocean_CO2 to compute transfer velocity for CO2
c
c taken from Watsons' code
c ---------------------------------------------------------------------
      REAL*8 FUNCTION sol_co2(pt,ps)
c-------------------------------------------------------------------
c
c     CO2 Solubility in seawater
c
c     pt:       temperature (degrees Celcius)
c     ps:       salinity    (o/oo)
c     sol_co2:  in mol/m3/pptv
c               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm
c-------------------------------------------------------------------

      REAL*8    pt,ps,ta,tk100,tk1002


      ta = (pt + 273.16d0)
      tk100 = ta*0.01d0
      tk1002 = tk100*tk100
      sol_co2 = exp(-162.8301d0 + 218.2968d0/tk100  +
     .          90.9241d0*log(tk100) - 1.47696d0*tk1002 +
     .          ps * (.025695d0 - .025225d0*tk100 +
     .          0.0049867d0*tk1002))

      END

      REAL*8 FUNCTION alpha_gas2(pt,ps)
c-------------------------------------------------------------------
c
c     CO2 Solubility in seawater
c
c     pt:       temperature (degrees Celcius)
c     ps:       salinity    (o/oo)
c     sol_co2:  in mol/m3/pptv
c               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm
c     Scaled by 1024.5 (density of sea water) and 10^6 to get
c     alpha_gas in  mol/m3/uatm
c-------------------------------------------------------------------

      REAL*8    pt,ps,sol_CO2

      alpha_gas2=sol_CO2(pt,ps) * 1.d-6 * 1024.5  !mol,CO2/m3/uatm
      
      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      REAL*8 FUNCTION sc_co2(t)
c---------------------------------------------------
c     CO2 Schmidt number 
c     as a function of temperature. 
c
c     t: temperature (degree Celcius)
c---------------------------------------------------
      IMPLICIT NONE 
      REAL*8  a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)
      REAL*8 t
c
      sc_co2 = 2073.1 - 125.62*t + 3.6276*t**2 - 0.043219*t**3
c
      RETURN 
      END 

