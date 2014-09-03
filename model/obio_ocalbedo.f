#include "rundeck_opts.h"      

      subroutine obio_ocalbedo(wind,solz,bocvn,xocvn,chl,
     .                         rod,ros,hycgr,vrbos,i,j)

***********************************************************************
***** this routine is used by both atmosphere and ocean at each (i,j)
***** where implicitly it is assumed that
***** wind,solz,bocvn,xocvn,chl are in the atmos gird if hycgr=.false. 
***** wind,solz,                are in the ocean gird if hycgr=.true.
***** when this routine is called from within ocean, it does not compute
***** albedo coefficients. Those are computed on the atmos grid.
***** when this routine is called from within ocean, it does not pass
***** reflectances
***********************************************************************

c  Computes ocean surface albedo from solar zenith angle (solz)
c  and wind speed (wind, m/s).
c  Albedo is provided as direct (albd) and diffuse (albs).
c  Derive surface reflectance as a function of solz and wind
c  Includes spectral dependence of foam reflectance derived from Frouin
c  et al., 1996 (JGR)
      USE CONSTANT, only : radian
      USE FILEMANAGER
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
#if (defined CHL_from_SeaWIFs) || (defined OBIO_RAD_coupling)
!wfac does not depend on (i,j) thus indept of grid choice
      USE RAD_COM, only : wfac  
#elif (defined  CHL_from_OBIO) ||  (defined TRACERS_OceanBiology)
      USE obio_incom, only : wfac 
#endif
#if (defined CHL_from_OBIO) || (defined TRACERS_OceanBiology)
      USE obio_incom, only : lam
#endif
      USE MODEL_COM,  only : nstep=>itime

      implicit none

      integer nl,i,j
      real*8 cn,rof,rosps,rospd,rtheta
      real*8 sintr,rthetar,rmin,rpls,sinrmin,sinrpls,tanrmin
      real*8 tanrpls,sinp,tanp,a,b

      real*8, intent(in)  :: wind, solz, chl
      real*8 :: sunz
 
      real*8, dimension(6), intent(out) :: bocvn
      real*8, dimension(6), intent(out) :: xocvn

      real*8 :: sum1, sum2, part_sum
      logical :: obio_reflectance, res
      integer, parameter :: nlt=33
      real*8, dimension(nlt) :: refl

!!!!!!!!!!  Boris' part !!!!!!!!!!!!!!!!!
!@sum Those are the weights which were obtained by normalizing solar flux
!@sum for the Lamda's range 0 - 4000 nm. We used Landau fitting function for it.
!@sum They are used to be used for getting 6 band approximation based on 33
!@sum Watson Gregg band calculations
!@sum band_6(j) = Sum(band_33(i)*weight(i))/Sum(part_sum(i))
C**** WHY IS WEIGHT ONLY DECLARED TO BE 31 AND NOT 33?
      real*8 :: weight(31) = (/0.0158378,0.0201205,0.0241885,0.0277778,
     .                       0.0307124,0.0329082,0.0343586,0.0351143,
     .                       0.0352609,0.0349008,0.0341389,0.0330742,
     .                       0.0317941,0.0303725,0.0288696,0.0273329,
     .                       0.0500921,0.0643897,0.0686573,0.0532013,
     .                       0.0416379,0.0330341,0.0265929,0.0217156,
     .                       0.0179725,0.0150596,0.0127618,0.0158128,
     .                       0.0232875,0.0313132,0.0184843/)

C**** gband is the distribution of the 31 bands for the 6 band GISS code
      integer :: gband(7) = (/ 1, 18, 19, 23, 26, 30, 32 /)
c      real*8 :: part_sum(6) = (/0.526854,0.0643897,0.196531,0.066281,
c     .                        0.066922,0.0497974/)
        real*8 :: lam8(nlt)
#ifdef CHL_from_SeaWIFs
       integer:: lam(nlt) = (/ 250, 325, 350, 375, 400
     .                       , 425, 450, 475, 500, 525
     .                       , 550, 575, 600, 625, 650
     .                       , 675, 700, 725, 775, 850
     .                       , 950,1050,1150,1250,1350
     .                       ,1450,1550,1650,1750,1900,2200,2900,3700/)
#endif
!!!!!!!!!! end Boris' part !!!!!!!!!!!!!!!!!
      logical vrbos,hycgr

      real*8 :: roair, rn, rod(nlt),ros(nlt)
      integer :: ngiss

      rn = 1.341d0  ! index of refraction of pure seawater
      roair = 1.2D3 ! density of air g/m3  SHOULD BE INTERACTIVE?

      sunz=acos(solz)/radian  !in degs

c  Foam and diffuse reflectance
      if (wind .gt. 4.0) then
        if (wind .le. 7.0) then
          cn = 6.2D-4 + 1.56D-3/wind
          rof = roair*cn*2.2D-5*wind*wind - 4.0D-4
        else
          cn = 0.49D-3 + 0.065D-3*wind
          rof = (roair*cn*4.5D-5 - 4.0D-5)*wind*wind
        endif
        rosps = 0.057d0
      else
        rof = 0.0
        rosps = 0.066d0
      endif
      
c  Direct
c   Fresnel reflectance for sunz < 40, wind < 2 m/s
      if (sunz .lt. 40.0 .or. wind .lt. 2.0) then
        if (sunz .eq. 0.0) then
          rospd = 0.0211d0
        else
          rtheta = sunz*radian
          sintr = sin(rtheta)/rn
          rthetar = asin(sintr)
          rmin = rtheta - rthetar
          rpls = rtheta + rthetar
          sinrmin = sin(rmin)
          sinrpls = sin(rpls)
          tanrmin = tan(rmin)
          tanrpls = tan(rpls)
          sinp = (sinrmin*sinrmin)/(sinrpls*sinrpls)
          tanp = (tanrmin*tanrmin)/(tanrpls*tanrpls)
          rospd = 0.5*(sinp + tanp)
        endif
      else
       !Empirical fit otherwise
        a = 0.0253d0
        b = -7.14D-4*wind + 0.0618d0
        rospd = a*exp(b*(sunz-40.0))
      endif

c  Reflectance totals
      do nl = 1,nlt
        ros(nl) = rosps + rof*wfac(nl)
        rod(nl) = rospd + rof*wfac(nl)
      enddo

      !lam is integer, lam8 is real8
      lam8=float(lam)

      if(vrbos.and..not.hycgr)then
      do nl=1,nlt
          write(*,'(a,i5,9d12.4)')'ocalbedo: A', 
     .    nl,lam8(nl),wfac(nl),wind,sunz,rof,rospd,rosps,rod(nl),ros(nl)
      enddo
      endif
      if(vrbos.and.hycgr)then
      do nl=1,nlt
          write(*,'(a,i5,8d12.4)')'ocalbedo O:', 
     .    nl,wfac(nl),wind,sunz,rof,rosps,rospd,rod(nl),ros(nl)
      enddo
      endif

      if (hycgr) return  !we do not compute albedo coefs
                         !from within ocean, but from atmos
      
!!!!!!!!!! Boris' part !!!!!!!!!!!!!!!!!
C**** get chlorophyll term
      ! function obio_reflectance calculates reflectance 
      ! as a function of chl and wavelength (lam)

      res = obio_reflectance(refl,chl,lam8,nlt,i,j)
 
!  transition between band33 and band6 approximation

! loop over giss radiation bands
      do ngiss=1,6
        sum1 = 0.0
        sum2 = 0.0
        do nl=gband(ngiss), gband(ngiss+1)-1 
          if (refl(nl).lt.0) refl(nl) = 0.0
#ifdef CHL0
          ros(nl) = ros(nl) 
#else
          ros(nl) = ros(nl) + refl(nl)
#endif
          sum1 = sum1+weight(nl)*rod(nl)
          sum2 = sum2+weight(nl)*ros(nl)
        enddo
        part_sum=sum(weight(gband(ngiss):gband(ngiss+1)-1))
        xocvn(ngiss) = sum1/part_sum
        bocvn(ngiss) = sum2/part_sum

      if (xocvn(ngiss).ge.1. .or. bocvn(ngiss).ge.1.) then
         print*, 'XOCVN/BOCVN greater than 1 at ngiss,i,j=',
     .                ngiss,i,j
         do nl=gband(ngiss), gband(ngiss+1)-1
           write(*,'(i5,4e12.4)')
     .     nl,rod(nl),ros(nl),refl(nl),weight(gband(ngiss))
         enddo
         stop
      endif
      enddo
    
!!!!!!!!!! end Boris' part !!!!!!!!!!!!!!!!!
      return
      end subroutine obio_ocalbedo
