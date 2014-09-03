#include "rundeck_opts.h"
!     MODULE CLOUD_DROP_PRED
!     USE CLOUDS_COM
!     END MODULE CLOUD_DROP_PRED

      SUBROUTINE GET_CC_CDNC(L,AIRM,DXYPJ,PL,TL,DSS,MCDNL1,MCDNO1)
      USE CLOUDS_COM
      USE TRACER_COM
      USE CONSTANT,only:mb2kg,by3 ,avog,bygasc,RGAS
      IMPLICIT NONE
      real*8 AIRM,EXPL,EXPO,WCDNO,WCDNL,rho
     *,MCDNL1,MCDNO1,amass,tams,smturb,DXYPJ,PL,TL
      real*8 SSM1,SSM2,SSM3,SSM4,SSM5,SSMAL,SSMAO,SSML,SSMO
      real*8 SSMD1,SSMD2,SSMD3
      integer, PARAMETER :: nt=17+ntm_soa/2
      real*8,dimension(nt)::DSS,DSU
      integer L,n

      do n = 1,nt
        DSU(n)=1.d-10                         
c       write(6,*)"In SM ROUTIN",DSS(n),n
      end do
      SSMAL=0.d0
      SSMAO=0.d0

C** add in terms for AMASS from other program to get aerosol mass conc.
C** amass is airmass in kg 
      amass=AIRM*mb2kg*DXYPJ  
C** This is air density in kg/m3 
      rho=1d2*PL/(RGAS*TL)
C*** DSU gives you aerosol mass in  kg/m3, DSU is in kg of species
C*** DSS/amass is mass mixing ratio of aerosol (kg/kg)
      tams=1.d0/amass*rho
      do n = 1,nt 
        DSU(n) =DSS(n)*tams
C** Special case if not including dust-sulfate hetchem reactions
        if (n.gt.9 .and. n.lt.18) then
          if (DSS(n).eq.1.d-10) DSU(n)=0.d0
        endif 
      enddo

C** Converting aerosol mass (kg/m3) to number (m-3) from Leon Rotstayn, based on IPCC Table 5
c Factor of 132.14/32.06 converts from sulfur to ammonium sulfate
c 1.69e11 converts from mass (kg/m3) to number conc. (/cm3) for assumed lognormal dist.
c 1.21e11 converts from hydrophilic mass (kg/m3) to number conc.(/cm3) for
c assumed lognormal distribution for carbonaceous aerosols.
 
      SSM1 = 1.69d11*(132.14d0/32.06d0)*DSU(1)         !all sulfate 
      SSM2=(DSU(2)/2169.d0)/(0.004189d0*(.44d0**3)   ) !seasalt 0.1-1 um range
      SSM4 = 1.21d11*((DSU(4))+ DSU(6))                !OCIA,BCIA: aged (1 day efolding time)
#ifdef TRACERS_AEROSOLS_SOA
     .      +1.21d11*(DSU(18)+DSU(19)+DSU(20)+DSU(21)) !SOA
#endif  /* TRACERS_AEROSOLS_SOA */
      SSM5 = 1.21d11*((0.8*DSU(5))+(0.6*DSU(7)))       !OCB,BCB: 80% and 60% as hydrophillic) 
!     SSMD1=(DSU(10)/2.5d3)/(0.004189d0*(0.75d0**3.))  !dust(clay) coated w/ sulf
!     SSMD2=(DSU(11)/2.65d3)/(0.004189d0*(2.2d0**3.))  !dust(silt1) coated w/ sulfate
!     SSMD3=(DSU(12)/2.65d3)/(0.004189d0*(4.4d0**3.))  !dust(silt2) coated w/ sulfate

C*** Old way of getting Na in cm-3 as from Lohmann et al. 1999, JGR
c     SSM1=DSU(1)/1769.d0                           !sulfate
c     SSM2=(DSU(2)/2169.)/(0.004189*(.44**3.))      !seasalt 0.1-1 um range
c     SSM3=(0.1d0*DSU(3)/2169.)/(0.004189*(1.7**3.))!seasalt 1-4 um range    
c     SSM4=((1.0d0*DSU(4))+(1.0d0*DSU(6)))/1000.d0   !OCI and BCI aged use 100% & 100% 
c     SSM5=((0.8d0*DSU(5))+(0.6d0*DSU(7)))/1000.d0   !OCB and BCB      use 80%  & 60%
   
c     SSML=SSMAL/(0.004189d0*(.050d0**3.d0))
c     SSMO=SSMAO/(0.004189d0*(.085d0**3.d0)) + SSM2
c     IF(SSML.le.100.d0) SSML=100.d0
c     IF(SSMO.le.50.d0) SSMO=50.d0
C** Land Na (cm-3) is from sulfate+OC+BC and is SSMAL
C** Ocean Na (cm-3) is from sulfate+OC+BC+Seasalt and is SSMAO
      SSMAL=SSM1+SSM4+SSM5 !+ (SSMD1+SSMD2+SSMD3)
      SSMAO=SSMAL+SSM2             
      IF(SSMAL.le.100.d0) SSMAL=100.d0
      IF(SSMAO.le.50.d0) SSMAO=50.d0

C** We use data from Texas based on Segal et al. 2004 from Leon R.
      MCDNL1 = 174.8d0 + (1.51d0*SSMAL**0.886d0)
      MCDNO1 = -29.6d0 + (4.917d0*SSMAO**0.694d0)
c     write(6,*)"CDNC for MC Clds",MCDNL1,MCDNO1,SSMAL,SSMAO,L
c     if (MCDNL1.gt.6000.d0.or.MCDNO1.gt.6000.d0)
c    *write(6,*) "Nc1st",MCDNL1,MCDNO1,DSU(1),DSU(2),DSU(4),DSU(5),
c    *DSU(6),DSU(7),L,SSMAL,SSMAO,SSM1,SSM2,SSM4,SSM5

      RETURN
      END SUBROUTINE GET_CC_CDNC 

      SUBROUTINE GET_CDNC(L,LHX,WCONST,WMUI,AIRM,WMX,DXYPJ,
     *FCLD,CAREA,CLDSAVL,DSS,SMFPML,PL,TL,OLDCDO,OLDCDL,
     *VVEL,SME,DSU,CDNL1,CDNO1)
      USE CLOUDS_COM
      USE TRACER_COM
      USE CONSTANT,only:mb2kg,LHE,LHS,RGAS
      IMPLICIT NONE
      real*8 CAREA,CLDSAVL,AIRM,WMX,SMFPML,OLDCDO,OLDCDL,VVEL
     *,SME,rho,PL,TL

      integer, PARAMETER :: nt=17+ntm_soa/2
      real*8,dimension(nt)::DSS,DSU
      real*8 EXPL,EXPO,WCDNO,WCDNL,CDNO0,CDNL0,
     *CCLD0,CCLD1,DCLD,dfn,CDNL1,CDNO1,amass,tams,smalphaf
     *,FCLD,WCONST,LHX,WMUI,DXYPJ
      real*8 SSM1,SSM2,SSM3,SSM4,SSM5,SSMAL,SSMAO,SSML,SSMO
      real*8 SSMD1,SSMD2,SSMD3
      
      integer L,n

      do n = 1,nt
        DSU(n)=1.d-10                       
c       write(6,*)"In SM LSS ROUTIN",DSS(n),n
      end do
      SSMAL=0.d0
      SSMAO=0.d0

!add in terms for AMASS from other program to get aerosol mass conc.
      amass=AIRM*mb2kg*DXYPJ  
C** This is air density in kg/m3
      rho=1d2*PL/(RGAS*TL)
C*** DSU gives you aerosol mass in  kg/m3, DSU is in kg of species
C*** DSS/amass is mass mixing ratio of aerosol (kg/kg)
      tams=1.d0/amass*rho
      do n = 1,nt
         DSU(n) =DSS(n)*tams
C** Special case if not including dust-sulfate hetchem reactions
        if (n.gt.9 .and. n.lt.18) then
          if (DSS(n).eq.1.d-10) DSU(n)=0.d0
        endif 
      enddo

C** use CTEI effect for CDNC as in Menon et al. 2002 JAS
      smalphaf=(1.d0+ 2.d0*SMFPML)*0.5d0 
C** Converting aerosol mass (kg/m3) to number (m-3) from Leon Rotstayn, based on IPCC Table 5
c Factor of 132.14/32.06 converts from sulfur to ammonium sulfate
c 1.69e11 converts from mass (kg/m3) to number conc. (/cm3) for assumed lognormal dist.
c 1.21e11 converts from hydrophilic mass (kg/m3) to number conc.(/cm3) for
c assumed lognormal distribution for carbonaceous aerosols.

      SSM1 = 1.69d11*(132.14d0/32.06d0)*DSU(1)            !all sulfate
      SSM2=(DSU(2)/2169.d0)/(0.004189d0*(.44d0**3)   )    !seasalt 0.1-1 um range
      SSM4 = 1.21d11*((DSU(4))+ DSU(6))               !OCIA,BCIA: aged (1 day efolding time)
#ifdef TRACERS_AEROSOLS_SOA
     .      +1.21d11*(DSU(18)+DSU(19)+DSU(20)+DSU(21)) !SOA
#endif  /* TRACERS_AEROSOLS_SOA */
      SSM5 = 1.21d11*((0.8*DSU(5))+(0.6*DSU(7)))      !OCB,BCB: 80% and 60% as hydrophillic)
c     SSMD1=(DSU(10)/2.5d3)/(0.004189d0*(0.75d0**3.))  !dust(clay) coated w/ sulf
c     SSMD2=(DSU(11)/2.65d3)/(0.004189d0*(2.2d0**3.))  !dust(silt1) coated w/ sulfate
c     SSMD3=(DSU(12)/2.65d3)/(0.004189d0*(4.4d0**3.))  !dust(silt2) coated w/ sulfate

C** Old way based on Lohmann et al. 1999, JGR  !Converting aerosol mass to number
c     SSM1=DSU(1)/1769.d0                       !sulfate
c     SSM2=(DSU(2)/2169.)/(0.004189*(.44**3.))  !seasalt 0.1-1 um range
c     SSM3=(0.1d0*DSU(3)/2169.)/(0.004189*(1.7**3.))  !seasalt 1-4 um range    
c     SSM4=((1.0d0*DSU(4))+(1.0d0*DSU(6)))/1000.d0   !OCI and BCI aged use 100% & 100% 
c     SSM5=((0.8d0*DSU(5))+(0.6d0*DSU(7)))/1000.d0   !OCB and BCB      use 80%  & 60%
c     SSMAL=SSM1+SSM4+SSM5
c     SSMAO=SSMAL              !add in larger size when including TKE effects
c     SSML=SSMAL/(0.004189d0*(.050d0**3.d0))
c     SSMO=SSMAO/(0.004189d0*(.085d0**3.d0)) + SSM2

C** Land Na (cm-3) is from sulfate+OC+BC and is SSMAL
C** Ocean Na (cm-3) is from sulfate+OC+BC+Seasalt and is SSMAO
      SSMAL=SSM1+SSM4+SSM5 !+ (SSMD1+SSMD2+SSMD3)
      SSMAO=SSMAL+SSM2
      IF(SSMAL.le.100.d0) SSMAL=100.d0
      IF(SSMAO.le.50.d0) SSMAO=50.d0
C** Here we use Gultepe's paramet for CDNC = f(Na)
      EXPL=(298.d0*log10(SSMAL))-595.d0
      EXPO=(162.d0*log10(SSMAO))-273.d0

      IF (EXPO.LT.10.d0) EXPO=10.d0
      IF (EXPL.LT.10.d0) EXPL=10.d0
      WCDNO= EXPO*smalphaf
      WCDNL= EXPL*smalphaf
      CDNL0=OLDCDL   !term initialised to 10. CLOUDS_DRV                  
      CDNO0=OLDCDO   
C** Using the new CDNC scheme where we calculate it as a function of
C** gas phase sulfate and cloud area changes
      CCLD0 = CLDSAVL   !cld frac from previous time step   was 3D
      CCLD1 = 1.d0-CAREA  !CLDSS(L)cld frac from present time step 
      DCLD = CCLD1-CCLD0     ! cloud fraction change
      dfn = 0.0d0
      IF(LHX.EQ.LHE.AND.WMX/FCLD.GE.WCONST*1.d-3) dfn=0.1d0  
      IF(LHX.EQ.LHS.AND.WMX/FCLD.GE.WMUI*1.d-3) dfn=0.1d0

C** If previous time step is clear sky
      if (CCLD0.eq.0.d0) then
        CDNL1 = WCDNL        !EXPL1 = N from present time step
        CDNO1 = WCDNO
C** If previous time step is cloudy then depending on cld frac change
      elseif (DCLD.le.0.d0) then
        CDNL1= CDNL0 - dfn*CDNL0  !EXPL =  N from previous time step
        CDNO1= CDNO0 - dfn*CDNO0  !dfn=fraction of cloud area raining
      elseif (DCLD.gt.0.d0) then
        CDNL1 = (((CDNL0*CCLD0)+(WCDNL*DCLD))/CCLD1) - dfn*CDNL0
        CDNO1 = (((CDNO0*CCLD0)+(WCDNO*DCLD))/CCLD1) - dfn*CDNO0
      endif

      IF (CDNL1.le.10.d0) CDNL1=10.d0
      IF (CDNO1.le.10.d0) CDNO1=10.d0
c     if(CDNL1.gt.2000.d0.or.CDNO1.gt.2000.d0) 
c     write(6,*) "Nc1st_ST",CDNL1,CDNO1,DSU(1),DSU(2),DSU(4),DSU(5),
c    *DSU(6),DSU(7),L,SSMAL,SSMAO,SSM1,SSM2,SSM4,SSM5,smalphaf 
c    *,EXPO,EXPL,dfn,DCLD,CCLD0,CDNL0,CDNO0 
c     STOP "CLD_AEROSOLS_Menon_TC"
      RETURN
      
      END SUBROUTINE GET_CDNC 


      SUBROUTINE GET_QAUT(L,PL,TL,FCLD,WMX,SCDNCW,RCLD,RHOW,r6,r6c,
     *QCRIT,QAUT)
      USE TRACER_COM
      USE CONSTANT,only:TWOPI,GRAV,by6,by3,RGAS
      IMPLICIT NONE
      real*8 TL,WMX,SCDNCW,QAUT,RHOW,FCLD,RCLD,rho,r3c
      real*8 QCRIT,rcr,dynvis,PL,qcr,QAU,GAMA2
      real*8 epsi,epsis,betu6,betd6,bet6,cwc,r6c,r6,QAUT1,QAUT2
      real*8, PARAMETER :: bcon=1.15d23,  ! s-1 condensation rate const.
     &                     kap2=1.9d11    ! cm-3 s-1
      integer L

C**sm Tripoli and Cotton (1980), adapted from Jones Tech rpt.
C     using Qcrit for QAUT
C** This is air density in kg/m3 
       rho=1d2*PL/(RGAS*TL)
c     dynvis=2.4784d0*1.d-6*TL/(TL+393.16d0)
c     rcr=7.0d0
c     QCRIT=(2.d0/3.d0)*TWOPI*1.d9*((rcr*1.d-06)**3.d0)*SCDNCW/
c    *(.001d0*rho)
c
c     QAUT=(0.104d0*GRAV*0.55d0*((.001d0*rho)**(4.d0/3.d0))*
c    * ((WMX/(FCLD+1.d-30))**(7.d0/3.d0)))/
c    * (dynvis*((SCDNCW*1.d09)**(1.d0/3.d0)))

c     if(QAUT.le.0.) write(6,*)"QCR",QAUT,SCDNCW,WMX,l
C***Using Berrys 1967 scheme
c      GAMA2 = 0.35
c      QAUT=(GAMA2*rho*(WMX/(FCLD+1.d-20))**2)/(120.d0+((5.7d0
c    *    *SCDNCW*(FCLD+1.d-30))/(rho*WMX)))

C*** Use R6AUTO from Rotstayn and Liu (2005) GRL...
C***Using dispersion parameter for alhpa of 0.003 (liu and Daum, 2002, Nature)
       epsi=1.d0-0.7d0*exp(-0.003d0*SCDNCW)
       epsis=epsi**2
       betu6=(1.d0+3.d0*epsis)*(1.d0+4.d0*epsis)*(1.d0+5.d0*epsis)
       betd6=(1.d0+epsis)*(1.d0+2.d0*epsis)
       bet6 = (betu6/betd6)
       cwc = WMX/(FCLD+1.d-30)*1.d3*rho    ! gives in-cld water content in g/m3
       r6=RCLD*bet6   ! check to make sure it is in um
       r6c=4.09d-4*(bcon**by6)*(SCDNCW**by6)/(cwc**by3)
       r3c=r6c/(bet6**by6)*1.d-6   ! gives r in meters 
       QCRIT = (2.d0*TWOPI/3.d0)*rhow*r3c**3*SCDNCW*1.d6/rho   !Qcrit unitless: r,N in m
       QAUT1=(3.d0/(2.d0*TWOPI*rhow))**2
       QAUT2=kap2*bet6*(cwc**3)/SCDNCW
C** Units of kg m^-3 s^-1
       QAU=QAUT1*QAUT2*1.d-9
       QAUT=QAU/rho     ! QAUT is in s-1
C*** In main CLOUDS2.f if qc is gt. qcri start Qaut otherwise it is 0        
c     if (RCLD.gt.5.d0)
c     write(6,*)"QAUT",r6c,r6,RCLD,cwc,WMX,FCLD,SCDNCW,QCRIT,QAUT
      RETURN
    
      END SUBROUTINE GET_QAUT

      SUBROUTINE GET_CDNC_UPD(L,LHX,WCONST,WMUI,WMX,FCLD,
     *CLDSSL,CLDSAVL,VVEL,SME,DSU,SMFPML,OLDCDO,OLDCDL,CDNL1,CDNO1)
      USE CLOUDS_COM
      USE CONSTANT,only:LHE,LHS
      USE TRACER_COM
      IMPLICIT NONE
      real*8 ::CLDSSL,CLDSAVL,WMX
     *,SMFPML,OLDCDO,OLDCDL,VVEL,SME
      real*8 EXPL,EXPO,WCDNO,WCDNL,CDNO0,CDNL0,
     *CCLD0,CCLD1,DCLD,dfn,CDNL1,CDNO1,smalfaf,FCLD
     *,LHX,WMUI,WCONST
      integer, PARAMETER :: nt=17+ntm_soa/2
      real*8,dimension(nt)::DSU

      real*8 SSM1,SSM2,SSM3,SSM4,SSM5,SSMAL,SSMAO,SSML,SSMO
      real*8 SSMD1,SSMD2,SSMD3

      integer L

      SSMAL=0.d0
      SSMAO=0.d0
      smalfaf=(1.d0+ 2.d0*SMFPML)*0.5d0   !SMFPM was 3D
C** Converting aerosol mass (kg/m3) to number (m-3) from Leon Rotstayn, based on IPCC Table 5
c Factor of 132.14/32.06 converts from sulfur to ammonium sulfate
c 1.69e11 converts from mass (kg/m3) to number conc. (/cm3) for assumed lognormal dist.
c 1.21e11 converts from hydrophilic mass (kg/m3) to number conc.(/cm3) for
c assumed lognormal distribution for carbonaceous aerosols.

      SSM1 = 1.69d11*(132.14d0/32.06d0)*DSU(1)            !all sulfate
      SSM2=(DSU(2)/2169.d0)/(0.004189d0*(.44d0**3)   )    !seasalt 0.1-1 um range
      SSM4 = 1.21d11*((DSU(4))+ DSU(6))               !OCIA,BCIA: aged (1 day efolding time)
#ifdef TRACERS_AEROSOLS_SOA
     .      +1.21d11*(DSU(18)+DSU(19)+DSU(20)+DSU(21)) !SOA
#endif  /* TRACERS_AEROSOLS_SOA */
      SSM5 = 1.21d11*((0.8*DSU(5))+(0.6*DSU(7)))      !OCB,BCB: 80% and 60% as hydrophillic)
c     SSMD1=(DSU(10)/2.5d3)/(0.004189d0*(0.75d0**3.))  !dust(clay) coated w/ sulf
c     SSMD2=(DSU(11)/2.65d3)/(0.004189d0*(2.2d0**3.))  !dust(silt1) coated w/ sulfate
c     SSMD3=(DSU(12)/2.65d3)/(0.004189d0*(4.4d0**3.))  !dust(silt2) coated w/ sulfate

C** Old way based on Lohmann et al. 1999, JGR  !Converting aerosol mass to number
C** Converting aerosol mass to number
c     SSM1=DSU(1)/1769.d0                       !sulfate
c     SSM2=(DSU(2)/2169.)/(0.004189*(.44**3.))  !seasalt 0.1-1 um range
c     SSM3=(DSU(3)/2169.)/(0.004189*(1.7**3.))  !seasalt 1-4 um range    
c     SSM4=((1.0d0*DSU(4))+(1.0d0*DSU(6)))/1000.d0   !OCI and BCI aged use 100% & 1000% 
c     SSM5=((0.8d0*DSU(5))+(0.6d0*DSU(7)))/1000.d0   !OCB and BCB      use 80%  & 60%
c     SSMAL=SSM1+SSM4+SSM5
c     SSMAO=SSMAL              !add in larger size when including TKE effects
c     SSML=SSMAL/(0.004189d0*(.050d0**3.d0))
c     SSMO=SSMAO/(0.004189d0*(.085d0**3.d0)) + SSM2

C** Land Na (cm-3) is from sulfate+OC+BC and is SSMAL
C** Ocean Na (cm-3) is from sulfate+OC+BC+Seasalt and is SSMAO
      SSMAL=SSM1+SSM4+SSM5 !+ (SSMD1+SSMD2+SSMD3)
      SSMAO=SSMAL+SSM2

      IF(SSMAL.le.100.d0) SSMAL=100.d0
      IF(SSMAO.le.50.d0) SSMAO=50.d0
C** Here we use Gultepe's paramet for CDNC = f(Na)
      EXPL=(298.d0*log10(SSMAL))-595.d0
      EXPO=(162.d0*log10(SSMAO))-273.d0

      IF (EXPO.LT.10.d0) EXPO=10.d0
      IF (EXPL.LT.10.d0) EXPL=10.d0
      WCDNO= EXPO*smalfaf
      WCDNL= EXPL*smalfaf
      CDNL0=OLDCDL
      CDNO0=OLDCDO 

C** Using the new CDNC scheme where we calculate it as a function of
C** gas phase sulfate and cloud area changes
      CCLD0 = CLDSAVL   ! cld frac from previous time step  was 3D
      CCLD1 = CLDSSL    !cld frac from present time step  was 3D
      DCLD = CCLD1-CCLD0   ! cloud fraction change
      dfn = 0.d0
      IF(LHX.EQ.LHE.AND.WMX/FCLD.GE.WCONST*1.d-3) dfn=0.1d0
      IF(LHX.EQ.LHS.AND.WMX/FCLD.GE.WMUI*1.d-3) dfn=0.1d0

C** If previous time step is clear sky
      if (CCLD0.eq.0.d0) then
        CDNL1 = WCDNL        !EXPL1 = N from present time step
        CDNO1 = WCDNO
C** If previous time step is cloudy then depending on cld frac change
      elseif (DCLD.le.0.d0) then
        CDNL1= CDNL0 - dfn*CDNL0  !EXPL =  N from previous time step
        CDNO1= CDNO0 - dfn*CDNO0  !dfn=fraction of cloud area raining
      elseif (DCLD.gt.0.d0) then
        CDNL1 = (((CDNL0*CCLD0)+(WCDNL*DCLD))/CCLD1) - dfn*CDNL0
        CDNO1 = (((CDNO0*CCLD0)+(WCDNO*DCLD))/CCLD1) - dfn*CDNO0
      endif
c     if(CDNL1.gt.1600.d0.or.CDNO1.gt.1600.d0) 
c     write(6,*) "Nc2nd",CDNL1,CDNO1,DSU(1),DSU(2),DSU(4),DSU(5),
c    *DSU(6),DSU(7),L,SSMAL,SSMAO,SSM1,SSM2,SSM4,SSM5,smalfaf 
c    *,EXPO,EXPL,dfn,DCLD,CCLD0,CDNL0,CDNO0 
c     STOP "CLD_AEROSOLS_Menonc"
      IF (CDNL1.le.10.d0) CDNL1=10.d0
      IF (CDNO1.le.10.d0) CDNO1=10.d0

      RETURN

      END SUBROUTINE GET_CDNC_UPD
