#include "rundeck_opts.h"
!  NOTE: Fastj photolysis scheme obtained from Oliver Wild (UCI) for
!  incorporation into the GISS 4x5 GCM with 10 tracer chemistry.
!  jlg. Sept. 98.
 
      subroutine photoj(nslon,nslat)
!@sum photoj (jv_trop.f): FAST J-Value code, troposphere only
!@+   (mjprather 6/96)
!@+   uses special wavelength quadrature spectral data (jv_spec.dat)
!@+   that includes only 289 nm - 800 nm  (later a single 205 nm add-on)
!@+   uses special compact Mie code based on Feautrier/Auer/Prather ver.
!@auth Oliver Wild/Lee Grenfell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on fastj000_M23p.f from model II)
!@calls CTM_ADJ,INT_PROF,PRTATM,JVALUE,JRATET

!**** GLOBAL parameters and variables:
      use domain_decomp_atm,only : grid,get
      use model_com, only    : LM,LS1
      use constant, only     : radian
      use dynamics, only     : ltropo
      use trchem_sHINDELL_com, ONLY: sza,tfastj,jfastj,jppj,zj,
     &            which_trop,szamax,U0,NCFASTJ,iprn,jprn,prnrts
 
      IMPLICIT NONE
 
!**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var i,j,k dummy loop variables
!@var NCFASTJ Number of levels in atmosphere
!@var maxl trop/strat boundary level
      integer, intent(IN) :: nslon, nslat
      integer   :: i,j,k,maxl
      logical   :: jay
      integer   :: J_0, J_1

      call get(grid, J_STRT=J_0,  J_STOP=J_1)
      jay = (nslat >= J_0 .and. nslat <= J_1)

#ifdef SHINDELL_STRAT_CHEM
      call stop_model('SHINDELL_STRAT_CHEM use TRCHEM_fastj2',255)
#endif
      select case(which_trop)
      case(0); maxl=ltropo(nslon,nslat)
      case(1); maxl=ls1-1
      case default; call stop_model('which_trop problem 3',255)
      end select
 
      zj(:,:)=0.d0
      jfastj(:,:)=0.d0
      ncfastj=2*LM+2 ! number of levels for clear-sky conditions.
      u0=dcos(sza*radian)
 
      if(SZA <= szamax) then
        call ctm_adj(nslon,nslat)       !define pres & aerosol profiles
        call int_prof(nslon,nslat,maxl) !define T & O3 profiles
        if(prnrts .and. nslon==iprn .and. nslat==jprn)
     &  call prtatm(nslon,nslat,jay)    ! Print out atmosphere
        call jvalue(nslon,nslat,maxl)   ! Calculate Actinic flux
        call jratet(nslon,nslat,maxl)   ! Calculate photolysis rates
        jfastj(:,:)=zj(:,:)             ! photolysis rates returned
      end if ! sza
      
      return
      end subroutine photoj
 
 
      subroutine ctm_adj(nslon,nslat)
!@sum ctm_adj to adapt the model input for the photolysis code; set up
!@+   the appropriate pressure and aerosol profiles. The pressure at the
!@+   bottom of CTM box (I,J,L) is given by   etaa(L) + etab(L)*P(I,J)
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
 
!**** GLOBAL parameters and variables:
      use model_com, only          : im,jm,lm
      use rad_com, only            : alb
      use trchem_sHINDELL_com, only: rcloudfj,jndlev,ncfastj,aer,
     & zfastj,jaddlv,jaddto,pfastj,rflect,odtmp,odsum,odmax,luselb,
     & nlbatm,zlbatm,iprn,jprn
 
      IMPLICIT NONE
 
!**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var i,j dummy variables
!@var MIEDX type of aerosol profile. See code for details.
!@var ASTRAT0 dummy variable
      integer, intent(IN) :: nslon, nslat
      integer             :: i,j,miedx
      real*8              :: astrat0

!-----------------------------------------------------------------------
!  Customise aerosol profile and scattering index
!
!     ASTRAT0  Extinction (/cm) always defined at 1000 nm = 1 micron
!     MAER1    Lowest level to put aerosol in (inclusive)
!     MAER2    Highest level to put aerosol in (inclusive)
!     MIEDX    Type of aerosol scattering, currently 6 set up
!
!-----------------------------------------------------------------------
! 1 Rayle 300 1.00 400 1.00  600 1.00  999 1.00 Rayleigh phase fn
! 2 IsoEX 300 1.00 400 1.00  600 1.00  999 1.00 isotropic exact opt.dep.
! 3 IsoEq 300 0.20 400 0.20  600 0.20  999 1.00 isotropic equiv.opt.dep.
! 4 SBkgd 310 2.77 400 2.38  600 1.58  999 .680 bckgrdsulfte,r/s=.09/.6
! 5 SVolc 310 2.69 400 2.55  600 2.18  999 1.46 volcnicsulfte,r/s=.08/.8
! 6 W_Cld 300 2.05 400 2.08  600 2.10  999 2.15 wtr cloud,r/gam=10./0.15
!-----------------------------------------------------------------------
! USER DEFINES THIS VARIABLE:
      miedx = 6

! Pressure Levels:
      j=0
      do i=1,lm
        j=2*i
        jndlev(i) = j ! centre of level we want J rate at
      enddo
 
! Associated altitudes
      zfastj(1) = 0.d0
      do i=2,ncfastj
        zfastj(i)= 16.d0*1.d+5*log10(pfastj(1)/pfastj(i))
      enddo
 
! Zero level indices
      do i=1,ncfastj
        jaddlv(i) = 0
        jaddto(i) = 0
      enddo
 
! Set surface albedo
      rflect = dble(1.- alb(nslon,nslat,1))
      rflect = max(0.d0,min(1.d0,rflect))
 
! Scale optical depths as appropriate
      odsum = 0.d0
      do i=1,lm
        odtmp(i) = rcloudfj(i,nslon,nslat)
        odsum    = odsum + odtmp(i)
      enddo
      if(odsum > odmax) then
        odsum = odmax/odsum
        do i=1,lm
          odtmp(i) = odtmp(i)*odsum
        enddo
      endif
 
! Lower boundary treatment; if cloud surface, zero ODs and set albedo
      nlbatm = 1
      if(luselb) then
        do i=1,lm
          if(odtmp(i) > zlbatm) nlbatm = (2*i)+1
          odtmp(i)=0.d0
        enddo
        if(nlbatm > 1) rflect = 0.6d0     ! Effective cloud albedo
      endif
 
! Write the AER profile directly - interpolate for boundary points...
! Or perhaps redefine if optical depths are provided directly
 
      aer(1:ncfastj) = 0.d0
 
      do i=1,LM
       if(odtmp(i) /= 0.d0) then
        astrat0 = odtmp(i)/(zfastj(jndlev(i)+1)-zfastj(jndlev(i)-1))
        aer(jndlev(i)-1) = aer(jndlev(i)-1) + astrat0*0.5d0
        aer(jndlev(i))   = aer(jndlev(i))   + astrat0
        aer(jndlev(i)+1) = aer(jndlev(i)+1) + astrat0*0.5d0
       endif
      enddo
 
      return
      end subroutine ctm_adj
 
 
      subroutine int_prof(nslon,nslat,maxl)
!@sum int_prof to interpolate T and O3 onto model grid. Currently only a
!@+   linear interpolation Oliver (23/05/97). Calculate the total number
!@+   density and the ozone profile.
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
 
!**** GLOBAL parameters and variables:
      use model_com, only: jm, month=>jmon
      USE GEOM, only: lat2d_dg
      use trchem_sHINDELL_com, only: o3_fastj,cboltz,dlogp,o3j,tj,dbc,
     &                   oref,tref,bref,do3,dmfastj,ncfastj,pfastj
 
      IMPLICIT NONE
 
!**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var i,j,l,m dummy loop variables
!@var pstd Approximate pressures of levels for supplied data
!@var tmp1,tmp2,ydgrd,month,tmpfra temporary variables
!@var maxl LTROPO or LS1-1, depending on which_trop variable
      integer, intent(IN) :: nslon, nslat, maxl
      integer             :: i,j,l,m
      real*8, dimension(51) :: pstd
      real*8              :: tmp1,tmp2,ydgrd,tmpfra
 
      pstd(1) = 1000.d0
      ydgrd=lat2d_dg(nslon,nslat)
      do i=2,51
        pstd(i) = pstd(i-1)*dlogp
      enddo
 
! Select appropriate monthly and latitudinal profiles
      m = max(1,min(12,month))
      l = max(1,min(18,(int(ydgrd)+99)/10))
 
! Linear interpolation
      do i=1,ncfastj
        if(pfastj(i) >= pstd(1)) then
          tmpfra = (pfastj(i)-pstd(2))/(pstd(1)-pstd(2))
          tj(i) = (tref(1,l,m)-tref(2,l,m))*tmpfra + tref(2,l,m)
          o3j(i) = (oref(1,l,m)-oref(2,l,m))*tmpfra + oref(1,l,m)
          dbc(i) = (bref(1)-bref(2))*tmpfra + bref(1)
        else
          do j=2,51
            if(PFASTJ(i) >= pstd(j)) then
              tmpfra = (pfastj(i)-pstd(j))/(pstd(j-1)-pstd(j))
!  Black carbon (assume constant above 80 km)
              if(j <= 41) then
                dbc(i) = (bref(j-1)-bref(j))*tmpfra + bref(j)
              else
                dbc(i) = bref(41)
              end if
!  Temperature (assume constant above 80 km)
              if(j <= 41) then
                tj(i) = (tref(j-1,l,m)-tref(j,l,m))*tmpfra + tref(j,l,m)
              else
                tj(i) = tref(41,l,m)
              endif
!  Ozone (extrapolate above 60 km)
              if(j <= 31) then
                o3j(i)=(oref(j-1,l,m)-oref(j,l,m))*tmpfra+oref(j,l,m)
              else
                tmp1=oref(31,l,m)*(oref(31,l,m)/oref(30,l,m))**(j-31)
                tmp2=oref(31,l,m)*(oref(31,l,m)/oref(30,l,m))**(j-32)
                o3j(i)=(tmp2-tmp1)*tmpfra + tmp1
              endif
              go to 10
            endif
          enddo
          write(6,*)'pfastj,i,pstd(51)',pfastj(i),i,pstd(51)
          if(pfastj(i) < pstd(51))
     &    call stop_model('Need mesosphere data in int_prof.',255)
 10       continue
        endif
      enddo
 
! overwrite troposphere levels of climatological O3 with GISS GCM O3:
      do i=1,2*maxl
        o3j(i)=o3_fastj(i)
      enddo
 
! Calculate total and O3 number densities
      do i=1,ncfastj
        dmfastj(i)  = pfastj(i)/(cboltz*tj(i))
        do3(i) = o3j(i)*1.d-6*dmfastj(i)
      enddo
 
      return
 1000 format(1x,i2,2x,f8.3,2x,f5.3,2x,f7.3,2x,f7.4)
      end subroutine int_prof
 
 
      SUBROUTINE JRATET(NSLON,NSLAT,maxl)
!@sum JRATET Calculate & print J-values. The loop in this routine
!@+   only covers the maxl levels actually needed by the CTM.
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only:FFF,TFASTJ,VALJ,NW1,NW2,jpdep,
     &                              TQQ,QQQ,zpdep,PFASTJ,jppj,zj,
     &                              jfacta,NJVAL,jind
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var i,j,k,l,jgas dummy loop variables
!@var QO2TOT,QO3TOT,QO31D,QO33P,QQQT,tfact temporary variables
!@var maxl the tropopause level
      INTEGER, INTENT(IN) :: nslon, nslat, maxl
      INTEGER i,j,k,l,jgas
      REAL*8 QO2TOT,QO3TOT,QO31D,QO33P,QQQT,tfact
      REAL*8 XSECO2,XSECO3,XSEC1D ! >>> FUNCTIONS <<<
C
      DO I=1,maxl
       VALJ(1) = 0.d0
       VALJ(2) = 0.d0
       VALJ(3) = 0.d0
       DO K=NW1,NW2                       ! Using model 'T's here
         QO2TOT= XSECO2(K,dble(TFASTJ(I)))
         VALJ(1) = VALJ(1) + QO2TOT*FFF(K,I)
         QO3TOT= XSECO3(K,dble(TFASTJ(I)))
         QO31D = XSEC1D(K,dble(TFASTJ(I)))*QO3TOT
         QO33P = QO3TOT - QO31D
         VALJ(2) = VALJ(2) + QO33P*FFF(K,I)
         VALJ(3) = VALJ(3) + QO31D*FFF(K,I)
       ENDDO
C
C------Calculate remaining J-values with T-dep X-sections
       DO J=4,NJVAL
         VALJ(J) = 0.d0
         TFACT = 0.d0
         L = jpdep(J)
         IF(TQQ(2,J).GT.TQQ(1,J)) TFACT = DMAX1(0.D0,DMIN1(1.D0,
     $        (TFASTJ(I)-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J)) ))
           DO K=NW1,NW2
           QQQT = QQQ(K,1,J-3) + (QQQ(K,2,J-3) - QQQ(K,1,J-3))*TFACT
           if(L.eq.0) then
             VALJ(J) = VALJ(J) + QQQT*FFF(K,I)
           else
           VALJ(J)=VALJ(J)+QQQT*FFF(K,I)*(1.d0+zpdep(K,L)*PFASTJ(2*i))
           endif
         ENDDO
       ENDDO

       DO jgas=1,jppj
         zj(i,jgas)=VALJ(jind(jgas))*jfacta(jgas)
       ENDDO
      ENDDO
C
      RETURN
      END SUBROUTINE JRATET
C
C
      SUBROUTINE PRTATM(NSLON,NSLAT)
!@sum PRTATM Print out atmosphere, calculate approx. columns.
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only:NLFASTJ,DO3,AER,DMFASTJ,DBC,ZZHT,
     &                              NCFASTJ,ZFASTJ,PFASTJ,TJ
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var i dummy variable
!@var COLO3,COLAX,COLO2,COLBC,ZKM,ZSTAR local variables for printing
      INTEGER, INTENT(IN) :: nslon, nslat
      INTEGER i
      REAL*4, DIMENSION(NLFASTJ) :: COLO3,COLAX,COLO2,COLBC
      REAL*4 ZKM,ZSTAR
C
C---Calculate columns, for diagnostic output only:
      COLO3(NCFASTJ) = DO3(NCFASTJ)*ZZHT
      COLAX(NCFASTJ) = AER(NCFASTJ)*ZZHT
      COLO2(NCFASTJ) = DMFASTJ(NCFASTJ)*0.20948*ZZHT
C     assuming 10 m2/g Black Carbon:
      COLBC(NCFASTJ) = DBC(NCFASTJ)*ZZHT*10.0
      DO I=NCFASTJ-1,1,-1
       COLO3(i)=COLO3(i+1)+(DO3(i)+DO3(i+1))*0.5*(ZFASTJ(i+1)-ZFASTJ(i))
       COLO2(i)=COLO2(i+1)+(DMFASTJ(i)+DMFASTJ(i+1))*0.5*
     $ (ZFASTJ(i+1)-ZFASTJ(i))*.20948
       COLAX(i)=COLAX(i+1)+(AER(i)+AER(i+1))*0.5*(ZFASTJ(i+1)-ZFASTJ(i))
       COLBC(i) = COLBC(i+1)+(DBC(i)+DBC(i+1))*0.5*
     $ (ZFASTJ(i+1)-ZFASTJ(i))*10.0
      ENDDO
      WRITE(6,1200) ' O3-column(DU)=',COLO3(1)/2.687E16,
     $              '  column aerosol @1000nm=',COLAX(1),
     $              '  column black carbon=',COLBC(1)
C---Print out atmosphere
      WRITE(6,1000)
      DO I=1,NCFASTJ
        ZKM = 1.d-5*ZFASTJ(I)
        ZSTAR = 16.d0*LOG10(1000.D0/PFASTJ(I))
        WRITE(6,1100) I,ZKM,ZSTAR,DMFASTJ(I),DO3(I),
     $  1.d6*DO3(I)/DMFASTJ(I),
     $  TJ(I),PFASTJ(I),AER(I),COLO3(I),COLAX(I),COLO2(I),COLBC(I)
      ENDDO
      RETURN
 1000 format(5X,'Zkm',3X,'Z*',8X,'M',8X,'O3',6X,'f-O3',5X,'T',7X,'P',
     $    7X,'AER-X',4X,'col-O3',2X,'col-AER',3X,'col-O2',3X,'col-BC')
 1100 format(1X,I2,0P,2F6.2,1P,2E10.3,0P,F7.3,F8.2,F10.4,1P,5E9.2)
 1200 format(A,F8.1,A,F8.4,A,1pE10.3)
      END SUBROUTINE PRTATM
C
C
      SUBROUTINE JVALUE(nslon,nslat,maxl)
!@sum JVALUE Calculate actinic flux at each level for current SZA value.
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
!@calls SPHERE,OPMIE,XSECO3,XSECO2
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: LM
      USE TRCHEM_Shindell_COM, only: NW1,NW2,FFF,WL,NCFASTJ,TJ,FL,
     &                          XQO2,XQO3
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var k,j dummy loop variables
!@var WAVE Effective wavelength of each wavelength bin
!@var AVGF   Attenuation of beam at each level for each wavelength
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var maxl the tropopause level
      INTEGER, INTENT(IN) :: nslon, nslat, maxl
      INTEGER j,k
      REAL*8 WAVE
      REAL*8 XSECO2,XSECO3 ! >>> FUNCTIONS <<<
      REAL*8, DIMENSION(LM) :: AVGF
C
      DO J=1,maxl
        DO K=NW1,NW2
         FFF(K,J) = 0.d0
        ENDDO
      ENDDO

C---Calculate spherical weighting functions
c
      CALL SPHERE
c
C---Loop over all wavelength bins
      DO K=NW1,NW2
        WAVE = WL(K)
        DO J=1,NCFASTJ
          XQO3(J) = XSECO3(K,dble(TJ(J)))
        ENDDO
        DO J=1,NCFASTJ
          XQO2(J) = XSECO2(K,dble(TJ(J)))
        ENDDO
c
C-----------------------------------------
        CALL OPMIE(K,WAVE,AVGF)
C-----------------------------------------
c
        DO J=1,maxl
          FFF(K,J) = FFF(K,J) + FL(K)*AVGF(J)
        ENDDO
      ENDDO
c
      RETURN
      END SUBROUTINE JVALUE
C
C
      FUNCTION XSECO3(K,TTT)
!@sum XSECO3  O3 Cross-sections for all processes interpolated across
!@+   3 temps
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
!@calls FLINT
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: TQQ,QO3
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var k passed index for wavelength bin
!@var TTT returned termperature profile
      INTEGER, INTENT(IN) :: k
      REAL*8 xseco3,FLINT ! >>> FUNCTIONS <<<
      real*8 TTT
c
      XSECO3  =
     F  FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2),QO3(K,1),QO3(K,2),QO3(K,3))
      RETURN
      END FUNCTION XSECO3
C
C
      FUNCTION XSEC1D(K,TTT)
!@sum XSEC1D  Quantum yields for O3 --> O2 + O(1D) interpolated across
!@+   3 temps
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
!@calls FLINT
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: TQQ,Q1D
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var k passed index for wavelength bin
!@var TTT returned termperature profile
      INTEGER, INTENT(IN) :: k
      REAL*8 xsec1d,FLINT ! >>> FUNCTIONS <<<
      real*8 TTT
c
      XSEC1D =
     F  FLINT(TTT,TQQ(1,3),TQQ(2,3),TQQ(3,3),Q1D(K,1),Q1D(K,2),Q1D(K,3))
      RETURN
      END FUNCTION XSEC1D
c
c
      FUNCTION XSECO2(K,TTT)
!@sum XSECO2 Cross-sections for O2 interpolated across 3 temps; No
!@+   S_R Bands yet!
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
!@calls FLINT
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: TQQ,QO2
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var k passed index for wavelength bin
!@var TTT returned termperature profile
      INTEGER, INTENT(IN) :: k
      REAL*8 xseco2,FLINT ! >>> FUNCTIONS <<<
      real*8 TTT
c
      XSECO2 =
     F  FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1),QO2(K,1),QO2(K,2),QO2(K,3))
      RETURN
      END FUNCTION XSECO2
c
c
      REAL*8 FUNCTION FLINT(TINT,T1,T2,T3,F1,F2,F3)
!@sum FLINT Three-point linear interpolation function
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
c
C**** GLOBAL parameters and variables:
C (none)
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var T1,T2,T3 passed temperature variables
!@var F1,F2,F3 passed X-section variables
!@var TINT returned temperature profile?
      REAL*8, INTENT(IN) :: T1,T2,T3,F1,F2,F3
      real*8 TINT
!!    REAL*8 FLINT ! >>> FUNCTIONS <<<

      IF (TINT .LE. T2)  THEN
        IF (TINT .LE. T1)  THEN
          FLINT  = F1
        ELSE
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        ENDIF
      ELSE
        IF (TINT .GE. T3)  THEN
          FLINT  = F3
        ELSE
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        ENDIF
      ENDIF
      RETURN
      END FUNCTION FLINT
c
c
      SUBROUTINE SPHERE
!@sum SPHERE Calculate spherical geometry; derive tangent heights, slant
!@+   path lengths and optical depth weighting for each layer. Not
!@+   called when SZA > 98 degrees.  Beyond 90 degrees, include treat-
!@+   ment of emergent beam (where tangent height is below altitude
!@+   J-value desired at).
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
c
C**** GLOBAL parameters and variables:
C
      USE CONSTANT, only: radius
      USE TRCHEM_Shindell_COM, only: U0,RZ,RQ,ZFASTJ,TANHT,nlbatm,WTAU,
     &                              NCFASTJ,ZZHT
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var airmas - statement function with dummy arguments Ux,HFASTJ
!@var xmu1,xmu2,xl,DIFF,WTING temp variables
!@var i,ii,j,k dummy loop index
      REAL*8 AIRMAS,Ux,HFASTJ,xmu1,xmu2,xl, diff, wting
      INTEGER i,ii,j,k
      AIRMAS(Ux,HFASTJ)=(1.0d0+HFASTJ)/SQRT(Ux*Ux+2.0d0*HFASTJ*(1.0d0-
     $       0.6817d0*EXP(-57.3d0*ABS(Ux)/SQRT(1.0d0+5500.d0*HFASTJ))/
     $                                  (1.0d0+0.625d0*HFASTJ)))
C
      RZ(1)=RADIUS*1.D2+ZFASTJ(1)
c
      DO II=2,NCFASTJ
        RZ(II) = RADIUS*1.D2 + ZFASTJ(II)
        RQ(II-1) = (RZ(II-1)/RZ(II))**2.
      END DO
      IF (U0.LT.0.0D0) THEN
        TANHT = RZ(nlbatm)/DSQRT(1.0D0-U0**2.)
      ELSE
        TANHT = RZ(nlbatm)
      ENDIF
C  Now go up from the surface calculating the weightings for each level
C  from the level we're at upwards - (WTAU(i,j),i=j,nc)
      DO 16 J=1,NCFASTJ
        DO K=1,NCFASTJ
          WTAU(K,J)=0.0D0
        END DO
C-----WEIGHTED LENGTHS TO SUN STARTING AT ZFASTJ(J);MU>0
        IF (RZ(J).LT.TANHT) GOTO 16
        XMU1=ABS(U0)
        DO I=J,NCFASTJ-1
          XMU2=DSQRT(max(0.d0,1.D0-RQ(I)*(1.D0-XMU1**2.)))
          XL=RZ(I+1)*XMU2-RZ(I)*XMU1
          WTAU(I,J)=WTAU(I,J)+XL*0.5D0
          WTAU(I+1,J)=WTAU(I+1,J)+XL*0.5D0
          XMU1=XMU2
        END DO
C-----SCALE HEIGHT AT TOP POINT
        WTAU(NCFASTJ,J)=WTAU(NCFASTJ,J) +
     &  ZZHT*AIRMAS(XMU1,ZZHT/(RADIUS*1.D2))
        IF (U0.GE.0.0D0) GOTO 16
C-----TWILIGHT CASE - Emergent Beam
        XMU1=ABS(U0)
        DO II=J,1,-1
          DIFF=RZ(II)*DSQRT(1.0D0-XMU1**2.)-RZ(II-1)
          IF (DIFF.LT.0.0D0) THEN
            XMU2=DSQRT(1.0D0-(1.0D0-XMU1**2.)/RQ(II-1))
            XL=ABS(RZ(II)*XMU1-RZ(II-1)*XMU2)
            WTAU(II,J)=WTAU(II,J)+XL
            WTAU(II-1,J)=WTAU(II-1,J)+XL
            XMU1=XMU2
          ELSE      ! Lowest level intersected by emergent beam
            XL=RZ(II)*XMU1*2.0D0
            WTING=DIFF/(RZ(II)-RZ(II-1))
            WTAU(II,J)=WTAU(II,J)+XL*0.5D0*(1.0D0+WTING)
            WTAU(II-1,J)=WTAU(II-1,J)+XL*0.5D0*(1.0D0-WTING)
            GOTO 16
          ENDIF
        END DO
   16 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE OPMIE(KW,WAVEL,FMEAN)
!@sum OPMIE Mie code for J's, only uses 8-term expansion, 4-Gauss pts
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
!@calls MIESCT
C ----------------------------------------------------------------------
C  Currently allow only one aerosol phase fn (ALL altitudes) to
C  be associated with extinction AER(1:NC) = aer ext @ 1000 nm
C  Pick Mie-wavelength with phase function and Qext:
c
C  MIEDX:1=Rayly 2=iso 30iso-equiv 4=bkgrd-sulf,5=volc-sulf,6=liq water
C    1   Rayleigh phase scattering-----
C    2   Isotropic scattering----------
C    3   Isotropic equivalent: take 0.2*Mie optical depth (at 1000 nm)
C    4   Mie phase for std background strat particles (r=0.09,sig=.6)
C    5   Mie phase fn for Volcanic (large) particles (r=0.08,sig=.8)
C    6   Mie phase fn for LARGE water drops (10 micron=r, gamma=0.150)
C----------------------------------------------------------------------
C  FUNCTION RAYLAY(WAVE)---RAYLEIGH CROSS-SECTION for wave > 170 nm
C       WSQI = 1.E6/(WAVE*WAVE)
C       REFRM1 = 1.0E-6*(64.328+29498.1/(146.-WSQI)+255.4/(41.-WSQI))
C       RAYLAY = 5.40E-21*(REFRM1*WSQI)**2
C---------------------------------------------------------------------
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: LM
      USE TRCHEM_Shindell_COM, only:XQO2,XQO3,N__,M__,MIEDX,QAAFASTJ,
     &                          ncFASTJ,jaddto,jaddlv,FTAU,NLBATM,XLTAU,
     &                          DO3,DMFASTJ,DBC,QBC,QRAYL,AER,DTAUDZ,
     &                          PIRAY,PIAER,ZZHT,TTAU,dtaumax,NLFASTJ,
     &                          RFLECT,MFIT,PAA,FJFASTJ,dpomega,NAA,
     &                          POMEGA,ZTAU,FZ,POMEGAJ,jndlv,jndlev,
     &                          ZREFL,ZU0,ZFASTJ,ZFLUX,WTAU,U0
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var KW passed wavelength bin index
!@var FMEAN    Mean actinic flux at desired levels
!@var WAVEL Effective wavelength of each wavelength bin
!@var i,j,k,ix dummy loop variables
!@var J1 local copy of NLBATM
!@var KM,QXMIE ?
!@var xlo2,xlo3,xlray,xlbc,xlaer {o2,o3,rayleigh,black carbon,aerosol}
!@+   scattering
!@var taudif,taudn,tauup ?
!@var zk        fractional increment in level
!@var dttau     change in ttau per increment    (linear, positive)
!@var ftaulog   change in ftau per increment    (exp, normally < 1)
!@var ND is 2*NCFASTJ + 3 - 2*J1

      REAL*8 QXMIE,WAVEL
      INTEGER, INTENT(IN)     :: KW
      REAL*8, DIMENSION(LM)   :: FMEAN
      INTEGER KM,j,i,k,ix,J1,ND
      REAL*8 xlo2,xlo3,xlray,xlbc,xlaer,
     &  zk,dttau,ftaulog,taudif,taudn,tauup
c
C---Pick nearest Mie wavelength, no interpolation--------------
                              KM=1
      IF( WAVEL .GT. 355.d0 ) KM=2
      IF( WAVEL .GT. 500.d0 ) KM=3
C     IF( WAVEL .GT. 800.d0 ) KM=4  !drop the 1000 nm wavelength
      MIEDX = MIN(NAA, MAX(1, MIEDX))
      QXMIE = QAAFASTJ(KM,MIEDX)/QAAFASTJ(4,MIEDX)
C---Reinitialize level arrays
      do j=1,ncFASTJ
        jaddlv(j)=0
        jaddto(j)=0
        ftau(j)=0.d0             ! Only start of array, note
      enddo
      do j=1,LM
        jndlv(j)=0
      enddo
C
C---Setup atmosphere, TTAU(1) is column optical depth at surface, TTAU(N
C   +1)=0
      J1 = NLBATM
      DO J=J1,NCFASTJ
        XLO3=DO3(J)*XQO3(J)
        XLO2=DMFASTJ(J)*XQO2(J)*0.20948d0
        XLBC=DBC(J)*QBC(KW)
        XLRAY=DMFASTJ(J)*QRAYL(KW)
C---For Mie code scale extinction at 1000 nm to wavelength WAVEL (QXMIE)
        XLAER=AER(J)*QXMIE
c
        DTAUDZ(J)=XLO3+XLO2+XLBC+XLRAY+XLAER
c
        PIRAY(J)=XLRAY/DTAUDZ(J)
        PIAER(J)=XLAER/DTAUDZ(J)
      ENDDO
      TTAU(NCFASTJ+1)=0.0D0
      TTAU(NCFASTJ)=DTAUDZ(NCFASTJ)*ZZHT
      jaddlv(NCFASTJ)=int(ttau(NCFASTJ)/dtaumax)
      DO J=NCFASTJ-1,J1,-1
        TAUDIF=(ZFASTJ(J+1)-ZFASTJ(J))*(DTAUDZ(J)+DTAUDZ(J+1))*0.5D0
        TTAU(J)=TTAU(J+1) + TAUDIF
        jaddlv(j)=int(taudif/dtaumax)
      ENDDO
c
c Calculate cumulative total and define levels we want J-values at.
c Sum upwards for levels, and then downwards for Mie code readjustments.
c
      jaddto(J1)=jaddlv(J1)
      do j=J1+1,ncFASTJ
        jaddto(j)=jaddlv(j)+jaddto(j-1)
      enddo
      if((jaddto(ncFASTJ)+ncFASTJ).gt.nlfastj) then
         write(6,*)' Too many levels in photolysis code: need'
         write(6,*) jaddto(ncFASTJ)+ncFASTJ
         WRITE(6,*)'but NLFASTJ dimensioned as ',NLFASTJ
         call stop_model('Too many levels in photolysis code.',255)
      endif
      do i=LM,1,-1
        if(jndlev(i).ne.1) jndlv(i)=jndlev(i)+jaddto(jndlev(i)-1)
      enddo
      jaddto(ncFASTJ)=jaddlv(ncFASTJ)
      do j=ncFASTJ-1,J1,-1
        jaddto(j)=jaddlv(j)+jaddto(j+1)
      enddo
c
C Calculate attenuated incident beam EXP(-TTAU/U0) and flux on surface
      DO J=J1,NCFASTJ
        IF (WTAU(J,J).GT.0.0D0) THEN
          XLTAU=0.0D0
          DO I=1,NCFASTJ
            XLTAU=XLTAU + DTAUDZ(I)*WTAU(I,J)
          ENDDO
          FTAU(J)=DEXP(-XLTAU)
        ELSE
          FTAU(J)=0.0D0
        ENDIF
      ENDDO
      FTAU(NCFASTJ+1)=1.D0
      IF (U0.GT.0.D0) THEN
        ZFLUX = U0*FTAU(J1)*RFLECT/(1.d0+RFLECT)
      ELSE
        ZFLUX = 0.d0
      ENDIF
C
C---------------------SET UP FOR MIE CODE-------------------------------
c
c  Transpose the ascending TTAU grid to a descending ZTAU grid.
c  Double the resolution - TTAU points become the odd points on the
c  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
c  Odd point added at top of grid for unattenuated beam   (Z='inf')
c
c        Surface:   TTAU(1)   now use ZTAU(2*NCFASTJ+1)
c        Top:       TTAU(NCFASTJ)  now use ZTAU(3)
c        Infinity:            now use ZTAU(1)
c
c  Mie scattering code only used from surface to level NCFASTJ
C-----------------------------------------------------------------
C
      ND = 2*NCFASTJ + 3 - 2*J1
c
C-------CLEAR THE POMEGA'S
      DO K=1,N__
        DO I=1,MFIT
          POMEGA(I,K) = 0.0D0
        ENDDO
      ENDDO
C
c Ascend through atmosphere transposing grid and adding extra points
c Define scattering phase function with mix of Rayleigh(1) & Mie(MIEDX)
c
c---Set up pomega for needed levels first
      do j=j1+1,ncFASTJ,2
        do i=1,MFIT
          pomegaj(i,j) = PIAER(J)*PAA(I,KM,MIEDX)+PIRAY(J)*PAA(I,KM,1)
        enddo
      enddo
c
c---Set up pomega for boundary levels
      do i=1,MFIT
        pomegaj(i,J1) = PIAER(J1)*PAA(I,KM,MIEDX)+PIRAY(J1)*PAA(I,KM,1)
      enddo
      do j=J1+2,ncFASTJ,2
        taudn = ttau(j-1)-ttau(j)
        tauup = ttau(j)-ttau(j+1)
        do i=1,MFIT
          pomegaj(i,j) = (pomegaj(i,j-1)*taudn +
     $                    pomegaj(i,j+1)*tauup) / (taudn+tauup)
        enddo
      enddo
      do i=1,MFIT
        pomegaj(i,ncFASTJ+1) = pomegaj(i,ncFASTJ)
      enddo
c
c Transpose profiles for standard model levels first
      do j=J1,ncFASTJ
        k = 2*(ncFASTJ+1-j)+2*jaddto(j)+1
        ztau(k)= ttau(j)
        fz(k)  = ftau(j)
        do i=1,MFIT
          pomega(i,k) = pomegaj(i,j)
        enddo
      enddo
c
c Add top-of-atmosphere point
      ztau(1) = 0.d0
      fz(1) = 1.d0
      do i=1,MFIT
        pomega(i,1) = pomega(i,3)
      enddo
c
C---------------------------------------------------------------------
c  Insert new levels, working downwards from the top of the atmosphere
c  to the surface (down in 'j', up in 'k'). This allows ztau and pomega
c  to be incremented linearly (in a +ve sense), and the flux fz to be
c  attenuated top-down (avoiding problems where lower level fluxes are
c  zero).
c
c    zk        fractional increment in level
c    dttau     change in ttau per increment    (linear, positive)
c    dpomega   change in pomega per increment  (linear)
c    ftaulog   change in ftau per increment    (exp, normally < 1)
c
C---------------------------------------------------------------------
      do j=ncFASTJ,J1,-1
        if(jaddlv(j).gt.0) then
          zk = 1.d0/(1.d0+dble(jaddlv(j)))
          dttau = (ttau(j)-ttau(j+1))*zk
          do i=1,MFIT
            dpomega(i) = (pomegaj(i,j)-pomegaj(i,j+1))*zk
          enddo
          if(ftau(j).eq.0.d0.or.ftau(j+1).eq.0.d0) then
            ftaulog=1.0d-05
          else
            ftaulog = ftau(j)/ftau(j+1)
            if(ftaulog.gt.1.d+150.or.ftaulog.lt.1.d-150) then
cc              ftaulog=exp(log(1.0d-150)*zk)
              ftaulog=1.0d-05
            else
              ftaulog=exp(log(ftaulog)*zk)
            endif
          endif
          k = 2*(ncFASTJ-j+jaddto(j)-jaddlv(j))+1   !  k at level j+1
c
          do ix=1,jaddlv(j)
            ztau(k+2) = ztau(k) + dttau
            fz(k+2) = fz(k)*ftaulog
            do i=1,MFIT
              pomega(i,k+2) = pomega(i,k) + dpomega(i)
            enddo
            k = k+2
          enddo
        endif
      enddo

c
C---Update total number of levels
      ncFASTJ=ncFASTJ+jaddto(J1)
      ND = 2*NCFASTJ + 3 - 2*J1
      if(nd.gt.N__) then
         write(6,'(a,a,i3,a,i3)')' Too many levels in photolysis code:',
     $            ' need ',nd,' but N__ dimensioned as ',N__
         call stop_model('Too many levels in photolysis code.',255)
      endif
c
C---Add boundary/ground layer to ensure no negative J's caused by
C---too large a TTAU-step in the 2nd-order lower b.c.
      ZTAU(ND+2) = 0.01d0 + ZTAU(ND)
      FZ(ND+2) = FZ(ND)
      DO I=1,MFIT
       POMEGA(I,ND+2)   = POMEGA(I,ND)
      ENDDO
      ND = ND+2
c
C---Fill in even points K
      DO K=2,ND-1,2
        ZTAU(K) = 0.5d0*(ZTAU(K-1)+ZTAU(K+1))
        FZ(K) = 0.5d0*(FZ(K-1)+FZ(K+1))
        DO I=1,MFIT
          POMEGA(I,K) = 0.5d0*(POMEGA(I,K-1)+POMEGA(I,K+1))
        ENDDO
      ENDDO
      ZU0 = U0
      ZREFL = RFLECT
c
c
C-----------------------------------------
      CALL MIESCT(ND)
C-----------------------------------------
c
c Accumulate attenuation for selected levels
      DO j=1,LM
        k=2*NCFASTJ+3-2*jndlv(j)
        if(k.gt.ND-2) then
          FMEAN(j) = 0.d0
        else
          FMEAN(j) = FJFASTJ(k)
        endif
      ENDDO
c
c Reset level parameters
      NCFASTJ = NCFASTJ - jaddto(J1)
c
      RETURN
 1000 format(1x,i3,3(2x,1pe10.4),1x,i3)
 1100 format(1x,a3,4(a9,2x))
 1200 format(1x,i3,11(1x,1pe9.3))
 1300 format(1x,50(i3))
      END SUBROUTINE OPMIE
c
C
       SUBROUTINE MIESCT(ND)
!@sum MIESCT radiative transfer code
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
!@calls LEGND0,BLKSLV
C-------------------------------------------------------------------
C   This is an adaption of the Prather rad transfer code, (mjp, 10/95)
C     Prather, 1974, Astrophys. J. 192, 787-792.
C         Sol'n of inhomogeneous Rayleigh scattering atmosphere.
C         (original Rayleigh w/ polarization)
C     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
C         Raman scattering in the atmospheres of the major planets.
C         (first use of anisotropic code)
C     Jacob, Gottlieb and Prather,89, J.G..Res., 94, 12975-13002.
C         Chemistry of a polluted cloudy boundary layer,
C         (documentation of extension to anisotropic scattering)
C
C    takes atmospheric structure and source terms from std J-code
C    ALSO limited to 4 Gauss points, only calculates mean field!
C
C   mean rad. field ONLY (M=1)
C   initialize variables FIXED/UNUSED in this special version:
C   FTOP=1.0=astrophys flux (unit of pi) at SZA, -ZU0, use for scaling
C   FBOT=0.0=ext isotropic flux on lower boundary
C   SISOTP=0.0=Spec Intensity of isotropic radiation incid from top
C
C   SUBROUTINES:  MIESCT              needs 'jv_mie.cmn'
C                 BLKSLV              needs 'jv_mie.cmn'
C                 GEN (ID)            needs 'jv_mie.cmn'
C                 LEGND0 (X,PL,N)
C                 MATIN4 (AFASTJ)
C-------------------------- --------------------------------------------
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only:NFASTJ,EMU,ZFLUX,ZU0,ZREFL,WTFASTJ,
     &                          MFASTJ,MFIT,PM0,PM,CMEQ1,FJFASTJ,FZ
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var i,im,id dummy loop variables
!@var see OPMIE ?
      INTEGER i,im,id,ND
c
C-------------------------- --------------------------------------------
C---fix scattering to 4 Gausss pts = 8-stream
C---solve eqn of R.T. only for first-order M=1
C      ZFLUX = (ZU0*FZ(ND)*ZREFL+FBOT)/(1.0d0+ZREFL)
      ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      DO I=1,NFASTJ
        CALL LEGND0 (EMU(I),PM0)
        DO IM=MFASTJ,MFIT
          PM(I,IM) = PM0(IM)
        ENDDO
      ENDDO
C
      CALL LEGND0 (-ZU0,PM0)
      DO IM=MFASTJ,MFIT
        PM0(IM) = CMEQ1*PM0(IM)
      ENDDO
C
      CALL BLKSLV(ND)
C
      DO ID=1,ND,2
        FJFASTJ(ID) = 4.0d0*FJFASTJ(ID) + FZ(ID)
      ENDDO
      RETURN
      END SUBROUTINE MIESCT
c
c
      SUBROUTINE BLKSLV(ND)
!@sum BLKSLV Solves the block tri-diagonal system:
!@+   A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
!@calls GEN,MATIN4
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: NFASTJ,RR2,BFASTJ,CC,DD,HFASTJ,
     &                          AFASTJ,C1,AAFASTJ,FJFASTJ,WTFASTJ,EMU
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var i,j,k,id dummy loop variables
!@var sum local variable for summing
!@var ND ? see OPMIE
      integer i, j, k, id, ND
      real*8  sum
c
C-----------UPPER BOUNDARY ID=1
      CALL GEN(1,ND)
      CALL MATIN4(BFASTJ)
      DO I=1,NFASTJ
        RR2(I,1) = 0.0d0
        DO J=1,NFASTJ
         SUM = 0.0d0
         DO K=1,NFASTJ
          SUM = SUM - BFASTJ(I,K)*CC(K,J)
         ENDDO
         DD(I,J,1) = SUM
         RR2(I,1) = RR2(I,1) + BFASTJ(I,J)*HFASTJ(J)
       ENDDO
      ENDDO
C----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
      DO ID=2,ND-1
        CALL GEN(ID,ND)
        DO I=1,NFASTJ
          DO J=1,NFASTJ
          BFASTJ(I,J) = BFASTJ(I,J) + AFASTJ(I)*DD(I,J,ID-1)
          ENDDO
          HFASTJ(I) = HFASTJ(I) - AFASTJ(I)*RR2(I,ID-1)
        ENDDO
        CALL MATIN4 (BFASTJ)
        DO I=1,NFASTJ
          RR2(I,ID) = 0.0d0
          DO J=1,NFASTJ
          RR2(I,ID) = RR2(I,ID) + BFASTJ(I,J)*HFASTJ(J)
          DD(I,J,ID) = - BFASTJ(I,J)*C1(J)
          ENDDO
        ENDDO
      ENDDO
C---------FINAL DEPTH POINT: ND
      CALL GEN(ND,ND)
      DO I=1,NFASTJ
        DO J=1,NFASTJ
          SUM = 0.0d0
          DO K=1,NFASTJ
          SUM = SUM + AAFASTJ(I,K)*DD(K,J,ND-1)
          ENDDO
        BFASTJ(I,J) = BFASTJ(I,J) + SUM
        HFASTJ(I) = HFASTJ(I) - AAFASTJ(I,J)*RR2(J,ND-1)
        ENDDO
      ENDDO
      CALL MATIN4 (BFASTJ)
      DO I=1,NFASTJ
        RR2(I,ND) = 0.0d0
        DO J=1,NFASTJ
        RR2(I,ND) = RR2(I,ND) + BFASTJ(I,J)*HFASTJ(J)
        ENDDO
      ENDDO
C-----------BACK SOLUTION
      DO ID=ND-1,1,-1
       DO I=1,NFASTJ
        DO J=1,NFASTJ
         RR2(I,ID) = RR2(I,ID) + DD(I,J,ID)*RR2(J,ID+1)
        ENDDO
       ENDDO
      ENDDO
C----------MEAN J & H
      DO ID=1,ND,2
       FJFASTJ(ID) = 0.0d0
       DO I=1,NFASTJ
        FJFASTJ(ID) = FJFASTJ(ID) + RR2(I,ID)*WTFASTJ(I)
       ENDDO
      ENDDO
      DO ID=2,ND,2
       FJFASTJ(ID) = 0.0d0
       DO I=1,NFASTJ
        FJFASTJ(ID) = FJFASTJ(ID) + RR2(I,ID)*WTFASTJ(I)*EMU(I)
       ENDDO
      ENDDO
      RETURN
      END SUBROUTINE BLKSLV
c
c
      SUBROUTINE GEN(ID,ND)
!@sum GEN Generate coefficient matrices for block tri-diagonal system:
!@ +  A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
!@calls
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: NFASTJ,MFASTJ,MFIT,POMEGA,PM,PM0,
     &                        HFASTJ,AFASTJ,FZ,WTFASTJ,SFASTJ,WFASTJ,U1,
     &                        EMU,C1,V1,ZTAU,CC,ZFLUX,AAFASTJ,BFASTJ,
     &                        ZREFL
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var id passed
!@var ND ? see OPMIE
!@var id0,id1,im,i,j,k,mstart dummy variables
!@var sum0,sum1,sum2,sum3,deltau,d1,d2,surfac dummy variables
      integer id, id0, id1, im, i, j, k, mstart, ND
      real*8  sum0, sum1, sum2, sum3, deltau, d1, d2, surfac, ttmp(8)
c
C---------------------------------------------
      IF(ID.EQ.1 .OR. ID.EQ.ND) THEN
C---------calculate generic 2nd-order terms for boundaries
       ID0 = ID
       ID1 = ID+1
       IF(ID.GE.ND) ID1 = ID-1
       DO 10 I=1,NFASTJ
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
        DO IM=MFASTJ,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        ENDDO
        DO IM=MFASTJ+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        ENDDO
         HFASTJ(I) = 0.5d0*(SUM0*FZ(ID0) + SUM2*FZ(ID1))
         AFASTJ(I) = 0.5d0*(SUM1*FZ(ID0) + SUM3*FZ(ID1))
        DO J=1,I
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
         DO IM=MFASTJ,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         ENDDO
         DO IM=MFASTJ+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         ENDDO
         SFASTJ(I,J) = - SUM2*WTFASTJ(J)
         SFASTJ(J,I) = - SUM2*WTFASTJ(I)
         WFASTJ(I,J) = - SUM1*WTFASTJ(J)
         WFASTJ(J,I) = - SUM1*WTFASTJ(I)
         U1(I,J) = - SUM3*WTFASTJ(J)
         U1(J,I) = - SUM3*WTFASTJ(I)
          SUM0 = 0.5d0*(SUM0 + SUM2)
         BFASTJ(I,J) = - SUM0*WTFASTJ(J)
         BFASTJ(J,I) = - SUM0*WTFASTJ(I)
        ENDDO
         SFASTJ(I,I) = SFASTJ(I,I) + 1.0d0
         WFASTJ(I,I) = WFASTJ(I,I) + 1.0d0
         U1(I,I) = U1(I,I) + 1.0d0
         BFASTJ(I,I) = BFASTJ(I,I) + 1.0d0
   10  CONTINUE
       DO I=1,NFASTJ
         SUM0 = 0.0d0
        DO J=1,NFASTJ
         SUM0 = SUM0 + SFASTJ(I,J)*AFASTJ(J)/EMU(J)
        ENDDO
        C1(I) = SUM0
       ENDDO
       DO I=1,NFASTJ
        DO J=1,NFASTJ
          SUM0 = 0.0d0
          SUM2 = 0.0d0
         DO K=1,NFASTJ
          SUM0 = SUM0 + SFASTJ(J,K)*WFASTJ(K,I)/EMU(K)
          SUM2 = SUM2 + SFASTJ(J,K)*U1(K,I)/EMU(K)
         ENDDO
         AFASTJ(J) = SUM0
         V1(J) = SUM2
        ENDDO
        DO J=1,NFASTJ
         WFASTJ(J,I) = AFASTJ(J)
         U1(J,I) = V1(J)
        ENDDO
       ENDDO
       IF (ID.EQ.1) THEN
C-------------upper boundary, 2nd-order, C-matrix is full (CC)
        DELTAU = ZTAU(2) - ZTAU(1)
        D2 = 0.25d0*DELTAU
        DO I=1,NFASTJ
          D1 = EMU(I)/DELTAU
          DO J=1,NFASTJ
           BFASTJ(I,J) = BFASTJ(I,J) + D2*WFASTJ(I,J)
           CC(I,J) = D2*U1(I,J)
          ENDDO
          BFASTJ(I,I) = BFASTJ(I,I) + D1
          CC(I,I) = CC(I,I) - D1
C         HFASTJ(I) = HFASTJ(I) + 2.0d0*D2*C1(I) + D1*SISOTP
          HFASTJ(I) = HFASTJ(I) + 2.0d0*D2*C1(I)
          AFASTJ(I) = 0.0d0
        ENDDO
       ELSE
C-------------lower boundary, 2nd-order, A-matrix is full (AAFASTJ)
        DELTAU = ZTAU(ND) - ZTAU(ND-1)
        D2 = 0.25d0*DELTAU
        SURFAC = 4.0d0*ZREFL/(1.0d0 + ZREFL)
        DO I=1,NFASTJ
          D1 = EMU(I)/DELTAU
          HFASTJ(I) = HFASTJ(I) - 2.0d0*D2*C1(I)
           SUM0 = 0.0d0
          DO J=1,NFASTJ
           SUM0 = SUM0 + WFASTJ(I,J)
          ENDDO
           SUM0 = D1 + D2*SUM0
           SUM1 = SURFAC*SUM0
          DO J=1,NFASTJ
          BFASTJ(I,J)=BFASTJ(I,J)+D2*WFASTJ(I,J)-SUM1*EMU(J)*WTFASTJ(J)
          ENDDO
          BFASTJ(I,I) = BFASTJ(I,I) + D1
          HFASTJ(I) = HFASTJ(I) + SUM0*ZFLUX
          DO J=1,NFASTJ
           AAFASTJ(I,J) = - D2*U1(I,J)
          ENDDO
           AAFASTJ(I,I) = AAFASTJ(I,I) + D1
           C1(I) = 0.0d0
        ENDDO
       ENDIF
C------------intermediate points:  can be even or odd, A & C diagonal
      ELSE
        DELTAU = ZTAU(ID+1) - ZTAU(ID-1)
        MSTART = MFASTJ + MOD(ID+1,2)
        DO I=1,NFASTJ
          AFASTJ(I) = EMU(I)/DELTAU
          C1(I) = -AFASTJ(I)
           SUM0 = 0.0d0
          DO IM=MSTART,MFIT,2
           TTMP(IM) = POMEGA(IM,ID)*PM(I,IM)
           SUM0 = SUM0 + TTMP(IM)*PM0(IM)
c           SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM0(IM)
          ENDDO
          HFASTJ(I) = SUM0*FZ(ID)
          DO J=1,I
            SUM0 = 0.0d0
           DO IM=MSTART,MFIT,2
            SUM0 = SUM0 + TTMP(IM)*PM(J,IM)
c            SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM(J,IM)
           ENDDO
            BFASTJ(I,J) =  - SUM0*WTFASTJ(J)
            BFASTJ(J,I) =  - SUM0*WTFASTJ(I)
          ENDDO
          BFASTJ(I,I) = BFASTJ(I,I) + 1.0d0
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE GEN
c
c
      SUBROUTINE LEGND0(X,PL)
!@sum LEGND0 Calculates ORDINARY LEGENDRE fns of X (real) from:
!@+   P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
!@calls
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: MFIT
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var I dummy loop variable
!@var PL ?
!@var X passed EMU or -ZU0
!@var DEN denominator
      INTEGER I
      REAL*8, DIMENSION(MFIT) :: PL
      REAL*8 X,DEN
C---Always does PL(2) = P[1]
        PL(1) = 1.D0
        PL(2) = X
        DO I=3,MFIT
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.D0/DEN) - PL(I-2)*(1.d0-1.D0/DEN)
        ENDDO
      RETURN
      END SUBROUTINE LEGND0
c
c
      SUBROUTINE MATIN4(AFASTJ)
!@sum MATIN4 invert 4x4 matrix A(4,4) in place with L-U decomp
!@+   (mjp, old...)
!@auth Lee Grenfell/Oliver Wild (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3_fastjlg_ozone_M23)
!@calls
c
C**** GLOBAL parameters and variables:
C (none)
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var AFASTJ passed (actually BFASTJ)
      REAL*8 AFASTJ(4,4)
C
C---SETUP L AND U
      AFASTJ(2,1) = AFASTJ(2,1)/AFASTJ(1,1)
      AFASTJ(2,2) = AFASTJ(2,2)-AFASTJ(2,1)*AFASTJ(1,2)
      AFASTJ(2,3) = AFASTJ(2,3)-AFASTJ(2,1)*AFASTJ(1,3)
      AFASTJ(2,4) = AFASTJ(2,4)-AFASTJ(2,1)*AFASTJ(1,4)
      AFASTJ(3,1) = AFASTJ(3,1)/AFASTJ(1,1)
      AFASTJ(3,2) = (AFASTJ(3,2)-AFASTJ(3,1)*AFASTJ(1,2))/AFASTJ(2,2)
      AFASTJ(3,3) = AFASTJ(3,3)-AFASTJ(3,1)*AFASTJ(1,3)-AFASTJ(3,2)
     $*AFASTJ(2,3)
      AFASTJ(3,4) = AFASTJ(3,4)-AFASTJ(3,1)*AFASTJ(1,4)-AFASTJ(3,2)
     $*AFASTJ(2,4)
      AFASTJ(4,1) = AFASTJ(4,1)/AFASTJ(1,1)
      AFASTJ(4,2) = (AFASTJ(4,2)-AFASTJ(4,1)*AFASTJ(1,2))/AFASTJ(2,2)
      AFASTJ(4,3) = (AFASTJ(4,3)-AFASTJ(4,1)*AFASTJ(1,3)-AFASTJ(4,2)
     $*AFASTJ(2,3))/AFASTJ(3,3)
      AFASTJ(4,4) = AFASTJ(4,4)-AFASTJ(4,1)*AFASTJ(1,4)-AFASTJ(4,2)*
     $AFASTJ(2,4)-AFASTJ(4,3)*AFASTJ(3,4)
C---INVERT L
      AFASTJ(4,3) = -AFASTJ(4,3)
      AFASTJ(4,2) = -AFASTJ(4,2)-AFASTJ(4,3)*AFASTJ(3,2)
      AFASTJ(4,1) = -AFASTJ(4,1)-AFASTJ(4,2)*AFASTJ(2,1)-AFASTJ(4,3)
     $*AFASTJ(3,1)
      AFASTJ(3,2) = -AFASTJ(3,2)
      AFASTJ(3,1) = -AFASTJ(3,1)-AFASTJ(3,2)*AFASTJ(2,1)
      AFASTJ(2,1) = -AFASTJ(2,1)
C---INVERT U
      AFASTJ(4,4) = 1.D0/AFASTJ(4,4)
      AFASTJ(3,4) = -AFASTJ(3,4)*AFASTJ(4,4)/AFASTJ(3,3)
      AFASTJ(3,3) = 1.D0/AFASTJ(3,3)
      AFASTJ(2,4) = -(AFASTJ(2,3)*AFASTJ(3,4)+AFASTJ(2,4)*
     $AFASTJ(4,4))/AFASTJ(2,2)
      AFASTJ(2,3) = -AFASTJ(2,3)*AFASTJ(3,3)/AFASTJ(2,2)
      AFASTJ(2,2) = 1.D0/AFASTJ(2,2)
      AFASTJ(1,4) = -(AFASTJ(1,2)*AFASTJ(2,4)+AFASTJ(1,3)*AFASTJ(3,4)
     $+AFASTJ(1,4)*AFASTJ(4,4))/AFASTJ(1,1)
      AFASTJ(1,3) = -(AFASTJ(1,2)*AFASTJ(2,3)+AFASTJ(1,3)*
     $AFASTJ(3,3))/AFASTJ(1,1)
      AFASTJ(1,2) = -AFASTJ(1,2)*AFASTJ(2,2)/AFASTJ(1,1)
      AFASTJ(1,1) = 1.D0/AFASTJ(1,1)
C---MULTIPLY (U-INVERSE)*(L-INVERSE)
      AFASTJ(1,1) = AFASTJ(1,1)+AFASTJ(1,2)*AFASTJ(2,1)+AFASTJ(1,3)
     $*AFASTJ(3,1)+AFASTJ(1,4)*AFASTJ(4,1)
      AFASTJ(1,2) = AFASTJ(1,2)+AFASTJ(1,3)*AFASTJ(3,2)+AFASTJ(1,4)
     $*AFASTJ(4,2)
      AFASTJ(1,3) = AFASTJ(1,3)+AFASTJ(1,4)*AFASTJ(4,3)
      AFASTJ(2,1) = AFASTJ(2,2)*AFASTJ(2,1)+AFASTJ(2,3)*AFASTJ(3,1)
     $+AFASTJ(2,4)*AFASTJ(4,1)
      AFASTJ(2,2) = AFASTJ(2,2)+AFASTJ(2,3)*AFASTJ(3,2)+AFASTJ(2,4)
     $*AFASTJ(4,2)
      AFASTJ(2,3) = AFASTJ(2,3)+AFASTJ(2,4)*AFASTJ(4,3)
      AFASTJ(3,1) = AFASTJ(3,3)*AFASTJ(3,1)+AFASTJ(3,4)*AFASTJ(4,1)
      AFASTJ(3,2) = AFASTJ(3,3)*AFASTJ(3,2)+AFASTJ(3,4)*AFASTJ(4,2)
      AFASTJ(3,3) = AFASTJ(3,3)+AFASTJ(3,4)*AFASTJ(4,3)
      AFASTJ(4,1) = AFASTJ(4,4)*AFASTJ(4,1)
      AFASTJ(4,2) = AFASTJ(4,4)*AFASTJ(4,2)
      AFASTJ(4,3) = AFASTJ(4,4)*AFASTJ(4,3)
      RETURN
      END SUBROUTINE MATIN4
c
