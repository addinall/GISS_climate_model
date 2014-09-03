#include "rundeck_opts.h"

      subroutine OCN_mesosc

      USE MODEL_COM,  only : nstep=>itime,itimei
     .                    ,JMON,jhour,nday,jdate,jday
     . ,iyear1,jdendofm,jyear,aMON,dtsrc
     . ,xlabel,lrunid
      USE CONSTANT,   only : grav,omega,sday
      USE OCEANR_DIM, only : ogrid
      USE OCEANRES,   only : idm=>imo,jdm=>jmo,kdm=>lmo,dzo
      USE OFLUXES,    only : oRSI,oAPRESS
      USE OCEAN,      only : ZOE=>ZE,g0m,s0m,mo,dxypo,focean,lmm
     .                      ,oLON_DG,oLAT_DG,uo,vo,sinpo
      USE KPP_COM,    only : kpl

      USE GM_COM, only: RHOX, RHOY
      USE ODIAG, only: oij=>oij_loc,oijl=>oijl_loc
     .            ,ij_eke,ij_rd,ijl_ueddy,ijl_veddy,ijl_n2

      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      use TimerPackage_mod


      implicit none

      integer, parameter :: itest=97,jtest=128
      integer i,j,k,l,ndepi
      integer i_0,i_1,j_0,j_1

      Real*8,External   :: VOLGSP,TEMGSP,TEMGS
      real*8  g,s,pres
      real*8   p1d(kdm+1),dp1d(kdm),temp1d(kdm),saln1d(kdm)
      real*8   rho_water,amld_cgs,Rd
      real*8   K0,Ustar(kdm),Ustar_star(kdm)
      real*8   Vstar(kdm),Vstar_star(kdm)
      real*8   z_cm(kdm+1),dens_cgs(kdm)
     .        ,uvel_cgs(kdm),vvel_cgs(kdm)
     .        ,drhodz_cgs(kdm),coriol,n2(kdm)
     .        ,drhodx_cgs(kdm),drhody_cgs(kdm)

      logical vrbos

      real*8 map(288,180)
 

      if (nstep.eq.0) return

      i_0=ogrid%I_STRT
      i_1=ogrid%I_STOP
      j_0=ogrid%J_STRT
      j_1=ogrid%J_STOP


      do 1000 j=j_0,j_1
      do 1000 i=i_0,i_1
      IF(FOCEAN(I,J).gt.0.) THEN

      vrbos=.false.
c     vrbos=.true.
      if (i.eq.itest.and.j.eq.jtest) vrbos=.true.

! max depth cell
      ndepi=lmm(i,j)

!change units
       pres = oAPRESS(i,j)    !surface atm. pressure
       do k=1,lmm(i,j)
         pres=pres+MO(I,J,k)*GRAV*.5
         g=G0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
         s=S0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
         temp1d(k)=TEMGSP(g,s,pres)    !in situ   temperature
         saln1d(k)=s*1000.             !convert to psu (eg. ocean mean salinity=35psu)
         rho_water = 1d0/VOLGSP(g,s,pres)
         dp1d(k)=MO(I,J,K)/rho_water   !local thickenss of each layer in meters
         !add missing part of density to get to the bottom of the layer
         !now pres is at the bottom of the layer
         pres=pres+MO(I,J,k)*GRAV*.5
      enddo

      p1d(1)=0.
      do k=2,kdm+1
      p1d(k)=p1d(k-1)+dp1d(k-1)
      enddo
      do k=1,kdm
       z_cm(k) = p1d(k)*100.d0                !depth in cm
       dens_cgs(k) = MO(i,j,k)/max(1.d0,dp1d(k)) *1.d-3 
       uvel_cgs(k) = uo(i,j,k) * 100.d0
       vvel_cgs(k) = vo(i,j,k) * 100.d0
      enddo

      z_cm(kdm+1) = p1d(kdm+1)*100.d0                !depth in cm
      coriol = 2d0*omega*sinpo(j)      !evaluated at tracer points
!drhodz_cgs
      do k=1,kdm-1
         if (dens_cgs(k+1).ne.0.d0.and.dens_cgs(k).ne.0.d0) then
           drhodz_cgs(k) = (dens_cgs(k+1) - dens_cgs(k))          ! at (i,j,k+1/2)   interface
     .                 /(max(1.d0,z_cm(k+1)+z_cm(k)))/2.
         elseif (dens_cgs(k+1).eq.0.d0.
     *       and.dens_cgs(k).ne.0.d0.
     *       and.dens_cgs(max(1,k-1)).ne.0.d0) then
           drhodz_cgs(k)=drhodz_cgs(max(1,k-1))
         else
           drhodz_cgs(k) = 1.d0
         endif
      enddo
      do k=1,kdm
         drhodx_cgs(k)=rhox(i,j,k)*1.d-3/100.d0
         drhody_cgs(k)=rhoy(i,j,k)*1.d-3/100.d0
      enddo
! compute Brunt-Vaisala
      do k=1,kdm-1
        if (drhodz_cgs(k).ne.1d30) then
        n2(k) = GRAV* drhodz_cgs(k)*100.d0     !no minus sign here, because drhodz defined as rho(k+1)-rho(k)
     .            /(max(1.d0,dens_cgs(k+1)+dens_cgs(k))/2.)   !at (i,j,k+1/2)   interface
        else
        n2(k) = 1.d30
        endif
      enddo
! mixed layer depth
      amld_cgs = 0.d0
      do k=1,kpl(i,j)
      if (dp1d(k).ne.1.d30) then
          amld_cgs = amld_cgs + dp1d(k)*100.d0      ! at (i,j,k) the middle of last layer in MLD
      endif
      enddo
      if (amld_cgs.gt.p1d(ndepi)*100.d0) amld_cgs=p1d(ndepi)*100.d0

      call mesoscales1d(kdm,ndepi,z_cm,
     .      dens_cgs,uvel_cgs,vvel_cgs,
     .      n2,drhodx_cgs,drhody_cgs,drhodz_cgs,coriol,amld_cgs
     .     ,Rd,K0,Ustar,Vstar,Ustar_star,Vstar_star,i,j)

!     if (vrbos) then
!     write(*,'(a,5i5,14e12.4)')'MESOSCALES1:',
!    .   nstep,i,j,k,ndepi,z_cm(k),dens_cgs(k),uvel_cgs(k)
!    .  ,vvel_cgs(k),n2(k),drhodx_cgs(k),drhody_cgs(k)
!    .  ,drhodz_cgs(k),coriol,amld_cgs,Rd,K0,ustar(k),vstar(k)
!     endif
c     if (nstep.eq.1) then
c     if (vrbos) then
c     do k=1,kdm
c     write(*,'(a,5i5,14e12.4)')'MESOSCALES:',
c    .   nstep,i,j,k,ndepi,z_cm(k),dens_cgs(k),uvel_cgs(k)
c    .  ,vvel_cgs(k),n2(k),drhodx_cgs(k),drhody_cgs(k)
c    .  ,drhodz_cgs(k),coriol,amld_cgs,Rd,K0,ustar(k),vstar(k)
c     enddo
c     endif
c     endif

       OIJ(I,J,IJ_eke)  = OIJ(I,J,IJ_eke) + K0      ! eddy kinetic energy, ocean
       OIJ(I,J,IJ_rd )  = OIJ(I,J,IJ_rd ) + Rd      ! Rossby radius of deformation
       DO k=1,kdm
          OIJL(I,J,k,IJL_n2   )= OIJL(I,J,k,IJL_n2   ) + n2(k)    ! brunt vaisala squared
          OIJL(I,J,k,IJL_ueddy)= OIJL(I,J,k,IJL_ueddy) + ustar(k) ! ustar, eddy induced velocity (Canuto)
          OIJL(I,J,k,IJL_veddy)= OIJL(I,J,k,IJL_veddy) + vstar(k) ! vstar, eddy induced velocity (Canuto)
       ENDDO

c     map(i,j)=map(i,j)+K0

      endif    !focean

 1000 continue

c     IF(nstep.eq.47) THEN
c     WRITE(*,*)'output K0 map'
c     OPEN(1010,FILE='map-k0.dat')
c     DO i=1,288
c       DO j=1,180
c         WRITE(1010,9050)i,j,map(i,j)
c       ENDDO
c     ENDDO
c9050 FORMAT(2I4,1(1pe12.3))
c     ENDIF


      end subroutine OCN_mesosc

      subroutine mesoscales1d(km,ndepi,z_cm,
     .      rhoi,UI,VI,n2,dxrho,dyrho,dzrho,f,ml
     .     ,rdm,K02,Ustar,Vstar,Ustar_star,Vstar_star,igrid,jgrid) 

        IMPLICIT NONE

       real*8, parameter :: huge=1.d30
        INTEGER k,km
        INTEGER i0,j0,kmli,kmli5,kn2i
        REAL*8 z_cm(km)
C
        integer igrid,jgrid
        integer ndepi
        REAL*8 Ustar(km),Ustar_star(km)
        REAL*8 Vstar(km),Vstar_star(km)
        REAL*8 KAPPAM
        REAL*8 SIGMAT,f,ml,yt,rd,rdm,pi
        REAL*8 n2m,
     *    rrho(km)
        REAL*8,ALLOCATABLE :: LXsm(:),LYsm(:),DZLXsm(:),DZLYsm(:)
        REAL*8 uve,vve,tem,sal
c       PARAMETER(SIGMAT=0.72)
        PARAMETER(SIGMAT=1.)
c       PARAMETER(SIGMAT=4.)
        INTEGER kkpmax,kkpmaxp5,kkpmax2,kkpmax2p5,kQmax,kkpint,kdbmax,
     *    kn2max,kkpmax2p
        REAL*8 DZLXB1AVE,DZLYB1AVE,kpmax,zkpmax,kpmax2,Qmax,dbmax,n2max,
     *    kpmax2p
        REAL*8 UI(km),VI(km)
        REAL*8 UB1AVE,VB1AVE,UL,VL
        REAL*8 rhoi(km),dbmaptop,dbmapml,
     *    dbmaphalf,map
        REAL*8 frd,n2(km),n2a(km),lr
        REAL*8 UD,VD,UD1(km),VD1(km)
        REAL*8 K0,K02,Qmap,
     *    PHIM,CK,K03,Cs
        REAL*8 dxrho(km),dyrho(km),dzrho(km)
        REAL*8 DAML,DBML
        REAL*8, ALLOCATABLE :: n2t(:),n2t1(:),z(:),zm(:),b1(:),b1t(:),
     *   U(:),V(:),rho(:),n2tr(:),
     *   LX(:),LY(:),DZLX(:),DZLY(:),
     *   UINT(:),VINT(:),UT(:),VT(:),
     *   M1X(:),M1Y(:),M2X(:),M2Y(:),
     *   dxb(:),dyb(:),KP(:),KINT(:),db(:),
     *   zint(:),UREV(:),VREV(:),KPINT(:),
     *   b1rev(:),DZLXREV(:),DZLYREV(:),
     *   UTINT(:),VTINT(:),DZU(:),DZV(:),
     *   F1X(:),F2X(:),F1Y(:),F2Y(:),KP2(:),KINT2(:),KP2INT(:),KP2p(:),
     *   Q(:),QINT(:),OMEGA1(:),OMEGA2(:),KV(:),KVSTAR(:),sx(:),sy(:),
     *   KP3(:),KP3INT(:),KINT3(:)
        INTEGER im,kmi
        REAL*8 Nm,factor
        REAL*8 zw(km),zwt(km),zt(km)
C
      pi=4.*atan(1.)
C
      if (ndepi.eq.0) goto 10

      do k=1,km
         zw(k) = z_cm(k)
         zt(k) =-z_cm(k)
      enddo
      zwt(1)=0.
      do k=2,km
        zwt(k) = zw(k-1)
      enddo

      kmi = ndepi

      kmli=km+1   !outside bounds
      DO k=1,kmi-1
        IF(zwt(k).le.ml.and.ml.lt.zwt(k+1)) kmli=k+1
      ENDDO
C
      kmli=min(kmli,kmi)
C
      ALLOCATE(n2t(kmi),n2t1(kmi),z(kmi),zm(kmi),b1(kmi),b1t(kmi),
     * U(kmi),V(kmi),rho(kmi),n2tr(kmi),
     * LX(kmi),LY(kmi),DZLX(kmi),DZLY(kmi),
     * UINT(kmi),VINT(kmi),UT(kmi),VT(kmi),
     * M1X(kmi),M1Y(kmi),M2X(kmi),M2Y(kmi),
     * dxb(kmi),dyb(kmi),KP(kmi),db(kmi),
     * zint(kmi),UREV(kmi),VREV(kmi),KPINT(kmi),
     * b1rev(kmi),DZLXREV(kmi),DZLYREV(kmi),
     * UTINT(kmi),VTINT(kmi),DZU(kmi),DZV(kmi),KP3(kmi),KP3INT(kmi),
     * F1X(kmi),F1Y(kmi),F2X(kmi),F2Y(kmi),KP2(kmi),KP2INT(kmi),
     * QINT(kmi),OMEGA1(kmi),OMEGA2(kmi),KV(kmi),KVSTAR(kmi),
     * sx(kmi),sy(kmi),KP2p(kmi))
C-- M1: Calculate B1, frd and rd
      n2t(1)=n2(1)
C     n2t(1)=max(n2t(1),1.d-8)
      DO k=2,kmi
        n2t(k)=n2(k-1)
C       n2t(k)=max(n2t(k),1.d-8)
      ENDDO
      DO k=1,kmi-1
        n2tr(k)=0.5*(n2(k)+n2(k+1))
      ENDDO
      n2tr(kmi)=n2(kmi-1)
C
      z=-zw(1:kmi)
      zm=-z
C     n2t=((7.d-3)**2)*exp(2.*z/100000.)
      n2t1=n2t
C
c     CALL baroclin1st(zm,n2t1,kmi,b1,im,Nm,frd)
      CALL baroclin1st(zm,n2(:),kmi,b1,im,Nm,frd,igrid,jgrid)
      rd=dabs(frd/f)
c     rdm=rd
      DO k=1,kmi-1
        b1t(k)=0.5*(b1(k)+b1(k+1))
      ENDDO
      b1t(kmi)=b1(kmi)
      b1=b1t
      U=UI(1:kmi)
      V=VI(1:kmi)
      rho=rhoi(1:kmi)
C-- Tony - 12/03/11 - Calulate the Rhines Scale and replace rd by lr between +/-30 degrees
c     lr=dsqrt(dsqrt(U(1)**2+V(1)**2)/0.95d-13)
c     IF(jgrid.ge.60.AND.jgrid.le.120) rd=lr
      rdm=rd
C--
C
C-- TONY - bug correction - 07/27/10
c     LX=-dxrho(1:kmi)/(-n2tr*rho/981.)
c     LY=-dyrho(1:kmi)/(-n2tr*rho/981.)
c
      CALL d1sym(kmi,z,U,DZU)
      CALL d1sym(kmi,z,V,DZV)
c
      do k=1,kmi
         if (n2tr(k).ne.0.d0)then 
                LX(k)=-(f/n2tr(k))*DZV(k)
                LY(k)=+(f/n2tr(k))*DZU(k)
         else
                LX(k)=0.d0
                LY(k)=0.d0
         endif
      enddo

c     do k=1,kmi
c        if (dzrho(k).ne.0.d0)then      !Natassa
c               sx(k)=-dxrho(k)/dzrho(k)
c               sy(k)=-dyrho(k)/dzrho(k)
c        else
c               sx=0.d0
c               sy=0.d0
c        endif
c     enddo
c     LX=sx
c     LY=sy
C
      UREV=0.
      VREV=0.
      DO k=1,kmi
        zint(k)=zt(kmi-k+1)
        UREV(k)=U(kmi-k+1)
        VREV(k)=V(kmi-k+1)
        b1rev(k)=b1(kmi-k+1)
      ENDDO
C
      z=zt(1:kmi)
      CALL B1AVERAGE(UREV,zint,b1rev,kmi,UB1AVE,kmli+0)
      CALL B1AVERAGE(VREV,zint,b1rev,kmi,VB1AVE,kmli+0)
C
      CALL d1sym(kmi,z,LX,DZLX)
      CALL d1sym(kmi,z,LY,DZLY)

c     if(i.eq.  1.and.j.eq.120) then
c     do k=1,kmi
c       write(*,'(a,i5,6e12.4)')'LX:',
c    .          k,-dxrho(k),-n2tr(k),rho(k),LX(k),z(k)
c    .           ,dzlx(k)
c     enddo
c     endif
c     if(i.eq.  4.and.j.eq.120) then
c     do k=1,kmi
c       write(*,'(a,i5,6e12.4)')'LX:',
c    .          k,-dxrho(k),-n2tr(k),rho(k),LX(k),z(k)
c    .           ,dzlx(k)
c     enddo
c     endif

      DO k=1,kmi
        DZLXREV(k)=DZLX(kmi-k+1)
        DZLYREV(k)=DZLY(kmi-k+1)
      ENDDO

c     if(i.eq.  1.and.j.eq.120)
c    .          write(*,'(a,3e12.4,i5,e12.4)')'DZLXB1AVE:', 
c    .          DZLXREV(1),zint(1),b1rev(1),kmi,DZLXB1AVE
c     if(i.eq.  4.and.j.eq.120)
c    .          write(*,'(a,3e12.4,i5,e12.4)')'DZLXB1AVE:', 
c    .          DZLXREV(1),zint(1),b1rev(1),kmi,DZLXB1AVE

      CALL B1AVERAGE(DZLXREV,zint,b1rev,kmi,DZLXB1AVE,kmli+0)
      CALL B1AVERAGE(DZLYREV,zint,b1rev,kmi,DZLYB1AVE,kmli+0)
C
C-- Calculate UT and VT
      DO k=1,kmi
        CALL dqtfg(zint(k:kmi),UREV(k:kmi),UINT(k:kmi),kmi-k+1)
        if (zint(k).ne.0.d0) then        !Natassa
        UTINT(k)=-(1./zint(k))*UINT(kmi)
        CALL dqtfg(zint(k:kmi),VREV(k:kmi),VINT(k:kmi),kmi-k+1)
        VTINT(k)=-(1./zint(k))*VINT(kmi)
        else
        UTINT(k)=0.d0
        VTINT(k)=0.d0
        endif
      ENDDO
      UTINT(kmi)=0.
      VTINT(kmi)=0.
      DO k=1,kmi
        UT(k)=UTINT(kmi-k+1)
        VT(k)=VTINT(kmi-k+1)
      ENDDO
C
C-- Calculate F1
      if (rd.ne.0.d0) then                !Natassa
         F1X = DZLXB1AVE + (1.+1./SIGMAT)*(1./f)*(1./rd**2)*(-VB1AVE)
      else
          F1X = 0.d0
      endif

c     if(i.eq.  1.and.j.eq.120)write(*,'(a,5e12.4)')
c    .     'f1x:',DZLXB1AVE,SIGMAT,f,rd,VB1AVE
c     if(i.eq.  4.and.j.eq.120)
c    .     write(*,'(a,5e12.4)')'f1x:',DZLXB1AVE,SIGMAT,f,rd,VB1AVE

      if (rd.ne.0.d0) then                !Natassa
          F1Y = DZLYB1AVE + (1.+1./SIGMAT)*(1./f)*(1./rd**2)*(+UB1AVE)
      else
          F1Y = 0.d0
      endif
c     F1X=0.
c     F1Y=0.
C-- Calculate F2
      if (rd.ne.0.d0) then                !Natassa
      F2X = (1./f)*(1./rd**2)*((1./SIGMAT)*(+VT)+V)
      F2Y = (1./f)*(1./rd**2)*((1./SIGMAT)*(-UT)-U)
      else
      F2X = 0.d0
      F2Y = 0.d0
      endif
c     F2X=0.
c     F2Y=0.
C-- Calculate DXB and DYB
      dxb=-981.*dxrho(1:kmi)/rho(1:kmi)
      dyb=-981.*dyrho(1:kmi)/rho(1:kmi)
C-- Calculate new K
c     CK=3.95
      CK=15.
      KP2=z*((F1X+F2X)*dxb+(F1Y+F2Y)*dyb)
c     KP2=-z*((F1X+F2X)*sx+(F1Y+F2Y)*sy)*n2tr
      DO k=1,kmi
        KP2INT(k)=KP2(kmi-k+1)
        QINT(k)=z(kmi-k+1)*n2tr(kmi-k+1)
      ENDDO
      kkpmax2=1
      kpmax2=0.
      DO k=1,kmi
        IF(dabs(KP2(k)).gt.kpmax2) THEN
          kpmax2=dabs(KP2(k))
          kkpmax2=k
        ENDIF
      ENDDO
C
      kkpint=max(1,kkpmax2)
c     kkpint=max(1,kmli)
c     kkpint=max(1,max(kmli,kkpmax2))

      ALLOCATE(KINT2(kkpint))
      CALL dqtfg(zint(kmi-(kkpint-1):kmi),KP2INT(kmi-(kkpint-1):kmi),
     *           KINT2,kkpint)
C--- Calculate Q
      ALLOCATE(Q(kmli+1))
      CALL dqtfg(zint(max(kmi-kmli,1):kmi),QINT(max(kmi-kmli,1):kmi),
     *           Q,kmli+1)
      if (rd.ne.0d0.and.ml.ne.0d0) then       !Natassa
      Q = (2.*CK/(ml*(f**2)*rd**2))*Q
      else
      Q = 0.d0
      endif
      Qmap=Q(kmli+1)
      if (rd.ne.0d0.and.z(kkpint).ne.0.d0) then       !Natassa
        K02=-(CK/(1.+Q(kmli+1)))*(rd**2)*(1./-z(kkpint))
     *         *KINT2(kkpint)
      else
        K02=0.d0
      endif

c     if (igrid.eq.10) write(*,*)jgrid,f
c     if (igrid.eq.215.and.jgrid.eq.73) then
c     write(*,*)'MESOSCALES:',igrid,jgrid,kmi,kmli,DZLXB1AVE,DZLYB1AVE,
c    .  K02
c     DO k=1,kmi
c     write(*,'(a,3i5,17e12.4)')'MESOSCALES:',
c    .   igrid,jgrid,k,z_cm(k),K02,LX(k),DZLX(k),LY(k),DZLY(k),
c    .   dxrho(k),dyrho(k),n2tr(k),rho(k),b1(k),sx(k),sy(k),
c    .   rd/1.d5,lr/1.d5,f,dsqrt(U(1)**2+V(1)**2)
c     ENDDO
c     endif
c     if (igrid.eq.216.and.jgrid.eq.73) then
c     write(*,*)
c     write(*,*)'MESOSCALES:',igrid,jgrid,kmi,kmli,DZLXB1AVE,DZLYB1AVE,
c    .  K02
c     DO k=1,kmi
c     write(*,'(a,3i5,17e12.4)')'MESOSCALES:',
c    .   igrid,jgrid,k,z_cm(k),K02,LX(k),DZLX(k),LY(k),DZLY(k),
c    .   dxrho(k),dyrho(k),n2tr(k),rho(k),b1(k),sx(k),sy(k),
c    .   rd/1.d5,lr/1.d5,f,dsqrt(U(1)**2+V(1)**2)
c     ENDDO
c     endif

      if (K02.lt.0d0)K02=0.d0     !Natassa
      KAPPAM = rdm * sqrt(K02)

c     if (i.eq.1.and.j.eq.120) then
c     do k=1,kmli
c     write(*,'(a,4i5,7e12.4)')'K02:',k,kkpint
c    *     ,CK,Q(kmli+1),rdm,z(kkpint),KINT2(kkpint)
c    *     ,K02,KAPPAM
c     enddo
c     endif

      DO k=1,kmli
        Ustar(k) =  KAPPAM * (F1X(k)+F2X(k))
        Vstar(k) =  KAPPAM * (F1Y(k)+F2Y(k))
        if ((dxb(k)**2+dyb(k)**2).ne.0.d0) then    !Natassa
        factor = (Ustar(k) * dxb(k) 
     *           +Vstar(k) * dyb(k))
     *         / (dxb(k)**2+dyb(k)**2)
        else
        factor = 0.d0
        endif
        Ustar_star(k) = Ustar(k) 
     *                    - factor * dxb(k)
        Vstar_star(k) = Vstar(k) 
     *                    - factor * dyb(k)

c     if (i.eq.1.and.j.eq.120) then
c     write(*,'(a,3i5,9e12.4)')'KAPPAM:',k,
c    *        rdm,
c    *        KAPPAM,F1X(k),F2X(k),Ustar(k),
c    *        dxb(k),dyb(k),factor,Ustar_star(k)
c     endif
c     if (i.eq.  4.and.j.eq.120) then
c     write(*,'(3i5,9e12.4)')k,
c    *        rdm,
c    *        KAPPAM,F1X(k),F2X(k),Ustar(k),
c    *        dxb(k),dyb(k),factor,Ustar_star(k)
c     endif
      ENDDO


C-- MKE
c     map=0.5*(U(1)**2+V(1)**2)

      DEALLOCATE(KINT2,Q)

C-- Calculate OMEGA1 and OMEGA2
      OMEGA1=z*((F1X+F2X)*LX+(F1Y+F2Y)*LY)

C-- Calculate KV and KVSTAR
      KVSTAR=rd*dsqrt(K02)*OMEGA1 

C
      DEALLOCATE(n2t,n2t1,z,zm,b1,b1t,
     *  U,V,rho,n2tr,
     *  LX,LY,DZLX,DZLY,
c    *  LXsm,LYsm,DZLXsm,DZLYsm,
     *  UINT,VINT,UT,VT,
     *  M1X,M1Y,M2X,M2Y,
     *  dxb,dyb,KP,
c    *  KINT,
     *  db,
     *  zint,UREV,VREV,KPINT,
     *  b1rev,DZLXREV,DZLYREV,
     *  UTINT,VTINT,DZU,DZV,sx,sy,KP3,KP3INT,KP2p,
     *  F1X,F1Y,F2X,F2Y,KP2,KP2INT,QINT,OMEGA1,OMEGA2,KV,KVSTAR)
C
 10   CONTINUE
 
      END  subroutine mesoscales1d


C*****
      SUBROUTINE B1AVERAGE(F,Z,B1,K,FAVE,KML)
      INTEGER K,KML
      REAL*8 F(K),Z(K),B1(K),FAVE,FB1AVE(K-KML+1),
     *  B1AVE(K-KML+1),a02
C
      a02=0.03
c     a02=0.5
C
      CALL dqtfg(Z(1:K-KML+1),F(1:K-KML+1)
     *           *dsqrt((a02+B1(1:K-KML+1)**2)/(1.+a02)),
     *           FB1AVE,K-KML+1)     
      CALL dqtfg(Z(1:K-KML+1),dsqrt((a02+B1(1:K-KML+1)**2)/(1.+a02)),
     *           B1AVE,K-KML+1)
C     WRITE(*,*),"B1AVE",FB1AVE/B1AVE
      FAVE=0.
      IF(B1AVE(K-KML+1).NE.0.)
     *FAVE=FB1AVE(K-KML+1)/B1AVE(K-KML+1)
C
      RETURN
      END
C*****
        SUBROUTINE d1sym(n,x,y,dyodx)
C080626 Input n,x(n),y(n)               Output dyodx(n)
C       Calculates the numerical ordinary derivative of y with respect to x, dyodx,
C       using the symmetrical approach of considering the nearest neighbor points
C       on both sides and assuming the derivative is changing linearly.
C       The ratio of differences on each side thus give the midpoint derivatives
C       on the respective sides and the averages of the midpoint derivatives
C       weighted each by the other midpoint's  distance give the derivative at the point.
C       Note for equal spacing same as ratio of neighbor differences with each other.
C       At top and bottom revert to ratio of differences with lone nearest neighbor.
        IMPLICIT NONE
        INTEGER n
        REAL*8 x(n),y(n),dyodx(n)
        REAL*8 dyodxm,dyodxp,dxm,dxp,dym,dyp

        INTEGER i


        DO i=1,n
           IF(i.EQ.1) THEN
              dyodx(i) = ((y(i+1)-y(i))/(x(i+1)-x(i)))
           ELSE IF(i.EQ.n) THEN
              dyodx(i) = ((y(i)-y(i-1))/(x(i)-x(i-1)))
           ELSE
              dxm = x(i)-x(i-1)
              dxp = x(i+1)-x(i)
              dym = y(i)-y(i-1)
              dyp = y(i+1)-y(i)
              dyodxm = dym/dxm
              dyodxp = dyp/dxp
              dyodx(i) = (dxm*dyodxp + dxp*dyodxm)/(dxm+dxp)
           END IF
        END DO

        RETURN
        END
C*****
C080625-30AH Subroutine to calculate the first baroclinic mode in CD WKB approximation
C       from an input N^2 profile, artificially treating N^2 as zero where it is less than
C       zero since its square and fourth roots appear in the WKB approximation formula,
C       B_1(z)=[N(z)/N_max]^1/2 cos[(f r_d)^-1 \Integral_z^z_max N(z) dz (where z_max was
C       called "z_0" by CD and take (f r_d) = {\Integral -H^z_max N(z) dz \over pi},
C       but with B_1(z) kept at 1 above z_max as formula wasn't intended for that range,
C       corrected and adapted from Cheng's 080522 program that calls dqtfg , named B1b.f:
!@sum B1b.f based on B1a.f but only reads form one data file
!@+   and ignors the last line in the data file
!@    5-22-2008

        SUBROUTINE baroclin1st(z,n2,m,ba1,im,Nm,frd,igrid,jgrid)
C       Inputs:
C               z(m)    !Depth                                                  [m]
C               n2(m)   !Square of Brunt Vaisala frequency                      [s^-2]
C
C       Outputs:
C               ba1(m)  !First baroclinic mode
C               im      !Index of maximum N on column
C               Nm      !Maximum N on column                                    [s^-2]
C               frd     !(f r_d)                                                [m/s]
C
C
C       Internal:
C               lifout  !Logical switch set .TRUE. iff diagnostic output here

      implicit real*8 (a-h,n,o-z)
C080625AH
        integer igrid,jgrid
        INTEGER m       ! Number of ocean levels
        REAL*8 ba1(m)   ! Array for first baroclinic mode as calculated by this routine
        LOGICAL lifout
        PARAMETER(lifout=.TRUE.)
      dimension n2(m)   !Array for N^2
      dimension z(m),n(m),ni(m) !Arrays for depth, SQRT(N^2), and z integral of N^2
        REAL*8 frd      !"f r_d" calculated in this routine
C******AH
C080625AH Work variables introduced.
        INTEGER im
C******AH
C080630AH Depth integrations are only done from maximum N level, "the pycnocline depth",
C       to the bottom so points above the depth of Maximum N are excluded unlike in B1b.f.
C       Introduce new arrays for depth and SQRT(N^2) which start at "pycnocline depth".
        REAL*8 zsubpyc(m),nsubpyc(m)
        INTEGER isubpyc(m),msubpyc
C******AH


      pi=acos(-1.d0)
      a=.1


C080625AH Canuto's WKB approximation formula contains (N^2)^1/2 and (N^2)^(1/4)
C       and therefore CANNOT HANDLE NEGATIVE N^2. Artificially set N^2=0 .
        DO i=1,m
           n2(i)=MAX(0.D0,n2(i))
        END DO
C******AH


      do i=1,m
         n(i)=sqrt(n2(i))
      end do

      ! find Nm = N_max
      Nm=N(1)
C080625AH Write out maximum N and depth when lifout=.TRUE. .
      im=1
      do i=2,m
         IF(N(i).GT.Nm) im=i
         Nm=max(N(i),Nm)
      end do
      IF(lifout) THEN
c       IF(ik.eq.177.and.jk.eq.163) THEN
c       write(*,*) "Nm=",Nm
c       WRITE(*,*) "im=",im
c       ENDIF
      END IF
C******AH

C080630AH Fill z & N arrays which go from max. N level down and integrate N on same range.
        DO i=im,m
           isubpyc(i)=i+1-im
           zsubpyc(isubpyc(i))=z(i)
           nsubpyc(isubpyc(i))=n(i)
        END DO
        msubpyc=m+1-im

      call dqtfg(zsubpyc,nsubpyc,ni,msubpyc)
C******AH

      do i=1,m ! z in meters
C80630AH Integrate N only over subpycnocline to find "f r_d". Keep B_1=1 above pycnocline.
         frd=ni(msubpyc)/pi
         IF(i.LT.im) THEN
           B1=1.D0
         ELSE
           if(frd.ne.0.d0) then           !Natassa
           arg=1./frd*ni(isubpyc(i))
           B1=(nsubpyc(isubpyc(i))/nm)**.5 * cos(arg)
           else
           B1=1.D0
           endif
         END IF
C******AH
C080625AH Store first baroclinic mode values in an array.
         ba1(i)=B1
C        Gam=1./(1+a) * (a+B1**2)
C        DMW=1./4.+3./4.*exp(-abs(z(i))/500.)
C        write(23,'(9e14.6)') z(i),Gam,sqrt(Gam)
C &        ,n2(i)/n2(1),DMW
C******AH
c        IF(ik.eq.177.and.jk.eq.163.and.i.eq.m) THEN
c          WRITE(*,*)'b1:ba1=',i,ba1(i)
c          WRITE(*,*)'b1:nsubpyc(isubpyc(i))=',nsubpyc(isubpyc(i))
c          WRITE(*,*)'b1:nm=',nm
c          WRITE(*,*)'b1:cos(arg)=',cos(arg)
c          WRITE(*,*)'b1:arg=',arg
c          WRITE(*,*)'b1:frd=',frd
c          WRITE(*,*)'b1:ni(isubpyc(i))=',ni(isubpyc(i))
c        ENDIF
         end do

c     IF(im.eq.m) WRITE(*,*)'b1:im,m,frd=',ik,jk,im,m,frd

      end
C*****
      subroutine dqtfg(x,y,z,n)
      !@sum integrate y from x(1) to x(i) and store it in z(i)
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n),z(n)
      sum2=0.
      if(n-1)4,3,1
    1 do 2 i=2,n
      sum1=sum2
      sum2=sum2+.5d0*(x(i)-x(i-1))*(y(i)+y(i-1))
    2 z(i-1)=sum1
    3 z(n)=sum2
    4 return
      end
