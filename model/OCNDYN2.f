#include "rundeck_opts.h"

!#define LUS_VERT_ADV

c Will add more documentation if this version becomes the modelE default.

      SUBROUTINE OCEANS
C****
      USE CONSTANT, only : rhows,grav
      USE MODEL_COM, only : modd5s,msurf,itime
      USE OCEANRES, only : NOCEAN
      USE OCEAN, only : im,jm,lmo,ndyno,mo,g0m,gxmo,gymo,gzmo,
     *    s0m,sxmo,symo,szmo,dts,dtofs,dto,dtolf,mdyno,msgso,
     *    ogeoz,ogeoz_sv,opbot,ze,lmm,imaxj, UO,VO,VONP,IVNP, ! VOSP,IVSP,
     *    OBottom_drag,OCoastal_drag,uod,vod,lmu,lmv
      USE OCEAN, only :
     &     nbyzm,nbyzu,nbyzv, i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      USE OCEAN_DYN, only : mmi,smu,smv,smw
      USE DOMAIN_DECOMP_1D, only : get, AM_I_ROOT, halo_update,
     &     south,north, hasSouthPole, hasNorthPole
      USE OCEANR_DIM, only : grid=>ogrid
      USE ODIAG, only : oijl=>oijl_loc,oij=>oij_loc,
     *    ijl_mo,ijl_g0m,ijl_s0m,  ijl_gflx, ijl_sflx, ijl_mfw2,
     *    ijl_mfu,ijl_mfv,ijl_mfw, ijl_ggmfl,ijl_sgmfl,ij_ssh,ij_pb

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : t_qlimit,ntm
      USE OCEAN, only : trmo,txmo,tymo,tzmo
      Use ODIAG, Only: toijl=>toijl_loc,
     *               toijl_conc,toijl_tflx,toijl_gmfl
#endif
#ifdef TRACERS_OceanBiology
      USE obio_com, only: gather_chl
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      USE obio_com, only: gather_pCO2
#endif

      IMPLICIT NONE
      Integer*4 I,J,L,N,NS,NST,NO,NEVEN ; real*8 now
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     &     MO1,MO2, UO1,UO2,UOD1,UOD2, VO1,VO2,VOD1,VOD2
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     OPBOT1,OPBOT2
      real*8 :: relfac,dt_odiff
      real*8, parameter :: byno=1./nocean

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0H,J_1H, J_0S,J_1S
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &     J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &     J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)

C***  Get the data from the atmospheric grid to the ocean grid
      call AG2OG_oceans

C***  Interpolate DYNSI outputs to the ocean grid
C***  (at present, only the ice-ocean stress is used)
      call IG2OG_oceans

c-------------------------------------------------------------------
c Begin ocean-processors-only code region
      ocean_processors_only: if(grid%have_domain) then
c-------------------------------------------------------------------

      OGEOZ_SV(:,:)=OGEOZ(:,:)

C**** Apply surface fluxes to ocean
      CALL GROUND_OC
         CALL CHECKO('GRNDOC')

C**** Apply ice/ocean and air/ocean stress to ocean
      CALL OSTRES2
         CALL CHECKO('OSTRES')
         CALL TIMER (NOW,MSURF)
         IF (MODD5S == 0) CALL DIAGCO (11)

C**** Apply ocean vertical mixing
      CALL OCONV
         CALL CHECKO('OCONV ')

C**** Apply bottom and coastal drags
      if (OBottom_drag  == 1) CALL OBDRAG2
      if (OCoastal_drag == 1) CALL OCOAST

C**** Add ocean biology
#ifdef TRACERS_OceanBiology
      call obio_model
      call gather_chl
#ifdef TRACERS_GASEXCH_ocean
      call gather_pco2
#endif
      IF (MODD5S.EQ.0) CALL DIAGCO (5)
#endif

         CALL TIMER (NOW,MSGSO)

c relax UOD,VOD toward 4-pt avgs of UO,VO
      call halo_update(grid,uo,from=north)
      call halo_update(grid,vo,from=south)
      relfac = .005d0
      do l=1,lmo
        if(hasNorthPole(grid)) then
          j = jm-1
          call polevel(uo(1,j_0h,l),vo(1,j_0h,l),l)
          i=1
          if(l.le.lmv(i,j)) then
            uod(i,j,l) = (1.-relfac)*uod(i,j,l) + relfac*.25*(
     &           uo(im,j,l)+uo(i,j,l)+2.*uo(i,j+1,l)
     &           )
          endif
          do n=1,nbyzv(j,l)
            do i=max(2,i1yzv(n,j,l)),i2yzv(n,j,l)
              uod(i,j,l) = (1.-relfac)*uod(i,j,l) + relfac*.25*(
     &             uo(i-1,j,l)+uo(i,j,l)+2.*uo(i,j+1,l)
     &             )
            enddo
          enddo
        endif
        do j=j_0s,min(jm-2,j_1s)
          i=1
          if(l.le.lmv(i,j)) then
            uod(i,j,l) = (1.-relfac)*uod(i,j,l) + relfac*.25*(
     &           uo(im,j,l)+uo(i,j,l)+uo(im,j+1,l)+uo(i,j+1,l)
     &           )
          endif
          do n=1,nbyzv(j,l)
            do i=max(2,i1yzv(n,j,l)),i2yzv(n,j,l)
              uod(i,j,l) = (1.-relfac)*uod(i,j,l) + relfac*.25*(
     &             uo(i-1,j,l)+uo(i,j,l)+uo(i-1,j+1,l)+uo(i,j+1,l)
     &             )
            enddo
          enddo
        enddo
        do j=j_0s,j_1s
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),min(im-1,i2yzu(n,j,l))
              vod(i,j,l) = (1.-relfac)*vod(i,j,l) + relfac*.25*(
     &             vo(i,j-1,l)+vo(i+1,j-1,l)+vo(i,j,l)+vo(i+1,j,l)
     &             )
            enddo
          enddo
          i=im
          if(l.le.lmu(i,j)) then
            vod(i,j,l) = (1.-relfac)*vod(i,j,l) + relfac*.25*(
     &           vo(i,j-1,l)+vo(1,j-1,l)+vo(i,j,l)+vo(1,j,l)
     &           )
          endif
        enddo
      enddo

      DO L=1,LMO
        MO1(:,:,L) = 0
        UO1(:,:,L) = 0
        VO1(:,:,L) = 0
        UOD1(:,:,L) = 0
        VOD1(:,:,L) = 0
        MO2(:,:,L) = 0
        UO2(:,:,L) = 0
        VO2(:,:,L) = 0
        UOD2(:,:,L) = 0
        VOD2(:,:,L) = 0
      EndDo

C****
C**** Integrate Ocean Dynamics
C****
      Do NO=1,NOCEAN

c check which of these are necessary
      call halo_update(grid,mo)
      call halo_update(grid,uo)
      call halo_update(grid,vo)

      CALL ODHORZ0

      do l=1,lmo
        SMU(:,J_0:J_1,L) = 0
        SMV(:,J_0H:J_1,L) = 0
c        SMW(:,J_0:J_1H,L) = 0 ! not summed
        do j=max(1,j_0h),min(jm,j_1h)
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              mo1(i,j,l) = mo(i,j,l)
              mo2(i,j,l) = mo(i,j,l)
            enddo
          enddo
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),i2yzu(n,j,l)
              uo1(i,j,l) = uo(i,j,l)
              uo2(i,j,l) = uo(i,j,l)
              vod1(i,j,l) = vod(i,j,l)
              vod2(i,j,l) = vod(i,j,l)
            enddo
          enddo
          do n=1,nbyzv(j,l)
            do i=i1yzv(n,j,l),i2yzv(n,j,l)
              vo1(i,j,l) = vo(i,j,l)
              vo2(i,j,l) = vo(i,j,l)
              uod1(i,j,l) = uod(i,j,l)
              uod2(i,j,l) = uod(i,j,l)
            enddo
          enddo
        enddo
        if(hasNorthPole(grid)) then
          mo1(:,jm,l) = mo(:,jm,l)
          mo2(:,jm,l) = mo(:,jm,l)
          uo1(:,jm,l) = uo(:,jm,l)
          uo2(:,jm,l) = uo(:,jm,l)
          vo1(:,jm,l) = vo(:,jm,l)
          vo2(:,jm,l) = vo(:,jm,l)
        endif
      enddo
      do j=max(1,j_0h),min(jm,j_1h)
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            opbot1(i,j) = opbot(i,j)
            opbot2(i,j) = opbot(i,j)
          enddo
        enddo
      enddo

c
c advance the short-timestep horizontal dynamics
c
c initialize the odd state
      Call ODHORZ(MO ,UO ,VO ,UOD ,VOD ,OPBOT ,
     &            MO2,UO2,VO2,UOD2,VOD2,OPBOT2, DTOFS,.false.)
      Call ODHORZ(MO2,UO2,VO2,UOD2,VOD2,OPBOT2,
     &            MO1,UO1,VO1,UOD1,VOD1,OPBOT1, DTO,.false.)
c loop over the leapfrog steps
      neven = NDYNO / (2*NOCEAN)
      do n=1,neven
c update the even state
        Call ODHORZ(MO1,UO1,VO1,UOD1,VOD1,OPBOT1,
     &              MO ,UO ,VO ,UOD ,VOD ,OPBOT , DTOLF,.true.)
        if(n == neven) exit ! no need to further update the odd state
c update the odd state
        Call ODHORZ(MO ,UO ,VO ,UOD ,VOD ,OPBOT ,
     &              MO1,UO1,VO1,UOD1,VOD1,OPBOT1, DTOLF,.false.)
      enddo

c is this still needed?
      if(j_1 == JM) UO(IVNP,JM,:) = VONP(:) ! not needed if Mod(IM,4)=0

c
c long-timestep vertical redistribution of mass (and momentum)
c
      Call OFLUXV

c
c long-timestep advection of potential enthalpy, salt, and tracers
c
      CALL OADVT2 (G0M,GXMO,GYMO,GZMO,DTOLF,.FALSE.
     *        ,OIJL(1,J_0H,1,IJL_GFLX))
      CALL OADVT2 (S0M,SXMO,SYMO,SZMO,DTOLF,.TRUE.
     *        ,OIJL(1,J_0H,1,IJL_SFLX))
#ifdef TRACERS_OCEAN
      DO N=1,NTM
        CALL OADVT2(TRMO(1,J_0H,1,N),TXMO(1,J_0H,1,N)
     *       ,TYMO(1,J_0H,1,N),TZMO(1,J_0H,1,N),DTOLF,t_qlimit(n)
     *       ,TOIJL(1,J_0H,1,TOIJL_TFLX,N))
      END DO
#endif
        CALL CHECKO ('OADVT ')

c
c diagnostics
c
      DO L=1,LMO
        do j=j_0s,j_1s
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),i2yzu(n,j,l)
              OIJL(I,J,L,IJL_MFU) = OIJL(I,J,L,IJL_MFU) + SMU(I,J,L)
            enddo
          enddo
        enddo
        do j=j_0,j_1
          do n=1,nbyzv(j,l)
            do i=i1yzv(n,j,l),i2yzv(n,j,l)
              OIJL(I,J,L,IJL_MFV) = OIJL(I,J,L,IJL_MFV) + SMV(I,J,L)
            enddo
          enddo
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              OIJL(I,J,L,IJL_MFW) = OIJL(I,J,L,IJL_MFW) + SMW(I,J,L)
              OIJL(I,J,L,IJL_MFW2)= OIJL(I,J,L,IJL_MFW2)+
     *             SMW(I,J,L)*SMW(I,J,L)
              OIJL(I,J,L,IJL_MO)  = OIJL(I,J,L,IJL_MO) +  MO(I,J,L)*byno
              OIJL(I,J,L,IJL_G0M) = OIJL(I,J,L,IJL_G0M) +G0M(I,J,L)*byno
              OIJL(I,J,L,IJL_S0M) = OIJL(I,J,L,IJL_S0M) +S0M(I,J,L)*byno
            enddo
          enddo
        enddo
      ENDDO

      do j=j_0,j_1
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            OIJ(I,J,IJ_SSH) = OIJ(I,J,IJ_SSH) + OGEOZ(I,J)*byno
            OIJ(I,J,IJ_PB)  = OIJ(I,J,IJ_PB)  +
     &           (OPBOT(I,J)-ZE(LMM(I,J))*RHOWS*GRAV)*byno
          enddo
        enddo
      enddo

#ifdef TRACERS_OCEAN
        DO N=1,NTM
          DO L=1,LMO
            TOIJL(:,:,L,TOIJL_CONC,N)=TOIJL(:,:,L,TOIJL_CONC,N)
     *           +TRMO(:,:,L,N)*byno
          END DO
        END DO
#endif

      call gather_ocean_straits()

      IF(AM_I_ROOT()) THEN
C****
C**** Acceleration and advection of tracers through ocean straits
C****
          CALL STPGF(DTS/NOCEAN)
          CALL STADV(DTS/NOCEAN)
          CALL CHECKO_serial ('STADV0')
        IF (NO .EQ. NOCEAN) THEN
          CALL STCONV
          CALL STBDRA
        END IF
      END IF
      call scatter_ocean_straits()
      call BCAST_straits (0)
        CALL CHECKO ('STADV ')

      ENDDO  !  End of Do-loop NO=1,NOCEAN

        CALL TIMER (NOW,MDYNO)
        IF (MODD5S == 0) CALL DIAGCO (12)

c
c recalculate vbar etc. for ocean physics
c
      CALL ODHORZ0

C**** Apply Wajowicz horizontal diffusion to UO and VO ocean currents
C**** every 3 hours
      dt_odiff = 3.*3600.
      if(mod(itime,int(dt_odiff/dts)).eq.0) then
        CALL ODIFF(dt_odiff)
      endif
c      CALL OABFILx ! binomial filter
c      CALL OABFILy ! binomial filter
c      CALL CHECKO ('ODIFF0')

C**** Apply GM + Redi tracer fluxes
      CALL GMKDIF
      CALL GMFEXP(G0M,GXMO,GYMO,GZMO,.FALSE.,OIJL(1,J_0H,1,IJL_GGMFL))
      CALL GMFEXP(S0M,SXMO,SYMO,SZMO,.TRUE. ,OIJL(1,J_0H,1,IJL_SGMFL))
#ifdef TRACERS_OCEAN
      DO N = 1,NTM
        CALL GMFEXP(TRMO(1,J_0H,1,N),TXMO(1,J_0H,1,N),TYMO(1,J_0H,1,N),
     *    TZMO(1,J_0H,1,N),t_qlimit(n),TOIJL(1,J_0H,1,TOIJL_GMFL,N))
      END DO
#endif
      CALL CHECKO ('GMDIFF')

#ifdef TRACERS_OCEAN
      CALL OC_TDECAY(DTS)
#ifdef TRACERS_AGE_OCEAN
      CALL OCN_TR_AGE(DTS)
#endif
#endif

#ifdef OCN_Mesoscales
      CALL OCN_mesosc
#endif
      CALL TIMER (NOW,MSGSO)

c-------------------------------------------------------------------
c End ocean-processors-only code region
c-------------------------------------------------------------------

      else ! no ocean domain, call the timers anyway
        CALL TIMER (NOW,MSURF)
        CALL TIMER (NOW,MSGSO)
        CALL TIMER (NOW,MDYNO)
        CALL TIMER (NOW,MSGSO)
      endif ocean_processors_only

C***  Get the data from the ocean grid to the atmospheric grid
      CALL TOC2SST
      call OG2AG_oceans

C***  Interpolate ocean surface velocity to the DYNSI grid
      call OG2IG_uvsurf

      RETURN
      END SUBROUTINE OCEANS

      Subroutine OFLUXV
      use constant, only : grav
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM,LMOU=>LMU,LMOV=>LMV,
     &     ZOE=>ZE,dZO, DTS, DXYPO, BYDXYPO, OPRESS, dtolf
      Use OCEAN, only : opbot,mo,uo,vo
      USE OCEAN, only :
     &     nbyzm,nbyzu,nbyzv,
     &     i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      Use OCEAN_DYN, Only: SMW
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH,SOUTH
      USE OCEANR_DIM, only : grid=>ogrid
      Implicit None
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: msum
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     mb,mtmp,mwtmp
      Real*8 :: mwfac,mfinal
      Integer*4 I,J,L,LM, J1,JN,J1P,J1H,JNP,JNQ,JNH, N
      Logical QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      J1H = Max(J1-1,1)     !    1      4     JM-4   Halo minimum
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNQ = Min(JN,JM-2)    !    4      8     JM-2   Exclude NP,NP-1
      JNH = Min(JN+1,JM)
      QNP = JN==JM          !    F      F      T
C****

      do l=1,lmo
        do j=j1p,jn
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              mb(i,j,l) = mo(i,j,l)
            enddo
          enddo
        enddo
        if(qnp) then
          j = jm
          do i=2,im
            mb(i,j,l) = mo(i,j,l)
          enddo
        endif
      enddo

      do j=j1p,jn !jnh
        mwfac = dxypo(j)/dtolf
        do n=1,nbyzm(j,2)
          do i=i1yzm(n,j,2),i2yzm(n,j,2)
            msum(i,j) = (opbot(i,j)-opress(i,j))/grav
            lm=lmom(i,j)
            mfinal = msum(i,j)*dzo(1)/zoe(lm)
c            mw(i,j) = mwfac*(mo(i,j,1)-mfinal)
            smw(i,j,1) = mwfac*(mo(i,j,1)-mfinal)
            mo(i,j,1) = mfinal
          enddo
        enddo
      enddo

      do l=2,lmo-1

      do j=j1p,jn !jnh
        mwfac = dxypo(j)/dtolf
        do n=1,nbyzm(j,l+1)
          do i=i1yzm(n,j,l+1),i2yzm(n,j,l+1)
            lm=lmom(i,j)
            mfinal = msum(i,j)*dzo(l)/zoe(lm)
c            mw(i,j) = mw(i,j) + mwfac*(mo(i,j,l)-mfinal)
            smw(i,j,l) = smw(i,j,l-1) + mwfac*(mo(i,j,l)-mfinal)
            mo(i,j,l) = mfinal
          enddo
        enddo
      enddo

      enddo ! l

c
c update bottom layer mass
c
      do j=j1p,jn !jnh
        do n=1,nbyzm(j,2)
          do i=i1yzm(n,j,2),i2yzm(n,j,2)
            lm=lmom(i,j)
            mfinal = msum(i,j)*dzo(lm)/zoe(lm)
            mo(i,j,lm) = mfinal
          enddo
        enddo
      enddo

      if(qnp)  then ! fill pole
        j = jm
        do l=1,lmom(1,j)
          do i=2,im
            mo(i,j,l) = mo(1,j,l)
            smw(i,j,l) = smw(1,j,l)
          enddo
        enddo
      endif

c
c for now, vertical advection of U,V uses the simplest upstream scheme
c
      call halo_update(grid,smw,from=north)
      call halo_update(grid,mb,from=north)
      do l=1,lmo
        do j=j1p,jnp
          mwfac = bydxypo(j)
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),min(i2yzu(n,j,l),im-1)
              mtmp(i,j,l) = .5*(mb(i,j,l)+mb(i+1,j,l))
              mwtmp(i,j,l) = .5*(smw(i,j,l)+smw(i+1,j,l))*mwfac
            enddo
            i = im
            if(l <= lmou(i,j)) then
              mtmp(i,j,l) = .5*(mb(i,j,l)+mb(1,j,l))
              mwtmp(i,j,l) = .5*(smw(i,j,l)+smw(1,j,l))*mwfac
            endif
          enddo
        enddo
      enddo
      do j=j1p,jnp
        do n=1,nbyzu(j,1)
          do i=i1yzu(n,j,1),i2yzu(n,j,1)
            mwtmp(i,j,lmou(i,j)) = 0.
          enddo
        enddo
      enddo
      call oadvuz(uo,mtmp,mwtmp,dtolf,j1p,jnp,nbyzu,i1yzu,i2yzu)
      do l=1,lmo
        do j=j1p,jnp
          do n=1,nbyzv(j,l)
            do i=i1yzv(n,j,l),i2yzv(n,j,l)
              mtmp(i,j,l) = .5*(mb(i,j,l)+mb(i,j+1,l))
              mwtmp(i,j,l) =
     &             .5*(smw(i,j,l)*bydxypo(j)+smw(i,j+1,l)*bydxypo(j+1))
            enddo
          enddo
        enddo
      enddo
      do j=j1p,jnp
        do n=1,nbyzv(j,1)
          do i=i1yzv(n,j,1),i2yzv(n,j,1)
            mwtmp(i,j,lmov(i,j)) = 0.
          enddo
        enddo
      enddo
      call oadvuz(vo,mtmp,mwtmp,dtolf,j1p,jnp,nbyzv,i1yzv,i2yzv)

      Return
      End Subroutine OFLUXV

      Subroutine OPFIL2(X,L,JMIN,JMAX)
C****
C**** OPFIL smoothes X in the zonal direction by reducing coefficients
C**** of its Fourier series for high wave numbers near the poles.
C****
      Use CONSTANT, Only: TWOPI
      Use OCEAN, Only: IM,JM,LMO,J1O, DLON, DXP=>DXPO, DYP=>DYPO,
     *  JMPF=>J40S  !  greatest J in SH where polar filter is applied
      Use FILEMANAGER, Only: OPENUNIT, CLOSEUNIT

      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
C****
      Real*8,Intent(InOut) :: X(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      Integer*4, Intent(In) :: L,JMIN,JMAX
      Integer*4,Parameter :: IMz2=IM/2  !  highest possible wave number
      Integer*4,Save ::
     *  JMEX,    !  exclude unfiltered latitudes, store J=JM-1 in JMEX
     *  NBASM,   !  maximum number of ocean basins at any latitude
     *  INDM,    !  number of entries in array REDUCO
     *  NMIN(JM) !  minimum wave number for Fourier smoothing
      Real*8,Save :: SMOOTH(IMz2,JM)  !  multiplies Fourier coefficient
C****
      Integer*2,Save,Allocatable,Dimension(:,:) ::
     *  NBAS    !  number of ocean basins for given layer and latitude
      Integer*2,Save,Allocatable,Dimension(:,:,:) ::
     *  IMINm1, !  western most cell in basin minus 1
     *  IWIDm1  !  number of cells in basin minus 1
      Integer*4,Save,Allocatable,Dimension(:,:) ::
     *  INDEX   !  index to concatenated matrices
      Real*4,Save,Allocatable,Dimension(:) ::
     *  REDUCO  !  concatenation of reduction matrices
      Real*4,Allocatable,Dimension(:) ::
     *  REDUCO_glob  !  concatenation of reduction matrices
C****
      Character*80 TITLE
      Integer*4 I,I1,INDX,IWm2,IWm1, J,JA,JX, K, LL, N,NB, IU_AVR
      integer :: i2,k1,k2
      integer, parameter :: hwid_max=15
      integer :: indx_min,indx_max
      Real*8 AN(0:IMz2), !  Fourier cosine coefficients
     *       BN(0:IMz2), !  Fourier sine coefficients
     *          Y(IM*2), !  original copy of X that wraps around IDL
     *       DRAT, REDUC
C****
      Integer*4,Save :: IFIRST=1, J1,JN,J1P,JNP
      If (IFIRST == 0)  GoTo 100
      IFIRST = 0
      JMEX = 2*JMPF-1
C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
C**** Calculate SMOOTHing factor for longitudes without land cells
      Do 30 J=J1O,JM-1
      DRAT = DXP(J)/DYP(3)
      Do 20 N=IMz2,1,-1
      SMOOTH(N,J) = IMz2*DRAT/N
   20 If (SMOOTH(N,J) >= 1)  GoTo 30
C     N = 0
   30 NMIN(J) = N+1
C**** Read in reduction contribution matrices from disk.  Only keep
C**** the matrices needed for the latitudes on this processor.
      Call OPENUNIT ('AVR',IU_AVR,.True.,.True.)
      Read (IU_AVR) TITLE,NBASM,INDM
      Allocate (IMINm1(NBASM,LMO,J1O:JMEX), NBAS(LMO,J1O:JMEX),
     *          IWIDm1(NBASM,LMO,J1O:JMEX), INDEX(IM,2:JMPF),
     *          REDUCO_glob(INDM))
      Read (IU_AVR) TITLE,NBAS,IMINm1,IWIDm1,INDEX,REDUCO_glob
      Call CLOSEUNIT (IU_AVR)
      Write (6,*) 'Read from unit',IU_AVR,': ',TITLE
      indx_min = indm+1; indx_max = -1
      do J=Max(J1O,J1-1),Min(JN+1,JM-1)
        JX=J  ;  If(J > JMPF) JX=J+2*JMPF-JM
        JA=J  ;  If(J > JMPF) JA=JM+1-J
        If (JA > JMPF)  cycle   !  skip latitudes J=JMPF+1,JM-JMPF
        do LL=1,LMO
          If (IWIDm1(1,LL,JX) >= IM)  cycle
          do NB=1,NBAS(LL,JX)
            IWm1 = IWIDm1(NB,LL,JX)
            indx_min = min(indx_min,INDEX(IWm1,JA)+1)
            indx_max = max(indx_max,INDEX(IWm1,JA)+IWm1**2)
          enddo
        enddo
      enddo
      if(indx_min.le.indx_max) then
        allocate(reduco(indx_min:indx_max))
        reduco(indx_min:indx_max) = reduco_glob(indx_min:indx_max)
      endif
      deallocate(reduco_glob)
 100  CONTINUE
C****
C**** Loop over J.  JX = eXclude unfiltered latitudes
C****               JA = Absolute latitude
C****
      Do J=Max(J1O,JMIN),JMAX !Max(J1O,J1P),JNP
      JX=J  ;  If(J > JMPF) JX=J+2*JMPF-JM
      JA=J  ;  If(J > JMPF) JA=JM+1-J
      If (JA > JMPF)  cycle  !  skip latitudes J=JMPF+1,JM-JMPF

      If (IWIDm1(1,L,JX) >= IM)  then
C****
C**** No land cells at this latitude and layer,
C**** perform standard polar filter
C****
        Call OFFT (X(1,J),AN,BN)
        Do N=NMIN(J),IMz2-1
          AN(N) = AN(N)*SMOOTH(N,J)
          BN(N) = BN(N)*SMOOTH(N,J)
        enddo
        AN(IMz2) = AN(IMz2)*SMOOTH(IMz2,J)
        Call OFFTI (AN,BN,X(1,J))

      else
C****
C**** Land cells exist at this latitude and layer, loop over ocean
C**** basins.
C****
        Do NB=1,NBAS(L,JX)
          I1   = IMINm1(NB,L,JX) + 1
          IWm2 = IWIDm1(NB,L,JX) - 1
          I2   = i1+iwm2
          INDX = INDEX(IWm2+1,JA)

          If (I2 <= IM)  then
C**** Ocean basin does not wrap around the IDL.
C**** Copy X to temporary array Y and filter X in place.
            do i=i1,i2
              y(i) = x(i,j)
            enddo
            do i=i1,i2
              reduc = 0
              k1 = max(i1,i-hwid_max)
              k2 = min(i2,i+hwid_max)
              indx = indx + k1-i1
              do k=k1,k2
                indx=indx+1
                reduc = reduc + reduco(indx)*y(k)
              enddo
              indx = indx + i2-k2
              x(i,j) = x(i,j) - reduc
            enddo

          else
C**** Ocean basin wraps around the IDL.
C**** Copy X to temporary array Y and filter X in place.
            do i=i1,im
              y(i) = x(i,j)
            enddo
            do i=im+1,i2
              y(i) = x(i-im,j)
            enddo
            do i=i1,im
              reduc = 0
              k1 = max(i1,i-hwid_max)
              k2 = min(i2,i+hwid_max)
              indx = indx + k1-i1
              do k=k1,k2
                indx=indx+1
                reduc = reduc + reduco(indx)*y(k)
              enddo
              indx = indx + i2-k2
              x(i,j) = x(i,j) - reduc
            enddo
            do i=1,i2-im
              reduc = 0
              k1 = max(i1,im+i-hwid_max)
              k2 = min(i2,im+i+hwid_max)
              indx = indx + k1-i1
              do k=k1,k2
                indx=indx+1
                reduc = reduc + reduco(indx)*y(k)
              enddo
              indx = indx + i2-k2
              x(i,j) = x(i,j) - reduc
            enddo
          endif     ! wraparound or not
        enddo       ! end loop over basins
      endif         ! fft vs matrix method
      enddo         ! end loop over j

      Return
      EndSubroutine OPFIL2

      Subroutine ODHORZ(MOH,UOH,VOH,UODH,VODH,OPBOTH,
     &                  MO ,UO ,VO ,UOD ,VOD ,OPBOT, DT,qeven)
      Use CONSTANT, Only: GRAV,omega
      Use OCEAN, Only: IM,JM,LMO,
     *                 LMOM=>LMM,LMOU=>LMU,LMOV=>LMV,
     *                 mZSOLID=>HOCEAN,OGEOZ,
     &     SINVO,SINPO,DXPO,DYPO,DXYVO,DXVO,DYVO,DXYPO
      USE OCEAN, only :
     &     nbyzm,nbyzu,nbyzv,nbyzc,
     &     i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv, i1yzc,i2yzc
      Use OCEAN_DYN, Only: SMU,SMV,VBAR,dZGdP,MMI
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH, SOUTH
     &     ,GLOBALMAX
      USE OCEANR_DIM, only : grid=>ogrid
      Implicit None
      Real*8, Intent(In),
     &     Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     &     MOH,UOH,VOH,UODH,VODH
      Real*8, Intent(In),
     &     Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: OPBOTH
      Real*8, Intent(Inout),
     &  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     &     MO,UO,VO,UOD,VOD
      Real*8, Intent(Inout),
     &     Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: OPBOT
      Real*8,Intent(In) :: DT
      Logical, Intent(In) :: qeven
c
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     ZG,P,PDN,DH,USMOOTH,PGFX,PGFY,MU,MV,UA,VA,VORT,KE
      Real*8 corofj,mufac,mvfac,bydx,bydy,mmid,xeven,
     &     convij,convfac,dp,pgf4pt,pgfac,uq,vq,uasmooth
      Integer*4 I,J,L, J1,JN,J1P,JNP,J1H,JNH, j1a, jminpfu,jmaxpfu, N
      integer :: smallmo_loc,smallmo
      Logical QNP

C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNH = Min(JN+1,JM)    !    5      9     JM     Halo maximum
      QNP = JN==JM          !    F      F      T
      JMINpfu = max(2,J1-1)
      JMAXpfu = min(JNH,JM-1)
      j1h = max(1,grid%j_strt_halo)
      j1a = grid%j_strt_halo

      if(qeven) then
        xeven = 1.
      else
        xeven = 0.
      endif

c
c initialize pressure and geopotential at the ocean bottom
c
      do j=j1h,jnh
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            pdn(i,j) = opboth(i,j)
            OGEOZ(I,J) = - mZSOLID(I,J)*GRAV
          enddo
        enddo
      enddo

! zero only once with bottom-up looping. see whether these are
! all necessary
      usmooth(:,:) = 0.
      mu(:,:) = 0.
      mv(:,:) = 0.
      pgfx(:,:) = 0.
      pgfy(:,:) = 0.
      vort(:,:) = 0.

      smallmo_loc = 0
c
c loop over layers, starting at the ocean bottom
c
      do l=lmo,1,-1

C**** Apply polar filter to West-East velocity
      if(any(nbyzu(jminpfu:jmaxpfu,l).gt.0)) then
        Do J=jminpfu,jmaxpfu
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),i2yzu(n,j,l)
              usmooth(I,J) = UOH(I,J,L)
            enddo
          enddo
        enddo
        Call OPFIL2(usmooth,L,jminpfu,jmaxpfu)
      endif
      if(qnp) then
        usmooth(:,jm) = uoh(:,jm,l)
      endif

c
c calculate pressure, geopotential, and kinetic energy
c
      do j=j1h,jnh
        i = 1
        if(l <= lmom(i,j)) then
          dp = moh(i,j,l)*grav
          dh(i,j) = moh(i,j,l)*vbar(i,j,l)
          p(i,j) = pdn(i,j) - .5*dp
          zg(i,j) = ogeoz(i,j) + dp*.5*dzgdp(i,j,l)
          pdn(i,j) = pdn(i,j) - dp
          ogeoz(i,j) = ogeoz(i,j) + dh(i,j)*grav
        endif
        do n=1,nbyzm(j,l)
          do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
            dp = moh(i,j,l)*grav
            dh(i,j) = moh(i,j,l)*vbar(i,j,l)
            p(i,j) = pdn(i,j) - .5*dp
            zg(i,j) = ogeoz(i,j) + dp*.5*dzgdp(i,j,l)
            pdn(i,j) = pdn(i,j) - dp
            ogeoz(i,j) = ogeoz(i,j) + dh(i,j)*grav
          enddo
        enddo
      enddo
      do j=j1,jnh
        i = 1
        if(l <= lmom(i,j)) then
          uasmooth = .5*(usmooth(im,j)+usmooth(i,j))
          ua(i,j) = .5*(uoh(im,j,l)+uoh(i,j,l))
          va(i,j) = .5*(voh(i,j-1,l)+voh(i,j,l))
          ke(i,j) = .5*(ua(i,j)*uasmooth + va(i,j)**2)
        endif
        do n=1,nbyzm(j,l)
          do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
            uasmooth = .5*(usmooth(i-1,j)+usmooth(i,j))
            ua(i,j) = .5*(uoh(i-1,j,l)+uoh(i,j,l))
            va(i,j) = .5*(voh(i,j-1,l)+voh(i,j,l))
            ke(i,j) = .5*(ua(i,j)*uasmooth + va(i,j)**2)
          enddo
        enddo
      enddo

c
c fill pole
c
      if(qnp .and. l <= lmom(1,jm)) then
        j = jm
        do i=2,im
          dh(i,j) = dh(1,j)
          p(i,j) = p(1,j)
          zg(i,j) = zg(1,j)
          ke(i,j) = ke(1,j)
          ua(i,j) = .5*(usmooth(i-1,j)+usmooth(i,j))
        enddo
      endif

c
c store, and smooth, the east-west pressure gradient
c
      DO J=J1P,JMAXpfu
        mufac = .5*dypo(j)
        bydx = 1d0/dxpo(j)
        do n=1,nbyzu(j,l)
          do i=i1yzu(n,j,l),min(im-1,i2yzu(n,j,l))
            mmid = MOH(I,J,L)+MOH(I+1,J,L)
            MU(I,J) = mufac*usmooth(i,j)*mmid
            SMU(I,J,L) = SMU(I,J,L) + MU(I,J)*xeven
            pgfx(i,j) =
     &           (ZG(I,J)-ZG(I+1,J))+(P(I,J)-P(I+1,J))*
     &           (dH(I,J)+dH(I+1,J))/mmid
            pgfx(i,j) = pgfx(i,j)*bydx
          enddo
        enddo
        I=IM
        if(l <= lmou(i,j)) then
          mmid = MOH(I,J,L)+MOH(1,J,L)
          MU(I,J) = mufac*usmooth(i,j)*mmid
          SMU(I,J,L) = SMU(I,J,L) + MU(I,J)*xeven
          pgfx(i,j) =
     &         (ZG(I,J)-ZG(1,J))+(P(I,J)-P(1,J))*
     &         (dH(I,J)+dH(1,J))/mmid
          pgfx(i,j) = pgfx(i,j)*bydx
        endif
      ENDDO
      if(any(nbyzu(j1p:JMAXpfu,l).gt.0)) then
        Call OPFIL2(pgfx,l,J1P,JMAXpfu)
      endif

      do j=j1p-1,jnp
        bydy = 1d0/dyvo(j)
        if(j == jm-1) bydy = bydy*2./3.
        do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            mmid = (moh(i,j,l)+moh(i,j+1,l))
            pgfy(i,j) = (ZG(I,J)-ZG(I,J+1))+(P(I,J)-P(I,J+1))*
     &             (dH(I,J)+dH(I,J+1))/mmid
            pgfy(i,j) = pgfy(i,j)*bydy
          enddo
        enddo
      enddo

c
c calculate vorticity at cell corners
c
      do j=j1p-1,jnp  ! merge with pgfy loop
        do n=1,nbyzc(j,l)
          do i=i1yzc(n,j,l),min(im-1,i2yzc(n,j,l))
            vort(i,j) = (dxpo(j)*usmooth(i,j)-dxpo(j+1)*usmooth(i,j+1)
     &                  +dyvo(j)*(voh(i+1,j,l)-voh(i,j,l)))/dxyvo(j)
          enddo
        enddo
        i=im
        vort(i,j) = (dxpo(j)*usmooth(i,j)-dxpo(j+1)*usmooth(i,j+1)
     &              +dyvo(j)*(voh(1,j,l)-voh(i,j,l)))/dxyvo(j)
      enddo

c
c update UO, VOD
c
      Do J=J1P,JNP
        bydx = 1d0/dxpo(j)
        corofj = 2.*omega*sinpo(j)
        do n=1,nbyzu(j,l)
          do i=i1yzu(n,j,l),min(im-1,i2yzu(n,j,l))
            vq = .25*(va(i,j)+va(i+1,j))*(vort(i,j-1)+vort(i,j))
            uo(i,j,l) = uo(i,j,l) + dt*(
     &           pgfx(i,j) +(ke(i,j)-ke(i+1,j))*bydx
     &           +vodh(i,j,l)*corofj + vq)
            pgf4pt=.25*(pgfy(i,j)+pgfy(i+1,j)+pgfy(i,j-1)+pgfy(i+1,j-1))
            vod(i,j,l) = vod(i,j,l) + dt*(pgf4pt - uoh(i,j,l)*corofj)
          enddo
        enddo
        i = im
        if(l <= lmou(i,j)) then
          vq = .25*(va(i,j)+va(1,j))*(vort(i,j-1)+vort(i,j))
          uo(i,j,l) = uo(i,j,l) + dt*(
     &         pgfx(i,j) +(ke(i,j)-ke(1,j))*bydx
     &         +vodh(i,j,l)*corofj + vq)
          pgf4pt = .25*(pgfy(i,j)+pgfy(1,j)+pgfy(i,j-1)+pgfy(1,j-1))
          vod(i,j,l) = vod(i,j,l) + dt*(pgf4pt - uoh(i,j,l)*corofj)
        endif
      enddo

c
c update VO, UOD
c

      do j=j1h,jnp
        bydy = 1d0/dyvo(j)
        if(j == jm-1) bydy = bydy*2./3.
        corofj = 2.*omega*sinvo(j)
        mvfac = .5*dxvo(j)
        pgfac = .25
        if(j == jm-1) pgfac = .5
        i = 1
        if(l <= lmov(i,j)) then
          mmid = (moh(i,j,l)+moh(i,j+1,l))
          mv(i,j) = mvfac*voh(i,j,l)*mmid
          smv(i,j,l) = smv(i,j,l) + mv(i,j)*xeven
          uq = .25*(ua(i,j)+ua(i,j+1))*(vort(im,j)+vort(i,j))
          vo(i,j,l) = vo(i,j,l) + dt*(
     &         pgfy(i,j) +(ke(i,j)-ke(i,j+1))*bydy
     &         -uodh(i,j,l)*corofj -uq)
          pgf4pt = pgfac*
     &         (pgfx(im,j)+pgfx(i,j)+pgfx(im,j+1)+pgfx(i,j+1))
          uod(i,j,l) = uod(i,j,l) + dt*(pgf4pt + voh(i,j,l)*corofj)
        endif
        do n=1,nbyzv(j,l)
          do i=max(2,i1yzv(n,j,l)),i2yzv(n,j,l)
            mmid = (moh(i,j,l)+moh(i,j+1,l))
            mv(i,j) = mvfac*voh(i,j,l)*mmid
            smv(i,j,l) = smv(i,j,l) + mv(i,j)*xeven
            uq = .25*(ua(i,j)+ua(i,j+1))*(vort(i-1,j)+vort(i,j))
            vo(i,j,l) = vo(i,j,l) + dt*(
     &           pgfy(i,j) +(ke(i,j)-ke(i,j+1))*bydy
     &           -uodh(i,j,l)*corofj -uq)
            pgf4pt = pgfac*
     &           (pgfx(i-1,j)+pgfx(i,j)+pgfx(i-1,j+1)+pgfx(i,j+1))
            uod(i,j,l) = uod(i,j,l) + dt*(pgf4pt + voh(i,j,l)*corofj)
          enddo
        enddo
      enddo

c
c update polar velocities
c
      call polevel(uo(1,j1a,l),vo(1,j1a,l),l)


c
c update MO
c
      Do J=J1P,JNP
        convfac = dt/dxypo(j)
        i = 1
        if(l <= lmom(i,j)) then
          convij = convfac*(MU(IM ,J)-MU(I,J) + MV(I,J-1)-MV(I,J))
          mo(i,j,l) = mo(i,j,l) + convij
          opbot(i,j) = opbot(i,j) + convij*grav
          if(mo(i,j,l)*dxypo(j) < 0.75d0*mmi(i,j,l)) then
            write(6,*) 'small mo',i,j,l,
     &           mo(i,j,l)*dxypo(j)/mmi(i,j,l)
            smallmo_loc = 1
          endif
        endif
        do n=1,nbyzm(j,l)
          do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
            convij = convfac*(MU(I-1,J)-MU(I,J) +MV(I,J-1)-MV(I,J))
            mo(i,j,l) = mo(i,j,l) + convij
            opbot(i,j) = opbot(i,j) + convij*grav
            if(mo(i,j,l)*dxypo(j) < 0.75d0*mmi(i,j,l)) then
              write(6,*) 'small mo',i,j,l,
     &             mo(i,j,l)*dxypo(j)/mmi(i,j,l)
              smallmo_loc = 1
            endif
          enddo
        enddo
      enddo
      If (QNP)  then
        j = jm
        convij = dt*sum(mv(:,jm-1))/(im*dxypo(jm))
        mo(1,j,l) = mo(1,j,l) + convij
        opbot(1,j) = opbot(1,j) + convij*grav
        do i=2,im
          mo(i,j,l) = mo(1,j,l)
        enddo
      endif

      enddo ! end loop over layers

      CALL GLOBALMAX(grid, smallmo_loc, smallmo)
      if(smallmo .gt. 0) call stop_model('small mo',255)

c emergency vertical regrid
c      Call OFLUXV(opbot,mo,qeven)

c fill the halos of just-updated quantities
      call halo_update(grid,mo)
      call halo_update(grid,opbot)
      call halo_update(grid,uo)
      call halo_update(grid,vo)

      Return
      End Subroutine ODHORZ

      subroutine polevel(u,v,l)
      Use OCEAN, Only: IM,JM,COSIC,SINIC, SINU,COSU
      USE OCEAN, only :
     &     nbyzm,nbyzu,nbyzv,
     &     i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      USE OCEANR_DIM, only : grid=>ogrid
      use DOMAIN_DECOMP_1D, only: hasNorthPole
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: u,v
      integer, intent(in) :: l
      integer :: i,j,n
      real*8 :: unp,vnp
c
c calculate U and V at the north pole
c
      if(hasNorthPole(grid)) then
        unp = 0.
        vnp = 0.
        j = jm-1
        do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            unp = unp - sinic(i)*v(i,j)
            vnp = vnp + cosic(i)*v(i,j)
          enddo
        enddo
        unp = unp*2/im
        vnp = vnp*2/im
c        vonp(l) = vnp
        do i=1,im
          u(i,jm) = unp*cosu(i)  + vnp*sinu(i)
          v(i,jm) = vnp*cosic(i) - unp*sinic(i)
        enddo
      endif
      return
      end subroutine polevel

      Subroutine ODHORZ0
C****
C**** OPGF0 calculates GUP, GDN, SUP and SDN inside each ocean cell,
C**** which will be kept constant during the dynamics of momentum
C****
      use constant, only : grav
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, MO,UO,VO,OPBOT,OPRESS,
     *                 G0M,GZM=>GZMO, S0M,SZM=>SZMO, DXYPO
      USE OCEAN, only : nbyzm,i1yzm,i2yzm
      Use OCEAN_DYN, Only: MMI, GUP,GDN, SUP,SDN, VBAR, dZGdP
      use ocean_dyn, only : dh3d=>dh ! for interface with other ocean routines
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
      Real*8,Parameter :: z12eH=.28867513d0  !  z12eH = 1/SQRT(12)
      Integer*4 I,J,L, J1,J1H,JN,J1A,JNH, N
      Logical :: QNP
      Real*8,External   :: VOLGSP
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     *          P
      Real*8 :: PUP,PDN, VUP,VDN
C****
C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      J1H = Max(J1-1,1)     !    1      4     JM-4   Halo minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      JNH = Min(JN+1,JM)    !    5      9     JM     Halo maximum
      QNP = JN==JM          !    F      F      T
      j1a = grid%j_strt_halo

C****
      Call HALO_UPDATE (GRID, G0M)
      Call HALO_UPDATE (GRID, GZM)
      Call HALO_UPDATE (GRID, S0M)
      Call HALO_UPDATE (GRID, SZM)

C**** Calculate pressure by integrating from the top down
      Call HALO_UPDATE (GRID,OPRESS)
      l=1
      do j=j1h,jnh
        do n=1,nbyzm(j,l)
          do i=i1yzm(n,j,l),i2yzm(n,j,l)
            opbot(i,j) = OPRESS(I,J)
          enddo
        enddo
      enddo
      do l=1,lmo
      do j=j1h,jnh
        do n=1,nbyzm(j,l)
          do i=i1yzm(n,j,l),i2yzm(n,j,l)
            P(I,J,L) = opbot(i,j)  + MO(I,J,L)*GRAV*.5
            opbot(i,j) = opbot(i,j) + MO(I,J,L)*GRAV
          enddo
        enddo
      enddo
      enddo

C****
      do l=lmo,1,-1
        Do J=J1H,JNH
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              MMI(I,J,L) = MO(I,J,L)*DXYPO(J)
              GUP(I,J,L) = (G0M(I,J,L) - 2*z12eH*GZM(I,J,L))/MMI(I,J,L)
              GDN(I,J,L) = (G0M(I,J,L) + 2*z12eH*GZM(I,J,L))/MMI(I,J,L)
              SUP(I,J,L) = (S0M(I,J,L) - 2*z12eH*SZM(I,J,L))/MMI(I,J,L)
              SDN(I,J,L) = (S0M(I,J,L) + 2*z12eH*SZM(I,J,L))/MMI(I,J,L)
              PUP = P(I,J,L) - MO(I,J,L)*GRAV*z12eH
              PDN = P(I,J,L) + MO(I,J,L)*GRAV*z12eH
              VUP = VOLGSP (GUP(I,J,L),SUP(I,J,L),PUP)
              VDN = VOLGSP (GDN(I,J,L),SDN(I,J,L),PDN)
              dZGdP(I,J,L) = VUP*(.5-z12eH) + VDN*(.5+z12eH)
              VBAR(I,J,L) = (VUP + VDN)*.5
              DH3D(I,J,L) = MO(I,J,L)*VBAR(I,J,L)
            enddo
          enddo
        enddo

C**** Copy to all longitudes at north pole
        If (QNP .and. l <= lmom(1,jm)) Then
          DH3D(2:IM,JM,L) = DH3D(1,JM,L)
          VBAR(2:IM,JM,L) = VBAR(1,JM,L)
          dZGdP(2:IM,JM,L) = dZGdP(1,JM,L)
          MO(2:IM,JM,L) = MO(1,JM,L)
        EndIf

c initialize polar velocities
        call polevel(uo(1,j1a,l),vo(1,j1a,l),l)

      enddo

      Return
      End Subroutine ODHORZ0

      SUBROUTINE OADVT2 (RM,RX,RY,RZ,DT,QLIMIT, OIJL)
!@sum  OADVT advects tracers using the linear upstream scheme.
C****
C**** Input:  MB (kg) = mass before advection
C****          DT (s) = time step
C****       MU (kg/s) = west to east mass flux
C****       MV (kg/s) = south to north mass flux
C****       MW (kg/s) = downward vertical mass flux
C****          QLIMIT = whether slope limitations should be used
C**** Output: RM (kg) = tracer mass
C****   RX,RY,RZ (kg) = first moments of tracer mass
C****       OIJL (kg) = diagnostic accumulation of tracer mass flux
C****
      USE OCEAN, only : im,jm,lmo
      USE OCEAN_DYN, only : mb=>mmi,smu,smv,smw

!      use domain_decomp_1d, only : grid, get
      use domain_decomp_1d, only : get
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      REAL*8, INTENT(INOUT),     DIMENSION
     *     (IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: RM,RX,RY,RZ
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO,3) :: OIJL
      REAL*8,
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: MA
      INTEGER I,J,L, J_0H
      LOGICAL, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT

      logical :: HAVE_NORTH_POLE  ! ,HAVE_SOUTH_POLE
         call get (grid, J_STRT_HALO=J_0H)
         call get (grid, HAVE_NORTH_POLE=HAVE_NORTH_POLE)
C        call get (grid, HAVE_SOUTH_POLE=HAVE_SOUTH_POLE)

C****
C**** Load mass after advection from mass before advection
C****
      MA = MB
C****
C**** Advect the tracer using the Slopes Scheme
C****
      CALL OADVTX2(RM,RX,RY,RZ,MA,SMU,5d-1*DT,QLIMIT,OIJL(1,J_0H,1,1))
      CALL OADVTY2(RM,RX,RY,RZ,MA,SMV,     DT,QLIMIT,OIJL(1,J_0H,1,2))
      CALL OADVTZ2(RM,RX,RY,RZ,MA,SMW,     DT,QLIMIT,OIJL(1,J_0H,1,3))
      CALL OADVTX2(RM,RX,RY,RZ,MA,SMU,5d-1*DT,QLIMIT,OIJL(1,J_0H,1,1))

C**** Fill in values at the poles
      if (HAVE_NORTH_POLE) then
      DO 20 L=1,LMO
      DO 20 I=1,IM
      RM(I,JM,L) = RM(1,JM,L)
   20 RZ(I,JM,L) = RZ(1,JM,L)
      end if
C     if (HAVE_SOUTH_POLE) then
C     DO 30 L=1,LMO
C     DO 30 I=1,IM
C     RM(I, 1,L) = RM(IM,1,L)
C  30 RZ(I, 1,L) = RZ(IM,1,L)
C     end if

      RETURN
      END SUBROUTINE OADVT2

      Subroutine OADVTX2(RM,RX,RY,RZ,MM,MU,DT,QLIMIT,OIJL)
C****
!@sum   OADVTX advects tracer in X direction via Linear Upstream Scheme
C****
C**** If QLIMIT is true, gradients are limited to prevent mean tracer
C**** from becoming negative; Abs[A(I)] must not exceed 1.
C****
C**** Input: DT (s) = time step
C****     MU (kg/s) = west to east mass flux
C****        QLIMIT = whether slope limitations should be used
C**** Input and Output: RM (kg) = tracer mass
C****             RX,RY,RZ (kg) = first moments of tracer mass
C****                   MM (kg) = ocean mass
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, LMOU=>LMU
      USE OCEAN, only : nbyzu,i1yzu,i2yzu, nbyzm,i1yzm,i2yzm
      Use OCEANR_DIM, Only: oGRID
      Implicit None
C**** Interface variables
      Logical,Intent(In) :: QLIMIT
      Real*8,   Intent(In) :: DT
      Real*8,Intent(In),
     *  Dimension(IM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO, LMO) ::
     *  MU
      Real*8,Intent(InOut),
     *  Dimension(IM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO, LMO) ::
     *  RM,RX,RY,RZ, MM,OIJL
C**** Local variables
      Integer :: I,J,L, J1P,JNP, N,NC,NCOURANT
      real*8, dimension(im) :: mudt
      Real*8 :: AM,A,FM,FX,FY,FZ,MMnew, zCOURANT,mcheck,courmax,rxlimit
      Real*8 :: AMim1,FMim1,FXim1,FYim1,FZim1,
     &          AM_im,FM_im,FX_im,FY_im,FZ_im

C**** Extract domain decomposition band parameters
      J1P = Max (oGRID%J_STRT, 2)     !  Exclude south pole
      JNP = Min (oGRID%J_STOP, JM-1)  !  Exclude north pole

      if(qlimit) then
        rxlimit = 1d0
      else
        rxlimit = 0d0
      endif

C****
C**** Loop over layers and latitudes
C****
      Do L=1,LMO
      Do J=J1P,JNP
        if(nbyzu(j,l) == 0) cycle


c
c adjust mass fluxes so that courant numbers are never greater than 1
c
        courmax = 0.
        mudt(1:2) = mu(1:2,j,l)*dt
        mudt(im) = mu(im,j,l)*dt
        i=1
        if(l <= lmou(i,j)) then
          if(mudt(i).ge.0.) then
            mcheck = mm(i,j,l)+min(0d0,mudt(im)-mudt(i))
          else
            mcheck = -mm(i+1,j,l)-min(0d0,mudt(i)-mudt(i+1))
          endif
          courmax = max(courmax,mudt(i)/mcheck)
        endif
        do n=1,nbyzu(j,l)
          i = i1yzu(n,j,l)
          if(i.gt.1) mudt(i-1) = 0. ! zero mudt at west boundary
          mudt(i) = mu(i,j,l)*dt
          do i=max(2,i1yzu(n,j,l)),min(i2yzu(n,j,l),im-1)
            mudt(i+1) = mu(i+1,j,l)*dt
            if(mudt(i).ge.0.) then
              mcheck = mm(i,j,l)+min(0d0,mudt(i-1)-mudt(i))
            else
              mcheck = -mm(i+1,j,l)-min(0d0,mudt(i)-mudt(i+1))
            endif
            courmax = max(courmax,mudt(i)/mcheck)
          enddo
        enddo
        i = im
        if(l <= lmou(i,j)) then
          if(mudt(i).ge.0.) then
            mcheck = mm(i,j,l)+min(0d0,mudt(i-1)-mudt(i))
          else
            mcheck = -mm(1,j,l)-min(0d0,mudt(i)-mudt(1))
          endif
          courmax = max(courmax,mudt(i)/mcheck)
        endif
        if(courmax.gt.1.) then
          ncourant = 1+int(courmax)
          write(6,*) 'ncourant ',j,l,ncourant
          zcourant = 1d0/ncourant
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),i2yzu(n,j,l)
              mudt(i) = mudt(i)*zcourant
            enddo
          enddo
        else
          ncourant = 1
        endif

        do nc=1,ncourant

c
c calculate fluxes at the dateline
c
          i = im
          if(l <= lmou(i,j)) then
            AM = mudt(I)
            if(am.ge.0.) then   ! mass flux is positive or zero
              A  = AM / MM(I,J,L)
              rx(i,j,l) = rx(i,j,l)-rxlimit*sign(min(0d0,
     &             rm(i,j,l)-abs(rx(i,j,l))),rx(i,j,l))
              FM = A * (RM(I,J,L) + (1-A)*RX(I,J,L))
              FX = AM * (A*A*RX(I,J,L) - 3*FM)
              FY = A*RY(I,J,L)
              FZ = A*RZ(I,J,L)
            else                ! mass flux is negative
              A  = AM / MM(1,J,L)
              rx(1,j,l) = rx(1,j,l)-rxlimit*sign(min(0d0,
     &             rm(1,j,l)-abs(rx(1,j,l))),rx(1,j,l))
              FM = A * (RM(1,J,L) - (1+A)*RX(1,J,L))
              FX = AM * (A*A*RX(1,J,L) - 3*FM)
              FY = A*RY(1,J,L)
              FZ = A*RZ(1,J,L)
            endif
          else
            am = 0.
            fm = 0.
            fx = 0.
            fy = 0.
            fz = 0.
          endif
          amim1 = am
          fmim1 = fm
          fxim1 = fx
          fyim1 = fy
          fzim1 = fz
          am_im = am
          fm_im = fm
          fx_im = fx
          fy_im = fy
          fz_im = fz

c
c loop over basins
c
          do n=1,nbyzm(j,l)
            if(i1yzm(n,j,l)>1 .and. i1yzm(n,j,l)==i2yzm(n,j,l)) cycle
            do i=i1yzm(n,j,l),min(i2yzm(n,j,l),im-1)
              AM = mudt(I)
              if(am.ge.0.) then ! mass flux is positive or zero
                A  = AM / MM(I,J,L)
                rx(i,j,l) = rx(i,j,l)-rxlimit*sign(min(0d0,
     &               rm(i,j,l)-abs(rx(i,j,l))),rx(i,j,l))
                FM = A * (RM(I,J,L) + (1-A)*RX(I,J,L))
                FX = AM * (A*A*RX(I,J,L) - 3*FM)
                FY = A*RY(I,J,L)
                FZ = A*RZ(I,J,L)
              else              ! mass flux is negative
                A  = AM / MM(I+1,J,L)
                rx(i+1,j,l) = rx(i+1,j,l)-rxlimit*sign(min(0d0,
     &               rm(i+1,j,l)-abs(rx(i+1,j,l))),rx(i+1,j,l))
                FM = A * (RM(I+1,J,L) - (1+A)*RX(I+1,J,L))
                FX = AM * (A*A*RX(I+1,J,L) - 3*FM)
                FY = A*RY(I+1,J,L)
                FZ = A*RZ(I+1,J,L)
              endif
              MMnew = MM(I,J,L) + (AMim1-AM)
              RM(I,J,L) = RM(I,J,L) +  (FMim1-FM)
              RX(I,J,L) = (RX(I,J,L)*MM(I,J,L) + (FXim1-FX) +
     &             3d0*((AMim1+AM)*RM(I,J,L)-MM(I,J,L)*(FMim1+FM)))
     &             / MMnew
              RY(I,J,L) = RY(I,J,L) + (FYim1-FY)
              RZ(I,J,L) = RZ(I,J,L) + (FZim1-FZ)
              MM(I,J,L) = MMnew
              OIJL(I,J,L) = OIJL(I,J,L) + FM
              amim1 = am
              fmim1 = fm
              fxim1 = fx
              fyim1 = fy
              fzim1 = fz
            enddo
          enddo

c
c update at dateline
c
          i = im
          if(l <= lmom(i,j)) then
            am = am_im
            fm = fm_im
            fx = fx_im
            fy = fy_im
            fz = fz_im
            MMnew = MM(I,J,L) + (AMim1-AM)
            RM(I,J,L) = RM(I,J,L) +  (FMim1-FM)
            RX(I,J,L) = (RX(I,J,L)*MM(I,J,L) + (FXim1-FX) +
     &           3d0*((AMim1+AM)*RM(I,J,L)-MM(I,J,L)*(FMim1+FM)))
     &           / MMnew
            RY(I,J,L) = RY(I,J,L) + (FYim1-FY)
            RZ(I,J,L) = RZ(I,J,L) + (FZim1-FZ)
            MM(I,J,L) = MMnew
            OIJL(I,J,L) = OIJL(I,J,L) + FM
          endif
        enddo ! nc
      enddo ! j
      enddo ! l

      Return
      EndSubroutine OADVTX2

      Subroutine OADVTY2(RM,RX,RY,RZ,MO,MV,DT,QLIMIT,OIJL)
C****
!@sum   OADVTY advects tracer in Y direction via Linear Upstream Scheme
c****
c**** oadvty advects tracers in the south to north direction using the
c**** linear upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mv (kg/s) = north-south mo flux, positive northward
c****        qlimit = whether slope limitations should be used
c****
c**** input/output:
c****     rm     (kg) = tracer mass
c****   rx,ry,rz (kg) = 1st moments of tracer mass
c****     mo     (kg) = ocean mass
c****
C****

      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, LMOV=>LMV
      USE OCEAN, only : nbyzm,i1yzm,i2yzm, nbyzv,i1yzv,i2yzv
      Use OCEANR_DIM, Only: grid=>oGRID
      use DOMAIN_DECOMP_1D, only : halo_update, 
     &     hasSouthPole, hasNorthPole
      Implicit None
C**** Interface variables
      Logical,Intent(In) :: QLIMIT
      Real*8,   Intent(In) :: DT
      Real*8,Intent(In),
     *  Dimension(IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO, LMO) ::
     *  MV
      Real*8,Intent(InOut),
     *  Dimension(IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO, LMO) ::
     *  RM,RX,RY,RZ, MO,OIJL
C**** Local variables
      logical :: qnp
      Integer :: I,J,L, J1H,J1P,JNP,JN, N, imin,imax
      Real*8 :: BM,B,FM,FX,FY,FZ,Mnew,rylimit
      real*8, dimension(im) :: BMjm1,FMjm1,FXjm1,FYjm1,FZjm1

C**** Extract domain decomposition band parameters
      J1H = Max (GRID%J_STRT-1, 1)   !  Exclude south pole
      J1P = Max (GRID%J_STRT, 2)     !  Exclude south pole
      JNP = Min (GRID%J_STOP, JM-1)  !  Exclude north pole
      JN = GRID%J_STOP
      QNP = hasNorthPole(grid)

      if(qlimit) then
        rylimit = 1d0
      else
        rylimit = 0d0
      endif

      call halo_update(grid,mo)
      call halo_update(grid,rm)
      call halo_update(grid,rx)
      call halo_update(grid,ry)
      call halo_update(grid,rz)

      do l=1,lmo

        do i=1,im
          bmjm1(i) = 0.
          fmjm1(i) = 0.
          fxjm1(i) = 0.
          fyjm1(i) = 0.
          fzjm1(i) = 0.
        enddo

        if(qnp) then
          j = jm
          do i=2,im
            mo(i,j,l) = mo(1,j,l)
            rm(i,j,l) = rm(1,j,l)
            rz(i,j,l) = rz(1,j,l)
          enddo
          do i=1,im
            rx(i,j,l) = 0.
            ry(i,j,l) = 0.
          enddo
        endif

        j = j1h
        do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            bm = mv(i,j,l)*dt
            if(bm.ge.0.) then   ! mass flux is positive or zero
              B  = BM / MO(I,J,L)
              ry(i,j,l) = ry(i,j,l)-rylimit*sign(min(0d0,
     &             rm(i,j,l)-abs(ry(i,j,l))),ry(i,j,l))
              FM = B * (RM(I,J,L) + (1-B)*RY(I,J,L))
              FY = BM * (B*B*RY(I,J,L) - 3*FM)
              FX = B*RX(I,J,L)
              FZ = B*RZ(I,J,L)
            else                ! mass flux is negative
              B  = BM / MO(i,j+1,L)
              ry(i,j+1,l) = ry(i,j+1,l)-rylimit*sign(min(0d0,
     &             rm(i,j+1,l)-abs(ry(i,j+1,l))),ry(i,j+1,l))
              FM = B * (RM(i,j+1,L) - (1+B)*RY(i,j+1,L))
              FY = BM * (B*B*RY(i,j+1,L) - 3*FM)
              FX = B*RX(i,j+1,L)
              FZ = B*RZ(i,j+1,L)
            endif
            bmjm1(i) = bm
            fmjm1(i) = fm
            fxjm1(i) = fx
            fyjm1(i) = fy
            fzjm1(i) = fz
          enddo
        enddo

        do j=j1p,jn
          do n=1,nbyzm(j,l)
            imin=i1yzm(n,j,l)
            imax=i2yzm(n,j,l)
            if(j == jm) imax=im
            do i=imin,imax
              bm = mv(i,j,l)*dt
              if(bm.ge.0.) then ! mass flux is positive or zero
                B  = BM / MO(I,J,L)
                ry(i,j,l) = ry(i,j,l)-rylimit*sign(min(0d0,
     &               rm(i,j,l)-abs(ry(i,j,l))),ry(i,j,l))
                FM = B * (RM(I,J,L) + (1-B)*RY(I,J,L))
                FY = BM * (B*B*RY(I,J,L) - 3*FM)
                FX = B*RX(I,J,L)
                FZ = B*RZ(I,J,L)
              else              ! mass flux is negative
                B  = BM / MO(i,j+1,L)
                ry(i,j+1,l) = ry(i,j+1,l)-rylimit*sign(min(0d0,
     &               rm(i,j+1,l)-abs(ry(i,j+1,l))),ry(i,j+1,l))
                FM = B * (RM(i,j+1,L) - (1+B)*RY(i,j+1,L))
                FY = BM * (B*B*RY(i,j+1,L) - 3*FM)
                FX = B*RX(i,j+1,L)
                FZ = B*RZ(i,j+1,L)
              endif
              Mnew = MO(I,J,L) + (BMjm1(i)-BM)
              RM(I,J,L) = RM(I,J,L) +  (FMjm1(i)-FM)
              RY(I,J,L) = (RY(I,J,L)*MO(I,J,L) + (FYjm1(i)-FY) + 3d0*
     &             ((BMjm1(i)+BM)*RM(I,J,L)-MO(I,J,L)*(FMjm1(i)+FM)))
     &             / Mnew
              RX(I,J,L) = RX(I,J,L) + (FXjm1(i)-FX)
              RZ(I,J,L) = RZ(I,J,L) + (FZjm1(i)-FZ)
              MO(I,J,L) = Mnew
              OIJL(I,J,L) = OIJL(I,J,L) + FM
              bmjm1(i) = bm
              fmjm1(i) = fm
              fxjm1(i) = fx
              fyjm1(i) = fy
              fzjm1(i) = fz
            enddo
          enddo
        enddo

c
c average the pole
c
        if(qnp .and. l<=lmom(1,jm)) then
          j = jm
          mo(:,j,l) = sum(mo(:,j,l))/im
          rm(:,j,l) = sum(rm(:,j,l))/im
          rz(:,j,l) = sum(rz(:,j,l))/im
          do i=1,im
            rx(i,j,l) = 0.
            ry(i,j,l) = 0.
          enddo
        endif

      enddo ! l

      return
      end subroutine oadvty2

      SUBROUTINE OADVTZ2(RM,RX,RY,RZ,MO,MW,DT,QLIMIT,OIJL)
C****
C**** OADVTZ advects tracers in the vertical direction using the
C**** linear upstream scheme.  If QLIMIT is true, the gradients are
C**** limited to prevent the mean tracer from becoming negative.
C****
C****
C**** Input: DT (s) = time step
C****     MW (kg/s) = downward vertical mass flux
C**** Input and Output: RM (kg) = tracer mass
C****             RX,RY,RZ (kg) = first moments of tracer mass
C****                    M (kg) = ocean mass
C****
      USE OCEAN, only : im,jm,lmo,lmm
      USE OCEAN, only : nbyzm,i1yzm,i2yzm
      use domain_decomp_1d, only : get
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     *  RM,RX,RY,RZ, OIJL, MO
      REAL*8, INTENT(IN),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: MW
      REAL*8, INTENT(IN) :: DT
      LOGICAL, INTENT(IN) :: QLIMIT
c
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     CMUP,FMUP,FXUP,FYUP,FZUP,FMUP_ctr
      REAL*8 :: CM,C,FM,FX,FY,FZ,MNEW,rzlim,no_rzlim
      real*8 :: fm_ctr,r_edge,wtdn,rdn,rup,rm_lus
c Applied when qlimit is true, edgmax is the maximum allowed ratio of
c r_edge to gridbox-mean r.  edgmax = 2 can be used if there is never
c flow out of the bottom and top edges simultaneously.
      real*8, parameter :: edgmax=1.5
      INTEGER :: I,J,L,N,LDN

      INTEGER :: J_0,J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      if(qlimit) then
        rzlim = 1d0
      else
        rzlim = 0d0
      endif
      no_rzlim = 1d0-rzlim

      do j=j_0,j_1
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            cmup(i,j) = 0.
            fmup(i,j) = 0.
            fxup(i,j) = 0.
            fyup(i,j) = 0.
            fzup(i,j) = 0.
            fmup_ctr(i,j) = 0.
          enddo
        enddo
      enddo

      do l=1,lmo
        do j=j_0,j_1
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              CM = DT*MW(I,J,L)
              ldn = min(l+1,lmm(i,j))
              wtdn = mo(i,j,l)/(mo(i,j,l)+mo(i,j,ldn)) ! use dzo instead?
              rdn = rm(i,j,ldn)/mo(i,j,ldn)
              rup = rm(i,j,l  )/mo(i,j,l  )
              R_edge = wtdn*rdn+(1.-wtdn)*rup
              if(cm.ge.0.) then ! mass flux is downward or zero
                C  = CM/MO(I,J,L)
                IF(C.GT.1d0)  WRITE (6,*) 'C>1:',I,J,L,C,MO(I,J,L)
                rz(i,j,l) = rz(i,j,l)-rzlim*sign(min(0d0,
     &               rm(i,j,l)-abs(rz(i,j,l))),rz(i,j,l))
                FM = C*(RM(I,J,L)+(1d0-C)*RZ(I,J,L))
                FX = C*RX(I,J,L)
                FY = C*RY(I,J,L)
                FZ = CM*(C*C*RZ(I,J,L)-3d0*FM)
                R_edge = R_edge*no_rzlim+rzlim*min(R_edge,edgmax*rup)
                fm_ctr = cm*(c*rup + (1.-c)*r_edge)
              else              ! mass flux is upward
                C  = CM/MO(I,J,L+1)
                IF(C.LT.-1d0)  WRITE (6,*) 'C<-1:',I,J,L,C,MO(I,J,L+1)
                rz(i,j,l+1) = rz(i,j,l+1)-rzlim*sign(min(0d0,
     &               rm(i,j,l+1)-abs(rz(i,j,l+1))),rz(i,j,l+1))
                FM = C*(RM(I,J,L+1)-(1d0+C)*RZ(I,J,L+1))
                FX = C*RX(I,J,L+1)
                FY = C*RY(I,J,L+1)
                FZ = CM*(C*C*RZ(I,J,L+1)-3d0*FM)
                R_edge = R_edge*no_rzlim+rzlim*min(R_edge,edgmax*rdn)
                fm_ctr = cm*(-c*rdn + (1.+c)*r_edge)
              endif
#ifdef LUS_VERT_ADV
              fm_ctr = fm
#endif
              RM_lus = RM(I,J,L) + (FMUP(I,J)-FM)
              RM(I,J,L) = RM(I,J,L) + (FMUP_ctr(I,J)-FM_ctr)
              mnew = MO(I,J,L) + CMUP(I,J)-CM
              RZ(I,J,L) = (RZ(I,J,L)*MO(I,J,L) + (FZUP(I,J)-FZ) +3d0*
     &             ((CMUP(I,J)+CM)*RM_lus-MO(I,J,L)*(FMUP(I,J)+FM)))
     &             / mnew
              MO(I,J,L) = mnew
              cmup(i,j) = cm
              fmup(i,j) = fm
              fmup_ctr(i,j) = fm_ctr
              fzup(i,j) = fz
              RX(I,J,L) = RX(I,J,L) + (FXUP(I,J)-FX)
              fxup(i,j) = fx
              RY(I,J,L) = RY(I,J,L) + (FYUP(I,J)-FY)
              fyup(i,j) = fy
              OIJL(I,J,L) = OIJL(I,J,L) + FM_ctr
            enddo
          enddo
        enddo
      enddo

      return
      END SUBROUTINE OADVTZ2

      SUBROUTINE OADVUZ(R,M,MW,DT,jmin,jmax,nbyz,i1yz,i2yz)
c simplest upstream scheme, for vertical direction
      USE OCEAN, only : im,jm,lmo
      USE OCEAN, only :  nbyzmax
      use domain_decomp_1d, only : get
      USE OCEANR_DIM, only : grid=>ogrid
      IMPLICIT NONE
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: R,M
      REAL*8, INTENT(IN),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: MW
      REAL*8, INTENT(IN) :: DT
      INTEGER, INTENT(IN) :: jmin,jmax
      INTEGER, INTENT(IN),
     &     DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: nbyz
      INTEGER, INTENT(IN),
     &     DIMENSION(nbyzmax,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     &     i1yz,i2yz
c
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     CMUP,FMUP
      REAL*8 :: CM,FM,MNEW
      INTEGER :: I,J,L,N
      do j=jmin,jmax
        do n=1,nbyz(j,1)
          do i=i1yz(n,j,1),i2yz(n,j,1)
            cmup(i,j) = 0.
            fmup(i,j) = 0.
          enddo
        enddo
      enddo
      do l=1,lmo
        do j=jmin,jmax
          do n=1,nbyz(j,l)
            do i=i1yz(n,j,l),i2yz(n,j,l)
              CM = DT*MW(I,J,L)
              if(cm.ge.0.) then ! mass flux is downward or zero
                FM = CM*R(I,J,L)
              else              ! mass flux is upward
                FM = CM*R(I,J,L+1)
              endif
              mnew = M(I,J,L) + CMUP(I,J)-CM
              R(I,J,L) = (R(I,J,L)*M(I,J,L) + (FMUP(I,J)-FM))/mnew
              cmup(i,j) = cm
              fmup(i,j) = fm
              M(I,J,L) = mnew
            enddo
          enddo
        enddo
      enddo
      return
      END SUBROUTINE OADVUZ

      SUBROUTINE OSTRES2
!@sum OSTRES applies the atmospheric surface stress over open ocean
!@sum and the sea ice stress to the layer 1 ocean velocities
!@auth Gary Russell
!@ver  1.0

      USE OCEAN, only : IMO=>IM,JMO=>JM
     *     , IVNP, UO,VO, UOD,VOD, MO,DXYSO,DXYNO,DXYVO
     *     , LMU,LMV, COSIC,SINIC

      USE DOMAIN_DECOMP_1D, only : get, halo_update, north,south

      USE OCEANR_DIM, only : ogrid

      USE OFLUXES, only : oDMUA,oDMVA, oDMUI,oDMVI

      IMPLICIT NONE
      INTEGER I,J,IM1,IP1

C****
C**** All stress now defined for whole box, not just ocn or ice fraction
C**** FLUXCB  DMUA(1)  U momentum downward into open ocean (kg/m*s)
C****         DMVA(1)  V momentum downward into open ocean (kg/m*s)
C****         DMUA(2,JM,1)  polar atmo. mass slowed to zero (kg/m**2)
C****         DMUI     U momentum downward from sea ice (kg/m*s)
C****         DMVI     V momentum downward from sea ice (kg/m*s)

      integer :: J_0, J_1, J_0S, J_1S  ; logical :: have_north_pole

      call get (ogrid, J_STRT=J_0, J_STOP=J_1,
     *                 J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     *                 have_north_pole=have_north_pole)

C****
C**** Surface stress is applied to U component
C****
      DO J=J_0S,J_1S
      I=IMO
      DO IP1=1,IMO
        IF(LMU(I,J).gt.0.)  UO(I,J,1) = UO(I,J,1) +
     *       (oDMUA(I,J,1) + oDMUA(IP1,J,1) + 2d0*oDMUI(I,J)) /
     *       (  MO(I,J,1) +   MO(IP1,J,1))
        I=IP1
      END DO
      END DO
      if (have_north_pole) then
        UO(IMO ,JMO,1) = UO(IMO ,JMO,1) + oDMUA(1,JMO,1)/MO(1,JMO,1)
        UO(IVNP,JMO,1) = UO(IVNP,JMO,1) + oDMVA(1,JMO,1)/MO(1,JMO,1)
      end if
      call halo_update(ogrid, odmvi, from=south)
      DO J=J_0S,J_1S
        I=IMO
        DO IP1=1,IMO
          IF(LMU(I,J).gt.0.) THEN
            VOD(I,J,1) = VOD(I,J,1) + (
     *           oDMVA(I,J,1)+oDMVA(IP1,J,1)
     *        +.5*(oDMVI(I,J-1)+oDMVI(IP1,J-1)+oDMVI(I,J)+oDMVI(IP1,J))
     *           )/(MO(I,J,1)+MO(IP1,J,1))
          ENDIF
          I=IP1
        END DO
      END DO
C****
C**** Surface stress is applied to V component
C****
      call halo_update(ogrid, odmva, from=north)
      call halo_update(ogrid,    mo, from=north)
      DO J=J_0S,min(J_1S,JMO-2)
      DO I=1,IMO
        IF(LMV(I,J).GT.0.)  VO(I,J,1) = VO(I,J,1) +
     *       (oDMVA(I,J  ,1)*DXYNO(J) + oDMVA(I,J+1,1)*DXYSO(J+1)
     *      + oDMVI(I,J)*DXYVO(J))  !!  2d0*oDMVI(I,J)*DXYVO(J) - error
     * / (MO(I,J,1)*DXYNO(J) + MO(I,J+1,1)*DXYSO(J+1))
      END DO
      END DO
C**** Surface stress is applied to V component at the North Pole
      if (have_north_pole) then
      DO I=1,IMO
        VO(I,JMO-1,1) = VO(I,JMO-1,1) +
     *    (oDMVA(I,JMO-1,1)*DXYNO(JMO-1)+
     *    (oDMVA(1,JMO,1)*COSIC(I) - oDMUA(1,JMO,1)*SINIC(I))*DXYSO(JMO)
     *   + oDMVI(I,JMO-1)*DXYVO(JMO-1)) /
     *  (MO(I,JMO-1,1)*DXYNO(JMO-1) + MO(I,JMO,1)*DXYSO(JMO))
      END DO
      end if
      call halo_update(ogrid, odmua, from=north)
      call halo_update(ogrid, odmui, from=north)
      DO J=J_0S,min(J_1S,JMO-2)
        IM1=IMO
        DO I=1,IMO
          IF(LMV(I,J).GT.0.) THEN
            UOD(I,J,1) = UOD(I,J,1) + (
     *           oDMUA(I,J,1)+oDMUA(I,J+1,1)
     *     +.5*(oDMUI(IM1,J)+oDMUI(I,J)+oDMUI(IM1,J+1)+oDMUI(I,J+1))
     *           )/(MO(I,J,1)+MO(I,J+1,1))
          ENDIF
          IM1=I
        END DO
      END DO
      RETURN
      END SUBROUTINE OSTRES2

      Subroutine OBDRAG2
!@sum  OBDRAG exerts a drag on the Ocean Model's bottom layer
!@auth Gary Russell
!@ver  2010/01/08
      Use OCEAN, Only: IM,JM,LMO,IVNP,J1O, MO,UO,VO,UOD,VOD, LMU,LMV,
     *                 DTS, COSI=>COSIC,SINI=>SINIC
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, SOUTH,NORTH
      Use OCEANR_DIM,       Only: oGRID

      Implicit None
      REAL*8, PARAMETER :: BDRAGX=1d0, SDRAGX=1d-1
      Integer*4 I,IP1,J,L,J1,JN,JNP
      Real*8    WSQ

C**** Define decomposition band parameters
      J1 = oGRID%J_STRT
      JN = oGRID%J_STOP
      JNP = JN  ;  If(JN==JM) JNP = JM-1
      
      Call HALO_UPDATE (oGRID, MO, From=NORTH)
C****
C**** Reduce ocean current at east edges of cells
C**** UO = UO*(1-x/y)  is approximated by  UO*y/(y+x)  for stability
C**** 
      Do J=Max(J1O,J1),JNP
        I=IM
        DO IP1=1,IM
          IF(LMU(I,J) > 0) THEN
            L=LMU(I,J)
            WSQ = UO(I,J,L)**2 + VOD(I,J,L)**2 + 1d-20
            UO(I,J,L) = UO(I,J,L) * (MO(I,J,L)+MO(IP1,J,L)) /
     *           (MO(I,J,L)+MO(IP1,J,L) + DTS*BDRAGX*SQRT(WSQ)*2d0)
            VOD(I,J,L) = VOD(I,J,L) * (MO(I,J,L)+MO(IP1,J,L)) /
     *           (MO(I,J,L)+MO(IP1,J,L) + DTS*BDRAGX*SQRT(WSQ)*2d0)
          ENDIF
          I=IP1
        ENDDO
      ENDDO
C****
C**** Reduce ocean current at north edges of cells
C****
      Do J=Max(J1O,J1),JNP
        DO I=1,IM
          IF(LMV(I,J) > 0) THEN
            L=LMV(I,J)   
            WSQ = VO(I,J,L)**2 + UOD(I,J,L)**2 + 1d-20
            VO(I,J,L) = VO(I,J,L) * (MO(I,J,L)+MO(I,J+1,L)) /
     *           (MO(I,J,L)+MO(I,J+1,L) + DTS*BDRAGX*SQRT(WSQ)*2d0) 
            UOD(I,J,L) = UOD(I,J,L) * (MO(I,J,L)+MO(I,J+1,L)) /
     *           (MO(I,J,L)+MO(I,J+1,L) + DTS*BDRAGX*SQRT(WSQ)*2d0)
          ENDIF
        ENDDO
      ENDDO
      RETURN
C****
      END Subroutine OBDRAG2
