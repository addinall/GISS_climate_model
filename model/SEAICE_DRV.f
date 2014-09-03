#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif
!@sum  SEAICE_DRV contains drivers for SEAICE related routines
!@auth Gavin Schmidt
!@ver  1.0
!@cont PRECIP_SI,GROUND_SI

      SUBROUTINE PRECIP_SI
!@sum  PRECIP_SI driver for applying precipitation to sea ice fraction
!@auth Original Development team
!@ver  1.0
!@calls seaice:prec_si
      USE CONSTANT, only : teeny,grav,tf,bylhm
      USE MODEL_COM, only : im,jm,fland,itoice,itlkice,focean
     *     ,p,ptop
#ifdef SCM
      USE MODEL_COM, only : I_TARG,J_TARG
      USE SCMCOM, only : iu_scm_prt,SCM_SURFACE_FLAG,ATSKIN
#endif
      USE GEOM, only : imaxj,axyp,byaxyp
      USE FLUXES, only : runpsi,prec,eprec,srunpsi,gtemp,apress,fwsim
     *     ,gtempr,erunpsi
#ifdef TRACERS_WATER
     *     ,trprec,trunpsi,gtracer
#endif
      USE SEAICE, only : prec_si, ace1i, lmi,xsi,debug
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi,flag_dsws,pond_melt
#ifdef TRACERS_WATER
     *     ,ntm,trsi
      USE TRDIAG_COM, only: taijn=>taijn_loc, tij_icocflx
#endif
      USE DIAG_COM, only : aij=>aij_loc,jreg,ij_f0oi,j_imelt,j_smelt
     *     ,j_hmelt,ij_fwio,ij_htio,ij_stio,ij_sisnwf,ij_sitopmlt
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_ATM, only : GET, GLOBALSUM
      IMPLICIT NONE

      REAL*8, DIMENSION(LMI) :: HSIL,TSIL,SSIL
      REAL*8 SNOW,MSI2,PRCP,ENRGP,RUN0,POICE,SRUN0,ERUN0
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: TRSIL
      REAL*8, DIMENSION(NTM) :: TRUN0,TRPRCP
#endif
      INTEGER I,J,JR,ITYPE
      LOGICAL WETSNOW

C****
C**** Extract useful local domain parameters from "grid"
C****
      integer :: J_0, J_1, J_0H, J_1H ,I_0,I_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE


      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Initialize work array
      DO J=J_0, J_1
      DO I=I_0,IMAXJ(J)
        JR=JREG(I,J)
      POICE=    RSI(I,J) *(1.-FLAND(I,J))
      RUNPSI(I,J)=0
      SRUNPSI(I,J)=0
      ERUNPSI(I,J)=0
#ifdef TRACERS_WATER
      TRUNPSI(:,I,J)=0
#endif
      IF (POICE.gt.0) THEN

        IF (FOCEAN(I,J).gt.0) THEN
          ITYPE=ITOICE
        ELSE
          ITYPE=ITLKICE
        END IF
        PRCP=PREC(I,J)
        ENRGP=EPREC(I,J)      ! energy of precip
        SNOW=SNOWI(I,J)
        MSI2=MSI(I,J)
        HSIL(:) = HSI(:,I,J)      ! sea ice temperatures
        SSIL(:) = SSI(:,I,J)      ! sea ice salt
#ifdef TRACERS_WATER
        TRSIL(:,:)=TRSI(:,:,I,J)  ! sea ice tracers
        TRPRCP(:)=TRPREC(:,I,J)*BYAXYP(I,J)   ! tracer in precip
#endif

        AIJ(I,J,IJ_F0OI)=AIJ(I,J,IJ_F0OI)+ENRGP*POICE
        AIJ(I,J,IJ_SISNWF)=AIJ(I,J,IJ_SISNWF)-MIN(ENRGP*BYLHM,0d0)*POICE

C**** CALL SUBROUTINE FOR CALCULATION OF PRECIPITATION OVER SEA ICE

        CALL PREC_SI(SNOW,MSI2,HSIL,TSIL,SSIL,PRCP,ENRGP,RUN0,SRUN0
     *       ,ERUN0,
#ifdef TRACERS_WATER
     *       TRSIL,TRPRCP,TRUN0,
#endif
     *       WETSNOW)

        SNOWI(I,J)  =SNOW
        RUNPSI(I,J) =RUN0
        SRUNPSI(I,J)=SRUN0
        ERUNPSI(I,J)=ERUN0
        HSI(:,I,J)=HSIL(:)
        SSI(:,I,J)=SSIL(:)
#ifdef TRACERS_WATER
        TRSI(:,:,I,J)=TRSIL(:,:)
        TRUNPSI(:,I,J)=TRUN0(:) ! tracer in runoff
#endif
        FLAG_DSWS(I,J)=FLAG_DSWS(I,J).or.WETSNOW
C**** reset flag if there was fresh snow (i.e. prcp but no rain!)
        IF (.not. WETSNOW .and. PRCP.gt.0.) FLAG_DSWS(I,J)=.FALSE.
C**** pond_melt accumulation
        pond_melt(i,j)=pond_melt(i,j)+0.3d0*RUN0

C**** set gtemp array
        MSI(I,J)=MSI2
        GTEMP(1:2,2,I,J)=TSIL(1:2)
        GTEMPR(2,I,J)   =TSIL(1)+TF
#ifdef SCM
        if (I.eq.I_TARG.and.J.eq.J_TARG) then
           if (SCM_SURFACE_FLAG.ge.1) then
               GTEMP(1,2,I,J) = ATSKIN
               GTEMP(2,2,I,J) = ATSKIN
               GTEMPR(2,I,J) = ATSKIN + TF
           endif
        endif
#endif
#ifdef TRACERS_WATER
        GTRACER(:,2,I,J) = TRSIL(:,1)/(XSI(1)*(SNOW+ACE1I)-SSIL(1))
#endif
        FWSIM(I,J) = RSI(I,J)*(ACE1I+SNOW+MSI2-SUM(SSIL(1:LMI)))

C**** Accumulate diagnostics for ice fraction
       CALL INC_AJ(I,J,ITYPE,J_IMELT,RUN0 *POICE)
       CALL INC_AJ(I,J,ITYPE,J_SMELT,SRUN0*POICE)
       CALL INC_AJ(I,J,ITYPE,J_HMELT,ERUN0*POICE)
       AIJ(I,J,IJ_SITOPMLT)=AIJ(I,J,IJ_SITOPMLT)+RUN0*POICE
       IF (FOCEAN(I,J).gt.0) THEN
          AIJ(I,J,IJ_FWIO)=AIJ(I,J,IJ_FWIO)+(RUN0-SRUN0)*POICE
          AIJ(I,J,IJ_HTIO)=AIJ(I,J,IJ_HTIO)+ERUN0*POICE
          AIJ(I,J,IJ_STIO)=AIJ(I,J,IJ_STIO)+SRUN0*POICE
#ifdef TRACERS_WATER
          TAIJN(I,J,TIJ_ICOCFLX,:)=TAIJN(I,J,TIJ_ICOCFLX,:)
     *                    +TRUN0(:)*POICE
#endif
       END IF

C**** Accumulate regional diagnostics
        CALL INC_AREG(I,J,JR,J_IMELT,RUN0 *POICE)
        CALL INC_AREG(I,J,JR,J_SMELT,SRUN0*POICE)
        CALL INC_AREG(I,J,JR,J_HMELT,ERUN0*POICE)

      END IF

C**** Calculate pressure anomaly at surface
      APRESS(I,J) = 100.*(P(I,J)+PTOP-1013.25d0)+
     *     RSI(I,J)*(SNOWI(I,J)+ACE1I+MSI(I,J))*GRAV

      END DO
      END DO

      IF (HAVE_SOUTH_POLE) APRESS(2:IM,1)  = APRESS(1,1)
      IF (HAVE_NORTH_POLE) APRESS(2:IM,JM) = APRESS(1,JM)
C****
      END SUBROUTINE PRECIP_SI

      SUBROUTINE UNDERICE
!@sum  underice calculates basal fluxes under sea and lake ice
!@+    saves the resulting fluxes
!@auth Gavin Schmidt
!@ver  1.0
!@calls iceocean,icelake
      USE CONSTANT, only : rhow,rhows,omega,rhoi,shw
      USE MODEL_COM, only : im,jm,focean,dtsrc,qcheck,kocean
      USE GEOM, only : sinlat2d,imaxj,axyp
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm, trname
#endif
      USE SEAICE, only : lmi,xsi,icelake,iceocean,ac2oim,alpha
     *     ,tfrez,debug,Ti,dEidTi,alami
      USE SEAICE_COM, only : msi,hsi,ssi,rsi
#ifdef TRACERS_WATER
     *     ,trsi
#endif
      USE FLUXES, only : fmsi_io,fhsi_io,fssi_io,ui2rho,gtemp,sss,mlhc
#ifdef TRACERS_WATER
     *     ,ftrsi_io,gtracer
#endif
      USE LAKES_COM, only : tlake,mwl,flake,gml,mldlk
#ifdef TRACERS_WATER
     *     ,trlake
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_ATM, only : GET
      USE TimerPackage_mod, only: startTimer => start
      USE TimerPackage_mod, only: stopTimer => stop
      IMPLICIT NONE
      INTEGER I,J
      REAL*8 coriol,ustar,Tm,Sm,Si,Tic,dh,mflux,hflux,sflux,fluxlim
     *     ,mlsh   !,mfluxmax
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM) :: Trm,Tri,trflux,tralpha
#ifdef TRACERS_SPECIAL_O18
      REAL*8 fracls
      INTEGER N
#endif
#endif
      integer :: J_0, J_1 ,I_0,I_1

      call startTimer('UNDERICE()')
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0, J_1
C**** Coriolis parameter (on tracer grid)
        DO I=I_0,IMAXJ(J)
          coriol = ABS(2.*OMEGA*SINLAT2D(I,J))
          IF ((FOCEAN(I,J)+FLAKE(I,J))*RSI(I,J).gt.0) THEN
C**** Set mixed layer conditions
            Tm = GTEMP(1,1,I,J)
#ifdef TRACERS_WATER
            Trm(:)=GTRACER(:,1,I,J)
            Tri(:)=TRSI(:,LMI,I,J)/(XSI(LMI)*MSI(I,J)-SSI(LMI,I,J))
#ifdef TRACERS_SPECIAL_O18
            do n=1,ntm
              tralpha(n)=fracls(n)
            end do
#else
            tralpha(:)=1.
#endif
#endif
            dh = 0.5*(XSI(LMI)*MSI(I,J))/RHOI
c            mfluxmax = (MSI(I,J)-AC2OIM)/dtsrc
            IF (FOCEAN(I,J).gt.0) THEN
C**** Ice lowest layer conditions
              Si = 1d3*SSI(LMI,I,J)/(XSI(LMI)*MSI(I,J))
              Tic = Ti(HSI(LMI,I,J)/(XSI(LMI)*MSI(I,J)),Si)
              IF (KOCEAN.ge.1) THEN
C**** should we calculate ocean rho(Tm,Sm) here?
                Ustar = MAX(5d-4,SQRT(UI2rho(I,J)/RHOWS))
                Sm = SSS(I,J)
                mlsh = MLHC(I,J)
                call iceocean(Tic,Si,Tm,Sm,dh,Ustar,Coriol,dtsrc,mlsh,
#ifdef TRACERS_WATER
     *               Tri,Trm,trflux,tralpha,
#endif
     *               mflux,sflux,hflux)
              ELSE ! for fixed SST assume freezing temp at base,implicit
                hflux=alami(Tic,Si)*(Tic-tfrez(sss(i,j)))/(dh+alpha
     *               *dtsrc*alami(Tic,Si)/(dEidTi(Tic,Si)*XSI(LMI)*MSI(I
     *               ,J)))
                mflux=0.
                sflux=0.
#ifdef TRACERS_WATER
                trflux= 0.
#endif
              END IF
            ELSE   ! for lakes (no salinity so solve directly)
              Tic =Ti(HSI(LMI,I,J)/(XSI(LMI)*MSI(I,J)),0d0)
              mlsh=SHW*MLDLK(I,J)*RHOW
              sflux = 0.
              call icelake(Tic,Tm,dh,dtsrc,mlsh,
#ifdef TRACERS_WATER
     *             Tri,Trm,trflux,tralpha,
#endif
     *             mflux,hflux)

C**** Limit lake-to-ice flux if lake is too shallow (< 40cm)
              IF (MWL(I,J).lt.0.4d0*RHOW*FLAKE(I,J)*AXYP(I,J)) THEN
                FLUXLIM=-GML(I,J)/(DTSRC*FLAKE(I,J)*AXYP(I,J))
                IF (hflux.lt.FLUXLIM) hflux = FLUXLIM
                if (mflux.lt.0) then
                  mflux = 0.
#ifdef TRACERS_WATER
                  trflux= 0.
#endif
                end if
                if (qcheck) print*,"Flux limiting",I,J,MWL(I,J)/
     *               (RHOW*FLAKE(I,J)*AXYP(I,J)),FLUXLIM*DTSRC
              END IF
            END IF
            FMSI_IO(I,J) = mflux*dtsrc   ! positive down
            FHSI_IO(I,J) = hflux*dtsrc
            FSSI_IO(I,J) = sflux*dtsrc
#ifdef TRACERS_WATER
            FTRSI_IO(:,I,J)=trflux(:)*dtsrc
#endif
          ELSE
            FMSI_IO(I,J) = 0.
            FHSI_IO(I,J) = 0.
            FSSI_IO(I,J) = 0.
#ifdef TRACERS_WATER
            FTRSI_IO(:,I,J)=0.
#endif
          END IF
        END DO
      END DO

C****
      call stopTimer('UNDERICE()')
      RETURN
      END SUBROUTINE UNDERICE

      SUBROUTINE MELT_SI
!@sum  MELT_SI driver for lateral melt of sea ice
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
!@calls SEAICE:SIMELT
      USE CONSTANT, only : sday,TF
      USE MODEL_COM, only : im,jm,kocean,focean,itoice,itlkice ! ,itime
     *     ,itocean,itlake,dtsrc                               ! ,nday
#ifdef SCM
      USE MODEL_COM, only : I_TARG,J_TARG
      USE SCMCOM, only : iu_scm_prt,SCM_SURFACE_FLAG,ATSKIN
#endif
      USE GEOM, only : axyp,imaxj
      USE DIAG_COM, only : j_imelt,j_hmelt,j_smelt,jreg,aij=>aij_loc,
     *     ij_fwio,ij_htio,ij_stio,ij_sigrlt
      USE SEAICE, only : simelt,tfrez,xsi,Ti,ace1i,debug
      USE SEAICE_COM, only : rsi,hsi,msi,lmi,snowi,ssi
#ifdef TRACERS_WATER
      USE SEAICE_COM, only : trsi,ntm
      USE TRDIAG_COM, only: taijn=>taijn_loc, tij_icocflx
#endif
      USE LAKES_COM, only : flake
      USE FLUXES, only : sss,melti,emelti,smelti,gtemp,gtempr,mlhc,fwsim
#ifdef TRACERS_WATER
     *     ,trmelti,gtracer
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_ATM, only : GET, GLOBALSUM
      USE TimerPackage_mod, only: startTimer => start
      USE TimerPackage_mod, only: stopTimer => stop
      IMPLICIT NONE
      REAL*8, DIMENSION(LMI) :: HSIL,TSIL,SSIL
      REAL*8 MSI2,ROICE,SNOW,ENRGUSED,RUN0,SALT,POCEAN,TFO
     *     ,PWATER,Tm,DT,ENRGMAX
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: TRSIL
      REAL*8, DIMENSION(NTM) :: TRUN0
#endif
      INTEGER I,J,ITYPE,JR,ITYPEO
      integer :: J_0, J_1 ,I_0,I_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      call startTimer('MELT_SI()')
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** CALCULATE LATERAL MELT (ALSO ELIMINATE SMALL AMOUNTS)
C**** EVERY PHYSICS TIME STEP
      DT=DTsrc
      DO J=J_0, J_1
        DO I=I_0,IMAXJ(J)
          PWATER=FOCEAN(I,J)+FLAKE(I,J)
          POCEAN=FOCEAN(I,J)
          RUN0=0. ; ENRGUSED=0. ; SALT=0.
#ifdef TRACERS_WATER
          TRUN0(:) = 0.
#endif
C**** Call simelt if (lake and v. small ice) or (q-flux ocean, some ice)
C**** now include lat melt for lakes and any RSI < 1
          IF ( (RSI(I,J).lt.1. .and. (FLAKE(I,J).gt.0 .and. RSI(I,J).gt
     *         .0)) .or. (KOCEAN.ge.1.and.POCEAN*RSI(I,J).gt.0) ) THEN
            JR=JREG(I,J)
            IF (POCEAN.gt.0) THEN
              ITYPE =ITOICE
              ITYPEO=ITOCEAN
              TFO = tfrez(sss(i,j))
            ELSE
              ITYPE =ITLKICE
              ITYPEO=ITLAKE
              TFO = 0.
            END IF
            Tm=GTEMP(1,1,I,J)
            ROICE=RSI(I,J)
            MSI2=MSI(I,J)
            SNOW=SNOWI(I,J)     ! snow mass
            HSIL(:)= HSI(:,I,J) ! sea ice enthalpy
            SSIL(:)= SSI(:,I,J) ! sea ice salt
            ENRGMAX= MAX(Tm-TFO,0d0)*MLHC(I,J) ! max energy for melt
#ifdef TRACERS_WATER
            TRSIL(:,:)=TRSI(:,:,I,J) ! tracer content of sea ice
#endif
            CALL SIMELT(DT,ROICE,SNOW,MSI2,HSIL,SSIL,POCEAN,Tm,TFO,TSIL
#ifdef TRACERS_WATER
     *           ,TRSIL,TRUN0
#endif
     *           ,ENRGMAX,ENRGUSED,RUN0,SALT)

C**** accumulate diagnostics
            AIJ(I,J,IJ_SIGRLT)=AIJ(I,J,IJ_SIGRLT)-RUN0*PWATER
            IF (FOCEAN(I,J).gt.0) THEN
              AIJ(I,J,IJ_FWIO)=AIJ(I,J,IJ_FWIO)+(RUN0-SALT)*PWATER
              AIJ(I,J,IJ_HTIO)=AIJ(I,J,IJ_HTIO)-ENRGUSED*PWATER
              AIJ(I,J,IJ_STIO)=AIJ(I,J,IJ_STIO)+SALT*PWATER
#ifdef TRACERS_WATER
              TAIJN(I,J,TIJ_ICOCFLX,:)=TAIJN(I,J,TIJ_ICOCFLX,:)
     *                                        +TRUN0(:)*PWATER
#endif
            END IF

           CALL INC_AJ(I,J,ITYPE,J_HMELT,-ENRGUSED*ROICE*PWATER)
           CALL INC_AJ(I,J,ITYPE,J_SMELT,     SALT*ROICE*PWATER)
           CALL INC_AJ(I,J,ITYPE,J_IMELT,     RUN0*ROICE*PWATER)
           CALL INC_AJ(I,J,ITYPEO,J_HMELT,-ENRGUSED*(1.-ROICE)*PWATER)
           CALL INC_AJ(I,J,ITYPEO,J_SMELT,     SALT*(1.-ROICE)*PWATER)
           CALL INC_AJ(I,J,ITYPEO,J_IMELT,     RUN0*(1.-ROICE)*PWATER)
           CALL INC_AREG(I,J,JR,J_HMELT,-ENRGUSED*PWATER)
           CALL INC_AREG(I,J,JR,J_SMELT,     SALT*PWATER)
           CALL INC_AREG(I,J,JR,J_IMELT,     RUN0*PWATER)
C**** Update prognostic sea ice variables + correction for rad. fluxes
            if (roice.gt.rsi(i,j)) ! ice from ocean
     *          call RESET_SURF_FLUXES(I,J,1,2,RSI(I,J),ROICE)
            if (roice.lt.rsi(i,j)) ! ocean from ice
     *          call RESET_SURF_FLUXES(I,J,2,1,1.-RSI(I,J),1.-ROICE)
C****
            RSI(I,J)=ROICE
            MSI(I,J)=MSI2
            SNOWI(I,J)=SNOW
            HSI(:,I,J)=HSIL(:)
            SSI(:,I,J)=SSIL(:)
#ifdef TRACERS_WATER
            TRSI(:,:,I,J)=TRSIL(:,:)
#endif
          END IF
C**** Save fluxes (in kg, J etc.), positive into ocean
          MELTI(I,J) = RUN0*PWATER*AXYP(I,J)
          EMELTI(I,J)=-ENRGUSED*PWATER*AXYP(I,J)
          SMELTI(I,J)= SALT*PWATER*AXYP(I,J)
#ifdef TRACERS_WATER
          TRMELTI(:,I,J)=TRUN0(:)*PWATER*AXYP(I,J)
#endif
C**** Reset some defaults if all ice is gone
          IF (RSI(I,J).eq.0) THEN
            GTEMP(1,2,I,J)=Ti(HSI(1,I,J)/(XSI(1)*ACE1I),1d3*SSI(1,I,J
     *           )/(XSI(1)*ACE1I))
            GTEMP(2,2,I,J)=Ti(HSI(2,I,J)/(XSI(2)*ACE1I),1d3*SSI(2,I,J
     *           )/(XSI(2)*ACE1I))
            GTEMPR(2,I,J) = GTEMP(1,2,I,J) + TF
#ifdef SCM
            if (I.eq.I_TARG.and.J.eq.J_TARG) then
              if (SCM_SURFACE_FLAG.ge.1) then
                  GTEMP(1,2,I,J) = ATSKIN
                  GTEMP(2,2,I,J) = ATSKIN
                  GTEMPR(2,I,J) = ATSKIN + TF
              endif
            endif
#endif
#ifdef TRACERS_WATER
            GTRACER(:,2,I,J) = 0.
#endif
          END IF
C****
        END DO
      END DO
C****
      call stopTimer('MELT_SI()')
      RETURN
      END SUBROUTINE MELT_SI

      SUBROUTINE GROUND_SI
!@sum  GROUND_SI driver for applying surface + base fluxes to sea ice
!@auth Gary Russell/Gavin Schmidt
!@ver  2010/11/12
!@calls SEAICE:SEA_ICE
      USE CONSTANT, only : grav,rhows,rhow,sday
      USE MODEL_COM, only : im,jm,dtsrc,fland,focean
     *     ,itoice,itlkice,p,ptop,jhour,jday
      USE GEOM, only : imaxj,axyp
      USE FLUXES, only : e0,e1,evapor,runosi,erunosi,srunosi,solar
     *     ,fmsi_io,fhsi_io,fssi_io,apress,gtemp,sss
#ifdef TRACERS_WATER
     *     ,ftrsi_io,trevapor,trunosi,gtracer
#ifdef TRACERS_DRYDEP
     *     ,trdrydep
#endif
#endif
      USE SEAICE, only : sea_ice,ssidec,lmi,xsi,ace1i,qsfix,debug
     *     ,snowice, snow_ice, rhos, Ti
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi,pond_melt,flag_dsws
#ifdef TRACERS_WATER
     *     ,trsi,ntm
      USE TRDIAG_COM, only: taijn=>taijn_loc, tij_icocflx
#endif
      USE LAKES_COM, only : mwl,gml,flake
      USE DIAG_COM, only : aij=>aij_loc,jreg,ij_rsoi,ij_msi
     *     ,j_imelt,j_hmelt,j_smelt,j_rsnow,ij_rsit,ij_rsnw,ij_snow
     *     ,ij_mltp,ij_zsnow,ij_fwio,ij_htio,ij_stio,ij_sntosi,ij_tsice
     *     ,ij_sihc,ij_sigrcg,ij_sitopmlt,ij_sibotmlt
     *     ,IJ_MSNFLOOD,IJ_HSNFLOOD
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_ATM, only : GET, GLOBALSUM
      USE TimerPackage_mod, only: startTimer => start
      USE TimerPackage_mod, only: stopTimer => stop
      IMPLICIT NONE

      REAL*8, DIMENSION(LMI) :: HSIL,SSIL
      REAL*8 SNOW,ROICE,MSI2,F0DT,F1DT,EVAP,SROX(2)
     *     ,FMOC,FHOC,FSOC,POICE,PWATER,SCOVI
      REAL*8 MFLUX,HFLUX,SFLUX,RUN,ERUN,SRUN,MELT12
      REAL*8 MSNWIC,HSNWIC,SSNWIC,SM,TM,Ti1,DSNOW
      INTEGER I,J,JR,ITYPE
      LOGICAL WETSNOW
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: trsil
      REAL*8, DIMENSION(NTM) :: trflux,ftroc,trevap,trrun,trsnwic,trm
     *     ,tralpha
#endif
      integer :: J_0, J_1, J_0H, J_1H ,I_0,I_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      call startTimer('GROUND_SI()')
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT = J_0,     J_STOP = J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      debug=.false.

C**** Initialize work array
      DO J=J_0, J_1
      DO I=I_0,IMAXJ(J)
c      debug=(i.eq.7.and.j.eq.42).or.(i.eq.30.and.j.eq.45).or.
c     *        (i.eq.1.and.j.eq.46)
c      debug=i.eq.40.and.j.eq.41
      PWATER=FOCEAN(I,J)+FLAKE(I,J)   ! 1.-FLAND(I,J)
      ROICE=RSI(I,J)
      POICE=ROICE*PWATER
      JR=JREG(I,J)
      SOLAR(3,I,J)=0
      RUNOSI(I,J)=0
      ERUNOSI(I,J)=0
      SRUNOSI(I,J)=0
#ifdef TRACERS_WATER
      TRUNOSI(:,I,J) = 0.
#endif
      IF (POICE.gt.0) THEN

        F0DT=E0(I,J,2) ! heat flux to the top ice surface (J/m^2)
        F1DT=E1(I,J,2) ! heat flux between 1st and 2nd ice layer (J/m^2)
        EVAP=EVAPOR(I,J,2) ! evaporation/dew at the ice surface (kg/m^2)
        SROX(1)=SOLAR(2,I,J) ! solar radiation absrbd by sea ice (J/m^2)
        FMOC=fmsi_io(i,j)  ! mass flux at base (kg/m^2)
        FSOC=fssi_io(i,j)  ! salt flux at base (kg/m^2)
        FHOC=fhsi_io(i,j)  ! heat flux at base (J/M^2)
        SNOW= SNOWI(I,J)  ! snow mass (kg/m^2)
        MSI2= MSI(I,J)
        HSIL(:) = HSI(:,I,J)  ! sea ice enthalpy
        SSIL(:) = SSI(:,I,J)  ! sea ice salt
        WETSNOW=FLAG_DSWS(I,J)  ! wetness of snow
        Tm=GTEMP(1,1,I,J)    ! ocean mixed layer temperature (C)
#ifdef TRACERS_WATER
        TREVAP(:) = TREVAPOR(:,2,I,J)
#ifdef TRACERS_DRYDEP
     *       -trdrydep(:,2,i,j)
#endif
        FTROC(:)  = ftrsi_io(:,i,j)
        TRSIL(:,:)= TRSI(:,:,I,J)
        Trm(:)=GTRACER(:,1,I,J)
        Tralpha(:)=1d0  ! no fractionation for snow ice formation
#endif
        IF (FOCEAN(I,J).gt.0) THEN
          ITYPE=ITOICE
          Sm=SSS(I,J)           ! ocean mixed layer salinity (psu)
        ELSE
          ITYPE=ITLKICE
          Sm=0.                 ! lakes always fresh
        END IF

        AIJ(I,J,IJ_RSOI) =AIJ(I,J,IJ_RSOI) +POICE
        AIJ(I,J,IJ_MSI) =AIJ(I,J,IJ_MSI) + (ACE1I+MSI2)*POICE

c        if (debug) write(6,'(A,2I4,4F11.6)') "si0",i,j,SNOW,1d3*SSIL(1)
c     $       /(XSI(1)*(ACE1I+SNOW)),1d3*SSIL(2)/(XSI(2)*(ACE1I+SNOW))
c     $       ,1d3*(SSIL(1)+SSIL(2))/ACE1I
c        if (debug) print*,"si0",i,j,SNOW,ROICE,HSIL,SSIL,MSI2,F0DT,F1DT,
c     *       EVAP,SROX,SNOW+ACE1I-SSIL(1)-SSIL(2) !, TRSIL(1,1)+TRSIL(1,2)
        CALL SEA_ICE(DTSRC,SNOW,ROICE,HSIL,SSIL,MSI2,F0DT,F1DT,EVAP,SROX
#ifdef TRACERS_WATER
     *       ,TRSIL,TREVAP,FTROC,TRRUN
#endif
     *       ,FMOC,FHOC,FSOC,RUN,ERUN,SRUN,WETSNOW,MELT12)
c        if (debug)  write(6,'(A,2I4,4F11.6)') "si1",i,j,SNOW,1d3*SSIL(1)
c     $       /(XSI(1)*(ACE1I+SNOW)),1d3*SSIL(2)/(XSI(2)*(ACE1I+SNOW))
c     $       ,1d3*(SSIL(1)+SSIL(2))/ACE1I
c        if (debug) print*,"si1",i,j,SNOW,HSIL,SSIL,MSI2,FMOC,FHOC,FSOC,
c     *       RUN,ERUN,SRUN,WETSNOW,MELT12,SNOW
c     $       +ACE1I-SSIL(1)-SSIL(2) !, TRSIL(1,1)+TRSIL(1,2)

C**** Decay sea ice salinity
        if (FOCEAN(I,J).gt.0) then
          CALL SSIDEC(SNOW,MSI2,HSIL,SSIL,DTsrc,MELT12,
#ifdef TRACERS_WATER
     *         TRSIL,TRFLUX,
#endif
     *         MFLUX,HFLUX,SFLUX)
        else
          MFLUX=0. ; SFLUX=0. ; HFLUX=0.
#ifdef TRACERS_WATER
          TRFLUX = 0.
#endif
        end if

c        if (debug)  write(6,'(A,2I4,4F11.6)') "si2",i,j,SNOW,1d3*SSIL(1)
c     $       /(XSI(1)*(ACE1I+SNOW)),1d3*SSIL(2)/(XSI(2)*(ACE1I+SNOW))
c     $       ,1d3*(SSIL(1)+SSIL(2))/ACE1I
c        if (debug) print*,"si2",i,j,SNOW,HSIL,SSIL,MSI2,MFLUX,HFLUX,
c     *       SFLUX,SNOW+ACE1I-SSIL(1)-SSIL(2) !, TRSIL(1,1)+TRSIL(1,2)
C**** Calculate snow-ice possibility
        if (snow_ice .eq. 1 .and. FOCEAN(I,J).gt.0) then
          call snowice(Tm,Sm,SNOW,MSI2,HSIL,SSIL,qsfix,
#ifdef TRACERS_WATER
     *         Trm,Tralpha,TRSIL,TRSNWIC,
#endif
     *         MSNWIC,HSNWIC,SSNWIC,DSNOW)
        else
          MSNWIC=0. ; SSNWIC=0. ; HSNWIC=0. ; DSNOW=0.
#ifdef TRACERS_WATER
          TRSNWIC = 0.
#endif
        end if
c        if (debug)  write(6,'(A,2I4,4F11.6)') "si3",i,j,SNOW,1d3*SSIL(1)
c     $       /(XSI(1)*(ACE1I+SNOW)),1d3*SSIL(2)/(XSI(2)*(ACE1I+SNOW))
c     $       ,1d3*(SSIL(1)+SSIL(2))/ACE1I
c        if (debug) print*,"si3",i,j,SNOW,HSIL,SSIL,MSI2, MSNWIC,HSNWIC,
c     *       SSNWIC,SNOW+ACE1I-SSIL(1)-SSIL(2) ! ,TRSIL(1,1)+TRSIL(1,2),

C**** RESAVE PROGNOSTIC QUANTITIES
        SNOWI(I,J)=SNOW
        HSI(:,I,J)=HSIL(:)
        SSI(:,I,J)=SSIL(:)
        MSI(I,J) = MSI2
#ifdef TRACERS_WATER
        TRSI(:,:,I,J) = TRSIL(:,:)
#endif
        FLAG_DSWS(I,J)=WETSNOW
        Ti1 = Ti(HSIL(1)/(XSI(1)*(SNOW+ACE1I)),1d3*SSIL(1)/(XSI(1)*(SNOW
     *       +ACE1I)))

C**** pond_melt accumulation
        pond_melt(i,j)=pond_melt(i,j)+0.3d0*MELT12
        pond_melt(i,j)=MIN(pond_melt(i,j),0.5*(MSI2+SNOW+ACE1I))

C**** decay is slow if there is some melting, faster otherwise
        if (MELT12.gt.0) then   ! 30 day decay
          pond_melt(i,j)=pond_melt(i,j)*(1.-dtsrc/(30.*sday))
        else                    ! 10 day decay
          pond_melt(i,j)=pond_melt(i,j)*(1.-dtsrc/(10.*sday))
        end if

C**** saftey valve to ensure that melt ponds eventually disappear (Ti<-10)
        if (Ti1 .lt.-10.) pond_melt(i,j)=0.  ! refreeze

c        if (debug) write(6,'(A,4I4,4F11.6)') "ponds",i,j,jday,jhour,
c     *       melt12,pond_melt(i,j),ti(HSIL(1)/(XSI(1)*(SNOW+ACE1I)),
c     *         1d3*SSIL(1)/(XSI(1)*(SNOW+ACE1I))),
c     *         1d3*(SSIL(1)+SSIL(2))/ACE1I

C**** Net fluxes to ocean
        RUNOSI(I,J) = FMOC + RUN  + MFLUX + MSNWIC
        ERUNOSI(I,J)= FHOC + ERUN + HFLUX + HSNWIC
        SRUNOSI(I,J)= FSOC + SRUN + SFLUX + SSNWIC
        SOLAR(3,I,J)= SROX(2)
#ifdef TRACERS_WATER
        TRUNOSI(:,I,J) = FTROC(:) + TRRUN(:) + TRFLUX(:) + TRSNWIC(:)
#endif

C**** ACCUMULATE DIAGNOSTICS
          SCOVI=0.
C**** snow cover diagnostic now matches that seen by the radiation
          IF (SNOW.GT.0) SCOVI=MIN(1d0,SNOW/(RHOS*0.1d0))*POICE

          AIJ(I,J,IJ_RSNW)=AIJ(I,J,IJ_RSNW)+SCOVI
          AIJ(I,J,IJ_SNOW)=AIJ(I,J,IJ_SNOW)+SNOW*POICE
          AIJ(I,J,IJ_RSIT)=AIJ(I,J,IJ_RSIT)+POICE
          AIJ(I,J,IJ_MLTP)=AIJ(I,J,IJ_MLTP)+pond_melt(i,j)*POICE
          AIJ(I,J,IJ_ZSNOW)=AIJ(I,J,IJ_ZSNOW)+POICE*SNOW/RHOS
          AIJ(I,J,IJ_SNTOSI)=AIJ(I,J,IJ_SNTOSI)+POICE*(DSNOW-MSNWIC)
          AIJ(I,J,IJ_TSICE)=AIJ(I,J,IJ_TSICE)+Ti1*POICE
          AIJ(I,J,IJ_SIHC)=AIJ(I,J,IJ_SIHC)+SUM(HSIL(:))*POICE
          AIJ(I,J,IJ_SITOPMLT)=AIJ(I,J,IJ_SITOPMLT)+POICE*(RUN+MFLUX)
          AIJ(I,J,IJ_MSNFLOOD) = AIJ(I,J,IJ_MSNFLOOD) - POICE*MSNWIC
          AIJ(I,J,IJ_HSNFLOOD) = AIJ(I,J,IJ_MSNFLOOD) - POICE*HSNWIC
          IF (FMOC.lt.0) THEN   ! define as congelation growth
            AIJ(I,J,IJ_SIGRCG)=AIJ(I,J,IJ_SIGRCG)-POICE*FMOC
          ELSE                  ! basal melt
            AIJ(I,J,IJ_SIBOTMLT)=AIJ(I,J,IJ_SIBOTMLT)+POICE*FMOC
          END IF

          IF (FOCEAN(I,J).gt.0) THEN
            AIJ(I,J,IJ_FWIO)=AIJ(I,J,IJ_FWIO)+(RUNOSI(I,J)-SRUNOSI(I,J))
     *           *POICE
            AIJ(I,J,IJ_HTIO)=AIJ(I,J,IJ_HTIO)+ERUNOSI(I,J)*POICE
            AIJ(I,J,IJ_STIO)=AIJ(I,J,IJ_STIO)+SRUNOSI(I,J)*POICE
#ifdef TRACERS_WATER
            TAIJN(I,J,TIJ_ICOCFLX,:)=TAIJN(I,J,TIJ_ICOCFLX,:)
     *                                       +TRUNOSI(:,I,J)*POICE
#endif
          END IF

          CALL INC_AJ(I,J,ITYPE,J_RSNOW,SCOVI)
          CALL INC_AJ(I,J,ITYPE,J_IMELT,(FMOC+ RUN+MFLUX+MSNWIC)*POICE)
          CALL INC_AJ(I,J,ITYPE,J_HMELT,(FHOC+ERUN+HFLUX+HSNWIC)*POICE)
          CALL INC_AJ(I,J,ITYPE,J_SMELT,(FSOC+SRUN+SFLUX+SSNWIC)*POICE)

          CALL INC_AREG(I,J,JR,J_RSNOW,SCOVI)
          CALL INC_AREG(I,J,JR,J_IMELT,(FMOC+ RUN+MFLUX+MSNWIC)*POICE)
          CALL INC_AREG(I,J,JR,J_HMELT,(FHOC+ERUN+HFLUX+HSNWIC)*POICE)
          CALL INC_AREG(I,J,JR,J_SMELT,(FSOC+SRUN+SFLUX+SSNWIC)*POICE)

      END IF
C**** set total atmopsheric pressure anomaly in case needed by ocean
      APRESS(I,J) = 100.*(P(I,J)+PTOP-1013.25d0)+RSI(I,J)
     *     *(SNOWI(I,J)+ACE1I+MSI(I,J))*GRAV

      END DO
      END DO

      IF (HAVE_SOUTH_POLE) APRESS(2:IM,1)  = APRESS(1,1)
      IF (HAVE_NORTH_POLE) APRESS(2:IM,JM) = APRESS(1,JM)
C****
      call stopTimer('GROUND_SI()')
      END SUBROUTINE GROUND_SI

      SUBROUTINE FORM_SI
!@sum  FORM_SI driver for adding new sea ice
!@auth Original Development team
!@ver  1.0
!@calls seaice:addice
      USE CONSTANT, only : tf
      USE MODEL_COM, only : im,jm,focean,kocean,fland
     *     ,itocean,itoice,itlake,itlkice,itime
#ifdef SCM
      USE MODEL_COM, only : I_TARG,J_TARG
      USE SCMCOM, only : iu_scm_prt,SCM_SURFACE_FLAG,ATSKIN
#endif
      USE GEOM, only : imaxj,axyp
#ifdef TRACERS_WATER
      USE TRACER_COM, only : itime_tr0,tr_wd_type,nWater,nPART
      USE TRDIAG_COM, only : taijn=>taijn_loc, tij_icocflx, tij_seaice
#endif
      USE DIAG_COM, only : aij=>aij_loc,jreg,j_rsi,j_ace1,j_ace2,j_snow
     *     ,j_smelt,j_imelt,j_hmelt,ij_tsi,ij_ssi1,ij_ssi2,j_implh
     *     ,j_implm,ij_smfx,ij_fwio,ij_htio,ij_stio,ij_sigrfr,ij_sigrcg
      USE SEAICE, only : ace1i,addice,lmi,fleadoc,fleadlk,xsi,debug
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi
#ifdef TRACERS_WATER
     *     ,trsi,ntm
#endif
      USE FLUXES, only : dmsi,dhsi,dssi,gtemp,fwsim,gtempr
#ifdef TRACERS_WATER
     *     ,dtrsi,gtracer
#endif
      USE LAKES_COM, only : flake
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_ATM, only : GET, GLOBALSUM
      IMPLICIT NONE

      REAL*8, DIMENSION(LMI) :: HSIL,TSIL,SSIL
      REAL*8 SNOW,ROICE,MSI2,ENRGFO,ACEFO,ACEFI,ENRGFI,SALTO,SALTI
     *     ,POICE,PWATER,FLEAD,POCEAN,DMIMP,DHIMP,DSIMP
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: trsil
      REAL*8, DIMENSION(NTM) :: tro,tri,dtrimp
#endif
      LOGICAL QFIXR
      INTEGER I,J,JR,ITYPE,ITYPEO,N

C****
C**** Extract useful local domain parmeters from "grid"
C****
      integer :: J_0, J_1, J_0H, J_1H ,I_0,I_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE   ,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE   )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      debug=.false.

C**** Initialize work array
      DO J=J_0, J_1
      DO I=I_0,IMAXJ(J)

c         debug=i.eq.40.and.j.eq.41

      PWATER=FOCEAN(I,J)+FLAKE(I,J)
      ROICE=RSI(I,J)
      POICE=ROICE*PWATER
      POCEAN=(1.-ROICE)*PWATER
      JR=JREG(I,J)
      IF (PWATER.gt.0) THEN

        SNOW= SNOWI(I,J)      ! snow mass (kg/m^2)
        MSI2= MSI(I,J)
        HSIL(:) = HSI(:,I,J)      ! sea ice enthalpy
        SSIL(:) = SSI(:,I,J)      ! sea ice salt
#ifdef TRACERS_WATER
        TRSIL(:,:)= TRSI(:,:,I,J)
#endif

        IF (FOCEAN(I,J).gt.0) THEN
          FLEAD=FLEADOC
          ITYPE=ITOICE
          ITYPEO=ITOCEAN
          IF (KOCEAN.ge.1) THEN
            QFIXR=.FALSE.
          ELSE
            QFIXR=.TRUE.
          END IF
        ELSE
          FLEAD=FLEADLK
          ITYPE=ITLKICE
          ITYPEO=ITLAKE
          QFIXR=.FALSE.
        END IF

        ACEFO=DMSI(1,I,J)
        ACEFI=DMSI(2,I,J)
        ENRGFO=DHSI(1,I,J)
        ENRGFI=DHSI(2,I,J)
        SALTO=DSSI(1,I,J)
        SALTI=DSSI(2,I,J)
#ifdef TRACERS_WATER
        TRO(:) = DTRSI(:,1,I,J)
        TRI(:) = DTRSI(:,2,I,J)
#endif
C**** ice formation diagnostics on the atmospheric grid
! define open ocean ice formation as frazil ice growth
        AIJ(I,J,IJ_SIGRFR)=AIJ(I,J,IJ_SIGRFR)+POCEAN*ACEFO
! define under ice formation as congelation ice growth
        AIJ(I,J,IJ_SIGRCG)=AIJ(I,J,IJ_SIGRCG)+POICE*ACEFI

        IF (FOCEAN(I,J).gt.0) THEN
          AIJ(I,J,IJ_FWIO)=AIJ(I,J,IJ_FWIO) - POCEAN*(ACEFO-SALTO)
     *         - POICE*(ACEFI-SALTI)
          AIJ(I,J,IJ_HTIO)=AIJ(I,J,IJ_HTIO) - POCEAN* ENRGFO
     *         - POICE* ENRGFI
          AIJ(I,J,IJ_STIO)=AIJ(I,J,IJ_STIO) - POCEAN* SALTO
     *         - POICE* SALTI
#ifdef TRACERS_WATER
          TAIJN(I,J,TIJ_ICOCFLX,:)=TAIJN(I,J,TIJ_ICOCFLX,:)
     *                    - POCEAN*TRO(:) - POICE*TRI(:)
#endif
        END IF
C**** open ocean diagnostics
        CALL INC_AJ(I,J,ITYPEO,J_SMELT,-SALTO *POCEAN)
        CALL INC_AJ(I,J,ITYPEO,J_HMELT,-ENRGFO*POCEAN)
        CALL INC_AJ(I,J,ITYPEO,J_IMELT,-ACEFO *POCEAN)
C**** Ice-covered ocean diagnostics
        CALL INC_AJ(I,J,ITYPE,J_SMELT,-SALTI *POICE)
        CALL INC_AJ(I,J,ITYPE,J_HMELT,-ENRGFI*POICE)
        CALL INC_AJ(I,J,ITYPE,J_IMELT,-ACEFI *POICE)
C**** regional diagnostics
        CALL INC_AREG(I,J,JR,J_IMELT,-(ACEFO *POCEAN+ACEFI *POICE))
        CALL INC_AREG(I,J,JR,J_HMELT,-(ENRGFO*POCEAN+ENRGFI*POICE))
        CALL INC_AREG(I,J,JR,J_SMELT,-(SALTO *POCEAN+SALTI *POICE))

        CALL ADDICE (SNOW,ROICE,HSIL,SSIL,MSI2,TSIL,ENRGFO,ACEFO,ACEFI,
     *       ENRGFI,SALTO,SALTI,
#ifdef TRACERS_WATER
     *       TRSIL,TRO,TRI,DTRIMP,
#endif
     *       DMIMP,DHIMP,DSIMP,FLEAD,QFIXR)

C**** RESAVE PROGNOSTIC QUANTITIES
        SNOWI(I,J) = SNOW
        MSI(I,J)=MSI2
        HSI(:,I,J) = HSIL(:)
        SSI(:,I,J) = SSIL(:)
#ifdef TRACERS_WATER
        TRSI(:,:,I,J) = TRSIL(:,:)
#endif
        IF (.not. QFIXR) THEN
          if (roice.gt.rsi(i,j)) ! ice from ocean
     *         call RESET_SURF_FLUXES(I,J,1,2,RSI(I,J),ROICE)
          if (roice.lt.rsi(i,j)) ! ocean from ice
     *         call RESET_SURF_FLUXES(I,J,2,1,1.-RSI(I,J),1.-ROICE)
          RSI(I,J)=ROICE
        ELSE
C**** save implicit mass-flux diagnostics
          AIJ(I,J,IJ_SMFX)=AIJ(I,J,IJ_SMFX)+ROICE*DMIMP
          CALL INC_AJ(I,J,ITOICE,J_IMPLM,-(DMIMP-DSIMP)*POICE)
          CALL INC_AJ(I,J,ITOICE,J_IMPLH,       -DHIMP *POICE)
        END IF
C**** set gtemp array
        GTEMP(1:2,2,I,J)=TSIL(1:2)
        GTEMPR(2,I,J)   =TSIL(1)+TF
#ifdef SCM
        if (I.eq.I_TARG.and.J.eq.J_TARG) then
            if (SCM_SURFACE_FLAG.ge.1) then
                GTEMP(1,2,I,J) = ATSKIN
                GTEMP(2,2,I,J) = ATSKIN
                GTEMPR(2,I,J) = ATSKIN + TF
            endif
        endif
#endif
#ifdef TRACERS_WATER
        GTRACER(:,2,I,J)=TRSIL(:,1)/(XSI(1)*(SNOW+ACE1I)-SSIL(1))
#endif
        FWSIM(I,J) = RSI(I,J)*(ACE1I+SNOW+MSI2-SUM(SSIL(1:LMI)))

C**** ACCUMULATE DIAGNOSTICS
        CALL INC_AJ(I,J,ITYPE,J_RSI ,      POICE)
        CALL INC_AJ(I,J,ITYPE,J_ACE1,ACE1I*POICE)
        CALL INC_AJ(I,J,ITYPE,J_ACE2, MSI2*POICE)
        CALL INC_AJ(I,J,ITYPE,J_SNOW, SNOW*POICE)
        IF (JR.ne.24) THEN
          CALL INC_AREG(I,J,JR,J_RSI ,      POICE)
          CALL INC_AREG(I,J,JR,J_SNOW, SNOW*POICE)
          CALL INC_AREG(I,J,JR,J_ACE1,ACE1I*POICE)
          CALL INC_AREG(I,J,JR,J_ACE2, MSI2*POICE)
        END IF
        AIJ(I,J,IJ_TSI)=AIJ(I,J,IJ_TSI)+
     *       POICE*(XSI(3)*TSIL(3)+XSI(4)*TSIL(4))
        AIJ(I,J,IJ_SSI1)=AIJ(I,J,IJ_SSI1)+POICE*(SSIL(1)+SSIL(2))/ACE1I
        AIJ(I,J,IJ_SSI2)=AIJ(I,J,IJ_SSI2)+POICE*(SSIL(3)+SSIL(4))
     *       /MSI(I,J)

#ifdef TRACERS_WATER
C**** Save sea ice tracer amount
      do n=1,ntm
        if (itime_tr0(n).le.itime .and. (tr_wd_TYPE(n).eq.nWater .or.
     *      tr_wd_TYPE(n).eq.nPART)) then
          taijn(i,j,tij_seaice,n)=taijn(i,j,tij_seaice,n)+
     *         POICE*sum(trsil(n,:))
        end if
      end do
#endif
      if (TSIL(1).lt.-100.) then
         write(6,*) "Seaice: T < -100. i,j,TSI = ",i,j,TSIL(1:LMI)
         call stop_model("Seaice too cold after ADDICE",255)
      end if

      END IF

      END DO
      END DO

C**** replicate ice values at the poles
      IF (HAVE_NORTH_POLE) THEN
        DO I=2,IM
          RSI(I,JM)=RSI(1,JM)
          MSI(I,JM)=MSI(1,JM)
          HSI(:,I,JM)=HSI(:,1,JM)
          SSI(:,I,JM)=SSI(:,1,JM)
          SNOWI(I,JM)=SNOWI(1,JM)
          GTEMP(1:2,2,I,JM)=GTEMP(1:2,2,1,JM)
          GTEMPR(2,I,JM)=GTEMPR(2,1,JM)
#ifdef TRACERS_WATER
          TRSI(:,:,I,JM) = TRSI(:,:,1,JM)
          GTRACER(:,2,I,JM) = GTRACER(:,2,1,JM)
#endif
          FWSIM(I,JM) = FWSIM(1,JM)
        END DO
      END IF
      IF (HAVE_SOUTH_POLE) THEN
        DO I=2,IM
          RSI(I,1)=RSI(1,1)
          MSI(I,1)=MSI(1,1)
          HSI(:,I,1)=HSI(:,1,1)
          SSI(:,I,1)=SSI(:,1,1)
          SNOWI(I,1)=SNOWI(1,1)
          GTEMP(1:2,2,I,1)=GTEMP(1:2,2,1,1)
          GTEMPR(2,I,1)=GTEMPR(2,1,1)
#ifdef TRACERS_WATER
          TRSI(:,:,I,1) = TRSI(:,:,1,1)
          GTRACER(:,2,I,1) = GTRACER(:,2,1,1)
#endif
          FWSIM(I,1) = FWSIM(1,1)
        END DO
      END IF
C****
      END SUBROUTINE FORM_SI

      SUBROUTINE vflx_OCEAN
!@sum  vflx_OCEAN saves quantities for OHT calculations
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm,focean
      USE DIAG_COM, only : oa
      USE SEAICE_COM, only : hsi,snowi
      USE FLUXES, only : fwsim
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_ATM, only : GET
      IMPLICIT NONE
      INTEGER I,J
      integer :: I_0, I_1, J_0, J_1
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C****
C****       DATA SAVED IN ORDER TO CALCULATE OCEAN TRANSPORTS
C****
C****       1  SNOWOI (INSTANTANEOUS AT NOON GMT)
C****       2  FWSIM  (INSTANTANEOUS AT NOON GMT)
C****       3  HSIT   (INSTANTANEOUS AT NOON GMT)
C****
      DO J=J_0, J_1
        DO I=I_0, I_1
          IF (FOCEAN(I,J).gt.0) THEN
            OA(I,J,1)=SNOWI(I,J)
            OA(I,J,2)=FWSIM(I,J)
            OA(I,J,3)=SUM(HSI(:,I,J))
          END IF
        END DO
      END DO

      RETURN
C****
      END SUBROUTINE vflx_OCEAN

      SUBROUTINE init_ice(iniOCEAN,istart)
!@sum  init_ice initialises ice arrays
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : rhows,tf
      USE MODEL_COM, only : im,jm,kocean,focean,flake0
#ifdef SCM
      USE MODEL_COM, only : I_TARG,J_TARG
      USE SCMCOM, only : iu_scm_prt,SCM_SURFACE_FLAG,ATSKIN
#endif
      USE SEAICE, only : xsi,ace1i,ac2oim,ssi0,tfrez,oi_ustar0,silmfac
     *     ,lmi,snow_ice,Ti,Ei,seaice_thermo
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi,pond_melt,flag_dsws
#ifdef TRACERS_WATER
     *     ,trsi,ntm
#endif
      USE FLUXES, only : gtemp,ui2rho,fwsim,msicnv,gtempr
#ifdef TRACERS_WATER
     *     ,gtracer
#endif
      USE DIAG_COM, only : npts,icon_OMSI,icon_OHSI,icon_OSSI,icon_LMSI
     *     ,icon_LHSI,conpt0
      USE PARAM
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_ATM, only : GET
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE. , iniOCEAN
      CHARACTER CONPT(NPTS)*10
      INTEGER I,J,istart
      REAL*8 MSI1,TFO
      integer :: I_0, I_1, J_0, J_1
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** set up a default ice-ocean stress field. This can be changed by
C**** adjusting oi_ustar0 in the parameter list. If ice dynamics
C**** is used, this is overwritten.
      call sync_param("oi_ustar0",oi_ustar0)
      UI2rho = rhows*(oi_ustar0)**2

C**** Adjust degree of lateral melt by changing silmfac
C**** Default is 2.5d-8, but could be changed by a factor of 2.
      call sync_param("silmfac",silmfac)

C**** Decide whether snow_ice formation is allowed
      call sync_param("snow_ice",snow_ice)

C**** Define the ice thermodynamics (SI or BP)
      call sync_param("seaice_thermo",seaice_thermo)

C**** clean up ice fraction/sea ice salinity possibly incorrect in I.C.
      if (istart.lt.9) then
      DO J=J_0, J_1
      DO I=I_0, I_1
        IF (FOCEAN(I,J)+FLAKE0(I,J).eq.0 .and. RSI(i,j).gt.0) RSI(I,J)=0
        IF (RSI(I,J).gt.0 .and. FLAKE0(I,J).gt.0) SSI(:,I,J)=0.
      END DO
      END DO
      end if

      IF (KOCEAN.EQ.0.and.iniOCEAN) THEN
C****   set defaults for no ice case
        DO J=J_0, J_1
        DO I=I_0, I_1
          IF (RSI(I,J).le.0) THEN
            MSI1        =ACE1I
            MSI(I,J)    =AC2OIM
            SNOWI(I,J)  =0.
            IF (FOCEAN(I,J).gt.0) THEN
              SSI(1:2,I,J)=SSI0*XSI(1:2)*ACE1I
              SSI(3:4,I,J)=SSI0*XSI(3:4)*AC2OIM
              TFO = -1.87d0  ! reasonable value, doesn't really matter
              HSI(1:2,I,J)=Ei(TFO,1d3*SSI0)*XSI(1:2)*ACE1I
              HSI(3:4,I,J)=Ei(TFO,1d3*SSI0)*XSI(3:4)*AC2OIM
            ELSE
              SSI(:,I,J)  = 0.
              TFO = 0.
              HSI(1:2,I,J)=Ei(TFO,0d0)*XSI(1:2)*ACE1I
              HSI(3:4,I,J)=Ei(TFO,0d0)*XSI(3:4)*AC2OIM
            END IF
#ifdef TRACERS_WATER
            TRSI(:,:,I,J)=0.
#endif
            pond_melt(i,j) = 0.
            flag_dsws(i,j) = .FALSE.
          END IF

        END DO
        END DO
      END IF
C**** set GTEMP etc. array for ice
      DO J=J_0, J_1
      DO I=I_0, I_1
        MSI1=SNOWI(I,J)+ACE1I
        GTEMP(1,2,I,J)=Ti(HSI(1,I,J)/(XSI(1)*MSI1),1d3*SSI(1,I,J
     *       )/(XSI(1)*MSI1))
        GTEMP(2,2,I,J)=Ti(HSI(2,I,J)/(XSI(2)*MSI1),1d3*SSI(2,I,J
     *       )/(XSI(2)*MSI1))
        GTEMPR(2,I,J) = GTEMP(1,2,I,J)+TF
#ifdef SCM
        if (I.eq.I_TARG.and.J.eq.J_TARG) then
            if (SCM_SURFACE_FLAG.ge.1) then
                GTEMP(1,2,I,J) = ATSKIN
                GTEMP(2,2,I,J) = ATSKIN
                GTEMPR(2,I,J) = ATSKIN + TF
            endif
        endif
#endif
#ifdef TRACERS_WATER
        GTRACER(:,2,I,J) = TRSI(:,1,I,J)/(XSI(1)*MSI1-SSI(1,I,J))
#endif
        FWSIM(I,J) = RSI(I,J)*(MSI1+MSI(I,J)-SUM(SSI(1:LMI,I,J)))
        MSICNV(I,J)=0.   ! always initialise to zero
      END DO
      END DO

C**** Set conservation diagnostics for ice mass, energy, salt
      CONPT=CONPT0
      CONPT(3)="LAT. MELT" ; CONPT(4)="PRECIP"
      CONPT(5)="THERMO"    ; CONPT(6)="ADVECT"
      CONPT(8)="OCN FORM"
      QCON=(/ F, F, T, T, T, T, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT,"OICE MAS","(KG/M^2)        ",
     *     "(10**-9 KG/SM^2)",1d0,1d9,icon_OMSI)
      QCON=(/ F, F, T, T, T, T, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT,"OICE ENR","(10**6 J/M^2)   ",
     *     "(10**-3 W/M^2)  ",1d-6,1d3,icon_OHSI)
      QCON=(/ F, F, T, T, T, T, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT,"OICE SLT","(10**-3 KG/M^2) ",
     *     "(10**-12KG/SM^2)",1d3,1d12,icon_OSSI)
      CONPT(8)="LK FORM"
      QCON=(/ F, F, T, T, T, F, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT,"LKICE MS","(KG/M^2)        ",
     *     "(10**-9 KG/SM^2)",1d0,1d9,icon_LMSI)
      QCON=(/ F, F, T, T, T, F, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT,"LKICE EN","(10**6 J/M^2)   ",
     *     "(10**-3 W/M^2)  ",1d-6,1d3,icon_LHSI)
C****
      END SUBROUTINE init_ice

      SUBROUTINE conserv_OMSI(ICE)
!@sum  conserv_MSI calculates total amount of snow and ice over ocean
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,focean
      USE GEOM, only : imaxj
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : rsi,msi,snowi
      USE DOMAIN_DECOMP_ATM, only : GRID, GET
      IMPLICIT NONE
!@var ICE total ocean snow and ice mass (kg/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: ICE
      INTEGER I,J

c**** Extract useful domain information from grid
      INTEGER J_0, J_1, J_0H, J_1H ,I_0,I_1
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(GRID, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H ,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE    ,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE    )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        ICE(I,J)=RSI(I,J)*(MSI(I,J)+ACE1I+SNOWI(I,J))*FOCEAN(I,J)
      END DO
      END DO
      IF (HAVE_SOUTH_POLE) ICE(2:im,1) =ICE(1,1)
      IF (HAVE_NORTH_POLE) ICE(2:im,JM)=ICE(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_OMSI

      SUBROUTINE conserv_OHSI(EICE)
!@sum  conserv_HSI calculates total ice energy over ocean
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,focean
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : rsi,hsi
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      IMPLICIT NONE
!@var EICE total ocean snow and ice energy (J/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: EICE
      INTEGER I,J

c**** Extract useful domain information from grid
      INTEGER J_0, J_1, J_0H, J_1H ,I_0,I_1
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(GRID, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H ,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE    ,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE    )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        EICE(I,J)=RSI(I,J)*FOCEAN(I,J)*SUM(HSI(:,I,J))
      END DO
      END DO
      IF (HAVE_SOUTH_POLE) EICE(2:im,1) =EICE(1,1)
      IF (HAVE_NORTH_POLE) EICE(2:im,JM)=EICE(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_OHSI

      SUBROUTINE conserv_OSSI(SALT)
!@sum  conserv_SSI calculates total amount of salt in ocean ice
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,focean
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : rsi,ssi,lmi
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      IMPLICIT NONE
!@var SALT total salt in ocean ice (kg/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: SALT
      INTEGER I,J

c**** Extract useful domain information from grid
      INTEGER J_0, J_1, J_0H, J_1H ,I_0,I_1
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(GRID, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H ,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE    ,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE    )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        IF (FOCEAN(I,J).gt.0) THEN
          SALT(I,J)=FOCEAN(I,J)*RSI(I,J)*SUM(SSI(:,I,J))
        ELSE
          SALT(I,J)=0
        END IF
      END DO
      END DO
      IF (HAVE_SOUTH_POLE) SALT(2:im,1) =SALT(1,1)
      IF (HAVE_NORTH_POLE) SALT(2:im,JM)=SALT(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_OSSI

      SUBROUTINE conserv_LMSI(ICE)
!@sum  conserv_LMSI calculates total amount of snow and ice over lakes
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE GEOM, only : imaxj
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : rsi,msi,snowi
      USE LAKES_COM, only : flake
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      IMPLICIT NONE
!@var ICE total lake snow and ice mass (kg/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: ICE
      INTEGER I,J

c**** Extract useful domain information from grid
      INTEGER J_0, J_1, J_0H, J_1H ,I_0,I_1
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(GRID, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H ,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE    ,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE    )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        ICE(I,J)=RSI(I,J)*(MSI(I,J)+ACE1I+SNOWI(I,J))*FLAKE(I,J)
      END DO
      END DO
      IF (HAVE_SOUTH_POLE) ICE(2:im,1) =ICE(1,1)
      IF (HAVE_NORTH_POLE) ICE(2:im,JM)=ICE(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_LMSI

      SUBROUTINE conserv_LHSI(EICE)
!@sum  conserv_LHSI calculates total ice energy over lakes
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : rsi,hsi
      USE LAKES_COM, only : flake
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      IMPLICIT NONE
!@var EICE total lake snow and ice energy (J/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: EICE
      INTEGER I,J

c**** Extract useful domain information from grid
      INTEGER J_0, J_1, J_0H, J_1H ,I_0,I_1
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(GRID, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H ,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE    ,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE    )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        EICE(I,J)=RSI(I,J)*FLAKE(I,J)*SUM(HSI(:,I,J))
      END DO
      END DO
      IF (HAVE_SOUTH_POLE) EICE(2:im,1) =EICE(1,1)
      IF (HAVE_NORTH_POLE) EICE(2:im,JM)=EICE(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_LHSI

      SUBROUTINE daily_ice
!@sum daily_ice performs ice processes that are needed everyday
!@auth Gavin Schmidt
      USE CONSTANT, only : tf
      USE MODEL_COM, only : jm
#ifdef SCM
      USE MODEL_COM, only : I_TARG,J_TARG
      USE SCMCOM, only : iu_scm_prt, SCM_SURFACE_FLAG,ATSKIN
#endif
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : msi,hsi,ssi,rsi,snowi
#ifdef TRACERS_WATER
     *     ,trsi
#endif
      USE SEAICE, only : ace1i,xsi,lmi,Ti
      USE FLUXES, only : gtemp,fwsim,gtempr
#ifdef TRACERS_WATER
     &     ,gtracer
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_ATM, only : GET
      IMPLICIT NONE
      INTEGER I,J, J_0, J_1 ,I_0,I_1
      REAL*8 MSI1
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0, J_1
      DO I=I_0,IMAXJ(J)
C**** set GTEMP etc. array for ice (to deal with daily_lake changes)
        MSI1=SNOWI(I,J)+ACE1I
        GTEMP(1,2,I,J)=Ti(HSI(1,I,J)/(XSI(1)*MSI1),1d3*SSI(1,I,J
     *       )/(XSI(1)*MSI1))
        GTEMP(2,2,I,J)=Ti(HSI(2,I,J)/(XSI(2)*MSI1),1d3*SSI(2,I,J
     *       )/(XSI(2)*MSI1))
        GTEMPR(2,I,J) = GTEMP(1,2,I,J)+TF
#ifdef SCM
        if (I.eq.I_TARG.and.J.eq.J_TARG) then
            if (SCM_SURFACE_FLAG.ge.1) then
                GTEMP(1,2,I,J) = ATSKIN
                GTEMP(2,2,I,J) = ATSKIN
                GTEMPR(2,I,J) = ATSKIN + TF
            endif
        endif
#endif
#ifdef TRACERS_WATER
        GTRACER(:,2,I,J) = TRSI(:,1,I,J)/(XSI(1)*MSI1-SSI(1,I,J))
#endif
        FWSIM(I,J) = RSI(I,J)*(MSI1+MSI(I,J)-SUM(SSI(1:LMI,I,J)))
      END DO
      END DO

      RETURN
      END SUBROUTINE daily_ice

