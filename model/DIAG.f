#include "rundeck_opts.h"
#ifdef SKIP_TRACER_DIAGS
#undef TRACERS_SPECIAL_O18
#endif

#ifdef CUBED_SPHERE
#define SLP_FROM_T1
#endif

!@sum  DIAG ModelE diagnostic calculations
!@auth G. Schmidt/J. Lerner/R. Ruedy/M. Kelley
!@ver  1.0
C**** AJ(J,N)  (ZONAL SUM OVER LONGITUDE AND TIME)
C****   See j_defs for contents
C****                                                             IDACC
C****
C**** CONTENTS OF AJL(J,L,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   See jl_defs for contents
C****
C**** CONTENTS OF ASJL(J,L,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   See jls_defs for contents
C****
C**** CONTENTS OF AIJ(I,J,N)  (SUM OVER TIME OF)
C****   See ij_defs for contents
C****
C**** CONTENTS OF AIL(I,L,N)  (SUM OVER TIME OF)
C****   See il_defs for contents
C****
C**** CONTENTS OF IDACC(N), NUMBER OF ACCUMULATION TIMES OF
C****   1  SOURCE TERMS  (dt: DTSRC)
C****   2  RADIATION SOURCE TERMS  (dt: NRAD*DTsrc)
C****   3  SURFACE INTERACTION SOURCE TERMS  (dt: NDASF*DTsrc+DTsurf)
C****   4  QUANTITIES IN DIAGA  (dt: NDAA*DTsrc+2*DTdyn)
C****   5  ENERGY NUMBERS IN DIAG4  (DETERMINED BY NDA4)
C****   6  KINETIC ENERGY IN DIAG5 FROM DYN'CS (dt: NDA5K*DTsrc+2*DTdyn)
C****   7  ENERGY IN DIAG5 FROM DYNAMICS  (dt: NDA5D*DTsrc)
C****   8  ENERGY IN DIAG5 FROM SOURCES  (DETERMINED BY NDA5S)
C****   9  WAVE ENERGY IN DIAG7  (dt: 12 HOURS)
C****  10  ENERGY IN DIAG5 FROM FILTER  (DT: NFILTR*DTsrc)
C****  11  NOT USED
C****  12  ALWAYS =1 (UNLESS SEVERAL RESTART FILES WERE ACCUMULATED)
C****

      MODULE DIAG_LOC
!@sum DIAG_LOC is a local module for some saved diagnostic calculations
!@auth Gavin Schmidt
      USE MODEL_COM, only : im,jm,lm
      IMPLICIT NONE
      SAVE
C**** Variables passed from DIAGA to DIAGB
!@var W,TX vertical velocity and in-situ temperature calculations
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: W
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TX

C**** Some local constants
!@var JET, LDEX model levels for various pressures
!@var LUPA,LDNA shorthand for above/below levels
!@var PMO,PLO,PM,PL some shorthand pressure level
      INTEGER :: JET
      INTEGER, DIMENSION(3) :: LDEX
      REAL*8, DIMENSION(LM) :: LUPA,LDNA
      REAL*8, DIMENSION(LM) :: PMO,PLO
      REAL*8, DIMENSION(LM+1) :: PM,PL

      END MODULE DIAG_LOC

      SUBROUTINE ALLOC_DIAG_LOC(grid)
      USE DOMAIN_DECOMP_ATM, only : GET
      USE DOMAIN_DECOMP_ATM, only : DIST_GRID
      USE MODEL_COM, only : lm
      USE DIAG_LOC, only  : W,TX
      IMPLICIT NONE
      LOGICAL, SAVE :: init=.false.
      INTEGER :: I_0H,I_1H,J_0H,J_1H
      INTEGER :: IER
      TYPE(DIST_GRID) :: grid

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = GRID%I_STRT_HALO
      I_1H = GRID%I_STOP_HALO

      ALLOCATE( W(I_0H:I_1H, J_0H:J_1H, LM),
     &         TX(I_0H:I_1H, J_0H:J_1H, LM),
     &     STAT = IER)

      !hack hack hack!
      TX(:,:,:) = 0.d0

      RETURN
      END SUBROUTINE ALLOC_DIAG_LOC

      SUBROUTINE DIAGA
!@sum  DIAGA accumulate various diagnostics during dynamics
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,kapa,lhe,lhs,sha,bygrav,tf
     *     ,rvap,gamd,teeny,undef,radius,omega,kg2mb,mair   
      USE MODEL_COM, only : im,jm,lm,ls1,idacc,ptop
     *     ,pmtop,psfmpt,mdyn,mdiag,sig,sige,dsig,zatmo,WM,ntype,ftype
     *     ,u,v,t,p,q,lm_req,req_fac_m,pmidl00
      USE GEOM, only : sinlat2d,coslat2d,axyp,imaxj,ddy_ci,ddy_cj,
     &     lon2d_dg,byaxyp
      USE RAD_COM, only : rqt
      USE DIAG_COM, only : ia_dga,jreg,
     *     aijl=>aijl_loc
     *     ,aij=>aij_loc,ij_dtdp,ij_phi1k,ij_pres,ij_slpq,ij_presq
     *     ,ij_slp,ij_t850,ij_t500,ij_t300,ij_t100,ij_q850,ij_q500
     *     ,ij_rh700,ij_t700,ij_q700,ij_q100,ij_rh100
     *     ,ij_RH1,ij_RH850,ij_RH500,ij_RH300,ij_qm,ij_q300,ij_ujet
     *     ,ij_vjet,j_tx1,j_tx,j_qp,j_dtdjt,j_dtdjs,j_dtdgtr,j_dtsgst
     &     ,ijl_dp,ijk_dp,ijl_u,ijl_v,ijl_w,ijk_tx,ijk_q,ijk_rh
     *     ,j_rictr,j_rostr,j_ltro,j_ricst,j_rosst,j_lstr,j_gamm,j_gam
     *     ,j_gamc,lstr,kgz_max,pmb,ght,ple
     *     ,jl_dtdyn,jl_dpa
     *     ,jl_epacwt,jl_uepac,jl_vepac,jl_wepac
     *     ,jl_wpacwt,jl_uwpac,jl_vwpac,jl_wwpac
     *     ,jk_dpwt,jk_tx,jk_hght,jk_q,jk_rh,jk_cldh2o
     *     ,jk_cldwtr,jk_cldice
     *     ,ij_p850,z_inst,rh_inst,t_inst,plm,ij_p1000,ij_p925,ij_p700
     *     ,ij_p600,ij_p500,ijl_templ,ijl_gridh,ijl_husl,ijl_zL
#ifdef TRACERS_SPECIAL_Shindell
     *     ,o_inst,x_inst,n_inst,m_inst
#endif
#ifdef TES_LIKE_DIAGS
     *     ,t_more,q_more,kgz_max_more,PMBmore
#ifdef TRACERS_SPECIAL_Shindell
     *     ,o_more,x_more,n_more,m_more
#endif
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE TRACER_COM, only: trm,mass2vol,n_CO,n_Ox,n_NOx
      USE TRCHEM_Shindell_COM, only : mNO2
#endif
      USE DYNAMICS, only : pk,pek,phi,pmid,pdsig,plij, SD,pedn,am
     &     ,ua=>ualij,va=>valij,wcp
      USE PBLCOM, only : tsavg
      USE CLOUDS_COM, only : svlhx
      USE DIAG_LOC, only : w,tx,jet
      USE DOMAIN_DECOMP_ATM, only : GET, GRID, HALO_UPDATE
      USE GETTIME_MOD
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     TX_TROP,TX_STRAT
      REAL*8, DIMENSION(LM_REQ) :: TRI

      REAL*8, PARAMETER :: ONE=1.,P1000=1000.
      INTEGER :: I,IM1,J,K,L,JR,
     &     IP1,LR,IT
      REAL*8 THBAR ! external
      REAL*8 ::
     &     BBYGV,BYSDSG,DLNP01,DLNP12,DLNP23,DBYSD,
     &     DXYPJ,
     *     ESEPS,GAMC,GAMM,GAMX,
     &     PDN,PE,PHI_REQ,PIJ,pfact,chemL,chemLm1,
     *     PKE,PL,PRT,W2MAX,RICHN,
     *     ROSSN,ROSSL,BYFCOR,BYBETA,BYBETAFAC,NH,SS,THETA,
     *     TZL,X,TIJK,QIJK,DTXDY
      LOGICAL qpress,qabove
      INTEGER nT,nQ,nRH
      REAL*8, PARAMETER :: EPSLON=1.

      REAL*8 QSAT, SLP, PS, ZS, TS_SLP, QLH, begin

      real*8, dimension(lm+1) :: pecp,pedge
      real*8, dimension(lm) :: dpwt,txdp,phidp,qdp,rhdp,wmdp,rh
      real*8, dimension(lm) :: wmliqdp,wmfrzdp
      integer, parameter :: lmxmax=2*lm
      integer :: lmx
      real*8, dimension(lmxmax) :: dpx
      integer, dimension(lmxmax) :: lmod,lcp

      INTEGER :: J_0, J_1, J_0S, J_1S, I_0,I_1,I_0H,I_1H
     &     ,IM1S,IP1S,IP1E
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GETTIME(BEGIN)

      CALL GET(grid, J_STRT=J_0,         J_STOP=J_1,
     &               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP
      I_0H = GRID%I_STRT_HALO
      I_1H = GRID%I_STOP_HALO

c
c get winds on the atmospheric primary grid
c
      call recalc_agrid_uv

      IDACC(ia_dga)=IDACC(ia_dga)+1

      BYSDSG=1./(1.-SIGE(LM+1))
      DLNP01=LOG(pmidl00(lm)/PLM(LM+1))
      DLNP12=LOG(REQ_FAC_M(1)/REQ_FAC_M(2))  ! LOG(.75/.35)
      DLNP23=LOG(REQ_FAC_M(2)/REQ_FAC_M(3))  ! LOG(.35/.1)
C****
C**** FILL IN HUMIDITY AND SIGMA DOT ARRAYS AT THE POLES
C****
      IF(HAVE_SOUTH_POLE) THEN
        DO L=1,LM
          DO I=2,IM
            Q(I,1,L)=Q(1,1,L)
          END DO
        END DO
      ENDIF        ! HAVE_SOUTH_POLE
      IF(HAVE_NORTH_POLE) THEN
        DO L=1,LM
          DO I=2,IM
            Q(I,JM,L)=Q(1,JM,L)
          END DO
        END DO
      ENDIF        ! HAVE_NORTH_POLE
C****
C**** CALCULATE PK AND TX, THE REAL TEMPERATURE
C****
      IF(HAVE_SOUTH_POLE) THEN
        DO L=1,LM
          TX(1,1,L)=T(1,1,L)*PK(L,1,1)
          DO I=2,IM
            T(I,1,L)=T(1,1,L)
            TX(I,1,L)=TX(1,1,L)
          END DO
        END DO
      ENDIF        ! HAVE_SOUTH_POLE
      IF(HAVE_NORTH_POLE) THEN
        DO L=1,LM
          TX(1,JM,L)=T(1,JM,L)*PK(L,1,JM)
          DO I=2,IM
            T(I,JM,L)=T(1,JM,L)
            TX(I,JM,L)=TX(1,JM,L)
          END DO
        END DO
      ENDIF          ! HAVE_NORTH_POLE

      DO L=1,LM
      DO J=J_0S,J_1S
        DO I=I_0,I_1
          TX(I,J,L)=T(I,J,L)*PK(L,I,J)
        END DO
      END DO
      END DO


C****
C**** J LOOPS FOR ALL PRIMARY GRID ROWS
C****
      DO J=J_0,J_1
C**** NUMBERS ACCUMULATED FOR A SINGLE LEVEL
        DO I=I_0,IMAXJ(J)
          DXYPJ=AXYP(I,J)
          JR=JREG(I,J)
          DO IT=1,NTYPE
            CALL INC_AJ(I,J,IT,J_TX1,(TX(I,J,1)-TF)*FTYPE(IT,I,J))
          END DO
          CALL INC_AREG(I,J,JR,J_TX1,(TX(I,J,1)-TF))
          PS=P(I,J)+PTOP
          ZS=BYGRAV*ZATMO(I,J)
          AIJ(I,J,IJ_PRES)=AIJ(I,J,IJ_PRES)+ PS
#ifdef SLP_FROM_T1
          TS_SLP=T(I,J,1)*PEK(1,I,J) ! todo: check if tmom(mz) helps
#else
          TS_SLP=TSAVG(I,J)
#endif
          AIJ(I,J,IJ_SLP)=AIJ(I,J,IJ_SLP)+SLP(PS,TS_SLP,ZS)-P1000
C**** calculate pressure diags including water
          PS=PS+SUM((Q(I,J,:)+WM(I,J,:))*AM(:,I,J))*kg2mb
          AIJ(I,J,IJ_PRESQ)=AIJ(I,J,IJ_PRESQ)+ PS
          AIJ(I,J,IJ_SLPQ)=AIJ(I,J,IJ_SLPQ)+SLP(PS,TS_SLP,ZS)-P1000

          AIJ(I,J,IJ_RH1)=AIJ(I,J,IJ_RH1)+Q(I,J,1)/QSAT(TX(I,J,1),LHE,
     *        PMID(1,I,J))
        END DO
c        APJ(J,1)=APJ(J,1)+PI(J)
#ifdef TES_LIKE_DIAGS
C**** Calculate T and Q at specific millibar levels for TES diags
C**** Follows logic for geopotential section following this...
        do I=I_0,IMAXJ(J)
          K=1
          L=1
          q_more(:,i,j) = undef ; t_more(:,i,j) = undef
#ifdef TRACERS_SPECIAL_Shindell
          o_more(:,i,j) = undef ; x_more(:,i,j) = undef
          n_more(:,i,j) = undef ; m_more(:,i,j) = undef
#endif
 1720     L=L+1
          pdn=pmid(L-1,I,J)
          pl=pmid(L,I,J)
          if (PMBmore(K)<pl .AND. L<LM) goto 1720
 1740     continue
          qabove = pmbmore(k).le.pedn(l-1,i,j)
          pfact=(PMBmore(K)-PL)/(PDN-PL)
          if(ABS(TX(I,J,L)-TX(I,J,L-1))>=epslon)then
            if(qabove) TIJK=(TX(I,J,L)-TF
     *      +(TX(I,J,L-1)-TX(I,J,L))*LOG(PMBmore(K)/PL)/LOG(PDN/PL))
          else
            if(qabove) TIJK=TX(I,J,L)-TF
          endif
          if (qabove) then
            QIJK=Q(I,J,L)+(Q(I,J,L-1)-Q(I,J,L))*pfact
            q_more(K,I,J)=QIJK
            t_more(K,I,J)=TIJK
#ifdef TRACERS_SPECIAL_Shindell
              chemL=1.d6*trm(i,j,L,n_Ox)*mass2vol(n_Ox)/
     &        (am(L,i,j)*axyp(i,j))
              chemLm1=1.d6*trm(i,j,L-1,n_Ox)*mass2vol(n_Ox)/
     &        (am(L-1,i,j)*axyp(i,j))
            o_more(K,I,J)= chemL+(chemLm1-chemL)*pfact
              chemL=1.d6*trm(i,j,L,n_NOx)*mass2vol(n_NOx)/
     &        (am(L,i,j)*axyp(i,j))
              chemLm1=1.d6*trm(i,j,L-1,n_NOx)*mass2vol(n_NOx)/
     &        (am(L-1,i,j)*axyp(i,j))
            x_more(K,I,J)= chemL+(chemLm1-chemL)*pfact
              chemL=1.d6*mNO2(i,j,L)
              chemLm1=1.d6*mNO2(i,j,L-1)
            n_more(K,I,J)= chemL+(chemLm1-chemL)*pfact
              chemL=1.d6*trm(i,j,L,n_CO)*mass2vol(n_CO)/
     &        (am(L,i,j)*axyp(i,j))
              chemLm1=1.d6*trm(i,j,L-1,n_CO)*mass2vol(n_CO)/
     &        (am(L-1,i,j)*axyp(i,j))
            m_more(K,I,J)= chemL+(chemLm1-chemL)*pfact
#endif
          endif
          if(K < KGZ_max_more) then
            K=K+1
            if(PMBmore(K)<pl .and. L<LM) goto 1720
            goto 1740
          endif
        enddo ! I
#endif

C**** CALCULATE GEOPOTENTIAL HEIGHTS AT SPECIFIC MILLIBAR LEVELS
        DO I=I_0,IMAXJ(J)
          K=1
          L=1
          rh_inst(:,i,j) = undef ; t_inst(:,i,j) = undef
          z_inst(:,i,j) = undef
#ifdef TRACERS_SPECIAL_Shindell
          o_inst(:,i,j) = undef ; x_inst(:,i,j) = undef
          n_inst(:,i,j) = undef ; m_inst(:,i,j) = undef
#endif
 172      L=L+1
          PDN=PMID(L-1,I,J)
          PL=PMID(L,I,J)
          IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 172
C**** Select pressure levels on which to save temperature and humidity
C**** Use masking for 850 mb temp/humidity
 174      qpress = .false.
          qabove = pmb(k).le.pedn(l-1,i,j)
          SELECT CASE (NINT(PMB(K)))
          CASE (850)            ! 850 mb
            nT = IJ_T850 ; nQ = IJ_Q850 ; nRH = IJ_RH850 ; qpress=.true.
            if (.not. qabove) qpress = .false.
            if (qpress) aij(i,j,ij_p850) = aij(i,j,ij_p850) + 1.
          CASE (700)            ! 700 mb
            nT = IJ_T700 ; nQ = IJ_Q700 ; nRH = IJ_RH700 ; qpress=.true.
          CASE (500)            ! 500 mb
            nT = IJ_T500 ; nQ = IJ_Q500 ; nRH = IJ_RH500 ; qpress=.true.
          CASE (300)            ! 300 mb
            nT = IJ_T300 ; nQ = IJ_Q300 ; nRH = IJ_RH300 ; qpress=.true.
          CASE (100)            ! 100 mb
            nT = IJ_T100 ; nQ = IJ_Q100 ; nRH = IJ_RH100 ; qpress=.true.
          END SELECT
C**** calculate geopotential heights + temperatures
          IF (ABS(TX(I,J,L)-TX(I,J,L-1)).GE.EPSLON) THEN
            BBYGV=(TX(I,J,L-1)-TX(I,J,L))/(PHI(I,J,L)-PHI(I,J,L-1))
            AIJ(I,J,IJ_PHI1K-1+K)=AIJ(I,J,IJ_PHI1K-1+K)+(PHI(I,J,L)
     *           -TX(I,J,L)*((PMB(K)/PL)**(RGAS*BBYGV)-1.)/BBYGV-GHT(K)
     *           *GRAV)
            IF (qabove) then
              TIJK=(TX(I,J,L)-TF
     *           +(TX(I,J,L-1)-TX(I,J,L))*LOG(PMB(K)/PL)/LOG(PDN/PL))
              Z_inst(K,I,J)=(PHI(I,J,L)
     *           -TX(I,J,L)*((PMB(K)/PL)**(RGAS*BBYGV)-1.)/BBYGV-GHT(K)
     *             *GRAV)
            END IF
          ELSE
            AIJ(I,J,IJ_PHI1K-1+K)=AIJ(I,J,IJ_PHI1K-1+K)+(PHI(I,J,L)
     *           -RGAS*TX(I,J,L)*LOG(PMB(K)/PL)-GHT(K)*GRAV)
            IF (qabove) then
              TIJK=TX(I,J,L)-TF
              Z_inst(K,I,J)=(PHI(I,J,L)
     *             -RGAS*TX(I,J,L)*LOG(PMB(K)/PL)-GHT(K)*GRAV)
            END IF
          END IF
          if (qabove) then
            QIJK=Q(I,J,L)+(Q(I,J,L-1)-Q(I,J,L))*(PMB(K)-PL)/(PDN-PL)
            RH_inst(K,I,J)=QIJK/qsat(TIJK+TF,LHE,PMB(K))
            T_inst(K,I,J) =TIJK
#ifdef TRACERS_SPECIAL_Shindell
            pfact=(PMB(K)-PL)/(PDN-PL)
              chemL=1.d6*trm(i,j,L,n_Ox)*mass2vol(n_Ox)/
     &        (am(L,i,j)*axyp(i,j))
              chemLm1=1.d6*trm(i,j,L-1,n_Ox)*mass2vol(n_Ox)/
     &        (am(L-1,i,j)*axyp(i,j))
            o_inst(K,I,J)= chemL+(chemLm1-chemL)*pfact
              chemL=1.d6*trm(i,j,L,n_NOx)*mass2vol(n_NOx)/
     &        (am(L,i,j)*axyp(i,j))
              chemLm1=1.d6*trm(i,j,L-1,n_NOx)*mass2vol(n_NOx)/
     &        (am(L-1,i,j)*axyp(i,j))
            x_inst(K,I,J)= chemL+(chemLm1-chemL)*pfact
              chemL=1.d6*mNO2(i,j,L)
              chemLm1=1.d6*mNO2(i,j,L-1)
            n_inst(K,I,J)= chemL+(chemLm1-chemL)*pfact
              chemL=1.d6*trm(i,j,L,n_CO)*mass2vol(n_CO)/
     &        (am(L,i,j)*axyp(i,j))
              chemLm1=1.d6*trm(i,j,L-1,n_CO)*mass2vol(n_CO)/
     &        (am(L-1,i,j)*axyp(i,j))
            m_inst(K,I,J)= chemL+(chemLm1-chemL)*pfact
#endif
            if (qpress) then
              AIJ(I,J,nT)=AIJ(I,J,nT)+TIJK
              AIJ(I,J,nQ)=AIJ(I,J,nQ)+QIJK
              if (PMB(K).ge.500) then  ! w.r.t. water
                AIJ(I,J,nRH)=AIJ(I,J,nRH)+QIJK/qsat(TIJK+TF,LHe,PMB(K))
              else                     ! w.r.t ice above 500mb
                AIJ(I,J,nRH)=AIJ(I,J,nRH)+QIJK/qsat(TIJK+TF,LHs,PMB(K))
              end if
            end if
          end if
C****
          IF (K.LT.KGZ_max) THEN
            K=K+1
            IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 172
            GO TO 174
          END IF
C**** BEGIN AMIP
          IF((P(I,J)+PTOP).LT.1000.)AIJ(I,J,IJ_P1000)=
     *      AIJ(I,J,IJ_P1000)+1.
          IF((P(I,J)+PTOP).LT.925.)AIJ(I,J,IJ_P925)=AIJ(I,J,IJ_P925)+1.
          IF((P(I,J)+PTOP).LT.700.)AIJ(I,J,IJ_P700)=AIJ(I,J,IJ_P700)+1.
          IF((P(I,J)+PTOP).LT.600.)AIJ(I,J,IJ_P600)=AIJ(I,J,IJ_P600)+1.
          IF((P(I,J)+PTOP).LT.500.)AIJ(I,J,IJ_P500)=AIJ(I,J,IJ_P500)+1.
C**** END AMIP
        END DO
      END DO

C**** ACCUMULATION OF TEMP., POTENTIAL TEMP., Q, AND RH
      DO J=J_0,J_1
        DO L=1,LM
          DBYSD=DSIG(L)*BYSDSG
          DO I=I_0,IMAXJ(J)
            DXYPJ=AXYP(I,J)
            JR=JREG(I,J)
            PIJ=PLIJ(L,I,J)
            aijl(i,j,l,ijl_dp) = aijl(i,j,l,ijl_dp) + pdsig(l,i,j)
            call inc_ajl(i,j,l,jl_dpa,pdsig(l,i,j))
c ajl(jl_dtdyn) was incremented by -t(i,j,l) before dynamics
            call inc_ajl(i,j,l,jl_dtdyn,tx(i,j,l)*pdsig(l,i,j))
            AIJ(I,J,IJ_QM)=AIJ(I,J,IJ_QM)+Q(I,J,L)*AM(L,I,J)
            aijl(i,j,L,ijl_tempL)=aijl(i,j,L,ijl_tempL)+TX(i,j,L)
            aijl(i,j,L,ijl_husL)=aijl(i,j,L,ijl_husL)+Q(i,j,L)
#ifdef HTAP_LIKE_DIAGS
            aijl(i,j,L,ijl_gridH)=aijl(i,j,L,ijl_gridH)+
     &      rgas/grav*TX(i,j,L)*log(pedn(l,i,j)/pedn(L+1,i,j))
#endif
            aijl(i,j,L,ijl_zL)=aijl(i,j,L,ijl_zL)+phi(i,j,l)/grav
            DO IT=1,NTYPE
              CALL INC_AJ(I,J,IT,J_TX,(TX(I,J,L)-TF)*FTYPE(IT,I,J)*
     *             DBYSD)
              CALL INC_AJ(I,J,IT,J_QP,(Q(I,J,L)+WM(I,J,L))*PIJ*DSIG(L)
     *             *FTYPE(IT,I,J))
            END DO
            CALL INC_AREG(I,J,JR,J_QP,(Q(I,J,L)+WM(I,J,L))*PIJ*DSIG(L))
            CALL INC_AREG(I,J,JR,J_TX,(TX(I,J,L)-TF)*DBYSD)
          END DO
        END DO
      END DO


C****
C**** STATIC STABILITIES: TROPOSPHERIC AND STRATOSPHERIC
C****
      DO J=J_0,J_1
C**** OLD TROPOSPHERIC STATIC STABILITY
        DO I=I_0,IMAXJ(J)
          DXYPJ=AXYP(I,J)
          JR=JREG(I,J)
          SS=(T(I,J,LS1-1)-T(I,J,1))/(PHI(I,J,LS1-1)-PHI(I,J,1)+teeny)
          DO IT=1,NTYPE
            CALL INC_AJ(I,J,IT,J_DTDGTR,SS*FTYPE(IT,I,J))
          END DO
          CALL INC_AREG(I,J,JR,J_DTDGTR,SS)
          AIJ(I,J,IJ_DTDP)=AIJ(I,J,IJ_DTDP)+SS
        END DO
C**** OLD STRATOSPHERIC STATIC STABILITY (USE LSTR as approx 10mb)
        DO I=I_0,IMAXJ(J)
          JR=JREG(I,J)
          SS=(T(I,J,LSTR)-T(I,J,LS1-1))/((PHI(I,J,LSTR)-PHI(I,J,LS1-1))
     *         +teeny)
          DO IT=1,NTYPE
            CALL INC_AJ(I,J,IT,J_DTSGST,SS*FTYPE(IT,I,J))
          END DO
          CALL INC_AREG(I,J,JR,J_DTSGST,SS)
        END DO

C****
C**** NUMBERS ACCUMULATED FOR THE RADIATION EQUILIBRIUM LAYERS
C****
        DO I=I_0,IMAXJ(J)
          DO LR=1,LM_REQ
            TRI(LR)=RQT(LR,I,J)
            call inc_asjl(i,j,lr,1,RQT(LR,I,J)-TF)
          END DO
          PHI_REQ=PHI(I,J,LM)
          PHI_REQ=PHI_REQ+RGAS*.5*(TX(I,J,LM)+TRI(1))*DLNP01
          call inc_asjl(i,j,1,2,PHI_REQ)
          PHI_REQ=PHI_REQ+RGAS*.5*(TRI(1)+TRI(2))*DLNP12
          call inc_asjl(i,j,2,2,PHI_REQ)
          PHI_REQ=PHI_REQ+RGAS*.5*(TRI(2)+TRI(3))*DLNP23
          call inc_asjl(i,j,3,2,PHI_REQ)
        END DO
      END DO

C****
C**** RICHARDSON NUMBER , ROSSBY NUMBER , RADIUS OF DEFORMATION
C****

c
c Spatial mean Rossby number is computed using the spatial mean of
c column-maximum wind speeds.
c
c Spatial mean Rossby radius is computed using the spatial mean of
c min(N*H/f,sqrt(N*H/beta)) where H is the either the depth of
c the troposphere or the lower stratosphere. (N*H)**2 is taken as
c log(theta2/theta1)*grav*H where theta2,theta1 are the potential
c temperatures at the top/bottom of the troposphere or lower stratosphere.
c f is the coriolis parameter and beta its latitudinal derivative.
c
c The Richardson number is not being computed.  It is set to 99 or 999
c to flag that it is missing.  For the large vertical ranges of these
c calculations, it should be replaced by something like the ratio of
c available mechanical energy to the work required to homogenize
c potential temperature.
c

      X=RGAS*LHE*LHE/(SHA*RVAP)
      bybetafac = radius/(2.*omega)

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)

        byfcor = 1d0/(2.*omega*abs(sinlat2d(i,j))+teeny)
        bybeta = bybetafac/(coslat2d(i,j)+teeny)
c troposphere
        w2max = maxval(ua(1:ls1-1,i,j)**2+va(1:ls1-1,i,j)**2)+teeny
        rossn = sqrt(w2max)*byfcor
        nh = sqrt((phi(i,j,ls1-1)-phi(i,j,1))*
     &       log(t(i,j,ls1-1)/t(i,j,1)))
        rossl = min(nh*byfcor,sqrt(nh*bybeta))
        richn = 99d0 ! for now

c lapse rates
        gamx = (tx(i,j,1)-tx(i,j,ls1-1))/
     &       (phi(i,j,ls1-1)-phi(i,j,1))
        gamm=0.
        do l=1,ls1-1
          tzl=tx(i,j,l)
          prt=(sig(l)*p(i,j)+ptop)*rgas*tzl
          eseps=qsat(tzl,lhe,one)
          gamm=gamm+dsig(l)*(prt+lhe*eseps)/(prt+x*eseps/tzl)
        end do

        do it=1,ntype
          call inc_aj(i,j,it,j_rostr,rossn*ftype(it,i,j))
          call inc_aj(i,j,it,j_ltro,rossl*ftype(it,i,j))
          call inc_aj(i,j,it,j_rictr,richn*ftype(it,i,j))

          call inc_aj(i,j,it,j_gam ,gamx*ftype(it,i,j))
          call inc_aj(i,j,it,j_gamm,gamm*ftype(it,i,j))

        end do

c stratosphere
        w2max = maxval(ua(ls1:lstr,i,j)**2+va(ls1:lstr,i,j)**2)
        rossn = sqrt(w2max)*byfcor
        nh = sqrt((phi(i,j,lstr)-phi(i,j,ls1))*
     &       log(t(i,j,lstr)/t(i,j,ls1)))
        rossl = min(nh*byfcor,sqrt(nh*bybeta))
        richn = 999d0 ! for now
        do it=1,ntype
          call inc_aj(i,j,it,j_rosst,rossn*ftype(it,i,j))
          call inc_aj(i,j,it,j_lstr,rossl*ftype(it,i,j))
          call inc_aj(i,j,it,j_ricst,richn*ftype(it,i,j))
        end do

      ENDDO
      ENDDO

C****
C**** NORTHWARD GRADIENT OF TEMPERATURE: TROPOSPHERIC AND STRATOSPHERIC
C**** GAMC, THE DYNAMICALLY DETERMINED LAPSE RATE IN THE EXTRATROPICS
C****
      do j=j_0,j_1
      do i=i_0,i_1
        tx_trop(i,j) = 0.
        tx_strat(i,j) = 0.
      enddo
      enddo
      do l=1,ls1-1
      do j=j_0,j_1
      do i=i_0,i_1
        tx_trop(i,j) = tx_trop(i,j) + tx(i,j,l)*dsig(l)
      enddo
      enddo
      enddo
      do l=ls1,lstr
      do j=j_0,j_1
      do i=i_0,i_1
        tx_strat(i,j) = tx_strat(i,j) + tx(i,j,l)*dsig(l)
      enddo
      enddo
      enddo
      do j=j_0,j_1
      do i=i_0,i_1
        tx_trop(i,j) = tx_trop(i,j)/(SIGE(1)-SIGE(LS1))
        tx_strat(i,j) = tx_strat(i,j)/(SIGE(LS1)-SIGE(LSTR+1)+1d-12)
      enddo
      enddo

      call halo_update(grid, tx_trop)
      call halo_update(grid, tx_strat)

      if(i_0h.lt.i_0) then      ! halo cells exist in i direction
        im1s=i_0-1; ip1s=i_0+1; ip1e=i_1+1
      else                      ! periodic
        im1s=im-1; ip1s=1; ip1e=im
      endif
      do j=j_0s,j_1s
        im1=im1s
        i=im1s+1
        do ip1=ip1s,ip1e
          dtxdy = (tx_trop(ip1,j)-tx_trop(im1,j))*ddy_ci(i,j)
     &          + (tx_trop(i,j+1)-tx_trop(i,j-1))*ddy_cj(i,j)
          gamc = gamd+grav*radius*dtxdy*sinlat2d(i,j)/
     &         (rgas*tx_trop(i,j)*(coslat2d(i,j)+.001))
          do it=1,ntype
            call inc_aj(i,j,it,j_dtdjt,dtxdy*ftype(it,i,j))
            call inc_aj(i,j,it,j_gamc,gamc*ftype(it,i,j))
          enddo
          dtxdy = (tx_strat(ip1,j)-tx_strat(im1,j))*ddy_ci(i,j)
     &          + (tx_strat(i,j+1)-tx_strat(i,j-1))*ddy_cj(i,j)
          do it=1,ntype
            call inc_aj(i,j,it,j_dtdjs,dtxdy*ftype(it,i,j))
          enddo
          im1=i
          i=ip1
        enddo
      enddo

      DO J=J_0S,J_1S
      DO I=I_0,I_1
        AIJ(I,J,IJ_UJET)=AIJ(I,J,IJ_UJET)+UA(JET,I,J)
        AIJ(I,J,IJ_VJET)=AIJ(I,J,IJ_VJET)+VA(JET,I,J)
      ENDDO
      ENDDO

C****
C**** CONVERT VERTICAL WINDS TO UNITS PROPORTIONAL TO M/S
C****
      DO L=1,LM-1
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        PIJ=PLIJ(L,I,J)
        PE=SIGE(L+1)*PIJ+PTOP
        PKE=PE**KAPA
        THETA=THBAR(T(I,J,L+1),T(I,J,L))
        W(I,J,L)=SD(I,J,L)*THETA*PKE/PE
      END DO
       if(have_south_pole .and. J==1) W(2:IM,J,L) = W(1,J,L)
       if(have_north_pole .and. J==JM) W(2:IM,J,L) = W(1,J,L)
      END DO
      END DO

c
c accumulate AIJL: U,V,W on model layers
c

      DO L=1,LM
      DO J=J_0S,J_1S
      DO I=I_0,I_1
        AIJL(I,J,L,IJL_U) = AIJL(I,J,L,IJL_U) + UA(L,I,J)
        AIJL(I,J,L,IJL_V) = AIJL(I,J,L,IJL_V) + VA(L,I,J)
      ENDDO
      ENDDO
      ENDDO ! L
      DO L=1,LM-1
      DO J=J_0,J_1
      DO I=I_0,I_1
        AIJL(I,J,L,IJL_W) = AIJL(I,J,L,IJL_W) + BYAXYP(I,J)*
#ifdef CUBED_SPHERE
     &       WCP(I,J,L)
#else
     &       W(I,J,L)
#endif
      ENDDO
      ENDDO
      ENDDO

C****
C**** CERTAIN HORIZONTAL WIND AVERAGES
C****
      do j=j_0,j_1
      do i=i_0,imaxj(j)

        if(lon2d_dg(i,j).ge.-135. .and. lon2d_dg(i,j).le.-110.) then
c east pacific
          do l=1,lm
            call inc_ajl(i,j,l,jl_epacwt,pdsig(l,i,j))
            call inc_ajl(i,j,l,jl_uepac, pdsig(l,i,j)*ua(l,i,j))
            call inc_ajl(i,j,l,jl_vepac, pdsig(l,i,j)*va(l,i,j))
          enddo
          do l=1,lm-1 ! using pdsig for vertical velocity weight
            call inc_ajl(i,j,l,jl_wepac,
     &           pdsig(l,i,j)*w(i,j,l)*byaxyp(i,j))
          enddo
        elseif(lon2d_dg(i,j).ge.150.) then
c west pacific
          do l=1,lm
            call inc_ajl(i,j,l,jl_wpacwt,pdsig(l,i,j))
            call inc_ajl(i,j,l,jl_uwpac, pdsig(l,i,j)*ua(l,i,j))
            call inc_ajl(i,j,l,jl_vwpac, pdsig(l,i,j)*va(l,i,j))
          enddo
          do l=1,lm-1 ! using pdsig for vertical velocity weight
            call inc_ajl(i,j,l,jl_wwpac,
     &           pdsig(l,i,j)*w(i,j,l)*byaxyp(i,j))
          enddo
        endif

      enddo
      enddo

c
c constant-pressure diagnostics
c
      pecp(2:lm+1) = ple(1:lm)
      pecp(1) = 1d30 ! ensure that all column mass is included
      qlh=lhe
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        pedge(:) = pedn(:,i,j)
        call get_dx_intervals(pedge,lm,pecp,lm,dpx,lmod,lcp,lmx,lmxmax)
        do l=1,lm
          dpwt(l) = 0d0
          txdp(l) = 0d0
          phidp(l) = 0d0
          qdp(l) = 0d0
          rhdp(l) = 0d0
          wmdp(l) = 0d0
          wmliqdp(l) = 0d0
          wmfrzdp(l) = 0d0
          rh(l) = q(i,j,l)/min(1d0,QSAT(TX(I,J,L),QLH,pmid(l,i,j)))
        enddo
        do l=1,lmx
          dpwt(lcp(l))  = dpwt(lcp(l))  + dpx(l)
          txdp(lcp(l))  = txdp(lcp(l))  + dpx(l)*(tx(i,j,lmod(l))-tf)
          phidp(lcp(l)) = phidp(lcp(l)) + dpx(l)*phi(i,j,lmod(l))
          qdp(lcp(l))   = qdp(lcp(l))   + dpx(l)*q(i,j,lmod(l))
          rhdp(lcp(l))  = rhdp(lcp(l))  + dpx(l)*rh(lmod(l))
          wmdp(lcp(l))  = wmdp(lcp(l))  + dpx(l)*wm(i,j,lmod(l))
          if( svlhx(lmod(l),i,j) == lhe)
     *      wmliqdp(lcp(l)) = wmliqdp(lcp(l)) + dpx(l)*wm(i,j,lmod(l))
          if( svlhx(lmod(l),i,j) == lhs)
     *      wmfrzdp(lcp(l)) = wmfrzdp(lcp(l)) + dpx(l)*wm(i,j,lmod(l))
        enddo
        do l=1,lm
          aijl(i,j,l,ijk_dp) = aijl(i,j,l,ijk_dp) + dpwt(l)
          aijl(i,j,l,ijk_tx) = aijl(i,j,l,ijk_tx) + txdp(l)
          aijl(i,j,l,ijk_q)  = aijl(i,j,l,ijk_q)  + qdp(l)
          aijl(i,j,l,ijk_rh) = aijl(i,j,l,ijk_rh) + rhdp(l)
          call inc_ajl(i,j,l,jk_dpwt,  dpwt(l))
          call inc_ajl(i,j,l,jk_tx,    txdp(l))
          call inc_ajl(i,j,l,jk_hght,  phidp(l))
          call inc_ajl(i,j,l,jk_q,     qdp(l))
          call inc_ajl(i,j,l,jk_rh,    rhdp(l))
          call inc_ajl(i,j,l,jk_cldh2o,wmdp(l))
          call inc_ajl(i,j,l,jk_cldice,wmfrzdp(l))
          call inc_ajl(i,j,l,jk_cldwtr,wmliqdp(l))
        enddo
      enddo
      enddo

C**** ACCUMULATE TIME USED IN DIAGA
      CALL TIMEOUT(BEGIN,MDIAG,MDYN)
      RETURN

      ENTRY DIAGA0
c increment ajl(jl_dtdyn) by -t before dynamics.
c ajl(jl_dtdyn) will be incremented by +t after the dynamics, giving
c the tendency.

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

      DO L=1,LM
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        call inc_ajl(i,j,l,jl_dtdyn,-t(i,j,l)*pk(l,i,j)*pdsig(l,i,j))
      END DO
      END DO
      END DO
      RETURN
C****
      END SUBROUTINE DIAGA

      subroutine get_dx_intervals(
     &     xesrc,nsrc,xedst,ndst,dx,lsrc,ldst,nxchng,nmax)
c
c This routine returns information needed for conservative
c remapping in 1 dimension.
c Given two lists of points xesrc(1:nsrc+1) and xedst(1:ndst+1)
c along the "x" axis, merge the two lists and calculate the
c distance increments dx(1:nxchng) separating the points
c in the merged list.  The i_th increment lies in the intervals
c xesrc(lsrc(i):lsrc(i)+1) and xedst(ldst(i):ldst(i)+1).
c Points in xedst lying outside the range xesrc(1:nsrc+1) are
c not included.
c
      implicit none
      integer :: nsrc,ndst,nxchng,nmax
      real*8, dimension(nsrc+1) :: xesrc
      real*8, dimension(ndst+1) :: xedst
      real*8, dimension(nmax) :: dx
      integer, dimension(nmax) :: lsrc,ldst
      integer :: ls,ld
      real*8 :: xlast,s
      s = sign(1d0,xesrc(2)-xesrc(1))
      ld = 1
      xlast = max(s*xesrc(1),s*xedst(1))
      do while(s*xedst(ld).le.xlast)
        ld = ld + 1
      enddo
      ls = 2
      do nxchng=1,nmax
        lsrc(nxchng) = ls-1
        ldst(nxchng) = ld-1
        if(s*xedst(ld).le.s*xesrc(ls)) then
          if(xedst(ld).eq.xesrc(ls)) ls = ls + 1
          dx(nxchng) = s*xedst(ld)-xlast
          xlast = s*xedst(ld)
          ld = ld + 1
        else
          dx(nxchng) = s*xesrc(ls)-xlast
          xlast = s*xesrc(ls)
          ls = ls + 1
        endif
        if(ld.gt.ndst+1) exit
        if(ls.gt.nsrc+1) exit
      enddo
      return
      end subroutine get_dx_intervals

      SUBROUTINE DIAGCA (M)
!@sum  DIAGCA Keeps track of the conservation properties of angular
!@+    momentum, kinetic energy, mass, total potential energy and water
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : mdiag,itime
#ifdef TRACERS_ON
      USE TRACER_COM, only: itime_tr0,ntm  !xcon
#endif
      USE DIAG_COM, only : icon_AM,icon_KE,icon_MS,icon_TPE
     *     ,icon_WM,icon_LKM,icon_LKE,icon_EWM,icon_WTG,icon_HTG
     *     ,icon_OMSI,icon_OHSI,icon_OSSI,icon_LMSI,icon_LHSI,icon_MLI
     *     ,icon_HLI,icon_MICB,icon_HICB,title_con
      !USE SOIL_DRV, only: conserv_WTG,conserv_HTG
      IMPLICIT NONE
!@var M index denoting from where DIAGCA is called
      INTEGER, INTENT(IN) :: M
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1  INITIALIZE CURRENT QUANTITY
C****   2  AFTER DYNAMICS
C****   3  AFTER CONDENSATION
C****   4  AFTER RADIATION
C****   5  AFTER PRECIPITATION
C****   6  AFTER LAND SURFACE (INCL. RIVER RUNOFF)
C****   7  AFTER FULL SURFACE INTERACTION
C****   8  AFTER FILTER
C****   9  AFTER OCEAN DYNAMICS (from MAIN)
C****  10  AFTER DAILY
C****  11  AFTER OCEAN DYNAMICS (from ODYNAM)
C****  12  AFTER OCEAN SUB-GRIDSCALE PHYS
C****
#ifndef SCM
      EXTERNAL conserv_AM,conserv_KE,conserv_MS,conserv_PE
     *     ,conserv_WM,conserv_EWM,conserv_LKM,conserv_LKE,conserv_OMSI
     *     ,conserv_OHSI,conserv_OSSI,conserv_LMSI,conserv_LHSI
     *     ,conserv_MLI,conserv_HLI,conserv_WTG,conserv_HTG
     *     ,conserv_MICB,conserv_HICB
      real*8 NOW
      INTEGER NT

C**** ATMOSPHERIC ANGULAR MOMENTUM
      CALL conserv_DIAG(M,conserv_AM,icon_AM)

C**** ATMOSPHERIC KINETIC ENERGY
      CALL conserv_DIAG(M,conserv_KE,icon_KE)

C**** ATMOSPHERIC MASS
      CALL conserv_DIAG(M,conserv_MS,icon_MS)

C**** ATMOSPHERIC TOTAL POTENTIAL ENERGY
      CALL conserv_DIAG(M,conserv_PE,icon_TPE)

C**** ATMOSPHERIC TOTAL WATER MASS
      CALL conserv_DIAG(M,conserv_WM,icon_WM)

C**** ATMOSPHERIC TOTAL WATER ENERGY
      CALL conserv_DIAG(M,conserv_EWM,icon_EWM)

C**** LAKE MASS AND ENERGY
      CALL conserv_DIAG(M,conserv_LKM,icon_LKM)
      CALL conserv_DIAG(M,conserv_LKE,icon_LKE)

C**** OCEAN ICE MASS, ENERGY, SALT
      CALL conserv_DIAG(M,conserv_OMSI,icon_OMSI)
      CALL conserv_DIAG(M,conserv_OHSI,icon_OHSI)
      CALL conserv_DIAG(M,conserv_OSSI,icon_OSSI)

C**** LAKE ICE MASS, ENERGY
      CALL conserv_DIAG(M,conserv_LMSI,icon_LMSI)
      CALL conserv_DIAG(M,conserv_LHSI,icon_LHSI)

C**** GROUND WATER AND ENERGY
      CALL conserv_DIAG(M,conserv_WTG,icon_WTG)
      CALL conserv_DIAG(M,conserv_HTG,icon_HTG)

C**** LAND ICE MASS AND ENERGY
      CALL conserv_DIAG(M,conserv_MLI,icon_MLI)
      CALL conserv_DIAG(M,conserv_HLI,icon_HLI)

C**** ICEBERG MASS AND ENERGY
      CALL conserv_DIAG(M,conserv_MICB,icon_MICB)
      CALL conserv_DIAG(M,conserv_HICB,icon_HICB)

C**** OCEAN CALLS ARE DEALT WITH SEPARATELY
      CALL DIAGCO (M)

#ifdef TRACERS_ON
C**** Tracer calls are dealt with separately
      do nt=1,ntm
        CALL DIAGTCA(M,NT)
      end do
#endif
C****
      CALL TIMER (NOW,MDIAG)

#endif /* not SCM */

      RETURN
      END SUBROUTINE DIAGCA

      SUBROUTINE conserv_DIAG (M,CONSFN,ICON)
!@sum  conserv_DIAG generic routine keeps track of conserved properties
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE GEOM, only : j_budg, j_0b, j_1b, imaxj
      USE DIAG_COM, only : consrv=>consrv_loc,nofm, jm_budg,wtbudg
      USE DOMAIN_DECOMP_ATM, only : GET, GRID
      IMPLICIT NONE
!@var M index denoting from where routine is called
      INTEGER, INTENT(IN) :: M
!@var ICON index for the quantity concerned
      INTEGER, INTENT(IN) :: ICON
!@var CONSFN external routine that calculates total conserved quantity
      EXTERNAL CONSFN
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: TOTAL
      REAL*8, DIMENSION(JM_BUDG) :: TOTALJ
      INTEGER :: I,J,NM,NI
      INTEGER :: I_0,I_1,J_0,J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C**** NOFM contains the indexes of the CONSRV array where each
C**** change is to be stored for each quantity. If NOFM(M,ICON)=0,
C**** no calculation is done.
C**** NOFM(1,ICON) is the index for the instantaneous value.
      IF (NOFM(M,ICON).gt.0) THEN
C**** Calculate current value TOTAL
        CALL CONSFN(TOTAL)
        NM=NOFM(M,ICON)
        NI=NOFM(1,ICON)
C**** Calculate zonal sums
        TOTALJ(J_0B:J_1B)=0.
        DO J=J_0,J_1
          DO I=I_0,IMAXJ(J) ! not I_1 b/c latlon wtbudg differs at poles
            TOTALJ(J_BUDG(I,J)) = TOTALJ(J_BUDG(I,J)) + TOTAL(I,J)
     &           *WTBUDG(I,J)
          END DO
        END DO
C**** Accumulate difference from last time in CONSRV(NM)
        IF (M.GT.1) THEN
          DO J=J_0B,J_1B
            CONSRV(J,NM)=CONSRV(J,NM)+(TOTALJ(J)-CONSRV(J,NI))
          END DO
        END IF
C**** Save current value in CONSRV(NI)
        DO J=J_0B,J_1B
          CONSRV(J,NI)=TOTALJ(J)
        END DO
      END IF
      RETURN
C****
      END SUBROUTINE conserv_DIAG


      SUBROUTINE conserv_MS(RMASS)
!@sum  conserv_MA calculates total atmospheric mass
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,p,pstrat
      USE GEOM, only : imaxj
      USE DOMAIN_DECOMP_ATM, only : GET, GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: RMASS
      INTEGER :: I,J
      INTEGER :: J_0,J_1 ,I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(grid, J_STRT=J_0,    J_STOP=J_1,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****
C**** MASS
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        RMASS(I,J)=(P(I,J)+PSTRAT)*mb2kg
      END DO
      END DO
      IF(HAVE_SOUTH_POLE) RMASS(2:im,1) =RMASS(1,1)
      IF(HAVE_NORTH_POLE) RMASS(2:im,JM)=RMASS(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_MS


      SUBROUTINE conserv_PE(TPE)
!@sum  conserv_TPE calculates total atmospheric potential energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sha,mb2kg
      USE MODEL_COM, only : im,jm,lm,t,p,ptop,zatmo
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pk,pdsig
      USE DOMAIN_DECOMP_ATM, only : GET,GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: TPE
      INTEGER :: I,J,L
      INTEGER :: J_0,J_1,I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****
C**** TOTAL POTENTIAL ENERGY (J/m^2)
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        TPE(I,J)=0.
        DO L=1,LM
          TPE(I,J)=TPE(I,J)+T(I,J,L)*PK(L,I,J)*PDSIG(L,I,J)
        ENDDO
        TPE(I,J)=(TPE(I,J)*SHA+ZATMO(I,J)*(P(I,J)+PTOP))*mb2kg
      ENDDO
      ENDDO
      IF(HAVE_SOUTH_POLE) TPE(2:im,1) =TPE(1,1)
      IF(HAVE_NORTH_POLE) TPE(2:im,JM)=TPE(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_PE

      SUBROUTINE conserv_WM(WATER)
!@sum  conserv_WM calculates total atmospheric water mass
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,lm,wm,q
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pdsig
      USE DOMAIN_DECOMP_ATM, only : GET, GRID
      IMPLICIT NONE

      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: WATER
      INTEGER :: I,J,L
      INTEGER :: J_0,J_1,I_0,I_1
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

      CALL GET(GRID, J_STRT=J_0, J_STOP=J_1,
     &     HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****
C**** TOTAL WATER MASS (kg/m^2)
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        WATER(I,J) = 0.
        DO L=1,LM
          WATER(I,J)=WATER(I,J)+(Q(I,J,L)+WM(I,J,L))*PDSIG(L,I,J)
        ENDDO
        WATER(I,J)=WATER(I,J)*mb2kg
      ENDDO
      ENDDO
      IF (HAVE_SOUTH_POLE) WATER(2:im,1) = WATER(1,1)
      IF (HAVE_NORTH_POLE) WATER(2:im,JM)= WATER(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_WM


      SUBROUTINE conserv_EWM(EWATER)
!@sum  conserv_EWM calculates total atmospheric water energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg,shv,grav,lhe
      USE MODEL_COM, only : im,jm,lm,wm,t,q,p
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pdsig, pmid, pk
      USE CLOUDS_COM, only : svlhx
      USE DOMAIN_DECOMP_ATM, only : GET, GRID
      IMPLICIT NONE
      REAL*8, PARAMETER :: HSCALE = 7.8d0 ! km
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: EWATER
      INTEGER :: I,J,L
      INTEGER :: J_0,J_1,I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE
      REAL*8 EL!,W

      CALL GET(GRID, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****
C**** TOTAL WATER ENERGY (J/m^2)
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        EWATER(I,J) = 0.
        DO L=1,LM
c this calculation currently only calculates latent heat
c          W =(Q(I,J,L)+WM(I,J,L))*PDSIG(L,I,J)*mb2kg
          EL=(Q(I,J,L)*LHE+WM(I,J,L)*(LHE-SVLHX(L,I,J)))*PDSIG(L,I,J)
          EWATER(I,J)=EWATER(I,J)+EL !+W*(SHV*T(I,J,L)*PK(L,I,J)+GRAV
!     *           *HSCALE*LOG(P(I,J)/PMID(L,I,J)))
        ENDDO
        EWATER(I,J)=EWATER(I,J)*mb2kg
      ENDDO
      ENDDO
      IF(HAVE_SOUTH_POLE) EWATER(2:im,1) = EWATER(1,1)
      IF(HAVE_NORTH_POLE) EWATER(2:im,JM)= EWATER(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_EWM

      SUBROUTINE DIAG4A
C****
C**** THIS ROUTINE PRODUCES A TIME HISTORY OF ENERGIES
C****
      USE MODEL_COM, only : im,istrat,IDACC
      USE DIAG_COM, only : energy,speca,ned
      IMPLICIT NONE

      INTEGER :: I,IDACC5,N,NM

#ifdef SCM
      return
#endif

      IF (IDACC(4).LE.0.OR.IDACC(7).LE.0) RETURN
      NM=1+IM/2
C****
C**** LOAD ENERGIES INTO TIME HISTORY ARRAY
C****
      IDACC5=IDACC(5)+1
      IF (IDACC5.GT.100) RETURN
      DO I=0,1+ISTRAT  ! loop over number of 'spheres'
        ENERGY(1+NED*I,IDACC5)=SPECA(1,19,1+4*I)   ! SH
        ENERGY(2+NED*I,IDACC5)=SPECA(1,19,2+4*I)   ! NH
        ENERGY(5+NED*I,IDACC5)=SPECA(2,19,2+4*I)   ! NH wave 1
        ENERGY(6+NED*I,IDACC5)=SPECA(3,19,2+4*I)   ! NH wave 2
        ENERGY(7+NED*I,IDACC5)=SPECA(1,20,1+4*I)
        ENERGY(8+NED*I,IDACC5)=SPECA(1,20,2+4*I)
        DO N=2,NM
        ENERGY( 3+NED*I,IDACC5)=ENERGY( 3+10*I,IDACC5)+SPECA(N,19,1+4*I)
        ENERGY( 4+NED*I,IDACC5)=ENERGY( 4+10*I,IDACC5)+SPECA(N,19,2+4*I)
        ENERGY( 9+NED*I,IDACC5)=ENERGY( 9+10*I,IDACC5)+SPECA(N,20,1+4*I)
        ENERGY(10+NED*I,IDACC5)=ENERGY(10+10*I,IDACC5)+SPECA(N,20,2+4*I)
        END DO
      END DO
      IDACC(5)=IDACC5
      RETURN
C****
      END SUBROUTINE DIAG4A

#ifndef CACHED_SUBDD

      module subdaily
!@sum SUBDAILY defines variables associated with the sub-daily diags
!@auth Gavin Schmidt
      use domain_decomp_atm, only: get,grid,am_i_root
      USE MODEL_COM, only : im,jm,lm,itime,itime0,nday,iyear1,jyear
     &     ,jmon,jday,jdate,jhour,dtsrc,xlabel,jdpery,JDendOfM,lrunid
      USE FILEMANAGER, only : openunit, closeunit, nameunit
      use ghy_com, only: gdeep,gsaveL,ngm
      USE DIAG_COM, only : kgz_max,pmname,P_acc,PM_acc,R_acc
#ifdef TES_LIKE_DIAGS
     *                    ,kgz_max_more,pmnamemore
#endif
      USE PARAM
#ifdef CALCULATE_FLAMMABILITY
      use flammability_com, only : raP_acc
#endif
#ifdef TRACERS_ON
      use rad_com, only: nTracerRadiaActive,tracerRadiaActiveFlag
      use tracer_com
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only : sOx_acc,sNOx_acc,sCO_acc
     &     ,l1Ox_acc,l1NO2_acc
#endif
#ifdef TRACERS_ON
      use trdiag_com, only: trcsurf,trcSurfByVol,trcSurfMixR_acc
     &     ,trcSurfByVol_acc
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
     &     ,sPM2p5_acc,sPM10_acc,l1PM2p5_acc,l1PM10_acc
     &     ,csPM2p5_acc,csPM10_acc
#endif
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only : BE7D_acc,BE7W_acc
#endif
#endif
#ifdef TRACERS_WATER
      USE TRDIAG_COM, only : trp_acc, tre_acc
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      use tracers_dust, only: dustDiagSubdd_acc,dust_names,n_soilDust
      use trdust_drv, only: accSubddDust
#endif
      IMPLICIT NONE
      SAVE
!@var kddmax maximum number of sub-daily diags output files
      INTEGER, PARAMETER :: kddmax = 55
!@var kdd total number of sub-daily diags
      INTEGER :: kdd
!@var kddunit total number of sub-daily files
      INTEGER :: kddunit
!@var namedd array of names of sub-daily diags
      CHARACTER*10, DIMENSION(kddmax) :: namedd
!@var iu_subdd array of unit numbers for sub-daily diags output
      INTEGER, DIMENSION(kddmax) :: iu_subdd
!@var subddt = subdd + subdd1,2,3 = all variables for sub-daily diags
      CHARACTER*320 :: subddt = " "
!@dbparam subdd string contains variables to save for sub-daily diags
!@dbparam subdd1 additional string of variables for sub-daily diags
!@dbparam subdd2 additional string of variables for sub-daily diags
!@dbparam subdd3 additional string of variables for sub-daily diags
!@dbparam subdd4 additional string of variables for sub-daily diags
C**** Note: for longer string increase MAX_CHAR_LENGTH in PARAM
      CHARACTER*64 :: subdd="SLP", 
     & subdd1=" ", subdd2=" ", subdd3=" ", subdd4=" "
!@dbparam Nsubdd: DT_save_SUBDD =  Nsubdd*DTsrc sub-daily diag freq.
      INTEGER :: Nsubdd = 0
!@dbparam LmaxSUBDD: the max L when writing "ALL" levels
      INTEGER :: LmaxSUBDD = LM
!@var lst level strings
      character*2, dimension(lm) :: lst

      character*14, private :: adate_sv

!@var LmaxSUBDD_array array for three-dimensional fields for subdd diagnostics
!@var ngm_array three-dimensional array for subdd ground diagnostics
      real(kind=8),allocatable,dimension(:,:,:) :: LmaxSUBDD_array
     &     ,ngm_array
!@var kgz_max_suffixes array of names for subdd-diagnostic on pressure levels
      character(len=8),allocatable,dimension(:) :: kgz_max_suffixes
!@var kgz_max_array three-dimensional array for diagnostics on pressure levels
      real(kind=8),allocatable,dimension(:,:,:) :: kgz_max_array
#ifdef TES_LIKE_DIAGS
!@var kgz_max_suffixes array of names for subdd-diagnostic on pressure levels
      character(len=8),allocatable,dimension(:) :: kgz_max_more_suffixes
!@var kgz_max_array three-dimensional array for diagnostics on more pressure levels
      real(kind=8),allocatable,dimension(:,:,:) :: kgz_max_more_array
#endif
#ifdef TRACERS_ON
!@var rTrname array with tracer names for subdd radiation diagnostics
      character(len=len(trname(1))),allocatable,dimension(:) :: rTrname
!@var TRACER_array tracer array for subdd diagnostics
!@var rTRACER_array tracer array for subdd radiation diagnostic
      real(kind=8),allocatable,dimension(:,:,:) :: TRACER_array
     &     ,rTRACER_array
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
!@var dust3d_array three-dimensional soil dust array for subdd diagnostics
!@var dust4d_array four-dimensional soil dust array for subdd diagnostics
      real(kind=8),allocatable,dimension(:,:,:) :: dust3d_array
      real(kind=8),allocatable,dimension(:,:,:,:) :: dust4d_array
#endif

#ifdef NEW_IO_SUBDD
      private :: write_2d,write_3d,write_4d,in_subdd_list
     &     ,def_global_attr_subdd,def_xy_coord_subdd,time_subdd
     &     ,get_calendarstring,get_referencetime_for_netcdf
     &     ,def_time_coord_subdd,write_xy_coord_subdd
     &     ,write_time_coord_subdd
      interface write_subdd
      module procedure write_2d
      module procedure write_3d
      module procedure write_4d
      end interface
#endif

      contains

      subroutine init_subdd(aDATE)
!@sum init_subdd initialise sub daily diags and position files
!@auth Gavin Schmidt
      implicit none
      character*14, intent(in) :: adate
      integer :: i,j,k,l,kunit,i_0h,i_1h,j_0h,j_1h

      call get(grid,i_strt_halo=i_0h,i_stop_halo=i_1h,j_strt_halo=j_0h
     &     ,j_stop_halo=j_1h)

      adate_sv = adate

      call sync_param( "subdd" ,subdd)
      call sync_param( "subdd1" ,subdd1)
      call sync_param( "subdd2" ,subdd2)
      call sync_param( "subdd3" ,subdd3)
      call sync_param( "subdd4" ,subdd4)
      call sync_param( "Nsubdd",Nsubdd)
      call sync_param( "LmaxSUBDD",LmaxSUBDD)

      if (nsubdd.ne.0) then
C**** combine strings subdd, subdd1...4:
        subddt=trim(subdd)//' '//
     &  trim(subdd1)//' '//trim(subdd2)//' '//trim(subdd3)//' '//subdd4
C**** calculate how many names
        k=0
        i=1
 10     j=index(subddt(i:len(subddt))," ")
        if (j.gt.1) then
          k=k+1
          i=i+j
        else
          i=i+1
        end if
        if (i.lt.len(subddt)) goto 10
        kdd=k
        if (kdd.gt.kddmax) call stop_model
     *       ("Increase kddmax: No. of sub-daily diags too big",255)

C**** make array of names
        read(subddt,*) namedd(1:kdd)

#ifndef NEW_IO_SUBDD
C**** open units and position
        call open_subdd(aDATE)

C**** position correctly
        do kunit=1,kddunit
          call io_POS(iu_SUBDD(kunit),Itime,im*jm,Nsubdd)
        end do
#endif

      end if

C**** define lst
      do l=1,lm
        if (l.lt.10) write(lst(l)(1:2),'(I1,1X)') l
        if (l.ge.10) write(lst(l)(1:2),'(I2)') l
      end do

C**** initialise special subdd accumulation
#ifdef TRACERS_COSMO
      BE7W_acc=0.
      BE7D_acc=0.
#endif
#ifdef TRACERS_WATER
      TRP_acc=0.
      TRE_acc=0.
#endif

      allocate(LmaxSUBDD_array(i_0h:i_1h,j_0h:j_1h,LmaxSUBDD))
      allocate(ngm_array(i_0h:i_1h,j_0h:j_1h,ngm))
      allocate(kgz_max_suffixes(kgz_max))
      allocate(kgz_max_array(i_0h:i_1h,j_0h:j_1h,kgz_max))
#ifdef TES_LIKE_DIAGS
      allocate(kgz_max_more_suffixes(kgz_max_more))
      allocate(kgz_max_more_array(i_0h:i_1h,j_0h:j_1h,kgz_max_more))
#endif
#ifdef TRACERS_ON
      allocate(rTrname(nTracerRadiaActive))
      allocate(rTRACER_array(i_0h:i_1h,j_0h:j_1h,nTracerRadiaActive))
      allocate(TRACER_array(i_0h:i_1h,j_0h:j_1h,ntm))
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      allocate(dust3d_array(i_0h:i_1h,j_0h:j_1h,ntm_dust))
      allocate(dust4d_array(i_0h:i_1h,j_0h:j_1h,LmaxSUBDD,ntm_dust))
#endif

      return
      end subroutine init_subdd

      subroutine open_subdd(aDATE)
!@sum open_subdd opens sub daily diag files
!@auth Gavin Schmidt
      implicit none
      character*14, intent(in) :: adate
      character*12 name
      integer :: k,kunit,kk

      adate_sv = adate

#if (defined NEW_IO_SUBDD)
      return
#endif

      kunit=0
      do k=1,kdd
C**** Some names have more than one unit associated (i.e. "ZALL")
        if (namedd(k)(len_trim(namedd(k))-2:len_trim(namedd(k))).eq.
     *       "ALL") then
          select case (namedd(k)(1:1))
          case ("U","V","W","C","D","o","B","n","t","q","z","r","m","x")
            ! velocities/tracers on model layers
            kunit=kunit+1
            write(name,'(A1,A3,A7)') namedd(k)(1:1),'ALL',aDATE(1:7)
            call openunit(name,iu_SUBDD(kunit),.true.,.false.)
          case ("Z", "R") ! heights, rel hum PMB levels
            do kk=1,kgz_max
              kunit=kunit+1
              call openunit(namedd(k)(1:1)//trim(PMNAME(kk))//
     *             aDATE(1:7),iu_SUBDD(kunit),.true.,.false.)
            end do
          case ("T", "Q", ! temps, spec hum PMB levels
     &          "O", "X", "M", "N")! Ox, NOx, CO, NO2
#ifdef TES_LIKE_DIAGS
            kunit=kunit+1
            call openunit(namedd(k)(1:1)//'TES'//
     *      aDATE(1:7),iu_SUBDD(kunit),.true.,.false.)
#else
            do kk=1,kgz_max
              kunit=kunit+1
              call openunit(namedd(k)(1:1)//trim(PMNAME(kk))//
     *             aDATE(1:7),iu_SUBDD(kunit),.true.,.false.)
            end do
#endif
          end select
          select case (namedd(k)(1:2))
          case ("GT","GW","GI")
            ! soil variables on soil levels
            kunit=kunit+1
            write(name,'(A2,A3,A7)') namedd(k)(1:2),'ALL',aDATE(1:7)
            call openunit(name,iu_SUBDD(kunit),.true.,.false.)
          end select
        else                    ! single file per name
          kunit=kunit+1
          call openunit(trim(namedd(k))//aDATE(1:7),iu_SUBDD(kunit),
     *         .true.,.false.)
        endif
      end do
      kddunit=kunit
C****
      return
      end subroutine open_subdd

      subroutine reset_subdd(aDATE)
!@sum reset_subdd resets sub daily diag files
!@auth Gavin Schmidt
      implicit none
      character*14, intent(in) :: adate

      adate_sv = adate

#if (defined NEW_IO_SUBDD)
      return
#endif

      if (nsubdd.ne.0) then
C**** close and re-open units
        call closeunit ( iu_SUBDD(1:kddunit) )
        call open_subdd( aDATE )
      end if
C****
      return
      end subroutine reset_subdd

c accSubdd
      subroutine accSubdd
!@sum  accSubdd accumulates variables for subdaily diagnostics
!@auth Jan Perlwitz

      implicit none

      integer :: i_0,i_1,j_0,j_1,i,j,n

      call get(grid,i_strt=i_0,i_stop=i_1,j_strt=j_0,j_stop=j_1)

#ifdef TRACERS_ON
!$OMP PARALLEL DO PRIVATE(i,j,n)
      do n=1,ntm
        do j=j_0,j_1
          do i=i_0,i_1
            trcSurfMixR_acc(i,j,n)=trcSurfMixR_acc(i,j,n)+trcSurf(i,j,n)
            trcSurfByVol_acc(i,j,n)=trcSurfByVol_acc(i,j,n)
     &           +trcSurfByVol(i,j,n)
          end do
        end do
      end do
!$OMP END PARALLEL DO
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      call accSubddDust(dustDiagSubdd_acc) ! in TRDUST_DRV.f
#endif

      return
      end subroutine accSubdd

c get_subdd
      subroutine get_subdd
!@sum get_SUBDD saves variables at sub-daily frequency
!@+   every ABS(NSUBDD)
!@+   Note that TMIN,TMAX,{ ,c,t,ct}AOD, are only output once/day.
!@+   If there is a choice between outputting pressure levels or
!@+   model levels, use lower case for the model levels:
!@+   Current options: SLP, PS, SAT, PREC, QS, LCLD, MCLD, HCLD, PTRO
!@+                    QLAT, QSEN, SWD, SWU, LWD, LWU, LWT, SWT, STX, STY,
!@+                    ICEF, SNOWD, TCLD, SST, SIT, US, VS, TMIN, TMAX
!@+                    MCP, SNOWC, RS, GT1, GTD, GW0, GWD, GI0, GID,
!@+                    GTALL, GWALL, GIALL (on soil levels)
!@+                    {L,M,H}CLDI,CTPI,TAUI (ISCCP quantities)
!@+                    LGTN, c2gLGTN (lightning flashs/cloud-to-ground)
!@+                    TRP*, TRE* (water tracers only)
!@+                    Z*, R*, T*, Q*  (on any fixed pressure level)
!@+                    z*, r*, t*, q*  (on any model level)
!@+                    U*, V*, W*, C*  (on any model level only)
!@+                    O*, X*, M*, N*  (Ox,NOx,CO,NO2 on fixed pres lvl)
!@+                    o*, x*, m*, n*  (Ox,NOx,CO,NO2 on any model lvl)
!@+                    oAVG  (SFC Ox time-average ppbv)
!@+                    nxAVG (SFC NOx time-average ppbv)
!@+                    cAVG (SFC CO time-average ppbv)
!@+                    oAVG1,nAVG1 (L=1 Ox and NO2 time-average ppbv)
!@+                    PM2p5, PM10 (SFC time-average PM2.5 and PM10 ppmm)
!@+                    PM2p51,PM101(L=1 time-average PM2.5 and PM10 ppmm)
!@+                    cPM2p5,cPM10 (SFC time-average PM2.5, PM10 kg/m3)
!@+                    NO2col NO2 column amount, instant., (kg/m2)
!@+                    D*          (HDO on any model level)
!@+                    B*          (BE7 on any model level)
!@+                    SO4, RAPR
!@+                    7BEW, 7BED, BE7ATM
!@+                    CTEM,CD3D,CI3D,CL3D,CDN3D,CRE3D,CLWP  ! aerosol
!@+                    TAUSS,TAUMC,CLDSS,CLDMC,MCCTP
!@+                    SO4_d1,SO4_d2,SO4_d3,   ! het. chem
!@+                    Clay, Silt1, Silt2, Silt3  ! dust
!@+                    TrSMIXR surface mixing ratio for all tracers [kg/kg]
!@+                    TrSCONC surface concentration for all tracers [kg/m^3]
!@+                    DuEMIS soil dust aerosol emission flux [kg/m^2/s]
!@+                    DuEMIS2 soil dust aerosol emission flux [kg/m^2/s]
!@+                            from cubed wind speed (only diagnostic)
!@+                    DuDEPTURB turbulent depo of soil dust aer [kg/m^2/s]
!@+                    DuDEPGRAV grav settling of soil dust aerosols [kg/m^2/s]
!@+                    DuDEPWET wet deposition of soil dust aerosols [kg/m^2/s]
!@+                    DuLOAD soil dust aer load of atmospheric column [kg/m^2]
!@+                    DuCONC three-dimensional soil dust concentrations [kg/m^3]
!@+                    DuSMIXR surface mix ratio of soil dust aerosols [kg/kg]
!@+                    DuSCONC surface conc of soil dust aerosols [kg/m^3]
!@+                    DuAOD dust aer opt depth daily avg [1]
!@+                    DuCSAOD clear sky dust aer opt depth daily avg [1]
!@+                    AOD aer opt dep (1,nTracerRadiaActive in rad code) daily avg
!@+                    tAOD aer opt dep (sum 1,nTracerRadiaActive) daily avg
!@+                    ctAOD and cAOD are clr-sky versions of tAOD/AOD
!@+                    ictAOD clr-sky sum AOD, 'instantaneous', 3D
!@+                    itAAOD all-sky sum AOD ext-scat band6, instant., 3D
!@+                    FRAC land fractions over 6 types
!@+                    RNFT total runoff over land surface
!@+
!@+   More options can be added as extra cases in this routine
!@auth Gavin Schmidt/Reto Ruedy
      USE CONSTANT, only : grav,rgas,bygrav,bbyg,gbyrb,sday,tf,mair,sha
     *     ,lhe,rhow,undef,stbo,bysha
      USE MODEL_COM, only : lm,p,ptop,zatmo,u,v,focean,flice,t,q
      USE GEOM, only : imaxj,axyp,byaxyp
      USE PBLCOM, only : tsavg,qsavg,usavg,vsavg
      USE CLOUDS_COM, only : llow,lmid,lhi,cldss,cldmc,taumc,tauss,fss
     *           ,svlat,svlhx
#ifdef CLD_AER_CDNC
     *           ,cdn3d,cre3d,clwp
#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
     *           ,ctem,cd3d,ci3d,cl3d
#endif
      USE DYNAMICS, only : ptropo,am,byam,wsave,pk,phi,pmid
      USE FLUXES, only : prec,dmua,dmva,tflux1,qflux1,uflux1,vflux1
     *     ,gtemp,gtempr
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only : mNO2,sOx_acc,sNOx_acc,sCO_acc
     *     ,l1Ox_acc,l1NO2_acc,save_NO2column
#endif
#if (defined TRACERS_SPECIAL_Shindell) || (defined CALCULATE_LIGHTNING)
      USE LIGHTNING, only : saveC2gLightning,saveLightning
#endif
      USE SEAICE_COM, only : rsi,snowi
      USE LANDICE_COM, only : snowli
      USE LAKES_COM, only : flake
      USE GHY_COM, only : snowe,fearth,wearth,aiearth,soil_surf_moist
      USE RAD_COM, only : trhr,srhr,srdn,salb,cfrac,cosz1
     &     ,tausumw,tausumi,fsrdif,difnir
#ifdef TRACERS_ON
     & ,ttausv_sum,ttausv_sum_cs,ttausv_count,ttausv_save,ttausv_cs_save
     & ,aerAbs6SaveInst
#endif
      USE DIAG_COM, only : z_inst,rh_inst,t_inst,tdiurn,pmb,lname_strlen
     * ,isccp_diags,saveHCLDI,saveMCLDI,saveLCLDI,saveCTPI,saveTAUI
     * ,saveSCLDI,saveTCLDI,saveMCCLDTP
#ifdef TRACERS_SPECIAL_Shindell
     * ,o_inst,n_inst,m_inst,x_inst
#endif
#ifdef TES_LIKE_DIAGS
     * ,t_more,q_more
#ifdef TRACERS_SPECIAL_Shindell
     * ,o_more,n_more,m_more,x_more
#endif
#endif

      IMPLICIT NONE
      REAL*4, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: DATA
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: DATAR8
      INTEGER :: I,J,K,L,kp,ks,kunit,n,n1,nc
      REAL*8 POICE,PEARTH,PLANDI,POCEAN,QSAT,PS,SLP, ZS,TAUL
      INTEGER :: J_0,J_1,J_0S,J_1S,I_0,I_1
      LOGICAL :: polefix,have_south_pole,have_north_pole,skip
      INTEGER :: DAY_OF_MONTH ! for daily averages

!@var qinstant flag whether output is instanteneous or accumulated
!@+            /averaged over nsubdd internal time units
      logical :: qinstant
!@var units_of_data units of data for netcdf output
      character(len=24) :: units_of_data
!@var long_name long name for netcdf output
      character(len=lname_strlen+len(kgz_max_suffixes)) :: long_name

      DAY_OF_MONTH = (1+ITIME-ITIME0)/NDAY

      CALL GET(GRID,J_STRT=J_0, J_STOP=J_1,
     &              J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE=have_south_pole,
     &               HAVE_NORTH_POLE=have_north_pole)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

      datar8 = 0.d0
      data = 0.

      kunit=0
C**** depending on namedd string choose what variables to output
      nameloop: do k=1,kdd

        qinstant = .true.
        units_of_data = 'not yet set in get_subdd'
        long_name = 'not yet set in get_subdd'

C**** simple diags (one record per file)
        select case (namedd(k))
        case ("SLP")            ! sea level pressure (mb)
          do j=J_0,J_1
          do i=I_0,imaxj(j)
            ps=(p(i,j)+ptop)
            zs=bygrav*zatmo(i,j)
            datar8(i,j)=slp(ps,tsavg(i,j),zs)
            units_of_data = '10^2 Pa'
            long_name = 'Sea Level Pressure'
          end do
          end do
        case ("PS")             ! surface pressure (mb)
          datar8=p+ptop
            units_of_data = '10^2 Pa'
            long_name = 'Surface Pressure'
        case ("SAT")            ! surf. air temp (C)
          datar8=tsavg-tf
            units_of_data = 'C'
            long_name = 'Surface Air Temperature'
        case ("US")             ! surf. u wind (m/s)
          datar8=usavg
            units_of_data = 'm/s'
            long_name = 'U Component of Surface Air Velocity'
        case ("VS")             ! surf. v wind (m/s)
          datar8=vsavg
            units_of_data = 'm/s'
            long_name = 'V Component of Surface Air Velocity'
        case ("SST")            ! sea surface temp (C)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (FOCEAN(I,J)+FLAKE(I,J).gt.0) then
                datar8(i,j)=GTEMP(1,1,i,j)
              else
                datar8(i,j)=undef
              end if
            end do
          end do
          units_of_data = 'C'
          long_name = 'Sea Surface Temperature'
        case ("SIT")       ! surface sea/lake ice temp (C)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J)).gt.0) then
                datar8(i,j)=GTEMP(1,2,i,j)
              else
                datar8(i,j)=undef
              end if
            end do
          end do
          units_of_data = 'C'
          long_name = 'Surface Sea/Lake Ice Temperature'
        case ("LIT")       ! surface land ice temp (C)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (FLICE(I,J).gt.0) then
                datar8(i,j)=GTEMP(1,3,i,j)
              else
                datar8(i,j)=undef
              end if
            end do
          end do
          units_of_data = 'C'
          long_name = 'Surface Land Ice Temperature'
        case ("GT1")      ! level 1 ground temp (LAND) (C)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                datar8(i,j)=gtemp(1,4,i,j)
              else
                datar8(i,j)=undef
              end if
            end do
          end do
          units_of_data = 'C'
          long_name = 'Level 1 Ground Temperature, Land'
        case ("GTD")   ! avg levels 2-6 ground temp (LAND) (C)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                datar8(i,j)=gdeep(i,j,1)
              else
                datar8(i,j)=undef
              end if
            end do
          end do
          units_of_data = 'C'
          long_name = 'Average Levels 2-6 Ground Temperature, Land'
        case ("GWD")  ! avg levels 2-6 ground liq water (m)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                datar8(i,j)=gdeep(i,j,2)
              else
                datar8(i,j)=undef
              end if
            end do
          end do
          units_of_data = 'm'
          long_name = 'Average Levels 2-6 Ground Liquid Water, Land'
        case ("GID")  ! avg levels 2-6 ground ice (m liq. equiv.)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                datar8(i,j)=gdeep(i,j,3)
              else
                datar8(i,j)=undef
              end if
            end do
          end do
          units_of_data = 'm liq. equiv.'
          long_name = 'Average Levels 2-6 Ground Ice, Land'
        case ("GW0")  ! ground lev 1 + canopy liq water (m)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                datar8(i,j)=wearth(i,j)
              else
                datar8(i,j)=undef
              end if
            end do
          end do
          units_of_data = 'm'
          long_name = 'Ground Level 1 + Canopy Liquid Water, Land'
        case ("GI0")  ! ground lev 1 + canopy ice (m liq. equiv.)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                datar8(i,j)=aiearth(i,j)
              else
                datar8(i,j)=undef
              end if
            end do
          end do
          units_of_data = 'm liq. equiv.'
          long_name = 'Ground Level 1 + Canopy Ice'
        case ("QS")             ! surf spec humidity (kg/kg)
          datar8=qsavg
          units_of_data = 'kg/kg'
        case ("RS")             ! surf rel humidity
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              datar8(i,j)=qsavg(i,j)/qsat(tsavg(i,j),lhe,p(i,j)+ptop)
     &             *100.d0
            enddo
          enddo
          units_of_data = '%'
          long_name = 'Surface Relative Humidity'
        case ("PREC")           ! precip (mm/day)
c          datar8=sday*prec/dtsrc
          datar8=sday*P_acc/(Nsubdd*dtsrc) ! accum over Nsubdd steps
          P_acc=0.
          units_of_data = 'mm/day'
          long_name = 'Precipitation'
          qinstant = .false.
        case ("RNFT")           ! total runoff (mm/day)
          datar8=sday*R_acc/(Nsubdd*dtsrc) ! accum over Nsubdd steps
          R_acc=0.
          units_of_data = 'mm/day'
          long_name = 'Total runoff'
          qinstant = .false.
#ifdef CALCULATE_FLAMMABILITY
        case ("RAPR")   !running avg precip (mm/day)
          datar8=sday*raP_acc/(Nsubdd*dtsrc) ! accum over Nsubdd steps
          raP_acc=0.
          units_of_data = 'mm/day'
          long_name = 'Running Average of Precipitation'
          qinstant = .false.
#endif
#ifdef TRACERS_SPECIAL_Shindell
        case ("oAVG")   ! Nsubdd-step average SFC Ox tracer (ppbv)
          datar8=sOx_acc/real(Nsubdd) ! accum over Nsubdd steps, already in ppbv
          sOx_acc=0.
          units_of_data = 'ppbv'
          long_name = 'Average Surface Ox Tracer'
          qinstant = .false.
        case ("nxAVG")   ! Nsubdd-step average SFC NOx tracer (ppbv)
          datar8=sNOx_acc/real(Nsubdd) ! accum over Nsubdd steps, already in ppbv
          sNOx_acc=0.
          units_of_data = 'ppbv'
          long_name = 'Average Surface NOx Tracer'
          qinstant = .false.
        case ("cAVG")   ! Nsubdd-step average SFC CO tracer (ppbv)
          datar8=sCO_acc/real(Nsubdd) ! accum over Nsubdd steps, already in ppbv
          sCO_acc=0.
          units_of_data = 'ppbv'
          long_name = 'Average Surface CO Tracer'
          qinstant = .false.
        case ("oAVG1")  ! Nsubdd-step average L=1 Ox tracer (ppbv)
          datar8=l1Ox_acc/real(Nsubdd) ! accum over Nsubdd steps, already in ppbv
          l1Ox_acc=0.
          units_of_data = 'ppbv'
          long_name = 'Average Level 1 Ox Tracer'
          qinstant = .false.
        case ("nAVG1")  ! Nsubdd-step average L=1 NO2 (ppbv)
          datar8=l1NO2_acc/real(Nsubdd) ! accum over Nsubdd steps, already in ppbv
          l1NO2_acc=0.
          units_of_data = 'ppbv'
          long_name = 'Average Level 1 NO2'
          qinstant = .false.
        case ("NO2col") ! instantaneous NO2 column amount (kg/m2)
          datar8=save_NO2column
          units_of_data = 'kg/m^2'
          long_name = 'NO2 Column Amount'
#endif /* TRACERS_SPECIAL_Shindell */

        case ("MCP")       ! moist conv precip (mm/day)
          datar8=sday*PM_acc/(Nsubdd*dtsrc) ! accum over Nsubdd steps
          PM_acc=0.
          units_of_data = 'mm/day'
          long_name = 'Moist Convective Precipitation'
          qinstant = .false.
#ifdef TRACERS_WATER
        case ("TRP1")
          datar8=sday*TRP_acc(1,:,:)/(Nsubdd*dtsrc*axyp(:,:)) ! accum over Nsubdd steps
          TRP_acc(1,:,:)=0.
          units_of_data = 'kg/(s m^2)'
          qinstant = .false.
        case ("TRE1")
          datar8=sday*TRE_acc(1,:,:)/(Nsubdd*dtsrc*axyp(:,:)) ! accum over Nsubdd steps
          TRE_acc(1,:,:)=0.
          units_of_data = 'kg/(s m^2)'
          qinstant = .false.
        case ("TRP2")
          datar8=sday*TRP_acc(2,:,:)/(Nsubdd*dtsrc*axyp(:,:)) ! accum over Nsubdd steps
          TRP_acc(2,:,:)=0.
          units_of_data = 'kg/(s m^2)'
          qinstant = .false.
        case ("TRE2")
          datar8=sday*TRE_acc(2,:,:)/(Nsubdd*dtsrc*axyp(:,:)) ! accum over Nsubdd steps
          TRE_acc(2,:,:)=0.
          units_of_data = 'kg/(s m^2)'
          qinstant = .false.
        case ("TRP3")
          datar8=sday*TRP_acc(3,:,:)/(Nsubdd*dtsrc*axyp(:,:)) ! accum over Nsubdd steps
          TRP_acc(3,:,:)=0.
          units_of_data = 'kg/(s m^2)'
          qinstant = .false.
        case ("TRE3")
          datar8=sday*TRE_acc(3,:,:)/(Nsubdd*dtsrc*axyp(:,:)) ! accum over Nsubdd steps
          TRE_acc(3,:,:)=0.
          units_of_data = 'kg/(s m^2)'
          qinstant = .false.
#endif
        case ("SNOWD")     ! snow depth (w.e. mm)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              POICE=RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J))
              PEARTH=FEARTH(I,J)
              PLANDI=FLICE(I,J)
              datar8(i,j)=1d3*(SNOWI(I,J)*POICE+SNOWLI(I,J)*PLANDI
     &             +SNOWE(I,J)*PEARTH)/RHOW
            end do
          end do
          units_of_data = 'w.e. mm'
          long_name = 'Snow Depth'
        case ("SNOWC")     ! snow cover (fraction of grid)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              datar8(i,j)=0.d0
              POICE=RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J))
              if(SNOWI(I,J) > 0.)datar8(i,j)=datar8(i,j)+POICE
              PEARTH=FEARTH(I,J)
              if(SNOWE(I,J) > 0.)datar8(i,j)=datar8(i,j)+PEARTH
              PLANDI=FLICE(I,J)
              if(SNOWLI(I,J) > 0.)datar8(i,j)=datar8(i,j)+PLANDI
              datar8(i,j)=min(1.0,datar8(i,j))
            end do
          end do
          units_of_data = 'fraction of grid area'
          long_name = 'Snow Cover'
        case ("QLAT")           ! latent heat (W/m^2)
          datar8=qflux1*lhe
          units_of_data = 'W/m^2'
          long_name = 'Latent Heat'
        case ("QSEN")           ! sensible heat flux (W/m^2)
          datar8=tflux1*sha
          units_of_data = 'W/m^2'
          long_name = 'Sensible Heat Flux'
        case ("SWD")            ! solar downward flux at surface (W/m^2)
          datar8=srdn*cosz1       ! multiply by instant cos zenith angle
          units_of_data = 'W/m^2'
          long_name = 'Solar Downward Flux at Surface'
        case ("SWU")            ! solar upward flux at surface (W/m^2)
! estimating this from the downward x albedo, since that's already saved
          datar8=srdn*(1.-salb)*cosz1
          units_of_data = 'W/m^2'
          long_name = 'Solar Upward Flux at Surface'
        case ("LWD")            ! LW downward flux at surface (W/m^2)
          datar8=TRHR(0,:,:)
          units_of_data = 'W/m^2'
          long_name = 'Longwave Downward Flux at Surface'
        case ("LWU")            ! LW upward flux at surface (W/m^2)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              POCEAN=(1.-RSI(I,J))*(FOCEAN(I,J)+FLAKE(I,J))
              POICE=RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J))
              PEARTH=FEARTH(I,J)
              PLANDI=FLICE(I,J)
              datar8(i,j)=STBO*(POCEAN*GTEMPR(1,I,J)**4+
     *             POICE *GTEMPR(2,I,J)**4+PLANDI*GTEMPR(3,I,J)**4+
     *             PEARTH*GTEMPR(4,I,J)**4)
            end do
          end do
          units_of_data = 'W/m^2'
          long_name = 'Longwave Upward Flux at Surface'
        case ("SWDF")           ! SW downward diffuse flux at surface (W/m^2)
          datar8=FSRDIF(:,:)+DIFNIR(:,:)
          units_of_data = 'W/m^2'
          long_name = 'Solar Downward Diffuse Flux at Surface'
        case ("LWT")            ! LW upward flux at TOA (P1) (W/m^2)
          do j=J_0,J_1     ! sum up all cooling rates + net surface emission
            do i=I_0,imaxj(j)
              POCEAN=(1.-RSI(I,J))*(FOCEAN(I,J)+FLAKE(I,J))
              POICE=RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J))
              PEARTH=FEARTH(I,J)
              PLANDI=FLICE(I,J)
              datar8(i,j)=-SUM(TRHR(0:LM,I,J))+
     *             STBO*(POCEAN*GTEMPR(1,I,J)**4+
     *             POICE *GTEMPR(2,I,J)**4+PLANDI*GTEMPR(3,I,J)**4+
     *             PEARTH*GTEMPR(4,I,J)**4)
            end do
          end do
          units_of_data = 'W/m^2'
          long_name = 'Longwave Upward Flux at Top of Atmosphere'
        case ("SWT")            ! SW net flux at TOA (P1) (W/m^2)
          do j=J_0,J_1     ! sum up all heating rates + surface absorption
            do i=I_0,imaxj(j)
              datar8(i,j)=SUM(SRHR(0:LM,I,J))*cosz1(i,j)
            end do
          end do
          units_of_data = 'W/m^2'
          long_name = 'Solar Net Flux at Top of Atmosphere'
        case ("ICEF")           ! ice fraction over open water (%)
          datar8=RSI*100.
          units_of_data = '%'
          long_name = 'Ice Fraction Over Open Water'
        case ("STX")            ! E-W surface stress (N/m^2)
          datar8=uflux1
          units_of_data = 'N/m^2'
          long_name = 'East-West Surface Stress'
        case ("STY")            ! N-S surface stress (N/m^2)
          datar8=vflux1
          units_of_data = 'N/m^2'
          long_name = 'North-South Surface Stress'
        case ("LCLD")           ! low level cloud cover (%)
          datar8=0.               ! Warning: these can be greater >100!
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              do l=1,llow
                datar8(i,j)=datar8(i,j)+(cldss(l,i,j)+cldmc(l,i,j))
              end do
              datar8(i,j)=datar8(i,j)*100.d0/real(llow,kind=8)
            end do
          end do
          units_of_data = '%'
          long_name = 'Low Level Cloud Cover'
        case ("MCLD")           ! mid level cloud cover (%)
          datar8=0.               ! Warning: these can be greater >100!
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              do l=llow+1,lmid
                datar8(i,j)=datar8(i,j)+(cldss(l,i,j)+cldmc(l,i,j))
              end do
              datar8(i,j)=datar8(i,j)*100.d0/real(lmid-llow,kind=8)
            end do
          end do
          units_of_data = '%'
          long_name = 'Mid Level Cloud Cover'
        case ("HCLD")           ! high level cloud cover (%)
          datar8=0.               ! Warning: these can be greater >100!
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              do l=lmid+1,lhi
                datar8(i,j)=datar8(i,j)+(cldss(l,i,j)+cldmc(l,i,j))
              end do
              datar8(i,j)=datar8(i,j)*100.d0/real(lhi-lmid,kind=8)
            end do
          end do
          units_of_data = '%'
          long_name = 'High Level Cloud Cover'
        case ("TCLD")           ! total cloud cover (%) (As seen by rad)
          datar8=cfrac*100.d0
          units_of_data = '%'
          long_name = 'Total Cloud Cover (as seen by rad)'
        case ("PTRO")           ! tropopause pressure (mb)
          datar8 = ptropo
          units_of_data = '10^2 Pa'
          long_name = 'Tropopause Pressure'

! attempting here {L,M,H}CLDI,CTPI,TAUI,TCLDI (ISCCP quantities):
! please use critically:
        case ("HCLDI")      ! HIGH LEVEL CLOUDINESS (ISCCP)
          if (isccp_diags.eq.1) then
            do j=J_0,J_1
              do i=I_0,imaxj(j)
                if(saveSCLDI(i,j)==0.)then
                  datar8(i,j) = undef
                else
                  datar8(i,j) = 100.d0*saveHCLDI(i,j)!/1.
                endif
              end do
            end do
          else
            datar8=undef
          end if
          units_of_data = '%'
          long_name = 'high level cloudiness isccp'
        case ("MCLDI")      ! MID LEVEL CLOUDINESS (ISCCP)
          if (isccp_diags.eq.1) then
            do j=J_0,J_1
              do i=I_0,imaxj(j)
                if(saveSCLDI(i,j)==0.)then
                  datar8(i,j) = undef
                else
                  datar8(i,j) = 100.d0*saveMCLDI(i,j)!/1.
                endif
              end do
            end do
          else
            datar8=undef
          end if
          units_of_data = '%'
          long_name = 'mid level cloudiness isccp'
        case ("LCLDI")      ! LOW LEVEL CLOUDINESS (ISCCP)
          if (isccp_diags.eq.1) then
            do j=J_0,J_1
              do i=I_0,imaxj(j)
                if(saveSCLDI(i,j)==0.)then
                  datar8(i,j) = undef
                else
                  datar8(i,j) = 100.d0*saveLCLDI(i,j)!/1.
                endif
              end do
            end do
          else
            datar8=undef
          end if
          units_of_data = '%'
          long_name = 'low level cloudiness isccp'
        case ("CTPI")      ! CLOUD TOP PRESSURE (ISCCP)
          if (isccp_diags.eq.1) then
            do j=J_0,J_1
              do i=I_0,imaxj(j)
                if(saveTCLDI(i,j)==0.)then
                  datar8(i,j) = undef
                else
                  datar8(i,j) = saveCTPI(i,j)!/1.
                endif
              end do
            end do
          else
            datar8=undef
          end if
          units_of_data = 'mb'
          long_name = 'cloud top pressure isccp'
        case ("TAUI")      ! CLOUD OPTICAL DEPTH (ISCCP)
          if (isccp_diags.eq.1) then
            do j=J_0,J_1
              do i=I_0,imaxj(j)
                if(saveTCLDI(i,j)==0.)then
                  datar8(i,j) = undef
                else
                  datar8(i,j) = saveTAUI(i,j)!/1.
                endif
              end do
            end do
          else
            datar8=undef
          end if
          units_of_data = '1'
          long_name = 'cloud optical depth isccp'
        case ("MCCTP")   ! MOIST CONVECTIVE CLOUD TOP PRESSURE
          datar8=saveMCCLDTP
          units_of_data = 'mb'
          long_name = 'moist convective cloud top pressure'
        case ("TAUSUMW")   ! WATER CLOUD OPTICAL DEPTH
          datar8=tausumw
          units_of_data = '1'
          long_name = 'water cloud optical depth, vertical sum'
        case ("TAUSUMI")   ! ICE CLOUD OPTICAL DEPTH
          datar8=tausumi
          units_of_data = '1'
          long_name = 'ice cloud optical depth, vertical sum'
#if (defined TRACERS_SPECIAL_Shindell) || (defined CALCULATE_LIGHTNING)
        case ("LGTN")  ! lightning flash rate (flash/m2/s)
          datar8 = saveLightning
          units_of_data = 'flash/m^2/s'
          long_name = 'Lightning Flash Rate'
        case ("c2gLGTN")!cloud-to-ground lightning flash rate(flash/m2/s)
          datar8 = saveC2gLightning
          units_of_data = 'flash/m^2/s'
          long_name = 'Cloud to Ground Lightning Flash Rate'
#endif /* TRACERS_SPECIAL_Shindell or CALCULATE_LIGHTNING*/
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
        case ("PM2p5") ! Nsubdd-step avg SFC PM2.5 (ppmm)
           datar8=sPM2p5_acc/real(Nsubdd)
           sPM2p5_acc=0.
          units_of_data = 'ppmm'
          long_name = 'Surface Particulate Matter <= 2.5 um'
          qinstant = .false.
        case ("PM10") ! Nsubdd-step avg SFC PM10 (ppmm)
           datar8=sPM10_acc/real(Nsubdd)
           sPM10_acc=0.
          units_of_data = 'ppmm'
          long_name = 'Surface Particulate Matter <= 10 um'
          qinstant = .false.
        case ("PM2p51") ! Nsubdd-step avg L=1 PM2.5 (ppmm)
           datar8=l1PM2p5_acc/real(Nsubdd)
           l1PM2p5_acc=0.
          units_of_data = 'ppmm'
          long_name = 'Layer 1 Particulate Matter <= 2.5 um'
          qinstant = .false.
        case ("PM101") ! Nsubdd-step avg L=1 PM10 (ppmm)
           datar8=l1PM10_acc/real(Nsubdd)
           l1PM10_acc=0.
          units_of_data = 'ppmm'
          long_name = 'Layer 1 Particulate Matter <= 10 um'
          qinstant = .false.
        case ("cPM2p5") ! Nsubdd-step avg SFC PM2.5 (kg/m3)
           datar8=csPM2p5_acc/real(Nsubdd)
           csPM2p5_acc=0.
          units_of_data = 'kg/m^3'
          long_name = 'Surface Particulate Matter <= 2.5 um'
          qinstant = .false.
        case ("cPM10") ! Nsubdd-step avg SFC PM10 (kg/m3)
           datar8=csPM10_acc/real(Nsubdd)
           csPM10_acc=0.
          units_of_data = 'kg/m^3'
          long_name = 'Surface Particulate Matter <= 10 um'
          qinstant = .false.
#endif /* (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) */

#ifdef TRACERS_AEROSOLS_Koch
        case ("SO4")      ! sulfate in L=1
          datar8=trm(:,:,1,n_SO4)
          units_of_data = 'kg'
          long_name = 'Layer 1 Sulfate Mass'
#ifdef TRACERS_HETCHEM
          datar8 = datar8 + trm(:,:,1,n_SO4_d1) + trm(:,:,1,n_SO4_d2) +
     &         trm(:,:,1,n_SO4_d3)
          long_name =
     &         'Layer 1 Mass of Sulfate + Dust Coated with Sulfate'
#endif
#endif
#ifdef CLD_AER_CDNC
        case ("CLWP")             !LWP (kg m-2)
          datar8=clwp
          units_of_data = 'kg/m^2'
          long_name = 'Cloud Liquid Water Path'
#endif
#ifdef TRACERS_COSMO
        case ("7BEW")
          datar8=Be7w_acc
          Be7w_acc=0.
          units_of_data = 'kg/m^2'
          qinstant = .false.

        case ("7BED")
          datar8=Be7d_acc
          Be7d_acc=0.
          units_of_data = 'kg/m^2'
          qinstant = .false.

        case ("7BES")
          datar8=1.d6*trcsurf(:,:,n_Be7)   ! 10^-6 kg/kg
          units_of_data = '10^-6 kg/kg'
#endif
        case ("SMST") ! near surface soil moisture (kg/m^3)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                datar8(i,j)=soil_surf_moist(i,j)
              else
                datar8(i,j)=undef
              end if
            end do
          end do
          units_of_data = 'kg/m^3'
          long_name = 'Near Surface Soil Moisture'

        case default
          goto 10
        end select
        kunit=kunit+1
        polefix=.true.
        data=datar8
        call write_data(data,kunit,polefix)
#ifdef NEW_IO_SUBDD
        call write_subdd(trim(namedd(k)),datar8,polefix,units_of_data
     &       ,long_name=long_name,qinstant=qinstant)
#endif
        cycle
 10     continue

c**** variables written at end of day only
        select case (namedd(k))
        case ("TMIN","TMAX","RHMIN","RHMAX","WSMAX")
          kunit=kunit+1
          if (mod(itime+1,Nday).ne.0) cycle ! except at end of day
          polefix=.true.
          select case (namedd(k))
          case ("TMIN")         ! min daily temp (C)
            datar8=tdiurn(:,:,9)-tf
            units_of_data = 'C'
            long_name = 'Minimum Daily Surface Temperature'
            qinstant=.false.
          case ("TMAX")         ! max daily temp (C)
            datar8=tdiurn(:,:,6)-tf
            units_of_data = 'C'
            long_name = 'Maximum Daily Surface Temperature'
            qinstant=.false.
          case ("RHMIN")         ! min relative humidity (%)
            datar8=tdiurn(:,:,10)
            units_of_data = '%'
            long_name = 'Minimum Daily Surface Relative Humidity'
            qinstant=.false.
          case ("RHMAX")         ! max relative humidity (%)
            datar8=tdiurn(:,:,11)
            units_of_data = '%'
            long_name = 'Maximum Daily Surface Relative Humidity'
            qinstant=.false.
          case ("WSMAX")         ! max surface wind speed (m/s)
            datar8=tdiurn(:,:,12)
            units_of_data = 'M/S'
            long_name = 'Maximum Daily Surface Wind Speed'
            qinstant=.false.
          end select
          data=datar8
          call write_data(data,kunit,polefix)
#ifdef NEW_IO_SUBDD
          call write_subdd(trim(namedd(k)),datar8,polefix,units_of_data
     &         ,long_name=long_name,record=day_of_month,qinstant=.false.
     &         )
#endif
          cycle
        end select

C**** diags on soil levels
        select case (namedd(k)(1:2))
        case ("GT","GW","GI") ! soil temperature, water, ice
          if (namedd(k)(3:5) .eq. "ALL") then
            kunit=kunit+1
            do ks=1,ngm
              select case (namedd(k)(1:2))
              case ("GT")        ! soil temperature lvls 1-6, land (C)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    if (fearth(i,j).gt.0) then
                      datar8(i,j)=gsaveL(i,j,ks,1)
                    else
                      datar8(i,j)=undef
                    end if
                  end do
                end do
                units_of_data = 'C'
                long_name = 'Soil Temperature Layers 1-6, Land'
              case ("GW")        ! ground wetness lvls 1-6, land (m)
                ! 8/13/10: for RELATIVE wetness, edit giss_LSM/GHY.f
                ! and activate the corresponding lines where wtr_L is set
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    if (fearth(i,j).gt.0) then
                      datar8(i,j)=gsaveL(i,j,ks,2)
                    else
                      datar8(i,j)=undef
                    end if
                  end do
                end do
                units_of_data = 'm'
                long_name = 'Ground Wetness Layers 1-6, Land'
              case ("GI")  ! ground ice lvls 1-6, land (m, liq equiv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    if (fearth(i,j).gt.0) then
                      datar8(i,j)=gsaveL(i,j,ks,3)
                    else
                      datar8(i,j)=undef
                    end if
                  end do
                end do
                units_of_data = 'liq. equiv. m'
                long_name = 'Ground Ice Layers 1-6, Land'
              end select
              polefix=.true.
              ngm_array(:,:,ks)=datar8
              data=datar8
              call write_data(data,kunit,polefix)
            end do
#ifdef NEW_IO_SUBDD
            call write_subdd(trim(namedd(k)),ngm_array,polefix
     &           ,units_of_data,long_name=long_name,positive='down')
#endif
            cycle
          end if
        end select

C**** diags on fixed pressure levels or velocity
        select case (namedd(k)(1:1))
        case ("Z","R","T","Q",  ! heights, rel/spec humidity or temp
     &        "O","X","M","N")  ! Ox, NOx, CO, NO2
C**** get pressure level
          do kp=1,kgz_max
            if (namedd(k)(2:5) .eq. PMNAME(kp)) then
              kunit=kunit+1
              select case (namedd(k)(1:1))
              case ("Z")        ! geopotential heights
                datar8=z_inst(kp,:,:)
                units_of_data = 'm'
                long_name = 'Geopotential Height at '//trim(PMNAME(kp))
     &               //' hPa'
              case ("R")        ! relative humidity (wrt water)
                datar8=rh_inst(kp,:,:)*100.d0
                units_of_data = '%'
                long_name = 'Relative Humidity at '// trim(PMNAME(kp))//
     &               ' hPa'
              case ("Q")        ! specific humidity
                do j=J_0,J_1
                do i=I_0,imaxj(j)
                  datar8(i,j)=rh_inst(kp,i,j)*qsat(t_inst(kp,i,j)+tf,lhe
     *                 ,PMB(kp))
                end do
                end do
                units_of_data = 'kg/kg'
                long_name = 'Specific Humidity at ' //trim(PMNAME(kp))//
     &               ' hPa'
              case ("T")        ! temperature (C)
                datar8=t_inst(kp,:,:)
                units_of_data = 'C'
                long_name = 'Temperature at '//trim(PMNAME(kp))//' hPa'
#ifdef TRACERS_SPECIAL_Shindell
              case ("O")        ! Ox  tracer (ppmv)
                datar8=O_inst(kp,:,:)
                units_of_data = 'ppmv'
                long_name = 'Ox tracer at '//trim(PMNAME(kp))//' hPa'
              case ("X")        ! NOx  tracer (ppmv)
                datar8=X_inst(kp,:,:)
                units_of_data = 'ppmv'
                long_name = 'NOx tracer at '//trim(PMNAME(kp))//' hPa'
              case ("M")        ! CO  tracer (ppmv)
                datar8=M_inst(kp,:,:)
                units_of_data = 'ppmv'
                long_name = 'CO tracer at '//trim(PMNAME(kp))//' hPa'
              case ("N")        ! NO2 non-tracer (ppmv)
                datar8=N_inst(kp,:,:)
                units_of_data = 'ppmv'
                long_name = 'NO2 at '//trim(PMNAME(kp))//' hPa'
#endif /* TRACERS_SPECIAL_Shindell */
              end select
              polefix=.true.
              data=datar8
              call write_data(data,kunit,polefix)
#ifdef NEW_IO_SUBDD
              call write_subdd(trim(namedd(k)),datar8,polefix
     &             ,units_of_data,long_name=long_name)
#endif
              cycle nameloop
            end if
          end do


#ifdef TES_LIKE_DIAGS
          if (namedd(k)(2:4) .eq. "ALL") then
            kunit=kunit+1
            do kp=1,kgz_max_more
              select case (namedd(k)(1:1))
              case ("Q")        ! specific humidity
                datar8=q_more(kp,:,:)
                units_of_data = 'kg/kg'
                long_name = 'Specific Humidity'
              case ("T")        ! temperature (C)
                datar8=t_more(kp,:,:)
                units_of_data = 'C'
                long_name = 'Temperature'
#ifdef TRACERS_SPECIAL_Shindell
              case ("O")        ! Ox  tracer (ppmv)
                datar8=o_more(kp,:,:)
                units_of_data = 'ppmv'
                long_name = 'Ox Tracer'
              case ("X")        ! NOx tracer (ppmv)
                datar8=x_more(kp,:,:)
                units_of_data = 'ppmv'
                long_name = 'NOx Tracer'
              case ("M")        ! CO  tracer (ppmv)
                datar8=m_more(kp,:,:)
                units_of_data = 'ppmv'
                long_name = 'CO Tracer'
              case ("N")        ! NO2 non-tracer (ppmv)
                datar8=n_more(kp,:,:)
                units_of_data = 'ppmv'
                long_name = 'NO2 (not a tracer)'
#endif /* TRACERS_SPECIAL_Shindell */
              end select
              polefix=.true.
              kgz_max_more_suffixes(kp)=trim(PMNAMEmore(kp))//'_hPa'
              kgz_max_more_array(:,:,kp)=datar8
              data=datar8
              call write_data(data,kunit,polefix)
            end do
#ifdef NEW_IO_SUBDD
            call write_subdd(trim(namedd(k)),kgz_max_more_array,polefix
     &           ,units_of_data,long_name=long_name,suffixes
     &           =kgz_max_more_suffixes,positive='down')
#endif
            cycle
          end if
#endif /* TES_LIKE_DIAGS */

          if (namedd(k)(2:4) .eq. "ALL") then
            do kp=1,kgz_max
              kunit=kunit+1
              select case (namedd(k)(1:1))
              case ("Z")        ! geopotential heights
                datar8=z_inst(kp,:,:)
                units_of_data = 'm'
                long_name = 'Geopotential Height'
              case ("R")        ! relative humidity (wrt water)
                datar8=rh_inst(kp,:,:)*100.d0
                units_of_data = '%'
                long_name = 'Relative Humidity'

#ifndef TES_LIKE_DIAGS /* note NOT defined */
              case ("Q")        ! specific humidity
                do j=J_0,J_1
                do i=I_0,imaxj(j)
                  datar8(i,j)=rh_inst(kp,i,j)*qsat(t_inst(kp,i,j)+tf,lhe
     *                 ,PMB(kp))
                end do
                end do
                units_of_data = 'kg/kg'
                long_name = 'Specific Humidity'
              case ("T")        ! temperature (C)
                datar8=t_inst(kp,:,:)
                units_of_data = 'C'
                long_name = 'Temperature'
#ifdef TRACERS_SPECIAL_Shindell
              case ("O")        ! Ox  tracer (ppmv)
                datar8=o_inst(kp,:,:) 
                units_of_data = 'ppmv'
                long_name = 'Ox Tracer'
              case ("X")        ! NOx tracer (ppmv)
                datar8=x_inst(kp,:,:) 
                units_of_data = 'ppmv'
                long_name = 'NOx Tracer'
              case ("M")        ! CO  tracer (ppmv)
                datar8=m_inst(kp,:,:) 
                units_of_data = 'ppmv'
                long_name = 'CO Tracer'
              case ("N")        ! NO2 non-tracer (ppmv)
                datar8=n_inst(kp,:,:) 
                units_of_data = 'ppmv'
                long_name = 'NO2 (not a tracer)'
#endif /* TRACERS_SPECIAL_Shindell */
#endif /* NOT defined TES_LIKE_DIAGS */
              end select
              polefix=.true.
              kgz_max_suffixes(kp)=trim(PMNAME(kp))//'_hPa'
              kgz_max_array(:,:,kp)=datar8
              data=datar8
              call write_data(data,kunit,polefix)
            end do
#ifdef NEW_IO_SUBDD
            call write_subdd(trim(namedd(k)),kgz_max_array,polefix
     &           ,units_of_data,long_name=long_name,suffixes
     &           =kgz_max_suffixes,positive='down')
#endif
            cycle
          end if

C**** diagnostics on model levels
        case ("U","V","W","C","o","B","D","x","t","q","z","r","m","n")
             ! velocity/clouds/tracers, temp,spec.hum.,geo.ht
          if (namedd(k)(2:4) .eq. "ALL") then
            kunit=kunit+1
            do kp=1,LmaxSUBDD
              skip = .false.
              select case (namedd(k)(1:1))
              case ("t")        ! temperature (C)
                if(have_south_pole) datar8(1:im,1)=
     &          t(1,1,kp)*pk(kp,1, 1)-tf
                if(have_north_pole) datar8(1:im,jm)=
     &          t(1,jm,kp)*pk(kp,1,jm)-tf
                datar8(:,J_0S:J_1S)=
     &          t(:,J_0S:J_1S,kp)*pk(kp,:,J_0S:J_1S)-tf
                units_of_data = 'C'
                long_name = 'Temperature'
              case ("r")        ! relative humidity
                if(have_south_pole) datar8(1:im, 1)=q(1,1,kp)/
     &          qsat(t(1,1,kp)*pk(kp,1,1),lhe,pmid(kp,1,1))
                if(have_north_pole) datar8(1:im,jm)=q(1,jm,kp)/
     &          qsat(t(1,jm,kp)*pk(kp,1,jm),lhe,pmid(kp,1,jm))
                do j=J_0S,J_1S; do i=i_0,i_1
                  datar8(i,j)=q(i,j,kp)/qsat(t(i,j,kp)*pk(kp,i,j),
     &            lhe,pmid(kp,i,j))*100.d0
                enddo         ; enddo
                units_of_data = '%'
                long_name = 'Relative Humidity'
              case ("q")        ! specific humidity
                if(have_south_pole) datar8(1:im, 1)=q(1, 1,kp)
                if(have_north_pole) datar8(1:im,jm)=q(1,jm,kp)
                datar8(:,J_0S:J_1S)=q(:,J_0S:J_1S,kp)
                units_of_data = 'kg/kg'
                long_name = 'Specific Humidity'
              case ("z")        ! geopotential height
                if(have_south_pole) datar8(1:im, 1)=phi(1, 1,kp)
                if(have_north_pole) datar8(1:im,jm)=phi(1,jm,kp)
                datar8(:,J_0S:J_1S)=phi(:,J_0S:J_1S,kp)
                units_of_data = 'm'
              case ("U")        ! E-W velocity
                datar8=u(:,:,kp)
                units_of_data = 'm/s'
                long_name = 'U-Velocity'
              case ("V")        ! N-S velocity
                datar8=v(:,:,kp)
                units_of_data = 'm/s'
                long_name = 'V-Velocity'
              case ("W")        ! vertical velocity
                if(kp<lm) then
                  datar8=wsave(:,:,kp)
                units_of_data = 'Pa/s'
                long_name = 'Vertical Velocity'
                else
                  datar8=undef
                  skip = .true.
                end if
              case ("C")        ! estimate of cloud optical depth
                datar8=(1.-fss(kp,:,:))*taumc(kp,:,:)+fss(kp,:,:)
     *               *tauss(kp,:,:)
                units_of_data = '1'
                long_name = 'Estimate of Cloud Optical Depth'
#ifdef TRACERS_SPECIAL_Shindell
              case ("o")                ! Ox ozone tracer (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1.d6*trm(i,j,kp,n_Ox)*mass2vol(n_Ox)/
     *                   (am(kp,i,j)*axyp(i,j))
                  end do
                end do
                units_of_data = 'ppmv'
                long_name = 'Ox Ozone Tracer'
              case ("x")                ! NOx tracer (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1.d6*trm(i,j,kp,n_NOx)*mass2vol(n_NOx)/
     *                   (am(kp,i,j)*axyp(i,j))
                  end do
                end do
                units_of_data = 'ppmv'
                long_name = 'NOx Tracer'
              case ("m")                ! CO tracer (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1.d6*trm(i,j,kp,n_CO)*mass2vol(n_CO)/
     *                   (am(kp,i,j)*axyp(i,j))
                  end do
                end do
                units_of_data = 'ppmv'
                long_name = 'CO Tracer'
              case ("n")                ! NO2 (not a tracer) (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1.d6*mNO2(i,j,kp)
                  end do
                end do
                units_of_data = 'ppmv'
                long_name = 'NO2 (not a tracer)'
#endif
#ifdef TRACERS_COSMO
              case ("B")                ! Be7 tracer
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1.d6*trm(i,j,kp,n_Be7)* mass2vol(n_Be7)/
     *                   (am(kp,i,j)*axyp(i,j))
                  end do
                end do
                units_of_data = 'ppmv'
                long_name = 'Be7 Tracer'
#endif
#ifdef TRACERS_SPECIAL_O18
              case ("D")                ! HDO tracer (permil)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1d3*(trm(i,j,kp,n_HDO)/
     *                   (trm(i,j,kp,n_water)*trw0(n_HDO))-1.)
                  end do
                end do
                units_of_data = 'permil'
                long_name = 'HDO Tracer'
#endif
              end select
              polefix=(namedd(k)(1:1).ne."U".and.namedd(k)(1:1).ne."V")
              LmaxSUBDD_array(:,:,kp)=datar8
              data=datar8
              if(.not.skip) call write_data(data,kunit,polefix)
            end do
#ifdef NEW_IO_SUBDD
            call write_subdd(trim(namedd(k)),LmaxSUBDD_array,polefix
     &           ,units_of_data,long_name=long_name,positive='up')
#endif
            cycle
          end if
C**** get model level
          do l=1,lm
            if (trim(namedd(k)(2:5)) .eq. lst(l)) then
              kunit=kunit+1
              select case (namedd(k)(1:1))
              case ("t")        ! temperature (C)
                if(have_south_pole) datar8(1:im,1)=
     &               t(1,1,l)*pk(l,1, 1)-tf
                if(have_north_pole) datar8(1:im,jm)=
     &               t(1,jm,l)*pk(l,1,jm)-tf
                datar8(:,J_0S:J_1S)=
     &               t(:,J_0S:J_1S,l)*pk(l,:,J_0S:J_1S)-tf
                units_of_data = 'C'
                long_name = 'Temperature at Level '//trim(lst(l))
              case ("r")        ! relative humidity
                if(have_south_pole) datar8(1:im, 1)=q(1, 1,l)/
     &               qsat(t(1,1,l)*pk(l,1,1),lhe,pmid(l,1,1))
                if(have_north_pole) datar8(1:im,jm)=q(1,jm,l)/
     &               qsat(t(1,jm,l)*pk(l,1,jm),lhe,pmid(l,1,jm))
                do j=J_0S,J_1S; do i=i_0,i_1
                  datar8(i,j)=q(i,j,l)/qsat(t(i,j,l)*pk(l,i,j),
     &                 lhe,pmid(l,i,j))*100.d0
                enddo
              enddo
                units_of_data = '%'
                long_name = 'Relative Humidity at Level '//trim(lst(l))
              case ("q")        ! specific humidity
                if(have_south_pole) datar8(1:im, 1)=q(1, 1,l)
                if(have_north_pole) datar8(1:im,jm)=q(1,jm,l)
                datar8(:,J_0S:J_1S)=q(:,J_0S:J_1S,l)
                units_of_data = 'kg/kg'
                long_name = 'Specific Humidity at Level '//trim(lst(l))
              case ("z")        ! geopotential height
                if(have_south_pole) datar8(1:im, 1)=phi(1, 1,l)
                if(have_north_pole) datar8(1:im,jm)=phi(1,jm,l)
                datar8(:,J_0S:J_1S)=phi(:,J_0S:J_1S,l)
                units_of_data = 'm'
                long_name = 'Geopotential Height at Level '//trim(lst(l)
     &               )
              case ("U")        ! U velocity
                datar8=u(:,:,l)
                units_of_data = 'm/s'
                long_name = 'U-Velocity at Level '//trim(lst(l))
              case ("V")        ! V velocity
                datar8=v(:,:,l)
                units_of_data = 'm/s'
                long_name = 'V-Velocity at Level '//trim(lst(l))
              case ("W")        ! W velocity
                datar8=wsave(:,:,l)
                units_of_data = 'Pa/s'
                long_name = 'Vertical Velocity at Level '//trim(lst(l))
              case ("C")        ! estimate of cloud optical depth
                datar8=(1.-fss(l,:,:))*taumc(l,:,:)+fss(l,:,:)
     *               *tauss(l,:,:)
                units_of_data = '1'
                long_name = 'Estimate of Cloud Optical Depth at Level '
     &               //trim(lst(l))
#ifdef TRACERS_SPECIAL_Shindell
              case ("o")                ! Ox ozone tracer (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1.d6*trm(i,j,l,n_Ox)*mass2vol(n_Ox)/
     *                   (am(l,i,j)*axyp(i,j))
                  end do
                end do
                units_of_data = 'ppmv'
                long_name = 'Ox Ozone Tracer at Level '//trim(lst(l))
              case ("x")                ! NOx tracer (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1.d6*trm(i,j,l,n_NOx)*mass2vol(n_NOx)/
     *                   (am(l,i,j)*axyp(i,j))
                  end do
                end do
                units_of_data = 'ppmv'
                long_name = 'NOx Tracer at Level '//trim(lst(l))
              case ("m")                ! CO tracer (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1.d6*trm(i,j,l,n_CO)*mass2vol(n_CO)/
     *                   (am(l,i,j)*axyp(i,j))
                  end do
                end do
                units_of_data = 'ppmv'
                long_name = 'CO Tracer at Level '//trim(lst(l))
              case ("n")                ! NO2 (not a tracer) (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1.d6*mNO2(i,j,l)
                  end do
                end do
                units_of_data = 'ppmv'
                long_name = 'NO2 (not a tracer) at Level '//trim(lst(l))
#endif
#ifdef TRACERS_COSMO
              case ("B")                ! Be7 tracer
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1.d6*trm(i,j,l,n_Be7)* mass2vol(n_Be7)/
     *                   (am(l,i,j)*axyp(i,j))
                  end do
                end do
                units_of_data = 'ppmv'
                long_name = 'Be7 Tracer at Level '//trim(lst(l))
#endif
#ifdef TRACERS_SPECIAL_O18
              case ("D")                ! HDO tracer (permil)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    datar8(i,j)=1d3*(trm(i,j,l,n_HDO)/(trm(i,j,l,n_water
     *                   )*trw0(n_HDO))-1.)
                  end do
                end do
                units_of_data = 'permil'
                long_name = 'HDO Tracer at Level '//trim(lst(l))
#endif
              end select
              polefix=(namedd(k)(1:1).ne."U".and.namedd(k)(1:1).ne."V")
              data=datar8
              call write_data(data,kunit,polefix)
#ifdef NEW_IO_SUBDD
              call write_subdd(trim(namedd(k)),datar8,polefix
     &             ,units_of_data,long_name=long_name)
#endif
              cycle nameloop
            end if
          end do
        end select

C**** Additional diags - multiple records per file
        select case (namedd(k))

C**** cases using all levels up to LmaxSUBDD
          case ("SO2", "SO4", "SO4_d1", "SO4_d2", "SO4_d3", "Clay",
     *         "Silt1", "Silt2", "Silt3", "CTEM", "CL3D", "CI3D", "CD3D"
     *         , "CLDSS", "CLDMC", "CDN3D", "CRE3D", "TAUSS", "TAUMC",
     *         "RADHEAT","CLWP","itAOD","ictAOD","itAAOD")
          kunit=kunit+1
          do l=1,LmaxSUBDD
            select case(namedd(k))

#ifdef TRACERS_ON
c***** 3D i(c)tAOD instantaneous sum over tracers of aerosol opt depth
c***** (keep in mind that depending on nRAD and NSUBDD, this could be
c***** "instantaneous" is a relative term.)
            case ("itAOD","ictAOD")   !tot aero(+dust,etc) opt dep, inst.
              if (any(tracerRadiaActiveFlag)) then
                datar8=0.
                do n=1,ntm        ! sum over rad code tracers is used
                  if(tracerRadiaActiveFlag(n))then
                    select case(namedd(k))
                    case('itAOD')
                      datar8=datar8+ttausv_save(:,:,n,L)
                      units_of_data='1'
                      long_name = 'Total All Sky Aerosol Optical Depth'
                    case('ictAOD')
                      datar8=datar8+ttausv_cs_save(:,:,n,L)
                      units_of_data='1'
                      long_name =
     &                     'Total Clear Sky Aerosol Optical Depth'
                    end select
                  end if
                end do
              else
                write(6,*) 'Warning: No radiatively active tracers'
                write(6,*) ' ',trim(namedd(k)),' not written'
              end if
c***** 3D itAAOD instantaneous sum over tracers of aerosol opt depth
c***** Band 6 extinction-scatter (so absorption)
c***** (keep in mind that depending on nRAD and NSUBDD, this could be
c***** "instantaneous" is a relative term.)
            case ("itAAOD")   !tot abs aero opt dep, inst.
              datar8=aerAbs6SaveInst(:,:,L)
              units_of_data='1'
              long_name = 'Total All Sky Aerosol Optical Depth'
#endif /*TRACERS_ON*/


#ifdef TRACERS_HETCHEM
            case ("SO2")
              datar8=trm(:,:,l,n_SO2)
              units_of_data = 'kg'
              long_name = 'Sulfur Dioxide Mass'
            case ("SO4")
              datar8=trm(:,:,l,n_SO4)
              units_of_data = 'kg'
              long_name = 'Sulfate Mass'
            case ("SO4_d1")
              datar8=trm(:,:,l,n_SO4_d1)
              units_of_data = 'kg'
              long_name = 'Mass of Sulfate Coated With 0.1 to 1 um Dust'
            case ("SO4_d2")
              datar8= trm(:,:,l,n_SO4_d2)
              units_of_data = 'kg'
              long_name = 'Mass of Sulfate Coated With 1 to 2 um Dust'
            case ("SO4_d3")
              datar8=trm(:,:,l,n_SO4_d3)
              units_of_data = 'kg'
              long_name = 'Mass of Sulfate Coated With 2 to 4 um Dust'
            case ("Clay")
              datar8=trm(:,:,l,n_Clay)
              units_of_data = 'kg'
              long_name = 'Mass of Sulfate Coated With 4 to 8 um Dust'
            case ("Silt1")
              datar8=trm(:,:,l,n_Silt1)
              units_of_data = 'kg'
              long_name = 'Mass of Silt 1 to 2 um'
            case ("Silt2")
              datar8=trm(:,:,l,n_Silt2)
              units_of_data = 'kg'
              long_name = 'Mass of Silt 2 to 4 um'
            case ("Silt3")
              datar8=trm(:,:,l,n_Silt3)
              units_of_data = 'kg'
              long_name = 'Mass of Silt 4 to 8 um'
#endif
            case ("CLDSS")
              datar8=100.d0*cldss(l,:,:) ! Cld cover LS(%)
              units_of_data = '%'
              long_name = 'Large Scale Cloud Cover'
            case ("CLDMC")
              datar8=100.d0*cldmc(l,:,:) ! Cld cover MC(%)
              units_of_data = '%'
              long_name = 'Moist Convective Cloud Cover'
            case ("TAUSS")
              datar8=tauss(l,:,:) ! LS cld tau
              units_of_data = '1'
              long_name = 'Large Scale Cloud Optical Depth'
            case ("TAUMC")
              datar8=taumc(l,:,:) ! MC cld tau
              units_of_data = '1'
              long_name = 'Moist Convective Cloud Optical Depth'
            case ("RADHEAT")
              datar8=(SRHR(L,:,:)*COSZ1(:,:)+TRHR(L,:,:))*
     &             SDAY*bysha*byam(l,:,:)
              units_of_data = 'K/day'
              long_name = 'Radiative Heating Rate'
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
            case ("CTEM")
              datar8=ctem(l,:,:) ! cld temp (K) at cld top
              units_of_data = 'K'
              long_name = 'Cloud Temperature at Cloud Top'
            case ("CL3D")
              datar8=cl3d(l,:,:) ! cld LWC (kg m-3)
              units_of_data = 'kg/m^3'
              long_name = 'Liquid Cloud Water'
            case ("CI3D")
              datar8=ci3d(l,:,:) ! cld IWC (kg m-3)
              units_of_data = 'kg/m^3'
            case ("CD3D")
              datar8=cd3d(l,:,:) ! cld thickness (m)
              units_of_data = 'm'
              long_name = 'Cloud Thickness'
#endif
#ifdef CLD_AER_CDNC
            case ("CDN3D")
              datar8=cdn3d(l,:,:) ! cld CDNC (cm^-3)
              units_of_data = 'cm^3'
            case ("CRE3D")
              datar8=1.d-6*cre3d(l,:,:) ! cld Reff (m)
              units_of_data = 'm'
#endif
            end select
            polefix=.true.
            LmaxSUBDD_array(:,:,l)=datar8
            data=datar8
            call write_data(data,kunit,polefix)
          end do
#ifdef NEW_IO_SUBDD
          call write_subdd(trim(namedd(k)),LmaxSUBDD_array,polefix
     &         ,units_of_data,long_name=long_name,positive='up')
#endif
          cycle

#ifdef TRACERS_ON
c***** for (c)tAOD the sum over tracers of aerosol optical depth
          case ("tAOD","ctAOD")   !tot aero(+dust,etc) opt dep, daily avg
            kunit=kunit+1
            if (any(tracerRadiaActiveFlag)) then
              polefix=.true.
              if(mod(itime+1,Nday).ne.0) cycle ! except at end of day
              if(ttausv_count==0.)call stop_model('ttausv_count=0',255)
              datar8=0.
              do n=1,ntm        ! sum over rad code tracers is used
                if(tracerRadiaActiveFlag(n))then
                  select case(namedd(k))
                  case('tAOD')
                    datar8=datar8+ttausv_sum(:,:,n)
                    units_of_data='1'
                    long_name = 'Total All Sky Aerosol Optical Depth'
                  case('ctAOD')
                    datar8=datar8+ttausv_sum_cs(:,:,n)
                    units_of_data='1'
                    long_name = 'Total Clear Sky Aerosol Optical Depth'
                  end select
                end if
              end do
              datar8=datar8/ttausv_count
              data=datar8
              call write_data(data,kunit,polefix)
#ifdef NEW_IO_SUBDD
              call write_subdd(trim(namedd(k)),datar8,polefix
     &             ,units_of_data,long_name=long_name,record
     &             =day_of_month,qinstant=.false.)
#endif
            else
              write(6,*) 'Warning: No radiatively active tracers'
              write(6,*) ' ',trim(namedd(k)),' not written'
            end if
            cycle

C**** for (c)AOD multiple tracers are written to one file:
          case ('AOD','cAOD')!aerosol opt dep daily avg (all/clear sky)
            kunit=kunit+1
            if (any(tracerRadiaActiveFlag)) then
              polefix=.true.
              if(mod(itime+1,Nday).ne.0) cycle ! except at end of day
              if(ttausv_count==0.)call stop_model('ttausv_count=0',255)
              nc=0
              do n=1,ntm
                if(tracerRadiaActiveFlag(n))then
                  nc=nc+1
                  select case(namedd(k))
                  case('AOD')
                    datar8=ttausv_sum(:,:,n)/ttausv_count
                    units_of_data='1'
                    long_name = 'All Sky Aerosol Optical Depth of'
                  case('cAOD')
                    datar8=ttausv_sum_cs(:,:,n)/ttausv_count
                    units_of_data='1'
                    long_name = 'Clear Sky Aerosol Optical Depth of'
                  end select
                  rTrname(nc)=trname(n)
                  rTRACER_array(:,:,nc)=datar8
                  data=datar8
                  call write_data(data,kunit,polefix)
                end if
              end do
#ifdef NEW_IO_SUBDD
              call write_subdd(trim(namedd(k)),rTRACER_array,polefix
     &             ,units_of_data,long_name=long_name,record
     &             =day_of_month,suffixes=rTrname ,qinstant=.false.)
#endif
            else
              write(6,*) 'Warning: No radiatively active tracers'
              write(6,*) ' ',trim(namedd(k)),' not written'
            end if
            cycle

c**** Mixing ratio for all tracers at surface [kg/kg]
          case('TrSMIXR')
            kunit=kunit+1
            polefix=.true.
            do n=1,ntm
!$OMP PARALLEL DO PRIVATE(i,j)
              do j=j_0,j_1
                do i=i_0,i_1
                  trcSurfMixR_acc(i,j,n)=trcSurfMixR_acc(i,j,n)
     &                 /real(Nsubdd,kind=8)
                  datar8(i,j)=trcSurfMixR_acc(i,j,n)
                  trcSurfMixR_acc(i,j,n)=0.D0
                end do
              end do
!$OMP END PARALLEL DO
              units_of_data='kg/kg'
              long_name = 'Mixing Ratio at Surface of'
              TRACER_array(:,:,n)=datar8
              data=datar8
              call write_data(data,kunit,polefix)
            end do
#ifdef NEW_IO_SUBDD
            call write_subdd(trim(namedd(k)),TRACER_array,polefix
     &           ,units_of_data,long_name=long_name,suffixes=trname
     &           ,qinstant=.false.)
#endif
            cycle

c**** Concentration for all tracers at surface [kg/m^3]
          case('TrSCONC')
            kunit=kunit+1
            polefix=.true.
            do n=1,ntm
!$OMP PARALLEL DO PRIVATE(i,j)
              do j=j_0,j_1
                do i=i_0,i_1
                  trcSurfByVol_acc(i,j,n)=trcSurfByVol_acc(i,j,n)
     &                 /real(Nsubdd,kind=8)
                  datar8(i,j)=trcSurfByVol_acc(i,j,n)
                  trcSurfByVol_acc(i,j,n)=0.D0
                end do
              end do
!$OMP END PARALLEL DO
              units_of_data='kg/m^3'
              long_name = 'Concentration at Surface of'
              TRACER_array(:,:,n)=datar8
              data=datar8
              call write_data(data,kunit,polefix)
            end do
#ifdef NEW_IO_SUBDD
            call write_subdd(trim(namedd(k)),TRACER_array,polefix
     &           ,units_of_data,long_name=long_name,suffixes=trname
     &           ,qinstant=.false.)
#endif
            cycle
#endif /*TRACERS_ON*/

C**** Write land,lake,ocean,and ice fractions to one file
C**** Possbily will remove at some point, since kinda
C**** overkill, but useful now for EPA down-scaling:
          case ('FRAC') ! land, ocean, lake fractions incl ice+ice-free
            kunit=kunit+1
            polefix=.true.
            data=RSI(:,:)*FOCEAN(:,:)         ! ocean ice
            call write_data(data,kunit,polefix)
            data=RSI(:,:)*FLAKE(:,:)          ! lake ice
            call write_data(data,kunit,polefix)
            data=(1.-RSI(:,:))*FOCEAN(:,:)    ! ice-free ocean
            call write_data(data,kunit,polefix)
            data=(1.-RSI(:,:))*FLAKE(:,:)     ! ice-free lake
            call write_data(data,kunit,polefix)
            data=FLICE(:,:)                   ! land ice
            call write_data(data,kunit,polefix)
            data=FEARTH(:,:)                  ! ice-free land
            call write_data(data,kunit,polefix)
            cycle

C**** cases where multiple records go to one file for dust

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
        case ('DuEMIS','DuEMIS2','DuSMIXR','DuSCONC','DuLOAD')
          kunit=kunit+1
          do n=1,Ntm_dust
            n1=n_soilDust+n-1
C**** first set: no 'if' tests
            select case (namedd(k))

            CASE ('DuEMIS')     ! Dust emission flux [kg/m^2/s]
!$OMP PARALLEL DO PRIVATE(i,j)
              do j=j_0,j_1
                do i=i_0,i_1
                  datar8(i,j)=dustDiagSubdd_acc%dustEmission(i,j,n)
     &                 /real(Nsubdd,kind=8)
                  dustDiagSubdd_acc%dustEmission(i,j,n)=0.D0
                end do
              end do
!$OMP END PARALLEL DO
              units_of_data='kg/(s*m^2)'
              long_name = 'Emission of'
#ifdef TRACERS_DUST
            CASE ('DuEMIS2')    ! Dust emission flux 2 (diag. var. only) [kg/m^2/s]
!$OMP PARALLEL DO PRIVATE(i,j)
              do j=j_0,j_1
                do i=i_0,i_1
                  datar8(i,j)=dustDiagSubdd_acc%dustEmission2(i,j,n)
     &                 /real(Nsubdd,kind=8)
                  dustDiagSubdd_acc%dustEmission2(i,j,n)=0.D0
                end do
              end do
!$OMP END PARALLEL DO
              units_of_data='kg/(s*m^2)'
              long_name =
     &             'Emission According to Cubic Formula (diag only) of'
#endif
            CASE ('DuSMIXR')      ! Mixing ratio of dust tracers at surface [kg/kg]
!$OMP PARALLEL DO PRIVATE(i,j)
              do j=j_0,j_1
                do i=i_0,i_1
                  datar8(i,j)=dustDiagSubdd_acc%dustSurfMixR(i,j,n)
     &                 /real(Nsubdd,kind=8)
                  dustDiagSubdd_acc%dustSurfMixR(i,j,n)=0.D0
                end do
              end do
!$OMP END PARALLEL DO
              units_of_data='kg/kg'
              long_name = 'Mixing Ratio at Surface of'
            CASE ('DuSCONC')  ! Concentration of dust tracers at surface [kg/m^3]
!$OMP PARALLEL DO PRIVATE(i,j)
              do j=j_0,j_1
                do i=i_0,i_1
                  datar8(i,j)=dustDiagSubdd_acc%dustSurfConc(i,j,n)
     &                 /real(Nsubdd,kind=8)
                  dustDiagSubdd_acc%dustSurfConc(i,j,n)=0.D0
                end do
              end do
!$OMP END PARALLEL DO
              units_of_data='kg/m^3'
              long_name = 'Concentration at Surface of'
            CASE ('DuLOAD')     ! Dust load [kg/m^2]
!$OMP PARALLEL DO PRIVATE(i,j,l)
              do l=1,lm
                do j=j_0,j_1
                  do i=i_0,i_1
                    dustDiagSubdd_acc%dustMass(i,j,l,n)
     &                   =dustDiagSubdd_acc%dustMass(i,j,l,n)
     &                   /real(Nsubdd,kind=8)
                  end do
                end do
              end do
!$OMP END PARALLEL DO
              datar8=0.D0
!$OMP PARALLEL DO PRIVATE(i,j,l)
              do j=j_0,j_1
                do i=i_0,i_1
                  do l=1,lm
                    datar8(i,j)=datar8(i,j)+dustDiagSubdd_acc%dustMass(i
     &                   ,j,l,n)
                  end do
                  datar8(i,j)=datar8(i,j)*byaxyp(i,j)
                  dustDiagSubdd_acc%dustMass(i,j,:,n)=0.D0
                end do
              end do
!$OMP END PARALLEL DO
              units_of_data='kg/m^2'
              long_name = 'Mass in Atmospheric Column of'
            end select
            polefix=.true.
            dust3d_array(:,:,n)=datar8
            data=datar8
            call write_data(data,kunit,polefix)
          end do
#ifdef NEW_IO_SUBDD
          call write_subdd(trim(namedd(k)),dust3d_array,polefix
     &         ,units_of_data,long_name=long_name,suffixes=dust_names
     &         ,qinstant=.false.)
#endif
          cycle

        case('DuCONC')      ! Dust concentration [kg/m^3]
          kunit=kunit+1
          polefix=.true.
          do n=1,Ntm_dust
            do l=1,LmaxSUBDD
!$OMP PARALLEL DO PRIVATE(i,j)
              do j=j_0,j_1
                do i=i_0,i_1
                  datar8(i,j) = dustDiagSubdd_acc%dustConc(i,j,l,n)
     &                 *byaxyp(i,j)/real(Nsubdd,kind=8)
                  dustDiagSubdd_acc%dustConc(i,j,l,n)=0.d0
                end do
              end do
!$OMP END PARALLEL DO
              dust4d_array(:,:,l,n) = datar8
              data=datar8
              call write_data(data,kunit,polefix)
            end do
            units_of_data='kg/m^3'
            long_name = 'Concentration of'
          end do
#ifdef NEW_IO_SUBDD
          call write_subdd(trim(namedd(k)),dust4d_array,polefix
     &         ,units_of_data,long_name=long_name,suffixes=dust_names
     &         ,qinstant=.false.,positive='up')
#endif
          cycle

C**** other dust special cases

          case ("DuAOD","DuCSAOD") !tot dust aero opt dep, daily avg
            kunit=kunit+1
            if(mod(itime+1,Nday).ne.0) cycle ! except at end of day
            if(ttausv_count==0.)call stop_model('ttausv_count=0',255)
            datar8=0.
            do n=1,Ntm_dust
              n1=n_soilDust+n-1
              select case(namedd(k))
              case('DuAOD')
                datar8(:,:)=ttausv_sum(:,:,n1)
                units_of_data='1'
                long_name = 'All Sky Optical Depth of'
              case('DuCSAOD')
                datar8(:,:)=ttausv_sum_cs(:,:,n1)
                units_of_data='1'
                long_name = 'Clear Sky Optical Depth of'
              end select
              datar8=datar8/ttausv_count
              polefix=.true.
              dust3d_array(:,:,n) = datar8
              data=datar8
              call write_data(data,kunit,polefix)
            end do
#ifdef NEW_IO_SUBDD
            call write_subdd(trim(namedd(k)),dust3d_array,polefix
     &           ,units_of_data,long_name=long_name,record=day_of_month
     &           ,suffixes=dust_names ,qinstant=.false.)
#endif
            cycle

#ifdef TRACERS_DRYDEP
          case ('DuDEPTURB')        ! Turb. deposition flux of dust tracers [kg/m^2/s]
          kunit=kunit+1
          do n=1,Ntm_dust
            n1=n_soilDust+n-1
            if (dodrydep(n1)) then
!$OMP PARALLEL DO PRIVATE(i,j)
              do j=j_0,j_1
                do i=i_0,i_1
                  datar8(i,j)=dustDiagSubdd_acc%dustDepoTurb(i,j,n)
     &                 /real(Nsubdd,kind=8)
                  dustDiagSubdd_acc%dustDepoTurb(i,j,n)=0.D0
                end do
              end do
!$OMP END PARALLEL DO
              polefix=.true.
              dust3d_array(:,:,n)=datar8
              data=datar8
              call write_data(data,kunit,polefix)
            end if
            units_of_data='kg/(s*m^2)'
            long_name = 'Turbulent Deposition of'
          end do
#ifdef NEW_IO_SUBDD
          call write_subdd(trim(namedd(k)),dust3d_array,polefix
     &         ,units_of_data,long_name=long_name,suffixes=dust_names
     &         ,qinstant=.false.)
#endif
          cycle

          case ('DuDEPGRAV')      ! Gravit. settling flux of dust tracers [kg/m^2/s]
          kunit=kunit+1
          do n=1,Ntm_dust
            n1=n_soilDust+n-1
            IF (dodrydep(n1)) THEN
!$OMP PARALLEL DO PRIVATE(i,j)
              do j=j_0,j_1
                do i=i_0,i_1
                  datar8(i,j)=dustDiagSubdd_acc%dustDepoGrav(i,j,n)
     &                 /real(Nsubdd,kind=8)
                  dustDiagSubdd_acc%dustDepoGrav(i,j,n)=0.D0
                end do
              end do
!$OMP END PARALLEL DO
              polefix=.true.
              dust3d_array(:,:,n)=datar8
              data=datar8
              call write_data(data,kunit,polefix)
            end if
            units_of_data='kg/(s*m^2)'
            long_name = 'Gravitational Settling of'
          end do
#ifdef NEW_IO_SUBDD
          call write_subdd(trim(namedd(k)),dust3d_array,polefix
     &         ,units_of_data,long_name=long_name,suffixes=dust_names
     &         ,qinstant=.false.)
#endif
          cycle
#endif /*TRACERS_DRYDEP*/

          case ('DuDEPWET')         ! Wet deposition flux of dust tracers [kg/m^2/s]
          kunit=kunit+1
          do n=1,Ntm_dust
            n1=n_soilDust+n-1
#ifdef TRACERS_WATER
            if (dowetdep(n1)) then
#endif
!$OMP PARALLEL DO PRIVATE(i,j)
              do j=j_0,j_1
                do i=i_0,i_1
                  datar8(i,j)=dustDiagSubdd_acc%dustMassInPrec(i,j,n)
     &                 *byaxyp(i,j)/Dtsrc/real(Nsubdd,kind=8)
                  dustDiagSubdd_acc%dustMassInPrec(i,j,n)=0.D0
                end do
              end do
!$OMP END PARALLEL DO
#ifdef TRACERS_WATER
            end if
#endif
            polefix=.true.
            dust3d_array(:,:,n)=datar8
            data=datar8
            call write_data(data,kunit,polefix)
            units_of_data='kg/(s*m^2)'
            long_name = 'Wet Deposition of'
          end do
#ifdef NEW_IO_SUBDD
          call write_subdd(trim(namedd(k)),dust3d_array,polefix
     &         ,units_of_data,long_name=long_name,suffixes=dust_names
     &         ,qinstant=.false.)
#endif
          cycle
#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM*/

C**** this prevents tokens that are not caught from messing up the file data
        case default
          kunit=kunit+1
        end select

      end do nameloop ! end of  k=1,kdd loop
c****
      return
      end subroutine get_subdd

      subroutine write_data(data,kunit,polefix)
!@sum write out subdd data array with optional pole fix
      use domain_decomp_1d, only : grid,get,writei_parallel,
     &     hasSouthPole, hasNorthPole

      implicit none

      real*4, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: data
      integer kunit
      logical :: polefix

#ifdef NEW_IO_SUBDD
      return
#endif

c**** fix polar values
      if (polefix) then
        if(hassouthpole(grid)) data(2:im,1) =data(1,1)
        if(hasnorthpole(grid)) data(2:im,jm)=data(1,jm)
      end if
      call writei_parallel(grid,iu_subdd(kunit),
     *     nameunit(iu_subdd(kunit)),data,itime+1)

      end subroutine write_data

#ifdef NEW_IO_SUBDD

c def_global_attr_subdd
      subroutine def_global_attr_subdd(fid,q24,qinst)
!@sum def_global_attr_subdd defines global attributes in subdd output files
!@auth Jan Perlwitz

      use pario, only: write_attr
      use domain_decomp_atm, only: grid

      implicit none

      integer,intent(in) :: fid
      logical,intent(in) :: q24,qinst

      character(len=3) :: cnsubdd
      character(len=37) :: method

      call write_attr(grid,fid,'global','xlabel',xlabel)
      if (qinst) then
        method = 'instantaneous value for current ITU'
      else
        if (q24) then
          method = 'daily'
        else
          method = 'accumulated/averaged over nsubdd ITUs'
        end if
      end if
      call write_attr(grid,fid,'global','calculation_method',method)
      write(cnsubdd,'(i3)') nsubdd
      call write_attr(grid,fid,'global','nsubdd',cnsubdd)

      return
      end subroutine def_global_attr_subdd

c def_xy_coord_subdd
      subroutine def_xy_coord_subdd(fid)
!@sum def_xy_coord_subdd defines x,y-coordinates in subdd output files
!@auth Jan Perlwitz

#ifdef CUBED_SPHERE
      use geom, only: lon2d_dg,lat2d_dg,lonbds,latbds
#else
      use geom, only: lat_dg,lon_dg
#endif
      use geom, only: axyp
      use pario, only: defvar,write_attr
      use domain_decomp_atm, only: grid

      implicit none

      integer,intent(in) :: fid

      logical :: qr4
      character(len=8) :: dist_x,dist_y

#ifdef CUBED_SPHERE
      dist_x='dist_im'
      dist_y='dist_jm'
      qr4=.false.
      call defvar(grid,fid,lon2d_dg,'lon('//trim(dist_x)//','//
     &     trim(dist_y)//')',r4_on_disk=qr4)
      call write_attr(grid,fid,'lon','units','degrees_east')
      call write_attr(grid,fid,'lon','long_name','Longitude')
      call write_attr(grid,fid,'lon','bounds','lonbds')
      call defvar(grid,fid,lat2d_dg,'lat('//trim(dist_x)//','//
     &     trim(dist_y)//')',r4_on_disk=qr4)
      call write_attr(grid,fid,'lat','units','degrees_north')
      call write_attr(grid,fid,'lat','long_name','Latitude')
      call write_attr(grid,fid,'lat','bounds','latbds')
      call defvar(grid,fid,lonbds,'lonbds(four,'//trim(dist_x)//','//
     &     trim(dist_y)//')',r4_on_disk=qr4)
      call write_attr(grid,fid,'lonbds','units','degrees_east')
      call write_attr(grid,fid,'lonbds','long_name'
     &     ,'Longitude Boundaries')
      call defvar(grid,fid,latbds,'latbds(four,'//trim(dist_x)//','//
     &     trim(dist_y)//')',r4_on_disk=qr4)
      call write_attr(grid,fid,'latbds','units','degrees_north')
      call write_attr(grid,fid,'latbds','long_name'
     &     ,'Latitude Boundaries')
#else
      dist_x='dist_lon'
      dist_y='dist_lat'
      qr4=.true.
      call defvar(grid,fid,lon_dg(:,1),'lon(lon)',r4_on_disk=qr4)
      call write_attr(grid,fid,'lon','units','degrees_east')
      call write_attr(grid,fid,'lon','long_name','Longitude')
      call defvar(grid,fid,lat_dg(:,1),'lat(lat)',r4_on_disk=qr4)
      call write_attr(grid,fid,'lat','units','degrees_north')
      call write_attr(grid,fid,'lat','long_name','Latitude')
#endif
      call defvar(grid,fid,axyp,'axyp('//trim(dist_x)//','//trim(dist_y)
     &     //')',r4_on_disk=qr4)
      call write_attr(grid,fid,'axyp','units','m^2')
      call write_attr(grid,fid,'axyp','long_name','Gridcell Area')

      return
      end subroutine def_xy_coord_subdd

c time_subdd
      real(kind=8) function time_subdd(q24,rec)
!@sum time_subdd calculates value of time coordinate 'time' for subdd output
!@auth Jan Perlwitz

      implicit none

      integer,intent(in) :: rec
      logical,intent(in) :: q24

      if (q24) then
        time_subdd = real((jyear - iyear1)*jdpery + jday - 1,kind=8)
      else
        time_subdd = real((jyear - iyear1)*jdpery + JDendOfM(jmon - 1)
     &       ,kind=8)*24. + (rec - 1)*nsubdd*dtsrc/3600.
      end if

      return
      end function time_subdd

c get_calendarstring
      subroutine get_calendarstring(qinstant,q24,year,mon,day,hour,itu
     &     ,calendarstring)
!@sum get_calendarstring provides calendarstring from year, month, day,
!@+                      hour, and internal time unit depending on flags
!@auth Jan Perlwitz

      implicit none

!@var q24 flag is true for output every 24 hours
!@var qinstant flag is true for instantaneous variables
      logical,intent(in) :: q24,qinstant
!@var year,mon,day,hour,itu year,months,day,internal time unit
      integer,intent(in) :: year,mon,day,hour,itu
!@var calendar string with date and time
!@+            format: 'year-mm-dd' or 'year-mm-dd hr:mm'
      character(len=*),intent(out) :: calendarstring

!@var modelDate string variable for model date
      character(len=10) :: modelDate
!@var hour_of_day hour of day
      integer :: hour_of_day
!@var minute_of_hour minute of hour
      integer :: minute_of_hour
!@var time_of_day string variable for time of day
      character(len=5) :: time_of_day

      real(kind=8) :: hrs

      modelDate = '0000-00-00'
      write(modelDate(1:4),'(i4)') year
      if (mon <= 9) then
        write(modelDate(7:7),'(i1)') mon
      else
        write(modelDate(6:7),'(i2)') mon
      end if
      if (day <= 9) then
        write(modelDate(10:10),'(i1)') day
      else
        write(modelDate(9:10),'(i2)') day
      end if
      if (q24 .and. .not. qinstant) then
        time_of_day = ''
      else
        time_of_day = '00:00'
      end if
      if (qinstant) then
        hour_of_day = hour
        minute_of_hour = int(mod((real(itu,kind=8) + 0.5)*dtsrc/3600.,
     &       1.)*60.)
      else if (.not. q24) then
        hrs = (real(itu,kind=8) + 1. - nsubdd/2.)*dtsrc/3600.
        hour_of_day = int(mod(hrs,24.))
        minute_of_hour = int(mod(hrs,1.)*60.)
      end if
      if (qinstant .or. (.not. qinstant .and. .not. q24)) then
        if (hour_of_day <= 9) then
          write(time_of_day(2:2),'(i1)') hour_of_day
        else
          write(time_of_day(1:2),'(i2)') hour_of_day
        end if
        if (minute_of_hour <= 9) then
          write(time_of_day(5:5),'(i1)') minute_of_hour
        else
          write(time_of_day(4:5),'(i2)') minute_of_hour
        end if
      end if
      write(calendarstring,'(a10,a1,a5)') modelDate,' ',time_of_day

      return
      end subroutine get_calendarstring

c get_referencetime_for_netcdf
      subroutine get_referencetime_for_netcdf(qinstant,q24,itu
     &     ,referencetime)
!@sum get_referencetime_for_netcdf provides string with reference time of time
!@+                                coordinate for high frequency netcdf-output
!@auth Jan Perlwitz

      implicit none

!@var q24 flag is true for output every 24 hours
!@var qinstant flag is true for instantaneous variables
      logical,intent(in) :: q24,qinstant
!@var itu internal time unit
      integer,intent(in) :: itu
!@var referencetime string with reference time of time coordinate for netcdf
      character(len=*),intent(out) :: referencetime

!@var calendarstring calendar string variable with date and time
      character(len=50) :: calendarstring

      integer :: year1,mon1,day1,jdate1,hour1
      character(len=4) :: amon1

      call getdte(itu,nday,iyear1,year1,mon1,day1,jdate1,hour1,amon1)
      call get_calendarstring(qinstant,q24,year1,mon1,jdate1,hour1,itu
     &     ,calendarstring)
      if (q24) then
        referencetime = 'days since '//trim(calendarstring)
      else
        referencetime = 'hours since '//trim(calendarstring)//' UTC'
      end if

      return
      end subroutine get_referencetime_for_netcdf

c def_time_coord_subdd
      subroutine def_time_coord_subdd(fid,q24,qinst,time,calendarstring)
!@sum def_time_coord_subdd defines time coordinates in subdd output files
!@auth Jan Perlwitz

      use pario, only: defvar,write_attr,set_record_dimname
      use domain_decomp_atm, only: grid

      implicit none

      integer,intent(in) :: fid
      real(kind=8),intent(in) :: time
      logical,intent(in) :: q24,qinst
      character(len=*),intent(in) :: calendarstring

      integer :: length_of_calendarstring
      character(len=23) :: relation_to_itime
      character(len=50) :: itimeUnits,referencetime

      call set_record_dimname(grid,fid,'time')
      call get_referencetime_for_netcdf(qinst,q24,nsubdd-1
     &     ,referencetime)
      call defvar(grid,fid,itime,'itime',with_record_dim=.true.)
      write(itimeUnits,'(a11,i4,a16)') 'ITUs since ',iyear1
     &     ,'-01-01 00:00 UTC'
      call write_attr(grid,fid,'itime','units',itimeUnits)
      call write_attr(grid,fid,'itime','long_name'
     &     ,'Internal Time Unit (ITU)')
      call defvar(grid,fid,time,'time',with_record_dim=.true.
     &     ,r4_on_disk=.true.)
      call write_attr(grid,fid,'time','calendar','noleap')
      call write_attr(grid,fid,'time','units',referencetime)
      call write_attr(grid,fid,'time','long_name'
     &     ,'Time Coordinate Relative to IYEAR1/01/01')
      if (qinst) then
        relation_to_itime = 'midpoint of last ITU'
      else if (.not. q24) then
        relation_to_itime = 'midpoint of nsubdd ITUs'
      else
        relation_to_itime = 'all ITUs of day'
      end if
      call write_attr(grid,fid,'time','relation_to_itime',
     &     relation_to_itime)
      call defvar(grid,fid,calendarstring,
     &     'calendar(length_of_calendarstring)',with_record_dim=.true.)
      call write_attr(grid,fid,'calendar','units','UTC')
      call write_attr(grid,fid,'calendar','long_name'
     &     ,'Date/Time of Midpoint')
      call write_attr(grid,fid,'calendar','relation_to_itime',
     &     relation_to_itime)

      return
      end subroutine def_time_coord_subdd

c write_xy_coord_subdd
      subroutine write_xy_coord_subdd(fid)
!@sum write_xy_coord_subdd writes x,y-coordinates to subdd output files
!@auth Jan Perlwitz

#ifdef CUBED_SPHERE
      use geom, only: lon2d_dg,lat2d_dg,lonbds,latbds
#else
      use geom, only: lat_dg,lon_dg
#endif
      use geom, only: axyp
      use pario, only: write_data,write_dist_data
      use domain_decomp_atm, only: grid

      implicit none

      integer,intent(in) :: fid

#ifdef CUBED_SPHERE
      call write_dist_data(grid,fid,'lon',lon2d_dg)
      call write_dist_data(grid,fid,'lat',lat2d_dg)
      call write_dist_data(grid,fid,'lonbds',lonbds,jdim=3)
      call write_dist_data(grid,fid,'latbds',latbds,jdim=3)
#else
      call write_data(grid,fid,'lon',lon_dg(:,1))
      call write_data(grid,fid,'lat',lat_dg(:,1))
#endif
      call write_dist_data(grid,fid,'axyp',axyp)

      return
      end subroutine write_xy_coord_subdd

c write_time_coord_subdd
      subroutine write_time_coord_subdd(fid,rec,time,qinst,q24
     &     ,calendarstring)
!@sum write_time_coord_subdd writes time coordinates to subdd output files
!@auth Jan Perlwitz

      use model_com, only: itime,jdate,jhour,jmon,jyear
      use pario, only: write_data
      use domain_decomp_atm, only: grid

      implicit none

      integer,intent(in) :: fid,rec
      real(kind=8),intent(in) :: time
      logical,intent(in) :: q24,qinst

      character(len=*),intent(out) :: calendarstring

      call write_data(grid,fid,'itime',itime+1,record=rec)
      call write_data(grid,fid,'time',time,record=rec)
      call get_calendarstring(qinst,q24,jyear,jmon,jdate,jhour,itime
     &     ,calendarstring)
      call write_data(grid,fid,'calendar',calendarstring,record=rec)

      return
      end subroutine write_time_coord_subdd

c write_2d
      subroutine write_2d(qtyname,data,polefix,units_of_data,long_name
     &     ,record,positive,qinstant)
!@sum write_2d high frequency netcdf-output of two-dimensional arrays
!@auth M. Kelley and Jan Perlwitz

      use pario, only : par_open,par_close,par_enddef,defvar,
     &     write_data,write_dist_data,write_attr
      use domain_decomp_atm, only : grid,hasNorthPole,hasSouthPole

      implicit none

      character(len=*),intent(in) :: qtyname
      real*8,dimension(grid%i_strt_halo:,grid%j_strt_halo:) :: data
      logical,intent(in) :: polefix
      character(len=*),intent(in) :: units_of_data
      character(len=*),intent(in),optional :: long_name
      integer, intent(in), optional :: record
      character(len=*),intent(in),optional :: positive
      logical,intent(in),optional :: qinstant

      integer :: fid,rec
      real(kind=8) :: time
      logical :: q24,qinst,qr4
      character(len=80) :: qname
      character(len=80) :: fname
      character(len=9) :: cform
      character(len=16) :: calendarstring
      character(len=8) :: dist_x,dist_y

      if(.not. in_subdd_list(qtyname)) return

      qinst = .true.
      if (present(qinstant)) qinst = qinstant

      fname = trim(qtyname)//aDATE_sv(1:7)//'.'//xlabel(1:lrunid)//'.nc'
      if(present(record)) then
        rec = record
      else
        rec = (1+itime-itime0)/nsubdd
      endif
      if (present(record) .or. mod(nsubdd*dtsrc/3600.,24.) == 0.) then
        q24 = .true.
      else
        q24 = .false.
      end if
      time = time_subdd(q24,rec)

      if(rec==1) then ! define this output file
        fid = par_open(grid,trim(fname),'create')
        call def_global_attr_subdd(fid,q24,qinst)
        call def_xy_coord_subdd(fid)
        call def_time_coord_subdd(fid,q24,qinst,time,calendarstring)
#ifdef CUBED_SPHERE
        dist_x='dist_im'
        dist_y='dist_jm'
        qr4=.true.
#else
        dist_x='dist_lon'
        dist_y='dist_lat'
        qr4=.true.
#endif
        call defvar(grid,fid,data,trim(qtyname)//'('//trim(dist_x)//','
     &       //trim(dist_y)//')',with_record_dim=.true.,r4_on_disk=qr4)
        call write_attr(grid,fid,trim(qtyname),'units',units_of_data)
        if (present (long_name)) then
          if (long_name /= 'not yet set in get_subdd') call
     &         write_attr(grid,fid,trim(qtyname),'long_name',long_name)
        end if

        call par_enddef(grid,fid)
        call write_xy_coord_subdd(fid)
      else
        fid = par_open(grid,trim(fname),'write')
      endif

      call write_time_coord_subdd(fid,rec,time,qinst,q24,calendarstring)

      if(polefix) then
        if(hasSouthPole(grid)) data(2:im,1) = data(1,1)
        if(hasNorthPole(grid)) data(2:im,jm) = data(1,jm)
      endif
      call write_dist_data(grid,fid, trim(qtyname), data, record=rec)
      call par_close(grid,fid)

      return
      end subroutine write_2d

c write_3d
      subroutine write_3d(qtyname,data,polefix,units_of_data,long_name
     &     ,record,suffixes,positive,qinstant)
!@sum write_3d high frequency netcdf-output of three-dimensional arrays
!@auth M. Kelley and Jan Perlwitz

      use pario, only : par_open,par_close,par_enddef,defvar,
     &     write_data,write_dist_data,write_attr
      use domain_decomp_atm, only : grid,hasNorthPole,hasSouthPole

      implicit none

      character(len=*),intent(in) :: qtyname
      real*8,dimension(grid%i_strt_halo:,grid%j_strt_halo:,:) :: data
      logical,intent(in) :: polefix
      character(len=*),intent(in) :: units_of_data
      character(len=*),intent(in),optional :: long_name
      integer, intent(in), optional :: record
      character(len=*), intent(in), optional :: suffixes(:)
      character(len=*),intent(in),optional :: positive
      logical,intent(in),optional :: qinstant

      integer :: fid,l,rec,level,ivertical,l5
      real(kind=8) :: time
      real(kind=8),dimension(size(data,dim=3)) :: vertical
      logical :: q24,qinst,qr4
      character(len=80) :: qname
      character(len=80) :: fname
      character(len=9) :: cform
      character(len=50) :: lname
      character(len=16) :: calendarstring
      character(len=8) :: dist_x,dist_y
      character(len=5) :: c5

      if(.not. in_subdd_list(qtyname)) return

      if (present(suffixes)) then
        if (size(suffixes) /= size(data,3)) call stop_model
     &       ('write_3d: bad sizes',255)
      end if

      qinst = .true.
      if (present(qinstant)) qinst=qinstant

      fname = trim(qtyname)//aDATE_sv(1:7)//'.'//xlabel(1:lrunid)//'.nc'
      if(present(record)) then
        rec = record
      else
        rec = (1+itime-itime0)/nsubdd
      endif
      if (present(record) .or. mod(nsubdd*dtsrc/3600.,24.) == 0.) then
        q24 = .true.
      else
        q24 = .false.
      end if
      time = time_subdd(q24,rec)

      if(rec==1) then ! define this output file
        fid = par_open(grid,trim(fname),'create')
        call def_global_attr_subdd(fid,q24,qinst)
        call def_xy_coord_subdd(fid)

c define vertical coordinates
        if (present(positive)) then ! define vertical coordinates
          call defvar(grid,fid,vertical,'level(level)',r4_on_disk=
     &         .true.)
          if (.not. (positive == 'up' .or. positive == 'down')) call
     &         stop_model ('In subdd: Wrong netcdf positive-attribute'
     &         ,255)
          if (.not. present(suffixes)) then ! define sigma levels or soil layers
            if (positive == 'down') then
              call write_attr(grid,fid,'level','units','layer')
              call write_attr(grid,fid,'level','long_name'
     &             ,'Ground Layer')
           else
              call write_attr(grid,fid,'level','units','sigma_level')
              call write_attr(grid,fid,'level','long_name'
     &             ,'Sigma Level')
            end if
            do l=1,size(data,dim=3)
              vertical(l) = real(l,kind=8)
            end do
          else if (qtyname(2:4) == 'ALL') then ! define vertical pressure coordinates
            call write_attr(grid,fid,'level','units','10^2 Pa')
            call write_attr(grid,fid,'level','long_name'
     &           ,'Pressure Level')
            do l=1,size(data,3)
              if (index(suffixes(l)(1:4),'.') == 0) then
                read(suffixes(l)(1:min(4,index(suffixes(l),'_')-1))
     &               ,'(i4)') ivertical
                vertical(l) = ivertical
              else
                read(suffixes(l)(1:min(4,index(suffixes(l)(1:4),'_')-1))
     &               ,'(f4.2)') vertical(l)
              end if
            end do
          end if
          call write_attr(grid,fid,'level','positive',trim(positive))
        end if

        call def_time_coord_subdd(fid,q24,qinst,time,calendarstring)

c define physical variable
#ifdef CUBED_SPHERE
        dist_x='dist_im'
        dist_y='dist_jm'
        qr4=.true.
#else
        dist_x='dist_lon'
        dist_y='dist_lat'
        qr4=.true.
#endif
        if ((.not. present(positive)) .or. (present(suffixes) .and.
     &       qtyname(2:4) /= 'ALL')) then
                                ! define 2-dim fields with different names
          do l=1,size(data,3)
            if (qtyname(2:4) == 'ALL') then
              qname = qtyname(1:1)//'_'//trim(suffixes(l))
            else
              qname = trim(qtyname)//'_'//trim(suffixes(l))
            end if
            call defvar(grid,fid,data(:,:,l),trim(qname)//'('//
     &           trim(dist_x)//','//trim(dist_y)//')',with_record_dim=
     &           .true.,r4_on_disk=qr4)
            call write_attr(grid,fid,trim(qname),'units',units_of_data)
            if (present(long_name)) then
              if (long_name /= 'not yet set in get_subdd') then
                lname = trim(long_name)//' '//trim(suffixes(l))
                call write_attr(grid,fid,trim(qname),'long_name',lname)
              end if
            end if
          end do
        else ! define 3-dimensional fields
          l5 = min(len_trim(qtyname),5)
          c5 = ''
          c5(1:l5) = qtyname(1:l5)
          if (c5(2:4) == 'ALL' .or. c5(3:5) == 'ALL') then
            qname = qtyname(1:index(qtyname,'ALL')-1)
          else
            qname = qtyname
          end if
          call defvar(grid,fid,data,trim(qname)//'('//trim(dist_x)//','
     &         //trim(dist_y)//',level)',with_record_dim=.true.
     &         ,r4_on_disk=qr4)
          call write_attr(grid,fid,trim(qname),'units',units_of_data)
          if (present(long_name)) then
            if (long_name /= 'not yet set in get_subdd') call
     &           write_attr(grid,fid,trim(qname),'long_name',long_name)
          end if
        endif

        call par_enddef(grid,fid)
        call write_xy_coord_subdd(fid)
        if (present(positive)) call write_data(grid,fid,'level',vertical
     &       )
      else
        fid = par_open(grid,trim(fname),'write')
      endif

      call write_time_coord_subdd(fid,rec,time,qinst,q24,calendarstring)

c write physical variable
      if(polefix) then
        do l=1,size(data,3)
          if(hasSouthPole(grid)) data(2:im,1,l) = data(1,1,l)
          if(hasNorthPole(grid)) data(2:im,jm,l) = data(1,jm,l)
        enddo
      endif
      if ((.not. present(positive)) .or. (present(suffixes) .and.
     &     qtyname(2:4) /= 'ALL')) then
        do l=1,size(data,3)
          if (qtyname(2:4) == 'ALL') then
            qname = qtyname(1:1)//'_'//trim(suffixes(l))
          else
            qname = trim(qtyname)//'_'//trim(suffixes(l))
          end if
          call write_dist_data(grid,fid,trim(qname),data(:,:,l),record
     &         =rec)
        enddo
      else
        l5 = min(len_trim(qtyname),5)
        c5 = ''
        c5(1:l5) = qtyname(1:l5)
        if (c5(2:4) == 'ALL' .or. c5(3:5) == 'ALL') then
          qname = qtyname(1:index(qtyname,'ALL')-1)
        else
          qname = qtyname
        end if
        call write_dist_data(grid,fid,trim(qname),data,record=rec)
      endif
      call par_close(grid,fid)

      return
      end subroutine write_3d

c write_4d
      subroutine write_4d(qtyname,data,polefix,units_of_data,long_name
     &     ,record,suffixes,positive,qinstant)
!@sum write_4d high frequency netcdf-output of four-dimensional arrays
!@auth M. Kelley and Jan Perlwitz

      use pario, only : par_open,par_close,par_enddef,defvar,
     &     write_data,write_dist_data,write_attr
      use domain_decomp_atm, only : grid,hasNorthPole,hasSouthPole

      implicit none

      character(len=*),intent(in) :: qtyname
      real*8, dimension(grid%i_strt_halo:,grid%j_strt_halo:,:,:) :: data
      logical,intent(in) :: polefix
      character(len=*),intent(in) :: units_of_data
      character(len=*),intent(in),optional :: long_name
      integer, intent(in), optional :: record
      character(len=*), intent(in), optional :: suffixes(:)
      character(len=*),intent(in),optional :: positive
      logical,intent(in),optional :: qinstant

      integer :: fid,l,n,rec,level
      real(kind=8),dimension(size(data,dim=3)) :: vertical
      real(kind=8) :: time
      logical :: q24,qinst,qr4
      character(len=80) :: qname
      character(len=80) :: fname
      character(len=9) :: cform
      character(len=50) :: lname
      character(len=16) :: calendarstring
      character(len=8) :: dist_x,dist_y

      if(.not. in_subdd_list(qtyname)) return

      if (present(suffixes)) then
        if (size(suffixes) /= size(data,4)) call stop_model
     &       ('write_4d: bad sizes',255)
      end if

      qinst = .true.
      if (present(qinstant)) qinst = qinstant

      fname = trim(qtyname)//aDATE_sv(1:7)//'.'//xlabel(1:lrunid)//'.nc'
      if(present(record)) then
        rec = record
      else
        rec = (1+itime-itime0)/nsubdd
      endif
      if (present(record) .or. mod(nsubdd*dtsrc/3600.,24.) == 0.) then
        q24 = .true.
      else
        q24 = .false.
      end if
      time = time_subdd(q24,rec)

      if(rec==1) then ! define this output file
        fid = par_open(grid,trim(fname),'create')
        call def_global_attr_subdd(fid,q24,qinst)
        call def_xy_coord_subdd(fid)
        do l=1,size(data,dim=3)
          vertical(l) = real(l,kind=8)
        end do
        call defvar(grid,fid,vertical,'level(level)',r4_on_disk=.true.)

c define vertical coordinates
        if (present(positive)) then
          if (.not. (positive == 'up' .or. positive == 'down')) then
            call stop_model ('In subdd: Wrong netcdf positive-attribute'
     &           ,255)
          else
            if (positive == 'down') then
              call write_attr(grid,fid,'level','units','layer')
              call write_attr(grid,fid,'level','long_name'
     &             ,'Ground Layer')
            else
              call write_attr(grid,fid,'level','units','sigma_level')
              call write_attr(grid,fid,'level','long_name'
     &             ,'Sigma Level')
            end if
            call write_attr(grid,fid,'level','positive',trim(positive))
          end if
        end if

        call def_time_coord_subdd(fid,q24,qinst,time,calendarstring)

c define physical variable
#ifdef CUBED_SPHERE
        dist_x='dist_im'
        dist_y='dist_jm'
        qr4=.true.
#else
        dist_x='dist_lon'
        dist_y='dist_lat'
        qr4=.true.
#endif
        if (present(suffixes)) then
          do n=1,size(data,4)
            qname = trim(qtyname)//'_'//trim(suffixes(n))
            call defvar(grid,fid,data(:,:,:,n),trim(qname)//'('//
     &           trim(dist_x)//','//trim(dist_y)//',level)'
     &           ,with_record_dim=.true.,r4_on_disk=qr4)
            call write_attr(grid,fid,trim(qname),'units',units_of_data)
            if (present(long_name)) then
              if (long_name /= 'not yet set in get_subdd') then
                lname = trim(long_name)//' '//trim(suffixes(n))
                call write_attr(grid,fid,trim(qname),'long_name',lname)
              end if
            end if
          end do
        else
          call defvar(grid,fid,data,trim(qtyname)//'
     &         ('// trim(dist_x)//','//trim(dist_y)//',level,ntm)'
     &         ,with_record_dim=.true.,r4_on_disk=qr4)
          call write_attr(grid,fid,trim(qtyname),'units',units_of_data)
          if (present(long_name)) then
            if (long_name /= 'not yet set in get_subdd') call
     &           write_attr(grid,fid,trim(qtyname),'long_name',long_name
     &           )
          end if
        endif

        call par_enddef(grid,fid)
        call write_xy_coord_subdd(fid)
        call write_data(grid,fid,'level',vertical)
      else
        fid = par_open(grid,trim(fname),'write')
      endif

      call write_time_coord_subdd(fid,rec,time,qinst,q24,calendarstring)

c write physical variable
      if(polefix) then
        do n=1,size(data,4)
          do l=1,size(data,3)
            if(hasSouthPole(grid)) data(2:im,1,l,n) = data(1,1,l,n)
            if(hasNorthPole(grid)) data(2:im,jm,l,n) = data(1,jm,l,n)
          end do
        end do
      end if
      if(present(suffixes)) then
        do n=1,size(data,4)
          qname = trim(qtyname)//'_'//trim(suffixes(n))
          call write_dist_data(grid,fid,trim(qname),data(:,:,:,n),
     &         record=rec)
        end do
      else
        call write_dist_data(grid,fid, trim(qtyname), data, record=rec)
      end if
      call par_close(grid,fid)

      return
      end subroutine write_4d

c in_subdd_list
      logical function in_subdd_list(qtyname)
!@sum in_subdd_list tests presence of qtyname in subdd list
!@auth M. Kelley

      implicit none

      character(len=*) :: qtyname
      integer kq
      do kq=1,kdd
        if(trim(qtyname).eq.trim(namedd(kq))) then
          in_subdd_list = .true.
          return
        endif
      enddo
      in_subdd_list = .false.
      end function in_subdd_list

#endif /* NEW_IO_SUBDD */

      end module subdaily

#endif /* not using CACHED_SUBDD */

#ifdef CUBED_SPHERE
      subroutine get_vorticity(vortl)
      use model_com, only : lm
      use domain_decomp_atm, only : grid,get
      use cs2ll_utils, only : uv_derivs_cs_agrid
      use strat, only : dfm_type ! temporarily borrowing
      use constant, only : radius
      use dynamics, only : ualij,valij
      implicit none
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) :: vortl
      integer :: l, i_0,i_1,j_0,j_1
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: ul,vl
      call get(grid, i_strt=i_0,i_stop=i_1, j_strt=j_0,j_stop=j_1)
      do l=1,lm
        ul(:,:) = ualij(l,:,:)
        vl(:,:) = valij(l,:,:)
        call uv_derivs_cs_agrid(grid,dfm_type,ul,vl,vort=vortl(:,:,l))
        vortl(i_0:i_1,j_0:j_1,l) = vortl(i_0:i_1,j_0:j_1,l)/radius
      enddo
      return
      end subroutine get_vorticity
#else
      subroutine get_vorticity(avt)
      use model_com, only : im,jm,lm,u,v
      use domain_decomp_atm, only : grid,hassouthpole,hasnorthpole,
     *     halo_update,north
      use geom, only : dxv,dyp,bydxyp
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) ::
     *     avt
      integer :: i,j,l
      if(hassouthpole(grid)) avt(:, 1,:)=0.
      if(hasnorthpole(grid)) avt(:,jm,:)=0.
      call halo_update(grid,u,from=north)
      call halo_update(grid,v,from=north)
      do l=1,lm
        do j=grid%j_strt_skp,grid%j_stop_skp
          i=1
            avt(i,j,l)=
     *         (((u(i,j,l)+u(im,j,l))/2.*DXV(J)-
     *         (u(i,j+1,l)+u(im,j+1,l))/2.*DXV(J+1))
     *         +((v(i,j,l)+v(i,j+1,l))/2.-(v(im,j,l)+
     *         v(im,j+1,l))/2.)*DYP(J))*BYDXYP(J)
          do i=2,im
            avt(i,j,l)=
     *           (((u(i,j,l)+u(i-1,j,l))/2.*DXV(J)-
     *           (u(i,j+1,l)+u(i-1,j+1,l))/2.*DXV(J+1))
     *           +((v(i,j,l)+v(i,j+1,l))/2.-(v(i-1,j,l)+
     *           v(i-1,j+1,l))/2.)*DYP(J))*BYDXYP(J)
          end do
        end do
      end do
      return
      end subroutine get_vorticity
#endif  /* CUBED_SPHERE */

      subroutine ahourly
!@sum ahourly saves instantaneous variables at sub-daily frequency
!@+   for diurnal cycle diagnostics
!@auth Reha Cakmur/Jan Perlwitz

      USE MODEL_COM, only : u,v,t,p,q,jdate,jhour,ptop,sig
      USE CONSTANT, only : bygrav
      USE domain_decomp_atm, ONLY : am_i_root,get,globalsum,grid
      USE GEOM, only : imaxj,axyp,byaxyp
      USE DYNAMICS, only : phi,wsave,pek,byam
      USE rad_com,ONLY : cosz1,srnflb_save,trnflb_save,ttausv_save,
     &     ttausv_cs_save
      USE diag_com,ONLY : adiurn_dust,ndiupt,ndiuvar,lmax_dd2,ijdd
     &     ,adiurn=>adiurn_loc
#ifndef NO_HDIURN
     &     ,hdiurn=>hdiurn_loc
#endif
#ifdef TRACERS_DUST
     *     ,idd_u1,idd_v1,idd_uv1,idd_t1,idd_qq1,idd_p1,idd_w1,idd_phi1
     *     ,idd_sr1,idd_tr1,idd_load1,idd_conc1,idd_tau1,idd_tau_cs1
#endif
#ifdef TRACERS_ON
      USE TRACER_COM, only : trm
#ifdef TRACERS_DUST
     &     ,Ntm_dust,n_clay
#endif
#endif

      IMPLICIT NONE

      INTEGER :: i,j,ih,ihm,kr,n,n1
      REAL*8 :: psk
      INTEGER,PARAMETER :: n_idxd=14*lmax_dd2
      INTEGER :: idxd(n_idxd)
      REAL*8 :: tmp(NDIUVAR)

C****   define local grid
      INTEGER J_0, J_1, I_0, I_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL get(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

#ifdef TRACERS_DUST
      IF (adiurn_dust == 1) THEN

        idxd=(/
     *       (idd_u1+i-1,i=1,lmax_dd2), (idd_v1+i-1,i=1,lmax_dd2),
     *       (idd_uv1+i-1,i=1,lmax_dd2), (idd_t1+i-1,i=1,lmax_dd2),
     *       (idd_qq1+i-1,i=1,lmax_dd2), (idd_p1+i-1,i=1,lmax_dd2),
     *       (idd_w1+i-1,i=1,lmax_dd2), (idd_phi1+i-1,i=1,lmax_dd2),
     *       (idd_sr1+i-1,i=1,lmax_dd2), (idd_tr1+i-1,i=1,lmax_dd2),
     *       (idd_load1+i-1,i=1,lmax_dd2), (idd_conc1+i-1,i=1,lmax_dd2),
     *       (idd_tau1+i-1,i=1,lmax_dd2), (idd_tau_cs1+i-1,i=1,lmax_dd2)
     *       /)

      END IF
#endif

      ih=jhour+1
      ihm=ih+(jdate-1)*24
!$OMP PARALLEL DO PRIVATE(i,j,kr,n,psk,n1,tmp)
!$OMP*   SCHEDULE(DYNAMIC,2)
      do j=j_0,j_1
      do i=I_0,imaxj(j)
      psk=pek(1,i,j)
      do kr=1,ndiupt
        if(i.eq.ijdd(1,kr).and.j.eq.ijdd(2,kr)) then
#ifdef TRACERS_DUST
          IF (adiurn_dust == 1) THEN
            tmp=0.D0

            tmp(idd_u1:idd_u1+lmax_dd2-1)=u(i,j,1:lmax_dd2)
            tmp(idd_v1:idd_v1+lmax_dd2-1)=v(i,j,1:lmax_dd2)
            tmp(idd_uv1:idd_uv1+lmax_dd2-1)=sqrt(u(i,j,1:lmax_dd2)*u(i
     *           ,j,1:lmax_dd2)+v(i,j,1:lmax_dd2)*v(i,j,1:lmax_dd2))
            tmp(idd_t1:idd_t1+lmax_dd2-1)=t(i,j,1:lmax_dd2)*psk
            tmp(idd_qq1:idd_qq1+lmax_dd2-1)=q(i,j,1:lmax_dd2)
            tmp(idd_p1:idd_p1+lmax_dd2-1)=p(i,j)*sig(1:lmax_dd2)+ptop
            tmp(idd_w1:idd_w1+lmax_dd2-1)=wsave(i,j,1:lmax_dd2)
            tmp(idd_phi1:idd_phi1+lmax_dd2-1)=phi(i,j,1:lmax_dd2)*bygrav
            tmp(idd_sr1:idd_sr1+lmax_dd2-1)=srnflb_save(i,j,1:lmax_dd2)
     *           *cosz1(i,j)
            tmp(idd_tr1:idd_tr1+lmax_dd2-1)=trnflb_save(i,j,1:lmax_dd2)

            DO n=1,Ntm_dust
              n1=n_clay+n-1

              tmp(idd_load1:idd_load1+lmax_dd2-1)
     *             =tmp(idd_load1:idd_load1+lmax_dd2-1)+trm(i,j
     *             ,1:lmax_dd2,n1)*byaxyp(i,j)
              tmp(idd_conc1:idd_conc1+lmax_dd2-1)
     *             =tmp(idd_conc1:idd_conc1+lmax_dd2-1)+trm(i,j
     *             ,1:lmax_dd2,n1)*byam(1,i,j)*byaxyp(i,j)
              tmp(idd_tau1:idd_tau1+lmax_dd2-1)=tmp(idd_tau1:idd_tau1
     *             +lmax_dd2-1)+ttausv_save(i,j,n1,1:lmax_dd2)
              tmp(idd_tau_cs1:idd_tau_cs1+lmax_dd2-1)
     *             =tmp(idd_tau_cs1:idd_tau_cs1+lmax_dd2-1)
     *             +ttausv_cs_save(i,j,n1,1:lmax_dd2)

            END DO

            ADIURN(idxd(:),kr,ih)=ADIURN(idxd(:),kr,ih)+tmp(idxd(:))
#ifndef NO_HDIURN
            HDIURN(idxd(:),kr,ihm)=HDIURN(idxd(:),kr,ihm)+tmp(idxd(:))
#endif

          END IF
#endif
        endif
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end subroutine ahourly

      module msu_wts_mod
      implicit none
      save
      integer, parameter :: nmsu=200 , ncols=4
      real*8 plbmsu(nmsu),wmsu(ncols,nmsu)
      contains
      subroutine read_msu_wts
      use filemanager
      integer n,l,iu_msu
c**** read in the MSU weights file
      call openunit('MSU_wts',iu_msu,.false.,.true.)
      do n=1,4
        read(iu_msu,*)
      end do
      do l=1,nmsu
        read(iu_msu,*) plbmsu(l),(wmsu(n,l),n=1,ncols)
      end do
      call closeunit(iu_msu)
      end subroutine read_msu_wts
      end module msu_wts_mod

      subroutine diag_msu(pland,ts,tlm,ple,tmsu2,tmsu3,tmsu4)
!@sum diag_msu computes MSU channel 2,3,4 temperatures as weighted means
!@auth Reto A Ruedy (input file created by Makiko Sato)
      USE MODEL_COM, only : lm
      use msu_wts_mod
      implicit none
      real*8, intent(in) :: pland,ts,tlm(lm),ple(lm+1)
      real*8, intent(out) :: tmsu2,tmsu3,tmsu4

      real*8 tlmsu(nmsu),tmsu(ncols)
      real*8 plb(0:lm+2),tlb(0:lm+2)
      integer l

c**** find edge temperatures (assume continuity and given means)
      tlb(0)=ts ; plb(0)=plbmsu(1) ; tlb(1)=ts
      plb(1:lm+1) = ple(1:lm+1)
      do l=1,lm
        tlb(l+1)=2*tlm(l)-tlb(l)
      end do
      tlb(lm+2)=tlb(lm+1) ; plb(lm+2)=0.
      call vntrp1 (lm+2,plb,tlb, nmsu-1,plbmsu,tlmsu)
c**** find MSU channel 2,3,4 temperatures
      tmsu(:)=0.
      do l=1,nmsu-1
        tmsu(:)=tmsu(:)+tlmsu(l)*wmsu(:,l)
      end do
      tmsu2 = (1-pland)*tmsu(1)+pland*tmsu(2)
      tmsu3 = tmsu(3)
      tmsu4 = tmsu(4)

      return
      end subroutine diag_msu

      SUBROUTINE init_DIAG(istart,num_acc_files)
!@sum  init_DIAG initializes the diagnostics
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sday,kapa,undef
      USE MODEL_COM, only : lm,Itime,ItimeI,Itime0,pmtop,nfiltr,jhour
     *     ,jdate,jmon,amon,jyear,jhour0,jdate0,jmon0,amon0,jyear0,idacc
     *     ,ioread_single,xlabel,iowrite_single,iyear1,nday,dtsrc,dt
     *     ,nmonav,ItimeE,lrunid,focean,pednl00,pmidl00,lm_req
      USE GEOM, only : axyp,imaxj,lon2d_dg,lat2d_dg
      USE GEOM, only : lonlat_to_ij
      USE SEAICE_COM, only : rsi
      USE LAKES_COM, only : flake
      USE DIAG_COM, only : TSFREZ => TSFREZ_loc
      USE DIAG_COM, only : NPTS, NAMDD, NDIUPT, IJDD,LLDD, ISCCP_DIAGS
      USE DIAG_COM, only : monacc, acc_period, keyct, KEYNR, PLE
      USE DIAG_COM, only : PLM, p1000k, icon_AM, NOFM
      USE DIAG_COM, only : PLE_DN, icon_KE, NSUM_CON, IA_CON, SCALE_CON
      USE DIAG_COM, only : TITLE_CON, PSPEC, LSTR, NSPHER, KLAYER
      USE DIAG_COM, only : ISTRAT, kgz, pmb, kgz_max
      USE DIAG_COM, only : TF_DAY1, TF_LAST, TF_LKON, TF_LKOFF
      USE DIAG_COM, only : name_consrv, units_consrv, lname_consrv
      USE DIAG_COM, only : CONPT0, icon_MS, icon_TPE, icon_WM, icon_EWM
      USE DIAG_COM, only : nreg,jreg,titreg,namreg,sarea_reg
      USE diag_com,ONLY : adiurn_dust,adiurn_loc,areg_loc,aisccp_loc
     &     ,consrv_loc
#ifndef NO_HDIURN
     &     ,hdiurn_loc
#endif
      USE diag_com,only : lh_diags
#ifdef TES_LIKE_DIAGS
      USE DIAG_COM, only : kgz_max_more,KGZmore,pmbmore
#endif
      USE DIAG_LOC
      USE PARAM
      USE FILEMANAGER
      USE DOMAIN_DECOMP_ATM, only: GRID,GET,WRITE_PARALLEL,
     &     AM_I_ROOT,GLOBALSUM
      use msu_wts_mod
      IMPLICIT NONE
      integer, intent(in) :: istart,num_acc_files
      INTEGER I,J,L,K,KL,n,ioerr,months,years,mswitch,ldate
     *     ,jday0,jday,moff,kb,l850,l300,l50
      REAL*8 PLE_tmp
      CHARACTER CONPT(NPTS)*10
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
      INTEGER :: J_0,J_1, I_0,I_1
!@var out_line local variable to hold mixed-type output for parallel I/O
      character(len=300) :: out_line
      character(len=80) :: filenm
!@var iu_REG unit number for regions file
      INTEGER iu_REG
#ifdef CUBED_SPHERE
#define ASCII_REGIONS
#endif
#ifdef ASCII_REGIONS
C***  regions defined as rectangles in an ASCII file
      integer, dimension(23) :: NRECT
      character*4, dimension(23,6) :: CORLON,CORLAT
      real*8, dimension(23,6) :: DCORLON,DCORLAT   !lat-lon coordinates of rect. corners
      integer :: ireg,irect,icorlon,icorlat
      real*8::lon,lat
#else
c dummy global array to read special-format regions file lacking a
c a parallelized i/o routine that understands it
      integer :: jreg_glob(im,jm)
#endif
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     area_part

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****   READ SPECIAL REGIONS
#ifndef ASCII_REGIONS
      call openunit("REG",iu_REG,.true.,.true.)
      READ(iu_REG) TITREG,JREG_glob,NAMREG
      jreg(i_0:i_1,j_0:j_1) = jreg_glob(i_0:i_1,j_0:j_1)
#else
c**** Regions are defined in reg.txt input file as union of rectangles
c**** independant of resolution & grid type

      call openunit("REG",iu_REG,.false.,.true.)
      READ(iu_REG,'(A80)') TITREG !read title
      READ (iu_REG,'(I2)') (NRECT(I), I=1,23 ) !#of rectangles per region
      READ (iu_REG,'(A4,1X,A4)') (NAMREG(1,I),NAMREG(2,I),I=1,23) !Read region name

c**** Read cordinates of rectangles
c**** (NWcorner long, NW corner lat, SE corner long, SE corner lat)(1:23)
c**** 0555 = no rectangle
      READ (iu_REG,'(11(A4,1X),A4)') (
     &       CORLON(I,1),CORLAT(I,1),
     &       CORLON(I,2),CORLAT(I,2),
     &       CORLON(I,3),CORLAT(I,3),
     &       CORLON(I,4),CORLAT(I,4),
     &       CORLON(I,5),CORLAT(I,5),
     &       CORLON(I,6),CORLAT(I,6), I=1,23)

c**** Convert to integer
      do i=1,23
        do j=1,6
          read(corlon(i,j),'(I4)') icorlon
          dcorlon(i,j) =icorlon
          read(corlat(i,j),'(I4)') icorlat
          dcorlat(i,j) =icorlat
        enddo
      enddo

c**** determine the region to which each cell belongs
      JREG(:,:)=24
      do j=j_0,j_1
        do i=i_0,i_1
          lon=lon2d_dg(i,j)
          lat=lat2d_dg(i,j)
          do ireg=1,23
            do irect=1,NRECT(ireg)
              if ( lat > dcorlat(ireg,2*irect-1) .or.
     &             lat < dcorlat(ireg,2*irect  ) )  cycle
              if(dcorlon(ireg,2*irect-1)<dcorlon(ireg,2*irect)) then
                if( lon < dcorlon(ireg,2*irect-1) .or.
     &              lon > dcorlon(ireg,2*irect  ) ) cycle
              else ! wraparound
                if( lon < dcorlon(ireg,2*irect-1) .and.
     &              lon > dcorlon(ireg,2*irect  ) ) cycle
              endif
              JREG(i,j)=ireg
            enddo
          enddo
        enddo
      enddo
#endif
      IF (AM_I_ROOT()) then
        WRITE(6,*) ' read REGIONS from unit ',iu_REG,': ',TITREG
      endif
      call closeunit(iu_REG)
c
c calculate the areas of the special regions
c
#ifndef SCM
      do n=1,nreg
        where(jreg(i_0:i_1,j_0:j_1).eq.n)
          area_part(i_0:i_1,j_0:j_1) = axyp(i_0:i_1,j_0:j_1)
        elsewhere
          area_part(i_0:i_1,j_0:j_1) = 0d0
        end where
        call globalsum(grid,area_part,sarea_reg(n),all=.true.)
      enddo
#endif

C**** Initialse diurnal diagnostic locations (taken from the 4x5 res)
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      NAMDD =
     &   (/'AUSD', 'MWST', 'SAHL', 'EPAC', 'AF01',
     &     'AF02', 'AF03', 'AF04', 'AF05', 'ASA1',
     &     'AF06', 'AME1', 'AF07', 'AF08', 'AF09',
     &     'ARAB', 'AF10', 'AF11', 'ASA2', 'AUS1',
     &     'AUS2', 'AUS3', 'AF12', 'AF13', 'ASA3',
     &     'AF14', 'ASA4', 'AUS4', 'AME2', 'AF15',
     &     'AME3', 'AF16', 'AF17', 'AF18' /)
      LLDD = RESHAPE( (/
     &   132.5,-26.,  -97.5, 42.,    2.5, 14.,
     &  -117.5, -2.,   32.5, 18.,   22.5, 30.,
     &    -2.5, 26.,  -12.5, 26.,   12.5, 18.,
     &    62.5, 38.,   17.5, 18.,  -67.5,-42.,
     &    -7.5, 26.,   22.5, 26.,   27.5, 30.,
     &    42.5, 34.,   12.5, 14.,   -2.5, 22.,
     &    62.5, 42.,  137.5,-30.,  137.5,-26.,
     &   127.5,-30.,   22.5, 18.,   -2.5, 18.,
     &   102.5, 42.,    7.5, 30.,   72.5, 26.,
     &   132.5,-30.,  -67.5,-38.,    2.5, 26.,
     &   -67.5,-46.,   22.5, 22.,    2.5, 22.,
     &    17.5, 30.  /),(/2,NDIUPT/))
#else
c defaults for diurnal diagnostics
      NAMDD = (/ 'PITT', 'SGPL', 'NYCC', 'HOUS' /)
c Pittsburgh, Southern Great Plains, NYC, Houston
      LLDD = RESHAPE( (/
c        Longitude, Latitude
     &      -79.9,  40.4,
     &      -97.5,  36.4,
     &      -74.0,  40.8,
     &      -95.4,  29.8
     &     /),(/2,4/))
#endif
      call sync_param( "LLDD", LLDD(1:2,1), 2*NDIUPT )
      do n=1,ndiupt
        call lonlat_to_ij(lldd(1,n),ijdd(1,n))
      enddo

      call sync_param( "NAMDD", NAMDD, NDIUPT )
#ifndef CUBED_SPHERE
c if people still want to specify dd points as ij, let them
      call sync_param( "IJDD", IJDD(1:2,1), 2*NDIUPT )
#endif

      call sync_param( "isccp_diags",isccp_diags)
      call sync_param( "adiurn_dust",adiurn_dust)
      call sync_param( "lh_diags",lh_diags)

      IF(ISTART.LT.1) THEN  ! initialize for post-processing
        call getdte(Itime0,Nday,Iyear1,Jyear0,Jmon0,Jday0,Jdate0,Jhour0
     *       ,amon0)
        call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour
     *       ,amon)
        months=1 ; years=monacc(jmon0) ; mswitch=0 ; moff=0 ; kb=jmon0
        do kl=jmon0+1,jmon0+11
          k = kl
          if (k.gt.12) k=k-12
          if (monacc(k).eq.years) then
            months=months+1
          else if (monacc(k).ne.0) then
C****            write(6,*) 'uneven period:',monacc
            CALL WRITE_PARALLEL(monacc, UNIT=6, format=
     &                          "('uneven period:',12I5)")
            call stop_model( 'uneven period', 255 )
          end if
          if(monacc(k).ne.monacc(kb)) mswitch = mswitch+1
          if(mswitch.eq.2) moff = moff+1
          kb = k
        end do
        if (mswitch.gt.2) then
C****          write(6,*) 'non-consecutive period:',monacc
            CALL WRITE_PARALLEL(monacc, UNIT=6, format=
     &                          "('non-consecutive period:',12I5)")
          call stop_model( 'non-consecutive period', 255 )
        end if
        call aPERIOD (JMON0,JYEAR0,months,years,moff, acc_period,Ldate)
        if (num_acc_files.gt.1) then  ! save the summed acc-file
          write(out_line,*) num_acc_files,' files are summed up'
          CALL WRITE_PARALLEL(TRIM(out_line), UNIT=6)
          keyct=1 ; KEYNR=0
          XLABEL(128:132)='     '
          XLABEL(120:132)=acc_period(1:3)//' '//acc_period(4:Ldate)
C****          write(6,*) XLABEL
          CALL WRITE_PARALLEL(XLABEL, UNIT=6)
          filenm=acc_period(1:Ldate)//'.acc'//XLABEL(1:LRUNID)
          call io_rsf (filenm,Itime,iowrite_single,ioerr)
        end if
        ItimeE = -1
        close (6)
        open(6,file=acc_period(1:Ldate)//'.'//XLABEL(1:LRUNID)//'.PRT',
     *       FORM='FORMATTED')
      END IF

C**** Initialize certain arrays used by more than one print routine
      DO L=1,LM
        PLE(L)   =pednl00(l+1)
        PLE_DN(L)=pednl00(l)
        PLM(L)   =pmidl00(l)
      END DO
      PLM(LM+1:LM+LM_REQ)=pmidl00(lm+1:lm+lm_req)

      p1000k=1000.0**kapa

C**** Initialise some local constants (replaces IFIRST constructions)
C**** From DIAGA:
      DO L=1,LM
        LUPA(L)=L+1
        LDNA(L)=L-1
      END DO
      LDNA(1)=1
      LUPA(LM)=LM

C**** From DIAGB (PM, PMO are fixed, PL,PLO will vary)
      PM(1)=1200.   ! ensures below surface for extrapolation
      DO L=2,LM+1
        PL(L)=pednl00(l)
        PM(L)=pednl00(l)
      END DO
      DO L=1,LM
        PLO(L)=pmidl00(l)
        PMO(L)=.5*(PM(L)+PM(L+1))
      END DO

C**** From DIAG7A
      L850=LM
      L300=LM
      L50=LM
      DO L=LM-1,1,-1
        PLE_tmp=.25*(PEDNL00(L)+2.*PEDNL00(L+1)+PEDNL00(L+2))
        IF (PLE_tmp.LT.850.) L850=L
        IF (PLE_tmp.LT.300.) L300=L
        IF (PLE_tmp.LT.250.) JET=L
        IF (PLE_tmp.LT.50.) L50=L
      END DO
C      WRITE (6,888) JET
      CALL WRITE_PARALLEL(JET, UNIT=6, format=
     & "(' JET WIND LEVEL FOR DIAG',I3)")
C 888  FORMAT (' JET WIND LEVEL FOR DIAG',I3)
C****      WRITE (6,889) L850,L300,L50
      WRITE (out_line,889) L850,L300,L50
      CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
 889  FORMAT (' LEVELS FOR WIND WAVE POWER DIAG  L850=',I3,
     *     ' L300=',I3,' L50=',I3)
      LDEX(1)=L850
      LDEX(2)=L300
      LDEX(3)=L50

C**** Initialize conservation diagnostics
C**** NCON=1:25 are special cases: Angular momentum and kinetic energy
      icon_AM=1
      NOFM(:,icon_AM) = (/  1, 8, 0, 0, 0, 0, 9,10, 0,11, 0, 0/)
      icon_KE=2
      NOFM(:,icon_KE) = (/ 13,20,21, 0, 0, 0,22,23, 0,24, 0, 0/)
      NSUM_CON(1:25) = (/-1,-1,-1,-1,-1,-1,-1,12,12,12,12, 0,
     *                   -1,-1,-1,-1,-1,-1,-1,25,25,25,25,25, 0/)
      IA_CON(1:25) =   (/12, 1, 1, 1, 1, 1, 1, 7, 8,10, 9,12,
     *                   12, 1, 1, 1, 1, 1, 1, 7, 8, 8,10, 9,12/)
      SCALE_CON(1)              = 1d-9
      SCALE_CON((/2,3,4,5,6,7,8,9/))= 1d-2/DTSRC
      SCALE_CON(10)              = 1d-2/(NFILTR*DTSRC)
      SCALE_CON(11)             = 2d-2/SDAY
      SCALE_CON((/12,25/))      = 1.
      SCALE_CON(13)             = 1d-3
      SCALE_CON((/14,15,16,17,18,19,20,21,22/)) = 1d3/DTSRC
      SCALE_CON(23)             = 1d3/(NFILTR*DTSRC)
      SCALE_CON(24)             = 2d3/SDAY
      TITLE_CON(1:25) = (/
     *  ' INSTANTANE AM (10**9 J*S/M^2)  ',
     *  '     DELTA AM BY ADVECTION      ',
     *  '     DELTA AM BY CORIOLIS FORCE ',
     *  '     DELTA AM BY PRESSURE GRAD  ',
     *  '     DELTA AM BY STRATOS DRAG   ',
     *  '     DELTA AM BY UV FILTER      ',
     *  '     DELTA AM BY GW DRAG        ',
     *  ' CHANGE OF AM BY DYNAMICS       ',
     *  ' CHANGE OF AM BY SURF FRIC+TURB ',
     *  ' CHANGE OF AM BY FILTER         ',
     *  ' CHANGE OF AM BY DAILY RESTOR   ',
     *  ' SUM OF CHANGES (10**2 J/M^2)   ',
     *  '0INSTANTANEOUS KE (10**3 J/M^2) ',
     *  '     DELTA KE BY ADVECTION      ',
     *  '     DELTA KE BY CORIOLIS FORCE ',
     *  '     DELTA KE BY PRESSURE GRAD  ',
     *  '     DELTA KE BY STRATOS DRAG   ',
     *  '     DELTA KE BY UV FILTER      ',
     *  '     DELTA KE BY GW DRAG        ',
     *  ' CHANGE OF KE BY DYNAMICS       ',
     *  ' CHANGE OF KE BY MOIST CONVEC   ',
     *  ' CHANGE OF KE BY SURF + DC/TURB ',
     *  ' CHANGE OF KE BY FILTER         ',
     *  ' CHANGE OF KE BY DAILY RESTOR   ',
     *  ' SUM OF CHANGES (10**-3 W/M^2)  '/)
      name_consrv(1:25) = (/
     *     'inst_AM   ','del_AM_ADV','del_AM_COR','del_AM_PRE',
     *     'del_AM_STR','del_AM_UVF','del_AM_GWD','chg_AM_DYN'
     *     ,'chg_AM_SUR','chg_AM_FIL','chg_AM_DAI','sum_chg_AM'
     *     ,'inst_KE   ','del_KE_ADV','del_KE_COR','del_KE_PRE'
     *     ,'del_KE_STR','del_KE_UVF','del_KE_GWD','chg_KE_DYN'
     *     ,'chg_KE_MOI','chg_KE_SUR','del_KE_FIL','chg_KE_DAI'
     *     ,'sum_chg_KE'/)
      units_consrv(1)    ="10**9 J*S/M^2"
      units_consrv(2:12) ="10**2 J/M^2"
      units_consrv(13)   ="10**3 J/M^2"
      units_consrv(14:24)="10**-3 W/M^2"
      lname_consrv(1:25)=TITLE_CON(1:25)
C**** To add a new conservation diagnostic:
C****    i) Add 1 to NQUANT, and increase KCON in DIAG_COM.f
C****   ii) Set up a QCON, and call SET_CON to allocate array numbers,
C****       set up scales, titles, etc. The icon_XX index must be
C****       declared in DIAG_COM.f for the time being
C**** QCON denotes when the conservation diags should be done
C**** 1:NPTS ==> DYN,   COND,   RAD,   PREC,   LAND,  SURF,
C****            FILTER,STRDG/OCEAN, DAILY, OCEAN1, OCEAN2,
C****  iii) Write a conserv_XYZ routine that returns the zonal average
C****       of your quantity
C****   iv) Add a line to DIAGCA that calls conserv_DIAG (declared
C****       as external)
C****    v) Note that the conserv_XYZ routine, and call to SET_CON
C****       should be in the driver module for the relevant physics

C**** Set up atmospheric component conservation diagnostics
      CONPT=CONPT0
C**** Atmospheric mass
      QCON=(/ T, F, F, F, F, F, T, F, T, F, F/)
      CALL SET_CON(QCON,CONPT,"MASS    ","(KG/M^2)        ",
     *     "(10**-8 KG/SM^2)",1d0,1d8,icon_MS)
C**** Atmospheric total potential energy
      CONPT(8)="SURF+TURB" ; CONPT(6)="KE DISSIP"
      QCON=(/ T, T, T, F, F, T, T, T, F, F, F/)
      CALL SET_CON(QCON,CONPT,"TPE     ","(10**5 J/M^2)   ",
     *     "(10**-2 W/M^2)  ",1d-5,1d2,icon_TPE)
C**** Atmospheric water mass
      CONPT(6)="SURF+TURB"
      QCON=(/ T, T, F, F, F, T, F, F, F, F, F/)
      CALL SET_CON(QCON,CONPT,"ATM WAT ","(10**-2 KG/M^2) ",
     *     "(10**-8 KG/SM^2)",1d2,1d8,icon_WM)
C**** Atmospheric water latent heat (at some point should include
C**** sensible + potential energy associated with water mass as well)
      QCON=(/ T, T, F, F, F, T, F, F, F, F, F/)
      CALL SET_CON(QCON,CONPT,"ENRG WAT","(10**3 J/M^2)   ",
     *     "(10**-2 W/M^2)  ",1d-3,1d2,icon_EWM)

C**** Initialize layering for spectral diagnostics
C**** add in epsilon=1d-5 to avoid roundoff mistakes
      KL=1
      DO L=1,LM
        IF (PEDNL00(L+1)+1d-5.lt.PSPEC(KL) .and.
     *      PEDNL00(L)  +1d-5.gt.PSPEC(KL)) THEN
          IF (KL.eq.2) LSTR = L  ! approx. 10mb height
          KL=KL+1
        END IF
        KLAYER(L)=4*(KL-1)+1
      END DO
      IF (KL*4 .gt. NSPHER) THEN
        CALL WRITE_PARALLEL("Inconsistent definitions of stratosphere:"
     &                       ,UNIT=6)
        CALL WRITE_PARALLEL("Adjust PSPEC, ISTRAT so that KL*4 = NSPHER"
     &                       ,UNIT=6)
        WRITE(out_line,*) "ISTRAT,PSPEC,NSPHER,KL=",
     &                     ISTRAT,PSPEC,NSPHER,KL
        CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
        call stop_model(
     *    "Stratospheric definition problem for spectral diags.",255)
      END IF

C**** Calculate the max number of geopotential heights
      do k=1,kgz
        if (pmb(k).le.pmtop) exit
        kgz_max = k
      end do
#ifdef TES_LIKE_DIAGS
      do k=1,KGZmore
        if (pmbmore(k).le.pmtop) exit
        kgz_max_more = k
      end do
#endif
      CALL WRITE_PARALLEL(" Geopotential height diagnostics at (mb): ",
     &                      UNIT=6)
      CALL WRITE_PARALLEL(PMB(1:kgz_max), UNIT=6, format="(20F9.3)")

c**** Initialize acc-array names, units, idacc-indices
      call def_acc

C**** Ensure that diagnostics are reset at the beginning of the run
      IF (Itime.le.ItimeI .and. ISTART.gt.0) THEN
        call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour
     *       ,amon)
        CALL reset_DIAG(0)
C**** Initiallise ice freeze diagnostics at beginning of run
        DO J=J_0,J_1
          DO I=I_0,IMAXJ(J)
            TSFREZ(I,J,TF_DAY1)=365.
            TSFREZ(I,J,TF_LAST)=365.
            IF (FOCEAN(I,J)+FLAKE(I,J).gt.0) then
              IF (RSI(I,J).gt.0) then
                TSFREZ(I,J,TF_LKON) = JDAY-1
                TSFREZ(I,J,TF_LKOFF) = JDAY
              ELSE
                TSFREZ(I,J,TF_LKON) = JDAY
                TSFREZ(I,J,TF_LKOFF) = undef
              END IF
            ELSE
              TSFREZ(I,J,TF_LKON) = undef
              TSFREZ(I,J,TF_LKOFF) = undef
            END IF
          END DO
        END DO
        CALL daily_DIAG
      END IF
c
c zero out certain non-distributed arrays
c
      areg_loc = 0
      aisccp_loc = 0
      consrv_loc = 0
      adiurn_loc = 0
#ifndef NO_HDIURN
      hdiurn_loc = 0
#endif

c
c read MSU weighting functions for diagnostics
c
      call read_msu_wts

      RETURN
      END SUBROUTINE init_DIAG


      SUBROUTINE reset_DIAG(isum)
!@sum  reset_DIAG resets/initializes diagnostics
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : Itime,iyear1,nday,kradia,
     *     Itime0,jhour0,jdate0,jmon0,amon0,jyear0,idacc,u
      USE DIAG_COM
      USE PARAM
      USE DOMAIN_DECOMP_ATM, only: grid,am_i_root
      IMPLICIT NONE
      INTEGER :: isum !@var isum if =1 preparation to add up acc-files
      INTEGER jd0

      IDACC(1:12)=0
      if (kradia.gt.0) then
        AFLX_ST = 0.
        if (isum.eq.1) return
        go to 100
      end if

      if(am_i_root()) then
        aj = 0; ajl = 0; asjl = 0; agc = 0
      endif
      AJ_loc=0    ; AREG_loc=0; AREG=0
      AJL_loc=0  ; ASJL_loc=0   ; AIJ_loc=0
      AIJL_loc=0   ; ENERGY=0 ; CONSRV = 0; CONSRV_loc=0
      SPECA=0 ; ATPE=0 ; WAVE=0 ; AGC_loc=0   ; AIJK_loc=0
#ifndef NO_HDIURN
      HDIURN=0; HDIURN_loc=0
#endif
      ADIURN=0 ; ADIURN_loc=0; AISCCP=0; AISCCP_loc=0
#ifdef TRACERS_ON
      call reset_trdiag
#endif
      call reset_ODIAG(isum)  ! ocean diags if required
      call reset_icdiag       ! ice dynamic diags if required

      if (isum.eq.1) return ! just adding up acc-files

      AIJ_loc(:,:,IJ_TMNMX)=1000. ; IDACC(12)=1

#ifndef CUBED_SPHERE
      CALL EPFLXI (U)  ! strat
#endif

  100 Itime0=Itime
      call getdte(Itime0,Nday,Iyear1,Jyear0,Jmon0,Jd0,
     *     Jdate0,Jhour0,amon0)

      RETURN
      END SUBROUTINE reset_DIAG


      SUBROUTINE daily_DIAG
!@sum  daily_DIAG resets diagnostics at beginning of each day
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : undef
      USE MODEL_COM, only : im,jm,jday,focean
      USE GEOM, only : imaxj,lat2d
      USE SEAICE_COM, only : rsi
      USE LAKES_COM, only : flake
      USE GHY_COM, only : fearth
      USE DIAG_COM, only : aij=>aij_loc
     *     ,ij_lkon,ij_lkoff,ij_lkice,tsfrez=>tsfrez_loc,tdiurn
     *     ,tf_lkon,tf_lkoff,tf_day1,tf_last
      USE DOMAIN_DECOMP_ATM, only : GRID,GET,am_i_root
#ifdef TRACERS_ON
      USE RAD_COM,only: ttausv_sum,ttausv_sum_cs,ttausv_count
#endif
      IMPLICIT NONE
      INTEGER I,J
      INTEGER :: J_0, J_1, I_0,I_1

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C**** INITIALIZE SOME ARRAYS AT THE BEGINNING OF SPECIFIED DAYS
      IF (JDAY.EQ.32) THEN
        DO J=J_0,J_1
        DO I=I_0,I_1
          if(lat2d(i,j).gt.0.) then
            TSFREZ(I,J,TF_DAY1)=JDAY
          else
            TSFREZ(I,J,TF_LAST)=JDAY
          endif
        ENDDO
        ENDDO
      ELSEIF (JDAY.EQ.213) THEN
        DO J=J_0,J_1
        DO I=I_0,I_1
          if(lat2d(i,j).lt.0.) then
            TSFREZ(I,J,TF_DAY1)=JDAY
          endif
        ENDDO
        ENDDO
      END IF
C**** set and initiallise freezing diagnostics
C**** Note that TSFREZ saves the last day of no-ice and some-ice.
C**** The AIJ diagnostics are set once a year (zero otherwise)
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          if(lat2d(i,j).lt.0.) then
C**** initialize/save South. Hemi. on Feb 28
            IF (JDAY.eq.59 .and. TSFREZ(I,J,TF_LKOFF).ne.undef) THEN
              AIJ(I,J,IJ_LKICE)=1.
              AIJ(I,J,IJ_LKON) =MOD(NINT(TSFREZ(I,J,TF_LKON)) +307,365)
              AIJ(I,J,IJ_LKOFF)=MOD(NINT(TSFREZ(I,J,TF_LKOFF))+306,365)
     *             +1
              IF (RSI(I,J).gt.0) THEN
                TSFREZ(I,J,TF_LKON) = JDAY-1
              ELSE
                TSFREZ(I,J,TF_LKOFF) = undef
              END IF
            END IF
          ELSE
C**** initiallise/save North. Hemi. on Aug 31
C**** Note that for continuity across the new year, the julian days
C**** are counted from Sep 1 (NH only).
            IF (JDAY.eq.243 .and. TSFREZ(I,J,TF_LKOFF).ne.undef) THEN
              AIJ(I,J,IJ_LKICE)=1.
              AIJ(I,J,IJ_LKON) =MOD(NINT(TSFREZ(I,J,TF_LKON)) +123,365)
              AIJ(I,J,IJ_LKOFF)=MOD(NINT(TSFREZ(I,J,TF_LKOFF))+122,365)
     *             +1
              IF (RSI(I,J).gt.0) THEN
                TSFREZ(I,J,TF_LKON) = JDAY-1
              ELSE
                TSFREZ(I,J,TF_LKOFF) = undef
              END IF
            END IF
          END IF
C**** set ice on/off days
          IF (FOCEAN(I,J)+FLAKE(I,J).gt.0) THEN
            IF (RSI(I,J).eq.0.and.TSFREZ(I,J,TF_LKOFF).eq.undef)
     *           TSFREZ(I,J,TF_LKON)=JDAY
            IF (RSI(I,J).gt.0) TSFREZ(I,J,TF_LKOFF)=JDAY
          END IF
        END DO
      END DO

C**** INITIALIZE SOME ARRAYS AT THE BEGINNING OF EACH DAY
      DO J=J_0,J_1
         DO I=I_0,I_1
            TDIURN(I,J,1)= 1000.
            TDIURN(I,J,2)=-1000.
            TDIURN(I,J,3)= 1000.
            TDIURN(I,J,4)=-1000.
            TDIURN(I,J,5)=    0.
            TDIURN(I,J,6)=-1000.
            TDIURN(I,J,7)=-1000.
            TDIURN(I,J,8)=-1000.
            TDIURN(I,J,9)= 1000.
            TDIURN(I,J,10)= 1000.
            TDIURN(I,J,11)=-1000.
            TDIURN(I,J,12)=-1000.
            IF (FEARTH(I,J).LE.0.) THEN
               TSFREZ(I,J,TF_DAY1)=365.
               TSFREZ(I,J,TF_LAST)=365.
            END IF
#ifdef TRACERS_ON
            ttausv_sum(I,J,:)=0.d0
            ttausv_sum_cs(I,J,:)=0.d0
#endif
         END DO
      END DO
#ifdef TRACERS_ON
      ttausv_count=0.d0
#endif
      END SUBROUTINE daily_DIAG


      SUBROUTINE SET_CON(QCON,CONPT,NAME_CON,INST_UNIT,SUM_UNIT,INST_SC
     *     ,CHNG_SC,ICON)
!@sum  SET_CON assigns conservation diagnostic array indices
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sday
      USE MODEL_COM, only : dtsrc,nfiltr
      USE DIAG_COM, only : kcon,nquant,npts,title_con,scale_con,nsum_con
     *     ,nofm,ia_con,kcmx,ia_d5d,ia_d5s,ia_filt,ia_12hr,ia_inst
     *     ,name_consrv,lname_consrv,units_consrv
      USE DOMAIN_DECOMP_ATM, only : WRITE_PARALLEL
      IMPLICIT NONE
!@var QCON logical variable sets where conservation diags are saved
      LOGICAL, INTENT(IN),DIMENSION(NPTS) :: QCON
!@var CONPT names for points where conservation diags are saved
      CHARACTER*10, INTENT(IN),DIMENSION(NPTS) :: CONPT
!@var INST_SC scale for instantaneous value
      REAL*8, INTENT(IN) :: INST_SC
!@var CHNG_SC scale for changes
      REAL*8, INTENT(IN) :: CHNG_SC
!@var NAME_CON name of conservation quantity
      CHARACTER*8, INTENT(IN) :: NAME_CON
!@var sname name of conservation quantity (no spaces)
      CHARACTER*8 :: sname
!@var INST_UNIT string for unit for instant. values
      CHARACTER*16, INTENT(IN) :: INST_UNIT
!@var SUM_UNIT string for unit for summed changes
      CHARACTER*16, INTENT(IN) :: SUM_UNIT
!@var ICON index for the conserved quantity
      INTEGER, INTENT(OUT) :: ICON
!@var out_line local variable to hold mixed-type output for parallel I/O
      character(len=300) :: out_line
!@var CONPT_us CONPT with blanks replaced by underscores
      CHARACTER*10 :: CONPT_us
      CHARACTER*40 :: clean_str
      INTEGER NI,NM,NS,N,k
      INTEGER, SAVE :: NQ = 2   ! first 2 special cases AM + KE

      NQ=NQ+1
      IF (NQ.gt.NQUANT) THEN
        WRITE(out_line,*)
     *       "Number of conserved quantities larger than NQUANT"
     *       ,NQUANT,NQ
        CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
        call stop_model("Change NQUANT in diagnostic common block",255)
      END IF
C**** make nice netcdf names
      sname=trim(clean_str(name_CON))
C****
      NI=KCMX+1
      NOFM(1,NQ) = NI
      TITLE_CON(NI) = "0INSTANT "//TRIM(NAME_CON)//" "//TRIM(INST_UNIT)
      SCALE_CON(NI) = INST_SC
      name_consrv(NI) ="inst_"//sname
      lname_consrv(NI) = "INSTANT "//TRIM(NAME_CON)
      units_consrv(NI) = INST_UNIT
      IA_CON(NI) = ia_inst
      NM=NI
      DO N=1,NPTS
        IF (QCON(N)) THEN
          NM = NM + 1
          NOFM(N+1,NQ) = NM
          TITLE_CON(NM) = " CHANGE OF "//TRIM(NAME_CON)//" BY "//
     *         CONPT(N)
          CONPT_us = CONPT(N)
          do k=1,len_trim(CONPT_us)
            if (CONPT_us(k:k).eq." ") CONPT_us(k:k)="_"
          end do
          name_consrv(NM) ="chg_"//trim(sname)//"_"//TRIM(CONPT_us)
          lname_consrv(NM) = TITLE_CON(NM)
          units_consrv(NM) = SUM_UNIT
          SELECT CASE (N)
          CASE (1)
            SCALE_CON(NM) = CHNG_SC/DTSRC
            IA_CON(NM) = ia_d5d
          CASE (2,3,4,5,6,8,10,11)
            SCALE_CON(NM) = CHNG_SC/DTSRC
            IA_CON(NM) = ia_d5s
          CASE (7)
            SCALE_CON(NM) = CHNG_SC/(NFILTR*DTSRC)
            IA_CON(NM) = ia_filt
          CASE (9)
            SCALE_CON(NM) = CHNG_SC*2./SDAY
            IA_CON(NM) = ia_12hr
          END SELECT
        ELSE
          NOFM(N+1,NQ) = 0
        END IF
      END DO
      NS=NM+1
      IF (NS.gt.KCON) THEN
        WRITE(out_line,*)
     *      "KCON not large enough for extra conserv diags",
     *       KCON,NI,NM,NQ,NS,NAME_CON
        CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
        call stop_model("Change KCON in diagnostic common block",255)
      END IF
      TITLE_CON(NS) = " SUM OF CHANGES "//TRIM(SUM_UNIT)
      name_consrv(NS) ="sum_chg_"//trim(sname)
      lname_consrv(NS) = " SUM OF CHANGES OF "//TRIM(NAME_CON)
      units_consrv(NS) = SUM_UNIT
      SCALE_CON(NS) = 1.
      IA_CON(NS) = ia_inst
      NSUM_CON(NI) = -1
      NSUM_CON(NI+1:NS-1) = NS
      NSUM_CON(NS) = 0
      KCMX=NS
      ICON=NQ
      RETURN
C****
      END SUBROUTINE set_con

      SUBROUTINE UPDTYPE
!@sum UPDTYPE updates FTYPE array to ensure correct budget diagnostics
!@auth Gavin Schmidt
      USE MODEL_COM, only : im,jm,focean,flice,itocean
     *     ,itoice,itlandi,itearth,itlake,itlkice,ftype
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : rsi
      USE LAKES_COM, only : flake
      USE GHY_COM, only : fearth
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      IMPLICIT NONE
      INTEGER I,J
      INTEGER :: J_0,J_1,I_0,I_1

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          FTYPE(ITOICE ,I,J)=FOCEAN(I,J)*RSI(I,J)
          FTYPE(ITOCEAN,I,J)=FOCEAN(I,J)-FTYPE(ITOICE,I,J)
          FTYPE(ITLKICE,I,J)=FLAKE(I,J)*RSI(I,J)
          FTYPE(ITLAKE ,I,J)=FLAKE(I,J)-FTYPE(ITLKICE,I,J)
C**** set land components of FTYPE array. Summation is necessary for
C**** cases where Earth and Land Ice are lumped together
          FTYPE(ITLANDI,I,J)=0.
          FTYPE(ITEARTH,I,J)=FEARTH(I,J)
          FTYPE(ITLANDI,I,J)=FTYPE(ITLANDI,I,J)+FLICE(I,J)
        END DO
      END DO
      RETURN
C****
      END SUBROUTINE UPDTYPE

      subroutine calc_derived_aij
!@sum Calculate derived lat/lon diagnostics prior to printing
!@auth Group
      USE CONSTANT, only : grav,rgas,bygrav,tf,teeny
      USE DOMAIN_DECOMP_ATM, only : GRID,SUMXPE,AM_I_ROOT
      USE Domain_decomp_1d, only: hasSouthPole, hasNorthPole
      USE MODEL_COM, only : idacc,zatmo,fearth0,flice,focean,lm,pmtop,
     &     im,jm
      USE GEOM, only : imaxj,axyp,lat2d,areag
      USE DIAG_COM, only : aij=>aij_loc,tsfrez=>tsfrez_loc,
     &     kaij,hemis_ij,jgrid_ij,
     &     aijl=>aijl_loc,ia_ij,ia_src,ia_inst,ia_dga,tf_last,tf_day1,
     *     ij_topo, ij_wsmn, ij_wsdir, ij_jet, ij_jetdir, ij_grow,
     *     ij_netrdp, ij_albp, ij_albg, ij_albv,   ij_pwater, ij_lk,
     *     ij_fland, ij_dzt1, ij_albgv, ij_clrsky, ij_pocean, ij_ts,
     *     ij_RTSE, ij_HWV, ij_PVS,
     &     IJ_TRNFP0,IJ_SRNFP0,IJ_TRSUP,IJ_TRSDN,IJ_EVAP,IJ_QS,IJ_PRES,
     &     IJ_SRREF,IJ_SRVIS,IJ_SRINCP0,IJ_SRINCG,IJ_SRNFG,IJ_PHI1K,
     &     IJ_US,IJ_VS,IJ_UJET,IJ_VJET,IJ_CLDCV,IJ_TATM,IJK_DP,IJK_TX,
     &     IJ_MSU2,IJ_MSU3,IJ_MSU4,KGZ_MAX,GHT,PMB,
     &     ij_TminC,ij_TmaxC,ij_TDcomp,
     *     ij_swaerrf,ij_lwaerrf,ij_swaersrf,ij_lwaersrf,ij_swaerabs,
     *     ij_lwaerabs,ij_swaerrfnt,ij_lwaerrfnt,ij_swaersrfnt,
     *     ij_lwaersrfnt,ij_swaerabsnt,ij_lwaerabsnt
      IMPLICIT NONE
      INTEGER :: I,J,L,K,K1,K2,N,KHEM
      INTEGER :: J_0,J_1,I_0,I_1
      REAL*8 :: SCALEK
      real*8 :: ts,pland,tlm(lm),ple(lm+1),dp,tmsu2,tmsu3,tmsu4
      real*8, dimension(2,kaij) :: shnh_loc,shnh

      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP
      J_0 = GRID%J_STRT
      J_1 = GRID%J_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)

        k = ij_topo
        aij(i,j,k) = zatmo(i,j)*bygrav*idacc(ia_inst)

        k = ij_fland
        aij(i,j,k) = (fearth0(i,j)+flice(i,j))*idacc(ia_src)

        k = ij_netrdp
        aij(i,j,k) = aij(i,j,ij_trnfp0)+aij(i,j,ij_srnfp0)

        k = ij_grow
        aij(i,j,k) = (tsfrez(i,j,tf_last)-tsfrez(i,j,tf_day1))
     &       *idacc(ia_inst)

        k = ij_rtse
        aij(i,j,k) = aij(i,j,ij_trsup) - aij(i,j,ij_trsdn)

        k = ij_hwv
        aij(i,j,k) = aij(i,j,ij_evap)

        k = ij_pvs
        aij(i,j,k) = aij(i,j,ij_qs)*
     &       aij(i,j,ij_pres)/idacc(ia_ij(ij_pres))

        k = ij_albv
        aij(i,j,k) = aij(i,j,ij_srref)

        k = ij_albgv
        aij(i,j,k) = aij(i,j,ij_srvis)

        k = ij_albp
        aij(i,j,k) = aij(i,j,ij_srincp0)-aij(i,j,ij_srnfp0)*
     &       idacc(ia_ij(ij_srincp0))/idacc(ia_ij(ij_srnfp0))

        k = ij_albg
        aij(i,j,k) = aij(i,j,ij_srincg)-aij(i,j,ij_srnfg)*
     &       idacc(ia_ij(ij_srincg))/idacc(ia_ij(ij_srnfg))

        do k=ij_dzt1,ij_dzt1+kgz_max-2
          k1 = k-ij_dzt1+1  ; k2 = ij_phi1k + k1
          scalek = 1./(rgas*log(pmb(k1)/pmb(k1+1)))
          aij(i,j,k) = (scalek*(ght(k1+1)-ght(k1))*grav-tf)*
     &         idacc(ia_ij(ij_phi1k))
     &         +scalek*(aij(i,j,k2)-aij(i,j,k2-1))
        enddo

        k = ij_jet
        aij(i,j,k) = sqrt(aij(i,j,ij_ujet)**2+aij(i,j,ij_vjet)**2)

        k = ij_wsmn
        aij(i,j,k) = sqrt(aij(i,j,ij_us)**2+aij(i,j,ij_vs)**2)

        k = ij_wsdir
        aij(i,j,k) = atan2(aij(i,j,ij_us)+teeny,aij(i,j,ij_vs)+teeny)
     &       *idacc(ia_inst)

        k = ij_jetdir
        aij(i,j,k)=atan2(aij(i,j,ij_ujet)+teeny,aij(i,j,ij_vjet)+teeny)
     &       *idacc(ia_inst)

        k = ij_clrsky
        aij(i,j,k) = idacc(ia_ij(ij_cldcv))-aij(i,j,ij_cldcv)

        k = ij_pocean
        aij(i,j,k) = idacc(ia_ij(k))*focean(i,j)

        k = ij_pwater
        aij(i,j,k) = idacc(ia_ij(k))*focean(i,j)+aij(i,j,ij_lk)

        k = ij_tatm
        aij(i,j,k) = sum(aijl(i,j,:,ijk_tx))

        k = ij_swaerabs
        do n=1,8
        aij(i,j,k+n-1)=aij(i,j,ij_swaerrf+n-1)-aij(i,j,ij_swaersrf+n-1)
        end do

        k = ij_lwaerabs
        do n=1,8
        aij(i,j,k+n-1)=aij(i,j,ij_lwaerrf+n-1)-aij(i,j,ij_lwaersrf+n-1)
        end do

        k = ij_swaerabsnt
        aij(i,j,k) = aij(i,j,ij_swaerrfnt)-aij(i,j,ij_swaersrfnt)

        k = ij_lwaerabsnt
        aij(i,j,k) = aij(i,j,ij_lwaerrfnt)-aij(i,j,ij_lwaersrfnt)

        k = ij_TminC
        aij(i,j,k) = aij(i,j,ij_TmaxC) - aij(i,j,ij_TDcomp)

C**** Find MSU channel 2,3,4 temperatures (simple lin.comb. of Temps)
        pland = fearth0(i,j)+flice(i,j)
        ts = aij(i,j,ij_ts)/idacc(ia_ij(ij_ts))
        ple(lm+1) = pmtop
        do l=lm,1,-1
          dp = aijl(i,j,l,ijk_dp)
          ple(l) = ple(l+1)+dp/idacc(ia_dga)
          tlm(l) = ts
          if(dp.gt.0.) tlm(l)=aijl(i,j,l,ijk_tx)/dp
        enddo
        call diag_msu(pland,ts,tlm,ple,tmsu2,tmsu3,tmsu4)
        aij(i,j,ij_msu2) = tmsu2*idacc(ia_inst)
        aij(i,j,ij_msu3) = tmsu3*idacc(ia_inst)
        aij(i,j,ij_msu4) = tmsu4*idacc(ia_inst)

      ENDDO
      ENDDO

c
c fill poles
c
      if(hassouthpole(grid)) then
        do k=1,kaij
          if(jgrid_ij(k).eq.1) aij(2:im,1,k) = aij(1,1,k)
        enddo
      endif
      if(hasnorthpole(grid)) then
        do k=1,kaij
          if(jgrid_ij(k).eq.1) aij(2:im,jm,k) = aij(1,jm,k)
        enddo
      endif

c
c compute hemispheric/global means
c
      shnh_loc = 0.
      do k=1,kaij
        if(jgrid_ij(k).ne.1) cycle ! only do primary grid for now
        do j=j_0,j_1
        do i=i_0,i_1
          khem = 1
          if(lat2d(i,j).ge.0.) khem = 2
          shnh_loc(khem,k) = shnh_loc(khem,k) + axyp(i,j)*aij(i,j,k)
        enddo
        enddo
      enddo
      call sumxpe(shnh_loc,shnh)
      if(am_i_root()) then
        do k=1,kaij
          shnh(:,k) = 2.*shnh(:,k)/areag
          hemis_ij(1,1:2,k) = shnh(1:2,k)
          hemis_ij(1,3,k) = .5*(shnh(1,k)+shnh(2,k))
        enddo
      endif

      return
      end subroutine calc_derived_aij

      SUBROUTINE DIAGJ_PREP
      USE DOMAIN_DECOMP_ATM, ONLY : AM_I_ROOT
      USE CONSTANT, only : teeny
      USE MODEL_COM, only : dtsrc,idacc,ntype
      USE DIAG_COM, only : jm_budg,
     &     aj,ntype_out,wt=>wtj_comp,aj_out,areg,areg_out,
     &     nreg,kaj,j_albp0,j_srincp0,j_albg,j_srincg,
     &     j_srabs,j_srnfp0,j_srnfg,j_trnfp0,j_hsurf,j_trhdt,j_trnfp1,
     *     j_hatm,j_rnfp0,j_rnfp1,j_srnfp1,j_rhdt,j_hz1,j_prcp,j_prcpss,
     *     j_prcpmc,j_hz0,j_implh,j_shdt,j_evhdt,j_eprcp,j_erun,
     *     j_hz2,j_ervr,
     *     ia_src,ia_rad,ia_inst,
     &     sarea=>sarea_reg,
     &     hemis_j,dxyp_budg,
     &     consrv,hemis_consrv,kcon,nsum_con,scale_con,ia_con
      IMPLICIT NONE
      REAL*8 :: A1BYA2,hemfac
      INTEGER :: J,JR,J1,J2,K,M,IT

      IF(.NOT.AM_I_ROOT()) RETURN

C**** CALCULATE THE DERIVED QUANTTIES
      A1BYA2=IDACC(ia_src)/(IDACC(ia_rad)+teeny)
      DO JR=1,23
        AREG(JR,J_SRABS) =AREG(JR,J_SRNFP0)-AREG(JR,J_SRNFG)
        AREG(JR,J_RNFP0) =AREG(JR,J_SRNFP0)+AREG(JR,J_TRNFP0)
        AREG(JR,J_RNFP1) =AREG(JR,J_SRNFP1)+AREG(JR,J_TRNFP1)
        AREG(JR,J_ALBP0)=AREG(JR,J_SRINCP0)-AREG(JR,J_SRNFP0)
        AREG(JR,J_ALBG)=AREG(JR,J_SRINCG)-AREG(JR,J_SRNFG)
        AREG(JR,J_RHDT)  =A1BYA2*AREG(JR,J_SRNFG)*DTSRC+AREG(JR,J_TRHDT)
        AREG(JR,J_PRCP)  =AREG(JR,J_PRCPSS)+AREG(JR,J_PRCPMC)
        AREG(JR,J_HZ0)=AREG(JR,J_RHDT)+AREG(JR,J_SHDT)+
     *                 AREG(JR,J_EVHDT)+AREG(JR,J_EPRCP)
        AREG(JR,J_HZ1)=AREG(JR,J_HZ0)+AREG(JR,J_ERVR)
        AREG(JR,J_HZ2)=AREG(JR,J_HZ1)-AREG(JR,J_ERUN)-AREG(JR,J_IMPLH)
      END DO
      DO J=1,JM_BUDG
      DO IT=1,NTYPE
        AJ(J,J_SRABS ,IT)=AJ(J,J_SRNFP0,IT)-AJ(J,J_SRNFG,IT)
        AJ(J,J_RNFP0 ,IT)=AJ(J,J_SRNFP0,IT)+AJ(J,J_TRNFP0,IT)
        AJ(J,J_RNFP1 ,IT)=AJ(J,J_SRNFP1,IT)+AJ(J,J_TRNFP1,IT)
        AJ(J,J_ALBP0,IT)=AJ(J,J_SRINCP0,IT)-AJ(J,J_SRNFP0,IT)
        AJ(J,J_ALBG,IT)=AJ(J,J_SRINCG,IT)-AJ(J,J_SRNFG,IT)
        AJ(J,J_RHDT  ,IT)=A1BYA2*AJ(J,J_SRNFG,IT)*DTSRC+AJ(J,J_TRHDT,IT)
        AJ(J,J_PRCP  ,IT)=AJ(J,J_PRCPSS,IT)+AJ(J,J_PRCPMC,IT)
        AJ(J,J_HZ0,IT)=AJ(J,J_RHDT,IT)+AJ(J,J_SHDT,IT)+
     *                 AJ(J,J_EVHDT,IT)+AJ(J,J_EPRCP,IT)
        AJ(J,J_HZ1,IT)=AJ(J,J_HZ0,IT)+AJ(J,J_ERVR,IT)
        AJ(J,J_HZ2,IT)=AJ(J,J_HZ1,IT)-AJ(J,J_ERUN,IT)-AJ(J,J_IMPLH,IT)
      END DO
      END DO

c
c construct the composites of surface types
c
      aj_out = 0d0
      do m=1,ntype_out
        do it=1,ntype
          aj_out(:,:,m) = aj_out(:,:,m) + aj(:,:,it)*wt(m,it)
        enddo
      enddo


C**** CALCULATE SUM OF CONSRV CHANGES
C**** LOOP BACKWARDS SO THAT INITIALISATION IS DONE BEFORE SUMMATION!
      do k=kcon,1,-1
        if(nsum_con(k).eq.0) then
          consrv(:,k)=0.
        elseif(nsum_con(k).gt.0) then
          consrv(:,nsum_con(k)) = consrv(:,nsum_con(k)) +
     &         consrv(:,k)*scale_con(k)*
     &         max(1,idacc(ia_inst))/(idacc(ia_con(k))+1d-20)
        endif
      enddo

c
c compute hemispheric and global means
c
      hemfac = 2./sum(dxyp_budg)
      do m=1,ntype_out
      do k=1,kaj
        j1 = 1; j2 = jm_budg/2
        hemis_j(1,k,m) = hemfac*sum(aj_out(j1:j2,k,m)*dxyp_budg(j1:j2))
        j1 = jm_budg/2+1; j2 = jm_budg
        hemis_j(2,k,m) = hemfac*sum(aj_out(j1:j2,k,m)*dxyp_budg(j1:j2))
        hemis_j(3,k,m) = .5*(hemis_j(1,k,m)+hemis_j(2,k,m))
      enddo
      enddo
      do k=1,kcon
        j1 = 1; j2 = jm_budg/2
        hemis_consrv(1,k) =
     &       hemfac*sum(consrv(j1:j2,k)*dxyp_budg(j1:j2))
        j1 = jm_budg/2+1; j2 = jm_budg
        hemis_consrv(2,k) =
     &       hemfac*sum(consrv(j1:j2,k)*dxyp_budg(j1:j2))
        hemis_consrv(3,k) = .5*(hemis_consrv(1,k)+hemis_consrv(2,k))
      enddo

c
c scale areg by area
c
      do jr=1,nreg
        areg_out(jr,:) = areg(jr,:)/sarea(jr)
      enddo

      RETURN
      END SUBROUTINE DIAGJ_PREP

      subroutine diagjl_prep
      use model_com, only : lm,lm_req,do_gwdrag
      use domain_decomp_atm, only : am_i_root
      use diag_com, only : kajl,jm_budg,
     &     ajl,asjl,jl_srhr,jl_trcr,jl_rad_cool,
     &     jl_sumdrg,jl_dumtndrg,jl_dushrdrg,
     &     jl_mcdrgpm10,jl_dumcdrgm10,jl_dumcdrgp10,
     &     jl_mcdrgpm20,jl_dumcdrgm20,jl_dumcdrgp20,
     &     jl_mcdrgpm40,jl_dumcdrgm40,jl_dumcdrgp40,
     &     jl_dudfmdrg,jl_dudtsdif,
     &     dxyp_budg,hemis_jl,vmean_jl
      implicit none
      integer :: j,j1,j2,l,k,lr,n
      real*8 :: hemfac

      if(.not.am_i_root()) return

      do j=1,jm_budg
        do lr=1,lm_req
          asjl(j,lr,5)=asjl(j,lr,3)+asjl(j,lr,4)
        enddo
        do l=1,lm
          ajl(j,l,jl_rad_cool)=ajl(j,l,jl_srhr)+ajl(j,l,jl_trcr)
        enddo
      enddo

      if(do_gwdrag) then
        n = jl_mcdrgpm10
        ajl(:,:,n) = ajl(:,:,jl_dumcdrgm10)+ajl(:,:,jl_dumcdrgp10)
        n = jl_mcdrgpm20
        ajl(:,:,n) = ajl(:,:,jl_dumcdrgm20)+ajl(:,:,jl_dumcdrgp20)
        n = jl_mcdrgpm40
        ajl(:,:,n) = ajl(:,:,jl_dumcdrgm40)+ajl(:,:,jl_dumcdrgp40)
        n = jl_sumdrg
        ajl(:,:,n) = sum( ajl(:,:,
     &       (/
     &       jl_dumtndrg, jl_dushrdrg,
     &       jl_mcdrgpm10, jl_mcdrgpm20, jl_mcdrgpm40,
     &       jl_dudfmdrg, jl_dudtsdif
     &       /)              ), dim=3)
      endif

c
c compute hemispheric/global means and vertical sums
c
      hemfac = 2./sum(dxyp_budg)
      do k=1,kajl
        do l=1,lm
          j1 = 1; j2 = jm_budg/2
          hemis_jl(1,l,k) = hemfac*sum(ajl(j1:j2,l,k)*dxyp_budg(j1:j2))
          j1 = jm_budg/2+1; j2 = jm_budg
          hemis_jl(2,l,k) = hemfac*sum(ajl(j1:j2,l,k)*dxyp_budg(j1:j2))
          hemis_jl(3,l,k) = .5*(hemis_jl(1,l,k)+hemis_jl(2,l,k))
        enddo
        vmean_jl(1:jm_budg,1,k) = sum(ajl(:,:,k),dim=2)
        vmean_jl(jm_budg+1:jm_budg+3,1,k) = sum(hemis_jl(:,:,k),dim=2)
      enddo

      return
      end subroutine diagjl_prep

      subroutine diag_isccp_prep
c calculates the denominator array for ISCCP histograms
      use diag_com, only : aij=>aij_loc,ij_scldi,nisccp,wisccp
      use clouds_com, only : isccp_reg2d
      use domain_decomp_atm, only : grid,sumxpe
      use geom, only : axyp
      implicit none
      real*8 wisccp_loc(nisccp)
      integer :: i,j,n,j_0,j_1,i_0,i_1
      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop
      wisccp_loc(:) = 0.
      do j=j_0,j_1
      do i=i_0,i_1
        n=isccp_reg2d(i,j)
        if(n.gt.0)
     &       wisccp_loc(n)=wisccp_loc(n)+aij(i,j,ij_scldi)*axyp(i,j)
      enddo
      enddo
      call sumxpe(wisccp_loc,wisccp)
      return
      end subroutine diag_isccp_prep

      SUBROUTINE VNTRP1 (KM,P,AIN,  LMA,PE,AOUT)
C**** Vertically interpolates a 1-D array
C**** Input:       KM = number of input pressure levels
C****            P(K) = input pressure levels (mb)
C****          AIN(K) = input quantity at level P(K)
C****             LMA = number of vertical layers of output grid
C****           PE(L) = output pressure levels (mb) (edges of layers)
C**** Output: AOUT(L) = output quantity: mean between PE(L-1) & PE(L)
C****
      implicit none
      integer, intent(in) :: km,lma
      REAL*8, intent(in)  :: P(0:KM),AIN(0:KM),    PE(0:LMA)
      REAL*8, intent(out) :: AOUT(LMA)

      integer k,k1,l
      real*8 pdn,adn,pup,aup,psum,asum

C****
      PDN = PE(0)
      ADN = AIN(0)
      K=1
C**** Ignore input levels below ground level pe(0)=p(0)
      IF(P(1).GT.PE(0)) THEN
         DO K1=2,KM
         K=K1
         IF(P(K).LT.PE(0)) THEN  ! interpolate to ground level
           ADN=AIN(K)+(AIN(K-1)-AIN(K))*(PDN-P(K))/(P(K-1)-P(K))
           GO TO 300
         END IF
         END DO
         STOP 'VNTRP1 - error - should not get here'
      END IF
C**** Integrate - connecting input data by straight lines
  300 DO 330 L=1,LMA
      ASUM = 0.
      PSUM = 0.
      PUP = PE(L)
  310 IF(P(K).le.PUP)  GO TO 320
      PSUM = PSUM + (PDN-P(K))
      ASUM = ASUM + (PDN-P(K))*(ADN+AIN(K))/2.
      PDN  = P(K)
      ADN  = AIN(K)
      K=K+1
      IF(K.LE.KM) GO TO 310
      stop 'VNTRP1 - should not happen'
C****
  320 AUP  = AIN(K) + (ADN-AIN(K))*(PUP-P(K))/(PDN-P(K))
      PSUM = PSUM + (PDN-PUP)
      ASUM = ASUM + (PDN-PUP)*(ADN+AUP)/2.
      AOUT(L) = ASUM/PSUM
      PDN = PUP
  330 ADN = AUP
C****
      RETURN
      END subroutine vntrp1

      subroutine calc_derived_acc_atm
      use diag_com, only : isccp_diags
      implicit none
      call gather_zonal_diags
      call collect_scalars
      call calc_derived_aij
      call calc_derived_aijk
      call diagj_prep
      call diagjl_prep
      call diaggc_prep
      call diag_river_prep
      if(isccp_diags.eq.1) call diag_isccp_prep
      return
      end subroutine calc_derived_acc_atm
