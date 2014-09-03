#include "rundeck_opts.h"
#define JJ(J) (J)-J_0H+1

      module ATMDYN
      implicit none

      contains

      SUBROUTINE init_ATMDYN
      return
      end SUBROUTINE init_ATMDYN

      SUBROUTINE DYNAM
      USE MODEL_COM, only : im,lm,t,p,q,ls1,NSTEPSCM,I_TARG,J_TARG
      USE SCMDIAG, only : dTfrc,dqfrc,dTtot,dqtot 
      USE SOMTQ_COM, only : tmom,mz
      USE DOMAIN_DECOMP_1D, only : grid
      USE DYNAMICS, only : PMID,PEDN,PUA,PVA,SDA     

      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &     TZ,PIJL

      INTEGER L

      do L=1,LM
         dTtot(L) = T(I_TARG,J_TARG,L)
         dqtot(L) = Q(I_TARG,J_TARG,L)
         dTfrc(L) = T(I_TARG,J_TARG,L)
         dqfrc(L) = Q(I_TARG,J_TARG,L)
      enddo

      do L=1,LM
         PUA(:,:,L) = 0.
         PVA(:,:,L) = 0.
         SDA(:,:,L) = 0.
      ENDDO

      call pass_SCMDATA

      CALL CALC_PIJL(LM,P,PIJL)
      CALL CALC_AMPK(LM)

      call SCM_FORCN

      CALL tq_zmom_init(T,Q,PMID,PEDN)

      DO L=1,LM
         TZ(:,:,L)  = TMOM(MZ,:,:,L)
      ENDDO

      CALL PGF_SCM(T,TZ,PIJL)

      call FCONV
     
c
c     calculate change in T and Q due to large scale forcings
      do L=1,LM
         dTfrc(L) = T(I_TARG,J_TARG,L)-dTfrc(L)
         dqfrc(L) = Q(I_TARG,J_TARG,L)-dqfrc(L)
      enddo

      return
      END SUBROUTINE DYNAM

      SUBROUTINE SCM_FORCN
c     apply advective forcings from ARM Variational analysis to T and Q

      USE RESOLUTION , only : LM
      USE MODEL_COM , only : P,T,Q,PTOP,SIG,NSTEPSCM,DTSRC        
     &                ,I_TARG,J_TARG
      USE DYNAMICS, only : PK
      USE CONSTANT , only : KAPA 

      USE SCMCOM , only : SG_HOR_TMP_ADV, SG_VER_S_ADV, SG_HOR_Q_ADV,
     &           SG_VER_Q_ADV,SG_T,SG_Q,iu_scm_prt, 
     &           SCM_RELAX_FORCING_FLAG,SCM_DEL_T, SCM_DEL_Q

      IMPLICIT NONE

      REAL*8  dTdtls,dqdtls,tadv,dTrel,dqrel
      real*8  deltnorel,tnorel
c     tadv = advective time scale for relaxation toward observed values
c            run with  tadv = 3 hours  (in secs=10800)
      parameter (tadv=10800)

      INTEGER L

      do L = 1,LM
         T(I_TARG,J_TARG,L) = T(I_TARG,J_TARG,L)*PK(L,I_TARG,J_TARG) 
c        write(iu_scm_prt,100) NSTEPSCM,L,T(I_TARG,J_TARG,L),
c    &        SG_T(L)
ccc           ,Q(I_TARG,J_TARG,L)*1000.0,SG_Q(L)*1000.0
 100     format(1x,'FORCN before  NSTEP L T TOBS ',i5,i5,f10.3,f10.3)
      enddo

      do L = 1,LM
         deltnorel = 0.0
         Tnorel = 0.0
c        write(iu_scm_prt,*) 'tadvs ',L,SG_HOR_TMP_ADV(L),
c    *                       SG_VER_S_ADV(L)
         if (SCM_RELAX_FORCING_FLAG.eq.0) then
             SCM_DEL_T(L)=SG_HOR_TMP_ADV(L)*DTSRC+SG_VER_S_ADV(L)*DTSRC
         else
c            calculate delta T with relaxation toward observed value
c            write(iu_scm_prt,109) L,SG_T(L),T(I_TARG,J_TARG,L)
 109         format(1x,' before frc    L SGT Tmodel ',i5,2(f10.3))
             dTdtls = SG_HOR_TMP_ADV(L)+SG_VER_S_ADV(L)
             dTrel = (SG_T(L)-T(I_TARG,J_TARG,L))/tadv
             SCM_DEL_T(L) = dTdtls + dTrel
             SCM_DEL_T(L) = SCM_DEL_T(L)*DTSRC
             deltnorel = dTdtls*DTSRC
c            write(iu_scm_prt,110) dTdtls,dTrel
 110         format(1x,'wth relaxed frcing dTdtls dTrel ',2(f12.8))
         endif
         Tnorel = T(I_TARG,J_TARG,L) + deltnorel
         T(I_TARG,J_TARG,L) = T(I_TARG,J_TARG,L) + SCM_DEL_T(L)
c        write(iu_scm_prt,111) L,SCM_DEL_T(L),T(I_TARG,J_TARG,L),
c    &              deltnorel,Tnorel
 111     format(1x,' after frc L delt Tmodel ',i5,2(f10.3),
     &             ' deltnorel tnorel ',2(f10.3))
      enddo   
      do L = 1,LM
c        write(iu_scm_prt,*) 'qadvs ',L,SG_HOR_Q_ADV(L),SG_VER_Q_ADV(L)
         if (SCM_RELAX_FORCING_FLAG.eq.0) then
             SCM_DEL_Q(L)=SG_HOR_Q_ADV(L)*DTSRC+SG_VER_Q_ADV(L)*DTSRC
         else
c            calculate delta q with relaxation toward oberved value
             dqdtls = SG_HOR_Q_ADV(L)+SG_VER_Q_ADV(L)
             dqrel = (SG_Q(L)-Q(I_TARG,J_TARG,L))/tadv
             SCM_DEL_Q(L) = dqdtls + dqrel
             SCM_DEL_Q(L) = SCM_DEL_Q(L)*DTSRC
         endif
         Q(I_TARG,J_TARG,L) = Q(I_TARG,J_TARG,L) + SCM_DEL_Q(L)
         if (Q(I_TARG,J_TARG,L).lt.0.0) then
            write(99,151) NSTEPSCM,I_TARG,J_TARG,L,Q(I_TARG,J_TARG,L)
 151        format(1x,'SCM_FORCN NSTEP  I_TARG J_TARG L Q ',
     &               4(i5),f10.7)
            SCM_DEL_Q(L) = -Q(I_TARG,J_TARG,L)
            Q(I_TARG,J_TARG,L) = 0.0
         endif
      enddo
    
      do L = 1,LM
c        write(iu_scm_prt,200) NSTEPSCM,L,T(I_TARG,J_TARG,L)
cc   &               ,Q(I_TARG,J_TARG,L)*1000.0
 200     format(1x,'FORCN after NSTEP L T ',i5,i5,f10.3)
         T(I_TARG,J_TARG,L) = T(I_TARG,J_TARG,L)/PK(L,I_TARG,J_TARG)
      enddo

      RETURN

      END SUBROUTINE SCM_FORCN 

  
      SUBROUTINE FCONV
C*****
C     for single column model
C     compute CONV=Horizontal Mass Convergence
C     as filled in subroutine AFLUX in the GCM for use in
C     the CONDSE and MSTCNV Subroutines
C     Use the Wind Divergence from the ARM data
C     CONV = Wind Divergence*dSigma*P*DelArea
C
C     NOTE:    Wind Divergence is calculated for the area of the
C              ARM site. Therefore we need to take into account the
C              difference between the GCM grid box area and the ARM
C              Site.   Oklahoma site (SGP)  300 x 365 KM = 109500KM**2
C                      GCM 2 x 2.5 degrees (smaller for SGP)
C                          ~ 222.63 * 2223.42 = 49739.01 KM**2
C                     SGP/GCM = 2.2
C                     
c              Note: for NSA  domain for the variational analysis
c                    is  230KM (longitudinal) x 100KM (latitudinal)
c                          230x100 = 23000
c                    GCM 2x2.5 degree grid box ~ 21266
c               area/box = (sin(q1)-sin(q2))*2(pi)R**2/144
c                        72 degrees-71degrees
c                    ARMFAC = NSA/GCM = 23000/21266 ~ 1.08
c   
c              What about for TWP site ? ? ?
c
c * * *    now set ARMFAC in SCM_COM   so not buried here
c
c
  
      USE RESOLUTION , only : LM

      USE MODEL_COM , only : P,DSIG, I_TARG, J_TARG   
      USe GEOM , only : AXYP
   
      USE SCMCOM , only : SG_WINDIV, SG_CONV, ARMFAC, iu_scm_prt
   
      IMPLICIT NONE

      integer L,ifirst

cc    now set in SCM_COM.f   ALLOC_SCM()
cc    DATA ARMFAC/1.0/
cc    DATA ARMFAC/2.2/
cc    DATA ARMFAC/1.08/
      
      DATA ifirst/0/

c     want to fill SD (IDUM,JDUM)  check out 

c     check grid box size
      if (ifirst.eq.0) then
          write(iu_scm_prt,100) I_TARG,J_TARG,AXYP(I_TARG,J_TARG),
     &                          ARMFAC
 100      format(1x,'for I_TARG J_TARG    AXYP  ARMFAC ',i5,i5,
     &          f15.2,f10.4)
          ifirst = 1
      endif

      DO L=1,LM
         SG_CONV(L) = SG_WINDIV(L)*DSIG(L)*P(I_TARG,J_TARG)
     &                 *AXYP(I_TARG,J_TARG)*ARMFAC
      ENDDO

      return

      end SUBROUTINE FCONV  

      SUBROUTINE SDRAG(DT1)
      REAL*8, INTENT(IN) :: DT1 
      return
      END SUBROUTINE SDRAG 


      SUBROUTINE PGF_SCM (T,SZ,P)
!@SCM-version    For SCM need to calculate geopotential height. 
!                Remove other calculations.
!@sum  PGF Adds pressure gradient forces to momentum
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,kapa,bykapa,bykapap1,bykapap2
      USE MODEL_COM, only : im,jm,lm,ls1,mrch,dsig,psfmpt,sige,ptop
     *     ,zatmo,sig,modd5k,bydsig
     &     ,do_polefix,I_TARG,J_TARG   
      USE GEOM, only : imaxj,dxyv,dxv,dyv,dxyp,dyp,dxp,acor,acor2
      USE DYNAMICS, only : gz,pu,pit,phi,spa,dut,dvt
      USE DOMAIN_DECOMP_1D, Only : grid, GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      USE DOMAIN_DECOMP_1D, only : haveLatitude
      USE SCMCOM, only : iu_scm_prt
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM):: T
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO):: FD,RFDUX
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *  P, SZ

      REAL*8 PKE(LS1:LM+1)
      REAL*8 PIJ,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP,DP,P0,X
     *     ,BYDP
      REAL*8 TZBYDP,FLUX,FDNP,FDSP,RFDU,PHIDN,FACTOR
      INTEGER I,J,L,IM1,IP1,IPOLE  !@var I,J,IP1,IM1,L,IPOLE loop variab.
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)
C****
      DO L=LS1,LM+1
        PKE(L)=(PSFMPT*SIGE(L)+PTOP)**KAPA
      END DO
C****
C**** VERTICAL DIFFERENCING
C****
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=LS1,LM
      SPA(:,:,L)=0.
      END DO
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO PRIVATE(I,J,L,DP,P0,PIJ,PHIDN,TZBYDP,X,
!$OMP*             BYDP,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP)
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
        PIJ=P(I,J,1)
        PDN=PIJ+PTOP
        PKDN=PDN**KAPA
        PHIDN=ZATMO(I,J)
C**** LOOP OVER THE LAYERS
        DO L=1,LM
          PKPDN=PKDN*PDN
          PKPPDN=PKPDN*PDN
          IF(L.GE.LS1) THEN
            DP=DSIG(L)*PSFMPT
            BYDP=1./DP
            P0=SIG(L)*PSFMPT+PTOP
            TZBYDP=2.*SZ(I,J,L)*BYDP
            X=T(I,J,L)+TZBYDP*P0
            PUP=SIGE(L+1)*PSFMPT+PTOP
            PKUP=PKE(L+1)
            PKPUP=PKUP*PUP
            PKPPUP=PKPUP*PUP
          ELSE
            DP=DSIG(L)*PIJ
            BYDP=1./DP
            P0=SIG(L)*PIJ+PTOP
            TZBYDP=2.*SZ(I,J,L)*BYDP
            X=T(I,J,L)+TZBYDP*P0
            PUP=SIGE(L+1)*PIJ+PTOP
            PKUP=PUP**KAPA
            PKPUP=PKUP*PUP
            PKPPUP=PKPUP*PUP
C****   CALCULATE SPA, MASS WEIGHTED THROUGHOUT THE LAYER
            SPA(I,J,L)=RGAS*((X+TZBYDP*PTOP)*(PKPDN-PKPUP)*BYKAPAP1
     *      -X*PTOP*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *      *BYDP
          END IF
C**** CALCULATE PHI, MASS WEIGHTED THROUGHOUT THE LAYER
          PHI(I,J,L)=PHIDN+RGAS*(X*PKDN*BYKAPA-TZBYDP*PKPDN*BYKAPAP1
     *      -(X*(PKPDN-PKPUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *      *BYDP*BYKAPAP1)
C**** CALULATE PHI AT LAYER TOP (EQUAL TO BOTTOM OF NEXT LAYER)
          PHIDN=PHIDN+RGAS*(X*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPDN-PKPUP)
     *     *BYKAPAP1)
          PDN=PUP
          PKDN=PKUP
        END DO
      END DO
      END DO
!$OMP END PARALLEL DO
C**** SET POLAR VALUES FROM THOSE AT I=1
      IF (haveLatitude(grid, J=1)) THEN
        DO L=1,LM
          SPA(2:IM,1,L)=SPA(1,1,L)
          PHI(2:IM,1,L)=PHI(1,1,L)
        END DO
      END IF
      IF (haveLatitude(grid, J=JM)) THEN
        DO L=1,LM
          SPA(2:IM,JM,L)=SPA(1,JM,L)
          PHI(2:IM,JM,L)=PHI(1,JM,L)
        END DO
      END IF

!$OMP  PARALLEL DO PRIVATE(L)
      DO L=1,LM
        GZ(:,:,L)=PHI(:,:,L)
      END DO
c     do L=1,LM
c        write(iu_scm_prt,*) 'PGF_SCM  L GZ ',L,GZ(I_TARG,J_TARG,L)
c     enddo
!$OMP END PARALLEL DO
C****
C
      RETURN
      END SUBROUTINE PGF_SCM


c     SUBROUTINE AFLUX (U,V,PIJL)
c     END SUBROUTINE AFLUX


      SUBROUTINE FILTER
      return
      END SUBROUTINE FILTER

C**** Dummy routines

      SUBROUTINE COMPUTE_DYNAM_AIJ_DIAGNOSTICS( PUA,PVA,dt)
!@sum COMPUTE_DYNAM_AIJ_DIAGNOSTICS Dummy
      use DOMAIN_DECOMP_1D, only: grid

      real*8, intent(in) :: PUA(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: PVA(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: dt

      return
      END SUBROUTINE COMPUTE_DYNAM_AIJ_DIAGNOSTICS

      end module ATMDYN

      subroutine DIAG5D(M5,NDT,DUT,DVT)
c     dummy of subroutine
      USE MODEL_COM, only : im,imh,jm,lm
      USE DOMAIN_DECOMP_1D, only : GRID
      IMPLICIT NONE
      INTEGER :: M5,NDT
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        DUT,DVT
      return
      end subroutine DIAG5D

      SUBROUTINE DIAG5F(UX,VX)
c     dummy of subroutine
      USE MODEL_COM, only : im,imh,jm,lm
      USE DOMAIN_DECOMP_1D, only : GRID
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        UX,VX
      return
      end subroutine DIAG5F

      subroutine DIAGCD(grid,M,UX,VX,DUT,DVT,DT1)
c     dummy of subroutine

      USE MODEL_COM, only : im,jm,lm,fim,mdiag,mdyn
      USE DOMAIN_DECOMP_1D, only : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid
!@var M index denoting from where DIAGCD is called
      INTEGER, INTENT(IN) :: M
!@var DT1 current time step
      REAL*8, INTENT(IN) :: DT1
!@var UX,VX current velocities
      REAL*8, INTENT(IN),
     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        UX,VX
!@var DUT,DVT current momentum changes
      REAL*8, INTENT(IN),
     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        DUT,DVT

      return
      end subroutine DIAGCD

      SUBROUTINE conserv_AM(AM)
      USE DOMAIN_DECOMP_ATM, only : GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: AM
c empty - not needed yet.
      am = 0d0
      return
      end subroutine conserv_AM

      SUBROUTINE conserv_KE(RKE)
!@sum  conserv_KE calculates A-grid column-sum atmospheric kinetic energy,
!@sum  multiplied by cell area
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,dsig,ls1,p,u,v,psfmpt
      USE GEOM, only : dxyn,dxys,dxyv
      USE DOMAIN_DECOMP_1D, only : GET, CHECKSUM, HALO_UPDATE, GRID
      USE DOMAIN_DECOMP_1D, only : SOUTH
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: RKE
      INTEGER :: I,IP1,J,L
      INTEGER :: J_0STG,J_1STG
      REAL*8 :: PSJ,PSIJ

      CALL GET(grid, J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)

C****
C**** KINETIC ENERGY ON B GRID
C****

      CALL HALO_UPDATE(grid, P, FROM=SOUTH)
      DO J=J_0STG,J_1STG
      PSJ=(2.*PSFMPT*DXYV(J))
      I=IM
      DO IP1=1,IM
        PSIJ=(P(I,J-1)+P(IP1,J-1))*DXYN(J-1)+(P(I,J)+P(IP1,J))*DXYS(J)
        RKE(I,J)=0.
        DO L=1,LS1-1
          RKE(I,J)=RKE(I,J)+
     &         (U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*DSIG(L)
        END DO
        RKE(I,J)=RKE(I,J)*PSIJ
        DO L=LS1,LM
          RKE(I,J)=RKE(I,J)+
     &         (U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*DSIG(L)*PSJ
        END DO
        RKE(I,J)=0.25*RKE(I,J)*mb2kg
        I=IP1
      END DO
      END DO

c
c move to A grid
c
      call regrid_btoa_ext(rke)

      RETURN
C****
      END SUBROUTINE conserv_KE

      SUBROUTINE calc_kea_3d(kea)
!@sum  calc_kea_3d calculates square of wind speed on the A grid
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,byim,u,v
c      USE GEOM, only : ravps,ravpn
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, GRID
      USE DOMAIN_DECOMP_1D, only : NORTH
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: KEA

      RETURN

      END SUBROUTINE calc_kea_3d

      subroutine recalc_agrid_uv
!@sum Computes u_a,v_a from u and v
!@var u x-component at secondary grids (B_grid)
!@var v y-component at secondary grids (B_grid)
!@var u_a x-component at primary grids (A_grid)
!@var v_a y-component at primary grids (A_grid)
!@auth Ye Cheng
!@ver  1.0

      USE MODEL_COM, only : im,jm,lm,u,v
      USE DYNAMICS, only : ua=>ualij,va=>valij
      USE DOMAIN_DECOMP_1D, only : grid,get,NORTH, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP_1D, only : halo_update
      USE GEOM, only : imaxj,idij,idjj,kmaxj,rapj,cosiv,siniv
      implicit none
      real*8, dimension(im) :: ra
      integer, dimension(im) :: idj
      real*8 :: HEMI,u_t,v_t,rak,ck,sk,uk,vk
      integer :: i,j,l,k,idik,idjk,kmax

      integer :: J_0S,J_1S
      logical :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      call get(grid, J_STRT_SKP=J_0S,   J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE    )
!     polar boxes

C**** Update halos of U and V
      CALL HALO_UPDATE(grid,u, from=NORTH)
      CALL HALO_UPDATE(grid,v, from=NORTH)

      if (HAVE_SOUTH_POLE) then
        J=1
        KMAX=KMAXJ(J)
        HEMI=-1.
!$OMP  PARALLEL DO PRIVATE (I,L,u_t,v_t,K,IDIK,IDJK,RAK,ck,sk,uk,vk)
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              RAK=RAPJ(K,J)
              ck=cosiv(k)
              sk=siniv(k)
              uk=u(idik,idjk,L)
              vk=v(idik,idjk,L)
              u_t=u_t+rak*(uk*ck-hemi*vk*sk)
              v_t=v_t+rak*(vk*ck+hemi*uk*sk)
            END DO
            ua(l,i,j)=u_t
            va(l,i,j)=v_t
          END DO
        END DO
!$OMP  END PARALLEL DO
      end if              !south pole
!
      if (HAVE_NORTH_POLE) then
        J=JM
        KMAX=KMAXJ(J)
        HEMI=1.
!$OMP  PARALLEL DO PRIVATE (I,L,u_t,v_t,K,IDIK,IDJK,RAK,ck,sk,uk,vk)
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              RAK=RAPJ(K,J)
              ck=cosiv(k)
              sk=siniv(k)
              uk=u(idik,idjk,L)
              vk=v(idik,idjk,L)
              u_t=u_t+rak*(uk*ck-hemi*vk*sk)
              v_t=v_t+rak*(vk*ck+hemi*uk*sk)
            END DO
            ua(l,i,j)=u_t
            va(l,i,j)=v_t
          END DO
        END DO
!$OMP  END PARALLEL DO
      end if                !north pole

!     non polar boxes
C**** Update halos of u and v. (Needed bcs. IDJJ(3:4,J_1S)=J_1S+1)
C     ---> done by calling routine...

c      CALL HALO_UPDATE(grid, u, FROM=NORTH)
c      CALL HALO_UPDATE(grid, v, FROM=NORTH)
!$OMP  PARALLEL DO PRIVATE (J,I,L,u_t,v_t,K,KMAX,IDJ,RA,IDIK,IDJK,RAK)
!$OMP*    SCHEDULE(DYNAMIC,2)
      DO J=J_0S,J_1S
        KMAX=KMAXJ(J)
        DO K=1,KMAX
          IDJ(K)=IDJJ(K,J)
          RA(K)=RAPJ(K,J)
        END DO
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJ(K)
              RAK=RA(K)
              u_t=u_t+u(IDIK,IDJK,L)*RAK
              v_t=v_t+v(IDIK,IDJK,L)*RAK
            END DO
            ua(l,i,j)=u_t
            va(l,i,j)=v_t
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO
C****
      return
      end subroutine recalc_agrid_uv

      subroutine regrid_atov_1d(u_a,v_a,uv1d)
      USE MODEL_COM, only : im,jm
      USE DOMAIN_DECOMP_1D, only : grid,halo_update,SOUTH
      USE GEOM, only : rapvs,rapvn,cosiv,siniv
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo)  ::
     &          u_a,v_a
      real*8, dimension(2*im*(1+grid%j_stop_stgr-grid%j_strt_stgr)),
     &        intent(out) :: uv1d
      real*8 :: hemi
      integer :: i,ip1,j,n
      real*8, dimension(im) :: usouth,vsouth,unorth,vnorth
      integer :: j_0stg, j_1stg
      j_0stg = grid%j_strt_stgr
      j_1stg = grid%j_stop_stgr
      CALL HALO_UPDATE(grid,U_A,from=SOUTH)
      CALL HALO_UPDATE(grid,V_A,from=SOUTH)
      j=j_0stg-1
      if(j.eq.1) then
        hemi = -1.
        usouth(:)=2.*(u_a(1,j)*cosiv(:)+v_a(1,j)*siniv(:)*hemi)
        vsouth(:)=2.*(v_a(1,j)*cosiv(:)-u_a(1,j)*siniv(:)*hemi)
      else
        i=im
        do ip1=1,im
          usouth(i)=(u_a(i,j)+u_a(ip1,j))
          vsouth(i)=(v_a(i,j)+v_a(ip1,j))
          i=ip1
        enddo
      endif
      n = 0
      do j=j_0stg,j_1stg
        if(j.lt.jm) then
          i=im
          do ip1=1,im
            unorth(i)=(u_a(i,j)+u_a(ip1,j))
            vnorth(i)=(v_a(i,j)+v_a(ip1,j))
            i=ip1
          enddo
        else
          hemi = +1.
          unorth(:)=2.*(u_a(1,j)*cosiv(:)+v_a(1,j)*siniv(:)*hemi)
          vnorth(:)=2.*(v_a(1,j)*cosiv(:)-u_a(1,j)*siniv(:)*hemi)
        endif
        do i=1,im
          n = n + 1
          uv1d(n) = rapvn(j-1)*usouth(i)+rapvs(j)*unorth(i)
          n = n + 1
          uv1d(n) = rapvn(j-1)*vsouth(i)+rapvs(j)*vnorth(i)
          usouth(i) = unorth(i)
          vsouth(i) = vnorth(i)
        enddo
      enddo
      return
      end subroutine regrid_atov_1d
      subroutine get_nuv(nuv)
      use model_com, only : im
      USE DOMAIN_DECOMP_1D, only : GRID
      implicit none
      integer :: nuv
      nuv = 2*im*(1+grid%j_stop_stgr-grid%j_strt_stgr)
      return
      end subroutine get_nuv
      subroutine get_vpkey_of_n(n,vpkey)
      implicit none
      integer :: n,vpkey
      vpkey = 1+(n-1)/2
      return
      end subroutine get_vpkey_of_n
      subroutine get_regrid_info_for_n(n,ilist,jlist,wts,nnbr)
      use model_com, only : im
      use geom, only : rapvn,rapvs
      implicit none
      integer :: n
      integer, dimension(4) :: ilist,jlist
      real*8, dimension(4) :: wts
      integer :: nnbr
      integer :: iv,jv,ivp1
      call get_ivjv_of_n(n,iv,jv)
      nnbr = 4
      ivp1 = iv+1 - im*(iv/im)
      ilist(1:4) = (/ iv, ivp1, iv, ivp1 /)
      jlist(1:4) = (/ jv-1, jv-1, jv, jv /)
      wts(1:4) = (/ rapvn(jv-1), rapvn(jv-1), rapvs(jv), rapvs(jv) /)
      return
      end subroutine get_regrid_info_for_n
      subroutine get_uv_of_n(n,uv)
      use model_com, only : im,lm,u,v
      use domain_decomp_1d, only : am_i_root
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      integer :: iv,jv
      call get_ivjv_of_n(n,iv,jv)
      if(mod(n,2).eq.1) then
        uv(1:lm) = u(iv,jv,1:lm)
      else
        uv(1:lm) = v(iv,jv,1:lm)
      endif
      return
      end subroutine get_uv_of_n
      subroutine store_uv_of_n(n,uv)
      use model_com, only : im,lm,u,v
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      integer :: iv,jv
      call get_ivjv_of_n(n,iv,jv)
      if(mod(n,2).eq.1) then
        u(iv,jv,1:lm) = uv(1:lm)
      else
        v(iv,jv,1:lm) = uv(1:lm)
      endif
      return
      end subroutine store_uv_of_n
      subroutine get_ivjv_of_n(n,iv,jv)
      use model_com, only : im
      USE DOMAIN_DECOMP_1D, only : GRID
      implicit none
      integer :: n
      integer :: iv,jv
      integer :: nv,njm1
      nv = 1+(n-1)/2
      njm1 = (nv-1)/im
      jv = grid%j_strt_stgr + njm1
      iv = nv - njm1*im
      return
      end subroutine get_ivjv_of_n

      subroutine replicate_uv_to_agrid(ur,vr,k,ursp,vrsp,urnp,vrnp)
      USE MODEL_COM, only : im,jm,lm,u,v
      USE DOMAIN_DECOMP_1D, only : GRID, hasNorthPole, hasSouthPole
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     UR,VR
      real*8, dimension(im,lm) :: ursp,vrsp,urnp,vrnp
      integer :: i,j,l
      integer :: J_0S,J_1S
      if(k.ne.4)
     &     call stop_model('incorrect k in replicate_uv_to_agrid',255)
      J_0S = GRID%J_STRT_SKP
      J_1S = GRID%J_STOP_SKP
#ifdef SCM
      do j=j_0s,j_1s
      do i=1,im
      do l=1,lm
        ur(1,l,i,j) = u(i,j,l)
        vr(1,l,i,j) = v(i,j,l)
        ur(2,l,i,j) = u(i,j,l)
        vr(2,l,i,j) = v(i,j,l)
        ur(3,l,i,j) = u(i,j,l)
        vr(3,l,i,j) = v(i,j,l)
        ur(4,l,i,j) = u(i,j,l)
        vr(4,l,i,j) = v(i,j,l)
      enddo ! l
      enddo ! i
      enddo ! j
#else
      do j=j_0s,j_1s
      do i=2,im
      do l=1,lm
        ur(1,l,i,j) = u(i-1,j  ,l)
        vr(1,l,i,j) = v(i-1,j  ,l)
        ur(2,l,i,j) = u(i  ,j  ,l)
        vr(2,l,i,j) = v(i  ,j  ,l)
        ur(3,l,i,j) = u(i-1,j+1,l)
        vr(3,l,i,j) = v(i-1,j+1,l)
        ur(4,l,i,j) = u(i  ,j+1,l)
        vr(4,l,i,j) = v(i  ,j+1,l)
      enddo ! l
      enddo ! i
      i = 1
      do l=1,lm
        ur(1,l,i,j) = u(im ,j  ,l)
        vr(1,l,i,j) = v(im ,j  ,l)
        ur(2,l,i,j) = u(i  ,j  ,l)
        vr(2,l,i,j) = v(i  ,j  ,l)
        ur(3,l,i,j) = u(im ,j+1,l)
        vr(3,l,i,j) = v(im ,j+1,l)
        ur(4,l,i,j) = u(i  ,j+1,l)
        vr(4,l,i,j) = v(i  ,j+1,l)
      enddo ! l
      enddo ! j
#endif
      if(hasSouthPole(grid)) then
        ursp(:,:) = u(:,2,:)
        vrsp(:,:) = v(:,2,:)
      endif
      if(hasNorthPole(grid)) then
        urnp(:,:) = u(:,jm,:)
        vrnp(:,:) = v(:,jm,:)
      endif
      return
      end subroutine replicate_uv_to_agrid

      subroutine avg_replicated_duv_to_vgrid(du,dv,k,
     &     dusp,dvsp,dunp,dvnp)
      USE MODEL_COM, only : im,jm,lm,u,v
      USE DOMAIN_DECOMP_1D, only : GRID, HALO_UPDATE_BLOCK,SOUTH
      use Domain_Decomp_1d, only: hasNorthPole, hasSouthPole
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     DU,DV
      real*8, dimension(im,lm) :: dusp,dvsp,dunp,dvnp
      integer :: i,j,l
      integer :: J_0STG,J_1STG
      if(k.ne.4) call stop_model(
     &     'incorrect k in avg_replicated_duv_to_vgrid',255)
      J_0STG = GRID%J_STRT_STGR
      J_1STG = GRID%J_STOP_STGR
      CALL HALO_UPDATE_BLOCK(GRID,DU,FROM=SOUTH)
      CALL HALO_UPDATE_BLOCK(GRID,DV,FROM=SOUTH)
c
c copy circumpolar data into the appropriate spots in du,dv
c
      if(hasSouthPole(grid)) then
        j=1
        do i=2,im
        do l=1,lm
          du(3,l,i,j) = dusp(i-1,l)
          du(4,l,i,j) = dusp(i  ,l)
          dv(3,l,i,j) = dvsp(i-1,l)
          dv(4,l,i,j) = dvsp(i  ,l)
        enddo
        enddo
        i=1
        do l=1,lm
          du(3,l,i,j) = dusp(im ,l)
          du(4,l,i,j) = dusp(i  ,l)
          dv(3,l,i,j) = dvsp(im ,l)
          dv(4,l,i,j) = dvsp(i  ,l)
        enddo
c compensate for the factor of 2 in ravj(1).  change ravj(1) later.
        du(:,:,:,j) = du(:,:,:,j)*.5
        dv(:,:,:,j) = dv(:,:,:,j)*.5
      endif
      if(hasNorthPole(grid)) then
        j=jm
        do i=2,im
        do l=1,lm
          du(1,l,i,j) = dunp(i-1,l)
          du(2,l,i,j) = dunp(i  ,l)
          dv(1,l,i,j) = dvnp(i-1,l)
          dv(2,l,i,j) = dvnp(i  ,l)
        enddo
        enddo
        i=1
        do l=1,lm
          du(1,l,i,j) = dunp(im ,l)
          du(2,l,i,j) = dunp(i  ,l)
          dv(1,l,i,j) = dvnp(im ,l)
          dv(2,l,i,j) = dvnp(i  ,l)
        enddo
c compensate for the factor of 2 in ravj(jm).  change ravj(jm) later.
        du(:,:,:,j) = du(:,:,:,j)*.5
        dv(:,:,:,j) = dv(:,:,:,j)*.5
      endif
c
c now do the averaging
c
      do j=j_0stg,j_1stg
      do i=1,im-1
      do l=1,lm
        u(i,j,l)=u(i,j,l)+ !.25*
     &       (du(4,l,i,j-1)+du(3,l,i+1,j-1)+du(2,l,i,j)+du(1,l,i+1,j))
        v(i,j,l)=v(i,j,l)+ !.25*
     &       (dv(4,l,i,j-1)+dv(3,l,i+1,j-1)+dv(2,l,i,j)+dv(1,l,i+1,j))
      enddo ! l
      enddo ! i
      i = im
      do l=1,lm
        u(i,j,l)=u(i,j,l)+ !.25*
     &       (du(4,l,i,j-1)+du(3,l,1,j-1)+du(2,l,i,j)+du(1,l,1,j))
        v(i,j,l)=v(i,j,l)+ !.25*
     &       (dv(4,l,i,j-1)+dv(3,l,1,j-1)+dv(2,l,i,j)+dv(1,l,1,j))
      enddo ! l
      enddo ! j
      return
      end subroutine avg_replicated_duv_to_vgrid

      SUBROUTINE regrid_btoa_3d(x)
      USE MODEL_COM, only : im,jm,lm,byim
c      USE GEOM, only : ravps,ravpn
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, GRID
      USE DOMAIN_DECOMP_1D, only : NORTH, hasNorthPole, hasSouthPole
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: X
      INTEGER :: I,IM1,J,L
      REAL*8 :: XIM1J,XIJ
      call halo_update(grid,x,from=north)
      DO L=1,LM
      if(hasSouthPole(grid)) then
        x(:,1,l) = sum(x(:,2,l))*byim
      endif
      DO J=GRID%J_STRT_SKP,GRID%J_STOP_SKP
      IM1=IM
      XIM1J = x(im1,j,l)
      DO I=1,IM
        XIJ = x(i,j,l)
        X(I,J,L)=.25*(XIM1J+XIJ+
     &       x(im1,j+1,l)+x(i,j+1,l))
c        X(I,J,L)=(
c     &       ravps(j)*(x(im1,j,l)+x(i,j,l))
c     &      +ravpn(j)*(x(im1,j+1,l)+x(i,j+1,l))
        XIM1J = XIJ
        IM1=I
      ENDDO
      ENDDO
      if(hasNorthPole(grid)) then
        x(:,jm,l) = sum(x(:,jm,l))*byim
      endif
      ENDDO
      RETURN
      END SUBROUTINE regrid_btoa_3d

      subroutine regrid_btoa_ext(x)
c regrids scalar x_bgrid*dxyv -> x_agrid*dxyp
      USE MODEL_COM, only : im,jm,byim
      USE GEOM, only : rapvs,rapvn,dxyp,dxyv
      USE DOMAIN_DECOMP_1D, only : GET, HALO_UPDATE, GRID, NORTH
      use Domain_Decomp_1d, only: hasNorthPole, hasSouthPole
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: X
      INTEGER :: I,IM1,J
      INTEGER :: J_0S,J_1S
      REAL*8 :: XIM1J,XIJ
      CALL GET(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)
      call halo_update(grid,x,from=north)
      if(hasSouthPole(grid)) then
        X(:,1) = SUM(X(:,2))*BYIM*(DXYP(1)/DXYV(2))
      endif
      DO J=J_0S,J_1S
      IM1=IM
      XIM1J = X(IM1,J)
      DO I=1,IM
        XIJ = X(I,J)
c        X(I,J) = .25*(XIM1J+X(I,J)+X(IM1,J+1)+X(I,J+1))
        X(I,J) = (
     &       (XIM1J+X(I,J))*RAPVS(J)
     &      +(X(IM1,J+1)+X(I,J+1))*RAPVN(J) )
        XIM1J = XIJ
        IM1 = I
      ENDDO
      ENDDO
      if(hasNorthPole(grid)) then
        X(:,JM) = SUM(X(:,JM))*BYIM*(DXYP(JM)/DXYV(JM))
      endif
      return
      end subroutine regrid_btoa_ext

      SUBROUTINE QDYNAM
      return
      END SUBROUTINE QDYNAM
