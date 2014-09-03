!     SCM_DIAG_COM.f    
!     save diagnostics for run of MODELE SCM

      MODULE SCMDIAG  

      USE RESOLUTION , only : LM 
      USE SCMCOM, only : NRESET
      USE DIAG_COM , only : npres,ntau
      USE CLOUDS_COM, only : ncol

      IMPLICIT NONE

      real*8   CLCVSS(LM),CLCVMC(LM),CLTHCK(LM),CSIZE(LM,2),   
     *         EFFRAD(LM),TAUSSC(LM),TAUMCC(LM),CUMFLX(LM,2,LM),
     *         CUMHET(LM),CUMOST(LM),PRCSS,PRCMC,EVPFLX,SHFLX,
     *         SOILMS,CLDFLG(LM),DWNFLX(LM,2,LM),RHC(LM)
      real*8   SG_RH(LM),SG_CLDWAT(LM),SG_CLDICE(LM)
      real*8   clsav(LM)
      real*8 SRDFLBTOP,SRNFLBTOP,SRUFLBTOP,SRUFLBBOT,
     *       TRUFLBTOP,TRDFLBTOP,SRDFLBBOT,SRNFLBBOT,
     *       TRUFLBBOT,TRDFLBBOT,CSSRNTOP,CSTRUTOP,
     *       CSSRNBOT,CSTRNBOT,CSSRDBOT,TRNFLBBOT
      real*8 SCM_PBL_HGT
      real*8 SRFHRLCOL(LM),TRFCRLCOL(LM)
c     
c     isccp diagnostics
c
      real*8 isccp_boxptop(ncol),isccp_boxtau(ncol)
      real*8 isccp_ctp,isccp_tauopt,isccp_lowcld,isccp_midcld,
     &       isccp_highcld,isccp_fq(ntau,npres)
      integer isccp_sunlit
 
      real*8 LHPSAV(LM),LHPMC(LM),PRESAV(LM),PREMC(LM)
c
c     cumulus updraft speed
      real*8   WCUSCM(LM,2),WCUALL(LM,2,LM)
      real*8   WCUDEEP(LM,2)
C--- Added by J.W. starting ---C
      real*8   MPLUMESCM(LM,2),MPLUMEALL(LM,2,LM)
      real*8   MPLUMEDEEP(LM,2)
      real*8   ENTSCM(LM,2),ENTALL(LM,2,LM)
      real*8   ENTDEEP(LM,2)
      real*8   DETRAINDEEP(LM,2,LM)
C--- Added by J.W. ending ---C
c
c     precipitating and non-precipitating convective condensate for
c     deep convection 
      real*8 PRCCDEEP(LM,2,LM),NPRCCDEEP(LM,2,LM),TPALL(LM,2,LM)   
      real*8 PRCCGRP(LM,2,LM),PRCCICE(LM,2,LM)
     
c     condensate for all convection
      real*8 mccond(LM,2,LM)
c
      real*8 SCM_LWP_MC,SCM_LWP_SS,SCM_IWP_MC,SCM_IWP_SS,
     &       SCM_WM_MC(LM),SCM_SVWMXL(LM)

c     turbulent momentum fluxes,effective cloud albedo at surface
      real*8 scm_tmf_u,scm_tmf_v,scm_ECAS
      real*8 scm_Ssrf,scm_Hsrf
      real*8 scm_CAPE,scm_CIN
      real*8 scm_FDIF_vis,scm_FDIF_nir,scm_FDIR_vis,scm_FDIR_nir
c
c     dt and dq over time step from different processes
      real*8 dTtot(LM),dqtot(LM),dTfrc(LM),dqfrc(LM),dTrad(LM),
     &       dTradlw(LM),dTradsw(LM)
      real*8 dTHmc(LM),dqmc(LM),dTHbl(LM),dqbl(LM),dTHss(LM),dqss(LM)

c
c     convectice plume diagnostics
      real*8 PLUME_MAX(2,LM), PLUME_MIN(2,LM)
c     moist static energy
      real*8 SCM_GZ(LM)
      real*8 SCM_H(LM),SCM_HSAT(LM)
      real*8 ARM_H(LM),ARM_HSAT(LM)



      integer BBNSTEPSCM(NRESET),BBISVDATE(NRESET)

      END MODULE SCMDIAG  
