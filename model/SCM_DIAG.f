c     SCM_DIAG.f
c     save diagnostics for run of MODELE SCM

      SUBROUTINE  SCM_DIAG  


      USE RESOLUTION, only : LM
      USE MODEL_COM , only :  p,u,v,t,q,wm,NSTEPSCM,sige,sig,
     &                        I_TARG,J_TARG,dtsrc,PTOP,NSTEP
      USE CLOUDS_COM, only : SVLHX,SVLAT,RHSAV,CLDSAV,tauss,taumc,
     &                cldss,cldmc,csizmc,csizss,ncol
      use DIAG_COM, only : npres,ntau,isccp_press,isccp_tau
      USE SCMCOM
      USE SCMDIAG
      USE RAD_COM, only : srhr,trhr
      USE PBLCOM, only : TSAVG,WSAVG,QSAVG,USAVG,VSAVG
      USE FLUXES, only : GTEMP,UFLUX1,VFLUX1 
C     USE DYNAMICS, only : PK 
C--- Added by J.W. starting ---C
      USE DYNAMICS, only : PK,GZ
C--- Added by J.W. ending ---C
      USE CONSTANT, only : SHA, GRAV, kapa, RGAS, LHE 
      USE GEOM, only : axyp 
      USE FILEMANAGER, only : openunit,closeunit
      

      IMPLICIT NONE


C     ALERT!!! LX=57 is for a 53-layer model (i.e., LX=LM+4)
C              LX=40          23-layer model 
C     If LX is changed, need to change data statements in 
C     BLOCK DATA RADPAR that initialize PLB and HLB
CCC   PARAMETER (LX=57,LMOX=12*(1998-1881),MO3X=12*(2050-1950) )
c     PARAMETER (LX=40,LMOX=12*(1998-1881),MO3X=12*(2050-1850+1) )
c     REAL*4 TROAER,VDBCSU,TDUST,EPLMHC,UVLEAN   ! ,FVEG11,FOLGIZ     
c     REAL*4 ATAU98,SIZE98,HTF498,O3CLIM         !  only for offline use
c     CHARACTER*80 ATITLE,VTITLE,DTITLE,TITLE
C
c     COMMON/RADCM1/V6ALB(11,4,7),TAUWC0,TAUIC0,EPSCON,RO3COL,FULGAS(13)
c    A             ,FRAYLE,FALGAE,FCLDTR,FCLDSR,PTLISO,TLGRAD,FGOLDH(13)
c    B             ,WETTRA,WETSRA,FSXAER(5),FTXAER(5),FZSRA(6),FEMTRA(6)
c    C             ,KWVCON,KEEPAL,KEEPRH,KEEP10,KCNORM,KCLDEP,ICE012,NV
c    D             ,MADVES,MRELAY,MOZONE,KO3LON,NO3COL,NORMS0,KSOLAR
c    E             ,KTREND,NTRACE,ITR(8),NL,NLP,MLAT46,MLON72,LASTVC
c    X             ,HLB(LX),RHL(LX),TRACER(LX,8),AERTAU(LX),ZOICE
c    X             ,S00WM2,RATLS0,S0,XXJDAY,JYEARR,JDAYR,ISPARE(98)
c    X             ,DMOICE,DMLICE,TRSLCR,TRDFSL,TRUFSL
C     INPUT DATA
c     COMMON/RADCM2/
c    F              PLB(LX),        TLB(LX),TLT(LX),TLM(LX),U0GAS(LX,12)
c    G             ,ULGAS(LX,12),SHL(LX)
c    H             ,TAUWC(LX),TAUIC(LX),SIZEWC(LX),SIZEIC(LX),CLDEPS(LX)
c    I          ,POCEAN,PEARTH,POICE,PLICE,AGESN(3),SNOWE,SNOWOI,SNOWLI
c    J             ,TGO,TGE,TGOI,TGLI,TSL,WMAG,WEARTH,      FSPARE(998)
c    K                              ,COSZ,PVT(11),BXA(153),SRBXAL(15,2)
c    L                              ,JLAT,ILON
C     OUTPUT DATA
c    M             ,TRDFLB(LX),TRUFLB(LX),TRNFLB(LX),TRFCRL(LX)
c    N             ,SRDFLB(LX),SRUFLB(LX),SRNFLB(LX),SRFHRL(LX),SRSLHR
c    O             ,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR,SRXATM(4)
c    P             ,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR,ALBNIR,FSRNFG(4)
c    Q             ,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR,FTRUFG(4)
c    R             ,TRDFGW,TRUFGW,TRUFTW,BTEMPW,              DTRUFG(4)
c    S             ,TRSLTS,TRSLTG,TRSLWV,TRSLBS,TTRUFG,LBOTCL,LTOPCL
c    X             ,Z0(12)
c
c-------------------------------------------------------------------------------

C
C 
C             Record Layout  
C 
C             NSTEPSCM            Time Stamp - note --- add date/time stamp also    
C             P                   Column Pressure (mb)
C             SG_P       (LM)     Pressure at sigma layers (mb)
C             SG_HGT     (LM)     Height at sigma layers (m)
C             T          (LM)     Temperature (K) 
C             Q          (LM)     Specific Humidity (Kg/Kg)
C             SG_RH      (LM)     Relative Humidity (wrt water)
C             WMCOL      (LM)     Cloud Water Content (Kg/Kg)
C             SVLHXCOL   (LM)     Liquid/Ice Flag (SS) save Latent Heats (j/Kg)
C             SVLATCOL   (LM)     Liquid/Ice Flag (MC) save Latent Heats (j/Kg)
C             SG_CLDWAT  (LM)     Cloud Water Mixing Ratio
C             SG_CLDICE  (LM)     Cloud Ice Mixing Ratio
C             CLCVSS     (LM)     Cloud Cover SS (by area) 
C             CLCVMC     (LM)     Cloud Cover MC 
C             SG_U       (LM)     ARM u-wind (m/s)
C             SG_V       (LM)     ARM v-wind (m/s)
C             CSIZE      (LM,2)   Particle Size (10**-06m)     1=mc,2=ss 
C             EFFRAD     (LM)     Effective Radius (10**-06m)
C             CUMFLXCOL  (LM)     Cumulus Mass Flux (kg/m**2 /s) 
C             DWNFLXCOL  (LM)     Downdraft Cloud Mass Flux (kg/m**2 /s)
C             TAUSSC     (LM)     Cloud Optical Thickness - SS 
C             TAUMCC     (LM)     Cloud Optical Thickness - MC 
C             SG_T       (LM)     Arm T interpolated to Sig levels (K)
C             SG_Q       (LM)     ARM Q interpolated to Sig levels (g/KG)
C             SG_OMEGA   (LM)     ARM Omega (mb/hour)
C             SG_HOR_TMP_ADV(LM)  ARM Horizontal Temperature Advection(K/s) 
C             SG_VER_S_ADV(LM)    ARM Vertical S Advection (K/s)
C             SG_HOR_Q_ADV(LM)    ARM Horizontal Q Advection (kg/kg/s) 
C             SG_VER_Q_ADV(LM)    ARM Vertiacl Q Advection (kg/Kg/s)
C             SG_ARSCL(LM)        ARM  ARSCL CLD AMOUNT (%)
C             dTtot(LM)           dT/dt over time step (K/day)
C             dqtot(LM)           dq/dt over time step (kg/kg /day)
C             dTfrc(LM)           dT/dt over time step from FORCN (K/day)
C             dqfrc(LM)           dq/dt over time step from FORCN (kg/kg /day)
C             dTrad(LM)           dT/dt over time step from radiation (K/day)
C             dTHmc(LM)           dT/dt over time step from mstcnv (K/day)
C             dqmc(LM)            dq/dt over time step from mstcnv (kg/kg/day)
C             dTHbl(LM)           dT/dt over time step from boundary layer (srf flxs + aturb) (K/day)
C             dqbl(LM)            dq/dt over time step from boundary layer (srf flxs + aturb) (kg/kg/day)
C             dTHss(LM)           dT/dt over time step from large scale clouds (K/day) 
C             dqss(LM)            dq/dt over time step from large scale clouds (kg/kg/day) 
C             dTradlw(LM)         dT/dt over time step from radiation (lw) (K/day)
C             dTradsw(LM)         dT/dt over time step from radiation (sw) (K/day)
C             GZPRT(LM)           geopotential height
C             TSKIN               GTEMP(1,4,itarg,jtarg) - Skin temperature (C)
C             TSURF               Ts  -  Surface Air T (K)  (TSAVG)
C             QSURF               near surface   specific humidity (kg/kg)
C             USURF               surface winds  x direction (m/s)
C             VSURF               surface winds  y direction (m/s)
C             EVPFLX              Evaporation Flux 
C             SHFLX               Sensible Heat Flux 
C             PRCSS               Precip - Large Scale SS (mm/hour) 
C             PRCMC               Precip - Convective (mm/hour)
C             SCM_PBL_HGT         height of the top of the boundary layer  (m)
C             scm_tmf_u           turbulent momentum flux x-direction
C             scm_tmf_v           turbulent momentum flux y-direction
C             scm_ECAS            Effective Cloud Albedo at surface
C             scm_Ssrf            near surface Dry Static Energy (kJ/kg)
C             scm_Hsrf            near surface Moist Static Energy (kJ/kg)
C             scm_CAPE            Convective Available Potential Energy
C             scm_CIN             Convective Inhibition
C             scm_FDIR_vis        Direct vis solar radiative flux (W/m**2)
C             scm_FDIF_vis        Diffuse vis solar radiative flux (W/m**2)
C             scm_FDIR_nir        Direct near-ir solar radiative flux (W/m**2)
C             scm_FDIF_nir        Diffuse near-ir solar radiative flux (W/m**2)
C             SRDFLBTOP           INC SW on Top of Atmos (W/m**2) 
C             SRNFLBTOP           SW ABSORBTION BELOW P0 (W/m**2)
C             TRUFLBTOP           UPward LW at P0  (W/m**2)
C             SRDFLBBOT           SW INC on Z0  (W/m**2)
C             SRNFLBBOT           SW ABS on Z0  (W/m**2)
C             TRUFLBBOT           Upward LW at Z0  (W/m**2)
C             TRDFLBBOT           LW INC on Z0  (W/m**2)
C             SRUFLBBOT           Short Wave radiation up at z0 (W/m2)
C             SRUFLBTOP           Short Wave radiation up at p0 (W/m2)
C             TRDFLBTOP           Long Wave radiation down at p0 (W/m2)
C             TRNFLBBOT           LW surface net radiation (W/m2)
C             CSSRNTOP            clear sky toa sw net downward rad (W/m2)
C             CSTRUTOP            clear sky toa lw up rad (W/m2)
C             CSSRNBOT            clear sky srf sw net downward rad (W/m2)
C             CSTRNBOT            clear sky srf lw net upward rad(W/m2)
C             CSSRDBOT            clear sky srf sw downward rad (W/m2)
C             SRFHRLCOL(1-LM)*COSZ1  SW Heating 
C             TRFCRLCOL (1-LM)    LW Heating 
C             APREC               ARM precipitation (mm/hour)
C             ATSAIR              ARM Surface Temperature Air (C)
C             ATSKIN              ARM Surface Skin Temperature (C)
C             ARHSAIR             ARM Surface Relative Humidity (%)
C             ALH                 ARM Latent heat upwards (W/m2)
C             ASH                 ARM Sensible heat upwards (W/m2)
C             AMEANPS             ARM Area mean Suface Pressure (mb)
C
C             isccp record layout
C
C             isccp_sunlit        ISCCP   1-day 0-night
C             isccp_ctp           ISCCP   cloud top pressure
C             isccp_tauopt        ISCCP   optical thickness
C             isccp_lowcld        ISCCP   low cloud fraction
C             isccp_midcld        ISCCP   mid cloud fraction
C             isccp_highcld       ISCCP   high cloud fraction
C             isccp_fq(ntau,npres)ISCCP  fraction of model grid box
C                                     covered by each of the 49 ISCCP D level cloud
C                                     types
C             isccp_boxtau(ncol)  ISCCP optical thickness in each column
C             isccp_boxptop(ncol) ISCCP cloud top pressure (mb) in each column
C            
C              
C
CC    some diags not currently written out
C             PRESAV     (LM)     SCM Large Scale Precip by layer (kg/kg)
C             LHPSAV     (LM)     SCM Phase of Precip for SS
C             PREMC      (LM)     SCM MSTCNV Precip by layer (kg/kg)
C             LHPMC      (LM)     SCM Phase of Precip for MC
C             CLTHCK     (LM)     Cloud Thickness
C             CUMHET     (LM)     Cumulus Heating  (10**14 W)
C             CUMOST     (LM)     Cumulus Moistening (10**14 W)
C             SOILMS              Soil Moisture
C             CLSAV(LM)           SCM cloud fraction SS (by volume)
C             CLDFLG(LM)          SCM Cloud flag from Radia-outcome of Rand(0,1)
C             RHC(LM)             SCM Relative Humidity saved after Cloud
C                                     routines (Qs(over water))
C             ALWP                ARM MWR cloud liquid water path (cm)
C             ADWDT               ARM d(Column_H2O)/dt (mm/hr)
C             ADWADV              ARM Column_H2O_Advection)_(mm/hr)
C             ATLWUP              ARM TOA LW Up (W/m**2)
C             ATSWDN              ARM TOA SW Dn (W/m**2)
C             ATSWIN              ARM TOA SW IN (W/m**2)
C--- Added by J.W. starting ---C
C             ENTSCM(LM,2)        SCM Entrainment rate for 2 two plume types
C             ENTDEEP(LM,2)       SCM Entrainment rate for deep convection - 2 plumes
C             MPLUMESCM(LM,2)        SCM Mass flux for 2 two plume types
C             MPLUMEDEEP(LM,2)       SCM Mass flux for deep convection - 2 plumes
C             DETRAINDEEP(LM,2,LM)  SCM Detrained convective condensate for Deep conv
C--- Added by J.W. ending ---C
C             WCUSCM(LM,2)        SCM Cumulus updraft speed for 2 two plume types
C             WCUDEEP(LM,2)       SCM  Cumulus updraft speed for deep convection - 2 plumes
C             PRCCDEEP(LM,2,LM)   SCM Precipitating convective condensate for Deep conv
C             NPRCCDEEP(LM,2,LM)  SCM Non-Precipitating conv condensate for Deep conv
C             TPALL(LM,2,LM)      SCM Plume Temperature for Deep Conv
C             MCCOND(LM,2,LM)     SCM convective condensate for deep and shallow
C             PRCCGRP(LM,2,LM)    SCM deep convective condensate graupel
C             PRCCICE(LM,2,LM)    SCM deep convective condensate ICE
C             SCM_LWP_MC          SCM MC liquid water path (kg/m2)
C             SCM_IWP_MC          SCM MC ice water path (kg/m2)
C             SCM_LWP_SS          SCM SS liquid water path (kg/m2)
C             SCM_IWP_SS          SCM SS ice water path (kg/m2)
C             SCM_WM_MC(LM)       SCM Cloud water for moist convective clouds  kg/kg
C------------------------------------------
C
C
      real*8 GZPRT(LM)
      real*8 TPRT(LM), QPRT(LM) ,TSURF, TSKIN, WMCOL(LM),
     &       QSURF,ZSURF,USURF,VSURF
      real*8 TDIFF,QDIFF
      real*8 PCOL, SVLHXCOL(LM),SVLATCOL(LM)    
      real*8 CUMFLXCOL(LM),DWNFLXCOL(LM)
      real*8 daysec,pk1000
      real*8 tt,tf,tr,tmc,tss,tbl,ZE,DZ(LM),HL(LM)

c************************************************************************
c     need to calculate CAPE/CIN since they are not currently being
c     calculated in the model
c     need to call cape_buoy_calc() subroutine
c     need to read data set psadilookup.dat
c
      INTEGER*4 mkzh
      REAL*4  cpp(LM),cpt(LM),cpq(LM),cph(LM),elev
      real*4  ccape,ccin,ghtpari,kmax
      real*4  buoy(150),buoy_dt(150),zrel(150)
c
c************************************************************************
      INTEGER L,LMIN,IC,IU
      INTEGER IPLUM,IPL,IPLUMSV      
      INTEGER IDEBUG
 
      DATA  daysec/86400./,ZSURF/10.0/

      real*8, external :: qsat

      if (NSTEPSCM.eq.0) then
          call openunit('scm.save.sige',iu,.true.,.false.)
          WRITE(iu) SIGE
          call closeunit(iu)
          call openunit('scm.output',iu_scm_diag,.true.,.false.)
      endif

      pk1000 = 1000.**kapa
c     terrain height at central facility in m
      elev = ARM_ELEV

      PCOL = P(I_TARG,J_TARG)
      TSURF = TSAVG(I_TARG,J_TARG)
      QSURF = QSAVG(I_TARG,J_TARG)
      USURF = USAVG(I_TARG,J_TARG)
      VSURF = VSAVG(I_TARG,J_TARG)
      TSKIN = GTEMP(1,4,I_TARG,J_TARG)    !!!  if ocean point not (1,4,I,J)
      scm_tmf_u = uflux1(I_TARG,J_TARG)
      scm_tmf_v = vflux1(I_TARG,J_TARG)
      if (SRDFLBTOP.gt.0.0) then
          scm_ECAS = 1.0 - (SRDFLBBOT/CSSRDBOT)
      else
          scm_ECAS = 0.0
      endif
      scm_Ssrf = (SHA*TSURF+GRAV*ZSURF)/1000.0
      scm_Hsrf = (SHA*TSURF+GRAV*ZSURF+LHE*QSURF)/1000.0
      
      do L = 1,LM
         TPRT(L) = T(I_TARG,J_TARG,L)*PK(L,I_TARG,J_TARG) 
         QPRT(L) = Q(I_TARG,J_TARG,L)
         WMCOL(L) = WM(I_TARG,J_TARG,L)
         SVLHXCOL(L) = SVLHX(L,I_TARG,J_TARG)
         SVLATCOL(L) = SVLAT(L,I_TARG,J_TARG)
         CLCVSS(L) = CLDSS(L,I_TARG,J_TARG)
         CLCVMC(L) = CLDMC(L,I_TARG,J_TARG)
         CLSAV(L) = CLDSAV(L,I_TARG,J_TARG)   
         TAUSSC(L) = TAUSS(L,I_TARG,J_TARG)
         TAUMCC(L) = TAUMC(L,I_TARG,J_TARG)
         GZPRT(L) = GZ(I_TARG,J_TARG,L)
ccc Now use temp  in K/day and q still in kg/kg (not potential temp)
         dTtot(L) = PK(L,I_TARG,J_TARG)*dTtot(L)*daysec/dtsrc
         dqtot(L) = dqtot(L)*daysec/dtsrc
         dTfrc(L) = PK(L,I_TARG,J_TARG)*dTfrc(L)*daysec/dtsrc
         dqfrc(L) = dqfrc(L)*daysec/dtsrc
         dTrad(L) = PK(L,I_TARG,J_TARG)*dTrad(L)*daysec/dtsrc
         dTHmc(L) = PK(L,I_TARG,J_TARG)*dTHmc(L)*daysec/dtsrc
         dqmc(L) = dqmc(L)*daysec/dtsrc
         dTHbl(L) = PK(L,I_TARG,J_TARG)*dTHbl(L)*daysec/dtsrc
         dqbl(L) = dqbl(L)*daysec/dtsrc
         dTHss(L) = PK(L,I_TARG,J_TARG)*dTHss(L)*daysec/dtsrc
         dqss(L) = dqss(L)*daysec/dtsrc
         dTradlw(L) = PK(L,I_TARG,J_TARG)*dTradlw(L)*daysec/dtsrc
         dTradsw(L) = PK(L,I_TARG,J_TARG)*dTradsw(L)*daysec/dtsrc
c        write(iu_scm_prt,18) NSTEPSCM,L,
c    *      dTtot(L),dTfrc(L),dTrad(L),
c    *      dTHmc(L),dTHbl(L),dTHss(L) 
c  18    format(1x,'N L dT tot frc rad mc bl ss ',i5,i5,6(f12.3))
      enddo      
c
c calculate SG_HGT and SG_RH
      ZE = 0.d0
      do L=1,LM
         DZ(L) = ((SGE_P(L)-SGE_P(L+1))/SG_P(L))*((RGAS/GRAV)*TPRT(L))
         SG_HGT(L) = ZE + DZ(L)/2.0
         SG_RH(L) = QPRT(L)/QSAT(TPRT(L),LHE,SG_P(L))
c        write(iu_scm_prt,22) L,SG_P(L),TPRT(L),ZE,SG_HGT(L),
c    &        QPRT(L),SG_RH(L)
c22      format(1x,'l p t ze z ',i5,f10.2,f10.2,f12.2,f12.2,
c    &        '  Q RH ',f10.6,f10.5)
         ZE = ZE + DZ(L)
      enddo
      do L=1,LM
         SG_CLDICE(L)=0.0
         SG_CLDWAT(L)=0.0
         if (SVLHXCOL(L).gt.2400000.and.SVLHXCOL(L).lt.2800000.) then
             SG_CLDWAT(L)=WMCOL(L)
         elseif (SVLHXCOL(L).gt.2800000.) then
             SG_CLDICE(L)=WMCOL(L)
         endif
c        write(iu_scm_prt,23) L,WMCOL(L)*1000.0,SVLHXCOL(L),
c    &               SG_CLDWAT(L)*1000.0,SG_CLDICE(L)*1000.
c23      format(1x,'L WM SVLHX  cld wat ice ',i5,f10.5,f15.0,2(f10.5))
      enddo
      CUMFLXCOL = 0.d0
      DWNFLXCOL = 0.d0
      do LMIN = 1,LM
         do IC = 1,2
            do L=1,LM
               CUMFLXCOL(L) = CUMFLXCOL(L) + CUMFLX(L,IC,LMIN)
               DWNFLXCOL(L) = DWNFLXCOL(L) + DWNFLX(L,IC,LMIN)
c              if (CUMFLX(L,IC,LMIN).gt.0.d0.or.
c    &             DWNFLX(L,IC,LMIN).gt.0.d0)
c    &         write(iu_scm_prt,81) LMIN,IC,L,CUMFLX(L,IC,LMIN),
c    &            DWNFLX(L,IC,LMIN)
c81            format(1x,'DIAG  LMIN IC L CUMFLX DWNFLX ',i5,i5,i5,
c    &            2(f12.6))
            enddo
         enddo
      enddo
c     do L = 1,LM
c        if (CUMFLXCOL(L).gt.0.d0.or.DWNFLXCOL(L).gt.0.d0)
c    &       write(iu_scm_prt,82) L,CUMFLXCOL(L),DWNFLXCOL(L)
c82      format(1x,'DIAG  L CUMFLX DWNFLX ',i5,f12.6,f12.6)
c     enddo

C--- Added by J.W. starting ---C
c     do L=1,LM
c        ENTSCM(L,1) = 0.0
c        ENTSCM(L,2) = 0.0
c        MPLUMESCM(L,1) = 0.0
c        MPLUMESCM(L,2) = 0.0
c        WCUSCM(L,1) = 0.0
c        WCUSCM(L,2) = 0.0
c     enddo
C--- Added by J.W. ending ---C

C--- Added by J.W. starting ---C
c     do ic=1,2
c        IPLUM = 0
c        do LMIN = 1,LM
c           do L=1,LM
c              if (WCUALL(L,ic,LMIN).ne.0.0d0) then
c                  IPLUM = LMIN 
cc                 write(iu_scm_prt,24) ic,lmin,iplum
cc24               format(1x,
cc   *               'WCUPLUM found   ic lmin iplum ',i5,i5,i5)
c                  go to 25
c              endif
c           enddo
c        enddo
c25      continue
c        if (IPLUM.gt.0) then
c            do L=1,LM
c               WCUSCM(L,ic) = WCUALL(L,ic,IPLUM) 
c               MPLUMESCM(L,ic) = MPLUMEALL(L,ic,IPLUM)
c               ENTSCM(L,ic) = ENTALL(L,ic,IPLUM)
c            enddo
c        endif
c     enddo
cc    do ic=1,2
cc       do l=1,LM 
cc          write(iu_scm_prt,26) ic,l,wcudeep(l,ic)
c26         format(1x,'ic l wcudeep ',i5,i5,f10.4)
cc       enddo
cc    enddo
C--- Added by J.W. ending ---C

C     before writing out diagnostics convert cumulus diagnostics
c     do L = 1,LM
c        write(iu_scm_prt,80) L,CUMHET(L),CUMOST(L)
c 80     format(1x,'L  het mst ',i5,2(f12.3))    
c        CUMHET(L) = CUMHET(L)*10.E-13*SHA*AXYP(I_TARG,J_TARG)/(GRAV*DTSRC)
c        CUMOST(L) = CUMOST(L)*10.E-13*SHA*AXYP(I_TARG,J_TARG)/(GRAV*DTSRC)
c     enddo


c******************************************************************************
c     do calculation of CAPE and CIN
c******************************************************************************
      mkzh = LM
      do L=1,LM
         cpp(L) = SG_P(LM-L+1)
         cpt(L) = TPRT(LM-L+1)
         cpq(L) = QPRT(LM-L+1)
         cph(L) = SG_HGT(LM-L+1)
c        write(iu_scm_prt,80) L,cpp(L),cpt(L),cpq(L)*1000.,cph(L)
c80      format(1x,'for cape   L P T Q H ',i5,3(f10.4),f12.2)
      enddo
      scm_CAPE = 0.0
      scm_CIN = 0.0
cc    call cape_buoy_calc(cpp,cpt,cpq,cph,elev,mkzh,
cc   *        ccape,ccin,buoy,buoy_dt,zrel,ghtpari,kmax)
cc    scm_CAPE = ccape
cc    scm_CIN = ccin
cc    write(iu_scm_prt,90) scm_CAPE,scm_CIN
c90   format(1x,'CAPE CIN ',2(f10.5))

C
ccc   WRITE(iu_scm_diag) NSTEPSCM,PCOL,TPRT,QPRT,TSURF,TSKIN,CLCVSS,
ccc  *           CLCVMC,CLTHCK,WMCOL,SVLHXCOL,SVLATCOL,CSIZE,EFFRAD,
ccc  *           CUMFLX,CUMHET, CUMOST,SRDFLBTOP,SRNFLBTOP,TRUFLBTOP,
ccc  *           SRDFLBBOT,SRNFLBBOT,TRUFLBBOT,TRDFLBBOT,PRCSS,PRCMC,
ccc  *           EVPFLX,SHFLX,SOILMS,SRFHRLCOL,
ccc  *           TRFCRLCOL,TAUSSC,TAUMCC,SG_P,SG_T,SG_Q,
ccc  *           SG_OMEGA,APREC,ALH,ASH,AMEANPS,ATSAIR,ATSKIN,
ccc  *           ARHSAIR,SG_U,SG_V,SG_HOR_TMP_ADV,
ccc  *           SG_VER_S_ADV,SG_HOR_Q_ADV,SG_VER_Q_ADV,CLSAV,
ccc  *           CLDFLG,DWNFLX,RHC,ALWP,ADWDT,ADWADV,ATLWUP,
ccc  *           ATSWDN,ATSWIN,SG_ARSCL,PRESAV,PREMC,LHPSAV,LHPMC,
ccc  *           WCUSCM,WCUDEEP,PRCCDEEP,NPRCCDEEP,TPALL,MCCOND,
ccc  *           PRCCGRP,PRCCICE,MPLUMESCM,MPLUMEDEEP
ccc  *           ,GZPRT,ENTSCM,ENTDEEP,DETRAINDEEP
ccc  *           ,SCM_LWP_MC,SCM_IWP_MC,SCM_LWP_SS,SCM_IWP_SS
ccc  *           ,SCM_WM_MC,SRUFLBTOP,SRUFLBBOT,TRDFLBTOP 
ccc  *           ,dTtot,dqtot,dTfrc,dqfrc,dTrad


      WRITE(iu_scm_diag) NSTEPSCM,ASTIME,PCOL,SG_P,SG_HGT,TPRT,QPRT,
     *   SG_RH,WMCOL,SVLHXCOL,SVLATCOL,SG_CLDWAT,SG_CLDICE,CLCVSS,
     *   CLCVMC,SG_U,SG_V,
     *   CSIZE,EFFRAD,CUMFLXCOL,DWNFLXCOL,TAUSSC,TAUMCC,
     *   SG_T,SG_Q,SG_OMEGA,SG_HOR_TMP_ADV,SG_VER_S_ADV,SG_HOR_Q_ADV,
     *   SG_VER_Q_ADV,SG_ARSCL,dTtot,dqtot,dTfrc,dqfrc,dTrad,dTHmc,
     *   dqmc,dTHbl,dqbl,dTHss,dqss,dTradlw,dTradsw,GZPRT,TSKIN,TSURF,
     *   QSURF,USURF,VSURF,EVPFLX,SHFLX,PRCSS,PRCMC,SCM_PBL_HGT,
     *   scm_tmf_u,scm_tmf_v,scm_ECAS,scm_Ssrf,scm_Hsrf,scm_CAPE,
     *   scm_CIN,scm_FDIR_vis,scm_FDIF_vis,scm_FDIR_nir,
     *   scm_FDIF_nir,SRDFLBTOP,SRNFLBTOP,TRUFLBTOP,
     *   SRDFLBBOT,SRNFLBBOT,TRUFLBBOT,TRDFLBBOT,SRUFLBBOT,
     *   SRUFLBTOP,TRDFLBTOP,TRNFLBBOT,CSSRNTOP,CSTRUTOP,
     *   CSSRNBOT,CSTRNBOT,CSSRDBOT,SRFHRLCOL,TRFCRLCOL,APREC,
     *   ATSAIR,ATSKIN,
     *   ALH,ASH,AMEANPS,PLUME_MIN,PLUME_MAX,SCM_H,SCM_HSAT,
     *   ARM_H,ARM_HSAT,isccp_sunlit,isccp_ctp,
     *   isccp_tauopt,isccp_lowcld,isccp_midcld,isccp_highcld,
     *   isccp_fq,isccp_boxtau,isccp_boxptop

      write(iu_scm_prt,99) NSTEPSCM,ASTIME
 99   format(//1x,'END OF TIME STEP      NSTEPSCM ASTIME ',
     &      i6,f10.4)
      write(iu_scm_prt,100) NSTEPSCM,PCOL+PTOP,TSURF,TSKIN 
100   format(1x,'NSTEP P tsurfair tskin  ',i5,3(f10.4))
      write(iu_scm_prt,105) TSURF,QSURF,USURF,VSURF
105   format(1x,'near surf  t q u v ',f10.2,f10.6,2(f10.4)) 
      write(iu_scm_prt,110) PRCSS,PRCMC
110   format(1x,'PRCSS MC ',2(f10.4))
      write(iu_scm_prt,120) SRNFLBBOT,SRNFLBTOP,SRDFLBBOT,SRDFLBTOP, 
     *              TRUFLBTOP,TRUFLBBOT,TRDFLBBOT
120   format(1x,'rad  sw nbot top dbot top ',4(f10.2),
     *       '    lw utop bot dbot ',3(f10.2)) 
      write(iu_scm_prt,125) EVPFLX,SHFLX,ALH,ASH
125   format(1x,'EVPFLX SHFLX ',2(f12.6),'  ALH ASH ',2(f10.4))

c
      IDEBUG = 0

      do L=1,LM
         write(iu_scm_prt,140) L,SG_P(L),TPRT(L),SG_T(L),
     +                 QPRT(L)*1000.0,SG_Q(L)*1000.0,    
     +                 WMCOL(L)*1000.0,TAUSSC(L),TAUMCC(L),
     +                 CLCVSS(L)*100.,CLCVMC(L)*100.,SG_ARSCL(L)
 140     format(1x,i3,f8.2,' T ',2(f7.2),'  Q ',2(f7.3),'  WM ',
     +          f7.4,'  tauss mc ',2(f8.3),'  camtss mc ',
     +          2(f7.3),' arscl ',f6.2)
      enddo 

      write(iu_scm_prt,151) isccp_sunlit,isccp_ctp,isccp_tauopt,
     &       isccp_lowcld,isccp_midcld,isccp_highcld
151   format(1x,'ISCCP diags   sunlit cldptop tau  ',I4,f8.2,f8.2,
     &           '    low mid high ',3(f8.3))
c     do ic = 1,ncol
c        write(iu_scm_prt,152) ic,isccp_boxptop(ic),isccp_boxtau(ic)
c152     format(1x,' ic   ptop   tau  ',i5,f9.2,f10.4)
c     enddo
      write(iu_scm_prt,153) (isccp_tau(ic),ic=1,7)
 153  format(1x,10x,7(f10.4))
      do L = 1,npres
         write(iu_scm_prt,154) isccp_press(L),
     &                          (isccp_fq(ic,L),ic=1,ntau)
 154     format(1x,I10,7(f10.4))
      enddo

      RETURN 

      END SUBROUTINE SCM_DIAG 

      SUBROUTINE CLDRESET

      USE RESOLUTION, only : LM
      USE MODEL_COM, only : I_TARG,J_TARG,WM,NSTEPSCM
      USE DYNAMICS, only : PTOLD
      USE CLOUDS_COM, only : TTOLD,QTOLD,SVLHX,RHSAV,CLDSAV,
     *                       CLDSAV1
      USE SCMCOM, only : NRESET,NRAMP,CBPTOLD,CBTTOLD,CBQTOLD,
     *                   CBWM,CBSVLHX,CBRHSAV,CBCLDSAV,
     *                   CBCLDSAV1,iu_scm_prt
      USE SCMDIAG, only : BBNSTEPSCM

      IMPLICIT NONE

      INTEGER L,IRAMP

      IRAMP=NRESET-NRAMP


      write(iu_scm_prt,*) 'enter cldreset   nstepscm iramp BBNSTEP',
     &                     nstepscm,iramp,BBNSTEPSCM(iramp)

      PTOLD(I_TARG,J_TARG) = CBPTOLD(iramp)

      do L=1,LM
         TTOLD(L,I_TARG,J_TARG) = CBTTOLD(L,iramp)
         QTOLD(L,I_TARG,J_TARG) = CBQTOLD(L,iramp)
         WM(I_TARG,J_TARG,L) = CBWM(L,iramp)
         SVLHX(L,I_TARG,J_TARG) = CBSVLHX(L,iramp)
         RHSAV(L,I_TARG,J_TARG) = CBRHSAV(L,iramp)
         CLDSAV(L,I_TARG,J_TARG) = CBCLDSAV(L,iramp)
         CLDSAV1(L,I_TARG,J_TARG) = CBCLDSAV1(L,iramp)
         write(iu_scm_prt,200) L,TTOLD(L,I_TARG,J_TARG),
     &           QTOLD(L,I_TARG,J_TARG)*1000.,WM(I_TARG,J_TARG,L)*1000.,
     *           SVLHX(L,I_TARG,J_TARG),
     *           CLDSAV(L,I_TARG,J_TARG)*100.,CLDSAV1(L,I_TARG,J_TARG)
 200     format(1x,'overwriting  L T Q WM LHX CLDSAV ',i5,f10.2,
     &           f10.5,f10.5,f11.0,2(f10.4))
      enddo

      return

      END SUBROUTINE CLDRESET




c*******get the crrspding hgt prof to the buoy prof*******************c
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine cape_buoy_calc(prs,tmk,qvp,ght,ter,mkzh,
     &                          cape,cin,buoy,buoy_dt,zrel,ghtpari,kmax)
c
c   If i3dflag=1, this routine calculates CAPE and CIN (in m**2/s**2,
c   or J/kg) for every grid point in the entire 3D domain (treating
c   each grid point as a parcel).  If i3dflag=0, then it
c   calculates CAPE and CIN only for the parcel with max theta-e in
c   the column, (i.e. something akin to Colman's MCAPE).  By "parcel",
c   we mean a 500-m deep parcel, with actual temperature and moisture
c   averaged over that depth.
c
c   In the case of i3dflag=0,
c   CAPE and CIN are 2D fields that are placed in the k=mkzh slabs of
c   the cape and cin arrays.  Also, if i3dflag=0, LCL and LFC heights
c   are put in the k=mkzh-1 and k=mkzh-2 slabs of the cin array.
c
      USE SCMCOM, only : iu_scm_prt
      USE FILEMANAGER, only : openunit,closeunit
      real*4,dimension(mkzh):: prs,tmk,qvp,ght,prsf
      real*4 kmax,ghtpari,cape,cin,ter,lcl,lfc
      integer iucape,nthte,nprs
c
      real*4, dimension(150):: buoy,buoy_dt,zrel,benaccum
c
c     include 'comconst'
      dimension psadiprs(150),psadithte(150),psaditmk(150,150)
      common /const/ rgas,rgasmd,grav,cp,cpmd,ds,dskm,
     &   gamma,gammamd,thtecon1,thtecon2,thtecon3,eps,ezero,
     &   xlhc0,xlhctd,xlhf,celkel,abscoef,ussalr,eslcon1,eslcon2,
     &   esicon1,esicon2,flmin,frmax,fbmin,ftmax,dsc,dskmc,
     &   tlclc1,tlclc2,tlclc3,tlclc4,
     &   xlatc,xlonc,rhoice,rhowat,mdateb,mhourb,rhourb,
     &   miycors,mjxcors,refrat,yicorn,xjcorn,iice,rpd,pi,
     &   rktpmps,true1,true2,abscoefi,nproj,sclht,rmsg,
     &   psadiprs,psadithte,psaditmk,ilandset,iplevdata
c
      kmax=0.
      ghtpari=0.
      buoy=0.
      buoy_dt=0.
      zrel=0.
      benaccum=0.
      cape=0.
      cin=0.
c
c   Define constants.  Many are taken from Bolton (1980, MWR 108,1046-1053).
c
      rgas=287.04  !J/K/kg
      rgasmd=.608   ! rgas_moist=rgas*(1.+rgasmd*qvp)
      cp=1004.     ! J/K/kg  Note: not using Bolton's value of 1005.7
      cp_b=1005.7     ! J/K/kg  Note: Bolton's value of 1005.7
      cpmd=.887   ! cp_moist=cp*(1.+cpmd*qvp)
      gamma=rgas/cp
      gamma_b=rgas/cp_b  ! Note: using Bolton's cp of 1005.7
      gammamd=rgasmd-cpmd  ! gamma_moist=gamma*(1.+gammamd*qvp)
      grav=9.81           ! m/s**2
      sclht=rgas*256./grav   ! 256 K is avg. trop. temp. from USSA.
      eps=0.622
      ezero=6.112  ! hPa
      xlhc0=3.1484e6   ! J/kg
      xlhctd=2370.  !
      xlhf=3.34e5
      rktpmps=1.94
      celkel=273.15
      eslcon1=17.67
      eslcon2=29.65
      esicon1=22.514
      esicon2=6.15e3
      thtecon1=3376. ! K
      thtecon2=2.54
      thtecon3=.81
      tlclc1=2840.
      tlclc2=3.5
      tlclc3=4.805
      tlclc4=55.
      rhoice=917.
      rhowat=1000.
      pi=4.*atan(1.)
      rpd=pi/180.
      abscoef=.145      ! cloud water absorption coefficient in m^2/g
      abscoefi=.272     ! cloud ice absorption coefficient in m^2/g
      ussalr=.0065      ! deg C per m
      rmsg=9.0e+9       ! indicates missing data or specification

      do k=1,mkzh
         if (k.eq.mkzh) then
c           if (iplevdata.ge.4) then ! terrain-following data
c              prsf(k)=sfp
c           else ! pressure-level data
               prsf(k)=.5*(3.*prs(k)-prs(k-1))
c           endif
         else
            prsf(k)=.5*(prs(k+1)+prs(k))
         endif
      enddo



c
c   Set up lookup table for getting temperature on a pseudoadiabat.
c   (Borrow the unit number for the stationlist, just for the moment.)
c
ccccccc fname='psadilookup.dat'
ccc   iucape=8
c
      call openunit ('SCMcape',iucape,.false.,.true.)
c     open(unit=iucape,file='psadilookup.dat',form='formatted',
c    &     status='old')
      do i=1,14
         read(iucape,*)
      enddo
      read(iucape,172) nthte,nprs
 172  format(1x,i5,i5)
c     write(iu_scm_prt,*) 'CAPE  nthte nprs ',nthte,nprs
      if (nthte.ne.150.or.nprs.ne.150) then
         write(iu_scm_prt,*)
     &      'Number of pressure or theta_e levels in lookup table'
         write(iu_scm_prt,*) 'file not = 150.  Check lookup table file.'
         stop
      endif
      read(iucape,173) (psadithte(jt),jt=1,nthte)
      read(iucape,173) (psadiprs(ip),ip=1,nprs)
      read(iucape,173) ((psaditmk(ip,jt),ip=1,nprs),jt=1,nthte)
 173  format(5e15.7)
      call closeunit(iucape)
ccc   close(iucape)


         cape=0.
         cin=0.
c
c
c      Find parcel with max theta-e in lowest 3 km AGL.
c
         ethmax=-1.
         do k=mkzh,1,-1
            if (ght(k)-ter.lt.3000.) then
               q=max(qvp(k),1.e-15)
               t=tmk(k)
               p=prs(k)
               e=q*p/(eps+q)
               tlcl=tlclc1/(log(t**tlclc2/e)-tlclc3)+tlclc4
c              eth=t*(1000./p)**(gamma*(1.+gammamd*q))*
c    &            exp((thtecon1/tlcl-thtecon2)*q*(1.+thtecon3*q))
cccc      to be exactly consistently with bolton's formula
               eth=t*(1000./p)**(gamma_b*(1.+gammamd*q))*
     &            exp((thtecon1/tlcl-thtecon2)*q*(1.+thtecon3*q))
               if (eth.gt.ethmax) then
                  klev=k
                  ethmax=eth
               endif
            endif
         enddo
         kpar=klev
c
c      Establish average properties of that parcel
c         (over depth of approximately davg meters)
c
c         davg=.1
         eps=0.622
         davg=500.
         pavg=davg*prs(kpar)*grav/
     &      (rgas*tmk(kpar)*(eps+qvp(kpar))/(eps*(1.+qvp(kpar))))

         p2=min(prs(kpar)+.5*pavg,prsf(mkzh))
         p1=p2-pavg
         totthe=0.
         totqvp=0.
         totprs=0.
         do k=mkzh,2,-1
            if (prsf(k).le.p1) goto 35
            if (prsf(k-1).ge.p2) goto 34
            p=prs(k)
            pup=prsf(k)
            pdn=prsf(k-1)
            q=max(qvp(k),1.e-15)
            th=tmk(k)*(1000./p)**(gamma*(1.+gammamd*q))
            pp1=max(p1,pdn)
            pp2=min(p2,pup)
            if (pp2.gt.pp1) then
               deltap=pp2-pp1
               totqvp=totqvp+q*deltap
               totthe=totthe+th*deltap
               totprs=totprs+deltap
            endif
 34         continue
         enddo
 35      continue
         qvppari=totqvp/totprs
         tmkpari=(totthe/totprs)*(prs(kpar)/1000.)**
     &      (gamma*(1.+gammamd*qvp(kpar)))

c
c
c   Calculate temperature and moisture properties of parcel
c     (Note, qvppari and tmkpari already calculated above for 2D case.)
c
      prspari=prs(kpar)
      ghtpari=ght(kpar)
      gammam=gamma*(1.+gammamd*qvppari)
      cpm=cp*(1.+cpmd*qvppari)
c
      e=max(1.e-20,qvppari*prspari/(eps+qvppari))
      tlcl=tlclc1/(log(tmkpari**tlclc2/e)-tlclc3)+tlclc4
c     ethpari=tmkpari*(1000./prspari)**(gamma*(1.+gammamd*qvppari))*
c    &   exp((thtecon1/tlcl-thtecon2)*qvppari*
c    &   (1.+thtecon3*qvppari))
ccccc to be exactly consistent with bolton's formula
      ethpari=tmkpari*(1000./prspari)**(gamma_b*(1.+gammamd*qvppari))*
     &   exp((thtecon1/tlcl-thtecon2)*qvppari*
     &   (1.+thtecon3*qvppari))
      zlcl=ghtpari+(tmkpari-tlcl)/(grav/cpm)
c
c   Calculate buoyancy and relative height of lifted parcel at
c   all levels, and store in bottom up arrays.  Add a level at the LCL,
c   and at all points where buoyancy is zero.
c
      kk=0 ! for arrays that go bottom to top
      write(0,*) 'inside cape_buoy_calc   kk ',kk
      ilcl=0
      if (ghtpari.ge.zlcl) then
c
c      initial parcel already saturated or supersaturated.
c
         ilcl=2
         klcl=1
      endif
      do k=kpar,1,-1
 33      kk=kk+1                ! for arrays that go bottom to top
         write(0,*) 'inside k loop  k,kk ',k,kk
         if (ght(k).lt.zlcl) then ! model level is below LCL
            qvplift=qvppari
            tmklift=tmkpari-grav/cpm*(ght(k)-ghtpari)
            tvenv=tmk(k)*(eps+qvp(k))/(eps*(1.+qvp(k)))
            tvlift=tmklift*(eps+qvplift)/(eps*(1.+qvplift))
            ghtlift=ght(k)
         elseif (ght(k).ge.zlcl.and.ilcl.eq.0) then
c
c         This model level and previous model level straddle the LCL,
c         so first create a new level in the bottom-up array, at the LCL.
c
            tmklift=tlcl
            qvplift=qvppari
            facden=ght(k)-ght(k+1)
            fac1=(zlcl-ght(k+1))/facden
            fac2=(ght(k)-zlcl)/facden
            tmkenv=tmk(k+1)*fac2+tmk(k)*fac1
            qvpenv=qvp(k+1)*fac2+qvp(k)*fac1
            tvenv=tmkenv*(eps+qvpenv)/(eps*(1.+qvpenv))
            tvlift=tmklift*(eps+qvplift)/(eps*(1.+qvplift))
            ghtlift=zlcl
            ilcl=1
         else
            tmklift=tonpsadiabat(ethpari,prs(k))
            eslift=ezero*exp(eslcon1*(tmklift-celkel)/
     &         (tmklift-eslcon2))
            if ((tmklift-eslcon2).lt.0.0) eslift = 0.0
c           write(iu_scm_prt,*) 'tvlift else eslift exps  ',eslcon1,
c    &                         (tmklift-celkel),(tmklift-eslcon2)
c           write(iu_scm_prt,*) 'tvlift else eslift ',k,kk,ezero,
c    &         eslcon1,tmklift,celkel,eslcon2,eslift
            qvplift=eps*eslift/(prs(k)-eslift)
c           write(iu_scm_prt,*) 'tvlift else qvplift ',k,kk,eps,eslift,
c    &           prs(K),qvplift
            tvenv=tmk(k)*(eps+qvp(k))/(eps*(1.+qvp(k)))
            tvlift=tmklift*(eps+qvplift)/(eps*(1.+qvplift))
c           write(iu_scm_prt,*) 'tvlift else ',k,kk,tvlift,eps,qvplift
            ghtlift=ght(k)
         endif
c        write(iu_scm_prt,*) 'kk tvlift tvenv tvenv ',kk,tvlift,
c    &                        tvenv
         buoy(kk)=grav*(tvlift-tvenv)/tvenv  ! buoyancy
         write(0,*) 'calc buoy(kk)   ',kk,buoy(kk)
         buoy_dt(kk)=tvlift-tvenv            ! buoyancy
         zrel(kk)=ghtlift-ghtpari
c         if (buoy(kk)*buoy(kk-1).lt.0.0) then      ! changed April 26, 2006
         if (kk.gt.1) then
             if (buoy(kk)*buoy(kk-1).lt.0.0) then
c
c            Parcel ascent curve crosses sounding curve, so create a new level
c            in the bottom-up array at the crossing.
c
               kk=kk+1
               buoy(kk)=buoy(kk-1)
               write(0,*) 'calc buoynewlevel ',kk,buoy(kk)
               buoy_dt(kk)=buoy_dt(kk-1)
               zrel(kk)=zrel(kk-1)
               buoy(kk-1)=0.
               buoy_dt(kk-1)=0.
               if (kk.gt.2) then
                 zrel(kk-1)=zrel(kk-2)+
     &            buoy(kk-2)/(buoy(kk-2)-buoy(kk))*(zrel(kk)-zrel(kk-2))
               endif
             endif
         endif
         if (ilcl.eq.1) then
            klcl=kk
            ilcl=2
            goto 33
         endif
      enddo
      kmax=kk
      write(0,*) 'kmax = ',kmax
c     if (kmax.gt.150) then
c        write(iu_scm_prt,*)'in capecalc3d: kmax got too big. kmax=',kmax
c        stop
c     endif
c
c   Get the accumulated buoyant energy from the parcel's starting
c   point, at all levels up to the top level.
c
      benaccum(1)=0.0
      benamin=9e9
      ikmax = kmax
      do k=2,ikmax
         dz=zrel(k)-zrel(k-1)
         benaccum(k)=benaccum(k-1)+.5*dz*(buoy(k-1)+buoy(k))
         if (benaccum(k).lt.benamin) then
            benamin=benaccum(k)
         endif
      enddo
c
c     Determine equilibrium level (EL), which we define as the highest
c     level of non-negative buoyancy above the LCL. Note, this may be
c     the top level if the parcel is still buoyant there.
c
      ikmax = kmax
      do k=ikmax,klcl,-1
         if (buoy(k).ge.0.) then
            kel=k   ! k of equilibrium level
            goto 50
         endif
      enddo
c
c   If we got through that loop, then there is no non-negative
c   buoyancy above the LCL in the sounding.  In these situations,
c   both CAPE and CIN will be set to -0.1 J/kg.  Also, where CAPE is
c   non-zero, CAPE and CIN will be set to a minimum of +0.1 J/kg, so
c   that the zero contour in either the CIN or CAPE fields will
c   circumscribe regions of non-zero CAPE.
c
      cape=-0.1
      cin=-0.1
      klfc=kmax
c
      goto 102
c
 50   continue
c
c   If there is an equilibrium level, then CAPE is positive.  We'll
c   define the level of free convection (LFC) as the point below the
c   EL, but at or above the LCL, where accumulated buoyant energy is a
c   minimum.  The net positive area (accumulated buoyant energy) from
c   the LFC up to the EL will be defined as the CAPE, and the net
c   negative area (negative of accumulated buoyant energy) from the
c   parcel starting point to the LFC will be defined as the convective
c   inhibition (CIN).
c
c   First get the LFC according to the above definition.
c
      benamin=9e9
      klfc=kmax
c     write(iu_scm_prt,*) 'klcl kel ',klcl,kel
      do k=klcl,kel
c        write(iu_scm_prt,*) 'k = ',k
         if (benaccum(k).lt.benamin) then
            benamin=benaccum(k)
            klfc=k
         endif
      enddo
c
c   Now we can assign values to cape and cin
c
      cape=max(benaccum(kel)-benamin,0.1)
      cin=max(-benamin,0.1)
c
c   CIN is uninteresting when CAPE is small (< 100 J/kg), so set
c   CIN to -.1 in that case.
c
      if ( cape .lt. 100. ) cin = -0.1
 102  continue

      lcl=zrel(klcl)+ghtpari-ter   ! meters AGL
      lfc=zrel(klfc)+ghtpari-ter   ! meters AGL

      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      function tonpsadiabat(thte,prs)
c
c   This function gives the temperature (in K) on a moist adiabat
c   (specified by thte in K) given pressure in hPa.  It uses a
c   lookup table, with data that was generated by the Bolton (1980)
c   formula for theta_e.
c
c     include 'comconst'
      dimension psadiprs(150),psadithte(150),psaditmk(150,150)
      common /const/ rgas,rgasmd,grav,cp,cpmd,ds,dskm,
     &   gamma,gammamd,thtecon1,thtecon2,thtecon3,eps,ezero,
     &   xlhc0,xlhctd,xlhf,celkel,abscoef,ussalr,eslcon1,eslcon2,
     &   esicon1,esicon2,flmin,frmax,fbmin,ftmax,dsc,dskmc,
     &   tlclc1,tlclc2,tlclc3,tlclc4,
     &   xlatc,xlonc,rhoice,rhowat,mdateb,mhourb,rhourb,
     &   miycors,mjxcors,refrat,yicorn,xjcorn,iice,rpd,pi,
     &   rktpmps,true1,true2,abscoefi,nproj,sclht,rmsg,
     &   psadiprs,psadithte,psaditmk,ilandset,iplevdata
c
c     First check if pressure is less than min pressure in lookup table.
c     If it is, assume parcel is so dry that the given theta-e value can
c     be interpretted as theta, and get temperature from the simple dry
c     theta formula.
c
      if (prs.lt.psadiprs(150)) then
         tonpsadiabat=thte*(prs/1000.)**gamma
         return
      endif
c
c   Otherwise, look for the given thte/prs point in the lookup table.
c
      do jtch=1,150-1
         if (thte.ge.psadithte(jtch).and.thte.lt.psadithte(jtch+1)) then
            jt=jtch
            goto 213
         endif
      enddo
      jt=-1
 213  continue
      do ipch=1,150-1
         if (prs.le.psadiprs(ipch).and.prs.gt.psadiprs(ipch+1)) then
            ip=ipch
            goto 215
         endif
      enddo
      ip=-1
 215  continue
      if (jt.eq.-1.or.ip.eq.-1) then
         write(iu_scm_prt,*)
     &      'Outside of lookup table bounds. prs,thte=',prs,thte
         stop
      endif
      fracjt=(thte-psadithte(jt))/(psadithte(jt+1)-psadithte(jt))
      fracjt2=1.-fracjt
      fracip=(psadiprs(ip)-prs)/(psadiprs(ip)-psadiprs(ip+1))
      fracip2=1.-fracip
      if (psaditmk(ip,jt).gt.1e9.or.psaditmk(ip+1,jt).gt.1e9.or.
     &    psaditmk(ip,jt+1).gt.1e9.or.psaditmk(ip+1,jt+1).gt.1e9) then
         write(iu_scm_prt,*)
     &      'Tried to access missing tmperature in lookup table.'
         write(iu_scm_prt,*)
     &      'Prs and Thte probably unreasonable. prs,thte=',prs,thte
         stop
      endif
      tonpsadiabat=fracip2*fracjt2*psaditmk(ip  ,jt  )+
     &       fracip *fracjt2*psaditmk(ip+1,jt  )+
     &       fracip2*fracjt *psaditmk(ip  ,jt+1)+
     &       fracip *fracjt *psaditmk(ip+1,jt+1)
c
      return
      end
