#include "rundeck_opts.h"

!#define ROUGHL_HACK

      module PBL_DRV
      use SOCPBL, only : t_pbl_args, xdelt
      implicit none

      private

      public t_pbl_args, pbl, xdelt

      contains

      SUBROUTINE PBL(I,J,ITYPE,PTYPE,pbl_args)
!@sum  PBL calculate pbl profiles for each surface type
!@+        Contains code common for all surfaces
!@auth Greg. Hartke/Ye Cheng
!@ver  1.0
!@var DDMS downdraft mass flux in kg/(m^2 s), (i,j)
!@var TDN1 downdraft temperature in K, (i,j)
!@var QDN1 downdraft humidity in kg/kg, (i,j)

      USE CONSTANT, only :  rgas,grav,omega2,deltx,teeny
      USE MODEL_COM, only : t,q,u,v,ls1
#ifdef SCM
      USE MODEL_COM, only : I_TARG,J_TARG
      USE SCMCOM, only : iu_scm_prt,SCM_SURFACE_FLAG
      USE SCMDIAG, only : SCM_PBL_HGT
#endif
      USE GEOM, only : sinlat2d
      USE DYNAMICS, only : pmid,pk,pedn,pek
     &    ,DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
     &    ,ua=>ualij,va=>valij
      USE CLOUDS_COM, only : ddm1
      USE CLOUDS_COM, only : DDMS,TDN1,QDN1,DDML
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,trdn1
#ifdef TRACERS_DRYDEP
     &    ,trradius,trpdens,tr_mm
#endif
#endif
#ifdef TRACERS_AMP
     & ,AMP_MODES_MAP,AMP_NUMB_MAP,ntmAMP
      USE AMP_AEROSOL, only : DIAM, AMP_dens,AMP_TR_MM
      USE AERO_SETUP,  only : CONV_DPAM_TO_DGN
#endif

      use SOCPBL, only : npbl=>n, zgs, advanc
      USE PBLCOM
      use QUSDEF, only : mz
      use SOMTQ_COM, only : tmom


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: I,J  !@var I,J grid point
      INTEGER, INTENT(IN) :: ITYPE  !@var ITYPE surface type
      REAL*8, INTENT(IN) :: PTYPE  !@var PTYPE percent surface type
      type (t_pbl_args) :: pbl_args

      REAL*8, parameter :: dbl_max=3000., dbl_max_stable=500. ! meters
      real*8, parameter :: S1byG1=.57735d0
      REAL*8 Ts

#ifdef TRACERS_ON
      integer nx,n
#endif
c
      REAL*8 ztop,zpbl,pl1,tl1,pl,tl,tbar,thbar,zpbl1,coriol
      REAL*8 qtop,utop,vtop,ufluxs,vfluxs,tfluxs,qfluxs,psitop,psisrf
      INTEGER LDC,L,k
!@var uocean,vocean ocean/ice velocities for use in drag calulation
!@var evap_max maximal evaporation from unsaturated soil
!@var  fr_sat fraction of saturated soil
!@var ZS1    = height of the first model layer (m)
!@var TGV    = virtual potential temperature of the ground (K)
!@+            (if xdelt=0, TGV is the actual temperature)
!@var TKV    = virtual potential temperature of first model layer (K)
!@+            (if xdelt=0, TKV is the actual temperature)
!@var WS     = magnitude of the surface wind (m/s)
!@var PSI    = angular diff. btw geostrophic and surface winds (rads)
!@var WG     = magnitude of the geostrophic wind (m/s)
!@var HEMI   = 1 for northern hemisphere, -1 for southern hemisphere
!@var TG     = bulk ground temperature (K)
!@var ELHX   = latent heat for saturation humidity (J/kg)
!@var dskin  = skin-bulk SST difference (C)
!@VAR QSOL   = solar heating (W/m2)
      real*8 zs1,psi,hemi
!@var POLE   = .TRUE. if at the north or south pole, .FALSE. otherwise
c      logical pole

!**** the following is output from advance (mostly passed through pbl_args)
!@var US     = x component of surface wind, positive eastward (m/s)
!@var VS     = y component of surface wind, positive northward (m/s)
!@var WSGCM  = magnitude of the GCM surface wind - ocean currents (m/s)
!@var WSPDF  = mean surface wind calculated from PDF of wind speed (m/s)
!@var WS     = magn. of GCM surf wind - ocean curr + buoyancy + gust (m/s)
!@var TSV    = virtual potential temperature of the surface (K)
!@+            (if xdelt=0, TSV is the actual temperature)
!@var QS     = surface value of the specific moisture
!@var DBL    = boundary layer height (m)
!@var KMS    = momentum transport coefficient at ZGS (m**2/s)
!@var KHS    = heat transport coefficient at ZGS (m**2/s)
!@var KHQ    = moist transport coefficient at ZGS (m**2/s)
!@var USTAR  = friction speed (square root of momentum flux) (m/s)
!@var CM     = drag coefficient (dimensionless surface momentum flux)
!@var CH     = Stanton number   (dimensionless surface heat flux)
!@var CQ     = Dalton number    (dimensionless surface moisture flux)
!@var z0m   = roughness length for momentum,
!@+           prescribed for itype=3,4 but computed for itype=1,2 (m)
!@var z0h   = roughness length for temperature (m)
!@var z0q   = roughness length for water vapor (m)
!@var UG     = eastward component of the geostrophic wind (m/s)
!@var VG     = northward component of the geostrophic wind (m/s)
!@var MDF    = downdraft mass flux (m/s)
!@var WINT   = integrated surface wind speed over sgs wind distribution
      real*8 :: dbl,kms,kqs,cm,ch,cq,z0m,z0h,z0q,ug,vg,w2_1,mdf
!@var dtdt_gcm temp. tendency from processes other than turbulence (K/s)
      real*8 ::  dpdxr,dpdyr,dpdxr0,dpdyr0,dtdt_gcm
      real*8 ::  mdn  ! ,mup
      real*8, dimension(npbl) :: upbl,vpbl,tpbl,qpbl
      real*8, dimension(npbl-1) :: epbl
#if defined(TRACERS_ON)
!@var  tr local tracer profile (passive scalars)
      real*8, dimension(npbl,pbl_args%ntx) :: tr
      real*8, dimension(ntm) :: trnradius,trndens,trnmm
#endif

ccc extract data needed in driver from the pbl_args structure
      zs1 = pbl_args%zs1
      hemi = pbl_args%hemi
c      pole = pbl_args%pole

      ! Redelsperger et al. 2000, eqn(13), J. Climate, 13, 402-421
      ! tprime,qprime are the pertubation of t and q due to gustiness

      ! pick up one of the following two expressions for gusti

      ! for down draft:
      mdn=max(DDMS(i,j), -0.07d0)
      pbl_args%gusti=log(1.-600.4d0*mdn-4375.*mdn*mdn)

      ! for up draft:
      ! mup=min(DDMS(i,j), 0.1d0)
      ! pbl_args%gusti=log(1.+386.6d0*mup-1850.*mup*mup)

C        ocean and ocean ice are treated as rough surfaces
C        roughness lengths from Brutsaert for rough surfaces

      IF (ITYPE.GT.2) THEN
        Z0M=ROUGHL(I,J)           ! 30./(10.**ROUGHL(I,J))
      ENDIF
      ztop=zgs+zs1  ! zs1 is calculated before pbl is called
      IF (pbl_args%TKV.EQ.pbl_args%TGV)
     &     pbl_args%TGV = 1.0001d0*pbl_args%TGV

      ! FIND THE PBL HEIGHT IN METERS (DBL) AND THE CORRESPONDING
      ! GCM LAYER (L) AT WHICH TO COMPUTE UG AND VG.
      ! LDC IS THE LAYER TO WHICH DRY CONVECTION/TURBULENCE MIXES

c       IF (TKV.GE.TGV) THEN
c         ! ATMOSPHERE IS STABLE WITH RESPECT TO THE GROUND
c         ! DETERMINE VERTICAL LEVEL CORRESPONDING TO HEIGHT OF PBL:
c         ! WHEN ATMOSPHERE IS STABLE, CAN COMPUTE DBL BUT DO NOT
c         ! KNOW THE INDEX OF THE LAYER.
c         ustar=ustar_pbl(itype,i,j)
c         DBL=min(0.3d0*USTAR/OMEGA2,dbl_max_stable)
c         if (dbl.le.ztop) then
c           dbl=ztop
c           L=1
c         else
c           ! FIND THE VERTICAL LEVEL NEXT HIGHER THAN DBL AND
c           ! COMPUTE Ug and Vg THERE:
c           zpbl=ztop
c           pl1=pmid(1,i,j)         ! pij*sig(1)+ptop
c           tl1=t(i,j,1)*(1.+xdelt*q(i,j,1))*pk(1,i,j)
c           do l=2,ls1
c             pl=pmid(l,i,j)        !pij*sig(l)+ptop
c             tl=t(i,j,l)*(1.+xdelt*q(i,j,l))*pk(l,i,j) !virtual,absolute
c             tbar=thbar(tl1,tl)
c             zpbl=zpbl-(rgas/grav)*tbar*(pl-pl1)/(pl1+pl)*2.
c             if (zpbl.ge.dbl) exit
c             pl1=pl
c             tl1=tl
c           end do
c         endif

c     ELSE
        ! ATMOSPHERE IS UNSTABLE WITH RESPECT TO THE GROUND
        ! LDC IS THE LEVEL TO WHICH DRYCNV/ATURB MIXES.
        ! FIND DBL FROM LDC.  IF BOUNDARY
        ! LAYER HEIGHT IS LESS THAN DBL_MAX, ASSIGN LDC TO L, OTHERWISE
        ! MUST FIND INDEX FOR NEXT MODEL LAYER ABOVE 3 KM:

        LDC=max(int(DCLEV(I,J)+.5d0),1)
        IF (LDC.EQ.0) LDC=1
        if (ldc.eq.1) then
          dbl=ztop
          l=1
        else
          zpbl=ztop
          pl1=pmid(1,i,j)                             ! pij*sig(1)+ptop
          tl1=t(i,j,1)*(1.+xdelt*q(i,j,1))*pk(1,i,j)  ! expbyk(pl1)
          zpbl1=ztop
          do l=2,ldc
            pl=pmid(l,i,j)                            ! pij*sig(l)+ptop
            tl=t(i,j,l)*(1.+xdelt*q(i,j,l))*pk(l,i,j) ! expbyk(pl)
            tbar=thbar(tl1,tl)
            zpbl=zpbl-(rgas/grav)*tbar*(pl-pl1)/(pl1+pl)*2.
            if (zpbl.ge.dbl_max) then
              zpbl=zpbl1
              exit
            endif
            pl1=pl
            tl1=tl
            zpbl1=zpbl
          end do
          l=min(l,ldc)
          dbl=zpbl
        endif

c     ENDIF

      coriol=sinlat2d(i,j)*omega2
      qtop=q(i,j,1)

      utop = ua(1,i,j)
      vtop = va(1,i,j)
      ug   = ua(L,i,j)
      vg   = va(L,i,j)

#ifdef SCM
      if (I.eq.I_TARG.and.J.eq.J_TARG) then
c         write(iu_scm_prt,'(a33,i4,i4,2(f8.3),i5,2(f8.3),f10.2)')
c    &            'in PBL j i utop vtop L ug vg dbl ',
c    &             j,i,utop,vtop,L,ug,vg,dbl
          SCM_PBL_HGT = dbl
      endif
#endif
      upbl(:)=uabl(:,itype,i,j)
      vpbl(:)=vabl(:,itype,i,j)
      tpbl(:)=tabl(:,itype,i,j)
      qpbl(:)=qabl(:,itype,i,j)
      epbl(1:npbl-1)=eabl(1:npbl-1,itype,i,j)

#ifdef TRACERS_ON
      do nx=1,pbl_args%ntx
        tr(:,nx)=trabl(:,pbl_args%ntix(nx),itype,i,j)
      end do

      do n = 1,ntm
#ifdef TRACERS_DRYDEP
           trnradius(n) = trradius(n)
           trndens(n)   = trpdens(n)
           trnmm(n)     = tr_mm(n)
#endif
#ifdef TRACERS_AMP
       if (n.le.ntmAMP) then
        if(AMP_MODES_MAP(n).gt.0) then
         if(DIAM(i,j,1,AMP_MODES_MAP(n)).gt.0.) then
          if(AMP_NUMB_MAP(n).eq. 0) then    ! Mass
        trnradius(n)=0.5*DIAM(i,j,1,AMP_MODES_MAP(n))
          else                              ! Number
        trnradius(n)=0.5*DIAM(i,j,1,AMP_MODES_MAP(n))
     +               *CONV_DPAM_TO_DGN(AMP_MODES_MAP(n))
          endif

           call AMPtrdens(i,j,1,n)
           call AMPtrmass(i,j,1,n)

          trndens(n) =AMP_dens(i,j,1,AMP_MODES_MAP(n))
          trnmm(n)   =AMP_TR_MM(i,j,1,AMP_MODES_MAP(n))
        endif   
        endif   
       endif 
#endif  
      enddo
#endif

      cm=cmgs(itype,i,j)
      ch=chgs(itype,i,j)
      cq=cqgs(itype,i,j)
      dpdxr  = DPDX_BY_RHO(i,j)
      dpdyr  = DPDY_BY_RHO(i,j)
      dpdxr0 = DPDX_BY_RHO_0(i,j)
      dpdyr0 = DPDY_BY_RHO_0(i,j)

      mdf = ddm1(i,j)

!!! put some results from above to pbl_args
      pbl_args%dbl = dbl
      pbl_args%ug = ug
      pbl_args%vg = vg
      pbl_args%wg = sqrt(ug*ug+vg*vg)
      pbl_args%cm = cm
      pbl_args%ch = ch
      pbl_args%cq = cq

#ifdef USE_PBL_E1
      pbl_args%ddml_eq_1=.false.
#else
      pbl_args%ddml_eq_1=DDML(i,j).eq.1
#endif

      ! if ddml_eq_1=.false.,
      ! i.e., either USE_PBL_E1 or DDML(i,j) is not 1,
      ! then tdns,qdns,tprime,qprime are not in use

      if (pbl_args%ddml_eq_1) then
        pbl_args%tdns=TDN1(i,j)*pek(1,i,j)/pk(1,i,j)
        pbl_args%qdns=QDN1(i,j)
#ifdef TRACERS_ON
        do nx=1,pbl_args%ntx
          pbl_args%trdn1(nx)=TRDN1(pbl_args%ntix(nx),i,j)
        end do
#endif
      else
        pbl_args%tdns=0.d0
        pbl_args%qdns=0.d0
#ifdef TRACERS_ON
        pbl_args%trdn1(:)=0.
#endif
      endif

      dtdt_gcm = (pbl_args%tkv - t1_after_aturb(i,j)*pek(1,i,j))/
     &     pbl_args%dtsurf
      call advanc( pbl_args,coriol,utop,vtop,qtop,ztop,mdf
     &     ,dpdxr,dpdyr,dpdxr0,dpdyr0
     &     ,dtdt_gcm,u1_after_aturb(i,j),v1_after_aturb(i,j)
     &     ,i,j,itype
     &     ,kms,kqs,z0m,z0h,z0q,w2_1,ufluxs,vfluxs,tfluxs,qfluxs
     &     ,upbl,vpbl,tpbl,qpbl,epbl
#if defined(TRACERS_ON)
     &     ,tr,ptype,trnradius,trndens,trnmm
#endif
     &     )

      uabl(:,itype,i,j)=upbl(:)
      vabl(:,itype,i,j)=vpbl(:)
      tabl(:,itype,i,j)=tpbl(:)
      qabl(:,itype,i,j)=qpbl(:)
      eabl(1:npbl-1,itype,i,j)=epbl(1:npbl-1)
#ifdef TRACERS_ON
      do nx=1,pbl_args%ntx
        trabl(:,pbl_args%ntix(nx),itype,i,j)=tr(:,nx)
      end do
#endif

      cmgs(itype,i,j)=pbl_args%cm
      chgs(itype,i,j)=pbl_args%ch
      cqgs(itype,i,j)=pbl_args%cq
      ipbl(itype,i,j)=1  ! ipbl is used in subroutine init_pbl

      psitop=atan2(vg,ug+teeny)
      psisrf=atan2(pbl_args%vs,pbl_args%us+teeny)
      psi   =psisrf-psitop
      ustar_pbl(itype,i,j)=pbl_args%ustar
C ******************************************************************
      TS=pbl_args%TSV/(1.+pbl_args%QSRF*xdelt)
      if ( ts.lt.152d0 .or. ts.gt.423d0 ) then
        write(6,*) 'PBL: Ts bad at',i,j,' itype',itype,ts
        if (ts.gt.1d3) call stop_model("PBL: Ts out of range",255)
        if (ts.lt.50d0) call stop_model("PBL: Ts out of range",255)
      end if
#ifdef SCM
      if (SCM_SURFACE_FLAG.eq.1.and.
     &                   I.eq.I_TARG.and.J.eq.J_TARG) then
c         write(iu_scm_prt,'(a26,2(i5),5(f10.4))')
c    &       'SCM I J WSAVG TS QS US VS ',
c    &       I,J,WSAVG(I,J),TSAVG(I,J),QSAVG(I,J)*1000.,
c    &       USAVG(I,J),VSAVG(I,J)
      else
         WSAVG(I,J)=WSAVG(I,J)+pbl_args%WS*PTYPE
         TSAVG(I,J)=TSAVG(I,J)+TS*PTYPE
  !       if(itype.ne.4) QSAVG(I,J)=QSAVG(I,J)+QSRF*PTYPE
         QSAVG(I,J)=QSAVG(I,J)+pbl_args%QSRF*PTYPE
         USAVG(I,J)=USAVG(I,J)+pbl_args%US*PTYPE
         VSAVG(I,J)=VSAVG(I,J)+pbl_args%VS*PTYPE
      endif
#else
      WSAVG(I,J)=WSAVG(I,J)+pbl_args%WS*PTYPE
      TSAVG(I,J)=TSAVG(I,J)+TS*PTYPE
  !    if(itype.ne.4) QSAVG(I,J)=QSAVG(I,J)+QSRF*PTYPE
      QSAVG(I,J)=QSAVG(I,J)+pbl_args%QSRF*PTYPE
      USAVG(I,J)=USAVG(I,J)+pbl_args%US*PTYPE
      VSAVG(I,J)=VSAVG(I,J)+pbl_args%VS*PTYPE
#endif
      TAUAVG(I,J)=TAUAVG(I,J)+pbl_args%CM*pbl_args%WS*pbl_args%WS*PTYPE
      uflux(I,J)=uflux(I,J)+ufluxs*PTYPE
      vflux(I,J)=vflux(I,J)+vfluxs*PTYPE
      tflux(I,J)=tflux(I,J)+tfluxs*PTYPE
      qflux(I,J)=qflux(I,J)+qfluxs*PTYPE

      tgvAVG(I,J)=tgvAVG(I,J)+pbl_args%tgv*PTYPE
      qgAVG(I,J)=qgAVG(I,J)+pbl_args%qg_aver*PTYPE
      w2_l1(I,J)=w2_l1(I,J)+w2_1*PTYPE

ccc put drive output data to pbl_args structure
      pbl_args%psi = psi ! maybe also should be moved to ADVANC
                         ! or completely otside of PBL* ?

      RETURN
      END SUBROUTINE PBL

      end module PBL_DRV

      subroutine init_pbl(inipbl)
c -------------------------------------------------------------
c These routines include the array ipbl which indicates if the
c  computation for a particular ITYPE was done last time step.
c Sets up the initialization of wind, temperature, and moisture
c  fields in the boundary layer. The initial values of these
c  fields are obtained by solving the static equations of the
c  Level 2 model. This is used when starting from a restart
c  file that does not have this data stored.
c -------------------------------------------------------------
      USE FILEMANAGER
      USE PARAM
      USE CONSTANT, only : lhe,lhs,tf,omega2,deltx
      USE MODEL_COM
      USE GEOM, only : imaxj,sinlat2d
!      USE SOCPBL, only : dpdxr,dpdyr,dpdxr0,dpdyr0

      USE SOCPBL, only : npbl=>n,zgs,inits,XCDpbl,ccoeff0,skin_effect
     &     ,xdelt
      USE GHY_COM, only : fearth
      USE PBLCOM
      USE DOMAIN_DECOMP_ATM, only : GRID, GET, READT_PARALLEL
      USE DOMAIN_DECOMP_1D, only : WRITET_PARALLEL
      USE DYNAMICS, only : pmid,pk,pedn,pek
     &    ,DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
     &    ,ua=>ualij,va=>valij
      USE SEAICE_COM, only : rsi,snowi
      USE FLUXES, only : gtemp
#ifdef USE_ENT
      use ent_mod, only: ent_get_exports
      use ent_com, only : entcells
#endif


      IMPLICIT NONE

C**** ignore ocean currents for initialisation.
      real*8, parameter :: uocean=0.,vocean=0.
!@var inipbl whether to init prog vars
      logical, intent(in) :: inipbl
!@var iu_CDN unit number for roughness length input file
      integer :: iu_CDN
      integer :: ilong  !@var ilong  longitude identifier
      integer :: jlat   !@var jlat  latitude identifier
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,4) ::
     *                                                      tgvdat

      integer :: itype  !@var itype surface type
      integer i,j,k,lpbl !@var i,j,k loop variable
      real*8 pland,pwater,plice,psoil,poice,pocean,
     *     ztop,elhx,coriol,tgrndv,pij,ps,psk,qgrnd
     *     ,utop,vtop,qtop,ttop,zgrnd,cm,ch,cq,ustar
      real*8 qsat
      real*8 ::  dpdxr,dpdyr,dpdxr0,dpdyr0
      real*8, dimension(npbl) :: upbl,vpbl,tpbl,qpbl
      real*8, dimension(npbl-1) :: epbl
      real*8 ug,vg
      real*8, allocatable :: buf(:,:)
      real*8 :: canopy_height, fv
      integer, save :: roughl_from_file = 0

      integer :: I_1, I_0, J_1, J_0
      integer :: I_1H, I_0H, J_1H, J_0H

       character*80 :: titrrr
       real*8 rrr(im,grid%J_STRT_HALO:grid%J_STOP_HALO)

        titrrr = "roughness length over land"
        rrr = 0.

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               J_STRT=J_0,       J_STOP=J_1)

      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

C things to be done regardless of inipbl

      call sync_param( 'roughl_from_file', roughl_from_file )
!!#if ( ! defined ROUGHL_HACK ) || ( ! defined USE_ENT )
      if ( roughl_from_file .ne. 0 ) then
        allocate ( buf(I_0H:I_1H, J_0H:J_1H) )
        call openunit("CDN",iu_CDN,.TRUE.,.true.)
        CALL READT_PARALLEL(grid,iu_CDN,NAMEUNIT(iu_CDN),buf,1)
        call closeunit(iu_CDN)
        roughl(:,:)=30./(10.**buf(:,:))
        deallocate ( buf )
      endif
!!#endif

      call sync_param( 'XCDpbl', XCDpbl )
      call sync_param( 'skin_effect', skin_effect )

      do j=J_0,J_1
        do i=I_0,I_1
C**** fix roughness length for ocean ice that turned to land ice
          if (snowi(i,j).lt.-1.and.flice(i,j).gt.0)
     &         roughl(i,j)=30./(10.**1.84d0)
          if (fland(i,j).gt.0.and.roughl(i,j) .gt. 29.d0) then
            print*,"Roughness length not defined for i,j",i,j
     *           ,roughl(i,j),fland(i,j),flice(i,j)
            print*,"Setting to .01"
            roughl(i,j)=1d-2
          end if
        end do
      end do

      call ccoeff0
      call getztop(zgs,ztop)

      if(.not.inipbl) return

      do j=J_0,J_1
      do i=I_0,I_1
        pland=fland(i,j)
        pwater=1.-pland
        plice=flice(i,j)
        psoil=fearth(i,j)
        poice=rsi(i,j)*pwater
        pocean=pwater-poice
        if (pocean.le.0.) then
          tgvdat(i,j,1)=0.
        else
          tgvdat(i,j,1)=gtemp(1,1,i,j)+TF
        end if
        if (poice.le.0.) then
          tgvdat(i,j,2)=0.
        else
          tgvdat(i,j,2)=gtemp(1,2,i,j)+TF
        end if
        if (plice.le.0.) then
          tgvdat(i,j,3)=0.
        else
          tgvdat(i,j,3)=gtemp(1,3,i,j)+TF
        end if
        if (psoil.le.0.) then
          tgvdat(i,j,4)=0.
        else
          tgvdat(i,j,4)=gtemp(1,4,i,j)+TF
        end if
      end do
      end do

      do itype=1,4
        if ((itype.eq.1).or.(itype.eq.4)) then
          elhx=lhe
        else
          elhx=lhs
        endif

        do j=J_0,J_1
          jlat=j
          do i=I_0,imaxj(j)
            coriol=sinlat2d(i,j)*omega2
            tgrndv=tgvdat(i,j,itype)
            if (tgrndv.eq.0.) then
              ipbl(itype,i,j)=0
              go to 200
            endif
            ilong=i
            pij=p(i,j)
            ps=pedn(1,i,j)    !pij+ptop
            psk=pek(1,i,j)    !expbyk(ps)
            qgrnd=qsat(tgrndv,elhx,ps)

            utop = ua(1,i,j)
            vtop = va(1,i,j)
            qtop=q(i,j,1)
            ttop=t(i,j,1)*(1.+qtop*xdelt)*psk
            t1_after_aturb(i,j) = ttop/psk
            u1_after_aturb(i,j) = utop
            v1_after_aturb(i,j) = vtop

            zgrnd=.1d0 ! formal initialization
            if (itype.gt.2) zgrnd=roughl(i,j) !         30./(10.**roughl(i,j))

            if (itype.gt.2) rrr(i,j) = zgrnd

            dpdxr  = DPDX_BY_RHO(i,j)
            dpdyr  = DPDY_BY_RHO(i,j)
            dpdxr0 = DPDX_BY_RHO_0(i,j)
            dpdyr0 = DPDY_BY_RHO_0(i,j)
#ifdef SCM
            utop = u(i,j,1)
            vtop = v(i,j,1)
            ug = utop
            vg = vtop
#endif
            call inits(tgrndv,qgrnd,zgrnd,zgs,ztop,utop,vtop,
     2                 ttop,qtop,coriol,cm,ch,cq,ustar,
     3                 uocean,vocean,ilong,jlat,itype
     &                 ,dpdxr,dpdyr,dpdxr0,dpdyr0
     &                 ,upbl,vpbl,tpbl,qpbl,epbl,ug,vg)
            cmgs(itype,i,j)=cm
            chgs(itype,i,j)=ch
            cqgs(itype,i,j)=cq

            do lpbl=1,npbl
              uabl(lpbl,itype,i,j)=upbl(lpbl)
              vabl(lpbl,itype,i,j)=vpbl(lpbl)
              tabl(lpbl,itype,i,j)=tpbl(lpbl)
              qabl(lpbl,itype,i,j)=qpbl(lpbl)
            end do

            do lpbl=1,npbl-1
              eabl(lpbl,itype,i,j)=epbl(lpbl)
            end do

            ipbl(itype,i,j)=1
            ustar_pbl(itype,i,j)=ustar

 200      end do
        end do
      end do

      !write(981) titrrr,rrr
#ifndef SCM
#ifndef CUBED_SPHERE
      call WRITET_PARALLEL(grid,981,"fort.981",rrr,titrrr)
#endif
#endif

      return
 1000 format (1x,//,1x,'completed initialization, itype = ',i2,//)
      end subroutine init_pbl

      subroutine loadbl
!@sum loadbl initiallise boundary layer calc each surface time step
!@auth Ye Cheng
c ----------------------------------------------------------------------
c             This routine checks to see if ice has
c              melted or frozen out of a grid box.
c
c For ITYPE=1 (ocean; melted ocean ice since last time step):
c  If there was no computation made for ocean at the last time step,
c  this time step may start from ocean ice result. If there was no
c  ocean nor ocean ice computation at the last time step, nothing
c  need be done. Also deals with newly created lake (from land)
c
c For ITYPE=2 (ocean ice; frozen from ocean since last time step):
c  If there was no computation made for ocean ice at the last time step,
c  this time step may start from ocean result. If there was no
c  ocean nor ocean ice computation at the last time step, nothing
c  need be done.
c
c For ITYPE=3 (land ice; frozen on land since last time step):
c  If there was no computation made for land ice at the last time step,
c  this time step may start from land result. If there was no
c  land ice nor land computation at the last time step, nothing
c  need be done.
c
c For ITYPE=4 (land; melted land ice since last time step):
c  If there was no computation made for land at the last time step,
c  this time step may start from land ice result. If there was no
c  land nor land ice computation at the last time step, nothing
c  need be done. Also deal with newly created earth (from lake)
c
c In the current version of the GCM, there is no need to check the
c  land or land ice components of the grid box for ice formation and
c  melting because pland and plice are fixed. The source code to do
c  this is retained and deleted in the update deck in the event this
c  capability is added in future versions of the model.
c ----------------------------------------------------------------------
      USE MODEL_COM
      USE GEOM, only : imaxj
      USE DOMAIN_DECOMP_ATM, only : GRID, GET
      USE PBLCOM, only : ipbl,wsavg,tsavg,qsavg,usavg,vsavg,tauavg
     &     ,uflux,vflux,tflux,qflux,tgvavg,qgavg,w2_l1
#ifdef SCM
      USE SCMCOM, only : iu_scm_prt,SCM_SURFACE_FLAG
#endif
      IMPLICIT NONE
      integer i,j  !@var i,j loop variable

      integer :: J_1, J_0, I_1, I_0
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
        do i=I_0,imaxj(j)

c ******* itype=1: Ocean

          if (ipbl(1,i,j).eq.0) then
            if (ipbl(2,i,j).eq.1) then
              call setbl(2,1,i,j)
            elseif (ipbl(4,i,j).eq.1) then ! initialise from land
              call setbl(4,1,i,j)
            endif
          endif

c ******* itype=2: Ocean ice

          if (ipbl(2,i,j).eq.0) then
            if (ipbl(1,i,j).eq.1) call setbl(1,2,i,j)
          endif

c ******* itype=3: Land ice

          if (ipbl(3,i,j).eq.0) then
            if (ipbl(4,i,j).eq.1) call setbl(4,3,i,j)
          endif

c ******* itype=4: Land

          if (ipbl(4,i,j).eq.0) then
            if (ipbl(3,i,j).eq.1) then
              call setbl(3,4,i,j)
            elseif (ipbl(1,i,j).eq.1) then
              call setbl(1,4,i,j)
            endif
          endif

C**** initialise some pbl common variables
#ifdef SCM
         IF (SCM_SURFACE_FLAG.eq.1.and.
     &                    I.eq.I_TARG.and.J.eq.J_TARG) then
c            write(iu_scm_prt,8888) I_TARG,J_TARG,
c    &           TSAVG(I_TARG,J_TARG),QSAVG(I_TARG,J_TARG),
c    &           USAVG(I_TARG,J_TARG),VSAVG(I_TARG,J_TARG),
c    &           WSAVG(I_TARG,J_TARG)
c8888        format(1x,'SCM PBL init  TS QS US VS WS ',i5,i5,5(f10.6))
         else
             WSAVG(I,J)=0.
             TSAVG(I,J)=0.
             QSAVG(I,J)=0.
             USAVG(I,J)=0.
             VSAVG(I,J)=0.
         endif
#else
          WSAVG(I,J)=0.
          TSAVG(I,J)=0.
          QSAVG(I,J)=0.
          USAVG(I,J)=0.
          VSAVG(I,J)=0.
#endif
          TAUAVG(I,J)=0.
          TGVAVG(I,J)=0.
          QGAVG(I,J)=0.
          w2_l1(I,J)=0.

          uflux(I,J)=0.
          vflux(I,J)=0.
          tflux(I,J)=0.
          qflux(I,J)=0.

          ipbl(:,i,j) = 0       ! - will be set to 1s when pbl is called

        end do
      end do

      return
      end subroutine loadbl

      subroutine setbl(itype_in,itype_out,i,j)
!@sum setbl initiallise bl from another surface type for one grid box
!@auth Ye Cheng
      USE PBLCOM, only : npbl,uabl,vabl,tabl,qabl,eabl,cmgs,chgs,cqgs
     *     ,ipbl,ustar_pbl
#ifdef TRACERS_ON
     *     ,trabl
#endif
      IMPLICIT NONE
      integer, INTENT(IN) :: itype_in,itype_out,i,j
      integer lpbl  !@var lpbl loop variable

      do lpbl=1,npbl-1
        uabl(lpbl,itype_out,i,j)=uabl(lpbl,itype_in,i,j)
        vabl(lpbl,itype_out,i,j)=vabl(lpbl,itype_in,i,j)
        tabl(lpbl,itype_out,i,j)=tabl(lpbl,itype_in,i,j)
        qabl(lpbl,itype_out,i,j)=qabl(lpbl,itype_in,i,j)
        eabl(lpbl,itype_out,i,j)=eabl(lpbl,itype_in,i,j)
      end do
      uabl(npbl,itype_out,i,j)=uabl(npbl,itype_in,i,j)
      vabl(npbl,itype_out,i,j)=vabl(npbl,itype_in,i,j)
      tabl(npbl,itype_out,i,j)=tabl(npbl,itype_in,i,j)
      qabl(npbl,itype_out,i,j)=qabl(npbl,itype_in,i,j)
#ifdef TRACERS_ON
      trabl(:,:,itype_out,i,j)=trabl(:,:,itype_in,i,j)
#endif
      cmgs(itype_out,i,j)=cmgs(itype_in,i,j)
      chgs(itype_out,i,j)=chgs(itype_in,i,j)
      cqgs(itype_out,i,j)=cqgs(itype_in,i,j)
      ustar_pbl(itype_out,i,j)=ustar_pbl(itype_in,i,j)

      return
      end subroutine setbl

      subroutine getztop(zgs,ztop)
!@sum  getztop computes the value of ztop which is the height in meters
!@+  of the first GCM layer from the surface.
!@+  This subroutine only needs to be called when the BL fields require
!@+  initialization.
!@+  This form for z1 = zgs + zs1 (in terms of GCM parameters) yields an
!@+  average value for zs1. The quantity theta was computed on the
!@+  assumption of zs1=200 m from the original 9-layer model (actually
!@+  was misconstrued as z1 = 200m when it should have been zs1 = 200m)
!@+  and is then applied to all vertical resolutions.
!@auth Greg. Hartke/Ye Cheng
!@var zgs The height of the surface layer.
!@var ztop The height of the top of the BL simulation domain.
!@+   Corresponds to averaged height of the middle of first model layer.

      USE CONSTANT, only : rgas,grav
      USE MODEL_COM, only : pednl00,psf
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: ZGS
      REAL*8, INTENT(OUT) :: ZTOP
      real*8, parameter :: theta=269.0727251d0

      ztop=zgs+0.5d0*(pednl00(1)-pednl00(2))*rgas*theta/(grav*psf)

      return
      end subroutine getztop

      SUBROUTINE CHECKPBL(SUBR)
!@sum  CHECKPBL Checks whether PBL data are reasonable
!@auth Original Development Team
!@ver  1.0
      USE DOMAIN_DECOMP_ATM, only : GRID, GET
      USE PBLCOM, only : wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg
     *     ,ustar_pbl,uflux,vflux,tflux,qflux,tgvavg,qgavg,w2_l1
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

      integer :: I_1, I_0, J_1, J_0, njpol
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, I_STRT=I_0, I_STOP=I_1,
     *               J_STRT=J_0, J_STOP=J_1)
      njpol = grid%J_STRT_SKP-grid%J_STRT

C**** Check for NaN/INF in boundary layer data
      CALL CHECK3B(wsavg(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'wsavg')
      CALL CHECK3B(tsavg(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'tsavg')
      CALL CHECK3B(qsavg(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'qsavg')
      CALL CHECK3B(dclev(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'dclev')
      CALL CHECK3B(usavg(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'usavg')
      CALL CHECK3B(vsavg(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'vsavg')
      CALL CHECK3B(tauavg(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'tauavg')
      CALL CHECK3C(ustar_pbl(:,I_0:I_1,J_0:J_1),4,I_0,I_1,J_0,J_1,NJPOL,
     &     SUBR,'ustar')

      CALL CHECK3B(uflux(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'uflux')
      CALL CHECK3B(vflux(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'vflux')
      CALL CHECK3B(tflux(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'tflux')
      CALL CHECK3B(qflux(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'qflux')

      CALL CHECK3B(tgvavg(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'tgvavg')
      CALL CHECK3B(qgavg(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'qgavg')
      CALL CHECK3B(w2_l1(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'w2_l1')

      END SUBROUTINE CHECKPBL

