      module soilbgc

!@sum Routines to simulate soil biogeochemistry:
!@sum microbial dynamics, C & N pools, respiration and N fluxes.

      use ent_const
      use ent_types
      use ent_pfts
      implicit none

      public soil_bgc

      contains
      
!***********************************************************************      
      subroutine soil_bgc(dtsec, pp)
      use patches, only : print_Tpool

      real*8, intent(in) :: dtsec  !main ent time step (s)
      type(patch),pointer :: pp
      !----Local----------
      real*8 :: Soilmoist(N_CASA_LAYERS) !soil volumetric moisture, depth-structured -PK
      integer :: ivt                     !ivt = pft
      real*8 :: Soiltemp(N_CASA_LAYERS)  !soil temperature (C), depth-structured
      real*8 :: clayfrac                 !fractional clay content in soil
      real*8 :: sandfrac                 !fractional clay content in soil
      real*8 :: siltfrac                 !fractional clay content in soil
      real*8 :: Tpool(PTRACE,NPOOLS,N_CASA_LAYERS)  !total plant and soil C,N pools (gC/m2), depth-structured
      real*8 :: Cflux                    !total soil C flux to atm (gC/m2/s) **C flux, NOT CO2!**

      ! do nothing if no vegetation
      if ( .not. ASSOCIATED(pp%tallest) ) return

      ivt = pp%tallest%pft     
!      Soilmoist = pp%Soilmoist    !soil moist varies by patch  !**omit for now** -PK 7/6/07
      Soilmoist(:) = pp%cellptr%Soilmoist(:) 
      Soiltemp(:) = pp%cellptr%Soiltemp(:)  !soil temp, texture vary by cell
!       print *, __FILE__,__LINE__,'soiltemp=',soiltemp !***test*** -PK 7/24/07  
!       print *, __FILE__,__LINE__,'soilmoist=',soilmoist !***test*** -PK 7/24/07  
      clayfrac = pp%cellptr%soil_texture(3) !in GHY.f, texture order is sand,loam,clay,peat(+bedrock) -PK 7/13/06
      sandfrac = pp%cellptr%soil_texture(1)
! use siltfrac = 1 - (clayfrac + sandfrac) or 0.4*"loam" for now -PK 6/14/06
      siltfrac = 0.4d0*pp%cellptr%soil_texture(2)  !hack to allow use of GCM soil textures -PK  
      Tpool(:,:,:) = pp%Tpool(:,:,:) !Added - NYK 7/27/06
!       print *, __FILE__,__LINE__,'Tpool=',Tpool(CARBON,:,:) !***test*** -PK 7/24/07  
      
      call casa_bgfluxes(dtsec, Soilmoist,  ivt    
     &                 , Soiltemp, clayfrac, sandfrac, siltfrac
     &                 , Tpool, Cflux)

      !* Convert Cflux from gC/m2/s to kgC/m2/s 
      !* and assign patch soil_resp, Tpool 
      pp%Soil_resp = Cflux*1.d-3 
      pp%Tpool(:,:,:) = Tpool(:,:,:)
!       print *, __FILE__,__LINE__,'pp%Tpool=',pp%Tpool(CARBON,:,:) !***test*** -PK 7/24/07  

      end subroutine soil_bgc

!***********************************************************************
      subroutine casa_bgfluxes(dtsec, Soilmoist, ivt
     &                 ,Soiltemp, clayfrac, sandfrac, siltfrac 
     &                 ,Tpool,Cflux)

      implicit none

! ------------------------ input/output variables -----------------
! input 
      real*8,intent(in) :: dtsec           !main ent time step (s)
      real*8,intent(in) :: Soilmoist(N_CASA_LAYERS)  !soil moisture, depth-structured
      integer,intent(in) :: ivt            !ivt = pp%tallest%pft (assigned in soilbgc)
      real*8,intent(in) :: Soiltemp(N_CASA_LAYERS)   !soil temperature, depth-structured
      real*8,intent(in) :: clayfrac        !fractional clay content in soil
      real*8,intent(in) :: sandfrac        !fractional sand content in soil
      real*8,intent(in) :: siltfrac        !fractional silt content in soil
      
! i/o
      real*8,intent(inout) :: Tpool(PTRACE,NPOOLS,N_CASA_LAYERS)  !total plant and soil C,N pools, depth-structured
! output
      real*8,intent(out) :: Cflux      !total respiration flux to atm (gC/m2/s) !main output from this routine -PK 6/15/06

! ------------------------ local variables ------------------------
      integer ::  n,m
      integer ::  ipool, i
      real*8 :: Cdead(N_CASA_LAYERS) 
!decomp coefs      
      real*8 :: fact_soilmic(N_PFT)
      real*8 :: fact_slow(N_PFT)
      real*8 :: fact_passive(N_PFT)
      real*8 :: eff(NRESP_POOLS)          
      real*8 :: frac_donor(NRESP_POOLS)   
      real*8 :: kdt(N_PFT,NPOOLS)
!met variables and functions
      real*8,dimension(N_CASA_LAYERS) :: atmp, bgmoist, bgtemp
      real*8,dimension(N_CASA_LAYERS) :: Wlim    !water limitation function for CASA (func of next 2)
      real*8 :: watopt  !"optimal" water content for ET
      real*8 :: watdry  !water content when ET stops (~wilting point)
      !watopt, watdry are funcs of next 3 (which are funcs of clayfrac,sandfrac -- see lsmtci.F) 
      real*8 :: watsat  !saturated volumetric soil water content (porosity)
      real*8 :: smpsat  !soil matric potential at saturation (mm)
      real*8 :: bch     !clapp and hornberger "b"

      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS)    !amt C lost per pool (gC/m2)
      real*8 :: Resp(PTRACE,NPOOLS,N_CASA_LAYERS)     !amt C lost to atm per pool (gC/m2), used to calculate Cflux
      real*8 :: poolsum                 !accumulates Resp

!!! HACK  set Closs to something so that 1:NLIVE are initialized
!!! please fix! ??
      Closs(:,:,:) = 0.d0
      
! define rate coefs fact_soilmic, fact_slow, fact_passive (from casatci.F) -PK 5/25/06
        if(ivt.eq.CROPS)then      
          fact_soilmic(ivt) = 1.25d0
          fact_slow(ivt)    = 1.5d0
          fact_passive(ivt) = 1.5d0
        else
          fact_soilmic(ivt) = 1.d0
          fact_slow(ivt)    = 1.d0
          fact_passive(ivt) = 1.d0
        end if
        if(ivt.eq.0)then
          fact_soilmic(ivt) = 0.d0
          fact_slow(ivt)    = 0.d0
          fact_passive(ivt) = 0.d0
        end if

* Maximum RATE CONSTANTS FOR EACH POOL SCALED TO LENGTH OF TIME STEP
* For small delta_t, kdt   is the same as annK*delta_t
*iyf:  Consider dM/dt =  - M/tau
*iyf:  Analytic solution:  M(t) = M0 exp (-t/tau)
*iyf:  Integrate Flux*dt from t=(n-1)*dt to t=n*dt:
*iyf:  integral = M(n*dt) - M[(n-1)*dt]
*iyf:           = - M[(n-1)*dt] {1 - exp [-dt/tau]}
*iyf:    approx = M[(n-1)*dt] * [dt/tau]
*iyf:  variable kdt was previously Krate in CASA

** NOTE: kdt will be used for WOOD and dead pools only
**       so no need to worry about adding stressCD to annK here
      do n = 1, NPOOLS
        do m = 1, N_PFT    
          kdt(m,n)= 1.d0 - (exp(-annK(m,n)*dtsec))
        enddo
      enddo  


*---Step 1: heterotrophic (microbial) respiration -------------------------
C.. Initialize respiration fluxes each timestep
C.. Note: these should be over dead pools only (see resp_pool_index)
*iyf:  Resp in unit of gC/m2/timestep

      Resp(:,:,:) = 0.d0

*---Step 1a: TEMPERATURE AND MOISTURE CONSTRAINTS ON DECOMP 
! TEMPERATURE DEPENDENCE
      !* Original CASA.
            bgtemp(:) = (Q10 ** ((Soiltemp(:) - 30.d0) / 10.d0))  !original CASA function -PK
      !* Exponential Arrhenius temperature response (close to Q10 but smaller at 0 Celsius).
!      do i = 1,N_CASA_LAYERS
!        if (Soiltemp(i).le.-33.d0) then
!          bgtemp(i) = 0.d0
!        else
!          bgtemp(i) = exp(308.56d0*(1/63.15d0 - 1/(Soiltemp(i)+33.15))) !exp(308.56d0*(1/(KELVIN+30.d0-240.d0) - 1/(Soiltemp(:)+KELVIN-240.d0)
!        endif
!      enddo
      !* Function f(Tsoil) from DelGrosso et al. (Biogeoch. 73, 2005)**  -PK 2/07
        !allows for variable Q10 rather than fixed 
!            bgtemp(:) = 0.56d0+
!     &      (1.46d0*atan(PI*0.0309d0*(Soiltemp(:)-15.7d0))) / PI

      !* Linear fit to Del Grosso ftemp. Min=0.125 @Soiltemp=0, Max=1.0 @Soiltemp=30 degC- NK
!      bgtemp = max(0.125d0, 
!     &     0.125d0 + (1.d0-0.125d0)/(30.d0-0.d0)*Soiltemp(:)) !linear, no upper cap -NK
!      bgtemp = max(0.125d0, min(1.d0,
!     &     0.125d0 + (1.d0-0.125d0)/(30.d0-0.d0)*Soiltemp(:))) !linear, Del Grosso intercept -NK
!      bgtemp = max(0.d0, min(1.d0,
!     &     0.d0 + (1.d0-0.d0)/(30.d0-0.d0)*Soiltemp(:))) !linear, zero intercept at freezing - NK
!      bgtemp = max(0.046d0, min(1.d0,
!     &     0.046d0 + (1.d0-0.046d0)/(30.d0-0.d0)*Soiltemp(:))) !linear, log-fitted intercept to original Del Grosso data -NK
      
      !* S-function fit to Del Grosso. - NK
!      bgtemp=1.15d0*(1.d0/(1.d0+EXP(-0.14d0*(Soiltemp(:)-17.d0))))

! MOISTURE DEPENDENCE
* mimic calculation of bevap in surphy.F to get Wlim
* but use Soilmoist,Soiltemp instead of h2osoi,tsoi 
*   watdry = water content when evapotranspiration stops = wp
*   watdry x 0.5d0 = rough estimate of hygroscopic point for microbial respiration.
! Equations take soil textures in percents rather than fractions.
            watsat = 0.489d0 - 0.00126d0*(sandfrac*100.d0)
            smpsat = -10.d0 * ( 10.d0**(1.88d0-
     &           (0.0131d0*sandfrac*100.d0)) )
            bch = 2.91d0 + 0.159d0*(clayfrac*100.d0)
            watdry = 0.5d0*watsat * (-316230.d0/smpsat) ** (-1.d0/bch)
            watopt = watsat * (-158490.d0/smpsat) ** (-1.d0/bch)
!! 03/11/21 note: there are no limits on Wlim, except if Soiltemp < 0
!!Wlim ultimately used to get total C loss per pool,  
!!to both atm and other pools (see step 1b below) -PK 6/8/06
           do n=1,N_CASA_LAYERS
            if (Soiltemp(n) .gt. 0.d0) then
!               Wlim(n) = min( max(Soilmoist(n)-watdry,0.d0) /  !original CASA function -PK
!     &                   (watopt-watdry), 1.d0)
        !**function RWC from DelGrosso et al., 2005** -PK 2/07
!               Wlim(n) = (Soilmoist(n)-watdry)/(watopt - watdry)
               Wlim(n) = max(0.0d0,
     &             (Soilmoist(n)-watdry)/(watsat - watdry)) !Made this REW instead of Wlim - NK
            else
               Wlim = 0.01d0
            end if
           end do

!           bgmoist(:) = 0.25d0 + 0.75d0*Wlim(:)  !original CASA function -PK
        !**functions f(RWC), Rh from DelGrosso et al. (Biogeoch. 73, 2005)** -PK 2/07
!           bgmoist(:) = 5.d0 *
!     &                 (0.287d0+(atan(PI*0.009d0*(Wlim(:)-17.47d0)))/PI)
           bgmoist(:) = min(1.d0,0.01d0 + 
     &                     (1.d0-0.01d0)/(0.7-0.d0)*Wlim(:)) !linear - NK
           atmp(:) = bgtemp(:) * bgmoist(:)

*---Step 1b: DETERMINE loss of C FROM EACH DEAD POOL (donor) PER TIMESTEP
*iyf:  Closs is the amount of carbon each pool loses in gC/m2/timestep.
*iyf:  A fraction of Closs is transferred to another pool (receiver pool), 
*iyf:  the remainder to the atm (Resp).
*iyf:  The distribution of Closs is done in subroutine casa_respire, Step 1c. 
          do i=1,N_CASA_LAYERS
            do n = 1, NDEAD
             ipool = NLIVE + n
             Cdead(i) = Tpool(CARBON,ipool,i) * kdt(ivt,ipool) * atmp(i)
             Closs(CARBON,ipool,i) = Cdead(i)
            end do  

** adjust pools
            Closs(CARBON,SURFSTR,i) = Closs(CARBON,SURFSTR,i)
     &                  * lignineffect(ivt)
            Closs(CARBON,SOILSTR,i) = Closs(CARBON,SOILSTR,i)
     &                  * lignineffect(ivt)
            Closs(CARBON,SOILMIC,i) = Closs(CARBON,SOILMIC,i)
     &                  * (1.d0-(0.75d0*(siltfrac+clayfrac)))  !see Potter et al 1993 for ref -PK 6/8/06
     &                  * fact_soilmic(ivt)
            Closs(CARBON,SLOW,i) = Closs(CARBON,SLOW,i)
     &                  * fact_slow(ivt)
            Closs(CARBON,PASSIVE,i) = Closs(CARBON,PASSIVE,i)
     &                  * fact_passive(ivt)
!       print *, __FILE__,__LINE__,'closs=',closs(CARBON,:,:) !***test*** -PK 7/24/07  

*iyf:  limits on loss from dead pools.
ciyf - No need to track limits on loss rate
ciyf (only need to track limits on inventories)
          !do i=1,N_CASA_LAYERS
           do n = NLIVE+1,NPOOLS
            Closs(CARBON,n,i)=MIN(Closs(CARBON,n,i),Tpool(CARBON,n,i))
           end do
           
          end do  !N_CASA_LAYERS
 
*---Step 1c:  SOM C DECOMPOSITION
*iyf:  update Tpool's:  C is transferred between donor and receiver pools 
*iyf:  Resp is the amount of C to atm, in gC/m2/timestep

* MICROBIAL EFFICIENCIES FOR PARTICULAR FLOWS
        eff( 1) =  0.45    ! SLOW,PASSIVE
        eff( 2) =  0.45    ! SLOW,SOILMIC
        eff( 3) =  0.40    ! SURFMET,SURFMIC
        eff( 4) =  0.40    ! SURFSTR,SURFMIC
        eff( 5) =  0.70    ! SURFSTR,SLOW
        eff( 6) =  0.45    ! SOILMET,SOILMIC
        eff( 7) =  0.45    ! SOILSTR,SOILMIC
        eff( 8) =  0.70    ! SOILSTR,SLOW
        eff( 9) =  0.40    ! CWD,SURFMIC
        eff(10) =  0.70    ! CWD,SLOW
        eff(11) =  0.40    ! SURFMIC,SLOW
        eff(12) =  0.85 - (0.68 * (siltfrac+clayfrac))  ! SOILMIC,PASSIVE 
        eff(13) =  0.85 - (0.68 * (siltfrac+clayfrac))  ! SOILMIC,SLOW 
        eff(14) =  0.45    ! PASSIVE,SOILMIC
* EXTRA RESPIRATION TRANSFER EFFICIENCIES
        frac_donor( 1) =  0.003 + (0.009*clayfrac)
        frac_donor( 2) =  1.0 - frac_donor(1)
        frac_donor( 3) =  1.0
        frac_donor( 4) =  1.0 - structurallignin(ivt)
        frac_donor( 5) =  structurallignin(ivt)
        frac_donor( 6) =  1.0
        frac_donor( 7) =  1.0 - structurallignin(ivt)
        frac_donor( 8) =  structurallignin(ivt)
        frac_donor( 9) =  1.0 - woodligninfract
        frac_donor(10) =  woodligninfract
        frac_donor(11) =  1.0
        frac_donor(12) =  0.003 + (0.032*clayfrac)
        frac_donor(13) =  1.0 - frac_donor(12)
        frac_donor(14) =  1.0

! casa_respire determines amount of C transferred from each "donor" pool and 
! where it goes ("receiver" pool), in addition to total C respired (Resp-->Cflux) -PK 6/15/06
      call casa_respire(ivt ,eff ,frac_donor
     &                 ,Closs, Resp, Tpool)
     
*---Step 2-----------------------------------------------------------
* CALCULATE NITROGEN POOLS !***keep this for future use*** -PK 8/23/06                                 
!            do n = 1,NPOOLS
!              Tpool(Nitrogen,n)=Tpool(Carbon,n)/CNratio(n)
!            enddo   

!--- Step 2.5 ------------------------------------------------------
! **calculate vertical C transport from each pool (assume closed at top&bottom)
! and adjust each pool accordingly** -PK 5/07
      if (N_CASA_LAYERS == 2)  call vertCtransport(dtsec, Tpool)  !only applies when multiple layers present

*--- Step 3 --------------------------------------------------------
*  Get Total C Fluxes to atm = Sum over respiring (dead) pools/dtsec
*  Note: Resp for live pools should be zero
*iyf:  Cflux in gC/m2/s
         poolsum = 0.0
         do i=1,N_CASA_LAYERS
          do n = NLIVE+1,NPOOLS
            poolsum=poolsum+Resp(CARBON,n,i)
          end do  
         end do  !over dead pools -PK
         Cflux=poolsum/dtsec           !total C flux to atm (gC/m2/s)

C.. mod 03/03/11 change closs/nloss to be g/m2/sec instead of g/m2/timestep
         Closs(CARBON,:,:)=Closs(CARBON,:,:)/dtsec
!         Closs(NITROGEN,:,:)=Closs(NITROGEN,:,:)/dtsec  !not used (yet) -PK  

      return
      end subroutine casa_bgfluxes
      
!***********************************************************************
      subroutine casa_respire(ivt ,eff ,frac_donor
     &                       ,Closs, Resp, Tpool)

      implicit none

! ------------------------ input/output variables -----------------
! input
      integer,intent(in) :: ivt            !ivt = pp%tallest%pft (assigned in soilbgc)
      real*8,intent(in) :: eff(NRESP_POOLS)        !decomp coefs
      real*8,intent(in) :: frac_donor(NRESP_POOLS) !decomp coefs
      real*8,intent(in) :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS)    !amt C lost per pool (gC/m2)
      
! i/o
      real*8,intent(inout) :: Tpool(PTRACE,NPOOLS,N_CASA_LAYERS)  !total plant and soil C,N pools
      
! output
      real*8,intent(out) :: Resp(PTRACE,NPOOLS,N_CASA_LAYERS)  !C lost to atm per pool (gC/m2)-->Cflux in casa_bgfluxes (step 3)
      
! ------------------------ local variables ------------------------
      integer ::  resp_pool_index(2,NRESP_POOLS)   
      integer ::  irtype,n
      integer ::  donor_pool
      integer ::  recvr_pool
      real*8 :: Out(N_CASA_LAYERS)

! ----------------------------------------------------------
! set values of parameters/constants used in CASA Respiration
! these are in order in which respiration is called
      resp_pool_index = reshape ( 
     1              (/SLOW     ,PASSIVE,
     2              SLOW     ,SOILMIC,
     3              SURFMET  ,SURFMIC,
     4              SURFSTR  ,SURFMIC,
     5              SURFSTR  ,SLOW   ,
     6              SOILMET  ,SOILMIC,
     7              SOILSTR  ,SOILMIC,
     8              SOILSTR  ,SLOW   ,
     9              CWD      ,SURFMIC,
     a              CWD      ,SLOW   ,
     b              SURFMIC  ,SLOW   ,
     c              SOILMIC  ,PASSIVE,
     d              SOILMIC  ,SLOW   ,
     e              PASSIVE  ,SOILMIC/),
     &              (/2,NRESP_POOLS/)   )

! Loop over all respiring pools
      do n=1,N_CASA_LAYERS
       do irtype = 1, NRESP_POOLS
         donor_pool = resp_pool_index(1,irtype)
         recvr_pool = resp_pool_index(2,irtype)
         Out(n)  = Closs(CARBON,donor_pool,n) * frac_donor(irtype)
         Tpool(CARBON,donor_pool,n) = Tpool(CARBON,donor_pool,n)
     &                                  - Out(n)
         Tpool(CARBON,recvr_pool,n) = Tpool(CARBON,recvr_pool,n)  
     &                                  + (Out(n) * eff(irtype))
         Resp(CARBON,donor_pool,n) =  Resp(CARBON,donor_pool,n) 
     &                                  + Out(n) * (1.d0 - eff(irtype))

!** 01/10/01 make sure donor pool does not fall below zero
        if(Tpool(CARBON,donor_pool,n).le.0.0d0)then
            Tpool(CARBON,donor_pool,n)=0.d0
        end if
        
       end do  !irtype - over all 14 respiring pools
       
      end do  !N_CASA_LAYERS

      return
      end subroutine casa_respire
      
!***********************************************************************
      subroutine vertCtransport(dtsec, Tpool)
      !calculates inter-layer C transport (based on Baisden et al., GBC 16, 2002) -PK 5/07
      implicit none
      
! input 
      real*8,intent(in) :: dtsec           !main ent time step (s)
! i/o
      real*8,intent(inout) :: Tpool(PTRACE,NPOOLS,N_CASA_LAYERS)  !total plant and soil C,N pools
! ------------------------ local variables ------------------------
      integer ::  i, n, numlayers
      !coefs for vertical C transport (based on Baisden et al., GBC 16, 2002) -PK 5/07        
      !layer depths(mm) for use w/vertical transport coefs. **need to know these a priori**  -PK 5/07
      real*8 :: dz(N_CASA_LAYERS) !move to ent_const at some point? -PK 7/07   
      !rate coefficients 
      real*8, dimension(N_CASA_LAYERS) :: vfast, vslow, vpass
      
      !**temporary workaround to avoid compile trap and allow for N_CASA_LAYERS=2** -PK 7/9/07
        numlayers = 2
      !might increase dz(2) to 1700 mm, i.e. put bottom of column at 2 m -PK
        dz(1) = 300.d0  !layer 1 = 300 mm
#ifdef NCASA2
        dz(2) = 700.d0  !layer 2 = 700 mm   
#endif
        
      !calculate 1st-order transport coefs in appropriate units (unitless) 
      do i=1,N_CASA_LAYERS
        vfast(i) = 4.d0/dz(i)/SECPY*dtsec   !v1=4.0 mm/yr for fast pool of annual grassland from Baisden et al., 2002 
        vslow(i) = 0.5d0/dz(i)/SECPY*dtsec  !v2=0.5 mm/yr for slow pool ...
        vpass(i) = 0.4d0/dz(i)/SECPY*dtsec  !v3=0.4 mm/yr for passive pool ...
      end do
      
      !adjust TPOOL to reflect vertical (inter-layer) C transport
      do n = NLIVE+1,NPOOLS
        !top(1st) layer -- loss term only
        if (n<=7 .OR. n==9 .OR. n==10) then  !for SURFMET,SURFSTR,SOILMET,SOILSTR,SURFMIC,SOILMIC  
          TPOOL(CARBON,n,1) =  (1.d0-vfast(1)) * TPOOL(CARBON,n,1) 
        else if (n==8 .OR. n==11) then  !for CWD,SLOW
          TPOOL(CARBON,n,1) =  (1.d0-vslow(1)) * TPOOL(CARBON,n,1)
        else if (n==12) then  !for PASSIVE
          TPOOL(CARBON,n,1) =  (1.d0-vpass(1)) * TPOOL(CARBON,n,1)
        end if
        !second layer -- gain from layer 1 only
        if (n==6 .OR. n==7 .OR. n==10) then  !for SOILMET,SOILSTR,SOILMIC  
          TPOOL(CARBON,n,numlayers) = TPOOL(CARBON,n,numlayers) 
     &                       + vfast(1) * TPOOL(CARBON,n,1)
        else if (n==8 .OR. n==11) then  !for CWD,SLOW
          TPOOL(CARBON,n,numlayers) = TPOOL(CARBON,n,numlayers)
     &                       + vslow(1) * TPOOL(CARBON,n,1)
        else if (n==12) then  !for PASSIVE
          TPOOL(CARBON,n,numlayers) = TPOOL(CARBON,n,numlayers)
     &                       + vpass(1) * TPOOL(CARBON,n,1)
        end if
      end do  !loop through all 9 soil C ('dead') pools
       
      end subroutine vertCtransport
      
!***********************************************************************

      end module soilbgc
