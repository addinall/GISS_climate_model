#include "rundeck_opts.h"
      PROGRAM GISS_modelE
!@sum  MAIN GISS modelE main time-stepping routine
!@auth Original Development Team
!@ver  2009/05/11 (Based originally on B399)
      USE FILEMANAGER, only : openunit,closeunit,nameunit
      USE TIMINGS, only : ntimemax,ntimeacc,timing,timestr
      USE PARAM
      USE PARSER
      USE MODEL_COM
      USE DOMAIN_DECOMP_1D, ONLY : init_app,AM_I_ROOT,ESMF_BCAST
      USE DOMAIN_DECOMP_ATM, ONLY : grid,init_grid,sumxpe
      use domain_decomp_atm, only : writei8_parallel
      USE DYNAMICS
      USE RANDOM
      USE GETTIME_MOD
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM, only: mtrace
#ifdef TRAC_ADV_CPU
      USE TRACER_COM, only: mtradv
#endif
#endif
      USE DIAG_COM, only : ia_src,ia_d5s,ia_d5d,ia_filt
     &     ,oa,monacc,koa,acc_period
      USE SOIL_DRV, only: daily_earth, ground_e
#ifdef CACHED_SUBDD
      USE SUBDD_MOD, only : write_daily_files
#else
      USE SUBDAILY, only : nsubdd,init_subdd,get_subdd,reset_subdd
     &     ,accSubdd
#endif
      USE DIAG_SERIAL, only : print_diags
#ifdef BLK_2MOM
      USE mo_bulk2m_driver_gcm, ONLY: init_bulk2m_driver
      USE mo_bulk2m_driver_gcm, ONLY: cleanup_bulk2m_driver
#endif
#ifdef USE_MPP
      USE fms_mod,         only : fms_init, fms_end
#endif
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: fv_core_wrapper
      USE FV_INTERFACE_MOD, only: Initialize
      USE FV_INTERFACE_MOD, only: Run
      USE FV_INTERFACE_MOD, only: Checkpoint
      USE FV_INTERFACE_MOD, only: Finalize
      USE FV_INTERFACE_MOD, only: Compute_Tendencies
      USE FV_INTERFACE_MOD, only: init_app_clock
c$$$      USE MODEL_COM, only: clock
      USE ESMF_MOD, only: ESMF_Clock
      USE ESMF_CUSTOM_MOD, Only: vm => modelE_vm
#endif
#ifndef CUBED_SPHERE
      USE ATMDYN, only : DYNAM,SDRAG
     &     ,FILTER, COMPUTE_DYNAM_AIJ_DIAGNOSTICS
#endif
#ifdef BLK_2MOM
#ifdef TRACERS_AMP
      USE AERO_CONFIG, ONLY: NMODES
#endif
#endif
#ifdef SCM
      USE SCMCOM , only : SG_CONV,SCM_SAVE_T,SCM_SAVE_Q,SCM_DEL_T,
     &    SCM_DEL_Q,iu_scm_prt,iu_scm_diag,iu_scm_seed,jrandscm,
     &    IFLRESET,IRESET,NRESET,NRAMP
      USE SCMDIAG, only : dTtot,dqtot
#endif
      use TimerPackage_mod, only: startTimer => start
      use TimerPackage_mod, only: stopTimer => stop
      !use soil_drv, only : conserv_wtg, conserv_htg

      IMPLICIT NONE

      INTEGER K,M,MSTART,MNOW,MODD5D,months,ioerr,Ldate,istart
      INTEGER iu_VFLXO,iu_ODA
      INTEGER :: MDUM = 0

      character(len=80) :: filenm

      REAL*8, DIMENSION(NTIMEMAX) :: PERCENT
      REAL*8, DIMENSION(0:NTIMEMAX) ::TIMING_glob
      REAL*8 start,now, DTIME,TOTALT

      CHARACTER aDATE*14
      CHARACTER*8 :: flg_go='___GO___'      ! green light
      CHARACTER*8 YYYYMMDD
      integer :: iflag=1
      external sig_stop_model
      logical :: start9

C**** Command line options
      CHARACTER*32 :: ifile
      logical :: coldRestart=.false.
      integer :: iu_IFILE
      real :: lat_min=-90.,lat_max=90.,longt_min=0.,longt_max=360.
      real*8 :: tloopbegin, tloopend
#ifdef USE_FVCORE
      Character(Len=*), Parameter :: fv_config = 'fv_config.rc'
      Type (FV_CORE_WRAPPER) :: fv
      Type (ESMF_CLOCK) :: clock
      character(len=80) :: fv_fname, fv_dfname
      character(len=1)  :: suffix
#endif
#ifdef BLK_2MOM
      LOGICAL      :: ldummy=.false.
      INTEGER      :: il0,jl0,kl0,nm0
      CHARACTER*20 :: bname='kuku.txt'
      CHARACTER*15 :: sname='MODELE_mainV3: '
#endif
      integer :: I,J,L,I_0,I_1,J_0,J_1
      real*8 :: initialTotalEnergy, finalTotalEnergy
      real*8 :: gettotalenergy ! external for now
#ifdef USE_SYSUSAGE
      integer :: i_su, max_su=3
#endif

#ifdef USE_SYSUSAGE
      do i_su=0,max_su
        call sysusage(i_su,0)
      enddo
#endif

C****
C**** Processing command line options
C****
      call read_options( coldRestart, ifile )
C****
C**** Reading rundeck (I-file) options
C****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      call parse_params(iu_IFILE)
      call closeunit(iu_IFILE)

#ifdef USE_MPP
      call fms_init( )
#endif

      call init_app()
      call initializeDefaultTimers()

#ifdef SCM
      call sync_param( "J_TARG", J_TARG )
      call init_grid(grid, im, jm, lm, j_scm=j_targ)
      call openunit("scm.prt",iu_scm_prt,.false.,.false.)
#else
c initialize the atmospheric domain decomposition
c for now, CREATE_CAP is only relevant to the cubed sphere grid
      call init_grid(grid, im, jm, lm, CREATE_CAP=.true.)

#endif

      I_0 = GRID%I_STRT; I_1 = GRID%I_STOP
      J_0 = GRID%J_STRT; J_1 = GRID%J_STOP

#ifndef ADIABATIC


#ifdef TRACERS_ON
#ifdef RUNTIME_NTM
! allocation of tracer arrays in physics modules needs to know NTM
      call read_tracer_config
#endif
#endif

#endif /* ADIABATIC */
      call alloc_drv()

#if !defined(ADIABATIC) || defined( CUBED_SPHERE)
C****
C**** INITIALIZATIONS
C****
         CALL TIMER (NOW,MDUM)

#ifdef BLK_2MOM
C initialize microphysics
c        print *,sname,'im        = ',im
c        print *,sname,'jm        = ',jm
c        print *,sname,'lm        = ',lm
c        print *,sname,'dtsrc        = ',dtsrc
         kl0=12;il0=im;jl0=jm ;il0=1;jl0=1;kl0=1
         nm0 = 1 ! or whatever, put a correct value here
#ifdef TRACERS_AMP
         nm0=NMODES
#endif
c        print *,sname,'il0       = ',il0
c        print *,sname,'jl0       = ',jl0
c        print *,sname,'kl0       = ',kl0
c        print *,sname,'nm0       = ',nm0
        ldummy = init_bulk2m_driver(dtsrc,il0,jl0,kl0,nm0,bname)
        if(ldummy) then
          print *,sname,'BLK Initialization is completed...'
        else
          call stop_model("BLK Initialization is not completed: ",255)
        endif
c        print *,sname,'Before:istart,ifile = ',istart,ifile
c        print *,sname,'Before:im,jm        = ',im,jm
#endif

#endif /* ADIABATIC */


C**** Read input/ic files
#ifdef USE_FVCORE
      CALL INPUT (istart,ifile,coldRestart,clock)
#else
      CALL INPUT (istart,ifile,coldRestart)
#endif

#ifdef SCM
c***************************************************************************
c   for multiple runs of the SCM - changing seed to random number generator
c   read in seed from file
c     call openunit("scmseed.data",iu_scm_seed,.false.,.true.)
c     rewind(iu_scm_seed)
c     read(iu_scm_seed,'(1x,i15)') jrandscm
c     CALL RINIT(jrandscm)
c     write(iu_scm_prt,*) 'random number generator initialized with ',
c    &     jrandscm
c***************************************************************************
#endif

#if !defined(ADIABATIC) || defined( CUBED_SPHERE)

C**** Set run_status to "run in progress"
      if(istart > 0) call write_run_status("Run in progress...",1)

C****
C**** If run is already done, just produce diagnostic printout
C****
      IF (Itime.GE.ItimeE.and.Kradia.le.0) then ! includes ISTART<1 case
        if (ItimeE.gt.0) then
          months=(Jyear-Jyear0)*JMperY + JMON-JMON0
          call aPERIOD (JMON0,JYEAR0,months,1,0, acc_period,Ldate)
        end if
        call print_diags(0)
        if(istart < 1) then
          call stop_model ('Finished post-processing',-1)
        else
          call stop_model ('The run has already completed',13)
        end if
        ! no output files are affected
      END IF

      IF (AM_I_ROOT()) Then
         open(3,file='flagGoStop',form='FORMATTED',status='REPLACE')
         write (3,'(A8)') flg_go
         close (3)
      END IF
      call sys_signal( 15, sig_stop_model )  ! works only on single CPU
         START=NOW
         DO M=1,NTIMEACC
           START= START-TIMING(M)
         END DO
C**** INITIALIZE TIME PARAMETERS
      NSTEP=(Itime-ItimeI)*NIdyn
#ifdef SCM
      NSTEPSCM=ITIME-ITIMEI
#endif
         MODD5K=1000
      CALL DAILY(.false.)                  ! not end_of_day
      CALL daily_RAD(.false.)
      if (istart.le.9) call reset_diag(0)
      if(Kradia==10) call daily_OCEAN(.false.) ! to test OCLIM
      if (Kradia.le.0) then
        CALL daily_EARTH(.false.)          ! not end_of_day
        CALL daily_OCEAN(.false.)          ! not end_of_day
        CALL CALC_AMPK(LS1-1)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
        CALL daily_tracer(.false.)
#endif
           if (kradia.le.0) CALL CHECKT ('INPUT ')
      end if
      CALL UPDTYPE

#endif

#ifdef USE_FVCORE
C****
C**** Initialize FV dynamical core (ESMF component) if requested
C**** For restarts/continuations, FV state and import files are
C**** assumed to have been copied to the appropriate names by INPUT().
C**** For cold starts the FV interface code generates the required files.
C****
      Call Initialize(fv, istart, vm, grid%esmf_grid, clock,fv_config)
#endif

      if (AM_I_ROOT())
     *   WRITE (6,'(A,11X,A4,I5,A5,I3,A4,I3,6X,A,I4,I10)')
     *   '0NASA/GISS Climate Model (re)started',
     *   'Year',JYEAR,aMON,JDATE,', Hr',JHOUR,
     *   'Internal clock: DTsrc-steps since 1/1/',Iyear1,ITIME

         CALL TIMER (NOW,MELSE)
C****
C**** Open and position output history files if needed
C****
C**** Files for an accumulation period (1-12 months)
      write(aDATE(1:7),'(a3,I4.4)') aMON0(1:3),Jyear0
      if (Kvflxo.ne.0) then
        call openunit('VFLXO'//aDATE(1:7),iu_VFLXO,.true.,.false.)
        call io_POS(iu_VFLXO,Itime,2*im*jm*koa,Nday) ! real*8-dim -> 2*
      end if
#ifndef CACHED_SUBDD
C**** Initiallise file for sub-daily diagnostics, controlled by space-
C**** separated string segments in SUBDD,SUBDD{1,2,3,4} in the rundeck
      call init_subdd(aDATE)
#endif
      call sys_flush(6)

C****
C**** MAIN LOOP
C****
      call gettime(tloopbegin)
      start9 = .false.
      if (istart==9) start9=.true.
      call startTimer('Main Loop')
      DO WHILE (Itime.lt.ItimeE)
#ifdef USE_SYSUSAGE
        call sysusage(0,1)
#endif


#if !defined( ADIABATIC ) || defined( CUBED_SPHERE)

c$$$         call test_save(__LINE__, itime)
C**** Every Ndisk Time Steps (DTsrc), starting with the first one,
C**** write restart information alternately onto 2 disk files
      if (Ndisk > 0) then
        IF (MOD(Itime-ItimeI,Ndisk).eq.0.or.start9) THEN
          start9 = .false.
          CALL RFINAL (IRAND)
          call set_param( "IRAND", IRAND, 'o' )
          call io_rsf(rsf_file_name(KDISK),Itime,iowrite,ioerr)
#if defined( USE_FVCORE )
          fv_fname='fv.'   ; write(fv_fname(4:4),'(i1)') kdisk
          fv_dfname='dfv.' ; write(fv_dfname(5:5),'(i1)') kdisk
          call Checkpoint(fv, clock, fv_fname, fv_dfname)
#endif
          if (AM_I_ROOT())
     *      WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)')
     *     '0Restart file written on fort.',KDISK,'Year',
     *     JYEAR,aMON,JDATE,', Hr',JHOUR,'  Internal clock time:',ITIME
          KDISK=3-KDISK
          CALL TIMER (NOW,MELSE)
        END IF
      end if
C**** THINGS THAT GET DONE AT THE BEGINNING OF EVERY DAY
      IF (MOD(Itime,NDAY).eq.0) THEN
C**** INITIALIZE SOME DIAG. ARRAYS AT THE BEGINNING OF SPECIFIED DAYS
        if (kradia.le.0) call daily_DIAG
C**** THINGS THAT GET DONE AT THE BEGINNING OF EVERY MONTH
        IF ( JDAY.eq.1+JDendOfM(Jmon-1) ) then
          write(aDATE(1:7),'(a3,I4.4)') aMON(1:3),Jyear
C**** THINGS THAT GET DONE AT THE BEGINNING OF EVERY ACC.PERIOD
          months=(Jyear-Jyear0)*JMperY + JMON-JMON0
          if ( months.ge.NMONAV ) then
            call reset_DIAG(0)
            if (Kvflxo.ne.0) then
              call closeunit( iu_VFLXO )
              call openunit('VFLXO'//aDATE(1:7),iu_VFLXO,.true.,.false.)
            end if
#ifndef CACHED_SUBDD
C**** reset sub-daily diag files
            call reset_subdd(aDATE)
#endif
          end if   !  beginning of acc.period
        END IF     !  beginning of month
      END IF       !  beginning of day
#ifdef CACHED_SUBDD
      call set_subdd_period()
#endif
C****
C**** INTEGRATE DYNAMIC TERMS (DIAGA AND DIAGB ARE CALLED FROM DYNAM)
C****
      if(Kradia>9) go to 100 ! to test daily/monthly procedures fast
      CALL CHECKT ('DYNAM0')
      if (kradia.le.0) then                   ! full model,kradia le 0
         MODD5D=MOD(Itime-ItimeI,NDA5D)

         IF (MODD5D.EQ.0) IDACC(ia_d5d)=IDACC(ia_d5d)+1
         IF (MODD5D.EQ.0) CALL DIAG5A (2,0)
         IF (MODD5D.EQ.0) CALL DIAGCA (1)

      PTOLD = P ! save for clouds
C**** Initialize pressure for mass fluxes used by tracers and Q
      PS (:,:)   = P(:,:)

C**** Initialise total energy (J/m^2)
      initialTotalEnergy = getTotalEnergy()

#ifdef SCM
      NSTEPSCM = ITIME-ITIMEI
      write(0,*) 'NSTEPSCM ',NSTEPSCM
      do L=1,LM
         SCM_SAVE_T(L) = T(I_TARG,J_TARG,L)
         SCM_SAVE_Q(L) = Q(I_TARG,J_TARG,L)
      enddo
c     do L=1,LM
c        write(iu_scm_prt,'(a13,i3,4(f9.3))')
c    &              'before dynam ',
c    &               L,T(I_TARG,J_TARG,L)*PK(L,I_TARG,J_TARG),
c    &               Q(I_TARG,J_TARG,L)*1000.0,
c    &               U(I_TARG,J_TARG,L),V(I_TARG,J_TARG,L)
c     enddo
#endif

! ADIABATIC
#endif
! ADIABATIC

#ifdef USE_SYSUSAGE
      call sysusage(1,1)
#endif
      call startTimer('Atm. Dynamics')

#ifndef USE_FVCORE

      CALL DYNAM()
#else

      ! Using FV instead
        IF (MOD(Itime-ItimeI,NDAA).eq.0) CALL DIAGA0

      call Run(fv, clock)

#ifndef CUBED_SPHERE
      CALL SDRAG (DTsrc)
#endif
        if (MOD(Itime-ItimeI,NDAA).eq.0) THEN
          call DIAGA
          call DIAGB
#ifdef CUBED_SPHERE
          call EPFLUX
#endif
        endif
#endif /* USE_FVCORE */

#if !defined( ADIABATIC ) || defined( CUBED_SPHERE)
C**** This fix adjusts thermal energy to conserve total energy TE=KE+PE
C**** Currently energy is put in uniformly weighted by mass
      finalTotalEnergy = getTotalEnergy()
      call addEnergyAsDiffuseHeat(finalTotalEnergy - initialTotalEnergy)
#ifndef CUBED_SPHERE
      call COMPUTE_DYNAM_AIJ_DIAGNOSTICS(PUA, PVA, DT)
#endif
#ifdef SCM
       do L=1,LM
          CONV(I_TARG,J_TARG,L) = SG_CONV(L)
       enddo
#endif
      SD_CLOUDS(:,:,:) = CONV(:,:,:)
      call COMPUTE_WSAVE
C**** Scale WM mixing ratios to conserve liquid water
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LS1-1
      DO J=J_0,J_1
      DO I=I_0,I_1
        WM(I,J,L)=WM(I,J,L)* (PTOLD(I,J)/P(I,J))
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO
      CALL QDYNAM  ! Advection of Q by integrated fluxes
      call stopTimer('Atm. Dynamics')
         CALL TIMER (NOW,MDYN)
#ifdef TRACERS_ON
      CALL TrDYNAM   ! tracer dynamics
#ifdef TRAC_ADV_CPU
         CALL TIMER (NOW,MTRADV)
#else
         CALL TIMER (NOW,MTRACE)
#endif
#endif
#ifdef USE_SYSUSAGE
         call sysusage(1,2)
#endif

#ifdef CACHED_SUBDD
      call get_subdd_vinterp_coeffs ! doing after dynamics
#endif

C****
C**** Calculate tropopause level and pressure
C****
      CALL CALC_TROP
C**** calculate some dynamic variables for the PBL
#ifndef SCM
      CALL PGRAD_PBL
#endif
C**** calculate zenith angle for current time step
      CALL CALC_ZENITH_ANGLE

         CALL CHECKT ('DYNAM ')
         CALL TIMER (NOW,MSURF)
         IF (MODD5D.EQ.0) CALL DIAG5A (7,NIdyn)
         IF (MODD5D.EQ.0) CALL DIAGCA (2)
         IF (MOD(Itime,NDAY/2).eq.0) CALL DIAG7A
C****
C**** INTEGRATE SOURCE TERMS
C****

c calculate KE before atmospheric column physics
         call calc_kea_3d(kea)

#ifdef CUBED_SPHERE
c GWDRAG, SDRAG considered as column physics so that their KE
c dissipation gets included in the KE->PE adjustment
      CALL GWDRAG
      CALL SDRAG (DTsrc)
#endif

         IDACC(ia_src)=IDACC(ia_src)+1
         MODD5S=MOD(Itime-ItimeI,NDA5S)
         IF (MODD5S.EQ.0) IDACC(ia_d5s)=IDACC(ia_d5s)+1
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAG5A (1,0)
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAGCA (1)

C**** FIRST CALL MELT_SI SO THAT TOO SMALL ICE FRACTIONS ARE REMOVED
C**** AND ICE FRACTION CAN THEN STAY CONSTANT UNTIL END OF TIMESTEP
      CALL MELT_SI
         CALL UPDTYPE
         CALL TIMER (NOW,MSURF)
C**** CONDENSATION, SUPER SATURATION AND MOIST CONVECTION
      CALL CONDSE
         CALL CHECKT ('CONDSE')
         CALL TIMER (NOW,MCNDS)
         IF (MODD5S.EQ.0) CALL DIAG5A (9,NIdyn)
         IF (MODD5S.EQ.0) CALL DIAGCA (3)
      end if                                  ! full model,kradia le 0
C**** RADIATION, SOLAR AND THERMAL
      MODRD=MOD(Itime-ItimeI,NRAD)
      if (kradia.le.0. or. MODRD.eq.0) then
         CALL RADIA
         if (kradia.le.0) CALL CHECKT ('RADIA ')
      end if
         CALL TIMER (NOW,MRAD)
      if (kradia.le.0) then                    ! full model,kradia le 0
         IF (MODD5S.EQ.0) CALL DIAG5A (11,NIdyn)
         IF (MODD5S.EQ.0) CALL DIAGCA (4)
C****
C**** SURFACE INTERACTION AND GROUND CALCULATION
C****
C**** NOTE THAT FLUXES ARE APPLIED IN TOP-DOWN ORDER SO THAT THE
C**** FLUXES FROM ONE MODULE CAN BE SUBSEQUENTLY APPLIED TO THAT BELOW
C****
C**** APPLY PRECIPITATION TO SEA/LAKE/LAND ICE
      call startTimer('Surface')
      call startTimer('Precip')
      CALL PRECIP_SI
      CALL PRECIP_LI
C**** APPLY PRECIPITATION AND RUNOFF TO LAKES/OCEANS
#ifdef IRRIGATION_ON
C**** CHECK FOR IRRIGATION POSSIBILITY
      CALL IRRIG_LK
#endif
      CALL PRECIP_LK
      CALL PRECIP_OC
         CALL TIMER (NOW,MSURF)
         CALL CHECKT ('PRECIP')
      call stopTimer('Precip')
#ifdef TRACERS_ON
C**** Calculate non-interactive tracer surface sources and sinks
         call set_tracer_2Dsource
         CALL TIMER (NOW,MTRACE)
#endif
C**** CALCULATE SURFACE FLUXES AND EARTH
      CALL SURFCE
      call stopTimer('Surface')
         CALL CHECKT ('SURFCE')
         CALL TIMER (NOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (5)
#ifdef CALCULATE_FLAMMABILITY
      call flammability_drv
#endif
C**** CALCULATE ICE DYNAMICS
      CALL DYNSI
C**** CALCULATE BASE ICE-OCEAN/LAKE FLUXES
      CALL UNDERICE
C**** APPLY SURFACE/BASE FLUXES TO SEA/LAKE ICE
      CALL GROUND_SI
C**** APPLY SURFACE FLUXES TO LAND ICE
      CALL GROUND_LI
         CALL CHECKT ('GRNDSI')
C**** APPLY FLUXES TO LAKES AND DETERMINE ICE FORMATION
      CALL GROUND_LK
         CALL CHECKT ('GRNDLK')
         CALL TIMER (NOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (6)
C**** CALCULATE RIVER RUNOFF FROM LAKE MASS
      CALL RIVERF
      CALL GROUND_E    ! diagnostic only - should be merged with EARTH
C**** APPLY FLUXES TO OCEAN, DO OCEAN DYNAMICS AND CALC. ICE FORMATION
#ifdef USE_SYSUSAGE
      call sysusage(2,1)
#endif
      call startTimer('OCEANS')
      CALL OCEANS
      call stopTimer('OCEANS')
#ifdef USE_SYSUSAGE
      call sysusage(2,2)
#endif
         CALL CHECKT ('OCEANS')
C**** APPLY ICE FORMED IN THE OCEAN/LAKES TO ICE VARIABLES
      CALL FORM_SI
         CALL CHECKT ('FORMSI')
C**** IF ATURB is used in rundeck then this is a dummy call
C**** CALCULATE DRY CONVECTION ABOVE PBL
      CALL ATM_DIFFUS (2,LM-1,dtsrc)
         CALL CHECKT ('DRYCNV')
         CALL TIMER (NOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (9)
C**** ADVECT ICE
      CALL ADVSI
      CALL ADVSI_DIAG ! needed to update qflux model, dummy otherwise
         CALL CHECKT ('ADVSI ')
C**** UPDATE DIAGNOSTIC TYPES
         CALL UPDTYPE
C**** ADD DISSIPATED KE FROM COLUMN PHYSICS CALCULATION BACK AS LOCAL HEAT
      CALL DISSIP ! uses kea calculated before column physics
         CALL CHECKT ('DISSIP')
         CALL TIMER (NOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (7)
         IF (MODD5S.EQ.0) CALL DIAG5A (12,NIdyn)

#ifdef CUBED_SPHERE
      IDACC(ia_filt)=IDACC(ia_filt)+1 ! prevent /0
#else
C**** SEA LEVEL PRESSURE FILTER
      IF (MFILTR.GT.0.AND.MOD(Itime-ItimeI,NFILTR).EQ.0) THEN
           IDACC(ia_filt)=IDACC(ia_filt)+1
           IF (MODD5S.NE.0) CALL DIAG5A (1,0)
           CALL DIAGCA (1)
           CALL FILTER
           CALL CHECKT ('FILTER')
           CALL TIMER (NOW,MDYN)
           CALL DIAG5A (14,NFILTR*NIdyn)
           CALL DIAGCA (8)
      END IF
#endif
#ifdef TRACERS_ON
#ifdef CUBED_SPHERE
! Reinitialize instantaneous consrv qtys (every timestep since
! DIAGTCA is called every timestep for 3D sources)
      CALL DIAGCA (1) ! was not called w/ SLP filter
#endif
C**** 3D Tracer sources and sinks
C**** Tracer gravitational settling for aerosols
      CALL TRGRAV
C**** Tracer radioactive decay (and possible source)
      CALL TDECAY
C**** Calculate 3D tracers sources and sinks

      call tracer_3Dsource
C**** Accumulate tracer distribution diagnostics
      CALL TRACEA
         CALL TIMER (NOW,MTRACE)
         CALL CHECKT ('T3DSRC')
#endif
      end if                                  ! full model,kradia le 0
C****
C**** WRITE SUB-DAILY DIAGNOSTICS EVERY NSUBDD hours
C****
#ifdef CACHED_SUBDD
      call accum_subdd
#else
      if (Nsubdd.ne.0) then
        call accSubdd
        if (mod(Itime+1,Nsubdd).eq.0) call get_subdd
      end if
#endif
#ifdef TRACERS_DUST
      call ahourly
#endif
C****
C**** UPDATE Internal MODEL TIME AND CALL DAILY IF REQUIRED
C****
  100 Itime=Itime+1                       ! DTsrc-steps since 1/1/Iyear1
      Jhour=MOD(Itime*24/NDAY,24)         ! Hour (0-23)
      Nstep=Nstep+NIdyn                   ! counts DT(dyn)-steps

      IF (MOD(Itime,NDAY).eq.0) THEN      ! NEW DAY
         call startTimer('Daily')
      if (kradia.gt.0) then               ! radiative forcing run
        CALL DAILY(.false.)
        if(Kradia<10)  CALL daily_RAD(.true.)
        if(Kradia<10)  call daily_EARTH(.false.)
        if(Kradia==10) CALL daily_OCEAN(.true.) ! to test OCLIM
        months=(Jyear-Jyear0)*JMperY + JMON-JMON0
      else                                ! full model, kradia le 0
#ifdef CACHED_SUBDD
        if(write_daily_files) then
          write(yyyymmdd,'(i4,i2.2,i2.2)') jyear,jmon,jdate
          filenm=yyyymmdd//'.subdd'//XLABEL(1:LRUNID)
          call write_subdd_accfile (filenm)
        endif
#endif
           CALL DIAG5A (1,0)
           CALL DIAGCA (1)
        CALL DAILY(.true.)                 ! end_of_day
        CALL daily_RAD(.true.)
        months=(Jyear-Jyear0)*JMperY + JMON-JMON0
           CALL TIMER (NOW,MELSE)

        call daily_LAKE
        call daily_EARTH(.true.)           ! end_of_day

        call daily_OCEAN(.true.)           ! end_of_day
        call daily_ICE
        call daily_LI
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
        call daily_tracer(.true.)
           CALL TIMER (NOW,MTRACE)
#endif
           CALL CHECKT ('DAILY ')
           CALL TIMER (NOW,MSURF)
           CALL DIAG5A (16,NDAY*NIdyn)
           CALL DIAGCA (10)
        call sys_flush(6)
      end if   ! kradia: full model (or rad.forcing run)
      CALL UPDTYPE
      call stopTimer('Daily')
      END IF   !  NEW DAY

#ifdef USE_FVCORE
      Call Compute_Tendencies(fv)
#endif

      if (kradia.le.0) then   ! full model
C****
C**** WRITE INFORMATION FOR OHT CALCULATION EVERY 24 HOURS
C****
      IF (Kvflxo.EQ.0.) OA(:,:,4:KOA)=0. ! to prepare for future saves
      IF (Kvflxo.NE.0.) THEN
         IF (MOD(Itime,NDAY).eq.0) THEN
c            call pack_data (grid, OA, OA_glob)
c            if (am_I_root()) call WRITEI8 (iu_vflxo,Itime,OA_glob,im*jm*koa)
           call writei8_parallel(grid,iu_vflxo,
     &          nameunit(iu_vflxo),oa,Itime)
C**** ZERO OUT INTEGRATED QUANTITIES
            OA(:,:,4:KOA)=0.
         ELSEIF (MOD(Itime,NDAY/2).eq.0) THEN
            call vflx_OCEAN
         END IF
         CALL TIMER (NOW,MELSE)
      END IF
C****
C**** CALL DIAGNOSTIC ROUTINES
C****
        call startTimer('Diagnostics')
#ifdef SCM
c*****call scm diagnostics every time step
      do L=1,LM
         dTtot(L) = T(I_TARG,J_TARG,L)-dTtot(L)
         dqtot(L) = Q(I_TARG,J_TARG,L)-dqtot(L)
      enddo
      call scm_diag
#endif

      IF (MOD(Itime-ItimeI,NDA4).EQ.0) CALL DIAG4A ! at hr 23 E-history
C**** PRINT CURRENT DIAGNOSTICS (INCLUDING THE INITIAL CONDITIONS)
      IF (NIPRNT.GT.0) THEN
        acc_period='PARTIAL      '
#ifdef NEW_IO
        filenm='PARTIAL.acc'//XLABEL(1:LRUNID)
        call io_rsf (filenm,Itime,iowrite_single,ioerr)
#endif
        call print_diags(1)
        NIPRNT=NIPRNT-1
        call set_param( "NIPRNT", NIPRNT, 'o' )
      END IF

      end if   ! full model ; kradia le 0

C**** THINGS TO DO BEFORE ZEROING OUT THE ACCUMULATING ARRAYS
C**** (after the end of a diagn. accumulation period)
      IF (MOD(Itime,NDAY).eq.0) then
      IF (months.ge.NMONAV .and. JDAY.eq.1+JDendOfM(JMON-1)) then

C**** PRINT DIAGNOSTIC TIME AVERAGED QUANTITIES
        call aPERIOD (JMON0,JYEAR0,months,1,0, aDATE(1:12),Ldate)
        acc_period=aDATE(1:12)
        WRITE (aDATE(8:14),'(A3,I4.4)') aMON(1:3),JYEAR
        if (kradia.le.0) call print_diags(0)
C**** SAVE ONE OR BOTH PARTS OF THE FINAL RESTART DATA SET
        IF (KCOPY.GT.0) THEN
C**** KCOPY > 0 : SAVE THE DIAGNOSTIC ACCUM ARRAYS IN SINGLE PRECISION
          monacc = 0
          do k=JMON0,JMON0+NMONAV-1
            m = k
            if(m.gt.12) m = m-12
            monacc(m) = 1
          end do
          filenm=aDATE(1:7)//'.acc'//XLABEL(1:LRUNID)
          call io_rsf (filenm,Itime,iowrite_single,ioerr)
#ifdef CACHED_SUBDD
          if(.not.write_daily_files) then
            filenm=aDATE(1:7)//'.subdd'//XLABEL(1:LRUNID)
            call write_subdd_accfile (filenm)
          endif
#endif
C**** KCOPY > 1 : ALSO SAVE THE RESTART INFORMATION
          IF (KCOPY.GT.1) THEN
            CALL RFINAL (IRAND)
            call set_param( "IRAND", IRAND, 'o' )
            filenm='1'//aDATE(8:14)//'.rsf'//XLABEL(1:LRUNID)
            call io_rsf(filenm,Itime,iowrite_mon,ioerr)
#if defined( USE_FVCORE )
            fv_fname  = '1'//aDATE(8:14)//'.fv'//XLABEL(1:LRUNID)
            fv_dfname = '1'//aDATE(8:14)//'.dfv'//XLABEL(1:LRUNID)
            call Checkpoint(fv, clock, fv_fname, fv_dfname)
#endif
          END IF
C**** KCOPY > 2 : ALSO SAVE THE OCEAN DATA TO INITIALIZE DEEP OCEAN RUNS
          IF (KCOPY.GT.2) THEN
            If (AM_I_ROOT())
     *           call openunit(aDATE(1:7)//'.oda'//XLABEL(1:LRUNID)
     *           ,iu_ODA,.true.,.false.)
            call io_oda(iu_ODA,Itime,iowrite,ioerr)
            IF (AM_I_ROOT()) call closeunit(iu_ODA)
          END IF
        END IF

C**** PRINT AND ZERO OUT THE TIMING NUMBERS
        CALL TIMER (NOW,MDIAG)
        CALL SUMXPE(TIMING, TIMING_glob, increment=.true.)
        if (am_i_root()) then
        TOTALT=SUM(TIMING_glob(1:NTIMEACC))  ! over all processors
        DO M=1,NTIMEACC
          PERCENT(M) = 100d0*TIMING_glob(M)/(TOTALT+.00001)
        END DO
        TOTALT=SUM(TIMING(1:NTIMEACC)) ! on the root processor
        TOTALT=TOTALT/60.       ! seconds -> minutes
        DTIME = NDAY*TOTALT/(Itime-Itime0)        ! minutes/day
        WRITE (6,'(/A,F7.2,A,/(8(A13,F5.1/))//)')
     *   '0TIME',DTIME,'(MINUTES) ',(TIMESTR(M),PERCENT(M),M=1,NTIMEACC)
        end if
        TIMING = 0
        START= NOW

      END IF ! beginning of accumulation period
      END IF ! beginning of day

C**** CPU TIME FOR CALLING DIAGNOSTICS
      call stopTimer('Diagnostics')
      CALL TIMER (NOW,MDIAG)
C**** TEST FOR TERMINATION OF RUN
ccc
      IF (MOD(Itime,Nssw).eq.0) then
       IF (AM_I_ROOT()) then
        flg_go = '__STOP__'     ! stop if flagGoStop if missing
        iflag=0
        open(3,file='flagGoStop',form='FORMATTED',status='OLD',err=210)
        read (3,'(A8)',end=210) flg_go
        close (3)
 210    continue
        IF (flg_go .eq. '___GO___') iflag=1
        call ESMF_BCAST( grid, iflag)
       else
        call ESMF_BCAST( grid, iflag)
        if (iflag .eq. 1) flg_go = '___GO___'
        if (iflag .eq. 0) flg_go = '__STOP__'
       end if
      endif
      IF (flg_go.ne.'___GO___' .or. stop_on) THEN
C**** Flag to continue run has been turned off
         WRITE (6,'("0Flag to continue run has been turned off.")')
         EXIT
      END IF

c$$$      call test_save(__LINE__, itime-1)

! ADIABATIC
#else
      Itime=Itime+1                       ! DTsrc-steps since 1/1/Iyear1
#endif
! ADIABATIC
#ifdef USE_SYSUSAGE
      call sysusage(0,2)
#endif
      END DO
      call stopTimer('Main Loop')

      call gettime(tloopend)
      if (AM_I_ROOT())
     *     write(*,*) "Time spent in the main loop in seconds:",
     *     tloopend-tloopbegin
C****
C**** END OF MAIN LOOP
C****

#if !defined( ADIABATIC ) || defined( CUBED_SPHERE)
C**** ALWAYS PRINT OUT RSF FILE WHEN EXITING
      CALL RFINAL (IRAND)
      call set_param( "IRAND", IRAND, 'o' )
      call io_rsf(rsf_file_name(KDISK),Itime,iowrite,ioerr)
#endif

#ifdef USE_FVCORE
         fv_fname='fv.' ; write(fv_fname(4:4),'(i1)') kdisk
         fv_dfname='dfv.' ; write(fv_dfname(5:5),'(i1)') kdisk
         call Finalize(fv, clock, fv_fname, fv_dfname)
#endif

#if !defined( ADIABATIC ) || defined( CUBED_SPHERE)
      if (AM_I_ROOT()) then
      WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)')
     *  '0Restart file written on fort.',KDISK,'Year',JYEAR,
     *     aMON,JDATE,', Hr',JHOUR,'  Internal clock time:',ITIME
      end if

C**** RUN TERMINATED BECAUSE IT REACHED TAUE (OR SS6 WAS TURNED ON)
#ifdef SCM
      call closeunit(iu_scm_prt)
      call closeunit(iu_scm_diag)
#endif
#endif

#ifdef USE_SYSUSAGE
      do i_su=0,max_su
        call sysusage(i_su,3)
      enddo
#endif

      IF (AM_I_ROOT())
     *   WRITE (6,'(/////4(1X,33("****")/)//,A,I8
     *             ///4(1X,33("****")/))')
     *  ' PROGRAM TERMINATED NORMALLY - Internal clock time:',ITIME

      IF (Itime.ge.ItimeE) then
         call reportProfile((itimee-itimei)*dtSRC)
         if (AM_I_ROOT()) call print_unused_param(6)
         CALL stop_model (
     &     'Terminated normally (reached maximum time)',13)
      END IF

#ifndef ADIABATIC
#ifdef BLK_2MOM
      ldummy=cleanup_bulk2m_driver()
#endif
#endif

      CALL stop_model ('Run stopped with sswE',12)  ! voluntary stop
#ifdef USE_MPP
      call fms_end( )
#endif

      contains

      subroutine initializeDefaultTimers()
      use TimerPackage_mod

      call initialize()
      call addTimer('Main Loop')
      call addTimer(' Atm. Dynamics')
      call addTimer(' MELT_SI()')
      call addTimer(' CONDSE()')
      call addTimer(' RADIA()')
      call addTimer(' Surface')
      call addTimer('  Precip')
      call addTimer('  SURFCE()')
      call addTimer(' DYNSI()')
      call addTimer(' UNDERICE()')
      call addTimer(' GROUND_SI()')
      call addTimer(' GROUND_LI()')
      call addTimer(' GROUND_LK()')
      call addTimer(' RIVERF()')
      call addTimer(' OCEANS')
      call addTimer(' Daily')
      call addTimer(' Diagnostics')

      end subroutine initializeDefaultTimers

      subroutine reportProfile(elapsedTimeInSeconds)
      use TimerPackage_mod
      use DOMAIN_DECOMP_1D, only: AM_I_ROOT

#if  defined(USE_ESMF) || defined(USE_MPP)
      include 'mpif.h'
#endif
      real*8, intent(in) :: elapsedTimeInSeconds
      character(len=MAX_RECORD_LENGTH), pointer :: lines(:)
      integer :: i

      integer, parameter :: SECONDS_PER_MINUTE = 60
      integer, parameter :: MINUTES_PER_DAY = 24*60
      type (ProfileReport_type) :: report
      type (ReportColumn_type) :: column
      real*8 :: totalTime

      call finalize()
      report = newReport()
      call addColumn(report, newColumn(NAME_COLUMN, fieldWidth=20))

      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth=7)
      totalTime = getInclusiveTime(getTimer('main'))
      call setPrecision(column, 2)
      call setScale(column, 100/totalTime, ' %')
      call addColumn(report, column)

      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth=11)
      call setPrecision(column, 5)
      call setScale(column, MINUTES_PER_DAY/elapsedTimeInSeconds,
     &     'min/day')
      call addColumn(report, column)

      column = newColumn(EXCLUSIVE_TIME_COLUMN, fieldWidth=11)
      call setPrecision(column, 5)
      call setScale(column, MINUTES_PER_DAY/elapsedTimeInSeconds,
     &     'min/day')
      call addColumn(report, column)

      call addColumn(report, newColumn(TRIP_COUNTS_COLUMN,fieldWidth=6))

      column = newColumn(MAXIMUM_TIME_COLUMN, fieldWidth=10)
      call setPrecision(column, 6)
      call addColumn(report, column)

      call addColumn(report, newColumn(MAX_PROCESS_COLUMN,fieldWidth=4))

      column = newColumn(MINIMUM_TIME_COLUMN, fieldWidth=10)
      call setPrecision(column, 6)
      call addColumn(report, column)
      call addColumn(report, newColumn(MIN_PROCESS_COLUMN,fieldWidth=4))

      column = newColumn(AVERAGE_TIME_COLUMN, fieldWidth=10)
      call setPrecision(column, 6)
      call addColumn(report, column)

#if  defined(USE_ESMF) || defined(USE_MPP)
      lines => generateParallelReport(report, MPI_COMM_WORLD)
#else
      lines => generateReport(report)
#endif

      if (AM_I_ROOT()) then
         do i = 1, size(lines)
            write(*,'(a)') trim(lines(i))
         end do
      end if
      deallocate(lines)
      call delete(report)

      end subroutine reportProfile

      end program GISS_modelE


      subroutine sig_stop_model
      USE MODEL_COM, only : stop_on
      implicit none
      stop_on = .true.
      end subroutine sig_stop_model


      subroutine init_Model
!@sum This program reads most of parameters from the database (DB)
!@+   get_param( "A", X ) reads parameter A into variable X
!@+   if "A" is not in the database, it will generate an error
!@+   message and stop
!@+   sync_param( "B", Y ) reads parameter B into variable Y
!@+   if "B" is not in the database, then Y is unchanged and its
!@+   value is saved in the database as "B" (here sync = synchronize)
      USE MODEL_COM, only : JM,LM,NIPRNT,MFILTR,NFILTR,NRAD
     *     ,NDASF,NDA4,NDA5S,NDA5K,NDA5D,NDAA,Kvflxo,kradia
     *     ,NMONAV,Ndisk,Nssw,KCOPY,KOCEAN,NIsurf,iyear1
     $     ,LS1,IRAND,ItimeI,PSTRAT,UOdrag,USE_UNR_DRAG
     $     ,X_SDRAG,C_SDRAG,LSDRAG,P_SDRAG,LPSDRAG,PP_SDRAG,ang_sdrag
     $     ,P_CSDRAG,CSDRAGL,Wc_Jdrag,wmax,VSDRAGL,COUPLED_CHEM,dt
     *     ,DT_XUfilter,DT_XVfilter,DT_YVfilter,DT_YUfilter,QUVfilter
     &     ,do_polefix,pednl00,pmidl00,ij_debug
      USE RAD_COM, only : variable_orb_par,orb_par_year_bp,orb_par
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      USE PARAM
      implicit none
      INTEGER L,LCSDRAG

C**** Rundeck parameters:
      call sync_param( "NMONAV", NMONAV )
      call sync_param( "NIPRNT", NIPRNT )
      call sync_param( "DT_XVfilter", DT_XVfilter )
      call sync_param( "DT_XUfilter", DT_XUfilter )
      call sync_param( "DT_YVfilter", DT_YVfilter )
      call sync_param( "DT_YUfilter", DT_YUfilter )
      call sync_param( "MFILTR", MFILTR )
      call sync_param( "USE_UNR_DRAG", USE_UNR_DRAG )
      call sync_param( "X_SDRAG", X_SDRAG, 2 )
      call sync_param( "C_SDRAG", C_SDRAG )
      call sync_param( "P_CSDRAG", P_CSDRAG )
      call sync_param( "P_SDRAG", P_SDRAG )
      call sync_param( "PP_SDRAG", PP_SDRAG )
      call sync_param( "ANG_SDRAG", ANG_SDRAG )
      call sync_param( "Wc_Jdrag", Wc_Jdrag )
      call sync_param( "VSDRAGL", VSDRAGL, lm-ls1+1 )
      call sync_param( "wmax", wmax )
      call sync_param( "do_polefix", do_polefix )
      call sync_param( "NDASF", NDASF )
      call sync_param( "NDA4", NDA4 ) !!
      call sync_param( "NDA5S", NDA5S ) !!
      call sync_param( "NDA5K", NDA5K ) !!
      call sync_param( "NDA5D", NDA5D ) !!
      call sync_param( "NDAA", NDAA ) !!
      call sync_param( "NFILTR", NFILTR ) !!
      call sync_param( "NRAD", NRAD ) !!
      call sync_param( "Kvflxo", Kvflxo ) !!
      call sync_param( "Ndisk", Ndisk )
      call sync_param( "Nssw", Nssw )
      call sync_param( "KCOPY", KCOPY )
      call sync_param( "KOCEAN", KOCEAN )
      call sync_param( "KRADIA", KRADIA )
      call sync_param( "NIsurf", NIsurf )
      call sync_param( "UOdrag", UOdrag )
      call sync_param( "IRAND", IRAND )
      call sync_param( "COUPLED_CHEM", COUPLED_CHEM )
      call sync_param( "ij_debug",ij_debug , 2)
      call sync_param( "variable_orb_par", variable_orb_par )
      call sync_param( "orb_par_year_bp", orb_par_year_bp )
      call sync_param( "orb_par", orb_par, 3 )

C**** Parameters derived from Rundeck parameters:

C**** Calculate levels for application of SDRAG: LSDRAG,LPSDRAG->LM i.e.
C**** all levels above and including P_SDRAG mb (PP_SDRAG near poles)
C**** If P is the edge between 2 levels, take the higher level.
C**** Also find CSDRAGL, the coefficients of C_Sdrag as a function of L

      LSDRAG=LM ; LPSDRAG=LM ; LCSDRAG=LM ; CSDRAGL=C_SDRAG
      DO L=1,LM
        IF (PEDNL00(L+1)-1d-5.lt.P_SDRAG .and.
     *      PEDNL00(L)  +1d-5.gt.P_SDRAG)         LSDRAG=L
        IF (PEDNL00(L+1)-1d-5.lt.PP_SDRAG .and.
     *      PEDNL00(L)  +1d-5.gt.PP_SDRAG)        LPSDRAG=L
        IF (PEDNL00(L+1)-1d-5.lt.P_CSDRAG .and.
     *      PEDNL00(L)  +1d-5.gt.P_CSDRAG)        LCSDRAG=L
      END DO
      DO L=LCSDRAG,LSDRAG-1
         CSDRAGL(L) = C_SDRAG + max( 0.d0 , (X_SDRAG(1)-C_SDRAG) *
     *     LOG(P_CSDRAG/(PMIDL00(L))) / LOG(P_CSDRAG/P_SDRAG) )
      END DO
      if (AM_I_ROOT()) then
         WRITE(6,*) "Levels for  LSDRAG =",LSDRAG ,"->",LM
         WRITE(6,*) "Levels for LPSDRAG =",LPSDRAG,"->",LM," near poles"
         WRITE(6,*) "C_SDRAG coefficients:",CSDRAGL(LS1:LSDRAG-1)
      end if

C**** Determine if FLTRUV is called.
      QUVfilter = .false.
      if (DT_XUfilter>0. .or. DT_XVfilter>0. .or.
     *    DT_YUfilter>0. .or. DT_YVfilter>0.)  QUVfilter = .true.
      if (QUVfilter) then
         if (DT_XUfilter > 0. .and. DT_XUfilter < DT) then
             DT_XUfilter = DT
             WRITE(6,*) "DT_XUfilter too small; reset to :",DT_XUfilter
         end if
         if (DT_XVfilter > 0. .and. DT_XVfilter < DT) then
             DT_XVfilter = DT
             WRITE(6,*) "DT_XVfilter too small; reset to :",DT_XVfilter
         end if
         if (DT_YUfilter > 0. .and. DT_YUfilter < DT) then
             DT_YUfilter = DT
             WRITE(6,*) "DT_YUfilter too small; reset to :",DT_YUfilter
         end if
         if (DT_YVfilter > 0. .and. DT_YVfilter < DT) then
             DT_YVfilter = DT
             WRITE(6,*) "DT_YVfilter too small; reset to :",DT_YVfilter
         end if
      end if
c Warn if polar fixes requested for a model not having a half polar box
c     if(do_polefix.eq.1 .and. jm.ne.46) then
c        do_polefix = 0
c        write(6,*) 'Polar fixes are currently applicable only to'//
c    &           'models having a half polar box; no fixes applied'
c     endif
      RETURN
C****
      end subroutine init_Model

#ifdef USE_FVCORE
      SUBROUTINE INPUT (istart,ifile,coldRestart,clock)
#else
      SUBROUTINE INPUT (istart,ifile,coldRestart)
#endif
C****
C**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE
C**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS
C****
      USE FILEMANAGER, only : openunit,closeunit,nameunit
      USE TIMINGS, only : timing,ntimeacc
      USE PARAM
      !USE PARSER
      USE CONSTANT, only : grav,kapa,sday,by3,twopi
      USE MODEL_COM, only : im,jm,lm,wm,u,v,t,p,q,fearth0,fland
     *     ,focean,flake0,flice,hlake,zatmo,plbot,sig,dsig,sige,kradia
     *     ,bydsig,xlabel,lrunid,nmonav,qcheck,irand,ptop
     *     ,nisurf,nidyn,nday,dt,dtsrc,kdisk,jmon0,jyear0
     *     ,iyear1,itime,itimei,itimee
     *     ,ls1,psfmpt,pstrat,idacc,jyear,jmon,jday,jdate,jhour
     *     ,aMONTH,jdendofm,jdpery,aMON,aMON0,ioread,irerun
     *     ,ioread_single,irsfic,irsficnt,iowrite_single,ioreadnt
     *     ,irsficno,mdyn,mcnds,mrad,msurf,mdiag,melse,Itime0,Jdate0
     *     ,Jhour0,rsf_file_name,lm_req
     *     ,pl00,aml00,pednl00,pdsigl00,pmidl00,byaml00,coupled_chem
     *     ,USE_UNR_DRAG
#ifdef BLK_2MOM
     * ,wmice
#endif

#ifdef SCM
      USE MODEL_COM, only : I_TARG,J_TARG
      USE SCMCOM, only : iu_scm_prt
#endif
      USE SOMTQ_COM, only : mz,tmom,qmom
#ifdef CUBED_SPHERE
       use GEOM, only : geom_cs,imaxj
#else
       USE GEOM, only : geom_b,imaxj
#endif
      USE DIAG_ZONAL, only : imlon
      USE RANDOM
      USE RAD_COM, only : rqt,cloud_rad_forc
      USE DYNAMICS, only : pk,pmid,pedn,ualij,valij
      USE CLOUDS_COM, only : ttold,qtold,svlhx,rhsav,cldsav
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM,only: MTRACE,NTM,TRNAME,daily_z
#ifdef TRACERS_SPECIAL_Shindell
     *     ,mchem
#endif
#ifdef TRAC_ADV_CPU
      USE TRACER_COM,only: MTRADV
#endif
#endif
#ifdef TRACERS_AMP
      USE AERO_CONFIG
      USE AERO_COAG
      USE AERO_INIT
      USE AERO_SETUP
      USE AERO_SUBS
      USE AERO_NPF
      USE AERO_DIAM
      USE AMP_AEROSOL
#endif
      USE DIAG_COM, only : acc_period,monacc
     &  ,hr_in_day,iwrite,jwrite,itwrite,kdiag,qdiag,qdiag_ratios,oa
      USE PBLCOM
     &     , only : wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,ustar_pbl
     &  ,egcm,w2gcm,tgvavg,qgavg
      USE LAKES_COM, only : flake
      USE GHY_COM, only : fearth
      USE SOIL_DRV, only: init_LSM
      USE DOMAIN_DECOMP_1D, only : HERE
      USE DOMAIN_DECOMP_ATM, only : grid, GET, AM_I_ROOT
      USE DOMAIN_DECOMP_ATM, only : HALO_UPDATE,READT_PARALLEL
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: init_app_clock
      USE CONSTANT, only : hrday
      USE ESMF_MOD, only: ESMF_Clock
#endif
#ifndef CUBED_SPHERE
      USE ATMDYN, only : init_ATMDYN
#endif
#ifdef IRRIGATION_ON
      use irrigate_crop, only : init_irrigate
#endif
cddd#ifdef USE_ENT
cddd      USE ENT_DRV, only : init_module_ent
cddd#endif
      IMPLICIT NONE
      CHARACTER(*) :: ifile
      logical, intent(in) :: coldRestart
!@var iu_AIC,iu_TOPO unit numbers for input files
      INTEGER iu_AIC,iu_TOPO,iu_IFILE
!@var num_acc_files number of acc files for diag postprocessing
      INTEGER I,J,L,K,LID1,LID2,IM1,NOFF,ioerr,num_acc_files
!@dbparam init_topog_related : set = 1 if IC and topography are incompatible
      integer :: init_topog_related = 0
!@dbparam do_IC_fixups : set = 1 if IC are to be checked/corrected
      integer :: do_IC_fixups = 0
!@nlparam HOURI,DATEI,MONTHI,YEARI        start of model run
!@nlparam TIMEE,HOURE,DATEE,MONTHE,YEARE,IHOURE   end of model run
!@var  IHRI,IHOURE start and end of run in hours (from 1/1/IYEAR1 hr 0)
      INTEGER ::   HOURI=0 , DATEI=1, MONTHI=1, YEARI=-1, IHRI=-1,
     *    TIMEE=-1,HOURE=0 , DATEE=1, MONTHE=1, YEARE=-1, IHOURE=-1,
!@nlparam ISTART  postprocessing(-1)/start(1-8)/restart(>8)  option
!@nlparam IRANDI  random number seed to perturb init.state (if>0)
     *             ISTART, IRANDI=0
      REAL*8 TIJL,CDM,TEMP,X
      INTEGER ItimeX,IhrX, LMR
#ifdef SCM
      REAL*4 SCMDATEB,SCMDATEE
      INTEGER ISCMMONTHI
#endif

!@ egcm_init_max maximum initial vaule of egcm
      real*8, parameter :: egcm_init_max=0.5d0

      LOGICAL :: redoGH = .FALSE.,iniPBL = .FALSE., inilake = .FALSE.,
     &           iniSNOW = .FALSE.  ! true = restart from "no snow" rsf
     &           ,iniOCEAN = .FALSE.
cddd#ifdef USE_ENT
cddd     &     ,iniENT = .FALSE.
cddd#endif
#ifdef USE_FVCORE
      type (ESMF_Clock) :: clock
      integer :: minti,minte
      character(len=1) :: suffix
#endif
      CHARACTER NLREC*80,filenm*100,RLABEL*132
      NAMELIST/INPUTZ/ ISTART,IRANDI
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK,QDIAG,KDIAG,QDIAG_RATIOS
     *     ,IHOURE, TIMEE,HOURE,DATEE,MONTHE,YEARE,IYEAR1
C****    List of parameters that are disregarded at restarts
     *     ,        HOURI,DATEI,MONTHI,YEARI
      NAMELIST/INPUTZ_cold/ ISTART,IRANDI
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK,QDIAG,KDIAG,QDIAG_RATIOS
     *     ,IHOURE, TIMEE,HOURE,DATEE,MONTHE,YEARE,IYEAR1
C****    List of parameters that are disregarded at restarts
     *     ,        HOURI,DATEI,MONTHI,YEARI
      integer ISTART_kradia, istart_fixup
      character*132 :: bufs
      integer :: kdisk_sv  ! for delayed reads from checkpoint file

      integer :: nij_before_j0,nij_after_j1,nij_after_i1
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
CCCC      INTEGER :: stdin ! used to read 'I' file
      integer, parameter :: MAXLEN_RUNID = 32
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

#ifdef USE_ESMF
      write(6,*) 'mpi-zone',J_0,' - ',J_1
#endif

C****
C**** Default setting for ISTART : restart from latest save-file (10)
C****
      ISTART=10
      num_acc_files=0
C****
C**** Set dependent vertical resolution variables
C****
      SIGE(:) = (PLbot(:)-PTOP)/PSFMPT
      SIG(:)  = (sige(1:lm)+sige(2:lm+1))*0.5d0
      DSIG(:) =  sige(1:lm)-sige(2:lm+1)
      byDSIG  =  1./DSIG
#ifdef CUBED_SPHERE
      call geom_cs
#else
C**** CALCULATE SPHERICAL GEOMETRY
      CALL GEOM_B
#endif
C**** Calculate default vertical arrays (including rad. eq. layers)
      LMR=LM+LM_REQ
      CALL CALC_VERT_AMP(PSFMPT,LMR,PL00,AML00,PDSIGL00,PEDNL00,PMIDL00)
      BYAML00(:)=1./AML00(:)
C****
C**** default settings for prog. variables etc
C****
      TEMP=250.
      TSAVG(:,:)=TEMP
      U(:,:,:)=0.
      V(:,:,:)=0.
      T(:,:,:)=TEMP  ! will be changed to pot.temp later
      Q(:,:,:)=3.D-6
      P(:,:)=PSFMPT
C**** Advection terms for first and second order moments
      TMOM(:,:,:,:)=0.
      QMOM(:,:,:,:)=0.
C**** Auxiliary clouds arrays
      RHSAV (:,:,:)=.85d0
      CLDSAV(:,:,:)=0.
      SVLHX (:,:,:)=0.
      WM    (:,:,:)=0.
#ifdef BLK_2MOM
      WMICE (:,:,:)=0.
#endif

C****    Ocean info saved for ocean heat transport calculations
         OA = 0.
C**** All diagn. are enabled unless KDIAG is changed in the rundeck
      KDIAG(1:12)=0
      KDIAG(13)=9
C**** Set global default timing descriptions
C**** Other speciality descriptions can be added/used locally
      NTIMEACC = 0
      CALL SET_TIMER("ATMOS. DYNAM",MDYN)
      CALL SET_TIMER("CONDENSATION",MCNDS)
      CALL SET_TIMER("   RADIATION",MRAD)
      CALL SET_TIMER("     SURFACE",MSURF)
      CALL SET_TIMER(" DIAGNOSTICS",MDIAG)
      CALL SET_TIMER("       OTHER",MELSE)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      CALL SET_TIMER("     TRACERS",MTRACE)
#ifdef TRAC_ADV_CPU
      CALL SET_TIMER(" TRACER ADV.",MTRADV)
#endif
#endif
#ifdef TRACERS_SPECIAL_Shindell
      CALL SET_TIMER("   CHEMISTRY",MCHEM)
#endif
C****
C**** Set some documentary parameters in the database
C****
      call set_param("IM",IM)
      call set_param("JM",JM)
      call set_param("LM",LM)
      call set_param("LS1",LS1)
      call set_param("PLBOT",Plbot,LM+1)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      call set_param("NTM",NTM)
      call set_param("TRNAME",TRNAME,ntm)
#endif
C****
C**** Print Header and Label (2 lines) from rundeck
C****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      if (AM_I_ROOT()) WRITE (6,'(A,40X,A/)') '0','GISS CLIMATE MODEL'
      READ(iu_IFILE,'(A80)') XLABEL(1:80),NLREC
      NOFF=0
      IF (XLABEL(73:80).EQ.'        ') NOFF=8   ! for 72-column rundecks
      XLABEL(81-NOFF:132)=NLREC(1:52+NOFF)
      if (AM_I_ROOT()) WRITE (6,'(A,A/)') '0',XLABEL
      RLABEL = XLABEL !@var RLABEL rundeck-label

C****
C**** Print preprocessing options (if any are defined)
C****
      IF(AM_I_ROOT()) call print_and_check_PPopts

C****
C**** Read parameters from the rundeck to database and namelist
C****
      !call parse_params(iu_IFILE)
      ! skip "&&PARAMETERS" section
      do
        read( iu_IFILE, *, err=910, end=910 ) bufs
        if ( bufs == '&&END_PARAMETERS' ) exit
      enddo

      READ (iu_IFILE,NML=INPUTZ,ERR=900)
      if (coldRestart) READ (iu_IFILE,NML=INPUTZ_cold,ERR=900)
      call closeunit(iu_IFILE)

C**** Get those parameters which are needed in this subroutine
      if(is_set_param("DTsrc"))  call get_param( "DTsrc", DTsrc )
      if(is_set_param("DT"))     call get_param( "DT", DT )
      if(is_set_param("NIsurf")) call get_param( "NIsurf", NIsurf ) !
      if(is_set_param("IRAND"))  call get_param( "IRAND", IRAND )
      if(is_set_param("NMONAV")) call get_param( "NMONAV", NMONAV )
      if(is_set_param("Kradia")) call get_param( "Kradia", Kradia )
      if(is_set_param("coupled_chem"))
     &call get_param( "coupled_chem", coupled_chem )

C***********************************************************************
C****                                                               ****
C****        Post-process one or more ACC-files : ISTART < 1        ****
C****                                                               ****
C***********************************************************************
      if (istart.le.0) then
        call reset_diag(1)
        monacc = 0
        do
          call nextarg(filenm, 0)
          if ( filenm == "" ) exit ! end of args
          call io_rsf(filenm,itime,ioread_single,ioerr)
          num_acc_files = num_acc_files + 1
        end do
        GO TO 500
      end if

      if (istart.ge.9 .or. Kradia.gt.0) go to 400
C***********************************************************************
C****                                                               ****
C****                  INITIAL STARTS - ISTART: 1 to 8              ****
C****                                                               ****
C****   Current settings: 1 - from defaults                         ****
C****                     2 - from observed data                    ****
C****                     3 - so far unused                         ****
C****                     4 - from coupled model M-file - reset ocn ****
C****                     5 - tracer run from M-file w/o tracers    ****
C****                     6 - pred.ocn run from M-file w/o ocn data ****
C****                     7 - from mod. II' M-file - reset snow/ocn ****
C****                     8 - from current model M-file - no resets ****
C****                                                               ****
C***********************************************************************
C****
C**** Set quantities that are derived from the namelist parameters
C****
!@var NDAY=(1 day)/DTsrc : even integer; adjust DTsrc later if necessary
      NDAY = 2*NINT(.5*SDAY/DTsrc)

C**** Get Start Time; at least YearI HAS to be specified in the rundeck
      IF (YearI.lt.0) then
        IF (AM_I_ROOT())
     *   WRITE(6,*) 'Please choose a proper start year yearI, not',yearI
        call stop_model('INPUT: yearI not provided',255)
      END IF
      IF (Iyear1.lt.0) Iyear1 = yearI
      IhrI = HourI +
     +  HR_IN_DAY*(dateI-1 + JDendofM(monthI-1) + JDperY*(yearI-Iyear1))
      ITimeI = IhrI*NDAY/HR_IN_DAY ! internal clock counts DTsrc-steps
      Itime=ItimeI
      IF (IhrI.lt.0) then
        IF (AM_I_ROOT())
     *  WRITE(6,*) 'Improper start time OR Iyear1=',Iyear1,' > yearI;',
     *  ' yearI,monthI,dateI,hourI=',yearI,monthI,dateI,hourI
        call stop_model(
     &       'INPUT: Improper start date or base year Iyear1',255)
      END IF
C**** Check the vertical layering defined in RES_ (is sige(ls1)=0 ?)
      IF (SIGE(LS1).ne.0.) then
        if (AM_I_ROOT())
     *       write(6,*) 'bad vertical layering: ls1,sige(ls1)',
     &       ls1,sige(ls1)
        call stop_model('INPUT: ls1 incorrectly set in RES_',255)
      END IF
C****
C**** Get Ground conditions from a separate file - ISTART=1,2
C****

      IF (ISTART.LE.2) THEN
    
C**** Set flag to initialise topography-related variables
        init_topog_related = 1
cddd#ifdef USE_ENT
cddd        iniENT = .TRUE.
cddd#endif
        if (istart.eq.1) redogh=.true.

C**** Read in ground initial conditions
        call read_ground_ic() ! code moved to IORSF

      END IF

C****
C**** Get primary Atmospheric data from NMC tapes - ISTART=2
C****
      IF (ISTART.EQ.2) THEN
C**** Use title of first record to get the date and make sure  ???
C**** it is consistent with IHRI (at least equal mod 8760)     ???
C****            not yet implemented but could easily be done  ???

C**** open atmospheric initial conditions file
        call openunit("AIC",iu_AIC,.true.,.true.)

        XLABEL(1:80)='Observed atmospheric data from NMC tape'
Csoon   READ (iu_AIC) XLABEL(1:80)

        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),P,1) ! Psurf
        DO J=J_0,J_1
          DO I=I_0,I_1
            P(I,J)=P(I,J)-PTOP                        ! Psurf -> P
          END DO
        END DO
        DO L=1,LM
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),U(:,:,L),1) ! U
        END DO
        DO L=1,LM
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),V(:,:,L),1) ! V
        END DO
        DO L=1,LM
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),T(:,:,L),1) ! Temperature
        END DO
        DO L=1,LM  ! alternatively, only read in L=1,LS1 ; skip rest
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),Q(:,:,L),1) ! Q
        END DO
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),TSAVG,1)  ! Tsurf

C**** Close "AIC"
        call closeunit(iu_AIC)

      END IF

C****
C**** Derive other data from primary data if necessary - ISTART=1,2
C****                                                    currently
      IF (ISTART.LE.2) THEN

#if defined(SCM) || defined(CUBED_SPHERE)
c in these cases, assume input U/V are on the A grid
        DO J=J_0,J_1
        DO I=I_0,I_1
          ualij(:,i,j) = u(i,j,:)
          valij(:,i,j) = v(i,j,:)
        END DO
        END DO
#else
c assume input U/V are on the B grid.  Need to calculate A-grid winds.
        call recalc_agrid_uv
c the latlon version of recalc_agrid_uv does not fill the poles.
c replicate polar data to avoid compiler traps in INPUT only.
        if(have_south_pole) then
          ualij(1,2:im,1) = ualij(1,1,1)
          valij(1,2:im,1) = valij(1,1,1)
        endif
        if(have_north_pole) then
          ualij(1,2:im,jm) = ualij(1,1,jm)
          valij(1,2:im,jm) = valij(1,1,jm)
        endif
#endif

        do j=j_0,j_1
        do i=i_0,i_1
          usavg(i,j) = ualij(1,i,j)
          vsavg(i,j) = valij(1,i,j)
          wsavg(i,j) = sqrt(usavg(i,j)**2 + vsavg(i,j)**2)
        enddo
        enddo

        CDM=.001d0

#ifdef SCM
c      enter SCM part of INPUT
       call sync_param( "I_TARG",I_TARG)
       call sync_param( "J_TARG",J_TARG)
       write(0,*) 'I/J Targets set ',I_TARG,J_TARG
       write(iu_scm_prt,*) 'I/J Targets set ',I_TARG,J_TARG

c      write(iu_scm_prt,*) 'before scm inputs L u v '
c      do L=1,LM
c         write(iu_scm_prt,'(a6,i5,2(f9.3))') 'l u v ',l,
c    &         u(i_targ,j_targ,l),v(i_targ,j_targ,l)
c      enddo

       if ((I_TARG.lt.1 .or. I_TARG.gt. 144) .or.
     &      (J_TARG.lt.2 .or. J_TARG.gt.89)) then
             write(iu_scm_prt,*)
     &             'Invalid grid coordinates for selected box ',
     &             I_TARG,J_TARG
             STOP 100
       endif
!      read scm data and initialize model
!      note:  usavg,vsavg and wsavg filled from here
       SCMDATEB = 0.0
       SCMDATEE = 0.0
       ISCMMONTHI = monthi
       call init_scmdata(SCMDATEB,SCMDATEE,ISCMMONTHI)
c      do L=1,LM
c         write(iu_scm_prt,'(a6,i5,2(f10.4))') 'l u v ',l,
c    &         u(i_targ,j_targ,l),v(i_targ,j_targ,l)
c      enddo
#endif

        CALL CALC_AMPK(LM)

        DO J=J_0,J_1
        DO I=I_0,I_1
C**** SET SURFACE MOMENTUM TRANSFER TAU0
          TAUAVG(I,J)=CDM*WSAVG(I,J)**2
C**** SET LAYER THROUGH WHICH DRY CONVECTION MIXES TO 1
          DCLEV(I,J)=1.
C**** SET SURFACE SPECIFIC HUMIDITY FROM FIRST LAYER HUMIDITY
          QSAVG(I,J)=Q(I,J,1)
          QGAVG(I,J)=Q(I,J,1)
          TGVAVG(I,J)=T(I,J,1)
C**** SET RADIATION EQUILIBRIUM TEMPERATURES FROM LAYER LM TEMPERATURE
          DO K=1,LM_REQ
            RQT(K,I,J)=T(I,J,LM)
          END DO
C**** REPLACE TEMPERATURE BY POTENTIAL TEMPERATURE
          DO L=1,LM
            T(I,J,L)=T(I,J,L)/PK(L,I,J)
            TTOLD(L,I,J)=T(I,J,L)
            QTOLD(L,I,J)=Q(I,J,L)
          END DO
C**** initialize egcm to be used in ATURB.f
          DO L=1,LM
            egcm(l,i,j)=egcm_init_max/(float(l)**2)
            w2gcm(l,i,j)=egcm(l,i,j)*2.*by3
          END DO
        END DO
        END DO
C**** Initialize surface friction velocity
        DO J=J_0,J_1
        DO I=I_0,I_1
          USTAR_pbl(:,I,J)=WSAVG(I,J)*SQRT(CDM)
        END DO
        END DO
C**** INITIALIZE VERTICAL SLOPES OF T,Q
        call tq_zmom_init(t,q,PMID,PEDN)

      END IF

C****
C**** I.C from possibly older/incomplete MODEL OUTPUT, ISTART=3-8
C****
      SELECT CASE (ISTART)
      CASE (3)               ! just general hints - not to be used as is
C**** Read what's there and substitute rest as needed (as above)
C**** To be implemented as needed. Sometimes it is safer to
C**** combine the ground layers into 2 layers (top 10cm and rest) and
C**** set   redoGH  to .true.  (after major changes in the GH code or
C**** after changing to a new horizontal grid)
C     redoGH=.TRUE.
C**** Set flag to initialise pbl/snow variables if obsolete or missing
C     iniPBL=.TRUE.  ; iniSNOW = .TRUE.
C**** if iniPBL is to be used, you must set A-grid winds as well
C     call recalc_agrid_uv
        go to 890            !  not implemented; stop model
C****
C**** I.C FROM FULL MODEL RESTART FILE (but re-initialise ocean)
C****
      CASE (4)
        call io_rsf("AIC",IhrX,irsficno,ioerr)
        if (ioerr.eq.1) goto 800
        iniOCEAN = .TRUE. ! read in ocean ic
C****
C**** I.C FROM FULL MODEL RESTART FILE (but no tracers)
C****
      CASE (5)             ! this model's rsf file, no tracers
        call io_rsf("AIC",IhrX,irsficnt,ioerr)
        if (ioerr.eq.1) goto 800
C****
C**** I.C FROM RESTART FILE that may not match land-ocean mask  ISTART=6
C****
      CASE (6)             ! converted model II' (B399) format (no snow)
        call io_rsf("AIC",IhrX,irsficno,ioerr)
        if (ioerr.eq.1) goto 800
        init_topog_related = 1
        redoGH=.TRUE.
C****
C**** I.C FROM RESTART FILE WITH almost COMPLETE DATA    ISTART=7
C****
      CASE (7)             ! converted model II' (B399) format (no snow)
        call io_rsf("AIC",IhrX,irsfic,ioerr)
        if (ioerr.eq.1) goto 800
        iniSNOW = .TRUE.      ! extract snow data from first soil layer
        iniOCEAN = .TRUE. ! read in ocean ic
        inilake = .TRUE. ! use flake from topog file
C****
C****   Data from current type of RESTART FILE           ISTART=8
C****
      CASE (8)  ! no need to read SRHR,TRHR,FSF,TSFREZ,diag.arrays
        call io_rsf("AIC",IhrX,irsfic,ioerr)
        if (ioerr.eq.1) goto 800
      END SELECT

C**** Check consistency of starting time
      IF (ISTART.ge.3) THEN
        IF( (MOD(IHRI-IHRX,8760).ne.0) ) THEN
         WRITE (6,*) ' Difference in hours between ',
     *       'Starting date and Data date:',MOD(IHRI-IHRX,8760)
         WRITE (6,*) 'Please change HOURI,DATEI,MONTHI'
         call stop_model('INPUT: start date inconsistent with data',255)
        ENDIF
      END IF

C**** Set flags to initialise some variables related to topography
      call sync_param( "init_topog_related", init_topog_related )
      IF (init_topog_related == 1) then
        iniOcean=.true.
        iniLAKE=.TRUE.
        iniPBL=.TRUE.
        iniSNOW = .TRUE.        ! extract snow data from first soil layer
        do_IC_fixups = 1        ! new default, not necessarily final
      endif

      CALL CALC_AMPK(LM)
C****
!**** IRANDI seed for random perturbation of initial conditions (if/=0):
C****        tropospheric temperatures changed by at most 1 degree C
      IF (IRANDI.NE.0) THEN
        CALL RINIT (IRANDI)
        DO L=1,LS1-1
        call burn_random(nij_before_j0(J_0))
        DO J=J_0,J_1
        call burn_random((I_0-1))
        DO I=I_0,I_1
           TIJL=T(I,J,L)*PK(L,I,J)-1.+2*RANDU(X)
           T(I,J,L)=TIJL/PK(L,I,J)
        END DO
        call burn_random(nij_after_i1(I_1))
        END DO
        call burn_random(nij_after_j1(J_1))
        END DO
        IF (AM_I_ROOT())
     *       WRITE(6,*) 'Initial conditions were perturbed !!',IRANDI
      END IF

      IF (AM_I_ROOT())
     *     WRITE(6,'(A,i3,1x,a4,i5,a3,i3,3x,a,i2/" ",a)')
     *  '0Model started on',datei,aMONTH(monthi),yeari,' Hr',houri,
     *  'ISTART =',ISTART,XLABEL(1:80)    ! report input file label
      XLABEL = RLABEL                     ! switch to rundeck label

      GO TO 600
C***********************************************************************
C****                                                               ****
C****                  RESTARTS: ISTART > 8                         ****
C****                                                               ****
C****   Current settings: 9 - from own model M-file                 ****
C****                    10 - from later of fort.1 or fort.2        ****
C****                    11 - from fort.1                           ****
C****                    12 - from fort.2                           ****
C****               13 & up - from earlier of fort.1 or fort.2      ****
C****                                                               ****
C***********************************************************************
  400 SELECT CASE (ISTART)
C****
C****   DATA FROM end-of-month RESTART FILE     ISTART=9
C****        mainly used for REPEATS and delayed EXTENSIONS
      CASE (1:9)                      !  diag.arrays are not read in
        if(istart.eq.9) call io_rsf("AIC",Itime,irerun,ioerr)
        if(istart.le.8) then         !  initial start of rad.forcing run
          call openunit("AIC",iu_AIC,.true.,.true.)
          call io_label(iu_AIC,Itime,ItimeX,irerun,ioerr)
          if (Kradia.gt.0) call io_rad (iu_AIC,irsfic,ioerr)
          call closeunit(iu_AIC)
        end if
        if (ioerr.eq.1) goto 800
        WRITE (6,'(A,I2,A,I11,A,A/)') '0Model restarted; ISTART=',
     *    ISTART,', TIME=',Itime,' ',XLABEL(1:80) ! sho input file label
        XLABEL = RLABEL                        ! switch to rundeck label

        CALL CALC_AMPK(LM)
C****
!**** IRANDI seed for random perturbation of current state (if/=0)
C****        tropospheric temperatures are changed by at most 1 degree C
        IF (IRANDI.ne.0 .and. Kradia.le.0) THEN
          CALL RINIT (IRANDI)
          DO L=1,LS1-1
          DO J=J_0,J_1
          DO I=I_0,I_1
             TIJL=T(I,J,L)*PK(L,I,J)-1.+2*RANDU(X)
             T(I,J,L)=TIJL/PK(L,I,J)
          END DO
          END DO
          END DO
          IF (AM_I_ROOT())
     *         WRITE(6,*) 'Current temperatures were perturbed !!',IRANDI
        END IF
        TIMING = 0
        GO TO 500
C****
C**** RESTART ON DATA SETS 1 OR 2, ISTART=10 or more
C****
C**** CHOOSE DATA SET TO RESTART ON
      CASE (10,13:)
         call find_later_rsf(kdisk)
         IF (ISTART.GE.13)     KDISK=3-KDISK
      CASE (11,12)
                               KDISK=ISTART-10
      END SELECT
  430 continue
#ifdef USE_FVCORE
        if(AM_I_ROOT()) then
          write(suffix,'(i1)') kdisk
          call system('cp  fv.'// suffix // ' dyncore_internal_restart')
          call system('cp dfv.'// suffix // ' tendencies_checkpoint')
        endif
#endif
      CALL HERE(__FILE__//'::io_rsf',__LINE__ + 10000*KDISK)
      call io_rsf(rsf_file_name(KDISK),Itime,ioread,ioerr)
      kdisk_sv = kdisk
      if (ioerr.eq.1) then
         if (istart.gt.10) go to 850  ! no 2nd chance if istart/=10
         KDISK=3-KDISK                ! try the earlier restart file
         WRITE (6,'(A,I1,A,I1)')
     *     ' Read Error on fort.',3-kdisk,' trying fort.',kdisk
         ISTART=110
         go to 430
      end if
      if (AM_I_ROOT())
     * WRITE (6,'(A,I2,A,I11,A,A/)') '0RESTART DISK READ, UNIT',
     *   KDISK,', Time=',Itime,' ',XLABEL(1:80)

C**** Switch KDISK if the other file is (or may be) bad (istart>10)
C****     so both files will be fine after the next write execution
      IF (istart.gt.10) KDISK=3-KDISK
C**** Keep KDISK after reading from the later restart file, so that
C****     the same file is overwritten first; in case of trouble,
C****     the earlier restart file will still be available

  500 CONTINUE
C**** Get parameters we just read from rsf file. Only those
C**** parameters which we need in "INPUT" should be extracted here.
      if(is_set_param("DTsrc"))  call get_param( "DTsrc", DTsrc )
      if(is_set_param("DT"))     call get_param( "DT", DT )
      if(is_set_param("NMONAV")) call get_param( "NMONAV", NMONAV )
      if(is_set_param("Kradia")) call get_param( "Kradia", Kradia )

C***********************************************************************
C****                                                              *****
C****       INITIAL- AND RESTARTS: Final Initialization steps      *****
C****                                                              *****
C***********************************************************************
  600 CONTINUE

C**** initialize Lrunid (length of the identifying part of XLABEL)
C****
      lid1 = INDEX(XLABEL,'(') -1
      if (lid1.lt.1) lid1=MAXLEN_RUNID+1
      lid2 = INDEX(XLABEL,' ') -1
      if (lid2.lt.1) lid2=MAXLEN_RUNID+1
      LRUNID = min(lid1,lid2)
      IF (LRUNID.gt.MAXLEN_RUNID) call stop_model
     *     ('INPUT: Rundeck name too long. Shorten to 32 char or less'
     *     ,255)

C**** Update ItimeE only if YearE or IhourE is specified in the rundeck
C****
      if(timee.lt.0) timee=houre*nday/HR_IN_DAY
      IF(yearE.ge.0) ItimeE = (( (yearE-iyear1)*JDperY +
     *  JDendofM(monthE-1)+dateE-1 )*HR_IN_DAY )*NDAY/HR_IN_DAY + TIMEE
C**** Alternate (old) way of specifying end time
      if(IHOURE.gt.0) ItimeE=IHOURE*NDAY/HR_IN_DAY

C**** Check consistency of DTsrc (with NDAY) and dt (with NIdyn)
      if (is_set_param("DTsrc") .and. nint(sday/DTsrc).ne.NDAY) then
        if (AM_I_ROOT())
     *        write(6,*) 'DTsrc=',DTsrc,' has to stay at/be set to',SDAY/NDAY
        call stop_model('INPUT: DTsrc inappropriately set',255)
      end if
      DTsrc = SDAY/NDAY
      call set_param( "DTsrc", DTsrc, 'o' )   ! copy DTsrc into DB

      NIdyn=nint(dtsrc/dt)
#ifndef USE_FVCORE
C**** NIdyn=dtsrc/dt(dyn) has to be a multiple of 2
C****
      if(istart>0) then
        NIdyn = 2*nint(.5*dtsrc/dt)
        if (is_set_param("DT") .and. nint(DTsrc/dt).ne.NIdyn) then
          if (AM_I_ROOT())
     *        write(6,*) 'DT=',DT,' has to be changed to',DTsrc/NIdyn
          call stop_model('INPUT: DT inappropriately set',255)
        end if
      end if
#else
C**** need a clock to satisfy ESMF interfaces
      call getdte(itimei,nday,iyear1,YEARI,MONTHI,jday,DATEI,HOURI,amon)
      MINTI = nint(mod( mod(Itimei*hrday/Nday,hrday) * 60d0, 60d0))
      call getdte(itimee,nday,iyear1,YEARE,MONTHE,jday,DATEE,HOURE,amon)
      MINTE = nint(mod( mod(Itimee*hrday/Nday,hrday) * 60d0, 60d0))
      clock = init_app_clock( (/ YEARI, MONTHI, DATEI, HOURI, MINTI,0/),
     &                        (/ YEARE, MONTHE, DATEE, HOURE, MINTE,0/),
     &             interval = int(dt) )
C**** copy FV restart files to expected names
      if(AM_I_ROOT() .and. istart.gt.2 .and. istart.lt.10) then
        call system('cp AICfv  dyncore_internal_restart')
        call system('cp AICdfv tendencies_checkpoint')
      endif
#endif
      DT = DTsrc/NIdyn
      call set_param( "DT", DT, 'o' )         ! copy DT into DB

C**** NMONAV has to be 1(default),2,3,4,6,12, i.e. a factor of 12
      if (NMONAV.lt.1 .or. MOD(12,NMONAV).ne.0) then
        write (6,*) 'NMONAV has to be 1,2,3,4,6 or 12, not',NMONAV
        call stop_model('INPUT: nmonav inappropriately set',255)
      end if
      if (AM_I_ROOT())
     *     write (6,*) 'Diag. acc. period:',NMONAV,' month(s)'

C**** Updating Parameters: If any of them changed beyond this line
C**** use set_param(.., .., 'o') to update them in the database (DB)

C**** Get the rest of parameters from DB or put defaults to DB
      call init_Model

C**** Set julian date information
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
      call getdte(Itime0,Nday,iyear1,Jyear0,Jmon0,J,Jdate0,Jhour0,amon0)

C****
C**** READ IN TIME-INDEPENDENT ARRAYS
C****
      if (Kradia.le.0) then   !  full model
        CALL CALC_AMPK(LM)
      end if  ! full model: Kradia le 0

C**** READ IN LANDMASKS AND TOPOGRAPHIC DATA
C**** Note that FLAKE0 is read in only to provide initial values
C**** Actual array is set from restart file.
      call openunit("TOPO",iu_TOPO,.true.,.true.)

      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),FOCEAN,1) ! Ocean fraction
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),FLAKE0,1) ! Orig. Lake fraction
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),FEARTH0,1) ! Earth frac. (no LI)
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),FLICE ,1) ! Land ice fraction
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),ZATMO ,1) ! Topography
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),HLAKE ,2) ! Lake Depths
      ZATMO(:,J_0:J_1) = ZATMO(:,J_0:J_1)*GRAV                  ! Geopotential

C**** Deal with single -> double precision problems and potential
C**** ocean/lake inconsistency. Adjust FLAKE0 and FLICE if necessary.
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        IF (FOCEAN(I,J).gt.0) THEN
          FLAND(I,J)=1.-FOCEAN(I,J) ! Land fraction if focean>0
          IF (FLAKE0(I,J).gt.0) THEN
            WRITE(6,*) "Ocean and lake cannot co-exist in same grid box"
     *       ,i,j,FOCEAN(I,J),FLAKE0(I,J)
            FLAKE0(I,J)=0
          END IF
        ELSEIF (FLAKE0(I,J).gt.0) THEN
          FLAND(I,J)=1.-FLAKE0(I,J)  ! for initialization only
        ELSE
          FLAND(I,J)=1.              ! for initialization only
        END IF
C**** Ensure that no round off error effects land with ice and earth
        IF (FLICE(I,J)-FLAND(I,J).gt.-1d-4 .and. FLICE(I,J).gt.0) THEN
          FLICE(I,J)=FLAND(I,J)
        END IF
      END DO
      END DO
      call closeunit(iu_TOPO)

      CALL HALO_UPDATE(GRID, FOCEAN)
      CALL HALO_UPDATE(GRID, FLAKE0)
      CALL HALO_UPDATE(GRID, FEARTH0)
      CALL HALO_UPDATE(GRID, ZATMO)

C**** Check polar uniformity
      if(have_south_pole) then
        do i=2,im
          if (zatmo(i,1).ne.zatmo(1,1)) then
            print*,"Polar topography not uniform, corrected",i,1
     *           ,zatmo(i,1),zatmo(1,1)
            zatmo(i,1)=zatmo(1,1)
          end if
        end do
      end if
      if (have_north_pole) then
        do i=2,im
          if (zatmo(i,jm).ne.zatmo(1,jm)) then
            print*,"Polar topography not uniform, corrected",i,jm
     *           ,zatmo(i,jm),zatmo(1,jm)
            zatmo(i,jm)=zatmo(1,jm)
          end if
        end do
      end if

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** Initialise tracer parameters and diagnostics
C**** MUST be before other init routines
      call init_tracer
#endif
#ifdef TRACERS_ON
      if(istart.le.2) then
        call COMPUTE_GZ(p,t,tmom(mz,:,:,:),daily_z)
        daily_z = daily_z/grav
      endif
#endif

      call sync_param ("do_IC_fixups", do_IC_fixups)

!!! hack:
      istart_fixup = istart
      if(istart==8 .and. do_IC_fixups==1 ) istart_fixup = 9

C**** Initialise some modules before finalising Land/LI mask
C**** Initialize ice
      CALL init_ice(iniOCEAN,istart_fixup)

C**** Initialize lake variables (including river directions)
      CALL init_LAKES(inilake,istart_fixup)

C**** Initialize ice dynamics code (if required)
      CALL init_icedyn(iniOCEAN)

C**** Initialize ocean variables
C****  KOCEAN = 1 => ocean heat transports/max. mixed layer depths
C****  KOCEAN = 0 => RSI/MSI factor
      CALL init_OCEAN(iniOCEAN,istart_fixup)

C**** Initialize land ice (must come after oceans)
      CALL init_LI(istart_fixup)
C**** Make sure that constraints are satisfied by defining FLAND/FEARTH
C**** as residual terms.
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
!!      FLAND(I,J)=1.-FOCEAN(I,J)  !! already set if FOCEAN>0
        IF (FOCEAN(I,J).le.0) THEN
          FLAND(I,J)=1
          IF (FLAKE(I,J).gt.0) FLAND(I,J)=1.-FLAKE(I,J)
        END IF
        FEARTH(I,J)=FLAND(I,J)-FLICE(I,J) ! Earth fraction
      END DO
      END DO

      If (HAVE_SOUTH_POLE) Then
         FLAND(2:IM,1)=FLAND(1,1)
         FEARTH(2:IM,1)=FEARTH(1,1)
         FLICE(2:IM,1)=FLICE(1,1)
      End If
      If (HAVE_NORTH_POLE) Then
         FLAND(2:IM,JM)=FLAND(1,JM)
         FEARTH(2:IM,JM)=FEARTH(1,JM)
         FLICE(2:IM,JM)=FLICE(1,JM)
      End If

C****
C**** INITIALIZE GROUND HYDROLOGY ARRAYS (INCL. VEGETATION)
C**** Recompute Ground hydrology data if redoGH (new soils data)
C****
!!! hack: make sure that ISTART_kradia==0 if Kradia>0
!!! do we need it ? I.A.
      ISTART_kradia = ISTART
      if ( Kradia.gt.0 ) ISTART_kradia = 0
      CALL init_LSM(DTsrc/NIsurf,redoGH,iniSNOW,inilake,ISTART_kradia)
cddd#ifdef USE_ENT
cddd      CALL init_module_ent(iniENT, Jday, Jyear, FOCEAN) !!! FEARTH)
cddd#endif
#ifdef IRRIGATION_ON
      call init_irrigate()
#endif
      if (Kradia.gt.0) then   !  radiative forcing run
        CALL init_RAD(istart)
        if(istart.lt.0) CALL init_DIAG(0,num_acc_files) !post-processing
        if (AM_I_ROOT()) Then
           WRITE (6,INPUTZ)
           call print_param( 6 )
           WRITE (6,'(A14,4I4)') "IM,JM,LM,LS1=",IM,JM,LM,LS1
           WRITE (6,*) "PLbot=",PLbot
        end if
        if(istart.lt.0)
     &       CALL stop_model ('Terminated normally, istart<0',13)
        return
      end if                  !  Kradia>0; radiative forcing run

#if !defined(ADIABATIC) || defined(CUBED_SPHERE)

C**** Initialize pbl (and read in file containing roughness length data)
#ifndef CUBED_SPHERE /* until a better solution is found */
      if (iniPBL) call recalc_agrid_uv   ! PBL needs A-grid winds
#endif
      if(istart.gt.0) CALL init_pbl(iniPBL)
C****
C**** Initialize the use of gravity wave drag diagnostics
C****
      CALL init_GWDRAG
C
#ifdef CALCULATE_FLAMMABILITY
      CALL init_flammability
#endif
C**** Initialize nudging
#ifdef NUDGE_ON
      if (istart.gt.0) CALL NUDGE_INIT
#endif
#ifdef TRACERS_AMP
      CALL SETUP_CONFIG
      CALL SETUP_SPECIES_MAPS
      CALL SETUP_DP0
      CALL SETUP_AERO_MASS_MAP
      CALL SETUP_COAG_TENSORS
      CALL SETUP_DP0
      CALL SETUP_KIJ
      CALL SETUP_EMIS
      CALL SETUP_KCI
      CALL SETUP_NPFMASS
      CALL SETUP_DIAM
      CALL SETUP_RAD
#endif

#ifndef SCM
#ifndef CUBED_SPHERE
      if (USE_UNR_DRAG==1) CALL init_UNRDRAG
#endif
#endif

C****
      if(istart.gt.0) CALL RINIT (IRAND)
c Note on FFT initialization: IMLON is defined by the diag_zonal module,
c not by the resolution module.  IMLON==IM for a latlon grid.
      CALL FFT0 (IMLON)  ! CALL FFT0(IM)
      CALL init_CLD
      CALL init_DIAG(istart,num_acc_files) ! initialize for accumulation
      CALL UPDTYPE
      if(istart.gt.0) CALL init_QUS(grid,im,jm,lm)
#ifndef CUBED_SPHERE
      if(istart.gt.0) CALL init_ATMDYN
#endif
      CALL init_RAD(istart)
      if (AM_I_ROOT()) then
         WRITE (6,INPUTZ)
         call print_param( 6 )
         WRITE (6,'(A7,12I6)') "IDACC=",(IDACC(I),I=1,12)
         WRITE (6,'(A14,4I4)') "IM,JM,LM,LS1=",IM,JM,LM,LS1
         WRITE (6,*) "PLbot=",PLbot
      end if
#ifdef CACHED_SUBDD
      ! Initialize subdaily diagnostics
      call parse_subdd
      call reset_cached_subdd
      if(istart.ge.10) then
        call read_subdd_rsf(trim(rsf_file_name(kdisk_sv))//'.nc')
      endif
#endif

#endif

C****
      RETURN
C****
C**** TERMINATE BECAUSE OF IMPROPER PICK-UP
C****
  800 WRITE (6,'(A,I4/" ",A)')
     *  '0ERROR ENCOUNTERED READING AIC ISTART=', ISTART,XLABEL(1:80)
      call stop_model('INPUT: READ ERROR FOR AIC',255)
  850 WRITE (6,'(A)')
     *  '0ERRORS ON BOTH RESTART DATA SETS. TERMINATE THIS JOB'
      call stop_model('INPUT: ERRORS ON BOTH RESTART FILES',255)
  890 WRITE (6,'(A,I5)') '0INCORRECT VALUE OF ISTART',ISTART
      call stop_model('INPUT: ISTART-SPECIFICATION INVALID',255)
  900 write (6,*) 'Error in NAMELIST parameters'
      call stop_model('Error in NAMELIST parameters',255)
  910 write (6,*) 'Error readin I-file'
      call stop_model('Error reading I-file',255)
      END SUBROUTINE INPUT

      SUBROUTINE DAILY(end_of_day)
!@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
!@auth Original Development Team
!@ver  1.0
!@calls constant:orbit, calc_ampk, getdte
      USE MODEL_COM, only : im,jm,lm,ls1,ptop,psf,p,q
     *     ,itime,itimei,iyear1,nday,jdpery,jdendofm
     *     ,jyear,jmon,jday,jdate,jhour,aMON,aMONTH,ftype,ntype
      USE GEOM, only : areag,axyp,imaxj,lat2d
      USE DYNAMICS, only : byAM
      USE RADPAR, only : ghgam,ghgyr2,ghgyr1
      USE RAD_COM, only : RSDIST,COSD,SIND,COSZ_day,SUNSET,
     *     dh2o,H2ObyCH4,ghg_yr,
     *     omegt,obliq,eccn,omegt_def,obliq_def,eccn_def,
     *     variable_orb_par,orb_par_year_bp
#ifdef TRACERS_WATER
      USE TRACER_COM, only: trm,tr_wd_type,nwater,tr_H2ObyCH4,itime_tr0
     *     ,ntm
#endif
      USE DIAG_COM, only : aj=>aj_loc,j_h2och4
      USE DOMAIN_DECOMP_ATM, only : grid, GET, GLOBALSUM, AM_I_ROOT,
     &         WRITE_PARALLEL
c      USE ATMDYN, only : CALC_AMPK
      use RAD_COSZ0, only : daily_cosz
      IMPLICIT NONE
      REAL*8 :: DELTAP,PBAR,SMASS,SUNLON,SUNLAT,LAM,xCH4,xdH2O,EDPY,
     $          VEDAY,PYEAR
      REAL*8 :: CMASS(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                grid%J_STRT_HALO:grid%J_STOP_HALO)
      INTEGER i,j,l,iy,it
      character(len=300) :: out_line
      LOGICAL, INTENT(IN) :: end_of_day
#ifdef TRACERS_WATER
      INTEGER n
#endif
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Tasks to be done at end of day and at each start or restart
C****
C**** CALCULATE THE DAILY CALENDAR
C****
      call getdte(Itime,Nday,iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)

C**** CALCULATE SOLAR ANGLES AND ORBIT POSITION
C**** This is for noon (GMT) for new day.

C**** The orbital calculation will need to vary depending on the kind
C**** of calendar adopted (i.e. a generic 365 day year, or a transient
C**** calendar including leap years etc.).  For transient calendars the
C**** JDAY passed to orbit needs to be adjusted to represent the number
C**** of days from Jan 1 2000AD.
c      EDPY=365.2425d0, VEDAY=79.3125d0  ! YR 2000AD
c      JDAY => JDAY + 365 * (JYEAR-2000) + appropriate number of leaps
C**** Default calculation (no leap, VE=Mar 21 hr 0)
c      EDPY=365d0 ; VEDAY=79d0           ! Generic year
C**** PMIP calculation (no leap, VE=Mar 21 hr 12)
      EDPY=365d0 ; VEDAY=79.5d0           ! Generic year
C**** Update orbital parameters at start of year
      if (variable_orb_par == 1.and.JDAY == 1) then
        pyear = JYEAR - orb_par_year_bp ! bp=before present model year
        call orbpar(pyear, eccn, obliq, omegt)
        if (am_I_root()) then
          write(6,*) 'Set orbital parameters for year ',pyear,' (CE)'
          if (orb_par_year_bp.ne.0) write(6,*) 'offset by',
     *      orb_par_year_bp,' years from model year'
          write(6,*) "   Eccentricity: ",eccn
          write(6,*) "   Obliquity (degs): ",obliq
          write(6,*) "   Precession (degs from ve): ",omegt
        end if
      end if

      CALL ORBIT (OBLIQ,ECCN,OMEGT,VEDAY,EDPY,REAL(JDAY,KIND=8)-.5
     *     ,RSDIST,SIND,COSD,SUNLON,SUNLAT,LAM)
      call daily_cosz(sind,cosd,cosz_day,sunset)

      IF (.not.(end_of_day.or.itime.eq.itimei)) RETURN

C**** Tasks to be done at end of day and at initial starts only
C****
C**** THE GLOBAL MEAN PRESSURE IS KEPT CONSTANT AT PSF MILLIBARS
C****
C**** CALCULATE THE CURRENT GLOBAL MEAN PRESSURE
#ifndef SCM
      DO J=J_0,J_1
      DO I=I_0,I_1
        CMASS(I,J)=P(I,J)*AXYP(I,J)
      END DO
      END DO
      if(have_south_pole) cmass(2:im,1)=cmass(1,1)
      if(have_north_pole) cmass(2:im,jm)=cmass(1,jm)
      CALL GLOBALSUM(grid, CMASS, SMASS, ALL=.TRUE.)
      PBAR=SMASS/AREAG+PTOP
C**** CORRECT PRESSURE FIELD FOR ANY LOSS OF MASS BY TRUNCATION ERROR
C****   except if it was just done (restart from itime=itimei)
      DELTAP=PSF-PBAR
      if(itime.eq.itimei .and. abs(deltap).lt.1.d-10) return
      P=P+DELTAP

      CALL CALC_AMPK(LS1-1)

      if (AM_I_ROOT()) then
         IF (ABS(DELTAP).gt.1d-6)
     *      WRITE (6,'(A25,F10.6/)') '0PRESSURE ADDED IN GMP IS',DELTAP
      end if
#endif

      IF (.not.end_of_day) RETURN

C**** Tasks to be done at end of day only
      if (H2ObyCH4.gt.0) then
C****   Add obs. H2O generated by CH4(*H2ObyCH4) using a 2 year lag
        iy = jyear - 2 - ghgyr1 + 1
        if (ghg_yr.gt.0) iy = ghg_yr - 2 - ghgyr1 + 1
        if (iy.lt.1) iy=1
        if (iy.gt.ghgyr2-ghgyr1+1) iy=ghgyr2-ghgyr1+1
        xCH4=ghgam(3,iy)*H2ObyCH4
c        If (AM_I_ROOT())
c     &    write(6,*) 'add in stratosphere: H2O gen. by CH4(ppm)=',xCH4

        do l=1,lm
        do j=J_0,J_1
        do i=I_0,imaxj(j)
#ifdef CUBED_SPHERE
          call lat_interp_qma(lat2d(i,j),l,jmon,xdH2O)
#else
          xdH2O = dH2O(j,l,jmon)
#endif
          q(i,j,l)=q(i,j,l)+xCH4*xdH2O*byAM(l,i,j)
#ifdef TRACERS_WATER
C**** Add water to relevant tracers as well
          do n=1,ntm
            if (itime_tr0(n).le.itime) then
              select case (tr_wd_type(n))
              case (nWater)    ! water: add CH4-sourced water to tracers
                trm(i,j,l,n) = trm(i,j,l,n) +
     +                tr_H2ObyCH4(n)*xCH4*xdH2O*axyp(i,j)
              end select
            end if
          end do
#endif
          do it=1,ntype
            call inc_aj(i,j,it,j_h2och4,xCH4*xdH2O*ftype(it,i,j))
          end do
        end do
        end do
        If (HAVE_NORTH_POLE) q(2:im,jm,l)=q(1,jm,l)
        If (HAVE_SOUTH_POLE) q(2:im, 1,l)=q(1, 1,l)
#ifdef TRACERS_WATER
        do n=1,ntm
          If (HAVE_SOUTH_POLE) trm(2:im, 1,l,n)=trm(1, 1,l,n)
          If (HAVE_NORTH_POLE) trm(2:im,jm,l,n)=trm(1,jm,l,n)
        end do
#endif
        end do
      end if

      RETURN
      END SUBROUTINE DAILY

      SUBROUTINE CHECKT (SUBR)
!@sum  CHECKT Checks arrays for NaN/INF and reasonablness
!@auth Original Development Team
!@ver  1.0

C**** CHECKT IS TURNED ON BY SETTING QCHECK=.TRUE. IN NAMELIST
C**** REMEMBER TO SET QCHECK BACK TO .FALSE. AFTER THE ERRORS ARE
C**** CORRECTED.
      USE CONSTANT, only : tf
      USE MODEL_COM
      USE DYNAMICS, only : pk
      USE DOMAIN_DECOMP_ATM, only : grid, GET, AM_I_ROOT
      USE soil_drv, only : checke
      IMPLICIT NONE
      INTEGER I,J,L
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0H, J_1H, I_0,I_1, I_0H,I_1H, njpol
      INTEGER :: I_0STG,I_1STG,J_0STG,J_1STG
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     *     J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
c      I_0STG = grid%I_STRT_STGR
c      I_1STG = grid%I_STOP_STGR
      I_0STG = I_0
      I_1STG = I_1
      J_0STG = grid%J_STRT_STGR
      J_1STG = grid%J_STOP_STGR
      njpol = grid%J_STRT_SKP-grid%J_STRT

      IF (QCHECK) THEN
C**** Check all prog. arrays for Non-numbers
        CALL CHECK3B(U(I_0STG:I_1STG,J_0STG:J_1STG,:),
     &       I_0STG,I_1STG,J_0STG,J_1STG,0,LM,SUBR,'u     ')
        CALL CHECK3B(V(I_0STG:I_1STG,J_0STG:J_1STG,:),
     &       I_0STG,I_1STG,J_0STG,J_1STG,0,LM,SUBR,'v     ')
        CALL CHECK3B(T(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'t     ')
        CALL CHECK3B(Q(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'q     ')
        CALL CHECK3B(P(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &       SUBR,'p     ')
        CALL CHECK3B(WM(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'wm    ')
#ifdef BLK_2MOM
        CALL CHECK3B(WMICE(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'wmice    ')
#endif

        DO J=J_0,J_1
        DO I=I_0,I_1
          IF (Q(I,J,1).gt.1d-1)print*,SUBR," Q BIG ",i,j,Q(I,J,1:LS1)
          IF (T(I,J,1)*PK(1,I,J)-TF.gt.50.) print*,SUBR," T BIG ",i,j
     *         ,T(I,J,1:LS1)*PK(1:LS1,I,J)-TF
        END DO
        END DO
        DO L=1,LM
        DO J=J_0,J_1
        DO I=I_0,I_1
          IF (Q(I,J,L).lt.0.) then
            print*,"After ",SUBR," Q < 0 ",i,j,Q(I,J,L)
            call stop_model('Q<0 in CHECKT',255)
          END IF
          IF (WM(I,J,L).lt.0.) then
            print*,"After ",SUBR," WM < 0 ",i,j,WM(I,J,L)
            call stop_model('WM<0 in CHECKT',255)
          END IF
#ifdef BLK_2MOM
          IF (WMICE(I,J,L).lt.0.) then
            print*,"After ",SUBR," WMICE < 0 ",i,j,WMICE(I,J,L)
            call stop_model('WMICE<0 in CHECKT',255)
          END IF
#endif
        END DO
        END DO
        END DO
C**** Check PBL arrays
        CALL CHECKPBL(SUBR)
C**** Check Ocean arrays
        CALL CHECKO(SUBR)
C**** Check Ice arrays
        CALL CHECKI(SUBR)
C**** Check Lake arrays
        CALL CHECKL(SUBR)
C**** Check Earth arrays
        CALL CHECKE(SUBR)
C**** Check Land Ice arrays
        CALL CHECKLI(SUBR)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** check tracers
        CALL CHECKTR(SUBR)
#endif
      END IF

      RETURN
      END SUBROUTINE CHECKT


      subroutine read_options( coldRestart, ifile )
!@sum reads options from the command line (for now only one option)
!@auth I. Aleinov
!@ver 1.0
      implicit none
      character(*),intent(inout)  :: ifile
      logical, intent(inout) :: coldRestart
      character*80 arg,arg1

      do
        call nextarg( arg, 1 )
        if ( arg == "" ) exit          ! end of args
        select case (arg)
        case ("-i")
          call nextarg( arg1, 0 )
          ifile=arg1
        case ("-cold-restart")
          coldRestart = .true.
        ! new options can be included here
        case default
          print *,'Unknown option specified: ', arg
          print *,'Aborting...'
          call stop_model("Unknown option on a command line",255)
        end select
      enddo

      return
      end subroutine read_options


      subroutine print_and_check_PPopts
!@sum prints preprocessor options in english and checks some
!@+  interdependencies. (moved from subroutine INPUT).
!@+  Called by root thread only.
!@ver 1.0

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      write(6,*) 'This program includes tracer code'
#endif
#ifdef TRACERS_WATER
      write(6,*) '...and water tracer code'
#ifndef TRACERS_ON
      call stop_model(
     &' Water tracers need TRACERS_ON as well as TRACERS_WATER',255)
#endif
#endif
#ifdef TRACERS_OCEAN
      write(6,*) '...and ocean tracer code'
#endif
#ifdef TRACERS_SPECIAL_O18
      write(6,*) '...and water isotope code'
#ifndef TRACERS_WATER
      call stop_model('Water isotope tracers need TRACERS_WATER '//
     *'as well as TRACERS_SPECIAL_O18',255)
#endif
#endif
#ifdef TRACERS_SPECIAL_Lerner
      write(6,*) '...and Jean/David tracers and chemistry'
#endif
#ifdef TRACERS_GASEXCH_ocean
      write(6,*) '          '
      write(6,*) '...and Natassa Romanou air-sea GAS EXCHANGE'
#ifdef TRACERS_OceanBiology
      write(6,*) '          '
      write(6,*) '...and Natassa Romanou/Watson Gregg ocean biology '
#endif
#ifdef TRACERS_GASEXCH_ocean_CFC
      write(6,*) '****CFC flux across air/sea interface****'
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      write(6,*) '****CO2 flux across air/sea interface****'
#endif
#endif
#ifdef TRACERS_SPECIAL_Shindell
      write(6,*) '...and Drew Shindell tracers and chemistry'
#endif
#ifdef TRACERS_TERP
#ifdef TRACERS_SPECIAL_Shindell
      write(6,*) '...and Terpenes tracer'
#else
      call stop_model('Terpenes tracer needs tropo chemistry',255)
#endif
#endif  /* TRACERS_TERP */
#ifdef CALCULATE_FLAMMABILITY
      write(6,*) '...and calculating sfc veg flammability'
#endif
#ifdef CALCULATE_LIGHTNING
      write(6,*) '...and calculating lightning flash rate'
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
      write(6,*) '...and dynamic biomass burning srcs by flammability'
#endif
#ifdef TRACERS_AEROSOLS_Koch
      write(6,*) '...and Dorothy Koch aerosols'
#endif
#ifdef TRACERS_AMP
      write(6,*) '...and aerosol microphysics'
#endif
#ifdef TRACERS_AEROSOLS_SOA
#ifdef TRACERS_SPECIAL_Shindell
      write(6,*) '...and secondary organic aerosols'
#else
      call stop_model('SOA version needs tropo chemistry',255)
#endif
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef SOA_DIAGS
#ifdef TRACERS_AEROSOLS_SOA
      write(6,*) '...and additional SOA diagnostics'
#else
      call stop_model('SOA_DIAGS needs TRACERS_AEROSOLS_SOA',255)
#endif  /* TRACERS_AEROSOLS_SOA */
#endif  /* SOA_DIAGS */
#ifdef TRACERS_DRYDEP
      write(6,*) '...and tracer dry deposition'
#endif
#ifdef EDGAR_HYDE_SOURCES
      write(6,*) '...and EDGAR HYDE sources instead of GISS'
#endif
#ifdef SHINDELL_STRAT_CHEM
      write(6,*) '...and Drew Shindell stratospheric chemistry'
#endif
#ifdef SHINDELL_STRAT_EXTRA
      write(6,*) '...and Drew Shindell extra strat tracers'
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      write(6,*) '...and interactive CH4 wetlands emissions'
#endif
#ifdef NUDGE_ON
      write(6,*) '...and nudging of meteorology'
#endif
#ifdef HTAP_LIKE_DIAGS
      write(6,*) '...and HTAP set of diagnostics'
#endif
#ifdef ACCMIP_LIKE_DIAGS
      write(6,*) '...and ACCMIP set of diagnostics'
#ifndef SHINDELL_STRAT_EXTRA
      call stop_model
     & ('SHINDELL_STRAT_EXTRA should be on for ACCMIP_LIKE_DIAGS',255)
#endif
#ifndef HTAP_LIKE_DIAGS
      call stop_model
     & ('HTAP_LIKE_DIAGS should be on for ACCMIP_LIKE_DIAGS',255)
#endif
#endif /* ACCMIP_LIKE_DIAGS */

      return
      end subroutine print_and_check_PPopts
