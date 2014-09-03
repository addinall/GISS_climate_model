#include "rundeck_opts.h"
      MODULE MODEL_COM
!@sum  MODEL_COM Main model variables, independent of resolution
!@auth Original Development Team
!@ver  1.0
      USE RESOLUTION, only : im,jm,lm,ls1,kep,istrat,
     *     psf,pmtop,ptop,psfmpt,pstrat,plbot
!      USE DOMAIN_DECOMP_ATM, only : grid
!AOO      USE ESMF_CUSTOM_MOD, ONLY: FIELD
c$$$#ifdef USE_FVCORE
c$$$      USE ESMF_MOD, only: esmf_clock
c$$$#endif
      IMPLICIT NONE
! just to make all compilers happy (should check later)
#if ! defined(COMPILER_G95)
      SAVE
#endif

!@param IMH half the number of latitudinal boxes
      INTEGER, PARAMETER :: IMH=IM/2
!@param IVNP,IVSP V at north/south pole is stored in U(IVNP,JM)/U(IVSP,1)
      INTEGER, PARAMETER :: IVNP = IM/4 , IVSP = 3*IM/4
!@param FIM,BYIM real values related to number of long. grid boxes
      REAL*8, PARAMETER :: FIM=IM, BYIM=1./FIM
!@param JEQ grid box zone around or immediately north of the equator
      INTEGER, PARAMETER :: JEQ=1+JM/2

      CHARACTER*132 XLABEL !@var XLABEL=runID+brief description of run
      INTEGER :: LRUNID    !@var Run name stored in XLABEL(1:LRUNID)
      INTEGER :: LHEAD=15  !@var length of crucial beg of module_headers

!@var LM_REQ Extra number of radiative equilibrium layers
      INTEGER, PARAMETER :: LM_REQ=3
!@var REQ_FAC/REQ_FAC_M factors for REQ layer pressures
      REAL*8, PARAMETER, DIMENSION(LM_REQ-1) ::
     *     REQ_FAC=(/ .5d0, .2d0 /)               ! edge
      REAL*8, PARAMETER, DIMENSION(LM_REQ) ::
     *     REQ_FAC_M=(/ .75d0, .35d0, .1d0 /),    ! mid-points
     *     REQ_FAC_D=(/ .5d0,  .3d0,  .2d0 /)     ! delta

!**** Vertical resolution dependent variables (set in INPUT)
!@var SIGE sigma levels at layer interfaces (1)
!!!!  Note:   sige(1)=1,  sige(ls1)=0,  sige(lm+1)=-pstrat/psfmpt
      REAL*8, DIMENSION(LM+1) :: SIGE
!@var SIG,DSIG,byDSIG mid point, depth, 1/depth of sigma levels (1)
      REAL*8, DIMENSION(LM) ::
     &     SIG,    ! = (sige(1:lm)+sige(2:lm+1))*0.5d0,
     &     DSIG,   ! =  sige(1:lm)-sige(2:lm+1),
     &     byDSIG  ! =  1./DSIG

!@var PL00, PMIDL00, PDSIGL00, AML00 press (mb), mid-pressure (mb),
!@+        mass (kg/m2) for mean profile
!@var PEDNL00 edge pressure for mean profile (mb)
      REAL*8, DIMENSION(LM+LM_REQ) :: PL00, PMIDL00, PDSIGL00, AML00,
     *     BYAML00
      REAL*8, DIMENSION(LM+LM_REQ+1) :: PEDNL00

!**** Model control parameters:
!@dbparam KOCEAN: if 0 => specified, if 1 => predicted ocean
!@dbparam MFILTR: if 1 => SLP, if 2 => T, if 3 => SLP&T is filtered
      integer :: KOCEAN = 1, MFILTR = 1
!@dbparam COUPLED_CHEM: if 0 => uncoupled, if 1 => coupled
      integer :: COUPLED_CHEM = 0
!@var ij_debug: if i > 0, print out some extra info on bad ij box
      integer, dimension(2) :: ij_debug = (/ 0 , 1 /)

!**** Controls on FLTRUV (momentum/velocity filter)
!@dbparam DT_XUfilter dU is multiplied by dt/DT_XUfilter in E-W
!@dbparam DT_XVfilter dV is multiplied by dt/DT_XVfilter in E-W
!@dbparam DT_YUfilter dU is multiplied by dt/DT_YUfilter in N-S
!@dbparam DT_YVfilter dV is multiplied by dt/DT_YVfilter in N-S
      REAL*8 :: DT_XUfilter=0. ! U-filter is NOT used in E-W direction
      REAL*8 :: DT_XVfilter=0. ! V-filter is NOT used in E-W direction
      REAL*8 :: DT_YUfilter=0. ! U-filter is NOT used in N-S direction
      REAL*8 :: DT_YVfilter=0. ! V-filter is NOT used in N-S direction
!     Target Coordinates for SCM
      INTEGER*4 :: I_TARG,J_TARG   !TWP I=125,J=39  set targets in parameter list
      INTEGER*4 :: NSTEPSCM=0      !Time step counter for SCM
!@var QUVfilter: True if any of DT_[XY][UV]filter are not=0
      LOGICAL :: QUVfilter
!@dbparam ang_uv =1 to conserve ang mom in UVfilter
      INTEGER :: ang_uv = 1 ! UV-filter conserves ang mom

!@dbparam X_SDRAG.  SDRAG ~X_SDRAG(1)+X_SDRAG(2)*wind_magnitude
      REAL*8, DIMENSION(2) :: X_SDRAG = (/2.5D-4,2.5D-5/)
!@dbparam C_SDRAG.  SDRAG=C_SDRAG (const.) above PTOP
      REAL*8 :: C_SDRAG = 2.5D-5
      REAL*8, DIMENSION(LS1:LM) :: CSDRAGL
!@dbparam P_CSDRAG pressure level above which const.drag is increased
      REAL*8 :: P_CSDRAG=0.
!@dbparam P(P)_SDRAG pressure level above which SDRAG is applied (mb)
      REAL*8 :: P_SDRAG=0., PP_SDRAG = 1.d0 ! (PP_... near poles)
!@var L(P)SDRAG level above which SDRAG is applied (near pole)
      INTEGER :: LSDRAG=LM, LPSDRAG=LM  ! non-polar, polar limit
!@var ANG_SDRAG if =1: ang.momentum lost by SDRAG is added in below PTOP
      INTEGER :: ANG_SDRAG=1  ! default: SDRAG does conserve ang.mom
!@dbparam Wc_JDRAG critical velocity for J.Hansen/Judith Perlwitz drag
      REAL*8 :: Wc_JDRAG=30.d0  !  if 0.: no JDRAG-feature in Sdrag
!@dbparam wmax imposed limit for stratospheric winds (m/s) in SDRAG
      real*8 :: wmax=200.d0
!@dbparam do_polefix if =1 : u,v tendencies are corrected near the pole
      INTEGER :: do_polefix=1     ! default is to enable corrections
!@dbparam VSDRAGL layer dependent tuning factor for stratospheric drag
!@+   (not =1 e.g. if used with explicit grav.wave drag scheme)
      real*8 :: VSDRAGL(LS1:LM) = 1d0
!@dbparam  USE_UNR_DRAG   if 1 =>SDRAG is turned off and GWD is applied.
!@+    if 0 => SDRAG is kept intact and alternative GWD is not employed.
      INTEGER :: USE_UNR_DRAG=0  ! default: SDRAG is kept intact.

!**** Diagnostic control parameters
!@dbparam KCOPY: if 1 => acc, if 2 => +rsf, if 3 => +od are saved
!@dbparam NMONAV number of months in a diagnostic accuml. period
!@dbparam Kvflxo if 1 => vert.fluxes into ocean are saved daily
!@dbparam Kradia if -1 save data for, if 1|2 do   inst|adj forcing run
!@dbparam NIPRNT number of instantaneous initial printouts
      integer :: KCOPY=2, NMONAV=1, Kvflxo=0, Kradia=0,iu_rad, NIPRNT=1

C**** (Simplified) Calendar Related Terms
!@param JDperY,JMperY    number of days,months per year
!@var   JDendOfM(0:12)   last Julian day in month
!@var   JDmidOfM(0:13)   middle Julian day in month
      integer, PARAMETER :: JDPERY = 365, JMPERY = 12
      integer :: JDendOfM(0:JMPERY) = (
     *     /0,31,59,90,120,151,181,212,243,273,304,334,365/)
      integer :: JDmidOfM(0:JMPERY+1) = (
     *     /-15,16,45,75,106,136,167,197,228,259,289,320,350,381/)

!@var AMON,AMONTH(0:12)  (3-4 letter) names for current,all months
!@var AMON0  (3-4 letter) name of first month of the current acc-period
      CHARACTER*4 :: AMON='none',AMON0='none', AMONTH(0:12) = (/'IC  ',
     *  'JAN ','FEB ','MAR ','APR ','MAY ','JUNE',
     *  'JULY','AUG ','SEP ','OCT ','NOV ','DEC '/)

!@var NDAY and IYEAR1 relate CALENDAR TIME and INTERNAL TIME Itime :
!@var NDAY number of Internal Time Units per day (1 ITU = DTsrc sec)
!@nlparam IYEAR1  year 1 of internal clock (Itime=0 to 365*NDAY)
      INTEGER :: NDAY,IYEAR1=-1   !@var relate internal to calendar time

!@var ITIME current time in ITUs (1 ITU = DTsrc sec, currently 1 hour)
!@var JDAY,JMON,JDATE,JYEAR,JHOUR current Julian day,month,day,year,hour
      INTEGER :: Itime,JDAY,JMON,JDATE,JYEAR,JHOUR
!@var ItimeI,ItimeE   time at start,end of run
!@var Itime0          time at start of current accumulation period
!@var JMON0,JDATE0,JYEAR0,JHOUR0 date-info about Itime0 (beg.of acc.per)
      INTEGER :: ItimeI,ItimeE,   Itime0,JMON0,JDATE0,JYEAR0,JHOUR0

!@var ESMF clock required for some interfaces
c$$$#ifdef USE_FVCORE
c$$$      Type (ESMF_CLOCK) :: clock
c$$$#endif
!@dbparam DTSRC source time step (s)   = 1 ITU
      REAL*8 :: DTsrc = 3600.
!@dbparam DT (atmospheric) dynamics time step (s)
      REAL*8 :: DT    =  450.         ! DT = DTdyn_atm

C**** Time step related multipliers:  N... NI...
C**** general rule:   DTxxx = Nxxx*DTsrc  and  DTxxx = DTsrc/NIxxx
C**** except that the time steps related to NDAa, NDA5k, NDAsf are
C**** slightly larger, to sample all points within the cycle

!@var     NIdyn:  DT atm_dyn  =  DTsrc/NIdyn     (NIdyn=DTsrc/DT)
!@dbparam NIsurf: DT_Surface  =  DTsrc/NIsurf
!@dbparam NRad:   DT_Rad      =  NRad*DTsrc
!@dbparam NFILTR: DT_filter   =  NFILTR*DTsrc
      INTEGER :: NIdyn, NIsurf = 2, NRad = 5 , NFILTR = 1

!@dbparam Ndisk:  DT_saversf    =  Ndisk *DTsrc fort.1/fort.2 saves
!@dbparam Nssw:   DT_checkSsw   =  Nssw  *DTsrc
      INTEGER :: NDisk = 24, Nssw = 1

!@dbparam NDAA:   DT_DiagA    =  NDAA*DTsrc + 2*DT(dyn)
!@dbparam NDA5k:  DT_Diag5k   =  NDA5k*DTsrc + 2*DT(dyn) SpAnal KE
!@dbparam NDA5d:  DT_Diag5d   =  NDA5d*DTsrc     Consrv  SpAnal dyn
!@dbparam NDA5s:  DT_Diag5s   =  NDA5s*DTsrc     Consrv  SpAnal src
!@dbparam NDASf:  DT_DiagSrfc =  NDASf*DTsrc + DTsrc/NIsurf
!@dbparam NDA4:   DT_Diag4    =  NDA4 *DTsrc   Energy history
      INTEGER :: NDAa=7, NDA5d=1, NDA5k=7, NDA5s=1, NDASf=1, NDA4=24

!**** Accounting variables
!@dbparam IRAND last seed used by rand.number generator
!@var KDISK next rsf (fort.)1 or 2 to be written to
!@var NSTEP number of dynamics steps since start of run
!@var MRCH  flags position in dynamics cycle (>0 fw, <0 bw step)
      INTEGER :: IRAND=123456789, KDISK=1, NSTEP,MRCH
!@param rsf_file_name names of restart files
      CHARACTER(6), PARAMETER :: rsf_file_name(2)=(/'fort.1','fort.2'/)
!@var MODRD,MODD5K,MODD5S: if MODxxx=0 do xxx, else skip xxx
      INTEGER :: MODRD, MODD5K, MODD5S
!@var MDYN,MCNDS,MRAD,MSURF,MDIAG,MELSE timing-indices
      INTEGER  MDYN,MCNDS,MRAD,MSURF,MDIAG,MELSE
!@param NSAMPL number of diagnostic sampling schemes
      INTEGER, PARAMETER :: NSAMPL = 12
!@var IDACC(NSAMPL) counters for diagn. accumulations
      INTEGER, DIMENSION(NSAMPL) :: IDACC


!**** IO read/write flags used by the io_xyz routines
!@param IOWRITE Flag used for writing normal restart files
!@param IOWRITE_SINGLE Flag used for saving diags in single precision
!@param IOWRITE_MON Flag used for saving restart part only (no diags)
!@param IOREAD Flag used for reading in (composite) restart files
!@param IOREADNT Flag used for reading in restart files (w/o tracers)
!@param IRSFIC Flag used for reading in restart part to start NEW run
!@param IRSFICNT Flag used for reading restart (w/o tracers) for NEW run
!@param IRSFICNO Flag used for reading restart (w/o ocean) for NEW run
!@param IRERUN Flag used for reading in restart part to extend OLD run
      INTEGER, PARAMETER :: ioread=1,ioread_single=2,
     *     irerun=3,irsfic=4,irsficnt=5,ioreadnt=6,irsficno=7,
     *     ioread_nodiag=8,
     *     iowrite=-1,iowrite_single=-2,iowrite_mon=-3

!**** Main model prognostic variables
!@var U,V east-west, and north-south velocities (m/s)
!@var T potential temperature (referenced to 1 mb) (K)
!@var Q specific humidity (kg water vapor/kg air)
!@var WM cloud liquid water amount (kg water/kg air)
#ifdef BLK_2MOM
!@var WMICE cloud ice amount (kg water/kg air)
#endif
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: U
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: V
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: T
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: Q
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: WM
#ifdef BLK_2MOM
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: WMICE
#endif

!**** Boundary condition arrays:
!@var ZATMO,HLAKE Topography arrays: elevation (m), lake sill depth (m)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: ZATMO
      REAL*8, ALLOCATABLE,  DIMENSION(:,:)   :: HLAKE
!@var Fxx fraction of gridbox of type xx (land,ocean,...)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLAND
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FOCEAN
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLICE
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLAKE0
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FEARTH0

!@var WFCS water field capacity of first ground layer (kg/m2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: WFCS

!@var P surface pressure (hecto-Pascals - PTOP)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: P

C**** Define surface types (mostly used for weighting diagnostics)
!@param NTYPE number of different surface types
      INTEGER, PARAMETER :: NTYPE=6   ! orig = 3
!@var FTYPE fractions of each surface type
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: FTYPE

!@param ITxx indices of various types (used only when it matters)
      INTEGER, PARAMETER :: ITOCEAN=1, ITOICE=2, ITEARTH=3,
     *                      ITLANDI=4, ITLAKE=5, ITLKICE=6

C**** Variables specific for stratosphere and/or strat diagnostics
!@var DO_GWDRAG when true, prints Gravity Wave diagnostics
      LOGICAL :: DO_GWDRAG = .false.
!@var iDO_GWDRAG number if AIJ Gravity wave diagnostics
      INTEGER :: iDO_GWDRAG = 0

!@dbparam UOdrag parameter that decides whether ocean.ice velocities
!@+   feed into drag calculation in surface (default = 0)
      INTEGER :: UOdrag = 0

!@nlparam QCHECK TRUE for running diagnostic checks
      LOGICAL :: QCHECK = .FALSE.

!@var stop_on TRUE stops the model (set with "kill -15 PID)
      LOGICAL :: stop_on = .FALSE.

      END MODULE MODEL_COM

      SUBROUTINE ALLOC_MODEL_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
      USE RESOLUTION, ONLY : IM,JM,LM
      USE MODEL_COM, ONLY : NTYPE
      USE MODEL_COM, ONLY : ZATMO,HLAKE,FLAND,FOCEAN,FLICE,FLAKE0,
     *                      FEARTH0,WFCS,P,U,V,T,Q,WM,FTYPE
#ifdef BLK_2MOM
     *  ,WMICE
#endif
!AOO      USE ESMF_CUSTOM_MOD, ONLY: modelE_grid
!AOO      USE ESMF_CUSTOM_MOD, ONLY: ESMF_CELL_SFACE
!AOO      USE ESMF_CUSTOM_MOD, ONLY: ESMF_CELL_CENTER

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE(ZATMO(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(HLAKE(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(FLAND(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(FOCEAN(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(FLICE(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(FLAKE0(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(FEARTH0(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(WFCS(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(P(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(U(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(V(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(T(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(Q(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(WM(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)
#ifdef BLK_2MOM
      ALLOCATE(WMICE(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)
#endif


      ALLOCATE(FTYPE(NTYPE,I_0H:I_1H,J_0H:J_1H), STAT = IER)

! initialize non-initialized arrays
      FTYPE(:,:,:) = 0.d0

      END SUBROUTINE ALLOC_MODEL_COM


      MODULE TIMINGS
!@sum  TIMINGS contains variables for keeping track of computing time
!@auth Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE
      SAVE
!@param NTIMEMAX maximum number of possible time accumulators
      INTEGER, PARAMETER :: NTIMEMAX=12
!@var NTIMEACC actual number of time accumulators
      INTEGER :: NTIMEACC = 0
!@var TIMING array that holds timing info
      REAL*8, DIMENSION(0:NTIMEMAX) :: TIMING
!@var TIMESTR array that holds timing info description
      CHARACTER*12, DIMENSION(NTIMEMAX) :: TIMESTR

      END MODULE TIMINGS

      SUBROUTINE SET_TIMER(STR,MINDEX)
!@sum  SET_TIMER sets an index of TIMING for a particular description
!@auth Gavin Schmidt
!@ver  1.0
      USE TIMINGS
      IMPLICIT NONE
!@var STR string that describes timing accumulator
      CHARACTER*12, INTENT(IN) :: STR
!@var MINDEX index for that accumulator
      INTEGER, INTENT(OUT) :: MINDEX
      INTEGER N

C**** Check whether index has been set
      DO N=1,NTIMEACC
        IF (STR.EQ.TIMESTR(N)) THEN
          MINDEX=N
          RETURN
        END IF
      END DO
C**** Otherwise increase number of indexes
      NTIMEACC = NTIMEACC + 1
      IF (NTIMEACC.gt.NTIMEMAX) call stop_model(
     &     "Too many timing indices: increase NTIMEMAX",255)
      MINDEX = NTIMEACC
      TIMESTR(MINDEX) = STR
C****
      RETURN
      END SUBROUTINE SET_TIMER

      SUBROUTINE TIMER (NOW,MSUM)
!@sum  TIMER keeps track of elapsed CPU time in hundredths of seconds
!@auth Gary Russell
!@ver  1.0
      USE TIMINGS
      USE GETTIME_MOD
      IMPLICIT NONE
      REAL*8, INTENT(OUT) :: NOW     !@var NOW current CPU time (seconds)
      INTEGER, INTENT(INOUT) :: MSUM !@var MSUM index for running total
      REAL*8 :: INC                  !@var INC time since last call
      REAL*8, SAVE :: LAST = 0       !@var LAST  last CPU time
      REAL*8 :: CMAX                 !@var CMAX max.count before 0-reset

      CALL GETTIME(NOW, CMAX)
      INC  = NOW - LAST
      if(inc<0) inc=inc+cmax ! offset system_clock reset
      TIMING(MSUM)  = TIMING(MSUM) + INC
      LAST = NOW
      RETURN
      END SUBROUTINE TIMER

      SUBROUTINE TIMEOUT (BEGIN,MIN,MOUT)
!@sum  TIMEOUT redistributes timing info between counters
!@auth Gary Russell
!@ver  1.0
      USE TIMINGS
      USE GETTIME_MOD
      IMPLICIT NONE
!@var MBEGIN CPU time start of section (.01 s)
      REAL*8, INTENT(IN) :: BEGIN
      INTEGER, INTENT(INOUT) :: MIN  !@var MIN index to be added to
      INTEGER, INTENT(INOUT) :: MOUT !@var MOUT index to be taken from
      REAL*8 :: INC                  !@var INC time since MBEGIN
      REAL*8 :: NOW                  !@var NOW current CPU time (s)
      REAL*8 :: CMAX                 !@var CMAX max.count before 0-reset

      CALL GETTIME(NOW, CMAX)
      INC  = NOW - BEGIN
      if(inc<0) inc=inc+cmax ! offset system_clock reset
      TIMING(MIN)  = TIMING(MIN)  + INC
      TIMING(MOUT) = TIMING(MOUT) - INC
      RETURN
      END SUBROUTINE TIMEOUT

      SUBROUTINE io_label(kunit,it,itm,iaction,ioerr)
!@sum  io_label reads and writes label/parameters to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT
      USE TIMINGS, only : ntimemax,ntimeacc,timestr,timing
      USE PARAM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var it input/ouput value of hour
!@var itm maximum hour returned (different from it if post-processing)
      INTEGER, INTENT(INOUT) :: it,itm
!@var LABEL2 content of record 2
      CHARACTER*80 :: LABEL2
!@var NTIM1,TSTR1,TIM1 timing related dummy arrays
      INTEGER NTIM1,ITIM1(NTIMEMAX)
      REAL*8 ::      TIM1(NTIMEMAX)
      CHARACTER*12  TSTR1(NTIMEMAX)
!@var ITmin,ITmax minimal/maximal time in acc periods to be combined
      INTEGER, SAVE :: ITmax=-1, ITmin=-1 ! to protect against long runs
      INTEGER nd1,iy1,iti1,ite1,it01,im0,jm0,lm0,ls10

      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output to end-of-month restart file
        IF (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) it,XLABEL,nday,iyear1,itimei,itimee,
     *         itime0,NTIMEACC,TIMING(1:NTIMEACC),TIMESTR(1:NTIMEACC)
C**** doc line: basic model parameters
          write(label2,'(a13,4i4,a)') 'IM,JM,LM,LS1=',im,jm,lm,ls1,' '
          label2(74:80) = 'LABEL01'
          WRITE (kunit,err=10) LABEL2
C**** write parameters database here
          call write_param(kunit)
        END IF
      CASE (IOREAD:)          ! label always from input file
!****   Determine whether timing numbers were saved as integers or reals
        read(kunit) ; read(kunit) label2 ; rewind kunit
        if (label2(80:80)==' ') then         ! integers, convert to seconds
          READ (kunit,err=10) it,XLABEL,nd1,iy1,iti1,ite1,it01,
     *        NTIM1,ITIM1(1:NTIM1),TSTR1(1:NTIM1)
          tim1(1:NTIM1) = 1d-2*itim1(1:NTIM1)
        else                                 ! real*8
          READ (kunit,err=10) it,XLABEL,nd1,iy1,iti1,ite1,it01,
     *        NTIM1,TIM1(1:NTIM1),TSTR1(1:NTIM1)
        end if
C**** use doc-record to check the basic model parameters
        READ (kunit,err=10) LABEL2
        READ (label2,'(13x,4i4)',err=10) im0,jm0,lm0,ls10
        if (im.ne.im0.or.jm.ne.jm0.or.lm.ne.lm0.or.ls10.ne.ls1) then
          ioerr = 0   ! warning
        end if
        SELECT CASE (IACTION)   ! set model common according to iaction
        CASE (ioread)           ! parameters from rundeck & restart file
          call read_param(kunit,.false.)
          nday=nd1 ; itimei=iti1      ! changeable only at new starts
          itimee=ite1 ; itime0=it01   ! are changed later if appropriate
          if (iyear1.lt.0) iyear1=iy1 ! rarely changes on restarts
          NTIMEACC=NTIM1
          TIMESTR(1:NTIM1)=TSTR1(1:NTIM1)
          TIMING(1:NTIM1)=TIM1(1:NTIM1)
        CASE (IRSFIC,irsficnt,IRSFICNO) ! rundeck/defaults except label
          read(kunit,err=10)          ! skip parameters, dates
          it=it*24/nd1                ! switch itime to ihour
        CASE (IRERUN)           ! parameters from rundeck & restart file
          call read_param(kunit,.false.)
          nday=nd1 ; itimei=iti1      ! changeable only at new starts
          itimee=ite1 ; itime0=it01   ! is changed later if appropriate
          if (iyear1.lt.0) iyear1=iy1 ! rarely changes on restart/reruns
        CASE (IOREAD_SINGLE)    ! parameters/label from 1-many acc files
          call read_param(kunit,.false.)  ! use rundeck
          call sync_param( "kradia",kradia)
          nday=nd1 ; iyear1=iy1 ; itime0=it01
          NTIMEACC=NTIM1                 ! use timing from current file
          TIMESTR(1:NTIM1)=TSTR1(1:NTIM1)
          TIMING(1:NTIM1)=TIMING(1:NTIM1)+TIM1(1:NTIM1)

C**** keep track of min/max time over the combined diagnostic period
          if (it.gt.ITmax)                   ITmax = it
          if (ITmin.lt.0 .or. it01.lt.ITmin) ITmin = it01
          itime0 = ITmin

        END SELECT
      END SELECT
      itm = max(it,ITmax)
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_label

      SUBROUTINE io_model(kunit,iaction,ioerr)
!@sum  io_model reads and writes model variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM
      USE DOMAIN_DECOMP_1D, only: grid,PACK_DATA,UNPACK_DATA,AM_I_ROOT
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "MODEL01"
!@var U_glob Work array for parallel I/O
!@var V_glob Work array for parallel I/O
!@var T_glob Work array for parallel I/O
!@var Q_glob Work array for parallel I/O
!@var WM_glob Work array for parallel I/O
#ifdef BLK_2MOM
!@var WMICE_glob Work array for parallel I/O
#endif
!@var P_glob Work array for parallel I/O
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: !(IM,JM,LM)
     &     U_glob,V_glob,T_glob,Q_glob,WM_glob
#ifdef BLK_2MOM
     *,WMICE_glob
#endif
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: P_glob!(IM,JM)
      integer :: img, jmg, lmg

      MODULE_HEADER(lhead+1:80) = 'R8 dim(im,jm,lm):u,v,t, p(im,jm),'//
     *  ' dim(im,jm,lm):q,MliqW'

      if(am_i_root()) then
         img = IM
         jmg = JM
         lmg = LM
      else
         img = 1
         jmg = 1
         lmg = 1
      end if
      allocate(U_glob(img,jmg,lmg))
      allocate(V_glob(img,jmg,lmg))
      allocate(T_glob(img,jmg,lmg))
      allocate(P_glob(IM,JM))
      allocate(Q_glob(img,jmg,lmg))
      allocate(WM_glob(img,jmg,lmg))
#ifdef BLK_2MOM
      allocate(WMICE_glob(img,jmg,lmg))
#endif

      SELECT CASE (IACTION)
      CASE (:IOWRITE) ! output to end-of-month restart file
        CALL PACK_DATA(grid, U, U_GLOB)
        CALL PACK_DATA(grid, V, V_GLOB)
        CALL PACK_DATA(grid, T, T_GLOB)
        CALL PACK_DATA(grid, Q, Q_GLOB)
        CALL PACK_DATA(grid, WM, WM_GLOB)
#ifdef BLK_2MOM
        CALL PACK_DATA(grid, WMICE, WMICE_GLOB)
#endif
        CALL PACK_DATA(grid, P, P_GLOB)
        IF (AM_I_ROOT())
     &    WRITE (kunit,err=10) MODULE_HEADER,U_glob,V_glob,T_glob,
     &                         P_glob,Q_glob,WM_glob
#ifdef BLK_2MOM
     &   ,WMICE_glob
#endif
      CASE (IOREAD:)          ! input from restart file
        if ( AM_I_ROOT() ) then
          READ (kunit,err=10) HEADER,U_glob,V_glob,T_glob,
     &                           P_glob,Q_glob,WM_glob
#ifdef BLK_2MOM
     &   ,WMICE_glob
#endif
          IF (HEADER(1:LHEAD).ne.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        end if
        CALL UNPACK_DATA(grid, U_GLOB, U)
        CALL UNPACK_DATA(grid, V_GLOB, V)
        CALL UNPACK_DATA(grid, T_GLOB, T)
        CALL UNPACK_DATA(grid, Q_GLOB, Q)
        CALL UNPACK_DATA(grid, WM_GLOB, WM)
#ifdef BLK_2MOM
        CALL UNPACK_DATA(grid, WMICE_GLOB, WMICE)
#endif
        CALL UNPACK_DATA(grid, P_GLOB, P)
      END SELECT
      call freemem
      RETURN
 10   IOERR=1
      call freemem
      RETURN
      contains
      subroutine freemem
        deallocate(U_glob)
        deallocate(V_glob)
        deallocate(T_glob)
        deallocate(P_glob)
        deallocate(Q_glob)
        deallocate(WM_glob)
#ifdef BLK_2MOM
        deallocate(WMICE_glob)
#endif
      end subroutine freemem
      END SUBROUTINE io_model

#ifdef NEW_IO
ccc was not sure where to dump these routines ... IA
      module conserv_diags
      implicit none

      contains
      subroutine declare_conserv_diags( grid, fid, name_dims )
      use domain_decomp_atm, only : dist_grid, get
      use pario, only : defvar
      implicit none
      type (dist_grid), intent(in) :: grid
      integer ::  fid
      character(len=*) :: name_dims
      integer :: i_0h, i_1h, j_0h, j_1h
      integer :: ier
      real*8, allocatable :: buf(:,:)
      call get( grid, j_strt_halo=j_0h, j_stop_halo=j_1h,
     &     i_strt_halo=i_0h, i_stop_halo=i_1h )
      allocate( buf(i_0h:i_1h,j_0h:j_1h), stat=ier)
      call defvar(grid, fid, buf, trim(name_dims))     
      deallocate( buf )
      end subroutine declare_conserv_diags

      subroutine dump_conserv_diags( grid, fid, name, conserv )
      use domain_decomp_atm, only : dist_grid, get
      use pario, only : write_dist_data
      implicit none
      type (dist_grid), intent(in) :: grid
      integer ::  fid
      character(len=*) :: name
      external :: conserv
      integer :: i_0h, i_1h, j_0h, j_1h
      integer :: ier
      real*8, allocatable :: buf(:,:)
      call get( grid, j_strt_halo=j_0h, j_stop_halo=j_1h,
     &     i_strt_halo=i_0h, i_stop_halo=i_1h )
      allocate( buf(i_0h:i_1h,j_0h:j_1h), stat=ier)
      if(grid%have_domain) call conserv(buf)
      call write_dist_data(grid,fid,trim(name),buf)
      deallocate( buf )
      end subroutine dump_conserv_diags

      end module conserv_diags

      subroutine def_rsf_model(fid)
!@sum  def_rsf_model defines U,V,T,P,Q,WM array structure in restart files
!@auth M. Kelley
!@ver  beta
      use model_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use conserv_diags
      implicit none
      integer fid   !@var fid file id
      character(len=20) :: ijlstr
      ijlstr='(dist_im,dist_jm,lm)'
      call defvar(grid,fid,u,'u'//ijlstr)
      call defvar(grid,fid,v,'v'//ijlstr)
      call defvar(grid,fid,t,'t'//ijlstr)
      call defvar(grid,fid,q,'q'//ijlstr)
      call defvar(grid,fid,wm,'wm'//ijlstr)
      call defvar(grid,fid,p,'p(dist_im,dist_jm)')
#ifdef BLK_2MOM
      call defvar(grid,fid,wmice,'wmice'//ijlstr)
#endif
      call declare_conserv_diags( grid, fid, 'watmo(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'ekatmo(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'epatmo(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'ewatmo(dist_im,dist_jm)' )
      return
      end subroutine def_rsf_model

      subroutine new_io_model(fid,iaction)
!@sum  new_io_model read/write U,V,T,P,Q,WM arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com
      use domain_decomp_atm, only: grid
      use pario, only : write_dist_data,read_dist_data
      use conserv_diags
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      external conserv_WM, conserv_KE, conserv_PE, conserv_EWM
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 'u', u)
        call write_dist_data(grid, fid, 'v', v)
        call write_dist_data(grid, fid, 't', t)
        call write_dist_data(grid, fid, 'p', p)
        call write_dist_data(grid, fid, 'q', q)
        call write_dist_data(grid, fid, 'wm', wm)
#ifdef BLK_2MOM
        call write_dist_data(grid, fid, 'wmice', wmice)
#endif
        call dump_conserv_diags( grid, fid, 'watmo', conserv_WM )
        call dump_conserv_diags( grid, fid, 'ekatmo', conserv_KE )
        call dump_conserv_diags( grid, fid, 'epatmo', conserv_PE )
        call dump_conserv_diags( grid, fid, 'ewatmo', conserv_EWM  )
      case (ioread)             ! input from restart file
        call read_dist_data(grid, fid, 'u', u)
        call read_dist_data(grid, fid, 'v', v)
        call read_dist_data(grid, fid, 't', t)
        call read_dist_data(grid, fid, 'p', p)
        call read_dist_data(grid, fid, 'q', q)
        call read_dist_data(grid, fid, 'wm', wm)
#ifdef BLK_2MOM
        call read_dist_data(grid, fid, 'wmice', wmice)
#endif
      end select
      return
      end subroutine new_io_model
#endif /* NEW_IO */

      subroutine getdte(It,Nday,Iyr0,Jyr,Jmn,Jd,Jdate,Jhour,amn)
!@sum  getdte gets julian calendar info from internal timing info
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : hrday
      USE MODEL_COM, only : JDperY,JDendOfM,amonth
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: It,Nday,Iyr0
      INTEGER, INTENT(OUT) :: Jyr,Jmn,Jd,Jdate,Jhour
      CHARACTER*4, INTENT(OUT) :: amn

      Jyr=Iyr0+It/(Nday*JDperY)
      Jd=1+It/Nday-(Jyr-Iyr0)*JDperY
      Jmn=1
      do while (Jd.GT.JDendOfM(Jmn))
        Jmn=Jmn+1
      end do
      Jdate=Jd-JDendOfM(Jmn-1)
      Jhour=mod(It*hrday/Nday,hrday)
      amn=amonth(Jmn)

      return
      end subroutine getdte

      SUBROUTINE CALC_VERT_AMP(P0,LMAX,PL,AM,PDSIG,PEDN,PMID)
!@sum  CALC_VERT_AMPK calculates air mass and pressure vertical arrays
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : bygrav
      USE MODEL_COM, only : lm,ls1,dsig,sig,sige,ptop,psfmpt,lm_req
     *     ,req_fac,req_fac_m,req_fac_d,pmtop
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: P0 !@var P0 surface pressure (-PTOP) (mb)
      INTEGER, INTENT(IN) :: LMAX !@var LMAX max level for calculation
!@var AM mass at each level (kg/m2)
!@var PDSIG pressure interval at each level (mb)
!@var PMID mid-point pressure (mb)
      REAL*8, INTENT(OUT), DIMENSION(LMAX) :: AM,PDSIG,PMID,PL
!@var PEDN edge pressure (top of box) (mb)
      REAL*8, INTENT(OUT), DIMENSION(LMAX+1) :: PEDN
      INTEGER :: L  !@var L  loop variables

C**** Calculate air mass, layer pressures
C**** Note that only layers LS1 and below vary as a function of surface
C**** pressure.
C**** Note Air mass is calculated in (kg/m^2)

      DO L=1,LS1-1
        PL(L)   = P0
        PDSIG(L)= P0*DSIG(L)
        PMID(L) = SIG(L)*P0+PTOP
        PEDN(L) = SIGE(L)*P0+PTOP
        AM  (L) = PDSIG(L)*1d2*BYGRAV
      END DO
      DO L=LS1,MIN(LMAX,LM)
        PL(L)   = PSFMPT
        PDSIG(L)= PSFMPT*DSIG(L)
        PMID(L) = SIG(L)*PSFMPT+PTOP
        PEDN(L) = SIGE(L)*PSFMPT+PTOP
        AM  (L) = PDSIG(L)*1d2*BYGRAV
      END DO
      IF (LMAX.ge.LM) PEDN(LM+1) = SIGE(LM+1)*PSFMPT+PTOP
C**** Rad. equ. layers if necessary (only PEDN,AM,PMID)
      IF (LMAX.eq.LM+LM_REQ) THEN
        PMID(LM+1:LM+LM_REQ) = REQ_FAC_M(1:LM_REQ)*PMTOP
          AM(LM+1:LM+LM_REQ) = REQ_FAC_D(1:LM_REQ)*PMTOP*1d2*BYGRAV
        PEDN(LM+2:LM+LM_REQ) = REQ_FAC(1:LM_REQ-1)*PEDN(LM+1)
        PEDN(LM+LM_REQ+1)=0.    ! 1d-5  ! why not zero?
      END IF

      RETURN
      END SUBROUTINE CALC_VERT_AMP

