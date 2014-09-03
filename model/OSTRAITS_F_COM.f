#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_WATER
#endif

      MODULE STRAITS
C****
!@sum  STRAITS ocean strait related variables
!@+    RESOLUTION DEPENDENT: This version is for 144x90 - F
!@auth Gary Russell/Gavin Schmidt
!@ver  2009/05/26
C****
      USE SEAICE, only : lmi
      USE OCEAN, only : lmo
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
      USE OCN_TRACER_COM, only : ntm
#endif
      IMPLICIT NONE

      SAVE
C**** These values are highly resolution dependent
C****
C****     Strait           From         To       LM    Width
C****     ------           ----         --       --    -----
C****  1  Fury & Hecla   38,80 ES    40,80 WN     2    20000
C****  2  Nares          45,85 EN    47,86 WS     5     5000
C****  3  Belle Isle     49,71 EN    51,71 WS     4    40000
C****  4  Gibraltar      70,63 EN    72,64 WS     5    25000
C****  5  Irish Sea      70,72 NW    70,74 ES     3    30000
C****  6  Kattegat       77,74 EN    78,73 WS     2    60000
C****  7  Bosporous      83,65 EN    84,66 WS     2     6000
C****  8  White Sea      88,78 EN    89,79 WS     2    40000
C****  9  Bab al Mandab  89,53 ES    90,52 WN     6    25000
C**** 10  Hormuz         94,58 ES    96,58 WN     2   100000  ! 50000
C**** 11  Malacca       112,48 EN   114,46 WS     2    50000
C**** 12  Soya-kaikyo   128,68 EN   129,69 WS     2    40000
C****
      INTEGER, PARAMETER :: NMST=12  !@param NMST no. of ocean straits
!@var MMST mass of water in strait (kg)
!@var MUST mass flux of water in strait (kg/s)
!@var G0MST,GXMST,GYMST pot. enthalpy of water in strait (+ moments) (J)
!@var S0MST,SXMST,SYMST salinity of water in strait (+ moments) (kg)
      REAL*8, DIMENSION(LMO,NMST) :: MMST,MUST,G0MST,GXMST,GZMST,S0MST
     *     ,SXMST,SZMST

!@var QTYE workspace holding values of QTY at endpoints of straits
      Real*8 :: OPRESE(2,NMST)
      Real*8, Dimension(2,NMST,LMO) ::
     *       MOE, G0ME,GXME,GYME,GZME, S0ME,SXME,SYME,SZME

!@var WIST width of strait (m)
!@var DIST distance along strait (m)
!@var DISTPG distance between centre points of adjoining ocean boxes (m)
      REAL*8, DIMENSION(NMST) :: DIST,DISTPG,
     *     WIST = (/  2d4,  5d3,  4d4, 2.5d4,   3d4,  6d4,
     *                6d3,  4d4,2.5d4,   1d5,   5d4,  4d4 /)

!@var XST,YST local coordinates [-1,1] for strait entry/exit points
      REAL*8, DIMENSION(NMST,2) :: XST = RESHAPE( (/
     *     .6d0,   0d0,   1d0,  .8d0,  .6d0,  .8d0,  .8d0,  .8d0,
     *     .8d0,   1d0,  .8d0,  .8d0, -.6d0,  -1d0,  -1d0, -.8d0,
     *      0d0, -.8d0, -.8d0, -.8d0, -.8d0,  -1d0, -.8d0, -.8d0
     *     /), (/NMST,2/) ),
     *     YST = RESHAPE( (/
     *     .8d0,   1d0,   0d0,  .8d0,   1d0, -.8d0,  .8d0,  .8d0,
     *    -.8d0,   0d0, -.8d0,  .8d0,  .8d0,   0d0,  .6d0, -.8d0,
     *     -1d0,  .8d0, -.8d0, -.8d0,  .8d0,   0d0,  .8d0, -.8d0
     *     /), (/NMST,2/) )

!@var IST,JST i,j coordinates of ends of straits
      INTEGER, DIMENSION(NMST,2) ::
     *     IST = RESHAPE( (/
     *     38, 45, 49, 70, 70, 77, 83, 88, 89, 94,112,128,
     *     40, 47, 51, 72, 70, 78, 84, 89, 90, 96,114,129/),
     *     (/NMST,2/) ),
     *     JST = RESHAPE( (/
     *     80, 85, 71, 63, 72, 74, 65, 78, 53, 58, 48, 68,
     *     80, 86, 71, 64, 74, 73, 66, 79, 52, 58, 46, 69/),
     *     (/NMST,2/) )

!@var kn2: Where kn2>0, the index pair k,n [k=1,2;n=1,nmst] corresponds
!@+        to the same i,j endpoint as index pair kn2(1:2,k,n).
!@+        Currently, an i,j endpoint can be shared by only 2 straits.
      integer, dimension(2,2,nmst) :: kn2

!@var LMST no. of levels in strait
      INTEGER, DIMENSION(NMST) :: LMST = (/
     *      2,  5,  4,  5,  3,  2,  2,  2,  6,  2,  2,  2/)

!@var name_st Names of straits
      CHARACTER*20, DIMENSION(NMST) :: NAME_ST = (/
     *     'Fury & Hecla  ', 'Nares         ', 'Belle Isle    ',
     *     'Gibraltar     ', 'Irish Sea     ', 'Kattegat      ',
     *     'Bosporous     ', 'White Sea     ', 'Bab al Mandab ',
     *     'Hormuz        ', 'Malacca       ', 'Soya-kaikyo   '/)

!@var RSIST Sea ice fraction in strait
!@var RSIXST Center of sea ice in strait (m)
!@var MSIST Mass of ice within strait (kg)
!@var HSIST Enthalpy of ice within strait (J)
!@var SSIST Salinity of ice within strait (kg)
      REAL*8, DIMENSION(NMST) :: RSIST,RSIXST
      REAL*8, DIMENSION(2,NMST) :: MSIST
      REAL*8, DIMENSION(LMI,NMST) :: HSIST,SSIST

!@param USIFAC ratio of strait sea ice velocity to current
      REAL*8 :: USIFAC = 0.1d0   ! used to be 1. (too much)

#ifdef TRACERS_OCEAN
!@var TRMST,TXMST,TZMST tracer amount in strait (+ moments) (kg)
      REAL*8, DIMENSION(LMO,NMST,NTM) :: TRMST, TXMST, TZMST
!@var TRME,TXME,TYME,TZME tracers at the endpoints of straits
      Real*8, Dimension(2,NMST,LMO,NTM) :: TRME,TXME,TYME,TZME
#endif
#ifdef TRACERS_WATER
!@var TRSIST tracer amount in with strait (kg)
      REAL*8, DIMENSION(NTM,LMI,NMST) :: TRSIST
#endif

      END MODULE STRAITS
