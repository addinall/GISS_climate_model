#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_WATER
#endif
      MODULE STRAITS
!@sum  STRAITS ocean strait related variables
!@+    RESOLUTION DEPENDENT: This version for 72x46 - M, 9ky topography
!@+    Fury+Hecla blocked.
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm_atm=>ntm
#endif
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : ntm
#endif
      USE SEAICE, only : lmi
      USE OCEAN, only : lmo
      IMPLICIT NONE
      SAVE
C**** These values are highly resolution dependent
C****
C****     Strait           From         To       LM    Width
C****     ------           ----         --       --    -----
C****  x  Fury & Hecla   19,42 ES    20,40 WN     2    20000
C****  1  Nares          22,43 EN    24,44 WS     5     5000
C****  2  Gibraltar      35,32 EN    37,33 WS     5    25000
C****  3  English        36,36 EN    37,37 WS     2    35000
C****  4  Kattegat       38,38 EN    40,38 WS     2    60000
C****  5  Bosporous      42,33 EN    43,34 WS     2     6000
C****  6  Red Sea        44,29 ES    45,28 WN     6   250000
C****  7  Bab al Mandab  45,28 ES    46,27 WN     6    25000
C****  8  Hormuz         47,30 ES    49,29 WN     2   100000  ! 50000
C****  9  Malacca        56,25 EN    58,24 WS     3    50000
C**** 10  Korea          62,32 EN    63,33 WS     4   170000
C**** 11  Soya-kaikyo    64,34 EN    65,35 WS     2    40000
C**** 11alt Tsugaru      64,34 EN    66,34 WS     4    20000
C****
      INTEGER, PARAMETER :: NMST=11  !@param NMST no. of ocean straits
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
     *     WIST = (/         5d3, 2.5d4, 3.5d4,   6d4, 6d3,
     *              2.5d5, 2.5d4,   1d5,   5d4, 1.7d5, 4d4/)

!@var XST,YST local coordinates [-1,1] for strait entry/exit points
      REAL*8, DIMENSION(NMST,2) :: XST = RESHAPE( (/
     *            6d-1, 6d-1,  1d0,  1d0,  0d0, 6d-1, 6d-1,  1d0, 6d-1,
     *     0d0  , 6d-1,
     *           -8d-1,-8d-1,  0d0,-6d-1,-8d-1,-8d-1, -1d0,-7d-1, -1d0,
     *     -6d-1, -1d0 /), (/NMST,2/) ),
     *     YST = RESHAPE( (/
     *            8d-1, 8d-1,  0d0,  0d0,  1d0,-8d-1,-8d-1,  0d0,-8d-1,
     *      1d0 , 8d-1,
     *           -6d-1,-6d-1, -1d0,-8d-1,-6d-1, 6d-1,  0d0, 7d-1,  0d0,
     *     -8d-1,  0d0 /), (/NMST,2/) )

!@var IST,JST i,j coordinates of ends of straits
      INTEGER, DIMENSION(NMST,2) ::
     *     IST = RESHAPE( (/
     *         22, 35, 36, 38, 42, 44, 45, 47, 56, 62, 64,
     *         24, 37, 37, 40, 43, 45, 46, 49, 58, 63, 65/),
     *     (/NMST,2/) ),
     *     JST = RESHAPE( (/
     *         43, 32, 36, 38, 33, 29, 28, 30, 25, 32, 34,
     *         44, 33, 37, 38, 34, 28, 27, 29, 24, 33, 35/),
     *     (/NMST,2/) )

!@var kn2: Where kn2>0, the index pair k,n [k=1,2;n=1,nmst] corresponds
!@+        to the same i,j endpoint as index pair kn2(1:2,k,n).
!@+        Currently, an i,j endpoint can be shared by only 2 straits.
      integer, dimension(2,2,nmst) :: kn2

!@var LMST no. of levels in strait
      INTEGER, DIMENSION(NMST) :: LMST = (/
     *         5,  5,  2,  2,  2,  6,  6,  2,  3,  4,  2/)

!@var name_st Names of straits
      CHARACTER*20, DIMENSION(NMST) :: NAME_ST = (/
     *                       'Nares         ', 'Gibraltar     ',
     *     'English       ', 'Kattegat      ', 'Bosporous     ',
     *     'Red Sea       ', 'Bab al Mandab ', 'Hormuz        ',
     *     'Malacca       ', 'Korea         ', 'Soya-kaikyo   '/)

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
      REAL*8, DIMENSION(NTM_ATM,LMI,NMST) :: TRSIST
#endif

      END MODULE STRAITS
