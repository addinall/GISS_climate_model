#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_WATER
#endif
      MODULE STRAITS
!@sum  STRAITS ocean strait related variables
!@+    RESOLUTION DEPENDENT: This version is for 72x46 - M
!@auth Gary Russell/Gavin Schmidt/Allegra LeGrande
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
C**** LGM 21 kyr
C****     Strait           From         To       LM    Width
C****     ------           ----         --       --    -----
C****  X  Fury & Hecla  Above ground/under ice
C****  X  Nares         Above ground/under ice
C****  1  Gibraltar      35,32 EN    37,33 WS     5    25000
C****  X  English       Above ground/under ice
C****  X  Kattegat      Above ground/under ice
C****  X  Bosporous     Black Sea -> Lake
C****  2  Red Sea        44,29 ES    45,28 WN     5   250000
C****  3  Bab al Mandab  45,28 ES    46,27 WN     2    25000
C****  X  Hormuz         Above ground
C****  X  Malacca        Above ground
C****  4  Korea          63,31 EN    63,33 WS     4    17000
C****                   deep salty (Tsushima gone)     85000
C****  X  Soya-kaikyo    Above ground

      INTEGER, PARAMETER :: NMST=4  !@param NMST no. of ocean straits
      integer, parameter :: n_STgbrltr=1 ,n_STrs=2,n_STbam=3,n_STkor=4
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
     *     WIST = (/  2.5d4, 2.5d5, 2.5d4, 1.7d5/)

!@var XST,YST local coordinates [-1,1] for strait entry/exit points
      REAL*8, DIMENSION(NMST,2) :: XST = RESHAPE( (/
     *     6d-1,  6d-1, 6d-1, -1d0,
     *     -8d-1,-8d-1, -1d0, -6d1 /), (/NMST,2/) ),
     *     YST = RESHAPE( (/
     *     8d-1,  -8d-1,-8d-1, 1d0,
     *     -6d-1, 6d-1,  0d0, -8d-1 /), (/NMST,2/) )

!@var IST,JST i,j coordinates of ends of straits
      INTEGER, DIMENSION(NMST,2) ::
     *     IST = RESHAPE( (/
     *     35, 44, 45, 63,
     *     37, 45, 46, 63 /), (/NMST,2/) )
     *     , JST = RESHAPE( (/
     *     32, 29, 28, 31,
     *     33, 28, 27, 33 /), (/NMST,2/) )

!@var kn2: Where kn2>0, the index pair k,n [k=1,2;n=1,nmst] corresponds
!@+        to the same i,j endpoint as index pair kn2(1:2,k,n).
!@+        Currently, an i,j endpoint can be shared by only 2 straits.
      integer, dimension(2,2,nmst) :: kn2

!@var LMST no. of levels in strait
      INTEGER, DIMENSION(NMST) :: LMST = (/
     *     5,   5,  2,  4 /) ! 180m,144,17,100 ... given 120m slr to 0k

!@var name_st Names of straits
      CHARACTER*20, DIMENSION(NMST) :: NAME_ST = (/
     *     'Gibraltar     ', 'Red Sea       '
     *     ,'Bab al Mandab ', 'Korea         '/)

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
