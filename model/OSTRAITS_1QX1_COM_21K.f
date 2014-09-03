#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_WATER
#endif

      Module STRAITS
C****
!@sum   STRAITS ocean strait related variables
!@+     RESOLUTION DEPENDENT: This version is for 288x180 - S
!@auth  Gary Russell/Gavin Schmidt
!@ver   2009/05/26
C****
      Use OCEAN,  Only: LMO
      Use SEAICE, Only: LMI

#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm_atm=>ntm
#endif
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : ntm
#endif
      Implicit None
      Save
C****
C****     Strait              From         To         LM    Width
C****     ------              ----         --         --    -----
C****  X  Dolphin & Union   51,160 ES    53,159  N     3    32000
C****  X  Dease             57,159 EN    61,159 WN     3    25000
C****  X  Fury & Hecla      76,160 EN    79,160 W      2    15000
C****  X  Nares             91,171 EN    93,172 WS     6    30000
C****  1  Gibraltar        139,126 EN   141,126 WS     7    15000
C****  X  English          145,141 EN   146,142 WS     2    35000
C****  X  Dardanelles      165,130 EN   167,131 WS     2     1000
C****  X  Bosporous        167,131 EN   168,132 WS     2     2000
C****  2  Bab al Mandab    179,104 ES   180,103 WN     5    25000
C****  X  Malacca          225, 93 ES   227, 92 WN     2    40000
C****  X  Salat Sunda      229, 84 EN   230, 85 WS     2    25000
C****
C**** Fixed parameters of straits
C****      NMST    = number of ocean straits
C****  IST(NMST,2) = I index of west and east ocean cells
C****  JST(NMST,2) = J index of west and east ocean cells
C**** LMST(NMST)   = number of layers in strait
C****  XST(NMST,2) = X location of strait end in ocean cell [-1:1]
C****  YST(NMST,2) = Y location of strait end in ocean cell [-1:1]
C**** WIST(NMST)   = width of strait (m)
C****
      Integer*4,Parameter :: NMST =2,
     *  IST(NMST,2) = Reshape ( (/
     *     139,179,
     *     141,180/),(/NMST,2/) ),
     *  JST(NMST,2) = Reshape ( (/
     *    126,104,
     *    126,103/),(/NMST,2/) ),
     *  LMST(NMST) = (/ 6, 2/)
      Real*8,Parameter ::
     *  XST(NMST,2) = Reshape ( (/
     *      .7d0, .7d0,
     *     -.7d0,-.7d0 /), (/ NMST,2 /) ),
     *  YST(NMST,2) = Reshape ( (/
     *      .7d0,-.7d0,
     *     -.7d0, .7d0/), (/ NMST,2 /) ),
     *  WIST(NMST) = (/ 15d3, 25d3 /)
      Character*20,Parameter :: NAME_ST(NMST) = (/
     *     'Gibraltar   ', 'Bab al Mandab' /)

      Real*8 ::
     *   MMST(LMO,NMST), !  mass of water in strait (kg)
     *   MUST(LMO,NMST), !  mass flux of water in strait (kg/s)
     *  G0MST(LMO,NMST), !  pot. enthalpy of water in strait (J)
     *  GXMST(LMO,NMST), !  west-east gradient of pot. enthalpy (J)
     *  GZMST(LMO,NMST), !  down. vert. gradient of pot. enthalpy (J)
     *  S0MST(LMO,NMST), !  mass of salt in strait (kg)
     *  SXMST(LMO,NMST), !  west-east gradient of salt (kg)
     *  SZMST(LMO,NMST), !  down. vert. gradient of salt (kg)
     *       DIST(NMST), !  fixed length of strait (m)
     *     DISTPG(NMST), !  fixed distance between ocean cell centers
     *      RSIST(NMST), !  horizontal sea ice cover in strait
     *     RSIXST(NMST), !  center of sea ice cover in strait
     *    MSIST(2,NMST), !  2 mass layers of sea ice in strait (kg)
     *  HSIST(LMI,NMST), !  LMI layers of heat content in strait (J)
     *  SSIST(LMI,NMST), !  LMI layers of salt in strait (kg)
     *    USIFAC = .1d0  !  ratio of strait sea ice velocity to current

!@var QTYE workspace holding values of QTY at endpoints of straits
      Real*8 :: OPRESE(2,NMST)
      Real*8, Dimension(2,NMST,LMO) ::
     *       MOE, G0ME,GXME,GYME,GZME, S0ME,SXME,SYME,SZME

!@var kn2: Where kn2>0, the index pair k,n [k=1,2;n=1,nmst] corresponds
!@+        to the same i,j endpoint as index pair kn2(1:2,k,n).
!@+        Currently, an i,j endpoint can be shared by only 2 straits.
      integer, dimension(2,2,nmst) :: kn2

#ifdef TRACERS_OCEAN
      Real*8 TRMST(LMO,NMST,NTM),TXMST(LMO,NMST,NTM),TZMST(LMO,NMST,NTM)
!@var TRME,TXME,TYME,TZME tracers at the endpoints of straits
      Real*8, Dimension(2,NMST,LMO,NTM) :: TRME,TXME,TYME,TZME
#endif
#ifdef TRACERS_WATER
      Real*8 TRSIST(NTM_ATM,LMI,NMST)
#endif

!@param NMST no. of ocean straits
!@var IST,JST i,j coordinates of ends of straits
!@var LMST no. of levels in strait
!@var XST,YST local coordinates [-1,1] for strait entry/exit points
!@var WIST width of strait (m)
!@var NAME_ST Names of straits
!@var MMST mass of water in strait (kg)
!@var MUST mass flux of water in strait (kg/s)
!@var G0MST,GXMST,GZMST pot. enthalpy of water in strait (+ moments) (J)
!@var S0MST,SXMST,SZMST salinity of water in strait (+ moments) (kg)
!@var DIST distance along strait (m)
!@var DISTPG distance between centre points of adjoining ocean boxes (m)
!@var RSIST Sea ice fraction in strait
!@var RSIXST Center of sea ice in strait (m)
!@var MSIST Mass of ice within strait (kg)
!@var HSIST Enthalpy of ice within strait (J)
!@var SSIST Salinity of ice within strait (kg)
!@param USIFAC ratio of strait sea ice velocity to current
!@var TRMST,TXMST,TZMST tracer amount in strait (+ moments) (kg)
!@var TRSIST tracer amount in with strait (kg)
      EndModule STRAITS
