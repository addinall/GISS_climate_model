!@sum  RES_F79 Resolution info for 79 layer, 2x2.5 model
!@+    (top at .1 mb, no GWDRAG)
!@auth Original Development Team
!@ver  1.0

      MODULE RESOLUTION
!@sum  RESOLUTION contains horiz/vert resolution variables
!@auth Original Development Team
!@ver  1.0
      IMPLICIT NONE
      SAVE
!@var IM,JM longitudinal and latitudinal number of grid boxes
!@var LM number of vertical levels
!@var LS1 Layers LS1->LM: constant pressure levels, L<LS1: sigma levels
      INTEGER, PARAMETER :: IM=144,JM=90,LM=79, LS1=24

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
      REAL*8, PARAMETER :: PSF=984.d0, PTOP = 150.d0, PMTOP = .1d0
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
      REAL*8, PARAMETER :: PSFMPT=PSF-PTOP, PSTRAT=PTOP-PMTOP

!@var PLbot pressure levels at bottom of layers (mb)
      REAL*8, PARAMETER, DIMENSION(LM+1) :: PLbot = (/
     *     PSF,   964d0, 942d0, 917d0, 890d0, 860d0, 825d0,        ! Pbot L=1,..   
     *     785d0, 740d0, 692d0, 642d0, 591d0, 539d0, 489d0,        !      L=...    
     *     441d0, 396d0, 354d0, 316d0, 282d0, 251d0, 223d0,  
     *     197d0, 173d0,                                    
     *     PTOP,                                                   !      L=LS1    
     *            128d0,       108d0,       100d0, 95.804824d0,    !      L=...    
     *      91.785644d0, 87.935075d0, 84.246044d0, 80.711775d0, 
     *      77.325774d0, 74.081822d0, 70.973960d0, 67.996477d0,  
     *      65.143906d0, 62.411005d0, 59.792753d0, 57.284342d0,  
     *      54.881164d0, 52.578802d0, 50.373029d0, 48.259792d0,  
     *      46.235209d0, 44.295561d0, 42.437285d0, 40.656966d0,
     *           38.9d0,     37.16d0,     35.43d0,     33.71d0,  
     *             32d0,      30.3d0,     28.61d0,     26.93d0,
     *          25.26d0,      23.6d0,     21.95d0,     20.31d0,  
     *          18.69d0,     17.09d0,     15.51d0,     13.96d0,  
     *          12.44d0,     10.99d0,      9.62d0,      8.36d0, 
     *           7.12d0,       5.9d0,       4.7d0,      3.59d0, 
     *          2.697d0,      1.99d0,       1.4d0,     0.955d0, 
     *            0.6d0,     0.346d0,      0.19d0,     PMTOP /)    !      L=..,LM+1  
 
C**** KEP depends on whether stratos. EP flux diagnostics are calculated
C**** If dummy EPFLUX is used set KEP=0, otherwise KEP=21
!@param KEP number of lat/height E-P flux diagnostics
      INTEGER, PARAMETER :: KEP=0

C**** Based on model top, determine how much of stratosphere is resolved
C****         PMTOP >= 10 mb,    ISTRAT = 0
C**** 1 mb <= PMTOP <  10 mb,    ISTRAT = 1
C****         PMTOP <   1 mb,    ISTRAT = 2
      INTEGER, PARAMETER :: ISTRAT = 2

      END MODULE RESOLUTION

C**** The vertical resolution also determines whether
C**** stratospheric wave drag will be applied or not.
C**** Hence also included here are some dummy routines for non-strat
C**** models.

      SUBROUTINE DUMMY_STRAT
!@sum DUMMY dummy routines for non-stratospheric models
C**** Dummy routines in place of STRATDYN
      ENTRY init_GWDRAG
      ENTRY GWDRAG
      ENTRY VDIFF
      ENTRY io_strat
      ENTRY alloc_strat_com
C**** Dummy routines in place of STRAT_DIAG (EP flux calculations)
C**** Note that KEP=0 is set to zero above for the dummy versions.
      ENTRY EPFLUX
      ENTRY EPFLXI
      ENTRY EPFLXP
C****
      RETURN
      END SUBROUTINE DUMMY_STRAT

