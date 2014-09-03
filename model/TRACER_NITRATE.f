#include "rundeck_opts.h"
      SUBROUTINE EQSAM_DRV
      USE TRACER_COM
      USE AEROSOL_SOURCES, only: NH3_src_con, NH3_src_cyc
     & ,off_HNO3,off_SS

      USE MODEL_COM, only : im,jm,lm     ! dimensions
     $                     ,t            ! potential temperature (C)
     $                     ,q            ! saturatered pressure
     $                     ,dtsrc
      USE GEOM, only: axyp,BYAXYP
      USE CONSTANT,   only:  lhe,mair,gasc    ! latent heat of evaporation at 0 C       
      USE FLUXES, only: tr3Dsource
      USE DYNAMICS,   only: pmid,pk,byam,am   ! midpoint pressure in hPa (mb)
!                                             and pk is t mess up factor
!                                             BYAM  1/Air mass (m^2/kg)
      USE DOMAIN_DECOMP_ATM,only: GRID, GET

      IMPLICIT NONE
      ! Call parameters for the EQSAM thermodynamic model. 
      INTEGER, PARAMETER :: NCA  = 11    ! fixed number of input variables
      INTEGER, PARAMETER :: NCO  = 36    ! fixed number of output variables
      INTEGER, PARAMETER :: IOPT =  1    ! =1 selects the metastable (wet) state and history
!     INTEGER, PARAMETER :: IOPT =  2    ! =2 selects the solid      (dry) state and history
      INTEGER, PARAMETER :: LOOP =  1    ! only a single time step done
      INTEGER, PARAMETER :: IMAX =  1    ! only a single time step done

      REAL*4 :: YI(IMAX,NCA)            ! [umol/m^3] for chemical species - input
      REAL*4 :: YO(IMAX,NCO)            ! [umol/m^3] for chemical species - output
      REAL*8 :: yM,yS
      INTEGER:: j,l,i,J_0, J_1,n,I_0,I_1
C**** functions
      REAL*8 :: QSAT, AVOL
      LOGICAL, SAVE :: NO_SS = .TRUE.

      CALL GET(grid, J_STRT =J_0, J_STOP =J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

#ifndef  TRACERS_SPECIAL_Shindell
      CALL READ_OFFHNO3(OFF_HNO3)
#endif


      do n=1,ntm
      if (trname(n).eq.'seasalt1')    NO_SS=.FALSE.
      enddo   

      if (NO_SS) CALL READ_OFFSS(OFF_SS)
      YI(1,:) = 0.d0
      YO(1,:) = 0.d0

      DO L=1,LM                            
      DO J=J_0,J_1                               
      DO I=I_0,I_1
! meteo
      YI(1,1) = pk(l,i,j)*t(i,j,l)           !should be in [K]
      YI(1,2) = q(i,j,l)/QSAT (pk(l,i,j)*t(i,j,l),lhe,pmid(l,i,j)) ! rH [0-1]
      YI(1,11)= pmid(l,i,j)                  ! p in [hPa]
! calculate M and set fixed ratios for O2 & H2:
! chem trm in [kg/gb] -> molecule/cm3
! conversion molecule/cm3 -> umol/m3: 1.d6 * 1.d6 /6.022e23
      yS =  1.d6 * 1.d6 /6.022e23
      yM    = YI(1,11)/(YI(1,1)*1.38d-19) * yS
      YI(1,3) = trm(i,j,l,n_NH3)*yM*mass2vol(n_NH3)*
     *     BYAXYP(I,J)*BYAM(L,I,J) ! NH3  (g) + NH4+  (p)   [umol/m^3 air]
      YI(1,3) = YI(1,3) + (trm(i,j,l,n_NH4)*yM*mass2vol(n_NH4)*
     *     BYAXYP(I,J)*BYAM(L,I,J)) ! NH3  (g) + NH4+  (p)   [umol/m^3 air]
c#ifdef TRACERS_HETCHEM
c      YI(1,4) = (trm(i,j,l,n_SO4)+trm(i,j,l,n_SO4_d1)+
c     *          trm(i,j,l,n_SO4_d2)+trm(i,j,l,n_SO4_d3))
c     *         * yM*mass2vol(n_SO4)*
c     *         BYAXYP(I,J)*BYAM(L,I,J)    ! H2SO4    + SO4-- (p)   [umol/m^3 air]
c#else
      YI(1,4) = trm(i,j,l,n_SO4)*yM*mass2vol(n_SO4)*
     *     BYAXYP(I,J)*BYAM(L,I,J) ! H2SO4    + SO4-- (p)   [umol/m^3 air]
c#endif
#ifdef  TRACERS_SPECIAL_Shindell
      YI(1,5) = trm(i,j,l,n_HNO3)*yM*mass2vol(n_HNO3)*
     *     BYAXYP(I,J)*BYAM(L,I,J)   ! HNO3 (g) + NO3-  (p)   [umol/m^3 air]
#else 
!off-line HNO3
      YI(1,5) = off_HNO3(i,j,l)*yM!*(mair/63.018)    ! HNO3 (g)   [umol/m^3 air]
#endif
      YI(1,5) =YI(1,5)+ (trm(i,j,l,n_NO3p)*yM*mass2vol(n_NO3p)*
     *     BYAXYP(I,J)*BYAM(L,I,J) ) ! HNO3 (g) + NO3-  (p)   [umol/m^3 air]
      
      if (NO_SS) then
! estimated sea salt = NaCl
      YI(1,6) = off_SS(i,j,l)*0.5 ! off_SS [kg/kg]
     *         *yM*(mair/23.)*0.1 ! Na+ (ss  + xsod) (a)   [umol/m^3 air]

      YI(1,7) = off_SS(i,j,l)*0.5
     *         *yM*(mair/36.5)*0.1  ! HCl  (g) + Cl-   (p)   [umol/m^3 air]
      else
        YI(1,6) = (trm(i,j,l,n_seasalt1)+ trm(i,j,l,n_seasalt2))*0.5
     *       *yM*(mair/23.)*
     *       BYAXYP(I,J)*BYAM(L,I,J) *0.1 ! Na+ (ss  + xsod) (a)   [umol/m^3 air]
        YI(1,7) = (trm(i,j,l,n_seasalt1)+ trm(i,j,l,n_seasalt2))*0.5
     *       *yM*(mair/36.5)*
     *       BYAXYP(I,J)*BYAM(L,I,J)*0.1 ! HCl  (g) + Cl-   (p)   [umol/m^3 air]
      endif
#ifdef  TRACERS_DUST
! estimated after Trochkine et al. 2003, dust = 10% Ca + 10% K + 20% Mg +[60% (Na, Al, Si, Fe)]
      YI(1,8) = (trm(i,j,l,n_Clay)+trm(i,j,l,n_Silt1)+
     *         trm(i,j,l,n_Silt2)+trm(i,j,l,n_Silt3))*0.1
     *         *yM*(mair/39.1)*
     *         BYAXYP(I,J)*BYAM(L,I,J)   ! K+   (p) from Dust     [umol/m^3 air]
      YI(1,9) = (trm(i,j,l,n_Clay)+trm(i,j,l,n_Silt1)+
     *         trm(i,j,l,n_Silt2)+trm(i,j,l,n_Silt3))*0.1 
     *         *yM*(mair/40.)*
     *         BYAXYP(I,J)*BYAM(L,I,J) ! Ca++ (p) from Dust     [umol/m^3 air]
      YI(1,10)= (trm(i,j,l,n_Clay)+trm(i,j,l,n_Silt1)+
     *         trm(i,j,l,n_Silt2)+trm(i,j,l,n_Silt3))*0.2
     *         *yM*(mair/24.3)*
     *         BYAXYP(I,J)*BYAM(L,I,J) ! Mg++ (p) from Dust     [umol/m^3 air]  
#endif


      call EQSAM_V03D(YI,YO,NCA,NCO,IOPT,LOOP,IMAX,66)

      YO(1,12) = MAX(YO(1,12),1.D-30)
      YO(1,20) = MAX(YO(1,20),1.D-30)
      YO(1,10) = MAX(YO(1,10),1.D-30)
      YO(1,19) = MAX(YO(1,19),1.D-30)
      YO(1,9)  = MAX(YO(1,9),1.D-30)
! Nitrate production   
      tr3Dsource(i,j,l,1,n_NO3p)= ((YO(1,20)
     *        /(yM*mass2vol(n_NO3p)* BYAXYP(I,J)*BYAM(L,I,J)) )
     *        -trm(i,j,l,n_NO3p)) /dtsrc
! Ammonia residual
      tr3Dsource(i,j,l,1,n_NH3) = ((YO(1,10)
     *        /(yM*mass2vol(n_NH3)* BYAXYP(I,J)*BYAM(L,I,J)) )
     *        -trm(i,j,l,n_NH3)) /dtsrc
! Ammonium production
      tr3Dsource(i,j,l,1,n_NH4) = ((YO(1,19) 
     *        /(yM*mass2vol(n_NH4)* BYAXYP(I,J)*BYAM(L,I,J)) )
     *        -trm(i,j,l,n_NH4))/dtsrc
! Aerosol Water [ug/m3]

c aval [m3/gb] mass of air pro m3      
c      AVOL = am(l,i,j)*axyp(i,j)/mair*1000.d0*gasc*
c     *      (pk(l,i,j)*t(i,j,l)) /(pmid(l,i,j)     *100.)
c      trm(i,j,l,n_AW) = YO(1,12) *AVOL * 1.d-9

#ifdef  TRACERS_SPECIAL_Shindell
! Nitric Acid residual
      tr3Dsource(i,j,l,3,n_HNO3) =((YO(1,9)
     *        /(yM*mass2vol(n_HNO3)* BYAXYP(I,J)*BYAM(L,I,J)) )
     *        -trm(i,j,l,n_HNO3))/dtsrc
#endif

      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE EQSAM_DRV
