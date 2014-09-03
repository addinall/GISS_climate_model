#include "rundeck_opts.h"

#ifdef TRACERS_ON
!@sum  TRACERS: generic tracer routines used for all tracers
!@+    Routines included:
!@+      Generic diags: set_generic_tracer_diags
!@+      Apply previously set sources: apply_tracer_sources
!@+      Radioactive Decay: tdecay
!@+      Gravitaional Settling: trgrav
!@+      Check routine: checktr
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      SUBROUTINE set_generic_tracer_diags
!@sum set_generic_tracer_diags init trace gas attributes and diagnostics
!@auth J. Lerner
!@calls sync_param
      USE CONSTANT, only: mair,sday
      USE MODEL_COM, only: dtsrc,nisurf
      USE DIAG_COM, only: ia_src,ia_12hr,ir_log2,ir_0_71
      USE TRACER_COM
      USE TRDIAG_COM
      USE PARAM
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      implicit none
      integer :: l,k,n,nd
      character*10, DIMENSION(NTM) :: CMR,CMRWT
      logical :: T=.TRUE. , F=.FALSE.
      character*50 :: unit_string

C**** Get factor to convert from mass to volume mr if necessary
      do n=1,ntm
        if (to_volume_MixRat(n) .eq.1) then
          MMR_to_VMR(n) = mass2vol(n)
          cmr(n) = 'V/V air'
        else
          MMR_to_VMR(n) = 1.d0
          cmr(n) = 'kg/kg air'
        endif
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) then
          cmrwt(n) = 'per mil'
          cmr(n) = 'per mil'     ! this overrides to_volume_ratio
          MMR_to_VMR(n) = 1.
        else
          cmrwt(n) = 'kg/m^2'
        end if
#endif
      end do
C****
C**** TAJLN(J,L,KQ,N)  (SUM OVER LONGITUDE AND TIME OF)
C****
C**** jlq_power Exponent associated with a physical process
C****      (for printing tracers).   (ntm_power+jlq_power)
C**** jls_power Exponent associated with a source/sink (for printing)
C**** Defaults for JLN
      scale_jlq(:) = 1./DTsrc
      jgrid_jlq(:) = 1
      ia_jlq(:) = ia_src

C**** Tracer concentration
      do n=1,ntm
        k = 1        ! <<<<< Be sure to do this
        jlnt_conc = k
        sname_jln(k,n) = trim(trname(n))//'_CONCENTRATION'
        lname_jln(k,n) = trim(trname(n))//' CONCENTRATION'
        jlq_power(k) = 0.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),cmr(n))
        scale_jlq(k) = 1.d0
        scale_jln(n) = MMR_to_VMR(n)
#ifdef TRACERS_WATER
        if (to_per_mil(n).gt.0) units_jln(k,n) = unit_string(0,cmr(n))
#endif
C**** Tracer mass
        k = k + 1
        jlnt_mass = k
        sname_jln(k,n) = trim(trname(n))//'_MASS'
        lname_jln(k,n) = trim(trname(n))//' MASS'
        jlq_power(k) = 4
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP)
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k)+13
     *       ,'kg')
#else
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k)
     *       ,'kg/m^2')
#endif
        scale_jlq(k) = 1.d0
#ifdef TRACERS_WATER
C****   TRACER CONCENTRATION IN CLOUD WATER
        k = k + 1
        jlnt_cldh2o = k
        sname_jln(k,n) = trim(trname(n))//'_WM_CONC'
        lname_jln(k,n) = trim(trname(n))//' CLOUD WATER CONCENTRATION'
        jlq_power(k) = 4
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k)
     *       ,'kg/kg water')
        scale_jlq(k) = 1.d0
        if (to_per_mil(n).gt.0) units_jln(k,n) = unit_string(0,cmr(n))
#endif
C**** Physical processes affecting tracers
C****   F (TOTAL NORTHWARD TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_nt_tot = k
        sname_jln(k,n) = 'tr_nt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL NORTHWARD TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11
        jgrid_jlq(k) = 2
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   STM/SM (MEAN MERIDIONAL N.T. OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_nt_mm = k
        sname_jln(k,n) = 'tr_nt_mm_'//trname(n)
        lname_jln(k,n) = 'NORTHWARD TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10
        jgrid_jlq(k) = 2
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   F (TOTAL VERTICAL TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_tot = k
        sname_jln(k,n) = 'tr_vt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL VERTICAL TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   STM/SM (MEAN MERIDIONAL V.T. OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_mm = k
        sname_jln(k,n) = 'tr_vt_mm_'//trname(n)
        lname_jln(k,n) = 'VERTICAL TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   TMBAR-TM (CHANGE OF TRACER MASS BY MOIST CONVEC)(kg)
        k = k + 1
        jlnt_mc = k
        sname_jln(k,n) = 'tr_mc_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY MOIST CONVECTION'
        jlq_power(k) = 10
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   TMBAR-TM (CHANGE OF TRACER MASS BY Large-scale CONDENSE)  (kg)
        k = k + 1
        jlnt_lscond = k
        sname_jln(k,n) = 'tr_lscond'//trname(n)
        lname_jln(k,n) ='CHANGE OF '//
     &     trim(trname(n))//' MASS BY LARGE-SCALE CONDENSE'
        jlq_power(k) = 10.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   TMBAR-TM (CHANGE OF TRACER MASS BY DRY CONVEC)  (kg)
        k = k + 1
        jlnt_turb = k
        sname_jln(k,n) = 'tr_turb_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY TURBULENCE/DRY CONVECTION'
        jlq_power(k) = 10
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')


        if (k.gt. ktajl) then
           if (AM_I_ROOT()) write (6,*)
     &          'tjl_defs: Increase ktajl=',ktajl,' to at least ',k
           call stop_model('ktajl too small',255)
        end if
      end do

C**** CONTENTS OF TAIJLN(I,J,LM,N)  (SUM OVER TIME OF)
C****        TML (M*M * KG TRACER/KG AIR)
C**** Set defaults that are true for all tracers and layers
      do n=1,ntm
        ijtm_power(n) = ntm_power(n)+4 !n for integrated mass
        ijtc_power(n) = ntm_power(n)+1 !n for concentration
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) ijtc_power(n) = 0
#endif
      end do
C**** Tracer concentrations (TAIJLN)
      do n=1,ntm
        write(sname_ijt(n),'(a)') trim(TRNAME(n))
        write(lname_ijt(n),'(a)') trim(TRNAME(n))
        if (to_conc(n).eq.1) then   ! diag in kg/m3
          units_ijt(n) = unit_string(ijtc_power(n),'kg/m3')
          scale_ijt(n) = 10.**(-ijtc_power(n))
        else ! mixing ratio
          units_ijt(n) = unit_string(ijtc_power(n),cmr(n))
          scale_ijt(n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
        end if
      end do

C**** AIJN
C****     1  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER/KG AIR)
C****     2  TRS (SURFACE TRACER CONC.) (M*M * KG TRACER/KG AIR)
C****     3  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER)
      do n=1,ntm
C**** Summation of mass over all layers
      k = 1        ! <<<<< Be sure to do this
      tij_mass = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Total_Mass'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Total Mass'
        units_tij(k,n) = unit_string(ijtm_power(n),'kg/m^2')
        scale_tij(k,n) = 10.**(-ijtm_power(n))
C**** Average concentration over layers
      k = k+1
      tij_conc = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Average'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Average'
        units_tij(k,n) = unit_string(ijtc_power(n),cmr(n))
        scale_tij(k,n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) then
          denom_tij(k,n)=n_Water
          scale_tij(k,n)=1.
        endif
#endif
C**** Surface concentration
      k = k+1
      tij_surf = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_At_Surface'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' At Surface'
        units_tij(k,n) = unit_string(ijtc_power(n),cmr(n))
        scale_tij(k,n)=MMR_to_VMR(n)*10.**(-ijtc_power(n))/
     *                 REAL(NIsurf,KIND=8)
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) then
          denom_tij(k,n)=n_Water
          scale_tij(k,n)=1.
        endif
#endif
C**** Surface concentration by volume (units kg/m^3)
      k = k+1
      tij_surfbv = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *     '_byVol_At_Surface'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' By Vol At Surface'
        units_tij(k,n) = unit_string(ijtc_power(n),'kg/m^3')
        scale_tij(k,n)=MMR_to_VMR(n)*10.**(-ijtc_power(n))/
     *                 REAL(NIsurf,KIND=8)
#ifdef TRACERS_WATER
C**** the following diagnostics are set assuming that the particular
C**** tracer exists in water.
C**** Tracers in precipitation (=Wet deposition)
      k = k+1
      tij_prec = k
      if (dowetdep(n)) then
        if (to_per_mil(n) .eq.1) then
          write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_prec'
          write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *         ' in Precip'
          units_tij(k,n)=cmrwt(n)
          scale_tij(k,n)=1.
          denom_tij(k,n)=n_Water
        else
          write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_wet_dep'
          write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *         ' Wet Deposition'
          units_tij(k,n)=unit_string(ijtc_power(n)-5,trim(cmrwt(n))
     *         //'/s')
          scale_tij(k,n)=10.**(-ijtc_power(n)+5)/dtsrc
        end if
      end if
C**** Tracers in evaporation
      k = k+1
      tij_evap = k
      if (tr_wd_type(n).eq.nWater) then
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_evap'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in Evaporation'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(ijtc_power(n),cmrwt(n))
          scale_tij(k,n)=1.
          denom_tij(k,n)=n_Water
        else
          units_tij(k,n)=unit_string(ijtc_power(n),trim(cmrwt(n))//'/s')
          scale_tij(k,n)=10.**(-ijtc_power(n))/dtsrc
        end if
      endif
C**** Tracers in river runoff (two versions - for inflow and outflow)
      k = k+1
      tij_rvr = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_rvr'
        lname_tij(k,n) = trim(TRNAME(n))//' in River Inflow'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
      k = k+1
      tij_rvro = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_rvro'
        lname_tij(k,n) = trim(TRNAME(n))//' in River Outflow'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracers in iceberg runoff 
      k = k+1
      tij_icb = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_icb'
        lname_tij(k,n) = trim(TRNAME(n))//' in Iceberg Inflow'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracers in sea ice
      k = k+1
      tij_seaice = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_ice'
        lname_tij(k,n) = trim(TRNAME(n))//' in Sea Ice'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
          denom_tij(k,n)=n_Water
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,cmrwt(n))
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
          dname_tij(k,n) = 'oicefr'
c        denom_tij(k,n)=n_Water ! if kg/kg units for non-water-isotopes
        end if
C**** Tracers conc. in ground component (ie. water or ice surfaces)
      k = k+1
      tij_grnd = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_at_Grnd'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' at Ground'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
          denom_tij(k,n)=n_Water
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)/REAL(NIsurf,KIND=8)
        end if
C**** Tracers conc. in lakes (layer 1)
      k = k+1
      tij_lk1 = k
        sname_tij(k,n) = trim(TRNAME(n))//'_Lake1'
        lname_tij(k,n) = trim(TRNAME(n))//' Lakes layer 1'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracers conc. in lakes (layer 2)
      k = k+1
      tij_lk2 = k
        sname_tij(k,n) = trim(TRNAME(n))//'_Lake2'
        lname_tij(k,n) = trim(TRNAME(n))//' Lakes layer 2'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracers conc. in soil water
      k = k+1
      tij_soil = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_Soil'
        lname_tij(k,n) = trim(TRNAME(n))//' Soil Water'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracers conc. in land snow water
      k = k+1
      tij_snow = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_Snow'
        lname_tij(k,n) = trim(TRNAME(n))//' Land Snow Water'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracer ice-ocean flux
      k = k+1
      tij_icocflx = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_ic_oc_flx'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Ice-Ocean Flux'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          denom_tij(k,n)=n_Water
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)-5,'kg/m^2/s')
          scale_tij(k,n)=10.**(-ijtc_power(n)+5)/DTsrc
        end if
C**** Tracers integrated E-W atmospheric flux
      k = k+1
      tij_uflx = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_uflx'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' E-W Atmos Flux'
        units_tij(k,n)=unit_string(ijtc_power(n)+10,'kg/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)-10)/DTsrc
C**** Tracers integrated N-S atmospheric flux
      k = k+1
      tij_vflx = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_vflx'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' N-S Atmos Flux'
        units_tij(k,n)=unit_string(ijtc_power(n)+10,'kg/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)-10)/DTsrc
C**** Tracers integrated E-W sea ice flux
      k = k+1
      tij_tusi = k
        sname_tij(k,n) = trim(TRNAME(n))//'_tusi'
        lname_tij(k,n) = trim(TRNAME(n))//' E-W Ice Flux'
        units_tij(k,n) = unit_string(ntrocn(n),'kg/s')
        scale_tij(k,n) = (10.**(-ntrocn(n)))/DTsrc
C**** Tracers integrated N-S sea ice flux
      k = k+1
      tij_tvsi = k
        sname_tij(k,n) = trim(TRNAME(n))//'_tvsi'
        lname_tij(k,n) = trim(TRNAME(n))//' N-S Ice Flux'
        units_tij(k,n) = unit_string(ntrocn(n),'kg/s')
        scale_tij(k,n) = (10.**(-ntrocn(n)))/DTsrc
#endif
#ifdef TRACERS_DRYDEP
C**** Tracers dry deposition flux.
      k = k+1
      tij_drydep = k
      if (dodrydep(n)) then
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_dry_dep'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Dry Deposition'
        units_tij(k,n)=unit_string(ijtc_power(n)-5,'kg/m^2/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)+5)/DTsrc
      end if
      k = k+1
      tij_gsdep = k
      if (trradius(n).gt.0) then
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_gs_dep'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Gravitational Settling'
        units_tij(k,n)=unit_string(ijtc_power(n)-5,'kg/m^2/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)+5)/DTsrc
      end if
#endif

      if (k .gt. ktaij) then
        if (AM_I_ROOT()) write (6,*)
     &   'tij_defs: Increase ktaij=',ktaij,' to at least ',k
        call stop_model('ktaij too small',255)
      end if

      end do

c
c Collect denominator short names for later use
c
      do n=1,ntm
      do k=1,ktaij
        nd = denom_tij(k,n)
        if(nd.gt.0) dname_tij(k,n) = sname_tij(k,nd)
      enddo
      enddo

      RETURN
      END SUBROUTINE set_generic_tracer_diags


      SUBROUTINE apply_tracer_2Dsource(dtstep)
!@sum apply_tracer_2Dsource adds surface sources to tracers
!@auth Jean Lerner/Gavin Schmidt
      USE MODEL_COM, only : jm
      USE GEOM, only : imaxj
      USE QUSDEF, only : mz,mzz
      USE TRACER_COM, only : ntm,trm,trmom,ntsurfsrc,ntisurfsrc,trname
      USE FLUXES, only : trsource,trflux1,trsrfflx
      USE TRDIAG_COM, only : taijs=>taijs_loc
      USE TRDIAG_COM, only : ijts_source,jls_source,itcon_surf
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, GET
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: dtstep
      INTEGER n,ns,naij,najl,j,i
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO
     *     ,grid%J_STRT_HALO:grid%J_STOP_HALO) :: dtracer
      REAL*8 ftr1

      INTEGER :: J_0, J_1, I_0, I_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      
C**** This is tracer independent coding designed to work for all
C**** surface sources.
C**** Note that tracer flux is added to first layer either implicitly
C**** in ATURB or explicitly in 'apply_fluxes_to_atm' call in SURFACE.

      do n=1,ntm
        trflux1(:,:,n) = 0.
C**** Non-interactive sources
        do ns=1,ntsurfsrc(n)
C**** diagnostics
          naij = ijts_source(ns,n)
          IF (naij > 0) THEN
          taijs(:,:,naij) = taijs(:,:,naij) + trsource(:,:,ns,n)*dtstep
          ENDIF
          najl = jls_source(ns,n)
          IF (najl > 0) THEN
            DO J=J_0,J_1
              DO I=I_0,imaxj(j)
                call inc_tajls(i,j,1,najl,trsource(i,j,ns,n)*dtstep)
              END DO
            END  DO
          END IF
          DO J=J_0,J_1
            do i=i_0,imaxj(j)
              dtracer(i,j)=trsource(i,j,ns,n)*dtstep
            end do
          end do
          if (itcon_surf(ns,n).gt.0)
     *         call DIAGTCB(dtracer,itcon_surf(ns,n),n)
C**** trflux1 is total flux into first layer
          trflux1(:,:,n) = trflux1(:,:,n)+trsource(:,:,ns,n)
        end do

       do j=j_0,j_1
C**** modify vertical moments (only from non-interactive sources)
c this is disabled until vertical moments are also modified
c during the vertical transport between layer 1 and the others
c        trmom( mz,:,j,1,n) = trmom( mz,:,j,1,n)-1.5*trflux1(:,j,n)
c     *       *dtstep
c        trmom(mzz,:,j,1,n) = trmom(mzz,:,j,1,n)+0.5*trflux1(:,j,n)
c     *       *dtstep

C**** Accumulate interactive sources as well
        trflux1(:,j,n) = trflux1(:,j,n)+trsrfflx(:,j,n)
       end do

C**** Technically speaking the vertical moments should be modified here
C**** as well. But for consistency with water vapour we only modify
C**** moments for dew.
        do j=J_0,J_1
          do i=i_0,imaxj(j)
            if (trsrfflx(i,j,n).lt.0 .and. trm(i,j,1,n).gt.0) then
              ftr1=min(1d0,-trsrfflx(i,j,n)*dtstep/trm(i,j,1,n))
              trmom(:,i,j,1,n)=trmom(:,i,j,1,n)*(1.-ftr1)
            end if
          end do
        end do
      end do
C****
      RETURN
      END SUBROUTINE apply_tracer_2Dsource

      MODULE apply3d
!@sum apply3d is used simply so that I can get optional arguments
!@+   to work. If anyone can some up with something neater, let me know.
      CONTAINS
      SUBROUTINE apply_tracer_3Dsource( ns , n , momlog )
!@sum apply_tracer_3Dsource adds 3D sources to tracers
!@auth Jean Lerner/Gavin Schmidt
      USE CONSTANT, only : teeny
      USE MODEL_COM, only : jm,im,lm,dtsrc
      USE GEOM, only : imaxj,byaxyp,lat2d_dg,lon2d_dg
      USE QUSDEF, only: nmom
      USE TRACER_COM, only : ntm,trm,trmom,trname,alter_sources,
     * num_regions,reg_S,reg_N,reg_E,reg_W,ef_FACT3d,num_tr_sectors
     * ,tr_sect_index3D,num_tr_sectors3d
      USE FLUXES, only : tr3Dsource
      USE TRDIAG_COM, only : jls_3Dsource,itcon_3Dsrc
     *     ,ijts_3Dsource,taijs=>taijs_loc
      USE DOMAIN_DECOMP_ATM, only : GRID, GET, am_i_root
      IMPLICIT NONE
!@var MOM true (default) if moments are to be modified
      logical, optional, intent(in) :: momlog
      integer, intent(in) :: n,ns
      real*8 fred(grid%i_strt:grid%i_stop),
     &     dtrm(grid%i_strt_halo:grid%i_stop_halo,
     &          grid%j_strt_halo:grid%j_stop_halo,lm),
     &     dtrml(lm),eps

      logical :: domom
      integer najl,i,j,l,naij,kreg,nsect,nn
      INTEGER :: J_0, J_1, I_0, I_1
      integer :: mask(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop)
      real*8 :: coef(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop)

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Ensure that this is a valid tracer and source
      if (n.eq.0 .or. ns.eq.0) then
         return
      end if
C**** parse options
      domom=.true.
      if( present(momlog) ) then
        domom=momlog
      end if
C**** This is tracer independent coding designed to work for all
C**** 3D sources.
C**** Modify tracer amount, moments, and diagnostics
      najl = jls_3Dsource(ns,n)
      naij = ijts_3Dsource(ns,n)

C**** apply tracer source alterations if requested in rundeck:
      if (alter_sources) then
        do kreg = 1, num_regions
          do j = j_0, j_1
            do i = i_0, imaxj(j)
              if (lat2d_dg(i,j) >= REG_S(kreg) .and. 
     &            lat2d_dg(i,j) <= REG_N(kreg) .and.
     &            lon2d_dg(i,j) >= REG_W(kreg) .and.
     &            lon2d_dg(i,j) < REG_E(kreg) ) then
                 mask(i,j) = 1
              else
                 mask(i,j) = 0
              end if
            end do
          end do

          do nsect = 1, num_tr_sectors3d(n,ns)
            nn = tr_sect_index3D(n, ns, nsect)
            if (ef_FACT3d(nn,kreg) > -1.e20) then
              do j = j_0, j_1
                 do i = i_0, imaxj(j)
                    coef(i,j) = 
     &                   mask(i,j) * ef_Fact3d(nn,kreg) + 
     &                   (1-mask(i,j))
                 end do
              end do
              do l = 1, lm
                do j = j_0, j_1
                  do i = i_0, imaxj(j)
                    tr3Dsource(i,j,l,ns,n) = tr3Dsource(i,j,l,ns,n) *
     &                    coef(i,j)
        
                  end do
                end do
              end do
            end if
          end do
        end do
      end if

      eps = tiny(trm(i_0,j_0,1,n))
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE (L,I,J,fred)
!$OMP&  SHARED(trm, dtrm, tr3Dsource, dtsrc, domom, naij, taijs, 
!$OMP&         imaxj, eps, j_0, j_1, i_0, ns, n, trmom)
      do l=1,lm
      do j=j_0,j_1
        do i=i_0,imaxj(j)
          dtrm(i,j,l) = tr3Dsource(i,j,l,ns,n)*dtsrc
C**** calculate fractional loss and update tracer mass
          fred(i) = max(0.,1.+min(0.,dtrm(i,j,l))/(trm(i,j,l,n)+eps))
          trm(i,j,l,n) = trm(i,j,l,n)+dtrm(i,j,l)
          if(fred(i).le.1d-16) trm(i,j,l,n) = 0.
        end do
        if(domom .and. any(fred.lt.1.)) then
          do i=i_0,imaxj(j)
            trmom(:,i,j,l,n) = trmom(:,i,j,l,n)*fred(i)
          enddo
        endif
      end do
      if (naij.gt.0) then
        do j=j_0,j_1
          do i=i_0,imaxj(j)
            taijs(i,j,naij) = taijs(i,j,naij) + dtrm(i,j,l)
          end do
        end do
      end if
      end do ! l
!$OMP END PARALLEL DO

      if(jls_3Dsource(ns,n) > 0) then
        do j=j_0,j_1
          do i=i_0,imaxj(j)
            dtrml(:) = dtrm(i,j,:)
            call inc_tajls_column(i,j,1,lm,lm,najl,dtrml)
          enddo
        enddo
      endif

      if (itcon_3Dsrc(ns,n).gt.0)
     *  call DIAGTCA(itcon_3Dsrc(ns,n),n)

C****
      RETURN
      END SUBROUTINE apply_tracer_3Dsource

      end module apply3d

      SUBROUTINE TDECAY
!@sum TDECAY decays radioactive tracers every source time step
!@auth Gavin Schmidt/Jean Lerner
      USE MODEL_COM, only : im,jm,lm,itime,dtsrc
      USE FLUXES, only : tr3Dsource
      USE GEOM, only : imaxj
      USE TRACER_COM, only : ntm,trm,trmom,trdecay,itime_tr0,n_Pb210,
     &     trname, n_Rn222
#ifdef TRACERS_WATER
     *     ,trwm
      USE SEAICE_COM, only : trsi
      USE LAKES_COM, only : trlake
      USE LANDICE_COM, only : trlndi,trsnowli
      USE GHY_COM, only : tr_w_ij,tr_wsn_ij
#endif
      USE TRDIAG_COM, only : jls_decay,itcon_decay
      USE DOMAIN_DECOMP_ATM, only : GRID, GET
      IMPLICIT NONE
      real*8, save, dimension(ntm) :: expdec = 1.
      real*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: told
      logical, save :: ifirst=.true.
      integer n,najl,j,l,i
      integer :: J_0, J_1, I_0, I_1

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      if (ifirst) then
        do n=1,ntm
          if (trdecay(n).gt.0.0) expdec(n)=exp(-trdecay(n)*dtsrc)
        end do
        ifirst = .false.
      end if

      do n=1,ntm
        if (trdecay(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** Atmospheric decay
          told(:,:,:)=trm(:,:,:,n)

#ifdef TRACERS_WATER
     *               +trwm(:,:,:,n)
          trwm(:,:,:,n)   = expdec(n)*trwm(:,:,:,n)
#endif
          if (trname(n) .eq. "Rn222" .and. n_Pb210.gt.0) then
            tr3Dsource(:,:,:,1,n_Pb210)= trm(:,:,:,n)*(1-expdec(n))*210.
     *           /222./dtsrc
          end if

          trm(:,:,:,n)    = expdec(n)*trm(:,:,:,n)
          trmom(:,:,:,:,n)= expdec(n)*trmom(:,:,:,:,n)

#ifdef TRACERS_WATER
C**** Note that ocean tracers are dealt with by separate ocean code.
C**** Decay sea ice tracers
          trsi(n,:,:,:)   = expdec(n)*trsi(n,:,:,:)
C**** ...lake tracers
          trlake(n,:,:,:) = expdec(n)*trlake(n,:,:,:)
C**** ...land surface tracers
          tr_w_ij(n,:,:,:,:) = expdec(n)*tr_w_ij(n,:,:,:,:)
          tr_wsn_ij(n,:,:,:,:)= expdec(n)*tr_wsn_ij(n,:,:,:,:)
          trsnowli(n,:,:) = expdec(n)*trsnowli(n,:,:)
          trlndi(n,:,:)   = expdec(n)*trlndi(n,:,:)
#endif
C**** atmospheric diagnostics
          najl = jls_decay(n)
          do l=1,lm
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              call inc_tajls(i,j,l,najl,trm(i,j,l,n)
#ifdef TRACERS_WATER
     *             +trwm(i,j,l,n)
#endif
     *             -told(i,j,l))
            enddo
          enddo
          enddo

          call DIAGTCA(itcon_decay(n),n)
        end if
      end do
C****
      return
      end subroutine tdecay


      SUBROUTINE TRGRAV
!@sum TRGRAV gravitationally settles particular tracers
!@auth Gavin Schmidt/Reha Cakmur
      USE CONSTANT, only : grav,deltx,lhe,rgas,visc_air
      USE MODEL_COM, only : im,jm,lm,itime,dtsrc,zatmo,t,q
      USE GEOM, only : imaxj,byaxyp
      USE SOMTQ_COM, only : mz,mzz,mzx,myz,zmoms
      USE DYNAMICS, only : gz,pmid,pk
      USE TRACER_COM, only : ntm,trm,trmom,itime_tr0,trradius
     *     ,trname,trpdens
#ifdef TRACERS_AMP
     *     ,AMP_MODES_MAP,AMP_NUMB_MAP,ntmAMP
      USE AMP_AEROSOL, only : DIAM, AMP_dens
      USE AERO_SETUP,  only : CONV_DPAM_TO_DGN
#endif
      USE TRDIAG_COM, only : jls_grav
      USE DOMAIN_DECOMP_ATM, only : GRID, GET
      IMPLICIT NONE
      real*8 :: stokevdt,press,fgrfluxd,qsat,vgs,tr_radius,tr_dens,temp
      real*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: told,airden,visc,rh
     *     ,gbygz
      real*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO) :: fluxd, fluxu
      integer n,najl,i,j,l
      integer :: J_0, J_1, I_0, I_1
      logical :: hydrate

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Calculate some tracer independent arrays      
C**** air density + relative humidity (wrt water) + air viscosity
      do l=1,lm
        do j=J_0,J_1
          do i=I_0,imaxj(j)
            press=pmid(l,i,j)
            temp=pk(l,i,j)*t(i,j,l)
            airden(i,j,l)=100.d0*press/(rgas*temp*(1.+q(i,j,l)*deltx))
            rh(i,j,l)=q(i,j,l)/qsat(temp,lhe,press)
            visc(i,j,l)=visc_air(temp)
            if (l.eq.1) then
              gbygz(i,j,l)=0.
            else
              gbygz(i,j,l)=grav/(gz(i,j,l)-gz(i,j,l-1))
            end if
          end do
        end do
      end do

C**** Gravitational settling
!$OMP  PARALLEL DO DEFAULT (NONE)
!$OMP&   SHARED(J_0,J_1,I_0,IMAXJ,trm,trpdens,trradius,dtsrc,
!$OMP&          airden,rh,visc,gbygz,trmom,itime,itime_tr0,jls_grav )
!$OMP&   PRIVATE (n,l,i,j,stokevdt,fgrfluxd,tr_dens,tr_radius,
!$OMP&      hydrate, fluxd, fluxu, told, najl)
      do n=1,ntm
        if (trradius(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** need to hydrate the sea salt before determining settling
          hydrate = (trname(n).eq.'seasalt1'.or.trname(n).eq.'seasalt2')

          fluxd=0.
          do l=lm,1,-1          ! loop down
            do j=J_0,J_1
              do i=I_0,imaxj(j)

C*** save original tracer mass
                told(i,j,l)=trm(i,j,l,n)
C**** set incoming flux from previous level
                fluxu(i,j)=fluxd(i,j)

C**** set particle properties
                tr_dens = trpdens(n)
                tr_radius = trradius(n)

#ifdef TRACERS_AMP
       if (n.le.ntmAMP) then
        if(AMP_MODES_MAP(n).gt.0) then
        if(DIAM(i,j,l,AMP_MODES_MAP(n)).gt.0.) then
        if(AMP_NUMB_MAP(n).eq. 0) then    ! Mass
        tr_radius=0.5*DIAM(i,j,l,AMP_MODES_MAP(n))
        else                              ! Number
        tr_radius=0.5*DIAM(i,j,l,AMP_MODES_MAP(n))
     +            *CONV_DPAM_TO_DGN(AMP_MODES_MAP(n))
        endif

        call AMPtrdens(i,j,l,n)
        tr_dens =AMP_dens(i,j,l,AMP_MODES_MAP(n))
        endif   
        endif   
       endif 
#endif  

C**** calculate stoke's velocity (including possible hydration effects
C**** and slip correction factor)
                stokevdt=dtsrc*vgs(airden(i,j,l),rh(i,j,l),tr_radius
     *               ,tr_dens,visc(i,j,l),hydrate)

C**** Calculate height differences using geopotential
                fgrfluxd=stokevdt*gbygz(i,j,l) 
                fluxd(i,j) = trm(i,j,l,n)*fgrfluxd ! total flux down
                trm(i,j,l,n) = trm(i,j,l,n)*(1.-fgrfluxd)+fluxu(i,j)
                if (1.-fgrfluxd.le.1d-16) trm(i,j,l,n) = fluxu(i,j)
                trmom(zmoms,i,j,l,n)=trmom(zmoms,i,j,l,n)*(1.-fgrfluxd)
              end do
            end do
          end do
          najl = jls_grav(n)
          IF (najl > 0) THEN
            do l=1,lm
              do j=J_0,J_1
                do i=I_0,imaxj(j)
                  call inc_tajls(i,j,l,najl,trm(i,j,l,n)-told(i,j,l))
                enddo
              enddo
            enddo
          END IF
        end if
      end do
!$OMP END PARALLEL DO

C****

      return
      end subroutine trgrav

      REAL*8 FUNCTION vgs(airden,rh1,tr_radius,tr_dens,visc,hydrate)
!@sum vgs returns settling velocity for tracers (m/s)
!@auth Gavin Schmidt/Reha Cakmur
      USE CONSTANT, only : by3,pi,gasc,avog,rt2,deltx
     *     ,mair,grav
      IMPLICIT NONE
      real*8, intent(in) ::  airden,rh1,tr_radius,tr_dens,visc
      logical, intent(in) :: hydrate
      real*8  wmf,frpath
      real*8, parameter :: dair=3.65d-10 !m diameter of air molecule
C**** coefficients for slip factor calculation
      real*8, parameter :: s1=1.247d0, s2=0.4d0, s3=1.1d0
C**** coefficients for hydration radius if required 
      real*8, parameter :: c1=0.7674d0, c2=3.079d0, c3=2.573d-11,
     *     c4=-1.424d0
      real*8 dens, rad, rh

C**** Determine if hydration of aerosols is required
      if (hydrate) then
        rh=max(0.01d0,min(rh1,0.99d0))
C**** hydrated radius (Gerber, 1988)
        rad=(c1*tr_radius**(c2)/(c3*tr_radius**(c4)-log10(rh))
     *       + tr_radius**3)**by3
C**** hydrated density
        dens=((rad**3 - tr_radius**3)*1000.d0
     *       + tr_radius**3*tr_dens)/rad**3
      else                      ! take dry values
        dens = tr_dens
        rad = tr_radius 
      end if

C**** calculate stoke's velocity
      vgs=2.*grav*dens*rad**2/(9.*visc)

C**** slip correction factor
c wmf is the additional velocity if the particle size is small compared
c   to the mean free path of the air; important in the stratosphere
      frpath=1d-3*mair/(pi*rt2*avog*airden*(dair)**2.)
      wmf=frpath/tr_radius*(s1+s2*exp(-s3*tr_radius/frpath))
      vgs=(1.d0+wmf)*vgs
C****
      return
      end function vgs

      SUBROUTINE read_monthly_sources(iu,jdlast,tlca,tlcb,data,
     *  frac,imon)
!@sum Read in monthly sources and interpolate to current day
!@+   Calling routine must have the lines:
!@+      real*8 tlca(im,jm,nm),tlcb(im,jm,nm)
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      data jdlast /0/
!@+      save jdlast,tlca,tlcb,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
!@auth Jean Lerner and others
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP_ATM, only : GRID, GET, AM_I_ROOT
      USE DOMAIN_DECOMP_ATM, only : READT_PARALLEL, REWIND_PARALLEL
      USE MODEL_COM, only: jday,im,jm,idofm=>JDmidOfM
      implicit none
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     tlca,tlcb,data
      real*8 :: frac
      integer imon,iu,jdlast

      integer :: J_0, J_1, I_0, I_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      CALL GET(grid, I_STRT=I_0, I_STOP=I_1)

      if (jdlast.EQ.0) then ! NEED TO READ IN FIRST MONTH OF DATA
        imon=1          ! imon=January
        if (jday.le.16)  then ! JDAY in Jan 1-15, first month is Dec
          CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),tlca,12)
          CALL REWIND_PARALLEL( iu )
        else            ! JDAY is in Jan 16 to Dec 16, get first month
  120     imon=imon+1
          if (jday.gt.idofm(imon) .AND. imon.le.12) go to 120
          CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),tlca,imon-1)
          if (imon.eq.13)  CALL REWIND_PARALLEL( iu )
        end if
      else              ! Do we need to read in second month?
        if (jday.ne.jdlast+1) then ! Check that data is read in daily
          if (jday.ne.1 .OR. jdlast.ne.365) then
            if (AM_I_ROOT()) write(6,*)
     *      'Incorrect values in Tracer Source:JDAY,JDLAST=',JDAY,JDLAST
            call stop_model('stopped in TRACERS.f',255)
          end if
          imon=imon-12  ! New year
          go to 130
        end if
        if (jday.le.idofm(imon)) go to 130
        imon=imon+1     ! read in new month of data
        tlca(I_0:I_1,J_0:J_1) = tlcb(I_0:I_1,J_0:J_1)
        if (imon.eq.13) CALL REWIND_PARALLEL( iu  )
      end if
      CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),tlcb,1)
  130 continue
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
      data(I_0:I_1,J_0:J_1) = tlca(I_0:I_1,J_0:J_1)*frac + 
     & tlcb(I_0:I_1,J_0:J_1)*(1.-frac)
      return
      end subroutine read_monthly_sources
#endif

      subroutine checktr(subr)
!@sum  CHECKTR Checks whether atmos tracer variables are reasonable
!@auth Gavin Schmidt
!@ver  1.0
#ifdef TRACERS_ON
      USE CONSTANT, only : teeny
      USE MODEL_COM, only : ls1,im,jm,lm,q,wm
      USE GEOM, only : axyp,imaxj
      USE SOMTQ_COM, only : qmom
      USE DYNAMICS, only : am
      USE FLUXES, only : gtracer
      USE TRACER_COM
      USE DOMAIN_DECOMP_ATM, ONLY: GRID, GET, AM_I_ROOT
      IMPLICIT NONE
      LOGICAL QCHECKT
      INTEGER I,J,L,N,m, imax,jmax,lmax
      REAL*8 relerr, errmax,errsc,tmax,amax,qmax,wmax,twmax,qmomax(nmom)
     *     ,tmomax(nmom)
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
      INTEGER :: J_0, J_1, nj, I_0,I_1
      CALL GET(GRID, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      nj = J_1 - J_0 + 1

      CALL CHECK4(gtracer(1,1,1,J_0),NTM,4,IM,nJ,SUBR,'GTRACE')
      do n=1,ntm
        CALL CHECK4(trmom(:,:,J_0:J_1,:,n),NMOM,IM,nJ,LM,SUBR,
     *       'X'//trname(n))
        CALL CHECK3(trm(:,J_0:J_1,:,n),IM,nJ,LM,SUBR,trname(n))
#ifdef TRACERS_WATER
        CALL CHECK3(trwm(:,J_0:J_1,:,n),IM,nJ,LM,SUBR,'WM'//trname(n))
#endif

C**** check for negative tracer amounts (if t_qlimit is set)
        if (t_qlimit(n)) then
          QCHECKT=.false.
          do l=1,lm
          do j=j_0,j_1
          do i=i_0,imaxj(j)
            if (trm(i,j,l,n).lt.0) then
              if (AM_I_ROOT())
     *         write(6,*) "Negative mass for ",trname(n),i,j,l,trm(i,j,l
     *             ,n)," after ",SUBR,"."
              QCHECKT=.true.
            end if
          end do
          end do
          end do
          if (QCHECKT)
     &         call stop_model("CHECKTR: Negative tracer amount",255)
        end if

C**** check whether air mass is conserved

        if (trname(n).eq.'Air') then
          errmax = 0. ; lmax=1 ; imax=I_0 ; jmax=J_0; tmax=0. ; amax=0.
          do l=1,lm
          do j=j_0,j_1
          do i=i_0,imaxj(j)
            relerr=abs(trm(i,j,l,n)-am(l,i,j)*axyp(i,j))/
     *           (am(l,i,j)*axyp(i,j))
            if (relerr.gt.errmax) then
              lmax=l ; imax=i ; jmax=j ; errmax=relerr
              tmax=trm(i,j,l,n) ; amax=am(l,i,j)*axyp(i,j)
            end if
          end do
          end do
          end do
          print*,"Relative error in air mass after ",trim(subr),":",imax
     *         ,jmax,lmax,errmax,tmax,amax
        end if

#ifdef TRACERS_WATER
        if (trname(n).eq.'Water') then
          errmax = 0. ; lmax=1 ; imax=I_0 ; jmax=J_0
          tmax=0. ; twmax=0. ; qmax=0. ; wmax=0.
          tmomax = 0. ; qmomax = 0. 
          do l=1,lm
          do j=j_0,j_1
          do i=i_0,imaxj(j)
            errsc=(q(i,j,l)+sum(abs(qmom(:,i,j,l))))*am(l,i,j)*axyp(i,j)
            if (errsc.eq.0.) errsc=1.
            relerr=abs(trm(i,j,l,n)-q(i,j,l)*am(l,i,j)*axyp(i,j))/errsc
            if (wm(i,j,l).gt.0 .and. trwm(i,j,l,n).gt.1.) relerr
     *           =max(relerr,(trwm(i,j,l,n)-wm(i,j,l)*am(l,i,j)*axyp(i,j
     *           ))/(wm(i,j,l)*am(l,i,j)*axyp(i,j)))
            if ((wm(i,j,l).eq.0 .and.trwm(i,j,l,n).gt.1) .or. (wm(i,j,l)
     *           .gt.teeny .and.trwm(i,j,l,n).eq.0))
     *           print*,"Liquid water mismatch: ",subr,i,j,l,trwm(i,j,l
     *           ,n),wm(i,j,l)*am(l,i,j)*axyp(i,j)
            do m=1,nmom
              relerr=max(relerr,(trmom(m,i,j,l,n)-qmom(m,i,j,l)*am(l,i,j
     *             )*axyp(i,j))/errsc)
            end do
            if (relerr.gt.errmax) then
              lmax=l ; imax=i ; jmax=j ; errmax=relerr
              tmax=trm(i,j,l,n) ; qmax=q(i,j,l)*am(l,i,j)*axyp(i,j)
              twmax=trwm(i,j,l,n) ; wmax=wm(i,j,l)*am(l,i,j)*axyp(i,j)
              tmomax(:)=trmom(:,i,j,l,n)
              qmomax(:)=qmom(:,i,j,l)*am(l,i,j)*axyp(i,j)
            end if
          end do
          end do
          end do
          print*,"Relative error in water mass after ",trim(subr),":"
     *         ,imax,jmax,lmax,errmax,tmax,qmax,twmax,wmax,(tmomax(m)
     *         ,qmomax(m),m=1,nmom) 
        end if
#endif
      end do
#endif
      return
      end subroutine checktr


      SUBROUTINE io_tracer(kunit,iaction,ioerr)
!@sum  io_tracer reads and writes tracer variables to file
!@auth Jean Lerner
!@ver  1.0
#ifdef TRACERS_ON
      USE MODEL_COM, only: ioread,iowrite,irsfic,irsficno,irerun,lhead
     &     ,coupled_chem
      USE DOMAIN_DECOMP_1D, only : grid,AM_I_ROOT,PACK_DATA,UNPACK_DATA
     &     ,PACK_BLOCK, UNPACK_BLOCK, PACK_COLUMN
     &     ,UNPACK_COLUMN, esmf_bcast, get
      USE TRACER_COM
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,
     &     yROR,yXO2,yAldehyde,yXO2N,yRXPAR,ss,JPPJ,ydms,yso2,sulfate
     &     ,acetone, sOx_acc,sNOx_acc,sCO_acc,l1Ox_acc,l1NO2_acc
#ifdef SHINDELL_STRAT_CHEM
     &     ,SF3,SF2,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2
#endif
#ifdef INTERACTIVE_WETLANDS_CH4 
      use TRACER_SOURCES, only: day_ncep,DRA_ch4,sum_ncep,PRS_ch4,
     &     HRA_ch4,iday_ncep,i0_ncep,iHch4,iDch4,i0ch4,first_ncep,
     &     first_mod,max_days,nra_ncep,nra_ch4,maxHR_ch4,avg_model,
     &     avg_ncep
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE fluxes,ONLY : nstype,pprec,pevap
      USE tracers_dust,ONLY : hbaij,ricntd
      use trdust_drv, only: io_trDust
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      USE AEROSOL_SOURCES, only : snosiz
#endif
      USE param, only : sync_param
      use trdiag_com, only: trcSurfMixR_acc,trcSurfByVol_acc
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
     &     ,sPM2p5_acc,sPM10_acc,l1PM2p5_acc,l1PM10_acc
     &     ,csPM2p5_acc,csPM10_acc
#endif

      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records

      REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: TRMOM_GLOB
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: TRM_GLOB
#ifdef TRACERS_WATER
     &     ,TRWM_GLOB
#endif
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Aijl_glob
  
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: snosiz_glob
#endif
#ifdef TRACERS_SPECIAL_Shindell
      REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: ss_glob
      real(kind=8),allocatable,dimension(:,:) :: sOx_acc_glob,
     & sNOx_acc_glob, sCO_acc_glob, l1Ox_acc_glob, l1NO2_acc_glob
#ifdef INTERACTIVE_WETLANDS_CH4 
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE ::
     &     rHch4,rDch4,r0ch4,rfirst_mod
      REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE ::
     &     day_ncep_glob,DRA_ch4_glob,HRA_ch4_glob
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Rijch4_glob
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Rijncep_glob
#endif
#endif     
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      REAL*8,DIMENSION(Im,Jm) :: pprec_glob,ricntd_glob,hbaij_glob
      REAL*8,DIMENSION(Im,Jm,Nstype) :: pevap_glob
#endif
      real(kind=8),allocatable,dimension(:,:,:) :: trcSurfMixR_acc_glob
     &     ,trcSurfByVol_acc_glob
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
      real(kind=8),allocatable,dimension(:,:) :: sPM2p5_acc_glob
     &     ,sPM10_acc_glob,l1PM2p5_acc_glob,l1PM10_acc_glob
     &     ,csPM2p5_acc_glob,csPM10_acc_glob
#endif
      INTEGER :: ITM,ITM1,ITM2
#ifdef TRACERS_WATER
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRACERW01"
#else
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRACER01"
#endif

      INTEGER :: J_0, J_1, J_1H, J_0H, I_0, I_1
      INTEGER :: img, jmg

      CALL GET(grid,I_STRT=I_0,I_STOP=I_1, 
     &              J_STRT=J_0,J_STOP=J_1, 
     &         J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

      if(am_i_root()) then
         img = IM
         jmg = JM
      else
         img = 1
         jmg = 1
      end if

      ! update database parameters which are needed for reading
#ifdef TRACERS_ON
      call sync_param( "coupled_chem", coupled_chem )
#endif

      allocate(
     &     TRM_GLOB(img,jmg,LM)
     &     ,TRMOM_GLOB(NMOM,img,jmg,LM)
#ifdef TRACERS_WATER
     &    ,TRWM_GLOB(img,jmg,LM)
#endif
     &    ,Aijl_glob(img,jmg,LM)
     &     )

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      allocate( snosiz_glob(img,jmg) )
#endif

#ifdef TRACERS_SPECIAL_Shindell
      allocate(
     &    ss_glob(JPPJ,LM,img,jmg)
     &    ,sOx_acc_glob(img,jmg)
     &    ,sNOx_acc_glob(img,jmg)
     &    ,sCO_acc_glob(img,jmg)
     &    ,l1Ox_acc_glob(img,jmg) 
     &    ,l1NO2_acc_glob(img,jmg))
#ifdef INTERACTIVE_WETLANDS_CH4
      allocate(
     &     day_ncep_glob(img,jmg,max_days,nra_ncep)
     &    ,DRA_ch4_glob(img,jmg,max_days,nra_ch4)
     &    ,HRA_ch4_glob(img,jmg,maxHR_ch4,nra_ch4)
     &    ,Rijncep_glob(img,jmg,nra_ncep)
     &    ,Rijch4_glob(img,jmg,nra_ch4) )
      allocate(
     &     rfirst_mod(IM,J_0H:J_1H,nra_ch4) 
     &    ,rHch4(IM,J_0H:J_1H,nra_ch4)
     &    ,rDch4(IM,J_0H:J_1H,nra_ch4)
     &    ,r0ch4(IM,J_0H:J_1H,nra_ch4) )
#endif
#endif

      allocate(trcSurfMixR_acc_glob(im,jm,ntm)
     &        ,trcSurfByVol_acc_glob(im,jm,ntm))
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
      allocate(sPM2p5_acc_glob(im,jm)
     &        ,sPM10_acc_glob(im,jm)
     &        ,l1PM2p5_acc_glob(im,jm)
     &        ,l1PM10_acc_glob(im,jm)
     &        ,csPM2p5_acc_glob(im,jm)
     &        ,csPM10_acc_glob(im,jm) )
#endif

      SELECT CASE (IACTION)

      CASE (:IOWRITE) ! output to end-of-month restart file
      
        DO ITM=1,NTM
          CALL PACK_DATA(grid, TRM(:,:,:,ITM), TRM_GLOB)
          CALL PACK_COLUMN(grid, TRmom(:,:,:,:,ITM),TRmom_GLOB)
#ifdef TRACERS_WATER
          CALL PACK_DATA(grid, TRWM(:,:,:,ITM), TRWM_GLOB)
#endif
          IF (AM_I_ROOT()) THEN
            write (MODULE_HEADER(lhead+1:80),'(a,a,i1,a,a)')
     &           trname(itm),
     &           ' R8 TRM(im,jm,lm),TRmom(',
     &           NMOM,
     &           ',im,jm,lm)'
#ifdef TRACERS_WATER
     &           ,',trwm(im,jm,lm)'
#endif
            WRITE (kunit,err=10) MODULE_HEADER,TRM_glob,TRmom_glob
#ifdef TRACERS_WATER
     *           ,TRWM_glob
#endif
          END IF                !only root processor writes
        END DO

c not yet        header='For tracer 3D emissions: daily_z(i,j,l)'
c not yet        call pack_data(grid,daily_z,aijl_glob)
c not yet        if(am_i_root()) write(kunit,err=10) header,aijl_glob

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
       CALL pack_data(grid,pprec,pprec_glob)
       CALL pack_data(grid,pevap,pevap_glob)
       CALL pack_data(grid,hbaij,hbaij_glob)
       CALL pack_data(grid,ricntd,ricntd_glob)
       header='For dust tracers: hbaij,ricntd,pprec,pevap'
       IF (am_i_root()) WRITE(kunit,ERR=10) header,hbaij_glob,
     &      ricntd_glob,pprec_glob,pevap_glob
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
       header='BC_ALB: snosiz(i,j)'
        call pack_data(grid,snosiz(:,:),snosiz_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,snosiz_glob
#endif

#ifdef TRACERS_SPECIAL_Shindell       
       header='TRACERS_SPECIAL_Shindell: ss(JPPJ,l,i,j)'
        call pack_block(grid,ss(:,:,:,:),ss_glob(:,:,:,:))
        if(am_i_root())write(kunit,err=10)header,ss_glob
       header='TRACERS_SPECIAL_Shindell: yNO3(i,j,l)'
        call pack_data(grid,yNO3,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: pHOx(i,j,l)'
        call pack_data(grid,pHOx,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: pNOx(i,j,l)'
        call pack_data(grid,pNOx,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: pOx(i,j,l)'
        call pack_data(grid,pOx,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: yCH3O2(i,j,l)'
        call pack_data(grid,yCH3O2,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: yC2O3(i,j,l)'
        call pack_data(grid,yC2O3,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: yROR(i,j,l)'
        call pack_data(grid,yROR,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: yXO2(i,j,l)'
        call pack_data(grid,yXO2,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: yXO2N(i,j,l)'
        call pack_data(grid,yXO2N,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: yAldehyde(i,j,l)'
        call pack_data(grid,yAldehyde,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: yRXPAR(i,j,l)'
        call pack_data(grid,yRXPAR,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: ydms(i,j,l)'
        call pack_data(grid,ydms,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: ySO2(i,j,l)'
        call pack_data(grid,ySO2,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: sulfate(i,j,l)'
        call pack_data(grid,sulfate,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='TRACERS_SPECIAL_Shindell: acetone(i,j,l)'
        call pack_data(grid,acetone,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       if(coupled_chem == 1)then
         header='TRACERS_SPECIAL_Shindell: oh_live(i,j,l)'
          call pack_data(grid,oh_live,Aijl_glob)
          if(am_i_root())write(kunit,err=10)header,Aijl_glob
         header='TRACERS_SPECIAL_Shindell: no3_live(i,j,l)'
          call pack_data(grid,no3_live,Aijl_glob)
          if(am_i_root())write(kunit,err=10)header,Aijl_glob
       endif
#ifdef SHINDELL_STRAT_CHEM
       header='SHINDELL_STRAT_CHEM: SF3(i,j,l)'
        call pack_data(grid,SF3,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='SHINDELL_STRAT_CHEM: SF2(i,j,l)'
        call pack_data(grid,SF2,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='SHINDELL_STRAT_CHEM: pClOx(i,j,l)'
        call pack_data(grid,pClOx,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='SHINDELL_STRAT_CHEM: pClx(i,j,l)'
        call pack_data(grid,pClx,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='SHINDELL_STRAT_CHEM: pOClOx(i,j,l)'
        call pack_data(grid,pOClOx,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='SHINDELL_STRAT_CHEM: pBrOx(i,j,l)'
        call pack_data(grid,pBrOx,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='SHINDELL_STRAT_CHEM: yCl2(i,j,l)'
        call pack_data(grid,yCl2,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
       header='SHINDELL_STRAT_CHEM: yCl2O2(i,j,l)'
        call pack_data(grid,yCl2O2,Aijl_glob)
        if(am_i_root())write(kunit,err=10)header,Aijl_glob
#endif
#ifdef INTERACTIVE_WETLANDS_CH4 
       header='INTERACTIVE_WETLANDS_CH4: day_ncep(i,j,days,#raN)'
        do itm=1,max_days
         do itm2=1,nra_ncep
           call pack_data
     &     (grid,day_ncep(:,:,itm,itm2),day_ncep_glob(:,:,itm,itm2))
         end do
        end do
        if(am_i_root())write(kunit,err=10)header,day_ncep_glob
       header='INTERACTIVE_WETLANDS_CH4: dra_ch4(i,j,days,#raC)'
        do itm=1,max_days
         do itm2=1,nra_ch4
           call pack_data
     &     (grid,dra_ch4(:,:,itm,itm2),dra_ch4_glob(:,:,itm,itm2))
         end do
        end do
        if(am_i_root())write(kunit,err=10)header,dra_ch4_glob
       header='INTERACTIVE_WETLANDS_CH4: sum_ncep(i,j,#raN)'
        call pack_data(grid,sum_ncep(:,:,:),Rijncep_glob(:,:,:))
        if(am_i_root())write(kunit,err=10)header,Rijncep_glob
       header='INTERACTIVE_WETLANDS_CH4: prs_ch4(i,j,#raC)'
        call pack_data(grid,prs_ch4(:,:,:),Rijch4_glob(:,:,:))
        if(am_i_root())write(kunit,err=10)header,Rijch4_glob
       header='INTERACTIVE_WETLANDS_CH4: HRA_ch4(i,j,mxHc,#raC)'
        do itm=1,maxHR_ch4
          do itm2=1,nra_ch4
            call pack_data
     &      (grid,HRA_ch4(:,:,itm,itm2),HRA_ch4_glob(:,:,itm,itm2))
          end do
        end do
        if(am_i_root())write(kunit,err=10)header,HRA_ch4_glob
       header='INTERACTIVE_WETLANDS_CH4: i0ch4(i,j,#raC) (real)'
        r0ch4(I_0:I_1,J_0:J_1,:)=REAL(i0ch4(I_0:I_1,J_0:J_1,:))
        call pack_data(grid,r0ch4(:,:,:),Rijch4_glob(:,:,:))
        if(am_i_root())write(kunit,err=10)header,Rijch4_glob
       header='INTERACTIVE_WETLANDS_CH4: iDch4(i,j,#raC) (real)'
        rDch4(I_0:I_1,J_0:J_1,:)=REAL(iDch4(I_0:I_1,J_0:J_1,:))
        call pack_data(grid,rDch4(:,:,:),Rijch4_glob(:,:,:))
        if(am_i_root())write(kunit,err=10)header,Rijch4_glob
       header='INTERACTIVE_WETLANDS_CH4: iHch4(i,j,#raC) (real)'
        rHch4(I_0:I_1,J_0:J_1,:)=REAL(iHch4(I_0:I_1,J_0:J_1,:))
        call pack_data(grid,rHch4(:,:,:),Rijch4_glob(:,:,:))
        if(am_i_root())write(kunit,err=10)header,Rijch4_glob
       header='INTERACTIVE_WETLANDS_CH4: first_mod(i,j,#raC) (real)'
        rfirst_mod(I_0:I_1,J_0:J_1,:)=REAL(first_mod(I_0:I_1,J_0:J_1,:))
        call pack_data(grid,rfirst_mod(:,:,:),Rijch4_glob(:,:,:))
        if(am_i_root())write(kunit,err=10)header,Rijch4_glob
       header='INTERACTIVE_WETLANDS_CH4: avg_model(i,j,#raC)'
        call pack_data(grid,avg_model(:,:,:),Rijch4_glob(:,:,:))
        if(am_i_root())write(kunit,err=10)header,Rijch4_glob
       header='INTERACTIVE_WETLANDS_CH4: avg_ncep(i,j,#raN)'
        call pack_data(grid,avg_ncep(:,:,:),Rijncep_glob(:,:,:))
        if(am_i_root())write(kunit,err=10)header,Rijncep_glob
       header='INTERACTIVE_WETLANDS_CH4: iday_ncep,i0_ncep,first_ncep'
        if(am_i_root())write(kunit,err=10)
     &  header,iday_ncep,i0_ncep,first_ncep
#endif /* INTERACTIVE_WETLANDS_CH4 */
       call pack_data(grid,sOx_acc,sOx_acc_glob)
       call pack_data(grid,sNOx_acc,sNOx_acc_glob)
       call pack_data(grid,sCO_acc,sCO_acc_glob)
       call pack_data(grid,l1Ox_acc,l1Ox_acc_glob)
       call pack_data(grid,l1NO2_acc,l1NO2_acc_glob)
#endif /* TRACERS_SPECIAL_Shindell */

       call pack_data(grid,trcSurfMixR_acc,trcSurfMixR_acc_glob)
       call pack_data(grid,trcSurfByVol_acc,trcSurfByVol_acc_glob)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
       call pack_data(grid,sPM2p5_acc,sPM2p5_acc_glob)
       call pack_data(grid,sPM10_acc,sPM10_acc_glob)
       call pack_data(grid,l1PM2p5_acc,l1PM2p5_acc_glob)
       call pack_data(grid,l1PM10_acc,l1PM10_acc_glob)
       call pack_data(grid,csPM2p5_acc,csPM2p5_acc_glob)
       call pack_data(grid,csPM10_acc,csPM10_acc_glob)
#endif
       header='accumulation arrays for subdd diagnostics for tracers'
       if (am_i_root()) write(kunit,err=10) header,trcSurfMixR_acc_glob
     &      ,trcSurfByVol_acc_glob
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
     &      ,sPM2p5_acc_glob,sPM10_acc_glob,l1PM2p5_acc_glob
     &      ,l1PM10_acc_glob,csPM2p5_acc_glob,csPM10_acc_glob
#endif
#ifdef TRACERS_SPECIAL_Shindell
     &      ,sOx_acc_glob,sNOx_acc_glob,sCO_acc_glob,l1Ox_acc_glob
     &      ,l1NO2_acc_glob
#endif

      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread,irerun,irsfic,irsficno) ! restarts
 
          DO ITM=1,NTM
            if ( AM_I_ROOT() ) then
              READ (kunit,err=10) HEADER,TRM_glob,TRmom_glob
#ifdef TRACERS_WATER
     *             ,TRWM_glob
#endif
              IF (HEADER(1:lhead).ne.MODULE_HEADER(1:lhead)) THEN
                PRINT*,"Discrepancy in module version ",HEADER,
     &               MODULE_HEADER
                GO TO 10
              ENDIF
            endif               ! AM_I_ROOT
C**** ESMF: Copy global data into the corresponding local (distributed) arrays.
            CALL UNPACK_DATA  (grid,   TRM_GLOB, TRM(:,:,:,  itm))
            CALL UNPACK_COLUMN(grid, TRMOM_GLOB, TRmom(:,:,:,:,itm))
#ifdef TRACERS_WATER
            CALL UNPACK_DATA  (grid,  TRWM_GLOB, TRWM(:,:,:,  itm))
#endif
          ENDDO

c not yet          if(am_i_root()) read(kunit,err=10) header,aijl_glob
c not yet          call unpack_data(grid,aijl_glob,daily_z)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
          IF (am_i_root()) READ(kunit,ERR=10) header,hbaij_glob,
     &         ricntd_glob,pprec_glob,pevap_glob
          CALL unpack_data(grid,hbaij_glob,hbaij)
          CALL unpack_data(grid,ricntd_glob,ricntd)
          CALL unpack_data(grid,pprec_glob,pprec)
          CALL unpack_data(grid,pevap_glob,pevap)
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
          if(am_i_root())read(kunit,err=10)header,snosiz_glob
          call unpack_data(grid,snosiz_glob(:,:),snosiz(:,:))
#endif

#ifdef TRACERS_SPECIAL_Shindell       
          if(am_i_root())read(kunit,err=10)header,ss_glob
          call unpack_block(grid,ss_glob(:,:,:,:),ss(:,:,:,:))
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,yNO3)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,pHOx)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,pNOx)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,pOx)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,yCH3O2)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,yC2O3)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,yROR)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,yXO2)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,yXO2N)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,yAldehyde)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,yRXPAR)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,ydms)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,ySO2)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,sulfate)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,acetone)
          if(coupled_chem == 1)then
            if(am_i_root())read(kunit,err=10)header,Aijl_glob
            call unpack_data(grid,Aijl_glob,oh_live)
            if(am_i_root())read(kunit,err=10)header,Aijl_glob
            call unpack_data(grid,Aijl_glob,no3_live)
          endif
#ifdef SHINDELL_STRAT_CHEM
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,SF3)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,SF2)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,pClOx)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,pClx)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,pOClOx)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,pBrOx)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,yCl2)
          if(am_i_root())read(kunit,err=10)header,Aijl_glob
          call unpack_data(grid,Aijl_glob,yCl2O2)
#endif
#ifdef INTERACTIVE_WETLANDS_CH4 
          if(am_i_root())read(kunit,err=10)header,day_ncep_glob
          do itm=1,max_days ;do itm2=1,nra_ncep
            call unpack_data
     &      (grid,day_ncep_glob(:,:,itm,itm2),day_ncep(:,:,itm,itm2))
          end do            ;end do
          if(am_i_root())read(kunit,err=10)header,dra_ch4_glob
          do itm=1,max_days ;do itm2=1,nra_ch4
            call unpack_data
     &      (grid,dra_ch4_glob(:,:,itm,itm2),dra_ch4(:,:,itm,itm2))
          end do            ;end do
          if(am_i_root())read(kunit,err=10)header,Rijncep_glob   
          call unpack_data(grid,Rijncep_glob(:,:,:),sum_ncep(:,:,:))
          if(am_i_root())read(kunit,err=10)header,Rijch4_glob 
          call unpack_data(grid,Rijch4_glob(:,:,:),prs_ch4(:,:,:))
          if(am_i_root())read(kunit,err=10)header,HRA_ch4_glob
          do itm=1,maxHR_ch4 ;do itm2=1,nra_ch4
            call unpack_data
     &      (grid,HRA_ch4_glob(:,:,itm,itm2),HRA_ch4(:,:,itm,itm2))
          end do             ;end do
          if(am_i_root())read(kunit,err=10)header,Rijch4_glob
          call unpack_data(grid,Rijch4_glob(:,:,:),r0ch4(:,:,:))
          i0ch4(I_0:I_1,J_0:J_1,:)=NINT(r0ch4(I_0:I_1,J_0:J_1,:))
          if(am_i_root())read(kunit,err=10)header,Rijch4_glob
          call unpack_data(grid,Rijch4_glob(:,:,:),rDch4(:,:,:)) 
          iDch4(I_0:I_1,J_0:J_1,:)=NINT(rDch4(I_0:I_1,J_0:J_1,:))
          if(am_i_root())read(kunit,err=10)header,Rijch4_glob
          call unpack_data(grid,Rijch4_glob(:,:,:),rHch4(:,:,:)) 
          iHch4(I_0:I_1,J_0:J_1,:)=NINT(rHch4(I_0:I_1,J_0:J_1,:))
          if(am_i_root())read(kunit,err=10)header,Rijch4_glob
          call unpack_data(grid,Rijch4_glob(:,:,:),rfirst_mod(:,:,:))
          first_mod(I_0:I_1,J_0:J_1,:)=
     &    NINT(rfirst_mod(I_0:I_1,J_0:J_1,:))
          if(am_i_root())read(kunit,err=10)header,Rijch4_glob
          call unpack_data(grid,Rijch4_glob(:,:,:),avg_model(:,:,:))
          if(am_i_root())read(kunit,err=10)header,Rijncep_glob
          call unpack_data(grid,Rijncep_glob(:,:,:),avg_ncep(:,:,:))
          if(am_i_root())read(kunit,err=10)
     &    header,iday_ncep,i0_ncep,first_ncep
C**** ESMF: Broadcast all non-distributed read arrays.
          call ESMF_BCAST( grid, iday_ncep )  
          call ESMF_BCAST( grid, i0_ncep   )  
          call ESMF_BCAST( grid, first_ncep)
#endif
#endif

          if (am_i_root()) read(kunit,err=10) header
     &         ,trcSurfMixR_acc_glob,trcSurfByVol_acc_glob
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
     &         ,sPM2p5_acc_glob,sPM10_acc_glob,l1PM2p5_acc_glob
     &         ,l1PM10_acc_glob,csPM2p5_acc_glob,csPM10_acc_glob
#endif
#ifdef TRACERS_SPECIAL_Shindell
     &         ,sOx_acc_glob,sNOx_acc_glob,sCO_acc_glob,l1Ox_acc_glob
     &         ,l1NO2_acc_glob
#endif

          call unpack_data(grid,trcSurfMixR_acc_glob,trcSurfMixR_acc)
          call unpack_data(grid,trcSurfByVol_acc_glob,trcSurfByVol_acc)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
          call unpack_data(grid,sPM2p5_acc_glob,sPM2p5_acc)
          call unpack_data(grid,sPM10_acc_glob,sPM10_acc)
          call unpack_data(grid,l1PM2p5_acc_glob,l1PM2p5_acc)
          call unpack_data(grid,l1PM10_acc_glob,l1PM10_acc)
          call unpack_data(grid,csPM2p5_acc_glob,csPM2p5_acc)
          call unpack_data(grid,csPM10_acc_glob,csPM10_acc)
#endif
#ifdef TRACERS_SPECIAL_Shindell
          call unpack_data(grid,sOx_acc_glob,sOx_acc)
          call unpack_data(grid,sNOx_acc_glob,sNOx_acc)
          call unpack_data(grid,sCO_acc_glob,sCO_acc)
          call unpack_data(grid,l1Ox_acc_glob,l1Ox_acc)
          call unpack_data(grid,l1NO2_acc_glob,l1NO2_acc)
#endif

        END SELECT
      END SELECT

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      call io_trDust(kunit,iaction)
#endif

      call freemem
      RETURN

 10   IOERR=1
      call freemem
      RETURN

      contains

      subroutine freemem
      deallocate(TRM_GLOB,TRMOM_GLOB
#ifdef TRACERS_WATER
     &     ,TRWM_GLOB
#endif
     &     )

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      deallocate(snosiz_glob)
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
      deallocate(sPM2p5_acc_glob, sPM10_acc_glob, l1PM2p5_acc_glob,
     & l1PM10_acc_glob, csPM2p5_acc_glob, csPM10_acc_glob)
#endif

      deallocate(Aijl_glob)
#ifdef TRACERS_SPECIAL_Shindell
      deallocate(ss_glob,sOx_acc_glob,sNOx_acc_glob,sCO_acc_glob,
     & l1Ox_acc_glob,l1NO2_acc_glob)
#ifdef INTERACTIVE_WETLANDS_CH4 
      deallocate(day_ncep_glob,DRA_ch4_glob,HRA_ch4_glob,Rijch4_glob,
     & Rijncep_glob,rfirst_mod,rHch4,rDch4,r0ch4)
#endif
#endif
      end subroutine freemem

#endif
      END SUBROUTINE io_tracer

#ifdef NEW_IO
#ifdef TRACERS_ON /* only declare NEW_IO routines when needed */
      subroutine def_rsf_tracer(fid)
!@sum  def_rsf_tracer defines tracer array structure in restart files
!@auth M. Kelley
!@ver  beta
      use model_com, only : coupled_chem
      use tracer_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,
     &yROR,yXO2,yAldehyde,yXO2N,yRXPAR,ss,ydms,yso2,sulfate
     &,acetone,sOx_acc,sNOx_acc,sCO_acc,l1Ox_acc,l1NO2_acc
#ifdef SHINDELL_STRAT_CHEM
     &,SF3,SF2,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2
#endif
#ifdef INTERACTIVE_WETLANDS_CH4 
      use TRACER_SOURCES, only: day_ncep,DRA_ch4,sum_ncep,PRS_ch4,
     & HRA_ch4,iday_ncep,i0_ncep,iHch4,iDch4,i0ch4,first_ncep,first_mod,
     & avg_model,avg_ncep
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE fluxes,ONLY : pprec,pevap
      USE tracers_dust,ONLY : hbaij,ricntd
      use trdust_drv, only: def_rsf_trdust
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      USE AEROSOL_SOURCES, only : snosiz
#endif
      use trdiag_com, only: trcSurfMixR_acc,trcSurfByVol_acc
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
     &     ,sPM2p5_acc,sPM10_acc,l1PM2p5_acc,l1PM10_acc
     &     ,csPM2p5_acc,csPM10_acc
#endif

      implicit none
      integer fid   !@var fid file id
      integer :: n
      character(len=80) :: compstr
      character(len=20) :: ijldims
      ijldims='(dist_im,dist_jm,lm)'
      do n=1,ntm
        call defvar(grid,fid,trm(:,:,:,n),
     &       'trm_'//trim(trname(n))//ijldims)
        call defvar(grid,fid,trmom(:,:,:,:,n),
     &       'trmom_'//trim(trname(n))//'(nmom,dist_im,dist_jm,lm)')
#ifdef TRACERS_WATER
        call defvar(grid,fid,trwm(:,:,:,n),
     &       'trwm_'//trim(trname(n))//ijldims)
#endif
      enddo

#ifdef CUBED_SPHERE
c daily_z is currently only needed for CS
      call defvar(grid,fid,daily_z,'daily_z'//ijldims)
#endif

#ifdef TRACERS_SPECIAL_Shindell       
      compstr='TRACERS_SPECIAL_Shindell'
      call defvar(grid,fid,ss,'ss(JPPJ,lm,dist_im,dist_jm)'
     &     ,defby=compstr)
      call defvar(grid,fid,yNO3,'yNO3'//ijldims,defby=compstr)
      call defvar(grid,fid,pHOx,'pHOx'//ijldims,defby=compstr)
      call defvar(grid,fid,pNOx,'pNOx'//ijldims,defby=compstr)
      call defvar(grid,fid,pOx ,'pOx'//ijldims,defby=compstr)
      call defvar(grid,fid,yCH3O2,'yCH3O2'//ijldims,defby=compstr)
      call defvar(grid,fid,yC2O3,'yC2O3'//ijldims,defby=compstr)
      call defvar(grid,fid,yROR,'yROR'//ijldims,defby=compstr)
      call defvar(grid,fid,yXO2,'yXO2'//ijldims,defby=compstr)
      call defvar(grid,fid,yXO2N,'yXO2N'//ijldims,defby=compstr)
      call defvar(grid,fid,yAldehyde,'yAldehyde'//ijldims,
     &     defby=compstr)
      call defvar(grid,fid,yRXPAR,'yRXPAR'//ijldims,defby=compstr)
      call defvar(grid,fid,ydms,'ydms'//ijldims,defby=compstr)
      call defvar(grid,fid,ySO2,'ySO2'//ijldims,defby=compstr)
      call defvar(grid,fid,sulfate,'sulfate'//ijldims,
     &     defby=compstr)
      call defvar(grid,fid,acetone,'acetone'//ijldims,
     &     defby=compstr)
      if(coupled_chem == 1) then
        call defvar(grid,fid,oh_live,'oh_live'//ijldims,
     &       defby=compstr)
        call defvar(grid,fid,no3_live,'no3_live'//ijldims,
     &       defby=compstr)
      endif
#ifdef SHINDELL_STRAT_CHEM
      compstr='SHINDELL_STRAT_CHEM'
      call defvar(grid,fid,SF3,'SF3'//ijldims,defby=compstr)
      call defvar(grid,fid,SF2,'SF2'//ijldims,defby=compstr)
      call defvar(grid,fid,pClOx,'pClOx'//ijldims,defby=compstr)
      call defvar(grid,fid,pClx,'pClx'//ijldims,defby=compstr)
      call defvar(grid,fid,pOClOx,'pOClOx'//ijldims,defby=compstr)
      call defvar(grid,fid,pBrOx,'pBrOx'//ijldims,defby=compstr)
      call defvar(grid,fid,yCl2,'yCl2'//ijldims,defby=compstr)
      call defvar(grid,fid,yCl2O2,'yCl2O2'//ijldims,defby=compstr)
#endif
#ifdef INTERACTIVE_WETLANDS_CH4 
      compstr='INTERACTIVE_WETLANDS_CH4'
      call defvar(grid,fid,day_ncep,
     &     'day_ncep(dist_im,dist_jm,max_days,nra_ncep)')
      call defvar(grid,fid,dra_ch4,
     &     'dra_ch4(dist_im,dist_jm,max_days,nra_ch4)')
      call defvar(grid,fid,sum_ncep,
     &     'sum_ncep(dist_im,dist_jm,nra_ncep)')
      call defvar(grid,fid,prs_ch4,
     &     'prs_ch4(dist_im,dist_jm,nra_ch4)')
      call defvar(grid,fid,HRA_ch4,
     &     'HRA_ch4(dist_im,dist_jm,maxHR_ch4,nra_ch4)')
      call defvar(grid,fid,i0ch4,
     &     'i0ch4(dist_im,dist_jm,nra_ch4)')
      call defvar(grid,fid,iDch4,
     &     'iDch4(dist_im,dist_jm,nra_ch4)')
      call defvar(grid,fid,iHch4,
     &     'iHch4(dist_im,dist_jm,nra_ch4)')
      call defvar(grid,fid,first_mod,
     &     'first_mod(dist_im,dist_jm,nra_ch4)')
      call defvar(grid,fid,avg_model,
     &     'avg_model(dist_im,dist_jm,nra_ch4)')
      call defvar(grid,fid,avg_ncep, 
     &     'avg_ncep(dist_im,dist_jm,nra_ncep)')
      call defvar(grid,fid,iday_ncep,'iday_ncep(nra_ncep)')
      call defvar(grid,fid,i0_ncep,'i0_ncep(nra_ncep)')
      call defvar(grid,fid,first_ncep,'first_ncep(nra_ncep)')
#endif /* INTERACTIVE_WETLANDS_CH4 */
      call defvar(grid,fid,sOx_acc,'sOx_acc(dist_im,dist_jm)')
      call defvar(grid,fid,sNOx_acc,'sNOx_acc(dist_im,dist_jm)')
      call defvar(grid,fid,sCO_acc,'sCO_acc(dist_im,dist_jm)')
      call defvar(grid,fid,l1Ox_acc,'l1Ox_acc(dist_im,dist_jm)')
      call defvar(grid,fid,l1NO2_acc,'l1NO2_acc(dist_im,dist_jm)')
#endif /* TRACERS_SPECIAL_Shindell */

      call defvar(grid,fid,trcSurfMixR_acc
     &     ,'trcSurfMixR_acc(dist_im,dist_jm,Ntm)')
      call defvar(grid,fid,trcSurfByVol_acc
     &     ,'trcSurfByVol_acc(dist_im,dist_jm,Ntm)')
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
      call defvar(grid,fid,sPM2p5_acc,'sPM2p5_acc(dist_im,dist_jm)')
      call defvar(grid,fid,sPM10_acc,'sPM10_acc(dist_im,dist_jm)')
      call defvar(grid,fid,l1PM2p5_acc,'l1PM2p5_acc(dist_im,dist_jm)')
      call defvar(grid,fid,l1PM10_acc,'l1PM10_acc(dist_im,dist_jm)')
      call defvar(grid,fid,csPM2p5_acc,'csPM2p5_acc(dist_im,dist_jm)')
      call defvar(grid,fid,csPM10_acc,'csPM10_acc(dist_im,dist_jm)')
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      compstr='TRACERS_DUST||TRACERS_MINERALS||TRACERS_QUARZHEM'
      call defvar(grid,fid,hbaij,'hbaij(dist_im,dist_jm)')
      call defvar(grid,fid,ricntd,'ricntd(dist_im,dist_jm)')
      call defvar(grid,fid,pprec,'pprec(dist_im,dist_jm)')
      call defvar(grid,fid,pevap,'pevap(dist_im,dist_jm,nstype)')
      call def_rsf_trdust(fid)
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      call defvar(grid,fid,snosiz,'snosiz(dist_im,dist_jm)')
#endif

      return
      end subroutine def_rsf_tracer

      subroutine new_io_tracer(fid,iaction)
!@sum  new_io_tracer read/write tracer arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite,coupled_chem
      use tracer_com
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,
     &yROR,yXO2,yAldehyde,yXO2N,yRXPAR,ss,ydms,yso2,sulfate
     &,acetone,sOx_acc,sNOx_acc,sCO_acc,l1Ox_acc,l1NO2_acc
#ifdef SHINDELL_STRAT_CHEM
     &,SF3,SF2,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2
#endif
#ifdef INTERACTIVE_WETLANDS_CH4 
      use TRACER_SOURCES, only: day_ncep,DRA_ch4,sum_ncep,PRS_ch4,
     & HRA_ch4,iday_ncep,i0_ncep,iHch4,iDch4,i0ch4,first_ncep,first_mod,
     & avg_model,avg_ncep
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE fluxes,ONLY : pprec,pevap
      USE tracers_dust,ONLY : hbaij,ricntd
      use trdust_drv, only: new_io_trdust
#endif
      use trdiag_com, only: trcSurfMixR_acc,trcSurfByVol_acc
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
     &     ,sPM2p5_acc,sPM10_acc,l1PM2p5_acc,l1PM10_acc
     &     ,csPM2p5_acc,csPM10_acc
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      USE AEROSOL_SOURCES, only : snosiz
#endif
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data,read_data,
     & write_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: n
      select case (iaction)
      case (iowrite)            ! output to restart file
        do n=1,ntm
          call write_dist_data(grid,fid, 'trm_'//trim(trname(n)),
     &         trm(:,:,:,n))
          call write_dist_data(grid,fid, 'trmom_'//trim(trname(n)),
     &         trmom(:,:,:,:,n), jdim=3)
#ifdef TRACERS_WATER
          call write_dist_data(grid,fid, 'trwm_'//trim(trname(n)),
     &         trwm(:,:,:,n))
#endif
        enddo

#ifdef CUBED_SPHERE
c daily_z is currently only needed for CS
        call write_dist_data(grid,fid,'daily_z',daily_z)
#endif

#ifdef TRACERS_SPECIAL_Shindell       
        call write_dist_data(grid,fid,'ss',ss,jdim=4)
        call write_dist_data(grid,fid,'yNO3',yNO3)
        call write_dist_data(grid,fid,'pHOx',pHOx)
        call write_dist_data(grid,fid,'pNOx',pNOx)
        call write_dist_data(grid,fid,'pOx ',pOx)
        call write_dist_data(grid,fid,'yCH3O2',yCH3O2)
        call write_dist_data(grid,fid,'yC2O3',yC2O3)
        call write_dist_data(grid,fid,'yROR',yROR)
        call write_dist_data(grid,fid,'yXO2',yXO2)
        call write_dist_data(grid,fid,'yXO2N',yXO2N)
        call write_dist_data(grid,fid,'yAldehyde',yAldehyde)
        call write_dist_data(grid,fid,'yRXPAR',yRXPAR)
        call write_dist_data(grid,fid,'ydms',ydms)
        call write_dist_data(grid,fid,'ySO2',ySO2)
        call write_dist_data(grid,fid,'sulfate',sulfate)
        call write_dist_data(grid,fid,'acetone',acetone)
        if(coupled_chem == 1) then
          call write_dist_data(grid,fid,'oh_live',oh_live)
          call write_dist_data(grid,fid,'no3_live',no3_live)
        endif
#ifdef SHINDELL_STRAT_CHEM
        call write_dist_data(grid,fid,'SF3',SF3)
        call write_dist_data(grid,fid,'SF2',SF2)
        call write_dist_data(grid,fid,'pClOx',pClOx)
        call write_dist_data(grid,fid,'pClx',pClx)
        call write_dist_data(grid,fid,'pOClOx',pOClOx)
        call write_dist_data(grid,fid,'pBrOx',pBrOx)
        call write_dist_data(grid,fid,'yCl2',yCl2)
        call write_dist_data(grid,fid,'yCl2O2',yCl2O2)
#endif
#ifdef INTERACTIVE_WETLANDS_CH4 
        call write_dist_data(grid,fid,'day_ncep',day_ncep)
        call write_dist_data(grid,fid,'dra_ch4',dra_ch4)
        call write_dist_data(grid,fid,'sum_ncep',sum_ncep)
        call write_dist_data(grid,fid,'prs_ch4',prs_ch4)
        call write_dist_data(grid,fid,'HRA_ch4',HRA_ch4)
        call write_dist_data(grid,fid,'i0ch4',i0ch4)
        call write_dist_data(grid,fid,'iDch4',iDch4)
        call write_dist_data(grid,fid,'iHch4',iHch4)
        call write_dist_data(grid,fid,'first_mod',first_mod)
        call write_dist_data(grid,fid,'avg_model',avg_model)
        call write_dist_data(grid,fid,'avg_ncep',avg_ncep)
        call write_data(grid,fid,'iday_ncep',iday_ncep)
        call write_data(grid,fid,'i0_ncep',i0_ncep)
        call write_data(grid,fid,'first_ncep',first_ncep)
#endif /* INTERACTIVE_WETLANDS_CH4 */
        call write_dist_data(grid,fid,'sOx_acc',sOx_acc)
        call write_dist_data(grid,fid,'sNOx_acc',sNOx_acc)
        call write_dist_data(grid,fid,'sCO_acc',sCO_acc)
        call write_dist_data(grid,fid,'l1Ox_acc',l1Ox_acc)
        call write_dist_data(grid,fid,'l1NO2_acc',l1NO2_acc)
#endif /* TRACERS_SPECIAL_Shindell */

        call write_dist_data(grid,fid,'trcSurfMixR_acc',trcSurfMixR_acc)
        call write_dist_data(grid,fid,'trcSurfByVol_acc'
     &       ,trcSurfByVol_acc)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
        call write_dist_data(grid,fid,'sPM2p5_acc',sPM2p5_acc)
        call write_dist_data(grid,fid,'sPM10_acc',sPM10_acc)
        call write_dist_data(grid,fid,'l1PM2p5_acc',l1PM2p5_acc)
        call write_dist_data(grid,fid,'l1PM10_acc',l1PM10_acc)
        call write_dist_data(grid,fid,'csPM2p5_acc',csPM2p5_acc)
        call write_dist_data(grid,fid,'csPM10_acc',csPM10_acc)
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
        call write_dist_data(grid,fid,'hbaij',hbaij)
        call write_dist_data(grid,fid,'ricntd',ricntd)
        call write_dist_data(grid,fid,'pprec',pprec)
        call write_dist_data(grid,fid,'pevap',pevap)
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
        call write_dist_data(grid,fid,'snosiz',snosiz)
#endif

      case (ioread)            ! input from restart file
        do n=1,ntm
          call read_dist_data(grid,fid, 'trm_'//trim(trname(n)),
     &         trm(:,:,:,n))
          call read_dist_data(grid,fid, 'trmom_'//trim(trname(n)),
     &         trmom(:,:,:,:,n), jdim=3)
#ifdef TRACERS_WATER
          call read_dist_data(grid,fid, 'trwm_'//trim(trname(n)),
     &         trwm(:,:,:,n))
#endif
        enddo

#ifdef CUBED_SPHERE
c daily_z is currently only needed for CS
        call read_dist_data(grid,fid,'daily_z',daily_z)
#endif

#ifdef TRACERS_SPECIAL_Shindell       
        call read_dist_data(grid,fid,'ss',ss,jdim=4)
        call read_dist_data(grid,fid,'yNO3',yNO3)
        call read_dist_data(grid,fid,'pHOx',pHOx)
        call read_dist_data(grid,fid,'pNOx',pNOx)
        call read_dist_data(grid,fid,'pOx ',pOx)
        call read_dist_data(grid,fid,'yCH3O2',yCH3O2)
        call read_dist_data(grid,fid,'yC2O3',yC2O3)
        call read_dist_data(grid,fid,'yROR',yROR)
        call read_dist_data(grid,fid,'yXO2',yXO2)
        call read_dist_data(grid,fid,'yXO2N',yXO2N)
        call read_dist_data(grid,fid,'yAldehyde',yAldehyde)
        call read_dist_data(grid,fid,'yRXPAR',yRXPAR)
        call read_dist_data(grid,fid,'ydms',ydms)
        call read_dist_data(grid,fid,'ySO2',ySO2)
        call read_dist_data(grid,fid,'sulfate',sulfate)
        call read_dist_data(grid,fid,'acetone',acetone)
        if(coupled_chem == 1) then
          call read_dist_data(grid,fid,'oh_live',oh_live)
          call read_dist_data(grid,fid,'no3_live',no3_live)
        endif
#ifdef SHINDELL_STRAT_CHEM
        call read_dist_data(grid,fid,'SF3',SF3)
        call read_dist_data(grid,fid,'SF2',SF2)
        call read_dist_data(grid,fid,'pClOx',pClOx)
        call read_dist_data(grid,fid,'pClx',pClx)
        call read_dist_data(grid,fid,'pOClOx',pOClOx)
        call read_dist_data(grid,fid,'pBrOx',pBrOx)
        call read_dist_data(grid,fid,'yCl2',yCl2)
        call read_dist_data(grid,fid,'yCl2O2',yCl2O2)
#endif
#ifdef INTERACTIVE_WETLANDS_CH4 
        call read_dist_data(grid,fid,'day_ncep',day_ncep)
        call read_dist_data(grid,fid,'dra_ch4',dra_ch4)
        call read_dist_data(grid,fid,'sum_ncep',sum_ncep)
        call read_dist_data(grid,fid,'prs_ch4',prs_ch4)
        call read_dist_data(grid,fid,'HRA_ch4',HRA_ch4)
        call read_dist_data(grid,fid,'i0ch4',i0ch4)
        call read_dist_data(grid,fid,'iDch4',iDch4)
        call read_dist_data(grid,fid,'iHch4',iHch4)
        call read_dist_data(grid,fid,'first_mod',first_mod)
        call read_dist_data(grid,fid,'avg_model',avg_model)
        call read_dist_data(grid,fid,'avg_ncep',avg_ncep)
        call read_data(grid,fid,'iday_ncep',iday_ncep,
     &       bcast_all=.true.)
        call read_data(grid,fid,'i0_ncep',i0_ncep,
     &       bcast_all=.true.)
        call read_data(grid,fid,'first_ncep',first_ncep,
     &       bcast_all=.true.)
#endif /* INTERACTIVE_WETLANDS_CH4 */
        call read_dist_data(grid,fid,'sOx_acc',sOx_acc)
        call read_dist_data(grid,fid,'sNOx_acc',sNOx_acc)
        call read_dist_data(grid,fid,'sCO_acc',sCO_acc)
        call read_dist_data(grid,fid,'l1Ox_acc',l1Ox_acc)
        call read_dist_data(grid,fid,'l1NO2_acc',l1NO2_acc)
#endif /* TRACERS_SPECIAL_Shindell */

        call read_dist_data(grid,fid,'trcSurfMixR_acc',trcSurfMixR_acc)
        call read_dist_data(grid,fid,'trcSurfByVol_acc'
     &       ,trcSurfByVol_acc)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
        call read_dist_data(grid,fid,'sPM2p5_acc',sPM2p5_acc)
        call read_dist_data(grid,fid,'sPM10_acc',sPM10_acc)
        call read_dist_data(grid,fid,'l1PM2p5_acc',l1PM2p5_acc)
        call read_dist_data(grid,fid,'l1PM10_acc',l1PM10_acc)
        call read_dist_data(grid,fid,'csPM2p5_acc',csPM2p5_acc)
        call read_dist_data(grid,fid,'csPM10_acc',csPM10_acc)
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
        call read_dist_data(grid,fid,'hbaij',hbaij)
        call read_dist_data(grid,fid,'ricntd',ricntd)
        call read_dist_data(grid,fid,'pprec',pprec)
        call read_dist_data(grid,fid,'pevap',pevap)
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
        call read_dist_data(grid,fid,'snosiz',snosiz)
#endif

      end select

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      call new_io_trdust(fid,iaction)
#endif

      return
      end subroutine new_io_tracer
#endif /* TRACERS_ON */
#endif /* NEW_IO */

      subroutine num_srf_sources(nt,nsrc,checkname)
!@sum reads headers from emission files to return
!@+ source names and determine the number of sources
!@+ from the number of files in the rundeck of the form:
!@+ trname_##. Then assigns each source to sector(s),
!@+ based on definitions in the rundeck.
!@auth Greg Faluvegi

      use TRACER_COM, only : ntsurfsrcmax,trname,
     & num_tr_sectors,n_max_sect,tr_sect_index,num_sectors,
     & tr_sect_name, sect_name
      USE DOMAIN_DECOMP_ATM, only: GRID, GET, write_parallel
      USE FILEMANAGER, only: openunit,closeunit
      USE PARAM, only : sync_param
    
      implicit none

!@var nsrc number of source to define ntsurfsrc(n)
      integer, intent(out) :: nsrc
      integer, intent(in) :: nt
      integer :: n,iu,i,j,nsect,nn
      character*2 :: fnum
      character*80 :: fname
      character*32 :: pname
      character*124 :: tr_sectors_are
      character(len=300) :: out_line
      logical :: qexist,checkname

! loop through potential number of surface sources, checking if
! those files exist. If they do, obtain the source name by reading
! the header. If not, the number of sources for this tracer has 
! been reached.

      nsrc=0
      loop_n: do n=1,ntsurfsrcmax
        if(n < 10) then ; write(fnum,'(a1,I1)')'0',n
        else ; write(fnum,'(I2)')n ; endif
        fname=trim(trname(nt))//'_'//fnum
        inquire(file=fname,exist=qexist)
        if(qexist) then
          nsrc=nsrc+1
          call openunit(fname,iu,.true.)
          call read_emis_header(nt,n,iu,checkname)
          call closeunit(iu)
! -- begin sector  stuff --
          tr_sectors_are = ' '
          pname=trim(trim(fname)//'_sect')
          call sync_param(pname,tr_sectors_are)
          num_tr_sectors(nt,n)=0
          i=1
          do while(i < len(tr_sectors_are))
            j=index(tr_sectors_are(i:len(tr_sectors_are))," ")
            if (j > 1) then
              num_tr_sectors(nt,n)=num_tr_sectors(nt,n)+1
              i=i+j
            else
              i=i+1
            end if
          enddo
          if(num_tr_sectors(nt,n) > n_max_sect) call stop_model
     &    ("num_tr_sectors problem",255)
          if(num_tr_sectors(nt,n) > 0)then
            read(tr_sectors_are,*)
     &      tr_sect_name(nt,n,1:num_tr_sectors(nt,n))
            do nsect=1,num_tr_sectors(nt,n)
              tr_sect_index(nt,n,nsect)=0
              loop_nn: do nn=1,num_sectors
                if(trim(tr_sect_name(nt,n,nsect)) ==
     &          trim(sect_name(nn))) then
                  tr_sect_index(nt,n,nsect)=nn
                  exit loop_nn
                endif
              enddo loop_nn
            enddo
          endif
! -- end sector stuff --
        else
          exit loop_n
        endif
      enddo loop_n

! and make sure there isn't a skip:

      n=n+1
      if(n < 10) then ; write(fnum,'(a1,I1)')'0',n
      else ; write(fnum,'(I2)')n ; endif
      fname=trim(trname(nt))//'_'//fnum
      inquire(file=fname,exist=qexist)
      if(qexist) then
        write(out_line,*)'problem with num_srf_sources'
        call write_parallel(trim(out_line))
        call stop_model(trim(out_line),255)
      endif

      end subroutine num_srf_sources


      subroutine read_sfc_sources(n,nsrc,xyear,xday,checkname)
!@sum reads surface (2D generally non-interactive) sources
!@auth Jean Lerner/Greg Faluvegi
      USE MODEL_COM, only: itime,im,jm
      USE DOMAIN_DECOMP_ATM, only: GRID, GET,
     &     readt_parallel, write_parallel
      USE FILEMANAGER, only: openunit,closeunit, nameunit
      USE TRACER_COM,only:itime_tr0,trname,sfc_src,ntm,ntsurfsrcmax,
     & freq,nameT,ssname,ty_start,ty_end,delTyr
   
      implicit none
      
      integer :: iu,ns,k,ipos,kx,kstep=10
      integer, intent(in) :: nsrc,n
      character*80 :: fname
      character*2 :: fnum
      character(len=300) :: out_line
      logical,dimension(ntm,ntsurfsrcmax) :: ifirst2=.true.
      real*8 :: alpha
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     & sfc_a,sfc_b
      integer, intent(IN) :: xyear, xday
      logical :: checkname
      
      save ifirst2

      INTEGER :: J_1, J_0, J_0H, J_1H, I_0, I_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      CALL GET(grid, I_STRT=I_0, I_STOP=I_1)
      call GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      if (itime < itime_tr0(n)) return
      if (nsrc <= 0) return
   
      do ns=1,nsrc

! open file and read its header:

        if(ns < 10) then
          write(fnum,'(a1,I1)')'0',ns
        else
          write(fnum,'(I2)')ns
        endif
        fname=trim(trname(n))//'_'//fnum
        call openunit(fname,iu,.true.)
        call read_emis_header(n,ns,iu,checkname)

! now read the data: (should be in kg/m2/s please)
 
! -------------- non-transient emissions ----------------------------!
        if(ty_start(n,ns)==ty_end(n,ns))then 

          select case(freq(n,ns))
          case('a')        ! annual file, only read first time
            if(ifirst2(n,ns)) then
              call readt_parallel(grid,iu,fname,sfc_src(:,:,n,ns),1)
              write(out_line,*)trim(nameT(n,ns)),
     &        ' ',trim(ssname(n,ns)),' ann source read.'
              call write_parallel(trim(out_line))
              ifirst2(n,ns) = .false.
            endif
          case('m')        ! monthly file, interpolate to now
            call read_mon_src_2(n,ns,iu,sfc_src(:,:,n,ns),xyear,xday)
            ifirst2(n,ns) = .false. ! needed?
          end select

! --------------- transient emissions -------------------------------!
        else                
          select case(freq(n,ns))
          case('a')        ! annual file, only read first time + new steps
            kstep=delTyr(n,ns)
            ipos=1
            alpha=0.d0 ! before start year, use start year value
            kx=ty_start(n,ns) ! just for printing
            if(xyear>ty_end(n,ns).or.
     &      (xyear==ty_end(n,ns).and.xday>=183))then
              alpha=1.d0 ! after end year, use end year value     
              ipos=(ty_end(n,ns)-ty_start(n,ns))/kstep
              kx=ty_end(n,ns)-kstep
            endif
            do k=ty_start(n,ns),ty_end(n,ns)-kstep,kstep
!should do!   if(xyear==k .and. xday==183)ifirst2(n,ns)=.true.
              if(xyear>k .or. (xyear==k.and.xday>=183)) then
                if(xyear<k+kstep.or.(xyear==k+kstep.and.xday<183))then
                  ipos=1+(k-ty_start(n,ns))/kstep ! (integer artithmatic)
                  alpha=(365.d0*(0.5+real(xyear-1-k))+xday) / 
     &                  (365.d0*real(kstep))
                  kx=k
                  exit
                endif
              endif
            enddo
!should do! if(ifirst2(n,ns)) then
              call readt_parallel(grid,iu,fname,sfc_a(:,:),ipos)
              call readt_parallel(grid,iu,fname,sfc_b(:,:),1)
!should do! endif
            sfc_src(I_0:I_1,J_0:J_1,n,ns)=sfc_a(I_0:I_1,J_0:J_1)*
     &      (1.d0-alpha)+sfc_b(I_0:I_1,J_0:J_1)*alpha

            write(out_line,'(a,1X,a,a4,F9.4,a16,I4,a8,I4)')
     &      trim(nameT(n,ns)),trim(ssname(n,ns)),' at ',
     &      100.d0*alpha,'% of period mid ',kx,' to mid ',kx+kstep
            call write_parallel(trim(out_line))
            ifirst2(n,ns) = .false.

          case('m')        ! monthly file, interpolate to now
            call read_mon_src_2(n,ns,iu,sfc_src(:,:,n,ns),xyear,xday)
            ifirst2(n,ns) = .false. ! needed?
          end select

        endif
        call closeunit(iu)

      enddo

      return
      end subroutine read_sfc_sources


      subroutine read_emis_header(n,ns,iu,checkname)
!@sum read_emis_header reads the emissions file's header and
!@+   reports back the meta-data. 
!@auth Greg Faluvegi
     
      use TRACER_COM, only : trname,nameT
      USE DOMAIN_DECOMP_ATM, only: GRID, GET, write_parallel
      USE FILEMANAGER, only: openunit,closeunit

      implicit none

      integer, intent(in) :: iu,n,ns
      integer :: error
      character*80 :: message,header
      character(len=300) :: out_line
      logical :: checkname

      error=0
      read(iu)header
      if(len_trim(header) < 1)error=1

      call parse_header(n,ns,trim(header),error)
  
      select case(error)
      case default ! nothing
      case(1) ; message='read_emis_header: missing header'
      case(2) ; message='read_emis_header: problem with freq'
      case(3) ; message='read_emis_header: a and m are choices for freq'
      case(4) ; message='read_emis_header: problem with tracer name'
      case(5) ; message='read_emis_header: tracer name mismatch'
      case(6) ; message='read_emis_header: problem with source'
      case(7) ; message='read_emis_header: problem with res'
      case(8) ; message='read_emis_header: M and F are choices for res'
      case(9) ; message='read_emis_header: transient years seem wrong'
      case(10); message='read_emis_header: delTyr(e.g. kstep) is zero'
      case(11); message='read_emis_header: trans yrs step/years suspect'
      end select
      if(error > 0) then
        if(error == 5 .and. .not. checkname)then
          continue
        else
          write(out_line,*) trim(header)
          call write_parallel(trim(out_line))
          write(out_line,*) trim(message)
          call write_parallel(trim(out_line))
          call stop_model('problem reading emissions',255)
        endif
      endif

      end subroutine read_emis_header


      subroutine parse_header(n,ns,str,error)
!@sum parse_header gets the informantion from emissions file's
!@+  header and reports back the meta-data. 
!@auth Greg Faluvegi
     
      use TRACER_COM, only : trname,freq,nameT,res,ssname,Tyears,
     & ty_start,ty_end,delTyr
      use StringUtilities_mod, only : toLowerCase
      implicit none

      integer, intent(in) :: n,ns
      integer :: n1,n2,error
      character*(*) str

      n1 = scan( str,'=')         ! freq
      if(str(1:n1-1) /= 'freq')error=2
      read(str(n1+1:n1+1),*)freq(n,ns)
      if(freq(n,ns) /= 'a' .and. freq(n,ns) /= 'm')error=3

      str = str(n1+3:)            ! tracer name
      n1 = scan( str, '=')
      if(str(1:n1-1) /= 'name')error=4
      n2 = scan( str,' ')
      read(str(n1+1:n2-1),*)nameT(n,ns)
      if(toLowerCase(trim(trname(n))) /=
     &   toLowerCase(trim(nameT(n,ns))) )error=5

      str = str(n2+1:)            ! source name
      n1 = scan( str,'=')
      if(str(1:n1-1) /= 'source')error=6
      n2 = scan( str,' ')
      read(str(n1+1:n2-1),*)ssname(n,ns)

      str = str(n2+1:)            ! resolution
      n1 = scan( str,'=')
      if(str(1:n1-1) /= 'res')error=7
      read(str(n1+1:n1+1),*)res(n,ns)
#ifndef CUBED_SPHERE
      if(res(n,ns) /= 'M' .and. res(n,ns) /= 'F')error=8
#endif

      n1 = scan( str, ' ')
      str = str(n1+1:)            ! optional transient years
      n1 = scan( str, '=')
      if(str(1:n1-1) /= 'y')then
        ty_start(n,ns)=0; ty_end(n,ns)=0
      else 
        read(str(n1+1:n1+9),*)Tyears(n,ns)
        read(Tyears(n,ns)(1:4),'(I4)')ty_start(n,ns)
        read(Tyears(n,ns)(6:9),'(I4)')ty_end(n,ns)
        if(ty_start(n,ns) /= ty_end(n,ns))then
          if(ty_start(n,ns) < 0 .or. ty_start(n,ns) > 3000)error=9
          if(ty_end(n,ns)   < 0 .or. ty_end(n,ns)   > 3000)error=9
        endif
      endif

      delTyr(n,ns)=10   ! default=decades for backwards compatability
      str = str(n1+9+1:)
      n1 = scan( str, '=')     ! optional transient slice step.
      if(n1 > 0) then
        if(str(n1-3:n1-1) == 'del')then
          str = str(n1+1:)
          n2 = scan( str,' ')
          read(str(1:n2-1),*)delTyr(n,ns)
        endif
      endif

      if(delTyr(n,ns)==0)then
        error=10
      else
        ! check for integer number of slices:
        if(MOD((ty_end(n,ns)-ty_start(n,ns)),delTyr(n,ns))/=0.)error=11
      endif

      end subroutine parse_header


      SUBROUTINE read_mon_src_2(n,ns,iu,data,xyear,xday)
!@sum Read in monthly sources and interpolate to current day
!@+   Calling routine must have the lines:
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      data jdlast /0/
!@+      save jdlast,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
!@+ Note: this started out very similar to read_monthly_sources, 
!@+ but I had problems with that routine. This one is less effecient,
!@+ in that it doesn't keep file open, and has to read monthly files
!@+ each day. My intent is to eventually fix this problem and merge
!@+ with read_monthly_sources. --gsf
!@auth Greg Faluvegi, Jean Lerner and others

      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP_ATM, only : GRID,GET,AM_I_ROOT,write_parallel,
     & READT_PARALLEL, REWIND_PARALLEL, BACKSPACE_PARALLEL
      USE MODEL_COM, only: im,jm,idofm=>JDmidOfM
      USE TRACER_COM, only: ssname,nameT,ty_start,ty_end,delTyr

      implicit none

      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     tlca,tlcb,data
      real*8 :: frac,alpha
      integer ::  imon,iu,n,ns,ipos,k,nn,kx,kstep=10
      character*80 :: junk
      character(len=300) :: out_line
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     & sfc_a,sfc_b
      integer, intent(IN) :: xyear, xday

      integer :: J_0, J_1, I_0, I_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      CALL GET(grid, I_STRT=I_0, I_STOP=I_1)

! -------------- non-transient emissions ----------------------------!
      if(ty_start(n,ns)==ty_end(n,ns))then
        imon=1
        if (xday <= 16)  then ! xDAY in Jan 1-15, first month is Dec
          call readt_parallel(grid,iu,nameunit(iu),tlca,12)
          call rewind_parallel( iu );if (AM_I_ROOT()) read( iu ) junk
        else            ! xDAY is in Jan 16 to Dec 16, get first month
          do while(xday > idofm(imon) .AND. imon <= 12)
            imon=imon+1
          enddo
          call readt_parallel(grid,iu,nameunit(iu),tlca,imon-1)
          if (imon == 13)then
            call rewind_parallel( iu );if (AM_I_ROOT()) read( iu ) junk
          endif
        end if
        call readt_parallel(grid,iu,nameunit(iu),tlcb,1)

c****   Interpolate two months of data to current day
        frac = float(idofm(imon)-xday)/(idofm(imon)-idofm(imon-1))
        data(I_0:I_1,J_0:J_1)=tlca(I_0:I_1,J_0:J_1)*frac + 
     &  tlcb(I_0:I_1,J_0:J_1)*(1.-frac)
        write(out_line,*)
     &  trim(nameT(n,ns)),' ',trim(ssname(n,ns)),' interp to now ',frac
        call write_parallel(trim(out_line))
! --------------- transient emissions -------------------------------!
      else 
        kstep=delTyr(n,ns)
        ipos=1
        alpha=0.d0 ! before start year, use start year value
        kx=ty_start(n,ns) ! just for printing
        if(xyear>ty_end(n,ns).or.
     &  (xyear==ty_end(n,ns).and.xday>=183))then
          alpha=1.d0 ! after end year, use end year value
          ipos=(ty_end(n,ns)-ty_start(n,ns))/kstep
          kx=ty_end(n,ns)-kstep
        endif
        do k=ty_start(n,ns),ty_end(n,ns)-kstep,kstep
          if(xyear>k .or. (xyear==k.and.xday>=183)) then
            if(xyear<k+kstep .or. (xyear==k+kstep.and.xday<183))then
              ipos=1+(k-ty_start(n,ns))/kstep ! (integer artithmatic)
              alpha=real(xyear-k)/real(kstep)
              kx=k
              exit
            endif
          endif
        enddo
!
! read the two necessary months from the first decade:
!
        imon=1
        if (xday <= 16)  then ! xDAY in Jan 1-15, first month is Dec
         call readt_parallel(grid,iu,nameunit(iu),tlca,(ipos-1)*12+12)
         do nn=1,12; call backspace_parallel(iu); enddo
        else            ! xDAY is in Jan 16 to Dec 16, get first month
         do while(xday > idofm(imon) .AND. imon <= 12)
           imon=imon+1
         enddo
         call  readt_parallel
     &   (grid,iu,nameunit(iu),tlca,(ipos-1)*12+imon-1)
         if (imon == 13)then
           do nn=1,12; call backspace_parallel(iu); enddo               
         endif
        end if
        call readt_parallel(grid,iu,nameunit(iu),tlcb,1)
c****   Interpolate two months of data to current day
        frac = float(idofm(imon)-xday)/(idofm(imon)-idofm(imon-1))
        sfc_a(I_0:I_1,J_0:J_1)=tlca(I_0:I_1,J_0:J_1)*frac + 
     &  tlcb(I_0:I_1,J_0:J_1)*(1.-frac)
        call rewind_parallel( iu );if (AM_I_ROOT()) read( iu ) junk

        ipos=ipos+1
        imon=1
        if (xday <= 16)  then ! xDAY in Jan 1-15, first month is Dec
         call readt_parallel(grid,iu,nameunit(iu),tlca,(ipos-1)*12+12)
         do nn=1,12; call backspace_parallel(iu); enddo
        else            ! xDAY is in Jan 16 to Dec 16, get first month
         do while(xday > idofm(imon) .AND. imon <= 12)
           imon=imon+1
         enddo
         call readt_parallel
     &   (grid,iu,nameunit(iu),tlca,(ipos-1)*12+imon-1)
         if (imon == 13)then
           do nn=1,12; call backspace_parallel(iu); enddo
         endif
        end if
        call readt_parallel(grid,iu,nameunit(iu),tlcb,1)
c****   Interpolate two months of data to current day
        frac = float(idofm(imon)-xday)/(idofm(imon)-idofm(imon-1))
        sfc_b(I_0:I_1,J_0:J_1)=tlca(I_0:I_1,J_0:J_1)*frac + 
     &  tlcb(I_0:I_1,J_0:J_1)*(1.-frac)

! now interpolate between the two time periods:

        data(I_0:I_1,J_0:J_1)=sfc_a(I_0:I_1,J_0:J_1)*(1.d0-alpha) + 
     &  sfc_b(I_0:I_1,J_0:J_1)*alpha

        write(out_line,'(a,1X,a,a4,F9.4,a21,I4,a13,I4,a22,F9.4)')
     &  trim(nameT(n,ns)),trim(ssname(n,ns)),' at ',100.d0*alpha,
     &  '% of period this day ',kx,' to this day ',kx+kstep,
     &  ' and monthly fraction=',frac
        call write_parallel(trim(out_line))
 
      endif

      return
      end subroutine read_mon_src_2


      subroutine setup_emis_sectors_regions
!@sum setup_emis_sectors_regions reads from the rundeck the 
!@+ geographic regions and sectors associated with tracer
!@+ emissions and saves names.
!@+ Also reads the factors associated with each sector and
!@+ region. Output IJ map of regions.
!@auth Greg Faluvegi

      use TRACER_COM, only : n_max_sect,reg_n,reg_s,reg_e,reg_w,
     & n_max_reg,alter_sources,ef_REG_IJ,
     & ef_fact,num_regions,num_sectors,sect_name
      USE DOMAIN_DECOMP_ATM, only:GRID,GET,AM_I_ROOT,writet_parallel
      USE GEOM, only: lat2d_dg, lon2d_dg, imaxj
      USE FILEMANAGER, only: openunit,closeunit,nameunit
      USE PARAM, only : sync_param

      implicit none

      integer :: i,j,n,iu
      character*80 :: title
      character*2 :: fnum
      character*124 :: sectors_are,regions_are

      INTEGER J_0, J_1, I_0, I_1
      CALL GET(grid, J_STRT=J_0,J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      sectors_are=' '
      regions_are=' '
      call sync_param("sectors_are",sectors_are)
      call sync_param("regions_are",regions_are)

! see how many sectors there are, save names in array:

      num_sectors=0
      i=1
      do while(i < len(sectors_are))
        j=index(sectors_are(i:len(sectors_are))," ")
        if (j > 1) then
          num_sectors=num_sectors+1
          i=i+j
        else
          i=i+1
        end if
      enddo
      if (num_sectors > n_max_sect) call stop_model
     &("n_max_sect must be increased",255)
      if(num_sectors > 0 ) read(sectors_are,*)
     & sect_name(1:num_sectors)

! see how many regions there are, save names in array:

      num_regions=0
      i=1
      do while(i < len(regions_are))
        j=index(regions_are(i:len(regions_are))," ")
        if (j > 1) then
          num_regions=num_regions+1
          i=i+j
        else
          i=i+1
        end if
      enddo
      if (num_regions > n_max_reg) call stop_model
     &("n_max_reg must be increased",255)
      if(num_regions>0) then
        REG_N(:)=-1.e30 ! initialize to undefined
        REG_S(:)=-1.e30; REG_E(:)=-1.e30; REG_W(:)=-1.e30
        call sync_param("REG_N",REG_N,num_regions)
        call sync_param("REG_S",REG_S,num_regions)
        call sync_param("REG_E",REG_E,num_regions)
        call sync_param("REG_W",REG_W,num_regions)
      end if 

! read the actual emission altering factors:
      do n=1,num_sectors
        if(n < 10) then ; write(fnum,'(a1,I1)')'0',n
        else ; write(fnum,'(I2)')n ; endif
        ef_fact(n,:)=1.d0
        call sync_param("SECT_"//fnum,ef_fact(n,:),num_regions)
      enddo

! check if there is any altering requested:
      alter_sources = .false.
      do i=1,num_sectors; do j=1,num_regions
        if(ef_fact(i,j)/=1.0 .and. ef_fact(i,j)/=-1.d30)
     &  alter_sources = .true.
      enddo ; enddo

! write out regions to GISS-format IJ file, each restart:
      if(alter_sources)then
        ef_REG_IJ(I_0:I_1,J_0:J_1)=0.d0
        do j=J_0,J_1; do i=i_0,imaxj(j); do n=1,num_regions
          if(lat2d_dg(i,j) >= REG_S(n) .and. lat2d_dg(i,j)
     &    <= REG_N(n) .and. lon2d_dg(i,j) >= REG_W(n)
     &    .and. lon2d_dg(i,j) < REG_E(n) ) ef_REG_IJ(i,j)=
     &    max(ef_REG_IJ(i,j),dble(n))
        enddo; enddo; enddo
        title='Regions defined in rundeck for altering tracer sources'
        call openunit('EF_REG',iu,.true.)
        call WRITET_PARALLEL(grid,iu,nameunit(iu),ef_REG_IJ,title)
        call closeunit(iu)
      endif

      end subroutine setup_emis_sectors_regions

#ifdef CACHED_SUBDD
      subroutine tijh_defs(arr,nmax,decl_count)
! 2D tracer outputs (model horizontal grid).
! Each tracer output must be declared separately (no bundling).
      use model_com, only : dtsrc,nday
      use subdd_mod, only : info_type
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      implicit none
      integer :: nmax,decl_count
      type(info_type) :: arr(nmax)

      decl_count = 0

! NONE AT THE MOMENT:
!     arr(next()) = info_type_(
!    &  sname = 'someName',
!    &  lname = 'some tracer output field',
!    &  units = 'kg/m2'
!    &     )
!
      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine tijh_defs

      subroutine tijlh_defs(arr,nmax,decl_count)
! 3D tracer outputs (model horizontal grid and layers).
! Each tracer output must be declared separately (no bundling).
      use model_com, only : dtsrc,nday
      use subdd_mod, only : info_type
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      use tracer_com, only : ntm,trname
      use trdiag_com, only : to_volume_MixRat
      implicit none
      integer :: nmax,decl_count
      integer :: n
      character*80 :: unitString
      type(info_type) :: arr(nmax)

      decl_count = 0

      ! First, diagnostics available for all tracers:
      do n=1,ntm

        ! 3D mixing ratios (SUBDD string is just tracer name):
        if (to_volume_MixRat(n) .eq.1) then
          unitString='mole species / mole air'
        else
          unitString='kg species / kg air'
        endif
        arr(next()) = info_type_(
     &    sname = trim(trname(n)),
     &    lname = trim(trname(n))//' mixing ratio',
     &    units = trim(unitString)
     &    )

      end do ! tracers loop

      ! Other tracer diags on model levels:

#ifdef TRACERS_SPECIAL_Shindell
      arr(next()) = info_type_(
     &  sname = 'MRNO2', ! because not a tracer
     &  lname = 'NO2 mixing ratio',
     &  units = 'mole species / mole air'
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'MRNO', ! because not a tracer
     &  lname = 'NO mixing ratio',
     &  units = 'mole species / mole air'
     &  )
#endif /* TRACERS_SPECIAL_Shindell */

      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine tijlh_defs


      subroutine accumCachedTracerSUBDDs

      use domain_decomp_atm, only : grid
      use model_com, only : LM
      use geom, only : byaxyp
      use dynamics, only : am
      use tracer_com, only : ntm,trname,trm,mass2vol
      use trdiag_com, only : to_volume_MixRat
      use subdd_mod, only : subdd_groups,subdd_type,subdd_ngroups
     &     ,inc_subdd,find_groups, LmaxSUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups)
      type(subdd_type), pointer :: subdd
      integer :: L, n, k
!     real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
!    &                  grid%j_strt_halo:grid%j_stop_halo) :: sddarr2d
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,
     &                  LM                               ) :: sddarr3d
      real*8 :: convert

      call find_groups('taijlh',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        do k=1,subdd%ndiags
          ntm_loop: do n=1,ntm
            ! tracer 3D mixing ratios (SUBDD names are just tracer name):
            if(trim(trname(n)).eq.trim(subdd%name(k))) then
              if (to_volume_MixRat(n) .eq.1) then
                convert=mass2vol(n)
              else
                convert=1.d0
              endif
              do L=1,LmaxSUBDD
                sddarr3d(:,:,L) = 
     &          trm(:,:,L,n)*convert*byaxyp(:,:)/am(L,:,:)
              end do
              call inc_subdd(subdd,k,sddarr3d)
              exit ntm_loop
            end if
         end do ntm_loop
        enddo ! k
      enddo ! igroup

      ! Eventual 2D diags can go below (uncomment sddarr2d declaration
      ! above, if needed)...
!     call find_groups('taijh',grpids,ngroups)
!     do igrp=1,ngroups
!     subdd => subdd_groups(grpids(igrp))
!     do k=1,subdd%ndiags
!     select case (subdd%name(k))
!      ...
!     end select
!     enddo ! k
!     enddo ! igroup

      end subroutine accumCachedTracerSUBDDs
#endif /* CACHED_SUBDD */
