module dyn_core_mod

    use mpp_domains_mod,  only: CGRID_NE, DGRID_NE, mpp_get_boundary,   &
                                mpp_update_domains
    use fv_mp_mod,        only: domain, isd, ied, jsd, jed, is, ie, js, je
    use fv_control_mod,   only: hord_mt, hord_vt, hord_tm, hord_ze, n_sponge,  &
                                dddmp, dddm4, d_ext, m_grad_p, a2b_ord, master
    use nh_core_mod,      only: Riem_Solver_C, Riem_Solver, update_dz_c, update_dz_d
    use sw_core_mod,      only: c_sw, d_sw, divergence_corner
    use a2b_edge_mod,     only: a2b_ord2, a2b_ord4
    use fv_grid_tools_mod,   only: rdx, rdy, rdxc, dxc, dyc, rdyc, dx, dy, area, rarea, grid_type
    use fv_grid_utils_mod,   only: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n,  &
                                   ec1, ec2, en1, en2, da_min_c

!4TAF#ifdef SW_DYNAMICS
!4TAF     use test_cases_mod,  only: test_case, case9_forcing1, case9_forcing2
!4TAF#endif

     use fv_timing_mod, only: timing_on, timing_off
!    use fv_diagnostics_mod, only: prt_maxmin
     use mpp_parameter_mod, only: CORNER

implicit none
private

public :: dyn_core


contains

!-----------------------------------------------------------------------
!     dyn_core :: FV Lagrangian dynamics driver
!-----------------------------------------------------------------------
 
 subroutine dyn_core(npx, npy, npz, ng, bdt, n_split, cp, akap, grav, hydrostatic,  &
                     u,  v,  w, delz, pt, delp, pe, pk, phis, omga, ptop, pfull, ua, va, & 
                     uc, vc, mfx, mfy, cx, cy, pem, delzc, peln, uniform_ppm, time_total)
    integer, intent(IN) :: npx
    integer, intent(IN) :: npy
    integer, intent(IN) :: npz
    integer, intent(IN) :: ng
    integer, intent(IN) :: n_split
    real   , intent(IN) :: bdt
    real   , intent(IN) :: cp, akap, grav
    real   , intent(IN) :: ptop
    logical, intent(IN) :: uniform_ppm
    logical, intent(IN) :: hydrostatic
    real, intent(in) :: pfull(npz)
    real, intent(inout) :: u(   isd:ied  ,jsd:jed+1,npz)  ! D grid zonal wind (m/s)
    real, intent(inout) :: v(   isd:ied+1,jsd:jed  ,npz)  ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   isd:ied  ,jsd:jed  ,npz)  ! vertical vel. (m/s)
    real, intent(inout) :: delz(is :ie   ,js :je   ,npz)  ! delta-height (m)
    real, intent(inout) :: pt(  isd:ied  ,jsd:jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(isd:ied  ,jsd:jed  ,npz)  ! pressure thickness (pascal)
    real, intent(IN), optional:: time_total  ! total time (seconds) since start

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout):: phis(isd:ied,jsd:jed)      ! Surface geopotential (g*Z_surf)
    real, intent(inout):: pe(is-1:ie+1, npz+1,js-1:je+1)  ! edge pressure (pascal)
    real, intent(out) :: pem(is-1:ie+1, npz+1,js-1:je+1)
    real, intent(out):: peln(is:ie,npz+1,js:je)           ! ln(pe)
    real, intent(inout):: pk(is:ie,js:je, npz+1)        ! pe**kappa

!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, intent(inout):: omga(isd:ied,jsd:jed,npz)    ! Vertical pressure velocity (pa/s)
    real, intent(inout):: uc(isd:ied+1,jsd:jed  ,npz)  ! (uc, vc) are mostly used as the C grid winds
    real, intent(inout):: vc(isd:ied  ,jsd:jed+1,npz)
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: ua, va

! The Flux capacitors: accumulated Mass flux arrays
    real, intent(inout)::  mfx(is:ie+1, js:je,   npz)
    real, intent(inout)::  mfy(is:ie  , js:je+1, npz)
! Accumulated Courant number arrays
    real, intent(inout)::  cx(is:ie+1, jsd:jed, npz)
    real, intent(inout)::  cy(isd:ied ,js:je+1, npz)
! Work:
    real, intent(  out):: delzc(is:ie,js:je,npz)  ! 

! Auto 1D & 2D arrays:
    real wbuffer(npy+2,npz)
    real sbuffer(npx+2,npz)
! ----   For external mode:
    real divg2(is:ie+1,js:je+1)
    real wk(isd:ied,jsd:jed)
!-------------------------------------
! Allocatable 3D
    real, allocatable:: delpc(:,:,:)
    real, allocatable::   ptc(:,:,:)
    real, allocatable::   pkc(:,:,:)
    real, allocatable::   pk3(:,:,:)
    real, allocatable::    gz(:,:,:)
    real, allocatable::    zh(:,:,:)
    real, allocatable::    ut(:,:,:)
    real, allocatable::    vt(:,:,:)
    real, allocatable:: crx(:,:,:), xfx(:,:,:)
    real, allocatable:: cry(:,:,:), yfx(:,:,:)
    real, allocatable:: divg_d(:,:,:)

    integer :: hord_m, hord_v, hord_t
    integer :: i,j,k, it
    integer :: ism1, iep1, jsm1, jep1, is2, js2, imax1, jmax1
    integer ieb1, jeb1
    real    :: dt, dt2, rdt, rgrav
    real    :: d_divg
    logical :: do_omega, msg_done, last_step, dord4
    real ptk

    ptk  = ptop ** akap

      allocate( delpc(isd:ied, jsd:jed  ,npz  ) )
      allocate(   ptc(isd:ied, jsd:jed  ,npz  ) )
      allocate(    gz(isd:ied, jsd:jed  ,npz+1) )
      allocate(   pkc(isd:ied, jsd:jed  ,npz+1) )

      allocate( ut(isd:ied, jsd:jed, npz) )
      allocate( vt(isd:ied, jsd:jed, npz) )

! For advection of zh in the D core
      allocate( crx(is :ie+1, jsd:jed,  npz) )
      allocate( xfx(is :ie+1, jsd:jed,  npz) )
      allocate( cry(isd:ied,  js :je+1, npz) )
      allocate( yfx(isd:ied,  js :je+1, npz) )

      if ( dddm4> 0. ) then
           dord4 = .true.
           allocate ( divg_d(isd:ied+1,jsd:jed+1,npz) )
      else
           dord4 = .false.
      endif

      if ( .not. hydrostatic ) then 
           allocate( zh(isd:ied, jsd:jed, npz) )
           if ( m_grad_p==0 ) allocate ( pk3(isd:ied,jsd:jed,npz+1) )
           msg_done = .false.
!          if ( max(npx,npy)>360 )  msg_done = .true.    ! reduce the buffer size
      else
           msg_done = .true.
      endif
                                 call timing_on('COMM_TOTAL')
      if ( npz>1 )   &
      call mpp_update_domains(  pt, domain, complete=.false.)
      call mpp_update_domains(delp, domain, complete=msg_done)

                                 call timing_off('COMM_TOTAL')

      dt  = bdt / real(n_split)
      dt2 = 0.5*dt
      rdt = 1./dt
      rgrav = 1./grav

! Indexes:
      ism1 = is - 1;  iep1 = ie + 1
      jsm1 = js - 1;  jep1 = je + 1
      is2  = max(2,is)
      js2  = max(2,js)
      imax1 = min(npx-1,ie+1)
      jmax1 = min(npy-1,je+1)

! Empty the "flux capacitors"
      mfx(:,:,:) = 0.;  mfy(:,:,:) = 0.
       cx(:,:,:) = 0.;   cy(:,:,:) = 0.

!-----------------------------------------------------
  do it=1,n_split
!    if(master) write(*,*) it
!-----------------------------------------------------
     if ( .not. hydrostatic ) then
        do j=js,je
           do i=is,ie
              zh(i,j,npz) = phis(i,j)*rgrav - delz(i,j,npz)
           enddo
           do k=npz-1,1,-1
              do i=is,ie
                 zh(i,j,k) = zh(i,j,k+1) - delz(i,j,k)
              enddo
           enddo
        enddo
                                 call timing_on('COMM_TOTAL')
        call mpp_update_domains(zh, domain, complete=.false.)
        call mpp_update_domains(w,  domain, complete=.true.)
                                call timing_off('COMM_TOTAL')
     endif

!4TAF#ifdef SW_DYNAMICS
!4TAF      do_omega  = .false.
!4TAF      if (test_case>1) then
!4TAF      if (test_case==9) call case9_forcing1(phis, time_total)
!4TAF#else
      if ( it==n_split ) then
!$omp parallel do private (i, j, k)
      do j=jsm1,jep1
         do i=ism1,iep1
            pem(i,1,j) = ptop
         enddo
         do k=1,npz
            do i=ism1,iep1
               pem(i,k+1,j) = pem(i,k,j) + delp(i,j,k)
            enddo
         enddo
      enddo
           do_omega  = .true.
      else
           do_omega  = .false.
      endif
!4TAF#endif

      if ( it==n_split ) then
           last_step = .true.
      else
           last_step = .false.
      endif

                                                     call timing_on('c_sw')
!$omp parallel do default(shared) private(i,j,k)
      do k=1,npz
         call c_sw(delpc(isd,jsd,k), delp(isd,jsd,k),  ptc(isd,jsd,k),  &
                      pt(isd,jsd,k),    u(isd,jsd,k),    v(isd,jsd,k),  &
                       w(isd,jsd,k),   uc(isd,jsd,k),   vc(isd,jsd,k),  &
                      ua(isd,jsd,k),   va(isd,jsd,k), omga(isd,jsd,k),  &
                      ut(isd,jsd,k),   vt(isd,jsd,k), dt2, hydrostatic, dord4 )
! on output omga is updated w
      enddo
                                                     call timing_off('c_sw')

      if ( hydrostatic ) then
           call geopk(ptop, pe, peln, delpc, pkc, gz, phis, ptc, npz, akap, .false., .false., .true.)
      else
           call update_dz_c(is,   ie, js, je,  npz,    ng,    &
                            area, zh, ut, vt, delz, delzc, gz)
                                               call timing_on('Riem_C')
           call Riem_Solver_C( dt2,   is,  ie,   js,   je,   npz,   ng,   &
                               akap,  cp,  ptop, phis, omga, delzc, ptc,  &
                               delpc, gz,  pkc,  1 )
                                               call timing_off('Riem_C')
! pkc is full non-hydro pressure
                                               call timing_on('COMM_TOTAL')
           call mpp_update_domains(pkc, domain,  whalo=1, ehalo=1,     &
                                        shalo=1, nhalo=1, complete=.false.)
           call mpp_update_domains(gz , domain,  whalo=1, ehalo=1,     &
                                        shalo=1, nhalo=1, complete=.true.)
                                               call timing_off('COMM_TOTAL')
      endif

!-----------------------------------------
! Update time-centered winds on the C-Grid
!-----------------------------------------
!4TAF#ifdef FIX_C_BOUNDARY 
!4TAF      ieb1 = ie;     jeb1 = je
!4TAF#else   
      ieb1 = ie+1;   jeb1 = je+1
!4TAF#endif

!$omp parallel do default(shared) private(i,j,k, wk)
      do k=1,npz
         if ( hydrostatic ) then
              do j=jsm1,jeb1
                 do i=ism1,ieb1
                    wk(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
                 enddo
              enddo
         else
              do j=jsd,jed
                 do i=isd,ied
                       wk(i,j)   = delpc(i,j,k)
                    delpc(i,j,k) =  delp(i,j,k)   ! Save delp for update_dz_d
                 enddo
              enddo
         endif

         do j=js,je
            do i=is,ieb1
               uc(i,j,k) = uc(i,j,k) + dt2*rdxc(i,j) / (wk(i-1,j)+wk(i,j)) *   &
                      ( (gz(i-1,j,k+1)-gz(i,j,k  ))*(pkc(i,j,k+1)-pkc(i-1,j,k))  &
                      + (gz(i-1,j,k) - gz(i,j,k+1))*(pkc(i-1,j,k+1)-pkc(i,j,k)) )
            enddo
         enddo
         do j=js,jeb1
            do i=is,ie
               vc(i,j,k) = vc(i,j,k) + dt2*rdyc(i,j) / (wk(i,j-1)+wk(i,j)) *   &
                      ( (gz(i,j-1,k+1)-gz(i,j,k  ))*(pkc(i,j,k+1)-pkc(i,j-1,k))  &
                      + (gz(i,j-1,k) - gz(i,j,k+1))*(pkc(i,j-1,k+1)-pkc(i,j,k)) )
            enddo
         enddo
      enddo

!--------------------------------------------------------------------------------------------
!4TAF#ifdef FIX_C_BOUNDARY 
!4TAF                                                                 call timing_on('COMM_TOTAL')
!4TAF      call mpp_get_boundary(uc, vc, domain, wbufferx=wbuffer, ebufferx=uc(ie+1,js:je,1:npz), &
!4TAF                          sbuffery=sbuffer, nbuffery=vc(is:ie,je+1,1:npz), gridtype=CGRID_NE )
!4TAF                                                                 call timing_off('COMM_TOTAL')
!4TAF#endif
!--------------------------------------------------------------------------------------------
                                                     call timing_on('COMM_TOTAL')
      call mpp_update_domains( uc, vc, domain, gridtype=CGRID_NE, complete=.true.)
                                                     call timing_off('COMM_TOTAL')
!4TAF#ifdef SW_DYNAMICS
!4TAF      if (test_case==9) call case9_forcing2(phis)
!4TAF      else
!4TAF         last_step = .true.
!4TAF      endif !test_case>1
!4TAF#endif

    if ( dord4 ) then
!                                           call timing_on('  divg_d')
         call divergence_corner(u, v, ua, va, divg_d, npz)
         if ( grid_type/=4 ) then
                                            call timing_on('COMM_TOTAL')
              call mpp_update_domains(divg_d, domain, position=CORNER)
                                            call timing_off('COMM_TOTAL')
         endif
!                                           call timing_off('  divg_d')
    endif

                                                     call timing_on('d_sw')
!$omp parallel do default(shared) private(i,j,k, d_divg, hord_m, hord_v, hord_t, wk)
    do k=1,npz
       hord_m = hord_mt
       hord_t = hord_tm
       hord_v = hord_vt
       if ( n_sponge==-1 .or. npz==1 ) then
! Constant divg damping coefficient:
           d_divg = dddmp
       else
#ifndef NEW_CODE
           if( k <= n_sponge .and. npz>16 ) then
! Apply first order scheme for damping the sponge layer
!              hord_m = 1
               hord_v = 1
               hord_t = 1
               d_divg = min(0.24, 6.*dddmp)   ! 0.25 is the stability limit
               d_divg = max(0.02,   d_divg)
           elseif( k == n_sponge+1 .and. npz>16 ) then
               d_divg = min(0.24, 4.*dddmp)
               d_divg = max(0.01 ,  d_divg)
           else
               d_divg = min(0.24, dddmp*(1.-3.*tanh(0.1*log(pfull(k)/pfull(npz)))))
           endif
#else
           if( k <= n_sponge .and. npz>16 ) then
! Apply first order scheme for damping the sponge layer
!              hord_m = 1
               d_divg = min(0.24, 20.*dddmp)   ! 0.25 is the stability limit
               d_divg = max(0.02,   d_divg)
           elseif( k == n_sponge+1 .and. npz>16 ) then
               d_divg = min(0.24, 10.*dddmp)
               d_divg = max(0.01 ,  d_divg)
           elseif( k == n_sponge+2 .and. npz>16 ) then
               d_divg = min(0.24, 4.*dddmp)
               d_divg = max(0.01 ,  d_divg)
           else
               d_divg = min(0.24, dddmp*(1.-3.*tanh(0.1*log(pfull(k)/pfull(npz)))))
           endif
#endif
       endif

!--- external mode divergence damping ---
       if ( d_ext > 0. )  &
            call a2b_ord2(delp(isd,jsd,k), wk, npx, npy, is,    &
                          ie, js, je, ng, .false.)

       call d_sw( vt(isd,jsd,k), delp(isd,jsd,k), ptc(isd,jsd,k),  pt(isd,jsd,k), &
                   u(isd,jsd,k),    v(isd,jsd,k),   w(isd,jsd,k),  uc(isd,jsd,k), &
                  vc(isd,jsd,k),   ua(isd,jsd,k),  va(isd,jsd,k), divg_d(isd,jsd,k),  &
                 mfx(is, js, k),  mfy(is, js, k),  cx(is, jsd,k),  cy(isd,js, k), &
                 crx(is, jsd,k),  cry(isd,js, k), xfx(is, jsd,k), yfx(isd,js, k), &
                 dt, hord_m, hord_v, hord_t, d_divg, dddm4, hydrostatic, uniform_ppm )

       if ( d_ext > 0. ) then
            do j=js,jep1
               do i=is,iep1
                  ptc(i,j,k) = wk(i,j)
               enddo
            enddo
       endif
    enddo         
                                                     call timing_off('d_sw')

    if ( d_ext > 0. ) then
!$omp parallel do default(shared) private(i,j,k)
! Barotropic mode:
          d_divg = d_ext * da_min_c
          do j=js,jep1
              do i=is,iep1
                    wk(i,j) = ptc(i,j,1)
                 divg2(i,j) = wk(i,j)*vt(i,j,1)
              enddo
              do k=2,npz
                 do i=is,iep1
                       wk(i,j) =    wk(i,j) + ptc(i,j,k)
                    divg2(i,j) = divg2(i,j) + ptc(i,j,k)*vt(i,j,k)
                 enddo
              enddo
              do i=is,iep1
                 divg2(i,j) = d_divg*divg2(i,j)/wk(i,j)
              enddo
          enddo
    else
        divg2 = 0.
        vt = 0.
    endif
                               call timing_on('COMM_TOTAL')
    call mpp_update_domains(  pt, domain, complete=.false.)
    call mpp_update_domains(delp, domain, complete=.true.)
                             call timing_off('COMM_TOTAL')
     if ( hydrostatic ) then
          call geopk(ptop, pe, peln, delp, pkc, gz, phis, pt, npz, akap, last_step, .false., .false.)
     else
                                            call timing_on('UPDATE_DZ')
          call update_dz_d(hord_tm, is, ie, js, je, npz, ng, npx, npy, area,  &
                           zh, crx, cry, xfx, yfx, delz, delzc, delpc, n_sponge)
                                            call timing_off('UPDATE_DZ')
                                                          call timing_on('Riem_D')
!-----------------------------------------------------------
! mgrad_p = 1: pkc is full pressure
! mgrad_p = 0: pkc is non-hydrostatic perturbation pressure
!-----------------------------------------------------------
          call Riem_Solver(dt,   is,   ie,   js,   je, npz,  ng,  &
                           akap, cp,   ptop, phis, peln, w,  delz,      &
                           pt,   delp, gz,   pkc,  pk, pe, last_step, m_grad_p)
                                                 call timing_off('Riem_D')

                                       call timing_on('COMM_TOTAL')
          if ( m_grad_p==0 ) then
             do k=1,npz+1
                do j=js,je
                   do i=is,ie
                      pk3(i,j,k) = pk(i,j,k)
                   enddo
                enddo
             enddo
             if ( a2b_ord==4 ) then
                call mpp_update_domains(pk3, domain, whalo=2, ehalo=2,   &
                                        shalo=2, nhalo=2, complete=.false.)
!               call mpp_update_domains(pk3, domain, complete=.false.)
             else
                call mpp_update_domains(pk3, domain, whalo=1, ehalo=1,   &
                                        shalo=1, nhalo=1, complete=.false.)
             endif
          endif

          if ( a2b_ord==4 ) then
               call mpp_update_domains(pkc, domain, whalo=2, ehalo=2,   &
                                       shalo=2, nhalo=2, complete=.false.)
               call mpp_update_domains(gz , domain, whalo=2, ehalo=2,   &
                                       shalo=2, nhalo=2, complete=.true.)
!              call mpp_update_domains(pkc, domain, complete=.true.)
!              call mpp_update_domains(gz , domain, complete=.true.)
          else
               call mpp_update_domains(pkc, domain, whalo=1, ehalo=1,   &
                                       shalo=1, nhalo=1, complete=.false.)
               call mpp_update_domains(gz , domain, whalo=1, ehalo=1,   &
                                       shalo=1, nhalo=1, complete=.true.)
          endif
                                       call timing_off('COMM_TOTAL')
     endif    ! end hydro case


!4TAF#ifdef SW_DYNAMICS
!4TAF      if (test_case > 1) then
!4TAF#else
      if ( last_step .and. hydrostatic ) then
!$omp parallel do private(i, j, k)
           do k=1,npz+1
              do j=js,je
                 do i=is,ie
                    pk(i,j,k) = pkc(i,j,k)
                 enddo
              enddo
           enddo
      endif
    
      if ( do_omega ) then
!------------------------------
! Compute time tendency
!------------------------------
         do k=1,npz
            do j=js,je
               do i=is,ie
                  omga(i,j,k) = (pe(i,k+1,j) - pem(i,k+1,j)) * rdt 
               enddo
            enddo
         enddo
!------------------------------
! Compute the "advective term"
!------------------------------
         call adv_pe(ua, va, pem, omga, npx, npy,  npz, ng)
      endif
!4TAF#endif

      if ( .not.hydrostatic .and. m_grad_p == 0 ) then
           call two_grad_p(u, v, pkc, gz, delp, pk3, divg2, dt, ng, npx, npy,   &
                           npz, ptk)  
      else
           call one_grad_p(u, v, pkc, gz, divg2, delp, dt, ng, npx, npy, npz,   &
                           ptop, ptk, hydrostatic)  
      endif

                                                                call timing_on('COMM_TOTAL')
!     if( last_step ) then
      if( last_step .and. grid_type<4 ) then
! Prevent accumulation of rounding errors at overlapped domain edges:
          call mpp_get_boundary(u, v, domain, wbuffery=wbuffer, ebuffery=v(ie+1,js:je,1:npz),  &
                            sbufferx=sbuffer, nbufferx=u(is:ie,je+1,1:npz), gridtype=DGRID_NE )
      else
          call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
      endif
                                                                call timing_off('COMM_TOTAL')
!4TAF#ifdef SW_DYNAMICS
!4TAF      endif
!4TAF#endif

!-----------------------------------------------------
  enddo   ! time split loop
!-----------------------------------------------------

   deallocate(    gz )
   deallocate(   pkc )
   deallocate( delpc )
   deallocate(   ptc )
   deallocate(    ut )
   deallocate(    vt )
   deallocate(   crx )
   deallocate(   xfx )
   deallocate(   cry )
   deallocate(   yfx )
   if ( dord4 ) deallocate ( divg_d )

   if ( .not. hydrostatic ) then
         deallocate( zh )
         if ( m_grad_p==0 ) deallocate ( pk3 )
   endif

 end subroutine dyn_core



 subroutine two_grad_p(u, v, pkc, gz, delp, pk3, divg2, dt, ng, npx, npy, npz, ptk)  

    integer, intent(IN) :: ng, npx, npy, npz
    real,    intent(IN) :: dt, ptk
    real,    intent(in) :: divg2(is:ie+1, js:je+1)
    real, intent(inout) ::  delp(isd:ied, jsd:jed, npz)
    real, intent(inout) ::   pkc(isd:ied, jsd:jed, npz+1)  ! perturbation pressure
    real, intent(inout) ::   pk3(isd:ied, jsd:jed, npz+1)  ! p**kappa
    real, intent(inout) ::    gz(isd:ied, jsd:jed, npz+1)  ! g * zh
    real, intent(inout) ::     u(isd:ied,  jsd:jed+1,npz) 
    real, intent(inout) ::     v(isd:ied+1,jsd:jed,  npz)
! Local:
    real wk1(isd:ied, jsd:jed)
    real  wk(is: ie+1,js: je+1)
    integer iep1, jep1
    integer i,j,k

    iep1 = ie + 1
    jep1 = je + 1

    do j=js,jep1
       do i=is,iep1
          pkc(i,j,1) = 0.
          pk3(i,j,1) = ptk
       enddo
    enddo

    do k=1,npz+1

       if ( k/=1 ) then
         if ( a2b_ord==4 ) then
           call a2b_ord4(pkc(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
           call a2b_ord4(pk3(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
         else
           call a2b_ord2(pkc(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
           call a2b_ord2(pk3(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
         endif
       endif

       if ( a2b_ord==4 ) then
           call a2b_ord4( gz(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
       else
           call a2b_ord2( gz(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng, .true.)
       endif
    enddo

    do k=1,npz

       if ( a2b_ord==4 ) then
            call a2b_ord4(delp(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng)
       else
            call a2b_ord2(delp(isd,jsd,k), wk1, npx, npy, is, ie, js, je, ng)
       endif

       do j=js,jep1
          do i=is,iep1
             wk(i,j) = pk3(i,j,k+1) - pk3(i,j,k)
          enddo
       enddo

       do j=js,jep1
          do i=is,ie
!------------------
! Perturbation term:
!------------------
             u(i,j,k) = u(i,j,k) + dt/(wk1(i,j)+wk1(i+1,j)) *   &
                   ((gz(i,j,k+1)-gz(i+1,j,k))*(pkc(i+1,j,k+1)-pkc(i,j,k)) &
                  + (gz(i,j,k)-gz(i+1,j,k+1))*(pkc(i,j,k+1)-pkc(i+1,j,k)))
!-----------------
! Hydrostatic term
!-----------------
!4TAF#ifdef NO_EXT_DAMP
!4TAF             u(i,j,k) = rdx(i,j)*(u(i,j,k) + dt/(wk(i,j)+wk(i+1,j)) *   &
!4TAF#else
             u(i,j,k) = rdx(i,j)*(divg2(i,j)-divg2(i+1,j)+u(i,j,k) + dt/(wk(i,j)+wk(i+1,j)) *      &
!4TAF#endif
                   ((gz(i,j,k+1)-gz(i+1,j,k))*(pk3(i+1,j,k+1)-pk3(i,j,k)) &
                  + (gz(i,j,k)-gz(i+1,j,k+1))*(pk3(i,j,k+1)-pk3(i+1,j,k))))
          enddo
       enddo

       do j=js,je
          do i=is,iep1
!------------------
! Perturbation term:
!------------------
             v(i,j,k) = v(i,j,k) + dt/(wk1(i,j)+wk1(i,j+1)) *   &
                   ((gz(i,j,k+1)-gz(i,j+1,k))*(pkc(i,j+1,k+1)-pkc(i,j,k)) &
                  + (gz(i,j,k)-gz(i,j+1,k+1))*(pkc(i,j,k+1)-pkc(i,j+1,k)))
!-----------------
! Hydrostatic term
!-----------------
!4TAF#ifdef NO_EXT_DAMP
!4TAF             v(i,j,k) = rdy(i,j)*(v(i,j,k) + dt/(wk(i,j)+wk(i,j+1)) *   &
!4TAF#else
             v(i,j,k) = rdy(i,j)*(divg2(i,j)-divg2(i,j+1)+v(i,j,k) + dt/(wk(i,j)+wk(i,j+1)) *      &
!4TAF#endif
                   ((gz(i,j,k+1)-gz(i,j+1,k))*(pk3(i,j+1,k+1)-pk3(i,j,k)) &
                  + (gz(i,j,k)-gz(i,j+1,k+1))*(pk3(i,j,k+1)-pk3(i,j+1,k))))
          enddo
       enddo
    enddo    ! end k-loop

 end subroutine two_grad_p



 subroutine one_grad_p(u, v, pkc, gz, divg2, delp, dt, ng, npx, npy, npz,  &
                       ptop, ptk,  hydrostatic)  

    integer, intent(IN) :: ng, npx, npy, npz
    real,    intent(IN) :: dt, ptop, ptk
    logical, intent(in) :: hydrostatic
    real,    intent(in) :: divg2(is:ie+1,js:je+1)
    real, intent(inout) ::   pkc(isd:ied,  jsd:jed  ,npz+1)
    real, intent(inout) ::    gz(isd:ied,  jsd:jed  ,npz+1)
    real, intent(inout) ::  delp(isd:ied,  jsd:jed  ,npz)
    real, intent(inout) ::     u(isd:ied  ,jsd:jed+1,npz) 
    real, intent(inout) ::     v(isd:ied+1,jsd:jed  ,npz)
! Local:
    real, dimension(isd:ied,jsd:jed):: wk
    real top_value
    integer :: iep1, jep1
    integer i,j,k

    iep1 = ie + 1
    jep1 = je + 1


    if ( hydrostatic ) then
! pkc is pe**kappa if hydrostatic
         top_value = ptk
    else
! pkc is full pressure if non-hydrostatic
         top_value = ptop
    endif

    do j=js,jep1
       do i=is,iep1
          pkc(i,j,1) = top_value
       enddo
    enddo

    do k=2,npz+1
       if ( a2b_ord==4 ) then
         call a2b_ord4(pkc(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       else
         call a2b_ord2(pkc(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       endif
    enddo

    do k=1,npz+1
       if ( a2b_ord==4 ) then
         call a2b_ord4( gz(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       else
         call a2b_ord2( gz(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng, .true.)
       endif
    enddo

    do k=1,npz

       if ( hydrostatic ) then
            do j=js,jep1
               do i=is,iep1
                  wk(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
               enddo
            enddo
       else
         if ( a2b_ord==4 ) then
            call a2b_ord4(delp(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng)
         else
            call a2b_ord2(delp(isd,jsd,k), wk, npx, npy, is, ie, js, je, ng)
         endif
       endif

       do j=js,jep1
          do i=is,ie
!4TAF#ifdef NO_EXT_DAMP
!4TAF             u(i,j,k) = rdx(i,j)*(u(i,j,k) + dt/(wk(i,j)+wk(i+1,j)) *   &
!4TAF#else
             u(i,j,k) = rdx(i,j)*(divg2(i,j)-divg2(i+1,j)+u(i,j,k) + dt/(wk(i,j)+wk(i+1,j)) *      &
!4TAF#endif
                   ((gz(i,j,k+1)-gz(i+1,j,k))*(pkc(i+1,j,k+1)-pkc(i,j,k)) &
                  + (gz(i,j,k)-gz(i+1,j,k+1))*(pkc(i,j,k+1)-pkc(i+1,j,k))))
          enddo
       enddo
       do j=js,je
          do i=is,iep1
!4TAF#ifdef NO_EXT_DAMP
!4TAF             v(i,j,k) = rdy(i,j)*(v(i,j,k) + dt/(wk(i,j)+wk(i,j+1)) *   &
!4TAF#else
             v(i,j,k) = rdy(i,j)*(divg2(i,j)-divg2(i,j+1)+v(i,j,k) + dt/(wk(i,j)+wk(i,j+1)) *      &
!4TAF#endif
                   ((gz(i,j,k+1)-gz(i,j+1,k))*(pkc(i,j+1,k+1)-pkc(i,j,k)) &
                  + (gz(i,j,k)-gz(i,j+1,k+1))*(pkc(i,j,k+1)-pkc(i,j+1,k))))
          enddo
       enddo
    enddo    ! end k-loop

 end subroutine one_grad_p



 subroutine geopk(ptop, pe, peln, delp, pk, gz, hs, pt, km, akap, last_call, dp_check, CG)

     integer, intent(IN) :: km
     real   , intent(IN) :: akap, ptop
     real   , intent(IN) :: hs(isd:ied,jsd:jed)
     real, intent(INOUT), dimension(isd:ied,jsd:jed,km):: pt, delp
     logical, intent(IN) :: last_call, dp_check, CG
! !OUTPUT PARAMETERS
     real, intent(OUT), dimension(isd:ied,jsd:jed,km+1):: gz, pk
     real, intent(OUT)::   pe(is-1:ie+1,km+1,js-1:je+1)
     real, intent(out):: peln(is:ie,km+1,js:je)           ! ln(pe)
! !DESCRIPTION:
!    Calculates geopotential and pressure to the kappa.
! Local:
     real  p1d(is-2:ie+2)
     real logp(is-2:ie+2)
     real ptk, dp, dpmin
     integer i, j, k
     integer ifirst, ilast
     integer jfirst, jlast

     dpmin = 0.01*ptop
     ptk  = ptop ** akap

          ifirst = is-1
          jfirst = js-1

!4TAF#ifdef FIX_C_BOUNDARY 
!4TAF     if ( CG ) then   ! C-Grid
!4TAF          ilast  = ie
!4TAF          jlast  = je
!4TAF     else
!4TAF          ilast  = ie+1
!4TAF          jlast  = je+1
!4TAF     endif
!4TAF#else
          ilast  = ie+1
          jlast  = je+1
!4TAF#endif

     if ( .not. CG .and. a2b_ord==4 ) then   ! D-Grid
          ifirst = is-2; ilast = ie+2
          jfirst = js-2; jlast = je+2
     endif

!$omp parallel do default(shared) private(i,j,k, p1d, dp)
     do 2000 j=jfirst,jlast

        do i=ifirst, ilast
           p1d(i) = ptop
           pk(i,j,1) = ptk
           gz(i,j,km+1) = hs(i,j)
        enddo

        if( last_call .and. j>(js-2) .and. j<(je+2) ) then
           do i=max(ifirst,is-1), min(ilast,ie+1) 
              pe(i,1,j) = ptop
           enddo
        endif

!4TAF#ifdef DP_CHECK
!4TAF        if( dp_check ) then
!4TAF          do k=1, km-1
!4TAF             do i=ifirst, ilast
!4TAF              if(delp(i,j,k) < dpmin) then
! Remap from below and mix pt
!4TAF                dp = dpmin - delp(i,j,k)
!4TAF                pt(i,j,k) = (pt(i,j,k)*delp(i,j,k) + pt(i,j,k+1)*dp) / dpmin
!4TAF                delp(i,j,k) = dpmin
!4TAF                delp(i,j,k+1) = delp(i,j,k+1) - dp
!4TAF              endif
!4TAF            enddo
!4TAF          enddo

! Bottom (k=km):
!4TAF          do i=ifirst, ilast
!4TAF            if(delp(i,j,km) < dpmin) then
! Remap from above and mix pt
!4TAF              dp = dpmin - delp(i,j,km)
!4TAF              pt(i,j,km) = (pt(i,j,km)*delp(i,j,km) + pt(i,j,km-1)*dp)/dpmin
!4TAF              delp(i,j,km) = dpmin
!4TAF              delp(i,j,km-1) = delp(i,j,km-1) - dp
!4TAF            endif
!4TAF          enddo
!4TAF        endif
!4TAF#endif

! Top down
        do k=2,km+1
          do i=ifirst, ilast
               p1d(i)  = p1d(i) + delp(i,j,k-1)
!            pk(i,j,k) = p1d(i) ** akap
! Optimized form:
              logp(i)  = log(p1d(i))
             pk(i,j,k) = exp( akap*logp(i) )
          enddo

          if( last_call .and. j>(js-2) .and. j<(je+2) ) then
             do i=max(ifirst,is-1), min(ilast,ie+1) 
                pe(i,k,j) = p1d(i)
             enddo
             if ( j>=js .and. j<=je ) then
                  do i=is,ie
                     peln(i,k,j) = logp(i)
                  enddo
             endif
          endif

        enddo

! Bottom up
        do k=km,1,-1
           do i=ifirst, ilast
              gz(i,j,k) = gz(i,j,k+1) + pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
           enddo
        enddo
2000  continue

 end subroutine geopk

 
 subroutine adv_pe(ua, va, pem, om, npx, npy, npz, ng)

 integer, intent(in) :: npx, npy, npz, ng
! Contra-variant wind components:
 real, intent(in), dimension(isd:ied,jsd:jed,npz):: ua, va
! Pressure at edges:
 real, intent(in) :: pem(is-1:ie+1,1:npz+1,js-1:je+1)
 real, intent(inout) :: om(isd:ied,jsd:jed,npz)

! Local:
 real, dimension(is:ie,js:je):: ut, vt
 real v3(3,is:ie,js:je)

 real pin(isd:ied,jsd:jed)
 real  pb(isd:ied,jsd:jed)

 real grad(3,is:ie,js:je)
 real pdx(3,is:ie,js:je+1)
 real pdy(3,is:ie+1,js:je)
 integer :: i,j,k, n

!$omp parallel do private (i, j, k, n, pdx, pdy, pin, pb, ut, vt, grad, v3)
 do k=1,npz
    if ( k==npz ) then
       do j=js,je
          do i=is,ie
             ut(i,j) = ua(i,j,npz)
             vt(i,j) = va(i,j,npz)
          enddo
       enddo
    else
       do j=js,je
          do i=is,ie
             ut(i,j) = 0.5*(ua(i,j,k)+ua(i,j,k+1))
             vt(i,j) = 0.5*(va(i,j,k)+va(i,j,k+1))
          enddo
       enddo
    endif

! Compute Vect wind:
    do j=js,je
       do i=is,ie
          do n=1,3
             v3(n,i,j) = ut(i,j)*ec1(n,i,j) + vt(i,j)*ec2(n,i,j) 
          enddo
       enddo
    enddo

    do j=js-1,je+1
       do i=is-1,ie+1
          pin(i,j) = pem(i,k+1,j)
       enddo
    enddo

! Compute pe at 4 cell corners:
    call a2b_ord2(pin, pb, npx, npy, is, ie, js, je, ng)


    do j=js,je+1
       do i=is,ie
          do n=1,3
             pdx(n,i,j) = (pb(i,j)+pb(i+1,j))*dx(i,j)*en1(n,i,j)
          enddo
       enddo
    enddo
    do j=js,je
       do i=is,ie+1
          do n=1,3
             pdy(n,i,j) = (pb(i,j)+pb(i,j+1))*dy(i,j)*en2(n,i,j)
          enddo
       enddo
    enddo

! Compute grad (pe) by Green's theorem
    do j=js,je
       do i=is,ie
          do n=1,3
             grad(n,i,j) = pdx(n,i,j+1) - pdx(n,i,j) - pdy(n,i,j) + pdy(n,i+1,j)
          enddo
       enddo
    enddo

! Compute inner product: V3 * grad (pe)
       do j=js,je
          do i=is,ie
             om(i,j,k) = om(i,j,k) + 0.5*rarea(i,j)*(v3(1,i,j)*grad(1,i,j) +   &
                         v3(2,i,j)*grad(2,i,j) + v3(3,i,j)*grad(3,i,j))
          enddo
       enddo
 enddo

 end subroutine adv_pe 

end module dyn_core_mod
