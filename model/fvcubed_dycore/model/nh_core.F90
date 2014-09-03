module nh_core_mod

! Notes:
! Using k_top=2 to treat the top layer hydrostatically so that delz will
! be computed using hydrostatic balance (instead of the update by
! advection of height using extrapolated winds at the model top)
!
! To do list:
! include moisture effect in pt
!------------------------------

   use constants_mod,  only: rdgas, grav
   use fv_control_mod,     only: m_split, quick_p_c, quick_p_d, uniform_ppm,   &
                             k_top, m_riem,  master
   use tp_core_mod,         only: fv_tp_2d, copy_corners
!  use fv_timing_mod,   only: timing_on, timing_off

   implicit none
   private

   public Riem_Solver, Riem_Solver_C, update_dz_c, update_dz_d
   real, parameter:: dz_max = -0.5               ! (meters)

CONTAINS 

  subroutine update_dz_c(is, ie, js, je, km, ng, area,   &
                         zh, ut, vt, dz_in, dz_out, wk)
! !INPUT PARAMETERS:
  integer, intent(in):: is, ie, js, je, ng, km
  real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: ut, vt, zh
  real, intent(in ):: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(in ):: dz_in (is:ie,js:je,km) 
  real, intent(out):: dz_out(is:ie,js:je,km) 
  real, intent(out):: wk(is-ng:ie+ng,js-ng:je+ng,km+1)  ! work array
! Local Work array:
  real, dimension(is:ie+1,js:je  ):: xfx, fx
  real, dimension(is:ie  ,js:je+1):: yfx, fy
  integer  i, j, k

  end subroutine update_dz_c



  subroutine update_dz_d(hord, is, ie, js, je, km, ng, npx, npy, area,    &
                         zh, crx, cry, xfx, yfx, delz, wk, delp, n_sponge)

  integer, intent(in):: is, ie, js, je, ng, km, npx, npy
  integer, intent(in):: hord, n_sponge
  real, intent(in)   :: area(is-ng:ie+ng,js-ng:je+ng)
  real, intent(inout) ::  zh(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout) ::delz(is:ie,js:je,km)
  real, intent(inout) ::delp(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout), dimension(is:ie+1,js-ng:je+ng,km):: crx, xfx
  real, intent(inout), dimension(is-ng:ie+ng,js:je+1,km):: cry, yfx
  real, intent(  out) ::   wk(is:ie,js:je,km)  ! work array
!-----------------------------------------------------
! Local array:
  real, dimension(is:   ie+1, js-ng:je+ng):: crx_adv, xfx_adv
  real, dimension(is-ng:ie+ng,js:   je+1 ):: cry_adv, yfx_adv
  real, dimension(is:ie+1,js:je  ):: fx
  real, dimension(is:ie  ,js:je+1):: fy
  real  delx(is:ie+1,km), dely(is-ng:ie+ng,km)
  real :: ra_x(is:ie,js-ng:je+ng)
  real :: ra_y(is-ng:ie+ng,js:je)
!--------------------------------------------------------------------
  integer  i, j, k, iord, isd, ied, jsd, jed, lm


  end subroutine update_dz_d


  subroutine Riem_Solver_C(dt,   is,  ie,   js, je, km,   ng,  &
                           akap, cp,  ptop, hs, w,  delz, pt,  &
                           delp, gz,  pk,   ip)

   integer, intent(in):: is, ie, js, je, ng, km
   integer, intent(in):: ip       ! ip==1 pk is full pressure
   real, intent(in):: dt,  akap, cp, ptop
   real, intent(in):: delz(is:ie,js:je,km)
   real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng,km):: pt, delp
   real, intent(in)::       hs(is-ng:ie+ng,js-ng:je+ng)
   real, intent(inout):: w(is-ng:ie+ng,js-ng:je+ng,km)
! OUTPUT PARAMETERS 
   real, intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz, pk
! Local:
  real, dimension(is:ie,km  ):: pm, dm, dz2
  real, dimension(is:ie,km+1):: pem, pk2
  real gama, rgrav, ptk
  integer i, j, k
  integer m_split_c

  end subroutine Riem_Solver_C


  subroutine Riem_Solver(dt,   is,   ie,   js, je, km, ng,    &
                         akap, cp,   ptop, hs, peln, w,  delz, pt,  &
                         delp, gz,   pkc, pk, pe, last_call, ip)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: gz: grav*height at edges
!        pe: full     hydrostatic pressure
!       pkc: full non-hydrostatic pressure
!--------------------------------------------
   integer, intent(in):: is, ie, js, je, km, ng
   integer, intent(in):: ip      ! ip==0 pkc is perturbation pressure
   real, intent(in):: dt         ! the BIG horizontal Lagrangian time step
   real, intent(in):: akap, cp, ptop
   real, intent(in):: hs(is-ng:ie+ng,js-ng:je+ng)
   logical, intent(in):: last_call
   real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km):: w, delp, pt
   real, intent(inout):: delz(is:ie,js:je,km)
   real, intent(out), dimension(is-ng:ie+ng,js-ng:je+ng,km+1):: gz, pkc
   real, intent(out):: pk(is:ie,js:je,km+1)
   real, intent(out):: pe(is-1:ie+1,km+1,js-1:je+1)
   real, intent(out):: peln(is:ie,km+1,js:je)           ! ln(pe)
! Local:
  real, dimension(is:ie,km):: pm, dm, dz2
  real :: pem(is:ie,km+1)
  real gama, rgrav, ptk
  integer i, j, k

  end subroutine Riem_Solver


  subroutine Riem_3D(ns, bdt, is, ie, js, je, ng, j, km, cp, gama, cappa, p3, dm2,    &
                     pm2, w, dz2, pt, quick_p, c_core, ktop, iad)

  integer, intent(in):: ns, is, ie, js, je, ng,  km, j
  integer, intent(in):: iad      ! time step scheme 
  integer, intent(in):: ktop     ! starting layer for non-hydrostatic dynamics
                                 ! 1: All non-hydrostatic
                                 ! 2: top sponge layer is hydrostatic
  real,    intent(in):: bdt, cp, gama, cappa
  real,    intent(in), dimension(is:ie,km):: dm2, pm2
  logical, intent(in):: quick_p       ! fast algorithm for pressure
  logical, intent(in):: c_core
  real, intent(in  ) :: pt (is-ng:ie+ng,js-ng:je+ng,km)
! IN/OUT:
  real, intent(inout):: dz2(is:ie,km)
  real, intent(inout)::   w(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(out  )::  p3(is-ng:ie+ng,js-ng:je+ng,km+1)
! --- Local 1D copyies -----
#ifdef USE_2D
  real, dimension(km,is:ie):: t2, p2, pt2
#else
  real, dimension(km):: c2, p2, pt2
#endif
  real, dimension(km):: r_p, r_n, rden, dz, dm, wm, dts, pdt
  real, dimension(km+1):: m_bot, m_top, r_bot, r_top, time_left, pe1, pbar, wbar

  real, parameter:: dzmx = 0.5*dz_max
  real    :: dt, rdt, grg, z_frac, t_left
  real    :: a1, b1, g2, rcp
  real    :: seq(ns)       ! time stepping sequence
  integer :: k2(km+1)
  integer :: i, k, n, ke, kt, k0, k1, k3

 end subroutine Riem_3D

 subroutine time_sequence ( iad, ns, bdt, tseq )
 integer, intent(in) :: iad, ns
 real,    intent(in) :: bdt
 real, intent(out):: tseq(ns)
! local
  integer :: seq(ns)
  integer :: n, nstep
  real :: sdt

 end subroutine time_sequence


 subroutine edge_profile(q1, q2, i1, i2, j1, j2, j, km, limiter, delp)
! Optimized for wind profile reconstruction:
! Developer: S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2, j1, j2
 integer, intent(in):: j, km
 integer, intent(in):: limiter
 real, intent(in), optional::  delp(i1:i2,km)  ! layer thickness
 real, intent(inout), dimension(i1:i2,j1:j2,km):: q1, q2
!-----------------------------------------------------------------------
 real, dimension(i1:i2,km+1):: qe1, qe2  ! edge values
 real   d4(i1:i2)
 real  gam(i1:i2,km)
 real  gak(km)
 real  bet, a_bot, gratio, r2o3, r4o3, xt1, xt2
 integer i, k

 end subroutine edge_profile


 subroutine top_edge(p, qe, is, ie, js, je, km)
! Constant grid spacing:
  integer, intent(in) :: is, ie, js, je, km
  real,    intent(in) ::  p(is:ie,js:je,km)
  real, intent(out)::    qe(is:ie,js:je)
!----------------------------------------
  real dq1, dq2
  real a3, b2, sc
  integer i,j

! three-cell parabolic subgrid distribution at model top
! three-cell PP-distribution
! Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
! a3 = a / 3
! b2 = b / 2

 end subroutine top_edge

 subroutine fix_dz(is, ie, js, je, km, ng, dz)
   integer,  intent(in):: is, ie, km
   integer,  intent(in):: js, je, ng
   real,  intent(inout):: dz(is:ie, js-ng:je+ng,km)
   integer i, j, k
   logical modified

end subroutine fix_dz
!-----------------------------------------------------------------------

#ifdef DEEP_ATM
  subroutine rotate_uvw(dt, im, jm, km, jfirst, jlast, ng_d, ng_s, g_d,  &
                        ua, va, u, v, w, du)
  integer, intent(in):: im, jm, km
  integer, intent(in):: ng_d, ng_s
  integer, intent(in):: jfirst       ! first latitude of the subdomain
  integer, intent(in):: jlast        ! last latitude of the subdomain
  real, intent(in)::  dt
  real, intent(in)::  g_d(im,jfirst:jlast)

  real, intent(inout):: ua(im,jfirst:jlast,km)  ! a-grid u-Wind (m/s)
  real, intent(inout):: va(im,jfirst:jlast,km)  ! a-grid v-Wind (m/s)
  real, intent(inout)::  u(im,jfirst-ng_d:jlast+ng_s,km)  ! u-Wind (m/s)
  real, intent(inout)::  v(im,jfirst-ng_d:jlast+ng_d,km)  ! v-Wind (m/s)
  real, intent(inout)::  w(im,jfirst-ng_d:jlast+ng_d,km)  ! w-Wind (m/s)
  real, intent(out):: du(im,jfirst-1:jlast,km)

! Local:
  real wo(0:im)
  integer i,j,k
  
  end subroutine rotate_uvw
#endif

end module nh_core_mod
