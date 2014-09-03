#include "rundeck_opts.h"

#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
!@sum  OCN_TRACER: tracer-dependent routines for GISS Ocean tracers
!@+    Routines included:
!@+      Those that MUST EXIST for all tracers:
!@+        Tracer initialisation + sources: tracer_ic_ocean
!@+
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      subroutine tracer_ic_ocean
!@sum tracer_ic_ocean initialise ocean tracers
!@auth Gavin Schmidt
!@ver 1.0
      USE MODEL_COM, only: itime
      USE OCN_TRACER_COM, only : itime_tr0, ntm, trname, trw0, n_age
#ifdef TRACERS_SPECIAL_O18
      USE TRACER_COM, only : water_tracer_ic
#endif

      USE OCEAN, only : im,jm,lmo,focean,dxypo,mo,lmm,imaxj,oXYP
#ifdef TRACERS_OCEAN
     *     ,trmo,txmo,tymo,tzmo,s0m,sxmo,symo,szmo,oc_tracer_mean
     *     ,trmo_glob,s0m_glob,mo_glob
#endif
      USE SEAICE, only : xsi,lmi
      USE STRAITS, only : nmst,msist,ssist
#ifdef TRACERS_OCEAN
     *     ,lmst,ist,jst,xst,yst,mmst,s0mst,sxmst,szmst,trmst,txmst
     *     ,tzmst
#endif
#ifdef TRACERS_WATER
     *     ,trsist
#endif
      USE FILEMANAGER, only : openunit,closeunit
      USE DOMAIN_DECOMP_1D, only : get, haveLatitude, GLOBALSUM,
     *     ESMF_BCAST
      USE OCEANR_DIM, only : grid=>ogrid
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT, pack_data, unpack_data
      USE OCEAN, only : scatter_ocean, gather_ocean

      IMPLICIT NONE
      integer n,i,j,l,nst,i1,j1,i2,j2,ll
#ifdef TRACERS_SPECIAL_O18
      integer iu_O18ic,ip1,im1
      character*80 title
      real*4, dimension(im,jm,lmo) :: t0m4,tzm4
      real*8 fwm,afac
#endif
      real*8 t01,t02,trsum,wsum,tratio,frac_tr

      real*8 :: OTRACJ(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) 
      INTEGER :: J_0S, J_1S, J_0, J_1, J_0H, J_1H
      CALL GET(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S, J_STRT = J_0,
     *     J_STOP = J_1, J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)

C**** Note that only sea ice related arrays are initialised if
C**** only TRACERS_WATER is true.

      do n=1,ntm

        if (trname(n).eq.'Age') n_age=n

        if (itime.eq.itime_tr0(n)) then
        select case (trname(n)(1:6))

        case default
#ifdef TRACERS_OCEAN
          trmo(:,:,:,n)=0.
          txmo(:,:,:,n)=0.
          tymo(:,:,:,n)=0.
          tzmo(:,:,:,n)=0.
C**** straits
          if (am_i_root()) then
            trmst(:,:,n)=0.
            txmst(:,:,n)=0.
            tzmst(:,:,n)=0.
          end if
          CALL ESMF_BCAST(grid, trmst)
          CALL ESMF_BCAST(grid, txmst)
          CALL ESMF_BCAST(grid, tzmst)

#endif

#ifdef TRACERS_WATER
          if (am_i_root()) trsist(:,:,n)=0.
          CALL ESMF_BCAST(grid, trsist)
#endif

#if (defined TRACERS_OCEAN) && (defined TRACERS_ZEBRA)
        case ('zebraL')

           read(trname(n)(7:8),'(I2)') ll ! gets level from name
           do l=1,lmo
             do j=J_0,J_1
               do i=1,im
                 if (l.eq.ll .and. l.le.lmm(i,j)) then
                   trmo(i,j,l,n) = mo(i,j,l)*dxypo(j) ! set conc=1 for l=ll
                 else
                   trmo(i,j,l,n) = 0.   ! zero otherwise
                 end if
               end do
             end do
           end do

           if (am_i_root()) then
             do nst=1,nmst
               do l=1,lmst(nst)
                 if (l.eq.ll) then
                   trmst(nst,l,n)= mmst(nst,l)
                 else
                   trmst(nst,l,n)= 0.
                 endif
               end do
             end do
           end if

! set all gradients to zero initially
           txmo(:,:,:,n) = 0; tymo(:,:,:,n)= 0. ; tzmo(:,:,:,n)=0.
           if (am_i_root()) then
             txmst(:,:,n)=0. ; tzmst(:,:,n)=0.
          endif
          CALL ESMF_BCAST(grid, trmst)
          CALL ESMF_BCAST(grid, txmst)
          CALL ESMF_BCAST(grid, tzmst)
#endif

        case ('Water', 'H2O18', 'HDO', 'H2O17' )
#if (defined TRACERS_OCEAN) && (defined TRACERS_SPECIAL_O18)
C**** Open ic file for isotope tracers
          if(water_tracer_ic.eq.1)
     *         call openunit("H2O18ic",iu_O18ic,.true.,.true.)
#endif

#ifdef TRACERS_OCEAN
#ifndef TRACERS_SPECIAL_O18
C**** main ocean variabiles and gradients
          do j=J_0,J_1
          do i=1,im
            do l=1,lmm(i,j)
              trmo(i,j,l,n)=trw0(n)*(mo(i,j,l)*dxypo(j)-s0m(i,j,l))
              txmo(i,j,l,n)=-trw0(n)*sxmo(i,j,l)
              tymo(i,j,l,n)=-trw0(n)*symo(i,j,l)
              tzmo(i,j,l,n)=-trw0(n)*szmo(i,j,l)
            end do
          end do
          end do
#else
C**** read in initial conditions for isotopes
C**** search through for correct tracer (since there is no guarantee
C**** that they will be in same order as tracer indices).
C**** data are now in 'per mil'
          if(water_tracer_ic.eq.1) then
            rewind (iu_O18ic)
 10         read  (iu_O18ic,err=800,end=810) title,t0m4,tzm4
            if (index(title,trim(trname(n))).eq.0) goto 10
            write (6,*) 'Read from H2O18ic: ',title
            call closeunit(iu_O18ic)
          else
            t0m4(:,:,:)=0.
            tzm4(:,:,:)=0.
c            if(n.eq.n_Water) t0m4(:,:,:)=1.
          endif

C**** Turn per mil data into mass ratios (using current standard)
          t0m4(:,:,:)=(t0m4(:,:,:)*1d-3+1.)*trw0(n)
          tzm4(:,:,:)=0.    ! tzm4(:,:,:)*1d-3*trw0(n) corrupted?
C****
          do l=1,lmo
            txmo(:,J_0:J_1,l,n) = 0.
            tymo(:,J_0:J_1,l,n) = 0.
C**** Define East-West horizontal gradients
            im1=im-1
            i=im
            do j=J_0S,J_1S
              do ip1=1,im
                if (lmm(i,j).ge.l) then
                  if (lmm(im1,j).ge.l) then
                    if (lmm(ip1,j).ge.l) then
                      txmo(i,j,l,n)=2.5d-1*(t0m4(ip1,j,l)-t0m4(im1,j,l))
                    else
                      txmo(i,j,l,n)=  5d-1*(t0m4(  i,j,l)-t0m4(im1,j,l))
                    end if
                  else
                    if (lmm(ip1,j).ge.l)
     *                   txmo(i,j,l,n)=5d-1*(t0m4(ip1,j,l)-t0m4(i,j,l))
                  end if
                end if
                im1=i
                i=ip1
              end do
            end do
C**** Define North-South horizontal gradients
            do j=J_0S,J_1S
              do i=1,im
                if (lmm(i,j).ge.l)  then
                  if (lmm(i,j-1).ge.l)  then
                    if (lmm(i,j+1).ge.l)  then
                      tymo(i,j,l,n)=2.5d-1*(t0m4(i,j+1,l)-t0m4(i,j-1,l))
                    else
                      tymo(i,j,l,n)=  5d-1*(t0m4(i,  j,l)-t0m4(i,j-1,l))
                    end if
                  else
                    if (lmm(i,j+1).ge.l)
     *                   tymo(i,j,l,n)=5d-1*(t0m4(i,j+1,l)-t0m4(i,j,l))
                  end if
                end if
              end do
            end do
C**** Multiply ratios by freshwater mass
            do j=J_0,J_1
            do i=1,im
              if (l.le.lmm(i,j)) then
                fwm = mo(i,j,l)*dxypo(j)-s0m(i,j,l)
                trmo(i,j,l,n)=t0m4(i,j,l)*fwm
                txmo(i,j,l,n)=txmo(i,j,l,n)*fwm-sxmo(i,j,l)*t0m4(i,j,l)
                tymo(i,j,l,n)=tymo(i,j,l,n)*fwm-symo(i,j,l)*t0m4(i,j,l)
                tzmo(i,j,l,n)=tzm4(i,j,l)  *fwm-szmo(i,j,l)*t0m4(i,j,l)
              else
                trmo(i,j,l,n)=0.
                txmo(i,j,l,n)=0.
                tymo(i,j,l,n)=0.
                tzmo(i,j,l,n)=0.
              end if
            end do
            end do
          end do
#endif

C**** Initiallise strait values based on adjacent ocean boxes
          call gather_ocean(1)  ! mo,g0m,gx-zmo,s0m,sx-zmo,trmo,tx-zmo

          if(am_I_root()) then
          do nst=1,nmst
            i1=ist(nst,1)
            j1=jst(nst,1)
            i2=ist(nst,2)
            j2=jst(nst,2)
            do l=1,lmst(nst)
              t01=trmo_glob(i1,j1,l,n)/(mo_glob(i1,j1,l)*dxypo(j1)
     *             -s0m_glob(i1,j1,l))
              t02=trmo_glob(i2,j2,l,n)/(mo_glob(i2,j2,l)*dxypo(j2)
     *             -s0m_glob(i2,j2,l))
              trmst(l,nst,n) = 5d-1*(mmst(l,nst)-s0mst(l,nst))*(t01+t02)
              txmst(l,nst,n) = -5d-1*sxmst(l,nst)*(t01+t02)
              tzmst(l,nst,n) = -5d-1*szmst(l,nst)*(t01+t02)
            end do
            do l=lmst(nst)+1,lmo
              trmst(l,nst,n) = 0.
              txmst(l,nst,n) = 0.
              tzmst(l,nst,n) = 0.
            end do
#ifdef TRACERS_WATER
            trsist(n,1:2,nst) = trw0(n)*(msist(1,nst)*xsi(1:2)
     *           -ssist(1:2,nst))
            trsist(n,3:lmi,nst)=trw0(n)*(msist(2,nst)*xsi(3:lmi)
     *           -ssist(3:lmi,nst))
#endif
          end do
          end if

          call bcast_straits(.false.) ! bcst tracers

C**** Balance tracers so that average concentration is TRW0 
C**** or oc_tracer_mean
          
          CALL CONSERV_OTR(OTRACJ,N)
          CALL GLOBALSUM(grid, OTRACJ, trsum, ALL=.true.)

          if (AM_I_ROOT()) then
            if (oc_tracer_mean(n).ne.-999.) then
              tratio=trw0(n)*(oc_tracer_mean(n)*1d-3+1.)
            else
              tratio=trw0(n)
            end if
            
            if (trname(n).eq.'Water') then
              wsum = trsum
              frac_tr = -999.
            else
              frac_tr = tratio/(trsum/wsum)
              write(6,*) "Average oceanic tracer concentration ",
     *             trname(n),(trsum/(wsum*trw0(n))-1d0)*1d3,frac_tr
            end if
          end if
          call ESMF_BCAST(grid,  frac_tr)
         
          if (frac_tr.ne.-999.) then    
            do l=1,lmo
              do j=J_0,J_1
                do i=1,imaxj(j)
                  trmo(i,j,l,n) = trmo(i,j,l,n) * frac_tr
                  if (frac_tr.lt.1.) then
                    TXMO(I,J,L,n)=frac_tr*TXMO(I,J,L,n)
                    TYMO(I,J,L,n)=frac_tr*TYMO(I,J,L,n)
                    TZMO(I,J,L,n)=frac_tr*TZMO(I,J,L,n)
                  end if
                end do
              end do
            end do
C**** straits
            if (am_i_root()) then
              DO NST=1,NMST
                DO L=1,LMST(NST)
                  TRMST(L,NST,N)=frac_tr*TRMST(L,NST,N)
                  if (frac_tr.lt.1) then
                    TXMST(L,NST,N)=frac_tr*TXMST(L,NST,N)
                    TZMST(L,NST,N)=frac_tr*TZMST(L,NST,N)
                  end if
                END DO
              END DO
            end if
            CALL ESMF_BCAST(grid, TRMST)
            CALL ESMF_BCAST(grid, TXMST)
            CALL ESMF_BCAST(grid, TZMST)

C**** Check
            CALL CONSERV_OTR(OTRACJ,N)
            CALL GLOBALSUM(grid, OTRACJ, trsum, ALL=.true.)

            if (AM_I_ROOT()) then
              tratio=(trsum/(wsum*trw0(n))-1.)*1000.
              write(6,*) "New ocean tracer mean: ",tratio
     *             ,oc_tracer_mean(n)
            end if
          end if
#endif
C****
        end select
        write(6,*) trname(n)," tracer initialised in ocean"
        end if
      end do

C**** ensure that atmospheric arrays are properly updated (i.e. gtracer)
      CALL TOC2SST

      return
 800  write(6,*) "Error reading input file H2O18ic"
      call stop_model('stopped in OCN_TRACER.f',255)
 810  write(6,*) "Tracer ",trname(n)," not found in file H2O18ic"
      call stop_model('stopped in OCN_TRACER.f',255)
      end subroutine tracer_ic_ocean
#endif

#ifdef TRACERS_OCEAN
      SUBROUTINE OC_TDECAY(DTS)
!@sum OC_TDECAY decays radioactive tracers in ocean
!@auth Gavin Schmidt/Jean Lerner
      USE MODEL_COM, only : itime
      USE OCN_TRACER_COM, only : ntm,trdecay,itime_tr0
      USE OCEAN, only : trmo,txmo,tymo,tzmo
      IMPLICIT NONE
      real*8, save, dimension(ntm) :: expdec = 1.
      real*8, intent(in) :: dts
      logical, save :: ifirst=.true.
      integer n

      if (ifirst) then               
        do n=1,ntm
          if (trdecay(n).gt.0.0) expdec(n)=exp(-trdecay(n)*dts)
        end do
        ifirst = .false.
      end if

      do n=1,ntm
        if (trdecay(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** Oceanic decay
          trmo(:,:,:,n)   = expdec(n)*trmo(:,:,:,n)
          txmo(:,:,:,n)   = expdec(n)*txmo(:,:,:,n)
          tymo(:,:,:,n)   = expdec(n)*tymo(:,:,:,n)
          tzmo(:,:,:,n)   = expdec(n)*tzmo(:,:,:,n)
        end if
      end do
C****
      return
      end subroutine oc_tdecay

      SUBROUTINE OCN_TR_AGE(DTS)
!@sum OCN_TR_AGE age tracers in ocean
!@auth Gavin Schmidt/Natassa Romanou
      USE CONSTANT, only : sday
      USE MODEL_COM, only : itime,JDperY
      USE OCN_TRACER_COM, only : n_age
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo, imaxj, focean,
     *     lmm, lmo

      USE DOMAIN_DECOMP_1D, only : get
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      real*8, intent(in) :: dts
      real*8 age_inc
      integer i,j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)

C**** at each time step set surface tracer conc=0 and add 1 below
C**** this is mass*age (kg*year)
C**** age=1/(JDperY*24*3600) in years
      age_inc=dts/(JDperY*SDAY)
      DO L=1,LMO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            if (l.le.lmm(i,j)) then
              if (L.eq.1) then
                TRMO(I,J,1,n_age)=0 ; TXMO(I,J,1,n_age)=0 
                TYMO(I,J,1,n_age)=0 ; TZMO(I,J,1,n_age)=0
              else
                TRMO(I,J,L,n_age)= TRMO(I,J,L,n_age) +
     +                age_inc * MO(I,J,L) * oXYP(I,J)
              end if
            end if
          ENDDO
        ENDDO
      ENDDO
C****
      return
      end subroutine ocn_tr_age
#endif

      SUBROUTINE DIAGTCO (M,NT0)
!@sum  DIAGTCO Keeps track of the conservation properties of ocean tracers
!@auth Gary Russell/Gavin Schmidt/Jean Lerner
!@ver  1.0
      USE OCEAN, only : IMO=>IM, oJ_BUDG, oJ_0B, oJ_1B
      USE DIAG_COM, only : jm_budg
      USE TRDIAG_COM, only: tconsrv=>tconsrv_loc,nofmt,title_tcon,
     *     natmtrcons,ntmxcon
      USE DOMAIN_DECOMP_1D, only : GET
      USE OCEANR_DIM, only : oGRID
      IMPLICIT NONE
!@var M index denoting which process changed the tracer
      INTEGER, INTENT(IN) :: m
!@var NT0 index denoting tracer number, NT is index in tconsrv 
      INTEGER, INTENT(IN) :: nt0
      INTEGER :: nt
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(IMO,
     &                  oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) :: TOTAL
      REAL*8, DIMENSION(JM_BUDG) :: TOTALJ
      INTEGER :: nm,ni
      INTEGER :: I, J, J_0, J_1, I_0, I_1

      CALL GET(ogrid, J_STRT=J_0, J_STOP=J_1)
      I_0 = oGRID%I_STRT
      I_1 = oGRID%I_STOP
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1,2...12:  See DIAGCA in DIAG.f
C****   13+ AFTER Sources and Sinks
C****
C**** NOFMT contains the indexes of the TCONSRV array where each
C**** change is to be stored for each quantity. If NOFMT(M,NT)=0,
C**** no calculation is done.
C**** NOFMT(1,NT) is the index for the instantaneous value.
      nt=nt0+natmtrcons
      if (nofmt(m,nt).gt.0) then
C**** Calculate current value TOTAL (kg)
        call conserv_otr(total,nt0)
        nm=nofmt(m,nt)
        ni=nofmt(1,nt)

C**** Calculate zonal sums
        totalj(oj_0b:oj_1b)=0.
        do j=j_0,j_1
        do i=i_0,i_1
          totalj(oj_budg(i,j)) = totalj(oj_budg(i,j)) + total(i,j)
        end do
        end do

c**** Accumulate difference from last time in TCONSRV(NM)
        if (m.gt.1) then
          tconsrv(oJ_0b:oJ_1b,nm,nt) =
     &         tconsrv(oJ_0b:oJ_1b,nm,nt)+(totalj(oJ_0b:oJ_1b)
     *         -tconsrv(oJ_0b:oJ_1b,ni,nt)) 
        end if
C**** Save current value in TCONSRV(NI)
        tconsrv(oJ_0b:oJ_1b,ni,nt)=totalj(oJ_0b:oJ_1b)
      end if
      return
      end subroutine diagtco

      SUBROUTINE conserv_OTR(OTR,ITR)
!@sum  conserv_OTR calculates zonal ocean tracer (kg) on ocean grid
!@auth Gavin Schmidt
!@ver  1.0
      USE OCEAN, only : IMO=>IM,JMO=>JM,oXYP,LMM, imaxj,trmo
      USE STRAITS, only : nmst,jst,ist,lmst,trmst
      USE DOMAIN_DECOMP_1D, only : GET
      USE OCEANR_DIM, only : ogrid

      IMPLICIT NONE
!@var OTR zonal ocean tracer (kg)
      REAL*8, DIMENSION(IMO,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) :: OTR 
      INTEGER, INTENT(IN) :: ITR
      INTEGER I,J,L,N

      INTEGER :: J_0, J_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(ogrid, J_STRT=J_0,    J_STOP=J_1,
     &  HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          OTR(I,J) = SUM(TRMO(I,J,:LMM(I,J),ITR))
        END DO
      END DO
      if (HAVE_SOUTH_POLE) OTR(2:IMO,1)   = OTR(1,1)
      if (HAVE_NORTH_POLE) OTR(2:IMO,JMO) = OTR(1,JMO)
C**** include straits variables
      DO N=1,NMST
        I=IST(N,1)
        J=JST(N,1)
        If (J >= J_0 .and. J <= J_1)
     &       OTR(I,J)=OTR(I,J)+SUM(TRMST(:LMST(N),N,ITR))
      END DO
C****
      RETURN
      END SUBROUTINE conserv_OTR

      SUBROUTINE SET_TCONO(QCON,CONPT,NAME_CON,QSUM,INST_UNIT,SUM_UNIT
     *     ,INST_SC,CHNG_SC, itr0)
!@sum  SET_TCONO assigns ocean conservation diagnostic array indices
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only: sday
      USE MODEL_COM, only: dtsrc,nfiltr
      USE DIAG_COM, only: npts,ia_d5d,ia_d5s,ia_filt,ia_12hr,ia_src
     *     ,conpt0
      USE TRDIAG_COM, only: ktcon,title_tcon,scale_tcon,nsum_tcon
     *     ,nofmt,ia_tcon,name_tconsrv,lname_tconsrv,units_tconsrv
     *     ,ntcons,npts,natmtrcons
      IMPLICIT NONE
!@var QCON denotes at which points conservation diags are saved
      LOGICAL, INTENT(IN),DIMENSION(npts) :: QCON
!@var QSUM sets whether each diag is included in final sum
!@+   should be zero for diags set using DIAGTCB (i.e. using difference)
      LOGICAL, INTENT(IN),DIMENSION(npts) :: QSUM
      LOGICAL, DIMENSION(npts) :: QSUM_CON   ! local version
!@var INST_SC scale for instantaneous value
      REAL*8, INTENT(IN) :: INST_SC
!@var CHNG_SC scale for changes
      REAL*8, INTENT(IN) :: CHNG_SC
!@var NAME_CON name of conservation quantity
      CHARACTER*8, INTENT(IN) :: NAME_CON
!@var sname name of conservation quantity (no spaces)
      CHARACTER*8 :: sname
!@var INST_UNIT string for unit for instant. values
      CHARACTER*20, INTENT(IN) :: INST_UNIT
!@var SUM_UNIT string for unit for summed changes
      CHARACTER*20, INTENT(IN) :: SUM_UNIT
!@var ITR index for the tracer
      INTEGER, INTENT(IN) :: ITR0
!@var CONPT0_sname like CONPT0 but without spaces
      CHARACTER*10, DIMENSION(npts) :: CONPT0_sname, CONPT
      CHARACTER*11 CHGSTR
      CHARACTER*40 clean_str
      INTEGER NI,NM,NS,N,k,itr

C**** make nice netcdf names
      sname=trim(clean_str(name_con))
      do n=1,npts
         conpt0_sname(n) = trim(clean_str(conpt(n)))
      enddo
C****
      NI=1
      itr=itr0+natmtrcons
      NOFMT(1,itr) = NI
      TITLE_TCON(NI,itr) =
     *     "0INSTANT "//TRIM(NAME_CON)//" "//TRIM(INST_UNIT)
      SCALE_TCON(NI,itr) = INST_SC
      name_tconsrv(NI,itr) ="inst_oc_"//sname
      lname_tconsrv(NI,itr) = "INSTANT "//TRIM(NAME_CON)
      units_tconsrv(NI,itr) = INST_UNIT
      NSUM_TCON(NI,itr) = -1
      IA_TCON(NI,itr) = 12
      NM=NI
      DO N=2,npts+1
        IF (QCON(N-1)) THEN
          NM = NM + 1
          NOFMT(N,itr) = NM
          QSUM_CON(NM)=.FALSE.
          IF (QSUM(N-1)) QSUM_CON(NM)=.TRUE.
          CHGSTR=" CHANGE OF "
          if (n.le.npts+1) then
            TITLE_TCON(NM,itr) = CHGSTR//TRIM(NAME_CON)//" BY "//
     *         CONPT(N-1)
            name_tconsrv(NM,itr) =
     *           "chg_oc_"//trim(sname)//"_"//TRIM(CONPT0_sname(N-1))
c          else
c            IF (.not. QSUM(N-1)) CHGSTR="     DELTA "
c            TITLE_TCON(NM,itr) = CHGSTR//TRIM(NAME_CON)//" BY "//
c     *           CONPTs(N-npts-1)
c            name_tconsrv(NM,itr) =
c     *           "chg_"//trim(sname)//"_"//TRIM(CONPTs_sname(N-npts-1))
          end if
          lname_tconsrv(NM,itr) = TITLE_TCON(NM,itr)
          units_tconsrv(NM,itr) = SUM_UNIT
          SELECT CASE (N)
          CASE (2)
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_d5d
          CASE (3,4,5,6,7,9,11,12)
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_d5s
          CASE (8)
            SCALE_TCON(NM,itr) = CHNG_SC/(NFILTR*DTSRC)
            IA_TCON(NM,itr) = ia_filt
          CASE (10)
            SCALE_TCON(NM,itr) = CHNG_SC*2./SDAY
            IA_TCON(NM,itr) = ia_12hr
          CASE (13:)   ! special tracer sources
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_src
          END SELECT
        ELSE
          NOFMT(N,itr) = 0
        END IF
      END DO
      NS=NM+1
      IF (NS.gt.KTCON) THEN
        WRITE(6,*) "KTCON not large enough for extra conserv diags",
     *       KTCON,NI,NM,NS,NAME_CON
        call stop_model(
     &       "Change KTCON in tracer diagnostic common block",255)
      END IF
      DO NM=NI+1,NS-1
        NSUM_TCON(NM,itr) = -1
        IF (QSUM_CON(NM)) NSUM_TCON(NM,itr) = NS
      END DO
      TITLE_TCON(NS,itr) = " SUM OF CHANGES "//TRIM(SUM_UNIT)
      name_Tconsrv(NS,itr) ="sum_chg_oc_"//trim(sname)
      lname_Tconsrv(NS,itr) = " SUM OF CHANGES OF "//TRIM(NAME_CON)
      units_Tconsrv(NS,itr) = SUM_UNIT
      SCALE_TCON(NS,itr) = 1.
      IA_TCON(NS,itr) = 12
      NSUM_TCON(NS,itr) = 0
      RETURN
      END SUBROUTINE set_tcono


