      MODULE TRACER_ADV
!@sum MODULE TRACER_ADV arrays needed for tracer advection

      USE MODEL_COM, ONLY : IM,JM,LM,byim
      SAVE
      INTEGER, PARAMETER :: ncmax=10
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: NSTEPX1, NSTEPX2
      INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: NSTEPZ
      INTEGER NSTEPY(LM,NCMAX), NCYC
C**** zonal mean diags
      REAL*8,  ALLOCATABLE, DIMENSION(:,:)   :: sfbm,sbm,sbf,
     *                                          sfcm,scm,scf
C**** vertically integrated fluxes
      REAL*8,  ALLOCATABLE, DIMENSION(:,:)   :: safv,sbfv

      contains

      SUBROUTINE AADVQ (RM,RMOM,QLIMIT,tname)
!@sum  AADVQ advection driver
!@auth G. Russell, modified by Maxwell Kelley
!@+        Jean Lerner modified this for tracers in mass units
c****
c**** AADVQ advects tracers using the Quadradic Upstream Scheme.
c****
c**** input:
c****  pu,pv,sd (kg/s) = east-west,north-south,vertical mass fluxes
c****      qlimit = whether moment limitations should be used
C****         DT (s) = time step
c****
c**** input/output:
c****     rm = tracer mass
c****   rmom = moments of tracer mass
c****     ma (kg) = fluid mass
c****
      USE MODEL_COM, only : fim
      USE DOMAIN_DECOMP_1D, only : GRID, GET

      USE QUSCOM, ONLY : MFLX,nmom
      USE DYNAMICS, ONLY: pu=>pua, pv=>pva, sd=>sda, mb,ma
      IMPLICIT NONE

      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: rm
      REAL*8, dimension(nmom,im,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: rmom
      logical, intent(in) :: qlimit
      character*8 tname          !tracer name
      integer :: I,J,L,n,nx

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0H, J_1H
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** Fill in values at the poles
      if (HAVE_SOUTH_POLE) then
!$OMP  PARALLEL DO PRIVATE (I,L,N)
      do l=1,lm
         do i=2,im
           rm(i,1 ,l) = rm(1,1 ,l)
           do n=1,nmom
             rmom(n,i,1 ,l) = rmom(n,1,1 ,l)
           enddo
         enddo
      enddo
!$OMP  END PARALLEL DO
      endif
      if (HAVE_NORTH_POLE) then
!$OMP  PARALLEL DO PRIVATE (I,L,N)
        do l=1,lm
          do i=2,im
            rm(i,jm,l) = rm(1,jm,l)
            do n=1,nmom
              rmom(n,i,jm,l) = rmom(n,1,jm,l)
            enddo
          enddo
        enddo
!$OMP  END PARALLEL DO
      endif
C****
C**** Load mass after advection from mass before advection
C****
ccc   ma(:,:,:) = mb(:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MA(:,:,L) = MB(:,:,L)
      ENDDO
!$OMP  END PARALLEL DO
C****
C**** Advect the tracer using the quadratic upstream scheme
C****
C**** loop over cycles
      do n=1,ncyc
ccc   mflx(:,:,:)=pu(:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MFLX(:,:,L) = PU(:,:,L)
      ENDDO
!$OMP  END PARALLEL DO

      call aadvqx (rm,rmom,ma,mflx,qlimit,tname,nstepx1(J_0H,1,n),
     &    safv)

ccc   mflx(:,:,:)=pv(:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MFLX(:,:,L) = PV(:,:,L)
      ENDDO
!$OMP  END PARALLEL DO

      call aadvqy (rm,rmom,ma,mflx,qlimit,tname,nstepy(1,n),
     &    sbf,sbm,sfbm,sbfv)

ccc   mflx(:,:,:)=sd(:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MFLX(:,:,L) = SD(:,:,L)
      ENDDO
!$OMP  END PARALLEL DO

      call aadvqz (rm,rmom,ma,mflx,qlimit,tname,nstepz(1,n),
     &    scf,scm,sfcm)

ccc   mflx(:,:,:)=pu(:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MFLX(:,:,L) = PU(:,:,L)
      ENDDO
!$OMP  END PARALLEL DO

      call aadvqx (rm,rmom,ma,mflx,qlimit,tname,nstepx2(J_0H,1,n),
     *     safv)
      end do

C**** deal with vertical polar box diagnostics outside ncyc loop
      if (HAVE_SOUTH_POLE) then
!$OMP  PARALLEL DO PRIVATE (L)
        do l=1,lm-1
          sfcm(1 ,l) = fim*sfcm(1 ,l)
          scm (1 ,l) = fim*scm (1 ,l)
          scf (1 ,l) = fim*scf (1 ,l)
        end do
!$OMP  END PARALLEL DO
      endif
      if (HAVE_NORTH_POLE) then
!$OMP  PARALLEL DO PRIVATE (L)
        do l=1,lm-1
          sfcm(jm,l) = fim*sfcm(jm,l)
          scm (jm,l) = fim*scm (jm,l)
          scf (jm,l) = fim*scf (jm,l)
        end do
!$OMP  END PARALLEL DO
      endif

      return
      end SUBROUTINE AADVQ

      SUBROUTINE AADVQ0(DT)
!@sum AADVQ0 initialises advection of tracer.
!@+   Decide how many cycles to take such that mass does not become
!@+   too small during any of the operator splitting steps of each cycle
!@auth Maxwell Kelley
c****
C**** The MA array space is temporarily put to use in this section
      USE DYNAMICS, ONLY: mu=>pua, mv=>pva, mw=>sda, mb, ma
      USE DOMAIN_DECOMP_1D, ONLY : GRID, GET, GLOBALSUM, HALO_UPDATE
      USE DOMAIN_DECOMP_1D, ONLY : NORTH, SOUTH, AM_I_ROOT
      USE QUSCOM, ONLY : IM,JM,LM
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: DT
      INTEGER :: i,j,l,n,nc,im1,nbad,nbad_loc
      REAL*8 :: byn,ssp,snp

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0H, J_1H
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

ccc   mu(:,:,:) = mu(:,:,:)*(.5*dt)
ccc   mv(:,J_0:J_1S,:) = mv(:,2:jm,:)*dt
ccc   mv(:,jm,:) = 0.
ccc   mw(:,:,1:lm-1) = mw(:,:,1:lm-1)*(-dt)
C
!$OMP  PARALLEL DO PRIVATE (J,L)
      DO L=1,LM
         IF (HAVE_SOUTH_POLE) MU(:,1,L) = 0.
         DO J=J_0S,J_1S
            MU(:,J,L) = MU(:,J,L)*(.5*DT)
         ENDDO
         IF (HAVE_NORTH_POLE) MU(:,JM,L) = 0.
      ENDDO
!$OMP  END PARALLEL DO
C
      CALL HALO_UPDATE(grid, MV, FROM=NORTH)
!$OMP  PARALLEL DO PRIVATE (J,L,I)
        DO L=1,LM
          DO J=J_0H,J_1S
            DO I=1,IM
              MV(I,J,L) = MV(I,J+1,L)*DT
            END DO
          END DO
          IF (HAVE_NORTH_POLE) MV(:,JM,L) = 0.
        END DO
!$OMP  END PARALLEL DO
C
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM-1
         MW(:,:,L) = MW(:,:,L)*(-DT)
      ENDDO
!$OMP  END PARALLEL DO
C
c     for some reason mu is not zero at the poles...
ccc   mu(:,1,:) = 0.
ccc   mu(:,jm,:) = 0.
C**** Set things up
      nbad = 1
      ncyc = 0
      do while(nbad.gt.0)
      ncyc = ncyc + 1
      byn = 1./ncyc
      nbad_loc = 0
!$OMP  PARALLEL DO PRIVATE (L)
      do l=1,lm
         ma(:,:,l) = mb(:,:,l)
      enddo
!$OMP  END PARALLEL DO
      do nc=1,ncyc

C****     1/2 x-direction
!$OMP  PARALLEL DO PRIVATE (I,IM1,J,L)
!$OMP* REDUCTION(+:NBAD_LOC)
        lloopx1: do l=1,lm
        do j=J_0,J_1
          im1 = im
          do i=1,im
            ma(i,j,l) = ma(i,j,l) + (mu(im1,j,l)-mu(i,j,l))*byn
            if (ma(i,j,l).lt.0.5*mb(i,j,l)) then
               nbad_loc = nbad_loc + 1
c               exit lloopx1 ! saves time in single-processor mode
            endif
            im1 = i
          end do
        end do
        end do lloopx1
!$OMP  END PARALLEL DO
        CALL GLOBALSUM(grid, nbad_loc, nbad, all=.true.)
        IF(NBAD.GT.0) exit ! nc loop
        nbad_loc = nbad

C****         y-direction
!$OMP  PARALLEL DO PRIVATE (I,J,L,SSP,SNP)
!$OMP* REDUCTION(+:NBAD_LOC)
        lloopy: do l=1,lm              !Interior
        do j=J_0S,J_1S
        do i=1,im
          ma(i,j,l) = ma(i,j,l) + (mv(i,j-1,l)-mv(i,j,l))*byn
          if (ma(i,j,l).lt.0.5*mb(i,j,l)) then
             nbad_loc = nbad_loc + 1
c             exit lloopy ! saves time in single-processor mode
          endif
        end do
        end do
        if (HAVE_SOUTH_POLE) then
           ssp = sum(ma(:, 1,l)-mv(:,   1,l)*byn)*byim
           ma(:,1 ,l) = ssp
           if (ma(1,1,l).lt.0.5*mb(1,1,l)) then
              nbad_loc = nbad_loc + 1
c              exit lloopy ! saves time in single-processor mode
           endif
        endif
        if (HAVE_NORTH_POLE) then
           snp = sum(ma(:,jm,l)+mv(:,jm-1,l)*byn)*byim
           ma(:,jm,l) = snp
           if (ma(1,jm,l).lt.0.5*mb(1,jm,l)) then
              nbad_loc = nbad_loc + 1
c              exit lloopy ! saves time in single-processor mode
           endif
        endif
        end do lloopy
!$OMP  END PARALLEL DO
        CALL GLOBALSUM(grid, nbad_loc, nbad, all=.true.)
        IF(NBAD.GT.0) exit ! nc loop
        nbad_loc = nbad
C****         z-direction
!$OMP  PARALLEL DO PRIVATE (I,J,L)
!$OMP* REDUCTION(+:NBAD_LOC)
        lloopz: do l=1,lm
        if(l.eq.1) then ! lowest layer
        do j=J_0,J_1
        do i=1,im
          ma(i,j,l) = ma(i,j,l)-mw(i,j,l)*byn
          if (ma(i,j,l).lt.0.5*mb(i,j,l)) then
             nbad_loc = nbad_loc + 1
c             exit lloopz ! saves time in single-processor mode
          endif
        end do
        end do
        else if(l.eq.lm) then ! topmost layer
        do j=J_0,J_1
        do i=1,im
          ma(i,j,l) = ma(i,j,l)+mw(i,j,l-1)*byn
          if (ma(i,j,l).lt.0.5*mb(i,j,l)) then
             nbad_loc = nbad_loc + 1
c             exit lloopz ! saves time in single-processor mode
          endif
        end do
        end do
        else ! interior layers
        do j=J_0,J_1
        do i=1,im
          ma(i,j,l) = ma(i,j,l)+(mw(i,j,l-1)-mw(i,j,l))*byn
          if (ma(i,j,l).lt.0.5*mb(i,j,l)) then
             nbad_loc = nbad_loc + 1
c             exit lloopz ! saves time in single-processor mode
          endif
        end do
        end do
        endif
        end do lloopz
!$OMP  END PARALLEL DO
        CALL GLOBALSUM(grid, nbad_loc, nbad, all=.true.)
        IF(NBAD.GT.0) exit ! nc loop
        nbad_loc=nbad
C****     1/2 x-direction
!$OMP  PARALLEL DO PRIVATE (I,IM1,J,L)
!$OMP* REDUCTION(+:NBAD_LOC)
        lloopx2: do l=1,lm
        do j=J_0,J_1
          im1 = im
          do i=1,im
            ma(i,j,l) = ma(i,j,l) + (mu(im1,j,l)-mu(i,j,l))*byn
            if (ma(i,j,l).lt.0.5*mb(i,j,l)) then
               nbad_loc = nbad_loc + 1
c               exit lloopx2 ! saves time in single-processor mode
            endif
            im1 = i
          end do
        end do
        end do lloopx2
!$OMP  END PARALLEL DO
        CALL GLOBALSUM(grid, nbad_loc, nbad, all=.true.)
        IF(NBAD.GT.0) exit ! nc loop
        nbad_loc=nbad

      end do ! nc loop

      if(ncyc.ge.10) then
         if (AM_I_ROOT()) then
            write(6,*) 'stop: ncyc=10 in AADVQ0'
            call stop_model('AADVQ0: ncyc>=10',11)
         end if
      end if
      enddo ! while(nbad.gt.0)
      if(ncyc.gt.2) then
         if (AM_I_ROOT()) write(6,*) 'AADVQ0: ncyc>2',ncyc
      end if
C**** Divide the mass fluxes by the number of cycles
      byn = 1./ncyc
ccc   mu(:,:,:)=mu(:,:,:)*byn
ccc   mv(:,J_0:J_1S,:)=mv(:,J_0:J_1S,:)*byn
ccc   mw(:,:,:)=mw(:,:,:)*byn
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         mu(:,:,l)=mu(:,:,l)*byn
      ENDDO
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         mv(:,J_0:J_1S,l)=mv(:,J_0:J_1S,l)*byn
      ENDDO
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         mw(:,:,l)=mw(:,:,l)*byn
      ENDDO
!$OMP  END PARALLEL DO
C****
C**** Decide how many timesteps to take by computing Courant limits
C****
ccc   MA(:,:,:) = MB(:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MA(:,:,L) = MB(:,:,L)
      ENDDO
!$OMP  END PARALLEL DO
      do n=1,ncyc
        call xstep (MA,nstepx1(J_0H,1,n))
        call ystep (MA,nstepy(1,n))
        call zstep (MA,nstepz(1,n))
        call xstep (MA,nstepx2(J_0H,1,n))
      end do
      RETURN
  900 format (1x,a,3i4,f10.4,i5)
  910 format (1x,a,i4,f10.4,i5)
      END subroutine AADVQ0

      end MODULE TRACER_ADV


      subroutine aadvQx(rm,rmom,mass,mu,qlimit,tname,nstep,
     *     safv)
!@sum  AADVQX advection driver for x-direction
!@auth Maxwell Kelley; modified by J. Lerner
c****
c**** aadvtQ advects tracers in the west to east direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mu (kg) = west-east mass flux, positive eastward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      use QUSDEF
      USE DOMAIN_DECOMP_1D, only : GRID, GET, GLOBALSUM
ccc   use QUSCOM, only : im,jm,lm, xstride,am,f_i,fmom_i
      use QUSCOM, only : im,jm,lm, xstride
      implicit none
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) ::
     &                                         rm,mass,mu
      REAL*8, dimension(nmom,im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) ::
     &                                         rmom
      logical ::  qlimit
      REAL*8  AM(IM), F_I(IM), FMOM_I(NMOM,IM)
      REAL*8, intent(inout),
     *        dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: safv
      character*8 tname
      integer :: nstep(GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm)
      integer :: i,j,l,ierr,nerr,ns,ICKERR,ICKERR_LOC

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      REAL*8 :: safvl(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

c**** loop over layers and latitudes
      ICKERR_LOC=0
      safvl = 0
!$OMP  PARALLEL DO PRIVATE (J,L,NS,AM,F_I,FMOM_I,IERR,NERR)
!$OMP* SHARED(IM,QLIMIT,XSTRIDE)
!$OMP* REDUCTION(+:ICKERR_LOC)
      do l=1,lm
      do j=J_0S,J_1S
      am(:) = mu(:,j,l)/nstep(j,l)
c****
c**** call 1-d advection routine
c****
      do ns=1,nstep(j,l)
      call adv1d(rm(1,j,l),rmom(1,1,j,l), f_i,fmom_i, mass(1,j,l),
     &        am, im, qlimit,xstride,xdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvQx: i,j,l=",nerr,j,l,' ',tname
        if (ierr.eq.2) write(6,*) "Error in qlimit: abs(a) > 1"
        if (ierr.eq.2) ICKERR_LOC=ICKERR_LOC+1
      end if

! store tracer flux in safv array
      safvl(:,j,l) = safvl(:,j,l) + f_i(:) 

      enddo ! ns
      enddo ! j
      enddo ! l
!$OMP  END PARALLEL DO
      safv = safv + sum(safvl,3)

      If (HAVE_NORTH_POLE) safv(:,jm) = 0. ! no horizontal flux at poles
      If (HAVE_SOUTH_POLE) safv(:,1) = 0.  
C
      CALL GLOBALSUM(grid, ICKERR_LOC, ICKERR, all=.true.)
      IF(ICKERR.GT.0)  CALL stop_model('Stopped in aadvQx',11)
C
      return
c****
      end subroutine aadvQx


      subroutine aadvQy(rm,rmom,mass,mv,qlimit,tname,nstep,
     &   sbf,sbm,sfbm, sbfv)
!@sum  AADVQY advection driver for y-direction
!@auth Maxwell Kelley; modified by J. Lerner
c****
c**** aadvQy advects tracers in the south to north direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mv (kg) = north-south mass flux, positive northward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      USE DOMAIN_DECOMP_1D, only : GRID, GET, GLOBALSUM
      use DOMAIN_DECOMP_1D, only : AM_I_ROOT
      USE DOMAIN_DECOMP_1D, ONLY : TRANSP, TRANSPOSE_COLUMN
      use CONSTANT, only : teeny
      use QUSDEF
ccc   use QUSCOM, only : im,jm,lm, ystride,bm,f_j,fmom_j, byim
      use QUSCOM, only : im,jm,lm, ystride,               byim
      implicit none
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: 
     &                                         rm,mass,mv
      REAL*8, dimension(nmom,im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: 
     &                                         rmom
      logical ::  qlimit
      REAL*8, intent(inout), 
     &        dimension(GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: 
     &                                         sfbm,sbm,sbf
      REAL*8, intent(inout),
     *        dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: sbfv
      character*8 tname
      integer :: i,j,l,ierr,ns,nstep(lm),ICKERR, ICKERR_LOC
      integer :: err_loc(3)
      REAL*8, DIMENSION(LM) :: m_sp,m_np,rm_sp,rm_np,rzm_sp,rzm_np,
     &     rzzm_sp,rzzm_np
      REAL*8, DIMENSION(LM) :: dm_sp,dm_np,drm_sp,drm_np,drzm_sp,
     &     drzm_np,drzzm_sp, drzzm_np

      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: fqv
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: 
     &     F_J,BM
      REAL*8  FMOM_J(NMOM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)

      INTEGER :: nm

c****Get relevant local distributed parameters
      INTEGER J_0, J_1
      INTEGER J_0H, J_1H
      INTEGER J_0S, J_1S
      LOGICAL :: idx(LM)
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1, 
     *               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               J_STRT_SKP=J_0S,  J_STOP_SKP=J_1S,
     *               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     *               HAVE_NORTH_POLE=HAVE_NORTH_POLE)

c**** loop over layers
      ICKERR_LOC=0
c**** loop over timesteps
        do ns=1,maxval(nstep)

!$OMP  PARALLEL DO PRIVATE (I,J,L)
!$OMP* SHARED(JM,QLIMIT,YSTRIDE,ns)
!$OMP* REDUCTION(+:ICKERR_LOC)
      do l=1,lm
        if (ns == 1)  fqv(:,:,l)=0
        if (ns > nstep(l)) cycle

c**** scale polar boxes to their full extent
          If (HAVE_SOUTH_POLE) THEN
            mass(:,1,l)=mass(:,1,l)*im
            m_sp(l) = mass(1,1 ,l)
            rm(:,1,l)=rm(:,1,l)*im
            rm_sp(l) = rm(1,1 ,l)
            do i=1,im
              rmom(:,i,1 ,l)=rmom(:,i,1 ,l)*im
            enddo
            rzm_sp(l)  = rmom(mz ,1,1 ,l)
            rzzm_sp(l) = rmom(mzz,1,1 ,l)
          End If

          If (HAVE_NORTH_POLE) THEN
            mass(:,jm,l)=mass(:,jm,l)*im
            m_np(l) = mass(1,jm,l)
            rm(:,jm,l)=rm(:,jm,l)*im
            rm_np(l) = rm(1,jm,l)
            do i=1,im
              rmom(:,i,jm,l)=rmom(:,i,jm,l)*im
            enddo
            rzm_np(l)  = rmom(mz ,1,jm,l)
            rzzm_np(l) = rmom(mzz,1,jm,l)
          End IF

c****
c**** load 1-dimensional arrays
c****
          bm (:,:,l) = mv(:,:,l)/nstep(l)
          If (HAVE_SOUTH_POLE) rmom(ihmoms,:,1,l)=0
          If (HAVE_NORTH_POLE) THEN
            bm(:,jm,l) = 0.
            rmom(ihmoms,:,jm,l) = 0.
          End IF
      enddo  ! end loop over levels
!$OMP  END PARALLEL DO

c****
c**** call 1-d advection routine
c****
      idx =(ns <= nstep)
        call advection_1D_custom( rm(1,j_0h,1), rmom(1,1,j_0h,1), 
     &       f_j(1,j_0h,1),fmom_j(1,1,j_0h,1), mass(1,j_0h,1),
     &       bm(1,j_0h,1),j_1-j_0+1,LM,idx,
     &     qlimit,ystride,ydir,ierr,err_loc)

      if (ierr.gt.0) then
        write(6,*) "Error in aadvQy: i,j,l=",err_loc,' ',tname
        if (ierr.eq.2) write(6,*) "Error in qlimit: abs(b) > 1"
        if (ierr.eq.2) ICKERR_LOC=ICKERR_LOC+1
      end if

!$OMP  PARALLEL DO PRIVATE (I,J,L)
!$OMP* SHARED(JM,QLIMIT,YSTRIDE,ns)
!$OMP* REDUCTION(+:ICKERR_LOC)
      do l = 1,LM
        if (ns > nstep(l)) cycle

! horizontal moments are zero at pole
        IF (HAVE_SOUTH_POLE) rmom(ihmoms,:,1, l) = 0
        IF (HAVE_NORTH_POLE) rmom(ihmoms,:,jm,l) = 0.
c     sbfijl(i,:,l) = sbfijl(i,:,l)+f_j(:)

      fqv(:,j_0:j_1,l) = fqv(:,j_0:j_1,l) + f_j(:,j_0:j_1,l)  !store tracer flux in fqv array
      If (HAVE_NORTH_POLE) fqv(:,jm,l) = 0.       ! play it safe


c**** average and unscale polar boxes
      if (HAVE_SOUTH_POLE) then
        mass(:,1 ,l) = (m_sp(l) + sum(mass(:,1 ,l)-m_sp(l)))*byim
        rm(:,1 ,l) = (rm_sp(l) + sum(rm(:,1 ,l)-rm_sp(l)))*byim
        rmom(mz ,:,1 ,l) =
     &       (rzm_sp(l) + sum(rmom(mz ,:,1 ,l)-rzm_sp(l) ))*byim
        rmom(mzz,:,1 ,l) =
     &       (rzzm_sp(l)+ sum(rmom(mzz,:,1 ,l)-rzzm_sp(l)))*byim
      end if   !SOUTH POLE

      if (HAVE_NORTH_POLE) then
        mass(:,jm,l) = (m_np(l) + sum(mass(:,jm,l)-m_np(l)))*byim
        rm(:,jm,l) = (rm_np(l) + sum(rm(:,jm,l)-rm_np(l)))*byim
        rmom(mz ,:,jm,l) =
     &       (rzm_np(l) + sum(rmom(mz ,:,jm,l)-rzm_np(l) ))*byim
        rmom(mzz,:,jm,l) =
     &       (rzzm_np(l)+ sum(rmom(mzz,:,jm,l)-rzzm_np(l)))*byim
      end if  !NORTH POLE

      enddo  ! end loop over levels
!$OMP  END PARALLEL DO

      end do ! time step

      do l=1,lm
        do j=J_0,J_1S           ! zonal mean diagnostics
          sfbm(j,l) = sfbm(j,l) + sum(fqv(:,j,l)/(mv(:,j,l)+teeny))
          sbm (j,l) = sbm (j,l) + sum(mv(:,j,l))
          sbf (j,l) = sbf (j,l) + sum(fqv(:,j,l))
        end do
      end do

      do l=1,lm                 ! vert. integrated diagnostics
        do j=J_0,J_1S
            sbfv (:,j) = sbfv (:,j) + fqv(:,j,l)
        end do
      end do
C
      CALL GLOBALSUM(grid, ICKERR_LOC, ICKERR, all=.true.)
      IF(ICKERR.NE.0)  call stop_model('Stopped in aadvQy',11)
C

      return
c****
      end subroutine aadvQy


      subroutine aadvQz(rm,rmom,mass,mw,qlimit,tname,nstep,
     &  scf,scm,sfcm)
!@sum  AADVQZ advection driver for z-direction
!@auth Maxwell Kelley; modified by J. Lerner
c****
c**** aadvQz advects tracers in the upward vertical direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mw (kg) = vertical mass flux, positive upward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      use CONSTANT, only : teeny
      use GEOM, only : imaxj
      USE DOMAIN_DECOMP_1D, only : GRID, GET
      use QUSDEF
ccc   use QUSCOM, only : im,jm,lm, zstride,cm,f_l,fmom_l
      use QUSCOM, only : im,jm,lm, zstride
      implicit none
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) ::
     &                                         rm,mass,mw
      REAL*8, dimension(nmom,im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) ::
     &                                         rmom
      INTEGER, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: nstep
      REAL*8, intent(inout),
     &               dimension(GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) ::
     &                                         sfcm,scm,scf
      logical ::  qlimit
      REAL*8, dimension(lm) :: fqw
      character*8 tname
      REAL*8  CM(LM),F_L(LM),FMOM_L(NMOM,LM)
      integer :: i,j,l,ierr,nerr,ns,ICKERR,ICKERR_LOC

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

c**** loop over latitudes and longitudes
      ICKERR_LOC=0.
!$OMP  PARALLEL DO PRIVATE (I,J,L,NS,CM,F_L,FMOM_L,FQW,IERR,NERR)
!$OMP* SHARED(LM,QLIMIT,ZSTRIDE)
!$OMP* REDUCTION(+:ICKERR_LOC)
      do j=J_0,J_1
      do i=1,imaxj(j)
      fqw(:) = 0.
      cm(:) = mw(i,j,:)/nstep(i,j)
      cm(lm)= 0.
c****
c**** call 1-d advection routine
c****
      do ns=1,nstep(i,j)
      call adv1d(rm(i,j,1),rmom(1,i,j,1),f_l,fmom_l,mass(i,j,1),
     &        cm,lm,qlimit,zstride,zdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvQz: i,j,l=",i,j,nerr,' ',tname
        if (ierr.eq.2) write(6,*) "Error in qlimit: abs(c) > 1"
        if (ierr.eq.2) ICKERR_LOC=ICKERR_LOC+1
      end if
      fqw(:)  = fqw(:) + f_l(:) !store tracer flux in fqw array
      enddo ! ns
      do l=1,lm-1   !diagnostics
        sfcm(j,l) = sfcm(j,l) + fqw(l)/(mw(i,j,l)+teeny)
        scm (j,l) = scm (j,l) + mw(i,j,l)
        scf (j,l) = scf (j,l) + fqw(l)
      enddo
      enddo ! i
      if (j.eq.1.or.j.eq.jm) then
        do l=1,lm
          do i=2,im
            rm(i,j,l)=rm(1,j,l)
            rmom(:,i,j,l)=rmom(:,1,j,l)
            mass(i,j,l)=mass(1,j,l)
          end do
        end do
      end if
      enddo ! j
!$OMP  END PARALLEL DO
C
      IF(ICKERR_LOC.GT.0)  call stop_model('Stopped in aadvQz',11)
C
      return
c****
      end subroutine aadvQz



      SUBROUTINE XSTEP (M,NSTEPX)
!@sum XSTEP determines the number of X timesteps for tracer dynamics
!@+    using Courant limits
!@auth J. Lerner and M. Kelley
!@ver  1.0
      USE DOMAIN_DECOMP_1D, ONLY : GRID, GET, GLOBALSUM
      USE QUSCOM, ONLY : IM,JM,LM,byim
      USE DYNAMICS, ONLY: mu=>pua
      IMPLICIT NONE
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: m
      REAL*8, dimension(im) :: a,am,mi
      integer, dimension(GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: nstepx
      integer :: l,j,i,ip1,im1,nstep,ns,ICKERR,ICKERR_LOC
      REAL*8 :: courmax

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** Decide how many timesteps to take by computing Courant limits
C
      ICKERR_LOC = 0
!$OMP  PARALLEL DO PRIVATE (I,IP1,IM1,J,L,NSTEP,NS,COURMAX,A,AM,MI)
!$OMP* REDUCTION(+:ICKERR_LOC)
      DO 420 L=1,LM
      DO 420 J=J_0S,J_1S
      nstep=0
      courmax = 2.
      do while(courmax.gt.1.)
        nstep = nstep+1   !(1+int(courmax))
        am(:) = mu(:,j,l)/nstep
        mi(:) = m (:,j,l)
        courmax = 0.
        do ns=1,nstep
          i = im
          do ip1=1,im
            if(am(i).gt.0.) then
               a(i) = am(i)/mi(i)
               courmax = max(courmax,+a(i))
            else
               a(i) = am(i)/mi(ip1)
               courmax = max(courmax,-a(i))
            endif
          i = ip1
          enddo  ! ip1=1,im
C**** Update air mass
          im1 = im
          do i=1,im
            mi(i) = mi(i)+am(im1)-am(i)
          im1 = i
          enddo
        enddo    ! ns=1,nstep
        if(nstep.ge.20) write(6,*) 'aadvqx: nstep.ge.20'
        if(nstep.ge.20)  then
           write(6,*) 'aadvqx: j,l,nstep,courmax=',j,l,nstep,courmax
           courmax=-1.
           ICKERR_LOC=ICKERR_LOC+1
        end if
      enddo      ! while(courmax.gt.1.)
C**** Correct air mass
      M(:,J,L) = MI(:)
      NSTEPX(J,L) = NSTEP
c     if(nstep.gt.2 .and. nx.eq.1)
c    *  write(6,'(a,3i3,f7.4)')
c    *  'aadvqx: j,l,nstep,courmax=',j,l,nstep,courmax
  420 CONTINUE
!$OMP  END PARALLEL DO
C
      ! Ensure all processes agree on stop criteria
      CALL GLOBALSUM(grid, ICKERR_LOC, ICKERR, all=.true.)
      IF(ICKERR.GT.0)  call stop_model('Stopped in XSTEP',11)
C
      RETURN
      END SUBROUTINE XSTEP


      SUBROUTINE YSTEP (M,NSTEPY)
!@sum YSTEP determines the number of Y timesteps for tracer dynamics
!@+    using Courant limits
!@auth J. Lerner and M. Kelley
!@ver  1.0
      USE DOMAIN_DECOMP_1D, ONLY : GRID, GET, HALO_UPDATE, NORTH, SOUTH
      USE DOMAIN_DECOMP_1D, ONLY : GLOBALSUM, GLOBALMAX
      USE QUSCOM, ONLY : IM,JM,LM,byim
      USE DYNAMICS, ONLY: mv=>pva
      IMPLICIT NONE
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: m
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: mij
      REAL*8, dimension(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: b,bm
      integer, dimension(LM) :: nstepy
      integer :: jprob,iprob,nstep,ns,i,j,l,ICKERR,ICKERR_LOC
      REAL*8 :: courmax,courmax_loc, byn,sbms,sbmn

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** decide how many timesteps to take (all longitudes at this level)
      ICKERR_LOC=0
      CALL HALO_UPDATE(grid, mv, FROM=SOUTH)
      CALL HALO_UPDATE(grid,  m, FROM=NORTH)
c$$$!$OMP  PARALLEL DO PRIVATE (I,J,L,NS,NSTEP,COURMAX,BYN,B,BM,MIJ,
c$$$!$OMP*          COURMAX_LOC,IPROB,JPROB,SBMS,SBMN)
c$$$!$OMP* REDUCTION(+:ICKERR_LOC)
      DO 440 L=1,LM
C**** Scale poles
      if (HAVE_SOUTH_POLE) m(:, 1,l) =   m(:, 1,l)*im !!!!! temporary
      if (HAVE_NORTH_POLE) m(:,jm,l) =   m(:,jm,l)*im !!!!! temporary
C**** begin computation
      nstep=0
      courmax = 2.
      do while(courmax.gt.1.)
        nstep = nstep+1   !(1+int(courmax_loc))
        byn = 1./nstep
        courmax_loc = 0.
        mij(:,:) = m(:,:,l)
        do ns=1,nstep
          do j=J_0,J_1S
          do i=1,im
            bm(j) = mv(i,j,l)*byn
            if(bm(j).gt.0.) then
               b(j) = bm(j)/mij(i,j)
               courmax_loc = max(courmax_loc,+b(j))
               if (courmax_loc.eq.b(j)) jprob = j
               if (courmax_loc.eq.b(j)) iprob = i
            else
               b(j) = bm(j)/mij(i,j+1)
               courmax_loc = max(courmax_loc,-b(j))
               if (courmax_loc.eq.-b(j)) jprob = j
               if (courmax_loc.eq.-b(j)) iprob = i
            endif
          enddo
          enddo
          ! Get courmax across all PEs
          CALL GLOBALMAX(grid, courmax_loc, courmax)
C**** Update air mass at poles
          if (HAVE_SOUTH_POLE) then
            sbms = sum(mv(:,   1,l))*byn
            mij(:, 1) = mij(:, 1)-sbms
          endif
          if (HAVE_NORTH_POLE) then
            sbmn = sum(mv(:,jm-1,l))*byn
            mij(:,jm) = mij(:,jm)+sbmn
          endif
C**** Update air mass in the interior
          do j=J_0S,J_1S
            mij(:,j) = mij(:,j)+(mv(:,j-1,l)-mv(:,j,l))*byn
          enddo
        enddo    ! ns=1,nstep
        if(nstep.ge.20) then
          If (courmax == courmax_loc) THEN ! I have the worst case
           write(6,*) 'courmax=',courmax,l,iprob,jprob
           write(6,*) 'aadvqy: nstep.ge.20'
          endif
          ICKERR_LOC=ICKERR_LOC+1
          courmax = -1.
        endif
      enddo      ! while(courmax.gt.1.)
C**** Correct air mass
      m(:,:,l) = mij(:,:)
      NSTEPY(L) = nstep
c       if(nstep.gt.1. and. nTRACER.eq.1) write(6,'(a,2i3,f7.4)')
c    *    'aadvqy: l,nstep,courmax=',l,nstep,courmax
C**** Unscale poles
      if (HAVE_SOUTH_POLE) m(:, 1,l) =   m(:, 1,l)*byim !!! undo temporary
      if (HAVE_NORTH_POLE) m(:,jm,l) =   m(:,jm,l)*byim !!! undo temporary
  440 CONTINUE
c$$$!$OMP  END PARALLEL DO
C
      CALL GLOBALSUM(grid, ICKERR_LOC, ICKERR, all=.true.)
      IF(ICKERR.GT.0)  call stop_model('Stopped in YSTEP',11)
C
      RETURN
      END SUBROUTINE YSTEP


      SUBROUTINE ZSTEP (M,NSTEPZ)
!@sum ZSTEP determines the number of Z timesteps for tracer dynamics
!@+    using Courant limits
!@auth J. Lerner and M. Kelley
!@ver  1.0
      USE DOMAIN_DECOMP_1D, ONLY : GRID, GET, GLOBALSUM
      USE QUSCOM, ONLY : IM,JM,LM,byim
      USE DYNAMICS, ONLY: mw=>sda
      IMPLICIT NONE
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: m
      REAL*8, dimension(lm) :: ml
      REAL*8, dimension(0:lm) :: c,cm
      integer, dimension(im*(GRID%J_STOP_HALO-GRID%J_STRT_HALO+1)) :: 
     &                                                         nstepz
      integer :: nstep,ns,l,i,j,ICKERR,ICKERR_LOC
      REAL*8 :: courmax,byn

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** decide how many timesteps to take
      ICKERR_LOC=0
!$OMP  PARALLEL DO PRIVATE (I,J,L,NS,NSTEP,COURMAX,BYN,C,CM,ML)
!$OMP* REDUCTION(+:ICKERR_LOC)
      DO J=J_0,J_1
      DO I=1,IM
      nstep=0
      courmax = 2.
      do while(courmax.gt.1.)
        nstep = nstep+1   !(1+int(courmax))
        byn = 1./nstep
        cm(1:lm) = mw(i,j,1:lm)*byn
        ml(:)  = m(i,j,:)
        CM(LM)= 0. ! VERY IMPORTANT TO SET THIS TO ZERO
        CM( 0)= 0. ! VERY IMPORTANT TO SET THIS TO ZERO
        courmax = 0.
        do ns=1,nstep
          do l=1,lm-1
            if(cm(l).gt.0.) then
               c(l) = cm(l)/ml(l)
               courmax = max(courmax,+c(l))
            else
               c(l) = cm(l)/ml(l+1)
               courmax = max(courmax,-c(l))
            endif
          enddo
          do l=1,lm
            ml(l) = ml(l)+(cm(l-1)-cm(l))
          enddo
        enddo    ! ns=1,nstep
        if(nstep.ge.20) write(6,*) 'aadvqz: nstep.ge.20'
        if(nstep.ge.20)  then
           write(6,*)  'aadvqz: nstep.ge.20'
           ICKERR_LOC=ICKERR_LOC+1
           courmax = -1.
        end if
      enddo      ! while(courmax.gt.1.)
C**** Correct air mass
      m(i,j,:) = ml(:)
cgsfc      NSTEPZ(I+IM*(J-1)) = NSTEP
      NSTEPZ(I+IM*(J-J_0H)) = NSTEP
c     if(nstep.gt.1 .and. nTRACER.eq.1) write(6,'(a,2i7,f7.4)')
c    *   'aadvqz: i,j,nstep,courmax=',i,j,nstep,courmax
      END DO
      END DO
!$OMP  END PARALLEL DO
C
      CALL GLOBALSUM(grid, ICKERR_LOC, ICKERR, all=.true.)
      IF(ICKERR.GT.0)  call stop_model('Stopped in ZSTEP',11)
C
      RETURN
      END SUBROUTINE ZSTEP

      SUBROUTINE ALLOC_TRACER_ADV(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE TRACER_ADV
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID, GET
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE(  NSTEPX1(J_0H:J_1H,LM,NCMAX),
     *           NSTEPX2(J_0H:J_1H,LM,NCMAX),
     *            NSTEPZ(IM*(J_1H-J_0H+1),NCMAX) )
  
      ALLOCATE( sfbm(J_0H:J_1H,LM),
     *           sbm(J_0H:J_1H,LM),
     *           sbf(J_0H:J_1H,LM),
     *          sfcm(J_0H:J_1H,LM),
     *           scm(J_0H:J_1H,LM),
     *           scf(J_0H:J_1H,LM) )

      ALLOCATE( safv(IM,J_0H:J_1H),
     *          sbfv(IM,J_0H:J_1H) )

      END SUBROUTINE ALLOC_TRACER_ADV


