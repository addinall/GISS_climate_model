      program agcstat
!@sum agcstat calculates AGC output fields defined with respect
!@+   to time-averaged quantities (standing eddies etc.)
!@auth M. Kelley
      implicit none
      real*4, dimension(:), allocatable :: lat_dg,plm,ple,lons,dxv,dxyp
     &     ,fcor,cosv
      real*4, dimension(:,:), allocatable :: xjl,lats_dg
      character(len=80) :: run_name,time_period,gcfile,ijfile,ijkfile
      character(len=30) :: units
      character(len=80) :: lname,title
      character(len=80) :: vname
      character(len=20) :: latname
      integer :: nargs,j,k,l,im,jm,jm_,lm,kgz
      real*4, parameter :: missing=-1.e30
      include 'netcdf.inc'
      integer :: gc_fid,ij_fid,ijk_fid,vid,dids(7)
      integer :: status,varid,ij_nvars
      integer :: plm_id,ple_id, lat_id,lat2_id
      real*4, dimension(:,:,:), allocatable :: phi
      real*4, dimension(:,:), allocatable :: ampltd,phase
      real*4, dimension(:), allocatable :: sinkx,coskx
      character(len=1) :: strk
      real*4 :: ak,bk,ampfac
      real*4 :: pi,radian,grav,radius
      real*4, dimension(:,:,:), allocatable :: dpuv,u3d,v3d,dse
      real*4, dimension(:,:), allocatable :: seke,seuv,sedv,dpa,
     &     temp,theta,u,n_sqr,dalphadp,dx,dqdy,dqdybyu,nbyf_sqr,
     &     by2h_sqr,refr_ind,nptsavg,pmid
      real*4, dimension(:), allocatable :: vort
      real*4 :: psum,ubar,vbar,dbar
      real*4 :: tf,rgas,teeny,omega,dlat,
     &     bydy,dtheta,dlnthetadp,dlnthetadlnp,hscale,
     &     dp,dp2,dxavg,MBYRCOS_sqr
      real*4 :: zeros_hemis(10000)
      integer :: ldn,lup
      INTEGER, PARAMETER, DIMENSION(5) :: MW=(/1,2,3,6,9/)
      logical :: have_aijk

      zeros_hemis = 0.
      omega = 7.292e-5
      tf = 273.16
      rgas = 287.05
      teeny = 1.e-20
      grav = 9.81
      pi = acos(-1.)
      radian = pi/180.
      radius = 6371000.

      nargs = iargc()
      if(nargs.ne.2) then
        write(6,*) 'usage: agcstat run_name time_period'
        write(6,*) '  e.g. agcstat Exyz JAN1901'
        write(6,*) '       agcstat Exyz ANN1901-1910'
        stop
      endif

      call getarg(1,run_name)
      call getarg(2,time_period)

      gcfile=trim(time_period)//'.agc'//trim(run_name)//'.nc'
      ijfile=trim(time_period)//'.aij'//trim(run_name)//'.nc'
      ijkfile=trim(time_period)//'.aijk'//trim(run_name)//'.nc'

      write(6,*) 'Searching for the following input files:'
      write(6,*) '   ',trim(gcfile)
      write(6,*) '   ',trim(ijfile)
      write(6,*) '   ',trim(ijkfile),' (optional)'

c
c open input files
c
      call handle_err(nf_open(gcfile,nf_write,gc_fid),
     &     'opening '//trim(gcfile))
      call handle_err(nf_open(ijfile,nf_nowrite,ij_fid),
     &     'opening '//trim(ijfile))
      status = nf_open(ijkfile,nf_nowrite,ijk_fid)
      have_aijk = status==nf_noerr
      if(.not.have_aijk) write(6,*)
     &     'no aijk-file: skipping standing eddy calculations'

      status = nf_inq_nvars(ij_fid,ij_nvars)

      latname = 'lat'
      call get_dimsize(gc_fid,'lat',jm)
      call get_dimsize(gc_fid,'plm',lm)
      call get_dimsize(gc_fid,'pgz',kgz)

      call get_dimsize(ij_fid,'lon',im)
      call get_dimsize(ij_fid,'lat',jm_)
      if(jm.ne.jm_) stop 'mismatched jmlat in agc and aij'

      status = nf_inq_dimid(gc_fid,'lat',lat_id)
      lat2_id = -1
      status = nf_inq_dimid(gc_fid,'lat2',lat2_id)
      status = nf_inq_dimid(gc_fid,'plm',plm_id)
      status = nf_inq_dimid(gc_fid,'ple',ple_id)

c
c allocate workspace
c
      allocate(lat_dg(jm),lats_dg(jm,2),xjl(jm,lm))
      allocate(plm(lm),ple(lm),lons(im))
      allocate(dxv(jm),dxyp(jm),fcor(jm),cosv(jm))

c
c read geometry
c
      dlat = pi/jm
      call get_var_real(gc_fid,'lat',lats_dg(1,1))
      if(lat2_id.ge.0)
     &     call get_var_real(gc_fid,'lat2',lats_dg(1,2))
      call get_var_real(gc_fid,'plm',plm)
      call get_var_real(gc_fid,'ple',ple)
      call get_var_real(ij_fid,'lon',lons)
      cosv = cos(lats_dg(:,2)*radian)
      dxv = 2.*pi*radius*cosv
      dxyp(1) = 2.*pi*radius*radius*(sin(radian*lats_dg(2,2))+1.)
      do j=2,jm-1
        dxyp(j) = 2.*pi*radius*radius*
     &       (sin(radian*lats_dg(j+1,2))-sin(radian*lats_dg(j,2)))
        fcor(j) = -(cos(radian*lats_dg(j+1,2))
     &             -cos(radian*lats_dg(j,2)))*2.*omega/dlat
      enddo
      dxyp(jm) = dxyp(1)
      fcor(jm) = cos(radian*lats_dg(jm,2))*2.*omega/dlat
      fcor(1) = -fcor(jm)

c
c amplitude and phase of geopotential heights.  reading would be
c simpler if 3D data were output as a 3D array rather than as a
c collection of 2D arrays
c
      allocate(phi(im,jm,kgz),ampltd(jm,kgz),phase(jm,kgz))
      allocate(sinkx(im),coskx(im))
      l = 0
      do varid=1,ij_nvars
        vname = ''
        status = nf_inq_varname(ij_fid,varid,vname)
        if(vname(1:4).ne.'phi_') cycle
        if(index(vname,'_hemis').gt.0) cycle
        l = l + 1
        if(l.gt.kgz) exit
        status = nf_get_var_real(ij_fid,varid,phi(1,1,l))
      enddo
      if(l.ne.kgz) then
        write(6,*) 'Differing numbers of geopot. levels in '//
     &     'agc and aij'
        stop
      endif
      ampfac = 2./real(im)
      do k=1,4
        sinkx = sin(k*lons*radian)
        coskx = cos(k*lons*radian)
        do l=1,kgz
          do j=1,jm
            ak = sum(coskx*phi(:,j,l))
            bk = sum(sinkx*phi(:,j,l))
            ampltd(j,l) = sqrt(ak*ak+bk*bk)*ampfac
            phase(j,l) = -atan2(bk,ak)/radian
          enddo
        enddo
        write(strk,'(i1)') k
        vname = 'phi_amp_wave'//strk
        call put_var_real(gc_fid,trim(vname),ampltd)
        call put_var_real(gc_fid,trim(vname)//'_hemis',zeros_hemis)
        vname = 'phi_phase_wave'//strk
        call put_var_real(gc_fid,trim(vname),phase)
        call put_var_real(gc_fid,trim(vname)//'_hemis',zeros_hemis)
      enddo

c
c standing eddies
c
      if(have_aijk) then
      allocate(dpuv(im,jm,lm),u3d(im,jm,lm),v3d(im,jm,lm),dse(im,jm,lm))
      allocate(seke(jm,lm),seuv(jm,lm),sedv(jm,lm))
      call get_var_real(ijk_fid,'dpb',dpuv)
      call get_var_real(ijk_fid,'ub',u3d)
      call get_var_real(ijk_fid,'vb',v3d)
      call get_var_real(ijk_fid,'dse',dse)
      where(dpuv.eq.0.)
        u3d = 0.
        v3d = 0.
        dse = 0.
      end where
      do l=1,lm
        do j=1,jm
          psum = sum(dpuv(:,j,l))+teeny
          ubar = sum(dpuv(:,j,l)*u3d(:,j,l))/psum
          vbar = sum(dpuv(:,j,l)*v3d(:,j,l))/psum
          dbar = sum(dpuv(:,j,l)*dse(:,j,l))/psum
          seke(j,l) = .5*(
     &         sum(dpuv(:,j,l)*(u3d(:,j,l)**2+v3d(:,j,l)**2))/psum
     &         -ubar*ubar -vbar*vbar)
          seuv(j,l) = sum(dpuv(:,j,l)*(u3d(:,j,l)*v3d(:,j,l)))/psum
     &         -ubar*vbar
          sedv(j,l) = (sum(dpuv(:,j,l)*(dse(:,j,l)*v3d(:,j,l)))/psum
     &         -dbar*vbar)*dxv(j)*(100./grav)
        enddo
      enddo
      call put_var_real(gc_fid,'stand_eddy_ke',seke)
      call put_var_real(gc_fid,'stand_eddy_ke_hemis',zeros_hemis)
      call put_var_real(gc_fid,'nt_u_stand_eddy',seuv)
      call put_var_real(gc_fid,'nt_u_stand_eddy_hemis',zeros_hemis)
      call put_var_real(gc_fid,'nt_dse_stand_eddy',sedv)
      call put_var_real(gc_fid,'nt_dse_stand_eddy_hemis',zeros_hemis)
      deallocate(dpuv,u3d,v3d,dse)
      deallocate(seke,seuv,sedv)
      endif

C****
C**** D/DY OF Q-G POTENTIAL VORTICITY AND REFRACTION INDICES
C****

      allocate(dpa(jm,lm),temp(jm,lm),theta(jm,lm),u(jm,lm),
     &     n_sqr(jm,lm),dalphadp(jm,lm),dx(jm,lm),dqdy(jm,lm),
     &     dqdybyu(jm,lm),nbyf_sqr(jm,lm),refr_ind(jm,lm),
     &     by2h_sqr(jm,lm),nptsavg(jm,lm),pmid(jm,lm))
      allocate(vort(jm))

      call get_var_real(gc_fid,'dp_cp1',dpa)
      call get_var_real(gc_fid,'npts_avg1',nptsavg)
      call get_var_real(gc_fid,'temp',temp)
      call get_var_real(gc_fid,'pot_temp',theta)
      call get_var_real(gc_fid,'u',u)

      ! shift u to secondary latitudes if necessary
      status = nf_inq_varid(gc_fid,'u',vid)
      status = nf_inq_vardimid(gc_fid,vid,dids)
      if(dids(1).eq.lat_id) then
        where(dpa.eq.0.) u=0.
        do l=1,lm
          do j=jm,2,-1
            u(j,l) = .5*(u(j-1,l)+u(j,l))
          enddo
        enddo
      endif

      do l=1,lm
        pmid(:,l) = ple(l) + .5*dpa(:,l)*im/(nptsavg(:,l)+teeny)
      enddo

      where(dpa.eq.0.)
        temp=0.
        theta=0.
      end where
      temp = temp + tf

      bydy = jm/(radius*pi)

      dalphadp = 0.
      dx = 0.
      N_Sqr = 0.

      do j=1,jm
        do l=lm,1,-1
          if(dpa(j,l).lt.teeny) exit
          lup=min(l+1,lm)
          ldn=max(l-1,1)
          if(dpa(j,ldn).lt.teeny) ldn=l
          dtheta=theta(j,lup)-theta(j,ldn)
          if(abs(dtheta).lt.teeny) dtheta=sign(teeny,dtheta)
          dp=pmid(j,ldn)-pmid(j,lup)+teeny
          dlnthetadp = dtheta/(theta(j,l)*dp)
          dlnthetadlnp = pmid(j,l)*dlnthetadp
          hscale = rgas*temp(j,l)/grav
          dx(j,l) = -fcor(j)*pmid(j,l)/(temp(j,l)*dlnthetadp)
          dalphadp(j,l) = -(temp(j,lup)/pmid(j,lup)
     &                     -temp(j,ldn)/pmid(j,ldn))/dp
          n_sqr(j,l) = dlnthetadlnp*grav/hscale
        enddo
      enddo

      do l=1,lm
        vort(1)=-u(2,l)*dxv(2)/dxyp(1) +fcor(1)
        do j=2,jm-1
          vort(j)=(u(j,l)*dxv(j)-u(j+1,l)*dxv(j+1))/dxyp(j)+fcor(j)
        enddo
        vort(jm)=u(jm,l)*dxv(jm)/dxyp(jm) +fcor(jm)
        do j=jm,2,-1
          dp2=dpa(j,l)+dpa(j-1,l)+teeny
          by2h_sqr(j,l)=.25*(grav/(rgas*.5*(temp(j-1,l)+temp(j,l))))**2
          nbyf_sqr(j,l) =
     &         ((n_sqr(j-1,l)*dpa(j-1,l)+n_sqr(j,l)*dpa(j,l))/dp2) *
     &         (2./(fcor(j-1)*fcor(j-1)+fcor(j)*fcor(j)))
          ubar=u(j,l)
          if(abs(ubar).lt.teeny) ubar=sign(teeny,ubar)
          dxavg = (dx(j,l)*dpa(j,l)+dx(j-1,l)*dpa(j-1,l))/dp2
          dqdy(j,l)=bydy*(vort(j)-vort(j-1)
     &         +(dalphadp(j,l)-dalphadp(j-1,l))*dxavg )
          dqdybyu(j,l)=dqdy(j,l)/ubar
        enddo
      enddo

      call put_var_real(gc_fid,'del_qgpv',dqdy)
      call put_var_real(gc_fid,'del_qgpv_hemis',zeros_hemis)

      refr_ind = 0.
      do k=1,5
        do j=2,jm
          MBYRCOS_sqr = (MW(K)/(RADIUS*COSV(J)))**2
          do l=1,lm
            refr_ind(j,l) =
     &           Nbyf_Sqr(J,L)*(DqDybyU(J,L)-MBYRCOS_sqr)-By2h_Sqr(J,L)
          enddo
        enddo
        refr_ind = refr_ind * 1.e8
        write(strk,'(i1)') mw(k)
        vname = 'refr_ind_wave'//strk
        call put_var_real(gc_fid,trim(vname),refr_ind)
        call put_var_real(gc_fid,trim(vname)//'_hemis',zeros_hemis)
      enddo

c
c deallocate workspace
c
      deallocate(lat_dg,lats_dg,xjl)
      deallocate(plm,ple,lons)
      deallocate(dxv,dxyp,fcor,cosv)

      deallocate(phi,ampltd,phase)
      deallocate(sinkx,coskx)

      deallocate(dpa,temp,theta,u,
     &     n_sqr,dalphadp,dx,dqdy,
     &     dqdybyu,nbyf_sqr,refr_ind,
     &     by2h_sqr,nptsavg,pmid)
      deallocate(vort)

c
c close files
c
      status = nf_close(gc_fid)
      status = nf_close(ij_fid)

      end program agcstat
