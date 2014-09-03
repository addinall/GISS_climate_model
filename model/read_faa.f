      subroutine read_daily_faa
      use resolution, only : lm ! just for temporary declaration of emiss_5d
      use resolution, only : im,jm ! just for debug write
      use model_com, only : JMON,JDATE,JYEAR
      use domain_decomp_atm, only : get,grid,am_i_root
      use geom, only : lonlat_to_ij
      use domain_decomp_atm, only : pack_data ! just for debug write
      use faa_com
      use tracer_com, only : daily_z   !altitude of model levels (m), updated once per day
!@+   and used for vertical distribution of 3D emissions
      implicit none
c
      character(len=5), dimension(nqty) :: qtynames
      real*8, dimension(:,:), allocatable :: emiss_2d_sum_root
      character(len=80) :: title_debug
c
      include 'netcdf.inc'
      integer :: fid,did,vid,status
      character(len=80) :: infile
      character(len=2) :: mon_str,day_str
      character(len=4) :: year_str
      integer :: ij(2),i,j,l,lz,n,hr,ncells,qty
      real*4, dimension(:,:), allocatable :: lons,lats,alts,emiss1d
      real*8 :: ll(2)
      integer, dimension(:,:), allocatable :: iofn,jofn,lofn,i2ofn
      integer :: i_0,i_1,j_0,j_1
      call get(grid, i_strt=i_0,i_stop=i_1, j_strt=j_0,j_stop=j_1)

      qtynames(1) = 'ifuel' 
      qtynames(2) = 'iCO' 
      qtynames(3) = 'iHC'       
      qtynames(4) = 'iNOx' 
      qtynames(5) = 'ipmnv' 
      qtynames(6) = 'ipmfo' 
c
c determine the input filename
c
      write(mon_str,'(i2.2)') jmon
      write(day_str,'(i2.2)') jdate
      write(year_str,'(i4)') jyear
      infile='FAA_DIR/faa_emissions_'//mon_str//'_'//
     &     day_str//'_'//'2006.nc'
      if(am_i_root()) write(6,*) 'opening FAA file ',trim(infile)
c
c open the input file and determine model i,j,l for each data point
c
      call handle_err(nf_open(trim(infile),nf_nowrite,fid),
     &     'problem opening '//trim(infile))
      call handle_err(nf_inq_dimid(fid,'ncells',did),
     &     'ncells dim missing in '//trim(infile))
      status = nf_inq_dimlen(fid,did,ncells)
      allocate(lons(ncells,24),lats(ncells,24),
     &         alts(ncells,24),emiss1d(ncells,24))
      allocate(iofn(ncells,24),jofn(ncells,24),i2ofn(ncells,24),
     &         lofn(ncells,24))
      call get_var_real(fid,'lon',lons)
      call get_var_real(fid,'lat',lats)
      call get_var_real(fid,'alt',alts)
      do hr=1,24
        do n=1,ncells
          lofn(n,hr) = -999
          if(alts(n,hr).lt.0.) cycle ! missing data
          ll = (/ lons(n,hr)-.1, lats(n,hr) /)
          call lonlat_to_ij(ll,ij)
          if(ij(2) < j_0 .or. ij(2) > j_1 .or.
     &       ij(1) < i_0 .or. ij(1) > i_1) then
            cycle ! data is not in the domain of this proessor
          endif
          iofn(n,hr) = ij(1)
          jofn(n,hr) = ij(2)
          ll(1) = lons(n,hr)+.1d0
          call lonlat_to_ij(ll,ij)
          i2ofn(n,hr) = ij(1)
c determine which model level corresponds to the emission altitude
C is there a less expensive way to do this?
          loop_l: do lz=1,LM
          if(alts(n,hr) <= daily_z(iofn(n,hr),jofn(n,hr),lz)) then
          lofn(n,hr) = lz
          exit loop_l
          endif
          enddo loop_l
        enddo
      enddo
c
c read emissions and place into model i,j,l array
c
      emiss_5d = 0.
      do qty=1,nqty
        call get_var_real(fid,trim(qtynames(qty)),emiss1d)
        do hr=1,24
          do n=1,ncells
            l = lofn(n,hr)
            if(l < 1) cycle ! missing data or not on this PE
            j = jofn(n,hr)
            i = iofn(n,hr)
            emiss_5d(i,j,l,hr,qty) = emiss_5d(i,j,l,hr,qty)
     &           + .5*emiss1d(n,hr)
            i = i2ofn(n,hr)
            emiss_5d(i,j,l,hr,qty) = emiss_5d(i,j,l,hr,qty)
     &           + .5*emiss1d(n,hr)
          enddo
        enddo
      enddo

c
c close the input file
c
      status = nf_close(fid)
c
c deallocate workspace
c
      deallocate(lons,lats,alts,emiss1d)
      deallocate(iofn,jofn,lofn,i2ofn)
c
c debug write to fort.743
c
c      emiss_2d_sum = sum(sum(emiss_5d(:,:,:,:,1),4),3)
c      if(am_i_root()) allocate(emiss_2d_sum_root(im,jm))
c      call pack_data(grid,emiss_2d_sum,emiss_2d_sum_root)
c      if(am_i_root()) then
c        title_debug = trim(qtynames(1))//' vertical+daily sum'
c        write(743) title_debug,real(emiss_2d_sum_root,kind=4)
c        deallocate(emiss_2d_sum_root)
c      endif
c      return
      contains
      subroutine handle_err(rc,errstr)
      integer :: rc
      character(len=*) :: errstr
      if(rc.ne.nf_noerr) call stop_model(errstr,255)
      end subroutine handle_err

      subroutine get_var_real(fid,var_name,var)
      implicit none
      integer :: fid
      character(len=*) :: var_name
      real*4 :: var(1)
      integer :: rc,varid
      rc = nf_inq_varid(fid,var_name,varid)
      if(rc.ne.nf_noerr) then
        write(6,*) 'nonexistent variable ',trim(var_name)
        stop
      else
        rc = nf_get_var_real(fid,varid,var)
      endif
      return
      end subroutine get_var_real
      end subroutine read_daily_faa
