      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer, intent(in), optional :: jdim
#ifndef SERIAL_MODE
      real*8, dimension(:), allocatable :: arr_nohalo
      integer :: nl,nk,jd,nd,vid,rc,arrsize
      integer*8, dimension(7) :: srt,cnt
      rc = nfmpi_inq_varid(fid,trim(varname),vid)
      if(rc.ne.nf_noerr) then
        if(grid%am_i_globalroot) write(6,*) 'variable ',
     &       trim(varname),' not found in input file - skipping'
	return
      endif
      if(present(jdim)) then
        jd=jdim
      else
        jd=2
      endif
      nd = size(shape(arr))
      cnt(1:nd) = shape(arr)
      srt(1:nd) = 1
      srt(jd) = grid%js
      srt(jd-1:jd) = (/ grid%is, grid%js /)
      cnt(jd-1:jd) = 1+(/ grid%ie-grid%is, grid%je-grid%js /)
      srt(nd+1) = grid%tile
      cnt(nd+1) = 1
      if(grid%have_domain) then
        call get_nlnk(arr,jd,nl,nk)
      else
        srt=1; cnt=0   ! empty domain
      endif
      arrsize = product(cnt(1:nd))
      allocate(arr_nohalo(arrsize))
      rc = nfmpi_get_vara_double_all(fid,vid,srt,cnt,arr_nohalo)
      if(arrsize.gt.0) then
        call copy_from_1D( arr_nohalo,arr, nl,nk,
     &       grid%isd,grid%ied,grid%jsd,grid%jed,
     &       grid%is,grid%ie,grid%js,grid%je )
      endif
      deallocate(arr_nohalo)
#else
c call read_data when in serial mode
      call read_data(grid,fid,varname,arr)
#endif
