      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer, intent(in), optional :: jdim
      logical, intent(in), optional :: no_xdim
      integer :: jd,nd,jm,vid,rc
      integer, dimension(7) :: srt,cnt
      real*8, dimension(:), allocatable :: arr1d
#ifndef SERIAL_MODE
      integer :: nl,nk,nsend,ierr,src_pe,tag,modpe
      integer :: mpi_status(mpi_status_size)
      logical :: has_xdim
#endif
      if(present(jdim)) then
        jd=jdim
      else
        jd=2
      endif
      nd = size(shape(arr))
      srt(:) = 1
      cnt(1:nd) = shape(arr)
      jm = grid%npy
      if(grid%am_i_globalroot) then
        rc = nf_inq_varid(fid,trim(varname),vid)
        if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &       trim(varname),' not found in output file - stopping'
      endif
      call stoprc(rc,nf_noerr)
#ifndef SERIAL_MODE
      tag=1024
      if(present(no_xdim)) then
        has_xdim = .not.no_xdim
      elseif(jd.eq.1) then
        has_xdim = .false.
      else
        has_xdim = .true.
      endif
c if latlon model, call routines better suited for very thin
c zonal strips if they are available
      if(grid%ntiles.eq.1 .and. has_xdim .and. jd.ne.4 .and.
     &     (jd.eq.2 .or. jd.lt.nd)) then
        call par_write_jdecomp_optimized(grid,fid,vid,arr,jd)
        return
      endif
      if(grid%am_i_rowroot) then
        call get_nlnk(arr,jd,nl,nk)
        nsend = nl*nk*grid%maxnj
        if(has_xdim) then
          nsend = nsend*grid%npx
        elseif(jd.gt.1) then
          nsend = nsend*size(arr,jd-1)
        endif
        allocate(arr1d(nsend))
        if(has_xdim) cnt(jd-1) = grid%npx
        cnt(nd+1) = 1
      endif
      if(has_xdim) then
        call pack_row(grid,arr,arr1d,jdim=jd)
      else
        call pack_row_no_xdim(grid,arr,arr1d,jd)
      endif
      if(grid%am_i_globalroot) then
        do src_pe=0,grid%nproc-grid%nprocx,grid%nprocx
          if(src_pe.gt.0) then
            call mpi_recv(arr1d,nsend,MPI_DOUBLE_PRECISION,src_pe,
     &           tag,MPI_COMM_WORLD,mpi_status,ierr)
          endif
          modpe=mod(src_pe,grid%nproc_tile)+1
          srt(jd) = grid%jsr(modpe)
          cnt(jd) = 1+grid%jer(modpe)-grid%jsr(modpe)
          srt(nd+1) = 1 + src_pe/grid%nproc_tile
          rc = nf_put_vara_double(fid,vid,srt,cnt,arr1d)
        enddo
      elseif(grid%am_i_rowroot) then
        call mpi_send(arr1d,nsend,MPI_DOUBLE_PRECISION,0,
     &       tag,MPI_COMM_WORLD,ierr)
      endif
      if(grid%am_i_rowroot) deallocate(arr1d)
      call mpi_barrier(MPI_COMM_WORLD, ierr)
#else
      if(size(arr,jd).eq.jm .or. size(arr,jd).eq.1) then
        if(size(arr,jd).eq.1) then      ! SCM hack: initialize global
          allocate(arr1d(size(arr)*jm)) ! domain to zero
          arr1d = 0.
          rc = nf_put_var_double(fid,vid,arr1d)
          deallocate(arr1d)
        endif
        srt(jd-1) = grid%i_strt
        srt(jd)   = grid%j_strt
        rc = nf_put_vara_double(fid,vid,srt,cnt,arr)
      else
        call stop_model('unrecognized jdim size in do_par_write_nc',255)
      endif
#endif
