      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: varname
      integer, intent(in), optional :: record
      integer :: vid,did,rc,nd
      integer*8, dimension(7) :: srt,cnt
      integer*8 :: nrecs8
      rc = nfmpi_inq_varid(fid,trim(varname),vid)
      if(grid%am_i_globalroot .and. rc.ne.nf_noerr) then
        write(6,*) 'variable ',
     &       trim(varname),' not found in output file - stopping'
      endif
      call stoprc(rc,nf_noerr)
      if(present(record)) then
        nrecs8 = 0
        rc = nfmpi_inq_unlimdim(fid,did)
        rc = nfmpi_inq_dimlen(fid,did,nrecs8)
        if(record.le.0 .or. nrecs8+1.lt.record) then
          if(grid%am_i_globalroot) write(6,*)
     &         'error in record dim spec. for variable ',trim(varname)
          call stoprc(0,1)
        endif
      endif
      rc = nfmpi_begin_indep_data(fid)
      if(grid%am_i_globalroot) then
        if(present(record)) then
          nd = size(shape(arr))
          if(nd.gt.0) then
            cnt(1:nd) = shape(arr)
            srt(1:nd) = 1
          endif
          cnt(nd+1) = 1
          srt(nd+1) = record
          rc = nfmpi_put_vara_double(fid,vid,srt,cnt,arr)
        else
          rc = nfmpi_put_var_double(fid,vid,arr)
        endif
      endif
      rc = nfmpi_end_indep_data(fid)
