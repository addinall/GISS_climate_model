      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      logical, intent(in), optional :: bcast_all
      integer :: vid,rc,ierr,nd,arrsize,shp(7)
      logical :: bcast_all_
      bcast_all_ = .false.
      if(present(bcast_all)) bcast_all_ = bcast_all
      rc = nfmpi_inq_varid(fid,trim(varname),vid)
      if(rc.ne.nf_noerr) then
        if(grid%am_i_globalroot) write(6,*) 'variable ',
     &       trim(varname),' not found in input file - skipping'
        return
      endif
      if(bcast_all_) then
        rc = nfmpi_get_var_double_all(fid,vid,arr)
      else
        rc = nfmpi_begin_indep_data(fid)
        if(grid%am_i_globalroot) rc = nfmpi_get_var_double(fid,vid,arr)
        rc = nfmpi_end_indep_data(fid)
      endif
