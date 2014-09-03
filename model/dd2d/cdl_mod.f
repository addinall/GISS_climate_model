      module cdl_mod
      implicit none
      private

      integer, parameter :: cdl_strlen=100
      integer, parameter :: ndims_max=100
      integer, parameter :: ndata_max=1000
      integer, parameter :: nlines_max=20000
      character(len=4) :: indent='    '
      character(len=8) :: indent2='        '

      public :: cdl_type
      type cdl_type
        character(len=30) :: name ! name of this cdl for error reporting
        integer :: nlines     ! current number of lines
        integer :: ndims,ncoordlines,nvarlines,ndatalines
        character(len=cdl_strlen), dimension(nlines_max) :: text
        character(len=cdl_strlen), dimension(ndims_max) :: dims
        character(len=cdl_strlen), dimension(3*ndims_max) :: coords
        character(len=cdl_strlen), dimension(nlines_max) :: vars
        character(len=cdl_strlen), dimension(ndata_max) :: datavalues
      end type cdl_type

      public :: cdl_strlen,init_cdl_type,print_cdl,add_dim,add_coord,
     &     add_var,add_varline,add_dataline,add_vardata,add_unlimdim,
     &     copy_dims,copy_coord_vars,
     &     copy_data,merge_cdl,assemble_cdl,defvar_cdl,write_cdl
      

      interface add_vardata
        module procedure add_vardata_r8_1d
        module procedure add_vardata_int_1d
        module procedure add_vardata_1d_array_of_strings
      end interface add_vardata

      contains

      subroutine init_cdl_type(cdl_name,cdl)
      character(len=*), intent(in) :: cdl_name
      type(cdl_type), intent(inout) :: cdl
      cdl%name = cdl_name
      cdl%ndims = 0
      cdl%ncoordlines = 0
      cdl%nvarlines = 0
      cdl%ndatalines = 0
      return
      end subroutine init_cdl_type

      subroutine add_dim(cdl,dimname,dimsize)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: dimname
      integer, intent(in) :: dimsize
      character(len=10) :: dimstr
      integer :: k
      write(dimstr,'(i10)') dimsize
      k = cdl%ndims + 1
      cdl%dims(k) = indent//trim(dimname)//' = '//
     &     trim(adjustl(dimstr))//' ;'
      cdl%ndims  = k
      return
      end subroutine add_dim

      subroutine add_unlimdim(cdl,dimname)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: dimname
      integer :: k
      k = cdl%ndims + 1
      cdl%dims(k) = indent//trim(dimname)//
     &     ' = UNLIMITED ; // (0 currently)'
      cdl%ndims  = k
      return
      end subroutine add_unlimdim

      subroutine copy_dims(cdl_in,cdl_out)
      type(cdl_type), intent(in) :: cdl_in
      type(cdl_type), intent(inout) :: cdl_out
      integer :: k,kout
      kout = cdl_out%ndims
      do k=1,cdl_in%ndims
        kout = kout + 1
        cdl_out%dims(kout) = cdl_in%dims(k)
      enddo
      cdl_out%ndims = kout
      return
      end subroutine copy_dims

      subroutine add_coord(cdl,coordname,coordsize,
     &     units,long_name,coordvalues)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: coordname
      character(len=*), intent(in), optional :: units,long_name
      integer, intent(in) :: coordsize
      real*8, dimension(coordsize), intent(in), optional :: coordvalues
      integer :: k
      call add_dim(cdl,coordname,coordsize)
      k = cdl%ncoordlines + 1
      cdl%coords(k) = indent//'float '//trim(coordname)//'('//
     &     trim(coordname)//') ;'
      if(present(units)) then
        k = k + 1
        cdl%coords(k) = indent2//trim(coordname)//':units = "'//
     &     trim(units)//'" ;'
      endif
      if(present(long_name)) then
        k = k + 1
        cdl%coords(k) = indent2//trim(coordname)//':long_name = "'//
     &     trim(long_name)//'" ;'
      endif
      cdl%ncoordlines  = k
      if(present(coordvalues)) then
        call add_vardata(cdl,trim(coordname),coordvalues)
      endif
      return
      end subroutine add_coord

      subroutine copy_coord_vars(cdl_in,cdl_out)
      type(cdl_type), intent(in) :: cdl_in
      type(cdl_type), intent(inout) :: cdl_out
      integer :: k,kout
      kout = cdl_out%ncoordlines
      do k=1,cdl_in%ncoordlines
        kout = kout + 1
        cdl_out%coords(kout) = cdl_in%coords(k)
      enddo
      cdl_out%ncoordlines = kout
      return
      end subroutine copy_coord_vars

      subroutine copy_data(cdl_in,cdl_out)
      type(cdl_type), intent(in) :: cdl_in
      type(cdl_type), intent(inout) :: cdl_out
      integer :: k,kout
      kout = cdl_out%ndatalines
      do k=1,cdl_in%ndatalines
        kout = kout + 1
        cdl_out%datavalues(kout) = cdl_in%datavalues(k)
      enddo
      cdl_out%ndatalines = kout
      return
      end subroutine copy_data

      subroutine add_var(cdl,varstr,units,long_name)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varstr
      character(len=*), intent(in), optional :: units,long_name
      character(len=cdl_strlen) :: varname,tmpstr
      integer :: k,n1,n2
      k = cdl%nvarlines + 1
      tmpstr = adjustl(varstr)
      cdl%vars(k) = indent//tmpstr
      n1 = index(tmpstr,' ')+1
      n2 = index(tmpstr,'(')-1
      varname = tmpstr(n1:n2)
      if(present(units)) then
        if(len_trim(units).gt.0) then
          k = k + 1
          cdl%vars(k) = indent2//trim(varname)//':units = "'//
     &         trim(units)//'" ;'
        endif
      endif
      if(present(long_name)) then
        k = k + 1
        cdl%vars(k) = indent2//trim(varname)//':long_name = "'//
     &     trim(long_name)//'" ;'
      endif
      cdl%nvarlines  = k
      return
      end subroutine add_var

      subroutine add_varline(cdl,varstr)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varstr
      integer :: k
      k = cdl%nvarlines + 1
      cdl%vars(k) = indent2//trim(varstr)
      cdl%nvarlines  = k
      return
      end subroutine add_varline

      subroutine add_dataline(cdl,varstr)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varstr
      integer :: k
      k = cdl%ndatalines + 1
      cdl%datavalues(k) = indent//trim(varstr)
      cdl%ndatalines  = k
      return
      end subroutine add_dataline

      subroutine add_vardata_r8_1d(cdl,varname,values,fmtstr)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varname
      real*8, dimension(:), intent(in) :: values
      character(len=*), intent(in), optional :: fmtstr
      integer :: k,line,i1,i2,pos,varsize
      character(len=32) :: fmtstr_
      integer, parameter :: npl=6 ! 6 values per line
      if(present(fmtstr)) then
        fmtstr_ = '(6('//trim(fmtstr)//',","))'
      else
        fmtstr_ = '(6(1pe13.5,","))'
      endif
      varsize = size(values)
      k = cdl%ndatalines + 1
      cdl%datavalues(k) = indent//trim(varname)//' = '
      do line=1,(varsize+npl-1)/npl
        i1 = 1 + npl*(line-1)
        i2 = min(varsize,i1+npl-1)
        k = k + 1
        cdl%datavalues(k) = ''
        write(cdl%datavalues(k),fmtstr_) values(i1:i2)
        if(i2.eq.varsize) then
          pos = len_trim(cdl%datavalues(k))
          cdl%datavalues(k)(pos:pos) = ';'
        endif
      enddo
      cdl%ndatalines  = k
      return
      end subroutine add_vardata_r8_1d

      subroutine add_vardata_int_1d(cdl,varname,values)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varname
      integer, dimension(:), intent(in) :: values
      integer :: k,line,i1,i2,pos,varsize
      integer, parameter :: npl=8 ! 8 values per line
      varsize = size(values)
      k = cdl%ndatalines + 1
      cdl%datavalues(k) = indent//trim(varname)//' = '
      do line=1,(varsize+npl-1)/npl
        i1 = 1 + npl*(line-1)
        i2 = min(varsize,i1+npl-1)
        k = k + 1
        cdl%datavalues(k) = ''
        write(cdl%datavalues(k),'(8(i10,","))') values(i1:i2)
        if(i2.eq.varsize) then
          pos = len_trim(cdl%datavalues(k))
          cdl%datavalues(k)(pos:pos) = ';'
        endif
      enddo
      cdl%ndatalines  = k
      return
      end subroutine add_vardata_int_1d

      subroutine add_vardata_1d_array_of_strings(cdl,varname,strings)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varname
      character(len=*), dimension(:), intent(in) :: strings(:)
      integer :: k,line,varsize
      character*1 punct
      varsize = size(strings)
      k = cdl%ndatalines + 1
      cdl%datavalues(k) = indent//trim(varname)//' = '
      punct = ','
      do line=1,varsize
        k = k + 1
        if(line.eq.varsize) punct=';'
        cdl%datavalues(k) = ''
        write(cdl%datavalues(k),'(a)')
     &       '"'//strings(line)//'"'//punct
      enddo
      cdl%ndatalines  = k
      return
      end subroutine add_vardata_1d_array_of_strings

      subroutine assemble_cdl(cdl)
      type(cdl_type), intent(inout) :: cdl
      integer :: k,n
      cdl%text = ''
      cdl%text(1) = 'netcdf xxx { '
      cdl%text(2) = 'dimensions:  '
      k = 2
      do n=1,cdl%ndims
        k = k + 1
        cdl%text(k) = cdl%dims(n)
      enddo
      k = k + 1
      cdl%text(k) = 'variables:  '
      do n=1,cdl%ncoordlines
        k = k + 1
        cdl%text(k) = cdl%coords(n)
      enddo
      do n=1,cdl%nvarlines
        k = k + 1
        cdl%text(k) = cdl%vars(n)
      enddo
      if(cdl%ndatalines.gt.0) then
        k = k + 1
        cdl%text(k) = 'data:  '
        do n=1,cdl%ndatalines
          k = k + 1
          cdl%text(k) = cdl%datavalues(n)
        enddo
      endif
      k = k + 1
      cdl%text(k) = '}'
      cdl%nlines = k
      return
      end subroutine assemble_cdl

      subroutine merge_cdl(cdl1,cdl2,cdl)
      type(cdl_type), intent(in) :: cdl1,cdl2
      type(cdl_type), intent(inout) :: cdl
      integer :: k,n
c dims
      k = 0
      do n=1,cdl1%ndims
        k = k + 1
        cdl%dims(k) = cdl1%dims(n)
      enddo
      do n=1,cdl2%ndims
        k = k + 1
        cdl%dims(k) = cdl2%dims(n)
      enddo
      cdl%ndims = k
c coords
      k = 0
      do n=1,cdl1%ncoordlines
        k = k + 1
        cdl%coords(k) = cdl1%coords(n)
      enddo
      do n=1,cdl2%ncoordlines
        k = k + 1
        cdl%coords(k) = cdl2%coords(n)
      enddo
      cdl%ncoordlines = k
c vars
      k = 0
      do n=1,cdl1%nvarlines
        k = k + 1
        cdl%vars(k) = cdl1%vars(n)
      enddo
      do n=1,cdl2%nvarlines
        k = k + 1
        cdl%vars(k) = cdl2%vars(n)
      enddo
      cdl%nvarlines = k
c data
      k = 0
      do n=1,cdl1%ndatalines
        k = k + 1
        cdl%datavalues(k) = cdl1%datavalues(n)
      enddo
      do n=1,cdl2%ndatalines
        k = k + 1
        cdl%datavalues(k) = cdl2%datavalues(n)
      enddo
      cdl%ndatalines = k
      return
      end subroutine merge_cdl

      subroutine print_cdl(cdl)
      type(cdl_type), intent(in) :: cdl
      integer :: k
      do k=1,cdl%nlines
        write(6,*) cdl%text(k)
      enddo
      end subroutine print_cdl

      subroutine defvar_cdl(grid,fid,cdl,varstr)
      use dd2d_utils, only : dist_grid
      use pario, only : defvar
      type(dist_grid) :: grid
      integer :: fid
      type(cdl_type) :: cdl
      character(len=*) :: varstr
      call assemble_cdl(cdl)
      call defvar(grid,fid,cdl%text(1:cdl%nlines),trim(varstr))
      return
      end subroutine defvar_cdl

      subroutine write_cdl(grid,fid,varstr,cdl)
      use dd2d_utils, only : dist_grid
      use pario, only : write_data
      type(dist_grid) :: grid
      integer :: fid
      type(cdl_type) :: cdl
      character(len=*) :: varstr
      call assemble_cdl(cdl)
      call write_data(grid,fid,trim(varstr),cdl%text(1:cdl%nlines))
      return
      end subroutine write_cdl

      end module cdl_mod
