subroutine cube2latlon(npx, npy, nlon, nlat, data_cs, data_ll)

 use ESMF_Mod
 use MAPL_ConstantsMod, only : pi=> MAPL_PI
 use fv_grid_utils_mod, only : gnomonic_grids
 use fv_grid_tools_mod, only : mirror_grid
 use CUB2LATLON_mod,    only : init_latlon_grid, &
                           read_c2l_weight,  write_c2l_weight,         &
                           get_c2l_weight,   do_c2l_interpolation
 use GHOST_CUBSPH_mod,   only: B_grid, A_grid, ghost_cubsph_update

 implicit none

 integer, intent(in) :: npx, npy, nlon, nlat
 real, dimension(npx , npy ), intent(in ) :: data_cs
 real, dimension(nlon, nlat), intent(out) :: data_ll

 integer :: ntiles=6
 integer :: npts
 integer :: ndims=2
 real(ESMF_KIND_R8), allocatable :: xs(:,:), ys(:,:)
 real(ESMF_KIND_R8), allocatable :: grid_global(:,:,:,:)
 real(ESMF_KIND_R8), allocatable :: sph_corner(:,:,:,:)

 real(ESMF_KIND_R8), allocatable :: xlon(:), ylat(:)

 real(ESMF_KIND_R8), save, dimension(:,:,:), allocatable :: c2l_weight                 
 integer           , save, dimension(:,:,:), allocatable :: c2l_index
 real(ESMF_KIND_R8), save, dimension(:,:,:,:), allocatable  :: elon_cubsph, elat_cubsph 
 real(ESMF_KIND_R8), save, dimension(:,:,:),   allocatable  :: elon_latlon, elat_latlon 
 logical, save :: do_init=.true.

 real(ESMF_KIND_R8) :: varmisval=1.e25
 real(ESMF_KIND_R8), allocatable :: var_cubsph(:,:,:)
 real(ESMF_KIND_R8), allocatable :: var_latlon(:,:)

 logical :: found
 integer :: grid_type = 0

 integer :: i,j,n,l,itile, j1,j2

  npts = npx+1

 if (do_init) then

  !--------------------------------------------------------------------!
  ! check for existing c2l_res files                                   !
  !--------------------------------------------------------------------!
  allocate(c2l_index(3,nlon,nlat),c2l_weight(4,nlon,nlat))
  allocate(elon_cubsph(3,0:npts,0:npts,ntiles), elon_latlon(3,nlon,nlat), &
           elat_cubsph(3,0:npts,0:npts,ntiles), elat_latlon(3,nlon,nlat))
  call read_c2l_weight(c2l_index, c2l_weight, nlon, nlat, npts, npts, ntiles, &
                       elon_cubsph, elat_cubsph, elon_latlon, elat_latlon, found)

  if (.not. found) then

  !--------------------------------------------------------------------!
  ! initialize cubed sphere grid                                       !
  !--------------------------------------------------------------------!
  allocate( xs(npts,npts) )
  allocate( ys(npts,npts) )
  allocate( grid_global(npts,npts,ndims,ntiles) )
  allocate( sph_corner(ndims,0:npts+1,0:npts+1,ntiles) )
  call gnomonic_grids(grid_type, npts-1, xs, ys)
  do j=1,npts
     do i=1,npts
        grid_global(i,j,1,1) = xs(i,j)
        grid_global(i,j,2,1) = ys(i,j)
     enddo
  enddo
! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi]
  call mirror_grid(grid_global, 0, npts, npts, 2, 6)

  do n=1,ntiles
     do j=1,npts
        do i=1,npts
!---------------------------------
! Shift the corner away from Japan
!---------------------------------
#ifdef SHIFT_WEST
! This will result in the corner close to east coast of China
           grid_global(i,j,1,n) = grid_global(i,j,1,n) - pi/18.
#endif
           if ( grid_global(i,j,1,n) < 0. )              &
                grid_global(i,j,1,n) = grid_global(i,j,1,n) + 2.*pi
           if (ABS(grid_global(i,j,1,1)) < 1.e-10) grid_global(i,j,1,1) = 0.0
           if (ABS(grid_global(i,j,2,1)) < 1.e-10) grid_global(i,j,2,1) = 0.0
        enddo
     enddo
  enddo
!---------------------------------
! Clean Up Corners
!---------------------------------
  grid_global(  1,1:npts,:,2)=grid_global(npts,1:npts,:,1)
  grid_global(  1,1:npts,:,3)=grid_global(npts:1:-1,npts,:,1)
  grid_global(1:npts,npts,:,5)=grid_global(1,npts:1:-1,:,1)
  grid_global(1:npts,npts,:,6)=grid_global(1:npts,1,:,1)
  grid_global(1:npts,  1,:,3)=grid_global(1:npts,npts,:,2)
  grid_global(1:npts,  1,:,4)=grid_global(npts,npts:1:-1,:,2)
  grid_global(npts,1:npts,:,6)=grid_global(npts:1:-1,1,:,2)
  grid_global(  1,1:npts,:,4)=grid_global(npts,1:npts,:,3)
  grid_global(  1,1:npts,:,5)=grid_global(npts:1:-1,npts,:,3)
  grid_global(npts,1:npts,:,3)=grid_global(1,1:npts,:,4)
  grid_global(1:npts,  1,:,5)=grid_global(1:npts,npts,:,4)
  grid_global(1:npts,  1,:,6)=grid_global(npts,npts:1:-1,:,4)
  grid_global(  1,1:npts,:,6)=grid_global(npts,1:npts,:,5)

  sph_corner(1,1:npts,1:npts,:) = grid_global(:,:,1,:)
  sph_corner(2,1:npts,1:npts,:) = grid_global(:,:,2,:)
  !------------------------------------------------------------------!
  ! do halo update                                                   !
  !------------------------------------------------------------------!
  do l=1,ntiles
     sph_corner(1:2,0     ,0     ,l)=0.
     sph_corner(1:2,npts+1,0     ,l)=0.
     sph_corner(1:2,0     ,npts+1,l)=0.
     sph_corner(1:2,npts+1,npts+1,l)=0.
     call ghost_cubsph_update(sph_corner(1,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1, 1, &
                              1, ntiles, 1, 1, l, B_grid)
     call ghost_cubsph_update(sph_corner(2,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1, 1, &
                              1, ntiles, 1, 1, l, B_grid)
  enddo
  deallocate ( xs )
  deallocate ( ys )
  deallocate ( grid_global )

  !--------------------------------------------------------------------!
  ! initialize latlon grid                                             !
  !--------------------------------------------------------------------!
  allocate(xlon(nlon), ylat(nlat))
  call init_latlon_grid(xlon, ylat, nlon, nlat)

  !--------------------------------------------------------------------!
  ! calculate weights for bilinear interpolation                       !
  ! from cubed sphere to latlon grid                                   !
  !--------------------------------------------------------------------! 
     call get_c2l_weight(sph_corner, npts, npts, ntiles, xlon, ylat, nlon, nlat,  &
                         c2l_index, c2l_weight, elon_cubsph, elat_cubsph,       &
                         elon_latlon, elat_latlon)
     call write_c2l_weight(c2l_index, c2l_weight, nlon, nlat, npts, npts, ntiles, &
                           elon_cubsph, elat_cubsph, elon_latlon, elat_latlon)

  deallocate ( xlon )
  deallocate ( ylat )
  deallocate ( sph_corner )

  endif ! found
  do_init = .false.
 endif ! do_init

  !--------------------------------------------------------------------!
  ! perform interpolation                                              !
  !--------------------------------------------------------------------!
  allocate ( var_cubsph(0:npts,0:npts,ntiles) )
  allocate ( var_latlon(nlon,nlat) )
  do itile=1,ntiles
     j1 = (npts-1)*(itile-1) + 1
     j2 = (npts-1)*(itile-1) + npts-1
     var_cubsph(1:npts-1,1:npts-1,itile)=data_cs(:,j1:j2)
  enddo
  do itile=1,ntiles
     call ghost_cubsph_update(var_cubsph, 0, npts, 0, npts, 1, 1, ntiles,  &
                              1, 1, itile, A_grid)
     call do_c2l_interpolation(var_cubsph(:,:,itile), 0, npts, 0, npts, 1, itile, 1, &
                               var_latlon, nlon, nlat, c2l_index, c2l_weight,                      &
                               .true., varmisval, .true.)
  enddo
  data_ll = var_latlon
  deallocate ( var_cubsph )
  deallocate ( var_latlon )

! NOTE: CAUTION: Code assumes interpolation is the same throughout execution !!!
! SAVE Variables to avoid costly Initialization of weights
!  deallocate( c2l_index )
!  deallocate( c2l_weight )
!  deallocate( elon_cubsph )
!  deallocate( elon_latlon )
!  deallocate( elat_cubsph )
!  deallocate( elat_latlon )
! SAVE Variables to avoid costly Initialization of weights

end subroutine cube2latlon

