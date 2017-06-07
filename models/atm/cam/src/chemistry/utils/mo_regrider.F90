      module mo_regrider
!---------------------------------------------------------------------
!	... General horizontal regriding
!---------------------------------------------------------------------

      use shr_kind_mod, only : r8 => shr_kind_r8
      use abortutils,   only : endrun
      use cam_logfile,  only : iulog

      implicit none

      integer, private, parameter  :: max_conv = 25
      integer, private             :: conv_cnt = 0
      real(r8), private, parameter :: max_diff = 1.e-4_r8

      save

      type grid_conv

      integer :: nfrom_lats, nto_lats, &
                 nfrom_lons, nto_lons, &
                 latl, latu, lati, &
                 lonl, lonu, loni, &
                 from_min_lat, from_max_lat, &                ! global indicies
                 wings                                        ! number of ghost wings
      integer, pointer, dimension(:) :: interp_lats, &
                                        interp_lons
      real(r8)    :: max_lat, min_lat
      real(r8)    :: max_lon, min_lon
      real(r8), pointer, dimension(:)    :: lat_del, lon_del
      real(r8), pointer, dimension(:)    :: from_lats, to_lats
      real(r8), pointer, dimension(:)    :: from_lons, to_lons

      logical :: do_lon, do_lat
      logical :: from_lats_mono_pos, from_lons_mono_pos
      logical :: active

      end type grid_conv

      type(grid_conv), private :: converter(max_conv)

      public  :: regrid_inti, regrid_3d, regrid_2d, regrid_1d, regrid_lat_limits, regrid_diagnostics
      private

      contains

      integer function regrid_inti( from_nlats, to_nlats, &
                                    from_nlons, to_nlons, &
                                    from_lons,  to_lons, &
                                    from_lats,  to_lats, &
                                    wing_cnt, &
                                    do_lons,    do_lats )
!---------------------------------------------------------------------
!	... Determine indicies and deltas for transform
!           Note : it is assumed that the latitude and longitude
!                  arrays are monotonic
!--------------------------------------------------------------------
      implicit none

!---------------------------------------------------------------------
!	... Dummy args
!---------------------------------------------------------------------
      integer, intent(in)  :: from_nlats, to_nlats, &
                              from_nlons, to_nlons, &
                              wing_cnt                                        ! number of wing terms
      real(r8), intent(in) :: from_lats(from_nlats), to_lats(to_nlats), &
                              from_lons(from_nlons), to_lons(to_nlons)
      logical, optional, intent(in) :: do_lons, do_lats

!---------------------------------------------------------------------
!	... Local variables
!---------------------------------------------------------------------
      integer :: from_lat, to_lat
      integer :: from_lon, to_lon
      integer :: astat
      integer :: i, j, jglb, m
      integer, dimension(1) :: max_ind, min_ind
      real(r8) :: target_lat, target_lon, pi, r2d
      logical  :: match, check_lats, check_lons, lat_xform, lon_xform

!---------------------------------------------------------------------
!	... Check if dimension transform is required or requested
!---------------------------------------------------------------------
      check_lats = .not. present( do_lats )
      if( .not. check_lats ) then
         check_lats = do_lats
      end if
      check_lons = .not. present( do_lons )
      if( .not. check_lons ) then
         check_lons = do_lons
      end if


!---------------------------------------------------------------------
!	... No transform requested; leave
!---------------------------------------------------------------------
      if( .not. check_lats .and. .not. check_lons ) then
	 regrid_inti = -1
	 return
      end if

!---------------------------------------------------------------------
!	... Check to see if from lat grid == to lat grid
!---------------------------------------------------------------------
      if( check_lats ) then
         lat_xform = from_nlats /= to_nlats
         if( .not. lat_xform ) then
            do j = 1,to_nlats
	       if( abs( from_lats(j) - to_lats(j) ) > max_diff*abs( from_lats(j) ) ) then
	          lat_xform = .true.
	          exit
	       end if
	    end do
         end if
      else
	 lat_xform = .false.
      end if
!---------------------------------------------------------------------
!	... Check to see if from lon grid == to lon grid
!---------------------------------------------------------------------
      if( check_lons ) then
         lon_xform = from_nlons /= to_nlons
         if( .not. lon_xform ) then
            do i = 1,to_nlons
	       if( abs( from_lons(i) - to_lons(i) ) > max_diff*abs( from_lons(i) ) ) then
	          lon_xform = .true.
	          exit
	       end if
	    end do
         end if
      else
	 lon_xform = .false.
      end if
!---------------------------------------------------------------------
!	... No transform necessary; leave
!---------------------------------------------------------------------
      if( .not. lat_xform .and. .not. lon_xform ) then
	 regrid_inti = 0
	 return
      end if

!---------------------------------------------------------------------
!	... Check for match with existing transform
!---------------------------------------------------------------------
      if( conv_cnt > 0 ) then
         do m = 1,conv_cnt
	    if( wing_cnt /= converter(m)%wings ) then
	       match = .false.
	       cycle
	    end if
            if( lat_xform .and. converter(m)%do_lat ) then
	       match = converter(m)%nfrom_lats == from_nlats
	       if( match ) then
	          do j = 1,from_nlats
	             if( ABS( converter(m)%from_lats(j) - from_lats(j) ) > max_diff*ABS( converter(m)%from_lats(j) ) ) then
		        match = .false.
		        exit
		     end if
	          end do
	       end if
	       if( .not. match ) then
	          cycle
	       end if
	       match = converter(m)%nto_lats == to_nlats
	       if( match ) then
	          do j = 1,to_nlats
	             if( abs( converter(m)%to_lats(j) - to_lats(j) ) > max_diff*abs( converter(m)%to_lats(j) ) ) then
		        match = .false.
		        exit
		     end if
	          end do
	       end if
            else if( lat_xform .eqv. converter(m)%do_lat ) then
	       match = .true.
	    else
	       match = .false.
	    end if
	    if( .not. match ) then
	       cycle
	    end if
            if( lon_xform .and. converter(m)%do_lon ) then
	       match = converter(m)%nfrom_lons == from_nlons
	       if( match ) then
	          do i = 1,from_nlons
	             if( abs( converter(m)%from_lons(i) - from_lons(i) ) > max_diff*abs( converter(m)%from_lons(i) ) )  then
		        match = .false.
		        exit
		     end if
	          end do
	       end if
	       if( .not. match ) then
	          cycle
	       end if
	       match = converter(m)%nto_lons == to_nlons
	       if( match ) then
	          do i = 1,to_nlons
	             if( ABS( converter(m)%to_lons(i) - to_lons(i) ) > max_diff*ABS( converter(m)%to_lons(i) ) ) then
		        match = .false.
		        exit
		     end if
	          end do
	       end if
            else if( lon_xform .eqv. converter(m)%do_lon ) then
	       match = .true.
	    else
	       match = .false.
	    end if
	    if( match ) then
	       exit
	    end if
	 end do
      else
         match = .false.
      end if
      if( match ) then
	 regrid_inti = m
	 return
      else
!---------------------------------------------------------------------
!	... Check for conversion count
!---------------------------------------------------------------------
         if( conv_cnt >= max_conv ) then
	    write(iulog,*) 'REGRID_INTI: Reached max conversion count of ',max_conv
	    regrid_inti = -2
	    return
         end if
         conv_cnt = conv_cnt + 1
      end if

      converter(conv_cnt)%do_lat = lat_xform
      converter(conv_cnt)%do_lon = lon_xform
      converter(conv_cnt)%wings  = wing_cnt

!---------------------------------------------------------------------
!	... New transform; store grids
!---------------------------------------------------------------------
      write(iulog,*) 'REGRID_INTI: Diagnostics for transform index = ',conv_cnt
      write(iulog,'(1x,''REGRID_INTI: from_nlats, to_nlats = '',2i6)') from_nlats, to_nlats
      write(iulog,'(1x,''REGRID_INTI: from_nlons, to_nlons = '',2i6)') from_nlons, to_nlons
      write(iulog,*) 'REGRID_INTI: lat_xform, lon_xform = ',lat_xform,lon_xform
      if( lat_xform ) then
         allocate( converter(conv_cnt)%from_lats(from_nlats),stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'REGRID_INTI: Failed to allocate from_lats array'
	    call endrun
         end if
         converter(conv_cnt)%from_lats(:) = from_lats(:)
         allocate( converter(conv_cnt)%to_lats(to_nlats),stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'REGRID_INTI: Failed to allocate to_lats array'
	    call endrun
         end if
         write(iulog,*) 'REGRID_INTI: size of to_lats = ',size( to_lats )
         write(iulog,*) 'REGRID_INTI: size of converter(conv_cnt)%to_lats = ',size( converter(conv_cnt)%to_lats )
         converter(conv_cnt)%to_lats(:) = to_lats(:)
#ifdef DEBUG
	 pi  = 4._r8 * atan( 1._r8 )
	 r2d = 180._r8/pi
         write(iulog,*) 'REGRID_INTI: to_lats (deg):'
         write(iulog,'(10F8.3)') to_lats(:to_nlats)*r2d
#endif
      end if
      if( lon_xform ) then
         allocate( converter(conv_cnt)%from_lons(from_nlons),stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'REGRID_INTI: Failed to allocate from_lons array'
	    call endrun
         end if
         converter(conv_cnt)%from_lons(:) = from_lons(:)
         allocate( converter(conv_cnt)%to_lons(to_nlons),stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'REGRID_INTI: Failed to allocate to_lons array'
	    call endrun
         end if
         write(iulog,*) 'REGRID_INTI: size of to_lons = ',size( to_lons )
         write(iulog,*) 'REGRID_INTI: size of converter(conv_cnt)%to_lons = ',size( converter(conv_cnt)%to_lons )
#ifdef DEBUG
         write(iulog,*) 'REGRID_INTI: to_lons (deg):'
         write(iulog,'(10F8.3)') to_lons(:to_nlons)*r2d
#endif
         converter(conv_cnt)%to_lons(:) = to_lons(:)
      end if

!---------------------------------------------------------------------
!	... Set "module" variables
!---------------------------------------------------------------------
      converter(conv_cnt)%nfrom_lats = from_nlats
      converter(conv_cnt)%nto_lats = to_nlats
      converter(conv_cnt)%nfrom_lons = from_nlons
      converter(conv_cnt)%nto_lons = to_nlons
      write(iulog,*) 'REGRID_INTI: size of to_lons = ',size( to_lons )

      if( converter(conv_cnt)%do_lat ) then
         max_ind(:) = maxloc( from_lats(:) )
         min_ind(:) = minloc( from_lats(:) )
         converter(conv_cnt)%max_lat = from_lats(max_ind(1))
         converter(conv_cnt)%min_lat = from_lats(min_ind(1))
         if( max_ind(1) >= min_ind(1) ) then
	    converter(conv_cnt)%latl = 1
	    converter(conv_cnt)%latu = from_nlats
	    converter(conv_cnt)%lati = 1
	    converter(conv_cnt)%from_lats_mono_pos = .true.
         else
	    converter(conv_cnt)%latl = from_nlats
	    converter(conv_cnt)%latu = 1
	    converter(conv_cnt)%lati = -1
	    converter(conv_cnt)%from_lats_mono_pos = .false.
         end if
      end if

      if( converter(conv_cnt)%do_lon ) then
         max_ind(:) = maxloc( from_lons(:) )
         min_ind(:) = minloc( from_lons(:) )
         converter(conv_cnt)%max_lon = from_lons(max_ind(1))
         converter(conv_cnt)%min_lon = from_lons(min_ind(1))
         if( max_ind(1) >= min_ind(1) ) then
	    converter(conv_cnt)%lonl = 1
	    converter(conv_cnt)%lonu = from_nlons
	    converter(conv_cnt)%loni = 1
	    converter(conv_cnt)%from_lons_mono_pos = .true.
         else
	    converter(conv_cnt)%lonl = from_nlons
	    converter(conv_cnt)%lonu = 1
	    converter(conv_cnt)%loni = -1
	    converter(conv_cnt)%from_lons_mono_pos = .false.
         end if
      end if

      if( converter(conv_cnt)%do_lat ) then
!---------------------------------------------------------------------
!	... Allocate interpolation latitude indicies
!---------------------------------------------------------------------
         allocate( converter(conv_cnt)%interp_lats(to_nlats),stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'REGRID_INTI: Failed to allocate interp lats array'
	    call endrun
         end if
!---------------------------------------------------------------------
!	... Allocate interpolation latitude deltas
!---------------------------------------------------------------------
         allocate( converter(conv_cnt)%lat_del(to_nlats),stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'REGRID_INTI: Failed to allocate lat del array'
	    call endrun
         end if
!---------------------------------------------------------------------
!	... Set interpolation latitude indicies and deltas
!---------------------------------------------------------------------
         do to_lat = 1,converter(conv_cnt)%nto_lats
            target_lat = to_lats(to_lat)
	    if( target_lat <= converter(conv_cnt)%min_lat ) then
	       converter(conv_cnt)%lat_del(to_lat) = 0._r8
	       converter(conv_cnt)%interp_lats(to_lat) = converter(conv_cnt)%latl
	    else if( target_lat >= converter(conv_cnt)%max_lat ) then
	       converter(conv_cnt)%lat_del(to_lat) = 1._r8
	       if( converter(conv_cnt)%from_lats_mono_pos ) then
	          converter(conv_cnt)%interp_lats(to_lat) = converter(conv_cnt)%latu - 1
	       else
	          converter(conv_cnt)%interp_lats(to_lat) = converter(conv_cnt)%latu + 1
	       end if
	    else
	       do from_lat = converter(conv_cnt)%latl,converter(conv_cnt)%latu,converter(conv_cnt)%lati
	          if( target_lat < from_lats(from_lat) ) then
		     j = from_lat - converter(conv_cnt)%lati
	             converter(conv_cnt)%interp_lats(to_lat) = min( converter(conv_cnt)%nfrom_lats,max( 1,j ) )
	             converter(conv_cnt)%lat_del(to_lat) = &
                                        (target_lat - from_lats(j))/(from_lats(from_lat) - from_lats(j))
	             exit
	          end if
	       end do
            end if
         end do
!        jglb = 1
!        converter(conv_cnt)%from_max_lat = min( converter(conv_cnt)%nfrom_lats, &
!                                                converter(conv_cnt)%interp_lats(jglb) + converter(conv_cnt)%lati )
         converter(conv_cnt)%from_max_lat = converter(conv_cnt)%nfrom_lats
! converter(conv_cnt)%from_min_lat = max( 1,converter(conv_cnt)%interp_lats(jglb) )
	 converter(conv_cnt)%from_min_lat = 1
      end if

      if( converter(conv_cnt)%do_lon ) then
!---------------------------------------------------------------------
!	... Allocate interpolation longitude indicies
!---------------------------------------------------------------------
         allocate( converter(conv_cnt)%interp_lons(to_nlons),stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'REGRID_INTI: Failed to allocate interp lon array'
	    call endrun
         end if
!---------------------------------------------------------------------
!	... Allocate interpolation longitude deltas
!---------------------------------------------------------------------
         allocate( converter(conv_cnt)%lon_del(to_nlons),stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'REGRID_INTI: Failed to allocate lon del array'
	    call endrun
         end if
!---------------------------------------------------------------------
!	... Set interpolation longitude indicies and deltas
!---------------------------------------------------------------------
         do to_lon = 1,converter(conv_cnt)%nto_lons
            target_lon = to_lons(to_lon)
	    if( target_lon <= converter(conv_cnt)%min_lon ) then
	       converter(conv_cnt)%lon_del(to_lon) = 0._r8
	       converter(conv_cnt)%interp_lons(to_lon) = converter(conv_cnt)%lonl
	    else if( target_lon >= converter(conv_cnt)%max_lon ) then
	       converter(conv_cnt)%lon_del(to_lon) = 1._r8
	       if( converter(conv_cnt)%from_lons_mono_pos ) then
	          converter(conv_cnt)%interp_lons(to_lon) = converter(conv_cnt)%lonu - 1
	       else
	          converter(conv_cnt)%interp_lons(to_lon) = converter(conv_cnt)%lonu + 1
	       end if
	    else
	       do from_lon = converter(conv_cnt)%lonl,converter(conv_cnt)%lonu,converter(conv_cnt)%loni
	          if( target_lon < from_lons(from_lon) ) then
		     i = from_lon - converter(conv_cnt)%loni
	             converter(conv_cnt)%interp_lons(to_lon) = min( converter(conv_cnt)%nfrom_lons,max( 1,i ) )
	             converter(conv_cnt)%lon_del(to_lon) = &
                                        (target_lon - from_lons(i))/(from_lons(from_lon) - from_lons(i))
	             exit
	          end if
	       end do
            end if
         end do
      end if

      regrid_inti = conv_cnt

      end function regrid_inti

      subroutine regrid_3d( from_field, to_field, index, to_lat_min, to_lat_max, do_poles, scaling )
!--------------------------------------------------------------------
!	... Regrid horizontal data
!           Note: this subroutine works on latitude "slices"
!--------------------------------------------------------------------

      use ppgrid, only : pver
      use scamMod, only: single_column
      use dycore, only : dycore_is
      implicit none

!--------------------------------------------------------------------
!	... Dummy args
!--------------------------------------------------------------------
      integer, intent(in) :: index
      integer, intent(in) :: to_lat_min, to_lat_max                     ! globals
      real(r8), optional, intent(in)    :: scaling
      real(r8), intent(in)  :: from_field(:,:,:)
      real(r8), intent(out) :: to_field(:,:,:)
      logical, optional, intent(in) :: do_poles

!--------------------------------------------------------------------
!	... Local variables
!--------------------------------------------------------------------
      integer   :: i, j, k, astat
      integer   :: nlons, nlats, offset
      real(r8)  :: temp, rnlons
      real(r8), allocatable :: wrk(:,:)

      if(dycore_is('UNSTRUCTURED')) then
         call endrun('unstructured grid regriding not supported in mo_regrider')
      end if


      nlats = to_lat_max - to_lat_min + 1
      if( index /= 0 ) then
!--------------------------------------------------------------------
!	... Check index for validity
!--------------------------------------------------------------------
         if( index < 1 .or. index > conv_cnt ) then
	    write(iulog,'(''REGRID_3D: '',3x,'' is out of range'')') index
	    call endrun
         end if
!--------------------------------------------------------------------
!	... Allocate work array
!--------------------------------------------------------------------
	 nlons = converter(index)%nto_lons
         allocate( wrk(nlons,nlats),stat=astat )
         if( astat /= 0 ) then
	    write(iulog,*) 'REGRID_3D: Failed to allocate work array'
	    call endrun
         end if
!--------------------------------------------------------------------
!	... Latitude interp
!--------------------------------------------------------------------
         do k = 1,pver
            call regrid_horiz( from_field(:,:,k), wrk, to_lat_min, to_lat_max, index )
	    do j = 1,nlats
	       do i = 1,nlons
	          to_field(i,k,j) = wrk(i,j)
	       end do
	    end do
         end do
         deallocate( wrk )
      else
!--------------------------------------------------------------------
!	... Transparent transform
!--------------------------------------------------------------------
	 nlons  = size(from_field,dim=1)
	 offset = (size(from_field,dim=2) - nlats)/2
         do k = 1,pver
	    do j = 1,nlats
	       do i = 1,nlons
	          to_field(i,k,j) = from_field(i,j+offset,k)
	       end do
	    end do
	 end do
      end if

      if( present(scaling) ) then
         if( scaling /= 1._r8 ) then
	    do j = 1,nlats
               do k = 1,pver
	          do i = 1,nlons
		     to_field(i,k,j) = scaling * to_field(i,k,j)
	          end do
	       end do
	    end do
         end if
      end if

      if( present(do_poles).and..not.single_column ) then
         if( do_poles ) then
	    rnlons = 1._r8/nlons
            do k = 1,pver
	       temp = sum( to_field(:,k,2) )*rnlons
	       to_field(:,k,1) = temp
	    end do
            do k = 1,pver
	       temp = sum( to_field(:,k,nlats-1) )*rnlons
	       to_field(:,k,nlats) = temp
	    end do
         end if
      end if

      end subroutine regrid_3d

      subroutine regrid_2d( from_field, to_field, index, to_lat_min, to_lat_max, do_poles, scaling )
!--------------------------------------------------------------------
!	... Regrid horizontal data
!           Note: this subroutine works on horizontal "slices"
!--------------------------------------------------------------------

      use scamMod, only: single_column
      implicit none

!--------------------------------------------------------------------
!	... Dummy args
!--------------------------------------------------------------------
      integer, intent(in) :: index
      integer, intent(in) :: to_lat_min, to_lat_max                    ! globals
      real(r8), optional, intent(in)    :: scaling
      real(r8), intent(in)    :: from_field(:,:)
      real(r8), intent(inout) :: to_field(:,:)
      logical, optional, intent(in) :: do_poles

!--------------------------------------------------------------------
!	... Local variables
!--------------------------------------------------------------------
      integer :: i, j
      integer :: nlons, nlats, offset
      real(r8)    :: temp, rnlons

      nlats = to_lat_max - to_lat_min + 1
      if( index /= 0 ) then
!--------------------------------------------------------------------
!	... Check index for validity
!--------------------------------------------------------------------
         if( index < 1 .or. index > conv_cnt ) then
	    write(iulog,'(''REGIRD_2D: '',3x,'' is out of range'')') index
	    call endrun
         end if
         call regrid_horiz( from_field, to_field, to_lat_min, to_lat_max, index )
	 nlons  = converter(index)%nto_lons
      else
!--------------------------------------------------------------------
!	... Transparent transform
!--------------------------------------------------------------------
	 nlons  = size(from_field,dim=1)
	 offset = (size(from_field,dim=2) - nlats)/2
         do j = 1,nlats
            do i = 1,nlons
	       to_field(i,j) = from_field(i,j+offset)
	    end do
	 end do
      end if

      if( present(scaling) ) then
         if( scaling /= 1._r8 ) then
	    do j = 1,nlats
	       do i = 1,nlons
	          to_field(i,j) = scaling * to_field(i,j)
	       end do
	    end do
         end if
      end if

      if( present(do_poles).and..not.single_column ) then
         if( do_poles ) then
	    rnlons = 1._r8/nlons
	    temp = sum( to_field(:,2) )*rnlons
	    to_field(:,1) = temp
	    temp = sum( to_field(:,nlats-1) )*rnlons
	    to_field(:,nlats) = temp
         end if
      end if

      end subroutine regrid_2d

      subroutine regrid_horiz( from_field, to_field, latl, latu, index )
!--------------------------------------------------------------------
!	... Regrid horizontal data
!           Note: this subroutine works on horizontal "slices"
!--------------------------------------------------------------------

      implicit none

!--------------------------------------------------------------------
!	... Dummy args
!--------------------------------------------------------------------
      integer, intent(in)   :: latl, latu                                  ! globals
      integer, intent(in)   :: index
      real(r8), intent(in)  :: from_field(:,:)
      real(r8), intent(out) :: to_field(:,:)

!--------------------------------------------------------------------
!	... Local variables
!--------------------------------------------------------------------
      integer :: j, ji, ji1, jmin, jmax, jloc, jl, ju
      integer :: i, ii, ii1
      real(r8)    :: wrk(converter(index)%nto_lons,converter(index)%from_min_lat:converter(index)%from_max_lat)

      jmax = size( from_field,dim=2 )
      jl   = converter(index)%from_min_lat
      ju   = min( jl + jmax - 1,converter(index)%from_max_lat )
!--------------------------------------------------------------------
!	... First longitude interp
!--------------------------------------------------------------------
      if( converter(index)%do_lon ) then
         do i = 1,converter(index)%nto_lons
            ii = converter(index)%interp_lons(i)
	    ii1 = ii + converter(index)%loni
	    wrk(i,jl:ju) = from_field(ii,1:jmax) &
                           + converter(index)%lon_del(i) * (from_field(ii1,1:jmax) - from_field(ii,1:jmax))
         end do
      else
	 wrk(:,jl:ju) = from_field(:,1:jmax)
      end if
!--------------------------------------------------------------------
!	... Then latitude interp
!--------------------------------------------------------------------
      if( converter(index)%do_lat ) then
	 jmin = max( 1,latl )
	 jmax = min( converter(index)%nto_lats,latu )
         do j = jmin,jmax
	    ji = converter(index)%interp_lats(j)
	    ji1 = ji + converter(index)%lati
	    jloc = j - jmin + 1
	    to_field(:,jloc) = wrk(:,ji) &
                             + converter(index)%lat_del(j) * (wrk(:,ji1) - wrk(:,ji))
         end do
      else
	 to_field(:,:) = wrk(:,:)
      end if

      end subroutine regrid_horiz

      subroutine regrid_1d( from_field, to_field, index, scaling, do_lat, to_lat_min, to_lat_max, do_lon, to_lon_min, to_lon_max )
!--------------------------------------------------------------------
!	... Regrid horizontal data
!           Note: this subroutine works on a horizontal "line"
!--------------------------------------------------------------------

      implicit none

!--------------------------------------------------------------------
!	... Dummy args
!--------------------------------------------------------------------
      integer, intent(in) :: index
      integer, optional, intent(in) :: to_lat_min, to_lat_max
      integer, optional, intent(in) :: to_lon_min, to_lon_max
      real(r8), intent(in)    :: from_field(:)
      real(r8), intent(out)   :: to_field(:)
      real(r8), optional, intent(in) :: scaling
      logical, optional, intent(in) :: do_lat, do_lon

!--------------------------------------------------------------------
!	... Local variables
!--------------------------------------------------------------------
      integer :: j, ji, ji1, offset
      integer :: jmin, jmax, jloc
      integer :: i, ii, ii1
      integer :: nlons, nlats
      integer :: astat
      real(r8), allocatable :: wrk(:)

!--------------------------------------------------------------------
!	... Check index for validity
!--------------------------------------------------------------------
      if( index < 0 .or. index > conv_cnt ) then
	 write(iulog,'(''REGIRD_1D: '',3x,'' is out of range'')') index
	 call endrun
      end if
      if( present(do_lat) ) then
         if( index /= 0 ) then
!--------------------------------------------------------------------
!	... Latitude interp
!--------------------------------------------------------------------
	    nlats = to_lat_max - to_lat_min + 1
	    allocate( wrk(nlats),stat=astat )
	    if( astat /= 0 ) then
	       write(iulog,*) 'REGRID_1D: Failed to allocate wrk space'
	       call endrun
	    end if
	    jmin = max( 1,to_lat_min )
	    jmax = min( converter(index)%nto_lats,to_lat_max )
            offset = 1 - converter(index)%from_min_lat
            do j = jmin,jmax
               ji = converter(index)%interp_lats(j) + offset
	       ji1 = ji + converter(index)%lati
	       jloc = j - jmin + 1
	       wrk(jloc) = from_field(ji) &
                           + converter(index)%lat_del(j) * (from_field(ji1) - from_field(ji))
            end do
            to_field(1:nlats) = wrk(1:nlats)
	 else
	    nlats = size( from_field )
	    to_field(1:nlats) = from_field(1:nlats)
	 end if
	 if( present(scaling) ) then
	    if( scaling /= 1._r8 ) then
               to_field(1:nlats) = scaling*to_field(1:nlats)
	    end if
	 end if
      else if( present(do_lon ) ) then
         if( index /= 0 ) then
!--------------------------------------------------------------------
!	... Check dimensions
!--------------------------------------------------------------------
	    if( .not. converter(index)%do_lon ) then
	       write(iulog,*) 'REGRID_1D: Requesting lon interp; not set in intialization'
	       call endrun
	    end if
	    if( size( from_field ) /= converter(index)%nfrom_lons ) then
	       write(iulog,*) 'REGRID_1D: Input field does not match module dimension'
	       call endrun
	    end if
	    if( size( to_field ) /= converter(index)%nto_lons ) then
	       write(iulog,*) 'REGRID_1D: Output field does not match module dimension'
	       call endrun
	    end if
!--------------------------------------------------------------------
!	... Lontitude interp
!--------------------------------------------------------------------
	    nlons = converter(index)%nto_lons
	    allocate( wrk(nlons),stat=astat )
	    if( astat /= 0 ) then
	       write(iulog,*) 'REGRID_1D: Failed to allocate wrk space'
	       call endrun
	    end if
            do i = 1,nlons
               ii = converter(index)%interp_lons(i)
	       ii1 = ii + converter(index)%loni
	       wrk(i) = from_field(ii) &
                        + converter(index)%lon_del(i) * (from_field(ii1) - from_field(ii))
            end do
            to_field(1:nlons) = wrk(1:nlons)
	 else
	    nlons = size( from_field )
	    to_field(1:nlons) = from_field(1:nlons)
	 end if
	 if( present(scaling) ) then
	    if( scaling /= 1._r8 ) then
               to_field(1:nlons) = scaling*to_field(1:nlons)
	    end if
	 end if
      end if
      if( allocated( wrk ) ) then
         deallocate( wrk )
      end if

      end subroutine regrid_1d

      function regrid_lat_limits( index )
!--------------------------------------------------------------------
!	... Return the from latitude limits
!--------------------------------------------------------------------

      implicit none

!--------------------------------------------------------------------
!	... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: index

!--------------------------------------------------------------------
!	... Function declaration
!--------------------------------------------------------------------
      integer :: regrid_lat_limits(2)

      regrid_lat_limits(:2) = (/ converter(index)%from_min_lat, converter(index)%from_max_lat /)

      end function regrid_lat_limits

      subroutine regrid_diagnostics( index )

      implicit none

      integer, intent(in) :: index

      if( index /= 0 ) then
!--------------------------------------------------------------------
!	... Check index for validity
!--------------------------------------------------------------------
         if( index < 1 .or. index > conv_cnt ) then
	    write(iulog,'(''REGRID_DIAGNOSTICS: '',3x,'' is out of range'')') index
	    call endrun
         end if
	 write(iulog,*) ' '
	 write(iulog,*) '-------------------------------------------------------------------------------------------'
	 write(iulog,*) 'Regrid diagnostics for index ',index
	 write(iulog,*) 'Lon, Lat xform = ',converter(index)%do_lon,converter(index)%do_lat
	 if( converter(index)%do_lat ) then
	    write(iulog,*) 'Number from lats, to lats = ',converter(index)%nfrom_lats,converter(index)%nto_lats
	    write(iulog,*) 'latl, latu, lati = ',converter(index)%latl,converter(index)%latu,converter(index)%lati
	    write(iulog,*) 'From min,max lat = ',converter(index)%from_min_lat,converter(index)%from_max_lat
	    write(iulog,*) 'Wing count = ',converter(index)%wings
	    write(iulog,*) 'From lats monotically increasing = ',converter(index)%from_lats_mono_pos
	    write(iulog,*) 'Lat interp indicies'
	    write(iulog,'(10i5)') converter(index)%interp_lats(:)
	    write(iulog,*) 'Lat interp delta'
	    write(iulog,'(1p,5e22.15)') converter(index)%lat_del(:)
	 end if
	 if( converter(index)%do_lon ) then
	    write(iulog,*) ' '
	    write(iulog,*) 'Number from lons, to lons = ',converter(index)%nfrom_lons,converter(index)%nto_lons
	    write(iulog,*) 'lonl, lonu, loni = ',converter(index)%lonl,converter(index)%lonu,converter(index)%loni
	    write(iulog,*) 'From lons monotically increasing = ',converter(index)%from_lons_mono_pos
	    write(iulog,*) 'Lon interp indicies'
	    write(iulog,'(10i5)') converter(index)%interp_lons(:)
	    write(iulog,*) 'Lon interp delta size = ',size( converter(index)%lon_del )
	    write(iulog,'(1p,5e22.15)') converter(index)%lon_del(:)
	 end if
	 write(iulog,*) '-------------------------------------------------------------------------------------------'
	 write(iulog,*) ' '
      end if

      end subroutine regrid_diagnostics

      end module mo_regrider
