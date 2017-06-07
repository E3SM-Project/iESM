

      module mo_sulf
!---------------------------------------------------------------
!	... Annual cycle for sulfur
!---------------------------------------------------------------

      use shr_kind_mod, only : r8 => shr_kind_r8

      use abortutils,   only : endrun
      use wrap_nf
      use cam_logfile,  only : iulog

      implicit none

      private
      public  :: sulf_inti, set_sulf_time, sulf_interp

      save

      integer :: time_cnt = 0, &
                 last, &
                 next, &
                 sulf_nlevs
      real(r8) :: dels
      real(r8) :: sulf_P0, sp0tmp(1)
      real(r8), allocatable :: times(:)
      real(r8), allocatable :: psi(:,:,:), &
                               sulfatei(:,:,:,:), &
                               sulf_hyam(:), &
                               sulf_hybm(:)
      logical :: cyclical = .true.
      logical :: read_sulf = .true.

      contains 

      subroutine sulf_inti( sulf_file )
!-----------------------------------------------------------------------
! 	... Open netCDF file containing annual sulfur data.  Initialize
!           arrays with the data to be interpolated to the current time.
!
!           It is assumed that the time coordinate is increasing
!           and represents calendar days; range = [1.,366.).
!-----------------------------------------------------------------------

      use mo_regrider,   only : regrid_inti, regrid_lat_limits, regrid_2d
      use phys_grid,     only : clat_p, clon_p
      use ioFileMod,     only : getfil
      use spmd_utils,    only : masterproc
      use physconst,     only : pi
#ifdef SPMD
      use mpishorthand,  only : mpicom, mpir8, mpiint
#endif
      use mo_chem_utls,  only : get_spc_ndx, get_rxt_ndx
      use dyn_grid,      only : get_dyn_grid_parm
      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: sulf_file

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: gndx = 0, &
                 ncid = 0, &
                 nlon, &
                 nlat
      integer  :: dimid_lat, dimid_lon, dimid_lev
      integer  :: &
        ierr, &              ! return code
        times_id, &
        vid, k, month, &
	jl, ju
      integer  ::  jlim_in(2)
      integer  ::  start(4), count(4)
      integer  ::  plon, plat
      real(r8) ::  d2r
      real(r8), allocatable :: lat(:), lon(:), &
           psi_in(:,:,:), &
           sulfatei_in(:,:,:,:)
      character(len=256) :: locfn

      integer :: ndxs(5), so4_ndx

      ndxs(1) = get_rxt_ndx( 'usr_N2O5_aer' )
      ndxs(2) = get_rxt_ndx( 'usr_NO3_aer' )
      ndxs(3) = get_rxt_ndx( 'usr_NO2_aer' )
      ndxs(4) = get_rxt_ndx( 'usr_HO2_aer' )
      ndxs(5) = get_rxt_ndx( 'het1' )
      so4_ndx = get_spc_ndx('SO4')

      read_sulf = any( ndxs > 0) .and. (so4_ndx < 0)

      if ( .not. read_sulf ) return
      plon = get_dyn_grid_parm('plon')
      plat = get_dyn_grid_parm('plat')
Masterproc_only : &
      if( masterproc ) then
!-----------------------------------------------------------------------
!     	... Open netcdf file
!-----------------------------------------------------------------------
         call getfil (sulf_file, locfn, 0)
         call wrap_open (trim(locfn), NF_NOWRITE, ncid)

!-----------------------------------------------------------------------
!     	... Inquire about file
!-----------------------------------------------------------------------
      ierr = nf_inq_unlimdim( ncid, times_id )
      if( ierr /= NF_NOERR ) then
         write(iulog,*) 'sulf_inti : Failed to get unlimited dimension id'
	 call endrun
      end if
      call wrap_inq_dimlen( ncid, times_id, time_cnt )
      if( time_cnt /= 12 ) then
         write(iulog,*) 'sulf_inti : There are ',time_cnt,' times; expecting 12'
	 call endrun
      end if

!-----------------------------------------------------------------------
!     	... Allocate space for time coordinate data
!-----------------------------------------------------------------------
      allocate( times(time_cnt), stat = ierr )
      if( ierr /= 0 ) then
	 write(iulog,*) 'sulf_inti : Failed to allocate times array; error = ',ierr
	 call endrun
      end if
!-----------------------------------------------------------------------
!     	... Get time coordinate
!-----------------------------------------------------------------------
      call wrap_inq_varid( ncid, 'time', vid )
      call wrap_get_var_realx( ncid, vid, times )

      write(iulog,*) 'sulf_inti: times' 
      write(iulog,'(1p,5e21.13)') times

!-----------------------------------------------------------------------
!     	... Get vertical coordinate
!-----------------------------------------------------------------------
      call wrap_inq_dimid( ncid, 'lev', dimid_lev )
      call wrap_inq_dimlen( ncid, dimid_lev, sulf_nlevs )
      allocate( sulf_hyam(sulf_nlevs), stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'sulf_inti: hyam allocation error = ',ierr
         call endrun
      end if
      allocate( sulf_hybm(sulf_nlevs), stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'sulf_inti: hybm allocation error = ',ierr
         call endrun
      end if
      call wrap_inq_varid( ncid, 'hyam', vid )
      call wrap_get_var_realx( ncid, vid, sulf_hyam )
      call wrap_inq_varid( ncid, 'hybm', vid )
      call wrap_get_var_realx( ncid, vid, sulf_hybm )
      call wrap_inq_varid( ncid, 'P0', vid )
      call wrap_get_var_realx( ncid, vid, sP0tmp )
      sulf_P0=sP0tmp(1)
      write(iulog,*) 'sulf_inti: P0 = ',sulf_P0

!-----------------------------------------------------------------------
!     	... Get latitude and longitude
!-----------------------------------------------------------------------
      call wrap_inq_dimid( ncid, 'lat', dimid_lat )
      call wrap_inq_dimlen( ncid, dimid_lat, nlat )
      allocate( lat(nlat), stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'sulf_inti: lat allocation error = ',ierr
         call endrun
      end if
      call wrap_inq_varid( ncid, 'lat', vid )
      call wrap_get_var_realx( ncid, vid, lat )
      d2r = pi/180._r8
      lat(:nlat) = lat(:nlat) * d2r
 
      call wrap_inq_dimid( ncid, 'lon', dimid_lon )
      call wrap_inq_dimlen( ncid, dimid_lon, nlon )
      allocate( lon(nlon), stat=ierr )
      if( ierr /= 0 ) then
         write(iulog,*) 'sulf_inti: lon allocation error = ',ierr
         call endrun
      end if
      call wrap_inq_varid( ncid, 'lon', vid )
      call wrap_get_var_realx( ncid, vid, lon )
      lon(:nlon) = lon(:nlon) * d2r

!-----------------------------------------------------------------------
!     	... Get grid interp limits
!-----------------------------------------------------------------------
      gndx = regrid_inti( nlat, plat, &
                          nlon, plon, &
                          lon,  clon_p, &
                          lat,  clat_p, &
                          0, &
                          do_lons=.true.,do_lats=.true. )
      if( ierr /= 0 ) then
         write(iulog,*) 'sulf_inti: Failed to deallocate lat,lon; ierr = ',ierr
         call endrun
      end if
      if( gndx /= 0 )then
         jlim_in = regrid_lat_limits( gndx )
      else
         jlim_in = (/ 1, nlat /)
      end if
      jl = 1
      ju = plat

      write(iulog,'(''sulf_inti: gndx='',i2,'', grid limits = '',2i4,'', jl,ju='',2i4)') &
         gndx, jlim_in, jl, ju

      deallocate( lat, lon, stat=ierr )

      allocate( psi_in(nlon,jlim_in(1):jlim_in(2),time_cnt), stat=ierr)
      if( ierr /= 0 ) then
         write(iulog,*) 'sulf_inti: psi_in allocation error = ',ierr
         call endrun
      end if
      allocate( sulfatei_in(nlon,sulf_nlevs,jlim_in(1):jlim_in(2),time_cnt), stat=ierr)
      if( ierr /= 0 ) then
         write(iulog,*) 'sulf_inti: sulfatei_in allocation error = ',ierr
         call endrun
      end if
!-----------------------------------------------------------------------
!	... Read the sulfate inputs
!-----------------------------------------------------------------------
      call wrap_inq_varid( ncid, 'PS', vid )
      start(:3) = (/ 1, jlim_in(1), 1 /)
      count(:3) = (/ nlon, jlim_in(2) - jlim_in(1) + 1, time_cnt /)
      call wrap_get_vara_realx( ncid, vid, start(:3), count(:3), psi_in )

      call wrap_inq_varid( ncid, 'SULFATE', vid )
      start(:) = (/ 1, 1, jlim_in(1), 1 /)
      count(:) = (/ nlon, sulf_nlevs, jlim_in(2) - jlim_in(1) + 1, time_cnt /)
      call wrap_get_vara_realx( ncid, vid, start, count, sulfatei_in )
!-----------------------------------------------------------------------
!	... Regrid sulfate inputs
!-----------------------------------------------------------------------
      allocate( psi(plon,plat,time_cnt ), stat=ierr)
      if( ierr /= 0 ) then
         write(iulog,*) 'sulf_inti: psi allocation error = ',ierr
         call endrun
      end if
      allocate( sulfatei(plon,sulf_nlevs,plat,time_cnt), stat=ierr)
      if( ierr /= 0 ) then
         write(iulog,*) 'sulf_inti: sulfatei allocation error = ',ierr
         call endrun
      end if
      do month = 1,time_cnt
         call regrid_2d( psi_in(:,jlim_in(1):jlim_in(2),month), psi(:,:,month), gndx, &
	                 jl, ju, do_poles=.true. )
         do k = 1,sulf_nlevs
            call regrid_2d( sulfatei_in(:,k,jlim_in(1):jlim_in(2),month), sulfatei(:,k,:,month), gndx, &
	                    jl, ju, do_poles=.true. )
	 end do
      end do

      deallocate( psi_in, sulfatei_in )

      call wrap_close( ncid )

      end if Masterproc_only

#ifdef SPMD
      call mpibarrier( mpicom )
      call mpibcast( time_cnt, 1, mpiint, 0, mpicom )
      call mpibcast( sulf_nlevs, 1, mpiint, 0, mpicom )
      if( .not. masterproc ) then
         allocate( sulf_hyam(sulf_nlevs), stat=ierr )
         if( ierr /= 0 ) then
            write(iulog,*) 'sulf_inti: hyam allocation error = ',ierr
            call endrun
         end if
         allocate( sulf_hybm(sulf_nlevs), stat=ierr )
         if( ierr /= 0 ) then
            write(iulog,*) 'sulf_inti: hybm allocation error = ',ierr
            call endrun
         end if
         allocate( times(time_cnt), stat = ierr )
         if( ierr /= 0 ) then
	    write(iulog,*) 'sulf_inti : Failed to allocate times array; error = ',ierr
	    call endrun
         end if
         allocate( psi(plon,plat,time_cnt ), stat=ierr)
         if( ierr /= 0 ) then
            write(iulog,*) 'sulf_inti: psi allocation error = ',ierr
            call endrun
         end if
         allocate( sulfatei(plon,sulf_nlevs,plat,time_cnt ), stat=ierr)
         if( ierr /= 0 ) then
            write(iulog,*) 'sulf_inti: sulfatei allocation error = ',ierr
            call endrun
         end if
      end if
      call mpibcast( sulf_P0, 1, mpir8, 0, mpicom )
      call mpibcast( times, time_cnt, mpir8, 0, mpicom )
      call mpibcast( sulf_hyam, sulf_nlevs, mpir8, 0, mpicom )
      call mpibcast( sulf_hybm, sulf_nlevs, mpir8, 0, mpicom )
      call mpibcast( psi, plon*plat*time_cnt, mpir8, 0, mpicom )
      call mpibcast( sulfatei, plon*sulf_nlevs*plat*time_cnt, mpir8, 0, mpicom )
#endif

      end subroutine sulf_inti

      subroutine set_sulf_time( calday )
!--------------------------------------------------------------------
!	... Check and set time interpolation indicies
!--------------------------------------------------------------------

      implicit none

!--------------------------------------------------------------------
!	... Dummy args
!--------------------------------------------------------------------
      real(r8), intent(in)    :: calday

!--------------------------------------------------------------------
!	... Local variables
!--------------------------------------------------------------------
      real(r8), parameter :: dayspy = 365._r8
      integer  ::  m, upper
      real(r8) ::  numer, denom

      if ( .not. read_sulf ) return

!--------------------------------------------------------------------
!	... Setup the time interpolation
!--------------------------------------------------------------------
      if( calday < times(1) ) then
	 next = 1
	 last = 12
      else
	 if( times(12) < dayspy ) then
	    upper = 12
	 else
	    upper = 11
	 end if
         do m = upper,1,-1
	    if( calday >= times(m) ) then
	       exit
	    end if
         end do
	 last = m
	 next = mod( m,12 ) + 1
      end if
      numer = calday - times(last)
      denom = times(next) - times(last)
      if( numer < 0._r8 ) then
	    numer = dayspy + numer
      end if
      if( denom < 0._r8 ) then
	    denom = dayspy + denom
      end if
      dels = max( min( 1._r8,numer/denom ),0._r8 )

      end subroutine set_sulf_time

      subroutine sulf_interp( lonndx, latndx, pmid, ccm_sulf, ncol )
!-----------------------------------------------------------------------
! 	... Time interpolate sulfatei to current time
!-----------------------------------------------------------------------

      use ppgrid, only : pcols, pver

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)   :: ncol                    ! columns in chunk
      integer, intent(in)   :: lonndx(pcols)           ! longitude indicies in chunk
      integer, intent(in)   :: latndx(pcols)           ! latitude  indicies in chunk
      real(r8), intent(in)  :: pmid(pcols,pver)        ! midpoint pressure ( pascals )
      real(r8), intent(out) :: ccm_sulf(ncol,pver)     ! output sulfate

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer     ::  i, ii, jj, k, ku, kl             ! working indicies
      real(r8)    ::  delp                             ! pressure delta
      real(r8)    ::  pinterp                          ! interpolation pressure
      real(r8)    ::  sulf_press(ncol,sulf_nlevs,2)
      real(r8)    ::  zint_sulf(ncol,pver,2)

      if ( .not. read_sulf ) return

!-----------------------------------------------------------------------
! 	... Form MOZ1 pressures at enclosing months
!-----------------------------------------------------------------------
      do k = 1,sulf_nlevs
         do i = 1,ncol
	    ii = lonndx(i)
	    jj = latndx(i)
	    sulf_press(i,k,1) = sulf_hyam(k)*sulf_P0 + sulf_hybm(k)*psi(ii,jj,last)
	    sulf_press(i,k,2) = sulf_hyam(k)*sulf_P0 + sulf_hybm(k)*psi(ii,jj,next)
         end do
      end do

!-----------------------------------------------------------------------
! 	... Vertical interpolation of images data to model levels
!	    Note: images data only up to 50mb
!-----------------------------------------------------------------------
Column_loop : &
      do i = 1,ncol
	 ii = lonndx(i)
	 jj = latndx(i)
Level_loop : &
	 do k = 1,pver
	    pinterp = pmid(i,k)
	    if( pinterp < maxval( sulf_press(i,1,:2) ) ) then
	       zint_sulf(i,k,1) = 0._r8
	       zint_sulf(i,k,2) = 0._r8
	    else
	       if( pinterp > sulf_press(i,sulf_nlevs,1) ) then
		  zint_sulf(i,k,1) = sulfatei(ii,sulf_nlevs,jj,last)
	       else
	          do ku = 2,sulf_nlevs
		     if( pinterp <= sulf_press(i,ku,1) ) then
			kl = ku - 1
			delp =  log( pinterp/sulf_press(i,kl,1) ) &
                             /  log( sulf_press(i,ku,1)/sulf_press(i,kl,1) )
			zint_sulf(i,k,1) = sulfatei(ii,kl,jj,last) &
                                           + delp * (sulfatei(ii,ku,jj,last) - sulfatei(ii,kl,jj,last))
		        exit
		     end if
	          end do
	       end if
	       if( pinterp > sulf_press(i,sulf_nlevs,2) ) then
		  zint_sulf(i,k,2) = sulfatei(ii,sulf_nlevs,jj,next)
		  cycle
	       else
	          do ku = 2,sulf_nlevs
		     if( pinterp <= sulf_press(i,ku,2) ) then
			kl = ku - 1
			delp = log( pinterp/sulf_press(i,kl,2) ) &
                             / log( sulf_press(i,ku,2)/sulf_press(i,kl,2) )
			zint_sulf(i,k,2) = sulfatei(ii,kl,jj,next) &
                                           + delp * (sulfatei(ii,ku,jj,next) - sulfatei(ii,kl,jj,next))
		        exit
		     end if
	          end do
	       end if
	    end if
	 end do Level_loop
      end do Column_loop

!-----------------------------------------------------------------------
!     	... Linear time interpolation     
!-----------------------------------------------------------------------
      do k = 1,pver
	 ccm_sulf(:,k) = zint_sulf(:,k,1) + dels*(zint_sulf(:,k,2) - zint_sulf(:,k,1))
      end do

      end subroutine sulf_interp

      end module mo_sulf
