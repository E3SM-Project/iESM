!=================================================================================
! input data necessary to drive radiation 
! Francis Vitt -- Created 15 Dec 2009
!=================================================================================
module rad_data_input

  use abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use shr_kind_mod, only : r8 => shr_kind_r8,r4 => shr_kind_r4, cl=>shr_kind_cl, cs=>shr_kind_cs
  use ppgrid,       only : pcols, pver, pverp, begchunk, endchunk
  use cam_logfile,  only : iulog
  use pio,          only : file_desc_t

  use radiation_data

  use physics_types, only: physics_state
  use physics_buffer,only: physics_buffer_desc, pbuf_get_chunk
  use radconstants,     only: nradgas, gaslist

  implicit none

  save

  private

  public :: init_rad_data_input
  public :: get_rad_data_input
  public :: close_rad_data_input

  integer, public :: ntimes
  integer, public :: nlons
  integer, public :: nlats
  integer, public :: nlevs
  integer, public :: nilevs

  integer, public, allocatable :: dates(:)
  integer, public, allocatable :: secs(:)

  integer, public :: dtime
  integer, public :: ref_ymd
  integer, public :: ref_tod

  character(len=cs), public :: calendar  ! Calendar type
  type(file_desc_t) :: piofile

  type data_ptr_3d
     real(r8), pointer :: array(:,:)
  endtype data_ptr_3d
  type data_ptr_2d
     real(r8), pointer :: array(:)
  endtype data_ptr_2d

  type data_iptr_2d
     integer, pointer :: array(:)
  endtype data_iptr_2d

  ! rad constituents names

  integer :: ngas, naer
  character(len=64), allocatable :: gasnames(:)
  character(len=64), allocatable :: aernames(:)

contains

  !=================================================================================
  !=================================================================================
  subroutine init_rad_data_input( infile )

    use ioFileMod, only: getfil
    use cam_pio_utils, only: cam_pio_openfile
    use pio,       only: PIO_NOCLOBBER, pio_inq_dimid, pio_inq_dimlen
    use pio,       only: pio_inq_varid, pio_get_var, pio_get_att

    use shr_const_mod, only: SHR_CONST_CDAY
    use shr_string_mod,only: shr_string_toUpper

    use pmgrid, only: plon,plat,plev !,plevp,plnlv
    use rad_constituents, only: rad_cnst_get_info

    implicit none

    character(len=*), intent(in) :: infile


    character(len=cl) :: filen
    integer :: id, ierr
    integer :: nvars, vid, ndims

    character(len=cs) :: vname
    character(len=cs) :: time_units
    character(len=10) :: ref_date_str
    character(len=4) :: ref_year_str
    character(len=2) :: ref_mon_str
    character(len=2) :: ref_day_str
    character(len=2) :: ref_hr_str
    character(len=2) :: ref_min_str
    character(len=2) :: ref_sec_str
    character(len=8) :: ref_time_str
    integer :: ref_year, ref_mon, ref_day, ref_hr, ref_min, ref_sec

    !
    ! open file and get fileid
    !
    call getfil( infile, filen, 0 )
    call cam_pio_openfile( piofile, filen, PIO_NOCLOBBER)

    if(masterproc) write(iulog,*)'init_rad_data_input: opened: ',trim(filen)

    ierr = pio_inq_dimid( piofile, 'time', id )
    ierr = pio_inq_dimlen( piofile, id, ntimes )

    ierr = pio_inq_dimid( piofile, 'lon', id )
    ierr = pio_inq_dimlen( piofile, id, nlons )

    ierr = pio_inq_dimid( piofile, 'lat', id )
    ierr = pio_inq_dimlen( piofile, id, nlats )

    ierr = pio_inq_dimid( piofile, 'lev', id )
    ierr = pio_inq_dimlen( piofile, id, nlevs )

    ierr = pio_inq_dimid( piofile, 'ilev', id )
    ierr = pio_inq_dimlen( piofile, id, nilevs )

    allocate( dates(ntimes), secs(ntimes) )

    ierr = pio_inq_varid( piofile, 'date',  id  )
    ierr = pio_get_var( piofile, id, dates )

    ierr = pio_inq_varid( piofile, 'datesec',  id  )
    ierr = pio_get_var( piofile, id, secs )

    ierr = pio_inq_varid( piofile, 'time',  id  )
    ierr = pio_get_att( piofile, id, 'units', time_units )
    ierr = pio_get_att( piofile, id, 'calendar', calendar )

    calendar = shr_string_toUpper(calendar)

    ref_date_str = time_units(12:21)
    ref_year_str = ref_date_str(1:4)
    ref_mon_str = ref_date_str(6:7)
    ref_day_str = ref_date_str(9:10)

    read(ref_year_str,fmt='(I4)') ref_year
    read(ref_mon_str, fmt='(I2)') ref_mon
    read(ref_day_str, fmt='(I2)') ref_day

    ref_ymd = ref_year*10000 + ref_mon*100 + ref_day

    ref_time_str = time_units(23:)
    ref_hr_str = ref_time_str(1:2)
    ref_min_str = ref_time_str(4:5)
    ref_sec_str = ref_time_str(7:8)

    read(ref_hr_str, fmt='(I2)') ref_hr
    read(ref_min_str,fmt='(I2)') ref_min
    read(ref_sec_str,fmt='(I2)') ref_sec

    ref_tod = ref_hr*3600 + ref_min*60 + ref_sec

    ierr = pio_inq_varid( piofile, 'mdt',  id  )
    ierr = pio_get_var( piofile, id, dtime )

    if (plon /= nlons) then
       call endrun('init_rad_data_input plon /= nlons')
    endif
    if (plat /= nlats) then
       call endrun('init_rad_data_input plat /= nlats')
    endif
    if (plev /= nlevs) then
       call endrun('init_rad_data_input plev /= nlevs')
    endif

    call rad_cnst_get_info( 0, ngas=ngas, naero=naer )
    allocate(gasnames(ngas))
    allocate(aernames(naer))

    call rad_cnst_get_info( 0, gasnames=gasnames, aernames=aernames )

  end subroutine init_rad_data_input

  !=================================================================================
  !=================================================================================
  subroutine get_rad_data_input( phys_state, pbuf2d, cam_in, landm, recno )

    use physics_types,    only: physics_state
    use camsrfexch,       only: cam_in_t
    use physics_buffer,   only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use constituents,     only: cnst_get_ind

    implicit none

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_state), target, intent(inout) :: phys_state(begchunk:endchunk)
    type(cam_in_t),      target, intent(inout) :: cam_in(begchunk:endchunk)
    real(r8),            target, intent(inout) :: landm(pcols,begchunk:endchunk)
    integer,             intent(in)    :: recno

    ! local vars

    type(physics_buffer_desc), pointer :: pbuf(:)
    type(data_ptr_3d) :: cld_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: concld_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: cldfsnow_ptrs(begchunk:endchunk)

    type(data_ptr_3d) :: rel_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: rei_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: dei_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: des_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: mu_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: lambdac_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: iciwp_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: iclwp_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: icswp_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: rel_fn_ptrs(begchunk:endchunk)

    type(data_ptr_3d) :: watvap_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: watliq_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: watice_ptrs(begchunk:endchunk)

    type(data_ptr_3d) :: zi_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: pint_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: lnpint_ptrs(begchunk:endchunk)

    type(data_ptr_3d) :: temp_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: pdel_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: pdeldry_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: lnpmid_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: pmid_ptrs(begchunk:endchunk)

    type(data_ptr_2d) :: landm_ptrs(begchunk:endchunk)

    type(data_ptr_2d) :: lndfrac_ptrs(begchunk:endchunk)
    type(data_ptr_2d) :: icefrac_ptrs(begchunk:endchunk)
    type(data_ptr_2d) :: snowhland_ptrs(begchunk:endchunk)
    type(data_ptr_2d) :: asdir_ptrs(begchunk:endchunk)
    type(data_ptr_2d) :: asdif_ptrs(begchunk:endchunk)
    type(data_ptr_2d) :: aldir_ptrs(begchunk:endchunk)
    type(data_ptr_2d) :: aldif_ptrs(begchunk:endchunk)

    type(data_iptr_2d) :: nmxrgn_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: pmxrgn_ptrs(begchunk:endchunk)

    type(data_ptr_3d) :: cldemis_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: cldtau_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: cicewp_ptrs(begchunk:endchunk)
    type(data_ptr_3d) :: cliqwp_ptrs(begchunk:endchunk)

    type(data_ptr_2d) :: lwup_ptrs(begchunk:endchunk)
    type(data_ptr_2d) :: ts_ptrs(begchunk:endchunk)

    integer :: i, c, ncol, itim

    integer :: ixcldice   ! cloud ice water index
    integer :: ixcldliq   ! cloud liquid water index

    ! get index of (liquid+ice) cloud water
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

    ! phys buffer time index
    itim = pbuf_old_tim_idx()

    ! setup the data pointers
    do c=begchunk,endchunk

       pbuf => pbuf_get_chunk(pbuf2d, c)

       watvap_ptrs (c)%array => phys_state(c)%q(:,:,1)
       watliq_ptrs (c)%array => phys_state(c)%q(:,:,ixcldliq)
       watice_ptrs (c)%array => phys_state(c)%q(:,:,ixcldice)

       zi_ptrs     (c)%array => phys_state(c)%zi
       pint_ptrs   (c)%array => phys_state(c)%pint
       lnpint_ptrs (c)%array => phys_state(c)%lnpint

       temp_ptrs   (c)%array => phys_state(c)%t
       pdel_ptrs   (c)%array => phys_state(c)%pdel
       pdeldry_ptrs(c)%array => phys_state(c)%pdeldry
       lnpmid_ptrs (c)%array => phys_state(c)%lnpmid
       pmid_ptrs   (c)%array => phys_state(c)%pmid

       landm_ptrs    (c)%array => landm(:,c) 
       lndfrac_ptrs  (c)%array => cam_in(c)%landfrac
       icefrac_ptrs  (c)%array => cam_in(c)%icefrac
       snowhland_ptrs(c)%array => cam_in(c)%snowhland
       asdir_ptrs    (c)%array => cam_in(c)%asdir
       asdif_ptrs    (c)%array => cam_in(c)%asdif
       aldir_ptrs    (c)%array => cam_in(c)%aldir
       aldif_ptrs    (c)%array => cam_in(c)%aldif
       lwup_ptrs     (c)%array => cam_in(c)%lwup
       ts_ptrs       (c)%array => cam_in(c)%ts

       call pbuf_get_field(pbuf, cld_ifld, cld_ptrs   (c)%array, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
       call pbuf_get_field(pbuf, concld_ifld, concld_ptrs(c)%array, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
       call pbuf_get_field(pbuf, rel_ifld,   rel_ptrs(c)%array )
       call pbuf_get_field(pbuf, rei_ifld,   rei_ptrs(c)%array )

       if (mg_microphys) then
          call pbuf_get_field(pbuf, dei_ifld,   dei_ptrs(c)%array )
          call pbuf_get_field(pbuf, des_ifld,   des_ptrs(c)%array )
          call pbuf_get_field(pbuf, mu_ifld,    mu_ptrs     (c)%array )
          call pbuf_get_field(pbuf, lambdac_ifld, lambdac_ptrs(c)%array )
          call pbuf_get_field(pbuf, iciwp_ifld, iciwp_ptrs  (c)%array )
          call pbuf_get_field(pbuf, iclwp_ifld, iclwp_ptrs  (c)%array )
          call pbuf_get_field(pbuf, icswp_ifld, icswp_ptrs  (c)%array )
          call pbuf_get_field(pbuf, rel_fn_ifld, rel_fn_ptrs (c)%array )
          call pbuf_get_field(pbuf, cldfsnow_ifld, cldfsnow_ptrs(c)%array, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
       else
          call pbuf_get_field(pbuf, nmxrgn_ifld,   nmxrgn_ptrs(c)%array )
          call pbuf_get_field(pbuf, pmxrgn_ifld,   pmxrgn_ptrs(c)%array )
          call pbuf_get_field(pbuf, cldemis_ifld,  cldemis_ptrs(c)%array )
          call pbuf_get_field(pbuf, cldtau_ifld,   cldtau_ptrs(c)%array )
          call pbuf_get_field(pbuf, cicewp_ifld,   cicewp_ptrs(c)%array )
          call pbuf_get_field(pbuf, cliqwp_ifld,   cliqwp_ptrs(c)%array )
       endif

    enddo

    ! get the 2D data

    call get_data2d( landm_fldn,  recno, landm_ptrs )
    call get_data2d( lndfrc_fldn, recno, lndfrac_ptrs )
    call get_data2d( icefrc_fldn, recno, icefrac_ptrs )
    call get_data2d( snowh_fldn,  recno, snowhland_ptrs )
    call get_data2d( asdir_fldn,  recno, asdir_ptrs )
    call get_data2d( asdif_fldn,  recno, asdif_ptrs )
    call get_data2d( aldir_fldn,  recno, aldir_ptrs )
    call get_data2d( aldif_fldn,  recno, aldif_ptrs )
    call get_data2d( lwup_fldn,   recno, lwup_ptrs )
    call get_data2d( ts_fldn,     recno, ts_ptrs )

    ! get the 3D data

    call get_data3d( cld_fldn,    'lev', pver, recno, cld_ptrs )
    call get_data3d( concld_fldn, 'lev', pver, recno, concld_ptrs )
    call get_data3d( rel_fldn,    'lev', pver, recno, rel_ptrs )
    call get_data3d( rei_fldn,    'lev', pver, recno, rei_ptrs )

    if (mg_microphys) then
       call get_data3d( dei_fldn,    'lev', pver, recno, dei_ptrs )
       call get_data3d( des_fldn,    'lev', pver, recno, des_ptrs )
       call get_data3d( mu_fldn,     'lev', pver, recno, mu_ptrs )
       call get_data3d( lambdac_fldn,'lev', pver, recno, lambdac_ptrs )
       call get_data3d( iciwp_fldn,  'lev', pver, recno, iciwp_ptrs )
       call get_data3d( iclwp_fldn,  'lev', pver, recno, iclwp_ptrs )
       call get_data3d( icswp_fldn,  'lev', pver, recno, icswp_ptrs )
       call get_data3d( rel_fn_fldn, 'lev', pver, recno, rel_fn_ptrs )
       call get_data3d( cldfsnow_fldn,'lev',pver, recno, cldfsnow_ptrs )
    else
       call get_idata2d( nmxrgn_fldn, recno, nmxrgn_ptrs )
       call get_data3d( pmxrgn_fldn, 'ilev',pverp,recno, pmxrgn_ptrs )
       call get_data3d( cldemis_fldn,'lev', pver, recno, cldemis_ptrs )
       call get_data3d( cldtau_fldn, 'lev', pver, recno, cldtau_ptrs )
       call get_data3d( cicewp_fldn, 'lev', pver, recno, cicewp_ptrs )
       call get_data3d( cliqwp_fldn, 'lev', pver, recno, cliqwp_ptrs )
    endif

    call get_data3d( watvap_fldn, 'lev', pver, recno, watvap_ptrs )
    call get_data3d( watliq_fldn, 'lev', pver, recno, watliq_ptrs )
    call get_data3d( watice_fldn, 'lev', pver, recno, watice_ptrs )
    call get_data3d( temp_fldn,   'lev', pver, recno, temp_ptrs )
    call get_data3d( pdel_fldn,   'lev', pver, recno, pdel_ptrs )
    call get_data3d( pdeldry_fldn,'lev', pver, recno, pdeldry_ptrs )
    call get_data3d( pmid_fldn,   'lev', pver, recno, pmid_ptrs )

    call get_data3d( zint_fldn,   'ilev',pverp,recno, zi_ptrs )
    call get_data3d( pint_fldn,   'ilev',pverp,recno, pint_ptrs )

    !$OMP PARALLEL DO PRIVATE (c,i,ncol)
    do c=begchunk,endchunk
       ncol = phys_state(c)%ncol
       if (all(pmid_ptrs(c)%array(:ncol,:) > 0._r8 )) then
          do i=1,ncol
             lnpmid_ptrs(c)%array(i,:) = log(pmid_ptrs(c)%array(i,:))
             lnpint_ptrs(c)%array(i,:) = log(pint_ptrs(c)%array(i,:))
          enddo
       endif
    enddo

    call get_rad_cnst_data(phys_state, pbuf2d, recno)

  end subroutine get_rad_data_input

  !================================================================================================
  !================================================================================================
  subroutine close_rad_data_input
    use pio, only: pio_closefile
    implicit none

    deallocate( dates, secs )

    call pio_closefile( piofile )

  end subroutine close_rad_data_input

  !================================================================================================
  ! Private routines
  !================================================================================================

  !================================================================================================
  !================================================================================================
  subroutine get_data3d(infld_name, lev_name, nlev, recno, chunk_ptrs)

    character(len=*),    intent(in)    :: infld_name
    character(len=*),    intent(in)    :: lev_name
    integer,             intent(in)    :: nlev
    integer,             intent(in)    :: recno
    type(data_ptr_3d),   intent(inout) :: chunk_ptrs(begchunk:endchunk)

    real(r8), allocatable :: data (:,:,:)

    integer :: c, ncol

    allocate( data (pcols, nlev,  begchunk:endchunk) )

    data = read_data_3d( infld_name, lev_name, nlev, recno )
    do c=begchunk,endchunk
       chunk_ptrs(c)%array(:,:) = data(:,:,c)
    enddo

    deallocate( data )

  end subroutine get_data3d

  !================================================================================================
  !================================================================================================
  subroutine get_data2d(infld_name, recno, chunk_ptrs)

    character(len=*),    intent(in)    :: infld_name
    integer,             intent(in)    :: recno
    type(data_ptr_2d),   intent(inout) :: chunk_ptrs(begchunk:endchunk)

    real(r8), allocatable :: data (:,:)

    integer :: c, ncol

    allocate( data (pcols,  begchunk:endchunk) )

    data = read_data_2d( infld_name, recno )
    do c=begchunk,endchunk
       chunk_ptrs(c)%array(:) = data(:,c)
    enddo

    deallocate( data )

  end subroutine get_data2d

  !================================================================================================
  !================================================================================================
  subroutine get_idata2d(infld_name, recno, chunk_ptrs)

    character(len=*),    intent(in)    :: infld_name
    integer,             intent(in)    :: recno
    type(data_iptr_2d),   intent(inout) :: chunk_ptrs(begchunk:endchunk)

    real(r8), allocatable :: data (:,:)

    integer :: c, ncol

    allocate( data (pcols,  begchunk:endchunk) )

    data = read_data_2d( infld_name, recno )
    do c=begchunk,endchunk
       chunk_ptrs(c)%array(:) = int(data(:,c))
    enddo

    deallocate( data )

  end subroutine get_idata2d

  !================================================================================================
  !================================================================================================
  subroutine get_rad_cnst_data(state, pbuf2d, recno)

    use physics_types,    only: physics_state
    use physics_buffer,   only: physics_buffer_desc

    implicit none

    ! Arguments
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_state), intent(inout) :: state(begchunk:endchunk)
    integer,             intent(in)    :: recno

    ! Local variables
    integer :: i, ncol
    integer :: idx
    character(len=32) :: name

    !-----------------------------------------------------------------------------

    ! read in mixing ratio of rad constituents 

    do i = 1,ngas
       call read_rad_gas_data(gasnames(i), i, state, pbuf2d, recno )
    enddo

    do i = 1,naer
       call read_rad_aer_data(aernames(i), i, state, pbuf2d, recno )
    enddo

  end subroutine get_rad_cnst_data

  !=================================================================================
  !=================================================================================
  subroutine read_rad_gas_data(name, idx, state, pbuf2d, recno )
    use rad_constituents, only: rad_cnst_get_gas

    character(len=32),   intent(in)    :: name
    integer,             intent(in)    :: idx
    type(physics_state),target, intent(inout) :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
    integer,             intent(in)    :: recno

    type(physics_buffer_desc), pointer, dimension(:) :: phys_buffer_chunk
    character(len=32) :: radname
    integer            :: c, ncol
    real(r8)          :: mass(pcols,pver,begchunk:endchunk)
    real(r8), pointer :: mmr(:,:)

    radname = 'rad_'//trim(name)
    mass = read_data_3d( radname, 'lev', pver, recno )

    do c = begchunk,endchunk
       ncol = state(c)%ncol
       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
       call rad_cnst_get_gas(0, gaslist(idx), state(c), phys_buffer_chunk, mmr)
       mmr(:ncol,:) = mass(:ncol,:,c)
    enddo

  end subroutine read_rad_gas_data

  !=================================================================================
  !=================================================================================
  subroutine read_rad_aer_data(name, idx, state, pbuf2d, recno )
    use rad_constituents, only: rad_cnst_get_aer_mmr

    character(len=32),   intent(in)    :: name
    integer,             intent(in)    :: idx
    type(physics_state),target, intent(inout) :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
    integer,             intent(in)    :: recno

    type(physics_buffer_desc), pointer, dimension(:) :: phys_buffer_chunk
    character(len=32) :: radname
    integer :: c, ncol
    real(r8)          :: mass(pcols,pver,begchunk:endchunk)
    real(r8), pointer :: mmr(:,:)

    radname = 'rad_'//trim(name)
    mass = read_data_3d( radname, 'lev', pver, recno )

    do c = begchunk,endchunk
       ncol = state(c)%ncol
       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
       call rad_cnst_get_aer_mmr(0, idx, state(c), phys_buffer_chunk, mmr)
       mmr(:ncol,:) = mass(:ncol,:,c)
    enddo

  end subroutine read_rad_aer_data

  !=================================================================================
  !=================================================================================
  function read_data_2d( fldname, recno )
    use ncdio_atm,        only: infld

    implicit none

    character(len=*), intent(in) :: fldname
    integer,          intent(in) :: recno

    logical  :: found
    real(r8) :: read_data_2d(pcols,begchunk:endchunk)

    call infld( fldname, piofile, 'lon','lat', 1,pcols, begchunk,endchunk, &
         read_data_2d, found, grid_map='PHYS',timelevel=recno)

    if (.not.found) then
       call endrun('rad_data_input.read_data_2d: did not find '// trim(fldname))
    endif

  endfunction read_data_2d

  !=================================================================================
  !=================================================================================
  function read_data_3d( fldname, vertname, vertsize, recno )
    use ncdio_atm,        only: infld

    implicit none

    character(len=*), intent(in) :: fldname
    character(len=*), intent(in) :: vertname
    integer,          intent(in) :: vertsize
    integer,          intent(in) :: recno

    logical  :: found
    real(r8) :: read_data_3d(pcols,vertsize,begchunk:endchunk)

    call infld( fldname, piofile, 'lon',vertname,'lat', 1,pcols, 1,vertsize, begchunk,endchunk, &
         read_data_3d, found, grid_map='PHYS',timelevel=recno)

    if (.not.found) then
       call endrun('rad_data_input.read_data_3d: did not find '// trim(fldname))
    endif

  endfunction read_data_3d

end module rad_data_input
