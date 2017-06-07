!-------------------------------------------------------------------
! Manages reading and interpolation of prescribed aerosol concentrations.
! This places the concentration fields in the physics buffer so that
! radiation package can access them.
!
! This has been generalized so that the field names in the data files
! and the field names in the physics buffer can be arbitrary.
!
! The prescribed_aero_specifier namelist variable specifies a list of 
! variable names of the concentration fields in the netCDF dataset (ncdf_fld_name)
! and the corresponding names used in the physics buffer:
!
! prescribed_aero_specifier = 'pbuf_name1:ncdf_fld_name1','pbuf_name2:ncdf_fld_name2', ...
!
! If there is no ":" then the specified name is used as both the 
! pbuf_name and ncdf_fld_name
!
! Created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_aero

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld, trfile
  use cam_logfile,  only : iulog

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: prescribed_aero_init
  public :: prescribed_aero_adv
  public :: write_prescribed_aero_restart
  public :: read_prescribed_aero_restart
  public :: has_prescribed_aero
  public :: prescribed_aero_register
  public :: init_prescribed_aero_restart
  public :: prescribed_aero_readnl

  logical :: has_prescribed_aero = .false.

  integer, parameter, public :: N_AERO = 50

  integer :: number_flds

  character(len=256) :: filename = ' '
  character(len=256) :: filelist = ' '
  character(len=256) :: datapath = ' '
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  character(len=32)  :: specifier(N_AERO) = ''

  ! prescribed aerosol names 
  character(len=16), allocatable :: pbuf_names(:)

  integer :: aero_cnt

contains

!-------------------------------------------------------------------
! registers aerosol fields to the phys buffer
!-------------------------------------------------------------------
  subroutine prescribed_aero_register()

    use ppgrid,         only: pver,pcols
    
    use physics_buffer, only : pbuf_add_field, dtype_r8
    integer :: i,idx

    if (has_prescribed_aero) then
       do i = 1,aero_cnt
          call pbuf_add_field(pbuf_names(i),'physpkg',dtype_r8,(/pcols,pver/),idx)
       enddo
    endif

  endsubroutine prescribed_aero_register

!-------------------------------------------------------------------
! parses prescribed_aero_specifier namelist option
!-------------------------------------------------------------------
  subroutine prescribed_aero_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    implicit none

    integer :: ndx, istat, i
    
    if ( has_prescribed_aero ) then
       if ( masterproc ) then
          write(iulog,*) 'aero is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    file%in_pbuf = .true.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)
        
    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if( number_flds < 1 ) then
       if ( masterproc ) then
          write(iulog,*) 'There are no prescribed aerosols'
          write(iulog,*) ' '
       endif
       return
    end if
    
    do i = 1,number_flds
       call addfld(trim(fields(i)%fldnam)//'_D',trim(fields(i)%units), pver, 'A', 'prescribed aero', phys_decomp )
    enddo

  end subroutine prescribed_aero_init

!-------------------------------------------------------------------
! reads namelist options
!-------------------------------------------------------------------
subroutine prescribed_aero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_aero_readnl'

   character(len=32)  :: prescribed_aero_specifier(N_AERO)
   character(len=256) :: prescribed_aero_file
   character(len=256) :: prescribed_aero_filelist
   character(len=256) :: prescribed_aero_datapath
   character(len=32)  :: prescribed_aero_type
   logical            :: prescribed_aero_rmfile
   integer            :: prescribed_aero_cycle_yr
   integer            :: prescribed_aero_fixed_ymd
   integer            :: prescribed_aero_fixed_tod
   integer :: i,k

   namelist /prescribed_aero_nl/ &
      prescribed_aero_specifier, &
      prescribed_aero_file,      &
      prescribed_aero_filelist,  &
      prescribed_aero_datapath,  &
      prescribed_aero_type,      &
      prescribed_aero_rmfile,    &
      prescribed_aero_cycle_yr,  &
      prescribed_aero_fixed_ymd, &
      prescribed_aero_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_aero_specifier= specifier
   prescribed_aero_file     = filename
   prescribed_aero_filelist = filelist
   prescribed_aero_datapath = datapath
   prescribed_aero_type     = datatype
   prescribed_aero_rmfile   = rmv_file
   prescribed_aero_cycle_yr = cycle_yr
   prescribed_aero_fixed_ymd= fixed_ymd
   prescribed_aero_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_aero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_aero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_aero_specifier,len(prescribed_aero_specifier(1))*N_AERO,     mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_file,     len(prescribed_aero_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_filelist, len(prescribed_aero_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_datapath, len(prescribed_aero_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_type,     len(prescribed_aero_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_aero_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_aero_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_aero_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier  = prescribed_aero_specifier
   filename   = prescribed_aero_file
   filelist   = prescribed_aero_filelist
   datapath   = prescribed_aero_datapath
   datatype   = prescribed_aero_type
   rmv_file   = prescribed_aero_rmfile
   cycle_yr   = prescribed_aero_cycle_yr
   fixed_ymd  = prescribed_aero_fixed_ymd
   fixed_tod  = prescribed_aero_fixed_tod

   ! Turn on prescribed aerosols if user has specified an input dataset.
   has_prescribed_aero = len_trim(filename) > 0 

   if ( .not. has_prescribed_aero) return

   ! determine which prescibed aerosols are specified
   aero_cnt = 0
   cnt_loop: do i = 1,N_AERO
      if ( len_trim(specifier(i))==0 ) exit cnt_loop
      aero_cnt = aero_cnt+1
   end do cnt_loop

   has_prescribed_aero = aero_cnt>0
   if ( .not. has_prescribed_aero) return

   allocate(pbuf_names(aero_cnt))

   do i = 1,aero_cnt
      k = scan( specifier(i),':')
      if (k>1) then
         pbuf_names(i) = trim(specifier(i)(:k-1)) 
      else
         pbuf_names(i) = trim(specifier(i))
      endif
   enddo

end subroutine prescribed_aero_readnl

!-------------------------------------------------------------------
! advances the aerosol fields to the current time step
!-------------------------------------------------------------------
  subroutine prescribed_aero_adv( state, pbuf2d )

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer :: c,ncol, i
    real(r8),pointer :: outdata(:,:)

    if( .not. has_prescribed_aero ) return

    call advance_trcdata( fields, file, state, pbuf2d )
    
    ! set the tracer fields with the correct units
    do i = 1,number_flds

!$OMP PARALLEL DO PRIVATE (C, NCOL, OUTDATA,PBUF_CHNK)
       do c = begchunk,endchunk
          ncol = state(c)%ncol
          pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
          call pbuf_get_field(pbuf_chnk, fields(i)%pbuf_ndx, outdata)
          call outfld( trim(fields(i)%fldnam)//'_D', outdata(:ncol,:), ncol, state(c)%lchnk )
       enddo
    enddo

  end subroutine prescribed_aero_adv

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine init_prescribed_aero_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'prescribed_aero', piofile, file )

  end subroutine init_prescribed_aero_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_aero_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_aero_restart

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine read_prescribed_aero_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'prescribed_aero', piofile, file )

  end subroutine read_prescribed_aero_restart

end module prescribed_aero
