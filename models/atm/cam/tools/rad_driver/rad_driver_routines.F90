!=================================================================================
! utility routines for the offline radiation driver 
! Francis Vitt -- Created 15 Dec 2009
!=================================================================================
module rad_driver_routines

  use shr_kind_mod,    only: r8=>SHR_KIND_R8, cl=>SHR_KIND_CL, cs=>SHR_KIND_CS

  implicit none
  private

  public :: rad_driver_init
  public :: rad_driver_final

  ! private data:

  character(len=cl) :: rad_drv_infile = 'rad_drv_infile'
  character(len=cs) :: rad_drv_case = 'RAD_DRIVER'

contains

  subroutine rad_driver_init()
    use filenames,       only: caseid
    use rad_data_input,  only: init_rad_data_input

    implicit none

    !
    ! Read namelist
    !
    call rad_driver_readnl('atm_in')
    caseid = rad_drv_case

    call init_rad_data_input( rad_drv_infile )

  end subroutine rad_driver_init

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

  subroutine rad_driver_final
    use shr_mpi_mod,     only: shr_mpi_finalize
    use rad_data_input,  only: close_rad_data_input
    implicit none
  
    call close_rad_data_input()
    call shr_mpi_finalize('rad_driver_final')

  endsubroutine rad_driver_final

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

  subroutine init_clock()
    use time_manager, only: timemgr_init, tmgr_dtime=>dtime

    use rad_data_input,  only: calendar, dtime, ref_ymd, ref_tod, data_dtime=>dtime, dates, secs, ntimes

    implicit none

    logical :: perpetual_run    ! If in perpetual mode or not
    integer :: perpetual_ymd    ! Perpetual date (YYYYMMDD)
    integer :: start_ymd        ! Start date (YYYYMMDD)
    integer :: start_tod        ! Start time of day (sec)
    integer :: stop_ymd         ! Stop date (YYYYMMDD)
    integer :: stop_tod         ! Stop time of day (sec)

    start_ymd = dates(1)
    start_tod = secs(1)
    stop_ymd = dates(ntimes)
    stop_tod = secs(ntimes)+data_dtime
    perpetual_run =  .false.
    perpetual_ymd =  0

    if (trim(calendar)=='NOLEAP') calendar = 'NO_LEAP'

    call timemgr_init( &
         calendar_in=calendar, &
         start_ymd=start_ymd, &
         start_tod=start_tod, &
         ref_ymd=ref_ymd, &
         ref_tod=ref_tod, &
         stop_ymd=stop_ymd, &
         stop_tod=stop_tod, &
         perpetual_run=perpetual_run, &
         perpetual_ymd=perpetual_ymd  &
         )

     tmgr_dtime = dtime
  endsubroutine init_clock

!--------------------------------------------------------------------------------

  subroutine rad_driver_readnl( nlfile )
    use namelist_utils,only: find_group_name
    use units,         only: getunit, freeunit
    use abortutils,    only: endrun

    implicit none


    ! arguments
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars
    integer :: unitn, ierr

    namelist /rad_drv_nl/ rad_drv_infile, rad_drv_case

    unitn = getunit()
    open( unitn, file=trim(nlfile), status='old' )
    call find_group_name(unitn, 'rad_drv_nl', status=ierr)
    if (ierr == 0) then
       read(unitn, rad_drv_nl, iostat=ierr)
       if (ierr /= 0) then
          call endrun('rad_driver_readnl: ERROR reading namelist')
       end if
    end if
    close(unitn)
    call freeunit(unitn)

  endsubroutine rad_driver_readnl

endmodule rad_driver_routines
