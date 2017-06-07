module error_messages

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! General purpose routines for issuing error messages.
   ! 
   ! Author: B. Eaton
   ! 
   !----------------------------------------------------------------------- 
   use abortutils,  only: endrun
   use cam_logfile, only: iulog

   implicit none
   save
   private
   public :: &
      alloc_err,      &! Issue error message after non-zero return from an allocate statement.
      handle_err,     &! Issue error message after non-zero return from anything
      handle_ncerr     ! Handle error returns from netCDF library procedures.

!##############################################################################
contains
!##############################################################################

   subroutine alloc_err( istat, routine, name, nelem )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Issue error message after non-zero return from an allocate statement.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      integer, intent(in) ::&
         istat           ! status from allocate statement
      character(len=*), intent(in) ::&
         routine,       &! routine that called allocate
         name            ! name of array
      integer, intent(in) ::&
         nelem           ! number of elements attempted to allocate
      !-----------------------------------------------------------------------

      if ( istat .ne. 0 ) then
         write(iulog,*)'ERROR trying to allocate memory in routine: ' &
                   //trim(routine)
         write(iulog,*)'  Variable name: '//trim(name)
         write(iulog,*)'  Number of elements: ',nelem
         call endrun ('ALLOC_ERR')
      end if

      return

   end subroutine alloc_err

!##############################################################################

   subroutine handle_err( istat, msg )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Issue error message after non-zero return from anything.
      !
      ! Author: T. Henderson
      !----------------------------------------------------------------------- 

      integer,          intent(in) :: istat  ! status, zero = "no error"
      character(len=*), intent(in) :: msg    ! error message to print
      !-----------------------------------------------------------------------

      if ( istat .ne. 0 ) then
         call endrun (trim(msg))
      end if

      return

   end subroutine handle_err

!##############################################################################

   subroutine handle_ncerr( ret, mes, line )
      
      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Check netCDF library function return code.  If error detected 
      ! issue error message then abort.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

!-----------------------------------------------------------------------
     use netcdf
!-----------------------------------------------------------------------

      integer, intent(in) ::&
         ret                 ! return code from netCDF library routine
      character(len=*), intent(in) ::&
         mes                 ! message to be printed if error detected
      integer, intent(in), optional :: line
      !-----------------------------------------------------------------------

      if ( ret .ne. NF90_NOERR ) then
         if(present(line)) then
            write(iulog,*) mes, line
         else	
            write(iulog,*) mes
         end if
         write(iulog,*) nf90_strerror( ret )
         call endrun ('HANDLE_NCERR')
      endif

      return

   end subroutine handle_ncerr

!##############################################################################

end module error_messages
