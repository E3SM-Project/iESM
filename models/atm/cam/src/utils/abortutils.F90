
module abortutils

   implicit none
   private
   save

   public :: endrun

CONTAINS

subroutine endrun (msg)
!-----------------------------------------------------------------------
! Purpose:
!
! Abort the model for abnormal termination
!
!-----------------------------------------------------------------------
#if (defined SPMD) 
   use mpishorthand, only: MPI_COMM_WORLD
#endif
#if defined(Darwin)
   use ifport
#endif
   use shr_sys_mod,  only: shr_sys_flush
   use cam_logfile,  only: iulog

   ! Arguments
   character(len=*), intent(in), optional :: msg    ! string to be printed

   ! Local variables
   integer :: ierr
   !-----------------------------------------------------------------------

   if (present (msg)) then
      write(iulog,*)'ENDRUN:', msg
   else
      write(iulog,*)'ENDRUN: called without a message string'
   end if

#if defined(NEC_SX)
   call mesput("ENDRUN", len("ENDRUN"), 1)
#elif defined(AIX)
   close(5)    ! needed to prevent batch jobs from hanging in xl__trbk
#if !defined(BGL) && !defined(BGP)
   call xl__trbk()
#endif
#endif

   call shr_sys_flush( iulog )

#if (defined SPMD) 
   ! passing an argument of 1 to mpi_abort will lead to a STOPALL output 
   ! error code of 257
   call mpi_abort (MPI_COMM_WORLD, 1, ierr)  
#else

#if defined(Darwin)
   call abort("Aborting ...")
#else
   call abort
#endif

#endif

end subroutine endrun
end module abortutils
