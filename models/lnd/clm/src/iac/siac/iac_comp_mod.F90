
module iac_comp_mod
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: iac_comp_mod
!
!  Interface of the integrated assessment component in CCSM
!
! !DESCRIPTION:
!
! !USES:
  use iac_fields_mod

  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public :: iac_init_mod               ! clm initialization
  public :: iac_run_mod                ! clm run phase
  public :: iac_final_mod              ! clm finalization/cleanup

! !PUBLIC DATA MEMBERS: None

! !REVISION HISTORY:
! Author: T Craig


! !PRIVATE DATA MEMBERS:

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: iac_init_mod

! !INTERFACE:
  subroutine iac_init_mod( EClock, cdata, iaci, iaco)

! !DESCRIPTION:
! Initialize interface for iac

! !USES:
    implicit none

! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real*8, pointer :: iaci(:,:)
    real*8, pointer :: iaco(:,:)

! !LOCAL VARIABLES:


! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

!  allocate(iaci(iac_iaci_nflds,iac_data_size))
!  allocate(iaco(iac_iaco_nflds,iac_data_size))
!  iaci = 0.0
!  iaco = 0.0

  cdata%l(iac_cdatal_iac_present) = .false.
  cdata%l(iac_cdatal_iac_prognostic) = .false.

  end subroutine iac_init_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: iac_run_mod

! !INTERFACE:
  subroutine iac_run_mod( EClock, cdata, iaci, iaco)

! !DESCRIPTION:
! Run interface for iac

! !USES:
    implicit none

! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real*8, pointer :: iaci(:,:)
    real*8, pointer :: iaco(:,:)

! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='(iac_run_mod)'


! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

  ! nothing to do 
  end subroutine iac_run_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: iac_final_mod

! !INTERFACE:
  subroutine iac_final_mod( )

! !DESCRIPTION:
! Finalize iac model

!------------------------------------------------------------------------------

   implicit none
! !ARGUMENTS:

! !REVISION HISTORY:
! Author: T Craig

!EOP
!---------------------------------------------------------------------------

  ! nothing to do 
  end subroutine iac_final_mod

!====================================================================================

end module iac_comp_mod
