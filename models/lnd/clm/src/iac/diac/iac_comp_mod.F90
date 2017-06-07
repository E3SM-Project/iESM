
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
  use glm2iac_mod
  use shr_cal_mod
  use shr_file_mod
  use shr_sys_mod
  use netcdf

  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public :: iac_init_mod               ! clm initialization
  public :: iac_run_mod                ! clm run phase
  public :: iac_final_mod              ! clm finalization/cleanup

! !PUBLIC DATA MEMBERS: None

  integer, parameter :: glm_data_size   = iac_glm_nx*iac_glm_ny  ! should be set by glm

! !REVISION HISTORY:
! Author: T Craig


! !PRIVATE DATA MEMBERS:

  integer,save :: iulog
  character(len=512) :: glm2iac_glmofile = 'unknown'
  real*8, pointer :: glmi(:,:)
  real*8, pointer :: glmo(:,:)

  namelist /iacnml/   &
       glm2iac_glmofile

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
    character(len=*),parameter :: subname='(iac_init_mod)'
    integer :: iac_data_size
    integer :: nunit, ier
    character(len=128) :: casename


! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

  nunit = shr_file_getUnit()
  open(nunit,file="iac_in",status="old",action="read")
  read(nunit, iacnml, iostat=ier)
  if (ier /= 0) then
     write(iulog,*)'error: iacnml namelist input resulted in error code ',ier
     call shr_sys_abort(subname//' ERROR: iacnml error')
  endif
  close(nunit)
  call shr_file_freeUnit(nunit)

  cdata%l(iac_cdatal_iac_present) = .true.
  cdata%l(iac_cdatal_iac_prognostic) = .false.

  call iac_fields_init

  casename = trim(cdata%c(iac_cdatac_casename))
  iac_data_size = cdata%i(iac_cdatai_iac_size)
  allocate(iaci(iac_iaci_nflds,iac_data_size))
  allocate(iaco(iac_iaco_nflds,iac_data_size))
  iaci = iac_spval
  iaco = iac_spval
!
! This is from glm_init_mod - need a 'data' glm
!
  cdata%i(iac_cdatai_glm_nx) = iac_glm_nx
  cdata%i(iac_cdatai_glm_ny) = iac_glm_ny
  cdata%i(iac_cdatai_glm_size) = iac_glm_nx * iac_glm_ny
  
  allocate(glmi(iac_glmi_nflds,glm_data_size))
  allocate(glmo(iac_glmo_nflds,glm_data_size))
  
  glmi = iac_spval
  glmo = iac_spval

  call glm2iac_init_mod (EClock, cdata, glmo,  iaco)

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
    character(len=256) :: hfile,vname
    integer :: i,j,ij,h,k
    integer :: yyyy,mm,dd
    integer :: iacymd,iactod,iacymd_hold
    integer :: ncidglmo,dimid,varid,nmode,n,ierr,ncid
    integer :: dimidglm(2)
    integer :: start3(3),count3(3)
    real*8, allocatable :: array3(:,:,:)
    real*8, allocatable :: array2(:,:)
    character(len=128) :: casename


! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

  iulog = cdata%i(iac_cdatai_logunit)

  ! iac time
  iacymd = EClock(iac_EClock_ymd)
  iactod = EClock(iac_Eclock_tod)
  call shr_cal_date2ymd(iacymd,yyyy,mm,dd)

  ! compute "alarms" 0 = off, 1 = on
  EClock(iac_EClock_Agcam) = 0
  EClock(iac_EClock_Aglm)  = 0
  EClock(iac_EClock_AclmC) = 0
  EClock(iac_EClock_Agcamsetden) = 0

  if (iactod == EClock(iac_EClock_dt)) then   ! first timestep of day
     if (dd==1 .and. mm==1) then
        EClock(iac_EClock_Aglm)  = 1   ! every year
     end if
  endif

#ifdef DEBUG
  write(iulog,*) trim(subname),'current model date1 ',iacymd,iactod
  write(iulog,*) trim(subname),'current model date2 ',yyyy,mm,dd
  write(iulog,*) trim(subname),'current model alarm ', &
       EClock(iac_Eclock_Agcam),EClock(iac_Eclock_Aglm),EClock(iac_Eclock_AclmC)
#endif

  if (EClock(iac_EClock_Aglm) == 1) then
     iacymd_hold = EClock(iac_Eclock_ymd)
     EClock(iac_Eclock_ymd) = EClock(iac_Eclock_ymd) + 10000
#ifdef DEBUG
     write(iulog,*) trim(subname),'calling glm2iac_run',EClock(iac_EClock_ymd),EClock(iac_EClock_tod)
#endif

     call shr_cal_date2ymd(EClock(iac_EClock_ymd),yyyy,mm,dd)

     if (yyyy.gt.2005) then 
          write(iulog,*) trim(subname),'Data Mode IAC: Year ',yyyy,' is out of range'
          write(iulog,*) trim(subname),'Data Mode IAC is hardwired to work between the years 1850 to 2005'
          call shr_sys_abort(subname//' ERROR: iacnml error')
     end if
     ierr = nf90_open(trim(glm2iac_glmofile),nf90_nowrite,ncidglmo)
     call iac_ncerr(ierr,'open')

     allocate(array3(cdata%i(iac_cdatai_glm_nx),cdata%i(iac_cdatai_glm_ny),1))
     array3=0
     start3=1
     start3(3) = yyyy-1700+1
     count3(1) = cdata%i(iac_cdatai_glm_nx)
     count3(2) = cdata%i(iac_cdatai_glm_ny)
     count3(3) = 1
     ierr = nf90_inq_varid(ncidglmo,'gcrop',varid)
     call iac_ncerr(ierr,'inq_varid gcrop')
     ierr = nf90_get_var(ncidglmo,varid,array3,start3,count3)
     call iac_ncerr(ierr,'get_var gcrop')
     ij = 0
     do j = 1,cdata%i(iac_cdatai_glm_ny)
        do i = 1,cdata%i(iac_cdatai_glm_nx)
           ij = ij + 1
           glmo(iac_glmo_gcrop,ij) = array3(i,j,1)
        enddo
     enddo
     ierr = nf90_inq_varid(ncidglmo,'gpast',varid)
     call iac_ncerr(ierr,'inq_varid gpast')
     ierr = nf90_get_var(ncidglmo,varid,array3,start3,count3)
     call iac_ncerr(ierr,'get var gpast')
     ij = 0
     do j = 1,cdata%i(iac_cdatai_glm_ny)
        do i = 1,cdata%i(iac_cdatai_glm_nx)
           ij = ij + 1
           glmo(iac_glmo_gpast,ij) = array3(i,j,1)
        enddo
     enddo
     ierr = nf90_inq_varid(ncidglmo,'gothr',varid)
     call iac_ncerr(ierr,'inq var gothr')
     ierr = nf90_get_var(ncidglmo,varid,array3,start3,count3)
     call iac_ncerr(ierr,'get var gothr')
     ij = 0
     do j = 1,cdata%i(iac_cdatai_glm_ny)
        do i = 1,cdata%i(iac_cdatai_glm_nx)
           ij = ij + 1
           glmo(iac_glmo_gothr,ij) = array3(i,j,1)
        enddo
     enddo
     ierr = nf90_inq_varid(ncidglmo,'gsecd',varid)
     call iac_ncerr(ierr,'inq var gsecd')
     ierr = nf90_get_var(ncidglmo,varid,array3,start3,count3)
     call iac_ncerr(ierr,'get var gsecd')
     ij = 0
     do j = 1,cdata%i(iac_cdatai_glm_ny)
        do i = 1,cdata%i(iac_cdatai_glm_nx)
           ij = ij + 1
           glmo(iac_glmo_gsecd,ij) = array3(i,j,1)
        enddo
     enddo
! Harvest variables are rates and are offset by 1 year from the state
! ie read harvest rates for 2005 and states for 2006

     start3(3) = yyyy-1700
     
     ierr = nf90_inq_varid(ncidglmo,'gfvh1',varid)
     call iac_ncerr(ierr,'inq var gfvh1')
     ierr = nf90_get_var(ncidglmo,varid,array3,start3,count3)
     call iac_ncerr(ierr,'get var gfvh1')
     ij = 0
     do j = 1,cdata%i(iac_cdatai_glm_ny)
        do i = 1,cdata%i(iac_cdatai_glm_nx)
           ij = ij + 1
           glmo(iac_glmo_gfvh1,ij) = array3(i,j,1)
        enddo
     enddo
     ierr = nf90_inq_varid(ncidglmo,'gfvh2',varid)
     call iac_ncerr(ierr,'inq var gfvh2')
     ierr = nf90_get_var(ncidglmo,varid,array3,start3,count3)
     call iac_ncerr(ierr,'get var gfvh2')
     ij = 0
     do j = 1,cdata%i(iac_cdatai_glm_ny)
        do i = 1,cdata%i(iac_cdatai_glm_nx)
           ij = ij + 1
           glmo(iac_glmo_gfvh2,ij) = array3(i,j,1)
        enddo
     enddo
     ierr = nf90_inq_varid(ncidglmo,'gfsh1',varid)
     call iac_ncerr(ierr,'inq var gfsh1')
     ierr = nf90_get_var(ncidglmo,varid,array3,start3,count3)
     call iac_ncerr(ierr,'get var gfsh1')
     ij = 0
     do j = 1,cdata%i(iac_cdatai_glm_ny)
        do i = 1,cdata%i(iac_cdatai_glm_nx)
           ij = ij + 1
           glmo(iac_glmo_gfsh1,ij) = array3(i,j,1)
        enddo
     enddo
     ierr = nf90_inq_varid(ncidglmo,'gfsh2',varid)
     call iac_ncerr(ierr,'inq var gfsh2')
     ierr = nf90_get_var(ncidglmo,varid,array3,start3,count3)
     call iac_ncerr(ierr,'get var gfsh2')
     ij = 0
     do j = 1,cdata%i(iac_cdatai_glm_ny)
        do i = 1,cdata%i(iac_cdatai_glm_nx)
           ij = ij + 1
           glmo(iac_glmo_gfsh2,ij) = array3(i,j,1)
        enddo
     enddo
     ierr = nf90_inq_varid(ncidglmo,'gfsh3',varid)
     call iac_ncerr(ierr,'inq var gfsh3')
     ierr = nf90_get_var(ncidglmo,varid,array3,start3,count3)
     call iac_ncerr(ierr,'get var gfsh3')
     ij = 0
     do j = 1,cdata%i(iac_cdatai_glm_ny)
        do i = 1,cdata%i(iac_cdatai_glm_nx)
           ij = ij + 1
           glmo(iac_glmo_gfsh3,ij) = array3(i,j,1)
        enddo
     enddo
     ierr = nf90_close(ncidglmo)
     call iac_ncerr(ierr,'close')
     
     deallocate(array3)

     call iac_diag(' glmo: ',glmo)
     call glm2iac_run_mod(EClock, cdata, glmo, iaco)
     EClock(iac_Eclock_ymd) = iacymd_hold
!-------- history file ---------------------------------------------------
!jt     write(hfile,'(a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)') trim(casename)//'.iac.hi.',yyyy,'-',mm,'-',dd,'-',iactod,'.nc'
     write(hfile,'(a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)') 'DIAC1.iac.hi.',yyyy,'-',mm,'-',dd,'-',iactod,'.nc'
     write(iulog,*) trim(subname),' writing history file ',trim(hfile)
     
     nmode = ior(NF90_CLOBBER,NF90_64BIT_OFFSET)
     ierr = nf90_create(trim(hfile),nmode,ncid)
     call iac_ncerr(ierr,'create')
     ierr = nf90_put_att(ncid,NF90_GLOBAL,'missing_value',iac_spval)
     call iac_ncerr(ierr,'putatt_missval')
     
     ierr = nf90_def_dim(ncid,'glm_nx' ,cdata%i(iac_cdatai_glm_nx) ,dimidglm(1))
     call iac_ncerr(ierr,'defdim_glmnx')
     ierr = nf90_def_dim(ncid,'glm_ny' ,cdata%i(iac_cdatai_glm_ny) ,dimidglm(2))
     call iac_ncerr(ierr,'defdim_glmny')

     do n = 1,size(glmo,dim=1)
        write(vname,'(a,i2.2,a)') 'glmo',n,'_'//trim(iac_glmo_fld_names(n))
        ierr = nf90_def_var(ncid,vname,NF90_DOUBLE,dimidglm,varid)
        call iac_ncerr(ierr,'defvar_'//trim(vname))
     enddo

     ierr = nf90_enddef(ncid)
     call iac_ncerr(ierr,'enddef')

     allocate(array2(cdata%i(iac_cdatai_glm_nx),cdata%i(iac_cdatai_glm_ny)))
     do n = 1,size(glmo,dim=1)
        ij = 0
        do j = 1,cdata%i(iac_cdatai_glm_ny)
           do i = 1,cdata%i(iac_cdatai_glm_nx)
              ij = ij+1
              array2(i,j) = glmo(n,ij)
           enddo
        enddo
        write(vname,'(a,i2.2,a)') 'glmo',n,'_'//trim(iac_glmo_fld_names(n))
        ierr = nf90_inq_varid(ncid,vname,varid)
        call iac_ncerr(ierr,'inqvar_'//trim(vname))
        ierr = nf90_put_var(ncid,varid,array2)
        call iac_ncerr(ierr,'putvar_'//trim(vname))
     enddo
     deallocate(array2)
     
     ierr = nf90_close(ncid)
     call iac_ncerr(ierr,'close')
     
     !-------- end history file ------------
     
  endif ! Aglm = 1

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

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: iac_diag

! !INTERFACE:
  subroutine iac_diag(string,array)

! !DESCRIPTION:
! iac array diagnostic

!------------------------------------------------------------------------------

   implicit none
! !ARGUMENTS:
    character(len=*) :: string
    real*8, pointer :: array(:,:)

! !LOCAL VARIABLES:
    integer :: nj,j
    character(len=*),parameter :: subname='(iac_diag)'

! !REVISION HISTORY:
! Author: T Craig

!EOP
!---------------------------------------------------------------------------

     nj = size(array,dim=1)
     do j = 1,nj
        write(iulog,'(2a,i3,2f13.6)') trim(subname)//' ',trim(string),j, &
              minval(array(j,:)),maxval(array(j,:))
     enddo

  end subroutine iac_diag

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: iac_ncerr

! !INTERFACE:
  subroutine iac_ncerr(ierr,str)

! !DESCRIPTION:
! iac netcdf error diagnostic

!------------------------------------------------------------------------------

   implicit none
! !ARGUMENTS:
    integer :: ierr
    character(len=*),optional :: str

! !LOCAL VARIABLES:
    integer :: nj,j
    character(len=128) :: lstr
    character(len=*),parameter :: subname='(iac_ncerr)'

! !REVISION HISTORY:
! Author: T Craig

!EOP
!---------------------------------------------------------------------------

     lstr = ' '
     if (present(str)) then
        lstr = trim(str)
     endif

     if (ierr /= NF90_NOERR) then
        write(iulog,*) trim(subname),':',trim(lstr),':',trim(nf90_strerror(ierr))
     endif

  end subroutine iac_ncerr

!====================================================================================

end module iac_comp_mod
