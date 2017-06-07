module microp_driver

!-------------------------------------------------------------------------------------------------------
!
! Driver for CAM microphysics parameterizations
!
!-------------------------------------------------------------------------------------------------------

use shr_kind_mod,  only: r8 => shr_kind_r8
use ppgrid,        only: pcols, pver
use physics_types, only: physics_state, physics_ptend, physics_tend, &
                         physics_ptend_copy, physics_ptend_sum
use physics_buffer,only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, physics_buffer_desc, pbuf_add_field, &
                         dtype_r8, physics_buffer_desc
use phys_control,  only: phys_getopts

use cldwat2m_macro,only: ini_macro
use microp_aero,   only: microp_aero_init, microp_aero_run
use micro_mg_cam,  only: micro_mg_cam_readnl, micro_mg_cam_register, &
                         micro_mg_cam_implements_cnst, micro_mg_cam_init_cnst, &
                         micro_mg_cam_init, micro_mg_cam_tend
use cam_logfile,   only: iulog
use abortutils,    only: endrun
use perf_mod,      only: t_startf, t_stopf

implicit none
private
save

public :: &
   microp_driver_readnl,          &
   microp_driver_register,        &
   microp_driver_init_cnst,       &
   microp_driver_implements_cnst, &
   microp_driver_init,            &
   microp_driver_tend

character(len=16)  :: microp_scheme   ! Microphysics scheme

! physics buffer indices
integer ::  ast_idx
integer ::  prec_str_idx, snow_str_idx, prec_sed_idx, snow_sed_idx
integer ::  prec_pcw_idx, snow_pcw_idx

!===============================================================================
contains
!===============================================================================

subroutine microp_driver_readnl(nlfile)

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Read in namelist for microphysics scheme
   !-----------------------------------------------------------------------

   call phys_getopts(microp_scheme_out=microp_scheme)

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_readnl(nlfile)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_readnl:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_readnl

subroutine microp_driver_register

   ! Register microphysics constituents and fields in the physics buffer.
   !-----------------------------------------------------------------------


   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_register()
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_register:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_register

!===============================================================================

function microp_driver_implements_cnst(name)

   ! Return true if specified constituent is implemented by the
   ! microphysics package

   character(len=*), intent(in) :: name        ! constituent name
   logical :: microp_driver_implements_cnst    ! return value

   ! Local workspace
   integer :: m
   !-----------------------------------------------------------------------

   microp_driver_implements_cnst = .false.

   select case (microp_scheme)
   case ('MG')
      microp_driver_implements_cnst = micro_mg_cam_implements_cnst(name)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_implements_cnst:: unrecognized microp_scheme')
   end select

end function microp_driver_implements_cnst

!===============================================================================

subroutine microp_driver_init_cnst(name, q, gcid)

   ! Initialize the microphysics constituents, if they are
   ! not read from the initial file.

   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,          intent(in)  :: gcid(:)  ! global column id
   !-----------------------------------------------------------------------

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_init_cnst(name, q, gcid)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_init_cnst:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_init_cnst

!===============================================================================

subroutine microp_driver_init(pbuf2d)

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! Initialize the microphysics parameterizations
   !-----------------------------------------------------------------------

   call ini_macro()
   call microp_aero_init()

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_init(pbuf2d)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_init:: unrecognized microp_scheme')
   end select


   ! Get indices for pbuf fields
   ast_idx      = pbuf_get_index('AST')
   prec_str_idx = pbuf_get_index('PREC_STR')
   snow_str_idx = pbuf_get_index('SNOW_STR')
   prec_sed_idx = pbuf_get_index('PREC_SED')
   snow_sed_idx = pbuf_get_index('SNOW_SED')
   prec_pcw_idx = pbuf_get_index('PREC_PCW')
   snow_pcw_idx = pbuf_get_index('SNOW_PCW')


end subroutine microp_driver_init

!===============================================================================

subroutine microp_driver_tend( &
             state, ptend, dtime, pbuf, cmeliq)

   ! Call the microphysics parameterization run methods.

   ! Input arguments

   type(physics_state), intent(in)    :: state       ! State variables
   type(physics_ptend), intent(out)   :: ptend       ! Package tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in)  :: dtime                    ! Timestep
   real(r8), intent(in)  :: cmeliq(pcols,pver)       ! Rate of cond-evap of liq within the cloud

   ! Local variables

   type(physics_ptend) :: ptend_aero     ! tendencies from aerosol microphysics
   type(physics_ptend) :: ptend_micro    ! tendencies from cloud microphysics

   integer :: lchnk
   integer :: ncol

   real(r8), pointer :: prec_str(:)          ! [Total] Sfc flux of precip from stratiform [ m/s ] 
   real(r8), pointer :: snow_str(:)          ! [Total] Sfc flux of snow from stratiform   [ m/s ]
   real(r8), pointer :: prec_sed(:)          ! Surface flux of total cloud water from sedimentation
   real(r8), pointer :: snow_sed(:)          ! Surface flux of cloud ice from sedimentation
   real(r8), pointer :: prec_pcw(:)          ! Sfc flux of precip from microphysics [ m/s ]
   real(r8), pointer :: snow_pcw(:)          ! Sfc flux of snow from microphysics [ m/s ]

   ! Output from microp_aero_run for aerosol actication
   real(r8) :: naai(pcols,pver)      !ice nucleation number
   real(r8) :: naai_hom(pcols,pver)  !ice nucleation number (homogeneous)
   real(r8) :: npccn(pcols,pver)     !liquid activation number tendency
   real(r8) :: rndst(pcols,pver,4)
   real(r8) :: nacon(pcols,pver,4)

   ! Physics buffer fields
   integer itim
   real(r8), pointer :: ast(:,:)       ! Relative humidity cloud fraction
   real(r8) :: alst_mic(pcols,pver)
   real(r8) :: aist_mic(pcols,pver)
   !======================================================================

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Microphysics assumes 'liquid stratus frac = ice stratus frac 
   !                      = max( liquid stratus frac, ice stratus frac )'.
   itim = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx,      ast, start=(/1,1,itim/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, prec_str_idx, prec_str)
   call pbuf_get_field(pbuf, snow_str_idx, snow_str)
   call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
   call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
   call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
   call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)

   alst_mic(:ncol,:pver) = ast(:ncol,:pver)
   aist_mic(:ncol,:pver) = ast(:ncol,:pver)

   ! calculate aerosol activation (naai for ice, npccn for liquid) and 
   ! dust size (rndst) and number (nacon) for contact nucleation

   call t_startf('microp_aero_run')
   call microp_aero_run(state, ptend_aero, dtime, pbuf, alst_mic, &
                        aist_mic, naai, naai_hom, npccn, rndst,   &
                        nacon)
   call t_stopf('microp_aero_run')

   ! Call MG Microphysics

   select case (microp_scheme)
   case ('MG')
      call t_startf('microp_mg_tend')
      call micro_mg_cam_tend(state, ptend_micro, dtime, pbuf, cmeliq, &
                             naai, naai_hom, npccn, rndst, nacon, &
                             prec_str, snow_str, prec_sed, snow_sed, &
                             prec_pcw, snow_pcw)
      call t_stopf('microp_mg_tend')
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_tend:: unrecognized microp_scheme')
   end select

   ! combine aero and micro tendencies
   call physics_ptend_copy(ptend_micro, ptend)
   call physics_ptend_sum(ptend_aero, ptend, state)

end subroutine microp_driver_tend

end module microp_driver
