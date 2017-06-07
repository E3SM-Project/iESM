!===========================================================================
! Combine several reactions into one pseudo reaction to correct the 
! photolysis rate J(O1D) to incorporate the effect of the other reactions. 
!
! Creator: Philip Cameron-Smith
!===========================================================================

module llnl_O1D_to_2OH_adj

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

  private
  public :: O1D_to_2OH_adj, O1D_to_2OH_adj_init

  integer :: jo1d_ndx

contains
!===========================================================================

!===========================================================================
!===========================================================================
  subroutine O1D_to_2OH_adj_init
    use mo_chem_utls, only : get_rxt_ndx
    implicit none

    jo1d_ndx  = get_rxt_ndx( 'jo1d' )

  end subroutine O1D_to_2OH_adj_init

!===========================================================================
!===========================================================================
  subroutine O1D_to_2OH_adj( p_rate, inv, m, ncol, tfld )

    use chem_mods,    only : nfs, phtcnt, rxntot, nfs !PJC added rxntot, nfs
    use ppgrid,       only : pcols, pver              !PJC added pcols
    use mo_setinv,    only : n2_ndx, o2_ndx, h2o_ndx  !PJC

    implicit none

    !--------------------------------------------------------------------
    ! ... dummy arguments
    !--------------------------------------------------------------------
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: inv(ncol,pver,nfs)
    real(r8), intent(in) :: m(ncol,pver)
    real(r8), intent(inout) :: p_rate(ncol,pver,rxntot)
    real(r8), intent(in)    :: tfld(pcols,pver)               ! midpoint temperature (K)

    !--------------------------------------------------------------------
    ! ... local variables
    !--------------------------------------------------------------------
    integer :: k
    real(r8) :: im(ncol)
    real(r8) :: n2_rate(ncol,pver)
    real(r8) :: o2_rate(ncol,pver)
    real(r8) :: h2o_rate(ncol,pver)

    real(r8), parameter :: x1 = 2.15e-11_r8
    real(r8), parameter :: x2 = 3.30e-11_r8
    real(r8), parameter :: x3 = 1.63e-10_r8
    real(r8), parameter :: y1 = 110.0_r8
    real(r8), parameter :: y2 =  55.0_r8
    real(r8), parameter :: y3 =  60.0_r8

    n2_rate(:,:)  = x1 * Exp ( y1 / tfld(:ncol,:)) * inv(:,:,n2_ndx)
    o2_rate(:,:)  = x2 * Exp ( y2 / tfld(:ncol,:)) * inv(:,:,o2_ndx)
    h2o_rate(:,:) = x3 * Exp ( y3 / tfld(:ncol,:)) * inv(:,:,h2o_ndx)

    p_rate(:,:,jo1d_ndx) = p_rate(:,:,jo1d_ndx) *   &
                          (h2o_rate(:,:) / (h2o_rate(:,:) + n2_rate(:,:) + o2_rate(:,:)))

  end subroutine O1D_to_2OH_adj

end module llnl_O1D_to_2OH_adj
