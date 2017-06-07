





      module mo_lin_matrix

      private
      public :: linmat

      contains

      subroutine linmat01( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)

         mat(34) = -( rxt(1) + het_rates(1) )
         mat(28) = rxt(3)

         mat(66) = -( rxt(15) + rxt(16)*y(8) + het_rates(2) )
         mat(37) = 2.000_r8*rxt(1)
         mat(11) = 2.000_r8*rxt(2)
         mat(24) = rxt(6)

         mat(43) = -( het_rates(3) )
         mat(19) = 2.000_r8*rxt(4)
         mat(23) = rxt(6)
         mat(64) = rxt(16)*y(8)

         mat(9) = -( rxt(2) + het_rates(4) )

         mat(51) = -( het_rates(5) )
         mat(29) = rxt(3)

         mat(27) = -( rxt(3) + rxt(23) + het_rates(6) )

         mat(1) = -( het_rates(7) )
         mat(26) = .500_r8*rxt(23)

         mat(18) = -( rxt(4) + rxt(5) + het_rates(9) )
         mat(21) = rxt(6)

         mat(75) = -( het_rates(10) )
         mat(67) = rxt(15)

         mat(22) = -( rxt(6) + het_rates(11) )

         mat(3) = -( het_rates(15) )

         mat(7) = -( het_rates(14) )

         mat(2) = -( het_rates(13) )

         mat(12) = -( het_rates(12) )


      end subroutine linmat01

      subroutine linmat( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)

      call linmat01( mat, y, rxt, het_rates )

      end subroutine linmat

      end module mo_lin_matrix
