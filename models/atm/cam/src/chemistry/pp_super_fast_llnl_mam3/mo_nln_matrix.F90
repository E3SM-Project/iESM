      module mo_nln_matrix
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: nlnmat
      contains
      subroutine nlnmat01( mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat(34) = -(rxt(7)*y(2) + rxt(8)*y(3) + rxt(12)*y(5) + rxt(28)*y(12))
         mat(63) = -rxt(7)*y(1)
         mat(42) = -rxt(8)*y(1)
         mat(49) = -rxt(12)*y(1)
         mat(14) = -rxt(28)*y(1)
         mat(66) = -(rxt(7)*y(1) + rxt(9)*y(3) + rxt(11)*y(4) + rxt(14)*y(6) + rxt(17) &
                      *y(9) + rxt(19)*y(11) + (rxt(24) + rxt(25)) * y(15) + rxt(26) &
                      *y(14) + rxt(27)*y(12))
         mat(37) = -rxt(7)*y(2)
         mat(45) = -rxt(9)*y(2)
         mat(11) = -rxt(11)*y(2)
         mat(30) = -rxt(14)*y(2)
         mat(20) = -rxt(17)*y(2)
         mat(24) = -rxt(19)*y(2)
         mat(5) = -(rxt(24) + rxt(25)) * y(2)
         mat(8) = -rxt(26)*y(2)
         mat(16) = -rxt(27)*y(2)
         mat(37) = mat(37) + rxt(8)*y(3)
         mat(45) = mat(45) + rxt(8)*y(1) + rxt(13)*y(5)
         mat(52) = rxt(13)*y(3)
         mat(43) = -(rxt(8)*y(1) + rxt(9)*y(2) + 4._r8*rxt(10)*y(3) + rxt(13)*y(5) &
                      + rxt(18)*y(10))
         mat(35) = -rxt(8)*y(3)
         mat(64) = -rxt(9)*y(3)
         mat(50) = -rxt(13)*y(3)
         mat(72) = -rxt(18)*y(3)
         mat(35) = mat(35) + rxt(7)*y(2) + .060_r8*rxt(28)*y(12)
         mat(64) = mat(64) + rxt(7)*y(1) + rxt(11)*y(4) + rxt(17)*y(9)
         mat(10) = rxt(11)*y(2)
         mat(50) = mat(50) + rxt(21)*y(10)
         mat(19) = rxt(17)*y(2)
         mat(72) = mat(72) + rxt(21)*y(5) + 1.600_r8*rxt(22)*y(10)
         mat(15) = .060_r8*rxt(28)*y(1)
         mat(9) = -(rxt(11)*y(2))
         mat(58) = -rxt(11)*y(4)
         mat(39) = 2.000_r8*rxt(10)*y(3)
         mat(51) = -(rxt(12)*y(1) + rxt(13)*y(3) + rxt(21)*y(10))
         mat(36) = -rxt(12)*y(5)
         mat(44) = -rxt(13)*y(5)
         mat(73) = -rxt(21)*y(5)
         mat(27) = -(rxt(14)*y(2))
         mat(62) = -rxt(14)*y(6)
         mat(33) = rxt(12)*y(5)
         mat(41) = rxt(13)*y(5)
         mat(48) = rxt(12)*y(1) + rxt(13)*y(3) + rxt(21)*y(10)
         mat(70) = rxt(21)*y(5)
         mat(54) = rxt(14)*y(6)
         mat(26) = rxt(14)*y(2)
         mat(18) = -(rxt(17)*y(2))
         mat(60) = -rxt(17)*y(9)
         mat(32) = .870_r8*rxt(28)*y(12)
         mat(60) = mat(60) + rxt(20)*y(11)
         mat(47) = rxt(21)*y(10)
         mat(68) = rxt(21)*y(5) + 4.000_r8*rxt(22)*y(10)
         mat(21) = rxt(20)*y(2)
         mat(13) = .870_r8*rxt(28)*y(1)
         mat(75) = -(rxt(18)*y(3) + rxt(21)*y(5) + 4._r8*rxt(22)*y(10))
         mat(46) = -rxt(18)*y(10)
         mat(53) = -rxt(21)*y(10)
         mat(38) = 1.860_r8*rxt(28)*y(12)
         mat(67) = rxt(19)*y(11)
         mat(25) = rxt(19)*y(2)
         mat(17) = 1.860_r8*rxt(28)*y(1)
         mat(22) = -((rxt(19) + rxt(20)) * y(2))
         mat(61) = -(rxt(19) + rxt(20)) * y(11)
         mat(40) = rxt(18)*y(10)
         mat(69) = rxt(18)*y(3)
         mat(3) = -((rxt(24) + rxt(25)) * y(2))
         mat(56) = -(rxt(24) + rxt(25)) * y(15)
         mat(7) = -(rxt(26)*y(2))
         mat(57) = -rxt(26)*y(14)
         mat(57) = mat(57) + (rxt(24)+.750_r8*rxt(25))*y(15)
         mat(4) = (rxt(24)+.750_r8*rxt(25))*y(2)
         mat(55) = rxt(26)*y(14)
         mat(6) = rxt(26)*y(2)
         mat(12) = -(rxt(27)*y(2) + rxt(28)*y(1))
         mat(59) = -rxt(27)*y(12)
         mat(31) = -rxt(28)*y(12)
      end subroutine nlnmat01
      subroutine nlnmat_finit( mat, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat( 1) = lmat( 1)
         mat( 2) = lmat( 2)
         mat( 3) = mat( 3) + lmat( 3)
         mat( 7) = mat( 7) + lmat( 7)
         mat( 9) = mat( 9) + lmat( 9)
         mat( 11) = mat( 11) + lmat( 11)
         mat( 12) = mat( 12) + lmat( 12)
         mat( 18) = mat( 18) + lmat( 18)
         mat( 19) = mat( 19) + lmat( 19)
         mat( 21) = mat( 21) + lmat( 21)
         mat( 22) = mat( 22) + lmat( 22)
         mat( 23) = lmat( 23)
         mat( 24) = mat( 24) + lmat( 24)
         mat( 26) = mat( 26) + lmat( 26)
         mat( 27) = mat( 27) + lmat( 27)
         mat( 28) = lmat( 28)
         mat( 29) = lmat( 29)
         mat( 34) = mat( 34) + lmat( 34)
         mat( 37) = mat( 37) + lmat( 37)
         mat( 43) = mat( 43) + lmat( 43)
         mat( 51) = mat( 51) + lmat( 51)
         mat( 64) = mat( 64) + lmat( 64)
         mat( 66) = mat( 66) + lmat( 66)
         mat( 67) = mat( 67) + lmat( 67)
         mat( 75) = mat( 75) + lmat( 75)
         mat( 65) = 0._r8
         mat( 71) = 0._r8
         mat( 74) = 0._r8
         mat( 1) = mat( 1) - dti
         mat( 2) = mat( 2) - dti
         mat( 3) = mat( 3) - dti
         mat( 7) = mat( 7) - dti
         mat( 9) = mat( 9) - dti
         mat( 12) = mat( 12) - dti
         mat( 18) = mat( 18) - dti
         mat( 22) = mat( 22) - dti
         mat( 27) = mat( 27) - dti
         mat( 34) = mat( 34) - dti
         mat( 43) = mat( 43) - dti
         mat( 51) = mat( 51) - dti
         mat( 66) = mat( 66) - dti
         mat( 75) = mat( 75) - dti
      end subroutine nlnmat_finit
      subroutine nlnmat( mat, y, rxt, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
      call nlnmat01( mat, y, rxt )
      call nlnmat_finit( mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
