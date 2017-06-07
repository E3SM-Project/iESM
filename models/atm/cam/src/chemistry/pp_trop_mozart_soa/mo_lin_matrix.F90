




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

         mat(760) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(863) = .890_r8*rxt(8) + .890_r8*rxt(9)
         mat(252) = rxt(41)

         mat(83) = -( rxt(43) + rxt(44) + rxt(45) + rxt(46)*y(11) + rxt(57)*y(4) &
                      + rxt(58)*y(4) + rxt(73)*y(15) + het_rates(3) )
         mat(735) = rxt(2)

         mat(248) = -( rxt(41) + het_rates(2) )
         mat(740) = rxt(3)
         mat(822) = rxt(5)
         mat(84) = rxt(43) + rxt(44)

         mat(729) = -( het_rates(5) )
         mat(830) = rxt(5) + .500_r8*rxt(208)
         mat(862) = .110_r8*rxt(8) + .110_r8*rxt(9)
         mat(88) = 2.000_r8*rxt(58)*y(4)

         mat(834) = -( rxt(5) + rxt(208) + het_rates(6) )
         mat(97) = rxt(6) + rxt(65)
         mat(228) = rxt(7)
         mat(866) = .890_r8*rxt(8) + .890_r8*rxt(9)
         mat(138) = .660_r8*rxt(10) + .660_r8*rxt(11) + rxt(71)
         mat(217) = .600_r8*rxt(19) + rxt(109)
         mat(260) = rxt(20) + rxt(176)
         mat(393) = rxt(30)

         mat(867) = -( rxt(8) + rxt(9) + rxt(207) + het_rates(7) )
         mat(98) = rxt(6) + rxt(65)
         mat(139) = .330_r8*rxt(10) + .330_r8*rxt(11)
         mat(218) = .400_r8*rxt(19)

         mat(226) = -( rxt(7) + het_rates(8) )
         mat(96) = 2.000_r8*rxt(206)
         mat(841) = rxt(207)
         mat(820) = .500_r8*rxt(208)

         mat(135) = -( rxt(10) + rxt(11) + rxt(71) + het_rates(9) )

         mat(95) = -( rxt(6) + rxt(65) + rxt(206) + het_rates(10) )

         mat(684) = -( rxt(47)*y(11) + rxt(72)*y(15) + rxt(81)*y(16) + rxt(82)*y(16) &
                      + rxt(217)*y(116) + rxt(218)*y(117) + het_rates(12) )
         mat(227) = rxt(7)
         mat(137) = .330_r8*rxt(10) + .330_r8*rxt(11)
         mat(148) = rxt(12)
         mat(65) = 2.000_r8*rxt(15)
         mat(211) = rxt(17)
         mat(202) = rxt(18)
         mat(457) = .330_r8*rxt(21) + .330_r8*rxt(22)
         mat(164) = rxt(24)
         mat(144) = rxt(25)
         mat(133) = rxt(26)
         mat(73) = rxt(29)
         mat(297) = rxt(37)
         mat(128) = rxt(38)
         mat(173) = rxt(39)
         mat(225) = rxt(40)
         mat(87) = 2.000_r8*rxt(45) + rxt(46)*y(11) + .750_r8*rxt(73)*y(15)
         mat(829) = .500_r8*rxt(208)

         mat(604) = -( rxt(216) + het_rates(13) )
         mat(136) = .660_r8*rxt(10) + .660_r8*rxt(11) + rxt(71)
         mat(147) = rxt(12)
         mat(463) = 2.000_r8*rxt(13)
         mat(412) = rxt(16)
         mat(210) = rxt(17)
         mat(456) = .670_r8*rxt(21) + .670_r8*rxt(22)
         mat(163) = rxt(24)
         mat(143) = rxt(25)
         mat(470) = rxt(28)
         mat(391) = rxt(30)
         mat(272) = rxt(31)
         mat(448) = rxt(32)
         mat(286) = 2.000_r8*rxt(33)
         mat(245) = .560_r8*rxt(35)
         mat(231) = 2.000_r8*rxt(36)
         mat(296) = .900_r8*rxt(37)
         mat(224) = rxt(40)
         mat(188) = rxt(86)
         mat(102) = rxt(93) + rxt(94)
         mat(86) = rxt(46)*y(11) + .400_r8*rxt(73)*y(15)
         mat(683) = rxt(47)*y(11) + rxt(81)*y(16) + rxt(82)*y(16) + rxt(217)*y(116) &
                      + rxt(218)*y(117)

         mat(63) = -( rxt(15) + het_rates(14) )
         mat(560) = .500_r8*rxt(216)

         mat(790) = -( het_rates(17) )
         mat(414) = rxt(16)
         mat(203) = rxt(18)
         mat(215) = .400_r8*rxt(19)
         mat(539) = .300_r8*rxt(23)
         mat(310) = rxt(27)
         mat(687) = rxt(72)*y(15)
         mat(89) = .750_r8*rxt(73)*y(15)

         mat(145) = -( rxt(12) + het_rates(18) )

         mat(462) = -( rxt(13) + rxt(14) + het_rates(19) )
         mat(146) = rxt(12)
         mat(209) = rxt(17)
         mat(452) = .670_r8*rxt(21) + .670_r8*rxt(22)
         mat(132) = rxt(26)
         mat(389) = rxt(30)
         mat(268) = .690_r8*rxt(31)
         mat(446) = rxt(32)
         mat(285) = rxt(33)
         mat(295) = .100_r8*rxt(37)
         mat(187) = rxt(86)
         mat(101) = 2.000_r8*rxt(94)
         mat(85) = .250_r8*rxt(73)*y(15)

         mat(262) = -( het_rates(20) )

         mat(109) = -( het_rates(21) )

         mat(150) = -( het_rates(22) )

         mat(99) = -( rxt(93) + rxt(94) + het_rates(23) )

         mat(179) = -( het_rates(24) )

         mat(197) = -( het_rates(25) )

         mat(284) = -( rxt(33) + het_rates(26) )
         mat(100) = rxt(93)

         mat(49) = -( het_rates(27) )

         mat(350) = -( het_rates(28) )
         mat(194) = rxt(34)

         mat(160) = -( rxt(24) + het_rates(29) )

         mat(411) = -( rxt(16) + het_rates(30) )
         mat(207) = rxt(17)
         mat(162) = rxt(24)
         mat(294) = .400_r8*rxt(37)
         mat(127) = rxt(38)

         mat(813) = -( het_rates(31) )
         mat(216) = .600_r8*rxt(19) + rxt(109)
         mat(459) = .670_r8*rxt(21) + .670_r8*rxt(22)
         mat(540) = .300_r8*rxt(23)
         mat(134) = rxt(26)
         mat(311) = rxt(27)
         mat(472) = rxt(28)
         mat(450) = rxt(32)
         mat(196) = rxt(34)
         mat(247) = .130_r8*rxt(35)
         mat(129) = rxt(38)

         mat(200) = -( rxt(18) + het_rates(32) )

         mat(398) = -( het_rates(33) )
         mat(527) = .700_r8*rxt(23)

         mat(46) = -( het_rates(34) )

         mat(360) = -( het_rates(35) )

         mat(140) = -( rxt(25) + het_rates(36) )

         mat(313) = -( het_rates(37) )

         mat(205) = -( rxt(17) + het_rates(38) )

         mat(307) = -( rxt(27) + het_rates(39) )
         mat(141) = .820_r8*rxt(25)
         mat(291) = .250_r8*rxt(37)
         mat(220) = .100_r8*rxt(40)

         mat(420) = -( het_rates(40) )

         mat(130) = -( rxt(26) + het_rates(41) )

         mat(34) = -( het_rates(42) )

         mat(118) = -( het_rates(43) )

         mat(37) = -( het_rates(47) )

         mat(335) = -( het_rates(48) )

         mat(289) = -( rxt(37) + het_rates(49) )

         mat(192) = -( rxt(34) + het_rates(44) )
         mat(288) = .800_r8*rxt(37)

         mat(300) = -( het_rates(45) )

         mat(125) = -( rxt(38) + het_rates(46) )

         mat(372) = -( het_rates(50) )

         mat(516) = -( het_rates(51) )

         mat(266) = -( rxt(31) + het_rates(52) )

         mat(533) = -( rxt(23) + het_rates(53) )
         mat(271) = .402_r8*rxt(31)
         mat(223) = rxt(40)

         mat(451) = -( rxt(21) + rxt(22) + het_rates(54) )
         mat(267) = .288_r8*rxt(31)
         mat(222) = rxt(40)

         mat(498) = -( het_rates(55) )

         mat(113) = -( het_rates(56) )

         mat(549) = -( het_rates(57) )
         mat(257) = rxt(20) + rxt(176)
         mat(455) = .330_r8*rxt(21) + .330_r8*rxt(22)

         mat(165) = -( het_rates(58) )

         mat(445) = -( rxt(32) + het_rates(59) )

         mat(469) = -( rxt(28) + het_rates(60) )
         mat(244) = .180_r8*rxt(35)
         mat(172) = .450_r8*rxt(39)

         mat(482) = -( het_rates(61) )

         mat(71) = -( rxt(29) + het_rates(62) )

         mat(274) = -( het_rates(63) )

         mat(433) = -( het_rates(64) )

         mat(219) = -( rxt(40) + het_rates(65) )

         mat(55) = -( het_rates(66) )

         mat(60) = -( het_rates(67) )

         mat(235) = -( het_rates(68) )

         mat(168) = -( rxt(39) + het_rates(69) )

         mat(74) = -( het_rates(78) )

         mat(243) = -( rxt(35) + het_rates(79) )
         mat(171) = .900_r8*rxt(39)

         mat(230) = -( rxt(36) + het_rates(80) )
         mat(242) = .130_r8*rxt(35)
         mat(169) = .450_r8*rxt(39)

         mat(40) = -( het_rates(70) )

         mat(80) = -( het_rates(71) )

         mat(1) = -( het_rates(72) )

         mat(2) = -( het_rates(73) )

         mat(43) = -( het_rates(74) )

         mat(92) = -( het_rates(75) )

         mat(3) = -( het_rates(76) )

         mat(4) = -( het_rates(77) )

         mat(212) = -( rxt(19) + rxt(109) + het_rates(81) )

         mat(174) = -( het_rates(82) )

         mat(254) = -( rxt(20) + rxt(176) + het_rates(83) )

         mat(321) = -( het_rates(84) )

         mat(388) = -( rxt(30) + het_rates(85) )

         mat(53) = -( het_rates(100) )

         mat(104) = -( het_rates(101) )

         mat(5) = -( het_rates(102) )

         mat(32) = -( het_rates(103) )

         mat(6) = -( het_rates(104) )


      end subroutine linmat01

      subroutine linmat02( mat, y, rxt, het_rates )
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

         mat(7) = -( het_rates(105) )

         mat(8) = -( het_rates(90) )

         mat(9) = -( het_rates(91) )

         mat(10) = -( het_rates(92) )

         mat(11) = -( het_rates(93) )

         mat(12) = -( het_rates(94) )

         mat(13) = -( het_rates(95) )

         mat(14) = -( het_rates(96) )

         mat(15) = -( het_rates(97) )

         mat(16) = -( het_rates(98) )

         mat(17) = -( het_rates(99) )

         mat(18) = -( rxt(209) + het_rates(86) )

         mat(20) = -( het_rates(87) )
         mat(19) = rxt(209)

         mat(21) = -( rxt(215) + het_rates(88) )

         mat(23) = -( het_rates(89) )
         mat(22) = rxt(215)

         mat(66) = -( het_rates(118) )

         mat(157) = -( het_rates(119) )

         mat(186) = -( rxt(86) + het_rates(120) )

         mat(24) = -( het_rates(106) )

         mat(25) = -( het_rates(107) )

         mat(26) = -( het_rates(108) )

         mat(27) = -( het_rates(109) )

         mat(28) = -( het_rates(110) )

         mat(29) = -( het_rates(111) )

         mat(30) = -( het_rates(112) )

         mat(31) = -( het_rates(113) )


      end subroutine linmat02

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
      call linmat02( mat, y, rxt, het_rates )

      end subroutine linmat

      end module mo_lin_matrix
