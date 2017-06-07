




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

         mat(863) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(776) = rxt(71)

         mat(773) = -( rxt(71) + het_rates(2) )
         mat(860) = rxt(3)
         mat(1108) = rxt(5)
         mat(823) = rxt(6)
         mat(96) = rxt(8)
         mat(733) = rxt(10)
         mat(875) = rxt(46)
         mat(61) = rxt(49)
         mat(894) = rxt(56)
         mat(321) = rxt(74) + rxt(75)
         mat(73) = rxt(102)

         mat(317) = -( rxt(74) + rxt(75) + rxt(77)*y(18) + rxt(78)*y(4) + rxt(79)*y(4) &
                      + rxt(80)*y(12) + rxt(81)*y(12) + rxt(82)*y(12) + rxt(260)*y(41) &
                      + rxt(261)*y(42) + rxt(262)*y(43) + rxt(263)*y(44) + rxt(264)*y(45) &
                      + rxt(265)*y(40) + rxt(266)*y(48) + rxt(267)*y(47) + het_rates(3) )
         mat(839) = rxt(2)
         mat(870) = rxt(45)

         mat(72) = -( rxt(102) + het_rates(5) )
         mat(1072) = rxt(5)

         mat(1117) = -( rxt(5) + het_rates(6) )
         mat(832) = rxt(6) + .500_r8*rxt(253)
         mat(98) = rxt(8)
         mat(740) = rxt(11)
         mat(75) = rxt(102)
         mat(327) = 2.000_r8*rxt(79)*y(4)

         mat(825) = -( rxt(6) + rxt(253) + het_rates(7) )
         mat(97) = rxt(7) + rxt(112)
         mat(409) = rxt(9)
         mat(734) = rxt(10)
         mat(163) = rxt(13) + rxt(121)
         mat(196) = .600_r8*rxt(22) + rxt(164)
         mat(256) = rxt(23) + rxt(210)
         mat(416) = rxt(33)
         mat(459) = rxt(54)
         mat(250) = rxt(59)

         mat(1008) = -( rxt(94)*y(18) + rxt(122)*y(12) + rxt(129)*y(17) + rxt(130)*y(17) &
                      + rxt(311)*y(39) + rxt(312)*y(46) + rxt(313)*y(44) + rxt(314)*y(40) &
                 + het_rates(20) )
         mat(411) = rxt(9)
         mat(165) = rxt(12)
         mat(172) = rxt(14)
         mat(270) = rxt(20)
         mat(212) = rxt(21)
         mat(524) = .330_r8*rxt(24) + .330_r8*rxt(25)
         mat(102) = rxt(27)
         mat(139) = rxt(28)
         mat(146) = rxt(29)
         mat(124) = rxt(32)
         mat(306) = rxt(40)
         mat(115) = rxt(41)
         mat(152) = rxt(42)
         mat(190) = rxt(43)
         mat(882) = rxt(44)
         mat(295) = 2.000_r8*rxt(47)
         mat(427) = rxt(51)
         mat(378) = rxt(57)
         mat(830) = .500_r8*rxt(253)
         mat(325) = rxt(77)*y(18) + rxt(80)*y(12)

         mat(732) = -( rxt(10) + rxt(11) + rxt(252) + het_rates(8) )
         mat(95) = rxt(7) + rxt(8) + rxt(112)
         mat(162) = rxt(12)
         mat(195) = .400_r8*rxt(22)
         mat(455) = rxt(53)
         mat(248) = rxt(58)

         mat(407) = -( rxt(9) + het_rates(9) )
         mat(94) = 2.000_r8*rxt(251) + 2.000_r8*rxt(315) + 2.000_r8*rxt(321) &
                      + 2.000_r8*rxt(326)
         mat(715) = rxt(252)
         mat(812) = .500_r8*rxt(253)
         mat(451) = rxt(316) + rxt(322) + rxt(327)
         mat(246) = rxt(317) + rxt(325) + rxt(328)

         mat(161) = -( rxt(12) + rxt(13) + rxt(121) + het_rates(10) )

         mat(93) = -( rxt(7) + rxt(8) + rxt(112) + rxt(251) + rxt(315) + rxt(321) &
                      + rxt(326) + het_rates(11) )

         mat(676) = -( het_rates(13) )
         mat(465) = rxt(19)
         mat(209) = rxt(21)
         mat(194) = .400_r8*rxt(22)
         mat(566) = .300_r8*rxt(26)
         mat(357) = rxt(30)
         mat(320) = rxt(80)*y(12)
         mat(997) = rxt(122)*y(12)
         mat(788) = rxt(274)*y(12)

         mat(167) = -( rxt(14) + het_rates(14) )

         mat(76) = -( het_rates(114) )

         mat(32) = -( het_rates(115) )

         mat(918) = -( rxt(15) + rxt(16) + het_rates(16) )
         mat(171) = rxt(14)
         mat(269) = rxt(20)
         mat(523) = .670_r8*rxt(24) + .670_r8*rxt(25)
         mat(145) = rxt(29)
         mat(417) = rxt(33)
         mat(241) = .690_r8*rxt(34)
         mat(513) = rxt(35)
         mat(274) = rxt(36)
         mat(305) = .100_r8*rxt(40)
         mat(203) = rxt(137)
         mat(66) = 2.000_r8*rxt(160)
         mat(324) = rxt(81)*y(12) + rxt(82)*y(12)

         mat(428) = -( rxt(84) + het_rates(19) )
         mat(168) = rxt(14)
         mat(907) = 2.000_r8*rxt(15)
         mat(871) = rxt(44) + 2.000_r8*rxt(46)
         mat(691) = rxt(52)
         mat(318) = rxt(81)*y(12)
         mat(980) = rxt(94)*y(18) + rxt(130)*y(17)
         mat(786) = rxt(269)*y(18)

         mat(1069) = -( rxt(259) + het_rates(21) )
         mat(166) = rxt(13) + rxt(121)
         mat(470) = rxt(19)
         mat(271) = rxt(20)
         mat(525) = .670_r8*rxt(24) + .670_r8*rxt(25)
         mat(103) = rxt(27)
         mat(140) = rxt(28)
         mat(542) = rxt(31)
         mat(419) = rxt(33)
         mat(243) = rxt(34)
         mat(515) = rxt(35)
         mat(276) = 2.000_r8*rxt(36)
         mat(217) = .560_r8*rxt(38)
         mat(220) = 2.000_r8*rxt(39)
         mat(307) = .900_r8*rxt(40)
         mat(191) = rxt(43)
         mat(433) = rxt(84)
         mat(205) = rxt(137)
         mat(67) = rxt(159) + rxt(160)
         mat(326) = rxt(77)*y(18) + rxt(81)*y(12)
         mat(1009) = rxt(129)*y(17) + rxt(311)*y(39) + rxt(314)*y(40)
         mat(800) = rxt(310)*y(39)

         mat(290) = -( rxt(47) + het_rates(22) )
         mat(1031) = .500_r8*rxt(259)

         mat(879) = -( rxt(44) + rxt(45) + rxt(46) + het_rates(136) )
         mat(1005) = rxt(94)*y(18) + rxt(122)*y(12) + rxt(311)*y(39) + rxt(312)*y(46) &
                      + rxt(313)*y(44) + rxt(314)*y(40)

         mat(793) = -( rxt(269)*y(18) + rxt(274)*y(12) + rxt(310)*y(39) + het_rates(25) )
         mat(63) = 2.000_r8*rxt(48)
         mat(22) = 2.000_r8*rxt(50)
         mat(425) = rxt(51)
         mat(698) = rxt(52)
         mat(458) = rxt(53)
         mat(82) = rxt(55)
         mat(322) = 3.000_r8*rxt(260)*y(41) + 2.000_r8*rxt(261)*y(42) &
                      + 3.000_r8*rxt(262)*y(43) + rxt(263)*y(44) + 4.000_r8*rxt(264)*y(45)
         mat(1002) = rxt(311)*y(39) + 3.000_r8*rxt(312)*y(46) + rxt(313)*y(44)

         mat(62) = -( rxt(48) + het_rates(26) )

         mat(750) = -( het_rates(27) )
         mat(60) = rxt(49)
         mat(456) = rxt(54)
         mat(21) = 2.000_r8*rxt(285)

         mat(59) = -( rxt(49) + het_rates(28) )

         mat(20) = -( rxt(50) + rxt(285) + het_rates(29) )

         mat(694) = -( rxt(52) + het_rates(30) )
         mat(789) = rxt(269)*y(18) + rxt(274)*y(12) + 2.000_r8*rxt(310)*y(39)

         mat(421) = -( rxt(51) + het_rates(31) )
         mat(452) = rxt(316) + rxt(322) + rxt(327)

         mat(453) = -( rxt(53) + rxt(54) + rxt(316) + rxt(322) + rxt(327) + het_rates(32) &
       )

         mat(80) = -( rxt(55) + het_rates(33) )

         mat(527) = -( het_rates(34) )
         mat(81) = rxt(55)
         mat(890) = rxt(56)
         mat(372) = rxt(57)
         mat(247) = rxt(58)
         mat(319) = rxt(265)*y(40) + rxt(266)*y(48) + rxt(267)*y(47)
         mat(989) = rxt(314)*y(40)

         mat(899) = -( rxt(56) + het_rates(35) )
         mat(251) = rxt(59)

         mat(231) = -( het_rates(36) )

         mat(371) = -( rxt(57) + het_rates(37) )
         mat(245) = rxt(317) + rxt(325) + rxt(328)

         mat(244) = -( rxt(58) + rxt(59) + rxt(317) + rxt(325) + rxt(328) + het_rates(38) &
       )

         mat(260) = -( het_rates(15) )

         mat(89) = -( het_rates(49) )

         mat(125) = -( het_rates(50) )

         mat(64) = -( rxt(159) + rxt(160) + het_rates(51) )

         mat(174) = -( het_rates(52) )

         mat(286) = -( het_rates(53) )

         mat(273) = -( rxt(36) + het_rates(54) )
         mat(65) = rxt(159)

         mat(35) = -( het_rates(55) )

         mat(346) = -( het_rates(56) )
         mat(182) = rxt(37)

         mat(99) = -( rxt(27) + het_rates(57) )

         mat(463) = -( rxt(19) + het_rates(58) )
         mat(266) = rxt(20)
         mat(101) = rxt(27)
         mat(303) = .400_r8*rxt(40)
         mat(113) = rxt(41)

         mat(644) = -( het_rates(59) )
         mat(193) = .600_r8*rxt(22) + rxt(164)
         mat(520) = .670_r8*rxt(24) + .670_r8*rxt(25)
         mat(565) = .300_r8*rxt(26)
         mat(143) = rxt(29)
         mat(356) = rxt(30)
         mat(537) = rxt(31)
         mat(512) = rxt(35)
         mat(183) = rxt(37)
         mat(216) = .130_r8*rxt(38)
         mat(114) = rxt(41)


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

         mat(207) = -( rxt(21) + het_rates(60) )

         mat(437) = -( het_rates(61) )
         mat(559) = .700_r8*rxt(26)

         mat(39) = -( het_rates(62) )

         mat(396) = -( het_rates(63) )

         mat(135) = -( rxt(28) + het_rates(64) )

         mat(361) = -( het_rates(65) )

         mat(264) = -( rxt(20) + het_rates(66) )

         mat(354) = -( rxt(30) + het_rates(67) )
         mat(136) = .820_r8*rxt(28)
         mat(302) = .250_r8*rxt(40)
         mat(186) = .100_r8*rxt(43)

         mat(499) = -( het_rates(68) )

         mat(141) = -( rxt(29) + het_rates(69) )

         mat(23) = -( het_rates(70) )

         mat(104) = -( het_rates(71) )

         mat(26) = -( het_rates(75) )

         mat(332) = -( het_rates(76) )

         mat(298) = -( rxt(40) + het_rates(77) )

         mat(180) = -( rxt(37) + het_rates(72) )
         mat(297) = .800_r8*rxt(40)

         mat(309) = -( het_rates(73) )

         mat(111) = -( rxt(41) + het_rates(74) )

         mat(380) = -( het_rates(78) )

         mat(602) = -( het_rates(79) )

         mat(236) = -( rxt(34) + het_rates(80) )

         mat(563) = -( rxt(26) + het_rates(81) )
         mat(239) = .402_r8*rxt(34)
         mat(189) = rxt(43)

         mat(516) = -( rxt(24) + rxt(25) + het_rates(82) )
         mat(237) = .288_r8*rxt(34)
         mat(188) = rxt(43)

         mat(582) = -( het_rates(83) )

         mat(116) = -( het_rates(84) )

         mat(619) = -( het_rates(85) )
         mat(254) = rxt(23) + rxt(210)
         mat(519) = .330_r8*rxt(24) + .330_r8*rxt(25)

         mat(132) = -( het_rates(86) )

         mat(510) = -( rxt(35) + het_rates(87) )

         mat(536) = -( rxt(31) + het_rates(88) )
         mat(215) = .180_r8*rxt(38)
         mat(151) = .450_r8*rxt(42)

         mat(549) = -( het_rates(89) )

         mat(121) = -( rxt(32) + het_rates(90) )

         mat(277) = -( het_rates(91) )

         mat(486) = -( het_rates(92) )

         mat(185) = -( rxt(43) + het_rates(93) )

         mat(46) = -( het_rates(94) )

         mat(51) = -( het_rates(95) )

         mat(224) = -( het_rates(96) )

         mat(147) = -( rxt(42) + het_rates(97) )

         mat(68) = -( het_rates(98) )

         mat(213) = -( rxt(38) + het_rates(99) )
         mat(148) = .900_r8*rxt(42)

         mat(218) = -( rxt(39) + het_rates(100) )
         mat(214) = .130_r8*rxt(38)
         mat(149) = .450_r8*rxt(42)

         mat(192) = -( rxt(22) + rxt(164) + het_rates(101) )

         mat(153) = -( het_rates(102) )

         mat(252) = -( rxt(23) + rxt(210) + het_rates(103) )

         mat(473) = -( het_rates(104) )

         mat(413) = -( rxt(33) + het_rates(105) )

         mat(44) = -( het_rates(107) )

         mat(84) = -( het_rates(108) )

         mat(29) = -( het_rates(109) )

         mat(1) = -( het_rates(110) )

         mat(2) = -( het_rates(106) )

         mat(54) = -( het_rates(116) )

         mat(157) = -( het_rates(117) )

         mat(200) = -( rxt(137) + het_rates(118) )

         mat(3) = -( het_rates(119) )

         mat(4) = -( het_rates(120) )

         mat(5) = -( het_rates(121) )

         mat(6) = -( het_rates(122) )

         mat(7) = -( het_rates(123) )

         mat(8) = -( het_rates(124) )

         mat(9) = -( het_rates(125) )

         mat(10) = -( het_rates(126) )

         mat(11) = -( het_rates(127) )

         mat(12) = -( het_rates(128) )

         mat(13) = -( het_rates(129) )

         mat(14) = -( het_rates(130) )

         mat(15) = -( het_rates(131) )

         mat(16) = -( het_rates(132) )

         mat(17) = -( het_rates(133) )

         mat(18) = -( het_rates(134) )

         mat(19) = -( het_rates(135) )


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
