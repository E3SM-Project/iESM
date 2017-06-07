




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

         mat(1113) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(1077) = rxt(71)

         mat(1076) = -( rxt(71) + het_rates(2) )
         mat(1112) = rxt(3)
         mat(895) = rxt(5)
         mat(926) = rxt(6)
         mat(113) = rxt(8)
         mat(756) = rxt(10)
         mat(770) = rxt(46)
         mat(76) = rxt(49)
         mat(1054) = rxt(56)
         mat(341) = rxt(74) + rxt(75)
         mat(86) = rxt(102)

         mat(332) = -( rxt(74) + rxt(75) + rxt(77)*y(18) + rxt(78)*y(4) + rxt(79)*y(4) &
                      + rxt(80)*y(12) + rxt(81)*y(12) + rxt(82)*y(12) + rxt(260)*y(41) &
                      + rxt(261)*y(42) + rxt(262)*y(43) + rxt(263)*y(44) + rxt(264)*y(45) &
                      + rxt(265)*y(40) + rxt(266)*y(48) + rxt(267)*y(47) + het_rates(3) )
         mat(1085) = rxt(2)
         mat(758) = rxt(45)

         mat(83) = -( rxt(102) + het_rates(5) )
         mat(851) = rxt(5)

         mat(890) = -( rxt(5) + het_rates(6) )
         mat(921) = rxt(6) + .500_r8*rxt(253)
         mat(111) = rxt(8)
         mat(753) = rxt(11)
         mat(84) = rxt(102)
         mat(339) = 2.000_r8*rxt(79)*y(4)

         mat(922) = -( rxt(6) + rxt(253) + het_rates(7) )
         mat(112) = rxt(7) + rxt(112)
         mat(412) = rxt(9)
         mat(754) = rxt(10)
         mat(174) = rxt(13) + rxt(121)
         mat(212) = .600_r8*rxt(22) + rxt(164)
         mat(261) = rxt(23) + rxt(210)
         mat(433) = rxt(33)
         mat(472) = rxt(54)
         mat(284) = rxt(59)

         mat(1012) = -( rxt(94)*y(18) + rxt(122)*y(12) + rxt(129)*y(17) + rxt(130)*y(17) &
                      + rxt(311)*y(39) + rxt(312)*y(46) + rxt(313)*y(44) + rxt(314)*y(40) &
                 + het_rates(20) )
         mat(413) = rxt(9)
         mat(175) = rxt(12)
         mat(181) = rxt(14)
         mat(278) = rxt(20)
         mat(245) = rxt(21)
         mat(542) = .330_r8*rxt(24) + .330_r8*rxt(25)
         mat(139) = rxt(27)
         mat(187) = rxt(28)
         mat(155) = rxt(29)
         mat(134) = rxt(32)
         mat(322) = rxt(40)
         mat(125) = rxt(41)
         mat(161) = rxt(42)
         mat(206) = rxt(43)
         mat(768) = rxt(44)
         mat(300) = 2.000_r8*rxt(47)
         mat(439) = rxt(51)
         mat(390) = rxt(57)
         mat(923) = .500_r8*rxt(253)
         mat(340) = rxt(77)*y(18) + rxt(80)*y(12)

         mat(749) = -( rxt(10) + rxt(11) + rxt(252) + het_rates(8) )
         mat(110) = rxt(7) + rxt(8) + rxt(112)
         mat(171) = rxt(12)
         mat(210) = .400_r8*rxt(22)
         mat(470) = rxt(53)
         mat(283) = rxt(58)

         mat(409) = -( rxt(9) + het_rates(9) )
         mat(109) = 2.000_r8*rxt(251) + 2.000_r8*rxt(315) + 2.000_r8*rxt(321) &
                      + 2.000_r8*rxt(326)
         mat(733) = rxt(252)
         mat(908) = .500_r8*rxt(253)
         mat(466) = rxt(316) + rxt(322) + rxt(327)
         mat(281) = rxt(317) + rxt(325) + rxt(328)

         mat(170) = -( rxt(12) + rxt(13) + rxt(121) + het_rates(10) )

         mat(108) = -( rxt(7) + rxt(8) + rxt(112) + rxt(251) + rxt(315) + rxt(321) &
                      + rxt(326) + het_rates(11) )

         mat(714) = -( het_rates(13) )
         mat(480) = rxt(19)
         mat(242) = rxt(21)
         mat(209) = .400_r8*rxt(22)
         mat(604) = .300_r8*rxt(26)
         mat(372) = rxt(30)
         mat(335) = rxt(80)*y(12)
         mat(1005) = rxt(122)*y(12)
         mat(1121) = rxt(274)*y(12)

         mat(176) = -( rxt(14) + het_rates(14) )

         mat(87) = -( het_rates(113) )

         mat(44) = -( het_rates(114) )

         mat(842) = -( rxt(15) + rxt(16) + het_rates(16) )
         mat(180) = rxt(14)
         mat(277) = rxt(20)
         mat(541) = .670_r8*rxt(24) + .670_r8*rxt(25)
         mat(154) = rxt(29)
         mat(432) = rxt(33)
         mat(269) = .690_r8*rxt(34)
         mat(532) = rxt(35)
         mat(294) = rxt(36)
         mat(321) = .100_r8*rxt(40)
         mat(218) = rxt(137)
         mat(94) = 2.000_r8*rxt(160)
         mat(338) = rxt(81)*y(12) + rxt(82)*y(12)

         mat(443) = -( rxt(84) + het_rates(19) )
         mat(177) = rxt(14)
         mat(836) = 2.000_r8*rxt(15)
         mat(759) = rxt(44) + 2.000_r8*rxt(46)
         mat(582) = rxt(52)
         mat(333) = rxt(81)*y(12)
         mat(987) = rxt(94)*y(18) + rxt(130)*y(17)
         mat(1118) = rxt(269)*y(18)

         mat(823) = -( rxt(259) + het_rates(21) )
         mat(173) = rxt(13) + rxt(121)
         mat(483) = rxt(19)
         mat(276) = rxt(20)
         mat(540) = .670_r8*rxt(24) + .670_r8*rxt(25)
         mat(138) = rxt(27)
         mat(186) = rxt(28)
         mat(558) = rxt(31)
         mat(431) = rxt(33)
         mat(268) = rxt(34)
         mat(531) = rxt(35)
         mat(293) = 2.000_r8*rxt(36)
         mat(226) = .560_r8*rxt(38)
         mat(228) = 2.000_r8*rxt(39)
         mat(320) = .900_r8*rxt(40)
         mat(205) = rxt(43)
         mat(445) = rxt(84)
         mat(217) = rxt(137)
         mat(93) = rxt(159) + rxt(160)
         mat(337) = rxt(77)*y(18) + rxt(81)*y(12)
         mat(1008) = rxt(129)*y(17) + rxt(311)*y(39) + rxt(314)*y(40)
         mat(1124) = rxt(310)*y(39)

         mat(296) = -( rxt(47) + het_rates(22) )
         mat(793) = .500_r8*rxt(259)

         mat(763) = -( rxt(44) + rxt(45) + rxt(46) + het_rates(151) )
         mat(1007) = rxt(94)*y(18) + rxt(122)*y(12) + rxt(311)*y(39) + rxt(312)*y(46) &
                      + rxt(313)*y(44) + rxt(314)*y(40)

         mat(1133) = -( rxt(269)*y(18) + rxt(274)*y(12) + rxt(310)*y(39) + het_rates(25) &
       )
         mat(78) = 2.000_r8*rxt(48)
         mat(49) = 2.000_r8*rxt(50)
         mat(442) = rxt(51)
         mat(596) = rxt(52)
         mat(476) = rxt(53)
         mat(103) = rxt(55)
         mat(342) = 3.000_r8*rxt(260)*y(41) + 2.000_r8*rxt(261)*y(42) &
                      + 3.000_r8*rxt(262)*y(43) + rxt(263)*y(44) + 4.000_r8*rxt(264)*y(45)
         mat(1017) = rxt(311)*y(39) + 3.000_r8*rxt(312)*y(46) + rxt(313)*y(44)

         mat(77) = -( rxt(48) + het_rates(26) )

         mat(1033) = -( het_rates(27) )
         mat(75) = rxt(49)
         mat(474) = rxt(54)
         mat(48) = 2.000_r8*rxt(285)

         mat(74) = -( rxt(49) + het_rates(28) )

         mat(47) = -( rxt(50) + rxt(285) + het_rates(29) )

         mat(585) = -( rxt(52) + het_rates(30) )
         mat(1120) = rxt(269)*y(18) + rxt(274)*y(12) + 2.000_r8*rxt(310)*y(39)

         mat(436) = -( rxt(51) + het_rates(31) )
         mat(467) = rxt(316) + rxt(322) + rxt(327)

         mat(468) = -( rxt(53) + rxt(54) + rxt(316) + rxt(322) + rxt(327) + het_rates(32) &
       )

         mat(101) = -( rxt(55) + het_rates(33) )

         mat(545) = -( het_rates(34) )
         mat(102) = rxt(55)
         mat(1043) = rxt(56)
         mat(387) = rxt(57)
         mat(282) = rxt(58)
         mat(334) = rxt(265)*y(40) + rxt(266)*y(48) + rxt(267)*y(47)
         mat(996) = rxt(314)*y(40)

         mat(1053) = -( rxt(56) + het_rates(35) )
         mat(285) = rxt(59)

         mat(246) = -( het_rates(36) )

         mat(386) = -( rxt(57) + het_rates(37) )
         mat(280) = rxt(317) + rxt(325) + rxt(328)

         mat(279) = -( rxt(58) + rxt(59) + rxt(317) + rxt(325) + rxt(328) + het_rates(38) &
       )

         mat(251) = -( het_rates(15) )

         mat(104) = -( het_rates(49) )

         mat(140) = -( het_rates(50) )

         mat(91) = -( rxt(159) + rxt(160) + het_rates(51) )

         mat(189) = -( het_rates(52) )

         mat(287) = -( het_rates(53) )

         mat(292) = -( rxt(36) + het_rates(54) )
         mat(92) = rxt(159)

         mat(50) = -( het_rates(55) )

         mat(361) = -( het_rates(56) )
         mat(197) = rxt(37)

         mat(135) = -( rxt(27) + het_rates(57) )

         mat(478) = -( rxt(19) + het_rates(58) )
         mat(273) = rxt(20)
         mat(137) = rxt(27)
         mat(318) = .400_r8*rxt(40)
         mat(123) = rxt(41)

         mat(665) = -( het_rates(59) )
         mat(208) = .600_r8*rxt(22) + rxt(164)
         mat(537) = .670_r8*rxt(24) + .670_r8*rxt(25)
         mat(603) = .300_r8*rxt(26)
         mat(152) = rxt(29)
         mat(371) = rxt(30)
         mat(555) = rxt(31)
         mat(530) = rxt(35)
         mat(198) = rxt(37)
         mat(225) = .130_r8*rxt(38)
         mat(124) = rxt(41)


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

         mat(240) = -( rxt(21) + het_rates(60) )

         mat(452) = -( het_rates(61) )
         mat(597) = .700_r8*rxt(26)

         mat(54) = -( het_rates(62) )

         mat(416) = -( het_rates(63) )

         mat(182) = -( rxt(28) + het_rates(64) )

         mat(376) = -( het_rates(65) )

         mat(271) = -( rxt(20) + het_rates(66) )

         mat(369) = -( rxt(30) + het_rates(67) )
         mat(183) = .820_r8*rxt(28)
         mat(317) = .250_r8*rxt(40)
         mat(201) = .100_r8*rxt(43)

         mat(488) = -( het_rates(68) )

         mat(150) = -( rxt(29) + het_rates(69) )

         mat(35) = -( het_rates(70) )

         mat(114) = -( het_rates(71) )

         mat(38) = -( het_rates(75) )

         mat(347) = -( het_rates(76) )

         mat(313) = -( rxt(40) + het_rates(77) )

         mat(195) = -( rxt(37) + het_rates(72) )
         mat(312) = .800_r8*rxt(40)

         mat(324) = -( het_rates(73) )

         mat(121) = -( rxt(41) + het_rates(74) )

         mat(395) = -( het_rates(78) )

         mat(640) = -( het_rates(79) )

         mat(263) = -( rxt(34) + het_rates(80) )

         mat(601) = -( rxt(26) + het_rates(81) )
         mat(266) = .402_r8*rxt(34)
         mat(204) = rxt(43)

         mat(534) = -( rxt(24) + rxt(25) + het_rates(82) )
         mat(264) = .288_r8*rxt(34)
         mat(203) = rxt(43)

         mat(620) = -( het_rates(83) )

         mat(126) = -( het_rates(84) )

         mat(682) = -( het_rates(85) )
         mat(257) = rxt(23) + rxt(210)
         mat(538) = .330_r8*rxt(24) + .330_r8*rxt(25)

         mat(147) = -( het_rates(86) )

         mat(528) = -( rxt(35) + het_rates(87) )

         mat(554) = -( rxt(31) + het_rates(88) )
         mat(224) = .180_r8*rxt(38)
         mat(160) = .450_r8*rxt(42)

         mat(567) = -( het_rates(89) )

         mat(131) = -( rxt(32) + het_rates(90) )

         mat(303) = -( het_rates(91) )

         mat(515) = -( het_rates(92) )

         mat(200) = -( rxt(43) + het_rates(93) )

         mat(61) = -( het_rates(94) )

         mat(66) = -( het_rates(95) )

         mat(233) = -( het_rates(96) )

         mat(156) = -( rxt(42) + het_rates(97) )

         mat(79) = -( het_rates(98) )

         mat(222) = -( rxt(38) + het_rates(99) )
         mat(157) = .900_r8*rxt(42)

         mat(227) = -( rxt(39) + het_rates(100) )
         mat(223) = .130_r8*rxt(38)
         mat(158) = .450_r8*rxt(42)

         mat(207) = -( rxt(22) + rxt(164) + het_rates(101) )

         mat(162) = -( het_rates(102) )

         mat(255) = -( rxt(23) + rxt(210) + het_rates(103) )

         mat(501) = -( het_rates(104) )

         mat(428) = -( rxt(33) + het_rates(105) )

         mat(59) = -( het_rates(107) )

         mat(96) = -( het_rates(108) )

         mat(41) = -( het_rates(109) )

         mat(1) = -( het_rates(106) )

         mat(69) = -( het_rates(115) )

         mat(166) = -( het_rates(116) )

         mat(215) = -( rxt(137) + het_rates(117) )

         mat(2) = -( het_rates(118) )

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

         mat(20) = -( het_rates(136) )

         mat(21) = -( het_rates(137) )

         mat(22) = -( het_rates(138) )

         mat(23) = -( het_rates(139) )

         mat(24) = -( het_rates(140) )

         mat(25) = -( het_rates(141) )

         mat(26) = -( het_rates(142) )

         mat(27) = -( het_rates(143) )

         mat(28) = -( het_rates(144) )

         mat(29) = -( het_rates(145) )

         mat(30) = -( het_rates(146) )

         mat(31) = -( het_rates(147) )

         mat(32) = -( het_rates(148) )

         mat(33) = -( het_rates(149) )

         mat(34) = -( het_rates(150) )


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
