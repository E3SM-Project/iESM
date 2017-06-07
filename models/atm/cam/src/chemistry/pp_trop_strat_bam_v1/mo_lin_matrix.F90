




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

         mat(788) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(917) = rxt(71)

         mat(920) = -( rxt(71) + het_rates(2) )
         mat(791) = rxt(3)
         mat(838) = rxt(5)
         mat(970) = rxt(6)
         mat(96) = rxt(8)
         mat(737) = rxt(10)
         mat(985) = rxt(46)
         mat(55) = rxt(49)
         mat(939) = rxt(56)
         mat(354) = rxt(74) + rxt(75)
         mat(73) = rxt(102)

         mat(348) = -( rxt(74) + rxt(75) + rxt(77)*y(18) + rxt(78)*y(4) + rxt(79)*y(4) &
                      + rxt(80)*y(12) + rxt(81)*y(12) + rxt(82)*y(12) + rxt(262)*y(41) &
                      + rxt(263)*y(42) + rxt(264)*y(43) + rxt(265)*y(44) + rxt(266)*y(45) &
                      + rxt(267)*y(40) + rxt(268)*y(48) + rxt(269)*y(47) + het_rates(3) )
         mat(768) = rxt(2)
         mat(977) = rxt(45)

         mat(71) = -( rxt(102) + het_rates(5) )
         mat(799) = rxt(5)

         mat(836) = -( rxt(5) + het_rates(6) )
         mat(968) = rxt(6) + .500_r8*rxt(253)
         mat(95) = rxt(8)
         mat(735) = rxt(11)
         mat(72) = rxt(102)
         mat(352) = 2.000_r8*rxt(79)*y(4)

         mat(972) = -( rxt(6) + rxt(253) + het_rates(7) )
         mat(97) = rxt(7) + rxt(112)
         mat(408) = rxt(9)
         mat(738) = rxt(10)
         mat(170) = rxt(13) + rxt(121)
         mat(195) = .600_r8*rxt(22) + rxt(164)
         mat(264) = rxt(23) + rxt(210)
         mat(416) = rxt(33)
         mat(441) = rxt(54)
         mat(258) = rxt(59)

         mat(1078) = -( rxt(94)*y(18) + rxt(122)*y(12) + rxt(129)*y(17) + rxt(130)*y(17) &
                      + rxt(313)*y(39) + rxt(314)*y(46) + rxt(315)*y(44) + rxt(316)*y(40) &
                 + het_rates(20) )
         mat(410) = rxt(9)
         mat(172) = rxt(12)
         mat(177) = rxt(14)
         mat(241) = rxt(20)
         mat(210) = rxt(21)
         mat(526) = .330_r8*rxt(24) + .330_r8*rxt(25)
         mat(102) = rxt(27)
         mat(139) = rxt(28)
         mat(144) = rxt(29)
         mat(123) = rxt(32)
         mat(305) = rxt(40)
         mat(114) = rxt(41)
         mat(158) = rxt(42)
         mat(190) = rxt(43)
         mat(989) = rxt(44)
         mat(294) = 2.000_r8*rxt(47)
         mat(425) = rxt(51)
         mat(376) = rxt(57)
         mat(974) = .500_r8*rxt(253)
         mat(356) = rxt(77)*y(18) + rxt(80)*y(12)

         mat(733) = -( rxt(10) + rxt(11) + rxt(252) + het_rates(8) )
         mat(94) = rxt(7) + rxt(8) + rxt(112)
         mat(168) = rxt(12)
         mat(194) = .400_r8*rxt(22)
         mat(438) = rxt(53)
         mat(255) = rxt(58)

         mat(406) = -( rxt(9) + het_rates(9) )
         mat(93) = 2.000_r8*rxt(251) + 2.000_r8*rxt(317) + 2.000_r8*rxt(323) &
                      + 2.000_r8*rxt(328)
         mat(717) = rxt(252)
         mat(956) = .500_r8*rxt(253)
         mat(434) = rxt(318) + rxt(324) + rxt(329)
         mat(253) = rxt(319) + rxt(327) + rxt(330)

         mat(167) = -( rxt(12) + rxt(13) + rxt(121) + het_rates(10) )

         mat(92) = -( rxt(7) + rxt(8) + rxt(112) + rxt(251) + rxt(317) + rxt(323) &
                      + rxt(328) + het_rates(11) )

         mat(678) = -( het_rates(13) )
         mat(464) = rxt(19)
         mat(208) = rxt(21)
         mat(193) = .400_r8*rxt(22)
         mat(568) = .300_r8*rxt(26)
         mat(344) = rxt(30)
         mat(351) = rxt(80)*y(12)
         mat(1067) = rxt(122)*y(12)
         mat(1086) = rxt(276)*y(12)

         mat(173) = -( rxt(14) + het_rates(14) )

         mat(75) = -( het_rates(129) )

         mat(31) = -( het_rates(130) )

         mat(1117) = -( rxt(15) + rxt(16) + het_rates(16) )
         mat(178) = rxt(14)
         mat(242) = rxt(20)
         mat(527) = .670_r8*rxt(24) + .670_r8*rxt(25)
         mat(145) = rxt(29)
         mat(418) = rxt(33)
         mat(250) = .690_r8*rxt(34)
         mat(517) = rxt(35)
         mat(275) = rxt(36)
         mat(306) = .100_r8*rxt(40)
         mat(205) = rxt(137)
         mat(66) = 2.000_r8*rxt(160)
         mat(358) = rxt(81)*y(12) + rxt(82)*y(12)

         mat(427) = -( rxt(84) + het_rates(19) )
         mat(174) = rxt(14)
         mat(1103) = 2.000_r8*rxt(15)
         mat(978) = rxt(44) + 2.000_r8*rxt(46)
         mat(693) = rxt(52)
         mat(349) = rxt(81)*y(12)
         mat(1050) = rxt(94)*y(18) + rxt(130)*y(17)
         mat(1084) = rxt(271)*y(18)

         mat(897) = -( rxt(261) + het_rates(21) )
         mat(169) = rxt(13) + rxt(121)
         mat(466) = rxt(19)
         mat(239) = rxt(20)
         mat(524) = .670_r8*rxt(24) + .670_r8*rxt(25)
         mat(101) = rxt(27)
         mat(137) = rxt(28)
         mat(541) = rxt(31)
         mat(415) = rxt(33)
         mat(248) = rxt(34)
         mat(515) = rxt(35)
         mat(273) = 2.000_r8*rxt(36)
         mat(216) = .560_r8*rxt(38)
         mat(218) = 2.000_r8*rxt(39)
         mat(304) = .900_r8*rxt(40)
         mat(189) = rxt(43)
         mat(429) = rxt(84)
         mat(201) = rxt(137)
         mat(65) = rxt(159) + rxt(160)
         mat(353) = rxt(77)*y(18) + rxt(81)*y(12)
         mat(1073) = rxt(129)*y(17) + rxt(313)*y(39) + rxt(316)*y(40)
         mat(1092) = rxt(312)*y(39)

         mat(289) = -( rxt(47) + het_rates(22) )
         mat(865) = .500_r8*rxt(261)

         mat(988) = -( rxt(44) + rxt(45) + rxt(46) + het_rates(49) )
         mat(1077) = rxt(94)*y(18) + rxt(122)*y(12) + rxt(313)*y(39) + rxt(314)*y(46) &
                      + rxt(315)*y(44) + rxt(316)*y(40)

         mat(1098) = -( rxt(271)*y(18) + rxt(276)*y(12) + rxt(312)*y(39) + het_rates(25) &
       )
         mat(62) = 2.000_r8*rxt(48)
         mat(21) = 2.000_r8*rxt(50)
         mat(426) = rxt(51)
         mat(706) = rxt(52)
         mat(444) = rxt(53)
         mat(81) = rxt(55)
         mat(357) = 3.000_r8*rxt(262)*y(41) + 2.000_r8*rxt(263)*y(42) &
                      + 3.000_r8*rxt(264)*y(43) + rxt(265)*y(44) + 4.000_r8*rxt(266)*y(45)
         mat(1079) = rxt(313)*y(39) + 3.000_r8*rxt(314)*y(46) + rxt(315)*y(44)

         mat(61) = -( rxt(48) + het_rates(26) )

         mat(751) = -( het_rates(27) )
         mat(54) = rxt(49)
         mat(439) = rxt(54)
         mat(20) = 2.000_r8*rxt(287)

         mat(53) = -( rxt(49) + het_rates(28) )

         mat(19) = -( rxt(50) + rxt(287) + het_rates(29) )

         mat(696) = -( rxt(52) + het_rates(30) )
         mat(1087) = rxt(271)*y(18) + rxt(276)*y(12) + 2.000_r8*rxt(312)*y(39)

         mat(420) = -( rxt(51) + het_rates(31) )
         mat(435) = rxt(318) + rxt(324) + rxt(329)

         mat(436) = -( rxt(53) + rxt(54) + rxt(318) + rxt(324) + rxt(329) + het_rates(32) &
       )

         mat(79) = -( rxt(55) + het_rates(33) )

         mat(529) = -( het_rates(34) )
         mat(80) = rxt(55)
         mat(932) = rxt(56)
         mat(371) = rxt(57)
         mat(254) = rxt(58)
         mat(350) = rxt(267)*y(40) + rxt(268)*y(48) + rxt(269)*y(47)
         mat(1059) = rxt(316)*y(40)

         mat(940) = -( rxt(56) + het_rates(35) )
         mat(257) = rxt(59)

         mat(230) = -( het_rates(36) )

         mat(370) = -( rxt(57) + het_rates(37) )
         mat(252) = rxt(319) + rxt(327) + rxt(330)

         mat(251) = -( rxt(58) + rxt(59) + rxt(319) + rxt(327) + rxt(330) + het_rates(38) &
       )

         mat(267) = -( het_rates(15) )

         mat(88) = -( het_rates(50) )

         mat(124) = -( het_rates(51) )

         mat(63) = -( rxt(159) + rxt(160) + het_rates(52) )

         mat(147) = -( het_rates(53) )

         mat(276) = -( het_rates(54) )

         mat(272) = -( rxt(36) + het_rates(55) )
         mat(64) = rxt(159)

         mat(34) = -( het_rates(56) )

         mat(333) = -( het_rates(57) )
         mat(181) = rxt(37)

         mat(98) = -( rxt(27) + het_rates(58) )

         mat(462) = -( rxt(19) + het_rates(59) )
         mat(237) = rxt(20)
         mat(100) = rxt(27)
         mat(302) = .400_r8*rxt(40)
         mat(112) = rxt(41)

         mat(646) = -( het_rates(60) )
         mat(192) = .600_r8*rxt(22) + rxt(164)
         mat(522) = .670_r8*rxt(24) + .670_r8*rxt(25)
         mat(567) = .300_r8*rxt(26)
         mat(142) = rxt(29)
         mat(343) = rxt(30)
         mat(539) = rxt(31)
         mat(514) = rxt(35)
         mat(182) = rxt(37)
         mat(215) = .130_r8*rxt(38)
         mat(113) = rxt(41)


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

         mat(206) = -( rxt(21) + het_rates(61) )

         mat(448) = -( het_rates(62) )
         mat(561) = .700_r8*rxt(26)

         mat(38) = -( het_rates(63) )

         mat(395) = -( het_rates(64) )

         mat(134) = -( rxt(28) + het_rates(65) )

         mat(360) = -( het_rates(66) )

         mat(235) = -( rxt(20) + het_rates(67) )

         mat(341) = -( rxt(30) + het_rates(68) )
         mat(135) = .820_r8*rxt(28)
         mat(301) = .250_r8*rxt(40)
         mat(185) = .100_r8*rxt(43)

         mat(472) = -( het_rates(69) )

         mat(140) = -( rxt(29) + het_rates(70) )

         mat(22) = -( het_rates(71) )

         mat(103) = -( het_rates(72) )

         mat(25) = -( het_rates(76) )

         mat(319) = -( het_rates(77) )

         mat(297) = -( rxt(40) + het_rates(78) )

         mat(179) = -( rxt(37) + het_rates(73) )
         mat(296) = .800_r8*rxt(40)

         mat(308) = -( het_rates(74) )

         mat(110) = -( rxt(41) + het_rates(75) )

         mat(379) = -( het_rates(79) )

         mat(604) = -( het_rates(80) )

         mat(243) = -( rxt(34) + het_rates(81) )

         mat(565) = -( rxt(26) + het_rates(82) )
         mat(246) = .402_r8*rxt(34)
         mat(188) = rxt(43)

         mat(518) = -( rxt(24) + rxt(25) + het_rates(83) )
         mat(244) = .288_r8*rxt(34)
         mat(187) = rxt(43)

         mat(584) = -( het_rates(84) )

         mat(115) = -( het_rates(85) )

         mat(621) = -( het_rates(86) )
         mat(261) = rxt(23) + rxt(210)
         mat(521) = .330_r8*rxt(24) + .330_r8*rxt(25)

         mat(131) = -( het_rates(87) )

         mat(512) = -( rxt(35) + het_rates(88) )

         mat(538) = -( rxt(31) + het_rates(89) )
         mat(214) = .180_r8*rxt(38)
         mat(157) = .450_r8*rxt(42)

         mat(551) = -( het_rates(90) )

         mat(120) = -( rxt(32) + het_rates(91) )

         mat(280) = -( het_rates(92) )

         mat(499) = -( het_rates(93) )

         mat(184) = -( rxt(43) + het_rates(94) )

         mat(45) = -( het_rates(95) )

         mat(50) = -( het_rates(96) )

         mat(223) = -( het_rates(97) )

         mat(153) = -( rxt(42) + het_rates(98) )

         mat(67) = -( het_rates(99) )

         mat(212) = -( rxt(38) + het_rates(100) )
         mat(154) = .900_r8*rxt(42)

         mat(217) = -( rxt(39) + het_rates(101) )
         mat(213) = .130_r8*rxt(38)
         mat(155) = .450_r8*rxt(42)

         mat(191) = -( rxt(22) + rxt(164) + het_rates(102) )

         mat(159) = -( het_rates(103) )

         mat(259) = -( rxt(23) + rxt(210) + het_rates(104) )

         mat(485) = -( het_rates(105) )

         mat(412) = -( rxt(33) + het_rates(106) )

         mat(43) = -( het_rates(112) )

         mat(83) = -( het_rates(113) )

         mat(1) = -( het_rates(114) )

         mat(28) = -( het_rates(115) )

         mat(2) = -( het_rates(116) )

         mat(3) = -( het_rates(117) )

         mat(4) = -( het_rates(111) )

         mat(5) = -( rxt(254) + het_rates(107) )

         mat(7) = -( het_rates(108) )
         mat(6) = rxt(254)

         mat(8) = -( rxt(260) + het_rates(109) )

         mat(10) = -( het_rates(110) )
         mat(9) = rxt(260)

         mat(56) = -( het_rates(131) )

         mat(163) = -( het_rates(132) )

         mat(199) = -( rxt(137) + het_rates(133) )

         mat(11) = -( het_rates(118) )

         mat(12) = -( het_rates(119) )

         mat(13) = -( het_rates(120) )

         mat(14) = -( het_rates(121) )

         mat(15) = -( het_rates(122) )

         mat(16) = -( het_rates(123) )

         mat(17) = -( het_rates(124) )

         mat(18) = -( het_rates(125) )


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
