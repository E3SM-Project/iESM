




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

         mat(608) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(640) = .890_r8*rxt(8) + .890_r8*rxt(9)
         mat(223) = rxt(41)

         mat(61) = -( rxt(43) + rxt(44) + rxt(45) + rxt(46)*y(11) + rxt(57)*y(4) &
                      + rxt(58)*y(4) + rxt(73)*y(15) + het_rates(3) )
         mat(585) = rxt(2)

         mat(221) = -( rxt(41) + het_rates(2) )
         mat(590) = rxt(3)
         mat(750) = rxt(5)
         mat(62) = rxt(43) + rxt(44)

         mat(547) = -( het_rates(5) )
         mat(755) = rxt(5) + .500_r8*rxt(202)
         mat(638) = .110_r8*rxt(8) + .110_r8*rxt(9)
         mat(64) = 2.000_r8*rxt(58)*y(4)

         mat(761) = -( rxt(5) + rxt(202) + het_rates(6) )
         mat(71) = rxt(6) + rxt(65)
         mat(213) = rxt(7)
         mat(644) = .890_r8*rxt(8) + .890_r8*rxt(9)
         mat(126) = .660_r8*rxt(10) + .660_r8*rxt(11) + rxt(71)
         mat(209) = .600_r8*rxt(19) + rxt(109)
         mat(232) = rxt(20) + rxt(176)
         mat(366) = rxt(30)

         mat(641) = -( rxt(8) + rxt(9) + rxt(201) + het_rates(7) )
         mat(70) = rxt(6) + rxt(65)
         mat(124) = .330_r8*rxt(10) + .330_r8*rxt(11)
         mat(206) = .400_r8*rxt(19)

         mat(210) = -( rxt(7) + het_rates(8) )
         mat(69) = 2.000_r8*rxt(200)
         mat(621) = rxt(201)
         mat(749) = .500_r8*rxt(202)

         mat(123) = -( rxt(10) + rxt(11) + rxt(71) + het_rates(9) )

         mat(68) = -( rxt(6) + rxt(65) + rxt(200) + het_rates(10) )

         mat(740) = -( rxt(47)*y(11) + rxt(72)*y(15) + rxt(81)*y(16) + rxt(82)*y(16) &
                      + rxt(209)*y(86) + rxt(210)*y(87) + het_rates(12) )
         mat(212) = rxt(7)
         mat(125) = .330_r8*rxt(10) + .330_r8*rxt(11)
         mat(131) = rxt(12)
         mat(52) = 2.000_r8*rxt(15)
         mat(219) = rxt(17)
         mat(169) = rxt(18)
         mat(430) = .330_r8*rxt(21) + .330_r8*rxt(22)
         mat(96) = rxt(24)
         mat(111) = rxt(25)
         mat(137) = rxt(26)
         mat(56) = rxt(29)
         mat(270) = rxt(37)
         mat(102) = rxt(38)
         mat(146) = rxt(39)
         mat(183) = rxt(40)
         mat(66) = 2.000_r8*rxt(45) + rxt(46)*y(11) + .750_r8*rxt(73)*y(15)
         mat(760) = .500_r8*rxt(202)

         mat(812) = -( rxt(208) + het_rates(13) )
         mat(127) = .660_r8*rxt(10) + .660_r8*rxt(11) + rxt(71)
         mat(132) = rxt(12)
         mat(447) = 2.000_r8*rxt(13)
         mat(390) = rxt(16)
         mat(220) = rxt(17)
         mat(431) = .670_r8*rxt(21) + .670_r8*rxt(22)
         mat(97) = rxt(24)
         mat(112) = rxt(25)
         mat(439) = rxt(28)
         mat(367) = rxt(30)
         mat(242) = rxt(31)
         mat(423) = rxt(32)
         mat(260) = 2.000_r8*rxt(33)
         mat(202) = .560_r8*rxt(35)
         mat(187) = 2.000_r8*rxt(36)
         mat(271) = .900_r8*rxt(37)
         mat(184) = rxt(40)
         mat(157) = rxt(86)
         mat(75) = rxt(93) + rxt(94)
         mat(67) = rxt(46)*y(11) + .400_r8*rxt(73)*y(15)
         mat(742) = rxt(47)*y(11) + rxt(81)*y(16) + rxt(82)*y(16) + rxt(209)*y(86) &
                      + rxt(210)*y(87)

         mat(51) = -( rxt(15) + het_rates(14) )
         mat(764) = .500_r8*rxt(208)

         mat(577) = -( het_rates(17) )
         mat(385) = rxt(16)
         mat(167) = rxt(18)
         mat(205) = .400_r8*rxt(19)
         mat(508) = .300_r8*rxt(23)
         mat(291) = rxt(27)
         mat(736) = rxt(72)*y(15)
         mat(65) = .750_r8*rxt(73)*y(15)

         mat(128) = -( rxt(12) + het_rates(18) )

         mat(442) = -( rxt(13) + rxt(14) + het_rates(19) )
         mat(129) = rxt(12)
         mat(218) = rxt(17)
         mat(426) = .670_r8*rxt(21) + .670_r8*rxt(22)
         mat(135) = rxt(26)
         mat(362) = rxt(30)
         mat(237) = .690_r8*rxt(31)
         mat(420) = rxt(32)
         mat(258) = rxt(33)
         mat(268) = .100_r8*rxt(37)
         mat(153) = rxt(86)
         mat(74) = 2.000_r8*rxt(94)
         mat(63) = .250_r8*rxt(73)*y(15)

         mat(243) = -( het_rates(20) )

         mat(82) = -( het_rates(21) )

         mat(117) = -( het_rates(22) )

         mat(72) = -( rxt(93) + rxt(94) + het_rates(23) )

         mat(159) = -( het_rates(24) )

         mat(170) = -( het_rates(25) )

         mat(257) = -( rxt(33) + het_rates(26) )
         mat(73) = rxt(93)

         mat(23) = -( het_rates(27) )

         mat(323) = -( het_rates(28) )
         mat(175) = rxt(34)

         mat(93) = -( rxt(24) + het_rates(29) )

         mat(384) = -( rxt(16) + het_rates(30) )
         mat(216) = rxt(17)
         mat(95) = rxt(24)
         mat(267) = .400_r8*rxt(37)
         mat(100) = rxt(38)

         mat(664) = -( het_rates(31) )
         mat(207) = .600_r8*rxt(19) + rxt(109)
         mat(429) = .670_r8*rxt(21) + .670_r8*rxt(22)
         mat(511) = .300_r8*rxt(23)
         mat(136) = rxt(26)
         mat(292) = rxt(27)
         mat(436) = rxt(28)
         mat(421) = rxt(32)
         mat(176) = rxt(34)
         mat(200) = .130_r8*rxt(35)
         mat(101) = rxt(38)

         mat(165) = -( rxt(18) + het_rates(32) )

         mat(371) = -( het_rates(33) )
         mat(500) = .700_r8*rxt(23)

         mat(26) = -( het_rates(34) )

         mat(333) = -( het_rates(35) )

         mat(108) = -( rxt(25) + het_rates(36) )

         mat(281) = -( het_rates(37) )

         mat(214) = -( rxt(17) + het_rates(38) )

         mat(289) = -( rxt(27) + het_rates(39) )
         mat(109) = .820_r8*rxt(25)
         mat(264) = .250_r8*rxt(37)
         mat(179) = .100_r8*rxt(40)

         mat(393) = -( het_rates(40) )

         mat(133) = -( rxt(26) + het_rates(41) )

         mat(29) = -( het_rates(42) )

         mat(86) = -( het_rates(43) )

         mat(32) = -( het_rates(47) )

         mat(308) = -( het_rates(48) )

         mat(262) = -( rxt(37) + het_rates(49) )

         mat(173) = -( rxt(34) + het_rates(44) )
         mat(261) = .800_r8*rxt(37)

         mat(273) = -( het_rates(45) )

         mat(98) = -( rxt(38) + het_rates(46) )

         mat(345) = -( het_rates(50) )

         mat(489) = -( het_rates(51) )

         mat(235) = -( rxt(31) + het_rates(52) )

         mat(506) = -( rxt(23) + het_rates(53) )
         mat(240) = .402_r8*rxt(31)
         mat(182) = rxt(40)

         mat(424) = -( rxt(21) + rxt(22) + het_rates(54) )
         mat(236) = .288_r8*rxt(31)
         mat(181) = rxt(40)

         mat(471) = -( het_rates(55) )

         mat(103) = -( het_rates(56) )

         mat(828) = -( het_rates(57) )
         mat(234) = rxt(20) + rxt(176)
         mat(432) = .330_r8*rxt(21) + .330_r8*rxt(22)

         mat(138) = -( het_rates(58) )

         mat(418) = -( rxt(32) + het_rates(59) )

         mat(434) = -( rxt(28) + het_rates(60) )
         mat(199) = .180_r8*rxt(35)
         mat(145) = .450_r8*rxt(39)

         mat(455) = -( het_rates(61) )

         mat(54) = -( rxt(29) + het_rates(62) )

         mat(247) = -( het_rates(63) )

         mat(406) = -( het_rates(64) )

         mat(178) = -( rxt(40) + het_rates(65) )

         mat(38) = -( het_rates(66) )

         mat(43) = -( het_rates(67) )

         mat(190) = -( het_rates(68) )

         mat(141) = -( rxt(39) + het_rates(69) )

         mat(57) = -( het_rates(70) )

         mat(198) = -( rxt(35) + het_rates(71) )
         mat(144) = .900_r8*rxt(39)

         mat(185) = -( rxt(36) + het_rates(72) )
         mat(197) = .130_r8*rxt(35)
         mat(142) = .450_r8*rxt(39)

         mat(203) = -( rxt(19) + rxt(109) + het_rates(73) )

         mat(147) = -( het_rates(74) )

         mat(227) = -( rxt(20) + rxt(176) + het_rates(75) )

         mat(294) = -( het_rates(76) )

         mat(361) = -( rxt(30) + het_rates(77) )

         mat(36) = -( het_rates(79) )

         mat(77) = -( het_rates(80) )

         mat(21) = -( het_rates(81) )

         mat(1) = -( het_rates(82) )

         mat(2) = -( het_rates(83) )

         mat(3) = -( het_rates(78) )

         mat(46) = -( het_rates(88) )

         mat(113) = -( het_rates(89) )

         mat(152) = -( rxt(86) + het_rates(90) )

         mat(4) = -( het_rates(91) )

         mat(5) = -( het_rates(92) )

         mat(6) = -( het_rates(93) )

         mat(7) = -( het_rates(94) )


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

         mat(8) = -( het_rates(95) )

         mat(9) = -( het_rates(96) )

         mat(10) = -( het_rates(97) )

         mat(11) = -( het_rates(98) )

         mat(12) = -( het_rates(99) )

         mat(13) = -( het_rates(100) )

         mat(14) = -( het_rates(101) )

         mat(15) = -( het_rates(102) )

         mat(16) = -( het_rates(103) )

         mat(17) = -( het_rates(104) )

         mat(18) = -( het_rates(105) )

         mat(19) = -( het_rates(106) )

         mat(20) = -( het_rates(107) )


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
