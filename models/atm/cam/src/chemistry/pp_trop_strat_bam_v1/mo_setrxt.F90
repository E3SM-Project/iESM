
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,76) = 2.20e-10_r8
      rate(:,:,77) = 1.10e-10_r8
      rate(:,:,78) = 4.90e-11_r8
      rate(:,:,79) = 6.70e-11_r8
      rate(:,:,80) = 1.125e-10_r8
      rate(:,:,81) = 3.00e-11_r8
      rate(:,:,82) = 7.50e-12_r8
      rate(:,:,86) = 7.20e-11_r8
      rate(:,:,87) = 6.90e-12_r8
      rate(:,:,88) = 1.60e-12_r8
      rate(:,:,92) = 1.80e-12_r8
      rate(:,:,95) = 1.80e-12_r8
      rate(:,:,116) = 1.00e-11_r8
      rate(:,:,117) = 2.20e-11_r8
      rate(:,:,118) = 3.50e-12_r8
      rate(:,:,135) = 4.5e-13_r8
      rate(:,:,144) = 7.00e-13_r8
      rate(:,:,147) = 2.00e-13_r8
      rate(:,:,148) = 6.80e-14_r8
      rate(:,:,157) = 1.00e-12_r8
      rate(:,:,159) = 1.00e-14_r8
      rate(:,:,161) = 1.00e-11_r8
      rate(:,:,162) = 1.10e-11_r8
      rate(:,:,165) = 4.00e-14_r8
      rate(:,:,182) = 3.00e-12_r8
      rate(:,:,185) = 6.80e-13_r8
      rate(:,:,186) = 5.40e-11_r8
      rate(:,:,198) = 2.40e-12_r8
      rate(:,:,201) = 1.40e-11_r8
      rate(:,:,204) = 5.00e-12_r8
      rate(:,:,216) = 2.40e-12_r8
      rate(:,:,220) = 1.40e-11_r8
      rate(:,:,222) = 2.40e-12_r8
      rate(:,:,224) = 3.50e-12_r8
      rate(:,:,225) = 4.50e-11_r8
      rate(:,:,232) = 2.40e-12_r8
      rate(:,:,242) = 3.00e-12_r8
      rate(:,:,243) = 1.00e-11_r8
      rate(:,:,250) = 2.10e-6_r8
      rate(:,:,254) = 7.10e-6_r8
      rate(:,:,260) = 7.10e-6_r8
      rate(:,:,262) = 1.70e-10_r8
      rate(:,:,263) = 1.20e-10_r8
      rate(:,:,264) = 1.50e-10_r8
      rate(:,:,265) = 7.20e-11_r8
      rate(:,:,266) = 2.84e-10_r8
      rate(:,:,267) = 1.80e-10_r8
      rate(:,:,268) = 9.60e-11_r8
      rate(:,:,269) = 4.10e-11_r8
      rate(:,:,290) = 1.70e-13_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,72) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,74) = 2.10e-11_r8 * exp( 115._r8 * itemp(:,:) )
      rate(:,:,75) = 3.20e-11_r8 * exp( 70._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 100._r8 * itemp(:,:) )
      rate(:,:,83) = 7.70e-11_r8 * exp_fac(:,:)
      rate(:,:,103) = 2.10e-11_r8 * exp_fac(:,:)
      rate(:,:,85) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,89) = 2.20e-11_r8 * exp( 120._r8 * itemp(:,:) )
      rate(:,:,90) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,91) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,106) = 3.50e-12_r8 * exp_fac(:,:)
      rate(:,:,300) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,94) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,97) = 7.80e-13_r8 * exp( -1050._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,98) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,125) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,134) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,149) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,172) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,176) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,181) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,193) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,202) = 2.30e-11_r8 * exp_fac(:,:)
      rate(:,:,218) = 1.52e-11_r8 * exp_fac(:,:)
      rate(:,:,230) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,241) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,249) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,308) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,99) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,101) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,213) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,102) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,104) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,280) = 2.70e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,107) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,309) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,108) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,110) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,283) = 3.00e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 170._r8 * itemp(:,:) )
      rate(:,:,115) = 1.50e-11_r8 * exp_fac(:,:)
      rate(:,:,273) = 1.80e-11_r8 * exp_fac(:,:)
      rate(:,:,120) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,122) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 300._r8 * itemp(:,:) )
      rate(:,:,123) = 2.80e-12_r8 * exp_fac(:,:)
      rate(:,:,178) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,124) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,126) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,127) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,128) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,315) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,131) = 5.00e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,132) = 1.90e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,133) = 2.90e-12_r8 * exp( -345._r8 * itemp(:,:) )
      rate(:,:,136) = 9.7e-15_r8 * exp( 625._r8 * itemp(:,:) )
      rate(:,:,137) = 2.4e12_r8 * exp( -7000._r8 * itemp(:,:) )
      rate(:,:,138) = 2.6e-12_r8 * exp( 265._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 700._r8 * itemp(:,:) )
      rate(:,:,139) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,146) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,170) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,175) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,179) = 8.60e-13_r8 * exp_fac(:,:)
      rate(:,:,192) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,199) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,217) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,223) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,229) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,233) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,240) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,248) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,141) = 8.70e-12_r8 * exp( -1070._r8 * itemp(:,:) )
      rate(:,:,143) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:,:) )
      rate(:,:,145) = 2.60e-12_r8 * exp( 365._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,150) = 5.60e-12_r8 * exp_fac(:,:)
      rate(:,:,152) = 8.10e-12_r8 * exp_fac(:,:)
      rate(:,:,278) = 7.40e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1900._r8 * itemp(:,:) )
      rate(:,:,151) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,167) = 6.50e-15_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._r8 * itemp(:,:) )
      rate(:,:,154) = 4.30e-13_r8 * exp_fac(:,:)
      rate(:,:,205) = 4.30e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,155) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,156) = 2.50e-12_r8 * exp_fac(:,:)
      rate(:,:,180) = 7.10e-13_r8 * exp_fac(:,:)
      rate(:,:,206) = 2.00e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 180._r8 * itemp(:,:) )
      rate(:,:,158) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,169) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,174) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,187) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,191) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,215) = 4.40e-12_r8 * exp_fac(:,:)
      rate(:,:,228) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,239) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,247) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,160) = 1.60e11_r8 * exp( -4150._r8 * itemp(:,:) )
      rate(:,:,163) = 6.90e-12_r8 * exp( -230._r8 * itemp(:,:) )
      rate(:,:,168) = 4.60e-13_r8 * exp( -1156._r8 * itemp(:,:) )
      rate(:,:,171) = 3.75e-13_r8 * exp( -40._r8 * itemp(:,:) )
      rate(:,:,173) = 1.00e-11_r8 * exp( -665._r8 * itemp(:,:) )
      rate(:,:,183) = 8.40e-13_r8 * exp( 830._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1860._r8 * itemp(:,:) )
      rate(:,:,184) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,226) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,188) = 4.13e-12_r8 * exp( 452._r8 * itemp(:,:) )
      rate(:,:,189) = 7.52e-16_r8 * exp( -1521._r8 * itemp(:,:) )
      rate(:,:,190) = 2.30e-12_r8 * exp( -170._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 175._r8 * itemp(:,:) )
      rate(:,:,194) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,227) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,195) = 4.40e-15_r8 * exp( -2500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,196) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,197) = 1.30e-13_r8 * exp_fac(:,:)
      rate(:,:,203) = 5.30e-12_r8 * exp_fac(:,:)
      rate(:,:,221) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,231) = 2.7e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 400._r8 * itemp(:,:) )
      rate(:,:,200) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,219) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,234) = 5.00e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 530._r8 * itemp(:,:) )
      rate(:,:,207) = 4.60e-12_r8 * exp_fac(:,:)
      rate(:,:,208) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,212) = 2.54e-11_r8 * exp( 410._r8 * itemp(:,:) )
      rate(:,:,214) = 3.03e-12_r8 * exp( -446._r8 * itemp(:,:) )
      rate(:,:,235) = 1.30e-12_r8 * exp( 640._r8 * itemp(:,:) )
      rate(:,:,236) = 1.90e-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,238) = 1.70e-12_r8 * exp( 352._r8 * itemp(:,:) )
      rate(:,:,244) = 1.2e-11_r8 * exp( 444._r8 * itemp(:,:) )
      rate(:,:,245) = 1.e-15_r8 * exp( -732._r8 * itemp(:,:) )
      rate(:,:,246) = 1.2e-12_r8 * exp( 490._r8 * itemp(:,:) )
      rate(:,:,256) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,258) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,259) = 1.70e-12_r8 * exp( -710._r8 * itemp(:,:) )
      rate(:,:,270) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,271) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,272) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,274) = 4.10e-11_r8 * exp( -450._r8 * itemp(:,:) )
      rate(:,:,275) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,276) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,277) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,279) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,299) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,307) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,281) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,306) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,284) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,285) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,288) = 2.60e-12_r8 * exp( -350._r8 * itemp(:,:) )
      rate(:,:,289) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,291) = 2.50e-12_r8 * exp( -130._r8 * itemp(:,:) )
      rate(:,:,292) = 3.00e-12_r8 * exp( -500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,293) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,296) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,298) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,294) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,295) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,297) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      rate(:,:,301) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,302) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,305) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,304) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,310) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,311) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      rate(:,:,312) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,313) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,314) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      rate(:,:,316) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 4.70e-11_r8 * itemp(:,:)**0.2_r8
      call jpl( rate(1,1,84), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,93), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.28e-33_r8
      kinf(:,:) = 9.30e-15_r8 * itemp(:,:)**(-4.42_r8)
      call jpl( rate(1,1,96), m, 0.8_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,105), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,109), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,111), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,113), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,119), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,129), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.5e-30_r8
      kinf(:,:) = 8.3e-13_r8 * itemp(:,:)**(-2.0_r8)
      call jpl( rate(1,1,140), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.00e-28_r8 * itemp(:,:)**0.8_r8
      kinf(:,:) = 8.80e-12_r8
      call jpl( rate(1,1,142), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.50e-29_r8 * itemp(:,:)**6.5_r8
      kinf(:,:) = 1.10e-11_r8 * itemp(:,:)
      call jpl( rate(1,1,153), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,166), m, 0.5_r8, ko, kinf, n )

      ko(:,:) = 8.e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.e-11_r8
      call jpl( rate(1,1,211), m, .5_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,282), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 2.0e-12_r8 * itemp(:,:)**2.4_r8
      call jpl( rate(1,1,286), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,303), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)


      end subroutine setrxt_hrates

      end module mo_setrxt
