
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

      rate(:,:,77) = 8.00e-14_r8
      rate(:,:,78) = 3.90e-17_r8
      rate(:,:,81) = 4.20e-13_r8
      rate(:,:,82) = 8.50e-2_r8
      rate(:,:,83) = 1.30e-16_r8
      rate(:,:,85) = 1.00e-20_r8
      rate(:,:,86) = 2.58e-04_r8
      rate(:,:,93) = 1.20e-10_r8
      rate(:,:,94) = 1.70e-10_r8
      rate(:,:,95) = 1.20e-10_r8
      rate(:,:,96) = 1.50e-10_r8
      rate(:,:,97) = 7.20e-11_r8
      rate(:,:,98) = 2.84e-10_r8
      rate(:,:,99) = 1.80e-10_r8
      rate(:,:,100) = 9.60e-11_r8
      rate(:,:,101) = 4.10e-11_r8
      rate(:,:,102) = 1.125e-10_r8
      rate(:,:,103) = 3.00e-11_r8
      rate(:,:,104) = 7.50e-12_r8
      rate(:,:,105) = 1.10e-10_r8
      rate(:,:,106) = 1.50e-10_r8
      rate(:,:,107) = 1.50e-10_r8
      rate(:,:,108) = 5.00e-12_r8
      rate(:,:,109) = 7.00e-13_r8
      rate(:,:,124) = 1.00e-11_r8
      rate(:,:,125) = 2.20e-11_r8
      rate(:,:,126) = 3.50e-12_r8
      rate(:,:,134) = 5.80e-16_r8
      rate(:,:,141) = 7.20e-11_r8
      rate(:,:,142) = 6.90e-12_r8
      rate(:,:,143) = 1.60e-12_r8
      rate(:,:,147) = 1.80e-12_r8
      rate(:,:,150) = 1.80e-12_r8
      rate(:,:,175) = 1.70e-13_r8
      rate(:,:,222) = 1.e-10_r8
      rate(:,:,223) = 4.4e-10_r8
      rate(:,:,224) = 4.e-10_r8
      rate(:,:,225) = 2.e-10_r8
      rate(:,:,226) = 1.e-12_r8
      rate(:,:,227) = 6.e-11_r8
      rate(:,:,228) = 5.e-16_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,75) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,79) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,80) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,84) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,87) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,88) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:,89) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:,90) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,91) = 6.70e-11_r8 * exp_fac(:,:)
      rate(:,:,92) = 4.70e-11_r8 * exp_fac(:,:)
      rate(:,:,110) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:,111) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,112) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,165) = 2.70e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,114) = 3.50e-12_r8 * exp_fac(:,:)
      rate(:,:,146) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,185) = 1.70e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,115) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,194) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,116) = 5.20e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,118) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,168) = 3.00e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 170._r8 * itemp(:,:) )
      rate(:,:,123) = 1.50e-11_r8 * exp_fac(:,:)
      rate(:,:,158) = 1.80e-11_r8 * exp_fac(:,:)
      rate(:,:,128) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,130) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,131) = 2.80e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,132) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,133) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,151) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,193) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,135) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,136) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,200) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,140) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,144) = 2.20e-11_r8 * exp( 120._r8 * itemp(:,:) )
      rate(:,:,145) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,149) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,152) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,154) = 1.40e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      rate(:,:,155) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,156) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,157) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,159) = 4.10e-11_r8 * exp( -450._r8 * itemp(:,:) )
      rate(:,:,160) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,161) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,162) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      rate(:,:,163) = 7.40e-12_r8 * exp( 270._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,164) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,184) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,192) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,166) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,191) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,169) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,170) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,173) = 2.60e-12_r8 * exp( -350._r8 * itemp(:,:) )
      rate(:,:,174) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,176) = 2.50e-12_r8 * exp( -130._r8 * itemp(:,:) )
      rate(:,:,177) = 3.00e-12_r8 * exp( -500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,178) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,181) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,183) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,179) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,180) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,182) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      rate(:,:,186) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,187) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,190) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,189) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,195) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,196) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      rate(:,:,197) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,198) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,199) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      rate(:,:,201) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,113), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,117), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,119), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,121), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,127), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,137), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 4.7e-11_r8 * itemp(:,:)**0.2_r8
      call jpl( rate(1,1,139), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.6e-11_r8
      call jpl( rate(1,1,148), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,167), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 2.0e-12_r8 * itemp(:,:)**2.4_r8
      call jpl( rate(1,1,171), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,188), m, 0.6_r8, ko, kinf, n )

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

      rate(:,:kbot,77) = 8.00e-14_r8
      rate(:,:kbot,78) = 3.90e-17_r8
      rate(:,:kbot,82) = 8.50e-2_r8
      rate(:,:kbot,83) = 1.30e-16_r8
      rate(:,:kbot,85) = 1.00e-20_r8
      rate(:,:kbot,86) = 2.58e-04_r8
      rate(:,:kbot,108) = 5.00e-12_r8
      rate(:,:kbot,109) = 7.00e-13_r8
      rate(:,:kbot,142) = 6.90e-12_r8
      rate(:,:kbot,222) = 1.e-10_r8
      rate(:,:kbot,223) = 4.4e-10_r8
      rate(:,:kbot,224) = 4.e-10_r8
      rate(:,:kbot,225) = 2.e-10_r8
      rate(:,:kbot,226) = 1.e-12_r8
      rate(:,:kbot,227) = 6.e-11_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,75) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,79) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,80) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,84) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,87) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,88) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,89) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,110) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,111) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,114) = 3.50e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,146) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,115) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:kbot,116) = 5.20e-12_r8 * exp( 210._r8 * itemp(:,:) )
      rate(:,:kbot,140) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,144) = 2.20e-11_r8 * exp( 120._r8 * itemp(:,:) )
      rate(:,:kbot,145) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,151) = 3.00e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,152) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)







      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 4.7e-11_r8 * itemp(:,:)**0.2_r8
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,139) = wrk(:,:)





      end subroutine setrxt_hrates

      end module mo_setrxt
