
      module mo_sim_dat

      private
      public :: set_sim_dat

      contains

      subroutine set_sim_dat

      use chem_mods,   only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass
      use chem_mods,   only : diag_map
      use chem_mods,   only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use chem_mods,   only : pht_alias_lst, pht_alias_mult
      use chem_mods,   only : extfrc_lst, inv_lst, slvd_lst
      use abortutils,  only : endrun
      use mo_tracname, only : solsym
      use chem_mods,   only : frc_from_dataset
      use shr_kind_mod,only : r8 => shr_kind_r8
      use cam_logfile, only : iulog

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios

      clscnt(:) = (/      0,     0,     0,     5,     0 /)

      cls_rxt_cnt(:,4) = (/      0,     5,     0,     5 /)

      solsym(:  5) = (/ 'CH4             ','N2O             ','CFC11           ','CFC12           ','H2O             ' /)

      solsym(: 5) = (/ 'CH4     ','N2O     ','CFC11   ','CFC12   ','H2O     ' /)
      adv_mass(: 5) = (/ 16.0405998_r8, 44.0128784_r8, 137.367508_r8, 120.913208_r8, 18.0142002_r8 /)

      clsmap(:  5,4) = (/    1,   2,   3,   4,   5 /)

      permute(:  5,4) = (/    1,   2,   3,   4,   5 /)

      diag_map(:  5) = (/    1,   3,   4,   5,   6 /)

      inv_lst(:  3) = (/ 'M               ', 'N2              ', 'O2              ' /)

      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_map; error = ',ios
         call endrun
      end if
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'ch4_loss        ', 'n2o_loss        ', 'cfc11_loss      ', 'cfc12_loss      ', &
                                     'lyman_alpha     ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
