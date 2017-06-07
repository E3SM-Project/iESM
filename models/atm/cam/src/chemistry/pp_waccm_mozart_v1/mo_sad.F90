

      module mo_sad

      use shr_kind_mod,  only : r8 => shr_kind_r8
      use physconst,     only : pi
      use ppgrid,        only : pcols, pver
      use m_sad_data,    only : a, b
      use cam_logfile,   only : iulog
      

      implicit none

      private
      public  :: sad_inti
      public  :: sad_strat_calc
      public  :: sad_top
      public  :: sad_setopts
      public  :: sad_defaultopts

      save

      real(r8), parameter :: one_thrd = 1._r8/3._r8
      real(r8), parameter :: two_thrd = 2._r8/3._r8
      real(r8), parameter :: four_pi  = 4._r8*pi

      integer :: sad_top
      integer :: sad_topp

      logical :: rad_feedback = .false.
      integer :: h2so4_ndx
      integer :: h2so4_r_g_ndx  !index to geometric mean radius of wet sulfuric aerosol

    contains
        
      subroutine sad_defaultopts( strat_aero_feedback_out )
        implicit none
        logical, intent(out), optional :: strat_aero_feedback_out

        if ( present(strat_aero_feedback_out) ) then
           strat_aero_feedback_out = rad_feedback
        endif

      end subroutine sad_defaultopts

      subroutine sad_setopts( strat_aero_feedback_in )

        use physics_buffer, only : pbuf_add_field, dtype_r8

        implicit none
        logical, intent(in), optional :: strat_aero_feedback_in

        if ( present(strat_aero_feedback_in) ) then
           rad_feedback = strat_aero_feedback_in
        endif
        if ( rad_feedback ) then
           call pbuf_add_field('H2SO4M',         'global',dtype_r8,(/pcols,pver/), h2so4_ndx) 
           call pbuf_add_field('VOLC_RAD_GEOM' , 'global',dtype_r8,(/pcols,pver/), h2so4_r_g_ndx) 
        endif
      end subroutine sad_setopts

      subroutine sad_inti(pbuf2d)
!----------------------------------------------------------------------
!     ... initialize the sad module
!----------------------------------------------------------------------

      use time_manager, only : is_first_step
      use phys_grid,    only : get_ncols_p, get_lon_all_p, get_lat_all_p
      use cam_history,  only : addfld, phys_decomp
      use ref_pres,     only : pref_mid
      use physics_buffer, only : physics_buffer_desc, pbuf_set_field

      implicit none

      type(physics_buffer_desc), pointer :: pbuf2d(:,:)

!---------------------------------------------------------------------- 
!	... Local variables
!---------------------------------------------------------------------- 
      integer  ::  k


!---------------------------------------------------------------------- 
!	... find level where etamids are all > 1 hPa
!---------------------------------------------------------------------- 
      sad_top = 0
      do k = pver,1,-1
	 if( pref_mid(k) < .001_r8 ) then
             sad_top = k
             exit
         end if
      end do
      sad_topp = sad_top + 1
      write(iulog,*) 'sad_inti: sad capped at level ',sad_top
      write(iulog,*) '          whose midpoint is ',(pref_mid(sad_topp))*1.e3_r8,' hPa'

      if (rad_feedback .and. is_first_step()) then
         call pbuf_set_field(pbuf2d, h2so4_ndx,  0.0_r8)
         call pbuf_set_field(pbuf2d, h2so4_r_g_ndx,  1.0_r8)
      endif

      call addfld( 'H2SO4M_C',  'ug/m3 ',   pver, 'I', 'chemical sulfate aerosol mass',  phys_decomp )

      if ( rad_feedback ) then
         call addfld( 'SAD_SULFR', 'cm2/cm3 ', pver, 'I', 'radiative sulfate aerosol SAD',  phys_decomp )
         call addfld( 'RAD_SULFR', 'cm ',      pver, 'I', 'radiative sad sulfate',          phys_decomp )
         call addfld( 'H2SO4MMR',  'kg/kg',    pver, 'I', 'radiative sulfate aerosol mmr',  phys_decomp )
         call addfld( 'H2SO4M_R',  'ug/m3 ',   pver, 'I', 'radiative sulfate aerosol mass', phys_decomp )
         call addfld( 'VOLC_RAD_GEOM', 'm',    pver, 'I', 'geometric mean radius of wet aerosol', phys_decomp)
      endif

      end subroutine sad_inti
!===============================================================================
! ROUTINE
!   sad_strat_calc
!
!   Date...
!     14 October 1999
!
!   Programmed by...
!     Douglas E. Kinnison
!   Modified by
!     Stacy Walters
!     1 November 1999
!   Modified by 
!     Doug Kinnison
!     1 September 2004; Condensed phase H2O passed in from CAM
!     2 November 2004; New treatment of denitrificatoin (NAT)
!    14 November 2004; STS mode of operation.
!    27 March    2008; Using original NAT approach.
!
! DESCRIPTION
!
!     This routine has the logic to derive the surface area density for
!     three types of aerosols: Sulfate; Nitric Acid Trihydrate (NAT);
!     and ICE. The surface area density is stored in sad_strat(3). The
!     first, second, and third dimensions are SULFATE, NAT, and ICE SAD
!     respectively. The effective radius of each particle is also stored
!     in radius_strat(3).
!
!     NOTE1: For Sulfate and H2O ICE
!     The Surface Area and Volume Densities are calculated from the
!     second and third moment of the LOG Normal Distribution. For an example
!     see Binkowski et al., JGR, 100, 26191-26209, 1995. The Volume Density
!     is substituted into the SAD equation so that the SAD is dependent on
!     the # of particles cm-3, the width of the distribution (sigma), and
!     the volume density of the aerosol. For the ternary solution calculation
!     the total sulfate mass is derived from the SAGEII SAD data. This approach
!     has been previously used in Considine et al., JGR, 1999. The thermodynamic
!     models used in this routine are from A. Tabazedeh.
!
!     NOTE2: For NAT, there are two modes centered at: 1) 0.5um radius; 2) 6.5um
!     radius. This is based on ER2 measurements during NASA SOLVE, see Fahey et al., 
!     Science, 291, 1026-1031, 2001. . The number density of the large mode is set to 2.3e-4
!     Particle cm-3 and the condensed phase HNO3 mass assigned to this mode to 
!     produce particles with a mode radius of 6.5 microns. Any additional
!     HNO3 is then assigned to the smaller mode. Denitrification occurs on the
!     larger mode. The HNO3 in the smaller mode is assigned to STS aerosols.
!
! ARGUMENTS
!   INPUT:
!     hno3_gas		Nitric Acid gas phase abundance (mole fraction)
!     hno3_cond(3)	Nitric Acid condensed phase abundance (mole fraction)
!                       (1) in STS; (2) in 0.5um NAT; (3) 6.5um NAT
!     h2o_cond		Water condensed phase           (mole fraction)
!     h2o_gas           Water gas-phase abundance       (mole fraction)

!     sage_sad		SAGEII surface area density     (cm2 cm-3)
!     m                 Airdensity                      (molecules cm-3)
!     press             Pressure                        (hPa, i.e., mb)

!
!   OUTPUT:
!
!	hno3_gas     = Gas-phase HNO3... Used in chemical solver. 
!       hno3_cond(1) = Condensed HNO3 from STS... Not used in mo_aero_settling
!       hno3_cond(2) = Condensed HNO3 from NAT, sm mode... Not used in mo_aero_settling
!       hno3_cond(3) = Condensed HNO3 from NAT, lg mode... USED in mo_aero_settling
!
!       SAD_strat(1) = Sulfate Aerosol... Used in heterogeneous rate module 
!       SAD_strat(2) = NAT small mode.... Used in heterogeneous rate module
!       SAD_strat(3) = NAT large mode.... Not used in heterogeneous rate module.
!       SAD_strat(4) = Water-Ice......... Used in heterogeneous rate module
!
!       RAD_strat(1) = Sulfate Aerosol... Used in heterogeneous rate module
!       RAD_strat(2) = NAT small mode.... Not used in mo_aero_settling (0.5um)
!       RAD_strat(3) = NAT large mode.... Used in mo_aero_settling (6.5um)
!       RAD_strat(4) = Water-Ice......... Not used in heterogeneous rate module
!
!   NOTE1: The sum of hno3_cond(1-3) will be added to hno3_gas for advection of HNO3 in
!         WACCM3.
!
!   NOTE2: This routine does not partition H2O.
!
!
! ROUTINES CallED (in and below this routine):
!
!	sad_strat_calc
!         nat_sat_temp ...... derives the NAT saturation temp
!
!         sulfate_sad_calc .. Calculates the sulfate SAD;  HNO3, H2O cond. phase
!         calc_radius_lbs ... Calculates the radius for a H2SO4 Binary soln. (T>200K)
!         sad_to_h2so4 ...... Derives H2SO4 abundance (micrograms m-3)
!                             from SAGEII SAD.
!         density............ A. Tabazedeh binary thermo model
!         equil.............. A. Tabazedeh ternary thermo. model
!
!         nat_sad_calc....... Calculates the NAT SAD; HNO3, H2O cond. phase
!         nat_cond........... Derives the NAT HNO3 gas/cond partitioning
!
!         ice_sad_calc....... derives the ICE SAD and H2O gas/cond partitioning
!===============================================================================

      subroutine sad_strat_calc( lchnk, m, press, temper, hno3_gas, &
                                 hno3_cond_arg, h2o_gas, h2o_cond, sad_sage, radius_strat_arg, &
                                 sad_strat_arg, ncol, pbuf )

      use cam_history, only : outfld
      use physics_buffer, only : physics_buffer_desc

      implicit none

!-------------------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------------------
      integer, intent(in)     ::  lchnk                      ! chnk id
      integer, intent(in)     ::  ncol                       ! columns in chunk
      real(r8), intent(in)    ::  m(ncol,pver)               ! Air density (molecules cm-3)
      real(r8), intent(in)    ::  sad_sage(ncol,pver)        ! SAGEII surface area density (cm2 aer. cm-3 air)
      real(r8), intent(in)    ::  press(ncol,pver)           ! Pressure, hPa
      real(r8), intent(in)    ::  temper(pcols,pver)         ! Temperature (K)
      real(r8), intent(inout) ::  h2o_gas(ncol,pver)         ! H2O gas-phase (mole fraction)
      real(r8), intent(inout) ::  h2o_cond(ncol,pver)        ! H2O condensed phase (mole fraction)
      real(r8), intent(inout) ::  hno3_gas(ncol,pver)        ! HNO3 condensed phase (mole fraction)
      real(r8), intent(inout) ::  hno3_cond_arg(ncol,pver,2)     ! HNO3 condensed phase (mole fraction)
      real(r8), intent(out)   ::  radius_strat_arg(ncol,pver,3)  ! Radius of Sulfate, NAT_sm, NAT_lg, and ICE (cm)
      real(r8), intent(out)   ::  sad_strat_arg(ncol,pver,3)     ! Surface area density of Sulfate, NAT_sm, NAT_lg, ICE (cm2 cm-3)

      real(r8) ::  hno3_cond(ncol,pver,3)     ! HNO3 condensed phase (mole fraction)
      real(r8) ::  radius_strat(ncol,pver,4)  ! Radius of Sulfate, NAT_sm, NAT_lg, and ICE (cm)
      real(r8) ::  sad_strat(ncol,pver,4)     ! Surface area density of Sulfate, NAT_sm, NAT_lg, ICE (cm2 cm-3)

!-------------------------------------------------------------------------------
!	... local variables
!-------------------------------------------------------------------------------
      real(r8), parameter :: temp_floor = 0._r8

      integer  ::  i, k, n
      integer  ::  dims(1)
      real(r8) ::  hno3_total(ncol,pver)         ! HNO3 total = gas-phase + condensed
      real(r8) ::  h2o_total(ncol,pver)          ! H2O total  = gas-phase + condensed
      real(r8) ::  radius_lbs(ncol,pver)         ! Radius of Liquid Binary Sulfate (cm)
      real(r8) ::  radius_sulfate(ncol,pver)     ! Radius of Sulfate aerosol (cm)
      real(r8) ::  radius_nat_sm(ncol,pver)      ! Radius of NAT aerosol at 0.5um  (cm)
      real(r8) ::  radius_nat_lg(ncol,pver)      ! Radius of NAT aerosol at 6.5um  (cm)
      real(r8) ::  radius_ice(ncol,pver)         ! Radius of ICE aerosol     (cm)
      real(r8) ::  sad_nat_sm(ncol,pver)         ! SAD of NAT aerosol at 0.5um (cm2 cm-3)
      real(r8) ::  sad_nat_lg(ncol,pver)         ! SAD of NAT aerosol at 6.5um (cm2 cm-3)
      real(r8) ::  sad_sulfate(ncol,pver)        ! SAD of Sulfate aerosol    (cm2 cm-3)
      real(r8) ::  sad_ice(ncol,pver)            ! SAD of ICE aerosol        (cm2 cm-3)
      real(r8) ::  tsat_nat(ncol,pver)           ! Temperature for NAT saturation
      real(r8) ::  h2o_avail(ncol,pver)          ! H2O temporary arrays
      real(r8) ::  hno3_avail(ncol,pver)         ! HNO3 temporary array
      real(r8) ::  hno3_gas_nat(ncol,pver)       ! HNO3 after call to NAT routines
      real(r8) ::  hno3_gas_sulf(ncol,pver)      ! HNO3 after call to STS routines
      real(r8) ::  hno3_cond_nat_sm(ncol,pver)   ! HNO3 condensed after call to NAT, small mode
      real(r8) ::  hno3_cond_nat_lg(ncol,pver)   ! HNO3 condensed after call to NAT, large mode
      real(r8) ::  hno3_cond_sulf(ncol,pver)     ! HNO3 condensed after call to STS routines
      real(r8) ::  temp(pcols,pver)              ! wrk temperature array
      real(r8) ::  h2so4m(ncol,pver)             ! wrk array

      logical  ::  z_val(ncol)
      logical  ::  mask(ncol,pver)
      logical  ::  mask1(ncol,pver)
      logical  ::  mask2(ncol,pver)

      type(physics_buffer_desc), pointer :: pbuf(:)

      hno3_cond(:,:,:) = 0._r8
      radius_strat(:,:,:) = 0._r8
      sad_strat(:,:,:) = 0._r8

      if ( rad_feedback ) then
         call sad_strat_calc_rad( lchnk, ncol, m, press, temper, hno3_gas, &
                                  hno3_cond, h2o_gas, h2o_cond, sad_sage,  &
                                  pbuf )
      endif

!----------------------------------------------------------------------
!     ... initialize to zero
!----------------------------------------------------------------------
      do n = 1,4
         do k = 1,pver
            radius_strat(:,k,n) = 0._r8
            sad_strat(:,k,n)    = 0._r8
         end do
      end do

!----------------------------------------------------------------------
!     ... limit temperature
!----------------------------------------------------------------------
      do k = 1,pver
         h2so4m(:ncol,k)  = 0._r8
         temp(:ncol,k)    = max( temp_floor,temper(:ncol,k) )
      end do

!----------------------------------------------------------------------
!     ... total HNO3 and H2O gas and condensed phases
!----------------------------------------------------------------------
      do k = sad_topp,pver
         hno3_total(:,k) = hno3_gas(:,k) + hno3_cond(:,k,1) + hno3_cond(:,k,2) + hno3_cond(:,k,3)
	 h2o_total(:,k)  = h2o_gas(:,k)  + h2o_cond(:,k)
      end do

!----------------------------------------------------------------------
!     ... Derive NAT Saturation Temperature
!     ... using total H2O - assuming NAT forms before ICE
!----------------------------------------------------------------------
      call nat_sat_temp( ncol, hno3_total, h2o_total, press, tsat_nat )

!----------------------------------------------------------------------
!     ... Set SAD to SAGEII if Temperature is GT 200K and Return
!     ... mask2 = true  .... T > 200K or SAD_SULF < 1e-15 or P >2hPa or P <300hPa
!     ... mask1 = false .... T <= TSAT_NAT
!     ... mask  = false .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
         mask2(:,k) = temp(:ncol,k) > 200._r8 .or. sad_sage(:,k) <= 1.e-15_r8 &
                      .or. press(:ncol,k) < 2._r8 .or. press(:ncol,k) > 300._r8
      end do

sage_sad : &
      if( any( mask2(:,sad_topp:pver) ) ) then
         do k = sad_topp,pver
            where( mask2(:,k) )
	       sad_strat(:,k,1)    = sad_sage(:,k)
	       sad_strat(:,k,2)    = 0._r8
	       sad_strat(:,k,3)    = 0._r8
               sad_strat(:,k,4)    = 0._r8
            endwhere
         end do
!----------------------------------------------------------------------
!     ... Calculate Liquid Binary Sulfate Radius
!----------------------------------------------------------------------
         call calc_radius_lbs( ncol, mask2, sad_sage, radius_lbs )
         do k = sad_topp,pver
            where( mask2(:,k) )
               radius_strat(:,k,1)     = radius_lbs(:,k)
	       radius_strat(:,k,2)     = 0._r8
	       radius_strat(:,k,3)     = 0._r8
               radius_strat(:,k,4)     = 0._r8
	       hno3_gas(:,k)           = hno3_total(:,k)
               hno3_cond(:,k,1)        = 0._r8
	       hno3_cond(:,k,2)        = 0._r8
	       hno3_cond(:,k,3)        = 0._r8
            endwhere
         end do
         if( all( mask2(:,sad_topp:pver) ) ) then
            call outfld( 'H2SO4M_C', h2so4m(:ncol,:), ncol, lchnk )

            hno3_cond_arg(:,:,1) = hno3_cond(:,:,1)
            hno3_cond_arg(:,:,2) = hno3_cond(:,:,3)

            radius_strat_arg(:,:,1) = radius_strat(:,:,1)
            radius_strat_arg(:,:,2) = radius_strat(:,:,3)
            radius_strat_arg(:,:,3) = radius_strat(:,:,4)

            sad_strat_arg(:,:,1) = sad_strat(:,:,1)
            sad_strat_arg(:,:,2) = sad_strat(:,:,3)
            sad_strat_arg(:,:,3) = sad_strat(:,:,4)

	    return
         end if
      end if sage_sad
!----------------------------------------------------------------------
!     ... Logic for deriving sad for sulfate, nat, and ice aerosol
!         Ice formation occurs here if condensed phase H2O exists (from CAM).
!         mask2 = false.... T > 200K or SAD_SULF < 1e-15 or P >2hPa or P <300hPa
!         mask1 = true .... T <= TSAT_NAT
!         mask  = true .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
	 do i = 1,ncol
	    if( .not. mask2(i,k) ) then
              mask(i,k) = h2o_cond(i,k) > 0._r8
	    else
	       mask(i,k) = .false.
	    end if
         end do
      end do

all_sad : &
      if( any( mask(:,sad_topp:pver) ) ) then
         do k = sad_topp,pver
            where( mask(:,k) )
	       h2o_avail(:,k) = h2o_cond(:,k)
            endwhere
         end do

         call ice_sad_calc( ncol, press, temp, m, h2o_avail, &
			    sad_ice, radius_ice, mask )

         do k = sad_topp,pver
            where( mask(:,k) )
               h2o_avail(:,k)  = h2o_gas(:,k)
               hno3_avail(:,k) = hno3_total(:,k)
            endwhere
	    if( any(mask(:,k)) ) then
	       where( mask(:,k) )
	          z_val(:) = hno3_avail(:,k) == 0._r8
	       elsewhere
	          z_val(:) = .false.
	       endwhere
	       if( any( z_val(:) ) ) then
	          write(iulog,*) 'sad_strat_calc: Before NAT_SAD_CALC_1 has zero hno3_avail at lchnk,k = ',lchnk,k
	       end if
	    end if
         end do

         call nat_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, m, &
			    hno3_gas_nat, hno3_cond_nat_sm, hno3_cond_nat_lg, &
                            sad_nat_sm, sad_nat_lg, &
			    radius_nat_sm, radius_nat_lg, mask )

!----------------------------------------------------------------------
!  	... H2O Gas available for forming STS
!           The HNO3 avail for the STS aerosol is the amount in the gas
!           phase, plus any additional HNO3 not contained in the large
!           NAT mode at 6.5um.
!----------------------------------------------------------------------
         do k = sad_topp,pver
            where( mask(:,k) )
               h2o_avail(:,k)   = h2o_gas(:,k)
               hno3_avail(:,k)  = hno3_gas_nat(:,k)+hno3_cond_nat_sm(:,k)
            endwhere
	    if( any(mask(:,k)) ) then
	       where( mask(:,k) )
	          z_val(:) = hno3_avail(:,k) == 0._r8
	       elsewhere
	          z_val(:) = .false.
	       endwhere
	       if( any( z_val(:) ) ) then
	          write(iulog,*) 'sad_strat_calc: After NAT_SAD_CALC_1 has zero hno3_avail at lchnk,k = ',lchnk,k
	       end if
	    end if
         end do

         call sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, &
                                sad_sage, m, hno3_gas_sulf, hno3_cond_sulf, sad_sulfate, &
                                radius_sulfate, mask, lchnk, 1, h2so4m, .true. )
         do k = sad_topp,pver
            where( mask(:,k) )
               sad_strat(:,k,1)    = sad_sulfate(:,k)
               sad_strat(:,k,2)    = 0.0_r8
               sad_strat(:,k,3)    = sad_nat_lg(:,k)
               sad_strat(:,k,4)    = sad_ice(:,k)
               radius_strat(:,k,1) = radius_sulfate(:,k)
               radius_strat(:,k,2) = 0.0_r8
               radius_strat(:,k,3) = radius_nat_lg(:,k)
               radius_strat(:,k,4) = radius_ice(:,k)
               hno3_gas(:,k)       = hno3_gas_sulf(:,k)
               hno3_cond(:,k,1)    = hno3_cond_sulf(:,k)  ! de-NOx the atmosphere
               hno3_cond(:,k,2)    = 0.0_r8
               hno3_cond(:,k,3)    = hno3_cond_nat_lg(:,k) ! de-NOY the atmosphere
            endwhere
         end do
      end if all_sad
!----------------------------------------------------------------------
!  	... Saturation of NAT occurs here it temp LE NAT saturation temp.
!           There is no condensed phase H2O (mask is false)
!           mask2 = false .... T > 200K or SAD_SULF < 1e-15 or P >2hPa or P <300hPa
!           mask1 = true  .... T <= TSAT_NAT
!           mask  = false .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
	 do i = 1,ncol
	    if( .not. mask2(i,k) .and. .not. mask(i,k) ) then
               mask1(i,k)   = temp(i,k) <= tsat_nat(i,k)
	    else
	       mask1(i,k) = .false.
	    end if
         end do
      end do

nat_sad : &
      if( any( mask1(:,sad_topp:pver) ) ) then

         do k = sad_topp,pver
            where( mask1(:,k) )
               h2o_avail(:,k)  = h2o_gas(:,k)
               hno3_avail(:,k) = hno3_total(:,k)
            endwhere
	    if( any(mask1(:,k)) ) then
	       where( mask1(:,k) )
	          z_val(:) = hno3_avail(:,k) == 0._r8
	       elsewhere
	          z_val(:) = .false.
	       endwhere
	       if( any( z_val(:) ) ) then
	          write(iulog,*) 'sad_strat_calc: Before NAT_SAD_CALC_2 has zero hno3_avail at lchnk,k = ',lchnk,k
	       end if
	    end if
         end do

         call nat_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, m, &
			    hno3_gas_nat, hno3_cond_nat_sm, hno3_cond_nat_lg, &
                            sad_nat_sm, sad_nat_lg, &
			    radius_nat_sm, radius_nat_lg, mask1 )
!----------------------------------------------------------------------
!  	... H2O Gas available for forming STS
!           The HNO3 avail for the STS aerosol is the amount in the gas
!           phase, plus any additional HNO3 not contained in the large
!           NAT mode at 6.5um.
!----------------------------------------------------------------------
         do k = sad_topp,pver
            where( mask1(:,k) )
               h2o_avail(:,k)    = h2o_gas(:,k)
               hno3_avail(:,k)   = hno3_gas_nat(:,k)+hno3_cond_nat_sm(:,k)
            endwhere
	    if( any(mask1(:,k)) ) then
	       where( mask1(:,k) )
	          z_val(:) = hno3_avail(:,k) == 0._r8
	       elsewhere
	          z_val(:) = .false.
	       endwhere
	       if( any( z_val(:) ) ) then
	          write(iulog,*) 'sad_strat_calc: After NAT_SAD_CALC_2 has zero hno3_avail at lchnk,k = ',lchnk,k
	       end if
	    end if
         end do
         call sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, &
                                sad_sage, m, hno3_gas_sulf, hno3_cond_sulf, sad_sulfate, &
                                radius_sulfate, mask1, lchnk, 2, h2so4m, .true. )
         do k = sad_topp,pver
            where( mask1(:,k) )
               sad_strat(:,k,1)    = sad_sulfate(:,k)
               sad_strat(:,k,2)    = 0._r8
               sad_strat(:,k,3)    = sad_nat_lg(:,k)
               sad_strat(:,k,4)    = 0._r8
               radius_strat(:,k,1) = radius_sulfate(:,k)
               radius_strat(:,k,2) = 0._r8
               radius_strat(:,k,3) = radius_nat_lg(:,k)
               radius_strat(:,k,4) = 0._r8
               hno3_gas(:,k)       = hno3_gas_sulf(:,k)
               hno3_cond(:,k,1)    = hno3_cond_sulf(:,k)  !de-NOx atmosphere
               hno3_cond(:,k,2)    = 0._r8
               hno3_cond(:,k,3)    = hno3_cond_nat_lg(:,k) !de-NOy atmosphere
            endwhere
         end do
      end if nat_sad
!----------------------------------------------------------------------
!  	... Sulfate Ternary (only)
!           No condensed phase H2O (mask = false)
!           No condensed phase NAT (mask1 = false)
!           mask2 = false .... T > 200K or SAD_SULF < 1e-15 or P >2hPa or P <300hPa
!           mask1 = false .... T <= TSAT_NAT
!           mask  = false .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
         mask(:,k) = .not. mask(:,k) .and. .not. mask1(:,k) .and. .not. mask2(:,k)
      end do

sulfate_sad : &
      if( any( mask(:,sad_topp:pver) ) ) then
         do k = sad_topp,pver
            where( mask(:,k) )
               h2o_avail(:,k)  = h2o_gas(:,k)
               hno3_avail(:,k) = hno3_total(:,k)
            endwhere
	    if( any(mask(:,k)) ) then
	       where( mask(:,k) )
	          z_val(:) = hno3_avail(:,k) == 0._r8
	       elsewhere
	          z_val(:) = .false.
	       endwhere
	       if( any( z_val(:) ) ) then
	          write(iulog,*) 'sad_strat_calc: Before NAT_SAD_CALC_3 has zero hno3_avail at lchnk,k = ',lchnk,k
	       end if
	    end if
         end do

         call sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, &
                                sad_sage, m, hno3_gas_sulf, hno3_cond_sulf, sad_sulfate, &
                                radius_sulfate, mask, lchnk, 3, h2so4m, .true. )
         do k = sad_topp,pver
            where( mask(:,k) )
               sad_strat(:,k,1)    = sad_sulfate(:,k)
               sad_strat(:,k,2)    = 0._r8
               sad_strat(:,k,3)    = 0._r8
               sad_strat(:,k,4)    = 0._r8
               radius_strat(:,k,1) = radius_sulfate(:,k)
               radius_strat(:,k,2) = 0._r8
               radius_strat(:,k,3) = 0._r8
               radius_strat(:,k,4) = 0._r8
               hno3_gas(:,k)       = hno3_gas_sulf(:,k)
               hno3_cond(:,k,1)    = hno3_cond_sulf(:,k)
               hno3_cond(:,k,2)    = 0._r8
               hno3_cond(:,k,3)    = 0._r8
            endwhere
         end do
      end if sulfate_sad

      call outfld( 'H2SO4M_C', h2so4m(:ncol,:), ncol, lchnk )
      
      hno3_cond_arg(:,:,1) = hno3_cond(:,:,1)
      hno3_cond_arg(:,:,2) = hno3_cond(:,:,3)

      radius_strat_arg(:,:,1) = radius_strat(:,:,1)
      radius_strat_arg(:,:,2) = radius_strat(:,:,3)
      radius_strat_arg(:,:,3) = radius_strat(:,:,4)

      sad_strat_arg(:,:,1) = sad_strat(:,:,1)
      sad_strat_arg(:,:,2) = sad_strat(:,:,3)
      sad_strat_arg(:,:,3) = sad_strat(:,:,4)

      end subroutine sad_strat_calc

      subroutine sad_strat_calc_rad( lchnk, ncol, m, press, temper, &
                                     hno3_gas_in, hno3_cond_in, h2o_gas_in, h2o_cond_in, sad_sage, &
                                     pbuf )

      use physconst,   only : avogad,boltz,mwdry
      use cam_history, only : outfld
      use physics_buffer, only : physics_buffer_desc, pbuf_get_field
      use cam_pio_utils, only : fillvalue

      implicit none

!-------------------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------------------
      integer, intent(in)     ::  lchnk                      ! chnk id
      integer, intent(in)     ::  ncol                       ! columns in chunk
      real(r8), intent(in)    ::  m(ncol,pver)               ! Air density (molecules cm-3)
      real(r8), intent(in)    ::  sad_sage(ncol,pver)        ! SAGEII surface area density (cm2 aer. cm-3 air)
      real(r8), intent(in)    ::  press(ncol,pver)           ! Pressure, hPa
      real(r8), intent(in)    ::  temper(pcols,pver)         ! Temperature (K)
      real(r8), intent(in)    ::  h2o_gas_in(ncol,pver)      ! H2O gas-phase (mole fraction)
      real(r8), intent(in)    ::  h2o_cond_in(ncol,pver)     ! H2O condensed phase (mole fraction)
      real(r8), intent(in)    ::  hno3_gas_in(ncol,pver)     ! HNO3 condensed phase (mole fraction)
      real(r8), intent(in)    ::  hno3_cond_in(ncol,pver,3)  ! HNO3 condensed phase (mole fraction)
 
!-------------------------------------------------------------------------------
!	... local variables
!-------------------------------------------------------------------------------
      real(r8), parameter :: temp_floor = 0._r8

      integer  ::  i, k, n
      integer  ::  dims(1)
      real(r8) ::  hno3_total(ncol,pver)         ! HNO3 total = gas-phase + condensed
      real(r8) ::  h2o_total(ncol,pver)          ! H2O total  = gas-phase + condensed
      real(r8) ::  radius_lbs(ncol,pver)         ! Radius of Liquid Binary Sulfate (cm)
      real(r8) ::  radius_sulfate(ncol,pver)     ! Radius of Sulfate aerosol (cm)
      real(r8) ::  radius_nat_sm(ncol,pver)      ! Radius of NAT aerosol at 0.5um  (cm)
      real(r8) ::  radius_nat_lg(ncol,pver)      ! Radius of NAT aerosol at 6.5um  (cm)
      real(r8) ::  radius_ice(ncol,pver)         ! Radius of ICE aerosol     (cm)
      real(r8) ::  sad_nat_sm(ncol,pver)         ! SAD of NAT aerosol at 0.5um (cm2 cm-3)
      real(r8) ::  sad_nat_lg(ncol,pver)         ! SAD of NAT aerosol at 6.5um (cm2 cm-3)
      real(r8) ::  sad_sulfate(ncol,pver)        ! SAD of Sulfate aerosol    (cm2 cm-3)
      real(r8) ::  sad_ice(ncol,pver)            ! SAD of ICE aerosol        (cm2 cm-3)
      real(r8) ::  tsat_nat(ncol,pver)           ! Temperature for NAT saturation
      real(r8) ::  h2o_avail(ncol,pver)          ! H2O temporary arrays
      real(r8) ::  hno3_avail(ncol,pver)         ! HNO3 temporary array
      real(r8) ::  hno3_gas_nat(ncol,pver)       ! HNO3 after call to NAT routines
      real(r8) ::  hno3_gas_sulf(ncol,pver)      ! HNO3 after call to STS routines
      real(r8) ::  hno3_cond_nat_sm(ncol,pver)   ! HNO3 condensed after call to NAT, small mode
      real(r8) ::  hno3_cond_nat_lg(ncol,pver)   ! HNO3 condensed after call to NAT, large mode
      real(r8) ::  hno3_cond_sulf(ncol,pver)     ! HNO3 condensed after call to STS routines
      real(r8) ::  temp(pcols,pver)              ! wrk temperature array
      real(r8) ::  h2o_gas(ncol,pver)            ! H2O gas-phase (mole fraction)
      real(r8) ::  h2o_cond(ncol,pver)           ! H2O condensed phase (mole fraction)
      real(r8) ::  hno3_gas(ncol,pver)           ! HNO3 condensed phase (mole fraction)
      real(r8) ::  hno3_cond(ncol,pver,3)        ! HNO3 condensed phase (mole fraction)
      real(r8) ::  radius_strat(ncol,pver,4)     ! Radius of Sulfate, NAT_sm, NAT_lg, and ICE (cm)
      real(r8) ::  sad_strat(ncol,pver,4)        ! Surface area density of Sulfate, NAT_sm, NAT_lg, ICE (cm2 cm-3)

      logical  ::  z_val(ncol)
      logical  ::  mask(ncol,pver)
      logical  ::  mask1(ncol,pver)
      logical  ::  mask2(ncol,pver)
      real(r8) ::  tmpfld(ncol,pver)

      real(r8), pointer, dimension(:,:) :: h2so4m
      real(r8), pointer, dimension(:,:) :: r_g_sulfate
      real(r8) dens_fctr      ! = 1.e-11 * boltz * avogad / mwdry

      type(physics_buffer_desc), pointer :: pbuf(:)
      dens_fctr = 1.e-11_r8 * boltz * avogad / mwdry

      call pbuf_get_field(pbuf, h2so4_ndx, h2so4m )
      call pbuf_get_field(pbuf, h2so4_r_g_ndx, r_g_sulfate )
   
!----------------------------------------------------------------------
!     ... initialize to zero
!----------------------------------------------------------------------
      do n = 1,4
         do k = 1,pver
            radius_strat(:,k,n) = 0._r8
            sad_strat(:,k,n)    = 0._r8
         end do
      end do

!----------------------------------------------------------------------
!     ... initialize h2so4m, copy incoming variables to local space
!----------------------------------------------------------------------
      do k = 1,pver
         h2so4m(:ncol,k)  = 0._r8
         r_g_sulfate(:ncol,k) = 1._r8 
         h2o_gas(:ncol,k) = h2o_gas_in(:ncol,k)
         h2o_cond(:ncol,k) = h2o_cond_in(:ncol,k)
         hno3_gas(:ncol,k) = hno3_gas_in(:ncol,k)
      end do

      do n = 1,3
         do k = 1,pver
            hno3_cond(:ncol,k,n) = hno3_cond_in(:ncol,k,n)
         end do
      end do
!----------------------------------------------------------------------
!     ... limit temperature
!----------------------------------------------------------------------
      do k = 1,pver
         temp(:ncol,k) = max( temp_floor,temper(:ncol,k) )
      end do

!----------------------------------------------------------------------
!     ... total HNO3 and H2O gas and condensed phases
!----------------------------------------------------------------------
      do k = sad_topp,pver
         hno3_total(:,k) = hno3_gas(:,k) + hno3_cond(:,k,1) + hno3_cond(:,k,2) + hno3_cond(:,k,3)
	 h2o_total(:,k)  = h2o_gas(:,k)  + h2o_cond(:,k)
      end do

!----------------------------------------------------------------------
!     ... Derive NAT Saturation Temperature
!     ... using total H2O - assuming NAT forms before ICE
!----------------------------------------------------------------------
      call nat_sat_temp( ncol, hno3_total, h2o_total, press, tsat_nat )

!----------------------------------------------------------------------
!     ... Set SAD to SAGEII primary mask
!     ... mask2 = true  .... 180 <= T && T <= 350K && SAD_SULF > 1e-15 && P >= 2hPa && P <= 300hPa
!     ... mask1 = false .... T <= TSAT_NAT
!     ... mask  = false .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
         mask2(:,k) = 180._r8 <= temp(:ncol,k) .and. temp(:ncol,k) <= 350._r8 &
                      .and. sad_sage(:,k) > 1.e-15_r8 &
                      .and. press(:ncol,k) >= 2._r8 .and. press(:ncol,k) <= 300._r8
      end do
!----------------------------------------------------------------------
!     ... Logic for deriving sad for sulfate, nat, and ice aerosol
!         Ice formation occurs here if condensed phase H2O exists (from CAM).
!         mask2 = false.... T > 200K or SAD_SULF < 1e-15 or P >2hPa or P <300hPa
!         mask1 = true .... T <= TSAT_NAT
!         mask  = true .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
	 do i = 1,ncol
	    if( mask2(i,k) ) then
               mask(i,k) = h2o_cond(i,k) > 0._r8
	    else
	       mask(i,k) = .false.
	    end if
         end do
      end do

all_sad : &
      if( any( mask(:,sad_topp:pver) ) ) then
         do k = sad_topp,pver
            where( mask(:,k) )
	       h2o_avail(:,k) = h2o_cond(:,k)
            endwhere
         end do

         call ice_sad_calc( ncol, press, temp, m, h2o_avail, &
			    sad_ice, radius_ice, mask )

         do k = sad_topp,pver
            where( mask(:,k) )
               h2o_avail(:,k)  = h2o_gas(:,k)
               hno3_avail(:,k) = hno3_total(:,k)
            endwhere
	    if( any(mask(:,k)) ) then
	       where( mask(:,k) )
	          z_val(:) = hno3_avail(:,k) == 0._r8
	       elsewhere
	          z_val(:) = .false.
	       endwhere
	       if( any( z_val(:) ) ) then
	          write(iulog,*) 'sad_strat_calc: Before NAT_SAD_CALC_1 has zero hno3_avail at lchnk,k = ',lchnk,k
	       end if
	    end if
         end do

         call nat_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, m, &
			    hno3_gas_nat, hno3_cond_nat_sm, hno3_cond_nat_lg, &
                            sad_nat_sm, sad_nat_lg, &
			    radius_nat_sm, radius_nat_lg, mask )

!----------------------------------------------------------------------
!  	... H2O Gas available for forming STS
!           The HNO3 avail for the STS aerosol is the amount in the gas
!           phase, plus any additional HNO3 not contained in the large
!           NAT mode at 6.5um.
!----------------------------------------------------------------------
         do k = sad_topp,pver
            where( mask(:,k) )
               h2o_avail(:,k)   = h2o_gas(:,k)
               hno3_avail(:,k)  = hno3_gas_nat(:,k)+hno3_cond_nat_sm(:,k)
            endwhere
	    if( any(mask(:,k)) ) then
	       where( mask(:,k) )
	          z_val(:) = hno3_avail(:,k) == 0._r8
	       elsewhere
	          z_val(:) = .false.
	       endwhere
	       if( any( z_val(:) ) ) then
	          write(iulog,*) 'sad_strat_calc_rad: After NAT_SAD_CALC_1 has zero hno3_avail at lchnk,k = ',lchnk,k
	       end if
	    end if
         end do

         call sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, &
                                sad_sage, m, hno3_gas_sulf, hno3_cond_sulf, sad_sulfate, &
                                radius_sulfate, mask, lchnk, 1, h2so4m(:ncol,:), .false., r_g_sulfate=r_g_sulfate(:ncol,:) )
         do k = sad_topp,pver
            where( mask(:,k) )
               sad_strat(:,k,1)    = sad_sulfate(:,k)
               sad_strat(:,k,2)    = 0.0_r8
               sad_strat(:,k,3)    = sad_nat_lg(:,k)
               sad_strat(:,k,4)    = sad_ice(:,k)
               radius_strat(:,k,1) = radius_sulfate(:,k)
               radius_strat(:,k,2) = 0.0_r8
               radius_strat(:,k,3) = radius_nat_lg(:,k)
               radius_strat(:,k,4) = radius_ice(:,k)
               hno3_gas(:,k)       = hno3_gas_sulf(:,k)
               hno3_cond(:,k,1)    = hno3_cond_sulf(:,k)  ! de-NOx the atmosphere
               hno3_cond(:,k,2)    = 0.0_r8
               hno3_cond(:,k,3)    = hno3_cond_nat_lg(:,k) ! de-NOY the atmosphere
            endwhere
         end do
      end if all_sad
!----------------------------------------------------------------------
!  	... Saturation of NAT occurs here it temp LE NAT saturation temp.
!           There is no condensed phase H2O (mask is false)
!           mask2 = false .... T > 200K or SAD_SULF < 1e-15 or P >2hPa or P <300hPa
!           mask1 = true  .... T <= TSAT_NAT
!           mask  = false .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
	 do i = 1,ncol
	    if( mask2(i,k) .and. .not. mask(i,k) ) then
               mask1(i,k)   = temp(i,k) <= tsat_nat(i,k)
	    else
	       mask1(i,k) = .false.
	    end if
         end do
      end do

nat_sad : &
      if( any( mask1(:,sad_topp:pver) ) ) then
         do k = sad_topp,pver
            where( mask1(:,k) )
               h2o_avail(:,k)  = h2o_gas(:,k)
               hno3_avail(:,k) = hno3_total(:,k)
            endwhere
	    if( any(mask1(:,k)) ) then
	       where( mask1(:,k) )
	          z_val(:) = hno3_avail(:,k) == 0._r8
	       elsewhere
	          z_val(:) = .false.
	       endwhere
	       if( any( z_val(:) ) ) then
	          write(iulog,*) 'sad_strat_calc: Before NAT_SAD_CALC_2 has zero hno3_avail at lchnk,k = ',lchnk,k
	       end if
	    end if
         end do

         call nat_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, m, &
			    hno3_gas_nat, hno3_cond_nat_sm, hno3_cond_nat_lg, &
                            sad_nat_sm, sad_nat_lg, &
			    radius_nat_sm, radius_nat_lg, mask1 )
!----------------------------------------------------------------------
!  	... H2O Gas available for forming STS
!           The HNO3 avail for the STS aerosol is the amount in the gas
!           phase, plus any additional HNO3 not contained in the large
!           NAT mode at 6.5um.
!----------------------------------------------------------------------
         do k = sad_topp,pver
            where( mask1(:,k) )
               h2o_avail(:,k)    = h2o_gas(:,k)
               hno3_avail(:,k)   = hno3_gas_nat(:,k)+hno3_cond_nat_sm(:,k)
            endwhere
	    if( any(mask1(:,k)) ) then
	       where( mask1(:,k) )
	          z_val(:) = hno3_avail(:,k) == 0._r8
	       elsewhere
	          z_val(:) = .false.
	       endwhere
	       if( any( z_val(:) ) ) then
	          write(iulog,*) 'sad_strat_calc_rad: After NAT_SAD_CALC_2 has zero hno3_avail at lchnk,k = ',lchnk,k
	       end if
	    end if
         end do
         call sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, &
                                sad_sage, m, hno3_gas_sulf, hno3_cond_sulf, sad_sulfate, &
                                radius_sulfate, mask1, lchnk, 2, h2so4m(:ncol,:), .false., r_g_sulfate=r_g_sulfate(:ncol,:) )
         do k = sad_topp,pver
            where( mask1(:,k) )
               sad_strat(:,k,1)    = sad_sulfate(:,k)
               sad_strat(:,k,2)    = 0._r8
               sad_strat(:,k,3)    = sad_nat_lg(:,k)
               sad_strat(:,k,4)    = 0._r8
               radius_strat(:,k,1) = radius_sulfate(:,k)
               radius_strat(:,k,2) = 0._r8
               radius_strat(:,k,3) = radius_nat_lg(:,k)
               radius_strat(:,k,4) = 0._r8
               hno3_gas(:,k)       = hno3_gas_sulf(:,k)
               hno3_cond(:,k,1)    = hno3_cond_sulf(:,k)  !de-NOx atmosphere
               hno3_cond(:,k,2)    = 0._r8
               hno3_cond(:,k,3)    = hno3_cond_nat_lg(:,k) !de-NOy atmosphere
            endwhere
         end do
      end if nat_sad
!----------------------------------------------------------------------
!  	... Sulfate Ternary (only)
!           No condensed phase H2O (mask = false)
!           No condensed phase NAT (mask1 = false)
!           mask2 = false .... T > 200K or SAD_SULF < 1e-15 or P >2hPa or P <300hPa
!           mask1 = false .... T <= TSAT_NAT
!           mask  = false .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
         mask(:,k) = mask2(:,k) .and. .not. mask(:,k) .and. .not. mask1(:,k)
      end do

sulfate_sad : &
      if( any( mask(:,sad_topp:pver) ) ) then
         do k = sad_topp,pver
            where( mask(:,k) )
               h2o_avail(:,k)  = h2o_gas(:,k)
               hno3_avail(:,k) = hno3_total(:,k)
            endwhere
	    if( any(mask(:,k)) ) then
	       where( mask(:,k) )
	          z_val(:) = hno3_avail(:,k) == 0._r8
	       elsewhere
	          z_val(:) = .false.
	       endwhere
	       if( any( z_val(:) ) ) then
	          write(iulog,*) 'sad_strat_calc_rad: Before NAT_SAD_CALC_3 has zero hno3_avail at lchnk,k = ',lchnk,k
	       end if
	    end if
         end do

         call sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, &
                                sad_sage, m, hno3_gas_sulf, hno3_cond_sulf, sad_sulfate, &
                                radius_sulfate, mask, lchnk, 3, h2so4m(:ncol,:), .false., r_g_sulfate=r_g_sulfate(:ncol,:) )
         do k = sad_topp,pver
            where( mask(:,k) )
               sad_strat(:,k,1)    = sad_sulfate(:,k)
               sad_strat(:,k,2)    = 0._r8
               sad_strat(:,k,3)    = 0._r8
               sad_strat(:,k,4)    = 0._r8
               radius_strat(:,k,1) = radius_sulfate(:,k)
               radius_strat(:,k,2) = 0._r8
               radius_strat(:,k,3) = 0._r8
               radius_strat(:,k,4) = 0._r8
               hno3_gas(:,k)       = hno3_gas_sulf(:,k)
               hno3_cond(:,k,1)    = hno3_cond_sulf(:,k)
               hno3_cond(:,k,2)    = 0._r8
               hno3_cond(:,k,3)    = 0._r8
            endwhere
         end do
      end if sulfate_sad

      call outfld( 'SAD_SULFR', sad_strat(:ncol,:,1), ncol, lchnk )
      call outfld( 'RAD_SULFR', radius_strat(:ncol,:,1), ncol, lchnk )
      call outfld( 'H2SO4M_R', h2so4m(:ncol,:), ncol, lchnk )

      ! need mmr and geometric mean radius for radiation
      ! assume aerosol is 75% H2SO4, 25% H2O  hence factor of (4/3) in conversion below
      !  The assumption of 75%/25% weight ratio is part of assumptions in optics
      ! convert: (micrograms dry H2SO4) / m^3 --> kg (wet aerosol) / kg (dry air)
      ! geometric_mean_radius = effective_radius * exp(-5*sigma^2/2)

      do k = 1,pver
!         h2so4m(:ncol,k) = (4._r8/3._r8)*h2so4m(:ncol,k) * dens_fctr * temper(:ncol,k) / press(:ncol,k)
! leave this unchanged for waccm_mozart_v1
         h2so4m(:ncol,k) = h2so4m(:ncol,k) * dens_fctr * temper(:ncol,k) / press(:ncol,k)
      enddo

      call outfld( 'H2SO4MMR', h2so4m(:ncol,:), ncol, lchnk )
      tmpfld(:,:) = fillvalue
      where ( mask(:ncol,:) ) 
         tmpfld(:ncol,:) = r_g_sulfate(:ncol,:)
      endwhere
      call outfld( 'VOLC_RAD_GEOM', tmpfld(:ncol,:), ncol, lchnk)

      end subroutine sad_strat_calc_rad

      subroutine nat_sat_temp( ncol, hno3_total, h2o_avail, press, tsat_nat )

      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      real(r8), intent(in)  :: press(ncol,pver)
      real(r8), intent(in)  :: h2o_avail(ncol,pver)
      real(r8), intent(in)  :: hno3_total(ncol,pver)
      real(r8), intent(out) :: tsat_nat(ncol,pver)

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      real(r8), parameter :: ssrNAT = 10.0_r8
      real(r8), parameter :: ssrNATi = .1_r8
      real(r8), parameter :: aa = -2.7836_r8, &
                             bb = -0.00088_r8, &
                             cc = 38.9855_r8, &
                             dd = -11397.0_r8, &
                             ee = 0.009179_r8
      integer  :: k
      real(r8) :: bbb(ncol)                     ! temporary variable
      real(r8) :: ccc(ncol)                     ! temporary variable
      real(r8) :: wrk(ncol)                     ! temporary variable
      real(r8) :: tmp(ncol)                     ! temporary variable
      real(r8) :: phno3(ncol)                   ! hno3 partial pressure
      real(r8) :: ph2o(ncol)                    ! h2o  partial pressure

!----------------------------------------------------------------------
!     	... Derive HNO3 and H2O partial pressure (torr)
!          where: 0.7501 = 760/1013.
!----------------------------------------------------------------------
      do k = sad_topp,pver
         bbb(:)   = press(:,k) * .7501_r8
         phno3(:) = hno3_total(:,k) * bbb(:)
         ph2o(:)  = h2o_avail(:,k)  * bbb(:)
	 where( phno3(:) > 0._r8 )
!----------------------------------------------------------------------
!     	... Calculate NAT Saturation Threshold Temperature
!           Hanson and Mauersberger: GRL, Vol.15, 8, p855-858, 1988.
!           Substitute m(T) and b(T) into Equation (1). Rearrange and
!           solve quadratic eqation. 
!----------------------------------------------------------------------
            tmp(:) = log10( ph2o(:) )
            wrk(:) = 1._r8 / (ee + bb*tmp(:))
            bbb(:) = (aa*tmp(:) - log10( phno3(:)*ssrNATi ) + cc) * wrk(:)
            ccc(:) = dd *wrk(:)
            tsat_nat(:,k) = .5_r8 * (-bbb(:) + sqrt( bbb(:)*bbb(:) - 4._r8*ccc(:) ))
	 elsewhere
            tsat_nat(:,k) = 0._r8
	 endwhere
      end do

      end subroutine nat_sat_temp

      subroutine calc_radius_lbs( ncol, mask, sad_sage, radius_lbs )

      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      real(r8), intent(in)  :: sad_sage(ncol,pver)
      real(r8), intent(out) :: radius_lbs(ncol,pver)
      logical, intent(in)   :: mask(ncol,pver)

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      integer  :: k
      real(r8) :: lbs_vol_dens(ncol)       ! Vol Density (cm3 aer / cm3 air)

!----------------------------------------------------------------------
!     	... parameters
!----------------------------------------------------------------------
      real(r8), parameter :: lbs_part_dens = 10._r8, sigma_lbs = 1.6_r8

!----------------------------------------------------------------------
!     	... calculate the volume density (cm3 aerosol / cm3 air)
!     	    calculate the mean radius for binary soln
!----------------------------------------------------------------------
      do k = sad_topp,pver
	 where( mask(:,k) )
            lbs_vol_dens(:) = ((sad_sage(:,k)**1.5_r8)/3._r8)/sqrt( four_pi*lbs_part_dens ) &
                              *exp( 1.5_r8*(log( sigma_lbs ))**2 )
            radius_lbs(:,k) = (3._r8*lbs_vol_dens(:)/(four_pi*lbs_part_dens))**one_thrd &
                              *exp( -1.5_r8*(log( sigma_lbs ))**2 )
	 endwhere
      end do

      end subroutine calc_radius_lbs

      subroutine ice_sad_calc( ncol, press, temp, m, h2o_avail, &
			       sad_ice, radius_ice, mask )

      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      real(r8), intent(in)  :: press(ncol,pver)
      real(r8), intent(in)  :: temp(pcols,pver)
      real(r8), intent(in)  :: m(ncol,pver)
      real(r8), intent(in)  :: h2o_avail(ncol,pver)
      real(r8), intent(out) :: sad_ice(ncol,pver)
      real(r8), intent(out) :: radius_ice(ncol,pver)
      logical, intent(in)   :: mask(ncol,pver)

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      real(r8), parameter :: &
                 avo_num = 6.02214e23_r8, &
                 aconst  = -2663.5_r8, &
                 bconst  = 12.537_r8, &
                 ice_mass_dens = 1._r8, &
                 ice_part_dens = 1.e-3_r8, &
                 mwh2o         = 18._r8, &
                 sigma_ice     = 1.6_r8, &
                 ice_dens_aer  = ice_mass_dens / (mwh2o/avo_num), &
                 ice_dens_aeri = 1._r8/ice_dens_aer

      integer  :: k
      real(r8) :: h2o_cond_ice(ncol)              ! Condensed phase H2O (from CAM)
      real(r8) :: voldens_ice(ncol)               ! Volume Density, um3 cm-3

      do k = sad_topp,pver
	 where( mask(:,k) )
!----------------------------------------------------------------------
!     .... Convert condensed phase to molecules cm-3 units
!----------------------------------------------------------------------
	   h2o_cond_ice(:) = h2o_avail(:,k) * m(:,k)
!----------------------------------------------------------------------
!     .... ICE volume density .....
!----------------------------------------------------------------------
           voldens_ice(:) = h2o_cond_ice(:)*ice_dens_aeri
!----------------------------------------------------------------------
!     .... Calculate the SAD from log normal distribution .....
!----------------------------------------------------------------------
           sad_ice(:,k) = (four_pi*ice_part_dens)**one_thrd &
                         *(3._r8*voldens_ice(:))**two_thrd &
                         *exp( -(log( sigma_ice ))**2 )
!----------------------------------------------------------------------
!    .... Calculate the radius from log normal distribution .....
!----------------------------------------------------------------------
           radius_ice(:,k) = (3._r8*h2o_cond_ice(:) &
                              /(ice_dens_aer*four_pi*ice_part_dens))**one_thrd &
                             *exp( -1.5_r8*(log( sigma_ice ))**2 )
         endwhere
      end do

      end subroutine ice_sad_calc

      subroutine sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, &
                                   sad_sage, m, hno3_gas, hno3_cond, sad_sulfate, &
                                   radius_sulfate, mask, lchnk, flag, h2so4m, is_chem, r_g_sulfate )

      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      integer, intent(in)   :: lchnk, flag
      real(r8), intent(in)  :: temp(pcols,pver)
      real(r8), intent(in)  :: press(ncol,pver)
      real(r8), intent(in)  :: m(ncol,pver)
      real(r8), intent(in)  :: h2o_avail(ncol,pver)
      real(r8), intent(in)  :: hno3_avail(ncol,pver)
      real(r8), intent(in)  :: sad_sage(ncol,pver)
      real(r8), intent(out) :: hno3_gas(ncol,pver)             ! Gas-phase HNO3, mole fraction
      real(r8), intent(out) :: hno3_cond(ncol,pver)            ! Condensed phase HNO3, mole fraction
      real(r8), intent(out) :: sad_sulfate(ncol,pver)   
      real(r8), intent(out) :: radius_sulfate(ncol,pver)
      real(r8), intent(inout) :: h2so4m(ncol,pver)             ! mass per volume, micro grams m-3
      logical, intent(in)   :: is_chem                         ! chemistry calc switch
      logical, intent(in)   :: mask(ncol,pver)

      real(r8), optional, intent(out) :: r_g_sulfate(ncol,pver)          ! geometric mean radius

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      real(r8), parameter :: t_limit = 200._r8
      real(r8), parameter :: avo_num = 6.02214e23_r8
      real(r8), parameter :: mwh2so4 = 98.076_r8
      real(r8), parameter :: sigma_sulfate     = 1.6_r8
      real(r8), parameter :: sulfate_part_dens = 10._r8

      integer  :: i, k
      integer  :: cnt1, cnt2
      real(r8) :: ratio
      real(r8) :: vals(2)
      real(r8) :: h2so4_aer_dens(ncol,pver)       ! grams cm-3 solution
      real(r8) :: h2so4_cond(ncol,pver)           ! Condensed H2SO4 (moles cm-3 air)
      real(r8) :: sulfate_vol_dens(ncol,pver)     ! Volume Density, um3 cm-3
      real(r8) :: wtf(ncol,pver)                  ! weight fraction of H2SO4 in ternary soln
      real(r8) :: wts(ncol,pver)                  ! weight percent of ternary solution
      real(r8) :: packer(ncol*pver)
      logical  :: do_equil                        ! local mask
      logical  :: mask_lt(ncol,pver)              ! local temperature mask
      logical  :: maskx(ncol,pver)

      do k = sad_topp,pver
         mask_lt(:,k)  = mask(:,k)
      end do
!----------------------------------------------------------------------
!  	... derive H2SO4 (micro grams / m3) from SAGEII SAD
!----------------------------------------------------------------------
      call sad2h2so4( h2o_avail, press, sad_sage, temp, sulfate_vol_dens, & 
                      h2so4_aer_dens, h2so4m, mask, ncol )

!----------------------------------------------------------------------
!  	... limit h2so4m
!----------------------------------------------------------------------
      if( is_chem ) then
         do k = sad_topp,pver
            do i = 1,ncol
               if( mask(i,k) ) then
                  if( h2so4m(i,k) <= 0._r8 ) then
                     h2so4m(i,k) = 1.e-4_r8
                  end if
               end if
            end do
         end do
      else
        do k = sad_topp,pver
            where( mask(:ncol,k) )
               mask_lt(:ncol,k) = temp(:ncol,k) < t_limit
            end where
         end do
      end if
      if( is_chem ) then
         do_equil = .true.
      else
         do_equil = any( mask_lt(:,sad_topp:pver) )
      end if

!----------------------------------------------------------------------
!     .... Calculate the ternary soln volume density
!----------------------------------------------------------------------
      if( do_equil ) then
!        cnt1 = count( mask(:,sad_topp:pver) )
!        cnt2 = count( mask_lt(:,sad_topp:pver) )
!        ratio = 100._r8*real(cnt2,r8)/real(cnt1,r8)   
!        write(iulog,*) 'sulfate_sad_calc: calling equil with is_chem,flag,count,percent = ',is_chem,flag,count(mask_lt(:,sad_topp:pver)),ratio
         call equil( temp, h2so4m, hno3_avail, h2o_avail, press, & 
                     hno3_cond, h2so4_cond, wts, mask_lt, ncol, &
                     lchnk, flag, is_chem )

         do k = sad_topp,pver
	    where( mask_lt(:,k) )
!----------------------------------------------------------------------
!     .... convert h2o, hno3 from moles cm-3 air to molecules cm-3 air
!----------------------------------------------------------------------
               hno3_cond(:,k) = min( hno3_cond(:,k),hno3_avail(:,k) )
               hno3_gas(:,k)  = hno3_avail(:,k) - hno3_cond(:,k)
!----------------------------------------------------------------------
!     .... Derive ternary volume density (cm3 aer / cm3 air)
!----------------------------------------------------------------------
               wtf(:,k) = .01_r8* wts(:,k)
               sulfate_vol_dens(:,k) = h2so4_cond(:,k)*mwh2so4/(wtf(:,k)*h2so4_aer_dens(:,k))
	    endwhere
         end do
!        packer(:cnt2) = pack( sulfate_vol_dens, mask_lt )
!        write(iulog,'(10(1p,g15.7))') packer(:cnt2)
!        maskx(:,sad_topp:pver) = mask_lt(:,sad_topp:pver) .and. sulfate_vol_dens(:,sad_topp:pver) > 0._r8
!        vals(1) = minval( sulfate_vol_dens(:,sad_topp:pver), maskx(:,sad_topp:pver) )
!        vals(2) = maxval( sulfate_vol_dens(:,sad_topp:pver), maskx(:,sad_topp:pver) )
!        write(iulog,*) 'sulfate_sad_calc: min,max sulfate_vol_dens = ',vals(:)
      end if
      do k = sad_topp,pver
	 where( mask(:,k) )
!----------------------------------------------------------------------
!     .... Calculate the SAD (assuming ternary solution)
!----------------------------------------------------------------------
            sad_sulfate(:,k) = (four_pi*sulfate_part_dens)**one_thrd &
                               *(3._r8*sulfate_vol_dens(:,k))**two_thrd &
                               *exp( -(log( sigma_sulfate ))**2 )
!----------------------------------------------------------------------
!     .... Calculate the radius (assuming ternary solution) (in cm?)
!----------------------------------------------------------------------
            radius_sulfate(:,k) = (3._r8*sulfate_vol_dens(:,k) &
                                   /(four_pi*sulfate_part_dens))**one_thrd &
                                  *exp( -1.5_r8*(log( sigma_sulfate ))**2 )
	 endwhere

!----------------------------------------------------------------------
!     .... Calculate the geometric mean radius (assuming ternary solution) 
!          r_g = r_eff*exp(-5*ln(sigma)^2/2)  
!    convert to cm to meters (1/100)
!----------------------------------------------------------------------
         if (present(r_g_sulfate)) then
            where( mask(:,k) )
               r_g_sulfate(:,k) = (3._r8*sulfate_vol_dens(:,k) &
                    /(four_pi*sulfate_part_dens))**one_thrd &
                    *exp( -4._r8*(log( sigma_sulfate ))**2 ) &
                    /100._r8
            endwhere
         endif

      end do

      end subroutine sulfate_sad_calc


      subroutine nat_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, m, &
			       hno3_gas, hno3_cond_sm, hno3_cond_lg, sad_nat_sm, &
                               sad_nat_lg, radius_nat_sm, radius_nat_lg, mask )

      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      real(r8), intent(in)  :: press     (ncol,pver)
      real(r8), intent(in)  :: m         (ncol,pver)
      real(r8), intent(in)  :: temp      (pcols,pver)
      real(r8), intent(in)  :: h2o_avail (ncol,pver)
      real(r8), intent(in)  :: hno3_avail(ncol,pver)
      real(r8), intent(out) :: hno3_cond_sm(ncol,pver)             ! HNO3 in condensed phase (mole fraction)
      real(r8), intent(out) :: hno3_cond_lg(ncol,pver)             ! HNO3 in condensed phase (mole fraction)
      real(r8), intent(out) :: hno3_gas  (ncol,pver)               ! HNO3 in gas-phase (mole fraction)
      real(r8), intent(out) :: sad_nat_sm(ncol,pver)   
      real(r8), intent(out) :: sad_nat_lg(ncol,pver)   
      real(r8), intent(out) :: radius_nat_sm(ncol,pver)   
      real(r8), intent(out) :: radius_nat_lg(ncol,pver)
      logical, intent(in)   :: mask(ncol,pver)                     ! grid mask

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      integer  :: k, i
      real(r8) :: nat_dens_condphase(ncol, pver)      ! Condensed phase NAT, molec cm-3
      real(r8) :: voldens_nat       (ncol, pver)      ! Volume Density,  um3 cm-3
      real(r8) :: hno3_cond_total   (ncol, pver)       ! Total Condensed phase HNO3 

!----------------------------------------------------------------------
!     ... parameters
!----------------------------------------------------------------------
      real(r8), parameter :: avo_num          = 6.02214e23_r8, &
                             nat_mass_dens    = 1.6_r8, &
                             nat_part_dens_lg = 1.0e-2_r8, &
                             mwnat            = 117._r8, &
                             sigma_nat        = 1.6_r8, &
                             nat_dens_aer     = nat_mass_dens / (mwnat/avo_num), &
                             nat_dens_aeri    = 1._r8/nat_dens_aer

!----------------------------------------------------------------------
!     ... Derive HNO3 paritioning (call A. Tabazedeh routine for NAT)
!----------------------------------------------------------------------
      call nat_cond( ncol, press, temp, h2o_avail, hno3_avail, &
		     hno3_gas, hno3_cond_total, mask )

      do k = sad_topp,pver
         do i = 1,ncol
masked :   if( mask(i,k) ) then
            
!----------------------------------------------------------------------
!     .... Set Condensed phase for return arguments
!----------------------------------------------------------------------
            hno3_cond_sm(i,k) = 0.0_r8
            hno3_cond_lg(i,k) = hno3_cond_total(i,k)

!----------------------------------------------------------------------
!     .... Calculated Condensed Phase NAT (i.e. HNO3) in
!            molecules cm-3 of air units.
!----------------------------------------------------------------------
            nat_dens_condphase(i,k) = hno3_cond_total(i,k) * m(i,k)

!----------------------------------------------------------------------
!     .... Calculate the Volume Density
!----------------------------------------------------------------------
            voldens_nat(i,k) = nat_dens_condphase(i,k) * nat_dens_aeri

!----------------------------------------------------------------------
!     .... Calculate the SAD from log normal distribution
!     .... Assuming sigma and nad_part_dens (# particles per cm3 of air)
!----------------------------------------------------------------------
            sad_nat_lg(i,k) = (four_pi*nat_part_dens_lg)**(one_thrd) &
                          *(3._r8*voldens_nat(i,k))**(two_thrd) &
                          *exp( -(log( sigma_nat )**2 ))

            sad_nat_sm(i,k) = 0.0_r8

!----------------------------------------------------------------------
!     .... Calculate the radius of NAT from log normal distribution
!     .... Assuming sigma and nat_part_dens_lg (# particles per cm3 
!     .... of air)
!----------------------------------------------------------------------
            radius_nat_lg(i,k) = (3._r8*nat_dens_condphase(i,k) &
                        /(nat_dens_aer*four_pi*nat_part_dens_lg))**(one_thrd) &
                        *exp( -1.5_r8*(log( sigma_nat ))**2 )

            radius_nat_sm(i,k) = 0.0_r8

           end if masked
         end do
      end do

      end subroutine nat_sad_calc

      subroutine nat_cond( ncol, press, temp, h2o_avail, hno3_avail, &
			   hno3_gas, hno3_cond, mask )

      implicit none

!----------------------------------------------------------------------
! 	... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      real(r8), intent(in)  :: press(ncol,pver)
      real(r8), intent(in)  :: temp(pcols,pver)
      real(r8), intent(in)  :: h2o_avail(ncol,pver)
      real(r8), intent(in)  :: hno3_avail(ncol,pver)
      real(r8), intent(out) :: hno3_gas(ncol,pver)
      real(r8), intent(out) :: hno3_cond(ncol,pver)
      logical, intent(in)   :: mask(ncol,pver)

!----------------------------------------------------------------------
! 	... local variables
!----------------------------------------------------------------------
      real(r8), parameter ::  aa = -2.7836_r8,  &
                              bb = -0.00088_r8, &
                              cc = 38.9855_r8,  &
                              dd = -11397.0_r8, &
                              ee = 0.009179_r8

      integer  :: i, k
      real(r8) :: bt                                 ! temporary variable
      real(r8) :: mt                                 ! temporary variable
      real(r8) :: t                                  ! temporary variable
      real(r8) :: logPhno3                           ! temporary variable
      real(r8) :: phno3                              ! hno3 partial pressure
      real(r8) :: ph2o                               ! h2o  partial pressure
      real(r8) :: phno3_eq                           ! partial pressure above NAT
      real(r8) :: wrk      
      
      do k = sad_topp,pver
	 do i = 1,ncol
!----------------------------------------------------------------------
!     .... Derive HNO3 and H2O partial pressure (torr)
!          where: 0.7501 = 760/1013.
!----------------------------------------------------------------------        
	    if( mask(i,k) ) then
               wrk   = press(i,k) * .7501_r8
               phno3 = hno3_avail(i,k) * wrk
               ph2o  = h2o_avail(i,k)  * wrk
!----------------------------------------------------------------------
!     Calculating the temperature coefficients for the variation of HNO3
!     and H2O vapor pressure (torr) over a trihydrate solution of HNO3/H2O
!     The coefficients are taken from Hanson and Mauersberger: 
!     GRL, Vol.15, 8, p855-858, 1988.
!----------------------------------------------------------------------
	       t   = temp(i,k)
               bt  = cc + dd/t + ee*t
               mt  = aa + bb*t
  
               logphno3 = mt*log10( ph2o ) + bt
               phno3_eq = 10._r8**logphno3

	       if( phno3 > phno3_eq ) then
                  wrk            = 1._r8 / wrk
                  hno3_cond(i,k) = (phno3 - phno3_eq) * wrk
                  hno3_gas(i,k)  = phno3_eq * wrk
	       else
                  hno3_cond(i,k) = 0._r8
                  hno3_gas(i,k)  = hno3_avail(i,k)
	       end if
	    end if
         end do
      end do

      end subroutine nat_cond

      subroutine sad2h2so4( h2o, press, sad_sage, temp, lbs_vol_dens, &
                            h2so4_aer_dens, h2so4m, mask, ncol )

      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol                                    ! columns in chunk
      real(r8), intent(in)  :: h2o(ncol,pver)                          ! h2o mole fraction
      real(r8), intent(in)  :: sad_sage(ncol,pver)                     ! sad from SAGEII cm2 aer, cm-3 air
      real(r8), intent(in)  :: press(ncol,pver)                        ! pressure (hPa)
      real(r8), intent(in)  :: temp(pcols,pver)                        ! temperature (K)
      real(r8), intent(inout) :: h2so4m(ncol,pver)                       ! microgram/m**3 of air,
      real(r8), intent(out) :: h2so4_aer_dens(ncol,pver)               ! units: grams / cm3-aerosol
      real(r8), intent(out) :: lbs_vol_dens(ncol,pver)                 ! cm3 aer / cm3 air
      logical, intent(in)   :: mask(ncol,pver)                         ! activation mask

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      real(r8), parameter :: lbs_part_dens = 10._r8
      real(r8), parameter :: sigma_lbs     = 1.6_r8
      real(r8), parameter :: t_floor       = 180._r8

      integer  :: i, k, l
      real(r8) :: wts0   
      real(r8) :: p                         ! pressure, torr
      real(r8) :: tr                        ! inverse temperature
      real(r8) :: c(6)

!----------------------------------------------------------------------
!     	... Calculate the volume density (cm3 aerosol / cm3 air)
!----------------------------------------------------------------------
      do k = sad_topp,pver
	 do i = 1,ncol
	    if( mask(i,k) ) then
               lbs_vol_dens(i,k) = ((sad_sage(i,k)**1.5_r8)/3._r8)/sqrt( four_pi*lbs_part_dens ) &
                                   *exp( 1.5_r8*(log( sigma_lbs ))**2 )
!----------------------------------------------------------------------
!     	... calculate Molality from Tabazadeh EQUIL routine (binary soln)
!----------------------------------------------------------------------
!     	... DEK, added a minimum to temperature
!----------------------------------------------------------------------
               p   = h2o(i,k) * press(i,k) * .7501_r8
               tr  = 1._r8 / max( t_floor,temp(i,k) )
             
	       do l = 1,6
                  c(l) = exp( a(1,l) + tr*(a(2,l) + tr*(a(3,l) + tr*(a(4,l) + tr*a(5,l)))) )
	       end do
!----------------------------------------------------------------------
!     	... H2SO4/H2O pure weight percent and molality (moles gram-1)
!----------------------------------------------------------------------
               wts0 = max( 0._r8,c(1) + p*(-c(2) + p*(c(3) + p*(-c(4) + p*(c(5) - p*c(6))))) )
!----------------------------------------------------------------------
!     	... derive weight fraction for density routine
!----------------------------------------------------------------------
               wts0 = .01_r8 *wts0
!----------------------------------------------------------------------
!     	... calculate the binary solution density, grams / cm3-aerosol
!----------------------------------------------------------------------
               h2so4_aer_dens(i,k) = max( 0._r8,density( temp(i,k), wts0 ) )
!----------------------------------------------------------------------
!     	... calculate the H2SO4 micrograms m-3 abundance for binary soln
!----------------------------------------------------------------------
               h2so4m(i,k) = lbs_vol_dens(i,k)*h2so4_aer_dens(i,k)*wts0*1.e12_r8
	    end if
         end do
      end do
   
      end subroutine sad2h2so4

!======================================================================
!
! ROUTINE
!   EQUIL
!
!   Date...
!     7 October 1999
!
!   Programmed by...
!     A. Tabazadeh
!
! DESCRIPTION
!	Ternary solution routine
!
! ARGUMENTS
!
!....  INPUT:
!
!        H2SO4m    = microgram/m**3 of air
!        HNO3r     = mole fraction
!        H2Or      = mole fraction
!        PTOTAL    = hPa
!
!....  Output
!
!        Cwater    = Total moles of liguid water / cm**3 of air
!        hno3_cond = HNO3 Condensed phase (mole fraction)
!        CH2SO4    = Total H2SO4 moles / cm**3 of air
!        WTS       = Weight percent of H2SO4 in the ternary aerosol
!
!======================================================================

      subroutine equil( temper, h2so4m, hno3_avail, h2o_avail, press, &
                        hno3_cond, ch2so4, wts, mask, ncol, &
                        lchnk, flag, is_chem )
!----------------------------------------------------------------------
!                       Written by Azadeh Tabazadeh (1993)
!           (modified from EQUISOLV -- M. Z. Jacobson -- see below) 
!            NASA Ames Research Center , Tel. (415) 604 - 1096
!
!       This program solves the equilibrium composition for the ternary  
!       system of H2SO4/HNO3/H2O under typical stratospheric conditions. 
!       The formulation of this work is described by Tabazadeh, A.,      
!       Turco, R. P., and Jacobson, M. Z. (1994), "A model for studying  
!       the composition and chemical effects of stratospheric aerosols," 
!       J. Geophys. Res., 99, 12,897, 1994.        *
!
!       The solution mechanism for the equilibrium equations is des-     
!       cribed by Jacobson, M. Z., Turco, R. P., and Tabazadeh, A.       
!       (1994), "Simulating Equilibrium within aerosols and non-equil-   
!       ibrium between gases and aerosols," J. Geophys. Res., in review. 
!       The mechanism is also codified in the fortran program, EQUISOLV, 
!       by M.Z. Jacobson (1991-3). EQUISOLV solves any number of         
!       gas / liquid / solid / ionic equilibrium equations simultan-     
!       eously and includes treatment of the water equations and act-    
!       ivity coefficients. The activity coeffients currently in         
!       EQUISOLV are valid for tropospheric temperatures. The acitiv-    
!       ities listed here are valid for stratospheric temperatures only. 
!
!	DEFINING PARAMETERS
!
!       *NOTE*	  Solver parameters including, F, Z, QN, QD, and deltaX
!                 are described in Jacobson et al.
!
!	PTOTAL	= Total atmospheric pressure in mb
!	H2SO4m	= Total mass of H2SO4 (microgram/m**3 of air)
!	HNO3r	= HNO3 mixing ratio
!	H2Or	= H2O mixing ratio
!	P	    = Partial pressure of water in units of torr
!	pures	= molality for a pure H2SO4/H2O system
!	puren	= molality for a pure HNO3/H2O sytem
!	WTS0	= weight percent of H2SO4 in a pure H2SO4/H2O system
!	WTN0	= weight percent of HNO3 in a pure HNO3/H2O system
!	WTS	    = weight percent of H2SO4 in the ternary aerosol
!	WTN 	= weight percent of HNO3 in the ternary aerosol
!	PHNO3	= HNO3 vapor pressure over the ternary system in atm
!	HNO3	= HNO3 vapor concentration over the ternary system (#/cm3)
!	CH2SO4	= Total H2SO4 moles / cm**3 of air
!	CHNO3	= Total HNO3 moles / cm**3 of air
!	CHplus	= Total H+ moles / cm**3 0f air
!	CPHNO3	= Total moles of HNO3 gas / cm**3 of air
!	CNO3	= Total moles of NO3- / cm**3 0f air
!	Cwater	= Total moles of liguid water / cm**3 of air
!	KS  	= Solubility constant for dissolution of HNO3 in
!		      water ( HNO3(gas) === H+(aq) + NO3- (aq) )
!	nm  	= HNO3 molality at the STREN of the ternary solution
!	sm  	= H2SO4 molality at the STREN of the ternary solution
!	molHNO3	= Equilibrium molality of HNO3 in the ternary solution
!	molH2SO4= Equilibrium molality of H2SO4 in the ternary solution
!     STREN   = ionic strenght for the ternary solutin, which in
!		      this case is = 3 * molH2SO4 + molHNO3
!	acts	= Pure mean binary activity coefficient for the H2SO4/
!		      H2O system evaluated at the STREN of the ternary system
!	actn	= Pure mean binary activity coefficient for the HNO3/
!		      H2O system evaluated at the STREN of the ternary system
!     ymix    = Mixed binary activity coefficient for the HNO3/H2O in
!		      the ternary solution
!----------------------------------------------------------------------

      use abortutils, only : endrun

      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: lchnk
      integer, intent(in)   :: flag
      integer, intent(in)   :: ncol                         ! columns in chunk
      real(r8), intent(in)  :: h2so4m(ncol,pver)    
      real(r8), intent(in)  :: hno3_avail(ncol,pver)    
      real(r8), intent(in)  :: h2o_avail(ncol,pver)    
      real(r8), intent(in)  :: press(ncol,pver)
      real(r8), intent(in)  :: temper(pcols,pver)
      real(r8), intent(out) :: hno3_cond(ncol,pver)    
      real(r8), intent(out) :: ch2so4(ncol,pver)    
      real(r8), intent(out) :: wts(ncol,pver)
      logical, intent(in)   :: is_chem
      logical, intent(in)   :: mask(ncol,pver)              ! activation mask

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      integer, parameter  :: itermax = 50
      real(r8), parameter :: con_lim  = .00005_r8
      real(r8), parameter :: t0       = 298.15_r8
      real(r8), parameter :: ks0      = 2.45e6_r8
      real(r8), parameter :: lower_delx = 1.e-10_r8
      real(r8), parameter :: upper_delx = .98_r8
      real(r8), parameter :: con_crit_chem = 5.e-5_r8

      integer  :: i, iter, k, l, nstep
      real(r8) :: reduction_factor
      real(r8) :: p
      real(r8) :: tr
      real(r8) :: wts0
      real(r8) :: wtn0
      real(r8) :: pures
      real(r8) :: puren
      real(r8) :: chno3
      real(r8) :: chplus
      real(r8) :: cno3
      real(r8) :: wrk
      real(r8) :: z, num, den
      real(r8) :: deltax
      real(r8) :: chplusnew
      real(r8) :: cno3new
      real(r8) :: stren
      real(r8) :: sm
      real(r8) :: actn
      real(r8) :: acts
      real(r8) :: nm
      real(r8) :: ks
      real(r8) :: lnks
      real(r8) :: lnks0
      real(r8) :: mixyln
      real(r8) :: molhno3
      real(r8) :: molh2so4
      real(r8) :: wrk_h2so4
      real(r8) :: cphno3new
      real(r8) :: con_val
      real(r8) :: t, t1, t2, f, f1, f2, ymix, hplus, wtotal, wtn, ratio 
      real(r8) :: con_crit
      real(r8) :: h2o_cond(ncol,pver)
      real(r8) :: fratio(0:itermax)
      real(r8) :: delx(0:itermax)
      real(r8) :: delz(0:itermax)
      real(r8) :: c(12)
      real(r8) :: d(13:22)
      logical  :: interval_set
      logical  :: positive
      logical  :: converged

      lnks0 = log( ks0 )
      if( is_chem ) then
         con_crit = con_crit_chem
      else
         con_crit = con_crit_chem
      end if
Level_loop : &
      do k = sad_topp,pver
Column_loop : &
	 do i = 1,ncol
	    if( mask(i,k) ) then
               p = h2o_avail(i,k) * press(i,k) * .7501_r8
!----------------------------------------------------------------------
!	Calculating the molality for pure binary systems of H2SO4/H2O
!	and HNO3/H2O at a given temperature and water vapor pressure
!	profile (relative humiditiy). Water activities were used to
!	calculate the molalities as described in Tabazadeh et al. (1994).
!----------------------------------------------------------------------
               t  = max( 180._r8,temper(i,k) )
               tr = 1._r8/t
               do l = 1,12
                  c(l) = exp( a(1,l) + tr*(a(2,l) + tr*(a(3,l) + tr*(a(4,l) + tr*a(5,l)))) )
               end do
!----------------------------------------------------------------------
!	... H2SO4/H2O pure weight percent and molality
!----------------------------------------------------------------------
               wts0  = max( 0._r8,c(1) + p*(-c(2) + p*(c(3) + p*(-c(4) + p*(c(5) - p*c(6))))) )
               pures = (wts0 * 1000._r8)/(100._r8 - wts0)
               pures = pures / 98._r8
!----------------------------------------------------------------------
!	... HNO3/H2O pure weight percent and molality
!----------------------------------------------------------------------
               puren = c(7) + p*(-c(8) + p*(c(9) + p*(-c(10) + p*(c(11) - p*c(12)))))
!              wtn0 = (puren * 6300._r8) /(puren * 63._r8 + 1000._r8)
!----------------------------------------------------------------------
!	The solving scheme is described both in Jacobson et al. and Tabazadeh
!	et al.. Assumptions:
!	(1) H2SO4 is present only in the aqueous-phase
!	(2) H2SO4 and HNO3 in solution are fully dissocated into H+
!	    SO42- and NO3-
!	(3) PHNO3 + NO3- = constant
!----------------------------------------------------------------------
	       ch2so4(i,k) = (h2so4m(i,k)*1.e-12_r8) / 98._r8
	       if( pures > 0._r8 ) then
	          wrk_h2so4 = (1000._r8*ch2so4(i,k))/(pures*18._r8)
	       else
	          wrk_h2so4 = 0._r8
	       end if
	       chno3 = 1.2029e-5_r8 * press(i,k) * tr * hno3_avail(i,k)
	       do l = 13,22
                  d(l) = b(1,l) + t*(b(2,l) + t*(b(3,l) + t*(b(4,l) + t*b(5,l))))
	       end do
!----------------------------------------------------------------------
!	Note that KS depends only on the temperature
!----------------------------------------------------------------------
               t1	= (t - t0)/(t*t0)
               t2	= t0/t - 1._r8 - log( t0/t )
               lnks     = lnks0 - 8792.3984_r8 * t1  - 16.8439_r8 * t2
               ks	= exp( lnks )

	       converged = .false.
!----------------------------------------------------------------------
!	Setting up initial guesses for the equations above.  Note that
!	for the initial choices the mass and the charge must be conserved.
!----------------------------------------------------------------------
               delx(0)  = .5_r8
               z        = .5_r8
               delz(0)  = .5_r8
               fratio(0) = 0._r8
               reduction_factor = .1_r8
               interval_set = .false.
Iter_loop :    do iter = 1,itermax
!----------------------------------------------------------------------
!	Cwater is the water equation as described in Tabazadeh et
!	al. and Jacobson et al.
!----------------------------------------------------------------------
                  cno3new   = chno3 * delx(iter-1)
                  cphno3new = chno3 * (1._r8 - delx(iter-1))
		  if( puren > 0._r8 ) then
		     t1 = (1000._r8*cno3new)/(puren*18._r8)
		  else
		     t1 = 0._r8
		  end if
	          h2o_cond(i,k) = t1 + wrk_h2so4
		  if( h2o_cond(i,k) > 0._r8 ) then
                     wrk      = 1.e3_r8 / (18._r8 * h2o_cond(i,k))
                     molhno3  = cno3new * wrk
                     molh2so4 = ch2so4(i,k) * wrk
		  else
                     molhno3  = 0._r8
                     molh2so4 = 0._r8
		  end if
                  stren	= molhno3 + 3._r8 * molh2so4
!----------------------------------------------------------------------
!	(1) Calculate the activity of H2SO4 at a given STREN
!----------------------------------------------------------------------
                  sm	= stren/3._r8
                  acts = d(13) + sm*(d(14) + sm*(d(15) + sm*d(16)))
!----------------------------------------------------------------------
!	(2) Calculate the activity for HNO3 at a given STREN
!----------------------------------------------------------------------
                  nm	= stren
                  actn 	= d(17) + nm*(d(18) + nm*(d(19) + nm*(d(20) + nm*(d(21) + nm*d(22)))))
!----------------------------------------------------------------------
!	(3) Calculate the mixed activity coefficient for HNO3 at STREN
!	    as described by Tabazadeh et al.
!----------------------------------------------------------------------
                  f1	 = 2._r8 * (molh2so4 + molhno3) * actn
                  f2	 = 2.25_r8 * molh2so4 * acts
                  mixyln = (f1 + f2) / (2._r8 * stren)
                  ymix	 = exp( mixyln )
                  hplus	 = 2._r8 * molh2so4 + molhno3
                  num = ymix**2 * hplus * molhno3
                  den = 1000._r8 * cphno3new * .0820578_r8 * t * ks
		  if( chno3 == 0._r8 ) then
	             converged = .true.
		     exit Iter_loop
		  end if
!----------------------------------------------------------------------
!       the denominator
!       Calculate the ratio F, check convergence
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!       Calculate the ratio F and reset the deltaX (see Jacobson et al.)
!----------------------------------------------------------------------
                  f = num / den
                  fratio(iter) = abs( f ) - 1._r8
		  con_val      = abs( f - 1._r8 )
		  if( con_val <= con_lim ) then
		     converged = .true.
		     exit Iter_loop
                  end if
!----------------------------------------------------------------------
!       non-convergence; setup next iterate
!----------------------------------------------------------------------
                  if( interval_set ) then
                     z = reduction_factor * z
                     delz(iter) = z
		     if( f > 1._r8 ) then
                        deltax = -z
		     else
                        deltax = z
		     end if
                     delx(iter) = delx(iter-1) + deltax
                  else
                     if( iter == 1 ) then
                        if( fratio(iter) >= 1._r8 ) then
                           positive = .false.
                        else
                           positive = .true.
                        end if
                     end if
                     if( fratio(iter)*fratio(iter-1) < 0._r8 ) then
                        interval_set = .true.
                        reduction_factor = .5_r8
                        delx(iter) = .5_r8*(delx(iter-1) + delx(iter-2))
                        z = .5_r8*abs( delx(iter-1) - delx(iter-2) )
                     else
                        if( .not. positive ) then
                           delx(iter) = reduction_factor * delx(iter-1)
                        else
                           delx(iter) = reduction_factor + delx(iter-1)
                           if( delx(iter) > upper_delx ) then
                              delx(iter) = .5_r8
                              interval_set = .true.
                              reduction_factor = .5_r8
                           end if
                        end if
                     end if
		  end if
               end do Iter_loop

               wtotal   = molhno3 * 63._r8 + molh2so4 * 98._r8 + 1000._r8
               wts(i,k) = (molh2so4 * 9800._r8) / wtotal
	       if( cno3new /= 0._r8 .or. cphno3new /= 0._r8 ) then
                  ratio	= max( 0._r8,min( 1._r8,cno3new/(cphno3new + cno3new) ) )
                  hno3_cond(i,k) = ratio*hno3_avail(i,k)
	       else
                  hno3_cond(i,k) = 0._r8
	       end if
	       if( .not. converged ) then
		  write(iulog,*) 'equil: Failed to converge @ is_chem,flag,lchnk,i,k,f = ',is_chem,flag,lchnk,i,k,f
!	          write(iulog,'(''equil: temper = '',z16)') temper(i,k)
		  write(iulog,*) '       h2o_avail,hno3_avail,p,t = ',h2o_avail(i,k),hno3_avail(i,k),press(i,k),temper(i,k)
		  write(iulog,*) '       molhno3,molh2so4,h2o_cond,hno3_cond = ',molhno3,molh2so4,h2o_cond(i,k),hno3_cond(i,k)
                  if( con_val > .05_r8 ) then
                     write(iulog,*) ' '
                     write(iulog,*) 'equil; diagnostics at lchnk, flag, i, k, iter = ',lchnk,flag,i,k,iter
                     write(iulog,*) 'equil; fratio'
                     write(iulog,'(5(1pg15.7))') fratio(1:iter)
                     write(iulog,*) ' '
                     write(iulog,*) 'equil; delx'
                     write(iulog,'(5(1pg15.7))') delx(0:iter-1)
                     write(iulog,*) ' '
                     write(iulog,*) 'equil; delz'
                     write(iulog,'(5(1pg15.7))') delz(0:iter-1)
                     write(iulog,*) ' '
	  	     call endrun('mo_sad.equil: ERROR 2 -- did not converge')
                  end if
	       end if
            end if
         end do Column_loop
      end do Level_loop

      end subroutine equil

!======================================================================
!
!
! ROUTINE
!   DENSITY
!
!   Date...
!     7 October 1999
!
!   Programmed by...
!     A. Tabazadeh
!
! DESCRIPTION
!     Calculates the density (g cm-3) of a binary sulfate solution.
!
! ARGUMENTS
!   INPUT
!      T           Temperature
!      w           Weight fraction
!
!   OUTPUT
!        den       Density of the Binary Solution (g cm-3)
!
!======================================================================
       
      function density( temp, w )

      implicit none

!----------------------------------------------------------------------
!	... Dummy arguments
!----------------------------------------------------------------------
      real(r8), intent(in) :: temp, w

!----------------------------------------------------------------------
!	... Function declarations
!----------------------------------------------------------------------
      real(r8) :: density

!----------------------------------------------------------------------
!	... Local variables
!----------------------------------------------------------------------
      real(r8), parameter :: a9 = -268.2616e4_r8, a10 = 576.4288e3_r8

      real(r8) :: a0, a1, a2, a3, a4, a5, a6, a7 ,a8
      real(r8) :: c1, c2, c3, c4

!----------------------------------------------------------------------
!	... Temperature variables
!----------------------------------------------------------------------
      c1 = temp - 273.15_r8
      c2 = c1**2
      c3 = c1*c2
      c4 = c1*c3
!----------------------------------------------------------------------
!	Polynomial Coefficients
!----------------------------------------------------------------------
      a0 = 999.8426_r8 + 334.5402e-4_r8*c1 - 569.1304e-5_r8*c2
      a1 = 547.2659_r8 - 530.0445e-2_r8*c1 + 118.7671e-4_r8*c2 + 599.0008e-6_r8*c3
      a2 = 526.295e1_r8 + 372.0445e-1_r8*c1 + 120.1909e-3_r8*c2 - 414.8594e-5_r8*c3 + 119.7973e-7_r8*c4
      a3 = -621.3958e2_r8 - 287.7670_r8*c1 - 406.4638e-3_r8*c2 + 111.9488e-4_r8*c3 + 360.7768e-7_r8*c4
      a4 = 409.0293e3_r8 + 127.0854e1_r8*c1 + 326.9710e-3_r8*c2 - 137.7435e-4_r8*c3 - 263.3585e-7_r8*c4
      a5 = -159.6989e4_r8 - 306.2836e1_r8*c1 + 136.6499e-3_r8*c2 + 637.3031e-5_r8*c3
      a6 = 385.7411e4_r8 + 408.3717e1_r8*c1 - 192.7785e-3_r8*c2
      a7 = -580.8064e4_r8 - 284.4401e1_r8*c1
      a8 = 530.1976e4_r8 + 809.1053_r8*c1
!----------------------------------------------------------------------
!	... Summation
!----------------------------------------------------------------------
      density = .001_r8*(a0 + w*(a1 + w*(a2 + w*(a3 + w*(a4 + w*(a5 + w*(a6 + w*(a7 + w*(a8 + w*(a9 + w*a10))))))))))

      end function density

      end module mo_sad
