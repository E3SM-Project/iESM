!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ecosys_mod

!BOP
! !MODULE: ecosys_mod
!
! !DESCRIPTION:
!
!  Multispecies ecosystem based on Doney et al. 1996, Moore et al., 2002
!  Based on POP Global NCAR Nitrogen Ecosystem Model
!  version 0.0 (June 15th, 1998) from S.C. Doney.
!  Based on Doney et al., 1996 model.
!  Climate and Global Dynamics, NCAR
!  (doney@whoi.edu)
!
!  Version 1.0
!  Multispecies, multiple limiting nutrient version of ecosystem
!  based on mixed layer model of Moore et al.(2002).  Implemented here with
!  fixed elemental ratios and including only the diatoms and small
!  phytoplankton, with a parameterization of calcification,
!  by Keith Lindsay and Keith Moore, Fall 2001 - Spring 2002.
!  Calcification parameterization based on Moore et al. 2002.
!
!  Version 2.0, January 2003
!    Adds diazotrophs as a phytoplankton group, (based on Moore et al., 2002a)
!    Allows for variable fe/C for all phytoplankton groups
!     Allows for variable si/C for the diatoms
!     Adds explicit tracers for DON, DOP, DOFe
!     variable remin length scale for detrital soft POM and bSi f(temperature)
!     Extensive modifications to iron scavenging parameterization
!     Addition of a sedimentary dissolved iron source,
!        (implemented in ballast code as excess remin in bottom cell)
!        coded by J.K. Moore, (jkmoore@uci.edu)
!
!   Version 2.01. March 2003
!     corrected O2 bug
!     corrected grazing parameter z_grz bug at depth
!     dust dissolution at depth releases iron,
!     increased length scale for dust diss., increased hard fraction dust
!     no deep ocean reduction in scavenging rates,
!     increase bSi OC/ballast ratio 0.3 -> 0.35,
!     corrected bug in diazotroph photoadaptation, and diat and sp adapatation
!
!   Version 2.02.
!     corrected bug in Fe_scavenge (units for dust), May 2003
!     changed C/N/P ratios to 117/16/1 (Anderson & Sarmiento, 1994)
!
!   Version 2.03., July 2003
!     Remin of DOM no longer temperature dependent,
!     new iron scavenging parameterization added,
!     some dissolution of hard fraction of ballast materials added
!
!   Version 2.1, September 2003
!     modfied iron scavenging and dust dissolution at depth
!
!   Version 2.11, March 2004
!     fixed bug in iron scavenging code, replace dust and POC flux_in w/ flux_out
!
!   Version 2.12, April 2004 - Final version for GBC paper revision,
!     (Questions/comments, Keith Moore - jkmoore@uci.edu
!
!   References
!   Doney, S.C., Glover, D.M., Najjar, R.G., 1996. A new coupled, one-dimensional
!   biological-physical model for the upper ocean: applications to the JGOFS
!   Bermuda Time-Series Study (BATS) site. Deep-Sea Res. II, 43: 591-624.
!
!   Moore, JK, Doney, SC, Kleypas, JA, Glover, DM, Fung, IY, 2002. An intermediate
!   complexity marine ecosystem model for the global domain. Deep-Sea Res. II, 49:
!   403-462.
!
!   Moore, JK, Doney, SC, Glover, DM, Fung, IY, 2002. Iron cycling and nutrient
!   limitation patterns in surface waters of the world ocean. Deep-Sea Res. II,
!   49: 463-507.

! !REVISION HISTORY:

!  SVN:$Id:  $

!-----------------------------------------------------------------------
!  variables/subroutines/function used from other modules
!  The following are used extensively in this ecosys, so are used at
!  the module level. The use statements for variables that are only needed
!  locally are located at the module subprogram level.
!-----------------------------------------------------------------------

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_GridHorzMod
   use POP_FieldMod
   use POP_HaloMod

   use kinds_mod
   use constants
   use communicate
   use broadcast
   use global_reductions
   use blocks
   use domain_size
   use domain
   use exit_mod
   use prognostic
   use grid
   use io
   use io_types
   use io_tools
   use tavg
   use timers
   use passive_tracer_tools
   use named_field_mod
   use forcing_tools
   use time_management
   use ecosys_parms
   use registry
   use named_field_mod
   use co2calc
#ifdef CCSMCOUPLED
   use POP_MCT_vars_mod
   use shr_strdata_mod
#endif

! !INPUT PARAMETERS:
!-----------------------------------------------------------------------
!  include ecosystem parameters
!  all variables from this modules have a parm_ prefix
!-----------------------------------------------------------------------

   implicit none
   save
   private

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   public :: &
      ecosys_tracer_cnt,            &
      ecosys_init,                  &
      ecosys_tracer_ref_val,        &
      ecosys_set_sflux,             &
      ecosys_tavg_forcing,          &
      ecosys_set_interior,          &
      ecosys_write_restart,         &
      ecosys_qsw_distrb_const

!-----------------------------------------------------------------------
!  module variables required by forcing_passive_tracer
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      ecosys_tracer_cnt = 24

!-----------------------------------------------------------------------
!  flags controlling which portion of code are executed
!  usefull for debugging
!-----------------------------------------------------------------------

  logical (log_kind) :: &
     lsource_sink, &
     lflux_gas_o2, &
     lflux_gas_co2,&
     locmip_k1_k2_bug_fix

  logical (log_kind), dimension(:,:,:), allocatable :: &
     LAND_MASK

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      po4_ind          =  1,  & ! dissolved inorganic phosphate
      no3_ind          =  2,  & ! dissolved inorganic nitrate
      sio3_ind         =  3,  & ! dissolved inorganic silicate
      nh4_ind          =  4,  & ! dissolved ammonia
      fe_ind           =  5,  & ! dissolved inorganic iron
      o2_ind           =  6,  & ! dissolved oxygen
      dic_ind          =  7,  & ! dissolved inorganic carbon
      alk_ind          =  8,  & ! alkalinity
      doc_ind          =  9,  & ! dissolved organic carbon
      spC_ind          = 10,  & ! small phytoplankton carbon
      spChl_ind        = 11,  & ! small phytoplankton chlorophyll
      spCaCO3_ind      = 12,  & ! small phytoplankton caco3
      diatC_ind        = 13,  & ! diatom carbon
      diatChl_ind      = 14,  & ! diatom chlorophyll
      zooC_ind         = 15,  & ! zooplankton carbon
      spFe_ind         = 16,  & ! small phytoplankton iron
      diatSi_ind       = 17,  & ! diatom silicon
      diatFe_ind       = 18,  & ! diatom iron
      diazC_ind        = 19,  & ! diazotroph carbon
      diazChl_ind      = 20,  & ! diazotroph Chlorophyll
      diazFe_ind       = 21,  & ! diazotroph iron
      don_ind          = 22,  & ! dissolved organic nitrogen
      dofe_ind         = 23,  & ! dissolved organic iron
      dop_ind          = 24     ! dissolved organic phosphorus

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(ecosys_tracer_cnt) :: &
      ind_name_table

!   type(ind_name_pair), dimension(ecosys_tracer_cnt) :: &
!      ind_name_table = (/ &
!      ind_name_pair(po4_ind,         'PO4'), &
!      ind_name_pair(no3_ind,         'NO3'), &
!      ind_name_pair(sio3_ind,        'SiO3'), &
!      ind_name_pair(nh4_ind,         'NH4'), &
!      ind_name_pair(fe_ind,          'Fe'), &
!      ind_name_pair(o2_ind,          'O2'), &
!      ind_name_pair(dic_ind,         'DIC'), &
!      ind_name_pair(alk_ind,         'ALK'), &
!      ind_name_pair(doc_ind,         'DOC'), &
!      ind_name_pair(spC_ind,         'spC'), &
!      ind_name_pair(spChl_ind,       'spChl'), &
!      ind_name_pair(spCaCO3_ind,     'spCaCO3'), &
!      ind_name_pair(diatC_ind,       'diatC'), &
!      ind_name_pair(diatChl_ind,     'diatChl'), &
!      ind_name_pair(zooC_ind,        'zooC'), &
!      ind_name_pair(spFe_ind,        'spFe'), &
!      ind_name_pair(diatSi_ind,      'diatSi'), &
!      ind_name_pair(diatFe_ind,      'diatFe'), &
!      ind_name_pair(diazC_ind,       'diazC'), &
!      ind_name_pair(diazChl_ind,     'diazChl'), &
!      ind_name_pair(diazFe_ind,      'diazFe'), &
!      ind_name_pair(don_ind,         'DON'), &
!      ind_name_pair(dofe_ind,        'DOFe'), &
!      ind_name_pair(dop_ind,         'DOP') /)
!
!-----------------------------------------------------------------------
!  options for forcing of gas fluxes
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      gas_flux_forcing_iopt_drv   = 1,   &
      gas_flux_forcing_iopt_file  = 2,   &
      atm_co2_iopt_const          = 1,   &
      atm_co2_iopt_drv_prog       = 2,   &
      atm_co2_iopt_drv_diag       = 3

   integer (int_kind) :: &
      gas_flux_forcing_iopt,             &
      atm_co2_iopt

   real (r8)       :: &
      atm_co2_const  ! value of atmospheric co2 (ppm, dry-air, 1 atm)

   character(char_len) :: &
      gas_flux_forcing_file    ! file containing gas flux forcing fields

!-----------------------------------------------------------------------

   real (r8) :: &
      rest_time_inv_surf,  & ! inverse restoring timescale at surface
      rest_time_inv_deep,  & ! inverse restoring timescale at depth
      rest_z0,             & ! shallow end of transition regime
      rest_z1                ! deep end of transition regime

   type(tracer_read) :: &
      po4_rest,            & ! restoring data for PO4
      no3_rest,            & ! restoring data for NO3
      sio3_rest,           & ! restoring data for SiO3
      gas_flux_fice,       & ! ice fraction for gas fluxes
      gas_flux_ws,         & ! wind speed for gas fluxes
      gas_flux_ap,         & ! atmospheric pressure for gas fluxes
      fesedflux_input        ! namelist input for iron_flux

!-----------------------------------------------------------------------
!  module variables related to ph computations
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:), allocatable, target :: &
      PH_PREV,        & ! computed ph from previous time step
      IRON_PATCH_FLUX   ! localized iron patch flux

   real (r8), dimension(:,:,:), allocatable :: &
      dust_FLUX_IN      ! dust flux not stored in STF since dust is not prognostic

   real (r8), dimension(:,:,:,:), allocatable :: &
      PH_PREV_3D        ! computed pH_3D from previous time step

!-----------------------------------------------------------------------
!  restoring climatologies for nutrients
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      lrest_po4,  & ! restoring on po4 ?
      lrest_no3,  & ! restoring on no3 ?
      lrest_sio3    ! restoring on sio3 ?

   real (r8), dimension(km) :: &
      nutr_rest_time_inv ! inverse restoring time scale for nutrients (1/secs)

   real (r8), dimension(:,:,:,:), allocatable, target :: &
      PO4_CLIM, NO3_CLIM, SiO3_CLIM

   real (r8), dimension(:,:,:,:), allocatable, target :: &
      FESEDFLUX

   character(char_len) :: &
      nutr_rest_file               ! file containing nutrient fields

!maltrud variable restoring
   logical (log_kind) :: &
      lnutr_variable_restore       ! geographically varying nutrient restoring

   character(char_len) :: &
      nutr_variable_rest_file,   & ! file containing variable restoring info
      nutr_variable_rest_file_fmt  ! format of file containing variable restoring info

   real (r8), dimension(:,:,:), allocatable, target :: &
      NUTR_RESTORE_RTAU            ! inverse restoring timescale for variable
                                   ! interior restoring

   integer (int_kind), dimension(:,:,:), allocatable :: &
      NUTR_RESTORE_MAX_LEVEL       ! maximum level for applying variable
                                   ! interior restoring

   real (r8), dimension(:,:,:,:), allocatable :: &
      INTERP_WORK                  ! temp array for interpolate_forcing output

   type(forcing_monthly_every_ts) :: &
      dust_flux,                 & ! surface dust flux
      iron_flux,                 & ! iron component of surface dust flux
      fice_file,                 & ! ice fraction, if read from file
      xkw_file,                  & ! a * wind-speed ** 2, if read from file
      ap_file                      ! atmoshperic pressure, if read from file

   character(char_len) :: &
      ndep_data_type               ! type of ndep forcing

   type(forcing_monthly_every_ts) :: &
      nox_flux_monthly,          & ! surface NOx species flux, added to nitrate pool
      nhy_flux_monthly             ! surface NHy species flux, added to ammonium pool

   integer (int_kind) :: &
      ndep_shr_stream_year_first, & ! first year in stream to use
      ndep_shr_stream_year_last,  & ! last year in stream to use
      ndep_shr_stream_year_align    ! align ndep_shr_stream_year_first with this model year

   integer (int_kind), parameter :: &
      ndep_shr_stream_var_cnt = 2, & ! number of variables in ndep shr_stream
      ndep_shr_stream_no_ind  = 1, & ! index for NO forcing
      ndep_shr_stream_nh_ind  = 2    ! index for NH forcing

   character(char_len) :: &
      ndep_shr_stream_file          ! file containing domain and input data

   real (r8) :: &
      ndep_shr_stream_scale_factor  ! unit conversion factor

#ifdef CCSMCOUPLED
   type(shr_strdata_type) :: ndep_sdat ! input data stream for ndep
#endif

!-----------------------------------------------------------------------
!  derived type for implicit handling of sinking particulate matter
!-----------------------------------------------------------------------

   type sinking_particle
      real (r8) :: &
         diss,        & ! dissolution length for soft subclass
         gamma,       & ! fraction of production -> hard subclass
         mass,        & ! mass of 1e9 base units in g
         rho            ! QA mass ratio of POC to this particle class

      real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
         sflux_in,    & ! incoming flux of soft subclass (base units/cm^2/sec)
         hflux_in,    & ! incoming flux of hard subclass (base units/cm^2/sec)
         prod,        & ! production term (base units/cm^3/sec)
         sflux_out,   & ! outgoing flux of soft subclass (base units/cm^2/sec)
         hflux_out,   & ! outgoing flux of hard subclass (base units/cm^2/sec)
         remin          ! remineralization term (base units/cm^3/sec)
   end type sinking_particle

!-----------------------------------------------------------------------
!  define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_ECOSYS_IFRAC, &! tavg id for ice fraction
      tavg_ECOSYS_IFRAC_2,&! tavg id for duplicate of ice fraction
      tavg_ECOSYS_XKW,   &! tavg id for xkw
      tavg_ECOSYS_XKW_2, &! tavg id for duplicate of xkw
      tavg_ECOSYS_ATM_PRESS, &! tavg id for atmospheric pressure
      tavg_PV_O2,        &! tavg id for o2 piston velocity
      tavg_SCHMIDT_O2,   &! tavg id for O2 schmidt number
      tavg_O2SAT,        &! tavg id for O2 saturation
      tavg_CO2STAR,      &! tavg id for co2star
      tavg_DCO2STAR,     &! tavg id for dco2star
      tavg_pCO2SURF,     &! tavg id for surface pco2
      tavg_DpCO2,        &! tavg id for delta pco2
      tavg_DpCO2_2,      &! tavg id for duplicate of delta pco2
      tavg_PV_CO2,       &! tavg id for co2 piston velocity
      tavg_SCHMIDT_CO2,  &! tavg id for co2 schmidt number
      tavg_DIC_GAS_FLUX, &! tavg id for dic flux
      tavg_DIC_GAS_FLUX_2,&! tavg id for duplicate of dic flux
      tavg_O2_GAS_FLUX_2,&! tavg id for duplicate of O2 flux
      tavg_PH,           &! tavg id for surface pH
      tavg_ATM_CO2,      &! tavg id for atmospheric CO2
      tavg_IRON_FLUX,    &! tavg id for dust flux
      tavg_DUST_FLUX,    &! tavg id for dust flux
      tavg_NOx_FLUX,     &! tavg id for nox flux
      tavg_NHy_FLUX       ! tavg id for nhy flux

!-----------------------------------------------------------------------
!  define tavg id for nonstandard 2d fields
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_O2_ZMIN,      &! tavg id for vertical minimum of O2
      tavg_O2_ZMIN_DEPTH  ! tavg id for depth of vertical minimum of O2

!-----------------------------------------------------------------------
!  define tavg id for nonstandard 3d fields
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_O2_PRODUCTION,&! tavg id for o2 production
      tavg_O2_CONSUMPTION,&! tavg id for o2 consumption
      tavg_AOU,          &! tavg id for AOU
      tavg_PO4_RESTORE,  &! tavg id for po4 restoring
      tavg_NO3_RESTORE,  &! tavg id for no3 restoring
      tavg_SiO3_RESTORE, &! tavg id for sio3 restoring
      tavg_PAR_avg,      &! tavg id for available radiation avg over mixed layer
      tavg_POC_FLUX_IN,  &! tavg id for poc flux into cell
      tavg_POC_PROD,     &! tavg id for poc production
      tavg_POC_REMIN,    &! tavg id for poc remineralization
      tavg_POC_ACCUM,    &! tavg id for poc accumulation
      tavg_CaCO3_FLUX_IN,&! tavg id for caco3 flux into cell
      tavg_CaCO3_PROD,   &! tavg id for caco3 production
      tavg_CaCO3_REMIN,  &! tavg id for caco3 remineralization
      tavg_SiO2_FLUX_IN, &! tavg id for sio2 flux into cell
      tavg_SiO2_PROD,    &! tavg id for sio2 production
      tavg_SiO2_REMIN,   &! tavg id for sio2 remineralization
      tavg_dust_FLUX_IN, &! tavg id for dust flux into cell
      tavg_dust_REMIN,   &! tavg id for dust remineralization
      tavg_P_iron_FLUX_IN, &! tavg id for p_iron flux into cell
      tavg_P_iron_PROD,    &! tavg id for p_iron production
      tavg_P_iron_REMIN,   &! tavg id for p_iron remineralization
      tavg_graze_sp,     &! tavg id for small phyto grazing
      tavg_graze_diat,   &! tavg id for diatom grazing
      tavg_graze_diaz,   &! tavg id for diaz grazing
      tavg_graze_TOT,    &! tavg id for total grazing
      tavg_sp_loss,      &! tavg id for small phyto loss
      tavg_diat_loss,    &! tavg id for diatom loss
      tavg_diaz_loss,    &! tavg id for diaz loss
      tavg_zoo_loss,     &! tavg id for zooplankton loss
      tavg_sp_agg,       &! tavg id for small phyto aggregate
      tavg_diat_agg       ! tavg id for diatom aggregate

!-----------------------------------------------------------------------
!  define tavg id for MORE nonstandard 3d fields
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_photoC_sp,            &! tavg id for small phyto C fixation
      tavg_CaCO3_form,           &! tavg id for CaCO3 formation
      tavg_photoC_diat,          &! tavg id for diatom C fixation
      tavg_photoC_diaz,          &! tavg id for diaz C fixation
      tavg_photoC_TOT,           &! tavg id for total C fixation
      tavg_photoC_sp_zint,       &! tavg id for small phyto C fixation vertical integral
      tavg_CaCO3_form_zint,      &! tavg id for CaCO3 formation vertical integral
      tavg_photoC_diat_zint,     &! tavg id for diatom C fixation vertical integral
      tavg_photoC_diaz_zint,     &! tavg id for diaz C fixation vertical integral
      tavg_photoC_TOT_zint,      &! tavg id for total C fixation vertical integral
      tavg_photoC_NO3_sp,        &! tavg id for small phyto C fixation from NO3
      tavg_photoC_NO3_sp_zint,   &! tavg id for small phyto C fixation from NO3 vertical integral
      tavg_photoC_NO3_diat,      &! tavg id for diatom C fixation from NO3
      tavg_photoC_NO3_diat_zint, &! tavg id for diatom C fixation from NO3 vertical integral
      tavg_photoC_NO3_diaz,      &! tavg id for diaz C fixation from NO3
      tavg_photoC_NO3_diaz_zint, &! tavg id for diaz C fixation from NO3 vertical integral
      tavg_photoC_NO3_TOT,       &! tavg id for total C fixation from NO3
      tavg_photoC_NO3_TOT_zint    ! tavg id for total C fixation from NO3 vertical integral

!-----------------------------------------------------------------------
!  define tavg id for MORE nonstandard 3d fields
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_DOC_prod,       &! tavg id for doc production
      tavg_DOC_remin,      &! tavg id for doc remineralization
      tavg_DON_prod,       &! tavg id for don production
      tavg_DON_remin,      &! tavg id for don remineralization
      tavg_DOFe_prod,      &! tavg id for dofe production
      tavg_DOFe_remin,     &! tavg id for dofe remineralization
      tavg_DOP_prod,       &! tavg id for dop production
      tavg_DOP_remin,      &! tavg id for dop remineralization
      tavg_Fe_scavenge,    &! tavg id for iron scavenging
      tavg_Fe_scavenge_rate,   &! tavg id for iron scavenging rate
      tavg_sp_N_lim,       &! tavg id for small phyto N limitation
      tavg_sp_Fe_lim,      &! tavg id for small phyto Fe limitation
      tavg_sp_PO4_lim,     &! tavg id for small phyto PO4 limitation
      tavg_sp_light_lim,   &! tavg id for small phyto light limitation
      tavg_diat_N_lim,     &! tavg id for diatom N limitation
      tavg_diat_Fe_lim,    &! tavg id for diatom Fe limitation
      tavg_diat_PO4_lim,   &! tavg id for diatom PO4 limitation
      tavg_diat_SiO3_lim,  &! tavg id for diatom SiO3 limitation
      tavg_diat_light_lim, &! tavg id for diatom light limitation
      tavg_diaz_Fe_lim,    &! tavg id for diaz Fe limitation
      tavg_diaz_P_lim,     &! tavg id for diaz PO4 limitation
      tavg_diaz_light_lim, &! tavg id for diaz light limitation
      tavg_diaz_Nfix,      &! tavg id for diaz N fixation
      tavg_bSi_form,       &! tavg id for diatom Si uptake
      tavg_NITRIF,         &! tavg id for nitrification
      tavg_DENITRIF,       &! tavg id for denitrification
      tavg_photoFe_sp,     &! tavg id for small phyto Fe uptake
      tavg_photoFe_diat,   &! tavg id for diatom Fe uptake
      tavg_photoFe_diaz     ! tavg id for diaz Fe uptake

   integer (int_kind) :: &
      tavg_CO3,            &! tavg id for 3D carbonate ion
      tavg_HCO3,           &! tavg id for 3D bicarbonate ion
      tavg_H2CO3,          &! tavg id for 3D carbonic acid
      tavg_pH_3D,          &! tavg id for 3D pH
      tavg_co3_sat_calc,   &! tavg id for co3 concentration at calcite saturation
      tavg_zsatcalc,       &! tavg id for calcite saturation depth
      tavg_co3_sat_arag,   &! tavg id for co3 concentration at aragonite saturation
      tavg_zsatarag         ! tavg id for aragonite saturation depth

!-----------------------------------------------------------------------
!  define array for holding flux-related quantities that need to be time-averaged
!  this is necessary since the forcing routines are called before tavg flags
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable :: &
      ECO_SFLUX_TAVG

!-----------------------------------------------------------------------
!  average surface tracer value related variables
!  used as reference value for virtual flux computations
!-----------------------------------------------------------------------

   logical (log_kind), dimension(ecosys_tracer_cnt) :: &
      vflux_flag                ! which tracers get virtual fluxes applied

   integer (int_kind) :: &
      comp_surf_avg_flag        ! time flag id for computing average
                                ! surface tracer values

   real (r8), dimension(ecosys_tracer_cnt) :: &
      surf_avg                  ! average surface tracer values

   logical (log_kind) :: &
      ecosys_qsw_distrb_const

!-----------------------------------------------------------------------
!  iron patch fertilization
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      liron_patch               ! flag for iron patch fertilization

   character(char_len) :: &
      iron_patch_flux_filename  ! file containing name of iron patch file

   integer (int_kind) :: &
      iron_patch_month          !  integer month to add patch flux

!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ecosys_shr_strdata_advance_timer,        &
      ecosys_comp_CO3terms_timer,              &
      ecosys_interior_timer,                   &
      ecosys_sflux_timer

!-----------------------------------------------------------------------
!  named field indices
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      totChl_surf_nf_ind = 0,    & ! total chlorophyll in surface layer
      sflux_co2_nf_ind   = 0,    & ! air-sea co2 gas flux
      atm_co2_nf_ind     = 0       ! atmospheric co2

!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      PAR_out           ! photosynthetically available radiation (W/m^2)

!-----------------------------------------------------------------------

   real (r8), parameter :: &
      phlo_surf_init = 7.0_r8, & ! low bound for surface ph for no prev soln
      phhi_surf_init = 9.0_r8, & ! high bound for surface ph for no prev soln
      phlo_3d_init = 6.0_r8,   & ! low bound for subsurface ph for no prev soln
      phhi_3d_init = 9.0_r8,   & ! high bound for subsurface ph for no prev soln
      del_ph = 0.20_r8           ! delta-ph for prev soln

!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: ecosys_init
! !INTERFACE:

 subroutine ecosys_init(init_ts_file_fmt, read_restart_filename, &
                        tracer_d_module, TRACER_MODULE, tadvect_ctype, &
                        errorCode)

! !DESCRIPTION:
!  Initialize ecosys tracer module. This involves setting metadata, reading
!  the module namelist, setting initial conditions, setting up forcing,
!  and defining additional tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

  !type (tracer_field), dimension(ecosys_tracer_cnt), intent(inout) :: &
   type (tracer_field), dimension(:)                , intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

  !real (r8), dimension(nx_block,ny_block,km,ecosys_tracer_cnt,3,max_blocks_clinic), &
   real (r8), dimension(:,:,:,:,:,:), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

  !character (char_len), dimension(ecosys_tracer_cnt), intent(out) :: &
   character (char_len), dimension(:), intent(out) :: &
      tadvect_ctype     ! advection method for ecosys tracers

   integer (POP_i4), intent(out) :: &
      errorCode

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'ecosys_mod:ecosys_init'

   character(char_len) :: &
      init_ecosys_option,        & ! option for initialization of bgc
      init_ecosys_init_file,     & ! filename for option 'file'
      init_ecosys_init_file_fmt, & ! file format for option 'file'
      comp_surf_avg_freq_opt,    & ! choice for freq of comp_surf_avg
      gas_flux_forcing_opt,      & ! option for forcing gas fluxes
      atm_co2_opt,               & ! option for atmospheric co2 concentration
      ecosys_tadvect_ctype         ! advection method for ecosys tracers

   type(tracer_read), dimension(ecosys_tracer_cnt) :: &
      tracer_init_ext              ! namelist variable for initializing tracers

   type(tracer_read) :: &
      dust_flux_input,           & ! namelist input for dust_flux
      iron_flux_input,           & ! namelist input for iron_flux
      nox_flux_monthly_input,    & ! namelist input for nox_flux_monthly
      nhy_flux_monthly_input       ! namelist input for nhy_flux_monthly

   logical (log_kind) :: &
      default,                   & ! arg to init_time_flag
      lnml_found,                & ! Was ecosys_nml found ?
      lmarginal_seas               ! Is ecosystem active in marginal seas ?

   integer (int_kind) :: &
      n,                         & ! index for looping over tracers
      k,                         & ! index for looping over depth levels
      l,                         & ! index for looping over time levels
      ind,                       & ! tracer index for tracer name from namelist
      iblock,                    & ! index for looping over blocks
      nml_error                    ! namelist i/o error flag

   integer (int_kind) :: &
      freq_opt, freq,            & ! args for init_time_flag
      comp_surf_avg_freq_iopt,   & ! choice for freq of comp_surf_avg
      comp_surf_avg_freq           ! choice for freq of comp_surf_avg

   logical (log_kind) :: &
      use_nml_surf_vals            ! do namelist surf values override values from restart file

   logical (log_kind) :: &
      lecovars_full_depth_tavg     ! should ecosystem vars be written full depth

!-----------------------------------------------------------------------
!  values to be used when comp_surf_avg_freq_opt==never
!-----------------------------------------------------------------------

   real (r8) :: &
      surf_avg_dic_const, surf_avg_alk_const

   namelist /ecosys_nml/ &
      init_ecosys_option, init_ecosys_init_file, tracer_init_ext, &
      init_ecosys_init_file_fmt, &
      dust_flux_input, iron_flux_input, fesedflux_input, &
      ndep_data_type, nox_flux_monthly_input, nhy_flux_monthly_input, &
      ndep_shr_stream_year_first, ndep_shr_stream_year_last, &
      ndep_shr_stream_year_align, ndep_shr_stream_file, &
      ndep_shr_stream_scale_factor, &
      gas_flux_forcing_opt, gas_flux_forcing_file, &
      gas_flux_fice, gas_flux_ws, gas_flux_ap, &
      lrest_po4, lrest_no3, lrest_sio3, &
      rest_time_inv_surf, rest_time_inv_deep, rest_z0, rest_z1, &
      nutr_rest_file, po4_rest, no3_rest, sio3_rest, &
      comp_surf_avg_freq_opt, comp_surf_avg_freq,  &
      use_nml_surf_vals, surf_avg_dic_const, surf_avg_alk_const, &
      ecosys_qsw_distrb_const, lmarginal_seas, &
      lsource_sink, lflux_gas_o2, lflux_gas_co2, locmip_k1_k2_bug_fix, &
      lnutr_variable_restore, nutr_variable_rest_file,  &
      nutr_variable_rest_file_fmt,atm_co2_opt,atm_co2_const, &
      ecosys_tadvect_ctype, &
      liron_patch,iron_patch_flux_filename,iron_patch_month, &
      lecovars_full_depth_tavg

   character (char_len) :: &
      ecosys_restart_filename  ! modified file name for restart file

   real (r8), dimension (nx_block,ny_block) :: WORK

!-----------------------------------------------------------------------
!  initialize name table 
!-----------------------------------------------------------------------

   errorCode = POP_Success

   ind_name_table( 1) = ind_name_pair(po4_ind,     'PO4')
   ind_name_table( 2) = ind_name_pair(no3_ind,     'NO3')
   ind_name_table( 3) = ind_name_pair(sio3_ind,    'SiO3')
   ind_name_table( 4) = ind_name_pair(nh4_ind,     'NH4')
   ind_name_table( 5) = ind_name_pair(fe_ind,      'Fe')
   ind_name_table( 6) = ind_name_pair(o2_ind,      'O2')
   ind_name_table( 7) = ind_name_pair(dic_ind,     'DIC')
   ind_name_table( 8) = ind_name_pair(alk_ind,     'ALK')
   ind_name_table( 9) = ind_name_pair(doc_ind,     'DOC')
   ind_name_table(10) = ind_name_pair(spC_ind,     'spC')
   ind_name_table(11) = ind_name_pair(spChl_ind,   'spChl')
   ind_name_table(12) = ind_name_pair(spCaCO3_ind, 'spCaCO3')
   ind_name_table(13) = ind_name_pair(diatC_ind,   'diatC')
   ind_name_table(14) = ind_name_pair(diatChl_ind, 'diatChl')
   ind_name_table(15) = ind_name_pair(zooC_ind,    'zooC')
   ind_name_table(16) = ind_name_pair(spFe_ind,    'spFe')
   ind_name_table(17) = ind_name_pair(diatSi_ind,  'diatSi')
   ind_name_table(18) = ind_name_pair(diatFe_ind,  'diatFe')
   ind_name_table(19) = ind_name_pair(diazC_ind,   'diazC')
   ind_name_table(20) = ind_name_pair(diazChl_ind, 'diazChl')
   ind_name_table(21) = ind_name_pair(diazFe_ind,  'diazFe')
   ind_name_table(22) = ind_name_pair(don_ind,     'DON')
   ind_name_table(23) = ind_name_pair(dofe_ind,    'DOFe')
   ind_name_table(24) = ind_name_pair(dop_ind,     'DOP')

!-----------------------------------------------------------------------
!  initialize forcing_monthly_every_ts variables
!-----------------------------------------------------------------------

   call init_forcing_monthly_every_ts(dust_flux)
   call init_forcing_monthly_every_ts(iron_flux)
   call init_forcing_monthly_every_ts(fice_file)
   call init_forcing_monthly_every_ts(xkw_file)
   call init_forcing_monthly_every_ts(ap_file)
   call init_forcing_monthly_every_ts(nox_flux_monthly)
   call init_forcing_monthly_every_ts(nhy_flux_monthly)

!-----------------------------------------------------------------------
!  initialize ecosystem parameters
!-----------------------------------------------------------------------

   call ecosys_parms_init

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   do n = 1, ecosys_tracer_cnt
      tracer_d_module(n)%short_name = ind_name_table(n)%name
      if (any(n == (/ spChl_ind, diatChl_ind, diazChl_ind /))) then
         tracer_d_module(n)%units      = 'mg/m^3'
         tracer_d_module(n)%tend_units = 'mg/m^3/s'
         tracer_d_module(n)%flux_units = 'mg/m^3 cm/s'
      else if (n == alk_ind) then
         tracer_d_module(n)%units      = 'meq/m^3'
         tracer_d_module(n)%tend_units = 'meq/m^3/s'
         tracer_d_module(n)%flux_units = 'meq/m^3 cm/s'
      else
         tracer_d_module(n)%units      = 'mmol/m^3'
         tracer_d_module(n)%tend_units = 'mmol/m^3/s'
         tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
      endif
   end do
   tracer_d_module(po4_ind    )%long_name='Dissolved Inorganic Phosphate'
   tracer_d_module(no3_ind    )%long_name='Dissolved Inorganic Nitrate'
   tracer_d_module(sio3_ind   )%long_name='Dissolved Inorganic Silicate'
   tracer_d_module(nh4_ind    )%long_name='Dissolved Ammonia'
   tracer_d_module(fe_ind     )%long_name='Dissolved Inorganic Iron'
   tracer_d_module(o2_ind     )%long_name='Dissolved Oxygen'
   tracer_d_module(dic_ind    )%long_name='Dissolved Inorganic Carbon'
   tracer_d_module(alk_ind    )%long_name='Alkalinity'
   tracer_d_module(doc_ind    )%long_name='Dissolved Organic Carbon'
   tracer_d_module(spC_ind    )%long_name='Small Phytoplankton Carbon'
   tracer_d_module(spChl_ind  )%long_name='Small phytoplankton Chlorophyll'
   tracer_d_module(spCaCO3_ind)%long_name='Small Phytoplankton CaCO3'
   tracer_d_module(diatC_ind  )%long_name='Diatom Carbon'
   tracer_d_module(diatChl_ind)%long_name='Diatom Chlorophyll'
   tracer_d_module(zooC_ind   )%long_name='Zooplankton Carbon'
   tracer_d_module(spFe_ind   )%long_name='Small Phytoplankton Iron'
   tracer_d_module(diatSi_ind )%long_name='Diatom Silicon'
   tracer_d_module(diatFe_ind )%long_name='Diatom Iron'
   tracer_d_module(diazC_ind  )%long_name='Diazotroph Carbon'
   tracer_d_module(diazChl_ind)%long_name='Diazotroph Chlorophyll'
   tracer_d_module(diazFe_ind )%long_name='Diazotroph Iron'
   tracer_d_module(don_ind    )%long_name='Dissolved Organic Nitrogen'
   tracer_d_module(dofe_ind   )%long_name='Dissolved Organic Iron'
   tracer_d_module(dop_ind    )%long_name='Dissolved Organic Phosphorus'

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_ecosys_option = 'unknown'
   init_ecosys_init_file = 'unknown'
   init_ecosys_init_file_fmt = 'bin'

   gas_flux_forcing_opt  = 'drv'
   gas_flux_forcing_file = 'unknown'

   gas_flux_fice%filename     = 'unknown'
   gas_flux_fice%file_varname = 'FICE'
   gas_flux_fice%scale_factor = c1
   gas_flux_fice%default_val  = c0
   gas_flux_fice%file_fmt     = 'bin'

   gas_flux_ws%filename     = 'unknown'
   gas_flux_ws%file_varname = 'XKW'
   gas_flux_ws%scale_factor = c1
   gas_flux_ws%default_val  = c0
   gas_flux_ws%file_fmt     = 'bin'

   gas_flux_ap%filename     = 'unknown'
   gas_flux_ap%file_varname = 'P'
   gas_flux_ap%scale_factor = c1
   gas_flux_ap%default_val  = c0
   gas_flux_ap%file_fmt     = 'bin'

   lrest_po4      = .false.
   lrest_no3      = .false.
   lrest_sio3     = .false.
   nutr_rest_file = 'unknown'

   rest_time_inv_surf = c0
   rest_time_inv_deep = c0
   rest_z0            = c1000
   rest_z1            = c2 * c1000

   po4_rest%filename     = 'unknown'
   po4_rest%file_varname = tracer_d_module(po4_ind)%short_name
   po4_rest%scale_factor = c1
   po4_rest%default_val  = c0
   po4_rest%file_fmt     = 'bin'

   no3_rest%filename     = 'unknown'
   no3_rest%file_varname = tracer_d_module(no3_ind)%short_name
   no3_rest%scale_factor = c1
   no3_rest%default_val  = c0
   no3_rest%file_fmt     = 'bin'

   sio3_rest%filename     = 'unknown'
   sio3_rest%file_varname = tracer_d_module(sio3_ind)%short_name
   sio3_rest%scale_factor = c1
   sio3_rest%default_val  = c0
   sio3_rest%file_fmt     = 'bin'

!maltrud variable restoring
   lnutr_variable_restore      = .false.
   nutr_variable_rest_file     = 'unknown'
   nutr_variable_rest_file_fmt = 'bin'

   dust_flux_input%filename     = 'unknown'
   dust_flux_input%file_varname = 'dust_flux'
   dust_flux_input%scale_factor = c1
   dust_flux_input%default_val  = c0
   dust_flux_input%file_fmt     = 'bin'

   iron_flux_input%filename     = 'unknown'
   iron_flux_input%file_varname = 'iron_flux'
   iron_flux_input%scale_factor = c1
   iron_flux_input%default_val  = c0
   iron_flux_input%file_fmt     = 'bin'

   fesedflux_input%filename     = 'unknown'
   fesedflux_input%file_varname = 'FESEDFLUXIN'
   fesedflux_input%scale_factor = c1
   fesedflux_input%default_val  = c0
   fesedflux_input%file_fmt     = 'bin'

   ndep_data_type              = 'monthly-calendar'

   nox_flux_monthly_input%filename     = 'unknown'
   nox_flux_monthly_input%file_varname = 'nox_flux'
   nox_flux_monthly_input%scale_factor = c1
   nox_flux_monthly_input%default_val  = c0
   nox_flux_monthly_input%file_fmt     = 'bin'

   nhy_flux_monthly_input%filename     = 'unknown'
   nhy_flux_monthly_input%file_varname = 'nhy_flux'
   nhy_flux_monthly_input%scale_factor = c1
   nhy_flux_monthly_input%default_val  = c0
   nhy_flux_monthly_input%file_fmt     = 'bin'

   ndep_shr_stream_year_first = 1
   ndep_shr_stream_year_last  = 1
   ndep_shr_stream_year_align = 1
   ndep_shr_stream_file       = 'unknown'
   ndep_shr_stream_scale_factor = c1

   do n = 1,ecosys_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'
   end do

   lmarginal_seas        = .true.
   lsource_sink          = .true.
   lflux_gas_o2          = .true.
   lflux_gas_co2         = .true.
   locmip_k1_k2_bug_fix  = .true.

   comp_surf_avg_freq_opt        = 'never'
   comp_surf_avg_freq            = 1
   use_nml_surf_vals             = .false.
   surf_avg_dic_const            = 1944.0_r8
   surf_avg_alk_const            = 2225.0_r8

   ecosys_qsw_distrb_const  = .true.

   liron_patch              = .false.
   iron_patch_flux_filename = 'unknown_iron_patch_filename'
   iron_patch_month         = 1

   atm_co2_opt   = 'const'
   atm_co2_const = 280.0_r8

   ecosys_tadvect_ctype = 'base_model'

   lecovars_full_depth_tavg = .false.

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=ecosys_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'ecosys_nml not found')
      call exit_POP(sigAbort, 'ERROR : stopping in '/&
                           &/ subname)
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' ecosys:'
      write(stdout,blank_fmt)
      write(stdout,*) ' ecosys_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,ecosys_nml)
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_ecosys_option, master_task)
   call broadcast_scalar(init_ecosys_init_file, master_task)
   call broadcast_scalar(init_ecosys_init_file_fmt, master_task)

   call broadcast_scalar(gas_flux_forcing_opt, master_task)
   if (trim(gas_flux_forcing_opt) == 'drv') then
      gas_flux_forcing_iopt = gas_flux_forcing_iopt_drv
   else if (trim(gas_flux_forcing_opt) == 'file') then
      gas_flux_forcing_iopt = gas_flux_forcing_iopt_file
   else
      call document(subname, 'gas_flux_forcing_opt', gas_flux_forcing_opt)
      call exit_POP(sigAbort, 'unknown gas_flux_forcing_opt')
   endif

   call broadcast_scalar(gas_flux_forcing_file, master_task)

   call broadcast_scalar(gas_flux_fice%filename, master_task)
   call broadcast_scalar(gas_flux_fice%file_varname, master_task)
   call broadcast_scalar(gas_flux_fice%scale_factor, master_task)
   call broadcast_scalar(gas_flux_fice%default_val, master_task)
   call broadcast_scalar(gas_flux_fice%file_fmt, master_task)

   fice_file%input = gas_flux_fice

   call broadcast_scalar(gas_flux_ws%filename, master_task)
   call broadcast_scalar(gas_flux_ws%file_varname, master_task)
   call broadcast_scalar(gas_flux_ws%scale_factor, master_task)
   call broadcast_scalar(gas_flux_ws%default_val, master_task)
   call broadcast_scalar(gas_flux_ws%file_fmt, master_task)

   xkw_file%input = gas_flux_ws

   call broadcast_scalar(gas_flux_ap%filename, master_task)
   call broadcast_scalar(gas_flux_ap%file_varname, master_task)
   call broadcast_scalar(gas_flux_ap%scale_factor, master_task)
   call broadcast_scalar(gas_flux_ap%default_val, master_task)
   call broadcast_scalar(gas_flux_ap%file_fmt, master_task)

   ap_file%input = gas_flux_ap

   call broadcast_scalar(lrest_po4, master_task)
   call broadcast_scalar(lrest_no3, master_task)
   call broadcast_scalar(lrest_sio3, master_task)
   call broadcast_scalar(nutr_rest_file, master_task)

   call broadcast_scalar(rest_time_inv_surf, master_task)
   call broadcast_scalar(rest_time_inv_deep, master_task)
   call broadcast_scalar(rest_z0, master_task)
   call broadcast_scalar(rest_z1, master_task)

   call broadcast_scalar(po4_rest%filename, master_task)
   call broadcast_scalar(po4_rest%file_varname, master_task)
   call broadcast_scalar(po4_rest%scale_factor, master_task)
   call broadcast_scalar(po4_rest%default_val, master_task)
   call broadcast_scalar(po4_rest%file_fmt, master_task)

   call broadcast_scalar(no3_rest%filename, master_task)
   call broadcast_scalar(no3_rest%file_varname, master_task)
   call broadcast_scalar(no3_rest%scale_factor, master_task)
   call broadcast_scalar(no3_rest%default_val, master_task)
   call broadcast_scalar(no3_rest%file_fmt, master_task)

   call broadcast_scalar(sio3_rest%filename, master_task)
   call broadcast_scalar(sio3_rest%file_varname, master_task)
   call broadcast_scalar(sio3_rest%scale_factor, master_task)
   call broadcast_scalar(sio3_rest%default_val, master_task)
   call broadcast_scalar(sio3_rest%file_fmt, master_task)

!maltrud variable restoring
   call broadcast_scalar(lnutr_variable_restore, master_task)
   call broadcast_scalar(nutr_variable_rest_file, master_task)
   call broadcast_scalar(nutr_variable_rest_file_fmt, master_task)

   call broadcast_scalar(dust_flux_input%filename, master_task)
   call broadcast_scalar(dust_flux_input%file_varname, master_task)
   call broadcast_scalar(dust_flux_input%scale_factor, master_task)
   call broadcast_scalar(dust_flux_input%default_val, master_task)
   call broadcast_scalar(dust_flux_input%file_fmt, master_task)

   dust_flux%input = dust_flux_input

   call broadcast_scalar(iron_flux_input%filename, master_task)
   call broadcast_scalar(iron_flux_input%file_varname, master_task)
   call broadcast_scalar(iron_flux_input%scale_factor, master_task)
   call broadcast_scalar(iron_flux_input%default_val, master_task)
   call broadcast_scalar(iron_flux_input%file_fmt, master_task)

   iron_flux%input = iron_flux_input

   call broadcast_scalar(fesedflux_input%filename, master_task)
   call broadcast_scalar(fesedflux_input%file_varname, master_task)
   call broadcast_scalar(fesedflux_input%scale_factor, master_task)
   call broadcast_scalar(fesedflux_input%default_val, master_task)
   call broadcast_scalar(fesedflux_input%file_fmt, master_task)

   call broadcast_scalar(ndep_data_type, master_task)

   call broadcast_scalar(nox_flux_monthly_input%filename, master_task)
   call broadcast_scalar(nox_flux_monthly_input%file_varname, master_task)
   call broadcast_scalar(nox_flux_monthly_input%scale_factor, master_task)
   call broadcast_scalar(nox_flux_monthly_input%default_val, master_task)
   call broadcast_scalar(nox_flux_monthly_input%file_fmt, master_task)

   nox_flux_monthly%input = nox_flux_monthly_input

   call broadcast_scalar(nhy_flux_monthly_input%filename, master_task)
   call broadcast_scalar(nhy_flux_monthly_input%file_varname, master_task)
   call broadcast_scalar(nhy_flux_monthly_input%scale_factor, master_task)
   call broadcast_scalar(nhy_flux_monthly_input%default_val, master_task)
   call broadcast_scalar(nhy_flux_monthly_input%file_fmt, master_task)

   nhy_flux_monthly%input = nhy_flux_monthly_input

   call broadcast_scalar(ndep_shr_stream_year_first, master_task)
   call broadcast_scalar(ndep_shr_stream_year_last, master_task)
   call broadcast_scalar(ndep_shr_stream_year_align, master_task)
   call broadcast_scalar(ndep_shr_stream_file, master_task)
   call broadcast_scalar(ndep_shr_stream_scale_factor, master_task)

   do n = 1,ecosys_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
   end do

   call broadcast_scalar(comp_surf_avg_freq_opt, master_task)
   call broadcast_scalar(comp_surf_avg_freq, master_task)
   call broadcast_scalar(use_nml_surf_vals, master_task)
   call broadcast_scalar(surf_avg_dic_const, master_task)
   call broadcast_scalar(surf_avg_alk_const, master_task)

   call broadcast_scalar(ecosys_qsw_distrb_const, master_task)

   call broadcast_scalar(lmarginal_seas, master_task)
   call broadcast_scalar(lsource_sink, master_task)
   call broadcast_scalar(lflux_gas_o2, master_task)
   call broadcast_scalar(lflux_gas_co2, master_task)
   call broadcast_scalar(locmip_k1_k2_bug_fix, master_task)

   call broadcast_scalar(liron_patch, master_task)
   call broadcast_scalar(iron_patch_flux_filename, master_task)
   call broadcast_scalar(iron_patch_month, master_task)

   call broadcast_scalar(atm_co2_opt, master_task)
   call broadcast_scalar(atm_co2_const, master_task)

   call broadcast_scalar(ecosys_tadvect_ctype, master_task)
   tadvect_ctype = ecosys_tadvect_ctype

   call broadcast_scalar(lecovars_full_depth_tavg, master_task)

!-----------------------------------------------------------------------
!  set variables immediately dependent on namelist variables
!-----------------------------------------------------------------------

   select case (comp_surf_avg_freq_opt)
   case ('never')
      comp_surf_avg_freq_iopt = freq_opt_never
   case ('nyear')
      comp_surf_avg_freq_iopt = freq_opt_nyear
   case ('nmonth')
      comp_surf_avg_freq_iopt = freq_opt_nmonth
   case default
      call document(subname, 'comp_surf_avg_freq_opt', comp_surf_avg_freq_opt)
      call exit_POP(sigAbort, 'unknown comp_surf_avg_freq_opt')
   end select

  call init_time_flag('ecosys_comp_surf_avg', comp_surf_avg_flag, &
     default=.false., freq_opt=comp_surf_avg_freq_iopt,  &
     freq=comp_surf_avg_freq, owner='ecosys_init')

   select case (atm_co2_opt)
   case ('const')
      atm_co2_iopt = atm_co2_iopt_const
   case ('drv_prog')
      atm_co2_iopt = atm_co2_iopt_drv_prog
   case ('drv_diag')
      atm_co2_iopt = atm_co2_iopt_drv_diag
   case default
      call document(subname, 'atm_co2_opt', atm_co2_opt)
      call exit_POP(sigAbort, 'unknown atm_co2_opt')
   end select

!-----------------------------------------------------------------------
!  namelist consistency checking
!-----------------------------------------------------------------------

   if (use_nml_surf_vals .and. comp_surf_avg_freq_iopt /= freq_opt_never) then
      call document(subname, 'use_nml_surf_vals', use_nml_surf_vals)
      call document(subname, 'comp_surf_avg_freq_opt', comp_surf_avg_freq_opt)
      call exit_POP(sigAbort, 'use_nml_surf_vals can only be .true. if ' /&
                           &/ ' comp_surf_avg_freq_opt is never')
   endif

!-----------------------------------------------------------------------
!  initialize virtual flux flag array
!-----------------------------------------------------------------------

   vflux_flag = .false.
   vflux_flag(dic_ind) = .true.
   vflux_flag(alk_ind) = .true.

!-----------------------------------------------------------------------
!  allocate various ecosys allocatable module variables
!-----------------------------------------------------------------------

   allocate( PH_PREV(nx_block,ny_block,max_blocks_clinic) )
   allocate( dust_FLUX_IN(nx_block,ny_block,max_blocks_clinic) )
   allocate( PH_PREV_3D(nx_block,ny_block,km,max_blocks_clinic) )
   PH_PREV_3D = c0

!-----------------------------------------------------------------------
!  allocate and initialize LAND_MASK
!-----------------------------------------------------------------------

   allocate( LAND_MASK(nx_block,ny_block,nblocks_clinic) )

   if (lmarginal_seas) then
      LAND_MASK = REGION_MASK /= 0
   else
      LAND_MASK = REGION_MASK > 0
   endif

!-----------------------------------------------------------------------
!  initialize tracers
!-----------------------------------------------------------------------

   select case (init_ecosys_option)

   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid')

      ecosys_restart_filename = char_blank

      if (init_ecosys_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read ecosys from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         ecosys_restart_filename = read_restart_filename
         init_ecosys_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         ecosys_restart_filename = trim(init_ecosys_init_file)

      endif

      call rest_read_tracer_block(init_ecosys_init_file_fmt, &
                                  ecosys_restart_filename,   &
                                  tracer_d_module,           &
                                  TRACER_MODULE)

      call read_field(init_ecosys_init_file_fmt, &
                      ecosys_restart_filename,   &
                      'PH_SURF', PH_PREV)

      if (use_nml_surf_vals) then
         surf_avg = c0
         surf_avg(dic_ind) = surf_avg_dic_const
         surf_avg(alk_ind) = surf_avg_alk_const
      else
         call extract_surf_avg(init_ecosys_init_file_fmt, &
                               ecosys_restart_filename)
      endif

      call eval_time_flag(comp_surf_avg_flag) ! evaluates time_flag(comp_surf_avg_flag)%value via time_to_do

      if (check_time_flag(comp_surf_avg_flag)) &
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:), &
                            TRACER_MODULE(:,:,1,:,curtime,:))

   case ('file', 'ccsm_startup')
      call document(subname, 'ecosystem vars being read from separate files')

      call file_read_tracer_block(init_ecosys_init_file_fmt, &
                                  init_ecosys_init_file,     &
                                  tracer_d_module,           &
                                  ind_name_table,            &
                                  tracer_init_ext,           &
                                  TRACER_MODULE)

      if (n_topo_smooth > 0) then
         do n = 1, ecosys_tracer_cnt
            do k=1,km
               call fill_points(k,TRACER_MODULE(:,:,k,n,oldtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'ecosys_init: error in fill points for tracers(oldtime)')
                  return
               endif

               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'ecosys_init: error in fill points for tracers(newtime)')
                  return
               endif

            enddo
         enddo
      endif

      PH_PREV = c0

      if (use_nml_surf_vals) then
         surf_avg = c0
         surf_avg(dic_ind) = surf_avg_dic_const
         surf_avg(alk_ind) = surf_avg_alk_const
      else
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:), &
                            TRACER_MODULE(:,:,1,:,curtime,:))
      endif

   case default
      call document(subname, 'init_ecosys_option', init_ecosys_option)
      call exit_POP(sigAbort, 'unknown init_ecosys_option')

   end select

!-----------------------------------------------------------------------
!  register Chl field for short-wave absorption
!  apply land mask to tracers
!  set Chl field for short-wave absorption
!-----------------------------------------------------------------------

   call named_field_register('model_chlorophyll', totChl_surf_nf_ind)

   !$OMP PARALLEL DO PRIVATE(iblock,n,k,WORK)
   do iblock=1,nblocks_clinic
      do n = 1,ecosys_tracer_cnt
         do k = 1,km
            where (.not. LAND_MASK(:,:,iblock) .or. k > KMT(:,:,iblock))
               TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
               TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            end where
         end do
      end do

      WORK = max(c0,p5*(TRACER_MODULE(:,:,1,spChl_ind,oldtime,iblock) + &
                        TRACER_MODULE(:,:,1,spChl_ind,curtime,iblock))) + &
             max(c0,p5*(TRACER_MODULE(:,:,1,diatChl_ind,oldtime,iblock) + &
                        TRACER_MODULE(:,:,1,diatChl_ind,curtime,iblock))) + &
             max(c0,p5*(TRACER_MODULE(:,:,1,diazChl_ind,oldtime,iblock) + &
                        TRACER_MODULE(:,:,1,diazChl_ind,curtime,iblock)))
      call named_field_set(totChl_surf_nf_ind, iblock, WORK)
   enddo
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  timer init
!-----------------------------------------------------------------------

   call get_timer(ecosys_comp_CO3terms_timer, 'comp_CO3terms', &
                  nblocks_clinic, distrb_clinic%nprocs)
   call get_timer(ecosys_interior_timer, 'ECOSYS_INTERIOR', &
                  nblocks_clinic, distrb_clinic%nprocs)
   call get_timer(ecosys_sflux_timer, 'ECOSYS_SFLUX',1, &
                  distrb_clinic%nprocs)
   if (ndep_data_type == 'shr_stream') then
      call get_timer(ecosys_shr_strdata_advance_timer, &
                     'ecosys_shr_strdata_advance',1, distrb_clinic%nprocs)
   endif

!-----------------------------------------------------------------------
!  call other initialization subroutines
!-----------------------------------------------------------------------

   call ecosys_init_tavg
   call ecosys_init_sflux
   call ecosys_init_interior_restore

!-----------------------------------------------------------------------
!  set lfull_depth_tavg flag for short-lived ecosystem tracers
!-----------------------------------------------------------------------

   tracer_d_module(spC_ind    )%lfull_depth_tavg = lecovars_full_depth_tavg
   tracer_d_module(spChl_ind  )%lfull_depth_tavg = lecovars_full_depth_tavg
   tracer_d_module(spCaCO3_ind)%lfull_depth_tavg = lecovars_full_depth_tavg
   tracer_d_module(diatC_ind  )%lfull_depth_tavg = lecovars_full_depth_tavg
   tracer_d_module(diatChl_ind)%lfull_depth_tavg = lecovars_full_depth_tavg
   tracer_d_module(zooC_ind   )%lfull_depth_tavg = lecovars_full_depth_tavg
   tracer_d_module(spFe_ind   )%lfull_depth_tavg = lecovars_full_depth_tavg
   tracer_d_module(diatSi_ind )%lfull_depth_tavg = lecovars_full_depth_tavg
   tracer_d_module(diatFe_ind )%lfull_depth_tavg = lecovars_full_depth_tavg
   tracer_d_module(diazC_ind  )%lfull_depth_tavg = lecovars_full_depth_tavg
   tracer_d_module(diazChl_ind)%lfull_depth_tavg = lecovars_full_depth_tavg
   tracer_d_module(diazFe_ind )%lfull_depth_tavg = lecovars_full_depth_tavg

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_init

!***********************************************************************
!BOP
! !IROUTINE: extract_surf_avg
! !INTERFACE:

 subroutine extract_surf_avg(init_ecosys_init_file_fmt, &
                             ecosys_restart_filename)

! !DESCRIPTION:
!  Extract average surface values from restart file.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      init_ecosys_init_file_fmt, & ! file format (bin or nc)
      ecosys_restart_filename      ! file name for restart file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   type (datafile) ::&
      restart_file    ! io file descriptor

   integer (int_kind) :: &
      n               ! tracer index

   character (char_len) :: &
      short_name      ! tracer name temporaries

!-----------------------------------------------------------------------

   surf_avg = c0

   restart_file = construct_file(init_ecosys_init_file_fmt, &
                                 full_name=trim(ecosys_restart_filename), &
                                 record_length=rec_type_dbl, &
                                 recl_words=nx_global*ny_global)

   do n = 1, ecosys_tracer_cnt
      if (vflux_flag(n)) then
         short_name = 'surf_avg_' /&
                   &/ ind_name_table(n)%name
         call add_attrib_file(restart_file, trim(short_name), surf_avg(n))
      endif
   end do

   call data_set(restart_file, 'open_read')

   do n = 1, ecosys_tracer_cnt
      if (vflux_flag(n)) then
         short_name = 'surf_avg_' /&
                   &/ ind_name_table(n)%name
         call extract_attrib_file(restart_file, trim(short_name), surf_avg(n))
      endif
   end do

   call data_set (restart_file, 'close')

   call destroy_file (restart_file)

!-----------------------------------------------------------------------
!EOC

 end subroutine extract_surf_avg

!***********************************************************************
!BOP
! !IROUTINE: ecosys_init_tavg
! !INTERFACE:

 subroutine ecosys_init_tavg

! !DESCRIPTION:
!  call define_tavg_field for nonstandard tavg fields
!
! !REVISION HISTORY:
!  same as module
!

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      var_cnt             ! how many tavg variables are defined

!-----------------------------------------------------------------------
!  2D fields related to surface fluxes
!-----------------------------------------------------------------------

   var_cnt = 0

   call define_tavg_field(tavg_ECOSYS_IFRAC,'ECOSYS_IFRAC',2,          &
                          long_name='Ice Fraction for ecosys fluxes',  &
                          units='fraction', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ECOSYS_IFRAC_2,'ECOSYS_IFRAC_2',2,      &
                          long_name='Ice Fraction for ecosys fluxes',  &
                          units='fraction', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_ECOSYS_XKW,'ECOSYS_XKW',2,              &
                          long_name='XKW for ecosys fluxes',           &
                          units='cm/s', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ECOSYS_XKW_2,'ECOSYS_XKW_2',2,          &
                          long_name='XKW for ecosys fluxes',           &
                          units='cm/s', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_ECOSYS_ATM_PRESS,'ECOSYS_ATM_PRESS',2,  &
                          long_name='Atmospheric Pressure for ecosys fluxes', &
                          units='atmospheres', grid_loc='2110',        &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_PV_O2,'PV_O2',2,                        &
                          long_name='PV_O2',                           &
                          units='cm/s', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_SCHMIDT_O2,'SCHMIDT_O2',2,              &
                          long_name='O2 Schmidt Number',               &
                          units='none', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_O2SAT,'O2SAT',2,                        &
                          long_name='O2 Saturation',                   &
                          units='mmol/m^3', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CO2STAR,'CO2STAR',2,                    &
                          long_name='CO2 Star',                        &
                          units='mmol/m^3', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DCO2STAR,'DCO2STAR',2,                  &
                          long_name='D CO2 Star',                      &
                          units='mmol/m^3', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_pCO2SURF,'pCO2SURF',2,                  &
                          long_name='surface pCO2',                    &
                          units='ppmv', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DpCO2,'DpCO2',2,                        &
                          long_name='D pCO2',                          &
                          units='ppmv', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DpCO2_2,'DpCO2_2',2,                    &
                          long_name='D pCO2',                          &
                          units='ppmv', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_PV_CO2,'PV_CO2',2,                      &
                          long_name='CO2 Piston Velocity',             &
                          units='cm/s', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_SCHMIDT_CO2,'SCHMIDT_CO2',2,            &
                          long_name='CO2 Schmidt Number',              &
                          units='none', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DIC_GAS_FLUX,'FG_CO2',2,                &
                          long_name='DIC Surface Gas Flux',            &
                          units='mmol/m^3 cm/s', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DIC_GAS_FLUX_2,'FG_CO2_2',2,            &
                          long_name='DIC Surface Gas Flux',            &
                          units='mmol/m^3 cm/s', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_PH,'PH',2,                              &
                          long_name='Surface pH',                      &
                          units='none', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ATM_CO2,'ATM_CO2',2,                    &
                          long_name='Atmospheric CO2',                 &
                          units='ppmv', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_O2_GAS_FLUX_2,'STF_O2_2',2,             &
                          long_name='Dissolved Oxygen Surface Flux',   &
                          units='mmol/m^3 cm/s', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

!-----------------------------------------------------------------------
 
   allocate(ECO_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
   ECO_SFLUX_TAVG = c0
 
!-----------------------------------------------------------------------
!  The follow surface flux related fields are not with the above because
!  their values are directly available from STF in ecosys_tavg_forcing
!  and do not need to be stored in ECO_SFLUX_TAVG.
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_IRON_FLUX,'IRON_FLUX',2,                &
                          long_name='Iron Flux',                       &
                          units='mmol/m^2/s', grid_loc='2110',        &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_DUST_FLUX,'DUST_FLUX',2,                &
                          long_name='Dust Flux',                       &
                          units='g/cm^2/s', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_NOx_FLUX,'NOx_FLUX',2,                  &
                          long_name='Flux of NOx from Atmosphere',     &
                          units='nmol/cm^2/s', grid_loc='2110',        &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_NHy_FLUX,'NHy_FLUX',2,                  &
                          long_name='Flux of NHy from Atmosphere',     &
                          units='nmol/cm^2/s', grid_loc='2110',        &
                          coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!  nonstandard 2D fields
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_O2_ZMIN,'O2_ZMIN',2,                    &
                          long_name='Vertical Minimum of O2',          &
                          units='mmol/m^3', grid_loc='2111',           &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_O2_ZMIN_DEPTH,'O2_ZMIN_DEPTH',2,        &
                          long_name='Depth of Vertical Minimum of O2', &
                          units='cm', grid_loc='2111',                 &
                          coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!  nonstandard 3D fields
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_O2_PRODUCTION,'O2_PRODUCTION',3,        &
                          long_name='O2 Production',                   &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_O2_CONSUMPTION,'O2_CONSUMPTION',3,      &
                          long_name='O2 Consumption',                  &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_AOU,'AOU',3,                            &
                          long_name='Apparent O2 Utilization ',        &
                          units='mmol/m^3', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_PO4_RESTORE,'PO4_RESTORE',3,            &
                          long_name='PO4 Restoring',                   &
                          units='mmol/m^3', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_NO3_RESTORE,'NO3_RESTORE',3,            &
                          long_name='NO3 Restoring',                   &
                          units='mmol/m^3', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_SiO3_RESTORE,'SiO3_RESTORE',3,          &
                          long_name='SiO3 Restoring',                  &
                          units='mmol/m^3', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_PAR_avg,'PAR_avg',3,                    &
                          long_name='PAR Average over Model Cell',     &
                          units='w/m^2', grid_loc='3114',              &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_POC_FLUX_IN,'POC_FLUX_IN',3,            &
                          long_name='POC Flux into Cell',              &
                          units='mmol/m^3 cm/s', grid_loc='3111',      &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_POC_PROD,'POC_PROD',3,                  &
                          long_name='POC Production',                  &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_POC_REMIN,'POC_REMIN',3,                &
                          long_name='POC Remineralization',            &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_POC_ACCUM,'POC_ACCUM',3,                &
                          long_name='POC Accumulation',                &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CaCO3_FLUX_IN,'CaCO3_FLUX_IN',3,        &
                          long_name='CaCO3 flux into cell',            &
                          units='mmol/m^3 cm/s', grid_loc='3111',      &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CaCO3_PROD,'CaCO3_PROD',3,              &
                          long_name='CaCO3 Production',                &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_CaCO3_REMIN,'CaCO3_REMIN',3,            &
                          long_name='CaCO3 Remineralization',          &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_SiO2_FLUX_IN,'SiO2_FLUX_IN',3,          &
                          long_name='SiO2 Flux into Cell',             &
                          units='mmol/m^3 cm/s', grid_loc='3111',      &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_SiO2_PROD,'SiO2_PROD',3,                &
                          long_name='SiO2 Production',                 &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_SiO2_REMIN,'SiO2_REMIN',3,              &
                          long_name='SiO2 Remineralization',           &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_dust_FLUX_IN,'dust_FLUX_IN',3,          &
                          long_name='Dust Flux into Cell',             &
                          units='ng/s/m^2', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_dust_REMIN,'dust_REMIN',3,              &
                          long_name='Dust Remineralization',           &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_P_iron_FLUX_IN,'P_iron_FLUX_IN',3,      &
                          long_name='P_iron Flux into Cell',           &
                          units='mmol/m^3 cm/s', grid_loc='3111',      &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_P_iron_PROD,'P_iron_PROD',3,            &
                          long_name='P_iron Production',               &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_P_iron_REMIN,'P_iron_REMIN',3,          &
                          long_name='P_iron Remineralization',         &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_graze_sp,'graze_sp',3,                  &
                          long_name='Small Phyto Grazing',             &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_graze_diat,'graze_diat',3,              &
                          long_name='Diatom Grazing',                  &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_graze_diaz,'graze_diaz',3,              &
                          long_name='Diazotroph Grazing',              &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_graze_TOT,'graze_TOT',3,                &
                          long_name='Total Grazing',                   &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_sp_loss,'sp_loss',3,                    &
                          long_name='Small Phyto Loss',                &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diat_loss,'diat_loss',3,                &
                          long_name='Diatom Loss',                     &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diaz_loss,'diaz_loss',3,                &
                          long_name='Diaz Loss',                       &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_zoo_loss,'zoo_loss',3,                  &
                          long_name='Zooplankton Loss',                &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_sp_agg,'sp_agg',3,                      &
                          long_name='Small Phyto Aggregate',           &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diat_agg,'diat_agg',3,                  &
                          long_name='Diatom Aggregate',                &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_photoC_sp,'photoC_sp',3,                &
                          long_name='Small Phyto C Fixation',          &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_CaCO3_form,'CaCO3_form',3,              &
                          long_name='CaCO3 Formation',                 &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_photoC_diat,'photoC_diat',3,            &
                          long_name='Diatom C Fixation',               &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_photoC_diaz,'photoC_diaz',3,            &
                          long_name='Diaz C Fixation',                 &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_photoC_TOT,'photoC_TOT',3,              &
                          long_name='Total C Fixation',                &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_photoC_sp_zint,'photoC_sp_zint',2,      &
                          long_name='Small Phyto C Fixation Vertical Integral',&
                          units='mmol/m^3 cm/s', grid_loc='2111',      &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_CaCO3_form_zint,'CaCO3_form_zint',2,    &
                          long_name='CaCO3 Formation Vertical Integral',&
                          units='mmol/m^3 cm/s', grid_loc='2111',      &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_photoC_diat_zint,'photoC_diat_zint',2,  &
                          long_name='Diatom C Fixation Vertical Integral',&
                          units='mmol/m^3 cm/s', grid_loc='2111',      &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_photoC_diaz_zint,'photoC_diaz_zint',2,  &
                          long_name='Diaz C Fixation Vertical Integral',&
                          units='mmol/m^3 cm/s', grid_loc='2111',      &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_photoC_TOT_zint,'photoC_TOT_zint',2,    &
                          long_name='Total C Fixation Vertical Integral',&
                          units='mmol/m^3 cm/s', grid_loc='2111',      &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_photoC_NO3_sp,'photoC_NO3_sp',3,        &
                          long_name='Small Phyto C Fixation from NO3',      &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_photoC_NO3_sp_zint,'photoC_NO3_sp_zint',2,&
                          long_name='Small Phyto C Fixation from NO3 Vertical Integral',&
                          units='mmol/m^3 cm/s', grid_loc='2111',      &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_photoC_NO3_diat,'photoC_NO3_diat',3,    &
                          long_name='Diatom C Fixation from NO3',           &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_photoC_NO3_diat_zint,'photoC_NO3_diat_zint',2,&
                          long_name='Diatom C Fixation from NO3 Vertical Integral',&
                          units='mmol/m^3 cm/s', grid_loc='2111',      &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_photoC_NO3_diaz,'photoC_NO3_diaz',3,    &
                          long_name='Diaz C Fixation from NO3',             &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_photoC_NO3_diaz_zint,'photoC_NO3_diaz_zint',2,&
                          long_name='Diaz C Fixation from NO3 Vertical Integral',&
                          units='mmol/m^3 cm/s', grid_loc='2111',      &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_photoC_NO3_TOT,'photoC_NO3_TOT',3,      &
                          long_name='Total C Fixation from NO3',            &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_photoC_NO3_TOT_zint,'photoC_NO3_TOT_zint',2,&
                          long_name='Total C Fixation from NO3 Vertical Integral',&
                          units='mmol/m^3 cm/s', grid_loc='2111',      &
                          coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!  MORE nonstandard 3D fields
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_DOC_prod,'DOC_prod',3,                  &
                          long_name='DOC Production',                  &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DOC_remin,'DOC_remin',3,                &
                          long_name='DOC Remineralization',            &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DON_prod,'DON_prod',3,                  &
                          long_name='DON Production',                  &
                          units='mmol/m^3/s', grid_loc='3114',          &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_DON_remin,'DON_remin',3,                &
                          long_name='DON Remineralization',            &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_DOFe_prod,'DOFe_prod',3,                &
                          long_name='DOFe Production',                 &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_DOFe_remin,'DOFe_remin',3,              &
                          long_name='DOFe Remineralization',           &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_DOP_prod,'DOP_prod',3,                  &
                          long_name='DOP Production',                  &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_DOP_remin,'DOP_remin',3,                &
                          long_name='DOP Remineralization',            &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_Fe_scavenge,'Fe_scavenge',3,            &
                          long_name='Iron Scavenging',                 &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_Fe_scavenge_rate,'Fe_scavenge_rate',3,  &
                          long_name='Iron Scavenging Rate',            &
                          units='1/y', grid_loc='3111',                &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_sp_N_lim,'sp_N_lim',3,                  &
                          long_name='Small Phyto N Limitation',        &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_sp_Fe_lim,'sp_Fe_lim',3,                &
                          long_name='Small Phyto Fe Limitation',       &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_sp_PO4_lim,'sp_PO4_lim',3,              &
                          long_name='Small Phyto PO4 Limitation',      &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_sp_light_lim,'sp_light_lim',3,          &
                          long_name='Small Phyto Light Limitation',    &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diat_N_lim,'diat_N_lim',3,              &
                          long_name='Diatom N Limitation',             &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diat_Fe_lim,'diat_Fe_lim',3,            &
                          long_name='Diatom Fe Limitation',            &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diat_PO4_lim,'diat_PO4_lim',3,          &
                          long_name='Diatom PO4 Limitation',           &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diat_SiO3_lim,'diat_SiO3_lim',3,        &
                          long_name='Diatom SiO3 Limitation',          &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diat_light_lim,'diat_light_lim',3,      &
                          long_name='Diatom Light Limitation',         &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diaz_Fe_lim,'diaz_Fe_lim',3,            &
                          long_name='Diaz Fe Limitation',              &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diaz_P_lim,'diaz_P_lim',3,              &
                          long_name='Diaz PO4 Limitation',             &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diaz_light_lim,'diaz_light_lim',3,      &
                          long_name='Diaz LIGHT Limitation',           &
                          units='none', grid_loc='3114',               &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_diaz_Nfix,'diaz_Nfix',3,                &
                          long_name='Diaz N Fixation',                 &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_bSi_form,'bSi_form',3,                  &
                          long_name='Diatom Si Uptake',                &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_NITRIF,'NITRIF',3,                      &
                          long_name='Nitrification',                   &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DENITRIF,'DENITRIF',3,                  &
                          long_name='Denitrification',                 &
                          units='mmol/m^3/s', grid_loc='3111',         &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_photoFe_sp,'photoFe_sp',3,              &
                          long_name='Small Phyto Fe Uptake',           &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_photoFe_diat,'photoFe_diat',3,          &
                          long_name='Diatom Fe Uptake',                &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_photoFe_diaz,'photoFe_diaz',3,          &
                          long_name='Diaz Fe Uptake',                  &
                          units='mmol/m^3/s', grid_loc='3114',         &
                          coordinates='TLONG TLAT z_t_150m time')

   call define_tavg_field(tavg_CO3,'CO3',3,                            &
                          long_name='Carbonate Ion Concentration',     &
                          units='mmol/m^3', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_HCO3,'HCO3',3,                          &
                          long_name='Bicarbonate Ion Concentration',   &
                          units='mmol/m^3', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_H2CO3,'H2CO3',3,                        &
                          long_name='Carbonic Acid Concentration',     &
                          units='mmol/m^3', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_pH_3D,'pH_3D',3,                        &
                          long_name='pH',                              &
                          units='none', grid_loc='3111',               &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_co3_sat_calc,'co3_sat_calc',3,          &
                          long_name='CO3 concentration at calcite saturation', &
                          units='mmol/m^3', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_zsatcalc,'zsatcalc',2,                  &
                          long_name='Calcite Saturation Depth',        &
                          units='cm', grid_loc='2110',                 &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_co3_sat_arag,'co3_sat_arag',3,          &
                          long_name='CO3 concentration at aragonite saturation', &
                          units='mmol/m^3', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_zsatarag,'zsatarag',2,                  &
                          long_name='Aragonite Saturation Depth',      &
                          units='cm', grid_loc='2110',                 &
                          coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_init_tavg

!***********************************************************************
!BOP
! !IROUTINE: ecosys_set_interior
! !INTERFACE:

 subroutine ecosys_set_interior(k, TEMP_OLD, TEMP_CUR, SALT_OLD, SALT_CUR, &
    TRACER_MODULE_OLD, TRACER_MODULE_CUR, DTRACER_MODULE, this_block)

! !DESCRIPTION:
!  Compute time derivatives for ecosystem state variables
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k                   ! vertical level index

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      TEMP_OLD,          &! old potential temperature (C)
      TEMP_CUR,          &! current potential temperature (C)
      SALT_OLD,          &! old salinity (msu)
      SALT_CUR            ! current salinity (msu)

  !real (r8), dimension(nx_block,ny_block,km,ecosys_tracer_cnt), intent(in) :: &
   real (r8), dimension(:,:,:,:), intent(in) :: &
      TRACER_MODULE_OLD, &! old tracer values
      TRACER_MODULE_CUR   ! current tracer values

   type (block), intent(in) :: &
      this_block          ! block info for the current block

! !OUTPUT PARAMETERS:

  !real (r8), dimension(nx_block,ny_block,ecosys_tracer_cnt), intent(out) :: &
   real (r8), dimension(:,:,:), intent(out) :: &
      DTRACER_MODULE      ! computed source/sink terms

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: &
      subname = 'ecosys_mod:ecosys_set_interior'

   real (r8), parameter :: &
      epsC      = 1.00e-8, & ! small C concentration (mmol C/m^3)
      epsTinv   = 3.17e-8, & ! small inverse time scale (1/year) (1/sec)
      epsnondim = 1.00e-6    ! small non-dimensional number (non-dim)

   type(sinking_particle), save :: &
      POC,            & ! base units = nmol C
      P_CaCO3,        & ! base units = nmol CaCO3
      P_SiO2,         & ! base units = nmol SiO2
      dust,           & ! base units = g
      P_iron            ! base units = nmol Fe

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), save :: &
      QA_dust_def,    & ! incoming deficit in the QA(dust) POC flux
      ZSATCALC,       & ! Calcite Saturation Depth
      ZSATARAG,       & ! Aragonite Saturation Depth
      CO3_CALC_ANOM_km1,&! CO3 concentration above calcite saturation at k-1
      CO3_ARAG_ANOM_km1 ! CO3 concentration above aragonite saturation at k-1

   real (r8), dimension(nx_block,ny_block) :: &
      TEMP,           & ! local copy of model TEMP
      SALT,           & ! local copy of model SALT
      DIC_loc,        & ! local copy of model DIC
      ALK_loc,        & ! local copy of model ALK
      PO4_loc,        & ! local copy of model PO4
      NO3_loc,        & ! local copy of model NO3
      SiO3_loc,       & ! local copy of model SiO3
      NH4_loc,        & ! local copy of model NH4
      Fe_loc,         & ! local copy of model Fe
      O2_loc,         & ! local copy of model O2
      DOC_loc,        & ! local copy of model DOC
      spC_loc,        & ! local copy of model spC
      spChl_loc,      & ! local copy of model spChl
      spCaCO3_loc,    & ! local copy of model spCaCO3
      diatC_loc,      & ! local copy of model diatC
      diatChl_loc,    & ! local copy of model diatChl
      zooC_loc,       & ! local copy of model zooC
      spFe_loc,       & ! local copy of model spFe
      diatSi_loc,     & ! local copy of model diatSi
      diatFe_loc,     & ! local copy of model diatFe
      diazC_loc,      & ! local copy of model diazC
      diazChl_loc,    & ! local copy of model diazChl
      diazFe_loc,     & ! local copy of model diazFe
      DON_loc,        & ! local copy of model DON
      DOFe_loc,       & ! local copy of model DOFe
      DOP_loc           ! local copy of model DOP

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1,WORK2,WORK3 ! temporaries

   real (r8) :: &
      C_loss_thres,   & ! bio-C threshold at which losses go to zero (mmol C/m^3)
      f_loss_thres      ! fraction of grazing loss reduction at depth

   real (r8), dimension(nx_block,ny_block) :: &
      PAR_in,         & ! photosynthetically available radiation (W/m^2)
      KPARdz,         & ! PAR adsorption coefficient (non-dim)
      PAR_avg,        & ! average PAR over mixed layer depth (W/m^2)
      DOC_prod,       & ! production of DOC (mmol C/m^3/sec)
      DOC_remin,      & ! remineralization of DOC (mmol C/m^3/sec)
      NITRIF,         & ! nitrification (NH4 -> NO3) (mmol N/m^3/sec)
      DENITRIF,       & ! nitrification (NO3 -> N2) (mmol N/m^3/sec)
      RESTORE           ! restoring terms for nutrients (mmol ./m^3/sec)

   real (r8), dimension(nx_block,ny_block) :: &
      z_umax,         & ! max. zoo growth rate on sp at local T (1/sec)
      diat_umax,      & ! max. zoo growth rate on diatoms at local T (1/sec)
      z_mort,         & ! zoo respiration loss, (1/sec/((mmol C/m3))
      C_loss_diaz,    & ! bio-C threshold at which losses go to zero (mmol C/m^3)
      z_mort2,        & ! zoo quad mort rate, tohigherlevels (1/sec/((mmol C/m3))
      diaz_umax         ! max. zoo growth rate on diazotrophs at local T (1/sec)


   real (r8), dimension(nx_block,ny_block) :: &
      thetaC_sp,      & ! local Chl/C ratio in small phyto (mg Chl/mmol C)
      thetaC_diat,    & ! local Chl/C ratio in diatoms (mg Chl/mmol C)
      QCaCO3,         & ! small phyto CaCO3/C ratio (mmol CaCO3/mmol C)
      Tfunc,          & ! temp response function GD98 (non-dim)
      VNO3_sp,        & ! small phyto NO3 uptake rate (non-dim)
      VNH4_sp,        & ! small phyto NH4 uptake rate (non-dim)
      VNtot_sp,       & ! small phyto total N uptake rate (non-dim)
      VFeC_sp,        & ! ??? small phyto C-specific iron uptake (non-dim)
      VPO4_sp,        & ! ??? small phyto C-specific po4 uptake (non-dim)
      f_nut,          & ! nut limitation factor, modifies C fixation (non-dim)
      PCmax,          & ! max value of PCphoto at temperature TEMP (1/sec)
      PCphoto_sp,     & ! small C-specific rate of photosynth. (1/sec)
      photoC_sp,      & ! small phyto C-fixation (mmol C/m^3/sec)
      NO3_V_sp,       & ! nitrate uptake by small phyto (mmol NO3/m^3/sec)
      NH4_V_sp,       & ! ammonium uptake by small phyto (mmol NH4/m^3/sec)
      VNC_sp,         & ! small phyto C-specific N uptake rate (mmol N/mmol C/sec)
      pChl,           & ! Chl synth. regulation term (mg Chl/mmol N)
      photoacc_sp,    & ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
      CaCO3_PROD,     & ! prod. of CaCO3 by small phyto (mmol CaCO3/m^3/sec)
      VNO3_diat,      & ! diatom nitrate uptake rate (non-dim)
      VNH4_diat,      & ! diatom ammonium uptake rate (non-dim)
      VNtot_diat,     & ! diatom total N uptake rate (non-dim)
      VFeC_diat,      & ! diatom C-specific iron uptake (non-dim)
      VPO4_diat,      & ! diatom C-specific PO4 uptake (non-dim)
      VSiO3_diat,     & ! C-specific SiO3 uptake for diatoms (non-dim)
      PCphoto_diat,   & ! diatom C-specific rate of photosynth. (1/sec)
      photoC_diat,    & ! diatom C-fixation (mmol C/m^3/sec)
      NO3_V_diat,     & ! nitrate uptake by diatoms (mmol NO3/m^3/sec)
      NH4_V_diat,     & ! ammonium uptake by diatoms (mmol NH4/m^3/sec)
      VNC_diat,       & ! diatom C-specific N uptake rate (mmol N/mmol C/sec)
      photoacc_diat,  & ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
      reduceV,        & ! factor in nutrient uptake (mmol C/m^3)^2
      graze_sp,       & ! grazing rate on small phyto (mmol C/m^3/sec)
      graze_sp_zoo,   & ! graze_sp routed to zoo (mmol C/m^3/sec)
      graze_sp_poc,   & ! graze_sp routed to poc (mmol C/m^3/sec)
      graze_sp_doc,   & ! graze_sp routed to doc (mmol C/m^3/sec)
      graze_sp_dic      ! graze_sp routed to dic (mmol C/m^3/sec)

   real (r8), dimension(nx_block,ny_block) :: & ! max of 39 continuation lines
      graze_diat,     & ! grazing rate on diatoms (mmol C/m^3/sec)
      graze_diat_zoo, & ! graze_diat routed to zoo (mmol C/m^3/sec)
      graze_diat_poc, & ! graze_diat routed to poc (mmol C/m^3/sec)
      graze_diat_doc, & ! graze_diat routed to doc (mmol C/m^3/sec)
      graze_diat_dic, & ! graze_diat routed to dic (mmol C/m^3/sec)
      Pprime,         & ! used to limit phyto mort at low biomass (mmol C/m^3)
      sp_loss,        & ! small phyto non-grazing mort (mmol C/m^3/sec)
      sp_loss_poc,    & ! sp_loss routed to poc (mmol C/m^3/sec)
      sp_loss_doc,    & ! sp_loss routed to doc (mmol C/m^3/sec)
      sp_loss_dic,    & ! sp_loss routed to dic (mmol C/m^3/sec)
      sp_agg,         & ! small phyto agg loss (mmol C/m^3/sec)
      diat_loss,      & ! diatom non-grazing mort (mmol C/m^3/sec)
      diat_loss_poc,  & ! diat_loss routed to poc (mmol C/m^3/sec)
      diat_loss_doc,  & ! diat_loss routed to doc (mmol C/m^3/sec)
      diat_loss_dic,  & ! diat_loss routed to dic (mmol C/m^3/sec)
      diat_agg,       & ! diatom agg (mmol C/m^3/sec)
      f_zoo_detr,     & ! frac of zoo losses into large detrital pool (non-dim)
      Fe_scavenge,    & ! loss of dissolved iron, scavenging (mmol Fe/m^3/sec)
      Zprime,         & ! used to limit zoo mort at low biomass (mmol C/m^3)
      zoo_loss,       & ! mortality & higher trophic grazing on zooplankton (mmol C/m^3/sec)
      zoo_loss_doc,   & ! zoo_loss routed to doc (mmol C/m^3/sec)
      zoo_loss_dic,   & ! zoo_loss routed to dic (mmol C/m^3/sec)
      light_lim,      & ! light limitation factor
      Qsi,            & ! Diatom initial Si/C ratio (mmol Si/mmol C)
      gQsi,           & ! diatom Si/C ratio for growth (new biomass)
      Qfe_sp,         & ! small phyto init fe/C ratio (mmolFe/mmolC)
      gQfe_sp,        & ! small phyto fe/C for growth
      Qfe_diat,       & ! diatom init fe/C ratio
      gQfe_diat,      & ! diatom fe/C ratio for growth
      Qfe_diaz,       & ! diazotrophs init fe/C ratio
      gQfe_diaz         ! diazotroph fe/C ratio for new growth

   real (r8), dimension(nx_block,ny_block) :: & ! max of 39 continuation lines
      PCphoto_diaz,   & ! diazotroph C-specific rate of photosynth. (1/sec)
      photoC_diaz,    & ! diazotroph C-fixation (mmol C/m^3/sec)
      Vfec_diaz,      & ! diazotroph C-specific iron uptake (non-dim)
      Vpo4_diaz,      & ! diazotroph C-specific po4 uptake (non-dim)
      photoacc_diaz,  & ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
      Vnc_diaz,       & ! diazotroph C-specific N uptake rate (mmol N/mmol C/sec)
      diaz_Nexcrete,  & ! diazotroph fixed N excretion
      diaz_Nfix,      & ! diazotroph total Nitrogen fixation (mmol N/m^3/sec)
      thetaC_diaz,    & ! local Chl/C ratio in diazotrophs (mg Chl/mmol C)
      photoFe_diaz,   & ! iron uptake by diazotrophs (mmolFe/m^3/sec)
      photoFe_diat,   & ! iron uptake by diatoms
      photoFe_sp,     & ! iron uptake by small phyto
      photoSi_diat,   & ! silicon uptake by diatoms (mmol Si/m^3/sec)
      remaining_diazP,& ! used in routing P from diazotroph losses
      diaz_loss,      & ! diazotroph non-grazing mort (mmol C/m^3/sec)
      diaz_loss_doc,  & ! mortality routed to DOM pool
      diaz_loss_poc,  & ! mortalitiy routed to POC pool
      diaz_loss_dic,  & ! mortality routed to remin
      diaz_loss_dop,  & ! P from mort routed to DOP pool
      diaz_loss_dip,  & ! P from mort routed to remin
                        ! P from diaz losses must be routed differently than
                        ! other elements to ensure that sinking detritus and
                        ! zooplankton pools get their fixed P/C ratios, the
                        ! remaining P is split evenly between DOP and PO4
      graze_diaz,     & ! grazing rate on diazotrophs (mmol C/m^3/sec)
      graze_diaz_poc, & ! grazing routed to sinking detr (mmol C/m^3/sec)
      graze_diaz_doc, & ! grazing routed to DOC (mmol C/m^3/sec)
      graze_diaz_dic, & ! grazing routed to remin (mmol C/m^3/sec)
      graze_diaz_zoo, & ! grazing routed to new zoo biomass
      DON_remin,      & ! portion of DON remineralized
      DOFe_remin,     & ! portion of DOFe remineralized
      DOP_remin,      & ! portion of DOP remineralized
      DOM_remin,      & ! fraction of DOM remineralized at current TEMP
      Fe_scavenge_rate  ! annual scavenging rate of iron as % of ambient

   real (r8), dimension(nx_block,ny_block) :: & ! max of 39 continuation lines
      DON_prod,       & ! production of dissolved organic N
      DOFe_prod,      & ! produciton of dissolved organic Fe
      DOP_prod,       & ! production of dissolved organic P
      po4_V_sp,       & ! sp uptake of po4 (mmol po4/m3/s)
      po4_V_diaz,     & ! diaz uptake of po4 (mmol po4/m3/s)
      dop_V_sp,       & ! sp uptake of dop (mmol dop/m3/s)
      dop_V_diaz,     & ! diaz uptake of dop (mmol dop/m3/s)
      cks,            & ! constant used in Fe quota modification
      cksi,           & ! constant used in Si quota modification
      VNO3_diaz,      & ! relative NO3 uptake by diazotrophs
      VNH4_diaz,      & ! relative NH4 uptake by diazotrophs
      VNtot_diaz,     & ! total relative N uptake by diazotrophs
      photoNO3_diaz,  & ! total nitrate uptake by diazotrophs
      photoNH4_diaz,  & ! total NH4 uptake by diazotrophs
      O2_PRODUCTION,  & ! O2 production
      O2_CONSUMPTION    ! O2 consumption

   real (r8), dimension(nx_block,ny_block) :: &
      CO3,            &! carbonate ion
      HCO3,           &! bicarbonate ion
      H2CO3,          &! carbonic acid
      OMEGA_CALC,     &! solubility ratio for aragonite
      OMEGA_ARAG       ! solubility ratio for calcite

   integer (int_kind) :: &
      bid,            & ! local_block id
      kk                ! index for looping over k levels

!-----------------------------------------------------------------------

   bid = this_block%local_id

   call timer_start(ecosys_interior_timer, block_id=bid)

   DTRACER_MODULE = c0

!-----------------------------------------------------------------------
!  exit immediately if computations are not to be performed
!-----------------------------------------------------------------------

   if (.not. lsource_sink) then
      call timer_stop(ecosys_interior_timer, block_id=bid)
      return
   endif

!-----------------------------------------------------------------------
!  create local copies of model tracers
!  treat negative values as zero
!  apply mask to local copies
!-----------------------------------------------------------------------

   TEMP         = p5*(TEMP_OLD + TEMP_CUR)
   SALT         = p5*(SALT_OLD + SALT_CUR)*salt_to_ppt

   DIC_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,dic_ind) + &
                              TRACER_MODULE_CUR(:,:,k,dic_ind)))
   ALK_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,alk_ind) + &
                              TRACER_MODULE_CUR(:,:,k,alk_ind)))
   PO4_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,po4_ind) + &
                              TRACER_MODULE_CUR(:,:,k,po4_ind)))
   NO3_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,no3_ind) + &
                              TRACER_MODULE_CUR(:,:,k,no3_ind)))
   SiO3_loc     = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,sio3_ind) + &
                              TRACER_MODULE_CUR(:,:,k,sio3_ind)))
   NH4_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,nh4_ind) + &
                              TRACER_MODULE_CUR(:,:,k,nh4_ind)))
   Fe_loc       = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,fe_ind) + &
                              TRACER_MODULE_CUR(:,:,k,fe_ind)))
   O2_loc       = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,o2_ind) + &
                              TRACER_MODULE_CUR(:,:,k,o2_ind)))
   DOC_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,doc_ind) + &
                              TRACER_MODULE_CUR(:,:,k,doc_ind)))
   spC_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,spC_ind) + &
                              TRACER_MODULE_CUR(:,:,k,spC_ind)))
   spChl_loc    = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,spChl_ind) + &
                              TRACER_MODULE_CUR(:,:,k,spChl_ind)))
   spCaCO3_loc  = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,spCaCO3_ind) + &
                              TRACER_MODULE_CUR(:,:,k,spCaCO3_ind)))
   diatC_loc    = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,diatC_ind) + &
                              TRACER_MODULE_CUR(:,:,k,diatC_ind)))
   diatChl_loc  = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,diatChl_ind) + &
                              TRACER_MODULE_CUR(:,:,k,diatChl_ind)))
   zooC_loc     = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,zooC_ind) + &
                              TRACER_MODULE_CUR(:,:,k,zooC_ind)))
   spFe_loc     = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,spFe_ind) + &
                              TRACER_MODULE_CUR(:,:,k,spFe_ind)))
   diatSi_loc   = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,diatSi_ind) + &
                              TRACER_MODULE_CUR(:,:,k,diatSi_ind)))
   diatFe_loc   = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,diatFe_ind) + &
                              TRACER_MODULE_CUR(:,:,k,diatFe_ind)))
   diazC_loc    = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,diazC_ind) + &
                              TRACER_MODULE_CUR(:,:,k,diazC_ind)))
   diazChl_loc  = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,diazChl_ind) + &
                              TRACER_MODULE_CUR(:,:,k,diazChl_ind)))
   diazFe_loc   = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,diazFe_ind) + &
                              TRACER_MODULE_CUR(:,:,k,diazFe_ind)))
   DON_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,don_ind) + &
                              TRACER_MODULE_CUR(:,:,k,don_ind)))
   DOFe_loc     = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,dofe_ind) + &
                              TRACER_MODULE_CUR(:,:,k,dofe_ind)))
   DOP_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,dop_ind) + &
                              TRACER_MODULE_CUR(:,:,k,dop_ind)))

   where (.not. LAND_MASK(:,:,bid) .or. k > KMT(:,:,bid))
      PO4_loc      = c0
      NO3_loc      = c0
      SiO3_loc     = c0
      NH4_loc      = c0
      Fe_loc       = c0
      O2_loc       = c0
      DOC_loc      = c0
      spC_loc      = c0
      spChl_loc    = c0
      spCaCO3_loc  = c0
      diatC_loc    = c0
      diatChl_loc  = c0
      zooC_loc     = c0
      spFe_loc     = c0
      diatSi_loc   = c0
      diatFe_loc   = c0
      diazC_loc    = c0
      diazChl_loc  = c0
      diazFe_loc   = c0
      DON_loc      = c0
      DOFe_loc     = c0
      DOP_loc      = c0
   end where

!-----------------------------------------------------------------------
!  If any phyto box are zero, set others to zeros.
!-----------------------------------------------------------------------

   where (spC_loc == c0 .or. spChl_loc == c0 .or. spFe_loc == c0)
      spC_loc     = c0
      spChl_loc   = c0
      spCaCO3_loc = c0
      spFe_loc    = c0
   end where

   where (diatC_loc == c0 .or. diatChl_loc == c0 .or. diatFe_loc == c0 &
             .or. diatSi_loc == c0)
      diatC_loc   = c0
      diatChl_loc = c0
      diatFe_loc  = c0
      diatSi_loc  = c0
   end where

   where (diazC_loc == c0 .or. diazChl_loc == c0 .or. diazFe_loc == c0)
      diazC_loc   = c0
      diazChl_loc = c0
      diazFe_loc  = c0
   end where

!-----------------------------------------------------------------------
!  set local variables, with incoming ratios
!-----------------------------------------------------------------------

   thetaC_sp   = spChl_loc   / (spC_loc + epsC)
   thetaC_diat = diatChl_loc / (diatC_loc + epsC)
   thetaC_diaz = diazChl_loc / (diazC_loc + epsC)
   Qsi         = diatSi_loc  / (diatC_loc + epsC)
   Qfe_diat    = diatFe_loc  / (diatC_loc + epsC)
   Qfe_sp      = spFe_loc    / (spC_loc + epsC)
   Qfe_diaz    = diazFe_loc  / (diazC_loc + epsC)
   where (Qsi > gQsi_max) Qsi = gQsi_max

!-----------------------------------------------------------------------
!  DETERMINE NEW ELEMENTAL RATIOS FOR GROWTH (NEW BIOMASS)
!-----------------------------------------------------------------------

   gQsi      = gQsi_0
   gQfe_diat = gQfe_diat_0
   gQfe_sp   = gQfe_sp_0
   gQfe_diaz = gQfe_diaz_0
   cks       = 10._r8
   cksi      = 10._r8

!-----------------------------------------------------------------------
!  Modify these initial ratios under low ambient iron conditions
!-----------------------------------------------------------------------

   where (Fe_loc < cks * parm_diat_kfe)
      gQfe_diat = max((gQfe_diat * Fe_loc /(cks * parm_diat_kfe)), &
                      gQfe_diat_min)
   end where

   where ((Fe_loc < cksi * parm_diat_kfe) .and. (Fe_loc > c0) .and. &
          (SiO3_loc > (cksi * parm_diat_kSiO3)))
      gQsi = min((gQsi*cksi*parm_diat_kfe/Fe_loc), gQsi_max)
   end where

   where (Fe_loc == c0)
      gQsi = gQsi_max
   end where

   where (Fe_loc < cks * parm_sp_kfe)
      gQfe_sp = max((gQfe_sp*Fe_loc/(cks * parm_sp_kfe)), &
                    gQfe_sp_min)
   end where

   where (Fe_loc < cks * diaz_kFe)
      gQfe_diaz = max((gQfe_diaz*Fe_loc/(cks * diaz_kFe)), &
                      gQfe_diaz_min)
   end where

!-----------------------------------------------------------------------
!  Modify the initial si/C ratio under low ambient Si conditions
!-----------------------------------------------------------------------

   where (SiO3_loc < (cksi * parm_diat_kSiO3))
      gQsi = max((gQsi*SiO3_loc/(cksi * parm_diat_kSiO3)), &
                 gQsi_min)
   end where

!-----------------------------------------------------------------------
!  various k==1 initializations
!
!  0.45   fraction of incoming SW -> PAR (non-dim)
!-----------------------------------------------------------------------

   if (k == 1) then
      call init_particulate_terms(POC, P_CaCO3, P_SiO2, dust, P_iron, &
                                  QA_dust_def, dust_FLUX_IN, this_block)
   endif

!-----------------------------------------------------------------------
!  QCaCO3 is the percentage of sp organic matter which is associated
!  with coccolithophores
!-----------------------------------------------------------------------

   QCaCO3 = spCaCO3_loc / (spC_loc + epsC)
   where (QCaCO3 > QCaCO3_max) QCaCO3 = QCaCO3_max

!-----------------------------------------------------------------------
!  compute PAR related quantities
!
!  0.03e-2   atten. coeff. per unit chlorophyll (1/cm/(mg Chl/m^3))
!  0.04e-2   atten. coeff. for water (1/cm)
!-----------------------------------------------------------------------

   PAR_in = PAR_out(:,:,bid)
   where (.not. LAND_MASK(:,:,bid) .or. k > KMT(:,:,bid)) PAR_in = c0

   if (partial_bottom_cells) then
      KPARdz = (k_chl * (spChl_loc + diatChl_loc + &
                         diazChl_loc) + k_h2o) * DZT(:,:,k,bid)
   else
      KPARdz = (k_chl * (spChl_loc + diatChl_loc + &
                         diazChl_loc) + k_h2o) * dz(k)
   endif

   PAR_out(:,:,bid) = PAR_in * exp(-KPARdz)
   PAR_avg = PAR_in * (c1 - exp(-KPARdz)) / KPARdz

!-----------------------------------------------------------------------
!  Tref = 30.0 reference temperature (deg. C)
!
!  Using q10 formulation with Q10 value of 2.0 (Doney et al., 1996).
!-----------------------------------------------------------------------

   Tfunc = Q_10**(((TEMP + T0_Kelvin)-(Tref + T0_Kelvin))/c10)

!-----------------------------------------------------------------------
!  modify growth mort rates by Tfunc
!-----------------------------------------------------------------------

   z_umax    = parm_z_umax_0    * Tfunc
   diat_umax = parm_diat_umax_0 * Tfunc
   z_mort2   = parm_z_mort2_0   * Tfunc
   z_mort    = parm_z_mort_0    * Tfunc
   diaz_umax = parm_diaz_umax_0 * Tfunc

   DOM_remin = parm_sd_remin_0

!-----------------------------------------------------------------------
!  Get relative nutrient uptake rates for phytoplankton,
!  min. relative uptake rate modifies C fixation in the manner
!  that the min. cell quota does in GD98.
!-----------------------------------------------------------------------

   VNO3_sp = (NO3_loc / parm_sp_kNO3) / &
      (c1 + (NO3_loc / parm_sp_kNO3) + (NH4_loc / parm_sp_kNH4))

   VNH4_sp = (NH4_loc / parm_sp_kNH4) / &
      (c1 + (NO3_loc / parm_sp_kNO3) + (NH4_loc / parm_sp_kNH4))

   VNtot_sp = VNO3_sp + VNH4_sp

   call accumulate_tavg_field(VNtot_sp, tavg_sp_N_lim,bid,k)

!-----------------------------------------------------------------------
!  get relative Fe uptake by phytoplankton
!  get relative P uptake rates for phytoplankton
!-----------------------------------------------------------------------

   VFeC_sp = Fe_loc / (Fe_loc + parm_sp_kFe)

   call accumulate_tavg_field(VFeC_sp, tavg_sp_Fe_lim,bid,k)

!-----------------------------------------------------------------------
!  Partition diaz P uptake in same manner to no3/nh4 partition
!-----------------------------------------------------------------------

   VPO4_sp = PO4_loc / (PO4_loc + parm_sp_kPO4)

   call accumulate_tavg_field(VPO4_sp, tavg_sp_PO4_lim,bid,k)

!-----------------------------------------------------------------------
!  Small Phytoplankton C-fixation - given light and relative uptake rates
!  determine small phyto nutrient limitation factor for carbon fixation
!-----------------------------------------------------------------------

   f_nut = min(VNtot_sp, VFeC_sp)
   f_nut = min(f_nut, VPO4_sp)

!-----------------------------------------------------------------------
!  get small phyto photosynth. rate, phyto C biomass change, photoadapt
!-----------------------------------------------------------------------

   PCmax = PCref * f_nut * Tfunc

   light_lim = (c1 - exp((-c1 * parm_alphaChlsp * thetaC_sp * PAR_avg) / &
                         (PCmax + epsTinv)))
   PCphoto_sp = PCmax * light_lim

   call accumulate_tavg_field(light_lim, tavg_sp_light_lim,bid,k)

   photoC_sp = PCphoto_sp * spC_loc

!-----------------------------------------------------------------------
!  Get nutrient uptakes by small phyto based on calculated C fixation
!  total N uptake (VNC_sp) is used in photoadaption
!-----------------------------------------------------------------------

   where (VNtot_sp > c0)
      NO3_V_sp = (VNO3_sp / VNtot_sp) * photoC_sp * Q
      NH4_V_sp = (VNH4_sp / VNtot_sp) * photoC_sp * Q
      VNC_sp = PCphoto_sp * Q
   elsewhere
      NO3_V_sp = c0
      NH4_V_sp = c0
      VNC_sp = c0
      photoC_sp = c0
   end where

   po4_V_sp = photoC_sp * Qp
   photoFe_sp = photoC_sp * gQfe_sp

   call accumulate_tavg_field(photoFe_sp, tavg_photoFe_sp,bid,k)

!-----------------------------------------------------------------------
!  calculate pChl, (used in photoadapt., GD98)
!  2.3   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
!  GD 98 Chl. synth. term
!-----------------------------------------------------------------------

   WORK1 = parm_alphaChlsp * thetaC_sp * PAR_avg
   where (WORK1 > c0)
      pChl = thetaN_max_sp * PCphoto_sp / WORK1
      photoacc_sp = (pChl * VNC_sp / thetaC_sp) * spChl_loc
   elsewhere
      photoacc_sp = c0
   end where

!-----------------------------------------------------------------------
!  CaCO3 Production, parameterized as function of small phyto production
!  decrease CaCO3 as function of nutrient limitation decrease CaCO3 prod
!  at low temperatures increase CaCO3 prod under bloom conditions
!  maximum calcification rate is 40% of primary production
!-----------------------------------------------------------------------

   CaCO3_PROD = f_prod_sp_CaCO3 * photoC_sp
   CaCO3_PROD = CaCO3_PROD * f_nut

   where (TEMP < CaCO3_temp_thres1)  &
      CaCO3_PROD = CaCO3_PROD * max((TEMP-CaCO3_temp_thres2), c0) / &
                   (CaCO3_temp_thres1-CaCO3_temp_thres2)

   where (spC_loc > CaCO3_sp_thres)  &
      CaCO3_PROD = min((CaCO3_PROD*spC_loc/CaCO3_sp_thres), &
                       (f_photosp_CaCO3*photoC_sp))

   call accumulate_tavg_field(CaCO3_PROD, tavg_CaCO3_form,bid,k)

   if (accumulate_tavg_now(tavg_CaCO3_form_zint)) then
      if (partial_bottom_cells) then
         WORK1 = DZT(:,:,k,bid) * CaCO3_PROD
      else
         WORK1 = dz(k) * CaCO3_PROD
      endif
      call accumulate_tavg_field(WORK1, tavg_CaCO3_form_zint,bid,k)
   endif

!-----------------------------------------------------------------------
!  Relative uptake rates for diatoms nitrate is VNO3, ammonium is VNH4
!-----------------------------------------------------------------------

   VNO3_diat = (NO3_loc / parm_diat_kNO3) / &
      (c1 + (NO3_loc / parm_diat_kNO3) + (NH4_loc / parm_diat_kNH4))

   VNH4_diat = (NH4_loc / parm_diat_kNH4) / &
      (c1 + (NO3_loc / parm_diat_kNO3) + (NH4_loc / parm_diat_kNH4))

   VNtot_diat = VNO3_diat + VNH4_diat

   call accumulate_tavg_field(VNtot_diat, tavg_diat_N_lim,bid,k)

!-----------------------------------------------------------------------
!  get relative Fe uptake by diatoms
!  get relative P uptake rates for diatoms
!  get relative SiO3 uptake rate for diatoms
!-----------------------------------------------------------------------

   VFeC_diat = Fe_loc / (Fe_loc + parm_diat_kFe)

   call accumulate_tavg_field(VFeC_diat, tavg_diat_Fe_lim,bid,k)

   VPO4_diat = PO4_loc / (PO4_loc + parm_diat_kPO4)

   call accumulate_tavg_field(VPO4_diat, tavg_diat_PO4_lim,bid,k)

   VSiO3_diat = SiO3_loc / (SiO3_loc + parm_diat_kSiO3)

   call accumulate_tavg_field(VSiO3_diat, tavg_diat_SiO3_lim,bid,k)

!-----------------------------------------------------------------------
!  Diatom carbon fixation and photoadapt.
!  determine diatom nutrient limitation factor for carbon fixation
!-----------------------------------------------------------------------

   f_nut = min(VNtot_diat, VFeC_diat)
   f_nut = min(f_nut, VSiO3_diat)
   f_nut = min(f_nut, VPO4_diat)

!-----------------------------------------------------------------------
!  get diatom photosynth. rate, phyto C biomass change, photoadapt
!-----------------------------------------------------------------------

   PCmax = PCref * f_nut * Tfunc

   light_lim = (c1 - exp((-c1 * parm_alphaChl * thetaC_diat * PAR_avg) / &
                         (PCmax + epsTinv)))
   PCphoto_diat = PCmax * light_lim

   call accumulate_tavg_field(light_lim, tavg_diat_light_lim,bid,k)

   photoC_diat = PCphoto_diat * diatC_loc

!-----------------------------------------------------------------------
!  Get nutrient uptake by diatoms based on C fixation
!-----------------------------------------------------------------------

   where (VNtot_diat > c0)
      NO3_V_diat = (VNO3_diat / VNtot_diat) * photoC_diat * Q
      NH4_V_diat = (VNH4_diat / VNtot_diat) * photoC_diat * Q
      VNC_diat = PCphoto_diat * Q
   elsewhere
      NO3_V_diat = c0
      NH4_V_diat = c0
      VNC_diat = c0
      photoC_diat = c0
   end where

   photoFe_diat = photoC_diat * gQfe_diat
   photoSi_diat = photoC_diat * gQsi

   call accumulate_tavg_field(photoFe_diat, tavg_photoFe_diat,bid,k)

   call accumulate_tavg_field(photoSi_diat, tavg_bSi_form,bid,k)

!-----------------------------------------------------------------------
!  calculate pChl, (used in photoadapt., GD98)
!  3.0   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
!  GD 98 Chl. synth. term
!-----------------------------------------------------------------------

   WORK1 = parm_alphaChl * thetaC_diat * PAR_avg
   where (WORK1 > c0)
      pChl = thetaN_max_diat * PCphoto_diat / WORK1
      photoacc_diat = (pChl * VNC_diat / thetaC_diat) * diatChl_loc
   elsewhere
      photoacc_diat = c0
   end where

!-----------------------------------------------------------------------
!  get relative Fe uptake by diazotrophs
!  get relative P uptake rates for diazotrophs
!-----------------------------------------------------------------------

   Vfec_diaz = Fe_loc/(Fe_loc + diaz_kFe)

   call accumulate_tavg_field(Vfec_diaz, tavg_diaz_Fe_lim,bid,k)

   Vpo4_diaz = PO4_loc / (PO4_loc + diaz_kPO4)

   call accumulate_tavg_field(Vpo4_diaz, tavg_diaz_P_lim,bid,k)

   f_nut = min(Vpo4_diaz, Vfec_diaz)

!-----------------------------------------------------------------------
!  get diazotroph photosynth. rate, phyto C biomass change
!-----------------------------------------------------------------------

   PCmax = PCrefDiaz * f_nut * Tfunc
   where (TEMP < diaz_temp_thres) PCmax = c0

   light_lim = (c1 - exp((-c1 * parm_alphaDiaz * thetaC_diaz * PAR_avg) / &
                         (PCmax + epsTinv)))
   PCphoto_diaz = PCmax * light_lim

   call accumulate_tavg_field(light_lim, tavg_diaz_light_lim,bid,k)

   photoC_diaz = PCphoto_diaz * diazC_loc

!-----------------------------------------------------------------------
!  Get Fe and po4 uptake by diazotrophs based on C fixation
!-----------------------------------------------------------------------

   po4_V_diaz = photoC_diaz * Qp_diaz
   photoFe_diaz = photoC_diaz * gQfe_diaz

   call accumulate_tavg_field(photoFe_diaz, tavg_photoFe_diaz,bid,k)

!-----------------------------------------------------------------------
!  Relative uptake rates for diazotrophs nitrate is VNO3, ammonium is VNH4
!-----------------------------------------------------------------------

   VNO3_diaz = (NO3_loc / diaz_kNO3) / &
      (c1 + (NO3_loc / diaz_kNO3) + (NH4_loc / diaz_kNH4))

   VNH4_diaz = (NH4_loc / diaz_kNH4) / &
      (c1 + (NO3_loc / diaz_kNO3) + (NH4_loc / diaz_kNH4))

   VNtot_diaz = c1

!------------------------------------------------------------------------
!  Compute inoranic N uptake by diazotrophs.
!------------------------------------------------------------------------

   WORK1 = photoC_diaz * Q
   photoNO3_diaz = VNO3_diaz / VNtot_diaz *WORK1
   photoNH4_diaz = VNH4_diaz / VNtot_diaz *WORK1

!-----------------------------------------------------------------------
!  Get N fixation by diazotrophs based on C fixation,
!  Diazotrophs fix more than they need then 30% is excreted
!-----------------------------------------------------------------------

   diaz_Nfix     = (WORK1 * r_Nfix_photo) - photoNO3_diaz - photoNH4_diaz
   diaz_Nexcrete = diaz_Nfix + photoNO3_diaz + photoNH4_diaz - WORK1

   Vnc_diaz = PCphoto_diaz * Q

   call accumulate_tavg_field(diaz_Nfix, tavg_diaz_Nfix,bid,k)

!-----------------------------------------------------------------------
!  calculate pChl, (used in photoadapt., GD98)
!  3.4   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
!  GD 98 Chl. synth. term
!-----------------------------------------------------------------------

   WORK1 = parm_alphaDiaz * thetaC_diaz * PAR_avg
   where (WORK1 > c0)
      pChl = thetaN_max_diaz * PCphoto_diaz / WORK1
      photoacc_diaz = (pChl * Vnc_diaz / thetaC_diaz) * diazChl_loc
   elsewhere
      photoacc_diaz = c0
   end where

!-----------------------------------------------------------------------
!  CALCULATE GRAZING AND OTHER MORT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  calculate the loss threshold interpolation factor
!-----------------------------------------------------------------------

   if (zt(k) > thres_z1) then
      if (zt(k) < thres_z2) then
         f_loss_thres = (thres_z2 - zt(k))/(thres_z2 - thres_z1)
      else
         f_loss_thres = c0
      endif
   else
      f_loss_thres = c1
   endif

!-----------------------------------------------------------------------
!  0.001 small phytoplankton threshold C concentration (mmol C/m^3)
!  get small phyto loss(in C units)
!  small phyto agg loss
!  get grazing rate (graze_sp) on small phyto  (in C units)
!-----------------------------------------------------------------------

   C_loss_thres = f_loss_thres * loss_thres_sp

   Pprime = max(spC_loc - C_loss_thres, c0)

   sp_loss = sp_mort * Pprime * Tfunc

   sp_agg = min((sp_agg_rate_max * dps) * Pprime, sp_mort2 * Pprime * Pprime)

   reduceV = Pprime * Pprime
   graze_sp = z_umax * zooC_loc * (reduceV / (reduceV + parm_z_grz))

!-----------------------------------------------------------------------
!  routing of graze_sp & sp_loss
!  sp_agg all goes to POC
!  currently assumes that 33% of grazed caco3 is remineralized
!  if z_ingest ever changes, coefficients on routing grazed sp must change!
!  min.%C routed to POC from grazing for ballast requirements = 0.4 * Qcaco3
!  min.%C routed from sp_loss = 0.59 * QCaCO3, or P_CaCO3%rho
!-----------------------------------------------------------------------

   graze_sp_zoo = z_ingest * graze_sp
   graze_sp_poc = graze_sp * max((caco3_poc_min * QCaCO3),  &
                                 min(max(0.05_r8,(spc_poc_fac * Pprime)),&
                                     f_graze_sp_poc_lim))

   graze_sp_doc = f_graze_sp_doc * graze_sp
   graze_sp_dic = f_graze_sp_dic * graze_sp - graze_sp_poc

   sp_loss_poc = QCaCO3 * sp_loss
   sp_loss_doc = (c1 - parm_labile_ratio) * (sp_loss - sp_loss_poc)
   sp_loss_dic = parm_labile_ratio * (sp_loss - sp_loss_poc)

!-----------------------------------------------------------------------
!  0.01 small diatom threshold C concentration (mmol C/m^3)
!  get diatom loss(in C units)
!  Diatom agg loss, min. 1%/day
!  get grazing rate (graze_diat) on diatoms  (in C units)
!-----------------------------------------------------------------------

   C_loss_thres = f_loss_thres * loss_thres_diat

   Pprime = max(diatC_loc - C_loss_thres, c0)

   diat_loss = diat_mort * Pprime * Tfunc

   diat_agg = min((diat_agg_rate_max * dps) * Pprime, diat_mort2 * Pprime * Pprime)
   diat_agg = max((diat_agg_rate_min * dps) * Pprime, diat_agg)

!-----------------------------------------------------------------------
!  Lower z_grz term for diatoms and diazotrophs, larger, more mobile predators
!-----------------------------------------------------------------------

   reduceV = Pprime * Pprime
   graze_diat = diat_umax * zooC_loc * &
      (reduceV / (reduceV + 0.7_r8))

!-----------------------------------------------------------------------
!  routing of graze_diat & diat_loss
!  diat_agg all goes to POC
!  NOTE: if z_ingest is changed, coeff.s for poc,doc and dic must change!
!-----------------------------------------------------------------------

   graze_diat_zoo = z_ingest * graze_diat
   graze_diat_poc = f_graze_diat_poc * graze_diat
   graze_diat_doc = f_graze_diat_doc * graze_diat
   graze_diat_dic = f_graze_diat_dic * graze_diat

   diat_loss_poc = f_diat_loss_poc * diat_loss
   diat_loss_doc = (c1-parm_labile_ratio) * f_diat_loss_dc * diat_loss
   diat_loss_dic = parm_labile_ratio * f_diat_loss_dc * diat_loss

!-----------------------------------------------------------------------
!  0.03 small diazotroph threshold C concentration (mmol C/m^3)
!  Lower value used at temperatures < 16 deg. C, negligible biomass
!  get diazotroph loss(in C units)
!  get grazing rate (graze_diaz) on diazotrophs  (in C units)
!  no aggregation loss for diazotrophs
!-----------------------------------------------------------------------

   C_loss_diaz = f_loss_thres * loss_thres_diaz

   where (TEMP .lt. diaz_temp_thres) &
      C_loss_diaz = f_loss_thres * loss_thres_diaz2

   Pprime = max(diazC_loc - C_loss_diaz, c0)

   diaz_loss = diaz_mort * Pprime * Tfunc

   reduceV = Pprime * Pprime
   graze_diaz = diaz_umax * zooC_loc * (reduceV / (reduceV + parm_z_grz))

!-----------------------------------------------------------------------
!  routing of graze_diaz & diaz_loss
!-----------------------------------------------------------------------

   graze_diaz_zoo = f_graze_diaz_zoo * graze_diaz
   graze_diaz_poc = f_graze_diaz_poc * graze_diaz
   graze_diaz_doc = f_graze_diaz_doc * graze_diaz
   graze_diaz_dic = f_graze_diaz_dic * graze_diaz

   diaz_loss_poc = f_diaz_loss_poc * diaz_loss
   diaz_loss_doc = (c1 - parm_labile_ratio) * (diaz_loss - diaz_loss_poc)
   diaz_loss_dic = parm_labile_ratio * (diaz_loss - diaz_loss_poc)

!-----------------------------------------------------------------------
!  Note as diazotrophs have different Qp, we must route enough P into zoopl
!  and sinking detritus pools to fill their fixed p/C ratios.  The left over
!  P (remaining_diazP) is split between DOP and DIP pools
!-----------------------------------------------------------------------

   remaining_diazP =  ((graze_diaz + diaz_loss) * Qp_diaz) - &
                      ((graze_diaz_poc+graze_diaz_zoo) * Qp)
   diaz_loss_dop = (c1 - parm_labile_ratio) * remaining_diazP
   diaz_loss_dip = parm_labile_ratio * remaining_diazP

!-----------------------------------------------------------------------
!  get fractional factor for routing of zoo losses, based on food supply
!  more material is routed to large detrital pool when diatoms eaten
!-----------------------------------------------------------------------

   f_zoo_detr = (f_diat_zoo_detr * (graze_diat + epsC * epsTinv) + &
                 f_sp_zoo_detr * (graze_sp + epsC * epsTinv) + &
                 f_diaz_zoo_detr * (graze_diaz + epsC * epsTinv)) / &
                (graze_diat + graze_sp + graze_diaz + c3 * epsC * epsTinv)

!-----------------------------------------------------------------------
!  0.01 small zoo threshold C concentration (mmol C/m^3)
!  zoo losses
!-----------------------------------------------------------------------

   C_loss_thres = f_loss_thres * loss_thres_zoo

   Zprime = max(zooC_loc - C_loss_thres, c0)

   zoo_loss = z_mort2 * Zprime**1.4_r8 + z_mort * Zprime

   zoo_loss_doc = (c1 - parm_labile_ratio) * (c1 - f_zoo_detr) * zoo_loss
   zoo_loss_dic = parm_labile_ratio * (c1 - f_zoo_detr) * zoo_loss

!-----------------------------------------------------------------------
!  compute terms for DOM
!-----------------------------------------------------------------------

   DOC_prod = sp_loss_doc + graze_sp_doc + zoo_loss_doc + diat_loss_doc &
              + graze_diat_doc + diaz_loss_doc + graze_diaz_doc

   DON_prod = DOC_prod * Q
   DOP_prod = (sp_loss_doc + graze_sp_doc + zoo_loss_doc + diat_loss_doc &
               + graze_diat_doc) * Qp + diaz_loss_dop
   DOFe_prod = (zoo_loss_doc * Qfe_zoo) &
               + (Qfe_sp * (graze_sp_doc + sp_loss_doc)) &
               + (Qfe_diat * (graze_diat_doc + diat_loss_doc)) &
               + (Qfe_diaz * (graze_diaz_doc + diaz_loss_doc))

   DOC_remin = DOC_loc * DOM_remin
   DON_remin = DON_loc * DOM_remin
   DOFe_remin = DOFe_loc * DOM_remin
   DOP_remin = DOP_loc * DOM_remin

!-----------------------------------------------------------------------
!  large detritus C
!-----------------------------------------------------------------------

   POC%prod(:,:,bid) = sp_agg + graze_sp_poc + sp_loss_poc + &
                       f_zoo_detr * zoo_loss + diat_loss_poc + &
                       diat_agg + graze_diat_poc + graze_diaz_poc

!-----------------------------------------------------------------------
!  large detrital CaCO3
!  33% of CaCO3 is remin when phyto are grazed
!-----------------------------------------------------------------------

   P_CaCO3%prod(:,:,bid) = ((c1 - f_graze_CaCO3_REMIN) * graze_sp + &
                            sp_loss + sp_agg) * QCaCO3

!-----------------------------------------------------------------------
!  large detritus SiO2
!  grazed diatom SiO2, 60% is remineralized
!-----------------------------------------------------------------------

   P_SiO2%prod(:,:,bid) = ((c1 - f_graze_si_remin) * graze_diat + &
                           diat_agg  + f_diat_loss_poc * diat_loss) * Qsi

   dust%prod(:,:,bid) = c0

!-----------------------------------------------------------------------
!  Compute iron scavenging :
!  1) compute in terms of loss per year per unit iron (%/year/fe)
!  2) scale by sinking POMx6 + Dust + bSi + CaCO3 flux
!  3) increase scavenging at higher iron (>0.6nM)
!  4) convert to net loss per second
!-----------------------------------------------------------------------

   Fe_scavenge_rate = Fe_scavenge_rate0

   Fe_scavenge_rate = Fe_scavenge_rate * &
      ((POC%sflux_out(:,:,bid) + POC%hflux_out(:,:,bid)) * 72.06_r8 + &
       (P_CaCO3%sflux_out(:,:,bid) + P_CaCO3%hflux_out(:,:,bid)) * P_CaCO3%mass + &
       (P_SiO2%sflux_out(:,:,bid)+P_SiO2%hflux_out(:,:,bid))*P_SiO2%mass + &
       (dust%sflux_out(:,:,bid) + dust%hflux_out(:,:,bid)) * dust_fescav_scale)

   where (Fe_loc > Fe_scavenge_thres1) &
      Fe_scavenge_rate = Fe_scavenge_rate + &
                         (Fe_loc - Fe_scavenge_thres1) * fe_max_scale2

   Fe_scavenge = yps * Fe_loc * Fe_scavenge_rate

   P_iron%prod(:,:,bid) = &
      ((sp_agg + graze_sp_poc + sp_loss_poc) * Qfe_sp) &
      + (zoo_loss * f_zoo_detr * Qfe_zoo) &
      + ((diat_agg + graze_diat_poc + diat_loss_poc) * Qfe_diat) &
      + (graze_diaz_poc * Qfe_diaz) + (f_fescav_P_iron * Fe_scavenge)

   call compute_particulate_terms(k, POC, P_CaCO3, P_SiO2, dust, P_iron, &
                                  QA_dust_def, TEMP, O2_loc, this_block)

!-----------------------------------------------------------------------
!  nitrate & ammonium
!  nitrification in low light
!  use exponential decay of PAR across model level to compute taper factor
!-----------------------------------------------------------------------

   if (lrest_no3) then
      if (lnutr_variable_restore) then
         RESTORE = NUTR_RESTORE_RTAU(:,:,bid) * &
                      merge((NO3_CLIM(:,:,k,bid) - NO3_loc), &
                            c0, k <= NUTR_RESTORE_MAX_LEVEL(:,:,bid))
      else
         RESTORE = (NO3_CLIM(:,:,k,bid) - NO3_loc) * nutr_rest_time_inv(k)
      endif
   else
      RESTORE = c0
   endif

   call accumulate_tavg_field(RESTORE, tavg_NO3_RESTORE,bid,k)

   where (PAR_out(:,:,bid) < parm_nitrif_par_lim)
      NITRIF = parm_kappa_nitrif * NH4_loc
      where (PAR_in > parm_nitrif_par_lim)
         NITRIF = NITRIF * log(PAR_out(:,:,bid) / parm_nitrif_par_lim) / (-KPARdz)
      end where
   elsewhere
      NITRIF = c0
   end where

   call accumulate_tavg_field(NITRIF, tavg_NITRIF,bid,k)

!-----------------------------------------------------------------------
!  Compute denitrification under low O2 conditions
!-----------------------------------------------------------------------

   WORK1 = ((parm_o2_min+parm_o2_min_delta) - O2_loc) / parm_o2_min_delta
   WORK1 = min(max(WORK1,c0),c1)
   where (NO3_loc .lt. parm_no3_min)
      WORK1 = WORK1 * (NO3_loc / parm_no3_min)
   end where
   DENITRIF = WORK1 * (DOC_remin + POC%remin(:,:,bid)) / denitrif_C_N

   call accumulate_tavg_field(DENITRIF, tavg_DENITRIF,bid,k)

!-----------------------------------------------------------------------
!  nitrate & ammonium
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,no3_ind) = RESTORE + NITRIF - NO3_V_diat - &
      NO3_V_sp - DENITRIF - photoNO3_diaz

   DTRACER_MODULE(:,:,nh4_ind) = -NH4_V_diat - NH4_V_sp - NITRIF + &
        Q * (zoo_loss_dic + sp_loss_dic + graze_sp_dic + diat_loss_dic + &
        graze_diat_dic + POC%remin(:,:,bid) + diaz_loss_dic + &
        graze_diaz_dic) + DON_remin + diaz_Nexcrete - photoNH4_diaz

!-----------------------------------------------------------------------
!  dissolved iron
!-----------------------------------------------------------------------

     DTRACER_MODULE(:,:,fe_ind) = P_iron%remin(:,:,bid)  &
       + (Qfe_zoo * zoo_loss_dic) + DOFe_remin  &
       + (Qfe_sp * (sp_loss_dic + graze_sp_dic))     &
       + (Qfe_diat * (diat_loss_dic + graze_diat_dic))     &
       + (Qfe_diaz * (diaz_loss_dic + graze_diaz_dic))     &
       + graze_diaz_zoo * (Qfe_diaz-Qfe_zoo)     &
       + graze_diat_zoo * (Qfe_diat-Qfe_zoo)     &
       + graze_sp_zoo * (Qfe_sp-Qfe_zoo)     &
       - photoFe_sp - photoFe_diat - photoFe_diaz - Fe_scavenge

!-----------------------------------------------------------------------
!  dissolved SiO3
!-----------------------------------------------------------------------

   if (lrest_sio3) then
      if (lnutr_variable_restore) then
         RESTORE = NUTR_RESTORE_RTAU(:,:,bid) * &
                      merge((SiO3_CLIM(:,:,k,bid) - SiO3_loc), &
                            c0, k <= NUTR_RESTORE_MAX_LEVEL(:,:,bid))
      else
         RESTORE = (SiO3_CLIM(:,:,k,bid) - SiO3_loc) * nutr_rest_time_inv(k)
      endif
   else
      RESTORE = c0
   endif

   call accumulate_tavg_field(RESTORE, tavg_SiO3_RESTORE,bid,k)

   DTRACER_MODULE(:,:,sio3_ind) = RESTORE + P_SiO2%remin(:,:,bid) + &
      Qsi * (f_graze_si_remin * graze_diat +f_diat_loss_dc * diat_loss) &
      - photoSi_diat

!-----------------------------------------------------------------------
!  phosphate
!-----------------------------------------------------------------------

   if (lrest_po4) then
      if (lnutr_variable_restore) then
         RESTORE = NUTR_RESTORE_RTAU(:,:,bid) * &
                      merge((PO4_CLIM(:,:,k,bid) - PO4_loc), &
                            c0, k <= NUTR_RESTORE_MAX_LEVEL(:,:,bid))
      else
         RESTORE = (PO4_CLIM(:,:,k,bid) - PO4_loc) * nutr_rest_time_inv(k)
      endif
   else
      RESTORE = c0
   endif

   call accumulate_tavg_field(RESTORE, tavg_PO4_RESTORE,bid,k)

   DTRACER_MODULE(:,:,po4_ind) = RESTORE + (Qp * (POC%remin(:,:,bid) + &
      zoo_loss_dic + sp_loss_dic + graze_sp_dic + diat_loss_dic + &
      graze_diat_dic - photoC_diat)) &
      + DOP_remin + diaz_loss_dip - po4_V_sp - po4_V_diaz

!-----------------------------------------------------------------------
!  small phyto Carbon
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,spC_ind) = photoC_sp - graze_sp - sp_loss - sp_agg

!-----------------------------------------------------------------------
!  small phyto Chlorophyll
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,spChl_ind) = photoacc_sp - &
      thetaC_sp * (graze_sp + sp_loss + sp_agg)

!-----------------------------------------------------------------------
!  small phytoplankton CaCO3
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,spCaCO3_ind) = CaCO3_PROD - &
      (graze_sp + sp_loss + sp_agg) * QCaCO3

!-----------------------------------------------------------------------
!  diatom Carbon
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,diatC_ind) = photoC_diat - graze_diat - &
      diat_loss - diat_agg

!-----------------------------------------------------------------------
!  diatom Chlorophyll
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,diatChl_ind) = photoacc_diat - &
      thetaC_diat * (graze_diat + diat_loss + diat_agg)

!-----------------------------------------------------------------------
!  zoo Carbon
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,zooC_ind) = graze_sp_zoo + graze_diat_zoo &
      + graze_diaz_zoo - zoo_loss

!-----------------------------------------------------------------------
!  dissolved organic Matter
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,doc_ind) = DOC_prod - DOC_remin

   DTRACER_MODULE(:,:,don_ind) = DON_prod - DON_remin

   DTRACER_MODULE(:,:,dop_ind) = DOP_prod - DOP_remin

   DTRACER_MODULE(:,:,dofe_ind) = DOFe_prod - DOFe_remin

!-----------------------------------------------------------------------
!  small phyto Fe
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,spFe_ind) =  photoFe_sp &
      - (Qfe_sp * (graze_sp+sp_loss+sp_agg))

!-----------------------------------------------------------------------
!  Diatom Fe
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,diatFe_ind) =  photoFe_diat &
      - (Qfe_diat * (graze_diat+diat_loss+diat_agg))

!-----------------------------------------------------------------------
!  Diatom Si
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,diatSi_ind) =  photoSi_diat &
      - (Qsi * (graze_diat+diat_loss+diat_agg))

!-----------------------------------------------------------------------
!  Diazotroph C
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,diazC_ind) =  photoC_diaz - graze_diaz - diaz_loss

!-----------------------------------------------------------------------
!  diazotroph Chlorophyll
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,diazChl_ind) = photoacc_diaz - &
      thetaC_diaz * (graze_diaz + diaz_loss)

!-----------------------------------------------------------------------
!  Diazotroph Fe
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,diazFe_ind) =  photoFe_diaz &
      - (Qfe_diaz * (graze_diaz + diaz_loss))

!-----------------------------------------------------------------------
!   dissolved inorganic Carbon
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,dic_ind) = &
      DOC_remin + POC%remin(:,:,bid) + P_CaCO3%remin(:,:,bid) + &
      f_graze_CaCO3_REMIN * graze_sp * QCaCO3 + zoo_loss_dic + sp_loss_dic + &
      graze_sp_dic + diat_loss_dic + graze_diat_dic - photoC_sp - photoC_diat &
      - CaCO3_PROD + graze_diaz_dic + diaz_loss_dic - photoC_diaz

!-----------------------------------------------------------------------
!  alkalinity
!-----------------------------------------------------------------------

   DTRACER_MODULE(:,:,alk_ind) = -DTRACER_MODULE(:,:,no3_ind) + &
      DTRACER_MODULE(:,:,nh4_ind) + &
      c2 * (P_CaCO3%remin(:,:,bid) + f_graze_CaCO3_REMIN * graze_sp * QCaCO3 &
            - CaCO3_PROD)

!-----------------------------------------------------------------------
!  oxygen
!-----------------------------------------------------------------------

!  DTRACER_MODULE(:,:,o2_ind) = (photoC_sp + photoC_diat) / &
!     parm_Red_D_C_O2 + photoC_diaz/parm_Red_D_C_O2_diaz

!  WORK1 = (O2_loc - parm_o2_min) / parm_o2_min_delta
!  WORK1 = min(max(WORK1,c0),c1)
!  DTRACER_MODULE(:,:,o2_ind) = DTRACER_MODULE(:,:,o2_ind) + WORK1 &
!     * ((- POC%remin(:,:,bid) - DOC_remin - zoo_loss_dic - sp_loss_dic  &
!        - graze_sp_dic - diat_loss_dic - graze_diat_dic   &
!        - graze_diaz_dic - diaz_loss_dic) / parm_Red_P_C_O2)

   DTRACER_MODULE(:,:,o2_ind) = c0

   O2_PRODUCTION = c0

   where (photoC_diat > c0)
      O2_PRODUCTION = O2_PRODUCTION + photoC_diat * &
         ((NO3_V_diat/(NO3_V_diat+NH4_V_diat)) / parm_Red_D_C_O2 + &
          (NH4_V_diat/(NO3_V_diat+NH4_V_diat)) / parm_Remin_D_C_O2)
   end where

   where (photoC_sp > c0)
      O2_PRODUCTION = O2_PRODUCTION + photoC_sp * &
         ((NO3_V_sp/(NO3_V_sp+NH4_V_sp)) / parm_Red_D_C_O2 + &
          (NH4_V_sp/(NO3_V_sp+NH4_V_sp)) / parm_Remin_D_C_O2)
   end where

   where (photoC_diaz > c0)
      O2_PRODUCTION = O2_PRODUCTION + photoC_diaz * &
         ((photoNO3_diaz/(photoNO3_diaz+photoNH4_diaz+diaz_Nfix)) / parm_Red_D_C_O2 + &
          (photoNH4_diaz/(photoNO3_diaz+photoNH4_diaz+diaz_Nfix)) / parm_Remin_D_C_O2 + &
          (diaz_Nfix/(photoNO3_diaz+photoNH4_diaz+diaz_Nfix)) / parm_Red_D_C_O2_diaz)
   end where

   WORK1 = (O2_loc - parm_o2_min) / parm_o2_min_delta
   WORK1 = min(max(WORK1,c0),c1)
   O2_CONSUMPTION = WORK1 * &
      ((POC%remin(:,:,bid) + DOC_remin + zoo_loss_dic + &
        sp_loss_dic + graze_sp_dic + diat_loss_dic + graze_diat_dic + &
        graze_diaz_dic + diaz_loss_dic)/ parm_Remin_D_C_O2 + (c2*NITRIF))

   DTRACER_MODULE(:,:,o2_ind) = O2_PRODUCTION - O2_CONSUMPTION

!-----------------------------------------------------------------------
!  various tavg/history variables
!-----------------------------------------------------------------------

   if (k == 1) then
      if (accumulate_tavg_now(tavg_O2_ZMIN) .or. accumulate_tavg_now(tavg_O2_ZMIN_DEPTH)) then
         ! WORK1 = O2 at this level
         ! WORK2 = vertical min of O2
         ! WORK3 = depth of min

         kk = 1
         WORK1 = p5*(TRACER_MODULE_OLD(:,:,kk,o2_ind) + &
                     TRACER_MODULE_CUR(:,:,kk,o2_ind))
         WORK2 = WORK1
         WORK3 = zt(kk)

         do kk = 2,km
            WORK1 = p5*(TRACER_MODULE_OLD(:,:,kk,o2_ind) + &
                        TRACER_MODULE_CUR(:,:,kk,o2_ind))
            where (kk <= KMT(:,:,bid) .and. (WORK1 < WORK2))
               WORK2 = WORK1
               WORK3 = zt(kk)
            endwhere
         end do

         call accumulate_tavg_field(WORK2, tavg_O2_ZMIN,bid,k)

         call accumulate_tavg_field(WORK3, tavg_O2_ZMIN_DEPTH,bid,k)
      endif
   endif

   call accumulate_tavg_field(O2_PRODUCTION, tavg_O2_PRODUCTION,bid,k)

   call accumulate_tavg_field(O2_CONSUMPTION, tavg_O2_CONSUMPTION,bid,k)

   if (accumulate_tavg_now(tavg_AOU)) then
      WORK1 = O2SAT(TEMP, SALT, &
                    LAND_MASK(:,:,bid) .and. (k .le. KMT(:,:,bid)))
      WORK1 = WORK1 - p5*(TRACER_MODULE_OLD(:,:,k,o2_ind) + &
                          TRACER_MODULE_CUR(:,:,k,o2_ind))
      call accumulate_tavg_field(WORK1, tavg_AOU,bid,k)
   endif

   call accumulate_tavg_field(PAR_avg, tavg_PAR_avg,bid,k)

   call accumulate_tavg_field(graze_sp, tavg_graze_sp,bid,k)

   call accumulate_tavg_field(graze_diat, tavg_graze_diat,bid,k)

   call accumulate_tavg_field(graze_diaz, tavg_graze_diaz,bid,k)

   if (accumulate_tavg_now(tavg_graze_TOT)) then
      WORK1 = graze_sp + graze_diat + graze_diaz
      call accumulate_tavg_field(WORK1, tavg_graze_TOT,bid,k)
   endif

   call accumulate_tavg_field(sp_loss, tavg_sp_loss,bid,k)

   call accumulate_tavg_field(diat_loss, tavg_diat_loss,bid,k)

   call accumulate_tavg_field(diaz_loss, tavg_diaz_loss,bid,k)

   call accumulate_tavg_field(zoo_loss, tavg_zoo_loss,bid,k)

   call accumulate_tavg_field(sp_agg, tavg_sp_agg,bid,k)

   call accumulate_tavg_field(diat_agg, tavg_diat_agg,bid,k)

   call accumulate_tavg_field(photoC_sp, tavg_photoC_sp,bid,k)

   call accumulate_tavg_field(photoC_diat, tavg_photoC_diat,bid,k)

   call accumulate_tavg_field(photoC_diaz, tavg_photoC_diaz,bid,k)

   if (accumulate_tavg_now(tavg_photoC_TOT)) then
      WORK1 = photoC_sp + photoC_diat + photoC_diaz
      call accumulate_tavg_field(WORK1, tavg_photoC_TOT,bid,k)
   endif

   if (accumulate_tavg_now(tavg_photoC_sp_zint)) then
      if (partial_bottom_cells) then
         WORK1 = DZT(:,:,k,bid) * photoC_sp
      else
         WORK1 = dz(k) * photoC_sp
      endif
      call accumulate_tavg_field(WORK1, tavg_photoC_sp_zint,bid,k)
   endif

   if (accumulate_tavg_now(tavg_photoC_diat_zint)) then
      if (partial_bottom_cells) then
         WORK1 = DZT(:,:,k,bid) * photoC_diat
      else
         WORK1 = dz(k) * photoC_diat
      endif
      call accumulate_tavg_field(WORK1, tavg_photoC_diat_zint,bid,k)
   endif

   if (accumulate_tavg_now(tavg_photoC_diaz_zint)) then
      if (partial_bottom_cells) then
         WORK1 = DZT(:,:,k,bid) * photoC_diaz
      else
         WORK1 = dz(k) * photoC_diaz
      endif
      call accumulate_tavg_field(WORK1, tavg_photoC_diaz_zint,bid,k)
   endif

   if (accumulate_tavg_now(tavg_photoC_TOT_zint)) then
      if (partial_bottom_cells) then
         WORK1 = DZT(:,:,k,bid) * (photoC_sp + photoC_diat + photoC_diaz)
      else
         WORK1 = dz(k) * (photoC_sp + photoC_diat + photoC_diaz)
      endif
      call accumulate_tavg_field(WORK1, tavg_photoC_TOT_zint,bid,k)
   endif

   if (accumulate_tavg_now(tavg_photoC_NO3_sp) .or. &
       accumulate_tavg_now(tavg_photoC_NO3_sp_zint)) then
      where (VNtot_sp > c0)
         WORK1 = (VNO3_sp / VNtot_sp) * photoC_sp
      elsewhere
         WORK1 = c0
      end where

      call accumulate_tavg_field(WORK1, tavg_photoC_NO3_sp,bid,k)

      if (accumulate_tavg_now(tavg_photoC_NO3_sp_zint)) then
         if (partial_bottom_cells) then
            WORK1 = DZT(:,:,k,bid) * WORK1
         else
            WORK1 = dz(k) * WORK1
         endif
         call accumulate_tavg_field(WORK1, tavg_photoC_NO3_sp_zint,bid,k)
      endif
   endif

   if (accumulate_tavg_now(tavg_photoC_NO3_diat) .or. &
       accumulate_tavg_now(tavg_photoC_NO3_diat_zint)) then
      where (VNtot_diat > c0)
         WORK1 = (VNO3_diat / VNtot_diat) * photoC_diat
      elsewhere
         WORK1 = c0
      end where

      call accumulate_tavg_field(WORK1, tavg_photoC_NO3_diat,bid,k)

      if (accumulate_tavg_now(tavg_photoC_NO3_diat_zint)) then
         if (partial_bottom_cells) then
            WORK1 = DZT(:,:,k,bid) * WORK1
         else
            WORK1 = dz(k) * WORK1
         endif
         call accumulate_tavg_field(WORK1, tavg_photoC_NO3_diat_zint,bid,k)
      endif
   endif

   if (accumulate_tavg_now(tavg_photoC_NO3_diaz) .or. &
       accumulate_tavg_now(tavg_photoC_NO3_diaz_zint)) then
      where (VNtot_diaz > c0)
         WORK1 = (VNO3_diaz / VNtot_diaz) * photoC_diaz
      elsewhere
         WORK1 = c0
      end where

      call accumulate_tavg_field(WORK1, tavg_photoC_NO3_diaz,bid,k)

      if (accumulate_tavg_now(tavg_photoC_NO3_diaz_zint)) then
         if (partial_bottom_cells) then
            WORK1 = DZT(:,:,k,bid) * WORK1
         else
            WORK1 = dz(k) * WORK1
         endif
         call accumulate_tavg_field(WORK1, tavg_photoC_NO3_diaz_zint,bid,k)
      endif
   endif

   if (accumulate_tavg_now(tavg_photoC_NO3_TOT) .or. &
       accumulate_tavg_now(tavg_photoC_NO3_TOT_zint)) then
      where (VNtot_sp > c0)
         WORK1 = WORK1 + (VNO3_sp / VNtot_sp) * photoC_sp
      elsewhere
         WORK1 = c0
      end where
      where (VNtot_diat > c0)
         WORK1 = WORK1 + (VNO3_diat / VNtot_diat) * photoC_diat
      end where
      where (VNtot_diaz > c0)
         WORK1 = WORK1 + (VNO3_diaz / VNtot_diaz) * photoC_diaz
      end where

      call accumulate_tavg_field(WORK1, tavg_photoC_NO3_TOT,bid,k)

      if (accumulate_tavg_now(tavg_photoC_NO3_TOT_zint)) then
         if (partial_bottom_cells) then
            WORK1 = DZT(:,:,k,bid) * WORK1
         else
            WORK1 = dz(k) * WORK1
         endif
         call accumulate_tavg_field(WORK1, tavg_photoC_NO3_TOT_zint,bid,k)
      endif
   endif

   call accumulate_tavg_field(DOC_prod, tavg_DOC_prod,bid,k)

   call accumulate_tavg_field(DOC_remin, tavg_DOC_remin,bid,k)

   call accumulate_tavg_field(DON_prod, tavg_DON_prod,bid,k)

   call accumulate_tavg_field(DON_remin, tavg_DON_remin,bid,k)

   call accumulate_tavg_field(DOP_prod, tavg_DOP_prod,bid,k)

   call accumulate_tavg_field(DOP_remin, tavg_DOP_remin,bid,k)

   call accumulate_tavg_field(DOFe_prod, tavg_DOFe_prod,bid,k)

   call accumulate_tavg_field(DOFe_remin, tavg_DOFe_remin,bid,k)

   call accumulate_tavg_field(Fe_scavenge, tavg_Fe_scavenge,bid,k)

   call accumulate_tavg_field(Fe_scavenge_rate, tavg_Fe_scavenge_rate,bid,k)

   if (accumulate_tavg_now(tavg_CO3) .or. &
       accumulate_tavg_now(tavg_HCO3) .or. &
       accumulate_tavg_now(tavg_H2CO3) .or. &
       accumulate_tavg_now(tavg_pH_3D) .or. &
       accumulate_tavg_now(tavg_zsatcalc) .or. &
       accumulate_tavg_now(tavg_zsatarag)) then

      where (PH_PREV_3D(:,:,k,bid) /= c0)
         WORK1 = PH_PREV_3D(:,:,k,bid) - del_ph
         WORK2 = PH_PREV_3D(:,:,k,bid) + del_ph
      elsewhere
         WORK1 = phlo_3d_init
         WORK2 = phhi_3d_init
      end where

      call timer_start(ecosys_comp_CO3terms_timer, block_id=bid)
      call comp_CO3terms(bid, k, LAND_MASK(:,:,bid) .and. k <= KMT(:,:,bid), &
                         TEMP, SALT, DIC_loc, ALK_loc, PO4_loc, SiO3_loc, &
                         WORK1, WORK2, WORK3, H2CO3, HCO3, CO3)
      call timer_stop(ecosys_comp_CO3terms_timer, block_id=bid)
      call accumulate_tavg_field(CO3, tavg_CO3,bid,k)
      call accumulate_tavg_field(HCO3, tavg_HCO3,bid,k)
      call accumulate_tavg_field(H2CO3, tavg_H2CO3,bid,k)
      call accumulate_tavg_field(WORK3, tavg_pH_3D,bid,k)

      PH_PREV_3D(:,:,k,bid) = WORK3
   endif

   if (accumulate_tavg_now(tavg_co3_sat_calc) .or. &
       accumulate_tavg_now(tavg_zsatcalc) .or. &
       accumulate_tavg_now(tavg_co3_sat_arag) .or. &
       accumulate_tavg_now(tavg_zsatarag)) then
      call comp_co3_sat_vals(k, LAND_MASK(:,:,bid) .and. k <= KMT(:,:,bid), &
                             TEMP, SALT, WORK1, WORK2)
      call accumulate_tavg_field(WORK1, tavg_co3_sat_calc,bid,k)
      if (accumulate_tavg_now(tavg_zsatcalc)) then
         if (k == 1) then
            ! set to -1, i.e. depth not found yet,
            ! if mask == .true. and surface supersaturated to -1
            ZSATCALC(:,:,bid) = merge(-c1, c0, LAND_MASK(:,:,bid) .and. CO3 > WORK1)
         else
            where (ZSATCALC(:,:,bid) == -c1 .and. CO3 <= WORK1)
               ZSATCALC(:,:,bid) = zt(k-1) + (zt(k)-zt(k-1)) * &
                  CO3_CALC_ANOM_km1(:,:,bid) / (CO3_CALC_ANOM_km1(:,:,bid) - (CO3 - WORK1))
            endwhere
            where (ZSATCALC(:,:,bid) == -c1 .and. KMT(:,:,bid) == k)
               ZSATCALC(:,:,bid) = zw(k)
            endwhere
         endif
         CO3_CALC_ANOM_km1(:,:,bid) = CO3 - WORK1
         if (k == km) then
            call accumulate_tavg_field(ZSATCALC(:,:,bid), tavg_zsatcalc,bid,k)
         endif
      endif
      call accumulate_tavg_field(WORK2, tavg_co3_sat_arag,bid,k)
      if (accumulate_tavg_now(tavg_zsatarag)) then
         if (k == 1) then
            ! set to -1, i.e. depth not found yet,
            ! if mask == .true. and surface supersaturated to -1
            ZSATARAG(:,:,bid) = merge(-c1, c0, LAND_MASK(:,:,bid) .and. CO3 > WORK2)
         else
            where (ZSATARAG(:,:,bid) == -c1 .and. CO3 <= WORK2)
               ZSATARAG(:,:,bid) = zt(k-1) + (zt(k)-zt(k-1)) * &
                  CO3_ARAG_ANOM_km1(:,:,bid) / (CO3_ARAG_ANOM_km1(:,:,bid) - (CO3 - WORK2))
            endwhere
            where (ZSATARAG(:,:,bid) == -c1 .and. KMT(:,:,bid) == k)
               ZSATARAG(:,:,bid) = zw(k)
            endwhere
         endif
         CO3_ARAG_ANOM_km1(:,:,bid) = CO3 - WORK2
         if (k == km) then
            call accumulate_tavg_field(ZSATARAG(:,:,bid), tavg_zsatarag,bid,k)
         endif
      endif
   endif

   call timer_stop(ecosys_interior_timer, block_id=bid)

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_set_interior

!***********************************************************************
!BOP
! !IROUTINE: init_particulate_terms
! !INTERFACE:

 subroutine init_particulate_terms(POC, P_CaCO3, P_SiO2, dust, P_iron, &
       QA_dust_def, NET_DUST_IN, this_block)

! !DESCRIPTION:
!  Set incoming fluxes (put into outgoing flux for first level usage).
!  Set dissolution length, production fraction and mass terms.
!
!  The first 6 arguments are intent(inout) in
!  order to preserve contents on other blocks.

! !INPUT/OUTPUT PARAMETERS:

   type(sinking_particle), intent(inout) :: &
      POC,          & ! base units = nmol C
      P_CaCO3,      & ! base units = nmol CaCO3
      P_SiO2,       & ! base units = nmol SiO2
      dust,         & ! base units = g
      P_iron          ! base units = nmol Fe

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(inout) :: &
      QA_dust_def     ! incoming deficit in the QA(dust) POC flux

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      NET_DUST_IN     ! dust flux

   type (block), intent(in) :: &
      this_block      ! block info for the current block

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bid                 ! local_block id

   character (char_len) :: &
      tracer_data_label   ! label for what is being updated

   character (char_len), dimension(:), allocatable :: &
      tracer_data_names   ! short names for input data fields

   integer (int_kind), dimension(:), allocatable :: &
      tracer_bndy_loc,  & ! location and field type for ghost
      tracer_bndy_type    !    cell updates

!-----------------------------------------------------------------------
!  parameters, from Armstrong et al. 2000
!
!  July 2002, length scale for excess POC and bSI modified by temperature
!  Value given here is at Tref of 30 deg. C, JKM
!-----------------------------------------------------------------------

    POC%diss      = 21000.0_r8  ! diss. length (cm), modified by TEMP
    POC%gamma     = c0          ! not used
    POC%mass      = 12.01_r8    ! molecular weight of POC
    POC%rho       = c0          ! not used

    P_CaCO3%diss  = 80000.0_r8  ! diss. length (cm)
    P_CaCO3%gamma = 0.55_r8     ! prod frac -> hard subclass
    P_CaCO3%mass  = 100.09_r8   ! molecular weight of CaCO
    P_CaCO3%rho   = 0.05_r8 * P_CaCO3%mass / POC%mass
                                       ! QA mass ratio for CaCO3
                                       ! This ratio is used in ecos_set_interior

    P_SiO2%diss   = 21000.0_r8  ! diss. length (cm), modified by TEMP
    P_SiO2%gamma  = 0.25_r8     ! prod frac -> hard subclass
    P_SiO2%mass   = 60.08_r8    ! molecular weight of SiO2
    P_SiO2%rho    = 0.05_r8 * P_SiO2%mass / POC%mass
                                       ! QA mass ratio for SiO2

    dust%diss     = 60000.0_r8  ! diss. length (cm)
    dust%gamma    = 0.97_r8     ! prod frac -> hard subclass
    dust%mass     = 1.0e9_r8    ! base units are already grams
    dust%rho      = 0.05_r8 * dust%mass / POC%mass
                                       ! QA mass ratio for dust

    P_iron%diss   = 60000.0_r8  ! diss. length (cm) - not used
    P_iron%gamma  = c0          ! prod frac -> hard subclass - not used
    P_iron%mass   = c0          ! not used
    P_iron%rho    = c0          ! not used

!-----------------------------------------------------------------------
!  Set incoming fluxes
!-----------------------------------------------------------------------

    bid = this_block%local_id

    P_CaCO3%sflux_out(:,:,bid) = c0
    P_CaCO3%hflux_out(:,:,bid) = c0

    P_SiO2%sflux_out(:,:,bid) = c0
    P_SiO2%hflux_out(:,:,bid) = c0

    if (dust_flux%has_data) then
       dust%sflux_out(:,:,bid) = (c1 - dust%gamma) * NET_DUST_IN(:,:,bid)
       dust%hflux_out(:,:,bid) = dust%gamma * NET_DUST_IN(:,:,bid)
    else
       dust%sflux_out(:,:,bid) = c0
       dust%hflux_out(:,:,bid) = c0
    endif

    P_iron%sflux_out(:,:,bid) = c0
    P_iron%hflux_out(:,:,bid) = c0

!-----------------------------------------------------------------------
!  Hard POC is QA flux and soft POC is excess POC.
!-----------------------------------------------------------------------

    POC%sflux_out(:,:,bid) = c0
    POC%hflux_out(:,:,bid) = c0

!-----------------------------------------------------------------------
!  Compute initial QA(dust) POC flux deficit.
!-----------------------------------------------------------------------

    QA_dust_def(:,:,bid) = dust%rho * &
       (dust%sflux_out(:,:,bid) + dust%hflux_out(:,:,bid))

!-----------------------------------------------------------------------
!EOC

 end subroutine init_particulate_terms

!***********************************************************************
!BOP
! !IROUTINE: compute_particulate_terms
! !INTERFACE:

 subroutine compute_particulate_terms(k, POC, P_CaCO3, P_SiO2, dust, P_iron, &
       QA_dust_def, TEMP, O2_loc, this_block)

! !DESCRIPTION:
!  Compute outgoing fluxes and remineralization terms. Assumes that
!  production terms have been set. Incoming fluxes are assumed to be the
!  outgoing fluxes from the previous level.
!
!  It is assumed that there is no production of dust.
!
!  Instantaneous remineralization in the bottom cell is implemented by
!  setting the outgoing flux to zero.
!
!  For POC, the hard subclass is the POC flux qualitatively associated
!  with the ballast flux. The soft subclass is the excess POC flux.
!
!  Remineralization for the non-iron particulate pools is computing
!  by first computing the outgoing flux and then computing the
!  remineralization from conservation, i.e.
!     flux_in - flux_out + prod * dz - remin * dz == 0.
!
!  For iron, remineralization is first computed from POC remineralization
!  and then flux_out is computed from conservation. If the resulting
!  flux_out is negative or should be zero because of the sea floor, the
!  remineralization is adjusted.
!  Note: all the sinking iron is in the P_iron%sflux pool, hflux Fe not
!        explicitly tracked, it is assumed that total iron remin is
!        proportional to total POC remin.
!
!  Based upon Armstrong et al. 2000
!
!  July 2002, added temperature effect on remin length scale of
!  excess POC (all soft POM& Iron) and on SiO2.
!  new variable passed into ballast, Tfunc, main Temperature function
!  computed in ecosystem routine.  scaling factor for dissolution
!  of excess POC, Fe, and Bsi now varies with location (f(temperature)).
!
!  Added diffusive iron flux from sediments at depths < 1100m,
!  based on Johnson et al., 1999, value of 5 umolFe/m2/day,
!      this value too high, using 2 umolFe/m2/day here
!
!  Allow hard fraction of ballast to remin with long length scale 40,000m
!     thus ~ 10% of hard ballast remins over 4000m water column.
!
!  Sinking dust flux is decreased by assumed instant solubility/dissolution
!     at ocean surface from the parm_Fe_bioavail.
!
!  Modified to allow different Q10 factors for soft POM and bSI remin,
!  water TEMP is now passed in instead of Tfunc (1/2005, JKM)

! !USES:

#ifdef CCSMCOUPLED
   use shr_sys_mod, only: shr_sys_abort
#endif

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k ! vertical model level

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      TEMP,         & ! temperature for scaling functions
      O2_loc          ! dissolved oxygen used to change/increase POC%diss

   type (block), intent(in) :: &
      this_block      ! block info for the current block

! !INPUT/OUTPUT PARAMETERS:

   type(sinking_particle), intent(inout) :: &
      POC,          & ! base units = nmol C
      P_CaCO3,      & ! base units = nmol CaCO3
      P_SiO2,       & ! base units = nmol SiO2
      dust,         & ! base units = g
      P_iron          ! base units = nmol Fe

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(inout) :: &
      QA_dust_def     ! incoming deficit in the QA(dust) POC flux

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8) :: poc_diss ! diss. length used (cm)

   character(*), parameter :: &
      subname = 'ecosys_mod:compute_particulate_terms'

   real (r8), dimension(nx_block,ny_block) :: &
      WORK,               & ! temporary for summed quantities to be averaged
      TfuncP,             & ! temperature scaling from soft POM remin
      TfuncS,             & ! temperature scaling from soft POM remin
      DECAY_CaCO3,        & ! scaling factor for dissolution of CaCO3
      DECAY_dust,         & ! scaling factor for dissolution of dust
      DECAY_Hard,         & ! scaling factor for dissolution of Hard Ballast
      DECAY_HardDust        ! scaling factor for dissolution of Hard dust

   real (r8) :: &
      decay_POC_E,        & ! scaling factor for dissolution of excess POC
      decay_SiO2,         & ! scaling factor for dissolution of SiO2
      POC_PROD_avail,     & ! POC production available for excess POC flux
      new_QA_dust_def,    & ! outgoing deficit in the QA(dust) POC flux
      dz_loc, dzr_loc       ! dz, dzr at a particular i,j location

   integer (int_kind) :: &
      i, j,               & ! loop indices
      bid                   ! local_block id

   logical (log_kind) :: &
      poc_error             ! POC error flag

!-----------------------------------------------------------------------
!  incoming fluxes are outgoing fluxes from previous level
!-----------------------------------------------------------------------

   bid = this_block%local_id

   P_CaCO3%sflux_in(:,:,bid) = P_CaCO3%sflux_out(:,:,bid)
   P_CaCO3%hflux_in(:,:,bid) = P_CaCO3%hflux_out(:,:,bid)

   P_SiO2%sflux_in(:,:,bid) = P_SiO2%sflux_out(:,:,bid)
   P_SiO2%hflux_in(:,:,bid) = P_SiO2%hflux_out(:,:,bid)

   dust%sflux_in(:,:,bid) = dust%sflux_out(:,:,bid)
   dust%hflux_in(:,:,bid) = dust%hflux_out(:,:,bid)

   POC%sflux_in(:,:,bid) = POC%sflux_out(:,:,bid)
   POC%hflux_in(:,:,bid) = POC%hflux_out(:,:,bid)

   P_iron%sflux_in(:,:,bid) = P_iron%sflux_out(:,:,bid)
   P_iron%hflux_in(:,:,bid) = P_iron%hflux_out(:,:,bid)

!-----------------------------------------------------------------------
!  compute decay factors
!-----------------------------------------------------------------------

   if (partial_bottom_cells) then
      DECAY_CaCO3    = exp(-DZT(:,:,k,bid) / P_CaCO3%diss)
      DECAY_dust     = exp(-DZT(:,:,k,bid) / dust%diss)
      DECAY_Hard     = exp(-DZT(:,:,k,bid) / 4.0e6_r8)
      DECAY_HardDust = exp(-DZT(:,:,k,bid) / 1.2e7_r8)
   else
      DECAY_CaCO3    = exp(-dz(k) / P_CaCO3%diss)
      DECAY_dust     = exp(-dz(k) / dust%diss)
      DECAY_Hard     = exp(-dz(k) / 4.0e6_r8)
      DECAY_HardDust = exp(-dz(k) / 1.2e7_r8)
   endif

!----------------------------------------------------------------------
!   Tref = 30.0 reference temperature (deg. C)
!
!   Using q10 formulation with Q10 value of 1.1 soft POM (TfuncP) and
!       a Q10 value of 2.5 soft bSi (TfuncS)
!-----------------------------------------------------------------------

!
!  NOTE: Turning off temperature effect on POM lengthscale, three instances
!    of TfuncP below have been removed, see comment lines.
!-------------------------------------------------------------------------
!  TfuncP = 1.1_r8**(((TEMP + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)

   TfuncS = 2.5_r8**(((TEMP + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)

   poc_error = .false.

   do j = 1,ny_block
      do i = 1,nx_block

         if (LAND_MASK(i,j,bid) .and. k <= KMT(i,j,bid)) then

            if (partial_bottom_cells) then
               dz_loc = DZT(i,j,k,bid)
            else
               dz_loc = dz(k)
            endif
            dzr_loc = c1 / dz_loc

!-----------------------------------------------------------------------
!  increase POC diss where O2 concentrations are low
!-----------------------------------------------------------------------

!            if (O2_loc(i,j) > c0) then
!               poc_diss = min((POC%diss * 2.0_r8), (POC%diss / &
!                  (O2_loc(i,j) / (O2_loc(i,j) + 5.0_r8))))
!            else
!               poc_diss = POC%diss * 2.0_r8
!            endif

             poc_diss = POC%diss

!-----------------------------------------------------------------------
!  decay_POC_E and decay_SiO2 set locally, modified by Tfunc
!-----------------------------------------------------------------------

!            decay_POC_E = exp(-dz_loc / (poc_diss / TfuncP(i,j)))

            decay_POC_E = exp(-dz_loc / poc_diss)
            decay_SiO2  = exp(-dz_loc / (P_SiO2%diss / TfuncS(i,j)))

!-----------------------------------------------------------------------
!  Set outgoing fluxes for non-iron pools.
!  The outoing fluxes for ballast materials are from the
!  solution of the coresponding continuous ODE across the model
!  level. The ODE has a constant source term and linear decay.
!  It is assumed that there is no sub-surface dust production.
!-----------------------------------------------------------------------

            P_CaCO3%sflux_out(i,j,bid) = P_CaCO3%sflux_in(i,j,bid) * DECAY_CaCO3(i,j) + &
               P_CaCO3%prod(i,j,bid) * ((c1 - P_CaCO3%gamma) * (c1 - DECAY_CaCO3(i,j)) &
                  * P_CaCO3%diss)

            P_CaCO3%hflux_out(i,j,bid) = P_CaCO3%hflux_in(i,j,bid) * DECAY_Hard(i,j) + &
               P_CaCO3%prod(i,j,bid) * (P_CaCO3%gamma * dz_loc)

            P_SiO2%sflux_out(i,j,bid) = P_SiO2%sflux_in(i,j,bid) * decay_SiO2 + &
               P_SiO2%prod(i,j,bid) * ((c1 - P_SiO2%gamma) * (c1 - decay_SiO2) &
                  * (P_SiO2%diss / TfuncS(i,j)))

            P_SiO2%hflux_out(i,j,bid) = P_SiO2%hflux_in(i,j,bid) * DECAY_Hard(i,j) + &
               P_SiO2%prod(i,j,bid) * (P_SiO2%gamma * dz_loc)

            dust%sflux_out(i,j,bid) = dust%sflux_in(i,j,bid) * DECAY_dust(i,j)

            dust%hflux_out(i,j,bid) = dust%hflux_in(i,j,bid) * DECAY_HardDust(i,j)

!-----------------------------------------------------------------------
!  Compute how much POC_PROD is available for deficit reduction
!  and excess POC flux after subtracting off fraction of non-dust
!  ballast production from net POC_PROD.
!-----------------------------------------------------------------------

            POC_PROD_avail = POC%prod(i,j,bid) - &
               P_CaCO3%rho * P_CaCO3%prod(i,j,bid) - &
               P_SiO2%rho * P_SiO2%prod(i,j,bid)

!-----------------------------------------------------------------------
!  Check for POC production bounds violations
!-----------------------------------------------------------------------

            if (POC_PROD_avail < c0) then
               poc_error = .true.
            endif

!-----------------------------------------------------------------------
!  Compute 1st approximation to new QA_dust_def, the QA_dust
!  deficit leaving the cell. Ignore POC_PROD_avail at this stage.
!-----------------------------------------------------------------------

            if (QA_dust_def(i,j,bid) > 0) then
               new_QA_dust_def = QA_dust_def(i,j,bid) * &
                  (dust%sflux_out(i,j,bid) + dust%hflux_out(i,j,bid)) / &
                  (dust%sflux_in(i,j,bid) + dust%hflux_in(i,j,bid))
            else
               new_QA_dust_def = c0
            endif

!-----------------------------------------------------------------------
!  Use POC_PROD_avail to reduce new_QA_dust_def.
!-----------------------------------------------------------------------

            if (new_QA_dust_def > c0) then
               new_QA_dust_def = new_QA_dust_def - POC_PROD_avail * dz_loc
               if (new_QA_dust_def < c0) then
                  POC_PROD_avail = -new_QA_dust_def * dzr_loc
                  new_QA_dust_def = c0
               else
                  POC_PROD_avail = c0
               endif
            endif

            QA_dust_def(i,j,bid) = new_QA_dust_def

!-----------------------------------------------------------------------
!  Compute outgoing POC fluxes. QA POC flux is computing using
!  ballast fluxes and new_QA_dust_def. If no QA POC flux came in
!  and no production occured, then no QA POC flux goes out. This
!  shortcut is present to avoid roundoff cancellation errors from
!  the dust%rho * dust_flux_out - QA_dust_def computation.
!  Any POC_PROD_avail still remaining goes into excess POC flux.
!-----------------------------------------------------------------------

            if (POC%hflux_in(i,j,bid) == c0 .and. POC%prod(i,j,bid) == c0) then
               POC%hflux_out(i,j,bid) = c0
            else
               POC%hflux_out(i,j,bid) = P_CaCO3%rho * &
                  (P_CaCO3%sflux_out(i,j,bid) + P_CaCO3%hflux_out(i,j,bid)) + &
                  P_SiO2%rho * &
                  (P_SiO2%sflux_out(i,j,bid) + P_SiO2%hflux_out(i,j,bid)) + &
                  dust%rho * &
                  (dust%sflux_out(i,j,bid) + dust%hflux_out(i,j,bid)) - &
                  new_QA_dust_def
               POC%hflux_out(i,j,bid) = max(POC%hflux_out(i,j,bid), c0)
            endif

            POC%sflux_out(i,j,bid) = POC%sflux_in(i,j,bid) * decay_POC_E + &
               POC_PROD_avail *((c1 - decay_POC_E) * &
               poc_diss)

!               (poc_diss / TfuncP(i,j)))

!-----------------------------------------------------------------------
!  Compute remineralization terms. It is assumed that there is no
!  sub-surface dust production.
!-----------------------------------------------------------------------

            P_CaCO3%remin(i,j,bid) = P_CaCO3%prod(i,j,bid) + &
               ((P_CaCO3%sflux_in(i,j,bid) - P_CaCO3%sflux_out(i,j,bid)) + &
               (P_CaCO3%hflux_in(i,j,bid) - P_CaCO3%hflux_out(i,j,bid))) * dzr_loc

            P_SiO2%remin(i,j,bid) = P_SiO2%prod(i,j,bid) + &
               ((P_SiO2%sflux_in(i,j,bid) - P_SiO2%sflux_out(i,j,bid)) + &
               (P_SiO2%hflux_in(i,j,bid) - P_SiO2%hflux_out(i,j,bid))) * dzr_loc

            POC%remin(i,j,bid) = POC%prod(i,j,bid) + &
               ((POC%sflux_in(i,j,bid) - POC%sflux_out(i,j,bid)) + &
               (POC%hflux_in(i,j,bid) - POC%hflux_out(i,j,bid))) * dzr_loc

            dust%remin(i,j,bid) = &
               ((dust%sflux_in(i,j,bid) - dust%sflux_out(i,j,bid)) + &
               (dust%hflux_in(i,j,bid) - dust%hflux_out(i,j,bid))) * dzr_loc

!-----------------------------------------------------------------------
!  Compute iron remineralization and flux out.
!-----------------------------------------------------------------------

            if (POC%sflux_in(i,j,bid) + POC%hflux_in(i,j,bid) == c0) then
               P_iron%remin(i,j,bid) = (POC%remin(i,j,bid) * parm_Red_Fe_C)
            else
               P_iron%remin(i,j,bid) = (POC%remin(i,j,bid) * &
                  (P_iron%sflux_in(i,j,bid) + P_iron%hflux_in(i,j,bid)) / &
                  (POC%sflux_in(i,j,bid) + POC%hflux_in(i,j,bid))) + &
                  (P_iron%sflux_in(i,j,bid) * 3.0e-6_r8)
            endif

            P_iron%sflux_out(i,j,bid) = P_iron%sflux_in(i,j,bid) + dz_loc * &
               ((c1 - P_iron%gamma) * P_iron%prod(i,j,bid) - P_iron%remin(i,j,bid))

            if (P_iron%sflux_out(i,j,bid) < c0) then
               P_iron%sflux_out(i,j,bid) = c0
               P_iron%remin(i,j,bid) = P_iron%sflux_in(i,j,bid) * dzr_loc + &
                  (c1 - P_iron%gamma) * P_iron%prod(i,j,bid)
            endif

!-----------------------------------------------------------------------
!  Compute iron release from dust remin/dissolution
!
!  dust remin gDust = 0.035 / 55.847 * 1.0e9 = 626712.0 nmolFe
!                      gFe     molFe     nmolFe
!  Also add in Fe source from sediments if applicable to this cell.
!-----------------------------------------------------------------------


            P_iron%remin(i,j,bid) = P_iron%remin(i,j,bid) &
               + dust%remin(i,j,bid) * dust_to_Fe &
               + (FESEDFLUX(i,j,k,bid) * dzr_loc)

!maltrud what is up here--jkm does not have this
            P_iron%hflux_out(i,j,bid) = P_iron%hflux_in(i,j,bid)

         else
            P_CaCO3%sflux_out(i,j,bid) = c0
            P_CaCO3%hflux_out(i,j,bid) = c0
            P_CaCO3%remin(i,j,bid) = c0

            P_SiO2%sflux_out(i,j,bid) = c0
            P_SiO2%hflux_out(i,j,bid) = c0
            P_SiO2%remin(i,j,bid) = c0

            dust%sflux_out(i,j,bid) = c0
            dust%hflux_out(i,j,bid) = c0
            dust%remin(i,j,bid) = c0

            POC%sflux_out(i,j,bid) = c0
            POC%hflux_out(i,j,bid) = c0
            POC%remin(i,j,bid) = c0

            P_iron%sflux_out(i,j,bid) = c0
            P_iron%hflux_out(i,j,bid) = c0
            P_iron%remin(i,j,bid) = c0
         endif

!-----------------------------------------------------------------------
!  Remineralize everything in bottom cell.
!-----------------------------------------------------------------------

         if (LAND_MASK(i,j,bid) .and. k == KMT(i,j,bid)) then
            P_CaCO3%remin(i,j,bid) = P_CaCO3%remin(i,j,bid) + &
               (P_CaCO3%sflux_out(i,j,bid) + P_CaCO3%hflux_out(i,j,bid)) * dzr_loc
            P_CaCO3%sflux_out(i,j,bid) = c0
            P_CaCO3%hflux_out(i,j,bid) = c0

            P_SiO2%remin(i,j,bid) = P_SiO2%remin(i,j,bid) + &
               (P_SiO2%sflux_out(i,j,bid) + P_SiO2%hflux_out(i,j,bid)) * dzr_loc
            P_SiO2%sflux_out(i,j,bid) = c0
            P_SiO2%hflux_out(i,j,bid) = c0

!           dust%remin(i,j,bid) = dust%remin(i,j,bid) + &
!              (dust%sflux_out(i,j,bid) + dust%hflux_out(i,j,bid)) * dzr_loc
!           dust%sflux_out(i,j,bid) = c0
!           dust%hflux_out(i,j,bid) = c0

            POC%remin(i,j,bid) = POC%remin(i,j,bid) + &
               (POC%sflux_out(i,j,bid) + POC%hflux_out(i,j,bid)) * dzr_loc
            POC%sflux_out(i,j,bid) = c0
            POC%hflux_out(i,j,bid) = c0

            P_iron%remin(i,j,bid) = P_iron%remin(i,j,bid) + &
               (P_iron%sflux_out(i,j,bid) + P_iron%hflux_out(i,j,bid)) * dzr_loc
            P_iron%sflux_out(i,j,bid) = c0
            P_iron%hflux_out(i,j,bid) = c0

         endif

      end do
   end do

#ifdef CCSMCOUPLED
    if (poc_error) then
       call shr_sys_abort(subname /&
                       &/ ': mass ratio of ballast ' /&
                       &/ 'production exceeds POC production')
    endif
#endif

!-----------------------------------------------------------------------
!  Set tavg variables.
!-----------------------------------------------------------------------

   if (accumulate_tavg_now(tavg_POC_FLUX_IN)) then
      WORK = POC%sflux_in(:,:,bid) + POC%hflux_in(:,:,bid)
      call accumulate_tavg_field(WORK, tavg_POC_FLUX_IN,bid,k)
   endif

   call accumulate_tavg_field(POC%prod(:,:,bid), tavg_POC_PROD,bid,k)

   call accumulate_tavg_field(POC%remin(:,:,bid), tavg_POC_REMIN,bid,k)

   if (accumulate_tavg_now(tavg_CaCO3_FLUX_IN)) then
      WORK = P_CaCO3%sflux_in(:,:,bid) + P_CaCO3%hflux_in(:,:,bid)
      call accumulate_tavg_field(WORK, tavg_CaCO3_FLUX_IN,bid,k)
   endif

   call accumulate_tavg_field(P_CaCO3%prod(:,:,bid), tavg_CaCO3_PROD,bid,k)

   call accumulate_tavg_field(P_CaCO3%remin(:,:,bid), tavg_CaCO3_REMIN,bid,k)

   if (accumulate_tavg_now(tavg_SiO2_FLUX_IN)) then
      WORK = P_SiO2%sflux_in(:,:,bid) + P_SiO2%hflux_in(:,:,bid)
      call accumulate_tavg_field(WORK, tavg_SiO2_FLUX_IN,bid,k)
   endif

   call accumulate_tavg_field(P_SiO2%prod(:,:,bid), tavg_SiO2_PROD,bid,k)

   call accumulate_tavg_field(P_SiO2%remin(:,:,bid), tavg_SiO2_REMIN,bid,k)

   if (accumulate_tavg_now(tavg_dust_FLUX_IN)) then
      WORK = dust%sflux_in(:,:,bid) + dust%hflux_in(:,:,bid)
      call accumulate_tavg_field(WORK, tavg_dust_FLUX_IN,bid,k)
   endif

   call accumulate_tavg_field(dust%remin(:,:,bid), tavg_dust_REMIN,bid,k)

   if (accumulate_tavg_now(tavg_P_iron_FLUX_IN)) then
      WORK = P_iron%sflux_in(:,:,bid) + P_iron%hflux_in(:,:,bid)
      call accumulate_tavg_field(WORK, tavg_P_iron_FLUX_IN,bid,k)
   endif

   call accumulate_tavg_field(P_iron%prod(:,:,bid), tavg_P_iron_PROD,bid,k)

   call accumulate_tavg_field(P_iron%remin(:,:,bid), tavg_P_iron_REMIN,bid,k)

!-----------------------------------------------------------------------
!EOC

 end subroutine compute_particulate_terms

!***********************************************************************
!BOP
! !IROUTINE: ecosys_init_sflux
! !INTERFACE:

 subroutine ecosys_init_sflux

! !DESCRIPTION:
!  Initialize surface flux computations for ecosys tracer module.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'ecosys_mod:ecosys_init_sflux'

   logical (log_kind) :: &
      luse_INTERP_WORK     ! does INTERP_WORK need to be allocated

   integer (int_kind) :: &
      n,                 & ! index for looping over tracers
      iblock               ! index for looping over blocks

   real (r8), dimension (nx_block,ny_block) :: WORK

   real (r8), dimension (nx_block,ny_block,12,max_blocks_clinic), target :: &
      WORK_READ            ! temporary space to read in fields

!-----------------------------------------------------------------------

   luse_INTERP_WORK = .false.

!-----------------------------------------------------------------------
!  read gas flux forcing (if required)
!-----------------------------------------------------------------------

   if ((lflux_gas_o2 .or. lflux_gas_co2) .and. &
       gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then

      luse_INTERP_WORK = .true.

!-----------------------------------------------------------------------
!  first, read ice file
!-----------------------------------------------------------------------

      allocate(fice_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
      if (trim(fice_file%input%filename) == 'unknown') &
         fice_file%input%filename = gas_flux_forcing_file

      call read_field(fice_file%input%file_fmt, &
                      fice_file%input%filename, &
                      fice_file%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         fice_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            fice_file%DATA(:,:,iblock,1,n) = c0
         fice_file%DATA(:,:,iblock,1,n) = &
            fice_file%DATA(:,:,iblock,1,n) * fice_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(fice_file%data_time, &
                              fice_file%data_inc, fice_file%interp_type, &
                              fice_file%data_next, fice_file%data_time_min_loc, &
                              fice_file%data_update, fice_file%data_type)

!-----------------------------------------------------------------------
!  next, read piston velocity file
!-----------------------------------------------------------------------

      allocate(xkw_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
      if (trim(xkw_file%input%filename) == 'unknown') &
         xkw_file%input%filename = gas_flux_forcing_file

      call read_field(xkw_file%input%file_fmt, &
                      xkw_file%input%filename, &
                      xkw_file%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         xkw_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            xkw_file%DATA(:,:,iblock,1,n) = c0
         xkw_file%DATA(:,:,iblock,1,n) = &
            xkw_file%DATA(:,:,iblock,1,n) * xkw_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(xkw_file%data_time, &
                              xkw_file%data_inc, xkw_file%interp_type, &
                              xkw_file%data_next, xkw_file%data_time_min_loc, &
                              xkw_file%data_update, xkw_file%data_type)

!-----------------------------------------------------------------------
!  last, read atmospheric pressure file
!-----------------------------------------------------------------------

      allocate(ap_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
      if (trim(ap_file%input%filename) == 'unknown') &
         ap_file%input%filename = gas_flux_forcing_file

      call read_field(ap_file%input%file_fmt, &
                      ap_file%input%filename, &
                      ap_file%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         ap_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            ap_file%DATA(:,:,iblock,1,n) = c0
         ap_file%DATA(:,:,iblock,1,n) = &
            ap_file%DATA(:,:,iblock,1,n) * ap_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(ap_file%data_time, &
                              ap_file%data_inc, ap_file%interp_type, &
                              ap_file%data_next, ap_file%data_time_min_loc, &
                              ap_file%data_update, ap_file%data_type)

    endif

!-----------------------------------------------------------------------
!  load dust flux fields (if required)
!-----------------------------------------------------------------------

   if (trim(dust_flux%input%filename) /= 'none' .and. &
       trim(dust_flux%input%filename) /= 'unknown') then

      luse_INTERP_WORK = .true.

      allocate(dust_flux%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
      if (trim(dust_flux%input%filename) == 'unknown') &
         dust_flux%input%filename = gas_flux_forcing_file

      call read_field(dust_flux%input%file_fmt, &
                      dust_flux%input%filename, &
                      dust_flux%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         dust_flux%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            dust_flux%DATA(:,:,iblock,1,n) = c0
         dust_flux%DATA(:,:,iblock,1,n) = &
            dust_flux%DATA(:,:,iblock,1,n) * dust_flux%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(dust_flux%data_time, &
                              dust_flux%data_inc, dust_flux%interp_type, &
                              dust_flux%data_next, dust_flux%data_time_min_loc, &
                              dust_flux%data_update, dust_flux%data_type)

      dust_flux%has_data = .true.
   else
      dust_flux%has_data = .false.
   endif

!-----------------------------------------------------------------------
!  load iron flux fields (if required)
!-----------------------------------------------------------------------

   if (trim(iron_flux%input%filename) /= 'none' .and. &
       trim(iron_flux%input%filename) /= 'unknown') then

      luse_INTERP_WORK = .true.

      allocate(iron_flux%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
      if (trim(iron_flux%input%filename) == 'unknown') &
         iron_flux%input%filename = gas_flux_forcing_file

      call read_field(iron_flux%input%file_fmt, &
                      iron_flux%input%filename, &
                      iron_flux%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         iron_flux%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            iron_flux%DATA(:,:,iblock,1,n) = c0
         iron_flux%DATA(:,:,iblock,1,n) = &
            iron_flux%DATA(:,:,iblock,1,n) * iron_flux%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(iron_flux%data_time, &
                              iron_flux%data_inc, iron_flux%interp_type, &
                              iron_flux%data_next, iron_flux%data_time_min_loc, &
                              iron_flux%data_update, iron_flux%data_type)

      iron_flux%has_data = .true.
   else
      iron_flux%has_data = .false.
   endif

!-----------------------------------------------------------------------
!  load nox & noy flux fields (if required)
!-----------------------------------------------------------------------

   if (trim(ndep_data_type) /= 'none' .and. &
       trim(ndep_data_type) /= 'monthly-calendar' .and. &
       trim(ndep_data_type) /= 'shr_stream') then
      call document(subname, 'ndep_data_type', ndep_data_type)
      call exit_POP(sigAbort, 'unknown ndep_data_type')
   endif

   if (trim(ndep_data_type) == 'monthly-calendar' .and. &
       trim(nox_flux_monthly%input%filename) /= 'none' .and. &
       trim(nox_flux_monthly%input%filename) /= 'unknown') then

      luse_INTERP_WORK = .true.

      allocate(nox_flux_monthly%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
      if (trim(nox_flux_monthly%input%filename) == 'unknown') &
         nox_flux_monthly%input%filename = gas_flux_forcing_file

      call read_field(nox_flux_monthly%input%file_fmt, &
                      nox_flux_monthly%input%filename, &
                      nox_flux_monthly%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         nox_flux_monthly%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            nox_flux_monthly%DATA(:,:,iblock,1,n) = c0
         nox_flux_monthly%DATA(:,:,iblock,1,n) = &
            nox_flux_monthly%DATA(:,:,iblock,1,n) * nox_flux_monthly%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(nox_flux_monthly%data_time, &
                              nox_flux_monthly%data_inc, nox_flux_monthly%interp_type, &
                              nox_flux_monthly%data_next, nox_flux_monthly%data_time_min_loc, &
                              nox_flux_monthly%data_update, nox_flux_monthly%data_type)

      nox_flux_monthly%has_data = .true.
   else
      nox_flux_monthly%has_data = .false.
   endif

   if (trim(ndep_data_type) == 'monthly-calendar' .and. &
       trim(nhy_flux_monthly%input%filename) /= 'none' .and. &
       trim(nhy_flux_monthly%input%filename) /= 'unknown') then

      luse_INTERP_WORK = .true.

      allocate(nhy_flux_monthly%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
      if (trim(nhy_flux_monthly%input%filename) == 'unknown') &
         nhy_flux_monthly%input%filename = gas_flux_forcing_file

      call read_field(nhy_flux_monthly%input%file_fmt, &
                      nhy_flux_monthly%input%filename, &
                      nhy_flux_monthly%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         nhy_flux_monthly%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            nhy_flux_monthly%DATA(:,:,iblock,1,n) = c0
         nhy_flux_monthly%DATA(:,:,iblock,1,n) = &
            nhy_flux_monthly%DATA(:,:,iblock,1,n) * nhy_flux_monthly%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(nhy_flux_monthly%data_time, &
                              nhy_flux_monthly%data_inc, nhy_flux_monthly%interp_type, &
                              nhy_flux_monthly%data_next, nhy_flux_monthly%data_time_min_loc, &
                              nhy_flux_monthly%data_update, nhy_flux_monthly%data_type)

      nhy_flux_monthly%has_data = .true.
   else
      nhy_flux_monthly%has_data = .false.
   endif

!-----------------------------------------------------------------------

#ifndef CCSMCOUPLED
   if (trim(ndep_data_type) == 'shr_stream') then
      call document(subname, 'ndep_data_type', ndep_data_type)
      call exit_POP(sigAbort, &
         'shr_stream option only supported when CCSMCOUPLED is defined')
   endif
#endif

!-----------------------------------------------------------------------
!  allocate space for interpolate_forcing
!-----------------------------------------------------------------------

   if (luse_INTERP_WORK) &
      allocate(INTERP_WORK(nx_block,ny_block,max_blocks_clinic,1))

!-----------------------------------------------------------------------
!  load iron PATCH flux fields (if required)
!-----------------------------------------------------------------------

   if (liron_patch) then

!maltrud iron patch
!  assume patch file has same normalization and format as deposition file

      allocate(IRON_PATCH_FLUX(nx_block,ny_block,max_blocks_clinic))

      if (trim(iron_flux%input%filename) == 'unknown') &
         iron_flux%input%filename = gas_flux_forcing_file

      call read_field(iron_flux%input%file_fmt, &
                      iron_flux%input%filename, &
                      iron_patch_flux_filename, &
                      IRON_PATCH_FLUX)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         where (.not. LAND_MASK(:,:,iblock)) &
            IRON_PATCH_FLUX(:,:,iblock) = c0
         iron_flux%DATA(:,:,iblock,1,n) = &
            IRON_PATCH_FLUX(:,:,iblock) * iron_flux%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

   endif

!-----------------------------------------------------------------------
!  register and set
!     fco2, the air-sea co2 gas flux
!-----------------------------------------------------------------------

   call named_field_register('SFLUX_CO2', sflux_co2_nf_ind)
   !$OMP PARALLEL DO PRIVATE(iblock,WORK)
   do iblock=1,nblocks_clinic
      WORK = c0
      call named_field_set(sflux_co2_nf_ind, iblock, WORK)
   end do
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  verify running coupled if gas fluxes use coupler forcing
!-----------------------------------------------------------------------

   if ((lflux_gas_o2 .or. lflux_gas_co2) .and. &
       (gas_flux_forcing_iopt == gas_flux_forcing_iopt_drv .or. &
        atm_co2_iopt == atm_co2_iopt_drv_prog .or. &
        atm_co2_iopt == atm_co2_iopt_drv_diag) .and. &
       .not. registry_match('lcoupled')) then
      call exit_POP(sigAbort, 'ecosys_init: ecosys module requires the ' /&
                           &/ 'flux coupler when gas_flux_forcing_opt=drv')
   endif

!-----------------------------------------------------------------------
!  get named field index for atmospheric CO2, if required
!-----------------------------------------------------------------------

   if (lflux_gas_co2 .and. atm_co2_iopt == atm_co2_iopt_drv_prog) then
      call named_field_get_index('ATM_CO2_PROG', atm_co2_nf_ind, &
                                 exit_on_err=.false.)
      if (atm_co2_nf_ind == 0) then
         call exit_POP(sigAbort, 'ecosys_init: ecosys module requires ' /&
                              &/ 'atmopsheric CO2 from the flux coupler ' /&
                              &/ 'and it is not present')
      endif
   endif

   if (lflux_gas_co2 .and. atm_co2_iopt == atm_co2_iopt_drv_diag) then
      call named_field_get_index('ATM_CO2_DIAG', atm_co2_nf_ind, &
                                 exit_on_err=.false.)
      if (atm_co2_nf_ind == 0) then
         call exit_POP(sigAbort, 'ecosys_init: ecosys module requires ' /&
                              &/ 'atmopsheric CO2 from the flux coupler ' /&
                              &/ 'and it is not present')
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: ecosys_init_interior_restore
! !INTERFACE:

 subroutine ecosys_init_interior_restore

! !DESCRIPTION:
!  Initialize interior restoring computations for ecosys tracer module.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                   & ! index for looping over tracers
      k,                   & ! index for looping over levels
      i,j,                 & ! index for looping over horiz. dims.
      iblock                 ! index for looping over blocks

   real (r8) :: &
      subsurf_fesed          ! sum of subsurface fesed values

!-----------------------------------------------------------------------
!  initialize restoring timescale (if required)
!-----------------------------------------------------------------------

   if (lrest_po4 .or. lrest_no3 .or. lrest_sio3) then

      if (lnutr_variable_restore) then

         allocate(NUTR_RESTORE_RTAU(nx_block,ny_block,max_blocks_clinic))
         allocate(NUTR_RESTORE_MAX_LEVEL(nx_block,ny_block,max_blocks_clinic))

         call read_field(nutr_variable_rest_file_fmt, &
                         nutr_variable_rest_file,     &
                         'NUTR_RESTORE_MAX_LEVEL',    &
                         NUTR_RESTORE_RTAU)

         NUTR_RESTORE_MAX_LEVEL = nint(NUTR_RESTORE_RTAU)

         call read_field(nutr_variable_rest_file_fmt, &
                         nutr_variable_rest_file,     &
                         'NUTR_RESTORE_RTAU',         &
                         NUTR_RESTORE_RTAU)

         NUTR_RESTORE_RTAU = NUTR_RESTORE_RTAU/seconds_in_day ! convert days to secs

      else

         do k = 1,km
            if (zt(k) < rest_z0) then
               nutr_rest_time_inv(k) = rest_time_inv_surf
            else if (zt(k) > rest_z1) then
               nutr_rest_time_inv(k) = rest_time_inv_deep
            else if (rest_z1 == rest_z0) then
               nutr_rest_time_inv(k) = rest_time_inv_surf + p5 * &
                    (rest_time_inv_deep - rest_time_inv_surf)
            else
               nutr_rest_time_inv(k) = rest_time_inv_surf + &
                    (zt(k) - rest_z0) / (rest_z1 - rest_z0) * &
                    (rest_time_inv_deep - rest_time_inv_surf)
            endif
         end do

      endif  !  variable restoring

   endif  ! restoring

!-----------------------------------------------------------------------
!  load restoring fields (if required)
!-----------------------------------------------------------------------

   if (lrest_po4) then

      allocate(PO4_CLIM(nx_block,ny_block,km,max_blocks_clinic))

      call read_field(po4_rest%file_fmt, &
                      po4_rest%filename, &
                      po4_rest%file_varname, &
                      PO4_CLIM)

      do iblock=1,nblocks_clinic
         do k = 1, km
            where (.not. LAND_MASK(:,:,iblock) .or. k > KMT(:,:,iblock)) &
               PO4_CLIM(:,:,k,iblock) = c0
            PO4_CLIM(:,:,k,iblock) = PO4_CLIM(:,:,k,iblock) * po4_rest%scale_factor
         enddo
      end do

   endif

   if (lrest_no3) then

      allocate(NO3_CLIM(nx_block,ny_block,km,max_blocks_clinic))

      call read_field(no3_rest%file_fmt, &
                      no3_rest%filename, &
                      no3_rest%file_varname, &
                      NO3_CLIM)

      do iblock=1,nblocks_clinic
         do k = 1, km
            where (.not. LAND_MASK(:,:,iblock) .or. k > KMT(:,:,iblock)) &
               NO3_CLIM(:,:,k,iblock) = c0
            NO3_CLIM(:,:,k,iblock) = NO3_CLIM(:,:,k,iblock) * no3_rest%scale_factor
         enddo
      end do

   endif

   if (lrest_sio3) then

      allocate(SiO3_CLIM(nx_block,ny_block,km,max_blocks_clinic))

      call read_field(sio3_rest%file_fmt, &
                      sio3_rest%filename, &
                      sio3_rest%file_varname, &
                      SiO3_CLIM)

      do iblock=1,nblocks_clinic
         do k = 1, km
            where (.not. LAND_MASK(:,:,iblock) .or. k > KMT(:,:,iblock)) &
               SiO3_CLIM(:,:,k,iblock) = c0
            SiO3_CLIM(:,:,k,iblock) = SiO3_CLIM(:,:,k,iblock) * sio3_rest%scale_factor
         enddo
      end do

   endif

!-----------------------------------------------------------------------
!  load fesedflux
!  add subsurface positives to 1 level shallower, to accomodate overflow pop-ups
!-----------------------------------------------------------------------

   allocate(FESEDFLUX(nx_block,ny_block,km,max_blocks_clinic))

   call read_field(fesedflux_input%file_fmt, &
                   fesedflux_input%filename, &
                   fesedflux_input%file_varname, &
                   FESEDFLUX)

   do iblock=1,nblocks_clinic
      do j=1,ny_block
      do i=1,nx_block
         if (KMT(i,j,iblock) > 0 .and. KMT(i,j,iblock) < km) then
            subsurf_fesed = c0
            do k=KMT(i,j,iblock)+1,km
               subsurf_fesed = subsurf_fesed + FESEDFLUX(i,j,k,iblock)
            enddo
            FESEDFLUX(i,j,KMT(i,j,iblock),iblock) = FESEDFLUX(i,j,KMT(i,j,iblock),iblock) + subsurf_fesed
         endif
      enddo
      enddo

      do k = 1, km
         where (.not. LAND_MASK(:,:,iblock) .or. k > KMT(:,:,iblock)) &
            FESEDFLUX(:,:,k,iblock) = c0
         FESEDFLUX(:,:,k,iblock) = FESEDFLUX(:,:,k,iblock) * fesedflux_input%scale_factor
      enddo
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_init_interior_restore

!***********************************************************************
!BOP
! !IROUTINE: ecosys_set_sflux
! !INTERFACE:

 subroutine ecosys_set_sflux(SHF_QSW_RAW, SHF_QSW, &
                             U10_SQR,IFRAC,PRESS,SST,SSS, &
                             SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)
   use shr_pio_mod, only : shr_pio_getiotype, shr_pio_getiosys
   use POP_IOUnitsMod, only: inst_name

! !DESCRIPTION:
!  Compute surface fluxes for ecosys tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      SHF_QSW_RAW,  &! penetrative solar heat flux, from coupler (degC*cm/s)
      SHF_QSW,      &! SHF_QSW used by physics, may have diurnal cylce imposed (degC*cm/s)
      U10_SQR,      &! 10m wind speed squared (cm/s)**2
      IFRAC,        &! sea ice fraction (non-dimensional)
      PRESS,        &! sea level atmospheric pressure (dyne/cm**2)
      SST,          &! sea surface temperature (C)
      SSS            ! sea surface salinity (psu)

  !real (r8), dimension(nx_block,ny_block,ecosys_tracer_cnt,max_blocks_clinic), &
   real (r8), dimension(:,:,:,:), &
      intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

! !INPUT/OUTPUT PARAMETERS:

  !real (r8), dimension(nx_block,ny_block,ecosys_tracer_cnt,max_blocks_clinic), &
   real (r8), dimension(:,:,:,:), &
      intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'ecosys_mod:ecosys_set_sflux'

   logical (log_kind) :: first_call = .true.

   type (block) :: &
      this_block      ! block info for the current block

   integer (int_kind) :: &
      i,j,iblock,n, & ! loop indices
      mcdate,sec,   & ! date vals for shr_strdata_advance
      errorCode       ! errorCode from HaloUpdate call

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      IFRAC_USED,   & ! used ice fraction (non-dimensional)
      XKW_USED,     & ! portion of piston velocity (cm/s)
      AP_USED,      & ! used atm pressure (converted from dyne/cm**2 to atm)
      IRON_FLUX_IN    ! iron flux

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      SHR_STREAM_WORK

   real (r8), dimension(nx_block,ny_block) :: &
      XKW_ICE,      & ! common portion of piston vel., (1-fice)*xkw (cm/s)
      SCHMIDT_USED, & ! used Schmidt number
      PV,           & ! piston velocity (cm/s)
      O2SAT_1atm,   & ! O2 saturation @ 1 atm (mmol/m^3)
      O2SAT_USED,   & ! used O2 saturation (mmol/m^3)
      XCO2,         & ! atmospheric co2 conc. (dry-air, 1 atm)
      FLUX            ! tracer flux (nmol/cm^2/s)

   real (r8), dimension(nx_block) :: &
      PHLO,         & ! lower bound for ph in solver
      PHHI,         & ! upper bound for ph in solver
      PH_NEW,       & ! computed PH from solver
      DIC_ROW,      & ! row of DIC values for solver
      ALK_ROW,      & ! row of ALK values for solver
      PO4_ROW,      & ! row of PO4 values for solver
      SiO3_ROW,     & ! row of SiO3 values for solver
      CO2STAR_ROW,  & ! CO2STAR from solver
      DCO2STAR_ROW, & ! DCO2STAR from solver
      pCO2SURF_ROW, & ! pCO2SURF from solver
      DpCO2_ROW       ! DpCO2 from solver

   character (char_len) :: &
      tracer_data_label,       & ! label for what is being updated
      ndep_shr_stream_fldList

   character (char_len), dimension(1) :: &
      tracer_data_names          ! short names for input data fields

   integer (int_kind), dimension(1) :: &
      tracer_bndy_loc,         & ! location and field type for ghost
      tracer_bndy_type           !    cell updates

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1, WORK2 ! temporaries for averages

   real (r8) :: scalar_temp

!-----------------------------------------------------------------------
!  local parameters
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      xkw_coeff = 8.6e-9_r8     ! a = 0.31 cm/hr s^2/m^2 in (s/cm)

!-----------------------------------------------------------------------

   call timer_start(ecosys_sflux_timer)

!-----------------------------------------------------------------------

   if (check_time_flag(comp_surf_avg_flag))  &
      call comp_surf_avg(SURF_VALS_OLD,SURF_VALS_CUR)

!-----------------------------------------------------------------------
!  fluxes initially set to 0
!  set Chl field for short-wave absorption
!  store incoming shortwave in PAR_out field, converting to W/m^2
!-----------------------------------------------------------------------

   scalar_temp = f_qsw_par / hflux_factor

   !$OMP PARALLEL DO PRIVATE(iblock,WORK1)
   do iblock = 1, nblocks_clinic
      STF_MODULE(:,:,:,iblock) = c0

      WORK1 = max(c0,p5*(SURF_VALS_OLD(:,:,spChl_ind,iblock) + &
                         SURF_VALS_CUR(:,:,spChl_ind,iblock))) + &
              max(c0,p5*(SURF_VALS_OLD(:,:,diatChl_ind,iblock) + &
                         SURF_VALS_CUR(:,:,diatChl_ind,iblock))) + &
              max(c0,p5*(SURF_VALS_OLD(:,:,diazChl_ind,iblock) + &
                         SURF_VALS_CUR(:,:,diazChl_ind,iblock)))
      call named_field_set(totChl_surf_nf_ind, iblock, WORK1)

      if (ecosys_qsw_distrb_const) then
         PAR_out(:,:,iblock) = SHF_QSW_RAW(:,:,iblock)
      else
         PAR_out(:,:,iblock) = SHF_QSW(:,:,iblock)
      endif

      where (LAND_MASK(:,:,iblock))
         PAR_out(:,:,iblock) = max(c0, scalar_temp * PAR_out(:,:,iblock))
      elsewhere
         PAR_out(:,:,iblock) = c0
      end where
   enddo
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  Interpolate gas flux forcing data if necessary
!-----------------------------------------------------------------------

   if ((lflux_gas_o2 .or. lflux_gas_co2) .and. &
        gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then
      if (thour00 >= fice_file%data_update) then
         tracer_data_names = fice_file%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'Ice Fraction'
         call update_forcing_data(fice_file%data_time,      &
            fice_file%data_time_min_loc,  fice_file%interp_type,    &
            fice_file%data_next,          fice_file%data_update,    &
            fice_file%data_type,          fice_file%data_inc,       &
            fice_file%DATA(:,:,:,:,1:12), fice_file%data_renorm,    &
            tracer_data_label,            tracer_data_names,        &
            tracer_bndy_loc,              tracer_bndy_type,         &
            fice_file%filename,           fice_file%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK, &
         fice_file%DATA(:,:,:,:,1:12), &
         fice_file%data_time,         fice_file%interp_type, &
         fice_file%data_time_min_loc, fice_file%interp_freq, &
         fice_file%interp_inc,        fice_file%interp_next, &
         fice_file%interp_last,       0)
      IFRAC_USED = INTERP_WORK(:,:,:,1)

      if (thour00 >= xkw_file%data_update) then
         tracer_data_names = xkw_file%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'Piston Velocity'
         call update_forcing_data(xkw_file%data_time,      &
            xkw_file%data_time_min_loc,  xkw_file%interp_type,    &
            xkw_file%data_next,          xkw_file%data_update,    &
            xkw_file%data_type,          xkw_file%data_inc,       &
            xkw_file%DATA(:,:,:,:,1:12), xkw_file%data_renorm,    &
            tracer_data_label,           tracer_data_names,       &
            tracer_bndy_loc,             tracer_bndy_type,        &
            xkw_file%filename,           xkw_file%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK,     &
         xkw_file%DATA(:,:,:,:,1:12), &
         xkw_file%data_time,         xkw_file%interp_type, &
         xkw_file%data_time_min_loc, xkw_file%interp_freq, &
         xkw_file%interp_inc,        xkw_file%interp_next, &
         xkw_file%interp_last,       0)
      XKW_USED = INTERP_WORK(:,:,:,1)

      if (thour00 >= ap_file%data_update) then
         tracer_data_names = ap_file%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'Atmospheric Pressure'
         call update_forcing_data(ap_file%data_time,    &
            ap_file%data_time_min_loc,  ap_file%interp_type,    &
            ap_file%data_next,          ap_file%data_update,    &
            ap_file%data_type,          ap_file%data_inc,       &
            ap_file%DATA(:,:,:,:,1:12), ap_file%data_renorm,    &
            tracer_data_label,          tracer_data_names,      &
            tracer_bndy_loc,            tracer_bndy_type,       &
            ap_file%filename,           ap_file%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK, &
         ap_file%DATA(:,:,:,:,1:12), &
         ap_file%data_time,         ap_file%interp_type, &
         ap_file%data_time_min_loc, ap_file%interp_freq, &
         ap_file%interp_inc,        ap_file%interp_next, &
         ap_file%interp_last,       0)
      AP_USED = INTERP_WORK(:,:,:,1)

   endif

!-----------------------------------------------------------------------
!  calculate gas flux quantities if necessary
!-----------------------------------------------------------------------

   if (lflux_gas_o2 .or. lflux_gas_co2) then

      !$OMP PARALLEL DO PRIVATE(iblock,j,XKW_ICE,SCHMIDT_USED,PV,O2SAT_USED, &
      !$OMP                     O2SAT_1atm,FLUX,XCO2,PHLO,PHHI,DIC_ROW,ALK_ROW, &
      !$OMP                     PO4_ROW,SiO3_ROW,PH_NEW,CO2STAR_ROW, &
      !$OMP                     DCO2STAR_ROW,pCO2SURF_ROW,DpCO2_ROW)

      do iblock = 1, nblocks_clinic

!-----------------------------------------------------------------------
!  Apply OCMIP ice fraction mask when input is from a file.
!-----------------------------------------------------------------------

         if (gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then
            where (IFRAC_USED(:,:,iblock) < 0.2000_r8) &
               IFRAC_USED(:,:,iblock) = 0.2000_r8
            where (IFRAC_USED(:,:,iblock) > 0.9999_r8) &
               IFRAC_USED(:,:,iblock) = 0.9999_r8
         endif

         if (gas_flux_forcing_iopt == gas_flux_forcing_iopt_drv) then
            IFRAC_USED(:,:,iblock) = IFRAC(:,:,iblock)
            where (IFRAC_USED(:,:,iblock) < c0) IFRAC_USED(:,:,iblock) = c0
            where (IFRAC_USED(:,:,iblock) > c1) IFRAC_USED(:,:,iblock) = c1
            XKW_USED(:,:,iblock) = xkw_coeff * U10_SQR(:,:,iblock)
            AP_USED(:,:,iblock) = PRESS(:,:,iblock)
         endif

!-----------------------------------------------------------------------
!  assume PRESS is in cgs units (dyne/cm**2) since that is what is
!    required for pressure forcing in barotropic
!  want units to be atmospheres
!  convertion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
!  convertion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8
!-----------------------------------------------------------------------

         AP_USED(:,:,iblock) = PRESS(:,:,iblock) / 101.325e+4_r8

         ECO_SFLUX_TAVG(:,:,1,iblock) = IFRAC_USED(:,:,iblock)
         ECO_SFLUX_TAVG(:,:,2,iblock) = XKW_USED(:,:,iblock)
         ECO_SFLUX_TAVG(:,:,3,iblock) = AP_USED(:,:,iblock)

!-----------------------------------------------------------------------
!  Compute XKW_ICE. XKW is zero over land, so XKW_ICE is too.
!-----------------------------------------------------------------------

         XKW_ICE = (c1 - IFRAC_USED(:,:,iblock)) * XKW_USED(:,:,iblock)

!-----------------------------------------------------------------------
!  compute O2 flux
!-----------------------------------------------------------------------

         if (lflux_gas_o2) then
            SCHMIDT_USED = SCHMIDT_O2(SST(:,:,iblock), LAND_MASK(:,:,iblock))

            O2SAT_1atm = O2SAT(SST(:,:,iblock),SSS(:,:,iblock),   &
                               LAND_MASK(:,:,iblock))

            where (LAND_MASK(:,:,iblock))
               PV = XKW_ICE * SQRT(660.0_r8 / SCHMIDT_USED)
               O2SAT_USED = AP_USED(:,:,iblock) * O2SAT_1atm
               FLUX = PV * (O2SAT_USED - p5*(SURF_VALS_OLD(:,:,o2_ind,iblock) + &
                                             SURF_VALS_CUR(:,:,o2_ind,iblock)))
               STF_MODULE(:,:,o2_ind,iblock) = FLUX
            elsewhere
               O2SAT_USED = c0
               STF_MODULE(:,:,o2_ind,iblock) = c0
            end where

            ECO_SFLUX_TAVG(:,:,4,iblock) = PV
            ECO_SFLUX_TAVG(:,:,5,iblock) = SCHMIDT_USED
            ECO_SFLUX_TAVG(:,:,6,iblock) = O2SAT_USED
            ECO_SFLUX_TAVG(:,:,16,iblock) = STF_MODULE(:,:,o2_ind,iblock)

         endif  ! lflux_gas_o2

!-----------------------------------------------------------------------
!  compute CO2 flux, computing disequilibrium one row at a time
!-----------------------------------------------------------------------

         if (lflux_gas_co2) then

            SCHMIDT_USED = SCHMIDT_CO2(SST(:,:,iblock), LAND_MASK(:,:,iblock))

            where (LAND_MASK(:,:,iblock))
               PV = XKW_ICE * SQRT(660.0_r8 / SCHMIDT_USED)
            elsewhere
               PV = c0
            end where

!-----------------------------------------------------------------------
!  Set XCO2
!-----------------------------------------------------------------------

            select case (atm_co2_iopt)
            case (atm_co2_iopt_const)
               XCO2 = atm_co2_const
            case (atm_co2_iopt_drv_prog, atm_co2_iopt_drv_diag)
               call named_field_get(atm_co2_nf_ind, iblock, XCO2)
            end select

            do j = 1,ny_block
               where (PH_PREV(:,j,iblock) /= c0)
                  PHLO = PH_PREV(:,j,iblock) - del_ph
                  PHHI = PH_PREV(:,j,iblock) + del_ph
               elsewhere
                  PHLO = phlo_surf_init
                  PHHI = phhi_surf_init
               end where

               DIC_ROW = p5*(SURF_VALS_OLD(:,j,dic_ind,iblock) + &
                             SURF_VALS_CUR(:,j,dic_ind,iblock))
               ALK_ROW = p5*(SURF_VALS_OLD(:,j,alk_ind,iblock) + &
                             SURF_VALS_CUR(:,j,alk_ind,iblock))
               PO4_ROW = p5*(SURF_VALS_OLD(:,j,po4_ind,iblock) + &
                             SURF_VALS_CUR(:,j,po4_ind,iblock))
               SiO3_ROW = p5*(SURF_VALS_OLD(:,j,sio3_ind,iblock) + &
                              SURF_VALS_CUR(:,j,sio3_ind,iblock))

               call co2calc_row(iblock, j, LAND_MASK(:,j,iblock), locmip_k1_k2_bug_fix, &
                                SST(:,j,iblock), SSS(:,j,iblock), &
                                DIC_ROW, ALK_ROW, PO4_ROW, SiO3_ROW, &
                                PHLO, PHHI, PH_NEW, XCO2(:,j), &
                                AP_USED(:,j,iblock), CO2STAR_ROW, &
                                DCO2STAR_ROW, pCO2SURF_ROW, DpCO2_ROW)

               PH_PREV(:,j,iblock) = PH_NEW

               FLUX(:,j) = PV(:,j) * DCO2STAR_ROW

               ECO_SFLUX_TAVG(:,j, 7,iblock) = CO2STAR_ROW
               ECO_SFLUX_TAVG(:,j, 8,iblock) = DCO2STAR_ROW
               ECO_SFLUX_TAVG(:,j, 9,iblock) = pCO2SURF_ROW
               ECO_SFLUX_TAVG(:,j,10,iblock) = DpCO2_ROW

            end do

!-----------------------------------------------------------------------
!  set air-sea co2 gas flux named field, converting units from
!  nmol/cm^2/s (positive down) to kg CO2/m^2/s (positive down)
!-----------------------------------------------------------------------

            call named_field_set(sflux_co2_nf_ind, iblock, 44.0e-8_r8 * FLUX)

            STF_MODULE(:,:,dic_ind,iblock) = STF_MODULE(:,:,dic_ind,iblock) + FLUX

            ECO_SFLUX_TAVG(:,:,11,iblock) = PV
            ECO_SFLUX_TAVG(:,:,12,iblock) = SCHMIDT_USED
            ECO_SFLUX_TAVG(:,:,13,iblock) = FLUX
            ECO_SFLUX_TAVG(:,:,14,iblock) = PH_PREV(:,:,iblock)
            ECO_SFLUX_TAVG(:,:,15,iblock) = XCO2

         endif  !  lflux_gas_co2

      enddo
      !$OMP END PARALLEL DO

   endif  ! lflux_gas_o2 .or. lflux_gas_co2

!-----------------------------------------------------------------------
!  calculate iron and dust fluxes if necessary
!-----------------------------------------------------------------------

   if (iron_flux%has_data) then
      if (thour00 >= iron_flux%data_update) then
         tracer_data_names = iron_flux%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'Iron Flux'
         call update_forcing_data(iron_flux%data_time,    &
            iron_flux%data_time_min_loc,  iron_flux%interp_type,    &
            iron_flux%data_next,          iron_flux%data_update,    &
            iron_flux%data_type,          iron_flux%data_inc,       &
            iron_flux%DATA(:,:,:,:,1:12), iron_flux%data_renorm,    &
            tracer_data_label,            tracer_data_names,        &
            tracer_bndy_loc,              tracer_bndy_type,         &
            iron_flux%filename,           iron_flux%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK,     &
         iron_flux%DATA(:,:,:,:,1:12), &
         iron_flux%data_time,         iron_flux%interp_type, &
         iron_flux%data_time_min_loc, iron_flux%interp_freq, &
         iron_flux%interp_inc,        iron_flux%interp_next, &
         iron_flux%interp_last,       0)
      if (liron_patch .and. imonth == iron_patch_month) then
         IRON_FLUX_IN = INTERP_WORK(:,:,:,1) + IRON_PATCH_FLUX
      else
         IRON_FLUX_IN = INTERP_WORK(:,:,:,1)
      endif
   else
      IRON_FLUX_IN = c0
   endif

   IRON_FLUX_IN = IRON_FLUX_IN * parm_Fe_bioavail

   STF_MODULE(:,:,fe_ind,:) = IRON_FLUX_IN

   if (dust_flux%has_data) then
      if (thour00 >= dust_flux%data_update) then
         tracer_data_names = dust_flux%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'Dust Flux'
         call update_forcing_data(dust_flux%data_time,    &
            dust_flux%data_time_min_loc,  dust_flux%interp_type,    &
            dust_flux%data_next,          dust_flux%data_update,    &
            dust_flux%data_type,          dust_flux%data_inc,       &
            dust_flux%DATA(:,:,:,:,1:12), dust_flux%data_renorm,    &
            tracer_data_label,            tracer_data_names,        &
            tracer_bndy_loc,              tracer_bndy_type,         &
            dust_flux%filename,           dust_flux%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK, &
         dust_flux%DATA(:,:,:,:,1:12),    &
         dust_flux%data_time,         dust_flux%interp_type, &
         dust_flux%data_time_min_loc, dust_flux%interp_freq, &
         dust_flux%interp_inc,        dust_flux%interp_next, &
         dust_flux%interp_last,       0)
      dust_FLUX_IN = INTERP_WORK(:,:,:,1)

!-----------------------------------------------------------------------
!  Reduce surface dust flux due to assumed instant surface dissolution
!-----------------------------------------------------------------------

      dust_FLUX_IN = dust_FLUX_IN * (c1 - parm_Fe_bioavail)
   else
      dust_FLUX_IN = c0
   endif

!-----------------------------------------------------------------------
!  calculate nox and nhy fluxes if necessary
!-----------------------------------------------------------------------

   if (nox_flux_monthly%has_data) then
      if (thour00 >= nox_flux_monthly%data_update) then
         tracer_data_names = nox_flux_monthly%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'NOx Flux'
         call update_forcing_data(nox_flux_monthly%data_time,    &
            nox_flux_monthly%data_time_min_loc,  nox_flux_monthly%interp_type, &
            nox_flux_monthly%data_next,          nox_flux_monthly%data_update, &
            nox_flux_monthly%data_type,          nox_flux_monthly%data_inc,    &
            nox_flux_monthly%DATA(:,:,:,:,1:12), nox_flux_monthly%data_renorm, &
            tracer_data_label,                   tracer_data_names,            &
            tracer_bndy_loc,                     tracer_bndy_type,             &
            nox_flux_monthly%filename,           nox_flux_monthly%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK,     &
         nox_flux_monthly%DATA(:,:,:,:,1:12), &
         nox_flux_monthly%data_time,         nox_flux_monthly%interp_type, &
         nox_flux_monthly%data_time_min_loc, nox_flux_monthly%interp_freq, &
         nox_flux_monthly%interp_inc,        nox_flux_monthly%interp_next, &
         nox_flux_monthly%interp_last,       0)
      STF_MODULE(:,:,no3_ind,:) = INTERP_WORK(:,:,:,1)
   endif

   if (nhy_flux_monthly%has_data) then
      if (thour00 >= nhy_flux_monthly%data_update) then
         tracer_data_names = nhy_flux_monthly%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'NHy Flux'
         call update_forcing_data(nhy_flux_monthly%data_time,    &
            nhy_flux_monthly%data_time_min_loc,  nhy_flux_monthly%interp_type, &
            nhy_flux_monthly%data_next,          nhy_flux_monthly%data_update, &
            nhy_flux_monthly%data_type,          nhy_flux_monthly%data_inc,    &
            nhy_flux_monthly%DATA(:,:,:,:,1:12), nhy_flux_monthly%data_renorm, &
            tracer_data_label,                   tracer_data_names,            &
            tracer_bndy_loc,                     tracer_bndy_type,             &
            nhy_flux_monthly%filename,           nhy_flux_monthly%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK,     &
         nhy_flux_monthly%DATA(:,:,:,:,1:12), &
         nhy_flux_monthly%data_time,         nhy_flux_monthly%interp_type, &
         nhy_flux_monthly%data_time_min_loc, nhy_flux_monthly%interp_freq, &
         nhy_flux_monthly%interp_inc,        nhy_flux_monthly%interp_next, &
         nhy_flux_monthly%interp_last,       0)
      STF_MODULE(:,:,nh4_ind,:) = INTERP_WORK(:,:,:,1)
   endif

#ifdef CCSMCOUPLED
   if (trim(ndep_data_type) == 'shr_stream') then
      if (first_call) then

         ndep_shr_stream_fldList = ' '
         do n = 1, ndep_shr_stream_var_cnt
            if (n == ndep_shr_stream_no_ind) &
               ndep_shr_stream_fldList = trim(ndep_shr_stream_fldList) /&
                  &/ 'NOy_deposition'
            if (n == ndep_shr_stream_nh_ind) &
               ndep_shr_stream_fldList = trim(ndep_shr_stream_fldList) /&
                  &/ 'NHx_deposition'
            if (n < ndep_shr_stream_var_cnt) &
               ndep_shr_stream_fldList = trim(ndep_shr_stream_fldList) /&
                  &/ ':'
         end do

         call shr_strdata_create(ndep_sdat,name='ndep data',                   &
                                 mpicom=POP_communicator,                      &
                                 compid=POP_MCT_OCNID,                         &
                                 gsmap=POP_MCT_gsMap_o, ggrid=POP_MCT_dom_o,   &
                                 nxg=nx_global, nyg=ny_global,                 &
                                 yearFirst=ndep_shr_stream_year_first,         &
                                 yearLast=ndep_shr_stream_year_last,           &
                                 yearAlign=ndep_shr_stream_year_align,         &
                                 offset=0,                                     &
                                 domFilePath='',                               &
                                 domFileName=ndep_shr_stream_file,             &
                                 domTvarName='time',                           &
                                 domXvarName='TLONG', domYvarName='TLAT',      &
                                 domAreaName='TAREA', domMaskName='KMT',       &
                                 FilePath='',                                  &
                                 FileName=(/trim(ndep_shr_stream_file)/),      &
                                 fldListFile=ndep_shr_stream_fldList,          &
                                 fldListModel=ndep_shr_stream_fldList,         &
                                 pio_subsystem=shr_pio_getiosys(inst_name),     &
                                 pio_iotype=shr_pio_getiotype(inst_name),       &
                                 fillalgo='none', mapalgo='none')
         if (my_task == master_task) then
            call shr_strdata_print(ndep_sdat)
         endif
         first_call = .false.
      endif

      mcdate = iyear*10000 + imonth*100 + iday
      sec = isecond + 60 * (iminute + 60 * ihour)

      call timer_start(ecosys_shr_strdata_advance_timer)
      call shr_strdata_advance(ndep_sdat, mcdate, sec, &
                               POP_communicator, 'ndep')
      call timer_stop(ecosys_shr_strdata_advance_timer)

      !
      ! process NO3 flux, store results in SHR_STREAM_WORK array
      ! instead of directly into STF_MODULE
      ! to avoid argument copies in HaloUpdate calls
      !
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            SHR_STREAM_WORK(i,j,iblock) = &
               ndep_sdat%avs(1)%rAttr(ndep_shr_stream_no_ind,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(SHR_STREAM_WORK,POP_haloClinic, &
                          POP_gridHorzLocCenter,          &
                          POP_fieldKindScalar, errorCode, &
                          fillValue = 0.0_POP_r8)
      if (errorCode /= POP_Success) then
         call exit_POP(sigAbort, subname /&
            &/ ': error updating halo for Ndep fields')
      endif

      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1, nblocks_clinic
         where (LAND_MASK(:,:,iblock))
            STF_MODULE(:,:,no3_ind,iblock) = &
               ndep_shr_stream_scale_factor * SHR_STREAM_WORK(:,:,iblock)
         elsewhere
            STF_MODULE(:,:,no3_ind,iblock) = c0
         endwhere
      enddo
      !$OMP END PARALLEL DO

      !
      ! process NH4 flux, store results in SHR_STREAM_WORK array
      ! instead of directly into STF_MODULE
      ! to avoid argument copies in HaloUpdate calls
      !
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            SHR_STREAM_WORK(i,j,iblock) = &
               ndep_sdat%avs(1)%rAttr(ndep_shr_stream_nh_ind,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(SHR_STREAM_WORK,POP_haloClinic, &
                          POP_gridHorzLocCenter,          &
                          POP_fieldKindScalar, errorCode, &
                          fillValue = 0.0_POP_r8)
      if (errorCode /= POP_Success) then
         call exit_POP(sigAbort, subname /&
            &/ ': error updating halo for Ndep fields')
      endif

      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1, nblocks_clinic
         where (LAND_MASK(:,:,iblock))
            STF_MODULE(:,:,nh4_ind,iblock) = &
               ndep_shr_stream_scale_factor * SHR_STREAM_WORK(:,:,iblock)
         elsewhere
            STF_MODULE(:,:,nh4_ind,iblock) = c0
         endwhere
      enddo
      !$OMP END PARALLEL DO

   endif
#endif

   call timer_stop(ecosys_sflux_timer)

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_set_sflux

!*****************************************************************************
!BOP
! !IROUTINE: SCHMIDT_O2
! !INTERFACE:

 function SCHMIDT_O2(SST, LAND_MASK)

! !DESCRIPTION:
!  Compute Schmidt number of O2 in seawater as function of SST
!  where LAND_MASK is true. Give zero where LAND_MASK is false.
!
!  ref : Keeling et al, Global Biogeochem. Cycles, Vol. 12,
!        No. 1, pp. 141-163, March 1998
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(in) :: SST

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: SCHMIDT_O2

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      a = 1638.0_r8, &
      b = 81.83_r8, &
      c = 1.483_r8, &
      d = 0.008004_r8

!-----------------------------------------------------------------------

   where (LAND_MASK)
      SCHMIDT_O2 = a + SST * (-b + SST * (c + SST * (-d)))
   elsewhere
      SCHMIDT_O2 = c0
   end where

!-----------------------------------------------------------------------
!EOC

 end function SCHMIDT_O2

!*****************************************************************************
!BOP
! !IROUTINE: O2SAT
! !INTERFACE:

 function O2SAT(SST, SSS, LAND_MASK)

! !DESCRIPTION:
!
!  Computes oxygen saturation concentration at 1 atm total pressure
!  in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
!  in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
!
!  FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
!  THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
!
!  *** NOTE: THE "A_3*TS^2" TERM (IN THE PAPER) IS INCORRECT. ***
!  *** IT SHOULD NOT BE THERE.                                ***
!
!  O2SAT IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
!  0 permil <= S <= 42 permil
!  CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil,
!  O2SAT = 282.015 mmol/m^3
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SST, & ! sea surface temperature (C)
      SSS    ! sea surface salinity (psu)

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK

! !OUTPUT PARAMETERS:

    real (r8), dimension(nx_block,ny_block) :: O2SAT

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block) :: TS

!-----------------------------------------------------------------------
!  coefficients in expansion
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      a_0 = 2.00907_r8, &
      a_1 = 3.22014_r8, &
      a_2 = 4.05010_r8, &
      a_3 = 4.94457_r8, &
      a_4 = -2.56847E-1_r8, &
      a_5 = 3.88767_r8, &
      b_0 = -6.24523E-3_r8, &
      b_1 = -7.37614E-3_r8, &
      b_2 = -1.03410E-2_r8, &
      b_3 = -8.17083E-3_r8, &
      c_0 = -4.88682E-7_r8

   where (LAND_MASK)
      TS = log( ((T0_Kelvin+25.0_r8) - SST) / (T0_Kelvin + SST) )

      O2SAT = exp(a_0+TS*(a_1+TS*(a_2+TS*(a_3+TS*(a_4+TS*a_5)))) + &
         SSS*( (b_0+TS*(b_1+TS*(b_2+TS*b_3))) + SSS*c_0 ))
   elsewhere
      O2SAT = c0
   end where

!-----------------------------------------------------------------------
!  Convert from ml/l to mmol/m^3
!-----------------------------------------------------------------------

   O2SAT = O2SAT / 0.0223916_r8

!-----------------------------------------------------------------------
!EOC

 end function O2SAT

!*****************************************************************************
!BOP
! !IROUTINE: SCHMIDT_CO2
! !INTERFACE:

 function SCHMIDT_CO2(SST, LAND_MASK)

! !DESCRIPTION:
!  Compute Schmidt number of CO2 in seawater as function of SST
!  where LAND_MASK is true. Give zero where LAND_MASK is false.
!
!  ref : Wanninkhof, J. Geophys. Res, Vol. 97, No. C5,
!  pp. 7373-7382, May 15, 1992
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(in) :: SST

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: SCHMIDT_CO2

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      a = 2073.1_r8, &
      b = 125.62_r8, &
      c = 3.6276_r8, &
      d = 0.043219_r8

   where (LAND_MASK)
      SCHMIDT_CO2 = a + SST * (-b + SST * (c + SST * (-d)))
   elsewhere
      SCHMIDT_CO2 = c0
   end where

!-----------------------------------------------------------------------
!EOC

 end function SCHMIDT_CO2

!*****************************************************************************
!BOP
! !IROUTINE: ecosys_tavg_forcing
! !INTERFACE:

 subroutine ecosys_tavg_forcing(STF_MODULE)

! !DESCRIPTION:
!  Accumulate non-standard forcing related tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

 !real (r8), dimension(nx_block,ny_block,ecosys_tracer_cnt,max_blocks_clinic), &
  real (r8), dimension(:,:,:,:), &
     intent(in) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock              ! block loop index

!-----------------------------------------------------------------------
!  multiply IRON, DUST fluxes by mpercm (.01) to convert from model
!    units (cm/s)(mmol/m^3) to mmol/s/m^2
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)

   do iblock = 1,nblocks_clinic

      call accumulate_tavg_field(STF_MODULE(:,:,fe_ind,iblock)*mpercm, &
                                 tavg_IRON_FLUX,iblock,1)
      call accumulate_tavg_field(dust_FLUX_IN(:,:,iblock)*mpercm,  &
                                    tavg_DUST_FLUX,iblock,1)
      call accumulate_tavg_field(STF_MODULE(:,:,no3_ind,iblock),  &
                                 tavg_NOx_FLUX,iblock,1)
      call accumulate_tavg_field(STF_MODULE(:,:,nh4_ind,iblock),  &
                                 tavg_NHy_FLUX,iblock,1)

!-----------------------------------------------------------------------
!  now do components of flux calculations saved in hack array
!-----------------------------------------------------------------------

      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,1,iblock),  &
                                 tavg_ECOSYS_IFRAC,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,1,iblock),  &
                                 tavg_ECOSYS_IFRAC_2,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,2,iblock),  &
                                 tavg_ECOSYS_XKW,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,2,iblock),  &
                                 tavg_ECOSYS_XKW_2,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,3,iblock),  &
                                 tavg_ECOSYS_ATM_PRESS,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,4,iblock),  &
                                 tavg_PV_O2,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,5,iblock),  &
                                 tavg_SCHMIDT_O2,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,6,iblock),  &
                                 tavg_O2SAT,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,7,iblock),  &
                                 tavg_CO2STAR,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,8,iblock),  &
                                 tavg_DCO2STAR,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,9,iblock),  &
                                 tavg_pCO2SURF,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,10,iblock),  &
                                 tavg_DpCO2,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,10,iblock),  &
                                 tavg_DpCO2_2,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,11,iblock),  &
                                 tavg_PV_CO2,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,12,iblock),  &
                                 tavg_SCHMIDT_CO2,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,13,iblock),  &
                                 tavg_DIC_GAS_FLUX,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,13,iblock),  &
                                 tavg_DIC_GAS_FLUX_2,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,14,iblock),  &
                                 tavg_PH,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,15,iblock),  &
                                 tavg_ATM_CO2,iblock,1)
      call accumulate_tavg_field(ECO_SFLUX_TAVG(:,:,16,iblock),  &
                                 tavg_O2_GAS_FLUX_2,iblock,1)

   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_tavg_forcing

!*****************************************************************************
!BOP
! !IROUTINE: comp_surf_avg
! !INTERFACE:

 subroutine comp_surf_avg(SURF_VALS_OLD,SURF_VALS_CUR)

! !DESCRIPTION:
!  compute average surface tracer values
!
!  avg = sum(SURF_VAL*TAREA) / sum(TAREA)
!  with the sum taken over ocean points only
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

 !real (r8), dimension(nx_block,ny_block,ecosys_tracer_cnt,max_blocks_clinic), &
  real (r8), dimension(:,:,:,:), &
      intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,       & ! tracer index
      iblock,  & ! block index
      dcount,  & ! diag counter
      ib,ie,jb,je

   real (r8), dimension(max_blocks_clinic,ecosys_tracer_cnt) :: &
      local_sums ! array for holding block sums of each diagnostic

   real (r8) :: &
      sum_tmp ! temp for local sum

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1, &! local work space
      TFACT   ! factor for normalizing sums

   type (block) :: &
      this_block ! block information for current block

!-----------------------------------------------------------------------

   local_sums = c0

!jw   !$OMP PARALLEL DO PRIVATE(iblock,this_block,ib,ie,jb,je,TFACT,n,WORK1)
   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      ib = this_block%ib
      ie = this_block%ie
      jb = this_block%jb
      je = this_block%je
      TFACT = TAREA(:,:,iblock)*RCALCT(:,:,iblock)

      do n = 1, ecosys_tracer_cnt
         if (vflux_flag(n)) then
            WORK1 = p5*(SURF_VALS_OLD(:,:,n,iblock) + &
                        SURF_VALS_CUR(:,:,n,iblock))*TFACT
            local_sums(iblock,n) = sum(WORK1(ib:ie,jb:je))
         endif
      end do
   end do
!jw   !$OMP END PARALLEL DO

   do n = 1, ecosys_tracer_cnt
      if (vflux_flag(n)) then
         sum_tmp = sum(local_sums(:,n))
         surf_avg(n) = global_sum(sum_tmp,distrb_clinic)/area_t
      endif
   end do

   if(my_task == master_task) then
      write(stdout,*)' Calculating surface tracer averages'
      do n = 1, ecosys_tracer_cnt
         if (vflux_flag(n)) then
            write(stdout,*) n, surf_avg(n)
         endif
      end do
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_surf_avg

!*****************************************************************************
!BOP
! !IROUTINE: ecosys_write_restart
! !INTERFACE:

 subroutine ecosys_write_restart(restart_file, action)

! !DESCRIPTION:
!  write auxiliary fields & scalars to restart files
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character(*), intent(in) :: action

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: restart_file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character (char_len) :: &
      short_name  ! tracer name temporaries

   type (io_dim) :: &
      i_dim, j_dim ! dimension descriptors

   integer (int_kind) :: n

   type (io_field_desc), save :: PH_SURF

!-----------------------------------------------------------------------

   if (trim(action) == 'add_attrib_file') then
      short_name = char_blank
      do n=1,ecosys_tracer_cnt
         if (vflux_flag(n)) then
            short_name = 'surf_avg_' /&
                      &/ ind_name_table(n)%name
            call add_attrib_file(restart_file,trim(short_name),surf_avg(n))
         endif
      end do
   endif

   if (trim(action) == 'define') then
      i_dim = construct_io_dim('i', nx_global)
      j_dim = construct_io_dim('j', ny_global)

      PH_SURF = construct_io_field('PH_SURF', i_dim, j_dim,     &
                   long_name='surface pH at current time',      &
                   units    ='pH', grid_loc ='2110',            &
                   field_loc = field_loc_center,                &
                   field_type = field_type_scalar,              &
                   d2d_array = PH_PREV)
      call data_set (restart_file, 'define', PH_SURF)
   endif

   if (trim(action) == 'write') then
      call data_set (restart_file, 'write', PH_SURF)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ecosys_write_restart

!*****************************************************************************
!BOP
! !IROUTINE: ecosys_tracer_ref_val
! !INTERFACE:

 function ecosys_tracer_ref_val(ind)

! !DESCRIPTION:
!  return reference value for tracers using virtual fluxes
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: ind

! !OUTPUT PARAMETERS:

   real (r8) :: ecosys_tracer_ref_val

!EOP
!BOC
!-----------------------------------------------------------------------

   if (vflux_flag(ind)) then
      ecosys_tracer_ref_val = surf_avg(ind)
   else
      ecosys_tracer_ref_val = c0
   endif

!-----------------------------------------------------------------------
!EOC

 end function ecosys_tracer_ref_val

!***********************************************************************

 end module ecosys_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
