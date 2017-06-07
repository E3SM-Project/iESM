MODULE ecosys_parms

  !-----------------------------------------------------------------------------
  !   This module manages the parameter variables for the module ecosys_mod.
  !   Most of the variables are not parameters in the Fortran sense. In the
  !   the Fortran sense, they are vanilla module variables.
  !
  !   This modules handles initializing the variables to default values and
  !   reading them from the namelist ecosys_parms. The values used are echoed
  !   to stdout for record keeping purposes.
  !
  !   CVS:$Id: ecosys_parms.F90 941 2006-05-12 21:36:48Z klindsay $
  !   CVS:$Name$
  !-----------------------------------------------------------------------------
  !   Modified to include parameters for diazotrophs, JKM  4/2002
  !-----------------------------------------------------------------------------
  !   variables/subroutines/function used from other modules
  !   The following are used extensively in this ecosys, so are used at
  !   the module level. The use statements for variables that are only needed
  !   locally are located at the module subprogram level.
  !-----------------------------------------------------------------------------

  USE exit_mod, ONLY : sigAbort, exit_POP
  USE communicate, ONLY : my_task, master_task
  USE constants, ONLY : c1
  USE kinds_mod
  USE io_tools, ONLY : document

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  !   public/private declarations
  !   all module variables are public and should have their values preserved
  !-----------------------------------------------------------------------------

  PUBLIC
  SAVE

  !-----------------------------------------------------------------------------
  !   floating point constants used across ecosystem module
  !-----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       spd = 86400.0_r8,    & ! number of seconds in a day
       dps = c1 / spd,            & ! number of days in a second
       yps = c1 / (365.0_r8*spd)! number of years in a second

  !-----------------------------------------------------------------------------
  !   Redfield Ratios, dissolved & particulate
  !-----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       parm_Red_D_C_P  = 117.0_r8,                 & ! carbon:phosphorus
       parm_Red_D_N_P  =  16.0_r8,                 & ! nitrogen:phosphorus
       parm_Red_D_O2_P = 170.0_r8,                 & ! oxygen:phosphorus
       parm_Remin_D_O2_P = 138.0_r8,               & ! oxygen:phosphorus
       parm_Red_P_C_P  = parm_Red_D_C_P,                 & ! carbon:phosphorus
       parm_Red_D_C_N  = parm_Red_D_C_P/parm_Red_D_N_P,  & ! carbon:nitrogen
       parm_Red_P_C_N  = parm_Red_D_C_N,                 & ! carbon:nitrogen
       parm_Red_D_C_O2 = parm_Red_D_C_P/parm_Red_D_O2_P, & ! carbon:oxygen
       parm_Remin_D_C_O2 = parm_Red_D_C_P/parm_Remin_D_O2_P, & ! carbon:oxygen
       parm_Red_P_C_O2 = parm_Red_D_C_O2,                & ! carbon:oxygen
       parm_Red_Fe_C   = 3.0e-6_r8,                & ! iron:carbon
       parm_Red_D_C_O2_diaz = parm_Red_D_C_P/150.0_r8! carbon:oxygen
                                                           ! for diazotrophs

  !----------------------------------------------------------------------------
  !   ecosystem parameters accessible via input file
  !----------------------------------------------------------------------------

  REAL(KIND=r8) :: &
       parm_Fe_bioavail,      & ! fraction of Fe flux that is bioavailable
       parm_o2_min,           & ! min O2 needed for prod & consump. (nmol/cm^3)
       parm_o2_min_delta,     & ! width of min O2 range (nmol/cm^3)
       parm_no3_min,          & ! min NO3 needed for denitrification (mmol/m^3)
       parm_kappa_nitrif,     & ! nitrification inverse time constant (1/sec)
       parm_nitrif_par_lim,   & ! PAR limit for nitrif. (W/m^2)
       parm_z_umax_0,         & ! max. zoo growth rate on sphyto at tref (1/sec)
       parm_diat_umax_0,      & ! max. zoo growth rate on diat at tref (1/sec)
       parm_z_mort_0,         & ! zoo linear mort rate (1/sec)
       parm_z_mort2_0,        & ! zoo quad mort rate (1/sec/((mmol C/m3))
       parm_sd_remin_0,       & ! small detrital remineralization rate (1/sec)
       parm_sp_kNO3,          & ! sphyto nitrate half sat. coef. (mmol N/m3)
       parm_diat_kNO3,        & ! diatom nitrate half sat. coef. (mmol N/m3)
       parm_sp_kNH4,          & ! sphyto ammonium half sat. coef. (mmol N/m3)
       parm_diat_kNH4,        & ! diatom ammonium half sat. coef. (mmol N/m3)
       parm_sp_kFe,           & ! sphyto iron half sat. coef. (nmol Fe/m3)
       parm_diat_kFe,         & ! diatom iron half sat. coef. (nmol Fe/m3)
       parm_diat_kSiO3,       & ! diatom si half sat. coef. (mmol SiO3/m3)
       parm_sp_kPO4,          & ! sphyto PO4 uptake (mmol P/m^3)
       parm_diat_kPO4,        & ! diatom PO4 uptate (mmol P/m^3)
       parm_z_grz,            & ! grazing coef. for small phyto (mmol C/m^3)
       parm_alphaChl,         & ! init slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
       parm_alphaChlsp,       & ! init slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
       parm_labile_ratio,     & ! fraction of loss to DOC that routed directly to DIC (non-dimensional)
       parm_alphaDiaz,        & ! chl. spec. init. slope of diaz. P_I curve
       parm_diaz_umax_0         ! max. zoo growth rate on diazotrophs at tref (1/sec)

  real(kind=r8), parameter ::      &
      PCref = 4.8_r8 * dps, & !max phyto C-spec. grth rate at tref (1/sec)
      sp_mort    = 0.15_r8  * dps, & !sphyto mort rate (1/sec)
      sp_mort2   = 0.0035_r8 * dps, & !sphyto quad. mort rate (1/sec/((mmol C/m3))
      diat_mort  = 0.15_r8   * dps, & !diatom mort rate (1/sec)
      diat_mort2 = 0.0035_r8 * dps, & !diatom quad mort rate (1/sec/((mmol C/m3))
      PCrefDiaz  = 0.70_r8  * dps,  & !max Diaz C-specific growth rate at tref (1/sec)
      diaz_mort  = 0.17_r8 * dps,  & !diaz mort rate (1/sec)
      diaz_kPO4  = 0.02_r8,      & !diaz half-sat. const. for P (diatom value)
      diaz_kFe   = 0.06e-3_r8       !diaz half-sat. const. for Fe

  !---------------------------------------------------------------------
  !     Misc. Rate constants
  !---------------------------------------------------------------------
  real(kind=r8), parameter :: &
       sp_agg_rate_max   = 0.75_r8, & !max agg. rate for small phyto (1/d)
       diat_agg_rate_max = 0.75_r8, & !max agg. rate for diatoms (1/d)
       diat_agg_rate_min = 0.01_r8,& !min agg. rate for diatoms (1/d)
       fe_scavenge_rate0 = 1.5_r8, & !base scavenging rate
       fe_scavenge_thres1 = 0.6e-3_r8,  & !upper thres. for Fe scavenging
       dust_fescav_scale  = 1.0e9,      & !dust scavenging scale factor
       fe_max_scale2      = 1000.0_r8,&   !unitless scaling coeff.
       f_fescav_P_iron    = 0.9_r8        !fraction of Fe scavenging 
                                          !        to particulate Fe

  !---------------------------------------------------------------------
  !     Compute iron remineralization and flux out.
  !     dust remin gDust = 0.035 gFe      mol Fe     1e9 nmolFe
  !                        --------- *  ---------- * ----------
  !			    gDust       55.847 gFe     molFe
  !
  !     dust_to_Fe          conversion - dust to iron (nmol Fe/g Dust) 
  !---------------------------------------------------------------------
  real(kind=r8), parameter :: &
       dust_to_Fe=0.035_r8/55.847_r8*1.0e9_r8
 
  !----------------------------------------------------------------------------
  !     Partitioning of phytoplankton growth, grazing and losses
  !
  !     All f_* variables are fractions and are non-dimensional
  !----------------------------------------------------------------------------

  real(kind=r8), parameter ::     &
      z_ingest         = 0.3_r8,  & !zoo ingestion coefficient (non-dim)
      caco3_poc_min    = 0.4_r8,  & !minimum proportionality between 
                                          !   QCaCO3 and grazing losses to POC 
                                          !   (mmol C/mmol CaCO3)
      spc_poc_fac      = 0.2_r8, & !small phyto grazing factor (1/mmolC)
      f_graze_sp_poc_lim = 0.25_r8, & 
      f_prod_sp_CaCO3  = 0.042_r8, & !fraction of sp prod. as CaCO3 prod.
      f_photosp_CaCO3  = 0.4_r8,  & !proportionality between small phyto 
                                          !    production and CaCO3 production
      f_graze_sp_doc   = 0.4_r8, & !fraction sm. phyto. grazing to DOC
      f_graze_sp_dic   = c1 - z_ingest - f_graze_sp_doc, & !fraction to DIC
      f_graze_diat_poc = 0.25_r8, & !fraction diatom grazing to POC
      f_graze_diat_doc = 0.35_r8, & !fraction diatom grazing to DOC
      f_graze_diat_dic = c1 - z_ingest - f_graze_diat_poc &
                          - f_graze_diat_doc, & !fraction diatom grazing to DIC
      f_diat_loss_poc  = 0.01_r8, &  !fraction diatom loss to POC
      f_diat_loss_dc   = c1-f_diat_loss_poc, & !fraction diatom loss to DOC
      f_graze_diaz_zoo = 0.21_r8, & !fraction diaz. grazing to zoo
      f_graze_diaz_poc = 0.0_r8, &  !fraction diaz grazing to POC
      f_graze_diaz_doc = 0.35_r8, & !fraction diaz grazing to DOC
      f_graze_diaz_dic = c1-f_graze_diaz_zoo-f_graze_diaz_poc &
                         - f_graze_diaz_doc, & !fraction diaz grazing to DIC
      f_diaz_loss_poc = 0.0_r8,& !fraction diaz loss to sinking pool
      f_sp_zoo_detr   = 0.1_r8,& !fraction of zoo losses to detrital 
                                          !  pool when eating sphyto
      f_diat_zoo_detr = 0.25_r8,& !fraction of zoo losses to detrital 
                                          !  pool when eating diatoms
      f_diaz_zoo_detr = 0.1_r8,& !fraction of zoo losses to detrital 
                                          !  pool when eating diaz
      f_graze_CaCO3_remin = 0.33_r8, & !fraction of spCaCO3 grazing 
                                             !          which is remin
      f_graze_si_remin    = 0.33_r8      !fraction of diatom Si grazing 
                                             !          which is remin

  !----------------------------------------------------------------------------
  !     fixed ratios
  !----------------------------------------------------------------------------
  real(kind=r8), parameter :: &
       r_Nfix_photo=1.3335_r8         ! N fix relative to C fix (non-dim)

  !-----------------------------------------------------------------------
  !     SET FIXED RATIOS for N/C, P/C, SiO3/C, Fe/C
  !     assumes C/N/P of 117/16/1 based on Anderson and Sarmiento, 1994
  !     for diazotrophs a N/P of 45 is assumed based on Letelier & Karl, 1998
  !-----------------------------------------------------------------------

  real(kind=r8), parameter ::  &
      Q             = 0.137_r8,  & !N/C ratio (mmol/mmol) of phyto & zoo
      Qp            = 0.00855_r8,& !P/C ratio (mmol/mmol) sphyto,diat,zoo
      Qp_diaz       = 0.002735_r8,& !diazotroph P/C ratio
      Qfe_zoo       = 3.0e-6_r8, & !zooplankton fe/C ratio
      gQsi_0        = 0.137_r8,  & !initial diatom Si/C ratio
      gQfe_diat_0   = 6.0e-6_r8, & !initial diatom fe/C ratio
      gQfe_sp_0     = 6.0e-6_r8, & !initial sphyto fe/C ratio
      gQfe_diaz_0   = 42.0e-6_r8,& !initial diaz. fe/C ratio
      gQfe_diat_min = 3.0e-6_r8, & !min diatom fe/C ratio
      gQsi_max      = 0.685_r8,  & !max diatom Si/C ratio
      gQsi_min      = 0.0685_r8, & !min diatom Si/C ratio
      gQfe_sp_min   = 3.0e-6_r8, & !min sphyto fe/C ratio
      gQfe_diaz_min = 14.0e-6_r8,& !min diaz fe/C ratio
      QCaCO3_max    = 0.4_r8,    & !max QCaCO3
      thetaN_max_sp   = 2.5_r8, & !sp max thetaN (Chl/N) (mg Chl/mmol N)
      thetaN_max_diat = 4.0_r8, & !diat max thetaN (Chl/N) (mg Chl/mmol N)
      thetaN_max_diaz = 2.5_r8, & !diaz max thetaN (Chl/N) (mg Chl/mmol N)
      ! carbon:nitrogen ratio for denitrification
      ! net removal of 120 mols NO3 for 117 mols C (136 = 120 + 16) 
      denitrif_C_N  = parm_Red_D_C_P/136.0_r8

  !----------------------------------------------------------------------------
  !     loss term threshold parameters, chl:c ratios
  !----------------------------------------------------------------------------

  real(kind=r8), parameter ::    &
      thres_z1          = 100.0e2_r8, & !threshold = C_loss_thres for z shallower than this (cm)
      thres_z2          = 200.0e2_r8, & !threshold = 0 for z deeper than this (cm)
      loss_thres_sp     = 0.001_r8, & !small phyto conc. where losses go to zero 
      loss_thres_diat   = 0.02_r8, & !diat conc. where losses go to zero
      loss_thres_zoo    = 0.2_r8,  & !zoo conc. where losses go to zero
      loss_thres_diaz   = 0.01_r8,  & !diaz conc. where losses go to zero
      loss_thres_diaz2  = 0.001_r8, & !diaz conc. thres at low temp
      diaz_temp_thres   = 17.0_r8,  & !Temp. where diaz conc thres drops
      CaCO3_temp_thres1 = 2.0_r8,   & !upper temp threshold for CaCO3 prod
      CaCO3_temp_thres2 = -2.0_r8,  & !lower temp threshold
      CaCO3_sp_thres    = 3.0_r8,   & ! bloom condition thres (mmolC/m3)
      diaz_kNO3         = 1.0_r8,   & ! diazotroph Ks for nitrate (mmolN/m3)
      diaz_kNH4         = 0.1_r8      ! diazotroph Ks for ammonimum (mmolN/m3)

  !---------------------------------------------------------------------
  !     attenuation coefficients for PAR and related parameters
  !---------------------------------------------------------------------
  real(kind=r8), parameter :: &
       k_chl = 0.03e-2_r8, & ! Chl atten. coeff. (1/cm/(mg Chl/m^3))
       k_h2o = 0.04e-2_r8, & ! water atten. coeff (1/cm)
       f_qsw_par = 0.45_r8   ! PAR fraction

  !---------------------------------------------------------------------
  !     Temperature parameters
  !---------------------------------------------------------------------
  real(kind=r8), parameter :: &
       Tref = 30.0_r8, & ! reference temperature (C)
       Q_10 = 2.0_r8     ! factor for temperature dependence (non-dim)

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  SUBROUTINE ecosys_parms_init

    USE io_types, ONLY: stdout, nml_in, nml_filename
    USE broadcast, ONLY : broadcast_scalar

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: subname = 'ecosys_parms:ecosys_parms_init'

    LOGICAL(KIND=log_kind) :: &
         lnml_found             ! Was ecosys_parms_nml found ?

    NAMELIST /ecosys_parms_nml/ &
         parm_Fe_bioavail, &
         parm_o2_min, &
         parm_o2_min_delta, &
         parm_no3_min, &
         parm_kappa_nitrif, &
         parm_nitrif_par_lim, &
         parm_z_umax_0, &
         parm_diat_umax_0, &
         parm_z_mort_0, &
         parm_z_mort2_0, &
         parm_sd_remin_0, &
         parm_sp_kNO3, &
         parm_diat_kNO3, &
         parm_sp_kNH4, &
         parm_diat_kNH4, &
         parm_sp_kFe, &
         parm_diat_kFe, &
         parm_diat_kSiO3, &
         parm_sp_kPO4, &
         parm_diat_kPO4, &
         parm_z_grz, &
         parm_alphaChl, &
         parm_alphaChlsp, &
         parm_labile_ratio, &
         parm_alphaDiaz, &
         parm_diaz_umax_0

    !---------------------------------------------------------------------------
    !   default namelist settings
    !---------------------------------------------------------------------------

    parm_Fe_bioavail    = 0.02_r8
    parm_o2_min         = 4.0_r8
    parm_o2_min_delta   = 2.0_r8
    parm_no3_min        = 110.0_r8
    parm_kappa_nitrif   = 0.06_r8 * dps       ! (= 1/( days))
    parm_nitrif_par_lim = 5.0_r8
    parm_z_umax_0       = 2.5_r8 * dps
    parm_diat_umax_0    = 1.95_r8 * dps     
    parm_z_mort_0       = 0.08_r8 * dps
    parm_z_mort2_0      = 0.42_r8 * dps
    parm_sd_remin_0     = 0.006667_r8 * dps       ! (= 1/(100 days))
    parm_sp_kNO3        = 0.5_r8
    parm_diat_kNO3      = 2.5_r8
    parm_sp_kNH4        = 0.01_r8
    parm_diat_kNH4      = 0.1_r8
    parm_sp_kFe         = 0.03e-3_r8
    parm_diat_kFe       = 0.08e-3_r8
    parm_diat_kSiO3     = 1.0_r8
    parm_sp_kPO4        = 0.01_r8
    parm_diat_kPO4      = 0.1_r8
    parm_z_grz          = 1.0_r8              
    parm_alphaChl       = 0.3_r8 * dps
    parm_alphaChlsp     = 0.34_r8 * dps
    parm_labile_ratio   = 0.65_r8
    parm_alphaDiaz      = 0.17_r8 * dps
    parm_diaz_umax_0    = 0.9_r8 * dps

    !---------------------------------------------------------------------------
    !   read in namelist
    !---------------------------------------------------------------------------

    IF (my_task == master_task) THEN
       lnml_found = .FALSE.
       OPEN(UNIT=nml_in, FILE=nml_filename, STATUS='OLD')
10     CONTINUE
       READ(UNIT=nml_in, NML=ecosys_parms_nml, ERR=10, END=20)
       CLOSE(UNIT=nml_in)
       lnml_found = .TRUE.
20     CONTINUE
    END IF

    CALL broadcast_scalar(lnml_found, master_task)
    IF (.NOT. lnml_found) THEN
       CALL document(subname, 'ecosys_parms_nml not found')
       CALL exit_POP(sigAbort, 'ERROR : stopping in ' // subname)
    END IF

    !---------------------------------------------------------------------------
    !   broadcast all namelist variables
    !---------------------------------------------------------------------------

    CALL broadcast_scalar(parm_Fe_bioavail, master_task)
    CALL broadcast_scalar(parm_o2_min, master_task)
    CALL broadcast_scalar(parm_o2_min_delta, master_task)
    CALL broadcast_scalar(parm_no3_min, master_task)
    CALL broadcast_scalar(parm_kappa_nitrif, master_task)
    CALL broadcast_scalar(parm_nitrif_par_lim, master_task)
    CALL broadcast_scalar(parm_z_umax_0, master_task)
    CALL broadcast_scalar(parm_diat_umax_0, master_task)
    CALL broadcast_scalar(parm_z_mort_0, master_task)
    CALL broadcast_scalar(parm_z_mort2_0, master_task)
    CALL broadcast_scalar(parm_sd_remin_0, master_task)
    CALL broadcast_scalar(parm_sp_kNO3, master_task)
    CALL broadcast_scalar(parm_diat_kNO3, master_task)
    CALL broadcast_scalar(parm_sp_kNH4, master_task)
    CALL broadcast_scalar(parm_diat_kNH4, master_task)
    CALL broadcast_scalar(parm_sp_kFe, master_task)
    CALL broadcast_scalar(parm_diat_kFe, master_task)
    CALL broadcast_scalar(parm_diat_kSiO3, master_task)
    CALL broadcast_scalar(parm_sp_kPO4, master_task)
    CALL broadcast_scalar(parm_diat_kPO4, master_task)
    CALL broadcast_scalar(parm_z_grz, master_task)
    CALL broadcast_scalar(parm_alphaChl, master_task)
    CALL broadcast_scalar(parm_alphaChlsp, master_task)
    CALL broadcast_scalar(parm_labile_ratio, master_task)
    CALL broadcast_scalar(parm_alphaDiaz, master_task)
    CALL broadcast_scalar(parm_diaz_umax_0, master_task)

    !---------------------------------------------------------------------------
    !   echo all namelist variables to stdout
    !---------------------------------------------------------------------------

    IF (my_task == master_task) THEN
       WRITE (stdout,*) '----------------------------------------'
       WRITE (stdout,*) '----- ecosys_parms namelist values -----'
       WRITE (stdout,*) 'parm_Fe_bioavail    = ', parm_Fe_bioavail
       WRITE (stdout,*) 'parm_o2_min         = ', parm_o2_min
       WRITE (stdout,*) 'parm_o2_min_delta   = ', parm_o2_min_delta
       WRITE (stdout,*) 'parm_no3_min        = ', parm_no3_min
       WRITE (stdout,*) 'parm_kappa_nitrif   = ', parm_kappa_nitrif
       WRITE (stdout,*) 'parm_nitrif_par_lim = ', parm_nitrif_par_lim
       WRITE (stdout,*) 'parm_z_umax_0       = ', parm_z_umax_0
       WRITE (stdout,*) 'parm_diat_umax_0    = ', parm_diat_umax_0
       WRITE (stdout,*) 'parm_z_mort_0       = ', parm_z_mort_0
       WRITE (stdout,*) 'parm_z_mort2_0      = ', parm_z_mort2_0
       WRITE (stdout,*) 'parm_sd_remin_0     = ', parm_sd_remin_0
       WRITE (stdout,*) 'parm_sp_kNO3        = ', parm_sp_kNO3
       WRITE (stdout,*) 'parm_diat_kNO3      = ', parm_diat_kNO3
       WRITE (stdout,*) 'parm_sp_kNH4        = ', parm_sp_kNH4
       WRITE (stdout,*) 'parm_diat_kNH4      = ', parm_diat_kNH4
       WRITE (stdout,*) 'parm_sp_kFe         = ', parm_sp_kFe
       WRITE (stdout,*) 'parm_diat_kFe       = ', parm_diat_kFe
       WRITE (stdout,*) 'parm_diat_kSiO3     = ', parm_diat_kSiO3
       WRITE (stdout,*) 'parm_sp_kPO4        = ', parm_sp_kPO4
       WRITE (stdout,*) 'parm_diat_kPO4      = ', parm_diat_kPO4
       WRITE (stdout,*) 'parm_z_grz          = ', parm_z_grz
       WRITE (stdout,*) 'parm_alphaChl       = ', parm_alphaChl
       WRITE (stdout,*) 'parm_alphaChlsp     = ', parm_alphaChlsp
       WRITE (stdout,*) 'parm_labile_ratio   = ', parm_labile_ratio
       WRITE (stdout,*) 'parm_alphaDiaz      = ', parm_alphaDiaz
       WRITE (stdout,*) 'parm_diaz_umax_0    = ', parm_diaz_umax_0
       WRITE (stdout,*) '----------------------------------------'
    END IF

  END SUBROUTINE ecosys_parms_init

  !*****************************************************************************

END MODULE ecosys_parms
