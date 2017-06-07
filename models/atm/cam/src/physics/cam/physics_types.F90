!-------------------------------------------------------------------------------
!physics data types module
!-------------------------------------------------------------------------------
module physics_types

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver
  use constituents, only: pcnst, qmin, cnst_name
  use geopotential, only: geopotential_dse
  use physconst,    only: zvir, gravit, cpair, rair, cpairv, rairv
  use dycore,       only: dycore_is
  use phys_grid,    only: get_ncols_p, get_rlon_all_p, get_rlat_all_p, get_gcol_all_p
  use cam_logfile,  only: iulog
  use abortutils,   only: endrun
  use phys_control, only: waccmx_is
  use shr_const_mod,only: shr_const_rwv

  implicit none
  private          ! Make default type private to the module

  logical, parameter :: adjust_te = .FALSE.

! Public types:

  public physics_state
  public physics_tend
  public physics_ptend
  
! Public interfaces

  public physics_update
  public physics_ptend_reset
  public physics_ptend_init
  public physics_state_set_grid
  public physics_dme_adjust  ! adjust dry mass and energy for change in water
                             ! cannot be applied to eul or sld dycores
  public physics_state_copy  ! copy a physics_state object
  public physics_ptend_copy  ! copy a physics_ptend object
  public physics_ptend_sum   ! accumulate physics_ptend objects
  public physics_tend_init   ! initialize a physics_tend object

  public set_state_pdry      ! calculate dry air masses in state variable
  public set_wet_to_dry
  public set_dry_to_wet
  public physics_type_alloc
!-------------------------------------------------------------------------------
  type physics_state
     integer                                     :: &
          lchnk,   &! chunk index
          ncol      ! number of active columns
     real(r8), dimension(pcols)                  :: &
          lat,     &! latitude (radians)
          lon,     &! longitude (radians)
          ps,      &! surface pressure
          psdry,   &! dry surface pressure
          phis,    &! surface geopotential
          ulat,    &! unique latitudes  (radians)
          ulon      ! unique longitudes (radians)
     real(r8), dimension(pcols,pver)             :: &
          t,       &! temperature (K)
          u,       &! zonal wind (m/s)
          v,       &! meridional wind (m/s)
          s,       &! dry static energy
          omega,   &! vertical pressure velocity (Pa/s) 
          pmid,    &! midpoint pressure (Pa) 
          pmiddry, &! midpoint pressure dry (Pa) 
          pdel,    &! layer thickness (Pa)
          pdeldry, &! layer thickness dry (Pa)
          rpdel,   &! reciprocal of layer thickness (Pa)
          rpdeldry,&! recipricol layer thickness dry (Pa)
          lnpmid,  &! ln(pmid)
          lnpmiddry,&! log midpoint pressure dry (Pa) 
          exner,   &! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
          uzm,     &! zonal wind for qbo (m/s)
          zm        ! geopotential height above surface at midpoints (m)

     real(r8), dimension(pcols,pver,pcnst) :: &
          q         ! constituent mixing ratio (kg/kg moist or dry air depending on type)

     real(r8), dimension(pcols,pver+1)           :: &
          pint,    &! interface pressure (Pa)
          pintdry, &! interface pressure dry (Pa) 
          lnpint,  &! ln(pint)
          lnpintdry,&! log interface pressure dry (Pa) 
          zi        ! geopotential height above surface at interfaces (m)

     real(r8), dimension(pcols)                  :: &
          te_ini,  &! vertically integrated total (kinetic + static) energy of initial state
          te_cur,  &! vertically integrated total (kinetic + static) energy of current state
          tw_ini,  &! vertically integrated total water of initial state
          tw_cur    ! vertically integrated total water of new state
     integer :: count ! count of values with significant energy or water imbalances
     integer, dimension(pcols) :: &
          latmapback, &! map from column to unique lat for that column
          lonmapback, &! map from column to unique lon for that column
          cid        ! unique column id
     integer :: ulatcnt, &! number of unique lats in chunk
                uloncnt   ! number of unique lons in chunk

#if (defined WACCM_PHYS)
     real(r8), dimension(pcols,pver)             :: &
          frontgf, &  ! frontogenesis function
          frontga     ! frontogenesis angle
#endif
  end type physics_state

!-------------------------------------------------------------------------------
  type physics_tend
     real(r8), dimension(pcols,pver)             :: dtdt, dudt, dvdt
     real(r8), dimension(pcols     )             :: flx_net
     real(r8), dimension(pcols)                  :: &
          te_tnd,  &! cumulative boundary flux of total energy
          tw_tnd    ! cumulative boundary flux of total water
  end type physics_tend

!-------------------------------------------------------------------------------
! This is for tendencies returned from individual parameterizations
  type physics_ptend
     character*24 :: name    ! name of parameterization which produced tendencies.

     logical ::             &
          ls,               &! true if dsdt is returned
          lu,               &! true if dudt is returned
          lv,               &! true if dvdt is returned
          lq(pcnst)          ! true if dqdt() is returned

     integer ::             &
          top_level,        &! top level index for which nonzero tendencies have been set
          bot_level          ! bottom level index for which nonzero tendencies have been set

     real(r8), dimension(pcols,pver)             :: &
          s,                &! heating rate (J/kg/s)
          u,                &! u momentum tendency (m/s/s)
          v                  ! v momentum tendency (m/s/s)
     real(r8), dimension(pcols,pver,pcnst) :: &
          q                  ! consituent tendencies (kg/kg/s)

! boundary fluxes
     real(r8), dimension(pcols) ::&
          hflux_srf,     &! net heat flux at surface (W/m2)
          hflux_top,     &! net heat flux at top of model (W/m2)
          taux_srf,      &! net zonal stress at surface (Pa)
          taux_top,      &! net zonal stress at top of model (Pa)
          tauy_srf,      &! net meridional stress at surface (Pa)
          tauy_top        ! net meridional stress at top of model (Pa)
     real(r8), dimension(pcols,pcnst) ::&
          cflx_srf,      &! constituent flux at surface (kg/m2/s)
          cflx_top        ! constituent flux top of model (kg/m2/s)

  end type physics_ptend


!===============================================================================
contains
!===============================================================================
  subroutine physics_type_alloc(phys_state, phys_tend, begchunk, endchunk)
    use infnan, only : inf
    implicit none
    type(physics_state), pointer :: phys_state(:)
    type(physics_tend), pointer :: phys_tend(:)
    integer, intent(in) :: begchunk, endchunk
    
    integer :: ierr, lchnk
    type(physics_state), pointer :: state
    type(physics_tend), pointer :: tend

    allocate(phys_state(begchunk:endchunk), stat=ierr)
    if( ierr /= 0 ) then
       write(iulog,*) 'physics_types: phys_state allocation error = ',ierr
       call endrun('physics_types: failed to allocate physics_state array')
    end if

    allocate(phys_tend(begchunk:endchunk), stat=ierr)
    if( ierr /= 0 ) then
       write(iulog,*) 'physics_types: phys_tend allocation error = ',ierr
       call endrun('physics_types: failed to allocate physics_tend array')
    end if

   ! Set chunk id, number of columns, and coordinates
    do lchnk = begchunk,endchunk
       state => phys_state(lchnk)
       tend => phys_tend(lchnk)
       state%lat(:) = inf
       state%lon(:) = inf
       state%ulat(:) = inf
       state%ulon(:) = inf
       state%ps(:) = inf
       state%psdry(:) = inf
       state%phis(:) = inf
       state%t(:,:) = inf
       state%u(:,:) = inf
       state%v(:,:) = inf
       state%s(:,:) = inf
       state%omega(:,:) = inf
       state%pmid(:,:) = inf
       state%pmiddry(:,:) = inf
       state%pdel(:,:) = inf
       state%pdeldry(:,:) = inf
       state%rpdel(:,:) = inf
       state%rpdeldry(:,:) = inf
       state%lnpmid(:,:) = inf
       state%lnpmiddry(:,:) = inf
       state%exner(:,:) = inf
       state%uzm(:,:) = inf
       state%zm(:,:) = inf
       state%q(:,:,:) = inf
      
       state%pint(:,:) = inf
       state%pintdry(:,:) = inf
       state%lnpint(:,:) = inf
       state%lnpintdry(:,:) = inf
       state%zi(:,:) = inf
      
       state%te_ini(:) = inf
       state%te_cur(:) = inf
       state%tw_ini(:) = inf
       state%tw_cur(:) = inf

       tend%dtdt(:,:) = inf
       tend%dudt(:,:) = inf
       tend%dvdt(:,:) = inf
       tend%flx_net(:) = inf
       tend%te_tnd(:) = inf
       tend%tw_tnd(:) = inf
    end do

  end subroutine physics_type_alloc
!===============================================================================
  subroutine physics_update(state, tend, ptend, dt)
!-----------------------------------------------------------------------
! Update the state and or tendency structure with the parameterization tendencies
!-----------------------------------------------------------------------
    use geopotential, only: geopotential_dse
    use constituents, only: cnst_get_ind, cnst_mw
    use scamMod,      only: scm_crm_mode, single_column
    use phys_control, only: phys_getopts
    use physconst,    only: physconst_update ! Routine which updates physconst variables (WACCM-X)

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies

    type(physics_state), intent(inout)  :: state   ! Physics state variables
    type(physics_tend ), intent(inout)  :: tend    ! Physics tendencies

    real(r8), intent(in) :: dt                     ! time step
!
!---------------------------Local storage-------------------------------
    integer :: i,k,m                               ! column,level,constituent indices
    integer :: ixcldice, ixcldliq                  ! indices for CLDICE and CLDLIQ
    integer :: ixnumice, ixnumliq
    integer :: ncol                                ! number of columns
    character*40 :: name    ! param and tracer name for qneg3

    integer :: ixo, ixo2, ixh, ixh2, ixn    ! indices for O, O2, H2, and N

    real(r8) :: zvirv(pcols,pver)                  ! Local zvir array pointer

    character(len=16) :: microp_scheme='RK'  ! microphysics scheme    
    !-----------------------------------------------------------------------

    ! The column radiation model does not update the state
    if(single_column.and.scm_crm_mode) return

    call phys_getopts(microp_scheme_out=microp_scheme)

    ncol = state%ncol

    ! Update u,v fields
    if(ptend%lu) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%u  (i,k) = state%u  (i,k) + ptend%u(i,k) * dt
             tend%dudt(i,k) = tend%dudt(i,k) + ptend%u(i,k)
          end do
       end do
    end if

    if(ptend%lv) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%v  (i,k) = state%v  (i,k) + ptend%v(i,k) * dt
             tend%dvdt(i,k) = tend%dvdt(i,k) + ptend%v(i,k)
          end do
       end do
    end if

   ! Update constituents, all schemes use time split q: no tendency kept
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    ! Check for number concentration of cloud liquid and cloud ice (if not present
    ! the indices will be set to -1)
    call cnst_get_ind('NUMICE', ixnumice, abort=.false.)
    call cnst_get_ind('NUMLIQ', ixnumliq, abort=.false.)

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
      call cnst_get_ind('H', ixh)
      call cnst_get_ind('H2', ixh2)
    endif
  
    do m = 1, pcnst
       if(ptend%lq(m)) then
          do k = ptend%top_level, ptend%bot_level
             do i = 1,ncol
                state%q(i,k,m) = state%q(i,k,m) + ptend%q(i,k,m) * dt
             end do
          end do

          ! now test for mixing ratios which are too small
          ! don't call qneg3 for number concentration variables
          if (m .ne. ixnumice  .and.  m .ne. ixnumliq) then
             name = trim(ptend%name) // '/' // trim(cnst_name(m))
             call qneg3(trim(name), state%lchnk, ncol, pcols, pver, m, m, qmin(m), state%q(1,1,m))
          else
             do k = ptend%top_level, ptend%bot_level
                do i = 1,ncol
                   ! checks for number concentration
                   state%q(i,k,m) = max(1.e-12_r8,state%q(i,k,m))
                   state%q(i,k,m) = min(1.e10_r8,state%q(i,k,m))
                end do
             end do
          end if

       end if

       !------------------------------------------------------------------------
       ! This is a temporary fix for the large H, H2 in WACCM-X
       !------------------------------------------------------------------------
       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
         if (m == ixh) then
           do k = ptend%top_level, ptend%bot_level
             do i = 1,ncol
                if (state%q(i,k,m) .gt. 0.01_r8) then
!                   write(*,*)'Large H in physics_update',i,k,m,state%q(i,k,m)
                   state%q(i,k,m) = 0.01_r8
                endif
             end do
           end do
         endif
         if (m == ixh2) then
           do k = ptend%top_level, ptend%bot_level
             do i = 1,ncol
                if (state%q(i,k,m) .gt. 6.e-5_r8) then
!                   write(*,*)'Large H2 in physics_update',i,k,m,state%q(i,k,m)
                   state%q(i,k,m) = 6.e-5_r8
                endif
             end do
           end do
         endif
       endif

    end do

    ! special tests for cloud liquid
    if (ixcldliq > 1) then
       if(ptend%lq(ixcldliq)) then
          if (ptend%name == 'stratiform' .or. ptend%name == 'cldwat'  ) then
#ifdef PERGRO
             where (state%q(:ncol,:pver,ixcldliq) < 1.e-12_r8)
                state%q(:ncol,:pver,ixcldliq) = 0._r8
             end where

             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' ) then
                where (state%q(:ncol,:pver,ixcldliq) < 1.e-12_r8)
                   state%q(:ncol,:pver,ixnumliq) = 0._r8
                end where
             end if
#endif
          else if (ptend%name == 'convect_deep') then
             where (state%q(:ncol,:pver,ixcldliq) < 1.e-36_r8)
                state%q(:ncol,:pver,ixcldliq) = 0._r8
             end where
             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' ) then
                where (state%q(:ncol,:pver,ixcldliq) < 1.e-36_r8)
                   state%q(:ncol,:pver,ixnumliq) = 0._r8
                end where
             end if
          end if
       end if
    end if

    ! special tests for cloud ice
    if (ixcldice > 1) then
       if(ptend%lq(ixcldice)) then
          if (ptend%name == 'stratiform' .or. ptend%name == 'cldwat'  ) then
#ifdef PERGRO
             where (state%q(:ncol,:pver,ixcldice) < 1.e-12_r8)
                state%q(:ncol,:pver,ixcldice) = 0._r8
             end where
             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' ) then
                where (state%q(:ncol,:pver,ixcldice) < 1.e-12_r8)
                   state%q(:ncol,:pver,ixnumice) = 0._r8
                end where
             end if
#endif
          else if (ptend%name == 'convect_deep') then
             where (state%q(:ncol,:pver,ixcldice) < 1.e-36_r8)
                state%q(:ncol,:pver,ixcldice) = 0._r8
             end where
             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' ) then
                where (state%q(:ncol,:pver,ixcldice) < 1.e-36_r8)
                   state%q(:ncol,:pver,ixnumice) = 0._r8
                end where
             end if
          end if
       end if
    end if

    !------------------------------------------------------------------------
    ! Get indices for molecular weights and call WACCM-X physconst_update
    !------------------------------------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
      call cnst_get_ind('O', ixo)
      call cnst_get_ind('O2', ixo2)
      call cnst_get_ind('N', ixn)             

      call physconst_update(state%q, state%t, &
	         cnst_mw(ixo), cnst_mw(ixo2), cnst_mw(ixh), cnst_mw(ixn), &
	                              ixo, ixo2, ixh, pcnst, state%lchnk, ncol)
    endif	  
   
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
      zvirv(:,:) = shr_const_rwv / rairv(:,:,state%lchnk) - 1._r8
    else
      zvirv(:,:) = zvir    
    endif

    !-------------------------------------------------------------------------------------------
    ! Update dry static energy(moved from above for WACCM-X so updating after cpairv update)
    !-------------------------------------------------------------------------------------------
    if(ptend%ls) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%s(i,k)   = state%s(i,k)   + ptend%s(i,k) * dt
             tend%dtdt(i,k) = tend%dtdt(i,k) + ptend%s(i,k)/cpairv(i,k,state%lchnk)
          end do
       end do
    end if

    ! Derive new temperature and geopotential fields if heating or water tendency not 0.
    if (ptend%ls .or. ptend%lq(1)) then
       call geopotential_dse(                                                                    &
            state%lnpint, state%lnpmid, state%pint  , state%pmid  , state%pdel  , state%rpdel  , &
            state%s     , state%q(1,1,1),state%phis , rairv(:,:,state%lchnk), gravit  , cpairv(:,:,state%lchnk), &
            zvirv    , state%t     , state%zi    , state%zm    , ncol         )
    end if

    ! Reset all parameterization tendency flags to false
    call physics_ptend_reset(ptend)

  end subroutine physics_update

!===============================================================================

  subroutine physics_ptend_sum(ptend, ptend_sum, state)
!-----------------------------------------------------------------------
! Add ptend fields to ptend_sum for ptend logical flags = .true.
! Where ptend logical flags = .false, don't change ptend_sum
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(in)     :: ptend   ! New parameterization tendencies
    type(physics_ptend), intent(inout)  :: ptend_sum   ! Sum of incoming ptend_sum and ptend
    type(physics_state), intent(in)     :: state   ! New parameterization tendencies

!---------------------------Local storage-------------------------------
    integer :: i,k,m                               ! column,level,constituent indices
    integer :: ncol                                ! number of columns

!-----------------------------------------------------------------------
    ncol = state%ncol
    

! Update u,v fields
    if(ptend%lu) then
       ptend_sum%lu = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%u(i,k) = ptend_sum%u(i,k) + ptend%u(i,k)
          end do
          ptend_sum%taux_srf(i) = ptend_sum%taux_srf(i) + ptend%taux_srf(i)
          ptend_sum%taux_top(i) = ptend_sum%taux_top(i) + ptend%taux_top(i)
       end do
    end if

    if(ptend%lv) then
       ptend_sum%lv = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%v(i,k) = ptend_sum%v(i,k) + ptend%v(i,k)
          end do
          ptend_sum%tauy_srf(i) = ptend_sum%tauy_srf(i) + ptend%tauy_srf(i)
          ptend_sum%tauy_top(i) = ptend_sum%tauy_top(i) + ptend%tauy_top(i)
       end do
    end if


    if(ptend%ls) then
       ptend_sum%ls = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%s(i,k) = ptend_sum%s(i,k) + ptend%s(i,k)
          end do
          ptend_sum%hflux_srf(i) = ptend_sum%hflux_srf(i) + ptend%hflux_srf(i)
          ptend_sum%hflux_top(i) = ptend_sum%hflux_top(i) + ptend%hflux_top(i)
       end do
    end if

! Update constituents
    do m = 1, pcnst
       if(ptend%lq(m)) then
          ptend_sum%lq(m) = .true.
          do i = 1,ncol
             do k = ptend%top_level, ptend%bot_level
                ptend_sum%q(i,k,m) = ptend_sum%q(i,k,m) + ptend%q(i,k,m)
             end do
             ptend_sum%cflx_srf(i,m) = ptend_sum%cflx_srf(i,m) + ptend%cflx_srf(i,m)
             ptend_sum%cflx_top(i,m) = ptend_sum%cflx_top(i,m) + ptend%cflx_top(i,m)
          end do
       end if
    end do


  end subroutine physics_ptend_sum

!===============================================================================

subroutine physics_ptend_copy(ptend, ptend_cp)

   !-----------------------------------------------------------------------
   ! Copy a physics_ptend object
   !-----------------------------------------------------------------------

   type(physics_ptend), intent(in)  :: ptend    ! ptend source
   type(physics_ptend), intent(out) :: ptend_cp ! copy of ptend

   !-----------------------------------------------------------------------

   ptend_cp%name  = ptend%name

   ptend_cp%ls = ptend%ls
   ptend_cp%lu = ptend%lu
   ptend_cp%lv = ptend%lv
   ptend_cp%lq = ptend%lq

   ptend_cp%top_level = ptend%top_level
   ptend_cp%bot_level = ptend%bot_level

   ptend_cp%s = ptend%s
   ptend_cp%u = ptend%u
   ptend_cp%v = ptend%v
   ptend_cp%q = ptend%q

   ptend_cp%hflux_srf = ptend%hflux_srf
   ptend_cp%hflux_top = ptend%hflux_top
   ptend_cp%taux_srf  = ptend%taux_srf
   ptend_cp%taux_top  = ptend%taux_top
   ptend_cp%tauy_srf  = ptend%tauy_srf
   ptend_cp%tauy_top  = ptend%tauy_top
   ptend_cp%cflx_srf  = ptend%cflx_srf
   ptend_cp%cflx_top  = ptend%cflx_top

end subroutine physics_ptend_copy

!===============================================================================

  subroutine physics_ptend_reset(ptend)
!-----------------------------------------------------------------------
! Reset the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
!-----------------------------------------------------------------------
    integer :: m             ! Index for constiuent
!-----------------------------------------------------------------------

    if(ptend%ls) then
       ptend%s = 0._r8
       ptend%hflux_srf = 0._r8
       ptend%hflux_top = 0._r8
    endif
    if(ptend%lu) then
       ptend%u = 0._r8
       ptend%taux_srf = 0._r8
       ptend%taux_top = 0._r8
    endif
    if(ptend%lv) then
       ptend%v = 0._r8
       ptend%tauy_srf = 0._r8
       ptend%tauy_top = 0._r8
    endif
    do m = 1, pcnst
       if(ptend%lq(m)) then
          ptend%q(:,:,m) = 0._r8
          ptend%cflx_srf(:,m) = 0._r8
          ptend%cflx_top(:,m) = 0._r8
       endif
    end do

    ptend%name  = "none"
    ptend%lq(:) = .FALSE.
    ptend%ls    = .FALSE.
    ptend%lu    = .FALSE.
    ptend%lv    = .FALSE.

    ptend%top_level = 1
    ptend%bot_level = pver

    return
  end subroutine physics_ptend_reset

!===============================================================================
  subroutine physics_ptend_init(ptend)
!-----------------------------------------------------------------------
! Initialize the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
!-----------------------------------------------------------------------
    ptend%name  = "none"
    ptend%lq(:) = .true.
    ptend%ls    = .true.
    ptend%lu    = .true.
    ptend%lv    = .true.

    call physics_ptend_reset(ptend)

    return
  end subroutine physics_ptend_init

!===============================================================================

  subroutine physics_state_set_grid(lchnk, phys_state)
!-----------------------------------------------------------------------
! Set the grid components of the physics_state object
!-----------------------------------------------------------------------

    integer,             intent(in)    :: lchnk
    type(physics_state), intent(inout) :: phys_state

    ! local variables
    integer  :: i, ncol
    real(r8) :: rlon(pcols)
    real(r8) :: rlat(pcols)
    !-----------------------------------------------------------------------

    ncol = get_ncols_p(lchnk)

    if(ncol<=0) then
       write(iulog,*) lchnk, ncol
       call endrun('physics_state_set_grid')
    end if

    call get_rlon_all_p(lchnk, ncol, rlon)
    call get_rlat_all_p(lchnk, ncol, rlat)
    phys_state%ncol  = ncol
    phys_state%lchnk = lchnk
    do i=1,ncol
       phys_state%lat(i) = rlat(i)
       phys_state%lon(i) = rlon(i)
    end do
    call init_geo_unique(phys_state,ncol)

  end subroutine physics_state_set_grid


  subroutine init_geo_unique(phys_state,ncol)
    integer,             intent(in)    :: ncol
    type(physics_state), intent(inout) :: phys_state
    logical :: match
    integer :: i, j, ulatcnt, uloncnt

    phys_state%ulat=-999._r8
    phys_state%ulon=-999._r8
    phys_state%latmapback=0
    phys_state%lonmapback=0
    match=.false.
    ulatcnt=0
    uloncnt=0
    match=.false.
    do i=1,ncol
       do j=1,ulatcnt
          if(phys_state%lat(i) .eq. phys_state%ulat(j)) then
             match=.true.
             phys_state%latmapback(i)=j
          end if
       end do
       if(.not. match) then
          ulatcnt=ulatcnt+1
          phys_state%ulat(ulatcnt)=phys_state%lat(i)
          phys_state%latmapback(i)=ulatcnt
       end if

       match=.false.
       do j=1,uloncnt
          if(phys_state%lon(i) .eq. phys_state%ulon(j)) then
             match=.true.
             phys_state%lonmapback(i)=j
          end if
       end do
       if(.not. match) then
          uloncnt=uloncnt+1
          phys_state%ulon(uloncnt)=phys_state%lon(i)
          phys_state%lonmapback(i)=uloncnt
       end if
       match=.false.

    end do
    phys_state%uloncnt=uloncnt
    phys_state%ulatcnt=ulatcnt

    call get_gcol_all_p(phys_state%lchnk,pcols,phys_state%cid)


  end subroutine init_geo_unique

!===============================================================================
  subroutine physics_dme_adjust(state, tend, qini, dt)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Adjust the dry mass in each layer back to the value of physics input state
    ! 
    ! Method: Conserve the integrated mass, momentum and total energy in each layer
    !         by scaling the specific mass of consituents, specific momentum (velocity)
    !         and specific total energy by the relative change in layer mass. Solve for
    !         the new temperature by subtracting the new kinetic energy from total energy
    !         and inverting the hydrostatic equation
    !
    !         The mass in each layer is modified, changing the relationship of the layer 
    !         interfaces and midpoints to the surface pressure. The result is no longer in 
    !         the original hybrid coordinate. 
    !
    !         This procedure cannot be applied to the "eul" or "sld" dycores because they
    !         require the hybrid coordinate.
    ! 
    ! Author: Byron Boville

    ! !REVISION HISTORY:
    !   03.03.28  Boville    Created, partly from code by Lin in p_d_adjust
    ! 
    !-----------------------------------------------------------------------

    use constituents, only : cnst_get_type_byind

    implicit none
    !
    ! Arguments
    !
    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    real(r8),            intent(in   ) :: qini(pcols,pver)    ! initial specific humidity
    real(r8),            intent(in   ) :: dt                  ! model physics timestep
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: lchnk         ! chunk identifier
    integer  :: ncol          ! number of atmospheric columns
    integer  :: i,k,m         ! Longitude, level indices
    real(r8) :: fdq(pcols)    ! mass adjustment factor
    real(r8) :: te(pcols)     ! total energy in a layer
    real(r8) :: utmp(pcols)   ! temp variable for recalculating the initial u values
    real(r8) :: vtmp(pcols)   ! temp variable for recalculating the initial v values

    real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer
    !
    !-----------------------------------------------------------------------
    ! verify that the dycore is FV
    if (.not. dycore_is('LR') ) return

    lchnk = state%lchnk
    ncol  = state%ncol

    ! adjust dry mass in each layer back to input value, while conserving
    ! constituents, momentum, and total energy
    do k = 1, pver

       ! adjusment factor is just change in water vapor
       fdq(:ncol) = 1._r8 + state%q(:ncol,k,1) - qini(:ncol,k)

       ! adjust constituents to conserve mass in each layer
       do m = 1, pcnst
          state%q(:ncol,k,m) = state%q(:ncol,k,m) / fdq(:ncol)
       end do

       if (adjust_te) then
          ! compute specific total energy of unadjusted state (J/kg)
          te(:ncol) = state%s(:ncol,k) + 0.5_r8*(state%u(:ncol,k)**2 + state%v(:ncol,k)**2) 

          ! recompute initial u,v from the new values and the tendencies
          utmp(:ncol) = state%u(:ncol,k) - dt * tend%dudt(:ncol,k)
          vtmp(:ncol) = state%v(:ncol,k) - dt * tend%dvdt(:ncol,k)
          ! adjust specific total energy and specific momentum (velocity) to conserve each
          te     (:ncol)   = te     (:ncol)     / fdq(:ncol)
          state%u(:ncol,k) = state%u(:ncol,k  ) / fdq(:ncol)
          state%v(:ncol,k) = state%v(:ncol,k  ) / fdq(:ncol)
          ! compute adjusted u,v tendencies
          tend%dudt(:ncol,k) = (state%u(:ncol,k) - utmp(:ncol)) / dt
          tend%dvdt(:ncol,k) = (state%v(:ncol,k) - vtmp(:ncol)) / dt

          ! compute adjusted static energy
          state%s(:ncol,k) = te(:ncol) - 0.5_r8*(state%u(:ncol,k)**2 + state%v(:ncol,k)**2)
       end if

! compute new total pressure variables
       state%pdel  (:ncol,k  ) = state%pdel(:ncol,k  ) * fdq(:ncol)
       state%pint  (:ncol,k+1) = state%pint(:ncol,k  ) + state%pdel(:ncol,k)
       state%lnpint(:ncol,k+1) = log(state%pint(:ncol,k+1))
       state%rpdel (:ncol,k  ) = 1._r8/ state%pdel(:ncol,k  )
    end do

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
      zvirv(:,:) = shr_const_rwv / rairv(:,:,state%lchnk) - 1._r8
    else
      zvirv(:,:) = zvir    
    endif

! compute new T,z from new s,q,dp
    if (adjust_te) then
       call geopotential_dse(state%lnpint, state%lnpmid  , state%pint ,  &
            state%pmid  , state%pdel    , state%rpdel,  &
            state%s     , state%q(1,1,1), state%phis , rairv(:,:,state%lchnk), &
	    gravit, cpairv(:,:,state%lchnk), zvirv, &
            state%t     , state%zi      , state%zm   , ncol)
    end if

  end subroutine physics_dme_adjust
!-----------------------------------------------------------------------

!===============================================================================
  subroutine physics_state_copy(state_in, state_out)
    
    use ppgrid,       only: pver, pverp
    use constituents, only: pcnst

    implicit none

    !
    ! Arguments
    !
    type(physics_state), intent(in) :: state_in
    type(physics_state), intent(out) :: state_out

    !
    ! Local variables
    !
    integer i, k, m, ncol


    ncol = state_in%ncol
    
    state_out%lchnk = state_in%lchnk
    state_out%ncol  = state_in%ncol
    state_out%count = state_in%count 

    do i = 1, ncol
       state_out%lat(i)    = state_in%lat(i)
       state_out%lon(i)    = state_in%lon(i)
       state_out%ps(i)     = state_in%ps(i)
       state_out%phis(i)   = state_in%phis(i)
       state_out%te_ini(i) = state_in%te_ini(i) 
       state_out%te_cur(i) = state_in%te_cur(i) 
       state_out%tw_ini(i) = state_in%tw_ini(i) 
       state_out%tw_cur(i) = state_in%tw_cur(i) 
    end do

    do k = 1, pver
       do i = 1, ncol
          state_out%t(i,k)         = state_in%t(i,k) 
          state_out%u(i,k)         = state_in%u(i,k) 
          state_out%v(i,k)         = state_in%v(i,k) 
          state_out%s(i,k)         = state_in%s(i,k) 
          state_out%omega(i,k)     = state_in%omega(i,k) 
          state_out%pmid(i,k)      = state_in%pmid(i,k) 
          state_out%pdel(i,k)      = state_in%pdel(i,k) 
          state_out%rpdel(i,k)     = state_in%rpdel(i,k) 
          state_out%lnpmid(i,k)    = state_in%lnpmid(i,k) 
          state_out%exner(i,k)     = state_in%exner(i,k) 
          state_out%zm(i,k)        = state_in%zm(i,k)
       end do
    end do

    do k = 1, pverp
       do i = 1, ncol
          state_out%pint(i,k)      = state_in%pint(i,k) 
          state_out%lnpint(i,k)    = state_in%lnpint(i,k) 
          state_out%zi(i,k)        = state_in% zi(i,k) 
       end do
    end do


       do i = 1, ncol
          state_out%psdry(i)  = state_in%psdry(i) 
       end do
       do k = 1, pver
          do i = 1, ncol
             state_out%lnpmiddry(i,k) = state_in%lnpmiddry(i,k) 
             state_out%pmiddry(i,k)   = state_in%pmiddry(i,k) 
             state_out%pdeldry(i,k)   = state_in%pdeldry(i,k) 
             state_out%rpdeldry(i,k)  = state_in%rpdeldry(i,k) 
          end do
       end do
       do k = 1, pverp
          do i = 1, ncol
             state_out%pintdry(i,k)   = state_in%pintdry(i,k)
             state_out%lnpintdry(i,k) = state_in%lnpintdry(i,k) 
          end do
       end do

    do m = 1, pcnst
       do k = 1, pver
          do i = 1, ncol
             state_out%q(i,k,m) = state_in%q(i,k,m) 
          end do
       end do
    end do

  end  subroutine physics_state_copy
!===============================================================================

  subroutine physics_tend_init(tend)
    
    implicit none
    
    !
    ! Arguments
    !
    type(physics_tend), intent(inout) :: tend

    !
    ! Local variables
    !

    tend%dtdt    = 0._r8
    tend%dudt    = 0._r8
    tend%dvdt    = 0._r8
    tend%flx_net = 0._r8
    tend%te_tnd  = 0._r8
    tend%tw_tnd  = 0._r8
    
end subroutine physics_tend_init

!===============================================================================

subroutine set_state_pdry (state,pdeld_calc)

  use ppgrid,  only: pver
  use pmgrid,  only: plev, plevp
  implicit none

  type(physics_state), intent(inout) :: state
  logical, optional, intent(in) :: pdeld_calc    !  .true. do calculate pdeld [default]
                                                 !  .false. don't calculate pdeld 
  integer ncol
  integer i, k
  logical do_pdeld_calc

  if ( present(pdeld_calc) ) then
     do_pdeld_calc = pdeld_calc
  else
     do_pdeld_calc = .true.
  endif
  
  ncol = state%ncol


  state%psdry(:ncol) = state%pint(:ncol,1)
  state%pintdry(:ncol,1) = state%pint(:ncol,1)

  if (do_pdeld_calc)  then
     do k = 1, pver
        state%pdeldry(:ncol,k) = state%pdel(:ncol,k)*(1._r8-state%q(:ncol,k,1))
     end do
  endif
  do k = 1, pver
     state%pintdry(:ncol,k+1) = state%pintdry(:ncol,k)+state%pdeldry(:ncol,k)
     state%pmiddry(:ncol,k) = (state%pintdry(:ncol,k+1)+state%pintdry(:ncol,k))/2._r8
     state%psdry(:ncol) = state%psdry(:ncol) + state%pdeldry(:ncol,k)
  end do

  state%rpdeldry(:ncol,:) = 1._r8/state%pdeldry(:ncol,:)
  state%lnpmiddry(:ncol,:) = log(state%pmiddry(:ncol,:))
  state%lnpintdry(:ncol,:) = log(state%pintdry(:ncol,:))

end subroutine set_state_pdry 

!===============================================================================

subroutine set_wet_to_dry (state)

  use constituents,  only: pcnst, cnst_type

  type(physics_state), intent(inout) :: state

  integer m, ncol
  
  ncol = state%ncol

  do m = 1,pcnst
     if (cnst_type(m).eq.'dry') then
        state%q(:ncol,:,m) = state%q(:ncol,:,m)*state%pdel(:ncol,:)/state%pdeldry(:ncol,:)
     endif
  end do

end subroutine set_wet_to_dry 

!===============================================================================

subroutine set_dry_to_wet (state)

  use constituents,  only: pcnst, cnst_type

  type(physics_state), intent(inout) :: state

  integer m, ncol
  
  ncol = state%ncol

  do m = 1,pcnst
     if (cnst_type(m).eq.'dry') then
        state%q(:ncol,:,m) = state%q(:ncol,:,m)*state%pdeldry(:ncol,:)/state%pdel(:ncol,:)
     endif
  end do

end subroutine set_dry_to_wet


end module physics_types
