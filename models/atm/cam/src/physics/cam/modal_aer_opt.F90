module modal_aer_opt

! parameterizes aerosol coefficients using chebychev polynomial
! parameterize aerosol radiative properties in terms of
! surface mode wet radius and wet refractive index

! Ghan and Zaveri, JGR 2007.

! uses Wiscombe's (1979) mie scattering code


use shr_kind_mod,      only: r8 => shr_kind_r8, shr_kind_cl
use ppgrid,            only: pcols, pver, pverp
use constituents,      only: pcnst
use spmd_utils,        only: masterproc
use phys_control,      only: cam_chempkg_is
use ref_pres,          only: top_lev=>trop_cloud_top_lev
use physconst,         only: rhoh2o, rga, rair
use radconstants,      only: nswbands, nlwbands, idx_sw_diag
use rad_constituents,  only: n_diag, rad_cnst_get_call_list, rad_cnst_get_info, rad_cnst_get_aer_mmr, &
                             rad_cnst_get_aer_props, rad_cnst_get_mode_props
use physics_types,     only: physics_state

use physics_buffer, only : pbuf_get_index,physics_buffer_desc, pbuf_get_field
use pio,               only: file_desc_t, var_desc_t, pio_inq_dimlen, pio_inq_dimid, pio_inq_varid, &
                             pio_get_var, pio_nowrite, pio_closefile
use cam_pio_utils,     only: cam_pio_openfile
use cam_history,       only: phys_decomp, addfld, add_default, outfld
use cam_history_support, only: fillvalue
use cam_logfile,       only: iulog
use perf_mod,          only: t_startf, t_stopf
use abortutils,        only: endrun

implicit none
private
save

public :: modal_aer_opt_readnl, modal_aer_opt_init, modal_aero_sw, modal_aero_lw


character(len=*), parameter :: unset_str = 'UNSET'

! Namelist variables:
character(shr_kind_cl)      :: modal_optics_file = unset_str   ! full pathname for modal optics dataset
character(shr_kind_cl)      :: water_refindex_file = unset_str ! full pathname for water refractive index dataset

! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
! in terms of refractive index and wet radius
integer, parameter :: ncoef=5, prefr=7, prefi=10

real(r8) :: xrmin, xrmax

! refractive index for water read in read_water_refindex
complex(r8) :: crefwsw(nswbands) ! complex refractive index for water visible
complex(r8) :: crefwlw(nlwbands) ! complex refractive index for water infrared

! physics buffer indices
integer :: dgnumwet_idx = -1
integer :: qaerwat_idx  = -1

!===============================================================================
CONTAINS
!===============================================================================

subroutine modal_aer_opt_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'modal_aer_opt_readnl'

   namelist /modal_aer_opt_nl/ water_refindex_file
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'modal_aer_opt_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, modal_aer_opt_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   call mpibcast(water_refindex_file, len(water_refindex_file), mpichar, 0, mpicom)
#endif


end subroutine modal_aer_opt_readnl

!===============================================================================

subroutine modal_aer_opt_init()

   use ioFileMod,        only: getfil

   ! Local variables

   integer  :: i, m
   real(r8) :: rmmin, rmmax       ! min, max aerosol surface mode radius treated (m)
   character(len=256) :: locfile

   logical :: call_list(0:n_diag)
   integer :: ilist, nmodes, m_ncoef, m_prefr, m_prefi
   integer :: errcode

   character(len=*), parameter :: routine='modal_aer_opt_init'
   !----------------------------------------------------------------------------

   rmmin = 0.01e-6_r8
   rmmax = 25.e-6_r8
   xrmin = log(rmmin)
   xrmax = log(rmmax)

   ! Check that dimension sizes in the coefficient arrays used to
   ! parameterize aerosol radiative properties are consistent between this
   ! module and the mode physprop files.
   call rad_cnst_get_call_list(call_list)
   do ilist = 0, n_diag
      if (call_list(ilist)) then
         call rad_cnst_get_info(ilist, nmodes=nmodes)
         do m = 1, nmodes
            call rad_cnst_get_mode_props(ilist, m, ncoef=m_ncoef, prefr=m_prefr, prefi=m_prefi)
            if (m_ncoef /= ncoef .or. m_prefr /= prefr .or. m_prefi /= prefi) then
               write(iulog,*) routine//': ERROR - file and module values do not match:'
               write(iulog,*) '   ncoef:', ncoef, m_ncoef
               write(iulog,*) '   prefr:', prefr, m_prefr
               write(iulog,*) '   prefi:', prefi, m_prefi
               call endrun(routine//': ERROR - file and module values do not match')
            end if
         end do
      end if
   end do

   ! Initialize physics buffer indices for dgnumwet and qaerwat.  Note the implicit assumption
   ! that the loops over modes in the optics calculations will use the values for dgnumwet and qaerwat
   ! that are set in the aerosol_wet_intr code.
   dgnumwet_idx = pbuf_get_index('DGNUMWET',errcode)
   if (errcode < 0) then
      call endrun(routine//' ERROR: cannot find physics buffer field DGNUMWET')
   end if
   qaerwat_idx  = pbuf_get_index('QAERWAT')
   if (errcode < 0) then
      call endrun(routine//' ERROR: cannot find physics buffer field QAERWAT')
   end if

   call getfil(water_refindex_file, locfile)
   call read_water_refindex(locfile)
   if (masterproc) write(iulog,*) "modal_aer_opt_init: read water refractive index file:", trim(locfile)


   ! Add diagnostic fields to history output.
   call addfld ('EXTINCT','/m  ',pver,    'A','Aerosol extinction',phys_decomp, flag_xyfill=.true.)
   call addfld ('ABSORB','/m  ',pver,    'A','Aerosol absorption',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODVIS','  ',1,    'A','Aerosol optical depth 550 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODABS','  ',1,    'A','Aerosol absorption optical depth 550 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODMODE1','  ',1,    'A','Aerosol optical depth 550 nm mode 1',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODMODE2','  ',1,    'A','Aerosol optical depth 550 nm mode 2',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODMODE3','  ',1,    'A','Aerosol optical depth 550 nm mode 3',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODDUST1','  ',1,    'A','Aerosol optical depth 550 nm model 1 from dust',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODDUST2','  ',1,    'A','Aerosol optical depth 550 nm model 2 from dust',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODDUST3','  ',1,    'A','Aerosol optical depth 550 nm model 3 from dust',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDEN1','kg/m2',1,    'A','Aerosol burden mode 1',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDEN2','kg/m2',1,    'A','Aerosol burden mode 2',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDEN3','kg/m2',1,    'A','Aerosol burden mode 3',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENDUST','kg/m2',1,    'A','Dust aerosol burden',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENSO4','kg/m2',1,    'A','Sulfate aerosol burden',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENPOM','kg/m2',1,    'A','POM aerosol burden',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENSOA','kg/m2',1,    'A','SOA aerosol burden',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENBC','kg/m2',1,    'A','Black carbon aerosol burden',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENSEASALT','kg/m2',1,    'A','Seasalt aerosol burden',phys_decomp, flag_xyfill=.true.)
   call addfld ('SSAVIS','  ',1,    'A','Aerosol singel-scatter albedo',phys_decomp, flag_xyfill=.true.)

   call add_default ('EXTINCT', 1, ' ')
   call add_default ('ABSORB', 1, ' ')
   call add_default ('AODVIS', 1, ' ')
   call add_default ('AODABS', 1, ' ')
   call add_default ('AODMODE1', 1, ' ')
   call add_default ('AODMODE2', 1, ' ')
   call add_default ('AODMODE3', 1, ' ')
   call add_default ('AODDUST1', 1, ' ')
   call add_default ('AODDUST2', 1, ' ')
   call add_default ('AODDUST3', 1, ' ')
   call add_default ('BURDEN1', 1, ' ')
   call add_default ('BURDEN2', 1, ' ')
   call add_default ('BURDEN3', 1, ' ')
   call add_default ('BURDENDUST', 1, ' ')
   call add_default ('BURDENSO4', 1, ' ')
   call add_default ('BURDENPOM', 1, ' ')
   call add_default ('BURDENSOA', 1, ' ')
   call add_default ('BURDENBC', 1, ' ')
   call add_default ('BURDENSEASALT', 1, ' ')
   call add_default ('SSAVIS', 1, ' ')

   if (cam_chempkg_is('trop_mam7').or.cam_chempkg_is('trop_strat_mam7')) then
      call addfld ('AODDUST4','  ',1,    'A','Aerosol optical depth 550 nm model 4 from dust',phys_decomp, flag_xyfill=.true.)
      call addfld ('AODDUST5','  ',1,    'A','Aerosol optical depth 550 nm model 5 from dust',phys_decomp, flag_xyfill=.true.)
      call addfld ('AODDUST6','  ',1,    'A','Aerosol optical depth 550 nm model 6 from dust',phys_decomp, flag_xyfill=.true.)
      call addfld ('AODDUST7','  ',1,    'A','Aerosol optical depth 550 nm model 7 from dust',phys_decomp, flag_xyfill=.true.)
      call add_default ('AODDUST4', 1, ' ')
      call add_default ('AODDUST5', 1, ' ')
      call add_default ('AODDUST6', 1, ' ')
      call add_default ('AODDUST7', 1, ' ')

      call addfld ('AODMODE4','  ',1,    'A','Aerosol optical depth 550 nm mode 4',phys_decomp, flag_xyfill=.true.)
      call addfld ('AODMODE5','  ',1,    'A','Aerosol optical depth 550 nm mode 5',phys_decomp, flag_xyfill=.true.)
      call addfld ('AODMODE6','  ',1,    'A','Aerosol optical depth 550 nm mode 6',phys_decomp, flag_xyfill=.true.)
      call addfld ('AODMODE7','  ',1,    'A','Aerosol optical depth 550 nm mode 7',phys_decomp, flag_xyfill=.true.)
      call addfld ('BURDEN4','kg/m2',1,    'A','Aerosol burden mode 4',phys_decomp, flag_xyfill=.true.)
      call addfld ('BURDEN5','kg/m2',1,    'A','Aerosol burden mode 5',phys_decomp, flag_xyfill=.true.)
      call addfld ('BURDEN6','kg/m2',1,    'A','Aerosol burden mode 6',phys_decomp, flag_xyfill=.true.)
      call addfld ('BURDEN7','kg/m2',1,    'A','Aerosol burden mode 7',phys_decomp, flag_xyfill=.true.)
      call add_default ('AODMODE4', 1, ' ')
      call add_default ('AODMODE5', 1, ' ')
      call add_default ('AODMODE6', 1, ' ')
      call add_default ('AODMODE7', 1, ' ')
      call add_default ('BURDEN4', 1, ' ')
      call add_default ('BURDEN5', 1, ' ')
      call add_default ('BURDEN6', 1, ' ')
      call add_default ('BURDEN7', 1, ' ')
   end if

end subroutine modal_aer_opt_init

!===============================================================================

subroutine modal_aero_sw(list_idx, state, pbuf, nnite, idxnite, &
                         tauxar, wa, ga, fa)

   ! calculates aerosol sw radiative properties

   integer,             intent(in) :: list_idx       ! index of the climate or a diagnostic list
   type(physics_state), intent(in) :: state          ! state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer,             intent(in) :: nnite          ! number of night columns
   integer,             intent(in) :: idxnite(nnite) ! local column indices of night columns

   real(r8), intent(out) :: tauxar(pcols,0:pver,nswbands) ! layer extinction optical depth
   real(r8), intent(out) :: wa(pcols,0:pver,nswbands)     ! layer single-scatter albedo
   real(r8), intent(out) :: ga(pcols,0:pver,nswbands)     ! asymmetry factor
   real(r8), intent(out) :: fa(pcols,0:pver,nswbands)     ! forward scattered fraction

   ! Local variables
   integer :: i, ifld, isw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec

   real(r8) :: mass(pcols,pver)        ! layer mass
   real(r8) :: air_density(pcols,pver) ! (kg/m3)

   real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specdens            ! species density (kg/m3)
   complex(r8), pointer :: specrefindex(:)     ! species refractive index
   character*32         :: spectype            ! species type

   real(r8), pointer :: dgnumwet(:,:)     ! number mode diameter
   real(r8), pointer :: qaerwat(:,:)      ! aerosol water (g/g)

   real(r8) :: sigma_logr_aer         ! geometric standard deviation of number distribution
   real(r8) :: radsurf(pcols,pver)    ! aerosol surface mode radius
   real(r8) :: logradsurf(pcols,pver) ! log(aerosol surface mode radius)
   real(r8) :: cheb(ncoef,pcols,pver)

   real(r8)    :: refr(pcols)     ! real part of refractive index
   real(r8)    :: refi(pcols)     ! imaginary part of refractive index
   complex(r8) :: crefin(pcols)   ! complex refractive index
   real(r8), pointer :: refrtabsw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitabsw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: extpsw(:,:,:,:) ! specific extinction
   real(r8), pointer :: abspsw(:,:,:,:) ! specific absorption
   real(r8), pointer :: asmpsw(:,:,:,:) ! asymmetry factor

   real(r8) :: vol(pcols)      ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)   ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: watervol(pcols) ! volume concentration of water in each mode (m3/kg)
   real(r8) :: wetvol(pcols)   ! volume concentration of wet mode (m3/kg)

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cext(pcols,ncoef), cabs(pcols,ncoef), casm(pcols,ncoef)
   real(r8) :: pext(pcols)     ! parameterized specific extinction (m2/kg)
   real(r8) :: specpext(pcols) ! specific extinction (m2/kg)
   real(r8) :: dopaer(pcols)   ! aerosol optical depth in layer
   real(r8) :: pabs(pcols)     ! parameterized specific absorption (m2/kg)
   real(r8) :: pasm(pcols)     ! parameterized asymmetry factor
   real(r8) :: palb(pcols)     ! parameterized single scattering albedo

   ! Diagnostics output for visible band only
   real(r8) :: extinct(pcols,pver)
   real(r8) :: absorb(pcols,pver)
   real(r8) :: aodvis(pcols)               ! extinction optical depth
   real(r8) :: aodabs(pcols)               ! absorption optical depth
   real(r8) :: ssavis(pcols)
   real(r8) :: dustvol(pcols)              ! volume concentration of dust in aerosol mode (m3/kg)
   real(r8) :: burden(pcols)
   real(r8) :: burdendust(pcols)
   real(r8) :: burdenso4(pcols)
   real(r8) :: burdenpom(pcols)
   real(r8) :: burdensoa(pcols)
   real(r8) :: burdenbc(pcols)
   real(r8) :: burdenseasalt(pcols)
   real(r8) :: aodmode(pcols)
   real(r8) :: dustaodmode(pcols)          ! dust aod in aerosol mode

   logical :: savaervis ! true if visible wavelength (0.55 micron)
   character(len=32) :: outname

   ! debug output
   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf            ! volume fraction of insoluble aerosol
   character(len=*), parameter :: subname = 'modal_aero_sw'
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8
   wa(:ncol,:,:)     = 0._r8
   ga(:ncol,:,:)     = 0._r8
   fa(:ncol,:,:)     = 0._r8

   ! zero'th layer does not contain aerosol
   tauxar(1:ncol,0,:)  = 0._r8
   wa(1:ncol,0,:)      = 0.925_r8
   ga(1:ncol,0,:)      = 0.850_r8
   fa(1:ncol,0,:)      = 0.7225_r8

   mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
   air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

   ! diagnostics for visible band summed over modes
   extinct(1:ncol,:)     = 0.0_r8
   absorb(1:ncol,:)      = 0.0_r8
   aodvis(1:ncol)        = 0.0_r8
   aodabs(1:ncol)        = 0.0_r8
   ssavis(1:ncol)        = 0.0_r8
   burdendust(:ncol)     = 0._r8
   burdenso4(:ncol)      = 0._r8
   burdenpom(:ncol)      = 0._r8
   burdensoa(:ncol)      = 0._r8
   burdenbc(:ncol)       = 0._r8
   burdenseasalt(:ncol)  = 0._r8

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   do m = 1, nmodes

      ! diagnostics for visible band for each mode
      burden(:ncol)       = 0._r8
      aodmode(1:ncol)     = 0.0_r8
      dustaodmode(1:ncol) = 0.0_r8

      ! pointers to physics buffer

      ! ***N.B.*** It is implicitly assumed that the loops over modes in the
      ! optics calculations will use the values for dgnumwet and qaerwat
      ! that are set in the aerosol_wet_intr code.

      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet, start=(/1,1,m/), kount=(/pcols,pver,1/) )
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat,  start=(/1,1,m/), kount=(/pcols,pver,1/) )

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtabsw=refrtabsw , &
         refitabsw=refitabsw, extpsw=extpsw, abspsw=abspsw, asmpsw=asmpsw)

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      ! calc size parameter for all columns
      call modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

      do isw = 1, nswbands
         savaervis = (isw .eq. idx_sw_diag)

         do k = top_lev, pver

            ! form bulk refractive index
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8
            dustvol(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec
               call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
               call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                           refindex_aer_sw=specrefindex, spectype=spectype)

               do i = 1, ncol
                  vol(i)      = specmmr(i,k)/specdens
                  dryvol(i)   = dryvol(i) + vol(i)
                  crefin(i)   = crefin(i) + vol(i)*specrefindex(isw)
               end do

               ! compute some diagnostics for visible band only
               if (savaervis) then

                  do i = 1, ncol
                     burden(i) = burden(i) + specmmr(i,k)*mass(i,k)
                  end do

                  if (trim(spectype) == 'dust') then
                     do i = 1, ncol
                        dustvol(i) = specmmr(i,k)/specdens
                        burdendust(i)=burdendust(i) + specmmr(i,k)*mass(i,k)
                     end do
                  end if
                  if (trim(spectype) == 'sulfate') then
                     do i = 1, ncol
                        burdenso4(i)=burdenso4(i) + specmmr(i,k)*mass(i,k)
                     end do
                  end if
                  if (trim(spectype) == 'black-c') then
                     do i = 1, ncol
                        burdenbc(i)=burdenbc(i) + specmmr(i,k)*mass(i,k)
                   end do
                  end if
                  if (trim(spectype) == 'p-organic') then
                     do i = 1, ncol
                        burdenpom(i)=burdenpom(i) + specmmr(i,k)*mass(i,k)
                      end do
                  end if
                  if (trim(spectype) == 's-organic') then
                     do i = 1, ncol
                        burdensoa(i)=burdensoa(i) + specmmr(i,k)*mass(i,k)
                     end do
                  end if
                  if (trim(spectype) == 'seasalt') then
                     do i = 1, ncol
                        burdenseasalt(i)=burdenseasalt(i) + specmmr(i,k)*mass(i,k)
                      end do
                  end if

               end if
            end do ! species loop

            do i = 1, ncol
               watervol(i) = qaerwat(i,k)/rhoh2o
               wetvol(i) = watervol(i) + dryvol(i)
               if (watervol(i) < 0._r8) then
                  if (abs(watervol(i)) .gt. 1.e-1_r8*wetvol(i)) then
                     write(iulog,'(a,2e10.2,a)') 'watervol,wetvol=', &
                        watervol(i), wetvol(i), ' in '//subname
                  end if
                  watervol(i) = 0._r8
                  wetvol(i) = dryvol(i)
               end if

               ! volume mixing
               crefin(i) = crefin(i) + watervol(i)*crefwsw(isw)
               crefin(i) = crefin(i)/max(wetvol(i),1.e-60_r8)
               refr(i)   = real(crefin(i))
               refi(i)   = abs(aimag(crefin(i)))
            end do

            ! call t_startf('binterp')

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(extpsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                         itab, jtab, ttab, utab, cext)
            call binterp(abspsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                         itab, jtab, ttab, utab, cabs)
            call binterp(asmpsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                         itab, jtab, ttab, utab, casm)

            ! call t_stopf('binterp')

            ! parameterized optical properties
            do i=1,ncol

               if (logradsurf(i,k) .le. xrmax) then
                  pext(i) = 0.5_r8*cext(i,1)
                  do nc = 2, ncoef
                     pext(i) = pext(i) + cheb(nc,i,k)*cext(i,nc)
                  enddo
                  pext(i) = exp(pext(i))
               else
                  pext(i) = 1.5_r8/(radsurf(i,k)*rhoh2o) ! geometric optics
               endif

               ! convert from m2/kg water to m2/kg aerosol
               specpext(i) = pext(i)
               pext(i) = pext(i)*wetvol(i)*rhoh2o
               pabs(i) = 0.5_r8*cabs(i,1)
               pasm(i) = 0.5_r8*casm(i,1)
               do nc = 2, ncoef
                  pabs(i) = pabs(i) + cheb(nc,i,k)*cabs(i,nc)
                  pasm(i) = pasm(i) + cheb(nc,i,k)*casm(i,nc)
               enddo
               pabs(i) = pabs(i)*wetvol(i)*rhoh2o
               pabs(i) = max(0._r8,pabs(i))
               pabs(i) = min(pext(i),pabs(i))

               palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)
               palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)

               dopaer(i) = pext(i)*mass(i,k)
            end do

            ! Save aerosol optical depth at longest visible wavelength
            ! sum over layers
            if (savaervis) then
               ! aerosol extinction (/m)
               do i = 1, ncol
                  extinct(i,k) = extinct(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                  absorb(i,k)  = absorb(i,k) + pabs(i)*air_density(i,k)
                  aodvis(i)    = aodvis(i) + dopaer(i)
                  aodabs(i)    = aodabs(i) + pabs(i)*mass(i,k)
                  aodmode(i)   = aodmode(i) + dopaer(i)
                  if (wetvol(i) > 1.e-40_r8) then
                     dustaodmode(i) = dustaodmode(i) + dopaer(i)*dustvol(i)/wetvol(i)
                  endif
                  ssavis(i) = ssavis(i) + dopaer(i)*palb(i)
               end do
            endif

            do i = 1, ncol

               if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 30._r8)) then

                  write(iulog,*) 'dopaer(', i, ',', k, ',', m, ',', lchnk, ')=', dopaer(i)
                  ! write(iulog,*) 'itab,jtab,ttab,utab=',itab(i),jtab(i),ttab(i),utab(i)
                  write(iulog,*) 'k=', k, ' pext=', pext(i), ' specext=', specpext(i)
                  write(iulog,*) 'wetvol=', wetvol(i), ' dryvol=', dryvol(i), ' watervol=', watervol(i)
                  ! write(iulog,*) 'cext=',(cext(i,l),l=1,ncoef)
                  ! write(iulog,*) 'crefin=',crefin(i)
                  write(iulog,*) 'nspec=', nspec
                  ! write(iulog,*) 'cheb=', (cheb(nc,m,i,k),nc=2,ncoef)
                  do l = 1, nspec
                     call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
                     call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                                 refindex_aer_sw=specrefindex)
                     volf = specmmr(i,k)/specdens
                     write(iulog,*) 'l=', l, 'vol(l)=', volf
                     write(iulog,*) 'isw=', isw, 'specrefindex(isw)=', specrefindex(isw)
                     write(iulog,*) 'specdens=', specdens
                  end do

                  nerr_dopaer = nerr_dopaer + 1
                  if (nerr_dopaer >= nerrmax_dopaer) then
                     ! write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                     ! call endrun('exit from '//subname)
                  end if

               end if
            end do

            do i=1,ncol
               tauxar(i,k,isw) = tauxar(i,k,isw) + dopaer(i)
               wa(i,k,isw)     = wa(i,k,isw)     + dopaer(i)*palb(i)
               ga(i,k,isw)     = ga(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)
               fa(i,k,isw)     = fa(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)*pasm(i)
            end do

         end do ! pver

      end do ! sw bands

      ! mode diagnostics
      ! The diagnostics are currently only output for the climate list.  Code mods will
      ! be necessary to provide output for the rad_diag lists.
      if (list_idx == 0) then
         do i = 1, nnite
            burden(idxnite(i))  = fillvalue
            aodmode(idxnite(i)) = fillvalue
            dustaodmode(idxnite(i)) = fillvalue
         end do

         write(outname,'(a,i1)') 'BURDEN', m
         call outfld(trim(outname), burden, pcols, lchnk)

         write(outname,'(a,i1)') 'AODMODE', m
         call outfld(trim(outname), aodmode, pcols, lchnk)

         write(outname,'(a,i1)') 'AODDUST', m
         call outfld(trim(outname), dustaodmode, pcols, lchnk)

      end if

   end do ! nmodes

   ! Output visible band diagnostics for quantities summed over the modes
   ! The diagnostics are currently only output for the climate list.  Code mods will
   ! be necessary to provide output for the rad_diag lists.
   if (list_idx == 0) then
      do i = 1, ncol
         if (aodvis(i) > 1.e-10_r8) then
            ssavis(i) = ssavis(i)/aodvis(i)
         else
            ssavis(i) = 0.925_r8
         endif
      end do

      do i = 1, nnite
         extinct(idxnite(i),:) = fillvalue
         absorb(idxnite(i),:)  = fillvalue
         aodvis(idxnite(i))    = fillvalue
         aodabs(idxnite(i))    = fillvalue
         ssavis(idxnite(i))    = fillvalue
      end do

      call outfld('EXTINCT',  extinct, pcols, lchnk)
      call outfld('ABSORB',   absorb,  pcols, lchnk)
      call outfld('AODVIS',   aodvis,  pcols, lchnk)
      call outfld('AODABS',   aodabs,  pcols, lchnk)
      call outfld('SSAVIS',   ssavis,  pcols, lchnk)
      call outfld('BURDENDUST',burdendust ,pcols,lchnk)
      call outfld('BURDENSO4' ,burdenso4  ,pcols,lchnk)
      call outfld('BURDENPOM' ,burdenpom  ,pcols,lchnk)
      call outfld('BURDENSOA' ,burdensoa  ,pcols,lchnk)
      call outfld('BURDENBC'  ,burdenbc   ,pcols,lchnk)
      call outfld('BURDENSEASALT',burdenseasalt ,pcols,lchnk)
   end if

end subroutine modal_aero_sw

!===============================================================================

subroutine modal_aero_lw(list_idx, state, pbuf, tauxar)

   ! calculates aerosol lw radiative properties

   integer,             intent(in)  :: list_idx ! index of the climate or a diagnostic list
   type(physics_state), intent(in)  :: state    ! state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(out) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth

   ! Local variables
   integer :: i, ifld, ilw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec

   real(r8), pointer :: dgnumwet(:,:)  ! wet number mode diameter (m)
   real(r8), pointer :: qaerwat(:,:)   ! aerosol water (g/g)

   real(r8) :: sigma_logr_aer          ! geometric standard deviation of number distribution
   real(r8) :: alnsg_amode
   real(r8) :: xrad(pcols)
   real(r8) :: cheby(ncoef,pcols,pver)  ! chebychef polynomials

   real(r8) :: mass(pcols,pver) ! layer mass

   real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specdens            ! species density (kg/m3)
   complex(r8), pointer :: specrefindex(:)     ! species refractive index

   real(r8) :: vol(pcols)       ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)    ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: wetvol(pcols)    ! volume concentration of wet mode (m3/kg)
   real(r8) :: watervol(pcols)  ! volume concentration of water in each mode (m3/kg)
   real(r8) :: refr(pcols)      ! real part of refractive index
   real(r8) :: refi(pcols)      ! imaginary part of refractive index
   complex(r8) :: crefin(pcols) ! complex refractive index
   real(r8), pointer :: refrtablw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitablw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: absplw(:,:,:,:) ! specific absorption

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cabs(pcols,ncoef)
   real(r8) :: pabs(pcols)      ! parameterized specific absorption (m2/kg)
   real(r8) :: dopaer(pcols)    ! aerosol optical depth in layer

   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf             ! volume fraction of insoluble aerosol

   character(len=*), parameter :: subname = 'modal_aero_lw'
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8

   ! dry mass in each cell
   mass(:ncol,:) = state%pdeldry(:ncol,:)*rga

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   do m = 1, nmodes

      ! pointers to physics buffer
      ! ***N.B.*** It is implicitly assumed that the loops over modes in the
      ! optics calculations will use the values for dgnumwet and qaerwat
      ! that are set in the aerosol_wet_intr code.

      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet, start=(/1,1,m/), kount=(/pcols,pver,1/) )
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat,  start=(/1,1,m/), kount=(/pcols,pver,1/) )

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtablw=refrtablw , &
         refitablw=refitablw, absplw=absplw)

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      ! calc size parameter for all columns
      ! this is the same calculation that's done in modal_size_parameters, but there
      ! some intermediate results are saved and the chebyshev polynomials are stored
      ! in a array with different index order.  Could be unified.
      do k = top_lev, pver
         do i = 1, ncol
            alnsg_amode = log( sigma_logr_aer )
            ! convert from number diameter to surface area
            xrad(i) = log(0.5_r8*dgnumwet(i,k)) + 2.0_r8*alnsg_amode*alnsg_amode
            ! normalize size parameter
            xrad(i) = max(xrad(i), xrmin)
            xrad(i) = min(xrad(i), xrmax)
            xrad(i) = (2*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
            ! chebyshev polynomials
            cheby(1,i,k) = 1.0_r8
            cheby(2,i,k) = xrad(i)
            do nc = 3, ncoef
               cheby(nc,i,k) = 2.0_r8*xrad(i)*cheby(nc-1,i,k)-cheby(nc-2,i,k)
            end do
         end do
      end do

      do ilw = 1, nlwbands

         do k = top_lev, pver

            ! form bulk refractive index. Use volume mixing for infrared
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec
               call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
               call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                           refindex_aer_lw=specrefindex)

               do i = 1, ncol
                  vol(i)    = specmmr(i,k)/specdens
                  dryvol(i) = dryvol(i) + vol(i)
                  crefin(i) = crefin(i) + vol(i)*specrefindex(ilw)
               end do
            end do

            do i = 1, ncol
               watervol(i) = qaerwat(i,k)/rhoh2o
               wetvol(i)   = watervol(i) + dryvol(i)
               if (watervol(i) < 0.0_r8) then
                  if (abs(watervol(i)) .gt. 1.e-1_r8*wetvol(i)) then
                     write(iulog,*) 'watervol,wetvol,dryvol=',watervol(i),wetvol(i),dryvol(i),' in '//subname
                  end if
                  watervol(i) = 0._r8
                  wetvol(i)   = dryvol(i)
               end if

               crefin(i) = crefin(i) + watervol(i)*crefwlw(ilw)
               if (wetvol(i) > 1.e-40_r8) crefin(i) = crefin(i)/wetvol(i)
               refr(i) = real(crefin(i))
               refi(i) = aimag(crefin(i))
            end do

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(absplw(:,:,:,ilw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtablw(:,ilw), refitablw(:,ilw), &
                         itab, jtab, ttab, utab, cabs)

            ! parameterized optical properties
            do i = 1, ncol
               pabs(i) = 0.5_r8*cabs(i,1)
               do nc = 2, ncoef
                  pabs(i) = pabs(i) + cheby(nc,i,k)*cabs(i,nc)
               end do
               pabs(i)   = pabs(i)*wetvol(i)*rhoh2o
               pabs(i)   = max(0._r8,pabs(i))
               dopaer(i) = pabs(i)*mass(i,k)
            end do

            do i = 1, ncol

               if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 20._r8)) then

                  write(iulog,*) 'dopaer(',i,',',k,',',m,',',lchnk,')=', dopaer(i)
                  write(iulog,*) 'k=',k,' pabs=', pabs(i)
                  write(iulog,*) 'wetvol=',wetvol(i),' dryvol=',dryvol(i),     &
                     ' watervol=',watervol(i)
                  write(iulog,*) 'cabs=', (cabs(i,l),l=1,ncoef)
                  write(iulog,*) 'crefin=', crefin(i)
                  write(iulog,*) 'nspec=', nspec
                  do l = 1,nspec
                     call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
                     call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                                 refindex_aer_lw=specrefindex)
                     volf = specmmr(i,k)/specdens
                     write(iulog,*) 'l=',l,'vol(l)=',volf
                     write(iulog,*) 'ilw=',ilw,' specrefindex(ilw)=',specrefindex(ilw)
                     write(iulog,*) 'specdens=',specdens
                  end do

                  nerr_dopaer = nerr_dopaer + 1
                  if (nerr_dopaer >= nerrmax_dopaer) then
                     write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                     call endrun()
                  end if

               end if
            end do

            do i = 1, ncol
               tauxar(i,k,ilw) = tauxar(i,k,ilw) + dopaer(i)
            end do

         end do ! k = top_lev, pver

      end do  ! nlwbands

   end do ! m = 1, nmodes

end subroutine modal_aero_lw

!===============================================================================
! Private routines
!===============================================================================

subroutine read_water_refindex(infilename)

   ! read water refractive index file and set module data

   character*(*), intent(in) :: infilename   ! modal optics filename

   ! Local variables

   integer            :: i, ierr
   type(file_desc_t)  :: ncid              ! pio file handle
   integer            :: did               ! dimension ids
   integer            :: dimlen            ! dimension lengths
   type(var_desc_t)   :: vid               ! variable ids
   real(r8) :: refrwsw(nswbands), refiwsw(nswbands) ! real, imaginary ref index for water visible
   real(r8) :: refrwlw(nlwbands), refiwlw(nlwbands) ! real, imaginary ref index for water infrared
   !----------------------------------------------------------------------------

   ! open file
   call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)

   ! inquire dimensions.  Check that file values match parameter values.

   ierr = pio_inq_dimid(ncid, 'lw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nlwbands) then
      write(iulog,*) 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwbands
      call endrun('read_modal_optics: bad lw_band value')
   endif

   ierr = pio_inq_dimid(ncid, 'sw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nswbands) then
      write(iulog,*) 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswbands
      call endrun('read_modal_optics: bad sw_band value')
   endif

   ! read variables
   ierr = pio_inq_varid(ncid, 'refindex_real_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refrwsw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refiwsw)

   ierr = pio_inq_varid(ncid, 'refindex_real_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refrwlw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refiwlw)

   ! set complex representation of refractive indices as module data
   do i = 1, nswbands
      crefwsw(i)  = cmplx(refrwsw(i), abs(refiwsw(i)))
   end do
   do i = 1, nlwbands
      crefwlw(i)  = cmplx(refrwlw(i), abs(refiwlw(i)))
   end do

   call pio_closefile(ncid)

end subroutine read_water_refindex

!===============================================================================

subroutine modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: sigma_logr_aer  ! geometric standard deviation of number distribution
   real(r8), intent(in)  :: dgnumwet(:,:)   ! aerosol wet number mode diameter (m)
   real(r8), intent(out) :: radsurf(:,:)    ! aerosol surface mode radius
   real(r8), intent(out) :: logradsurf(:,:) ! log(aerosol surface mode radius)
   real(r8), intent(out) :: cheb(:,:,:)

   integer  :: i, k, nc
   real(r8) :: alnsg_amode
   real(r8) :: explnsigma
   real(r8) :: xrad(pcols) ! normalized aerosol radius
   !-------------------------------------------------------------------------------

   alnsg_amode = log(sigma_logr_aer)
   explnsigma = exp(2.0_r8*alnsg_amode*alnsg_amode)

   do k = top_lev, pver
      do i = 1, ncol
         ! convert from number mode diameter to surface area
         radsurf(i,k) = 0.5_r8*dgnumwet(i,k)*explnsigma
         logradsurf(i,k) = log(radsurf(i,k))
         ! normalize size parameter
         xrad(i) = max(logradsurf(i,k),xrmin)
         xrad(i) = min(xrad(i),xrmax)
         xrad(i) = (2._r8*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
         ! chebyshev polynomials
         cheb(1,i,k) = 1._r8
         cheb(2,i,k) = xrad(i)
         do nc = 3, ncoef
            cheb(nc,i,k) = 2._r8*xrad(i)*cheb(nc-1,i,k)-cheb(nc-2,i,k)
         end do
      end do
   end do

end subroutine modal_size_parameters

!===============================================================================

      subroutine binterp(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)

!     bilinear interpolation of table
!
      implicit none
      integer im,jm,km,ncol
      real(r8) table(km,im,jm),xtab(im),ytab(jm),out(pcols,km)
      integer i,ix(pcols),ip1,j,jy(pcols),jp1,k,ic
      real(r8) x(pcols),dx,t(pcols),y(pcols),dy,u(pcols), &
             tu(pcols),tuc(pcols),tcu(pcols),tcuc(pcols)

      if(ix(1).gt.0)go to 30
      if(im.gt.1)then
        do ic=1,ncol
          do i=1,im
            if(x(ic).lt.xtab(i))go to 10
          enddo
   10     ix(ic)=max0(i-1,1)
          ip1=min(ix(ic)+1,im)
          dx=(xtab(ip1)-xtab(ix(ic)))
          if(abs(dx).gt.1.e-20_r8)then
             t(ic)=(x(ic)-xtab(ix(ic)))/dx
          else
             t(ic)=0._r8
          endif
	end do
      else
        ix(:ncol)=1
        t(:ncol)=0._r8
      endif
      if(jm.gt.1)then
        do ic=1,ncol
          do j=1,jm
            if(y(ic).lt.ytab(j))go to 20
          enddo
   20     jy(ic)=max0(j-1,1)
          jp1=min(jy(ic)+1,jm)
          dy=(ytab(jp1)-ytab(jy(ic)))
          if(abs(dy).gt.1.e-20_r8)then
             u(ic)=(y(ic)-ytab(jy(ic)))/dy
             if(u(ic).lt.0._r8.or.u(ic).gt.1._r8)then
                write(iulog,*) 'u,y,jy,ytab,dy=',u(ic),y(ic),jy(ic),ytab(jy(ic)),dy
             endif
          else
            u(ic)=0._r8
          endif
	end do
      else
        jy(:ncol)=1
        u(:ncol)=0._r8
      endif
   30 continue
      do ic=1,ncol
         tu(ic)=t(ic)*u(ic)
         tuc(ic)=t(ic)-tu(ic)
         tcuc(ic)=1._r8-tuc(ic)-u(ic)
         tcu(ic)=u(ic)-tu(ic)
         jp1=min(jy(ic)+1,jm)
         ip1=min(ix(ic)+1,im)
         do k=1,km
            out(ic,k)=tcuc(ic)*table(k,ix(ic),jy(ic))+tuc(ic)*table(k,ip1,jy(ic))   &
               +tu(ic)*table(k,ip1,jp1)+tcu(ic)*table(k,ix(ic),jp1)
	 end do
      enddo
      return
      end subroutine binterp

end module modal_aer_opt
