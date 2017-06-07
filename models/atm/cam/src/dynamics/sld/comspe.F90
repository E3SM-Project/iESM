
module comspe

!----------------------------------------------------------------------- 
! 
! Purpose: Spectral space arrays
! 
! Method: 
! 
! Author: CCM Core Group
! $Author$
! $Id$
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use infnan, only: nan, bigint
  use pmgrid, only: plev, plat
  use pspect, only: psp, pmax, pmmax, pspt, pmaxp

  implicit none

  private

  public vz, d, t, q, alps, hs, hsnm, dsnm, dnm, vznm, lnpstar
  public a0nm, bmnm, bpnm, atri
  public btri, ctri, ncutoff, nalp, numm, maxm, lpspt, locm, locrm, lnstart
  public ncoefi, nm, nco2, nstart, nlen, alp, dalp

  real(r8) :: vz(psp,plev)   = nan     ! Vorticity spectral coefficients
  real(r8) :: d(psp,plev)    = nan     ! Divergence spectral coefficients
  real(r8) :: t(psp,plev)    = nan     ! Temperature spectral coefficients
  real(r8) :: q(psp,plev)    = nan     ! Moisture     spectral coefficients
  real(r8) :: alps(psp)      = nan     ! Log-pressure spectral coefficients
  real(r8) :: hs(psp,plev)   = nan     ! hydrostatic matrix for "real" atmosphere
  real(r8) :: hsnm(psp,plev) = nan   ! vertical normal modes of "hs"
  real(r8) :: dsnm(psp,plev) = nan   ! vertical normal modes of "ds"
  real(r8) :: dnm(psp,plev)  = nan   ! vertical normal modes of "d"
  real(r8) :: vznm(psp,plev) = nan   ! vertical normal modes of "vz"
  real(r8) :: lnpstar(psp)   = nan   ! ln (Ps*) (SLD term; Ritchie & Tanguay, 1995)
  real(r8) :: a0nm(psp)      = nan   ! wave # coefs (use in vert normal mode space)
  real(r8) :: bmnm(psp)      = nan   ! wave # coefs (use in vert normal mode space)
  real(r8) :: bpnm(psp)      = nan   ! wave # coefs (use in vert normal mode space)
  real(r8) :: atri(psp)      = nan   ! wave # coefs (use in vert normal mode space)
  real(r8) :: btri(psp)      = nan   ! wave # coefs (use in vert normal mode space)
  real(r8) :: ctri(psp)      = nan   ! wave # coefs (use in vert normal mode space)

  integer :: ncutoff         = bigint   ! Break-even point for vector lengths in GRCALC
  integer :: nalp(pmax)      = bigint   ! Pointer into polynomial arrays
#if ( defined SPMD )
  integer :: maxm            = bigint  ! max number of Fourier wavenumbers per MPI task
  integer :: lpspt           = bigint  ! number of local spectral coefficients (NOT USED YET)
  integer, dimension(:), allocatable   :: numm
                                       ! number of Fourier wavenumbers owned per task
  integer, dimension(:,:), allocatable :: locm, locrm
                                       ! assignment of wavenumbers to MPI tasks
  integer, dimension(:), allocatable   :: lnstart 
                                       ! Starting indices for local spectral arrays (real) (NOT USED YET)
#else
  integer :: numm(0:0)       = pmmax
  integer :: maxm            = pmmax
  integer :: lpspt           = pspt
  integer :: locm(1:pmmax, 0:0) = bigint
  integer :: locrm(1:2*pmmax, 0:0) = bigint
  integer :: lnstart(1:pmmax) = bigint
#endif

  integer :: ncoefi(pmaxp)   = bigint   ! Pointer to start of coefficient diagonals
  integer :: nm(pmax)        = bigint   ! Number of coeffs stored on a given diagonal
  integer :: nco2(pmax)      = bigint   ! Complex form of ncoefi
  integer :: nstart(pmmax)   = bigint   ! Starting indices for spectral arrays (real)
  integer :: nlen(pmmax)     = bigint   ! Length vectors for spectral arrays

  real(r8) :: alp(pspt,plat/2)  = nan   ! Legendre polynomials
  real(r8) :: dalp(pspt,plat/2) = nan   ! Legendre polynomial derivatives

end module comspe
