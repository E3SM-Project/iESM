
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
  use infnan,       only: bigint
  use pmgrid,       only: plev, plat
  use pspect,       only: pmmax, pspt

  implicit none

  private

!
! $Id$
! $Author$
!

  public vz, d, t, alps, numm, maxm, lpspt, locm, locrm, lnstart
  public nstart, nlen, alp, dalp, lalp, ldalp

! Spectral space arrays
!
  real(r8), dimension(:,:), allocatable :: vz   ! Vorticity spectral coefficients
  real(r8), dimension(:,:), allocatable :: d    ! Divergence spectral coefficients
  real(r8), dimension(:,:), allocatable :: t    ! Temperature spectral coefficients
  real(r8), dimension(:), allocatable :: alps ! Log-pressure spectral coefficients

#if ( defined SPMD )
  integer :: maxm           = bigint  ! max number of Fourier wavenumbers per MPI task
  integer :: lpspt          = bigint  ! number of local spectral coefficients
  integer, dimension(:), allocatable   :: numm
                                      ! number of Fourier wavenumbers owned per task
  integer, dimension(:,:), allocatable :: locm, locrm
                                      ! assignment of wavenumbers to MPI tasks
  integer, dimension(:), allocatable   :: lnstart 
                                      ! Starting indices for local spectral arrays (real)
#else
  integer :: numm(0:0)      = pmmax
  integer :: maxm           = pmmax
  integer :: lpspt          = pspt
  integer :: locm(1:pmmax, 0:0) = bigint
  integer :: locrm(1:2*pmmax, 0:0) = bigint
  integer :: lnstart(1:pmmax) = bigint
#endif

  integer :: nstart(pmmax) = bigint   ! Starting indices for spectral arrays (real)
  integer :: nlen(pmmax)   = bigint   ! Length vectors for spectral arrays

  real(r8), dimension(:,:), allocatable :: alp  ! Legendre polynomials (pspt,plat/2)
  real(r8), dimension(:,:), allocatable :: dalp ! Legendre polynomial derivatives (pspt,plat/2)
!
  real(r8), dimension(:,:), allocatable :: lalp  ! local Legendre polynomials
  real(r8), dimension(:,:), allocatable :: ldalp ! local Legendre polynomial derivatives

end module comspe
