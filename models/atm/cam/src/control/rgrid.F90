
module rgrid

  use pmgrid, only: plat, plon
  use pspect, only: pmmax, pmax, ptrm
  use infnan, only: bigint

  implicit none

  integer :: nlon(plat)        = plon ! num longitudes per latitude
  integer :: beglatpair(pmmax) = bigint
  integer :: nmmax(plat/2)     = bigint
  integer :: wnummax(plat)     = ptrm ! cutoff Fourier wavenumber

  logical :: fullgrid                   ! true => no grid reduction towards poles
end module rgrid
