
module prognostics

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Prognostic variables held in-core for convenient access.
! q3 is specific humidity (water vapor) and other constituents.
! 
! Author: G. Grant
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev, beglat, endlat
   use infnan,       only: inf
   use constituents, only: pcnst


   implicit none

   private

   public ps, u3, v3, t3, q3, qminus, vort, div, dpsl, dpsm, dps, omga, phis, hadv, pdeld
   public n3, n3m1, n3m2, ptimelevels
   public initialize_prognostics
   public shift_time_indices

   integer, parameter :: ptimelevels = 3  ! number of time levels in the dycore
   integer :: n3   = 3
   integer :: n3m1 = 2
   integer :: n3m2 = 1

   real(r8), allocatable :: ps(:,:,:)
   real(r8), allocatable :: u3(:,:,:,:)
   real(r8), allocatable :: v3(:,:,:,:)
   real(r8), allocatable :: t3(:,:,:,:)
   real(r8), allocatable :: pdeld(:,:,:,:)
   real(r8), allocatable :: q3(:,:,:,:,:)
   real(r8), allocatable :: qminus(:,:,:,:)
   real(r8), allocatable :: hadv  (:,:,:,:)

   real(r8), allocatable :: vort(:,:,:,:)   ! vorticity
   real(r8), allocatable :: div(:,:,:,:)    ! divergence

   real(r8), allocatable :: dpsl(:,:)       ! longitudinal pressure gradient
   real(r8), allocatable :: dpsm(:,:)       ! meridional pressure gradient
   real(r8), allocatable :: dps(:,:)        ! pressure gradient
   real(r8), allocatable :: phis(:,:)       ! surface geopotential
   real(r8), allocatable :: omga(:,:,:)     ! vertical velocity

CONTAINS

   subroutine initialize_prognostics
!
! Purpose:  Allocate and initialize the prognostic arrays.
!

      allocate (ps    (plon           ,beglat:endlat    ,ptimelevels))
      allocate (u3    (plon,plev      ,beglat:endlat,ptimelevels))
      allocate (v3    (plon,plev      ,beglat:endlat,ptimelevels))
      allocate (t3    (plon,plev      ,beglat:endlat,ptimelevels))
      allocate (q3    (plon,plev,pcnst,beglat:endlat,ptimelevels))
      allocate (qminus(plon,plev,pcnst,beglat:endlat  ))
      allocate (hadv  (plon,plev,pcnst,beglat:endlat  ))

      allocate (vort  (plon,plev,beglat:endlat,ptimelevels))   
      allocate (div   (plon,plev,beglat:endlat,ptimelevels))    

      allocate (dpsl  (plon,beglat:endlat))        
      allocate (dpsm  (plon,beglat:endlat))        
      allocate (dps   (plon,beglat:endlat))         
      allocate (phis  (plon,beglat:endlat))        
      allocate (omga  (plon,plev,beglat:endlat))    
      allocate (pdeld (plon,plev,beglat:endlat,ptimelevels))

      ps(:,:,:)       = inf
      u3(:,:,:,:)     = inf
      v3(:,:,:,:)     = inf
      t3(:,:,:,:)     = inf
      pdeld(:,:,:,:)  = inf
      q3(:,:,:,:,:)   = inf
      qminus(:,:,:,:) = inf
      hadv  (:,:,:,:) = inf

      vort(:,:,:,:)   = inf
      div (:,:,:,:)   = inf

      dpsl  (:,:) = inf
      dpsm  (:,:) = inf
      dps   (:,:) = inf
      phis  (:,:) = inf
      omga  (:,:,:) = inf

      return
   end subroutine initialize_prognostics

   subroutine shift_time_indices
!
! Purpose: 
! Shift the indices that keep track of which index stores
! the relative times (current time, previous, time before previous etc).
!
      integer :: itmp

      itmp = n3m2

      n3m2 = n3m1
      n3m1 = n3
      n3   = itmp
   end subroutine shift_time_indices

end module prognostics
