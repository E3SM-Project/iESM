
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
   use pmgrid 
   use infnan,       only: nan
   use constituents, only: pcnst

   implicit none

   integer, parameter :: ptimelevels = 2  ! number of time levels in the dycore
   integer :: n3   = 2
   integer :: n3m1 = 1

   real(r8), allocatable :: ps(:,:,:)
   real(r8), allocatable :: u3(:,:,:,:)
   real(r8), allocatable :: v3(:,:,:,:)
   real(r8), allocatable :: t3(:,:,:,:)
   real(r8), allocatable :: q3(:,:,:,:,:)
   real(r8), allocatable :: qm1(:,:,:,:)
   real(r8), allocatable :: tarrsld(:,:,:)
   real(r8), allocatable :: parrsld(:,:,:)
   real(r8), allocatable :: etadot(:,:,:)
   
   real(r8), allocatable :: div(:,:,:,:)    ! divergence

   real(r8), allocatable :: urhs(:,:,:)
   real(r8), allocatable :: vrhs(:,:,:)
   real(r8), allocatable :: trhs(:,:,:)
   real(r8), allocatable :: prhs(:,:,:)
   
   real(r8), allocatable :: ql(:,:,:)
   real(r8), allocatable :: qm(:,:,:)
   real(r8), allocatable :: ed1(:,:,:)
   real(r8), allocatable :: tlm1(:,:,:)
   real(r8), allocatable :: tl(:,:,:)
   real(r8), allocatable :: tmm1(:,:,:)
   real(r8), allocatable :: tm(:,:,:)
   real(r8), allocatable :: omga(:,:,:)     ! vertical velocity
   real(r8), allocatable :: dpsl(:,:)       ! longitudinal pressure gradient
   real(r8), allocatable :: dpslm1(:,:)     ! longitudinal pressure gradient
   real(r8), allocatable :: dpslp1(:,:)     ! longitudinal pressure gradient
   real(r8), allocatable :: dpsm(:,:)       ! meridional pressure gradient
   real(r8), allocatable :: dpsmm1(:,:)     ! meridional pressure gradient
   real(r8), allocatable :: dpsmp1(:,:)     ! meridional pressure gradient
   real(r8), allocatable :: dps(:,:)        ! pressure gradient
   real(r8), allocatable :: phis(:,:)       ! surface geopotential
   real(r8), allocatable :: phisl(:,:)      ! surface geopotential
   real(r8), allocatable :: phism(:,:)      ! surface geopotential
   
CONTAINS

  subroutine initialize_prognostics
!
! Purpose: Allocate and initialize the prgnostic variables
!
    allocate (ps(plon                 ,beglat:endlat    ,ptimelevels))
    allocate (u3(plon,plev            ,beglat:endlat,ptimelevels))
    allocate (v3(plon,plev            ,beglat:endlat,ptimelevels))
    allocate (t3(plon,plev            ,beglat:endlat,ptimelevels))
    allocate (q3(plon,plev,pcnst,beglat:endlat,ptimelevels))
    allocate (qm1(plon,plev,pcnst,beglat:endlat))
    allocate (etadot(plon,plevp,beglat:endlat))
    allocate (tarrsld(plon,plev,beglat:endlat))
    allocate (parrsld(plon,plev,beglat:endlat))

    allocate (div   (plon,plev,beglat:endlat,ptimelevels))

    allocate (urhs(plon,plev,beglat:endlat))
    allocate (vrhs(plon,plev,beglat:endlat))
    allocate (trhs(plon,plev,beglat:endlat))
    allocate (prhs(plon,plev,beglat:endlat))

    allocate (ql     (plon,plev,beglat:endlat))    
    allocate (qm     (plon,plev,beglat:endlat))    
    allocate (ed1    (plon,plev,beglat:endlat))  
    allocate (tlm1   (plon,plev,beglat:endlat))   
    allocate (tl     (plon,plev,beglat:endlat))   
    allocate (tmm1   (plon,plev,beglat:endlat))   
    allocate (tm     (plon,plev,beglat:endlat))   
    allocate (omga   (plon,plev,beglat:endlat))    
    allocate (dpsl   (plon,beglat:endlat))        
    allocate (dpslm1 (plon,beglat:endlat))        
    allocate (dpslp1 (plon,beglat:endlat))        
    allocate (dpsm   (plon,beglat:endlat))        
    allocate (dpsmm1 (plon,beglat:endlat))        
    allocate (dpsmp1 (plon,beglat:endlat))        
    allocate (dps    (plon,beglat:endlat))         
    allocate (phis   (plon,beglat:endlat))        
    allocate (phisl  (plon,beglat:endlat))        
    allocate (phism  (plon,beglat:endlat))        

    ps(:,:,:)     = nan
    u3(:,:,:,:)   = nan
    v3(:,:,:,:)   = nan
    t3(:,:,:,:)   = nan
    q3(:,:,:,:,:) = nan
    qm1(:,:,:,:) = nan
    tarrsld(:,:,:) = nan
    parrsld(:,:,:) = nan
    etadot(:,:,:) = nan

    div(:,:,:,:) = nan

    urhs(:,:,:) = nan
    vrhs(:,:,:) = nan
    trhs(:,:,:) = nan
    prhs(:,:,:) = nan

    ql     (:,:,:) = nan
    qm     (:,:,:) = nan
    ed1    (:,:,:) = nan
    tlm1   (:,:,:) = nan
    tl     (:,:,:) = nan
    tmm1   (:,:,:) = nan
    tm     (:,:,:) = nan
    omga   (:,:,:) = nan
    dpsl   (:,:) = nan
    dpslm1 (:,:) = nan
    dpslp1 (:,:) = nan
    dpsm   (:,:) = nan
    dpsmm1 (:,:) = nan
    dpsmp1 (:,:) = nan
    dps    (:,:) = nan
    phis   (:,:) = nan
    phisl  (:,:) = nan
    phism  (:,:) = nan

    return
  end subroutine initialize_prognostics

  subroutine shift_time_indices
!
! Purpose: Shift the time indices that track the current and previous times.
!
     integer :: itmp

     itmp = n3m1
     n3m1 = n3
     n3   = itmp
   end subroutine shift_time_indices

end module prognostics
