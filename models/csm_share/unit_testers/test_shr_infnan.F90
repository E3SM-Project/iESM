 program test_shr_infnan
!
! Simple unit-test program for the shr_infnan_mod module.
!
! Erik Kluzek
!
! $Id: test_shr_sys.F90 7482 2007-11-07 20:54:58Z erik $
!
 use shr_kind_mod, only: r8 => shr_kind_r8
 use shr_kind_mod, only: r4 => shr_kind_r4
 use shr_infnan_mod
 use test_mod
  
 implicit none
 real(r8) :: x, zero
 real(r8) :: r8array(100), r82Darray(10,10), r83Darray(4,4,4)
 real(r8) :: r84Darray(3,3,3,3), r85Darray(2,2,2,2,2)
 real(r8), parameter :: inf    = O'0777600000000000000000'
 real(r8), parameter :: nan    = O'0777700000000000000000'
 real(r8), parameter :: nans   = O'0777610000000000000000'
 !real(r4), parameter :: spnan  = O'07770000000'
 real(r4), parameter :: spnan  = Z"7FC00000"
 real(r4), parameter :: spnans = Z"7FC10000"
 integer   :: count
 intrinsic :: count

 call test_init( 29 )

 print *, "Unit-tester for shr_infnan_mod"

 x    = 0.0
 zero = 0.0

 call test_is(       shr_infnan_isnan( nan ),    "Test that value set to nan is nan" )
 call test_is(       shr_infnan_isnan( nans ),   "Test that value set to nans is nan" )
 call test_is(       shr_infnan_isnan( spnan ),  "Test that value set to sp nan is nan" )
 call test_is(       shr_infnan_isnan( spnans ), "Test that value set to sp nans is nan" )
 call test_is( .not. shr_infnan_isnan( 1.0_r8 ), "Test that value set to one is NOT nan" )
 call test_is( .not. shr_infnan_isnan( 1.0_r4 ), "Test that value set to SP one is NOT nan" )
 call test_is( .not. shr_infnan_isnan( huge(x) ), "Test that value set to huge is NOT nan" )
 x    = 1.0/zero
 call test_is( .not. shr_infnan_isnan( x ),     "Test that 1/0 is NOT nan" )
 x    = -1.0/zero
 call test_is( .not. shr_infnan_isnan( x ),    "Test that -1/0 is NOT nan" )

 r8array(:)      = 1.0d00
 r8array(10)     = nan
 r8array(15)     = nan
 r82Darray(:,:)  = 1.0d00
 r82Darray(5,5)  = nan
 r82Darray(10,7) = nan
 r82Darray(7,9)  = nan
 r83Darray(:,:,:) = 1.0d00
 r83Darray(4,2,2) = nan
 r83Darray(3,1,2) = nan
 r83Darray(1,1,1) = nan
 r83Darray(1,1,4) = nan
 r84Darray(:,:,:,:) = 1.0d00
 r84Darray(3,2,2,1) = nan
 r84Darray(3,1,2,1) = nan
 r84Darray(1,1,1,1) = nan
 r84Darray(1,1,3,1) = nan
 r84Darray(1,2,3,1) = nan
 r85Darray(:,:,:,:,:) = 1.0d00
 r85Darray(1,2,2,1,1) = nan
 r85Darray(1,1,2,1,2) = nan
 r85Darray(1,1,1,2,1) = nan
 r85Darray(1,2,2,2,1) = nan
 r85Darray(1,2,1,1,2) = nan
 r85Darray(1,1,1,1,1) = nan
 call test_is(       any(shr_infnan_isnan( r8array )),      "Test that array with 2 nans is nan" )
 call test_is(       count(shr_infnan_isnan( r8array )),2,    "Test that there are 2 nans in that array" )
 call test_is(       any(shr_infnan_isnan( r82Darray )),    "Test that 2D array with 3 nans is nan" )
 call test_is(       count(shr_infnan_isnan( r82Darray )),3,  "Test that there are 3 nans in that array" )
 call test_is(       any(shr_infnan_isnan( r83Darray )),    "Test that 3D array with 4 nans is nan" )
 call test_is(       count(shr_infnan_isnan( r83Darray )),4,  "Test that there are 4 nans in that array" )
 call test_is(       any(shr_infnan_isnan( r84Darray )),    "Test that 4D array with 5 nans is nan" )
 call test_is(       count(shr_infnan_isnan( r84Darray )),5,  "Test that there are 5 nans in that array" )
 call test_is(       any(shr_infnan_isnan( r85Darray )),    "Test that 5D array with 6 nans is nan" )
 call test_is(       count(shr_infnan_isnan( r85Darray )),6,  "Test that there are 6 nans in that array" )
 call test_is(       shr_infnan_isposinf( inf ),   "Test that value set to inf is inf" )
 call test_is( .not. shr_infnan_isposinf( 1.0_r8 ), "Test that value set to one is NOT inf" )
 call test_is( .not. shr_infnan_isposinf( 1.0_r4 ), "Test that value set to SP one is NOT inf" )
 call test_is(       shr_infnan_isneginf( -inf ),   "Test that value set to -inf is -inf" )
 call test_is( .not. shr_infnan_isneginf( 1.0_r8 ), "Test that value set to one is NOT -inf" )
 call test_is( .not. shr_infnan_isneginf( 1.0_r4 ), "Test that value set to SP one is NOT -inf" )
 x = 1.0/zero
 call test_is(       shr_infnan_isposinf( x ),      "Test that 1/0 is inf" )
 x = -1.0/zero
 call test_is(       shr_infnan_isneginf( x ),     "Test that -1/0 is -inf" )

 x = -1.0
 call test_is(       shr_infnan_isnan( sqrt(x) ), "Test that sqrt-1 is nan" )
 call test_is(       shr_infnan_isnan( log(x) ),  "Test that log-1 is nan" )

 call test_final()

 end program test_shr_infnan
