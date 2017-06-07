module gamma_function

! Double precision version of gamma function.
! Upper incomplete gamma function ("igamma").

! This approximation if the gamma function is both fast and very
! accurate for an argument x such that 0 <= x <= 12. Above 12,
! the accuracy gradually decreases.

! Some comments and warnings:
! 
! 1) "gamma" here is pure, but not elemental. This allows it to be
!    passed as an argument to MG. No CAM code requires igamma as an
!    argument, so it is declared elemental for convenience.
! 
! 2) Fortran 2008 provides an elemental function named "gamma" as
!    well. A module or routine that uses gamma from here, but then
!    calls gamma with a TKR incompatible argument (e.g. a vector
!    argument) may actually end up calling the intrinsic, if the
!    compiler implements it. Compilers that do not implement the
!    intrinsic will simply throw an error in this case.
! 
! 3) "gamma" is believed to be safe from overflow and underflow, and
!    will return the largest possible value in case of overflow. It
!    also handles negative arguments. "igamma" also takes measures to
!    avoid overflow, but does not appear to have been written with
!    negative arguments in mind.

implicit none
private
save

! Public functions
public :: gamma
public :: igamma

! Private variables
integer, parameter :: r4 = selected_real_kind(6)  ! 4 byte real
integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real

! i and square root of pi
real(r8), parameter :: pi = 3.1415926535897932384626434E0_r8
real(r8), parameter :: sqrtpi = 0.9189385332046727417803297E0_r8

! "Data" statements must be outside of routine for it to be pure.

real(r8) :: C(7), P(8), Q(8)

!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
DATA P/-1.71618513886549492533811E+0_r8,2.47656508055759199108314E+1_r8,&
     -3.79804256470945635097577E+2_r8,6.29331155312818442661052E+2_r8,&
     8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8,&
     -3.61444134186911729807069E+4_r8,6.64561438202405440627855E+4_r8/
DATA Q/-3.08402300119738975254353E+1_r8,3.15350626979604161529144E+2_r8,&
     -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8,&
     2.25381184209801510330112E+4_r8,4.75584627752788110767815E+3_r8,&
     -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
DATA C/-1.910444077728E-03_r8,8.4171387781295E-04_r8, &
     -5.952379913043012E-04_r8,7.93650793500350248E-04_r8,&
     -2.777777777777681622553E-03_r8,8.333333333333333331554247E-02_r8,&
     5.7083835261E-03_r8/
!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
real(r8), parameter :: eps = epsilon(pi)
real(r8), parameter :: xinf = huge(pi)
! Assume IEEE double precision and hard-code approximate values for
! these.
real(r8), parameter :: xbig = 171.624_r8
real(r8), parameter :: xminin = 2.23E-308_r8

!------------------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------------------
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

pure function gamma(X)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------

  real(r8), intent(in) :: x
  real(r8) :: gamma
  real(r8) :: fact, res, sum, xden, xnum, y, y1, ysq, z

  integer :: i, n
  logical :: negative_odd

  negative_odd = .false.
  fact = 1._r8
  n = 0
  y = x
  if (y <= 0._r8) then
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
     y = -x
     y1 = aint(y)
     res = y - y1
     if (res /= 0._r8) then
        negative_odd = (y1 /= aint(y1*0.5_r8)*2._r8)
        fact = -pi/sin(pi*res)
        y = y + 1._r8
     else
        gamma = xinf
        return
     end if
  end if
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
  if (y < EPS) then
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
     if (y >= xminin) then
        res = 1._r8/y
     else
        gamma = xinf
        return
     end if
  elseif (y < 12._r8) then
     y1 = y
     if (y < 1._r8) then
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
        z = y
        y = y + 1._r8
     else
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
        n = int(y) - 1
        y = y - real(n, r8)
        z = y - 1._r8
     end if
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
     xnum = 0._r8
     xden = 1._r8
     do i=1,8
        xnum = (xnum+P(i))*z
        xden = xden*z + Q(i)
     end do
     res = xnum/xden + 1._r8
     if (y1 < y) then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
        res = res/y1
     elseif (y1 > y) then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
        do i = 1,n
           res = res*y
           y = y + 1._r8
        end do
     end if
  else
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
     if (y <= xbig) then
        ysq = y*y
        sum = C(7)
        do i=1,6
           sum = sum/ysq + C(i)
        end do
        sum = sum/y - y + sqrtpi
        sum = sum + (y-0.5_r8)*log(y)
        res = exp(sum)
     else
        gamma = xinf
        return
     end if
  end if
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
  if (negative_odd)  res = -res
  if (fact /= 1._r8) res = fact/res
  gamma = res
! ---------- LAST LINE OF GAMMA ----------
end function gamma

!! Incomplete Gamma function
!!
!! @author  Tianyi Fan
!! @version August-2010
real(r8) elemental function igamma(a, x)
  ! Upper incomplete gamma function.
  ! Modified for inclusion in this module and made
  ! pure elemental, September 2012
  
  real(r8), intent(in) ::      a
  real(r8), intent(in) ::      x
  
  ! local variable
  real(r8) :: xam, gin, s, r, t0
  integer  :: k
    
    
  xam = -x + a * log(x)
  
  if ((xam > 700.0_r8) .or. (a > xbig)) then
     ! Out of bounds
     ! Return "huge" value.
     igamma = xinf
     return

  else if (x == 0.0_r8) then           
     igamma = gamma(a)
  
  else if (x <= (1.0_r8 + a)) then
     s = 1.0_r8 / a
     r = s
        
     do  k = 1,60
        r = r * x / (a+k)
        s = s + r
           
        if (abs(r/s) < 1.0e-15_r8) exit
     end do
        
     gin = exp(xam) * s           
     igamma = gamma(a) - gin
        
  else
     t0 = 0.0_r8
    
     do k = 60,1,-1
        t0 = (k - a) / (1.0_r8 + k / (x + t0))
     end do
  
     igamma = exp(xam) / (x + t0)
  endif
    
end function igamma

end module gamma_function
