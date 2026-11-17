!< SELAH verification A5: the linear stability of every Runge–Kutta scheme, from the library's own tableaux.
program selah_vv_stability
!< SELAH verification A5: the linear stability of every Runge–Kutta scheme, from the library's own tableaux.
!<
!< For y' = lambda y a step multiplies y by R(z), z = lambda dt. For an explicit tableau R is the polynomial
!< R(z) = 1 + sum_{j=1..s} b^T A^{j-1} 1 z^j, computed here from butcher(scheme) in quadruple precision (QP; R8P where
!< the compiler has no quad kind). Asserted:
!<+ the first p+1 coefficients are those of exp(z), 1/j! (the order p, through the tall trees);
!<+ the limits on the negative real axis and on the imaginary axis match the closed forms where they are known:
!<  euler and ssp-22 -2 and 0; ssp-33 -2.512745 and sqrt(3); rk4 -2.785294 and 2 sqrt(2);
!<+ an SSP scheme with coefficient C is stable on the disc |z + C| <= C (checked on 720 points of its boundary), so
!<  its real-axis limit is at most -2C.
!< Report mode (SELAH_VV_DIR): |R| on a grid around each scheme's region (stability-<scheme>.dat), the coefficients
!< and the limits (stability-summary.dat).
use selah
use selah_test_tools, only : QP, check, failures
use selah_vv_cases
implicit none
real(QP), parameter   :: PI_Q = acos(-1._QP)
real(QP), parameter   :: TOL = max(1.e-9_QP, 1.e3_QP * real(epsilon(1._R_P), QP)) !< R4P coefficients: ~1e-4
real(QP), parameter   :: TOL_AXIS = 64._QP * real(epsilon(1._R_P), QP) !< |R| - 1 allowed on the imaginary axis: the
                                                                       !< rounding of the coefficients
type(butcher_tableau) :: bt
real(QP), allocatable :: coef(:)
real(QP)              :: xr, yi
integer               :: s, j, us
character(:), allocatable :: name

if (vv_report()) call vv_open('stability-summary.dat', us, &
                              'scheme stages order real_limit imag_limit ssp_coefficient coef_0..coef_s (R(z) = sum coef_j z^j)')
do s=1, size(RK_SCHEMES)
   name = trim(RK_SCHEMES(s))
   bt = butcher(name)
   coef = stability_polynomial(bt)
   call check(all([(abs(coef(j+1) - 1._QP / gamma(real(j+1, QP))) < TOL, j=0, bt%order)]), &
              name//': R(z) agrees with exp(z) up to z^p')
   xr = real_limit(coef) ; yi = imag_limit(coef)
   print '(A,2F12.6)', '        '//name//' real and imaginary limits:', real(xr, R8P), real(yi, R8P)
   select case(name)
   case('euler', 'ssp-22')
      ! no interval on the imaginary axis: |R(iy)| - 1 ~ y^4/8 (ssp-22), so the rounding opens one up to (8 TOL_AXIS)^(1/4)
      call check(abs(xr + 2._QP) < TOL .and. yi <= (8._QP * TOL_AXIS)**0.25_QP + 1.e-3_QP, name//': limits -2 and 0')
   case('ssp-33')
      call check(abs(xr + 2.5127453266183286_QP) < TOL .and. abs(yi - sqrt(3._QP)) < TOL, &
                 name//': limits -2.512745 and sqrt(3)')
   case('rk4')
      call check(abs(xr + 2.7852935634761548_QP) < TOL .and. abs(yi - 2._QP*sqrt(2._QP)) < TOL, &
                 name//': limits -2.785294 and 2 sqrt(2)')
   endselect
   if (bt%ssp_coefficient > 0._R_P) then
      call check(disc_stable(coef, real(bt%ssp_coefficient, QP)), name//': stable on the disc |z + C| <= C (SSP)')
      call check(xr <= -2._QP * real(bt%ssp_coefficient, QP) + TOL, name//': real limit <= -2C')
   endif
   if (vv_report()) then
      write(us, '(A,2I4,3ES24.15,*(ES24.15))') name, bt%stages, bt%order, real(xr, R8P), real(yi, R8P), &
            real(bt%ssp_coefficient, R8P), real(coef, R8P)
      call write_grid(name, coef, xr, yi)
   endif
enddo
if (vv_report()) close(us)
if (failures > 0) error stop 'selah_vv_stability: failed'
print '(A)', 'selah_vv_stability: ok'

contains
   function stability_polynomial(bt) result(c)
   !< c(j+1) = b^T A^{j-1} 1 (c(1) = 1): R(z) = sum_j c(j+1) z^j, in QP.
   type(butcher_tableau), intent(in) :: bt
   real(QP), allocatable             :: c(:), a(:,:), v(:)
   integer                           :: j

   a = real(bt%a, QP)
   allocate(c(bt%stages+1)) ; c(1) = 1._QP
   v = [(1._QP, j=1, bt%stages)]
   do j=1, bt%stages
      c(j+1) = sum(real(bt%b, QP) * v)
      v = matmul(a, v)
   enddo
   endfunction stability_polynomial

   pure function absr(c, z) result(r)
   !< |R(z)|, Horner.
   real(QP),    intent(in) :: c(:)
   complex(QP), intent(in) :: z
   real(QP)                :: r
   complex(QP)             :: p
   integer                 :: j

   p = cmplx(c(size(c)), 0._QP, QP)
   do j=size(c)-1, 1, -1
      p = p * z + c(j)
   enddo
   r = abs(p)
   endfunction absr

   function real_limit(c) result(x)
   !< The left end of the interval of the negative real axis that contains 0 and where |R| <= 1.
   real(QP), intent(in) :: c(:)
   real(QP)             :: x, lo, hi, mid
   integer              :: k

   x = 0._QP
   do while (absr(c, cmplx(x - 1.e-3_QP, 0._QP, QP)) <= 1._QP .and. x > -100._QP)
      x = x - 1.e-3_QP
   enddo
   lo = x - 1.e-3_QP ; hi = x                          ! |R(lo)| > 1, |R(hi)| <= 1
   do k=1, 200
      mid = 0.5_QP * (lo + hi)
      if (absr(c, cmplx(mid, 0._QP, QP)) <= 1._QP) then ; hi = mid ; else ; lo = mid ; endif
   enddo
   x = hi
   endfunction real_limit

   function imag_limit(c) result(y)
   !< The top end of the interval of the imaginary axis that contains 0 and where |R| <= 1 + TOL_AXIS (0 if none): the
   !< tolerance of the coefficients' rounding (rk4 is tangent to |R| = 1 at 0, and R4P coefficients cross it by ~1e-8).
   real(QP), intent(in) :: c(:)
   real(QP)             :: y, lo, hi, mid
   integer              :: k

   y = 0._QP
   do while (absr(c, cmplx(0._QP, y + 1.e-3_QP, QP)) <= 1._QP + TOL_AXIS .and. y < 100._QP)
      y = y + 1.e-3_QP
   enddo
   if (y == 0._QP) return
   lo = y ; hi = y + 1.e-3_QP
   do k=1, 200
      mid = 0.5_QP * (lo + hi)
      if (absr(c, cmplx(0._QP, mid, QP)) <= 1._QP + TOL_AXIS) then ; lo = mid ; else ; hi = mid ; endif
   enddo
   y = lo
   endfunction imag_limit

   function disc_stable(c, cssp) result(ok)
   !< |R| <= 1 on the boundary of the disc |z + C| <= C (720 points; the maximum principle covers the inside).
   real(QP), intent(in) :: c(:), cssp
   logical              :: ok
   integer              :: k

   ok = .true.
   do k=0, 719
      if (absr(c, cmplx(-cssp + cssp*cos(2*PI_Q*k/720), cssp*sin(2*PI_Q*k/720), QP)) > 1._QP + TOL) ok = .false.
   enddo
   endfunction disc_stable

   subroutine write_grid(name, c, xr, yi)
   !< |R| on a grid covering the region: x from 1.15 xr - 0.5 to 1, y symmetric beyond the region, 241 x 241.
   character(*), intent(in) :: name
   real(QP),     intent(in) :: c(:), xr, yi
   real(QP)                 :: x0, x1, y1, x, y
   integer                  :: u, i, k

   x0 = 1.15_QP * xr - 0.5_QP ; x1 = 1._QP ; y1 = max(0.75_QP * (x1 - x0), 1.6_QP * yi + 1._QP)
   call vv_open('stability-'//name//'.dat', u, 'x y |R(x + i y)|, 241 x 241, x fastest')
   do k=0, 240
      y = -y1 + 2._QP * y1 * k / 240
      do i=0, 240
         x = x0 + (x1 - x0) * i / 240
         write(u, '(3ES14.6)') real(x, R8P), real(y, R8P), real(min(absr(c, cmplx(x, y, QP)), 1.e3_QP), R8P)
      enddo
   enddo
   close(u)
   endsubroutine write_grid
endprogram selah_vv_stability
