!< SELAH, L1: the two-level particle pushers (Boris, relativistic Boris, Vay, Higuera–Cary).

#include "fundal.H"

module selah_pushers
!< SELAH, L1: the two-level particle pushers (Boris, relativistic Boris, Vay, Higuera–Cary).
!<
!< A staggered leapfrog for charged particles in given fields: the positions at t^n, the velocities (or the momenta
!< per unit mass u = gamma v) at t^(n-1/2), the fields E and B gathered by the caller at the positions x^n. One call
!< advances every particle by one step dt: v (or u) to t^(n+1/2), then x to t^(n+1). These are not Runge–Kutta stage
!< schemes (no tableau, no registers): a kernel of their own. The theory and the sources: `notes/companions-pushers-pml.md`.
!<
!<| Pusher | State | Defining equation (u bar = (u^{n+1/2} + u^{n-1/2})/2) | Properties |
!<|--------|-------|--------------------------------------------------------|------------|
!<| `PUSHER_BORIS` | v | (v^+ - v^-)/dt = (q/m)(E + v bar x B) | volume-preserving; \|v\| exact in pure B; ExB exact |
!<| `PUSHER_BORIS_RELATIVISTIC` | u | (u^+ - u^-)/dt = (q/m)(E + u bar x B / gamma(u^- + qh E)) | volume-preserving; wrong ExB velocity, fails E + v x B = 0 |
!<| `PUSHER_VAY` | u | (u^+ - u^-)/dt = (q/m)(E + (v^+ + v^-)/2 x B) | exact for E + v x B = 0; not volume-preserving |
!<| `PUSHER_HIGUERA_CARY` | u | (u^+ - u^-)/dt = (q/m)(E + u bar x B / gamma(u bar)) | volume-preserving; exact for E + v x B = 0 and the ExB drift |
!<
!< SI units: F = q (E + v x B), with q/m per particle from two components of the particle array (`iq`, `im`); c the
!< speed of light (relativistic pushers; a Gaussian-units caller passes B/c). The particle array holds a particle per
!< column (`LAYOUT_V_FIRST`, `p(ncomp, np)`, ADAM's host layout) or per row (`LAYOUT_V_LAST`, `p(np, ncomp)`, its device
!< layout); the positions at components ix..ix+2, the velocities or momenta at iv..iv+2. The fields `e`, `b` hold
!< three components per particle in the same layout. The first n particles are pushed; the rest are untouched.
use penf, only : I4P
use selah_kind_parameters, only : R_P
use selah_box, only : LAYOUT_V_FIRST, LAYOUT_V_LAST
implicit none
private
public :: PUSHER_BORIS, PUSHER_BORIS_RELATIVISTIC, PUSHER_VAY, PUSHER_HIGUERA_CARY
public :: selah_push

integer(I4P), parameter :: PUSHER_BORIS              = 1_I4P !< Boris (1970), non-relativistic.
integer(I4P), parameter :: PUSHER_BORIS_RELATIVISTIC = 2_I4P !< Boris, relativistic (u = gamma v).
integer(I4P), parameter :: PUSHER_VAY                = 3_I4P !< Vay (2008).
integer(I4P), parameter :: PUSHER_HIGUERA_CARY       = 4_I4P !< Higuera & Cary (2017).

contains
   subroutine selah_push(pusher, layout, n, dt, ix, iv, iq, im, p, e, b, c)
   !< Push the first n particles of p by one step dt in the fields e, b (gathered at their positions).
   integer(I4P), intent(in)           :: pusher   !< PUSHER_*.
   integer(I4P), intent(in)           :: layout   !< LAYOUT_V_FIRST (p(ncomp, np)) or LAYOUT_V_LAST (p(np, ncomp)).
   integer(I4P), intent(in)           :: n        !< Particles pushed (the first n).
   real(R_P),    intent(in)           :: dt       !< Time step.
   integer(I4P), intent(in)           :: ix, iv   !< First components of the position and of the velocity (momentum).
   integer(I4P), intent(in)           :: iq, im   !< Components of the charge and of the mass.
   real(R_P),    intent(inout)        :: p(1:,1:) !< Particles.
   real(R_P),    intent(in)           :: e(1:,1:) !< Electric field at the particles, (3, np) or (np, 3).
   real(R_P),    intent(in)           :: b(1:,1:) !< Magnetic field at the particles, (3, np) or (np, 3).
   real(R_P),    intent(in), optional :: c        !< Speed of light (relativistic pushers; default 1).
   real(R_P)                          :: ci, c2i  !< 1/c, 1/c^2.
   real(R_P)                          :: qh, tt, g, us, sg !< Per particle: q dt/(2m), scalars of the push.
   real(R_P)                          :: m1, m2, m3, t1, t2, t3, s1, s2, s3, w1, w2, w3 !< Per particle: vectors.
   integer(I4P)                       :: i        !< Particle.

   ci = 1._R_P ; if (present(c)) ci = 1._R_P / c
   c2i = ci * ci
   if (pusher < PUSHER_BORIS .or. pusher > PUSHER_HIGUERA_CARY) error stop 'selah push: unknown pusher'
   if (layout == LAYOUT_V_LAST) then
#define PP(k) p(i,k)
#define EE(k) e(i,k)
#define BB(k) b(i,k)
#include "selah_kernel_push.INC"
#undef PP
#undef EE
#undef BB
   else
#define PP(k) p(k,i)
#define EE(k) e(k,i)
#define BB(k) b(k,i)
#include "selah_kernel_push.INC"
#undef PP
#undef EE
#undef BB
   endif
   endsubroutine selah_push
endmodule selah_pushers
