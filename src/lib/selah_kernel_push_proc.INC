   ! One layout of the particle push of SELAH (selah_pushers) as a procedure of its own: SELAH_PROC is its name (undefined at
   ! the end), the macros PP(k), EE(k) and BB(k) select the layout (see selah_pushers). The dummy arguments are those of
   ! selah_push, without the layout and with ci = 1/c and c2i = 1/c^2 computed by the dispatcher.
   subroutine SELAH_PROC(pusher, n, dt, ix, iv, iq, im, p, e, b, ci, c2i)
   integer(I4P), intent(in)    :: pusher   !< PUSHER_*.
   integer(I4P), intent(in)    :: n        !< Particles pushed (the first n).
   real(R_P),    intent(in)    :: dt       !< Time step.
   integer(I4P), intent(in)    :: ix, iv   !< First components of the position and of the velocity (momentum).
   integer(I4P), intent(in)    :: iq, im   !< Components of the charge and of the mass.
   real(R_P),    intent(inout) :: p(1:,1:) !< Particles.
   real(R_P),    intent(in)    :: e(1:,1:) !< Electric field at the particles, (3, np) or (np, 3).
   real(R_P),    intent(in)    :: b(1:,1:) !< Magnetic field at the particles, (3, np) or (np, 3).
   real(R_P),    intent(in)    :: ci, c2i  !< 1/c, 1/c^2.
   real(R_P)                   :: qh, tt, g, us, sg !< Per particle: q dt/(2m), scalars of the push.
   real(R_P)                   :: m1, m2, m3, t1, t2, t3, s1, s2, s3, w1, w2, w3 !< Per particle: vectors.
   integer(I4P)                :: i        !< Particle.

#include "selah_kernel_push.INC"
   endsubroutine SELAH_PROC
#undef SELAH_PROC
