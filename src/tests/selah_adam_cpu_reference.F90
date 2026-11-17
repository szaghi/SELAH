!< SELAH tests: ADAM's CPU Runge–Kutta loops, transcribed, the oracle of the CPU envelope test.
module selah_adam_cpu_reference
!< SELAH tests: ADAM's CPU Runge–Kutta loops, transcribed, the oracle of the CPU envelope test.
!<
!< ADAM's `rk_object` methods (src/lib/common/adam_rk_object.F90 at ff879596: initialize_stages 516-556, assign_stage
!< 188-240, compute_stage 242-302, compute_stage_ls 304-403, update_q 574-648) need ADAM's field, grid and MPI
!< objects, so they cannot be copied whole. Their loops are transcribed here: the arithmetic statements and the loop
!< order are verbatim; adapted are `self%q_rk`/`self%alph`/`self%beta`/`self%ark`... (now arguments), the sizes (now
!< arguments instead of `field`/`self` components), the OpenMP `shared` lists (the arrays named instead of `self`), and
!< the first-use zeroing of `initialize_stages` (the caller allocates the stages zeroed). Layout V_FIRST, (v,i,j,k,b).
use penf, only : I4P, R8P
implicit none
private
public :: adam_cpu_initialize_stages
public :: adam_cpu_assign_stage
public :: adam_cpu_compute_stage
public :: adam_cpu_compute_stage_ls
public :: adam_cpu_update_q

contains
   subroutine adam_cpu_initialize_stages(ni, nj, nk, ngc, nv, blocks_number, q, q_rk)
   integer(I4P), intent(in)    :: ni, nj, nk, ngc, nv, blocks_number
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) :: q_rk(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:)
   integer(I4P)                :: i, j, k, b, v, s

   !$omp parallel do collapse(6) default(firstprivate) shared(q,q_rk)
   do s=lbound(q_rk,dim=6),ubound(q_rk,dim=6)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     q_rk(v,i,j,k,b,s) = q(v,i,j,k,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine adam_cpu_initialize_stages

   subroutine adam_cpu_assign_stage(ni, nj, nk, ngc, nv, blocks_number, s, q, q_rk, phi)
   integer(I4P), intent(in)           :: ni, nj, nk, ngc, nv, blocks_number, s
   real(R8P),    intent(in)           :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout)        :: q_rk(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:)
   real(R8P),    intent(in), optional :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P)                       :: all_solids
   integer(I4P)                       :: i, j, k, b, v

   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(5) default(firstprivate) shared(phi,q,q_rk)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     if (phi(all_solids,i,j,k,b) < 0._R8P) then
                        q_rk(v,i,j,k,b,s) = q(v,i,j,k,b)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   else
      !$omp parallel do collapse(5) default(firstprivate) shared(q,q_rk)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     q_rk(v,i,j,k,b,s) = q(v,i,j,k,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endif
   endsubroutine adam_cpu_assign_stage

   subroutine adam_cpu_compute_stage(ni, nj, nk, ngc, nv, blocks_number, s, dt, alph, q_rk, phi)
   integer(I4P), intent(in)           :: ni, nj, nk, ngc, nv, blocks_number, s
   real(R8P),    intent(in)           :: dt
   real(R8P),    intent(in)           :: alph(1:,1:)
   real(R8P),    intent(inout)        :: q_rk(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:)
   real(R8P),    intent(in), optional :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P)                       :: all_solids
   integer(I4P)                       :: i, j, k, b, v, ss

   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(5) default(firstprivate) shared(phi,q_rk,alph)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     do ss=1, s-1
                        if (phi(all_solids,i,j,k,b) < 0._R8P) then
                           q_rk(v,i,j,k,b,s) = q_rk(v,i,j,k,b,s) + dt * alph(s,ss) * q_rk(v,i,j,k,b,ss)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   else
      !$omp parallel do collapse(5) default(firstprivate) shared(q_rk,alph)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     do ss=1, s-1
                        q_rk(v,i,j,k,b,s) = q_rk(v,i,j,k,b,s) + dt * alph(s, ss) * q_rk(v,i,j,k,b,ss)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endif
   endsubroutine adam_cpu_compute_stage

   subroutine adam_cpu_compute_stage_ls(ni, nj, nk, ngc, nv, blocks_number, s, dt, ark, brk, crk, q_rk, dq, q, phi)
   integer(I4P), intent(in)           :: ni, nj, nk, ngc, nv, blocks_number, s
   real(R8P),    intent(in)           :: dt
   real(R8P),    intent(in)           :: ark(1:), brk(1:), crk(1:)
   real(R8P),    intent(in)           :: q_rk(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:)
   real(R8P),    intent(in)           :: dq(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout)        :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(in), optional :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P)                       :: all_solids
   integer(I4P)                       :: i, j, k, b, v

   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(5) default(firstprivate) shared(phi,q,dq,q_rk,ark,brk,crk)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     if (phi(all_solids,i,j,k,b) < 0._R8P) then
                        q(v,i,j,k,b) = ark(s) * q_rk(v,i,j,k,b,1) + &
                                       brk(s) * q(v,i,j,k,b) + dt * crk(s) * dq(v,i,j,k,b)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   else
      !$omp parallel do collapse(5) default(firstprivate) shared(q,dq,q_rk,ark,brk,crk)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     q(v,i,j,k,b) = ark(s) * q_rk(v,i,j,k,b,1) + &
                                    brk(s) * q(v,i,j,k,b) + dt * crk(s) * dq(v,i,j,k,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endif
   endsubroutine adam_cpu_compute_stage_ls

   subroutine adam_cpu_update_q(ni, nj, nk, ngc, nv, blocks_number, nrk, dt, beta, q_rk, q, dq, phi)
   integer(I4P), intent(in)              :: ni, nj, nk, ngc, nv, blocks_number, nrk
   real(R8P),    intent(in)              :: dt
   real(R8P),    intent(in)              :: beta(1:)
   real(R8P),    intent(in)              :: q_rk(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:)
   real(R8P),    intent(inout)           :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout), optional :: dq(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(in),    optional :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P)                          :: all_solids
   integer(I4P)                          :: i, j, k, b, v, s

   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(5) default(firstprivate) shared(phi,q,q_rk,beta)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     do s=1, nrk
                        if (phi(all_solids,i,j,k,b) < 0._R8P) then
                           q(v,i,j,k,b) = q(v,i,j,k,b) + dt * beta(s) * q_rk(v,i,j,k,b,s)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   else
      if (present(dq)) dq(:,:,:,:,1:blocks_number) = 0._R8P
      !$omp parallel do collapse(5) default(firstprivate) shared(q,dq,q_rk,beta)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     do s=1, nrk
                        if (present(dq)) dq(v,i,j,k,b) = dq(v,i,j,k,b) + beta(s)*q_rk(v,i,j,k,b,s)
                        q(v,i,j,k,b) = q(v,i,j,k,b) + dt * beta(s) * q_rk(v,i,j,k,b,s)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endif
   endsubroutine adam_cpu_update_q
endmodule selah_adam_cpu_reference
