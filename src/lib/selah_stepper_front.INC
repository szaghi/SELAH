   ! The rank-k fronts of the stepper (selah_stepper), included once per rank 1 to 4 with SELAH_FRONT_RANK, its
   ! assumed-shape dimensions SELAH_FRONT_DIMS, the specific names and SELAH_VIEW_BOUNDS. Each front checks that q has
   ! the field's rank and bounds, then passes q to the explicit-shape bridge: the field as the rank-5 array it is in
   ! memory, no copy when it is contiguous (selah_box: the dimensions as rank 5).
   subroutine SELAH_BEGIN_STEP(self, bx, q, mask, mask_comp)
   !< begin_step for a field of rank k < 5.
   class(stepper), intent(inout)        :: self                    !< Stepper.
   type(box),      intent(in)           :: bx                      !< Box of the field's rank.
   real(R_P),      intent(in)           :: q(SELAH_FRONT_DIMS)     !< Solution.
   real(R_P),      intent(in), optional :: mask(SELAH_FRONT_DIMS)  !< Mask (q's rank, its own components).
   integer(I4P),   intent(in), optional :: mask_comp               !< Mask component.
   integer(I4P)                         :: ubq(5), ubm(5)          !< Upper bounds of q and of the mask, as rank 5.

   call check_front(self, bx, SELAH_FRONT_RANK, shape(q, kind=I4P), ubq)
   ubm = ubq
   if (present(mask)) ubm = front_ub(bx, shape(mask, kind=I4P))
   call begin_step_x(self, bx, ubq, q, ubm, mask, mask_comp)
   endsubroutine SELAH_BEGIN_STEP

   subroutine SELAH_BEGIN_STAGE(self, k, bx, dt, q, mask, mask_comp, decay)
   !< begin_stage for a field of rank k < 5.
   class(stepper), intent(inout)        :: self                    !< Stepper.
   integer(I4P),   intent(in)           :: k                       !< Stage.
   type(box),      intent(in)           :: bx                      !< Box of the field's rank.
   real(R_P),      intent(in)           :: dt                      !< Time step.
   real(R_P),      intent(in)           :: q(SELAH_FRONT_DIMS)     !< Solution.
   real(R_P),      intent(in), optional :: mask(SELAH_FRONT_DIMS)  !< Mask (q's rank, its own components).
   integer(I4P),   intent(in), optional :: mask_comp               !< Mask component.
   real(R_P),      intent(in), optional :: decay(SELAH_FRONT_DIMS) !< Decay rate (q's shape).
   integer(I4P)                         :: ubq(5), ubm(5)          !< Upper bounds of q and of the mask, as rank 5.

   call check_front(self, bx, SELAH_FRONT_RANK, shape(q, kind=I4P), ubq)
   if (present(decay)) then
      if (any(shape(decay) /= shape(q))) error stop 'selah stepper: decay must have the shape of q'
   endif
   ubm = ubq
   if (present(mask)) ubm = front_ub(bx, shape(mask, kind=I4P))
   call begin_stage_x(self, k, bx, dt, ubq, q, ubm, mask, mask_comp, decay)
   endsubroutine SELAH_BEGIN_STAGE

   subroutine SELAH_END_STAGE(self, k, bx, dt, q, mask, mask_comp)
   !< end_stage for a field of rank k < 5.
   class(stepper), intent(inout)        :: self                    !< Stepper.
   integer(I4P),   intent(in)           :: k                       !< Stage.
   type(box),      intent(in)           :: bx                      !< Box of the field's rank.
   real(R_P),      intent(in)           :: dt                      !< Time step.
   real(R_P),      intent(inout)        :: q(SELAH_FRONT_DIMS)     !< Solution.
   real(R_P),      intent(in), optional :: mask(SELAH_FRONT_DIMS)  !< Mask (q's rank, its own components).
   integer(I4P),   intent(in), optional :: mask_comp               !< Mask component.
   integer(I4P)                         :: ubq(5), ubm(5)          !< Upper bounds of q and of the mask, as rank 5.

   call check_front(self, bx, SELAH_FRONT_RANK, shape(q, kind=I4P), ubq)
   ubm = ubq
   if (present(mask)) ubm = front_ub(bx, shape(mask, kind=I4P))
   call end_stage_x(self, k, bx, dt, ubq, q, ubm, mask, mask_comp)
   endsubroutine SELAH_END_STAGE

   subroutine SELAH_FINISH_STEP(self, bx, dt, q, dq, mask, mask_comp, decay)
   !< finish_step for a field of rank k < 5.
   class(stepper), intent(inout)           :: self                   !< Stepper.
   type(box),      intent(in)              :: bx                     !< Box of the field's rank.
   real(R_P),      intent(in)              :: dt                     !< Time step.
   real(R_P),      intent(inout)           :: q(SELAH_FRONT_DIMS)    !< Solution.
   real(R_P),      intent(inout), optional :: dq(SELAH_FRONT_DIMS)   !< Step residual (q's shape).
   real(R_P),      intent(in),    optional :: mask(SELAH_FRONT_DIMS) !< Mask (q's rank, its own components).
   integer(I4P),   intent(in),    optional :: mask_comp              !< Mask component.
   real(R_P),      intent(in),    optional :: decay(SELAH_FRONT_DIMS) !< Decay rate (q's shape).
   integer(I4P)                            :: ubq(5), ubm(5)         !< Upper bounds of q and of the mask, as rank 5.

   call check_front(self, bx, SELAH_FRONT_RANK, shape(q, kind=I4P), ubq)
   if (present(dq)) then
      if (any(shape(dq) /= shape(q))) error stop 'selah stepper: dq must have the shape of q'
   endif
   if (present(decay)) then
      if (any(shape(decay) /= shape(q))) error stop 'selah stepper: decay must have the shape of q'
   endif
   ubm = ubq
   if (present(mask)) ubm = front_ub(bx, shape(mask, kind=I4P))
   call finish_step_x(self, bx, dt, ubq, q, dq, ubm, mask, mask_comp, decay)
   endsubroutine SELAH_FINISH_STEP

   subroutine SELAH_STATE_VIEW(self, u)
   !< The stage state of a field of rank k < 5, with the field's bounds: null when stage_state() is.
   class(stepper), intent(in)      :: self                 !< Stepper.
   real(R_P), pointer, intent(out) :: u(SELAH_FRONT_DIMS)  !< Stage state.
   integer(I4P)                    :: slot                 !< Register slot.
   integer(I4P)                    :: l(5), h(5)           !< Bounds of the field.

   u => null()
   l = self%lbk ; h = self%ubk
   if (self%rank /= SELAH_FRONT_RANK) error stop 'selah stepper: the view does not have the rank of the field'
   slot = state_slot(self)
   if (slot > 0_I4P) u(SELAH_VIEW_BOUNDS) => self%regs(:,:,:,:,:,slot)
   endsubroutine SELAH_STATE_VIEW

   subroutine SELAH_RHS_VIEW(self, r)
   !< The residual register of a field of rank k < 5, with the field's bounds: null outside a stage window.
   class(stepper), intent(in)      :: self                 !< Stepper.
   real(R_P), pointer, intent(out) :: r(SELAH_FRONT_DIMS)  !< Residual register.
   integer(I4P)                    :: slot                 !< Register slot.
   integer(I4P)                    :: l(5), h(5)           !< Bounds of the field.

   r => null()
   l = self%lbk ; h = self%ubk
   if (self%rank /= SELAH_FRONT_RANK) error stop 'selah stepper: the view does not have the rank of the field'
   slot = rhs_slot(self)
   if (slot > 0_I4P) r(SELAH_VIEW_BOUNDS) => self%regs(:,:,:,:,:,slot)
   endsubroutine SELAH_RHS_VIEW

   subroutine SELAH_LOOKAHEAD(self, k, bx, dt, q, u, mask, mask_comp)
   !< lookahead for a field of rank k < 5.
   class(stepper), intent(in)           :: self                    !< Stepper.
   integer(I4P),   intent(in)           :: k                       !< Stage.
   type(box),      intent(in)           :: bx                      !< Box of the field's rank.
   real(R_P),      intent(in)           :: dt                      !< Time step.
   real(R_P),      intent(in)           :: q(SELAH_FRONT_DIMS)     !< Solution.
   real(R_P),      intent(inout)        :: u(SELAH_FRONT_DIMS)     !< Look-ahead state (q's shape).
   real(R_P),      intent(in), optional :: mask(SELAH_FRONT_DIMS)  !< Mask (q's rank, its own components).
   integer(I4P),   intent(in), optional :: mask_comp               !< Mask component.
   integer(I4P)                         :: ubq(5), ubm(5)          !< Upper bounds of q and of the mask, as rank 5.

   call check_front(self, bx, SELAH_FRONT_RANK, shape(q, kind=I4P), ubq)
   if (any(shape(u) /= shape(q))) error stop 'selah stepper: the look-ahead state must have the shape of q'
   ubm = ubq
   if (present(mask)) ubm = front_ub(bx, shape(mask, kind=I4P))
   call lookahead_x(self, k, bx, dt, ubq, q, u, ubm, mask, mask_comp)
   endsubroutine SELAH_LOOKAHEAD
