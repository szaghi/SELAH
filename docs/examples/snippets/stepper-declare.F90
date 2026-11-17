type(stepper)            :: st
type(box)                :: bx
real(R_P), allocatable   :: q(:,:,:,:,:)              ! the solution: the caller's array
real(R_P), pointer       :: u(:,:,:,:,:), r(:,:,:,:,:)
