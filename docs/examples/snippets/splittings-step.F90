call st%begin_step(bx, y)
do k=1, st%stages()
   call st%begin_stage(k, bx, dt, y)
   r => st%rhs_register()
   select case(st%stage_partition(k))                      ! the flow of stage k advances one partition:
   case(1) ; r(1,1,1,1,1) = y(2,1,1,1,1)                   !   q' = p  (its residual depends on p only)
   case(2) ; r(2,1,1,1,1) = -y(1,1,1,1,1)                  !   p' = -q (on q only)
   endselect
   call st%end_stage(k, bx, dt, y)                         ! the partition's components += c_k dt r
enddo
call st%finish_step(bx, dt, y)
