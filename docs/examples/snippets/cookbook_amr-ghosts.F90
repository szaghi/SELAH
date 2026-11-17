if (right) then
   call fill_ghosts(u)                            ! from the neighbours' stage state: right
else
   call fill_ghosts_from(u, q)                    ! from the neighbours' q^n: a stale exchange
endif
