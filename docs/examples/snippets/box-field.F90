! V_FIRST layout (v,i,j,k,b): the components innermost (the CPU layout); ghosts 1-NGC..N+NGC on i, j, k
allocate(x(1:NV, 1-NGC:NI+NGC, 1-NGC:NJ+NGC, 1-NGC:NK+NGC, 1:NB))
allocate(y, mold=x)
