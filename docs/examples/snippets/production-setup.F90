! the caller's pool: the field's bounds, the slots numbered from 1, as many as the form needs
allocate(pool(1:1, 0:NI+1, 0:2, 0:2, 1:NB, 1:registers_needed(SCHEME_SSP_54, FORM_BUTCHER)))
call st%initialize(SCHEME_SSP_54, FORM_BUTCHER, lb=lbound(q), ub=ubound(q), registers=pool, stages_on_q=.true.)
