!< SELAH expected-failure test: an additive splitting with backward substeps, not declared reversible, stops.
program selah_test_xfail_additive_backward
!< SELAH expected-failure test: an additive splitting with backward substeps, not declared reversible, stops.
!<
!< yoshida-4 runs each operator backward over part of the step (every real splitting of order 3 or more does): ill-posed
!< for an irreversible field step or flow (an upwind field step, cooling). `initialize` must stop the program unless the
!< caller declares both operators reversible. Exit status non-zero = pass (scripts/run_tests.sh, *_xfail_*).
use selah
implicit none
type(additive_splitting) :: sp

call sp%initialize(SCHEME_YOSHIDA_4)
print '(A)', 'selah_test_xfail_additive_backward: yoshida-4 accepted without reversible=.true. (wrong)'
endprogram selah_test_xfail_additive_backward
