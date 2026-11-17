!< SELAH verification: the report mode shared by the verification programs.
module selah_vv_cases
!< SELAH verification: the report mode shared by the verification programs.
!<
!< Every verification program asserts its bounds and exits non-zero on a violation, so that it runs in the test suite
!< (scripts/run_tests.sh). With the environment variable SELAH_VV_DIR set to a directory, the same program also writes
!< its data there (tables, curves, frames) for the documentation's tables, figures and animations (scripts/run_vv.sh,
!< docs/scripts/make_doc_figures.py): the report mode. The asserted bounds are the same in both modes.
implicit none
private
public :: vv_dir, vv_report, vv_open

contains
   function vv_dir() result(d)
   !< The report directory (SELAH_VV_DIR), '' outside the report mode.
   character(:), allocatable :: d
   character(512)            :: buf
   integer                   :: l, st

   call get_environment_variable('SELAH_VV_DIR', buf, l, st)
   d = ''
   if (st == 0 .and. l > 0) d = trim(buf)
   endfunction vv_dir

   function vv_report() result(r)
   !< True in the report mode.
   logical :: r

   r = len(vv_dir()) > 0
   endfunction vv_report

   subroutine vv_open(name, u, header)
   !< Open <SELAH_VV_DIR>/<name> for writing, with a '#' header line.
   character(*), intent(in)  :: name, header
   integer,      intent(out) :: u

   call execute_command_line('mkdir -p '//vv_dir())
   open(newunit=u, file=vv_dir()//'/'//name, status='replace', action='write')
   write(u, '(A)') '# '//header
   endsubroutine vv_open
endmodule selah_vv_cases
