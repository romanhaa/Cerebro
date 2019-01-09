
# ps 1.1.0

* New `ps_num_fds()` returns the number of open files/handles.

* New `ps_open_files()` lists all open files of a process.

* New `ps_interrupt()` interrupts a process. It sends a `SIGINT` signal on
  POSIX systems, and it can send CTRL+C or CTRL+BREAK events on Windows.

* New `ps_users()` lists users connected to the system.

* New `ps_mark_tree()`, `ps_find_tree()`, `ps_kill_tree()`,
  `with_process_cleanup()`: functions to mark and clean up child
  processes.

* New `CleanupReporter`, to be used with testthat: it checks for
  leftover child processes and open files in `test_that()` blocks.

# ps 1.0.0

First released version.
