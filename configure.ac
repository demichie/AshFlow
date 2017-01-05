# initial information about the project
AC_INIT([AshFlow], [1.0], [sameng@bgs.ac.uk])

# check if the source folder is available
AC_CONFIG_SRCDIR([src/flow_model.f90])

# check for fortran compiler
AC_PROG_FC

# automake initialization (mandatory) including a check for automake API version >= 1.9
AM_INIT_AUTOMAKE([1.9])

# files to generate via autotools (prepare .am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])


AC_PREFIX_DEFAULT([..])

# finally this generates the Makefiles etc. for the build
AC_OUTPUT
