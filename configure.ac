#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([TiDE], [2.0], [stermzsofi@freemail.hu])
AC_CONFIG_SRCDIR([src/tde_lcurve_assistant.cpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
#AC_CHECK_FUNCS([pow sqrt])
#Check c++17 standard
AX_CXX_COMPILE_STDCXX(17)

AC_CONFIG_FILES([Makefile src/Makefile])

#automake
AM_INIT_AUTOMAKE([foreign -Wall subdir-objects])
AC_OUTPUT
