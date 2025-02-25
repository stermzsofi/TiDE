#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([TiDE], [2.0], [stermzsofi@freemail.hu])
AC_CONFIG_SRCDIR([src/tde_lcurve_assistant.cpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

AC_LANG(C++)

# Checks for programs.

AC_PROG_CXX

# Checks for libraries.
#AC_CHECK_HEADER(boost/math/special_functions/bessel.hpp,[AC_DEFINE([HAVE_BOOST],[],["Bessel function will use from boost"])
#ac_cv_have_boost=yes],[ ac_cv_have_boost=no ])
# Checks for header files.
#AC_CHECK_HEADER(tr1/cmath,[AC_DEFINE([HAVE_TR1],[],["Bessel function will use from tr1"])
#	]ac_cv_have_tr,[])
AX_CXX_HEADER_STDCXX_TR1

BOOST_REQUIRE([1.23.0])
#BOOST_MATH
AX_CXX_BOOST_CHECK_BESSEL
# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AX_CXX_CHECK_BESSEL

#Check c++17 standard
AX_CXX_COMPILE_STDCXX(17,[],mandatory)

if test "$ax_cv_cxx_stdcxx_tr1" = no && test "$ax_cv_cxx_boost_bessel" = no && test "$ax_cv_cxx_bessel" = no ;then
	AC_MSG_ERROR([Your compiler misses ISO 29124:2010 and TR1. Please install a compiler that supports at least TR1 C++ extension or install boost C++ library.])
fi

#initialize libtool
AM_PROG_AR
LT_INIT([disable-static])
LT_LANG([C++])


AC_CONFIG_FILES([Makefile src/Makefile])

#automake
AM_INIT_AUTOMAKE([foreign -Wall subdir-objects])
AC_OUTPUT
