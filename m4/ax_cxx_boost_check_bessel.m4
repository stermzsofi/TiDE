AC_DEFUN([AX_CXX_BOOST_CHECK_BESSEL],[
	AC_CACHE_CHECK(for special math functions in boost,
		ax_cv_cxx_boost_bessel,
		[AC_LANG_PUSH([C++])
		ac_cxx_save=$CXX
		switch="-std=c++17"
		CXX="$CXX $switch"
		AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
		#include <boost/math/special_functions/bessel.hpp>
		
		namespace boost_bessel_test
		{
			double test = boost::math::cyl_bessel_i(0,2);
		}
		]], [])],
		[ax_cv_cxx_boost_bessel=yes],[ax_cv_cxx_boost_bessel=no])
		CXX=$ac_cxx_save
		AC_LANG_POP([C++])
		])
	if test "$ax_cv_cxx_boost_bessel" = yes; then
    AC_DEFINE(STDCXX_BOOST_SPEC_MATH,,[Special math functions are present in boost. ])
  fi		
])
