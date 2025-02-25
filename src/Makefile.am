bin_PROGRAMS = TiDE
TiDE_SOURCES = tde_lcurve_run.cpp tde_lcurve_assistant.cpp tde_lcurve_assistant.hpp tde_lcurve.cpp tde_lcurve.hpp tde_parameters.cpp tde_parameters.hpp Linear_interpol/linear_interpol.cpp Linear_interpol/linear_interpol.hpp Trapezoidal_rule/trapezoidal.cpp Trapezoidal_rule/trapezoidal.hpp
TiDE_CXXFLAGS = -std=c++17 -Wall -Wextra -g


lib_LTLIBRARIES = libtide.la

libtide_la_SOURCES = tde_lcurve_assistant.cpp tde_lcurve_assistant.hpp tde_lcurve.cpp tde_lcurve.hpp tde_parameters.cpp tde_parameters.hpp Linear_interpol/linear_interpol.cpp Linear_interpol/linear_interpol.hpp Trapezoidal_rule/trapezoidal.cpp Trapezoidal_rule/trapezoidal.hpp  Python_wrapper/tide_wrapper.cpp
include_HEADERS = Python_wrapper/tide_wrapper.h

