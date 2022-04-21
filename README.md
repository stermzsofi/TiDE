# TiDE

TiDE (TIdal Disruption Event) is a C++ code that computes the light curves or spectrum of tidal disruption events.

# Structure of TiDE and installation

The code is decomposed into four hpp and cpp file.
	• tde_assistant: contains the used physical constants, some unit conversions, parameter file and argument management, Planck function and a print help function of the code.
 
  • tde_lcurve: this is the principal file of the code. It contain the parameters, the two parts of the light curve routine (wind and disk part) and various classes for the computation of the light curves and spectra.
 
  • Trapezoidal_rule/trapezoidal: contains the integrator that uses the trapezoidal rule.
 
  • tde_lcurve_run: this is the only file which has not got a hpp file. This is a cpp file which contain the main function of the code.

To run the program one needs to compile all of these cpp files. For convenience, there is a Makefile that can be used to compile and link the different parts of the code together. IMPORTANT: ONE MUST USE AN AT LEAST C++17 COMPLIANT compiler. The code was tested with g++ version 9.4.0 compiler with -g -Wall -Wextra -std=c++17 flags, which meets this requirement (c++17 standard).

