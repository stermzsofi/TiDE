#ifndef TDE_ASSISTANT
#define TDE_ASSISTANT

#include "tde_lcurve.hpp"
//#include "tde_parameters.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <cmath>

//Declaration of Parameters struct and definition of Parameters_run: itt will contain the parameter settings of the ligt-curve
struct Parameters;
extern Parameters Parameters_run;

/*used constants*/

const double c = 3.e8;                          //light velocity in m/s
const double G = 6.673e-11;                     //Gravitational constant in m^3/(s^2*kg)
const double M_SUN = 1.9891*1e30;               //Sun's mass in kg
const double METER_AU = 1./149597870000.;       //Meter in AU
const double AU_METER = 1.5e11;                 //AU in m
const double h = 6.63e-34;                      //Planck constant in Js
const double k = 1.38e-23;                      //Boltzmann constant in J/K
const double PI = 3.14159;                      //Value of pi
const double hPERk = h/k;
const double hPERc2 = h/(c*c);
const double SIGMA = 5.67e-8;                    //Sigma in J*/(m^2*s*K^4)
const double PLANCK_CONST = 2.0*hPERc2;
const double a_const = 7.566 * 1e-16;            //a is the constant in the equation of radiationenergy, SI units
const double R_SUN_METER = 696340000;            //Radius of sun in m
const double YEAR_S = 31556926;
const double KAPPA_ES = 0.03977264;              //opacity of electron scattering in m^2/kg
const double s_to_day = 1./86400.;
const double YEAR_TO_DAY = 365.242199;
const double RHO_C_N3_CONSTPART = 12.9351;               //The constant part of rho_c if n =3
const double K_N3_CONSTPART = 4.40728e-12;          //The constant part of K if n=3
const double RHO_C_N3PER2_CONSTPART = 1.43017;      //The constant part of rho_c if n =3/2
const double K_N3PER2_CONSTPART = 2.51254e-11;      //The constant part of K if n=3/2
const double tpeak_relative_to_tmin_n3 = 5.77316;
const double tpeak_relative_to_tmin_n3per2 = 2.65911;
const double overflow_corrector = 1e-300;




enum startype : int {ms, wd, no_relationship};                              //startype enum: ms - main sequence, wd - white dwarf, no_relationship: no relationship - mass and radius of the star will independent

void chech_option(std::string name, std::string value, bool argument);      //check the options and set the needed parameter

void check_help(int argc, char* argv[]);

void read_parameter_file();                                                 //read parameters.txt file and set the parameters
void check_arguments(int argc, char* argv[]);                               //This function will check the options, set the parameters
double Planck_function_for_nu(double nu, double T);                         //calculate the Planck function to a specified frequency (nu [Hz], T [K])
double frequency_to_wavelength(double f);                                   //exchange frequency (Hz) to wavelength (nm)
double wavelength_to_frequency(double w);                                   //exchange wavelength (nm) to frequency (Hz)

void print_help();                                                          //print the help of the TDE program

bool IsDouble(const std::string& s);                                        //decide on a string that it is a double or not

#endif