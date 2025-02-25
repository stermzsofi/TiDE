

from ctypes import cdll
import ctypes
import numpy as np
from matplotlib import pyplot as plt
import os

libpath=os.getenv('TIDE_LIB')
if libpath is not None:
    lib = cdll.LoadLibrary(libpath+'/libtide.so')
    print("DEBUG: TIDE_LIB path is:"+libpath)
else:
    try:
        lib = cdll.LoadLibrary('libtide.so')
        print("DEBUG: No TIDE_LIB defined, try ldconfig defaults")
    except:
        lib = cdll.LoadLibrary(__file__[:-8]+'/../../../../src/.libs/libtide.so')
        print("DEBUG: backup to build directory in src.")
lib.get_parameters_tmin.restype=ctypes.c_double
lib.Lwind_at_nu_time.restype=ctypes.c_double
lib.Ldisk_at_nu_time.restype=ctypes.c_double
lib.CLdisk_at_time.restype=ctypes.c_double
lib.Cparameters.restype=ctypes.c_void_p
lib.Ldiffused_at_time.restype=ctypes.c_double

class Parameters(object):
    """
        Parameter class of TDE.
        This class handles most of the simulation data.

        Settable parameters and ther default values:
        -bh_M6 = 1: mass of the black hole in 10^6 M_sun units [float]
        -star_mstar = 1: mass of the star in M_sun units [float]
        -star_rstar = ms: radius of the star in R_sun units. ms MS: main sequence; wd WD: white dwarf [float or string]
        -star_politrop = default: politrop parameter of the star. default: setted based on the type of star. Can be 4per3 or 5per3 [string]
        -bh_eta = 0.1: radiative efficiency of the black hole [float]
        -eta_r = 0.0: the efficiency of the reprocessing [float]
        -fout = bifout: fout value in [0;1] range. bifout - build in calculation. [float or string]
        -fv = 1: constant in the velocity of the wind. must be geq 1 [float]
        -d = 0: distance. if 0, 1/(4 pi d^2) =1 will use [float]
        -i = 0: inclination of the disk [float]
        -diff_timescale = default: diffusion timescale. default or number [float or string]
        -beta_fulldisrupt = 1.85 #the limit of full disruption [float]
        -N = 1000 #Number of concentric rings used to calculate the accretion disk during numerical integration [int]
        -Mdotpeakcalc = "default" #Mdotpeak calculation method: default(not set), classic, L09_const, GR13
        -Mdotfb_calc = "L09" #Mdotfb calculation method: default (classic), L09, L09_const
        -tmin_calc = "default" #tmin calculation method: default (classic), with_rp, GR13
        -tstart = "default" #default (tmin) or number: the started time of the light curve [float or string]
        -tend=100: the end of the light curve [float]
        -dt=0.5: time step during light curve [float]

        Usable functions:
        -print_params(): print informations about C Parameters struct
        -param_init(): initialize the C parameters struct
        -

    """
    def __init__(self):
        self.param = lib.Cparameters() #create a C++ parameter class
        
        self.bh_M6 = 1.0 #mass of the black hole is 10^6 M_sun units [float]
        self.star_mstar = 1.0 #mass of the star in M_sun units [float]
        self.star_rstar = "ms" # mass of the star in R_sun units [float or string] (ms MS wd WD)
        self.star_rstar_copy = ""
        self.star_politrop = "default" # politrop index of the star: must be string: 4per3 or 5per3
        self.star_politrop_copy = ""
        self.bh_eta = 0.1 #radiative efficiency of the black hole [float]
        self.eta_r = 0.0 #the efficiency of the reprocessing [float]
        self.fout = "bifout" #bifout - build in calculation or a number in range [0;1] [float]
        self.fout_copy = ""
        self.fv = 1 # number larger or equal than 1 [float]
        self.d = 0 # distance: if 0, 1/(4 pi d^2) = 1 will use [float]
        self.i = 0 # inclination of the disk [float]
        self.diff_timescale = "default" # diffusion timescale: default or number [float]
        self.beta_fulldisrupt = 1.85 #the limit of full disruption [float]
        self.N = 1000 #Number of concentric rings used to calculate the accretion disk during numerical integration [int]
        self.Mdotpeakcalc = "default" #Mdotpeak calculation method: default(not set), classic, L09_const, GR13
        self.Mdotpeakcalc_copy = ""
        self.Mdotfb_calc = "L09" #Mdotfb calculation method: default (classic), L09, L09_const
        self.Mdotfb_calc_copy = ""
        self.tmin_calc = "default" #tmin calculation method: default (classic), with_rp, GR13
        self.tmin_calc_copy = ""
        self.tstart = "default" #default or number
        self.tstart_value = 0
        self.tend = 100
        self.dt = 0.5
        self.nu = 6.3e14

        self.param_init()

    #def read_param(self):
    #    pass

    def refresh_all_params(self):
        """
            Refresh all parameters (include calculation methods) az C parameters struct based on Python setted parameters
        """
        lib.param_set_allparams(self.param,ctypes.c_double(self.bh_M6), 
                        ctypes.c_double(self.star_mstar), ctypes.c_char_p(str(self.star_rstar).encode('utf-8')),
                        ctypes.c_char_p(self.star_politrop.encode('utf-8')),ctypes.c_double(self.bh_eta), 
                        ctypes.c_double(self.eta_r), ctypes.c_char_p(self.fout.encode('utf-8')),
                        ctypes.c_double(self.fv), ctypes.c_double(self.d), ctypes.c_double(self.i),
                        ctypes.c_char_p(self.diff_timescale.encode('utf-8')), ctypes.c_double(self.beta_fulldisrupt), ctypes.c_int(self.N), 
                        ctypes.c_char_p(self.Mdotpeakcalc.encode('utf-8')), ctypes.c_char_p(self.Mdotfb_calc.encode('utf-8')), ctypes.c_char_p(self.tmin_calc.encode('utf-8')),
                        ctypes.c_double(self.nu))
        lib.param_init(self.param)

    def refresh_params(self):
        """
            Refresh only the parameters (without calculation methods)
            at C parameters struct based on Python setted parameters
        """
        lib.param_set_params(self.param,ctypes.c_double(self.bh_M6), ctypes.c_double(self.star_mstar),
                        ctypes.c_double(self.bh_eta), ctypes.c_double(self.eta_r),
                        ctypes.c_double(self.fv), ctypes.c_double(self.d), ctypes.c_double(self.i),
                        ctypes.c_char_p(self.diff_timescale.encode('utf-8')), ctypes.c_double(self.beta_fulldisrupt),
                        ctypes.c_int(self.N), ctypes.c_double(self.nu))
        lib.param_refresh(self.param)

    def print_params(self):
        """
            Print some important values about C parameters struct
        """
        lib.print_params(self.param)

    def param_init(self):
        """
            Initialize the C parameters struct based on Python setted parameters
        """
        if (self.star_rstar != self.star_rstar_copy or 
            self.star_politrop != self.star_politrop_copy or 
            self.fout != self.fout_copy or
            self.Mdotpeakcalc != self.Mdotpeakcalc_copy or
            self.Mdotfb_calc != self.Mdotfb_calc_copy or
            self.tmin_calc != self.tmin_calc_copy):
            self.refresh_all_params()
            self.star_rstar_copy = self.star_rstar
            self.star_politrop_copy = self.star_politrop
            self.fout_copy = self.fout
            self.Mdotpeakcalc_copy = self.Mdotpeakcalc
            self.Mdotfb_calc_copy = self.Mdotfb_calc
            self.tmin_calc_copy = self.tmin_calc 
        else:
            self.refresh_params()
        if (self.tstart == "default"):
            self.set_tstart()
        else:
            self.tstart_value = self.tstart

    def parameters_calc_timedep_parameters(self, time):
        """
            Calculate and set the time dependent parameters at C parameters struct
        """
        lib.param_set_timedependent_pars(self.param, ctypes.c_double(time))

    def get_parameters_time(self):
        return lib.get_parameters_time(self.param)
    
    def set_tstart(self):
        """
            Set tstart value equal as tmin
        """
        self.tstart_value = lib.get_parameters_tmin(self.param)
        #print(lib.get_parameters_tmin(self.param))

class Wind_part_of_lc(object):
    """
        Python wind part of the light curve class

        -Init function parameter: Python Parameters class
    """
    def __init__(self, param):
        """
            Parameters: 
            -param: Python Parameters class
        """
        self.pythonpar = param
        self.Cpar=param.param
        self.wind_obj = lib.Cwind(self.Cpar)
    
    def wind_luminosity_at_time_nu(self, time, nu):
        """
            Calculate the wind luminosity at a given time and nu (frequency)
        """
        lum = lib.Lwind_at_nu_time(self.wind_obj,ctypes.c_double(time),ctypes.c_double(nu))
        return lum

    def wind_luminosity_at_time(self, time):
        """
            Calculate the wind luminosity at a given time
        """
        lum = lib.Lwind_at_nu_time(self.wind_obj,ctypes.c_double(time),ctypes.c_double(self.pythonpar.nu))
        return lum
    
class Disk_part_of_lc(object):
    """
        Python disk part of the light curve class

        -Init function parameter: Python Parameters class
    """
    def __init__(self, param):
        """
            Parameters: 
            -param: Python Parameters class
        """
        self.pythonpar = param
        self.Cpar=param.param
        self.disk_obj = lib.Cdisk(self.Cpar)

    def disk_luminosity_at_time_nu(self, time, nu):
        """
            Calculate the disk luminosity at a given time and nu (frequency)
        """
        lum = lib.Ldisk_at_nu_time(self.Cpar, self.disk_obj,ctypes.c_double(time),ctypes.c_double(nu))
        return lum

    def disk_luminosity_at_time(self, time):
        """
            Calculate the disk luminosity at a given time
        """
        lum = lib.CLdisk_at_time(self.Cpar, self.disk_obj,ctypes.c_double(time))
        return lum
    
class Diffused_luminosity(object):
    def __init__(self, param):
        self.Cpar = param.param
        self.diff_obj = lib.Cdiffused(self.Cpar)

    def calc_lum_at_time(self, time):
        lum = lib.Ldiffused_at_time(self.diff_obj,self.Cpar,ctypes.c_double(time))
        return lum
        
class Light_curve_of_tde(object):
    """
        Class to calculate the light curve of a TDE at a given time range.

        Functions:
            Init: 
            -param: Python Parameters class

            light_curve
            -Return arrays: times, sumlum, windlum, disklum
    """
    def __init__(self, param):
        self.Cpar = param.param
        self.PythonPar = param
        self.PythonDisk = Disk_part_of_lc(param)
        self.PythonWind = Wind_part_of_lc(param)
        self.PythonDiffusion = Diffused_luminosity(param)

    def light_curve(self):
        """
        Calculation of a light curve.

            Return arrays:
            -times: the calculated times
            -sumlum: the total calculated luminosity
            -windlum: the wind part of the luminosity
            -disklum: the disk part of the luminosity
        """
        #print(self.PythonPar.tstart, self.PythonPar.bh_M6)
        
        times=np.arange(self.PythonPar.tstart_value,self.PythonPar.tend,self.PythonPar.dt)
        windlumcalc = np.vectorize(self.PythonWind.wind_luminosity_at_time)
        windlum = windlumcalc(times)
        disklumcalc = np.vectorize(self.PythonDisk.disk_luminosity_at_time)
        disklum = disklumcalc(times)
        sumlum = windlum + disklum
        return times,sumlum, windlum, disklum
    


    def diffused_light_curve(self):
        """
        Calculate light curve with diffusion.
        Returns:
            times: time array
            difflum: Total luminosity with diffusion
        """
        lib.Ldiffused_reset(self.PythonDiffusion.diff_obj)
        times=np.arange(self.PythonPar.tstart_value,self.PythonPar.tend,self.PythonPar.dt)
        difflumcalc = np.vectorize(self.PythonDiffusion.calc_lum_at_time)
        difflum = difflumcalc(times)
        return times, difflum





def main():
    p = Parameters()
    #p.star_rstar = "ms"
    #p.Mdotfb_calc = "L09"
    #p.param_init()
    #print(p.tstart)
    #print("hello",float(lib.get_parameters_tmin(p.param)))
    #p.print_params()
    #print("ez:",type(p.param),p.param)
    #w = Wind_part_of_lc(p)
    #times=np.arange(p.tstart,p.tend,p.dt)
    #windcalc=np.vectorize(lambda x: w.wind_luminosity_at_time_nu(x,p.nu))
    #windlums=windcalc(times)
    #plt.plot(times,[w.wind_luminosity_at_time_nu(x,p.nu) for x in times])
    #plt.yscale('log')
    #plt.show()
    lc = Light_curve_of_tde(p)
    res = lc.light_curve()
    #plt.plot(res[0],res[1]*p.nu, label="sum")
    #plt.plot(res[0],res[2]*p.nu,label="wind")
    #plt.plot(res[0],res[3]*p.nu,label="disk")
    #print(lc.PythonPar.tstart_value)

    p.bh_M6 = 5.0
    #print("ez:",type(p.param),p.param)
    #p.param_init()
    #p.param_refresh()
    p.param_init()
    #p.print_params()
    #print(lc.PythonPar.tstart_value)
    #lc.PythonPar.print_params()
    #lc.PythonWind.pythonpar.print_params()
    res2 = lc.light_curve()
    #plt.plot(res2[0],res2[1]*p.nu, label="sum2")

    p.star_rstar = "wd"
    p.param_init()
    #p.print_params()
    res3 = lc.light_curve()
    #plt.plot(res3[0],res3[1]*p.nu, label="sum3")

    plt.legend()
    plt.yscale('log')
    plt.ylim((1e36,None))
    #plt.show()
    plt.clf()
    massrange=[0.5,1,10,15]
    r = light_curve_in_mass_range(massrange)
    for i in range(0,4):
        plt.plot(r[i][0],r[i][1]*p.nu,label=f"{i}")
    plt.legend()
    
    plt.ylim(bottom=1e36)
    plt.xscale('log')
    plt.yscale('log')
    #plt.show()

    p1 = Parameters()
    p1.tend=1000.0
    
    lc1  = Light_curve_of_tde(p1)
    difflum = lc1.diffused_light_curve()
    plt.clf()
    plt.plot(difflum[0],difflum[1]*p.nu,label="0")

    p1.diff_timescale = "6.68"
    p1.param_init()
    difflum1 = lc1.diffused_light_curve()
    plt.plot(difflum1[0],difflum1[1]*p.nu,label="1")
    plt.show()


if __name__=='__main__':
    main()
        
    
