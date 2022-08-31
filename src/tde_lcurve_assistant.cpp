#include "tde_lcurve_assistant.hpp"

Parameters Parameters_run;

void chech_option(std::string name, std::string value, bool argument)
{
    if(IsDouble(value))
    {
        if("tstart" == name)
        {
            Parameters_run.t_start = stod(value);
            Parameters_run.t_start_default = false;
        }
        else if("tend" == name)
        {
            Parameters_run.t_end = stod(value);
            Parameters_run.t_end_relative_to_tmin = false;
        }
        else if("tend_rt_tmin" == name)
        {
            Parameters_run.t_end = stod(value);
            Parameters_run.t_end_relative_to_tmin = true;
        }
        else if("dt" == name)
        {
            Parameters_run.dt = stod(value);
        }
        else if("M6" == name)
        {
            Parameters_run.bh.set_mass(stod(value));
        }
        else if("mstar" == name || "Mstar" == name)
        {
            Parameters_run.s.set_mass(stod(value));
        }
        else if("rstar" == name || "Rstar" == name)
        {
            Parameters_run.s.set_startype(value);
            Parameters_run.s.set_radius(stod(value));
        }
        else if("eta" == name)
        {
            Parameters_run.eta = stod(value);
        }
        else if("beta" == name)
        {
            Parameters_run.beta = stod(value);
        }
        else if("fout" == name)
        {
            Parameters_run.fout = stod(value);
        }
        else if("fv" == name)
        {
            Parameters_run.fv = stod(value);
        }
        else if("d" == name)
        {
            Parameters_run.d = stod(value);
        }
        else if("i" == name)
        {
            Parameters_run.i = stod(value);
        }
        else if("nu" == name)
        {
            Parameters_run.nu = stod(value);
        }
        else if("tdiff" == name)
        {
            Parameters_run.diffusion_timescale = stod(value);
        }
        else if("N" == name)
        {
            Parameters_run.N = stoi(value);
        }
        else if("nustart" == name)
        {
            Parameters_run.nu_start = stod(value);
        }
        else if("nuend" == name)
        {
            Parameters_run.nu_end = stod(value);
        }
        else if("dnu" == name)
        {
            Parameters_run.dnu = stod(value);
        }
        else if("time" == name)
        {
            Parameters_run.time = stod(value);
        }
        else if("beta_limit" == name)
        {
            Parameters_run.beta_fulldisrupt = stod(value);
        }
        else
        {
            if(argument)
            {
                throw std::runtime_error("ERROR! Unknown argument: -"+name);
            }
            else
            {
                throw std::runtime_error("ERROR! Error with parameter file line: "+name +" " +value);
            }
        }
    }
    else
    {
        if("rstar" == name)
        {
            Parameters_run.s.set_startype(value);
        }
        else if("politrop" == name)
        {
                Parameters_run.s.set_politrop(value);
        }
        else
        {
            if(argument)
            {
                throw std::runtime_error("ERROR! Error with argument: -"+name +" " +value);
            }
            else
            {
                throw std::runtime_error("ERROR! Error with parameter file line: "+name +value);
            }
        }
    }
}

void read_parameter_file()
{
    /*Open parameters.txt file*/
    std::ifstream input;
    input.open("parameters.txt");
    /*lines: it is a string vector. It will contain the lines of parameters.txt file*/
    std::vector<std::string> lines;
    /*Put the parameters.txt file lines to lines, string type vector*/
    if(input.is_open())
    {
        /*line: actual line of the file*/
        std::string line = "START";
        while(!input.eof())
        {
            std::getline(input, line);
            if(input.eof()) break;
            lines.push_back(line);
        }
    }
    else
    {
        throw std::runtime_error("ERROR! Couldn't open parameters.txt file\n");
    }

    /*One line contain two things: first is the name and second is the value*/
    std::string name;
    std::string value;
    /*Value should be a number most cases*/
    for(auto actual_line : lines)
    {
        std::istringstream act_line_stream(actual_line);
        std::getline(act_line_stream, name, '\t');
        std::getline(act_line_stream, value, '\t');
        chech_option(name, value, false);
    }
}

void check_help(int argc, char* argv[])
{
    std::string name;
    for(auto i = 1; i < argc; i++)
    {
        name = argv[i];
        if("-help" == name || "-h" == name)
        {
            print_help();
        }
    }
}

void check_arguments(int argc, char* argv[])
{
    std::string name, value;
    for(auto i = 1; i < argc; i++)
    {
        name = argv[i];
        /*set output type*/
        if("-b" == name || "-bolometric" == name)
        {
            Parameters_run.bolometric = true;
        }
        else if("-lc" == name)
        {
            Parameters_run.lcurve = true;
        }
        else if("-e" == name || "-extra" == name)
        {
            Parameters_run.extras = true;
        }
        else if("-s" == name || "-spectra" == name)
        {
            Parameters_run.spectra = true;
        }
        else if("-diffusion" == name)
        {
            Parameters_run.diffusion = true;
        }
        /*time dependent fout*/
        else if("-bi_fout" == name)
        {
            Parameters_run.change_foutcalculation(new Fout_Lodato_Rossi_eq28(Parameters_run));
            Parameters_run.is_bifout = true;
        }
        /*set type of Mdotfb*/
        else if("-Mdotfb_classic" == name)
        {
            Parameters_run.change_Mdotfbcalculation(new Mdotfb_classical_powerfunction(Parameters_run));
            if(i < argc - 1)
            {
                value = argv[i+1];
                if(IsDouble(value))
                {   
                    i++;
                    //Parameters_run.change_Mdotfbcalculation(new Mdotfb_L09_constrho_x(Parameters_run, stod(value)));
                    Parameters_run.calcMdotfb -> set_calc_ninf(value);
                }
                else if(value.front() != '-')
                {
                    i++;
                    Parameters_run.calcMdotfb -> set_calc_ninf(value);
                }
            }    
        }
        else if("-Mdotfb_L09_const" == name)
        {
            Parameters_run.change_Mdotfbcalculation(new Mdotfb_L09_constrho_x(Parameters_run));
            if(i < argc - 1)
            {
                value = argv[i+1];
                if(IsDouble(value))
                {   
                    i++;
                    //Parameters_run.change_Mdotfbcalculation(new Mdotfb_L09_constrho_x(Parameters_run, stod(value)));
                    Parameters_run.calcMdotfb -> set_calc_ninf(value);
                }
                else if(value.front() != '-')
                {
                    i++;
                    Parameters_run.calcMdotfb -> set_calc_ninf(value);
                }
            }    
        }
        else if("-Mdotfb_L09" == name)
        {
            Parameters_run.change_Mdotfbcalculation(new Mdotfb_L09_all(Parameters_run));
            Parameters_run.change_mdotpeakcalculation(new Mdotpeak_L09_all(Parameters_run));
        }
        /*set Mdotpeak type*/
        else if("-Mdotp_classic" == name)
        {
            Parameters_run.change_mdotpeakcalculation(new Mdotpeak_LodatoRossi(Parameters_run));
        }
        else if("-Mdotp_GR13" == name)
        {
            Parameters_run.change_mdotpeakcalculation(new Mdotpeak_GuillochonRamirez(Parameters_run));
        }
        else if("-Mdotp_L09_const" == name)
        {
            Parameters_run.change_mdotpeakcalculation(new Mdotpeak_L09_constrho(Parameters_run));
        }
        /*Set tmin type*/
        else if("-tmin_classic" == name)
        {
            Parameters_run.change_tmincalculation(new tmin_with_rt(Parameters_run));
        }
        else if("-tmin_with_rp" == name)
        {
            Parameters_run.change_tmincalculation(new tmin_Lodato(Parameters_run));
        }
        else if("-tmin_from_GR13" == name)
        {
            Parameters_run.change_tmincalculation(new tmin_from_GR13(Parameters_run));
        }
        /*Set a model family*/
        else if("-GR13" == name)
        {
            Parameters_run.change_tmincalculation(new tmin_from_GR13(Parameters_run));
            Parameters_run.change_mdotpeakcalculation(new Mdotpeak_GuillochonRamirez(Parameters_run));
            Parameters_run.change_Mdotfbcalculation(new Mdotfb_L09_constrho_x(Parameters_run));
            Parameters_run.calcMdotfb -> set_calc_ninf("GR13");
        }

        
        
        else if("-use_rph_limit" == name)
        {
            Parameters_run.use_rph_limit = true;
        }
        else
        {
            if(name.front() == '-')
            {
                name.erase(0,1);
                if(i != argc - 1)
                {
                    i++;
                    value = argv[i];
                    chech_option(name, value, true);
                }
                else
                {
                    throw std::runtime_error("ERROR! Error with argument: -"+name);
                }
            }
            else
            {
                throw std::runtime_error("ERROR! Unknown argument: "+name);
            }
        }
    }
}   

double Planck_function_for_nu(double nu, double T)
{
    double pl;
    pl = PLANCK_CONST*nu*nu*nu/(std::exp(hPERk*nu/T) -1);
    return pl;
}

//frequency in Hz, wavelength in nm
double frequency_to_wavelength(double f)
{
    double wavelength;
    wavelength = c/f;
    return wavelength*1e9;  //m to nm
}

//wavelength in nm, frequency in Hz
double wavelength_to_frequency(double w)
{
    double frequency;
    frequency = c/(w*1e-9);     //nm to m
    return frequency;
}

void print_help()
{
    std::cout << "---------------------TiDE code help----------------------" << std::endl << std::endl;

    std::cout << "TiDE program is for calculate TDE event light curves and spectrum." << std::endl << std::endl;

    std::cout << "Flags to create output file:" << std::endl;
    std::cout << "     -lc:                           Create lightcurve output file (lcurve.dat)" << std::endl;
    std::cout << "     -diffusion:                    Create lightcurve with diffusion (lcurve_diffusion.dat)" << std::endl;
    std::cout << "                                    WARNING: some cases during this process will throw a warning massage: 'WARNING: S is nan!'" << std::endl;
    std::cout << "                                    The process will create the full output file in this case too. For the reasons see the documentation." << std::endl;
    std::cout << "     -spectra (-s):                 Create spectra" << std::endl;
    std::cout << "     -bolometric (-b):              Create bolometric lightcurve" << std::endl;
    std::cout << "     -extra (-e):                   Use it together with lc flag. Will create an extras.dat file" << std::endl << std::endl;
    
    std::cout << "fout calculation methods:" << std::endl;
    std::cout << "     Constant fout value: this is the default case. No flag needed" << std::endl;
    std::cout << "     -bi_fout                       Time dependent (built-in) fout parameter" << std::endl << std::endl;

    std::cout << "Usable accretion rate models:" << std::endl;
    std::cout << "     -Mdotfb_classic                Set Mdotfb calculation classic power function (default)" << std::endl;
    std::cout << "     -Mdotfb_L09_const              Set Mdotfb calculation like in L09 paper with constant rho(R) case" << std::endl;
    std::cout << "     -Mdotfb_L09                    Set Mdotfb calculation like in L09 paper" << std::endl << std::endl;

    std::cout << "Usable tmin calculation models:" << std::endl;
    std::cout << "     -tmin_classic                  Set tmin calculation to the classic case with rt (default)" << std::endl;
    std::cout << "              For corrected tmin: set beta_d value with beta_limit flag" << std::endl;
    std::cout << "     -tmin_rp                       Set tmin calculation to old equation with rp" << std::endl;
    std::cout << "     -tmin_from_GR13                Set tmin from GR13 tpeak results. WARNING: use only if accretion rate model is L09" << std::endl << std::endl;

    std::cout << "Parameters file and settable parameters" << std::endl;
    std::cout << "parameters.txt file" << std::endl;
    std::cout << "     At parameters.txt use this sintax at all lines:      name value" << std::endl;
    std::cout << "     Separate these with tabulator" << std::endl;
    std::cout << "     Or: use -name value at flags" << std::endl;
    std::cout << "WARNING!  This file must exist in the same directory where you run TiDE." << std::endl;
    std::cout << "WARNING: Flags will overwrite the file settings!" << std::endl;
    std::cout << "Available settable parameters" << std::endl;
    std::cout << "     eta:             radiative efficiency of the black hole" << std::endl;
    std::cout << "     beta:            penetration factor" << std::endl;
    std::cout << "     fout:            factional mass of the debris that is supposed to leave the system via super-Eddington wind" << std::endl;
    std::cout << "     fv:              ratio between the velocity of the super-Eddington wind and the escape velocity at rL" << std::endl;
    std::cout << "     d:               distance of the event in meters" << std::endl;
    std::cout << "     i:               inclination of the disk" << std::endl;
    std::cout << "     tdiff:           photon diffusion timescale in days" << std::endl;
    std::cout << "     N:               Number of concentric rings used to calculate the accretion disk during the numerical integration" << std::endl;
    std::cout << "     M6:              Mass of the black hole in 10^6 M_sun units" << std::endl;
    std::cout << "     mstar (Mstar):   Mass of the star in M_sun units" << std::endl;
    std::cout << "     rstar (Rstar):   Radius of the star. you can use the name of the relation (MS, WD) or a number" << std::endl;
    std::cout << "     politrop:        default case: use 5/3 politrop for wd stars and 4/3 all other case" << std::endl;
    std::cout << "                      to use other: print 5per3 or 4per3 after switch" << std::endl;
    std::cout << "     tstart:          the moment of the beginning of the light curve" << std::endl;
    std::cout << "     tend:            the end time of the light curve " << std::endl;
    std::cout << "     tend_rt_tmin:    set tend parameter relative to tmin (tend = tmin + add value)" << std::endl;
    std::cout << "     dt:              time step of the light curve" << std::endl;
    std::cout << "     nu:              the frequency for computing the monochromatic light curve " << std::endl;
    std::cout << "     nustart:         the start of the frequency interval " << std::endl;
    std::cout << "     nuend:           the end of the frequency interval " << std::endl;
    std::cout << "     dnu:             the frequency step of the spectra " << std::endl;
    std::cout << "     time:            the epoch for the computed spectrum" << std::endl;
    std::cout << "     beta_limit:      parameter for corrected tmin calculation" << std::endl << std::endl;

    std::cout << "For more information: see the documentation of TiDE" << std::endl;
    exit(1);
}


bool IsDouble(const std::string& s)
{
    std::regex reg("^[+-]?[[:digit:]]+\\.?[[:digit:]]*[eE][+-]?[[:digit:]]+$|^[+-]?[[:digit:]]+\\.?[[:digit:]]*$");
    return std::regex_match(s,reg);
}



