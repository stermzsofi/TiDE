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
        else
        {
            if(argument)
            {
                throw std::runtime_error("Unknown argument: -"+name);
            }
            else
            {
                throw std::runtime_error("Error with parameter file line: "+name +" " +value);
            }
        }
    }
    else
    {
        if("rstar" == name)
        {
            Parameters_run.s.set_startype(value);
        }
        else
        {
            if(argument)
            {
                throw std::runtime_error("Error with argument: -"+name +" " +value);
            }
            else
            {
                throw std::runtime_error("Error with parameter file line: "+name +value);
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
        throw std::runtime_error("Couldn't open parameters.txt file\n");
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

void check_arguments(int argc, char* argv[])
{
    std::string name, value;
    for(auto i = 1; i < argc; i++)
    {
        name = argv[i];
        if("-b" == name || "-bolometric" == name)
        {
            Parameters_run.bolometric = true;
        }
        else if("-e" == name || "-extra" == name)
        {
            Parameters_run.extras = true;
        }
        else if("-lc" == name)
        {
            Parameters_run.lcurve = true;
        }
        else if("-s" == name || "-spectra" == name)
        {
            Parameters_run.spectra = true;
        }
        else if("-bi_fout" == name)
        {
            Parameters_run.change_foutcalculation(new Fout_Lodato_Rossi_eq28(Parameters_run));
            Parameters_run.is_bifout = true;
        }
        else if("-Mdot_Guillochon" == name)
        {
            Parameters_run.change_mdotpeakcalculation(new Mdotpeak_GuillochonRamirez(Parameters_run));
            if(i < argc - 1)
            {
                value = argv[i+1];
            }
            if("5per3" == value || "4per3" == value)
            {
                i++;
                Parameters_run.calcmdotpeak -> set_politrop(argv[i]);
            }
            else
            {
                Parameters_run.calcmdotpeak  -> set_default_politrop();
            }
        }
        else if("-diffusion" == name)
        {
            Parameters_run.diffusion = true;
        }
        else if("-use_rph_limit" == name)
        {
            Parameters_run.use_rph_limit = true;
        }
        else if("-help" == name || "-h" == name)
        {
            print_help();
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
                    throw std::runtime_error("Error with argument: -"+name);
                }
            }
            else
            {
                throw std::runtime_error("Unknown argument: "+name);
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

    std::cout << "Available flags:" << std::endl;
    std::cout << "     -lc:                           Make lightcurve" << std::endl;
    std::cout << "     -diffusion:                    Make lightcurve used diffusion modell" << std::endl;
    std::cout << "     -spectra (-s):                 Make spectra" << std::endl;
    std::cout << "     -extra (-e):                   Make extras.dat file (during lightcurve case)" << std::endl;
    std::cout << "     -bolometric (-b):              Make bolometric lightcurve" << std::endl;
    std::cout << "     -bi_fout                       Use built-in fout, time dependent fout parameter" << std::endl;
    std::cout << "     -Mdot_Guillochon               Use Mdot_peak calculation from Guillochon paper" << std::endl;
    std::cout << "                                    default case: use 5/3 politrop index for wd stars and 4/3 all other case" << std::endl;
    std::cout << "                                    to use other: print 5per3 or 4per3 after Mdot_Guillochon switch" << std::endl << std::endl;
    
    std::cout << "Parameters file and settable parameters" << std::endl;
    std::cout << "parameters.txt file" << std::endl;
    std::cout << "     At parameters.txt use this sintax at all lines:      name value" << std::endl;
    std::cout << "     Separate these with tabulator" << std::endl;
    std::cout << "     Or: use -name value at flags" << std::endl;
    std::cout << "WARNING: Flags will overwrite the file settings!" << std::endl;
    std::cout << "Available settable parameters" << std::endl;
    std::cout << "     eta" << std::endl;
    std::cout << "     beta" << std::endl;
    std::cout << "     fout" << std::endl;
    std::cout << "     fv" << std::endl;
    std::cout << "     d" << std::endl;
    std::cout << "     i" << std::endl;
    std::cout << "     tdiff" << std::endl;
    std::cout << "     N" << std::endl;
    std::cout << "     M6" << std::endl;
    std::cout << "     mstar (Mstar)" << std::endl;
    std::cout << "     rstar (Rstar): you can use the name of the relation (MS, WD) or a number" << std::endl;
    std::cout << "     tstart" << std::endl;
    std::cout << "     tend" << std::endl;
    std::cout << "     tend_rt_tmin: set tend parameter relative to tmin (tend = tmin + add value)" << std::endl;
    std::cout << "     dt" << std::endl;
    std::cout << "     nu" << std::endl;
    std::cout << "     nustart" << std::endl;
    std::cout << "     nuend" << std::endl;
    std::cout << "     dnu" << std::endl;
    std::cout << "     time" << std::endl << std::endl;

    std::cout << "For more information: see the documentation of TiDE" << std::endl;
    exit(1);
}


bool IsDouble(const std::string& s)
{
    std::regex reg("^[+-]?[[:digit:]]+\\.?[[:digit:]]*[eE][+-]?[[:digit:]]+$|^[+-]?[[:digit:]]+\\.?[[:digit:]]*$");
    return std::regex_match(s,reg);
}



