#include "tde_lcurve_assistant.hpp"
#include "tde_lcurve.hpp"

int main(int argc, char* argv[])
{
    try
    {
        //check the options
        if(argc == 1)
        {
            std::cout << "Not used any flags. Use -help option!" << std::endl;
            exit(1);
        }

        //read parameter file and check the flags
        read_parameter_file();
        check_arguments(argc, argv);
        
        //Check that use any flags which create an output file
        if(!Parameters_run.lcurve && !Parameters_run.spectra && !Parameters_run.bolometric && !Parameters_run.diffusion)
        {
            std::cout << "Not used any flag which creat output file. Use -help option!" << std::endl;
            exit(1);
        }
        
        //Init Parameters_run struct
        Parameters_run.init(0);

        //Make simple monochromatic lightcurve
        if(Parameters_run.lcurve)
        {
            Total_lightcurve total_lc(Parameters_run);
            Parameters_run.print_parameters();
            total_lc.calculate_lightcurve_timeintervall(Parameters_run.t_start, Parameters_run.t_end, Parameters_run.dt);
        }

        //Make a spectra
        if(Parameters_run.spectra)
        {
            Spectrum spectra(Parameters_run);
            spectra.calculate_spectrum_nuintervall(Parameters_run.nu_start, Parameters_run.nu_end, Parameters_run.dnu);
        }

        //Make bolometric lightcurve
        if(Parameters_run.bolometric)
        {
            Lbol Bolometric_lcurve(Parameters_run);
            Bolometric_lcurve.Lbol_time(Parameters_run.t_start, Parameters_run.t_end, Parameters_run.dt);
        }

        //Make monochromatic lightcurve with diffusion
        if(Parameters_run.diffusion)
        {
            Diffusion_luminosity_timeintervall diff_lum(Parameters_run);
            diff_lum.calc_diffusion_luminosity_timeintervall(Parameters_run.tmin + 0.1, Parameters_run.t_end, Parameters_run.dt);
        }
    }

    //catch error massages
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    return 0;
}