//#include "string.h"
#include "../tde_lcurve_assistant.hpp"
#include "../tde_parameters.hpp"

extern "C"
{
    //Parameters struct functions
    Parameters* Cparameters();
    void param_init(Parameters* par);
    void param_refreshtime(Parameters* par, double time);
    void param_refresh(Parameters* par);
    void param_set_allparams(Parameters* par, double bh_M6, double star_mstar, 
        char* star_rstar, char* star_politrop, double eta, double eta_r, char* fout,
        double fv, double d, double i, char* diff_timescale, 
        double beta_fulldisrupt, int N, 
        char* Mdotpeakcalc, char* Mdotfb_calc, char* tmin_calc,
        double nu);
    void param_set_params(Parameters* par, double bh_M6, double star_mstar, 
        double eta, double eta_r,
        double fv, double d, double i, char* diff_timescale, 
        double beta_fulldisrupt, int N, 
        double nu);
    void print_params(Parameters* par);
    void param_set_timedependent_pars(Parameters* par, double time);
    double get_parameters_time(Parameters* par);
    double get_parameters_tmin(Parameters* par);

    //Wind part of lc functions
    Wind_part_of_lightcurve* Cwind(Parameters* par);
    double Lwind_at_nu_time(Wind_part_of_lightcurve* w, double time, double nu);

    //Disk part of lc functions
    Disk_part_of_lightcurve* Cdisk(Parameters* par);
    double Ldisk_at_nu_time(Parameters* par, Disk_part_of_lightcurve* d, double time, double nu);
    double CLdisk_at_time(Parameters* par, Disk_part_of_lightcurve* d, double time);

    //Diffused light curve
    Diffusion_luminosity* Cdiffused(Parameters* par);
    void Ldiffused_reset(Diffusion_luminosity* difflum);
    double Ldiffused_at_time(Diffusion_luminosity* difflum, Parameters* par, double time);

}