#include "tide_wrapper.h"

/*
    Parameters struct
*/

Parameters* Cparameters()
{
    return new Parameters();
}

void param_init(Parameters* par)
{
    par->init();
}

void param_refreshtime(Parameters* par, double time)
{
    par->refresh(time);
}

void param_refresh(Parameters* par)
{
    par->s.init();
    par->bh.init();
    par->calculate_rt();
    par->calctmin->calc_tmin();
    //std::cout << "tim calculated as " << par->tmin << std::endl;
    par->refresh(par->tmin);
}

void param_set_allparams(Parameters* par, double bh_M6, double star_mstar, 
    char* star_rstar, char* star_politrop, double eta, double eta_r, char* fout, double fv,
double d, double i, char* diff_timescale, double beta_fulldisrupt, int N, 
char* Mdotpeakcalc, char* Mdotfb_calc, char* tmin_calc,
double nu)
{
    par->bh.set_mass(bh_M6);
    par->s.set_mass(star_mstar);
    par->s.set_startype(std::string(star_rstar));
    std::string mypolytrop(star_politrop);
    if(mypolytrop != "default")
    {
        par->s.set_politrop(mypolytrop);
    }
    par->eta = eta;
    par->eta_reprocessing = eta_r;
    std::string myfout(fout);
    if(IsDouble(myfout))
    {
        par->fout = stod(myfout);
    }
    else
    {
        par->change_foutcalculation(new Fout_Lodato_Rossi_eq28(*par));
        par->is_bifout = true;
        //std::cout << "fout done" << std::endl;
    }
    par->fv = fv;
    par->d = d;
    par->i = i;
    std::string my_diff_timescale(diff_timescale);
    if(IsDouble(my_diff_timescale))
    {
        par->diffusion_timescale = stod(my_diff_timescale);
    }
    par->beta_fulldisrupt = beta_fulldisrupt;
    par->N = N;
    std::string my_mdotpeakcalc(Mdotpeakcalc);
    if(my_mdotpeakcalc == "L09_const")
    {
        par->change_mdotpeakcalculation(new Mdotpeak_L09_constrho(*par));
    }
    else if(my_mdotpeakcalc == "GR13")
    {
        par->change_mdotpeakcalculation(new Mdotpeak_GuillochonRamirez(*par));
    }
    else if(my_mdotpeakcalc == "classical")
    {
        par->change_mdotpeakcalculation(new Mdotpeak_LodatoRossi(*par));
    }
    std::string my_mdotfbcalc(Mdotfb_calc);
    if(my_mdotfbcalc == "default"){}
    else if(my_mdotfbcalc == "L09")
    {
        par->change_Mdotfbcalculation(new Mdotfb_L09_all(*par));
        par->change_mdotpeakcalculation(new Mdotpeak_L09_all(*par));
    }
    else if(my_mdotfbcalc == "L09_const")
    {
        par->change_Mdotfbcalculation(new Mdotfb_L09_constrho_x(*par));
    }
    std::string my_tmincalc(tmin_calc);
    if(my_tmincalc == "with_rp")
    {
        par->change_tmincalculation(new tmin_Lodato(*par));
    }
    else if(my_tmincalc == "GR13")
    {
        par->change_tmincalculation(new tmin_from_GR13(*par));
    }
    par->nu = nu;
}

void param_set_params(Parameters* par, double bh_M6, double star_mstar, 
    double eta, double eta_r,
    double fv, double d, double i, char* diff_timescale, 
    double beta_fulldisrupt, int N, 
    double nu)
{
    par->bh.set_mass(bh_M6);
    par->s.set_mass(star_mstar);
    par->eta = eta;
    par->eta_reprocessing = eta_r;
    par->fv = fv;
    par->d = d;
    par->i = i;
    std::string my_diff_timescale(diff_timescale);
    if(IsDouble(my_diff_timescale))
    {
        par->diffusion_timescale = stod(my_diff_timescale);
    }
    par->beta_fulldisrupt = beta_fulldisrupt;
    par->N = N;
    par->nu = nu;
}

void print_params(Parameters* par)
{
    par->print_parameters();
}

void param_set_timedependent_pars(Parameters* par, double time)
{
    par->calculate_timedependent_parameters(time);
}

double get_parameters_time(Parameters* par)
{
    return par->time;
}

double get_parameters_tmin(Parameters* par)
{
    //std::cout << "Hali " << par->tmin << std::endl;
    return par->tmin;
}

/*
    Wind part of lightcurve
*/

Wind_part_of_lightcurve* Cwind(Parameters* par)
{
    return new Wind_part_of_lightcurve(*par);
}

double Lwind_at_nu_time(Wind_part_of_lightcurve* w, double time, double nu)
{
    if(w->par.time != time)
    {
        w->par.calculate_timedependent_parameters(time);
    }
    //w->par.nu = nu;
    w->par.calculate_timedependent_parameters(time);
    return w->calc_luminosity_at_nu(nu);
}

/*
    Disk part of lightcurve
*/

Disk_part_of_lightcurve* Cdisk(Parameters* par)
{
    return new Disk_part_of_lightcurve(*par);
}

double Ldisk_at_nu_time(Parameters* par, Disk_part_of_lightcurve* d, double time, double nu)
{
    if(par->time != time)
    {
        par->calculate_timedependent_parameters(time);
    }
    return d->calc_disk_part_at_t_nu(nu);
}

double CLdisk_at_time(Parameters* par, Disk_part_of_lightcurve* d, double time)
{
    if(par->time != time)
    {
        par->calculate_timedependent_parameters(time);
    }
    return d->calc_disk_part_at_t();
}

Diffusion_luminosity* Cdiffused(Parameters* par)
{
    return new Diffusion_luminosity(*par);
}

void Ldiffused_reset(Diffusion_luminosity* difflum)
{
    difflum->reset();
}

double Ldiffused_at_time(Diffusion_luminosity* difflum, Parameters* par, double time)
{
    par->calculate_timedependent_parameters(time);
    difflum->set_time(time);
    return difflum->calc_Lum_at_t();
}