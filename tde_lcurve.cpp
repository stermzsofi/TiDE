#include "tde_lcurve.hpp"

/*
    *******

        Objects derivation: Blackhole and Star

        Mass_to_radii derivated classes calc_radius functions
    
    *******
*/

/*Blackhole class functions*/

    /*mass: 10^6 M_sun units
      radius: m units*/
    Blackhole::Blackhole(double m6)
    {
        mass = m6;
        radius = calc_schwarzschild_radius();
    }

    /*Default constructor*/
    Blackhole::Blackhole()
    {
        mass = 1;
        calc_schwarzschild_radius();
    }

    /*in m*/
    double Blackhole::calc_schwarzschild_radius()
    {
        radius = 2.*G*M_SUN*mass*std::pow(10,6)/(c*c);
        return radius;
    }

    void Blackhole::init()
    {
        calc_schwarzschild_radius();
    }

/*Star class functions*/

    /*mass: M_sun units
      radius: R_sun units*/
    Star::Star(double m, double r, startype relation)
    {
        mass = m;
        radius = r;
        st = relation;
    }

    /*default constructor*/
    Star::Star()
    {
        mass = 1;
        radius = 1;
        st = no_relationship;
    }

    void Star::set_startype(std::string rad)
    {
        if(IsDouble(rad))
        {
            radius = std::stod(rad);
            st = no_relationship;
        }
        else
        {
            if(rad == "ms" || rad == "MS")
            {
                m_t_r = new Main_sequence_star;
                st = ms;
                radius = m_t_r->calc_radius(mass);
            }
            else if (rad == "wd" || rad == "WD")
            {
                m_t_r = new White_dwarf_star;
                st = wd;
                radius = m_t_r->calc_radius(mass);
            }
            else if(rad == "no_relationship")
            {
                st = no_relationship;
            }
            else
            {
                throw std::runtime_error("There is no type star with id: "+ rad +"\n");
            }
        }
    }

    void Star::init()
    {
        if(st != no_relationship)
        {
            radius = m_t_r->calc_radius(mass);
        }
    }

    void Star::print_startype()
    {
        if(st == no_relationship)
        {
            std::cout << "No relationship" << std::endl;
        }
        else if(st == wd)
        {
            std::cout << "White Dwarf" << std::endl;
        }
        else if(st == ms)
        {
            std::cout << "Main sequence" << std::endl;
        }
    }

/*Mass_to_radii classes and their calculations*/
    
    /*Main sequence star*/
    double Main_sequence_star::calc_radius(double mass)
    {
        return std::pow(mass, 0.6);
    }

    /*White dwarf star*/
    double White_dwarf_star::calc_radius(double mass)
    {
        return std::pow(mass, -1./3.)*1e-2;
    }


/*
    ******

        Parameters struct functions

    ******
*/

/*Parameters class functions*/

    /*Constructor*/
    Parameters::Parameters(int N_o, double e, double b, double fo, double fv_, double dist, double ink, double nu_, Star star, Blackhole blackhole)
    {
        N = N_o;
        eta = e;
        beta = b;
        fout = fo;
        fv = fv_;
        d = dist;
        i = ink;
        nu = nu_;
        s = star;
        bh = blackhole;
        calcfout = std::unique_ptr<Fout_calculation>(new No_fout_relation(*this));
        calcmdotpeak = std::unique_ptr<Mdot_peak_calculation>(new Mdotpeak_LodatoRossi(*this));

        if(d == 0)
        {
            oneper4pid2 = 1;
        }
        else
        {
            oneper4pid2 = 1./(4.*PI*d*d);
        }
        fv_check();
        calculate_rphtmin();
        calculate_tphtmin();
        calculate_tmin();
        calculate_rin();
        calculate_rout();
        calculate_rhalf();
        calculate_rL();
        calculate_mdot_edd();
        calcmdotpeak ->calc_mdotpeak();
        calculate_time_wind_limit();
        calculate_time_rph_limit();
        cosi = cos(i);
        diffusion_timescale = tmin*0.06;

        bolometric = false;
        extras = false;
        spectra = false;
        lcurve = false;
        diffusion = false;
        use_rph_limit = false;
    }

    /*Default constructor*/
    Parameters::Parameters()
    {
        N = 1000;
        eta = 0.1;
        beta = 1;
        fout = 0.1;
        fv = 1;
        d = 0;
        i = 0;
        nu = 6.3e14;
        oneper4pid2 = 1;
        calcfout = std::unique_ptr<Fout_calculation>(new No_fout_relation(*this));
        calcmdotpeak = std::unique_ptr<Mdot_peak_calculation>(new Mdotpeak_LodatoRossi(*this));

        fv_check();
        calculate_rphtmin();
        calculate_tphtmin();
        calculate_tmin();
        calculate_rin();
        calculate_rout();
        calculate_rhalf();
        calculate_rL();
        calculate_mdot_edd();
        calcmdotpeak -> calc_mdotpeak();
        calculate_time_wind_limit();
        calculate_time_rph_limit();
        cosi = 1;
        diffusion_timescale = tmin*0.06;


        t_start = tmin;
        t_end = 100;
        dt = 0.5;

        nu_start = wavelength_to_frequency(10000.);
        nu_end = wavelength_to_frequency(1.);
        dnu = (nu_end - nu_start)/10000.;
        time = tmin;

        lcurve = false;
        spectra = false;
        bolometric = false;
        extras = false;
        diffusion = false;
        use_rph_limit = false;

    }

    //Destructor
    Parameters::~Parameters(){}

    Parameters::Parameters(Parameters& copied)
    {
        N =    copied.N;
        eta =  copied.eta;
        beta = copied.beta;
        fout = copied.fout;
        fv =   copied.fv;
        d =    copied.d;
        i =    copied.i;
        diffusion_timescale = copied.diffusion_timescale;
        nu =   copied.nu;
        oneper4pid2 = copied.oneper4pid2;

        t_start = copied.t_start;
        t_end   = copied.t_end;
        dt      = copied.dt;
        nu      = copied.nu;

        /*Spectra intervalls and step*/
        nu_start = copied.nu_start;
        nu_end   = copied.nu_end;
        dnu      = copied.dnu;
        time     = copied.time;
        
        /*Derived quantities*/
        tmin     = copied.tmin;
        rphtmin  = copied.rphtmin;
        mdot_edd = copied.mdot_edd;
        tphtmin  = copied.tphtmin;
        rin      = copied.rin;
        rout     = copied.rout;
        rhalf    = copied.rhalf;
        rL       = copied.rL;
        time_wind_limit = copied.time_wind_limit;
        time_rph_limit = copied.time_rph_limit;
        cosi     = copied.cosi;
        Mdot_peak = copied.Mdot_peak;
        rpper3rs = copied.rpper3rs;
        tedge = copied.tedge;
        vwind = copied.vwind;

        /*Switch for bolometric lightcurve and extra informations*/
        bolometric = copied.bolometric;
        extras     = copied.extras;
        spectra    = copied.spectra;
        lcurve     = copied.lcurve;
        diffusion = copied.diffusion;
        use_rph_limit = copied.use_rph_limit;

        t_start_default = copied.t_start_default;
        t_end_relative_to_tmin = copied.t_end_relative_to_tmin;

        is_bifout = copied.is_bifout;

        s = copied.s;
        bh = copied.bh;
        calcfout.reset(copied.calcfout ->clone(*this));
        calcmdotpeak.reset(copied.calcmdotpeak -> clone(*this));
    }

    void Parameters::change_foutcalculation(Fout_calculation* new_fout)
    {
        calcfout.reset(new_fout);
    }

    void Parameters::change_mdotpeakcalculation(Mdot_peak_calculation* new_mdotpeak)
    {
        calcmdotpeak.reset(new_mdotpeak);
    }

    /*Functions used constructor*/
    void Parameters::fv_check()
	{
	    double vw, fv_max;
	    fv_max = c/(4.4e7*std::pow(beta, 0.5)*std::pow(bh.get_mass(),1./3.)*std::pow(s.get_radius(), -0.5)*std::pow(s.get_mass(),1./6.));
	    vw = 4.4e7*fv*std::pow(beta, 0.5)*std::pow(bh.get_mass(),1./3.)*std::pow(s.get_radius(), -0.5)*std::pow(s.get_mass(),1./6.);
        if(vw > c)
        {
            fv = fv_max;
        }
	}

    //It returns in day
       void Parameters::calculate_tmin()
       {
            tmin = 41.*std::pow(bh.get_mass(),0.5)*std::pow(s.get_mass(),-1.)*std::pow(beta,-3.)*std::pow(s.get_radius(),3./2.);
       }

       //in m unit
       void Parameters::calculate_rphtmin()
       {
            double rphtmin_newest;
            rphtmin_newest = KAPPA_ES * fout * Mdot_peak / (4 * M_PI * fv * c * sqrt(eta));
            
            rphtmin = rphtmin_newest;
       }

       //in K unit
       void Parameters::calculate_tphtmin()
       {
            double tltmin, tltmin4, tphtmin_newest,rl;
            
            tltmin4 = (Mdot_peak * eta * eta * eta * std::pow(c,6)) / (4. * M_PI * G * G * std::pow(bh.get_mass()*1e6 * M_SUN,2.) * SIGMA);

            tltmin = std::pow(tltmin4, 1./4.);
            rl = 1./(2.*eta)*bh.get_radius();
            tphtmin_newest = tltmin * std::pow(rphtmin/rl, -2./3.) *  std::pow(fout/fv, 1./3.);
            tphtmin = tphtmin_newest;
       }

       //in m
       void Parameters::calculate_rin()
       {
           rin = 3*bh.get_radius();
       }

       //in m
       void Parameters::calculate_rout()
       {
           double rt;
           rt = 0.47*std::pow(bh.get_mass()/s.get_mass(),1./3.)*AU_METER*s.get_radius();       //rt in m
           rout = 2*rt/beta;
       }

       void Parameters::calculate_rhalf()
       {
           rhalf = (rout - rin)/2.0 + rin;
       }

        // in m
       void Parameters::calculate_rL()
       {
           rL = bh.get_radius() / (2.0 * eta);
       }

        //in kg/s units
        void Parameters::calculate_mdot_edd()
        {
            mdot_edd = 1.3e21*bh.get_mass()*(0.1/eta);
        }

        void Parameters::calculate_rpper3rs()
        {
            double rp, rt;
            rt = std::pow(bh.get_mass()*1.e6/s.get_mass(), 1./3.) * s.get_radius() * R_SUN_METER; //in m
            rp = rt / beta; //in m
            rpper3rs = rp / (3. * bh.get_radius()); //dimensionless
        }

        //day unit
        void Parameters::calculate_tedge()
        {
            tedge = 1. * std::pow(fout, 3./8.) * std::pow(fv, -3./4.) * std::pow(bh.get_mass(), 5./8.) * std::pow(rpper3rs, 9./8.) * std::pow(s.get_mass(), 3./8.) * std::pow(s.get_radius(), -3./8.);
        }

        //in m/s
        void Parameters::calculate_vwind()
        {
            vwind = fv * std::sqrt(2.*G * bh.get_mass()*1e6 * M_SUN/(6.*bh.get_radius()*rpper3rs));
        }

        //in days
        void Parameters::calculate_time_wind_limit()
        {
            time_wind_limit = tmin * std::pow(Mdot_peak/mdot_edd, 3./5.);
        }

        void Parameters::calculate_time_rph_limit()
        {
            time_rph_limit = tmin * std::pow((rphtmin * 2. * eta)/bh.get_radius(), 3./5.);
        }

        void Parameters::print_parameters()
        {
            std::cout << "M6 = " << bh.get_mass() << "\tRS = " << bh.get_radius() << std::endl;
            std::cout << "Mstar = " << s.get_mass() << "\tRstar = " << s.get_radius() << std::endl;
            std::cout << "startype = ";
            s.print_startype();
            if(is_bifout)
            {
                std::cout << "fout: time dependent fout " << "\tfv = " << fv << "\teta = " << eta << "\tbeta = " << beta << std::endl;    
            }
            else
            {
                std::cout << "fout = " << fout << "\tfv = " << fv << "\teta = " << eta << "\tbeta = " << beta << std::endl;
            }
            std::cout << "d = " << d << "\ti = " << i << "\tnu = " << nu << "\tN = " << N << std::endl;
            std::cout << "tmin = " << tmin << "\ttstart = " << t_start << "\ttend = " << t_end << "\tdt = " << dt << std::endl; 
            std::cout << "rin = " << rin << "\trout = " << rout << std::endl;
            std::cout << "Mdot_peak = " << Mdot_peak << std::endl;
        }

        /*Comment line functions*/
        void Parameters::total_lc_commentline(std::ofstream& out)
        {
            out << "#M6 " << bh.get_mass() << std::endl;
            out << "#m* " << s.get_mass() << std::endl;
            out << "#x* " << s.get_radius() << std::endl;
            out << "#eta " << eta << std::endl;
            out << "#beta " << beta << std::endl;
            if(is_bifout)
            {
                out << "#used build-in fout calculation" << std::endl;
            }
            else
            {
                out << "#fout " << fout << std::endl;
            }
            out << "#fv " << fv << std::endl;
            out << "#d " << d << std::endl;
            out << "#i " << i << std::endl;
            out << "#nu " << nu << std::endl;
            out << "#lambda " << frequency_to_wavelength(nu) << std::endl;
            out << "#N " << N << std::endl;
            out << "#" << std::endl;
            out << "#time (day)" << "\t" << "Ltotal (erg/s/Hz)\t" << "Lwind (erg/s/Hz)" << "\t" << "Ldisk (erg/s/Hz)" << std::endl;
        }

        void Parameters::lbol_commentline(std::ofstream & out)
        {
            out << "#M6 " << bh.get_mass() << std::endl;
            out << "#m* " << s.get_mass() << std::endl;
            out << "#x* " << s.get_radius() << std::endl;
            out << "#eta " << eta << std::endl;
            out << "#beta " << beta << std::endl;
            if(is_bifout)
            {
                out << "#used build-in fout calculation" << std::endl;
            }
            else
            {
                out << "#fout " << fout << std::endl;
            }
            out << "#fv " << fv << std::endl;
            out << "#d " << d << std::endl;
            out << "#i " << i << std::endl;
            out << "#N " << N << std::endl;
            out << "#" << std::endl;
            out << "#time (day)" << "\t" << "Lbol (erg/s)" << std::endl;
        }

        void Parameters::extras_commentline(std::ofstream & out)
        {
            out << "#M6 " << bh.get_mass() << std::endl;
            out << "#m* " << s.get_mass() << std::endl;
            out << "#x* " << s.get_radius() << std::endl;
            out << "#eta " << eta << std::endl;
            out << "#beta " << beta << std::endl;
            if(is_bifout)
            {
                out << "#used build-in fout calculation" << std::endl;
            }
            else
            {
                out << "#fout " << fout << std::endl;
            }
            out << "#fv " << fv << std::endl;
            out << "#d " << d << std::endl;
            out << "#i " << i << std::endl;
            out << "#N " << N << std::endl;
            out << "#rhalf " << rhalf << std::endl; 
            out << "#" << std::endl;
            out << "#time (day)" << "\t" << "Photosperic radius (m)" << "\t" << "Photospheric temperature (K)" << "\t" << "T_disk(rhalf) (K)" << "\tfout" <<"\ttedge (day)" << std::endl;
        }

        void Parameters::spectra_commentline(std::ofstream & out)
        {
            out << "#M6 " << bh.get_mass() << std::endl;
            out << "#m* " << s.get_mass() << std::endl;
            out << "#x* " << s.get_radius() << std::endl;
            out << "#eta " << eta << std::endl;
            out << "#beta " << beta << std::endl;
            if(is_bifout)
            {
                out << "#used build-in fout calculation" << std::endl;
            }
            else
            {
                out << "#fout " << fout << std::endl;
            }
            out << "#fv " << fv << std::endl;
            out << "#d " << d << std::endl;
            out << "#i " << i << std::endl;
            out << "#N " << N << std::endl;
            out << "#time " << time << std::endl;
            out << "#" << std::endl;
            out << "#nu (Hz)" << "\tlambda (nm)" << "\tLwind (erg/s/Hz)" << "\tLdisk (erg/s/Hz)" << std::endl;
        }

        void Parameters::init(double dtime)
        {
            s.init();
            bh.init();
            calculate_tmin();
            calculate_mdot_edd();
            calcmdotpeak ->calc_mdotpeak();
            calcfout->calculate_fout(tmin);
            fv_check();
            calculate_rphtmin();
            calculate_tphtmin();
            calcfout->calculate_fout(dtime);
            calculate_rin();
            calculate_rout();
            calculate_rhalf();
            calculate_rL();
            calculate_rpper3rs();
            calculate_tedge();
            calculate_vwind();
            calculate_time_wind_limit();
            calculate_time_rph_limit();
            if(t_start_default)
            {
                t_start = tmin;
            }
            if(t_end_relative_to_tmin)
            {
                t_end += tmin;
            }
        }   

/*
    ********

        Fout_calculation functions: mother class: onli constructor and destructor

            No_fout_relation: derivated class: functions (clone, calculate_fout)

            Fout_Lodato_Rossi_eq28: derivated class functions, time dependent fout

    ********
*/

/*Fout_claculation classes and functions*/
Fout_calculation::Fout_calculation(){}
Fout_calculation::~Fout_calculation(){}

    /*No_realtion*/
    No_fout_relation::No_fout_relation(Parameters& _param_obj) : Fout_calculation(), param_obj(_param_obj) {}
    void No_fout_relation::calculate_fout([[maybe_unused]] double time){}
    Fout_calculation* No_fout_relation::clone(Parameters& where)
    {
        return new No_fout_relation(where);
    }
    No_fout_relation::~No_fout_relation(){}

    /*Fout_Lodato_Rossi_eq28*/
    Fout_Lodato_Rossi_eq28::Fout_Lodato_Rossi_eq28(Parameters& _param_obj) : Fout_calculation(), param_obj(_param_obj){}
    void Fout_Lodato_Rossi_eq28::calculate_fout(double time)
    {
        double fout, Mdot_rate, Mdotfb;
        Mdotfb = param_obj.Mdot_peak * std::pow(time/param_obj.tmin, -5./3.);
        Mdot_rate = Mdotfb / param_obj.mdot_edd;
        if(Mdot_rate < 2.188)
        {
            fout = 0.1;
        }
        else
        {
            fout = 2./M_PI * std::atan(1./7.5 * (Mdot_rate - 1.));
        }
        param_obj.fout = fout;
    }
    Fout_calculation* Fout_Lodato_Rossi_eq28::clone(Parameters& where)
    {
        return new Fout_Lodato_Rossi_eq28(where);
    }
    Fout_Lodato_Rossi_eq28::~Fout_Lodato_Rossi_eq28(){}

/*
    *********

        Mdot_peak calculation classes, its derivated classes functions 
        (Mdotpeak_LodatoRossi, Mdotpeak_GuillochonRamirez)

    *********
*/


/*Mdot_peak_calculation classes*/
    
    Mdot_peak_calculation::Mdot_peak_calculation(){}
    Mdot_peak_calculation::~Mdot_peak_calculation(){}
    
    Mdotpeak_LodatoRossi::Mdotpeak_LodatoRossi(Parameters& _param_obj) : param_obj(_param_obj){}
    
    /* kg/s */
    void Mdotpeak_LodatoRossi::calc_mdotpeak()
    {
        double mdotpeak;
        mdotpeak = 1./3.*param_obj.s.get_mass()*M_SUN/(param_obj.tmin*24.*60.*60.);
        param_obj.Mdot_peak = mdotpeak;
    }
    void Mdotpeak_LodatoRossi::set_politrop([[maybe_unused]]std::string gamma){}
    void Mdotpeak_LodatoRossi::set_default_politrop(){}
    Mdot_peak_calculation* Mdotpeak_LodatoRossi::clone(Parameters& where)
    {
        return new Mdotpeak_LodatoRossi(where);
    }
    Mdotpeak_LodatoRossi::~Mdotpeak_LodatoRossi(){}


    Mdotpeak_GuillochonRamirez::Mdotpeak_GuillochonRamirez(Parameters& _param_obj) : param_obj(_param_obj) {}
    double Mdotpeak_GuillochonRamirez::Mdotpeak_politrop4per3()
    {
        double agamma, exp_body, mdotpeak, beta;
        beta = param_obj.beta;
        exp_body = (27.261 - 27.516*beta + 3.8716 * beta * beta)/(1. - 3.2605 * beta - 1.3865 * beta * beta);
        agamma = std::exp(exp_body);
        mdotpeak = agamma * std::pow(param_obj.bh.get_mass(), -0.5) * std::pow(param_obj.s.get_mass(), 2.) * std::pow(param_obj.s.get_radius(), -1.5);
        mdotpeak = mdotpeak * M_SUN / YEAR_S; //Msun/year -> kg/s
        return mdotpeak;
    }

    double Mdotpeak_GuillochonRamirez::Mdotpeak_politrop5per3()
    {
        double agamma, exp_body, mdotpeak, beta;
        beta = param_obj.beta;
        exp_body = (10. - 17. * beta + 6. * beta * beta)/(1. - 0.47 * beta - 4.5 * beta * beta);
        agamma = std::exp(exp_body);
        mdotpeak = agamma * std::pow(param_obj.bh.get_mass(), -0.5) * std::pow(param_obj.s.get_mass(), 2.) * std::pow(param_obj.s.get_radius(), -1.5);
        mdotpeak = mdotpeak * M_SUN / YEAR_S; //Msun/year -> kg/s
        return mdotpeak;
    }

    void Mdotpeak_GuillochonRamirez::set_default_politrop()
    {
        if(param_obj.s.get_startype() == wd)
        {
            set_politrop("5per3");
        }
        else
        {
            set_politrop("4per3");
        }
    }

    void Mdotpeak_GuillochonRamirez::calc_mdotpeak()
    {
        if(politrop == "5per3")
        {
            param_obj.Mdot_peak = Mdotpeak_politrop5per3();
        }
        else
        {
            param_obj.Mdot_peak = Mdotpeak_politrop4per3();
        }
    }

    Mdot_peak_calculation* Mdotpeak_GuillochonRamirez::clone(Parameters& where)
    {
        return new Mdotpeak_GuillochonRamirez(where);
    }

    void Mdotpeak_GuillochonRamirez::set_politrop(std::string gamma)
    {
        politrop = gamma;
    }

    Mdotpeak_GuillochonRamirez::~Mdotpeak_GuillochonRamirez(){}


/*
    *********

        Here will start the different part of the light curve calculation: Wind and disk part

    *********
*/

/*
    *********Wind part of the light curve functions*********
*/

        /*Constructor*/
        Wind_part_of_lightcurve::Wind_part_of_lightcurve(Parameters& p) : par(p) {}
       

        /*Default constructor*/
        Wind_part_of_lightcurve::Wind_part_of_lightcurve() : par(*new Parameters){}

        /*Other functions*/

        double Wind_part_of_lightcurve::calc_temperature_ph(double t)
        {
            double tphtime;
            tphtime = par.tphtmin * std::pow(t/par.tmin, 25./36.);
            return tphtime;
        }

        //in m
        double Wind_part_of_lightcurve::calc_radius_ph(double t)
        {
            double rphtime;
            rphtime = par.rphtmin * std::pow(t/par.tmin,-5./3.);
            return rphtime;
        }

        //in K
        double Wind_part_of_lightcurve::calc_Tph_at_rL(double t)
        {
            double TL, TLtmin, TL4,Tph;
            TL4 = par.Mdot_peak * std::pow(par.eta, 3) * std::pow(c, 6)/ (4. * M_PI * G * G * std::pow(par.bh.get_mass()*1e6 * M_SUN, 2) * SIGMA);
            TLtmin = std::pow(TL4, 1./4.);
            TL = TLtmin * std::pow(t/par.tmin, -5./12.);
            Tph = TL * std::pow(par.fout/par.fv, 1./3.);
            return Tph;
        }

        double Wind_part_of_lightcurve::calc_luminosity_at_rldist_nu(double t, double nu)
        {
            double lnu;
            photospheric_temperature = calc_Tph_at_rL(t);
            photospheric_rad = par.rL;
            lnu = 4.*PI*PI*std::pow(photospheric_rad,2.)*Planck_function_for_nu(nu, photospheric_temperature);
            lnu = lnu * par.oneper4pid2;        //take into account of distance
            return lnu*std::pow(10.,7.);    //J -> erg
        }

        //will return in erg unit
        double Wind_part_of_lightcurve::calc_luminosity_at_nu(double t, double nu)
        {
            double lnu;
            photospheric_temperature = calc_temperature_ph(t);
            photospheric_rad = calc_radius_ph(t);
            lnu = 4.*PI*PI*std::pow(photospheric_rad,2.)*Planck_function_for_nu(nu, photospheric_temperature);
            lnu = lnu * par.oneper4pid2;        //take into account of distance
            return lnu*std::pow(10.,7.);    //J -> erg
        }

/*
    **********Disk part of the lightcurve*******
*/
    
    /*calc_Disk_at_r class functions*/

        /*Constructor*/
        calc_Disk_at_r::calc_Disk_at_r(Parameters& p) : par(p){}

        calc_Disk_at_r::calc_Disk_at_r() : par(*new Parameters){}
        
        void calc_Disk_at_r::set_time(double t)
        {
            time = t;
            mdot = calc_Mdot(time);
        }

        void calc_Disk_at_r::set_nu(double new_nu)
        {
            own_nu = new_nu;
        }
        
        //Retrurns kg/s unit
        double calc_Disk_at_r::calc_Mdot(double t)
        {
            double mdotp, mdotfb;
            mdotp = par.Mdot_peak;
            mdotfb = mdotp*std::pow(t/par.tmin,-5./3.);
            mdot = (1. - par.fout)*mdotfb;
            return mdot;
        }

        double calc_Disk_at_r::calc_Td(double r)
        {
            double Td, f, factor1, factor2, factor_help, Td4;
            f = 1. - std::pow(par.rin/r, 0.5);
            factor1 = 3.*G*par.bh.get_mass()*1.0e6*M_SUN*mdot*f/(8.*M_PI*r*r*r);
            factor_help = 1./4. + 3./2.*f*std::pow((mdot*par.bh.get_radius())/(par.eta*par.mdot_edd*r), 2);
	        factor2 = 1./2. + std::pow(factor_help, 0.5);
            Td4  = factor1/(factor2 * SIGMA);
            Td = std::pow(Td4, 1./4.);
            return Td;
        }

        double calc_Disk_at_r::calc_Td_at_rhalf()
        {
            return calc_Td(par.rhalf);
        }

        double calc_Disk_at_r::calc_F_at_r(double r, double nu)
        {
            double F;
            F = 2*M_PI*r*Planck_function_for_nu(nu, calc_Td(r));
            return F;
        }

        double calc_Disk_at_r::calc_Fx(double x)
        {
            return calc_F_at_r(x, own_nu);
        }


    /*Disk_part_of_lightcurve functions*/

        /*Constructor*/
        Disk_part_of_lightcurve::Disk_part_of_lightcurve(Parameters& p) : calc_disk_r(p)
        {
            calc_disk_r.set_nu(p.nu);
            sum_lum_disk = Trapezoidal_rule(calc_disk_r, calc_disk_r.par.rin, calc_disk_r.par.rout, calc_disk_r.par.N);
        }

        /*Default constructor*/
        Disk_part_of_lightcurve::Disk_part_of_lightcurve()
        {
            sum_lum_disk = Trapezoidal_rule(calc_disk_r, calc_disk_r.par.rin, calc_disk_r.par.rout, calc_disk_r.par.N);
        }

        /*Functions*/
        double Disk_part_of_lightcurve::calc_disk_part_at_t(double time)
        {
            double sum_disk = 0;
            calc_disk_r.set_time(time);
            sum_lum_disk.set_newfx(calc_disk_r);
            sum_lum_disk.set_ranges(calc_disk_r.par.rin, calc_disk_r.par.rout);
            if(calc_disk_r.par.rin < calc_disk_r.par.rout)
            {
                sum_disk = sum_lum_disk.calc_integral();
                sum_disk = sum_disk * 1e7 * calc_disk_r.par.cosi * calc_disk_r.par.oneper4pid2;  //J -> erg and take into account inclination and distance
            }
            else {sum_disk = 0;}
            return sum_disk;
        }

        double Disk_part_of_lightcurve::calc_disk_part_at_t_nu(double time, double nu)
        {
            calc_disk_r.set_nu(nu);
            return calc_disk_part_at_t(time);
        }

        void Disk_part_of_lightcurve::set_nu(double new_nu)
        {
            calc_disk_r.set_nu(new_nu);
        }

        double Disk_part_of_lightcurve::get_Tdisk_at_rhalf(double time)
        {
            calc_disk_r.set_time(time);
            return calc_disk_r.calc_Td_at_rhalf();
        }

/*
    **********

        Different light curve calculation classes

        Total_lightcurve:
            calculate the monochromatic light curve of a tde event at a given time intervall

        Diffusion calculation:
            To calculate diffusion, we need integral.
            - Lum_diffusion_fx class: derivated from Fx class
                It will calculate the needed integral
            - Diffusion_luminosity class: it will calculate the diffused luminosity at a given time
                To speeds up the integration time it uses the previous result of the integral (previous time step result)
            - Diffusion_luminosity_timeintervall class:
                it will calculate the monochromatic diffusion light curve of a tde event at a given time intervall

        Bolometric luminosity
            It is again an integral
            - Lsum class: derivated from Fx class
                It will calculate the integral of luminosities at different frequencies
            - Lbol class: 
                It will calculate the bolometric luminosity at a given time intervall

        Spectrum of a TDE event
            Spectrum class:
                Will calculate the spectrum of an event to a given time and given frequency intervall

    **********
*/

/*Total_lightcurve class functions*/

    /*Constructor*/
    Total_lightcurve::Total_lightcurve(Parameters& p) : par(p), wind_part(p), disk_part(p){}

    /*Functions*/
    void Total_lightcurve::output_commentline(std::ofstream & out)
    {
        par.total_lc_commentline(out);
    }

    void Total_lightcurve::print_output_before_rphlimit(double time, std::ofstream & out)
    {
        double lwind, ldisk;
        par.calcfout->calculate_fout(time);
        par.calculate_tedge();
        lwind = wind_part.calc_luminosity_at_nu(time, par.nu);
        ldisk = disk_part.calc_disk_part_at_t(time);
        out << time << "\t" << lwind + ldisk << "\t" << lwind << "\t" << ldisk << std::endl;
    }

    void Total_lightcurve::print_output_after_rphlimit(double time, std::ofstream & out)
    {
        double lwind, ldisk;
        par.calcfout->calculate_fout(time);
        par.calculate_tedge();
        lwind = wind_part.calc_luminosity_at_rldist_nu(time, par.nu);
        ldisk = disk_part.calc_disk_part_at_t(time);
        out << time << "\t" << lwind + ldisk << "\t" << lwind << "\t" << ldisk << std::endl;
    }

    void Total_lightcurve::calculate_lightcurve_timeintervall(double tmin, double tmax, double dt)
    {
        double t;//, lwind, ldisk;
        std::ofstream lum;
	    lum.open("lcurve.dat");
        output_commentline(lum);
        disk_part.set_nu(par.nu);

        if(par.use_rph_limit)
        {
            if(par.extras)
            {
                std::ofstream extra;
                extra.open("extras.dat");
                par.extras_commentline(extra);
                for(t = tmin; t <= par.time_rph_limit && t <= tmax; t += dt)
                {
                    print_output_before_rphlimit(t, lum);
                    extra << t << "\t" << wind_part.calc_radius_ph(t) << "\t" << wind_part.calc_temperature_ph(t) << "\t" << disk_part.get_Tdisk_at_rhalf(t) << "\t" << par.fout << "\t" << par.tedge << std::endl;
                }
                for(; t <= tmax; t+= dt)
                {
                    print_output_after_rphlimit(t, lum);
                    extra << t << "\t" << par.rL << "\t" << wind_part.calc_Tph_at_rL(t) << "\t" << disk_part.get_Tdisk_at_rhalf(t) << "\t" << par.fout << "\t" << par.tedge << std::endl;
                }
            }
            else
            {
                for(t = tmin; t <= par.time_rph_limit && t <= tmax; t += dt)
                {
                    print_output_before_rphlimit(t, lum);
                }
                for(; t <= tmax; t+= dt)
                {
                    print_output_after_rphlimit(t, lum);
                }
            }
        }
        else
        {
            if(par.extras)
            {
                std::ofstream extra;
                extra.open("extras.dat");
                par.extras_commentline(extra);
                for(t = tmin; t <= tmax; t+= dt)
                {
                    print_output_before_rphlimit(t, lum);
                    extra << t << "\t" << wind_part.calc_radius_ph(t) << "\t" << wind_part.calc_temperature_ph(t) << "\t" << disk_part.get_Tdisk_at_rhalf(t) << "\t" << par.fout << "\t" << par.tedge << std::endl;
                }
            }
            else
            {
                for(t = tmin; t <= tmax; t += dt)
                {
                    print_output_before_rphlimit(t, lum);
                }
            }
        }
        lum.close();
    }

/*Diffusion calculations*/
    /*Fx derivation class: Lum_diffusion_fx functions*/
        Lum_diffusion_fx::Lum_diffusion_fx(Parameters& p) : par(p), disk_part(p), wind_part(p)
        {}

        Lum_diffusion_fx::Lum_diffusion_fx() : par(*new Parameters), disk_part(par), wind_part(par)
        {}

        double Lum_diffusion_fx::L_inp(double t_offset)
        {
            double Linp;
            Linp = par.eta * par.Mdot_peak * std::pow(1 + t_offset/par.tmin, -5./3.) * c*c; //kg*m2/s3
            Linp *= 1e7; // -> to erg/s
            return Linp;
        }

        double Lum_diffusion_fx::L_inp_new(double t)
        {
            double l_wind, l_disk, Linp, t_offset;
            t_offset = t + par.tmin;
            par.calcfout->calculate_fout(t_offset);
            l_wind = wind_part.calc_luminosity_at_nu(t_offset, par.nu);
            l_disk = disk_part.calc_disk_part_at_t(t_offset);
            Linp = l_disk + l_wind;
            return Linp;
        }

        double Lum_diffusion_fx::calc_Lt_integral(double t)
        {
            double res;
            res = std::exp(t/par.diffusion_timescale) * L_inp_new(t);
            return res;
        }

        double Lum_diffusion_fx::calc_Fx(double x)
        {
            return calc_Lt_integral(x);
        }

    /*Diffusion_luminosity class: calculate diffused luminosity at t time*/

        Diffusion_luminosity::Diffusion_luminosity(Parameters& p) : par(p), diff_fx(p)
        {
            time = par.tmin;
            time_offset = time - par.tmin;
            diffusion = Trapezoidal_rule(diff_fx, 0., time_offset, 0.01, 1);
            previous_integral = 0.;
            previous_time = 0.;
        }

        Diffusion_luminosity::Diffusion_luminosity() : par(*new Parameters), diff_fx(par)
        {
            time = par.tmin;
            time_offset = time - par.tmin;
            diffusion = Trapezoidal_rule(diff_fx, 0., time_offset, 0.1, 1);
            previous_integral = 0.;
            previous_time = 0.;
        }

        void Diffusion_luminosity::set_time(double t)
        {
            time = t;
            time_offset = time - par.tmin;
            diffusion.set_xend(time_offset);
            diffusion.set_xstart(previous_time);
        }

        double Diffusion_luminosity::calc_Lum_at_t()
        {
            double res, integral;
            integral = diffusion.calc_integral_dx() + previous_integral;
            res = std::exp(-time_offset/par.diffusion_timescale) / par.diffusion_timescale * integral;
            previous_integral = integral;
            previous_time = time_offset;
            return res;
        }

        double Diffusion_luminosity::get_Linp(double t)
        {
            return diff_fx.L_inp_new(t);
        }

        double Diffusion_luminosity::get_integral(double t)
        {
            diffusion.set_xend(t - par.tmin);
            return diffusion.calc_integral_dx();
        }

    /*Diffusion calculations for a given time intervall*/

        Diffusion_luminosity_timeintervall::Diffusion_luminosity_timeintervall(Parameters& p) : par(p), lum_diff(p)
        {}

        Diffusion_luminosity_timeintervall::Diffusion_luminosity_timeintervall() : par(*new Parameters), lum_diff(par)
        {}

        void Diffusion_luminosity_timeintervall::calc_diffusion_luminosity_timeintervall(double tmin, double tmax, double dt)
        {
            double t, lum;

            std::ofstream lum_file;
	        lum_file.open("lcurve_diffusion.dat");
            lum_file << "#tdiff\t" << par.diffusion_timescale << std::endl;
            par.total_lc_commentline(lum_file);
            for(t = tmin; t < tmax; t += dt)
            {
                par.calcfout->calculate_fout(t);
                par.init(t);
                lum_diff.set_time(t);
                lum = lum_diff.calc_Lum_at_t();
                lum_file << t << "\t" << lum << std::endl;
            }
        }



    

/*Bolometric luminosity:*/
    /*Lsum class functions*/

        /*Constructor*/
        Lsum::Lsum(Parameters& p, double t) : par(p), wind_part(p), disk_part(p), time(t){}

        Lsum::Lsum(Parameters& p) : par(p), wind_part(p), disk_part(p)
        {
            time = 0;
        }

        Lsum::Lsum() : par(*new Parameters), wind_part(par), disk_part(par)
        {
            time = 0;
        }

        /*Functions*/
        void Lsum::set_time(double t)
        {
            time = t;
        }

        void Lsum::set_nu(double nu)
        {
            par.nu = nu;
        }

        void Lsum::comment(std::ofstream & out)
        {
            par.lbol_commentline(out);
        }

        double Lsum::calc_Lsum(double nu)
        {
            double Lwind, Ldisk, L_sum;
            set_nu(nu);
            Lwind = wind_part.calc_luminosity_at_nu(time, nu);
            Ldisk = disk_part.calc_disk_part_at_t(time);
            L_sum = Lwind + Ldisk;
            return L_sum;
        }

        double Lsum::calc_Fx(double x)
        {
            return calc_Lsum(x);
        }


    /*Lbol class functions*/

        /*Constructor*/
        Lbol::Lbol(Parameters& p) : par(p), Ls(p)
        {
            bolometric = Trapezoidal_rule(Ls, wavelength_to_frequency(2500), wavelength_to_frequency(100), 24);
        }

        double Lbol::calc_Lbol(double time)
        {
            double lbol;
            Ls.set_time(time);
            bolometric.set_newfx(Ls);
            lbol = bolometric.calc_integral();
            return lbol;
        }

        void Lbol::comment_line(std::ofstream & out)
        {
            par.lbol_commentline(out);
        }

        void Lbol::Lbol_time(double tstart, double tend, double dt)
        {
            std::ofstream bol;
            bol.open("L_bolometric.dat");
            comment_line(bol);
            double Lbol_now, time_now;
            for(time_now = tstart; time_now <= tend; time_now+=dt)
            {
                Lbol_now = calc_Lbol(time_now);
                bol << time_now << "\t" << Lbol_now << std::endl;
            }
        }

/*Spectrum class functions*/
    Spectrum::Spectrum(Parameters& p) : par(p), wind_part(p), disk_part(p){}

    void Spectrum::output_commentline(std::ofstream & out)
    {
        par.spectra_commentline(out);
    }

    void Spectrum::calculate_spectrum_nuintervall(double nu_min, double nu_max, double dnu)
    {
        double nu, lwind, ldisk;
        std::ofstream spectra;
	    spectra.open("spectra.dat");
        output_commentline(spectra);
        for(nu = nu_min; nu <= nu_max; nu+= dnu)
        {
            disk_part.set_nu(nu);
            lwind = wind_part.calc_luminosity_at_nu(par.time, nu);
            ldisk = disk_part.calc_disk_part_at_t(par.time);
            spectra << nu << "\t" << frequency_to_wavelength(nu) << "\t" << lwind << "\t" << ldisk << std::endl;
        }
    }
