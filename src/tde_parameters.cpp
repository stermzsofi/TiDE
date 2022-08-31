#include "tde_parameters.hpp"

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

    void Blackhole::print_settings()
    {
        std::cout << "Parameters of the black hole:" << std::endl;
        std::cout << "\tMass: " << mass << " 10^6 M_sun\tSchwarzschild radius: " << radius << " m" << std::endl; 
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

    Star::Star(double m, double r, startype relation, std::string _politrop)
    {
        mass = m;
        radius = r;
        st = relation;
        politrop = _politrop;
    }

    /*default constructor*/
    Star::Star()
    {
        mass = 1;
        radius = 1;
        st = no_relationship;
        politrop = "4per3";
        n = 3;
        rho_c = RHO_C_N3_CONSTPART;
        K = K_N3_CONSTPART;
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
                throw std::runtime_error("ERROR! There is no type star with id: "+ rad +"\n");
            }
        }
    }

    void Star::init()
    {
        if(st != no_relationship)
        {
            radius = m_t_r->calc_radius(mass);
        }
        if(!setted_politrop)
        {
            if(st == no_relationship || st == ms)
            {
                politrop = "4per3";
                n = 3;
                rho_c = RHO_C_N3_CONSTPART;
                K = K_N3_CONSTPART;
            }
            else if(st == wd)
            {
                politrop = "5per3";
                n = 1.5;
                rho_c = RHO_C_N3PER2_CONSTPART;
                K = K_N3PER2_CONSTPART;
            }
            
        }
        rho_c *= mass * M_SUN * std::pow(radius * R_SUN_METER,-3);
        K *= std::pow(radius * R_SUN_METER, 2) * std::pow(rho_c, 1.0 - 1./n);
    }

    void Star::print_startype()
    {
        std::cout << "\tUsed mass-radius relation: ";
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
        //std::cout << "gamma = " << politrop << "\tn = " << n << "\t Kconst = " << K << "\trhoc = " << rho_c << std::endl;
    }

    void Star::set_politrop(std::string _politrop)
    {
        if("4per3" == _politrop || "5per3" == _politrop)
        {
            politrop = _politrop;
            setted_politrop = true;
            if("4per3" == politrop)
            {
                n = 3;
                rho_c = RHO_C_N3_CONSTPART;
                K = K_N3_CONSTPART;
            }
            else
            {
                n = 1.5;
                rho_c = RHO_C_N3PER2_CONSTPART;
                K = K_N3PER2_CONSTPART;
            }
        }
        else
        {
            throw std::runtime_error("ERROR! Unknown politrop parameter: "+_politrop);
        }
        
    }

    void Star::print_settings()
    {
        std::cout << "Parameters of the star:" << std::endl;
        std::cout << "\tMass: " << mass << " M_sun\tradius: " << radius << " R_sun" << std::endl;
        print_startype();
        std::cout << "\tPolitrop: " << politrop << "\tn: " << n << std::endl;
        std::cout << "\tK: " << K << "\trho_c: " << rho_c << std::endl;
    }

    std::string Star::get_politrop()
    {
        return politrop;
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
        beta_fulldisrupt = 1.;
        calcfout = std::unique_ptr<Fout_calculation>(new No_fout_relation(*this));
        calcmdotpeak = std::unique_ptr<Mdot_peak_calculation>(new Mdotpeak_LodatoRossi(*this));
        calctmin = std::unique_ptr<tmin_calculation>(new tmin_with_rt(*this));
        calcMdotfb = std::unique_ptr<Mdotfb_calculation>(new Mdotfb_classical_powerfunction(*this));

        if(d == 0)
        {
            oneper4pid2 = 1;
        }
        else
        {
            oneper4pid2 = 1./(4.*PI*d*d);
        }
        calculate_rL();
        fv_check();
        calculate_rt();
        calculate_rphtmin();
        calculate_tphtmin();
        calctmin -> calc_tmin();
        calculate_rin();
        calculate_rout();
        calculate_rhalf();
        
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
        beta_fulldisrupt = 1.;
        calcfout = std::unique_ptr<Fout_calculation>(new No_fout_relation(*this));
        calcmdotpeak = std::unique_ptr<Mdot_peak_calculation>(new Mdotpeak_LodatoRossi(*this));
        calctmin = std::unique_ptr<tmin_calculation>(new tmin_with_rt(*this));
        calcMdotfb = std::unique_ptr<Mdotfb_calculation>(new Mdotfb_classical_powerfunction(*this));

        calculate_rL();
        fv_check();
        calculate_rt();
        calculate_rphtmin();
        calculate_tphtmin();
        calctmin -> calc_tmin();
        calculate_rin();
        calculate_rout();
        calculate_rhalf();
        
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
        beta_fulldisrupt = copied.beta_fulldisrupt;
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
        rt = copied.rt;
        rp = copied.rp;
        tmin     = copied.tmin;
        rphtmin  = copied.rphtmin;
        mdot_edd = copied.mdot_edd;
        tphtmin  = copied.tphtmin;
        rin      = copied.rin;
        rout     = copied.rout;
        rhalf    = copied.rhalf;
        rL       = copied.rL;
        vw       = copied.vw; 
        time_wind_limit = copied.time_wind_limit;
        time_rph_limit = copied.time_rph_limit;
        cosi     = copied.cosi;
        Mdot_peak = copied.Mdot_peak;

        Mdotfb_t = copied.Mdotfb_t;

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
        calctmin.reset(copied.calctmin -> clone(*this));
        calcMdotfb.reset(copied.calcMdotfb -> clone(*this));
    }

    void Parameters::change_foutcalculation(Fout_calculation* new_fout)
    {
        calcfout.reset(new_fout);
    }

    void Parameters::change_mdotpeakcalculation(Mdot_peak_calculation* new_mdotpeak)
    {
        calcmdotpeak.reset(new_mdotpeak);
    }

    void Parameters::change_tmincalculation(tmin_calculation* new_tmincalc)
    {
        calctmin.reset(new_tmincalc);
    }

    void Parameters::change_Mdotfbcalculation(Mdotfb_calculation* new_calcMdotfb)
    {
        calcMdotfb.reset(new_calcMdotfb);
    }

    /*Functions used constructor*/
    void Parameters::fv_check()
	{
	    double fv_max;
        fv_max = c/(std::sqrt((2.0*G*bh.get_mass()*1e6*M_SUN)/rL));
        if(fv > fv_max)
        {
            fv = fv_max;
        }
        vw = fv * std::sqrt((2.0*G*bh.get_mass()*1e6*M_SUN)/rL);
	}

    //rt and rp in m unit
    void Parameters::calculate_rt()
    {
        rt = s.get_radius() * std::pow((bh.get_mass() * 1e6)/s.get_mass(), 1./3.) * R_SUN_METER;
        rp = rt / beta;
        rt /= beta_fulldisrupt;
    }

       //in m unit
       void Parameters::calculate_rphtmin()
       {
            double _rphtmin;
            _rphtmin = KAPPA_ES * fout * Mdot_peak/(4.0 * M_PI * vw);
            rphtmin = _rphtmin;
       }

       //in K unit
       void Parameters::calculate_tphtmin()
       {
            
            double _tltmin4, _tltmin;
            double _tphtmin;
            _tltmin4 = (G * bh.get_mass()*1e6 * M_SUN * Mdot_peak)/(4.0 * M_PI * std::pow(rL,3)*SIGMA);
            _tltmin = std::pow(_tltmin4, 1./4.);
            _tphtmin = _tltmin * std::pow(rphtmin/rL, -2./3.) *  std::pow(fout/fv, 1./3.);
            tphtmin = _tphtmin;
       }

       //in m
       void Parameters::calculate_rin()
       {
           rin = 3*bh.get_radius();
       }

       //in m
       void Parameters::calculate_rout()
       {
           rout = 2.0 * rp;
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


        //ez az időpont limitek, ahol a be illetve ki kellene kapcsolni valamit. nem igazán használom őket
        //in days
        //ki kellene szedni, amúgy is csak a sima hatványfüggvény akkréciós rátára igaz
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
            bh.print_settings();
            s.print_settings();
            calcfout ->print_settings();
            std::cout << "fv = " << fv << "\teta = " << eta << "\tbeta = " << beta << std::endl;
            std::cout << "d = " << d << "\ti = " << i << "\tnu = " << nu << " Hz\tN = " << N << std::endl;
            calctmin ->print_setted_type();
            std::cout << "tstart = " << t_start << " day\ttend = " << t_end << " day\tdt = " << dt <<" day" << std::endl; 
            std::cout << "rin = " << rin << " m\trout = " << rout << " m" << std::endl;
            calcMdotfb -> print_setted_type();
            std::cout << "Mdot_peak = " << Mdot_peak << " kg/s" << std::endl;
            //std::cout << "ninf = " << calcMdotfb -> get_ninf() << std::endl;
            std::cout << "vwind = " << vw << " m/s" << std::endl;
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
            calculate_rt();
            calctmin -> calc_tmin();
            calculate_mdot_edd();
            calcMdotfb->init();
            calcmdotpeak ->calc_mdotpeak();
            calculate_timedependent_parameters(tmin);
            calculate_rL();
            fv_check();
            calculate_rphtmin();
            calculate_tphtmin();
            
            calculate_timedependent_parameters(dtime);
            calculate_rin();
            calculate_rout();
            calculate_rhalf();
            
            calculate_time_wind_limit();
            calculate_time_rph_limit();
            if(t_start_default)
            {
                t_start = tmin;
            }
            else
            {
                if(t_start < tmin)
                {
                    std::__throw_runtime_error("ERROR! t_start is less than tmin.");
                }
            }
            if(t_end_relative_to_tmin)
            {
                t_end += tmin;
            }
        } 

        void Parameters::calculate_timedependent_parameters(double t)
        {
            calcMdotfb -> calc_Mdotfb_at_t(t);
            calcfout -> calculate_fout();
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
    void No_fout_relation::calculate_fout(){}
    Fout_calculation* No_fout_relation::clone(Parameters& where)
    {
        return new No_fout_relation(where);
    }
    No_fout_relation::~No_fout_relation(){}
    void No_fout_relation::print_settings()
    {
        std::cout << "Used constant fout parameter: " << param_obj.fout << std::endl;
    }

    /*Fout_Lodato_Rossi_eq28*/
    Fout_Lodato_Rossi_eq28::Fout_Lodato_Rossi_eq28(Parameters& _param_obj) : Fout_calculation(), param_obj(_param_obj){}
    void Fout_Lodato_Rossi_eq28::calculate_fout()
    {
        double fout, Mdot_rate, Mdotfb;
        Mdotfb = param_obj.Mdotfb_t;
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

    void Fout_Lodato_Rossi_eq28::print_settings()
    {
        std::cout << "Used time dependent fout" << std::endl;
    }

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

    void Mdotpeak_GuillochonRamirez::calc_mdotpeak()
    {
        if(param_obj.s.get_politrop() == "5per3")
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

    Mdotpeak_GuillochonRamirez::~Mdotpeak_GuillochonRamirez(){}

    Mdotpeak_L09_constrho::Mdotpeak_L09_constrho(Parameters& _param_obj) : param_obj(_param_obj){}
    void Mdotpeak_L09_constrho::calc_mdotpeak()
    {
        double mdotpeak;
        mdotpeak = 1./2.*param_obj.s.get_mass()*M_SUN/(param_obj.tmin*24.*60.*60.);
        param_obj.Mdot_peak = mdotpeak;
    }
    Mdot_peak_calculation* Mdotpeak_L09_constrho::clone(Parameters& where)
    {
        return new Mdotpeak_L09_constrho(where);
    }
    Mdotpeak_L09_constrho::~Mdotpeak_L09_constrho(){}

    Mdotpeak_L09_all::Mdotpeak_L09_all(Parameters& _param_obj) : param_obj(_param_obj){}
    Mdotpeak_L09_all::~Mdotpeak_L09_all(){}
    Mdot_peak_calculation* Mdotpeak_L09_all::clone(Parameters& where)
    {
        return new Mdotpeak_L09_all(where);
    }
    void Mdotpeak_L09_all::calc_mdotpeak()
    {
        double tpeak;
        if("5per3" == param_obj.s.get_politrop())
        {
            tpeak = tpeak_relative_to_tmin_n3per2 * param_obj.tmin;
        }
        else
        {
            tpeak = tpeak_relative_to_tmin_n3 * param_obj.tmin;
        }
        param_obj.calcMdotfb ->calc_Mdotfb_at_t(tpeak);
        param_obj.Mdot_peak = param_obj.Mdotfb_t;
    }

    /*tmin calculation classes, its derivated classes functions*/

        tmin_calculation::tmin_calculation(Parameters& _param_obj) : param_obj(_param_obj){}
        tmin_calculation::~tmin_calculation(){}

        tmin_Lodato::tmin_Lodato(Parameters& _param_obj) : tmin_calculation(_param_obj){}
        tmin_Lodato::~tmin_Lodato(){}
        tmin_calculation* tmin_Lodato::clone(Parameters& where) 
        {
            return new tmin_Lodato(where);
        }

        void tmin_Lodato::calc_tmin()
        {
            double eq1,eq2;
            eq1 = (param_obj.rp/(param_obj.s.get_radius()*R_SUN_METER));
            eq2 = std::pow(param_obj.rp,3) / (G * param_obj.bh.get_mass()*1e6*M_SUN);
            param_obj.tmin = M_PI / std::sqrt(2) * std::pow(eq1,1.5) * std::sqrt(eq2) * s_to_day;
        }

        void tmin_Lodato::print_setted_type()
        {
            std::cout << "Used tmin calculation with rp distance. tmin = " << param_obj.tmin << " day" << std::endl;
        }

        tmin_with_rt::tmin_with_rt(Parameters& _param_obj) : tmin_calculation(_param_obj){}
        tmin_with_rt::~tmin_with_rt(){}
        tmin_calculation* tmin_with_rt::clone(Parameters& where)
        {
            return new tmin_with_rt(where);
        }
        void tmin_with_rt::calc_tmin()
        {
            double eq1,eq2;
            eq1 = (param_obj.rt/(param_obj.s.get_radius()*R_SUN_METER));
            eq2 = std::pow(param_obj.rt,3) / (G * param_obj.bh.get_mass()*1e6*M_SUN);
            param_obj.tmin = M_PI / std::sqrt(2) * std::pow(eq1,1.5) * std::sqrt(eq2) * s_to_day;
        }

        void tmin_with_rt::print_setted_type()
        {
            std::cout << "Used classical tmin calculation with rt distance. tmin = " << param_obj.tmin << " day" << std::endl;
        }

        tmin_from_GR13::tmin_from_GR13(Parameters& _param_obj) : tmin_calculation(_param_obj) {}
        tmin_from_GR13::~tmin_from_GR13(){}
        tmin_calculation* tmin_from_GR13::clone(Parameters& where)
        {
            return new tmin_from_GR13(where);
        }


        double tmin_from_GR13::B4per3()
        {
            double b = param_obj.beta;
            double b4per3;
            b4per3 = (-0.38670 + 0.57291 * std::sqrt(b) - 0.31231 * b)/ (1.0 - 1.2744 * std::sqrt(b) - 0.90053 * b);
            return b4per3;
        }

        double tmin_from_GR13::B5per3()
        {
            double b = param_obj.beta;
            double b5per3;
            b5per3 = (-0.30908 + 1.1804*std::sqrt(b) - 1.1764*b)/(1. + 1.3089*std::sqrt(b) - 4.1940*b);
            return b5per3;
        }

        void tmin_from_GR13::calc_tmin()
        {
            double B, tp_rt_tmin;
            if ("4per3" == param_obj.s.get_politrop())
            {
                B = B4per3();
                tp_rt_tmin = tpeak_relative_to_tmin_n3;
            }
            else
            {
                B = B5per3();
                tp_rt_tmin = tpeak_relative_to_tmin_n3per2;
            }
            param_obj.tmin = B / tp_rt_tmin * std::pow(param_obj.bh.get_mass(),0.5) * std::pow(param_obj.s.get_mass(), -1) * std::pow(param_obj.s.get_radius(),1.5) * YEAR_TO_DAY;
        }

        void tmin_from_GR13::print_setted_type()
        {
            std::cout << "Used tmin calculation from GR13 tpeak value. tmin = " << param_obj.tmin << " day" << std::endl;
        }


    n_inf_calculation::n_inf_calculation(Parameters& _param_obj) : param_obj(_param_obj) {}
    n_inf_calculation::~n_inf_calculation(){}

        n_inf_const::n_inf_const(Parameters& _param_obj) : n_inf_calculation(_param_obj) {}
        n_inf_const::n_inf_const(Parameters& _param_obj, double _ninf) : n_inf_calculation(_param_obj), ninf(_ninf) {}
        n_inf_const::~n_inf_const(){}
        n_inf_calculation* n_inf_const::clone(Parameters& where)
        {
            return new n_inf_const(where, ninf);
        }
        void n_inf_const::set_ninf(double _ninf)
        {
            ninf =_ninf;
        }
        double n_inf_const::get_ninf()
        {
            return ninf;
        }

        n_inf_GR13::n_inf_GR13(Parameters& _param_obj) : n_inf_calculation(_param_obj) {}
        n_inf_GR13::~n_inf_GR13(){}
        n_inf_calculation* n_inf_GR13::clone(Parameters& where)
        {
            return new n_inf_const(where);
        }
        double n_inf_GR13::D_4per3()
        {
            double D;
            double b = param_obj.beta;
            D = (-2.7332 + 6.9465*b - 3.2743*b*b - 0.84659*std::pow(b,3) + 0.56254*std::pow(b,4))/(1.0 - 2.3585*b + 0.47593*b*b + 0.96280*std::pow(b,3) - 0.37996*std::pow(b,4));
            return D;
        }
        double n_inf_GR13::D_5per3()
        {
            double D;
            double b = param_obj.beta;
            D = (-0.93653 + 11.109*b - 38.161*b*b + 50.418*std::pow(b,3) - 22.965*std::pow(b,4))/(1.0 - 8.6394*b + 26.012*b*b - 32.383*std::pow(b,3) + 14.350*std::pow(b,4));
            return D;
        }

        double n_inf_GR13::get_ninf()
        {
            if("4per3" == param_obj.s.get_politrop())
            {
                return D_4per3();
            }
            else
            {
                return D_5per3();
            }
        }



    Mdotfb_calculation::Mdotfb_calculation(Parameters& _param_obj) : param_obj(_param_obj) 
    {
        calc_ninf=std::unique_ptr<n_inf_calculation>(new n_inf_const(param_obj, -5./3.));
        n_inf = calc_ninf -> get_ninf();
    }
    Mdotfb_calculation::~Mdotfb_calculation(){}
    void Mdotfb_calculation::set_calc_ninf(std::string ninf)
    {
        if(IsDouble(ninf))
        {
            calc_ninf.reset(new n_inf_const(param_obj, stod(ninf)));
        }
        else if("GR13" == ninf)
        {
            calc_ninf.reset(new n_inf_GR13(param_obj));
        }
        else
        {
            throw std::runtime_error("ERROR! Unknown ninf calculation: " + ninf);
        }
        n_inf = calc_ninf -> get_ninf();
    }

        Mdotfb_classical_powerfunction::Mdotfb_classical_powerfunction(Parameters& _param_obj) : Mdotfb_calculation(_param_obj) {}
        Mdotfb_classical_powerfunction::~Mdotfb_classical_powerfunction(){}

        Mdotfb_calculation* Mdotfb_classical_powerfunction::clone(Parameters& where)
        {
            return new Mdotfb_classical_powerfunction(where);
        }

        void Mdotfb_classical_powerfunction::calc_Mdotfb_at_t(double t)
        {
            param_obj.Mdotfb_t = param_obj.Mdot_peak * std::pow(t/param_obj.tmin, n_inf);
        }

        void Mdotfb_classical_powerfunction::print_setted_type()
        {
            std::cout << "Mdotfb calculated as classical power function." << std::endl;
        }

        Mdotfb_L09_constrho_x::Mdotfb_L09_constrho_x(Parameters& _param_obj) : Mdotfb_calculation(_param_obj){}
        Mdotfb_L09_constrho_x::~Mdotfb_L09_constrho_x(){}

        Mdotfb_calculation* Mdotfb_L09_constrho_x::clone(Parameters& where)
        {
            return new Mdotfb_L09_constrho_x(where);
        }

        void Mdotfb_L09_constrho_x::calc_Mdotfb_at_t(double t)
        {
            param_obj.Mdotfb_t = param_obj.Mdot_peak * (1.0 - std::pow(t/param_obj.tmin, -4./3.)) * std::pow(t/param_obj.tmin, n_inf);
        }

        void Mdotfb_L09_constrho_x::print_setted_type()
        {
            std::cout << "Mdotfb calculated as L09 with constant rho" << std::endl;
        }

        Mdotfb_L09_all::Mdotfb_L09_all(Parameters& _param_obj) : Mdotfb_calculation(_param_obj) {}
        Mdotfb_L09_all::~Mdotfb_L09_all(){}

        Mdotfb_calculation* Mdotfb_L09_all::clone(Parameters& where)
        {
            return new Mdotfb_L09_all(where);
        }

        void Mdotfb_L09_all::calc_constpart()
        {
            const_part = param_obj.s.get_radius() * R_SUN_METER * param_obj.s.get_K() * (param_obj.s.get_n() + 1.) * std::pow(param_obj.s.get_rho_c(), 1.0/param_obj.s.get_n()) / (3.0 * G * param_obj.tmin * 24. * 60. * 60.);
        }

        void Mdotfb_L09_all::init()
        {
            calc_constpart();
            std::string needed_file;
            if(const char* path = std::getenv("TIDE_PATH"))
            {
                needed_file = path;
                needed_file += "/share/";
            }
            else
            {
                throw std::runtime_error("ERROR: TIDE_PATH environment variable not found");
            }
            if("4per3" == param_obj.s.get_politrop())
            {
                needed_file += "table_n3.dat";
                interpolate.init(needed_file);
            }
            else if("5per3" == param_obj.s.get_politrop())
            {
                //std::cout << "Opened file: star_structure/table_n3per2.dat" << std::endl;
                needed_file += "table_n3per2.dat";
                interpolate.init(needed_file);
            }
            interpolate.set_tmin(param_obj.tmin);
        }

        //kg/s kell legyen
        void Mdotfb_L09_all::calc_Mdotfb_at_t(double t)
        {
            param_obj.Mdotfb_t = const_part * std::pow(t/param_obj.tmin, -5./3.) * interpolate.interpolateAtX(t);
        }

        void Mdotfb_L09_all::print_setted_type()
        {
            std::cout << "Mdotfb calculated as L09." << std::endl;
        }