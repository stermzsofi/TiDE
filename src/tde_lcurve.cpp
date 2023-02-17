#include "tde_lcurve.hpp"




/*
    *********

        Here will start the different part of the light curve calculation: Wind and disk part

    *********
*/

/*
    *********Wind part of the light curve functions*********
*/

    wind_calculation_radius_temperature::wind_calculation_radius_temperature(Parameters& _param) : param(_param){}
    wind_calculation_radius_temperature::~wind_calculation_radius_temperature(){}
    void wind_calculation_radius_temperature::calc_rphconst()
    {
        rph_const_part = KAPPA_ES/(4.0*M_PI*param.vw);
    } 

        /*rph old*/
        rph_old::rph_old(Parameters& _param) : wind_calculation_radius_temperature(_param)
        {
            calc_rphconst();
        }
        rph_old::~rph_old(){}
        wind_calculation_radius_temperature* rph_old::clone(Parameters& where)
        {
            return new rph_old(where);
        }

        double rph_old::calc_radius()
        {
            double rph;
            rph = rph_const_part * param.fout * param.Mdotfb_t;
            return rph;
        }

        double rph_old::calc_temperature(double r)
        {
            double tph;
            /*ronda de próbáljuk meg*/
            /*double tau, Labs;
            tau = KAPPA_GAMMA * param.fout * param.Mdotfb_t /(4.0*M_PI*param.vw * param.rL);
            Labs = param.Ldiskbol * (1.0 - std::exp(-tau));*/
            tph = param.TL * std::pow(r/param.rL, -2./3.) * std::pow(param.fout/param.fv, 1./3.);
            return tph;
        }

        /*rtr*/
        rtr::rtr(Parameters& _param) : wind_calculation_radius_temperature(_param)
        {
            calc_rphconst();
        }
        rtr::~rtr(){}
        wind_calculation_radius_temperature* rtr::clone(Parameters& where)
        {
            return new rtr(where);
        }

        double rtr::calc_radius()
        {
            double r_tr;
            r_tr = rph_const_part * param.vw * param.fout * param.Mdotfb_t/(c);
            return r_tr;
        }

        double rtr::calc_temperature(double r)
        {
            double tph;
            tph = param.TL * std::pow(r/param.rL, -2./3.) * std::pow(param.fout/param.fv, 1./3.);
            return tph;
        }

        /*rtr_full*/
        rtr_full::rtr_full(Parameters& _param) : wind_calculation_radius_temperature(_param)
        {
            calc_rphconst();
        }
        rtr_full::~rtr_full(){}
        wind_calculation_radius_temperature* rtr_full::clone(Parameters& where)
        {
            return new rtr_full(where);
        }

        double rtr_full::calc_radius()
        {
            double A, rw, a, b, c_var, r_tr;
            A = rph_const_part * param.vw * param.fout * param.Mdotfb_t/(param.rL * c);
            rw = param.rL + (param.time - param.tmin) * DAY_TO_S * param.vw;
            a = A * param.rL / (rw*rw);
            b = -1.0 * (1.0 + 2.0* param.rL * A/rw);
            c_var = param.rL * (1.0 + A);
            r_tr = (-b - std::sqrt(b*b-4.0*a*c_var))/(2.0*a);
            return r_tr;
        }

        double rtr_full::calc_temperature(double r)
        {
            double tph;
            tph = param.TL * std::pow(r/param.rL, -2./3.) * std::pow(param.fout/param.fv, 1./3.);
            return tph;
        }

        /*rc_vs_rtr*/
        rc_vs_rtr::rc_vs_rtr(Parameters& _param) : wind_calculation_radius_temperature(_param)
        {
            calc_rphconst();
            calc_rc_const_part();
        }
        rc_vs_rtr::~rc_vs_rtr(){}
        wind_calculation_radius_temperature* rc_vs_rtr::clone(Parameters& where)
        {
            return new rc_vs_rtr(where);
        }

        void rc_vs_rtr::calc_rc_const_part()
        {
            rc_const_part = std::sqrt(3.0 * KAPPA_0 * KAPPA_ES) * std::pow(param.fv, 7./12.) * std::pow(param.rL, -7./6.) * std::pow(4.0*M_PI*param.vw, -1.5) * 0.5;
        }

        void rc_vs_rtr::calc_rtr()
        {
            double A, rw, a, b, c_var, r_tr;
            A = rph_const_part * param.vw * param.fout * param.Mdotfb_t/(param.rL * c);
            rw = param.rL + (param.time - param.tmin) * DAY_TO_S * param.vw;
            a = A * param.rL / (rw*rw);
            b = -1.0 * (1.0 + 2.0* param.rL * A/rw);
            c_var = param.rL * (1.0 + A);
            r_tr = (-b - std::sqrt(b*b-4.0*a*c_var))/(2.0*a);
            rtr = r_tr;
            /*double r_tr;
            r_tr = rph_const_part * param.vw * param.fout * param.Mdotfb_t/(c);
            rtr = r_tr;*/
        }

        void rc_vs_rtr::calc_rc()
        {
            calc_rtr();
            double A = rc_const_part * std::pow(param.fout, 11./12.) * std::pow(param.Mdotfb_t, 1.5) * std::pow(param.TL, -7./4.) * std::pow(rtr, 7./24.);
            rc = std::pow(A, 8./9.);
            //double randa_nagy_szam = std::pow(param.bh.get_mass()*1e6*M_SUN, 2./9.) * std::pow(param.Mdotfb_t, 53./54.) * std::pow(27,1./27.) /(std::pow(2.,2.888889) * std::pow(SIGMA,7./18.) * std::pow(KAPPA_ES, 8./27.) * std::pow(4.0 * M_PI, 53./54.)) * std::pow(KAPPA_0, 12./27.) * std::pow(G, 6./27.) * std::pow(param.fout, 2./27.);
            //randa_nagy_szam = std::pow(2.0,-53./27.) * std::pow(3.0, 4./9.) * std::pow(M_PI, -11./54.) * std::pow(SIGMA, 7./18.) * std::pow(KAPPA_0, 4./9.) * std::pow(KAPPA_ES, -8./27.) * std::pow(param.fout, 2./27.) * std::pow(param.Mdotfb_t, 11./54.) * std::pow(c, -23./27.) * std::pow(G * param.bh.get_mass()*1e6 * M_SUN, -7./27.);
            //std::cout << param.time << "\t" << rc << "\t" << rtr << "\t" << randa_nagy_szam << "\t" << randa_nagy_szam * std::pow(param.eta, -43./54.) << std::endl;
        }

        double rc_vs_rtr::calc_T_rtr()
        {
            double tph;
            tph = param.TL * std::pow(rtr/param.rL, -2./3.) * std::pow(param.fout/param.fv, 1./3.);
            return tph;
        }

        double rc_vs_rtr::calc_T_rc()
        {
            double tph, tc;
            tph = calc_T_rtr();
            tc = std::sqrt(rtr/rc) * tph;
            return tc;
        }

        double rc_vs_rtr::calc_radius()
        {
            calc_rc();
            if(rc > rtr)
            {
                std::cout << "Time " << param.time << " used rc" << std::endl;
                return rc;
            }
            else
            {
                return rtr;
            }
        }

        double rc_vs_rtr::calc_temperature([[maybe_unused]] double r)
        {
            if(rc > rtr)
            {
                //std::cout << "Used rc at time " << param.time << std::endl;
                return calc_T_rc();
            }
            else
            {
                return calc_T_rtr();
            }
        }




        /*Constructor*/
        Wind_part_of_lightcurve::Wind_part_of_lightcurve(Parameters& p) : par(p)
        {
            if(par.wind_radius == "rph")
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rph_old(par));
            }
            else if(par.wind_radius == "rtr")
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rtr(par));
            }
            else if(par.wind_radius == "rtr_full")
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rtr_full(par));
            }
            else
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rc_vs_rtr(par));
            }
        }
       

        /*Default constructor*/
        Wind_part_of_lightcurve::Wind_part_of_lightcurve() : par(*new Parameters)
        {
            if(par.wind_radius == "rph")
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rph_old(par));
            }
            else if(par.wind_radius == "rtr")
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rtr(par));
            }
            else if(par.wind_radius == "rtr_full")
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rtr_full(par));
            }
            else
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rc_vs_rtr(par));
            }
        }

        void Wind_part_of_lightcurve::init()
        {
            if(par.wind_radius == "rph")
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rph_old(par));
            }
            else if(par.wind_radius == "rtr")
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rtr(par));
            }
            else if(par.wind_radius == "rtr_full")
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rtr_full(par));
            }
            else
            {
                calc_radius_temp = std::unique_ptr<wind_calculation_radius_temperature>(new rc_vs_rtr(par));
            }
        }

        /*Other functions*/
        /*double Wind_part_of_lightcurve::calc_temperature_at_L()
        {
            double Tl4, Tl;
            Tl4 = G * par.bh.get_mass() * 1e6 * M_SUN * par.Mdotfb_t / (4. * M_PI * std::pow(par.rL, 3) * SIGMA);
            Tl = std::pow(Tl4, 1./4.);
            return Tl; 
        }*/

        /*double Wind_part_of_lightcurve::calc_temperature_ph()
        {
            double tphtime;
            temperature_at_L = calc_temperature_at_L();
            tphtime = temperature_at_L * std::pow(photospheric_rad/par.rL, -2./3.) * std::pow(par.fout/par.fv, 1./3.);
            return tphtime;
        }*/

        //in m
        /*double Wind_part_of_lightcurve::calc_radius_ph()
        {
            double rphtime;
            rphtime = KAPPA_ES/(4. * M_PI * par.vw) * par.fout * par.Mdotfb_t;
            photospheric_rad = rphtime;
            return rphtime;
        }*/

        //will return in erg unit
        double Wind_part_of_lightcurve::calc_luminosity_at_nu(double nu)
        {
            double lnu;
            //calc_radius_ph();
            //photospheric_temperature = calc_temperature_ph();
            wind_calculation_radius = calc_radius_temp ->calc_radius();
            wind_calculation_temperature = calc_radius_temp ->calc_temperature(wind_calculation_radius);
            //ide kell akkor beletenni Józsi féle reprocessinget
            //kell: tau=kappa*mstar/(4*pi*10*rT*rph)
            //Lrep = eps_rep * Ldiskbol*(1-exp(tau))
            //wind_calculation_temp = (wind_calculation_temp**4 + Lrep/(4*pi*rph**2*sigma))^(1/4)
            //std::cout << "wind temp before = " << wind_calculation_temperature << "\t";
            //double tau = KAPPA_ES * par.s.get_mass()*M_SUN/(40.0*M_PI*par.rt * wind_calculation_radius);
            double tau = KAPPA_ES * par.fout * par.Mdotfb_t/(4.0 * M_PI* par.rL * par.vw);
            double Lrep = par.eta_reprocessing * par.Ldisk_bolometric_at_t * (1.0 - std::exp(-tau));
            //std::cout << "Lrep = " << Lrep << "\t" << "1-exp(-tau) = " << (1.0 - std::exp(-tau)) << "\t" << "tau = " << tau << "\tLdiskbol = " << par.Ldisk_bolometric_at_t << "\t";
            wind_calculation_temperature = std::pow(std::pow(wind_calculation_temperature,4) + Lrep/(4.0*M_PI*wind_calculation_radius*wind_calculation_radius*SIGMA),0.25);
            //std::cout << "wind temp after = " << wind_calculation_temperature << std::endl;
            lnu = 4.*M_PI*M_PI*std::pow(wind_calculation_radius,2.)*Planck_function_for_nu(nu, wind_calculation_temperature);
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

        void calc_Disk_at_r::set_nu(double new_nu)
        {
            own_nu = new_nu;
        }
        
        //Retrurns kg/s unit
        double calc_Disk_at_r::calc_Mdot()
        {   
            mdot = (1. - par.fout)*par.Mdotfb_t;
            return mdot;
        }

        double calc_Disk_at_r::calc_Td(double r)
        {
            double Td, f, factor1, factor2, factor_help, Td4;
            calc_Mdot();
            f = 1. - std::pow(par.rin/r, 0.5);
            factor1 = 3.*G*par.bh.get_mass()*1.0e6*M_SUN*mdot*f/(8.*M_PI*r*r*r);
            factor_help = 1./4. + 3./2.*f*std::pow((mdot*par.bh.get_radius())/(par.eta*par.mdot_edd*r), 2);
	        factor2 = 1./2. + std::pow(factor_help, 0.5);
            Td4  = factor1/(factor2 * SIGMA);
            //std::cout  << par.fout << " " <<par.Mdotfb_t << " " <<"\nTd4 sigma not bol = " << Td4 * SIGMA;
            //std::cout << par.time << "\tmon\t" << mdot << "\t" << f << "\t" << factor1 << "\t" << factor_help << "\t" << factor2 << "\t" << Td4 << "\t" << r << std::endl; 
            Td = std::pow(Td4, 1./4.);
            if(Td > Tmax) 
            {   
                //std::cerr << "new tmax found when: " << Td << "> " << Tmax << "at r:" << r << std::endl;
                Tmax = Td;
                rmax_per_rin = r/par.rin;
            }
            return Td;
        }

        double calc_Disk_at_r::calc_Td_at_rhalf()
        {
            return calc_Td(par.rhalf);
        }

        double calc_Disk_at_r::calc_F_at_r(double r, double nu)
        {
            double F;
            //F = 2*M_PI*r*Planck_function_for_nu(nu, calc_Td(r));
            F = 4.0*M_PI*M_PI*r*Planck_function_for_nu(nu, calc_Td(r));
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
        double Disk_part_of_lightcurve::calc_disk_part_at_t()
        {
            double sum_disk = 0;
            /**ez lesz itt az új*/
            calc_disk_r.Tmax = 0.0;
            calc_disk_r.calc_Mdot();
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

        double Disk_part_of_lightcurve::calc_disk_part_at_t_nu(double nu)
        {
            calc_disk_r.set_nu(nu);
            return calc_disk_part_at_t();
        }

        void Disk_part_of_lightcurve::set_nu(double new_nu)
        {
            calc_disk_r.set_nu(new_nu);
        }

        double Disk_part_of_lightcurve::reprocessing_g()
        {
            double eps_g;
            double a_alap, b_alap, a_Tin,b_Tin,c_Tin,a_Tout,b_Tout,a_cover,b_cover,c_cover;
            double alap, Tin, Tout, cover;
            alap = std::log10(calc_disk_r.par.n)/std::log10(calc_disk_r.par.N_col) * std::log10(calc_disk_r.par.Rmax);
            Tin = std::log10(get_Tmax());
            Tout = std::log10(get_Tdisk_at_r(calc_disk_r.par.rout));
            cover = calc_disk_r.par.delta_omega/(4.0*M_PI);
            a_alap = -0.0004701;
            b_alap = 0.00372817;
            a_Tin = -0.00432228;
            b_Tin = 0.04835209;
            c_Tin = -0.13490034;
            a_Tout = -0.00052651;
            b_Tout = 0.00261462;
            a_cover = 1.23506820e+01;
            b_cover = -1.07035887e-01;
            c_cover = 1.71222168e-05;
            eps_g = a_alap*alap + b_alap + a_Tin * Tin*Tin + b_Tin * Tin + c_Tin + a_Tout * Tout + b_Tout + a_cover * cover*cover + b_cover*cover + c_cover;
            return eps_g * calc_disk_r.par.Ldisk_bolometric_at_t * 1e7; //1e7: J->erg
        }

        double Disk_part_of_lightcurve::get_Tdisk_at_rhalf()
        {
            calc_disk_r.calc_Mdot();
            return calc_disk_r.calc_Td_at_rhalf();
        }

        double Disk_part_of_lightcurve::get_Tdisk_at_r(double r)
        {
            calc_disk_r.calc_Mdot();
            return calc_disk_r.calc_Td(r);
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

    void Total_lightcurve::print_output(double time, std::ofstream& out)
    {
        double lwind, ldisk;
        par.calculate_timedependent_parameters(time);
        lwind = wind_part.calc_luminosity_at_nu(par.nu);
        ldisk = disk_part.calc_disk_part_at_t();
        out << time << "\t" << lwind + ldisk << "\t" << lwind << "\t" << ldisk << "\t" << disk_part.reprocessing_g() << std::endl;
    }

    void Total_lightcurve::calculate_lightcurve_timeintervall(double tmin, double tmax, double dt)
    {
        double t;
        std::ofstream lum;
	    lum.open("lcurve.dat");
        output_commentline(lum);
        disk_part.set_nu(par.nu);
        /*for (auto i = 0; i< 1000; i++)
        {
            double r = i*(par.rout - par.rin) /1000.0 + par.rin;
            std::cout << r << "\t" << disk_part.get_Tdisk_at_r(r) << std::endl;
        }*/
        if(par.extras)
        {
            std::ofstream extra;
            extra.open("extras.dat");
            par.extras_commentline(extra);
            for(t = tmin; t <= tmax; t+= dt)
            {
                print_output(t,lum);
                extra << t << "\t" << wind_part.get_phot_rad() << "\t" << wind_part.get_phot_temp() <<
                 "\t" << par.TL << "\t" << /*disk_part.get_Tdisk_at_rhalf(t)*/disk_part.get_Tdisk_at_rhalf() <<
                  "\t" << par.fout << "\t" << /*par.tedge << "\t" <<*/ par.Mdotfb_t << "\t" << par.Ldisk_bolometric_at_t << 
                  "\t" << disk_part.get_Tmax() << "\t" << disk_part.get_Tdisk_at_r(par.rout) << "\t" << disk_part.get_rmax() << 
                  "\t" << par.N_col << "\t" << par.n << "\t" << par.Rmax << std::endl;
            }
            extra.close();
        }
        else
        {
            for(t = tmin; t <= tmax; t += dt)
            {
                print_output(t, lum);
            }
        }
        lum.close();
    }


/*Diffusion calculations*/
    /*Fx derivation class: Lum_diffusion_fx functions*/
        Lum_diffusion_fx::Lum_diffusion_fx(Parameters& p) : par(p), disk_part(p), wind_part(p)
        {
            overflow_corrected = false;
        }

        Lum_diffusion_fx::Lum_diffusion_fx() : par(*new Parameters), disk_part(par), wind_part(par)
        {
            overflow_corrected = false;
        }

        void Lum_diffusion_fx::set_tend(double _t_end)
        {
            t_end = _t_end;
        }

        double Lum_diffusion_fx::L_inp(double t_offset)
        {
            double Linp;
            Linp = par.eta * par.Mdot_peak * std::pow(1 + t_offset/par.tmin, -5./3.) * c*c; //kg*m2/s3
            Linp *= 1e7; // -> to erg/s
            return Linp;
        }

        double Lum_diffusion_fx::L_inp_new(double t)
        {
            double l_wind, l_disk, l_reprocess_g, Linp, t_offset;
            t_offset = t + par.tmin;
            par.calculate_timedependent_parameters(t_offset);   
            l_wind = wind_part.calc_luminosity_at_nu(par.nu);
            l_disk = disk_part.calc_disk_part_at_t();
            if(par.my_rep /*&& par.nu > 6.2e14 && par.nu < 6.4e14*/) l_reprocess_g = disk_part.reprocessing_g()/6.3e14;
            else l_reprocess_g = 0.0;
            Linp = l_disk + l_wind + l_reprocess_g;
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
            diffusion_new = Quad_Trapezoidal(diff_fx, 0., time_offset, 0.01);
            previous_integral = 0.;
            previous_time = 1e-5;
            overflow = false;
        }

        Diffusion_luminosity::Diffusion_luminosity() : par(*new Parameters), diff_fx(par)
        {
            time = par.tmin;
            time_offset = time - par.tmin;
            diffusion = Trapezoidal_rule(diff_fx, 0., time_offset, 0.01, 1);
            diffusion_new = Quad_Trapezoidal(diff_fx, 0., time_offset, 0.7);
            previous_integral = 0.;
            previous_time = 1e-5;
            overflow = false;
        }

        void Diffusion_luminosity::set_time(double t)
        {
            time = t;
            time_offset = time - par.tmin;
            diffusion.set_xend(time_offset);
            diffusion_new.set_xend(time_offset);
            diffusion.set_xstart(previous_time);
            diffusion_new.set_xstart(previous_time);
            diff_fx.set_tend(time_offset);
        }

        double Diffusion_luminosity::calc_Lum_at_t()
        {
            double res, integral;
            if(!overflow)
            {
                try
                {
                    integral = diffusion_new.qtrap() + previous_integral;
                    previous_integral = integral;
                    integral = std::log(integral);
                    res = (-time_offset/par.diffusion_timescale) + integral;
                    res = std::exp(res) / par.diffusion_timescale;
                    previous_time = time_offset;
                }
                catch(const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    overflow = true;
                    res = diff_fx.L_inp_new(time_offset);
                }
                catch(char const* e)
                {
                    std::cerr << e << '\n';
                    overflow = true;
                    res = diff_fx.L_inp_new(time_offset);
                }
            }
            else
            {
                res = diff_fx.L_inp_new(time_offset);
            }
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

        void Diffusion_luminosity::reset()
        {
            previous_integral = 0.;
            previous_time = 1e-5;
            overflow = false;
        }

        void Diffusion_luminosity::set_prevtime_previntegral(double pt, double pi)
        {
            previous_time = pt;
            previous_integral = pi;
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
                par.calculate_timedependent_parameters(t);
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
            Lwind = wind_part.calc_luminosity_at_nu(par.nu);
            Ldisk = disk_part.calc_disk_part_at_t();
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
            lwind = wind_part.calc_luminosity_at_nu(par.nu);
            ldisk = disk_part.calc_disk_part_at_t();
            spectra << nu << "\t" << frequency_to_wavelength(nu) << "\t" << lwind << "\t" << ldisk << std::endl;
        }
    }
