#include "tde_lcurve.hpp"





/*
    *********

        Here will start the different part of the light curve calculation: Wind and disk part

    *********
*/

/*
    *********Wind part of the light curve functions*********
*/

        /*Constructor*/
        Wind_part_of_lightcurve::Wind_part_of_lightcurve(Parameters& p) : par(p) 
        {
            const_multiplier_at_rph = KAPPA_ES/(4. * M_PI * par.vw);
        }
       

        /*Default constructor*/
        Wind_part_of_lightcurve::Wind_part_of_lightcurve() : par(*new Parameters)
        {
            const_multiplier_at_rph = KAPPA_ES/(4. * M_PI * par.vw);
        }

        /*Other functions*/
        double Wind_part_of_lightcurve::calc_temperature_at_L()
        {
            double Tl4, Tl;
            Tl4 = G * par.bh.get_mass() * 1e6 * M_SUN * par.Mdotfb_t / (4. * M_PI * std::pow(par.rL, 3) * SIGMA);
            Tl = std::pow(Tl4, 1./4.);
            return Tl; 
        }

        double Wind_part_of_lightcurve::calc_temperature_ph()
        {
            double tphtime;
            temperature_at_L = calc_temperature_at_L();
            tphtime = temperature_at_L * std::pow(photospheric_rad/par.rL, -2./3.) * std::pow(par.fout/par.fv, 1./3.);
            return tphtime;
        }

        //in m
        double Wind_part_of_lightcurve::calc_radius_ph()
        {
            double rphtime;
            rphtime = const_multiplier_at_rph * par.fout * par.Mdotfb_t;
            photospheric_rad = rphtime;
            return rphtime;
        }

        //will return in erg unit
        double Wind_part_of_lightcurve::calc_luminosity_at_nu(double nu)
        {
            double lnu;
            calc_radius_ph();
            photospheric_temperature = calc_temperature_ph();
            lnu = 4.*M_PI*M_PI*std::pow(photospheric_rad,2.)*Planck_function_for_nu(nu, photospheric_temperature);
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
        double Disk_part_of_lightcurve::calc_disk_part_at_t()
        {
            double sum_disk = 0;
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

        double Disk_part_of_lightcurve::get_Tdisk_at_rhalf()
        {
            calc_disk_r.calc_Mdot();
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

    void Total_lightcurve::print_output(double time, std::ofstream& out)
    {
        double lwind, ldisk;
        par.calculate_timedependent_parameters(time);
        lwind = wind_part.calc_luminosity_at_nu(par.nu);
        ldisk = disk_part.calc_disk_part_at_t();
        out << time << "\t" << lwind + ldisk << "\t" << lwind << "\t" << ldisk << std::endl;
    }

    void Total_lightcurve::calculate_lightcurve_timeintervall(double tmin, double tmax, double dt)
    {
        double t;
        std::ofstream lum;
	    lum.open("lcurve.dat");
        output_commentline(lum);
        disk_part.set_nu(par.nu);
        if(par.extras)
        {
            std::ofstream extra;
            extra.open("extras.dat");
            par.extras_commentline(extra);
            for(t = tmin; t <= tmax; t+= dt)
            {
                print_output(t,lum);
                extra << t << "\t" << wind_part.get_phot_rad() << "\t" << wind_part.get_phot_temp() << "\t" << wind_part.calc_temperature_at_L() << "\t" << /*disk_part.get_Tdisk_at_rhalf(t)*/disk_part.get_Tdisk_at_rhalf() << "\t" << par.fout << "\t" << /*par.tedge << "\t" <<*/ par.Mdotfb_t << std::endl;
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
            double l_wind, l_disk, Linp, t_offset;
            t_offset = t + par.tmin;
            par.calculate_timedependent_parameters(t_offset);   
            l_wind = wind_part.calc_luminosity_at_nu(par.nu);
            l_disk = disk_part.calc_disk_part_at_t();
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
