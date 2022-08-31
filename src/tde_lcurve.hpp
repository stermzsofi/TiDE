#ifndef TDE_LCURVE
#define TDE_LCURVE

#include "tde_lcurve_assistant.hpp"
#include "Trapezoidal_rule/trapezoidal.hpp"
#include "star_structure/linear_interpol.hpp"
#include "tde_parameters.hpp"
#include <memory>


struct Parameters;

/* 
    *********

    Here will start the different part of the light curve calculation: Wind and disk part

    Wind_part_of_light curve
        calculate the wind part of the light curve
        need a parameters struct (par)
        calc_luminosity_at_nu: calculate the luminosity at t time and nu ferquency [erg/s/Hz]

    Disk part of light curve
        We have calculate the luminosity at different r radius of the disk than integral it
        For the integration: we need a class derivated from Fx
        - calc_Disk_at_r: Fx derivated class:
                calc_F_at_r function: luminosity at r radius and nu frequency
        - Disk_part_of_lightcurve: calculate the disk part of the light curve
                calc_disk_part_at_t_nu function: calculate the total luminosity of the disk
                at t time and nu freqency [erg/s/Hz]

    ********* 
*/

/*Wind part of the light curve*/
        class Wind_part_of_lightcurve
        {
            public:
                Wind_part_of_lightcurve(Parameters& p);
                Wind_part_of_lightcurve();
                double calc_temperature_at_L();
                double calc_temperature_ph();
                double calc_radius_ph();
                double calc_luminosity_at_nu(double nu);

                /*Get photospheric temperature and radius*/
                double get_phot_temp(){return photospheric_temperature;}
                double get_phot_rad(){return photospheric_rad;}

                Parameters& par;

            private:
                double photospheric_temperature;
                double photospheric_rad;
                double temperature_at_L;
                double const_multiplier_at_rph;
        };

/*Calculate disk part of the lightcurve:
    wee need two class
    first will be derived from Fx (definied at trapezoidal.hpp, we will call it calc_disk_at_r)
    second will use the first class and a trapezoidal_rule class (also definied at trapezoidal.hpp)*/

    class calc_Disk_at_r : public Fx
    {
        public:
            calc_Disk_at_r(Parameters& p);
            calc_Disk_at_r();

            void set_nu(double new_nu);
            double calc_Mdot();     //returns kg/s unit
            double calc_Td(double r);
            double calc_Td_at_rhalf();
            double calc_F_at_r(double r, double nu);

            double calc_Fx(double x);

            Parameters& par;
        
        private:
            double mdot;
            double own_nu;
    };

    class Disk_part_of_lightcurve
    {
        public:
            Disk_part_of_lightcurve(Parameters& p);
            Disk_part_of_lightcurve();
            double calc_disk_part_at_t();
            double calc_disk_part_at_t_nu(double nu);
            void set_nu(double new_nu);

            double get_Tdisk_at_rhalf();


        private:
            calc_Disk_at_r calc_disk_r;
            Trapezoidal_rule sum_lum_disk;
    };



/*
    ********
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
    
    ********
*/

/*Calculate the total lightcurve of the event*/
    /*It will use a Wind_part_of_lightcurve class and a Disk_part_of_lightcurve class*/

    class Total_lightcurve
    {
        public:
            Total_lightcurve(Parameters& p);
            void output_commentline(std::ofstream & out);
            void print_output(double time, std::ofstream& out);
            void calculate_lightcurve_timeintervall(double tmin, double tmax, double dt);
        private:
            Parameters& par;
            Wind_part_of_lightcurve wind_part;
            Disk_part_of_lightcurve disk_part;
    };

/*Diffusion at early times*/
    /*Again two classes:
    Lum_diffusion_fx: class derived from Fx, use at trapezoidal.hpp, will call it at Lum_diffusion class
    */

    class Lum_diffusion_fx : public Fx
    {
        public:
            Lum_diffusion_fx(Parameters& p);
            Lum_diffusion_fx();
            void set_tend(double _t_end);
            void overflow_on(){overflow_corrected = true;}
            bool get_overflow(){return overflow_corrected;}
            double calc_Lt_integral(double t);
            double calc_Fx(double x);
            double L_inp(double t);
            double L_inp_new(double t);
        private:
            
            Parameters &par;
            Disk_part_of_lightcurve disk_part;
            Wind_part_of_lightcurve wind_part;
            bool overflow_corrected = false;
            double t_end;
    };

    class Lum_diffusion_fx_after_rl_limit : public Fx
    {
        public:
            Lum_diffusion_fx_after_rl_limit(Parameters& p);
            Lum_diffusion_fx_after_rl_limit();
            double calc_Lt_integral(double t);
            double calc_Fx(double x);
            double L_inp(double t);
            double L_inp_new(double t);
        private:
            
            Parameters &par;
            Disk_part_of_lightcurve disk_part;
            Wind_part_of_lightcurve wind_part;
    };

    class Diffusion_luminosity
    {
        public:
            Diffusion_luminosity(Parameters& p);
            Diffusion_luminosity();
            void set_time(double t);
            double calc_Lum_at_t();
            double get_Linp(double t);
            double get_integral(double t);
        private:
            Parameters& par;
            Trapezoidal_rule diffusion;
            Quad_Trapezoidal diffusion_new;
            Lum_diffusion_fx diff_fx;
            double time;
            double time_offset;
            double previous_time;
            double previous_integral;
            bool overflow;
    };

    class Diffusion_luminosity_timeintervall
    {
        public:
            Diffusion_luminosity_timeintervall(Parameters& p);
            Diffusion_luminosity_timeintervall();
            void calc_diffusion_luminosity_timeintervall(double tmin, double tmax, double dt);
        private:
            Parameters& par;
            Diffusion_luminosity lum_diff;
    };


/*Calculate the bolometric luminosity*/
    /*We will again need two classes:
    first will be derived from Fx (definied at trapezoidal.hpp, we will call it Lsum)
    second will use the first class and a trapezoidal_rule class (also definied at trapezoidal.hpp)*/

    class Lsum : public Fx
    {
        public:
            Lsum(Parameters& p, double t);
            Lsum(Parameters& p);
            Lsum();

            void set_time(double t);
            void set_nu(double nu);
            void comment(std::ofstream &);
            double calc_Lsum(double nu);
            double calc_Fx(double x);
        private:
            Parameters& par;
            Wind_part_of_lightcurve wind_part;
            Disk_part_of_lightcurve disk_part;
            double time;

    };

    class Lbol
    {
        public:
            Lbol(Parameters& p);

            double calc_Lbol(double time);
            void comment_line(std::ofstream & out);
            void Lbol_time(double tstart, double tend, double dt);
        private:
            Parameters& par;
            Trapezoidal_rule bolometric;
            Lsum Ls;
    };

/*Calculate the spectrum */
    /*It will use a Wind_part_of_lightcurve class and a Disk_part_of_lightcurve class*/
    class Spectrum
    {
        public:
            Spectrum(Parameters& p);

            void output_commentline(std::ofstream & out);
            void calculate_spectrum_nuintervall(double nu_min, double nu_max, double dnu);
        private:
            Parameters& par;
            Wind_part_of_lightcurve wind_part;
            Disk_part_of_lightcurve disk_part;
    };

 #endif




    