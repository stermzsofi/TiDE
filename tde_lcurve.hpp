#ifndef TDE_LCURVE
#define TDE_LCURVE

#include "tde_lcurve_assistant.hpp"
#include "Trapezoidal_rule/trapezoidal.hpp"
#include <memory>

enum startype : int;

/*

    **********
        Here start the needed object: 
            -Object mother class
                -Blackhole derivated class
                -Star derivated class
            
            -Mass_to_radii mother class
                -Main_sequence_star derivated class
                -White_dwarf_star derivated class
    **********
*/

/*Object class: derivated will be Blackhole and Star*/
class Object
{
    public:
        double get_mass(){return mass;}
        double get_radius(){return radius;}
        void set_mass(double new_mass) {mass = new_mass;}
        void set_radius(double new_radius) {radius = new_radius;}
    protected:
        double mass;
        double radius;
};

    //Blackhole class: its radius will be the Schwarzschild radius
    class Blackhole : public Object
    {
        public:
            Blackhole(double m6);
            Blackhole();
            void init();
        private:
            double calc_schwarzschild_radius();
    };

    class Mass_to_Radii;
    
    //Star class
    class Star : public Object
    {
        public:
            Star(double m, double r, startype relation);
            Star();
            void set_startype(std::string rad);
            startype get_startype(){return st;}
            void print_startype();
            void init();
        private:
            Mass_to_Radii* m_t_r;
            startype st;
    };

/*Mass_to_Radii class: relationship between mass and radius*/
class Mass_to_Radii
{
    public:
        virtual double calc_radius(double) = 0;//mass-> radius calculate (will be definied at derivated classes)
};

    /*Main seqence star relationship*/
    class Main_sequence_star : public Mass_to_Radii
    {
        public:
            double calc_radius(double mass);
    };

    /*White dwarf star relationship*/
    class White_dwarf_star : public Mass_to_Radii
    {
        public:
            double calc_radius(double mass);
    };

class Fout_calculation;
class Mdot_peak_calculation;


/*
    *******

    Parameters struct
        This is the main struct for calculations
        Contain all settable parameters 
        and all derivated constants which needed and during lightcurve calculations
        The settings of the to be calculation light curves also there

    Fout_calculation class: mother class
        To set the fout parameter:
            independent or
            time dependent parameter, calculated from Mdotfb and Mdoteddington
        
        - No_fout_relation: derivated class, independent fout
        - Fout_Lodato_Rossi_eq28: derivated class, time dependent fout

    Mdot_peak_calculation: mother class
        To calculate Mdot at peak
        
        - Mdotpeak_LodatoRossi: derivated class, calculate Mdot peak form 
                                Lodato and Rossi 2011 paper eq 6
                                This is the default set
        - Mdotpeak_GuillochonRamirez: deruvated class, calculate Mdot peak from
                                      Guillochon and Ramirez paper
                                      There is 2 equation: 
                                         4/3 politrop index
                                         5/3 politrop index

    ********
*/


/*Other needed parameters, some simple calculations*/
struct Parameters
{
        Parameters(int N_o, double e, double b, double fo, double fv_, double dist, double ink, double nu_, Star star, Blackhole blackhole);
        Parameters();
        ~Parameters();
        Parameters(Parameters& copied);
        
        //Functions used in constructor
        void fv_check();  //check using fv, v_wind may is higher, than light velocity
        void calculate_rphtmin();
        void calculate_tphtmin();
        void calculate_tmin();
        void calculate_rin();
        void calculate_rout();
		void calculate_rhalf();
        void calculate_rL();
        void calculate_mdot_edd();  //in kg/s
        void calculate_rpper3rs();
        void calculate_tedge();
        void calculate_vwind();
        void calculate_time_wind_limit();
        void calculate_time_rph_limit();
    
        /*Print parametersfunction*/
        void print_parameters();

        /*commentline functions*/
        void total_lc_commentline(std::ofstream& out);
        void lbol_commentline(std::ofstream& out);
        void extras_commentline(std::ofstream& out);
        void spectra_commentline(std::ofstream& out);
        
        //Functions to set fout (build in fout or independent) and Mdot peak calculation (Lodato and Rossi or Guillochon)
        void change_foutcalculation(Fout_calculation* new_fout);
        void change_mdotpeakcalculation(Mdot_peak_calculation* new_mdotpeak);

        /*Parameters*/
        double eta;
        double beta;
        double fout;
        double fv;
        double d;
        double i;
        double diffusion_timescale;
        int N;
        /*needed classes*/
        Star s;
        Blackhole bh;
        std::unique_ptr<Fout_calculation> calcfout;
        std::unique_ptr<Mdot_peak_calculation> calcmdotpeak;

        /*Lightcurve intervalls and step*/
        double t_start;
        double t_end;
        double dt;
        double nu;

        /*Spectra intervalls and step*/
        double nu_start;
        double nu_end;
        double dnu;
        double time;
        
        /*Derived quantities*/
        double tmin;
        double rphtmin;
        double mdot_edd;
        double tphtmin;
        double rin;
        double rout;
        double rhalf;
        double rL;
        //time limit until wind will exist
        double time_wind_limit;
        //time limit while rph = rl
        double time_rph_limit;
        double oneper4pid2;
        double cosi;
        double Mdot_peak;
        double rpper3rs;
        double tedge;
        double vwind;

        /*Switch for bolometric lightcurve and extra informations*/
        bool bolometric;
        bool extras;
        bool spectra;
        bool lcurve;
        bool diffusion;
        bool use_rph_limit;

        /*Bool variable for t_start parameter: 
            default value, bool will true else it will false
            true: t_start = tmin
            false: t_start will give by user*/
        bool t_start_default = true;
        /*Bool variable for t_end parameter:
            if it is true, t_end will calculate relative to tmin (use tend_rt_tmin option)
            if it is false, t_end will give by the user (use tend option)
            */
        bool t_end_relative_to_tmin = true;

        /*Bool varioable for build-in fout calculation*/
        bool is_bifout = false;

        

        //Init function for Parameters struct
        void init(double dtime);
};

/*Class to set fout calculation*/
class Fout_calculation
{
    public:
        Fout_calculation();
        virtual void calculate_fout(double time) = 0;
        virtual Fout_calculation* clone(Parameters& where) = 0;
        virtual ~Fout_calculation();
};

    //Independent fout parameter: set by user
    class No_fout_relation : public Fout_calculation
    {
        public:
            No_fout_relation(Parameters& _param_obj);
            void calculate_fout([[maybe_unused]] double time);
            Fout_calculation* clone(Parameters& where);
            ~No_fout_relation();
        private:
            Parameters& param_obj;
    };

    //Build in fout parameter: calculate from equation was presented at Lodato and Rossi 2011 paper, eq 28
    class Fout_Lodato_Rossi_eq28 : public Fout_calculation
    {
        public:
            Fout_Lodato_Rossi_eq28(Parameters& _param_obj);
            void calculate_fout(double time);
            Fout_calculation* clone(Parameters& where);
            ~Fout_Lodato_Rossi_eq28();
        private:
            Parameters& param_obj;
    };

/*Mdot peak calculation*/
class Mdot_peak_calculation
{
    public:
        Mdot_peak_calculation();
        virtual void calc_mdotpeak() = 0;
        virtual Mdot_peak_calculation* clone(Parameters& where) = 0;
        virtual void set_default_politrop() = 0;
        virtual void set_politrop(std::string gamma) = 0;
        virtual ~Mdot_peak_calculation();
};

    //Lodati and Rossi calculation: 2011 paper, eq 6
    class Mdotpeak_LodatoRossi : public Mdot_peak_calculation
    {
        public:
            Mdotpeak_LodatoRossi(Parameters& _param_obj);
            void calc_mdotpeak();
            Mdot_peak_calculation* clone(Parameters& where);
            void set_default_politrop();
            void set_politrop([[maybe_unused]]std::string gamma);
            ~Mdotpeak_LodatoRossi();
        private:
            Parameters& param_obj;
    };

    //Guillochon equation to calculate Mdotpeak
    class Mdotpeak_GuillochonRamirez : public Mdot_peak_calculation
    {
        public:
            Mdotpeak_GuillochonRamirez(Parameters& _param_obj);
            void calc_mdotpeak();
            double Mdotpeak_politrop5per3();
            double Mdotpeak_politrop4per3();
            Mdot_peak_calculation* clone(Parameters& where);
            void set_default_politrop();
            void set_politrop(std::string gamma);
            ~Mdotpeak_GuillochonRamirez();
        private:
            Parameters& param_obj;
            std::string politrop;
    };



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
                double calc_temperature_ph(double t);
                double calc_radius_ph(double t);
                double calc_Tph_at_rL(double t);
                double calc_luminosity_at_rldist_nu(double t, double nu);
                double calc_luminosity_at_nu(double t, double nu);

                /*Get photospheric temperature and radius*/
                double get_phot_temp(){return photospheric_temperature;}
                double get_phot_rad(){return photospheric_rad;}

                Parameters& par;

            private:
                double photospheric_temperature;
                double photospheric_rad;
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

            void set_time(double t);
            void set_nu(double new_nu);
            double calc_Mdot(double t);     //returns kg/s unit
            double calc_Td(double r);
            double calc_Td_at_rhalf();
            double calc_F_at_r(double r, double nu);

            double calc_Fx(double x);

            Parameters& par;
        
        private:
            double time;
            double mdot;
            double own_nu;
    };

    class Disk_part_of_lightcurve
    {
        public:
            Disk_part_of_lightcurve(Parameters& p);
            Disk_part_of_lightcurve();
            double calc_disk_part_at_t(double time);
            double calc_disk_part_at_t_nu(double time, double nu);
            void set_nu(double new_nu);

            double get_Tdisk_at_rhalf(double time);


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
            void print_output_before_rphlimit(double time, std::ofstream & out);
            void print_output_after_rphlimit(double time, std::ofstream & out);


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
            double calc_Lt_integral(double t);
            double calc_Fx(double x);
            double L_inp(double t);
            double L_inp_new(double t);
        private:
            
            Parameters &par;
            Disk_part_of_lightcurve disk_part;
            Wind_part_of_lightcurve wind_part;
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
            Lum_diffusion_fx diff_fx;
            double time;
            double time_offset;
            double previous_time;
            double previous_integral;
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




    