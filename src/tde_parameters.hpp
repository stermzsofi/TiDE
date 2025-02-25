#ifndef TDE_PARAMETERS
#define TDE_PARAMETERS

#include "tde_lcurve_assistant.hpp"
#include "Trapezoidal_rule/trapezoidal.hpp"
#include "Linear_interpol/linear_interpol.hpp"
#include <memory>
#include <cstdlib>
#include <cmath>
#include <limits>

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
        virtual void print_settings() = 0;
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
            void print_settings();
            double get_RISCO();
        private:
            double calc_schwarzschild_radius();
            double RISCO;
            void calc_RISCO();
    };

    class Mass_to_Radii;
    
    //Star class
    class Star : public Object
    {
        public:
            Star(double m, double r, startype relation);
            Star(double m, double r, startype relation, std::string _politrop);
            Star();
            void set_startype(std::string rad);
            void set_politrop(std::string _politrop);
            startype get_startype(){return st;}
            void print_startype();
            double get_K() {return K;}
            double get_rho_c() {return rho_c;}
            double get_n() {return n;}
            std::string get_politrop();
            void init();
            void print_settings();
        private:
            Mass_to_Radii* m_t_r;
            startype st;
            std::string politrop;
            double n;
            double K;
            double rho_c;
            bool setted_politrop = false;

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
class tmin_calculation;
class Mdotfb_calculation;
class TL_calculation;
class rout_calculation;
struct Parameters;
class calc_Disk_at_r_for_bolometric : public Fx
    {
        public:
            calc_Disk_at_r_for_bolometric(Parameters& p);
            calc_Disk_at_r_for_bolometric();

            //void set_nu(double new_nu);
            double calc_Mdot();     //returns kg/s unit
            double calc_sigma_Td4(double r);
            //double calc_Td_at_rhalf();
            double calc_Lbolometric_at_r(double r);

            double calc_Fx(double x);

            Parameters& par;
        
        private:
            double mdot;
            //double own_nu;
    };

/*Fx_for_sigma class: to calculate the integral in the sigma*/
class Fx_for_sigma : public Fx
{
    public:
        Fx_for_sigma(Parameters& par);
        Fx_for_sigma();
        void set_w(double w_) {w = w_;}
        void set_w_Rc(double w_Rc_){w_Rc = w_Rc_;}
        void set_T(double T_){T = T_;}
        //void set_t(double t_){t = t_;}
        double calc_Fx(double x);
    private:
        Parameters& p;
        double w;
        double w_Rc;
        double T;
        //double t;
        double G(double z);     //Green function
        double St(double t);    //Mass source function
};

struct sigma
{
    sigma();
    sigma(Parameters& par);

    Parameters& p;
    Fx_for_sigma fx;
    Quad_Trapezoidal integ;
    double T;
    double w_R;
    double w_Rc;
    //double t;
    double R;

    double calc_T(double t);
    double calc_wR(double _R);
    double calc_sigma_R_t(double t, double _R);
};

class Fx_for_Redd : public Fx
{
    public:
        Fx_for_Redd();
        Fx_for_Redd(Parameters& par);
        void set_t(double _t) {t = _t;}
        double calc_Fx(double x);
        double last_fx;
    private:
        Parameters& p;
        sigma s;
        double t;
};

class dFx_for_Redd : public dFx
{
    public:
        dFx_for_Redd();
        dFx_for_Redd(Parameters & par, Fx_for_Redd & fx_);
        double calc_dFx(double x);
        void set_eps_f(double eps){eps_f = eps;}
    private:
        Parameters& p;
        Fx_for_Redd& fx;
        double eps_f;

};


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

/*Parameters struct sub structs*/
struct parameters_for_timedependent_rout
{
    double nu_c;    //valami viscosity: talán central?
    double alpha = 0.1;   //effective kinematic viscosity
    double HperR = 0.5;   //disk half thickness per distance from BH
    double Rc;      //circularization radius: 2rt/beta
    double n = 0.5;       //parameter in Green fuction
    double l = 3.0/4.0;       //parameter in Green function
};


/*Other needed parameters, some simple calculations*/
struct Parameters
{
        Parameters(int N_o, double e, double b, double fo, double fv_, double dist, double ink, double nu_, Star star, Blackhole blackhole);
        Parameters();
        ~Parameters();
        //Parameters(Parameters& copied);
        
        //Functions used in constructor and init
        void fv_check();  //check using fv, v_wind may is higher, than light velocity and also calculate the vwind velocity
        void calculate_rt();    //calculate rt tidal radius and also rp pericentrum distance
        void calculate_rphtmin();
        void calculate_tphtmin();
        void calculate_rin();
        void calculate_rout();
		void calculate_rhalf();
        void calculate_rL();
        void calculate_mdot_edd();  //in kg/s
        void calculate_time_wind_limit();
        void calculate_time_rph_limit();
        //void calculate_tpeak();
        void calculate_N_col();
        void calculate_n();
        void calculate_Rmax();
        void calculate_deltaomega();
    
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
        void change_tmincalculation(tmin_calculation* new_tmincalc);
        void change_Mdotfbcalculation(Mdotfb_calculation* new_calcMdotfb);
        void change_TLcalculation(TL_calculation* new_TLcalculation);
        void change_routcalculation(rout_calculation* new_routcalculation);

        /*Parameters*/
        double eta;
        double eta_reprocessing;
        double beta;
        double fout;
        double fv;
        double d;
        double i;
        double diffusion_timescale;
        int N;
        double beta_fulldisrupt;    //the real beta value of the full disruption
        /*needed classes*/
        Star s;
        Blackhole bh;
        std::unique_ptr<Fout_calculation> calcfout;
        std::unique_ptr<Mdot_peak_calculation> calcmdotpeak;
        std::unique_ptr<tmin_calculation> calctmin;
        std::unique_ptr<Mdotfb_calculation> calcMdotfb;
        std::unique_ptr<TL_calculation> calcTL;
        std::unique_ptr<rout_calculation> calcrout;

        /*Lightcurve intervalls and step*/
        double t_start;
        double t_end;
        double dt;
        double nu;

        /*Spectra intervalls and step*/
        double nu_start;
        double nu_end;
        double dnu;
        double time;    //also used as time dependent parameter when lightcurve is calculated
        
        /*Derived quantities*/
        double rt;      //tidal radius [m]
        double rp;      //pericentrum distance [m]
        double tmin;
        double rphtmin;
        double mdot_edd;
        double tphtmin;
        double rin;
        double rout;
        double rhalf;
        double rL;
        double vw;
        //double tpeak_per_tmin;
        //time limit until wind will exist
        //elvileg ezt a kettőt is ki kéne irtani majd
        double time_wind_limit;
        //time limit while rph = rl
        double time_rph_limit;
        double oneper4pid2;
        double cosi;
        double Mdot_peak;
        double delta_omega; //Deltaomega in Strubbe: in sr

        /*Needed classes and functions for calculate Ldisk bolometric luminosity at time*/
        calc_Disk_at_r_for_bolometric disk_bol_calc;
        Quad_Trapezoidal Ldiskbol;
        void Ldiskbol_init();
        void calc_Ldisk();
        

        /*Some time dependent values, but it useful to speed up the calculations*/
        double Mdotfb_t;
        double Ldisk_bolometric_at_t;
        double TL;
        double n; //number density : in cm^-3
        double N_col; //column density: in cm^-2
        double Rmax; //Rmax in Strubbe: in cm
        //double t;

        //Mdotfb_t, fout, time, TL
        void calculate_timedependent_parameters(double t);

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

        /**nagyon nagyon randa: ez ahhoz, hogy a cloudy-s reprocesst használja diffúziónál*/
        bool my_rep = false;

        /*Bool varioable for build-in fout calculation*/
        bool is_bifout = false;

        //Init function for Parameters struct
        void init(/*double dtime*/);
        void refresh(double dtime);
        void init_tmin();

        /*wind radius*/
        std::string wind_radius = "rph";

        parameters_for_timedependent_rout time_dep_rout;
};

/*Class to set fout calculation*/
class Fout_calculation
{
    public:
        Fout_calculation();
        virtual void calculate_fout() = 0;
        virtual Fout_calculation* clone(Parameters& where) = 0;
        virtual ~Fout_calculation();
        virtual void print_settings() = 0;
};

    //Independent fout parameter: set by user
    class No_fout_relation : public Fout_calculation
    {
        public:
            No_fout_relation(Parameters& _param_obj);
            void calculate_fout();
            Fout_calculation* clone(Parameters& where);
            ~No_fout_relation();
            void print_settings();
        private:
            Parameters& param_obj;
    };

    //Build in fout parameter: calculate from equation was presented at Lodato and Rossi 2011 paper, eq 28
    class Fout_Lodato_Rossi_eq28 : public Fout_calculation
    {
        public:
            Fout_Lodato_Rossi_eq28(Parameters& _param_obj);
            void calculate_fout();
            Fout_calculation* clone(Parameters& where);
            ~Fout_Lodato_Rossi_eq28();
            void print_settings();
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
        virtual ~Mdot_peak_calculation();
};

    //Lodati and Rossi calculation: 2011 paper, eq 6
    class Mdotpeak_LodatoRossi : public Mdot_peak_calculation
    {
        public:
            Mdotpeak_LodatoRossi(Parameters& _param_obj);
            void calc_mdotpeak();
            Mdot_peak_calculation* clone(Parameters& where);
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
            ~Mdotpeak_GuillochonRamirez();
        private:
            Parameters& param_obj;
    };

    class Mdotpeak_L09_constrho : public Mdot_peak_calculation
    {
        public:
            Mdotpeak_L09_constrho(Parameters& _param_obj);
            void calc_mdotpeak();
            Mdot_peak_calculation* clone(Parameters& where);
            ~Mdotpeak_L09_constrho();
        private:
            Parameters& param_obj;
    };

    class Mdotpeak_L09_all : public Mdot_peak_calculation
    {
        public:
            Mdotpeak_L09_all(Parameters& _param_obj);
            void calc_mdotpeak();
            Mdot_peak_calculation* clone(Parameters& where);
            ~Mdotpeak_L09_all();
        private:
            Parameters& param_obj;
    };

/*tmin calculation class*/
class tmin_calculation
{
    public:
        tmin_calculation(Parameters& _param_obj);
        virtual void calc_tmin() = 0;
        virtual tmin_calculation* clone(Parameters& where) = 0;
        virtual ~tmin_calculation();
        virtual void print_setted_type() = 0;
    protected:
        Parameters& param_obj;
};

    class tmin_Lodato : public tmin_calculation
    {
        public:
            tmin_Lodato(Parameters& _param_obj);
            tmin_calculation* clone(Parameters& where);
            ~tmin_Lodato();
            void calc_tmin();
            void print_setted_type();
    };

    class tmin_with_rt : public tmin_calculation
    {
        public:
            tmin_with_rt(Parameters& _param_obj);
            tmin_calculation* clone(Parameters& where);
            ~tmin_with_rt();
            void calc_tmin();
            void print_setted_type();
    };

    class tmin_from_GR13 : public tmin_calculation
    {
        public:
            tmin_from_GR13(Parameters& _param_obj);
            tmin_calculation* clone(Parameters& where);
            ~tmin_from_GR13();
            double B4per3();
            double B5per3();
            void calc_tmin();
            void print_setted_type();
    };

    class n_inf_calculation
    {
        public:
            n_inf_calculation(Parameters& _param_obj);
            virtual ~n_inf_calculation();
            virtual n_inf_calculation* clone(Parameters& where) = 0;
            virtual double get_ninf() = 0; 
        protected:
            Parameters& param_obj;
    };

        class n_inf_const : public n_inf_calculation
        {
            public:
                n_inf_const(Parameters& _param_obj);
                n_inf_const(Parameters& _param_obj, double _ninf);
                ~n_inf_const();
                void set_ninf(double _ninf);
                n_inf_calculation* clone(Parameters& where);
                double get_ninf();
            private:
                double ninf;
        };

        class n_inf_GR13 : public n_inf_calculation
        {
            public:
                n_inf_GR13(Parameters& _param_obj);
                ~n_inf_GR13();
                n_inf_calculation* clone(Parameters& where);
                double D_5per3();
                double D_4per3();
                double get_ninf();
        };

    class Mdotfb_calculation
    {
        public:
            Mdotfb_calculation(Parameters& _param_obj);
            ~Mdotfb_calculation();
            virtual Mdotfb_calculation* clone(Parameters& where) = 0;
            virtual void calc_Mdotfb_at_t(double t) = 0;
            virtual void init() = 0;
            virtual void refresh() = 0;
            virtual void print_setted_type() = 0;
            void set_calc_ninf(std::string ninf);
            double get_ninf(){return n_inf;}
        protected:
            Parameters& param_obj;
            std::unique_ptr<n_inf_calculation> calc_ninf;
            double n_inf;   //az idő hatványkitevője
    };

        //in this case Mdotfb = Mdotp * (t/tmin)^(n_inf)
        class Mdotfb_classical_powerfunction : public Mdotfb_calculation
        {
            public:
                Mdotfb_classical_powerfunction(Parameters& _param_obj);
                ~Mdotfb_classical_powerfunction();
                Mdotfb_calculation* clone(Parameters& where);
                void calc_Mdotfb_at_t(double t);
                void print_setted_type();
                void init() {}
                void refresh() {}
        };

        //Mdotfb = Mdotp * (1 - (t/tmin)^(-4/3)) * (t/tmin)^(n_inf)
        class Mdotfb_L09_constrho_x : public Mdotfb_calculation
        {
            public:
                Mdotfb_L09_constrho_x(Parameters& _param_obj);
                ~Mdotfb_L09_constrho_x();
                Mdotfb_calculation* clone(Parameters& where);
                void calc_Mdotfb_at_t(double t);
                void print_setted_type();
                void init(){}
                void refresh() {}
        };

        class Mdotfb_L09_all : public Mdotfb_calculation
        {
            public:
                Mdotfb_L09_all(Parameters& _param_obj);
                ~Mdotfb_L09_all();
                Mdotfb_calculation* clone(Parameters& where);
                void calc_Mdotfb_at_t(double t);
                void calc_constpart();
                void print_setted_type();
                void init();
                void refresh();
            private:
                double const_part;
                Interpolator interpolate;
        };

    /*TL calculation
        with reprocessing or not*/
    
    class TL_calculation
    {
        public:
            TL_calculation(Parameters& _param_obj);
            virtual double calc_TL([[maybe_unused]] double Ldisk) = 0;
            virtual TL_calculation* clone(Parameters& where) = 0;
            virtual ~TL_calculation();
        protected:
            Parameters& param_obj;
    };

        class TL_without_reprocessing : public TL_calculation
        {
            public:
                TL_without_reprocessing(Parameters& _param_obj);
                TL_calculation* clone(Parameters& where);
                ~TL_without_reprocessing();
                double calc_TL([[maybe_unused]] double Ldisk);
        };

        class TL_with_reprocessing : public TL_calculation
        {
            public:
                TL_with_reprocessing(Parameters& _param_obj);
                TL_calculation* clone(Parameters& where);
                ~TL_with_reprocessing();
                double calc_TL([[maybe_unused]] double Ldisk);
        };

    class rout_calculation
    {
        public:
            //rout_calculation();
            rout_calculation(Parameters& par);
            virtual double calc_rout([[maybe_unused]] double t) = 0;
        protected:
            Parameters& p;
    };

        class constant_rout : public rout_calculation
        {
            public:
                constant_rout(Parameters& par);
                double calc_rout( [[maybe_unused]] double t);
        };

        class time_dependent_rout : public rout_calculation
        {
            public:
                time_dependent_rout(Parameters& par);
                double calc_rout(double t);
            private:
                Fx_for_Redd fx;
                dFx_for_Redd dfx;
                Newthon_Raphson iter;
        };

    

#endif
