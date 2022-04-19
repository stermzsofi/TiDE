#include "trapezoidal.hpp"

/*
    **********

        Test classes: for test the operation of Trapezoidal_rule

    **********
*/


/**Functions for Line class*/
    /*Constructor*/
    Line::Line(double slope, double yintercept)
    {
        a = slope;
        b = yintercept;
    }

    /*calc_Fx function*/
    double Line::calc_Fx(double x)
    {
        double res;
        res = a*x + b;
        return res;
    }

/**Functions for Squarerootmy class**/
    /*Constructor*/
    Squarerootmy::Squarerootmy(double eltol)
    {
        a = eltol;
    }

    /*calc_Fx function*/
    double Squarerootmy::calc_Fx(double x)
    {
        double res;
        res = sqrt(x)+a;
        return res;
    }



    /*Power_function class functions*/
    Power_function::Power_function(double a, double b, double c)
    {
        A = a;
        B = b;
        C = c;
    }

    double Power_function::calc_Fx(double x)
    {
        double res;
        res = A * pow(x,B) + C;
        return res;
    }

    /*
        **********

            Trapezoidal_rule class functions: 
            It will calculate the integral with trapezoidal rule

        **********
    */

    /*Constructors*/
    Trapezoidal_rule::Trapezoidal_rule(Fx& fx, double xstart, double xend, double d_x, int n)
    {
        fv = &fx;
        x_start = xstart;
        x_end = xend;
        dx = d_x;
        partions = n;
    }

    Trapezoidal_rule::Trapezoidal_rule(Fx& fx, double xstart, double xend, int n)
    {
        fv = &fx;
        x_start = xstart;
        x_end = xend;
        partions = n;
    }

    void Trapezoidal_rule::set_ranges(double xstart, double xend)
    {
        x_start = xstart;
        x_end = xend;
    }

    void Trapezoidal_rule::set_xstart(double xstart)
    {
        x_start = xstart;
    }

    void Trapezoidal_rule::set_xend(double xend)
    {
        x_end = xend;
    }

    void Trapezoidal_rule::set_partions(int n)
    {
        partions = n;
    }

    void Trapezoidal_rule::set_newfx(Fx& fx)
    {
        fv = &fx;
    }

    /**Integrate part*/
    
    //This function will divide the intervall "partions" number of equivalent parts and return with the integral value
    double Trapezoidal_rule::calc_integral()
    {
        double integral = 0;
        double h = (x_end - x_start) / partions;
        double x0, x1, int_akt;
        for(int i = 1; i<= partions; i++)
        {
            x0 = x_start + (i-1)*h;
            x1 = x_start + i*h;
            int_akt = h*(fv ->calc_Fx(x1)+fv ->calc_Fx(x0))/2.0;
            integral += int_akt;
        }
        return integral;
    }

    //This function will goes through the intervall with dx steps and return with the vaule of the integral
    double Trapezoidal_rule::calc_integral_dx()
    {
        double integral = 0;
        double x0, x1, int_akt,x;
        for(x = x_start; x <= x_end; x += dx)
        {
            x0 = x;
            x1 = x + dx;
            int_akt = dx*(fv ->calc_Fx(x1)+fv ->calc_Fx(x0))/2.0;
            integral += int_akt;
        }
        return integral;
    }
