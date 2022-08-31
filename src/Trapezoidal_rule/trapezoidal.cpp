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

    Quad_Trapezoidal::Quad_Trapezoidal(Fx& _fv, double _x_start, double _x_end, double _eps) : fv(&_fv), x_start(_x_start), x_end(_x_end), eps(_eps)
    {
        n = 0;
    }

    double Quad_Trapezoidal::next()
    {
        double x, tnm, del;
        long double sum;
        int it, j;
        n++;
        if ( 1 == n )
        {
            s = 0.5 * (x_end - x_start) * (fv->calc_Fx(x_start) + fv->calc_Fx(x_end));
            return s;
        }
        else
        {
            for (it = 1, j = 1; j < n - 1; j++)
            {
                it <<=1;    //bit eltolása balra: ha jól gondolom, ezzel adja vissza, hogy hány részre van felolsztva a tartomány
            }
            tnm = it;
            del = (x_end - x_start)/tnm;
            x = x_start + 0.5 * del;
            for(sum = 0.0, j = 0; j < it; j++, x += del)
            {
                sum += fv->calc_Fx(x);
            }
            s = 0.5 * (s + (x_end - x_start) * sum/tnm);
            return s;
        }
    }

    double Quad_Trapezoidal::qtrap()
    {
        const int JMAX = 10;        //max number of steps: 2^20
        double olds = 0.0;
        for (int j = 0; j < JMAX; j++)
        {
            next();
            if (!std::isfinite(s)) throw("WARNING: S is nan!"); 
            if(j > 3)
            {
                if (std::abs(s - olds) < eps * std::abs(olds))
                {
                    n = 0;
                    //std::cout << "The needed number of j is: " << j << " corresponds to " << std::pow(2,j) << " step. The error is: " << std::abs(s)/std::abs(olds) << std::endl;
                    return s;
                }
                //std::cout << "The error is: " << std::abs(s)/std::abs(olds) << std::endl;
            }
            //if(j < JMAX - 1) 
            olds = s;
        }
        n = 0;
        throw("Too many steps in routine qtrap.");
        //std::cout << "Too many steps in routine qtrap. The error is: " << std::abs(s)/std::abs(olds) << std::endl;
        //return s;
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
        double x0, x1, int_akt,x,h;
        for(x = x_start; x <= x_end; x += dx)
        {
            x0 = x;
            x1 = x + dx;
            h=dx;
            if(x1 > x_end)
            {
                x1 = x_end;
                h=x_end-x0;
            }
            int_akt = h*(fv ->calc_Fx(x1)+fv ->calc_Fx(x0))/2.0;
            integral += int_akt;
        }
        /*if(x < x_end + dx)
        {
            x0 = x - dx;
            x1 = x_end;
            double h = x1-x0;
            int_akt = h*(fv ->calc_Fx(x1)+fv ->calc_Fx(x0))/2.0;
            integral += int_akt;
        }*/
        return integral;
    }
