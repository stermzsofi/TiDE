#ifndef TRAPEZOIDAL
#define TRAPEZOIDAL

#include <iostream>
#include <cmath>
#include <memory>


/*Fx class: this is a mother class. Children classes must be a calc_Fx function, which calculate the value of the function at x*/
class Fx
{
    public:
        virtual double calc_Fx(double x) = 0;
};

	/*Some classes to the test*/
	/*Line class. It derive from fx. It definie a line, f(x)=a*x+b*/
	class Line : public Fx
	{
		public:
			Line(double, double);
			double calc_Fx(double x);
		protected:
			double a;
			double b;
	};


	/*Squareroot class. It derive from fx. It definie sqrt(x)+a*/
	class Squarerootmy : public Fx
	{
		public:
			Squarerootmy(double);
			double calc_Fx(double x);
		private:
			double a;
	};

	/*Power_function class. f(x)=A*x^B+C*/

	class Power_function : public Fx
	{
		public:
			Power_function(double, double, double);
			double calc_Fx(double x);
		private:
			double A;
			double B;
			double C;
	};


/*Trapezoidal rule class*/
	/* It need an Fx derived class
	 * and 2 double (the start of the integrate and the end of it
	 * An int variable, which set how many part do you want to divide the integrate (partions)
	 * We can set the start and the end of the integrate, so we don't need many classes for different ranges
	 * calc_integral will return the integral for the current setups
	 * */
class Trapezoidal_rule
{
    public:
        Trapezoidal_rule(Fx& fx, double xstart, double xend, double d_x, int n);
		Trapezoidal_rule(Fx& fx, double xstart, double xend, int n);
		Trapezoidal_rule(){}
        void set_ranges(double sxtart, double xend);
        void set_xstart(double xstart);
        void set_xend(double xend);
        void set_partions(int n);
		void set_newfx(Fx& fx);
        double calc_integral();
		double calc_integral_dx();
    private:
        Fx* fv;
        double x_start;
        double x_end;
		double dx;
        int partions;
};

#endif
