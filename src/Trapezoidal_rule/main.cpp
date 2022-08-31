#include "trapezoidal.hpp"

/*Test file for test trapezoidal_rule class*/

int main()
{
    //Line l(2,5);
    Squarerootmy s(3);
    Power_function p(2, 3.2, 5);
    Trapezoidal_rule integral(s, 0,20,200);
    Quad_Trapezoidal fejlett_trapez(s,0,20,0.01);
    try
    {
        double int1_f = fejlett_trapez.qtrap();
        double int1_n = integral.calc_integral();
    
        //double int1 = integral.calc_integral();
        std::cout << "fejlett: " << int1_f << "\tnormal: " << int1_n << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    catch(char const* hiba)
    {
        std::cerr << hiba << std::endl;
    }
    
    
    
    return 0;
}
