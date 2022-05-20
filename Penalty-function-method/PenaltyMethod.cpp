#include "PenaltyMethod.h"

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
std::vector<double> Penaltymethod(Func_constraint* f) {

    std::cout << std::fixed;
    std::cout.precision(15);

    // set start random vector
    std::vector<double> random_vector = { fRand(-10, 10), fRand(-10, 10) };
    f->setC(random_vector);


    const double EPS = 0.000001;

    int index = 0;
    std::vector<double> x;
    do {
        index++;
        x = minGrad2P(f);
        std::cout << "g(x) = " << f->g(x) << std::endl;
        std::cout << "iter = " << index << std::endl;
        std::cout << "function calls = " << Func_constraint::val << std::endl;
        std::cout << "----" << std::endl;
        f->iter();
    } while (f->g(x) > EPS);

    return x;
}