#pragma once
#include "Matrix.h"


/*
    function interface
*/
class Func {
private:
    std::vector<double> c{ 1,1 };
public:
    static inline int val = 0;
    double f2P(double l) {
        std::vector<double> c_dc = c;
        Matrix* dc = *(*np::inv(H(c)) * *(grad(c))) * l;
        c_dc[0] = c[0] - dc->get(0, 0);
        c_dc[1] = c[1] - dc->get(0, 1);
        return f(c_dc) - f(c);
    }
    double f1P(double l) {
        std::vector<double> c_dc = c;
        Matrix* dc = *grad(c) * l;
        c_dc[0] = c[0] - dc->get(0, 0);
        c_dc[1] = c[1] - dc->get(0, 1);
        return f(c_dc);
    }
    void setC(std::vector<double> c) {
        this->c = c;
    }
    std::vector <double> getC() {
        return c;
    }

    /*
     function for opt
    */
    virtual double f(std::vector<double> x) = 0;
    /*
        gessian of function
    */
    virtual Matrix* H(std::vector<double> x) = 0;
    /*
        gradient of function
    */
    virtual Matrix* grad(std::vector<double> x) = 0;
};


/*
    function with constraint interface
*/
class Func_constraint : public Func {
public:
    /*
        constraint function
    */
    virtual double g(std::vector<double> val) = 0;

    virtual void iter() = 0;

    /*
        weight of constraint
    */
    virtual double getValC() = 0;
};

/*
    index of min value
*/
int argmin(std::vector<double> vals);

/*
    1D optimization
*/
double minValue1D(double (Func::* ptrfunc)(double), Func* f, double minX, double maxX);

/*
    Gradient method of the second order
*/
std::vector<double> minGrad2P(Func* f);
std::vector<double> minGrad2P(Func_constraint* f);

/*
    Gradient method of the first order
*/
std::vector<double> minGrad1P(Func* f);
