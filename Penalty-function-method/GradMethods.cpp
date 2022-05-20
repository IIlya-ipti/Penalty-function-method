#include "GradMethods.h"

int argmin(std::vector<double> vals) {
    double mn = vals[0];
    int index = 0;
    for (int i = 1; i < vals.size(); i++) {
        if (vals[i] < mn) {
            mn = vals[i];
            index = i;
        }
    }
    return index;
}

double minValue1D(double (Func::* ptrfunc)(double), Func* f, double minX, double maxX) {
    const double eps = 0.00001;

    // number of x
    const int n = 100;

    // array of func value
    std::vector<double> fs(n);

    // step of method
    double step, mn_val = (f->*ptrfunc)(minX);
    int index;

    do {
        step = (maxX - minX) / n;

        for (int i = 1; i <= n; i++) {
            fs[i - 1] = (f->*ptrfunc)(minX + step * i);
        }
        index = argmin(fs);
        minX = minX + step * (index - 1.0);
        maxX = minX + 2 * step;
    } while (abs(minX - maxX) > eps);

    return (minX + maxX) / 2;
}
std::vector<double> minGrad2P(Func* f) {
    double eps = 0.0001;
    double l;
    int steps = 0;
    std::vector<double> x = f->getC();
    Matrix* dx;
    Matrix* dv;
    do {
        f->setC(x);
        l = minValue1D(&(Func::f2P), f, 0.00000000000001, 0.99999999999999);
        dx = *(*np::inv(f->H(x)) * *(f->grad(x))) * l;
        dv = f->grad(x);
        x[0] = x[0] - dx->get(0, 0);
        x[1] = x[1] - dx->get(0, 1);
        steps += 1;
        if (steps == 10)break;
    } while (abs(std::max<double>(dv->get(0, 0), dv->get(0, 1))) > eps);
    f->setC(x);
    return x;
}
std::vector<double> minGrad2P(Func_constraint* f) {
    double eps = 0.0001;
    double l;
    int steps = 0;
    std::vector<double> x = f->getC();
    Matrix* dx;
    Matrix* dv;
    do {
        f->setC(x);
        l = 1 / f->getValC();
        dx = *(*np::inv(f->H(x)) * *(f->grad(x))) * l;
        dv = f->grad(x);
        x[0] = x[0] - dx->get(0, 0);
        x[1] = x[1] - dx->get(0, 1);
        steps += 1;
        if (steps == 10)break;
    } while (abs(std::max<double>(dv->get(0, 0), dv->get(0, 1))) > eps);
    f->setC(x);
    return x;
}
std::vector<double> minGrad1P(Func* f) {

    // constants
    double eps = 0.00001;

    // value of alpha k (learning rate)
    double l;

    int steps = 0;

    // start position
    std::vector<double> x = f->getC();


    Matrix* dx;
    do {
        steps += 1;
        f->setC(x);
        l = minValue1D(&(Func::f1P), f, 0.00000000000001, 0.99999999999999);
        dx = *f->grad(x) * l;
        x[0] = x[0] - dx->get(0, 0);
        x[1] = x[1] - dx->get(0, 1);
    } while (abs(std::max<double>(abs(dx->get(0, 0)), abs(dx->get(0, 1)))) > eps);

    return x;
}


