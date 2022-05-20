#include "GradMethods.h"
#include "Matrix.h"
#include "PenaltyMethod.h"



class f1 :public Func {

public:
    double f(std::vector<double> x) {
        return x[0] * x[0] + x[1] * x[1] + 2 * sin((x[0] - x[1]) / 2);
    }
    Matrix* H(std::vector<double> x) {
        Matrix* Hm = new Matrix(2, 2);
        Hm->set(0, 0, 2 - 0.5 * sin((x[0] - x[1]) / 2));
        Hm->set(1, 0, 0.5 * sin((x[0] - x[1]) / 2));
        Hm->set(0, 1, 0.5 * sin((x[0] - x[1]) / 2));
        Hm->set(1, 1, 2 - 0.5 * sin((x[0] - x[1]) / 2));
        return Hm;
    }
    Matrix* grad(std::vector<double> x) {
        Matrix* gr_M = new Matrix(1, 2);
        gr_M->set(0, 0, 2 * x[0] + cos((x[0] - x[1]) / 2));
        gr_M->set(0, 1, 2 * x[1] - cos((x[0] - x[1]) / 2));
        return gr_M;
    }
};

class f2 : public Func_constraint {
public:
    double c_val = 1;
    double f(std::vector<double> x) {
        val++;
        return x[0] * x[0] + x[1] * x[1] + c_val * std::pow(std::max(0.0, g(x)), 3);
    }
    double f_original(std::vector<double> x) {
        return x[0] * x[0] + x[1] * x[1];
    }
    double g(std::vector<double> x) {
        return (2 * x[0] + x[1] + 4);
    }
    Matrix* H(std::vector<double> x) {
        val++;
        Matrix* Hm = new Matrix(2, 2);
        Hm->set(0, 0, 2 + 6 * c_val * std::max(0.0, g(x)) * 4);
        Hm->set(1, 0, 6 * c_val * std::max(0.0, g(x)) * 2);
        Hm->set(0, 1, 6 * c_val * std::max(0.0, g(x)) * 2);
        Hm->set(1, 1, 2 + 6 * c_val * std::max(0.0, g(x)));
        return Hm;
    }
    Matrix* grad(std::vector<double> x) {
        val++;
        Matrix* gr_M = new Matrix(1, 2);
        gr_M->set(0, 0, 2 * x[0] + 3 * std::pow(c_val * std::max(0.0, g(x)), 2) * 2);
        gr_M->set(0, 1, 2 * x[1] + 3 * std::pow(c_val * std::max(0.0, g(x)), 2) * 1);
        return gr_M;
    }
    void iter() {
        c_val *= 10;
    };

    double getValC() {
        return c_val;
    }
};







int main() {

    f2* f = new f2;
    std::vector<double> x = Penaltymethod(f);
    std::cout << f->f(x) << std::endl;


    return 0;
}

