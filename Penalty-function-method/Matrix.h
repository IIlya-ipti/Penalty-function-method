#pragma once

#include "includes.h"

class Matrix {
private:
    int columns;
    int lines;
    std::vector<double> value;
public:
    Matrix(std::vector<double> val, int columns, int lines) {
        this->value = val;
        this->lines = lines;
        this->columns = columns;
    }
    Matrix() {
        this->lines = 0;
        this->columns = 0;
    };
    Matrix(int columns, int lines) {
        value = std::vector<double>(lines * columns);
        this->columns = columns;
        this->lines = lines;
    }
    std::vector<double> getVector() const {
        return value;
    }
    int numCol() const {
        return columns;
    }
    int numLines() const {
        return lines;
    }
    friend std::ostream& operator<< (std::ostream& out, const Matrix& matrix);
    friend Matrix* operator*(Matrix one, Matrix two);
    friend Matrix* operator*(const Matrix one, const double val);
    friend Matrix operator-(Matrix one, Matrix two);
    double get(int col, int lin) const {
        return value[lin * columns + col];
    }
    void set(int col, int lin, double val) {
        value[lin * columns + col] = val;
    }
};


namespace np {
    Matrix* concatenateX(const Matrix* one, Matrix* two);
    Matrix* inv(const Matrix* matr);
}