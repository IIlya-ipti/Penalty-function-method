#include "Matrix.h"



Matrix* np::concatenateX(const Matrix* one, Matrix* two) {
    std::vector<double> total_vec(one->getVector().size() + two->getVector().size());
    int column1 = one->numCol();
    int column2 = two->numCol();
    int lines1 = one->numLines();
    int lines2 = two->numLines();
    for (int i = 0; i < lines1; i++) {
        for (int j = 0; j < column1 + column2; j++) {
            if (j < column1) {
                total_vec[(column1 + column2) * i + j] = one->getVector()[i * column1 + j];
            }
            else {
                total_vec[(column1 + column2) * i + j] = two->getVector()[i * column2 + j - column1];
            }
        }
    }
    return new Matrix(total_vec, column1 + column2, lines1);
}
Matrix* np::inv(const Matrix* matr) {
    Matrix* ones = new Matrix(std::vector<double>(matr->numLines() * matr->numLines(), 0), matr->numLines(), matr->numLines());
    for (int i = 0; i < ones->numLines(); i++) {
        ones->set(i, i, 1);
    }
    Matrix* m = np::concatenateX(matr, ones);
    double diag, val;
    for (int i = 0; i < m->numLines(); i++) {
        diag = m->get(i, i);
        for (int j = 0; j < m->numCol(); j++) {
            m->set(j, i, m->get(j, i) / diag);
        }
        for (int j = i + 1; j < m->numLines(); j++) {
            val = m->get(i, j);
            for (int k = 0; k < m->numCol(); k++) {
                m->set(k, j, m->get(k, j) - val * m->get(k, i));
            }
        }
    }
    for (int j = 1; j < m->numLines(); j++) {
        for (int i = j; i < m->numLines(); i++) {
            val = m->get(i, j - 1);
            for (int k = i; k < m->numCol(); k++) {
                m->set(k, j - 1, m->get(k, j - 1) - val * m->get(k, i));
            }
        }
    }
    Matrix* mm = new Matrix(matr->numCol(), matr->numLines());
    for (int i = 0; i < m->numLines(); i++) {
        for (int j = mm->numCol(); j < m->numCol(); j++) {
            mm->set(j - mm->numCol(), i, m->get(j, i));
        }
    }
    return mm;
}
std::ostream& operator << (std::ostream& out, const Matrix& matrix) {
    out << "Matrix : " << "\n";
    for (int i = 0; i < matrix.numLines(); i++) {
        for (int j = 0; j < matrix.numCol(); j++) {
            out << matrix.get(j, i) << ' ';
        }
        out << '\n';
    }
    return out;
}
Matrix* operator*(const Matrix one, const Matrix two) {
    std::vector<double>new_matrix(one.numLines() * two.numCol());
    for (int z = 0; z < two.numCol(); z++) {
        for (int j = 0; j < one.numLines(); j++) {
            for (int i = 0; i < one.numCol(); i++) {
                new_matrix[two.numCol() * j + z] += one.get(i, j) * two.get(z, i);
            }
        }
    }
    return new Matrix(new_matrix, two.numCol(), one.numLines());
};
Matrix* operator*(const Matrix one, const double val) {
    std::vector<double>new_matrix(one.numLines() * one.numCol());
    for (int i = 0; i < new_matrix.size(); i++) {
        new_matrix[i] = one.getVector()[i] * val;
    }
    return new Matrix(new_matrix, one.numCol(), one.numLines());
};
Matrix operator-(Matrix one, Matrix two) {
    std::vector<double> res(one.numCol() * one.numLines());
    Matrix res_matrix = Matrix(res, one.numLines(), one.numCol());
    for (int i = 0; i < one.numCol(); i++) {
        for (int j = 0; j < one.numLines(); j++) {
            res_matrix.set(one.get(i, j) - two.get(i, j), i, j);
        }
    }
    return res_matrix;
}