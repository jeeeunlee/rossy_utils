#pragma once

#include <iostream>
#include <Eigen/Dense>


// solve LP, linprog(f,A,b)
// min f'x s.t. Ax<=>b
class HighsSparseMatrix;
class Highs;

namespace rossy_utils {
    double linprog(const Eigen::VectorXd& f,
                const Eigen::MatrixXd& A,
                const Eigen::VectorXd& b,
                Eigen::VectorXd& x);

    void EigenMatrix2PackedMat(const Eigen::MatrixXd& emat,
                    HighsSparseMatrix& packedMat);
    void printOptimizerMsgs(const Highs& highs);

}// namespace rossy_utils
