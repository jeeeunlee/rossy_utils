#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

class HighsHessian;

namespace rossy_utils {
    double qpprog(Eigen::MatrixXd& _G, Eigen::VectorXd& _g0,
                      const Eigen::MatrixXd& _CE, const Eigen::VectorXd& _ce0,
                      const Eigen::MatrixXd& _CI, const Eigen::VectorXd& _ci0,
                      Eigen::VectorXd& _x);
    double qpprog(Eigen::MatrixXd& G, Eigen::VectorXd& g0,
                const Eigen::MatrixXd& CI, const Eigen::VectorXd& ci0,
                Eigen::VectorXd& x);

    // add qpsolver from HiGHS later
    // f =   min 0.5 * x Q x + q x, where Q=Q'(symmetric, psd)
    // s.t. _A x + <= _b 
    // return f
    double qpprogHiGHS(const Eigen::MatrixXd & _Q, const Eigen::VectorXd & _q, 
                    const Eigen::MatrixXd & _A, const Eigen::VectorXd & _b, 
                    Eigen::VectorXd & _x);

    void EigenMatrix2Hessian(const Eigen::MatrixXd& H,
                            HighsHessian& hessian);

    // add qpsolver from OSQP
    float qpprogOSQPSparse(
        Eigen::SparseMatrix<float> &Q_float, 
        Eigen::VectorXf &q,
        Eigen::SparseMatrix<float> &A_float, 
        Eigen::VectorXf &b,
        Eigen::VectorXf &x);

    float qpprogOSQP(
        Eigen::MatrixXf &Q, 
        Eigen::VectorXf &q,
        Eigen::MatrixXf &A, 
        Eigen::VectorXf &b,
        Eigen::VectorXf &x);
    

}