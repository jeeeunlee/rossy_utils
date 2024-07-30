#pragma once

#include <iostream>
#include <Eigen/Dense>

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
    double qpprogHiGHS(const Eigen::MatrixXd & _Q, const Eigen::VectorXd & _q, 
                    const Eigen::MatrixXd & _A, const Eigen::VectorXd & _b, 
                    Eigen::VectorXd & _x);

    void EigenMatrix2Hessian(const Eigen::MatrixXd& H,
                            HighsHessian& hessian);
}