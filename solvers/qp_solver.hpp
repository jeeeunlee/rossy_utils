#pragma once

#include <iostream>
#include <Eigen/Dense>


namespace rossy_utils {
    double qpprog(Eigen::MatrixXd& _G, Eigen::VectorXd& _g0,
                      const Eigen::MatrixXd& _CE, const Eigen::VectorXd& _ce0,
                      const Eigen::MatrixXd& _CI, const Eigen::VectorXd& _ci0,
                      Eigen::VectorXd& _x);

    // add qpsolver from HiGHS later
}