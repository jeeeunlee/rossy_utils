#include "qp_solver.hpp"
#include "rossy_utils/thirdparty/goldfarb/QuadProg++.hh"

double rossy_utils::qpprog(Eigen::MatrixXd& G, Eigen::VectorXd& g0,
                      const Eigen::MatrixXd& CE, const Eigen::VectorXd& ce0,
                      const Eigen::MatrixXd& CI, const Eigen::VectorXd& ci0,
                      Eigen::VectorXd& x)
{
//     min 0.5 * x G x + g0 x
// s.t.
//     CE^T x + ce0 = 0
//     CI^T x + ci0 >= 0
 return solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
}