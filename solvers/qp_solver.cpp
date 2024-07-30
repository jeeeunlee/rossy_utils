#include "qp_solver.hpp"
#include "rossy_utils/thirdparty/goldfarb/QuadProg++.hh"
#include "lp_solver.hpp"
#include <Highs.h>
#include <cassert>

#define ZCE 1e-8
const double INF = std::numeric_limits<double>::infinity();

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

double rossy_utils::qpprog(Eigen::MatrixXd& G, Eigen::VectorXd& g0,
                      const Eigen::MatrixXd& CI, const Eigen::VectorXd& ci0,
                      Eigen::VectorXd& x)
{
//     min 0.5 * x G x + g0 x
// s.t.
//     CE^T x + ce0 = 0
//     CI^T x + ci0 >= 0
Eigen::MatrixXd CE = Eigen::MatrixXd::Zero(0,0);
Eigen::VectorXd ce0 = Eigen::VectorXd::Zero(0);
 return solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
}

double rossy_utils::qpprogHiGHS(const Eigen::MatrixXd & Q, 
                    const Eigen::VectorXd & q, 
                    const Eigen::MatrixXd & A, 
                    const Eigen::VectorXd & b, 
                    Eigen::VectorXd & x)
{
//   f =   min 0.5 * x Q x + q x, where Q=Q'(symmetric, psd)
// s.t.
//     _A x + <= _b 
// return f
    int Nvar = q.size();
    int Nconst = A.rows();

    HighsModel model;

    model.lp_.sense_ = ObjSense::kMinimize;
    model.lp_.offset_ = 0.;
    // # of variables
    model.lp_.num_col_ = Nvar;
    // # of constraints
    model.lp_.num_row_ = Nconst;  
    // set Q
    EigenMatrix2Hessian(Q, model.hessian_);
    // set q
    model.lp_.col_cost_.resize(Nvar); 
    Eigen::VectorXd::Map(&model.lp_.col_cost_[0], Nvar) = q;
    model.lp_.col_lower_.resize(Nvar);
    model.lp_.col_upper_.resize(Nvar);
    for(int i(0); i<Nvar; ++i){
        model.lp_.col_lower_[i] = -INF;
        model.lp_.col_upper_[i] = INF;
    }  

    // set constraints A, b
    EigenMatrix2PackedMat(A, model.lp_.a_matrix_);
    model.lp_.row_lower_.resize(Nconst);
    model.lp_.row_upper_.resize(Nconst);
    for(int i(0); i<Nconst; ++i){
        model.lp_.row_lower_[i] = -INF;
        model.lp_.row_upper_[i] = b(i);
    }

    // Create a Highs instance and solve the problem
    Highs highs;
    HighsStatus return_status;
    // return_status = highs.passHessian(model.hessian_);
    // assert(return_status==HighsStatus::kOk);
    assert(model.isQp());
    return_status = highs.passModel(model);    
    assert(return_status==HighsStatus::kOk);    
    return_status = highs.run();
    assert(return_status==HighsStatus::kOk);

    const HighsSolution& solution = highs.getSolution();
    x = Eigen::Map<const Eigen::VectorXd>(
        solution.col_value.data(), solution.col_value.size());
    const HighsInfo& info = highs.getInfo();
    return info.objective_function_value;
}

void rossy_utils::EigenMatrix2Hessian(const Eigen::MatrixXd& H,
    HighsHessian& hessian){
    // Assume Symmetric H = H'
    // hessian.format_ = HessianFormat::kTriangular; // 
    hessian.dim_ = H.cols();
    hessian.format_ = HessianFormat::kSquare; // 
    hessian.start_.clear();
    hessian.index_.clear();
    hessian.value_.clear();

    int idx(0);
    for(int c(0); c<H.cols(); ++c){
        hessian.start_.push_back(idx);    
        for(int r(0); r<H.rows(); ++r){                    
            if( std::abs(H(r,c)) > ZCE)
            {                
                hessian.index_.push_back(r);
                hessian.value_.push_back(H(r,c));
                idx ++;
            }
        }
    }
    hessian.start_.push_back(idx);
    // check
    // std::cout<<" hessian.start_ = " << std::endl;
    // for (auto &s : hessian.start_)    
    //     std::cout<< s << ", ";
    // std::cout<<std::endl;

    // std::cout<<" hessian.index_ = " << std::endl;
    // for (auto &s : hessian.index_)    
    //     std::cout<< s << ", ";
    // std::cout<<std::endl;

    // std::cout<<" hessian.value_ = " << std::endl;
    // for (auto &s : hessian.value_)    
    //     std::cout<< s << ", ";
    // std::cout<<std::endl;
}