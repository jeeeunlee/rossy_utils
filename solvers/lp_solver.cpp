#include "rossy_utils/solvers/lp_solver.hpp"
#include "Highs.h"
#include <cassert>

#define ZCE 1e-8
#define INF 1e6


namespace rossy_utils {
// solve LP, linprog(f,A,b)
// min f'x s.t. Ax<=b
double linprog(const Eigen::VectorXd& f,
                const Eigen::MatrixXd& A,
                const Eigen::VectorXd& b,
                Eigen::VectorXd& x){

    HighsModel model;
    model.lp_.sense_ = ObjSense::kMinimize;
    // model.lp_.sense_ = ObjSense::kMaximize;
    model.lp_.offset_ = 0.;

    int Nvar = f.size();
    int Nconst = A.rows();

    // cost
    model.lp_.col_cost_.resize(Nvar); 
    Eigen::VectorXd::Map(&model.lp_.col_cost_[0], Nvar) = f;
    std::cout<< "model.lp_.col_cost_ = " << std::endl;
    for(auto & d : model.lp_.col_cost_) std::cout << d << ", ";
    std::cout<<std::endl;
    // # of variables
    model.lp_.num_col_ = Nvar;
    // # of constraints
    model.lp_.num_row_ = Nconst;    

    model.lp_.col_lower_.resize(Nvar);
    model.lp_.col_upper_.resize(Nvar);
    for(int i(0); i<Nvar; ++i){
        model.lp_.col_lower_[i] = -INF;
        model.lp_.col_upper_[i] = INF;
    }

    model.lp_.row_lower_.resize(Nconst);
    model.lp_.row_upper_.resize(Nconst);
    for(int i(0); i<Nconst; ++i){
        model.lp_.row_lower_[i] = -INF;
        model.lp_.row_upper_[i] = b(i);
    }

    EigenMatrix2PackedMat(A, model.lp_.a_matrix_);

    // Create a Highs instance and solve the problem
    Highs highs;
    HighsStatus return_status = highs.passModel(model);
    assert(return_status==HighsStatus::kOk);    
    return_status = highs.run();
    assert(return_status==HighsStatus::kOk);
    // printOptimizerMsgs(highs);

    const HighsSolution& solution = highs.getSolution();
    x = Eigen::Map<const Eigen::VectorXd>(solution.col_value.data(), solution.col_value.size());

    const HighsInfo& info = highs.getInfo();

    return info.objective_function_value;
}

void EigenMatrix2PackedMat(const Eigen::MatrixXd& emat,
    HighsSparseMatrix& packedMat){

    packedMat.format_ = MatrixFormat::kColwise;
    packedMat.start_.clear();
    packedMat.index_.clear();
    packedMat.value_.clear();    

    int idx(0);
    for(int c(0); c<emat.cols(); ++c){
        packedMat.start_.push_back(idx);    
        for(int r(0); r<emat.rows(); ++r){                    
            if( std::abs(emat(r,c)) > ZCE)
            {                
                packedMat.index_.push_back(r);
                packedMat.value_.push_back(emat(r,c));
                idx ++;
            }            
        }        
    }
    packedMat.start_.push_back(idx);    
}

void printOptimizerMsgs(const Highs& highs){
    const HighsLp& lp = highs.getLp();
    // Get the model status
    const HighsModelStatus& model_status = highs.getModelStatus();
    std::cout << "Model status: " << highs.modelStatusToString(model_status) << std::endl;
    //  
    // Get the solution information
    const HighsInfo& info = highs.getInfo();
    std::cout << "Simplex iteration count: " << info.simplex_iteration_count << std::endl;
    std::cout << "Objective function value: " << info.objective_function_value << std::endl;
    std::cout << "Primal  solution status: " << highs.solutionStatusToString(info.primal_solution_status) << std::endl;
    std::cout << "Dual    solution status: " << highs.solutionStatusToString(info.dual_solution_status) << std::endl;
    std::cout << "Basis: " << highs.basisValidityToString(info.basis_validity) << std::endl;
    const bool has_values = info.primal_solution_status;
    const bool has_duals = info.dual_solution_status;
    const bool has_basis = info.basis_validity;

    const HighsSolution& solution = highs.getSolution();
    const HighsBasis& basis = highs.getBasis();
  //
  // Report the primal and solution values and basis
  for (int col=0; col < lp.num_col_; col++) {
    std::cout << "Column " << col;
    if (has_values) std::cout << "; value = " << solution.col_value[col];
    if (has_duals) std::cout << "; dual = " << solution.col_dual[col];
    if (has_basis) std::cout << "; status: " << highs.basisStatusToString(basis.col_status[col]);
    std::cout << std::endl;
  }
  for (int row=0; row < lp.num_row_; row++) {
    std::cout << "Row    " << row;
    if (has_values) std::cout << "; value = " << solution.row_value[row];
    if (has_duals) std::cout << "; dual = " << solution.row_dual[row];
    if (has_basis) std::cout << "; status: " << highs.basisStatusToString(basis.row_status[row]);
    std::cout << std::endl;
  }
}




// solve LP, linprog(f,A,b)
// min f'x s.t. lb<Ax<=ub
// double linprog(const Eigen::VectorXd& f,
//                 const Eigen::MatrixXd& A,
//                 const Eigen::VectorXd& lb,
//                 const Eigen::VectorXd& ub,
//                 Eigen::VectorXd& x){

//     HighsModel model;
//     // 
//     model.lp_.num_col_ = 2;
//     model.lp_.num_row_ = 3;
//     model.lp_.sense_ = ObjSense::kMinimize;
//     model.lp_.offset_ = 3;
//     model.lp_.col_cost_ = {1.0, 1.0};
//     model.lp_.col_lower_ = {0.0, 1.0};
//     model.lp_.col_upper_ = {4.0, 1.0e30};
//     model.lp_.row_lower_ = {-1.0e30, 5.0, 6.0};
//     model.lp_.row_upper_ = {7.0, 15.0, 1.0e30};

//     double ret;
//     return ret;    
// }
// // solve LP, linprog(f,A,b)
// // min f'x s.t. lb<Ax<=ub, l<x<u
// double linprog(const Eigen::VectorXd& f,
//                 const Eigen::MatrixXd& A,
//                 const Eigen::VectorXd& lb,
//                 const Eigen::VectorXd& ub,
//                 const Eigen::VectorXd& l,
//                 const Eigen::VectorXd& u,
//                 Eigen::VectorXd& x){

//     HighsModel model;
//     model.lp_.sense_ = ObjSense::kMinimize;
//     model.lp_.col_cost_.resize(f.size()); 
//     Eigen::VectorXd::Map(&model.lp_.col_cost_[0], f.size()) = f;

//     // 
//     model.lp_.num_col_ = 2;
//     model.lp_.num_row_ = 3;
    
//     model.lp_.offset_ = 3;
//     model.lp_.col_cost_ = {1.0, 1.0};
//     model.lp_.col_lower_ = {0.0, 1.0};
//     model.lp_.col_upper_ = {4.0, 1.0e30};
//     model.lp_.row_lower_ = {-1.0e30, 5.0, 6.0};
//     model.lp_.row_upper_ = {7.0, 15.0, 1.0e30};

//     double ret;
//     return ret;    
// }

} //namespace rossy_utils
