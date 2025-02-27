#include "qp_solver.hpp"
// QuadProg
#include "rossy_utils/thirdparty/goldfarb/QuadProg++.hh"
// highs qp solver
#include "lp_solver.hpp"
#include <Highs.h>
#include <cassert>
// osqp
#include <OsqpEigen/OsqpEigen.h>

// for benchmark
#include "rossy_utils/general/clock.hpp"




#define ZCE 1e-8
const double INF = std::numeric_limits<float>::infinity();

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

// f =   min   0.5 * x Q x + q x, where Q=Q'(symmetric, psd)
//       s.t.  _A x + <= _b 
// return f
double rossy_utils::qpprogHiGHS(const Eigen::MatrixXd & Q, 
                    const Eigen::VectorXd & q, 
                    const Eigen::MatrixXd & A, 
                    const Eigen::VectorXd & b, 
                    Eigen::VectorXd & x){

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

float rossy_utils::qpprogOSQPSparse(
        Eigen::SparseMatrix<float> &Q_sparse, 
        Eigen::VectorXf &q,
        Eigen::SparseMatrix<float> &A_sparse, 
        Eigen::VectorXf &b,
        Eigen::VectorXf &x) {
    int n = Q_sparse.rows(); // Number of variables
    int m = A_sparse.rows(); // Number of constraints
    x = Eigen::VectorXf::Zero(n);
    if (m != b.rows()) {
        std::cerr << "Error: A.rows() must match b.rows()!" << std::endl;
        return -1.0f;
    }
    Clock timer;
    timer.start();

    // Create an OSQP Solver instance
    OsqpEigen::Solver solver;

    // Set problem dimensions
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);
    // solver.settings()->setAbsoluteTolerance(1e-3);    // Relax tolerance
    // solver.settings()->setRelativeTolerance(1e-3);    // Relax tolerance
    // solver.settings()->setMaxIteraction(500);        // Limit iterations
    // solver.settings()->setAdaptiveRho(true);          // Enable adaptive penalty
    // solver.settings()->setPolish(false);              // Disable polishing
    // solver.settings()->setAdaptiveRhoInterval(25);    // Faster adaptation

    solver.data()->setNumberOfVariables(n);
    solver.data()->setNumberOfConstraints(m);
    timer.printElapsedMiliSec("OSQP: setting = ");

    // Define sparse matrices (OSQP requires sparse format)
    Eigen::SparseMatrix<double> Q_double = Q_sparse.cast<double>();
    Eigen::SparseMatrix<double> A_double = A_sparse.cast<double>();
    Eigen::VectorXd lower_bound = Eigen::VectorXd::Constant(m, -1e17);
    Eigen::VectorXd qdouble = q.cast<double>();
    Eigen::VectorXd bdouble = b.cast<double>();
    Eigen::VectorXd xdouble = x.cast<double>();
    timer.printElapsedMiliSec("OSQP: cast  = ");

    // Load data into solver   
    if (!solver.data()->setHessianMatrix(Q_double) ||
        !solver.data()->setGradient(qdouble) ||
        !solver.data()->setLinearConstraintsMatrix(A_double) ||
        !solver.data()->setLowerBound(lower_bound) ||
        !solver.data()->setUpperBound(bdouble)) {
        std::cout <<"##### OSQP set data failed!"  << std::endl;
        return -1.0f; 
    }

    // Initialize solver
    solver.clearSolver();
    if (!solver.initSolver()) {
        std::cerr << " ##### OSQP Solver initialization failed!" << std::endl;
        return -1.0f;
    }
    timer.printElapsedMiliSec("OSQP: init solver = ");

    // Solve the QP problem
    if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
        return -1.0f;
    }
    timer.printElapsedMiliSec("OSQP: solve = ");

    // Get the optimal solution (convert back to float)
    xdouble = solver.getSolution();
    x = xdouble.cast<float>();
    // Compute and return the optimal objective function value (float precision)
    float ret = (0.5f * x.transpose() * Q_sparse * x + q.transpose() * x).value();
    timer.printElapsedMiliSec("OSQP: get solution = ");
    return ret;
}



float rossy_utils::qpprogOSQP(
        Eigen::MatrixXf &Q, 
        Eigen::VectorXf &q,
        Eigen::MatrixXf &A, 
        Eigen::VectorXf &b,
        Eigen::VectorXf &x) {

    Clock timer;
    timer.start();
    Eigen::SparseMatrix<float> Q_sparse = Q.sparseView();
    Eigen::SparseMatrix<float> A_sparse = A.sparseView();
    Q_sparse.makeCompressed();
    A_sparse.makeCompressed();
    timer.printElapsedMiliSec("OSQP: dense to sparse  = ");

    return qpprogOSQPSparse(Q_sparse, q, A_sparse, b, x);
}



// float rossy_utils::qpprogOSQP(Eigen::MatrixXf &Q, Eigen::VectorXf &q,
//                  Eigen::MatrixXf &A, Eigen::VectorXf &b,
//                  Eigen::VectorXf &x) {
//     int n = Q.rows(); // Number of variables
//     int m = A.rows(); // Number of constraints
//     if (m != b.rows()) {
//         std::cerr << "Error: A.rows() must match b.rows()!" << std::endl;
//         return -1.0f;
//     }

//     // Create an OSQP Solver instance
//     OsqpEigen::Solver solver;
//     solver.clearSolver();

//     // Set problem dimensions
//     solver.settings()->setVerbosity(true);
//     solver.settings()->setWarmStart(false);
//     solver.data()->setNumberOfVariables(n);
//     solver.data()->setNumberOfConstraints(m);

//     std::cout <<" here 1 : m=" << m <<", n="<<n << std::endl;

//     // Define sparse matrices (OSQP requires sparse format)
//     Q = 0.5f * (Q + Q.transpose());

//     // Check if Q is PSD
//     Eigen::LLT<Eigen::MatrixXf> llt(Q);
//     if (llt.info() != Eigen::Success) {
//         std::cerr << "Error: Q is not positive semi-definite!" << std::endl;
//         std::cout<<"Q(8,8) = " << Q.topLeftCorner(8,8) << std::endl;
//         std::cout<<"Q(end,end) = " << Q.bottomRightCorner(1,1) << std::endl;
//         return -1.0f;   
//     }

//     // Check for NaN or Inf values
//     if (!Q.allFinite()) {
//         std::cerr << "Error: Q contains NaN or Inf values!" << std::endl;
//         return -1.0f;
//     }

//     // Ensure Q is square
//     if (Q.rows() != Q.cols()) {
//         std::cerr << "Error: Q must be square (n x n)!" << std::endl;
//         return -1.0f;
//     }

//     Eigen::SparseMatrix<float> Q_sparse = Q.sparseView();
//     Eigen::SparseMatrix<float> A_sparse = A.sparseView();
//     Eigen::VectorXf lower_bound = Eigen::VectorXf::Constant(m, -OSQP_INFTY);

//     std::cout <<" here 2 " << std::endl;

//     // Load data into solver
//     Q_sparse.makeCompressed();
//     A_sparse.makeCompressed();
//     if (!solver.data()->setHessianMatrix(Q_sparse)){
//         std::cout <<" setHessianMatrix "  << std::endl;
//         return -1.0f; 
//     }
//     std::cout <<" here 3 " << std::endl;
//     if(!solver.data()->setGradient(q)){
//         std::cout <<" setGradient "  << std::endl;
//         return -1.0f; 
//     }
//     std::cout <<" here 4 " << std::endl;
//     if(!solver.data()->setLinearConstraintsMatrix(A_sparse)){
//         std::cout <<" setLinearConstraintsMatrix "  << std::endl;
//         return -1.0f; 
//     }
//     std::cout <<" here 5 " << std::endl;
//     if(!solver.data()->setLowerBound(lower_bound)){
//         std::cout <<" setLowerBound "  << std::endl;
//         return -1.0f; 
//     }
//     std::cout <<" here 6 " << std::endl;
//     if(!solver.data()->setUpperBound(b)) {
//         std::cout <<" setUpperBound "  << std::endl;
//         return -1.0f; 
//     }
//     std::cout <<" here 7 " << std::endl;
//     std::cout << "is set = " << solver.data()->isSet() << std::endl;

//     std::cout << "m = " << solver.data()->getData()->m << std::endl;

//     // Initialize solver
//     if (!solver.initSolver()) {
//         std::cerr << "OSQP Solver initialization failed!" << std::endl;
//         return -1.0f;
//     }
//     std::cout <<" here 8 " << std::endl;

//     solver.clearSolver();
//     // Solve the QP problem
//     if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
//         return -1.0f;
//     }

//     // Get the optimal solution (convert back to float)
//     x = solver.getSolution();
    
//     // Compute and return the optimal objective function value (float precision)
//     float ret = (0.5f * x.transpose() * Q * x + q.transpose() * x).value();
//     return ret;
// }

