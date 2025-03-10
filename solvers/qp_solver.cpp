#include "qp_solver.hpp"
// QuadProg
#include "rossy_utils/thirdparty/goldfarb/QuadProg++.hh"
#include "rossy_utils/io/io_utilities.hpp"
// highs qp solver
#include "lp_solver.hpp"
#include <Highs.h>
#include <cassert>
// osqp
#include <OsqpEigen/OsqpEigen.h>
// qpoases
#include <qpOASES.hpp>

// for benchmark
#include "rossy_utils/general/clock.hpp"

#define ZCE 1e-8

namespace rossy_utils {
const double INF = std::numeric_limits<float>::infinity();

double qpprog(Eigen::MatrixXd& G, Eigen::VectorXd& g0,
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

double qpprog(Eigen::MatrixXd& G, Eigen::VectorXd& g0,
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
double qpprogHiGHS(const Eigen::MatrixXd & Q, 
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

void EigenMatrix2Hessian(const Eigen::MatrixXd& H,
    HighsHessian& hessian){
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
}

float qpprogQPOASES(const Eigen::MatrixXf& Q, 
                    const Eigen::VectorXf& q,
                    const Eigen::MatrixXf& A, 
                    const Eigen::VectorXf& b,
                    Eigen::VectorXf& x) {


    int nV = Q.cols();    // number of variables
    int nC = A.rows();    // number of constraints

    // Setup QProblem object
    qpOASES::QProblem problem(nV, nC);

    // Set Options
    qpOASES::Options options;
    options.setToMPC();
    options.printLevel = qpOASES::PL_NONE;
    problem.setOptions(options);

    // Convert Eigen matrices to raw arrays
    qpOASES::real_t* H = const_cast<qpOASES::real_t*>(Q.data());
    qpOASES::real_t* g = const_cast<qpOASES::real_t*>(q.data());
    qpOASES::real_t* A_data = const_cast<qpOASES::real_t*>(A.data());
    qpOASES::real_t* ub = const_cast<qpOASES::real_t*>(b.data());
    qpOASES::real_t* lb = new qpOASES::real_t[nC];
    for(int i = 0; i < nC; i++) lb[i] = - 1e-17;

    // Solve QP
    int nWSR = 1000;
    qpOASES::returnValue status = problem.init(H, g, A_data, NULL, NULL, lb, ub, nWSR);

    // Get solution
    qpOASES::real_t* xOpt = new qpOASES::real_t[nV];
    problem.getPrimalSolution(xOpt);
    
    // Copy solution to Eigen vector
    x = Eigen::Map<Eigen::VectorXf>(xOpt, nV);

    delete[] lb;
    delete[] xOpt;

    return (status == qpOASES::SUCCESSFUL_RETURN) ? 0.0f : -1.0f;
}


OSQPSolver::OSQPSolver():
n_(0),m_(0) {
    solver_ = new OsqpEigen::Solver();
    // Set default solver settings
    solver_->settings()->setVerbosity(false);
    solver_->settings()->setWarmStart(true);
    solver_->settings()->setAbsoluteTolerance(1e-3);
    solver_->settings()->setRelativeTolerance(1e-3);
    solver_->settings()->setMaxIteration(500);
    solver_->settings()->setPolish(false);

    // Settings for faster convergence
    // solver_->settings()->setAlpha(1.8);
    // solver_->settings()->setRho(0.05);
    // solver_->settings()->setSigma(1e-6);
    // solver_->settings()->setPolish(true);
    // solver_->settings()->setPolishRefineIter(4);
}

OSQPSolver::~OSQPSolver() {
    if (solver_ != nullptr) {
        delete solver_;
        solver_ = nullptr;
    }
}

void OSQPSolver::saveProblem(
        const Eigen::SparseMatrix<float>& Q_sparse,
        const Eigen::VectorXf& q,
        const Eigen::SparseMatrix<float>& A_sparse,
        const Eigen::VectorXf& b){
    Eigen::MatrixXf Q_sparse_float = Eigen::MatrixXf(Q_sparse);
    saveMatrix(Q_sparse_float, "Q_sparse");
    saveVector(q, "q");
    Eigen::MatrixXf A_sparse_float = Eigen::MatrixXf(A_sparse);
    saveMatrix(A_sparse_float, "A_sparse");
    saveVector(b, "b");   
}


float OSQPSolver::qpprogOSQPSparse(
        Eigen::SparseMatrix<float>& Q_sparse,
        Eigen::VectorXf& q,
        Eigen::SparseMatrix<float>& A_sparse,
        Eigen::VectorXf& b,
        Eigen::VectorXf& x) {
    // Clock timer;

    // Problem dimensions
    int n = Q_sparse.cols();
    int m = A_sparse.rows();
        
    // Set problem dimensions
    solver_->data()->setNumberOfVariables(n);
    solver_->data()->setNumberOfConstraints(m);
    Eigen::VectorXf lower_bound = Eigen::VectorXf::Constant(m, -1e17f);
    // timer.printElapsedMiliSec("OSQP: setting = ");

    // Load data into solver 
    if(n_==n || m_==m || solver_->isInitialized()){
        // set triplet from sparse matrix
    A_triplets_old_ = A_triplets_new_;
    setTripletfromSparseMatrix(A_sparse, A_triplets_new_);

    if(isPatternChanged(A_triplets_old_, A_triplets_new_)){
        // std::cout <<"##### OSQP sparsity pattern changed!"  << std::endl;
        solver_->data()->clearLinearConstraintsMatrix();
        if(!solver_->data()->setLinearConstraintsMatrix(A_sparse) ||
                !solver_->data()->setLowerBound(lower_bound) ||
                !solver_->data()->setUpperBound(b) ||
                !solver_->data()->setGradient(q)){
            std::cout <<"##### OSQP set data failed!"  << std::endl;
            return -1.0f;
        } 

        solver_->clearSolver();
        if (!solver_->initSolver()) {
            std::cerr << " ##### OSQP Solver initialization failed!" << std::endl;
            return -1.0f;
        } 
    }        
    else if (!solver_->updateGradient(q) ||
            !solver_->updateUpperBound(b) ||
            !solver_->updateLinearConstraintsMatrix(A_sparse)) {
        // Assume Q_sparse and lower bound will be same
        std::cout <<"##### OSQP update data failed!"  << std::endl;
        return -1.0f; 
    }

    }
    else{
        if(solver_->data()->isSet()){
            solver_->data()->clearLinearConstraintsMatrix();
            solver_->data()->clearHessianMatrix();
            solver_->clearSolver();  
        }
        if (!solver_->data()->setHessianMatrix(Q_sparse) ||
            !solver_->data()->setGradient(q) ||
            !solver_->data()->setLinearConstraintsMatrix(A_sparse) ||
            !solver_->data()->setLowerBound(lower_bound) ||
            !solver_->data()->setUpperBound(b)) {
            std::cout <<"##### OSQP set data failed!"  << std::endl;
            return -1.0f; 
        }
        // Initialize solver   
        solver_->clearSolver();
        if (!solver_->initSolver()) {
            std::cerr << " ##### OSQP Solver initialization failed!" << std::endl;
            return -1.0f;
        }
        // timer.printElapsedMiliSec("OSQP: init solver = ");
        n_ = n;
        m_ = m;

        setTripletfromSparseMatrix(A_sparse, A_triplets_new_);
    }
        

    // Solve the QP problem
    if (solver_->solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
        return -1.0f;
    }
    // timer.printElapsedMiliSec("OSQP: solve = ");

    // Get the optimal solution (convert back to float)
    x = solver_->getSolution();
    // Compute and return the optimal objective function value (float precision)
    float ret = (0.5f * x.transpose() * Q_sparse * x + q.transpose() * x).value();
    // timer.printElapsedMiliSec("OSQP: get solution = ");
    return ret;
}

float OSQPSolver::qpprogOSQP(
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
    // timer.printElapsedMiliSec("OSQP: dense to sparse  = ");

    return qpprogOSQPSparse(Q_sparse, q, A_sparse, b, x);
}

void OSQPSolver::setTripletfromSparseMatrix(
          const Eigen::SparseMatrix<float>& A, 
          std::vector<Eigen::Triplet<float>> &A_triplets){
    A_triplets.clear();
    for(int i=0; i<A.outerSize(); ++i){
        for(Eigen::SparseMatrix<float>::InnerIterator it(A,i); it; ++it){
            A_triplets.emplace_back(it.row(), it.col(), it.value());
        }
    }
}

bool OSQPSolver::isPatternChanged(
          const std::vector<Eigen::Triplet<float>>& A_old,
          const std::vector<Eigen::Triplet<float>>& A_new){
    if(A_old.size() != A_new.size()){
        return true;
    }
    for (int i = 0; i < A_new.size(); i++)
    {
        // check if the sparsity pattern is changed
        if ((A_new[i].row() != A_old[i].row())
                || (A_new[i].col() != A_old[i].col()))
            return true;
    }
    return false;
}

} // namespace rossy_utils