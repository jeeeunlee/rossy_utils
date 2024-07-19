#include "rossy_utils/math/cccb_spline.hpp"
#include "rossy_utils/math/math_utilities.hpp"
#include "cccb_spline.hpp"

CCCBSpline::CCCBSpline()
{
    cp_.clear();
    h_ = 1.;
    N_ = 0;
}

CCCBSpline::CCCBSpline(const double& pi, 
        const double& pf, 
        const std::vector<double>& cp_var_in,
        double h_in){
    h_ = h_in;
    setControlPoints(pi, pf, cp_var_in);     
}

CCCBSpline::CCCBSpline(const double& pi, 
        const double& pf, 
        const Eigen::VectorXd& cp_var_in,
        double h_in){
    h_ = h_in;

    // Eigen Vector 2 std vector
    std::vector<double> cp_var;
    cp_var.resize(cp_var_in.size());
    Eigen::VectorXd::Map(&cp_var[0], cp_var_in.size()) = cp_var_in;
    setControlPoints(pi, pf, cp_var);
}

CCCBSpline::~CCCBSpline()
{
    cp_.clear();
}

double CCCBSpline::evaluate(const double &t_in)
{
    int i = evaluateTimeInterval(t_in);
    double u;
    if (i<0) {
        i=0;
        u=0.;
    }
    else if (i>N_-1) {
        i = N_-1;
        u=1.;
    }
    else {
        u = t_in/h_ - (double)i;
    }
    double cp_1=cp_[i],cp_2=cp_[i+1],cp_3=cp_[i+2],cp_4=cp_[i+3];    
    double u2=u*u, u3=u*u*u;
    double ct = cp_1*(-u3/6. + u2/2. - u/2. + 1./6.)
                + cp_2*(u3/2. - u2 + 2./3.) 
                + cp_3*(-u3/2. + u2/2. + u/2. + 1./6.)
                + cp_4*(u3/6.);    
    return ct;
}

double CCCBSpline::evaluateFirstDerivative(const double &t_in)
{
    int i = evaluateTimeInterval(t_in);
    double u;
    if (i<0) {
        i=0; u=0.;
    }
    else if (i>N_-1) {
        i = N_-1; u=1.;
    }
    else {
        u = t_in/h_ - (double)i;
    }
    double cp_1=cp_[i],cp_2=cp_[i+1],cp_3=cp_[i+2],cp_4=cp_[i+3];    
    double u2=u*u;
    double dct = (  cp_1*(-u2 + 2.*u - 1.) 
                + cp_2*( 3.*u2 - 4.*u) 
                + cp_3*(- 3.*u2 + 2.*u + 1.) 
                + cp_4*u2 ) /(2*h_);                  
    return dct;
}

double CCCBSpline::evaluateSecondDerivative(const double &t_in)
{
    int i = evaluateTimeInterval(t_in);
    double u;
    if (i<0) {
        i=0; u=0.;
    }
    else if (i>N_-1) {
        i = N_-1; u=1.;
    }
    else {
        u = t_in/h_ - (double)i;
    }
    double cp_1=cp_[i],cp_2=cp_[i+1],cp_3=cp_[i+2],cp_4=cp_[i+3];    
    double ddct = (  cp_1*(-u+1) 
    + cp_2*(3*u - 2) 
    + cp_3*(-3*u + 1) 
    + cp_4*u) / h_/h_ ;            
    return ddct;
}

int CCCBSpline::evaluateTimeInterval(const double &t_in)
{
    // if i*h < t_in < (i+1)*h: return i = 0 ~ N_-1
    int i = (int)(t_in / h_); // 
    return i;
}


void CCCBSpline::setTimeDuration(const double &h_in)
{
    h_ = h_in;
}

void CCCBSpline::setControlPoints(const double & pi, 
    const double & pf, 
    const std::vector<double>& cp_in){
    // cp_var_in: cp1, cp2,...,
    // cp_ for cccbspline: pi,pi,pi, cp1,cp2,..., pf,pf,pf
    cp_.clear();
    cp_ = cp_in;
    cp_.insert(cp_.begin(), 3, pi);
    cp_.insert(cp_.end(), 3, pf);
    N_ = cp_.size()-3;
}

// Vector
CCCBSplineVec::CCCBSplineVec()
{
    initialize();
}

CCCBSplineVec::CCCBSplineVec(const Eigen::VectorXd & pi, 
    const Eigen::VectorXd & pf, 
    const std::vector<Eigen::VectorXd>& cp_var_in, 
    double h_in){
    h_ = h_in;
    setControlPoints(pi, pf, cp_var_in); 
}


CCCBSplineVec::~CCCBSplineVec()
{
    curves_.clear();
}

void CCCBSplineVec::initialize(){
    curves_.clear();
    h_ = 1.;
    N_ = 0;
    dim_ = 0;
}


int CCCBSplineVec::evaluateTimeInterval(const double & t_in)
{
    // if kv_[i] < t_in < kv_[i+1]: return i = 0 ~ N_-1
    int i = (int)(t_in / h_); // + kv_base_
    return i;
}

Eigen::VectorXd CCCBSplineVec::evaluate(const double &t_in)
{
    output_ = Eigen::VectorXd::Zero(dim_);
    for(int d(0); d<dim_; ++d){
        output_[d] = curves_[d].evaluate(t_in);
    }
    return output_;
}

Eigen::VectorXd CCCBSplineVec::evaluateFirstDerivative(const double &t_in)
{
    output_ = Eigen::VectorXd::Zero(dim_);
    for(int d(0); d<dim_; ++d){
        output_[d] = curves_[d].evaluateFirstDerivative(t_in);
    }
    return output_;
}

Eigen::VectorXd CCCBSplineVec::evaluateSecondDerivative(const double &t_in)
{
    output_ = Eigen::VectorXd::Zero(dim_);
    for(int d(0); d<dim_; ++d){
        output_[d] = curves_[d].evaluateSecondDerivative(t_in);
    }
    return output_;
}

void CCCBSplineVec::setTimeDuration(const double & h_in)
{
    h_ = h_in;
    for(auto & curve: curves_)
        curve.setTimeDuration(h_in);
}



void CCCBSplineVec::addControlPoints(const double &pi, 
        const double &pf, const Eigen::VectorXd &cp_in){
    
    curves_.push_back(CCCBSpline(pi,pf,cp_in));
    N_ = curves_[curves_.size()-1].getNumIntervals();
}

void CCCBSplineVec::setControlPoints(const Eigen::VectorXd & pi, 
    const Eigen::VectorXd & pf, 
    const std::vector<Eigen::VectorXd>& cp_in){
    assert(pi.size() == pf.size());
    dim_ = pi.size();

    // todo : check size
    if(curves_.size()!=dim_){
        curves_.clear();
        for(int d(0); d<dim_; ++d){
            curves_.push_back(CCCBSpline());
        }
    }
    std::vector<double> cp_in_d;
    for(int d(0); d<dim_; ++d){
        cp_in_d.clear();
        for(auto & cp : cp_in)
            cp_in_d.push_back(cp[d]);
        curves_[d].setControlPoints(pi[d], pf[d], cp_in_d);
        N_ = curves_[d].getNumIntervals();
    }    
}
