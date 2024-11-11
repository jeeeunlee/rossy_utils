#include "rossy_utils/math/cccb_spline.hpp"
#include "rossy_utils/math/math_utilities.hpp"
#include "cccb_spline.hpp"


template <typename Scalar>
CCCBSpline<Scalar>::CCCBSpline()
{
    cp_.clear();
    h_ = 1.;
    N_ = 0;
}

template <typename Scalar>
CCCBSpline<Scalar>::CCCBSpline(const Scalar& pi, 
        const Scalar& pf, 
        const std::vector<Scalar>& cp_var_in,
        Scalar h_in){
    h_ = h_in;
    setControlPoints(pi, pf, cp_var_in);     
}


template <typename Scalar>
CCCBSpline<Scalar>::CCCBSpline(const Scalar& pi, 
        const Scalar& pf, 
        const VectorX<Scalar>& cp_var_in,
        Scalar h_in){
    h_ = h_in;

    // Eigen Vector 2 std vector
    std::vector<Scalar> cp_var;
    cp_var.resize(cp_var_in.size());
    VectorX<Scalar>::Map(&cp_var[0], cp_var_in.size()) = cp_var_in;
    setControlPoints(pi, pf, cp_var);
}

template <typename Scalar>
CCCBSpline<Scalar>::~CCCBSpline()
{
    cp_.clear();
}

template <typename Scalar>
Scalar CCCBSpline<Scalar>::evaluate(const Scalar &t_in)
{
    int i = evaluateTimeInterval(t_in);
    Scalar u;
    if (i<0) {
        i=0;
        u=0.;
    }
    else if (i>N_-1) {
        i = N_-1;
        u=1.;
    }
    else {
        u = t_in/h_ - (Scalar)i;
    }
    Scalar cp_1=cp_[i],cp_2=cp_[i+1],cp_3=cp_[i+2],cp_4=cp_[i+3];    
    Scalar u2=u*u, u3=u*u*u;
    Scalar ct = cp_1*(-u3/6. + u2/2. - u/2. + 1./6.)
                + cp_2*(u3/2. - u2 + 2./3.) 
                + cp_3*(-u3/2. + u2/2. + u/2. + 1./6.)
                + cp_4*(u3/6.);    
    return ct;
}

template <typename Scalar>
Scalar CCCBSpline<Scalar>::evaluateFirstDerivative(const Scalar &t_in)
{
    int i = evaluateTimeInterval(t_in);
    Scalar u;
    if (i<0) {
        i=0; u=0.;
    }
    else if (i>N_-1) {
        i = N_-1; u=1.;
    }
    else {
        u = t_in/h_ - (Scalar)i;
    }
    Scalar cp_1=cp_[i],cp_2=cp_[i+1],cp_3=cp_[i+2],cp_4=cp_[i+3];    
    Scalar u2=u*u;
    Scalar dct = (  cp_1*(-u2 + 2.*u - 1.) 
                + cp_2*( 3.*u2 - 4.*u) 
                + cp_3*(- 3.*u2 + 2.*u + 1.) 
                + cp_4*u2 ) /(2*h_);                  
    return dct;
}

template <typename Scalar>
Scalar CCCBSpline<Scalar>::evaluateSecondDerivative(const Scalar &t_in)
{
    int i = evaluateTimeInterval(t_in);
    Scalar u;
    if (i<0) {
        i=0; u=0.;
    }
    else if (i>N_-1) {
        i = N_-1; u=1.;
    }
    else {
        u = t_in/h_ - (Scalar)i;
    }
    Scalar cp_1=cp_[i],cp_2=cp_[i+1],cp_3=cp_[i+2],cp_4=cp_[i+3];    
    Scalar ddct = (  cp_1*(-u+1) 
    + cp_2*(3*u - 2) 
    + cp_3*(-3*u + 1) 
    + cp_4*u) / h_/h_ ;            
    return ddct;
}

template <typename Scalar>
int CCCBSpline<Scalar>::evaluateTimeInterval(const Scalar &t_in)
{
    // if i*h < t_in < (i+1)*h: return i = 0 ~ N_-1
    int i = (int)(t_in / h_); // 
    return i;
}


template <typename Scalar>
void CCCBSpline<Scalar>::setTimeDuration(const Scalar &h_in)
{ h_ = h_in; }

template <typename Scalar>
void CCCBSpline<Scalar>::setControlPoints(const Scalar & pi, 
    const Scalar & pf, 
    const std::vector<Scalar>& cp_in){
    // cp_var_in: cp1, cp2,...,
    // cp_ for cccbspline: pi,pi,pi, cp1,cp2,..., pf,pf,pf
    cp_.clear();
    cp_ = cp_in;
    cp_.insert(cp_.begin(), 3, pi);
    cp_.insert(cp_.end(), 3, pf);
    N_ = cp_.size()-3;
}

// Vector
template <typename Scalar>
CCCBSplineVec<Scalar>::CCCBSplineVec()
{
    initialize();
}

template <typename Scalar>
CCCBSplineVec<Scalar>::CCCBSplineVec(const VectorX<Scalar> & pi, 
    const VectorX<Scalar> & pf, 
    const std::vector<VectorX<Scalar>>& cp_var_in, 
    Scalar h_in){
    h_ = h_in;
    setControlPoints(pi, pf, cp_var_in); 
}


template <typename Scalar>
CCCBSplineVec<Scalar>::~CCCBSplineVec()
{
    curves_.clear();
}

template <typename Scalar>
void CCCBSplineVec<Scalar>::initialize(){
    curves_.clear();
    h_ = 1.;
    N_ = 0;
    dim_ = 0;
}


template <typename Scalar>
int CCCBSplineVec<Scalar>::evaluateTimeInterval(const Scalar & t_in)
{
    // if kv_[i] < t_in < kv_[i+1]: return i = 0 ~ N_-1
    int i = (int)(t_in / h_); // + kv_base_
    return i;
}

template <typename Scalar>
VectorX<Scalar> CCCBSplineVec<Scalar>::evaluate(const Scalar &t_in)
{
    output_ = VectorX<Scalar>::Zero(dim_);
    for(int d(0); d<dim_; ++d){
        output_[d] = curves_[d].evaluate(t_in);
    }
    return output_;
}

template <typename Scalar>
VectorX<Scalar> CCCBSplineVec<Scalar>::evaluateFirstDerivative(const Scalar &t_in)
{
    output_ = VectorX<Scalar>::Zero(dim_);
    for(int d(0); d<dim_; ++d){
        output_[d] = curves_[d].evaluateFirstDerivative(t_in);
    }
    return output_;
}

template <typename Scalar>
VectorX<Scalar> CCCBSplineVec<Scalar>::evaluateSecondDerivative(const Scalar &t_in)
{
    output_ = VectorX<Scalar>::Zero(dim_);
    for(int d(0); d<dim_; ++d){
        output_[d] = curves_[d].evaluateSecondDerivative(t_in);
    }
    return output_;
}

template <typename Scalar>
void CCCBSplineVec<Scalar>::setTimeDuration(const Scalar & h_in)
{
    h_ = h_in;
    for(auto & curve: curves_)
        curve.setTimeDuration(h_in);
}



template <typename Scalar>
void CCCBSplineVec<Scalar>::addControlPoints(const Scalar &pi, 
        const Scalar &pf, const VectorX<Scalar> &cp_in){
    
    curves_.push_back(CCCBSpline(pi,pf,cp_in));
    N_ = curves_[curves_.size()-1].getNumIntervals();
}

template <typename Scalar>
void CCCBSplineVec<Scalar>::setControlPoints(const VectorX<Scalar> & pi, 
    const VectorX<Scalar> & pf, 
    const std::vector<VectorX<Scalar>>& cp_in){
    assert(pi.size() == pf.size());
    dim_ = pi.size();

    // todo : check size
    if(curves_.size()!=dim_){
        curves_.clear();
        for(int d(0); d<dim_; ++d){
            curves_.push_back(CCCBSpline<Scalar>());
        }
    }
    std::vector<Scalar> cp_in_d;
    for(int d(0); d<dim_; ++d){
        cp_in_d.clear();
        for(auto & cp : cp_in)
            cp_in_d.push_back(cp[d]);
        curves_[d].setControlPoints(pi[d], pf[d], cp_in_d);
        N_ = curves_[d].getNumIntervals();
    }    
}

template class CCCBSpline<double>;
template class CCCBSpline<float>;
template class CCCBSplineVec<double>;
template class CCCBSplineVec<float>;