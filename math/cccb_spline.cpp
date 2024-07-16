#include "rossy_utils/math/cccb_spline.hpp"
#include "rossy_utils/math_utilities.hpp"

CCCBSpline::CCCBSpline()
{
}

CCCBSpline::~CCCBSpline()
{
    cp_.clear();
    kv_.clear();
}

double CCCBSpline::evaluate(const double &t_in)
{
    return 0.0;
}

double CCCBSpline::evaluateFirstDerivative(const double &t_in)
{
    return 0.0;
}

double CCCBSpline::evaluateSecondDerivative(const double &t_in)
{
    return 0.0;
}

int CCCBSpline::evaluateTimeInterval(const double &t_in)
{
    // if kv_[i] < t_in < kv_[i+1]: return i = 0 ~ N_-1
    int i = (int)(t_in / h_) + kv_base_;
    return i;
}


void CCCBSpline::setTimeDuration(const double &h_in)
{
    h_ = h_in;
}

void CCCBSpline::resetKnotVectors()
{
    // reset kvs based on h_;
    kv_.resize(n_+5);
    for(int i{0}; i<n_+5; ++i) {
        kv_[i] = rossy_utils::CropValue((double)(i-kv_base_), 0., N_)*h_;
    }
}


CCCBSplineVec::CCCBSplineVec()
{
}

CCCBSplineVec::CCCBSplineVec(const Eigen::VectorXd &start_pos, const Eigen::VectorXd &start_vel, const Eigen::VectorXd &end_pos, const Eigen::VectorXd &end_vel, const double &duration)
{
}

CCCBSplineVec::~CCCBSplineVec()
{
}

void CCCBSplineVec::initialize(const Eigen::VectorXd &start_pos, const Eigen::VectorXd &start_vel, const Eigen::VectorXd &end_pos, const Eigen::VectorXd &end_vel, const double &duration)
{
}

Eigen::VectorXd CCCBSplineVec::evaluate(const double &t_in)
{
    return Eigen::VectorXd();
}

Eigen::VectorXd CCCBSplineVec::evaluateFirstDerivative(const double &t_in)
{
    return Eigen::VectorXd();
}

Eigen::VectorXd CCCBSplineVec::evaluateSecondDerivative(const double &t_in)
{
    return Eigen::VectorXd();
}
