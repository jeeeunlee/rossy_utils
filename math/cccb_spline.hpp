#ifndef CCCB_SPLINE_H
#define CCCB_SPLINE_H

#include <iostream>
#include <math.h>
#include <algorithm>
#include <Eigen/Dense>
#include <vector>

// CCCBspline: Clamped Cardinal Cubic B-spline interpolation

class CCCBSpline{
public:
	CCCBSpline();
	~CCCBSpline();

	double evaluate(const double & t_in);
	double evaluateFirstDerivative(const double & t_in);
	double evaluateSecondDerivative(const double & t_in);
    int evaluateTimeInterval(const double & t_in);

    void setTimeDuration(const double &h_in);

private:
    static const int p_ = 3; // degree = m-n-1 (Cubic spline with C2 continuity)
    int N_; // # of intervals
    int n_; // # of control points -1  = N+3
    int m_; // # of knots -1 = N+7

    // cp{-3}, cp{-2},..., cp{N-1}
    std::vector<double> cp_; // (n+1) control points
    static const int cp_base_ = 3; // cp{i} = cp_[i + cp_base_]
    
    // kv{-3}, kv{-2},..., kv{N+3}
    std::vector<double> kv_; // (m+1) knot vectors    
    static const int kv_base_ = 3; // kv{i} = kv_[i + kv_base_] = i*h if 0<=i<=N

	double h_;

    void resetKnotVectors();    
};



class CCCBSplineVec{
public:
	CCCBSplineVec();
	CCCBSplineVec(const Eigen::VectorXd & start_pos, const Eigen::VectorXd & start_vel, 
				   const Eigen::VectorXd & end_pos, const Eigen::VectorXd & end_vel, const double & duration);
	~CCCBSplineVec();
	
	void initialize(const Eigen::VectorXd & start_pos, const Eigen::VectorXd & start_vel, 
					const Eigen::VectorXd & end_pos, const Eigen::VectorXd & end_vel, const double & duration);
	Eigen::VectorXd evaluate(const double & t_in);
	Eigen::VectorXd evaluateFirstDerivative(const double & t_in);
	Eigen::VectorXd evaluateSecondDerivative(const double & t_in);

private:


	std::vector<CCCBSpline> curves;
 	Eigen::VectorXd output;
};

#endif 