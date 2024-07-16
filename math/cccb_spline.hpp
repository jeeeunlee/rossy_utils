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
	CCCBSpline(const double &pi, const double &pf, 
			const std::vector<double> &cp_var_in, double h_in=1.);
	~CCCBSpline();

	int evaluateTimeInterval(const double &t_in);
	double evaluate(const double &t_in);
	double evaluateFirstDerivative(const double &t_in);
	double evaluateSecondDerivative(const double &t_in);    

    void setTimeDuration(const double &h_in);
	void setControlPoints(const double &pi, 
						const double &pf,
						const std::vector<double> &cp_in);
	void resetKnotVectors();
	int getNumIntervals(){return N_;}

private:
	// p = m-n-1 = (N+6)-(N+2)-1 (Cubic spline with C2 continuity)
    static const int p_ = 3; // degree
    int N_; // # of intervals
	double h_; // time duration
    // cp{-3}, cp{-2},..., cp{N-1}
    std::vector<double> cp_; // (n+1)=N+3 control points    
    // kv{-3}, kv{-2},..., kv{N+3}
    std::vector<double> kv_; // (m+1)=N+7 knot vectors	
    
	int getNumCPs(){return (N_+2)+1; }
	int getNumKVs(){return (N_+6)+1; }	
};


class CCCBSplineVec{
public:
	CCCBSplineVec();
	CCCBSplineVec(const Eigen::VectorXd &pi, const Eigen::VectorXd &pf, 
			const std::vector<Eigen::VectorXd> &cp_var_in, double h_in=1.);
	~CCCBSplineVec();
	
	int evaluateTimeInterval(const double &t_in);	
	Eigen::VectorXd evaluate(const double &t_in);
	Eigen::VectorXd evaluateFirstDerivative(const double &t_in);
	Eigen::VectorXd evaluateSecondDerivative(const double &t_in);

	void setTimeDuration(const double &h_in);
	void setControlPoints(const Eigen::VectorXd &pi, 
						const Eigen::VectorXd &pf,
						const std::vector<Eigen::VectorXd> &cp_in);
	void resetKnotVectors();

	int getDim(){return dim_; }
	int getNumIntervals(){return N_;}

private:
	// p = m-n-1 = (N+6)-(N+2)-1 (Cubic spline with C2 continuity)
    static const int p_ = 3; // degree
    int N_; // # of intervals
	int dim_; // dim of vec
	double h_; // time duration

	std::vector<CCCBSpline> curves_; // size = dim_
 	Eigen::VectorXd output_;	
	
	int getNumCPs(){return (N_+2)+1; }
	int getNumKVs(){return (N_+6)+1; }
	
};

#endif 