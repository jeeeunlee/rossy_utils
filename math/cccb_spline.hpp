#ifndef CCCB_SPLINE_H
#define CCCB_SPLINE_H

#include <iostream>
#include <math.h>
#include <algorithm>
#include <Eigen/Dense>
#include <vector>
#include "rossy_utils/math/typedefs.h"

// CCCBspline: Clamped Cardinal Cubic B-spline interpolation

template <typename Scalar>
class CCCBSpline{
public:
	CCCBSpline();
	CCCBSpline(const Scalar &pi, const Scalar &pf, 
			const std::vector<Scalar> &cp_var_in, Scalar h_in=1.);
	CCCBSpline(const Scalar &pi, const Scalar &pf, 
			const VectorX<Scalar> &cp_var_in, Scalar h_in=1.);
	~CCCBSpline();

	int evaluateTimeInterval(const Scalar &t_in);
	Scalar evaluate(const Scalar &t_in);
	Scalar evaluateFirstDerivative(const Scalar &t_in);
	Scalar evaluateSecondDerivative(const Scalar &t_in);    

    void setTimeDuration(const Scalar &h_in);
	void setControlPoints(const Scalar &pi, 
						const Scalar &pf,
						const std::vector<Scalar> &cp_in);

	int getNumIntervals(){return N_;}
	Scalar getMotionPeriod(){return ((Scalar)N_) * h_;}

private:
	// p = m-n-1 = (N+6)-(N+2)-1 (Cubic spline with C2 continuity)
    static const int p_ = 3; // degree
    int N_; // # of intervals
	Scalar h_; // time duration

    // cp{-3}, cp{-2},..., cp{N-1}
    std::vector<Scalar> cp_; // (n+1)=N+3 control points
	// kv{-3}, cp{-2},..., kv{N+3}: N+7 knot vector
    
	int getNumCPs(){return (N_+2)+1; }
	int getNumKVs(){return (N_+6)+1; }	
};


template <typename Scalar>
class CCCBSplineVec{
public:
	CCCBSplineVec();
	CCCBSplineVec(const VectorX<Scalar> &pi, const VectorX<Scalar> &pf, 
			const std::vector<VectorX<Scalar>> &cp_var_in, Scalar h_in=1.);
	~CCCBSplineVec();

	void initialize();
	
	int evaluateTimeInterval(const Scalar &t_in);	
	VectorX<Scalar> evaluate(const Scalar &t_in);
	VectorX<Scalar> evaluateFirstDerivative(const Scalar &t_in);
	VectorX<Scalar> evaluateSecondDerivative(const Scalar &t_in);

	void setTimeDuration(const Scalar &h_in);
	void setControlPoints(const VectorX<Scalar> &pi, 
						const VectorX<Scalar> &pf,
						const std::vector<VectorX<Scalar>> &cp_in);
	void addControlPoints(const Scalar &pi, 
						const Scalar &pf, 
						const VectorX<Scalar> &cp_in);

	int getDim(){return dim_; }
	int getNumIntervals(){return N_;}
	Scalar getMotionPeriod(){return ((Scalar)N_) * h_;}

private:
	// p = m-n-1 = (N+6)-(N+2)-1 (Cubic spline with C2 continuity)
    static const int p_ = 3; // degree
    int N_; // # of intervals
	int dim_; // dim of vec
	Scalar h_; // time duration

	std::vector<CCCBSpline<Scalar>> curves_; // size = dim_
 	VectorX<Scalar> output_;	
	
	int getNumCPs(){return (N_+2)+1; }
	int getNumKVs(){return (N_+6)+1; }
	
};

#endif 