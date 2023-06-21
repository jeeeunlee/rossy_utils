
#include <math/hermite_quintic_splines.hpp>
#include <io/io_utilities.hpp>
#include <math.h>
#include <algorithm>

HermiteQuinticSplines::HermiteQuinticSplines(){
    initialize();
}

void HermiteQuinticSplines::initialize(){
    computed = false;
    ts.clear();
    ys.clear();   
    dys.clear();
    ddys.clear();
}

void HermiteQuinticSplines::push_back(double t, 
                                      double y,
                                      double ydot,
                                      double yddot){
    computed=false;
    ts.push_back(t);
    ys.push_back(y);
    dys.push_back(ydot);
    ddys.push_back(yddot);
}

// Hermite basis functions
// p(t) = p0
void HermiteQuinticSplines::compute(){
    // ts.size()=(n_wpts+1)
    n_wpts = ts.size()-1;

    v0s.clear();
    v1s.clear();
    a0s.clear();
    a1s.clear();

    double tk;
    for(int i(0); i<n_wpts; ++i){
        tk = ts[i+1] - ts[i];
        v0s.push_back(dys[i]*tk);
        v1s.push_back(dys[i+1]*tk);
        a0s.push_back(ddys[i]*tk*tk);
        a1s.push_back(ddys[i+1]*tk*tk);
    }
    

    computed=true;
}

double HermiteQuinticSplines::evaluate(const double & t_in){
    if(!computed) compute();
    int i = evaluateTimeInterval(t_in);
    if(i<0) return ys[0];    
    else if(i<n_wpts){
        double t = (t_in - ts[i])/(ts[i+1] - ts[i]);
        double h0 = 1. - 10.* t*t*t + 15.*t*t*t*t - 6.*t*t*t*t*t;
        double h1 = t - 6.*t*t*t + 8.*t*t*t*t - 3.*t*t*t*t*t;
        double h2 = 0.5*t*t - 1.5*t*t*t + 1.5*t*t*t*t - 0.5*t*t*t*t*t;
        double h3 = 0.5*t*t*t - 4.*t*t*t*t + 0.5*t*t*t*t*t;
        double h4 = -4.*t*t*t + 7.*t*t*t*t - 3.*t*t*t*t*t;
        double h5 = 10.*t*t*t - 15.*t*t*t*t + 6.*t*t*t*t*t;
        double s = ys[i]*h0 + v0s[i]*h1 + a0s[i]*h2 
            + a1s[i]*h3  + v1s[i]*h4 + ys[i+1]*h5;
        return s;        
    }
    else return ys[n_wpts];
}

double HermiteQuinticSplines::evaluateFirstDerivative(const double & t_in){
    if(!computed) compute();

    int i = evaluateTimeInterval(t_in);
    if(i<0) i=0;
    else if(i<n_wpts){}
    else i=n_wpts-1;

    double t = (t_in - ts[i])/(ts[i+1] - ts[i]);
    double sdot = 0.; // todo
    return sdot;     
}

double HermiteQuinticSplines::evaluateSecondDerivative(const double & t_in){
    if(!computed) compute();
    
    int i = evaluateTimeInterval(t_in);
    if(i<0) i=0;
    else if(i<n_wpts){}
    else i=n_wpts-1;
    double t = (t_in - ts[i])/(ts[i+1] - ts[i]);
    double sddot = 0.; // todo
    return sddot;
}

int HermiteQuinticSplines::evaluateTimeInterval(const double & t_in){    
    // if ts[i] < t_in < ts[i+1]: return i = 0 ~ (n_wpts-1)
    // if t_in<ts[0]: return i = -1, 
    // if ts[n_wpts]<t_in: return i = n_wpts 
    int i=-1;
    for(auto &ti : ts){ //ts.size()=(n_wpts+1)
        if(ti > t_in) break;
        else i++; // ti <= t_in        
    }
    return i;
}


/// -----------------------------------------

HQSpln4Vec::HQSpln4Vec(){
    initialize(0);
}

void HQSpln4Vec::initialize(int _dim){
    computed=false;

    dim = _dim;    
    curves.clear();
    for(int i(0); i<dim; ++i){
        curves.push_back( HermiteQuinticSplines() );
        curves[i].initialize();
    }
    output = Eigen::VectorXd::Zero(dim);
}
void HQSpln4Vec::push_back(double t, 
                  const Eigen::VectorXd& Y,
                  const Eigen::VectorXd& Ydot, 
                  const Eigen::VectorXd& Yddot){
    assert(Y.size() == dim);
    computed=false;    
    for(int i(0); i<dim; ++i){
        curves[i].push_back(t, Y(i), Ydot(i), Yddot(i));
    }
}
void HQSpln4Vec::compute(){
    computed = true;
    for(int i(0); i<dim; ++i)
        curves[i].compute();
}

Eigen::VectorXd HQSpln4Vec::evaluate(const double & t_in){
    output = Eigen::VectorXd::Zero(dim);
    for(int i(0); i<dim; ++i)
        output[i] = curves[i].evaluate(t_in);
    return output;
}   

Eigen::VectorXd HQSpln4Vec::evaluateFirstDerivative(const double & t_in){
    output = Eigen::VectorXd::Zero(dim);
    for(int i(0); i<dim; ++i)
        output[i] = curves[i].evaluateFirstDerivative(t_in);
    return output;
}

Eigen::VectorXd HQSpln4Vec::evaluateSecondDerivative(const double & t_in){
    output = Eigen::VectorXd::Zero(dim);
    for(int i(0); i<dim; ++i)
        output[i] = curves[i].evaluateSecondDerivative(t_in);
    return output;
}

