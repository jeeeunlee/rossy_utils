#pragma once

#include <stdio.h>
#include <Eigen/Dense>
#include "rossy_utils/io/io_utilities.hpp"
#include "rossy_utils/math/typedefs.h"
#include <iostream>


namespace rossy_utils {

// template <typename Scalar>
// using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
// template <typename Scalar>
// using VectorX = Eigen::Vector<Scalar, Eigen::Dynamic>;
// template <typename Scalar>
// using Vector3 = Eigen::Vector<Scalar, 3>;
// template <typename Scalar>
// using Quaternion = Eigen::Quaternion<Scalar>;
// template <typename Scalar>
// using AngleAxis = Eigen::AngleAxis<Scalar>;
// template <typename Scalar>
// using Isometry3 = Eigen::Transform<Scalar,3,Eigen::Isometry>;


template <typename Scalar>
MatrixX<Scalar> skew(const Vector3<Scalar>& w){
    // [[ 0, -3,  2],
    // [ 3,  0, -1],
    // [-2,  1,  0]]
    MatrixX<Scalar> Wx = MatrixX<Scalar>::Zero(3,3);
    Wx <<   0.0,   -w(2),  w(1),
           w(2),     0.0, -w(0),
          -w(1),    w(0),  0.0;
    return Wx;
};

template <typename Scalar>
MatrixX<Scalar> hStack(const MatrixX<Scalar>& a, const MatrixX<Scalar>& b) {
    if (a.rows()==0 || a.cols()==0)
        return b;
    if (b.rows()==0 || b.cols()==0)
        return a;
    if (a.rows() != b.rows()) {
        std::cout << "[hStack] Matrix Size is Wrong" << std::endl;
        exit(0);
    }

    MatrixX<Scalar> ab = MatrixX<Scalar>::Zero(a.rows(), a.cols() + b.cols());
    ab << a, b;
    return ab;
};

template <typename Scalar>
MatrixX<Scalar> hStack(const VectorX<Scalar>& a, const VectorX<Scalar>& b) {
    if (a.size()==0)
        return b;
    if (b.size()==0)
        return a;
    if (a.size() != b.size()) {
        std::cout << "[hStack] Vector Size is Wrong" << std::endl;
        exit(0);
    }
    MatrixX<Scalar> ab = MatrixX<Scalar>::Zero(a.size(), 2);
    ab << a, b;
    return ab;
};

template <typename Scalar>
MatrixX<Scalar> vStack(const MatrixX<Scalar>& a, const MatrixX<Scalar>& b) {
    if (a.rows()==0 || a.cols()==0)
        return b;
    if (b.rows()==0 || b.cols()==0)
        return a;
    if (a.cols() != b.cols()) {
        std::cout << "[vStack] Matrix Size is Wrong" << std::endl;
        exit(0);
    }
    MatrixX<Scalar> ab = MatrixX<Scalar>::Zero(a.rows() + b.rows(), a.cols());
    ab << a, b;
    return ab;
};

template <typename Scalar>
VectorX<Scalar> vStack(const VectorX<Scalar>& a, 
    const Eigen::DenseBase<MatrixX<Scalar>>& b) {
    return vStack((MatrixX<Scalar> const &)a,
                (MatrixX<Scalar> const &)b );
};

template <typename Scalar>
VectorX<Scalar> vStack(const VectorX<Scalar>& a, const VectorX<Scalar>& b) {
    if(a.size()==0) return b;
    if(b.size()==0) return a;
    VectorX<Scalar> ab = VectorX<Scalar>::Zero(a.size()+b.size());    
    // ab.head( a.size() ) = a; 
    // ab.tail( b.size() ) = b;
    ab << a, b;
    return ab;
};

template <typename Scalar>
VectorX<Scalar> vStack(
    const VectorX<Scalar>& a_, Scalar b_) {
    if(a_.size()==0) return VectorX<Scalar>::Constant(1, b_);
    VectorX<Scalar> ab = VectorX<Scalar>::Zero(a_.size()+1);    
    ab.head( a_.size() ) = a_; 
    ab(a_.size()) = b_;
    return ab;
};

template <typename Scalar>
MatrixX<Scalar> dStack(const MatrixX<Scalar>& a, const MatrixX<Scalar>& b) {
    // diagonally stack a,b -> [a 0; 0 b]
    if (a.rows()==0 || a.cols()==0)
        return b;
    if (b.rows()==0 || b.cols()==0)
        return a;
    MatrixX<Scalar> ab = MatrixX<Scalar>::Zero(a.rows() + b.rows(), a.cols() + b.cols());
    ab << a, MatrixX<Scalar>::Zero(a.rows(), b.cols()), 
        MatrixX<Scalar>::Zero(b.rows(), a.cols()), b;
    return ab;
};

template <typename Scalar>
VectorX<Scalar> MatrixtoVector(const MatrixX<Scalar>& a){
    // VectorX<Scalar> vec = VectorX<Scalar>::Zero(0);
    // for(int i(0); i<a.cols(); ++i){
    //     vec = vStack(vec, (MatrixX<Scalar>)(a.col(i)));
    // }
    // return vec;
    // colwise
    return Eigen::Map<const Eigen::VectorX<Scalar>>(a.data(), a.size());
};

template <typename Scalar>
MatrixX<Scalar> VectortoMatrix(const VectorX<Scalar>& a, int dim){
    int n=a.size();
    if( n%dim > 0) return a;
    MatrixX<Scalar> mat = MatrixX<Scalar>::Zero(dim , n/dim);
    for(int i(0); i<mat.cols(); ++i){
        mat.col(i) = a.segment(dim*i, dim);
    }
    return mat;
};

template <typename Scalar>
VectorX<Scalar> vector2EigenVector(
    const std::vector<Scalar>& vec, int k, int l){
    // from k, length l
    if(vec.size() < k+l){
        l = vec.size()-k;
        l = std::max(0, l);
    }
    VectorX<Scalar> ret = VectorX<Scalar>::Zero(l);
    for(int i(0); i<l; i++)
        ret(i) = vec[k+i];
    return ret;
};

template <typename Scalar>
VectorX<Scalar> vector2EigenVector(
    const std::vector<Scalar>& vec){
    return vector2EigenVector(vec,0,vec.size());
};

template <typename Scalar>
MatrixX<Scalar> vector2EigenMatrix(
        const std::vector<VectorX<Scalar>> &vOfv){
    int n = vOfv.size();
    if(n==0)
        return MatrixX<Scalar>::Zero(0,0);
    int dim = vOfv[0].size();
    MatrixX<Scalar> mat = MatrixX<Scalar>::Zero(dim,n); // dim x n
    for(int i(0); i<n; ++i){
        mat.col(i) = vOfv[i];
    }        
    return mat;
};

template <typename Scalar>
MatrixX<Scalar> kroneckerProduct(const MatrixX<Scalar> & A, const MatrixX<Scalar> & B)
{
    // Kronecker product A(m x n)⊗B(p x q) = C(pm x qn)
    // e.g. A⊗B = [ a11 B   a12 B   ... ]
    //            [ a21 B   a22 B   ... ]
    MatrixX<Scalar> ret = 
        MatrixX<Scalar>::Zero(B.rows()*A.rows(), B.cols()*A.cols());
    for(int i(0); i<A.rows(); ++i){
        for(int j(0); j<A.cols(); ++j){
            ret.block(i*B.rows(),j*B.cols(),B.rows(),B.cols()) = A(i,j)*B;
        }
    }
    return ret;
};

template <typename Scalar>
VectorX<Scalar> elementWiseDivisionExt(const VectorX<Scalar>& a, const VectorX<Scalar>& b){
    assert(a.size()%b.size()==0 && a.size()/b.size()>0);
    VectorX<Scalar> ret = VectorX<Scalar>::Zero(a.size());
    for(int i(0); i<a.size()/b.size();++i){
        ret.segment(i*b.size(),b.size()) = a.segment(i*b.size(),b.size()).array() / b.array();
    }
    return ret;
};

template <typename Scalar>
MatrixX<Scalar> deleteRow(const MatrixX<Scalar>& a_, int row_) {
    MatrixX<Scalar> ret = MatrixX<Scalar>::Zero(a_.rows() - 1, a_.cols());
    ret.block(0, 0, row_, a_.cols()) = a_.block(0, 0, row_, a_.cols());
    ret.block(row_, 0, ret.rows() - row_, a_.cols()) =
        a_.block(row_ + 1, 0, ret.rows() - row_, a_.cols());
    return ret;
};

template <typename Scalar>
Scalar getMaxRatioValue(
    const VectorX<Scalar>& val, 
    const VectorX<Scalar>& max){
    assert(val.size()==max.size());
    Scalar maxratio = 0.;
    for(int i(0); i<val.size(); ++i){
        maxratio = std::max(maxratio, abs(val(i)/max(i)));
    }
    return maxratio;
};  

template <typename Scalar>
Vector3<Scalar> convertQuatToExp(const Quaternion<Scalar>& q){
    AngleAxis aa = AngleAxis(q);
    return aa.axis() * aa.angle(); 
};

template <typename Scalar>
Vector3<Scalar> convertQuatToEulerAngles(const Quaternion<Scalar>& q) {
    Vector3<Scalar> angles;    //yaw pitch roll
    const auto x = q.x();
    const auto y = q.y();
    const auto z = q.z();
    const auto w = q.w();

    // roll (x-axis rotation)
    Scalar sinr_cosp = 2 * (w * x + y * z);
    Scalar cosr_cosp = 1 - 2 * (x * x + y * y);
    angles[2] = std::atan2(sinr_cosp, cosr_cosp);

    // pitch (y-axis rotation)
    Scalar sinp = 2 * (w * y - z * x);
    if (std::abs(sinp) >= 1)
        angles[1] = std::copysign(M_PI / 2, sinp); // use 90 degrees if out of range
    else
        angles[1] = std::asin(sinp);

    // yaw (z-axis rotation)
    Scalar siny_cosp = 2 * (w * z + x * y);
    Scalar cosy_cosp = 1 - 2 * (y * y + z * z);
    angles[0] = std::atan2(siny_cosp, cosy_cosp);
    return angles;
};

template <typename Scalar>
VectorX<Scalar> convertQuatDesToOriDes(const Quaternion<Scalar>& quat_in) {
  VectorX<Scalar> ori_out = VectorX<Scalar>::Zero(4);
  ori_out[0] = quat_in.w();
  ori_out[1] = quat_in.x();
  ori_out[2] = quat_in.y();
  ori_out[3] = quat_in.z();
  return ori_out;
};

template <typename Scalar>
void convertQuatDesToOriDes(const Quaternion<Scalar>& quat_in,  
                            VectorX<Scalar>& ori_out) {
  ori_out = VectorX<Scalar>::Zero(4);
  ori_out[0] = quat_in.w();
  ori_out[1] = quat_in.x();
  ori_out[2] = quat_in.y();
  ori_out[3] = quat_in.z();
};

template <typename Scalar>
void convertIsoToVec6d(const Isometry3<Scalar>& iso_in,
                    VectorX<Scalar>& vec_out){
    // Vec6d = (w,v) in se(3): 
    vec_out = VectorX<Scalar>::Zero(6);
    vec_out.tail(3) = iso_in.translation();
    auto ori = Eigen::Quaternion<Scalar>(iso_in.linear());   
    vec_out.head(3)=convertQuatToExp(ori);
};

template <typename Scalar>
void convertIsoToVec7d(const Isometry3<Scalar>& iso_in,
                    VectorX<Scalar>& vec_out){
    // Vec7d = (x,y,z, qw,qx,qy,qz) in R7
    vec_out = VectorX<Scalar>::Zero(7);
    vec_out.head(3) = iso_in.translation();
    auto ori = Eigen::Quaternion<Scalar>(iso_in.linear());   
    vec_out.tail(4)=convertQuatDesToOriDes(ori);
};

template <typename Scalar>
Scalar smooth_changing(Scalar ini, Scalar end, Scalar moving_duration,
                       Scalar curr_time) {
  Scalar ret;
  ret = ini + (end - ini) * 0.5 * (1 - cos(curr_time / moving_duration * M_PI));
  if (curr_time > moving_duration) {
    ret = end;
  }
  return ret;
};

template <typename Scalar>
Scalar smooth_changing_vel(Scalar ini, Scalar end, Scalar moving_duration,
                           Scalar curr_time) {
  Scalar ret;
  ret = (end - ini) * 0.5 * (M_PI / moving_duration) *
        sin(curr_time / moving_duration * M_PI);
  if (curr_time > moving_duration) {
    ret = 0.0;
  }
  return ret;
};

template <typename Scalar>
Scalar smooth_changing_acc(Scalar ini, Scalar end, Scalar moving_duration,
                           Scalar curr_time) {
  Scalar ret;
  ret = (end - ini) * 0.5 * (M_PI / moving_duration) *
        (M_PI / moving_duration) * cos(curr_time / moving_duration * M_PI);
  if (curr_time > moving_duration) {
    ret = 0.0;
  }
  return ret;
};

template <typename Scalar>
Scalar smoothing(Scalar ini, Scalar fin, Scalar rat) {
//   Scalar ret(0.);
  if (rat < 0) {
    return ini;
  } else if (rat > 1) {
    return fin;
  } else {
    return ini + (fin - ini) * rat;
  }
};

template <typename Scalar>
void getSinusoidTrajectory(Scalar initTime_, const VectorX<Scalar>& midPoint_,
                           const VectorX<Scalar>& amp_,
                           const VectorX<Scalar>& freq_, Scalar evalTime_,
                           VectorX<Scalar>& p_, VectorX<Scalar>& v_,
                           VectorX<Scalar>& a_) {
  int dim = midPoint_.size();
  p_ = VectorX<Scalar>::Zero(dim);
  v_ = VectorX<Scalar>::Zero(dim);
  a_ = VectorX<Scalar>::Zero(dim);
  for (int i = 0; i < dim; ++i) {
    p_[i] = amp_[i] * sin(2 * M_PI * freq_[i] * (evalTime_ - initTime_)) +
            midPoint_[i];
    v_[i] = amp_[i] * 2 * M_PI * freq_[i] *
            cos(2 * M_PI * freq_[i] * (evalTime_ - initTime_));
    a_[i] = -amp_[i] * 2 * M_PI * freq_[i] * 2 * M_PI * freq_[i] *
            sin(2 * M_PI * freq_[i] * (evalTime_ - initTime_));
  }
  Scalar smoothing_dur(1.0);
  if (evalTime_ < (initTime_ + smoothing_dur)) {
    for (int i = 0; i < dim; ++i) {
      v_[i] = 0. + v_[i] * (evalTime_ - initTime_) / smoothing_dur;
      a_[i] = 0. + a_[i] * (evalTime_ - initTime_) / smoothing_dur;
    }
  }
};

template <typename Scalar>
Scalar bind_half_pi(Scalar ang) {
    if (ang > M_PI / 2) {
        return ang - M_PI;
    }
    if (ang < -M_PI / 2) {
        return ang + M_PI;
    }
    return ang;
};

template <typename Scalar>
bool isInBoundingBox(const VectorX<Scalar>& val, const VectorX<Scalar>& lb,
                     const VectorX<Scalar>& ub) {
    int n = lb.size();
    bool ret(true);
    for (int i = 0; i < n; ++i) {
        if (lb[i] <= val[i] && val[i] <= ub[i]) {
        } else {
            rossy_utils::color_print(myColor::BoldMagneta, "Is not BoundingBox");
            std::cout << i << " th : lb = " << lb[i] << " val = " << val[i]
                      << " ub = " << ub[i] << std::endl;
            ret = false;
        }
    }
    return ret;
};

template <typename Scalar>
VectorX<Scalar> eulerIntegration(const VectorX<Scalar>& x,
                                 const VectorX<Scalar>& xdot, Scalar dt) {
    VectorX<Scalar> ret = x;
    ret += xdot * dt;
    return ret;
};

template <typename Scalar>
VectorX<Scalar> ScalarIntegration(const VectorX<Scalar>& q,
                                  const VectorX<Scalar>& alpha,
                                  const VectorX<Scalar>& alphad, Scalar dt) {
    VectorX<Scalar> ret = q;
    ret += alpha * dt + alphad * dt * dt * 0.5;
    return ret;
};

template <typename Scalar>
Scalar CropValue(Scalar value, Scalar min, Scalar max, std::string source) {
    assert(min < max);
    if (value > max) {
        printf("%s: %f is cropped to %f.\n", source.c_str(), value, max);
        value = max;
    }
    if (value < min) {
        printf("%s: %f is cropped to %f.\n", source.c_str(), value, min);
        value = min;
    }
    return value;
};

template <typename Scalar>
Scalar CropValue(Scalar value, Scalar min, Scalar max) {
    assert(min < max);
    if (value > max) {
        value = max;
    }
    if (value < min) {
        value = min;
    }
    return value;
};

template <typename Scalar>
VectorX<Scalar> CropVector(VectorX<Scalar> value, VectorX<Scalar> min,
                           VectorX<Scalar> max, std::string source) {
    assert(value.size() == min.size());
    assert(value.size() == max.size());
    int n_data = value.size();

    for (int i = 0; i < n_data; ++i) {
        if (value[i] > max[i]) {
            // printf("%s(%d): %f is cropped to %f\n", source.c_str(), i,
            // value[i], max[i]);
            value[i] = max[i];
        }
        if (value[i] < min[i]) {
            // printf("%s(%d): %f is cropped to %f\n", source.c_str(), i,
            // value[i], min[i]);
            value[i] = min[i];
        }
    }
    return value;
};

template <typename Scalar>
MatrixX<Scalar> CropMatrix(MatrixX<Scalar> value, MatrixX<Scalar> min,
                           MatrixX<Scalar> max, std::string source) {
    assert((value.cols() == min.cols()) && (value.cols() == max.cols()));
    assert((value.rows() == min.rows()) && (value.cols() == max.cols()));

    int n_row = value.rows();
    int n_cols = value.cols();

    for (int row_idx = 0; row_idx < n_row; ++row_idx) {
        for (int col_idx = 0; col_idx < n_cols; ++col_idx) {
            if (value(row_idx, col_idx) < min(row_idx, col_idx)) {
                // printf("%s(%d, %d): %f is cropped to %f\n", source.c_str(),
                // row_idx, col_idx, value(row_idx, col_idx), min(row_idx,
                // col_idx));
                value(row_idx, col_idx) = min(row_idx, col_idx);
            }
            if (value(row_idx, col_idx) > max(row_idx, col_idx)) {
                // printf("%s(%d, %d): %f is cropped to %f\n", source.c_str(),
                // row_idx, col_idx, value(row_idx, col_idx), max(row_idx,
                // col_idx));
                value(row_idx, col_idx) = max(row_idx, col_idx);
            }
        }
    }
    return value;
};

template <typename Scalar>
MatrixX<Scalar> GetRelativeMatrix(const MatrixX<Scalar> value,
                                  const MatrixX<Scalar> min,
                                  const MatrixX<Scalar> max) {
    assert((value.cols() == min.cols()) && (value.cols() == max.cols()));
    assert((value.rows() == min.rows()) && (value.cols() == max.cols()));

    MatrixX<Scalar> ret = value;
    for (int col_idx = 0; col_idx < value.cols(); ++col_idx) {
        for (int row_idx = 0; row_idx < value.rows(); ++row_idx) {
            Scalar width = max(row_idx, col_idx) - min(row_idx, col_idx);
            Scalar mid = (max(row_idx, col_idx) + min(row_idx, col_idx)) / 2.0;
            ret(row_idx, col_idx) =
                2.0 * (value(row_idx, col_idx) - mid) / width;
        }
    }
    return ret;
};

template <typename Scalar>
VectorX<Scalar> GetRelativeVector(const VectorX<Scalar> value,
                                  const VectorX<Scalar> min,
                                  const VectorX<Scalar> max) {
    assert((value.size() == min.size()) && (value.size() == max.size()));
    VectorX<Scalar> ret = value;

    for (int idx = 0; idx < value.size(); ++idx) {
        Scalar width = max(idx) - min(idx);
        Scalar mid = (max(idx) + min(idx)) / 2.0;
        ret(idx) = 2.0 * (value(idx) - mid) / width;
    }
    return ret;
};

}  // namespace rossy_utils
