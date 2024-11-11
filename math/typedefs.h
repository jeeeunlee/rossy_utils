#include <Eigen/Dense>

template <typename Scalar>
using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Scalar>
using VectorX = Eigen::Vector<Scalar, Eigen::Dynamic>;

template <typename Scalar>
using Vector3 = Eigen::Vector<Scalar, 3>;

template <typename Scalar>
using Quaternion = Eigen::Quaternion<Scalar>;

template <typename Scalar>
using AngleAxis = Eigen::AngleAxis<Scalar>;

template <typename Scalar>
using Isometry3 = Eigen::Transform<Scalar,3,Eigen::Isometry>;