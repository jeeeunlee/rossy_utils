#include "rossy_utils/robot_system/robot_system.hpp"
#include "rossy_utils/io/io_utilities.hpp"
#include <chrono>

template <typename Scalar>
RobotSystem<Scalar>::RobotSystem(const RobotSystem<Scalar>& robotsys) 
{
    rossy_utils::pretty_constructor(2, "Robot Model(temp)");
    urdf_file_ = robotsys.urdf_file_;       
    _initializeRobotInfo();
    // printRobotInfo();
}

template <typename Scalar>
RobotSystem<Scalar>::RobotSystem(const std::string& file): urdf_file_(file) {
    rossy_utils::pretty_constructor(1, "Robot Model");
    rossy_utils::color_print(myColor::BoldCyan, "|", false);
    rossy_utils::color_print(myColor::Green, "        path="+urdf_file_);
    _initializeRobotInfo();
    // printRobotInfo();
}

template <typename Scalar>
RobotSystem<Scalar>::~RobotSystem() {}

template <typename Scalar>
void RobotSystem<Scalar>::_initializeRobotInfo() {

    // set pinocchio model & data    
    // pinocchio::urdf::buildModel<Scalar>(urdf_file_, model_);
    pinocchio::Model pinocchio_model_double;
    pinocchio::urdf::buildModel(urdf_file_, pinocchio_model_double);
    model_ = pinocchio_model_double.cast<Scalar>();

    data_ = pinocchio::DataTpl<Scalar>(model_);
    n_q_ = model_.nq;
    n_qdot_ = model_.nv;
    n_dof_ = n_qdot_;

    q_ = VectorX<Scalar>::Zero(n_q_);
    qdot_ = VectorX<Scalar>::Zero(n_qdot_);
    qddot_ = VectorX<Scalar>::Zero(n_qdot_);

    // UPDATE link_idx_map_, joint_idx_map_
    for (pinocchio::FrameIndex i(0); // FrameIndex : size_t
        i < static_cast<pinocchio::FrameIndex>(model_.nframes); ++i) {
        if(model_.frames[i].type == pinocchio::FrameType::BODY){
            std::string frame_name = model_.frames[i].name;
            link_idx_map_[frame_name] = model_.getBodyId(frame_name);
            link_idx_map_inv_[model_.getBodyId(frame_name)] = frame_name;             
        }
    }
    for (pinocchio::JointIndex i(0);
        i < static_cast<pinocchio::JointIndex>(model_.njoints); ++i) {
        total_mass_ += model_.inertias[i].mass();
        std::string joint_name = model_.names[i];
        if (joint_name != "universe" && joint_name != "root_joint"){
            joint_idx_map_[joint_name] = model_.getJointId(joint_name); // i - 2; joint map excluding fixed joint
            joint_idx_map_inv_[model_.getJointId(joint_name)] = joint_name; // joint map excluding fixed joint
        }           
    }
    n_link_ = link_idx_map_.size();
    assert(n_dof_ == joint_idx_map_.size());
}

template <typename Scalar>
MatrixX<Scalar> RobotSystem<Scalar>::getMassMatrix() {
    // data_.M.triangularView<Eigen::StrictlyLower>() =
    //     data_.M.transpose().triangularView<Eigen::StrictlyLower>();
    data_.M = 0.5*(data_.M + data_.M.transpose());
    return data_.M;
}

template <typename Scalar>
MatrixX<Scalar> RobotSystem<Scalar>::getInvMassMatrix() {    
    data_.Minv = 0.5*(data_.Minv + data_.Minv.transpose());
    return data_.Minv;
}

template <typename Scalar>
VectorX<Scalar> RobotSystem<Scalar>::getCoriolisGravity() {
    return pinocchio::nonLinearEffects(model_, data_, q_, qdot_);
}

template <typename Scalar>
MatrixX<Scalar> RobotSystem<Scalar>::getCoriolisMatrix() {
    return computeCoriolisMatrix(model_, data_, q_, qdot_);
}

template <typename Scalar>
VectorX<Scalar> RobotSystem<Scalar>::getCoriolis() {
    return pinocchio::nonLinearEffects(model_, data_, q_, qdot_) -
         pinocchio::computeGeneralizedGravity(model_, data_, q_);
}

template <typename Scalar>
VectorX<Scalar> RobotSystem<Scalar>::getGravity() {
    return pinocchio::computeGeneralizedGravity(model_, data_, q_);
}

template <typename Scalar>
Isometry3<Scalar> RobotSystem<Scalar>::getBodyNodeIsometry(const std::string& name) {
    return this->getBodyNodeIsometry(link_idx_map_[name]); }

template <typename Scalar>
Isometry3<Scalar> RobotSystem<Scalar>::getBodyNodeIsometry(const int& link_idx) {
    Isometry3<Scalar> ret;
    const pinocchio::SE3Tpl<Scalar> trans =
        pinocchio::updateFramePlacement(model_, data_, link_idx);
    ret.template linear() = trans.rotation();
    ret.template translation() = trans.translation();
    return ret;
}

template <typename Scalar>
Eigen::Matrix<Scalar, 6, 1> RobotSystem<Scalar>::getBodyNodeSpatialVelocity(const std::string& name) {
    return this->getBodyNodeSpatialVelocity(link_idx_map_[name]); }

template <typename Scalar>
Eigen::Matrix<Scalar, 6, 1> RobotSystem<Scalar>::getBodyNodeSpatialVelocity(const int& link_idx) {
  Eigen::Matrix<Scalar, 6, 1> ret = Eigen::Matrix<Scalar, 6, 1>::Zero();
  pinocchio::MotionTpl<Scalar> fv = pinocchio::getFrameVelocity(
                                        model_, 
                                        data_, 
                                        link_idx, 
                                        pinocchio::LOCAL_WORLD_ALIGNED);
  ret.template head<3>() = fv.angular();
  ret.template tail<3>() = fv.linear();
  return ret;
}

template <typename Scalar>
MatrixX<Scalar> RobotSystem<Scalar>::getBodyNodeJacobian(const std::string& name) {
    return this->getBodyNodeJacobian(link_idx_map_[name]); }


template <typename Scalar>
MatrixX<Scalar> RobotSystem<Scalar>::getBodyNodeJacobian(const int& link_idx) {
    // Analytic Jacobian
  Eigen::Matrix<Scalar, 6, Eigen::Dynamic> jac =
      Eigen::Matrix<Scalar, 6, Eigen::Dynamic>::Zero(6, n_qdot_);
  pinocchio::getFrameJacobian(model_, data_, link_idx,
                              pinocchio::LOCAL_WORLD_ALIGNED, jac);
  Eigen::Matrix<Scalar, 6, Eigen::Dynamic> ret =
      Eigen::Matrix<Scalar, 6, Eigen::Dynamic>::Zero(6, n_qdot_);
  ret.topRows(3) = jac.bottomRows(3);
  ret.bottomRows(3) = jac.topRows(3);
  return ret;
}

template <typename Scalar>
VectorX<Scalar> RobotSystem<Scalar>::getBodyNodeJacobianDotQDot(const std::string& name) {
    return this->getBodyNodeJacobianDotQDot(link_idx_map_[name]); }

template <typename Scalar>
VectorX<Scalar> RobotSystem<Scalar>::getBodyNodeJacobianDotQDot(const int& link_idx) {
    // check 
    // pinocchio::Motion fa = pinocchio::getFrameAcceleration(
    //     model_, data_, link_idx, pinocchio::LOCAL_WORLD_ALIGNED);

    // pinocchio::forwardKinematics(model_, data_, q_, qdot_, 0 * qdot_);
    pinocchio::MotionTpl<Scalar> fa = pinocchio::getFrameClassicalAcceleration(
        model_, data_, link_idx, pinocchio::LOCAL_WORLD_ALIGNED);

    Eigen::Matrix<Scalar, 6, 1> ret = Eigen::Matrix<Scalar, 6, 1>::Zero();
    ret.template segment(0, 3) = fa.angular();
    ret.template segment(3, 3) = fa.linear();

    return ret;
}

template <typename Scalar>
Eigen::Matrix<Scalar, 6, 1> RobotSystem<Scalar>::getBodyNodeBodyVelocity(const std::string& name) {
    return this->getBodyNodeBodyVelocity(link_idx_map_[name]);
}

template <typename Scalar>
Eigen::Matrix<Scalar, 6, 1> RobotSystem<Scalar>::getBodyNodeBodyVelocity(const int& link_idx) {
  Eigen::Matrix<Scalar, 6, 1> ret = Eigen::Matrix<Scalar, 6, 1>::Zero();
  pinocchio::MotionTpl<Scalar> fv = pinocchio::getFrameVelocity(
                                    model_, data_, link_idx, pinocchio::LOCAL);
  ret.template head<3>() = fv.angular();
  ret.template tail<3>() = fv.linear();
  return ret;
}

template <typename Scalar>
MatrixX<Scalar> RobotSystem<Scalar>::getBodyNodeBodyJacobian(const std::string& name) {
    return this->getBodyNodeBodyJacobian(link_idx_map_[name]);
}

template <typename Scalar>
MatrixX<Scalar> RobotSystem<Scalar>::getBodyNodeBodyJacobian(const int& link_idx) {
  Eigen::Matrix<Scalar, 6, Eigen::Dynamic> jac =
      Eigen::Matrix<Scalar, 6, Eigen::Dynamic>::Zero(6, n_qdot_);
  pinocchio::getFrameJacobian(model_, data_, link_idx, pinocchio::LOCAL, jac);
  Eigen::Matrix<Scalar, 6, Eigen::Dynamic> ret =
      Eigen::Matrix<Scalar, 6, Eigen::Dynamic>::Zero(6, n_qdot_);
  ret.topRows(3) = jac.bottomRows(3);
  ret.bottomRows(3) = jac.topRows(3);
  return ret;
}

template <typename Scalar>
VectorX<Scalar> RobotSystem<Scalar>::getBodyNodeBodyJacobianDotQDot(const std::string& name) {
        return this->getBodyNodeBodyJacobianDotQDot(link_idx_map_[name]);
}

template <typename Scalar>
VectorX<Scalar> RobotSystem<Scalar>::getBodyNodeBodyJacobianDotQDot(const int& link_idx) {
    // pinocchio::forwardKinematics(model_, data_, q_, qdot_, 0 * qdot_);
    pinocchio::MotionTpl<Scalar> fa = pinocchio::getFrameClassicalAcceleration(
        model_, data_, link_idx, pinocchio::LOCAL);

    Eigen::Matrix<Scalar, 6, 1> ret = Eigen::Matrix<Scalar, 6, 1>::Zero();
    ret.template segment(0, 3) = fa.angular();
    ret.template segment(3, 3) = fa.linear();

    return ret;
}

template <typename Scalar>
int RobotSystem<Scalar>::getLinkIdx(const std::string& frame_name) {    
    return link_idx_map_[frame_name]; }

template <typename Scalar>
int RobotSystem<Scalar>::getJointIdx(const std::string& jointName) {    
    return joint_idx_map_[jointName]; }

template <typename Scalar>
std::string_view RobotSystem<Scalar>::getLinkName(const int& frame_idx) {    
    return link_idx_map_inv_[frame_idx]; }

template <typename Scalar>
std::string_view RobotSystem<Scalar>::getJointName(const int& joint_idx) {    
    return joint_idx_map_inv_[joint_idx]; }

template <typename Scalar>
void RobotSystem<Scalar>::updateSystem(const VectorX<Scalar> &joint_pos,
                                const VectorX<Scalar> &joint_vel) {
    // ASSUME FIXED BASE
    q_ = joint_pos;
    qdot_ = joint_vel;

    // update data
    _updateSystemData();
}

template <typename Scalar>
void RobotSystem<Scalar>::_updateSystemData(){
    pinocchio::crba(model_, data_, q_);
    pinocchio::forwardKinematics(model_, data_, q_, qdot_, 0 * qdot_);
    pinocchio::computeJointJacobians(model_, data_, q_);
}

template <typename Scalar>
void RobotSystem<Scalar>::printRobotInfo() {
    std::cout << " ==== Robot ====" << std::endl;
    std::cout << model_.name << std::endl;
    std::cout << " ==== Body Node ====" << std::endl;
    for (auto &[idx, name] : link_idx_map_inv_) {
        std::cout << "constexpr int " << name  << " = "
                  << std::to_string(idx) << ";" << std::endl;
    }
    std::cout << " ==== DoF ====" << std::endl;
    for (auto &[idx, name] : joint_idx_map_inv_) {
        std::cout << "constexpr int " << name << " = " 
                  << std::to_string(idx) << ";" << std::endl;
    }
    std::cout << " ==== Num ====" << std::endl;
    std::cout << "constexpr int n_bodynode = "
              << std::to_string(n_link_) << ";"
              << std::endl;
    std::cout << "constexpr int n_dof = " << std::to_string(n_qdot_) << ";"
              << std::endl;

    std::cout << " ==== Just info ====" << std::endl;
    std::cout << n_link_ <<", "  // 86
            << n_q_ << ", "  // 25 7 + 12 + 6
            << n_qdot_ <<", "  // 24  6 + 12 + 6
            << n_dof_ << std::endl; // 18 = 12 + 6

}

template class RobotSystem<double>;
template class RobotSystem<float>;

