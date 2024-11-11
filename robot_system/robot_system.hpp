#pragma once

#include <stdio.h>
#include <Eigen/Dense>
#include <string>

#include <pinocchio/algorithm/aba.hpp>
#include <pinocchio/algorithm/center-of-mass.hpp>
#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/fwd.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include "rossy_utils/math/typedefs.h"

template <typename Scalar>
class RobotSystem {
   protected:
    // pinocchio system
    pinocchio::ModelTpl<Scalar> model_;
    pinocchio::DataTpl<Scalar> data_;

    std::string urdf_file_;
    bool b_print_info_;

    // generalized coordinate configuration
    VectorX<Scalar> q_;
    VectorX<Scalar> qdot_;
    VectorX<Scalar> qddot_;

    // robot info
    Scalar total_mass_;
    int n_q_;
    int n_qdot_;
    int n_dof_;
    int n_link_;
    std::vector<int> idx_adof_;
    std::map<std::string, size_t> link_idx_map_;
    std::map<std::string, size_t> joint_idx_map_;
    std::map<size_t, std::string> link_idx_map_inv_;
    std::map<size_t, std::string> joint_idx_map_inv_; 

   public:
    RobotSystem(const RobotSystem<Scalar>& robotsys); // copier
    RobotSystem(const std::string& file);
    virtual ~RobotSystem(void);

    void printRobotInfo();
    std::string_view getRobotName(){return model_.name;}
    std::string_view getUrdfFile(){return urdf_file_;}

    // update fixed base system
    void updateSystem(const VectorX<Scalar> &joint_pos,
                    const VectorX<Scalar> &joint_vel);

    VectorX<Scalar> getQ() { return q_; };
    VectorX<Scalar> getQdot() { return qdot_; };
    VectorX<Scalar> getQddot() { return qddot_; };

    VectorX<Scalar> GetTorqueLowerLimits() { return -model_.effortLimit;}
    VectorX<Scalar> GetTorqueUpperLimits() { return model_.effortLimit;}
    VectorX<Scalar> GetPositionLowerLimits() { return model_.lowerPositionLimit;}
    VectorX<Scalar> GetPositionUpperLimits() { return model_.upperPositionLimit;} 
    VectorX<Scalar> GetVelocityLowerLimits() { return -model_.velocityLimit;}
    VectorX<Scalar> GetVelocityUpperLimits() { return model_.velocityLimit;}   
    
    Scalar getRobotMass() { return total_mass_; }
    int getNumDofs() { return n_qdot_; };
    int getNumBodyNodes() { return n_link_; };

    int getJointIdx(const std::string& joint_name);
    int getLinkIdx(const std::string& frame_name);
    std::string_view getLinkName(const int& frame_idx);
    std::string_view getJointName(const int& joint_idx); 

    MatrixX<Scalar> getMassMatrix();
    MatrixX<Scalar> getInvMassMatrix();
    VectorX<Scalar> getCoriolisGravity();
    VectorX<Scalar> getGravity();
    VectorX<Scalar> getCoriolis();
    MatrixX<Scalar> getCoriolisMatrix();

    Isometry3<Scalar> getBodyNodeIsometry(const std::string& name_);
    Eigen::Matrix<Scalar, 6, 1> getBodyNodeSpatialVelocity(const std::string& name_);    
    MatrixX<Scalar> getBodyNodeJacobian(const std::string& name_);
    VectorX<Scalar> getBodyNodeJacobianDotQDot(const std::string& name_);
    Eigen::Matrix<Scalar, 6, 1> getBodyNodeBodyVelocity(const std::string& name_);
    MatrixX<Scalar> getBodyNodeBodyJacobian(const std::string& name_);
    VectorX<Scalar> getBodyNodeBodyJacobianDotQDot(const std::string& name_);

    Isometry3<Scalar> getBodyNodeIsometry(const int& _bn_idx);
    Eigen::Matrix<Scalar, 6, 1> getBodyNodeSpatialVelocity(const int& _bn_idx);    
    MatrixX<Scalar> getBodyNodeJacobian(const int& _bn_idx);
    VectorX<Scalar> getBodyNodeJacobianDotQDot(const int& _bn_idx);
    Eigen::Matrix<Scalar, 6, 1> getBodyNodeBodyVelocity(const int& _bn_idx);
    MatrixX<Scalar> getBodyNodeBodyJacobian(const int& _bn_idx);
    VectorX<Scalar> getBodyNodeBodyJacobianDotQDot(const int& _bn_idx);


  private:
    void _initializeRobotInfo();
    void _updateSystemData();
};
