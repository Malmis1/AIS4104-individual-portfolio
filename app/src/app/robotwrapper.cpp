#include "app/robotwrapper.h"

#include <utility/math.h>
#include <utility/vectors.h>

using namespace AIS4104;

RobotWrapper::RobotWrapper(std::shared_ptr<threepp::Robot> robot, std::shared_ptr<Simulation::KinematicsSolver> solver)
    : m_tool_transform(Eigen::Matrix4d::Identity())
    , m_robot(std::move(robot))
    , m_solver(std::move(solver))
{
}

threepp::Robot& RobotWrapper::threepp_robot()
{
    return *m_robot;
}

std::shared_ptr<threepp::Robot> RobotWrapper::threepp_robot_ptr()
{
    return m_robot;
}

const threepp::Robot& RobotWrapper::threepp_robot() const
{
    return *m_robot;
}

uint8_t RobotWrapper::joint_count() const
{
    return m_solver->joint_count();
}

//TASK: Implement the function to calculate the joint positions for the desired tool pose
// a) Use m_tool_transform to calculate the flange pose required by m_solver.ik_solve()
// b) Use the m_solver.ik_solve() overload with the solution selector lambda to choose the most desirable IK solution.
Eigen::VectorXd RobotWrapper::ik_solve_pose(const Eigen::Matrix4d& eef_pose, const Eigen::VectorXd& j0) const
{
    return joint_positions();
}

Eigen::VectorXd RobotWrapper::ik_solve_flange_pose(const Eigen::Matrix4d& flange_pose, const Eigen::VectorXd& j0) const
{
    return m_solver->ik_solve(flange_pose, j0);
}

Eigen::Matrix4d RobotWrapper::tool_transform() const
{
    return m_tool_transform;
}

void RobotWrapper::set_tool_transform(Eigen::Matrix4d transform)
{
    m_tool_transform = std::move(transform);
}

//TASK: Calculate the pose of the end effector using forward kinematics;
// Relevant variables are m_solver and m_tool_transform.
Eigen::Matrix4d RobotWrapper::current_pose() const
{
    return Eigen::Matrix4d::Identity();
}

//TASK: Calculate the position of the end effector using forward kinematics.
// Relevant variables are m_solver and m_tool_transform (or possibly another function of RobotWrapper?).
Eigen::Vector3d RobotWrapper::current_position() const
{
    return Eigen::Vector3d::Zero();
}

//TASK: Calculate the orientation of the end effector using forward kinematics and m_solver (or rely on another function of RobotWrapper?).
Eigen::Vector3d RobotWrapper::current_orientation_zyx() const
{
    return Eigen::Vector3d::Zero();
}

Eigen::Matrix4d RobotWrapper::current_flange_pose() const
{
    return m_solver->fk_solve(joint_positions());
}

Eigen::Vector3d RobotWrapper::current_flange_position() const
{
    Eigen::Matrix4d tf = m_solver->fk_solve(joint_positions());
    return tf.block<3, 1>(0, 3);
}

Eigen::Vector3d RobotWrapper::current_flange_orientation_zyx() const
{
    Eigen::Matrix4d tf = m_solver->fk_solve(joint_positions());
    return utility::euler_zyx_from_rotation_matrix(tf.block<3, 3>(0, 0));
}

const Simulation::JointLimits& RobotWrapper::joint_limits() const
{
    return m_solver->joint_limits();
}

Eigen::VectorXd RobotWrapper::joint_positions() const
{
    return utility::to_eigen_vectord(m_robot->jointValues());
}

void RobotWrapper::set_joint_positions(const Eigen::VectorXd& joint_positions)
{
    m_robot->setJointValues(utility::to_std_vectorf(joint_positions));
}
