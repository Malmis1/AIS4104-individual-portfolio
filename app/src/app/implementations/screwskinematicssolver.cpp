#include "app/implementations/screwskinematicssolver.h"

#include <utility/math.h>

using namespace AIS4104;

ScrewsKinematicsSolver::ScrewsKinematicsSolver(Eigen::Matrix4d m, std::vector<Eigen::VectorXd> screws, Simulation::JointLimits limits)
    : ScrewsKinematicsSolver(std::move(m), std::move(screws), 4.e-3, 4.e-3, std::move(limits))
{
}

ScrewsKinematicsSolver::ScrewsKinematicsSolver(Eigen::Matrix4d m, std::vector<Eigen::VectorXd> space_screws, double v_e, double w_e, Simulation::JointLimits limits)
    : KinematicsSolver(std::move(limits))
    , m_ve(v_e)
    , m_we(w_e)
    , m_m(std::move(m))
    , m_screws(std::move(space_screws))
{
}

void ScrewsKinematicsSolver::set_epsilons(double v_e, double w_e)
{
    m_ve = v_e;
    m_we = w_e;
}

uint32_t ScrewsKinematicsSolver::joint_count() const
{
    return m_screws.size();
}

Eigen::Matrix4d ScrewsKinematicsSolver::fk_solve(const Eigen::VectorXd& joint_positions)
{
    // Equation (4.14) on page 140, MR 3rd print 2019
    Eigen::Matrix4d product = Eigen::Matrix4d::Identity();

    for (size_t i = 0; i < joint_positions.size(); i++)
    {
        product *= utility::matrix_exponential(m_screws[i], joint_positions[i]);
    }

    return product * m_m;
}

Eigen::VectorXd ScrewsKinematicsSolver::ik_solve(const Eigen::Matrix4d& t_sd, const Eigen::VectorXd& j0)
{
    return ik_solve(t_sd, j0, [&](const std::vector<Eigen::VectorXd>&) { return 0u; });
}

Eigen::VectorXd ScrewsKinematicsSolver::ik_solve(const Eigen::Matrix4d& t_sd, const Eigen::VectorXd& j0, const std::function<uint32_t(const std::vector<Eigen::VectorXd>&)>& solution_selector)
{
    // Inverse kinematics algorithm at bottom of page 228 and continuing on page 229, MR 3rd print 2019
    size_t iter = 0;
    size_t max_iter = 1000;

    Eigen::VectorXd joint_positions(6);
    Eigen::VectorXd previous_joint_positions(6);
    joint_positions = j0;
    previous_joint_positions = j0;

    Eigen::VectorXd result = joint_positions;
    bool crit = true;

    while ((iter < max_iter) && crit) {
        Eigen::Matrix4d t_sb_inv = fk_solve(previous_joint_positions).inverse();
        // Need to use static_case to specify which matrix_logarithm to use
        std::pair<Eigen::VectorXd, double> vb = utility::matrix_logarithm(static_cast<const Eigen::Matrix4d&>(t_sb_inv * t_sd));
        crit = (vb.first.block<3, 1>(0, 0).norm() > m_we) || (vb.first.block<3, 1>(3, 0).norm() > m_ve);

        Eigen::MatrixXd j_pinv(6, 6);
        // Built-in Eigen function for pseudoinverse
        j_pinv = body_jacobian(previous_joint_positions).completeOrthogonalDecomposition().pseudoInverse();

        Eigen::VectorXd temp_joint_positions = joint_positions;
        joint_positions = previous_joint_positions + j_pinv * vb.first;
        previous_joint_positions = temp_joint_positions;

        result = joint_positions;
        iter++;
    }

    return result;
}

std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ScrewsKinematicsSolver::space_chain()
{
    return { m_m, m_screws };
}

std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ScrewsKinematicsSolver::body_chain()
{
    auto [m, screws] = space_chain();

    // Below equation (4.16) on page 147, MR 3rd print 2019
    int count = joint_count();
    Eigen::MatrixXd adj_m = Eigen::MatrixXd::Identity(count, count);
    adj_m = utility::adjoint_matrix(m.inverse());

    std::vector<Eigen::VectorXd> b_screws;
    for (Eigen::VectorXd s : screws) {
        Eigen::VectorXd b = adj_m * s;
        b_screws.push_back(b);
    }

    return { m, b_screws };
}

Eigen::MatrixXd ScrewsKinematicsSolver::space_jacobian(const Eigen::VectorXd& current_joint_positions)
{
    auto [m, screws] = space_chain();

    // Equation (5.11) on page 178, MR 3rd print 2019
    int count = current_joint_positions.size();
    Eigen::Matrix4d t = Eigen::Matrix4d::Identity();
    Eigen::MatrixXd j = Eigen::MatrixXd::Identity(count, count);
    j.col(0) = screws[0];

    for (size_t i = 1; i < count; i++) {
        j.col(i) = utility::adjoint_matrix(t) * screws[i];
        t = t * utility::matrix_exponential(
            screws[i - 1].block<3, 1>(0, 0),
            screws[i - 1].block<3, 1>(3, 0),
            current_joint_positions[i - 1]);
    }

    return j;
}

Eigen::MatrixXd ScrewsKinematicsSolver::body_jacobian(const Eigen::VectorXd& current_joint_positions)
{
    auto [m, screws] = body_chain();

    // Equation (5.18) on page 183, MR 3rd print 2019
    int count = current_joint_positions.size();
    Eigen::Matrix4d t = Eigen::Matrix4d::Identity();
    Eigen::MatrixXd j = Eigen::MatrixXd::Identity(count, count);
    j.col(count - 1) = screws[count - 1];

    for (int i = count - 2; i >= 0; i--) {
        t = t * utility::matrix_exponential(
            -screws[i + 1].block<3, 1>(0, 0),
            -screws[i + 1].block<3, 1>(3, 0),
            current_joint_positions[i + 1]);
        j.col(i) = utility::adjoint_matrix(t) * screws[i];
    }

    return j;
}
