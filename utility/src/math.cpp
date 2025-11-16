#include "utility/math.h"
#include "utility/vectors.h"

#include <numbers>

namespace AIS4104::utility {

    bool double_equals(double a, double b, double epsilon = 0.000001)
    {
        // Taken from feedback
        return std::abs(a - b) < epsilon;
    }

    double cot(double x)
    {
        return 1 / (std::sin(x) / std::cos(x));
    }

    Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d& r)
    {
        // Equations on page 579 Section B.1.1, MR 3rd print 2019
        double alpha = std::atan2(r(1, 0), r(0, 0));
        double gamma = std::atan2(r(2, 1), r(2, 2));
        double beta;

        if (double_equals(r(2, 0), -1.0)) {
            // Equal to pi/2
            beta = std::acos(0.0);
        }
        else if (double_equals(r(2, 0), 1.0)) {
            // Equal to -pi/2
            beta = -std::acos(0.0);
        }
        else {
            beta = std::atan2(-r(2, 0), std::sqrt(std::pow(r(0, 0), 2) + std::pow(r(1, 0), 2)));
        }

        return Eigen::Vector3d{ alpha, beta, gamma } * rad_to_deg;
    }

    Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d& v)
    {
        // Equation (3.30) page 75, MR 3rd print 2019
        return Eigen::Matrix3d{ {0 , -v[2], v[1]},
                                {v[2], 0, -v[0]},
                                {-v[1], v[0], 0} };
    }

    Eigen::Vector3d from_skew_symmetric(const Eigen::Matrix3d& m)
    {
        // Based on values used in skew_symmetric function
        Eigen::Vector3d vec = { m(2,1), m(0,2), m(1,0) };

        return vec;
    }

    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix3d& r, const Eigen::Vector3d& p)
    {
        // Definition 3.20 on page 98, MR 3rd print 2019
        Eigen::MatrixXd adj_matrix(6, 6);
        adj_matrix << r, Eigen::Matrix3d::Zero(),
            skew_symmetric(p)* r, r;

        return adj_matrix;
    }

    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d& tf)
    {
        Eigen::Matrix3d r = tf.block<3, 3>(0, 0);
        Eigen::Vector3d p = tf.block<3, 1>(0, 3);

        return adjoint_matrix(r, p);
    }

    Eigen::VectorXd adjoint_map(const Eigen::VectorXd& twist, const Eigen::Matrix4d& tf)
    {
        // Definition 3.20 on page 98, MR 3rd print 2019
        Eigen::VectorXd adj_v(6, 1);
        adj_v << adjoint_matrix(tf) * twist;

        return adj_v;
    }

    Eigen::VectorXd twist(const Eigen::Vector3d& w, const Eigen::Vector3d& v)
    {
        // Equation (3.70) on page 96, MR 3rd print 2019
        Eigen::VectorXd vb(6, 1);
        vb << w, v;

        return vb;
    }

    Eigen::VectorXd twist(const Eigen::Vector3d& q, const Eigen::Vector3d& s, double h, double angular_velocity)
    {
        // Equation on page 101, MR 3rd print 2019
        Eigen::VectorXd vb(6, 1);
        Eigen::Vector3d w = s * angular_velocity;
        Eigen::Vector3d v = skew_symmetric(-s * angular_velocity) * q + h * s * angular_velocity;
        vb << w, v;

        return vb;
    }

    Eigen::Matrix4d twist_matrix(const Eigen::Vector3d& w, const Eigen::Vector3d& v)
    {
        // Equation (3.71) on page 96, MR 3rd print 2019
        Eigen::Matrix4d v_matrix;
        v_matrix << skew_symmetric(w), v,
            0, 0;

        return v_matrix;
    }

    Eigen::Matrix4d twist_matrix(const Eigen::VectorXd& twist)
    {
        const Eigen::Vector3d& w = twist.head<3>();
        const Eigen::Vector3d& v = twist.tail<3>();

        return twist_matrix(w, v);
    }

    Eigen::VectorXd screw_axis(const Eigen::Vector3d& w, const Eigen::Vector3d& v)
    {
        // Definition 3.24 on page 102, MR 3rd print 2019
        Eigen::VectorXd s(6, 1);
        s << w, v;

        return s;
    }

    Eigen::VectorXd screw_axis(const Eigen::Vector3d& q, const Eigen::Vector3d& s, double h)
    {
        // Equation on page 101, MR 3rd print 2019
        Eigen::VectorXd v(6, 1);
        v << s, -skew_symmetric(s) * q + h * s;

        return v;
    }

    Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d& w, double theta)
    {
        double radians = deg_to_rad * theta;

        // Equation (3.51) on page 82, MR 3rd print 2019
        return Eigen::Matrix3d::Identity() + std::sin(radians) * skew_symmetric(w)
            + (1.0 - std::cos(radians)) * (skew_symmetric(w) * skew_symmetric(w));
    }

    Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d& w, const Eigen::Vector3d& v, double theta)
    {
        Eigen::Matrix3d skew_w = skew_symmetric(w);
        Eigen::Matrix4d m_e;
        double radians = deg_to_rad * theta;

        // Proposition 3.25 on page 103, MR 3rd print 2019
        if (double_equals(w.norm(), 1.0)) {
            Eigen::Matrix3d r = matrix_exponential(w, theta);
            Eigen::Vector3d p = (Eigen::Matrix3d::Identity() * radians + (1 - std::cos(radians)) * skew_w + (radians - std::sin(radians)) * skew_w * skew_w) * v;
            m_e = transformation_matrix(r, p);
        }
        else if (double_equals(w.norm(), 0.0) && double_equals(v.norm(), 1.0)) {
            Eigen::Matrix3d r = Eigen::Matrix3d::Identity();
            Eigen::Vector3d p = v * theta;
            m_e = transformation_matrix(r, p);
        }

        return m_e;
    }

    Eigen::Matrix4d matrix_exponential(const Eigen::VectorXd& screw, double theta)
    {
        const Eigen::Vector3d& w = screw.head<3>();
        const Eigen::Vector3d& v = screw.tail<3>();

        return matrix_exponential(w, v, theta);
    }

    std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d& r)
    {
        Eigen::Vector3d w;
        double theta;
        double epsilon = 0.000001;

        // Algorithm from page 85 (Equations (3.58)-(3.61) on pages 85-86), MR 3rd print 2019
        if (r.isApprox(Eigen::Matrix3d::Identity(), epsilon)) {
            theta = 0.0;
        }
        else {
            // Built in function in Eigen for Equation (3.54) on page 84, MR 3rd print 2019
            double trr = r.trace();
            if (double_equals(trr, -1.0)) {
                theta = std::numbers::pi;
                // Equation (3.58) on page 85, MR 3rd print 2019
                w = (1.0 / std::sqrt(2.0 * (1.0 + r(2, 2))))
                    * Eigen::Vector3d{ r(0, 2), r(1, 2), 1.0 + r(2, 2) };
            }
            else {
                theta = std::acos((1.0 / 2.0) * (trr - 1.0));
                // Equations directly above Equation (3.53) on page 84, MR 3rd print 2019
                double w_1 = ((1.0 / (2.0 * std::sin(theta))) * (r(2, 1) - r(1, 2)));
                double w_2 = ((1.0 / (2.0 * std::sin(theta))) * (r(0, 2) - r(2, 0)));
                double w_3 = ((1.0 / (2.0 * std::sin(theta))) * (r(1, 0) - r(0, 1)));
                w = Eigen::Vector3d{ w_1, w_2, w_3 };
            }
        }

        return std::pair<Eigen::Vector3d, double>{ w, theta };
    }

    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix3d& r, const Eigen::Vector3d& p)
    {
        return matrix_logarithm(transformation_matrix(r, p));
    }

    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d& tf)
    {
        Eigen::Matrix3d r = tf.block<3, 3>(0, 0);
        Eigen::Vector3d p = tf.block<3, 1>(0, 3);
        Eigen::Vector3d w;
        Eigen::Vector3d v;
        double h = 0.0;
        double theta;
        double epsilon = 0.000001;

        // Algorithm in section 3.3.3.2 on page 104, MR 3rd print 2019
        if (r.isApprox(Eigen::Matrix3d::Identity(), epsilon)) {
            w = Eigen::Vector3d::Zero();
            v = p / p.norm();
            theta = p.norm();
        }
        else {
            std::pair<Eigen::Vector3d, double> m_l = matrix_logarithm(r);
            w = m_l.first;
            theta = m_l.second;
            Eigen::Matrix3d skew_w = skew_symmetric(w);
            v = ((1 / theta) * Eigen::Matrix3d::Identity() - (1.0 / 2.0) * skew_w
                + ((1.0 / theta) - (1.0 / 2.0) * cot(theta / 2.0)) * skew_w * skew_w) * p;
        }

        return std::pair<Eigen::VectorXd, double>{ screw_axis(w, v, h), theta };
    }

    Eigen::Matrix3d rotate_x(double radians)
    {
        Eigen::Matrix3d matrix;

        // Equations on page 72, MR 3rd print 2019
        matrix <<
            1, 0, 0,
            0, std::cos(radians), -std::sin(radians),
            0, std::sin(radians), std::cos(radians);

        return matrix;
    }

    Eigen::Matrix3d rotate_y(double radians)
    {
        Eigen::Matrix3d matrix;

        // Equations on page 72, MR 3rd print 2019
        matrix <<
            std::cos(radians), 0, std::sin(radians),
            0, 1, 0,
            -std::sin(radians), 0, std::cos(radians);

        return matrix;
    }

    Eigen::Matrix3d rotate_z(double radians)
    {
        Eigen::Matrix3d matrix;

        // Equations on page 72, MR 3rd print 2019
        matrix <<
            std::cos(radians), -std::sin(radians), 0,
            std::sin(radians), std::cos(radians), 0,
            0, 0, 1;

        return matrix;
    }

    Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d& x, const Eigen::Vector3d& y, const Eigen::Vector3d& z)
    {
        Eigen::Matrix3d matrix;

        // Equation (3.16) on page 65, MR 3rd print 2019
        matrix << x.normalized(), y.normalized(), z.normalized();

        return matrix;
    }

    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d& e)
    {
        // Appendix B.1 on page 577, MR 3rd print 2019
        Eigen::Matrix3d r_z = rotate_z(e[0]);
        Eigen::Matrix3d r_y = rotate_y(e[1]);
        Eigen::Matrix3d r_x = rotate_x(e[2]);

        return Eigen::Matrix3d::Identity() * r_z * r_y * r_x;
    }

    Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d& axis, double radians)
    {
        Eigen::Matrix3d matrix;
        double minus_cos = (1.0 - std::cos(radians));

        // Equations on page 72 and Equation (3.52) on page 84, MR 3rd print 2019
        matrix <<
            std::cos(radians) + std::pow(axis[0], 2) * minus_cos,
            axis[0] * axis[1] * minus_cos - axis[2] * std::sin(radians),
            axis[0] * axis[2] * minus_cos + axis[1] * std::sin(radians),

            axis[0] * axis[1] * minus_cos + axis[2] * std::sin(radians),
            std::cos(radians) + std::pow(axis[1], 2) * minus_cos,
            axis[1] * axis[2] * minus_cos - axis[0] * std::sin(radians),

            axis[0] * axis[2] * minus_cos - axis[1] * std::sin(radians),
            axis[1] * axis[2] * minus_cos + axis[0] * std::sin(radians),
            std::cos(radians) + std::pow(axis[2], 2) * minus_cos;

        return matrix;
    }

    Eigen::Matrix3d rotation_matrix(const Eigen::Matrix4d& tf)
    {
        // Equation (3.62) on page 87, MR 3rd print 2019
        return tf.block<3, 3>(0, 0);
    }

    Eigen::Matrix4d transformation_matrix(const Eigen::Vector3d& p)
    {
        return transformation_matrix(Eigen::Matrix3d::Identity(), p);
    }

    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d& r)
    {
        return transformation_matrix(r, Eigen::Vector3d::Zero());
    }

    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d& r, const Eigen::Vector3d& p)
    {
        Eigen::Matrix4d matrix = Eigen::Matrix4d::Zero();

        // Equation (3.62) on page 87, MR 3rd print 2019
        matrix.block<3, 3>(0, 0) = r;
        matrix.block<3, 1>(0, 3) = p;
        matrix(3, 3) = 1;

        return matrix;
    }
}
