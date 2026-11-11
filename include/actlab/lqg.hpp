#pragma once
#include <Eigen/Dense>

#include "actlab/controllers.hpp"

namespace actlab {

// ---------------------------------------------------------------- linear algebra utilities
// Zero-order-hold discretization: [Ad Bd; 0 I] = expm([A B; 0 0] * T).
// Matrix exponential by scaling and squaring with a 14-term Taylor series.
Eigen::MatrixXd expm(const Eigen::MatrixXd& M);
void c2d_zoh(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, double T, Eigen::MatrixXd& Ad,
             Eigen::MatrixXd& Bd);

// Discrete algebraic Riccati equation  P = A'PA - A'PB (R + B'PB)^-1 B'PA + Q
// solved with the structure-preserving doubling algorithm (quadratic convergence, typically
// < 30 iterations). Plain Riccati iteration needs tens of thousands of steps here because
// the closed-loop poles of a 5 kHz design sit very close to z = 1.
Eigen::MatrixXd dare(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& Q,
                     const Eigen::MatrixXd& R, int max_iter = 100, double tol = 1e-11,
                     int* iters = nullptr);

// ---------------------------------------------------------------- two-mass model
// Output-side coordinates: q_m = th_m / N, u = N * tau_m (motor torque seen at the output).
// x = [q_m, q_m', q_l, q_l', tau_d], tau_d = unknown link-side load (random walk in the KF).
struct TwoMassLinear {
  Eigen::Matrix<double, 5, 5> A;  // continuous
  Eigen::Matrix<double, 5, 1> B;
  Eigen::Matrix<double, 5, 1> G;  // known nonlinear link-side input (gravity, friction)
  static TwoMassLinear build(const NominalModel& m, double Jm_out, double bm_out, double d);
};

struct LqgReport {
  Eigen::VectorXd K;               // LQR gain on [q_m, q_m', q_l, q_l', int e_l]
  Eigen::MatrixXd L;               // Kalman gain (5x2)
  Eigen::VectorXcd poles_open;     // continuous-equivalent open-loop plant poles (4)
  Eigen::VectorXcd poles_ctrl;     // state-feedback closed loop (5)
  Eigen::VectorXcd poles_est;      // estimator error dynamics (5)
  int dare_iters_ctrl{0}, dare_iters_est{0};
  double synth_time_us{0};
};
LqgReport lqg_report(const LqgConfig& lc, const NominalModel& m, double f_vel);

std::unique_ptr<JointController> make_lqg(const LqgConfig& lc, const NominalModel& m, double f_vel,
                                          double f_pos);

}  // namespace actlab
