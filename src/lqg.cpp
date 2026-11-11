#include "actlab/lqg.hpp"
#include "actlab/output_rls.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>

#include <Eigen/Eigenvalues>

namespace actlab {

using Eigen::MatrixXd;

MatrixXd expm(const MatrixXd& M) {
  const double norm = M.cwiseAbs().rowwise().sum().maxCoeff();
  const int s = norm > 0.5 ? static_cast<int>(std::ceil(std::log2(norm / 0.5))) : 0;
  const MatrixXd X = M / std::pow(2.0, s);
  MatrixXd E = MatrixXd::Identity(M.rows(), M.cols());
  MatrixXd term = E;
  for (int k = 1; k <= 14; ++k) {
    term = term * X / static_cast<double>(k);
    E += term;
  }
  for (int i = 0; i < s; ++i) E = E * E;
  return E;
}

void c2d_zoh(const MatrixXd& A, const MatrixXd& B, double T, MatrixXd& Ad, MatrixXd& Bd) {
  const Eigen::Index n = A.rows(), m = B.cols();
  MatrixXd M = MatrixXd::Zero(n + m, n + m);
  M.topLeftCorner(n, n) = A * T;
  M.topRightCorner(n, m) = B * T;
  const MatrixXd E = expm(M);
  Ad = E.topLeftCorner(n, n);
  Bd = E.topRightCorner(n, m);
}

MatrixXd dare(const MatrixXd& A, const MatrixXd& B, const MatrixXd& Q, const MatrixXd& R,
              int max_iter, double tol, int* iters) {
  const Eigen::Index n = A.rows();
  const MatrixXd I = MatrixXd::Identity(n, n);
  MatrixXd Ak = A;
  MatrixXd Gk = B * R.ldlt().solve(B.transpose());
  MatrixXd Hk = Q;
  int k = 0;
  for (; k < max_iter; ++k) {
    const Eigen::PartialPivLU<MatrixXd> W(I + Gk * Hk);
    const MatrixXd WiA = W.solve(Ak);
    const MatrixXd WiG = W.solve(Gk);
    const MatrixXd A1 = Ak * WiA;
    const MatrixXd G1 = Gk + Ak * WiG * Ak.transpose();
    const MatrixXd H1 = Hk + Ak.transpose() * Hk * WiA;
    const double change = (H1 - Hk).norm();
    Ak = A1;
    Gk = 0.5 * (G1 + G1.transpose());
    Hk = 0.5 * (H1 + H1.transpose());
    if (change <= tol * Hk.norm()) break;
  }
  if (iters) *iters = k + 1;
  return Hk;
}

TwoMassLinear TwoMassLinear::build(const NominalModel& m, double Jm_out, double bm_out, double d) {
  TwoMassLinear s;
  s.A.setZero();
  s.B.setZero();
  s.G.setZero();
  const double k = m.k, Jl = m.Jl;
  s.A(0, 1) = 1.0;
  s.A(1, 0) = -k / Jm_out;
  s.A(1, 1) = -(bm_out + d) / Jm_out;
  s.A(1, 2) = k / Jm_out;
  s.A(1, 3) = d / Jm_out;
  s.A(2, 3) = 1.0;
  s.A(3, 0) = k / Jl;
  s.A(3, 1) = d / Jl;
  s.A(3, 2) = -k / Jl;
  s.A(3, 3) = -(d + m.b_out) / Jl;
  s.A(3, 4) = 1.0 / Jl;
  s.B(1) = 1.0 / Jm_out;
  s.G(3) = 1.0 / Jl;
  return s;
}

namespace {

class LqgController final : public JointController {
 public:
  using V5 = Eigen::Matrix<double, 5, 1>;
  using M5 = Eigen::Matrix<double, 5, 5>;

  LqgController(const LqgConfig& c, const NominalModel& m, double f_vel, double f_pos)
      : c_(c), m_(m), ts_(1.0 / f_vel),
        rls_(OutputSideRls::Vec(m.Jl + m.N * m.N * c.Jm_motor, m.b_out + m.N * m.N * c.bm_motor, m.tau_c_out, m.mgl),
             f_pos) {
    synthesize();
  }

  std::string name() const override { return c_.adaptive ? "lqg+rls" : "lqg"; }

  void reset(const Reference& r0, const Feedback& fb) override {
    xh_.setZero();
    xh_(0) = fb.th_m / m_.N;
    xh_(2) = fb.th_l;
    // Gravity is fed through the known-input channel, so the load state starts at zero.
    z_ = 0.0;
    (void)r0;
    rls_.reset(fb.th_l, m_.mgl * std::sin(fb.th_l));
  }

  double update(const Reference& r, const Feedback& fb, bool pos_tick) override {
    // ---- Kalman correction with both encoders
    const Eigen::Vector2d y(fb.th_m / m_.N, fb.th_l);
    const Eigen::Vector2d innov(y(0) - xh_(0), y(1) - xh_(2));
    xh_ += L_ * innov;

    if (pos_tick && c_.adaptive) adapt(fb.th_l);

    // ---- reference state and feedforward from the (possibly re-identified) model
    const double tau_s = expected_link_torque(m_, r);
    const double Jm_out = m_.N * m_.N * c_.Jm_motor, bm_out = m_.N * m_.N * c_.bm_motor;
    double qm_ref = r.pos + tau_s / m_.k;
    double u_ff = Jm_out * r.acc + bm_out * r.vel + tau_s;
    if (c_.use_dist_ff) {  // cancel the estimated load: extra torque plus the extra twist it needs
      qm_ref -= xh_(4) / m_.k;
      u_ff -= xh_(4);
    }
    e_link_ = fb.th_l - r.pos;
    z_prev_ = z_;
    z_ = std::clamp(z_ + ts_ * e_link_, -c_.link_integrator_limit_rad_s, c_.link_integrator_limit_rad_s);
    const Eigen::Matrix<double, 5, 1> e(xh_(0) - qm_ref, xh_(1) - r.vel, xh_(2) - r.pos,
                                        xh_(3) - r.vel, z_);
    const double u = u_ff - K_.dot(e);
    return u / m_.N;
  }

  void commit(double tau_unsat, double tau_sat) override {
    if (std::abs(tau_sat - tau_unsat) > 1e-9) z_ = z_prev_;  // conditional integration
    const double u = m_.N * tau_sat;
    tau_applied_ = tau_sat;
    // ---- Kalman prediction with the torque actually applied and known gravity
    const double g = -m_.mgl * std::sin(xh_(2));
    xh_ = Ad_ * xh_ + Bd_ * u + Gd_ * g;
  }

  double diag(int i) const override {
    switch (i) {
      case 0: return xh_(4);
      case 1: return m_.Jl;
      case 2: return m_.mgl;
      case 3: return xh_(3);
      case 4: return z_;
      default: return 0.0;
    }
  }

 private:
  void synthesize() {
    const double Jm_out = m_.N * m_.N * c_.Jm_motor, bm_out = m_.N * m_.N * c_.bm_motor;
    const TwoMassLinear sys = TwoMassLinear::build(m_, Jm_out, bm_out, c_.gearbox_d);

    // ---- estimator model (5 states incl. load torque), ZOH at the loop rate
    Eigen::MatrixXd Bc(5, 2);
    Bc << sys.B, sys.G;
    Eigen::MatrixXd Ad, Bd;
    c2d_zoh(sys.A, Bc, ts_, Ad, Bd);
    Ad_ = Ad;
    Bd_ = Bd.col(0);
    Gd_ = Bd.col(1);

    Eigen::Matrix<double, 2, 5> C;
    C.setZero();
    C(0, 0) = 1.0;
    C(1, 2) = 1.0;
    Eigen::Matrix<double, 5, 3> Bw;
    Bw.setZero();
    Bw(1, 0) = 1.0 / Jm_out;
    Bw(3, 1) = 1.0 / m_.Jl;
    Bw(4, 2) = 1.0;
    const Eigen::Vector3d qc(c_.q_tau_m, c_.q_tau_l, c_.q_dist);
    const Eigen::MatrixXd Qd = Bw * qc.asDiagonal() * Bw.transpose() * ts_;
    // Uniform quantization noise, LSB^2 / 12, with the motor encoder referred to the output. The
    // small floor keeps Rn positive definite if an encoder is configured as ideal.
    constexpr double kMinVarianceRad2 = 1e-13;
    const double lsb_m = 2.0 * kPi / std::pow(2.0, c_.encoder_motor_bits) / m_.N;
    const double lsb_l = 2.0 * kPi / std::pow(2.0, c_.encoder_link_bits);
    Eigen::Matrix2d Rn = Eigen::Vector2d(lsb_m * lsb_m / 12.0 + kMinVarianceRad2, lsb_l * lsb_l / 12.0 + kMinVarianceRad2).asDiagonal();
    const Eigen::MatrixXd P = dare(Ad.transpose(), C.transpose(), Qd, Rn, 100, 1e-11, &iters_est_);
    const Eigen::MatrixXd S = C * P * C.transpose() + Rn;
    L_ = P * C.transpose() * S.inverse();

    // ---- LQR on [4 plant states, integral of link error]
    Eigen::MatrixXd Ac = Eigen::MatrixXd::Zero(5, 5), Bu = Eigen::MatrixXd::Zero(5, 1);
    Ac.topLeftCorner(4, 4) = sys.A.topLeftCorner(4, 4);
    Ac(4, 2) = 1.0;
    Bu.topRows(4) = sys.B.topRows(4);
    Eigen::MatrixXd Adc, Bdc;
    c2d_zoh(Ac, Bu, ts_, Adc, Bdc);
    const Eigen::VectorXd qdiag = (Eigen::VectorXd(5) << c_.w_qm, c_.w_vm, c_.w_ql, c_.w_vl, c_.w_int).finished();
    const Eigen::MatrixXd Q = qdiag.asDiagonal();
    const Eigen::MatrixXd R = Eigen::MatrixXd::Constant(1, 1, c_.r_u);
    const Eigen::MatrixXd Pc = dare(Adc, Bdc, Q, R, 100, 1e-11, &iters_ctrl_);
    Adc_ = Adc;
    Bdc_ = Bdc;
    C_ = C;
    const Eigen::MatrixXd Kc = (R + Bdc.transpose() * Pc * Bdc).ldlt().solve(Bdc.transpose() * Pc * Adc);
    K_ = Kc.row(0).transpose();
    synth_Jl_ = m_.Jl;
  }

  // Plausibility bounds on what the estimate may write into the model, kg m^2 and N m.
  static constexpr double kMinLinkInertia = 0.1, kMaxLinkInertia = 5.0, kMaxGravityTorque = 60.0;
  // Re-synthesize LQR + KF when the inertia estimate has moved by more than this fraction.
  static constexpr double kResynthesisFraction = 0.05;

  void adapt(double th_l) {
    if (!rls_.update(th_l, m_.N * tau_applied_)) return;
    const auto& th = rls_.theta();
    const double Jm_out = m_.N * m_.N * c_.Jm_motor;
    m_.Jl = std::clamp(th(0) - Jm_out, kMinLinkInertia, kMaxLinkInertia);
    m_.mgl = std::clamp(th(3), 0.0, kMaxGravityTorque);
    if (std::abs(m_.Jl - synth_Jl_) > kResynthesisFraction * synth_Jl_) synthesize();
  }

  friend LqgReport actlab::lqg_report(const LqgConfig&, const NominalModel&, double);
  LqgConfig c_;
  NominalModel m_;
  double ts_;
  M5 Ad_;
  V5 Bd_, Gd_, K_;
  Eigen::Matrix<double, 5, 2> L_;
  V5 xh_ = V5::Zero();
  double z_{0}, z_prev_{0}, e_link_{0}, tau_applied_{0}, synth_Jl_{0};
  int iters_ctrl_{0}, iters_est_{0};
  Eigen::MatrixXd Adc_, Bdc_, C_;
  OutputSideRls rls_;
};

Eigen::VectorXcd to_continuous(const Eigen::VectorXcd& z, double T) {
  Eigen::VectorXcd s(z.size());
  for (Eigen::Index i = 0; i < z.size(); ++i) s(i) = std::log(z(i)) / T;
  return s;
}

}  // namespace

LqgReport lqg_report(const LqgConfig& lc, const NominalModel& m, double f_vel) {
  const auto t0 = std::chrono::steady_clock::now();
  LqgController c(lc, m, f_vel, 1000.0);
  const auto t1 = std::chrono::steady_clock::now();
  LqgReport r;
  r.synth_time_us = std::chrono::duration<double, std::micro>(t1 - t0).count();
  r.K = c.K_;
  r.L = c.L_;
  const double T = 1.0 / f_vel;
  const Eigen::MatrixXd Aol = c.Adc_.topLeftCorner(4, 4);
  r.poles_open = to_continuous(Eigen::EigenSolver<Eigen::MatrixXd>(Aol).eigenvalues(), T);
  const Eigen::MatrixXd Acl = c.Adc_ - c.Bdc_ * c.K_.transpose();
  r.poles_ctrl = to_continuous(Eigen::EigenSolver<Eigen::MatrixXd>(Acl).eigenvalues(), T);
  const Eigen::MatrixXd Ae = c.Ad_ - c.L_ * c.C_ * c.Ad_;
  r.poles_est = to_continuous(Eigen::EigenSolver<Eigen::MatrixXd>(Ae).eigenvalues(), T);
  r.dare_iters_ctrl = c.iters_ctrl_;
  r.dare_iters_est = c.iters_est_;
  return r;
}

std::unique_ptr<JointController> make_lqg(const LqgConfig& lc, const NominalModel& m, double f_vel,
                                          double f_pos) {
  return std::make_unique<LqgController>(lc, m, f_vel, f_pos);
}

}  // namespace actlab
