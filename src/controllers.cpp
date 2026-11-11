#include "actlab/controllers.hpp"
#include "actlab/lqg.hpp"
#include "actlab/output_rls.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace actlab {

double expected_link_torque(const NominalModel& m, const Reference& r) {
  constexpr double kCoulombSmoothingRadS = 0.005;  // same width as the plant's link friction
  return m.Jl * r.acc + m.b_out * r.vel + m.tau_c_out * std::tanh(r.vel / kCoulombSmoothingRadS) +
         m.mgl * std::sin(r.pos);
}

LinkCorrector::LinkCorrector(const LinkLoopConfig& c, double fs) : c_(c), ts_(1.0 / fs) {
  if (c_.notch) notch_ = Biquad<double>::notch(c_.notch_f0, c_.notch_zn, c_.notch_zd, fs);
}

void LinkCorrector::set_notch(double f0) {
  if (!c_.notch) return;
  c_.notch_f0 = f0;
  notch_ = Biquad<double>::notch(f0, c_.notch_zn, c_.notch_zd, 1.0 / ts_);
  notch_.reset(last_err_);  // start the new filter in steady state, no step transient
}

void LinkCorrector::update(double err_link) {
  last_err_ = err_link;
  const double e = notch_.step(err_link);
  // Conditional integration: stop integrating when the output is clamped in that direction.
  const double u = c_.kp * e + integ_;
  const bool clamped = std::abs(u) >= c_.limit && (u * e > 0.0);
  if (!clamped) integ_ = std::clamp(integ_ + ts_ * c_.ki * e, -c_.limit, c_.limit);
  corr_ = std::clamp(c_.kp * e + integ_, -c_.limit, c_.limit);
}

double PositionLoop::motor_ref(const Reference& r) const {
  // The gearbox twists by tau_s / k under load. Commanding the motor to r would leave the
  // link short by that deflection, so shift the motor target by the expected twist.
  const double defl = defl_ ? expected_link_torque(m_, r) / m_.k : 0.0;
  return m_.N * (r.pos + defl + link_correction());
}

void PositionLoop::update_link(const Reference& r, double th_l) {
  if (link_on_) link_.update(r.pos - th_l);
}

VelSetpoint PositionLoop::update(const Reference& r, const Feedback& fb) {
  update_link(r, fb.th_l);
  VelSetpoint sp;
  // Clamp just above the speed the bus voltage can support (~415 rad/s at 48 V, no load).
  // Without it a saturated position loop asks for speeds the drive cannot reach and the
  // nested integrators wind up into a bang-bang limit cycle.
  constexpr double kMaxMotorSpeedRadS = 420.0;
  sp.w_ref = std::clamp(kp_ * (motor_ref(r) - fb.th_m) + m_.N * r.vel, -kMaxMotorSpeedRadS, kMaxMotorSpeedRadS);
  sp.tau_ff = m_.Jm * m_.N * r.acc + expected_link_torque(m_, r) / m_.N;
  return sp;
}

VelocityLoop::VelocityLoop(const VelocityLoopConfig& c, double J_nominal, double fs)
    : c_(c), fs_(fs), J_(J_nominal) {
  const double wc = 2.0 * kPi * c_.bw_hz;
  const double kp = J_ * wc;
  const double ki = kp * c_.pi_zero_ratio * wc;
  pi_ = PI(kp, ki, 1.0 / fs_, ki / kp);
  notch_ = c_.notch ? Biquad<double>::notch(c_.notch_f0, c_.notch_zn, c_.notch_zd, fs_)
                    : Biquad<double>::identity();
  dob_lpf_ = LowPass1(c_.dob_hz, fs_);
}

void VelocityLoop::set_inertia(double J) {
  J_ = J;
  const double wc = 2.0 * kPi * c_.bw_hz;
  const double kp = J_ * wc;
  pi_.set_gains(kp, kp * c_.pi_zero_ratio * wc);
}

double VelocityLoop::update(const VelSetpoint& sp, double w_hat) {
  e_last_ = sp.w_ref - w_hat;
  w_last_ = w_hat;
  double u = pi_.output(e_last_) + sp.tau_ff;
  if (c_.dob) u -= d_hat_;
  return notch_.step(u);
}

void VelocityLoop::commit(double tau_unsat, double tau_sat) {
  pi_.commit(e_last_, tau_unsat, tau_sat);
  if (c_.dob) {
    // d_hat = Q(J s w - tau) written without differentiating w:
    // d_hat = J wq w - LPF(tau + J wq w), with Q = wq / (s + wq).
    const double wq = 2.0 * kPi * c_.dob_hz;
    d_hat_ = J_ * wq * w_last_ - dob_lpf_.step(tau_sat + J_ * wq * w_last_);
  }
}

namespace {

// ------------------------------------------------------------------------- cascade (+ adaptive)
class CascadeController final : public JointController {
 public:
  CascadeController(const ControllerConfig& c, const NominalModel& m, double f_vel, double f_pos)
      : cfg_(c),
        pos_(c.pos_kp_ratio * 2.0 * kPi * c.vel.bw_hz, m, c.deflection_comp, c.link, f_pos),
        vel_(c.vel, m.J_tot_motor(), f_vel),
        rls_(OutputSideRls::Vec(m.Jl + m.N * m.N * m.Jm, m.b_out, m.tau_c_out, m.mgl), f_pos) {}

  std::string name() const override {
    std::string n = "cascade";
    if (cfg_.link.enabled) n += cfg_.link.notch ? "+link(notch)" : "+link";
    if (cfg_.vel.notch) n += "+velnotch";
    if (cfg_.vel.dob) n += "+dob";
    if (cfg_.adaptive) n += "+rls";
    return n;
  }

  void reset(const Reference& r0, const Feedback& fb) override {
    sp_ = pos_.update(r0, fb);
    rls_.reset(fb.th_l, pos_.model().N * sp_.tau_ff);
  }

  double update(const Reference& r, const Feedback& fb, bool pos_tick) override {
    if (pos_tick) {
      if (cfg_.adaptive) adapt(fb);
      sp_ = pos_.update(r, fb);
    }
    return vel_.update(sp_, fb.w_m);
  }

  void commit(double tau_unsat, double tau_sat) override {
    vel_.commit(tau_unsat, tau_sat);
    tau_applied_ = tau_sat;
  }

  double diag(int i) const override {
    if (cfg_.adaptive && i >= 0 && i < 4) return rls_.theta()(i);
    if (i == 4) return vel_.dob_estimate();
    return 0.0;
  }

 private:
  // Plausibility bounds on what the estimate may write into the model, kg m^2 and N m.
  static constexpr double kMinOutputInertia = 0.2, kMaxOutputInertia = 5.0, kMinLinkInertia = 0.05;
  static constexpr double kMaxGravityTorque = 60.0;
  static constexpr double kNotchRecentreFraction = 0.05;

  void adapt(const Feedback& fb) {
    NominalModel& m = pos_.model();
    if (!rls_.update(fb.th_l, m.N * tau_applied_)) return;
    const auto& th = rls_.theta();
    const double J_eq = std::clamp(th(0), kMinOutputInertia, kMaxOutputInertia);
    m.Jl = std::max(kMinLinkInertia, J_eq - m.N * m.N * m.Jm);
    // Viscous and Coulomb regressors are nearly collinear on back-and-forth moves
    // (tanh(v/eps) ~ sign(v) and v always share a sign), so their individual estimates
    // wander even when the fit is good. Only inertia and gravity are applied.
    m.mgl = std::clamp(th(3), 0.0, kMaxGravityTorque);
    vel_.set_inertia(m.J_tot_motor());
    // The locked-motor mode moves with the payload: f = sqrt(k / Jl) / 2 pi. Re-centre the
    // link-loop notch when the estimate has moved by more than 5 %.
    const double f_new = std::sqrt(m.k / m.Jl) / (2.0 * kPi);
    if (std::abs(f_new - notch_f_) > kNotchRecentreFraction * notch_f_) {
      pos_.set_link_notch(f_new);
      notch_f_ = f_new;
    }
  }

  ControllerConfig cfg_;
  PositionLoop pos_;
  VelocityLoop vel_;
  VelSetpoint sp_{};
  double tau_applied_{0};
  double notch_f_{cfg_.link.notch_f0};
  OutputSideRls rls_;
};

// ------------------------------------------------------------------------- ADRC
// Linear ADRC (Gao 2003): treat everything except b0*u as a lumped "total disturbance" f,
// estimate it with a 3rd-order extended state observer and cancel it.
//   y'' = f + b0 u,  y = th_m / N (motor angle expressed at the output),  b0 = 1/(N J_tot)
class AdrcController final : public JointController {
 public:
  AdrcController(const ControllerConfig& c, const NominalModel& m, double f_vel, double f_pos)
      : cfg_(c), pos_(0.0, m, c.deflection_comp, c.link, f_pos), ts_(1.0 / f_vel) {
    b0_ = 1.0 / (m.N * m.J_tot_motor());
    wc_ = 2.0 * kPi * c.adrc_wc_hz;
    wo_ = 2.0 * kPi * c.adrc_wo_hz;
    notch_ = c.vel.notch ? Biquad<double>::notch(c.vel.notch_f0, c.vel.notch_zn, c.vel.notch_zd, f_vel)
                         : Biquad<double>::identity();
  }
  std::string name() const override { return "adrc"; }

  void reset(const Reference& r0, const Feedback& fb) override {
    z1_ = fb.th_m / pos_.model().N;
    z2_ = 0.0;
    // Start the disturbance estimate at the static gravity load so the first sample
    // does not see a large step.
    z3_ = -b0_ * expected_link_torque(pos_.model(), r0) / pos_.model().N;
  }

  double update(const Reference& r, const Feedback& fb, bool pos_tick) override {
    if (pos_tick) pos_.update_link(r, fb.th_l);
    const double y = fb.th_m / pos_.model().N;
    const double b1 = 3.0 * wo_, b2 = 3.0 * wo_ * wo_, b3 = wo_ * wo_ * wo_;
    const int n = 4;
    const double h = ts_ / n;
    for (int i = 0; i < n; ++i) {  // Euler sub-steps keep wo*h << 1
      const double e = y - z1_;
      const double d1 = z2_ + b1 * e;
      const double d2 = z3_ + b0_ * u_applied_ + b2 * e;
      const double d3 = b3 * e;
      z1_ += h * d1;
      z2_ += h * d2;
      z3_ += h * d3;
    }
    const double ry = pos_.motor_ref(r) / pos_.model().N;
    const double kp = wc_ * wc_, kd = 2.0 * wc_;
    const double u0 = kp * (ry - z1_) + kd * (r.vel - z2_) + r.acc;
    return notch_.step((u0 - z3_) / b0_);
  }

  void commit(double, double tau_sat) override { u_applied_ = tau_sat; }
  double diag(int i) const override { return i == 0 ? z3_ : 0.0; }

 private:
  ControllerConfig cfg_;
  PositionLoop pos_;  // only used for its deflection-compensated reference
  double ts_, b0_, wc_, wo_;
  double z1_{0}, z2_{0}, z3_{0}, u_applied_{0};
  Biquad<double> notch_;
};

// ------------------------------------------------------------------------- super-twisting SMC
// s = e' + lambda e. Super-twisting puts the discontinuity under an integral, so the
// torque command is continuous and the classic first-order SMC chatter mostly goes away.
class StsmcController final : public JointController {
 public:
  StsmcController(const ControllerConfig& c, const NominalModel& m, double f_vel, double f_pos)
      : cfg_(c), pos_(0.0, m, c.deflection_comp, c.link, f_pos), ts_(1.0 / f_vel) {
    b0_ = 1.0 / (m.N * m.J_tot_motor());
    notch_ = c.vel.notch ? Biquad<double>::notch(c.vel.notch_f0, c.vel.notch_zn, c.vel.notch_zd, f_vel)
                         : Biquad<double>::identity();
  }
  std::string name() const override { return "stsmc"; }

  void reset(const Reference& r0, const Feedback&) override {
    v_ = -b0_ * expected_link_torque(pos_.model(), r0) / pos_.model().N;
  }

  double update(const Reference& r, const Feedback& fb, bool pos_tick) override {
    if (pos_tick) pos_.update_link(r, fb.th_l);
    const double N = pos_.model().N;
    const double e = fb.th_m / N - pos_.motor_ref(r) / N;
    const double ed = fb.w_m / N - r.vel;
    const double s = ed + cfg_.smc_lambda * e;
    sgn_ = s / (std::abs(s) + cfg_.smc_phi);  // boundary layer only to tame encoder noise
    const double a = r.acc - cfg_.smc_lambda * ed - cfg_.smc_k1 * std::sqrt(std::abs(s)) * sgn_ + v_;
    // Expected (nominal) load is fed forward so the integral term only learns the mismatch.
    const double tau_ff = expected_link_torque(pos_.model(), r) / N;
    return notch_.step(a / b0_ + tau_ff - pos_.model().Jl * r.acc / N);
  }

  void commit(double tau_unsat, double tau_sat) override {
    // Freeze the twisting integrator while the actuator is saturated.
    if (std::abs(tau_sat - tau_unsat) < 1e-9) v_ += ts_ * (-cfg_.smc_k2 * sgn_);
  }
  double diag(int i) const override { return i == 0 ? v_ : 0.0; }

 private:
  ControllerConfig cfg_;
  PositionLoop pos_;
  double ts_, b0_;
  double v_{0}, sgn_{0};
  Biquad<double> notch_;
};

}  // namespace

std::unique_ptr<JointController> make_controller(const ControllerConfig& cfg,
                                                 const NominalModel& nominal, double f_vel,
                                                 double f_pos) {
  if (cfg.type == "cascade") return std::make_unique<CascadeController>(cfg, nominal, f_vel, f_pos);
  if (cfg.type == "lqg") return make_lqg(cfg.lqg, nominal, f_vel, f_pos);
  if (cfg.type == "adrc") return std::make_unique<AdrcController>(cfg, nominal, f_vel, f_pos);
  if (cfg.type == "stsmc") return std::make_unique<StsmcController>(cfg, nominal, f_vel, f_pos);
  throw std::invalid_argument("unknown controller type: " + cfg.type);
}

}  // namespace actlab
