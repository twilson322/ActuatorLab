#pragma once
#include <memory>
#include <string>

#include "actlab/filters.hpp"
#include "actlab/params.hpp"
#include "actlab/rls.hpp"

namespace actlab {

// What the control engineer believes about the joint. Deliberately separate from the
// PlantParams the simulator uses, so model mismatch is explicit.
struct NominalModel {
  double N{50}, kt{0.1}, Jm{6e-5}, Jl{0.6}, k{4000}, mgl{11.8}, b_out{0.05}, tau_c_out{0.5};
  double J_tot_motor() const { return Jm + Jl / (N * N); }  // rigid, reflected to motor
  static NominalModel from(const PlantParams& p) {
    return {p.joint.N, p.motor.kt(), p.motor.Jm, p.joint.Jl, p.joint.k,
            p.joint.mgl, p.joint.bl, p.joint.tau_c_l};
  }
};

struct Reference { double pos{0}, vel{0}, acc{0}; };  // link side
struct Feedback {
  double th_m{0}, w_m{0};  // motor side (encoder + velocity observer)
  double th_l{0};          // link encoder
};
struct VelSetpoint { double w_ref{0}, tau_ff{0}; };

// Expected spring (gearbox output) torque for a reference, used for deflection and
// torque feedforward. tanh knee matches a few mrad/s, not the plant's internal value.
double expected_link_torque(const NominalModel& m, const Reference& r);

// ---------------------------------------------------------------- link-side correction
// Dual-encoder outer loop. The motor-side loops cannot see gearbox wind-up caused by loads
// the model does not know about, so a slow PI on the *link* encoder error shifts the motor
// target. This loop is non-collocated: between motor and link sits the lightly damped
// locked-motor mode at sqrt(k/Jl) (the anti-resonance of the motor-side FRF), which caps the
// usable link gain. An optional notch on the link error targets exactly that mode.
struct LinkLoopConfig {
  bool enabled{false};
  double kp{0.0};      // rad of motor-reference shift (link units) per rad of link error
  double ki{0.0};      // 1/s
  bool notch{false};
  double notch_f0{13.0}, notch_zn{0.05}, notch_zd{0.7};
  double limit{0.05};  // clamp on the correction [rad]
};

class LinkCorrector {
 public:
  LinkCorrector() = default;
  LinkCorrector(const LinkLoopConfig& c, double fs);
  void update(double err_link);
  void set_notch(double f0);
  double correction() const { return corr_; }

 private:
  LinkLoopConfig c_{};
  double ts_{1e-3};
  Biquad<double> notch_ = Biquad<double>::identity();
  double integ_{0.0}, corr_{0.0}, last_err_{0.0};
};

// ---------------------------------------------------------------- position loop (host, 1 kHz)
class PositionLoop {
 public:
  PositionLoop(double kp, const NominalModel& m, bool deflection_comp,
               const LinkLoopConfig& link = {}, double f_pos = 1000.0)
      : kp_(kp), m_(m), defl_(deflection_comp), link_(link, f_pos), link_on_(link.enabled) {}
  VelSetpoint update(const Reference& r, const Feedback& fb);  // also steps the link loop
  void update_link(const Reference& r, double th_l);
  double motor_ref(const Reference& r) const;  // N * (r + expected deflection + link correction)
  NominalModel& model() { return m_; }
  const NominalModel& model() const { return m_; }
  double link_correction() const { return link_on_ ? link_.correction() : 0.0; }
  void set_link_notch(double f0) { link_.set_notch(f0); }

 private:
  double kp_;
  NominalModel m_;
  bool defl_;
  LinkCorrector link_;
  bool link_on_;
};

// ---------------------------------------------------------------- velocity loop (drive, 5 kHz)
struct VelocityLoopConfig {
  double bw_hz{40.0};          // target crossover of the rigid-body loop
  double pi_zero_ratio{0.2};   // PI zero at ratio * crossover
  bool notch{false};
  double notch_f0{29.0}, notch_zn{0.05}, notch_zd{0.6};
  bool dob{false};
  double dob_hz{15.0};
};

class VelocityLoop {
 public:
  VelocityLoop(const VelocityLoopConfig& c, double J_nominal, double fs);
  double update(const VelSetpoint& sp, double w_hat);  // returns unsaturated torque
  void commit(double tau_unsat, double tau_sat);       // anti-windup + DOB bookkeeping
  void set_inertia(double J);                          // gain scheduling hook
  double dob_estimate() const { return d_hat_; }

 private:
  VelocityLoopConfig c_;
  double fs_;
  double J_;
  PI pi_;
  Biquad<double> notch_;
  LowPass1 dob_lpf_;
  double e_last_{0}, w_last_{0}, d_hat_{0};
};

// ---------------------------------------------------------------- unified controller interface
class JointController {
 public:
  virtual ~JointController() = default;
  virtual std::string name() const = 0;
  virtual void reset(const Reference& r0, const Feedback& fb0) = 0;
  // Runs at f_vel. pos_tick marks the slower f_pos sub-rate.
  virtual double update(const Reference& r, const Feedback& fb, bool pos_tick) = 0;
  virtual void commit(double tau_unsat, double tau_sat) = 0;
  // Optional diagnostics for logging (estimated inertia etc.).
  virtual double diag(int /*i*/) const { return 0.0; }
};

// Two-mass LQG (see lqg.hpp). Lives here so ControllerConfig can carry it.
struct LqgConfig {
  // LQR weights on [q_m err, q_m' err, q_l err, q_l' err, integral of link err], R on u.
  double w_qm{0.0}, w_vm{0.0}, w_ql{4e6}, w_vl{2e3}, w_int{4e7}, r_u{1.0};
  // Kalman filter: process noise on motor/link accelerations [Nm]^2 and load random walk.
  double q_tau_m{0.5}, q_tau_l{4.0}, q_dist{4e4};
  bool use_dist_ff{true};
  bool adaptive{false};  // RLS on payload + periodic re-synthesis of LQR and KF
  double gearbox_d{2.0}, Jm_motor{6e-5}, bm_motor{2e-5};
  // Encoder resolutions the Kalman filter's quantization noise is computed from. They default to
  // the drive's encoders. They were hard-coded as 2^14 and 2^17 inside the filter design before,
  // so a different encoder in DriveParams would have left the filter tuned for the old one.
  int encoder_motor_bits{DriveParams{}.enc_motor_bits}, encoder_link_bits{DriveParams{}.enc_link_bits};
  double link_integrator_limit_rad_s{0.02};  // anti-windup clamp on the integral of link error
};

struct ControllerConfig {
  std::string type{"cascade"};  // cascade | adrc | stsmc | lqg
  VelocityLoopConfig vel{};
  double pos_kp_ratio{0.25};    // position P gain = ratio * velocity design bandwidth (rad/s)
  bool deflection_comp{true};
  LinkLoopConfig link{};
  bool adaptive{false};         // RLS on the output-side rigid model (cascade only)
  LqgConfig lqg{};
  // ADRC
  double adrc_wc_hz{8.0}, adrc_wo_hz{40.0};
  // Super-twisting SMC
  double smc_lambda{60.0}, smc_k1{25.0}, smc_k2{400.0}, smc_phi{2e-3};
};

std::unique_ptr<JointController> make_controller(const ControllerConfig& cfg,
                                                 const NominalModel& nominal, double f_vel,
                                                 double f_pos);

}  // namespace actlab
