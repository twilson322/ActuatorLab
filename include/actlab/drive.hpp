#pragma once
#include "actlab/filters.hpp"
#include "actlab/foc.hpp"
#include "actlab/params.hpp"
#include "actlab/plant.hpp"

namespace actlab {

// Second-order tracking observer (a PLL in disguise) that turns a quantized encoder angle
// into a usable speed estimate. A raw finite difference of a 14-bit encoder at 20 kHz has a
// quantization floor of 2*pi/16384 * 20000 = 7.7 rad/s, which is useless for a speed loop.
class VelocityObserver {
 public:
  VelocityObserver() = default;
  VelocityObserver(double bw_hz, double fs) : ts_(1.0 / fs) {
    const double w = 2.0 * kPi * bw_hz;
    l1_ = 2.0 * w;  // critically damped: zeta = 1
    l2_ = w * w;
  }
  void reset(double th) { th_ = th; w_ = 0.0; }
  void step(double th_meas) {
    const double e = th_meas - th_;
    th_ += ts_ * (w_ + l1_ * e);
    w_ += ts_ * l2_ * e;
  }
  double theta() const { return th_; }
  double omega() const { return w_; }

 private:
  double ts_{5e-5}, l1_{0}, l2_{0};
  double th_{0}, w_{0};
};

struct DriveConfig {
  double cur_bw_hz{1000.0};    // current-loop bandwidth (pole-zero cancellation design)
  double obs_bw_hz{400.0};     // speed observer bandwidth
  double vlimit_frac{0.96};    // fraction of the SVPWM linear limit Vdc/sqrt(3) to use
  bool delay_comp{true};       // advance the inverse-Park angle by 1.5 PWM periods
};

// The 20 kHz FOC current loop as it would run on the drive MCU:
// sense -> Clarke -> Park -> dq PI with decoupling -> circle limit -> inverse Park -> SVPWM.
// The computed voltage is applied in the *next* PWM period (one-sample transport delay).
class Drive {
 public:
  Drive(const MotorParams& nominal_motor, const DriveParams& dp, const DriveConfig& cfg,
        double f_cur);
  void reset(const Sensors& s);
  // Split into two calls so the speed loop can run between sensing and current control,
  // which is how an interrupt-driven drive orders the work inside one PWM period.
  void sense(const Sensors& s);
  foc::AlphaBeta control(double iq_ref);

  double id() const { return id_; }
  double iq() const { return iq_; }
  double w_m() const { return obs_.omega(); }
  double th_m() const { return th_m_; }
  bool saturated() const { return sat_; }
  double vq() const { return vq_; }

 private:
  MotorParams m_;
  DriveParams dp_;
  DriveConfig cfg_;
  double ts_;
  PI pi_d_, pi_q_;
  VelocityObserver obs_;
  double id_{0}, iq_{0}, th_m_{0}, vq_{0};
  bool sat_{false};
};

}  // namespace actlab
