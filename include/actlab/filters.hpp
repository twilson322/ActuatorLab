#pragma once
#include <algorithm>
#include <cmath>

namespace actlab {

inline constexpr double kPi = 3.14159265358979323846;

// Direct-form II transposed biquad. Templated so the same filter can run in float on an
// MCU target or double in simulation.
template <typename T = double>
class Biquad {
 public:
  Biquad() = default;
  Biquad(T b0, T b1, T b2, T a1, T a2) : b0_(b0), b1_(b1), b2_(b2), a1_(a1), a2_(a2) {}

  T step(T x) {
    const T y = b0_ * x + z1_;
    z1_ = b1_ * x - a1_ * y + z2_;
    z2_ = b2_ * x - a2_ * y;
    return y;
  }
  void reset(T x0 = T(0)) {
    // Initialize to the steady state for a constant input x0 (DC gain applied).
    const T dc = (b0_ + b1_ + b2_) / (T(1) + a1_ + a2_);
    const T y0 = dc * x0;
    z2_ = b2_ * x0 - a2_ * y0;
    z1_ = b1_ * x0 - a1_ * y0 + z2_;
  }

  // Generic second-order section H(s) = (s^2 + 2 zn w0 s + w0^2) / (s^2 + 2 zd w0 s + w0^2)
  // discretized with Tustin, prewarped so the notch centre lands exactly on f0.
  // Depth at f0 is zn/zd (e.g. zn=0.02, zd=0.4 -> -26 dB).
  static Biquad notch(double f0, double zeta_num, double zeta_den, double fs) {
    const double w0 = 2.0 * kPi * f0;
    const double K = w0 / std::tan(w0 / (2.0 * fs));  // prewarped Tustin constant
    const double K2 = K * K, w02 = w0 * w0;
    const double a0 = K2 + 2.0 * zeta_den * w0 * K + w02;
    const double b0 = (K2 + 2.0 * zeta_num * w0 * K + w02) / a0;
    const double b1 = (2.0 * w02 - 2.0 * K2) / a0;
    const double b2 = (K2 - 2.0 * zeta_num * w0 * K + w02) / a0;
    const double a1 = (2.0 * w02 - 2.0 * K2) / a0;
    const double a2 = (K2 - 2.0 * zeta_den * w0 * K + w02) / a0;
    return Biquad(T(b0), T(b1), T(b2), T(a1), T(a2));
  }

  // Second-order Butterworth-style low-pass, Tustin with prewarp.
  static Biquad lowpass(double fc, double zeta, double fs) {
    const double w0 = 2.0 * kPi * fc;
    const double K = w0 / std::tan(w0 / (2.0 * fs));
    const double K2 = K * K, w02 = w0 * w0;
    const double a0 = K2 + 2.0 * zeta * w0 * K + w02;
    return Biquad(T(w02 / a0), T(2.0 * w02 / a0), T(w02 / a0), T((2.0 * w02 - 2.0 * K2) / a0),
                  T((K2 - 2.0 * zeta * w0 * K + w02) / a0));
  }

  static Biquad identity() { return Biquad(T(1), T(0), T(0), T(0), T(0)); }

 private:
  T b0_{1}, b1_{0}, b2_{0}, a1_{0}, a2_{0};
  T z1_{0}, z2_{0};
};

// First-order low-pass via exact zero-order-hold discretization.
class LowPass1 {
 public:
  LowPass1() = default;
  LowPass1(double fc, double fs) : a_(std::exp(-2.0 * kPi * fc / fs)) {}
  double step(double x) {
    y_ = a_ * y_ + (1.0 - a_) * x;
    return y_;
  }
  void reset(double y0) { y_ = y0; }
  double value() const { return y_; }

 private:
  double a_{0.0};
  double y_{0.0};
};

// Discrete PI with back-calculation anti-windup. The caller computes the unsaturated
// output, applies whatever saturation is physically present (scalar clamp, voltage circle,
// current limit) and then calls commit() with both values so the integrator unwinds by
// exactly the amount of saturation.
class PI {
 public:
  PI() = default;
  PI(double kp, double ki, double ts, double kaw) : kp_(kp), ki_(ki), ts_(ts), kaw_(kaw) {}

  double output(double e) const { return kp_ * e + integ_; }
  void commit(double e, double u_unsat, double u_sat) {
    integ_ += ts_ * (ki_ * e + kaw_ * (u_sat - u_unsat));
  }
  void reset(double integ = 0.0) { integ_ = integ; }
  double integrator() const { return integ_; }
  void set_gains(double kp, double ki) { kp_ = kp; ki_ = ki; }
  double kp() const { return kp_; }
  double ki() const { return ki_; }

 private:
  double kp_{0}, ki_{0}, ts_{1e-3}, kaw_{0};
  double integ_{0};
};

}  // namespace actlab
