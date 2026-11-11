#pragma once
#include <cmath>

#include "actlab/filters.hpp"
#include "actlab/rls.hpp"

namespace actlab {

// Second-order state-variable filter H(s) = wf^2 / (s^2 + 2 z wf s + wf^2).
// Gives H*u, H*s*u and H*s^2*u from one input, which is how you get a usable
// acceleration regressor without differentiating a quantized encoder.
//
// The acceleration output is the mean of x2' over the last sample, (x2[k] - x2[k-1]) / T. It
// used to be x2' evaluated at the end of the sample, which contains the feedthrough term
// wf^2 u[k] of the newest input and so leads the other outputs by about half a sample. On a
// rigid joint with an exact model that biased the inertia estimate by +6.5 % and produced a
// viscous coefficient of 3.1 N m s/rad on a joint with none (tests/test_rls.cpp).
class SVF2 {
 public:
  SVF2(double f_hz, double zeta, double fs, int substeps = 4)
      : w_(2.0 * kPi * f_hz), z_(zeta), h_(1.0 / (fs * substeps)), n_(substeps) {}
  void reset(double u0) { x1_ = u0; x2_ = 0.0; acc_ = 0.0; }
  void step(double u) {
    const double x2_before = x2_;
    for (int i = 0; i < n_; ++i) {
      const double a = -w_ * w_ * x1_ - 2.0 * z_ * w_ * x2_ + w_ * w_ * u;
      x1_ += h_ * x2_;
      x2_ += h_ * a;
    }
    acc_ = (x2_ - x2_before) / (h_ * n_);
  }
  double y() const { return x1_; }
  double dy() const { return x2_; }
  double ddy() const { return acc_; }

 private:
  double w_, z_, h_;
  int n_;
  double x1_{0}, x2_{0}, acc_{0};
};

struct OutputRlsConfig {
  double filter_bandwidth_hz = 10.0;  // well below the 13 Hz mode, well above the motion
  double filter_damping = 0.8;
  double initial_covariance = 5.0;
  double forgetting = 0.998;           // per position-loop sample
  double covariance_trace_max = 1e4;
  double coulomb_smoothing_rad_s = 0.02;  // tanh width standing in for sign(velocity)
  // Only learn while moving. At rest the Coulomb term can sit anywhere in [-tc, tc] and would be
  // absorbed into the gravity estimate.
  double min_learning_speed_rad_s = 0.05;
};

// Online identification of the output-side rigid-joint model
//   N tau_m = J th'' + b th' + tc tanh(th' / eps) + mgl sin(th),   theta = [J, b, tc, mgl].
// Both sides of the regression go through the same SVF, so the filter's phase lag cancels.
// Shared by the adaptive cascade and the adaptive LQG, which apply different parts of theta.
class OutputSideRls {
 public:
  using Vec = RLS<4>::Vec;
  OutputSideRls(const Vec& theta0, double sample_rate_hz, const OutputRlsConfig& c = {})
      : c_(c),
        th_(c.filter_bandwidth_hz, c.filter_damping, sample_rate_hz),
        tau_(c.filter_bandwidth_hz, c.filter_damping, sample_rate_hz),
        friction_(c.filter_bandwidth_hz, c.filter_damping, sample_rate_hz),
        sine_(c.filter_bandwidth_hz, c.filter_damping, sample_rate_hz),
        rls_(theta0, c.initial_covariance, c.forgetting, c.covariance_trace_max) {}

  void reset(double th_link, double tau_output) {
    th_.reset(th_link);
    tau_.reset(tau_output);
    friction_.reset(0.0);
    sine_.reset(std::sin(th_link));
  }

  // One position-loop sample. Returns true if the estimate was updated.
  bool update(double th_link, double tau_output) {
    th_.step(th_link);
    tau_.step(tau_output);
    friction_.step(std::tanh(th_.dy() / c_.coulomb_smoothing_rad_s));
    sine_.step(std::sin(th_link));
    if (std::abs(th_.dy()) < c_.min_learning_speed_rad_s) return false;
    const Vec phi(th_.ddy(), th_.dy(), friction_.y(), sine_.y());
    rls_.update(phi, tau_.y());
    return true;
  }

  const Vec& theta() const { return rls_.theta(); }

 private:
  OutputRlsConfig c_;
  SVF2 th_, tau_, friction_, sine_;
  RLS<4> rls_;
};

}  // namespace actlab
