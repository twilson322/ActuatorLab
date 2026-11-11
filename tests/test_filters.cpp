#include <gtest/gtest.h>

#include <cmath>

#include "actlab/filters.hpp"

using namespace actlab;

namespace {
// Steady-state gain of a biquad at frequency f: RMS ratio over the second half of a 4 s run.
// (Peak picking fails near Nyquist because the samples rarely land on the crest.)
double measured_gain(Biquad<double> bq, double f, double fs) {
  double sy = 0.0, sx = 0.0;
  const int n = static_cast<int>(fs * 4.0);
  for (int k = 0; k < n; ++k) {
    const double x = std::sin(2.0 * kPi * f * k / fs);
    const double y = bq.step(x);
    if (k > n / 2) { sy += y * y; sx += x * x; }
  }
  return std::sqrt(sy / sx);
}
}  // namespace

TEST(Filters, NotchDepthAndPassband) {
  const double fs = 1000.0, f0 = 13.0;
  const auto n = Biquad<double>::notch(f0, 0.05, 0.7, fs);
  EXPECT_NEAR(measured_gain(n, f0, fs), 0.05 / 0.7, 0.01);  // depth = zn/zd at the centre
  EXPECT_NEAR(measured_gain(n, 0.5, fs), 1.0, 0.02);
  EXPECT_NEAR(measured_gain(n, 200.0, fs), 1.0, 0.02);
}

TEST(Filters, LowpassDcGainAndRolloff) {
  const double fs = 5000.0;
  const auto lp = Biquad<double>::lowpass(100.0, 0.707, fs);
  EXPECT_NEAR(measured_gain(lp, 2.0, fs), 1.0, 0.01);
  EXPECT_NEAR(measured_gain(lp, 100.0, fs), 0.707, 0.02);
  EXPECT_LT(measured_gain(lp, 1000.0, fs), 0.02);
}

TEST(Filters, BiquadResetIsSteadyState) {
  auto lp = Biquad<double>::lowpass(50.0, 0.7, 1000.0);
  lp.reset(3.0);
  for (int i = 0; i < 10; ++i) EXPECT_NEAR(lp.step(3.0), 3.0, 1e-12);
}

TEST(Filters, PiAntiWindupBoundsIntegrator) {
  PI pi(1.0, 50.0, 1e-3, 50.0);
  const double umax = 1.0;
  for (int i = 0; i < 5000; ++i) {  // large persistent error, actuator pinned at the limit
    const double e = 10.0, u = pi.output(e);
    pi.commit(e, u, std::clamp(u, -umax, umax));
  }
  // Back-calculation equilibrium: ki*e + kaw*(umax - kp*e - I) = 0  ->  I = 1 here.
  // Without anti-windup the integrator would have reached ki*e*5 s = 2500.
  EXPECT_NEAR(pi.integrator(), 1.0, 1e-6);
  // Error reverses: output must leave saturation within a few samples, not seconds.
  int k = 0;
  for (; k < 1000; ++k) {
    const double e = -0.5, u = pi.output(e);
    pi.commit(e, u, std::clamp(u, -umax, umax));
    if (u < umax) break;
  }
  EXPECT_LT(k, 5);
}
