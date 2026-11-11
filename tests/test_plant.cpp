#include <gtest/gtest.h>

#include <cmath>

#include "actlab/plant.hpp"

using namespace actlab;

// With every loss term removed the model conserves energy exactly, so whatever drift remains
// is RK4 truncation error. RK4 dissipates an oscillation of frequency w by ~(w h)^6 per step,
// i.e. drift per unit time scales with h^5: halving h must cut it by ~32x.
namespace {
double lossless_drift(double h) {
  PlantParams p;
  p.motor.R = 0.0;
  p.motor.bm = 0.0;
  p.motor.tau_c_m = 0.0;
  p.joint.d = 0.0;
  p.joint.bl = 0.0;
  p.joint.tau_c_l = 0.0;
  Plant plant(p);
  PlantState x;
  x.th_l = 0.02;  // released with the gearbox wound up by 20 mrad
  // (An earlier version also started with iq = 2 A. With R = 0 nothing limits the rotor and it
  // spun to 2200 rad/s, where w_e*h ~ 0.8 and RK4 is no longer in its asymptotic regime.)
  plant.set_state(x);
  plant.set_phase_voltage({0.0, 0.0});
  const double E0 = plant.energy();
  const int n = static_cast<int>(std::lround(1.0 / h));
  for (int k = 0; k < n; ++k) plant.step(h);  // 1 s
  return std::abs(plant.energy() / E0 - 1.0);
}
}  // namespace

TEST(Plant, EnergyConservedWhenLossless) {
  const double d1 = lossless_drift(25e-6), d2 = lossless_drift(12.5e-6);
  EXPECT_LT(d1, 1e-7);       // at the production step size
  EXPECT_LT(d2, d1 / 16.0);  // fourth-order integrator: drift collapses with h
}

TEST(Plant, StaticEquilibriumHolds) {
  PlantParams p;
  Plant plant(p);
  const double th0 = 0.7, tau_s = p.joint.mgl * std::sin(th0);
  PlantState x;
  x.th_l = th0;
  x.th_m = p.joint.N * (th0 + tau_s / p.joint.k);
  x.iq = tau_s / (p.joint.N * p.motor.kt());
  plant.set_state(x);
  // Hold the current with exactly R*iq in the q axis, expressed in the stationary frame.
  const double th_e = p.motor.pole_pairs * x.th_m;
  plant.set_phase_voltage(foc::inverse_park({0.0, p.motor.R * x.iq}, th_e));
  for (int k = 0; k < 4000; ++k) plant.step(25e-6);
  EXPECT_NEAR(plant.state().th_l, th0, 1e-6);
  EXPECT_NEAR(plant.state().w_l, 0.0, 1e-4);
}

TEST(Plant, LockedMotorModeFrequency) {
  // With the rotor clamped (huge Jm) the link rings at sqrt((k + mgl)/Jl) around th = 0.
  PlantParams p;
  p.motor.Jm = 1e6;
  p.joint.d = p.joint.bl = p.joint.tau_c_l = 0.0;
  Plant plant(p);
  PlantState x;
  x.th_l = 1e-3;
  plant.set_state(x);
  double t_prev = -1, period = 0;
  double prev = plant.state().th_l;
  int crossings = 0;
  for (int k = 0; k < 200000 && crossings < 3; ++k) {
    plant.step(25e-6);
    const double cur = plant.state().th_l;
    if (prev < 0 && cur >= 0) {
      const double t = k * 25e-6;
      if (t_prev > 0) period = t - t_prev;
      t_prev = t;
      ++crossings;
    }
    prev = cur;
  }
  const double f_expected = std::sqrt((p.joint.k + p.joint.mgl) / p.joint.Jl) / (2 * M_PI);
  EXPECT_NEAR(1.0 / period, f_expected, 0.05);
}
