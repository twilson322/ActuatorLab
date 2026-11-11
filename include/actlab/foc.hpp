#pragma once
#include <algorithm>
#include <array>
#include <cmath>

// Field-oriented control transforms. Amplitude-invariant convention:
// a balanced phase current of amplitude I maps to |i_alpha_beta| = I.
namespace actlab::foc {

struct AlphaBeta { double alpha{0}, beta{0}; };
struct DQ { double d{0}, q{0}; };
struct ABC { double a{0}, b{0}, c{0}; };

inline constexpr double kSqrt3 = 1.7320508075688772;

// Two measured phase currents are enough because ia + ib + ic = 0 (wye, no neutral).
inline AlphaBeta clarke(double ia, double ib) {
  return {ia, (ia + 2.0 * ib) / kSqrt3};
}

inline ABC inverse_clarke(const AlphaBeta& ab) {
  const double a = ab.alpha;
  const double b = -0.5 * ab.alpha + 0.5 * kSqrt3 * ab.beta;
  return {a, b, -a - b};
}

inline DQ park(const AlphaBeta& ab, double theta_e) {
  const double c = std::cos(theta_e), s = std::sin(theta_e);
  return {ab.alpha * c + ab.beta * s, -ab.alpha * s + ab.beta * c};
}

inline AlphaBeta inverse_park(const DQ& dq, double theta_e) {
  const double c = std::cos(theta_e), s = std::sin(theta_e);
  return {dq.d * c - dq.q * s, dq.d * s + dq.q * c};
}

// Space-vector PWM via min/max zero-sequence injection. Mathematically equivalent to
// the 6-sector formulation and avoids sector bookkeeping. Returns duty cycles in [0, 1].
// Linear range: |v_alpha_beta| <= Vdc / sqrt(3). Beyond that the duties clip.
inline std::array<double, 3> svpwm(const AlphaBeta& v, double vdc) {
  const ABC ph = inverse_clarke(v);
  const double vmax = std::max({ph.a, ph.b, ph.c});
  const double vmin = std::min({ph.a, ph.b, ph.c});
  const double v0 = -0.5 * (vmax + vmin);
  auto duty = [&](double vx) { return std::clamp(0.5 + (vx + v0) / vdc, 0.0, 1.0); };
  return {duty(ph.a), duty(ph.b), duty(ph.c)};
}

// Average phase-to-neutral voltages the inverter actually produces for given duties.
// The motor neutral floats, so the common-mode part of the leg voltages disappears.
inline AlphaBeta inverter_average(const std::array<double, 3>& duty, double vdc) {
  const double va = vdc * duty[0], vb = vdc * duty[1], vc = vdc * duty[2];
  const double vn = (va + vb + vc) / 3.0;
  return clarke(va - vn, vb - vn);
}

// Limit a dq voltage vector to the inscribed circle of the SVPWM hexagon.
inline DQ limit_circle(const DQ& v, double vmax) {
  const double mag = std::hypot(v.d, v.q);
  if (mag <= vmax || mag < 1e-12) return v;
  const double s = vmax / mag;
  return {v.d * s, v.q * s};
}

}  // namespace actlab::foc
