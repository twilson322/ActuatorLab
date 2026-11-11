#include "actlab/plant.hpp"

#include <cmath>

namespace actlab {

namespace {
// Smooth Coulomb friction. A hard sign() makes RK4 chatter around zero speed. The knee
// width matters for stability of the explicit integrator: the local stiffness is
// tau_c / (eps * J). Motor side uses 0.05 rad/s (motor runs at hundreds of rad/s),
// link side 5 mrad/s.
inline double coulomb(double tau_c, double w, double eps) { return tau_c * std::tanh(w / eps); }

inline double quantize(double x, double lsb) { return std::floor(x / lsb + 0.5) * lsb; }

PlantState axpy(const PlantState& x, const PlantState& k, double h) {
  return {x.id + h * k.id, x.iq + h * k.iq, x.th_m + h * k.th_m,
          x.w_m + h * k.w_m, x.th_l + h * k.th_l, x.w_l + h * k.w_l};
}
}  // namespace

Plant::Plant(const PlantParams& p, std::uint64_t seed) : p_(p), rng_(seed) {}

double Plant::spring_torque() const {
  const auto& j = p_.joint;
  return j.k * (x_.th_m / j.N - x_.th_l) + j.d * (x_.w_m / j.N - x_.w_l);
}

PlantState Plant::deriv(const PlantState& x) const {
  const auto& m = p_.motor;
  const auto& j = p_.joint;
  const double w_e = m.pole_pairs * x.w_m;
  const foc::DQ v = foc::park(v_ab_, m.pole_pairs * x.th_m);

  PlantState dx;
  dx.id = (v.d - m.R * x.id + w_e * m.L * x.iq) / m.L;
  dx.iq = (v.q - m.R * x.iq - w_e * m.L * x.id - w_e * m.lambda) / m.L;

  const double tau_m = m.kt() * x.iq;
  const double tau_s = j.k * (x.th_m / j.N - x.th_l) + j.d * (x.w_m / j.N - x.w_l);
  dx.th_m = x.w_m;
  // Coulomb friction is smoothed with tanh over 0.05 rad/s at the motor and 0.005 rad/s at the
  // link, so the right-hand side stays smooth for RK4.
  dx.w_m = (tau_m - m.bm * x.w_m - coulomb(m.tau_c_m, x.w_m, 0.05) - tau_s / j.N) / m.Jm;
  dx.th_l = x.w_l;
  dx.w_l = (tau_s - j.bl * x.w_l - coulomb(j.tau_c_l, x.w_l, 0.005) - j.mgl * std::sin(x.th_l) +
            tau_ext_) / j.Jl;
  return dx;
}

void Plant::step(double dt) {
  const PlantState k1 = deriv(x_);
  const PlantState k2 = deriv(axpy(x_, k1, 0.5 * dt));
  const PlantState k3 = deriv(axpy(x_, k2, 0.5 * dt));
  const PlantState k4 = deriv(axpy(x_, k3, dt));
  x_.id += dt / 6.0 * (k1.id + 2 * k2.id + 2 * k3.id + k4.id);
  x_.iq += dt / 6.0 * (k1.iq + 2 * k2.iq + 2 * k3.iq + k4.iq);
  x_.th_m += dt / 6.0 * (k1.th_m + 2 * k2.th_m + 2 * k3.th_m + k4.th_m);
  x_.w_m += dt / 6.0 * (k1.w_m + 2 * k2.w_m + 2 * k3.w_m + k4.w_m);
  x_.th_l += dt / 6.0 * (k1.th_l + 2 * k2.th_l + 2 * k3.th_l + k4.th_l);
  x_.w_l += dt / 6.0 * (k1.w_l + 2 * k2.w_l + 2 * k3.w_l + k4.w_l);
}

Sensors Plant::sample() {
  const auto& d = p_.drive;
  const double theta_e = p_.motor.pole_pairs * x_.th_m;
  const foc::ABC i = foc::inverse_clarke(foc::inverse_park({x_.id, x_.iq}, theta_e));
  const double adc_lsb = 2.0 * d.adc_range / std::pow(2.0, d.adc_bits);
  Sensors s;
  s.ia = quantize(i.a + d.i_noise * noise_(rng_), adc_lsb);
  s.ib = quantize(i.b + d.i_noise * noise_(rng_), adc_lsb);
  s.th_m = quantize(x_.th_m, 2.0 * M_PI / std::pow(2.0, d.enc_motor_bits));
  s.th_l = quantize(x_.th_l, 2.0 * M_PI / std::pow(2.0, d.enc_link_bits));
  return s;
}

double Plant::energy() const {
  const auto& m = p_.motor;
  const auto& j = p_.joint;
  const double defl = x_.th_m / j.N - x_.th_l;
  return 0.5 * m.Jm * x_.w_m * x_.w_m + 0.5 * j.Jl * x_.w_l * x_.w_l + 0.5 * j.k * defl * defl +
         j.mgl * (1.0 - std::cos(x_.th_l)) + 0.75 * m.L * (x_.id * x_.id + x_.iq * x_.iq);
}

}  // namespace actlab
