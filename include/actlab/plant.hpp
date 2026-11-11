#pragma once
#include <cstdint>
#include <random>

#include "actlab/foc.hpp"
#include "actlab/params.hpp"

namespace actlab {

// Continuous state of the electromechanical joint.
// Currents are in the rotor (dq) frame. th_m/w_m are motor-side mechanical angle/speed,
// th_l/w_l are the link angle/speed (output side of the compliant gearbox).
struct PlantState {
  double id{0}, iq{0};
  double th_m{0}, w_m{0};
  double th_l{0}, w_l{0};
};

// What the drive electronics can actually see: two phase currents through a noisy,
// quantized ADC and two quantized encoders.
struct Sensors {
  double ia{0}, ib{0};
  double th_m{0};  // motor encoder, mechanical rad
  double th_l{0};  // link encoder, rad
};

class Plant {
 public:
  explicit Plant(const PlantParams& p, std::uint64_t seed = 1);

  // The inverter output is held constant in the stationary frame between PWM updates
  // (zero-order hold), which is what a real average-model inverter does.
  void set_phase_voltage(const foc::AlphaBeta& v) { v_ab_ = v; }
  void set_external_torque(double tau) { tau_ext_ = tau; }
  void step(double dt);  // one RK4 step
  Sensors sample();

  const PlantState& state() const { return x_; }
  void set_state(const PlantState& x) { x_ = x; }
  PlantParams& params() { return p_; }
  const PlantParams& params() const { return p_; }

  double motor_torque() const { return p_.motor.kt() * x_.iq; }
  double spring_torque() const;
  double copper_power() const { return 1.5 * p_.motor.R * (x_.id * x_.id + x_.iq * x_.iq); }
  // Total stored energy (kinetic + spring + gravity + magnetic). Used by tests.
  double energy() const;

 private:
  PlantState deriv(const PlantState& x) const;

  PlantParams p_;
  PlantState x_{};
  foc::AlphaBeta v_ab_{};
  double tau_ext_{0.0};
  std::mt19937_64 rng_;
  std::normal_distribution<double> noise_{0.0, 1.0};
};

}  // namespace actlab
