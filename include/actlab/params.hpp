#pragma once
// Physical and drive parameters for a humanoid-class rotary joint:
// frameless PMSM -> 50:1 strain-wave gearbox (compliant) -> link with gravity load.
// All angles are mechanical radians unless named *_e (electrical).

namespace actlab {

struct MotorParams {
  double R = 0.15;          // phase resistance [ohm]
  double L = 1.0e-4;        // phase inductance, Ld = Lq (surface PMSM) [H]
  double pole_pairs = 14.0;
  double lambda = 0.1 / (1.5 * 14.0);  // flux linkage [Wb], chosen so Kt = 0.1 Nm/A
  double Jm = 6.0e-5;       // rotor inertia [kg m^2]
  double bm = 2.0e-5;       // rotor viscous friction [Nm s/rad]
  double tau_c_m = 0.01;    // rotor Coulomb friction [Nm]

  double kt() const { return 1.5 * pole_pairs * lambda; }
};

struct JointParams {
  double N = 50.0;          // gear ratio
  double k = 4000.0;        // gearbox torsional stiffness, output side [Nm/rad]
  double d = 2.0;           // gearbox damping, output side [Nm s/rad] (zeta ~ 0.05)
  double Jl = 0.6;          // link inertia about the joint [kg m^2]
  double bl = 0.05;         // output bearing viscous friction [Nm s/rad]
  double tau_c_l = 0.5;     // output Coulomb friction [Nm]
  double mgl = 11.8;        // gravity moment arm m*g*l [Nm]
};

struct DriveParams {
  double Vdc = 48.0;        // bus voltage [V]
  double i_max = 30.0;      // q-axis current limit [A]
  int adc_bits = 12;
  double adc_range = 40.0;  // +/- full scale [A]
  double i_noise = 0.05;    // current sensor noise std [A]
  int enc_motor_bits = 14;  // motor-side incremental encoder, counts per motor rev
  int enc_link_bits = 17;   // output absolute encoder
};

struct Rates {
  double f_sim = 40000.0;   // plant integration rate [Hz]
  double f_cur = 20000.0;   // current loop = PWM rate [Hz]
  double f_vel = 5000.0;    // velocity loop [Hz]
  double f_pos = 1000.0;    // position loop [Hz]
};

struct PlantParams {
  MotorParams motor;
  JointParams joint;
  DriveParams drive;
};

}  // namespace actlab
