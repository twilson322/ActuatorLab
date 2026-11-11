#pragma once
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "actlab/controllers.hpp"
#include "actlab/drive.hpp"
#include "actlab/plant.hpp"

namespace actlab {

// Piecewise quintic point-to-point moves with holds in between (zero vel/acc at both ends).
class Trajectory {
 public:
  struct Segment { double t0, t1, p0, p1; };
  explicit Trajectory(double p_init = 0.0) : p_init_(p_init) {}
  Trajectory& move(double t0, double t1, double p1);
  Reference eval(double t) const;

 private:
  double p_init_;
  std::vector<Segment> segs_;
};

struct Window { double t0{0}, t1{0}; };

struct SimConfig {
  Rates rates{};
  DriveConfig drive{};
  double duration{3.0};
  Trajectory traj{0.0};
  double theta0{0.0};  // initial link angle (plant starts at static equilibrium)
  // External link torque profile and a torque injected at the controller output
  // (used for chirp identification and loop-gain measurement).
  std::function<double(double)> tau_ext;
  std::function<double(double)> tau_inject;
  std::vector<std::pair<double, std::function<void(Plant&)>>> events;  // e.g. payload change
  // Metric windows
  Window motion{0.0, 0.0};
  Window hold{0.0, 0.0};
  double dist_time{-1.0};
  double dist_window{1.0};  // recovery is measured inside [dist_time, dist_time + window]
  double dist_settle_tol{2e-3};
  // Logging
  std::string log_path;  // empty = no log
  int log_decim{1};      // relative to f_vel
};

struct Metrics {
  double rms_err_motion{0}, max_err_motion{0}, ss_err{0};
  double dist_peak{0}, dist_recovery{0}, dist_iae{0};  // iae = integral |e| dt over the window
  double peak_iq{0}, copper_J{0}, chatter{0};
  double max_abs_err{0};
  bool stable{true};
  double diverge_time{-1.0};  // first time |e| > 0.5 rad (run stops there), -1 if never
};

Metrics run_sim(Plant& plant, JointController& ctrl, const SimConfig& cfg);

}  // namespace actlab
