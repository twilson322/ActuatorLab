#include "actlab/drive.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

#include "actlab/sim.hpp"
#include "actlab/telemetry.hpp"

namespace actlab {

// ============================================================================ Drive
Drive::Drive(const MotorParams& m, const DriveParams& dp, const DriveConfig& cfg, double f_cur)
    : m_(m), dp_(dp), cfg_(cfg), ts_(1.0 / f_cur), obs_(cfg.obs_bw_hz, f_cur) {
  // Pole-zero cancellation: PI zero at R/L cancels the electrical pole, leaving a
  // first-order closed loop with bandwidth wc (before delay effects).
  const double wc = 2.0 * kPi * cfg_.cur_bw_hz;
  const double kp = m_.L * wc, ki = m_.R * wc;
  pi_d_ = PI(kp, ki, ts_, ki / kp);
  pi_q_ = PI(kp, ki, ts_, ki / kp);
}

void Drive::reset(const Sensors& s) {
  obs_.reset(s.th_m);
  th_m_ = s.th_m;
}

void Drive::sense(const Sensors& s) {
  th_m_ = s.th_m;
  obs_.step(s.th_m);
  const foc::DQ i = foc::park(foc::clarke(s.ia, s.ib), m_.pole_pairs * s.th_m);
  id_ = i.d;
  iq_ = i.q;
}

foc::AlphaBeta Drive::control(double iq_ref) {
  const double w_e = m_.pole_pairs * obs_.omega();
  const double ed = 0.0 - id_, eq = iq_ref - iq_;
  // Cross-coupling and back-EMF decoupling feedforward.
  const double vd_u = pi_d_.output(ed) - w_e * m_.L * iq_;
  const double vq_u = pi_q_.output(eq) + w_e * (m_.L * id_ + m_.lambda);
  const double vmax = cfg_.vlimit_frac * dp_.Vdc / foc::kSqrt3;
  const foc::DQ v = foc::limit_circle({vd_u, vq_u}, vmax);
  sat_ = std::hypot(vd_u, vq_u) > vmax;
  pi_d_.commit(ed, vd_u, v.d);
  pi_q_.commit(eq, vq_u, v.q);
  vq_ = v.q;
  // The voltage is applied next period and held for a full period, so on average the rotor
  // will have advanced by 1.5 periods when it acts.
  const double th_apply = m_.pole_pairs * th_m_ + (cfg_.delay_comp ? 1.5 * ts_ * w_e : 0.0);
  return foc::inverter_average(foc::svpwm(foc::inverse_park(v, th_apply), dp_.Vdc), dp_.Vdc);
}

// ============================================================================ Trajectory
Trajectory& Trajectory::move(double t0, double t1, double p1) {
  const double p0 = segs_.empty() ? p_init_ : segs_.back().p1;
  segs_.push_back({t0, t1, p0, p1});
  return *this;
}

Reference Trajectory::eval(double t) const {
  double p = p_init_;
  for (const auto& s : segs_) {
    if (t < s.t0) return {p, 0.0, 0.0};
    if (t <= s.t1) {
      const double T = s.t1 - s.t0, tau = (t - s.t0) / T, d = s.p1 - s.p0;
      const double t3 = tau * tau * tau, t4 = t3 * tau, t5 = t4 * tau;
      return {s.p0 + d * (10 * t3 - 15 * t4 + 6 * t5),
              d * (30 * tau * tau - 60 * t3 + 30 * t4) / T,
              d * (60 * tau - 180 * tau * tau + 120 * t3) / (T * T)};
    }
    p = s.p1;
  }
  return {p, 0.0, 0.0};
}

// ============================================================================ run_sim
Metrics run_sim(Plant& plant, JointController& ctrl, const SimConfig& cfg) {
  const auto& P = plant.params();
  const double dt = 1.0 / cfg.rates.f_sim;
  const int sim_per_cur = static_cast<int>(std::lround(cfg.rates.f_sim / cfg.rates.f_cur));
  const int cur_per_vel = static_cast<int>(std::lround(cfg.rates.f_cur / cfg.rates.f_vel));
  const int vel_per_pos = static_cast<int>(std::lround(cfg.rates.f_vel / cfg.rates.f_pos));
  const double tau_max = P.motor.kt() * P.drive.i_max;

  // Start at static equilibrium: the gearbox already carries the gravity load.
  {
    const double th0 = cfg.theta0;
    const double tau_s = P.joint.mgl * std::sin(th0);
    PlantState x0;
    x0.th_l = th0;
    x0.th_m = P.joint.N * (th0 + tau_s / P.joint.k);
    x0.iq = tau_s / (P.joint.N * P.motor.kt());
    plant.set_state(x0);
  }

  Drive drive(P.motor, P.drive, cfg.drive, cfg.rates.f_cur);
  Sensors s = plant.sample();
  drive.reset(s);
  drive.sense(s);
  double iq_ref = plant.state().iq;
  foc::AlphaBeta pending = foc::inverse_park({0.0, P.motor.R * iq_ref}, P.motor.pole_pairs * plant.state().th_m);
  ctrl.reset(cfg.traj.eval(0.0), {s.th_m, 0.0, s.th_l});

  std::unique_ptr<CsvTelemetry> log;
  if (!cfg.log_path.empty()) {
    log = std::make_unique<CsvTelemetry>(
        cfg.log_path, std::vector<std::string>{"t", "ref", "th_l", "th_m_out", "err", "tau_cmd",
                                               "tau_inj", "iq", "w_m", "w_m_hat", "tau_ext",
                                               "d0", "d1", "d2", "d3", "d4"});
  }

  Metrics M;
  double sum_sq_motion = 0.0, sum_ss = 0.0, sum_chat = 0.0, tau_prev = 0.0;
  long n_motion = 0, n_ss = 0, n_chat = 0, cur_tick = 0, vel_tick = 0;
  double last_exceed = cfg.dist_time;
  std::size_t next_event = 0;
  const long n_steps = static_cast<long>(std::lround(cfg.duration * cfg.rates.f_sim));

  for (long k = 0; k < n_steps; ++k) {
    const double t = static_cast<double>(k) * dt;
    while (next_event < cfg.events.size() && t >= cfg.events[next_event].first) {
      cfg.events[next_event].second(plant);
      ++next_event;
    }
    const double tau_ext = cfg.tau_ext ? cfg.tau_ext(t) : 0.0;
    plant.set_external_torque(tau_ext);

    if (k % sim_per_cur == 0) {
      plant.set_phase_voltage(pending);
      s = plant.sample();
      drive.sense(s);

      if (cur_tick % cur_per_vel == 0) {
        const bool pos_tick = (vel_tick % vel_per_pos) == 0;
        const Reference r = cfg.traj.eval(t);
        const Feedback fb{drive.th_m(), drive.w_m(), s.th_l};
        const double inj = cfg.tau_inject ? cfg.tau_inject(t) : 0.0;
        const double tau_ctrl = ctrl.update(r, fb, pos_tick);
        const double tau_u = tau_ctrl + inj;
        const double tau_s = std::clamp(tau_u, -tau_max, tau_max);
        ctrl.commit(tau_u, tau_s);
        iq_ref = tau_s / P.motor.kt();

        // ------------------------------------------------ metrics (true link angle)
        const PlantState& x = plant.state();
        const double e = x.th_l - r.pos;
        if (!std::isfinite(e) || std::abs(e) > 0.5 || std::abs(x.w_m) > 2000.0) {
          M.stable = false;
          M.diverge_time = t;
          break;
        }
        M.max_abs_err = std::max(M.max_abs_err, std::abs(e));
        if (t >= cfg.motion.t0 && t <= cfg.motion.t1) {
          sum_sq_motion += e * e;
          ++n_motion;
          M.max_err_motion = std::max(M.max_err_motion, std::abs(e));
        }
        if (t >= cfg.hold.t0 && t <= cfg.hold.t1) {
          sum_ss += std::abs(e);
          ++n_ss;
        }
        if (cfg.dist_time >= 0 && t >= cfg.dist_time && t <= cfg.dist_time + cfg.dist_window) {
          M.dist_peak = std::max(M.dist_peak, std::abs(e));
          M.dist_iae += std::abs(e) / cfg.rates.f_vel;
          if (std::abs(e) > cfg.dist_settle_tol) last_exceed = t;
        }
        if (vel_tick > 0) {
          sum_chat += (tau_s - tau_prev) * (tau_s - tau_prev);
          ++n_chat;
        }
        tau_prev = tau_s;
        M.peak_iq = std::max(M.peak_iq, std::abs(x.iq));

        if (log && vel_tick % cfg.log_decim == 0) {
          LogRow row;
          const double vals[LogRow::kCols] = {t, r.pos, x.th_l, x.th_m / P.joint.N, e, tau_s, inj,
                                              x.iq, x.w_m, drive.w_m(), tau_ext, ctrl.diag(0),
                                              ctrl.diag(1), ctrl.diag(2), ctrl.diag(3), ctrl.diag(4)};
          std::copy(std::begin(vals), std::end(vals), row.v);
          log->log(row);
        }
        ++vel_tick;
      }
      pending = drive.control(iq_ref);
      ++cur_tick;
    }
    plant.step(dt);
    M.copper_J += plant.copper_power() * dt;
  }

  M.rms_err_motion = n_motion ? std::sqrt(sum_sq_motion / static_cast<double>(n_motion)) : 0.0;
  M.ss_err = n_ss ? sum_ss / static_cast<double>(n_ss) : 0.0;
  M.chatter = n_chat ? std::sqrt(sum_chat / static_cast<double>(n_chat)) : 0.0;
  M.dist_recovery = cfg.dist_time >= 0 ? last_exceed - cfg.dist_time : 0.0;
  return M;
}

}  // namespace actlab
