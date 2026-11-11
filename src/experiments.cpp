#include "actlab/experiments.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "actlab/lqg.hpp"
#include "actlab/telemetry.hpp"

namespace actlab::exp {

namespace {

std::map<std::string, double> load_kv(const std::string& path) {
  std::map<std::string, double> kv;
  std::ifstream f(path);
  std::string line;
  while (std::getline(f, line)) {
    const auto eq = line.find('=');
    if (eq == std::string::npos || line[0] == '#') continue;
    kv[line.substr(0, eq)] = std::stod(line.substr(eq + 1));
  }
  return kv;
}

// Logarithmic chirp, zero outside [t0, t0 + T].
std::function<double(double)> log_chirp(double amp, double f0, double f1, double t0, double T) {
  const double k = std::log(f1 / f0);
  return [=](double t) {
    if (t < t0 || t > t0 + T) return 0.0;
    const double tau = t - t0;
    const double phase = 2.0 * kPi * f0 * T / k * (std::exp(k * tau / T) - 1.0);
    // 0.2 s cosine fade in/out avoids broadband transients at the edges.
    const double w = std::min({1.0, tau / 0.2, (T - tau) / 0.2});
    return amp * w * std::sin(phase);
  };
}

void write_metrics_header(std::FILE* f) {
  std::fprintf(f,
               "controller,stable,rms_err_motion_mrad,max_err_motion_mrad,ss_err_mrad,"
               "dist_peak_mrad,dist_iae_mrad_s,dist_recovery_ms,peak_iq_A,copper_J,chatter_Nm\n");
}
void write_metrics_row(std::FILE* f, const std::string& name, const Metrics& m) {
  std::fprintf(f, "%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,%.3f,%.4f,%.5f\n", name.c_str(),
               m.stable ? 1 : 0, 1e3 * m.rms_err_motion, 1e3 * m.max_err_motion, 1e3 * m.ss_err,
               1e3 * m.dist_peak, 1e3 * m.dist_iae, 1e3 * m.dist_recovery, m.peak_iq, m.copper_J, m.chatter);
}

Metrics simulate(const PlantParams& pp, const ControllerConfig& cc, SimConfig sc,
                 std::uint64_t noise_seed = 1) {
  Plant plant(pp, noise_seed);
  const NominalModel nominal = NominalModel::from(PlantParams{});  // designer's belief
  auto ctrl = make_controller(cc, nominal, sc.rates.f_vel, sc.rates.f_pos);
  return run_sim(plant, *ctrl, sc);
}

}  // namespace

// ---------------------------------------------------------------------------- configs
namespace {
void apply_override(ControllerConfig& c, const std::string& key, double v) {
  if (key == "vel.bw_hz") c.vel.bw_hz = v;
  else if (key == "vel.pi_zero_ratio") c.vel.pi_zero_ratio = v;
  else if (key == "vel.dob") c.vel.dob = v != 0.0;
  else if (key == "vel.dob_hz") c.vel.dob_hz = v;
  else if (key == "vel.notch") c.vel.notch = v != 0.0;
  else if (key == "pos_kp_ratio") c.pos_kp_ratio = v;
  else if (key == "defl") c.deflection_comp = v != 0.0;
  else if (key == "link.enabled") c.link.enabled = v != 0.0;
  else if (key == "link.kp") c.link.kp = v;
  else if (key == "link.ki") c.link.ki = v;
  else if (key == "link.notch") c.link.notch = v != 0.0;
  else if (key == "link.notch_f0") c.link.notch_f0 = v;
  else if (key == "link.notch_zd") c.link.notch_zd = v;
  else if (key == "adrc.wc_hz") c.adrc_wc_hz = v;
  else if (key == "adrc.wo_hz") c.adrc_wo_hz = v;
  else if (key == "smc.lambda") c.smc_lambda = v;
  else if (key == "smc.k1") c.smc_k1 = v;
  else if (key == "smc.k2") c.smc_k2 = v;
  else if (key == "smc.phi") c.smc_phi = v;
  else if (key == "adaptive") c.adaptive = v != 0.0;
  else if (key == "lqg.w_qm") c.lqg.w_qm = v;
  else if (key == "lqg.w_vm") c.lqg.w_vm = v;
  else if (key == "lqg.w_ql") c.lqg.w_ql = v;
  else if (key == "lqg.w_vl") c.lqg.w_vl = v;
  else if (key == "lqg.w_int") c.lqg.w_int = v;
  else if (key == "lqg.r_u") c.lqg.r_u = v;
  else if (key == "lqg.q_tau_m") c.lqg.q_tau_m = v;
  else if (key == "lqg.q_tau_l") c.lqg.q_tau_l = v;
  else if (key == "lqg.q_dist") c.lqg.q_dist = v;
  else throw std::invalid_argument("unknown override key " + key);
}
}  // namespace

ControllerConfig named_config(const std::string& name, const Options& o) {
  const auto kv = load_kv(o.identified_path);
  // Without an identification file the notches fall back to the plant's design values. Say so:
  // a silently missing file would put every notch in the wrong place on a different plant.
  static std::atomic<bool> warned{false};  // named_config may run on the Monte Carlo workers
  if ((!kv.count("f_antires_hz") || !kv.count("f_res_hz")) && !warned.exchange(true))
    std::fprintf(stderr, "note: %s has no identified modes, using the design values 13 / 29 Hz\n", o.identified_path.c_str());
  const double f_ar = kv.count("f_antires_hz") ? kv.at("f_antires_hz") : 13.0;
  const double f_r = kv.count("f_res_hz") ? kv.at("f_res_hz") : 29.0;

  // Gains below are the grid-search winners from `actlab_sim tune` (results/tuning.csv).
  ControllerConfig c;
  c.vel.bw_hz = 40.0;
  c.pos_kp_ratio = 1.0;
  c.adrc_wc_hz = 10.0;
  c.adrc_wo_hz = 100.0;
  c.smc_lambda = 60.0;
  c.smc_k1 = 60.0;
  c.smc_k2 = 100.0;
  c.lqg.w_ql = 2e7;
  c.lqg.w_vl = 2e3;
  c.lqg.w_int = 1e7;
  c.lqg.q_tau_l = 0.2;
  c.lqg.q_dist = 4e3;
  c.vel.notch_f0 = f_r;
  c.link.notch_f0 = f_ar;
  LinkLoopConfig dual;
  dual.enabled = true;
  dual.kp = 1.0;
  dual.ki = 5.0;
  dual.notch = true;
  dual.notch_f0 = f_ar;

  if (name == "motor_only") {
    c.deflection_comp = false;
  } else if (name == "motor_defl") {
  } else if (name == "motor_defl_velnotch") {  // ablation: notch inside the collocated loop
    c.vel.notch = true;
  } else if (name == "motor_defl_dob") {  // ablation: disturbance observer on the speed loop
    c.vel.dob = true;
  } else if (name == "dual_pi_nonotch") {
    c.link = dual;
    c.link.notch = false;
  } else if (name == "dual_pi_notch") {
    c.link = dual;
  } else if (name == "dual_rls") {
    c.link = dual;
    c.adaptive = true;
  } else if (name == "lqg") {
    c.type = "lqg";
  } else if (name == "lqg_nodist") {
    c.type = "lqg";
    c.lqg.use_dist_ff = false;
  } else if (name == "lqg_rls") {
    c.type = "lqg";
    c.lqg.adaptive = true;
  } else if (name == "adrc_dual") {
    c.type = "adrc";
    c.link = dual;
  } else if (name == "stsmc_dual") {
    c.type = "stsmc";
    c.link = dual;
  } else {
    throw std::invalid_argument("unknown controller config " + name);
  }
  for (const auto& [k, v] : o.overrides) apply_override(c, k, v);
  return c;
}

std::vector<std::pair<std::string, ControllerConfig>> controller_suite(const Options& o) {
  std::vector<std::pair<std::string, ControllerConfig>> v;
  for (const char* n : {"motor_only", "motor_defl", "motor_defl_velnotch", "motor_defl_dob",
                        "dual_pi_notch", "adrc_dual", "stsmc_dual", "lqg", "lqg_nodist"})
    v.emplace_back(n, named_config(n, o));
  return v;
}

SimConfig tracking_scenario() {
  SimConfig sc;
  sc.duration = 2.6;
  sc.traj = Trajectory(0.0);
  sc.traj.move(0.2, 0.6, kPi / 2).move(1.0, 1.3, kPi / 4);
  sc.tau_ext = [](double t) { return (t >= 1.6 && t < 2.2) ? 25.0 : 0.0; };
  sc.dist_window = 0.6;
  sc.motion = {0.2, 1.3};
  sc.hold = {0.8, 1.0};
  sc.dist_time = 1.6;
  return sc;
}

// ---------------------------------------------------------------------------- experiments
namespace {

int exp_track(const Options& o) {
  std::FILE* f = std::fopen((o.out_dir + "/track_metrics.csv").c_str(), "w");
  write_metrics_header(f);
  for (const auto& [name, cc] : controller_suite(o)) {
    SimConfig sc = tracking_scenario();
    sc.log_path = o.out_dir + "/track_" + name + ".csv";
    sc.log_decim = 5;  // 1 kHz log
    const Metrics m = simulate(PlantParams{}, cc, sc);
    write_metrics_row(f, name, m);
    std::printf("%-20s stable=%d rms=%.3f mrad  max=%.3f mrad  dist_peak=%.2f mrad  "
                "recovery=%.0f ms  chatter=%.4f Nm\n",
                name.c_str(), m.stable, 1e3 * m.rms_err_motion, 1e3 * m.max_err_motion,
                1e3 * m.dist_peak, 1e3 * m.dist_recovery, m.chatter);
  }
  std::fclose(f);
  return 0;
}

// Open-loop plant FRF: a very soft position hold keeps the joint near 0.3 rad while a
// log chirp is added to the torque command. Python estimates tau -> w_m with Welch H1.
int exp_chirp_plant(const Options& o) {
  ControllerConfig cc;
  cc.vel.bw_hz = 3.0;
  SimConfig sc;
  sc.theta0 = 0.3;
  sc.traj = Trajectory(0.3);
  sc.duration = 26.0;
  sc.tau_inject = log_chirp(0.15, 2.0, 600.0, 1.0, 24.0);
  sc.log_path = o.out_dir + "/chirp_plant.csv";
  sc.log_decim = 1;  // 5 kHz
  const Metrics m = simulate(PlantParams{}, cc, sc);
  std::printf("chirp_plant done, stable=%d max|e|=%.2f mrad\n", m.stable, 1e3 * m.max_abs_err);
  return m.stable ? 0 : 1;
}

// Loop gain broken at the plant input: inject d, measure what the controller sends back.
int exp_loopgain(const Options& o) {
  for (const char* name : {"motor_defl", "motor_defl_velnotch", "dual_pi_notch", "adrc_dual", "stsmc_dual", "lqg"}) {
    SimConfig sc;
    sc.theta0 = 0.3;
    sc.traj = Trajectory(0.3);
    sc.duration = 32.0;
    sc.tau_inject = log_chirp(0.04, 1.0, 1000.0, 1.0, 30.0);
    sc.log_path = o.out_dir + "/loopgain_" + std::string(name) + ".csv";
    sc.log_decim = 1;
    const Metrics m = simulate(PlantParams{}, named_config(name, o), sc);
    std::printf("loopgain %-18s stable=%d\n", name, m.stable);
  }
  return 0;
}

// How far can the speed-loop bandwidth be pushed before the resonance takes over?
int exp_sweep(const Options& o) {
  std::FILE* f = std::fopen((o.out_dir + "/sweep.csv").c_str(), "w");
  std::fprintf(f, "notch,bw_hz,stable,rms_err_motion_mrad,dist_peak_mrad,chatter_Nm\n");
  for (int notch = 0; notch <= 1; ++notch) {
    for (double bw = 10.0; bw <= 120.0 + 1e-9; bw += 5.0) {
      ControllerConfig cc = named_config(notch ? "motor_defl_velnotch" : "motor_defl", o);
      cc.vel.bw_hz = bw;
      const Metrics m = simulate(PlantParams{}, cc, tracking_scenario());
      std::fprintf(f, "%d,%.1f,%d,%.4f,%.4f,%.5f\n", notch, bw, m.stable, 1e3 * m.rms_err_motion,
                   1e3 * m.dist_peak, m.chatter);
    }
  }
  std::fclose(f);
  std::printf("sweep written\n");
  return 0;
}

// Plant parameter Monte Carlo. Controllers keep the nominal model; the plant varies.
// Runs are independent, so a small thread pool pulls indices from an atomic counter and
// writes into pre-sized result slots (one writer per slot, no locking needed).
int exp_montecarlo(const Options& o) {
  struct Sample { PlantParams pp; };
  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  auto uni = [&](double a, double b) { return a + (b - a) * U(rng); };
  std::vector<Sample> samples(static_cast<std::size_t>(o.mc_samples));
  for (auto& s : samples) {
    PlantParams p;
    p.joint.Jl *= std::exp(uni(std::log(0.5), std::log(3.0)));
    p.joint.k *= uni(0.7, 1.3);
    p.joint.d *= uni(0.5, 2.0);
    const double fr = uni(0.5, 2.0);
    p.joint.tau_c_l *= fr;
    p.motor.tau_c_m *= fr;
    p.joint.bl *= uni(0.5, 2.0);
    p.joint.mgl *= uni(0.6, 1.6);
    p.motor.Jm *= uni(0.9, 1.1);
    s.pp = p;
  }
  const auto suite = controller_suite(o);
  const std::size_t n_jobs = samples.size() * suite.size();
  std::vector<Metrics> out(n_jobs);
  std::atomic<std::size_t> next{0};
  auto worker = [&] {
    for (std::size_t j = next.fetch_add(1); j < n_jobs; j = next.fetch_add(1)) {
      const auto& s = samples[j / suite.size()];
      const auto& cc = suite[j % suite.size()].second;
      out[j] = simulate(s.pp, cc, tracking_scenario(), 100 + j / suite.size());
    }
  };
  std::vector<std::thread> pool;
  for (int i = 0; i < std::max(1, o.threads); ++i) pool.emplace_back(worker);
  for (auto& t : pool) t.join();

  std::FILE* f = std::fopen((o.out_dir + "/montecarlo.csv").c_str(), "w");
  std::fprintf(f, "sample,controller,Jl,k,mgl,friction_scale,stable,rms_err_motion_mrad,"
                  "max_err_motion_mrad,dist_peak_mrad,dist_iae_mrad_s,dist_recovery_ms,chatter_Nm\n");
  for (std::size_t j = 0; j < n_jobs; ++j) {
    const auto& s = samples[j / suite.size()];
    const auto& m = out[j];
    std::fprintf(f, "%zu,%s,%.4f,%.1f,%.3f,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.2f,%.5f\n", j / suite.size(),
                 suite[j % suite.size()].first.c_str(), s.pp.joint.Jl, s.pp.joint.k,
                 s.pp.joint.mgl, s.pp.joint.tau_c_l / PlantParams{}.joint.tau_c_l, m.stable,
                 1e3 * m.rms_err_motion, 1e3 * m.max_err_motion, 1e3 * m.dist_peak,
                 1e3 * m.dist_iae, 1e3 * m.dist_recovery, m.chatter);
  }
  std::fclose(f);
  std::printf("montecarlo: %zu runs written\n", n_jobs);
  return 0;
}

// Payload change mid-run (e.g. the hand picks up a 4 kg object). The fixed-gain cascade
// keeps its nominal feedforward; the RLS variant re-identifies J and mgl and re-schedules.
int exp_adapt(const Options& o) {
  std::FILE* f = std::fopen((o.out_dir + "/adapt_metrics.csv").c_str(), "w");
  std::fprintf(f, "controller,window,rms_err_mrad,diverged_at_s\n");
  for (const char* name : {"dual_pi_notch", "dual_rls", "lqg", "lqg_rls"}) {
    SimConfig sc;
    sc.duration = 8.0;
    sc.traj = Trajectory(0.0);
    double t = 0.2;
    bool up = true;
    while (t + 0.5 < sc.duration) {
      sc.traj.move(t, t + 0.4, up ? 1.0 : 0.0);
      t += 0.5;
      up = !up;
    }
    sc.events.push_back({4.0, [](Plant& p) {
                           p.params().joint.Jl = 1.5;
                           p.params().joint.mgl = 24.0;
                         }});
    sc.log_path = o.out_dir + "/adapt_" + std::string(name) + ".csv";
    sc.log_decim = 5;
    // Score the second before and the last two seconds after the change separately.
    for (const auto& [label, w] : std::vector<std::pair<std::string, Window>>{
             {"before_3to4s", {3.0, 4.0}}, {"after_6to8s", {6.0, 8.0}}}) {
      SimConfig s2 = sc;
      s2.motion = w;
      s2.duration = w.t1;  // a later divergence must not hide the score of an earlier window
      if (label == "before_3to4s") s2.log_path.clear();
      const Metrics m = simulate(PlantParams{}, named_config(name, o), s2);
      if (m.stable) {
        std::fprintf(f, "%s,%s,%.4f,\n", name, label.c_str(), 1e3 * m.rms_err_motion);
        std::printf("adapt %-16s %-14s rms=%.3f mrad\n", name, label.c_str(), 1e3 * m.rms_err_motion);
      } else {
        std::fprintf(f, "%s,%s,,%.3f\n", name, label.c_str(), m.diverge_time);
        std::printf("adapt %-16s %-14s DIVERGED at t=%.3f s\n", name, label.c_str(), m.diverge_time);
      }
    }
  }
  std::fclose(f);
  return 0;
}

// Grid search per controller family on the nominal plant, same objective for everyone:
//   J = rms tracking error + 0.5 * peak disturbance deviation   [mrad]
//       + 2 * integral of |e| over the load window              [mrad s]
//       + a penalty once torque-command chatter exceeds 0.04 Nm (motor side, per 5 kHz step).
// The IAE term was added after the first tuning pass picked a link loop that rang for the
// whole load window: RMS and peak alone do not see sustained ringing.
// The winners are frozen into named_config(); the Monte Carlo then scores them on plants
// they were not tuned on.
int exp_tune(const Options& o) {
  std::FILE* f = std::fopen((o.out_dir + "/tuning.csv").c_str(), "w");
  std::fprintf(f, "family,params,stable,rms_mrad,dist_peak_mrad,dist_iae_mrad_s,recovery_ms,chatter_Nm,cost\n");
  auto cost = [](const Metrics& m) {
    if (!m.stable) return 1e9;
    return 1e3 * m.rms_err_motion + 0.5e3 * m.dist_peak + 2e3 * m.dist_iae +
           1e3 * std::max(0.0, m.chatter - 0.04);
  };
  struct Best { double J{1e18}; std::string params; };
  auto eval = [&](const std::string& family, const std::string& base,
                  const std::vector<std::pair<std::string, double>>& sets, Best& best) {
    Options oo = o;
    oo.overrides = sets;
    const Metrics m = simulate(PlantParams{}, named_config(base, oo), tracking_scenario());
    std::string ps;
    for (const auto& [k, v] : sets) ps += k + "=" + std::to_string(v).substr(0, 8) + " ";
    const double J = cost(m);
    std::fprintf(f, "%s,%s,%d,%.4f,%.4f,%.4f,%.1f,%.5f,%.4f\n", family.c_str(), ps.c_str(), m.stable,
                 1e3 * m.rms_err_motion, 1e3 * m.dist_peak, 1e3 * m.dist_iae, 1e3 * m.dist_recovery,
                 m.chatter, J);
    if (J < best.J) best = {J, ps};
  };
  Best b_motor, b_dual, b_dual_nn, b_adrc, b_smc, b_lqg;
  for (double bw : {10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 80.0})
    for (double kr : {0.1, 0.25, 0.5, 1.0})
      eval("motor_defl", "motor_defl", {{"vel.bw_hz", bw}, {"pos_kp_ratio", kr}}, b_motor);
  for (double kp : {0.0, 0.1, 0.2, 0.4, 0.6, 1.0})
    for (double ki : {2.0, 5.0, 10.0, 20.0, 40.0}) {
      eval("dual_pi_notch", "dual_pi_notch", {{"link.kp", kp}, {"link.ki", ki}}, b_dual);
      eval("dual_pi_nonotch", "dual_pi_nonotch", {{"link.kp", kp}, {"link.ki", ki}}, b_dual_nn);
    }
  for (double wc : {3.0, 6.0, 10.0, 15.0})
    for (double wo : {40.0, 60.0, 100.0, 150.0, 250.0})
      eval("adrc_dual", "adrc_dual", {{"adrc.wc_hz", wc}, {"adrc.wo_hz", wo}}, b_adrc);
  for (double l : {20.0, 60.0, 120.0, 200.0})
    for (double k1 : {15.0, 30.0, 60.0})
      for (double k2 : {50.0, 100.0, 300.0})
        eval("stsmc_dual", "stsmc_dual", {{"smc.lambda", l}, {"smc.k1", k1}, {"smc.k2", k2}}, b_smc);
  for (double wql : {4e6, 2e7, 1e8})
    for (double wvl : {2e3, 2e4})
      for (double wint : {4e7, 1e9})
        for (double qd : {4e3, 4e4})
          eval("lqg", "lqg", {{"lqg.w_ql", wql}, {"lqg.w_vl", wvl}, {"lqg.w_int", wint},
                              {"lqg.q_tau_l", 0.2}, {"lqg.q_dist", qd}}, b_lqg);
  std::fclose(f);
  std::printf("best motor_defl     J=%7.3f  %s\n", b_motor.J, b_motor.params.c_str());
  std::printf("best dual_pi_notch  J=%7.3f  %s\n", b_dual.J, b_dual.params.c_str());
  std::printf("best dual_pi_nonotch J=%7.3f %s\n", b_dual_nn.J, b_dual_nn.params.c_str());
  std::printf("best adrc_dual      J=%7.3f  %s\n", b_adrc.J, b_adrc.params.c_str());
  std::printf("best stsmc_dual     J=%7.3f  %s\n", b_smc.J, b_smc.params.c_str());
  std::printf("best lqg            J=%7.3f  %s\n", b_lqg.J, b_lqg.params.c_str());
  return 0;
}

}  // namespace

int run(const std::string& name, const Options& o) {
  if (name == "tune") return exp_tune(o);
  if (name.rfind("single:", 0) == 0) {  // single:<config> runs one tracking scenario
    const Metrics m = simulate(PlantParams{}, named_config(name.substr(7), o), tracking_scenario());
    std::printf("%s stable=%d rms=%.3f max=%.3f ss=%.3f dist_peak=%.2f mrad recovery=%.0f ms "
                "peak_iq=%.1f A chatter=%.4f Nm\n", name.c_str() + 7, m.stable,
                1e3 * m.rms_err_motion, 1e3 * m.max_err_motion, 1e3 * m.ss_err, 1e3 * m.dist_peak,
                1e3 * m.dist_recovery, m.peak_iq, m.chatter);
    return 0;
  }
  if (name == "lqr_info") {
    const ControllerConfig c = named_config("lqg", o);
    const LqgReport r = lqg_report(c.lqg, NominalModel::from(PlantParams{}), Rates{}.f_vel);
    auto pr = [](const char* label, const Eigen::VectorXcd& p) {
      std::printf("%s\n", label);
      for (Eigen::Index i = 0; i < p.size(); ++i) {
        const double wn = std::abs(p(i)), zeta = wn > 0 ? -p(i).real() / wn : 1.0;
        std::printf("   s = %10.2f %+10.2fj   f_n = %7.2f Hz  zeta = %.3f\n", p(i).real(),
                    p(i).imag(), wn / (2.0 * kPi), zeta);
      }
    };
    std::printf("LQR gain K = [%.4g %.4g %.4g %.4g %.4g]\n", r.K(0), r.K(1), r.K(2), r.K(3), r.K(4));
    std::printf("DARE iterations: controller %d, estimator %d. Full synthesis %.0f us\n",
                r.dare_iters_ctrl, r.dare_iters_est, r.synth_time_us);
    pr("open-loop plant poles:", r.poles_open);
    pr("LQR closed-loop poles:", r.poles_ctrl);
    pr("Kalman estimator poles:", r.poles_est);
    return 0;
  }
  if (name == "track") return exp_track(o);
  if (name == "chirp_plant") return exp_chirp_plant(o);
  if (name == "loopgain") return exp_loopgain(o);
  if (name == "sweep") return exp_sweep(o);
  if (name == "montecarlo") return exp_montecarlo(o);
  if (name == "adapt") return exp_adapt(o);
  throw std::invalid_argument("unknown experiment " + name);
}

}  // namespace actlab::exp
