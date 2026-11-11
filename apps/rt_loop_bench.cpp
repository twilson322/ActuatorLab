// Real-time loop harness: host-side 1 kHz position loop talking to an emulated drive.
//
//   [host thread, 1 kHz]  PositionLoop  --cmd-->  TripleBuffer  --> [drive thread, 1 kHz wall]
//                                       <--state-- TripleBuffer <--  plant + FOC + speed loop
//
// The drive thread advances the simulated actuator by 1 ms of simulated time per 1 ms of wall
// time (inner loops still run at 20 kHz / 5 kHz of simulated time inside that millisecond),
// the way an EtherCAT drive exchanges process data with a host once per cycle. Measured:
// host wake-up latency, compute time, deadline misses, command age at the drive, and watchdog
// trips (the drive falls back to a gravity-holding safe mode if commands stop arriving).
//
//   rt_loop_bench [--seconds S] [--stress N] [--fifo 0|1] [--out results/rt_bench.csv]
#include <pthread.h>
#include <sched.h>
#include <sys/mman.h>
#include <time.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

#include "actlab/concurrency.hpp"
#include "actlab/controllers.hpp"
#include "actlab/drive.hpp"
#include "actlab/plant.hpp"
#include "actlab/sim.hpp"
#include "actlab/telemetry.hpp"

using namespace actlab;

namespace {

constexpr long kPeriodNs = 1'000'000;  // 1 kHz

struct StateMsg {
  std::uint64_t seq{0};
  double t_sim{0};
  Feedback fb{};
  double th_l_true{0};
  bool safe_mode{false};
};
struct CmdMsg {
  std::uint64_t seq{0};
  std::uint64_t state_seq{0};  // which state sample this command was computed from
  VelSetpoint sp{};
};

std::int64_t now_ns() {
  timespec ts{};
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return static_cast<std::int64_t>(ts.tv_sec) * 1'000'000'000 + ts.tv_nsec;
}
void sleep_until_ns(std::int64_t t) {
  timespec ts{static_cast<time_t>(t / 1'000'000'000), static_cast<long>(t % 1'000'000'000)};
  while (clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &ts, nullptr) == EINTR) {}
}
bool try_fifo(int prio) {
  sched_param sp{};
  sp.sched_priority = prio;
  return pthread_setschedparam(pthread_self(), SCHED_FIFO, &sp) == 0;
}

double percentile(std::vector<double> v, double p) {
  if (v.empty()) return 0.0;
  std::sort(v.begin(), v.end());
  const std::size_t i = std::min(v.size() - 1, static_cast<std::size_t>(p * static_cast<double>(v.size())));
  return v[i];
}

}  // namespace

int main(int argc, char** argv) {
  double seconds = 10.0;
  int stress = 0;
  bool want_fifo = true;
  std::string out = "results/rt_bench.csv";
  for (int i = 1; i + 1 < argc; i += 2) {
    const std::string k = argv[i];
    if (k == "--seconds") seconds = std::atof(argv[i + 1]);
    else if (k == "--stress") stress = std::atoi(argv[i + 1]);
    else if (k == "--fifo") want_fifo = std::atoi(argv[i + 1]) != 0;
    else if (k == "--out") out = argv[i + 1];
  }
  const bool locked = mlockall(MCL_CURRENT | MCL_FUTURE) == 0;

  const PlantParams pp{};
  const Rates rates{};
  const NominalModel nominal = NominalModel::from(pp);
  TripleBuffer<StateMsg> state_buf;
  TripleBuffer<CmdMsg> cmd_buf;
  std::atomic<bool> run{true};
  std::atomic<int> fifo_granted{0};
  std::atomic<std::uint64_t> watchdog_trips{0};
  std::vector<double> cmd_age_ms;  // written only by the drive thread, read after join
  cmd_age_ms.reserve(static_cast<std::size_t>(seconds * 1000) + 16);

  Trajectory traj(0.0);
  for (double t = 0.2; t + 1.0 < seconds; t += 1.2) traj.move(t, t + 0.4, 1.0).move(t + 0.6, t + 1.0, 0.0);

  // ------------------------------------------------------------------ drive thread
  std::thread drive_thread([&] {
    if (want_fifo && try_fifo(80)) fifo_granted.fetch_add(1);
    Plant plant(pp, 3);
    Drive drive(pp.motor, pp.drive, DriveConfig{}, rates.f_cur);
    VelocityLoopConfig vc;
    vc.bw_hz = 40.0;
    VelocityLoop vel(vc, nominal.J_tot_motor(), rates.f_vel);
    const double dt = 1.0 / rates.f_sim, tau_max = pp.motor.kt() * pp.drive.i_max;
    Sensors s = plant.sample();
    drive.reset(s);
    drive.sense(s);
    foc::AlphaBeta pending{};
    double iq_ref = 0.0, t_sim = 0.0;
    std::uint64_t seq = 0, last_cmd_seq = 0;
    int stale = 0;
    bool safe = false;
    VelSetpoint sp{};
    long k = 0;
    std::int64_t next = now_ns();
    while (run.load(std::memory_order_acquire)) {
      // ---- receive the newest command (or count it as stale)
      if (cmd_buf.update()) {
        const CmdMsg& c = cmd_buf.read_buffer();
        if (c.seq != last_cmd_seq) {
          last_cmd_seq = c.seq;
          stale = 0;
          safe = false;
          sp = c.sp;
          cmd_age_ms.push_back(static_cast<double>(seq - c.state_seq));
        }
      } else if (++stale > 3 && !safe) {
        safe = true;  // watchdog: no fresh command for 3 cycles -> hold position against gravity
        watchdog_trips.fetch_add(1);
      }
      if (safe) {
        sp.w_ref = 0.0;
        sp.tau_ff = nominal.mgl * std::sin(s.th_l) / nominal.N;
      }
      // ---- advance 1 ms of simulated time with the inner loops at their own rates
      for (int i = 0; i < static_cast<int>(rates.f_sim / 1000.0); ++i, ++k) {
        if (k % 2 == 0) {
          plant.set_phase_voltage(pending);
          s = plant.sample();
          drive.sense(s);
          if ((k / 2) % 4 == 0) {
            const double tu = vel.update(sp, drive.w_m());
            const double ts = std::clamp(tu, -tau_max, tau_max);
            vel.commit(tu, ts);
            iq_ref = ts / pp.motor.kt();
          }
          pending = drive.control(iq_ref);
        }
        plant.step(dt);
        t_sim += dt;
      }
      // ---- publish state
      StateMsg& m = state_buf.write_buffer();
      m.seq = ++seq;
      m.t_sim = t_sim;
      m.fb = {drive.th_m(), drive.w_m(), s.th_l};
      m.th_l_true = plant.state().th_l;
      m.safe_mode = safe;
      state_buf.publish();
      next += kPeriodNs;
      sleep_until_ns(next);
    }
  });

  // ------------------------------------------------------------------ CPU stress (optional)
  std::vector<std::thread> stressors;
  for (int i = 0; i < stress; ++i)
    stressors.emplace_back([&] {
      volatile double x = 1.0;
      while (run.load(std::memory_order_relaxed)) x = std::sin(x) + 1.0;
    });

  // ------------------------------------------------------------------ host thread (this one)
  if (want_fifo && try_fifo(70)) fifo_granted.fetch_add(1);
  PositionLoop pos(2.0 * kPi * 40.0, nominal, true);  // same gains as motor_defl
  const long n_cycles = static_cast<long>(seconds * 1000.0);
  std::vector<double> wake_us, compute_us;
  wake_us.reserve(static_cast<std::size_t>(n_cycles));
  compute_us.reserve(static_cast<std::size_t>(n_cycles));
  long misses = 0;
  double sum_sq = 0.0;
  long n_err = 0;
  {
    CsvTelemetry log(out, {"cycle", "t_sim", "ref", "th_l", "err", "wake_us", "compute_us", "safe"});
    std::uint64_t cmd_seq = 0;
    std::int64_t next = now_ns() + kPeriodNs;
    for (long c = 0; c < n_cycles; ++c) {
      sleep_until_ns(next);
      const std::int64_t woke = now_ns();
      const double lat = 1e-3 * static_cast<double>(woke - next);
      if (woke - next > kPeriodNs) ++misses;
      state_buf.update();
      const StateMsg st = state_buf.read_buffer();
      const Reference r = traj.eval(st.t_sim);
      CmdMsg& cm = cmd_buf.write_buffer();
      cm.seq = ++cmd_seq;
      cm.state_seq = st.seq;
      cm.sp = pos.update(r, st.fb);
      cmd_buf.publish();
      const double comp = 1e-3 * static_cast<double>(now_ns() - woke);
      wake_us.push_back(lat);
      compute_us.push_back(comp);
      const double e = st.th_l_true - r.pos;
      if (st.t_sim > 0.1) { sum_sq += e * e; ++n_err; }
      LogRow row;
      const double v[8] = {static_cast<double>(c), st.t_sim, r.pos, st.th_l_true, e, lat, comp,
                           st.safe_mode ? 1.0 : 0.0};
      std::copy(std::begin(v), std::end(v), row.v);
      log.try_log(row);  // never blocks the control thread; drops are counted
      next += kPeriodNs;
      if (c == n_cycles - 1) {
        std::printf("telemetry rows dropped: %zu\n", log.dropped());
      }
    }
  }  // CsvTelemetry destructor drains and closes the file here
  run.store(false, std::memory_order_release);
  drive_thread.join();
  for (auto& t : stressors) t.join();

  std::printf("rt_loop_bench: %.1f s, stress threads=%d, SCHED_FIFO granted to %d/2 threads, "
              "mlockall=%s\n", seconds, stress, fifo_granted.load(), locked ? "yes" : "no");
  std::printf("  wake latency  p50 %.1f us  p99 %.1f us  p99.9 %.1f us  max %.1f us\n",
              percentile(wake_us, 0.5), percentile(wake_us, 0.99), percentile(wake_us, 0.999),
              *std::max_element(wake_us.begin(), wake_us.end()));
  std::printf("  compute time  p50 %.2f us  p99 %.2f us  max %.2f us\n", percentile(compute_us, 0.5),
              percentile(compute_us, 0.99), *std::max_element(compute_us.begin(), compute_us.end()));
  std::printf("  deadline misses (late by > 1 period): %ld of %ld\n", misses, n_cycles);
  std::printf("  command age at drive  p50 %.0f  p99 %.0f  max %.0f cycles\n",
              percentile(cmd_age_ms, 0.5), percentile(cmd_age_ms, 0.99),
              cmd_age_ms.empty() ? 0.0 : *std::max_element(cmd_age_ms.begin(), cmd_age_ms.end()));
  std::printf("  watchdog trips: %llu\n", static_cast<unsigned long long>(watchdog_trips.load()));
  std::printf("  link tracking RMS: %.3f mrad\n", 1e3 * std::sqrt(sum_sq / static_cast<double>(std::max(1L, n_err))));
  return 0;
}
