// Command-line entry point for the offline experiments.
//   actlab_sim <track|chirp_plant|loopgain|sweep|montecarlo|adapt|all> [--out DIR]
//              [--identified FILE] [--mc N] [--threads N] [--seed S]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>

#include "actlab/experiments.hpp"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <experiment> [--out DIR] [--identified FILE] [--mc N] "
                         "[--threads N] [--seed S] [--set key=value]...\n", argv[0]);
    return 2;
  }
  actlab::exp::Options o;
  const std::string name = argv[1];
  for (int i = 2; i + 1 < argc; i += 2) {
    const std::string k = argv[i], v = argv[i + 1];
    if (k == "--out") o.out_dir = v;
    else if (k == "--identified") o.identified_path = v;
    else if (k == "--mc") o.mc_samples = std::atoi(v.c_str());
    else if (k == "--threads") o.threads = std::atoi(v.c_str());
    else if (k == "--seed") o.seed = std::strtoull(v.c_str(), nullptr, 10);
    else if (k == "--set") {
      const auto eq = v.find('=');
      if (eq == std::string::npos) { std::fprintf(stderr, "--set expects key=value\n"); return 2; }
      o.overrides.emplace_back(v.substr(0, eq), std::stod(v.substr(eq + 1)));
    }
    else { std::fprintf(stderr, "unknown option %s\n", k.c_str()); return 2; }
  }
  try {
    const auto t0 = std::chrono::steady_clock::now();
    int rc = 0;
    if (name == "all") {
      for (const char* e : {"track", "loopgain", "sweep", "montecarlo", "adapt"})
        rc |= actlab::exp::run(e, o);
    } else {
      rc = actlab::exp::run(name, o);
    }
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("[%s] finished in %.1f s\n", name.c_str(), s);
    return rc;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
