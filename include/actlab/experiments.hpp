#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "actlab/controllers.hpp"
#include "actlab/sim.hpp"

namespace actlab::exp {

struct Options {
  std::string out_dir{"results"};
  std::string identified_path{"results/identified.txt"};
  int mc_samples{200};
  int threads{2};
  std::uint64_t seed{7};
  std::vector<std::pair<std::string, double>> overrides;  // --set key=value
};

// Named controller configurations used across experiments. The notch centre comes from
// the identification step (results/identified.txt) when that file exists.
std::vector<std::pair<std::string, ControllerConfig>> controller_suite(const Options& o);
ControllerConfig named_config(const std::string& name, const Options& o);

// The standard tracking + disturbance scenario all controllers are scored on.
SimConfig tracking_scenario();

int run(const std::string& name, const Options& o);

}  // namespace actlab::exp
