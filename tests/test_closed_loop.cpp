#include <gtest/gtest.h>

#include <Eigen/Eigenvalues>

#include "actlab/experiments.hpp"
#include "actlab/lqg.hpp"

using namespace actlab;

TEST(Lqg, ZohMatchesDoubleIntegrator) {
  Eigen::MatrixXd A(2, 2), B(2, 1), Ad, Bd;
  A << 0, 1, 0, 0;
  B << 0, 1;
  const double T = 0.01;
  c2d_zoh(A, B, T, Ad, Bd);
  EXPECT_NEAR(Ad(0, 1), T, 1e-14);
  EXPECT_NEAR(Bd(0, 0), 0.5 * T * T, 1e-14);
  EXPECT_NEAR(Bd(1, 0), T, 1e-14);
}

TEST(Lqg, DareSatisfiesRiccatiEquation) {
  const PlantParams pp;
  const NominalModel m = NominalModel::from(pp);
  const auto sys = TwoMassLinear::build(m, 0.15, 0.05, 2.0);
  Eigen::MatrixXd A = sys.A.topLeftCorner(4, 4), B = sys.B.topRows(4), Ad, Bd;
  c2d_zoh(A, B, 2e-4, Ad, Bd);
  Eigen::MatrixXd Q = Eigen::Vector4d(1e3, 1.0, 1e6, 1e2).asDiagonal();
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(1, 1);
  int it = 0;
  const Eigen::MatrixXd P = dare(Ad, Bd, Q, R, 100, 1e-12, &it);
  const Eigen::MatrixXd S = R + Bd.transpose() * P * Bd;
  const Eigen::MatrixXd res = Ad.transpose() * P * Ad -
                              Ad.transpose() * P * Bd * S.inverse() * Bd.transpose() * P * Ad + Q - P;
  EXPECT_LT(res.norm() / P.norm(), 1e-9);
  EXPECT_LT(it, 40);
  const Eigen::MatrixXd K = S.ldlt().solve(Bd.transpose() * P * Ad);
  const auto ev = Eigen::EigenSolver<Eigen::MatrixXd>(Ad - Bd * K).eigenvalues();
  for (int i = 0; i < ev.size(); ++i) EXPECT_LT(std::abs(ev(i)), 1.0);
}

TEST(ClosedLoop, TunedControllersTrackOnNominalPlant) {
  exp::Options o;
  o.identified_path = "nonexistent";  // use defaults, independent of results/
  for (const char* name : {"motor_defl", "dual_pi_notch", "lqg"}) {
    Plant plant(PlantParams{}, 1);
    auto c = make_controller(exp::named_config(name, o), NominalModel::from(PlantParams{}), 5000.0, 1000.0);
    const Metrics m = run_sim(plant, *c, exp::tracking_scenario());
    EXPECT_TRUE(m.stable) << name;
    EXPECT_LT(m.rms_err_motion, 1e-3) << name;  // < 1 mrad RMS on a 90 deg move
  }
}

TEST(ClosedLoop, LqgRejectsLoadFasterThanCascade) {
  exp::Options o;
  o.identified_path = "nonexistent";
  auto run = [&](const char* name) {
    Plant plant(PlantParams{}, 1);
    auto c = make_controller(exp::named_config(name, o), NominalModel::from(PlantParams{}), 5000.0, 1000.0);
    return run_sim(plant, *c, exp::tracking_scenario());
  };
  const Metrics casc = run("motor_defl"), lqg = run("lqg");
  EXPECT_LT(lqg.dist_peak, 0.5 * casc.dist_peak);
  EXPECT_LT(lqg.dist_recovery, 0.1);
}

TEST(Lqg, KalmanFilterFollowsTheEncoderResolution) {
  // The measurement noise used to be computed from hard-coded 14- and 17-bit encoders. A finer
  // link encoder must now make the filter trust the link measurement more.
  const NominalModel m = NominalModel::from(PlantParams{});
  LqgConfig coarse, fine;
  fine.encoder_link_bits = coarse.encoder_link_bits + 4;
  const LqgReport a = lqg_report(coarse, m, Rates{}.f_vel), b = lqg_report(fine, m, Rates{}.f_vel);
  EXPECT_GT(b.L.col(1).norm(), a.L.col(1).norm());  // column 1: gain on the link encoder
}
