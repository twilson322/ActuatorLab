#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "actlab/output_rls.hpp"
#include "actlab/rls.hpp"

using actlab::RLS;

TEST(Rls, ConvergesWithNoise) {
  std::mt19937 g(1);
  std::normal_distribution<double> n(0.0, 1.0);
  const Eigen::Vector3d truth(0.75, -2.0, 11.8);
  RLS<3> rls(Eigen::Vector3d::Zero(), 1e3, 1.0);
  for (int k = 0; k < 2000; ++k) {
    const Eigen::Vector3d phi(n(g), n(g), std::sin(0.01 * k));
    rls.update(phi, phi.dot(truth) + 0.01 * n(g));
  }
  EXPECT_LT((rls.theta() - truth).norm(), 0.01);
}

TEST(Rls, ForgettingTracksAParameterJump) {
  std::mt19937 g(2);
  std::normal_distribution<double> n(0.0, 1.0);
  Eigen::Vector2d truth(1.0, 2.0);
  RLS<2> rls(Eigen::Vector2d::Zero(), 1e3, 0.99);
  for (int k = 0; k < 4000; ++k) {
    if (k == 2000) truth << 3.0, -1.0;  // payload change
    const Eigen::Vector2d phi(n(g), n(g));
    rls.update(phi, phi.dot(truth));
  }
  EXPECT_LT((rls.theta() - truth).norm(), 1e-3);
}

TEST(Rls, CovarianceStaysBoundedWithoutExcitation) {
  RLS<2> rls(Eigen::Vector2d::Zero(), 1.0, 0.98, 1e4);
  for (int k = 0; k < 20000; ++k) rls.update(Eigen::Vector2d(1.0, 1.0), 2.0);  // rank-1 regressor
  EXPECT_LE(rls.covariance().trace(), 1e4 * (1.0 + 1e-9));
  EXPECT_TRUE(rls.theta().allFinite());
}

TEST(Rls, OutputSideEstimatorIdentifiesARigidJoint) {
  // A rigid joint driven by a known torque profile: the shared front end of the adaptive cascade
  // and the adaptive LQG must recover inertia and gravity (friction is left at zero here).
  const double J = 1.2, mgl = 15.0, fs = 1000.0;
  actlab::OutputSideRls est(actlab::OutputSideRls::Vec(0.6, 0.0, 0.0, 5.0), fs);
  double th = 0.0, w = 0.0;
  est.reset(th, 0.0);
  int learned = 0;
  for (int k = 0; k < 20000; ++k) {
    const double t = k / fs;
    const double tau = mgl * std::sin(th) + 8.0 * std::sin(2.0 * actlab::kPi * 0.7 * t) + 5.0 * std::sin(2.0 * actlab::kPi * 1.9 * t);
    constexpr int kSubsteps = 20;  // integrate the test plant finely so its error does not bias the fit
    for (int i = 0; i < kSubsteps; ++i) {
      w += (tau - mgl * std::sin(th)) / J / (fs * kSubsteps);
      th += w / (fs * kSubsteps);
    }
    learned += est.update(th, tau);
  }
  EXPECT_GT(learned, 10000);
  // With the end-of-sample acceleration output this gave J = 1.279 (+6.5 %) and b = 3.1.
  EXPECT_NEAR(est.theta()(0), J, 0.01 * J);
  EXPECT_NEAR(est.theta()(1), 0.0, 0.1);  // no viscous friction in this plant
  EXPECT_NEAR(est.theta()(3), mgl, 0.01 * mgl);
}
