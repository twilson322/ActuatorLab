#pragma once
#include <Eigen/Dense>

namespace actlab {

// Recursive least squares with exponential forgetting for y = phi^T theta + noise.
// Fixed-size Eigen types keep it allocation-free so it can live inside a control loop.
// Covariance is symmetrized every update and its trace is capped to stop the classic
// "covariance windup" when excitation disappears while forgetting is active.
template <int N>
class RLS {
 public:
  using Vec = Eigen::Matrix<double, N, 1>;
  using Mat = Eigen::Matrix<double, N, N>;

  RLS(const Vec& theta0, double p0, double lambda, double trace_max = 1e8)
      : theta_(theta0), P_(Mat::Identity() * p0), lambda_(lambda), trace_max_(trace_max) {}

  // Returns the a-priori prediction error.
  double update(const Vec& phi, double y) {
    const double err = y - phi.dot(theta_);
    const Vec Pphi = P_ * phi;
    const double denom = lambda_ + phi.dot(Pphi);
    const Vec K = Pphi / denom;
    theta_ += K * err;
    P_ = (P_ - K * Pphi.transpose()) / lambda_;
    P_ = 0.5 * (P_ + P_.transpose());
    const double tr = P_.trace();
    if (tr > trace_max_) P_ *= trace_max_ / tr;
    return err;
  }

  const Vec& theta() const { return theta_; }
  const Mat& covariance() const { return P_; }

 private:
  Vec theta_;
  Mat P_;
  double lambda_;
  double trace_max_;
};

}  // namespace actlab
