#include <gtest/gtest.h>

#include <cmath>

#include "actlab/foc.hpp"

using namespace actlab::foc;

TEST(Foc, ClarkeParkRoundTrip) {
  for (double th = -7.0; th < 7.0; th += 0.37) {
    const DQ dq{1.3, -2.1};
    const ABC abc = inverse_clarke(inverse_park(dq, th));
    EXPECT_NEAR(abc.a + abc.b + abc.c, 0.0, 1e-12);
    const DQ back = park(clarke(abc.a, abc.b), th);
    EXPECT_NEAR(back.d, dq.d, 1e-12);
    EXPECT_NEAR(back.q, dq.q, 1e-12);
  }
}

TEST(Foc, AmplitudeInvariant) {
  // Balanced currents of amplitude 5 A map to |i_ab| = 5.
  const double th = 0.8;
  const double ia = 5.0 * std::cos(th), ib = 5.0 * std::cos(th - 2.0 * M_PI / 3.0);
  const AlphaBeta ab = clarke(ia, ib);
  EXPECT_NEAR(std::hypot(ab.alpha, ab.beta), 5.0, 1e-12);
}

TEST(Foc, SvpwmReproducesVoltageInLinearRange) {
  const double vdc = 48.0, vmax = vdc / kSqrt3;
  for (double ang = 0.0; ang < 2 * M_PI; ang += 0.05) {
    for (double frac : {0.1, 0.5, 0.99}) {
      const AlphaBeta v{frac * vmax * std::cos(ang), frac * vmax * std::sin(ang)};
      const auto duty = svpwm(v, vdc);
      for (double d : duty) {
        EXPECT_GE(d, 0.0);
        EXPECT_LE(d, 1.0);
      }
      const AlphaBeta out = inverter_average(duty, vdc);
      EXPECT_NEAR(out.alpha, v.alpha, 1e-9);
      EXPECT_NEAR(out.beta, v.beta, 1e-9);
    }
  }
}

TEST(Foc, SvpwmSaturatesBeyondHexagon) {
  const double vdc = 48.0;
  const AlphaBeta v{2.0 * vdc, 0.0};
  const AlphaBeta out = inverter_average(svpwm(v, vdc), vdc);
  EXPECT_LT(std::hypot(out.alpha, out.beta), 2.0 * vdc / 3.0 + 1e-9);  // hexagon vertex
}

TEST(Foc, CircleLimit) {
  const DQ v = limit_circle({30.0, 40.0}, 10.0);
  EXPECT_NEAR(std::hypot(v.d, v.q), 10.0, 1e-12);
  EXPECT_NEAR(v.q / v.d, 40.0 / 30.0, 1e-12);
}
