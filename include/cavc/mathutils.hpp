#ifndef CAVC_MATHUTILS_HPP
#define CAVC_MATHUTILS_HPP

#include "internal/common.hpp"
#include <atomic>
#include <cmath>
#include <iterator>

namespace cavc {
namespace utils {
template <typename Real> constexpr Real defaultRealThreshold() { return Real(1e-8); }

template <typename Real> constexpr Real defaultRealPrecision() { return Real(1e-5); }

template <typename Real> constexpr Real defaultSliceJoinThreshold() { return Real(1e-4); }

template <typename Real> constexpr Real defaultOffsetThreshold() { return Real(1e-4); }

template <typename Real> struct EpsilonConfig {
  Real realThreshold;
  Real realPrecision;
  Real sliceJoinThreshold;
  Real offsetThreshold;
};

template <typename Real> constexpr EpsilonConfig<Real> defaultEpsilonConfig() {
  return {defaultRealThreshold<Real>(), defaultRealPrecision<Real>(),
          defaultSliceJoinThreshold<Real>(), defaultOffsetThreshold<Real>()};
}

template <typename Real> struct AtomicEpsilonConfig {
  std::atomic<Real> realThreshold;
  std::atomic<Real> realPrecision;
  std::atomic<Real> sliceJoinThreshold;
  std::atomic<Real> offsetThreshold;

  explicit AtomicEpsilonConfig(EpsilonConfig<Real> const &config)
      : realThreshold(config.realThreshold), realPrecision(config.realPrecision),
        sliceJoinThreshold(config.sliceJoinThreshold), offsetThreshold(config.offsetThreshold) {}
};

template <typename Real> AtomicEpsilonConfig<Real> &epsilonConfig() {
  static AtomicEpsilonConfig<Real> config(defaultEpsilonConfig<Real>());
  return config;
}

template <typename Real> EpsilonConfig<Real> getEpsilonConfig() {
  auto &config = epsilonConfig<Real>();
  return {config.realThreshold.load(std::memory_order_relaxed),
          config.realPrecision.load(std::memory_order_relaxed),
          config.sliceJoinThreshold.load(std::memory_order_relaxed),
          config.offsetThreshold.load(std::memory_order_relaxed)};
}

template <typename Real> void setEpsilonConfig(EpsilonConfig<Real> const &config) {
  CAVC_ASSERT(config.realThreshold > Real(0), "realThreshold must be > 0");
  CAVC_ASSERT(config.realPrecision > Real(0), "realPrecision must be > 0");
  CAVC_ASSERT(config.sliceJoinThreshold > Real(0), "sliceJoinThreshold must be > 0");
  CAVC_ASSERT(config.offsetThreshold > Real(0), "offsetThreshold must be > 0");
  auto &current = epsilonConfig<Real>();
  current.realThreshold.store(config.realThreshold, std::memory_order_relaxed);
  current.realPrecision.store(config.realPrecision, std::memory_order_relaxed);
  current.sliceJoinThreshold.store(config.sliceJoinThreshold, std::memory_order_relaxed);
  current.offsetThreshold.store(config.offsetThreshold, std::memory_order_relaxed);
}

template <typename Real> void resetEpsilonConfig() {
  setEpsilonConfig(defaultEpsilonConfig<Real>());
}

// absolute threshold to be used for comparing reals generally
template <typename Real> Real realThreshold() {
  return epsilonConfig<Real>().realThreshold.load(std::memory_order_relaxed);
}

// absolute threshold to be used for reals in common geometric computation (e.g. to check for
// singularities)
template <typename Real> Real realPrecision() {
  return epsilonConfig<Real>().realPrecision.load(std::memory_order_relaxed);
}

// absolute threshold to be used for joining slices together at end points
template <typename Real> Real sliceJoinThreshold() {
  return epsilonConfig<Real>().sliceJoinThreshold.load(std::memory_order_relaxed);
}

// absolute threshold to be used for pruning invalid slices for offset
template <typename Real> Real offsetThreshold() {
  return epsilonConfig<Real>().offsetThreshold.load(std::memory_order_relaxed);
}

template <typename Real> constexpr Real pi() { return Real(3.14159265358979323846264338327950288); }

template <typename Real> constexpr Real tau() { return Real(2) * pi<Real>(); }

template <typename Real> bool fuzzyEqual(Real x, Real y, Real epsilon = realThreshold<Real>()) {
  return std::abs(x - y) < epsilon;
}

template <typename Real>
bool fuzzyInRange(Real minValue, Real value, Real maxValue, Real epsilon = realThreshold<Real>()) {
  return (value + epsilon > minValue) && (value < maxValue + epsilon);
}

/// Normalize radius to be between 0 and 2PI, e.g. -PI/4 becomes 7PI/8 and 5PI becomes PI.
template <typename Real> Real normalizeRadians(Real angle) {
  if (angle >= Real(0) && angle <= tau<Real>()) {
    return angle;
  }

  return angle - std::floor(angle / tau<Real>()) * tau<Real>();
}

/// Returns the smaller difference between two angles, result is negative if angle2 < angle1.
template <typename Real> Real deltaAngle(Real angle1, Real angle2) {
  Real diff = normalizeRadians(angle2 - angle1);
  if (diff > pi<Real>()) {
    diff -= tau<Real>();
  }

  return diff;
}

/// Tests if angle is between a start and end angle (counter clockwise start to end, inclusive).
template <typename Real>
bool angleIsBetween(Real startAngle, Real endAngle, Real testAngle,
                    Real epsilon = realThreshold<Real>()) {
  Real endSweep = normalizeRadians(endAngle - startAngle);
  Real midSweep = normalizeRadians(testAngle - startAngle);

  return midSweep < endSweep + epsilon;
}

template <typename Real>
bool angleIsWithinSweep(Real startAngle, Real sweepAngle, Real testAngle,
                        Real epsilon = realThreshold<Real>()) {
  Real endAngle = startAngle + sweepAngle;
  if (sweepAngle < Real(0)) {
    return angleIsBetween(endAngle, startAngle, testAngle, epsilon);
  }

  return angleIsBetween(startAngle, endAngle, testAngle, epsilon);
}

/// Returns the solutions to for the quadratic equation -b +/- sqrt (b * b - 4 * a * c) / (2 * a).
template <typename Real>
std::pair<Real, Real> quadraticSolutions(Real a, Real b, Real c, Real discr) {
  // Function avoids loss in precision due to taking the difference of two floating point values
  // that are very near each other in value.
  // See:
  // https://math.stackexchange.com/questions/311382/solving-a-quadratic-equation-with-precision-when-using-floating-point-variables
  CAVC_ASSERT(fuzzyEqual(b * b - Real(4) * a * c, discr), "discriminate is not correct");
  Real sqrtDiscr = std::sqrt(discr);
  Real denom = Real(2) * a;
  Real sol1;
  if (b < Real(0)) {
    sol1 = (-b + sqrtDiscr) / denom;
  } else {
    sol1 = (-b - sqrtDiscr) / denom;
  }

  Real sol2 = (c / a) / sol1;

  return std::make_pair(sol1, sol2);
}

template <typename T> std::size_t nextWrappingIndex(std::size_t index, const T &container) {
  if (index == container.size() - 1) {
    return 0;
  }

  return index + 1;
}

template <typename T> std::size_t prevWrappingIndex(std::size_t index, const T &container) {
  if (index == 0) {
    return container.size() - 1;
  }

  return index - 1;
}
} // namespace utils
} // namespace cavc

#endif // CAVC_MATHUTILS_HPP
