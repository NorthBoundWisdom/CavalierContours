#include "benchmarkprofiles.h"
#include "cavc/polylinecombine.hpp"
#include <benchmark/benchmark.h>

struct CombineShiftedSetup {
  std::vector<cavccpp::Polyline<double>> shiftedProfiles;

  CombineShiftedSetup(TestProfile const &profile) {
    auto extents = cavccpp::getExtents(profile.pline);
    double halfWidth = (extents.xMax - extents.xMin) / 2.0;
    double halfHeight = (extents.yMax - extents.yMin) / 2.0;

    auto createShifted = [&](double angle) {
      auto shifted = profile.pline;
      double xOffset = halfWidth * std::cos(angle);
      double yOffset = halfHeight * std::sin(angle);
      cavccpp::translatePolyline(shifted, {xOffset, yOffset});
      return shifted;
    };

    std::size_t shiftedCount = 16;
    shiftedProfiles.reserve(shiftedCount);
    for (std::size_t i = 0; i < shiftedCount; ++i) {
      double angle = static_cast<double>(i) / static_cast<double>(shiftedCount) *
                     cavccpp::utils::tau<double>();
      shiftedProfiles.push_back(createShifted(angle));
    }
  }
};

static void combineShifted(CombineShiftedSetup const &setup, TestProfile const &profile) {
  for (auto const &shifted : setup.shiftedProfiles) {
    cavccpp::combinePolylines(profile.pline, shifted, cavccpp::PlineCombineMode::Union);
    cavccpp::combinePolylines(profile.pline, shifted, cavccpp::PlineCombineMode::Exclude);
    cavccpp::combinePolylines(profile.pline, shifted, cavccpp::PlineCombineMode::Intersect);
    cavccpp::combinePolylines(profile.pline, shifted, cavccpp::PlineCombineMode::XOR);
  }
}

CAVC_CREATE_BENCHMARKS(combine16Shifted, CombineShiftedSetup, combineShifted,
                       benchmark::kMicrosecond)
CAVC_CREATE_NO_ARCS_BENCHMARKS(combine16Shifted, CombineShiftedSetup, combineShifted, 0.01,
                               benchmark::kMicrosecond)

static void combineCoincident(NoSetup, TestProfile const &profile) {
  cavccpp::combinePolylines(profile.pline, profile.pline, cavccpp::PlineCombineMode::Union);
  cavccpp::combinePolylines(profile.pline, profile.pline, cavccpp::PlineCombineMode::Exclude);
  cavccpp::combinePolylines(profile.pline, profile.pline, cavccpp::PlineCombineMode::Intersect);
  cavccpp::combinePolylines(profile.pline, profile.pline, cavccpp::PlineCombineMode::XOR);
}

CAVC_CREATE_BENCHMARKS(combineCoincident, NoSetup, combineCoincident, benchmark::kMicrosecond)
CAVC_CREATE_NO_ARCS_BENCHMARKS(combineCoincident, NoSetup, combineCoincident, 0.01,
                               benchmark::kMicrosecond)

BENCHMARK_MAIN();
