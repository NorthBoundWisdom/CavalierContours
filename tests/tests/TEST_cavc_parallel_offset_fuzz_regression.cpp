#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <cavc/polyline.hpp>
#include <cavc/polylineoffset.hpp>
#include <gtest/gtest.h>

namespace {
using Pline = cavc::Polyline<double>;

bool allVertexesFinite(Pline const &pline) {
  return std::all_of(pline.vertexes().begin(), pline.vertexes().end(), [](auto const &v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.bulge());
  });
}

void expectBasicOffsetInvariants(Pline const &input, std::vector<Pline> const &results,
                                 std::string const &caseName) {
  ASSERT_FALSE(results.empty()) << caseName;
  for (std::size_t i = 0; i < results.size(); ++i) {
    auto const &result = results[i];
    SCOPED_TRACE(caseName + " result=" + std::to_string(i));
    EXPECT_EQ(result.isClosed(), input.isClosed());
    EXPECT_TRUE(allVertexesFinite(result));
    EXPECT_GT(cavc::getPathLength(result), cavc::utils::realThreshold<double>());
    if (result.isClosed()) {
      EXPECT_GE(result.size(), 2u);
      EXPECT_TRUE(std::isfinite(cavc::getArea(result)));
    } else {
      EXPECT_GE(result.size(), 2u);
    }
  }
}

Pline makeOpenMixedPolyline(std::uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> stepDist(0.5, 1.8);
  std::uniform_real_distribution<double> yDist(-0.8, 0.8);
  std::uniform_real_distribution<double> bulgeDist(-0.45, 0.45);

  Pline pline;
  double x = 0.0;
  double y = 0.0;
  for (int i = 0; i < 7; ++i) {
    double bulge = (i % 2 == 0) ? bulgeDist(rng) : 0.0;
    pline.addVertex(x, y, bulge);
    x += stepDist(rng);
    y += yDist(rng);
  }
  pline.lastVertex().bulge() = 0.0;
  return pline;
}

Pline makeConvexNgon(std::uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> radiusJitter(-0.45, 0.45);
  std::uniform_real_distribution<double> centerJitter(-0.2, 0.2);

  Pline pline;
  pline.isClosed() = true;
  constexpr int vertexCount = 9;
  double const cx = centerJitter(rng);
  double const cy = centerJitter(rng);
  for (int i = 0; i < vertexCount; ++i) {
    double const angle = static_cast<double>(i) * cavc::utils::tau<double>() /
                         static_cast<double>(vertexCount);
    double const radius = 5.0 + radiusJitter(rng);
    pline.addVertex(cx + radius * std::cos(angle), cy + radius * std::sin(angle), 0.0);
  }
  return pline;
}

Pline makeConcaveStar(std::uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> jitter(-0.15, 0.15);

  Pline pline;
  pline.isClosed() = true;
  constexpr int vertexCount = 10;
  for (int i = 0; i < vertexCount; ++i) {
    double const angle = static_cast<double>(i) * cavc::utils::tau<double>() /
                         static_cast<double>(vertexCount);
    double const radius = ((i % 2) == 0 ? 5.0 : 2.6) + jitter(rng);
    pline.addVertex(radius * std::cos(angle), radius * std::sin(angle), 0.0);
  }
  return pline;
}

Pline makeArcHeavyClosedPolyline(std::uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> radiusJitter(-0.25, 0.25);
  std::uniform_real_distribution<double> bulgeJitter(-0.08, 0.08);

  Pline pline;
  pline.isClosed() = true;
  constexpr int vertexCount = 8;
  for (int i = 0; i < vertexCount; ++i) {
    double const angle = static_cast<double>(i) * cavc::utils::tau<double>() /
                         static_cast<double>(vertexCount);
    double const radius = 4.0 + radiusJitter(rng);
    double const bulge = 0.18 + bulgeJitter(rng);
    pline.addVertex(radius * std::cos(angle), radius * std::sin(angle), bulge);
  }
  return pline;
}

Pline makeNearDegenerateOpenPolyline(std::uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> jitter(-1e-7, 1e-7);

  Pline pline;
  pline.addVertex(0.0, 0.0, 0.0);
  pline.addVertex(jitter(rng), jitter(rng), 0.0);
  pline.addVertex(1.0, 0.0, 0.25);
  pline.addVertex(2.0, 0.3, -0.2);
  pline.addVertex(3.1, -0.2, 0.0);
  return pline;
}

} // namespace

TEST(ParallelOffsetFuzzRegression, FixedSeedOpenJoinEndCapCorpusKeepsBasicInvariants) {
  constexpr std::array<std::uint32_t, 6> seeds = {17u, 29u, 43u, 71u, 113u, 191u};
  constexpr std::array<double, 4> offsets = {-0.65, -0.2, 0.2, 0.65};
  constexpr std::array<cavc::OffsetJoinType, 3> joinTypes = {
      cavc::OffsetJoinType::Round, cavc::OffsetJoinType::Miter, cavc::OffsetJoinType::Bevel};
  constexpr std::array<cavc::OffsetEndCapType, 3> endCapTypes = {
      cavc::OffsetEndCapType::Round, cavc::OffsetEndCapType::Square,
      cavc::OffsetEndCapType::Butt};

  for (std::uint32_t seed : seeds) {
    std::array<Pline, 2> inputs = {makeOpenMixedPolyline(seed), makeNearDegenerateOpenPolyline(seed)};
    for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex) {
      for (auto joinType : joinTypes) {
        for (auto endCapType : endCapTypes) {
          for (double offset : offsets) {
            cavc::ParallelOffsetOptions<double> options;
            options.joinType = joinType;
            options.endCapType = endCapType;
            options.miterLimit = 5.0;
            auto results = cavc::parallelOffset(inputs[inputIndex], offset, options);
            expectBasicOffsetInvariants(
                inputs[inputIndex], results,
                "open seed=" + std::to_string(seed) + " input=" + std::to_string(inputIndex));
          }
        }
      }
    }
  }
}

TEST(ParallelOffsetFuzzRegression, FixedSeedClosedJoinCorpusKeepsBasicInvariants) {
  constexpr std::array<std::uint32_t, 6> seeds = {23u, 37u, 59u, 89u, 131u, 197u};
  constexpr std::array<double, 4> offsets = {-0.5, -0.15, 0.15, 0.5};
  constexpr std::array<cavc::OffsetJoinType, 3> joinTypes = {
      cavc::OffsetJoinType::Round, cavc::OffsetJoinType::Miter, cavc::OffsetJoinType::Bevel};

  for (std::uint32_t seed : seeds) {
    std::array<Pline, 3> inputs = {makeConvexNgon(seed), makeConcaveStar(seed),
                                   makeArcHeavyClosedPolyline(seed)};
    for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex) {
      for (auto joinType : joinTypes) {
        for (double offset : offsets) {
          cavc::ParallelOffsetOptions<double> options;
          options.joinType = joinType;
          options.miterLimit = 5.0;
          auto results = cavc::parallelOffset(inputs[inputIndex], offset, options);
          expectBasicOffsetInvariants(
              inputs[inputIndex], results,
              "closed seed=" + std::to_string(seed) + " input=" + std::to_string(inputIndex));
        }
      }
    }
  }
}
