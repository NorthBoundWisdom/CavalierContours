#include "cavc/polylineoffsetislands.hpp"
#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <vector>

namespace {
using cavc::OffsetLoop;
using cavc::OffsetLoopSet;
using cavc::ParallelOffsetIslands;
using cavc::Polyline;
using cavc::StaticSpatialIndex;

Polyline<double> makeNgonLoop(std::size_t vertex_count, double radius, double center_x,
                              double center_y) {
  CAVC_ASSERT(vertex_count >= 3, "vertex_count must be >= 3");

  Polyline<double> loop;
  loop.isClosed() = true;
  for (std::size_t i = 0; i < vertex_count; ++i) {
    double angle =
        static_cast<double>(i) * cavc::utils::tau<double>() / static_cast<double>(vertex_count);
    loop.addVertex(center_x + radius * std::cos(angle), center_y + radius * std::sin(angle), 0.0);
  }
  return loop;
}

Polyline<double> makeAxisAlignedRect(double min_x, double min_y, double max_x, double max_y,
                                     bool clockwise) {
  Polyline<double> loop;
  loop.addVertex(min_x, min_y, 0.0);
  loop.addVertex(max_x, min_y, 0.0);
  loop.addVertex(max_x, max_y, 0.0);
  loop.addVertex(min_x, max_y, 0.0);
  loop.isClosed() = true;
  if (clockwise) {
    cavc::invertDirection(loop);
  }
  return loop;
}

OffsetLoop<double> makeOffsetLoop(Polyline<double> &&polyline, std::size_t parent_loop_index) {
  auto spatial_index = cavc::createApproxSpatialIndex(polyline);
  return {parent_loop_index, std::move(polyline), std::move(spatial_index)};
}

OffsetLoopSet<double> createGridOffsetLoopSet(std::size_t grid_size) {
  OffsetLoopSet<double> result;
  std::size_t loop_count = grid_size * grid_size;
  result.ccwLoops.reserve(loop_count);
  result.cwLoops.reserve(loop_count);

  double const pitch = 32.0;
  double const outer_size = 24.0;
  double const hole_margin = 5.0;

  for (std::size_t y = 0; y < grid_size; ++y) {
    for (std::size_t x = 0; x < grid_size; ++x) {
      std::size_t parent_index = y * grid_size + x;
      double min_x = static_cast<double>(x) * pitch;
      double min_y = static_cast<double>(y) * pitch;
      double max_x = min_x + outer_size;
      double max_y = min_y + outer_size;

      auto outer = makeAxisAlignedRect(min_x, min_y, max_x, max_y, false);
      auto hole = makeAxisAlignedRect(min_x + hole_margin, min_y + hole_margin, max_x - hole_margin,
                                      max_y - hole_margin, true);

      result.ccwLoops.push_back(makeOffsetLoop(std::move(outer), parent_index));
      result.cwLoops.push_back(makeOffsetLoop(std::move(hole), parent_index));
    }
  }

  return result;
}

struct BatchIndexSetup {
  std::vector<Polyline<double>> loops;

  BatchIndexSetup(std::size_t loop_count, std::size_t vertex_count) {
    loops.reserve(loop_count);
    std::size_t cols =
        static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(loop_count))));
    double const pitch = 96.0;

    for (std::size_t i = 0; i < loop_count; ++i) {
      std::size_t row = i / cols;
      std::size_t col = i % cols;
      double radius = 20.0 + static_cast<double>(i % 7);
      double center_x = static_cast<double>(col) * pitch;
      double center_y = static_cast<double>(row) * pitch;
      loops.push_back(makeNgonLoop(vertex_count, radius, center_x, center_y));
    }
  }
};

static void BM_createApproxSpatialIndicesScalar(benchmark::State &state) {
  std::size_t loop_count = static_cast<std::size_t>(state.range(0));
  std::size_t vertex_count = static_cast<std::size_t>(state.range(1));
  BatchIndexSetup setup(loop_count, vertex_count);

  state.counters["loopCount"] = static_cast<double>(loop_count);
  state.counters["vertexCount"] = static_cast<double>(vertex_count);

  for (auto _ : state) {
    (void)_;
    std::vector<StaticSpatialIndex<double>> indexes;
    indexes.reserve(setup.loops.size());
    for (auto const &loop : setup.loops) {
      indexes.push_back(cavc::createApproxSpatialIndex(loop));
    }
    benchmark::DoNotOptimize(indexes);
  }
}

static void BM_createApproxSpatialIndicesBatch(benchmark::State &state) {
  std::size_t loop_count = static_cast<std::size_t>(state.range(0));
  std::size_t vertex_count = static_cast<std::size_t>(state.range(1));
  BatchIndexSetup setup(loop_count, vertex_count);

  state.counters["loopCount"] = static_cast<double>(loop_count);
  state.counters["vertexCount"] = static_cast<double>(vertex_count);

  for (auto _ : state) {
    (void)_;
    auto indexes = cavc::createApproxSpatialIndices(setup.loops);
    benchmark::DoNotOptimize(indexes);
  }
}

static void BM_buildOffsetLoopTopologyLargeScale(benchmark::State &state) {
  std::size_t grid_size = static_cast<std::size_t>(state.range(0));
  OffsetLoopSet<double> input = createGridOffsetLoopSet(grid_size);

  state.counters["outerLoopCount"] = static_cast<double>(input.ccwLoops.size());
  state.counters["holeLoopCount"] = static_cast<double>(input.cwLoops.size());

  for (auto _ : state) {
    (void)_;
    auto topology = cavc::buildOffsetLoopTopology(input);
    benchmark::DoNotOptimize(topology);
  }
}

static void BM_parallelOffsetIslandsComputeLargeScale(benchmark::State &state) {
  std::size_t grid_size = static_cast<std::size_t>(state.range(0));
  OffsetLoopSet<double> input = createGridOffsetLoopSet(grid_size);
  ParallelOffsetIslands<double> algorithm;

  state.counters["outerLoopCount"] = static_cast<double>(input.ccwLoops.size());
  state.counters["holeLoopCount"] = static_cast<double>(input.cwLoops.size());

  for (auto _ : state) {
    (void)_;
    auto result = algorithm.compute(input, 1.0);
    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK(BM_createApproxSpatialIndicesScalar)
    ->Args({32, 16})
    ->Args({64, 32})
    ->Args({128, 64})
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_createApproxSpatialIndicesBatch)
    ->Args({32, 16})
    ->Args({64, 32})
    ->Args({128, 64})
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_buildOffsetLoopTopologyLargeScale)
    ->Arg(8)
    ->Arg(12)
    ->Arg(16)
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_parallelOffsetIslandsComputeLargeScale)
    ->Arg(4)
    ->Arg(8)
    ->Arg(12)
    ->Unit(benchmark::kMillisecond);
} // namespace

BENCHMARK_MAIN();
