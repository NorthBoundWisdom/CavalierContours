#include "cavc/polylineoffset.hpp"
#include <benchmark/benchmark.h>
#include <cmath>
#include <vector>

namespace {

cavc::Polyline<double> makeRect() {
  cavc::Polyline<double> p;
  p.isClosed() = true;
  p.addVertex(0.0, 0.0, 0.0);
  p.addVertex(20.0, 0.0, 0.0);
  p.addVertex(20.0, 8.0, 0.0);
  p.addVertex(0.0, 8.0, 0.0);
  return p;
}

cavc::Polyline<double> makeConcaveLineLoop() {
  cavc::Polyline<double> p;
  p.isClosed() = true;
  p.addVertex(0.0, 0.0, 0.0);
  p.addVertex(12.0, 0.0, 0.0);
  p.addVertex(12.0, 4.0, 0.0);
  p.addVertex(7.0, 4.0, 0.0);
  p.addVertex(7.0, 7.0, 0.0);
  p.addVertex(12.0, 7.0, 0.0);
  p.addVertex(12.0, 11.0, 0.0);
  p.addVertex(0.0, 11.0, 0.0);
  return p;
}

cavc::Polyline<double> makeArcHeavyLoop() {
  cavc::Polyline<double> p;
  p.isClosed() = true;
  p.addVertex(1.585, -0.404, 0.419);
  p.addVertex(1.856, -1.235, 0.125);
  p.addVertex(2.795, -2.210, 0.130);
  p.addVertex(4.185, -2.751, 0.773);
  p.addVertex(4.559, -1.233, 0.000);
  p.addVertex(5.090, -0.316, -0.312);
  return p;
}

cavc::Polyline<double> makeOpenArcPolyline() {
  cavc::Polyline<double> p;
  p.addVertex(0.0, 0.0, 0.25);
  p.addVertex(2.0, 0.4, 0.0);
  p.addVertex(4.5, -0.2, -0.3);
  p.addVertex(7.0, 0.0, 0.0);
  return p;
}

static void BM_roundJoinArcHeavyClosed(benchmark::State &state) {
  auto input = makeArcHeavyLoop();
  cavc::ParallelOffsetOptions<double> options;
  options.joinType = cavc::OffsetJoinType::Round;
  for (auto _ : state) {
    auto results = cavc::parallelOffset(input, 0.2, options);
    benchmark::DoNotOptimize(results);
  }
}

static void BM_miterJoinRectangleClosed(benchmark::State &state) {
  auto input = makeRect();
  cavc::ParallelOffsetOptions<double> options;
  options.joinType = cavc::OffsetJoinType::Miter;
  options.miterLimit = 5.0;
  for (auto _ : state) {
    auto results = cavc::parallelOffset(input, -1.0, options);
    benchmark::DoNotOptimize(results);
  }
}

static void BM_bevelJoinConcaveClosed(benchmark::State &state) {
  auto input = makeConcaveLineLoop();
  cavc::ParallelOffsetOptions<double> options;
  options.joinType = cavc::OffsetJoinType::Bevel;
  for (auto _ : state) {
    auto results = cavc::parallelOffset(input, -0.45, options);
    benchmark::DoNotOptimize(results);
  }
}

static void BM_squareEndCapOpenArc(benchmark::State &state) {
  auto input = makeOpenArcPolyline();
  cavc::ParallelOffsetOptions<double> options;
  options.endCapType = cavc::OffsetEndCapType::Square;
  for (auto _ : state) {
    auto results = cavc::parallelOffset(input, -0.25, options);
    benchmark::DoNotOptimize(results);
  }
}

static void BM_buttEndCapOpenArc(benchmark::State &state) {
  auto input = makeOpenArcPolyline();
  cavc::ParallelOffsetOptions<double> options;
  options.endCapType = cavc::OffsetEndCapType::Butt;
  for (auto _ : state) {
    auto results = cavc::parallelOffset(input, -0.25, options);
    benchmark::DoNotOptimize(results);
  }
}

} // namespace

BENCHMARK(BM_roundJoinArcHeavyClosed)->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_miterJoinRectangleClosed)->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_bevelJoinConcaveClosed)->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_squareEndCapOpenArc)->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_buttEndCapOpenArc)->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();
