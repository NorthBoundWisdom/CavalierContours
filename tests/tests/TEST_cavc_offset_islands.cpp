#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "cavc/polylineoffsetislands.hpp"

namespace {
using cavc::createApproxSpatialIndex;
using cavc::invertDirection;
using cavc::kNoParentOffsetLoop;
using cavc::OffsetLoop;
using cavc::OffsetLoopRole;
using cavc::OffsetLoopSet;
using cavc::OffsetLoopTopologyNode;
using cavc::ParallelOffsetIslands;
using cavc::Polyline;
using cavc::getArea;
using cavc::getPathLength;

struct LoopProperties {
  std::size_t vertexCount;
  double area;
  double pathLength;
  double minX;
  double minY;
  double maxX;
  double maxY;
};

OffsetLoop<double> makeAxisAlignedRectLoop(double min_x, double min_y, double max_x, double max_y,
                                           bool clockwise, std::size_t parent_index) {
  Polyline<double> loop;
  loop.addVertex(min_x, min_y, 0.0);
  loop.addVertex(max_x, min_y, 0.0);
  loop.addVertex(max_x, max_y, 0.0);
  loop.addVertex(min_x, max_y, 0.0);
  loop.isClosed() = true;
  if (clockwise) {
    invertDirection(loop);
  }

  auto spatial_index = createApproxSpatialIndex(loop);
  return {parent_index, std::move(loop), std::move(spatial_index)};
}

OffsetLoop<double> makeLoop(std::initializer_list<std::array<double, 3>> vertexes,
                            std::size_t parent_index = 0) {
  Polyline<double> loop;
  for (auto const &vertex : vertexes) {
    loop.addVertex(vertex[0], vertex[1], vertex[2]);
  }
  loop.isClosed() = true;

  auto spatial_index = createApproxSpatialIndex(loop);
  return {parent_index, std::move(loop), std::move(spatial_index)};
}

LoopProperties describeLoop(OffsetLoop<double> const &loop) {
  return {loop.polyline.size(), getArea(loop.polyline), getPathLength(loop.polyline),
          loop.spatialIndex.minX(), loop.spatialIndex.minY(), loop.spatialIndex.maxX(),
          loop.spatialIndex.maxY()};
}

void appendLoopByOrientation(OffsetLoopSet<double> &loop_set, OffsetLoop<double> loop) {
  if (getArea(loop.polyline) < 0.0) {
    loop_set.cwLoops.push_back(std::move(loop));
  } else {
    loop_set.ccwLoops.push_back(std::move(loop));
  }
}

void expectLoopPropertiesNear(LoopProperties const &actual, LoopProperties const &expected,
                              double epsilon = 1e-5) {
  EXPECT_EQ(actual.vertexCount, expected.vertexCount);
  EXPECT_NEAR(actual.area, expected.area, epsilon);
  EXPECT_NEAR(actual.pathLength, expected.pathLength, epsilon);
  EXPECT_NEAR(actual.minX, expected.minX, epsilon);
  EXPECT_NEAR(actual.minY, expected.minY, epsilon);
  EXPECT_NEAR(actual.maxX, expected.maxX, epsilon);
  EXPECT_NEAR(actual.maxY, expected.maxY, epsilon);
}

std::size_t findNodeIndex(std::vector<OffsetLoopTopologyNode<double>> const &topology,
                          OffsetLoopRole role, std::size_t source_index) {
  for (std::size_t i = 0; i < topology.size(); ++i) {
    if (topology[i].role == role && topology[i].sourceIndex == source_index) {
      return i;
    }
  }
  return kNoParentOffsetLoop;
}
} // namespace

TEST(OffsetIslands, BuildTopologyProducesStableOrderAndParentChild) {
  OffsetLoopSet<double> loop_set;
  // Unsorted input on purpose.
  loop_set.ccwLoops.push_back(makeAxisAlignedRectLoop(12.0, 4.0, 16.0, 8.0, false, 0)); // island
  loop_set.ccwLoops.push_back(makeAxisAlignedRectLoop(0.0, 0.0, 20.0, 20.0, false, 1)); // outer
  loop_set.cwLoops.push_back(makeAxisAlignedRectLoop(10.0, 2.0, 18.0, 10.0, true, 0));  // hole B
  loop_set.cwLoops.push_back(makeAxisAlignedRectLoop(2.0, 2.0, 8.0, 8.0, true, 1));     // hole A
  loop_set.cwLoops.push_back(makeAxisAlignedRectLoop(13.0, 5.0, 15.0, 7.0, true, 2));   // hole C

  std::vector<OffsetLoopTopologyNode<double>> topology = cavc::buildOffsetLoopTopology(loop_set);
  ASSERT_EQ(topology.size(), 5u);

  // Stable order is area-desc then bbox key.
  EXPECT_EQ(topology[0].role, OffsetLoopRole::Outer);
  EXPECT_EQ(topology[0].sourceIndex, 1u);
  EXPECT_EQ(topology[1].role, OffsetLoopRole::Hole);
  EXPECT_EQ(topology[1].sourceIndex, 0u);
  EXPECT_EQ(topology[2].role, OffsetLoopRole::Hole);
  EXPECT_EQ(topology[2].sourceIndex, 1u);
  EXPECT_EQ(topology[3].role, OffsetLoopRole::Outer);
  EXPECT_EQ(topology[3].sourceIndex, 0u);
  EXPECT_EQ(topology[4].role, OffsetLoopRole::Hole);
  EXPECT_EQ(topology[4].sourceIndex, 2u);

  std::size_t outer_index = findNodeIndex(topology, OffsetLoopRole::Outer, 1);
  std::size_t hole_b_index = findNodeIndex(topology, OffsetLoopRole::Hole, 0);
  std::size_t hole_a_index = findNodeIndex(topology, OffsetLoopRole::Hole, 1);
  std::size_t island_index = findNodeIndex(topology, OffsetLoopRole::Outer, 0);
  std::size_t hole_c_index = findNodeIndex(topology, OffsetLoopRole::Hole, 2);

  ASSERT_NE(outer_index, kNoParentOffsetLoop);
  ASSERT_NE(hole_b_index, kNoParentOffsetLoop);
  ASSERT_NE(hole_a_index, kNoParentOffsetLoop);
  ASSERT_NE(island_index, kNoParentOffsetLoop);
  ASSERT_NE(hole_c_index, kNoParentOffsetLoop);

  EXPECT_EQ(topology[outer_index].parentIndex, kNoParentOffsetLoop);
  EXPECT_EQ(topology[hole_b_index].parentIndex, outer_index);
  EXPECT_EQ(topology[hole_a_index].parentIndex, outer_index);
  EXPECT_EQ(topology[island_index].parentIndex, hole_b_index);
  EXPECT_EQ(topology[hole_c_index].parentIndex, island_index);
}

TEST(OffsetIslands, ComputeSortsOutputLoopsByStableGeometryOrder) {
  OffsetLoopSet<double> input;
  // Small then large to verify sorting is applied in compute result.
  input.ccwLoops.push_back(makeAxisAlignedRectLoop(0.0, 0.0, 2.0, 2.0, false, 0));
  input.ccwLoops.push_back(makeAxisAlignedRectLoop(10.0, 0.0, 20.0, 8.0, false, 1));

  ParallelOffsetIslands<double> algorithm;
  OffsetLoopSet<double> result = algorithm.compute(input, 0.5);

  ASSERT_EQ(result.ccwLoops.size(), 2u);
  ASSERT_TRUE(result.cwLoops.empty());

  double area0 = std::abs(cavc::getArea(result.ccwLoops[0].polyline));
  double area1 = std::abs(cavc::getArea(result.ccwLoops[1].polyline));
  EXPECT_GT(area0, area1);

  std::vector<OffsetLoopTopologyNode<double>> topology = cavc::buildOffsetLoopTopology(result);
  ASSERT_EQ(topology.size(), 2u);
  EXPECT_EQ(topology[0].role, OffsetLoopRole::Outer);
  EXPECT_EQ(topology[1].role, OffsetLoopRole::Outer);
  EXPECT_EQ(topology[0].parentIndex, kNoParentOffsetLoop);
  EXPECT_EQ(topology[1].parentIndex, kNoParentOffsetLoop);
}

TEST(OffsetIslands, ComputeHandlesIssue66SliceValidationRegression) {
  OffsetLoopSet<double> input;
  appendLoopByOrientation(input, makeLoop({
      {511.25220437557994, 328.84948025435654, 0.0},
      {561.2119896118824, 328.84948025435654, 0.0},
      {561.2119896118824, 363.8703101013724, 0.0},
      {511.25220437557994, 363.8703101013724, 0.0},
  }));
  appendLoopByOrientation(input, makeLoop({
      {540.0335350561843, 343.6169427142472, -0.2382488851276809},
      {537.4421349268171, 345.12517844750175, -0.009889532389405053},
      {537.3232102367999, 345.3220639672001, 0.0},
      {535.3578079577983, 348.7262405716385, 0.0},
      {535.32462892643, 348.7834560296831, -0.011646639385887355},
      {535.2271073347746, 348.9562631479, -0.25910007835503845},
      {535.2805330874, 352.1843242602999, 0.0},
      {537.0000084336, 355.1625429223, 0.0},
      {543.6691202685, 343.6113023833, 0.0},
      {540.2257285374734, 343.6113053474554, 0.0},
      {540.16153451323, 343.6120105305873, 0.0},
  }));
  appendLoopByOrientation(input, makeLoop({
      {535.4816659760771, 346.2417647657877, -0.23822264248219718},
      {535.4722003945319, 343.2448614622984, -0.009951542143624231},
      {535.3601905012999, 343.0416179385001, 0.0},
      {533.3951100416035, 339.6379987413748, 0.0},
      {533.3623248097243, 339.5809609408148, -0.011710653864622287},
      {533.2619655102294, 339.4110081813505, -0.11747560444569163},
      {532.1675755117166, 338.3268271911126, -0.13757287100703242},
      {530.4835288827071, 337.8423335382348, 0.0},
      {530.438959126, 337.8420348466, 0.0},
      {527.0000084336, 337.8420348466, 0.0},
      {533.6691202683, 349.3932753857, 0.0},
      {535.3908302312102, 346.4111802353502, 0.0},
      {535.4223097645596, 346.3552449203639, 0.0},
  }));

  ParallelOffsetIslands<double> algorithm;
  OffsetLoopSet<double> result = algorithm.compute(input, 0.8);

  std::vector<LoopProperties> actual;
  actual.reserve(result.ccwLoops.size() + result.cwLoops.size());
  for (auto const &loop : result.ccwLoops) {
    actual.push_back(describeLoop(loop));
  }
  for (auto const &loop : result.cwLoops) {
    actual.push_back(describeLoop(loop));
  }

  ASSERT_EQ(actual.size(), 2u);
  std::sort(actual.begin(), actual.end(), [](LoopProperties const &lhs, LoopProperties const &rhs) {
    return std::abs(lhs.area) > std::abs(rhs.area);
  });

  std::vector<LoopProperties> expected = {
      {4, 1616.2241538207163, 163.56123016663673, 512.0522043755799, 329.64948025435655,
       560.4119896118824, 363.0703101013724},
      {28, -148.47469897242397, 61.018056828113345, 526.2000084335999, 337.04203484659996,
       544.4691202685001, 355.96254292230003},
  };

  expectLoopPropertiesNear(actual[0], expected[0]);
  EXPECT_EQ(actual[1].vertexCount, expected[1].vertexCount);
  EXPECT_NEAR(actual[1].area, expected[1].area, 1e-5);
  EXPECT_NEAR(actual[1].pathLength, expected[1].pathLength, 1e-5);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
