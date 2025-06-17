#include <cmath>
#include <cstddef>

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

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
