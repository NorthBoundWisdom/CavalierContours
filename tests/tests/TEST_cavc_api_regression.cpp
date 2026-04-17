#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include <cavc/polyline.hpp>
#include <cavc/polylineintersects.hpp>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "c_api_include/cavaliercontours.h"
#include "c_api_test_helpers.hpp"
#include "cavc/intrlineseg2circle2.hpp"
#include "shape_offset_self_intersect_case.hpp"

namespace t = testing;

namespace {
struct PlineDeleter {
  void operator()(cavc_pline *pline) const {
    if (pline != nullptr) {
      cavc_pline_delete(pline);
    }
  }
};

struct PlineListDeleter {
  void operator()(cavc_pline_list *list) const {
    if (list != nullptr) {
      cavc_pline_list_delete(list);
    }
  }
};

struct TopologyDeleter {
  void operator()(cavc_offset_loop_topology *topology) const {
    if (topology != nullptr) {
      cavc_offset_loop_topology_delete(topology);
    }
  }
};

using PlinePtr = std::unique_ptr<cavc_pline, PlineDeleter>;
using PlineListPtr = std::unique_ptr<cavc_pline_list, PlineListDeleter>;
using SpatialIndexPtr = std::unique_ptr<cavc_spatial_index, void (*)(cavc_spatial_index *)>;
using TopologyPtr = std::unique_ptr<cavc_offset_loop_topology, TopologyDeleter>;

std::vector<cavc_vertex> readVertexes(cavc_pline const *pline) {
  uint32_t count = cavc_pline_vertex_count(pline);
  std::vector<cavc_vertex> result(count);
  if (count != 0) {
    cavc_pline_vertex_data(pline, result.data());
  }
  return result;
}

cavc::Polyline<cavc_real> toCppPolyline(cavc_pline const *pline) {
  cavc::Polyline<cavc_real> result;
  result.isClosed() = cavc_pline_is_closed(pline) != 0;

  const std::vector<cavc_vertex> vertexes = readVertexes(pline);
  result.vertexes().reserve(vertexes.size());
  for (cavc_vertex const &v : vertexes) {
    result.addVertex(v.x, v.y, v.bulge);
  }

  return result;
}

bool hasSelfIntersect(cavc_pline const *pline) {
  const cavc::Polyline<cavc_real> cpp_pline = toCppPolyline(pline);
  if (cpp_pline.size() < 2) {
    return false;
  }

  auto spatial_index = cavc::createApproxSpatialIndex(cpp_pline);
  std::vector<cavc::PlineIntersect<cavc_real>> intersects;
  cavc::allSelfIntersects(cpp_pline, intersects, spatial_index);
  return !intersects.empty();
}

std::vector<uint32_t> querySpatialIndexSorted(cavc_spatial_index const *spatial_index,
                                              cavc_real min_x, cavc_real min_y, cavc_real max_x,
                                              cavc_real max_y) {
  uint32_t count = cavc_spatial_index_query_count(spatial_index, min_x, min_y, max_x, max_y);
  std::vector<uint32_t> results(count);
  if (count != 0) {
    cavc_spatial_index_query(spatial_index, min_x, min_y, max_x, max_y, results.data());
    std::sort(results.begin(), results.end());
  }
  return results;
}

bool vertexListAllFinite(std::vector<cavc_vertex> const &vertexes) {
  return std::all_of(vertexes.begin(), vertexes.end(), [](cavc_vertex const &v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.bulge);
  });
}

bool closedVertexListsFuzzyEqual(std::vector<cavc_vertex> const &left,
                                 std::vector<cavc_vertex> const &right) {
  if (left.size() != right.size()) {
    return false;
  }
  if (left.empty()) {
    return true;
  }

  std::size_t const vertex_count = left.size();
  std::size_t start_index = 0;
  for (; start_index < vertex_count; ++start_index) {
    if (vertexesFuzzyEqual(left[0], right[start_index])) {
      break;
    }
  }
  if (start_index == vertex_count) {
    return false;
  }

  std::size_t index = start_index;
  for (std::size_t i = 0; i < vertex_count; ++i) {
    if (!vertexesFuzzyEqual(left[i], right[index])) {
      return false;
    }
    index = nextWrappingIndex(right, index);
  }
  return true;
}

std::vector<cavc_vertex> makeAxisAlignedRectLoopVertexes(cavc_real min_x, cavc_real min_y,
                                                         cavc_real max_x, cavc_real max_y,
                                                         bool clockwise) {
  std::vector<cavc_vertex> result = {
      {min_x, min_y, 0.0}, {max_x, min_y, 0.0}, {max_x, max_y, 0.0}, {min_x, max_y, 0.0}};
  if (clockwise) {
    reverseDirection(result);
  }
  return result;
}

std::vector<cavc_vertex> makeOffsetCase1Vertexes() {
  return {
      {0.0, 25.0, 1.0},
      {0.0, 0.0, 0.0},
      {2.0, 0.0, 1.0},
      {10.0, 0.0, -0.5},
      {8.0, 9.0, 0.374794619217547},
      {21.0, 0.0, 0.0},
      {23.0, 0.0, 1.0},
      {32.0, 0.0, -0.5},
      {28.0, 0.0, 0.5},
      {39.0, 21.0, 0.0},
      {28.0, 12.0, 0.5},
  };
}

std::vector<cavc_vertex> makeCollapsedCombVertexes() {
  return {
      {0.0, 0.0, 0.0},          {60.0, 0.0, 0.0},         {60.0, 20.0, 0.0},
      {59.0, 20.0, 0.0},        {59.0, 18.0, 0.0},        {57.0, 18.0, 0.0},
      {57.0, 20.0, 0.0},        {39.666666666667, 20.0, 0.0},
      {39.666666666667, 18.0, 0.0},                         {37.666666666667, 18.0, 0.0},
      {37.666666666667, 20.0, 0.0},                         {20.333333333333, 20.0, 0.0},
      {20.333333333333, 18.0, 0.0},                         {18.333333333333, 18.0, 0.0},
      {18.333333333333, 20.0, 0.0},                         {0.0, 20.0, 0.0},
  };
}

uint32_t findTopologyNodeIndex(cavc_offset_loop_topology const *topology,
                               cavc_offset_loop_role role, uint32_t source_index) {
  uint32_t node_count = cavc_offset_loop_topology_count(topology);
  for (uint32_t i = 0; i < node_count; ++i) {
    cavc_offset_loop_topology_node node = cavc_offset_loop_topology_get(topology, i);
    if (node.role == role && node.source_index == source_index) {
      return i;
    }
  }
  return CAVC_OFFSET_LOOP_NO_PARENT;
}
} // namespace

TEST(CApiRegression, SetIsClosedRoundTripDoesNotChangeVertices) {
  std::vector<cavc_vertex> vertexes = {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}};
  PlinePtr pline(plineFromVertexes(vertexes, false));

  EXPECT_EQ(cavc_pline_is_closed(pline.get()), 0);
  cavc_pline_set_is_closed(pline.get(), 1);
  EXPECT_EQ(cavc_pline_is_closed(pline.get()), 1);
  cavc_pline_set_is_closed(pline.get(), 0);
  EXPECT_EQ(cavc_pline_is_closed(pline.get()), 0);

  std::vector<cavc_vertex> actual = readVertexes(pline.get());
  EXPECT_THAT(actual, t::Pointwise(VertexEqual(), vertexes));
}

TEST(CApiRegression, NewWithNullVertexDataReservesCapacityAndStartsEmpty) {
  PlinePtr pline(cavc_pline_new(nullptr, 6, 1));

  ASSERT_NE(pline.get(), nullptr);
  EXPECT_EQ(cavc_pline_is_closed(pline.get()), 1);
  EXPECT_EQ(cavc_pline_vertex_count(pline.get()), 0u);
  EXPECT_GE(cavc_pline_capacity(pline.get()), 6u);

  cavc_pline_add_vertex(pline.get(), cavc_vertex{0.0, 0.0, 0.0});
  cavc_pline_add_vertex(pline.get(), cavc_vertex{2.0, 0.0, 0.0});
  EXPECT_EQ(cavc_pline_vertex_count(pline.get()), 2u);
}

TEST(CApiRegression, ClearLeavesCapacityAndAllowsReuse) {
  std::vector<cavc_vertex> vertexes = {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}};
  PlinePtr pline(plineFromVertexes(vertexes, false));

  cavc_pline_set_capacity(pline.get(), 16);
  const uint32_t capacity_before_clear = cavc_pline_capacity(pline.get());
  ASSERT_GE(capacity_before_clear, 16u);

  cavc_pline_clear(pline.get());
  EXPECT_EQ(cavc_pline_vertex_count(pline.get()), 0u);
  EXPECT_EQ(cavc_pline_capacity(pline.get()), capacity_before_clear);

  cavc_pline_add_vertex(pline.get(), cavc_vertex{5.0, 6.0, 0.0});
  cavc_pline_add_vertex(pline.get(), cavc_vertex{7.0, 8.0, 0.0});
  std::vector<cavc_vertex> expected = {{5.0, 6.0, 0.0}, {7.0, 8.0, 0.0}};
  std::vector<cavc_vertex> actual = readVertexes(pline.get());
  EXPECT_THAT(actual, t::Pointwise(VertexEqual(), expected));
}

TEST(CApiRegression, RemoveRangeWithZeroCountIsNoOp) {
  std::vector<cavc_vertex> vertexes = {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}};
  PlinePtr pline(plineFromVertexes(vertexes, false));

  cavc_pline_remove_range(pline.get(), 1, 0);
  std::vector<cavc_vertex> actual = readVertexes(pline.get());
  EXPECT_THAT(actual, t::Pointwise(VertexEqual(), vertexes));
}

TEST(CApiRegression, SetVertexDataNullClearsAndReserves) {
  std::vector<cavc_vertex> vertexes = {{1.0, 2.0, 0.1}, {3.0, 4.0, 0.2}, {5.0, 6.0, 0.3}};
  PlinePtr pline(plineFromVertexes(vertexes, false));

  cavc_pline_set_vertex_data(pline.get(), nullptr, 8);
  EXPECT_EQ(cavc_pline_vertex_count(pline.get()), 0u);
  EXPECT_GE(cavc_pline_capacity(pline.get()), 8u);
  EXPECT_EQ(cavc_pline_is_closed(pline.get()), 0);

  std::vector<cavc_vertex> sentinel(4, cavc_vertex{-1.0, -1.0, -1.0});
  std::vector<cavc_vertex> expected = sentinel;
  cavc_pline_vertex_data(pline.get(), sentinel.data());
  EXPECT_THAT(sentinel, t::Pointwise(VertexEqual(), expected));
}

TEST(CApiRegression, PlineListReleaseTransfersOwnership) {
  std::vector<cavc_vertex> square = {
      {0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}, {10.0, 5.0, 0.0}, {0.0, 5.0, 0.0}};
  PlinePtr pline(plineFromVertexes(square, true));

  cavc_pline_list *raw_results = nullptr;
  cavc_parallel_offset(pline.get(), -1.0, &raw_results, defaultParallelOffsetOptions());
  PlineListPtr results(raw_results);

  ASSERT_NE(results.get(), nullptr);
  ASSERT_EQ(cavc_pline_list_count(results.get()), 1u);

  cavc_pline *released_raw = cavc_pline_list_release(results.get(), 0);
  PlinePtr released(released_raw);
  EXPECT_NE(released.get(), nullptr);
  EXPECT_EQ(cavc_pline_list_count(results.get()), 0u);
  EXPECT_GT(cavc_pline_vertex_count(released.get()), 0u);

  cavc_pline_list_delete(results.release());
  EXPECT_GT(cavc_pline_vertex_count(released.get()), 0u);
  EXPECT_GT(cavc_get_path_length(released.get()), 0.0);
}

TEST(CApiRegression, ParallelOffsetDefaultOptionsUseRoundJoinAndRoundEndCap) {
  cavc_parallel_offset_options options = cavc_parallel_offset_default_options();
  EXPECT_EQ(options.may_have_self_intersects, 0);
  EXPECT_EQ(options.join_type, CAVC_OFFSET_JOIN_ROUND);
  EXPECT_EQ(options.end_cap_type, CAVC_OFFSET_END_CAP_ROUND);
  EXPECT_EQ(options.miter_limit, 4.0);
}

TEST(CApiRegression, ParallelOffsetStructuredOptionsCanEnableSelfIntersectHint) {
  std::vector<cavc_vertex> square = {
      {0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}, {10.0, 5.0, 0.0}, {0.0, 5.0, 0.0}};
  PlinePtr pline(plineFromVertexes(square, true));

  cavc_parallel_offset_options options = cavc_parallel_offset_default_options();
  options.may_have_self_intersects = 1;

  cavc_pline_list *raw_results = nullptr;
  cavc_parallel_offset(pline.get(), -1.0, &raw_results, options);
  PlineListPtr results(raw_results);

  ASSERT_NE(results.get(), nullptr);
  ASSERT_EQ(cavc_pline_list_count(results.get()), 1u);
  EXPECT_GT(cavc_get_area(cavc_pline_list_get(results.get(), 0)), 0.0);
}

TEST(CApiRegression, ParallelOffsetMiterJoinForLinePolyline) {
  std::vector<cavc_vertex> square = {
      {0.0, 0.0, 0.0}, {20.0, 0.0, 0.0}, {20.0, 10.0, 0.0}, {0.0, 10.0, 0.0}};
  PlinePtr pline(plineFromVertexes(square, true));

  cavc_parallel_offset_options options = cavc_parallel_offset_default_options();
  options.join_type = CAVC_OFFSET_JOIN_MITER;

  cavc_pline_list *raw_results = nullptr;
  cavc_parallel_offset(pline.get(), -2.0, &raw_results, options);
  PlineListPtr results(raw_results);

  ASSERT_NE(results.get(), nullptr);
  ASSERT_EQ(cavc_pline_list_count(results.get()), 1u);
  cavc_pline *offset = cavc_pline_list_get(results.get(), 0);

  EXPECT_EQ(cavc_pline_vertex_count(offset), 4u);
  std::vector<cavc_vertex> vertexes = readVertexes(offset);
  EXPECT_TRUE(std::all_of(vertexes.begin(), vertexes.end(),
                          [](cavc_vertex const &v) { return v.bulge == 0.0; }));

  EXPECT_NEAR(cavc_get_area(offset), 336.0, 1e-9);
  cavc_real min_x = 0.0;
  cavc_real min_y = 0.0;
  cavc_real max_x = 0.0;
  cavc_real max_y = 0.0;
  cavc_get_extents(offset, &min_x, &min_y, &max_x, &max_y);
  EXPECT_NEAR(min_x, -2.0, 1e-9);
  EXPECT_NEAR(min_y, -2.0, 1e-9);
  EXPECT_NEAR(max_x, 22.0, 1e-9);
  EXPECT_NEAR(max_y, 12.0, 1e-9);
}

TEST(CApiRegression, ParallelOffsetMiterLimitCanForceBevelEquivalentResult) {
  std::vector<cavc_vertex> square = {
      {0.0, 0.0, 0.0}, {20.0, 0.0, 0.0}, {20.0, 10.0, 0.0}, {0.0, 10.0, 0.0}};
  PlinePtr pline(plineFromVertexes(square, true));

  cavc_parallel_offset_options miter_options = cavc_parallel_offset_default_options();
  miter_options.join_type = CAVC_OFFSET_JOIN_MITER;
  miter_options.miter_limit = 1.0;

  cavc_parallel_offset_options bevel_options = cavc_parallel_offset_default_options();
  bevel_options.join_type = CAVC_OFFSET_JOIN_BEVEL;

  cavc_pline_list *raw_miter = nullptr;
  cavc_parallel_offset(pline.get(), -2.0, &raw_miter, miter_options);
  PlineListPtr miter_results(raw_miter);

  cavc_pline_list *raw_bevel = nullptr;
  cavc_parallel_offset(pline.get(), -2.0, &raw_bevel, bevel_options);
  PlineListPtr bevel_results(raw_bevel);

  ASSERT_EQ(cavc_pline_list_count(miter_results.get()), 1u);
  ASSERT_EQ(cavc_pline_list_count(bevel_results.get()), 1u);

  cavc_pline *miter_pline = cavc_pline_list_get(miter_results.get(), 0);
  cavc_pline *bevel_pline = cavc_pline_list_get(bevel_results.get(), 0);

  EXPECT_EQ(cavc_pline_vertex_count(miter_pline), 8u);
  EXPECT_EQ(cavc_pline_vertex_count(bevel_pline), 8u);

  std::vector<cavc_vertex> miter_vertexes = readVertexes(miter_pline);
  std::vector<cavc_vertex> bevel_vertexes = readVertexes(bevel_pline);
  EXPECT_THAT(miter_vertexes, t::Pointwise(VertexEqual(), bevel_vertexes));
  EXPECT_NEAR(cavc_get_area(miter_pline), 328.0, 1e-9);
}

TEST(CApiRegression, ParallelOffsetBevelJoinForLinePolyline) {
  std::vector<cavc_vertex> square = {
      {0.0, 0.0, 0.0}, {20.0, 0.0, 0.0}, {20.0, 10.0, 0.0}, {0.0, 10.0, 0.0}};
  PlinePtr pline(plineFromVertexes(square, true));

  cavc_parallel_offset_options options = cavc_parallel_offset_default_options();
  options.join_type = CAVC_OFFSET_JOIN_BEVEL;

  cavc_pline_list *raw_results = nullptr;
  cavc_parallel_offset(pline.get(), -2.0, &raw_results, options);
  PlineListPtr results(raw_results);

  ASSERT_NE(results.get(), nullptr);
  ASSERT_EQ(cavc_pline_list_count(results.get()), 1u);
  cavc_pline *offset = cavc_pline_list_get(results.get(), 0);

  EXPECT_EQ(cavc_pline_vertex_count(offset), 8u);
  std::vector<cavc_vertex> vertexes = readVertexes(offset);
  EXPECT_TRUE(std::all_of(vertexes.begin(), vertexes.end(),
                          [](cavc_vertex const &v) { return v.bulge == 0.0; }));

  EXPECT_NEAR(cavc_get_area(offset), 328.0, 1e-9);
  cavc_real min_x = 0.0;
  cavc_real min_y = 0.0;
  cavc_real max_x = 0.0;
  cavc_real max_y = 0.0;
  cavc_get_extents(offset, &min_x, &min_y, &max_x, &max_y);
  EXPECT_NEAR(min_x, -2.0, 1e-9);
  EXPECT_NEAR(min_y, -2.0, 1e-9);
  EXPECT_NEAR(max_x, 22.0, 1e-9);
  EXPECT_NEAR(max_y, 12.0, 1e-9);
}

TEST(CApiRegression, ParallelOffsetMiterJoinSupportsArcInvolvedJoins) {
  std::vector<cavc_vertex> arc_join_shape = {
      {0.0, 0.0, 0.0}, {4.0, 0.0, 0.5}, {8.0, 0.0, 0.0}, {8.0, 4.0, 0.0}, {0.0, 4.0, 0.0}};
  PlinePtr pline(plineFromVertexes(arc_join_shape, true));

  cavc_parallel_offset_options options = cavc_parallel_offset_default_options();
  options.join_type = CAVC_OFFSET_JOIN_MITER;
  options.miter_limit = 8.0;

  cavc_pline_list *raw_results = nullptr;
  cavc_parallel_offset(pline.get(), -0.5, &raw_results, options);
  PlineListPtr results(raw_results);

  ASSERT_NE(results.get(), nullptr);
  ASSERT_EQ(cavc_pline_list_count(results.get()), 1u);
  cavc_pline *offset = cavc_pline_list_get(results.get(), 0);

  std::vector<cavc_vertex> vertexes = readVertexes(offset);
  ASSERT_FALSE(vertexes.empty());

  bool has_arc_segment = false;
  for (cavc_vertex const &v : vertexes) {
    EXPECT_TRUE(std::isfinite(v.x));
    EXPECT_TRUE(std::isfinite(v.y));
    EXPECT_TRUE(std::isfinite(v.bulge));
    has_arc_segment = has_arc_segment || std::abs(v.bulge) > 1e-12;
  }
  EXPECT_TRUE(has_arc_segment);
}

TEST(CApiRegression, ArcArcOpposingDirectionTouchAtEndsReturnsSharedEndpoint) {
  using Vertex = cavc::PlineVertex<cavc_real>;
  using Point = cavc::Vector2<cavc_real>;

  auto expect_single_intersection = [](Vertex const &a1, Vertex const &a2, Vertex const &b1,
                                       Vertex const &b2, Point const &expected) {
    auto const result = cavc::intrPlineSegs(a1, a2, b1, b2);
    ASSERT_EQ(result.intrType, cavc::PlineSegIntrType::OneIntersect);
    EXPECT_TRUE(cavc::fuzzyEqual(result.point1, expected, cavc_real(1e-9)));
  };

  Vertex const v1(-189.0, -196.91384910249, 0.553407781718062);
  Vertex const v2(-170.999999999999, -225.631646989572, -0.553407781718061);
  Vertex const u1(-153.0, -196.91384910249, -0.553407781718061);
  Vertex const u2(-171.0, -225.631646989571, -0.553407781718061);
  Point const expected(-171.0, -225.631646989571);

  expect_single_intersection(v1, v2, u1, u2, expected);
  expect_single_intersection(u1, u2, v1, v2, expected);

  Vertex const reversed_u1(-171.0, -225.631646989571, 0.553407781718062);
  Vertex const reversed_u2(-153.0, -196.91384910249, -0.553407781718061);
  expect_single_intersection(v1, v2, reversed_u1, reversed_u2, expected);
}

TEST(CApiRegression, LineCircleTangentAtStartPointIsStable) {
  auto result = cavc::intrLineSeg2Circle2(
      cavc::Vector2<cavc_real>(161.28999999999999, 113.66500000000001),
      cavc::Vector2<cavc_real>(167.63999999999999, 113.66500000000001), 0.63499999999999801,
      cavc::Vector2<cavc_real>(161.28999999999999, 114.30000000000001));

  ASSERT_EQ(result.numIntersects, 1);
  EXPECT_NEAR(result.t0, 0.0, 1e-8);
}

TEST(CApiRegression, OverlappingPillEndsBooleanRegression) {
  std::vector<cavc_vertex> pline_a = {
      {113.1450199999994, 99.04090098302, 0.0},
      {113.1449999999994, 114.30000098302, 1.0},
      {111.6450000000006, 114.29999901698, 0.0},
      {111.6450200000006, 99.04089901698, 1.0},
  };
  std::vector<cavc_vertex> pline_b = {
      {113.145, 114.3, 0.0},
      {113.145, 117.475, 1.0},
      {111.645, 117.475, 0.0},
      {111.645, 114.3, 1.0},
  };

  PlinePtr a(plineFromVertexes(pline_a, true));
  PlinePtr b(plineFromVertexes(pline_b, true));

  auto expect_remaining = [&](int mode, std::vector<PolylineProperties> const &expected) {
    cavc_pline_list *remaining = nullptr;
    cavc_pline_list *subtracted = nullptr;
    cavc_combine_plines(a.get(), b.get(), mode, &remaining, &subtracted);
    ASSERT_EQ(cavc_pline_list_count(remaining), expected.size());
    ASSERT_EQ(cavc_pline_list_count(subtracted), 0u);

    std::vector<PolylineProperties> actual;
    actual.reserve(expected.size());
    for (uint32_t i = 0; i < expected.size(); ++i) {
      actual.emplace_back(cavc_pline_list_get(remaining, i));
    }

    ASSERT_THAT(actual, t::UnorderedPointwise(EqIgnoreSignOfArea(), expected));
    cavc_pline_list_delete(remaining);
    cavc_pline_list_delete(subtracted);
  };

  expect_remaining(0, {PolylineProperties(6, 29.418295867659253, 41.58058898041103, 111.645,
                                          98.29089999999995, 113.14502000000005, 118.225)});
  expect_remaining(2, {PolylineProperties(2, 1.7671458676433858, 4.712388980383754, 111.645,
                                          113.54999950849015, 113.14500000000015, 115.04999950848986)});
  expect_remaining(1, {PolylineProperties(4, -22.88864926275236, 35.230587997390565, 111.6450000000006,
                                          98.29089999999995, 113.14502000000005, 114.3)});
  expect_remaining(3, {PolylineProperties(4, -22.88864926275236, 35.230587997390565, 111.6450000000006,
                                          98.29089999999995, 113.14502000000005, 114.3),
                       PolylineProperties(4, 4.762500737263508, 11.062389963404218, 111.645,
                                          114.29999901698, 113.145, 118.225)});
}

TEST(CApiRegression, ParallelOffsetArcInvolvedMiterLimitCanForceBevelEquivalentResult) {
  std::vector<cavc_vertex> arc_join_shape = {
      {0.0, 0.0, 0.0}, {4.0, 0.0, 0.5}, {8.0, 0.0, 0.0}, {8.0, 4.0, 0.0}, {0.0, 4.0, 0.0}};
  PlinePtr pline(plineFromVertexes(arc_join_shape, true));

  cavc_parallel_offset_options miter_options = cavc_parallel_offset_default_options();
  miter_options.join_type = CAVC_OFFSET_JOIN_MITER;
  miter_options.miter_limit = 1.0;

  cavc_parallel_offset_options bevel_options = cavc_parallel_offset_default_options();
  bevel_options.join_type = CAVC_OFFSET_JOIN_BEVEL;

  cavc_pline_list *raw_miter = nullptr;
  cavc_parallel_offset(pline.get(), -0.5, &raw_miter, miter_options);
  PlineListPtr miter_results(raw_miter);

  cavc_pline_list *raw_bevel = nullptr;
  cavc_parallel_offset(pline.get(), -0.5, &raw_bevel, bevel_options);
  PlineListPtr bevel_results(raw_bevel);

  ASSERT_EQ(cavc_pline_list_count(miter_results.get()), 1u);
  ASSERT_EQ(cavc_pline_list_count(bevel_results.get()), 1u);

  std::vector<cavc_vertex> miter_vertexes =
      readVertexes(cavc_pline_list_get(miter_results.get(), 0));
  std::vector<cavc_vertex> bevel_vertexes =
      readVertexes(cavc_pline_list_get(bevel_results.get(), 0));
  ASSERT_FALSE(miter_vertexes.empty());
  ASSERT_EQ(miter_vertexes.size(), bevel_vertexes.size());

  std::vector<std::vector<cavc_vertex>> miter_vertex_sets = {miter_vertexes};
  std::vector<std::vector<cavc_vertex>> bevel_vertex_sets = {bevel_vertexes};
  EXPECT_THAT(miter_vertex_sets, t::Pointwise(VertexListsFuzzyEqual(true), bevel_vertex_sets));
}

TEST(CApiRegression, ParallelOffsetJoinArcHeavyStressMatrix) {
  struct JoinMatrixCase {
    const char *name;
    std::vector<cavc_vertex> vertexes;
    uint32_t min_bevel_extra_positive;
    uint32_t min_bevel_extra_negative;
  };

  std::vector<JoinMatrixCase> cases;
  cases.push_back(JoinMatrixCase{
      "high_curvature_arc_heavy",
      {
          {1.585, -0.404, 0.419},
          {1.856, -1.235, 0.125},
          {2.795, -2.210, 0.130},
          {4.185, -2.751, 0.773},
          {4.559, -1.233, 0.000},
          {5.090, -0.316, -0.312},
      },
      2,
      2,
  });
  cases.push_back(JoinMatrixCase{
      "near_tangent_turns",
      {
          {0.277, -0.934, 0.000},
          {0.616, -0.239, 0.125},
          {1.104, 0.666, 0.000},
          {2.029, 2.174, 0.000},
          {3.090, 0.843, 0.409},
          {3.904, 2.167, 0.000},
          {3.944, 3.010, -0.517},
      },
      2,
      1,
  });
  cases.push_back(JoinMatrixCase{
      "degenerate_repeated_vertex",
      {
          {1.092, 0.310, -0.847},
          {1.698, -0.231, 0.102},
          {1.698, -0.231, -0.400},
          {2.288, -0.255, 0.000},
          {3.848, 0.494, 0.515},
          {4.606, 0.328, 0.000},
          {4.888, 1.349, 0.374},
      },
      2,
      2,
  });

  auto run_case = [](JoinMatrixCase const &test_case, cavc_real delta, uint32_t min_bevel_extra) {
    PlinePtr pline(plineFromVertexes(test_case.vertexes, true));

    cavc_parallel_offset_options miter_high_options = cavc_parallel_offset_default_options();
    miter_high_options.join_type = CAVC_OFFSET_JOIN_MITER;
    miter_high_options.miter_limit = 8.0;

    cavc_parallel_offset_options miter_low_options = cavc_parallel_offset_default_options();
    miter_low_options.join_type = CAVC_OFFSET_JOIN_MITER;
    miter_low_options.miter_limit = 1.0;

    cavc_parallel_offset_options bevel_options = cavc_parallel_offset_default_options();
    bevel_options.join_type = CAVC_OFFSET_JOIN_BEVEL;

    cavc_parallel_offset_options round_options = cavc_parallel_offset_default_options();

    cavc_pline_list *raw_miter_high = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_miter_high, miter_high_options);
    PlineListPtr miter_high_results(raw_miter_high);

    cavc_pline_list *raw_miter_low = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_miter_low, miter_low_options);
    PlineListPtr miter_low_results(raw_miter_low);

    cavc_pline_list *raw_bevel = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_bevel, bevel_options);
    PlineListPtr bevel_results(raw_bevel);

    cavc_pline_list *raw_round = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_round, round_options);
    PlineListPtr round_results(raw_round);

    ASSERT_EQ(cavc_pline_list_count(miter_high_results.get()), 1u);
    ASSERT_EQ(cavc_pline_list_count(miter_low_results.get()), 1u);
    ASSERT_EQ(cavc_pline_list_count(bevel_results.get()), 1u);
    ASSERT_EQ(cavc_pline_list_count(round_results.get()), 1u);

    cavc_pline *miter_high_pline = cavc_pline_list_get(miter_high_results.get(), 0);
    cavc_pline *miter_low_pline = cavc_pline_list_get(miter_low_results.get(), 0);
    cavc_pline *bevel_pline = cavc_pline_list_get(bevel_results.get(), 0);
    cavc_pline *round_pline = cavc_pline_list_get(round_results.get(), 0);

    std::vector<cavc_vertex> miter_high_vertexes = readVertexes(miter_high_pline);
    std::vector<cavc_vertex> miter_low_vertexes = readVertexes(miter_low_pline);
    std::vector<cavc_vertex> bevel_vertexes = readVertexes(bevel_pline);
    std::vector<cavc_vertex> round_vertexes = readVertexes(round_pline);

    ASSERT_FALSE(miter_high_vertexes.empty());
    ASSERT_FALSE(miter_low_vertexes.empty());
    ASSERT_FALSE(bevel_vertexes.empty());
    ASSERT_FALSE(round_vertexes.empty());

    std::vector<std::vector<cavc_vertex>> miter_low_wrap = {miter_low_vertexes};
    std::vector<std::vector<cavc_vertex>> bevel_wrap = {bevel_vertexes};
    EXPECT_THAT(miter_low_wrap, t::Pointwise(VertexListsFuzzyEqual(true), bevel_wrap));

    const bool high_equals_bevel = closedVertexListsFuzzyEqual(miter_high_vertexes, bevel_vertexes);
    const bool high_equals_round = closedVertexListsFuzzyEqual(miter_high_vertexes, round_vertexes);
    EXPECT_TRUE(high_equals_round || !high_equals_bevel);
    if (!high_equals_bevel && !high_equals_round) {
      EXPECT_GE(bevel_vertexes.size(), miter_high_vertexes.size() + min_bevel_extra);
    }

    EXPECT_TRUE(vertexListAllFinite(miter_high_vertexes));
    EXPECT_TRUE(vertexListAllFinite(miter_low_vertexes));
    EXPECT_TRUE(vertexListAllFinite(bevel_vertexes));
    EXPECT_TRUE(vertexListAllFinite(round_vertexes));
    EXPECT_GT(cavc_get_path_length(miter_high_pline), 0.0);
    EXPECT_GT(cavc_get_path_length(miter_low_pline), 0.0);
    EXPECT_GT(cavc_get_path_length(bevel_pline), 0.0);
    EXPECT_GT(cavc_get_path_length(round_pline), 0.0);
  };

  for (auto const &test_case : cases) {
    SCOPED_TRACE(test_case.name);
    run_case(test_case, 0.2, test_case.min_bevel_extra_positive);
    run_case(test_case, -0.2, test_case.min_bevel_extra_negative);
  }
}

TEST(CApiRegression, ParallelOffsetOffsetCase1JoinModesKeepClosedOutputs) {
  const std::vector<cavc_vertex> offset_case1 = makeOffsetCase1Vertexes();

  auto run_case = [&](cavc_offset_join_type join_type, cavc_real delta) {
    PlinePtr pline(plineFromVertexes(offset_case1, true));
    cavc_parallel_offset_options options = cavc_parallel_offset_default_options();
    options.join_type = join_type;

    cavc_pline_list *raw_results = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_results, options);
    PlineListPtr results(raw_results);

    ASSERT_NE(results.get(), nullptr);
    ASSERT_GT(cavc_pline_list_count(results.get()), 0u);

    const uint32_t result_count = cavc_pline_list_count(results.get());
    for (uint32_t i = 0; i < result_count; ++i) {
      cavc_pline *offset = cavc_pline_list_get(results.get(), i);
      ASSERT_NE(offset, nullptr);

      EXPECT_EQ(cavc_pline_is_closed(offset), 1)
          << "join=" << join_type << ", delta=" << delta << ", index=" << i;
      EXPECT_GE(cavc_pline_vertex_count(offset), 3u);
      EXPECT_GT(cavc_get_path_length(offset), 0.0);

      const std::vector<cavc_vertex> vertexes = readVertexes(offset);
      EXPECT_TRUE(vertexListAllFinite(vertexes));
    }
  };

  run_case(CAVC_OFFSET_JOIN_ROUND, 1.0);
  run_case(CAVC_OFFSET_JOIN_ROUND, -1.0);
  run_case(CAVC_OFFSET_JOIN_MITER, 1.0);
  run_case(CAVC_OFFSET_JOIN_MITER, -1.0);
  run_case(CAVC_OFFSET_JOIN_BEVEL, 1.0);
  run_case(CAVC_OFFSET_JOIN_BEVEL, -1.0);
}

TEST(CApiRegression, ParallelOffsetOffsetCase1OpenMiterRoundEndCapNegativeRegression) {
  const std::vector<cavc_vertex> offset_case1 = makeOffsetCase1Vertexes();
  PlinePtr pline(plineFromVertexes(offset_case1, false));

  cavc_parallel_offset_options options = cavc_parallel_offset_default_options();
  options.join_type = CAVC_OFFSET_JOIN_MITER;
  options.end_cap_type = CAVC_OFFSET_END_CAP_ROUND;
  options.miter_limit = 4.0;

  cavc_pline_list *raw_results = nullptr;
  cavc_parallel_offset(pline.get(), -1.55, &raw_results, options);
  PlineListPtr results(raw_results);

  ASSERT_NE(results.get(), nullptr);
  ASSERT_EQ(cavc_pline_list_count(results.get()), 1u);

  cavc_pline *offset = cavc_pline_list_get(results.get(), 0);
  ASSERT_NE(offset, nullptr);
  EXPECT_EQ(cavc_pline_is_closed(offset), 0);
  EXPECT_GE(cavc_pline_vertex_count(offset), 3u);
  EXPECT_GT(cavc_get_path_length(offset), 0.0);

  const std::vector<cavc_vertex> vertexes = readVertexes(offset);
  ASSERT_FALSE(vertexes.empty());
  EXPECT_TRUE(vertexListAllFinite(vertexes));

  cavc_real min_x = 0.0;
  cavc_real min_y = 0.0;
  cavc_real max_x = 0.0;
  cavc_real max_y = 0.0;
  cavc_get_extents(offset, &min_x, &min_y, &max_x, &max_y);
  EXPECT_TRUE(std::isfinite(min_x));
  EXPECT_TRUE(std::isfinite(min_y));
  EXPECT_TRUE(std::isfinite(max_x));
  EXPECT_TRUE(std::isfinite(max_y));
  EXPECT_LT(min_x, max_x);
  EXPECT_LT(min_y, max_y);
}

TEST(CApiRegression, ParallelOffsetOffsetCase1ClosedJoinMatrixProducesSimpleLoops) {
  const std::vector<cavc_vertex> offset_case1 = makeOffsetCase1Vertexes();

  const std::array<cavc_real, 6> deltas = {0.2, 0.6, 1.0, -0.2, -0.6, -1.0};
  const std::array<cavc_offset_join_type, 3> join_types = {
      CAVC_OFFSET_JOIN_ROUND, CAVC_OFFSET_JOIN_MITER, CAVC_OFFSET_JOIN_BEVEL};

  for (cavc_offset_join_type join_type : join_types) {
    for (cavc_real delta : deltas) {
      SCOPED_TRACE(::testing::Message() << "join=" << join_type << ", delta=" << delta);

      PlinePtr pline(plineFromVertexes(offset_case1, true));
      cavc_parallel_offset_options options = cavc_parallel_offset_default_options();
      options.join_type = join_type;

      cavc_pline_list *raw_results = nullptr;
      cavc_parallel_offset(pline.get(), delta, &raw_results, options);
      PlineListPtr results(raw_results);

      ASSERT_NE(results.get(), nullptr);
      ASSERT_GT(cavc_pline_list_count(results.get()), 0u);

      const uint32_t result_count = cavc_pline_list_count(results.get());
      for (uint32_t i = 0; i < result_count; ++i) {
        cavc_pline *offset = cavc_pline_list_get(results.get(), i);
        ASSERT_NE(offset, nullptr);
        EXPECT_EQ(cavc_pline_is_closed(offset), 1);
        EXPECT_GE(cavc_pline_vertex_count(offset), 3u);
        EXPECT_GT(cavc_get_path_length(offset), 0.0);
        EXPECT_FALSE(hasSelfIntersect(offset));
        EXPECT_TRUE(vertexListAllFinite(readVertexes(offset)));
      }
    }
  }
}

TEST(CApiRegression, ParallelOffsetOffsetCase1OpenJoinEndCapMatrixProducesSimplePolylines) {
  const std::vector<cavc_vertex> offset_case1 = makeOffsetCase1Vertexes();

  const std::array<cavc_real, 4> deltas = {0.2, 0.6, -0.6, -1.55};
  const std::array<cavc_offset_join_type, 3> join_types = {
      CAVC_OFFSET_JOIN_ROUND, CAVC_OFFSET_JOIN_MITER, CAVC_OFFSET_JOIN_BEVEL};
  const std::array<cavc_offset_end_cap_type, 3> end_cap_types = {
      CAVC_OFFSET_END_CAP_ROUND, CAVC_OFFSET_END_CAP_SQUARE, CAVC_OFFSET_END_CAP_BUTT};

  for (cavc_offset_join_type join_type : join_types) {
    for (cavc_offset_end_cap_type end_cap_type : end_cap_types) {
      for (cavc_real delta : deltas) {
        SCOPED_TRACE(::testing::Message()
                     << "join=" << join_type << ", endcap=" << end_cap_type << ", delta=" << delta);

        PlinePtr pline(plineFromVertexes(offset_case1, false));
        cavc_parallel_offset_options options = cavc_parallel_offset_default_options();
        options.join_type = join_type;
        options.end_cap_type = end_cap_type;

        cavc_pline_list *raw_results = nullptr;
        cavc_parallel_offset(pline.get(), delta, &raw_results, options);
        PlineListPtr results(raw_results);

        ASSERT_NE(results.get(), nullptr);
        ASSERT_GT(cavc_pline_list_count(results.get()), 0u);

        const uint32_t result_count = cavc_pline_list_count(results.get());
        for (uint32_t i = 0; i < result_count; ++i) {
          cavc_pline *offset = cavc_pline_list_get(results.get(), i);
          ASSERT_NE(offset, nullptr);
          EXPECT_EQ(cavc_pline_is_closed(offset), 0);
          EXPECT_GE(cavc_pline_vertex_count(offset), 2u);
          EXPECT_GT(cavc_get_path_length(offset), 0.0);
          EXPECT_FALSE(hasSelfIntersect(offset));
          EXPECT_TRUE(vertexListAllFinite(readVertexes(offset)));
        }
      }
    }
  }
}

TEST(CApiRegression, ParallelOffsetSquareEndCapExtendsOpenLineEndpoints) {
  std::vector<cavc_vertex> line = {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}};
  PlinePtr pline(plineFromVertexes(line, false));

  cavc_parallel_offset_options options = cavc_parallel_offset_default_options();
  options.end_cap_type = CAVC_OFFSET_END_CAP_SQUARE;

  cavc_pline_list *raw_results = nullptr;
  cavc_parallel_offset(pline.get(), -1.0, &raw_results, options);
  PlineListPtr results(raw_results);

  ASSERT_NE(results.get(), nullptr);
  ASSERT_EQ(cavc_pline_list_count(results.get()), 1u);
  cavc_pline *offset = cavc_pline_list_get(results.get(), 0);

  EXPECT_EQ(cavc_pline_vertex_count(offset), 2u);
  std::vector<cavc_vertex> vertexes = readVertexes(offset);
  ASSERT_EQ(vertexes.size(), 2u);
  EXPECT_NEAR(vertexes[0].x, -1.0, 1e-9);
  EXPECT_NEAR(vertexes[0].y, -1.0, 1e-9);
  EXPECT_NEAR(vertexes[0].bulge, 0.0, 1e-12);
  EXPECT_NEAR(vertexes[1].x, 11.0, 1e-9);
  EXPECT_NEAR(vertexes[1].y, -1.0, 1e-9);
  EXPECT_NEAR(vertexes[1].bulge, 0.0, 1e-12);

  EXPECT_NEAR(cavc_get_path_length(offset), 12.0, 1e-9);
}

TEST(CApiRegression, ParallelOffsetSquareEndCapSupportsArcStartSegment) {
  std::vector<cavc_vertex> arc_start_open = {{0.0, 0.0, 1.0}, {2.0, 0.0, 0.0}, {4.0, 0.0, 0.0}};
  PlinePtr pline(plineFromVertexes(arc_start_open, false));

  cavc_parallel_offset_options round_options = cavc_parallel_offset_default_options();
  cavc_parallel_offset_options square_options = cavc_parallel_offset_default_options();
  square_options.end_cap_type = CAVC_OFFSET_END_CAP_SQUARE;

  cavc_pline_list *raw_round = nullptr;
  cavc_parallel_offset(pline.get(), -0.25, &raw_round, round_options);
  PlineListPtr round_results(raw_round);

  cavc_pline_list *raw_square = nullptr;
  cavc_parallel_offset(pline.get(), -0.25, &raw_square, square_options);
  PlineListPtr square_results(raw_square);

  ASSERT_EQ(cavc_pline_list_count(round_results.get()), 1u);
  ASSERT_EQ(cavc_pline_list_count(square_results.get()), 1u);

  cavc_pline *round_offset = cavc_pline_list_get(round_results.get(), 0);
  cavc_pline *square_offset = cavc_pline_list_get(square_results.get(), 0);
  std::vector<cavc_vertex> round_vertexes = readVertexes(round_offset);
  std::vector<cavc_vertex> square_vertexes = readVertexes(square_offset);
  ASSERT_FALSE(round_vertexes.empty());
  ASSERT_FALSE(square_vertexes.empty());

  auto endpoint_delta = [](cavc_vertex const &a, cavc_vertex const &b) {
    return std::hypot(a.x - b.x, a.y - b.y);
  };

  EXPECT_NEAR(endpoint_delta(square_vertexes.front(), round_vertexes.front()), 0.25, 1e-9);
  EXPECT_NEAR(endpoint_delta(square_vertexes.back(), round_vertexes.back()), 0.25, 1e-9);
  EXPECT_NEAR(cavc_get_path_length(square_offset) - cavc_get_path_length(round_offset), 0.5, 1e-8);
}

TEST(CApiRegression, ParallelOffsetSquareEndCapSupportsArcStartSegmentPositiveDelta) {
  std::vector<cavc_vertex> arc_start_open = {{0.0, 0.0, 1.0}, {2.0, 0.0, 0.0}, {4.0, 0.0, 0.0}};
  PlinePtr pline(plineFromVertexes(arc_start_open, false));

  cavc_parallel_offset_options round_options = cavc_parallel_offset_default_options();
  cavc_parallel_offset_options square_options = cavc_parallel_offset_default_options();
  square_options.end_cap_type = CAVC_OFFSET_END_CAP_SQUARE;

  cavc_pline_list *raw_round = nullptr;
  cavc_parallel_offset(pline.get(), 0.25, &raw_round, round_options);
  PlineListPtr round_results(raw_round);

  cavc_pline_list *raw_square = nullptr;
  cavc_parallel_offset(pline.get(), 0.25, &raw_square, square_options);
  PlineListPtr square_results(raw_square);

  ASSERT_EQ(cavc_pline_list_count(round_results.get()), 1u);
  ASSERT_EQ(cavc_pline_list_count(square_results.get()), 1u);

  cavc_pline *round_offset = cavc_pline_list_get(round_results.get(), 0);
  cavc_pline *square_offset = cavc_pline_list_get(square_results.get(), 0);
  std::vector<cavc_vertex> round_vertexes = readVertexes(round_offset);
  std::vector<cavc_vertex> square_vertexes = readVertexes(square_offset);
  ASSERT_FALSE(round_vertexes.empty());
  ASSERT_FALSE(square_vertexes.empty());

  auto endpoint_delta = [](cavc_vertex const &a, cavc_vertex const &b) {
    return std::hypot(a.x - b.x, a.y - b.y);
  };

  EXPECT_NEAR(endpoint_delta(square_vertexes.front(), round_vertexes.front()), 0.25, 1e-9);
  EXPECT_NEAR(endpoint_delta(square_vertexes.back(), round_vertexes.back()), 0.25, 1e-9);
  EXPECT_NEAR(cavc_get_path_length(square_offset) - cavc_get_path_length(round_offset), 0.5, 1e-8);
}

TEST(CApiRegression, ParallelOffsetSquareEndCapSupportsArcEndSegment) {
  std::vector<cavc_vertex> arc_end_open = {{0.0, 0.0, 0.0}, {2.0, 0.0, 1.0}, {4.0, 0.0, 0.0}};
  PlinePtr pline(plineFromVertexes(arc_end_open, false));

  cavc_parallel_offset_options round_options = cavc_parallel_offset_default_options();
  cavc_parallel_offset_options square_options = cavc_parallel_offset_default_options();
  square_options.end_cap_type = CAVC_OFFSET_END_CAP_SQUARE;

  cavc_pline_list *raw_round = nullptr;
  cavc_parallel_offset(pline.get(), -0.25, &raw_round, round_options);
  PlineListPtr round_results(raw_round);

  cavc_pline_list *raw_square = nullptr;
  cavc_parallel_offset(pline.get(), -0.25, &raw_square, square_options);
  PlineListPtr square_results(raw_square);

  ASSERT_EQ(cavc_pline_list_count(round_results.get()), 1u);
  ASSERT_EQ(cavc_pline_list_count(square_results.get()), 1u);

  cavc_pline *round_offset = cavc_pline_list_get(round_results.get(), 0);
  cavc_pline *square_offset = cavc_pline_list_get(square_results.get(), 0);
  std::vector<cavc_vertex> round_vertexes = readVertexes(round_offset);
  std::vector<cavc_vertex> square_vertexes = readVertexes(square_offset);
  ASSERT_FALSE(round_vertexes.empty());
  ASSERT_FALSE(square_vertexes.empty());

  auto endpoint_delta = [](cavc_vertex const &a, cavc_vertex const &b) {
    return std::hypot(a.x - b.x, a.y - b.y);
  };

  EXPECT_NEAR(endpoint_delta(square_vertexes.front(), round_vertexes.front()), 0.25, 1e-9);
  EXPECT_NEAR(endpoint_delta(square_vertexes.back(), round_vertexes.back()), 0.25, 1e-9);
  EXPECT_NEAR(cavc_get_path_length(square_offset) - cavc_get_path_length(round_offset), 0.5, 1e-8);
}

TEST(CApiRegression, ParallelOffsetSquareEndCapSupportsArcEndSegmentPositiveDelta) {
  std::vector<cavc_vertex> arc_end_open = {{0.0, 0.0, 0.0}, {2.0, 0.0, 1.0}, {4.0, 0.0, 0.0}};
  PlinePtr pline(plineFromVertexes(arc_end_open, false));

  cavc_parallel_offset_options round_options = cavc_parallel_offset_default_options();
  cavc_parallel_offset_options square_options = cavc_parallel_offset_default_options();
  square_options.end_cap_type = CAVC_OFFSET_END_CAP_SQUARE;

  cavc_pline_list *raw_round = nullptr;
  cavc_parallel_offset(pline.get(), 0.25, &raw_round, round_options);
  PlineListPtr round_results(raw_round);

  cavc_pline_list *raw_square = nullptr;
  cavc_parallel_offset(pline.get(), 0.25, &raw_square, square_options);
  PlineListPtr square_results(raw_square);

  ASSERT_EQ(cavc_pline_list_count(round_results.get()), 1u);
  ASSERT_EQ(cavc_pline_list_count(square_results.get()), 1u);

  cavc_pline *round_offset = cavc_pline_list_get(round_results.get(), 0);
  cavc_pline *square_offset = cavc_pline_list_get(square_results.get(), 0);
  std::vector<cavc_vertex> round_vertexes = readVertexes(round_offset);
  std::vector<cavc_vertex> square_vertexes = readVertexes(square_offset);
  ASSERT_FALSE(round_vertexes.empty());
  ASSERT_FALSE(square_vertexes.empty());

  auto endpoint_delta = [](cavc_vertex const &a, cavc_vertex const &b) {
    return std::hypot(a.x - b.x, a.y - b.y);
  };

  EXPECT_NEAR(endpoint_delta(square_vertexes.front(), round_vertexes.front()), 0.25, 1e-9);
  EXPECT_NEAR(endpoint_delta(square_vertexes.back(), round_vertexes.back()), 0.25, 1e-9);
  EXPECT_NEAR(cavc_get_path_length(square_offset) - cavc_get_path_length(round_offset), 0.5, 1e-8);
}

TEST(CApiRegression, ParallelOffsetButtEndCapMatchesRoundOnSimpleLine) {
  std::vector<cavc_vertex> line = {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}};
  PlinePtr pline(plineFromVertexes(line, false));

  cavc_parallel_offset_options round_options = cavc_parallel_offset_default_options();
  cavc_parallel_offset_options butt_options = cavc_parallel_offset_default_options();
  butt_options.end_cap_type = CAVC_OFFSET_END_CAP_BUTT;

  cavc_pline_list *raw_round = nullptr;
  cavc_parallel_offset(pline.get(), -1.0, &raw_round, round_options);
  PlineListPtr round_results(raw_round);

  cavc_pline_list *raw_butt = nullptr;
  cavc_parallel_offset(pline.get(), -1.0, &raw_butt, butt_options);
  PlineListPtr butt_results(raw_butt);

  ASSERT_EQ(cavc_pline_list_count(round_results.get()), 1u);
  ASSERT_EQ(cavc_pline_list_count(butt_results.get()), 1u);

  std::vector<cavc_vertex> round_vertexes =
      readVertexes(cavc_pline_list_get(round_results.get(), 0));
  std::vector<cavc_vertex> butt_vertexes = readVertexes(cavc_pline_list_get(butt_results.get(), 0));
  ASSERT_EQ(round_vertexes.size(), butt_vertexes.size());
  EXPECT_THAT(round_vertexes, t::Pointwise(VertexEqual(), butt_vertexes));
}

TEST(CApiRegression, ParallelOffsetButtEndCapCanDifferFromRoundOnMixedLineArcPolyline) {
  std::vector<cavc_vertex> mixed_open = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.5}, {2.0, 1.0, 0.0}};
  PlinePtr pline(plineFromVertexes(mixed_open, false));

  cavc_parallel_offset_options round_options = cavc_parallel_offset_default_options();
  cavc_parallel_offset_options butt_options = cavc_parallel_offset_default_options();
  butt_options.end_cap_type = CAVC_OFFSET_END_CAP_BUTT;

  auto run_offset = [&](cavc_real delta) {
    cavc_pline_list *raw_round = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_round, round_options);
    PlineListPtr round_results(raw_round);

    cavc_pline_list *raw_butt = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_butt, butt_options);
    PlineListPtr butt_results(raw_butt);

    ASSERT_EQ(cavc_pline_list_count(round_results.get()), 1u);
    ASSERT_EQ(cavc_pline_list_count(butt_results.get()), 1u);

    std::vector<cavc_vertex> round_vertexes =
        readVertexes(cavc_pline_list_get(round_results.get(), 0));
    std::vector<cavc_vertex> butt_vertexes =
        readVertexes(cavc_pline_list_get(butt_results.get(), 0));
    ASSERT_FALSE(round_vertexes.empty());
    ASSERT_FALSE(butt_vertexes.empty());

    std::vector<std::vector<cavc_vertex>> round_wrap = {round_vertexes};
    std::vector<std::vector<cavc_vertex>> butt_wrap = {butt_vertexes};
    EXPECT_THAT(round_wrap, t::Not(t::Pointwise(VertexListsFuzzyEqual(false), butt_wrap)));
    EXPECT_EQ(butt_vertexes.size(), round_vertexes.size() + 1);
  };

  run_offset(0.2);
  run_offset(-0.2);
}

TEST(CApiRegression, ParallelOffsetButtEndCapCanDifferFromRoundWhenArcStartsAtEndpoint) {
  std::vector<cavc_vertex> arc_start_open = {
      {0.0, 0.0, 0.5}, {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {3.0, 1.0, 0.0}};
  PlinePtr pline(plineFromVertexes(arc_start_open, false));

  cavc_parallel_offset_options round_options = cavc_parallel_offset_default_options();
  cavc_parallel_offset_options butt_options = cavc_parallel_offset_default_options();
  butt_options.end_cap_type = CAVC_OFFSET_END_CAP_BUTT;

  cavc_pline_list *raw_round = nullptr;
  cavc_parallel_offset(pline.get(), 0.5, &raw_round, round_options);
  PlineListPtr round_results(raw_round);

  cavc_pline_list *raw_butt = nullptr;
  cavc_parallel_offset(pline.get(), 0.5, &raw_butt, butt_options);
  PlineListPtr butt_results(raw_butt);

  ASSERT_EQ(cavc_pline_list_count(round_results.get()), 1u);
  ASSERT_EQ(cavc_pline_list_count(butt_results.get()), 1u);

  std::vector<cavc_vertex> round_vertexes =
      readVertexes(cavc_pline_list_get(round_results.get(), 0));
  std::vector<cavc_vertex> butt_vertexes = readVertexes(cavc_pline_list_get(butt_results.get(), 0));
  ASSERT_FALSE(round_vertexes.empty());
  ASSERT_FALSE(butt_vertexes.empty());

  std::vector<std::vector<cavc_vertex>> round_wrap = {round_vertexes};
  std::vector<std::vector<cavc_vertex>> butt_wrap = {butt_vertexes};
  EXPECT_THAT(round_wrap, t::Not(t::Pointwise(VertexListsFuzzyEqual(false), butt_wrap)));
  EXPECT_EQ(butt_vertexes.size(), round_vertexes.size() + 1);
}

TEST(CApiRegression, ParallelOffsetButtEndCapStressMatrix) {
  struct ButtMatrixCase {
    const char *name;
    std::vector<cavc_vertex> vertexes;
    uint32_t min_extra_positive;
    uint32_t min_extra_negative;
  };

  std::vector<ButtMatrixCase> cases;
  cases.push_back(ButtMatrixCase{
      "high_curvature",
      {
          {0.000, 0.000, -1.000},
          {1.726, -1.174, 0.000},
          {3.078, -0.618, 0.000},
          {3.734, -1.393, 0.050},
          {4.851, -0.804, 0.900},
          {6.329, 0.244, 0.500},
          {6.779, 1.285, 0.000},
      },
      1,
      1,
  });
  cases.push_back(ButtMatrixCase{
      "near_tangent",
      {
          {0.000, 0.000, -0.700},
          {0.443, 0.378, 0.080},
          {1.813, 0.591, -0.500},
          {2.839, 0.227, 0.000},
      },
      1,
      1,
  });
  cases.push_back(ButtMatrixCase{
      "multi_intersection_open",
      {
          {0.000, 0.000, 0.080},
          {0.714, -0.628, -0.700},
          {1.378, 0.844, -0.700},
          {3.126, -0.527, 0.000},
          {4.407, -0.858, 0.300},
          {5.060, -1.594, 0.000},
          {5.883, -0.695, -0.300},
          {6.388, -1.812, 0.000},
      },
      2,
      2,
  });

  auto run_case = [](ButtMatrixCase const &test_case, cavc_real delta, uint32_t min_extra_count) {
    PlinePtr pline(plineFromVertexes(test_case.vertexes, false));

    cavc_parallel_offset_options round_options = cavc_parallel_offset_default_options();
    cavc_parallel_offset_options butt_options = cavc_parallel_offset_default_options();
    butt_options.end_cap_type = CAVC_OFFSET_END_CAP_BUTT;

    cavc_pline_list *raw_round = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_round, round_options);
    PlineListPtr round_results(raw_round);

    cavc_pline_list *raw_butt = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_butt, butt_options);
    PlineListPtr butt_results(raw_butt);

    ASSERT_EQ(cavc_pline_list_count(round_results.get()), 1u);
    ASSERT_EQ(cavc_pline_list_count(butt_results.get()), 1u);

    cavc_pline *round_offset = cavc_pline_list_get(round_results.get(), 0);
    cavc_pline *butt_offset = cavc_pline_list_get(butt_results.get(), 0);
    std::vector<cavc_vertex> round_vertexes = readVertexes(round_offset);
    std::vector<cavc_vertex> butt_vertexes = readVertexes(butt_offset);
    ASSERT_FALSE(round_vertexes.empty());
    ASSERT_FALSE(butt_vertexes.empty());

    std::vector<std::vector<cavc_vertex>> round_wrap = {round_vertexes};
    std::vector<std::vector<cavc_vertex>> butt_wrap = {butt_vertexes};
    EXPECT_THAT(round_wrap, t::Not(t::Pointwise(VertexListsFuzzyEqual(false), butt_wrap)));
    EXPECT_GE(butt_vertexes.size(), round_vertexes.size() + min_extra_count);
    EXPECT_TRUE(vertexListAllFinite(round_vertexes));
    EXPECT_TRUE(vertexListAllFinite(butt_vertexes));
    EXPECT_GT(cavc_get_path_length(round_offset), 0.0);
    EXPECT_GT(cavc_get_path_length(butt_offset), 0.0);
  };

  for (auto const &test_case : cases) {
    SCOPED_TRACE(test_case.name);
    run_case(test_case, 0.2, test_case.min_extra_positive);
    run_case(test_case, -0.2, test_case.min_extra_negative);
  }
}

TEST(CApiRegression, ParallelOffsetEndCapDegenerateInputStressMatrix) {
  struct EndCapMatrixCase {
    const char *name;
    std::vector<cavc_vertex> vertexes;
    uint32_t min_square_extra_positive;
    uint32_t min_square_extra_negative;
  };

  std::vector<EndCapMatrixCase> cases;
  cases.push_back(EndCapMatrixCase{
      "repeated_endpoints",
      {
          {1.274, 0.067, -0.237},
          {1.274, 0.067, -0.488},
          {1.913, 0.833, -0.735},
          {3.039, 0.660, 0.160},
          {4.423, 0.165, -0.129},
          {4.423, 0.165, -0.129},
          {5.062, 0.863, 0.000},
          {6.105, 0.891, 0.000},
          {6.105, 0.891, -0.327},
      },
      1,
      1,
  });
  cases.push_back(EndCapMatrixCase{
      "near_zero_radius_arc",
      {
          {0.000, 0.000, 0.966},
          {0.133, 0.000, -0.911},
          {0.272, 0.000, 0.000},
          {1.272, 0.189, 0.000},
          {2.272, 0.802, -0.201},
          {3.272, -0.049, 0.000},
      },
      1,
      1,
  });
  cases.push_back(EndCapMatrixCase{
      "high_density_self_intersection_open",
      {
          {0.000, 0.000, 0.513},
          {1.647, 1.281, 0.000},
          {2.447, -1.499, 0.000},
          {3.167, 1.313, 0.118},
          {3.739, -0.908, 0.000},
          {4.688, 1.162, -0.771},
          {6.126, 0.000, 0.000},
      },
      1,
      1,
  });

  auto run_case = [](EndCapMatrixCase const &test_case, cavc_real delta,
                     uint32_t min_square_extra) {
    PlinePtr pline(plineFromVertexes(test_case.vertexes, false));

    cavc_parallel_offset_options round_options = cavc_parallel_offset_default_options();
    cavc_parallel_offset_options square_options = cavc_parallel_offset_default_options();
    square_options.end_cap_type = CAVC_OFFSET_END_CAP_SQUARE;
    cavc_parallel_offset_options butt_options = cavc_parallel_offset_default_options();
    butt_options.end_cap_type = CAVC_OFFSET_END_CAP_BUTT;

    cavc_pline_list *raw_round = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_round, round_options);
    PlineListPtr round_results(raw_round);

    cavc_pline_list *raw_square = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_square, square_options);
    PlineListPtr square_results(raw_square);

    cavc_pline_list *raw_butt = nullptr;
    cavc_parallel_offset(pline.get(), delta, &raw_butt, butt_options);
    PlineListPtr butt_results(raw_butt);

    ASSERT_EQ(cavc_pline_list_count(round_results.get()), 1u);
    ASSERT_EQ(cavc_pline_list_count(square_results.get()), 1u);
    ASSERT_EQ(cavc_pline_list_count(butt_results.get()), 1u);

    cavc_pline *round_pline = cavc_pline_list_get(round_results.get(), 0);
    cavc_pline *square_pline = cavc_pline_list_get(square_results.get(), 0);
    cavc_pline *butt_pline = cavc_pline_list_get(butt_results.get(), 0);

    std::vector<cavc_vertex> round_vertexes = readVertexes(round_pline);
    std::vector<cavc_vertex> square_vertexes = readVertexes(square_pline);
    std::vector<cavc_vertex> butt_vertexes = readVertexes(butt_pline);

    ASSERT_FALSE(round_vertexes.empty());
    ASSERT_FALSE(square_vertexes.empty());
    ASSERT_FALSE(butt_vertexes.empty());

    std::vector<std::vector<cavc_vertex>> round_wrap = {round_vertexes};
    std::vector<std::vector<cavc_vertex>> square_wrap = {square_vertexes};
    std::vector<std::vector<cavc_vertex>> butt_wrap = {butt_vertexes};
    EXPECT_THAT(square_wrap, t::Not(t::Pointwise(VertexListsFuzzyEqual(false), round_wrap)));
    EXPECT_THAT(butt_wrap, t::Not(t::Pointwise(VertexListsFuzzyEqual(false), round_wrap)));
    EXPECT_GE(square_vertexes.size(), round_vertexes.size() + min_square_extra);

    EXPECT_TRUE(vertexListAllFinite(round_vertexes));
    EXPECT_TRUE(vertexListAllFinite(square_vertexes));
    EXPECT_TRUE(vertexListAllFinite(butt_vertexes));
    EXPECT_GT(cavc_get_path_length(round_pline), 0.0);
    EXPECT_GT(cavc_get_path_length(square_pline), 0.0);
    EXPECT_GT(cavc_get_path_length(butt_pline), 0.0);
  };

  for (auto const &test_case : cases) {
    SCOPED_TRACE(test_case.name);
    run_case(test_case, 0.2, test_case.min_square_extra_positive);
    run_case(test_case, -0.2, test_case.min_square_extra_negative);
  }
}

TEST(CApiRegression, SpatialIndexClosedPolylineQueryAndItemCount) {
  std::vector<cavc_vertex> square = {
      {0.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {4.0, 2.0, 0.0}, {0.0, 2.0, 0.0}};
  PlinePtr pline(plineFromVertexes(square, true));
  SpatialIndexPtr index(cavc_spatial_index_create(pline.get()), cavc_spatial_index_delete);

  ASSERT_NE(index.get(), nullptr);
  EXPECT_EQ(cavc_spatial_index_item_count(index.get()), 4u);

  std::vector<uint32_t> right_hits = querySpatialIndexSorted(index.get(), 3.8, 0.5, 4.2, 1.5);
  std::vector<uint32_t> expected_right = {1u};
  EXPECT_THAT(right_hits, t::Pointwise(t::Eq(), expected_right));

  std::vector<uint32_t> all_hits = querySpatialIndexSorted(index.get(), -1.0, -1.0, 5.0, 3.0);
  std::vector<uint32_t> expected_all = {0u, 1u, 2u, 3u};
  EXPECT_THAT(all_hits, t::Pointwise(t::Eq(), expected_all));
}

TEST(CApiRegression, SpatialIndexOpenPolylineQueryAndItemCount) {
  std::vector<cavc_vertex> open_rect = {
      {0.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {4.0, 2.0, 0.0}, {0.0, 2.0, 0.0}};
  PlinePtr pline(plineFromVertexes(open_rect, false));
  SpatialIndexPtr index(cavc_spatial_index_create(pline.get()), cavc_spatial_index_delete);

  ASSERT_NE(index.get(), nullptr);
  EXPECT_EQ(cavc_spatial_index_item_count(index.get()), 3u);

  std::vector<uint32_t> right_hits = querySpatialIndexSorted(index.get(), 3.8, 0.5, 4.2, 1.5);
  std::vector<uint32_t> expected_right = {1u};
  EXPECT_THAT(right_hits, t::Pointwise(t::Eq(), expected_right));

  // Open polyline has no closing segment from last to first, so this query should miss.
  std::vector<uint32_t> left_hits = querySpatialIndexSorted(index.get(), -0.2, 0.5, 0.2, 1.5);
  EXPECT_TRUE(left_hits.empty());
}

TEST(CApiRegression, EmptyPolylineExtentsReturnInfiniteBounds) {
  PlinePtr pline(cavc_pline_new(nullptr, 0, 0));

  cavc_real min_x = 0.0;
  cavc_real min_y = 0.0;
  cavc_real max_x = 0.0;
  cavc_real max_y = 0.0;
  cavc_get_extents(pline.get(), &min_x, &min_y, &max_x, &max_y);

  EXPECT_TRUE(std::isinf(min_x));
  EXPECT_TRUE(std::isinf(min_y));
  EXPECT_TRUE(std::isinf(max_x));
  EXPECT_TRUE(std::isinf(max_y));
  EXPECT_GT(min_x, 0.0);
  EXPECT_GT(min_y, 0.0);
  EXPECT_LT(max_x, 0.0);
  EXPECT_LT(max_y, 0.0);
  EXPECT_DOUBLE_EQ(cavc_get_area(pline.get()), 0.0);
  EXPECT_DOUBLE_EQ(cavc_get_path_length(pline.get()), 0.0);
}

TEST(CApiRegression, CombinePlinesDoesNotMutateInputs) {
  std::vector<cavc_vertex> pline_a_vertexes = {
      {0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}, {10.0, 5.0, 0.0}, {0.0, 5.0, 0.0}};
  std::vector<cavc_vertex> pline_b_vertexes = {
      {5.0, 1.0, 0.0}, {12.0, 1.0, 0.0}, {12.0, 7.0, 0.0}, {5.0, 7.0, 0.0}};

  PlinePtr pline_a(plineFromVertexes(pline_a_vertexes, true));
  PlinePtr pline_b(plineFromVertexes(pline_b_vertexes, true));

  std::vector<cavc_vertex> before_a = readVertexes(pline_a.get());
  std::vector<cavc_vertex> before_b = readVertexes(pline_b.get());

  cavc_pline_list *raw_remaining = nullptr;
  cavc_pline_list *raw_subtracted = nullptr;
  cavc_combine_plines(pline_a.get(), pline_b.get(), 0, &raw_remaining, &raw_subtracted);
  PlineListPtr remaining(raw_remaining);
  PlineListPtr subtracted(raw_subtracted);

  ASSERT_NE(remaining.get(), nullptr);
  ASSERT_NE(subtracted.get(), nullptr);
  EXPECT_GE(cavc_pline_list_count(remaining.get()), 1u);

  std::vector<cavc_vertex> after_a = readVertexes(pline_a.get());
  std::vector<cavc_vertex> after_b = readVertexes(pline_b.get());
  EXPECT_THAT(after_a, t::Pointwise(VertexEqual(), before_a));
  EXPECT_THAT(after_b, t::Pointwise(VertexEqual(), before_b));
}

TEST(CApiRegression, InvertDirectionPreservesPathAndNegatesArea) {
  std::vector<cavc_vertex> original_vertexes = {
      {0.0, 0.0, 0.0}, {4.0, 0.0, 0.414213562373095}, {4.0, 4.0, 0.0}, {0.0, 4.0, 0.0}};
  PlinePtr pline(plineFromVertexes(original_vertexes, true));

  cavc_real original_area = cavc_get_area(pline.get());
  cavc_real original_path = cavc_get_path_length(pline.get());

  std::vector<cavc_vertex> expected_reversed = original_vertexes;
  reverseDirection(expected_reversed);

  cavc_pline_invert_direction(pline.get());
  std::vector<cavc_vertex> actual_reversed = readVertexes(pline.get());
  EXPECT_THAT(actual_reversed, t::Pointwise(VertexFuzzyEqual(), expected_reversed));
  EXPECT_NEAR(cavc_get_area(pline.get()), -original_area, TEST_EPSILON());
  EXPECT_NEAR(cavc_get_path_length(pline.get()), original_path, TEST_EPSILON());

  cavc_pline_invert_direction(pline.get());
  std::vector<cavc_vertex> roundtrip = readVertexes(pline.get());
  EXPECT_THAT(roundtrip, t::Pointwise(VertexFuzzyEqual(), original_vertexes));
}

TEST(CApiRegression, PruneSingularitiesRemovesDuplicatePositions) {
  std::vector<cavc_vertex> with_singularities = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.2}, {1.0, 0.0, 0.3},
                                                 {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}};
  PlinePtr input(plineFromVertexes(with_singularities, true));

  cavc_pline *raw_output = nullptr;
  cavc_pline_prune_singularities(input.get(), 1e-5, &raw_output);
  PlinePtr output(raw_output);

  ASSERT_NE(output.get(), nullptr);
  EXPECT_EQ(cavc_pline_vertex_count(input.get()), with_singularities.size());
  EXPECT_EQ(cavc_pline_vertex_count(output.get()), 4u);
  EXPECT_EQ(cavc_pline_is_closed(output.get()), 1);

  std::vector<cavc_vertex> expected = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.3}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};
  std::vector<cavc_vertex> actual = readVertexes(output.get());
  EXPECT_THAT(actual, t::Pointwise(VertexEqual(), expected));
}

TEST(CApiRegression, RemoveRedundantSimplifiesInPlace) {
  std::vector<cavc_vertex> redundant = {
      {0.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 2.0, 0.0}, {2.0, 2.0, 0.2}, {2.0, 3.0, 0.0}};
  PlinePtr pline(plineFromVertexes(redundant, false));

  cavc_pline_remove_redundant(pline.get(), 1e-5);

  std::vector<cavc_vertex> expected = {
      {0.0, 0.0, 0.0}, {2.0, 2.0, 0.2}, {2.0, 3.0, 0.0}};
  std::vector<cavc_vertex> actual = readVertexes(pline.get());
  EXPECT_THAT(actual, t::Pointwise(VertexFuzzyEqual(), expected));
  EXPECT_EQ(cavc_pline_is_closed(pline.get()), 0);
}

TEST(CApiRegression, ConvertArcsToLinesProducesLineSegments) {
  const cavc_real radius = 2.0;
  std::vector<cavc_vertex> circle = {{0.0, 0.0, 1.0}, {2.0 * radius, 0.0, 1.0}};
  PlinePtr input(plineFromVertexes(circle, true));

  cavc_pline *raw_output = nullptr;
  cavc_pline_convert_arcs_to_lines(input.get(), 0.01, &raw_output);
  PlinePtr output(raw_output);

  ASSERT_NE(output.get(), nullptr);
  EXPECT_EQ(cavc_pline_is_closed(output.get()), 1);
  EXPECT_GT(cavc_pline_vertex_count(output.get()), cavc_pline_vertex_count(input.get()));

  std::vector<cavc_vertex> output_vertexes = readVertexes(output.get());
  for (const cavc_vertex &v : output_vertexes) {
    EXPECT_NEAR(v.bulge, 0.0, TEST_EPSILON());
  }

  const cavc_real expected_area = PI() * radius * radius;
  EXPECT_NEAR(std::abs(cavc_get_area(output.get())), expected_area, 0.05 * expected_area);

  std::vector<cavc_vertex> input_after = readVertexes(input.get());
  EXPECT_THAT(input_after, t::Pointwise(VertexEqual(), circle));
}

TEST(CApiRegression, NormalizePrunesSingularitiesAndKeepsWinding) {
  std::vector<cavc_vertex> with_singularities = {{0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {2.0, 2.0, 0.2},
                                                 {2.0, 2.0, 0.3}, {2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  PlinePtr input(plineFromVertexes(with_singularities, true));

  const cavc_real area_before = cavc_get_area(input.get());

  cavc_pline *raw_output = nullptr;
  cavc_pline_normalize(input.get(), 1e-5, CAVC_WINDING_KEEP, &raw_output);
  PlinePtr output(raw_output);

  ASSERT_NE(output.get(), nullptr);
  EXPECT_EQ(cavc_pline_vertex_count(input.get()), with_singularities.size());
  EXPECT_EQ(cavc_pline_vertex_count(output.get()), 4u);
  EXPECT_LT(area_before, 0.0);
  EXPECT_LT(cavc_get_area(output.get()), 0.0);
}

TEST(CApiRegression, NormalizeCanForceCounterClockwiseWinding) {
  std::vector<cavc_vertex> clockwise = {
      {0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {2.0, 2.0, 0.0}, {2.0, 0.0, 0.0}};
  PlinePtr input(plineFromVertexes(clockwise, true));

  ASSERT_LT(cavc_get_area(input.get()), 0.0);
  cavc_pline *raw_output = nullptr;
  cavc_pline_normalize(input.get(), 1e-5, CAVC_WINDING_COUNTER_CLOCKWISE, &raw_output);
  PlinePtr output(raw_output);

  ASSERT_NE(output.get(), nullptr);
  EXPECT_GT(cavc_get_area(output.get()), 0.0);
  EXPECT_EQ(cavc_pline_vertex_count(output.get()), 4u);
}

TEST(CApiRegression, NormalizeCanForceClockwiseWinding) {
  std::vector<cavc_vertex> counter_clockwise = {
      {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 2.0, 0.0}, {0.0, 2.0, 0.0}};
  PlinePtr input(plineFromVertexes(counter_clockwise, true));

  ASSERT_GT(cavc_get_area(input.get()), 0.0);
  cavc_pline *raw_output = nullptr;
  cavc_pline_normalize(input.get(), 1e-5, CAVC_WINDING_CLOCKWISE, &raw_output);
  PlinePtr output(raw_output);

  ASSERT_NE(output.get(), nullptr);
  EXPECT_LT(cavc_get_area(output.get()), 0.0);
  EXPECT_EQ(cavc_pline_vertex_count(output.get()), 4u);
}

TEST(CApiRegression, PointContainmentClassifiesClosedPolyline) {
  std::vector<cavc_vertex> square = {
      {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 2.0, 0.0}, {0.0, 2.0, 0.0}};
  PlinePtr pline(plineFromVertexes(square, true));

  EXPECT_EQ(cavc_get_point_containment(pline.get(), cavc_point{1.0, 1.0}, 1e-6), CAVC_POINT_INSIDE);
  EXPECT_EQ(cavc_get_point_containment(pline.get(), cavc_point{3.0, 1.0}, 1e-6),
            CAVC_POINT_OUTSIDE);
  EXPECT_EQ(cavc_get_point_containment(pline.get(), cavc_point{0.0, 1.0}, 1e-6),
            CAVC_POINT_ON_BOUNDARY);
  EXPECT_EQ(cavc_get_point_containment(pline.get(), cavc_point{0.0, 0.0}, 1e-6),
            CAVC_POINT_ON_BOUNDARY);
}

TEST(CApiRegression, PointContainmentForOpenPolylineHasNoInside) {
  std::vector<cavc_vertex> seg = {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}};
  PlinePtr pline(plineFromVertexes(seg, false));

  EXPECT_EQ(cavc_get_point_containment(pline.get(), cavc_point{1.0, 0.0}, 1e-6),
            CAVC_POINT_ON_BOUNDARY);
  EXPECT_EQ(cavc_get_point_containment(pline.get(), cavc_point{1.0, 0.1}, 1e-6),
            CAVC_POINT_OUTSIDE);
  EXPECT_EQ(cavc_get_point_containment(pline.get(), cavc_point{1.0, 0.1}, 0.2),
            CAVC_POINT_ON_BOUNDARY);
}

TEST(CApiRegression, OffsetLoopTopologyBuildProducesStableOrderAndParentChild) {
  auto island_vertexes = makeAxisAlignedRectLoopVertexes(12.0, 4.0, 16.0, 8.0, false);
  auto outer_vertexes = makeAxisAlignedRectLoopVertexes(0.0, 0.0, 20.0, 20.0, false);
  auto hole_b_vertexes = makeAxisAlignedRectLoopVertexes(10.0, 2.0, 18.0, 10.0, true);
  auto hole_a_vertexes = makeAxisAlignedRectLoopVertexes(2.0, 2.0, 8.0, 8.0, true);
  auto hole_c_vertexes = makeAxisAlignedRectLoopVertexes(13.0, 5.0, 15.0, 7.0, true);

  PlinePtr island(plineFromVertexes(island_vertexes, true));
  PlinePtr outer(plineFromVertexes(outer_vertexes, true));
  PlinePtr hole_b(plineFromVertexes(hole_b_vertexes, true));
  PlinePtr hole_a(plineFromVertexes(hole_a_vertexes, true));
  PlinePtr hole_c(plineFromVertexes(hole_c_vertexes, true));

  cavc_pline const *ccw_loops[] = {island.get(), outer.get()};
  cavc_pline const *cw_loops[] = {hole_b.get(), hole_a.get(), hole_c.get()};

  TopologyPtr topology(
      cavc_offset_loop_topology_build(ccw_loops, 2, cw_loops, 3, static_cast<cavc_real>(1e-5)));
  ASSERT_NE(topology.get(), nullptr);
  ASSERT_EQ(cavc_offset_loop_topology_count(topology.get()), 5u);

  cavc_offset_loop_topology_node node0 = cavc_offset_loop_topology_get(topology.get(), 0);
  cavc_offset_loop_topology_node node1 = cavc_offset_loop_topology_get(topology.get(), 1);
  cavc_offset_loop_topology_node node2 = cavc_offset_loop_topology_get(topology.get(), 2);
  cavc_offset_loop_topology_node node3 = cavc_offset_loop_topology_get(topology.get(), 3);
  cavc_offset_loop_topology_node node4 = cavc_offset_loop_topology_get(topology.get(), 4);

  EXPECT_EQ(node0.role, CAVC_OFFSET_LOOP_ROLE_OUTER);
  EXPECT_EQ(node0.source_index, 1u);
  EXPECT_EQ(node1.role, CAVC_OFFSET_LOOP_ROLE_HOLE);
  EXPECT_EQ(node1.source_index, 0u);
  EXPECT_EQ(node2.role, CAVC_OFFSET_LOOP_ROLE_HOLE);
  EXPECT_EQ(node2.source_index, 1u);
  EXPECT_EQ(node3.role, CAVC_OFFSET_LOOP_ROLE_OUTER);
  EXPECT_EQ(node3.source_index, 0u);
  EXPECT_EQ(node4.role, CAVC_OFFSET_LOOP_ROLE_HOLE);
  EXPECT_EQ(node4.source_index, 2u);

  uint32_t outer_index = findTopologyNodeIndex(topology.get(), CAVC_OFFSET_LOOP_ROLE_OUTER, 1u);
  uint32_t hole_b_index = findTopologyNodeIndex(topology.get(), CAVC_OFFSET_LOOP_ROLE_HOLE, 0u);
  uint32_t hole_a_index = findTopologyNodeIndex(topology.get(), CAVC_OFFSET_LOOP_ROLE_HOLE, 1u);
  uint32_t island_index = findTopologyNodeIndex(topology.get(), CAVC_OFFSET_LOOP_ROLE_OUTER, 0u);
  uint32_t hole_c_index = findTopologyNodeIndex(topology.get(), CAVC_OFFSET_LOOP_ROLE_HOLE, 2u);

  ASSERT_NE(outer_index, CAVC_OFFSET_LOOP_NO_PARENT);
  ASSERT_NE(hole_b_index, CAVC_OFFSET_LOOP_NO_PARENT);
  ASSERT_NE(hole_a_index, CAVC_OFFSET_LOOP_NO_PARENT);
  ASSERT_NE(island_index, CAVC_OFFSET_LOOP_NO_PARENT);
  ASSERT_NE(hole_c_index, CAVC_OFFSET_LOOP_NO_PARENT);

  EXPECT_EQ(cavc_offset_loop_topology_get(topology.get(), outer_index).parent_index,
            CAVC_OFFSET_LOOP_NO_PARENT);
  EXPECT_EQ(cavc_offset_loop_topology_get(topology.get(), hole_b_index).parent_index, outer_index);
  EXPECT_EQ(cavc_offset_loop_topology_get(topology.get(), hole_a_index).parent_index, outer_index);
  EXPECT_EQ(cavc_offset_loop_topology_get(topology.get(), island_index).parent_index, hole_b_index);
  EXPECT_EQ(cavc_offset_loop_topology_get(topology.get(), hole_c_index).parent_index, island_index);
}

TEST(CApiRegression, OffsetLoopTopologySupportsOuterOnlyInput) {
  auto small_vertexes = makeAxisAlignedRectLoopVertexes(0.0, 0.0, 2.0, 2.0, false);
  auto large_vertexes = makeAxisAlignedRectLoopVertexes(10.0, 0.0, 20.0, 8.0, false);

  PlinePtr small(plineFromVertexes(small_vertexes, true));
  PlinePtr large(plineFromVertexes(large_vertexes, true));

  cavc_pline const *ccw_loops[] = {small.get(), large.get()};

  TopologyPtr topology(
      cavc_offset_loop_topology_build(ccw_loops, 2, nullptr, 0, static_cast<cavc_real>(1e-5)));
  ASSERT_NE(topology.get(), nullptr);
  ASSERT_EQ(cavc_offset_loop_topology_count(topology.get()), 2u);

  cavc_offset_loop_topology_node node0 = cavc_offset_loop_topology_get(topology.get(), 0);
  cavc_offset_loop_topology_node node1 = cavc_offset_loop_topology_get(topology.get(), 1);

  EXPECT_EQ(node0.role, CAVC_OFFSET_LOOP_ROLE_OUTER);
  EXPECT_EQ(node0.source_index, 1u);
  EXPECT_EQ(node0.parent_index, CAVC_OFFSET_LOOP_NO_PARENT);
  EXPECT_EQ(node1.role, CAVC_OFFSET_LOOP_ROLE_OUTER);
  EXPECT_EQ(node1.source_index, 0u);
  EXPECT_EQ(node1.parent_index, CAVC_OFFSET_LOOP_NO_PARENT);
}

TEST(CApiRegression, TolerancesSetGetResetRoundTrip) {
  cavc_tolerances defaults{};
  cavc_get_tolerances(&defaults);
  EXPECT_GT(defaults.real_threshold, 0.0);
  EXPECT_GT(defaults.real_precision, 0.0);
  EXPECT_GT(defaults.slice_join_threshold, 0.0);
  EXPECT_GT(defaults.offset_threshold, 0.0);

  cavc_tolerances custom{2e-7, 3e-6, 4e-5, 5e-5};
  cavc_set_tolerances(&custom);

  cavc_tolerances after_set{};
  cavc_get_tolerances(&after_set);
  EXPECT_EQ(after_set.real_threshold, custom.real_threshold);
  EXPECT_EQ(after_set.real_precision, custom.real_precision);
  EXPECT_EQ(after_set.slice_join_threshold, custom.slice_join_threshold);
  EXPECT_EQ(after_set.offset_threshold, custom.offset_threshold);

  cavc_reset_tolerances();
  cavc_tolerances after_reset{};
  cavc_get_tolerances(&after_reset);
  EXPECT_EQ(after_reset.real_threshold, defaults.real_threshold);
  EXPECT_EQ(after_reset.real_precision, defaults.real_precision);
  EXPECT_EQ(after_reset.slice_join_threshold, defaults.slice_join_threshold);
  EXPECT_EQ(after_reset.offset_threshold, defaults.offset_threshold);
}

TEST(CApiRegression, ParallelOffsetReportedShapeCaseReproducesSelfIntersection) {
  const std::vector<cavc_vertex> input_vertexes = makeShapeOffsetSelfIntersectInputVertexes();
  ASSERT_EQ(input_vertexes.size(), 249u);

  PlinePtr input(plineFromVertexes(input_vertexes, true));
  ASSERT_NE(input.get(), nullptr);
  EXPECT_LT(cavc_get_area(input.get()), 0.0);

  cavc_pline *raw_normalized = nullptr;
  cavc_pline_normalize(input.get(), 1e-5, CAVC_WINDING_COUNTER_CLOCKWISE, &raw_normalized);
  PlinePtr normalized(raw_normalized);
  ASSERT_NE(normalized.get(), nullptr);
  EXPECT_GT(cavc_get_area(normalized.get()), 0.0);

  cavc_pline_list *raw_results = nullptr;
  cavc_parallel_offset(normalized.get(), -0.125, &raw_results, defaultParallelOffsetOptions());
  PlineListPtr results(raw_results);

  ASSERT_NE(results.get(), nullptr);
  ASSERT_EQ(cavc_pline_list_count(results.get()), 1u);
  uint32_t open_count = 0;
  uint32_t self_intersecting_count = 0;
  for (uint32_t i = 0; i < cavc_pline_list_count(results.get()); ++i) {
    cavc_pline const *result = cavc_pline_list_get(results.get(), i);
    open_count += cavc_pline_is_closed(result) == 0 ? 1u : 0u;
    self_intersecting_count += hasSelfIntersect(result) ? 1u : 0u;
  }

  EXPECT_EQ(open_count, 0u);
  EXPECT_EQ(self_intersecting_count, 0u);
}

TEST(CApiRegression, ParallelOffsetCollapsedCombPrefersSimpleClosedLoop) {
  PlinePtr input(plineFromVertexes(makeCollapsedCombVertexes(), true));
  ASSERT_NE(input.get(), nullptr);

  cavc_pline_list *raw_results = nullptr;
  cavc_parallel_offset(input.get(), -1.0, &raw_results, defaultParallelOffsetOptions());
  PlineListPtr results(raw_results);

  ASSERT_NE(results.get(), nullptr);
  ASSERT_EQ(cavc_pline_list_count(results.get()), 1u);

  cavc_pline *result = cavc_pline_list_get(results.get(), 0);
  EXPECT_EQ(cavc_pline_is_closed(result), 1);
  EXPECT_FALSE(hasSelfIntersect(result));
  EXPECT_GT(cavc_get_area(result), 1000.0);
}

int main(int argc, char **argv) {
  t::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
