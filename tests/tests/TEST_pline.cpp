#include <gmock/gmock.h>

#include <gtest/gtest.h>

#include <cavc/mathutils.hpp>
#include <cavc/polyline.hpp>

#include "casebuilder.hpp"

using Polyline = cavc::Polyline<double>;
using PlineVertex = cavc::PlineVertex<double>;
using Vector2 = cavc::Vector2<double>;
using AABB = cavc::AABB<double>;

namespace {

// Helper function to check if two doubles are approximately equal
constexpr double EPSILON = 1e-9;
bool approxEqual(double a, double b, double epsilon = EPSILON) { return std::abs(a - b) < epsilon; }

bool approxEqual(const Vector2 &a, const Vector2 &b, double epsilon = EPSILON) {
  return approxEqual(a.x(), b.x(), epsilon) && approxEqual(a.y(), b.y(), epsilon);
}

bool approxEqual(const AABB &a, const AABB &b, double epsilon = EPSILON) {
  return approxEqual(a.xMin, b.xMin, epsilon) && approxEqual(a.yMin, b.yMin, epsilon) &&
         approxEqual(a.xMax, b.xMax, epsilon) && approxEqual(a.yMax, b.yMax, epsilon);
}

Polyline createFromVertices(const std::vector<PlineVertex> &vertices, bool isClosed = true) {
  Polyline pline;
  pline.isClosed() = isClosed;
  for (const auto &vertex : vertices) {
    pline.addVertex(vertex);
  }
  return pline;
}

} // namespace

// Basic Construction and Modification Tests

TEST(polyline, default_construction) {
  Polyline pline;
  EXPECT_EQ(pline.size(), 0);
  EXPECT_FALSE(pline.isClosed());
  EXPECT_TRUE(pline.vertexes().empty());
}

TEST(polyline, add_vertices) {
  Polyline pline;
  pline.addVertex(0.0, 0.0, 0.0);
  pline.addVertex(1.0, 0.0, 0.0);
  pline.addVertex(1.0, 1.0, 0.0);

  EXPECT_EQ(pline.size(), 3);
  EXPECT_TRUE(approxEqual(pline[0].pos(), Vector2{0.0, 0.0}));
  EXPECT_TRUE(approxEqual(pline[1].pos(), Vector2{1.0, 0.0}));
  EXPECT_TRUE(approxEqual(pline[2].pos(), Vector2{1.0, 1.0}));
}

TEST(polyline, add_vertex_from_plinevertex) {
  Polyline pline;
  PlineVertex vertex{2.0, 3.0, 0.5};
  pline.addVertex(vertex);

  EXPECT_EQ(pline.size(), 1);
  EXPECT_TRUE(approxEqual(pline[0].pos(), Vector2{2.0, 3.0}));
  EXPECT_TRUE(approxEqual(pline[0].bulge(), 0.5));
}

TEST(polyline, closed_open_state) {
  Polyline pline;
  EXPECT_FALSE(pline.isClosed());

  pline.isClosed() = true;
  EXPECT_TRUE(pline.isClosed());
}

TEST(polyline, last_vertex_access) {
  Polyline pline;
  pline.addVertex(1.0, 2.0, 0.0);
  pline.addVertex(3.0, 4.0, 0.5);

  EXPECT_TRUE(approxEqual(pline.lastVertex().pos(), Vector2{3.0, 4.0}));
  EXPECT_TRUE(approxEqual(pline.lastVertex().bulge(), 0.5));

  // Modify last vertex
  pline.lastVertex().bulge() = 1.0;
  EXPECT_TRUE(approxEqual(pline.lastVertex().bulge(), 1.0));
}

// Test with CaseBuilder patterns

TEST(polyline, simple_rectangle_case) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices);

  EXPECT_EQ(pline.size(), 4);
  EXPECT_TRUE(pline.isClosed());

  // Check vertices
  EXPECT_TRUE(approxEqual(pline[0].pos(), Vector2{0.0, 0.0}));
  EXPECT_TRUE(approxEqual(pline[1].pos(), Vector2{1.0, 0.0}));
  EXPECT_TRUE(approxEqual(pline[2].pos(), Vector2{1.0, 1.0}));
  EXPECT_TRUE(approxEqual(pline[3].pos(), Vector2{0.0, 1.0}));
}

TEST(polyline, positive_circle_case) {
  auto vertices = CaseBuilder::positiveCircle();
  Polyline pline = createFromVertices(vertices);

  EXPECT_EQ(pline.size(), 2);
  EXPECT_TRUE(approxEqual(pline[0].bulge(), 1.0));
  EXPECT_TRUE(approxEqual(pline[1].bulge(), 1.0));
}

TEST(polyline, negative_circle_case) {
  auto vertices = CaseBuilder::negativeCircle();
  Polyline pline = createFromVertices(vertices);

  EXPECT_EQ(pline.size(), 2);
  EXPECT_TRUE(approxEqual(pline[0].bulge(), -1.0));
  EXPECT_TRUE(approxEqual(pline[1].bulge(), -1.0));
}

TEST(polyline, figure_eight_case) {
  auto vertices = CaseBuilder::figureEightCase();
  Polyline pline = createFromVertices(vertices);

  EXPECT_EQ(pline.size(), 4);
  // Check that it has alternating positive and negative bulges
  EXPECT_TRUE(pline[0].bulge() > 0.0);
  EXPECT_TRUE(pline[1].bulge() > 0.0);
  EXPECT_TRUE(pline[2].bulge() < 0.0);
  EXPECT_TRUE(pline[3].bulge() < 0.0);
}

// Geometric Computation Tests

TEST(polyline, get_extents_empty_polyline) {
  Polyline pline;
  AABB extents = cavc::getExtents(pline);

  EXPECT_TRUE(std::isinf(extents.xMin) && extents.xMin > 0);
  EXPECT_TRUE(std::isinf(extents.yMin) && extents.yMin > 0);
  EXPECT_TRUE(std::isinf(extents.xMax) && extents.xMax < 0);
  EXPECT_TRUE(std::isinf(extents.yMax) && extents.yMax < 0);
}

TEST(polyline, get_extents_single_vertex) {
  Polyline pline;
  pline.addVertex(2.0, 3.0, 0.0);
  AABB extents = cavc::getExtents(pline);

  EXPECT_TRUE(approxEqual(extents.xMin, 2.0));
  EXPECT_TRUE(approxEqual(extents.yMin, 3.0));
  EXPECT_TRUE(approxEqual(extents.xMax, 2.0));
  EXPECT_TRUE(approxEqual(extents.yMax, 3.0));
}

TEST(polyline, get_extents_rectangle) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices);
  AABB extents = cavc::getExtents(pline);

  AABB expected{0.0, 0.0, 1.0, 1.0};
  EXPECT_TRUE(approxEqual(extents, expected));
}

TEST(polyline, get_extents_with_arcs) {
  auto vertices = CaseBuilder::positiveCircle();
  Polyline pline = createFromVertices(vertices);
  AABB extents = cavc::getExtents(pline);

  // Circle from (0,0) to (10,0) with bulge 1 creates a full circle
  // The circle should extend beyond the original points
  EXPECT_TRUE(extents.xMin <= 0.0);
  EXPECT_TRUE(extents.xMax >= 10.0);
  EXPECT_TRUE(extents.yMin < 0.0);
  EXPECT_TRUE(extents.yMax > 0.0);
}

TEST(polyline, get_area_open_polyline) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices, false); // Open polyline

  double area = cavc::getArea(pline);
  EXPECT_TRUE(approxEqual(area, 0.0)); // Open polylines have zero area
}

TEST(polyline, get_area_rectangle) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices);

  double area = cavc::getArea(pline);
  EXPECT_TRUE(approxEqual(area, 1.0)); // 1x1 rectangle

  auto reversed = CaseBuilder::reverseDirection(vertices);
  Polyline reversedPline = createFromVertices(reversed);
  double reversedArea = cavc::getArea(reversedPline);
  EXPECT_TRUE(approxEqual(reversedArea, -1.0)); // 1x1 rectangle
}

TEST(polyline, get_area_two_half_circle) {
  auto vertices = CaseBuilder::positiveCircle();
  Polyline pline = createFromVertices(vertices);

  double area = cavc::getArea(pline);
  EXPECT_TRUE(approxEqual(area, cavc::utils::pi<double>() * 25.0));

  auto reversed = CaseBuilder::reverseDirection(vertices);
  Polyline reversedPline = createFromVertices(reversed);
  double reversedArea = cavc::getArea(reversedPline);
  EXPECT_TRUE(approxEqual(reversedArea, -cavc::utils::pi<double>() * 25.0));
}

TEST(polyline, get_path_length_empty) {
  Polyline pline;
  double length = cavc::getPathLength(pline);
  EXPECT_TRUE(approxEqual(length, 0.0));
}

TEST(polyline, get_path_length_single_vertex) {
  Polyline pline;
  pline.addVertex(0.0, 0.0, 0.0);
  double length = cavc::getPathLength(pline);
  EXPECT_TRUE(approxEqual(length, 0.0));
}

TEST(polyline, get_path_length_line_segments) {
  Polyline pline;
  pline.addVertex(0.0, 0.0, 0.0);
  pline.addVertex(3.0, 0.0, 0.0);
  pline.addVertex(3.0, 4.0, 0.0);

  double length = cavc::getPathLength(pline);
  EXPECT_TRUE(approxEqual(length, 7.0)); // 3 + 4
}

TEST(polyline, get_path_length_with_arc) {
  auto vertices = CaseBuilder::quarterArcCase();
  Polyline pline = createFromVertices(vertices, false);

  double length = cavc::getPathLength(pline);
  // Quarter arc with radius 1 should have length π/2
  EXPECT_TRUE(approxEqual(length, cavc::utils::pi<double>() / 2.0, 1e-6));
}

// Winding Number Tests

TEST(polyline, get_winding_number_open_polyline) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices, false);

  int winding = cavc::getWindingNumber(pline, Vector2{0.5, 0.5});
  EXPECT_EQ(winding, 0); // Open polylines always return 0
}

TEST(polyline, get_winding_number_point_inside_rectangle) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices);

  int winding = cavc::getWindingNumber(pline, Vector2{0.5, 0.5});
  EXPECT_EQ(std::abs(winding), 1); // Point inside should have non-zero winding
}

TEST(polyline, get_winding_number_point_outside_rectangle) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices);

  int winding = cavc::getWindingNumber(pline, Vector2{2.0, 2.0});
  EXPECT_EQ(winding, 0); // Point outside should have zero winding
}

TEST(polyline, get_winding_number_circle) {
  auto vertices = CaseBuilder::positiveCircle();
  Polyline pline = createFromVertices(vertices);

  // Point inside circle
  int windingInside = cavc::getWindingNumber(pline, Vector2{5.0, 0.0});
  EXPECT_EQ(std::abs(windingInside), 1);

  // Point outside circle
  int windingOutside = cavc::getWindingNumber(pline, Vector2{20.0, 0.0});
  EXPECT_EQ(windingOutside, 0);
}

// Closest Point Tests

TEST(polyline, closest_point_single_vertex) {
  Polyline pline;
  pline.addVertex(1.0, 2.0, 0.0);

  cavc::ClosestPoint<double> cp(pline, Vector2{3.0, 4.0});
  EXPECT_EQ(cp.index(), 0);
  EXPECT_TRUE(approxEqual(cp.point(), Vector2{1.0, 2.0}));
  EXPECT_TRUE(approxEqual(cp.distance(), std::sqrt(8.0))); // Distance from (3,4) to (1,2)
}

TEST(polyline, closest_point_line_segment) {
  Polyline pline;
  pline.addVertex(0.0, 0.0, 0.0);
  pline.addVertex(2.0, 0.0, 0.0);

  cavc::ClosestPoint<double> cp(pline, Vector2{1.0, 1.0});
  EXPECT_EQ(cp.index(), 0);
  EXPECT_TRUE(approxEqual(cp.point(), Vector2{1.0, 0.0})); // Closest point on segment
  EXPECT_TRUE(approxEqual(cp.distance(), 1.0));
}

TEST(polyline, closest_point_to_vertex) {
  Polyline pline;
  pline.addVertex(0.0, 0.0, 0.0);
  pline.addVertex(2.0, 0.0, 0.0);
  pline.addVertex(2.0, 2.0, 0.0);

  cavc::ClosestPoint<double> cp(pline, Vector2{2.0, 0.0});
  EXPECT_TRUE(approxEqual(cp.point(), Vector2{2.0, 0.0}));
  EXPECT_TRUE(approxEqual(cp.distance(), 0.0));
}

// Transformation Tests

TEST(polyline, scale_polyline) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices);

  cavc::scalePolyline(pline, 2.0);

  EXPECT_TRUE(approxEqual(pline[0].pos(), Vector2{0.0, 0.0}));
  EXPECT_TRUE(approxEqual(pline[1].pos(), Vector2{2.0, 0.0}));
  EXPECT_TRUE(approxEqual(pline[2].pos(), Vector2{2.0, 2.0}));
  EXPECT_TRUE(approxEqual(pline[3].pos(), Vector2{0.0, 2.0}));

  // Bulge values should remain unchanged
  for (const auto &vertex : pline.vertexes()) {
    EXPECT_TRUE(approxEqual(vertex.bulge(), 0.0));
  }
}

TEST(polyline, translate_polyline) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices);

  Vector2 offset{3.0, 4.0};
  cavc::translatePolyline(pline, offset);

  EXPECT_TRUE(approxEqual(pline[0].pos(), Vector2{3.0, 4.0}));
  EXPECT_TRUE(approxEqual(pline[1].pos(), Vector2{4.0, 4.0}));
  EXPECT_TRUE(approxEqual(pline[2].pos(), Vector2{4.0, 5.0}));
  EXPECT_TRUE(approxEqual(pline[3].pos(), Vector2{3.0, 5.0}));
}

// Direction and Modification Tests

TEST(polyline, invert_direction) {
  Polyline pline;
  pline.addVertex(0.0, 0.0, 0.5);
  pline.addVertex(1.0, 0.0, -0.3);
  pline.addVertex(1.0, 1.0, 0.0);

  cavc::invertDirection(pline);

  // Vertices should be reversed
  EXPECT_TRUE(approxEqual(pline[0].pos(), Vector2{1.0, 1.0}));
  EXPECT_TRUE(approxEqual(pline[1].pos(), Vector2{1.0, 0.0}));
  EXPECT_TRUE(approxEqual(pline[2].pos(), Vector2{0.0, 0.0}));

  // Bulges should be negated and shifted
  EXPECT_TRUE(approxEqual(pline[0].bulge(), 0.3));  // was -0.3, negated
  EXPECT_TRUE(approxEqual(pline[1].bulge(), -0.5)); // was 0.5, negated
  EXPECT_TRUE(approxEqual(pline[2].bulge(), 0.0));  // was 0.0, negated
}

TEST(polyline, prune_singularities) {
  Polyline pline;
  pline.addVertex(0.0, 0.0, 0.0);
  pline.addVertex(0.0, 0.0, 0.5); // Duplicate position
  pline.addVertex(1.0, 0.0, 0.0);
  pline.addVertex(1.0, 1.0, 0.0);
  pline.addVertex(1.0, 1.0, 0.3); // Another duplicate

  Polyline pruned = cavc::pruneSingularities(pline, 1e-9);

  EXPECT_EQ(pruned.size(), 3);
  EXPECT_TRUE(approxEqual(pruned[0].pos(), Vector2{0.0, 0.0}));
  EXPECT_TRUE(approxEqual(pruned[0].bulge(), 0.5)); // Should keep last bulge value
  EXPECT_TRUE(approxEqual(pruned[1].pos(), Vector2{1.0, 0.0}));
  EXPECT_TRUE(approxEqual(pruned[2].pos(), Vector2{1.0, 1.0}));
  EXPECT_TRUE(approxEqual(pruned[2].bulge(), 0.3)); // Should keep last bulge value
}

TEST(polyline, prune_singularities_closed_polyline) {
  Polyline pline;
  pline.isClosed() = true;
  pline.addVertex(0.0, 0.0, 0.0);
  pline.addVertex(1.0, 0.0, 0.0);
  pline.addVertex(1.0, 1.0, 0.0);
  pline.addVertex(0.0, 0.0, 0.5); // Last vertex same as first

  Polyline pruned = cavc::pruneSingularities(pline, 1e-9);

  EXPECT_EQ(pruned.size(), 3); // Should remove duplicate last vertex
  EXPECT_TRUE(pruned.isClosed());
}

// Arc Conversion Tests

TEST(polyline, convert_arcs_to_lines_no_arcs) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices);

  Polyline converted = cavc::convertArcsToLines(pline, 0.1);

  EXPECT_EQ(converted.isClosed(), pline.isClosed());

  // Since there are no arcs, the result should have reasonable size
  EXPECT_GT(converted.size(), 0);

  // All vertices should have zero bulge after conversion
  for (const auto &vertex : converted.vertexes()) {
    EXPECT_TRUE(approxEqual(vertex.bulge(), 0.0));
  }
}

TEST(polyline, convert_arcs_to_lines_with_arc) {
  auto vertices = CaseBuilder::quarterArcCase();
  Polyline pline = createFromVertices(vertices, false);

  Polyline converted = cavc::convertArcsToLines(pline, 0.01);

  // Should have more vertices than original due to arc tessellation
  EXPECT_GT(converted.size(), pline.size());
  EXPECT_EQ(converted.isClosed(), pline.isClosed());

  // All vertices should have zero bulge
  for (const auto &vertex : converted.vertexes()) {
    EXPECT_TRUE(approxEqual(vertex.bulge(), 0.0));
  }

  // First and last vertices should be the same
  EXPECT_TRUE(approxEqual(converted[0].pos(), pline[0].pos()));
  EXPECT_TRUE(approxEqual(converted.lastVertex().pos(), pline.lastVertex().pos()));
}

// Spatial Index Tests

TEST(polyline, create_spatial_index) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices);

  auto spatialIndex = cavc::createApproxSpatialIndex(pline);

  // Just test that it doesn't crash and creates something
  // Check that we can query the spatial index
  std::vector<std::size_t> results;
  EXPECT_NO_THROW(spatialIndex.query(0.0, 0.0, 1.0, 1.0, results));
  EXPECT_GT(results.size(), 0); // Should find some segments
}

// Segment Iteration Tests

TEST(polyline, visit_segment_indices_open) {
  Polyline pline;
  pline.addVertex(0.0, 0.0, 0.0);
  pline.addVertex(1.0, 0.0, 0.0);
  pline.addVertex(1.0, 1.0, 0.0);
  pline.isClosed() = false;

  std::vector<std::pair<size_t, size_t>> segments;
  pline.visitSegIndices([&](size_t i, size_t j) {
    segments.emplace_back(i, j);
    return true;
  });

  EXPECT_EQ(segments.size(), 2);
  EXPECT_EQ(segments[0], std::make_pair(size_t(0), size_t(1)));
  EXPECT_EQ(segments[1], std::make_pair(size_t(1), size_t(2)));
}

TEST(polyline, visit_segment_indices_closed) {
  Polyline pline;
  pline.addVertex(0.0, 0.0, 0.0);
  pline.addVertex(1.0, 0.0, 0.0);
  pline.addVertex(1.0, 1.0, 0.0);
  pline.isClosed() = true;

  std::vector<std::pair<size_t, size_t>> segments;
  pline.visitSegIndices([&](size_t i, size_t j) {
    segments.emplace_back(i, j);
    return true;
  });

  EXPECT_EQ(segments.size(), 3);
  EXPECT_EQ(segments[0], std::make_pair(size_t(2), size_t(0))); // Last to first
  EXPECT_EQ(segments[1], std::make_pair(size_t(0), size_t(1)));
  EXPECT_EQ(segments[2], std::make_pair(size_t(1), size_t(2)));
}

TEST(polyline, visit_segment_indices_early_termination) {
  auto vertices = CaseBuilder::simpleRectangle();
  Polyline pline = createFromVertices(vertices);

  std::vector<std::pair<size_t, size_t>> segments;
  pline.visitSegIndices([&](size_t i, size_t j) {
    segments.emplace_back(i, j);
    return segments.size() < 2; // Stop after 2 segments
  });

  EXPECT_EQ(segments.size(), 2);
}

// Edge Cases

TEST(polyline, empty_polyline_basic_operations) {
  Polyline pline;

  // Basic query operations should not crash
  EXPECT_NO_THROW(cavc::getExtents(pline));
  EXPECT_NO_THROW(cavc::getArea(pline));
  EXPECT_NO_THROW(cavc::getPathLength(pline));
  EXPECT_NO_THROW(cavc::getWindingNumber(pline, Vector2{0.0, 0.0}));
}

TEST(polyline, empty_polyline_transformation_operations) {
  Polyline pline;

  // Transformation operations should not crash
  EXPECT_NO_THROW(cavc::scalePolyline(pline, 2.0));
  EXPECT_NO_THROW(cavc::translatePolyline(pline, Vector2{1.0, 1.0}));
}

TEST(polyline, empty_polyline_invert_direction) {
  Polyline pline;
  EXPECT_NO_THROW(cavc::invertDirection(pline));
}

TEST(polyline, empty_polyline_prune_singularities) {
  Polyline pline;
  EXPECT_NO_THROW(cavc::pruneSingularities(pline, 1e-9));
}

TEST(polyline, empty_polyline_convert_arcs_fixed) {
  Polyline pline;
  EXPECT_NO_THROW(cavc::convertArcsToLines(pline, 0.1));

  Polyline converted = cavc::convertArcsToLines(pline, 0.1);
  EXPECT_EQ(converted.size(), 0);
  EXPECT_EQ(converted.isClosed(), pline.isClosed());
}

TEST(polyline, single_vertex_operations) {
  Polyline pline;
  pline.addVertex(1.0, 2.0, 0.5);

  // These should work with single vertex
  EXPECT_NO_THROW(cavc::getExtents(pline));
  EXPECT_TRUE(approxEqual(cavc::getArea(pline), 0.0));
  EXPECT_TRUE(approxEqual(cavc::getPathLength(pline), 0.0));
  EXPECT_EQ(cavc::getWindingNumber(pline, Vector2{0.0, 0.0}), 0);
}

TEST(polyline, large_coordinates) {
  Polyline pline;
  pline.addVertex(1e6, 1e6, 0.0);
  pline.addVertex(1e6 + 1, 1e6, 0.0);
  pline.addVertex(1e6 + 1, 1e6 + 1, 0.0);
  pline.addVertex(1e6, 1e6 + 1, 0.0);
  pline.isClosed() = true;

  double area = cavc::getArea(pline);
  EXPECT_TRUE(approxEqual(area, 1.0, 1e-10));

  double length = cavc::getPathLength(pline);
  EXPECT_TRUE(approxEqual(length, 4.0, 1e-10));
}

TEST(polyline, very_small_coordinates) {
  Polyline pline;
  pline.addVertex(1e-6, 1e-6, 0.0);
  pline.addVertex(2e-6, 1e-6, 0.0);
  pline.addVertex(2e-6, 2e-6, 0.0);
  pline.addVertex(1e-6, 2e-6, 0.0);
  pline.isClosed() = true;

  double area = cavc::getArea(pline);
  EXPECT_TRUE(approxEqual(area, 1e-12, 1e-15));

  double length = cavc::getPathLength(pline);
  EXPECT_TRUE(approxEqual(length, 4e-6, 1e-12));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
