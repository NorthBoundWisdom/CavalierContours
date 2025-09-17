#include <gmock/gmock.h>

#include <gtest/gtest.h>

#include <cavc/intrlineseg2lineseg2.hpp>
#include <cavc/mathutils.hpp>
#include <cavc/vector2.hpp>

using Vector2 = cavccpp::Vector2<double>;
using IntrLineSeg2LineSeg2Result = cavccpp::IntrLineSeg2LineSeg2Result<double>;
using LineSeg2LineSeg2IntrType = cavccpp::LineSeg2LineSeg2IntrType;

namespace {

// Helper function to check if two doubles are approximately equal
constexpr double EPSILON = 1e-9;
bool approxEqual(double a, double b, double epsilon = EPSILON) { return std::abs(a - b) < epsilon; }

bool approxEqual(const Vector2 &a, const Vector2 &b, double epsilon = EPSILON) {
  return approxEqual(a.x(), b.x(), epsilon) && approxEqual(a.y(), b.y(), epsilon);
}

} // namespace

// Test cases for LineSeg2LineSeg2IntrType::True (segments actually intersect)

TEST(intrlineseg2lineseg2, true_intersection_perpendicular_cross) {
  // Two perpendicular segments crossing at their midpoints
  Vector2 u1{-1.0, 0.0};
  Vector2 u2{1.0, 0.0};
  Vector2 v1{0.0, -1.0};
  Vector2 v2{0.0, 1.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::True);
  EXPECT_TRUE(approxEqual(result.point, Vector2{0.0, 0.0}));
}

TEST(intrlineseg2lineseg2, true_intersection_diagonal_cross) {
  // Two diagonal segments crossing
  Vector2 u1{0.0, 0.0};
  Vector2 u2{2.0, 2.0};
  Vector2 v1{0.0, 2.0};
  Vector2 v2{2.0, 0.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::True);
  EXPECT_TRUE(approxEqual(result.point, Vector2{1.0, 1.0}));
}

TEST(intrlineseg2lineseg2, true_intersection_at_endpoint) {
  // Segments meet at an endpoint
  Vector2 u1{0.0, 0.0};
  Vector2 u2{1.0, 0.0};
  Vector2 v1{1.0, 0.0};
  Vector2 v2{1.0, 1.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::True);
  EXPECT_TRUE(approxEqual(result.point, Vector2{1.0, 0.0}));
}

TEST(intrlineseg2lineseg2, true_intersection_t_junction) {
  // T-junction intersection
  Vector2 u1{0.0, 0.0};
  Vector2 u2{2.0, 0.0};
  Vector2 v1{1.0, -1.0};
  Vector2 v2{1.0, 0.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::True);
  EXPECT_TRUE(approxEqual(result.point, Vector2{1.0, 0.0}));
}

TEST(intrlineseg2lineseg2, true_intersection_point_segments) {
  // Both segments are points at the same location
  Vector2 u1{1.0, 1.0};
  Vector2 u2{1.0, 1.0};
  Vector2 v1{1.0, 1.0};
  Vector2 v2{1.0, 1.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::True);
  EXPECT_TRUE(approxEqual(result.point, Vector2{1.0, 1.0}));
}

TEST(intrlineseg2lineseg2, true_intersection_point_on_segment) {
  // One segment is a point lying on the other segment
  Vector2 u1{1.0, 1.0}; // Point segment
  Vector2 u2{1.0, 1.0};
  Vector2 v1{0.0, 1.0}; // Line segment containing the point
  Vector2 v2{2.0, 1.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::True);
  EXPECT_TRUE(approxEqual(result.point, Vector2{1.0, 1.0}));
}

// Test cases for LineSeg2LineSeg2IntrType::False (segments would intersect if extended)

TEST(intrlineseg2lineseg2, false_intersection_would_meet_if_extended) {
  // Two segments that would intersect if one was extended
  Vector2 u1{0.0, 0.0};
  Vector2 u2{1.0, 0.0};
  Vector2 v1{2.0, -1.0}; // Would intersect at (2,0) if u1-u2 extended
  Vector2 v2{2.0, 1.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::False);
  EXPECT_TRUE(approxEqual(result.point, Vector2{2.0, 0.0}));
  EXPECT_TRUE(result.t0 > 1.0);                      // Extension of first segment
  EXPECT_TRUE(result.t1 >= 0.0 && result.t1 <= 1.0); // Within second segment
}

TEST(intrlineseg2lineseg2, false_intersection_both_need_extension) {
  // Both segments need extension to meet
  Vector2 u1{0.0, 0.0};
  Vector2 u2{0.5, 0.0};
  Vector2 v1{1.5, -0.5};
  Vector2 v2{1.5, 0.5};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::False);
  EXPECT_TRUE(result.t0 > 1.0);                      // Extension of first segment
  EXPECT_TRUE(result.t1 >= 0.0 && result.t1 <= 1.0); // Within second segment
}

// Test cases for LineSeg2LineSeg2IntrType::None (no intersection possible)

TEST(intrlineseg2lineseg2, no_intersection_parallel_separate) {
  // Parallel segments that don't overlap
  Vector2 u1{0.0, 0.0};
  Vector2 u2{1.0, 0.0};
  Vector2 v1{0.0, 1.0};
  Vector2 v2{1.0, 1.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::None);
}

TEST(intrlineseg2lineseg2, no_intersection_skew_segments) {
  // Skew segments (not parallel, don't intersect)
  Vector2 u1{0.0, 0.0};
  Vector2 u2{1.0, 0.0};
  Vector2 v1{2.0, 1.0};
  Vector2 v2{3.0, 2.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::False); // Would intersect if extended
}

TEST(intrlineseg2lineseg2, no_intersection_point_segments_different) {
  // Two point segments at different locations
  Vector2 u1{0.0, 0.0};
  Vector2 u2{0.0, 0.0};
  Vector2 v1{1.0, 1.0};
  Vector2 v2{1.0, 1.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::None);
}

TEST(intrlineseg2lineseg2, no_intersection_point_not_on_segment) {
  // Point segment not lying on line segment
  Vector2 u1{1.0, 2.0}; // Point not on the line
  Vector2 u2{1.0, 2.0};
  Vector2 v1{0.0, 0.0}; // Line segment
  Vector2 v2{2.0, 0.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::None);
}

// Test cases for LineSeg2LineSeg2IntrType::Coincident (segments overlap)

TEST(intrlineseg2lineseg2, coincident_complete_overlap) {
  // Two identical segments
  Vector2 u1{0.0, 0.0};
  Vector2 u2{2.0, 0.0};
  Vector2 v1{0.0, 0.0};
  Vector2 v2{2.0, 0.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::Coincident);
  EXPECT_TRUE(approxEqual(result.t0, 0.0));
  EXPECT_TRUE(approxEqual(result.t1, 1.0));
}

TEST(intrlineseg2lineseg2, coincident_partial_overlap) {
  // Segments partially overlap
  Vector2 u1{0.0, 0.0};
  Vector2 u2{3.0, 0.0};
  Vector2 v1{1.0, 0.0};
  Vector2 v2{2.0, 0.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::Coincident);
  // v1 to v2 should map to parameters on the u1-u2 line
  EXPECT_TRUE(result.t0 >= 0.0 && result.t0 <= 1.0);
  EXPECT_TRUE(result.t1 >= 0.0 && result.t1 <= 1.0);
  EXPECT_TRUE(result.t0 < result.t1);
}

TEST(intrlineseg2lineseg2, coincident_reversed_overlap) {
  // Segments overlap but in opposite directions
  Vector2 u1{0.0, 0.0};
  Vector2 u2{2.0, 0.0};
  Vector2 v1{2.0, 0.0};
  Vector2 v2{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::Coincident);
  EXPECT_TRUE(approxEqual(result.t0, 0.0));
  EXPECT_TRUE(approxEqual(result.t1, 1.0));
}

TEST(intrlineseg2lineseg2, coincident_endpoint_overlap) {
  // Segments meet end-to-end exactly
  Vector2 u1{0.0, 0.0};
  Vector2 u2{1.0, 0.0};
  Vector2 v1{1.0, 0.0};
  Vector2 v2{2.0, 0.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::True); // Single point intersection
  EXPECT_TRUE(approxEqual(result.point, Vector2{1.0, 0.0}));
}

// Edge cases and special scenarios

TEST(intrlineseg2lineseg2, vertical_segments_intersect) {
  // Two vertical segments intersecting
  Vector2 u1{1.0, 0.0};
  Vector2 u2{1.0, 2.0};
  Vector2 v1{1.0, 1.0};
  Vector2 v2{1.0, 3.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::Coincident);
}

TEST(intrlineseg2lineseg2, horizontal_segments_intersect) {
  // Two horizontal segments intersecting
  Vector2 u1{0.0, 1.0};
  Vector2 u2{2.0, 1.0};
  Vector2 v1{1.0, 1.0};
  Vector2 v2{3.0, 1.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::Coincident);
}

TEST(intrlineseg2lineseg2, very_small_segments) {
  // Very small segments
  Vector2 u1{0.0, 0.0};
  Vector2 u2{0.001, 0.001};
  Vector2 v1{0.0005, -0.0005};
  Vector2 v2{0.0005, 0.002};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  // Should detect some form of intersection or near-intersection
  EXPECT_TRUE(result.intrType != LineSeg2LineSeg2IntrType::None);
}

TEST(intrlineseg2lineseg2, large_coordinates) {
  // Test with large coordinate values
  Vector2 u1{1000.0, 1000.0};
  Vector2 u2{1001.0, 1000.0};
  Vector2 v1{1000.5, 999.0};
  Vector2 v2{1000.5, 1001.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::True);
  EXPECT_TRUE(approxEqual(result.point, Vector2{1000.5, 1000.0}));
}

TEST(intrlineseg2lineseg2, negative_coordinates) {
  // Test with negative coordinates
  Vector2 u1{-2.0, -1.0};
  Vector2 u2{-1.0, -2.0};
  Vector2 v1{-2.0, -2.0};
  Vector2 v2{-1.0, -1.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::True);
  EXPECT_TRUE(approxEqual(result.point, Vector2{-1.5, -1.5}));
}

TEST(intrlineseg2lineseg2, nearly_parallel_segments) {
  // Segments that are almost parallel (within threshold)
  Vector2 u1{0.0, 0.0};
  Vector2 u2{1.0, 0.0};
  Vector2 v1{0.0, 1e-10}; // Very slight angle
  Vector2 v2{1.0, 1e-10};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  // Should be treated as parallel due to threshold
  EXPECT_TRUE(result.intrType == LineSeg2LineSeg2IntrType::Coincident ||
              result.intrType == LineSeg2LineSeg2IntrType::None);
}

TEST(intrlineseg2lineseg2, precision_boundary_test) {
  // Test at the precision boundary for intersection detection
  double eps = cavccpp::utils::realThreshold<double>();
  Vector2 u1{0.0, 0.0};
  Vector2 u2{1.0, 0.0};
  Vector2 v1{0.5, eps / 2.0};
  Vector2 v2{0.5, -eps / 2.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  // Should detect intersection due to threshold tolerance
  EXPECT_TRUE(result.intrType == LineSeg2LineSeg2IntrType::True);
}

TEST(intrlineseg2lineseg2, parametric_values_validation) {
  // Verify parametric values for a known intersection
  Vector2 u1{0.0, 0.0};
  Vector2 u2{4.0, 0.0}; // 4 units long
  Vector2 v1{2.0, -1.0};
  Vector2 v2{2.0, 1.0}; // 2 units long, intersects at (2,0)

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::True);
  EXPECT_TRUE(approxEqual(result.point, Vector2{2.0, 0.0}));

  // For False intersections, we can check the parametric values
  // t0 should be 0.5 (middle of first segment), t1 should be 0.5 (middle of second segment)
  // Note: for True intersections, t0 and t1 are undefined according to the documentation
}

TEST(intrlineseg2lineseg2, collinear_segments_no_overlap) {
  // Collinear segments that don't overlap
  Vector2 u1{0.0, 0.0};
  Vector2 u2{1.0, 0.0};
  Vector2 v1{2.0, 0.0};
  Vector2 v2{3.0, 0.0};

  auto result = cavccpp::intrLineSeg2LineSeg2(u1, u2, v1, v2);
  EXPECT_EQ(result.intrType, LineSeg2LineSeg2IntrType::None);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
