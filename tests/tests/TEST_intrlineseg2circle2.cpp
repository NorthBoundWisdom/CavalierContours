#include <algorithm>
#include <cavc/intrlineseg2circle2.hpp>
#include <cavc/mathutils.hpp>
#include <cavc/vector2.hpp>
#include <gtest/gtest.h>

using Vector2 = cavccpp::Vector2<double>;
using IntrLineSeg2Circle2Result = cavccpp::IntrLineSeg2Circle2Result<double>;

namespace {

// Helper function to check if two doubles are approximately equal
constexpr double EPSILON = 1e-9;
bool approxEqual(double a, double b, double epsilon = EPSILON) { return std::abs(a - b) < epsilon; }

// Helper function to get point on segment using parametric t
Vector2 getPointOnSegment(const Vector2 &p0, const Vector2 &p1, double t) {
  return Vector2{p0.x() + t * (p1.x() - p0.x()), p0.y() + t * (p1.y() - p0.y())};
}

// Helper function to check if point is on circle
bool isPointOnCircle(const Vector2 &point, const Vector2 &center, double radius,
                     double epsilon = EPSILON) {
  double distSq = (point.x() - center.x()) * (point.x() - center.x()) +
                  (point.y() - center.y()) * (point.y() - center.y());
  return approxEqual(distSq, radius * radius, epsilon);
}

} // namespace

TEST(intrlineseg2circle2, no_intersection_segment_outside_circle) {
  // Horizontal line segment completely outside circle
  // The infinite line still intersects the circle, but not within the segment bounds
  Vector2 p0{3.0, 0.0};
  Vector2 p1{5.0, 0.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  // The infinite line y=0 intersects the circle at x=±1, but t values will be negative
  EXPECT_EQ(result.numIntersects, 2);
  // Both intersections should be outside the segment (t < 0)
  EXPECT_TRUE(result.t0 < 0.0 && result.t1 < 0.0);
}

TEST(intrlineseg2circle2, no_intersection_segment_inside_circle) {
  // Line segment completely inside circle
  // The infinite line still intersects the circle outside the segment
  Vector2 p0{0.2, 0.0};
  Vector2 p1{0.4, 0.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  // The infinite line y=0 intersects the circle at x=±1
  EXPECT_EQ(result.numIntersects, 2);
  // One intersection before segment (t < 0) and one after (t > 1)
  EXPECT_TRUE((result.t0 < 0.0 && result.t1 > 1.0) || (result.t1 < 0.0 && result.t0 > 1.0));
}

TEST(intrlineseg2circle2, one_intersection_tangent_horizontal) {
  // Horizontal line tangent to circle
  Vector2 p0{-2.0, 1.0};
  Vector2 p1{2.0, 1.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 1);
  EXPECT_TRUE(approxEqual(result.t0, 0.5)); // Tangent point at middle of segment

  Vector2 intersectionPoint = getPointOnSegment(p0, p1, result.t0);
  EXPECT_TRUE(isPointOnCircle(intersectionPoint, center, radius));
}

TEST(intrlineseg2circle2, one_intersection_tangent_vertical) {
  // Vertical line tangent to circle
  Vector2 p0{1.0, -2.0};
  Vector2 p1{1.0, 2.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 1);
  EXPECT_TRUE(approxEqual(result.t0, 0.5)); // Tangent point at middle of segment

  Vector2 intersectionPoint = getPointOnSegment(p0, p1, result.t0);
  EXPECT_TRUE(isPointOnCircle(intersectionPoint, center, radius));
}

TEST(intrlineseg2circle2, two_intersections_horizontal_through_center) {
  // Horizontal line passing through circle center
  Vector2 p0{-2.0, 0.0};
  Vector2 p1{2.0, 0.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);

  // Intersection points should be at t = 0.25 and t = 0.75 (at x = -1 and x = 1)
  EXPECT_TRUE(approxEqual(result.t0, 0.25) || approxEqual(result.t0, 0.75));
  EXPECT_TRUE(approxEqual(result.t1, 0.25) || approxEqual(result.t1, 0.75));
  EXPECT_FALSE(approxEqual(result.t0, result.t1));

  Vector2 intersection1 = getPointOnSegment(p0, p1, result.t0);
  Vector2 intersection2 = getPointOnSegment(p0, p1, result.t1);
  EXPECT_TRUE(isPointOnCircle(intersection1, center, radius));
  EXPECT_TRUE(isPointOnCircle(intersection2, center, radius));
}

TEST(intrlineseg2circle2, two_intersections_vertical_through_center) {
  // Vertical line passing through circle center
  Vector2 p0{0.0, -2.0};
  Vector2 p1{0.0, 2.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);

  Vector2 intersection1 = getPointOnSegment(p0, p1, result.t0);
  Vector2 intersection2 = getPointOnSegment(p0, p1, result.t1);
  EXPECT_TRUE(isPointOnCircle(intersection1, center, radius));
  EXPECT_TRUE(isPointOnCircle(intersection2, center, radius));
}

TEST(intrlineseg2circle2, two_intersections_diagonal) {
  // Diagonal line passing through circle
  Vector2 p0{-2.0, -2.0};
  Vector2 p1{2.0, 2.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);

  Vector2 intersection1 = getPointOnSegment(p0, p1, result.t0);
  Vector2 intersection2 = getPointOnSegment(p0, p1, result.t1);
  EXPECT_TRUE(isPointOnCircle(intersection1, center, radius));
  EXPECT_TRUE(isPointOnCircle(intersection2, center, radius));
}

TEST(intrlineseg2circle2, one_intersection_segment_starts_on_circle) {
  // Segment starts on circle boundary
  Vector2 p0{1.0, 0.0}; // On circle
  Vector2 p1{3.0, 0.0}; // Outside circle
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  // The infinite line intersects at x=±1, but we start at x=1
  EXPECT_EQ(result.numIntersects, 2);
  // One intersection at start (t=0) and one before the segment (t<0)
  EXPECT_TRUE((approxEqual(result.t0, 0.0) && result.t1 < 0.0) ||
              (approxEqual(result.t1, 0.0) && result.t0 < 0.0));
}

TEST(intrlineseg2circle2, one_intersection_segment_ends_on_circle) {
  // Segment ends on circle boundary
  Vector2 p0{3.0, 0.0}; // Outside circle
  Vector2 p1{1.0, 0.0}; // On circle
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  // The infinite line intersects at x=±1, we end at x=1
  EXPECT_EQ(result.numIntersects, 2);
  // One intersection at end (t=1) and one after the segment (t>1)
  EXPECT_TRUE((approxEqual(result.t0, 1.0) && result.t1 > 1.0) ||
              (approxEqual(result.t1, 1.0) && result.t0 > 1.0));
}

TEST(intrlineseg2circle2, two_intersections_segment_crosses_circle) {
  // Segment crosses circle from outside to outside
  Vector2 p0{-1.5, 0.0};
  Vector2 p1{1.5, 0.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);

  // Both t values should be between 0 and 1
  EXPECT_TRUE(result.t0 >= 0.0 && result.t0 <= 1.0);
  EXPECT_TRUE(result.t1 >= 0.0 && result.t1 <= 1.0);

  Vector2 intersection1 = getPointOnSegment(p0, p1, result.t0);
  Vector2 intersection2 = getPointOnSegment(p0, p1, result.t1);
  EXPECT_TRUE(isPointOnCircle(intersection1, center, radius));
  EXPECT_TRUE(isPointOnCircle(intersection2, center, radius));
}

TEST(intrlineseg2circle2, zero_length_segment_on_circle) {
  // Point (zero length segment) on circle
  Vector2 p0{1.0, 0.0};
  Vector2 p1{1.0, 0.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 1);
  EXPECT_TRUE(approxEqual(result.t0, 0.0));
}

TEST(intrlineseg2circle2, zero_length_segment_off_circle) {
  // Point (zero length segment) not on circle
  Vector2 p0{2.0, 0.0};
  Vector2 p1{2.0, 0.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 0);
}

TEST(intrlineseg2circle2, segment_chord_of_circle) {
  // Segment is a chord of the circle (both endpoints on circle)
  Vector2 p0{1.0, 0.0};
  Vector2 p1{0.0, 1.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);
  EXPECT_TRUE(approxEqual(result.t0, 0.0) || approxEqual(result.t0, 1.0));
  EXPECT_TRUE(approxEqual(result.t1, 0.0) || approxEqual(result.t1, 1.0));
}

TEST(intrlineseg2circle2, circle_not_at_origin) {
  // Circle centered at non-origin point
  Vector2 p0{1.0, 1.0};
  Vector2 p1{3.0, 1.0};
  double radius = 1.0;
  Vector2 center{2.0, 1.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);

  Vector2 intersection1 = getPointOnSegment(p0, p1, result.t0);
  Vector2 intersection2 = getPointOnSegment(p0, p1, result.t1);
  EXPECT_TRUE(isPointOnCircle(intersection1, center, radius));
  EXPECT_TRUE(isPointOnCircle(intersection2, center, radius));
}

TEST(intrlineseg2circle2, large_circle_small_segment) {
  // Large circle with small segment
  Vector2 p0{99.9, 0.0};
  Vector2 p1{100.1, 0.0};
  double radius = 100.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);
}

TEST(intrlineseg2circle2, very_small_circle) {
  // Very small circle
  Vector2 p0{-0.001, 0.0};
  Vector2 p1{0.001, 0.0};
  double radius = 0.0005;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  // Should detect tangency or very close intersections
  EXPECT_TRUE(result.numIntersects == 1 || result.numIntersects == 2);
}

TEST(intrlineseg2circle2, negative_coordinates) {
  // Test with negative coordinates
  Vector2 p0{-3.0, -1.0};
  Vector2 p1{-1.0, -1.0};
  double radius = 1.0;
  Vector2 center{-2.0, -1.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);

  Vector2 intersection1 = getPointOnSegment(p0, p1, result.t0);
  Vector2 intersection2 = getPointOnSegment(p0, p1, result.t1);
  EXPECT_TRUE(isPointOnCircle(intersection1, center, radius));
  EXPECT_TRUE(isPointOnCircle(intersection2, center, radius));
}

TEST(intrlineseg2circle2, extended_intersection_before_segment) {
  // Intersection would occur if line extended backward from p0
  Vector2 p0{0.5, 0.0}; // Inside circle
  Vector2 p1{2.0, 0.0}; // Outside circle
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);

  // Should have one intersection within segment (where it exits) and one before
  bool hasIntersectionInSegment =
      (result.t0 >= 0.0 && result.t0 <= 1.0) || (result.t1 >= 0.0 && result.t1 <= 1.0);
  bool hasIntersectionBefore = (result.t0 < 0.0) || (result.t1 < 0.0);
  EXPECT_TRUE(hasIntersectionInSegment && hasIntersectionBefore);
}

TEST(intrlineseg2circle2, extended_intersection_after_segment) {
  // Intersection would occur if line extended forward from p1
  Vector2 p0{-2.0, 0.0}; // Outside circle
  Vector2 p1{-0.5, 0.0}; // Inside circle
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);

  // Should have one intersection within segment (where it enters) and one after
  bool hasIntersectionInSegment =
      (result.t0 >= 0.0 && result.t0 <= 1.0) || (result.t1 >= 0.0 && result.t1 <= 1.0);
  bool hasIntersectionAfter = (result.t0 > 1.0) || (result.t1 > 1.0);
  EXPECT_TRUE(hasIntersectionInSegment && hasIntersectionAfter);
}

TEST(intrlineseg2circle2, precision_edge_case_near_tangent) {
  // Test near tangent condition with precision boundary
  double y = 1.0 - cavccpp::utils::realThreshold<double>() / 2.0;
  Vector2 p0{-2.0, y};
  Vector2 p1{2.0, y};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  // Should detect as either tangent (1 intersection) or very close intersections (2)
  EXPECT_TRUE(result.numIntersects == 1 || result.numIntersects == 2);
}

TEST(intrlineseg2circle2, zero_radius_circle) {
  // Point circle (zero radius)
  Vector2 p0{-1.0, 0.0};
  Vector2 p1{1.0, 0.0};
  double radius = 0.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 1);
  EXPECT_TRUE(approxEqual(result.t0, 0.5)); // Point at center of segment
}

TEST(intrlineseg2circle2, line_misses_circle_entirely) {
  // Line that doesn't intersect circle at all
  Vector2 p0{0.0, 2.0};
  Vector2 p1{1.0, 2.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 0);
}

TEST(intrlineseg2circle2, segment_entirely_before_intersections) {
  // Segment ends before any intersections would occur
  Vector2 p0{-5.0, 0.0};
  Vector2 p1{-3.0, 0.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);
  // Both intersections should be after the segment (t > 1)
  EXPECT_TRUE(result.t0 > 1.0 && result.t1 > 1.0);
}

TEST(intrlineseg2circle2, segment_entirely_after_intersections) {
  // Segment starts after all intersections
  Vector2 p0{3.0, 0.0};
  Vector2 p1{5.0, 0.0};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);
  // Both intersections should be before the segment (t < 0)
  EXPECT_TRUE(result.t0 < 0.0 && result.t1 < 0.0);
}

TEST(intrlineseg2circle2, verify_parametric_values_horizontal) {
  // Test specific parametric values for horizontal line through circle
  Vector2 p0{-3.0, 0.0};
  Vector2 p1{3.0, 0.0}; // 6 unit long segment
  double radius = 2.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);

  // Intersections at x = ±2, so t = (x - (-3))/6 = (x + 3)/6
  // For x = -2: t = 1/6 ≈ 0.1667
  // For x = 2: t = 5/6 ≈ 0.8333
  std::vector<double> tValues = {result.t0, result.t1};
  std::sort(tValues.begin(), tValues.end());

  EXPECT_TRUE(approxEqual(tValues[0], 1.0 / 6.0, 1e-6));
  EXPECT_TRUE(approxEqual(tValues[1], 5.0 / 6.0, 1e-6));
}

TEST(intrlineseg2circle2, verify_parametric_values_vertical) {
  // Test specific parametric values for vertical line through circle
  Vector2 p0{0.0, -3.0};
  Vector2 p1{0.0, 3.0}; // 6 unit long segment
  double radius = 2.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);

  // Intersections at y = ±2, so t = (y - (-3))/6 = (y + 3)/6
  // For y = -2: t = 1/6 ≈ 0.1667
  // For y = 2: t = 5/6 ≈ 0.8333
  std::vector<double> tValues = {result.t0, result.t1};
  std::sort(tValues.begin(), tValues.end());

  EXPECT_TRUE(approxEqual(tValues[0], 1.0 / 6.0, 1e-6));
  EXPECT_TRUE(approxEqual(tValues[1], 5.0 / 6.0, 1e-6));
}

TEST(intrlineseg2circle2, tangent_t_outside_segment) {
  // Tangent line where the tangent point is outside the segment
  Vector2 p0{2.0, 1.0}; // Start beyond tangent point
  Vector2 p1{4.0, 1.0}; // End even further
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 1);
  // Tangent point at (0,1) should give t = (0-2)/(4-2) = -1
  EXPECT_TRUE(approxEqual(result.t0, -1.0));
}

TEST(intrlineseg2circle2, intersection_points_validation) {
  // Verify that calculated intersection points actually lie on the circle
  Vector2 p0{-2.0, 0.5};
  Vector2 p1{2.0, 0.5};
  double radius = 1.0;
  Vector2 center{0.0, 0.0};

  auto result = cavccpp::intrLineSeg2Circle2(p0, p1, radius, center);
  EXPECT_EQ(result.numIntersects, 2);

  // Calculate actual intersection points and verify they're on the circle
  Vector2 intersection1 = getPointOnSegment(p0, p1, result.t0);
  Vector2 intersection2 = getPointOnSegment(p0, p1, result.t1);

  EXPECT_TRUE(isPointOnCircle(intersection1, center, radius));
  EXPECT_TRUE(isPointOnCircle(intersection2, center, radius));

  // Both points should have y = 0.5 (on the line)
  EXPECT_TRUE(approxEqual(intersection1.y(), 0.5));
  EXPECT_TRUE(approxEqual(intersection2.y(), 0.5));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
