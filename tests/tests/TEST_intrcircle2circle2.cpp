#include <cavc/intrcircle2circle2.hpp>
#include <cavc/mathutils.hpp>
#include <cavc/vector2.hpp>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using Vector2 = cavccpp::Vector2<double>;
using IntrCircle2Circle2Result = cavccpp::IntrCircle2Circle2Result<double>;
using Circle2Circle2IntrType = cavccpp::Circle2Circle2IntrType;

namespace {

// Helper function to check if two doubles are approximately equal
constexpr double EPSILON = 1e-9;
bool approxEqual(double a, double b, double epsilon = EPSILON) { return std::abs(a - b) < epsilon; }

bool approxEqual(const Vector2 &a, const Vector2 &b, double epsilon = EPSILON) {
  return approxEqual(a.x(), b.x(), epsilon) && approxEqual(a.y(), b.y(), epsilon);
}

} // namespace

TEST(intrcircle2circle2, coincident_circles) {
  // Two identical circles should be coincident
  auto result = cavccpp::intrCircle2Circle2(1.0, Vector2{0.0, 0.0}, 1.0, Vector2{0.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::Coincident);
}

TEST(intrcircle2circle2, no_intersect_far_apart) {
  // Two circles far apart should not intersect
  auto result = cavccpp::intrCircle2Circle2(1.0, Vector2{0.0, 0.0}, 1.0, Vector2{5.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::NoIntersect);
}

TEST(intrcircle2circle2, no_intersect_one_inside_other) {
  // Small circle inside large circle with no intersection
  auto result = cavccpp::intrCircle2Circle2(0.5, Vector2{0.0, 0.0}, 2.0, Vector2{0.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::NoIntersect);
}

TEST(intrcircle2circle2, external_tangent) {
  // Two circles touching externally (one intersection point)
  auto result = cavccpp::intrCircle2Circle2(1.0, Vector2{0.0, 0.0}, 1.0, Vector2{2.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::OneIntersect);
  EXPECT_TRUE(approxEqual(result.point1, Vector2{1.0, 0.0}));
}

TEST(intrcircle2circle2, internal_tangent) {
  // Large circle containing smaller circle, touching internally
  auto result = cavccpp::intrCircle2Circle2(2.0, Vector2{0.0, 0.0}, 1.0, Vector2{1.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::OneIntersect);
  EXPECT_TRUE(approxEqual(result.point1, Vector2{2.0, 0.0}));
}

TEST(intrcircle2circle2, two_intersects_horizontal) {
  // Two circles intersecting at two points (horizontal alignment)
  auto result = cavccpp::intrCircle2Circle2(1.0, Vector2{0.0, 0.0}, 1.0, Vector2{1.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::TwoIntersects);

  // The intersection points should be at (0.5, ±sqrt(3)/2)
  double expectedY = std::sqrt(3.0) / 2.0;
  Vector2 expected1{0.5, expectedY};
  Vector2 expected2{0.5, -expectedY};

  // Check that we get the expected points (order may vary)
  EXPECT_TRUE((approxEqual(result.point1, expected1) && approxEqual(result.point2, expected2)) ||
              (approxEqual(result.point1, expected2) && approxEqual(result.point2, expected1)));
}

TEST(intrcircle2circle2, two_intersects_vertical) {
  // Two circles intersecting at two points (vertical alignment)
  auto result = cavccpp::intrCircle2Circle2(1.0, Vector2{0.0, 0.0}, 1.0, Vector2{0.0, 1.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::TwoIntersects);

  // The intersection points should be at (±sqrt(3)/2, 0.5)
  double expectedX = std::sqrt(3.0) / 2.0;
  Vector2 expected1{expectedX, 0.5};
  Vector2 expected2{-expectedX, 0.5};

  EXPECT_TRUE((approxEqual(result.point1, expected1) && approxEqual(result.point2, expected2)) ||
              (approxEqual(result.point1, expected2) && approxEqual(result.point2, expected1)));
}

TEST(intrcircle2circle2, two_intersects_diagonal) {
  // Two circles intersecting at two points (diagonal alignment)
  auto result = cavccpp::intrCircle2Circle2(1.0, Vector2{0.0, 0.0}, 1.0, Vector2{1.0, 1.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::TwoIntersects);

  // Should have two intersection points
  EXPECT_FALSE(approxEqual(result.point1, result.point2));
}

TEST(intrcircle2circle2, different_radii_intersect) {
  // Circles with different radii intersecting
  auto result = cavccpp::intrCircle2Circle2(1.0, Vector2{0.0, 0.0}, 2.0, Vector2{2.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::TwoIntersects);

  // Should have two distinct intersection points
  EXPECT_FALSE(approxEqual(result.point1, result.point2));
}

TEST(intrcircle2circle2, different_radii_external_tangent) {
  // Circles with different radii touching externally
  auto result = cavccpp::intrCircle2Circle2(1.0, Vector2{0.0, 0.0}, 2.0, Vector2{3.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::OneIntersect);
  EXPECT_TRUE(approxEqual(result.point1, Vector2{1.0, 0.0}));
}

TEST(intrcircle2circle2, different_radii_internal_tangent) {
  // Circles with different radii touching internally
  auto result = cavccpp::intrCircle2Circle2(3.0, Vector2{0.0, 0.0}, 1.0, Vector2{2.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::OneIntersect);
  EXPECT_TRUE(approxEqual(result.point1, Vector2{3.0, 0.0}));
}

TEST(intrcircle2circle2, very_small_circles) {
  // Test with very small circles
  auto result = cavccpp::intrCircle2Circle2(0.001, Vector2{0.0, 0.0}, 0.001, Vector2{0.0015, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::TwoIntersects);
}

TEST(intrcircle2circle2, large_circles) {
  // Test with large circles
  auto result =
      cavccpp::intrCircle2Circle2(1000.0, Vector2{0.0, 0.0}, 1000.0, Vector2{1500.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::TwoIntersects);
}

TEST(intrcircle2circle2, zero_radius_first) {
  // Point (zero radius) and circle
  auto result = cavccpp::intrCircle2Circle2(0.0, Vector2{1.0, 0.0}, 1.0, Vector2{0.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::OneIntersect);
  EXPECT_TRUE(approxEqual(result.point1, Vector2{1.0, 0.0}));
}

TEST(intrcircle2circle2, zero_radius_second) {
  // Circle and point (zero radius)
  auto result = cavccpp::intrCircle2Circle2(1.0, Vector2{0.0, 0.0}, 0.0, Vector2{1.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::OneIntersect);
  EXPECT_TRUE(approxEqual(result.point1, Vector2{1.0, 0.0}));
}

TEST(intrcircle2circle2, both_zero_radius_same_point) {
  // Two points at same location
  auto result = cavccpp::intrCircle2Circle2(0.0, Vector2{1.0, 1.0}, 0.0, Vector2{1.0, 1.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::Coincident);
}

TEST(intrcircle2circle2, both_zero_radius_different_points) {
  // Two points at different locations
  auto result = cavccpp::intrCircle2Circle2(0.0, Vector2{0.0, 0.0}, 0.0, Vector2{1.0, 1.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::NoIntersect);
}

TEST(intrcircle2circle2, negative_coordinates) {
  // Test with negative coordinates
  auto result = cavccpp::intrCircle2Circle2(1.0, Vector2{-1.0, -1.0}, 1.0, Vector2{-1.0, 0.0});
  EXPECT_EQ(result.intrType, Circle2Circle2IntrType::TwoIntersects);
}

TEST(intrcircle2circle2, precision_edge_case) {
  // Test near the precision boundary
  double distance = 2.0 + cavccpp::utils::realThreshold<double>() / 2.0;
  auto result = cavccpp::intrCircle2Circle2(1.0, Vector2{0.0, 0.0}, 1.0, Vector2{distance, 0.0});
  // Should be close to tangent, might be OneIntersect or TwoIntersects depending on precision
  EXPECT_TRUE(result.intrType == Circle2Circle2IntrType::OneIntersect ||
              result.intrType == Circle2Circle2IntrType::TwoIntersects);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
