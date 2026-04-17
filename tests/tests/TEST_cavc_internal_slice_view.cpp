#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "cavc/internal/plinesliceview.hpp"
#include "c_api_test_helpers.hpp"

namespace t = testing;

namespace {
using cavc::PlineVertex;
using cavc::Polyline;
using cavc::Vector2;
using cavc::internal::PlineSliceViewData;
using cavc::internal::PlineSliceViewValidation;

std::vector<cavc_vertex> toVertexList(Polyline<double> const &pline) {
  std::vector<cavc_vertex> result;
  result.reserve(pline.size());
  for (std::size_t i = 0; i < pline.size(); ++i) {
    result.push_back({pline[i].x(), pline[i].y(), pline[i].bulge()});
  }
  return result;
}
} // namespace

TEST(InternalSliceView, CreateOnSingleSegmentLineProducesExpectedPolyline) {
  Polyline<double> source;
  source.addVertex(0.0, 0.0, 0.0);
  source.addVertex(10.0, 0.0, 0.0);

  auto view = PlineSliceViewData<double>::createOnSingleSegment(
      source, 0, PlineVertex<double>(2.0, 0.0, 0.0), Vector2<double>{7.0, 0.0});
  ASSERT_TRUE(view.has_value());
  EXPECT_EQ(view->debugValidateForSource(source), PlineSliceViewValidation::Valid);

  Polyline<double> materialized = view->toPolyline(source);
  std::vector<cavc_vertex> expected = {{2.0, 0.0, 0.0}, {7.0, 0.0, 0.0}};
  EXPECT_THAT(toVertexList(materialized), t::Pointwise(VertexFuzzyEqual(), expected));

  int visitedSegments = 0;
  view->visitSegments(source, [&](PlineVertex<double> const &v1, PlineVertex<double> const &v2) {
    ++visitedSegments;
    EXPECT_TRUE(vertexesFuzzyEqual({v1.x(), v1.y(), v1.bulge()}, expected[0]));
    EXPECT_TRUE(vertexesFuzzyEqual({v2.x(), v2.y(), v2.bulge()}, expected[1]));
    return true;
  });
  EXPECT_EQ(visitedSegments, 1);
}

TEST(InternalSliceView, CreateFromSlicePointsArcOnSingleSegmentTrimsBulgesCorrectly) {
  Polyline<double> source;
  source.addVertex(0.0, 0.0, 1.0);
  source.addVertex(10.0, 0.0, 0.0);

  Vector2<double> const center{5.0, 0.0};
  Vector2<double> const startPoint =
      cavc::pointOnCircle(5.0, center, 1.25 * cavc::utils::pi<double>());
  Vector2<double> const endPoint =
      cavc::pointOnCircle(5.0, center, 1.75 * cavc::utils::pi<double>());

  auto view =
      PlineSliceViewData<double>::createFromSlicePoints(source, startPoint, 0, endPoint, 0);
  ASSERT_TRUE(view.has_value());
  EXPECT_EQ(view->debugValidateForSource(source), PlineSliceViewValidation::Valid);

  auto startSplit = cavc::splitAtPoint(source[0], source[1], startPoint);
  auto expectedStart = cavc::splitAtPoint(startSplit.splitVertex, source[1], endPoint).updatedStart;
  EXPECT_NEAR(view->updatedStart.x(), expectedStart.x(), TEST_EPSILON());
  EXPECT_NEAR(view->updatedStart.y(), expectedStart.y(), TEST_EPSILON());
  EXPECT_NEAR(view->updatedStart.bulge(), expectedStart.bulge(), TEST_EPSILON());
  EXPECT_NEAR(view->updatedEndBulge, expectedStart.bulge(), TEST_EPSILON());

  Polyline<double> materialized = view->toPolyline(source);
  ASSERT_EQ(materialized.size(), 2u);
  EXPECT_NEAR(materialized[0].x(), expectedStart.x(), TEST_EPSILON());
  EXPECT_NEAR(materialized[0].y(), expectedStart.y(), TEST_EPSILON());
  EXPECT_NEAR(materialized[0].bulge(), expectedStart.bulge(), TEST_EPSILON());
  EXPECT_NEAR(materialized[1].x(), endPoint.x(), TEST_EPSILON());
  EXPECT_NEAR(materialized[1].y(), endPoint.y(), TEST_EPSILON());
  EXPECT_NEAR(materialized[1].bulge(), 0.0, TEST_EPSILON());
}

TEST(InternalSliceView, CreateFromSlicePointsWrapsClosedPolyline) {
  Polyline<double> source;
  source.isClosed() = true;
  source.addVertex(0.0, 0.0, 0.0);
  source.addVertex(10.0, 0.0, 0.0);
  source.addVertex(10.0, 10.0, 0.0);
  source.addVertex(0.0, 10.0, 0.0);

  auto view = PlineSliceViewData<double>::createFromSlicePoints(
      source, Vector2<double>{5.0, 10.0}, 2, Vector2<double>{5.0, 0.0}, 0);
  ASSERT_TRUE(view.has_value());
  EXPECT_EQ(view->debugValidateForSource(source), PlineSliceViewValidation::Valid);

  std::vector<cavc_vertex> expected = {
      {5.0, 10.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 0.0}, {5.0, 0.0, 0.0}};
  EXPECT_THAT(toVertexList(view->toPolyline(source)), t::Pointwise(VertexFuzzyEqual(), expected));
}

TEST(InternalSliceView, InvertedDirectionReversesVertexOrder) {
  Polyline<double> source;
  source.isClosed() = true;
  source.addVertex(0.0, 0.0, 0.0);
  source.addVertex(10.0, 0.0, 0.0);
  source.addVertex(10.0, 10.0, 0.0);
  source.addVertex(0.0, 10.0, 0.0);

  auto view = PlineSliceViewData<double>::createFromSlicePoints(
      source, Vector2<double>{5.0, 10.0}, 2, Vector2<double>{5.0, 0.0}, 0);
  ASSERT_TRUE(view.has_value());
  view->invertedDirection = true;

  std::vector<cavc_vertex> expected = {
      {5.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {5.0, 10.0, 0.0}};
  EXPECT_THAT(toVertexList(view->toPolyline(source)), t::Pointwise(VertexFuzzyEqual(), expected));
}

TEST(InternalSliceView, AppendToRemovesRepeatJoinVertex) {
  Polyline<double> source;
  source.addVertex(0.0, 0.0, 0.0);
  source.addVertex(10.0, 0.0, 0.0);
  source.addVertex(10.0, 10.0, 0.0);

  auto firstSlice = PlineSliceViewData<double>::createFromSlicePoints(
      source, Vector2<double>{0.0, 0.0}, 0, Vector2<double>{10.0, 0.0}, 0);
  auto secondSlice = PlineSliceViewData<double>::createFromSlicePoints(
      source, Vector2<double>{10.0, 0.0}, 0, Vector2<double>{10.0, 10.0}, 1);

  ASSERT_TRUE(firstSlice.has_value());
  ASSERT_TRUE(secondSlice.has_value());

  Polyline<double> target;
  firstSlice->appendTo(target, source);
  secondSlice->appendTo(target, source);

  std::vector<cavc_vertex> expected = {
      {0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}, {10.0, 10.0, 0.0}};
  EXPECT_THAT(toVertexList(target), t::Pointwise(VertexFuzzyEqual(), expected));
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
