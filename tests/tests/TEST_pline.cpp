#include <gmock/gmock.h>
#include <vector>

#include <gtest/gtest.h>

#include <cavc/polyline.hpp>

#include "casebuilder.hpp"

TEST(basic, basic_extent) {
  //   double minX;
  //   double minY;
  //   double maxX;
  //   double maxY;
  //   // Quad Circle
  //   std::vector<PlineVertex> plineAVertexes = CaseBuilder::quarterArcCase();
  //   auto plineA = cavc::Polyline::fromVertexes(plineAVertexes, true);
  //   plineA.getExtents(&minX, &minY, &maxX, &maxY);
  //   EXPECT_NEAR(minX, 0, 1e-15);
  //   EXPECT_NEAR(minY, -1, 1e-15);
  //   EXPECT_NEAR(maxX, 1, 1e-15);
  //   EXPECT_NEAR(maxY, 0, 1e-15);

  //   // Half circle cw
  //   std::vector<cavc_vertex> plineBVertexes = {{1.0, 0.00, -1.00}, {0.0, 0.00, 0.00}};
  //   cavc_pline *plineB1 = PolylineFactory::plineFromVertexes(plineBVertexes, true);
  //   cavc_get_extents(plineB1, &minX, &minY, &maxX, &maxY);
  //   EXPECT_NEAR(minX, 0, 1e-15);
  //   EXPECT_NEAR(minY, -0.5, 1e-15);
  //   EXPECT_NEAR(maxX, 1, 1e-15);
  //   EXPECT_NEAR(maxY, 0, 1e-15);

  //   // Half circle ccw
  //   plineBVertexes.front().bulge = 1;
  //   cavc_pline *plineB2 = PolylineFactory::plineFromVertexes(plineBVertexes, true);
  //   cavc_get_extents(plineB2, &minX, &minY, &maxX, &maxY);
  //   EXPECT_NEAR(minX, 0, 1e-15);
  //   EXPECT_NEAR(minY, 0, 1e-15);
  //   EXPECT_NEAR(maxX, 1, 1e-15);
  //   EXPECT_NEAR(maxY, 0.5, 1e-15);

  //   // Half circle ccw
  //   std::vector<cavc_vertex> plineCVertexes = {
  //       {0.0, 0.00, 1.00},
  //       {0.0, 1.00, 0.00},
  //   };
  //   cavc_pline *plineC1 = PolylineFactory::plineFromVertexes(plineCVertexes, true);
  //   cavc_get_extents(plineC1, &minX, &minY, &maxX, &maxY);
  //   EXPECT_NEAR(minX, 0, 1e-15);
  //   EXPECT_NEAR(minY, 0, 1e-15);
  //   EXPECT_NEAR(maxX, 0.5, 1e-15);
  //   EXPECT_NEAR(maxY, 1, 1e-15);

  //   // Half circle cw
  //   plineCVertexes.front().bulge = -1;
  //   cavc_pline *plineC2 = PolylineFactory::plineFromVertexes(plineCVertexes, true);
  //   cavc_get_extents(plineC2, &minX, &minY, &maxX, &maxY);
  //   EXPECT_NEAR(minX, -0.5, 1e-15);
  //   EXPECT_NEAR(minY, 0, 1e-15);
  //   EXPECT_NEAR(maxX, 0, 1e-15);
  //   EXPECT_NEAR(maxY, 1, 1e-15);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
