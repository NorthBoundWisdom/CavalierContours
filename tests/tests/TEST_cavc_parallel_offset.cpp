#include <gmock/gmock.h>
#include <vector>

#include <gtest/gtest.h>

#include "c_api_include/cavaliercontours.h"
#include "c_api_test_helpers.hpp"

namespace t = testing;

struct ParallelOffsetTestCase {
  std::string name;
  cavc_real offsetDelta;
  cavc_pline *pline;
  std::vector<PolylineProperties> expectedResult;

  ParallelOffsetTestCase(std::string p_name, cavc_real p_delta,
                         std::vector<cavc_vertex> const &p_plineVertexes, bool isClosed,
                         std::vector<PolylineProperties> p_expectedResult)
      : name(std::move(p_name)), offsetDelta(p_delta),
        pline(plineFromVertexes(p_plineVertexes, isClosed)),
        expectedResult(std::move(p_expectedResult)) {}
};
namespace {

[[maybe_unused]] std::ostream &operator<<(std::ostream &os, ParallelOffsetTestCase const &c) {
  os << "{ " << c.name << ", offsetDelta: " << c.offsetDelta << " }";
  return os;
}

std::vector<ParallelOffsetTestCase> createSpecificCases() {
  std::vector<ParallelOffsetTestCase> cases;
  {
    // offset arc just past line, in this case float epsilon thresholding can cause failures
    std::vector<cavc_vertex> vertexes = {{27.804688, 1, 0},
                                         {28.46842055794889, 0.3429054695163245, 0},
                                         {32.34577133994935, 0.9269762697003898, 0},
                                         {32.38116957207762, 1.451312562563487, 0},
                                         {31.5, 1, -0.31783751349740424},
                                         {30.79289310940682, 1.5, 0},
                                         {29.20710689059337, 1.5, -0.31783754777018053},
                                         {28.49999981323106, 1.00000000000007, 0.0}};
    std::vector<PolylineProperties> expectedResult;
    expectedResult.emplace_back(4, 0.094833810726263, 1.8213211761499, 31.533345690439,
                                0.90572346564886, 32.26949555256, 1.2817628453883);
    expectedResult.emplace_back(6, 1.7197931450343, 7.5140262005179, 28.047835685678,
                                0.44926177903859, 31.495431966272, 1.4);
    cases.emplace_back("offset_arc_just_past_line1", 0.1, std::move(vertexes), true,
                       std::move(expectedResult));
  }

  {
    // first vertex position is ontop of intersect with second segment
    std::vector<cavc_vertex> vertexes = {{27.804688, 1, 0},
                                         {27.804688, 0.75, 0},
                                         {32.195313, 0.75, 0},
                                         {32.195313, 1, 0},
                                         {31.5, 1, -0.3178375134974},
                                         {30.792893109407, 1.5, 0},
                                         {29.207106890593, 1.5, -0.31783754777018},
                                         {28.499999813231, 1.0000000000001, 0}};
    std::vector<PolylineProperties> expectedResult;
    expectedResult.emplace_back(4, 0.36247092523069, 3.593999211522, 29.16143806012, 1,
                                30.838561906052, 1.25);
    cases.emplace_back("intersect_ontop_first_vertex", 0.25, std::move(vertexes), true,
                       std::move(expectedResult));
  }

  {
    // collapsed rectangle
    std::vector<cavc_vertex> vertexes = {{0, 0, 0}, {120, 0, 0}, {120, 40, 0}, {0, 40, 0}};
    // expect no results
    std::vector<PolylineProperties> expectedResult;
    cases.emplace_back("collapsed_rectangle", 30.0, std::move(vertexes), true,
                       std::move(expectedResult));
  }

  return cases;
}
static std::vector<ParallelOffsetTestCase> specificCases = createSpecificCases();

std::vector<ParallelOffsetTestCase> createSimpleCases() {
  std::vector<ParallelOffsetTestCase> cases;

  // Rectangles
  {
    // closed rectangle offset inward
    std::vector<cavc_vertex> vertexes = {{0, 0, 0}, {20, 0, 0}, {20, 10, 0}, {0, 10, 0}};
    std::vector<PolylineProperties> expectedResult;
    expectedResult.emplace_back(4, 96, 44, 2, 2, 18, 8);
    cases.emplace_back("closed_rectangle_inward", 2, std::move(vertexes), true,
                       std::move(expectedResult));
  }

  {
    // open rectangle offset inward
    std::vector<cavc_vertex> vertexes = {{0, 0, 0}, {20, 0, 0}, {20, 10, 0}, {0, 10, 0}, {0, 0, 0}};
    std::vector<PolylineProperties> expectedResult;
    expectedResult.emplace_back(5, 0, 44, 2, 2, 18, 8);
    cases.emplace_back("open_rectangle_inward", 2, std::move(vertexes), false,
                       std::move(expectedResult));
  }

  {
    // closed rectangle offset outward
    std::vector<cavc_vertex> vertexes = {{0, 0, 0}, {20, 0, 0}, {20, 10, 0}, {0, 10, 0}};
    std::vector<PolylineProperties> expectedResult;
    expectedResult.emplace_back(8, 332.56637061436, 72.566370614359, -2, -2, 22, 12);
    cases.emplace_back("closed_rectangle_outward", -2, std::move(vertexes), true,
                       std::move(expectedResult));
  }

  {
    // open rectangle offset outward
    std::vector<cavc_vertex> vertexes = {{0, 0, 0}, {20, 0, 0}, {20, 10, 0}, {0, 10, 0}, {0, 0, 0}};
    std::vector<PolylineProperties> expectedResult;
    expectedResult.emplace_back(8, 0, 69.424777960769, -2, -2, 22, 12);
    cases.emplace_back("open_rectangle_outward", -2, std::move(vertexes), false,
                       std::move(expectedResult));
  }

  {
    // closed rectangle offset inward into coincident line
    std::vector<cavc_vertex> vertexes = {{0, 0, 0}, {20, 0, 0}, {20, 10, 0}, {0, 10, 0}};
    std::vector<PolylineProperties> expectedResult;
    expectedResult.emplace_back(2, 0, 20, 5, 5, 15, 5);
    cases.emplace_back("closed_rectangle_coincident", 5, std::move(vertexes), true,
                       std::move(expectedResult));
  }

  // Diamonds
  {
    // closed diamond offset inward
    std::vector<cavc_vertex> vertexes = {{-10, 0, 0}, {0, 10, 0}, {10, 0, 0}, {0, -10, 0}};
    std::vector<PolylineProperties> expectedResult;
    expectedResult.emplace_back(4, -17.157287525381, 16.568542494924, -2.9289321881345,
                                -2.9289321881345, 2.9289321881345, 2.9289321881345);
    cases.emplace_back("closed_diamond_inward", -5, std::move(vertexes), true,
                       std::move(expectedResult));
  }

  {
    // open diamond offset inward
    std::vector<cavc_vertex> vertexes = {
        {-10, 0, 0}, {0, 10, 0}, {10, 0, 0}, {0, -10, 0}, {-10, 0, 0}};
    std::vector<PolylineProperties> expectedResult;
    expectedResult.emplace_back(5, 0, 16.568542494924, -2.9289321881345, -2.9289321881345,
                                2.9289321881345, 2.9289321881345);
    cases.emplace_back("open_diamond_inward", -5, std::move(vertexes), false,
                       std::move(expectedResult));
  }

  {
    // closed diamond offset outward
    std::vector<cavc_vertex> vertexes = {{-10, 0, 0}, {0, 10, 0}, {10, 0, 0}, {0, -10, 0}};
    std::vector<PolylineProperties> expectedResult;
    expectedResult.emplace_back(8, -561.38252881436, 87.984469030822, -15, -15, 15, 15);
    cases.emplace_back("closed_diamond_outward", 5, std::move(vertexes), true,
                       std::move(expectedResult));
  }

  {
    // open diamond offset outward
    std::vector<cavc_vertex> vertexes = {
        {-10, 0, 0}, {0, 10, 0}, {10, 0, 0}, {0, -10, 0}, {-10, 0, 0}};
    std::vector<PolylineProperties> expectedResult;
    expectedResult.emplace_back(8, 0, 80.130487396847, -13.535533905933, -15, 15, 15);
    cases.emplace_back("open_diamond_outward", 5, std::move(vertexes), false,
                       std::move(expectedResult));
  }

  return cases;
}
} // namespace

static std::vector<ParallelOffsetTestCase> simpleCases = createSimpleCases();

class cavc_parallel_offsetTests : public t::TestWithParam<ParallelOffsetTestCase> {};

INSTANTIATE_TEST_SUITE_P(specifc_cases, cavc_parallel_offsetTests, t::ValuesIn(specificCases));

INSTANTIATE_TEST_SUITE_P(simple_cases, cavc_parallel_offsetTests, t::ValuesIn(simpleCases));

// Performs the offset and tests results of offset properties against expected properties
TEST_P(cavc_parallel_offsetTests, parallel_offset_test) {
  ParallelOffsetTestCase const &testCase = GetParam();
  cavc_pline_list *results;
  cavc_parallel_offset(testCase.pline, testCase.offsetDelta, &results,
                       defaultParallelOffsetOptions());
  ASSERT_EQ(cavc_pline_list_count(results), testCase.expectedResult.size());

  std::vector<PolylineProperties> resultsProperties;
  resultsProperties.reserve(testCase.expectedResult.size());

  for (uint32_t i = 0; i < testCase.expectedResult.size(); ++i) {
    cavc_pline *pline = cavc_pline_list_get(results, i);
    resultsProperties.emplace_back(pline);
  }

  ASSERT_THAT(resultsProperties, t::UnorderedPointwise(t::Eq(), testCase.expectedResult));

  cavc_pline_list_delete(results);
}

// Geometric properties should hold when polyline is reversed and offset
TEST_P(cavc_parallel_offsetTests, reversed_parallel_offset_test) {
  ParallelOffsetTestCase const &testCase = GetParam();
  // create the reversed polyline
  cavc_pline *revPline = createRevseredPline(testCase.pline);
  // negate delta to offset the same direction
  cavc_real delta = -testCase.offsetDelta;
  // sign of the area is also negated
  std::vector<PolylineProperties> expectedResult = testCase.expectedResult;
  for (auto &pp : expectedResult) {
    pp.area = -pp.area;
  }

  cavc_pline_list *results = nullptr;
  cavc_parallel_offset(revPline, delta, &results, defaultParallelOffsetOptions());
  ASSERT_EQ(cavc_pline_list_count(results), expectedResult.size());

  std::vector<PolylineProperties> resultsProperties;
  resultsProperties.reserve(expectedResult.size());

  for (uint32_t i = 0; i < expectedResult.size(); ++i) {
    cavc_pline *pline = cavc_pline_list_get(results, i);
    resultsProperties.emplace_back(pline);
  }

  ASSERT_THAT(resultsProperties, t::UnorderedPointwise(t::Eq(), expectedResult));

  cavc_pline_delete(revPline);
  cavc_pline_list_delete(results);
}

TEST(cavc_parallel_offsetTests, parallel_offset_does_not_modify_input_test) {
  ParallelOffsetTestCase const &testCase = simpleCases[0];
  std::vector<cavc_vertex> vertexesBefore(cavc_pline_vertex_count(testCase.pline));
  cavc_pline_vertex_data(testCase.pline, &vertexesBefore[0]);

  cavc_pline_list *results = nullptr;
  cavc_parallel_offset(testCase.pline, testCase.offsetDelta, &results,
                       defaultParallelOffsetOptions());

  std::vector<cavc_vertex> vertexesAfter(cavc_pline_vertex_count(testCase.pline));
  cavc_pline_vertex_data(testCase.pline, &vertexesAfter[0]);

  ASSERT_THAT(vertexesAfter, t::Pointwise(VertexEqual(), vertexesBefore));

  cavc_pline_list_delete(results);
}

TEST(cavc_parallel_offsetTests, parallel_offset_handles_repeat_position_input_test) {
  std::vector<cavc_vertex> vertexes = {{0, 0, 0}, {20, 0, 0}, {20, 0, 0}, {20, 10, 0}, {0, 10, 0}};
  cavc_pline *pline = plineFromVertexes(vertexes, true);

  cavc_pline_list *results = nullptr;
  cavc_parallel_offset(pline, -2, &results, defaultParallelOffsetOptions());

  ASSERT_EQ(cavc_pline_list_count(results), 1u);
  cavc_pline *resultPline = cavc_pline_list_get(results, 0);
  PolylineProperties actual(resultPline);
  PolylineProperties expected(8, 332.56637061436, 72.566370614359, -2, -2, 22, 12);
  ASSERT_EQ(actual, expected);

  cavc_pline_list_delete(results);
  cavc_pline_delete(pline);
}

TEST(cavc_parallel_offsetTests, parallel_offset_handles_collapsed_loop_input_test) {
  std::vector<cavc_vertex> vertexes = {{1, 0, 0}, {-1, 0, 0}};
  cavc_pline *pline = plineFromVertexes(vertexes, true);

  {
    cavc_pline_list *results = nullptr;
    cavc_parallel_offset(pline, 1.0, &results, defaultParallelOffsetOptions());

    ASSERT_EQ(cavc_pline_list_count(results), 1u);
    cavc_pline *resultPline = cavc_pline_list_get(results, 0);
    PolylineProperties actual(resultPline);
    PolylineProperties expected(4, -7.141592653589793, 10.283185307179586, -2, -1, 2, 1);
    ASSERT_EQ(actual, expected);

    cavc_pline_list_delete(results);
  }

  {
    cavc_pline_list *results = nullptr;
    cavc_parallel_offset(pline, -1.0, &results, defaultParallelOffsetOptions());

    ASSERT_EQ(cavc_pline_list_count(results), 1u);
    cavc_pline *resultPline = cavc_pline_list_get(results, 0);
    PolylineProperties actual(resultPline);
    PolylineProperties expected(4, 7.141592653589793, 10.283185307179586, -2, -1, 2, 1);
    ASSERT_EQ(actual, expected);

    cavc_pline_list_delete(results);
  }

  cavc_pline_delete(pline);
}

TEST(cavc_parallel_offsetTests, parallel_offset_self_intersect_hint_handles_repeated_offset_case) {
  std::vector<cavc_vertex> vertexes = {
      {2.0, 11.0, -0.6681786379192991},
      {2.7071067811865475, 9.292893218813452, 0.0},
      {-0.2928932188134524, 6.292893218813452, -0.6681786379192989},
      {-2.0, 7.0, 0.0},
      {-2.0, 15.0, -0.6681786379192989},
      {-0.2928932188134524, 15.707106781186548, 0.0},
      {2.7071067811865475, 12.707106781186548, -0.6681786379192991},
  };
  cavc_pline *pline = plineFromVertexes(vertexes, true);

  cavc_parallel_offset_options options = defaultParallelOffsetOptions();
  options.may_have_self_intersects = 1;

  cavc_pline_list *results = nullptr;
  cavc_parallel_offset(pline, 1.0, &results, options);

  ASSERT_EQ(cavc_pline_list_count(results), 1u);
  cavc_pline *resultPline = cavc_pline_list_get(results, 0);
  PolylineProperties actual(resultPline);
  PolylineProperties expected(7, -64.3633792727984, 31.14604709099094, -3.0, 5.0, 4.0, 17.0);
  ASSERT_EQ(actual, expected);

  cavc_pline_list_delete(results);
  cavc_pline_delete(pline);
}

TEST(cavc_parallel_offsetTests, parallel_offset_handles_overlap_result_regression) {
  std::vector<cavc_vertex> vertexes = {
      {0.0, 0.0, 0.0},
      {432.22004474869937, 0.0, 0.0},
      {432.22004474869937, -620.7191231042452, 0.0},
      {414.22004474869937, -620.7191231042452, 0.0},
      {414.22004474869937, -18.0, 0.0},
      {0.0, -18.0, 0.0},
  };
  cavc_pline *pline = plineFromVertexes(vertexes, true);

  cavc_pline_list *results = nullptr;
  cavc_parallel_offset(pline, -9.0, &results, defaultParallelOffsetOptions());

  ASSERT_EQ(cavc_pline_list_count(results), 1u);
  cavc_pline *resultPline = cavc_pline_list_get(results, 0);
  EXPECT_TRUE(cavc_pline_is_closed(resultPline));
  EXPECT_GE(cavc_pline_vertex_count(resultPline), 5u);
  EXPECT_NEAR(cavc_get_area(resultPline), -17.38274876480773, 1e-5);
  EXPECT_NEAR(cavc_get_path_length(resultPline), 2030.0155026470434, 1e-5);

  cavc_real min_x = 0.0;
  cavc_real min_y = 0.0;
  cavc_real max_x = 0.0;
  cavc_real max_y = 0.0;
  cavc_get_extents(resultPline, &min_x, &min_y, &max_x, &max_y);
  EXPECT_NEAR(min_x, 9.0, 1e-5);
  EXPECT_NEAR(min_y, -611.7191231042452, 1e-5);
  EXPECT_NEAR(max_x, 423.22004474869937, 1e-5);
  EXPECT_NEAR(max_y, -9.0, 1e-5);

  cavc_pline_list_delete(results);
  cavc_pline_delete(pline);
}

int main(int argc, char **argv) {
  t::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
