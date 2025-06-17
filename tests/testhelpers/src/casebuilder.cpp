#include "casebuilder.hpp"

#include <algorithm>

std::vector<PlineVertex> CaseBuilder::reverseDirection(std::vector<PlineVertex> const &vertices) {
  if (vertices.size() < 2) {
    return vertices;
  }
  std::vector<PlineVertex> result = vertices;
  std::reverse(result.begin(), result.end());
  double firstBulge = result[0].bulge();
  for (std::size_t i = 1; i < result.size(); ++i) {
    result[i - 1].bulge() = -result[i].bulge();
  }

  result.back().bulge() = -firstBulge;
  return result;
}

std::vector<PlineVertex> CaseBuilder::simpleRectangle() {
  std::vector<PlineVertex> result;
  result.reserve(4);

  result.push_back({0, 0, 0});
  result.push_back({1, 0, 0});
  result.push_back({1, 1, 0});
  result.push_back({0, 1, 0});

  return result;
}

std::vector<PlineVertex> CaseBuilder::offsetCase() {
  std::vector<PlineVertex> res;
  res.push_back({0, 25, 1});
  res.push_back({0, 0, 0});
  res.push_back({2, 0, 1});
  res.push_back({10, 0, -0.5});
  res.push_back({8, 9, 0.374794619217547});
  res.push_back({21, 0, 0});
  res.push_back({23, 0, 1});
  res.push_back({32, 0, -0.5});
  res.push_back({28, 0, 0.5});
  res.push_back({39, 21, 0});
  res.push_back({28, 12, 0.5});
  return res;
}

std::vector<PlineVertex> CaseBuilder::twoBadArcCase() {
  std::vector<PlineVertex> res;
  res.push_back({0, 0.5, -0.618033988749895});
  res.push_back({4, 0.5, -0.5});
  return res;
}

std::vector<PlineVertex> CaseBuilder::figureEightCase() {
  std::vector<PlineVertex> res;
  res.push_back({0, 0, 1});
  res.push_back({2, 0, 1});
  res.push_back({0, 0, -1});
  res.push_back({-2, 0, -1});
  return res;
}

std::vector<PlineVertex> CaseBuilder::complexSelfInterectCase() {
  std::vector<PlineVertex> res;
  res.push_back({0, 25, 1.995510113000107});
  res.push_back({6.684425313612419, 16.05394554133331, 0});
  res.push_back({-2.468698351579587, 15.973654982866707, 8.639429900629425});
  res.push_back({0.180890077818098, 14.60871548893456, 0.917148791358611});
  res.push_back({8, 9, 0.374794619217547});
  res.push_back({21, 0, 0});
  res.push_back({20.735273045266833, 4.331524005210202, 1.455901596158439});
  res.push_back({32, 0, -0.877081908564703});
  res.push_back({15.034643394138465, 9.55041030553898, 3.717958973262705});
  res.push_back({12.947088874006965, 14.126972138134981, 0});
  res.push_back({7.246459222878594, 20.9516696077957, 0.547663466075712});
  return res;
}

std::vector<PlineVertex> CaseBuilder::closedLineArcCase() {
  std::vector<PlineVertex> res;
  res.push_back({5, 5, 0});
  res.push_back({3, 9, 0});
  res.push_back({0, 10, 0});
  res.push_back({-4, 8, 0});
  res.push_back({-5, 5, 1});
  return res;
}

std::vector<PlineVertex> CaseBuilder::quarterArcCase() {
  std::vector<PlineVertex> res;
  res.push_back({1, 0, -0.414213562373095});
  res.push_back({0, -1, 0.00});
  return res;
}

std::vector<PlineVertex> CaseBuilder::positiveCircle() {
  std::vector<PlineVertex> data1;
  data1.push_back({0, 0, 1});
  data1.push_back({10, 0, 1});
  return data1;
}

std::vector<PlineVertex> CaseBuilder::negativeCircle() {
  std::vector<PlineVertex> data1;
  data1.push_back({0, 0, -1});
  data1.push_back({10, 0, -1});
  return data1;
}

std::vector<std::vector<PlineVertex>> CaseBuilder::simpleBoolCase() {
  std::vector<PlineVertex> data1;
  data1.push_back({0, 1, 1});
  data1.push_back({10, 1, 1});

  std::vector<PlineVertex> data2;
  data2.push_back({3, -10, 0});
  data2.push_back({6, -10, 0});
  data2.push_back({6, 10, 0});
  data2.push_back({3, 10, 0});

  return {data1, data2};
}
