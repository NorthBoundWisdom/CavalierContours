#ifndef CAVC_TESTHELPERS_HPP
#define CAVC_TESTHELPERS_HPP
#include <cmath>
#include <gmock/gmock.h>
#include <iostream>

#include <gtest/gtest.h>

#include "c_api_include/cavaliercontours.h"

constexpr inline cavc_real PI() { return 3.14159265358979323846264338327950288; }
constexpr inline cavc_real TEST_EPSILON() { return 1e-5; }
template <typename Real> inline bool fuzzyEqual(Real const &left, cavc_real const &right) {
  return std::abs(left - right) < TEST_EPSILON();
}

// type to hold summary properties of a polyline for test comparing, acts as a sort of geometric
// hash of the polyline, it is very unlikely that two polylines have the same PolylineProperties
// without being the same polyline, especially accidentally via generation in an algorithm
struct PolylineProperties {
  std::size_t vertexCount;
  cavc_real area;
  cavc_real pathLength;
  cavc_real minX;
  cavc_real minY;
  cavc_real maxX;
  cavc_real maxY;

  PolylineProperties(std::size_t p_vertexCount, cavc_real p_area, cavc_real p_pathLength,
                     cavc_real p_minX, cavc_real p_minY, cavc_real p_maxX, cavc_real p_maxY)
      : vertexCount(p_vertexCount), area(p_area), pathLength(p_pathLength), minX(p_minX),
        minY(p_minY), maxX(p_maxX), maxY(p_maxY) {}

  PolylineProperties(cavc_pline *pline) {
    vertexCount = cavc_pline_vertex_count(pline);
    area = cavc_get_area(pline);
    pathLength = cavc_get_path_length(pline);
    cavc_get_extents(pline, &minX, &minY, &maxX, &maxY);
  }
};

MATCHER(EqIgnoreSignOfArea, "") {
  auto const &left = std::get<0>(arg);
  auto const &right = std::get<1>(arg);
  return left.vertexCount == right.vertexCount &&
         fuzzyEqual(std::abs(left.area), std::abs(right.area)) &&
         fuzzyEqual(left.pathLength, right.pathLength) && fuzzyEqual(left.minX, right.minX) &&
         fuzzyEqual(left.minY, right.minY) && fuzzyEqual(left.maxX, right.maxX) &&
         fuzzyEqual(left.maxY, right.maxY);
}

// fuzzy equality operator== for testing
inline bool operator==(PolylineProperties const &left, PolylineProperties const &right) {
  return left.vertexCount == right.vertexCount && fuzzyEqual(left.area, right.area) &&
         fuzzyEqual(left.pathLength, right.pathLength) && fuzzyEqual(left.minX, right.minX) &&
         fuzzyEqual(left.minY, right.minY) && fuzzyEqual(left.maxX, right.maxX) &&
         fuzzyEqual(left.maxY, right.maxY);
}

inline std::ostream &operator<<(std::ostream &os, PolylineProperties const &p) {
  os << "{ vertexCount: " << p.vertexCount << ", area: " << p.area
     << ", pathLength: " << p.pathLength << ", minX: " << p.minX << ", minY: " << p.minY
     << ", maxX: " << p.maxX << ", maxY: " << p.maxY << " }";
  return os;
}

// Custom function to print differences
inline void PrintDiff(const PolylineProperties &expected, const PolylineProperties &actual) {
  if (expected.vertexCount != actual.vertexCount) {
    std::cout << "vertexCount: expected " << expected.vertexCount << ", actual "
              << actual.vertexCount << "\n";
  }
  if (!fuzzyEqual(expected.area, actual.area)) {
    std::cout << "area: expected " << expected.area << ", actual " << actual.area << "\n";
  }
  if (!fuzzyEqual(expected.pathLength, actual.pathLength)) {
    std::cout << "pathLength: expected " << expected.pathLength << ", actual " << actual.pathLength
              << "\n";
  }
  if (!fuzzyEqual(expected.minX, actual.minX)) {
    std::cout << "minX: expected " << expected.minX << ", actual " << actual.minX << "\n";
  }
  if (!fuzzyEqual(expected.minY, actual.minY)) {
    std::cout << "minY: expected " << expected.minY << ", actual " << actual.minY << "\n";
  }
  if (!fuzzyEqual(expected.maxX, actual.maxX)) {
    std::cout << "maxX: expected " << expected.maxX << ", actual " << actual.maxX << "\n";
  }
  if (!fuzzyEqual(expected.maxY, actual.maxY)) {
    std::cout << "maxY: expected " << expected.maxY << ", actual " << actual.maxY << "\n";
  }
}

inline bool vertexesFuzzyEqual(cavc_vertex const &left, cavc_vertex const &right) {
  return fuzzyEqual(left.x, right.x) && fuzzyEqual(left.y, right.y) &&
         fuzzyEqual(left.bulge, right.bulge);
}

template <typename Container>
inline std::size_t nextWrappingIndex(Container const &container, std::size_t index) {
  if (index == container.size() - 1) {
    return 0;
  }
  return index + 1;
}

MATCHER(VertexFuzzyEqual, "") {
  auto const &left = std::get<0>(arg);
  auto const &right = std::get<1>(arg);
  return vertexesFuzzyEqual(left, right);
}

MATCHER(VertexEqual, "") {
  auto const &left = std::get<0>(arg);
  auto const &right = std::get<1>(arg);
  return left.x == right.x && left.y == right.y && left.bulge == right.bulge;
}

MATCHER(PointFuzzyEqual, "") {
  auto const &left = std::get<0>(arg);
  auto const &right = std::get<1>(arg);
  return fuzzyEqual(left.x, right.x) && fuzzyEqual(left.y, right.y);
}

MATCHER_P(VertexListsFuzzyEqual, isClosed, "") {
  auto const &left = std::get<0>(arg);
  auto const &right = std::get<1>(arg);
  if (left.size() != right.size()) {
    *result_listener << "sizes of vertex lists do not match ";
    return false;
  }

  std::size_t const vertexCount = left.size();

  if (!isClosed) {
    // open polyline indexes much match up
    for (std::size_t i = 0; i < vertexCount; ++i) {
      if (!vertexesFuzzyEqual(left[i], right[i])) {
        *result_listener << "vertexes not equal at index: " << i << " ";
        return false;
      }
    }
    return true;
  }

  // vertexes may not have same indexes in the case of a closed polyline, find first matching and
  // start matching from there
  std::size_t startIndex = 0;
  for (; startIndex < vertexCount; ++startIndex) {
    if (vertexesFuzzyEqual(left[0], right[startIndex])) {
      break;
    }
  }

  if (startIndex == vertexCount) {
    *result_listener << "did not find matching vertex to start with ";
    return false;
  }

  *result_listener << " started matching at index: " << startIndex << " ";
  // first one already matched, start at next index
  std::size_t index = nextWrappingIndex(left, startIndex);
  for (std::size_t i = 1; i < vertexCount; ++i) {
    if (!vertexesFuzzyEqual(left[i], right[index])) {
      *result_listener << "vertexes not equal at index: " << index << " ";
      return false;
    }
    index = nextWrappingIndex(left, index);
  }

  return true;
}

inline std::ostream &operator<<(std::ostream &os, cavc_vertex const &v) {
  os << '[' << v.x << "," << v.y << "," << v.bulge << ']';
  return os;
}

inline std::ostream &operator<<(std::ostream &os, cavc_point const &p) {
  os << '[' << p.x << "," << p.y << ']';
  return os;
}

#endif // CAVC_TESTHELPERS_HPP
