#include "polylinefactory.hpp"
#include <algorithm>
#include <cmath>

static cavc_point pointOnCircle(cavc_real radius, cavc_point center, cavc_real angle) {
  return {center.x + radius * std::cos(angle), center.y + radius * std::sin(angle)};
}

inline static cavc_real PI() { return 3.14159265358979323846264338327950288; }

std::vector<cavc_vertex> PolylineFactory::createCircle(cavc_real radius, cavc_point center,
                                                       cavc_real vertexRotAngle, bool isCW) {
  std::vector<cavc_vertex> result;
  result.reserve(2);

  cavc_point point1 = pointOnCircle(radius, center, vertexRotAngle);
  cavc_point point2 = pointOnCircle(radius, center, vertexRotAngle + PI());
  cavc_real bulge = isCW ? -1.0 : 1.0;
  result.push_back({point1.x, point1.y, bulge});
  result.push_back({point2.x, point2.y, bulge});

  return result;
}

cavc_pline_ptr PolylineFactory::vertexesToPline(std::vector<cavc_vertex> const &vertexes,
                                                bool isClosed) {
  return cavc_pline_ptr(plineFromVertexes(vertexes, isClosed));
}

cavc_pline *PolylineFactory::plineFromVertexes(std::vector<cavc_vertex> const &vertexes,
                                               bool isClosed) {
  return cavc_pline_new(&vertexes[0], static_cast<uint32_t>(vertexes.size()), isClosed ? 1 : 0);
}

// reverses the direction of the polyline defined by vertexes
void PolylineFactory::reverseDirection(std::vector<cavc_vertex> &vertexes) {
  if (vertexes.size() < 2) {
    return;
  }

  std::reverse(vertexes.begin(), vertexes.end());
  cavc_real firstBulge = vertexes[0].bulge;
  for (std::size_t i = 1; i < vertexes.size(); ++i) {
    vertexes[i - 1].bulge = -vertexes[i].bulge;
  }

  vertexes.back().bulge = -firstBulge;
}

// create a reversed polyline from the given pline (caller must delete the created pline)
cavc_pline *PolylineFactory::createRevseredPline(cavc_pline *pline) {
  uint32_t count = cavc_pline_vertex_count(pline);
  std::vector<cavc_vertex> vertexes(count);
  cavc_pline_vertex_data(pline, &vertexes[0]);
  reverseDirection(vertexes);
  return cavc_pline_new(&vertexes[0], count, cavc_pline_is_closed(pline));
}
