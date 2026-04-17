#ifndef CAVC_INTRLINESEG2CIRCLE2_HPP
#define CAVC_INTRLINESEG2CIRCLE2_HPP
#include "vector2.hpp"

namespace cavc {
template <typename Real> struct IntrLineSeg2Circle2Result {
  // number of interescts found (0, 1, or 2)
  int numIntersects;
  // parametric value for first intersect (if numIntersects > 0) otherwise undefined
  Real t0;
  // parametric value for second intersect (if numintersects > 1) otherwise undefined
  Real t1;
};

// Gets the intersect between a segment and a circle, returning the parametric solution t to the
// segment equation P(t) = v1 + t * (v2 - v1) for t = 0 to t = 1, if t < 0 or t > 1 then intersect
// occurs only when extending the segment out past the points given (if t < 0 intersect nearest v1,
// if t > 0 then intersect nearest v2), intersects are "sticky" and "snap" to tangent points, e.g. a
// segment very close to being a tangent will be returned as a single intersect point
template <typename Real>
IntrLineSeg2Circle2Result<Real> intrLineSeg2Circle2(Vector2<Real> const &p0,
                                                    Vector2<Real> const &p1, Real radius,
                                                    Vector2<Real> const &circleCenter) {
  // This function solves for the line/circle intersects using the line equation rather than
  // directly solving for the parametric variable with the quadratic formula. This is more stable
  // in cases where the line is nearly vertical/horizontal.
  IntrLineSeg2Circle2Result<Real> result;
  Real dx = p1.x() - p0.x();
  Real dy = p1.y() - p0.y();
  Real h = circleCenter.x();
  Real k = circleCenter.y();
  Real const epsilon = utils::realThreshold<Real>();

  auto parametricFromPoint = [&](Vector2<Real> const &point) {
    if (std::abs(dx) < std::abs(dy)) {
      return (point.y() - p0.y()) / dy;
    }
    return (point.x() - p0.x()) / dx;
  };

  if (fuzzyEqual(p0, p1, epsilon)) {
    // v1 = v2, test if point is on the circle
    Real xh = (p0.x() + p1.x()) / Real(2) - h;
    Real yk = (p0.y() + p1.y()) / Real(2) - k;
    if (utils::fuzzyEqual(xh * xh + yk * yk, radius * radius, epsilon)) {
      result.numIntersects = 1;
      result.t0 = Real(0);
    } else {
      result.numIntersects = 0;
    }
  } else {
    Real p0ShiftX = p0.x() - h;
    Real p0ShiftY = p0.y() - k;
    Real p1ShiftX = p1.x() - h;
    Real p1ShiftY = p1.y() - k;

    Real a;
    Real b;
    Real c;
    if (std::abs(dx) < epsilon) {
      // vertical line, using average x value for fuzziness
      Real xPos = (p0ShiftX + p1ShiftX) / Real(2);
      a = Real(1);
      b = Real(0);
      c = -xPos;
    } else {
      Real m = dy / dx;
      a = m;
      b = Real(-1);
      c = p1ShiftY - m * p1ShiftX;
    }

    Real a2 = a * a;
    Real b2 = b * b;
    Real c2 = c * c;
    Real r2 = radius * radius;
    Real a2_b2 = a2 + b2;
    Real shortestDist = std::abs(c) / std::sqrt(a2_b2);

    if (shortestDist > radius + epsilon) {
      result.numIntersects = 0;
    } else {
      // adding h and k back to solution terms (shifting from origin back to real coordinates)
      Real x0 = -a * c / a2_b2 + h;
      Real y0 = -b * c / a2_b2 + k;

      if (utils::fuzzyEqual(shortestDist, radius, epsilon)) {
        result.numIntersects = 1;
        result.t0 = parametricFromPoint(Vector2<Real>(x0, y0));
      } else {
        Real d = r2 - c2 / a2_b2;
        Real mult = std::sqrt(std::abs(d / a2_b2));
        Real xSol1 = x0 + b * mult;
        Real xSol2 = x0 - b * mult;
        Real ySol1 = y0 - a * mult;
        Real ySol2 = y0 + a * mult;
        result.numIntersects = 2;
        result.t0 = parametricFromPoint(Vector2<Real>(xSol1, ySol1));
        result.t1 = parametricFromPoint(Vector2<Real>(xSol2, ySol2));
      }
    }
  }

  CAVC_ASSERT(result.numIntersects >= 0 && result.numIntersects <= 2, "invalid intersect count");
  return result;
}
} // namespace cavc

#endif // CAVC_INTRLINESEG2CIRCLE2_HPP
