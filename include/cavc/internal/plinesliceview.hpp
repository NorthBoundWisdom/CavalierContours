#ifndef CAVC_INTERNAL_PLINESLICEVIEW_HPP
#define CAVC_INTERNAL_PLINESLICEVIEW_HPP

#include "../polyline.hpp"
#include <optional>
#include <utility>

namespace cavc {
namespace internal {

enum class PlineSliceViewValidation {
  Valid,
  SourceHasNoSegments,
  OffsetOutOfRange,
  UpdatedStartNotOnSegment,
  EndPointNotOnSegment,
  EndPointOnFinalOffsetVertex,
  UpdatedBulgeDoesNotMatch,
};

template <typename Real> struct SlicePointAtIndex {
  std::size_t sIndex;
  Vector2<Real> pos;
};

template <typename Real> struct PlineSliceViewData {
  std::size_t startIndex;
  std::size_t endIndexOffset;
  PlineVertex<Real> updatedStart;
  Real updatedEndBulge;
  Vector2<Real> endPoint;
  bool invertedDirection = false;

  std::size_t vertexCount() const { return endIndexOffset + 2; }

  PlineVertex<Real> firstVertex(Polyline<Real> const &source) const { return getVertex(source, 0); }

  PlineVertex<Real> lastVertex(Polyline<Real> const &source) const {
    return getVertex(source, vertexCount() - 1);
  }

  Vector2<Real> firstPoint(Polyline<Real> const &source) const { return firstVertex(source).pos(); }

  Vector2<Real> lastPoint(Polyline<Real> const &source) const { return lastVertex(source).pos(); }

  template <typename Visitor>
  void visitSegments(Polyline<Real> const &source, Visitor &&visitor) const {
    if (vertexCount() < 2) {
      return;
    }

    PlineVertex<Real> prevVertex = getVertex(source, 0);
    for (std::size_t i = 1; i < vertexCount(); ++i) {
      PlineVertex<Real> currVertex = getVertex(source, i);
      if (!visitor(prevVertex, currVertex)) {
        return;
      }
      prevVertex = currVertex;
    }
  }

  void appendTo(Polyline<Real> &target, Polyline<Real> const &source,
                Real joinThreshold = utils::sliceJoinThreshold<Real>()) const {
    target.vertexes().reserve(target.size() + vertexCount());
    for (std::size_t i = 0; i < vertexCount(); ++i) {
      addOrReplaceIfSamePos(target, getVertex(source, i), joinThreshold);
    }
  }

  Polyline<Real> toPolyline(Polyline<Real> const &source,
                            Real joinThreshold = utils::sliceJoinThreshold<Real>()) const {
    Polyline<Real> result;
    result.isClosed() = false;
    result.vertexes().reserve(vertexCount());
    appendTo(result, source, joinThreshold);
    return result;
  }

  PlineSliceViewValidation debugValidateForSource(Polyline<Real> const &source) const {
    if (source.size() < 2) {
      return PlineSliceViewValidation::SourceHasNoSegments;
    }

    if (endIndexOffset > source.size()) {
      return PlineSliceViewValidation::OffsetOutOfRange;
    }

    Real const validationEps = Real(1e-5);
    Real const pointOnSegEps = Real(1e-3);
    auto pointIsOnSegment = [&](std::size_t segIndex, Vector2<Real> const &point) {
      PlineVertex<Real> const &v1 = source[segIndex];
      PlineVertex<Real> const &v2 = source[utils::nextWrappingIndex(segIndex, source)];
      if (fuzzyEqual(point, v1.pos(), pointOnSegEps) || fuzzyEqual(point, v2.pos(), pointOnSegEps)) {
        return true;
      }

      Vector2<Real> closestPoint = closestPointOnSeg(v1, v2, point, validationEps);
      return fuzzyEqual(closestPoint, point, pointOnSegEps);
    };

    if (!pointIsOnSegment(startIndex, updatedStart.pos())) {
      return PlineSliceViewValidation::UpdatedStartNotOnSegment;
    }

    std::size_t const endIndex = fwdWrappingIndex(startIndex, endIndexOffset, source);
    if (!pointIsOnSegment(endIndex, endPoint)) {
      return PlineSliceViewValidation::EndPointNotOnSegment;
    }

    if (fuzzyEqual(endPoint, source[endIndex].pos(), validationEps)) {
      return PlineSliceViewValidation::EndPointOnFinalOffsetVertex;
    }

    if (endIndexOffset == 0 &&
        !utils::fuzzyEqual(updatedEndBulge, updatedStart.bulge(), validationEps)) {
      return PlineSliceViewValidation::UpdatedBulgeDoesNotMatch;
    }

    return PlineSliceViewValidation::Valid;
  }

  static std::optional<PlineSliceViewData>
  createOnSingleSegment(Polyline<Real> const &source, std::size_t startIndex,
                        PlineVertex<Real> const &updatedStart, Vector2<Real> const &endIntersect,
                        Real posEqualEps = utils::realPrecision<Real>()) {
    (void)source;
    if (fuzzyEqual(updatedStart.pos(), endIntersect, posEqualEps)) {
      return std::nullopt;
    }

    PlineSliceViewData result{startIndex, 0, updatedStart, updatedStart.bulge(), endIntersect,
                              false};
    CAVC_ASSERT(result.debugValidateForSource(source) == PlineSliceViewValidation::Valid,
                "slice view data is invalid for source");
    return result;
  }

  static PlineSliceViewData fromEntirePolyline(Polyline<Real> const &source) {
    CAVC_ASSERT(source.size() >= 2, "source must have at least 2 vertexes to form slice view");

    PlineSliceViewData result;
    result.startIndex = 0;
    result.updatedStart = source[0];
    result.invertedDirection = false;

    if (source.isClosed()) {
      result.endIndexOffset = source.size() - 1;
      result.updatedEndBulge = source.lastVertex().bulge();
      result.endPoint = source[0].pos();
    } else {
      result.endIndexOffset = source.size() - 2;
      result.updatedEndBulge = source[source.size() - 2].bulge();
      result.endPoint = source.lastVertex().pos();
    }

    CAVC_ASSERT(result.debugValidateForSource(source) == PlineSliceViewValidation::Valid,
                "slice view data is invalid for source");
    return result;
  }

  static std::optional<PlineSliceViewData>
  createFromSlicePoints(Polyline<Real> const &source, Vector2<Real> const &startPoint,
                        std::size_t startIndex, Vector2<Real> const &endPoint,
                        std::size_t endIndex, Real posEqualEps = utils::realPrecision<Real>()) {
    CAVC_ASSERT(source.size() >= 2, "source must have at least 2 vertexes to form slice view");
    CAVC_ASSERT(startIndex <= endIndex || source.isClosed(),
                "open polyline slice start index must be <= end index");

    auto adjustedStart = [&]() {
      if (!source.isClosed() && startIndex >= endIndex) {
        return std::make_pair(startIndex, false);
      }

      std::size_t nextIndex = utils::nextWrappingIndex(startIndex, source);
      if (fuzzyEqual(source[nextIndex].pos(), startPoint, posEqualEps)) {
        return std::make_pair(nextIndex, true);
      }

      return std::make_pair(startIndex, false);
    }();

    startIndex = adjustedStart.first;
    bool startPointAtSegEnd = adjustedStart.second;

    auto forwardDistance = [&](std::size_t from, std::size_t to) {
      if (from <= to) {
        return to - from;
      }
      return source.size() - from + to;
    };

    std::size_t traverseCount;
    std::size_t const indexDist = forwardDistance(startIndex, endIndex);
    if (indexDist == 0 && source.isClosed() && !fuzzyEqual(startPoint, endPoint, posEqualEps)) {
      Vector2<Real> const &segStart = source[startIndex].pos();
      Real distToStartPoint = distSquared(segStart, startPoint);
      Real distToEndPoint = distSquared(segStart, endPoint);
      traverseCount = distToStartPoint < distToEndPoint ? 0 : source.size();
    } else {
      traverseCount = indexDist;
    }

    PlineVertex<Real> const &startV1 = source[startIndex];
    PlineVertex<Real> const &startV2 = source[utils::nextWrappingIndex(startIndex, source)];

    PlineVertex<Real> updatedStart;
    if (startPointAtSegEnd) {
      if (traverseCount == 0) {
        updatedStart = splitAtPoint(startV1, startV2, endPoint).updatedStart;
      } else {
        updatedStart = startV1;
      }
    } else {
      SplitResult<Real> startSplit = splitAtPoint(startV1, startV2, startPoint);
      updatedStart = startSplit.splitVertex;
      if (traverseCount == 0) {
        updatedStart = splitAtPoint(updatedStart, startV2, endPoint).updatedStart;
      }
    }

    if (traverseCount == 0) {
      return createOnSingleSegment(source, startIndex, updatedStart, endPoint, posEqualEps);
    }

    PlineSliceViewData result =
        createMultiSegment(source, startIndex, endPoint, endIndex, updatedStart, traverseCount,
                           posEqualEps);
    CAVC_ASSERT(result.debugValidateForSource(source) == PlineSliceViewValidation::Valid,
                "slice view data is invalid for source");
    return result;
  }

private:
  static std::size_t fwdWrappingIndex(std::size_t startIndex, std::size_t offset,
                                      Polyline<Real> const &source) {
    CAVC_ASSERT(startIndex < source.size(), "start index out of bounds");
    CAVC_ASSERT(offset <= source.size(), "offset wraps multiple times");
    std::size_t sum = startIndex + offset;
    if (sum < source.size()) {
      return sum;
    }
    return sum - source.size();
  }

  static PlineSliceViewData createMultiSegment(Polyline<Real> const &source, std::size_t startIndex,
                                               Vector2<Real> const &endIntersect,
                                               std::size_t intersectIndex,
                                               PlineVertex<Real> const &updatedStart,
                                               std::size_t traverseCount, Real posEqualEps) {
    CAVC_ASSERT(traverseCount != 0,
                "traverseCount must be > 0, use createOnSingleSegment for same-segment slices");

    PlineVertex<Real> const &currentVertex = source[intersectIndex];

    std::size_t endIndexOffset;
    Real updatedEndBulge;
    if (fuzzyEqual(endIntersect, currentVertex.pos(), posEqualEps)) {
      endIndexOffset = traverseCount - 1;
      if (endIndexOffset != 0) {
        updatedEndBulge = source[utils::prevWrappingIndex(intersectIndex, source)].bulge();
      } else {
        updatedEndBulge = updatedStart.bulge();
      }
    } else {
      std::size_t nextIndex = utils::nextWrappingIndex(intersectIndex, source);
      SplitResult<Real> split = splitAtPoint(currentVertex, source[nextIndex], endIntersect);
      endIndexOffset = traverseCount;
      updatedEndBulge = split.updatedStart.bulge();
    }

    return {startIndex, endIndexOffset, updatedStart, updatedEndBulge, endIntersect, false};
  }

  PlineVertex<Real> getVertex(Polyline<Real> const &source, std::size_t index) const {
    CAVC_ASSERT(index < vertexCount(), "index out of range for slice view");

    if (invertedDirection) {
      if (index == 0) {
        return PlineVertex<Real>(endPoint, -updatedEndBulge);
      }

      if (index < endIndexOffset) {
        std::size_t bulgeIndex = fwdWrappingIndex(startIndex, endIndexOffset - index, source);
        std::size_t i = utils::nextWrappingIndex(bulgeIndex, source);
        PlineVertex<Real> result = source[i];
        result.bulge() = -source[bulgeIndex].bulge();
        return result;
      }

      if (index == endIndexOffset) {
        std::size_t i = fwdWrappingIndex(startIndex, endIndexOffset - index + 1, source);
        PlineVertex<Real> result = source[i];
        result.bulge() = -updatedStart.bulge();
        return result;
      }

      return PlineVertex<Real>(updatedStart.pos(), Real(0));
    }

    if (index == 0) {
      return updatedStart;
    }

    if (index < endIndexOffset) {
      return source[fwdWrappingIndex(startIndex, index, source)];
    }

    if (index == endIndexOffset) {
      PlineVertex<Real> result = source[fwdWrappingIndex(startIndex, endIndexOffset, source)];
      result.bulge() = updatedEndBulge;
      return result;
    }

    return PlineVertex<Real>(endPoint, Real(0));
  }
};

template <typename Real> struct OffsetSliceRef {
  std::size_t intrStartIndex;
  PlineSliceViewData<Real> viewData;
};

template <typename Real> struct CombineSliceRef {
  PlineSliceViewData<Real> viewData;
  bool sourceIsPlineA;
  bool coincident;
};

template <typename Real> struct CoincidentSliceInfo {
  SlicePointAtIndex<Real> startPointOnA;
  SlicePointAtIndex<Real> endPointOnA;
  SlicePointAtIndex<Real> startPointOnB;
  SlicePointAtIndex<Real> endPointOnB;
  bool opposingDirection;
};

} // namespace internal
} // namespace cavc

#endif // CAVC_INTERNAL_PLINESLICEVIEW_HPP
