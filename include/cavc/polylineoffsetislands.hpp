#ifndef CAVC_POLYLINEOFFSETISLANDS_HPP
#define CAVC_POLYLINEOFFSETISLANDS_HPP
#include "polyline.hpp"
#include "polylinecombine.hpp"
#include "polylineoffset.hpp"
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <vector>

namespace cavc {
template <typename Real> struct OffsetLoop {
  std::size_t parentLoopIndex;
  Polyline<Real> polyline;
  StaticSpatialIndex<Real> spatialIndex;
};

template <typename Real> struct OffsetLoopSet {
  std::vector<OffsetLoop<Real>> ccwLoops;
  std::vector<OffsetLoop<Real>> cwLoops;
};

enum class OffsetLoopRole {
  Outer = 0,
  Hole = 1,
};

constexpr std::size_t kNoParentOffsetLoop = std::numeric_limits<std::size_t>::max();

template <typename Real> struct OffsetLoopTopologyNode {
  OffsetLoopRole role = OffsetLoopRole::Outer;
  std::size_t sourceIndex = 0;
  std::size_t parentIndex = kNoParentOffsetLoop;
};

template <typename Real> void sortOffsetLoopsStable(std::vector<OffsetLoop<Real>> &loops) {
  auto absArea = [](OffsetLoop<Real> const &loop) { return std::abs(getArea(loop.polyline)); };

  std::stable_sort(loops.begin(), loops.end(),
                   [&](OffsetLoop<Real> const &lhs, OffsetLoop<Real> const &rhs) {
                     Real lhsAbsArea = absArea(lhs);
                     Real rhsAbsArea = absArea(rhs);
                     if (lhsAbsArea != rhsAbsArea) {
                       return lhsAbsArea > rhsAbsArea;
                     }

                     if (lhs.spatialIndex.minX() != rhs.spatialIndex.minX()) {
                       return lhs.spatialIndex.minX() < rhs.spatialIndex.minX();
                     }
                     if (lhs.spatialIndex.minY() != rhs.spatialIndex.minY()) {
                       return lhs.spatialIndex.minY() < rhs.spatialIndex.minY();
                     }
                     if (lhs.spatialIndex.maxX() != rhs.spatialIndex.maxX()) {
                       return lhs.spatialIndex.maxX() < rhs.spatialIndex.maxX();
                     }
                     if (lhs.spatialIndex.maxY() != rhs.spatialIndex.maxY()) {
                       return lhs.spatialIndex.maxY() < rhs.spatialIndex.maxY();
                     }

                     return lhs.parentLoopIndex < rhs.parentLoopIndex;
                   });
}

template <typename Real> void sortOffsetLoopSetStable(OffsetLoopSet<Real> &loopSet) {
  sortOffsetLoopsStable(loopSet.ccwLoops);
  sortOffsetLoopsStable(loopSet.cwLoops);
}

template <typename Real>
std::vector<OffsetLoopTopologyNode<Real>>
buildOffsetLoopTopology(OffsetLoopSet<Real> const &loopSet,
                        Real boundaryEpsilon = utils::realPrecision<Real>()) {
  CAVC_ASSERT(boundaryEpsilon >= Real(0), "boundaryEpsilon must be >= 0");

  struct FlattenedLoopRef {
    OffsetLoopRole role;
    std::size_t sourceIndex;
    OffsetLoop<Real> const *loop;
    Real absArea;
  };

  std::vector<FlattenedLoopRef> flattened;
  flattened.reserve(loopSet.ccwLoops.size() + loopSet.cwLoops.size());
  for (std::size_t i = 0; i < loopSet.ccwLoops.size(); ++i) {
    flattened.push_back({OffsetLoopRole::Outer, i, &loopSet.ccwLoops[i],
                         std::abs(getArea(loopSet.ccwLoops[i].polyline))});
  }

  for (std::size_t i = 0; i < loopSet.cwLoops.size(); ++i) {
    flattened.push_back({OffsetLoopRole::Hole, i, &loopSet.cwLoops[i],
                         std::abs(getArea(loopSet.cwLoops[i].polyline))});
  }

  std::stable_sort(flattened.begin(), flattened.end(),
                   [](FlattenedLoopRef const &lhs, FlattenedLoopRef const &rhs) {
                     if (lhs.absArea != rhs.absArea) {
                       return lhs.absArea > rhs.absArea;
                     }

                     if (lhs.loop->spatialIndex.minX() != rhs.loop->spatialIndex.minX()) {
                       return lhs.loop->spatialIndex.minX() < rhs.loop->spatialIndex.minX();
                     }
                     if (lhs.loop->spatialIndex.minY() != rhs.loop->spatialIndex.minY()) {
                       return lhs.loop->spatialIndex.minY() < rhs.loop->spatialIndex.minY();
                     }
                     if (lhs.loop->spatialIndex.maxX() != rhs.loop->spatialIndex.maxX()) {
                       return lhs.loop->spatialIndex.maxX() < rhs.loop->spatialIndex.maxX();
                     }
                     if (lhs.loop->spatialIndex.maxY() != rhs.loop->spatialIndex.maxY()) {
                       return lhs.loop->spatialIndex.maxY() < rhs.loop->spatialIndex.maxY();
                     }

                     if (lhs.role != rhs.role) {
                       return lhs.role == OffsetLoopRole::Outer;
                     }

                     return lhs.sourceIndex < rhs.sourceIndex;
                   });

  std::vector<OffsetLoopTopologyNode<Real>> result(flattened.size());
  for (std::size_t i = 0; i < flattened.size(); ++i) {
    result[i].role = flattened[i].role;
    result[i].sourceIndex = flattened[i].sourceIndex;
  }

  auto pointInsideExpandedLoopBounds = [&](FlattenedLoopRef const &loop_ref,
                                           Vector2<Real> const &point) {
    Real const min_x = loop_ref.loop->spatialIndex.minX() - boundaryEpsilon;
    Real const min_y = loop_ref.loop->spatialIndex.minY() - boundaryEpsilon;
    Real const max_x = loop_ref.loop->spatialIndex.maxX() + boundaryEpsilon;
    Real const max_y = loop_ref.loop->spatialIndex.maxY() + boundaryEpsilon;
    return point.x() >= min_x && point.x() <= max_x && point.y() >= min_y && point.y() <= max_y;
  };

  for (std::size_t i = 0; i < flattened.size(); ++i) {
    auto const &child = flattened[i];
    if (child.loop->polyline.size() == 0) {
      continue;
    }

    Vector2<Real> samplePoint = child.loop->polyline[0].pos();
    std::size_t parentIndex = kNoParentOffsetLoop;
    Real parentAbsArea = std::numeric_limits<Real>::max();

    for (std::size_t j = 0; j < flattened.size(); ++j) {
      if (i == j) {
        continue;
      }

      auto const &candidate = flattened[j];
      if (candidate.role == child.role) {
        continue;
      }
      if (candidate.absArea <= child.absArea) {
        continue;
      }
      if (!pointInsideExpandedLoopBounds(candidate, samplePoint)) {
        continue;
      }

      PointContainment containment =
          getPointContainment(candidate.loop->polyline, samplePoint, boundaryEpsilon);
      if (containment == PointContainment::Outside) {
        continue;
      }

      if (candidate.absArea < parentAbsArea) {
        parentIndex = j;
        parentAbsArea = candidate.absArea;
      }
    }

    result[i].parentIndex = parentIndex;
  }

  return result;
}

template <typename Real> class ParallelOffsetIslands {
public:
  ParallelOffsetIslands() {}
  OffsetLoopSet<Real> compute(OffsetLoopSet<Real> const &input, Real offsetDelta);

private:
  // type to represent a slice point (intersect) on an OffsetLoop
  struct LoopSlicePoint {
    // intersect between loops
    PlineIntersect<Real> intr;
    bool noSliceAfterForIndex1;
  };

  // type to represent a set of slice points (intersects) between two loops
  struct SlicePointSet {
    // index of first loop involved
    std::size_t loopIndex1;
    // index of second loop involved
    std::size_t loopIndex2;
    // all of the slice points (intersects) between the two loops
    std::vector<LoopSlicePoint> slicePoints;
  };

  // get an offset loop by index i, maps to ccw then cw offset loops
  OffsetLoop<Real> &getOffsetLoop(std::size_t i) {
    return i < m_ccwOffsetLoops.size() ? m_ccwOffsetLoops[i]
                                       : m_cwOffsetLoops[i - m_ccwOffsetLoops.size()];
  }

  OffsetLoop<Real> const &getParentLoop(std::size_t i) {
    return i < m_inputSet->ccwLoops.size() ? m_inputSet->ccwLoops[i]
                                           : m_inputSet->cwLoops[i - m_inputSet->ccwLoops.size()];
  }

  void createOffsetLoops(const OffsetLoopSet<Real> &input, Real absDelta);
  void createOffsetLoopsIndex();
  void createSlicePoints();

  struct DissectionPoint {
    std::size_t otherLoopIndex;
    Vector2<Real> pos;
  };

  struct DissectedSlice {
    // open polyline representing the slice
    Polyline<Real> pline;
    // index of the loop the slice is from
    std::size_t sliceParentIndex;
    // index of the loop that intersected the parent loop to form the start of the slice
    std::size_t startLoopIndex;
    // index of the loop that intersected the parent loop to form the end of the slice
    std::size_t endLoopIndex;
  };

  bool pointOnOffsetValid(std::size_t skipIndex, Vector2<Real> const &pt, Real absDelta);
  void createSlicesFromLoop(std::size_t loopIndex, Real absDelta,
                            std::vector<DissectedSlice> &result);

  OffsetLoopSet<Real> const *m_inputSet;

  // counter clockwise offset loops, these surround the clockwise offset loops
  std::vector<OffsetLoop<Real>> m_ccwOffsetLoops;
  // clockwise (island) offset loops, these are surrounded by the counter clockwise loops
  std::vector<OffsetLoop<Real>> m_cwOffsetLoops;
  std::size_t totalOffsetLoopsCount() { return m_ccwOffsetLoops.size() + m_cwOffsetLoops.size(); }
  // spatial index of all the offset loops
  std::unique_ptr<StaticSpatialIndex<Real>> m_offsetLoopsIndex;
  using IndexPair = std::pair<std::size_t, std::size_t>;
  // set to keep track of already visited pairs of loops when finding intersects
  std::unordered_set<IndexPair, internal::IndexPairHash> m_visitedLoopPairs;
  // buffers to use for querying spatial indexes
  std::vector<std::size_t> m_queryStack;
  std::vector<std::size_t> m_queryResults;
  // slice point sets from intersects between loops
  std::vector<SlicePointSet> m_slicePointSets;
  // lookup used to get slice points for a particular loop index (holds indexes to sets in
  // m_slicePointSets)
  std::vector<std::vector<std::size_t>> m_slicePointsLookup;
  // dissection points used to form slices for a particular loop in createSlicesFromLoop
  std::unordered_map<std::size_t, std::vector<DissectionPoint>> m_loopDissectionPoints;
};

template <typename Real>
void ParallelOffsetIslands<Real>::createOffsetLoops(const OffsetLoopSet<Real> &input,
                                                    Real absDelta) {
  std::size_t const expected_loop_count = input.ccwLoops.size() + input.cwLoops.size();
  // create counter clockwise offset loops
  m_ccwOffsetLoops.clear();
  m_ccwOffsetLoops.reserve(expected_loop_count);
  m_cwOffsetLoops.clear();
  m_cwOffsetLoops.reserve(expected_loop_count);

  auto appendOffsetsWithBatchIndex = [&](std::vector<Polyline<Real>> &offsets,
                                         std::size_t parent_index,
                                         std::vector<OffsetLoop<Real>> &target) {
    if (offsets.size() == 0) {
      return;
    }

    std::vector<StaticSpatialIndex<Real>> indexes = createApproxSpatialIndices(offsets);
    CAVC_ASSERT(indexes.size() == offsets.size(), "batch index size mismatch");

    target.reserve(target.size() + offsets.size());
    for (std::size_t i = 0; i < offsets.size(); ++i) {
      target.push_back({parent_index, std::move(offsets[i]), std::move(indexes[i])});
    }
  };

  std::size_t parentIndex = 0;
  for (auto const &loop : input.ccwLoops) {
    auto offsets = parallelOffset(loop.polyline, absDelta);
    std::vector<Polyline<Real>> ccwOffsets;
    ccwOffsets.reserve(offsets.size());
    for (auto &offset : offsets) {
      // must check if orientation inverted (due to collapse of very narrow or small input)
      if (getArea(offset) < Real(0)) {
        continue;
      }
      ccwOffsets.push_back(std::move(offset));
    }

    appendOffsetsWithBatchIndex(ccwOffsets, parentIndex, m_ccwOffsetLoops);
    parentIndex += 1;
  }

  // create clockwise offset loops (note counter clockwise loops may result from outward offset)
  for (auto const &loop : input.cwLoops) {
    auto offsets = parallelOffset(loop.polyline, absDelta);
    std::vector<Polyline<Real>> cwOffsets;
    std::vector<Polyline<Real>> ccwOffsets;
    cwOffsets.reserve(offsets.size());
    ccwOffsets.reserve(offsets.size());
    for (auto &offset : offsets) {
      if (getArea(offset) < Real(0)) {
        cwOffsets.push_back(std::move(offset));
      } else {
        ccwOffsets.push_back(std::move(offset));
      }
    }

    appendOffsetsWithBatchIndex(cwOffsets, parentIndex, m_cwOffsetLoops);
    appendOffsetsWithBatchIndex(ccwOffsets, parentIndex, m_ccwOffsetLoops);
    parentIndex += 1;
  }
}

template <typename Real> void ParallelOffsetIslands<Real>::createOffsetLoopsIndex() {
  // create spatial index for all offset loop bounding boxes
  m_offsetLoopsIndex = std::make_unique<StaticSpatialIndex<Real>>(totalOffsetLoopsCount());
  for (auto const &posC : m_ccwOffsetLoops) {
    auto const &i = posC.spatialIndex;
    m_offsetLoopsIndex->add(i.minX(), i.minY(), i.maxX(), i.maxY());
  }

  for (auto const &negC : m_cwOffsetLoops) {
    auto const &i = negC.spatialIndex;
    m_offsetLoopsIndex->add(i.minX(), i.minY(), i.maxX(), i.maxY());
  }
  m_offsetLoopsIndex->finish();
}

template <typename Real> void ParallelOffsetIslands<Real>::createSlicePoints() {
  std::size_t const totalOffsetCount = totalOffsetLoopsCount();

  m_visitedLoopPairs.clear();
  m_visitedLoopPairs.reserve(totalOffsetCount * 2);
  m_slicePointSets.clear();
  m_slicePointSets.reserve(totalOffsetCount);
  m_slicePointsLookup.clear();
  m_queryResults.clear();
  m_queryResults.reserve(totalOffsetCount);

  // find all intersects between all offsets
  m_slicePointsLookup.resize(totalOffsetCount);
  PlineIntersectsResult<Real> intrsResults;
  for (std::size_t i = 0; i < totalOffsetCount; ++i) {
    auto const &loop1 = getOffsetLoop(i);
    auto const &index1 = loop1.spatialIndex;
    m_queryResults.clear();
    m_offsetLoopsIndex->query(index1.minX(), index1.minY(), index1.maxX(), index1.maxY(),
                              m_queryResults, m_queryStack);

    for (std::size_t j : m_queryResults) {
      // skip same index (no self intersects among the offset loops)
      if (i == j) {
        continue;
      }
      // skip reversed index order (would end up comparing the same loops)
      if (m_visitedLoopPairs.find({j, i}) != m_visitedLoopPairs.end()) {
        continue;
      }
      m_visitedLoopPairs.emplace(i, j);

      auto const &loop2 = getOffsetLoop(j);
      intrsResults.intersects.clear();
      intrsResults.coincidentIntersects.clear();
      // finding intersects
      findIntersects(loop1.polyline, loop2.polyline, index1, intrsResults);
      if (intrsResults.hasIntersects()) {
        m_slicePointSets.emplace_back();
        auto &slicePointSet = m_slicePointSets.back();
        slicePointSet.loopIndex1 = i;
        slicePointSet.loopIndex2 = j;
        for (auto &intr : intrsResults.intersects) {
          slicePointSet.slicePoints.push_back({std::move(intr), false});
        }

        // add coincident start and end points
        if (intrsResults.coincidentIntersects.size() != 0) {
          auto coinSliceResult = sortAndjoinCoincidentSlices(intrsResults.coincidentIntersects,
                                                             loop1.polyline, loop2.polyline);
          for (auto &sp : coinSliceResult.sliceStartPoints) {
            slicePointSet.slicePoints.push_back({std::move(sp), false});
          }
          for (auto &ep : coinSliceResult.sliceEndPoints) {
            slicePointSet.slicePoints.push_back({std::move(ep), true});
          }
        }

        m_slicePointsLookup[i].push_back(m_slicePointSets.size() - 1);
        m_slicePointsLookup[j].push_back(m_slicePointSets.size() - 1);
      }
    }
  }
}

template <typename Real>
bool ParallelOffsetIslands<Real>::pointOnOffsetValid(std::size_t skipIndex, const Vector2<Real> &pt,
                                                     Real absDelta) {
  // test distance against input polylines
  std::size_t const inputTotalCount = m_inputSet->ccwLoops.size() + m_inputSet->cwLoops.size();
  for (std::size_t i = 0; i < inputTotalCount; ++i) {
    if (i == skipIndex) {
      continue;
    }
    auto const &parentLoop = getParentLoop(i);
    if (!internal::pointValidForOffset(parentLoop.polyline, absDelta, parentLoop.spatialIndex, pt,
                                       m_queryStack)) {
      return false;
    }
  }

  return true;
}

template <typename Real>
void ParallelOffsetIslands<Real>::createSlicesFromLoop(std::size_t loopIndex, Real absDelta,
                                                       std::vector<DissectedSlice> &result) {
  OffsetLoop<Real> const &offsetLoop = getOffsetLoop(loopIndex);
  std::size_t const parentIndex = offsetLoop.parentLoopIndex;
  Polyline<Real> const &pline = offsetLoop.polyline;
  m_loopDissectionPoints.clear();
  for (auto const &setIndex : m_slicePointsLookup[loopIndex]) {
    auto const &set = m_slicePointSets[setIndex];
    bool isFirstIndex = loopIndex == set.loopIndex1;
    if (isFirstIndex) {
      for (auto const &p : set.slicePoints) {
        m_loopDissectionPoints[p.intr.sIndex1].push_back({set.loopIndex2, p.intr.pos});
      }
    } else {
      for (auto const &p : set.slicePoints) {
        m_loopDissectionPoints[p.intr.sIndex2].push_back({set.loopIndex1, p.intr.pos});
      }
    }
  }

  // sort points by distance from start vertex
  for (auto &kvp : m_loopDissectionPoints) {
    Vector2<Real> startPos = pline[kvp.first].pos();
    auto cmp = [&](DissectionPoint const &p1, DissectionPoint const &p2) {
      return distSquared(p1.pos, startPos) < distSquared(p2.pos, startPos);
    };
    std::sort(kvp.second.begin(), kvp.second.end(), cmp);
  }

  for (auto const &kvp : m_loopDissectionPoints) {
    // start index for the slice we're about to build
    std::size_t sIndex = kvp.first;
    // self intersect list for this start index
    std::vector<DissectionPoint> const &intrsList = kvp.second;

    const auto &firstSegStartVertex = pline[sIndex];
    std::size_t nextIndex = utils::nextWrappingIndex(sIndex, pline);
    const auto &firstSegEndVertex = pline[nextIndex];

    if (intrsList.size() != 1) {
      // build all the segments between the N intersects in siList (N > 1), skipping the first
      // segment (to be processed at the end)
      SplitResult<Real> firstSplit =
          splitAtPoint(firstSegStartVertex, firstSegEndVertex, intrsList[0].pos);
      auto prevVertex = firstSplit.splitVertex;
      for (std::size_t i = 1; i < intrsList.size(); ++i) {
        std::size_t const sliceStartIndex = intrsList[i - 1].otherLoopIndex;
        std::size_t const sliceEndIndex = intrsList[i].otherLoopIndex;
        SplitResult<Real> split = splitAtPoint(prevVertex, firstSegEndVertex, intrsList[i].pos);
        // update prevVertex for next loop iteration
        prevVertex = split.splitVertex;

        if (fuzzyEqual(split.updatedStart.pos(), split.splitVertex.pos(),
                       utils::realPrecision<Real>())) {
          continue;
        }

        auto sMidpoint = segMidpoint(split.updatedStart, split.splitVertex);
        if (!pointOnOffsetValid(parentIndex, sMidpoint, absDelta)) {
          // skip slice
          continue;
        }

        result.emplace_back();
        result.back().pline.addVertex(split.updatedStart);
        result.back().pline.addVertex(split.splitVertex);
        result.back().sliceParentIndex = loopIndex;
        result.back().startLoopIndex = sliceStartIndex;
        result.back().endLoopIndex = sliceEndIndex;
      }
    }

    // build the segment between the last intersect in instrsList and the next intersect found
    SplitResult<Real> split =
        splitAtPoint(firstSegStartVertex, firstSegEndVertex, intrsList.back().pos);

    DissectedSlice currSlice;
    currSlice.pline.addVertex(split.splitVertex);
    currSlice.sliceParentIndex = loopIndex;
    currSlice.startLoopIndex = intrsList.back().otherLoopIndex;

    std::size_t index = nextIndex;
    std::size_t loopCount = 0;
    const std::size_t maxLoopCount = pline.size();
    while (true) {
      if (loopCount++ > maxLoopCount) {
        CAVC_ASSERT(false, "Bug detected, should never loop this many times!");
        // break to avoid infinite loop
        break;
      }
      // add vertex
      internal::addOrReplaceIfSamePos(currSlice.pline, pline[index]);

      // check if segment that starts at vertex we just added has an intersect
      auto nextIntr = m_loopDissectionPoints.find(index);
      if (nextIntr != m_loopDissectionPoints.end()) {
        // there is an intersect, slice is done
        Vector2<Real> const &intersectPos = nextIntr->second[0].pos;

        // trim last added vertex and add final intersect position
        PlineVertex<Real> endVertex = PlineVertex<Real>(intersectPos, Real(0));
        std::size_t l_nextIndex = utils::nextWrappingIndex(index, pline);
        SplitResult<Real> l_split =
            splitAtPoint(currSlice.pline.lastVertex(), pline[l_nextIndex], intersectPos);
        currSlice.pline.lastVertex() = l_split.updatedStart;
        internal::addOrReplaceIfSamePos(currSlice.pline, endVertex);
        currSlice.endLoopIndex = nextIntr->second[0].otherLoopIndex;
        break;
      }
      // else there is not an intersect, increment index and continue
      index = utils::nextWrappingIndex(index, pline);
    }

    if (currSlice.pline.size() > 1) {
      auto sMidpoint =
          segMidpoint(currSlice.pline[currSlice.pline.size() - 2], currSlice.pline.lastVertex());
      if (!pointOnOffsetValid(parentIndex, sMidpoint, absDelta)) {
        // skip slice
        continue;
      }
      result.push_back(std::move(currSlice));
    }
  }
}

template <typename Real>
OffsetLoopSet<Real> ParallelOffsetIslands<Real>::compute(const OffsetLoopSet<Real> &input,
                                                         Real offsetDelta) {
  m_inputSet = &input;
  OffsetLoopSet<Real> result;
  Real absDelta = std::abs(offsetDelta);
  createOffsetLoops(input, absDelta);
  if (totalOffsetLoopsCount() == 0) {
    return result;
  }

  createOffsetLoopsIndex();
  createSlicePoints();

  std::vector<DissectedSlice> slices;
  std::size_t totalOffsetsCount = totalOffsetLoopsCount();
  slices.reserve(totalOffsetsCount * 2);

  std::vector<Polyline<Real>> resultSlices;
  resultSlices.reserve(totalOffsetsCount * 2);

  for (std::size_t i = 0; i < totalOffsetsCount; ++i) {
    if (m_slicePointsLookup[i].size() == 0) {
      // no intersects but still must test distance of one vertex position since it may be inside
      // another offset (completely eclipsed by island offset)
      auto &loop = getOffsetLoop(i);
      if (!pointOnOffsetValid(loop.parentLoopIndex, loop.polyline[0].pos(), absDelta)) {
        continue;
      }
      if (i < m_ccwOffsetLoops.size()) {
        result.ccwLoops.push_back(std::move(loop));
      } else {
        result.cwLoops.push_back(std::move(loop));
      }
      continue;
    }
    createSlicesFromLoop(i, absDelta, slices);
  }

  for (auto &slice : slices) {
    resultSlices.push_back(std::move(slice.pline));
  }

  std::vector<Polyline<Real>> stitched =
      internal::stitchOrderedSlicesIntoClosedPolylines(resultSlices);
  result.ccwLoops.reserve(result.ccwLoops.size() + stitched.size());
  result.cwLoops.reserve(result.cwLoops.size() + stitched.size());

  std::vector<Polyline<Real>> stitchedCcwLoops;
  std::vector<Polyline<Real>> stitchedCwLoops;
  stitchedCcwLoops.reserve(stitched.size());
  stitchedCwLoops.reserve(stitched.size());

  for (auto &r : stitched) {
    Real area = getArea(r);
    if (std::abs(area) < 1e-4) {
      continue;
    }
    if (area < Real(0)) {
      stitchedCwLoops.push_back(std::move(r));
    } else {
      stitchedCcwLoops.push_back(std::move(r));
    }
  }

  auto appendStitchedWithBatchIndex = [&](std::vector<Polyline<Real>> &loops,
                                          std::vector<OffsetLoop<Real>> &target) {
    if (loops.size() == 0) {
      return;
    }

    std::vector<StaticSpatialIndex<Real>> indexes = createApproxSpatialIndices(loops);
    CAVC_ASSERT(indexes.size() == loops.size(), "batch index size mismatch");

    target.reserve(target.size() + loops.size());
    for (std::size_t i = 0; i < loops.size(); ++i) {
      target.push_back({0, std::move(loops[i]), std::move(indexes[i])});
    }
  };

  appendStitchedWithBatchIndex(stitchedCcwLoops, result.ccwLoops);
  appendStitchedWithBatchIndex(stitchedCwLoops, result.cwLoops);

  sortOffsetLoopSetStable(result);
  return result;
}

} // namespace cavc

#endif // CAVC_POLYLINEOFFSETISLANDS_HPP
