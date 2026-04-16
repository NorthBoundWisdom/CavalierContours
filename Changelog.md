# Changelog

## Unreleased

- Add runtime-configurable geometry tolerances in `mathutils`:
  - `setEpsilonConfig`, `getEpsilonConfig`, `resetEpsilonConfig`
  - preserve existing defaults (`1e-8`, `1e-5`, `1e-4`, `1e-4`)
- Add point containment classification in C++ API:
  - `PointContainment` enum
  - `getPointContainment` (Outside/Inside/OnBoundary)
- Add polyline preprocessing in C++ API:
  - `ClosedPolylineWinding` enum
  - `normalizePolyline` (prune singularities + optional closed winding normalization)
  - `removeRedundant` and C API wrapper `cavc_pline_remove_redundant`
  - run parallel offset on redundancy-cleaned input to match the newer preprocessing flow
    for dense/near-degenerate real-world polylines
- Expand C API surface with low-risk geometry helpers:
  - `cavc_pline_invert_direction`
  - `cavc_pline_prune_singularities`
  - `cavc_pline_remove_redundant`
  - `cavc_pline_convert_arcs_to_lines`
  - `cavc_pline_normalize`
  - `cavc_spatial_index_create/delete/item_count/query_count/query`
  - `cavc_get_point_containment`
  - `cavc_get_tolerances`, `cavc_set_tolerances`, `cavc_reset_tolerances`
- Expand C API test coverage:
  - add `TEST_cavc_api_regression` for legacy mutating chains and new API regressions
  - add regression coverage for real-world self-intersecting shape-offset input cleanup and
    in-place `cavc_pline_remove_redundant`
  - fix path-length skip guard in `TEST_cavc_pline_function`
  - harden `plineFromVertexes` helper to safely support empty input
  - refactor `combine_with_self_invariants` into closed-only parameterized suite (remove stale
    runtime skips on open half-circle cases)
- Unify parallel offset options shape for join/end-cap controls:
  - C++: add `ParallelOffsetOptions`, `OffsetJoinType`, `OffsetEndCapType`
  - C API: replace offset bit flags with `cavc_parallel_offset_options`
  - add `cavc_parallel_offset_default_options`
- Clarify `miter_limit` validation semantics:
  - enforce `miter_limit >= 1` only when `join_type`/`joinType` is miter
- Improve runtime tolerance configuration reliability:
  - make tolerance reads/writes thread-safe for concurrent API usage
- Extend join behavior incrementally:
  - add `miter/bevel` for line-line joins in `parallelOffset`
  - add `miter_limit` option (C++ and C API) with line-line and arc-involved miter clipping
  - add arc-involved non-round joins using endpoint tangent intersection
    (collapsed-arc joins degrade to bevel)
  - add regression coverage for `CAVC_OFFSET_JOIN_MITER` and `CAVC_OFFSET_JOIN_BEVEL`, including
    arc-involved join scenarios
  - add arc-heavy join stress matrix coverage for high-curvature / near-tangent / degenerate
    repeated-vertex inputs (positive and negative deltas), including `miter_limit=1` bevel
    equivalence and high-limit miter divergence checks
- Extend end cap behavior incrementally:
  - add `square` end cap for open polylines with line start/end segments
    (endpoint extension + dual-slice cap-circle shift)
  - extend `square` end cap to arc start/end by preserving endpoint arc geometry and inserting
    explicit cap line segments when needed
  - separate `butt` from `round` for open polylines:
    `butt` uses endpoint-normal line clipping while `round` keeps endpoint-circle clipping
  - add regression coverage for `CAVC_OFFSET_END_CAP_SQUARE` on arc start/end for both
    positive and negative deltas
  - add regression coverage for `CAVC_OFFSET_END_CAP_BUTT` on arc start/end clip behavior
  - add stress matrix coverage for `CAVC_OFFSET_END_CAP_BUTT` on open polylines:
    high-curvature / near-tangent / multi-intersection, each with positive and negative offsets
  - add degenerate open-polyline end cap stress matrix across `round/square/butt`:
    repeated endpoints / near-zero-radius arcs / high-density self-intersection, each with positive
    and negative offsets
- Extend offset-islands engineering:
  - add stable geometry ordering for `ParallelOffsetIslands::compute` outputs
  - add loop role/tag and topology helpers:
    `OffsetLoopRole`, `OffsetLoopTopologyNode`, `buildOffsetLoopTopology`
  - expose topology output in C API with ownership-safe handle:
    `cavc_offset_loop_topology_build/delete/count/get`
  - add dedicated regression test target `TEST_cavc_offset_islands`
- Performance-oriented tuning (no behavior change):
  - reuse query buffers for open-polyline end-cap circle clipping in dual-slice offset path
  - pre-reserve frequently reused containers in offset-islands slicing/topology flow
  - add bounding-box early reject before expensive containment checks in topology parent search
  - add batch spatial-index helper `createApproxSpatialIndices`
  - switch `ParallelOffsetIslands` loop-index construction to batch path for
    offset loops and stitched loops
  - add `offsetislandsbenchmarks` for:
    `scalar vs batch index build`, large-scale `buildOffsetLoopTopology`,
    and large-scale `ParallelOffsetIslands::compute`

## 0.1.0

- Add option to enable clipper library for benchmarks: `CAVC_ENABLE_CLIPPER_BENCHMARKS`
- Rename variables whose declarations shadow local/previous variables
  - if parameter shadows the local variable, the local variable is renamed to `p_<old_name>`
  - if variable in function shadows the local variable, the parameter is renamed to `l_<old_name>`
