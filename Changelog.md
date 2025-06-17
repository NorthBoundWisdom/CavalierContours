# Changelog

## 0.1.1

- Add option to enable clipper library for benchmarks: `CAVC_ENABLE_CLIPPER_BENCHMARKS`. If enabled, clipper will be fetched and built automatically. default is disabled.
- Remove unused includes
- Rename of some variables declaration shadow the local/previous variable
  - if parameter shadows the local variable, the local variable is renamed to `p_<old_name>`
  - if variable in function shadows the local variable, the parameter is renamed to `l_<old_name>`
- Format code in tests
- Move some test utilities to `testhelpers` and add some new cases(with figures in testhelpers/assets)
- Renane of previous tests:
  - TEST_capi_pline -> TEST_capi_pline
  - TEST_capi_pline_function -> TEST_capi_pline_function
  - TEST_capi_parallel_offset -> TEST_capi_parallel_offset
  - TEST_capi_combine_plines -> TEST_capi_combine_plines
  - TEST_sample->TEST_pline
- Add new tests:
  - TEST_plinesegment
  - TEST_pline
- polyline: handle empty input in function convertArcsToLines



