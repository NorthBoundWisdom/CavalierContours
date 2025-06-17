# Changelog

## 0.1.1

- Add option to enable clipper library for benchmarks: `CAVC_ENABLE_CLIPPER_BENCHMARKS`. If enabled, clipper will be fetched and built automatically. default is disabled.
- Remove unused includes
- Rename of some variables declaration shadow the local/previous variable
  - if parameter shadows the local variable, the local variable is renamed to `p_<old_name>`
  - if variable in function shadows the local variable, the parameter is renamed to `l_<old_name>`
- Format code in tests
