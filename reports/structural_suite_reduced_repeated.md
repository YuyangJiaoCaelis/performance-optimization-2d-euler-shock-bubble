# Variant Suite Report

Performance case: `160 x 63`, `final_time = 0.0003`

Reported wall times are means over the timed repetitions. One warmup run per variant is discarded before timing, variants are executed in round-robin order, and confidence intervals use a normal 95% half-width estimate.

## Baseline verification

- Planar-shock density Linf y-invariance: `0.000000e+00`
- Planar-shock pressure Linf y-invariance: `0.000000e+00`
- Shock-bubble density Linf symmetry: `0.000000e+00`
- Shock-bubble pressure Linf symmetry: `0.000000e+00`

## Performance and regression table

| Variant | Mean time (s) | Std. dev. (s) | 95% CI (s) | Time / baseline | Density Linf vs baseline | Pressure Linf vs baseline | Description |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| baseline | 0.6036 | 0.0464 | 0.0303 | 1.000 | 0.000000e+00 | 0.000000e+00 | Cache-friendly baseline with std::array state storage and cached reconstruction. |
| vector_state | 11.0088 | 0.2016 | 0.1317 | 18.237 | 0.000000e+00 | 0.000000e+00 | Replace std::array state storage with std::vector. |
| state_by_value | 0.5884 | 0.0044 | 0.0029 | 0.975 | 0.000000e+00 | 0.000000e+00 | Pass individual state vectors by value into internal helpers. |
| grid_by_value | 0.6027 | 0.0359 | 0.0234 | 0.999 | 0.000000e+00 | 0.000000e+00 | Pass the full grid by value between sweep helpers. |
| cache_unfriendly | 0.6000 | 0.0256 | 0.0167 | 0.994 | 0.000000e+00 | 0.000000e+00 | Use a scrambled memory-access order when traversing the grid. |
| no_reconstruction_cache | 1.0021 | 0.1026 | 0.0670 | 1.660 | 0.000000e+00 | 0.000000e+00 | Recompute slope-limited half-step states at each interface instead of caching them. |
