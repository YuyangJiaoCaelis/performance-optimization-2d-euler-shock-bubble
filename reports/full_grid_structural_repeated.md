# Variant Suite Report

Performance case: `500 x 197`, `final_time = 0.0003`

Reported wall times are means over the timed repetitions. One warmup run per variant is discarded before timing, variants are executed in round-robin order, and confidence intervals use a normal 95% half-width estimate.

## Baseline verification

- Planar-shock density Linf y-invariance: `0.000000e+00`
- Planar-shock pressure Linf y-invariance: `0.000000e+00`
- Shock-bubble density Linf symmetry: `0.000000e+00`
- Shock-bubble pressure Linf symmetry: `0.000000e+00`

## Performance and regression table

| Variant | Mean time (s) | Std. dev. (s) | 95% CI (s) | Time / baseline | Density Linf vs baseline | Pressure Linf vs baseline | Description |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| baseline | 20.0898 | 0.1937 | 0.1265 | 1.000 | 0.000000e+00 | 0.000000e+00 | Cache-friendly baseline with std::array state storage and cached reconstruction. |
| state_by_value | 20.3471 | 0.3032 | 0.1981 | 1.013 | 0.000000e+00 | 0.000000e+00 | Pass individual state vectors by value into internal helpers. |
| grid_by_value | 20.3005 | 0.3511 | 0.2294 | 1.010 | 0.000000e+00 | 0.000000e+00 | Pass the full grid by value between sweep helpers. |
| cache_unfriendly | 20.2698 | 0.1740 | 0.1137 | 1.009 | 0.000000e+00 | 0.000000e+00 | Use a scrambled memory-access order when traversing the grid. |
| no_reconstruction_cache | 33.5086 | 0.2667 | 0.1742 | 1.668 | 0.000000e+00 | 0.000000e+00 | Recompute slope-limited half-step states at each interface instead of caching them. |
