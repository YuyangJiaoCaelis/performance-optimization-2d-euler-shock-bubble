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
| baseline | 19.9651 | 0.2639 | 0.2986 | 1.000 | 0.000000e+00 | 0.000000e+00 | Cache-friendly baseline with std::array state storage and cached reconstruction. |
| vector_state | 382.0406 | 1.3219 | 1.4959 | 19.135 | 0.000000e+00 | 0.000000e+00 | Replace std::array state storage with std::vector. |
