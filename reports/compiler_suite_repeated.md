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
| baseline | 20.2223 | 0.2566 | 0.1677 | 1.000 | 0.000000e+00 | 0.000000e+00 | Cache-friendly baseline with std::array state storage and cached reconstruction. |
| compiler_o1 | 21.5709 | 0.3242 | 0.2118 | 1.067 | 0.000000e+00 | 0.000000e+00 | Compiler sensitivity study with -O1. |
| compiler_o2 | 20.2149 | 0.4819 | 0.3148 | 1.000 | 0.000000e+00 | 0.000000e+00 | Compiler sensitivity study with -O2. |
| compiler_o3_generic | 20.1428 | 0.2366 | 0.1546 | 0.996 | 0.000000e+00 | 0.000000e+00 | Compiler sensitivity study with -O3 and no -march=native. |
