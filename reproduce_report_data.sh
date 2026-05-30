#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"

rm -rf build output reports
mkdir -p output/benchmark_series output/grid_sensitivity output/planar_shock_reference reports

make TARGET=shock_bubble_solver

./shock_bubble_solver \
  --nx 500 --ny 197 --final-time 3.0e-4 \
  --output-prefix output/default_run

./shock_bubble_solver \
  --nx 500 --ny 197 --final-time 3.0e-4 \
  --snapshot-interval 3.75e-05 \
  --output-prefix output/benchmark_series/baseline

./shock_bubble_solver \
  --nx 250 --ny 99 --final-time 3.0e-4 \
  --output-prefix output/grid_sensitivity/nx250

./shock_bubble_solver \
  --nx 1000 --ny 394 --final-time 3.0e-4 \
  --output-prefix output/grid_sensitivity/nx1000

./shock_bubble_solver \
  --case planar-shock \
  --nx 500 --ny 197 --final-time 1.0e-5 \
  --output-prefix output/planar_shock_reference/planar

python3 scripts/validate_outputs.py \
  output/planar_shock_reference/planar_0001.csv \
  --expect-y-invariant \
  > reports/planar_shock_validation.txt

python3 scripts/validate_outputs.py \
  output/default_run_0001.csv \
  > reports/final_baseline_validation.txt

python3 scripts/planar_shock_reference_check.py \
  output/planar_shock_reference/planar_0001.csv \
  --final-time 1.0e-5 \
  --csv-output reports/planar_shock_reference_check.csv \
  --txt-output reports/planar_shock_reference_check.txt

python3 scripts/run_variant_suite.py \
  --variants baseline,compiler_o1,compiler_o2,compiler_o3_generic \
  --nx 500 --ny 197 --final-time 3.0e-4 \
  --repetitions 9 --warmup-runs 1 \
  --output-root output/compiler_suite_repeated \
  --report-prefix reports/compiler_suite_repeated

python3 scripts/run_variant_suite.py \
  --variants baseline,vector_state,state_by_value,grid_by_value,cache_unfriendly,no_reconstruction_cache \
  --nx 160 --ny 63 --final-time 3.0e-4 \
  --repetitions 9 --warmup-runs 1 \
  --output-root output/structural_suite_reduced_repeated \
  --report-prefix reports/structural_suite_reduced_repeated

python3 scripts/run_variant_suite.py \
  --variants baseline,state_by_value,grid_by_value,cache_unfriendly,no_reconstruction_cache \
  --nx 500 --ny 197 --final-time 3.0e-4 \
  --repetitions 9 --warmup-runs 1 \
  --output-root output/full_grid_structural_repeated \
  --report-prefix reports/full_grid_structural_repeated

python3 scripts/run_variant_suite.py \
  --variants baseline,vector_state \
  --nx 500 --ny 197 --final-time 3.0e-4 \
  --repetitions 3 --warmup-runs 1 \
  --output-root output/full_grid_vector_repeated \
  --report-prefix reports/full_grid_vector_repeated

python3 scripts/generate_report_artifacts.py
