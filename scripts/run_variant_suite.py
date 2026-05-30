#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import statistics
from pathlib import Path
import subprocess
import sys
import time

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from snapshot_tools import load_snapshot
from validate_outputs import compare_snapshots, summarise_snapshot, symmetry_metrics, y_invariant_metrics


BASE_FLAGS_NATIVE = "-O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic"
BASE_FLAGS_GENERIC = "-O3 -std=c++17 -Wall -Wextra -Wpedantic"


@dataclass(frozen=True)
class Variant:
    name: str
    description: str
    cppflags: tuple[str, ...] = ()
    cxxflags: str = BASE_FLAGS_NATIVE


VARIANTS = {
    "baseline": Variant(
        name="baseline",
        description="Cache-friendly baseline with std::array state storage and cached reconstruction.",
    ),
    "compiler_o1": Variant(
        name="compiler_o1",
        description="Compiler sensitivity study with -O1.",
        cxxflags="-O1 -std=c++17 -Wall -Wextra -Wpedantic",
    ),
    "compiler_o2": Variant(
        name="compiler_o2",
        description="Compiler sensitivity study with -O2.",
        cxxflags="-O2 -std=c++17 -Wall -Wextra -Wpedantic",
    ),
    "compiler_o3_generic": Variant(
        name="compiler_o3_generic",
        description="Compiler sensitivity study with -O3 and no -march=native.",
        cxxflags=BASE_FLAGS_GENERIC,
    ),
    "vector_state": Variant(
        name="vector_state",
        description="Replace std::array state storage with std::vector.",
        cppflags=("-DUSE_VECTOR_STATE_STORAGE",),
    ),
    "state_by_value": Variant(
        name="state_by_value",
        description="Pass individual state vectors by value into internal helpers.",
        cppflags=("-DPASS_STATE_BY_VALUE",),
    ),
    "grid_by_value": Variant(
        name="grid_by_value",
        description="Pass the full grid by value between sweep helpers.",
        cppflags=("-DPASS_GRID_BY_VALUE",),
    ),
    "cache_unfriendly": Variant(
        name="cache_unfriendly",
        description="Use a scrambled memory-access order when traversing the grid.",
        cppflags=("-DNON_CACHE_FRIENDLY_ACCESS",),
    ),
    "no_reconstruction_cache": Variant(
        name="no_reconstruction_cache",
        description="Recompute slope-limited half-step states at each interface instead of caching them.",
        cppflags=("-DDISABLE_RECONSTRUCTION_CACHE",),
    ),
}


def parse_variant_names(requested: str) -> list[Variant]:
    if requested == "all":
        return [VARIANTS[name] for name in VARIANTS]

    names: list[str] = []
    for name in requested.split(","):
        key = name.strip()
        if not key:
            continue
        if key not in VARIANTS:
            raise ValueError(f"Unknown variant '{key}'.")
        if key not in names:
            names.append(key)
    if not names:
        raise ValueError("No variants selected.")

    if "baseline" in names:
        names.remove("baseline")
    names.insert(0, "baseline")

    return [VARIANTS[name] for name in names]


def remove_previous_outputs(prefix: Path) -> None:
    for path in prefix.parent.glob(f"{prefix.name}_*.csv"):
        path.unlink()


def build_variant(variant: Variant) -> Path:
    target = Path("build") / "variants" / variant.name / "shock_bubble_solver"
    cppflags = " ".join(variant.cppflags)

    command = [
        "make",
        "-B",
        f"TARGET={target}",
        f"CPPFLAGS={cppflags}",
        f"CXXFLAGS={variant.cxxflags}",
    ]
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
    return PROJECT_ROOT / target


def run_solver_once(
    binary: Path,
    output_prefix: Path,
    nx: int,
    ny: int,
    final_time: float,
    initial_case: str,
) -> tuple[Path, float]:
    remove_previous_outputs(output_prefix)

    command = [
        str(binary),
        "--case",
        initial_case,
        "--nx",
        str(nx),
        "--ny",
        str(ny),
        "--final-time",
        f"{final_time:.12g}",
        "--output-prefix",
        str(output_prefix),
    ]

    start = time.perf_counter()
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
    wall_time = time.perf_counter() - start

    snapshots = sorted(output_prefix.parent.glob(f"{output_prefix.name}_*.csv"))
    if not snapshots:
        raise RuntimeError(f"No snapshots were produced for prefix {output_prefix}.")
    return snapshots[-1], wall_time


def remove_snapshot_family(snapshot_path: Path) -> None:
    stem = snapshot_path.stem
    prefix = stem.rsplit("_", 1)[0]
    for path in snapshot_path.parent.glob(f"{prefix}_*.csv"):
        path.unlink()


def timing_summary(samples: list[float]) -> dict[str, float]:
    mean = statistics.fmean(samples)
    std = statistics.stdev(samples) if len(samples) > 1 else 0.0
    median = statistics.median(samples)
    minimum = min(samples)
    maximum = max(samples)
    ci95_halfwidth = 1.96 * std / (len(samples) ** 0.5) if len(samples) > 1 else 0.0
    cv_percent = 100.0 * std / max(mean, 1.0e-30)
    return {
        "wall_time_mean_s": mean,
        "wall_time_std_s": std,
        "wall_time_cv_percent": cv_percent,
        "wall_time_median_s": median,
        "wall_time_min_s": minimum,
        "wall_time_max_s": maximum,
        "ci95_halfwidth_s": ci95_halfwidth,
        "best_wall_time_s": minimum,
        "wall_time_s": mean,
    }


def run_round_robin_campaign(
    variants: list[Variant],
    binaries: dict[str, Path],
    output_root: Path,
    nx: int,
    ny: int,
    final_time: float,
    initial_case: str,
    repetitions: int,
    warmup_runs: int,
) -> dict[str, dict[str, object]]:
    results: dict[str, dict[str, object]] = {
        variant.name: {"timing_samples_s": [], "snapshot_path": None} for variant in variants
    }

    for warmup_index in range(warmup_runs):
        for variant in variants:
            prefix = output_root / f"{variant.name}_warmup{warmup_index:02d}"
            snapshot_path, _ = run_solver_once(
                binaries[variant.name],
                prefix,
                nx,
                ny,
                final_time,
                initial_case,
            )
            remove_snapshot_family(snapshot_path)

    for repetition_index in range(repetitions):
        for variant in variants:
            prefix = output_root / f"{variant.name}_run{repetition_index:02d}"
            snapshot_path, wall_time = run_solver_once(
                binaries[variant.name],
                prefix,
                nx,
                ny,
                final_time,
                initial_case,
            )
            results[variant.name]["timing_samples_s"].append(wall_time)
            if repetition_index == repetitions - 1:
                results[variant.name]["snapshot_path"] = snapshot_path
            else:
                remove_snapshot_family(snapshot_path)

    return results


def format_metric(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6e}"


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return path.name


def write_csv_report(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_markdown_report(
    path: Path,
    verification: dict[str, float],
    rows: list[dict[str, object]],
    performance_case: tuple[int, int, float],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as handle:
        handle.write("# Variant Suite Report\n\n")
        handle.write(
            f"Performance case: `{performance_case[0]} x {performance_case[1]}`, "
            f"`final_time = {performance_case[2]:.6g}`\n\n"
        )
        handle.write(
            "Reported wall times are means over the timed repetitions. "
            "One warmup run per variant is discarded before timing, variants are executed in round-robin order, "
            "and confidence intervals use a normal 95% half-width estimate.\n\n"
        )
        handle.write("## Baseline verification\n\n")
        handle.write(
            f"- Planar-shock density Linf y-invariance: `{verification['planar_density_linf']:.6e}`\n"
        )
        handle.write(
            f"- Planar-shock pressure Linf y-invariance: `{verification['planar_pressure_linf']:.6e}`\n"
        )
        handle.write(
            f"- Shock-bubble density Linf symmetry: `{verification['bubble_density_linf']:.6e}`\n"
        )
        handle.write(
            f"- Shock-bubble pressure Linf symmetry: `{verification['bubble_pressure_linf']:.6e}`\n\n"
        )
        handle.write("## Performance and regression table\n\n")
        handle.write(
            "| Variant | Mean time (s) | Std. dev. (s) | 95% CI (s) | Time / baseline | "
            "Density Linf vs baseline | Pressure Linf vs baseline | Description |\n"
        )
        handle.write(
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |\n"
        )
        for row in rows:
            handle.write(
                f"| {row['name']} | {row['wall_time_s']:.4f} | {row['wall_time_std_s']:.4f} | "
                f"{row['ci95_halfwidth_s']:.4f} | "
                f"{row['time_ratio']:.3f} | "
                f"{format_metric(row['density_linf_vs_baseline'])} | "
                f"{format_metric(row['pressure_linf_vs_baseline'])} | "
                f"{row['description']} |\n"
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build, run, time, and compare solver variants."
    )
    parser.add_argument(
        "--variants",
        default="all",
        help="Comma-separated variant names, or 'all'.",
    )
    parser.add_argument("--nx", type=int, default=120, help="Performance-study grid size in x.")
    parser.add_argument("--ny", type=int, default=48, help="Performance-study grid size in y.")
    parser.add_argument(
        "--final-time",
        type=float,
        default=4.0e-5,
        help="Performance-study final time.",
    )
    parser.add_argument("--verify-nx", type=int, default=80, help="Verification grid size in x.")
    parser.add_argument("--verify-ny", type=int, default=32, help="Verification grid size in y.")
    parser.add_argument(
        "--verify-final-time",
        type=float,
        default=1.0e-5,
        help="Verification final time for the planar-shock test.",
    )
    parser.add_argument(
        "--output-root",
        default="output/variant_suite",
        help="Root directory for generated snapshots.",
    )
    parser.add_argument(
        "--report-prefix",
        default="reports/variant_suite",
        help="Prefix for CSV and Markdown reports.",
    )
    parser.add_argument(
        "--repetitions",
        type=int,
        default=1,
        help="Number of timed runs per variant.",
    )
    parser.add_argument(
        "--warmup-runs",
        type=int,
        default=1,
        help="Number of unmeasured warmup runs per variant before the timed round-robin campaign.",
    )
    args = parser.parse_args()
    if args.repetitions < 1:
        raise ValueError("--repetitions must be at least 1.")
    if args.warmup_runs < 0:
        raise ValueError("--warmup-runs must be non-negative.")

    variants = parse_variant_names(args.variants)
    output_root = PROJECT_ROOT / args.output_root
    output_root.mkdir(parents=True, exist_ok=True)

    binaries = {variant.name: build_variant(variant) for variant in variants}
    baseline_variant = variants[0]
    baseline_binary = binaries[baseline_variant.name]

    planar_prefix = output_root / "baseline_planar"
    planar_snapshot_path, _ = run_solver_once(
        baseline_binary,
        planar_prefix,
        args.verify_nx,
        args.verify_ny,
        args.verify_final_time,
        "planar-shock",
    )
    planar_snapshot = load_snapshot(planar_snapshot_path)
    planar_invariance = y_invariant_metrics(planar_snapshot)

    campaign_results = run_round_robin_campaign(
        variants,
        binaries,
        output_root / "performance",
        args.nx,
        args.ny,
        args.final_time,
        "shock-bubble",
        args.repetitions,
        args.warmup_runs,
    )
    baseline_result = campaign_results[baseline_variant.name]
    baseline_snapshot_path = baseline_result["snapshot_path"]
    if baseline_snapshot_path is None:
        raise RuntimeError("Baseline campaign did not produce a final snapshot.")
    baseline_summary = timing_summary(baseline_result["timing_samples_s"])
    bubble_snapshot_path = baseline_snapshot_path
    bubble_snapshot = load_snapshot(bubble_snapshot_path)
    bubble_symmetry = symmetry_metrics(bubble_snapshot)

    verification_issues, _ = summarise_snapshot(planar_snapshot)
    bubble_issues, _ = summarise_snapshot(bubble_snapshot)
    if verification_issues or bubble_issues:
        raise RuntimeError(
            "Baseline verification failed: "
            + "; ".join(verification_issues + bubble_issues)
        )

    verification = {
        "planar_density_linf": planar_invariance["density"]["linf"],
        "planar_pressure_linf": planar_invariance["pressure"]["linf"],
        "bubble_density_linf": bubble_symmetry["density"]["linf"],
        "bubble_pressure_linf": bubble_symmetry["pressure"]["linf"],
    }

    rows: list[dict[str, object]] = []
    baseline_mean_time = float(baseline_summary["wall_time_mean_s"])
    rows.append(
        {
            "name": baseline_variant.name,
            "description": baseline_variant.description,
            **baseline_summary,
            "timing_samples_s": ";".join(f"{value:.9f}" for value in baseline_result["timing_samples_s"]),
            "repetitions": args.repetitions,
            "warmup_runs": args.warmup_runs,
            "time_ratio": 1.0,
            "density_linf_vs_baseline": 0.0,
            "pressure_linf_vs_baseline": 0.0,
            "velocity_x_linf_vs_baseline": 0.0,
            "velocity_y_linf_vs_baseline": 0.0,
            "binary": display_path(binaries[baseline_variant.name]),
            "snapshot": display_path(Path(bubble_snapshot_path)),
            "cppflags": " ".join(baseline_variant.cppflags),
            "cxxflags": baseline_variant.cxxflags,
        }
    )

    for variant in variants[1:]:
        result = campaign_results[variant.name]
        snapshot_path = result["snapshot_path"]
        if snapshot_path is None:
            raise RuntimeError(f"Variant {variant.name} did not produce a final snapshot.")
        summary = timing_summary(result["timing_samples_s"])
        snapshot = load_snapshot(snapshot_path)
        comparison = compare_snapshots(snapshot, bubble_snapshot)

        rows.append(
            {
                "name": variant.name,
                "description": variant.description,
                **summary,
                "timing_samples_s": ";".join(f"{value:.9f}" for value in result["timing_samples_s"]),
                "repetitions": args.repetitions,
                "warmup_runs": args.warmup_runs,
                "time_ratio": float(summary["wall_time_mean_s"]) / baseline_mean_time,
                "density_linf_vs_baseline": comparison["density"]["linf"],
                "pressure_linf_vs_baseline": comparison["pressure"]["linf"],
                "velocity_x_linf_vs_baseline": comparison["velocity_x"]["linf"],
                "velocity_y_linf_vs_baseline": comparison["velocity_y"]["linf"],
                "binary": display_path(binaries[variant.name]),
                "snapshot": display_path(Path(snapshot_path)),
                "cppflags": " ".join(variant.cppflags),
                "cxxflags": variant.cxxflags,
            }
        )

    report_prefix = PROJECT_ROOT / args.report_prefix
    write_csv_report(report_prefix.with_suffix(".csv"), rows)
    write_markdown_report(
        report_prefix.with_suffix(".md"),
        verification,
        rows,
        (args.nx, args.ny, args.final_time),
    )

    print("Baseline verification:")
    print(f"  Planar density Linf y-invariance: {verification['planar_density_linf']:.6e}")
    print(f"  Planar pressure Linf y-invariance: {verification['planar_pressure_linf']:.6e}")
    print(f"  Bubble density Linf symmetry: {verification['bubble_density_linf']:.6e}")
    print(f"  Bubble pressure Linf symmetry: {verification['bubble_pressure_linf']:.6e}")
    print()

    print("Variant results:")
    for row in rows:
        print(
            f"  {row['name']}: time = {row['wall_time_s']:.4f} +/- {row['wall_time_std_s']:.4f} s "
            f"(95% CI +/- {row['ci95_halfwidth_s']:.4f} s), "
            f"time/baseline = {row['time_ratio']:.3f}, "
            f"density Linf = {row['density_linf_vs_baseline']:.6e}, "
            f"pressure Linf = {row['pressure_linf_vs_baseline']:.6e}"
        )

    print()
    print(f"Wrote {report_prefix.with_suffix('.csv')}")
    print(f"Wrote {report_prefix.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
