#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import sys

from snapshot_tools import Snapshot, load_snapshot, relative_error_metrics


def field_metrics(differences: list[float], reference_scale: list[float]) -> dict[str, float]:
    abs_diff_sum = 0.0
    abs_ref_sum = 0.0
    diff_square_sum = 0.0
    ref_square_sum = 0.0
    max_diff = 0.0
    max_ref = 0.0

    for difference, reference in zip(differences, reference_scale):
        abs_difference = abs(difference)
        abs_reference = abs(reference)

        abs_diff_sum += abs_difference
        abs_ref_sum += abs_reference
        diff_square_sum += difference * difference
        ref_square_sum += reference * reference
        max_diff = max(max_diff, abs_difference)
        max_ref = max(max_ref, abs_reference)

    return {
        "l1": abs_diff_sum / max(abs_ref_sum, 1.0e-30),
        "l2": math.sqrt(diff_square_sum) / max(math.sqrt(ref_square_sum), 1.0e-30),
        "linf": max_diff / max(max_ref, 1.0e-30),
        "max_abs": max_diff,
    }


def summarise_snapshot(snapshot: Snapshot) -> tuple[list[str], list[str]]:
    issues: list[str] = []
    summary: list[str] = []

    density_values = snapshot.fields["density"]
    pressure_values = snapshot.fields["pressure"]
    velocity_x_values = snapshot.fields["velocity_x"]
    velocity_y_values = snapshot.fields["velocity_y"]

    finite_values = True
    for field_name, values in snapshot.fields.items():
        if any(not math.isfinite(value) for value in values):
            issues.append(f"{field_name}: found a non-finite value.")
            finite_values = False

    if finite_values:
        summary.append("All reported fields are finite.")

    minimum_density = min(density_values)
    minimum_pressure = min(pressure_values)
    if minimum_density <= 0.0:
        issues.append(f"density: minimum value is {minimum_density:.6e}.")
    if minimum_pressure <= 0.0:
        issues.append(f"pressure: minimum value is {minimum_pressure:.6e}.")

    summary.append(
        f"density range: [{minimum_density:.6e}, {max(density_values):.6e}]"
    )
    summary.append(
        f"pressure range: [{minimum_pressure:.6e}, {max(pressure_values):.6e}]"
    )
    summary.append(
        f"max |velocity_x|: {max(abs(value) for value in velocity_x_values):.6e}"
    )
    summary.append(
        f"max |velocity_y|: {max(abs(value) for value in velocity_y_values):.6e}"
    )

    if snapshot.ny > 1:
        coordinate_mismatch = max(
            abs(snapshot.y_values[j] + snapshot.y_values[snapshot.ny - 1 - j])
            for j in range(snapshot.ny // 2)
        )
        summary.append(f"top/bottom coordinate mismatch: {coordinate_mismatch:.6e}")

    return issues, summary


def symmetry_metrics(snapshot: Snapshot) -> dict[str, dict[str, float]]:
    pair_count = snapshot.ny // 2
    metrics: dict[str, dict[str, float]] = {}

    symmetric_fields = ("density", "pressure", "velocity_x")
    antisymmetric_fields = ("velocity_y",)

    for field_name in symmetric_fields:
        differences: list[float] = []
        scales: list[float] = []
        for j in range(pair_count):
            mirror_j = snapshot.ny - 1 - j
            for i in range(snapshot.nx):
                top = snapshot.field(field_name, i, j)
                bottom = snapshot.field(field_name, i, mirror_j)
                differences.append(top - bottom)
                scales.append(0.5 * (abs(top) + abs(bottom)))
        metrics[field_name] = field_metrics(differences, scales)

    for field_name in antisymmetric_fields:
        differences = []
        scales = []
        for j in range(pair_count):
            mirror_j = snapshot.ny - 1 - j
            for i in range(snapshot.nx):
                top = snapshot.field(field_name, i, j)
                bottom = snapshot.field(field_name, i, mirror_j)
                differences.append(top + bottom)
                scales.append(0.5 * (abs(top) + abs(bottom)))
        metrics[field_name] = field_metrics(differences, scales)

    total_y_momentum = 0.0
    for density, velocity_y in zip(snapshot.fields["density"], snapshot.fields["velocity_y"]):
        total_y_momentum += density * velocity_y
    metrics["global"] = {"total_y_momentum": total_y_momentum}

    return metrics


def y_invariant_metrics(snapshot: Snapshot) -> dict[str, dict[str, float]]:
    metrics: dict[str, dict[str, float]] = {}
    reference_row_indices = range(snapshot.nx)

    for field_name, values in snapshot.fields.items():
        reference_row = [values[index] for index in reference_row_indices]
        differences: list[float] = []
        scales: list[float] = []

        for j in range(snapshot.ny):
            row_start = j * snapshot.nx
            for i in range(snapshot.nx):
                current = values[row_start + i]
                reference = reference_row[i]
                differences.append(current - reference)
                scales.append(abs(reference))

        metrics[field_name] = field_metrics(differences, scales)

    return metrics


def compare_snapshots(candidate: Snapshot, reference: Snapshot) -> dict[str, dict[str, float]]:
    if candidate.nx != reference.nx or candidate.ny != reference.ny:
        raise ValueError("Candidate and baseline snapshots must use the same grid dimensions.")
    for candidate_x, reference_x in zip(candidate.x_values, reference.x_values):
        if abs(candidate_x - reference_x) > 1.0e-12 * max(1.0, abs(candidate_x), abs(reference_x)):
            raise ValueError("Candidate and baseline snapshots do not share the same x coordinates.")
    for candidate_y, reference_y in zip(candidate.y_values, reference.y_values):
        if abs(candidate_y - reference_y) > 1.0e-12 * max(1.0, abs(candidate_y), abs(reference_y)):
            raise ValueError("Candidate and baseline snapshots do not share the same y coordinates.")

    comparison: dict[str, dict[str, float]] = {}
    for field_name in sorted(candidate.fields):
        if field_name not in reference.fields:
            continue
        comparison[field_name] = relative_error_metrics(candidate.fields[field_name], reference.fields[field_name])
    return comparison


def print_metrics(title: str, metrics: dict[str, dict[str, float]]) -> None:
    print(title)
    for field_name, values in metrics.items():
        if field_name == "global":
            print(f"  {field_name}: total_y_momentum = {values['total_y_momentum']:.6e}")
            continue
        print(
            f"  {field_name}: L1 = {values['l1']:.6e}, "
            f"L2 = {values['l2']:.6e}, Linf = {values['linf']:.6e}, "
            f"max abs = {values['max_abs']:.6e}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check positivity, symmetry, and regression metrics for solver CSV snapshots."
    )
    parser.add_argument("snapshot", help="CSV snapshot produced by shock_bubble_solver.")
    parser.add_argument("--baseline", help="Reference snapshot for regression checks.")
    parser.add_argument(
        "--expect-y-invariant",
        action="store_true",
        help="Check that the solution stays identical in every y row. Use this with the planar-shock case.",
    )
    parser.add_argument(
        "--y-invariant-tolerance",
        type=float,
        default=None,
        help="Optional Linf tolerance for density and pressure in the y-invariance check.",
    )
    parser.add_argument(
        "--symmetry-tolerance",
        type=float,
        default=None,
        help="Optional Linf tolerance for density and pressure in the symmetry check.",
    )
    parser.add_argument(
        "--regression-tolerance",
        type=float,
        default=None,
        help="Optional Linf tolerance for density and pressure when comparing against --baseline.",
    )
    args = parser.parse_args()

    snapshot = load_snapshot(args.snapshot)

    print(f"Snapshot: {snapshot.path}")
    print(f"Grid: {snapshot.nx} x {snapshot.ny}")
    if "time" in snapshot.metadata:
        print(f"Time: {snapshot.metadata['time']}")
    if "step" in snapshot.metadata:
        print(f"Step: {snapshot.metadata['step']}")
    if "case" in snapshot.metadata:
        print(f"Case: {snapshot.metadata['case']}")
    print()

    issues, summary = summarise_snapshot(snapshot)
    print("Basic checks")
    for line in summary:
        print(f"  {line}")

    failed = False
    if issues:
        failed = True
        print("  Issues:")
        for issue in issues:
            print(f"    {issue}")
    print()

    symmetry = symmetry_metrics(snapshot)
    print_metrics("Symmetry check", symmetry)
    print()

    if args.symmetry_tolerance is not None:
        for field_name in ("density", "pressure"):
            if symmetry[field_name]["linf"] > args.symmetry_tolerance:
                failed = True
                print(
                    f"Symmetry tolerance failed for {field_name}: "
                    f"{symmetry[field_name]['linf']:.6e} > {args.symmetry_tolerance:.6e}"
                )
        if failed:
            print()

    if args.expect_y_invariant:
        invariance = y_invariant_metrics(snapshot)
        print_metrics("Y-invariance check", invariance)
        print()

        if args.y_invariant_tolerance is not None:
            for field_name in ("density", "pressure", "velocity_x", "velocity_y"):
                if invariance[field_name]["linf"] > args.y_invariant_tolerance:
                    failed = True
                    print(
                        f"Y-invariance tolerance failed for {field_name}: "
                        f"{invariance[field_name]['linf']:.6e} > {args.y_invariant_tolerance:.6e}"
                    )
            if failed:
                print()

    if args.baseline:
        baseline = load_snapshot(args.baseline)
        comparison = compare_snapshots(snapshot, baseline)
        print_metrics("Baseline comparison", comparison)
        print()

        if args.regression_tolerance is not None:
            for field_name in ("density", "pressure"):
                if comparison[field_name]["linf"] > args.regression_tolerance:
                    failed = True
                    print(
                        f"Regression tolerance failed for {field_name}: "
                        f"{comparison[field_name]['linf']:.6e} > {args.regression_tolerance:.6e}"
                    )
            if failed:
                print()

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
