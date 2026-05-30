#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from snapshot_tools import load_snapshot


def post_shock_state(gamma: float, ambient_pressure: float, air_density: float, mach: float) -> dict[str, float]:
    sound_speed = math.sqrt(gamma * ambient_pressure / air_density)
    shock_speed = mach * sound_speed
    mach_squared = mach * mach
    density_ratio = ((gamma + 1.0) * mach_squared) / ((gamma - 1.0) * mach_squared + 2.0)
    pressure_ratio = 1.0 + (2.0 * gamma / (gamma + 1.0)) * (mach_squared - 1.0)
    return {
        "density": air_density * density_ratio,
        "pressure": ambient_pressure * pressure_ratio,
        "velocity_x": shock_speed * (1.0 - 1.0 / density_ratio),
        "velocity_y": 0.0,
        "shock_speed": shock_speed,
    }


def field_metrics(values: list[float], reference: float) -> dict[str, float | str]:
    differences = [value - reference for value in values]
    abs_differences = [abs(value) for value in differences]
    l1_abs = sum(abs_differences) / len(abs_differences)
    linf_abs = max(abs_differences)

    reference_l1_scale = abs(reference)
    reference_linf_scale = abs(reference)
    l1_rel: float | str
    linf_rel: float | str
    if reference_l1_scale > 1.0e-30:
        l1_rel = l1_abs / reference_l1_scale
        linf_rel = linf_abs / reference_linf_scale
    else:
        l1_rel = ""
        linf_rel = ""

    return {
        "l1_abs": l1_abs,
        "linf_abs": linf_abs,
        "l1_rel": l1_rel,
        "linf_rel": linf_rel,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare planar-shock plateau regions against the exact Rankine-Hugoniot states."
    )
    parser.add_argument("snapshot", help="Planar-shock CSV snapshot produced by shock_bubble_solver.")
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--ambient-pressure", type=float, default=1.01325e5)
    parser.add_argument("--air-density", type=float, default=1.29)
    parser.add_argument("--mach", type=float, default=1.22)
    parser.add_argument("--shock-location", type=float, default=0.005)
    parser.add_argument(
        "--final-time",
        type=float,
        default=None,
        help="Override the final time in seconds. If omitted, use the snapshot metadata time.",
    )
    parser.add_argument(
        "--transition-cells",
        type=int,
        default=6,
        help="Exclude this many cells on each side of the exact shock location when defining plateau masks.",
    )
    parser.add_argument("--csv-output", type=Path, default=None)
    parser.add_argument("--txt-output", type=Path, default=None)
    args = parser.parse_args()

    snapshot = load_snapshot(args.snapshot)
    if args.final_time is None:
        if "time" not in snapshot.metadata:
            raise ValueError("Snapshot metadata do not contain time; pass --final-time explicitly.")
        final_time = float(snapshot.metadata["time"])
    else:
        final_time = args.final_time

    exact_post = post_shock_state(args.gamma, args.ambient_pressure, args.air_density, args.mach)
    exact_pre = {
        "density": args.air_density,
        "pressure": args.ambient_pressure,
        "velocity_x": 0.0,
        "velocity_y": 0.0,
    }

    if snapshot.nx < 2:
        raise ValueError("Need at least two x cells to infer grid spacing.")
    dx = snapshot.x_values[1] - snapshot.x_values[0]
    shock_position = args.shock_location + exact_post["shock_speed"] * final_time
    transition_half_width = args.transition_cells * dx

    left_indices: list[int] = []
    right_indices: list[int] = []
    for j in range(snapshot.ny):
        for i, x in enumerate(snapshot.x_values):
            flat_index = j * snapshot.nx + i
            if x <= shock_position - transition_half_width:
                left_indices.append(flat_index)
            elif x >= shock_position + transition_half_width:
                right_indices.append(flat_index)

    if not left_indices or not right_indices:
        raise ValueError("Plateau masks are empty; reduce --transition-cells or change the verification time.")

    rows: list[dict[str, float | int | str]] = []
    for region_name, indices, exact_state in (
        ("shocked", left_indices, exact_post),
        ("ambient", right_indices, exact_pre),
    ):
        for field_name in ("density", "pressure", "velocity_x", "velocity_y"):
            field_values = [snapshot.fields[field_name][index] for index in indices]
            metrics = field_metrics(field_values, float(exact_state[field_name]))
            rows.append(
                {
                    "region": region_name,
                    "field": field_name,
                    "reference_value": float(exact_state[field_name]),
                    "cell_count": len(indices),
                    "shock_position_mm": 1000.0 * shock_position,
                    "transition_half_width_mm": 1000.0 * transition_half_width,
                    "l1_abs": metrics["l1_abs"],
                    "linf_abs": metrics["linf_abs"],
                    "l1_rel": metrics["l1_rel"],
                    "linf_rel": metrics["linf_rel"],
                }
            )

    if args.csv_output is not None:
        args.csv_output.parent.mkdir(parents=True, exist_ok=True)
        with args.csv_output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    shocked_density = next(row for row in rows if row["region"] == "shocked" and row["field"] == "density")
    shocked_pressure = next(row for row in rows if row["region"] == "shocked" and row["field"] == "pressure")
    ambient_density = next(row for row in rows if row["region"] == "ambient" and row["field"] == "density")
    ambient_pressure = next(row for row in rows if row["region"] == "ambient" and row["field"] == "pressure")
    ambient_velocity_x = next(row for row in rows if row["region"] == "ambient" and row["field"] == "velocity_x")

    summary_lines = [
        "Planar-shock reference check",
        f"Snapshot: {snapshot.path}",
        f"Grid: {snapshot.nx} x {snapshot.ny}",
        f"Final time: {final_time:.12e} s",
        f"Exact shock position: {1000.0 * shock_position:.3f} mm",
        f"Excluded transition half-width: {1000.0 * transition_half_width:.3f} mm ({args.transition_cells} cells)",
        f"Shocked plateau cell count: {len(left_indices)}",
        f"Ambient plateau cell count: {len(right_indices)}",
        (
            "Density plateau relative Linf: "
            f"shocked = {100.0 * float(shocked_density['linf_rel']):.3f}%, "
            f"ambient = {100.0 * float(ambient_density['linf_rel']):.3f}%"
        ),
        (
            "Pressure plateau relative Linf: "
            f"shocked = {100.0 * float(shocked_pressure['linf_rel']):.3f}%, "
            f"ambient = {100.0 * float(ambient_pressure['linf_rel']):.3f}%"
        ),
        f"Ambient velocity_x absolute Linf: {float(ambient_velocity_x['linf_abs']):.6e} m/s",
    ]

    if args.txt_output is not None:
        args.txt_output.parent.mkdir(parents=True, exist_ok=True)
        args.txt_output.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    for line in summary_lines:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
