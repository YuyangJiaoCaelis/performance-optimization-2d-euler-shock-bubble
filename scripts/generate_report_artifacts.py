#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from snapshot_tools import load_snapshot


GAMMA = 1.4
AMBIENT_DENSITY = 1.29
HELIUM_DENSITY = 0.214
AMBIENT_PRESSURE = 1.01325e5
MACH = 1.22
BUBBLE_RADIUS = 0.025
INTERFACE_THRESHOLD = 0.5 * (AMBIENT_DENSITY + HELIUM_DENSITY)
INTERFACE_THRESHOLD_FRACTIONS = (0.45, 0.50, 0.55)
SOUND_SPEED = math.sqrt(GAMMA * AMBIENT_PRESSURE / AMBIENT_DENSITY)
TIME_SCALE = BUBBLE_RADIUS / (SOUND_SPEED * MACH)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPORTS_DIR = PROJECT_ROOT / "reports"
FIGURES_DIR = PROJECT_ROOT / "figures"
INPUTS_DIR = PROJECT_ROOT / "inputs"
PAPER_DIGITIZED_CSV = INPUTS_DIR / "paper_digitized_interfaces.csv"
PLANAR_REFERENCE_CHECK_CSV = REPORTS_DIR / "planar_shock_reference_check.csv"

PRODUCTION_SNAPSHOT = PROJECT_ROOT / "output" / "default_run_0001.csv"
BENCHMARK_SNAPSHOTS = sorted((PROJECT_ROOT / "output" / "benchmark_series").glob("baseline_*.csv"))
GRID_SNAPSHOTS = {
    "250x99": PROJECT_ROOT / "output" / "grid_sensitivity" / "nx250_0001.csv",
    "500x197": PRODUCTION_SNAPSHOT,
    "1000x394": PROJECT_ROOT / "output" / "grid_sensitivity" / "nx1000_0001.csv",
}

PAPER_REFERENCE_TSTAR = 3.0

GRID_SENSITIVITY_FIGURE = "Figure_1_grid_sensitivity_density.png"
PRODUCTION_FIELDS_FIGURE = "Figure_2_production_fields.png"
INTERFACE_TRACKING_FIGURE = "Figure_3_interface_tracking.png"


def load_paper_digitized_points() -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with PAPER_DIGITIZED_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                {
                    "source": row["source"],
                    "interface": row["interface"],
                    "time_star": float(row["time_star"]),
                    "position_mm": float(row["position_mm"]),
                }
            )
    if not rows:
        raise ValueError(f"No paper digitization points found in {PAPER_DIGITIZED_CSV}.")
    return rows


PAPER_DIGITIZED_POINTS = load_paper_digitized_points()


def centreline_index(snapshot) -> int:
    return min(range(snapshot.ny), key=lambda j: abs(snapshot.y_values[j]))


def structured_field(snapshot, field_name: str) -> np.ndarray:
    return np.asarray(snapshot.fields[field_name], dtype=float).reshape(snapshot.ny, snapshot.nx)


def centreline_values(snapshot, field_name: str) -> np.ndarray:
    data = structured_field(snapshot, field_name)
    return data[centreline_index(snapshot), :]


def threshold_crossings(snapshot, threshold: float = INTERFACE_THRESHOLD) -> list[float]:
    values = centreline_values(snapshot, "density")
    crossings: list[float] = []
    for i in range(snapshot.nx - 1):
        left = values[i]
        right = values[i + 1]
        if left == right:
            continue
        if (left - threshold) * (right - threshold) > 0.0:
            continue
        fraction = (threshold - left) / (right - left)
        x_position = snapshot.x_values[i] + fraction * (snapshot.x_values[i + 1] - snapshot.x_values[i])
        crossings.append(x_position)
    return crossings


def threshold_from_fraction(fraction: float) -> float:
    return HELIUM_DENSITY + fraction * (AMBIENT_DENSITY - HELIUM_DENSITY)


def nondimensional_time(physical_time: float) -> float:
    return physical_time / TIME_SCALE


def interface_series_for_threshold(threshold: float) -> tuple[list[dict[str, float]], float]:
    rows: list[dict[str, float]] = []
    initial_snapshot = load_snapshot(BENCHMARK_SNAPSHOTS[0])
    initial_crossings = threshold_crossings(initial_snapshot, threshold=threshold)
    initial_left = initial_crossings[0]

    for path in BENCHMARK_SNAPSHOTS:
        snapshot = load_snapshot(path)
        crossings = threshold_crossings(snapshot, threshold=threshold)
        if len(crossings) != 2:
            raise ValueError(f"Expected two centreline crossings in {path}, found {len(crossings)}.")

        time_s = float(snapshot.metadata["time"])
        rows.append(
            {
                "time_s": time_s,
                "time_star": nondimensional_time(time_s),
                "left_crossing_mm": 1000.0 * crossings[0],
                "right_crossing_mm": 1000.0 * crossings[1],
                "left_from_initial_left_mm": 1000.0 * (crossings[0] - initial_left),
                "right_from_initial_left_mm": 1000.0 * (crossings[1] - initial_left),
            }
        )

    return rows, initial_left


def save_interface_series() -> tuple[list[dict[str, float]], float]:
    rows, initial_left = interface_series_for_threshold(INTERFACE_THRESHOLD)

    output_path = REPORTS_DIR / "interface_proxy_series.csv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    return rows, initial_left


def save_interface_threshold_sensitivity() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []

    fig3_points = [row for row in PAPER_DIGITIZED_POINTS if row["source"] == "Fig. 3"]
    production_snapshot = load_snapshot(GRID_SNAPSHOTS["500x197"])
    fine_snapshot = load_snapshot(GRID_SNAPSHOTS["1000x394"])

    for fraction in INTERFACE_THRESHOLD_FRACTIONS:
        threshold = threshold_from_fraction(fraction)
        interface_rows, _ = interface_series_for_threshold(threshold)
        t3_row = min(interface_rows, key=lambda row: abs(row["time_star"] - PAPER_REFERENCE_TSTAR))

        jet_differences: list[float] = []
        downstream_differences: list[float] = []
        for point in fig3_points:
            matched_row = min(interface_rows, key=lambda row: abs(row["time_star"] - point["time_star"]))
            if point["interface"] == "jet":
                jet_differences.append(abs(matched_row["left_from_initial_left_mm"] - point["position_mm"]))
            elif point["interface"] == "downstream":
                downstream_differences.append(abs(matched_row["right_from_initial_left_mm"] - point["position_mm"]))

        production_crossings = [1000.0 * value for value in threshold_crossings(production_snapshot, threshold=threshold)]
        fine_crossings = [1000.0 * value for value in threshold_crossings(fine_snapshot, threshold=threshold)]

        rows.append(
            {
                "fraction_between_helium_and_air": fraction,
                "threshold_kg_m3": threshold,
                "t3_left_from_initial_left_mm": t3_row["left_from_initial_left_mm"],
                "t3_right_from_initial_left_mm": t3_row["right_from_initial_left_mm"],
                "fig3_jet_mean_abs_difference_mm": sum(jet_differences) / len(jet_differences),
                "fig3_downstream_mean_abs_difference_mm": sum(downstream_differences) / len(downstream_differences),
                "production_fine_left_shift_mm": abs(fine_crossings[0] - production_crossings[0]),
                "production_fine_right_shift_mm": abs(fine_crossings[1] - production_crossings[1]),
            }
        )

    with (REPORTS_DIR / "interface_threshold_sensitivity.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    central_row = next(row for row in rows if abs(row["fraction_between_helium_and_air"] - 0.50) < 1.0e-12)
    summary_lines = [
        "Interface-threshold sensitivity summary",
        (
            f"Threshold sweep fractions: {', '.join(f'{fraction:.2f}' for fraction in INTERFACE_THRESHOLD_FRACTIONS)} "
            "between helium and air density."
        ),
        (
            f"t* ~= 3 jet proxy range: {min(row['t3_left_from_initial_left_mm'] for row in rows):.3f} to "
            f"{max(row['t3_left_from_initial_left_mm'] for row in rows):.3f} mm "
            f"(central {central_row['t3_left_from_initial_left_mm']:.3f} mm)"
        ),
        (
            f"t* ~= 3 downstream proxy range: {min(row['t3_right_from_initial_left_mm'] for row in rows):.3f} to "
            f"{max(row['t3_right_from_initial_left_mm'] for row in rows):.3f} mm "
            f"(central {central_row['t3_right_from_initial_left_mm']:.3f} mm)"
        ),
        (
            f"Late-time Fig. 3 jet MAE range: {min(row['fig3_jet_mean_abs_difference_mm'] for row in rows):.3f} to "
            f"{max(row['fig3_jet_mean_abs_difference_mm'] for row in rows):.3f} mm"
        ),
        (
            f"Late-time Fig. 3 downstream MAE range: "
            f"{min(row['fig3_downstream_mean_abs_difference_mm'] for row in rows):.3f} to "
            f"{max(row['fig3_downstream_mean_abs_difference_mm'] for row in rows):.3f} mm"
        ),
        (
            f"Production-fine left-shift range: {min(row['production_fine_left_shift_mm'] for row in rows):.3f} to "
            f"{max(row['production_fine_left_shift_mm'] for row in rows):.3f} mm"
        ),
        (
            f"Production-fine right-shift range: {min(row['production_fine_right_shift_mm'] for row in rows):.3f} to "
            f"{max(row['production_fine_right_shift_mm'] for row in rows):.3f} mm"
        ),
    ]
    (REPORTS_DIR / "interface_threshold_sensitivity.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    return rows


def save_grid_metrics() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for grid_label, path in GRID_SNAPSHOTS.items():
        snapshot = load_snapshot(path)
        crossings = threshold_crossings(snapshot)
        density = centreline_values(snapshot, "density")
        velocity_x = centreline_values(snapshot, "velocity_x")

        min_index = int(np.argmin(density))
        max_velocity_index = int(np.argmax(velocity_x))

        rows.append(
            {
                "grid": grid_label,
                "left_crossing_mm": 1000.0 * crossings[0],
                "right_crossing_mm": 1000.0 * crossings[1],
                "min_density": float(density[min_index]),
                "min_density_x_mm": 1000.0 * snapshot.x_values[min_index],
                "max_velocity_x": float(velocity_x[max_velocity_index]),
                "max_velocity_x_x_mm": 1000.0 * snapshot.x_values[max_velocity_index],
            }
        )

    output_path = REPORTS_DIR / "grid_sensitivity_metrics.csv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    return rows


def interpolate_field_to_target(snapshot, field_name: str, target_x: np.ndarray, target_y: np.ndarray) -> np.ndarray:
    source_x = np.asarray(snapshot.x_values, dtype=float)
    source_y = np.asarray(snapshot.y_values, dtype=float)
    data = structured_field(snapshot, field_name)

    interpolated_x = np.vstack([np.interp(target_x, source_x, row) for row in data])
    interpolated_xy = np.vstack(
        [np.interp(target_y, source_y, interpolated_x[:, column]) for column in range(interpolated_x.shape[1])]
    ).T
    return interpolated_xy


def relative_metrics(candidate: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    difference = candidate - reference
    abs_reference = np.abs(reference)
    return {
        "l1": float(np.abs(difference).sum() / max(abs_reference.sum(), 1.0e-30)),
        "l2": float(np.sqrt((difference * difference).sum()) / max(np.sqrt((reference * reference).sum()), 1.0e-30)),
        "linf": float(np.abs(difference).max() / max(abs_reference.max(), 1.0e-30)),
        "abs_l1": float(np.abs(difference).sum() / difference.size),
        "abs_linf": float(np.abs(difference).max()),
    }


def save_convergence_metrics() -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    coarse = load_snapshot(GRID_SNAPSHOTS["250x99"])
    medium = load_snapshot(GRID_SNAPSHOTS["500x197"])
    fine = load_snapshot(GRID_SNAPSHOTS["1000x394"])

    target_x = np.asarray(coarse.x_values, dtype=float)
    target_y = np.asarray(coarse.y_values, dtype=float)

    field_rows: list[dict[str, float]] = []
    for field_name in ("density", "pressure"):
        coarse_field = structured_field(coarse, field_name)
        medium_on_coarse = interpolate_field_to_target(medium, field_name, target_x, target_y)
        fine_on_coarse = interpolate_field_to_target(fine, field_name, target_x, target_y)

        e_h = relative_metrics(coarse_field, medium_on_coarse)
        e_h2 = relative_metrics(medium_on_coarse, fine_on_coarse)

        field_rows.append(
            {
                "field": field_name,
                "eh_l1": e_h["l1"],
                "eh_l2": e_h["l2"],
                "eh_linf": e_h["linf"],
                "eh2_l1": e_h2["l1"],
                "eh2_l2": e_h2["l2"],
                "eh2_linf": e_h2["linf"],
                "order_l1": math.log(e_h["l1"] / e_h2["l1"], 2.0),
                "order_l2": math.log(e_h["l2"] / e_h2["l2"], 2.0),
                "order_linf": math.log(e_h["linf"] / e_h2["linf"], 2.0),
            }
        )

    interface_rows: list[dict[str, float]] = []
    coarse_crossings = [1000.0 * value for value in threshold_crossings(coarse)]
    medium_crossings = [1000.0 * value for value in threshold_crossings(medium)]
    fine_crossings = [1000.0 * value for value in threshold_crossings(fine)]

    for index, name in enumerate(("left", "right")):
        e_h = abs(coarse_crossings[index] - medium_crossings[index])
        e_h2 = abs(medium_crossings[index] - fine_crossings[index])
        interface_rows.append(
            {
                "interface": name,
                "coarse_mm": coarse_crossings[index],
                "medium_mm": medium_crossings[index],
                "fine_mm": fine_crossings[index],
                "eh_mm": e_h,
                "eh2_mm": e_h2,
                "order": math.log(e_h / e_h2, 2.0),
            }
        )

    with (REPORTS_DIR / "convergence_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(field_rows[0].keys()))
        writer.writeheader()
        writer.writerows(field_rows)

    with (REPORTS_DIR / "interface_convergence.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(interface_rows[0].keys()))
        writer.writeheader()
        writer.writerows(interface_rows)

    summary_lines = ["Three-level convergence summary"]
    for row in field_rows:
        summary_lines.append(
            (
                f"{row['field']}: L1 errors {row['eh_l1']:.6e} -> {row['eh2_l1']:.6e}, "
                f"observed order {row['order_l1']:.3f}; "
                f"L2 errors {row['eh_l2']:.6e} -> {row['eh2_l2']:.6e}, "
                f"observed order {row['order_l2']:.3f}"
            )
        )
    for row in interface_rows:
        summary_lines.append(
            (
                f"{row['interface']} interface: {row['coarse_mm']:.3f} -> {row['medium_mm']:.3f} -> {row['fine_mm']:.3f} mm, "
                f"observed order {row['order']:.3f}"
            )
        )
    (REPORTS_DIR / "convergence_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    return field_rows, interface_rows


def save_interface_richardson_summary(interface_convergence_rows: list[dict[str, float]]) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    refinement_ratio = 2.0
    safety_factor = 1.25

    for row in interface_convergence_rows:
        coarse = float(row["coarse_mm"])
        medium = float(row["medium_mm"])
        fine = float(row["fine_mm"])
        order = float(row["order"])

        extrapolated = fine + (fine - medium) / (refinement_ratio**order - 1.0)
        relative_error_medium = abs(extrapolated - medium) / extrapolated
        relative_error_fine = abs(extrapolated - fine) / extrapolated
        gci21 = safety_factor * abs((fine - medium) / fine) / (refinement_ratio**order - 1.0)
        gci32 = safety_factor * abs((medium - coarse) / medium) / (refinement_ratio**order - 1.0)
        asymptotic_ratio = gci32 / (refinement_ratio**order * gci21)

        rows.append(
            {
                "interface": row["interface"],
                "order": order,
                "extrapolated_mm": extrapolated,
                "production_relative_error": relative_error_medium,
                "fine_relative_error": relative_error_fine,
                "gci21": gci21,
                "gci32": gci32,
                "asymptotic_ratio": asymptotic_ratio,
            }
        )

    output_path = REPORTS_DIR / "interface_richardson_summary.csv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary_lines = ["Interface Richardson/GCI summary"]
    for row in rows:
        summary_lines.append(
            (
                f"{row['interface']}: observed order = {row['order']:.3f}, "
                f"extrapolated position = {row['extrapolated_mm']:.3f} mm, "
                f"production relative error = {100.0 * row['production_relative_error']:.2f}%, "
                f"fine relative error = {100.0 * row['fine_relative_error']:.2f}%, "
                f"asymptotic ratio = {row['asymptotic_ratio']:.3f}"
            )
        )
    (REPORTS_DIR / "interface_richardson_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    return rows


def save_field_richardson_summary(convergence_rows: list[dict[str, float]]) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    refinement_ratio = 2.0
    safety_factor = 1.25

    for row in convergence_rows:
        for norm in ("l1", "l2", "linf"):
            e32 = float(row[f"eh_{norm}"])
            e21 = float(row[f"eh2_{norm}"])
            order = float(row[f"order_{norm}"])
            denominator = refinement_ratio**order - 1.0
            gci21 = safety_factor * e21 / denominator
            gci32 = safety_factor * e32 / denominator

            rows.append(
                {
                    "field": row["field"],
                    "norm": norm.upper(),
                    "order": order,
                    "gci21": gci21,
                    "gci32": gci32,
                }
            )

    with (REPORTS_DIR / "field_richardson_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary_lines = ["Field-norm Richardson/GCI summary"]
    for field_name in ("density", "pressure"):
        subset = [row for row in rows if row["field"] == field_name and row["norm"] in ("L1", "L2")]
        summary_lines.append(
            field_name
            + ": "
            + "; ".join(
                (
                    f"{row['norm']}: order = {row['order']:.3f}, "
                    f"GCI21 = {100.0 * float(row['gci21']):.2f}%, "
                    f"GCI32 = {100.0 * float(row['gci32']):.2f}%"
                )
                for row in subset
            )
        )
    (REPORTS_DIR / "field_richardson_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    return rows


def save_paper_digitization() -> list[dict[str, float | str]]:
    output_path = REPORTS_DIR / "paper_digitized_interfaces.csv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source", "interface", "time_star", "position_mm"])
        writer.writeheader()
        writer.writerows(PAPER_DIGITIZED_POINTS)
    return PAPER_DIGITIZED_POINTS


def save_quantitative_summary(
    interface_rows: list[dict[str, float]],
    grid_rows: list[dict[str, float]],
    convergence_rows: list[dict[str, float]],
    interface_convergence_rows: list[dict[str, float]],
    interface_richardson_rows: list[dict[str, float]],
    field_richardson_rows: list[dict[str, float | str]],
    threshold_sensitivity_rows: list[dict[str, float]],
) -> None:
    t3_row = min(interface_rows, key=lambda row: abs(row["time_star"] - PAPER_REFERENCE_TSTAR))
    medium_row = next(row for row in grid_rows if row["grid"] == "500x197")
    fine_row = next(row for row in grid_rows if row["grid"] == "1000x394")

    medium_fine_left = abs(fine_row["left_crossing_mm"] - medium_row["left_crossing_mm"])
    medium_fine_right = abs(fine_row["right_crossing_mm"] - medium_row["right_crossing_mm"])

    paper_t3_points = [
        row for row in PAPER_DIGITIZED_POINTS if abs(row["time_star"] - PAPER_REFERENCE_TSTAR) < 0.05
    ]
    paper_jet_positions = [row["position_mm"] for row in paper_t3_points if row["interface"] == "jet"]
    paper_downstream_positions = [row["position_mm"] for row in paper_t3_points if row["interface"] == "downstream"]

    jet_difference_min = min(abs(t3_row["left_from_initial_left_mm"] - value) for value in paper_jet_positions)
    jet_difference_max = max(abs(t3_row["left_from_initial_left_mm"] - value) for value in paper_jet_positions)
    downstream_difference_min = min(
        abs(t3_row["right_from_initial_left_mm"] - value) for value in paper_downstream_positions
    )
    downstream_difference_max = max(
        abs(t3_row["right_from_initial_left_mm"] - value) for value in paper_downstream_positions
    )

    fig3_points = [row for row in PAPER_DIGITIZED_POINTS if row["source"] == "Fig. 3"]
    fig3_jet_differences: list[float] = []
    fig3_downstream_differences: list[float] = []
    for point in fig3_points:
        matched_row = min(interface_rows, key=lambda row: abs(row["time_star"] - point["time_star"]))
        if point["interface"] == "jet":
            fig3_jet_differences.append(abs(matched_row["left_from_initial_left_mm"] - point["position_mm"]))
        elif point["interface"] == "downstream":
            fig3_downstream_differences.append(abs(matched_row["right_from_initial_left_mm"] - point["position_mm"]))

    density_convergence = next(row for row in convergence_rows if row["field"] == "density")
    density_l1_richardson = next(
        row for row in field_richardson_rows if row["field"] == "density" and row["norm"] == "L1"
    )
    density_l2_richardson = next(
        row for row in field_richardson_rows if row["field"] == "density" and row["norm"] == "L2"
    )
    pressure_l1_richardson = next(
        row for row in field_richardson_rows if row["field"] == "pressure" and row["norm"] == "L1"
    )
    pressure_l2_richardson = next(
        row for row in field_richardson_rows if row["field"] == "pressure" and row["norm"] == "L2"
    )
    left_richardson = next(row for row in interface_richardson_rows if row["interface"] == "left")
    right_richardson = next(row for row in interface_richardson_rows if row["interface"] == "right")
    central_threshold_row = next(
        row for row in threshold_sensitivity_rows if abs(row["fraction_between_helium_and_air"] - 0.50) < 1.0e-12
    )
    planar_summary_line = None
    if PLANAR_REFERENCE_CHECK_CSV.exists():
        with PLANAR_REFERENCE_CHECK_CSV.open("r", newline="", encoding="utf-8") as handle:
            planar_rows = list(csv.DictReader(handle))
        shocked_density = next(row for row in planar_rows if row["region"] == "shocked" and row["field"] == "density")
        shocked_pressure = next(row for row in planar_rows if row["region"] == "shocked" and row["field"] == "pressure")
        ambient_density = next(row for row in planar_rows if row["region"] == "ambient" and row["field"] == "density")
        ambient_pressure = next(row for row in planar_rows if row["region"] == "ambient" and row["field"] == "pressure")
        ambient_velocity = next(row for row in planar_rows if row["region"] == "ambient" and row["field"] == "velocity_x")
        planar_summary_line = (
            "Planar-shock plateau check: "
            f"shocked density/pressure Linf = {100.0 * float(shocked_density['linf_rel']):.2f}%/"
            f"{100.0 * float(shocked_pressure['linf_rel']):.2f}%, "
            f"ambient density/pressure Linf = {100.0 * float(ambient_density['linf_rel']):.2f}%/"
            f"{100.0 * float(ambient_pressure['linf_rel']):.2f}%, "
            f"ambient velocity_x Linf = {float(ambient_velocity['linf_abs']):.3e} m/s"
        )

    summary_lines = [
        "Quantitative report evidence",
        f"Interface threshold: {INTERFACE_THRESHOLD:.3f} kg/m^3",
        f"Time scale r/(a_0 M_s): {TIME_SCALE:.12e} s",
        (
            "Centreline proxy at t* ~= 3: "
            f"left-from-initial-left = {t3_row['left_from_initial_left_mm']:.3f} mm, "
            f"right-from-initial-left = {t3_row['right_from_initial_left_mm']:.3f} mm"
        ),
        (
            "Digitized paper range at t* ~= 3: "
            f"jet = {min(paper_jet_positions):.1f} to {max(paper_jet_positions):.1f} mm, "
            f"downstream = {min(paper_downstream_positions):.1f} to {max(paper_downstream_positions):.1f} mm"
        ),
        (
            "Difference to digitized paper range: "
            f"jet proxy = {jet_difference_min:.3f} to {jet_difference_max:.3f} mm, "
            f"downstream proxy = {downstream_difference_min:.3f} to {downstream_difference_max:.3f} mm"
        ),
    ]
    if planar_summary_line is not None:
        summary_lines.append(planar_summary_line)
    summary_lines.extend(
        [
        (
            "Late-time Fig. 3 trajectory comparison: "
            f"jet mean abs difference = {sum(fig3_jet_differences) / len(fig3_jet_differences):.3f} mm "
            f"(max {max(fig3_jet_differences):.3f} mm), "
            f"downstream mean abs difference = {sum(fig3_downstream_differences) / len(fig3_downstream_differences):.3f} mm "
            f"(max {max(fig3_downstream_differences):.3f} mm)"
        ),
        (
            "Medium-fine interface shift at final time: "
            f"left = {medium_fine_left:.3f} mm, right = {medium_fine_right:.3f} mm"
        ),
        (
            "Threshold sensitivity (45%-55% between helium and air density): "
            f"t* ~= 3 jet proxy = {min(row['t3_left_from_initial_left_mm'] for row in threshold_sensitivity_rows):.3f} to "
            f"{max(row['t3_left_from_initial_left_mm'] for row in threshold_sensitivity_rows):.3f} mm "
            f"(central {central_threshold_row['t3_left_from_initial_left_mm']:.3f} mm), "
            f"downstream proxy = {min(row['t3_right_from_initial_left_mm'] for row in threshold_sensitivity_rows):.3f} to "
            f"{max(row['t3_right_from_initial_left_mm'] for row in threshold_sensitivity_rows):.3f} mm "
            f"(central {central_threshold_row['t3_right_from_initial_left_mm']:.3f} mm)"
        ),
        (
            "Three-level density self-convergence: "
            f"L1 order = {density_convergence['order_l1']:.3f}, "
            f"L2 order = {density_convergence['order_l2']:.3f}"
        ),
        (
            "Field-norm Richardson summary: "
            f"density GCI21 L1/L2 = {100.0 * float(density_l1_richardson['gci21']):.2f}%/"
            f"{100.0 * float(density_l2_richardson['gci21']):.2f}%; "
            f"pressure GCI21 L1/L2 = {100.0 * float(pressure_l1_richardson['gci21']):.2f}%/"
            f"{100.0 * float(pressure_l2_richardson['gci21']):.2f}%"
        ),
        (
            "Interface self-convergence orders: "
            + ", ".join(f"{row['interface']} = {row['order']:.3f}" for row in interface_convergence_rows)
        ),
        (
            "Interface Richardson summary: "
            f"left extrapolated = {left_richardson['extrapolated_mm']:.3f} mm "
            f"with production relative error {100.0 * left_richardson['production_relative_error']:.2f}%, "
            f"right extrapolated = {right_richardson['extrapolated_mm']:.3f} mm "
            f"with production relative error {100.0 * right_richardson['production_relative_error']:.2f}%"
        ),
        "Bagabir-Drikakis interface-location uncertainty reported around Fig. 3: +/-0.89 to +/-1.335 mm.",
        ]
    )

    output_path = REPORTS_DIR / "quantitative_report_evidence.txt"
    output_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


def make_production_fields_figure() -> None:
    snapshot = load_snapshot(PRODUCTION_SNAPSHOT)
    extent = [
        1000.0 * snapshot.x_values[0],
        1000.0 * snapshot.x_values[-1],
        1000.0 * snapshot.y_values[0],
        1000.0 * snapshot.y_values[-1],
    ]
    density = structured_field(snapshot, "density")
    pressure = structured_field(snapshot, "pressure") / 1.0e5

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    images = [
        axes[0].imshow(density, origin="lower", extent=extent, aspect="auto", cmap="cividis"),
        axes[1].imshow(pressure, origin="lower", extent=extent, aspect="auto", cmap="magma"),
    ]
    labels = [r"Density (kg m$^{-3}$)", "Pressure (bar)"]
    titles = ["Density at final time", "Pressure at final time"]

    for axis, image, label, title in zip(axes, images, labels, titles):
        axis.set_xlabel("x (mm)")
        axis.set_ylabel("y (mm)")
        axis.set_title(title)
        colourbar = fig.colorbar(image, ax=axis, shrink=0.9)
        colourbar.set_label(label)

    fig.savefig(FIGURES_DIR / PRODUCTION_FIELDS_FIGURE, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_interface_tracking_figure(interface_rows: list[dict[str, float]]) -> None:
    time_star = np.asarray([row["time_star"] for row in interface_rows], dtype=float)
    left_positions = np.asarray([row["left_from_initial_left_mm"] for row in interface_rows], dtype=float)
    right_positions = np.asarray([row["right_from_initial_left_mm"] for row in interface_rows], dtype=float)

    fig, ax = plt.subplots(figsize=(6.6, 4.3), constrained_layout=True)
    ax.plot(time_star, left_positions, marker="o", linewidth=1.8, markersize=4.5, label="Left centreline threshold")
    ax.plot(
        time_star,
        right_positions,
        marker="s",
        linewidth=1.8,
        markersize=4.5,
        label="Right centreline threshold",
    )
    fig3_points = [row for row in PAPER_DIGITIZED_POINTS if row["source"] == "Fig. 3"]
    fig12_points = [row for row in PAPER_DIGITIZED_POINTS if row["source"] == "Fig. 12b"]
    ax.scatter(
        [row["time_star"] for row in fig3_points],
        [row["position_mm"] for row in fig3_points],
        marker="D",
        s=40,
        color="#6a3d9a",
        label="Digitized paper points (Fig. 3)",
        zorder=4,
    )
    ax.scatter(
        [row["time_star"] for row in fig12_points],
        [row["position_mm"] for row in fig12_points],
        marker="X",
        s=48,
        color="#b15928",
        label="Digitized paper points (Fig. 12b)",
        zorder=4,
    )
    ax.axvline(PAPER_REFERENCE_TSTAR, linestyle="--", linewidth=1.0, color="0.45")
    ax.set_xlabel(r"Nondimensional time $t = T/(r/(a_0 M_s))$")
    ax.set_ylabel("Position from initial left interface (mm)")
    ax.set_title("Centreline interface proxies during shock-bubble interaction")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES_DIR / INTERFACE_TRACKING_FIGURE, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_grid_sensitivity_figure() -> None:
    fig, ax = plt.subplots(figsize=(6.6, 4.3), constrained_layout=True)
    styles = {
        "250x99": ("#d95f02", "--"),
        "500x197": ("#1b9e77", "-"),
        "1000x394": ("#7570b3", "-."),
    }

    for grid_label, path in GRID_SNAPSHOTS.items():
        snapshot = load_snapshot(path)
        x_values = 1000.0 * np.asarray(snapshot.x_values, dtype=float)
        density = centreline_values(snapshot, "density")
        colour, linestyle = styles[grid_label]
        ax.plot(x_values, density, color=colour, linestyle=linestyle, linewidth=1.8, label=grid_label)

    ax.axhline(INTERFACE_THRESHOLD, color="0.35", linestyle=":", linewidth=1.0, label=r"Interface threshold")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel(r"Centreline density (kg m$^{-3}$)")
    ax.set_title("Final-time centreline density on three grids")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES_DIR / GRID_SENSITIVITY_FIGURE, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    save_paper_digitization()
    interface_rows, _ = save_interface_series()
    threshold_sensitivity_rows = save_interface_threshold_sensitivity()
    grid_rows = save_grid_metrics()
    convergence_rows, interface_convergence_rows = save_convergence_metrics()
    interface_richardson_rows = save_interface_richardson_summary(interface_convergence_rows)
    field_richardson_rows = save_field_richardson_summary(convergence_rows)
    save_quantitative_summary(
        interface_rows,
        grid_rows,
        convergence_rows,
        interface_convergence_rows,
        interface_richardson_rows,
        field_richardson_rows,
        threshold_sensitivity_rows,
    )

    make_production_fields_figure()
    make_interface_tracking_figure(interface_rows)
    make_grid_sensitivity_figure()

    print(f"Wrote {REPORTS_DIR / 'paper_digitized_interfaces.csv'}")
    print(f"Wrote {REPORTS_DIR / 'interface_proxy_series.csv'}")
    print(f"Wrote {REPORTS_DIR / 'interface_threshold_sensitivity.csv'}")
    print(f"Wrote {REPORTS_DIR / 'interface_threshold_sensitivity.txt'}")
    print(f"Wrote {REPORTS_DIR / 'grid_sensitivity_metrics.csv'}")
    print(f"Wrote {REPORTS_DIR / 'convergence_metrics.csv'}")
    print(f"Wrote {REPORTS_DIR / 'interface_convergence.csv'}")
    print(f"Wrote {REPORTS_DIR / 'interface_richardson_summary.csv'}")
    print(f"Wrote {REPORTS_DIR / 'interface_richardson_summary.txt'}")
    print(f"Wrote {REPORTS_DIR / 'field_richardson_summary.csv'}")
    print(f"Wrote {REPORTS_DIR / 'field_richardson_summary.txt'}")
    print(f"Wrote {REPORTS_DIR / 'convergence_summary.txt'}")
    print(f"Wrote {REPORTS_DIR / 'quantitative_report_evidence.txt'}")
    print(f"Wrote {FIGURES_DIR / PRODUCTION_FIELDS_FIGURE}")
    print(f"Wrote {FIGURES_DIR / INTERFACE_TRACKING_FIGURE}")
    print(f"Wrote {FIGURES_DIR / GRID_SENSITIVITY_FIGURE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
