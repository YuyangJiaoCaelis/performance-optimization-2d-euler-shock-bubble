#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math


@dataclass
class Snapshot:
    path: Path
    metadata: dict[str, str]
    nx: int
    ny: int
    x_values: list[float]
    y_values: list[float]
    fields: dict[str, list[float]]

    def field(self, name: str, i: int, j: int) -> float:
        return self.fields[name][j * self.nx + i]


def _close(left: float, right: float, tolerance: float = 1.0e-12) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def _parse_metadata(line: str) -> tuple[str, str] | None:
    text = line[1:].strip()
    if "=" not in text:
        return None
    key, value = text.split("=", 1)
    return key.strip(), value.strip()


def load_snapshot(path: str | Path) -> Snapshot:
    snapshot_path = Path(path)
    metadata: dict[str, str] = {}
    header: list[str] | None = None
    columns: dict[str, list[float]] = {}

    with snapshot_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                parsed = _parse_metadata(line)
                if parsed is not None:
                    key, value = parsed
                    metadata[key] = value
                continue

            if header is None:
                header = [entry.strip() for entry in line.split(",")]
                columns = {name: [] for name in header}
                continue

            values = [float(entry) for entry in line.split(",")]
            if len(values) != len(header):
                raise ValueError(f"Row in {snapshot_path} has {len(values)} columns, expected {len(header)}.")

            for name, value in zip(header, values):
                columns[name].append(value)

    if header is None:
        raise ValueError(f"No tabular data found in {snapshot_path}.")

    if "x" not in columns or "y" not in columns:
        raise ValueError(f"{snapshot_path} does not contain x and y columns.")

    row_count = len(columns["x"])
    if row_count == 0:
        raise ValueError(f"{snapshot_path} does not contain any cell data.")

    nx = int(metadata["nx"]) if "nx" in metadata else 0
    ny = int(metadata["ny"]) if "ny" in metadata else 0

    if nx <= 0 or ny <= 0:
        first_y = columns["y"][0]
        nx = 0
        for y_value in columns["y"]:
            if not _close(y_value, first_y):
                break
            nx += 1
        if nx == 0 or row_count % nx != 0:
            raise ValueError(f"Could not infer a structured grid layout from {snapshot_path}.")
        ny = row_count // nx

    if nx * ny != row_count:
        raise ValueError(
            f"Metadata/grid shape mismatch in {snapshot_path}: nx * ny = {nx * ny}, rows = {row_count}."
        )

    x_values = columns["x"][:nx]
    y_values = [columns["y"][j * nx] for j in range(ny)]

    for j in range(ny):
        start = j * nx
        stop = start + nx
        row_x = columns["x"][start:stop]
        row_y = columns["y"][start:stop]

        for expected_x, actual_x in zip(x_values, row_x):
            if not _close(expected_x, actual_x):
                raise ValueError(f"Inconsistent x coordinates in row {j} of {snapshot_path}.")

        for actual_y in row_y:
            if not _close(y_values[j], actual_y):
                raise ValueError(f"Inconsistent y coordinates in row {j} of {snapshot_path}.")

    fields = {name: values for name, values in columns.items() if name not in {"x", "y"}}
    return Snapshot(snapshot_path, metadata, nx, ny, x_values, y_values, fields)


def relative_error_metrics(candidate: list[float], reference: list[float]) -> dict[str, float]:
    if len(candidate) != len(reference):
        raise ValueError("Candidate and reference arrays must have the same length.")

    abs_diff_sum = 0.0
    abs_ref_sum = 0.0
    diff_square_sum = 0.0
    ref_square_sum = 0.0
    max_diff = 0.0
    max_ref = 0.0

    for candidate_value, reference_value in zip(candidate, reference):
        difference = candidate_value - reference_value
        abs_diff = abs(difference)
        abs_ref = abs(reference_value)

        abs_diff_sum += abs_diff
        abs_ref_sum += abs_ref
        diff_square_sum += difference * difference
        ref_square_sum += reference_value * reference_value
        max_diff = max(max_diff, abs_diff)
        max_ref = max(max_ref, abs_ref)

    denominator_l1 = max(abs_ref_sum, 1.0e-30)
    denominator_l2 = max(math.sqrt(ref_square_sum), 1.0e-30)
    denominator_linf = max(max_ref, 1.0e-30)

    return {
        "l1": abs_diff_sum / denominator_l1,
        "l2": math.sqrt(diff_square_sum) / denominator_l2,
        "linf": max_diff / denominator_linf,
        "max_abs": max_diff,
    }
