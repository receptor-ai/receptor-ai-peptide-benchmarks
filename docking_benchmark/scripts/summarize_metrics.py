#!/usr/bin/env python3
"""Aggregate a metrics_results.csv (produced by score_docking.py) into summary
statistics: per-metric mean and the CAPRI class distribution."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

METRIC_COLUMNS = [
    "n_atoms",
    "backbone_rmsd", "heavy_atom_rmsd", "ca_rmsd", "sidechain_rmsd",
    "clash_native", "clash_pred",
    "fnat", "fnonnat", "LRMSD", "DockQ",
]

CAPRI_ORDER = ["Incorrect", "Acceptable", "Medium", "High"]


def load_rows(results_csv: Path) -> list[dict]:
    with results_csv.open(newline="") as f:
        return list(csv.DictReader(f))


def metric_stats(rows: list[dict], column: str) -> dict | None:
    raw = [float(row[column]) for row in rows if row[column] != ""]
    values = np.array(raw)
    values = values[np.isfinite(values)]  # some peptides have no sidechain atoms -> NaN
    if values.size == 0:
        return None
    return {"n": values.size, "mean": float(np.mean(values))}


def capri_distribution(rows: list[dict]) -> dict[str, int]:
    counts = {c: 0 for c in CAPRI_ORDER}
    for row in rows:
        counts[row["CAPRI_class"]] += 1
    return counts


def print_report(results_csv: Path, rows: list[dict]) -> None:
    print(f"{results_csv}  ({len(rows)} scored cases)\n")

    print(f"{'metric':<20} {'n':>5} {'mean':>10}")
    for column in METRIC_COLUMNS:
        stats = metric_stats(rows, column)
        if stats is None:
            print(f"{column:<20} {'--':>5}")
            continue
        print(f"{column:<20} {stats['n']:>5} {stats['mean']:>10.3f}")

    print("\nCAPRI class distribution:")
    counts = capri_distribution(rows)
    for capri_class in CAPRI_ORDER:
        n = counts[capri_class]
        pct = 100 * n / len(rows) if rows else 0.0
        print(f"  {capri_class:<12} {n:>5}  ({pct:>5.1f}%)")
    acceptable_or_better = sum(counts[c] for c in ("Acceptable", "Medium", "High"))
    pct = 100 * acceptable_or_better / len(rows) if rows else 0.0
    print(f"  {'>=Acceptable':<12} {acceptable_or_better:>5}  ({pct:>5.1f}%)")


def write_summary_csv(summary_csv: Path, rows: list[dict]) -> None:
    with summary_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "n", "mean"])
        for column in METRIC_COLUMNS:
            stats = metric_stats(rows, column)
            if stats is None:
                continue
            writer.writerow([column, stats["n"], stats["mean"]])
        writer.writerow([])
        writer.writerow(["CAPRI_class", "n", "pct"])
        counts = capri_distribution(rows)
        for capri_class in CAPRI_ORDER:
            n = counts[capri_class]
            pct = 100 * n / len(rows) if rows else 0.0
            writer.writerow([capri_class, n, round(pct, 1)])


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize a metrics_results.csv from score_docking.py.")
    parser.add_argument("results_csv", type=Path, help="per-case metrics CSV to summarize")
    args = parser.parse_args()

    rows = load_rows(args.results_csv)
    if not rows:
        print(f"No rows found in {args.results_csv}")
        return

    print_report(args.results_csv, rows)

    summary_csv = args.results_csv.with_name(f"{args.results_csv.stem}_summary.csv")
    write_summary_csv(summary_csv, rows)
    print(f"\nWrote summary -> {summary_csv}")


if __name__ == "__main__":
    main()
