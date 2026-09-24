#!/usr/bin/env python3
"""Recompute, directly from a shipped entry's coordinates, the same quality
statistics and pass/fail thresholds used to select every case in docking985 -
no external downloads needed, since receptor.pdb/peptide.pdb/meta.json already
carry everything the check needs. Point it at your own receptor.pdb +
peptide.pdb + meta.json (resolution_angstrom, is_nmr) to quality-screen a
candidate complex the same way before adding it to a benchmark of your own.

A peptide atom counts as "in contact" with the receptor if its nearest
receptor atom is within CONTACT_DISTANCE_THRESHOLD - the fraction of peptide
atoms in contact is the selection criterion, not the mean/std of per-atom
distances (a single solvent-exposed or flexible atom can drag a mean/std past
threshold on its own without the peptide actually being solvent-exposed).

Note on entries with more than one crystallographic peptide copy
(MANIFEST.csv's n_peptide_chain_copies > 1, ~40% of the dataset): receptor.pdb
ships only the receptor copy nearest the shipped peptide chain, not the full
deposited assembly. Selection itself was decided against the full assembly
(MANIFEST.csv's peptide_contact_fraction reflects that), so this script can
report a lower contact fraction - and even FAIL - on those cases when run
against the shipped files alone. That reflects the smaller receptor shipped
here, not a change to which entries belong in the benchmark. See
docking985/README.md.

Usage:
    python quality_filter.py docking985/entries/1B6J_entity2
    python quality_filter.py docking985/entries               # every entry, recheck
    python quality_filter.py docking985/entries --csv out.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import biotite.structure.io as strucio

CONTACT_DISTANCE_THRESHOLD = 5.0  # A: a peptide atom counts as "in contact" within this of the receptor
VDW_CLASH_FRACTION = 0.7          # a pair violates VDW if closer than this fraction of (r1 + r2)

MAX_PEPTIDE_LENGTH = 20
MAX_RESOLUTION = 2.5              # A; not applied to NMR structures
MAX_DIST_MIN_POCKET = 3.65        # A
MIN_PEPTIDE_CONTACT_FRACTION = 0.5
MAX_VDW_PER_ATOM = 0.1

VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
    "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98,
}


def vdw_radii(elements: np.ndarray) -> np.ndarray:
    return np.asarray([VDW_RADII.get(str(e).upper(), 1.70) for e in elements], dtype=np.float32)


def heavy_atoms(pdb_path: Path):
    atoms = strucio.load_structure(str(pdb_path), model=1)
    return atoms[atoms.element != "H"]


def distance_stats(peptide_coord: np.ndarray, peptide_element: np.ndarray,
                    protein_coord: np.ndarray, protein_element: np.ndarray) -> dict:
    peptide_vdw = vdw_radii(peptide_element)
    protein_vdw = vdw_radii(protein_element)

    min_per_peptide_atom = []
    global_min, global_max, n_vdw_violations = math.inf, 0.0, 0
    for i, coord in enumerate(peptide_coord):
        dist = np.linalg.norm(protein_coord - coord.reshape(1, 3), axis=1)
        min_per_peptide_atom.append(float(dist.min()))
        global_min = min(global_min, float(dist.min()))
        global_max = max(global_max, float(dist.max()))
        n_vdw_violations += int((dist < (peptide_vdw[i] + protein_vdw) * VDW_CLASH_FRACTION).sum())

    contact_fraction = float((np.asarray(min_per_peptide_atom) <= CONTACT_DISTANCE_THRESHOLD).mean())
    return {
        "dist_min_pocket": round(global_min, 3),
        "dist_max_pocket": round(global_max, 3),
        "peptide_contact_fraction": round(contact_fraction, 3),
        "n_vdw_violations": n_vdw_violations,
    }


def evaluate_entry(entry_dir: Path) -> dict:
    meta = json.loads((entry_dir / "meta.json").read_text())
    peptide = heavy_atoms(entry_dir / "peptide.pdb")
    receptor = heavy_atoms(entry_dir / "receptor.pdb")

    stats = distance_stats(peptide.coord, peptide.element, receptor.coord, receptor.element)
    peptide_length = len(set(zip(peptide.chain_id.tolist(), peptide.res_id.tolist())))
    vdw_per_atom = stats["n_vdw_violations"] / len(peptide) if len(peptide) else None

    resolution = meta.get("resolution_angstrom")
    is_nmr = bool(meta.get("is_nmr"))

    reasons = []
    if peptide_length > MAX_PEPTIDE_LENGTH:
        reasons.append(f"peptide_length={peptide_length} > {MAX_PEPTIDE_LENGTH}")
    if not is_nmr and (resolution is None or resolution > MAX_RESOLUTION):
        reasons.append(f"resolution={resolution} > {MAX_RESOLUTION} (or missing)")
    if stats["dist_min_pocket"] > MAX_DIST_MIN_POCKET:
        reasons.append(f"dist_min_pocket={stats['dist_min_pocket']} > {MAX_DIST_MIN_POCKET}")
    if stats["peptide_contact_fraction"] < MIN_PEPTIDE_CONTACT_FRACTION:
        reasons.append(f"peptide_contact_fraction={stats['peptide_contact_fraction']} < {MIN_PEPTIDE_CONTACT_FRACTION}")
    if vdw_per_atom is None or vdw_per_atom > MAX_VDW_PER_ATOM:
        reasons.append(f"vdw_violations_per_atom={vdw_per_atom} > {MAX_VDW_PER_ATOM}")

    return {
        "case_id": entry_dir.name,
        "n_peptide_atoms": len(peptide),
        "n_protein_atoms": len(receptor),
        "peptide_length": peptide_length,
        **stats,
        "vdw_violations_per_peptide_atom": None if vdw_per_atom is None else round(vdw_per_atom, 4),
        "passes_filter": not reasons,
        "fail_reasons": ";".join(reasons),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("path", type=Path, help="an entries/<case_id> dir, or an entries/ dir to check them all")
    parser.add_argument("--csv", type=Path, help="write per-entry results to this CSV")
    args = parser.parse_args()

    if (args.path / "meta.json").exists():
        entry_dirs = [args.path]
    else:
        entry_dirs = sorted(d for d in args.path.iterdir() if d.is_dir())

    results = []
    for entry_dir in entry_dirs:
        try:
            results.append(evaluate_entry(entry_dir))
        except Exception as error:
            print(f"{entry_dir.name}: SKIPPED - {error}")

    for r in results:
        status = "PASS" if r["passes_filter"] else f"FAIL ({r['fail_reasons']})"
        print(f"{r['case_id']:<20} {status}")

    n_pass = sum(r["passes_filter"] for r in results)
    print(f"\n{n_pass}/{len(results)} entries pass the quality filter")

    if args.csv and results:
        with args.csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            writer.writeheader()
            writer.writerows(results)
        print(f"Wrote {args.csv}")


if __name__ == "__main__":
    main()
