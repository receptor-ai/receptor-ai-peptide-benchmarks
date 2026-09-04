#!/usr/bin/env python3
"""Score a directory of docking outputs against the native poses shipped in
docking985/entries/. Reusable against any docking method: it only needs one
SDF per predicted pose, heavy atoms only, with the same atom count and atom
order as the corresponding entry's peptide.pdb (equivalently: peptide.sdf
with its explicit hydrogens stripped, e.g. via RDKit's Chem.RemoveHs) --
peptide.sdf itself is a docking *input* template, not the shape a predicted
pose is expected to come back in.

Name each output SDF so it can be matched back to a case: either the case_id
itself (optionally with a suffix, e.g. "1B6J_entity2.mymethod.sdf"), or the
bare pdb_id when that pdb_id has exactly one entry in the benchmark.

Usage:
    python score_docking.py path/to/docking_outputs/
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np
import biotite.structure.io as strucio

from metrics.structures import combine, load_native_peptide, load_pred_peptide_coord
from metrics.contact_weighted_rmsd import contact_weighted_rmsd
from metrics.rmsd_and_capri import backbone_rmsd, ca_rmsd, dockq_metrics, heavy_atom_rmsd, sidechain_rmsd
from metrics.utils import clashscore

ENTRIES_DIR = ROOT / "docking985" / "entries"
MANIFEST_CSV = ROOT / "docking985" / "MANIFEST.csv"

ROUND_NDIGITS = 3

FIELDNAMES = [
    "case_id", "n_atoms",
    "backbone_rmsd", "heavy_atom_rmsd", "ca_rmsd", "sidechain_rmsd",
    "contact_weighted_rmsd",
    "clash_native", "clash_pred",
    "fnat", "fnonnat", "LRMSD", "DockQ", "CAPRI_class",
]


def load_manifest_rows() -> list[dict]:
    with MANIFEST_CSV.open(newline="") as f:
        return list(csv.DictReader(f))


def unambiguous_pdb_id_lookup(rows: list[dict]) -> dict[str, str]:
    """pdb_id -> case_id, only where the pdb_id has exactly one entry (unambiguous)."""
    by_pdb_id: dict[str, list[dict]] = {}
    for row in rows:
        by_pdb_id.setdefault(row["pdb_id"], []).append(row)
    return {pdb_id: r[0]["case_id"] for pdb_id, r in by_pdb_id.items() if len(r) == 1}


def case_id_from_sdf(sdf_path: Path, case_ids: set[str], pdb_id_lookup: dict[str, str]) -> str | None:
    stem = sdf_path.stem
    if stem in case_ids:
        return stem
    for case_id in case_ids:
        if stem.startswith(case_id + "."):
            return case_id
    pdb_id = stem.split("_")[0].split(".")[0]
    return pdb_id_lookup.get(pdb_id)


def score_case(entry_dir: Path, output_sdf: Path) -> dict:
    receptor = strucio.load_structure(str(entry_dir / "receptor.pdb"))
    native_peptide = load_native_peptide(entry_dir / "peptide.pdb")
    peptide_pred = load_pred_peptide_coord(output_sdf)

    if len(native_peptide) != len(peptide_pred):
        raise ValueError(f"atom count mismatch: native={len(native_peptide)} pred={len(peptide_pred)}")

    rec_coord = receptor.coord
    pep_native_coord = native_peptide.coord

    complex_native = combine(receptor, native_peptide, pep_native_coord)
    complex_pred = combine(receptor, native_peptide, peptide_pred)

    row = {
        "case_id": entry_dir.name,
        "n_atoms": len(pep_native_coord),
        "backbone_rmsd": backbone_rmsd(native_peptide, peptide_pred, rec_coord, rec_coord.copy()),
        "heavy_atom_rmsd": heavy_atom_rmsd(native_peptide, peptide_pred, rec_coord, rec_coord.copy()),
        "ca_rmsd": ca_rmsd(native_peptide, peptide_pred, rec_coord, rec_coord.copy()),
        "sidechain_rmsd": sidechain_rmsd(native_peptide, peptide_pred, rec_coord, rec_coord.copy()),
        "contact_weighted_rmsd": contact_weighted_rmsd(pep_native_coord, peptide_pred, rec_coord, rec_coord.copy()),
        "clash_native": clashscore(complex_native),
        "clash_pred": clashscore(complex_pred),
    }
    row.update(dockq_metrics(receptor, native_peptide, peptide_pred))
    return row


def rounded(row: dict) -> dict:
    return {
        k: (round(float(v), ROUND_NDIGITS) if isinstance(v, (float, np.floating)) else v)
        for k, v in row.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outputs_dir", type=Path, help="folder of *.sdf docking outputs")
    args = parser.parse_args()

    results_csv = args.outputs_dir.parent / f"{args.outputs_dir.name}_metrics_results.csv"

    rows = load_manifest_rows()
    pdb_id_lookup = unambiguous_pdb_id_lookup(rows)
    case_ids = {row["case_id"] for row in rows}

    output_sdfs = sorted(args.outputs_dir.glob("*.sdf"))
    if not output_sdfs:
        print(f"No docking outputs found in {args.outputs_dir}")
        return

    results = []
    for sdf_path in output_sdfs:
        case_id = case_id_from_sdf(sdf_path, case_ids, pdb_id_lookup)
        entry_dir = ENTRIES_DIR / case_id if case_id else None
        if entry_dir is None or not entry_dir.exists():
            print(f"{sdf_path.name}: SKIPPED -- not in the benchmark, or ambiguous (name it <case_id>.sdf)")
            continue

        try:
            row = rounded(score_case(entry_dir, sdf_path))
        except Exception as error:
            print(f"{case_id}: SKIPPED -- {error}")
            continue

        results.append(row)
        print(f"{case_id:<20} DockQ={row['DockQ']} ({row['CAPRI_class']:<10}) "
              f"heavy_atom_RMSD={row['heavy_atom_rmsd']}  clash(native/pred)="
              f"{row['clash_native']}/{row['clash_pred']}")

    with results_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    print(f"\nScored {len(results)}/{len(output_sdfs)} available docking outputs -> {results_csv}")


if __name__ == "__main__":
    main()
