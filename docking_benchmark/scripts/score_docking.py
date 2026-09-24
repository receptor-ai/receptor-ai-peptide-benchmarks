#!/usr/bin/env python3
"""Score a directory of docking outputs against the native poses shipped in
docking985/entries/. Reusable against any docking method: it only needs one
SDF per predicted pose, heavy atoms only, with the same connectivity as the
corresponding entry's peptide.sdf (equivalently: peptide.sdf with its
explicit hydrogens stripped, e.g. via RDKit's Chem.RemoveAllHs) - peptide.sdf
itself is a docking *input* template, not the shape a predicted pose is
expected to come back in. Atom *order* does not need to match: predicted
atoms are matched to the native by substructure, since docking tools
routinely reorder atoms on the way out.

Name each output SDF so it can be matched back to a case: either the case_id
itself (optionally with a suffix, e.g. "1B6J_entity2.mymethod.sdf"), or the
bare pdb_id when that pdb_id has exactly one entry in the benchmark.

Usage:
    python score_docking.py path/to/docking_outputs/
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np
import biotite.structure.io as strucio
from rdkit import Chem

from metrics.structures import combine, load_native_peptide, load_pred_peptide_coord
from metrics.rmsd_and_capri import backbone_rmsd, ca_rmsd, dockq_metrics, heavy_atom_rmsd, sidechain_rmsd
from metrics.utils import clashscore

# checked once at import time (find_spec, not a real import) so a missing DockQ
# degrades to "no DockQ columns" instead of crashing every case's scoring
DOCKQ_AVAILABLE = importlib.util.find_spec("DockQ") is not None

ENTRIES_DIR = ROOT / "docking985" / "entries"
MANIFEST_CSV = ROOT / "docking985" / "MANIFEST.csv"

ROUND_NDIGITS = 3

FIELDNAMES = [
    "case_id", "n_atoms",
    "backbone_rmsd", "heavy_atom_rmsd", "ca_rmsd", "sidechain_rmsd",
    "clash_native", "clash_pred",
    "fnat", "fnonnat", "iRMSD", "LRMSD", "DockQ", "CAPRI_class",
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
    # some receptor.pdb files carry explicit hydrogens; clashscore() is heavy-atom-only
    # by design, so strip them or clash_native/clash_pred stop being comparable across cases
    receptor = receptor[receptor.element != "H"]
    native_peptide = load_native_peptide(entry_dir / "peptide.pdb")
    # RemoveAllHs, not RemoveHs - see load_pred_peptide_coord()
    native_mol = Chem.RemoveAllHs(Chem.SDMolSupplier(str(entry_dir / "peptide.sdf"), removeHs=False)[0])
    peptide_pred = load_pred_peptide_coord(output_sdf, native_mol)

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
        "clash_native": clashscore(complex_native),
        "clash_pred": clashscore(complex_pred),
    }
    if DOCKQ_AVAILABLE:
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

    if not DOCKQ_AVAILABLE:
        print("NOTE: DockQ not installed - fnat, fnonnat, LRMSD, DockQ, CAPRI_class will be blank. "
              "RMSD variants and clashscore are unaffected. `pip install DockQ` to get them.\n")

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
            print(f"{sdf_path.name}: SKIPPED - not in the benchmark, or ambiguous (name it <case_id>.sdf)")
            continue

        try:
            row = rounded(score_case(entry_dir, sdf_path))
        except MemoryError:
            print(f"{case_id}: SKIPPED - OUT OF MEMORY (large complex; not a data/scoring error)")
            continue
        except Exception as error:
            print(f"{case_id}: SKIPPED - {error}")
            continue

        results.append(row)
        dockq_part = f"DockQ={row['DockQ']} ({row['CAPRI_class']:<10}) " if DOCKQ_AVAILABLE else ""
        print(f"{case_id:<20} {dockq_part}"
              f"heavy_atom_RMSD={row['heavy_atom_rmsd']}  clash(native/pred)="
              f"{row['clash_native']}/{row['clash_pred']}")

    with results_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    print(f"\nScored {len(results)}/{len(output_sdfs)} available docking outputs -> {results_csv}")


if __name__ == "__main__":
    main()
