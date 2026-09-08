#!/usr/bin/env python3
"""Load the BMI-200 and BMI-MODES peptide–protein datasets.

This is a small example loader for the two datasets in this repository. It needs
only pandas. To read the mmCIF coordinates, install gemmi and use the paths that
this script returns.

Usage:
    python scripts/load.py            # print a summary of both datasets
    python scripts/load.py bmi200     # summary of BMI-200 only
    python scripts/load.py bmi_modes  # summary of BMI-MODES only

In your own code:
    from scripts.load import load_bmi200, load_bmi_modes
    entries = load_bmi200()           # one row per entry
    sites, poses = load_bmi_modes()   # one row per site, one row per pose
"""
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Always load the CSV files this way. keep_default_na=False stops the element
# symbol NA (sodium) from becoming a missing value.
_CSV = dict(low_memory=False, keep_default_na=False, na_values=[""])


def load_bmi200():
    """Return the 202 BMI-200 entries as one row per entry."""
    return pd.read_csv(os.path.join(ROOT, "bmi200", "MANIFEST.csv"), **_CSV)


def load_bmi_modes():
    """Return the BMI-MODES sites and poses as (sites, poses)."""
    base = os.path.join(ROOT, "bmi_modes")
    sites = pd.read_csv(os.path.join(base, "MANIFEST.csv"), **_CSV)
    poses = pd.read_csv(os.path.join(base, "poses.csv"), **_CSV)
    return sites, poses


def bmi200_complex_path(pdb_id):
    """Return the path to one BMI-200 complex.cif."""
    return os.path.join(ROOT, "bmi200", "entries", pdb_id, "complex.cif")


def bmi_modes_pose_path(site_id, pose):
    """Return the path to one BMI-MODES pose file (for example pose='pose_01')."""
    return os.path.join(ROOT, "bmi_modes", "sites", site_id, f"{pose}.cif")


def _summary_bmi200():
    e = load_bmi200()
    print("BMI-200")
    print(f"  entries: {len(e)}")
    print(f"  cyclic / linear: "
          f"{int(e.cell.str.startswith('cyclic').sum())} / "
          f"{int(e.cell.str.startswith('linear').sum())}")
    n = pd.to_numeric(e.n_res_observed, errors="coerce")
    print(f"  peptide length: {int(n.min())}-{int(n.max())} residues "
          f"(median {n.median():.0f})")
    print(f"  distinct receptors (UniProt): {e.receptor_uniprot.nunique()}")
    print(f"  example complex: {bmi200_complex_path(e.pdb_id.iloc[0])}")


def _summary_bmi_modes():
    sites, poses = load_bmi_modes()
    print("BMI-MODES")
    print(f"  sites: {len(sites)}")
    print(f"  poses: {len(poses)}")
    print(f"  distinct receptors (UniProt): {sites.uniprot.nunique()}")
    row = poses.iloc[0]
    print(f"  example pose: {bmi_modes_pose_path(row.site_id, row.pose)}")


def main(argv):
    which = argv[1] if len(argv) > 1 else "all"
    if which in ("all", "bmi200"):
        _summary_bmi200()
    if which in ("all", "bmi_modes"):
        _summary_bmi_modes()
    if which not in ("all", "bmi200", "bmi_modes"):
        print(__doc__)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
