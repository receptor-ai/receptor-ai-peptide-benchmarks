# docking985: a peptide–protein re-docking benchmark

docking985 is a curated dataset of 985 peptide–protein X-ray/NMR complexes, drawn from seven
published peptide-docking benchmarks and filtered to a consistent quality bar (resolution, peptide
size, and how well-packed the peptide is against its receptor). Every entry ships the native
(bound) structure, split into a receptor file and a peptide file that recombine atom for atom into
the complex, plus a bond-order-aware SDF of the peptide for use as a docking input.

This project is a benchmark, not a docking method: it ships structures and a scoring library, not
predicted poses. Bring your own docking method — of any kind, for any peptide representation — dock
the shipped peptides against the shipped receptors, and score the result against the native
structures with the tools in `scripts/`.

## Layout

```
docking_benchmark/
  README.md        this file
  LICENSE          CC-BY-4.0 - covers the docking985 dataset
  SHA256SUMS       a SHA-256 checksum for every delivered file
  requirements.txt minimum versions for scripts/ - see below
  scripts/         quality_filter.py, score_docking.py, summarize_metrics.py, and their own
                    LICENSE (MIT)
    metrics/       the scoring library (RMSD variants, DockQ/CAPRI, an approximate clashscore),
                    covered by scripts/LICENSE too
  docking985/      the dataset. See docking985/README.md.
```

## Scripts

- **`scripts/quality_filter.py`** — the filter used to select every entry, self-contained and
  runnable directly against a shipped `entries/<case_id>/` (or your own `receptor.pdb` +
  `peptide.pdb` + `meta.json`), no external downloads needed. Use it to quality-screen a candidate
  complex of your own before adding it to a benchmark.
- **`scripts/score_docking.py`** — scores a folder of predicted-pose SDFs (one per case, matched by
  filename to a `case_id`) against the native structures: backbone/CA/sidechain/heavy-atom RMSD,
  DockQ (fnat, fnonnat, iRMSD, LRMSD, CAPRI class), and an approximate clashscore. Hydrogens in a
  predicted-pose SDF are optional: the scorer strips them itself (RDKit's `Chem.RemoveAllHs`) and
  compares heavy atoms only, so a full-atom SDF (like `peptide.sdf` itself) scores the same as a
  heavy-atom one. Atom *order* need not match -
  predicted atoms are matched to the native by substructure, since docking tools routinely reorder
  atoms on the way out (PDBQT round-trips, OpenBabel, etc.). Note on `iRMSD`: the receptor is never
  docked here (only the peptide is), so `iRMSD`'s receptor side is always error-free by
  construction — it effectively reads as a local RMSD of the peptide's own interface atoms, not
  the joint receptor+ligand accuracy it measures in general protein-protein docking. Note on
  `LRMSD` vs `heavy_atom_rmsd`: these are not the same measurement. `LRMSD` comes from DockQ and is
  computed from **backbone atoms only** (it never looks at sidechains, so it can't distinguish a
  correct sidechain from a wrong one, symmetric or not). `heavy_atom_rmsd` (and `sidechain_rmsd`)
  are computed here, directly, over **all heavy atoms including sidechains**. Predicted atoms are
  matched to the native by enumerating every valid substructure match (not just the first one RDKit
  happens to find) and keeping whichever gives the lowest RMSD - this handles any local symmetry
  (a Val/Leu methyl pair, a Tyr/Phe ring, a symmetric non-standard group) automatically, since a
  "wrong" but chemically valid match is just another candidate that loses to a better one, without
  needing a list of known symmetric residue types. Cheap here (typically 16 matches per entry, at
  most 1024, the full dataset in about a second) precisely because it matches one specific molecule
  against a reordered copy of itself, not a general automorphism search across unrelated inputs.
  Neither `LRMSD` nor `heavy_atom_rmsd` substitutes for the other.
- **`scripts/summarize_metrics.py`** — aggregates a `score_docking.py` results CSV into per-metric
  means and the CAPRI class distribution.

```bash
python scripts/quality_filter.py docking985/entries/1B6J_entity2
python scripts/score_docking.py path/to/your_docking_outputs/
python scripts/summarize_metrics.py path/to/your_docking_outputs_metrics_results.csv
```

Dependencies (`requirements.txt`, minimum versions this package was built and tested against):
`biotite`, `numpy`, `rdkit`, and (for `score_docking.py`'s DockQ metrics only - fnat, fnonnat,
LRMSD, DockQ, CAPRI_class) [`DockQ`](https://github.com/bjornwallner/DockQ). Without it,
`score_docking.py` still runs and reports the RMSD variants and clashscore; it prints a note and
leaves the DockQ columns blank rather than failing.

```bash
pip install -r requirements.txt
```

## Provenance and licensing

Two separate licenses apply. The **`docking985` dataset** (structures and metadata) is under
[CC-BY-4.0](LICENSE). The **code** (`scripts/`, including `scripts/metrics/`) is under
[MIT](scripts/LICENSE) — use, modify, and redistribute it freely, including in closed-source tools.

Coordinates come from the RCSB PDB under CC0 1.0. Attribute the RCSB PDB when you use this
benchmark. See [`docking985/README.md`](docking985/README.md) for the source benchmarks each entry
was drawn from.

## How to cite

If you use docking985, please cite this benchmark:

> Receptor.AI, Inc. (2026). docking985: a peptide–protein re-docking benchmark (Version 1.0.0).
> https://github.com/receptor-ai/receptor-ai-peptide-benchmarks
