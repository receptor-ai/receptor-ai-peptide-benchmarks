# docking985: a peptide–protein docking benchmark

docking985 holds 985 peptide–protein X-ray/NMR complexes, curated from seven published
peptide-docking benchmarks and filtered to a consistent quality bar. Each entry ships the native
(bound) structure only, split into a receptor and a peptide file that recombine atom for atom into
the complex — there are no decoy or docked poses here. Bring your own docking method, generate
predicted poses for the shipped peptides, and score them against these native structures with the
tools in `../scripts/`.

```
docking985/
  README.md            this file
  DATA_DICTIONARY.md   definition of every column in MANIFEST.csv
  MANIFEST.csv          one row per entry, all metadata
  entries/<case_id>/    complex.pdb, receptor.pdb, peptide.pdb, peptide.sdf, meta.json
```

Every entry directory is self-contained: `meta.json` repeats that entry's `MANIFEST.csv` row, and
`peptide.pdb` plus `receptor.pdb` equals `complex.pdb`, atom for atom. Coordinates are plain PDB;
`peptide.sdf` carries the same atoms plus perceived bond orders and explicit hydrogens, for use as
a docking input. Load any of them with **biotite** (`biotite.structure.io.load_structure`) or
**RDKit** (the SDF). `peptide.sdf` is an input template; `../scripts/score_docking.py` compares
heavy atoms only and strips any hydrogens from a predicted SDF itself (RDKit's `Chem.RemoveAllHs`),
so a pose written in `peptide.sdf`'s full-atom form scores the same as a heavy-atom one.
Every `peptide.pdb` is heavy-atom only, but 92 of 985 `receptor.pdb`/`complex.pdb` files carry
explicit hydrogens (12 NMR + 80 X-ray, as deposited) — strip them yourself if your tooling assumes
heavy-atom-only receptors throughout.

## 1. Composition

| Property | Value |
|---|---|
| Entries | 985 |
| Distinct PDB entries | 955 (28 contain 2 independent peptide-binding sites and 1 — 4TR9 — contains 3, each kept as a separate entry) |
| Peptide length (observed residues) | 2–20, median 8 |
| Resolution (Å) | 0.75–2.50, median 1.90 (X-ray only) |
| Experimental method | 973 X-ray, 12 solution NMR |
| Distinct receptor UniProt accessions | 473 |
| Entries with >1 crystallographic peptide copy | 382 (see §3) |

### Source benchmarks

Every entry is traced back to the published benchmark(s) it was drawn from (`benchmarks` column,
`;`-separated where an entry appears in more than one):

| Source | Entries |
|---|---|
| CSPSet | 390 |
| ProtPep37_2021 | 342 |
| PepSet | 136 |
| PPDbench | 117 |
| LEADS-PEP | 53 |
| Vina_47 | 44 |
| PIPER-FlexPepDock | 41 |

## 2. Selection criteria

Every entry passed all five checks below, applied to the full deposited assembly (`../scripts/quality_filter.py`
documents the exact thresholds and recomputes them from a shipped entry's own coordinates):

- observed peptide length ≤ 20 residues;
- resolution ≤ 2.5 Å (not applied to the 12 NMR entries);
- closest peptide–receptor heavy-atom distance ≤ 3.65 Å;
- ≥ 50% of peptide heavy atoms within 5 Å of the receptor (`peptide_contact_fraction`) — chosen over
  the mean/std of per-atom distances, which a single solvent-exposed or flexible peptide atom can
  drag past a threshold on its own without the peptide actually being solvent-exposed;
- ≤ 0.1 severe Van-der-Waals violations per peptide heavy atom.

## 3. Multiple crystallographic peptide copies

382 of the 985 entries deposit more than one copy of the same peptide–receptor pair in the
asymmetric unit (`n_peptide_chain_copies` in `MANIFEST.csv`). For these, `entries/<case_id>/` keeps
a single representative pair: the first listed peptide chain, paired with whichever receptor chain
sits nearest to it. `dist_min_pocket`, `dist_max_pocket`, `peptide_contact_fraction` and the VDW
columns in `MANIFEST.csv` still reflect the original selection-time computation against the *full*
deposited assembly (every copy, every receptor chain) — they are not recomputed from the smaller,
single-copy `receptor.pdb` shipped here.

One consequence: running `quality_filter.py` against a shipped entry's own files can read a lower
`peptide_contact_fraction` than `MANIFEST.csv` reports, and can even fail the ≥0.5 threshold, for
entries where the peptide's contacts spread across more than one receptor copy. This happens for 58
of the 382 multi-copy entries. It reflects the smaller receptor shipped for docking, not a reason to
doubt the entry's inclusion — all 985 entries passed selection against the true, full assembly. A
further 4 entries fail the recheck for reasons unrelated to multi-copy receptor scope, from minor
differences between the original curation pipeline and `quality_filter.py`'s independent
reimplementation: `6ZCD_entity1` misses the 0.5 contact-fraction cutoff by 0.004; `4TR9_entity2` and
`4TR9_entity3` miss it by 0.021 and 0.026 respectively; `6D3Z_entity2` fails a different criterion
entirely (the VDW-violations-per-atom threshold, not contact fraction). `MANIFEST.csv` is always the
authoritative record of why an entry was selected.

## 4. Notes before you use this dataset

1. **Non-canonical and modified peptides.** Several source benchmarks (CSPSet in particular)
   include peptides built from non-standard monomers; these are kept as `HETATM` records in
   `peptide.pdb`/`complex.pdb` and are genuine peptide chemistry, not crystallization additives.
   `peptide_sequence` does **not** reliably mark these with `X` — a non-standard residue reads as
   `X` only when it has no recognized canonical parent; otherwise it silently reads as the parent's
   letter and the modification is invisible in the sequence string alone (e.g. 4YV9_entity2 carries
   9 non-canonical components — a D-amino acid, four N-methyl-leucines, an N-methyl-valine, BMT,
   sarcosine, an Abu — and its `peptide_sequence` is `ALLVTAGLVLA`, zero `X`). Use
   `noncanonical_components` (backbone-order CCD codes) for the real chemistry. This dataset does
   not identify whether the peptide is cyclic.
2. **`peptide.sdf`'s bond perception is automated** (RDKit from 3D coordinates, with an OpenBabel
   fallback) and was not manually reviewed per entry. RDKit only has a built-in template for the 20
   canonical residues; for anything else it guesses bond order from geometry, which used to
   mis-assign the backbone carbonyl (C=O read as C-O, then a spurious H added to balance the
   valence) on non-canonical residues - corrected for 324 entries (1048 bonds) by forcing the
   correct order wherever geometry unambiguously calls for a C=O. 5 entries (`2B9H_entity2`,
   `5N22_entity2`, `6T2F_entity2`, `7K2K_entity2`, `7K2O_entity2`) still carry the defect: each has
   a separate, unresolved proximity-bonding artifact (e.g. a spurious bond RDKit perceived between
   two atoms that aren't really bonded) that made a safe automatic fix impossible. Treat any other
   unusual bond order or formal charge on a non-standard residue as a possible perception artifact
   too, not necessarily deposited chemistry.
3. **No unbound (apo) receptor structures are shipped.** Only the bound complex is provided.
4. **This is a self/re-docking benchmark**: the receptor conformation shown is the one bound to
   this exact peptide. It is not a cross-docking or blind structure-prediction test.

## 5. Files and key columns

**Every column in `MANIFEST.csv` is defined in [`DATA_DICTIONARY.md`](DATA_DICTIONARY.md).**

```python
import pandas as pd
m = pd.read_csv("MANIFEST.csv")
paths = "entries/" + m.case_id + "/complex.pdb"
```
