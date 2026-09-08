# BMI-200: a diversity dataset of peptide–protein complexes

BMI-200 holds 202 peptide–protein X-ray complexes in four chemistry cells: cyclic or linear, crossed with
modified or unmodified. Every complex is metal-free.

```
bmi200/
  README.md            this file
  DATA_DICTIONARY.md   definition of every column in the tables below
  MANIFEST.csv         one row per entry, all metadata
  topology_audit.csv   ring topology measured from coordinates, one row per entry
  entries/<PDB_ID>/    complex.cif, peptide.cif, receptor.cif, meta.json
```

Every entry directory is self-contained. For every entry, `peptide.cif` plus `receptor.cif` equals
`complex.cif`, atom for atom. Load the coordinates with **gemmi**; `complex.cif` is the primary file.

## 1. Composition

| Cell | n | Receptor clusters | Length (observed) | Resolution (Å) |
|---|---|---|---|---|
| cyclic + unmodified | 50 | 23 | 7–20 | 1.00–2.00 |
| cyclic + modified | 52 | 46 | 5–20 | 0.96–3.10 |
| linear + unmodified | 50 | 42 | 5–20 | 1.00–2.00 |
| linear + modified | 50 | 43 | 5–20 | 0.85–2.00 |
| **total** | **202** | 130 | 5–20 (median 12) | 0.85–3.10 (median 1.78) |

- The dataset holds 102 cyclic and 100 linear peptides. This project split them by the ring topology
  measured from the coordinates.
- `modified` means the peptide contains 1 or more non-canonical amino acids.
- The dataset covers 21 receptor functional classes and 7 ring-closure chemistries. The dataset holds 126
  distinct receptors (UniProt accessions) across 130 clusters at 30% sequence identity.

## 2. Selection criteria

Every complex is metal-free and has a distinct peptide sequence. Any complex whose deposited structure
contained a metal is excluded. Monatomic halides are kept and flagged (`has_monatomic_halide`, 22 of 202).
Resolution ranges from 0.85 Å to 3.10 Å (median 1.78 Å). Within each cell, the entries are chosen for
diversity: they spread over peptide length, burial, secondary structure, ring-closure chemistry, and
receptor class, and reach the extreme values rather than the average.

## 3. Peptide chemistry

| Feature | Entries | How this project detects it |
|---|---|---|
| 1 or more D-amino acids | 25 | `n_d_residues > 0` |
| 1 or more cis amides | 33 | `n_cis_amides_corrected > 0` |
| a terminal cap | 51 | `peptide_caps` |
| 1 or more non-canonical residues (sequence token) | 101 | a parenthesised token in `peptide_seq_raw` |
| 1 or more non-canonical components (includes caps and linkers) | 117 | `has_noncanonical_component` |

The last two rows differ by definition. 101 peptides carry a non-standard amino-acid token. 117 peptides
carry any non-standard component, where a cap or a linker also counts. Filter on
`has_noncanonical_component` for chemistry rather than amino-acid vocabulary. The most D-substituted
peptide contains 15 D-residues. Helices are right-handed in 72 entries and left-handed in 2 entries. The
left-handed helix is the D-residue signature.

### Ring-closure chemistry (102 cyclic entries)

| Chemistry | n |
|---|---|
| disulfide | 31 |
| thioether or S–C bridge | 26 |
| head-to-tail backbone amide | 21 |
| amide or C–N bridge | 16 |
| carbon–carbon bridge | 6 |
| ester or C–O bridge | 2 |

`topology_audit.csv` gives the ring topology measured from the coordinates. The measured cyclic/linear
call matches the `cell` column.

## 4. Secondary structure

This project assigns secondary structure in two representations. The **bound** assignment counts every
backbone hydrogen bond, including bonds to the receptor. The **isolated** assignment counts only the
peptide's own backbone hydrogen bonds, on the same coordinates. This document reports both. The isolated
assignment is not a free-solution prediction. The isolated assignment measures how much of the bound
structure the peptide holds by itself.

| Category | bound | isolated | bound: cyclic | bound: linear |
|---|---|---|---|---|
| turn only | 51 | 60 | 37 | 14 |
| β (sheet) | 49 | 23 | 26 | 23 |
| α-helix | 34 | 38 | 19 | 15 |
| none / coil | 33 | 43 | 12 | 21 |
| PPII | 20 | 27 | 4 | 16 |
| 3₁₀-helix | 11 | 11 | 4 | 7 |
| helix + sheet | 4 | 0 | 0 | 4 |
| **total** | **202** | **202** | **102** | **100** |

- The cyclic arm is turn-dominated: 37 of the 51 turn-only entries are cyclic. The linear arm holds the
  coil, PPII, and mixed cases.
- 34 of 202 entries change category when this project removes the receptor. The loss concentrates in β.
  In this dataset, hydrogen bonds to the receptor hold most of the peptide β-structure.

The category is a coarse label. `none` and `turn_only` are residual categories: an entry in them can
still hold a short helical turn or a short β-strand that falls below the run-length threshold. For short
or nascent structure, read the per-residue string `pep_ss_string_bound` rather than the category alone.

## 5. Notes before you quote the dataset

1. **Receptor redundancy is deliberate.** The 202 entries span 130 receptor clusters at 30% sequence
   identity, so a receptor recurs up to 12 times. `cluster_rank` ranks the entries within each cluster.
   `cluster_rank == 1` selects 153 entries. `cluster_rank ≤ 2` selects 178 entries. Use these subsets to
   reduce receptor redundancy. Entries that share a receptor are not independent.
2. **Cross-dockable subset (161 of 202).** For these entries, a peptide-free structure of the same
   receptor exists in the PDB (`cross_dockable == 1`), so a cross-docking test is possible. This dataset
   does not ship that unbound receptor; fetch it from the PDB yourself. The remaining 41 entries support
   self-docking only.
3. **Gap-free-backbone subset (85 of 202).** `scoreable_core == 1` marks entries with a complete,
   gap-free peptide backbone. This is the subset on which a per-residue metric is meaningful.
4. **`n_res_observed` counts the modelled residues.** Use `n_res_observed` for any geometric measure.
   `peptide_seq_raw` is the designed (SEQRES) sequence and can differ from the observed residues.
5. **Antibody, designed-binder (24), and MHC (6) complexes** are a different binding-mode problem from a
   pocket. Report or exclude them separately when that distinction matters (`receptor_class`).

## 6. Files and key columns

**Every column in `MANIFEST.csv` and `topology_audit.csv` is defined in
[`DATA_DICTIONARY.md`](DATA_DICTIONARY.md).** This section lists only the columns most people start with.

Load `MANIFEST.csv` with `keep_default_na=False, na_values=['']`. This keeps the element `NA` (sodium)
from becoming null. Selected columns:

| Column | Meaning |
|---|---|
| `pdb_id` | PDB entry id |
| `cell` | chemistry cell (measured topology plus modified/unmodified) |
| `cyclomatic_number`, `ring_closure_class` | measured ring topology |
| `n_res_observed`, `peptide_seq_raw`, `peptide_ncaa`, `peptide_caps` | peptide identity |
| `has_noncanonical_component`, `noncanonical_components` | chemistry |
| `resolution_A`, `receptor_class`, `receptor_uniprot`, `receptor_seqid30_cluster` | context |
| `metal_free`, `has_monatomic_halide` | ion content |
| `ss_category_bound`, `ss_category_intrinsic` | secondary structure (isolated = intrinsic) |
| `pep_ss_string_bound`, `pep_ss_string_intrinsic` | per-residue SS strings |
| `cross_dockable`, `scoreable_core`, `cluster_rank` | usable-subset flags |
| `n_d_residues`, `n_cis_amides_corrected`, `helix_handedness` | chemistry and geometry |

```python
import pandas as pd
s = pd.read_csv('MANIFEST.csv', low_memory=False, keep_default_na=False, na_values=[''])
helical = s[s.ss_category_bound == 'helix']
paths   = 'entries/' + s.pdb_id + '/complex.cif'
```
