# BMI-MODES: a mode-discrimination dataset

BMI-MODES holds 31 sites, 73 poses, and 28 receptors. A **site** is one receptor pocket that holds two or
more experimentally supported poses of one peptide. Each pose is the same molecule, in the same chemical
environment, placed a different way. This dataset does not test whether a method finds the correct pose.
This dataset tests whether a method produces **all** the poses, keeps them **apart** as distinct modes,
and **ranks** them.

BMI-200 cannot support this test, because BMI-200 holds one pose per entry. The two datasets share 3 PDB
entries (`5OJR`, `7MX1`, `8IJ0`).

```
bmi_modes/
  README.md              this file
  MANIFEST.csv           one row per site (31)
  poses.csv              one row per pose (73)
  ss_per_pose.csv        per-pose secondary structure (73)
  pose_attribution.csv   one row per pose pair (59)
  sites/<site_id>/       pose_NN.cif, receptor_NN.cif (one per pose), receptor.cif, meta.json
```

## 1. What the dataset holds

A site holds 2 to 5 poses. All 31 sites resolve to exactly two binding modes.

| Poses in the site | Sites |
|--:|--:|
| 2 | 25 |
| 3 | 2 |
| 4 | 3 |
| 5 | 1 |

Each pose comes from one of three kinds of evidence. This project never pools the three kinds.

| Pose source | Sites | What it is |
|---|--:|---|
| separate ASU copies only | 19 | two fully occupied molecules in equivalent pockets of one crystal |
| altloc branches only | 7 | one molecule refined as two partially occupied placements |
| both separate copies and altloc branches | 4 | a site that contains poses of both kinds |
| repeat deposition | 1 | the same peptide and receptor in two independent crystals |
| **total** | **31** | |

## 2. Modes

A mode is a cluster of poses under single linkage at 3.0 Å on the **core** separation. The core separation
is the pose difference with the terminal paired residue trimmed off each end. This project trims the
terminus because whole-pose RMSD cannot separate two binding modes from one mode with a disordered
terminus. Single linkage means two poses are distinct modes only when every cross-mode pair is 3.0 Å apart
or more.

All 31 sites hold exactly two modes. Separations run from 3.1 Å to 21.1 Å. 12 sites are 10 Å apart or
more. 23 of 31 sites still read as two modes at a 5 Å threshold. `MANIFEST.csv` also carries mode counts
at 2.0 Å and 5.0 Å, so the sensitivity is visible.

### The attribution label: placement or conformation

For every pose pair, `pose_attribution.csv` records the **in-place** RMSD (as the crystal placed the
poses) and the **superposed** RMSD (the same atoms after a Kabsch fit, with placement removed). A pair far
apart in place but close after superposition is one conformation in two places. A pair far apart both ways
needs two conformations.

| Pose pairs by source | Pairs | One conformation, two places (superposed < 1.0 Å) |
|---|--:|--:|
| altloc branches | 32 | 21 |
| separate copies | 27 | 4 |
| **all pairs** | **59** | **25** |

Altloc pairs are usually one conformation placed twice (21 of 32). Separate-copy pairs usually differ in
conformation as well (only 4 of 27 are placement-only). Each pose has its own receptor conformation, so
the pocket Cα RMSD between a mode and the reference pose is also available (median 0.40 Å; 4 poses of
1.0 Å or more are induced-fit).

## 3. Site selection

Every delivered site meets these conditions:

- every pose overlaps the reference pocket (closest approach ≤ 1.5 Å) and shares at least 5 contact
  residues;
- the site holds 2 or more distinct modes at 3.0 Å single linkage on the core separation;
- every pose agrees on covalent adducts, coordinating metals, and pocket components;
- the mode split survives trimming the termini;
- every counted pose has mean backbone occupancy ≥ 0.25;
- every pose contacts the same copies of the same receptor entities;
- two separate copies do not occupy the same space in the crystal frame.

One site, `A_6RYF` (an HIV V3-loop peptide), met these conditions but this project removed it after
review. Its central five residues are unmodelled, so each pose's backbone splits into two disconnected
fragments.

## 4. Secondary structure, per pose

This project recomputes secondary structure per pose from the delivered coordinates. The **bound**
representation counts all backbone hydrogen bonds, including bonds to the receptor. The **isolated**
representation counts the peptide's own hydrogen bonds only. The CSV columns name the isolated
representation `intrinsic`. Per pose is the right granularity, because two poses of one peptide can hold
different secondary structure.

| Category | bound | isolated |
|---|--:|--:|
| turn only | 22 | 27 |
| β (sheet) | 19 | 0 |
| extended, unpaired | 15 | 27 |
| helix | 9 | 9 |
| PPII | 7 | 10 |
| none / coil | 1 | 0 |
| **total poses** | **73** | **73** |

- 20 of 73 poses change category when this project removes the receptor. Not one of the 19 β poses keeps
  β in isolation. Hydrogen bonds to the receptor hold the peptide β-structure.
- 12 of the 31 sites hold poses that differ in bound secondary structure (for example `A_5OJR`: β and
  helix; `A_4X34`: helix and turn). Group `ss_per_pose.csv` by `site_id` to find these sites.

> **Use `ss_per_pose.csv` for secondary structure.** A per-entry assignment is not reliable for this
> dataset. For five entries, the per-residue alternate-location choice mixes two opposite-running branches
> into one chain. The delivered pose coordinates are correct, because each pose takes one explicit branch.

## 5. Files delivered per pose, and the peptide-content convention

Per site, `sites/<site_id>/` contains, for each pose _NN_:

- `pose_NN.cif` — the pose (the peptide), in the common reference frame. `pose_01` is the reference.
- `receptor_NN.cif` — that pose's own receptor: its assembly copy, its conformation, and the metals and
  cofactors in its pocket, in the same frame.
- `receptor.cif` — the cleaned receptor for the frame entry.
- `meta.json` — the measurements plus the placement of each pose.

> **Peptide-content convention.** A pose file is the peptide polymer subchain only. Some depositions
> attach a covalent non-polymer moiety to the peptide, such as a synthetic staple or linker (`ZOY` in
> `7MX1`) or a glycan (`NAG` on `7VFL` and `7EC3`). This project ships that moiety in the pose's
> `receptor_NN.cif`, not inside the pose. BMI-200 treats a covalently bonded linker as part of the
> peptide, so the two datasets use different conventions for the 3 shared entries. Account for this when
> you compare a shared entry across the two datasets.

## 6. Files, columns, and usage

Load the CSV files with `keep_default_na=False, na_values=['']`. Selected columns:

- `MANIFEST.csv` (per site): `site_id`, `uniprot`, `receptor_class`, `seq`, `n_poses_counted`,
  `n_distinct_modes`, `n_modes_at_2A`, `n_modes_at_5A`, `max_separation_A`, `difference_type`,
  `mode_sizes`, `n_in_bmi200`.
- `poses.csv` (per pose): `site_id`, `pose`, `label`, `pdb_id`, `mode`, `separation_from_ref_A`,
  `receptor_file`, `pocket_ca_rmsd_vs_pose1_A`.
- `ss_per_pose.csv` (per pose): `ss_string_bound`, `ss_string_intrinsic`, `ss_category_bound`,
  `ss_category_intrinsic`, `topology_measured`, `n_hbonds_peptide_internal`, `n_hbonds_peptide_to_receptor`.
- `pose_attribution.csv` (per pose pair): `inplace_bb_rmsd_A`, `superposed_bb_rmsd_A`,
  `core_ca_separation_A`, `placement_only`, `receptor_pocket_rmsd_A`.

```python
import pandas as pd
man   = pd.read_csv('MANIFEST.csv',    low_memory=False, keep_default_na=False, na_values=[''])
poses = pd.read_csv('poses.csv',       low_memory=False, keep_default_na=False, na_values=[''])
ss    = pd.read_csv('ss_per_pose.csv', low_memory=False, keep_default_na=False, na_values=[''])
big_sep = man[man.difference_type == 'large']   # 12 sites, modes 10 Å apart or more
```

## 7. What the dataset does not tell you

- **Which pose is more populated.** Nothing here ranks poses by solution population. Score recovery of all
  modes, not their order.
- **A cyclic-arm statement.** The dataset holds one cyclic site.
- **Fine-grained rates.** 31 sites support statements such as "mode 2 was missed in _k_ of 31", not a
  percentage quoted to a fine precision.
- **Generalisation to structured peptides.** Half the poses hold no backbone hydrogen bond of their own.
  The dataset over-represents receptor-templated peptides.
- **Whether the modes found are all the modes.** These are the poses that were modelled.
