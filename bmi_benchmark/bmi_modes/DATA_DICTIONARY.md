# BMI-MODES data dictionary

This file explains every column in the BMI-MODES tables. Four tables ship:

- `MANIFEST.csv` — one row per **site** (30 rows).
- `poses.csv` — one row per **pose** (71 rows).
- `ss_per_pose.csv` — one row per pose (71 rows): secondary structure of that pose.
- `pose_attribution.csv` — one row per **pose pair** (58 rows): the pairwise comparison behind the mode call.

Load every table with `keep_default_na=False, na_values=['']`.

**The three units.** A **site** is one pocket on one receptor where the same peptide sequence is observed
in two or more different placements. A **pose** is one such observed placement (from a second copy in the
asymmetric unit, an alternate location, or a repeat deposition of the same peptide–receptor pair). A
**mode** is a group of poses that occupy the pocket the same way; two poses are in different modes when
their separation exceeds the clustering threshold. Every site in this set resolves into exactly two
distinct modes.

**Separation** is a backbone-RMSD-like distance between two poses in the shared receptor frame. The
**core** is the set of residues paired across the poses being compared; core separation is measured on
those residues only. Distances are in ångström (Å). A blank cell means "not applicable". Lists that run
over the site's poses are `;`-separated and `-` marks "none" for that pose.

**Site id format.** `<route>_<...>`. Route `A` (one deposition) gives `A_<PDB>` (for example `A_8ZVY`).
Route `B` (repeat deposition) gives `B_<UniProt>_<PDB>` (for example `B_P06873_3PTL`).

---

## MANIFEST.csv — one row per site

### Site identity

| Column | Meaning | Values |
|---|---|---|
| `site_id` | Site identifier (see format above) | e.g. `A_8ZVY`, `B_P06873_3PTL` |
| `site_type` | Where the poses come from | `one_deposition` (one PDB entry) or `repeat_deposition` (repeat depositions of the same pair) |
| `route` | Assembly route | `A_one_deposition` or `B_repeat_deposition` |
| `difference_type` | Qualitative magnitude of the difference between the two modes | `clear`, `substantial`, `large` |
| `redundancy_group` | Redundancy-group id (sites in one group would be redundant; each site here is its own group) | `RG001`… |
| `frame_pdb` | PDB entry used as the common coordinate frame | 4-character code |
| `n_in_bmi200` | 1 if `frame_pdb` is also a BMI-200 entry | 0 / 1 |
| `resolution_A` | Resolution of the frame structure | Å |
| `possible_lattice_contact` | `True` if the pocket may be a crystal-lattice contact | `True` / `False` |

### Receptor and peptide

| Column | Meaning | Values |
|---|---|---|
| `uniprot` | Receptor UniProt accession(s) | `;`-separated (`-` if none) |
| `receptor_one_line` | Receptor description | free text |
| `receptor_class` | Receptor functional class | 13 classes, e.g. `transcription_factor`, `protease` |
| `seq` | Peptide sequence (one-letter, with 3-letter tokens in parentheses for non-standard residues/caps) | e.g. `VSFN(FRD)PQITAA(NH2)` |
| `peptide_tokens` | The same sequence, tokenised | same form as `seq` |
| `cell` | Chemistry cell of the peptide | `linear_standard`, `linear_modified`, `cyclic_modified` |
| `max_token_diff` | Largest sequence-token difference between the site's poses | 0 (the poses must be the same sequence) |

### Pose and mode counts

| Column | Meaning | Values |
|---|---|---|
| `n_entries` | PDB entries contributing to the site | 1–2 |
| `n_poses_modelled` | Poses modelled at the site | 2–5 |
| `n_poses_in_pocket` | Poses that fall in the pocket | 2–5 |
| `n_poses_elsewhere` | Poses that landed away from the pocket | 0–1 |
| `n_poses_not_comparable` | Poses that could not be compared | 0 |
| `n_poses_counted` | Poses counted toward the mode determination | 2–5 |
| `n_distinct_modes` | Distinct binding modes at the site | 2 for every site (the defining property) |
| `n_modes_at_2A`, `n_modes_at_3A`, `n_modes_at_5A` | Modes resolved at a 2 / 3 / 5 Å separation threshold | 2 at 2–3 Å; at 5 Å some sites merge to 1 |
| `n_modes_wholepose_3A`, `n_modes_completelinkage_3A` | Modes at 3 Å by whole-pose clustering / complete-linkage clustering (a robustness cross-check) | 2 for every site |
| `mode_sizes` | Poses per mode | e.g. `3+2` (mode 1 has 3 poses, mode 2 has 2) |
| `mode_members` | Which poses belong to each mode | e.g. `C,FA,IA` and `FB,IB` |
| `n_pairs_measured` | Pose pairs measured | integer |
| `n_pairs_without_core` | Pose pairs with no shared core | 0 |
| `n_poses_symmetry_resolved` | Poses whose placement was resolved with a crystallographic symmetry operation | 0–2 |

### Separations

| Column | Meaning | Values |
|---|---|---|
| `max_separation_A`, `min_separation_A` | Largest / smallest pairwise pose separation | Å |
| `min_inter_mode_separation_A` | Smallest separation between poses of different modes | Å |
| `min_inter_mode_core_sep_A` | Same, measured on the core residues | Å |
| `core_separation_widest_pair_A` | Core separation of the widest-separated pair | Å |
| `n_core_residues` | Residues in the common core used for core separation | integer |
| `pose_pairing` | How poses were aligned before comparison | `residue_number` or `number_then_name` |
| `pairwise_separation_A` | Every pose-pair separation, as text | e.g. `D~C=21.05` |
| `pairwise_core_separation_A` | Every pose-pair core separation, with the core size | e.g. `D~C=19.67/8res` |
| `pose_sep_from_ref_A` | Each pose's separation from the reference pose | `;`-list, e.g. `0.0;21.05` |
| `min_pose_occupancy` | Smallest occupancy among the site's poses | 0–1 |

### Per-pose lists (one field per pose, `;`-separated across poses)

| Column | Meaning | Values |
|---|---|---|
| `pose_specs` | Each pose's source spec: subchain / altloc | e.g. `D/-;C/-` (`/-` = no altloc) |
| `pose_frames` | Frame alignment per pose: `-` for the reference, `X>Y` = receptor chain X superposed onto chain Y | e.g. `-;B>A` |
| `pose_modes` | Mode assigned to each pose | e.g. `1;2` |
| `pose_occupancy` | Occupancy of each pose | e.g. `1.0;1.0` |
| `pose_contact_entities` | Entity id(s) each pose contacts | e.g. `1;1` |
| `pose_env_adducts` | Covalent adducts near each pose (e.g. glycans) | e.g. `NAGx2;NAG` (`-` = none) |
| `pose_env_metals` | Metals near each pose | `-` = none |
| `pose_env_ions` | Ions near each pose | `-` = none |
| `pose_env_near` | Any component near each pose | codes (`-` = none) |
| `pose_env_additives` | Crystallisation additives near each pose | e.g. `SO4x2;SO4` (`-` = none) |
| `poses_elsewhere` | Description of any pose that landed away from the pocket | text (blank if none) |

---

## poses.csv — one row per pose

| Column | Meaning | Values |
|---|---|---|
| `site_id` | Site this pose belongs to | see `MANIFEST.csv` |
| `pose` | Pose number within the site | 1–5 |
| `label` | Pose label: `subchain/altloc`, or a PDB id for a repeat-deposition pose | e.g. `C/A`, `3PTL` |
| `pdb_id` | PDB entry the pose comes from | 4-character code |
| `source` | Provenance of the coordinates | e.g. `PDB 8ZVY subchain D`, `PDB 2AOJ subchain C altloc A`, `PDB 3PTL peptide` |
| `mode` | Mode assigned to this pose | 1 / 2 |
| `separation_from_ref_A` | Separation from the reference pose (pose 1) | Å |
| `chain_in_pdb` | Chain id given to the peptide in the shipped pose file | `P`, `Q`, `R`, `S`, `T` |
| `n_atoms`, `n_residues` | Heavy atoms / residues in the pose | integers |
| `occupancy` | Crystallographic occupancy of the pose | 0–1 |
| `env_adducts` | Covalent adducts near the pose | codes (`-` if none) |
| `env_near` | Any component near the pose | codes (`-` if none) |
| `env_additives` | Crystallisation additives near the pose | codes (`-` if none) |
| `receptor_file` | The receptor file that pairs with this pose (`receptor_NN.cif` in the site directory) | e.g. `receptor_01.cif` |
| `receptor_chains` | Receptor chains in that file | e.g. `A,B` |
| `receptor_residues` | Receptor residues in that file | integer |
| `receptor_sup_rmsd_A` | RMSD of the receptor superposition onto the common frame | Å |
| `pocket_ca_rmsd_vs_pose1_A` | Pocket Cα RMSD of this pose's receptor against pose 1's | Å |
| `pocket_n_res_paired` | Pocket residues paired in that comparison | integer (blank for pose 1) |

---

## ss_per_pose.csv — one row per pose

Secondary structure of each pose, from the `pep_ss` tool (`scripts/pep_ss.py`). The shipped CLI
reproduces the BMI-200 per-entry assignments (`--entry bmi200/entries/<id>`); the per-pose rows here were
produced by the same assignment engine run over the pose files. **bound** counts every
backbone hydrogen bond, including bonds to the receptor; **intrinsic** counts only the peptide's own
backbone hydrogen bonds. The per-residue string uses one character per residue, in backbone order:

```
H/h alpha-helix (right/left)   G/g 3-10 helix (right/left)   I pi-helix
E extended, bridged to RECEPTOR    e extended, bridged to another PEPTIDE segment   x extended, no partner
P/p polyproline II (L/D)   T beta-turn residue   S gamma-turn residue   C coil
X not assignable (backbone atom missing)   ? refused (not an alpha-amino-acid backbone)
```

| Column | Meaning | Values |
|---|---|---|
| `site_id`, `pose`, `label`, `pdb_id`, `mode` | Pose identity (as in `poses.csv`) | — |
| `pose_source` | Whether the pose is a second copy or an alternate location | `copies` / `altloc` |
| `n_residues` | Residues in the pose | integer |
| `n_assignable` | Residues that can be assigned (backbone complete and α-amino-acid) — the denominator for fractions | integer |
| `n_refused` | Residues refused (backbone is not an α-amino-acid backbone) | integer |
| `n_internal_chain_breaks` | Gaps in the modelled backbone | 0 for all poses |
| `backbone_order_ambiguous` | 1 if the backbone order could not be resolved | 0 for all poses |
| `topology_measured` | Ring topology of the pose | `linear` / `macrocycle` |
| `cyclomatic_number` | Independent rings | 0 / 1 |
| `ring_closure_class` | Ring-closure chemistry | `linear` or `head_to_tail_backbone_amide` |
| `head_to_tail_wrap` | 1 if head-to-tail cyclic | 0 / 1 |
| `ss_string_bound`, `ss_string_intrinsic` | Per-residue string (alphabet above) | e.g. `CCEEXxEEExCX` |
| `ss_category_bound` | Coarse category (bound) | `helix`, `sheet`, `ppii`, `turn_only`, `extended_unpaired`, `none_or_coil` |
| `ss_category_intrinsic` | Coarse category (intrinsic) | `helix`, `ppii`, `turn_only`, `extended_unpaired` |
| `n_strand_residues_bound`, `n_strand_residues_intrinsic` | β-strand residues (bound / intrinsic) | integers |
| `n_residues_receptor_induced` | Residues whose structure exists only because of receptor H-bonds | integer |
| `n_hbonds_peptide_internal`, `n_hbonds_peptide_to_receptor` | Backbone H-bonds inside the peptide / to the receptor | integers |
| `longest_helix_bound`, `longest_strand_bound`, `longest_ppii_bound` | Longest run of each type (bound) | residues |
| `helix_handedness` | Helix handedness | `right` / `none` |
| `n_d_residues` | D-amino acids in the pose | 0 for all poses |

---

## pose_attribution.csv — one row per pose pair

The pairwise comparisons that decide whether two poses are the same mode or different modes.

| Column | Meaning | Values |
|---|---|---|
| `site_id` | Site the pair belongs to | see `MANIFEST.csv` |
| `site_type` | `one_deposition` or `repeat_deposition` | — |
| `pose_source` | `copies` or `altloc` | — |
| `pose_a`, `pose_b` | The two poses compared | pose labels, e.g. `C/A` and `C/B` |
| `same_mode` | 1 if the pair is the same mode, 0 if different modes | 0 / 1 |
| `n_paired_bb_atoms` | Backbone atoms paired in the comparison | integer |
| `inplace_bb_rmsd_A` | Backbone RMSD as deposited (no superposition) — measures the placement difference | Å |
| `superposed_bb_rmsd_A` | Backbone RMSD after best-fit superposition — measures the conformation difference | Å |
| `core_ca_separation_A` | Cα separation on the shared core residues | Å |
| `n_core_residues` | Residues in the shared core | integer |
| `placement_only` | 1 if the two poses differ mainly in placement, not in conformation (low superposed RMSD, high in-place RMSD) | 0 / 1 |
| `receptor_pocket_rmsd_A` | Pocket RMSD between the two poses' receptors | Å |
| `uniprot` | Receptor UniProt accession(s) | `;`-separated |
| `seq` | Peptide sequence | as in `MANIFEST.csv` |
