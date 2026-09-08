# BMI-200 data dictionary

This file explains every column in the BMI-200 tables. Two tables ship:

- `MANIFEST.csv` — one row per entry (202 rows), all metadata.
- `topology_audit.csv` — one row per core entry (200 rows), ring topology measured from the coordinates.

Load either table with `keep_default_na=False, na_values=['']` so that the element symbol `NA` (sodium)
is not read as a missing value.

A blank cell means "not applicable" or "not measured" for that entry (for example, a ring-closure column
is blank for a linear peptide). Counts are integers; distances are in ångström (Å); fractions are 0–1.
Every value is measured from the deposited coordinates or read from the RCSB PDB / UniProt, not copied
from any deposited secondary-structure or topology annotation.

---

## MANIFEST.csv

### Identity and selection

| Column | Meaning | Values |
|---|---|---|
| `pdb_id` | RCSB PDB accession of the source structure | 4-character code |
| `set_membership` | Which selection stream the entry came from | `core_200` (stratified diversity core) or `requested_addition` (two extra entries folded in) |
| `selection_rule` | How the entry was chosen | `stratified_maxmin` or `requested_addition` |
| `selected_for_ss_stratum` | Secondary-structure stratum this entry was picked to fill | `helix`, `sheet`, `turn_only`, `ppii`, `helix_310`, `none`, `helix_sheet` |
| `cluster_rank` | Rank of the entry inside its receptor sequence cluster (see `receptor_seqid30_cluster`); 1 is the representative | 1–8. Use `== 1` or `<= 2` to reduce receptor redundancy |
| `quality_relaxed` | 1 if the quality gate was relaxed for this entry | 0 for the core; blank for the two additions |
| `peptide_entity` | PDB entry id plus the peptide entity number | e.g. `2YQ7_2` |
| `peptide_subchain` | mmCIF `label_asym_id` of the kept peptide copy | chain id |
| `peptide_auth_chain` | Author (deposited) chain id of the peptide | chain id |
| `peptide_n_chains` | Number of chains that make up the peptide | usually 1 |

### Peptide sequence and chemistry

| Column | Meaning | Values |
|---|---|---|
| `peptide_seq_raw` | Designed peptide sequence. Standard residues are one-letter; every non-standard residue, cap, or linker is a 3-letter code in parentheses | e.g. `(ACE)WIAQELREIGDKFNAYYA(NH2)` |
| `components_in_peptide_file` | All component codes in `peptide.cif`, in backbone order (N→C), `;`-separated | e.g. `ACE;TRP;ILE;…` |
| `n_components_in_peptide_file` | Number of components in `peptide.cif` (residues plus caps/linkers) | integer |
| `n_res_observed` | Number of residues modelled in the coordinates. **Use this for any geometric measure.** | 5–20 |
| `has_noncanonical_component` | 1 if the peptide contains any non-standard component (non-canonical residue, cap, or linker) | 0 / 1 |
| `noncanonical_components` | `;`-list of all non-standard component codes, caps and linkers included | e.g. `ACE;MK8;NH2` |
| `peptide_ncaa` | `;`-list of non-canonical **amino-acid** codes only (caps and linkers excluded) | e.g. `0EH;MK8` |
| `caps_in_coordinates` | Terminal caps present in the coordinates | e.g. `ACE;NH2` |
| `peptide_caps` | Terminal caps declared for the peptide entity | e.g. `ACE;NH2` |
| `n_d_residues` | Number of D-amino acids | 0–15 |
| `n_l_residues` | Number of L-amino acids | integer |
| `helix_handedness` | Handedness of any helix present | `right`, `left` (D-residue signature), `none` |

### Ring topology and bond geometry

| Column | Meaning | Values |
|---|---|---|
| `cell` | Chemistry cell: measured topology crossed with modification | `cyclic_standard`, `cyclic_modified`, `linear_standard`, `linear_modified` (`standard` = unmodified) |
| `cell_label_frozen` | The cell label as first assigned, before the coordinate re-measurement | same four values |
| `topology_label_disagrees` | 1 if `cell` differs from `cell_label_frozen` | 0 / 1 |
| `cyclic_measured`, `is_cyclic` | 1 if the peptide is cyclic by measured ring topology | 0 / 1 |
| `cyclomatic_number` | Number of independent rings in the peptide bond graph (0 = linear) | 0–4 |
| `ring_closure_class` | Chemistry of the ring closure | `disulfide`, `thioether_or_S_C_bridge`, `head_to_tail_backbone_amide`, `amide_or_C_N_bridge`, `ester_or_C_O_bridge`, `carbon_carbon_bridge`, `linear` |
| `topology_measured` | Ring topology from coordinates | `macrocycle` / `linear` |
| `topology_observed_in_coordinates` | Whether a ring-closing bond was seen | `macrocyclic` / `no_ring_bond_observed` |
| `topology_labels` | Free-text topology descriptors | e.g. `adjacent_residue_nonstandard_linker` |
| `n_ring_closures` | Number of ring-closing bonds | 0–4 |
| `ring_closure_max_A` | Longest ring-closure bond length | Å (blank if none) |
| `n_ring_closures_geometry_ok`, `n_ring_closures_geometry_failed` | Ring closures whose geometry validated / failed | integers |
| `has_macrocyclic_ring_closure` | 1 if a backbone-spanning (macrocyclic) closure exists | 0 / 1 |
| `macrocycle_closure_span_residues` | Residues spanned by the macrocycle closure | integer (0 if none) |
| `head_to_tail_wrap` | 1 if the peptide is head-to-tail cyclic (N-terminus bonded to C-terminus) | 0 / 1 |
| `ring_closure_omega_deg` | Omega dihedral of the ring-closure amide | degrees (blank if not a head-to-tail amide) |
| `ring_closure_omega_call` | cis/trans call for that amide | `trans` (blank if N/A) |
| `ring_closure_omega_c_n_A` | C–N distance of that amide | Å |
| `ring_closure_omega_reason` | Status of the omega measurement | `measured` (blank if N/A) |
| `n_peptide_receptor_covalent_links` | Covalent bonds between peptide and receptor | 0–2 |
| `n_unannotated_intrapeptide_ring_bonds` | Intramolecular ring bonds inferred from geometry but absent from the deposited `_struct_conn` | 0 / 1 |
| `unannotated_ring_bond_detail` | The inferred bond: atom pair and distance | e.g. `GLU154.CD-LYS158.NZ@1.336A` |
| `n_inferred_links` | Bonds inferred by distance (not in the deposition) | 0 / 1 |
| `n_bond_length_outliers` | Bonds outside the expected length window | integer |
| `worst_bond_deviation_A` | Largest signed bond-length deviation from ideal | Å |
| `backbone_c_n_min_A`, `backbone_c_n_max_A` | Range of backbone amide C–N bond lengths | Å |
| `n_backbone_c_n_outside_1.20_1.50A` | Backbone C–N bonds outside 1.20–1.50 Å | integer (0 for all shipped entries) |

### Amide (omega) geometry

| Column | Meaning | Values |
|---|---|---|
| `omega_pattern` | Per-bond omega string over the peptide | `t` trans, `c` cis, `?` undefined |
| `omega_pattern_adjacent_bonded` | Same, restricted to covalently adjacent residue pairs | as above |
| `n_cis_amides`, `n_twisted_amides` | Count of cis / twisted (non-planar) amides, raw | integers |
| `n_cis_amides_corrected`, `n_twisted_amides_corrected` | Same, after correcting the ring-closure amide call. **Prefer the corrected columns.** | integers |
| `n_cis_amides_adjacent_bonded`, `n_trans_amides_adjacent_bonded`, `n_twisted_amides_adjacent_bonded` | Omega counts restricted to bonded adjacent pairs | integers |
| `n_omega_calls_dropped_not_bonded` | Omega calls skipped because the residue pair is not covalently bonded | integer |

### Chirality

| Column | Meaning | Values |
|---|---|---|
| `n_chirality_disagreements` | Residues whose geometric chirality disagrees with the chemical-component reference | 0–2 |
| `n_chirality_conflicts_genuine` | Genuine chirality conflicts | 0–2 |
| `n_chirality_notes_achiral_or_undetermined` | Notes where chirality is achiral or cannot be determined | 0–2 |
| `chirality_conflict_detail` | Text detail of any conflict | e.g. `ALA9:geom_D_ccd_L@34.42deg` |

### Secondary structure

Two representations are given for every secondary-structure measure. **bound** counts every backbone
hydrogen bond, including bonds to the receptor. **intrinsic** counts only the peptide's own backbone
hydrogen bonds, on the same coordinates. The intrinsic assignment is not a free-solution prediction; it
shows how much of the bound structure the peptide holds by itself.

The per-residue strings (`pep_ss_string_*`) use one character per residue, in backbone order:

```
H right-handed alpha-helix      h left-handed alpha-helix (D)     G right 3-10 helix    g left 3-10 helix
I pi-helix                      E extended, bridged to RECEPTOR   e extended, bridged to another PEPTIDE segment
x extended, no bridge partner   P polyproline II (L)              p polyproline II (D)
T beta-turn central residue     S gamma-turn central residue      C coil
X not assignable (backbone atom missing)     ? refused (backbone is not an alpha-amino-acid backbone)
```

Because `E` means "bridged to the receptor", it appears only in the bound string; the intrinsic string
keeps `e` (peptide-to-peptide bridge) and turns receptor-only bridges into `x`.

| Column | Meaning | Values |
|---|---|---|
| `ss_category_bound`, `ss_category_intrinsic` | Coarse secondary-structure category | `helix`, `helix_310`, `sheet`, `ppii`, `turn_only`, `none`, `helix_sheet` |
| `ss_category` | Display label of the bound category | `helix`, `sheet`, `PPII`, `turn_only`, `none/coil`, `helix+sheet` |
| `pep_ss_string_bound`, `pep_ss_string_intrinsic` | Per-residue string (alphabet above) | e.g. `XHHHHHHHHHHHHHHHHHCX` |
| `longest_helix_any_bound` | Longest run of any helix (H/h/G/g/I) | residues |
| `longest_helix_alpha_bound`, `longest_helix_310_bound` | Longest α- / 3₁₀-helix run | residues |
| `longest_strand_bound`, `longest_ppii_bound` | Longest β-strand / PPII run | residues |
| `beta_turn_types_bound` | β-turn types present, with counts | e.g. `I:2;II_p:2` |
| `gamma_turn_types_bound` | γ-turn types present, with counts | e.g. `gamma_inverse:1` |
| `dssp_strand_partner`, `strand_partner` | DSSP β-strand pairing partner | `none`, `self`, `receptor`, `both` |
| `sheet_survives_receptor_removal` | 1 if the β-sheet remains after receptor H-bonds are removed | 0 / 1 (blank if no sheet) |
| `n_residues_receptor_induced` | Residues whose secondary structure exists only because of receptor H-bonds | integer |
| `n_hbonds_peptide_internal` | Peptide-internal backbone H-bonds | integer |
| `n_hbonds_peptide_to_receptor`, `n_hbonds_receptor_to_peptide` | Backbone H-bonds peptide→receptor and receptor→peptide | integers |
| `n_internal_chain_breaks` | Gaps in the modelled peptide backbone | integer |
| `dssp_complex_string` | DSSP string for the peptide inside the full complex | DSSP letters + `X`/`-` |
| `frac_helix`, `frac_strand`, `frac_turn`, `frac_ppii_dssp`, `frac_coil` | DSSP composition fractions (denominator = assignable residues) | 0–1 |

### Receptor

| Column | Meaning | Values |
|---|---|---|
| `receptor_class` | Functional class | 21 classes, e.g. `adaptor_scaffold_ppi`, `ubiquitin_system`, `protease` |
| `receptor_subclass` | Finer functional label | free text |
| `receptor_uniprot` | UniProt accession | blank if none |
| `receptor_pfam` | Pfam family id(s) | `;`-separated |
| `receptor_desc`, `receptor_one_line` | Protein description (short / one-line) | free text |
| `receptor_organism` | Source organism | free text |
| `receptor_ec_number` | Enzyme Commission number | e.g. `3.4.21.89` (blank if not an enzyme) |
| `receptor_len_seqres` | Receptor length (SEQRES residues) | integer |
| `receptor_seqid30_cluster` | Receptor cluster id at 30 % sequence identity (redundancy grouping) | e.g. `2933_30` |
| `receptor_annotation_source` | Source of the receptor annotation | e.g. `uniprot:Q07817` |
| `receptor_annotation_confidence` | Confidence of that annotation | `high`, `medium`, `low` |
| `receptor_annotation_note` | Note when annotation is partial or absent | free text (blank if none) |
| `peptide_free_receptor_available` | 1 if a peptide-free (apo) structure of this receptor exists in the PDB | 0 / 1 |
| `cross_dockable` | 1 if a cross-docking test is possible (apo receptor exists). **The apo structure is not shipped; fetch it from the PDB.** | 0 / 1 |
| `n_receptor_chains_kept`, `n_receptor_chains_dropped` | Receptor chains retained / removed during cleaning | integers |
| `receptor_subchains_kept` | Which receptor subchains were kept | `;`-list |

### Structure quality and validation

| Column | Meaning | Values |
|---|---|---|
| `resolution_A` | Crystallographic resolution | Å (0.85–3.10) |
| `resolution_band` | Resolution bin | `<=1.5`, `1.5-2.0`, `>2.5`, … |
| `deposit_date` | RCSB deposition date | `YYYY-MM-DD` |
| `space_group` | Crystal space group | e.g. `P 21 21 21` |
| `n_models` | Number of models in the file | 1 (single X-ray model) |
| `dep_R`, `dep_Rfree` | Deposited R and R-free | 0–1 |
| `dep_clashscore` | Deposited wwPDB clashscore | number |
| `dep_rama_outliers_pct`, `dep_rota_outliers_pct`, `dep_rsrz_outliers_pct` | Deposited Ramachandran / rotamer / RSRZ outlier percentages | % |
| `dep_heavy_clash_pairs_selection`, `dep_heavy_clash_pairs_entry`, `dep_heavy_clash_pairs_peptide` | Deposited heavy-atom clash pairs in the kept selection / whole entry / peptide | integers |
| `rscc_mean`, `rscc_min` | Real-space correlation coefficient (map-to-model fit) of the peptide, mean / minimum | 0–1 |
| `rscc_n_below_0.8` | Peptide residues with RSCC below 0.8 | integer |
| `n_res_with_rscc` | Peptide residues with an RSCC value | integer |
| `n_res_with_rsrz`, `rsrz_max` | Peptide residues with an RSRZ value / the maximum RSRZ | integer / number |
| `receptor_rscc_mean` | Mean RSCC of the receptor | 0–1 |
| `b_peptide_mean`, `b_receptor_mean`, `b_ratio` | Mean B-factor of peptide / receptor, and their ratio | Å² / Å² / ratio |
| `occ_mean`, `occ_min`, `occ_n_below_1` | Peptide occupancy mean / minimum / number of atoms below 1.0 | 0–1 / 0–1 / integer |
| `sidechain_completeness` | Fraction of expected peptide side-chain atoms present | 0–1 |
| `n_missing_backbone_atoms`, `n_missing_sidechain_atoms` | Missing peptide backbone / side-chain atoms | integers (backbone is 0 for all shipped entries) |
| `n_res_seqres`, `n_res_unmodelled`, `peptide_unmodelled_fraction` | Designed residues / unmodelled residues / their fraction | integers / 0–1 |
| `peptide_len_seqres` | Designed (SEQRES) peptide length | integer |
| `peptide_length_band` | Length bin | `5-8`, `9-12`, `13-16`, `17-20` |
| `clean_core_pass` | 1 if the entry passed the "clean core" quality gate | 0 / 1 |
| `clean_core_n_failed_tests`, `clean_core_fail_reasons` | Number of failed clean-core tests / their reasons | integer / text |
| `scoreable_core` | 1 if the peptide backbone is complete and gap-free (the subset where a per-residue metric is meaningful) | 0 / 1 |

### Clashes and steric overlap

| Column | Meaning | Values |
|---|---|---|
| `clash_intra_peptide`, `clash_intra_peptide_per1k` | Intra-peptide clashes, count / per 1000 atoms | integer / rate |
| `clash_intra_peptide_same_residue` | Intra-peptide clashes within one residue | integer |
| `clash_peptide_receptor`, `clash_peptide_receptor_per1k` | Peptide–receptor clashes, count / per 1000 atoms | integer / rate |
| `clash_all`, `clash_all_per1k` | All clashes, count / per 1000 atoms | integer / rate |
| `clash_peptide_hetero` | Peptide–heteroatom clashes | integer |
| `clash_peptide_receptor_nucleophile_like` | Peptide–receptor clashes that look like a nucleophile contact | integer |
| `clash_intra_peptide_hbond_excluded`, `clash_peptide_receptor_hbond_excluded` | Same clash counts, with hydrogen-bond contacts excluded | integers |
| `worst_overlap_A` | Worst steric overlap in the entry | Å |
| `worst_overlap_involving_peptide_A` | Worst overlap that involves the peptide | Å |
| `worst_overlap_involving_peptide_reason` | Why that value is what it is | `measured`, `no_clash_involves_the_peptide`, `peptide_clash_present_but_not_in…` |
| `n_flags`, `flags` | Number of QC flags and the `;`-list of flag strings | integer / text |

### Assembly, copies, alternate locations

| Column | Meaning | Values |
|---|---|---|
| `n_copies_asu` | Peptide copies in the asymmetric unit | integer |
| `copy_rule_step` | Rule that chose the kept copy | `single_eligible_copy`, `lowest_mean_b`, `highest_mean_occupancy`, `most_modelled_residues`, `most_modelled_heavy_atoms` |
| `copy_bb_rmsd_min_A`, `copy_bb_rmsd_max_A`, `copy_bb_rmsd_mean_A` | Backbone RMSD among copies | Å (blank if one copy) |
| `copies_collapse_at_1.0A` | `True` if all copies are within 1.0 Å (effectively identical) | `True` / `False` |
| `n_altloc_peptide_residues` | Peptide residues with alternate locations | integer |
| `peptide_backbone_altloc` | `True` if the peptide backbone has alternate locations | `True` / `False` |
| `altloc_branch_max_backbone_rmsd_A` | Largest backbone RMSD between altloc branches | Å (blank if none) |
| `n_altloc_other_residues` | Non-peptide residues with alternate locations | integer |

### Heteroatoms, metals, heavy-atom counts

| Column | Meaning | Values |
|---|---|---|
| `n_hetero_kept`, `n_hetero_dropped` | Heteroatom components kept / removed during cleaning | integers |
| `hetero_kept` | Which heteroatom components were kept | `;`-list (blank if none) |
| `n_waters_dropped` | Water molecules removed | integer |
| `additive_in_peptide_site` | Crystallisation additives that sat in the peptide site (and were removed) | `;`-list (blank if none) |
| `n_heavy_peptide`, `n_heavy_receptor`, `n_heavy_hetero` | Heavy-atom counts for peptide / receptor / kept heteroatoms | integers |
| `metal_free` | 1 for every entry (the set is metal-free by construction) | 1 |
| `has_monatomic_halide` | 1 if a monatomic halide ion is kept | 0 / 1 |
| `n_spec_deviations`, `spec_deviation_ids` | Documented build-specification deviations for the entry (see `meta.json` → `spec_deviation`) | integer / id text |
| `has_superposed_alternative_components` | 1 if two components are modelled superposed | 0 for all shipped entries |

### Interface

| Column | Meaning | Values |
|---|---|---|
| `bsa_total_A2` | Buried surface area at the peptide–receptor interface | Å² |
| `peptide_buried_fraction` | Fraction of the peptide surface that is buried | 0–1 |
| `n_contact_pairs` | Peptide–receptor atom contact pairs | integer |
| `n_receptor_contact_residues` | Receptor residues contacting the peptide | integer |
| `interface_annotation_agrees` | `True` if the measured interface agrees with the deposited annotation | `True` / `False` |
| `possible_lattice_contact` | `True` if the interface may be a crystal-lattice contact | `False` for all shipped entries |

### Provenance and files

| Column | Meaning | Values |
|---|---|---|
| `built_utc` | UTC time the entry's files were generated | ISO-8601 timestamp |
| `has_strict_clean_record` | 1 if the entry carries a strict-hetero-cleaning record in `meta.json` | 0 / 1 |
| `complex_cif_bytes`, `peptide_cif_bytes`, `receptor_cif_bytes` | Byte sizes of the three coordinate files | integers |

---

## topology_audit.csv

One row per core entry (200 rows). It cross-checks the deposited cyclic/linear label against the ring
topology measured from the coordinates, and records the detailed ring closures. The two folded-in entries
are in `MANIFEST.csv` but not in this audit table.

| Column | Meaning | Values |
|---|---|---|
| `pdb_id` | RCSB PDB accession | 4-character code |
| `cell` | Chemistry cell (as in `MANIFEST.csv`) | `cyclic_standard`, `cyclic_modified`, `linear_standard`, `linear_modified` |
| `label_says_cyclic` | 1 if the assigned label calls the peptide cyclic | 0 / 1 |
| `measured_cyclic` | 1 if the coordinates measure the peptide as cyclic | 0 / 1 |
| `agrees` | 1 if `label_says_cyclic` equals `measured_cyclic` | 1 for every row |
| `cyclomatic_measured` | Independent rings measured from the coordinates | 0–4 |
| `cyclomatic_corpus` | Independent rings recorded in the entry metadata | 0–4 |
| `ring_closure_class_corpus` | Ring-closure chemistry recorded in the metadata | same classes as `ring_closure_class` |
| `topology_measured_corpus` | Topology recorded in the metadata | `macrocycle` / `linear` |
| `closure_types_measured` | Ring-closure chemistries found in the coordinates | `;`-list (blank if linear) |
| `closure_detail` | Each ring closure: chemistry, the two atoms, and the distance | e.g. `disulfide CYS2.SG-CYS7.SG 2.03A` |
| `n_residues_in_file` | Components in `peptide.cif` (residues plus caps/linkers) | integer |
| `n_amino_acid_residues` | Amino-acid residues only | integer |
| `n_components` | Number of molecular components in `peptide.cif` | 1 (a single connected molecule) |
| `n_heavy` | Heavy atoms in the peptide | integer |
| `n_res_observed_manifest` | `n_res_observed` copied from `MANIFEST.csv` for cross-checking | integer |
