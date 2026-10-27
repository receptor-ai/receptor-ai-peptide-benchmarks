# Peptide–protein structural benchmarks: BMI-200 and BMI-MODES

BMI-200 and BMI-MODES are two curated datasets of peptide–protein X-ray complexes from the RCSB PDB. Both
datasets use uniform coordinates and carry per-entry metadata. Both datasets cover chemistry that most
peptide datasets omit:

- macrocycles;
- non-canonical amino acids;
- D-amino acids;
- N-methylation;
- staples;
- thioether links;
- head-to-tail cyclisation;
- caps.

This project measures the ring topology and the secondary structure from the coordinates. This project
does not copy these annotations from the PDB deposition.

| Property | BMI-200 | BMI-MODES |
|---|---|---|
| Type | Diversity dataset. 202 distinct peptide–protein complexes. | Mode-discrimination dataset. 29 pockets. Each pocket holds 2 or more experimentally supported poses of one peptide. |
| Unit | One entry is one pocket with one pose. | One site is one pocket with 2 to 4 poses. |
| Selection rule | Distinct peptide sequences, spread over topology, size, and structure. | One peptide sequence, placed two or more ways in one pocket. |
| Cases | 202 entries | 29 sites |
| Poses | 202. One per entry. | 66 |
| Distinct receptors (UniProt) | 131 | 26 |
| Distinct peptide sequences | 202 | 29 |
| Resolution (Å) | 0.85–3.10. Median 1.78. | 1.10–2.50. Median 1.72. |
| Cyclic / linear (measured) | 102 / 100 | 1 / 28 |

The two datasets use opposite selection rules. BMI-200 requires a distinct peptide sequence in every
entry. BMI-MODES requires the same peptide sequence in two or more poses at one site. The two datasets
are therefore separate datasets, not one dataset with a flag. The two datasets are almost disjoint. They
share 6 receptor UniProt accessions and 3 PDB entries (`5OJR`, `7MX1`, `8IJ0`).

Read [`bmi200/README.md`](bmi200/README.md) and [`bmi_modes/README.md`](bmi_modes/README.md) for the
full description of each dataset and its file layout.

## Layout

```
bmi_benchmark/
  README.md      this file
  SHA256SUMS     a SHA-256 checksum for every delivered file
  scripts/       an example loader (load.py), the pep_ss secondary-structure tool (pep_ss.py),
                 and the builder and checker of the derived peptide files (build_derived.py, verify.py)
  bmi200/        the diversity dataset. See bmi200/README.md.
  bmi_modes/     the mode-discrimination dataset. See bmi_modes/README.md.
```

To load either dataset, run `python scripts/load.py`. The script needs only pandas. To verify file
integrity, run `sha256sum -c SHA256SUMS` from `bmi_benchmark/`.

## What both datasets share

**Source.** Every complex comes from one deposited X-ray structure in the RCSB PDB. This project cleaned
every structure the same way:

- remove hydrogens;
- keep one alternate-location branch per residue;
- keep one peptide copy and its contacting receptor chains;
- remove waters and crystallisation additives;
- handle metals and monatomic ions explicitly.

**Files.** Every structure has two layers:

- **Deposited.** `peptide.cif` and `receptor.cif` (`pose_NN.cif` and `receptor_NN.cif` in BMI-MODES)
  hold the deposited atoms left after the cleaning above, unchanged, in mmCIF. Nothing is added to them.
- **Derived.** `peptide_fixed.cif` adds the heavy atoms that the deposited model lacks (Derived files,
  below). `peptide_fixed.sdf` holds the same atoms in the same order, with bond orders, every formal
  charge 0 and hydrogens implicit. `peptide_fixed_protonated.sdf` is the same molecule protonated as at
  pH 7.4, with explicit hydrogens (Protonation, below). In BMI-MODES the three files are
  `pose_NN_fixed.cif`, `pose_NN_fixed.sdf` and `pose_NN_fixed_protonated.sdf`.
  `scripts/build_derived.py` makes all three from the deposited file and `meta.json`.

Together, the peptide and receptor files hold the cleaned complex atom for atom. The SDFs are derived
files only, because an SDF of an incomplete deposition reads as different chemistry: wherever the
deposition leaves out the C-terminal O, it ends in an aldehyde.

**Connectivity.** Each `peptide.cif`/`pose_NN.cif` carries a `_chem_comp_bond` block (intra-residue
bonds from the RCSB CCD) and a `_struct_conn` block (inter-residue links, assigned from geometry and
valence-checked). The derived SDFs take their bonds from these two blocks, with the CCD bond orders,
and their stereochemistry from the coordinates. The bonds were verified bond-for-bond against an
independent BMI-200 assignment.

**Completeness.** `meta.json` distinguishes two cases of a residue missing heavy atoms:

- `incomplete_residues` — genuine crystallographic disorder (25 structures; mostly Lys/Arg/Glu/Asp/Ser
  tips). The derived files carry the missing atoms, rebuilt from the CCD and relaxed by a restrained
  GAFF2 minimisation with every crystal atom and the receptor frozen, so the crystallographic pose is
  unchanged. `meta.json` lists them under `derived` with origin `rebuilt`. The one exception is 7Y8D,
  whose incomplete component (JFF) is a covalent staple crosslinker with a disordered pendant that
  cannot be reliably rebuilt; its derived files keep JFF as deposited, and `meta.json` lists it under
  `derived.not_rebuilt`. `derived.rebuild` records the method and the minimisation of each rebuild.
  The derived files of 4K8Y also model Pro13, which the deposition leaves out, so that the
  head-to-tail ring closes. Its two C–N bonds are 1.40–1.42 Å, because the deposited F12 C…D14 N gap
  is 3.91 Å, against 3.03–3.18 Å in the SFTI-1 structures 1SFI, 3P8F and 4K1E.
- `covalent_junctions` — 12 structures where a carboxylate O, terminal OXT, or Cys-S looks missing but
  is a covalent-bond position (ester, isopeptide lactam, backbone amide, or thioether staple); the
  residue is fully modelled. In each record `atom` is the vacated position and `bonded_to` is its
  partner; `distance_A` is how far that vacated position sits from the partner atom (0.03–0.5 Å) — small
  because the partner occupies the vacated site, which is what confirms it is the covalent-bond position.
  It is not the bond length (the bond itself is in `_struct_conn`).

**Derived files.** In 129 of the 268 peptides and poses, the derived files add heavy atoms to the
deposited ones. No deposited atom moves. The added atoms are the rebuilt side chains of Completeness
(152 atoms, 4K8Y Pro13 included) and terminal atoms that the molecule has but the deposited model
lacks:

- `OXT`, the second carboxylate O of a true C-terminus (63 termini);
- a declared cap that the deposition did not model: `ACE` (3) or `NH2` (1);
- a stand-in cap for an unmodelled neighbouring residue: `NME` on a truncated C-terminus (36) or `ACE`
  on a truncated N-terminus (42). The stand-in supplies the amide that the real chain has there. It
  does not model the missing residue.

Terminal atoms are placed at ideal geometry from the deposited atoms, with no minimisation. Over all
added atoms, the closest non-bonded heavy-atom contact (topological distance > 3, peptide and receptor)
is 2.21 Å. 7ABT carries a stand-in at both ends although its deposited sequence is complete: its
peptide is an eight-residue window of a PR repeat that continues through the crystal, and ARG8 C and
PRO1 N lie 1.34 Å from the PRO1 N and ARG8 C of symmetry mates. Five termini stay as deposited. In 4MNV
(ARG10, the end of chain B) and 6MPZ the C-terminal carbon is bonded to the receptor (Ser195 OG; Cys27
SG, where the component `GLZ` is an aldehyde in the CCD). At the 3BU6 N-terminus, the A_4X34 pose 2
C-terminus and the A_7TO8 pose 1 N-terminus, no stand-in fits farther than 2.2 Å from every other atom.

`_bmi_atom_origin` in the CIF and the `atom.prop.origin` list in the SDFs give each atom's origin:
`deposited`, `rebuilt`, `oxt`, `declared_cap` or `stand_in_cap`, plus `hydrogen` in the protonated SDF.
`atom.prop.label` names each atom as chain/component/number/atom. `meta.json` lists every added atom,
with its coordinates, under `derived` (`derived__pose_NN` in BMI-MODES). `terminal_completion` in
`MANIFEST.csv` and `poses.csv` summarises the ends; every other column describes the deposited files.
Score a predicted pose against the deposited atoms (`peptide.cif`, or origin `deposited`). A stand-in
is not part of the deposited molecule, and the orientation of a stand-in ACE is a guess.

**Protonation.** `peptide_fixed.sdf` keeps every atom at formal charge 0 with the CCD bond orders.
`peptide_fixed_protonated.sdf` (`pose_NN_fixed_protonated.sdf`) is the molecule at pH 7.4, by fixed
group rules:

- ionised: carboxylic acids (C-terminus, Asp, Glu), both OH of a phosphate (pSer, pThr and pTyr are
  dianions), sulfates, sulfonates and sulfinic acids;
- protonated: aliphatic amines, whose N is bonded only to sp3 carbons (N-termini, Lys, Orn, N-methyl
  and N-alkyl amines), and guanidines and amidines (Arg and its analogues, acylguanidines included);
- neutral: everything else. His keeps its H on NE2, as in the CCD. Cys, Tyr, amides, anilines,
  pyridines and azoles are neutral. The side-chain amine of Dap (`DNP`) is neutral too, because its
  pKa in a peptide is about 6.3.

The net charge is the SDF property `net_charge` and the column `net_charge_ph7.4` in `MANIFEST.csv`
and `poses.csv`. It runs from −7 to +10 over the 268 peptides and poses. The rules use the pKa values
of model compounds. A residue whose pKa the site shifts, such as a buried Asp or a His in a salt
bridge, keeps the default state. Dimorphite-DL 2.0.2 at pH 7.4 gives the same states, with two
exceptions. It protonates tertiary amide N (Pro, N-methyl residues), deprotonates phenols and makes His
an imidazolate, and none of these is the majority state at pH 7.4. It also protonates the Dap
side-chain amine. RDKit places the hydrogens from the heavy atoms, without optimisation. In 76 of the
268 files a hydrogen sits within 1.6 Å of a non-bonded atom, where two heavy atoms are close
(`verify.py` lists them).

**Peptide–receptor bonds.** `receptor_links` in `meta.json` (`receptor_links__pose_NN`) lists the
covalent and disulfide bonds that the PDB entry records between the peptide and the receptor, with the
symmetry operators of both atoms: 14 records in 11 entries and 2 poses. 5F8T's bond is to a symmetry
mate.

**Checks.** `scripts/verify.py` checks every structure:

- the derived files are exactly what `build_derived.py` makes;
- the deposited atoms are unchanged, and the SDFs match the CIF;
- the protonated SDF carries the charges of the pH 7.4 rules and only the hydrogens they imply;
- every residue is complete against its CCD definition;
- no added atom clashes.

With `--pdb`, it also checks every deposited atom against the PDB entry, the completion of each chain
end against the deposition and the crystal, and `receptor_links`. Both scripts need gemmi, RDKit and
numpy. The current release passes every check. Among the warnings:

- four deposited residues have a Cα with the handedness opposite to their CCD definition: 8VJP Ala9
  and Ala13, 7Z4S Cys13, 5FC2 Arg8;
- three structures overlap a copy of themselves in the crystal, through a special position or an
  alternate orientation: 3FIV, A_5CS2, A_6FNZ.

**Topology.** This project classifies each peptide as cyclic or linear from the ring topology of its
coordinates, not from the deposited annotation. For 28 of the 202 BMI-200 entries the measured topology
differs from the label assigned at selection (`cell_label_frozen`); `cell` carries the measured one,
and `topology_audit.csv` reads `agrees == 1` on every row.

**Secondary structure.** This project assigns secondary structure in two representations. Each dataset's
README gives its own secondary-structure tables. The `pep_ss` tool that makes these assignments is in
`scripts/pep_ss.py`. It needs gemmi, numpy, and scipy. Annotate one entry with
`python scripts/pep_ss.py --entry bmi200/entries/<PDB_ID>`.

- **bound**: the peptide in its bound conformation. This assignment counts every backbone hydrogen bond,
  including bonds to the receptor. This assignment describes the structure in the deposited coordinates.
- **isolated**: the same coordinates with the receptor hydrogen bonds removed. This assignment counts
  only the peptide's own backbone hydrogen bonds. The isolated assignment is not a prediction of the free
  peptide in solution. Do not read it as one. The isolated assignment shows only how much of the bound
  structure the peptide holds by itself.

**Metal content.** No metal sits in any peptide binding pocket, in either dataset. BMI-200 is metal-free
outright: any complex whose deposited structure contained a metal was excluded, and BMI-200 keeps and
flags the monatomic halides it retains (`has_monatomic_halide`, 22 of 202). BMI-MODES is pocket-metal-free
rather than metal-free: 7 of its 29 sites keep a deposited assembly metal in the receptor, outside the
binding site (`A_2AOF`, `A_2ZNE`, `A_3UA7`, `A_5MTW`, `A_5N8E`, `A_6GQN`, `B_P06873_3PTL`). The nearest
metal-to-peptide approach is 6.9 Å, so no metal contacts or coordinates any pose, and every pose records
`pose_env_metals` = `-`. BMI-MODES contains no monatomic halides and has no `has_monatomic_halide` column.

## Provenance and licensing

Coordinates and crystallographic metadata come from the RCSB PDB under CC0 1.0. Receptor functional
annotation comes from UniProt under CC BY 4.0. This project computes the secondary structure. Attribute
the RCSB PDB and UniProt when you use these datasets.

## Reproducibility of the numbers

This project computes every count and table in the per-dataset READMEs directly from the delivered files.
The `MANIFEST.csv` and `poses.csv` tables and the coordinates are the only sources. The BMI-MODES
coordinates round-trip: the separations, the mode counts, and the per-pose chemistry checks re-derive
from the files alone.

## How to cite

If you use BMI-200 or BMI-MODES, please cite this benchmark. Machine-readable
metadata is in [`CITATION.cff`](CITATION.cff). Plain text:

> Receptor.AI, Inc. (2026). BMI-200 and BMI-MODES: peptide–protein binding-mode
> identification benchmarks (Version 1.0.0). https://github.com/receptor-ai/receptor-ai-peptide-benchmarks
