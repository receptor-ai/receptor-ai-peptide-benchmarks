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
| Type | Diversity dataset. 202 distinct peptide–protein complexes. | Mode-discrimination dataset. 30 pockets. Each pocket holds 2 or more experimentally supported poses of one peptide. |
| Unit | One entry is one pocket with one pose. | One site is one pocket with 2 to 5 poses. |
| Selection rule | Distinct peptide sequences, spread over topology, size, and structure. | One peptide sequence, placed two or more ways in one pocket. |
| Cases | 202 entries | 30 sites |
| Poses | 202. One per entry. | 71 |
| Distinct receptors (UniProt) | 131 | 27 |
| Distinct peptide sequences | 202 | 30 |
| Resolution (Å) | 0.85–3.10. Median 1.78. | 1.10–2.50. Median 1.76. |
| Cyclic / linear (measured) | 102 / 100 | 1 / 29 |

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
  scripts/       an example loader (load.py) and the pep_ss secondary-structure tool (pep_ss.py)
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

The coordinates are in mmCIF format. For every entry, `peptide.cif` plus `receptor.cif` equals
`complex.cif`, atom for atom.

**Connectivity.** Each `peptide.cif`/`pose_XX.cif` carries a `_chem_comp_bond` block (intra-residue
bonds from the RCSB CCD) and a `_struct_conn` block (inter-residue links, assigned from geometry and
valence-checked). A companion `.sdf` beside each holds the same heavy atoms with CCD bond orders and
crystal stereochemistry, for RDKit/OpenBabel. Verified bond-for-bond against an independent BMI-200
assignment.

**Completeness.** `meta.json` distinguishes two cases of a residue missing heavy atoms:

- `incomplete_residues` — genuine crystallographic disorder (25 structures; mostly Lys/Arg/Glu/Asp/Ser
  tips). Each has a companion `peptide_fixed.cif`/`.sdf` (`pose_XX_fixed.*` in BMI-MODES): the missing
  atoms are rebuilt from the CCD and relaxed by a restrained GAFF2 minimisation with every crystal atom
  and the receptor frozen, so the crystallographic pose is unchanged. `meta.json.fixed_geometry` lists
  the rebuilt atoms. The one exception is 7Y8D, whose incomplete component (JFF) is a covalent staple
  crosslinker with a disordered pendant that cannot be reliably rebuilt; it ships as deposited and is
  flagged `no_fixed_companion` in `meta.json`.
- `covalent_junctions` — 12 structures where a carboxylate O, terminal OXT, or Cys-S looks missing but
  is a covalent-bond position (ester, isopeptide lactam, backbone amide, or thioether staple); the
  residue is fully modelled. In each record `atom` is the vacated position and `bonded_to` is its
  partner; `distance_A` is how far that vacated position sits from the partner atom (0.03–0.5 Å) — small
  because the partner occupies the vacated site, which is what confirms it is the covalent-bond position.
  It is not the bond length (the bond itself is in `_struct_conn`).

**Topology.** This project classifies each peptide as cyclic or linear from the ring topology of its
coordinates, not from the deposited annotation. For 28 of the 202 BMI-200 entries the coordinate re-measurement overrode the initially assigned
topology label; the final labels all agree with the coordinates (`topology_audit.csv` reads
`agrees == 1` on every row).

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
rather than metal-free: 7 of its 30 sites keep a deposited assembly metal in the receptor, outside the
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
