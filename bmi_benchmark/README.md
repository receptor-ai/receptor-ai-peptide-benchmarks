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
| Type | Diversity dataset. 202 distinct peptide–protein complexes. | Mode-discrimination dataset. 31 pockets. Each pocket holds 2 or more experimentally supported poses of one peptide. |
| Unit | One entry is one pocket with one pose. | One site is one pocket with 2 to 5 poses. |
| Selection rule | Distinct peptide sequences, spread over topology, size, and structure. | One peptide sequence, placed two or more ways in one pocket. |
| Cases | 202 entries | 31 sites |
| Poses | 202. One per entry. | 73 |
| Distinct receptors (UniProt) | 126 | 28 |
| Distinct peptide sequences | 202 | 31 |
| Resolution (Å) | 0.85–3.10. Median 1.78. | 1.10–2.50. Median 1.80. |
| Cyclic / linear (measured) | 102 / 100 | 1 / 30 |

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

To load either dataset, run `python scripts/load.py`. The script needs only pandas.

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

**Topology.** This project classifies each peptide as cyclic or linear from the ring topology of its
coordinates, not from the deposited annotation. The deposited annotation disagrees with the coordinates
for 28 of the 202 BMI-200 entries.

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

**Metal content.** Every complex is metal-free. Any complex whose deposited structure contained a metal
is excluded. Each dataset keeps monatomic halides and flags them (`has_monatomic_halide`).

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
> identification benchmarks (Version 1.0.0). https://github.com/receptor-ai/peptide-benchmarks
