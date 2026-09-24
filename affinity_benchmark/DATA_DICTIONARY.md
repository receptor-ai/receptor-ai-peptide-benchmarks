# Data dictionary

This file explains every column of every file in this folder.

## `structures/<series>/affinity.csv`

All 11 files use the same 3 columns.

| Column | Meaning | Values |
|---|---|---|
| `peptide` | The peptide id, as the source paper names it | e.g. `pDIQ`, `SYN1327`, `12` |
| `measure` | What the paper measured | `IC50`, `Kd` or `Ki` |
| `p_potency` | The potency on a negative log10 scale | 3.17 to 11.66 |

`p_potency = -log10(potency in mol/L)`.

A higher `p_potency` means a stronger binder. To get a concentration:

```python
nM = 10 ** (9 - p_potency)
```

`p_potency` carries 4 decimal places. The measurements behind it carry 2 to 3
significant figures. The 4th decimal is not meaningful on its own. It is there
so the conversion back to a concentration is clean.

The metric of this benchmark is a rank correlation, so only the order of
`p_potency` matters. The log10 step cannot change that order.

| series | rows | `measure` |
|---|---|---|
| mcl1 | 12 | IC50 |
| mdm2 | 11 | IC50 |
| mdm2_new | 27 | Kd |
| mdmx | 11 | IC50 |
| fcrn | 72 | IC50 |
| keap1_circ | 14 | IC50 |
| keap1_lin | 7 | IC50 |
| calmodulin | 17 | Kd |
| src | 8 | IC50 |
| upa | 5 | Ki |
| upa/noncanonical | 7 | Ki |

`mdm2` and `mdmx` hold the same 11 peptides. Each file gives the potency
against its own target.

## `structures/<series>/alignment.tsv`

| Column | Meaning | Values |
|---|---|---|
| `peptide` | The peptide id. It joins to `affinity.csv` | e.g. `pDIQ` |
| `1`, `2`, `3`, ... | One column per sequence position | a residue code, or `-` for a gap |

The mutation map of an analog is the difference from the reference peptide.
The reference peptide is in `ref_peptide.pdb`. A residue code can name a
non-canonical residue, for example `Pen`, `Sar`, `NMeLeu`, `Nle`, `Abu` or
`Hse`. `NH2` and `Ace` are terminal caps.

Some files carry extra columns. They are the source table as it was supplied:

| Column | Where | Meaning |
|---|---|---|
| `IC50 (nM)` | mcl1, keap1_circ, keap1_lin | The potency as the paper prints it |
| `IC50 (MDM2) nM`, `IC50 (MDMX) nM` | mdm2, mdmx | The potency against both targets |
| `IC50 (μM)` | fcrn | The potency as the paper prints it, WITH the qualifier |
| `Mismatches` | fcrn | The number of positions that differ from `SYN1327` |
| `source_row` | fcrn | The row index in the supplied source table. Gaps mean dropped rows |

Prefer `affinity.csv`. These columns repeat it, in mixed units.

## Censored measurements — read this before you rank

Do not read a censored measurement as an exact value. Some papers report a
bound, and `affinity.csv` stores that bound. A censored peptide sits at its
bound, and its true potency is weaker.

| series | censored | what the paper prints |
|---|---|---|
| fcrn | 14 of 72 | `>125` (10 rows), `>250` (3), `>500` (1) |
| mcl1 | 2 of 12 | `>3,600` and `>10,000` |

The fcrn column `IC50 (μM)` in `alignment.tsv` keeps the qualifier, so you can
find those 14 rows. For mcl1 the qualifier is only in the paper.

The effect on a rank correlation is small, because a censored peptide stays
among the weakest either way. Worst case, if every censored row of a series is
treated as a tie: 0.03 on the fcrn rho and 0.007 on the mcl1 rho.

The uPA series is different. 4 peptides do not bind at all, above 1000 μM.
They are removed, not censored. They are in no file.

## `structures/<series>/ref_protein.pdb` and `ref_peptide.pdb`

The receptor and the peptide are in 2 separate files. Both files use chain A.
Neither file holds a water or an additive.

Each file starts with 3 `REMARK 999` lines. They name the PDB entry the
structure comes from:

```
REMARK 999 SOURCE PDB ENTRY: 3JZS
REMARK 999 CONTENT: the receptor of that entry
REMARK 999 These are not the deposited coordinates. See README section 2.
```

| series | PDB entry |
|---|---|
| mcl1 | 6VBX |
| mdm2 | 3JZS |
| mdmx | 3JZQ |
| mdm2_new | 3LNZ |
| fcrn | 3M17 |
| keap1_circ | 7K2S |
| keap1_lin | 7K2A |
| calmodulin | 2BBM |
| src | 1SPS |
| upa | 4X1Q |

Every receptor matches the amino-acid sequence of its entry at 100 %.

The calmodulin receptor also holds the 4 Ca2+ ions of 2BBM. They are `ATOM`
records with the residue name `CA`, and they are part of the binding mode. No
other series holds an ion or a cofactor.

The coordinates are prepared, not deposited. They hold the hydrogens and a
CONECT record for every bond that a standard residue template does not give.
They are protonated at pH 7.4.

## Load a file

```python
import pandas as pd
aff = pd.read_csv("structures/mdm2/affinity.csv")
aln = pd.read_csv("structures/mdm2/alignment.tsv", sep="\t")
df = aln.merge(aff, on="peptide")
```
