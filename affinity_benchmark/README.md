# Peptide scoring-function affinity benchmark

Rank a peptide analog series by predicted affinity. Every scoring method
reads the same poses.

This release holds the procedure and the inputs. It holds no results.

---

## 1. Design

- 10 series over 8 target proteins.
- 10 conformers for each complex. 5 independent repeats.
- Two variants of each system: **full** and **cropped**.
  - **Full** — the system given as a whole.
  - **Cropped** — cut to 12 A around the peptide.
- Metric: Spearman rho against the measured potency.
- A more negative score means a stronger predicted binder. The correlation is
  sign-flipped, so a positive rho is a correct ranking.

## 2. Procedure

1. Minimise the whole complex. Nothing restrained.
2. This receptor now serves every peptide of the series. Build one receptor
   per repeat.
3. Mutate the residues in place on it. Make 10 copies.
4. Minimise each copy in two stages:
   - Stage 1: only the new sidechains move.
   - Stage 2: the peptide and a 5 A protein shell move. The shell takes every
     whole residue that has an atom within 5 A of the peptide.
   - Hold the rest with a harmonic restraint, k = 10000 kJ/mol/nm2.
5. Score.
6. Repeat the whole procedure 5 times.

Every scoring method reads the pose of step 4. The score step minimises the
complex again inside a 12 A free shell. Measured on 40 poses: this moves the
structure less than 0.05 A. Step 4 is the step that relaxes the receptor.

Per repeat, take the best conformer of each peptide. Best = the lowest
predicted binding energy. Compute Spearman rho across the series. Report the
mean rho over the 5 repeats. The error bar is the standard deviation over the
repeats.

### 2.1 Reference states

| column | definition |
|---|---|
| `no_min` | complex - protein - peptide, at the bound geometry |
| `pept_min` | complex - protein - peptide minimised alone |
| `all_min` | complex - protein minimised - peptide minimised |

`pept_min` is the headline column. `src` uses `no_min`.

## 3. The datasets

Every series comes from published work. The last column gives the paper that
reports the peptides and the measured potency. We took no potency value from
another place.

| series | target | analogs in `affinity.csv` | potency | source |
|---|---|---|---|---|
| mcl1 | MCL1 | 12 | IC50 | [10.1021/acs.jmedchem.1c00005](https://doi.org/10.1021/acs.jmedchem.1c00005) |
| mdm2 | MDM2 | 11 | IC50 | [10.1074/jbc.M109.073056](https://doi.org/10.1074/jbc.M109.073056) |
| mdm2_new | MDM2 | 27 | Kd | [10.1016/j.jmb.2010.03.005](https://doi.org/10.1016/j.jmb.2010.03.005) |
| mdmx | MDMX | 11 | IC50 | [10.1074/jbc.M109.073056](https://doi.org/10.1074/jbc.M109.073056) |
| fcrn | FcRn | 72 | IC50 | [10.1016/j.bmc.2008.05.004](https://doi.org/10.1016/j.bmc.2008.05.004) |
| keap1_circ | KEAP1, cyclic | 14 | IC50 | [10.1021/jacs.0c09799](https://doi.org/10.1021/jacs.0c09799) |
| keap1_lin | KEAP1, linear | 7 | IC50 | [10.1021/jacs.0c09799](https://doi.org/10.1021/jacs.0c09799) |
| calmodulin | calmodulin | 17 | Kd | [10.1006/jmbi.1996.0229](https://doi.org/10.1006/jmbi.1996.0229) |
| src | Src SH2 | 8 | IC50 | [10.1371/journal.pone.0011215](https://doi.org/10.1371/journal.pone.0011215) |
| upa | uPA | 5 + 7 non-canonical | Ki | [10.1371/journal.pone.0115872](https://doi.org/10.1371/journal.pone.0115872) |

Notes on the sources:

- `mdm2` and `mdmx` hold the same peptides. One paper reports both potency sets.
- `keap1_circ` and `keap1_lin` come from one paper. That paper numbers the
  peptides 1 to 23 across the two series.
- `fcrn` names the parent peptide SYN746. SYN746 comes from
  [10.1073/pnas.0708960105](https://doi.org/10.1073/pnas.0708960105).
- `upa` takes 4 of its 12 measured peptides from two later papers:
  [10.1016/j.jmb.2015.08.005](https://doi.org/10.1016/j.jmb.2015.08.005) and
  [10.1016/j.biocel.2015.02.016](https://doi.org/10.1016/j.biocel.2015.02.016).
  All 3 papers measure Ki with the same chromogenic substrate assay. The values
  agree on the shared peptides, so they pool.

## 4. One receptor per series

The procedure holds one receptor conformation for the whole series. It ranks
analogs of one binding mode. It does not find a new binding mode.

Calmodulin shows why this matters. The analogs of that series can bind
different conformations of the target. A peptide that needs another receptor
conformation gets a wrong score. Check the binding mode of a series before you
trust the rank order.

## 5. Files

```
README.md             this file
DATA_DICTIONARY.md    every column of every file
CITATION.cff          machine-readable citation metadata
LICENSE               the licence of the data
SHA256SUMS            a SHA-256 checksum for every delivered file
structures/<series>/
  ref_protein.pdb     the reference receptor, the input of step 1
  ref_peptide.pdb     the reference peptide
  alignment.tsv       the peptide series, one row per analog
  affinity.csv        the measured potency
structures/upa/noncanonical/
  alignment.tsv       the 7 non-canonical analogs
  affinity.csv        their measured potency
```

[`DATA_DICTIONARY.md`](DATA_DICTIONARY.md) explains every column. Read it
before you load a file.

All 11 `affinity.csv` files use the same 3 columns:

```
peptide,measure,p_potency
pDIQ,IC50,7.0000
```

`measure` is `IC50`, `Kd` or `Ki`. `p_potency` is `-log10(potency in mol/L)`.
A higher value means a stronger binder. To get a concentration, use
`nM = 10 ** (9 - p_potency)`.

WARNING: Do not read 16 of the 191 measurements as exact values. The paper
gives a bound for them, and `affinity.csv` stores that bound. The true potency
is weaker. This affects 14 fcrn rows and 2 mcl1 rows.
[`DATA_DICTIONARY.md`](DATA_DICTIONARY.md) names them.

`alignment.tsv` gives one column per sequence position. The mutation map of an
analog is the difference from the reference peptide. It joins to `affinity.csv`
on `peptide`.

Each PDB starts with 3 `REMARK 999` lines. They name the PDB entry that the
structure comes from. Each PDB holds the hydrogens and a CONECT record for
every bond that a standard residue template does not give. The structures are
protonated at pH 7.4. They are prepared, not deposited, coordinates.

The residue parameters of the non-canonical residues UPA, UPD and UPG are
not part of this release. Build your own parameters for them.

## 6. How to cite

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff). Its
`references` block repeats the 11 source papers of section 3. Plain text:

> Receptor.AI, Inc. (2026). Peptide scoring-function affinity benchmark
> (Version 2.0.0). https://github.com/receptor-ai/receptor-ai-peptide-benchmarks
