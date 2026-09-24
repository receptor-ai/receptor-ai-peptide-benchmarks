# Peptide permeability benchmark — CycPeptMPDB

A benchmark dataset for the passive membrane permeability of peptides: a
curated reference set of PAMPA measurements for cyclic peptides, drawn from
CycPeptMPDB, to evaluate a permeability prediction against. One row gives one
peptide, its structure and one measured log P_app. A second file lists the
source paper of every measurement.

Every peptide carries its structure as SMILES, as HELM and as a monomer
sequence. No external download is necessary.

## Layout

```
permeability_benchmark/
  README.md                        this file
  DATA_DICTIONARY.md               every column of every CSV
  CITATION.cff                     machine-readable citation metadata
  LICENSE                          the licence of the curation
  SHA256SUMS                       a SHA-256 checksum for every delivered file
  inputs/                          the curated experimental data
```

## inputs/

| File | Content |
|---|---|
| `experimental.csv` | The structure and the curated measurement per peptide. 7027 rows over 42 source papers. |
| `pampa_sources.csv` | The source paper of each measurement, with its DOI. 42 rows. |

[`DATA_DICTIONARY.md`](DATA_DICTIONARY.md) explains every column of every
CSV.

## The pipeline

This data measures one workflow:

1. Take the peptide.
2. Generate a structure ensemble of it.
3. Score the ensemble with our software.
4. Compare the score with the measured log P_app of
   `inputs/experimental.csv`, as Pearson r or Spearman ρ.

## The source papers by tier

The tier groups a source paper by how well the workflow ranked its
peptides. GOOD is the best group, then MEDIUM, then FAIL.

| tier | source | n | what the library varies |
|---|---|---:|---|
| GOOD | 2015_Marelli | 10 | Enantiomeric cyclic-peptide pairs |
| GOOD | 2021_Lee | 5 | Cyclosporin-O derivatives (IHB/chameleonicity) |
| GOOD | 2020_Le Roux | 47 | Semipeptidic macrocycles, non-peptidic linkers |
| GOOD | 2022_Saunders | 11 | Backbone amides -> heterocycles |
| GOOD | 2022_Lee | 24 | Cyclosporin O + peptoid side chains |
| GOOD | 2015_Wang | 50 | Diverse cyclic-peptide natural products |
| GOOD | 2016_Hickey | 18 | Exocyclic amide macrocycles |
| GOOD | 2021_Golosov | 23 | Thioether cyclic-peptide scaffolds |
| GOOD | 2021_Wang | 24 | Cyclic decapeptide flexibility/lipophilicity scan |
| GOOD | 2016_Frost | 12 | Oxadiazole rings replacing amides |
| GOOD | 2022_Taechalertpaisarn | 52 | New side-chain-to-backbone H-bond scaffold |
| GOOD | 2015_Ahlbach | 31 | Diverse natural-product cyclic peptides |
| MEDIUM | 2020_Barlow | 26 | Prodrug masking of H-bond donors |
| MEDIUM | 2022_Tamura | 10 | Thiazoline ring-bridged macrocycles |
| MEDIUM | 2020_Townsend | 2881 | Combinatorial 6-/7-mer library: L/D stereochemistry, N-Me, peptoid, beta-AA substitutions on one fixed backbone |
| MEDIUM | 2021_Kelly | 1519 | Lariat scaffolds, position-scanning |
| MEDIUM | 2022_Bhardwaj | 133 | Computationally designed macrocycles |
| MEDIUM | 2021_Comeau | 42 | N-/C-methylation scan |
| FAIL | 2015_Bockus_2 | 17 | Small n |
| FAIL | 2020_Furukawa | 36 |  |
| FAIL | 2018_Lee | 6 | Small n |
| FAIL | 2006_Rezai_2 | 11 | Small n |
| FAIL | 2013_CHUGAI | 878 | A different assay format, at pH 6.5. Large n |
| FAIL | 2016_Furukawa | 680 | Peptomer (Calpha->N) scaffold. Large n |
| FAIL | 2017_Pye | 20 |  |
| FAIL | 2020_Hosono | 11 | Amide -> ester substitution |
| FAIL | 2018_Kaneda | 7 | Small n |
| FAIL | 2018_Naylor | 72 | PAMPA is the secondary assay here. The primary assay is RRCK |
| FAIL | 2006_Rezai_1 | 10 | Small n |
| FAIL | 2015_Bockus_1 | 15 | Small n |
| FAIL | 2019_Ono | 8 | Small n |

The table holds 31 of the 42 source papers of `inputs/experimental.csv`;
the remaining 11 are untiered.

`n` is the peptide count at the time of the work. `inputs/experimental.csv`
is a later curation snapshot, so a few counts differ. The largest
difference is 4 rows.

## Provenance and licensing

The experimental values come from CycPeptMPDB (Li et al., *J. Chem. Inf.
Model.* 2023, 63:2240, [10.1021/acs.jcim.2c01573](https://doi.org/10.1021/acs.jcim.2c01573)),
the public download `CycPeptMPDB_Peptide_All.csv`.
Attribute CycPeptMPDB and the original source papers when you use this
data. `inputs/pampa_sources.csv` gives the DOI of all 42 source papers.

The curation in `inputs/` is our work. It carries the CC-BY-4.0 licence
in [`LICENSE`](LICENSE).

## How to cite

If you use this benchmark, please cite it. Machine-readable metadata is in
[`CITATION.cff`](CITATION.cff). Plain text:

> Receptor.AI, Inc. (2026). Peptide passive-permeability benchmark on
> CycPeptMPDB (Version 1.0.0).
> https://github.com/receptor-ai/receptor-ai-peptide-benchmarks

## Warnings

Do not read the BLOD rows as measurements. The assay could not measure
them. CycPeptMPDB stores the placeholder -10. The curation removes these
rows.

The 7027 rows hold 6960 distinct structures. A peptide that 2 source
papers measure keeps one row per paper.
