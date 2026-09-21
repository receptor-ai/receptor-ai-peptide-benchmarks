# Cyclic-peptide PAMPA permeability data — CycPeptMPDB

This folder holds curated PAMPA measurements for cyclic peptides. One row
gives one peptide and one measured log P_app. A second file lists the
source paper of every measurement.

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
| `experimental.csv` | The curated measurement per peptide. 7027 rows over 42 source papers. |
| `pampa_sources.csv` | The source paper of each measurement. 36 rows. |

[`DATA_DICTIONARY.md`](DATA_DICTIONARY.md) explains every column of every
CSV.

## Provenance and licensing

The experimental values come from CycPeptMPDB (Li et al., *J. Chem. Inf.
Model.* 2023), the public download `CycPeptMPDB_Peptide_All.csv`.
Attribute CycPeptMPDB and the original source papers when you use this
data. `inputs/pampa_sources.csv` gives the DOI of every source paper.

The curation in `inputs/` is our work. It carries the CC-BY-4.0 licence
in [`LICENSE`](LICENSE).

## How to cite

If you use this data, please cite it. Machine-readable metadata is in
[`CITATION.cff`](CITATION.cff). Plain text:

> Receptor.AI, Inc. (2026). Cyclic-peptide PAMPA permeability data from
> CycPeptMPDB (Version 1.0.0).
> https://github.com/receptor-ai/receptor-ai-peptide-benchmarks

## Warnings

Do not read the BLOD rows as measurements. The assay could not measure
them. CycPeptMPDB stores the placeholder -10. The curation removes these
rows.

Do not trust a source paper with fewer than 20 peptides. A small series
gives a noisy correlation.

`inputs/experimental.csv` holds 42 source papers.
`inputs/pampa_sources.csv` lists 36 of them. The other 6 papers have no
row.

The assay conditions change the measured value. Compare two peptides of
one source paper. Do not pool two source papers without a check.
