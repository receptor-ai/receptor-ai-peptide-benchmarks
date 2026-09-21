# Peptide permeability benchmark — CycPeptMPDB

The model predicts passive permeability from one 3D peptide structure.
Every number below is the Pearson r
of the predicted `logp_raw` against the measured log P_app.

## Layout

```
permeability_benchmark/
  README.md                        this file
  DATA_DICTIONARY.md               every column of every CSV
  SHA256SUMS                       a SHA-256 checksum for every delivered file
  01_headline_results.csv          pooled r per endpoint and tier
  02_pampa_by_membrane_class.csv   pooled r per PAMPA membrane class
  03_pampa_per_source.csv          r per PAMPA source paper
  04_caco2_per_source.csv          r per Caco-2 source paper
  inputs/                          the curated experimental data
```

## 01_headline_results.csv — pooled r per endpoint and tier

| endpoint | pool | sources | n | pearson_r | signed_r2 |
|---|---|---:|---:|---:|---:|
| PAMPA | GOOD (r>0.5) | 12 | 307 | +0.666 | +0.443 |
| PAMPA | MEDIUM (0.2-0.5) | 6 | 4611 | +0.355 | +0.126 |
| PAMPA | GOOD+MEDIUM | 18 | 4918 | +0.432 | +0.187 |
| PAMPA | FAIL (r<0.2) | 13 | 1771 | -0.345 | -0.119 |
| Caco-2 (A->B only) | GOOD (r>0.5) | 2 | 58 | +0.703 | +0.494 |
| Caco-2 (A->B only) | MEDIUM (0.2-0.5) | 4 | 998 | +0.225 | +0.050 |
| Caco-2 (A->B only) | GOOD+MEDIUM | 6 | 1056 | +0.265 | +0.070 |
| Caco-2 (A->B only) | FAIL (r<0.2) | 4 | 111 | -0.016 | -0.000 |

Tiers: GOOD is r > 0.5. MEDIUM is 0.2 <= r <= 0.5. FAIL is r < 0.2.
A source appears only with 5 peptides or more.

## 02_pampa_by_membrane_class.csv — pooled r per membrane class

| merged_membrane_class | n_total | pooled_pearson_r | pooled_signed_r2 |
|---|---:|---:|---:|
| lecithin/dodecane (1-2.5%) | 5691 | +0.318 | +0.101 |
| Chugai biomimetic (10% lecithin + 0.5% cholesterol + bile acid) | 878 | -0.022 | -0.001 |
| DOPC + hexadecane trilayer (BD/Corning Gentest) | 57 | +0.708 | +0.502 |
| undocumented (CycPeptMPDB N.D.) | 38 | +0.405 | +0.164 |
| pure hexadecane | 26 | +0.553 | +0.306 |
| custom PC (Avanti 1.67% L-alpha-PC) | 10 | +0.437 | +0.191 |

The lecithin/dodecane class holds 95 % of the PAMPA rows. Its r of
+0.318 best shows the performance of the model on the main PAMPA membrane.

## 03_pampa_per_source.csv — PAMPA r per source paper

| tier | source | n | pearson_r | signed_r2 | what the library varies |
|---|---|---:|---:|---:|---|
| GOOD | 2015_Marelli | 10 | +0.950 | +0.902 | Enantiomeric cyclic-peptide pairs |
| GOOD | 2021_Lee | 5 | +0.850 | +0.723 | Cyclosporin-O derivatives (IHB/chameleonicity) |
| GOOD | 2020_Le Roux | 47 | +0.835 | +0.698 | Semipeptidic macrocycles, non-peptidic linkers |
| GOOD | 2022_Saunders | 11 | +0.829 | +0.688 | Backbone amides -> heterocycles |
| GOOD | 2022_Lee | 24 | +0.705 | +0.497 | Cyclosporin O + peptoid side chains |
| GOOD | 2015_Wang | 50 | +0.643 | +0.414 | Diverse cyclic-peptide natural products |
| GOOD | 2016_Hickey | 18 | +0.627 | +0.393 | Exocyclic amide macrocycles |
| GOOD | 2021_Golosov | 23 | +0.621 | +0.385 | Thioether cyclic-peptide scaffolds |
| GOOD | 2021_Wang | 24 | +0.606 | +0.367 | Cyclic decapeptide flexibility/lipophilicity scan |
| GOOD | 2016_Frost | 12 | +0.574 | +0.330 | Oxadiazole rings replacing amides |
| GOOD | 2022_Taechalertpaisarn | 52 | +0.564 | +0.318 | New side-chain-to-backbone H-bond scaffold |
| GOOD | 2015_Ahlbach | 31 | +0.544 | +0.296 | Diverse natural-product cyclic peptides |
| MEDIUM | 2020_Barlow | 26 | +0.476 | +0.227 | Prodrug masking of H-bond donors |
| MEDIUM | 2022_Tamura | 10 | +0.437 | +0.191 | Thiazoline ring-bridged macrocycles |
| MEDIUM | 2020_Townsend | 2881 | +0.330 | +0.109 | Combinatorial 6-/7-mer library: L/D stereochemistry, N-Me, peptoid, beta-AA substitutions on one fixed backbone |
| MEDIUM | 2021_Kelly | 1519 | +0.327 | +0.107 | Lariat scaffolds, position-scanning |
| MEDIUM | 2022_Bhardwaj | 133 | +0.292 | +0.085 | Computationally designed macrocycles |
| MEDIUM | 2021_Comeau | 42 | +0.275 | +0.075 | N-/C-methylation scan |
| FAIL | 2015_Bockus_2 | 17 | +0.160 | +0.026 | Small n |
| FAIL | 2020_Furukawa | 36 | +0.016 | +0.000 | |
| FAIL | 2018_Lee | 6 | +0.014 | +0.000 | Small n |
| FAIL | 2006_Rezai_2 | 11 | +0.008 | +0.000 | Small n |
| FAIL | 2013_CHUGAI | 878 | -0.022 | -0.001 | Non-bilayer biomimetic assay (10% lecithin, cholesterol, bile acid, pH 6.5). Outside the calibration range of the model. n is large, so this result is real |
| FAIL | 2016_Furukawa | 680 | -0.048 | -0.002 | Peptomer (Calpha->N) scaffold. The sign flip is systematic. n is large, so this result is real |
| FAIL | 2017_Pye | 20 | -0.113 | -0.013 | |
| FAIL | 2020_Hosono | 11 | -0.212 | -0.045 | Amide -> ester substitution |
| FAIL | 2018_Kaneda | 7 | -0.221 | -0.049 | Small n |
| FAIL | 2018_Naylor | 72 | -0.230 | -0.053 | PAMPA is the secondary assay here. The primary assay is RRCK |
| FAIL | 2006_Rezai_1 | 10 | -0.243 | -0.059 | Small n |
| FAIL | 2015_Bockus_1 | 15 | -0.266 | -0.071 | Small n |
| FAIL | 2019_Ono | 8 | -0.776 | -0.602 | Small n, strong anti-correlation |

Two large sources hold the PAMPA pool down: 2020_Townsend (n=2881) and
2021_Kelly (n=1519). Remove both and the pool falls to 518 rows, but r
rises to about +0.66.

## 04_caco2_per_source.csv — Caco-2 A->B r per source paper

| source | n | pearson_r | signed_r2 |
|---|---:|---:|---:|
| 2023_Ohta | 569 | +0.328 | +0.108 |
| 2018_CHUGAI | 374 | +0.225 | +0.051 |
| 2015_Wang | 50 | +0.711 | +0.506 |
| 2022_Bhardwaj | 40 | +0.126 | +0.016 |
| 2023_Tanada | 38 | +0.272 | +0.074 |
| 2024_Kage | 35 | -0.178 | -0.032 |
| 2015_Hewitt | 18 | +0.153 | +0.024 |
| 2018_Buckton | 18 | -0.860 | -0.740 |
| 2015_Bockus_2 | 17 | +0.322 | +0.104 |
| 2016_Furukawa | 8 | +0.854 | +0.730 |

## inputs/

| File | Content |
|---|---|
| `experimental.csv` | The curated experimental value per peptide. 8337 rows. |
| `pampa_source_membrane.csv` | The membrane composition per PAMPA source. 36 rows. |
| `caco2_source_direction.csv` | The transport direction per Caco-2 source. 21 rows. |

[`DATA_DICTIONARY.md`](DATA_DICTIONARY.md) explains every column of every
CSV.

## The metric definitions

Each table holds 3 numbers per source.

- `n` — the peptides with both a prediction and an experimental value.
- `pearson_r` — the Pearson r of `logp_raw` against the measured
  log P_app.
- `signed_r2` — `sign(r) * r^2`. The magnitude equals the r² of a
  per-source OLS refit. The sign shows the direction of the
  correlation.

Pearson r does not change under a linear remap of `logp_raw`. So r
measures the rank-order quality of the prediction. It is independent of
the empirical per-assay head.

## Reproduction

This folder holds the experimental data and the results. It does not
hold the predicted `logp_raw` per peptide. The per-peptide predictions
are not part of this release.

The scoring has no random step, so no seed is needed. The experimental
values come from the public download `CycPeptMPDB_Peptide_All.csv`
(19 MB), after the BLOD filter.

## Provenance and licensing

The experimental values come from CycPeptMPDB (Li et al., *J. Chem. Inf.
Model.* 2023), the public download `CycPeptMPDB_Peptide_All.csv`.
Attribute CycPeptMPDB and the original source papers when you use this
benchmark. `inputs/pampa_source_membrane.csv` and
`inputs/caco2_source_direction.csv` give the DOI of every source paper.

The curation in `inputs/` and the results in the 4 numbered CSVs are our
work. They carry the CC-BY-4.0 licence in [`LICENSE`](LICENSE).

## How to cite

If you use this benchmark, please cite it. Machine-readable metadata is
in [`CITATION.cff`](CITATION.cff). Plain text:

> Receptor.AI, Inc. (2026). Peptide passive-permeability benchmark on
> CycPeptMPDB (Version 1.0.0).
> https://github.com/receptor-ai/receptor-ai-peptide-benchmarks

## Warnings

Do not trust a source with fewer than 20 peptides. The r value is noisy.

Do not read the BLOD rows as measurements. The assay could not measure
them. CycPeptMPDB stores the placeholder -10. The filter removes these
rows before the scoring.

The Caco-2 conformers of Ohta, Tanada and Kage (658 peptides) come from
Balloon 1.8.4, not from RDKit. 642 of 658 runs succeeded.

Conformer sampling did not improve this result. An unpublished internal
test used CREST multi-conformer sampling. It gave no gain over the single
conformer. This release does not contain that test.

`inputs/experimental.csv` is a later curation snapshot than the 4 result
tables. Its row count differs from the `n` of the tables for 7 source and
assay pairs. The differences go in both directions. The largest is 15
rows (2023_Ohta, Caco-2). Use the `n` of the tables with the published r
values. Use `inputs/experimental.csv` as the current curated data.
