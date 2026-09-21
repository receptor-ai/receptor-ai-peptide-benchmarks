# Data dictionary

This file explains every column of every CSV in this benchmark.

## The result files

### `01_headline_results.csv` — 8 rows

| Column | Meaning | Values |
|---|---|---|
| `endpoint` | The assay endpoint | `PAMPA`, `Caco-2 (A->B only)` |
| `pool` | The tier, or the union of two tiers | `GOOD (r>0.5)`, `MEDIUM (0.2-0.5)`, `GOOD+MEDIUM`, `FAIL (r<0.2)` |
| `sources` | The source papers in the pool | integer |
| `n` | The peptides in the pool | integer |
| `pearson_r` | The Pearson r over the pooled rows | -1 to +1 |
| `signed_r2` | `sign(r) * r^2` | -1 to +1 |

### `02_pampa_by_membrane_class.csv` — 6 rows

| Column | Meaning | Values |
|---|---|---|
| `merged_membrane_class` | The PAMPA membrane class, with the full composition in the name. It maps one-to-one to `membrane_class` of `inputs/pampa_source_membrane.csv` | e.g. `lecithin/dodecane (1-2.5%)` |
| `n_total` | The peptides in the class | integer |
| `pooled_pearson_r` | The Pearson r over the pooled rows of the class | -1 to +1 |
| `pooled_signed_r2` | `sign(r) * r^2` | -1 to +1 |

### `03_pampa_per_source.csv` — 31 rows, one per PAMPA source paper

| Column | Meaning | Values |
|---|---|---|
| `tier` | The tier of the source | `GOOD`, `MEDIUM`, `FAIL` |
| `source` | The source paper id of CycPeptMPDB. Joins to `inputs/experimental.csv` and `inputs/pampa_source_membrane.csv` | e.g. `2020_Townsend` |
| `n` | The peptides with both a prediction and an experimental value | integer, 5 or more |
| `pearson_r` | The Pearson r for this source | -1 to +1 |
| `signed_r2` | `sign(r) * r^2` | -1 to +1 |
| `note_what_the_library_varies` | The chemistry that the source library changes, or the reason for a FAIL | free text, can be blank |

### `04_caco2_per_source.csv` — 10 rows, one per Caco-2 source paper

The columns `source`, `n`, `pearson_r` and `signed_r2` hold the same
values as in `03_pampa_per_source.csv`. There is no `tier` column and no
note column. `source` joins to `inputs/caco2_source_direction.csv`.

## The input files

### `inputs/experimental.csv` — 8337 rows, one per peptide and assay

| Column | Meaning | Values |
|---|---|---|
| `name` | The peptide id of CycPeptMPDB | e.g. `CycPeptMPDB_ID_1` |
| `assay` | The assay | `pampa` (7027 rows), `caco2` (1310 rows) |
| `logp_exp` | The measured log P_app, in log10 cm/s | -9.46 to -3.46 |
| `source` | The source paper id | e.g. `2006_Rezai_1` |

### `inputs/pampa_source_membrane.csv` — 36 rows, one per scored PAMPA source

| Column | Meaning | Values |
|---|---|---|
| `source` | The source paper id | e.g. `2020_Townsend` |
| `doi` | The link to the paper | URL |
| `paper_title` | The title of the paper | free text |
| `membrane_fine` | The membrane composition, as the paper states it | free text |
| `membrane_class` | The short membrane class. It maps one-to-one to `merged_membrane_class` of `02_pampa_by_membrane_class.csv` | `lecithin/dodecane`, `Chugai biomimetic`, `DOPC/hexadecane`, `custom PC`, `hexadecane only`, `undocumented` |

`inputs/experimental.csv` holds 42 PAMPA sources. The 6 sources that no
table scores have no row in this file.

### `inputs/caco2_source_direction.csv` — 21 rows, one per Caco-2 source

| Column | Meaning | Values |
|---|---|---|
| `source` | The source paper id | e.g. `2023_Ohta` |
| `doi` | The link to the paper or the patent | URL |
| `direction` | The transport direction | `AB`, `AB+BA_avg`, `unspecified_outsourced`, `unusable`, `unknown` |
| `evidence` | The quotation that fixes the direction | free text |
| `decision` | Whether the source enters the Caco-2 scoring | `keep`, `drop_averaged`, `drop_unverified` |

One source, `2018_Ramalho`, has no row in `inputs/experimental.csv`.

Load a file with:

```python
import pandas as pd
df = pd.read_csv("inputs/experimental.csv")
```
