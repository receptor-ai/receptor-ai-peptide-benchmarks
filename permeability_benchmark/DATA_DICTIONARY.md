# Data dictionary

This file explains every column of every CSV in this folder.

### `inputs/experimental.csv` — 7027 rows, one per peptide

| Column | Meaning | Values |
|---|---|---|
| `name` | The peptide id of CycPeptMPDB | e.g. `CycPeptMPDB_ID_1` |
| `assay` | The assay | always `pampa` |
| `logp_exp` | The measured log P_app, in log10 cm/s | -9.46 to -3.46 |
| `source` | The source paper id | e.g. `2006_Rezai_1` |

### `inputs/pampa_source_membrane.csv` — 36 rows, one per annotated source

| Column | Meaning | Values |
|---|---|---|
| `source` | The source paper id. It joins to `inputs/experimental.csv` | e.g. `2020_Townsend` |
| `doi` | The link to the paper | URL |
| `paper_title` | The title of the paper | free text |
| `membrane_fine` | The membrane composition, as the paper states it | free text |
| `membrane_class` | The short membrane class | `lecithin/dodecane`, `Chugai biomimetic`, `DOPC/hexadecane`, `custom PC`, `hexadecane only`, `undocumented` |

`inputs/experimental.csv` holds 42 source papers. The 6 papers that this
file does not list have no membrane annotation.

Load a file with:

```python
import pandas as pd
df = pd.read_csv("inputs/experimental.csv")
```
