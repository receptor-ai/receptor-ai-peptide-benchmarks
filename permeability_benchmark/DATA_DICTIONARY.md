# Data dictionary

This file explains every column of every CSV in this folder.

### `inputs/experimental.csv` — 7027 rows, one per peptide

| Column | Meaning | Values |
|---|---|---|
| `name` | The peptide id of CycPeptMPDB | e.g. `CycPeptMPDB_ID_1` |
| `assay` | The assay | always `pampa` |
| `logp_exp` | The measured log P_app, in log10 cm/s | -9.46 to -3.46 |
| `source` | The source paper id | e.g. `2006_Rezai_1` |
| `sequence` | The monomer sequence, as CycPeptMPDB gives it | e.g. `['dL', 'dL', 'L', 'dL', 'P', 'Y']` |
| `helm` | The structure in HELM notation | e.g. `PEPTIDE2{[dL].[dL].L.[dL].P.Y}$PEPTIDE2,PEPTIDE2,1:R1-6:R2$$$` |
| `smiles` | The structure as SMILES, with the stereochemistry | e.g. `CC(C)C[C@@H]1NC(=O)...` |

The 3 structure columns come from CycPeptMPDB without a change.

### `inputs/pampa_sources.csv` — 36 rows, one per listed source

| Column | Meaning | Values |
|---|---|---|
| `source` | The source paper id. It joins to `inputs/experimental.csv` | e.g. `2020_Townsend` |
| `doi` | The link to the paper | URL |
| `paper_title` | The title of the paper | free text |

`inputs/experimental.csv` holds 42 source papers. This file lists 36 of
them.

Load a file with:

```python
import pandas as pd
df = pd.read_csv("inputs/experimental.csv")
```
