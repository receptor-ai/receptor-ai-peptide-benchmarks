# docking985 data dictionary

This file explains every column in `MANIFEST.csv` (985 rows, one per entry). The same fields are
repeated in each entry's `meta.json`.

| Column | Meaning | Values |
|---|---|---|
| `case_id` | Unique entry id: `{pdb_id}_entity{peptide_entity_id}`. Matches the directory under `entries/` | e.g. `1B6J_entity2` |
| `pdb_id` | RCSB PDB accession of the source structure | 4-character code |
| `peptide_entity_id` | PDB entity number of the peptide chain(s), as deposited. Distinguishes independent peptide-binding sites within one PDB entry | integer |
| `benchmarks` | Published peptide-docking benchmark(s) this entry was drawn from, `;`-separated | `CSPSet`, `ProtPep37_2021`, `PepSet`, `PPDbench`, `LEADS-PEP`, `Vina_47`, `PIPER-FlexPepDock` |
| `uniprot_ids` | Receptor UniProt accession(s), `;`-separated | blank if unavailable |
| `experimental_method` | Deposited experimental method | `X-RAY DIFFRACTION`, `SOLUTION NMR` |
| `is_nmr` | 1 if `experimental_method` is an NMR method (the resolution filter does not apply) | `True` / `False` |
| `resolution_angstrom` | Crystallographic resolution | Å (blank for NMR) |
| `peptide_chain_id` | Author chain id of the single peptide copy shipped in `entries/<case_id>/peptide.pdb` | chain id |
| `n_peptide_chain_copies` | Number of crystallographic copies of this peptide in the deposited asymmetric unit. `peptide_chain_id` is the first of these; see docking985/README.md §3 for what this means for `receptor.pdb` when this is >1 | 1–28 |
| `peptide_length` | Number of residues in the shipped peptide chain (observed in the coordinates, not the designed/SEQRES length) | 2–20 |
| `peptide_sequence` | One-letter sequence of the shipped peptide chain. A non-standard residue reads as `X` only when it has no recognized canonical parent; otherwise it silently reads as the parent's letter (e.g. a D-amino acid or an N-methylated standard residue), which loses the modification, not just the identity. Use `noncanonical_components` for the real chemistry | e.g. `NXPIVX`; contrast `ALLVTAGLVLA` for 4YV9_entity2, which carries 9 non-canonical components and zero `X` |
| `noncanonical_components` | `;`-list of every non-standard component (residue, cap, or linker) in the shipped peptide, as its 3-letter CCD code, in backbone (N→C) order, one entry per occurrence (a repeated code appears more than once). Blank if the peptide is entirely standard | e.g. `DAL;MLE;MLE;MVA;BMT;ABA;SAR;MLE;MLE` |
| `n_protein_atoms` | Heavy-atom count of the shipped `receptor.pdb` | integer |
| `n_peptide_atoms` | Heavy-atom count of the shipped `peptide.pdb` (equivalently, `peptide.sdf` minus its explicit hydrogens) | integer |
| `dist_min_pocket` | Closest peptide-to-receptor heavy-atom distance, computed at selection time against the full deposited assembly | Å |
| `dist_max_pocket` | Farthest peptide-to-receptor heavy-atom distance, same computation | Å |
| `peptide_contact_fraction` | Fraction of peptide heavy atoms with a receptor heavy atom within 5 Å, same computation. The primary selection criterion (≥ 0.5 required) | 0–1 |
| `n_vdw_violations` | Count of peptide–receptor heavy-atom pairs closer than 0.7 × (sum of Van-der-Waals radii), same computation | integer |
| `vdw_violations_per_peptide_atom` | `n_vdw_violations` divided by the peptide's heavy-atom count at selection time (≤ 0.1 required) | ratio |

Load with:

```python
import pandas as pd
m = pd.read_csv("MANIFEST.csv")
```
