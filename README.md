# Receptor.AI Peptide Platform Benchmarks

Public benchmark datasets supporting the Receptor.AI Peptide Platform white paper. Each one
measures a different task in peptide design — docking, binding-mode identification, affinity
ranking, passive permeability — and each ships the curated experimental reference data that a
method is scored against, with the code needed to score it. None of them ship the predictions
of any particular method.

## Datasets

- [`docking_benchmark/`](docking_benchmark/) — **docking985**. 985 curated peptide–protein X-ray
  and NMR complexes. Every entry holds the native receptor and peptide coordinates and a
  bond-order-aware peptide SDF, with the pose-scoring library used for the RMSD, DockQ and CAPRI
  metrics. The dataset holds no predicted poses. See its
  [README](docking_benchmark/README.md).
- [`bmi_benchmark/`](bmi_benchmark/) — **BMI-200 and BMI-MODES**. 202 peptide–protein complexes
  with per-entry ring-topology, secondary-structure and crystallographic-quality metadata, and 31
  mode-discrimination sites that hold 73 experimentally supported poses. The receptor-class
  annotations are in `bmi200/MANIFEST.csv`. See its
  [README](bmi_benchmark/README.md).
- [`affinity_benchmark/`](affinity_benchmark/) — **peptide scoring-function affinity benchmark**.
  The reference receptor and peptide structures, the analog alignments and the measured potencies
  for 10 series and 191 peptides across 8 targets, with the full scoring procedure. The sampled
  conformers and the per-pose scores are not included. See its
  [README](affinity_benchmark/README.md).
- [`permeability_benchmark/`](permeability_benchmark/) — **peptide passive permeability**. Curated
  CycPeptMPDB PAMPA data: 7027 measurements over 42 source papers, each peptide carrying its monomer
  sequence, HELM and SMILES, with source-level quality-control results. Per-peptide predicted
  values are not included. See its
  [README](permeability_benchmark/README.md).

## Licensing

The datasets are under CC BY 4.0 — each benchmark folder carries its own `LICENSE`. The code
released with them is under [MIT](LICENSE): the docking985 pose-scoring library and
entry-selection filter, and the BMI-200 secondary-structure assignment tool.

Coordinates and crystallographic metadata in docking985, BMI-200 and BMI-MODES come from the RCSB
Protein Data Bank under CC0 1.0. Receptor functional annotation in BMI-200 and BMI-MODES comes
from UniProt under CC BY 4.0. The permeability measurements come from CycPeptMPDB (Li et al.,
*J. Chem. Inf. Model.* 2023). Attribute these sources when you use the datasets.

## What is not here

The platform components that these benchmarks were used to evaluate — QuorumMap, ForceFold,
ArtiDock-P, the MM/GBSA implementation and the permeability models — are proprietary software of
Receptor.AI, Inc., available under a commercial licence or a collaboration agreement. This
repository holds the public reference data, not the methods measured against it.

## Citation

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff); each dataset carries its own for
citing it on its own. To cite the collection:

> Receptor.AI, Inc. (2026). Receptor.AI Peptide Platform Benchmarks (Version 1.0.0).
> https://github.com/receptor-ai/receptor-ai-peptide-benchmarks
