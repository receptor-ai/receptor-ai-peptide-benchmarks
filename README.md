# Peptide Platform Benchmarks

Public benchmark datasets for Receptor.AI's peptide platform.

## Benchmarks

- [`bmi_benchmark/`](bmi_benchmark/) — binding-mode identification for peptides. Two
  curated datasets of peptide–protein X-ray complexes: BMI-200, a 202-entry diversity
  set; and BMI-MODES, 31 pockets that each hold 2 to 5 experimentally supported poses.
  Ring topology and secondary structure are measured from the coordinates. See its
  [README](bmi_benchmark/README.md).
- [`docking_benchmark/`](docking_benchmark/) — peptide–protein re-docking. docking985: 985 curated
  peptide–protein X-ray/NMR complexes drawn from seven published peptide-docking benchmarks, each
  shipped as a native structure plus a scoring library for evaluating predicted docking poses
  against it. See its [README](docking_benchmark/README.md).
- [`permeability_benchmark/`](permeability_benchmark/) — predict the passive membrane
  permeability of a cyclic peptide from one 3D structure. Scored per source paper against
  the curated CycPeptMPDB measurements: 31 PAMPA sources (7027 rows) and 10 Caco-2 sources
  (1310 rows). See its [README](permeability_benchmark/README.md).
