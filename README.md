# Peptide Platform Benchmarks

Public benchmark datasets for Receptor.AI's peptide platform.

## Benchmarks

- [`affinity_benchmark/`](affinity_benchmark/) — rank a peptide analog series by predicted
  affinity. 10 series over 8 target proteins, each in a full and a cropped variant. This
  release holds the procedure and the input structures; it holds no results. See its
  [README](affinity_benchmark/README.md).
- [`bmi_benchmark/`](bmi_benchmark/) — binding-mode identification for peptides. Two
  curated datasets of peptide–protein X-ray complexes: BMI-200, a 202-entry diversity
  set; and BMI-MODES, 31 pockets that each hold 2 to 5 experimentally supported poses.
  Ring topology and secondary structure are measured from the coordinates. See its
  [README](bmi_benchmark/README.md).
- [`docking_benchmark/`](docking_benchmark/) — peptide–protein re-docking. docking985: 985 curated
  peptide–protein X-ray/NMR complexes drawn from seven published peptide-docking benchmarks, each
  shipped as a native structure plus a scoring library for evaluating predicted docking poses
  against it. See its [README](docking_benchmark/README.md).
- [`permeability_benchmark/`](permeability_benchmark/) — curated PAMPA permeability
  measurements for cyclic peptides, taken from CycPeptMPDB. 7027 measured log P_app values
  over 42 source papers, with the DOI of 36 of those papers. See its
  [README](permeability_benchmark/README.md).
