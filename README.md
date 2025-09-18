# README

This repository stores all necessary input files and script to reproduce the results of the [BindFlow paper](https://example.com).

## Installation

```bash
micromamba env create -f https://raw.githubusercontent.com/ale94mleon/protein-ligand-benchmark/refs/heads/main/environment.yml -y

micromamba activate plbenchmark
python -m pip install seaborn scipy openpyxl git+https://github.com/ale94mleon/protein-ligand-benchmark.git
```

## Description of the data

```.
.
├── analysis  ->  (all scripts to reproduce the figures and analysis of the paper)
│   ├── color-code-latex-table.ipynb
│   ├── make-full-data-Chen2023.ipynb
│   ├── make-full-data-single-ff.py
│   ├── make-matrix.py
│   ├── mol-img-maker.ipynb
│   ├── significance-test.ipynb
│   ├── stats-contrib-A2A-validation.ipynb
│   ├── stats-contrib-all-horizontal.ipynb
│   ├── stats-contrib-CyclophilinD-validation.ipynb
│   ├── stats-contrib-p38-tyk2-validation.ipynb
│   └── utility.py
├── BindFlow-inputs ->  (Inputs to run BindFlow)
│   ├── README.md
│   ├── A2A
│   │   ├── config-fep.yml
│   │   ├── config-mmpbsa.yml
│   │   ├── executor-fep.py
│   │   ├── executor-mmpbsa.py
│   │   ├── inputs
│   │   └── make_ndx.py
│   ├── CyclophilinD
│   │   ├── config-fep.yml
│   │   ├── config-mmpbsa.yml
│   │   ├── executor-fep.py
│   │   ├── executor-mmpbsa.py
│   │   └── inputs
│   ├── mcl1
│   │   ├── config-fep.yml
│   │   ├── config-mmpbsa.yml
│   │   ├── executor-fep.py
│   │   ├── executor-mmpbsa.py
│   │   ├── executor.py
│   │   └── inputs
│   ├── p38
│   │   ├── config-fep.yml
│   │   ├── config-mmpbsa.yml
│   │   ├── executor-fep.py
│   │   ├── executor-mmpbsa.py
│   │   ├── inputs
│   │   └── make_ndx.py
│   ├── ptp1b
│   │   ├── config-fep.yml
│   │   ├── config-mmpbsa.yml
│   │   ├── executor-fep.py
│   │   ├── executor-mmpbsa.py
│   │   ├── inputs
│   │   └── make_ndx.py
│   ├── SAMPL6-OA
│   │   ├── config-fep.yml
│   │   ├── config-mmpbsa.yml
│   │   ├── executor-fep.py
│   │   ├── executor-mmpbsa.py
│   │   └── inputs
│   ├── thrombin
│   │   ├── config-fep.yml
│   │   ├── config-mmpbsa.yml
│   │   ├── executor-fep.py
│   │   ├── executor-mmpbsa.py
│   │   ├── inputs
│   │   └── make_ndx.py
│   └── tyk2
│       ├── config-fep.yml
│       ├── config-mmpbsa.yml
│       ├── executor-fep.py
│       ├── executor-mmpbsa.py
│       ├── inputs
│       └── make_ndx.py
├── data
│   ├── experimental  ->  (Experimental data collected from references)
│   │   ├── A2A
│   │   ├── CyclophilinD
│   │   ├── mcl1
│   │   ├── p38
│   │   ├── ptp1b
│   │   ├── SAMPL6-OA
│   │   ├── targets.yml
│   │   ├── thrombin
│   │   └── tyk2
│   └── simulation
│       ├── bindflow
│       │   ├── gather  ->  (Data gathered from all sets)
│       │   ├── individual  ->  (Individual final CSVs files for each set)
│       └── external-validation  ->  (Other studies that also compute binding free energies on the same data set. Used for comparison with BindFlow)
└── README.md
```

## Gathering data

1. Run `data/simulation/bindflow/gather/make-full-data.ipynb` --> `data/simulation/bindflow/gather/BindFlow.csv`
2. Run `make-stats.ipynb` --> `data/simulation/bindflow/gather/BindFlow-stats.csv`

## Analysis

All figures of the BindFlow paper can be exactly reproduced with the Python scrips or IPython-notebooks from the `analysis` directory.

## BindFlow inputs

The full benchmark campaign can be reproduced from the inputs on `BindFlow-inputs`. You will need to first install BindFlow in your HPC. For that follow the instructions on the [BidnFlow docs](https://bindlfow.readthedocs.io/en/latest/).

The Python version used is the one reported in the [BindFlow paper](https://example.com).
