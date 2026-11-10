# Combined thermodynamics and 13C metabolic flux analysis in yeast

Code and analyses for  combined thermodynamics and 13C metabolic flux analysis (13C-MFA) approach. The approach is applied to *Saccharomyces cerevisiae*, comparing the wild type (**WT**) with the reduced-glucose-uptake mutant strain **TM6**.

The workflow has four steps:

1. **Define the flux space.** A thermodynamic–stoichiometric model (TSM) is constrained with physiology (exchange fluxes) and intracellular metabolomics. Variability analysis gives bounds on fluxes, metabolite concentrations and reaction Gibbs energies.
2. **Sample the flux space.** Thermodynamically feasible flux distributions are sampled (PTA), and 1000 of them are drawn as starting points for 13C fitting.
3. **Fit the 13C labeling data.** Each sample is fitted to the labeling data with FreeFlux. The best-fitting samples give the final (bootstrapped) flux configurations.
4. **Interpret and validate.** Flux configurations, reaction directions (compared with Yeast9), prediction of independent isotopomer data (Shen et al. 2024), and metabolite concentrations of the best-fitting samples.

---

## Repository structure

```
.
├── 1_data/
│   ├── a_raw_data/           # external datasets as received (e.g. Shen 2024 isotopomer data)
│   ├── b_processed_files/    # model files, metabolomics, physiology and regression data
│   └── c_synthetic_data/     # derived inputs: ΔrG° mean and covariance, drawn flux samples,
│                             #   pathway mapping, Yeast9 model
├── 2_notebooks/              # analysis notebooks and plotting scripts (see 2_notebooks/README.md)
├── 3_src/                    # project code imported by the notebooks
│   ├── analysis/             # config modules, I/O and plotting helpers, distribution comparison
│   ├── helper_thermoflux/    # yeast model construction and data processing for thermo_flux
│   ├── HELPER_FREEFLUX_NF/   # wrappers around FreeFlux for the 13C fitting problems
│   └── freeflux/             # FreeFlux source used for 13C fitting
└── 4_results/
    ├── a_fluxspace/          # variability analysis bounds, reaction directions
    ├── b_sampling
    ├── c_freeflux_run	/     # flux samples per strain
    ├── d_predict_isotopes/   # simulated isotopomer distributions for the Shen 2024 comparison
    └── e_metconc_sampling    #metabolite concentration sampling results
```

### `1_data/`

Main inputs:

| File | Content |
|---|---|
| `b_processed_files/yeast_data_xiang_metabolites.xlsx` | Intracellular metabolomics for WT and TM6 |
| `b_processed_files/extended_regressiondata_Canelas_Xiang_yeast.xlsx` | Exchange-flux data used as physiology constraints |
| `b_processed_files/cobra_roundp3_ymodel.xml` | Stoichiometric model (SM) in SBML |
| `b_processed_files/yeast_v3.211_starting_curated_withprotons.xlsx` | Curated model spreadsheet (reactions, subsystems, EC numbers) |
| `c_synthetic_data/yeast-GEM9.xml` | Yeast9 genome-scale model, used to compare reaction directions |
| `a_raw_data/13C_isotopomerdata_Shen2024/` | Isotopomer data from 1,2-13C glucose labeling (Shen et al. 2024), used for validation |

### `2_notebooks/`

Jupyter notebooks and two Python scripts, numbered by workflow step (`1*` flux space, `2*` sampling, `3*` 13C fitting results, `4*` interpretation). Each notebook and the figure it produces are described in [`2_notebooks/README.md`](2_notebooks/README.md).

### `3_src/`

Code shared by the notebooks. The notebooks add `../3_src/` to the Python path and mostly access data through two configuration modules:

- `analysis/config_pta.py`: paths and loaders for the flux space and sampling results.
- `analysis/config_13c.py`: paths and loaders for the 13C fitting results (fitted fluxes, objective values, bootstrapped flux sets).

Other modules provide plotting configuration (`analysis/io_func.py`), analysis helpers (`analysis/an_func.py`), distribution comparison by overlap area and earth mover's distance (`analysis/DistribComp.py`), the yeast thermodynamic model builder (`helper_thermoflux/build_yeast_model.py`) and the 13C problem definitions (`HELPER_FREEFLUX_NF/`).

### `4_results/`

Intermediate and final results written by the notebooks and by the cluster runs. Figures are saved to the plot folders defined in the config modules.

---

## Computational steps run outside the notebooks

The heaviest steps were run on an HPC cluster. The notebooks prepare their inputs and analyze their outputs:

- **Variability analysis**: notebooks `1a`–`1d` write one Gurobi `.mps` problem per reaction; these are solved on the cluster and collected back into the notebooks.
- **Flux sampling (PTA)**: uses the ΔrG° mean and covariance written by `1c`; results are stored in `4_results/b3_flux_sampling/`.
- **13C fitting (FreeFlux)**: each of the 1000 drawn samples is fitted with several labeling strategies (1-13C, U-13C and combined).
- **Classic 13C-MFA confidence intervals** (NLP), used for comparison in `3a`.
- **Metabolite concentration sampling** for the best-fitting flux samples, analyzed in `4h`.

---

## Requirements

Python 3.11, with:

- `thermo_flux` and `equilibrator_api` (thermodynamic model and variability analysis)
- `gurobipy` (requires a Gurobi license)
- `cobra`
- `pta` (probabilistic thermodynamic analysis, for sampling)
- FreeFlux (included in `3_src/freeflux/`)
- `numpy`, `pandas`, `scipy`, `scikit-learn`, `umap-learn`
- `matplotlib`, `seaborn`, `cmap`

---

## Getting started

1. Clone the repository and install the dependencies (also see installation instructions files)
2. Set the project root in the config modules in `3_src/analysis/` (and `BASE_FOLDER` in the stage 1 notebooks) to your local path.
3. Run the notebooks in `2_notebooks/` in numerical order, from the `2_notebooks/` folder. Steps marked as cluster runs need their results in `4_results/` before the following notebooks can run.
