# 2_notebooks: analysis notebooks and scripts

This folder holds the notebooks and scripts for the thermodynamics-based 13C metabolic flux analysis of *S. cerevisiae*, comparing the wild type (**WT**) with the reduced-glucose-uptake strain (**TM6**). They cover four stages:

1. **Flux space definition** (`1*`): variability analysis (VA) on the thermodynamic–stoichiometric model (TSM) to get flux, ln(concentration) and ΔrG bounds, and comparison with the purely stoichiometric model (SM).
2. **Flux sampling** (`2*`): analysis and quality control of the flux samples drawn from the thermodynamically feasible space (PTA sampling).
3. **13C fitting results** (`3*`): comparison of the sample-based approach with classic 13C-MFA, and analysis of how the 13C fit reshapes the flux distributions.
4. **Biological interpretation and validation** (`4*`): flux splits and turnover, reaction directions compared with Yeast9, prediction of independent isotopomer data (Shen 2024), and metabolite concentrations of the best-fitting samples.

The numeric prefix gives the order in which the files should be run; the letter only orders files within a stage.

---

## Paths and configuration

The notebooks are run from `2_notebooks/` and import the project code from `../3_src/`. Input data are read from `1_data/` and results are written to `4_results/` (see the repository README for the full layout). Several stage 1 notebooks also set the project root explicitly as `BASE_FOLDER`.

Most stage 2, 3 and 4 files get their paths and data loaders from two config modules rather than hard-coding them:

| Module | Used by | Provides |
|---|---|---|
| `analysis.config_pta` (`cfp`) | 2a, 2b, 2c, 3b, 4b, 4h | `RESULTS_FOLDER`, `PLOTS_FOLDER`, `SYNDATA_FOLDER`, `config_prep(cond)` (model + VA flux bounds) |
| `analysis.config_13c` (`config`) | 2d, 3a, 3c–3h, 4a, 4b, 4d–4g | `model`, `net_flux_df(cond)`, `full_flux_df(cond)`, `obj_vals(cond)`, `bootstrap_flux_df(cond)`, folder constants |

Plot styling (palette and shared scatter/diagonal settings) always comes from `io_func.import_plotconfigs()`. WT is drawn with `palette[0]` (variants `palette[1]`, `palette[4]`) and TM6 with `palette[2]` (variants `palette[3]`, `palette[5]`).

---

## Input data

| File | Content |
|---|---|
| `1_data/b_processed_files/yeast_data_xiang_metabolites.xlsx` (sheet `Assembly for Joana`) | Intracellular metabolomics for WT and TM6 (mean, stdev); values below detection are set to 0.01 mM |
| `1_data/b_processed_files/extended_regressiondata_Canelas_Xiang_yeast.xlsx` | Exchange-flux data (`dvmeas`) used as physiology constraints; O2 uptake is ignored |
| `4_results/a_fluxspace/df_bounds_ex.csv` | Extracellular concentration ranges (from Canelas, widened to cover HCO3/CO2) |
| `1_data/b_processed_files/cobra_roundp3_ymodel.xml` | Stoichiometric-only COBRA model (SM) |
| `1_data/b_processed_files/yeast_v3.211_starting_curated_withprotons.xlsx` | Model spreadsheet, used for reaction subsystems |

The thermodynamic model itself is built in code with `helper_thermoflux.build_yeast_model.build_yeast_model()`.

---

## Stage 1: flux space (`1a`–`1f`)

### `1a`–`1d`: variability analysis per strain

Four near-identical notebooks, one per strain × constraint set:

| Notebook | Strain | Constraints | Output prefix | MPS folder |
|---|---|---|---|---|
| `1a_fluxspace_nomet_TM6_hpc2VA` | TM6 | Physiology only, TSM (P.) | `202601_13C_tf01_nomet` | `mpsVA_nomet/TM6/` |
| `1b_fluxspace_nomet_WT_hpc2VA` | WT | Physiology only, TSM (P.) | `202601_13C_tf01_nomet` | `mpsVA_nomet/WT/` |
| `1c_fluxspace_TM6_hpc2VA` | TM6 | Physiology + metabolome, TSM (P. + M.) | `202601_13C_tf01` | `mpsVA_cellconc/TM6/` |
| `1d_fluxspace_WT_hpc2VA` | WT | Physiology + metabolome, TSM (P. + M.) | `202601_13C_tf01` | `mpsVA_cellconc/WT/` |

The "metabolome" variants (`1c`, `1d`) add total cell concentration variables (`thermo_flux.solver.gurobi.total_cell_conc`) bounded by the measured metabolite levels ± CI·sd, on top of the shared steps below.

Each notebook runs the same sequence:

1. **Build the model and set concentration bounds.** Extracellular bounds come from `df_bounds_ex.csv`. Intracellular bounds are widened with the metabolomics data (CI = 5) using compartment volume fractions (WT: c 0.95 / m 0.05; TM6: c 0.90 / m 0.10). NAD(P)H, GTP/GDP and CTP/CDP data are excluded.
2. **Save model indices and ΔrG°.** `1a` and `1c` write `model211_indices.csv`; `1c` also writes `202601_13C_tf01_drG0mean.csv` and `202601_13C_tf01_drG0_cov_sqrt.csv` to `c_synthetic_data/` for the PTA sampling.
3. **Flux VA (run on the HPC).** Exchange fluxes are constrained to mean ± 3·sd, irreversible directions are fixed from the initial ΔrG bounds, then one Gurobi `.mps` problem per reaction is written to the MPS folder (`write_mps=True`). The VA itself runs on the cluster.
4. **Collect HPC results.** Per-reaction result files are read back from `4_results/a_fluxspace/no_metabolome/<cond>/` or `cellconc/<cond>/` and assembled into `<prefix>_flux_bounds_<cond>.csv`. Reactions that timed out are listed; in `1b` they get default bounds of [−500, 500].
5. **ln(concentration) VA** given the flux bounds → `<prefix>_lnconc_bounds_<cond>.csv`.
6. **ΔrG VA** given the flux and concentration bounds → `<prefix>_nocellconc_drG_bounds_<cond>.csv`. In `1b`/`1d`, a few reactions that did not finish are patched by hand with a one-sided bound of ±2000.

### `1e_fluxspace_analysereduction_nometabolome`

Compares flux ranges across the three model types (SM, TSM (P.), TSM (P. + M.)) for both strains.

- Counts forward / backward / reversible reactions from the VA bounds of each model.
- Runs COBRA FVA on the SM with the same exchange-flux constraints.
- Computes relative flux-range reduction (SM → TSM) and box volumes (product of all non-zero ranges).
- Lists the reactions whose ranges differ most between WT and TM6 (shown: PYRDC, FDH, PGI, G3PD2, ATPS3m).
- Writes `4_results/a_fluxspace/df_rxndirections.csv` (reaction directions in both strains).
- Figures: `Fig2a_boxvols_onethird`, `fig2c_densityplots_fluxrange_3models`, `fig2d_fluxbounds_scale_errorbar`, plus reversibility bar charts.

### `1f_stoichiosamples_TM6`

Applies the TM6 TSM (P. + M.) flux bounds to the model and draws 10,000 samples with COBRA's OptGP sampler (stoichiometry only). Output: `TM6_stoichio_fluxsamples_wbounds.csv`, used in `2b` for a supplementary comparison. Note: the code relies on a local edit to `cobra/sampling/hr_sampler.py` (line 198, no model copy).

---

## Stage 2: flux sampling (`2a`–`2d`)

These files read the PTA samples from `4_results/b3_flux_sampling/<cond>/<cond>fluxsamples.csv`. The sampling run itself is not part of these files.

### `2a_globalfluxplot`

- Draws the 1000 random samples used as starting points for the 13C fitting → `c_synthetic_data/2_fluxsamples_for13c/<cond>_drawn_samples.csv` (line commented out).
- Identifies growth-coupled reactions (setting the reaction to zero makes the model infeasible) → `b3_flux_sampling/biomasscoupled_reactions.csv`.
- Scatter of VA range vs. sampled credible-interval range for non-coupled, non-boundary reactions → `fig4d_<cond>_scatter_va_crI`.
- Exports the credible intervals to the manuscript supplement (`Supplementary/samplingCRI<cond>.csv`).

### `2b_fluxdistributions_plots`

- Histograms of sampled flux distributions with the VA bounds overlaid for five example reactions (ACGKm, GLCt+0, GAPD, G3PD2, ATPS3m) → `fig4f_<cond>_kdeplots_sampling`.
- Overlay of thermodynamic (PTA) vs. stoichiometric (`1f`) samples for SUCFUMtm-2-2, ETOHtm+0, ORNtm+1 and MEm → `fig_<cond>_samples_stoichiovthermo_overlay`.

### `2c_traceplots`

Sampler diagnostics on the ΔrG samples (`<cond>_samples.feather`, stored outside the project in `/usr/PhD/2601_13C/` because of their size): trace plots, autocorrelation plots and a normalized step-size ("mixing") metric per reaction → supplementary `202601_13C_<cond>_mixing`.

### `2d_count_revrxn_`

Counts reactions whose direction is undetermined at three levels: VA, the full sample set and the 1000 samples used as 13C input. The percentages in the bar chart are hard-coded from these counts → `fig4e_compare_percent_reversible`.

---

## Stage 3: 13C fitting results (`3a`–`3h`)

Fitting was run on the HPC with FreeFlux (single and combined labeling strategies: 1-13C, U-13C, and U+1-13C). These files load the fitted fluxes and objective values (residuals) through `config_13c`. "Accepted" samples are the best-fitting fraction, set with a percentile cutoff `nperc` (5 or 10 depending on the file); `bootstrap_flux_df` holds the bootstrapped best-fitting set.

### `3a_plotclassic_compare`

Compares flux uncertainty from classic 13C-MFA with the sample-based approach, for WT:

- Hessian-based CIs (`classic_13c/classic13C_fluxranges_hess.json`)
- NLP-based CIs, one min/max optimization per flux on the cluster (`classic_13c/NLPCI/`)
- VA bounds from stage 1
- 95% highest-density credible intervals (HDCI) from the accepted samples

Figures: `fig5e_comparisonclassic_multiplekdeplot_WT` (range distributions for the four methods) and `fig5a_comparisonclassic_3rxn_NLPCIrange_WT` (ICDHym, ALDD2ym, G6PDH2).

### `3b_paramtuning.py` (script)

Tuning of the number of optimizer restarts (4, 8, 16) for TM6 with U+1-13C labeling. Reads `param_tuning/n<N>/` results and plots residuals per sample, mean ± sd vs. number of restarts, and pairwise scatter between restart settings. Outputs go to `SuppFigures/SuppFig_paramtuning_*`. Run with `python 3b_paramtuning.py` (no arguments; condition and paths are set at the top).

### `3c_fitresults_res`

Residuals across labeling strategies (TM6): histogram of 1-13C, U-13C and U+1-13C residuals (`fig5b_residualshist_alllabels_tm6_noface`) and pairwise scatters (`fig5a_3subs_residualsscatters`).

### `3d_plotsmall_kde_new.py` (script)

Command-line tool for the small KDE inset plots of individual reactions:

```bash
python 3d_plotsmall_kde_new.py <condition> <rxn1,rxn2,...> <output_folder> [--plot_all] [--nperc 5 10]
```

`--plot_all` draws the KDE of 200 random samples from the full set; each `--nperc` value adds a curve in the strain's palette. Output: `13C_smallkde_<cond>_<rxn>_plotall<bool>.{png,svg}`.

### `3e_fitresults_dimreduc2cond`

Dimensionality reduction of the sampled flux configurations, colored by fitting residual. PCA is computed but the plotted embedding is PLS regression against the objective values, shown as a hexbin map → `fig5c_<cond>_normPLShex_2D`.

### `3f_fitresults_areanorm_bootstrap`

Quantifies how much the 13C fit changes each flux distribution, using `analysis.DistribComp.DistribComparer` to compare the full sample set with the bootstrapped best-fitting set:

- Overlap area of the aligned distributions (reported as 1 − overlap, "shape variation") and earth mover's distance (EMD).
- Example reaction (MDH) → `fig4_<cond>_overlaparea_example<rxn>_bootstrap`.
- Density of shape variation for both strains → `fig5d_<cond>overlaparea_bootstrap_2conds`.
- Mean shape variation and EMD per subsystem (mapping from `c_synthetic_data/pathways_mapping.csv`) → `fig5e_<cond>_allsubsystems_overlaparea_bootstrap`, `fig5Supp_<cond>_subsystems_EMD_bootstrap`.

### `3g_fitresults_dimreduc_vs_stoiciho`

Samples with and without thermodynamic constraints are compared using dimensionlity reduction (PLS). -> Appendix FIG

### `3h_fitresults_fluxunravelled_bootstrap`

Split violin plots of all samples vs. bootstrapped best-fitting samples for selected reactions (G3PD1, ATPASE, ComplexIII, PFK) → `fig5e_distribchange_violions_KSdiff_bootstrap_<cond>_small`.

### `3i_fitresults_pairwisedist`

Pairwise distances between net flux configurations within each strain, after normalizing to glucose uptake (GLCt+0) and scaling each reaction by its pooled standard deviation. Box plot → `pairwise_distances_boxplot_2strains`.

---

## Stage 4: interpretation and validation (`4a`–`4h`)

All stage 4 files work on the bootstrapped best-fitting flux sets (`config.bootstrap_flux_df(cond)`), usually for both strains.

### `4a_PPPflux`

Oxidative pentose phosphate pathway split, computed as −G6PDH2 / glc-D_EX for every bootstrapped WT sample and reported as a percentage of glucose uptake. No figure; the summary statistics are used in the text. (The markdown cell mentions PGI/(PGI + G6PDH2), but the code normalizes to glucose uptake.)

### `4b_fluxmap`

Median TM6 fluxes of the bootstrapped set for the reactions shown on the flux map: one list for glycolysis, PPP entry and fermentation (HEX1, G6PDH2, TKT1, PGI, GAPD, G3PD1, …, ALCD2m) and one for the TCA cycle and mitochondrial carriers (CS, CSm, ACONT, …, OAAtm-2). The values are read off the notebook output; the map itself is drawn outside Python.

### `4c_comparerxnvsYeast9`

Maps reactions of our model onto Yeast9 (`c_synthetic_data/yeast-GEM9.xml`) to compare predefined directions:

1. Matches by EC number (from the curated model spreadsheet).
2. Adds matches by gene/ORF, using `3_src/analysis/genes_mapping.csv`.
3. Keeps the mapped reactions that are irreversible in Yeast9 and prints both stoichiometries for manual checking of a subset (flagged as wrong direction: ACOAH, ASPTAm, OAAtm-2, PItm-2, SO4t-2).

Outputs (written to the working directory): `y9_mapping_irr_auto.csv` and `full_mapping_rid_rn.json`.

### `4d_fitdir_vs_Y9_bootstrap`

Assigns a direction (forward, backward, reversible) to every reaction in each strain from the sign of its bootstrapped fluxes, then adds the Yeast9 direction from three sources: a hand-curated dictionary in the notebook, `reaction_direction_check.json` and `reaction_direction_check_part2.json` (both built by hand from the `4c` output; 1 = forward, 0 = backward, 2 = reversible).

- Outputs: `rxn_finaldirections_2conds.csv` and `<FREEFLUXRUN_FOLDER>/df_bootstrap_compare_reaction_directions.csv` (both lines commented out).
- Ends with a table of reactions whose fitted direction contradicts an irreversible Yeast9 direction, for WT and TM6.

### `4e_zoom_fitdir_bootstrap`

Reads `df_bootstrap_compare_reaction_directions.csv` from `4d`, lists reactions that are reversible in our fit but constrained in Yeast9, and plots bounded KDEs of WT and TM6 bootstrapped fluxes for PYRtm-1, PItm-2, ASPTAm and SUCFUMtm-2-2 → `SuppFig_kderxn_wrong_predefineddir`.

### `4f_SHENsimulateranges_12`

Validation against independent labeling data: Shen et al. 2024, 50% [1,2-13C]glucose (`a_raw_data/13C_isotopomerdata_Shen2024/data_expmt2_<CENPK|Climited>.xlsx`; WT is compared with CEN.PK, TM6 with C-limited). The `_12` suffix refers to the 1,2-13C tracer.

1. Builds a FreeFlux problem with the 1,2-13C labeling strategy (`c13problem_class`), using the HPC null space (`nullspace_hpc.csv`).
2. For each of the 1000 fitted samples (`<cond>_fit2/optu_<cond>_combined_sample<N>.csv`), sets the total fluxes from the optimized free fluxes and simulates the MDVs of the metabolites measured by Shen → `d_predict_isotopes/<cond>_SHEN_intrametsMDVsall_12_GLC.csv`. Later cells reload the saved file instead of re-simulating.
3. Metabolites present in both cytosol and mitochondria are combined assuming equal concentrations and volume fractions of 0.95/0.05.
4. Compares simulated with measured isotopomer fractions for the bootstrapped samples (`bootstrap_chosensamples<cond>.csv`): regression plot with R² → `fig7a_bootstrap_regplot_simvsmeas_12GLCSHEN`; scatter of all vs. bootstrapped samples → `fig7a_scatter_simvsmeas_12GLCSHEN`; per-metabolite bar plots (glyc3p, g6p, fdp, s7p) → `fig6e_bars_residuals_<met>_12GLC`.

### `4g_met_turnover_bootstrap`

Producing and consuming fluxes of one metabolite (default `atp_c`; `nadph_c` is commented out), from the stoichiometric matrix (`b_processed_files/ymodel_stoichiometric_matrix_df.csv`) and the mean bootstrapped fluxes of each strain. Heatmap of the contributions above 1.5 for WT and TM6 → `SuppFig_<met>_fluxofmet_2conds`. Also prints the share of total ATP consumption taken by PPCK in WT.

### `4h_metconc`

Metabolite concentrations of the best-fitting flux samples, WT vs. TM6:

1. Loads ln(concentration) samples drawn for each bootstrapped flux sample (HDF5 files in `cfp.METSAMPLING_FOLDER/<cond>/`, 10^6 samples in total per strain). TM6 sample 663 is dropped.
2. Draws from the prior Gaussian ln(concentration) distribution (`prep.log_conc_mean`, `prep.log_conc_cov`) to get the initial ranges.
3. Assigns each metabolite to its majority pathway (via `pathways_mapping.csv`).
4. Compares WT and TM6 concentration distributions with `DistribComparer` (complement of overlap area), excluding biomass components and unconstrained metabolites.

Figures: `fig7_histplot_distcomp_wtvtm6` (distribution of the overlap metric), `fig7_distcomp_met<met>` (example: mal-L_m), `fig7_boostrap_lnconc_mean_wtvtm6` (mean ln conc WT vs. TM6, top 10 most different metabolites highlighted).

---

## Data flow summary

```
Input data ──► 1a–1d (VA, HPC) ──► *_flux_bounds_<cond>.csv ──► 1e (SM vs TSM comparison, Fig. 1)
                     │                        │
                     │                        └──► 1f (stoichiometric samples)
                     ▼
         drG0 mean / covariance ──► PTA sampling (external) ──► <cond>fluxsamples.csv
                                                                     │
                                         2a–2d (sampling analysis, Fig. 3) 
                                                                     ▼
                                         2a: 1000 drawn samples ──► FreeFlux fitting (HPC, external)
                                                                     │
                                         3a–3i (fit analysis, Fig. 4 and supplement) 
                                                                     ▼
                                         bootstrapped best-fitting set
                                           ├──► 4a, 4b, 4g (flux splits, flux map, ATP turnover)
                                           ├──► 4c → 4d → 4e (directions vs Yeast9)
                                           ├──► 4f (prediction of Shen 2024 isotopomers, Fig. EV)
                                           └──► ln conc sampling (external) ──► 4h (Fig. 6)
```
