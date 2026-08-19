# Methods — condition-specific metabolic modeling of zebrafish embryos under light regimes

## 1. Data

### 1.1 Transcriptomics
Bulk RNA-seq of whole zebrafish (*Danio rerio*) embryos at 24, 48, 72, 96, and 120 hpf under three light regimes: constant Blue light (BL), constant Dark (D), and 12h/12h Light-Dark cycle (LD). Four biological replicates per condition × timepoint (60 samples total). Raw gene-level counts (32,171 genes × 60 replicate columns) were used as input.

### 1.2 Reference metabolic model
The zebrafish genome-scale metabolic reconstruction `Benji_Wang_Zebrafish-GEM_v2.mat` (a Danio-annotated derivative of the Wang et al. 2021 zebrafish GEM [1], based on the Human-GEM / RAVEN framework) was used unmodified: 12,909 reactions, 8,505 metabolites, 2,714 genes, MARxxxxx identifiers. No cofactor pool modification and no metabolomics-derived exchange bounds were applied at extraction time.

Genes were mapped from ENSDARG identifiers to model gene symbols via Ensembl biomaRt [2] (release matching the zebrafish reference used for Human-GEM annotation), matched case-insensitively against `model.genes`, with paralogs aggregated by MAX (i.e., any paralog with high expression is treated as sufficient evidence for the function). The mapping is stored at `results/qc/gene_map_ensdarg_to_model.csv`.

## 2. Rationale for the thresholding pipeline

### 2.1 Exploratory investigation with TPM (abandoned)

An earlier iteration of this pipeline (retained in `_diagnostics/` for the audit trail; see also `phase3_transcriptomics/`) built gene expression signals in TPM space. That pipeline was:

1. **Gene lengths from Ensembl biomaRt** [2]. For each ENSDARG identifier in the raw count matrix, the sum of `exon_chrom_end − exon_chrom_start + 1` across the canonical transcript's non-overlapping exons was retrieved from the `drerio_gene_ensembl` mart (Ensembl release matched to the annotation used for read alignment). Genes without a length were dropped from TPM computation.
2. **Counts → TPM** using the standard Wagner et al. 2012 formulation [16]:
   ```
   rate_{i,j}      = count_{i,j} / length_i          (reads per base)
   TPM_{i,j}       = 10^6 × rate_{i,j} / sum_i(rate_{i,j})
   ```
   applied per sample j and per gene i.
3. **Median-of-ratios normalization** (DESeq2-style [17]) computed on the TPM matrix directly to remove between-sample composition bias:
   ```
   ref_i          = geometric_mean_j(TPM_{i,j})
   size_factor_j  = median_i( TPM_{i,j} / ref_i )
   norm_TPM_{i,j} = TPM_{i,j} / size_factor_j
   ```
4. **`log2(TPM + 1)`** applied element-wise to the normalized matrix.
5. **Thresholding attempt with GMM** on the ravelled matrix, then falling back to **localT2** [6] (per-gene threshold = the gene's mean across the 12 same-timepoint samples, clipped to the global Q25 / Q75 of nonzero `log2(TPM+1)`).

This exploratory pipeline is reproducible via `phase3_transcriptomics/s01_fetch_gene_lengths.R`, `s02_counts_to_tpm.py`, `s03_normalize.py`, `s04_gmm_fit.py`, and `s05_gene_thresholds.py`; the QC output is `results/qc_report.pdf` (per-stage distributions of raw / TPM / normalized-TPM / log2TPM) and the GMM probe is `results/cpm_gmm_probe.pdf`.

Three empirical observations from this investigation motivated the switch to the CPM-based pipeline described in § 3:

1. **`log2(TPM + 1)` is unimodal in whole-embryo bulk RNA-seq.** TPM's length normalization collapses the count-scale bimodal shape (short vs long genes) into a single expressed hump with a large zero spike. This makes it unsuitable for any threshold-based extractor.
2. **Zero-inflation dominates Gaussian mixture model (GMM) fits.** With ~43% zero-count genes, GMM on `log2(TPM + 1)` fits a delta at zero (σ ≈ 0.02, weight ≈ 30%) as one "component" and the entire real distribution as the other — an artefactual bimodal. Even after `filterByExpr`-filtering the same data, the fit remained a right-skewed unimodal (crossover threshold at the shoulder rather than a valley), confirming that the biological distribution is genuinely unimodal in this system.
3. **Filtering unreliable genes with `edgeR::filterByExpr` before ranking recovers a well-behaved distribution suitable for percentile thresholding.** No natural on/off valley exists in whole-embryo bulk — consistent with the biology of a mixed-tissue sample. Therefore *percentile*-based (rank-based) thresholding is more appropriate than model-fit (bimodality-based) thresholding.

The GMM bimodality probe (`s00b_cpm_gmm_probe.R`) directly compared five normalization variants — `log2(count + 1)`, raw `log2(CPM)`, `log2(CPM + 2)`, TMM-normalized `log2(CPM + 2)`, and our pre-existing `log2(TPM + 1)` — with and without pre-filtering. Every "bimodal" fit on the unfiltered inputs was an artefactual `σ ≈ 0.02` delta at zero + a wide catch-all Gaussian (BIC differences of ~50,000, meaningless because they were driven by the zero spike, not by genuine two-population structure). The filtered `log2(CPM + 2)` gave the only fit with two Gaussians of comparable width (σ ≈ 1.55 vs σ ≈ 2.45), but even that was fitting the right-tail skewness of a fundamentally unimodal distribution, so we adopted a percentile rule rather than a GMM crossover.

### 2.2 Decisions carried forward

The final choice — TMM-normalized `log2(CPM + 2)` with `filterByExpr` pre-filtering and a global 75th-percentile threshold — aligns with the guidance in Robaina Estévez & Nikoloski's 2022 review [4], the Opdam et al. 2017 systematic benchmark of GEM tailoring methods [5], and Richelle et al.'s 2019 assessment of transcriptomic integration decisions [6], all of which support simple global percentile thresholding for GIMME when the underlying distribution is not clearly bimodal, and the 2024 npj Systems Biology benchmark [10] that found between-sample normalization (TMM/RLE/GeTMM) outperforms within-sample methods (TPM/FPKM) for context-specific model extraction.

## 3. Threshold discovery pipeline (`s01_thresholds.R`)

The final normalization + thresholding pipeline is:

1. **DGEList wrapping**: raw counts loaded into an `edgeR::DGEList` [7,8].
2. **TMM normalization**: between-sample scale factors computed with `edgeR::calcNormFactors(method = "TMM")` [9]. TMM was chosen over TPM/FPKM (within-sample) methods based on the 2024 npj Systems Biology benchmark of RNA-seq normalization for context-specific metabolic models [10], which found between-sample methods (TMM, RLE, GeTMM) consistently outperform within-sample methods for downstream extraction accuracy.
3. **Filtering unreliable genes**: `edgeR::filterByExpr(dge, group = cond × hpf)` [11] dropped 18,196 of 32,171 genes (57%) whose counts never rose above sequencing noise across replicates within any of the 15 groups. Filtered genes are *not* discarded from downstream reasoning; they are re-introduced with a low-expression sentinel in step 5 (see § 4).
4. **Log-transform**: `edgeR::cpm(dge, log = TRUE, prior.count = 2)` yielded per-sample `log2(CPM + 2)`. The `+2` pseudo-count follows the edgeR EDA recommendation for smoothing low-count regions [11].
5. **Per-cell per-gene score**: for each of the 15 (condition × hpf) cells, gene scores were computed as the mean of `log2(CPM + 2)` across that cell's 4 replicates. Filtered genes retain `NA` at this stage.
6. **Global threshold** T = 75th percentile of the pooled per-cell means across all kept genes and all cells (`T = 4.882` log2-CPM units on this dataset). The 75th percentile follows Opdam et al. [5] as a widely-used default for GIMME-family extractors, and corresponds to designating the top ~25% of gene × cell measurements as "highly expressed."

Threshold sensitivity was reported at the 60/70/75/80/85/90 percentiles (`results/pipeline/threshold_sensitivity.csv`) to enable downstream sensitivity analyses.

## 4. Two-sentinel gene classification (`s03_gpr_penalty.m`)

Each model gene was assigned a per-cell expression score under a **two-sentinel** scheme that explicitly distinguishes "measured but low" from "not measured":

| gene state in this cell | assigned score | rationale |
|---|---|---|
| measured, `filterByExpr`-kept | actual `log2(CPM + 2)` mean | direct evidence |
| **measured but filtered** | **T − 5** (LOW sentinel) | evidence that the gene is at noise floor → penalize downstream |
| **model gene never measured** | **T + 5** (HIGH sentinel) | no evidence → protect downstream |

This distinction — motivated by mCADRE's three-tier Present/Marginal/Absent framework [12] — avoids the common pipeline bug in which filtered-out genes default to "unknown = protected," which would exempt entire reactions from penalty even when direct measurements indicated their genes were unexpressed.

Per-reaction expression and threshold values were computed by walking each reaction's `model.rules` GPR string with **min-of-AND, max-of-OR** aggregation (E-Flux convention). Paralog aggregation was MAX.

Per-reaction GIMME penalty:

    penalty_r = max(0,  T  −  expression_r)

Reactions with no GPR retained the HIGH sentinel default (no penalty).

## 5. Model extraction (`s04_gimme.m`)

For each of the 15 cells, condition-specific subnetworks were extracted using the COBRA Toolbox implementation of GIMME [3,13]:

    GIMME(model, exprRxns = −penalty, threshold = 0, obj_frac = 0.9)

where COBRA's internal penalty formula `pen_r = max(0, threshold − exprRxns_r)` reduces to `penalty_r` under our `exprRxns = −penalty, threshold = 0` substitution. The `obj_frac = 0.9` requires each extracted subnetwork to preserve at least 90% of the parent model's maximum biomass under its default bounds.

LP solutions were computed with Gurobi 12 via the COBRA Toolbox [13]. Extractions were run serially with per-cell error surfacing.

**Post-extraction validation** per cell:
- Biomass feasibility: `optimizeCbModel(m, 'max')` under a released biomass lower bound (`lb(MAR00021) := 0`); a cell was flagged infeasible if `stat ≠ 1` or `f ≤ 10⁻⁹`.
- Orphan metabolite count: metabolites appearing in exactly one reaction (`sum(S ≠ 0, 2) == 1`), which are a known artifact of GIMME's `|v| < tol` post-LP pruning and cannot be mass-balanced at steady state.

All 15 extracted models were biomass-feasible with post-release biomass max in 8.30–9.04 mmol biomass·gDW⁻¹·hr⁻¹ (parent-scale units), containing 3,394–3,573 reactions each. Orphan-metabolite counts were 842–902 per cell (consistent with GIMME's known pruning behavior [3,5]).

## 6. Downstream analyses

### 6.1 Structural analyses (`s06`, `s07`)
Reaction and gene presence matrices (parent-model rxn/gene × cell) were exported from the 15 extracted models. Subsystem composition was computed from `model.subSystems`. Presence-based Jaccard distance was used for hierarchical clustering (average linkage) and principal coordinates analysis (PCoA) of the 15 cells.

### 6.2 Condition contrast (`s11`, `s12`)
- **PERMANOVA** [14] via `vegan::adonis2` on Jaccard distance with condition + hpf effects (9,999 permutations).
- **Subsystem linear models**: `count ~ condition + hpf` per subsystem; nested F-test for the condition term; Benjamini-Hochberg FDR correction [15] across ~148 subsystems.
- **Planned contrasts**: light (BL+LD) vs dark, Blue-only vs cycled (BL vs LD), LD-specific vs {BL, D}.
- **Hypergeometric enrichment** on condition-unique reactions (present in exactly one condition at a given hpf), per (subsystem × condition), with BH-FDR correction within condition.

### 6.3 Gene and reaction essentiality (`s14`, `s15`, `s16`)
- **Gene essentiality**: for each of the ~840–920 active genes per cell, COBRA's `deleteModelGenes` [13] was used to compute the biomass-feasible flux under KO, respecting boolean GPR structure (AND-linked paralogs not knocked out unless all are removed). A gene was called essential if maximum biomass fell below 1% of the pre-KO baseline.
- **Reaction essentiality**: each reaction in each extracted model was KO'd by setting `[lb, ub] = [0, 0]`; same 1% threshold.
- **Differential essentiality**: genes / reactions essential in exactly one condition at a given hpf. Recurrent hits (essential in the same condition at ≥2 hpf) were highlighted as strong condition-specific dependencies.
- **Subsystem enrichment on differentially essential reactions**: hypergeometric test per (subsystem × condition), universe = reactions essential in any cell, BH-FDR corrected within condition.

### 6.4 Additional exploratory analyses (limitations noted)
Full FVA (`s10`, `fluxVariability(m, 90)` [13]) and per-reaction capacity probes (`s08`) were computed but yielded ranges pinned to the parent model's default ±1000 bounds because no physiological uptake constraints were applied at extraction time. Only structural signals (blocked reactions with FVA range = 0; ~560 per cell) carry biological information at this stage. Physiological constraints from the metabolomics dataset are the natural next step (`s17`).

## 7. Software and reproducibility

- **R 4.4.1** with `edgeR` 4.x, `mclust`, `vegan`, `ggplot2`, `dplyr`, `tidyr`, `pheatmap`, `RColorBrewer`, `jsonlite`, `R.matlab`, `proxy`.
- **MATLAB R2024a** with **COBRA Toolbox 3.x** [13] and **Gurobi 12** as LP solver.
- All scripts numbered s00–s17 in `gimme_transcriptomics_pipeline/`, intended to be executed in numerical order (R and MATLAB steps interleaved as documented in the folder README).
- Raw outputs stored under `results/pipeline/`; human-readable reports under `results/`. Diagnostic scripts from the initial thresholding investigation preserved in `_diagnostics/`.

## References

[1] Wang H, Marcišauskas S, Sánchez BJ, et al. "Genome-scale metabolic network reconstruction of model animals as a platform for translational research." *PNAS* 118 (2021): e2102344118.

[2] Durinck S, Spellman PT, Birney E, Huber W. "Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt." *Nature Protocols* 4 (2009): 1184-1191.

[3] Becker SA, Palsson BØ. "Context-specific metabolic networks are consistent with experiments." *PLoS Computational Biology* 4.5 (2008): e1000082. **[Original GIMME paper]**

[4] Robaina Estévez S, Nikoloski Z. "Guidelines for extracting biologically relevant context-specific metabolic models using gene expression data." *Metabolic Engineering* 74 (2022): 44-58. **[Standard-of-practice review; supports percentile thresholding when data is well-normalized and filter-then-threshold ordering]**

[5] Opdam S, Richelle A, Kellman B, Li S, Zielinski DC, Lewis NE. "A systematic evaluation of methods for tailoring genome-scale metabolic models." *Cell Systems* 4.3 (2017): 318-329. **[Benchmark supporting 75th percentile as a defensible global threshold for GIMME]**

[6] Richelle A, Joshi C, Lewis NE. "Assessing key decisions for transcriptomic data integration in biochemical networks." *PLoS Computational Biology* 15.7 (2019): e1007185. **[Introduces localT2 rule (tried and rejected here) and documents the global-vs-local threshold trade-off]**

[7] Robinson MD, McCarthy DJ, Smyth GK. "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." *Bioinformatics* 26.1 (2010): 139-140.

[8] Chen Y, Chen L, Lun ATL, Baldoni P, Smyth GK. "edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets." *Nucleic Acids Research* 53.2 (2025): gkaf018.

[9] Robinson MD, Oshlack A. "A scaling normalization method for differential expression analysis of RNA-seq data." *Genome Biology* 11 (2010): R25. **[TMM normalization]**

[10] Cornelius SP, Sadigh S, Guo Y, et al. "A benchmark of RNA-seq data normalization methods for transcriptome mapping on human genome-scale metabolic networks." *npj Systems Biology and Applications* 10 (2024): 42. **[Justifies between-sample normalization (TMM/RLE/GeTMM) over within-sample (TPM) for context-specific extraction]**

[11] Chen Y, Lun ATL, Smyth GK. "From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline." *F1000Research* 5 (2016): 1438. **[`filterByExpr` and `prior.count` guidance]**

[12] Wang Y, Eddy JA, Price ND. "Reconstruction of genome-scale metabolic models for 126 human tissues using mCADRE." *BMC Systems Biology* 6 (2012): 153. **[Present/Marginal/Absent evidence framework motivating the two-sentinel design]**

[13] Heirendt L, Arreckx S, Pfau T, et al. "Creation and analysis of biochemical constraint-based models using the COBRA Toolbox v.3.0." *Nature Protocols* 14.3 (2019): 639-702.

[14] Anderson MJ. "A new method for non-parametric multivariate analysis of variance." *Austral Ecology* 26 (2001): 32-46. **[PERMANOVA]**

[15] Benjamini Y, Hochberg Y. "Controlling the false discovery rate: a practical and powerful approach to multiple testing." *Journal of the Royal Statistical Society Series B* 57.1 (1995): 289-300.

[16] Wagner GP, Kin K, Lynch VJ. "Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples." *Theory in Biosciences* 131 (2012): 281-285. **[TPM definition]**

[17] Anders S, Huber W. "Differential expression analysis for sequence count data." *Genome Biology* 11 (2010): R106. **[DESeq median-of-ratios normalization]**
