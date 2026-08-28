# Pipeline Log — Zebrafish Condition-Specific GEMs

Chronological narrative of the modeling moves we've made and *why*, so future us (or a collaborator) can pick up the thread without excavating scripts. Complements `METHODS.md` (methods reference) and `README.md` (how to run).

Cells = 15: `{BL, D, LD} × {24, 48, 72, 96, 120 hpf}`
Parent model: `Benji_Wang_Zebrafish-GEM_v2.mat`

---

## CURRENT STATE (as of s51 — canonical method chosen; developmental, condition, and toxicological-stress axes done)

**Pipeline is complete through s51.** The s28–s39 branch (chosen-objective pFBA) has been superseded by objective-free FVA-envelope PCA (s40–s46). s47 uses metabolomics as a probe INTO the canonical model. s48 gives raw descriptive condition-change readouts. s49 (failed forced-flux CYP) → s50 (NADPH-tax CYP stress) → s51 (mechanistic overlay of susceptibility with subsystem envelopes and s44 condition changes).

**Canonical method for developmental biology:** fastCore (unbounded) with free-mode FVA envelope (optPercentage=0). GIMME_bounded free-mode is the orthogonal-extraction confirmation.

**Canonical method for toxicological/nutrient-stress analysis:** the metabolomics-**bounded** variants (GIMME_bounded, fastCore_bounded). Unbounded models cannot be meaningfully stress-tested because they can absorb any perturbation via unrestricted nutrient uptake. This is a separate methodological choice from the developmental-axis analysis — different scientific questions require different bounds.

### What each analysis dimension has established

- **Developmental axis (PC1):** yolk-catabolism → biosynthesis switch, recovered in every method + every objective floor. Cross-modality alignment: fastCore fPC1 ↔ tPC1 r = −0.981, fPC1 ↔ mPC1 r = +0.902. Metabolomics-lens (s47) confirms ~14 individual measured molecules each agree with the model at |r| > 0.9.
- **Condition axis:** no dominant PC captures BL/D/LD, but timepoint-stratified analysis (s44, s48) surfaces specific effects. Biggest measured change: FFA 22:1 at 120 hpf (LD−D = −5.65 log2 AUC). Biggest model subsystem change: peroxisomal β-ox even-chain FA at 120 hpf (LD−D = +1.85). Only unambiguously BL-specific effect: ROS detoxification at 48 hpf (BL−D = +0.60).
- **Measured-vs-model gap:** amino acids and nucleotides show large measured condition deltas at multiple timepoints without corresponding top model subsystem hits (which are dominated by fatty acid pathways). Something the model doesn't resolve yet.
- **Toxicological stress susceptibility (s50):** simulated CYP overexpression via a NADPH+O2 tax reaction. At discriminating stress level ~10,000 mmol/gDW/hr NADPH tax, GIMME_bounded produces a clean binary split: **late-embryo cells (BL_120, D_96, D_120, LD_96, LD_120) fully infeasible; all 24–72 hpf cells fully robust.** fastCore_bounded gives finer gradation with the same age direction. Two independent extraction methods reach the same conclusion. Unbounded fastCore shows no differentiation because unrestricted nutrient uptake absorbs any tax.
- **Mechanistic overlay (s51):** correlates per-cell susceptibility with per-cell subsystem envelope levels to identify which metabolic subsystems' activity predicts robustness vs susceptibility. Cross-references with s44 condition-change subsystems to find pathways that are both condition-plastic AND susceptibility-linked.

### Why the pivot happened

s39 used measured metabolomics exchange rxns as pFBA objectives on the bounded models. Two problems surfaced:
1. **Passthrough loophole.** `addDemandReaction` puts DM_ on the *extracellular* metabolite; the existing EX_ reaction can uptake and DM_ can drain the same molecule. LP shortcut: `EX = -100, DM = +100`, no internal biology used. Every s39 "capacity" number was partly an artifact of this shortcut, most severe on bounded models where EX bounds were symmetric ±100.
2. **Circularity.** Applying metabolomics-derived bounds and then using metabolomics rxns as objectives inflates any correlation between "flux features" and "measured metabolomics."

Stepping back and asking "what's the mathematically defensible object per cell?": the feasible flux polytope $P_i = \{v : S_i v = 0, \ell \le v \le u\}$. The two clean summaries with no arbitrary objective:
- **FVA envelope** = axis-aligned bounding box of $P_i$
- **Uniform sampling** (CHRR) = distribution of interior mass

FVA is cheap and deterministic; CHRR was reserved as a follow-up if FVA proved insufficient. Turns out FVA is enough.

### The new analysis chain (s40 → s46)

- **s40** — Objective-free FVA-envelope PCA. Per-cell feature = per-rxn `log10(1 + fva_range)`. Missing rxns encoded as 0 (absent = zero envelope). PCA per method, cross-modality correlations against transcriptomic PCA (positive control) and measured-metabolomic PCA.
- **s41 / s42** — Filled in missing FVA runs for fastCore (unbounded) and GIMME_bounded so all 4 methods have complete envelopes (parallel across cells, `ensure_solver()` pattern).
- **s43** — Pulled top ± loadings per method × PC (rxn-level + subsystem aggregate). Attached direction: negative loading ⇒ elevated in 24 hpf, positive ⇒ elevated in 48–120 hpf.
- **s44** — Full biological workup: (A) developmental trajectory of top-PC1 subsystems per condition, (B) three-way BL/D/LD condition comparison at each timepoint (per-rxn Δ + subsystem heatmaps).
- **s45** — Sensitivity FVA at two alternative objective floors (`optPercentage=0` = pure feasibility envelope; ATP-hydrolysis as objective at `optPercentage=90`). 120 jobs (4 methods × 15 cells × 2 modes) in a single flat parfor so 14 workers pick up next-available jobs across method/mode boundaries. ~76 min total.
- **s46** — Sensitivity PCA comparing biomass / free / atp modes. Key finding: **`free` and `atp` give identical PC1 (r = 1.000); biomass floor was compressing fastCore signal.** Rebuilt s43 & s44 to use free-mode features for fastCore variants (biomass mode kept for GIMME where the choice makes no difference).

### THE headline result (unbiased, no circularity)

**`fastCore` (unbounded, no biomass floor, no metabolomics bounds):**
- fPC1 ↔ transcriptomic PC1: **r = −0.981**
- fPC1 ↔ measured metabolomics PC1: **r = +0.902**

The developmental axis reconstructs both experimental modalities using nothing but network topology, transcriptomic-derived reaction inclusion, and the raw feasibility envelope. No chosen objective. No circular constraints.

### Biological narrative from s43 + s44

**PC1 = developmental time (yolk → organogenesis switch).** Consistent across all 4 methods:
- **24 hpf (yolk catabolism):** β-oxidation (mitochondrial), carnitine shuttle (ER + peroxisomal), linoleate metabolism, folate metabolism, cholesterol metabolism, BCAA catabolism.
- **48–120 hpf (biosynthetic surge):** cytosolic fatty acid activation, aminoacyl-tRNA biosynthesis, pentose phosphate pathway, pentose/glucuronate interconversions.

**Condition (BL/D/LD) effects are timepoint-specific:**
- **24 hpf:** peroxisomal carnitine shuttle **higher in both light conditions** (BL & LD > D). First light-vs-dark split.
- **48 hpf:** **ROS detoxification BL ≫ D > LD** (mean Δ = 1.7 log-units). Sharpest single condition effect in the dataset. Steady blue light drives oxidative-stress management; cycled light attenuates it.
- **72 hpf:** linoleate metabolism collapses under LD (BL−LD = +0.60).
- **96 hpf:** BL-specific PPP/glucuronate (BL−D = +0.81) and vitamin A metabolism (BL−LD = +0.40).
- **120 hpf:** Inversion — **LD-specific surge** in peroxisomal β-oxidation (mean Δ = 2.2), cholesterol biosynthesis, and linoleate metabolism. Circadian remodeling emerging at larval stages.

### File map for orientation

```
gimme_transcriptomics_pipeline/
  s01–s20g          transcriptomics → GIMME + fastCore → metabolomics bounds
  s21–s25           bounded model FVA/essentiality (superseded by s40+)
  s27               CHRR sampling (paused, laptop-crashing)
  s28–s39           _diagnostic_ — chosen-objective pFBA branch (s39 passthrough issue,
                    biomass-objective degeneracy in s34–s36). Kept for reference.
  s40               objective-free FVA-envelope PCA (canonical analysis)
  s41 / s42         FVA fill-in for fastCore + GIMME_bounded
  s43               PC loadings interpretation (rxn + subsystem)
  s44               developmental trajectory + 3-way condition contrast
  s45               FVA sensitivity: free (optPct=0) + ATP-hydrolysis objective
  s46               sensitivity PCA comparing biomass / free / atp modes
  results/pipeline/*.csv    all intermediates
  results/*.pdf     all reports
```

### Next steps — priority order

1. **Prune s28–s39 into `_diagnostics/`.** They're referenced in the log but no longer active analyses. Keeps the working directory scannable.
2. **Write-up.** Paper figure set from what's on disk:
   - **Fig 1** — s40 PCA panel (4 methods) + transcriptomic + metabolomic positive controls
   - **Fig 2** — s40/s46 correlation heatmaps (flux PCs vs reference PCs, all methods, robustness across objective floors)
   - **Fig 3** — s43 top-loading subsystems for PC1 (yolk-catabolism → biosynthesis axis)
   - **Fig 4** — s44 developmental trajectories + within-timepoint condition heatmaps
   - **Supplementary** — s46 sensitivity table (biomass / free / atp)
3. **Interpret D_96 fastCore outlier** (fastCore_bounded PC1 = +73.7 in s40; also flags in s44). Independent of the D_120 issue from s38.
4. **CHRR (optional).** Only if a reviewer asks for interior-mass sampling on top of envelope. Redesign for per-cell serial (5k samples, thin=100, ~2 hr total) — never re-attempt the batch approach that crashed the laptop in s27.
5. **Deferred / low priority:** pre-extraction vs post-extraction bounds (long-standing question, unlikely to change the story now).

---

## 1 — Transcriptomics normalization

**Move.** Raw counts (Salmon → tximport → edgeR) → filterByExpr (grouped by cell) → TMM normalization → `log2(CPM + prior=2)`.

**Why not TPM.** Initially explored TPM (gene-length divided, then median-of-ratios). Documented in `METHODS.md §2.1`. Abandoned because TMM+CPM is the community standard for cross-sample count normalization when downstream is threshold-based, and TPM's per-sample denominators break cross-sample comparability without extra normalization on top.

**Why filterByExpr on raw, before TMM.** filterByExpr's detection filter (min count in ≥N samples of a group) has to see raw counts to reject noise-floor genes; running it on already-normalized values distorts the group-count logic. TMM library-size scaling comes second, once the noise floor has been trimmed.

**Why log2(CPM + prior=2).** The +2 prior is a stabilizing pseudocount that keeps low-expression genes from producing huge negative values on log scale. Ends up being the input to any expression-driven thresholding.

**GMM alternative considered.** For picking LOW/HIGH sentinels via mixture components on the CPM histogram. Concluded that on filtered+TMM+log2 data, GMM effectively reproduces filterByExpr's bimodal split; two-sentinel percentile cuts are simpler and reproducible.

---

## 2 — GIMME extraction (transcriptomic layer)

**Move.** For each cell independently:
- GPR-to-reaction penalty: gene log2(CPM+2) → aggregate to reaction score using GPR (min for AND, sum for OR)
- Two sentinels: LOW threshold = 25th percentile of scored rxns, HIGH = 75th
- Score below LOW → penalized (proportional to distance below LOW); above HIGH → free; between → linear
- `s04_gimme.m` runs GIMME per cell with biomass as objective, tolerance 90%

**Output:** `extracted_{cell}.mat` (15 files). These are the "transcriptomically bounded" models — bounds inherited from parent, structure trimmed by transcript evidence.

**Why GIMME.** Well-known, protects transcript-supported rxns while allowing low-expression rxns in if biomass needs them. Model discretely differs across cells (~15% Jaccard between conditions).

**Weakness we accepted.** GIMME can prune rxns that leave the network flux-inconsistent if we later constrain harder. That's what motivated adding fastCore.

---

## 3 — fastCore extraction (structural integrity layer)

**Move.** `s18_fastcore.m` → `s18b_fastcore_strict.m` (current). fastcc first to get the flux-consistent subnetwork of the parent, then fastCore with a **core set** = utilities ∪ metabolomics exchanges ∪ transcript-high (score ≥ HIGH threshold) rxns. No-GPR reactions are included by default (evidence-neutral, protected).

**Output:** `extracted_fastcore_{cell}.mat` (15 files). Guarantees every rxn carries flux under some feasible LP — needed for CHRR sampling, FVA, and any downstream LP-based analysis.

**Why fastCore alongside GIMME.**
- GIMME captures transcript-driven biology (~15% Jaccard differentiation between conditions, condition-informative)
- fastCore gives structural integrity (~3–5% Jaccard between conditions, less differentiated but every rxn is functional)
- Two lenses on the same data; use each where it's stronger

**Metabolomics/utility coverage in fastCore models.**
- 101/105 metabolomics exchanges present in every cell (4 blocked in the parent's flux-consistent subnetwork)
- 18/24 utilities present (6 blocked E3 medium components have no internal consumers in parent)

**The "strict" iteration.** First strict cut excluded no-GPR rxns → biomass crushed to 0 in most cells. Fix: no-GPR rxns are protected (evidence-neutral) rather than excluded. Model sizes settled around ~6000 rxns.

---

## 4 — Metabolomics: geomean-relative bounds (foundational discretization step)

Metabolomics has been mostly used as a *bound-setting* mechanism, not yet as an analysis lens. This is important foundation for discrete modeling.

**Move.** `s20_apply_relative_metab_bounds.m` (fastCore) and `s20g_apply_relative_metab_bounds_gimme.m` (GIMME). "Option C, symmetric":

```
geomean_flux(rxn)     = geomean(|AUC-derived flux(rxn, cell)| across 15 cells)
fold(rxn, cell)       = |flux(rxn, cell)| / geomean_flux(rxn)
magnitude(rxn, cell)  = BASE_FLUX × fold(rxn, cell)
lb = -magnitude,  ub = +magnitude
```

**Why geomean, not one cell's baseline.** Any single-cell reference biases the whole panel toward that cell. Geomean over cells is symmetric, cell-democratic, and rewards the relative ordering that AUC measurements actually deliver.

**Why symmetric.** Metabolomics measures abundance, not transport direction. We do not know whether a metabolite is being taken up or secreted at that timepoint. Bounding as a "carrying capacity" in either direction lets the LP pick sign.

**Why absolute magnitude doesn't matter.** User's own words: "the biomass is a pseudo-reaction, as long as the concentrations are relative to one another, all is good." BASE_FLUX=10 is an arbitrary scale; only the fold differences carry information across cells. Original parent bounds were ±1000 — the geomean-relative bounds are much tighter but the *relative* ordering across cells is what's biologically meaningful.

**Output:** `extracted_fastcore_bounded_{cell}.mat`, `extracted_bounded_{cell}.mat`. Bounds review workbook: `results/s20_bounds_applied.xlsx` (LONG, MAGNITUDE, LB, UB, FOLD sheets, 101 rxns × 15 cells wide).

**Open question:** BASE_FLUX sensitivity. Parent was ±1000; we're at 10. Sweep {1, 10, 100, 1000} not yet done. User said "no don't patch" for now.

---

## 5 — Downstream on bounded models (already run)

Scripts `s21`–`s25`: capacity probes, FVA, gene/rxn essentiality, all on `extracted_fastcore_bounded_{cell}.mat`.
Report `s26_bounded_report.R` produced a PDF that felt unsatisfying — didn't surface biology, mostly recapitulated bound differences. Left on the shelf; may rework.

## 6 — CHRR sampling attempt (paused)

Scripts `s27_run_chrr.m`, `s27b_differential_flux.m`, `s27_chrr_report.R` built. First run crashed the laptop (30 parallel model copies + CHRR memory footprint). Back-burnered until we scale down (fewer chains, fewer workers, or serial).

---

## 7 — CURRENT MOVE: metabolomics as a **lens**, not a bound

New tack: use the 105 measured metabolomics exchanges to interrogate the *transcriptomically-bounded* fastCore models (before any metabolomics bounds are applied). Two scripts:

### `s28_metab_shadow_prices.m`
Max biomass, read metabolite shadow prices (LP duals) and reaction reduced costs, filter to the 105 metabolomics rxns.
- Non-zero |shadow price| ⇒ the model already treats this metabolite as biomass-limiting under transcriptomic constraints alone
- Correlate |shadow price| with measured AUC across cells ⇒ does the model "care" about mets we actually see?

**Output:** `s28_shadow_prices_long.csv` (cell × 105 rxns, dual values), `s28_biomass.csv`.

**Result (2026-08-21):** Every one of the 15 cells returned |shadow price| = 0 on all 101 in-model metabolomics exchanges. Interpretation:

- Under unbounded fastCore (parent ±1000, everything wide open), none of the 105 measured metabolites are currently biomass-limiting.
- Biomass in the 6–8 range is capped by non-metabolomics exchanges or internal stoichiometry, not by anything we measured.
- **The metabolomics data is off the LP's critical path unless we impose bounds.** This is the empirical argument *for* the s20 geomean-relative bounds — without them the measurements are invisible to the model.
- Biomass values: BL_24, D_24, LD_24 = 6.146; most 48–120 = 7.078; 72 hpf (BL/D) = 8.009. Three discrete plateaus suggest shared internal bottlenecks in biomass composition, not exchange limitations.

### `s29_single_source_screen.m`
For each measured met M, in each cell:
- **Mode A (permissive):** close only *other* metabolomics exchanges, allow M ± SOURCE_CAP, max biomass → M's marginal contribution above base medium
- **Mode B (harsh):** close all exchanges except M, max biomass → can M alone sustain the model?

Rank mets per cell by productivity. Correlate ranking to measured AUC.

**Output:** `s29_single_source_permissive_long.csv`, `s29_single_source_harsh_long.csv`, `s29_baseline_bio.csv`.

**Result (2026-08-21):**
- **Harsh mode:** every cell returned 0/101 nonzero. No measured metabolite can singly sustain biomass — consistent with the fact that the metabolomics panel is nucleosides / small biochemicals, not central C/N sources.
- **Permissive mode:** 0 or 1 hits per cell out of 101. Closing 100 of the measured mets barely moves biomass.
- Reinforces s28: unbounded fastCore models are effectively indifferent to the metabolomics measurements. To use these metabolites as a lens on biology, we have to **impose** them as bounds (as s20 does), not query them structurally.

### `s28g` and `s29g` — GIMME mirrors
Same probes run on `extracted_{cell}.mat` (GIMME extractions) instead of fastCore. Filenames prefixed `s28g_` / `s29g_`. Purpose: confirm whether the shadow-price = 0 pattern is a fastCore artifact or holds for the tighter transcript-driven extraction as well.

**Result (2026-08-21):**
- GIMME model sizes ~3500 rxns (half of fastCore's ~8500); only 56–73 of 105 measured mets survive extraction (fastCore kept 101/105)
- Baseline biomass higher and more continuous: 8.3–9.0 vs fastCore's 6.1/7.1/8.0 plateaus
- Single-source screens still 0/0 across all cells — no measured metabolite can singly sustain biomass in GIMME either
- **Shadow prices: exactly one hit.** `MAR09292` = **carnitine** (MAM02348e, HMDB0000062) has nonzero shadow price in 7/15 cells:
  - 24 hpf hits in all three conditions (BL, D, LD)
  - 120 hpf: no cells
  - Middle timepoints scattered: BL_48, D_72, D_96, LD_96
  - Shadow price ~-7.1e-4 in every hit — small but consistent negative sign; model is at a slight equilibrium tension around carnitine availability
- **Biological reading:** carnitine shuttles long-chain fatty acids into mitochondria for β-oxidation. The fact that this is the ONE measured metabolite the GIMME models feel tension on — and only at early/mid timepoints — is consistent with early-embryonic energy demand modulation. This is a genuine transcript-driven structural finding not created by any imposed metabolomics bound.
- **Method takeaway:** fastCore's structural bloat masks LP tensions that GIMME's transcript-tight extraction reveals. When the goal is "which measured mets is the model actually sensitive to," GIMME is the better lens; fastCore is better when the goal is "make everything carry flux for sampling."

## 8 — Per-pool biomass decomposition (s30 / s30g)

**Motivation.** Whole-biomass shadow prices lump ~32 precursor demands together and hide pool-specific tensions. Solution: probe each biomass precursor independently by installing a temporary demand rxn `DM_<met>: met ->` and maximizing it. Read max capacity + shadow prices on the 105 measured mets per (cell × precursor).

**Scripts.** `s30_biomass_pool_capacity.m` (fastCore), `s30g_biomass_pool_capacity_gimme.m` (GIMME). Both are parfor-across-cells with `ensure_solver()`. Each cell contributes 32 precursor probes × 105 measured-met shadow prices.

**Result (2026-08-21):**
- Same 32 biomass precursors in both extractions (identical MAR00021 structure).
- **fastCore: 0 hits across all 480 (cell × precursor) probes.** Every pool has slack — the bloated network absorbs any local demand without binding measured mets.
- **GIMME: 10 hits, all on carnitine (MAM02348e / MAR09292), tied to two specific pools:**
  - **MAM02392c = lipid droplet** — 4 hits (BL_24, BL_48, LD_24, LD_72), shadow price ~−0.038 (**50× stronger than the whole-biomass probe in s28g**)
  - **MAM01602c = cofactors and vitamins** — 6 hits (BL_24, BL_48, D_24, D_96, LD_24, LD_96), shadow price ~−0.0006

**Biological reading.**
- Lipid droplet biosynthesis requires LCFAs; carnitine shuttles LCFAs across membranes → the coupling is mechanistically expected and here it shows up as a hard LP constraint under GIMME's tight extraction.
- 24 hpf hits both pools in all three conditions → early-embryonic lipid mobilization is transcript-programmed to depend on the carnitine shuttle regardless of light.
- Middle-timepoint lipid-droplet hits are BL/LD only (not constant dark) → light conditions extend the lipid-droplet-carnitine coupling into mid-embryogenesis.
- 120 hpf: no hits in either pool → dependency dissolves as embryo matures.

**Method takeaway.** Per-pool probing amplifies the signal by 50× (from −7e-4 whole-biomass to −0.038 lipid-droplet). When looking for metabolite-specific biology, always decompose lumped objectives.

**Outputs.**
- `s30_pool_capacity.csv` / `s30g_pool_capacity.csv` — long: cell, precursor_met, biomass_stoich, max_capacity
- `s30g_pool_shadow_prices.csv` (fastCore had no hits, so no file)
- `results/s30_biomass_pool_report.pdf` — via `s30_biomass_pool_report.R`

## 9 — Biomass tree cascade (s31, s32, s33)

**Discovered the biomass bottleneck** by walking the entire precursor tree of MAR00021.

**Two parallel biomass systems in the model:**
- MAR00021 (flat, 32 direct precursors) — WORKS, this is the one we use
- MAR13082 + MAR10062–65 (pool-based, sub-pool assemblies) — DEAD (reactions exist but blocked in parent flux-consistent subnetwork; probably missing/broken tRNA-charging or phospholipid assembly steps)

**The 32 precursors of MAR00021 split into two classes:**
- **Class A (heavily cycled hubs):** all 20 amino acids + cholesterol. Turnover 4,000–74,000 vs biomass demand of 9.04. Peptide formation/degradation cycles (`MAR11xxx` dipeptides at ±1000) run in loops. Structurally cannot register LP tension.
- **Class B (narrow-neck):** 8 precursors — lipid droplet, cofactors and vitamins, phosphatidate-LD-TAG, CL pool, RNA, DNA, DNA-5-mC, glycogen. Turnover exactly 2× biomass demand. All the LP signal lives here.

**THE bottleneck (from s32 on parent):** `cofactors and vitamins (MAM01602c)` via MAR00022. Cap = 9.04 = biomass max. Under it, cytochrome-C (MAM01631m) sub-caps at 9.044.

**Per-cell cascade (s33):** all 30 extraction models cap on the same axis:
- fastCore: bottleneck at Level 1 (cofactor pool via MAR00022). Biomass values collapse to three plateaus (6.15 for all 24 hpf cells, 7.08 for most, 8.01 for BL_72 + D_72 only).
- GIMME: bottleneck at Level 2 (cytochrome-C via MAR04762). Biomass splits into two clusters — 9.04 for mid-embryogenesis (48–72 hpf) and 8.30 for edges (24 hpf + 96–120 hpf).

**Biological signal:**
- Bathtub curve in GIMME (24 and 96–120 hpf constrained, 48–72 hpf peak) → cytochrome-C sub-pathway is transcript-supported in mid-embryo, pinched at edges.
- LD_72 fastCore drops off the 8.01 plateau (stays at 7.08) while BL_72 + D_72 reach 8.01 → light-dark cycling at 72 hpf specifically constrains the cofactor pathway.

**Cofactor artifacts to filter:** naive lowest-cap ranking picks FAD, GTP, sedoheptulose-BP, tagatose-BP with cap = 0. These are cofactors that CYCLE (no net production) and cannot be drained by demand probes. Correct ranking = node with cap closest to biomass_max (not cap = 0). Filter added post-hoc in interpretation.

**Outputs.**
- `s32_biomass_tree_long.csv` + `s32_biomass_tree_dump.txt` — parent tree characterization
- `s33_cell_tree_long.csv` — per-cell trees (30 models × ~130 nodes = 3769 rows)
- `s33_biomass_per_cell.csv` — biomass max per (extraction, cell)
- `s33_bottleneck_ranking.csv` — top-10 lowest-cap narrow-neck nodes per cell (unfiltered; needs post-hoc filter of cap>0 and cap<999)

## 10 — Bounded biomass tree cascade with parent controls (s34)

**Why.** s33 characterized per-cell biomass bottlenecks on *unbounded* extractions. We also need to know whether the metabolomics bounds we imposed in s20/s20g change where the bottleneck lives. And we need a bounded-comparison-appropriate baseline — hence a parent model with all exchange bounds clamped to ±100 (matching the magnitude regime of the metabolomics bounds), so any drop in biomass on the bounded extractions can be attributed to *extraction pruning* vs *general exchange tightness*.

**Scripts.** `s34_biomass_tree_bounded.m` (MATLAB cascade on parent + parent_scaled100 + 15 fastCore_bounded + 15 GIMME_bounded); `s34_biomass_bottleneck_report.R` (two-axis PDF: axis 1 = timecourse 24–120 hpf, axis 2 = within-timepoint condition comparison).

## 11 — Tree feature-space PCA (s35)

**Why.** Collapsing each cell to biomass + one bottleneck was proof-of-concept; the actual data richness in a cascade is the full ~130-node feature vector per cell. Rather than reducing to one number, treat the tree as an input feature matrix and do PCA. If cells cluster by hpf or condition, the feature space is picking up biology.

**Scripts.** `s35_tree_feature_pca.R` — builds three feature matrices per cell (max_capacity_alone, flux_demanded, utilization = flux/cap), one PCA per method (own scale), plus loadings and a Kruskal-Wallis differential-node analysis within timepoint.

**Method note (post-hoc).** PCA on this feature space produced weak clustering. Reason: features are all downstream of the same biomass LP and dominated by the cofactor/cyt-C bottleneck axis, so effective dimensionality is low; LP degeneracy adds solver-tiebreaking noise on top. Not the right feature space for developmental biology. **This finding motivated s36.**

## 12 — Component-forward pFBA feature space (s36)

**Why.** Biomass optimization has three structural problems for feature extraction: (1) it always solves via the cofactor/cyt-C bottleneck, hiding perturbations elsewhere; (2) it's LP-degenerate with many equivalent flux vectors, so "flux at biomass optimum" is arbitrary tie-breaking; (3) 20 of 32 precursors are amino acids running futile cycles at ±1000, contributing zero LP signal. The fix: don't optimize biomass. Treat each of the 32 biomass precursors as its own objective. For each cell × precursor, install a temporary demand rxn, max it, then pFBA to get a *unique* minimum-flux solution. Reduce each pFBA flux vector to subsystem-level activity (sum |flux| per subsystem) to get an interpretable feature vector.

Feature space grows from ~100 dims (single-biomass tree) to ~1000 dims (32 precursors × ~100 subsystems), and every dimension is a distinct pFBA solve with no biomass bottleneck sitting in the middle.

**Scripts.** `s36_component_pfba_panel.m` (parallel pFBA panel on 62 models × 32 precursors = 1984 solves); `s36_component_pfba_report.R` (per-method PCA on the resulting feature matrix, labeled loadings).

## 13 — Transcriptomic vs flux PC correlation (s37)

**Why.** If the flux-based feature space is discovering biology, its top PCs should correlate with the transcriptomic PCs (which we know cluster developmentally and by condition). This is a **cross-modality validation** — do model-derived flux features align with the same axes of variation the transcriptomics identifies? If yes, the metabolic model is faithfully translating gene expression into flux behavior. If no, the model is decoupled from the transcript layer and something's off.

**Scripts.** `s37_transcriptomic_flux_correlation.R` — loads the per-cell log2(CPM+2) matrix, PCAs it, then correlates the top transcriptomic PCs against each method's top flux PCs (Pearson + Spearman). Reports which flux PC best aligns with which transcriptomic PC.

## 14 — Outlier drill-down: D_120 (s38)

**Why.** In s36, fastCore_bounded PC3 flagged D_120 as an extreme outlier (score −19.4 while all other cells were within ±13). Loadings pointed at cofactor × Folate/Ubiquinone/Prostaglandin + cholesterol × Heme degradation. Before writing up any biology from the PCA, we should check: is D_120 real biology (constant-dark constant-longest embryonic timepoint has a specific metabolic state), or is it a data/model artifact (extraction anomaly, bound scaling issue, one-off numerical instability)? Either way it needs to be characterized so we know whether to include or flag it.

**Scripts.** `s38_D_120_outlier.R` — pulls all D_120 features from s36, ranks by deviation from cross-cell mean, identifies distinct subsystems, and cross-references against BL_120 + LD_120 (same hpf) and D_24/48/72/96 (same condition) to isolate whether the signal is condition-specific, timepoint-specific, or truly a D_120-only pattern.

## 15 — Metabolomics-anchored pFBA panel (s39)

**Why.** s36 used biomass precursors as objectives; s39 uses **measured metabolomics exchange reactions** as objectives instead. Different question:
- **s36 asks:** what does the model do to make its own biomass?  (Self-referential to the model's objective.)
- **s39 asks:** what pathways does the model activate to handle each of the 105 metabolites we actually measured?  (Data-anchored to the experimental measurements.)

The output is a third feature space distinct from A (biomass precursors) and complementary to it. Where A's features are model-internal, s39's features are directly tied to the metabolites we can validate against.

**Cross-validation.** s37 correlated s36 flux PCs against transcriptomic PCs — that established the transcriptomic-flux modality bridge. s39 extends this: correlate the metabolomics-anchored flux PCs against BOTH (a) the transcriptomic PCs (same as s37) and (b) a fresh PCA on the **measured AUC-derived fluxes themselves**. If the model's flux features PCA correlates with the *actual* metabolomics data PCA, that's cross-validation from a totally independent source — model vs data, no transcript intermediate.

**Scripts.** `s39_metab_pfba_panel.m` (MATLAB, 62 models × 105 mets = 6510 pFBA solves in parallel); `s39_metab_pfba_report.R` (PCA per method + transcriptomic correlation + metabolomics correlation).

**Why this matters.** Answers a question we should have asked before bounding: do the metabolites we measured actually align with what the transcriptomic models consider metabolically important? If yes, the geomean bounds are well-motivated. If not, the two data layers are decoupled — which is *also* a finding worth having on record.

---

## Deferred / open

- **Pre-extraction vs post-extraction bounds:** apply metabolomics bounds to parent before extraction, then extract. Compare model sizes / biomass / subsystem coverage.
- **BASE_FLUX sensitivity sweep** ({1, 10, 100, 1000}).
- **CHRR redo** at smaller scale after we understand what memory was consumed.
- **s26 report rework** into a tighter biology-forward PDF.
- **s28/s29 report** — once s28+s29 finish, build an R script to plot cross-cell heatmaps, AUC correlations, and metabolite rankings.

---

## 16 — s40: FVA-envelope PCA (objective-free per-cell fingerprint)

**Why the pivot.** Everything from s28 onward chose an objective (biomass, biomass precursors, measured metabolites) and reported flux distributions conditional on that choice. Every choice is a modelling assumption. s39 additionally had a **passthrough loophole**: `DM_<met>` was added on the extracellular metabolite while `EX_<met>` was left reversible, so the LP could satisfy demand by uptake→drain without using any internal biology. Every s39 "capacity" was partly a shortcut, worst on bounded models where EX bounds were symmetric ±100.

Stepping back mathematically: each cell has a feasible flux polytope $P_i$. Two clean summaries with no arbitrary objective: **FVA envelope** (axis-aligned bounding box) and **uniform sampling** (interior distribution). FVA is cheap and deterministic. Try it first.

**Script.** `s40_fva_pca.R`: per-cell feature vector = per-rxn `log10(1 + fva_range)`, missing rxns encoded as 0 (topology + capacity combined). PCA per method, cross-modality correlation with transcriptomic PCs and measured-metabolomic PCs.

**Result.** Every method's PC1 aligns with both transcriptomics and metabolomics — objective-free.

---

## 17 — s41 / s42: fill-in FVA runs

**Why.** s10 covered GIMME unbounded; s23 covered fastCore_bounded. To have a complete 4-method panel for s40, ran FVA on fastCore (s41) and GIMME_bounded (s42). Same OPT_PCT=90, biomass floor released, parfor across cells, `ensure_solver()` on each worker.

---

## 18 — s43: PC-loading interpretation

**Why.** s40 established that PC1 is meaningful; s43 answers "meaningful *how*." Pulls top ± rxn loadings per method × PC (1..3), aggregated to subsystem, with cell-score bar so we know which sign of the loading corresponds to which cells (24 hpf are always negative on PC1 → negative-loaded pathways are 24 hpf-elevated).

**Script.** `s43_fva_pc_loadings.R`. After s46 (see below) the fastCore variants were switched to load from `fva_range_free_{method}.csv` because the biomass floor was compressing their signal.

**Finding.** PC1 = yolk-catabolism → biosynthesis switch. Reproduces textbook zebrafish embryo metabolism from FVA envelopes alone.

---

## 19 — s44: developmental trajectory + 3-way condition comparison

**Why.** s43 showed PC1 = developmental time, invariant across BL/D/LD. But PC2/PC3 didn't cleanly resolve condition. Question: is there no condition signal, or is it just hidden by the developmental gradient? s44 fixes each timepoint and asks "at hpf=X, which reactions/subsystems differ most between BL, D, LD?" — plus tracks the trajectory of top-PC1 subsystems over hpf, per condition.

**Script.** `s44_trajectory_and_condition.R`.

**Finding.** Condition signal is timepoint-specific:
- **24 hpf:** peroxisomal carnitine BL & LD > D
- **48 hpf:** ROS detox BL ≫ D > LD (sharpest single effect, mean Δ = 1.7)
- **72 hpf:** linoleate collapse under LD
- **96 hpf:** BL-specific PPP + vitamin A
- **120 hpf:** LD-specific peroxisomal β-ox + cholesterol biosynthesis surge (mean Δ = 2.2)

---

## 20 — s45: FVA sensitivity (biomass-free + ATP-anchored)

**Why.** s10/s23/s41/s42 all used `optPercentage=90` with biomass as objective — every envelope constrained to keep biomass ≥ 0.9 × max. Reviewers will (correctly) ask whether that choice drove PC1. Cheap sensitivity test: rerun FVA at `optPercentage=0` (no objective floor, pure feasibility envelope) and with ATP hydrolysis as the objective at `optPercentage=90` (energy-anchored, biologically defensible alternative).

**Script.** `s45_fva_sensitivity.m`. 4 methods × 15 cells × 2 modes = **120 jobs in one flat parfor** so 14 workers pull the next available job as soon as they finish (no barrier between methods or modes — critical for wall-time efficiency on a mixed job size). `MAR03964` for ATP in fastCore, `MAR06916` in GIMME (falls back to `DM_MAM01371c` demand rxn if neither exists). ~76 min total.

---

## 21 — s46: sensitivity PCA + free-mode adoption for fastCore

**Why.** Compare s40 (biomass) vs s45 (free) vs s45 (atp) at the PCA level. If PC1 loadings and cross-modality correlations survive across all three floors, biomass was innocent.

**Script.** `s46_sensitivity_pca.R`.

**Findings:**
- `free` and `atp` give **identical PC1** in all 4 methods (r = 1.000). ATP constraint is non-binding on top of pure feasibility.
- **GIMME variants perfectly robust** to objective choice (biomass↔free PC1 alignment r > 0.999). Biomass floor was a no-op.
- **fastCore variants: biomass floor was hurting the signal.** Releasing it pushes fastCore fPC1↔transcriptomic from −0.92 to **−0.98** and fPC1↔metabolomic from +0.76 to **+0.90**. fastCore_bounded metabolomic correlation goes +0.80 → **+0.95**.
- Consequence: s43 and s44 rebuilt to load `fva_range_free_{method}.csv` for fastCore variants (GIMME unchanged).

**Headline result to report (clean, no circularity):** fastCore unbounded, free-mode envelope → fPC1 vs transcriptomic PC1 r = **−0.981**, fPC1 vs measured metabolomics PC1 r = **+0.902**.

---

## 22 — s47: metabolomics as lens INTO the canonical model

**Why.** s46 settled the canonical method (fastCore unbounded, free-mode envelope). The bulk model↔metabolomics correlation is r = +0.90. That number is macroscopic — it doesn't tell you which measured molecules the model "sees" well vs poorly, nor does it tell you which model subsystems each measurement is coupled to. s47 flips the direction: uses the measured metabolomics as a **probe** into the model, rather than as a validation reference for the model's output.

**Three analyses in one script (`s47_metabolomics_lens.R`):**

- **Part A** — Metabolite × subsystem coupling map. For each of 176 measured metabolites (with model ex_rxn mapping) × 133 subsystems (n_rxn ≥ 3), Pearson correlation of concentration vs mean-subsystem-envelope across 15 cells. 20,300 pairs written to CSV.
- **Part B** — Metabolic PC anchor drill-down. Top-20 metabolites loading on measured mPC1 and mPC2 (deduped by molecule to remove LC-column variants). Envelope agreement per molecule: uses the ex_rxn if it exists in the model, falls back to the metabolite's top-correlated subsystem envelope from Part A when the metabolite has no exchange reaction.
- **Part C** — Metabolomically-filtered PCA rigor check. Restrict flux features to the 20 subsystems most linked to the top-CV measured metabolites (~416 of 12909 rxns). Compare filtered PC1 vs full PC1 and vs measured mPC1.

**Findings:**

- **Part A**: strongest single link is `glucose ↔ Linoleate metabolism` at r = −0.994. Multiple amino acids (Ornithine, Arginine, Alanine, Glycine, Threonine) couple to fatty acid subsystems (Linoleate, Arachidonic acid, FA desaturation) at |r| > 0.95. Consistent with amino acid ↔ lipid remodeling coupling during yolk→organogenesis transition.
- **Part B**: after dedup + subsystem fallback, ~14 of 20 top mPC1 loaders agree with the model at |r| > 0.9. The r = +0.90 macro correlation is carried by many molecules agreeing simultaneously — not a handful of well-fit outliers.
- **Part C**: filtered PC1 vs full PC1 r = 0.987. The developmental axis is carried by exactly the ~3% of reactions in metabolism visible to the assay — not hiding in unfalsifiable model territory.

**Implication for reporting.** Part A gives a defensible per-molecule and per-subsystem interpretation table. Part B upgrades the single macro correlation into ~15 individually-verified per-molecule correlations. Part C rebuts the "your correlation is coming from stuff you can't measure" objection.

---

## 23 — s48: descriptive condition-change readout (measured vs model, side-by-side)

**Why.** s47 finished the metabolomics-lens work on the developmental axis. Before layering biological interpretation on condition (BL/D/LD) differences, we wanted a purely descriptive readout: **at each hpf, what changed most in the measured data, and what changed most in the model subsystem envelope, side-by-side.** No interpretation, just the top hits and their signed contrasts (BL−D, LD−D, BL−LD). Grounds any condition biology claims in the raw data.

**Script.** `s48_condition_data_readout.R`. Reads `metabolomics_group_means.csv` (measured log2 AUC per cell) and `s44_subsystem_condition_deltas.csv` filtered to canonical method (fastCore free-mode). For each timepoint × top 12 measured metabolites × top 12 model subsystems, computes the three pairwise condition contrasts, outputs a grid of grouped bar plots + a per-hpf side-by-side table.

**Descriptive findings (no biology yet):**
- **96 hpf and 120 hpf** are the two timepoints with the largest condition deltas in the dataset, on both measured and model sides.
- **Largest single measured molecule change anywhere in the dataset:** FFA 22:1 at 120 hpf, BL−D=−3.91, LD−D=−5.65. FFA 22:2 shows the same-direction pattern.
- **Largest single measured non-lipid change:** picolinic acid at 96 hpf, BL−D=−2.54, LD−D=−3.58.
- **Largest single model subsystem change:** peroxisomal β-oxidation of even-chain FA at 120 hpf, LD−D=+1.85. Peroxisomal β-ox of unsaturated n-9 (+1.76) and cholesterol biosynthesis 3 (Kandustch-Russell) (+1.53) top the same timepoint.
- **Only clearly BL-specific model subsystem effect (not "light in general"):** ROS detoxification at 48 hpf, BL−D=+0.60, LD−D=0.00.
- **Coordinated nucleotide movements (AMP, UDP, UMP, GMP, adenine)** appear as measured top hits at 96 and 120 hpf, but do NOT surface as top model subsystems at those timepoints — a measured-vs-model gap worth noting.
- **Amino acid movements** (glutamic acid, N-acetylaspartate, 2-oxoglutarate) are prominent on the measured side across multiple timepoints but rarely top the model subsystem list, which is dominated by fatty acid subsystems.

**Outputs:**
- `results/s48_condition_data_readout.pdf` — per-timepoint side-by-side bar plots (measured left, model right, 3 contrasts per feature) + summary bars of peak and median-top-10 delta per timepoint + per-hpf tables.
- `results/pipeline/s48_measured_condition_deltas.csv` — full measured-side long table.
- `results/pipeline/s48_matched_side_by_side.csv` — top 12 measured + top 12 model per hpf in one flat file.

---

## 24 — s49: CYP overexpression susceptibility (failed formulation, kept for reference)

**Why we tried it.** Toxicological stressors upregulate CYP (cytochrome P450) family enzymes as a hallmark detox response. Question: are any of our 60 cell-specific models more/less able to sustain forced CYP flux — i.e., which developmental stages or light conditions produce metabolic states that would be more vulnerable to a toxicological insult?

**What we did.** Defined a "CYP reaction set" from the model's subsystems:
- **Narrow set**: Xenobiotics metabolism + Drug metabolism (~130–620 rxns per model, method-dependent)
- **Broad set**: adds Retinol, Steroid, Estrogen, Androgen, Glucocorticoid, Vitamin D, Eicosanoid, Arachidonic acid (all CYP-heavy subsystems) — ~250–900 rxns per model

Perturbation: for each rxn in the set, force `lb = FORCE * ub_baseline` (initial version) then reformulated to `lb = FORCE_total / n_cyp` (per-rxn distribution of total aggregate flux).

**Why it failed.** Both formulations were catastrophic. Every model went infeasible at the first nonzero FORCE level, regardless of scale (even `per_rxn = 1e-6 mmol/gDW/hr` broke it). Root cause: **many CYP reactions in extracted subnetworks are structurally dead** — the extraction kept the reaction because its GPR passed the transcript threshold, but its substrate isn't produced anywhere else in the network. Forcing ANY positive minimum flux on such a reaction is impossible, and there are dozens of such rxns per model. The sweep collapsed to "any cell has ≥1 dead CYP rxn?" (yes, all) instead of "how much CYP load can this cell handle?"

**What survives from s49.** Diagnostic knowledge: extracted-model CYP sets contain many structurally-dead reactions. Any forced-flux susceptibility analysis on these reactions is invalid. The bio_base numbers from the sweep at FORCE=0 are still valid and reproduce baseline biomass capacity per method x cell.

**Scripts.** `s49_cyp_susceptibility.m` + `s49_cyp_report.R`. Kept on disk as reference for the wrong approach.

---

## 25 — s50: CYP toxicological stress via NADPH-tax (the correct formulation)

**Why this reformulation.** Real CYP catalysis has a mechanistic biochemical cost: `R-H + O2 + NADPH + H+ -> R-OH + H2O + NADP+`. Every unit of CYP flux consumes NADPH and O2 regardless of which specific CYP is running. Instead of forcing arbitrary individual reactions, model the **aggregate biochemical cost** of the CYP system.

**Perturbation.** Add ONE new reaction to each model:
```
CYP_STRESS: NADPH + H+ + 0.5 O2 -> NADP+ + H2O
```
Mass-balanced (N, H, O, charge all check). Sweep `lb(CYP_STRESS)` from 0 to 200 mmol/gDW/hr. At each level, `max biomass` and `max ATP` under the constraint that the cell must run at least `lb` units of CYP-like NADPH oxidation.

**Why this is defensible.**
1. **No structural infeasibility.** The demand is on a single virtual reaction consuming intracellular NADPH + O2 — every viable model has some of both.
2. **Biologically real.** The NADPH+O2 → NADP++H2O half-reaction IS what CYP catalysis consumes; drug hydroxylation is the accompanying product. Modeling the cost side alone captures the metabolic burden without needing to specify which chemical is being detoxified.
3. **Graded dose response.** Cells with strong NADPH-producing capacity (PPP, malic enzyme, IDH1/2) and O2 delivery will handle high stress. Cells at their energetic edge collapse fast. This gives a real susceptibility rank.
4. **Model-agnostic.** Works identically on GIMME, fastCore, bounded, unbounded — one added reaction only.
5. **Cheap.** 60 models × ~10 stress levels × 2 LPs = 1200 LPs = seconds.

**Susceptibility metrics** (per cell x method):
- **biomass-loss AUC**: integral of `(1 - bio_frac)` across the stress sweep — overall susceptibility
- **LD50 analog**: stress level at which biomass_frac crosses 0.5 (linear interpolation between points) — half-max tolerance
- **stress_infeasible**: minimum stress level at which the model becomes LP-infeasible — tolerance ceiling
- **max_loss_bio**: fraction of baseline biomass capacity lost at max stress

**Scripts.** `s50_cyp_stress.m` (MATLAB, parfor across 60 models) + `s50_cyp_stress_report.R` (PDF: tradeoff curves per method, LD50 ranking, susceptibility heatmap, hpf x condition boxplots, ranked tables).

**Outputs:**
- `results/s50_cyp_stress.pdf` — full susceptibility report
- `results/pipeline/s50_cyp_stress_sweep.csv` — long table (method x cell x stress -> bio_frac, atp_frac)
- `results/pipeline/s50_susceptibility_rank.csv` — per-cell susceptibility metrics

**Findings after tuning the sweep.**

Initial sweep with STRESS_LEVELS ≤ 200 mmol/gDW/hr showed no effect — model has enormous baseline NADPH regeneration capacity. Widened to 1×10⁷; found a sharp cliff between 5000 and 20000. Root cause identified: the ub on the CYP_STRESS reaction was set to 1000 by `addReaction`, so any stress > 1000 was capped by MY constraint, not by the model. Raised ub to 1×10⁹ and re-swept. Now graded response with real between-cell variance in the 5000–15000 range.

**Susceptibility findings at stress = 10,000 mmol/gDW/hr NADPH tax (discriminating level):**

- **GIMME_bounded**: 5 late-embryo cells (BL_120, D_96, D_120, LD_96, LD_120) fully infeasible; all 24–72 hpf cells robust. Clean age-based binary split.
- **GIMME** (unbounded): D_120 and LD_120 fully infeasible; D_48, BL_120 mildly affected. Weaker signal because glucose/O2 uptake unrestricted.
- **fastCore_bounded**: graded 0.93–1.00; most susceptible are 72–120 hpf cells (LD_120, D_120, D_72, LD_72, BL_72).
- **fastCore** (unbounded): all cells at 1.000 — no differentiation. Unlimited nutrient inflow lets model absorb any tax.

**Cross-method concordance.** Both GIMME_bounded and fastCore_bounded flag the same age direction (late-embryo susceptibility). GIMME_bounded gives the binary split; fastCore_bounded gives finer gradation. Independent extraction methods reaching consistent conclusions.

**Key methodological insight.** **Nutrient inflow constraints from metabolomics are essential for meaningful toxicological stress simulation.** Unbounded extractions can absorb arbitrary stress by up-regulating uptake; the biology only emerges with the metabolomics-derived bounds. This is a real argument for using the bounded variants when studying stress response, even though we chose fastCore-unbounded as canonical for developmental-axis biology (s46).

**Note on bio_base vs susceptibility.** Baseline biomass capacity does NOT predict susceptibility. GIMME_bounded LD_48 has the lowest bio_base (1.97) yet is fully robust; LD_120 has higher bio_base (3.37) yet fully collapses. Susceptibility structure is orthogonal to raw capacity — it reflects the network's flexibility to redirect flux under stress, not its overall throughput.

---

## 26 — s51: mechanistic overlay linking CYP susceptibility to subsystem envelopes

**Why.** s50 identified WHICH cells are susceptible (late-embryo bounded models). s51 asks WHY: which subsystems' envelope levels predict per-cell susceptibility? This turns a per-cell susceptibility number into a per-subsystem interpretation.

**Approach.** For each method × subsystem: Pearson correlation of `mean subsystem envelope` (from s40, all 15 cells) with `bio_frac at the discriminating stress level` (from s50, all 15 cells). Positive correlation = "high envelope predicts robustness"; negative = "high envelope predicts susceptibility." The discriminating stress level is chosen per method as the level with maximum between-cell variance in bio_frac.

Then cross-reference with s44: does the susceptibility-driving subsystem list overlap with the condition-changing subsystem list? If a subsystem is BOTH condition-plastic and susceptibility-linked, that's the mechanistic connection between light/dark condition and toxicological vulnerability.

**Script.** `s51_cyp_mechanistic_overlay.R`.

**Outputs:**
- `results/s51_cyp_mechanism.pdf` — per-method subsystem ranking, susceptibility × condition-change overlay scatter, per-cell susceptibility bars, top ±12 drivers tables
- `results/pipeline/s51_susceptibility_subsystem_correlations.csv` — full method × subsystem × Pearson r long table
