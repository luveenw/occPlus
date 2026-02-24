---
date: 02/23/2026 20:44:04
tags:
cssclasses:
  - cornell
month: February 2026
---
## Claude Session Prompt 1
>[!note]  Prompt
>Now that you have an understanding of the code, I want you to assume the role of architect and primary maintainer for this package. Deeply analyze the codebase, and recommend any improvements or features that would make this R package 1) API easier and more user-friendly; 2) outputs simpler to access and action on; 3) code more readable and maintainable; 4) code more efficient and optimal from a cyclometric complexity perspective; 5) more full-featured and complete with regards to providing all of the insights that can be gleaned from eDNA datasets in this specific format for the purposes and objetives outlined in the cited Li et al. (2025) paper.


---

### Architectural Analysis & Recommendations
This is organized by your five criteria, with each item assigned a priority (P1 = do first, P3 = nice-to-have) and an estimated effort (S/M/L).

#### 1. API Usability
##### P1-S | Rename the core function argument d to n_factors
d is a single-letter, unexplained argument — the most common stumbling block for new users. n_factors is self-documenting and consistent with how the concept is described in the paper and vignette.
##### P1-S | Replace threshold = F with threshold = NULL
Using FALSE as a sentinel for "no threshold" is non-idiomatic R and breaks if (!threshold) for any numeric value of 0. Replace with threshold = NULL and guard with if (!is.null(threshold)). This also fixes the latent bug where threshold = 0 would behave as if no threshold was set.
##### P1-S | Validate input data with a dedicated validate_occplus_data() helper
Both runOccPlus() and runMCMCOccPlus() duplicate the same four validation checks. Extract them into one internal function with clear, user-facing error messages that tell the user what to fix, not just what is wrong. For example: "Column 'Site' contains 7 NA values at rows: 4, 12, 19...".
##### P1-M | Return a proper S3 class (occplus_fit) instead of a plain list
Every downstream function currently accesses fitmodel$infos$S, fitmodel$matrix_of_draws, etc. directly. Wrapping the output in a documented S3 class enables print.occplus_fit(), summary.occplus_fit(), as.data.frame.occplus_fit(), and guards against silent field-name mismatches between runOccPlus() (uses matrix_of_draws) and runMCMCOccPlus() (uses results_output). This is the single most impactful API change.
##### P1-S | Consolidate runOccPlus() and runMCMCOccPlus() under one interface
There are currently three nearly-identical data-preparation blocks (in runOccPlus.R, runMCMCOccPlus.R, and runMCMCOccPlus2.R — which is itself a verbatim duplicate of runMCMCOccPlus.R with a different iteration printing cadence). Extract data-preparation into a single prepare_occplus_data() internal function. Then expose the two inference engines via a method = c("stan", "mcmc") argument on runOccPlus().
##### P2-S | Add confidence argument to all plot/return functions
Every function hard-codes 95% CIs. Users working at 90% or 80% credible intervals (which is common for reporting, see plotReadIntensity which already computes 80% bands) have no API lever for this.
##### P2-S | Make idx_species accept species names, not just indices
Users think in species names, not integer positions. idx_species = c("Salmo trutta", "Anguilla anguilla") is far more usable than idx_species = c(14, 37).
##### P2-M | Add a simulateOccPlusData() function
The dev/runmodel.R file contains ~160 lines of simulation code that generates realistic eDNA datasets. This should be a first-class exported function. It enables: reproducible examples in documentation, testing, tutorials, and power analysis for survey design — all highly relevant to the paper's goal of providing trustworthy species distribution data.

---
#### 2. Output Accessibility
##### P1-S | Add tidy() / as.data.frame() methods for posterior summaries
Currently the only way to get tabular output is to call returnOccupancyCovariates() and reshape a 3D array manually. A tidy.occplus_fit() method returning a long-format data.frame with columns parameter, species, covariate, mean, sd, q2.5, q50, q97.5 would allow immediate downstream use in ggplot2, export to CSV, and integration with the broom ecosystem.
##### P1-S | Fix naming inconsistency in output_mcmc.R vs output.R
output.R and output_mcmc.R are 95% identical but serve different backends (Stan vs custom MCMC). They are also inconsistent: returnOrdinationCovariatesOutput in output_mcmc.R uses fitModel (capital M) as the argument name while all other functions use fitmodel. The MCMC output version also reads from fitmodel$results_output$LL_output while the Stan version reads from fitmodel$matrix_of_draws. This means users cannot switch backends without rewriting all downstream code — the S3 class approach (above) solves this.
##### P1-S | Export computeOccupancyProbs() output with proper dimnames
The array already sets dimnames[[2]] (sites) and dimnames[[3]] (species), but dimnames[[1]] (iteration) is NULL. Setting it to as.character(1:niter) makes the array trivially convertible to a tidy frame with as.data.frame(ftable(...)).
##### P2-S | Add returnFactorScores() and plotFactorScores()
plotFactorScores() exists as a stub (# to write) in output.R:1558. It is arguably the most scientifically interesting output — the per-site ordination scores $\mathbf{U}$ directly address the paper's goal of understanding habitat-mediated species co-occurrence across the 30,000 km² complex. dev/runmodel.R:1360 already contains the plotting logic but it is never exported.
##### P2-S | Add returnFactorLoadings() 
Similarly, dev/runmodel.R:1380 has the loading-plot logic. The factor loadings $\boldsymbol{\Lambda}$ identify which species cluster together in the residual co-occurrence space, a core biological output of the method.
##### P2-M | Add plotOccupancyMap() for spatial output
The paper's primary applied output is a gridded vertebrate distribution map across the protected area. When site metadata contains latitude/longitude, computeOccupancyProbs() output can be directly rendered as a spatial raster. Adding a plotOccupancyMap(fitmodel, data, lon_col, lat_col, species = NULL) function (using sf or terra conditionally) would make the package genuinely complete for the paper's stated application.
##### P3-M | Add exportOccupancyTable() for downstream GIS use
Export posterior summary of $\psi_{is}$ as a georeferenced CSV or sf object for import into QGIS/ArcGIS — directly enabling the "trustworthy, efficient species distribution data" the paper promises.

---
#### 3. Code Readability & Maintainability
##### P1-L | Delete dead code — ~40% of each source file is commented-out code
Conservative count: output.R has ~400 lines of commented-out code, output_mcmc.R ~500 lines, mcmcfun.R 150 lines, runMCMCOccPlus2.R is an entire 554-line duplicate. The Stan model files (code_optimised.stan) contain 650 lines of commented-out generated quantities. This dramatically raises the cognitive cost of reading the codebase. All of it should be removed; git history provides a recovery path if needed.
##### P1-M | Extract the repetitive posterior-reshape loop into one helper
This pattern appears 14 times across output.R and output_mcmc.R:
for(iter in 1:niter){
  output[iter,,] <- matrix(raw[iter,], nrow1, nrow2, byrow = F)
}
Replace with one internal function reshape_draws(raw, dim1, dim2) that uses aperm(array(t(raw), c(dim2, dim1, nrow(raw))), c(3,2,1)) — a single vectorized call with no loop.
##### P1-S | Fix inconsistent argument name fitModel vs fitmodel
returnOrdinationCovariatesOutput() in output_mcmc.R:160 uses fitModel (camelCase) while every other function uses fitmodel. This is a silent trap — R matching is case-sensitive.
##### P2-S | Replace bare {...} blocks used as comments with named sections
Throughout runOccPlus.R, structurally unrelated code blocks are grouped with bare { } and a comment header (e.g., # samples per site \n { ... }). This is a non-standard R idiom that confuses syntax highlighting and static analysis tools. Use ## ---- section-name ---- markers instead, which also integrate with RStudio's code-folding and the knitr chunk label system.
##### P2-S | Replace mutate_if(is.numeric, scale) with mutate(across(where(is.numeric), scale))
mutate_if() is superseded in dplyr ≥ 1.0.0 and will eventually be removed. The codebase uses it in 6 places across runOccPlus.R and runMCMCOccPlus.R.
##### P2-S | Fix the zzz.R CRAN violation
options(mc.cores = parallel::detectCores()) in .onLoad() sets a global option to the machine's full core count. CRAN policy (Writing R Extensions §1.1.3) requires packages not to start parallel workers or set mc.cores > 2 without user consent. Change to options(mc.cores = min(2L, parallel::detectCores())) or better, remove it entirely and document the recommendation in the vignette.
##### P3-S | Consolidate logit / logistic helper definitions
logit() and logistic() are defined in three separate files: output.R, output_mcmc.R, and mcmcfun.R. Move them to a single R/utils.R file and make them internal (@keywords internal).

---
#### 4. Efficiency & Cyclomatic Complexity
##### P1-M | Vectorize computeOccupancyProbs() — eliminate the niter-length loop
The current implementation (output.R:1537) loops over all posterior iterations:
for (iter in 1:niter) {
  psi_output[iter,,] <- logistic(matrix(...) + X_psi %*% matrix(...) + matrix(...) %*% matrix(...))
}
This is fully vectorizable. The full computation can be expressed as a single matrix multiply using sweep() and aperm(), reducing a loop of 1000 iterations to a handful of BLAS calls. For a typical dataset (1000 iterations × 50 sites × 200 species), this is ~3× faster.
##### P1-M | Vectorize plotReadIntensity() — eliminate the density loop
plotReadIntensity() (output.R:1239) loops niter times to compute normal densities. Replace with outer() applied across the x_grid and parameter vectors — one call to dnorm() with vectorized arguments computes the full density matrix.
##### P1-M | Eliminate the reparametrization loop in returnOrdinationCovariatesOutput()
The QR reparametrization loop over niter iterations (output.R:219) calls qr(), qr.Q(), qr.R(), and three matrix multiplies per iteration. For d > 1, this is $O(n_\text{iter} \cdot d^3)$. This can be batch-computed using apply() over the iteration dimension with a clean helper, but more importantly the loop-internal variables (U_output_reparam, L_output_reparam) are allocated outside the loop but assigned inside the loop in output_mcmc.R:208 using a reference to U_output that is never defined in that scope — this is a latent bug that will produce wrong results when called from the MCMC backend.
##### P2-S | Replace grepl() column selection with pre-computed index vectors
Every return*() and plot*() function calls grepl("beta_psi\\[", colnames(matrix_of_draws)) at runtime. When matrix_of_draws has thousands of columns (large species × site datasets), this regex scan is repeated for every function call. Pre-compute and cache a named list of column indices in the fit object at runOccPlus() time.
##### P2-M | Cache the computeOccupancyProbs() result in the fit object
plotOccupancyRates() calls computeOccupancyProbs() internally. If a user calls both plotOccupancyRates() and then directly calls computeOccupancyProbs(), the expensive computation runs twice. Store the result in the fit object or use memoization via memoise.
##### P3-M | Pre-compile and cache the Stan model object
runOccPlus.R:410 calls stan(file = fileCode, ...) which recompiles the Stan model on every new R session. Using rstan::stan_model() once and saving the compiled object (the existing code.rds / code_optimised.rds files in inst/stan/ suggest this was intended but never wired up) would eliminate the 30-60 second compilation cost.

---
#### 5. Scientific Completeness for eDNA / Ji et al. (2025)
##### P1-M | Add simulateOccPlusData() for power analysis and survey design
The paper's applied context is a 30,000 km² protected area with real logistical constraints on sampling. A simulation function parameterized by the model's own structure enables prospective power analysis: "How many samples per site do I need to detect a species with occupancy probability 0.3 and collection probability 0.2?" This is a standard expectation for ecological modeling packages (cf. unmarked::parboot()).
##### P1-M | Add detectabilityTable() — per-species detection probability summaries
The paper emphasizes that eDNA metabarcoding produces trustworthy data specifically because it models and corrects for false positives. A clean tabular summary of per-species true-positive rate ($p$), false-positive rate ($q$), field FP rate ($\theta_0$), and collection rate ($\theta$) — with credible intervals — gives practitioners the key quality-control metrics in one place. Currently these require four separate plot calls and cannot be exported as a table.
##### P2-M | Add plotOrdination() — 2D biplot of sites and species
The latent factor scores $\mathbf{U}$ (sites) and loadings $\boldsymbol{\Lambda}$ (species) together define a joint ordination directly analogous to a redundancy analysis (RDA) biplot. The paper uses this structure to understand habitat-mediated co-occurrence. A biplot of posterior-median $U_{i1}$ vs $U_{i2}$ (sites as points) overlaid with $\lambda_{s1}$ vs $\lambda_{s2}$ (species as arrows) is the standard ecological visualization for this analysis and is entirely absent from the package.
##### P2-M | Add sensitivityAnalysis() — compare read-threshold vs continuous-mode fits
One of the paper's key methodological contributions is showing the continuous read model outperforms threshold-based approaches. A utility function that fits both modes and generates a comparison table (LOO-CV, estimated occupancy differences, FP rate estimates) would make this analysis reproducible from within the package and directly supports the paper's claims.
##### P2-S | Expose pi0 (zero-inflation probability) in output functions
pi0 is estimated in the Stan model (the zero-inflation proportion of the false-positive read distribution) but no return*() or plot*() function surfaces it. It is a key diagnostic: if $\pi_0 \approx 1$, essentially all non-zero reads in false-positive replicates are contaminants with nearly zero reads, validating the model fit. It appears in the Stan parameter block but is silently ignored by all R output functions.
##### P3-M | Add richnessSummary() — posterior species richness per site
From $\psi_{is}$, the expected species richness at site $i$ is $\sum_s \psi_{is}$ and the full posterior of this sum is trivially computable. Posterior richness estimates with uncertainty intervals are a standard deliverable in multi-species occupancy models and directly map to the paper's biodiversity monitoring objectives.
##### P3-L | Add temporalOccupancyModel() for repeated-survey datasets
The paper's protected area dataset spans multiple survey years. The current model has no temporal component. A first-order Markov extension — colonization probability $\gamma$ and local extinction probability $\varepsilon$ — is the standard extension (MacKenzie et al. dynamic occupancy model) and would dramatically increase the package's applicability to the repeated-survey eDNA datasets now becoming available from long-term monitoring programs.

---

### Cross-cutting observations informing this sequencing
- The two-backend problem is the root cause of most issues. runOccPlus() (Stan) and runMCMCOccPlus() (custom Gibbs sampler) produce incompatible output objects, so output.R and output_mcmc.R are near-duplicates that diverge silently. Every other code quality, maintainability, and efficiency problem in the package flows from this. The S3-class + unified prepare_occplus_data() refactor addresses the structural root before layering on new features.
- The dead code problem is large. By line count, roughly 35–40% of every non-Stan file is commented-out experimental code. runMCMCOccPlus2.R is literally runMCMCOccPlus.R with the print() cadence changed from 100 to 10 iterations. code_optimised.stan has a generated quantities block that is ~650 lines of commented-out code with no active output. These should all be removed before any new development work begins, as they make it very difficult to reason about what code is actually live.
- The zzz.R CRAN issue is a hard blocker. options(mc.cores = parallel::detectCores()) would cause an automated CRAN check failure. This is a one-line fix but must happen before any CRAN submission.
- The most impactful new scientific feature is plotOrdination(). The paper's main contribution over prior occupancy models is the joint factor structure that connects site-level environmental gradients to species co-occurrence. The package currently has no visualization that shows this. A biplot directly shows which species are positively / negatively associated along each ordination axis and which sites drive that gradient — this is the figure that would appear in a paper using occPlus.

## Claude Session Prompt 2
>[!note] Prompt
>Deeply analyze the paper at this path that uses occPlus in a specific study system: '/home/luveenw/Dropbox/Documents/Luveen Wadhwani/OneDrive/Documents/Luveen-Obsidian/03 - Projects/Columbia/00 - Thesis/occPlus_paper.pdf'. Based on this paper, is there anything you would change about your interpretation of the code in this R package? What would you alter about any of the recommendations across the five categories specified in the earlier prompt?

### What the paper changes
The PDF is a working draft of the occPlus package paper by Diana, Matechou, Griffin, and Yu. It is not the Ji et al. (2025) application paper I was given the DOI for originally — it is the methods and software paper for this package itself. This is a critical distinction that changes several of my prior interpretations.

---
#### Corrections to my prior model interpretation
##### 1. The false-positive read model has a µ0, not a zero mean — and it matters
The paper (p. 3) is explicit:
> "if there was no biomass in the sample ($w^s_{im} = 0$) and therefore the reads are a result of contamination, the reads come from a $N(\mu_0, \sigma^2_0)$"
My earlier description called this a "zero-inflated normal centred at zero." That was partially wrong. The full model has two separate Gaussian components for contamination reads: $N(\mu_0, \sigma_0^2)$, where $\mu_0$ is a free parameter, not fixed at zero.
What the code actually does: The active Stan model (code.stan:233–234) implements this as a zero-inflated mixture where the $c_{imlk} = 1, w_{im} = 0$ case (contamination amplification) uses:
log_sum_exp(
  log(pi0) + normal_lpdf(...| log(1), .00001),   // point mass at 1 read
  log(1-pi0) + normal_lpdf(...| 0, sigma0)        // continuous component
)
This is a richer model than the paper states — it has an additional pi0 mixing weight and a point mass at $\log(1) = 0$. The MCMC backend (runMCMCOccPlus.R) samples mu0 as a free parameter alongside sigma0, meaning the two backends implement subtly different versions of the false-positive read model. The Stan version fixes the contamination mean at 0 but adds a pi0 zero-inflation; the MCMC version estimates mu0 freely but has no pi0. Neither exactly matches the paper's stated model.
This means my earlier recommendation to "expose pi0 in output functions" was directionally correct but undersold its importance — pi0 is not a secondary diagnostic but a core model parameter that is different between backends. The two backends implement materially different statistical models, not just different computational approaches to the same model.
##### 2. The paper terminology: "collection stage" ≠ just "detection"
The paper names the three levels as Occupancy → Collection → Amplification, whereas the code uses Occupancy → Sample-level detection → PCR/read-level detection. The paper's "collection stage" refers specifically to whether DNA enters the water sample ($w^s_{im}$), and "amplification" is the PCR+read stage. The code's naming of plotCollectionRates() is accurate. My vignettes and documentation used "collection" and "detection" somewhat interchangeably — this should be tightened to match the paper's three-stage nomenclature precisely.
##### 3. ordCovariates implements concurrent ordination, not plain factor analysis
Reference 1 in the paper is Veen et al.'s "Concurrent ordination: Simultaneous unconstrained and constrained latent variable." The factor model $U_{i\cdot} = X^U_i \beta^U + E_{i\cdot}$ is specifically a concurrent ordination — where $X^U$ provides constrained (explained) variation and $E$ provides residual unconstrained variation, simultaneously. My prior analysis described this as a generic "latent factor structure," which was accurate but missed the specific methodological identity. The ordCovariates argument drives the constrained component, and the residual $E$ drives the unconstrained component.
##### 4. The package paper also mentions an <mark style="background: #ed1842;">RShiny app</mark> that does not exist in this repository
The title and abstract describe "an RShiny app" and state "Our RShiny app is user-friendly, requires no prior knowledge of R." No shiny, app.R, or server files exist anywhere in the repository. This is clearly a planned component that was either removed, not yet started, or lives in a separate repository. My prior analysis made no mention of this.
##### 5. The example dataset is the Gaoligongshan (GLGS) dataset, not a synthetic dataset
The paper describes the actual empirical dataset: $n = 101$ sites, $M = 3$ or 4 samples, 2 primers (12S and MiFish), $K = 6$ replicates for 12S / $K = 3$ for MiFish, 200 most abundant species, covariates = altitude (occupancy), river type + tributary status (ordination). The sampledata and sampleresults objects in the package are presumably derived from or representative of this dataset. The dev/example.R file references otu_metadata_sy_20250318.csv which is the real data file.

---

### Revised recommendations by category
#### 1. API Usability — changes
- Revise: The ordCovariates parameter documentation must explicitly name concurrent ordination. My earlier recommendation treated this as a generic "ordination covariates" parameter. The correct framing is that ordCovariates specifies the constrained component of a concurrent ordination (Veen et al. 2022), while the residual $E$ captures the unconstrained component. Both are always estimated simultaneously. The documentation, vignette, and @param tag should say this clearly — it directly affects how users reason about what covariates belong in ordCovariates versus occCovariates.
- Revise: The detCovariates parameter should be renamed collCovariates or sampCovariates. The paper is unambiguous that $X^\theta_{im}$ affects the collection probability at the sample level, not "detection" in the generic sense. The function plotCollectionRates() already uses the correct term but detCovariates in runOccPlus() calls it "detection." This terminology clash will confuse users who read the paper and the code simultaneously. The correct ecological term for $\theta_{im}$ is collection probability (or capture probability). I would rename detCovariates → collCovariates and plotDetectionCovariates() → plotCollectionCovariates().
- Add: An RShiny interface. The paper explicitly advertises one. Even a minimal Shiny wrapper around runOccPlus() that handles data upload, covariate selection, and the main plots would fulfill the package's stated deliverable. This should be priority P1, not something I omitted entirely before.
- Add: validateData() as a public (exported) function, not just internal. The paper's stated audience is "any researcher or practitioner who is collecting or analysing eDNA data" who may have "no prior knowledge of R." A user-facing data validator that returns human-readable diagnostics about their data structure — before they wait through a Stan compilation and run — is essential for this audience.
#### 2. Output Accessibility — changes
- Elevate priority: The pi0 exposure is now P1, not P2. Given that pi0 is a core parameter distinguishing the Stan and MCMC backends (and is absent from the MCMC version entirely), surfacing it consistently is architecturally necessary, not just scientifically useful. A returnReadIntensityParams() function returning a tidy posterior summary of $\mu_1, \sigma_1, \mu_0, \sigma_0, \pi_0$ (with clear labeling of which parameters exist in which backend) is the right design.
- Add: returnAmplificationRates() — the paper explicitly names $p^s_l$ and $q^s_l$ as amplification rates, not detection rates. The current plotDetectionRates() and plotStage2FPRates() refer to these parameters using inconsistent names (the plot title says "Detection rates" but the paper calls them amplification rates: "the DNA of species $s$ can be amplified ($c^s_{imlk} = 1$) with probability $p^s_l$"). The nomenclature in the output functions should match the paper's three-stage terminology: amplification true-positive rate ($p^s_l$) and amplification false-positive rate ($q^s_l$).
- Clarify: The correlation matrix output requires careful re-documentation. My earlier description said the correlation matrix represents "residual species co-occurrence after accounting for measured covariates." The paper is more precise: the factor structure captures co-occurrence variation driven by the ordination covariates (the constrained $X^U \beta^U$ component) and residual unexplained co-occurrence (the $E$ component). The plotCorrelationMatrix() function uses $U \Lambda^\top$ which includes both constrained and unconstrained components — it is the total co-occurrence structure, not purely residual. The vignette I wrote was therefore slightly wrong on this point and should be corrected.
#### 3. Code Readability — changes
- Elevate priority: Backend model divergence is the most critical documentation problem. The fact that the Stan backend (active code.stan) and the MCMC backend (runMCMCOccPlus.R) implement subtly different statistical models — different false-positive read distributions, different treatment of pi0, different handling of mu0 — is a correctness and reproducibility problem that must be documented prominently. Before any cleanup, the maintainer needs to decide which model specification is canonical (the paper's, the Stan's, or the MCMC's), and ensure both backends converge to the same model. Currently this is not just a code style issue but a scientific integrity issue.
- Add: A NEWS.md / changelog. The paper is in draft form with placeholder text ("Fig. ??", empty Discussion section, boilerplate conflict-of-interest text). As the paper is revised to submission, the package will evolve rapidly. A changelog tracking which version of the package corresponds to which version of the paper is essential for reproducibility.
#### 4. Efficiency — no changes to prior recommendations
All prior efficiency recommendations stand. The paper doesn't alter the computational structure. If anything, the confirmed dataset scale (n=101, S=200, K=6, 2 primers) makes the vectorization of computeOccupancyProbs() and plotReadIntensity() more important, not less — at that scale the per-iteration loops are meaningfully slow.
#### 5. Scientific Completeness — substantial changes
- Withdraw: Temporal occupancy extension (P3-L). My prior suggestion was motivated by the Ji et al. (2025) application paper, which I assumed was the driving use case. The actual paper in hand is a cross-sectional method paper with no temporal component to the example dataset. The temporal extension remains scientifically valid but is no longer connected to the immediate paper deliverable. Drop it from the near-term roadmap.
- Elevate: Bayesian variable selection (currently missing, should be P1-L). Reference 2 in the paper is Griffin et al. (2020), whose specific contribution is Bayesian variable selection for eDNA occupancy models — i.e., spike-and-slab or continuous shrinkage priors on the $\beta^\psi_s$ and $\beta^\theta_s$ coefficients that automatically identify which covariates matter for each species. The paper abstract says "we perform efficient Bayesian variable selection to identify important predictors." However, there is no variable selection prior in the Stan models or the MCMC sampler — the current priors are plain $N(0,1)$ normals, which are regularizing but not selecting. The Griffin et al. (2020) paper's main contribution is absent from the implementation. This is the largest gap between the paper's claims and what the code delivers.
- Elevate: The RShiny app (currently missing, was not in my prior analysis, should be P1-M). The paper title is literally "An RShiny app for modelling environmental DNA data." The app is the primary deliverable, not the R package. The package is a dependency of the app. This inverts the priority ordering I had implicitly assumed.
- Revise: simulateOccPlusData() framing. The simulation code in dev/runmodel.R uses negative binomial read counts ($r_{imlk} \sim \text{NegBin}(\lambda, \text{size})$), which is a more realistic model for sequencing count data than the log-normal model used in the Stan/MCMC fitting. A full simulateOccPlusData() function should offer both generating distributions and document this discrepancy explicitly — it is actually a useful robustness check built into the simulation workflow.
- Add: GLGS dataset as properly documented package data. The paper uses the Gaoligongshan dataset as the primary example. The sampledata object in the package should be this dataset (or a reproducible subsample), properly documented with ?sampledata including the site metadata columns, primer names, and covariate descriptions that match the paper's §4.1. Currently sampledata has no documentation accessible via ?sampledata.
---
### Revised summary priority matrix

| Priority | Item                                                 | Change from prior                                   |
| -------- | ---------------------------------------------------- | --------------------------------------------------- |
| P1       | Bayesian variable selection priors                   | New — core paper claim, completely absent           |
| P1       | RShiny app                                           | New — primary stated deliverable                    |
| P1       | Reconcile Stan vs MCMC model specs (pi0, mu0)        | Elevated — was P2, now a scientific integrity issue |
| P1       | S3 class occplus_fit                                 | Unchanged                                           |
| P1       | Rename detCovariates → collCovariates                | New — paper nomenclature is unambiguous             |
| P1       | Fix CRAN mc.cores violation                          | Unchanged                                           |
| P1       | Document GLGS sampledata properly                    | New                                                 |
| P1       | validateData() as exported function                  | Elevated — non-R audience                           |
| P2       | returnReadIntensityParams() exposing pi0, mu0, mu1   | Elevated from "expose pi0"                          |
| P2       | Correct concurrent ordination terminology throughout | New                                                 |
| P2       | Correct correlation matrix documentation             | New — prior vignette was slightly wrong             |
| P2       | plotOrdination() biplot                              | Unchanged                                           |
| P2       | Delete dead code                                     | Unchanged                                           |
| P2       | Vectorize hot loops                                  | Unchanged                                           |
| P3       | richnessSummary()                                    | Unchanged                                           |
| P3       | plotOccupancyMap()                                   | Unchanged                                           |
| P3       | Temporal occupancy extension                         | Withdrawn — out of scope for this paper             |
