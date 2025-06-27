# BetaTREND 

* **BetaTREND** (Trait-Aligned Gradient Diagnostics and GCO Decomposition for Microbiome Beta Diversity): A lightweight R toolkit for trait-aligned gradient analysis of microbiome
β-diversity.*

---

## Overview

BetaTREND tackles a single question:

> **Does a continuous host or environmental trait structure variation in
> microbial community composition?**

The package delivers four complementary steps:

| Stage | Function | Purpose |
|-------|----------|---------|
| **1 Project**    | `gradient_projection()`      | Fit a 1-D gradient (least-squares projection) of a trait onto a PCoA embedding; return projection scores, direction vector, fit statistics, and an optional publication-quality plot. |
| **2 Diagnose**   | `evaluate_consistency_metrics()` | Quantify how well the 1-D projection preserves neighbourhood structure in the full PCoA space using trustworthiness, continuity, and Mantel correlation. |
| **3 Test**       | `cjtt_u_test()`             | Conduct the CJTT-U sign–sign U-statistic trend test for monotone association between beta-dispersion and trait, with analytic / jackknife / bootstrap variance, permutation option, and double FWL residualisation for confounders. |
| **4 Decompose**  | `gco_decomposition()`       | Classify each sample as **G-type** (gradient-driven), **C-type** (cluster-driven), **O-type** (outlier), or *Intermediate* using local neighbour consistency and an outlier filter. |

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("YiqianZhang1004/BetaTREND")
```

---

## Function Outline

* **`gradient_projection()`**  
  *Inputs*: distance matrix (`dist` or square matrix), numeric trait,  
  optional `ndim` (PCoA axes) and `figure`.  
  *Outputs*: projection scores, gradient vector, PCoA coordinates, adj R²,
  global *p*, optional plot, variance explained.

* **`evaluate_consistency_metrics()`**  
  *Inputs*: PCoA matrix, projection scores, neighbourhood size `k`,
  Mantel method.  
  *Outputs*: trustworthiness, continuity, Mantel *r* and *p*.

* **`cjtt_u_test()`**  
  *Inputs*: distance matrix, trait, optional confounder matrix,
  variance method, permutation flag, replicates.  
  *Outputs*: U-statistic, variance estimate, *z*-score, *p*-value.

* **`gco_decomposition()`**  
  *Inputs*: PCoA matrix, projection scores, `k`, thresholds `tau_H`, `tau_L`.  
  *Outputs*: data-frame of consistency scores, outlier flag, class label.

---

## Minimal Reproducible Example

Run the chunk below—results will print automatically.

```r
library(BetaTREND)
library(phyloseq)

set.seed(42)

## 1 ── Simulate data -------------------------------------------------
otu  <- matrix(abs(rnorm(100 * 50, 0, 3)), nrow = 50,
               dimnames = list(paste0("S", 1:50),
                               paste0("OTU", 1:100)))
phy  <- phyloseq(otu_table(otu, taxa_are_rows = FALSE))

trait <- seq(1, 50) + rnorm(50, 0, 2)           # continuous phenotype
dist_bray <- distance(phy, method = "bray")     # Bray–Curtis

## 2 ── Gradient projection ------------------------------------------
gp <- gradient_projection(dist_bray, trait, ndim = 5, figure = TRUE)

gp$adjusted_R2        # model fit
gp$p_value

## 3 ── Consistency diagnostics --------------------------------------
metrics <- evaluate_consistency_metrics(gp$PCoA_coordinates,
                                        gp$projection_scores,
                                        k = 10)

metrics$trustworthiness
metrics$continuity

## 4 ── Trend test ---------------------------------------------------
cjtt <- cjtt_u_test(dist_bray, trait, var_method = "jackknife")

cjtt$p_value

## 5 ── G – C – O decomposition --------------------------------------
gco <- gco_decomposition(gp$PCoA_coordinates,
                         gp$projection_scores,
                         k = 10)

table(gco$class)
```

---

