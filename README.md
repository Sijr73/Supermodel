

# Metabolic-Complexity-Evolvability&nbsp;: Code & Data Companion  
Reproducible pipeline for **“Metabolic complexity increases evolvability”**  
Claus J. Fritzemeier\*, Sajjad Ghaffarinasab\*, *et al.* (2025)

---

## 📑 Overview
This repo contains every script, notebook and raw result used in the manuscript.  
It automatically:

1. **Builds a pan-genome super-model** from 102 curated GSMs  
2. **Calculates environmental distances** via FBA + ARM-LP  
3. **Quantifies adaptability** (added-reaction distributions & power-law fits)  
4. **Analyses *E. coli* subsystems** involved in horizontal gene transfer  
5. **Generates publication-ready figures & tables** (Figs 1–3, S1–S10)

> **TL;DR** &nbsp;`bash run_pipeline.sh` → `results/` contains all figures, CSVs and intermediate models.

---

## 🗁 Repository Layout

```text
.
├── 01_make_supermodel/          # merge GSMs, balance mass, remove EGCs
├── 02_environmental_distance/   # FBA + ARM-LP for wet-lab & random media
├── 03_stats/                    # aggregation, regressions, main-text figures
├── 04_ecoli_subsystems/         # strain-restricted HGT & odds-ratio analysis
├── common/                      # generic biomass reactions, media files, utils
├── results/                     # ✓ auto-generated outputs (ignored by git)
├── environment.yml              # Conda spec (R + sybil, CPLEX, python libs)
├── run_pipeline.sh              # one-stop reproducibility driver
└── LICENSE
