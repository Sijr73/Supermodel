

# Metabolic-Complexity-Evolvability&nbsp;: Code & Data  
Reproducible pipeline for **“Metabolic complexity increases evolvability”**  

---

## 📑 Overview
This repo contains every script, notebook and raw result used in the manuscript.  
It automatically:

1. **Builds a pan-genome super-model** from 102 curated GSMs  
2. **Calculates environmental distances** via FBA + ARM-LP  
3. **Quantifies adaptability** (added-reaction distributions)  
4. **Analyses *E. coli* subsystems** involved in horizontal gene transfer  
5. **Generates publication-ready figures & tables** (Figs 1–3, S1–S10)



---

## 🗁 Repository Layout

```text
.
├── download/                    # download all the GSM models from BiGG database and it is saved in SourceData
├── modelCheck/                  # read and check all the GSMs and preparing for building the Supermodel
├── buildModel/                  # merge GSMs and analyze them
├── convertMedia/                # defining the culture media for all the models (wetlab+random)
├── massBalanceCheck/            # checking the mass balance of the models
├── imbalanceRemoval/            # removing the erroneous mass balances from the models
├── energyGeneratingCycleRemoval/ # removing erroneous energy-generating cycles (Running on HPC)
├── envirDist/                   # FBA + ARM-LP for wet-lab & random media (Running on HPC)
├── DATA/                        # including all the generated results from the analysis
├── paperPlots/                  # aggregation, regressions, main-text figures
├── makeFile.sh                  # one-stop reproducibility driver
└── LICENSE
```
## Contact
If you have any questions or run into issues, please open an issue in the repository or contact the repository maintainer.

## Citation
If you build upon this code or data, please cite:
