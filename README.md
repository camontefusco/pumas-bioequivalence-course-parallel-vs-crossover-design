# 📊 Student Project 1  
## Parallel vs Crossover Design Comparison in Bioequivalence Studies

This project compares **parallel** and **crossover** bioequivalence (BE) study designs using real clinical datasets analyzed with **Pumas** and **Bioequivalence.jl**.

The goal is to understand how **study design**, **variability**, and **sample size** affect:
- Bioequivalence outcomes
- Confidence interval precision
- Statistical power

---

## 🔬 Datasets Analyzed

| Design    | Dataset   | Subjects | Endpoint |
|-----------|-----------|----------|----------|
| Parallel  | FSL2015_5 | 60       | AUC      |
| Crossover | SLF2014_5 | 18       | AUC      |

Datasets were loaded from **PharmaDatasets.jl**.

---

## 🧪 Analyses Performed

- Descriptive statistics (mean, median, CV)
- Standard Average Bioequivalence (ABE)
- 90% confidence intervals (80–125%)
- Visual comparison of exposure distributions
- Paired subject analysis (crossover)

---

## 📈 Key Findings

- The **parallel study** demonstrated BE due to low variability.
- The **crossover study failed BE** because of extremely high within-subject variability.
- Crossover designs are **not automatically superior** — variability is decisive.
- Pilot data are essential for design justification.

---

## 📂 Repository Structure

