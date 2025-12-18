# 📊 Parallel vs Crossover Design Comparison in Bioequivalence Studies

![Pumas BE Course](https://img.shields.io/badge/Pumas-Bioequivalence-blue)
![Julia](https://img.shields.io/badge/Julia-1.9+-purple)
![License](https://img.shields.io/badge/License-CC%20BY--SA%204.0-green)
![Status](https://img.shields.io/badge/Status-Completed-success)

---

## 🧪 Project Overview

This repository contains **Student Project 1** from the **Pumas Bioequivalence Course**, focusing on a **comparative analysis of parallel and crossover study designs** using real bioequivalence datasets.

The project demonstrates how **study design, variability, and sample size** interact to influence:

- Bioequivalence outcomes  
- Confidence interval width  
- Statistical power  
- Regulatory success or failure  

Although crossover designs are often considered more efficient, this analysis shows that **design efficiency is highly dependent on variability**, and that crossover studies can fail under extreme within-subject variability.

---

## 🎯 Learning Objectives

By completing this project, the following objectives were addressed:

- Understand the **fundamental differences** between parallel and crossover designs
- Perform full **bioequivalence analysis** using `pumas_be`
- Evaluate the impact of **between- and within-subject variability**
- Compare **confidence interval precision** across designs
- Connect design choice to **power and planning considerations**
- Communicate results using **clear visualizations**

---

## 📦 Datasets Analyzed

| Design    | Dataset     | Subjects | Periods | Endpoint |
|---------- |------------ |----------|---------|----------|
| Parallel  | FSL2015_5   | 60       | 1       | AUC      |
| Crossover | SLF2014_5   | 18       | 2       | AUC      |

All datasets were loaded from **`PharmaDatasets.jl`** and analyzed using **Standard Average Bioequivalence (ABE)** with 90% confidence intervals and 80–125% acceptance limits.

---

## 🔬 Analyses Performed

### Task 1 — Data Exploration
- Descriptive statistics (mean, median, SD, CV%)
- Distribution comparison between designs
- Visualization of exposure patterns

### Task 2 — Bioequivalence Analysis
- GMR and 90% confidence interval estimation
- Pass/fail determination under standard ABE
- Interpretation of CI width and variability effects

### Task 3 — Power & Planning Implications
- Approximate planning power based on pilot CV
- Sample size requirements for 80% power
- Design efficiency comparison

---

## 📈 Visualizations

Figures are stored in the `figures/` directory and include:

- Parallel AUC histograms and boxplots
- Crossover AUC histograms and boxplots
- Paired subject plots (Test vs Reference) for crossover design

These plots visually demonstrate:
- Tight distributions under low variability
- CI inflation under high within-subject variability
- Loss of crossover efficiency when CV is extreme

---

## 🧠 Key Findings

- **Parallel design** succeeded due to low variability (CV ≈ 6%)
- **Crossover design failed** despite theoretical efficiency, due to extreme CV (>100%)
- Variability dominates design efficiency more than design type alone
- Crossover designs are not universally superior
- Pilot variability estimates are essential for design justification

---

## 🔗 Connection to Other Course Projects

This project connects directly to:

- **Project 2 — Reference Scaling & NTID**  
  Explains regulatory alternatives when variability is high.

- **Project 3 — Power & Sample Size Determination**  
  Quantifies how CV drives sample size requirements.

- **Project 4 — Tmax Nonparametric Analysis**  
  Highlights that not all endpoints are suitable for parametric inference.

Together, these projects form a **coherent bioequivalence analytics portfolio**.

---

## ▶️ How to Run

```bash
julia --project=. scripts/student_project_1_run_and_save.jl
```
This will:

- Load datasets  
- Run bioequivalence (BE) analyses  
- Generate figures  
- Write outputs to `outputs/`

---

## 📄 Reports & Documentation

- 📝 **Full analysis report:**  
  `docs/student_project_1_report.md`

- ❓ **Project questions:**  
  `docs/student_project_1_questions.md`

- 📊 **Execution outputs:**  
  `outputs/student_project_1_outputs.txt`

---

## 🧑‍🎓 Author

**Carlos Victor Montefusco Pereira**  
Pumas Bioequivalence Course
