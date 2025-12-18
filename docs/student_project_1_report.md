# 📊 Student Project 1

## Parallel vs Crossover Design Comparison in Bioequivalence Studies

**Course:** Pumas Bioequivalence Course  
**Author:** Carlos Victor Montefusco Pereira  
**Generated:** *via Julia / Pumas BE analysis*

---

## 1. Project Overview

This project compares **parallel** and **crossover** bioequivalence (BE) study designs using real datasets from published literature. The goal is to understand how **study design**, **variability**, and **sample size** jointly influence bioequivalence outcomes, statistical power, and confidence interval precision.

While crossover designs are often considered more efficient, this project demonstrates that **design efficiency critically depends on variability** and that crossover studies can fail when within-subject variability is extreme.

---

## 2. Learning Objectives

This project addresses the following objectives:

* Compare **parallel vs crossover** BE study designs  
* Perform full BE analysis using `pumas_be`  
* Understand how **variability affects confidence intervals**  
* Interpret BE outcomes under different designs  
* Connect design choice to **power and planning considerations**  
* Visualize exposure distributions and paired subject behavior  

---

## 3. Background Theory

### 3.1 Parallel Design

In a parallel design, each subject receives **only one formulation** (Test or Reference).

**Advantages**

* No carryover effects  
* Simple logistics  
* Suitable for long half-life drugs  

**Disadvantages**

* Requires larger sample sizes  
* Sensitive to between-subject variability  

---

### 3.2 Crossover Design

In a crossover design, each subject receives **both formulations** in separate periods.

**Advantages**

* Subjects act as their own controls  
* Reduced between-subject variability  
* Typically higher efficiency  

**Disadvantages**

* Requires washout  
* Vulnerable to carryover  
* Can fail under high within-subject variability  

---

## 4. Data and Methods

### 4.1 Datasets Used

| Design    | Dataset   | Subjects | Periods | Endpoint |
|-----------|-----------|----------|---------|----------|
| Parallel  | FSL2015_5 | 60       | 1       | AUC      |
| Crossover | SLF2014_5 | 18       | 2       | AUC      |

Both datasets were analyzed using **Standard Average Bioequivalence (ABE)** with 90% confidence intervals and 80–125% acceptance limits.

---

## 5. Results

### 5.1 Descriptive Statistics — Parallel Design (AUC)

| Formulation | Mean  | Median | SD   | CV (%) |
|-------------|-------|--------|------|--------|
| Reference   | 102.9 | 101.9  | 6.27 | ~6.1   |
| Test        | 112.3 | 115.2  | 6.56 | ~5.8   |

**Interpretation**

* Low variability  
* Tight distributions  
* Ideal conditions for a parallel BE study  

Visual inspection confirms:

* Narrow boxplots  
* Compact histograms  
* Stable exposure distributions  

---

### 5.2 Descriptive Statistics — Crossover Design (AUC)

| Formulation | Mean | Median | SD    | CV (%) |
|-------------|------|--------|-------|--------|
| Reference   | 5.23 | 4.91   | 2.67  | ~51    |
| Test        | 8.46 | 3.30   | 17.79 | ~210   |

**Interpretation**

* Extremely high variability  
* Strong skewness and outliers  
* Median differs substantially from mean  
* Test formulation particularly unstable  

Plots show:

* Long-tailed distributions  
* Extreme values dominating the scale  
* High within-subject inconsistency  

---

## 6. Bioequivalence Results

### 6.1 Parallel Design — AUC

* **GMR:** 109.2%  
* **90% CI:** [106.4%, 112.1%]  
* ✅ **Bioequivalence demonstrated**

**Reason**

* Low CV  
* Adequate sample size  
* Narrow confidence interval  

This represents a textbook successful parallel BE study.

---

### 6.2 Crossover Design — AUC

* **GMR:** 91.8%  
* **90% CI:** [55.7%, 151.4%]  
* ❌ **Bioequivalence not demonstrated**

**Reason**

* CV > 100%  
* Small sample size (18 subjects)  
* Extremely wide confidence interval  

⚠️ **Key Insight**

> Crossover designs do not guarantee success. When within-subject variability is extreme, crossover efficiency collapses.

---

## 7. Power and Planning Implications

### 7.1 Crossover Design (Pilot CV ≈ 51%)

* Approximate power ≈ 0  
* **~82 subjects required** to achieve 80% power  

### 7.2 Parallel Design (Pilot CV ≈ 6%)

* Approximate power ≈ 100%  
* **~8 subjects sufficient** for 80% power  

**Conclusion**

* Variability dominates design efficiency  
* A poorly behaved crossover may be less efficient than a stable parallel study  

This naturally motivates:

* **Project 2:** Reference-scaled BE (RSABE)  
* **Project 3:** Power and sample size planning  

---

## 8. Interpretation and Discussion

### 8.1 Design Choice Is Variability-Dependent

While crossover designs are generally preferred, this project demonstrates that:

* High within-subject variability can negate crossover advantages  
* Parallel designs can outperform crossover designs under favorable variability  
* Design choice must be justified using pilot variability estimates  

---

### 8.2 Visual Evidence

The paired crossover plot shows:

* Few subjects driving variability  
* Poor clustering around the identity line  
* Direct visual explanation for CI inflation  

---

## 9. Connection to Other Course Projects

* **Project 2 (RSABE & NTID):**  
  Explains how regulatory agencies handle highly variable crossover data.

* **Project 3 (Power Analysis):**  
  Demonstrates how CV drives sample size requirements.

* **Project 4 (Tmax Analysis):**  
  Shows that not all endpoints are suitable for parametric inference.

Together, these projects illustrate how **design, variability, endpoint choice, and regulation interact**.

---

## 10. Key Takeaways

* Crossover designs are not universally superior  
* Variability is the dominant driver of BE success  
* Low-variability parallel studies can outperform unstable crossover studies  
* Pilot data are essential for design justification  
* Study design is both a **statistical and strategic decision**  

---

## 11. Files Produced

* 📊 **Figures**
  * `figures/crossover_AUC_boxplot.png`
  * `figures/crossover_AUC_hist.png`
  * `figures/parallel_AUC_boxplot.png`
  * `figures/parallel_AUC_hist.png`
  * `figures/crossover_paired_AUC.png`

* 📄 **Outputs**
  * `outputs/student_project_1_outputs.txt`

* 📝 **Questions**
  * `docs/student_project_1_questions.md`

---

## 12. License

This project is shared under **CC BY-SA 4.0**, consistent with course materials.