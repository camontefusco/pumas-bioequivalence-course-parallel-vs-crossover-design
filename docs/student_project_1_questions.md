# Student Project 1 — Questions & What This Repo Answers

This file summarizes the **project questions** and points to where results appear in the repo.

## Task 1 — Data exploration

- **How many subjects are in each study?**  
  See `outputs/student_project_1_outputs.txt` ("N rows / N subjects").
- **What are the key variables?**  
  See `outputs/student_project_1_outputs.txt` ("dataset columns").
- **How do designs differ in organization?**  
  Parallel has one treatment per subject; crossover has multiple periods and uses within-subject comparisons.

## Task 2 — Bioequivalence analysis

- **Compare GMR and 90% CI between designs** (AUC/Cmax when available).  
  See `outputs/student_project_1_outputs.txt` ("Parallel BE" and "Crossover BE").
- **Which design gives tighter CI and why?**  
  Crossover typically yields tighter CIs because it controls between-subject variability.

## Task 3 — Power and sample size

- **What sample size is needed for ~80% power?**  
  See the "planning power" blocks in `outputs/student_project_1_outputs.txt` (didactic simulation-based approximation).
- **Why crossover is more efficient**: within-subject comparisons reduce unexplained variance → fewer subjects needed.

## Task 4 — Practical considerations

- Parallel design avoids carryover/washout but usually needs more participants.
- Crossover reduces sample size but requires washout and more complex conduct/analysis.

## Visualizations

- Distribution plots (hist/box) and paired plots (crossover) are stored in `figures_project1/`.

