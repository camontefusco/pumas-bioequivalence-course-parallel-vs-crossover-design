# scripts/student_project_1_run_and_save.jl
# Run from terminal (NOT inside the Julia REPL):
#   julia --project=. scripts/student_project_1_run_and_save.jl
#
# Student Project 1 — Parallel vs Crossover Design Comparison
# Datasets (PharmaDatasets):
#   - Parallel:  FSL2015_5
#   - Crossover: SLF2014_5
#
# Outputs:
#   - outputs/student_project_1_outputs.txt
#   - docs/student_project_1_report.md
#   - docs/student_project_1_questions.md
#   - figures_project1/*.png

using Pkg
using Dates
using Printf
using Statistics
using DataFrames
using Random

# Plotting
try
    using StatsPlots
catch
    @warn "StatsPlots not found. Installing..."
    Pkg.add("StatsPlots")
    using StatsPlots
end

# Pumas / BE / datasets
using Pumas
using PharmaDatasets
try
    using Bioequivalence
catch
    @warn "Bioequivalence not found. Installing..."
    Pkg.add("Bioequivalence")
    using Bioequivalence
end

# ----------------------------
# Paths (repo-friendly)
# ----------------------------
ROOT = dirname(@__DIR__)
OUTD = joinpath(ROOT, "outputs")
DOCD = joinpath(ROOT, "docs")
FIGD = joinpath(ROOT, "figures_project1")
mkpath.(String[OUTD, DOCD, FIGD])

OUTTXT = joinpath(OUTD, "student_project_1_outputs.txt")
OUTMD  = joinpath(DOCD, "student_project_1_report.md")
QMD    = joinpath(DOCD, "student_project_1_questions.md")

# ----------------------------
# Utilities
# ----------------------------
function save_block(io::IO, title::AbstractString, obj)
    println(io, "\n", "="^110)
    println(io, title)
    println(io, "="^110, "\n")
    show(io, MIME"text/plain"(), obj)
    println(io)
end

function try_run(io::IO, title::AbstractString, f::Function)
    try
        res = f()
        save_block(io, title, res)
        return (true, res)
    catch err
        println(io, "\n", "="^110)
        println(io, title, "  [ERROR]")
        println(io, "="^110, "\n")
        showerror(io, err, catch_backtrace())
        println(io)
        return (false, nothing)
    end
end

has_col(df::AbstractDataFrame, col::Symbol) = (String(col) in names(df))

function find_dataset(substr::AbstractString)
    ds = PharmaDatasets.datasets()
    hits = filter(x -> occursin(substr, x), ds)
    return isempty(hits) ? nothing : first(hits)
end

# ----------------------------
# Formulation derivation
# ----------------------------
"""
Ensure a :formulation column exists.

Supports:
- If :formulation exists -> keep it
- If :treatment exists -> use it
- If :trt exists -> use it
- If crossover with :sequence + :period -> derive from sequence[period]
Otherwise throws.
"""
function ensure_formulation!(df::DataFrame)
    if has_col(df, :formulation)
        df.formulation = String.(df.formulation)
        return df
    end
    if has_col(df, :treatment)
        df.formulation = String.(df.treatment)
        return df
    end
    if has_col(df, :trt)
        df.formulation = String.(df.trt)
        return df
    end
    if has_col(df, :sequence) && has_col(df, :period)
        seqs = String.(df.sequence)
        pers = Int.(df.period)
        f = Vector{String}(undef, nrow(df))
        for i in eachindex(seqs)
            s = seqs[i]
            p = pers[i]
            f[i] = (1 <= p <= lastindex(s)) ? string(s[p]) : "?"
        end
        df.formulation = f
        return df
    end
    error("Could not infer formulation. Columns: $(names(df))")
end

# ----------------------------
# Descriptives (mean/median/CV%)
# ----------------------------
cv_percent(x) = (mean(x) == 0 ? missing : 100 * std(x) / mean(x))

function summarize_endpoint(df::DataFrame, endpoint::Symbol)
    if !has_col(df, endpoint)
        return "Endpoint $(endpoint) not found."
    end
    g = groupby(df, :formulation)
    out = combine(g,
        endpoint => (x -> length(collect(skipmissing(x)))) => :n,
        endpoint => (x -> mean(collect(skipmissing(x)))) => :mean,
        endpoint => (x -> median(collect(skipmissing(x)))) => :median,
        endpoint => (x -> std(collect(skipmissing(x)))) => :sd,
        endpoint => (x -> cv_percent(collect(skipmissing(x)))) => :cv_pct,
    )
    sort!(out, :formulation)
    return out
end

# ----------------------------
# BE analysis wrapper
# ----------------------------
function run_be(df::DataFrame; endpoint::Symbol)
    # Standard ABE (TOST / 90% CI within 80-125)
    return pumas_be(df, StandardBioequivalenceCriterion; endpoint=endpoint)
end

# ----------------------------
# Simple planning-power approximation (simulation)
# ----------------------------
# “Power” here = P(90% CI fully within 80–125) under assumed CV and true GMR.
# This is not a replacement for a full validated planning tool; it’s a didactic approximation.
try
    using Distributions
catch
    @warn "Distributions not found. Installing..."
    Pkg.add("Distributions")
    using Distributions
end

logvar_from_cv(cv_pct) = log((cv_pct/100)^2 + 1)

function power_ci_sim(; n::Int, cv_pct::Real, gmr::Real=1.0, design::Symbol=:crossover,
                      alpha::Real=0.05, nsim::Int=12000, rng::AbstractRNG=Random.default_rng())
    μ  = log(gmr)
    s2 = logvar_from_cv(cv_pct)

    # Planning SE approximation:
    # - crossover: within-subject comparison ~ sqrt(2*s2/n)
    # - parallel:  between-subject ~ sqrt(4*s2/n)
    se = design == :crossover ? sqrt(2*s2/n) : sqrt(4*s2/n)

    dfree = max(n - 2, 1)
    tcrit = quantile(Distributions.TDist(dfree), 1 - alpha)
    lo = log(0.80); hi = log(1.25)

    ok = 0
    dist = Distributions.Normal(μ, se)
    for _ in 1:nsim
        est = rand(rng, dist)
        ci_lo = est - tcrit*se
        ci_hi = est + tcrit*se
        ok += (ci_lo >= lo && ci_hi <= hi) ? 1 : 0
    end
    return ok/nsim
end

function n_for_power(; target_power::Real, cv_pct::Real, design::Symbol=:crossover,
                     gmr::Real=1.0, n_min::Int=8, n_max::Int=220, step::Int=2, nsim::Int=8000)
    n0 = iseven(n_min) ? n_min : n_min + 1
    for n in n0:step:n_max
        p = power_ci_sim(n=n, cv_pct=cv_pct, gmr=gmr, design=design, nsim=nsim)
        if p >= target_power
            return (n=n, power=p)
        end
    end
    return (n=missing, power=missing)
end

# ----------------------------
# Plots
# ----------------------------
using StatsPlots  # ensure loaded
using Plots

function plot_endpoint_distributions(df::DataFrame, endpoint::Symbol, tag::String; outdir::String="figures")
    mkpath(outdir)

    # guard
    if !(String(endpoint) in names(df))
        return "Endpoint $(endpoint) not found."
    end
    if !("formulation" in names(df))
        error("Column :formulation missing. Columns: $(names(df))")
    end

    # explicit vectors (avoids Symbol-as-data errors)
    x = df[!, endpoint]
    g = df[!, :formulation]

    # histogram (by formulation)
    p1 = histogram(
        x;
        group=g,
        bins=:auto,
        xlabel=String(endpoint),
        ylabel="Frequency",
        title="$(tag): $(endpoint) distribution by formulation",
        legend=:topright
    )
    savefig(p1, joinpath(outdir, "$(tag)_$(endpoint)_hist.png"))

    # boxplot (by formulation)
    p2 = boxplot(
        g, x;
        xlabel="Formulation",
        ylabel=String(endpoint),
        title="$(tag): $(endpoint) boxplot (T vs R)"
    )
    savefig(p2, joinpath(outdir, "$(tag)_$(endpoint)_boxplot.png"))

    return (
        hist = joinpath(outdir, "$(tag)_$(endpoint)_hist.png"),
        box  = joinpath(outdir, "$(tag)_$(endpoint)_boxplot.png"),
    )
end


# Paired plot only for crossover (id has both R and T)
function plot_paired_crossover(df::DataFrame, endpoint::Symbol, tag::String)
    if !(has_col(df, :id) && has_col(df, endpoint) && has_col(df, :formulation))
        return "Missing columns for paired plot."
    end
    g = groupby(df, [:id, :formulation])
    subj = combine(g, endpoint => (x -> mean(collect(skipmissing(x)))) => :val)
    wide = unstack(subj, :formulation, :val)
    if !("R" in names(wide) && "T" in names(wide))
        return "Could not form R/T pairs for paired plot."
    end
    rename!(wide, "R" => :R, "T" => :T)
    keep = wide[.!ismissing.(wide.R) .&& .!ismissing.(wide.T), :]
    rr = Float64.(keep.R)
    tt = Float64.(keep.T)

    p = scatter(rr, tt,
        xlabel="Reference", ylabel="Test",
        title="$(tag): Paired subjects (mean per subject) — $(endpoint)",
        legend=false)
    mn = min(minimum(rr), minimum(tt))
    mx = max(maximum(rr), maximum(tt))
    plot!(p, [mn, mx], [mn, mx], lw=2)
    outp = joinpath(FIGD, "$(tag)_paired_$(endpoint).png")
    savefig(p, outp)
    return outp
end

# ----------------------------
# LOAD DATASETS
# ----------------------------
parallel_name  = find_dataset("FSL2015_5")
crossover_name = find_dataset("SLF2014_5")

parallel_df  = parallel_name  === nothing ? nothing : DataFrame(PharmaDatasets.dataset(parallel_name))
crossover_df = crossover_name === nothing ? nothing : DataFrame(PharmaDatasets.dataset(crossover_name))

if parallel_df !== nothing;  ensure_formulation!(parallel_df);  end
if crossover_df !== nothing; ensure_formulation!(crossover_df); end

# ----------------------------
# MAIN OUTPUT FILE
# ----------------------------
open(OUTTXT, "w") do io
    println(io, "Student Project 1 — Parallel vs Crossover Design Comparison")
    println(io, "Generated: ", Dates.now(), "\n")

    save_block(io, "Detected dataset names", (parallel=parallel_name, crossover=crossover_name))

    if parallel_df !== nothing
        save_block(io, "Parallel dataset columns", names(parallel_df))
        save_block(io, "Parallel N rows / N subjects", (rows=nrow(parallel_df), subjects=has_col(parallel_df,:id) ? length(unique(parallel_df.id)) : missing))
    else
        save_block(io, "Parallel dataset missing", "Could not locate FSL2015_5 in PharmaDatasets.datasets().")
    end

    if crossover_df !== nothing
        save_block(io, "Crossover dataset columns", names(crossover_df))
        save_block(io, "Crossover N rows / N subjects", (rows=nrow(crossover_df), subjects=has_col(crossover_df,:id) ? length(unique(crossover_df.id)) : missing))
    else
        save_block(io, "Crossover dataset missing", "Could not locate SLF2014_5 in PharmaDatasets.datasets().")
    end

    # ---- Task 1: descriptives
    if parallel_df !== nothing
        for ep in (:AUC, :Cmax)
            ok, res = try_run(io, "Task 1 — Parallel descriptives: $(ep)", () -> summarize_endpoint(parallel_df, ep))
        end
    end
    if crossover_df !== nothing
        for ep in (:AUC, :Cmax)
            ok, res = try_run(io, "Task 1 — Crossover descriptives: $(ep)", () -> summarize_endpoint(crossover_df, ep))
        end
    end

    # ---- Task 2: BE analysis
    if parallel_df !== nothing
        for ep in (:AUC, :Cmax)
            if has_col(parallel_df, ep)
                try_run(io, "Task 2 — Parallel BE: $(ep)", () -> run_be(parallel_df; endpoint=ep))
            end
        end
    end
    if crossover_df !== nothing
        for ep in (:AUC, :Cmax)
            if has_col(crossover_df, ep)
                try_run(io, "Task 2 — Crossover BE: $(ep)", () -> run_be(crossover_df; endpoint=ep))
            end
        end
    end

    # ---- Task 3: planning power + n for 80% (using CV from data, if possible)
    function extract_cvr_like(summary_df::DataFrame, formulation_ref::String="R")
        # use CV% of reference formulation row if available
        if !("formulation" in names(summary_df)) || !("cv_pct" in names(summary_df))
            return nothing
        end
        r = summary_df[summary_df.formulation .== formulation_ref, :]
        return nrow(r) == 1 ? r.cv_pct[1] : nothing
    end

    # Crossover
    if crossover_df !== nothing && has_col(crossover_df, :AUC)
        s = summarize_endpoint(crossover_df, :AUC)
        cv = s isa DataFrame ? extract_cvr_like(s) : nothing
        if cv !== nothing && cv !== missing
            p  = power_ci_sim(n=has_col(crossover_df,:id) ? length(unique(crossover_df.id)) : 24,
                              cv_pct=cv, design=:crossover, nsim=8000)
            n80 = n_for_power(target_power=0.80, cv_pct=cv, design=:crossover)
            save_block(io, "Task 3 — Crossover planning power (AUC, from pilot CV)", (cv_pct=cv, approx_power=p, n_for_80pct=n80))
        end
    end

    # Parallel
    if parallel_df !== nothing && has_col(parallel_df, :AUC)
        s = summarize_endpoint(parallel_df, :AUC)
        cv = s isa DataFrame ? extract_cvr_like(s) : nothing
        if cv !== nothing && cv !== missing
            p  = power_ci_sim(n=has_col(parallel_df,:id) ? length(unique(parallel_df.id)) : 48,
                              cv_pct=cv, design=:parallel, nsim=8000)
            n80 = n_for_power(target_power=0.80, cv_pct=cv, design=:parallel)
            save_block(io, "Task 3 — Parallel planning power (AUC, from pilot CV)", (cv_pct=cv, approx_power=p, n_for_80pct=n80))
        end
    end

    println(io, "\nDONE.\n")
end

# ----------------------------
# FIGURES
# ----------------------------
figlog = Dict{String,Any}()

if parallel_df !== nothing
    for ep in (:AUC, :Cmax)
        figlog["parallel_$ep"] = plot_endpoint_distributions(parallel_df, ep, "parallel")
    end
end

if crossover_df !== nothing
    for ep in (:AUC, :Cmax)
        figlog["crossover_$ep"] = plot_endpoint_distributions(crossover_df, ep, "crossover")
        figlog["crossover_paired_$ep"] = plot_paired_crossover(crossover_df, ep, "crossover")
    end
end

# ----------------------------
# DOCS: Questions + Report
# ----------------------------
open(QMD, "w") do io
    println(io, "# Student Project 1 — Questions & What This Repo Answers\n")
    println(io, "This file summarizes the **project questions** and points to where results appear in the repo.\n")

    println(io, "## Task 1 — Data exploration\n")
    println(io, "- **How many subjects are in each study?**  \n  See `outputs/student_project_1_outputs.txt` (\"N rows / N subjects\").")
    println(io, "- **What are the key variables?**  \n  See `outputs/student_project_1_outputs.txt` (\"dataset columns\").")
    println(io, "- **How do designs differ in organization?**  \n  Parallel has one treatment per subject; crossover has multiple periods and uses within-subject comparisons.\n")

    println(io, "## Task 2 — Bioequivalence analysis\n")
    println(io, "- **Compare GMR and 90% CI between designs** (AUC/Cmax when available).  \n  See `outputs/student_project_1_outputs.txt` (\"Parallel BE\" and \"Crossover BE\").")
    println(io, "- **Which design gives tighter CI and why?**  \n  Crossover typically yields tighter CIs because it controls between-subject variability.\n")

    println(io, "## Task 3 — Power and sample size\n")
    println(io, "- **What sample size is needed for ~80% power?**  \n  See the \"planning power\" blocks in `outputs/student_project_1_outputs.txt` (didactic simulation-based approximation).")
    println(io, "- **Why crossover is more efficient**: within-subject comparisons reduce unexplained variance → fewer subjects needed.\n")

    println(io, "## Task 4 — Practical considerations\n")
    println(io, "- Parallel design avoids carryover/washout but usually needs more participants.")
    println(io, "- Crossover reduces sample size but requires washout and more complex conduct/analysis.\n")

    println(io, "## Visualizations\n")
    println(io, "- Distribution plots (hist/box) and paired plots (crossover) are stored in `figures_project1/`.\n")
end

open(OUTMD, "w") do io
    println(io, "# Student Project 1 — Parallel vs Crossover Design Comparison\n")
    println(io, "- Generated: **$(Dates.now())**\n")

    println(io, "## Repository outputs\n")
    println(io, "- **Outputs (console-like):** `outputs/student_project_1_outputs.txt`")
    println(io, "- **Questions mapping:** `docs/student_project_1_questions.md`")
    println(io, "- **Figures:** `figures_project1/`\n")

    println(io, "## What we did\n")
    println(io, "1. Loaded a **parallel** dataset (`FSL2015_5`) and a **crossover** dataset (`SLF2014_5`) from PharmaDatasets.")
    println(io, "2. Computed **descriptive statistics** (mean, median, SD, CV%) by formulation for AUC/Cmax when available.")
    println(io, "3. Ran **Standard ABE** via `pumas_be` for endpoints present.")
    println(io, "4. Produced **visual comparisons** (histograms, boxplots; and paired plots for crossover).")
    println(io, "5. Included a **didactic planning-power approximation** to contrast design efficiency.\n")

    println(io, "## Key interpretation (student-level, regulator-aware)\n")
    println(io, "- **Crossover designs** are generally more statistically efficient because each subject serves as their own control.")
    println(io, "- **Parallel designs** are simpler to run (no washout / carryover concerns) but usually require larger sample sizes.")
    println(io, "- For BE, the operational target is often **≥80% power**, meaning a true-BE situation will pass about 8 out of 10 times under planning assumptions.\n")

    println(io, "## Figures produced\n")
    println(io, "See `figures_project1/` for endpoint distributions and paired crossover comparisons.\n")

    println(io, "## How to run\n")
    println(io, "```bash")
    println(io, "julia --project=. scripts/student_project_1_run_and_save.jl")
    println(io, "```\n")
end

println("✅ Project 1 done.")
println(" - ", OUTTXT)
println(" - ", OUTMD)
println(" - ", QMD)
println(" - figures in: ", FIGD)
