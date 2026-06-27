# FP Analyzer — User Guide

A Jupyter-based tool for analyzing fluorescence polarization (FP) saturation
binding experiments. Reads raw plate-reader output, applies instrument
corrections, fits direct-binding models, and exports publication-ready
results.

---

## 1. Overview

`fp_analyzer` is built around a single class, `FPAnalyzer`, which is
exposed through a tabbed `ipywidgets` interface (Load → QC → Fit →
Summary → Export). It expects:

- **a layout file** (CSV or XLSX) that maps each well of the plate to a
  sample type and metadata;
- **a data file** containing parallel and perpendicular intensities — either
  a custom CSV (`well, parallel, perpendicular`) or a native BMG Labtech
  FP export with per-cycle "Raw Data (parallel)" / "Raw Data
  (perpendicular)" blocks.

What the tool does, in order:

1. parses and merges layout + data;
2. computes the instrument *G*-factor from calibration wells;
3. subtracts buffer blanks and computes anisotropy *r* per well;
4. estimates the *Q*-factor (bound/free intensity ratio) per condition;
5. aggregates titration replicates into a saturation curve;
6. fits a binding model and reports K<sub>D</sub> with confidence
   intervals;
7. produces a Scatchard plot and a *Z*′-factor as quality diagnostics;
8. exports a ZIP bundle (CSV tables, fit reports, PNG figures).

What it does **not** do (see also §8 *Limitations*): it does not analyze
competitive displacement experiments, does not model non-specific
binding explicitly, and does not verify that the binding reaction has
reached equilibrium — the latter must be established by the experimenter
before data collection.

---

## 2. Theoretical background

This section is a brief operational summary of the physics needed to use
the tool sensibly. For full derivations see Lakowicz (2006) and Tranter
(2010); for the analytical framework adapted here, Roehrl et al. (2004)
and Jarmoskaite et al. (2020).

### 2.1 Why fluorescence polarization works for binding

When a fluorophore is excited with plane-polarized light, it preferentially
absorbs photons whose electric vector is aligned with the absorption
transition dipole, and emits with a polarization that depends on how much
the molecule has rotated during the excited-state lifetime τ<sub>F</sub>.
The relevant timescale is the rotational correlation time τ<sub>C</sub>,
which for a roughly spherical molecule of mass *M* in water is

```
τ_C ≈ 0.3 × M    (ns, with M in kDa)
```

A free, low-MW labelled probe (τ<sub>C</sub> ≪ τ<sub>F</sub>) tumbles fast
enough to depolarize the emission almost completely, while the same probe
bound to a much larger protein (τ<sub>C</sub> ≳ τ<sub>F</sub>) retains
much of its polarization. The change in anisotropy upon binding is the
signal exploited by the FP assay.

For typical fluorescein-class probes (τ<sub>F</sub> ≈ 4 ns), a useful
binding signal requires the protein to be at least ~10 kDa and the
labelling site to be reasonably rigid relative to the binding interface.
Floppy linkers between the dye and the binding epitope (the "propeller
effect") attenuate the bound-state anisotropy.

### 2.2 From raw intensities to anisotropy

Two emission intensities are recorded for each well: parallel
(*I*<sub>∥</sub>) and perpendicular (*I*<sub>⊥</sub>) to the excitation
polarization plane. Polarization *P* and anisotropy *r* are defined as

```
P = (I∥ − G·I⊥) / (I∥ + G·I⊥)
r = (I∥ − G·I⊥) / (I∥ + 2G·I⊥)
P = 3r / (2 + r)
```

The *G*-factor corrects for the fact that the parallel and perpendicular
detection channels have unequal sensitivities. **The tool uses anisotropy
internally**, not polarization — only *r* obeys a simple linear additivity
rule for mixtures of fluorescent species, which is essential for
extracting bound fractions from a two-state mixture.

**G-factor.** The pipeline derives *G* from dedicated calibration wells
containing the free fluorophore at the same concentration as the
labelled ligand. Given a user-supplied limiting anisotropy
*A*<sub>ref</sub> for the free probe (default 0.035), each calibration well contributes one *G* value:

```
G = I∥·(1 − A_ref) / [I⊥·(1 + 2·A_ref)]
```

The reported *G* is the mean across all calibration wells, with its
standard error. Typical values fall in the 0.7–1.3 range for healthy
instruments.

**Blank subtraction.** Buffer-only wells are averaged within each
experimental condition and subtracted from *I*<sub>∥</sub> and
*I*<sub>⊥</sub> before *r* is computed. This is important when the
labelled probe is dilute and the well-to-well buffer fluorescence is
non-negligible.

**Q-factor.** When the fluorescence quantum yield of the probe changes on
binding (quenching or enhancement), the contribution of the bound
species to the observed anisotropy is weighted accordingly. The
two-state mixing equation with a *Q*-factor (Roehrl et al., 2004) is

```
A_obs = [Q·F_SB·A_bound + (1 − F_SB)·A_free] / [1 − (1 − Q)·F_SB]
```

where *F*<sub>SB</sub> is the bound fraction, and
*Q* = *I*<sub>total,bound</sub> / *I*<sub>total,free</sub> is the ratio of
total fluorescence intensities of the bound vs. free probe. The pipeline
estimates *Q* empirically as the ratio of the mean
*I*<sub>∥</sub> + 2*G*·*I*<sub>⊥</sub> at the highest vs. lowest protein
concentration of the titration. *Q* = 1 reduces the equation to the
simple linear mixing form. Importantly, K<sub>D</sub> estimates are
robust to substantial errors in *Q* even when *A*<sub>bound</sub> is not
(Roehrl et al., 2004, Fig. 4B).

### 2.3 Binding regimes — which model to use

The single most common source of incorrect K<sub>D</sub> values in the
literature is fitting a hyperbolic equation to data taken in a regime
where it does not apply. Choose the model based on the relationship
between [L]<sub>T</sub> (the fixed labelled-ligand concentration) and the
expected K<sub>D</sub>:

| Regime                         | Recommended model | Notes |
|--------------------------------|-------------------|-------|
| [L]<sub>T</sub> ≪ K<sub>D</sub>           | One-site hyperbolic | "Binding regime"; hyperbolic and quadratic give the same answer. |
| [L]<sub>T</sub> ~ K<sub>D</sub>           | **Quadratic (default)** | Intermediate regime; quadratic is necessary, hyperbolic biases K<sub>D</sub> high. |
| [L]<sub>T</sub> ≫ K<sub>D</sub>           | None reliable     | "Titration regime"; only an upper bound on K<sub>D</sub> can be obtained. Lower [L]<sub>T</sub> if possible. |
| Sigmoidal, non-hyperbolic shape | Hill (diagnostic) | Indicates apparent cooperativity. Use the Hill *n*<sub>H</sub> as a phenomenological descriptor, not as evidence of a mechanism. |

**The three available models**, all expressed as the bound fraction
*F*<sub>SB</sub> of the labelled ligand:

```
One-site hyperbolic:   F_SB = [P]_T / (K_D + [P]_T)

Quadratic (default):   F_SB = { K_D + [L]_T + [P]_T
                              − √[(K_D + [L]_T + [P]_T)² − 4·[L]_T·[P]_T] }
                              / (2·[L]_T)

Hill:                  F_SB = [P]_T^nH / (K_D^nH + [P]_T^nH)
```

The fit is performed on the observed anisotropy via the *Q*-corrected
mixing equation in §2.2, using non-linear least squares (Levenberg–
Marquardt) as implemented in `lmfit`.

### 2.4 What good data look like

A clean saturation curve (i) starts on a low plateau at the lowest
protein concentrations close to *A*<sub>free</sub> measured from the
substrate-only wells, (ii) rises through a single sigmoidal transition,
(iii) reaches a clear high plateau (*A*<sub>bound</sub>) at the highest
protein concentrations, and (iv) is reproducible across at least three
technical replicates with small SEMs. The Scatchard plot
(*B*/[L]<sub>T</sub> vs. *B*) of well-behaved 1:1 binding is a straight
line; pronounced curvature suggests multiple sites, cooperativity, or
non-specific binding.

For screening applications, the *Z*′-factor (Zhang et al., 1999)

```
Z' = 1 − [3·(σ_pos + σ_neg)] / |μ_pos − μ_neg|
```

with positive control = saturated complex (highest [P]<sub>T</sub>
titration wells) and negative control = free probe (substrate-only
wells), summarizes assay robustness: *Z*′ > 0.5 is excellent,
0 < *Z*′ < 0.5 is marginal but usable, *Z*′ ≤ 0 is unreliable.

---

## 3. Input files

### 3.1 Layout file

CSV or XLSX with one row per well. Required columns:

| Column          | Type    | Notes |
|-----------------|---------|-------|
| `well`          | string  | e.g. `A1`, `H12`. Must be unique. |
| `condition`     | string  | Free-form label, e.g. `WT`, `mutant_K100A`. |
| `sample_type`   | string  | One of `titration`, `substrate`, `blank`, `fluorophore`. |
| `concentration` | numeric | Protein concentration in **nM**. Leave blank for non-titration wells. |
| `replicate`     | int     | Replicate index (1, 2, 3, …). Leave blank for non-titration wells. |

The four `sample_type` values map onto the well classes in §1:

- `titration` — labelled ligand + serial dilution of protein. Goes into
  the binding curve.
- `substrate` — labelled ligand only, no protein. Defines
  *A*<sub>free</sub> and the *Z*′ negative control.
- `blank` — assay buffer only. Subtracted from intensities, per
  condition.
- `fluorophore` — free fluorophore alone, used for *G*-factor
  calibration.

Example:

```csv
well,condition,sample_type,concentration,replicate
A1,WT,blank,,
A2,WT,fluorophore,,
A3,WT,substrate,,1
A4,WT,substrate,,2
B1,WT,titration,1000,1
B2,WT,titration,500,1
B3,WT,titration,250,1
...
```

### 3.2 Data file

**Custom CSV format** — three required columns:

```csv
well,parallel,perpendicular
A1,1234.5,1198.3
A2,9876.2,9543.1
...
```

**BMG Labtech native export** — the tool autodetects exports that contain
`Raw Data (parallel)` and `Raw Data (perpendicular)` plate blocks, and
parses all cycles.

---

## 4. Step-by-step workflow

The notebook has two code cells (run once) that define the class and the
UI. Everything else is done through the tabs.

### 4.1 Tab "Load"

Attach the layout and data files, set the labelled-ligand concentration
in nM (the constant component of the binding equilibrium), optionally
override *A*<sub>ref</sub> (default 0.035) or *G* (default: derived from
calibration wells), and click *Load*. Errors at this stage are usually
about layout/data mismatch (well naming, missing columns, unknown
`sample_type`); the message tells you exactly what to fix.

### 4.2 Tab "QC"

Inspect the per-condition QC report:

- *G*-factor and its standard error (one global value);
- *Q*-factor (per condition);
- *A*<sub>free</sub> with SE (per condition);
- mean blank intensities;
- number of titration points and replicates;
- *Z*′-factor (per condition);
- warnings (e.g. "*Q* deviates from 1 by > 0.3").

### 4.3 Tab "Fit"

For each condition, choose a model (default: `fp_quadratic`), optionally
override initial values or fix parameters, and run the fit. The fit
report (lmfit format) is shown along with the saturation plot overlaying
the data with SEM error bars and the best-fit curve. K<sub>D</sub> is
reported in nM with its standard error.

Per-parameter defaults:

- `Kd` — free, initial 50 nM;
- `A_bound` — free, initial 0.200;
- `A_free` — fixed to the value measured from substrate-only wells;
- `Q` — fixed to the empirical estimate from the titration extremes;
- `Dl` (= [L]<sub>T</sub>) — fixed to the value entered in the Load tab.

To explore sensitivity to assumptions, release `A_free` or `Q` from
their fixed values and refit; in well-behaved data K<sub>D</sub> should
not change much.

### 4.4 Tabs "Summary" and "Export"

Summary aggregates K<sub>D</sub> across conditions in a single table and
a comparative plot. Export bundles the results into a downloadable ZIP:
saturation data (CSV), Scatchard data (CSV), fit parameters (CSV), QC
report (TXT), and PNG figures of the saturation, Scatchard, and *G*
calibration plots.

---

## 5. Interpreting results

**A converged fit with small SE on K<sub>D</sub> and a low reduced χ²
is necessary but not sufficient.** Always cross-check:

- Does the fitted *A*<sub>bound</sub> match the high-concentration
  plateau visually? If the curve has not saturated, *A*<sub>bound</sub>
  and K<sub>D</sub> are correlated and both unreliable. Extend the
  titration to higher [P]<sub>T</sub>.
- Is the lowest [P]<sub>T</sub> close to *A*<sub>free</sub>? If not, the
  baseline is poorly defined and *A*<sub>free</sub> from substrate-only
  wells is more trustworthy than a free fit.
- Is K<sub>D</sub> > [L]<sub>T</sub> by at least a factor of ~3? If
  K<sub>D</sub> ≲ [L]<sub>T</sub> you are in or near the titration regime
  (§2.3) and the value should be treated as an upper limit. Lower
  [L]<sub>T</sub> in a follow-up experiment if the signal allows.
- Does the Scatchard plot look linear? Curvature is a flag for multi-site
  or non-specific binding (§6).

**Switching from hyperbolic to quadratic** is essentially free in terms
of code (just change the model selector) and is the right default. The
hyperbolic model is included only for backward compatibility and for the
[L]<sub>T</sub> ≪ K<sub>D</sub> regime.

**Hill *n*<sub>H</sub>** values significantly different from 1 (typically
|*n*<sub>H</sub> − 1| > 0.2 with non-overlapping error bars) suggest
apparent cooperativity but do **not** prove a specific mechanism;
treat as a hypothesis-generating observation.

---

## 6. Troubleshooting

**"G deviates significantly from 1."** Real *G* values are normally in
0.7–1.3. Larger deviations point to (i) a misconfigured polarizer or
emission filter, (ii) mismatch between the free-fluorophore concentration
in calibration wells and what the *A*<sub>ref</sub> default assumes, or
(iii) photobleaching during the read. Re-measure calibration wells fresh,
or override *G* manually after determining it independently.

**"Q deviates from 1 by more than 0.3."** A *Q* far from 1 means the
probe's quantum yield changes substantially on binding. This is not in
itself a problem — the fit accounts for it — but it can indicate that the
fluorophore is at or near the binding interface, in which case the
labelling site may be perturbing the interaction. Consider an
alternative labelling position. K<sub>D</sub> is robust to *Q*
mis-estimation; *A*<sub>bound</sub> is not.

**Fit does not converge.** Most common causes: (i) the initial
K<sub>D</sub> guess is far from reality (try a value near the visible
midpoint of the curve); (ii) the curve has not saturated, leaving
K<sub>D</sub> and *A*<sub>bound</sub> degenerate (extend the titration);
(iii) data are too noisy for the chosen model (try the simpler hyperbolic
form and inspect residuals).

**Scatchard is non-linear.** Possible causes, in rough order of
likelihood: (i) titration regime — repeat with lower [L]<sub>T</sub>; (ii)
non-specific binding to the wells, tubes, or another component (add 0.01
% Tween-20 or 0.1 mg/mL BSA); (iii) genuine multi-site or cooperative
binding — fit Hill as a diagnostic; (iv) protein aggregation at high
concentrations.

**Replicates disagree by more than the SEM suggests.** Check for
pipetting errors, edge effects on the plate (corner wells often give
slightly different reads), and protein activity loss during the assay
(re-run with fresh protein and shorter total assay time).

---

## 7. References

- Jameson, D. M. & Croney, J. C. (2003). Fluorescence polarization: past,
  present and future. *Comb. Chem. High Throughput Screen.* **6**,
  167–173.
- Jarmoskaite, I., AlSadhan, I., Vaidyanathan, P. P. & Herschlag, D.
  (2020). How to measure and evaluate binding affinities. *eLife*
  **9**, e57264.
- Lakowicz, J. R. (2006). *Principles of Fluorescence Spectroscopy*,
  3rd ed. Springer.
- Moerke, N. J. (2009). Fluorescence polarization (FP) assays for
  monitoring peptide–protein or nucleic acid–protein binding. *Curr.
  Protoc. Chem. Biol.* **1**, 1–15.
- Newville, M., Stensitzki, T., Allen, D. B. & Ingargiola, A. (2014).
  *LMFIT: Non-linear least-square minimization and curve-fitting for
  Python.* Zenodo. doi:10.5281/zenodo.11813.
- Roehrl, M. H. A., Wang, J. Y. & Wagner, G. (2004). A general framework
  for development and data analysis of competitive high-throughput
  screens for small-molecule inhibitors of protein–protein interactions
  by fluorescence polarization. *Biochemistry* **43**, 16056–16066.
- Tranter, G. E. (2010). Fluorescence polarization and anisotropy. In
  *Encyclopedia of Spectroscopy and Spectrometry*, 2nd ed., 625–627.
  Elsevier.
- Zhang, J. H., Chung, T. D. Y. & Oldenburg, K. R. (1999). A simple
  statistical parameter for use in evaluation and validation of high
  throughput screening assays. *J. Biomol. Screen.* **4**, 67–73.
