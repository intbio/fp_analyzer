# FP Analyzer

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/intbio/fp_analyzer/HEAD?urlpath=%2Fdoc%2Ftree%2Fnotebooks%2Ffp_analyzer.ipynb)

An interactive Jupyter notebook for analyzing **fluorescence polarization (FP)
saturation binding** experiments. Reads raw plate-reader output, applies
instrument corrections (G-factor, blank subtraction, Q-factor), fits direct
binding models, computes a Z′-factor for assay quality, and exports
publication-ready tables and figures.

## Quick start

**Option A — try it online (no install).** Click the Binder badge above. The
first launch takes 1–3 min while the environment is built (subsequent
launches are cached and start in seconds). The example data files are bundled with the notebook and load automatically, so you can click straight through. Additional example datasets are also available in the examples/ folder of the Binder session.

**Option B — run locally.**

```bash
git clone https://github.com/intbio/fp_analyzer
cd REPO
pip install -r requirements.txt
jupyter notebook notebooks/fp_analyzer.ipynb
```

In either case, run the cells from top to bottom. The notebook is laid out
as **six linear steps**: Setup → Load data → QC → Fit → Summary → Export.
Each step has a short instruction block followed by a code cell that
displays a widget panel.

## What it does

* Parses a layout file (CSV/XLSX) and a data file (custom CSV or native BMG
  Labtech FP export).
* Computes the instrument **G-factor** from dedicated calibration wells with
  free fluorophore of known reference anisotropy.
* Subtracts buffer **blanks** and computes per-well anisotropy *r* from
  raw parallel and perpendicular intensities.
* Estimates the **Q-factor** (bound/free total fluorescence ratio) per
  condition from the titration extremes.
* Aggregates titration replicates into a saturation curve.
* Fits one of three **direct-binding models**:
  * `fp_quadratic` (default) — exact closed-form solution accounting for
    ligand depletion (Roehrl et al. 2004; Jarmoskaite et al. 2020).
  * `one_site` — simple hyperbolic, valid only when [L\*] ≪ K<sub>D</sub>.
  * `hill` — phenomenological Hill model with coefficient n<sub>H</sub>.
* Reports K<sub>D</sub> with confidence intervals and reduced χ².
* Computes **Z′-factor** (Zhang 1999) per condition for assay quality.
* Generates a **Scatchard plot** as a visual diagnostic.
* Bundles results into a downloadable ZIP (CSV tables + PNG figures).

## What it does not do

* Does **not** fit competitive displacement experiments (only direct binding).
* Does **not** model non-specific binding explicitly — visual diagnosis via
  the Scatchard plot only.
* Does **not** verify that the binding reaction has reached equilibrium —
  this must be established by the experimenter prior to data collection
  (see Jarmoskaite et al. 2020 for guidance).

## Repository layout

```
fp-analyzer/
├── README.md                       # this file
├── USER_GUIDE.md                   # full user guide with theory
├── requirements.txt                # Python dependencies (used by Binder)
├── notebooks/
│   └── fp_analyzer.ipynb           # the interactive notebook
├── examples/
│   ├── example_layout.csv          # synthetic 384-well layout
│   ├── example_data.csv            # raw intensities, custom CSV format
│   ├── example_data_bmg.csv        # same data in BMG Labtech export format
│   ├── layout_dCas9.csv			# 384-well plate layout, dCas9 specific DNA binding
│   ├── data_dCas9.csv				# raw BMG Labtech export, dCas9 specific DNA binding
│   ├── layout_ns_dCas9.csv			# 384-well plate layout, dCas9 non-specific DNA binding
│   └── data_ns_dCas9.csv			# raw BMG Labtech export, dCas9 non-specific DNA binding
└── docs/
    └── methods.docx                # ready-to-paste Methods section for a paper
```

## Data format requirements

See **[USER_GUIDE.md §3](USER_GUIDE.md#3-input-files)** for the complete
specification of the layout and data file formats, including a worked example.

In short:

* **Layout file** — CSV/XLSX with columns `well, condition, sample_type,
  concentration, replicate`. Each well is one of four sample types:
  `titration`, `substrate`, `blank`, `fluorophore`.
* **Data file** — either a custom CSV with columns `well, parallel,
  perpendicular`, or a native BMG Labtech FP export (autodetected from
  `Raw Data (parallel)` / `Raw Data (perpendicular)` block headers).

Try the synthetic example in `examples/` first to see the expected
structure.

## Theory

A short theoretical primer covering FP basics (Perrin equation, G-factor,
anisotropy from raw intensities, Q-factor mixing, choice of binding model,
Z′-factor) is in **[USER_GUIDE.md §2](USER_GUIDE.md#2-theoretical-background)**.
A condensed Methods section ready for inclusion in a paper is in
`docs/methods.docx`.

Key references:

* Lakowicz JR (2006). *Principles of Fluorescence Spectroscopy*, 3rd ed.
  Springer.
* Tranter GE (2010). Fluorescence polarization and anisotropy. In
  *Encyclopedia of Spectroscopy and Spectrometry*, 2nd ed., 625–627.
* Jameson DM, Croney JC (2003). *Comb. Chem. High Throughput Screen.*
  **6**, 167.
* Moerke NJ (2009). *Curr. Protoc. Chem. Biol.* **1**, 1.
* Roehrl MHA, Wang JY, Wagner G (2004). *Biochemistry* **43**, 16056.
* Jarmoskaite I, AlSadhan I, Vaidyanathan PP, Herschlag D (2020). *eLife*
  **9**, e57264.
* Zhang JH, Chung TDY, Oldenburg KR (1999). *J. Biomol. Screen.* **4**, 67.
* Newville M, Stensitzki T, Allen DB, Ingargiola A (2014). LMFIT.
  doi:10.5281/zenodo.11813.

## License

MIT License

Copyright (c) 2026 Integrative Biology Group

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
