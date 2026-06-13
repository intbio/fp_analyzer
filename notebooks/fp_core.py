"""
fp_core.py — Fluorescence Polarization analysis core + interactive UI.

This module contains everything: the FPAnalyzer class, the binding models,
and the ipywidgets panels. The notebook only needs to import and call the
`show_*` functions; no analysis code is visible in the notebook itself.

Public API
----------
build_ui()                       -> displays the full tabbed interface
show_load_panel()                -> Step 2 (load data)
show_qc_panel()                  -> Step 3 (quality control)
show_fit_panel()                 -> Step 4 (fit binding model)
show_summary_panel()             -> Step 5 (summary)
show_export_panel()              -> Step 6 (export ZIP)
get_analyzer()                   -> returns the current FPAnalyzer instance
"""
from __future__ import annotations
import warnings, re, io, base64, zipfile, tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional, Union

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lmfit import Model, Parameters

import ipywidgets as widgets
from IPython.display import display, clear_output


# ══════════════════════════════════════════════════════════════════
# Exceptions
# ══════════════════════════════════════════════════════════════════
class FPError(Exception):   pass
class LayoutError(FPError): pass
class DataError(FPError):   pass
class MergeError(FPError):  pass
class FitError(FPError):    pass


# ══════════════════════════════════════════════════════════════════
# QC report
# ══════════════════════════════════════════════════════════════════
@dataclass
class QCReport:
    condition: str
    G: float; G_se: float; Q: float
    A_free: float; A_free_se: float
    blank_prl: float; blank_prp: float
    n_titration: int; n_replicates: int
    warnings: list = field(default_factory=list)

    def __str__(self):
        lines = [
            f"Condition      : {self.condition}",
            f"  G-factor     : {self.G:.4f} +/- {self.G_se:.4f}",
            f"  Q-factor     : {self.Q:.4f}",
            f"  A_free       : {self.A_free:.4f} +/- {self.A_free_se:.4f}",
            f"  Blank prl    : {self.blank_prl:.1f}",
            f"  Blank prp    : {self.blank_prp:.1f}",
            f"  Titration pts: {self.n_titration}",
            f"  Replicates   : {self.n_replicates}",
        ]
        if self.warnings:
            lines += [f"  !! {w}" for w in self.warnings]
        return "\n".join(lines)


# ══════════════════════════════════════════════════════════════════
# Binding models
# ══════════════════════════════════════════════════════════════════
def _model_fp_quadratic(x, Kd, A_bound, Q, Dl, A_free):
    """
    Quadratic binding equation for direct FP titration.
    x   = [protein]_total (variable)
    Dl  = [labelled ligand]_total (fixed)
    Returns observed anisotropy via the Q-corrected mixing equation
    (Roehrl et al. 2004, eq 6 and eq 39; Jarmoskaite et al. 2020, eq 5).
    """
    disc = np.maximum((Kd + Dl + x)**2 - 4*Dl*x, 0.0)
    F    = (Kd + Dl + x - np.sqrt(disc)) / (2*Dl)
    F    = np.clip(F, 0.0, 1.0)
    return (F*(Q*A_bound - A_free) + A_free) / (1 - (1-Q)*F)

def _model_one_site(x, Kd, A_free, A_bound):
    """Simple hyperbolic one-site binding (valid when [L*] << Kd)."""
    return A_free + (A_bound - A_free) * x / (Kd + x)

def _model_hill(x, Kd, A_free, A_bound, n):
    """Hill model with cooperativity coefficient n."""
    x = np.asarray(x, dtype=np.float64)
    # Guard against 0**n and overflow: work in a safe positive domain
    xs = np.clip(x, 1e-12, None)
    Kd = max(Kd, 1e-12)
    xn  = np.power(xs, n)
    Kdn = np.power(Kd, n)
    return A_free + (A_bound - A_free) * xn / (Kdn + xn)

_BUILTIN_MODELS = {
    'fp_quadratic': _model_fp_quadratic,
    'one_site':     _model_one_site,
    'hill':         _model_hill,
}

# Which parameters each model needs, and which are normally fixed from QC
_MODEL_PARAMS = {
    'fp_quadratic': {'free': ['Kd', 'A_bound'],
                     'fixed_from_qc': ['Q', 'Dl', 'A_free']},
    'one_site':     {'free': ['Kd', 'A_bound'],
                     'fixed_from_qc': ['A_free']},
    'hill':         {'free': ['Kd', 'A_bound', 'n'],
                     'fixed_from_qc': ['A_free']},
}


# ══════════════════════════════════════════════════════════════════
# FPAnalyzer
# ══════════════════════════════════════════════════════════════════
class FPAnalyzer:
    def __init__(self, layout_file, data_file, fixed_concentration,
                 A_ref=0.035, G=None, G_threshold=0.3, Q_threshold=0.3, cycle='last'):
        self.layout_path         = Path(layout_file)
        self.data_path           = Path(data_file)
        self.fixed_concentration = fixed_concentration
        self.A_ref               = A_ref
        self._G_override         = G
        self.G_threshold         = G_threshold
        self.Q_threshold         = Q_threshold
        self.cycle               = cycle
        self.G = self.G_se = None
        self.df_merged = self.df_saturation = self.df_scatchard = None
        self.qc_reports: dict = {}
        self._fit_results: dict = {}
        self._fit_curves:  dict = {}
        self._run_pipeline()

    def _run_pipeline(self):
        df_layout = self._parse_layout()
        df_data   = self._parse_data()
        self.df_merged = self._merge(df_layout, df_data)
        self._compute_G()
        self._compute_anisotropy()
        self._compute_qc()
        self._prepare_fitting_data()

    # ── Layout ─────────────────────────────────────────────────────
    def _parse_layout(self):
        if not self.layout_path.exists():
            raise LayoutError(f"Layout not found: {self.layout_path}")
        try:
            sfx = self.layout_path.suffix.lower()
            raw = pd.read_excel(self.layout_path) if sfx in ('.xlsx','.xls') \
                  else pd.read_csv(self.layout_path)
        except Exception as e:
            raise LayoutError(f"Cannot read layout: {e}") from e
        req  = {'well','condition','sample_type','concentration','replicate'}
        miss = req - set(raw.columns)
        if miss: raise LayoutError(f"Layout missing columns: {miss}")
        valid = {'titration','substrate','blank','fluorophore'}
        bad   = set(raw['sample_type'].dropna().unique()) - valid
        if bad: raise LayoutError(f"Unknown sample_type: {bad}")
        dups  = raw[raw.duplicated('well',keep=False)]['well'].unique()
        if len(dups): raise LayoutError(f"Duplicate wells: {list(dups)}")
        raw['concentration'] = pd.to_numeric(raw['concentration'], errors='coerce')
        if raw["replicate"].isna().all(): raw["replicate"] = 1
        else: raw["replicate"] = raw["replicate"].fillna(1).astype(int)
        return raw

    # ── Data ───────────────────────────────────────────────────────
    def _parse_data(self):
        if not self.data_path.exists():
            raise DataError(f"Data not found: {self.data_path}")
        fmt = self._detect_format()
        if fmt == 'custom':
            try:
                raw = pd.read_csv(self.data_path)
            except Exception as e:
                raise DataError(f"Cannot read data: {e}") from e
            df = self._validate_custom(raw)
        elif fmt == 'A':
            df = self._parse_format_a()
            df = self._select_cycle(df)
        else:
            raise DataError("Cannot detect format. Provide custom CSV "
                            "(well/parallel/perpendicular) or BMG Labtech FP export.")
        for col in ('parallel','perpendicular'):
            if df[col].isna().all():
                raise DataError(f"Column '{col}' is entirely NaN.")
        return df[['well','parallel','perpendicular']].copy()

    def _detect_format(self):
        try:
            with open(self.data_path) as f:
                head = [f.readline() for _ in range(15)]
        except Exception as e:
            raise DataError(f"Cannot read data: {e}") from e
        if any('Raw Data (parallel)' in l or 'Raw Data (perpendicular)' in l for l in head):
            return 'A'
        try:
            s = pd.read_csv(self.data_path, nrows=3)
            if {'well','parallel','perpendicular'}.issubset(s.columns):
                return 'custom'
        except Exception:
            pass
        return ''

    def _validate_custom(self, raw):
        miss = {'well','parallel','perpendicular'} - set(raw.columns)
        if miss: raise DataError(f"Custom data missing columns: {miss}")
        return raw

    _ROW_LABELS = list('ABCDEFGHIJKLMNOP')

    def _parse_plate_block(self, lines, start):
        hdr   = lines[start].rstrip('\n').split(',')
        while hdr and hdr[-1] == "": hdr.pop()
        ncols = len(hdr) - 1
        rows  = {}
        for j in range(1, 17):
            idx = start + j
            if idx >= len(lines): break
            cells = lines[idx].rstrip('\n').split(',')
            if cells[0] not in self._ROW_LABELS: break
            rows[cells[0]] = cells[1:ncols+1]
        if not rows: raise DataError(f"No plate data at line {start+1}")
        return pd.DataFrame(rows, index=hdr[1:]).T

    def _parse_format_a(self):
        with open(self.data_path) as f:
            lines = f.readlines()
        pairs, current = [], None
        def melt_block(df_b, val):
            return (df_b.reset_index().rename(columns={'index':'row'})
                    .melt(id_vars='row', var_name='column', value_name=val)
                    .assign(well=lambda d: d['row']+d['column']))
        for i, line in enumerate(lines):
            s = line.strip()
            m = re.match(r'^Cycle\s+(\d+)\s+\((.+?)\)', s)
            if m:
                current = {'cycle':int(m.group(1)),'time':m.group(2).strip(),
                           'had_header':True,'par':None,'perp':None}
                pairs.append(current); continue
            if re.match(r'^1\.\s*Raw Data \(parallel\)', s):
                if current is None or current['par'] is not None:
                    current = {'cycle':None,'time':None,'had_header':False,'par':None,'perp':None}
                    pairs.append(current)
                df_b = self._parse_plate_block(lines, i+2)
                df_l = melt_block(df_b,'parallel')
                df_l['parallel'] = pd.to_numeric(df_l['parallel'], errors='coerce')
                current['par'] = df_l
            if re.match(r'^2\.\s*Raw Data \(perpendicular\)', s):
                if current is None: continue
                df_b = self._parse_plate_block(lines, i+2)
                df_l = melt_block(df_b,'perpendicular')
                df_l['perpendicular'] = pd.to_numeric(df_l['perpendicular'], errors='coerce')
                current['perp'] = df_l
        if not pairs: raise DataError('Format A: no data blocks found.')
        used = sorted({p['cycle'] for p in pairs if p['had_header'] and p['cycle']})
        miss = sorted(set(range(1, max(used,default=0)+2)) - set(used))
        mi = 0
        for p in pairs:
            if not p['had_header'] or p['cycle'] is None:
                p['cycle'] = miss[mi] if mi < len(miss) else max(used)+mi+1
                mi += 1
        rows = []
        for p in pairs:
            if p['par'] is None or p['perp'] is None:
                warnings.warn(f"Incomplete pair cycle {p['cycle']}", UserWarning); continue
            mg = p['par'].merge(p['perp'][['well','perpendicular']], on='well', how='inner')
            mg['cycle'] = p['cycle']; mg['cycle_time'] = p['time']
            rows.append(mg)
        df = pd.concat(rows, ignore_index=True)
        df = df.dropna(subset=['parallel','perpendicular'], how='all').copy()
        df['cycle'] = df['cycle'].astype(int)
        return df

    def _select_cycle(self, df):
        if 'cycle' not in df.columns: return df
        cycles = sorted(df['cycle'].unique())
        if self.cycle == 'last':   return df[df['cycle'] == cycles[-1]].copy()
        elif self.cycle == 'mean': return df.groupby('well')[['parallel','perpendicular']].mean().reset_index()
        else:
            if int(self.cycle) not in cycles: raise DataError(f"Cycle {self.cycle} not in {cycles}")
            return df[df['cycle'] == int(self.cycle)].copy()

    # ── Merge ──────────────────────────────────────────────────────
    def _merge(self, layout, data):
        lw, dw = set(layout['well']), set(data['well'])
        issues = []
        if dw - lw: issues.append(f"In data not layout: {sorted(dw-lw)}")
        if lw - dw: issues.append(f"In layout not data: {sorted(lw-dw)}")
        if issues: raise MergeError('Well mismatch:\n  ' + '\n  '.join(issues))
        return layout.merge(data, on='well', how='inner')

    # ── G-factor ───────────────────────────────────────────────────
    def _compute_G(self):
        if self._G_override is not None:
            self.G = self._G_override; self.G_se = 0.0; return
        fluo = self.df_merged[self.df_merged['sample_type'] == 'fluorophore']
        if fluo.empty:
            warnings.warn("No fluorophore wells -- G set to 1.0", UserWarning)
            self.G = 1.0; self.G_se = 0.0; return
        A = self.A_ref
        G_vals    = fluo['parallel'] * (1-A) / (fluo['perpendicular'] * (1+2*A))
        self.G    = float(G_vals.mean())
        self.G_se = float(G_vals.sem()) if len(G_vals) > 1 else 0.0

    # ── Anisotropy ─────────────────────────────────────────────────
    def _compute_anisotropy(self):
        df = self.df_merged; G = self.G
        blank = (df[df['sample_type'] == 'blank']
                 .groupby('condition')[['parallel','perpendicular']].mean()
                 .rename(columns={'parallel':'blank_prl','perpendicular':'blank_prp'}))
        df = df.merge(blank, on='condition', how='left')
        df['blank_prl'] = df['blank_prl'].fillna(0.0)
        df['blank_prp'] = df['blank_prp'].fillna(0.0)
        df['Ipar_c']  = df['parallel']      - df['blank_prl']
        df['Iperp_c'] = df['perpendicular'] - df['blank_prp']
        denom = df['Ipar_c'] + G * df['Iperp_c']
        df['P'] = (df['Ipar_c'] - G*df['Iperp_c']) / denom.where(denom > 0)
        df['r'] = 2*df['P'] / (3 - df['P'])
        df['intensity'] = df['Ipar_c'] + 2*G*df['Iperp_c']
        self.df_merged = df

    # ── QC ─────────────────────────────────────────────────────────
    def _compute_qc(self):
        df = self.df_merged
        for cond, grp in df.groupby('condition'):
            warns = []
            G, G_se = self.G, self.G_se
            if abs(G - 1.0) > self.G_threshold:
                warns.append(f"G={G:.3f} deviates from 1 by >{self.G_threshold}")
            sub = grp[grp['sample_type'] == 'substrate']
            if sub.empty: sub = grp[grp['concentration'] == 0]
            A_free    = float(sub['r'].mean()) if not sub.empty else np.nan
            A_free_se = float(sub['r'].sem())  if len(sub) > 1  else 0.0
            tit = grp[grp['sample_type'] == 'titration'].dropna(subset=['concentration'])
            if not tit.empty:
                I_free  = tit.loc[tit['concentration']==tit['concentration'].min(),'intensity'].mean()
                I_bound = tit.loc[tit['concentration']==tit['concentration'].max(),'intensity'].mean()
                Q = I_bound/I_free if I_free and I_free > 0 else np.nan
            else:
                Q = np.nan
            if not np.isnan(Q) and abs(Q-1.0) > self.Q_threshold:
                warns.append(f"Q={Q:.3f} deviates from 1 by >{self.Q_threshold}")
            brows = grp[grp['sample_type'] == 'blank']
            b_prl = float(brows['parallel'].mean())      if not brows.empty else 0.0
            b_prp = float(brows['perpendicular'].mean()) if not brows.empty else 0.0
            n_tit  = tit['concentration'].nunique() if not tit.empty else 0
            n_reps = int(tit['replicate'].nunique()) if not tit.empty else 0
            for w in warns:
                warnings.warn(f"[{cond}] {w}", UserWarning, stacklevel=3)
            self.qc_reports[cond] = QCReport(
                condition=cond, G=G, G_se=G_se, Q=Q,
                A_free=A_free, A_free_se=A_free_se,
                blank_prl=b_prl, blank_prp=b_prp,
                n_titration=n_tit, n_replicates=n_reps, warnings=warns)

    # ── Saturation table (Scatchard now built post-fit, see get_scatchard) ──
    def _prepare_fitting_data(self):
        df  = self.df_merged
        tit = df[(df['sample_type']=='titration') & df['concentration'].notna()].copy()
        err_agg  = lambda x: ((x**2).sum()**0.5)/len(x)
        mean_agg = lambda x: x.mean()
        by_rep = (tit.groupby(['condition','concentration','replicate'], as_index=False)
                     .agg(r_mean=('r','mean'), r_sem=('r','sem'), r_std=('r','std')))
        sat = (by_rep.groupby(['condition','concentration'], as_index=False)
                     .agg(r=('r_mean',mean_agg), sem=('r_sem',err_agg),
                          std=('r_std',err_agg), n=('r_mean','count')))
        sat['r_mA']   = sat['r']   * 1000
        sat['sem_mA'] = sat['sem'] * 1000
        self.df_saturation = sat
        # df_scatchard is no longer precomputed; it is built post-fit on demand.
        self.df_scatchard = None

    def compute_z_prime(self, condition):
        """
        Z'-factor for assay quality (Zhang et al., 1999; Roehrl et al., 2004 eq 48).
            Z' = 1 - 3*(sigma_pos + sigma_neg) / |mu_pos - mu_neg|
        Positive control = titration wells at the highest [protein].
        Negative control = substrate-only wells.
        """
        df = self.df_merged
        cond_data = df[df['condition'] == condition]
        if cond_data.empty:
            raise FitError(f"Condition '{condition}' not found.")
        neg = cond_data[cond_data['sample_type'] == 'substrate']['r'].dropna()
        if len(neg) < 2:
            raise FitError(f"Z' for '{condition}': need >=2 substrate wells, found {len(neg)}.")
        tit = cond_data[(cond_data['sample_type'] == 'titration')
                        & cond_data['concentration'].notna()]
        if tit.empty:
            raise FitError(f"Z' for '{condition}': no titration wells found.")
        max_conc = tit['concentration'].max()
        pos = tit[tit['concentration'] == max_conc]['r'].dropna()
        if len(pos) < 2:
            raise FitError(f"Z' for '{condition}': need >=2 wells at max [P]={max_conc} nM, "
                           f"found {len(pos)}.")
        mu_pos, mu_neg       = float(pos.mean()), float(neg.mean())
        sigma_pos, sigma_neg = float(pos.std(ddof=1)), float(neg.std(ddof=1))
        gap = abs(mu_pos - mu_neg)
        z_prime = 1.0 - 3.0 * (sigma_pos + sigma_neg) / gap if gap > 0 else float('-inf')
        return {'condition': condition,
                'mu_pos': mu_pos, 'mu_neg': mu_neg,
                'sigma_pos': sigma_pos, 'sigma_neg': sigma_neg,
                'n_pos': int(len(pos)), 'n_neg': int(len(neg)),
                'Z_prime': z_prime,
                'pos_concentration_nM': float(max_conc)}

    # ── Fit (FIXED: A_free/A_bound now seeded for ALL models) ──────
    def fit(self, condition, model='fp_quadratic', params=None,
            bounds=None, fixed=None, method='leastsq'):
        if self.df_saturation is None or self.df_saturation.empty:
            raise FitError("No saturation data.")
        cond_data = self.df_saturation[self.df_saturation['condition'] == condition]
        if cond_data.empty:
            raise ValueError(f"Condition '{condition}' not found.")
        x = np.asarray(cond_data['concentration'].values, dtype=np.float64)
        y = np.asarray(cond_data['r'].values,             dtype=np.float64)
        if model not in _BUILTIN_MODELS:
            raise ValueError(f"Unknown model '{model}'.")
        model_func = _BUILTIN_MODELS[model]
        params = dict(params or {}); bounds = dict(bounds or {}); fixed = dict(fixed or {})
        qc = self.qc_reports.get(condition)

        # Sensible defaults derived from the data / QC, shared by all models
        A_free_def  = qc.A_free if (qc and not np.isnan(qc.A_free)) else float(y.min())
        A_bound_def = float(y.max())
        Q_def       = qc.Q if (qc and not np.isnan(qc.Q)) else 1.0

        if model == 'fp_quadratic':
            # Kd and A_bound float; Q, Dl, A_free fixed from QC unless overridden
            params.setdefault('Kd', 50.0)
            params.setdefault('A_bound', A_bound_def)
            fixed.setdefault('Q',      Q_def)
            fixed.setdefault('Dl',     self.fixed_concentration)
            fixed.setdefault('A_free', A_free_def)

        elif model == 'one_site':
            # Kd, A_bound float; A_free fixed from QC unless overridden
            params.setdefault('Kd', 50.0)
            params.setdefault('A_bound', A_bound_def)
            fixed.setdefault('A_free', A_free_def)

        elif model == 'hill':
            params.setdefault('Kd', 50.0)
            params.setdefault('A_bound', A_bound_def)
            params.setdefault('n', 1.0)
            fixed.setdefault('A_free', A_free_def)

        # If the caller fixed a param that also appears in params, drop it from params
        for k in list(params.keys()):
            if k in fixed:
                params.pop(k)

        # Sensible default bounds so the optimiser stays in physical territory.
        # Kd and n must stay positive; A_bound/A_free stay in anisotropy range.
        default_bounds = {
            'Kd':      (0.0, np.inf),
            'n':       (0.1, 10.0),
            'A_bound': (0.0, 1.0),
            'A_free':  (0.0, 1.0),
        }
        for k, db in default_bounds.items():
            bounds.setdefault(k, db)

        lm_model  = Model(model_func)
        lm_params = Parameters()
        for name, value in params.items():
            lo, hi = bounds.get(name, (-np.inf, np.inf))
            # Make sure the initial value sits inside the bounds
            if np.isfinite(lo): value = max(value, lo)
            if np.isfinite(hi): value = min(value, hi)
            lm_params.add(name, value=value, min=lo, max=hi, vary=True)
        for name, value in fixed.items():
            lm_params.add(name, value=float(value), vary=False)

        # Safety net: ensure every required model argument is present
        required = [p for p in model_func.__code__.co_varnames[:model_func.__code__.co_argcount]
                    if p != 'x']
        for p in required:
            if p not in lm_params:
                seed = {'A_free': A_free_def, 'A_bound': A_bound_def,
                        'Q': Q_def, 'Dl': self.fixed_concentration,
                        'Kd': 50.0, 'n': 1.0}.get(p, 1.0)
                lm_params.add(p, value=float(seed), vary=False)

        try:
            result = lm_model.fit(y, lm_params, x=x, method=method)
        except Exception as e:
            raise FitError(f"Fit failed for '{condition}': {e}") from e
        if not result.success:
            warnings.warn(f"[{condition}] Fit did not converge: {result.message}", UserWarning)
        self._fit_results[condition] = result
        x_fine = np.linspace(max(x.min(), 1e-9), x.max(), 400)
        y_fine = model_func(x_fine, **result.best_values)
        self._fit_curves[condition] = pd.DataFrame(
            {'condition':condition,'concentration':x_fine,'r':y_fine,'r_mA':y_fine*1000})

    # ── Post-fit Scatchard (FIXED: uses fitted A_bound and B/F) ────
    def get_scatchard(self, condition):
        """
        Build a post-fit Scatchard plot for one condition.

        This is a *diagnostic visualisation*, not an independent K_D estimate.
        It uses the fitted A_free and A_bound to convert each measured
        anisotropy into a bound-complex concentration [B], then plots B/F
        against B, where F = [L*]_0 - [B] is the free labelled-ligand
        concentration. A straight line is consistent with simple 1:1 binding;
        curvature suggests multiple sites, cooperativity, or non-specific
        binding.

        Returns a DataFrame: condition | concentration | B | F | B_over_F.
        Raises FitError if the condition has not been fitted yet.
        """
        if condition not in self._fit_results:
            raise FitError(f"No fit for '{condition}'. Run fit('{condition}', ...) first.")
        res = self._fit_results[condition]
        bv  = res.best_values
        A_free  = bv.get('A_free')
        A_bound = bv.get('A_bound')
        if A_free is None:
            A_free = self.qc_reports[condition].A_free if condition in self.qc_reports else None
        if A_free is None or A_bound is None:
            raise FitError(f"Cannot build Scatchard for '{condition}': "
                           f"missing fitted A_free / A_bound.")
        denom = A_bound - A_free
        if abs(denom) < 1e-12:
            raise FitError(f"A_bound - A_free ~ 0 for '{condition}'; Scatchard undefined.")

        Dl  = self.fixed_concentration
        sat = self.df_saturation[self.df_saturation['condition'] == condition]
        rows = []
        for _, row in sat.iterrows():
            frac = (row['r'] - A_free) / denom          # fraction bound (0..1)
            frac = float(np.clip(frac, 0.0, 1.0))
            B = frac * Dl                                # bound complex [nM]
            F = Dl - B                                   # free labelled ligand [nM]
            B_over_F = B / F if F > 1e-12 else np.nan
            rows.append({'condition': condition,
                         'concentration': row['concentration'],
                         'B': B, 'F': F, 'B_over_F': B_over_F})
        return pd.DataFrame(rows)

    # ── Accessors ──────────────────────────────────────────────────
    def fit_report(self, c):  return self._fit_results[c].fit_report()
    def get_fit_curve(self, c): return self._fit_curves[c].copy()
    def get_fit_params(self, c): return dict(self._fit_results[c].best_values)
    def get_fit_summary(self):
        rows = []
        for cond, res in self._fit_results.items():
            for pname, param in res.params.items():
                if param.vary:
                    rows.append({'condition':cond,'parameter':pname,
                                 'value':param.value,'stderr':param.stderr,'redchi':res.redchi})
        return pd.DataFrame(rows)
    def print_qc(self):
        print(f"G-factor (global): {self.G:.4f} +/- {self.G_se:.4f}\n")
        for r in self.qc_reports.values():
            print(r); print()
    def __repr__(self):
        n = self.df_saturation['condition'].nunique() if self.df_saturation is not None else 0
        return f"FPAnalyzer(conditions={n}, G={self.G:.3f}, Dl={self.fixed_concentration} nM, fits={list(self._fit_results)})"


# ══════════════════════════════════════════════════════════════════
# UI layer
# ══════════════════════════════════════════════════════════════════
COLORS = plt.cm.tab10.colors
_ana: Optional[FPAnalyzer] = None   # module-level current analyzer


def get_analyzer() -> Optional[FPAnalyzer]:
    """Return the FPAnalyzer instance created by the load panel."""
    return _ana


def _fig_to_widget(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=130, bbox_inches='tight')
    plt.close(fig); buf.seek(0)
    return widgets.Image(value=buf.read(), format='png')

def _save_bytes(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig); buf.seek(0)
    return buf.read()

def _styled_btn(desc, style='primary', icon=''):
    return widgets.Button(description=desc, button_style=style,
                          icon=icon, layout=widgets.Layout(width='auto'))

def _label(txt, bold=False):
    s = f'<b>{txt}</b>' if bold else txt
    return widgets.HTML(f'<span style="font-size:13px">{s}</span>')

def _condition_colors():
    if _ana is None: return {}
    conds = sorted(_ana.df_saturation['condition'].unique())
    return {c: COLORS[i % 10] for i, c in enumerate(conds)}


# ── Step 2: Load ───────────────────────────────────────────────────
def show_load_panel():
    upload_layout = widgets.FileUpload(accept='.csv,.xlsx', multiple=False,
                                       description='Layout file')
    upload_data   = widgets.FileUpload(accept='.csv,.xlsx,.txt', multiple=False,
                                       description='Data file')
    w_Dl    = widgets.FloatText(value=10.0, description='[L*] nM:',
                                style={'description_width':'80px'},
                                layout=widgets.Layout(width='200px'))
    w_Aref  = widgets.FloatText(value=0.035, description='A_ref:',
                                style={'description_width':'80px'},
                                layout=widgets.Layout(width='200px'))
    w_G     = widgets.FloatText(value=1.0, description='G override:',
                                style={'description_width':'90px'},
                                layout=widgets.Layout(width='200px'))
    w_G_use = widgets.Checkbox(value=False, description='Use G override',
                               layout=widgets.Layout(width='180px'))
    w_cycle = widgets.Dropdown(
        options=[('Last cycle','last'),('Mean of all','mean'),
                 ('Cycle 1',1),('Cycle 2',2),('Cycle 3',3),
                 ('Cycle 5',5),('Cycle 10',10)],
        value='last', description='Cycle:',
        style={'description_width':'80px'},
        layout=widgets.Layout(width='220px'))
    btn_load = _styled_btn('Load & process', 'success', 'upload')
    status   = widgets.Output()
    preview  = widgets.Output()

    def _save_upload(upload):
        if not upload.value: return None
        item = upload.value[0]
        fname   = item['name']    if isinstance(item, dict) else item.name
        content = item['content'] if isinstance(item, dict) else item.content
        content = bytes(content)
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=Path(fname).suffix)
        tmp.write(content); tmp.flush()
        return tmp.name

    def on_load(b):
        global _ana
        with status: clear_output(wait=True); display(_label('Loading...'))
        with preview: clear_output()
        try:
            lp = _save_upload(upload_layout)
            dp = _save_upload(upload_data)
            if lp is None: raise ValueError('Attach a layout file.')
            if dp is None: raise ValueError('Attach a data file.')
            G_val = w_G.value if w_G_use.value else None
            _ana = FPAnalyzer(layout_file=lp, data_file=dp,
                              fixed_concentration=w_Dl.value,
                              A_ref=w_Aref.value, G=G_val,
                              cycle=w_cycle.value)
            n_c = _ana.df_saturation['condition'].nunique()
            n_w = len(_ana.df_merged)
            with status:
                clear_output(wait=True)
                display(widgets.HTML(
                    f'<span style="color:green;font-size:13px">'
                    f'Loaded: {n_w} wells, {n_c} conditions. '
                    f'G = {_ana.G:.4f} +/- {_ana.G_se:.4f}</span>'))
            with preview:
                df = _ana.df_merged
                type_colors = {'titration':'#4878d0','substrate':'#6acc65',
                               'blank':'#d65f5f','fluorophore':'#ee854a'}
                fig, ax = plt.subplots(figsize=(10, 3.5))
                for _, row in df.iterrows():
                    w = row['well']
                    r_idx = ord(w[0]) - ord('A')
                    c_idx = int(w[1:]) - 1
                    col = type_colors.get(row['sample_type'], 'gray')
                    ax.add_patch(plt.Rectangle((c_idx, -r_idx), 0.9, 0.9, color=col, alpha=0.7))
                ax.set_xlim(-0.2, 24.2); ax.set_ylim(-16.5, 1.2)
                ax.set_xticks(np.arange(24)+0.45)
                ax.set_xticklabels(range(1,25), fontsize=7)
                ax.set_yticks([-i for i in range(16)])
                ax.set_yticklabels(list('ABCDEFGHIJKLMNOP'), fontsize=7)
                ax.set_title('Layout — plate map'); ax.axis('off')
                from matplotlib.patches import Patch
                legend = [Patch(color=v, label=k) for k,v in type_colors.items()]
                ax.legend(handles=legend, fontsize=8, loc='upper right')
                display(_fig_to_widget(fig))
        except Exception as e:
            with status:
                clear_output(wait=True)
                display(widgets.HTML(f'<span style="color:red">Error: {e}</span>'))

    btn_load.on_click(on_load)
    display(widgets.VBox([
        _label('Attach files and set parameters', bold=True),
        widgets.HBox([upload_layout, upload_data]),
        widgets.HTML('<hr style="margin:6px 0">'),
        _label('Experiment parameters'),
        widgets.HBox([w_Dl, w_Aref]),
        widgets.HBox([w_G, w_G_use, w_cycle]),
        widgets.HTML('<hr style="margin:6px 0">'),
        btn_load, status, preview,
    ]))


# ── Step 3: QC ─────────────────────────────────────────────────────
def show_qc_panel():
    btn_qc  = _styled_btn('Show QC', 'info', 'search')
    w_logx  = widgets.Checkbox(value=True, description='Log X axis',
                               layout=widgets.Layout(width='140px'))
    out_qc  = widgets.Output()

    def on_qc(b):
        with out_qc:
            clear_output(wait=True)
            if _ana is None:
                display(_label('Load data first (Step 2).')); return
            txt = f'G-factor (global): {_ana.G:.4f} +/- {_ana.G_se:.4f}\n\n'
            for r in _ana.qc_reports.values():
                txt += str(r) + '\n\n'
            display(widgets.Textarea(value=txt, rows=18,
                                      layout=widgets.Layout(width='98%')))
            # Z'-factor
            zp_rows = []
            for cond in sorted(_ana.df_saturation['condition'].unique()):
                try:
                    zp_rows.append(_ana.compute_z_prime(cond))
                except FitError:
                    zp_rows.append({'condition': cond, 'Z_prime': float('nan'),
                                    'mu_pos': float('nan'), 'mu_neg': float('nan'),
                                    'sigma_pos': float('nan'), 'sigma_neg': float('nan'),
                                    'n_pos': 0, 'n_neg': 0,
                                    'pos_concentration_nM': float('nan')})
            zp_df = pd.DataFrame(zp_rows)
            if not zp_df.empty:
                display(_label("Z'-factor:", bold=True))
                disp = zp_df[['condition','Z_prime','mu_pos','mu_neg',
                              'sigma_pos','sigma_neg','n_pos','n_neg',
                              'pos_concentration_nM']].copy()
                disp['Z_prime'] = disp['Z_prime'].round(3)
                for c in ['mu_pos','mu_neg','sigma_pos','sigma_neg']:
                    disp[c] = disp[c].round(4)
                disp['pos_concentration_nM'] = disp['pos_concentration_nM'].round(1)
                disp = disp.rename(columns={'pos_concentration_nM':'[P]_max nM'})
                display(disp)
            # G calibration
            fluo = _ana.df_merged[_ana.df_merged['sample_type'] == 'fluorophore']
            if not fluo.empty:
                A = _ana.A_ref
                G_vals = fluo['parallel']*(1-A) / (fluo['perpendicular']*(1+2*A))
                fig, ax = plt.subplots(figsize=(5, 3))
                ax.scatter(range(len(G_vals)), G_vals, zorder=3, color='#4878d0')
                ax.axhline(_ana.G, color='red', ls='--', label=f'G mean = {_ana.G:.4f}')
                ax.axhspan(_ana.G-_ana.G_se, _ana.G+_ana.G_se, alpha=0.15, color='red', label='+/-SE')
                ax.set(xlabel='Fluorophore well #', ylabel='G per well',
                       title='G-factor calibration')
                ax.legend(fontsize=8)
                display(_fig_to_widget(fig))
            # Titration overview (respects log/linear checkbox)
            cm = _condition_colors()
            fig, ax = plt.subplots(figsize=(8, 3.5))
            for cond in sorted(_ana.df_saturation['condition'].unique()):
                sat = _ana.df_saturation[_ana.df_saturation['condition'] == cond]
                ax.errorbar(sat['concentration'], sat['r_mA'], yerr=sat['sem_mA'],
                            fmt='o-', color=cm[cond], capsize=3,
                            label=cond, markersize=5)
            ax.set(xlabel='[Protein], nM', ylabel='Anisotropy, mA',
                   title='Titration curves (raw)')
            if w_logx.value:
                ax.set_xscale('log')
            ax.legend(fontsize=8)
            display(_fig_to_widget(fig))

    btn_qc.on_click(on_qc)
    display(widgets.VBox([widgets.HBox([btn_qc, w_logx]), out_qc]))


# ── Step 4: Fit ────────────────────────────────────────────────────
def show_fit_panel():
    mkw = lambda desc, val, w='180px': widgets.FloatText(
        value=val, description=desc,
        style={'description_width':'110px'}, layout=widgets.Layout(width=w))
    mkd = lambda desc, opts, val, w='220px': widgets.Dropdown(
        options=opts, value=val, description=desc,
        style={'description_width':'90px'}, layout=widgets.Layout(width=w))
    mkc = lambda desc, val: widgets.Checkbox(value=val, description=desc,
                                              layout=widgets.Layout(width='160px'))

    cond_sel   = mkd('Condition:', [], None, '240px')
    model_sel  = mkd('Model:', list(_BUILTIN_MODELS.keys()), 'fp_quadratic')
    method_sel = mkd('Method:', ['leastsq','least_squares','nelder','powell'], 'leastsq', '200px')
    w_logx     = mkc('Log X axis', True)
    w_Kd_init  = mkw('Kd init (nM):', 50.)
    w_Kd_min   = mkw('Kd min:', 0.)
    w_Kd_max   = mkw('Kd max:', 1e4)
    w_Ab_init  = mkw('A_bound init:', 0.20, '200px')
    w_Ab_fix   = mkc('Fix A_bound', False)
    w_Ab_min   = mkw('A_bound min:', 0., '200px')
    w_Ab_max   = mkw('A_bound max:', 0.5, '200px')
    w_nH_val   = mkw('Hill n init:', 1.0)
    w_nH_fix   = mkc('Fix n (Hill only)', False)
    w_nH_min   = mkw('Hill n min:', 0.1)
    w_nH_max   = mkw('Hill n max:', 5.0)
    btn_fit    = _styled_btn('Run fit', 'primary', 'play')
    btn_fitall = _styled_btn('Fit all conditions', 'warning', 'forward')
    out_fit    = widgets.Output()

    def refresh_conds():
        if _ana is not None:
            cond_sel.options = sorted(_ana.df_saturation['condition'].unique())
            if cond_sel.value is None and len(cond_sel.options):
                cond_sel.value = cond_sel.options[0]

    def _do_fit(condition):
        params = {'Kd': w_Kd_init.value}
        bounds = {'Kd': (w_Kd_min.value, w_Kd_max.value)}
        fixed  = {}
        if w_Ab_fix.value:
            fixed['A_bound'] = w_Ab_init.value
        else:
            params['A_bound'] = w_Ab_init.value
            bounds['A_bound'] = (w_Ab_min.value, w_Ab_max.value)
        if model_sel.value == 'hill':
            if w_nH_fix.value:
                fixed['n'] = w_nH_val.value
            else:
                params['n'] = w_nH_val.value
                bounds['n'] = (w_nH_min.value, w_nH_max.value)
        _ana.fit(condition=condition, model=model_sel.value,
                 params=params, bounds=bounds, fixed=fixed,
                 method=method_sel.value)

    def _show_result(cond):
        res   = _ana._fit_results[cond]
        col   = _condition_colors().get(cond, 'steelblue')
        sat   = _ana.df_saturation[_ana.df_saturation['condition'] == cond]
        curve = _ana.get_fit_curve(cond)
        lines = [f'Fit converged  |  redchi = {res.redchi:.5f}']
        for pname, param in res.params.items():
            if param.vary:
                se = f'+/- {param.stderr:.4f}' if param.stderr else '(no SE)'
                lines.append(f'  {pname:12s} = {param.value:.4f}  {se}')
        display(widgets.Textarea(value='\n'.join(lines), rows=6,
                                  layout=widgets.Layout(width='98%')))
        Kd    = res.params['Kd'].value
        Kd_se = res.params['Kd'].stderr or 0.
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        axes[0].errorbar(sat['concentration'], sat['r_mA'], yerr=sat['sem_mA'],
                         fmt='o', color=col, capsize=4, label='data', zorder=3)
        axes[0].plot(curve['concentration'], curve['r_mA'], '-', color=col,
                     label=f'fit  Kd = {Kd:.1f} +/- {Kd_se:.1f} nM')
        axes[0].set(xlabel='[Protein], nM', ylabel='Anisotropy, mA',
                    title=f'{cond} -- saturation')
        if w_logx.value:
            axes[0].set_xscale('log')
        axes[0].legend(fontsize=8)
        model_func = _BUILTIN_MODELS[model_sel.value]
        r_pred = model_func(np.asarray(sat['concentration'], dtype=np.float64),
                            **res.best_values)
        resid  = (np.asarray(sat['r'], dtype=np.float64) - r_pred) * 1000
        axes[1].bar(range(len(resid)), resid, color=col, alpha=0.7)
        axes[1].axhline(0, color='k', lw=0.8)
        axes[1].set(xlabel='Point #', ylabel='Residual (mA)', title='Residuals')
        plt.tight_layout()
        display(_fig_to_widget(fig))

    def on_fit(b):
        with out_fit:
            clear_output(wait=True)
            if _ana is None: display(_label('Load data first.')); return
            refresh_conds()
            try:
                _do_fit(cond_sel.value); _show_result(cond_sel.value)
            except Exception as e:
                display(widgets.HTML(f'<span style="color:red">Error: {e}</span>'))

    def on_fit_all(b):
        with out_fit:
            clear_output(wait=True)
            if _ana is None: display(_label('Load data first.')); return
            refresh_conds()
            for cond in sorted(_ana.df_saturation['condition'].unique()):
                try:
                    _do_fit(cond); display(_label(f'OK: {cond}'))
                except Exception as e:
                    display(_label(f'FAIL {cond}: {e}'))
            display(_label('Done. See Step 5 for the summary.', bold=True))

    btn_fit.on_click(on_fit)
    btn_fitall.on_click(on_fit_all)
    refresh_conds()
    display(widgets.VBox([
        _label('Fitting parameters', bold=True),
        widgets.HBox([cond_sel, model_sel, method_sel]),
        widgets.HBox([w_logx]),
        widgets.HTML('<hr style="margin:5px 0">'),
        _label('Kd'),
        widgets.HBox([w_Kd_init, w_Kd_min, w_Kd_max]),
        _label('A_bound'),
        widgets.HBox([w_Ab_init, w_Ab_fix, w_Ab_min, w_Ab_max]),
        _label('Hill coefficient (used only for the Hill model)'),
        widgets.HBox([w_nH_val, w_nH_fix, w_nH_min, w_nH_max]),
        widgets.HTML('<hr style="margin:5px 0">'),
        widgets.HBox([btn_fit, btn_fitall]),
        out_fit,
    ]))


# ── Step 5: Summary ────────────────────────────────────────────────
def show_summary_panel():
    btn_sum     = _styled_btn('Refresh', 'info', 'refresh')
    w_logx      = widgets.Checkbox(value=True, description='Log X axis',
                                   layout=widgets.Layout(width='140px'))
    scat_sel    = widgets.Dropdown(options=[], description='Scatchard:',
                                   style={'description_width':'90px'},
                                   layout=widgets.Layout(width='240px'))
    out_sum     = widgets.Output()

    def on_refresh(b):
        with out_sum:
            clear_output(wait=True)
            if _ana is None:
                display(_label('Load data first.')); return
            fitted = list(_ana._fit_results.keys())
            if not fitted:
                display(_label('No fits yet — run Step 4.')); return
            scat_sel.options = fitted
            if scat_sel.value not in fitted and fitted:
                scat_sel.value = fitted[0]

            df_sum  = _ana.get_fit_summary()
            kd_rows = df_sum[df_sum['parameter'] == 'Kd'].copy()
            if not kd_rows.empty:
                kd_rows['SE_%'] = (100 * kd_rows['stderr'] / kd_rows['value']).round(1)
                display(_label('Kd summary:', bold=True))
                display(kd_rows[['condition','value','stderr','SE_%','redchi']]
                        .rename(columns={'value':'Kd (nM)','stderr':'SE (nM)'})
                        .reset_index(drop=True))

            cm = _condition_colors()
            # Saturation curves (all fitted conditions)
            fig, ax = plt.subplots(figsize=(7, 5))
            for cond in fitted:
                col   = cm.get(cond, 'gray')
                sat   = _ana.df_saturation[_ana.df_saturation['condition'] == cond]
                curve = _ana.get_fit_curve(cond)
                Kd    = _ana.get_fit_params(cond)['Kd']
                ax.errorbar(sat['concentration'], sat['r_mA'], yerr=sat['sem_mA'],
                            fmt='o', color=col, capsize=3, markersize=5)
                ax.plot(curve['concentration'], curve['r_mA'], '-',
                        color=col, label=f'{cond}  Kd={Kd:.1f} nM')
            ax.set(xlabel='[Protein], nM', ylabel='Anisotropy, mA',
                   title='Saturation -- all fitted conditions')
            if w_logx.value:
                ax.set_xscale('log')
            ax.legend(fontsize=8)
            plt.tight_layout()
            display(_fig_to_widget(fig))

            # Post-fit Scatchard for the selected condition
            cond = scat_sel.value
            try:
                sc = _ana.get_scatchard(cond)
                col = cm.get(cond, 'steelblue')
                fig, ax = plt.subplots(figsize=(5, 4))
                ax.plot(sc['B'], sc['B_over_F'], 'o-', color=col)
                ax.set(xlabel='Bound complex B, nM', ylabel='B / F',
                       title=f'Post-fit Scatchard — {cond}')
                ax.text(0.02, 0.02,
                        'Diagnostic only (uses fitted A_bound/A_free; not an independent Kd)',
                        transform=ax.transAxes, fontsize=7, color='gray')
                plt.tight_layout()
                display(_fig_to_widget(fig))
            except FitError as e:
                display(_label(f'Scatchard unavailable: {e}'))

    btn_sum.on_click(on_refresh)
    scat_sel.observe(lambda ch: on_refresh(None) if ch['name'] == 'value' else None, names='value')
    display(widgets.VBox([widgets.HBox([btn_sum, w_logx, scat_sel]), out_sum]))


# ── Step 6: Export ─────────────────────────────────────────────────
def show_export_panel():
    chk_sat      = widgets.Checkbox(value=True, description='Saturation data (CSV)')
    chk_scat     = widgets.Checkbox(value=True, description='Scatchard data (CSV)')
    chk_params   = widgets.Checkbox(value=True, description='Fit parameters (CSV)')
    chk_qc       = widgets.Checkbox(value=True, description='QC report (TXT)')
    chk_zp       = widgets.Checkbox(value=True, description="Z'-factor (CSV)")
    chk_fig_sat  = widgets.Checkbox(value=True, description='Saturation plot (PNG)')
    chk_fig_scat = widgets.Checkbox(value=True, description='Scatchard plots (PNG)')
    chk_fig_qc   = widgets.Checkbox(value=True, description='G calibration (PNG)')
    w_logx       = widgets.Checkbox(value=True, description='Log X on saturation')
    btn_exp = _styled_btn('Export ZIP', 'warning', 'download')
    out_exp = widgets.Output()

    def on_export(b):
        with out_exp:
            clear_output(wait=True)
            if _ana is None:
                display(_label('Load data first.')); return
            zip_buf = io.BytesIO()
            cm = _condition_colors()
            fitted = list(_ana._fit_results.keys())
            with zipfile.ZipFile(zip_buf, 'w', zipfile.ZIP_DEFLATED) as zf:
                if chk_sat.value and _ana.df_saturation is not None:
                    zf.writestr('saturation_data.csv', _ana.df_saturation.to_csv(index=False))
                if chk_scat.value and fitted:
                    scat_all = pd.concat([_ana.get_scatchard(c) for c in fitted],
                                         ignore_index=True)
                    zf.writestr('scatchard_data.csv', scat_all.to_csv(index=False))
                if chk_params.value and fitted:
                    zf.writestr('fit_parameters.csv', _ana.get_fit_summary().to_csv(index=False))
                    reports = ''
                    for cond in fitted:
                        reports += f'=== {cond} ===\n' + _ana.fit_report(cond) + '\n\n'
                    zf.writestr('fit_reports.txt', reports)
                if chk_zp.value:
                    zp_rows = []
                    for cond in sorted(_ana.df_saturation['condition'].unique()):
                        try: zp_rows.append(_ana.compute_z_prime(cond))
                        except FitError: pass
                    if zp_rows:
                        zf.writestr('z_prime.csv', pd.DataFrame(zp_rows).to_csv(index=False))
                if chk_qc.value:
                    txt = f'G-factor: {_ana.G:.4f} +/- {_ana.G_se:.4f}\n\n'
                    for r in _ana.qc_reports.values():
                        txt += str(r) + '\n\n'
                    zf.writestr('qc_report.txt', txt)
                if chk_fig_sat.value:
                    fig, ax = plt.subplots(figsize=(8, 5))
                    all_conds = fitted if fitted else sorted(_ana.df_saturation['condition'].unique())
                    for cond in all_conds:
                        col = cm.get(cond,'gray')
                        sat = _ana.df_saturation[_ana.df_saturation['condition']==cond]
                        ax.errorbar(sat['concentration'], sat['r_mA'], yerr=sat['sem_mA'],
                                    fmt='o', color=col, capsize=3, markersize=5, label=f'{cond} data')
                        if cond in _ana._fit_results:
                            curve = _ana.get_fit_curve(cond)
                            Kd    = _ana.get_fit_params(cond)['Kd']
                            ax.plot(curve['concentration'], curve['r_mA'], '-',
                                    color=col, label=f'{cond} Kd={Kd:.1f} nM')
                    ax.set(xlabel='[Protein], nM', ylabel='Anisotropy, mA', title='Saturation binding')
                    if w_logx.value: ax.set_xscale('log')
                    ax.legend(fontsize=8); plt.tight_layout()
                    zf.writestr('saturation_plot.png', _save_bytes(fig))
                if chk_fig_scat.value and fitted:
                    for cond in fitted:
                        try:
                            sc = _ana.get_scatchard(cond)
                        except FitError:
                            continue
                        fig, ax = plt.subplots(figsize=(5, 4))
                        ax.plot(sc['B'], sc['B_over_F'], 'o-', color=cm.get(cond,'gray'))
                        ax.set(xlabel='Bound complex B, nM', ylabel='B / F',
                               title=f'Post-fit Scatchard — {cond}')
                        plt.tight_layout()
                        zf.writestr(f'scatchard_{cond}.png', _save_bytes(fig))
                if chk_fig_qc.value:
                    fluo = _ana.df_merged[_ana.df_merged['sample_type']=='fluorophore']
                    if not fluo.empty:
                        A = _ana.A_ref
                        G_vals = fluo['parallel']*(1-A)/(fluo['perpendicular']*(1+2*A))
                        fig, ax = plt.subplots(figsize=(5, 3.5))
                        ax.scatter(range(len(G_vals)), G_vals, color='#4878d0')
                        ax.axhline(_ana.G, color='red', ls='--', label=f'G = {_ana.G:.4f}')
                        ax.axhspan(_ana.G-_ana.G_se, _ana.G+_ana.G_se, alpha=0.15, color='red')
                        ax.set(xlabel='Well #', ylabel='G', title='G-factor calibration')
                        ax.legend(fontsize=8); plt.tight_layout()
                        zf.writestr('g_calibration.png', _save_bytes(fig))
            zip_buf.seek(0)
            zb    = zip_buf.read()
            b64   = base64.b64encode(zb).decode()
            kb    = len(zb)//1024
            html  = (f'<a href="data:application/zip;base64,{b64}" '
                     f'download="fp_results.zip" '
                     f'style="font-size:14px;font-weight:bold;">'
                     f'Download fp_results.zip ({kb} KB)</a>')
            display(widgets.HTML(html))

    btn_exp.on_click(on_export)
    display(widgets.VBox([
        _label('Select what to export:', bold=True),
        _label('Tables'),
        widgets.HBox([chk_sat, chk_scat, chk_params, chk_qc, chk_zp]),
        _label('Figures'),
        widgets.HBox([chk_fig_sat, chk_fig_scat, chk_fig_qc]),
        widgets.HBox([w_logx]),
        widgets.HTML('<hr style="margin:6px 0">'),
        btn_exp, out_exp,
    ]))


# ── Full tabbed interface (optional convenience) ───────────────────
def build_ui():
    """Display all six steps as a single tabbed widget."""
    out_load, out_qc, out_fit, out_sum, out_exp = (widgets.Output() for _ in range(5))
    with out_load: show_load_panel()
    with out_qc:   show_qc_panel()
    with out_fit:  show_fit_panel()
    with out_sum:  show_summary_panel()
    with out_exp:  show_export_panel()
    tab = widgets.Tab(children=[out_load, out_qc, out_fit, out_sum, out_exp])
    for i, t in enumerate(['Load','QC','Fit','Summary','Export']):
        tab.set_title(i, t)
    display(tab)
