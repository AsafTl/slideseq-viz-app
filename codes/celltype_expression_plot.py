"""
Cell-type expression renderer for the Slideseq viewer.

Two outputs:

1. Dotplot of expression per cell type (genes x RCTD cell types).
   - dot size  = fraction of singlet beads of that cell type with UMI > 0
   - dot color = mean UMI of that gene in that cell type
   Available per-puck and as a cumulative pooled dotplot across a chosen
   subset of pucks (raw bead pooling, no per-puck averaging).

2. Spatial highlight map. Singlet beads are colored by RCTD cell type, with
   beads matching `(first_type == chosen) & (gene UMI > 0)` rendered on top
   with full alpha + black ring + larger marker.

Reuses `multigene_neighborhood_plot._load_puck_data` for the count matrix,
gene index, and bead coordinates (already aligned to a row-ordered barcode
list). Only the small `_results.csv` is loaded fresh, cached per puck.
"""

from __future__ import annotations

import importlib.util
import threading
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from matplotlib_scalebar.scalebar import ScaleBar
import numpy as np
import pandas as pd

APP_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = APP_DIR.parent
WORKING_DIR = PROJECT_ROOT / "working"
OUTPUT_DIR = PROJECT_ROOT / "figs"
PUCK_FOLDER_PREFIX = "2024-04-22_"
RESULTS_FILE_PREFIX = "2024-04-22_"

# Shared filename date prefix for this view's PNG cache.
TODAY = "2026-04-28"

# Import the sibling renderer via importlib (its name is fine, but we want
# the same instance the Flask app loaded so the in-memory data cache is
# shared). The Flask app sets sys.modules["multigene_plot_v2"] = module
# at import time, so prefer that if available; otherwise load fresh.
import sys
_NEIGHBOR_MOD_NAME = "multigene_plot_v2"
if _NEIGHBOR_MOD_NAME in sys.modules:
    multigene_plot = sys.modules[_NEIGHBOR_MOD_NAME]
else:
    _spec = importlib.util.spec_from_file_location(
        _NEIGHBOR_MOD_NAME, APP_DIR / "multigene_neighborhood_plot.py"
    )
    multigene_plot = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(multigene_plot)
    sys.modules[_NEIGHBOR_MOD_NAME] = multigene_plot

PUCK_IDS = multigene_plot.PUCK_IDS
MICRONS_PER_PIXEL = multigene_plot.MICRONS_PER_PIXEL

# Stable 12-cell-type palette. Order matches the cell-type vocabulary in the
# RCTD results across all 6 melanoma pucks. Picked for distinguishability on
# a dark background (the existing puck tiles use a #1a1814 background).
CELL_TYPE_PALETTE: dict[str, str] = {
    "Melanoma":            "#e6194b",
    "Proliferating Cells": "#f58231",
    "Fibroblasts":         "#ffe119",
    "Macrophages":         "#3cb44b",
    "T cells":             "#4363d8",
    "B cells":             "#911eb4",
    "Endothelial Cells":   "#42d4f4",
    "Epithelial Cells":    "#f032e6",
    "Erythrocytes":        "#9a6324",
    "Neutrophils":         "#bcf60c",
    "Neurons":             "#fabebe",
    "Metaphocytes":        "#aaffc3",
}
UNASSIGNED_LABEL = "Unassigned"
UNASSIGNED_COLOR = "#7a746c"


# ---------- Per-puck cell-type call cache ----------

# np.ndarray[str], length = n_beads in puck (matrix row order).
_CALLS_CACHE: dict[str, np.ndarray] = {}
_CALLS_LOCKS: dict[str, threading.Lock] = {p: threading.Lock() for p in PUCK_IDS}


def _results_path(puck_id: str) -> Path:
    return (
        WORKING_DIR
        / f"{PUCK_FOLDER_PREFIX}{puck_id}"
        / f"{RESULTS_FILE_PREFIX}{puck_id}_results.csv"
    )


def _load_celltype_calls(puck_id: str) -> np.ndarray:
    """Singlet `first_type` per bead, aligned to `_PuckData.barcodes`.

    Doublets and rejects get `UNASSIGNED_LABEL`. Returns a fixed-length
    string array indexable by matrix row.
    """
    cached = _CALLS_CACHE.get(puck_id)
    if cached is not None:
        return cached

    lock = _CALLS_LOCKS.setdefault(puck_id, threading.Lock())
    with lock:
        cached = _CALLS_CACHE.get(puck_id)
        if cached is not None:
            return cached

        data = multigene_plot._load_puck_data(puck_id)
        bead_barcodes = data.barcodes  # list[str]

        df = pd.read_csv(_results_path(puck_id))
        # The first column is unnamed in the CSV; pandas gives it a name
        # like "Unnamed: 0". Normalize to "barcode".
        first_col = df.columns[0]
        df = df.rename(columns={first_col: "barcode"})
        # Keep only what we need.
        df = df[["barcode", "spot_class", "first_type"]]

        bc_to_call = {
            b: t
            for b, sc, t in zip(df["barcode"], df["spot_class"], df["first_type"])
            if sc == "singlet" and isinstance(t, str)
        }

        calls = np.array(
            [bc_to_call.get(b, UNASSIGNED_LABEL) for b in bead_barcodes],
            dtype=object,
        )
        _CALLS_CACHE[puck_id] = calls
        return calls


def list_celltypes(puck_id: str | None = None) -> list[str]:
    """Sorted distinct singlet cell types present in a puck (or pooled).

    Defaults to the union across all 6 pucks so the UI can populate one
    consistent dropdown.
    """
    seen: set[str] = set()
    pucks = [puck_id] if puck_id else PUCK_IDS
    for pid in pucks:
        try:
            calls = _load_celltype_calls(pid)
        except Exception:
            continue
        for t in calls:
            if t and t != UNASSIGNED_LABEL:
                seen.add(t)
    return sorted(seen)


# ---------- Dotplot computation ----------

def _per_celltype_stats(
    matrix_csc, gene_col_idx: int, calls: np.ndarray, celltypes: list[str]
):
    """For a single gene column, compute (frac_expressing, mean_expr) per celltype.

    Both arrays have shape (len(celltypes),). Empty cell types yield NaN.
    """
    n_beads = matrix_csc.shape[0]
    col_vec = matrix_csc.getcol(gene_col_idx)
    nz_rows = col_vec.indices
    nz_vals = col_vec.data.astype(np.float64, copy=False)

    # Dense per-bead UMI vector (fast, since n_beads ~ 26k).
    umi = np.zeros(n_beads, dtype=np.float64)
    umi[nz_rows] = nz_vals

    fracs = np.full(len(celltypes), np.nan, dtype=np.float64)
    means = np.full(len(celltypes), np.nan, dtype=np.float64)
    for i, ct in enumerate(celltypes):
        mask = calls == ct
        n = int(mask.sum())
        if n == 0:
            continue
        sub = umi[mask]
        fracs[i] = float((sub > 0).mean())
        means[i] = float(sub.mean())
    return fracs, means


def _dotplot_data(
    pucks: list[str], genes: list[str], celltypes: list[str]
):
    """Pool beads across pucks and compute (frac, mean) for each (gene, ct).

    Returns:
        frac:  np.ndarray shape (len(genes), len(celltypes))
        mean:  np.ndarray shape (len(genes), len(celltypes))
        n_per_ct: np.ndarray shape (len(celltypes),) - pooled bead counts
        missing_genes: list[str] of gene symbols not found in ANY puck
    """
    # Build per-puck (calls, matrix, gene_col_idx_for_each_gene).
    n_genes = len(genes)
    n_cts = len(celltypes)

    # Pooled UMI vectors (concatenate across pucks) and pooled cell-type calls.
    pooled_calls_chunks: list[np.ndarray] = []
    # For each gene, accumulate the umi vector per puck and concat at the end.
    pooled_umi_chunks: list[list[np.ndarray]] = [[] for _ in genes]
    gene_found_anywhere = [False] * n_genes

    for pid in pucks:
        data = multigene_plot._load_puck_data(pid)
        calls = _load_celltype_calls(pid)
        n = data.matrix_csc.shape[0]
        pooled_calls_chunks.append(calls)

        for gi, g in enumerate(genes):
            col = data.gene_lower_to_col.get(g.lower())
            if col is None:
                pooled_umi_chunks[gi].append(np.zeros(n, dtype=np.float64))
                continue
            gene_found_anywhere[gi] = True
            col_vec = data.matrix_csc.getcol(col)
            umi = np.zeros(n, dtype=np.float64)
            umi[col_vec.indices] = col_vec.data
            pooled_umi_chunks[gi].append(umi)

    pooled_calls = np.concatenate(pooled_calls_chunks)
    pooled_umis = [np.concatenate(chunks) for chunks in pooled_umi_chunks]

    frac = np.full((n_genes, n_cts), np.nan, dtype=np.float64)
    mean = np.full((n_genes, n_cts), np.nan, dtype=np.float64)
    n_per_ct = np.zeros(n_cts, dtype=np.int64)

    for ci, ct in enumerate(celltypes):
        ct_mask = pooled_calls == ct
        nct = int(ct_mask.sum())
        n_per_ct[ci] = nct
        if nct == 0:
            continue
        for gi in range(n_genes):
            sub = pooled_umis[gi][ct_mask]
            expressing = sub > 0
            frac[gi, ci] = float(expressing.mean())
            # Mean UMI restricted to expressing beads. Slide-seq is sparse
            # enough that mean-over-all-beads is dominated by zeros and
            # collapses the color axis; separating "% expressing" (size)
            # from "intensity when expressed" (color) is more informative.
            if expressing.any():
                mean[gi, ci] = float(sub[expressing].mean())
            # else: leave as NaN, drawn as empty cell.

    missing = [g for g, ok in zip(genes, gene_found_anywhere) if not ok]
    return frac, mean, n_per_ct, missing


def _draw_dotplot(ax, genes, celltypes, frac, mean, *, title: str):
    """Render a dotplot onto an existing axes. Returns the size scatter handle
    so the caller can add a colorbar via the figure."""
    n_genes = len(genes)
    n_cts = len(celltypes)

    # Scale fraction → marker area (s in scatter is in points^2).
    # Empty cells (NaN frac) get drawn as a tiny outline so the user can see
    # the cell exists but has no data.
    s_min, s_max = 8.0, 220.0
    sizes = np.where(np.isnan(frac), 0.0, frac) * (s_max - s_min) + s_min
    sizes[np.isnan(frac)] = 4.0

    # Color = mean expression. Use vmax = 95th percentile of finite values
    # so a single high outlier doesn't crush the contrast.
    finite = mean[np.isfinite(mean)]
    vmax = float(np.percentile(finite, 95)) if finite.size else 1.0
    if vmax <= 0:
        vmax = 1.0

    # Build coordinate grid.
    gx, cy = np.meshgrid(np.arange(n_genes), np.arange(n_cts), indexing="ij")
    sc = ax.scatter(
        gx.ravel(),
        cy.ravel(),
        s=sizes.ravel(),
        c=np.where(np.isnan(mean), 0.0, mean).ravel(),
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
        edgecolors="black",
        linewidths=0.4,
        alpha=0.95,
    )

    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(genes, rotation=30, ha="right", fontsize=10)
    ax.set_yticks(range(n_cts))
    ax.set_yticklabels(celltypes, fontsize=9)
    ax.set_xlim(-0.6, n_genes - 0.4)
    ax.set_ylim(-0.6, n_cts - 0.4)
    ax.invert_yaxis()
    ax.grid(True, linestyle=":", color="0.85", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=11)

    return sc, vmax, (s_min, s_max)


def _add_dotplot_legends(fig, sc, vmax: float, size_range: tuple[float, float]):
    """Add a color bar (mean expression) and a size legend (% expressing)."""
    cbar = fig.colorbar(sc, ax=fig.axes, fraction=0.04, pad=0.04, shrink=0.7)
    cbar.set_label("mean UMI (expressing beads)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Manual size legend in a side axes — drawn as separate scatter points.
    s_min, s_max = size_range
    handles = []
    for f in (0.05, 0.25, 0.5, 0.85):
        s = f * (s_max - s_min) + s_min
        handles.append(
            fig.axes[0].scatter([], [], s=s, c="white", edgecolor="black",
                                linewidth=0.4, label=f"{int(f*100)}%")
        )
    fig.axes[0].legend(
        handles=handles,
        title="% expressing",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        labelspacing=1.2,
        borderpad=0.6,
        fontsize=8,
        title_fontsize=9,
    )


# ---------- Public dotplot renderers ----------

def _slug_genes(genes: list[str]) -> str:
    return "_".join(genes)


def _short_puck(pid: str) -> str:
    # "Puck_240208_24" -> "24"
    return pid.rsplit("_", 1)[-1]


def _slug_pucks(pucks: list[str]) -> str:
    return "-".join(sorted(_short_puck(p) for p in pucks))


def make_dotplot(puck_id: str, genes: list[str]) -> Path | None:
    if puck_id not in PUCK_IDS:
        raise ValueError(f"unknown puck {puck_id!r}")
    if not genes:
        return None

    cts = list(CELL_TYPE_PALETTE.keys())
    frac, mean, n_per_ct, missing = _dotplot_data([puck_id], genes, cts)

    fig = Figure(figsize=(max(5.0, 1.6 + 0.9 * len(genes)), 5.5),
                 constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    title = f"{puck_id}: expression per cell type"
    if missing:
        title += f"\n(missing genes: {', '.join(missing)})"
    sc, vmax, sr = _draw_dotplot(ax, genes, cts, frac, mean, title=title)
    _add_dotplot_legends(fig, sc, vmax, sr)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUTPUT_DIR / f"{TODAY}_{puck_id}_{_slug_genes(genes)}_dotplot.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    return out


def make_cumulative_dotplot(puck_ids: list[str], genes: list[str]) -> Path | None:
    pucks = [p for p in puck_ids if p in PUCK_IDS]
    if not pucks or not genes:
        return None

    cts = list(CELL_TYPE_PALETTE.keys())
    frac, mean, n_per_ct, missing = _dotplot_data(pucks, genes, cts)

    n_total = int(n_per_ct.sum())
    short_list = ", ".join(_short_puck(p) for p in sorted(pucks))
    title = f"Cumulative: pucks {short_list} (n={n_total:,} singlet beads)"
    if missing:
        title += f"\n(missing genes: {', '.join(missing)})"

    fig = Figure(figsize=(max(6.0, 2.0 + 1.0 * len(genes)), 5.5),
                 constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    sc, vmax, sr = _draw_dotplot(ax, genes, cts, frac, mean, title=title)
    _add_dotplot_legends(fig, sc, vmax, sr)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUTPUT_DIR / (
        f"{TODAY}_cumulative_{_slug_pucks(pucks)}_{_slug_genes(genes)}_dotplot.png"
    )
    fig.savefig(out, dpi=180, bbox_inches="tight")
    return out


# ---------- Spatial highlight map ----------

def _ct_slug(ct: str) -> str:
    return ct.replace(" ", "-")


def make_celltype_spatial_plot(
    puck_id: str, gene: str, celltype: str
) -> Path | None:
    if puck_id not in PUCK_IDS:
        raise ValueError(f"unknown puck {puck_id!r}")
    if not gene or not celltype:
        return None

    data = multigene_plot._load_puck_data(puck_id)
    calls = _load_celltype_calls(puck_id)
    x, y = data.x, data.y
    n_beads = data.matrix_csc.shape[0]

    gene_col = data.gene_lower_to_col.get(gene.lower())
    if gene_col is None:
        umi_pos = np.zeros(n_beads, dtype=bool)
        gene_missing = True
    else:
        col_vec = data.matrix_csc.getcol(gene_col)
        umi_pos = np.zeros(n_beads, dtype=bool)
        umi_pos[col_vec.indices[col_vec.data > 0]] = True
        gene_missing = False

    chosen_color = CELL_TYPE_PALETTE.get(celltype, "#000000")
    chosen_mask = (calls == celltype)
    highlight_mask = chosen_mask & umi_pos
    n_highlight = int(highlight_mask.sum())
    n_chosen = int(chosen_mask.sum())

    fig = Figure(figsize=(9, 9), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)

    # Layer 1: doublets / rejects in dim gray, in the back.
    unassigned_mask = (calls == UNASSIGNED_LABEL)
    if unassigned_mask.any():
        ax.scatter(
            x[unassigned_mask], y[unassigned_mask],
            s=3, c=UNASSIGNED_COLOR, alpha=0.25, linewidths=0,
            zorder=1,
        )

    # Layer 2: each known cell type, faded.
    for ct, color in CELL_TYPE_PALETTE.items():
        m = (calls == ct) & ~highlight_mask  # leave highlights for layer 3
        if not m.any():
            continue
        ax.scatter(
            x[m], y[m], s=4, c=color, alpha=0.45, linewidths=0,
            zorder=2,
        )

    # Layer 3: highlighted beads on top.
    if n_highlight > 0:
        ax.scatter(
            x[highlight_mask], y[highlight_mask],
            s=90, c=chosen_color, alpha=1.0,
            edgecolors="black", linewidths=0.8,
            zorder=10,
        )

    ax.set_aspect("equal")
    pad = 50.0
    ax.set_xlim(x.min() - pad, x.max() + pad)
    ax.set_ylim(y.min() - pad, y.max() + pad)
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")
    title = (
        f"{puck_id}: {gene}+ within {celltype}  "
        f"(n={n_highlight} of {n_chosen})"
    )
    if gene_missing:
        title += "  [gene not in matrix]"
    ax.set_title(title, fontsize=11)

    legend_handles = [
        Patch(facecolor=color, edgecolor="none", label=ct)
        for ct, color in CELL_TYPE_PALETTE.items()
    ]
    legend_handles.append(
        Patch(facecolor=UNASSIGNED_COLOR, edgecolor="none",
              label=f"{UNASSIGNED_LABEL} (doublet/reject)")
    )
    legend_handles.append(
        Patch(facecolor=chosen_color, edgecolor="black", linewidth=0.6,
              label=f"highlight: {gene}+ {celltype}")
    )
    # Below the axes, in 2 rows so it never collides with the plot or title.
    ax.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
        fontsize=8,
        framealpha=0.9,
        ncol=5,
        handlelength=1.2,
        columnspacing=1.0,
        borderpad=0.4,
    )

    scalebar = ScaleBar(
        MICRONS_PER_PIXEL, units="um", location="lower right",
        length_fraction=0.15, box_alpha=0.7,
    )
    ax.add_artist(scalebar)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUTPUT_DIR / (
        f"{TODAY}_{puck_id}_{gene}_{_ct_slug(celltype)}_celltype_spatial.png"
    )
    fig.savefig(out, dpi=180)
    return out
