"""
Multi-gene neighborhood visualization across all 6 melanoma pucks (optimized).

Drop-in replacement for `multigene_neighborhood_plot.py`. Same CLI, same output
filenames (so the on-disk PNG cache is shared with the original), but much
faster when used as a library by the Flask viewer.

Key optimizations vs the original:
  - Per-puck data (sparse counts in CSC, coords, gene-symbol -> column index)
    is loaded once and held in process memory. Subsequent calls reuse it.
  - Gene lookup uses a precomputed dict (O(1)) instead of a linear scan, and
    only a single sparse column is materialised (not the whole dense column).
  - The matplotlib object-oriented API (`Figure`) is used instead of `pyplot`,
    so renders are thread-safe and can run in parallel from Flask workers.
  - Circle layer is drawn with `EllipseCollection` (array-based) rather than
    a `PatchCollection` of N Python `Circle` objects.
  - `prewarm_cache()` is exposed so the Flask app can preload all 6 pucks at
    startup in a background thread.

Usage (CLI is identical to the original):
    python codes/2026-04-25_multigene_neighborhood_plot.py fli1a enpp2 cd8a mitfa --radius-um 50
"""

from __future__ import annotations

import argparse
import gzip
import threading
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.collections import EllipseCollection
from matplotlib_scalebar.scalebar import ScaleBar
import numpy as np
import pandas as pd
import scipy.io

PROJECT_ROOT = Path(__file__).resolve().parent.parent
WORKING_DIR = PROJECT_ROOT / "working"
OUTPUT_DIR = PROJECT_ROOT / "figs"
PUCK_FOLDER_PREFIX = "2024-04-22_"

PUCK_IDS = [
    "Puck_240208_24",
    "Puck_240208_25",
    "Puck_240208_26",
    "Puck_240208_27",
    "Puck_240208_29",
    "Puck_240208_30",
]

# Order matters: 1st gene -> red, 2nd -> green, 3rd -> cyan, 4th -> purple.
GENE_COLORS = ["red", "green", "cyan", "purple"]

UMI_THRESHOLD = 0
RADIUS_UM = 50.0
MICRONS_PER_PIXEL = 0.65
# Match the original output prefix so the Flask cache lookup keeps working
# regardless of which renderer wrote the file. Keep this stable.
TODAY = "2026-04-24"


def puck_dir(puck_id: str) -> Path:
    return WORKING_DIR / f"{PUCK_FOLDER_PREFIX}{puck_id}"


# --------- Per-puck in-memory data cache ---------

class _PuckData:
    __slots__ = ("matrix_csc", "gene_lower_to_col", "x", "y")

    def __init__(self, matrix_csc, gene_lower_to_col, x, y):
        # CSC matrix, rows = beads (already aligned to coords), cols = genes.
        self.matrix_csc = matrix_csc
        self.gene_lower_to_col = gene_lower_to_col
        self.x = x
        self.y = y


_DATA_CACHE: dict[str, _PuckData] = {}
# One lock per puck so two threads asking for the same puck dedupe the load,
# but different pucks load in parallel.
_LOAD_LOCKS: dict[str, threading.Lock] = {p: threading.Lock() for p in PUCK_IDS}


def _load_puck_data(puck_id: str) -> _PuckData:
    cached = _DATA_CACHE.get(puck_id)
    if cached is not None:
        return cached
    lock = _LOAD_LOCKS.setdefault(puck_id, threading.Lock())
    with lock:
        cached = _DATA_CACHE.get(puck_id)
        if cached is not None:
            return cached

        pdir = puck_dir(puck_id)

        # The mtx is stored genes x barcodes; we want CSC over barcodes x genes
        # so getcol(gene_idx) is O(nnz_in_column).
        mtx = scipy.io.mmread(
            pdir / f"{puck_id}.matched.digital_expression_matrix.mtx"
        ).tocsr().T.tocsc()

        barcodes = pd.read_csv(
            pdir / f"{puck_id}.matched.digital_expression_barcodes.tsv",
            sep="\t", header=None
        )[0].tolist()
        features = pd.read_csv(
            pdir / f"{puck_id}.matched.digital_expression_features.tsv",
            sep="\t", header=None
        )
        gene_symbols = features[1].tolist()

        coord_path = pdir / "barcode_matching" / f"{puck_id}_barcode_xy.txt.gz"
        with gzip.open(coord_path, "rt") as fh:
            coords = pd.read_csv(fh, sep="\t", header=None,
                                 names=["barcode", "xcoord", "ycoord"])
        coords = coords.set_index("barcode")

        # Restrict to barcodes that have spatial coordinates (matches original
        # behaviour: counts.index.intersection(coords.index)).
        bc_to_row = {b: i for i, b in enumerate(barcodes)}
        common = [b for b in barcodes if b in coords.index]
        keep_rows = np.fromiter((bc_to_row[b] for b in common),
                                dtype=np.int64, count=len(common))
        mtx_filt = mtx[keep_rows, :]
        coords_filt = coords.loc[common]

        gene_lower_to_col: dict[str, int] = {}
        for i, s in enumerate(gene_symbols):
            if isinstance(s, str):
                # First-seen wins, mirroring `lower.index(target)` in the original.
                gene_lower_to_col.setdefault(s.lower(), i)

        data = _PuckData(
            matrix_csc=mtx_filt,
            gene_lower_to_col=gene_lower_to_col,
            x=coords_filt["xcoord"].to_numpy(),
            y=coords_filt["ycoord"].to_numpy(),
        )
        _DATA_CACHE[puck_id] = data
        return data


def prewarm_cache(puck_ids: list[str] = PUCK_IDS) -> None:
    """Eagerly load each puck's data so the first render after startup is fast."""
    for pid in puck_ids:
        try:
            _load_puck_data(pid)
        except Exception as exc:
            print(f"[prewarm] {pid}: {exc!r}")


def _gene_mask(data: _PuckData, gene: str, threshold: int) -> np.ndarray | None:
    col = data.gene_lower_to_col.get(gene.lower())
    if col is None:
        return None
    n_rows = data.matrix_csc.shape[0]
    col_vec = data.matrix_csc.getcol(col)  # CSC column slice, O(nnz)
    mask = np.zeros(n_rows, dtype=bool)
    nz_rows = col_vec.indices
    nz_vals = col_vec.data
    keep = nz_vals > threshold
    mask[nz_rows[keep]] = True
    return mask


def _add_circles(ax, xs, ys, diameter_data_units, color, alpha):
    if len(xs) == 0:
        return
    offsets = np.column_stack([xs, ys])
    coll = EllipseCollection(
        widths=diameter_data_units,
        heights=diameter_data_units,
        angles=0,
        units="x",  # widths interpreted in data x-units; aspect is equal
        offsets=offsets,
        offset_transform=ax.transData,
        facecolors=color,
        edgecolors=color,
        linewidths=0.6,
        alpha=alpha,
    )
    ax.add_collection(coll)


def make_plot(puck_id: str, genes: list[str], radius_um: float = RADIUS_UM,
              umi_threshold: int = UMI_THRESHOLD) -> Path | None:
    data = _load_puck_data(puck_id)
    n_beads = data.matrix_csc.shape[0]

    masks: list[np.ndarray] = []
    found: list[str] = []
    missing: list[str] = []
    for g in genes:
        m = _gene_mask(data, g, umi_threshold)
        if m is None:
            missing.append(g)
            masks.append(np.zeros(n_beads, dtype=bool))
        else:
            found.append(g)
            masks.append(m)

    if not found:
        print(f"[{puck_id}] none of the input genes found, skipping.")
        return None
    if missing:
        print(f"[{puck_id}] missing genes (drawn as empty): {missing}")

    x = data.x
    y = data.y
    radius_px = radius_um / MICRONS_PER_PIXEL
    diameter_px = 2 * radius_px

    any_mask = np.logical_or.reduce(masks)
    neg_mask = ~any_mask

    fig = Figure(figsize=(9, 9), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)

    ax.scatter(x[neg_mask], y[neg_mask], s=1, c="0.55", alpha=0.35,
               linewidths=0, zorder=1,
               label=f"non-expressing bead (n={int(neg_mask.sum())})")

    for gene, color, mask in zip(genes, GENE_COLORS, masks):
        n = int(mask.sum())
        if n == 0:
            ax.plot([], [], "o", color=color, markersize=6,
                    label=f"{gene} (n=0)")
            continue
        _add_circles(ax, x[mask], y[mask], diameter_px, color=color, alpha=0.15)
        ax.scatter(x[mask], y[mask], s=14, c=color,
                   edgecolors="black", linewidths=0.3, zorder=5,
                   label=f"{gene}+ bead + {radius_um:.0f}µm circle (n={n})")

    ax.set_aspect("equal")
    ax.set_xlim(x.min() - radius_px, x.max() + radius_px)
    ax.set_ylim(y.min() - radius_px, y.max() + radius_px)
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")
    gene_label = ", ".join(
        f"{g} [{c}]" for g, c in zip(genes, GENE_COLORS[:len(genes)])
    )
    ax.set_title(
        f"{puck_id}: expressing beads with {radius_um:.0f} µm proximity circles\n"
        f"{gene_label}; circle overlap = within-{radius_um:.0f}µm co-occurrence"
    )
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)

    scalebar = ScaleBar(MICRONS_PER_PIXEL, units="um", location="lower right",
                        length_fraction=0.15, box_alpha=0.7)
    ax.add_artist(scalebar)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    gene_slug = "_".join(genes)
    out = OUTPUT_DIR / (
        f"{TODAY}_{puck_id}_{gene_slug}_neighborhood_{int(radius_um)}um.png"
    )
    fig.savefig(out, dpi=200)
    print(f"[{puck_id}] wrote {out}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("genes", nargs="+",
                        help="1 to 4 gene symbols (1st=red, 2nd=green, 3rd=cyan, 4th=purple)")
    parser.add_argument("--radius-um", type=float, default=RADIUS_UM)
    parser.add_argument("--threshold", type=int, default=UMI_THRESHOLD,
                        help="UMI count strictly greater than this counts as 'expressing'")
    parser.add_argument("--pucks", nargs="*", default=PUCK_IDS,
                        help="Override which pucks to plot (default: all 6 melanoma pucks)")
    args = parser.parse_args()

    if len(args.genes) > 4:
        parser.error("supply at most 4 genes")

    for puck_id in args.pucks:
        try:
            make_plot(puck_id, args.genes, args.radius_um, args.threshold)
        except Exception as exc:
            print(f"[{puck_id}] FAILED: {exc!r}")


if __name__ == "__main__":
    main()
