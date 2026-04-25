"""
Local Flask backend for the Slideseq gene-neighborhood viewer (optimized).

Drop-in replacement for `slideseq_viz_app.py`. Same HTTP API, same HTML, same
on-disk PNG cache layout in `figs/`. The performance gains come from the new
renderer (`2026-04-25_multigene_neighborhood_plot.py`) plus a few server-side
changes:

  - Per-puck data is loaded once into process memory (renderer-side cache),
    so repeat renders for new gene sets don't re-read the .mtx files.
  - The global `_render_lock` is gone. The new renderer uses the matplotlib
    OO API (`Figure(...)` with no `pyplot`), which is thread-safe, so the 6
    parallel puck requests from the frontend now actually render in parallel.
  - Per-(puck, genes, radius) locks still dedupe identical concurrent requests.
  - All 6 pucks are pre-warmed in a background thread on startup, so the first
    user-visible render no longer pays the cold data-load cost.

Run:
    pip install flask
    python codes/2026-04-25_slideseq_viz_app.py
    # then open http://localhost:5000
"""

from __future__ import annotations

import importlib.util
import threading
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

from flask import Flask, jsonify, request, send_file, abort

APP_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = APP_DIR.parent
FIGS_DIR = PROJECT_ROOT / "figs"
HTML_PATH = APP_DIR / "slideseq_viz.html"
RENDERER_PATH = APP_DIR / "2026-04-25_multigene_neighborhood_plot.py"

# Filename starts with a date, so import via importlib rather than `import ...`.
_spec = importlib.util.spec_from_file_location("multigene_plot_v2", RENDERER_PATH)
multigene_plot = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(multigene_plot)

PUCK_IDS = multigene_plot.PUCK_IDS
RADIUS_UM = multigene_plot.RADIUS_UM
GENE_COLOR_NAMES = multigene_plot.GENE_COLORS

app = Flask(__name__)

# Per-(puck, genes, radius) lock so two clients asking for the exact same
# render dedupe instead of rendering twice. Different requests run in parallel.
_render_locks: dict[tuple, threading.Lock] = {}
_render_locks_guard = threading.Lock()


def _get_key_lock(key: tuple) -> threading.Lock:
    with _render_locks_guard:
        lk = _render_locks.get(key)
        if lk is None:
            lk = threading.Lock()
            _render_locks[key] = lk
        return lk


def _parse_genes(raw: str) -> list[str]:
    return [g.strip() for g in raw.split(",") if g.strip()][:4]


def _parse_radius(raw: str | None) -> float:
    if not raw:
        return RADIUS_UM
    try:
        v = float(raw)
    except ValueError:
        return RADIUS_UM
    return max(1.0, min(500.0, v))


def _cache_path(puck_id: str, genes: list[str], radius_um: float) -> Path | None:
    """Most recent cached PNG for this (puck, gene-list, radius), or None.

    Wildcard prefix matches PNGs written by either the original renderer
    (date prefix `2026-04-24_`) or the optimized one — they share the cache.
    """
    gene_slug = "_".join(genes)
    pattern = f"*_{puck_id}_{gene_slug}_neighborhood_{int(radius_um)}um.png"
    matches = sorted(FIGS_DIR.glob(pattern),
                     key=lambda p: p.stat().st_mtime, reverse=True)
    return matches[0] if matches else None


@app.route("/")
def index():
    return HTML_PATH.read_text(encoding="utf-8")


@app.route("/api/config")
def config():
    return jsonify({
        "puck_ids": PUCK_IDS,
        "gene_color_names": GENE_COLOR_NAMES,
        "radius_um": RADIUS_UM,
    })


@app.route("/api/render")
def render():
    genes = _parse_genes(request.args.get("genes", ""))
    puck = request.args.get("puck", "")
    radius_um = _parse_radius(request.args.get("radius"))
    if not genes:
        return jsonify({"error": "no genes provided"}), 400
    if puck not in PUCK_IDS:
        return jsonify({"error": f"unknown puck {puck!r}"}), 400

    cached = _cache_path(puck, genes, radius_um)
    if cached is None:
        key = (puck, tuple(genes), radius_um)
        with _get_key_lock(key):
            cached = _cache_path(puck, genes, radius_um)
            if cached is None:
                try:
                    cached = multigene_plot.make_plot(puck, genes, radius_um)
                except Exception as exc:
                    return jsonify({"error": repr(exc)}), 500
                if cached is None:
                    return jsonify({"error": "renderer returned no file"}), 500

    return jsonify({"image_url": f"/api/image/{cached.name}"})


def _safe_fig(filename: str) -> Path:
    if "/" in filename or "\\" in filename or ".." in filename:
        abort(400)
    path = FIGS_DIR / filename
    if not path.exists():
        abort(404)
    return path


@app.route("/api/image/<path:filename>")
def image(filename: str):
    return send_file(_safe_fig(filename), mimetype="image/png")


@app.route("/api/download/<path:filename>")
def download(filename: str):
    return send_file(_safe_fig(filename), mimetype="image/png",
                     as_attachment=True, download_name=filename)


if __name__ == "__main__":
    FIGS_DIR.mkdir(parents=True, exist_ok=True)
    # Load all 6 pucks into the renderer's in-memory cache before the first
    # user request. Done in a daemon thread so the server is responsive
    # immediately; per-puck load locks ensure a request that races the
    # pre-warm just blocks until that puck is ready.
    threading.Thread(target=multigene_plot.prewarm_cache, daemon=True).start()
    app.run(host="127.0.0.1", port=5000, debug=False, threaded=True)
