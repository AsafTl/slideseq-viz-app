# Slideseq Viz App

A small local web app for exploring **Slide-seq** spatial transcriptomics data
from 6 zebrafish melanoma pucks. The app has two tabs:

- **Neighborhood** — type up to 4 gene names and the app draws, for each puck,
  a map of where those genes are expressed. Each gene gets its own color, and
  a circle of adjustable radius is drawn around every expressing bead so you
  can see where genes co-occur in space.
- **Cell-type expression** — uses the RCTD per-bead cell-type calls to show
  either (a) a **dotplot** of how strongly your chosen genes are expressed in
  each cell type — one per puck, plus a cumulative one across a chosen subset
  of pucks — or (b) a **spatial highlight map** that colors the puck by cell
  type and highlights the beads of one chosen cell type that express one
  chosen gene.

The app runs **entirely on your own computer**. Nothing is uploaded anywhere.
The browser part is just the user interface; the actual rendering is done
locally by Python.

---

## What you'll need

1. **A copy of this code** (cloned from GitHub — instructions below).
2. **Python 3.10 or newer** installed on your computer.
3. **The Slide-seq data files** (not in this repo).
   Total size: roughly **700 MB** for all 6 pucks.
4. About 15 minutes for first-time setup. After that, starting the app
   takes 5 seconds.

You do **not** need to know any Python to use the app — but you'll be
typing a few commands into a terminal during setup.

---

## Step 1 — Install Python

Skip this step if you already have Python 3.10 or newer installed. To check,
open a terminal (see Step 2 for how) and type:

```
python --version
```

If you see `Python 3.10.x` or higher, you're good. If you see `Python 2.x`,
or a "command not found" error, install Python as below.

### Windows

1. Go to <https://www.python.org/downloads/windows/> and download the
   latest **Windows installer (64-bit)**.
2. **Important:** on the very first screen of the installer, tick the box
   that says **"Add python.exe to PATH"** before clicking Install. If you
   miss this, none of the commands below will work.
3. Click "Install Now" and wait.

### Mac

The version of Python that ships with macOS is too old. Install a fresh one:

- **Easiest option:** download the installer from
  <https://www.python.org/downloads/macos/> and run it.
- **Or, if you already use Homebrew:** `brew install python@3.12`

---

## Step 2 — Open a terminal

You'll run a few commands in a terminal.

- **Windows:** press the Start key, type `PowerShell`, and open
  **Windows PowerShell**. (Command Prompt also works, but the commands
  below assume PowerShell.)
- **Mac:** press `Cmd + Space`, type `Terminal`, and press Enter.

Keep this window open for the rest of setup.

---

## Step 3 — Clone the repository

"Cloning" just means downloading the code from GitHub. Pick **one** of the
two methods below.

### Method A — using Git (recommended if you already have Git installed)

```
git clone https://github.com/<asaf-or-org>/slideseq-viz-app.git
cd slideseq-viz-app
```

### Method B — download the ZIP (no Git needed)

1. Go to the repository page on GitHub.
2. Click the green **"Code"** button → **"Download ZIP"**.
3. Unzip the file somewhere convenient — e.g.
   - Windows: `C:\Users\<your-name>\Documents\slideseq-viz-app`
   - Mac: `~/Documents/slideseq-viz-app`
4. In your terminal, change into that folder. For example:
   - Windows: `cd C:\Users\<your-name>\Documents\slideseq-viz-app`
   - Mac: `cd ~/Documents/slideseq-viz-app`

From here on, every command in this README assumes your terminal is sitting
inside this project folder.

---

## Step 4 — Install the Python libraries

The app needs a handful of Python libraries (Flask for the web server,
matplotlib for plotting, etc.). They are all listed in
[`requirements.txt`](requirements.txt).

The cleanest way to install them is in a **virtual environment** ("venv") —
a private Python sandbox tied to this project, so it cannot interfere with
any other Python you might use.

### Windows (PowerShell)

```
python -m venv venv
venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

> **If PowerShell complains about "scripts being disabled":** run this once,
> then try the activate command again:
> ```
> Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
> ```

### Mac (Terminal)

```
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

After activation you should see `(venv)` at the start of your prompt — that
means the sandbox is active. The `pip install` step downloads several
hundred MB and may take a few minutes the first time.

> Whenever you open a fresh terminal later, **re-activate the venv first**
> (`venv\Scripts\Activate.ps1` on Windows, `source venv/bin/activate` on
> Mac) before running the app.

---

## Step 5 — Add the data files

The repository on GitHub contains **only the code**, not the actual
Slide-seq data.

You need 5 files per puck, for each of the 6 pucks:
`Puck_240208_24`, `Puck_240208_25`, `Puck_240208_26`, `Puck_240208_27`,
`Puck_240208_29`, `Puck_240208_30`.

The folder structure must look **exactly** like this — the app builds file
paths from these names, and a typo will cause it to fail to find the data:

```
slideseq-viz-app/
├── codes/
│   ├── slideseq_viz_app.py
│   ├── slideseq_viz.html
│   ├── multigene_neighborhood_plot.py
│   └── celltype_expression_plot.py
├── working/
│   ├── 2024-04-22_Puck_240208_24/
│   │   ├── Puck_240208_24.matched.digital_expression_matrix.mtx
│   │   ├── Puck_240208_24.matched.digital_expression_barcodes.tsv
│   │   ├── Puck_240208_24.matched.digital_expression_features.tsv
│   │   ├── 2024-04-22_Puck_240208_24_results.csv
│   │   └── barcode_matching/
│   │       └── Puck_240208_24_barcode_xy.txt.gz
│   ├── 2024-04-22_Puck_240208_25/
│   │   └── ... (same 5 files, with "25" instead of "24")
│   ├── 2024-04-22_Puck_240208_26/
│   │   └── ... (same 5 files, with "26" instead of "24")
│   ├── 2024-04-22_Puck_240208_27/
│   │   └── ... (same 5 files, with "27" instead of "24")
│   ├── 2024-04-22_Puck_240208_29/
│   │   └── ... (same 5 files, with "29" instead of "24")
│   └── 2024-04-22_Puck_240208_30/
│       └── ... (same 5 files, with "30" instead of "24")
│
├── requirements.txt
└── README.md
```

Notes:
- The `working/` folder and the 6 empty puck sub-folders are **already
  present in the repo** — you only need to drop the data files into the
  matching folders. (You'll see small `.gitkeep` placeholder files inside
  each empty folder; ignore them — they only exist so git keeps the
  folders around.)
- Each puck folder takes 4 files at the top level (`...mtx`,
  `...barcodes.tsv`, `...features.tsv`, `..._results.csv`) and 1 file
  inside the `barcode_matching/` sub-folder (`..._barcode_xy.txt.gz`).
- The `..._results.csv` is the per-bead RCTD output (cell-type assignment
  per bead). It's only used by the **Cell-type expression** tab; the
  Neighborhood tab works without it. If you don't have it, that tab will
  load empty — drop the file in and restart the app.
- It's fine if the puck folders contain extra files — the app only reads
  the 5 listed above and ignores everything else.


---

## Step 6 — Run the app

With the venv still active and your terminal sitting inside the project
folder:

```
python codes/slideseq_viz_app.py
```

You should see something like:

```
 * Running on http://127.0.0.1:5000
```

Open that URL in your browser (Chrome, Safari, Firefox, etc.). The app
loads with a search box for genes and a puck selector.

To stop the app, go back to the terminal and press `Ctrl + C`.

### How to use it

The app has two tabs at the top: **Neighborhood** and **Cell-type expression**.
Both render all 6 pucks at once in a 2×3 grid; click any tile to expand it
and download the PNG.

#### Neighborhood

1. Type 1 to 4 gene symbols in the box, separated by commas
   (e.g. `cd8a, mitfa, fli1a`). The first gene is shown in red, the
   second in green, the third in cyan, the fourth in purple.
2. Optionally change the proximity radius (default 50 µm).
3. Click **Visualize**. Each tile shows where those genes are expressed
   on its puck, with a circle around every expressing bead.

#### Cell-type expression

1. Type 1 to 4 gene symbols.
2. Pick a cell type from the dropdown (only used by the Spatial sub-mode).
3. Choose a sub-mode using the toggle: **Dotplot** or **Spatial**.
4. (Dotplot only) The checkboxes above the cumulative panel control which
   pucks are pooled into the cumulative dotplot. The 6 per-puck dotplots
   in the grid below are unaffected.
5. Click **Visualize**.
   - In **Dotplot** mode, a full-width "cumulative" dotplot appears above
     the per-puck grid. Cell types are on the Y axis, your genes on the X.
     Dot size = % of beads of that cell type with UMI > 0; dot color =
     mean UMI among the expressing beads (so a true marker shows up as a
     **large** *and* **bright** dot).
   - In **Spatial** mode, each tile shows the puck colored by RCTD cell
     type, with the beads of your chosen cell type that also express the
     first gene in your list highlighted on top.

---

## Troubleshooting

| Symptom | Likely fix |
|---|---|
| `python: command not found` (Mac) | Try `python3` instead of `python`. |
| `'python' is not recognized` (Windows) | Python wasn't added to PATH during install. Re-run the installer and tick "Add python.exe to PATH". |
| `ModuleNotFoundError: No module named 'flask'` | The venv isn't active, or `pip install -r requirements.txt` hasn't been run. Re-activate the venv (Step 4) and re-run the install. |
| App starts but says "no such file" when rendering | The data folder layout doesn't match Step 5 exactly — check the puck folder name (`2024-04-22_Puck_...`) and that `barcode_xy.txt.gz` is inside `barcode_matching/`. |
| `Address already in use` / port 5000 busy | Another program is using port 5000 (on Mac, AirPlay Receiver does this by default). Either turn that off, or edit the last line of `codes/slideseq_viz_app.py` and change `port=5000` to `port=5001`. |
| Browser shows the page but rendering hangs | Check the terminal — Python prints progress and any errors there. |


