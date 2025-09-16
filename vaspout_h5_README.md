# README: `vaspout_h5.py`

## Overview

`vaspout_h5.py` is a Python utility for **post-processing electronic structure data** from VASP’s `vaspout.h5` output files.
It provides a simple way to extract data and generate **publication-ready plots** of band structures (BS) and densities of states (DOS), directly from the command line or from within Python.

Key features:

* Extract **MAGMOM** data into VASP-compatible `INCAR` format.
* Parse **band structure** and **density of states** from `vaspout.h5`.
* Generate **customizable BS and DOS plots** as `.png` images.
* Export plot data in **gnuplot-compatible `.txt` files**.
* Includes a command-line interface (CLI) for quick post-processing.

---

## Requirements

* Python 3.8+
* Dependencies:

  * `h5py`
  * `numpy`
  * `matplotlib`
* Local modules:

  * `Plotting` (provides `BandStructurePlot`, `DensityOfStatesPlot`)
  * `InformationCollector`

⚠️ Users are expected to interact only with the `vaspout_h5` module. The other modules must be available in the project path.

---

## CLI Usage

### Band structure plot

```bash
python vaspout_h5.py band "1,2 Mo x p dxy" --kpaths "1 2 3 -4" --bands "3-10" --E0 fermi --save
```

### Density of states plot

```bash
python vaspout_h5.py dos "total" --E0 0 --save
```

### Arguments

* **Required**

  * `calc_type`: `band` or `dos`
  * `description`: Projection string (see below).
* **Optional**

  * `--kpaths`: Select ordering and orientation of k-point paths.
    Example: `"1 2 3 -4"` → use segments 1→2→3, then reverse segment 4.
  * `--bands`: Select band indices or ranges.
    Examples: `"5"` (only band 5), `"3-10"` (bands 3 through 10).
  * `--E0`: Reference energy (any `float` or `"fermi"`). By default it is set to the valence band maximum
  * `--gnuplot`: Export plot data to `.txt` file for gnuplot.
  * `--save`: Save `.png` plot instead of showing interactively.
  * `--bandnums`: Display band indices on BS plots.
  * `--color`, `--mult`, `--alpha`: Control marker color, size, and opacity.
  * `--folder`: Output directory (default: same as `vaspout.h5`).

---

## Syntax for `description`, `kpaths`, and `bands`

### Projection description (`description`)

The `description` string specifies **which atomic orbitals** to sum over when plotting.

* Format:

  ```
  "atom_indices atom_species direction orbitals"
  ```

* Example:

  ```
  "1,2 Mo x p dxy"
  ```

  means:

  * Atoms 1 and 2 of Molybdenum
  * Projection along the `x` direction
  * Orbitals: `p` and `dxy`

* Notes:

  * Atom indices are **1-based**.
  * You may combine multiple indices using commas: `"1,2,5"`.
  * Summation over **different atomic species is not supported** (unless the atomic species is not specified, in which case the summation is over every atom).
    E.g., `"Mo Te"` is **invalid**.
  * `"total"` can be used to plot the total DOS.
  * `"Mo"` will sum over all Molybdenum atoms, over all directions and over all orbitals
  * In spin polarized calculations if direction is not specified, the `"up"` direction is chosen

---

### K-path specification (`--kpaths`)

For band structure plots, the k-path order and direction can be customized.

* Example:

  ```
  --kpaths "2 3 -4 1"
  ```

  * Traverse segment 2 → 3 normally
  * Traverse segment 4 in reverse (negative sign)
  * Traverse segment 1

If omitted, the default k-path order from KPOINTS is used.

---

### Band selection (`--bands`)

Select which bands to display.

* Examples:

  * `"5"` → Only band 5.
  * `"3-10"` → Bands 3 through 10.
  * `"3,7,12"` → Specific bands 3, 7, and 12.
  * `"3,12-18,22"` → Specific bands 3, 22, and 12 through 18.

---

## Output

* **Plots**: Saved as `.png` images (default: `band_plot.png` or `dos_plot.png`).
* **Gnuplot data**: If `--gnuplot` is set, a `.txt` file is generated with the same base filename.
* Users can override output format and naming in custom scripts.

---

## Example Workflows

### 1. Extract MAGMOM

```python
from vaspout_h5 import vaspout_h5

calc = vaspout_h5("vaspout.h5")
calc.extract_magmom()
```

Generates a `MAGMOM` file compatible with `INCAR`.

---

### 2.1 Band structure plot (CLI)

```bash
python vaspout_h5.py band "1,2 Mo dxy" --bands "5-15" --E0 fermi
```

* Plots bands 5–15 projected onto the `dxy` orbital of atoms 1 and 2 of Mo.
* Fermi level set as energy reference.
* Saves as `band_plot.png`.

---

### 2.2 Band structure plot (Python API)

```python
from vaspout_h5 import vaspout_h5

plt.figure(figsize=(4,8))

calc = vaspout_h5("vaspout.h5")
calc.plot_BS("1,2 Mo dxy", bands="5-15", E0=fermi)

plt.tight_layout()
plt.savefig(f'{calc.ic.folder}band_plot')
```

* Gives the same output as 2.1

---

### 3. DOS plot (CLI)

```bash
python vaspout_h5.py dos "total" --E0 0 --save --gnuplot
```

* Plots the total DOS with zero as energy reference.
* Saves plot as `dos_plot.png` and gnuplot data as `dos_plot.txt`.

---

### 4. Combined BS + DOS (Python API)

```python
from vaspout_h5 import vaspout_h5, plot_BS_DOS

BS = vaspout_h5("vaspout.h5")
DOS = vaspout_h5("vaspout.h5")

plot_BS_DOS(BS, DOS, description="1,2 Mo dxy")
```

Produces a side-by-side BS and DOS plot, aligned to the same energy scale.