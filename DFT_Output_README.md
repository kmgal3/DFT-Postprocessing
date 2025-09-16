# Documentation for `DFT_Output.py`

This script provides tools for parsing, processing, and plotting data from **DFT (Density Functional Theory)** calculations, focusing on **VASP**, **Wien2k** and **DFTB+** output files.

It is intended for users familiar with the underlying physics and codes.

---

## General Conventions

* Band-/orbital-resolved data arrays are consistently shaped:

  ```
  (nkpts, nbands, 4)
  ```

  where the last axis is `[px, py, pz, z-term]`.

* Large text files (`PROCAR`, `WAVEDER`, etc.) are cached to `.npy` for efficiency.

* Plotting functions may use **hardcoded k-paths** (specific to the author’s system) and should be adapted for new use cases.

---

## Class: `DFT_Output`

### `__init__(self, path)`

Initialize the object with the given calculation directory.

**Parameters**

* `path : str`
  Path to directory containing VASP/Wien2k outputs.

**Notes**

* Initializes empty containers for band structure and momentum data.
* Does not perform file I/O — data must be loaded with `read_*` methods.

---

## File Readers

### `read_WAVEDER(self)`

Reads and parses the `WAVEDER` file.

**Returns**

* `P : ndarray, shape (nkpts, nbands, 4)`
  Parsed momentum/derivative data.

**Notes**

* First call parses raw text, then caches as `WAVEDER.npy`.
* Later calls load directly from `.npy`.

**Errors triggered**

* `FileNotFoundError` if `WAVEDER` is missing.
* `ValueError` if the file contents cannot be reshaped into `(nkpts, nbands, 4)`.

---

### `read_WAVEDERF(self)`

Reads and parses the `WAVEDERF` file.

**Returns**

* `P : ndarray, shape (nkpts, nbands, 4)`

**Notes**

* Structure consistent with `read_WAVEDER`.
* Cached to `WAVEDERF.npy`.

**Errors triggered**

* Same as `read_WAVEDER`.

---

### `read_procar(self)`

Reads the VASP `PROCAR` file containing orbital-resolved band weights.

**Returns**

* `bands : ndarray, shape (nkpts, nbands, norbitals)`

**Notes**

* On first run, parses raw text and saves `PROCAR.npy`.
* On later runs, loads `PROCAR.npy` directly.
* Designed to always use `.npy` for efficiency.

**Errors triggered**

* `FileNotFoundError` if `PROCAR` is missing.
* `ValueError` if parsing encounters mismatched k-point or band counts.

---

### `read_DOSCAR(self)`

Parses the VASP `DOSCAR` file.

**Returns**

* `E : ndarray`
  Energy mesh.
* `DOS : ndarray, shape (nE, norbitals)`
  Orbital-resolved density of states.

**Notes**

* Assumes standard VASP formatting.
* Units are as provided by VASP (not rescaled).

**Errors triggered**

* `FileNotFoundError` if `DOSCAR` is missing.
* `ValueError` if file structure deviates from expectations.

---

### `read_EIGENVAL(self)`

Reads eigenvalues and k-point information from `EIGENVAL`.

**Returns**

* `E : ndarray, shape (nkpts, nbands)`
  Band eigenvalues.
* `kpts : ndarray, shape (nkpts, 3)`
  K-point coordinates.

**Notes**

* Provides energy references for plotting band structures.

**Errors triggered**

* `FileNotFoundError` if `EIGENVAL` is missing.
* `ValueError` if header parsing fails.

---

### `mommat2up(self)`

Reads the Wien2k `mommat2up` file.

**Returns**

* `P : ndarray, shape (nkpts, nbands, 4)`

**Notes**

* Provides momentum matrix elements analogous to `WAVEDER`.

**Errors triggered**

* `FileNotFoundError` if `mommat2up` is missing.

---

## Utility Functions

### `addCHG(file1, file2, outfile)`

Adds two CHG (charge density) files.

**Parameters**

* `file1 : str`
* `file2 : str`
* `outfile : str`

**Notes**

* Assumes consistent grid sizes between `file1` and `file2`.
* Currently supports exactly **two** input CHG files.

**Errors triggered**

* `FileNotFoundError` if either input file is missing.
* `ValueError` if grid sizes differ.

---

## Analysis Methods

### `Optical_Selection_Rules(self, pol, Emin, Emax, nk, nb)`

Computes optical selection rules for transitions under a given polarization.

**Parameters**

* `pol : array-like, shape (3,)`
  Polarization vector.
* `Emin, Emax : float`
  Energy window for transitions.
* `nk : int`
  Number of k-points to include.
* `nb : int`
  Number of bands to include.

**Returns**

* `rules : ndarray`
  Selection rule weights for allowed transitions.

**Notes**

* Works on arrays shaped `(nkpts, nbands, 4)`.

**Errors triggered**

* `ValueError` if `pol` is not length 3.

---

### `L_Plot(self, Emin, Emax, plot)`

Plots angular momentum contributions within a given energy window.

**Parameters**

* `Emin, Emax : float`
* `plot : bool | str`

  * `True` → default plot
  * `'show'` → show interactive plot
  * `'save'` or filename → save figure

**Returns**

* None (creates or saves a plot).

**Errors triggered**

* `ValueError` if `plot` argument is not recognized.

---

### `band_plot(self, Emin, Emax)`

Plots band structure within an energy window.

**Parameters**

* `Emin, Emax : float`

**Returns**

* None (creates plot).

**Notes**

* Uses a **hardcoded k-path** (Γ–X–…) tailored to the author’s system.
* Must be adapted for other use cases.

**Errors triggered**

* `FileNotFoundError` if eigenvalue data not loaded.

---


# DFTB+ Band Structure Methods

### `read_band_out(file="band.out")`

Read band structure data from a **DFTB+ `band.out` file** and store it in the `DFT_Output` object.

#### Parameters

* **file** : `str`, optional
  Path to the `band.out` file. Default is `"band.out"`.

#### Returns

* **eval** : `ndarray, shape (nkpts, nbands)`
  Eigenvalues (band energies in eV) for each k-point.
* **occ** : `ndarray, shape (nkpts, nbands)`
  Band occupancies for each k-point.

#### Other Attributes Updated

* **self.eval** : stores band energies.
* **self.occ** : stores band occupancies.
* **self.kpoints** : list of dicts with metadata for each k-point block, e.g.

```python
[
  {"ik": 1, "spin": 1, "weight": 0.125},
  {"ik": 2, "spin": 1, "weight": 0.250},
  ...
]
```

#### Notes

The expected file format is:

```
KPT i spin j kweight w
 band_index   energy   occupation
 band_index   energy   occupation
 ...
```

* Each block starts with a line beginning with `KPT`.
* Band lines contain three columns: **index, eigenvalue (eV), occupancy**.

#### Example

```python
from DFT_Output import DFT_Output

dftb = DFT_Output("DFTB+")
evals, occs = dftb.read_band_out("band.out")
print(evals.shape)   # (nkpts, nbands)
```

---

### `band_plot_DFTB(ax=None, fermi_level=None, color="black", label="DFTB+")`

Plot a **DFTB+ band structure** loaded with `read_band_out`.

#### Parameters

* **ax** : `matplotlib.axes.Axes`, optional
  Axis object to draw on. If `None`, uses the current axis (`plt.gca()`).

* **fermi\_level** : `float`, optional
  Reference energy (in eV) to shift the bands.
  If `None`, automatically estimated as the **highest occupied energy** (max energy where occupancy > 0.5).

* **color** : `str`, default `"black"`
  Line color for the band structure.

* **label** : `str`, default `"DFTB+"`
  Title label for the plot.

#### Returns

* **ax** : `matplotlib.axes.Axes`
  The axis containing the band structure plot.

#### Notes

* Bands are plotted vs. the **k-point index** (no explicit k-path distances are reconstructed).
* The Fermi level is drawn as a dashed red horizontal line at **0 eV**.

#### Example

```python
import matplotlib.pyplot as plt
from DFT_Output import DFT_Output

dftb = DFT_Output("DFTB+")
dftb.read_band_out("band.out")

fig, ax = plt.subplots(figsize=(6,8))
dftb.band_plot_DFTB(ax=ax, color="blue")
plt.show()
```

---

## Workflow Summary

1. Initialize `DFT_Output(path)`.
2. Load files using `read_*`.
3. Analyze with `Optical_Selection_Rules` / `L_Plot`.
4. Plot results using `band_plot`.

---


# Workflow: VASP vs. DFTB+

The `DFT_Output` class provides a **unified interface** for reading and plotting band structures from different codes:

* **VASP** → uses `EIGENVAL` (with `read_EIGENVAL`)
* **DFTB+** → uses `band.out` (with `read_band_out`)

This makes analysis workflows nearly identical.

---

### 🔹 Example: VASP

```python
import matplotlib.pyplot as plt
from DFT_Output import DFT_Output

# Initialize for VASP
vasp = DFT_Output("VASP", dir="./")

# Read eigenvalues and occupancies
evals, occs = vasp.read_EIGENVAL()

# Plot band structure
fig, ax = plt.subplots(figsize=(6,8))
for ib in range(evals.shape[0]):  # loop over bands
    ax.plot(range(evals.shape[2]), evals[ib,0,:], color="black")

ax.axhline(0, color="red", linestyle="--", label="Fermi level")
ax.set_title("VASP Band Structure")
plt.show()
```

---

### 🔹 Example: DFTB+

```python
import matplotlib.pyplot as plt
from DFT_Output import DFT_Output

# Initialize for DFTB+
dftb = DFT_Output("DFTB+", dir="./")

# Read band.out (energies + occupancies)
evals, occs = dftb.read_band_out("band.out")

# Plot band structure
fig, ax = plt.subplots(figsize=(6,8))
dftb.band_plot_DFTB(ax=ax, color="blue")
plt.show()
```

---

### 🔹 Side-by-side comparison

| Task             | VASP                                | DFTB+                             |
| ---------------- | ----------------------------------- | --------------------------------- |
| Initialize       | `DFT_Output("VASP")`                | `DFT_Output("DFTB+")`             |
| Read bands       | `read_EIGENVAL("EIGENVAL")`         | `read_band_out("band.out")`       |
| Band energies in | `self.eval` (shape `(nbands,1,nk)`) | `self.eval` (shape `(nk,nbands)`) |
| Band occupancies | `self.occ`                          | `self.occ`                        |
| Plot bands       | manual or `band_plot(...)`          | `band_plot_DFTB(...)`             |
| Fermi level      | align to top of valence band        | guessed from occupancies (>0.5)   |

---

With this unified approach, you can use **the same analysis scripts** for both codes by just switching the `source` argument.

---