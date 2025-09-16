# README: Electronic Structure Analysis Toolkit

## Overview

This toolkit provides Python utilities for **post-processing and visualization of electronic structure data** from **DFT calculations** (VASP, Wien2k, DFTB+) and **Bethe–Salpeter Equation (BSE)** calculations.
It is designed to streamline workflows for extracting data, analyzing orbital/optical properties, and generating **publication-ready plots** of:

* Band structures
* Densities of states (DOS)
* Exciton contributions and fatbands
* Orbital and momentum-resolved properties

The toolkit is organized into three complementary modules:

* **`vaspout_h5.py`** → Process VASP’s `vaspout.h5` (band structure, DOS, MAGMOM extraction).
* **`DFT_Output.py`** → Unified reader/processor for VASP, Wien2k, and DFTB+ outputs.
* **`BSE_tools.py`** → Tools for analyzing VASP BSE calculations (exciton fatbands, orbital composition).

---

## Requirements

* Python **3.8+**
* Core dependencies:

  * `numpy`
  * `matplotlib`
  * `h5py`
* Local modules:

  * `Plotting` (band structure & DOS plotting classes)
  * `InformationCollector`

---

## Module Summaries

### 🔹 1. `vaspout_h5.py`

Utility for parsing VASP’s **HDF5 output** (`vaspout.h5`) and producing customizable plots.

**Features**

* Extract MAGMOM data (VASP-compatible `INCAR` format).
* Parse and plot band structures (BS) and densities of states (DOS).
* Generate `.png` plots or gnuplot-compatible `.txt` exports.
* Python API + CLI interface.

**Example CLI**

```bash
# Band structure projected on Mo dxy orbitals
python vaspout_h5.py band "1,2 Mo dxy" --bands "5-15" --E0 fermi --save

# Total DOS
python vaspout_h5.py dos "total" --E0 0 --save --gnuplot
```

**Example Python**

```python
from vaspout_h5 import vaspout_h5

calc = vaspout_h5("vaspout.h5")
calc.plot_BS("1,2 Mo dxy", bands="5-15", E0="fermi")
```

---

### 🔹 2. `DFT_Output.py`

General-purpose interface for **DFT outputs** from VASP, Wien2k, and DFTB+.

**Features**

* Readers for `WAVEDER`, `PROCAR`, `DOSCAR`, `EIGENVAL`, Wien2k `mommat2up`, and DFTB+ `band.out`.
* Caches large files (`.npy`) for efficiency.
* Provides band/orbital-resolved data in consistent array shapes.
* Analysis methods:

  * Optical selection rules
  * Angular momentum projections
  * Band structure plotting
  * CHG file combination

**Example Workflow (VASP)**

```python
from DFT_Output import DFT_Output

vasp = DFT_Output("VASP")
evals, kpts = vasp.read_EIGENVAL()
vasp.band_plot(Emin=-5, Emax=5)
```

**Example Workflow (DFTB+)**

```python
from DFT_Output import DFT_Output

dftb = DFT_Output("DFTB+")
evals, occs = dftb.read_band_out("band.out")
dftb.band_plot_DFTB(color="blue")
```

---

### 🔹 3. `BSE_tools.py`

Specialized tools for analyzing **Bethe–Salpeter Equation (BSE)** outputs (VASP `BSEFATBAND` + `PROCAR`).

**Features**

* Read and cache `PROCAR` projections.
* Process `BSEFATBAND` exciton contributions.
* Generate exciton fatbands (single bands or grids).
* Compute orbital composition of excitons (system-specific, Cr–S–Br style).

**Example Workflow**

```python
from BSE_tools import ReadPROCAR, ExcitonFatband

promat, energs, *_ = ReadPROCAR("PROCAR")
ExcitonFatband([1], i_band=5, i_kz=2, folder="./BSE/")
```

**Visualization**

* `ExcitonFatband`: exciton weights for a given band/kz slice.
* `FatbandGrid`: multi-band/kz exciton grids.
* `print_OrbComp`: orbital composition reports.

---

## Summary

1. **Electronic structure (VASP)** → use `vaspout_h5.py` or `DFT_Output.py`.
2. **Compare across codes (VASP, Wien2k, DFTB+)** → use `DFT_Output.py`.
3. **Exciton analysis (BSE)** → use `BSE_tools.py`.
4. **Visualization** → all modules provide plotting routines with `matplotlib`.

---

## Error Handling

* Missing files → `FileNotFoundError`.
* Unexpected format → `ValueError`.
* Array reshaping issues → `ValueError`.

---

## License & Contribution

This toolkit is intended for research use. Please adapt plotting routines (e.g., k-paths, orbital definitions) to your system as needed.
Contributions (extensions, generalizations, bug fixes) are welcome.

---
