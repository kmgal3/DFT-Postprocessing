# Documentation for `BSE_tools.py`

This script provides tools for processing and analyzing results of **Bethe–Salpeter Equation (BSE)** calculations, in particular VASP’s `BSEFATBAND` outputs and associated `PROCAR` projections.

It includes functionality for:

* Reading and caching `PROCAR` with orbital decompositions.
* Extracting and processing excitonic contributions from `BSEFATBAND`.
* Constructing exciton “fatband” plots across bands and k-points.
* Computing system-specific orbital contributions to excitons.

The script is tailored for **CrSBr-like systems** (with Cr–d, S–p, Br–p orbital decomposition) but can be adapted.

---

## General Conventions

* Large files are cached as `.npy` on first read for efficiency.
* Exciton indices (`l_exc`) are **integers** referencing excitons in `BSEFATBAND`.
* Spin channels are handled via:

  * `'up'` = spin-up only
  * `'down'` = spin-down only
  * `'sum'` = AFM-like (sum of both spin channels)
* Orbital contributions in `print_OrbComp` are **system-specific** (Cr, S, Br).

---

## Function Reference

### `ReadPROCAR(PROCAR)`

Reads a `PROCAR` file and caches results.

**Parameters**

* `PROCAR : str`
  Path to `PROCAR` file.

**Returns**

* `promat : ndarray, shape (nkpts, nbands, nprjct, nions, norbitals)`
* `energs : ndarray, shape (nkpts, nbands, nprjct)`
* `E_max : float` – max valence band energy
* `E_min : float` – min conduction band energy
* `E_gap : float` – band gap
* `kindexer : ndarray, shape (nkpts, 3)` – list of k-point vectors
* `weights : ndarray` – k-point weights (if available, otherwise empty placeholder)

**Notes**

* If `PROCAR.npy` exists, loads cached arrays instead of parsing text.
* Determines number of orbitals and projections automatically.

**Errors triggered**

* `FileNotFoundError` if `PROCAR` is missing.
* `ValueError` if parsing fails or reshaping arrays is inconsistent.

---

### `ProcessBSEFATBAND(l_exc, i_band, i_kpts=[None,None,None], c_kpts=[None,None,None], spins='sum', folder='./')`

Processes `BSEFATBAND` to extract exciton contributions for a given band.

**Parameters**

* `l_exc : list[int]` – exciton indices
* `i_band : int` – band index
* `i_kpts : list[int|None]` – indices of k-points along (kx, ky, kz); can restrict selection
* `c_kpts : list[float|None]` – coordinates of k-points (overrides indices if given)
* `spins : {'up','down','sum'}` – spin channel handling
* `folder : str` – path to `BSEFATBAND` and outputs

**Returns**

* `totmat : ndarray, shape (nkpts, 4)` – kx, ky, kz, weight

**Notes**

* Attempts to load small preprocessed `BSE_*` files.
* If not found, generates them by creating and running a `processing.sh` script.

**Errors triggered**

* `FileNotFoundError` if neither small files nor `BSEFATBAND` are available.

---

### `MakeBSEband(lexc, bands, spins, start_index=1, BSEfolder='./', PROCAR='./PROCAR')`

Writes band-resolved exciton weights into a file.

**Parameters**

* `lexc : list[int]` – exciton indices
* `bands : tuple(int,int)` – (min\_band, max\_band) range
* `spins : {'up','down','sum'}`
* `start_index : int` – starting k-point index (skip earlier ones)
* `BSEfolder : str` – path to BSE files
* `PROCAR : str` – path to PROCAR file

**Returns**

* None (writes `bseband_*` file to `BSEfolder`)

**Notes**

* Matches exciton weights with band energies from `PROCAR`.
* Output file contains columns: k\_idx, k\_x, k\_y, k\_z, energy, Abs(X)/W\_k.

---

### `ExcitonFatband(l_exc, i_band, i_kz, ax=None, folder='./')`

Generates exciton “fatband” data for a given kz slice.

**Parameters**

* `l_exc : list[int]` – exciton indices
* `i_band : int` – band index
* `i_kz : int` – kz index
* `ax : matplotlib.axes.Axes | None` – optional axes to plot on
* `folder : str`

**Returns**

* `totmat : ndarray, shape (nkpts, 3)` – kx, ky, weight

**Notes**

* If small BSE files are not present, generates them from `BSEFATBAND` using `processing.sh`.
* Saves plot as `BSE_<exc_list>_<band>_<kz>.png` if plotting is enabled.

---

### `FatbandGrid(l_exc, l_band, n_kz, quart=1, folder='./', plot_folder='./', figsize=(10,10))`

Creates a grid of exciton fatband plots across multiple bands and kz slices.

**Parameters**

* `l_exc : list[int]`
* `l_band : list[int]`
* `n_kz : int` – number of kz points
* `quart : {1,2,3,4,'all'}` – restricts plot quadrant
* `folder : str` – path to data
* `plot_folder : str` – output folder for plots
* `figsize : tuple(float,float)` – figure size

**Returns**

* None (saves grid plot to file)

---

### `print_OrbComp(lexc, iband, spins, nkz=3, folder='./', promat=[], kindexer={}, weights=[])`

Computes orbital composition of excitons for a given band.

**Parameters**

* `lexc : list[int]` – exciton indices
* `iband : int` – band index
* `spins : {'up','down','sum'}`
* `nkz : int` – number of kz slices
* `folder : str`
* `promat, kindexer, weights : ndarray` – optionally provide preloaded PROCAR data

**Returns**

* None (writes file `OrbComp_<exc>_<band>_<spin>` in folder)

**Notes**

* System-specific: assumes orbitals for Cr, S, Br with names hardcoded in the function.
* Contributions are normalized and reported as percentages.
* Output format: element + orbital + percentage.

---

## Error Triggers

* **Missing files**: All readers (`PROCAR`, `BSEFATBAND`, small BSE files).
* **Unexpected formatting**: Misaligned lines in `PROCAR` or `BSEFATBAND`.
* **Reshape failures**: If data arrays cannot be matched to expected dimensions.

---

## Workflow Summary

1. Run `ReadPROCAR` to prepare orbital-projected data.
2. Use `ProcessBSEFATBAND` or `ExcitonFatband` to extract excitonic weights.
3. Combine band structure with exciton weights using `MakeBSEband`.
4. Visualize results:

   * Single exciton/band → `ExcitonFatband`
   * Multiple → `FatbandGrid`
5. Compute orbital contributions with `print_OrbComp`.

---