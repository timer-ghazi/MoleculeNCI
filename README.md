# Molecule with NCIs Python Library

This repository provides two Python modules for representing molecules, detecting **covalent** and **non-covalent** interactions (NCIs), and calculating molecular geometry (distances, angles, dihedrals, etc.). 
---

## Table of Contents

1. [Repository Overview](#repository-overview)
2. [Getting Started](#getting-started)
3. [Data Structures and APIs](#data-structures-and-apis)
4. [Usage](#usage)
5. [Testing](#testing)
    - [Test Files and Organization](#test-files-and-organization)
    - [Running the Tests](#running-the-tests)
    - [What the Tests Check](#what-the-tests-check)
6. [Extending the Code](#extending-the-code)
7. [License](#license)
8. [Contact](#contact)

---

## Repository Overview

Your folder structure might look like this:

```
MoleculeNCI/
├─ molecule.py
├─ molecule_nci.py
├─ elements_table.py
├─ tests/
│   ├─ test_molecule_nci.py
│   └─ data/
│       ├─ CBr4_urea.xyz
│       ├─ CF3Br_H2O.xyz
│       ├─ ...
└─ README.md
```

**Key files**:

- **`molecule.py`**  
  Defines the base `Molecule` class with:
  - Methods to parse `.xyz` files  
  - Covalent bond detection (`detect_bonds`)  
  - Fragment detection  
  - Geometry utilities (distance, angle, dihedral)  
  - A `summary()` method for a human-readable overview  

- **`molecule_nci.py`**  
  Defines `MoleculeWithNCIs`, a subclass of `Molecule` that:
  - Detects various non-covalent interactions (NCIs) (hydrogen bonds, halogen/chalcogen bonds, steric clashes, etc.)  
  - Uses a registration system (`@register_nci`) for flexible addition of new NCI detection methods  
  - Offers methods to print or list NCIs, including geometry and whether they are intramolecular or intermolecular  

- **`elements_table.py`**  
  A helper module that provides chemical element data—covalent radii, van der Waals radii, etc. **We keep a local copy in this repository for convenience**, but the file is originally maintained in the [Elements-Table repo](https://github.com/timer-ghazi/Elements-Table). If you want to stay fully up to date with new data or bug fixes, you can download the latest version from there.

- **`tests/test_molecule_nci.py`**  
  A **pytest**-based test suite, referencing **`.xyz` files** under `tests/data/`. Verifies that the library detects the correct number of fragments and NCIs, and performs **fine-grained** geometry checks (distance/angle).

---

## Getting Started

1. **Clone or Download** this repository:

    ```bash
    git clone https://github.com/yourusername/NewMoleculeNCIs.git
    cd NewMoleculeNCIs
    ```

2. **Install Dependencies** (e.g., NumPy, pytest, etc.):

    ```bash
    pip install numpy pytest
    ```

3. **Directory Layout**  
   Make sure your `.xyz` files for testing are in `tests/data/`.

You can then run the scripts in place or proceed to testing (described below).

---

## Data Structures and APIs

### `Atom` Dataclass

Defined in `molecule.py`:
```python
@dataclass
class Atom:
    symbol: str
    x: float
    y: float
    z: float
    charge: float = 0.0
```

- A straightforward container for atomic data.  
- Coordinates are in Å by convention.

### `Molecule` Base Class

Core attributes:

- **`self.title`**: A string describing the molecule.  
- **`self.atoms`**: A list of `Atom` objects.  
- **`self.bond_matrix`**: A NumPy array (`shape = (n, n)`) for covalent bonding.  
- **`self.fragments`**: A dict mapping `fragment_id -> [atom_indices]` once `find_fragments()` is called.

Important Methods:

- **`from_xyz(cls, filename)`**: Class method to parse an XYZ file.  
- **`detect_bonds(tolerance=0.3)`**: Populate `bond_matrix` by comparing distances to sum of covalent radii.  
- **`find_fragments()`**: Identify connected components in `bond_matrix`.  
- **`distance(i, j)`**, **`angle(i, j, k)`**, **`dihedral(i, j, k, l)`**: Geometry calculations.  
- **`summary()`**: Return a string summarizing the molecule.

### `MoleculeWithNCIs` Subclass

Extends `Molecule` with attributes and methods for **non-covalent interactions** (NCIs):

- **`self.ncis`**: A dict mapping `(i, j) -> list of interaction dicts`. Each dict stores metadata (`type`, `distance`, `angle`, `intra_or_inter`, etc.).  
- **`detect_all_ncis(run_steric_clashes=False, debug=False)`**: Calls each registered NCI detection method in ascending order of priority.  
- **`@register_nci(priority=0)`**: Decorator to add custom NCI detection functions.  
- Built-in NCI detection methods:
  - **`detect_hydrogen_bonds(...)`**  
  - **`detect_sigma_hole_bonds(...)`**  
  - **`detect_steric_clashes(...)`**  
- **`list_interactions(...)`**: Returns a filtered list of detected NCIs.  
- **`get_fragment_number(atom_idx)`**: Returns a 1-based fragment ID for labeling intermolecular interactions.

---

## Usage

### 1) Command-Line Interface

Both `molecule.py` and `molecule_nci.py` can be run directly:

```bash
# Basic usage with molecule.py
python molecule.py <xyz_file>

# Usage with molecule_nci.py (can handle multiple files)
python molecule_nci.py <xyz_file1> [xyz_file2 ...] [--debug]
```

- The code prints a summary of the system and any NCIs it detects.

### 2) In Python Scripts / Notebooks

You can import these modules (assuming they’re discoverable by Python) and call them programmatically:

```python
from molecule_nci import MoleculeWithNCIs

mol = MoleculeWithNCIs.from_xyz("path/to/myfile.xyz")
mol.detect_bonds(tolerance=0.3)
mol.find_fragments()
mol.detect_all_ncis(run_steric_clashes=False, debug=False)

for (pair, info) in mol.list_interactions():
    print(f"Interaction between atoms {pair}, type={info['type']}, dist={info['distance']:.2f}")
```

---

## Testing

### Test Files and Organization

- **`tests/test_molecule_nci.py`**: Contains a **pytest** suite verifying the code’s functionality.
- **`tests/data/`**: Contains multiple `.xyz` files (e.g. `CBr4_urea.xyz`, `CF3Br_H2O.xyz`, etc.) that exhibit specific NCIs.  

We use a **path hack** (`sys.path.insert(...)`) in the test file so Python can locate the modules without needing a formal Python package or install.

### Running the Tests

From the **`NewMoleculeNCIs`** directory:

```bash
pytest -q
```

You should see a quick summary, such as:
```
...........
11 passed in 0.12s
```
indicating all tests passed.

### What the Tests Check

1. **Number of Fragments**: Ensures that `find_fragments()` yields the expected connected components.  
2. **NCI Counts and Types**: For each reference `.xyz` file, verifies that the correct number of NCIs and specific **interaction types** (e.g., `"halogen_bond"`, `"H-bond"`, `"chalcogen_bond"`) are found.  
3. **Fine-Grained Geometry Checks**: For certain interactions, the test compares the library’s reported **distance** and **angle** to a **reference** value with a small tolerance (±0.05 Å, ±2° or so).

**Example**: In `CBr4_urea.xyz`, the test expects:
- 2 fragments
- 1 halogen bond  
- Halogen bond distance ~2.75 Å and angle ~175.4° within a tolerance

---

## Extending the Code

You can easily **add new NCI detection** routines by writing a function in `MoleculeWithNCIs` and decorating it:

```python
from molecule_nci import MoleculeWithNCIs, register_nci

@register_nci(priority=50)
def detect_custom_interaction(mol, debug=False):
    # Custom logic to detect a brand-new NCI type...
    pass
```

If you want a more specialized molecule class, you can inherit from `MoleculeWithNCIs` and override or add new methods.

