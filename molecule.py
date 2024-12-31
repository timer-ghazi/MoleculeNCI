#!/usr/bin/env python3
# molecule.py
#
# A base Molecule class for reading atoms from an XYZ file, building
# a covalent adjacency matrix using data from elements_table.py,
# and providing geometry utilities (distance, angle, dihedral).
#
# This example includes:
#  1) An inline Atom dataclass.
#  2) The Molecule class with the methods described.
#  3) A command-line "main" entry point for testing.

import sys
import math
import numpy as np
from dataclasses import dataclass
# Make sure elements_table.py is in the same directory or installable path
from elements_table import Elements


@dataclass
class Atom:
    """
    A lightweight container for atomic data.
    """
    symbol: str
    x: float
    y: float
    z: float
    charge: float = 0.0


class Molecule:
    """
    A base class representing a molecular system with covalent bonding
    and fragment detection. Coordinates are in Å by convention.

    Key attributes:
        title (str): A descriptive name for the molecule.
        atoms (List[Atom]): List of Atom objects.
        bond_matrix (np.ndarray): 2D adjacency matrix of shape (n, n), storing bond information.
        fragments (Dict[int, List[int]]): Mapping of fragment_id -> list of atom indices.
    """

    def __init__(self, title: str = ""):
        self.title = title
        self.atoms = []               # type: List[Atom]
        self.bond_matrix = None       # type: np.ndarray
        self.fragments = {}           # type: Dict[int, List[int]]

    @classmethod
    def from_xyz(cls, filename: str) -> "Molecule":
        """
        Parse an XYZ file to create a Molecule instance.

        XYZ format reminder:
          1) First line = number_of_atoms
          2) Second line = optional comment (we store it in .title)
          3) Subsequent lines = symbol x y z

        Returns:
            Molecule object with self.atoms populated.
            bond_matrix is initialized to zeros, but not auto-detected yet.
        """
        with open(filename, 'r') as f:
            lines = f.readlines()

        # 1) number_of_atoms from first line
        num_atoms = int(lines[0].strip())
        # 2) comment line (optional)
        comment_line = lines[1].strip()

        # Create the Molecule
        mol = cls(title=comment_line)

        # 3) parse each subsequent line
        idx = 2
        for i in range(num_atoms):
            parts = lines[idx].split()
            idx += 1
            symbol = parts[0]
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
            atom = Atom(symbol, x, y, z)
            mol.atoms.append(atom)

        # Initialize bond_matrix
        n = len(mol.atoms)
        mol.bond_matrix = np.zeros((n, n), dtype=float)

        return mol

    def detect_bonds(self, tolerance: float = 0.3):
        """
        Fills self.bond_matrix by comparing interatomic distances
        to the sum of covalent radii (plus a tolerance).

        :param tolerance: Additional margin (in Å) on top of covalent radius sums.
                          Adjust as needed for borderline cases.
        """
        n = len(self.atoms)
        for i in range(n):
            for j in range(i + 1, n):
                dist_ij = self.distance(i, j)
                r1 = Elements.covalent_radius(self.atoms[i].symbol, order="single")
                r2 = Elements.covalent_radius(self.atoms[j].symbol, order="single")
                r_sum = r1 + r2 + tolerance

                if dist_ij <= r_sum:
                    # Simple approach: call it a single bond
                    self.bond_matrix[i, j] = 1
                    self.bond_matrix[j, i] = 1

                    # More advanced logic for double/triple bonds could go here

    def find_fragments(self):
        """
        Identify connected components ("fragments") in the bond_matrix.
        A DFS approach populates self.fragments = {frag_id: [atom_indices]}.
        """
        n = len(self.atoms)
        visited = set()
        frag_id = 0
        self.fragments = {}

        for start_atom in range(n):
            if start_atom not in visited:
                stack = [start_atom]
                connected = []
                while stack:
                    current = stack.pop()
                    if current not in visited:
                        visited.add(current)
                        connected.append(current)

                        # look for neighbors
                        for neigh in range(n):
                            if self.bond_matrix[current, neigh] > 0.0 and neigh not in visited:
                                stack.append(neigh)

                # store the connected fragment
                self.fragments[frag_id] = connected
                frag_id += 1

    def distance(self, i: int, j: int) -> float:
        r""" 
        Returns the Euclidean distance between atoms i and j (in ang).

        $$
        d_{ij} = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2}
        $$
        """
        ax, ay, az = self.atoms[i].x, self.atoms[i].y, self.atoms[i].z
        bx, by, bz = self.atoms[j].x, self.atoms[j].y, self.atoms[j].z
        return math.sqrt((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)

    def angle(self, i: int, j: int, k: int, degrees: bool = True) -> float:
        """
        Returns the angle at atom j formed by (i - j - k).
        If degrees=False, return radians.

        $$
        \\theta = \\cos^{-1} \\Bigl(
            \\frac{(\\mathbf{r}_i - \\mathbf{r}_j) \\cdot (\\mathbf{r}_k - \\mathbf{r}_j)}
                 {\\|\\mathbf{r}_i - \\mathbf{r}_j\\| \\, \\|\\mathbf{r}_k - \\mathbf{r}_j\\|}
        \\Bigr)
        $$
        """
        r_i = np.array([self.atoms[i].x, self.atoms[i].y, self.atoms[i].z])
        r_j = np.array([self.atoms[j].x, self.atoms[j].y, self.atoms[j].z])
        r_k = np.array([self.atoms[k].x, self.atoms[k].y, self.atoms[k].z])

        v_ji = r_i - r_j
        v_jk = r_k - r_j

        dot_val = np.dot(v_ji, v_jk)
        mag_ji = np.linalg.norm(v_ji)
        mag_jk = np.linalg.norm(v_jk)

        cos_theta = dot_val / (mag_ji * mag_jk)

        # numerical safety
        cos_theta = max(min(cos_theta, 1.0), -1.0)

        theta_radians = np.arccos(cos_theta)
        return np.degrees(theta_radians) if degrees else theta_radians

    def dihedral(self, i: int, j: int, k: int, l: int, degrees: bool = True) -> float:
        """
        Returns the dihedral angle formed by atoms (i - j - k - l).
        If degrees=False, returns radians.

        A standard approach uses cross products to find normal vectors
        and the sign of the torsion from the 'atan2' of those vectors.
        """
        r_i = np.array([self.atoms[i].x, self.atoms[i].y, self.atoms[i].z])
        r_j = np.array([self.atoms[j].x, self.atoms[j].y, self.atoms[j].z])
        r_k = np.array([self.atoms[k].x, self.atoms[k].y, self.atoms[k].z])
        r_l = np.array([self.atoms[l].x, self.atoms[l].y, self.atoms[l].z])

        b1 = r_i - r_j
        b2 = r_k - r_j
        b3 = r_l - r_k

        # normal vectors
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)

        # handle degenerate vectors (avoid /0)
        if np.linalg.norm(n1) == 0 or np.linalg.norm(n2) == 0:
            return 0.0  # or raise an exception

        # normalize
        n1 /= np.linalg.norm(n1)
        n2 /= np.linalg.norm(n2)

        # cross product to get sign
        m1 = np.cross(n1, b2 / np.linalg.norm(b2))

        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        angle_radians = np.arctan2(y, x)

        return np.degrees(angle_radians) if degrees else angle_radians

    def summary(self) -> str:
        """
        Returns a human-readable summary of the molecule:
        number of atoms, fragments, etc.
        """
        lines = []
        lines.append(f"Title: {self.title}")
        lines.append(f"Number of atoms: {len(self.atoms)}")
        if self.bond_matrix is not None:
            lines.append("Bond matrix is initialized.")
        if self.fragments:
            lines.append(f"Number of fragments: {len(self.fragments)}")
        else:
            lines.append("Fragments not yet determined.")
        return "\n".join(lines)


# Quick test / CLI usage

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python molecule.py <xyz_file>")
        sys.exit(1)

    xyz_file = sys.argv[1]

    # Build the Molecule from the XYZ file
    mol = Molecule.from_xyz(xyz_file)

    # Detect bonds
    mol.detect_bonds(tolerance=0.3)

    # Find connected components (fragments)
    mol.find_fragments()

    # Print a summary
    print(mol.summary())

    # --- New snippet: print bond matrix ---
    print("\nBond matrix:")
    n = len(mol.atoms)
    # Print header (atom indices), each index right-aligned in 2 characters
    header = "    " + " ".join([f"{j:>2}" for j in range(n)])
    print(header)
    
    # For each row, print the row index, then each column entry
    for i in range(n):
        row_str = " ".join(
            f"{int(mol.bond_matrix[i, j]):>2}"
            if mol.bond_matrix[i, j] != 0
            else f"{'•':>2}"
            for j in range(n)
        )
        print(f"{i:>2}  {row_str}")
    

    # Example geometry checks
    if n >= 2:
        d_01 = mol.distance(0, 1)
        print(f"\nDistance between atom 0 and 1: {d_01:.3f} Å")

    if n >= 4:
        torsion_angle = mol.dihedral(0, 1, 2, 3)
        print(f"Dihedral angle (0-1-2-3): {torsion_angle:.2f} degrees")

