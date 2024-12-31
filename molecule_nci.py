#!/usr/bin/env python3
# molecule_nci_debug.py
#
# A subclass "MoleculeWithNCIs" that detects non-covalent interactions
# (H-bonds, sigma-hole interactions, steric clashes), prints results
# in a formatted table, uses 1-based atom numbering, can handle multiple XYZ files
# from the command line, and supports debug output.
#
# NOW with an extra column: "Fragments" showing which fragment is donor vs. acceptor
# when the interaction is intermolecular.

import sys
import math
import numpy as np

from molecule import Molecule
from elements_table import Elements

_NCI_METHODS = []

def register_nci(priority=0):
    def decorator(func):
        _NCI_METHODS.append((priority, func))
        _NCI_METHODS.sort(key=lambda x: x[0])
        return func
    return decorator

# A unified set of acceptable "sigma-hole" acceptors for halogen/chalcogen
VALID_SIGMA_HOLE_ACCEPTORS = {"N", "O", "F", "S", "CL", "BR", "I"}

class MoleculeWithNCIs(Molecule):
    """
    Subclass of Molecule that adds Non-Covalent Interaction (NCI) detection.
    """
    def __init__(self, title: str = ""):
        super().__init__(title=title)
        self.ncis = {}  # dict of (i, j) -> list of interaction dicts

    def detect_all_ncis(self, run_steric_clashes: bool = False, debug: bool = False):
        if debug:
            print("[DEBUG] detect_all_ncis: run_steric_clashes =", run_steric_clashes)
        for priority, detect_func in _NCI_METHODS:
            if detect_func.__name__ == "detect_steric_clashes" and (not run_steric_clashes):
                if debug:
                    print("[DEBUG] Skipping steric_clashes.")
                continue
            if debug:
                print(f"[DEBUG] Calling {detect_func.__name__} (priority={priority})")
            detect_func(self, debug=debug)

    def detect_bonds(self, tolerance: float = 0.3, debug: bool = False):
        n = len(self.atoms)
        if self.bond_matrix is None:
            self.bond_matrix = np.zeros((n, n), dtype=int)

        if debug:
            print(f"[DEBUG] detect_bonds: tolerance={tolerance}")

        for i in range(n):
            for j in range(i+1, n):
                sym_i = self.atoms[i].symbol.upper()
                sym_j = self.atoms[j].symbol.upper()
                r_i = Elements.covalent_radius(sym_i)
                r_j = Elements.covalent_radius(sym_j)

                dist_ij = self.distance(i, j)
                cutoff = r_i + r_j + tolerance
                if debug:
                    print(f"[DEBUG] i={i}({sym_i}), j={j}({sym_j}), dist={dist_ij:.3f}, cutoff={cutoff:.3f}")

                if dist_ij < cutoff:
                    self.bond_matrix[i, j] = 1
                    self.bond_matrix[j, i] = 1
                    if debug:
                        print(f"  [DEBUG] -> BOND DETECTED between {i}({sym_i}) and {j}({sym_j})")

    @register_nci(priority=0)
    def detect_hydrogen_bonds(self,
                              max_dist=2.5,
                              angle_cutoff=160.0,
                              donors=("N", "O", "F"),
                              acceptors=("N", "O", "F"),
                              debug=False):
        n = len(self.atoms)
        if debug:
            print("[DEBUG] detect_hydrogen_bonds")

        for i in range(n):
            if self.atoms[i].symbol.upper() != "H":
                continue
            donor_idx = None
            if self.bond_matrix is not None:
                for d in range(n):
                    if d != i and self.bond_matrix[i, d] > 0:
                        if self.atoms[d].symbol.upper() in donors:
                            donor_idx = d
                            break
            if donor_idx is None:
                continue

            for j in range(n):
                if j == i or j == donor_idx:
                    continue
                if self.atoms[j].symbol.upper() not in acceptors:
                    continue

                dist_ij = self.distance(i, j)
                if dist_ij <= max_dist:
                    ang = self.angle(donor_idx, i, j)
                    if ang >= angle_cutoff:
                        pair = (min(i, j), max(i, j))
                        if pair not in self.ncis:
                            self.ncis[pair] = []
                        self.ncis[pair].append({
                            "type": "H-bond",
                            "distance": dist_ij,
                            "angle": ang,
                            "angle_atoms": (donor_idx, i, j),
                            "intra_or_inter": self._intra_or_inter(donor_idx, j)
                        })
                        if debug:
                            print(f"[DEBUG] H-bond: i={i}, donor={donor_idx}, j={j}, dist={dist_ij:.2f}, angle={ang:.1f}")

    @register_nci(priority=1)
    def detect_sigma_hole_bonds(self, debug=False):
        """
        Halogen + Chalcogen detection using one acceptor set: VALID_SIGMA_HOLE_ACCEPTORS.
        We check angle_min=120 (halogen) or 130 (chalcogen), distance < 0.9*(vdw_i+vdw_j).
        """
        groups = [
            {"label": "halogen_bond",   "symbols": {"CL", "BR", "I"},  "angle_min": 120.0},
            {"label": "chalcogen_bond", "symbols": {"S", "SE", "TE"},  "angle_min": 130.0},
        ]
        fraction_vdw = 0.9

        n = len(self.atoms)
        for group in groups:
            bond_type = group["label"]
            donor_symbols = group["symbols"]
            angle_min = group["angle_min"]

            if debug:
                print(f"[DEBUG] Checking {bond_type} with donors={donor_symbols}")

            for i in range(n):
                sym_i = self.atoms[i].symbol.upper()
                if sym_i not in donor_symbols:
                    continue

                neighbors_i = [r for r in range(n) if r != i and self.bond_matrix[i, r] > 0]
                vdw_i = Elements.vdw_radius(sym_i)

                for j in range(n):
                    if j == i or j in neighbors_i:
                        continue
                    if self.bond_matrix[i, j] > 0:
                        continue

                    sym_j = self.atoms[j].symbol.upper()
                    if sym_j not in VALID_SIGMA_HOLE_ACCEPTORS:
                        continue

                    dist_ij = self.distance(i, j)
                    vdw_j = Elements.vdw_radius(sym_j)
                    threshold = fraction_vdw * (vdw_i + vdw_j)
                    if dist_ij < threshold:
                        # check angles
                        best_angle = None
                        best_neighbor = None
                        for neighbor in neighbors_i:
                            ax = self.angle(neighbor, i, j)
                            if ax > angle_min:
                                best_angle = ax
                                best_neighbor = neighbor
                                break
                        if best_angle:
                            pair = (min(i, j), max(i, j))
                            if pair not in self.ncis:
                                self.ncis[pair] = []
                            self.ncis[pair].append({
                                "type": bond_type,
                                "distance": dist_ij,
                                "angle": best_angle,
                                "angle_atoms": (best_neighbor, i, j),  # R–X···A
                                "intra_or_inter": self._intra_or_inter(i, j)
                            })
                            if debug:
                                print(f"[DEBUG] {bond_type}: i={i}, j={j}, dist={dist_ij:.2f}, angle={best_angle:.1f}")

    @register_nci(priority=999)
    def detect_steric_clashes(self,
                              clash_tolerance=0.4,
                              only_hydrogen_clashes=True,
                              ignore_same_neighbor=True,
                              debug=False):
        n = len(self.atoms)
        if debug:
            print("[DEBUG] detect_steric_clashes")

        for i in range(n):
            for j in range(i+1, n):
                sym_i = self.atoms[i].symbol.upper()
                sym_j = self.atoms[j].symbol.upper()

                if only_hydrogen_clashes and (sym_i != "H" or sym_j != "H"):
                    continue

                if self.bond_matrix[i, j] > 0:
                    continue

                if ignore_same_neighbor and self._share_a_neighbor(i, j):
                    continue

                dist_ij = self.distance(i, j)
                vdw_i = Elements.vdw_radius(sym_i)
                vdw_j = Elements.vdw_radius(sym_j)
                if dist_ij < (vdw_i + vdw_j - clash_tolerance):
                    pair = (i, j)
                    existing = self.ncis.get((i, j), []) or self.ncis.get((j, i), [])
                    skip_clash = any(
                        x["type"] in ["H-bond", "halogen_bond", "chalcogen_bond"]
                        for x in existing
                    )
                    if not skip_clash:
                        if pair not in self.ncis:
                            self.ncis[pair] = []
                        self.ncis[pair].append({
                            "type": "steric_clash",
                            "distance": dist_ij,
                            "angle": None,
                            "angle_atoms": None,
                            "intra_or_inter": self._intra_or_inter(i, j)
                        })
                        if debug:
                            print(f"[DEBUG] steric_clash: i={i}, j={j}, dist={dist_ij:.2f}")

    def _share_a_neighbor(self, i: int, j: int) -> bool:
        if self.bond_matrix is None:
            return False
        n = len(self.atoms)
        neighbors_i = []
        neighbors_j = []
        for k in range(n):
            if k == i or k == j:
                continue
            if self.bond_matrix[i, k] > 0:
                neighbors_i.append(k)
            if self.bond_matrix[j, k] > 0:
                neighbors_j.append(k)
        return len(set(neighbors_i).intersection(set(neighbors_j))) > 0

    def _intra_or_inter(self, i: int, j: int) -> str:
        if not self.fragments:
            return "unknown"
        frag_i = frag_j = None
        for f_id, members in self.fragments.items():
            if i in members:
                frag_i = f_id
            if j in members:
                frag_j = f_id
        if frag_i is None or frag_j is None:
            return "unknown"
        return "intra" if (frag_i == frag_j) else "inter"

    def list_interactions(self, interaction_type=None, fragment_scope=None):
        results = []
        for pair, all_info in self.ncis.items():
            for info in all_info:
                if interaction_type and info["type"] != interaction_type:
                    continue
                if fragment_scope and info["intra_or_inter"] != fragment_scope:
                    continue
                results.append((pair, info))
        return results

    #
    # OPTIONAL: Utility to find which fragment an atom belongs to, returning a 1-based number
    #
    def get_fragment_number(self, atom_idx: int) -> int:
        """
        Return the fragment number (1-based) of 'atom_idx'.
        If not found, return 0 or -1
        """
        if not self.fragments:
            return 0
        for f_id, members in self.fragments.items():
            if atom_idx in members:
                return f_id + 1  # 1-based
        return 0


#
# Helper function to figure out donor vs. acceptor fragment
#
def get_fragments_for_interaction(mol: MoleculeWithNCIs, i: int, j: int, info: dict):
    """
    Returns a string describing which fragments are the donor and acceptor
    for an *intermolecular* interaction. For intramolecular, returns a single fragment.

    - If "intra", e.g. "Frag1"
    - If "inter", e.g. "Frag1->Frag2"
      * H-bond: donor is angle_atoms[0], acceptor is angle_atoms[2]
      * halogen/chalcogen: donor = 'i', acceptor = 'j'
      * clashes: no donor/acceptor concept, but we can do i->j anyway
    """
    scope = info.get("intra_or_inter", "unknown")
    if scope != "inter":
        # intramolecular
        frag = mol.get_fragment_number(i)
        return f"Frag{frag}"  # or check if j has same
    # "inter" -> we figure out which is donor vs. acceptor
    interaction_type = info["type"]

    # default: just do i->j
    donor_atom = i
    acceptor_atom = j

    if interaction_type == "H-bond":
        # angle_atoms = (donor_idx, H, acceptor)
        angle_atoms = info.get("angle_atoms", None)
        if angle_atoms and len(angle_atoms) == 3:
            donor_atom = angle_atoms[0]
            acceptor_atom = angle_atoms[2]

    elif interaction_type in ["halogen_bond", "chalcogen_bond"]:
        # angle_atoms = (neighbor, i, j), i is the halogen/chalcogen, j is acceptor
        # so donor = i, acceptor = j
        # (the default i->j is correct if i < j, but let's be explicit)
        donor_atom = info["angle_atoms"][1]
        acceptor_atom = info["angle_atoms"][2]
    # else: steric_clash or something else => just do i->j

    donor_frag = mol.get_fragment_number(donor_atom)
    acceptor_frag = mol.get_fragment_number(acceptor_atom)
    return f"Frag{donor_frag}->Frag{acceptor_frag}"


def print_interactions_table(mol: MoleculeWithNCIs, interactions):
    """
    Print a tabular summary of all interactions found, with:
      - Type
      - Pair (S1-O5, etc.)
      - Dist(Å)
      - Angle(°)
      - AngleAtoms (e.g. C1–S2–O6)
      - Fragments (e.g. 'Frag1->Frag2' for inter, or 'Frag1' if intra)
      - Scope
    """
    print(f"{'Type':<14} {'Pair':<10} {'Dist(Å)':>8} {'Angle(°)':>9} {'AngleAtoms':<15} {'Fragments':<12} {'Scope':>6}")
    print("-" * 100)

    for (i, j), info in interactions:
        typ   = info["type"]
        dist  = info.get("distance", None)
        angle = info.get("angle", None)
        angle_atoms = info.get("angle_atoms", None)
        scope = info.get("intra_or_inter", "unknown")

        sym_i = mol.atoms[i].symbol
        sym_j = mol.atoms[j].symbol
        pair_str = f"{sym_i}{i+1}-{sym_j}{j+1}"

        dist_str  = f"{dist:.2f}"   if dist  is not None else "N/A"
        angle_str = f"{angle:.1f}" if angle is not None else "N/A"

        angle_atoms_str = "N/A"
        if angle_atoms and len(angle_atoms) == 3:
            a1, a2, a3 = angle_atoms
            sym_a1 = mol.atoms[a1].symbol
            sym_a2 = mol.atoms[a2].symbol
            sym_a3 = mol.atoms[a3].symbol
            angle_atoms_str = f"{sym_a1}{a1+1}-{sym_a2}{a2+1}-{sym_a3}{a3+1}"

        # NEW: figure out donor/acceptor fragment labeling if it's inter
        fragments_str = get_fragments_for_interaction(mol, i, j, info)

        print(f"{typ:<14} {pair_str:<10} {dist_str:>8} {angle_str:>9} {angle_atoms_str:<15} {fragments_str:<12} {scope:>6}")


def main():
    args = sys.argv[1:]
    debug_mode = False
    if "--debug" in args:
        debug_mode = True
        args.remove("--debug")

    if len(args) < 1:
        print("Usage: python molecule_nci_debug.py <xyz_file1> [xyz_file2 ...] [--debug]")
        sys.exit(1)

    for xyz_file in args:
        print(f"\n=== Processing {xyz_file} ===")
        mol = MoleculeWithNCIs.from_xyz(xyz_file)

        mol.detect_bonds(tolerance=0.3, debug=debug_mode)
        mol.find_fragments()

        mol.detect_all_ncis(run_steric_clashes=False, debug=debug_mode)

        print("\n--- Molecule Summary ---")
        print(mol.summary())

        all_ncis = mol.list_interactions()
        if not all_ncis:
            print("\nNo NCIs found.")
        else:
            print("\n--- Non-Covalent Interactions Detected ---")
            print_interactions_table(mol, all_ncis)


if __name__ == "__main__":
    main()
