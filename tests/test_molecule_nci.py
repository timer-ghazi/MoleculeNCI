# tests/test_molecule_nci.py

import sys
import os

# This line inserts the parent directory of "tests" into sys.path,
# so Python can see molecule_nci.py, molecule.py, elements_table.py, etc.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pytest
from molecule_nci import MoleculeWithNCIs

###############################################################################
# 1) Basic Expectations: number of fragments, NCI types, total NCI count
###############################################################################
REFERENCE_DATA = {
    "CBr4_urea.xyz": {
        "num_fragments": 2,
        "nci_types": ["halogen_bond"],
        "num_ncis": 1
    },
    "CF3Br_H2O.xyz": {
        "num_fragments": 2,
        "nci_types": ["halogen_bond"],
        "num_ncis": 1
    },
    "CF3I_NH3.xyz": {
        "num_fragments": 2,
        "nci_types": ["halogen_bond"],
        "num_ncis": 1
    },
    "CF3I_acetone.xyz": {
        "num_fragments": 2,
        "nci_types": ["halogen_bond"],
        "num_ncis": 1
    },
    "CH3CH2COOH-HOCH2CH3.xyz": {
        "num_fragments": 2,
        "nci_types": ["H-bond"],
        "num_ncis": 1
    },
    "F2S_NH3.xyz": {
        "num_fragments": 2,
        "nci_types": ["chalcogen_bond"],
        "num_ncis": 1
    },
    "F2Se_H2O.xyz": {
        "num_fragments": 2,
        "nci_types": ["chalcogen_bond"],
        "num_ncis": 1
    },
    "F3C2S_HCN.xyz": {
        "num_fragments": 2,
        "nci_types": [],
        "num_ncis": 0
    },
    "FCl_HCN.xyz": {
        "num_fragments": 2,
        "nci_types": ["halogen_bond"],
        "num_ncis": 1
    },
    "HO(CH2)3COOH.xyz": {
        "num_fragments": 1,
        "nci_types": [],
        "num_ncis": 0
    },
    "SO2_NH3.xyz": {
        "num_fragments": 2,
        "nci_types": [],
        "num_ncis": 0
    },
}

###############################################################################
# 2) Fine-Grained Checks: expected distance, angle, etc. for a specific NCI.
#    Each list item includes:
#      - "type": "halogen_bond" | "H-bond" | "chalcogen_bond" | ...
#      - "pair": (i, j) in zero-based indices (i < j)
#      - "distance": expected distance
#      - "distance_tol": allowable deviation in Å
#      - "angle": expected angle in degrees
#      - "angle_tol": allowable deviation in degrees
###############################################################################
EXPLICIT_NCI_DATA = {
    "CBr4_urea.xyz": [
        {
            "type": "halogen_bond",
            "pair": (3, 7),          # from 1-based "Br4-O8"
            "distance": 2.75,
            "distance_tol": 0.05,
            "angle": 175.4,
            "angle_tol": 2.0
        }
    ],
    "CF3Br_H2O.xyz": [
        {
            "type": "halogen_bond",
            "pair": (4, 5),         # from 1-based "Br5-O6"
            "distance": 2.99,
            "distance_tol": 0.05,
            "angle": 178.1,
            "angle_tol": 2.0
        }
    ],
    "CF3I_NH3.xyz": [
        {
            "type": "halogen_bond",
            "pair": (4, 5),         # from 1-based "I5-N6"
            "distance": 2.91,
            "distance_tol": 0.05,
            "angle": 179.3,
            "angle_tol": 2.0
        }
    ],
    "CF3I_acetone.xyz": [
        {
            "type": "halogen_bond",
            "pair": (4, 7),         # from 1-based "I5-O8"
            "distance": 2.90,
            "distance_tol": 0.05,
            "angle": 179.2,
            "angle_tol": 2.0
        }
    ],
    "CH3CH2COOH-HOCH2CH3.xyz": [
        {
            "type": "H-bond",
            "pair": (9, 11),        # from 1-based "O10-H12"
            "distance": 1.75,
            "distance_tol": 0.05,
            "angle": 171.7,
            "angle_tol": 5.0
        }
    ],
    "F2S_NH3.xyz": [
        {
            "type": "chalcogen_bond",
            "pair": (0, 1),         # from 1-based "S1-N2"
            "distance": 2.47,
            "distance_tol": 0.05,
            "angle": 173.7,
            "angle_tol": 3.0
        }
    ],
    "F2Se_H2O.xyz": [
        {
            "type": "chalcogen_bond",
            "pair": (1, 3),         # from 1-based "Se2-O4"
            "distance": 2.56,
            "distance_tol": 0.05,
            "angle": 173.0,
            "angle_tol": 3.0
        }
    ],
    "FCl_HCN.xyz": [
        {
            "type": "halogen_bond",
            "pair": (1, 4),         # from 1-based "Cl2-N5"
            "distance": 2.55,
            "distance_tol": 0.05,
            "angle": 179.5,
            "angle_tol": 2.0
        }
    ],
    # Others have no interactions => no explicit checks
}

###############################################################################
# Test Function
###############################################################################
@pytest.mark.parametrize("xyz_file", REFERENCE_DATA.keys())
def test_nci_detection(xyz_file):
    """
    Test that MoleculeWithNCIs correctly identifies:
      1) The expected number of fragments
      2) The expected overall NCI types
      3) The expected total number of NCIs
      4) (Optionally) Specific distance & angle checks for certain interactions
    """
    data = REFERENCE_DATA[xyz_file]
    xyz_path = os.path.join("tests", "data", xyz_file)

    # 1) Build the molecule
    mol = MoleculeWithNCIs.from_xyz(xyz_path)

    # 2) Detect bonds and fragments
    mol.detect_bonds(tolerance=0.3)
    mol.find_fragments()

    # 3) Detect all NCIs (skip steric clashes for these tests)
    mol.detect_all_ncis(run_steric_clashes=False, debug=False)

    # Check number of fragments
    assert len(mol.fragments) == data["num_fragments"], (
        f"[{xyz_file}] Expected {data['num_fragments']} fragment(s), "
        f"found {len(mol.fragments)}"
    )

    # Retrieve all NCIs
    all_ncis = mol.list_interactions()

    # Check total number of NCIs
    assert len(all_ncis) == data["num_ncis"], (
        f"[{xyz_file}] Expected {data['num_ncis']} NCI(s), found {len(all_ncis)}"
    )

    # Check that the set of NCI types includes everything in data["nci_types"]
    found_types = {info["type"] for _, info in all_ncis}
    for expected_type in data["nci_types"]:
        assert expected_type in found_types, (
            f"[{xyz_file}] Expected NCI type '{expected_type}', but none found. "
            f"Detected types were: {found_types}"
        )

    # 4) Fine-grained checks: distances and angles
    if xyz_file in EXPLICIT_NCI_DATA:
        for expected_nci in EXPLICIT_NCI_DATA[xyz_file]:
            check_one_nci(mol, all_ncis, xyz_file, expected_nci)

def check_one_nci(mol, all_ncis, xyz_file, expected_nci):
    """
    Helper function that finds a specific NCI (by type + pair) and 
    checks distance and angle against expected values.
    """
    exp_type = expected_nci["type"]
    exp_pair = expected_nci["pair"]  # zero-based tuple (i, j)
    exp_dist = expected_nci["distance"]
    dist_tol = expected_nci["distance_tol"]
    exp_angle = expected_nci["angle"]
    angle_tol = expected_nci["angle_tol"]

    found_match = False

    for (pair, info) in all_ncis:
        # Pair is stored as (min, max), so let's compare sorted indices:
        if tuple(sorted(pair)) == exp_pair and info["type"] == exp_type:
            # Check distance
            dist = info.get("distance")
            assert dist is not None, (
                f"[{xyz_file}] Missing distance for {exp_type} NCI {exp_pair}"
            )
            diff_dist = abs(dist - exp_dist)
            assert diff_dist <= dist_tol, (
                f"[{xyz_file}] {exp_type} NCI {exp_pair} distance off by {diff_dist:.2f} Å; "
                f"expected ~{exp_dist}, got {dist:.2f}"
            )

            # Check angle
            angle = info.get("angle")
            assert angle is not None, (
                f"[{xyz_file}] Missing angle for {exp_type} NCI {exp_pair}"
            )
            diff_angle = abs(angle - exp_angle)
            assert diff_angle <= angle_tol, (
                f"[{xyz_file}] {exp_type} NCI {exp_pair} angle off by {diff_angle:.2f}°; "
                f"expected ~{exp_angle}, got {angle:.2f}"
            )

            found_match = True
            break

    # If we never found the match, fail the test
    assert found_match, (
        f"[{xyz_file}] Could not find NCI with type '{exp_type}' and pair {exp_pair}"
    )
