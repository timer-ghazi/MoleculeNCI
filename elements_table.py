#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dec 25 removed None's from the data, removed mulliken
"""
elements_table.py

Single-file Python module containing element data and a class `Elements`
that provides convenient classmethods for lookups, such as:

  - mass(symbol, unit="u")
  - vdw_radius(symbol, unit="Ang")
  - covalent_radius(symbol, order="single", source="cordero", unit="Ang")
  - electronegativity(symbol, scale="pauling")
  - atomic_number(symbol)
  - symbol(atomic_number)
  - name(symbol)
  - period(symbol)
  - group(symbol)
  - classification(symbol)
  - is_valid(symbol)
  - list_symbols()

Data is stored in `_ELEMENTS`, originating from combined sources with
all distances in Å and masses in atomic mass units (u).

Conversion factors are applied based on user-provided unit strings.
"""

# #############################################################################
# 1) Combined dictionary of element data (excerpts from combined_elements.py)
# #############################################################################

#
# Sources:
# 1. Van der Waals Radii:
#    - Source: https://en.wikipedia.org/wiki/Van_der_Waals_radius
#    - All values in Ångströms (Å)
#
# 2. Electronegativity:
#    - Three different scales:
#      * Pauling: dimensionless
#      * Mulliken: electron volts (eV)
#      * Allen: electron volts (eV)
#
# 3. Covalent Radii:
#    - Pyykko Parameters:
#      * Source: Pyykkö, P. & Atsumi, M.
#      * "Molecular Double-Bond Covalent Radii for Elements Li–E112"
#      * Chemistry: A European Journal 15 (46): 12770–12779
#      * DOI: dx.doi.org/10.1002/chem.200901472
#    - Cordero Parameters:
#      * Source: Cordero, B. et al.
#      * "Covalent radii revisited"
#      * Dalton Trans. (21): 2832–2838
#      * DOI: http://dx.doi.org/10.1039/b801115j
#    - All values in Ångströms (Å)
#
# 4. Atomic Masses:
#    - All values in unified atomic mass units (u)
#
# 5. Element Metadata:
#    - Includes: full names, atomic numbers, groups, periods, and classifications
#    - Classification categories: Nonmetal, Noble Gas, Alkali Metal,
#      Alkaline Earth Metal, Metalloid, Post-transition Metal,
#      Halogen, Transition Metal
#
# Data Organization:
# - Elements from Period 1-6
# - Main Group Elements followed by Transition Metals (Period 4 and 5)
# - All measurements in standard units (Å, eV, u)
#
# Created by combining multiple parameter files into a unified structure
# Last Updated: December 23, 2024
# Combined elements data from all source files

_ELEMENTS = {
    "H": {
        "atomic_number": 1,
        "name": "Hydrogen",
        "mass": 1.008,  # in u
        "vdw_radius": 1.2,  # in Å
        "bond_params": {
            "cordero": {'single': 0.31},
            "pyykko": {'single': 0.32},
        },
        "electronegativity": {
            "pauling": 2.2,
            "allen": 2.3,
        },
        "group": 1,
        "period": 1,
        "classification": "nonmetal",
    },
    "He": {
        "atomic_number": 2,
        "name": "Helium",
        "mass": 4.002602,  # in u
        "vdw_radius": 1.4,  # in Å
        "bond_params": {
            "cordero": {'single': 0.28},
            "pyykko": {'single': 0.46},
        },
        "electronegativity": {
            "allen": 4.16,
        },
        "group": 18,
        "period": 1,
        "classification": "noble gas",
    },
    "Li": {
        "atomic_number": 3,
        "name": "Lithium",
        "mass": 6.94,  # in u
        "vdw_radius": 1.82,  # in Å
        "bond_params": {
            "cordero": {'single': 1.28},
            "pyykko": {'single': 1.33, 'double': 1.24},
        },
        "electronegativity": {
            "pauling": 0.98,
            "allen": 0.912,
        },
        "group": 1,
        "period": 2,
        "classification": "alkali metal",
    },
    "Be": {
        "atomic_number": 4,
        "name": "Beryllium",
        "mass": 9.0122,  # in u
        "vdw_radius": 1.53,  # in Å
        "bond_params": {
            "cordero": {'single': 0.96},
            "pyykko": {'single': 1.02, 'double': 0.9, 'triple': 0.85},
        },
        "electronegativity": {
            "pauling": 1.57,
            "allen": 1.576,
        },
        "group": 2,
        "period": 2,
        "classification": "alkaline earth metal",
    },
    "B": {
        "atomic_number": 5,
        "name": "Boron",
        "mass": 10.81,  # in u
        "vdw_radius": 1.92,  # in Å
        "bond_params": {
            "cordero": {'single': 0.84},
            "pyykko": {'single': 0.85, 'double': 0.78, 'triple': 0.73},
        },
        "electronegativity": {
            "pauling": 2.04,
            "allen": 2.051,
        },
        "group": 13,
        "period": 2,
        "classification": "metalloid",
    },
    "C": {
        "atomic_number": 6,
        "name": "Carbon",
        "mass": 12.011,  # in u
        "vdw_radius": 1.7,  # in Å
        "bond_params": {
            "cordero": {'single': 0.76},
            "pyykko": {'single': 0.75, 'double': 0.67, 'triple': 0.6},
        },
        "electronegativity": {
            "pauling": 2.55,
            "allen": 2.544,
        },
        "group": 14,
        "period": 2,
        "classification": "nonmetal",
    },
    "N": {
        "atomic_number": 7,
        "name": "Nitrogen",
        "mass": 14.007,  # in u
        "vdw_radius": 1.55,  # in Å
        "bond_params": {
            "cordero": {'single': 0.71},
            "pyykko": {'single': 0.71, 'double': 0.6, 'triple': 0.54},
        },
        "electronegativity": {
            "pauling": 3.04,
            "allen": 3.066,
        },
        "group": 15,
        "period": 2,
        "classification": "nonmetal",
    },
    "O": {
        "atomic_number": 8,
        "name": "Oxygen",
        "mass": 15.999,  # in u
        "vdw_radius": 1.52,  # in Å
        "bond_params": {
            "cordero": {'single': 0.66},
            "pyykko": {'single': 0.63, 'double': 0.57, 'triple': 0.53},
        },
        "electronegativity": {
            "pauling": 3.44,
            "allen": 3.61,
        },
        "group": 16,
        "period": 2,
        "classification": "nonmetal",
    },
    "F": {
        "atomic_number": 9,
        "name": "Fluorine",
        "mass": 18.998403163,  # in u
        "vdw_radius": 1.47,  # in Å
        "bond_params": {
            "cordero": {'single': 0.57},
            "pyykko": {'single': 0.64, 'double': 0.59, 'triple': 0.53},
        },
        "electronegativity": {
            "pauling": 3.98,
            "allen": 4.193,
        },
        "group": 17,
        "period": 2,
        "classification": "halogen",
    },
    "Ne": {
        "atomic_number": 10,
        "name": "Neon",
        "mass": 20.1797,  # in u
        "vdw_radius": 1.54,  # in Å
        "bond_params": {
            "cordero": {'single': 0.58},
            "pyykko": {'single': 0.67, 'double': 0.96},
        },
        "electronegativity": {
            "allen": 4.787,
        },
        "group": 18,
        "period": 2,
        "classification": "noble gas",
    },
    "Na": {
        "atomic_number": 11,
        "name": "Sodium",
        "mass": 22.98976928,  # in u
        "vdw_radius": 2.27,  # in Å
        "bond_params": {
            "cordero": {'single': 1.66},
            "pyykko": {'single': 1.55, 'double': 1.6},
        },
        "electronegativity": {
            "pauling": 0.93,
            "allen": 0.869,
        },
        "group": 1,
        "period": 3,
        "classification": "alkali metal",
    },
    "Mg": {
        "atomic_number": 12,
        "name": "Magnesium",
        "mass": 24.305,  # in u
        "vdw_radius": 1.73,  # in Å
        "bond_params": {
            "cordero": {'single': 1.41},
            "pyykko": {'single': 1.39, 'double': 1.32, 'triple': 1.27},
        },
        "electronegativity": {
            "pauling": 1.31,
            "allen": 1.293,
        },
        "group": 2,
        "period": 3,
        "classification": "alkaline earth metal",
    },
    "Al": {
        "atomic_number": 13,
        "name": "Aluminum",
        "mass": 26.9815385,  # in u
        "vdw_radius": 1.84,  # in Å
        "bond_params": {
            "cordero": {'single': 1.21},
            "pyykko": {'single': 1.26, 'double': 1.13, 'triple': 1.11},
        },
        "electronegativity": {
            "pauling": 1.61,
            "allen": 1.613,
        },
        "group": 13,
        "period": 3,
        "classification": "post-transition metal",
    },
    "Si": {
        "atomic_number": 14,
        "name": "Silicon",
        "mass": 28.085,  # in u
        "vdw_radius": 2.1,  # in Å
        "bond_params": {
            "cordero": {'single': 1.11},
            "pyykko": {'single': 1.16, 'double': 1.07, 'triple': 1.02},
        },
        "electronegativity": {
            "pauling": 1.9,
            "allen": 1.916,
        },
        "group": 14,
        "period": 3,
        "classification": "metalloid",
    },
    "P": {
        "atomic_number": 15,
        "name": "Phosphorus",
        "mass": 30.973761998,  # in u
        "vdw_radius": 1.8,  # in Å
        "bond_params": {
            "cordero": {'single': 1.07},
            "pyykko": {'single': 1.11, 'double': 1.02, 'triple': 0.94},
        },
        "electronegativity": {
            "pauling": 2.19,
            "allen": 2.253,
        },
        "group": 15,
        "period": 3,
        "classification": "nonmetal",
    },
    "S": {
        "atomic_number": 16,
        "name": "Sulfur",
        "mass": 32.06,  # in u
        "vdw_radius": 1.8,  # in Å
        "bond_params": {
            "cordero": {'single': 1.05},
            "pyykko": {'single': 1.03, 'double': 0.94, 'triple': 0.95},
        },
        "electronegativity": {
            "pauling": 2.58,
            "allen": 2.589,
        },
        "group": 16,
        "period": 3,
        "classification": "nonmetal",
    },
    "Cl": {
        "atomic_number": 17,
        "name": "Chlorine",
        "mass": 35.45,  # in u
        "vdw_radius": 1.75,  # in Å
        "bond_params": {
            "cordero": {'single': 1.02},
            "pyykko": {'single': 0.99, 'double': 0.95, 'triple': 0.93},
        },
        "electronegativity": {
            "pauling": 3.16,
            "allen": 2.869,
        },
        "group": 17,
        "period": 3,
        "classification": "halogen",
    },
    "Ar": {
        "atomic_number": 18,
        "name": "Argon",
        "mass": 39.948,  # in u
        "vdw_radius": 1.88,  # in Å
        "bond_params": {
            "cordero": {'single': 1.06},
            "pyykko": {'single': 0.96, 'double': 1.07, 'triple': 0.96},
        },
        "electronegativity": {
            "allen": 3.242,
        },
        "group": 18,
        "period": 3,
        "classification": "noble gas",
    },
    "K": {
        "atomic_number": 19,
        "name": "Potassium",
        "mass": 39.0983,  # in u
        "vdw_radius": 2.75,  # in Å
        "bond_params": {
            "cordero": {'single': 2.03},
            "pyykko": {'single': 1.96, 'double': 1.93},
        },
        "electronegativity": {
            "pauling": 0.82,
            "allen": 0.734,
        },
        "group": 1,
        "period": 4,
        "classification": "alkali metal",
    },
    "Ca": {
        "atomic_number": 20,
        "name": "Calcium",
        "mass": 40.078,  # in u
        "vdw_radius": 2.31,  # in Å
        "bond_params": {
            "cordero": {'single': 1.76},
            "pyykko": {'single': 1.71, 'double': 1.47, 'triple': 1.33},
        },
        "electronegativity": {
            "pauling": 1.0,
            "allen": 1.034,
        },
        "group": 2,
        "period": 4,
        "classification": "alkaline earth metal",
    },
    "Sc": {
        "atomic_number": 21,
        "name": "Scandium",
        "mass": 44.955908,  # in u
        "vdw_radius": 2.11,  # in Å
        "bond_params": {
            "cordero": {'single': 1.7},
            "pyykko": {'single': 1.48, 'double': 1.16, 'triple': 1.14},
        },
        "electronegativity": {
            "pauling": 1.36,
            "allen": 1.19,
        },
        "group": 3,
        "period": 4,
        "classification": "transition metal",
    },
    "Ti": {
        "atomic_number": 22,
        "name": "Titanium",
        "mass": 47.867,  # in u
        "vdw_radius": 2.0,  # in Å
        "bond_params": {
            "cordero": {'single': 1.6},
            "pyykko": {'single': 1.36, 'double': 1.17, 'triple': 1.08},
        },
        "electronegativity": {
            "pauling": 1.54,
            "allen": 1.38,
        },
        "group": 4,
        "period": 4,
        "classification": "transition metal",
    },
    "V": {
        "atomic_number": 23,
        "name": "Vanadium",
        "mass": 50.9415,  # in u
        "vdw_radius": 2.0,  # in Å
        "bond_params": {
            "cordero": {'single': 1.53},
            "pyykko": {'single': 1.34, 'double': 1.12, 'triple': 1.06},
        },
        "electronegativity": {
            "pauling": 1.63,
            "allen": 1.53,
        },
        "group": 5,
        "period": 4,
        "classification": "transition metal",
    },
    "Cr": {
        "atomic_number": 24,
        "name": "Chromium",
        "mass": 51.9961,  # in u
        "vdw_radius": 2.0,  # in Å
        "bond_params": {
            "cordero": {'single': 1.39},
            "pyykko": {'single': 1.22, 'double': 1.11, 'triple': 1.03},
        },
        "electronegativity": {
            "pauling": 1.66,
            "allen": 1.65,
        },
        "group": 6,
        "period": 4,
        "classification": "transition metal",
    },
    "Mn": {
        "atomic_number": 25,
        "name": "Manganese",
        "mass": 54.938044,  # in u
        "vdw_radius": 2.0,  # in Å
        "bond_params": {
            "cordero": {'single': 1.5},
            "pyykko": {'single': 1.19, 'double': 1.05, 'triple': 1.03},
        },
        "electronegativity": {
            "pauling": 1.55,
            "allen": 1.75,
        },
        "group": 7,
        "period": 4,
        "classification": "transition metal",
    },
    "Fe": {
        "atomic_number": 26,
        "name": "Iron",
        "mass": 55.845,  # in u
        "vdw_radius": 2.0,  # in Å
        "bond_params": {
            "cordero": {'single': 1.42},
            "pyykko": {'single': 1.16, 'double': 1.09, 'triple': 1.02},
        },
        "electronegativity": {
            "pauling": 1.83,
            "allen": 1.8,
        },
        "group": 8,
        "period": 4,
        "classification": "transition metal",
    },
    "Co": {
        "atomic_number": 27,
        "name": "Cobalt",
        "mass": 58.933194,  # in u
        "vdw_radius": 2.0,  # in Å
        "bond_params": {
            "cordero": {'single': 1.38},
            "pyykko": {'single': 1.11, 'double': 1.03, 'triple': 0.96},
        },
        "electronegativity": {
            "pauling": 1.88,
            "allen": 1.84,
        },
        "group": 9,
        "period": 4,
        "classification": "transition metal",
    },
    "Ni": {
        "atomic_number": 28,
        "name": "Nickel",
        "mass": 58.6934,  # in u
        "vdw_radius": 1.63,  # in Å
        "bond_params": {
            "cordero": {'single': 1.24},
            "pyykko": {'single': 1.1, 'double': 1.01, 'triple': 1.01},
        },
        "electronegativity": {
            "pauling": 1.91,
            "allen": 1.88,
        },
        "group": 10,
        "period": 4,
        "classification": "transition metal",
    },
    "Cu": {
        "atomic_number": 29,
        "name": "Copper",
        "mass": 63.546,  # in u
        "vdw_radius": 1.4,  # in Å
        "bond_params": {
            "cordero": {'single': 1.32},
            "pyykko": {'single': 1.12, 'double': 1.15, 'triple': 1.2},
        },
        "electronegativity": {
            "pauling": 1.9,
            "allen": 1.85,
        },
        "group": 11,
        "period": 4,
        "classification": "transition metal",
    },
    "Zn": {
        "atomic_number": 30,
        "name": "Zinc",
        "mass": 65.38,  # in u
        "vdw_radius": 1.39,  # in Å
        "bond_params": {
            "cordero": {'single': 1.22},
            "pyykko": {'single': 1.18, 'double': 1.2},
        },
        "electronegativity": {
            "pauling": 1.65,
            "allen": 1.588,
        },
        "group": 12,
        "period": 4,
        "classification": "transition metal",
    },
    "Ga": {
        "atomic_number": 31,
        "name": "Gallium",
        "mass": 69.723,  # in u
        "vdw_radius": 1.87,  # in Å
        "bond_params": {
            "cordero": {'single': 1.22},
            "pyykko": {'single': 1.24, 'double': 1.17, 'triple': 1.21},
        },
        "electronegativity": {
            "pauling": 1.81,
            "allen": 1.756,
        },
        "group": 13,
        "period": 4,
        "classification": "post-transition metal",
    },
    "Ge": {
        "atomic_number": 32,
        "name": "Germanium",
        "mass": 72.63,  # in u
        "vdw_radius": 2.11,  # in Å
        "bond_params": {
            "cordero": {'single': 1.2},
            "pyykko": {'single': 1.21, 'double': 1.11, 'triple': 1.14},
        },
        "electronegativity": {
            "pauling": 2.01,
            "allen": 1.994,
        },
        "group": 14,
        "period": 4,
        "classification": "metalloid",
    },
    "As": {
        "atomic_number": 33,
        "name": "Arsenic",
        "mass": 74.921595,  # in u
        "vdw_radius": 1.85,  # in Å
        "bond_params": {
            "cordero": {'single': 1.19},
            "pyykko": {'single': 1.21, 'double': 1.14, 'triple': 1.06},
        },
        "electronegativity": {
            "pauling": 2.18,
            "allen": 2.211,
        },
        "group": 15,
        "period": 4,
        "classification": "metalloid",
    },
    "Se": {
        "atomic_number": 34,
        "name": "Selenium",
        "mass": 78.971,  # in u
        "vdw_radius": 1.9,  # in Å
        "bond_params": {
            "cordero": {'single': 1.2},
            "pyykko": {'single': 1.16, 'double': 1.07, 'triple': 1.07},
        },
        "electronegativity": {
            "pauling": 2.55,
            "allen": 2.424,
        },
        "group": 16,
        "period": 4,
        "classification": "nonmetal",
    },
    "Br": {
        "atomic_number": 35,
        "name": "Bromine",
        "mass": 79.904,  # in u
        "vdw_radius": 1.85,  # in Å
        "bond_params": {
            "cordero": {'single': 1.2},
            "pyykko": {'single': 1.14, 'double': 1.09, 'triple': 1.1},
        },
        "electronegativity": {
            "pauling": 2.96,
            "allen": 2.685,
        },
        "group": 17,
        "period": 4,
        "classification": "halogen",
    },
    "Kr": {
        "atomic_number": 36,
        "name": "Krypton",
        "mass": 83.798,  # in u
        "vdw_radius": 2.02,  # in Å
        "bond_params": {
            "cordero": {'single': 1.16},
            "pyykko": {'single': 1.17, 'double': 1.21, 'triple': 1.08},
        },
        "electronegativity": {
            "pauling": 3.0,
            "allen": 2.966,
        },
        "group": 18,
        "period": 4,
        "classification": "noble gas",
    },
    "Rb": {
        "atomic_number": 37,
        "name": "Rubidium",
        "mass": 85.4678,  # in u
        "vdw_radius": 3.03,  # in Å
        "bond_params": {
            "cordero": {'single': 2.2},
            "pyykko": {'single': 2.1, 'double': 2.02},
        },
        "electronegativity": {
            "pauling": 0.82,
            "allen": 0.706,
        },
        "group": 1,
        "period": 5,
        "classification": "alkali metal",
    },
    "Sr": {
        "atomic_number": 38,
        "name": "Strontium",
        "mass": 87.62,  # in u
        "vdw_radius": 2.49,  # in Å
        "bond_params": {
            "cordero": {'single': 1.95},
            "pyykko": {'single': 1.85, 'double': 1.57, 'triple': 1.39},
        },
        "electronegativity": {
            "pauling": 0.95,
            "allen": 0.963,
        },
        "group": 2,
        "period": 5,
        "classification": "alkaline earth metal",
    },
    "Y": {
        "atomic_number": 39,
        "name": "Yttrium",
        "mass": 88.90584,  # in u
        "vdw_radius": 2.19,  # in Å
        "bond_params": {
            "cordero": {'single': 1.9},
            "pyykko": {'single': 1.63, 'double': 1.3, 'triple': 1.24},
        },
        "electronegativity": {
            "pauling": 1.22,
            "allen": 1.12,
        },
        "group": 3,
        "period": 5,
        "classification": "transition metal",
    },
    "Zr": {
        "atomic_number": 40,
        "name": "Zirconium",
        "mass": 91.224,  # in u
        "vdw_radius": 2.06,  # in Å
        "bond_params": {
            "cordero": {'single': 1.75},
            "pyykko": {'single': 1.54, 'double': 1.27, 'triple': 1.21},
        },
        "electronegativity": {
            "pauling": 1.33,
            "allen": 1.32,
        },
        "group": 4,
        "period": 5,
        "classification": "transition metal",
    },
    "Nb": {
        "atomic_number": 41,
        "name": "Niobium",
        "mass": 92.90637,  # in u
        "vdw_radius": 1.98,  # in Å
        "bond_params": {
            "cordero": {'single': 1.64},
            "pyykko": {'single': 1.47, 'double': 1.25, 'triple': 1.16},
        },
        "electronegativity": {
            "pauling": 1.6,
            "allen": 1.41,
        },
        "group": 5,
        "period": 5,
        "classification": "transition metal",
    },
    "Mo": {
        "atomic_number": 42,
        "name": "Molybdenum",
        "mass": 95.95,  # in u
        "vdw_radius": 1.9,  # in Å
        "bond_params": {
            "cordero": {'single': 1.54},
            "pyykko": {'single': 1.38, 'double': 1.21, 'triple': 1.13},
        },
        "electronegativity": {
            "pauling": 2.16,
            "allen": 1.47,
        },
        "group": 6,
        "period": 5,
        "classification": "transition metal",
    },
    "Tc": {
        "atomic_number": 43,
        "name": "Technetium",
        "mass": 98.0,  # in u
        "vdw_radius": 1.83,  # in Å
        "bond_params": {
            "cordero": {'single': 1.47},
            "pyykko": {'single': 1.28, 'double': 1.2, 'triple': 1.1},
        },
        "electronegativity": {
            "pauling": 1.9,
            "allen": 1.51,
        },
        "group": 7,
        "period": 5,
        "classification": "transition metal",
    },
    "Ru": {
        "atomic_number": 44,
        "name": "Ruthenium",
        "mass": 101.07,  # in u
        "vdw_radius": 1.78,  # in Å
        "bond_params": {
            "cordero": {'single': 1.46},
            "pyykko": {'single': 1.25, 'double': 1.14, 'triple': 1.03},
        },
        "electronegativity": {
            "pauling": 2.2,
            "allen": 1.54,
        },
        "group": 8,
        "period": 5,
        "classification": "transition metal",
    },
    "Rh": {
        "atomic_number": 45,
        "name": "Rhodium",
        "mass": 102.9055,  # in u
        "vdw_radius": 1.73,  # in Å
        "bond_params": {
            "cordero": {'single': 1.42},
            "pyykko": {'single': 1.25, 'double': 1.1, 'triple': 1.06},
        },
        "electronegativity": {
            "pauling": 2.28,
            "allen": 1.56,
        },
        "group": 9,
        "period": 5,
        "classification": "transition metal",
    },
    "Pd": {
        "atomic_number": 46,
        "name": "Palladium",
        "mass": 106.42,  # in u
        "vdw_radius": 1.63,  # in Å
        "bond_params": {
            "cordero": {'single': 1.39},
            "pyykko": {'single': 1.2, 'double': 1.17, 'triple': 1.12},
        },
        "electronegativity": {
            "pauling": 2.2,
            "allen": 1.58,
        },
        "group": 10,
        "period": 5,
        "classification": "transition metal",
    },
    "Ag": {
        "atomic_number": 47,
        "name": "Silver",
        "mass": 107.8682,  # in u
        "vdw_radius": 1.72,  # in Å
        "bond_params": {
            "cordero": {'single': 1.45},
            "pyykko": {'single': 1.28, 'double': 1.39, 'triple': 1.37},
        },
        "electronegativity": {
            "pauling": 1.93,
            "allen": 1.87,
        },
        "group": 11,
        "period": 5,
        "classification": "transition metal",
    },
    "Cd": {
        "atomic_number": 48,
        "name": "Cadmium",
        "mass": 112.414,  # in u
        "vdw_radius": 1.58,  # in Å
        "bond_params": {
            "cordero": {'single': 1.44},
            "pyykko": {'single': 1.36, 'double': 1.44},
        },
        "electronegativity": {
            "pauling": 1.69,
            "allen": 1.521,
        },
        "group": 12,
        "period": 5,
        "classification": "transition metal",
    },
    "In": {
        "atomic_number": 49,
        "name": "Indium",
        "mass": 114.818,  # in u
        "vdw_radius": 1.93,  # in Å
        "bond_params": {
            "cordero": {'single': 1.42},
            "pyykko": {'single': 1.42, 'double': 1.36, 'triple': 1.46},
        },
        "electronegativity": {
            "pauling": 1.78,
            "allen": 1.656,
        },
        "group": 13,
        "period": 5,
        "classification": "post-transition metal",
    },
    "Sn": {
        "atomic_number": 50,
        "name": "Tin",
        "mass": 118.71,  # in u
        "vdw_radius": 2.17,  # in Å
        "bond_params": {
            "cordero": {'single': 1.39},
            "pyykko": {'single': 1.4, 'double': 1.3, 'triple': 1.32},
        },
        "electronegativity": {
            "pauling": 1.96,
            "allen": 1.824,
        },
        "group": 14,
        "period": 5,
        "classification": "post-transition metal",
    },
    "Sb": {
        "atomic_number": 51,
        "name": "Antimony",
        "mass": 121.76,  # in u
        "vdw_radius": 2.06,  # in Å
        "bond_params": {
            "cordero": {'single': 1.39},
            "pyykko": {'single': 1.4, 'double': 1.33, 'triple': 1.27},
        },
        "electronegativity": {
            "pauling": 2.05,
            "allen": 1.984,
        },
        "group": 15,
        "period": 5,
        "classification": "metalloid",
    },
    "Te": {
        "atomic_number": 52,
        "name": "Tellurium",
        "mass": 127.6,  # in u
        "vdw_radius": 2.06,  # in Å
        "bond_params": {
            "cordero": {'single': 1.38},
            "pyykko": {'single': 1.36, 'double': 1.28, 'triple': 1.21},
        },
        "electronegativity": {
            "pauling": 2.1,
            "allen": 2.158,
        },
        "group": 16,
        "period": 5,
        "classification": "metalloid",
    },
    "I": {
        "atomic_number": 53,
        "name": "Iodine",
        "mass": 126.90447,  # in u
        "vdw_radius": 1.98,  # in Å
        "bond_params": {
            "cordero": {'single': 1.39},
            "pyykko": {'single': 1.33, 'double': 1.29, 'triple': 1.25},
        },
        "electronegativity": {
            "pauling": 2.66,
            "allen": 2.359,
        },
        "group": 17,
        "period": 5,
        "classification": "halogen",
    },
    "Xe": {
        "atomic_number": 54,
        "name": "Xenon",
        "mass": 131.293,  # in u
        "bond_params": {
            "cordero": {'single': 1.4},
            "pyykko": {'single': 1.31, 'double': 1.35, 'triple': 1.22},
        },
        "electronegativity": {
            "allen": 2.582,
        },
        "group": 18,
        "period": 5,
        "classification": "noble gas",
    },
    "Cs": {
        "atomic_number": 55,
        "name": "Cesium",
        "mass": 132.90545196,  # in u
        "vdw_radius": 3.43,  # in Å
        "bond_params": {
            "cordero": {'single': 2.44},
            "pyykko": {'single': 2.32, 'double': 2.09},
        },
        "electronegativity": {
            "pauling": 0.79,
            "allen": 0.659,
        },
        "group": 1,
        "period": 6,
        "classification": "alkali metal",
    },
    "Ba": {
        "atomic_number": 56,
        "name": "Barium",
        "mass": 137.327,  # in u
        "vdw_radius": 2.68,  # in Å
        "bond_params": {
            "cordero": {'single': 2.15},
            "pyykko": {'single': 1.96, 'double': 1.61, 'triple': 1.49},
        },
        "electronegativity": {
            "pauling": 0.89,
            "allen": 0.881,
        },
        "group": 2,
        "period": 6,
        "classification": "alkaline earth metal",
    },
    "La": {
        "atomic_number": 57,
        "name": 'Lanthanum',
        "mass": 138.905,
        "group": 3,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.4,
        "electronegativity": {
            "pauling": 1.1,
        },
        "bond_params": {
            "cordero": {'single': 2.07},
            "pyykko": {'single': 1.8, 'double': 1.39, 'triple': 1.39},
        },
    },
    "Ce": {
        "atomic_number": 58,
        "name": 'Cerium',
        "mass": 140.116,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.35,
        "electronegativity": {
            "pauling": 1.12,
        },
        "bond_params": {
            "cordero": {'single': 2.04},
            "pyykko": {'single': 1.63, 'double': 1.37, 'triple': 1.31},
        },
    },
    "Pr": {
        "atomic_number": 59,
        "name": 'Praseodymium',
        "mass": 140.908,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.39,
        "electronegativity": {
            "pauling": 1.13,
        },
        "bond_params": {
            "cordero": {'single': 2.03},
            "pyykko": {'single': 1.76, 'double': 1.38, 'triple': 1.28},
        },
    },
    "Nd": {
        "atomic_number": 60,
        "name": 'Neodymium',
        "mass": 144.24,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.29,
        "electronegativity": {
            "pauling": 1.14,
        },
    },
    "Pm": {
        "atomic_number": 61,
        "name": 'Promethium',
        "mass": 145.0,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.36,
        "electronegativity": {
            "pauling": 1.13,
        },
    },
    "Sm": {
        "atomic_number": 62,
        "name": 'Samarium',
        "mass": 150.36,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.29,
        "electronegativity": {
            "pauling": 1.17,
        },
    },
    "Eu": {
        "atomic_number": 63,
        "name": 'Europium',
        "mass": 151.964,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.33,
        "electronegativity": {
            "pauling": 1.2,
        },
    },
    "Gd": {
        "atomic_number": 64,
        "name": 'Gadolinium',
        "mass": 157.25,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.37,
        "electronegativity": {
            "pauling": 1.2,
        },
        "bond_params": {
            "cordero": {'single': 1.96},
            "pyykko": {'single': 1.69, 'double': 1.35, 'triple': 1.32},
        },
    },
    "Tb": {
        "atomic_number": 65,
        "name": 'Terbium',
        "mass": 158.925,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.21,
        "electronegativity": {
            "pauling": 1.2,
        },
    },
    "Dy": {
        "atomic_number": 66,
        "name": 'Dysprosium',
        "mass": 162.5,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.29,
        "electronegativity": {
            "pauling": 1.22,
        },
    },
    "Ho": {
        "atomic_number": 67,
        "name": 'Holmium',
        "mass": 164.93,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.16,
        "electronegativity": {
            "pauling": 1.23,
        },
    },
    "Er": {
        "atomic_number": 68,
        "name": 'Erbium',
        "mass": 167.259,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.35,
        "electronegativity": {
            "pauling": 1.24,
        },
    },
    "Tm": {
        "atomic_number": 69,
        "name": 'Thulium',
        "mass": 168.934,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.27,
        "electronegativity": {
            "pauling": 1.25,
        },
    },
    "Yb": {
        "atomic_number": 70,
        "name": 'Ytterbium',
        "mass": 173.054,
        "group": None,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.42,
        "electronegativity": {
            "pauling": 1.1,
        },
    },
    "Lu": {
        "atomic_number": 71,
        "name": 'Lutetium',
        "mass": 174.967,
        "group": 3,
        "period": 6,
        "classification": 'lanthanoid',
        "vdw_radius": 2.21,
        "electronegativity": {
            "pauling": 1.27,
            "allen": 1.09,
        },
        "bond_params": {
            "cordero": {'single': 1.87},
            "pyykko": {'single': 1.62, 'double': 1.31, 'triple': 1.31},
        },
    },
    "Hf": {
        "atomic_number": 72,
        "name": 'Hafnium',
        "mass": 178.49,
        "group": 4,
        "period": 6,
        "classification": 'transition metal',
        "electronegativity": {
            "pauling": 1.3,
            "allen": 1.16,
        },
        "bond_params": {
            "cordero": {'single': 1.75},
            "pyykko": {'single': 1.52, 'double': 1.28, 'triple': 1.22},
        },
    },
    "Ta": {
        "atomic_number": 73,
        "name": 'Tantalum',
        "mass": 180.948,
        "group": 5,
        "period": 6,
        "classification": 'transition metal',
        "electronegativity": {
            "pauling": 1.5,
            "allen": 1.34,
        },
        "bond_params": {
            "cordero": {'single': 1.7},
            "pyykko": {'single': 1.46, 'double': 1.26, 'triple': 1.19},
        },
    },
    "W": {
        "atomic_number": 74,
        "name": 'Tungsten',
        "mass": 183.84,
        "group": 6,
        "period": 6,
        "classification": 'transition metal',
        "electronegativity": {
            "pauling": 2.36,
            "allen": 1.47,
        },
        "bond_params": {
            "cordero": {'single': 1.62},
            "pyykko": {'single': 1.37, 'double': 1.2, 'triple': 1.15},
        },
    },
    "Re": {
        "atomic_number": 75,
        "name": 'Rhenium',
        "mass": 186.207,
        "group": 7,
        "period": 6,
        "classification": 'transition metal',
        "electronegativity": {
            "pauling": 1.9,
            "allen": 1.6,
        },
        "bond_params": {
            "cordero": {'single': 1.51},
            "pyykko": {'single': 1.31, 'double': 1.19, 'triple': 1.1},
        },
    },
    "Os": {
        "atomic_number": 76,
        "name": 'Osmium',
        "mass": 190.23,
        "group": 8,
        "period": 6,
        "classification": 'transition metal',
        "electronegativity": {
            "pauling": 2.2,
            "allen": 1.65,
        },
        "bond_params": {
            "cordero": {'single': 1.44},
            "pyykko": {'single': 1.29, 'double': 1.16, 'triple': 1.09},
        },
    },
    "Ir": {
        "atomic_number": 77,
        "name": 'Iridium',
        "mass": 192.217,
        "group": 9,
        "period": 6,
        "classification": 'transition metal',
        "electronegativity": {
            "pauling": 2.2,
            "allen": 1.68,
        },
        "bond_params": {
            "cordero": {'single': 1.41},
            "pyykko": {'single': 1.22, 'double': 1.15, 'triple': 1.07},
        },
    },
    "Pt": {
        "atomic_number": 78,
        "name": 'Platinum',
        "mass": 195.084,
        "group": 10,
        "period": 6,
        "classification": 'transition metal',
        "electronegativity": {
            "pauling": 2.28,
            "allen": 1.72,
        },
        "bond_params": {
            "cordero": {'single': 1.36},
            "pyykko": {'single': 1.23, 'double': 1.12, 'triple': 1.1},
        },
    },
    "Au": {
        "atomic_number": 79,
        "name": 'Gold',
        "mass": 196.967,
        "group": 11,
        "period": 6,
        "classification": 'transition metal',
        "electronegativity": {
            "pauling": 2.54,
            "allen": 1.92,
        },
        "bond_params": {
            "cordero": {'single': 1.36},
            "pyykko": {'single': 1.24, 'double': 1.21, 'triple': 1.23},
        },
    },
    "Hg": {
        "atomic_number": 80,
        "name": 'Mercury',
        "mass": 200.59,
        "group": 12,
        "period": 6,
        "classification": 'transition metal',
        "electronegativity": {
            "pauling": 2.0,
            "allen": 1.765,
        },
    },
    "Tl": {
        "atomic_number": 81,
        "name": 'Thallium',
        "mass": 204.383,
        "group": 13,
        "period": 6,
        "classification": 'post-transition metal',
        "electronegativity": {
            "pauling": 1.62,
            "allen": 1.789,
        },
        "bond_params": {
            "cordero": {'single': 1.45},
            "pyykko": {'single': 1.44, 'double': 1.42, 'triple': 1.5},
        },
    },
    "Pb": {
        "atomic_number": 82,
        "name": 'Lead',
        "mass": 207.2,
        "group": 14,
        "period": 6,
        "classification": 'post-transition metal',
        "electronegativity": {
            "pauling": 2.33,
            "allen": 1.854,
        },
        "bond_params": {
            "cordero": {'single': 1.46},
            "pyykko": {'single': 1.44, 'double': 1.35, 'triple': 1.37},
        },
    },
    "Bi": {
        "atomic_number": 83,
        "name": 'Bismuth',
        "mass": 208.98,
        "group": 15,
        "period": 6,
        "classification": 'post-transition metal',
        "electronegativity": {
            "pauling": 2.02,
            "allen": 2.01,
        },
        "bond_params": {
            "cordero": {'single': 1.48},
            "pyykko": {'single': 1.51, 'double': 1.41, 'triple': 1.35},
        },
    },
    "Po": {
        "atomic_number": 84,
        "name": 'Polonium',
        "mass": 209.0,
        "group": 16,
        "period": 6,
        "classification": 'post-transition metal',
        "electronegativity": {
            "pauling": 2.0,
            "allen": 2.19,
        },
        "bond_params": {
            "cordero": {'single': 1.4},
            "pyykko": {'single': 1.45, 'double': 1.35, 'triple': 1.29},
        },
    },
    "At": {
        "atomic_number": 85,
        "name": 'Astatine',
        "mass": 210.0,
        "group": 17,
        "period": 6,
        "classification": 'halogen',
        "electronegativity": {
            "pauling": 2.2,
            "allen": 2.39,
        },
        "bond_params": {
            "cordero": {'single': 1.5},
            "pyykko": {'single': 1.47, 'double': 1.38, 'triple': 1.38},
        },
    },
    "Rn": {
        "atomic_number": 86,
        "name": 'Radon',
        "mass": 222.0,
        "group": 18,
        "period": 6,
        "classification": 'noble gas',
        "electronegativity": {
            "allen": 2.6,
        },
        "bond_params": {
            "cordero": {'single': 1.5},
            "pyykko": {'single': 1.42, 'double': 1.45, 'triple': 1.33},
        },
    },
    "Fr": {
        "atomic_number": 87,
        "name": 'Francium',
        "mass": 223.0,
        "group": 1,
        "period": 7,
        "classification": 'alkali metal',
        "electronegativity": {
            "pauling": 0.7,
            "allen": 0.67,
        },
    },
    "Ra": {
        "atomic_number": 88,
        "name": 'Radium',
        "mass": 226.0,
        "group": 2,
        "period": 7,
        "classification": 'alkaline earth metal',
        "electronegativity": {
            "pauling": 0.9,
            "allen": 0.89,
        },
        "bond_params": {
            "cordero": {'single': 2.21},
            "pyykko": {'single': 2.01, 'double': 1.73, 'triple': 1.59},
        },
    },
    "Ac": {
        "atomic_number": 89,
        "name": 'Actinium',
        "mass": 227.0,
        "group": 3,
        "period": 7,
        "classification": 'actinoid',
        "electronegativity": {
            "pauling": 1.1,
        },
        "bond_params": {
            "cordero": {'single': 2.15},
            "pyykko": {'single': 1.86, 'double': 1.53, 'triple': 1.4},
        },
    },
    "Th": {
        "atomic_number": 90,
        "name": 'Thorium',
        "mass": 232.038,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
        "electronegativity": {
            "pauling": 1.3,
        },
        "bond_params": {
            "cordero": {'single': 2.06},
            "pyykko": {'single': 1.75, 'double': 1.43, 'triple': 1.36},
        },
    },
    "Pa": {
        "atomic_number": 91,
        "name": 'Protactinium',
        "mass": 231.036,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
        "electronegativity": {
            "pauling": 1.5,
        },
        "bond_params": {
            "cordero": {'single': 2.0},
            "pyykko": {'single': 1.69, 'double': 1.38, 'triple': 1.29},
        },
    },
    "U": {
        "atomic_number": 92,
        "name": 'Uranium',
        "mass": 238.029,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
        "electronegativity": {
            "pauling": 1.38,
        },
        "bond_params": {
            "cordero": {'single': 1.96},
            "pyykko": {'single': 1.7, 'double': 1.34, 'triple': 1.18},
        },
    },
    "Np": {
        "atomic_number": 93,
        "name": 'Neptunium',
        "mass": 237.0,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
        "electronegativity": {
            "pauling": 1.36,
        },
        "bond_params": {
            "cordero": {'single': 1.9},
            "pyykko": {'single': 1.71, 'double': 1.36, 'triple': 1.16},
        },
    },
    "Pu": {
        "atomic_number": 94,
        "name": 'Plutonium',
        "mass": 244.0,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
        "electronegativity": {
            "pauling": 1.28,
        },
    },
    "Am": {
        "atomic_number": 95,
        "name": 'Americium',
        "mass": 243.0,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
    },
    "Cm": {
        "atomic_number": 96,
        "name": 'Curium',
        "mass": 247.0,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
    },
    "Bk": {
        "atomic_number": 97,
        "name": 'Berkelium',
        "mass": 247.0,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
    },
    "Cf": {
        "atomic_number": 98,
        "name": 'Californium',
        "mass": 251.0,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
    },
    "Es": {
        "atomic_number": 99,
        "name": 'Einsteinium',
        "mass": 252.0,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
    },
    "Fm": {
        "atomic_number": 100,
        "name": 'Fermium',
        "mass": 257.0,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
    },
    "Md": {
        "atomic_number": 101,
        "name": 'Mendelevium',
        "mass": 258.0,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
    },
    "No": {
        "atomic_number": 102,
        "name": 'Nobelium',
        "mass": 259.0,
        "group": None,
        "period": 7,
        "classification": 'actinoid',
    },
    "Lr": {
        "atomic_number": 103,
        "name": 'Lawrencium',
        "mass": 262.0,
        "group": 3,
        "period": 7,
        "classification": 'actinoid',
    },
    "Rf": {
        "atomic_number": 104,
        "name": 'Rutherfordium',
        "mass": 267.0,
        "group": 4,
        "period": 7,
        "classification": 'transition metal',
    },
    "Db": {
        "atomic_number": 105,
        "name": 'Dubnium',
        "mass": 268.0,
        "group": 5,
        "period": 7,
        "classification": 'transition metal',
    },
    "Sg": {
        "atomic_number": 106,
        "name": 'Seaborgium',
        "mass": 269.0,
        "group": 6,
        "period": 7,
        "classification": 'transition metal',
    },
    "Bh": {
        "atomic_number": 107,
        "name": 'Bohrium',
        "mass": 270.0,
        "group": 7,
        "period": 7,
        "classification": 'transition metal',
    },
    "Hs": {
        "atomic_number": 108,
        "name": 'Hassium',
        "mass": 269.0,
        "group": 8,
        "period": 7,
        "classification": 'transition metal',
    },
    "Mt": {
        "atomic_number": 109,
        "name": 'Meitnerium',
        "mass": 278.0,
        "group": 9,
        "period": 7,
        "classification": 'transition metal',
    },
    "Ds": {
        "atomic_number": 110,
        "name": 'Darmstadtium',
        "mass": 281.0,
        "group": 10,
        "period": 7,
        "classification": 'transition metal',
    },
    "Rg": {
        "atomic_number": 111,
        "name": 'Roentgenium',
        "mass": 282.0,
        "group": 11,
        "period": 7,
        "classification": 'transition metal',
    },
    "Cn": {
        "atomic_number": 112,
        "name": 'Copernicium',
        "mass": 285.0,
        "group": 12,
        "period": 7,
        "classification": 'transition metal',
    },
    "Nh": {
        "atomic_number": 113,
        "name": 'Nihonium',
        "mass": 286.0,
        "group": 13,
        "period": 7,
        "classification": 'post-transition metal',
    },
    "Fl": {
        "atomic_number": 114,
        "name": 'Flerovium',
        "mass": 289.0,
        "group": 14,
        "period": 7,
        "classification": 'post-transition metal',
    },
    "Mc": {
        "atomic_number": 115,
        "name": 'Moscovium',
        "mass": 290.0,
        "group": 15,
        "period": 7,
        "classification": 'post-transition metal',
    },
    "Lv": {
        "atomic_number": 116,
        "name": 'Livermorium',
        "mass": 293.0,
        "group": 16,
        "period": 7,
        "classification": 'post-transition metal',
    },
    "Ts": {
        "atomic_number": 117,
        "name": 'Tennessine',
        "mass": 294.0,
        "group": 17,
        "period": 7,
        "classification": 'halogen',
    },
    "Og": {
        "atomic_number": 118,
        "name": 'Oganesson',
        "mass": 294.0,
        "group": 18,
        "period": 7,
        "classification": 'noble gas',
        },

}

# #############################################################################
# 2) Unit normalization helpers
# #############################################################################

# Distances are stored in Å internally. For example, to convert from Å to pm,
# multiply by 100.

_DISTANCE_UNIT_ALIASES = {
    # Existing synonyms for Å
    "a": ("Ang", 1.0),
    "ang": ("Ang", 1.0),
    "angstrom": ("Ang", 1.0),
    "angstroms": ("Ang", 1.0),
    "ångström": ("Ang", 1.0),
    "å": ("Ang", 1.0),
    "Å": ("Ang", 1.0),

    # Picometers
    "pm": ("pm", 100.0),
    "picometer": ("pm", 100.0),
    "picometre": ("pm", 100.0),

    # Nanometers
    "nm": ("nm", 0.1),
    "nanometer": ("nm", 0.1),
    "nanometre": ("nm", 0.1),

    # Bohr (atomic units)
    "bohr": ("bohr", 1.889725989),
    "a0": ("bohr", 1.889725989),
    "au": ("bohr", 1.889725989),
    "atomic_unit": ("bohr", 1.889725989),

    # (Optional) Femtometers, micrometers, meters, etc.
    # "fm": ("fm", 1.0e5),  # 1 Å = 1e5 fm
    # "um": ("um", 1e-4),   # 1 Å = 1e-4 µm
    # "µm": ("um", 1e-4),
    # "m": ("m", 1e-10),    # 1 Å = 1e-10 m
    # "meter": ("m", 1e-10),
}
# Mass is stored in atomic mass units (u). 1 u ≈ 1 g/mol
_MASS_UNIT_ALIASES = {
    "u": ("u", 1.0),
    "amu": ("u", 1.0),
    "g/mol": ("g/mol", 1.0),
    "grams/mol": ("g/mol", 1.0),
}


def _normalize_symbol(symbol: str) -> str:
    """Normalize element symbol input: strip, capitalize first letter, check validity."""
    s = symbol.strip().capitalize()
    if s not in _ELEMENTS:
        raise KeyError(f"Unknown element symbol: '{symbol}'")
    return s


def _normalize_distance_unit(unit: str):
    """Returns (canonical_distance_unit_name, factor_from_Ang) for given distance unit."""
    key = unit.strip().lower()
    if key not in _DISTANCE_UNIT_ALIASES:
        raise KeyError(f"Unknown distance unit: '{unit}'")
    return _DISTANCE_UNIT_ALIASES[key]


def _normalize_mass_unit(unit: str):
    """Returns (canonical_mass_unit_name, factor_from_u) for given mass unit."""
    key = unit.strip().lower()
    if key not in _MASS_UNIT_ALIASES:
        raise KeyError(f"Unknown mass unit: '{unit}'")
    return _MASS_UNIT_ALIASES[key]


def _get_symbol_by_atomic_number(atomic_num: int) -> str:
    """Reverse lookup: atomic number -> symbol. Raises KeyError if not found."""
    for sym, data in _ELEMENTS.items():
        if data["atomic_number"] == atomic_num:
            return sym
    raise KeyError(f"No element found for atomic number {atomic_num}")


# #############################################################################
# 3) The class with @classmethod lookups
# #############################################################################

class Elements:
    @classmethod
    def is_valid(cls, symbol: str) -> bool:
        """Return True if `symbol` is a known element (case-insensitive)."""
        try:
            _normalize_symbol(symbol)
            return True
        except KeyError:
            return False

    @classmethod
    def list_symbols(cls):
        """Return a list of all element symbols in the internal dictionary."""
        return sorted(_ELEMENTS.keys(), key=lambda s: _ELEMENTS[s]["atomic_number"])

    @classmethod
    def atomic_number(cls, symbol: str) -> int:
        """Return atomic number of the given element symbol."""
        sym = _normalize_symbol(symbol)
        return _ELEMENTS[sym]["atomic_number"]

    @classmethod
    def symbol(cls, atomic_num: int) -> str:
        """Return element symbol for the given atomic number."""
        return _get_symbol_by_atomic_number(atomic_num)

    @classmethod
    def name(cls, symbol: str) -> str:
        """Return the full name of the element."""
        sym = _normalize_symbol(symbol)
        return _ELEMENTS[sym]["name"]

    @classmethod
    def mass(cls, symbol: str, unit: str = "u") -> float:
        """
        Return the atomic mass of the element in the specified unit.
        Internally stored in 'u'. (1 u = 1 g/mol)
        """
        sym = _normalize_symbol(symbol)
        mass_u = _ELEMENTS[sym]["mass"]
        _, factor = _normalize_mass_unit(unit)
        return mass_u * factor

    @classmethod
    def vdw_radius(cls, symbol: str, unit: str = "Ang") -> float:
        """
        Return the van der Waals radius for the given element symbol in
        the specified distance unit (e.g. 'Ang' or 'pm').
        """
        sym = _normalize_symbol(symbol)
        val = _ELEMENTS[sym]["vdw_radius"]
        if val is None:
            raise KeyError(f"No van der Waals radius available for '{sym}'")
        _, factor = _normalize_distance_unit(unit)
        return val * factor

    @classmethod
    def covalent_radius(cls,
                        symbol: str,
                        order: str = "single",
                        source: str = "cordero",
                        unit: str = "Ang") -> float:
        """
        Return the covalent radius from the specified source ('cordero' or 'pyykko')
        for the requested bond order ('single', 'double', 'triple') in the specified
        distance unit.
        """
        sym = _normalize_symbol(symbol)
        bond_data = _ELEMENTS[sym].get("bond_params", {})
        src_data = bond_data.get(source)
        if src_data is None:
            raise KeyError(f"No covalent radius data for source='{source}' in '{sym}'")

        if order not in src_data:
            raise KeyError(f"No covalent radius for bond order='{order}' in source='{source}' for '{sym}'")

        radius_A = src_data[order]  # in Å
        _, factor = _normalize_distance_unit(unit)
        return radius_A * factor

    @classmethod
    def electronegativity(cls, symbol: str, scale: str = "pauling") -> float:
        """
        Return the electronegativity of an element on the specified scale
        (pauling or allen). Raises KeyError if no data is found.
        """
        sym = _normalize_symbol(symbol)
        e_data = _ELEMENTS[sym].get("electronegativity", {})
        if scale not in e_data:
            raise KeyError(f"No electronegativity found for scale='{scale}' in '{sym}'")
        val = e_data[scale]
        if val is None:
            raise KeyError(f"No electronegativity value on '{scale}' scale for '{sym}' (data is None)")
        return val

    @classmethod
    def period(cls, symbol: str) -> int:
        """Return the period number for the given element symbol."""
        sym = _normalize_symbol(symbol)
        return _ELEMENTS[sym]["period"]

    @classmethod
    def group(cls, symbol: str) -> int:
        """Return the group number for the given element symbol."""
        sym = _normalize_symbol(symbol)
        return _ELEMENTS[sym]["group"]

    @classmethod
    def classification(cls, symbol: str) -> str:
        """
        Return the classification of the element (e.g. 'nonmetal', 'alkali metal',
        'noble gas', 'transition metal', etc.).
        """
        sym = _normalize_symbol(symbol)
        return _ELEMENTS[sym]["classification"]




# -----------------------------------------------------------------------------
# 4) Simple test suite using unittest
# -----------------------------------------------------------------------------
import unittest

class TestElements(unittest.TestCase):
    def test_is_valid(self):
        self.assertTrue(Elements.is_valid("H"))
        self.assertTrue(Elements.is_valid("he"))
        self.assertFalse(Elements.is_valid("foo"))

    def test_atomic_number(self):
        self.assertEqual(Elements.atomic_number("H"), 1)
        self.assertEqual(Elements.atomic_number("he"), 2)
        with self.assertRaises(KeyError):
            Elements.atomic_number("fakeSymbol")

    def test_symbol(self):
        self.assertEqual(Elements.symbol(1), "H")
        self.assertEqual(Elements.symbol(2), "He")
        with self.assertRaises(KeyError):
            Elements.symbol(9999)  # Non-existent atomic number

    def test_name(self):
        self.assertEqual(Elements.name("H"), "Hydrogen")
        self.assertEqual(Elements.name("hE"), "Helium")

    def test_mass(self):
        # Checking default units (u)
        self.assertAlmostEqual(Elements.mass("H"), 1.008, places=3)
        # Checking g/mol
        self.assertAlmostEqual(Elements.mass("H", "g/mol"), 1.008, places=3)
        # Unknown unit
        with self.assertRaises(KeyError):
            Elements.mass("H", "pounds")

    def test_vdw_radius(self):
        # Checking default units (Ang)
        self.assertAlmostEqual(Elements.vdw_radius("H"), 1.2, places=2)
        # Checking pm
        self.assertAlmostEqual(Elements.vdw_radius("H", "pm"), 120.0, places=1)
        # Element without VDW radius -> KeyError (if missing data in dictionary)

    def test_covalent_radius(self):
        # single bond, cordero for H
        self.assertAlmostEqual(Elements.covalent_radius("H"), 0.31, places=2)
        # single bond, pyykko for H
        self.assertAlmostEqual(Elements.covalent_radius("H", source="pyykko"), 0.32, places=2)
        # bond order not found
        with self.assertRaises(KeyError):
            Elements.covalent_radius("H", order="quadruple")

    def test_electronegativity(self):
        self.assertAlmostEqual(Elements.electronegativity("H", "pauling"), 2.2, places=2)
        with self.assertRaises(KeyError):
            Elements.electronegativity("H", "bogusScale")

    def test_period_and_group(self):
        self.assertEqual(Elements.period("Li"), 2)
        self.assertEqual(Elements.group("Li"), 1)
        self.assertEqual(Elements.period("He"), 1)
        self.assertEqual(Elements.group("He"), 18)

    def test_classification(self):
        self.assertEqual(Elements.classification("H"), "nonmetal")
        self.assertEqual(Elements.classification("he"), "noble gas")

if __name__ == "__main__":
    unittest.main()

