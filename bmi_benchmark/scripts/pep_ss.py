#!/usr/bin/env python
"""pep_ss.py -- coordinate-only secondary-structure assignment for peptides.

PER-RESIDUE ALPHABET (one character per residue of the peptide, in backbone order)
---------------------------------------------------------------------------------
    H   right-handed alpha-helix   -- two consecutive i->i+4 backbone H-bonds, positive twist
    h   LEFT-handed alpha-helix    -- same H-bond pattern, negative twist (D-peptides)
    G   right-handed 3-10 helix    -- two consecutive i->i+3 H-bonds, positive twist
    g   LEFT-handed 3-10 helix
    I   pi-helix                   -- two consecutive i->i+5 H-bonds (handedness reported
                                      separately; pi is rare enough that it is not split)
    E   extended and BRIDGED TO THE RECEPTOR backbone      <- receptor-induced by construction
    e   extended and BRIDGED TO ANOTHER PEPTIDE SEGMENT    <- intrinsic beta (hairpin, sheet)
    x   extended phi/psi, no bridge partner at all
    P   polyproline II, left-handed (the L-amino-acid form)
    p   polyproline II, right-handed (its mirror image, the D-amino-acid form)
    T   beta-turn central residue (position i+1 or i+2 of a typed four-residue turn)
    S   gamma-turn central residue (position i+1 of a typed three-residue turn)
    C   coil -- backbone complete and alpha, nothing above matched
    X   NOT ASSIGNABLE -- N, CA, C or O missing from the deposited model
    ?   REFUSED -- the backbone is not an alpha-amino-acid backbone (beta-amino acid,
        gamma-amino acid, non-amino-acid linker).  phi/psi and chirality are null for these
        residues and the reason is recorded.  A confident wrong answer is worse than a null.

Precedence, applied in this order:  helix (H h G g I) > bridge (E e) > turn (T S) >
PPII (P p) > extended-unpaired (x) > coil (C).  X and ? pre-empt everything.

X and ? residues are excluded from every fraction denominator; the denominator is
`n_assignable`, reported alongside every fraction.

THRESHOLDS.  Every number this module decides on is a module constant in the THRESH dict
below, each with the measurement or the citation that fixes it.  Nothing is hardcoded inline.

WHAT THIS TOOL DOES NOT DO.  Short version: it does not
predict what the free peptide would do, it does not assign anything to a residue whose
backbone is not an alpha-amino-acid backbone, it does not use side-chain H-bonds for
secondary structure, it cannot see a symmetry mate, and it cannot recover a residue that was
never modelled.

USAGE
    pep_ss.py --selftest                     ideal-geometry calibration, exits non-zero on failure
    pep_ss.py --entry 1HC9 [--json out.json] one entry, human-readable report + optional JSON
    pep_ss.py --entry 1HC9 --brief           one line
    pep_ss.py --refs                         print the derived reference table and exit
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

try:
    import gemmi
except ImportError:  # pragma: no cover
    gemmi = None
from scipy.spatial import cKDTree

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
VERSION = "pep_ss-1.1"
# 1.1, 2026-08-03: _bridges' two WIDE bridge clauses were the mirror image of Kabsch
#   & Sander -- a (donor, acceptor) versus K&S (acceptor, donor) convention inversion,
#   one cause, both clauses.  See the dated note inside _bridges and the new
#   _bridge_selftest.  Rows written by 1.0 carry the pre-correction beta assignment;
#   helix, turn, PPII, torsion, chirality, topology and H-bond columns are unaffected.
# 1.0: as delivered 2026-08-03.

# --------------------------------------------------------------------------- #
# THRESHOLDS.  Each entry carries the source that fixes it.  The measurements
# named "measured here" were made on this reference set by the calibration scripts; they are
# not quoted from anywhere.
# --------------------------------------------------------------------------- #
THRESH = {
    # ---- covalent connectivity -------------------------------------------- #
    # Amide C(i)-N(i+1) link.  Two measurements fix this, and the SECOND one is
    # the binding one because the peptides are the hard case.
    #   (a) standard-residue peptide bonds measured across a reference protein set:
    #       mean 1.332, sd 0.009, max 1.411 A; over the same chains the closest
    #       NON-bonded C...N pair anywhere is >= 2.6 A, zero pairs below 2.6.
    #   (b) EVERY C...N pair below 2.6 A across the peptide.cif files: a dense
    #       population from 1.15 to 1.75 A (99.90 %),
    #       then a sparse tail of 11 pairs spread over 1.75-2.35 A, then the
    #       non-bonded population resuming at 2.35 A and above (19 pairs).
    # 1.75 is the top of the dense population.  It matters: 4B8Y's head-to-tail
    # closure is a genuine amide modelled at 1.534 A and a 1.50 cutoff misses it,
    # turning a macrocycle into a linear chain.  Peptides in this reference set include
    # low-resolution and heavily modified components, and wwPDB checks no bond
    # geometry at all for non-standard components, so their amides are refined
    # more loosely than a receptor's -- which is exactly why (a) alone would have
    # been the wrong basis for this number.
    # Links outside [1.20, 1.45] are reported in `unusual_amide_links` so a
    # stretched or physically impossible one is visible rather than absorbed.
    "link_CN_max_A": 1.75,
    "link_CN_min_A": 1.05,
    "link_CN_normal_A": (1.20, 1.45),
    # Heavy-atom bond, generic: r_cov(i) + r_cov(j) + this.  Used for the
    # intra-residue bond graph (backbone-path length) only.  Covalent radii from
    # gemmi's element table.  0.35 puts the C-C cutoff at 1.87 and the C-N cutoff
    # at 1.82, above every real bond and below the shortest 1-3 contact (~2.4).
    "bond_slack_A": 0.35,
    # Bond from a backbone N to any heavy atom.  Measured over a set of
    # modified entries: distances from a peptide backbone N to any heavy atom
    # in the complex fall in 1.15-1.75 A (2179 pairs) and then in 2.05-2.8 A
    # (4900 pairs) with an EMPTY interval 1.75-2.05.  1.80 is the middle of that
    # gap.  Neighbours found in 1.70-2.00 are additionally flagged as borderline.
    "n_bond_max_A": 1.80,
    "n_bond_borderline_A": (1.70, 2.00),
    # ---- what counts as an alpha-amino-acid backbone ---------------------- #
    # Measured over a set of standard residues: N-CA 1.461 +- 0.013 (range
    # 1.387-1.540), CA-C 1.524 +- 0.012 (range 1.450-1.649).  The windows below
    # are ~+-10 sd wide, i.e. they will not reject a badly refined but genuine
    # alpha residue, and they still reject every real non-alpha case in the
    # reference set by a wide margin: 4BPJ/4BPI HR7 has N-CA 2.44 (gamma-amino acid
    # numbering), 4BPJ B3A/B3D/B3Q have CA-C > 1.9 (beta3 backbone), 7Z4S II7
    # has CA-C 4.04.
    "n_ca_range_A": (1.30, 1.65),
    "ca_c_range_A": (1.35, 1.75),
    # An alpha backbone has exactly two bonds on the covalent path N -> C
    # (N-CA-C).  A beta3 residue has three.  This is the PRIMARY test; the
    # distance windows above are the confirmation.
    "backbone_path_bonds": 2,
    # ---- hydrogen bonds --------------------------------------------------- #
    # Calibrated on i->i-4 backbone H-bonds inside phi/psi-defined helices
    # of >= 7 residues in a reference protein set:
    #     N...O  mean 3.017, sd 0.251, p95 3.386;  96.3 % <= 3.5 A
    #     N-H...O angle  mean 157.2, sd 10.0, p1 125.6;  99.5 % >= 120 deg
    # 3.5 A / 120 deg therefore retains ~96 % of genuine helical H-bonds.  The
    # 3.7 % it drops are helix termini and kinks where N...O runs to 4-5.6 A;
    # those are real distortions, not detection failures, and calling them
    # non-bonded shortens a helix by at most one residue at each end.
    "hb_NO_max_A": 3.50,
    "hb_NHO_min_deg": 120.0,
    # N-H bond length used to place the inferred hydrogen.  Neutron/X-ray amide
    # N-H is 1.01 A; the value only affects the angle, and the angle criterion is
    # insensitive to it (changing 1.01 -> 0.90 moves the mean N-H...O angle by
    # < 1 deg over the calibration set).
    "nh_length_A": 1.01,
    # DSSP's own Kabsch-Sander cutoff, reported for comparison, never decisive.
    "ks_energy_max_kcal": -0.50,
    # ---- helix / strand runs ---------------------------------------------- #
    # A helix needs two consecutive n-turns, which already implies >= 4 residues
    # for n=4.  These minima are for the AGGREGATE category only, never for the
    # per-residue string.  4: one full alpha turn (3.6 residues/turn).  2: the
    # shortest bridged stretch that can pair along its length; a single bridge is
    # one H-bond pair, reported as such but not called a strand.
    "helix_min_run": 4,
    "strand_min_run": 2,
    # PPII has 3.0 residues per turn, so 3 is one full turn.  Same value as
    # mkdssp's --min-pp-stretch default, so the two PPII opinions are comparable.
    "ppii_min_run": 3,
    # ---- turns ------------------------------------------------------------ #
    # Four-residue turn gate.  Derived, not quoted: the ideal builds below give
    # CA(i)-CA(i+3) = 5.0 (type I), 5.5 (type II), 5.0 (I'), 5.5 (II'), 6.1
    # (VIII); an ideal PPII trimer gives 9.4 and an ideal antiparallel strand
    # 10.0.  7.0 A separates the two populations with >= 0.9 A of margin on the
    # turn side and >= 2.4 A on the extended side.  This is also the value used
    # by the standard turn definition (Lewis/Richardson), which is why it is
    # familiar; here it is confirmed against our own ideal geometry.
    "turn_ca13_max_A": 7.00,
    # Tolerance on each of the four turn torsions.  +-30 deg is the Wilmot &
    # Thornton / Hutchinson & Thornton convention and CANNOT be derived from
    # geometry -- it is a convention about how much distortion still counts.  It
    # is applied strictly here (all four angles within 30 deg).  A second, looser
    # pass allows ONE of the four angles up to 45 deg, matching the other common
    # form of the same convention; both results are reported, never merged.
    "turn_tol_deg": 30.0,
    "turn_tol_loose_one_deg": 45.0,
    "gamma_tol_deg": 40.0,
    # ---- chirality -------------------------------------------------------- #
    # improper = dihedral(N, C, CB, CA).  Sign convention calibrated against real
    # coordinates, not asserted: over a set of standard L-residues the
    # improper is -33.6 +- 2.0 deg with zero positive
    # values, and 2WFJ DAL 1 (D-alanine of cyclosporin A, the only D residue in
    # that peptide) gives +33.3.  |improper| below the cutoff is a planar or
    # badly modelled centre and is refused.
    "chiral_L_sign": -1.0,
    "chiral_min_abs_deg": 10.0,
    # ---- geometry signals ------------------------------------------------- #
    # All of the following are DERIVED in refs() from ideal builds, not set here.
    # The only free choices are the acceptance half-widths, stated in refs().
    "axis_window": 5,          # CA positions per local helix-axis fit; 4 is the
                               # minimum that determines an axis, 5 gives one
                               # redundant point for the circle fit residual
    "axis_max_twist_deg": 150.0,  # above this the projected points collapse onto
                               # two clusters and the circle fit is meaningless
                               # (a beta strand is a 2-residue/turn "helix");
                               # return null rather than a number.
    "axis_max_circle_rms_A": 0.60,   # circle-fit residual above which the window
                                     # is not a regular helix at all
    # receptor residues within this of any peptide atom take part in bridge and
    # H-bond detection.  The H-bond criterion itself is 3.5 A; 8.0 leaves room
    # for the i-1/i+1 partners a bridge definition needs.
    "receptor_shell_A": 8.0,
    # CA-CA pre-filter on candidate bridge partners.  Measured over a set of
    # H-bonded backbone pairs: the CA-CA distance is
    # 5.65 +- , p99 6.67, max 6.79 A.  8.0 therefore cannot exclude a real
    # bridge; it only keeps the pair enumeration small.  The two-H-bond bridge
    # clause, not this number, decides whether a pair is a bridge.
    "bridge_ca_max_A": 8.0,
}

BB = ("N", "CA", "C", "O")

# --------------------------------------------------------------------------- #
# Ideal geometry.  Bond lengths and angles MEASURED on this reference set,
# not quoted:
#   N-CA  1.461 +- 0.013     N-CA-C   111.06 +- 2.28
#   CA-C  1.524 +- 0.012     CA-C-N   116.80 +- 1.17
#   C-N   1.332 +- 0.009     C-N-CA   121.39 +- 1.47
#   C-O   1.233 +- 0.008     CA-C-O   120.46 +- 0.91
# --------------------------------------------------------------------------- #
IDEAL = {
    "N_CA": 1.461, "CA_C": 1.524, "C_N": 1.332, "C_O": 1.233,
    "N_CA_C": 111.06, "CA_C_N": 116.80, "C_N_CA": 121.39, "CA_C_O": 120.46,
    # side-chain attachment at CA, also measured here, on a set of standard
    # non-Gly non-Pro residues:
    #   CA-CB 1.532 +- 0.012,  N-CA-CB 110.69 +- 1.25,
    #   dihedral(C, N, CA, CB) = -122.59 +- 2.24 for an L residue.
    # The last number is what makes the chirality selftest a real test: the ideal
    # L chain is BUILT from it and the assigner must then read L back out.
    "CA_CB": 1.532, "N_CA_CB": 110.69, "CB_DIH_L": -122.59,
}

# Reference (phi, psi) of the named conformations.  THESE ARE THE ONLY NUMBERS IN
# THIS MODULE TAKEN FROM THE LITERATURE, and they are inputs, not thresholds:
# every threshold derived from them is computed in refs() from the ideal chain
# they build.  Sources, all textbook:
#   alpha_R      (-57, -47)   Pauling & Corey 1951 alpha-helix
#   alpha_L      (+57, +47)   its exact mirror image
#   helix310     (-49, -26)   3-10 helix
#   helix_pi     (-55, -70)   pi-helix
#   beta_anti    (-139, 135)  antiparallel beta-sheet
#   beta_par     (-119, 113)  parallel beta-sheet
#   ppii         (-75, 145)   polyproline II
#   ppii_mirror  (+75, -145)  its exact mirror image
REF_TORSIONS = {
    "alpha_R":     (-57.0, -47.0),
    "alpha_L":     (57.0, 47.0),
    "helix310":    (-49.0, -26.0),
    "helix_pi":    (-55.0, -70.0),
    "beta_anti":   (-139.0, 135.0),
    "beta_par":    (-119.0, 113.0),
    "ppii":        (-75.0, 145.0),
    "ppii_mirror": (75.0, -145.0),
}

# Beta- and gamma-turn templates: (phi,psi) of i+1 and (phi,psi) of i+2.
# Standard turn typology (Venkatachalam 1968; Wilmot & Thornton 1988; Hutchinson
# & Thornton 1994).  These are definitions of named classes, so they are quoted,
# not derived -- but the CA(i)-CA(i+3) gate that admits a turn at all IS derived
# from the ideal chains these templates build (see THRESH["turn_ca13_max_A"]).
TURN_TEMPLATES = {
    "I":    ((-60.0, -30.0), (-90.0, 0.0)),
    "II":   ((-60.0, 120.0), (80.0, 0.0)),
    "I_p":  ((60.0, 30.0), (90.0, 0.0)),
    "II_p": ((60.0, -120.0), (-80.0, 0.0)),
    "VIII": ((-60.0, -30.0), (-120.0, 120.0)),
}
GAMMA_TEMPLATES = {"gamma_classic": (75.0, -64.0), "gamma_inverse": (-79.0, 69.0)}


# --------------------------------------------------------------------------- #
# vector helpers
# --------------------------------------------------------------------------- #
def _u(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def dihedral(p0, p1, p2, p3):
    """Signed dihedral p0-p1-p2-p3 in degrees, IUPAC sign convention."""
    if any(p is None for p in (p0, p1, p2, p3)):
        return None
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    n1 = np.linalg.norm(b1)
    if n1 == 0.0:
        return None
    b1 = b1 / n1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x, y = float(np.dot(v, w)), float(np.dot(np.cross(b1, v), w))
    if x == 0.0 and y == 0.0:
        return None
    return math.degrees(math.atan2(y, x))


def angle_deg(a, b, c):
    """Angle a-b-c in degrees."""
    if any(p is None for p in (a, b, c)):
        return None
    u, v = _u(a - b), _u(c - b)
    return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(u, v))))))


def dist(a, b):
    if a is None or b is None:
        return None
    return float(np.linalg.norm(a - b))


def dev(a, b):
    """Smallest absolute difference between two angles, degrees."""
    d = (a - b) % 360.0
    return min(d, 360.0 - d)


# --------------------------------------------------------------------------- #
# 1. IDEAL CHAIN BUILDER  (used only by --selftest / refs(); never on real data)
# --------------------------------------------------------------------------- #
def _place(a, b, c, r, theta, chi):
    """Natural-extension reference frame: position D with |CD|=r, angle BCD=theta,
    dihedral ABCD=chi.  theta and chi in degrees."""
    th, ch = math.radians(theta), math.radians(chi)
    bc = _u(c - b)
    n = _u(np.cross(b - a, bc))
    m = np.array([bc, np.cross(n, bc), n])
    d = np.array([-r * math.cos(th), r * math.sin(th) * math.cos(ch),
                  r * math.sin(th) * math.sin(ch)])
    return c + d @ m


def build_ideal(n_res, phi, psi, omega=180.0):
    """Build an ideal poly-peptide backbone.  Returns a list of dicts N/CA/C/O.

    phi/psi/omega may be scalars or per-residue sequences.  Ideal bond lengths
    and angles come from IDEAL, which was measured on this reference set.
    """
    def seq(v):
        return [float(v)] * n_res if np.isscalar(v) else [float(x) for x in v]
    phi, psi, omega = seq(phi), seq(psi), seq(omega)
    res = []
    # seed the first three atoms in a plane; the seed only sets the global frame
    N = np.array([0.0, 0.0, 0.0])
    CA = np.array([IDEAL["N_CA"], 0.0, 0.0])
    th = math.radians(180.0 - IDEAL["N_CA_C"])
    C = CA + IDEAL["CA_C"] * np.array([math.cos(th), math.sin(th), 0.0])
    res.append({"N": N, "CA": CA, "C": C})
    for i in range(1, n_res):
        p = res[-1]
        Nn = _place(p["N"], p["CA"], p["C"], IDEAL["C_N"], IDEAL["CA_C_N"], psi[i - 1])
        CAn = _place(p["CA"], p["C"], Nn, IDEAL["N_CA"], IDEAL["C_N_CA"], omega[i - 1])
        Cn = _place(p["C"], Nn, CAn, IDEAL["CA_C"], IDEAL["N_CA_C"], phi[i])
        res.append({"N": Nn, "CA": CAn, "C": Cn})
    # carbonyl O: anti to N(i+1) across the C, i.e. dihedral(N,CA,C,O) = psi + 180
    for i, r in enumerate(res):
        r["O"] = _place(r["N"], r["CA"], r["C"], IDEAL["C_O"], IDEAL["CA_C_O"],
                        psi[i] + 180.0)
    return res


# --------------------------------------------------------------------------- #
# 2. LOCAL HELIX AXIS
# --------------------------------------------------------------------------- #
def helix_params(cas):
    """Local helix parameters from >= 4 consecutive CA positions.

    Method, in full:
      * The bisector at an interior CA, b_i = (P[i-1]-P[i]) + (P[i+1]-P[i]), is
        perpendicular to the axis of a regular helix and points at it.  Stacking
        every interior bisector into a matrix, the axis direction is the null
        direction of that matrix -- the last right singular vector.  For exactly
        two bisectors this reduces to their cross product (Kahn 1989).
      * Sign is fixed so the axis points from the first CA to the last.
      * Project the CAs onto the plane perpendicular to the axis and fit a circle
        algebraically (Kasa).  Radius = fitted radius; rise per residue = mean
        step of the axial coordinate; twist per residue = mean signed step of the
        projected polar angle, positive for a right-handed helix.

    Returns a dict, or {"ok": False, "reason": ...} when the window is not a
    regular helix and a number would be a lie.
    """
    cas = np.asarray(cas, float)
    n = len(cas)
    if n < 4:
        return {"ok": False, "reason": "fewer_than_4_CA"}
    bis = np.array([(cas[i - 1] - cas[i]) + (cas[i + 1] - cas[i]) for i in range(1, n - 1)])
    nb = np.linalg.norm(bis, axis=1)
    if (nb < 1e-6).any():
        return {"ok": False, "reason": "degenerate_bisector_collinear_CA"}
    bis = bis / nb[:, None]
    _, s, vt = np.linalg.svd(bis)
    axis = _u(vt[-1])
    if np.dot(axis, cas[-1] - cas[0]) < 0:
        axis = -axis
    # in-plane frame
    e1 = _u(np.cross(axis, bis[0]))
    if np.linalg.norm(np.cross(axis, bis[0])) < 1e-6:
        return {"ok": False, "reason": "bisector_parallel_to_axis"}
    e2 = np.cross(axis, e1)
    origin = cas.mean(axis=0)
    rel = cas - origin
    z = rel @ axis
    xy = np.stack([rel @ e1, rel @ e2], axis=1)
    # Kasa circle fit: x^2+y^2 = a*x + b*y + c
    A = np.column_stack([xy[:, 0], xy[:, 1], np.ones(n)])
    b = (xy ** 2).sum(axis=1)
    try:
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    except np.linalg.LinAlgError:
        return {"ok": False, "reason": "circle_fit_failed"}
    cx, cy = sol[0] / 2.0, sol[1] / 2.0
    rr = sol[2] + cx * cx + cy * cy
    if rr <= 0:
        return {"ok": False, "reason": "circle_fit_negative_radius"}
    radius = math.sqrt(rr)
    d = np.sqrt(((xy - np.array([cx, cy])) ** 2).sum(axis=1))
    rms = float(np.sqrt(((d - radius) ** 2).mean()))
    ang = np.degrees(np.arctan2(xy[:, 1] - cy, xy[:, 0] - cx))
    steps = np.diff(ang)
    steps = (steps + 180.0) % 360.0 - 180.0
    twist = float(np.mean(steps))
    rise = float(np.mean(np.diff(z)))
    out = {"ok": True, "rise_per_residue_A": round(rise, 3),
           "radius_A": round(radius, 3), "twist_per_residue_deg": round(twist, 2),
           "residues_per_turn": (round(360.0 / twist, 3) if abs(twist) > 1e-6 else None),
           "circle_rms_A": round(rms, 3),
           "handedness": "right" if twist > 0 else "left",
           "axis": [round(float(x), 4) for x in axis],
           "bisector_planarity": round(float(s[1] / s[0]) if s[0] > 0 else 0.0, 4)}
    if abs(twist) > THRESH["axis_max_twist_deg"]:
        return {"ok": False, "reason": "twist_above_150deg_extended_not_helical",
                "twist_per_residue_deg": out["twist_per_residue_deg"],
                "rise_per_residue_A": out["rise_per_residue_A"]}
    if rms > THRESH["axis_max_circle_rms_A"]:
        return {"ok": False, "reason": "circle_fit_residual_above_tolerance",
                "circle_rms_A": out["circle_rms_A"]}
    return out


# --------------------------------------------------------------------------- #
# 3. DERIVED REFERENCE TABLE
# --------------------------------------------------------------------------- #
_REFS = None


def refs():
    """Derive every geometric reference value from the ideal builds.

    Nothing here is a quoted number: each entry is measured on a 12-residue ideal
    chain built from REF_TORSIONS with the reference-measured bond geometry.  The
    only free choices are the acceptance half-widths, which are stated in the
    'halfwidth' fields.
    """
    global _REFS
    if _REFS is not None:
        return _REFS
    out = {}
    for name, (phi, psi) in REF_TORSIONS.items():
        res = build_ideal(12, phi, psi)
        ca = np.array([r["CA"] for r in res])
        rec = {"phi": phi, "psi": psi,
               "ca_i_i1_A": round(float(np.mean([dist(ca[i], ca[i + 1])
                                                 for i in range(11)])), 3)}
        for k in (2, 3, 4, 5):
            vals = [dist(ca[i], ca[i + k]) for i in range(12 - k)]
            rec[f"ca_i_i{k}_A"] = round(float(np.mean(vals)), 3)
        hp = helix_params(ca[3:3 + THRESH["axis_window"]])
        rec["axis"] = hp
        # which n-turn H-bonds does this conformation close, by our own criterion?
        turns = {}
        for n in (3, 4, 5):
            hits = 0
            for i in range(12 - n):
                ok, _d, _a, _e = _hb_geom(res[i + n], res[i + n - 1], res[i])
                hits += 1 if ok else 0
            turns[f"n{n}_turns"] = hits
        rec["hbonds"] = turns
        out[name] = rec
    # turn templates
    for name, ((p1, s1), (p2, s2)) in TURN_TEMPLATES.items():
        # a 4-residue turn: phi/psi of i and i+3 are irrelevant to the CA span,
        # so build with the template pair in the middle and extended flanks
        res = build_ideal(4, [-120.0, p1, p2, -120.0], [120.0, s1, s2, 120.0])
        ca = [r["CA"] for r in res]
        out["turn_" + name] = {"ca_i_i3_A": round(dist(ca[0], ca[3]), 3),
                               "i1": (p1, s1), "i2": (p2, s2)}
    out["_halfwidths"] = {
        # Ramachandran acceptance boxes.  +-40 deg on phi and psi around the ideal
        # point.  Justification: the observed spread of phi/psi inside long
        # reference-set helices is sd 7-10 deg,
        # so +-40 is ~4 sd and admits distorted but genuine helical residues; it
        # is also narrow enough that the alpha_R box does not touch the beta or
        # PPII boxes.  Overlaps that DO remain are listed in out["_overlaps"].
        "rama_halfwidth_deg": 40.0,
        # helix-class separation on rise per residue.  The derived ideal rises are
        # printed in this table; the class boundary is placed at the midpoint of
        # each adjacent pair, so no rise threshold is chosen by hand.
        "rise_boundaries": "midpoints of the derived ideal rises",
    }
    _REFS = out
    return out


def rama_class(phi, psi, halfwidth=None):
    """Which reference conformations is (phi,psi) within halfwidth of?"""
    if phi is None or psi is None:
        return []
    hw = halfwidth if halfwidth is not None else refs()["_halfwidths"]["rama_halfwidth_deg"]
    hit = []
    for name, (rp, rs) in REF_TORSIONS.items():
        if dev(phi, rp) <= hw and dev(psi, rs) <= hw:
            hit.append(name)
    return hit


# --------------------------------------------------------------------------- #
# 4. HYDROGEN BOND GEOMETRY
# --------------------------------------------------------------------------- #
def infer_amide_H(n_pos, ca_pos, prev_c_pos):
    """Position of the backbone amide hydrogen, inferred from heavy atoms.

    A secondary amide N is sp2 with three substituents: CA, the preceding C, and
    H.  The H therefore lies in their plane, opposite the sum of the two
    heavy-atom bond directions.  This uses only the two heavy neighbours that are
    actually present, so it is correct whatever the residue is.

    Returns None when there is no preceding C -- a free N-terminal amine is
    rotatable and its H positions are not determined by the heavy atoms.
    """
    if n_pos is None or ca_pos is None or prev_c_pos is None:
        return None
    d = -(_u(prev_c_pos - n_pos) + _u(ca_pos - n_pos))
    if np.linalg.norm(d) < 1e-6:
        return None
    return n_pos + THRESH["nh_length_A"] * _u(d)


def _hb_geom(don, don_prev, acc):
    """H-bond test between backbone dicts.  don/don_prev/acc carry N/CA/C/O.

    Returns (ok, N...O distance, N-H...O angle, Kabsch-Sander energy).
    """
    N, CA = don.get("N"), don.get("CA")
    O, C = acc.get("O"), acc.get("C")
    if N is None or CA is None or O is None:
        return False, None, None, None
    d = dist(N, O)
    if d is None or d > THRESH["hb_NO_max_A"]:
        return False, d, None, None
    H = infer_amide_H(N, CA, (don_prev or {}).get("C"))
    if H is None:
        return False, d, None, None
    a = angle_deg(N, H, O)
    ok = a is not None and a >= THRESH["hb_NHO_min_deg"]
    e = None
    if C is not None:
        rON, rCH, rOH, rCN = dist(O, N), dist(C, H), dist(O, H), dist(C, N)
        if min(rON, rCH, rOH, rCN) > 0.1:
            e = 0.084 * 332.0 * (1.0 / rON + 1.0 / rCH - 1.0 / rOH - 1.0 / rCN)
    return ok, d, a, e


# --------------------------------------------------------------------------- #
# 5. STRUCTURE LOADING
# --------------------------------------------------------------------------- #
class Res:
    """One residue.  `name` is carried for reporting and is never read by the
    assignment logic."""
    __slots__ = ("name", "seqid", "icode", "auth_chain", "subchain", "het",
                 "atoms", "bb", "is_peptide", "idx", "path_bonds",
                 "n_neighbours", "donor", "chir", "phi", "psi", "omega",
                 "refuse", "flags", "alpha", "bt", "rama_used", "rama_mirrored")

    def __init__(self, res, chain_name, is_peptide):
        self.name = res.name
        self.seqid = res.seqid.num
        self.icode = res.seqid.icode.strip()
        self.auth_chain = chain_name
        self.subchain = res.subchain
        self.het = res.het_flag
        self.atoms = {}
        for at in res:
            if at.element.is_hydrogen:
                continue
            if at.name in self.atoms:      # keep first altloc only; the reference set is
                continue                   # already reduced to one per residue
            self.atoms[at.name] = (np.array([at.pos.x, at.pos.y, at.pos.z]),
                                   at.element.name)
        self.bb = {}
        for k in BB:
            if k in self.atoms:
                self.bb[k] = self.atoms[k][0]
        if "O" not in self.bb and "OXT" in self.atoms:
            self.bb["O"] = self.atoms["OXT"][0]     # C-terminal carboxylate
        self.is_peptide = is_peptide
        self.idx = None
        self.path_bonds = None
        self.n_neighbours = []
        self.donor = None
        self.chir = None
        self.phi = self.psi = self.omega = None
        self.refuse = None
        self.flags = []

    @property
    def label(self):
        """Human-readable identity, FOR REPORTING ONLY.  It contains the residue
        name, so it must never be used as a graph key -- use `key`."""
        return f"{self.name}{self.seqid}{self.icode}:{self.auth_chain}"

    @property
    def key(self):
        """Name-free identity, used for every graph key in this module (H-bond
        sets, bridge partners, neighbour lookups).  Chain plus number plus
        insertion code is unique within an entry and contains no residue name, so
        the audit for residue-name reads stays clean."""
        return f"{self.auth_chain}/{self.seqid}{self.icode}"

    @property
    def complete(self):
        return all(k in self.bb for k in BB)


def _read(path):
    st = gemmi.read_structure(path)
    st.setup_entities()
    return st


def load_entry(entry_dir):
    """Return (peptide residues, receptor/other residues, all heavy atoms).

    The peptide is taken from peptide.cif and the rest from complex.cif minus the
    peptide's atoms, so the partition is exactly the reference set partition and no
    heuristic decides what a peptide is.
    """
    ppath = os.path.join(entry_dir, "peptide.cif")
    cpath = os.path.join(entry_dir, "complex.cif")
    if not os.path.exists(ppath):
        raise FileNotFoundError(ppath)
    pst = _read(ppath)
    pep_keys, pep = set(), []
    for ch in pst[0]:
        for r in ch:
            pep.append(Res(r, ch.name, True))
            pep_keys.add((ch.name, r.seqid.num, r.seqid.icode.strip(), r.name))
    other = []
    if os.path.exists(cpath):
        cst = _read(cpath)
        for ch in cst[0]:
            for r in ch:
                if (ch.name, r.seqid.num, r.seqid.icode.strip(), r.name) in pep_keys:
                    continue
                other.append(Res(r, ch.name, False))
    return pep, other


# --------------------------------------------------------------------------- #
# 6. BACKBONE GRAPH: order, ring closure, chain breaks
# --------------------------------------------------------------------------- #
def amide_links(residues):
    """(i, j, distance) for every C(i)...N(j) below the amide link cutoff."""
    cs = [(i, r.bb["C"]) for i, r in enumerate(residues) if "C" in r.bb]
    ns = [(j, r.bb["N"]) for j, r in enumerate(residues) if "N" in r.bb]
    if not cs or not ns:
        return []
    tc = cKDTree(np.array([p for _, p in cs]))
    tn = cKDTree(np.array([p for _, p in ns]))
    out = []
    for a, b in tc.sparse_distance_matrix(tn, THRESH["link_CN_max_A"]).keys():
        i, j = cs[a][0], ns[b][0]
        if i == j:
            continue
        d = float(np.linalg.norm(cs[a][1] - ns[b][1]))
        if d < THRESH["link_CN_min_A"]:
            continue           # below this is not a bond, it is a modelling error
        out.append((i, j, round(d, 3)))
    return sorted(out)


def unusual_links(pep, links):
    """Amide links whose length is outside the normal range, reported so that a
    stretched closure or an impossible contact is visible in the record."""
    lo, hi = THRESH["link_CN_normal_A"]
    return [{"from_C_of": pep[i].label, "to_N_of": pep[j].label, "c_n_A": d,
             "kind": "short" if d < lo else "stretched"}
            for i, j, d in links if not (lo <= d <= hi)]


def order_backbone(pep):
    """Order the peptide residues by covalent connectivity and describe the
    topology.  Returns (ordered list of indices into pep, info dict).

    No residue numbering is trusted: the order comes from the C->N links alone.
    Numbering is only used, afterwards, to say how many residues are missing at a
    break.
    """
    n = len(pep)
    links = amide_links(pep)
    succ, pred = {}, {}
    extra = []
    for i, j, d in links:
        if i in succ:
            extra.append(("branched_C", pep[i].label, pep[j].label, d))
            continue
        if j in pred:
            extra.append(("branched_N", pep[i].label, pep[j].label, d))
            continue
        succ[i], pred[j] = (j, d), (i, d)
    starts = [i for i in range(n) if i not in pred]
    segments, seen = [], set()
    for s in sorted(starts):
        if s in seen:
            continue
        seg, cur = [], s
        while cur is not None and cur not in seen:
            seg.append(cur)
            seen.add(cur)
            cur = succ[cur][0] if cur in succ else None
        segments.append(seg)
    # anything left is on a cycle
    cyclic_segments = []
    for i in range(n):
        if i in seen:
            continue
        seg, cur = [], i
        while cur not in seen:
            seg.append(cur)
            seen.add(cur)
            cur = succ[cur][0]
        cyclic_segments.append(seg)
    is_cycle = bool(cyclic_segments) and not segments
    if cyclic_segments and segments:
        # a macrocycle with a tail; keep the cycle first, tails after
        pass
    order = []
    for seg in cyclic_segments + segments:
        order.extend(seg)
    for i in range(n):
        if i not in order:
            order.append(i)
    seg_of = {}
    for k, seg in enumerate(cyclic_segments + segments):
        for i in seg:
            seg_of[i] = k
    for i in range(n):
        seg_of.setdefault(i, -1)
    # describe the breaks between consecutive segments in the emitted order
    breaks, moieties = [], []
    allsegs = cyclic_segments + segments
    for k in range(len(allsegs) - 1):
        a, b = pep[allsegs[k][-1]], pep[allsegs[k + 1][0]]
        # A segment whose residues carry no CA at all is not a chain segment: it
        # is a covalent moiety (a farnesyl, a lipid, a non-amino-acid cap) that
        # the reference set keeps inside peptide.cif because it is bonded to the peptide.
        # Calling the join to it a "chain break" would invent a gap that does not
        # exist, so it is reported separately.
        if "CA" not in a.bb or "CA" not in b.bb:
            moieties.append({"between": a.label, "and": b.label,
                             "reason": "one side has no CA; not a backbone segment"})
            continue
        gap = None
        if a.auth_chain == b.auth_chain:
            gap = b.seqid - a.seqid - 1
        breaks.append({
            "after": a.label, "before": b.label,
            "n_unmodelled_residues_by_numbering": gap,
            "ca_ca_A": round(dist(a.bb["CA"], b.bb["CA"]), 3),
            "c_n_A": (round(dist(a.bb.get("C"), b.bb.get("N")), 3)
                      if "C" in a.bb and "N" in b.bb else None),
            "kind": ("unmodelled_gap" if gap and gap > 0 else
                     "sequential_numbering_but_not_bonded" if gap == 0 else
                     "separate_chain_or_unknown"),
        })
    info = {
        "n_residues": n,
        "n_amide_links": len(links),
        "is_macrocycle_by_backbone_amide": is_cycle,
        "neighbour_indices_wrap": is_cycle,
        "n_segments": len(allsegs),
        "segment_lengths": [len(s) for s in allsegs],
        "breaks": breaks,
        "unusual_amide_links": unusual_links(pep, links),
        # every break stops the neighbour walk, so all of them are counted here;
        # the second field is the subset that is an actual unmodelled gap, which
        # is what a consumer looking for missing density wants.
        "n_internal_chain_breaks": len(breaks),
        "n_internal_breaks_that_are_unmodelled_gaps":
            sum(1 for b in breaks if (b["n_unmodelled_residues_by_numbering"] or 0) > 0),
        "n_unmodelled_residues_at_internal_breaks":
            sum(b["n_unmodelled_residues_by_numbering"] or 0 for b in breaks),
        "non_backbone_moiety_joins": moieties,
        "branched_links": extra,
        # A branched amide link means one backbone C reaches two nitrogens (or one
        # N two carbons) within bonding distance.  The backbone order is then
        # genuinely ambiguous and the walk below picks the first link found.  The
        # branches are reported so the apply pass can exclude such an entry rather
        # than trust an arbitrary order.  7TO8 is the reference set example: six branches.
        "backbone_order_ambiguous": bool(extra),
        "closure_link": None,
        "note_terminal_truncation": (
            "residues unmodelled at the TERMINI leave no trace in the "
            "coordinates and are invisible here; take them from the manifest "
            "columns n_res_unmodelled / peptide_unmodelled_fraction"),
    }
    if is_cycle and cyclic_segments:
        seg = cyclic_segments[0]
        a, b = pep[seg[-1]], pep[seg[0]]
        info["closure_link"] = {
            "from_C_of": a.label, "to_N_of": b.label,
            "c_n_A": round(dist(a.bb["C"], b.bb["N"]), 3),
            "span_residues": len(seg),
            "kind": "head_to_tail_backbone_amide",
        }
    return order, info


def other_covalent_rings(pep, order, info):
    """Intra-peptide covalent links that are NOT backbone amides: disulfide,
    thioether, lactam, staple, ester.  These make a ring too, and a ring changes
    which residues are spatial neighbours even when the backbone stays linear.

    Detected from heavy-atom distances with element-aware covalent cutoffs.  No
    atom-name or residue-name list is used: any side-chain-to-side-chain or
    side-chain-to-backbone bond between residues that are not backbone
    neighbours is reported.
    """
    pos, key = [], []
    for i, r in enumerate(pep):
        for an, (p, el) in r.atoms.items():
            pos.append(p)
            key.append((i, an, el))
    if not pos:
        return []
    tree = cKDTree(np.array(pos))
    nbr = {}
    for i, j, d in info.get("_links", []):
        nbr[(i, j)] = True
    seqnbr = set()
    for a, b in zip(order, order[1:]):
        seqnbr.add((a, b))
        seqnbr.add((b, a))
    if info["is_macrocycle_by_backbone_amide"] and len(order) > 1:
        seqnbr.add((order[-1], order[0]))
        seqnbr.add((order[0], order[-1]))
    rad = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07, "SE": 1.20,
           "F": 0.57, "CL": 1.02, "BR": 1.20, "I": 1.39}
    out = []
    for a, b in tree.sparse_distance_matrix(tree, 2.6).keys():
        if a >= b:
            continue
        i, an, ei = key[a]
        j, bn, ej = key[b]
        if i == j:
            continue
        d = float(np.linalg.norm(pos[a] - pos[b]))
        cut = (rad.get(ei.upper(), 0.8) + rad.get(ej.upper(), 0.8)
               + THRESH["bond_slack_A"])
        if d > cut:
            continue
        if (i, j) in seqnbr and {an, bn} == {"C", "N"}:
            continue                      # the backbone amide itself
        out.append({"res_a": pep[i].label, "atom_a": an,
                    "res_b": pep[j].label, "atom_b": bn,
                    "dist_A": round(d, 3),
                    "backbone_neighbours": (i, j) in seqnbr,
                    "elements": f"{ei}-{ej}"})
    return out


# --------------------------------------------------------------------------- #
# 7. DONOR AVAILABILITY  (the fix for the DSSP defect)
# --------------------------------------------------------------------------- #
def donor_state(res, prev_res, env_tree, env_pos, env_key, env_res):
    """Does this residue's backbone N carry a hydrogen?

    Decided from heavy-atom connectivity only, because the reference set has no
    hydrogens.  Every heavy atom in the whole complex within the N-bond cutoff is
    counted, so an N-methyl, a proline-type ring closure, an N-terminal cap in a
    different residue and a side-chain-to-backbone-N crosslink are all caught by
    the same test.

    A neutral amide N is trivalent.  Bonded heavy neighbours:
        1  -> free terminal amine, NH2 / NH3+  : donor available, H NOT positioned
        2  -> secondary amide, N-H             : donor available, H positioned
        >=3 -> tertiary N, no hydrogen at all  : NO DONOR.  Refuse the H-bond.
    """
    if "N" not in res.bb:
        return {"available": None, "reason": "no_backbone_N_modelled",
                "n_heavy_neighbours": None, "neighbours": []}
    N = res.bb["N"]
    nb = []
    for k in env_tree.query_ball_point(N, THRESH["n_bond_max_A"]):
        ri, an = env_key[k]
        if env_res[ri] is res and an == "N":
            continue
        d = float(np.linalg.norm(N - env_pos[k]))
        if d < 0.4:
            continue
        nb.append({"res": env_res[ri].label, "atom": an, "dist_A": round(d, 3),
                   "same_residue": env_res[ri] is res,
                   "is_preceding_backbone_C": (prev_res is not None
                                               and env_res[ri] is prev_res
                                               and an == "C"),
                   "borderline": (THRESH["n_bond_borderline_A"][0] <= d
                                  <= THRESH["n_bond_borderline_A"][1])})
    nb.sort(key=lambda x: x["dist_A"])
    k = len(nb)
    if k >= 3:
        extra = [x for x in nb if not (x["same_residue"] and x["atom"] == "CA")
                 and not x["is_preceding_backbone_C"]]
        who = ", ".join(f"{x['atom']} of {x['res']} at {x['dist_A']} A" for x in extra) \
            or "unclear which neighbour is the third"
        return {"available": False,
                "reason": f"tertiary_N_{k}_heavy_substituents_no_amide_H ({who})",
                "n_heavy_neighbours": k, "neighbours": nb}
    if k == 2:
        return {"available": True, "reason": "secondary_amide_NH",
                "n_heavy_neighbours": 2, "neighbours": nb, "h_positioned": True}
    if k == 1:
        return {"available": True,
                "reason": "free_terminal_amine_NH2_or_NH3_rotatable_H_not_positioned",
                "n_heavy_neighbours": 1, "neighbours": nb, "h_positioned": False}
    return {"available": None, "reason": "isolated_N_no_heavy_neighbour_found",
            "n_heavy_neighbours": 0, "neighbours": nb}


# --------------------------------------------------------------------------- #
# 8. BACKBONE TYPE AND CHIRALITY
# --------------------------------------------------------------------------- #
def intra_bond_graph(res):
    """Heavy-atom bond graph inside one residue, from element-aware distances."""
    rad = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07, "SE": 1.20,
           "F": 0.57, "CL": 1.02, "BR": 1.20, "I": 1.39}
    names = list(res.atoms)
    g = {a: set() for a in names}
    for a in range(len(names)):
        for b in range(a + 1, len(names)):
            pa, ea = res.atoms[names[a]]
            pb, eb = res.atoms[names[b]]
            cut = (rad.get(ea.upper(), 0.8) + rad.get(eb.upper(), 0.8)
                   + THRESH["bond_slack_A"])
            if np.linalg.norm(pa - pb) <= cut:
                g[names[a]].add(names[b])
                g[names[b]].add(names[a])
    return g


def backbone_type(res):
    """Is this an alpha-amino-acid backbone?  Returns (ok, detail dict).

    PRIMARY test: the shortest covalent path from the backbone N to the backbone
    C, inside the residue, must be exactly two bonds (N-CA-C).  A beta3-amino
    acid gives three; a gamma residue whose CCD names the alpha carbon two bonds
    from N gives three; a non-amino-acid linker gives more or nothing.
    CONFIRMING test: N-CA and CA-C lengths inside the measured windows.
    """
    d = {"n_ca_A": None, "ca_c_A": None, "n_to_c_path_bonds": None,
         "path": None}
    if not all(k in res.bb for k in ("N", "CA", "C")):
        return False, dict(d, reason="incomplete_backbone_N_CA_C")
    d["n_ca_A"] = round(dist(res.bb["N"], res.bb["CA"]), 3)
    d["ca_c_A"] = round(dist(res.bb["CA"], res.bb["C"]), 3)
    g = intra_bond_graph(res)
    # BFS N -> C
    from collections import deque
    q, prev = deque([("N", 0)]), {"N": None}
    hops = None
    while q:
        a, h = q.popleft()
        if a == "C":
            hops = h
            break
        for b in sorted(g.get(a, ())):
            if b not in prev:
                prev[b] = a
                q.append((b, h + 1))
    d["n_to_c_path_bonds"] = hops
    if hops is not None:
        path, cur = [], "C"
        while cur is not None:
            path.append(cur)
            cur = prev[cur]
        d["path"] = "-".join(reversed(path))
    lo, hi = THRESH["n_ca_range_A"]
    lo2, hi2 = THRESH["ca_c_range_A"]
    ok_nca = lo <= d["n_ca_A"] <= hi
    ok_cac = lo2 <= d["ca_c_A"] <= hi2
    ok_path = hops == THRESH["backbone_path_bonds"]
    if ok_path and ok_nca and ok_cac:
        return True, dict(d, reason="alpha_amino_acid_backbone")
    bad = []
    if not ok_path:
        bad.append(f"N->C covalent path is {hops} bonds ({d['path']}), "
                   f"an alpha backbone has {THRESH['backbone_path_bonds']}")
    if not ok_nca:
        bad.append(f"N-CA {d['n_ca_A']} A outside {THRESH['n_ca_range_A']}")
    if not ok_cac:
        bad.append(f"CA-C {d['ca_c_A']} A outside {THRESH['ca_c_range_A']}")
    return False, dict(d, reason="not_an_alpha_backbone: " + "; ".join(bad))


def chirality(res, alpha_ok, bt):
    """D / L at the alpha carbon, from the improper dihedral N-C-CB-CA.

    Refuses -- returns a null value with a reason -- when the backbone is not
    alpha (this is the 4BPJ HR7 case: N sits 2.44 A from CA, so the improper is
    not a torsion about real bonds and the confident D call the build pass made
    is meaningless), when CA carries no substituent (glycine-like), when it
    carries two (Aib and every other quaternary centre: achiral at CA whatever
    the component dictionary says), or when the centre is too planar to call.
    """
    out = {"label": None, "improper_deg": None, "reason": None,
           "n_ca_substituents": None}
    if not alpha_ok:
        out["reason"] = "refused_non_alpha_backbone: " + bt["reason"]
        return out
    g = intra_bond_graph(res)
    subs = sorted(a for a in g.get("CA", ()) if a not in ("N", "C"))
    out["n_ca_substituents"] = len(subs)
    out["ca_substituent_atoms"] = subs
    if len(subs) == 0:
        out["reason"] = "achiral_no_substituent_at_CA_glycine_like"
        return out
    if len(subs) > 1:
        out["reason"] = (f"achiral_quaternary_CA_{len(subs)}_substituents "
                         f"({','.join(subs)})")
        return out
    imp = dihedral(res.bb["N"], res.bb["C"], res.atoms[subs[0]][0], res.bb["CA"])
    if imp is None:
        out["reason"] = "degenerate_improper"
        return out
    out["improper_deg"] = round(imp, 2)
    if abs(imp) < THRESH["chiral_min_abs_deg"]:
        out["reason"] = f"refused_near_planar_CA_|improper|={abs(imp):.2f}_deg"
        return out
    out["label"] = "L" if imp * THRESH["chiral_L_sign"] > 0 else "D"
    out["reason"] = "improper_dihedral_N-C-{}-CA".format(subs[0])
    return out


# --------------------------------------------------------------------------- #
# 9. THE ASSIGNER
# --------------------------------------------------------------------------- #
class Assignment:
    def __init__(self, entry_dir, pdb_id=None):
        self.entry_dir = entry_dir
        self.pdb_id = pdb_id or os.path.basename(os.path.normpath(entry_dir))
        self.pep, self.other = load_entry(entry_dir)
        self.notes = []
        self._prepare()

    # ---------------------------------------------------------------- setup
    def _prepare(self):
        pep, other = self.pep, self.other
        order, info = order_backbone(pep)
        self.order = order
        self.topo = info
        self.chain = [pep[i] for i in order]
        for k, r in enumerate(self.chain):
            r.idx = k
        self.wrap = info["neighbour_indices_wrap"]
        # segment id per position, so neighbour walks never cross a break
        self.segid = []
        s = 0
        for k in range(len(self.chain)):
            if k > 0 and self._broken_between(k - 1, k):
                s += 1
            self.segid.append(s)
        self.topo["crosslinks"] = other_covalent_rings(pep, order, info)
        # A break that a recorded NON-amide covalent bond spans is not an
        # unmodelled gap: it is a different chemistry joining the two residues --
        # an isoaspartyl CG-N link, a macrolactone ester, a thioester.  The
        # neighbour walk still stops there, because an i->i+4 hydrogen-bond
        # geometry across a backbone that is one atom longer is not comparable to
        # the alpha reference.  But the record must say which it is.
        labelpairs = {(c["res_a"], c["res_b"]) for c in self.topo["crosslinks"]}
        labelpairs |= {(b, a) for a, b in labelpairs}
        for b in self.topo["breaks"]:
            spans = [c for c in self.topo["crosslinks"]
                     if {c["res_a"], c["res_b"]} == {b["after"], b["before"]}]
            b["spanned_by_non_amide_covalent_link"] = spans or None
            if spans:
                b["kind"] = "non_amide_covalent_link_not_a_gap"
        self.topo["is_macrocycle_any_ring"] = bool(
            info["is_macrocycle_by_backbone_amide"]
            or any(not c["backbone_neighbours"] for c in self.topo["crosslinks"]))

        # environment tree over EVERY heavy atom of the complex, for N-substituent
        # detection.  Includes hetero and receptor, so a cap on another entity is
        # found.
        self._env_res, self._env_key, epos = [], [], []
        for r in pep + other:
            ri = len(self._env_res)
            self._env_res.append(r)
            for an, (p, _el) in r.atoms.items():
                self._env_key.append((ri, an))
                epos.append(p)
        self._env_pos = np.array(epos) if epos else np.zeros((0, 3))
        self._env_tree = cKDTree(self._env_pos) if len(epos) else None

        # receptor residues in the shell, ordered into chains by the same
        # connectivity rule
        self.rec = self._receptor_shell()

        # per-residue basics
        for r in self.chain:
            ok, bt = backbone_type(r)
            r.flags = []
            r.alpha = ok
            r.bt = bt
            r.refuse = None if ok else bt["reason"]
            r.chir = chirality(r, ok, bt)
        for k, r in enumerate(self.chain):
            prev = self._nbr(k, -1)
            r.donor = donor_state(r, self.chain[prev] if prev is not None else None,
                                  self._env_tree, self._env_pos, self._env_key,
                                  self._env_res)
        for r in self.rec:
            r.alpha, r.bt = backbone_type(r)
        self._rec_prev = self._receptor_prev_map()
        for r in self.rec:
            r.donor = donor_state(r, self._rec_prev.get(id(r)), self._env_tree,
                                  self._env_pos, self._env_key, self._env_res)
        self._torsions()
        self._hbonds()

    def _broken_between(self, a, b):
        """True if chain positions a and b are not joined by a backbone amide."""
        ra, rb = self.chain[a], self.chain[b]
        if "C" not in ra.bb or "N" not in rb.bb:
            return True
        d = dist(ra.bb["C"], rb.bb["N"])
        return not (THRESH["link_CN_min_A"] <= d <= THRESH["link_CN_max_A"])

    def _nbr(self, k, step):
        """Chain position k + step, walking the covalent backbone.

        WRAPS around a macrocycle closure -- this is the single most important
        line of the module for the 476 macrocyclic entries, because an
        implementation that stops at the end of the residue list silently loses
        every i->i+3 and i->i+4 relationship that crosses the closure.

        Returns None when a chain break intervenes (an unmodelled gap is a break,
        not a shortcut), and None when |step| >= the ring size, because in a ring
        of n residues position i+n is i itself and i+k for k>=n is not a distinct
        relationship at all.
        """
        n = len(self.chain)
        if n == 0:
            return None
        if self.wrap and abs(step) >= n:
            return None
        cur = k
        d = 1 if step > 0 else -1
        for _ in range(abs(step)):
            nxt = cur + d
            if nxt < 0 or nxt >= n:
                if not self.wrap:
                    return None
                nxt %= n
            # _broken_between(a, b) tests C of a against N of b, which is the
            # correct test for the wrap pair (last, first) as well as for any
            # adjacent pair, so no special case is needed here.
            a, b = (cur, nxt) if d > 0 else (nxt, cur)
            if self._broken_between(a, b):
                return None
            cur = nxt
        return cur

    def _ring_sep(self, i, j):
        """Separation between two peptide positions along the backbone, counted in
        residues and taking the short way round a macrocycle.  Returns a large
        number when the two positions are in different segments, because they are
        then not related along the chain at all."""
        if i == j:
            return 0
        n = len(self.chain)
        for step in range(1, n):
            if self._nbr(i, step) == j or self._nbr(i, -step) == j:
                return step
        return 10 ** 6

    def _receptor_shell(self):
        if self._env_tree is None or not self.pep:
            return []
        pp = np.array([p for r in self.pep for p, _ in r.atoms.values()])
        if len(pp) == 0:
            return []
        tree = cKDTree(pp)
        keep = []
        for r in self.other:
            if not r.atoms:
                continue
            q = np.array([p for p, _ in r.atoms.values()])
            if tree.query(q)[0].min() <= THRESH["receptor_shell_A"]:
                keep.append(r)
        return [r for r in keep if r.complete]

    def _receptor_prev_map(self):
        """Preceding backbone C partner for each shell residue, by connectivity.
        Searched over ALL other residues, not just the shell, so a residue at the
        shell edge still gets its real predecessor."""
        cs = [(r, r.bb["C"]) for r in self.other if "C" in r.bb]
        if not cs:
            return {}
        t = cKDTree(np.array([p for _, p in cs]))
        out = {}
        for r in self.rec:
            if "N" not in r.bb:
                continue
            for k in t.query_ball_point(r.bb["N"], THRESH["link_CN_max_A"]):
                if cs[k][0] is r:
                    continue
                out[id(r)] = cs[k][0]
                break
        return out

    # ---------------------------------------------------------------- torsions
    def _torsions(self):
        for k, r in enumerate(self.chain):
            p, n = self._nbr(k, -1), self._nbr(k, 1)
            rp = self.chain[p] if p is not None else None
            rn = self.chain[n] if n is not None else None
            if not r.alpha:
                r.phi = r.psi = r.omega = None
                r.flags.append("torsions_null_non_alpha_backbone")
                continue
            if rp is not None and "C" in rp.bb and r.alpha:
                r.phi = dihedral(rp.bb["C"], r.bb["N"], r.bb["CA"], r.bb["C"]) \
                    if all(x in r.bb for x in ("N", "CA", "C")) else None
            if rn is not None and "N" in rn.bb:
                r.psi = dihedral(r.bb["N"], r.bb["CA"], r.bb["C"], rn.bb["N"]) \
                    if all(x in r.bb for x in ("N", "CA", "C")) else None
                if "CA" in rn.bb:
                    r.omega = dihedral(r.bb["CA"], r.bb["C"], rn.bb["N"], rn.bb["CA"])
            # torsions of a residue whose neighbour is non-alpha are still real
            # torsions about real bonds, so they are kept; the neighbour's own
            # phi/psi are the ones that are null.
        # mirrored view: the (phi,psi) a D-residue would have if it were L
        for r in self.chain:
            if r.phi is None or r.psi is None:
                r.rama_used = (None, None)
                r.rama_mirrored = False
                continue
            if r.chir.get("label") == "D":
                r.rama_used = (-r.phi, -r.psi)
                r.rama_mirrored = True
            else:
                r.rama_used = (r.phi, r.psi)
                r.rama_mirrored = False

    # ---------------------------------------------------------------- H-bonds
    def _hbonds(self):
        """Every backbone N->O H-bond involving the peptide.

        Three sets are produced from the same geometry:
          hb          the honest set: donor availability enforced
          hb_nodonor  donor availability ignored (what DSSP assumes)
        Each bond records donor, acceptor, N...O, N-H...O, KS energy, whether
        both partners are in the peptide, and whether the donor's H position was
        determined by the heavy atoms.
        """
        # index every candidate donor (has N and CA) and acceptor (has O)
        don = [("P", k, r) for k, r in enumerate(self.chain) if "N" in r.bb and "CA" in r.bb]
        don += [("R", None, r) for r in self.rec if "N" in r.bb and "CA" in r.bb]
        acc = [("P", k, r) for k, r in enumerate(self.chain) if "O" in r.bb]
        acc += [("R", None, r) for r in self.rec if "O" in r.bb]
        if not don or not acc:
            self.hb, self.hb_nodonor, self.refused, self.near = [], [], [], []
            self.phantom_dssp = []
            return
        at = cKDTree(np.array([r.bb["O"] for _s, _k, r in acc]))
        prev_of = {}
        for k, r in enumerate(self.chain):
            p = self._nbr(k, -1)
            prev_of[id(r)] = self.chain[p] if p is not None else None
        for r in self.rec:
            prev_of[id(r)] = self._rec_prev.get(id(r))
        hb, hbn, refused, near = [], [], [], []
        for sd, kd, rd in don:
            for ai in at.query_ball_point(rd.bb["N"], THRESH["hb_NO_max_A"]):
                sa, ka, ra = acc[ai]
                if ra is rd:
                    continue
                if sd == "R" and sa == "R":
                    continue            # receptor-internal: not our business
                prev = prev_of.get(id(rd))
                if prev is not None and ra is prev:
                    # the carbonyl O of the immediately preceding residue is part
                    # of the same amide group as this N.  Its O sits ~2.25 A from
                    # the N by covalent geometry alone and it is never an H-bond
                    # acceptor from it.  Excluded explicitly rather than being
                    # left to the angle test, so that the near-miss diagnostic
                    # below counts real candidates and not one artefact per
                    # residue.
                    continue
                ok, d, a, e = _hb_geom(rd.bb, (prev.bb if prev else None) or {}, ra.bb)
                # the "as-if-DSSP" pass: place H from the same heavy atoms but do
                # not ask whether an H exists
                rec = {
                    "donor": rd.label, "donor_in": sd, "donor_pos": kd,
                    "donor_key": rd.key, "acceptor_key": ra.key,
                    "acceptor": ra.label, "acceptor_in": sa, "acceptor_pos": ka,
                    "N_O_A": round(d, 3) if d is not None else None,
                    "N_H_O_deg": round(a, 1) if a is not None else None,
                    "ks_energy_kcal": round(e, 3) if e is not None else None,
                    "class": ("peptide_internal" if sd == "P" and sa == "P"
                              else "peptide_to_receptor" if sd == "P"
                              else "receptor_to_peptide"),
                    "donor_h_positioned": bool(rd.donor.get("h_positioned")),
                    "donor_available": rd.donor.get("available"),
                    "donor_reason": rd.donor.get("reason"),
                }
                if not ok:
                    # within the distance criterion but failing the angle.  Kept
                    # as a count so that a threshold argument can be audited
                    # without re-running: these are the pairs a
                    # distance-only H-bond criterion would have accepted.
                    if d is not None and d <= THRESH["hb_NO_max_A"]:
                        near.append(dict(rec, failed="angle_or_no_H_position"))
                    continue
                hbn.append(rec)
                if rd.donor.get("available") is True:
                    hb.append(rec)
                else:
                    refused.append(rec)
        self.hb, self.hb_nodonor, self.refused, self.near = hb, hbn, refused, near
        self.phantom_dssp = self._phantom_by_dssp_energy(acc, at, prev_of)

    def _phantom_by_dssp_energy(self, acc, at, prev_of):
        """H-bonds DSSP WOULD score at a nitrogen that carries no hydrogen.

        The count of refused bonds above is measured against THIS module's
        criterion (3.5 A, 120 deg), which is stricter than DSSP's.  DSSP accepts
        any pair whose Kabsch-Sander electrostatic energy is below -0.5 kcal/mol,
        which reaches out to about 5 A at a favourable angle.  To state DSSP's
        exposure fairly rather than flatteringly, this scans donor-less nitrogens
        out to 5.5 A, places the hydrogen DSSP would have assumed, and applies
        DSSP's own energy cutoff.
        """
        out = []
        for k, rd in enumerate(self.chain):
            if rd.donor.get("available") is not False:
                continue
            if "N" not in rd.bb or "CA" not in rd.bb:
                continue
            prev = prev_of.get(id(rd))
            H = infer_amide_H(rd.bb["N"], rd.bb["CA"],
                              (prev.bb if prev else {}).get("C"))
            if H is None:
                continue
            for ai in at.query_ball_point(rd.bb["N"], 5.5):
                sa, ka, ra = acc[ai]
                if ra is rd or (prev is not None and ra is prev):
                    continue
                if "C" not in ra.bb:
                    continue
                N, O, C = rd.bb["N"], ra.bb["O"], ra.bb["C"]
                rON, rCH, rOH, rCN = dist(O, N), dist(C, H), dist(O, H), dist(C, N)
                if min(rON, rCH, rOH, rCN) < 0.1:
                    continue
                e = 0.084 * 332.0 * (1.0 / rON + 1.0 / rCH - 1.0 / rOH - 1.0 / rCN)
                if e <= THRESH["ks_energy_max_kcal"]:
                    out.append({"donor": rd.label, "acceptor": ra.label,
                                "acceptor_in": sa, "N_O_A": round(rON, 3),
                                "ks_energy_kcal": round(e, 3),
                                "donor_reason": rd.donor.get("reason"),
                                "meets_this_modules_criterion":
                                    rON <= THRESH["hb_NO_max_A"]
                                    and (angle_deg(N, H, O) or 0.0)
                                    >= THRESH["hb_NHO_min_deg"]})
        return out

    # ------------------------------------------------------------ SS per mode
    def _hbset(self, mode):
        """(donor,acceptor) position pairs visible in this mode, plus the raw
        list.  Positions are chain positions for the peptide and id() for the
        receptor."""
        src = self.hb_nodonor if mode == "nodonor" else self.hb
        out = []
        for h in src:
            if mode == "intrinsic" and h["class"] != "peptide_internal":
                continue
            out.append(h)
        return out

    def assign(self, mode="bound"):
        """Run the full cascade.  mode in {bound, intrinsic, nodonor}."""
        n = len(self.chain)
        code = ["C"] * n
        detail = [dict() for _ in range(n)]
        hbs = self._hbset(mode)

        # ---- n-turn table: nturn[(pos_i, n)] means CO(i) -> NH(i+n)
        nturn = set()
        pep_hb = set()
        for h in hbs:
            if h["donor_in"] == "P" and h["acceptor_in"] == "P":
                pep_hb.add((h["acceptor_pos"], h["donor_pos"]))
        # n = 2 is the gamma-turn relationship, 3/4/5 are the helical ones.
        for (i, j) in pep_hb:
            for nn in (2, 3, 4, 5):
                if self._nbr(i, nn) == j:
                    nturn.add((i, nn))

        # ---- helices: two consecutive n-turns.
        # A pair of consecutive n-turns at (i, i+1) closes one full helical turn
        # and marks residues i+1 .. i+n.  Nothing here is chosen by hand: n is 3,
        # 4 or 5 and the number of residues marked is n.
        helix_hits = {}
        for (i, nn) in sorted(nturn):
            if nn < 3:
                continue                  # n=2 is a gamma turn, not a helix
            j = self._nbr(i, 1)
            if j is None or (j, nn) not in nturn:
                continue
            span = [self._nbr(i, s) for s in range(1, nn + 1)]
            if any(s is None for s in span):
                continue
            for s in span:
                helix_hits.setdefault(s, set()).add(nn)
        # A residue covered by both a 4-turn and a 3-turn belongs to the 4-turn.
        # After that assignment a class can be left with a run shorter than its
        # own turn -- one residue calling itself a 3-10 helix in the middle of an
        # alpha-helix, which happens at a helix kink.  Such a residue is demoted
        # to its next available turn size, and out of the helix classes entirely
        # if it has none.  This is a consistency rule, not a threshold: a helix
        # of n residues per turn cannot be shorter than n residues.
        helix_n = {}
        pool = {k: set(v) for k, v in helix_hits.items()}
        for _ in range(4):
            helix_n = {}
            for k, ns in pool.items():
                if not ns:
                    continue
                helix_n[k] = 4 if 4 in ns else (3 if 3 in ns else 5)
            demoted = False
            for k, nn in list(helix_n.items()):
                run = [k]
                cur = k
                while True:
                    j = self._nbr(cur, 1)
                    if j is None or j in run or helix_n.get(j) != nn:
                        break
                    run.append(j)
                    cur = j
                cur = k
                while True:
                    j = self._nbr(cur, -1)
                    if j is None or j in run or helix_n.get(j) != nn:
                        break
                    run.append(j)
                    cur = j
                if len(run) < nn:
                    for p in run:
                        pool[p].discard(nn)
                    demoted = True
            if not demoted:
                break
        helix_hits = {k: pool[k] for k in pool if pool[k]}
        # Handedness is a property of the helical SEGMENT, not of one residue, so
        # it is decided once per segment: the local-axis fit over all of the
        # segment's CA positions when that fit succeeds, otherwise the sign of the
        # segment's mean phi (a right-handed alpha-helix has phi < 0 and its exact
        # mirror image has phi > 0, so the sign is not a fitted quantity).
        helix_class, helix_seg = {}, {}
        for k, nn in sorted(helix_n.items()):
            if k in helix_seg:
                continue
            run, cur = [k], k
            while True:
                j = self._nbr(cur, 1)
                if j is None or j in run or helix_n.get(j) != nn:
                    break
                run.append(j)
                cur = j
            cas = [self.chain[p].bb["CA"] for p in run
                   if "CA" in self.chain[p].bb]
            hp = helix_params(cas) if len(cas) >= 4 else {"ok": False,
                                                         "reason": "segment_shorter_than_4"}
            if hp.get("ok"):
                hand = hp["handedness"]
                src = "segment_axis_fit"
            else:
                phis = [self.chain[p].phi for p in run
                        if self.chain[p].phi is not None]
                hand = ("right" if np.mean(phis) < 0 else "left") if phis else None
                src = "sign_of_mean_phi"
            ch = {3: "G", 4: "H", 5: "I"}[nn]
            if hand == "left" and nn in (3, 4):
                ch = ch.lower()
            for p in run:
                helix_class[p] = ch
                helix_seg[p] = {"n_turn": nn, "handedness": hand,
                                "handedness_from": src,
                                "segment_positions": run,
                                "segment_axis": hp}

        # ---- bridges (DSSP's parallel/antiparallel bridge, wrap-aware)
        bridges = self._bridges(hbs)

        # ---- turns
        turns, gammas = self._turns(nturn)

        # ---- local axis per position
        axis = {}
        w = THRESH["axis_window"]
        for k in range(n):
            win = [self._nbr(k, s) for s in range(-(w // 2), w - w // 2)]
            if any(s is None for s in win):
                continue
            if any("CA" not in self.chain[s].bb for s in win):
                continue
            axis[k] = helix_params([self.chain[s].bb["CA"] for s in win])
            axis[k]["window_positions"] = win

        ppii = self._ppii()
        ext = self._extended()

        for k, r in enumerate(self.chain):
            d = detail[k]
            if not r.complete:
                code[k] = "X"
                d["reason"] = "incomplete_backbone_" + ",".join(
                    x for x in BB if x not in r.bb)
                continue
            if not r.alpha:
                code[k] = "?"
                d["reason"] = r.refuse
                continue
            if k in helix_class:
                code[k] = helix_class[k]
                d["helix"] = dict(helix_seg[k], n_turns_available=sorted(helix_hits[k]))
            elif k in bridges:
                partners = bridges[k]
                code[k] = "E" if any(p["partner_in"] == "R" for p in partners) else "e"
                d["bridge_partners"] = partners
            elif k in turns:
                code[k] = "T"
                d["turn"] = turns[k]
            elif k in gammas:
                code[k] = "S"
                d["gamma_turn"] = gammas[k]
            elif k in ppii:
                code[k] = "p" if ppii[k]["mirror"] else "P"
                d["ppii"] = ppii[k]
            elif k in ext:
                code[k] = "x"
                d["extended"] = ext[k]
            if k in axis:
                d["axis"] = axis[k]
        return {"mode": mode, "string": "".join(code), "detail": detail,
                "nturns": sorted(nturn), "bridges": bridges, "turns": turns,
                "gamma_turns": gammas, "ppii": ppii, "n_hbonds_visible": len(hbs)}

    def _bridges(self, hbs):
        """DSSP bridge definition, evaluated over peptide and receptor together.

        TWO CONVENTIONS ARE IN PLAY AND THEY ARE OPPOSITE.  Getting this wrong is
        the defect corrected on 2026-08-03 (see the dated note at the clause).

          * Kabsch & Sander write Hbond(i, j) for the electrostatic interaction
            between the C=O group of residue i and the N-H group of residue j.
            So in K&S notation the FIRST index ACCEPTS and the second DONATES.
          * `HB` in this method is a set of (donor_key, acceptor_key) pairs, built
            straight from _hbonds(): the FIRST element DONATES its N-H, the second
            ACCEPTS with its C=O.
          * Therefore  K&S Hbond(i, j)  is  HB[(j, i)]  here.  The tuple order is
            REVERSED between the two, always, with no exception.

        Kabsch & Sander's bridge definition, in their own notation:
            parallel      if [Hbond(i-1, j) and Hbond(j, i+1)]
                          or [Hbond(j-1, i) and Hbond(i, j+1)]
            antiparallel  if [Hbond(i, j) and Hbond(j, i)]
                          or [Hbond(i-1, j+1) and Hbond(j-1, i+1)]

        Translated into this module's (donor, acceptor) order, with am/ap = a-+1
        and bm/bp = b-+1, which is what the code below must contain:
            parallel      if [(b, am) and (ap, b)]  or  [(a, bm) and (bp, a)]
            antiparallel  if [(a, b) and (b, a)]    or  [(bp, am) and (ap, bm)]

        The narrow antiparallel clause is symmetric in the two tuples, so it reads
        the same under either convention; that is why only the two WIDE clauses
        were wrong.  Donor/acceptor are identified by `key`, not by residue name,
        so peptide and receptor residues can be mixed freely.
        """
        HB = set()
        for h in hbs:
            HB.add((h["donor_key"], h["acceptor_key"]))
        pos_of, label_of = {}, {}
        for k, r in enumerate(self.chain):
            pos_of[r.key] = ("P", k)
            label_of[r.key] = r.label
        for r in self.rec:
            pos_of.setdefault(r.key, ("R", None))
            label_of.setdefault(r.key, r.label)
        # neighbour lookup by label
        def lbl(sd, k, step):
            if sd == "P":
                j = self._nbr(k, step)
                return self.chain[j].key if j is not None else None
            return None
        rec_next, rec_prev = {}, {}
        for r in self.rec:
            p = self._rec_prev.get(id(r))
            if p is not None:
                rec_prev[r.key] = p.key
                rec_next[p.key] = r.key

        def nb(label, step):
            s = pos_of.get(label)
            if s is None:
                return None
            if s[0] == "P":
                return lbl("P", s[1], step)
            return rec_next.get(label) if step == 1 else rec_prev.get(label)

        # Candidate pairs.  It is NOT enough to test only pairs that are
        # themselves H-bonded: in an antiparallel sheet half the bridges are
        # "wide" pairs whose own N and O never touch, and their bridge is carried
        # entirely by the H-bonds of their neighbours.  Testing only H-bonded
        # pairs finds exactly the narrow pairs and produces the tell-tale
        # alternating E-x-E-x string.  So every peptide position is paired with
        # every residue whose CA is within bridge_ca_max_A, and the two-H-bond
        # clause then decides.  The distance is a pre-filter for speed only: a
        # reference set of H-bonded backbone pairs never exceeds 6.79 A CA-CA, so
        # 8.0 cannot exclude a real bridge.
        cand = set()
        pep_ca = [(k, r.bb["CA"]) for k, r in enumerate(self.chain) if "CA" in r.bb]
        all_ca = [(("P", k), r.bb["CA"]) for k, r in enumerate(self.chain)
                  if "CA" in r.bb]
        all_ca += [(("R", r.key), r.bb["CA"]) for r in self.rec if "CA" in r.bb]
        if pep_ca and all_ca:
            t = cKDTree(np.array([p for _, p in all_ca]))
            for k, p in pep_ca:
                for m in t.query_ball_point(p, THRESH["bridge_ca_max_A"]):
                    tag = all_ca[m][0]
                    if tag[0] == "P":
                        if self._ring_sep(k, tag[1]) < 3:
                            continue           # too close along the backbone to
                            # be a bridge: that is a turn, not a sheet
                        cand.add(tuple(sorted((self.chain[k].key,
                                               self.chain[tag[1]].key))))
                    else:
                        cand.add((self.chain[k].key, tag[1]))
        cand = sorted(cand)
        out = {}
        seen = set()
        for a, b in cand:
            if (a, b) in seen or (b, a) in seen:
                continue
            seen.add((a, b))
            am, ap = nb(a, -1), nb(a, 1)
            bm, bp = nb(b, -1), nb(b, 1)
            kind = None
            # CORRECTION 2026-08-03.  Both wide clauses used to be written with the
            # two H-bond tuples REVERSED -- they asked for
            #     parallel      [(am, b) and (b, ap)] or [(bm, a) and (a, bp)]
            #     antiparallel  ... or [(am, bp) and (bm, ap)]
            # which is the mirror image of Kabsch & Sander.  The cause was a single
            # convention inversion, not two typos: K&S Hbond(i, j) has i ACCEPTING
            # and j DONATING, while HB here is (donor, acceptor), so every K&S tuple
            # must be reversed on the way in.  The old docstring stated the clauses
            # in K&S notation and then declared HB to be (donor, acceptor) in the
            # same sentence, and the code followed the notation instead of the set.
            #
            # Measured against mkdssp 4.2.2's own declared bridges (199 reference set
            # entries, 10428 bridges, partners typed from the ladder-label case, so
            # the type is never inferred from the bonds under test): the old
            # parallel form matched 0 of 946 declared parallel bridges and invented
            # 133; the old antiparallel-wide form missed 2443 of 4721, failing at
            # every ladder end, bulge and irregular pair, and only succeeded inside
            # regular ladders where the flanking narrow pair is reciprocally bonded
            # so both tuple orders happen to exist.
            #
            # SCOPE OF THE 10428/10428 FIGURE -- qualified 2026-08-04 after audit.
            # Fed mkdssp's OWN H-bond table, these clause forms reproduce every bridge
            # mkdssp declares: 10428/10428 (an independent re-run got 10654/10654), zero
            # missed, zero wrong-sense, against 0/946 parallel recall for the old mirrored
            # form.  That is a statement about THE CLAUSE IN ISOLATION and nothing more.
            #
            # It is NOT a statement about this tool.  End-to-end, using pep_ss's own H-bond
            # criterion instead of mkdssp's, bridge agreement over 200 entries is
            # 180/212 = 0.849: 32 declared bridges missed, 7 invented.  Attribution of the
            # 32, measured: 0 are recovered by dropping the donor-availability check and
            # only 5 by dropping the 120 deg angle gate; the remaining 27 need H-bonds
            # beyond the 3.5 A N-O cutoff that DSSP accepts on energy alone (3GHE
            # ALA319->ARG306 is a near-miss at 113.6 deg; 4F27 GLU525:A->SER12:Q is at
            # DSSP -0.7 and is not even a distance candidate).
            #
            # So: the clause is now exactly right, and the residual beta gap against DSSP
            # is the H-BOND CRITERION, which is a separate and still-open question.  The
            # unqualified original sentence read as a tool-level claim it could not support
            # -- a failure mode where a benchmark claim outran its evidence, so corrected in place rather
            # than reworded away.
            #
            # Demonstrations, straight out of the mkdssp output, N-H(donor) -> O=C(acceptor):
            #   antiparallel wide, 1BE9 bridge Q B6 -- I A328, ladder 'A':
            #       N-H(T B7)   -> O=C(I A327)  -1.8   = (ap, bm)
            #       N-H(G A329) -> O=C(K B5)    -0.6   = (bp, am)
            #     the reversed pairs (am, bp) and (bm, ap) are absent from the file.
            #   parallel, 6XYX bridge T C1342 -- S A7, ladder 'a':
            #       N-H(T C1342) -> O=C(D A6)   -2.0   = (a, bm)
            #       N-H(Q A8)    -> O=C(T C1342) -2.5  = (bp, a)
            # DO NOT re-invert these.  _bridge_selftest() in section 10 pins them on
            # ideal two-strand sheets of both senses.
            if ((b, am) in HB and (ap, b) in HB) or ((a, bm) in HB and (bp, a) in HB):
                kind = "parallel"
            elif ((a, b) in HB and (b, a) in HB) or \
                 ((bp, am) in HB and (ap, bm) in HB):
                kind = "antiparallel"
            if kind is None:
                continue
            for x, y in ((a, b), (b, a)):
                sx = pos_of.get(x)
                sy = pos_of.get(y)
                if sx is None or sy is None or sx[0] != "P":
                    continue
                out.setdefault(sx[1], []).append(
                    {"partner": label_of.get(y, y), "partner_in": sy[0],
                     "kind": kind})
        return out

    def _turns(self, nturn):
        """Typed beta-turns and gamma-turns.  Wrap-aware."""
        n = len(self.chain)
        turns, gammas = {}, {}
        tol = THRESH["turn_tol_deg"]
        loose = THRESH["turn_tol_loose_one_deg"]
        for i in range(n):
            js = [self._nbr(i, s) for s in (1, 2, 3)]
            if any(j is None for j in js):
                continue
            i1, i2, i3 = js
            a, b = self.chain[i1], self.chain[i2]
            if a.phi is None or a.psi is None or b.phi is None or b.psi is None:
                continue
            ca13 = dist(self.chain[i].bb.get("CA"), self.chain[i3].bb.get("CA"))
            if ca13 is None or ca13 > THRESH["turn_ca13_max_A"]:
                continue
            devs_by_type = {}
            for name, ((p1, s1), (p2, s2)) in TURN_TEMPLATES.items():
                dd = [dev(a.phi, p1), dev(a.psi, s1), dev(b.phi, p2), dev(b.psi, s2)]
                devs_by_type[name] = dd
            strict = [t for t, dd in devs_by_type.items() if max(dd) <= tol]
            loosef = [t for t, dd in devs_by_type.items()
                      if sorted(dd)[-1] <= loose and sorted(dd)[-2] <= tol]
            pick = strict or loosef
            if not pick:
                continue
            best = min(pick, key=lambda t: sum(devs_by_type[t]))
            rec = {"type": best, "strict": best in strict,
                   "ca_i_i3_A": round(ca13, 3),
                   "i_i3_hbond": (i, 3) in nturn,
                   "max_torsion_deviation_deg": round(max(devs_by_type[best]), 1),
                   "candidates_strict": sorted(strict),
                   "candidates_loose_one_45": sorted(loosef),
                   "positions": [i, i1, i2, i3]}
            for p in (i1, i2):
                if p not in turns or sum(devs_by_type[best]) < \
                        sum(devs_by_type[turns[p]["type"]]):
                    turns[p] = rec
        for i in range(n):
            j1, j2 = self._nbr(i, 1), self._nbr(i, 2)
            if j1 is None or j2 is None:
                continue
            if (i, 2) not in nturn:
                continue
            a = self.chain[j1]
            if a.phi is None or a.psi is None:
                continue
            for name, (p, s) in GAMMA_TEMPLATES.items():
                if dev(a.phi, p) <= THRESH["gamma_tol_deg"] and \
                        dev(a.psi, s) <= THRESH["gamma_tol_deg"]:
                    gammas[j1] = {"type": name, "positions": [i, j1, j2],
                                  "i_i2_hbond": True,
                                  "deviation_deg": round(max(dev(a.phi, p),
                                                             dev(a.psi, s)), 1)}
                    break
        return turns, gammas

    def _ca_span(self, k):
        """CA(i-1)-CA(i+1), the direct per-residue measure of how extended the
        backbone is at position i.  Wrap-aware; None at a break.

        Derived ideal values: 3-10 helix 5.13, alpha-helix 5.48, PPII 6.54,
        parallel beta 6.56, antiparallel beta 6.94 A.  So this one distance
        separates helical from extended with ~1 A of margin, and it does so
        without any reference to phi/psi, which is what makes it usable on a
        residue whose torsions are unavailable.
        """
        a, b = self._nbr(k, -1), self._nbr(k, 1)
        if a is None or b is None:
            return None
        return dist(self.chain[a].bb.get("CA"), self.chain[b].bb.get("CA"))

    def _extended_cut(self):
        """The helical/extended boundary on CA(i-1)-CA(i+1), derived as the
        midpoint of the ideal alpha-helix and ideal parallel-beta values.  No
        number is chosen by hand."""
        R = refs()
        return 0.5 * (R["alpha_R"]["ca_i_i2_A"] + R["beta_par"]["ca_i_i2_A"])

    def _ppii(self):
        """PPII: runs of >= ppii_min_run residues that are locally EXTENDED and
        whose (phi,psi) sit in one of the two mirror-image PPII boxes.

        No residue name is involved: PPII does not require proline.  Chirality
        selects which mirror box is tested and the result records which matched.
        A run must be uniform in handedness -- a stretch that alternates P and p
        is not one PPII helix and is not reported as one.

        PPII and extended beta are NOT separable on backbone torsion or on CA
        distance: the derived ideal CA(i-1)-CA(i+1) is 6.54 for PPII and 6.56 for
        parallel beta.  The discriminator used here is the H-bond pattern -- a
        beta strand has a bridge partner, PPII has none -- applied by the
        precedence order in assign(), where E/e outrank P/p.
        """
        R = refs()
        hw = R["_halfwidths"]["rama_halfwidth_deg"]
        cut = self._extended_cut()
        hits = {}
        for k, r in enumerate(self.chain):
            if r.phi is None or r.psi is None:
                continue
            span = self._ca_span(k)
            if span is None or span < cut:
                continue
            m = None
            if dev(r.phi, REF_TORSIONS["ppii"][0]) <= hw and \
               dev(r.psi, REF_TORSIONS["ppii"][1]) <= hw:
                m = False
            elif dev(r.phi, REF_TORSIONS["ppii_mirror"][0]) <= hw and \
                    dev(r.psi, REF_TORSIONS["ppii_mirror"][1]) <= hw:
                m = True
            if m is None:
                continue
            hits[k] = {"mirror": m, "ca_im1_ip1_A": round(span, 3),
                       "extended_cut_A": round(cut, 3)}
        keep = {}
        for k in list(hits):
            run, cur = [k], k
            while True:
                j = self._nbr(cur, 1)
                if j is None or j in run or j not in hits:
                    break
                if hits[j]["mirror"] != hits[k]["mirror"]:
                    break
                run.append(j)
                cur = j
            if len(run) >= THRESH["ppii_min_run"]:
                for p in run:
                    if p not in keep or len(run) > keep[p]["run_length"]:
                        keep[p] = dict(hits[p], run_length=len(run))
        return keep

    def _extended(self):
        """Locally extended backbone, measured, with no torsion box involved.

        A residue is extended when CA(i-1)-CA(i+1) is at or above the derived
        helical/extended boundary.  Whether it is then called E, e, P, p or x is
        decided by the H-bond pattern in assign(): E/e when it has a bridge
        partner, P/p when it is part of a PPII run, x when it is neither -- an
        extended residue that pairs with nothing.
        """
        cut = self._extended_cut()
        out = {}
        for k in range(len(self.chain)):
            span = self._ca_span(k)
            if span is not None and span >= cut:
                out[k] = {"ca_im1_ip1_A": round(span, 3),
                          "extended_cut_A": round(cut, 3),
                          "nearest_rama_reference":
                              (rama_class(*self.chain[k].rama_used) or None)}
        return out

    # ---------------------------------------------------------------- report
    def report(self):
        R = refs()
        modes = {m: self.assign(m) for m in ("bound", "intrinsic", "nodonor")}
        n = len(self.chain)
        per = []
        for k, r in enumerate(self.chain):
            ax = modes["bound"]["detail"][k].get("axis")
            per.append({
                "pos": k,
                "residue": r.name, "seqid": r.seqid, "icode": r.icode,
                "auth_chain": r.auth_chain,
                "ss_bound": modes["bound"]["string"][k],
                "ss_intrinsic": modes["intrinsic"]["string"][k],
                "ss_nodonor": modes["nodonor"]["string"][k],
                "receptor_induced": (modes["bound"]["string"][k]
                                     != modes["intrinsic"]["string"][k]),
                "backbone_alpha": r.alpha,
                "backbone_detail": r.bt,
                "donor_available": r.donor.get("available"),
                "donor_reason": r.donor.get("reason"),
                "donor_n_heavy_neighbours": r.donor.get("n_heavy_neighbours"),
                "chirality": r.chir.get("label"),
                "chirality_improper_deg": r.chir.get("improper_deg"),
                "chirality_reason": r.chir.get("reason"),
                "phi": None if r.phi is None else round(r.phi, 1),
                "psi": None if r.psi is None else round(r.psi, 1),
                "omega": None if r.omega is None else round(r.omega, 1),
                "phi_psi_used_for_reference_match": (
                    None if r.rama_used[0] is None
                    else [round(r.rama_used[0], 1), round(r.rama_used[1], 1)]),
                "mirrored_for_D": r.rama_mirrored,
                "rama_reference_hits": rama_class(*r.rama_used),
                "ca_im1_ip1_A": (round(self._ca_span(k), 3)
                                 if self._ca_span(k) is not None else None),
                "ca_i_i2_A": self._cd(k, 2), "ca_i_i3_A": self._cd(k, 3),
                "ca_i_i4_A": self._cd(k, 4),
                "local_axis": ax,
                "flags": r.flags,
                "detail_bound": {kk: vv for kk, vv in
                                 modes["bound"]["detail"][k].items()
                                 if kk != "axis"},
            })
        agg = {m: self._agg(modes[m], n) for m in modes}
        out = {
            "tool": VERSION,
            "pdb_id": self.pdb_id,
            "entry_dir": self.entry_dir,
            "topology": self.topo,
            "n_residues": n,
            "n_assignable": sum(1 for p in per if p["ss_bound"] not in "X?"),
            "n_not_assignable_incomplete_backbone":
                sum(1 for p in per if p["ss_bound"] == "X"),
            "n_refused_non_alpha_backbone":
                sum(1 for p in per if p["ss_bound"] == "?"),
            "n_donor_unavailable": sum(1 for p in per
                                       if p["donor_available"] is False),
            # recorded to make the point that it is recorded and NOT used: the
            # assignment logic never reads het_flag, so a modified residue
            # deposited as HETATM cannot create a phantom chain break here.
            "n_residues_flagged_hetatm": sum(1 for r in self.chain
                                             if r.het == "H"),
            "chirality_counts": _count(p["chirality"] for p in per),
            "hbonds": {
                "n_total": len(self.hb),
                "n_peptide_internal": sum(1 for h in self.hb
                                          if h["class"] == "peptide_internal"),
                "n_peptide_to_receptor": sum(1 for h in self.hb
                                             if h["class"] == "peptide_to_receptor"),
                "n_receptor_to_peptide": sum(1 for h in self.hb
                                             if h["class"] == "receptor_to_peptide"),
                "n_refused_no_donor_hydrogen": len(self.refused),
                "n_within_distance_but_failing_angle": len(self.near),
                "n_phantom_by_dssp_energy_criterion": len(self.phantom_dssp),
                "phantom_by_dssp_energy_criterion": self.phantom_dssp,
                "refused": self.refused,
                "near_misses": self.near,
                "list": self.hb,
            },
            "dssp_defect_measured": {
                "n_phantom_hbonds_by_this_modules_criterion": len(self.refused),
                "n_phantom_hbonds_by_dssp_energy_criterion": len(self.phantom_dssp),
                "n_residues_changed_by_ignoring_donor_availability":
                    sum(1 for k in range(n)
                        if modes["bound"]["string"][k] != modes["nodonor"]["string"][k]),
                "string_bound": modes["bound"]["string"],
                "string_if_donors_ignored": modes["nodonor"]["string"],
            },
            "per_residue": per,
            "aggregate": agg,
            "references_derived": R,
            "thresholds": {k: v for k, v in THRESH.items()},
        }
        return out

    def _cd(self, k, step):
        j = self._nbr(k, step)
        if j is None:
            return None
        d = dist(self.chain[k].bb.get("CA"), self.chain[j].bb.get("CA"))
        return None if d is None else round(d, 3)

    def _agg(self, res, n):
        s = res["string"]
        den = sum(1 for c in s if c not in "X?")
        def frac(chars):
            if den == 0:
                return None
            return round(sum(1 for c in s if c in chars) / den, 4)
        return {
            "string": s,
            "n_assignable": den,
            "frac_helix_any": frac("HhGgI"),
            "frac_helix_alpha_right": frac("H"),
            "frac_helix_alpha_left": frac("h"),
            "frac_helix_310": frac("Gg"),
            "frac_helix_pi": frac("I"),
            "frac_strand_receptor_paired": frac("E"),
            "frac_strand_peptide_internal": frac("e"),
            "frac_extended_unpaired": frac("x"),
            "frac_ppii": frac("Pp"),
            "frac_turn_beta": frac("T"),
            "frac_turn_gamma": frac("S"),
            "frac_coil": frac("C"),
            "longest_helix_run": _run(s, "HhGgI", self.wrap),
            "longest_strand_run": _run(s, "Ee", self.wrap),
            "longest_ppii_run": _run(s, "Pp", self.wrap),
            # A refused ('?') or unassignable ('X') residue interrupts a run in
            # the strict count above.  For an alpha/beta-peptide foldamer such as
            # 4BPI or 4BPJ that is misleading: the helix really does continue
            # through the beta-amino acid, we simply cannot describe that residue.
            # So the bridging count, which steps over X and ? without resetting,
            # is reported alongside -- never instead.
            "longest_helix_run_bridging_unassignable":
                _run(s, "HhGgI", self.wrap, bridge="X?"),
            "longest_strand_run_bridging_unassignable":
                _run(s, "Ee", self.wrap, bridge="X?"),
            "category": _category(s, self.wrap),
            "category_bridging_unassignable": _category(s, self.wrap, bridge="X?"),
            # counted only where the residue is actually coded T / S.  A turn
            # window can also be found under a helix, where helix precedence
            # wins; counting those would inflate the turn census.
            "turn_types": _count(t["type"] for p, t in res["turns"].items()
                                 if s[p] == "T"),
            "gamma_turn_types": _count(t["type"] for p, t in
                                       res["gamma_turns"].items() if s[p] == "S"),
        }


def _count(it):
    d = {}
    for x in it:
        if x is None:
            continue
        d[x] = d.get(x, 0) + 1
    return dict(sorted(d.items()))


def _run(s, chars, wrap, bridge=""):
    """Longest run of chars.  Wraps for a macrocycle, because a helix or a strand
    in a head-to-tail cycle has no privileged start.  Characters in `bridge` are
    stepped over without resetting the run and without being counted."""
    if not s:
        return 0
    t = s + s if wrap else s
    best = cur = 0
    for c in t:
        if c in chars:
            cur += 1
            best = max(best, cur)
        elif c in bridge:
            continue
        else:
            cur = 0
    return min(best, len(s))


def _category(s, wrap, bridge=""):
    """Coarse per-peptide class, the same shape ss_assign.py uses so the two can
    be compared."""
    if _run(s, "HhGgI", wrap, bridge) >= THRESH["helix_min_run"]:
        return "helix"
    if _run(s, "Ee", wrap, bridge) >= THRESH["strand_min_run"]:
        return "sheet"
    if _run(s, "Pp", wrap, bridge) >= THRESH["ppii_min_run"]:
        return "ppii"
    if "T" in s or "S" in s:
        return "turn_only"
    if "x" in s:
        return "extended_unpaired"
    return "none_or_coil"


# --------------------------------------------------------------------------- #
# 10. SELF TEST / CALIBRATION ON IDEAL GEOMETRY
# --------------------------------------------------------------------------- #
class _IdealAssignment(Assignment):
    """Assignment over a synthetic ideal chain; no files, no receptor."""

    def __init__(self, residues, name="ideal", chir="L"):
        self.entry_dir = "(ideal)"
        self.pdb_id = name
        self.pep = []
        for i, r in enumerate(residues):
            fake = _FakeRes(r, i + 1, chir)
            self.pep.append(fake)
        self.other = []
        self.notes = []
        self._prepare()


class _FakeRes(Res):
    def __init__(self, bb, seqid, chir="L", chain="Z", is_peptide=True):
        self.name = "ALA"
        self.seqid = seqid
        self.icode = ""
        self.auth_chain = chain
        self.subchain = chain
        self.het = "A"
        # CB is BUILT from the measured L geometry -- CA-CB 1.532 A, N-CA-CB
        # 110.69 deg, dihedral(C,N,CA,CB) -122.59 deg -- so that the chirality
        # test in selftest() is a real round trip: build L, read L back.  chir="D"
        # mirrors only that dihedral, which is exactly what inverting the centre
        # means.
        N, CA, C = bb["N"], bb["CA"], bb["C"]
        sign = 1.0 if chir == "L" else -1.0
        CB = _place(C, N, CA, IDEAL["CA_CB"], IDEAL["N_CA_CB"],
                    sign * IDEAL["CB_DIH_L"])
        self.atoms = {"N": (N, "N"), "CA": (CA, "C"), "C": (C, "C"),
                      "O": (bb["O"], "O"), "CB": (CB, "C")}
        self.bb = {k: bb[k] for k in BB}
        self.is_peptide = is_peptide
        self.idx = None
        self.path_bonds = None
        self.n_neighbours = []
        self.donor = None
        self.chir = None
        self.phi = self.psi = self.omega = None
        self.refuse = None
        self.flags = []


def build_ideal_cycle(n=6, phi0=-70.0, psi0=130.0, seed=0):
    """Build a genuinely ring-closed ideal cyclic peptide.

    The torsions are relaxed by least squares until the amide that would be built
    from residue n-1's psi lands on residue 0's N and CA -- i.e. until the ring
    actually closes with the reference-measured amide geometry.  Bond lengths and
    angles are never adjusted; only phi and psi move.  Returns (residues, rms
    closure error in A) or (None, reason).
    """
    from scipy.optimize import least_squares
    rng = np.random.default_rng(seed)
    x0 = np.concatenate([np.full(n, phi0), np.full(n, psi0)]) \
        + rng.normal(0.0, 15.0, 2 * n)

    def resid(x):
        phi, psi = x[:n], x[n:]
        r = build_ideal(n + 2, np.concatenate([phi, phi[:2]]),
                        np.concatenate([psi, psi[:2]]))
        # residue n of that build is where residue 0 must sit
        return np.concatenate([r[n]["N"] - r[0]["N"], r[n]["CA"] - r[0]["CA"]])

    # the residual compares absolute positions, so first superpose: instead build
    # the open chain and compare the closing amide to the FIRST residue after
    # rigid-body alignment of residue n onto residue 0 is impossible -- so use an
    # internal-coordinate residual instead: the closing C-N length, the two
    # closing angles and the closing omega, all built from the open chain.
    def resid2(x):
        phi, psi = x[:n], x[n:]
        r = build_ideal(n, phi, psi)
        A, B = r[-1], r[0]
        return np.array([
            dist(A["C"], B["N"]) - IDEAL["C_N"],
            (angle_deg(A["CA"], A["C"], B["N"]) - IDEAL["CA_C_N"]) / 30.0,
            (angle_deg(A["C"], B["N"], B["CA"]) - IDEAL["C_N_CA"]) / 30.0,
            (dev(dihedral(A["CA"], A["C"], B["N"], B["CA"]) or 0.0, 180.0)) / 30.0,
        ])
    sol = least_squares(resid2, x0, xtol=1e-12, ftol=1e-12, max_nfev=20000)
    r = build_ideal(n, sol.x[:n], sol.x[n:])
    err = float(np.sqrt((resid2(sol.x) ** 2).mean()))
    if abs(dist(r[-1]["C"], r[0]["N"]) - IDEAL["C_N"]) > 0.05:
        return None, f"closure did not converge (C-N {dist(r[-1]['C'], r[0]['N']):.3f} A)"
    return r, err


def _cyclic_selftest():
    """Wrap-logic test on a genuinely ring-closed ideal cyclic hexapeptide.

    Under test: does order_backbone see a cycle rather than a linear chain, does
    `wrap` become True, does _nbr walk across the closure, and do the first
    residue's phi and the last residue's psi stop being null.  An implementation
    that assumes a linear chain fails every one of these, and it would fail them
    silently on the macrocyclic entries.
    """
    res, err = build_ideal_cycle(6)
    if res is None:
        return None, False, f"could not build a closed ring: {err}"
    a = _IdealAssignment(res, "cyclic_ideal")
    checks = {
        "cycle_detected": a.topo["is_macrocycle_by_backbone_amide"] is True,
        "wrap_true": a.wrap is True,
        "nbr_last_plus1_is_first": a._nbr(len(a.chain) - 1, 1) == 0,
        "nbr_first_minus1_is_last": a._nbr(0, -1) == len(a.chain) - 1,
        "phi_of_first_not_null": a.chain[0].phi is not None,
        "psi_of_last_not_null": a.chain[-1].psi is not None,
        "closure_recorded": a.topo["closure_link"] is not None,
        "no_breaks": a.topo["breaks"] == [],
    }
    bad = [k for k, v in checks.items() if not v]
    cn = dist(a.chain[-1].bb["C"], a.chain[0].bb["N"])
    return a, not bad, (f"closed ring built to C-N {cn:.3f} A (target "
                        f"{IDEAL['C_N']}); all {len(checks)} graph checks pass "
                        "(cycle detected, wrap on, neighbour walk crosses the "
                        "closure, terminal torsions computed through it)"
                        if not bad else f"C-N {cn:.3f} A; failed: " + ", ".join(bad))


def build_ideal_sheet(sense, offset_A, shift_A, n=10):
    """Two ideal beta strands in sheet register, from ONE ideal strand plus a rigid
    transform.  Returns (strand_A, strand_B) as lists of N/CA/C/O dicts.

    Frame of the reference strand: `ax` along the strand (CA[0] -> CA[-1]); `v` the
    in-plane direction an even-index residue's C=O points, taken perpendicular to
    ax; `w = ax x v` the sheet normal.  A PARALLEL partner is the same strand
    translated by offset_A along v and shift_A along ax.  An ANTIPARALLEL partner is
    additionally rotated 180 deg about w through the strand centroid, which reverses
    the chain direction and the carbonyl direction together and so turns the C=O
    groups back towards the reference strand.  Nothing is fitted or relaxed: the
    transform is closed-form, so the result is reproducible to the bit.

    Only alternate residues of a strand point their N-H and C=O at a given
    neighbour -- that is the beta pleat -- so the bridges come out alternating,
    exactly as they do in a real ladder.
    """
    A = build_ideal(n, *REF_TORSIONS["beta_anti"])
    ca = np.array([r["CA"] for r in A])
    ax = _u(ca[-1] - ca[0])
    k = 2 * ((n // 2) // 2)                    # an even-index interior residue
    d = A[k]["O"] - A[k]["C"]
    v = _u(d - np.dot(d, ax) * ax)
    w = np.cross(ax, v)
    R = (2.0 * np.outer(w, w) - np.eye(3)) if sense == "antiparallel" else np.eye(3)
    cen = ca.mean(axis=0)
    t = cen - R @ cen + offset_A * v + shift_A * ax
    B = [{q: R @ p + t for q, p in r.items()} for r in A]
    return A, B


class _IdealSheet(Assignment):
    """Ideal two-strand sheet: strand A is the peptide, strand B the receptor."""

    def __init__(self, sense, offset_A, shift_A, n=10):
        self.entry_dir = "(ideal)"
        self.pdb_id = "sheet_" + sense
        self.notes = []
        A, B = build_ideal_sheet(sense, offset_A, shift_A, n)
        self.pep = [_FakeRes(r, i + 1, chain="P") for i, r in enumerate(A)]
        self.other = [_FakeRes(r, i + 1, chain="R", is_peptide=False)
                      for i, r in enumerate(B)]
        self._prepare()


# Ideal-sheet geometries used by _bridge_selftest.  offset/shift are in A along the
# frame described in build_ideal_sheet.  `want_nE` and the reference string were
# MEASURED on 2026-08-03 with the corrected clause; the third column records what
# the pre-correction mirror-image clause produced on the same coordinates, which is
# what makes each case a discriminating test rather than a smoke test.
_SHEET_CASES = (
    # sense           offset shift  nE  string        nE under the old mirror clause
    ("parallel",        4.6,  0.0,   8, "CEEEEEEEEC", 0),
    ("antiparallel",    4.2,  1.6,   5, "CxxEEEEExC", 3),
    ("antiparallel",    4.0,  2.0,   5, "CxxEEEEExC", 3),
)


def _bridge_selftest():
    """Pin the DSSP bridge clause (_bridges) on ideal sheets of both senses.

    THIS IS THE TEST THAT DID NOT EXIST when the (donor, acceptor) versus
    Kabsch-Sander convention inversion corrected on 2026-08-03 went in, and its
    absence is why the mirror-image clause survived to four documents.  Measured on
    the same coordinates, the mirror clause loses every bridge of the parallel sheet
    (8 E residues -> 0, since it has no recall at all on parallel bridges) and loses
    the two WIDE bridges at the ends of the antiparallel ladder (5 -> 3), which is
    exactly the failure mode seen against mkdssp on the reference set.

    Returns a list of (ok, description) pairs.
    """
    out = []
    for sense, off, sh, want_nE, want_s, old_nE in _SHEET_CASES:
        a = _IdealSheet(sense, off, sh)
        r = a.assign("bound")
        s = r["string"]
        nE = sum(s.count(c) for c in "Ee")
        kinds = sorted({bp["kind"] for bl in r["bridges"].values() for bp in bl})
        ok = nE >= want_nE and kinds == [sense]
        out.append((ok, f"ideal {sense:12s} sheet (offset {off} A, shift {sh} A): "
                        f"{nE} E residues (want >={want_nE}, mirror clause gives "
                        f"{old_nE}), kinds {kinds or ['none']}, string {s} "
                        f"(reference {want_s})"))
    return out


def selftest(verbose=True):
    """Ideal-geometry calibration.  Returns (n_fail, lines)."""
    L, fails = [], 0
    R = refs()
    L.append("DERIVED REFERENCE TABLE (from ideal chains built with reference-measured "
             "bond geometry)")
    L.append(f"{'conformation':13s} {'phi':>7s} {'psi':>7s} {'i,i+1':>6s} {'i,i+2':>6s} "
             f"{'i,i+3':>6s} {'i,i+4':>6s} {'rise':>6s} {'radius':>7s} "
             f"{'twist':>7s} {'n/turn':>7s} {'hand':>5s}  n-turns(3/4/5)")
    for name in REF_TORSIONS:
        r = R[name]
        a = r["axis"]
        f = (lambda k, nd=2: ("%.*f" % (nd, a[k])) if a.get("ok") else "-")
        L.append(f"{name:13s} {r['phi']:7.1f} {r['psi']:7.1f} {r['ca_i_i1_A']:6.2f} "
                 f"{r['ca_i_i2_A']:6.2f} {r['ca_i_i3_A']:6.2f} {r['ca_i_i4_A']:6.2f} "
                 f"{f('rise_per_residue_A'):>6s} {f('radius_A'):>7s} "
                 f"{f('twist_per_residue_deg',1):>7s} {f('residues_per_turn',2):>7s} "
                 f"{(a.get('handedness') or '-'):>5s}  "
                 f"{r['hbonds']['n3_turns']}/{r['hbonds']['n4_turns']}/"
                 f"{r['hbonds']['n5_turns']}"
                 + ("" if a.get("ok") else f"   [axis: {a.get('reason')}]"))
    L.append("")
    L.append("TURN TEMPLATE CA(i)-CA(i+3), built from the template torsions:")
    for name in TURN_TEMPLATES:
        L.append(f"  {name:5s} {R['turn_' + name]['ca_i_i3_A']:5.2f} A")
    L.append(f"  gate = {THRESH['turn_ca13_max_A']} A;  ideal PPII trimer "
             f"{R['ppii']['ca_i_i3_A']:.2f} A, ideal antiparallel strand "
             f"{R['beta_anti']['ca_i_i3_A']:.2f} A")
    L.append("")

    L.append("RECOVERY OF IDEAL CONFORMATIONS BY THE FULL ASSIGNER")
    expect = {
        "alpha_R": ("H", 8), "alpha_L": ("h", 8), "helix310": ("G", 6),
        "beta_anti": ("x", 8), "beta_par": ("x", 8), "ppii": ("P", 8),
        "ppii_mirror": ("p", 8),
    }
    for name, (want, minimum) in expect.items():
        phi, psi = REF_TORSIONS[name]
        a = _IdealAssignment(build_ideal(14, phi, psi), name)
        s = a.assign("intrinsic")["string"]
        got = s.count(want)
        ok = got >= minimum
        fails += 0 if ok else 1
        L.append(f"  {'PASS' if ok else 'FAIL'} {name:13s} want >={minimum} '{want}'"
                 f"  got {got}  string {s}")
    # pi helix is marked I by the same rule; report it but do not fail on it,
    # because an ideal pi helix also satisfies the 4-turn criterion in places
    a = _IdealAssignment(build_ideal(14, *REF_TORSIONS["helix_pi"]), "helix_pi")
    L.append(f"  INFO helix_pi     string {a.assign('intrinsic')['string']}")

    # turn typology recovery: build each template turn with extended flanks and
    # check the assigner returns that type and only that type
    for name, ((p1, s1), (p2, s2)) in TURN_TEMPLATES.items():
        res = build_ideal(6, [-120.0, -120.0, p1, p2, -120.0, -120.0],
                          [130.0, 130.0, s1, s2, 130.0, 130.0])
        a = _IdealAssignment(res, "turn" + name)
        r = a.assign("intrinsic")
        got = {t["type"] for p, t in r["turns"].items() if r["string"][p] == "T"}
        ok = got == {name}
        fails += 0 if ok else 1
        L.append(f"  {'PASS' if ok else 'FAIL'} turn type {name:5s} recovered as "
                 f"{sorted(got) or ['nothing']}  string {r['string']}")
    for name, (p, s) in GAMMA_TEMPLATES.items():
        res = build_ideal(5, [-120.0, -120.0, p, -120.0, -120.0],
                          [130.0, 130.0, s, 130.0, 130.0])
        a = _IdealAssignment(res, name)
        r = a.assign("intrinsic")
        got = {t["type"] for t in r["gamma_turns"].values()}
        ok = got == {name} and "S" in r["string"]
        fails += 0 if ok else 1
        L.append(f"  {'PASS' if ok else 'FAIL'} {name} recovered as "
                 f"{sorted(got) or ['nothing']}  string {r['string']}")

    # chirality: build L, read L back; build the mirror, read D back
    for want in ("L", "D"):
        a = _IdealAssignment(build_ideal(8, *REF_TORSIONS["alpha_R"]), "chir", want)
        labs = {r.chir.get("label") for r in a.chain}
        imp = [r.chir.get("improper_deg") for r in a.chain]
        ok = labs == {want}
        fails += 0 if ok else 1
        L.append(f"  {'PASS' if ok else 'FAIL'} chirality round trip: built {want}, "
                 f"read {labs}, improper {min(imp):.2f}..{max(imp):.2f} deg "
                 f"(reference set L residues measure -32.94 +- 1.59)")
    # a quaternary CA (Aib-like: a second substituent) must be refused, not guessed
    res = build_ideal(6, *REF_TORSIONS["alpha_R"])
    a = _IdealAssignment(res, "aib")
    r = a.chain[2]
    r.atoms["CB2"] = (_place(r.bb["C"], r.bb["N"], r.bb["CA"], IDEAL["CA_CB"],
                             IDEAL["N_CA_CB"], -IDEAL["CB_DIH_L"]), "C")
    a._prepare()
    c = a.chain[2].chir
    ok = c["label"] is None and "quaternary" in (c["reason"] or "")
    fails += 0 if ok else 1
    L.append(f"  {'PASS' if ok else 'FAIL'} quaternary CA (Aib-like) refuses a "
             f"chirality call: label={c['label']} reason={c['reason']}")

    # open chain: the neighbour walk must NOT wrap
    a = _IdealAssignment(build_ideal(6, -70.0, 130.0), "openchain")
    ok = (a.wrap is False and a._nbr(5, 1) is None and a._nbr(0, -1) is None)
    fails += 0 if ok else 1
    L.append(f"  {'PASS' if ok else 'FAIL'} open chain: wrap={a.wrap}, "
             f"nbr(last,+1)={a._nbr(5, 1)}, nbr(first,-1)={a._nbr(0, -1)}")

    # head-to-tail macrocycle: build an ideal chain, then MOVE the last residue so
    # its C sits at amide distance from the first residue's N, and check that the
    # walk wraps and that residue 0 gets a phi and residue n-1 a psi.
    a, ok, why = _cyclic_selftest()
    fails += 0 if ok else 1
    L.append(f"  {'PASS' if ok else 'FAIL'} head-to-tail macrocycle: {why}")

    # bridge clause on ideal two-strand sheets of both senses.  Added 2026-08-03
    # with the convention-inversion correction in _bridges; see _bridge_selftest.
    for ok, why in _bridge_selftest():
        fails += 0 if ok else 1
        L.append(f"  {'PASS' if ok else 'FAIL'} bridge clause, {why}")

    # donor availability on the ideal chain: residue 0 is a free amine (1 heavy
    # neighbour), the rest are secondary amides (2)
    a = _IdealAssignment(build_ideal(6, *REF_TORSIONS["alpha_R"]), "donor")
    got = [r.donor["n_heavy_neighbours"] for r in a.chain]
    ok = got == [1, 2, 2, 2, 2, 2]
    fails += 0 if ok else 1
    L.append(f"  {'PASS' if ok else 'FAIL'} donor neighbour counts on an ideal open "
             f"chain: {got} (expect [1,2,2,2,2,2])")
    # N-methylation.  An ideal 10-residue alpha helix, then a methyl group on the
    # backbone N of residues 5 and 6 (0-based), which is where the amide H of the
    # i->i+4 bonds from residues 1 and 2 would have to be.  The donor must be
    # refused, and the helix must break -- which is precisely what DSSP cannot do.
    res = build_ideal(10, *REF_TORSIONS["alpha_R"])
    a = _IdealAssignment(res, "nme")
    for k in (5, 6):
        r = a.chain[k]
        prevC = a.chain[k - 1].bb["C"]
        CM = _place(a.chain[k - 1].bb["CA"], prevC, r.bb["N"], 1.47,
                    IDEAL["C_N_CA"], 180.0)
        r.atoms["CM"] = (CM, "C")
    a._prepare()
    ds = [a.chain[k].donor for k in (5, 6)]
    ok = all(d["available"] is False and d["n_heavy_neighbours"] == 3 for d in ds)
    fails += 0 if ok else 1
    L.append(f"  {'PASS' if ok else 'FAIL'} N-methyl on residues 6-7 refuses the "
             f"donor: available={[d['available'] for d in ds]}, "
             f"neighbours={[d['n_heavy_neighbours'] for d in ds]}")
    s_yes = a.assign("intrinsic")["string"]
    s_no = a.assign("nodonor")["string"]
    ok = s_yes != s_no
    fails += 0 if ok else 1
    L.append(f"  {'PASS' if ok else 'FAIL'} the same coordinates give "
             f"{s_yes!r} with the donor check and {s_no!r} without it "
             f"({sum(1 for x, y in zip(s_yes, s_no) if x != y)} residues differ) "
             "-- that difference IS the DSSP defect, in isolation")
    nph = len(a.refused)
    L.append(f"       {nph} H-bond(s) pass the distance and angle test but have no "
             f"donor hydrogen; DSSP would score every one of them")

    # beta-amino-acid refusal.  Built, not borrowed: insert a real methylene
    # between N and CA so the covalent path N->C is three bonds and N-CA is 2.4 A,
    # which is exactly the 4BPJ/4BPI HR7 situation the build pass got wrong.
    res = build_ideal(6, *REF_TORSIONS["alpha_R"])
    a = _IdealAssignment(res, "beta")
    r = a.chain[2]
    N, CA, C, O = r.bb["N"], r.bb["CA"], r.bb["C"], r.bb["O"]
    ax = _u(CA - N)
    perp = _u(np.cross(ax, C - CA))

    def rot(v, k, deg):
        """Rodrigues rotation of v about unit axis k."""
        t = math.radians(deg)
        return (v * math.cos(t) + np.cross(k, v) * math.sin(t)
                + k * np.dot(k, v) * (1 - math.cos(t)))
    u1 = ax
    CX = N + 1.46 * u1                       # the extra backbone methylene
    u2 = rot(u1, perp, 70.0)
    CA2 = CX + 1.53 * u2                     # the alpha carbon, now 2 bonds from N
    u3 = rot(u2, perp, -70.0)
    C2 = CA2 + 1.52 * u3
    r.atoms = {"N": (N, "N"), "CX": (CX, "C"), "CA": (CA2, "C"),
               "C": (C2, "C"), "O": (C2 + (O - C), "O")}
    r.bb["CA"], r.bb["C"], r.bb["O"] = CA2, C2, C2 + (O - C)
    ok_b, det = backbone_type(r)
    L.append(f"  {'FAIL' if ok_b else 'PASS'} synthetic beta backbone refused: "
             f"alpha={ok_b}  path={det['path']} ({det['n_to_c_path_bonds']} bonds)  "
             f"N-CA {det['n_ca_A']} A  CA-C {det['ca_c_A']} A")
    if ok_b:
        fails += 1
    else:
        L.append(f"       reason: {det['reason']}")
        c = chirality(r, ok_b, det)
        ok_c = c["label"] is None
        fails += 0 if ok_c else 1
        L.append(f"  {'PASS' if ok_c else 'FAIL'} and no chirality is guessed for "
                 f"it: label={c['label']}  reason={c['reason'][:90]}")

    if verbose:
        print("\n".join(L))
    return fails, L


# --------------------------------------------------------------------------- #
# 11. CLI
# --------------------------------------------------------------------------- #
def _fmt_entry(rep):
    o = []
    t = rep["topology"]
    o.append(f"{rep['pdb_id']}  {rep['n_residues']} residues, "
             f"{rep['n_assignable']} assignable")
    o.append(f"  topology: segments={t['n_segments']}{t['segment_lengths']}  "
             f"macrocycle_backbone={t['is_macrocycle_by_backbone_amide']}  "
             f"wrap={t['neighbour_indices_wrap']}  "
             f"any_ring={t['is_macrocycle_any_ring']}")
    if t["closure_link"]:
        c = t["closure_link"]
        o.append(f"    closure: C of {c['from_C_of']} -> N of {c['to_N_of']} "
                 f"at {c['c_n_A']} A, spanning {c['span_residues']} residues")
    for b in t["breaks"]:
        o.append(f"    break after {b['after']} before {b['before']}: {b['kind']}, "
                 f"{b['n_unmodelled_residues_by_numbering']} unmodelled by numbering, "
                 f"CA-CA {b['ca_ca_A']} A")
    for c in t["crosslinks"]:
        if not c["backbone_neighbours"]:
            o.append(f"    crosslink {c['atom_a']} of {c['res_a']} - "
                     f"{c['atom_b']} of {c['res_b']}  {c['dist_A']} A ({c['elements']})")
    seq = " ".join(f"{p['residue']}{p['seqid']}" for p in rep["per_residue"])
    o.append(f"  sequence (backbone order): {seq}")
    o.append(f"  ss bound      {rep['aggregate']['bound']['string']}")
    o.append(f"  ss intrinsic  {rep['aggregate']['intrinsic']['string']}")
    o.append(f"  ss if donors ignored (DSSP's assumption)  "
             f"{rep['dssp_defect_measured']['string_if_donors_ignored']}")
    o.append(f"  category bound={rep['aggregate']['bound']['category']}  "
             f"intrinsic={rep['aggregate']['intrinsic']['category']}")
    h = rep["hbonds"]
    o.append(f"  H-bonds: internal {h['n_peptide_internal']}, "
             f"peptide->receptor {h['n_peptide_to_receptor']}, "
             f"receptor->peptide {h['n_receptor_to_peptide']}, "
             f"REFUSED for lack of a donor H {h['n_refused_no_donor_hydrogen']}")
    for r in h["refused"]:
        o.append(f"    refused: {r['donor']} N-H...O {r['acceptor']}  "
                 f"N-O {r['N_O_A']} A  ({r['donor_reason'][:70]})")
    o.append(f"  chirality {rep['chirality_counts']}  "
             f"donor-unavailable {rep['n_donor_unavailable']}  "
             f"refused-non-alpha {rep['n_refused_non_alpha_backbone']}")
    o.append("  per residue:")
    o.append("    pos res      bnd int nod  phi     psi     omega   chir  imp    "
             "donor  span  CA(i,i+2/3/4)          rama")
    for p in rep["per_residue"]:
        def g(x, w=7, nd=1):
            return ("%*.*f" % (w, nd, x)) if x is not None else " " * (w - 1) + "-"
        o.append(f"    {p['pos']:3d} {p['residue']:>3s}{p['seqid']:<5d} "
                 f" {p['ss_bound']}   {p['ss_intrinsic']}   {p['ss_nodonor']}  "
                 f"{g(p['phi'])} {g(p['psi'])} {g(p['omega'])} "
                 f"{(p['chirality'] or '-'):>4s} {g(p['chirality_improper_deg'],6)} "
                 f"{('Y' if p['donor_available'] else 'N' if p['donor_available'] is False else '?'):>5s}  "
                 f"{g(p['ca_im1_ip1_A'], 5, 2)} "
                 f"{g(p['ca_i_i2_A'], 5, 2)} {g(p['ca_i_i3_A'], 5, 2)} "
                 f"{g(p['ca_i_i4_A'], 5, 2)}  "
                 f"{','.join(p['rama_reference_hits']) or '-'}")
        if p["ss_bound"] == "?":
            o.append(f"        REFUSED: {p['backbone_detail']['reason']}")
        if p["detail_bound"].get("turn"):
            tt = p["detail_bound"]["turn"]
            o.append(f"        turn type {tt['type']}"
                     f"{'' if tt['strict'] else ' (loose)'}  "
                     f"CA(i,i+3) {tt['ca_i_i3_A']} A  i->i+3 H-bond {tt['i_i3_hbond']}")
        if p["detail_bound"].get("bridge_partners"):
            for bp in p["detail_bound"]["bridge_partners"]:
                o.append(f"        bridge {bp['kind']} to {bp['partner']} "
                         f"({'receptor' if bp['partner_in'] == 'R' else 'peptide'})")
        ax = p["local_axis"]
        # printed only where a helix axis is a meaningful description of the
        # window; it is always present in the JSON
        if ax and ax.get("ok") and p["ss_bound"] in "HhGgIPp":
            o.append(f"        local axis rise {ax['rise_per_residue_A']} A  "
                     f"radius {ax['radius_A']} A  twist {ax['twist_per_residue_deg']} "
                     f"deg  {ax['residues_per_turn']} res/turn  {ax['handedness']}")
    return "\n".join(o)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--entry", help="PDB id under entries/, or a path to an entry dir")
    ap.add_argument("--json", help="write the full report as JSON here")
    ap.add_argument("--brief", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--refs", action="store_true")
    a = ap.parse_args(argv)
    if a.refs:
        print(json.dumps(refs(), indent=1, default=str))
        return 0
    if a.selftest:
        f, _ = selftest()
        print(f"\n{'SELFTEST OK' if f == 0 else f'SELFTEST FAILED: {f} checks'}")
        return 1 if f else 0
    if not a.entry:
        ap.error("give --entry, --selftest or --refs")
    d = a.entry if os.path.isdir(a.entry) else os.path.join(REPO_ROOT, "bmi200", "entries", a.entry)
    rep = Assignment(d).report()
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(rep, fh, indent=1, default=str)
    if a.brief:
        g = rep["aggregate"]
        print(f"{rep['pdb_id']}\t{g['bound']['string']}\t{g['intrinsic']['string']}\t"
              f"{g['bound']['category']}\t{g['intrinsic']['category']}\t"
              f"refused_hb={rep['hbonds']['n_refused_no_donor_hydrogen']}")
    else:
        print(_fmt_entry(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
