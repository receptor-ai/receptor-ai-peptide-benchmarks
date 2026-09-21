import numpy as np
import biotite.structure as struc
from rdkit import Chem

DEFAULT_VDW_RADIUS = 1.7

periodic_table = Chem.GetPeriodicTable()


def _rvdw(element: str) -> float:
    try:
        return periodic_table.GetRvdw(periodic_table.GetAtomicNumber(element.capitalize()))
    except RuntimeError:
        return DEFAULT_VDW_RADIUS


get_rvdw = np.vectorize(_rvdw)

METAL_ELEMENTS = {"NA", "MG", "K", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "LI", "AU"}
DISULFIDE_MAX_DIST = 2.5
METAL_COORDINATION_MAX_DIST = 2.8


def clashscore(atoms: struc.AtomArray, overlap_tolerance: float = 0.4, severity_weighted: bool = False) -> float:
    """
    Approximate MolProbity-style clashscore: serious VDW overlaps per 1000 atoms.
    Heavy-atom only (no explicit hydrogens); bonded pairs are approximated by
    excluding same-chain atoms within 1 residue of each other, plus disulfide
    bridges and metal-ion coordination bonds (both can span distant residues).

    severity_weighted: instead of counting each clashing pair as 1, sum how much
    each pair overlaps beyond the tolerance -- a single deep clash then weighs
    more than several pairs that barely cross the threshold.
    """
    coord = atoms.coord
    n = len(atoms)
    radii = get_rvdw(atoms.element)

    dist = np.linalg.norm(coord[:, None, :] - coord[None, :, :], axis=-1)
    vdw_sum = radii[:, None] + radii[None, :]
    clash_mask = dist < (vdw_sum - overlap_tolerance)

    bonded_like = (atoms.chain_id[:, None] == atoms.chain_id[None, :]) & \
                  (np.abs(atoms.res_id[:, None] - atoms.res_id[None, :]) <= 1)

    is_sg_cys = (atoms.res_name == "CYS") & (atoms.atom_name == "SG")
    disulfide = is_sg_cys[:, None] & is_sg_cys[None, :] & (dist < DISULFIDE_MAX_DIST)

    is_metal = np.isin(atoms.element, list(METAL_ELEMENTS))
    metal_coordination = (is_metal[:, None] | is_metal[None, :]) & (dist < METAL_COORDINATION_MAX_DIST)

    bonded_like = bonded_like | disulfide | metal_coordination
    upper_triangle = np.triu(np.ones((n, n), dtype=bool), k=1)

    real_clashes = clash_mask & upper_triangle & ~bonded_like

    if severity_weighted:
        severity = np.clip(vdw_sum - overlap_tolerance - dist, 0, None)
        score = severity[real_clashes].sum()
    else:
        score = real_clashes.sum()

    return float(score) / n * 1000
