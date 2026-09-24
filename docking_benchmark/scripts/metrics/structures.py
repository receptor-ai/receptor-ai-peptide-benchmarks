from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io as strucio
from rdkit import Chem

# crystallization buffer/cryoprotectant, not protein biology - already stripped from the
# shipped peptide.pdb/complex.pdb, matters only when re-deriving one from a raw extraction
JUNK_HETATM = {
    "SO4", "PO4", "CL", "NO3", "SO3", "CO3", "ACT", "FMT", "CIT", "TLA", "MLA", "MLI", "SIN",
    "AZI", "SCN", "IOD", "GOL", "EDO", "DMS", "IPA", "MPD", "PEG", "PGE", "PG4", "P6G", "1PE",
    "15P", "12P", "PG5", "BU3", "BU2", "DIO", "HEZ", "HEX", "MRD", "MRY", "DTD", "DTU", "DTT",
    "BME", "TRS", "MES", "EPE", "BTB", "CAC", "B3P", "DMF", "PDO", "TBU", "TRT", "CCN", "UNX",
    "LI", "AU", "OLA", "OLC", "CD", "K", "NI", "CO",
}
WATER = {"HOH", "WAT", "H2O"}


def load_native_peptide(pdb_path: Path, keep_chain: str | None = None) -> struc.AtomArray:
    """Heavy-atom peptide coordinates, ready to compare against a docked pose.
    entries/<case_id>/peptide.pdb already holds a single chain with water/buffer
    removed, so keep_chain is normally unnecessary; pass it only when loading
    directly from a raw, unfiltered multi-copy extraction.
    """
    atoms = strucio.load_structure(str(pdb_path), model=1)  # NMR files carry many models
    if keep_chain is not None:
        atoms = atoms[atoms.chain_id == keep_chain]
    atoms = atoms[~np.isin(atoms.res_name, list(JUNK_HETATM | WATER))]
    return atoms[atoms.element != "H"]


def _flatten_for_matching(mol: Chem.Mol) -> Chem.Mol:
    # GetSubstructMatches requires exact bond-order/charge equality by default, so a pose
    # with its own valid bond perception could fail to match. Only element and connectivity
    # matter for finding which atom is which.
    flat = Chem.RWMol(mol)
    for atom in flat.GetAtoms():
        atom.SetFormalCharge(0)
    for bond in flat.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
    return flat.GetMol()


def load_pred_peptide_coord(sdf_path: Path, native_mol: Chem.Mol) -> np.ndarray:
    """Coordinates of a predicted pose, reordered to match native_mol's atoms.
    Docking tools don't preserve atom order, and symmetric residues (Val/Leu
    methyls, Tyr/Phe rings, ...) can match more than one way - RDKit's first
    match may pair the wrong look-alikes, so we try them all and keep the
    lowest-RMSD one.
    """
    mol = Chem.RemoveAllHs(next(Chem.SDMolSupplier(str(sdf_path), removeHs=False)))
    # GetSubstructMatches only requires native_mol's atoms to appear somewhere in mol - it
    # does not care whether mol has extra atoms beyond the match, so a pose with a surplus
    # atom would otherwise match and get scored silently, with the extra atom just ignored
    if mol.GetNumAtoms() != native_mol.GetNumAtoms():
        raise ValueError(
            f"predicted pose in {sdf_path} has {mol.GetNumAtoms()} heavy atoms, "
            f"expected {native_mol.GetNumAtoms()} to match the native peptide"
        )
    matches = _flatten_for_matching(mol).GetSubstructMatches(
        _flatten_for_matching(native_mol), uniquify=False, maxMatches=10000
    )

    native_coord = native_mol.GetConformer().GetPositions()
    pred_coord = mol.GetConformer().GetPositions()

    best_coord, best_rmsd = None, float("inf")
    for match in matches:
        if len(match) != native_mol.GetNumAtoms():
            continue
        candidate = pred_coord[list(match)]
        rmsd = float(np.sqrt(np.mean(np.sum((candidate - native_coord) ** 2, axis=1))))
        if rmsd < best_rmsd:
            best_coord, best_rmsd = candidate, rmsd

    if best_coord is None:
        raise ValueError(
            f"could not match predicted pose atoms to the native peptide "
            f"(0/{native_mol.GetNumAtoms()} atoms matched) - "
            f"check the predicted SDF has the same connectivity as peptide.sdf"
        )
    return best_coord


def combine(receptor: struc.AtomArray, native_peptide: struc.AtomArray, peptide_coord: np.ndarray) -> struc.AtomArray:
    # re-use the native peptide's element/res_id/chain_id annotations - atom i in
    # peptide_coord is the same chemical atom as atom i in native_peptide, just moved
    peptide = native_peptide.copy()
    peptide.coord = peptide_coord
    return receptor + peptide
