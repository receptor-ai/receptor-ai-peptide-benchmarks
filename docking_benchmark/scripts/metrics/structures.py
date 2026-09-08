from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io as strucio
from rdkit import Chem

# crystallization buffer / cryoprotectant / solvent / phasing-agent components: not part of
# protein biology, dropped when a native peptide.pdb is re-derived from a raw extraction.
# The peptide.pdb/complex.pdb files shipped in docking985/entries/ already have these removed;
# this constant matters only if you re-run the extraction against your own raw structures.
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


def load_pred_peptide_coord(sdf_path: Path) -> np.ndarray:
    mol = next(Chem.SDMolSupplier(str(sdf_path), removeHs=False))
    return mol.GetConformer().GetPositions()


def combine(receptor: struc.AtomArray, native_peptide: struc.AtomArray, peptide_coord: np.ndarray) -> struc.AtomArray:
    # re-use the native peptide's element/res_id/chain_id annotations -- atom i in
    # peptide_coord is the same chemical atom as atom i in native_peptide, just moved
    peptide = native_peptide.copy()
    peptide.coord = peptide_coord
    return receptor + peptide
