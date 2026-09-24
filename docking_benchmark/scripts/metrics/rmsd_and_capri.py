import tempfile
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdb

BACKBONE_ATOMS = {"N", "CA", "C", "O"}


def _superposed_pred(receptor_native_coord, receptor_pred_coord, peptide_pred_coord):
    _, transformation = struc.superimpose(receptor_native_coord, receptor_pred_coord)
    return transformation.apply(peptide_pred_coord)


def _rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def _masked_rmsd(peptide_native, peptide_pred_coord, receptor_native_coord, receptor_pred_coord, mask):
    # no per-residue symmetry correction needed here: load_pred_peptide_coord()
    # (metrics/structures.py) already picks the atom correspondence, among every
    # valid substructure match, that minimizes overall RMSD - that subsumes
    # fixing up specific symmetric residue types (Tyr/Phe rings, etc.) as a
    # special case, so peptide_pred_coord is already the best assignment
    pred_fit = _superposed_pred(receptor_native_coord, receptor_pred_coord, peptide_pred_coord)
    return _rmsd(peptide_native.coord[mask], pred_fit[mask])


def backbone_rmsd(peptide_native, peptide_pred_coord, receptor_native_coord, receptor_pred_coord) -> float:
    mask = np.isin(peptide_native.atom_name, list(BACKBONE_ATOMS))
    return _masked_rmsd(peptide_native, peptide_pred_coord, receptor_native_coord, receptor_pred_coord, mask)


def heavy_atom_rmsd(peptide_native, peptide_pred_coord, receptor_native_coord, receptor_pred_coord) -> float:
    mask = np.ones(len(peptide_native), dtype=bool)
    return _masked_rmsd(peptide_native, peptide_pred_coord, receptor_native_coord, receptor_pred_coord, mask)


def ca_rmsd(peptide_native, peptide_pred_coord, receptor_native_coord, receptor_pred_coord) -> float:
    mask = peptide_native.atom_name == "CA"
    return _masked_rmsd(peptide_native, peptide_pred_coord, receptor_native_coord, receptor_pred_coord, mask)


def sidechain_rmsd(peptide_native, peptide_pred_coord, receptor_native_coord, receptor_pred_coord) -> float:
    mask = ~np.isin(peptide_native.atom_name, list(BACKBONE_ATOMS))
    return _masked_rmsd(peptide_native, peptide_pred_coord, receptor_native_coord, receptor_pred_coord, mask)


def capri_class(dockq_score: float) -> str:
    if dockq_score < 0.23:
        return "Incorrect"
    if dockq_score < 0.49:
        return "Acceptable"
    if dockq_score < 0.80:
        return "Medium"
    return "High"


def _free_chain_id(receptor: struc.AtomArray) -> str:
    used = set(receptor.chain_id)
    for letter in "LPQZYXWVUTSRO":
        if letter not in used:
            return letter
    raise ValueError("no free chain id for the peptide")


def _write_complex_pdb(receptor: struc.AtomArray, peptide: struc.AtomArray, peptide_coord: np.ndarray,
                        peptide_chain_id: str, path: Path) -> None:
    # receptor chain ids are kept as-is: relabeling multi-chain receptors onto a single
    # letter can collide residue numbers between the original chains
    pep = peptide.copy()
    pep.coord = peptide_coord
    pep.chain_id[:] = peptide_chain_id

    pdb_file = pdb.PDBFile()
    pdb_file.set_structure(receptor + pep)
    pdb_file.write(str(path))


def dockq_metrics(receptor: struc.AtomArray, peptide_native: struc.AtomArray, peptide_pred_coord: np.ndarray) -> dict:
    """Fnat, fnonnat, iRMSD, LRMSD, DockQ score and CAPRI class, via DockQ.
    Imports the DockQ package lazily - only this function needs it, not the
    RMSD variants or clashscore elsewhere in this module, so importing
    metrics.rmsd_and_capri (or calling backbone_rmsd/heavy_atom_rmsd/etc.)
    works without DockQ installed. Raises the normal ImportError, with a
    clearer message, only when this function is actually called.
    """
    try:
        import DockQ.DockQ as dockq
    except ImportError as error:
        raise ImportError(
            "dockq_metrics() needs the DockQ package (pip install DockQ) - "
            "not required for the RMSD variants or clashscore."
        ) from error

    peptide_chain_id = _free_chain_id(receptor)
    # DockQ's own loader drops HETATM-only chains (e.g. a lone glycan given its own chain
    # letter) - exclude those here too, or it KeyErrors looking them up internally
    real_chains = {c for c in set(receptor.chain_id) if np.any(~receptor.hetero[receptor.chain_id == c])}
    chain_map = {c: c for c in real_chains}
    chain_map[peptide_chain_id] = peptide_chain_id

    with tempfile.TemporaryDirectory() as tmp:
        native_path = Path(tmp) / "native.pdb"
        model_path = Path(tmp) / "model.pdb"
        _write_complex_pdb(receptor, peptide_native, peptide_native.coord, peptide_chain_id, native_path)
        _write_complex_pdb(receptor, peptide_native, peptide_pred_coord, peptide_chain_id, model_path)

        # small_molecule=True keeps HETATM residues - needed when the peptide is
        # entirely non-canonical (e.g. CSPSet) or it would vanish from DockQ's own
        # parsing entirely. Only turned on when actually needed: DockQ's small-molecule
        # path is dramatically slower (symmetry-corrected RMSD) for every case otherwise.
        peptide_is_all_hetatm = bool(np.all(peptide_native.hetero))
        native_structure = dockq.load_PDB(str(native_path), small_molecule=peptide_is_all_hetatm)
        model_structure = dockq.load_PDB(str(model_path), small_molecule=peptide_is_all_hetatm)
        result, _ = dockq.run_on_all_native_interfaces(
            model_structure, native_structure, chain_map=chain_map
        )

    # result may contain other native interfaces too (e.g. between two receptor
    # chains); keep only the one(s) actually involving the peptide chain
    peptide_interfaces = [v for k, v in result.items() if peptide_chain_id in k]
    if not peptide_interfaces:
        raise ValueError(f"no native interface found involving the peptide chain ({peptide_chain_id})")
    interface = max(peptide_interfaces, key=lambda v: v["DockQ"])

    # DockQ's "small molecule" scoring path (triggered when the peptide is entirely
    # non-canonical/HETATM residues, e.g. CSPSet) returns DockQ+LRMSD only - no
    # fnat/fnonnat/iRMSD, since it doesn't define residue-residue contacts the same way
    return {
        "fnat": interface.get("fnat"),
        "fnonnat": interface.get("fnonnat"),
        "iRMSD": interface.get("iRMSD"),
        "LRMSD": interface["LRMSD"],
        "DockQ": interface["DockQ"],
        "CAPRI_class": capri_class(interface["DockQ"]),
    }
