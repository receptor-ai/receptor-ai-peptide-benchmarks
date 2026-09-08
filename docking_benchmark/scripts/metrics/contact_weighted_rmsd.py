import numpy as np
import biotite.structure as struc


def contact_weighted_rmsd(
    peptide_native_coord: np.ndarray,
    peptide_pred_coord: np.ndarray,
    receptor_native_coord: np.ndarray,
    receptor_pred_coord: np.ndarray,
    lam: float = 2.5,
    power: float = 2.0,
) -> float:
    """
    RMSD between peptide_pred and peptide_native (after superimposing the
    receptor), where each peptide atom's error is weighted by how close that
    atom is to the receptor in the native structure -- atoms deep in the
    binding contact matter more than solvent-exposed ones.

    All coordinate arrays must have matching atom order (peptide_native_coord[i]
    <-> peptide_pred_coord[i], receptor_native_coord[j] <-> receptor_pred_coord[j]).
    """
    _, transformation = struc.superimpose(receptor_native_coord, receptor_pred_coord)
    peptide_pred_fit = transformation.apply(peptide_pred_coord)

    dist_to_receptor = np.linalg.norm(
        peptide_native_coord[:, None, :] - receptor_native_coord[None, :, :], axis=-1
    ).min(axis=1)
    weights = np.exp(-(dist_to_receptor / lam) ** power)

    sq_dev = np.sum((peptide_pred_fit - peptide_native_coord) ** 2, axis=1)
    return np.sqrt(np.sum(weights * sq_dev) / np.sum(weights))
