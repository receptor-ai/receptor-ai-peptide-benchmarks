#!/usr/bin/env python3
"""Check the coordinate files of BMI-200 and BMI-MODES (run before every release).

Offline checks, for every entry and pose:
  build      the derived files are exactly what build_derived.py makes from the deposited
             file and meta.json
  layout     *_fixed.cif holds every deposited atom unchanged and every added atom of
             meta.json with its origin; the SDFs hold the same heavy atoms in the same order
             and the same bonds; *_fixed.sdf has every formal charge 0; the protonated SDF
             adds hydrogens only, one more or one fewer than the neutral molecule implies
             for each unit of charge
  protonation  in the protonated SDF no carboxylic, phosphoric, sulfinic or sulfonic acid
             OH is left, every aliphatic amine (except a Dap side chain) and every guanidine
             or amidine is protonated, and every charge sits on one of those groups
  complete   every residue has the heavy atoms of its CCD definition (leaving atoms
             excepted) unless meta.json lists it under not_rebuilt
  chirality  every CCD stereocentre has the handedness of the CCD ideal coordinates
  added      bonds from added atoms; added atoms against the rest of the peptide and the
             receptor; hydrogens of the protonated SDF against non-bonded atoms

With --pdb, against the PDB entries (downloaded to --pdb-dir when missing):
  pdb        every deposited peptide and receptor atom is in the entry (a pose may be
             moved as a rigid body)
  ends       each chain end is completed as the deposition says: OXT only at a true
             C-terminus, a declared cap only where the sequence declares one, a stand-in only
             where residues are unobserved or the chain continues through the crystal
  lattice    peptide atoms within bonding distance of a crystal neighbour (the receptor,
             other chains, symmetry mates) that the files do not bond (warnings); of an atom
             with alternate locations only the one the files use counts, and for a pose taken
             from an altloc branch only that branch

Needs gemmi, RDKit and numpy.

Usage:
    python scripts/verify.py                   # offline checks, everything
    python scripts/verify.py --pdb 7ABT A_4X34 # all checks, some entries / sites
"""
import argparse
import csv
import json
import sys
import urllib.request
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_derived as bd  # noqa: E402

RDLogger.DisableLog("rdApp.*")
PDB_URL = "https://files.rcsb.org/download/{}.cif.gz"
CONTACT_MIN_A = 2.2   # an added heavy atom closer than this to a non-bonded heavy atom fails
H_CONTACT_A = 1.6     # a hydrogen closer than this to a non-bonded atom is reported
BOND_A = 2.0          # a crystal contact shorter than this is a bond
CAP_N = {"ACE", "FOR", "BOC", "CBZ", "PHQ", "BZ", "MYR", "PLM", "SIN", "ACY"}
CAP_C = {"NH2", "NME", "OME", "NH3", "ETA"}
# pH 7.4 protonation, written for the protonated molecule (aromaticity perceived)
ACID_LEFT = "[$([OX2H1][CX3]=O),$([OX2H1][PX4]=O),$([OX2H1][SX4](=O)=O),$([OX2H1][SX3]=O)]"
BASE_LEFT = ("[$([NX3;+0;!a;!$(N-[!#6;!#1]);!$(N-[#6;!X4])]),"
             "$([NX2;+0;!a]=[CX3;!a]([#7X3])[#6,#7])]")
DAP_AMINE = "[NX3;H2;+0]-[CH2]-[CX4]-[NX3]-[CX3]=O"
ANION = "[$([OX1-][CX3]=O),$([OX1-][PX4]=O),$([OX1-][SX4]=O),$([OX1-][SX3]=O)]"
CATION = ("[$([NX4+;!$(N~[!#6;!#1]);!$(N-[#6;!X4])]),"
          "$([NX3+;!a]=[CX3;!a]([#7X3])[#6,#7])]")


class Result:
    def __init__(self, name):
        self.name, self.errors, self.warnings = name, [], []

    def error(self, check, msg):
        self.errors.append(f"{check}: {msg}")

    def warn(self, check, msg):
        self.warnings.append(f"{check}: {msg}")


def ccd_def(code, ccd_dir, cache={}):
    """atoms {name: (element, leaving, stereo, reference xyz or None)}, neighbours {name: set}"""
    if code not in cache:
        block = bd.ccd_block(code, ccd_dir)
        atoms = {}
        cols = ["atom_id", "type_symbol", "?pdbx_leaving_atom_flag", "?pdbx_stereo_config",
                "?pdbx_model_Cartn_x_ideal", "?pdbx_model_Cartn_y_ideal",
                "?pdbx_model_Cartn_z_ideal", "?model_Cartn_x", "?model_Cartn_y", "?model_Cartn_z"]
        for r in block.find("_chem_comp_atom.", cols):
            get = [r.str(i) if r.has(i) else "?" for i in range(len(cols))]
            xyz = None
            for k in (4, 7):
                try:
                    xyz = np.array([float(get[k + j]) for j in range(3)])
                    break
                except ValueError:
                    continue
            atoms[get[0]] = (get[1].upper(), get[2].upper() == "Y", get[3].upper(), xyz)
        nb = {}
        for r in block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2"]):
            nb.setdefault(r.str(0), set()).add(r.str(1))
            nb.setdefault(r.str(1), set()).add(r.str(0))
        ctype = (block.find_value("_chem_comp.type") or "").strip("\"'").upper()
        cache[code] = atoms, nb, ctype
    return cache[code]


def read_sdf(path):
    mol = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)[0]
    return mol


# --------------------------------------------------------------------------- offline
def offline(item, ccd_dir, res):
    try:
        files, _ = bd.build(item, ccd_dir)
    except bd.BuildError as e:
        res.error("build", str(e))
        return None
    for p, text in files.items():
        if not p.exists():
            res.error("build", f"{p.name} is missing")
        elif p.read_text() != text:
            res.error("build", f"{p.name} differs from what build_derived.py makes")
    der = json.loads((item.dir / "meta.json").read_text())[item.meta_key]
    dblock = gemmi.cif.read(str(item.deposited)).sole_block()
    fblock = gemmi.cif.read(str(item.derived(".cif"))).sole_block()
    dtags, drows = bd.category(dblock, "_atom_site.")
    ftags, frows = bd.category(fblock, "_atom_site.")
    otags, orows = bd.category(fblock, bd.ORIGIN_CAT)
    if ftags != dtags or len(orows) != len(frows):
        res.error("layout", "_atom_site columns or the origin table do not match")
        return None
    fatoms = bd.atom_records(ftags, frows)
    idcol = dtags.index("id")
    dep_rows = {bd.value(r[idcol]): r for r in drows}
    seen, origin = set(), []
    for row, a, o in zip(frows, fatoms, orows):
        org, dep_id = bd.value(o[1]), bd.value(o[2])
        origin.append(org)
        if bd.value(o[0]) != a["id"]:
            res.error("layout", f"origin row {bd.value(o[0])} is not atom {a['id']}")
        if org == "deposited":
            d = dep_rows.get(dep_id)
            if d is None or dep_id in seen:
                res.error("layout", f"atom {a['id']}: deposited id {dep_id} unknown or used twice")
                continue
            seen.add(dep_id)
            if [x for i, x in enumerate(row) if i != idcol] != [x for i, x in enumerate(d) if i != idcol]:
                res.error("layout", f"deposited atom {dep_id} is changed")
        elif org not in bd.ORIGINS:
            res.error("layout", f"atom {a['id']}: unknown origin {org}")
    if seen != set(dep_rows):
        res.error("layout", f"{len(set(dep_rows) - seen)} deposited atoms are missing")
    labels = [bd.res_ref(*bd.rkey(a)) + "/" + a["name"] for a in fatoms]
    added = {lab: (a, o) for lab, a, o in zip(labels, fatoms, origin) if o != "deposited"}
    meta_added = {x["atom"]: x for x in der["added_atoms"]}
    if set(added) != set(meta_added):
        res.error("layout", f"added atoms differ from meta.json: {sorted(set(added) ^ set(meta_added))}")
    for lab, x in meta_added.items():
        if lab in added and (added[lab][1] != x["origin"]
                             or max(abs(p - q) for p, q in zip(added[lab][0]["xyz"], x["xyz"])) > 6e-4):
            res.error("layout", f"{lab}: origin or position differs from meta.json")

    # SDFs
    mol = read_sdf(item.derived(".sdf"))
    molp = read_sdf(item.derived(bd.PROTONATED))
    pname = item.derived(bd.PROTONATED).name
    if mol is None or molp is None:
        res.error("layout", "an SDF does not load in RDKit")
        return None
    xyz = mol.GetConformer().GetPositions()
    n = len(fatoms)
    if mol.GetNumAtoms() != n or any(
            at.GetSymbol().upper() != a["element"].upper() for at, a in zip(mol.GetAtoms(), fatoms)) \
            or np.abs(xyz - np.array([a["xyz"] for a in fatoms])).max() > 1e-3:
        res.error("layout", f"{item.derived('.sdf').name}: atoms differ from the CIF")
        return None
    if mol.GetProp("atom.prop.origin").split() != origin or mol.GetProp("atom.prop.label").split() != labels:
        res.error("layout", f"{item.derived('.sdf').name}: origin or label list differs from the CIF")
    q = [at.GetIdx() for at in mol.GetAtoms() if at.GetFormalCharge()]
    if q:
        res.error("layout", f"{item.derived('.sdf').name}: {len(q)} charged atoms")
    xp = molp.GetConformer().GetPositions()
    if molp.GetNumAtoms() < n or np.abs(xp[:n] - xyz).max() > 1e-4 or any(
            molp.GetAtomWithIdx(i).GetAtomicNum() != mol.GetAtomWithIdx(i).GetAtomicNum() for i in range(n)):
        res.error("layout", f"{pname} does not start with the heavy atoms of the fixed SDF")
        return None
    for at in list(molp.GetAtoms())[n:]:
        nbrs = at.GetNeighbors()
        if at.GetAtomicNum() != 1 or len(nbrs) != 1 or nbrs[0].GetIdx() >= n:
            res.error("layout", f"{pname} has more than hydrogens after the heavy atoms")
            break
    heavy_bonds = [b for b in molp.GetBonds() if b.GetBeginAtomIdx() < n and b.GetEndAtomIdx() < n]
    if len(heavy_bonds) != mol.GetNumBonds() or any(
            (b2 := mol.GetBondBetweenAtoms(b.GetBeginAtomIdx(), b.GetEndAtomIdx())) is None
            or b2.GetBondType() != b.GetBondType() for b in heavy_bonds):
        res.error("layout", f"{pname}: the heavy-atom bonds differ from the fixed SDF")
    for i in range(n):
        # protonation moves one H per unit of charge
        h = sum(1 for x in molp.GetAtomWithIdx(i).GetNeighbors() if x.GetAtomicNum() == 1)
        want = mol.GetAtomWithIdx(i).GetTotalNumHs() + molp.GetAtomWithIdx(i).GetFormalCharge()
        if h != want:
            res.error("layout", f"{labels[i]}: {h} explicit H in {pname}, {want} expected")
    if int(molp.GetProp("net_charge")) != Chem.GetFormalCharge(molp):
        res.error("layout", f"{pname}: net_charge is not the sum of the formal charges")

    # protonation at pH 7.4
    heavy = Chem.RemoveHs(molp)
    Chem.SetAromaticity(heavy)

    def hits(smarts):
        return {m[0] for m in heavy.GetSubstructMatches(Chem.MolFromSmarts(smarts))}

    charged = {a.GetIdx() for a in heavy.GetAtoms() if a.GetFormalCharge()}
    for what, idx in (("acid left un-ionised", hits(ACID_LEFT)),
                      ("base left unprotonated", hits(BASE_LEFT) - hits(DAP_AMINE)),
                      ("charge on no ionisable group", charged - hits(ANION) - hits(CATION))):
        if idx:
            res.error("protonation", f"{what}: {', '.join(labels[i] for i in sorted(idx))}")

    # completeness and chirality against the CCD
    not_rebuilt = {r["residue"]: set(r["missing_atoms"]) for r in der.get("not_rebuilt", [])}
    residues = {}
    for a in fatoms:
        residues.setdefault(bd.rkey(a), {})[a["name"]] = np.array(a["xyz"])
    index = {lab: i for i, lab in enumerate(labels)}

    def linked(ref, name):
        """the atom has a bond to another residue"""
        i = index.get(f"{ref}/{name}")
        return i is not None and any(bd.rkey(fatoms[x.GetIdx()]) != bd.rkey(fatoms[i])
                                     for x in mol.GetAtomWithIdx(i).GetNeighbors())

    for k, present in residues.items():
        ref = bd.res_ref(*k)
        atoms, nb, _ = ccd_def(k[1], ccd_dir)
        expected = {nm for nm, (el, leaving, _, _) in atoms.items() if el not in ("H", "D") and not leaving}
        missing = set()
        for x in expected - set(present):
            # an end atom whose only heavy neighbour is bonded to another residue was replaced
            # by that bond (a carboxyl O at an isopeptide, OXT of an acyl cap, Cys SG at a staple)
            heavy = [y for y in nb.get(x, ()) if atoms[y][0] not in ("H", "D")]
            if not (len(heavy) == 1 and linked(ref, heavy[0])):
                missing.add(x)
        unknown = set(present) - set(atoms)
        if unknown:
            res.error("complete", f"{ref}: atoms not in the CCD definition {sorted(unknown)}")
        if missing:
            if missing <= not_rebuilt.get(ref, set()):
                res.warn("complete", f"{ref}: {len(missing)} atoms missing (listed as not rebuilt)")
            else:
                res.error("complete", f"{ref}: missing {sorted(missing)}")
        for c, (el, _, stereo, cxyz) in atoms.items():
            if stereo not in ("R", "S") or c not in present or cxyz is None:
                continue
            hv = sorted(x for x in nb.get(c, ()) if atoms[x][0] not in ("H", "D")
                        and x in present and atoms[x][3] is not None)[:3]
            if len(hv) < 3:
                continue
            vol = [np.linalg.det(np.array([X[x] - X[c] for x in hv]))
                   for X in (present, {x: atoms[x][3] for x in hv + [c]})]
            if np.sign(vol[0]) != np.sign(vol[1]):
                res.warn("chirality", f"{ref}/{c}: handedness differs from the CCD ({stereo})")

    # added atoms: bonds, contacts
    topo = Chem.GetDistanceMatrix(mol)
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if origin[i] == "deposited" and origin[j] == "deposited":
            continue
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        junction = bd.rkey(fatoms[i]) != bd.rkey(fatoms[j])
        lo, hi = (1.28, 1.45) if junction else (1.1, 2.0)
        if not lo <= d <= hi:
            res.error("added", f"bond {labels[i]}-{labels[j]} {d:.2f} A")
    for i, o in enumerate(origin):
        if o == "oxt":
            ca, c, ox = (labels.index(labels[i].rsplit("/", 1)[0] + "/" + x) for x in ("CA", "C", "O"))
            d = float(np.linalg.norm(xyz[i] - xyz[c]))
            nrm = np.cross(xyz[ca] - xyz[c], xyz[ox] - xyz[c])
            off = abs(float((xyz[i] - xyz[c]) @ nrm / np.linalg.norm(nrm)))
            if not 1.20 <= d <= 1.30 or off > 0.15:
                res.error("added", f"{labels[i]}: C-OXT {d:.2f} A, {off:.2f} A off the carboxyl plane")
    rec_path = item.dir / ("receptor.cif" if item.stem == "peptide" else f"receptor_{item.stem[5:]}.cif")
    rec = np.array([a["xyz"] for a in bd.read_atoms(rec_path) if a["element"] not in ("H", "D")])
    new = [i for i, o in enumerate(origin) if o != "deposited"]
    worst = (99.0, "")
    for i in new:
        far = (topo[i] > 3)
        dp = np.linalg.norm(xyz[far] - xyz[i], axis=1)
        dr = np.linalg.norm(rec - xyz[i], axis=1)
        for d, what in ((dp, "peptide"), (dr, "receptor")):
            if len(d) and d.min() < worst[0]:
                worst = (float(d.min()), f"{labels[i]} to the {what}")
            if len(d) and d.min() < CONTACT_MIN_A:
                res.error("added", f"{labels[i]} is {d.min():.2f} A from a {what} atom")
    # hydrogens of the protonated SDF
    tp = Chem.GetDistanceMatrix(molp)
    hs = [at.GetIdx() for at in molp.GetAtoms() if at.GetAtomicNum() == 1]
    close = []
    for i in hs:
        far = tp[i] > 2
        d = np.linalg.norm(xp[far] - xp[i], axis=1)
        if len(d) and d.min() < H_CONTACT_A:
            close.append(float(d.min()))
    if close:
        res.warn("added", f"{len(close)} hydrogens within {H_CONTACT_A} A of a non-bonded atom "
                          f"(closest {min(close):.2f} A)")
    return dict(mol=mol, labels=labels, origin=origin, fatoms=fatoms, der=der,
                min_added_contact=worst)


# --------------------------------------------------------------------------- against the PDB
def fetch(url, dest):
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=120) as r:
            dest.write_bytes(r.read())


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return R, Qc - R @ Pc


def pdb_atoms(st):
    out = []
    for ch in st[0]:
        for r in ch:
            for at in r:
                out.append(dict(chain=ch.name, subchain=r.subchain, label_seq=r.label_seq, comp=r.name,
                                seq=r.seqid.num, icode=r.seqid.icode.strip(), name=at.name,
                                altloc=at.altloc if at.altloc != "\0" else "", el=at.element.name,
                                xyz=np.array(at.pos.tolist())))
    return out


def frame(dep, dep_label_seq, patoms, spec):
    """R, t taking the deposited file into the PDB frame (identity when it is already there)"""
    X = np.array([a["xyz"] for a in patoms])
    hit = 0
    for a in dep[:50]:
        d = np.linalg.norm(X - np.array(a["xyz"]), axis=1)
        hit += d.min() < 0.002
    if hit >= 0.9 * min(50, len(dep)):
        return np.eye(3), np.zeros(3), 0.0
    sub, alt = ((spec or "/") + "/").split("/")[:2]
    alt = "" if alt in ("-", "") else alt
    known = any(a["subchain"] == sub for a in patoms)
    n_poly = sum(s is not None for s in dep_label_seq)
    best = None
    for sc in ([sub] if known else sorted({a["subchain"] for a in patoms})):
        idx = {}
        for i, a in enumerate(patoms):
            if a["subchain"] == sc and a["label_seq"] is not None and (not alt or a["altloc"] in ("", alt)):
                if (a["label_seq"], a["name"]) not in idx or a["altloc"] == alt:
                    idx[(a["label_seq"], a["name"])] = i
        pairs = [(np.array(a["xyz"]), X[idx[(s, a["name"])]]) for a, s in zip(dep, dep_label_seq)
                 if s is not None and (s, a["name"]) in idx]
        # the subchain named by the pose label needs only a fit; another one must cover the pose
        if len(pairs) < (4 if known else max(4, 0.8 * n_poly)):
            continue
        P, Q = np.array([p for p, _ in pairs]), np.array([q for _, q in pairs])
        R, t = kabsch(P, Q)
        rms = float(np.sqrt(((P @ R.T + t - Q) ** 2).sum(1).mean()))
        if best is None or rms < best[2]:
            best = (R, t, rms)
    return best


def is_amino_acid(chem_comp_type):
    t = chem_comp_type.upper()
    return "PEPTIDE" in t and ("LINKING" in t or "TERMINUS" in t)


def receptor_links(block, where):
    """covalent and disulfide bonds that the PDB entry records between the peptide and the
    receptor; where maps (auth chain, seq, icode, comp, atom) -> (file role, atom label)"""
    cols = [f"{p}{c}" for p in ("1", "2") for c in (
        "_auth_asym_id", "_auth_comp_id", "_auth_seq_id", "_label_atom_id", "_symmetry")]
    cols = ["conn_type_id"] + [f"ptnr{c}" for c in cols] + [
        "?pdbx_ptnr1_PDB_ins_code", "?pdbx_ptnr2_PDB_ins_code", "?pdbx_dist_value"]
    out = []
    for r in block.find("_struct_conn.", cols):
        if not r.str(0).lower().startswith(("covale", "disulf")):
            continue
        ends = []
        for p in (0, 1):
            icode = r.str(11 + p) if r.has(11 + p) else "?"
            k = (r.str(1 + 5 * p), int(r.str(3 + 5 * p)), "" if icode in ("?", ".") else icode,
                 r.str(2 + 5 * p), r.str(4 + 5 * p))
            ends.append((where.get(k), r.str(5 + 5 * p)))
        roles = [e[0][0] if e[0] else None for e in ends]
        if sorted(roles, key=str) != ["peptide", "receptor"]:
            continue
        pep, rec = ends if roles[0] == "peptide" else ends[::-1]
        out.append(dict(type=r.str(0).lower(), peptide_atom=pep[0][1], receptor_atom=rec[0][1],
                        symmetry=[pep[1], rec[1]],
                        distance_A=float(r.str(13)) if r.has(13) and r.str(13) not in ("?", ".") else None))
    return out


def against_pdb(item, pdb_id, spec, pdb_dir, ccd_dir, off, res):
    path = Path(pdb_dir) / f"{pdb_id}.cif.gz"
    fetch(PDB_URL.format(pdb_id), path)
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    block = gemmi.cif.read(str(path)).sole_block()
    patoms = pdb_atoms(st)
    dblock = gemmi.cif.read(str(item.deposited)).sole_block()
    tags, rows = bd.category(dblock, "_atom_site.")
    dep = bd.atom_records(tags, rows)
    lseq = [None if bd.value(r[tags.index("label_seq_id")]) in (".", "?")
            else int(bd.value(r[tags.index("label_seq_id")])) for r in rows]
    fr = frame(dep, lseq, patoms, spec)
    if fr is None:
        res.error("pdb", f"cannot place {item.deposited.name} in {pdb_id}")
        return
    R, t, rms = fr
    if rms > 0.01:
        res.error("pdb", f"rigid fit to {pdb_id} leaves {rms:.3f} A")
    model = st[0]
    ns = gemmi.NeighborSearch(model, gemmi.UnitCell(), 3).populate(include_h=False)
    rec_path = item.dir / ("receptor.cif" if item.stem == "peptide" else f"receptor_{item.stem[5:]}.cif")
    pep_match, where, used_alt = [], {}, {}
    for role, path_, atoms in (("peptide", item.deposited, dep), ("receptor", rec_path, bd.read_atoms(rec_path))):
        lost = 0
        for a in atoms:
            p = R @ np.array(a["xyz"]) + t
            m = [x.to_cra(model) for x in ns.find_atoms(gemmi.Position(*p), "\0", radius=0.01)]
            m = [c for c in m if c.atom.name == a["name"] and c.residue.name == a["comp"]]
            lost += not m
            if role == "peptide":
                pep_match.append(m[0] if m else None)
            if m:
                c = m[0]
                where[(c.chain.name, c.residue.seqid.num, c.residue.seqid.icode.strip(), c.residue.name,
                       c.atom.name)] = (role, bd.res_ref(*bd.rkey(a)) + "/" + a["name"])
                if c.atom.altloc != "\0":
                    used_alt[(c.chain.name, c.residue.seqid.num, c.residue.seqid.icode)] = c.atom.altloc
        if lost:
            res.error("pdb", f"{lost} atoms of {path_.name} are not in {pdb_id}")
    # the peptide-receptor bonds of the entry are recorded in meta.json
    links = receptor_links(block, where)
    key = "receptor_links" + item.meta_key[len("derived"):]
    recorded = json.loads((item.dir / "meta.json").read_text()).get(key)
    if recorded is None:
        res.error("links", f"meta.json has no {key!r} record")
    elif sorted(map(json.dumps, recorded)) != sorted(map(json.dumps, links)):
        res.error("links", f"{key} in meta.json differs from the bonds {pdb_id} records")
    if off is None:
        return links
    # chain ends
    cell = st.cell if st.cell.is_crystal() else gemmi.UnitCell()
    st.setup_cell_images()
    nsx = gemmi.NeighborSearch(model, cell, 5).populate(include_h=False)
    self_res = {(c.chain.name, c.residue.seqid.num, c.residue.seqid.icode) for c in pep_match if c}
    chem_type = dict(zip(bd.column(block, "_chem_comp.", "id"), bd.column(block, "_chem_comp.", "type")))
    scheme = {}
    for r in block.find("_pdbx_poly_seq_scheme.", ["asym_id", "seq_id", "mon_id"]):
        scheme.setdefault(r.str(0), {}).setdefault(int(r.str(1)), r.str(2))
    mol, labels, origin, fatoms = off["mol"], off["labels"], off["origin"], off["fatoms"]
    tc = off["der"]["terminal_completion"]
    pose_alt = ((spec or "/") + "/").split("/")[1].strip("-")

    def contacts(i_dep):
        """bond-length contacts of a deposited atom with atoms outside the file"""
        p = gemmi.Position(*(R @ np.array(dep[i_dep]["xyz"]) + t))
        out = []
        for mk in nsx.find_atoms(p, "\0", radius=BOND_A):
            c = mk.to_cra(model)
            if c.residue.name in ("HOH", "DOD"):
                continue
            # an alternate location counts only when the files use it (or, for a residue
            # outside the files, when it is the branch the pose comes from)
            alt = c.atom.altloc
            if alt != "\0" and alt != used_alt.get(
                    (c.chain.name, c.residue.seqid.num, c.residue.seqid.icode), pose_alt or alt):
                continue
            img = cell.find_nearest_image(p, c.atom.pos, gemmi.Asu.Any) if cell.is_crystal() else None
            d = img.dist() if img else p.dist(c.atom.pos)
            sym = bool(img) and not (img.sym_idx == 0 and img.pbc_shift == (0, 0, 0))
            if d < BOND_A and not ((c.chain.name, c.residue.seqid.num, c.residue.seqid.icode) in self_res and not sym):
                out.append(dict(d=round(d, 2), atom=c.atom.name, sym=sym, entity=c.residue.entity_id,
                                what=f"{c.chain.name}:{c.residue.name}{c.residue.seqid.num}.{c.atom.name}",
                                copy=(c.residue.subchain,) + ((img.sym_idx, tuple(img.pbc_shift)) if sym else ())))
        return out

    fidx = {bd.akey(a): i for i, a in enumerate(fatoms)}
    dep_ent = {c.residue.entity_id for c in pep_match if c}
    ext_all = [contacts(i) for i in range(len(dep))]
    # a copy (symmetry mate or other chain) that four or more peptide atoms touch overlaps the
    # peptide: a special position or an alternative orientation, not a bond
    touching = {}
    for i, es in enumerate(ext_all):
        for e in es:
            touching.setdefault(e["copy"], set()).add(i)
    overlap = {k for k, v in touching.items() if len(v) >= 4}
    for k in sorted(overlap, key=str):
        res.warn("lattice", f"{len(touching[k])} peptide atoms within {BOND_A} A of "
                            f"{'a symmetry mate' if len(k) > 1 else 'another chain'} of subchain {k[0]} "
                            f"(overlapping copy)")
    for i, es in enumerate(ext_all):
        for e in es:
            if e["copy"] not in overlap:
                res.warn("lattice", f"{labels[fidx[bd.akey(dep[i])]]} {e['d']} A from {e['what']}"
                                    f"{' (symmetry mate)' if e['sym'] else ''}")
    by_chain = {}
    for i, (a, s) in enumerate(zip(dep, lseq)):
        if s is not None and pep_match[i] is not None:
            by_chain.setdefault(pep_match[i].residue.subchain, []).append((s, i))
    for sub, members in by_chain.items():
        seqs = sorted({s for s, _ in members})
        for side, s_end, atom in (("N", seqs[0], "N"), ("C", seqs[-1], "C")):
            i_end = next((i for s, i in members if s == s_end and dep[i]["name"] == atom), None)
            if i_end is None:
                continue
            a_end = dep[i_end]
            fi = fidx[bd.akey(a_end)]
            nb = [x.GetIdx() for x in mol.GetAtomWithIdx(fi).GetNeighbors()
                  if bd.rkey(fatoms[x.GetIdx()]) != bd.rkey(a_end)]
            nb_dep = [x for x in nb if origin[x] == "deposited"]
            nb_new = [x for x in nb if origin[x] != "deposited"]
            seq = scheme.get(pep_match[i_end].residue.subchain, {})
            flank = [seq[k] for k in sorted(seq) if (k < s_end if side == "N" else k > s_end)]
            ext = [e for e in ext_all[i_end] if e["copy"] not in overlap]
            if nb_dep:
                cls = "bonded"
            elif any(not e["sym"] and e["entity"] not in dep_ent for e in ext):
                cls = "bonded to the receptor"
            elif any(e["entity"] in dep_ent and e["d"] < 1.6 and e["atom"] == ("C" if side == "N" else "N")
                     for e in ext):
                cls = "continues through the crystal"
            elif not flank:
                cls = "terminus"
            elif all(c in (CAP_N if side == "N" else CAP_C) or not is_amino_acid(chem_type.get(c, ""))
                     for c in flank):
                cls = "declared cap, unobserved"
            else:
                cls = "truncated"
            got = sorted({origin[x] for x in nb_new})
            if side == "C":
                oxt = fidx.get(bd.akey(a_end)[:4] + ("OXT",))
                if oxt is not None and origin[oxt] == "oxt":
                    got.append("oxt")
            note = f"{side}:not added" in tc
            ref = bd.res_ref(*bd.rkey(a_end))
            if cls in ("bonded", "bonded to the receptor"):
                ok = not got
            elif cls == "terminus":
                ccd_oxt = "OXT" in ccd_def(a_end["comp"], ccd_dir)[0]
                dep_oxt = any(bd.akey(a) == bd.akey(a_end)[:4] + ("OXT",) for a in dep)
                ok = got == (["oxt"] if side == "C" and ccd_oxt and not dep_oxt else [])
            elif cls == "declared cap, unobserved":
                ok = got == ["declared_cap"] or (not got and note)
            else:
                ok = got == ["stand_in_cap"] or (not got and note)
            msg = f"{ref} {side}-end: {cls}; added {got or 'nothing'}"
            if ok:
                if cls in ("continues through the crystal", "bonded to the receptor"):
                    res.warn("ends", msg + " " + "; ".join(f"{e['what']} {e['d']} A" for e in ext))
            else:
                res.error("ends", msg)
    return links


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("items", nargs="*", help="BMI-200 PDB ids, BMI-MODES poses (SITE#NN) or sites")
    ap.add_argument("--root", type=Path, default=bd.ROOT, help="the bmi_benchmark directory")
    ap.add_argument("--ccd", type=Path, default=Path.home() / ".cache" / "bmi_benchmark" / "ccd",
                    help="directory of CCD component files (default: %(default)s)")
    ap.add_argument("--pdb", action="store_true", help="also check against the PDB entries")
    ap.add_argument("--pdb-dir", type=Path, default=Path.home() / ".cache" / "bmi_benchmark" / "pdb",
                    help="directory of PDB entries (default: %(default)s)")
    ap.add_argument("--json", type=Path, help="write every result here")
    ap.add_argument("-q", "--quiet", action="store_true", help="print errors only")
    args = ap.parse_args(argv)
    want = set(args.items)
    todo = [it for it in bd.items(args.root) if not want or it.name in want or it.name.split("#")[0] in want]
    if want and not todo:
        ap.error("nothing matches " + " ".join(sorted(want)))
    with open(args.root / "bmi_modes" / "poses.csv", newline="") as f:
        poses = {f"{r['site_id']}#{int(r['pose']):02d}": r for r in csv.DictReader(f)}
    results, n_err, n_warn = [], 0, 0
    for it in todo:
        res = Result(it.name)
        try:
            off = offline(it, args.ccd, res)
            if args.pdb:
                row = poses.get(it.name)
                against_pdb(it, row["pdb_id"] if row else it.name, row["label"] if row else None,
                            args.pdb_dir, args.ccd, off, res)
        except Exception as e:  # noqa: BLE001
            res.error("crash", f"{type(e).__name__}: {e}")
        n_err += bool(res.errors)
        n_warn += bool(res.warnings)
        results.append(dict(item=it.name, errors=res.errors, warnings=res.warnings,
                            min_added_contact=(off or {}).get("min_added_contact")))
        if res.errors or (res.warnings and not args.quiet):
            print(f"{it.name}: {len(res.errors)} errors, {len(res.warnings)} warnings")
            for m in res.errors:
                print("   ERROR", m)
            if not args.quiet:
                for m in res.warnings:
                    print("   warning", m)
    if args.json:
        args.json.write_text(json.dumps(results, indent=1))
    print(f"{len(todo)} checked: {n_err} with errors, {n_warn} with warnings")
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main())
