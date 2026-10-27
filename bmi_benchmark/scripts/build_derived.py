#!/usr/bin/env python3
"""Build the derived peptide files of BMI-200 and BMI-MODES from the deposited ones.

The deposited peptide files (peptide.cif, pose_NN.cif) are never edited. Next to each one
this script writes three derived files (pose_NN_fixed* in BMI-MODES):

    peptide_fixed.cif              every deposited atom unchanged, plus the heavy atoms listed
                                   in meta.json; the origin of each atom is in _bmi_atom_origin
    peptide_fixed.sdf              the same heavy atoms in the same order, every formal charge
                                   0, hydrogens implicit
    peptide_fixed_protonated.sdf   the same heavy atoms protonated as at pH 7.4, with explicit
                                   hydrogens appended after the heavy atoms

What to add is read from meta.json ("derived" of a BMI-200 entry, "derived__pose_NN" of a
BMI-MODES site): added_residues, added_atoms with their coordinates, and added_links. This
script places no atom itself.

Bonds come from the file's _chem_comp_bond (each component's bonds from the RCSB Chemical
Component Dictionary, Kekule bond orders), its covalent and disulfide _struct_conn links,
and the peptide bond C(i)-N(i+1) between consecutive residues of a chain when C-N < 2.0 A.
In peptide_fixed.sdf every atom is neutral and carries the hydrogens its valence implies.

peptide_fixed_protonated.sdf ionises each group whose pKa lies on the ionised side of 7.4
(PROTONATION below): carboxylic, sulfinic and sulfonic acids are anions, phosphate OH
groups are O- (pSer, pThr and pTyr are dianions), aliphatic amines are NH3+/NH2+/NH+ (Lys,
N-termini, N-methyl and N-alkyl amines, dimethyl-Lys) and guanidines are guanidinium (Arg
and its analogues). Everything else stays neutral: His (pKa 6.0-6.5, H on NE2, the CCD
Kekule form), Cys, Tyr, amides, pyridines, anilines, azoles. The side-chain amine of Dap
(DNP) stays neutral too: its pKa is about 6.3 in a peptide. Hydrogens are placed by RDKit
from the heavy atoms.

Needs gemmi and RDKit. A component the deposited file does not define (an added cap) is
read from --ccd, and downloaded from RCSB when it is not there.

Usage:
    python scripts/build_derived.py                   # every entry and pose
    python scripts/build_derived.py 1VDN A_2AOF#02    # some entries, poses or whole sites
    python scripts/build_derived.py --check           # compare with the files on disk
"""
import argparse
import csv
import io
import json
import re
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import gemmi
from rdkit import Chem, RDLogger
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
CCD_URL = "https://files.rcsb.org/ligands/download/{}.cif"
BOND_ORDER = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE}
LINKS = (gemmi.ConnectionType.Covale, gemmi.ConnectionType.Disulf)
PEPTIDE_BOND_MAX_A = 2.0
# origin of an atom; "hydrogen" is used in the protonated SDF only
ORIGINS = ("deposited", "rebuilt", "oxt", "declared_cap", "stand_in_cap")
ORIGIN_CAT = "_bmi_atom_origin."
PROTONATED = "_protonated.sdf"

# Protonation at pH 7.4: (group, SMARTS, index of the atom that takes the charge, charge).
# The SMARTS are matched on the neutral molecule with aromaticity perceived. pKa values are
# those of model compounds and peptides.
AMINE = "[NX3;+0;!a;!$(N-[!#6;!#1]);!$(N-[#6;!X4])]"  # N bonded to sp3 carbons only
PROTONATION = (
    ("carboxylic acid", "[CX3](=O)[OX2H1]", 2, -1),            # C-terminus 3.3, Asp 3.9, Glu 4.3
    ("phosphate OH", "[PX4](=O)[OX2H1]", 2, -1),               # pKa 1-2 and 5.6-6.1: both ionised
    ("sulfate or sulfonate OH", "[SX4](=O)(=O)[OX2H1]", 3, -1),  # pKa < 1
    ("sulfinic acid OH", "[SX3](=O)[OX2H1]", 2, -1),           # pKa ~2
    ("aliphatic amine", AMINE, 0, 1),                          # Lys 10.5, N-terminus 7.7
    ("guanidine or amidine",                                   # Arg 12.5, acylguanidine 8.2
     "[NX2;+0;!a;!$(N-[!#6;!#1]);!$(N-C#N)]=[CX3;!a]([#7X3;!a])[#6,#7]", 0, 1),
)
NOT_IONISED = (
    # Dap: an acylated alpha-N lowers the beta-amine pKa to ~6.3 (Lan et al., ChemBioChem 2010)
    ("Dap side-chain amine", "[NX3;H2;+0]-[CH2]-[CX4]-[NX3]-[CX3]=O", 0),
)


class BuildError(Exception):
    pass


@dataclass
class Item:
    name: str      # 1VDN, or A_2AOF#02 for a pose
    dir: Path
    stem: str      # peptide, or pose_02
    meta_key: str  # derived, or derived__pose_02
    title: str     # entries/1VDN/peptide, or sites/A_2AOF/pose_02

    @property
    def deposited(self):
        return self.dir / f"{self.stem}.cif"

    def derived(self, suffix):
        return self.dir / f"{self.stem}_fixed{suffix}"


def items(root=ROOT):
    with open(root / "bmi200" / "MANIFEST.csv", newline="") as f:
        for row in csv.DictReader(f):
            pid = row["pdb_id"]
            yield Item(pid, root / "bmi200" / "entries" / pid, "peptide", "derived",
                       f"entries/{pid}/peptide")
    with open(root / "bmi_modes" / "poses.csv", newline="") as f:
        for row in csv.DictReader(f):
            site, nn = row["site_id"], f"{int(row['pose']):02d}"
            yield Item(f"{site}#{nn}", root / "bmi_modes" / "sites" / site, f"pose_{nn}",
                       f"derived__pose_{nn}", f"sites/{site}/pose_{nn}")


# --------------------------------------------------------------------------- CIF helpers
def value(tok):
    """CIF token -> string; '?' and '.' stay as they are"""
    return tok if tok in ("?", ".") else gemmi.cif.as_string(tok)


def token(s):
    """string -> CIF token; '?' and '.' stay as they are"""
    s = str(s)
    return s if s in ("?", ".") else gemmi.cif.quote(s)


def category(block, cat):
    """(tags, rows) of an mmCIF category, the rows as raw CIF tokens"""
    t = block.find_mmcif_category(cat)
    tags = [tag[len(cat):] for tag in t.tags]
    return tags, [[row[i] for i in range(len(tags))] for row in t]


def put_category(block, cat, tags, rows):
    """replace a category (in place) or add it at the end"""
    loop = block.init_mmcif_loop(cat, tags)
    for r in rows:
        loop.add_row(r)


def add_rows(block, cat, new):
    """append rows (dicts tag -> string) to a category, creating it when absent"""
    tags, rows = category(block, cat)
    tags = tags or list(new[0])
    rows += [[token(d.get(t, "?")) for t in tags] for d in new]
    put_category(block, cat, tags, rows)


def column(block, cat, tag):
    tags, rows = category(block, cat)
    return [value(r[tags.index(tag)]) for r in rows] if tag in tags else []


def res_ref(chain, comp, seq, icode=""):
    return f"{chain}/{comp}/{seq}{icode}"


def parse_ref(ref):
    """'B/ACE/0' -> ('B', 'ACE', 0, ''); 'B/ARG/8/OXT' -> ('B', 'ARG', 8, '', 'OXT')"""
    parts = ref.split("/")
    m = re.fullmatch(r"(-?\d+)([A-Za-z]?)", parts[2]) if len(parts) in (3, 4) else None
    if not m:
        raise BuildError(f"bad reference {ref!r}")
    return (parts[0], parts[1], int(m[1]), m[2]) + tuple(parts[3:])


def atom_records(tags, rows):
    """parsed _atom_site rows"""
    out = []
    for raw in rows:
        d = dict(zip(tags, (value(x) for x in raw)))
        icode = d.get("pdbx_PDB_ins_code", "?")
        out.append(dict(id=d["id"], chain=d["auth_asym_id"], comp=d["label_comp_id"],
                        seq=int(d["auth_seq_id"]), icode="" if icode in ("?", ".") else icode,
                        name=d["label_atom_id"], element=d["type_symbol"],
                        xyz=tuple(float(d[k]) for k in ("Cartn_x", "Cartn_y", "Cartn_z"))))
    return out


def read_atoms(path):
    return atom_records(*category(gemmi.cif.read(str(path)).sole_block(), "_atom_site."))


def rkey(a):
    return a["chain"], a["comp"], a["seq"], a["icode"]


def akey(a):
    return a["chain"], a["comp"], a["seq"], a["icode"], a["name"]


def ccd_block(code, ccd_dir):
    path = Path(ccd_dir) / f"{code}.cif"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(CCD_URL.format(code), timeout=60) as r:
            path.write_bytes(r.read())
    return gemmi.cif.read(str(path)).sole_block()


# --------------------------------------------------------------------------- build
def build(item, ccd_dir):
    """-> ({path: text} of the three derived files, summary)"""
    meta = json.loads((item.dir / "meta.json").read_text())
    der = meta.get(item.meta_key)
    if der is None:
        raise BuildError(f"meta.json has no {item.meta_key!r} record")
    doc = gemmi.cif.read(str(item.deposited))
    block = doc.sole_block()
    tags, rows = category(block, "_atom_site.")
    col = {t: i for i, t in enumerate(tags)}
    dep = atom_records(tags, rows)

    # residues in file order; the deposited atoms of a residue are contiguous
    order, members = [], {}
    for i, a in enumerate(dep):
        k = rkey(a)
        if k not in members:
            members[k] = []
            order.append(k)
        elif order[-1] != k:
            raise BuildError(f"atoms of {res_ref(*k)} are not contiguous")
        members[k].append(("dep", i))
    new_res = {}
    for r in der["added_residues"]:
        k = parse_ref(r["residue"])
        where, anchor = r["anchor"].split()
        ak = parse_ref(anchor)
        if k in members or ak not in members or where not in ("before", "after"):
            raise BuildError(f"cannot add {r['residue']} {r['anchor']}")
        order.insert(order.index(ak) + (where == "after"), k)
        members[k] = []
        new_res[k] = r
    for a in der["added_atoms"]:
        *k, name = parse_ref(a["atom"])
        k = tuple(k)
        if k not in members:
            raise BuildError(f"{a['atom']}: residue not in the file")
        if a["origin"] not in ORIGINS[1:]:
            raise BuildError(f"{a['atom']}: unknown origin {a['origin']!r}")
        if any((dep[x]["name"] if src == "dep" else x["name"]) == name for src, x in members[k]):
            raise BuildError(f"{a['atom']}: already in the file")
        members[k].append(("new", dict(a, name=name)))

    out_rows, origin, atoms = [], [], []
    for k in order:
        first = next((rows[x] for src, x in members[k] if src == "dep"), None)
        for src, x in members[k]:
            n = str(len(out_rows) + 1)
            if src == "dep":
                raw = list(rows[x])
                origin.append(("deposited", raw[col["id"]]))
                atoms.append(dep[x])
            else:
                if first is not None:   # an atom added to a deposited residue
                    base = dict(zip(tags, first))
                else:                   # an atom of an added residue
                    r = new_res[k]
                    base = dict(group_PDB=r["group"], label_comp_id=token(k[1]),
                                label_asym_id=token(r["label_asym"]),
                                label_entity_id=token(r["entity"]),
                                label_seq_id=token(r["label_seq"]),
                                pdbx_PDB_ins_code=token(k[3] or "?"), auth_seq_id=str(k[2]),
                                auth_asym_id=token(k[0]), pdbx_PDB_model_num="1")
                xyz = dict(zip(("Cartn_x", "Cartn_y", "Cartn_z"), (f"{v:.3f}" for v in x["xyz"])))
                base.update(type_symbol=token(x["element"]), label_atom_id=token(x["name"]),
                            label_alt_id=".", occupancy="1", B_iso_or_equiv=token(x["b"]),
                            pdbx_formal_charge="?", **xyz)
                raw = [base.get(t, "?") for t in tags]
                origin.append((x["origin"], "."))
                atoms.append(dict(chain=k[0], comp=k[1], seq=k[2], icode=k[3], name=x["name"],
                                  element=x["element"],
                                  xyz=tuple(round(float(v), 3) for v in x["xyz"])))
            raw[col["id"]] = n
            out_rows.append(raw)
    put_category(block, "_atom_site.", tags, out_rows)

    # the categories that describe what the atom records refer to
    ents = set(column(block, "_entity.", "id"))
    asyms = set(column(block, "_struct_asym.", "id"))
    comps = set(column(block, "_chem_comp.", "id"))
    for k, r in new_res.items():
        if r["entity"] not in ents:
            add_rows(block, "_entity.", [dict(id=r["entity"], type=r["entity_type"])])
            ents.add(r["entity"])
        if r["label_asym"] not in asyms:
            add_rows(block, "_struct_asym.", [dict(id=r["label_asym"], entity_id=r["entity"])])
            asyms.add(r["label_asym"])
        if k[1] not in comps:
            add_rows(block, "_chem_comp.", [dict(id=k[1], type=".")])
            comps.add(k[1])
    symbols = set(column(block, "_atom_type.", "symbol"))
    for el in sorted({a["element"] for a in atoms} - symbols):
        add_rows(block, "_atom_type.", [dict(symbol=el)])
    defined = set(column(block, "_chem_comp_bond.", "comp_id"))
    ordinals = [int(x) for x in column(block, "_chem_comp_bond.", "pdbx_ordinal") if x.isdigit()]
    for comp in sorted({a["comp"] for a in atoms} - defined):
        tags_b = ["atom_id_1", "atom_id_2", "value_order", "?pdbx_aromatic_flag"]
        t = ccd_block(comp, ccd_dir).find("_chem_comp_bond.", tags_b)
        new = []
        for r in t:
            ordinals.append(max(ordinals, default=0) + 1)
            new.append(dict(comp_id=comp, atom_id_1=r.str(0), atom_id_2=r.str(1),
                            value_order=r.str(2), pdbx_aromatic_flag=r.str(3) if r.has(3) else "N",
                            pdbx_ordinal=str(ordinals[-1])))
        if new:
            add_rows(block, "_chem_comp_bond.", new)

    # added links
    index = {akey(a): i for i, a in enumerate(atoms)}
    sc_tags, _ = category(block, "_struct_conn.")
    if not sc_tags:
        sc_tags = ["id", "conn_type_id"] + [f"ptnr{p}_{c}" for p in (1, 2) for c in (
            "label_asym_id", "label_comp_id", "label_seq_id", "label_atom_id", "auth_asym_id",
            "auth_seq_id")] + ["pdbx_dist_value"]
    ids = column(block, "_struct_conn.", "id")
    new = []
    for link in der["added_links"]:
        ends = [index.get(parse_ref(link[f"atom{p}"])) for p in (1, 2)]
        if None in ends:
            raise BuildError(f"link {link}: atom not in the file")
        num = [int(m[1]) for x in ids for m in [re.fullmatch(re.escape(link["type"]) + r"(\d+)", x)]
               if m]
        ids.append(f"{link['type']}{max(num, default=0) + 1}")
        d = dict(id=ids[-1], conn_type_id=link["type"],
                 pdbx_dist_value=f"{dist(atoms[ends[0]], atoms[ends[1]]):.3f}")
        for p, i in zip((1, 2), ends):
            r = dict(zip(tags, (value(v) for v in out_rows[i])))
            d.update({f"ptnr{p}_label_asym_id": r["label_asym_id"],
                      f"ptnr{p}_label_comp_id": r["label_comp_id"],
                      f"ptnr{p}_label_seq_id": r["label_seq_id"],
                      f"ptnr{p}_label_atom_id": r["label_atom_id"],
                      f"ptnr{p}_auth_asym_id": r["auth_asym_id"],
                      f"ptnr{p}_auth_comp_id": r["label_comp_id"],
                      f"ptnr{p}_auth_seq_id": r["auth_seq_id"],
                      f"ptnr{p}_auth_atom_id": r["label_atom_id"],
                      f"pdbx_ptnr{p}_PDB_ins_code": r["pdbx_PDB_ins_code"],
                      f"ptnr{p}_symmetry": "1_555"})
        new.append({t: d.get(t, "?") for t in sc_tags})
    if new:
        add_rows(block, "_struct_conn.", new)
        types = column(block, "_struct_conn_type.", "id")
        if types:
            for t in sorted({link["type"] for link in der["added_links"]} - set(types)):
                add_rows(block, "_struct_conn_type.", [dict(id=t)])
    put_category(block, ORIGIN_CAT, ["atom_site_id", "origin", "deposited_atom_site_id"],
                 [[str(i + 1), o, d] for i, (o, d) in enumerate(origin)])
    cif_text = doc.as_string()

    # the molecule
    bonds, unresolved = bond_table(block, atoms)
    if unresolved:
        raise BuildError(f"_struct_conn rows that name no atom of the file: {unresolved}")
    mol = make_mol(atoms, bonds)
    labels = [res_ref(*rkey(a)) + "/" + a["name"] for a in atoms]
    origins = [o for o, _ in origin]
    molh = Chem.AddHs(protonate(mol), addCoords=True)
    Chem.AssignStereochemistryFrom3D(molh)
    labels_h = labels + [labels[a.GetNeighbors()[0].GetIdx()].rsplit("/", 1)[0] + "/H"
                         for a in molh.GetAtoms() if a.GetIdx() >= len(atoms)]
    origins_h = origins + ["hydrogen"] * (molh.GetNumAtoms() - len(atoms))
    charge = Chem.GetFormalCharge(molh)
    files = {item.derived(".cif"): cif_text,
             item.derived(".sdf"): sdf_text(mol, item.title + "_fixed", labels, origins),
             item.derived(PROTONATED): sdf_text(molh, item.title + "_fixed_protonated",
                                                labels_h, origins_h, net_charge=charge)}
    summary = dict(deposited=len(dep), added=len(atoms) - len(dep), heavy=len(atoms),
                   hydrogens=molh.GetNumAtoms() - len(atoms), bonds=mol.GetNumBonds(),
                   charge=charge)
    return files, summary


def dist(a, b):
    return sum((p - q) ** 2 for p, q in zip(a["xyz"], b["xyz"])) ** 0.5


def bond_table(block, atoms):
    """{(i, j): bond type} from _chem_comp_bond, _struct_conn and the peptide bonds"""
    index = {akey(a): i for i, a in enumerate(atoms)}
    by_res, order = {}, []
    for i, a in enumerate(atoms):
        if rkey(a) not in by_res:
            by_res[rkey(a)] = {}
            order.append(rkey(a))
        by_res[rkey(a)][a["name"]] = i
    comp_bonds = {}
    for r in block.find("_chem_comp_bond.", ["comp_id", "atom_id_1", "atom_id_2", "value_order"]):
        comp_bonds.setdefault(r.str(0), []).append((r.str(1), r.str(2), r.str(3).upper()))
    bonds = {}
    for k in order:
        names = by_res[k]
        if len(names) > 1 and k[1] not in comp_bonds:
            raise BuildError(f"no _chem_comp_bond rows for {k[1]}")
        for a1, a2, o in comp_bonds.get(k[1], []):
            if a1 in names and a2 in names:
                if o not in BOND_ORDER:
                    raise BuildError(f"{k[1]} {a1}-{a2}: bond order {o}")
                bonds[tuple(sorted((names[a1], names[a2])))] = BOND_ORDER[o]
    # _struct_conn with auth columns has no insertion code: then the atom must be unique
    # without one
    loose = {}
    for (chain, comp, seq, _, name), i in index.items():
        loose.setdefault((chain, comp, seq, name), []).append(i)
    unresolved = []
    for c in gemmi.make_structure_from_block(block).connections:
        if c.type not in LINKS:
            continue
        ends = []
        for p in (c.partner1, c.partner2):
            k = (p.chain_name, p.res_id.name, p.res_id.seqid.num, p.atom_name)
            i = index.get(k[:3] + (p.res_id.seqid.icode.strip(),) + k[3:])
            if i is None and len(loose.get(k, [])) == 1:
                i = loose[k][0]
            ends.append(i)
        if None in ends:
            unresolved.append(c.name)
        else:
            bonds.setdefault(tuple(sorted(ends)), Chem.BondType.SINGLE)
    for prev, cur in zip(order, order[1:]):
        i, j = by_res[prev].get("C"), by_res[cur].get("N")
        if prev[0] == cur[0] and i is not None and j is not None \
                and dist(atoms[i], atoms[j]) < PEPTIDE_BOND_MAX_A:
            bonds.setdefault((i, j), Chem.BondType.SINGLE)
    return bonds, unresolved


def make_mol(atoms, bonds):
    """heavy-atom molecule, every formal charge 0, hydrogens implicit, stereo from 3D"""
    rw = Chem.RWMol()
    for a in atoms:
        rw.AddAtom(Chem.Atom(gemmi.Element(a["element"]).name))
    for (i, j), bt in sorted(bonds.items()):
        rw.AddBond(i, j, bt)
    mol = rw.GetMol()
    conf = Chem.Conformer(len(atoms))
    conf.Set3D(True)
    for i, a in enumerate(atoms):
        conf.SetAtomPosition(i, Point3D(*a["xyz"]))
    mol.AddConformer(conf, assignId=True)
    try:
        # keep the CCD Kekule bond orders as they are
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL
                         ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    except Exception as e:  # noqa: BLE001
        raise BuildError(f"RDKit cannot sanitize the molecule: {e}") from e
    Chem.AssignStereochemistryFrom3D(mol)
    return mol


def protonate(mol):
    """the molecule as at pH 7.4 (PROTONATION, NOT_IONISED), hydrogens implicit"""
    q = Chem.Mol(mol)
    Chem.SetAromaticity(q)
    keep = {m[i] for _, smarts, i in NOT_IONISED
            for m in q.GetSubstructMatches(Chem.MolFromSmarts(smarts))}
    rw = Chem.RWMol(mol)
    for _, smarts, i, charge in PROTONATION:
        for m in q.GetSubstructMatches(Chem.MolFromSmarts(smarts)):
            if m[i] not in keep:
                rw.GetAtomWithIdx(m[i]).SetFormalCharge(charge)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out, Chem.SanitizeFlags.SANITIZE_ALL
                         ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    except Exception as e:  # noqa: BLE001
        raise BuildError(f"RDKit cannot sanitize the protonated molecule: {e}") from e
    Chem.AssignStereochemistryFrom3D(out)
    return out


def sdf_text(mol, title, labels, origins, net_charge=None):
    m = Chem.Mol(mol)
    m.SetProp("_Name", title)
    for a, lab, org in zip(m.GetAtoms(), labels, origins):
        a.SetProp("label", lab)
        a.SetProp("origin", org)
    Chem.CreateAtomStringPropertyList(m, "label")
    Chem.CreateAtomStringPropertyList(m, "origin")
    if net_charge is not None:
        m.SetProp("net_charge", str(net_charge))
    buf = io.StringIO()
    w = Chem.SDWriter(buf)
    w.write(m)
    w.close()
    return buf.getvalue()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("items", nargs="*", help="BMI-200 PDB ids, BMI-MODES poses (SITE#NN) or sites")
    ap.add_argument("--root", type=Path, default=ROOT, help="the bmi_benchmark directory")
    ap.add_argument("--ccd", type=Path, default=Path.home() / ".cache" / "bmi_benchmark" / "ccd",
                    help="directory of CCD component files (default: %(default)s)")
    ap.add_argument("--check", action="store_true",
                    help="compare with the derived files on disk instead of writing them")
    args = ap.parse_args(argv)
    want = set(args.items)
    todo = [it for it in items(args.root)
            if not want or it.name in want or it.name.split("#")[0] in want]
    if want and not todo:
        ap.error("nothing matches " + " ".join(sorted(want)))
    n_bad = 0
    for it in todo:
        try:
            files, summary = build(it, args.ccd)
        except BuildError as e:
            print(f"{it.name}: ERROR {e}")
            n_bad += 1
            continue
        if args.check:
            stale = [p.name for p, text in files.items() if not p.exists() or p.read_text() != text]
            n_bad += bool(stale)
            print(f"{it.name}: {'DIFFERS ' + ' '.join(stale) if stale else 'ok'}")
            continue
        for p, text in files.items():
            p.write_text(text)
        print(f"{it.name}: " + " ".join(f"{k}={v}" for k, v in summary.items()))
    print(f"{len(todo)} built, {n_bad} {'differ or ' if args.check else ''}failed")
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
