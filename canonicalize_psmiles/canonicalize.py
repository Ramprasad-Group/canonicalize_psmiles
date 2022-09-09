from __future__ import annotations

import logging
import pprint
import re
from typing import Dict, Tuple

from rdkit import Chem
from rdkit.Chem.rdmolops import GetSymmSSSR


def canonicalize(psmiles) -> str:
    """Canonicalize PSMILES string

    (1) Unify ([\*]CCO[\*] -> [\*]COO[\*]); apply multiple times

    (2) Reduce ([\*]CCOCCO[\*] -> [\*]CCO[\*]); apply multiple times

    Returns:
        str: canonicalized SMILES string
    """

    # Apply unify once
    new_ps = unify(psmiles)

    new_ps = add_brackets(new_ps)

    # recursively apply reduce
    for num in range(20):
        new = reduce_multiplication(new_ps)
        new = add_brackets(new)
        if new == new_ps:
            break
        new_ps = new

    return add_brackets(new_ps)


def get_connection_info(mol=None, symbol="*") -> Dict:
    """Get connection information a PSMILES string.

    Args:
        mol (Chem.RWMol, optional): _description_. Defaults to None.
        symbol (str, optional): _description_. Defaults to "*".

    Returns:
        Dict: Dictionary containing information of the mol
    """   

    ret_dict = {}

    stars_indices, stars_type, all_symbols, all_index = [], [], [], []
    for star_idx, atom in enumerate(mol.GetAtoms()):
        all_symbols.append(atom.GetSymbol())
        all_index.append(atom.GetIdx())
        if symbol in atom.GetSymbol():
            stars_indices.append(star_idx)
            stars_type.append(atom.GetSmarts())

    stars_bond = mol.GetBondBetweenAtoms(stars_indices[0], stars_indices[1])
    if stars_bond:
        stars_bond = stars_bond.GetBondType()

    ret_dict["symbols"] = all_symbols
    ret_dict["index"] = all_index

    ret_dict["star"] = {
        "index": stars_indices,
        "atom_type": stars_type,
        "bond_type": stars_bond,
    }

    # multiple neighbors are possible
    neighbor_indices = [
        [x.GetIdx() for x in mol.GetAtomWithIdx(stars_indices[0]).GetNeighbors()],
        [x.GetIdx() for x in mol.GetAtomWithIdx(stars_indices[1]).GetNeighbors()],
    ]

    neighbors_type = [
        [mol.GetAtomWithIdx(x).GetSmarts() for x in neighbor_indices[0]],
        [mol.GetAtomWithIdx(x).GetSmarts() for x in neighbor_indices[1]],
    ]

    # Bonds between stars and neighbors
    neighbor_bonds = [
        [
            mol.GetBondBetweenAtoms(stars_indices[0], x).GetBondType()
            for x in neighbor_indices[0]
        ],
        [
            mol.GetBondBetweenAtoms(stars_indices[1], x).GetBondType()
            for x in neighbor_indices[1]
        ],
    ]
    s_path = None
    if neighbor_indices[0][0] != neighbor_indices[1][0]:
        s_path = Chem.GetShortestPath(
            mol, neighbor_indices[0][0], neighbor_indices[1][0]
        )

    ret_dict["neighbor"] = {
        "index": neighbor_indices,
        "atom_type": neighbors_type,
        "bond_type": neighbor_bonds,
        "path": s_path,
    }

    # Stereo info
    stereo_info = []
    for b in mol.GetBonds():
        bond_type = b.GetStereo()
        if bond_type != Chem.rdchem.BondStereo.STEREONONE:
            idx = [b.GetBeginAtomIdx(), b.GetEndAtomIdx()]
            neigh_idx = b.GetStereoAtoms()
            stereo_info.append(
                {
                    "bond_type": bond_type,
                    "atom_idx": idx,
                    "bond_idx": b.GetIdx(),
                    "neighbor_idx": list(neigh_idx),
                }
            )

    ret_dict["stereo"] = stereo_info

    # Ring info
    ring_info = mol.GetRingInfo()
    ret_dict["atom_rings"] = ring_info.AtomRings()
    ret_dict["bond_rings"] = ring_info.BondRings()

    return ret_dict


def reduce_multiplication(psmiles: str) -> str:
    """Reduces PSMILES manifolds.
    Apply multiple times to reduce to smallest PSMILES.

    E.g.,

    [\*]CCC[\*] --> [\*]C[\*]

    [\*]COCO[\*] --> [\*]CO[\*]

    [\*]COCCOCCOCCOC[\*] -> [\*]COCCOC[\*]

    Returns:
        str: reduced PSMILES
    """
    sm = psmiles

    # skip if rings is polymer
    if re.findall("\d", sm):
        return psmiles

    # Handle the case of PE
    pe_test = list(set(sm.replace("[*]", "").upper()))
    if len(pe_test) == 1 and pe_test[0] == "C":
        return "[*]C[*]"

    # All other cases
    d = {}
    for sublen in range(1, int(len(sm) / 2)):
        for i in range(0, len(sm) - sublen):
            sub = sm[i : i + sublen]
            cnt = sm.count(sub)
            if cnt >= 2 and sub not in d:
                if "[" not in sub and "]" not in sub and "*" not in sub:
                    d[sub] = cnt
    if not d:
        return psmiles

    longest_string = sorted(d, key=lambda k: len(k), reverse=True)[0]
    check_s = sm.replace("[*]", "").replace(longest_string, "")

    if not check_s:
        sm = sm.replace(longest_string, "", d[longest_string] - 1)

    return sm


def unify(psmiles: str) -> str:
    """Unify PSMILES strings

    E.g.,

    [\*]COC[\*] --> [\*]COC[\*]

    [\*]OCC[\*] --> [\*]COC[\*]

    [\*]COC[\*] --> [\*]COC[\*]

    Returns:
        str: unified polymer SMILES string
    """

    # Make atom indices visible if in DEBUG mode
    if logging.DEBUG >= logging.root.level:
        from rdkit.Chem.Draw import IPythonConsole

        IPythonConsole.drawOptions.addAtomIndices = True
        # IPythonConsole.drawOptions.addBondIndices = True

    # (1) Get connection info
    mol = get_mol(psmiles)

    info = get_connection_info(mol)
    logging.debug(f"(1) Set labels and get connection dict \n {pprint.pformat(info)}")

    # SPECIAL CASES:
    # (1.1) if *C* return *C*
    if psmiles == "[*]C[*]":
        logging.debug("Found [*]C[*]; returning [*]C[*]")
        return Chem.MolToSmiles(mol)

    # (1.2) if neighbors are already connected PS or '[*]/C=C(\[*])C(C)(C)C'
    if info["neighbor"]["path"] and len(info["neighbor"]["path"]) == 2:
        logging.debug(f"Neighbors are already connected. Return unified SMILES here.")
        return Chem.MolToSmiles(mol)

    # (1.3) if *C(*)(F)F: neighbors are the same
    if info["neighbor"]["index"][0][0] == info["neighbor"]["index"][1][0]:
        logging.debug("Neighbors are the same. Returning unified SMILES here")
        return Chem.MolToSmiles(mol)

    mol.GetAtomWithIdx(info["neighbor"]["index"][0][0]).SetProp("atomLabel", "N1")
    mol.GetAtomWithIdx(info["neighbor"]["index"][1][0]).SetProp("atomLabel", "N2")
    if logging.DEBUG >= logging.root.level:
        nb_display(mol)

    # (2) Connect neighbors with bond type of neighbors
    logging.debug(
        f"(2) Add bond between neighbors {info['neighbor']['index']}; bond type {info['neighbor']['bond_type'][0][0]} {info['neighbor']['bond_type'][1][0]} "
    )

    mol.AddBond(
        info["neighbor"]["index"][0][0],
        info["neighbor"]["index"][1][0],
        info["neighbor"]["bond_type"][0][0],
    )
    if logging.DEBUG >= logging.root.level:
        nb_display(mol)

    # (3) Remove stars and bonds
    logging.debug(f"(3) Remove stars and bonds")
    mol.RemoveBond(info["star"]["index"][0], info["neighbor"]["index"][0][0])
    mol.RemoveBond(info["star"]["index"][1], info["neighbor"]["index"][1][0])
    mol.RemoveAtom(info["star"]["index"][1])
    mol.RemoveAtom(info["star"]["index"][0])
    if logging.DEBUG >= logging.root.level:
        nb_display(mol)

    # (4) Canonicalize
    # We have a cyclic molecule now that we can canonicalize
    # Because it is cyclic, it will always give the same SMILES
    logging.debug("(4) Canonicalize")
    Chem.Kekulize(mol, clearAromaticFlags=True)

    sm = Chem.MolToCXSmiles(mol)
    mol = Chem.RWMol(Chem.MolFromSmiles(sm))

    if logging.DEBUG >= logging.root.level:
        nb_display(mol)

    # (5) Find bond to break
    # Kekulize before breaking the bond (find conjugated bonds)
    Chem.Kekulize(mol, clearAromaticFlags=True)

    # Find the bond and get bond type
    b0_break, b1_break = get_index_to_break(mol, info)
    btype_removed = mol.GetBondBetweenAtoms(b0_break, b1_break).GetBondType()
    if btype_removed == Chem.rdchem.BondType.AROMATIC:
        conjugated = mol.GetBondBetweenAtoms(b0_break, b1_break).GetIsConjugated()
        if conjugated:
            btype_removed = Chem.rdchem.BondType.DOUBLE
        else:
            btype_removed = Chem.rdchem.BondType.SINGLE

    logging.debug(
        f"(5) Break bond at atom idx {b0_break} and {b1_break}. Bond type was {btype_removed}"
    )

    mol.RemoveBond(b0_break, b1_break)
    if logging.DEBUG >= logging.root.level:
        nb_display(mol)

    # (6) Add stars
    logging.debug(
        f"(6) Add stars at index {b0_break} and {b1_break} "
        f"with bond type {btype_removed}. "
        f"Renumber atoms."
    )
    idx = mol.AddAtom(Chem.AtomFromSmarts("*"))
    mol.AddBond(b0_break, idx, btype_removed)

    idx = mol.AddAtom(Chem.AtomFromSmarts("*"))
    mol.AddBond(b1_break, idx, btype_removed)

    # Renumber atoms
    sm = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(sm)

    if logging.DEBUG >= logging.root.level:
        nb_display(mol)

    return Chem.MolToSmiles(mol)


def get_index_to_break(mol: Chem.RWMol, info: Dict) -> Tuple[int, int]:
    """Get indices of two connected atoms that can be broken.
    Indices cannot be part of an original ring.

    Args:
        mol (Chem.RWMol): [description]
        info (Dict): [description]

    Returns:
        Tuple[int, int]: [description]
    """
    # Find index of N1 and N2

    n_idx = []
    for atom in mol.GetAtoms():
        if "atomLabel" in atom.GetPropsAsDict() and "N" in atom.GetProp("atomLabel"):
            n_idx.append(atom.GetIdx())

    # Get the original bonds to make sure we don't break them
    # (1) break N1-N2 bond
    # (2) find rings
    # (3) add bond back
    bnd = mol.GetBondBetweenAtoms(n_idx[0], n_idx[1]).GetBondType()

    mol.RemoveBond(n_idx[0], n_idx[1])
    # mol.ClearComputedProps()
    GetSymmSSSR(mol)
    original_atom_rings = mol.GetRingInfo().AtomRings()
    original_bond_rings = mol.GetRingInfo().BondRings()

    # Add bond back
    mol.AddBond(n_idx[0], n_idx[1], bnd)
    # mol.ClearComputedProps()
    GetSymmSSSR(mol)
    new_atom_rings = mol.GetRingInfo().AtomRings()
    new_bond_rings = mol.GetRingInfo().BondRings()

    # Find rings with N1 and N2 that have the computed length
    ring_length_to_break = len(info["neighbor"]["path"])

    rings_match, bond_rings_match = [], []
    for atom_ring, bond_ring in zip(new_atom_rings, new_bond_rings):
        if (
            set(atom_ring).issuperset(set(n_idx))
            and len(atom_ring) == ring_length_to_break
        ):
            rings_match.append(atom_ring)
            bond_rings_match.append(bond_ring)

    logging.debug(
        f"Found {len(rings_match)} rings that contain N1 and N2 and have length {ring_length_to_break}"
    )

    if len(rings_match) == 0:
        # Remove ring_length_to_break criteria
        rings_match, bond_rings_match = [], []
        for atom_ring, bond_ring in zip(new_atom_rings, new_bond_rings):
            if set(atom_ring).issuperset(set(n_idx)):
                rings_match.append(atom_ring)
                bond_rings_match.append(bond_ring)

        logging.debug(f"Found {len(rings_match)} rings that contain N1 and N2")

    def bond_in_other_rings(bond_index):
        # check rings that do not match
        # for _ring in bond_rings_no_match:
        # if bond_index in _ring:
        # return True
        # Check the original rings; they must stay intact
        for _ring in original_bond_rings:
            if bond_index in _ring:
                return True
        return False

    def get_break_idx(ring_to_break):
        # start to search with lower index first
        sorted_ring = sorted(ring_to_break)

        logging.debug(f"Searching in ring (sorted index): {sorted_ring}")
        for idx_1 in sorted_ring:

            # get neighbors of idx_1
            possible_neighbors = [
                x.GetIdx() for x in mol.GetAtomWithIdx(int(idx_1)).GetNeighbors()
            ]

            # Sort possible neighbors, rdkit does note guarantee the order
            possible_neighbors = sorted(possible_neighbors)

            for idx_2 in possible_neighbors:
                bond_idx = mol.GetBondBetweenAtoms(idx_1, idx_2).GetIdx()
                # print(idx_1, idx_2, bond_idx, bond_in_other_rings(bond_idx))

                # idx_2 cannot be in original rings but must be in sorted_ring
                if idx_2 in sorted_ring and not bond_in_other_rings(bond_idx):

                    logging.debug(
                        f"Index {idx_1} and {idx_2} are in ring and connected with bond {bond_idx}. "
                        f"Bond {bond_idx} is not in other rings. Break."
                    )
                    return idx_1, idx_2

        raise UserWarning("Canonicalization failed.")

    # We need to check all ring, sort by size
    logging.debug(f"Possible rings to break: {rings_match}")

    # TODO: There might be multiple rings that can be broken; do we need to sort here?
    # print(rings_match)
    rings_match = sorted(rings_match, key=sum)
    # print(rings_match)

    for ring_to_break in rings_match:
        break_idx_1, break_idx_2 = get_break_idx(ring_to_break)
        if break_idx_1 is not None and break_idx_2 is not None:
            # We found indices
            break

    # if ring has only 3 atoms break at N1 and N2
    if len(ring_to_break) == 3 and len(info["symbols"]) > 5:
        break_idx_1 = n_idx[0]
        break_idx_2 = n_idx[1]

    return break_idx_1, break_idx_2


def get_mol(psmiles) -> Chem.RWMol:
    """Returns a RDKit mol object.

    Note:
        In jupyter notebooks, this function draws the SMILES string

    Returns:
        Chem.MolFromSmiles: Mol object
    """
    return Chem.RWMol(Chem.MolFromSmiles(psmiles))


def nb_display(mol):
    """Helper function to display polymer for debug"""
    print(f"SMILES: {Chem.MolToCXSmiles(mol)}")
    try:
        display(mol)
    except:
        pass


def add_brackets(psmiles: str) -> str:
    # convert * to [*]
    stars_no_bracket = re.findall(r"(?<!\[)\*(?!\])", psmiles)
    if len(stars_no_bracket) == 2:
        psmiles = psmiles.replace("*", "[*]")
    return psmiles
