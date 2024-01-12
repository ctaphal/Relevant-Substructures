import statistics
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from rdkit.Chem.rdchem import RWMol
from rdkit.Chem.rdchem import GetPeriodicTable

# smile = 'CC(=O)N[C@@H](CO)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCCCN)C(=O)N1CCC[C@H]1C(=O)O'
# smile = 'C[C@@H](NC(=O)[C@@H](N)Cc1ccc(O)cc1)C(=O)N[C@@H](Cc1ccccc1)C(=O)NCC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CO)C(N)=O'
smile = 'CSCC[C@H](NC(=O)[C@H](Cc1ccccc1)n1c(=O)[nH]c2ccccc2c1=O)C(=O)N1CCC[C@@H]1C(=O)O'

def flatten(matrix):
    flat_list = []
    for row in matrix:
        flat_list.extend(row)
    return flat_list

def find_largest_repeated_submolecule(m):
    matches = []
    for i, bond in enumerate(m.GetBonds()):
        m_split = Chem.FragmentOnBonds(m, [bond.GetIdx()])
        fragments = Chem.GetMolFrags(m_split, asMols=True, sanitizeFrags=False)

        if len(fragments) != 2: continue
        m1, m2 = fragments

        m1_atoms = m1.GetNumAtoms()
        m2_atoms = m2.GetNumAtoms()
        # TODO: tune the constant here
        if max(m1_atoms, m2_atoms) / min(m1_atoms, m2_atoms) >= 5: continue

        mcs_result = Chem.MolFromSmarts(rdFMCS.FindMCS([m1, m2]).smartsString)
        # TODO: make chirality an input to the function
        matches_ = m.GetSubstructMatches(mcs_result, useChirality=True)

        if len(matches_) != 0:
            match = matches_[0]
            matches.append(match)

    if not matches:
        return None

    median_size = statistics.median(len(mat) for mat in matches)

    def rank_match(mat):
        size = len(mat)
        atomic_weights = [m.GetAtomWithIdx(a).GetMass() for a in mat]
        total_atomic_weight = sum(atomic_weights)
        average_atomic_weight = total_atomic_weight / len(atomic_weights)

        return size + average_atomic_weight

    # for match in matches:
    #     print(match, rank_match(match))

    best_match = max(matches, key=rank_match)

    match_mol = RWMol(m)
    idxs_to_remove = [
        atom.GetIdx() for atom in match_mol.GetAtoms() if atom.GetIdx() not in best_match
    ]

    for idx in sorted(idxs_to_remove, reverse=True):
        match_mol.RemoveAtom(idx)

    return (matches, match_mol)

def iterate_submolecules(mol):
    submolecules = [mol]
    highlights = []
    while True:
        result = find_largest_repeated_submolecule(submolecules[-1])
        if not result:
            break

        highlights.append(result[0])
        submolecules.append(result[1])

    # last submolecule gets no highlighting
    highlights.append([])

    return submolecules, highlights

def draw_submolecules(submolecules, highlights):
    for i, (m, h) in enumerate(zip(submolecules, highlights)):
        Draw.MolToFile(m, f"test{i}.png", size=(1000, 1000), highlightAtoms=flatten(h))

m = Chem.MolFromSmiles(smile)
draw_submolecules(*iterate_submolecules(m))
