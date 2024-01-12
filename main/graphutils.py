from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from rdkit.Chem.rdchem import RWMol
from rdkit.Chem.rdchem import GetPeriodicTable

def convertMolecule(molecule):
    return Chem.MolFromSmiles(molecule)

def flatten(matrix):
    flat_list = []
    for row in matrix:
        flat_list.extend(row)
    return flat_list

def average_tuples(*ts):
    return tuple(sum(es) / len(es) for es in zip(*ts))

def find_best_repeated_submolecule(m, original_molecule):
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
            matches.append(matches_)

    if not matches:
        return None

    # measures the proportion of the original molecule that the match covers
    def rank_match(mat):
        size = len(mat[0])
        num_matches = len(mat)

        return size

    best_match = max(matches, key=rank_match)

    # checks if the submolecule is too small
    if len(best_match[0]) <= 2:
        return None

    match_mol = RWMol(m)
    idxs_to_remove = [
        atom.GetIdx() for atom in match_mol.GetAtoms() if atom.GetIdx() not in best_match[0]
    ]

    for idx in sorted(idxs_to_remove, reverse=True):
        match_mol.RemoveAtom(idx)

    return (best_match, match_mol)

def iterate_submolecules(mol):
    submolecules = [mol]
    highlights = []
    while True:
        result = find_best_repeated_submolecule(submolecules[-1], mol)
        if not result:
            break

        highlights.append(result[0])
        submolecules.append(result[1])

    # last submolecule gets no highlighting
    highlights.append([])

    # manually add the special case of drawing the smallest found submolecule on the original
    smallest_submol = submolecules[-1]
    submolecules.append(mol)
    highlights.append(mol.GetSubstructMatches(smallest_submol))

    return submolecules, highlights

stock_colors = [(0,0,1,0.5), (0,1,0,0.5), (0,1,1,0.5), (1,0,0,0.5), (1,0,1,0.5), (1,1,0,0.5)]

def draw_submolecules(submolecules, highlights):
    images = []
    for i, (m, highlight) in enumerate(zip(submolecules, highlights)):
        drawer = Draw.rdMolDraw2D.MolDraw2DCairo(1000, 1000)

        colors = {}

        for j, match in enumerate(highlight):
            for atom in match:
                color = stock_colors[j % len(stock_colors)]
                if atom in colors:
                    colors[atom] = average_tuples(color, colors[atom])
                else:
                    colors[atom] = color

        drawer.DrawMolecule(
            m,
            highlightAtoms=flatten(highlight),
            highlightAtomColors=colors,
        )
        drawer.FinishDrawing()
        images.append(drawer.GetDrawingText())

    return images

def nth_smallest_submolecule(n, submolecules, highlights):
    if n >= len(submolecules):
        return None
    drawer = Draw.rdMolDraw2D.MolDraw2DCairo(1000, 1000)

    m = submolecules[0]
    submol = submolecules[-n-1]
    highlight = m.GetSubstructMatches(submol)

    drawer = Draw.rdMolDraw2D.MolDraw2DCairo(1000, 1000)

    colors = {}

    for j, match in enumerate(highlight):
        for atom in match:
            color = stock_colors[j % len(stock_colors)]
            if atom in colors:
                colors[atom] = average_tuples(color, colors[atom])
            else:
                colors[atom] = color

    drawer.DrawMolecule(
        m,
        highlightAtoms=flatten(highlight),
        highlightAtomColors=colors,
    )
    drawer.FinishDrawing()
    annotated_molecule_image = drawer.GetDrawingText()

    drawer = Draw.rdMolDraw2D.MolDraw2DCairo(1000, 1000)
    drawer.DrawMolecule(submol)
    drawer.FinishDrawing()
    submol_image = drawer.GetDrawingText()

    return Chem.MolToSmiles(submol), submol_image, annotated_molecule_image

if __name__ == "__main__":
    m = Chem.MolFromSmiles("CC(C)C[C@H](NC(=O)CNC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](C)NC(=O)[C@@H]1CCCN1C(=O)[C@H](CCCCN)NC(=O)CNC(=O)[C@@H]1CCCN1C(=O)[C@@H]1CCCN1C(=O)[C@@H]1CCCN1C(=O)[C@H](CCC(=O)O)NC(=O)CN)C(=O)N[C@H](C(=O)O)C(C)C")
    submol = nth_smallest_submolecule(2, *iterate_submolecules(m))


    print(submol[0])
    with open(f"submol.png", "wb") as f:
        f.write(submol[1])

    with open(f"annotated_molecule.png", "wb") as f:
        f.write(submol[2])
