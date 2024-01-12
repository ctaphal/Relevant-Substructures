from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import FunctionalGroups

def subunitIdent(smilesIn):

    # get BRICS decomposition substructures (of input molecule) in list1
    m = Chem.MolFromSmiles(smilesIn)
    list1 = list(BRICS.BRICSDecompose(m))

    # Overview of following code:
    # take functional groups in hierarchy (node objects) and get their SMARTS sequence values (an attribute of the node)
    # then use HasSubstructMatch function to compare input molecule with all of the standard substructures in the hierarchy
    # to identify substructure matches and append them to the BRICS decomposition substructures
    hierarchy = FunctionalGroups.BuildFuncGroupHierarchy()
    SMARTS = []
    for node in hierarchy:
        SMARTS.append(node.smarts)

    list2 = []
    for SMART in SMARTS:
        sub = Chem.MolFromSmarts(SMART)
        print(Chem.MolToSmiles(sub))
        isSubstruct = m.HasSubstructMatch(sub)
        print(isSubstruct)
        if (isSubstruct==True):
            list2.append(Chem.MolToSmiles(sub))

    return list1+list2

#print(subunitIdent('CN1CCC[C@H]1c2cccnc2'))
