from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import FunctionalGroups

def subunitIdent(smilesIn):
    m = Chem.MolFromSmiles(smilesIn)
    list1 = list(BRICS.BRICSDecompose(m))
    list2 = list(FunctionalGroups.CreateMolFingerprint(m, FunctionalGroups.BuildFuncGroupHierarchy()))
    all_subunits = list1 + list2
    return list2 #change so that it returns all_subunits after fixing code for li

print(subunitIdent('CCCOCc1cc(c2ncccc2)ccc1'))