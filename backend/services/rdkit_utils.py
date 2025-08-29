from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, DataStructs
import numpy as np

def compute_fingerprints(mol):
    morgan = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    maccs = Chem.MACCSkeys.GenMACCSKeys(mol)
    return {
        "morgan": morgan,
        "maccs": maccs.ToBitString()
    }

def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return {
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotatable": Lipinski.NumRotatableBonds(mol),
        "qed": Descriptors.qed(mol)
    }