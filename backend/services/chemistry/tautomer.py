# backend/services/chemistry/tautomer.py
from rdkit import Chem
from rdkir.Chem import rdMolDescriptors

def generate_tautomers(smiles: str) -> list:
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []
    # Используем TautomerEnumerator
    from rdkit.Chem.MolStandardize.tautomer import TautomerEnumerator
    enumerator = TautomerEnumerator()
    tautomers = enumerator.Enumerate(mol)
    return [Chem.MolToSmiles(taut) for taut in set(tautomers)]