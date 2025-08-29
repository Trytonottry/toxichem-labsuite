# backend/services/chemistry/reactions.py
from rdkit import Chem
from rdkit.Chem import AllChem, rdRetroPaths

def predict_metabolism(smiles: str) -> list:
    """Простое предсказание метаболических реакций (окисление, глюкуронидирование)"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []

    # Пример: окисление ароматического кольца
    reactions = [
        {"type": "Oxidation", "product": "C1=CC=C(O)C=C1", "enzyme": "CYP3A4"},
        {"type": "Glucuronidation", "product": "C1=CC=C(OC2C(C(C(C2O)O)O)O)C=C1", "enzyme": "UGT1A1"}
    ]
    return reactions