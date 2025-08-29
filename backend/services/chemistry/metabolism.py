# backend/services/chemistry/metabolism.py
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from typing import List, Dict
import re

# Правила CYP450 (упрощённые)
METABOLIC_RULES = [
    {"enzyme": "CYP3A4", "reaction": "Aromatic oxidation", "pattern": "c1ccccc1", "product": "c1cc(O)cccc1"},
    {"enzyme": "CYP2D6", "reaction": "N-Dealkylation", "pattern": "N(C)C", "product": "N(C)O"},
    {"enzyme": "CYP1A2", "reaction": "O-Demethylation", "pattern": "COc1ccccc1", "product": "Oc1ccccc1"},
    {"enzyme": "UGT1A1", "reaction": "Glucuronidation", "pattern": "O", "product": "OC1C(C(C(C1O)O)O)O"},
]

def predict_metabolites(smiles: str, max_steps: int = 2) -> List[Dict]:
    """
    Предсказывает возможные метаболиты с указанием фермента и пути
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []

    metabolites = []
    queue = [(smiles, 0, "Parent compound", "—")]

    while queue:
        current_smiles, depth, reaction, enzyme = queue.pop(0)
        if depth >= max_steps:
            continue

        mol = Chem.MolFromSmiles(current_smiles)
        if not mol:
            continue

        for rule in METABOLIC_RULES:
            pattern = Chem.MolFromSmarts(rule["pattern"])
            if pattern and mol.HasSubstructMatch(pattern):
                # Простая замена (в реальности — реакции с помощью RDKit Chemical Reactions)
                product_smiles = current_smiles.replace(
                    re.escape(rule["pattern"]), rule["product"], 1
                ) if rule["pattern"] in current_smiles else current_smiles

                if product_smiles != current_smiles:
                    metabolites.append({
                        "smiles": product_smiles,
                        "reaction": rule["reaction"],
                        "enzyme": rule["enzyme"],
                        "generation": depth + 1,
                        "parent": current_smiles
                    })
                    if depth + 1 < max_steps:
                        queue.append((product_smiles, depth + 1, rule["reaction"], rule["enzyme"]))

    return metabolites