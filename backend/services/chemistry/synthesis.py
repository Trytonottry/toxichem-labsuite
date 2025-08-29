# backend/services/chemistry/synthesis.py
from rdkit import Chem
from rdkit.Chem import BRICS, Recap
from typing import Dict, List

def retrosynthesis_brics(smiles: str) -> List[Dict]:
    """
    Разбивает молекулу на фрагменты по BRICS-связям
    Полезно для анализа синтетической доступности
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []

    # BRICS-декомпозиция
    brics_bonds = list(BRICS.FindBRICSBonds(mol))
    fragments = []

    for bond in brics_bonds:
        frag_mol = Chem.FragmentOnBonds(mol, [bond[0]], addDummies=True)
        frag_smiles = Chem.MolToSmiles(frag_mol)
        fragments.append({
            "type": "BRICS",
            "bond_idx": bond[0],
            "fragments": frag_smiles.split("."),
            "comment": f"Breakable bond for synthesis"
        })

    return fragments

def recap_decomposition(smiles: str) -> List[Dict]:
    """
    Использует RECAP для генерации синтетических путей
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []

    tree = Recap.RecapDecompose(mol)
    routes = []

    for child in tree.children.values():
        routes.append({
            "level": "primary",
            "smiles": Chem.MolToSmiles(child.mol),
            "origin": "RECAP",
            "comment": "Possible building block"
        })

    return routes