# backend/services/pubchem.py
import requests
import json
from typing import Dict, Optional

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def get_compound_by_smiles(smiles: str) -> Optional[Dict]:
    """Получает данные о соединении по SMILES из PubChem"""
    try:
        # Поиск CID по SMILES
        response = requests.get(
            f"{PUBCHEM_URL}/compound/smiles/{smiles}/cids",
            params={"format": "json"}
        )
        if response.status_code != 200:
            return None
        data = response.json()
        cid = data.get("IdentifierList", {}).get("CID", [None])[0]
        if not cid:
            return None

        # Получение свойств
        props = requests.get(
            f"{PUBCHEM_URL}/compound/cid/{cid}/property/MolecularWeight,LogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount",
            headers={"Accept": "application/json"}
        ).json()

        # Получение структуры в SDF
        sdf = requests.get(
            f"{PUBCHEM_URL}/compound/cid/{cid}/SDF",
            headers={"Accept": "chemical/x-mdl-sdfile"}
        ).text

        return {
            "cid": cid,
            "properties": props.get("PropertyTable", {}).get("Properties", [{}])[0],
            "sdf": sdf
        }
    except Exception as e:
        print(f"PubChem error: {e}")
        return None