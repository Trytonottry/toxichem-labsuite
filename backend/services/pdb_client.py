# backend/services/pdb_client.py
import requests

def fetch_pdb_structure(pdb_id: str) -> str:
    """Получает структуру белка в формате PDB"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise ValueError(f"PDB {pdb_id} not found")