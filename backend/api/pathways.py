# backend/api/pathways.py
import requests

def get_pathways_by_compound(kegg_id: str):
    url = f"http://rest.kegg.jp/get/{kegg_id}/json"
    return requests.get(url).json()

def get_reactome_pathways(smiles: str):
    # Поиск по химической структуре
    url = "https://reactome.org/ContentService/search/chemical/"
    resp = requests.get(url, params={"query": smiles})
    return resp.json()