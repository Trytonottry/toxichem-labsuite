# backend/services/kegg_client.py
import requests

def search_kegg_compound(name: str) -> list:
    """Поиск соединения в KEGG"""
    url = f"http://rest.kegg.jp/find/compound/{name}"
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.strip().split("\n")
        return [{"id": line.split("\t")[0], "name": line.split("\t")[1]} for line in lines]
    return []

def get_kegg_pathways_by_compound(kegg_id: str) -> list:
    """Получает пути, в которых участвует соединение"""
    url = f"http://rest.kegg.jp/link/pathway/{kegg_id}"
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.strip().split("\n")
        return [{"pathway_id": line.split("\t")[1], "name": line.split("\t")[1].split(" ")[-1]} for line in lines]
    return []