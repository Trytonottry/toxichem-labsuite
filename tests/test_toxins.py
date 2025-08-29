# tests/test_toxins.py
def test_create_toxin(client, test_db):
    response = client.post("/api/toxins/", json={
        "name": "Caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "ld50": 192
    })
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == "Caffeine"
    assert data["molecular_weight"] > 100
    assert data["logp"] < 1.0

def test_search_toxin(client, test_db):
    client.post("/api/toxins/", json={
        "name": "Nicotine",
        "smiles": "C1CCN(C1)C2=CC=CC=C2",
        "ld50": 50
    })
    response = client.get("/api/toxins/search/?smiles=C1CCN(C1)C2=CC=CC=C2")
    assert response.status_code == 200
    assert len(response.json()) > 0
    assert "Nicotine" in response.json()[0]["name"]