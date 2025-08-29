# tests/test_api.py
def test_predict_toxicity():
    response = client.post("/predict/toxicity", json={"smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"})
    assert response.status_code == 200
    data = response.json()
    assert "hERG" in data["predictions"]