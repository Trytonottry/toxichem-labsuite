# tests/test_predict.py
def test_predict_toxicity(client):
    response = client.post("/api/predict/toxicity", json={
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    })
    assert response.status_code == 200
    data = response.json()
    assert "predictions" in data
    assert "hERG" in data["predictions"]
    assert "risk" in data["predictions"]["hERG"]