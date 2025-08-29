# backend/api/predict.py
from fastapi import APIRouter
from pydantic import BaseModel
import joblib
import numpy as np
from services.rdkit_utils import calculate_properties
from fastapi import APIRouter, Body

router = APIRouter()

# Загружаем модели
models = {
    "hERG": joblib.load("ml_models/hERG_model.pkl"),
    "hepatotox": joblib.load("ml_models/hepato_model.pkl"),
    "ames": joblib.load("ml_models/ames_model.pkl")
}

class PredictRequest(BaseModel):
    smiles: str

@router.post("/predict/toxicity")
def predict_toxicity(req: PredictRequest):
    props = calculate_properties(req.smiles)
    if not props:
        return {"error": "Invalid SMILES"}

    X = np.array([[
        props["mw"], props["logp"], props["hbd"],
        props["hba"], props["tpsa"], props["rotatable"]
    ]])

    results = {}
    for name, model in models.items():
        pred = model.predict(X)[0]
        prob = model.predict_proba(X)[0].max()
        results[name] = {"risk": bool(pred), "confidence": float(prob)}

    return {"smiles": req.smiles, "predictions": results}

@router.post(
    "/predict/toxicity",
    summary="Predict key toxicity endpoints",
    description="""
    Predicts:
    - hERG cardiotoxicity
    - Hepatotoxicity
    - Mutagenicity (Ames test)
    
    Uses QSAR models trained on ChEMBL and ToxRefDB.
    """,
    responses={
        200: {
            "description": "Toxicity predictions with confidence",
            "content": {
                "application/json": {
                    "example": {
                        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                        "predictions": {
                            "hERG": {"risk": False, "confidence": 0.94},
                            "hepatotox": {"risk": True, "confidence": 0.87},
                            "ames": {"risk": False, "confidence": 0.91}
                        }
                    }
                }
            }
        }
    }
)
def predict_toxicity(req: PredictRequest = Body(..., example={
    "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
})):
    # ... (как раньше)
    return {"smiles": req.smiles, "predictions": results}