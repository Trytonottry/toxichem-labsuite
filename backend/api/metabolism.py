# backend/api/metabolism.py
from fastapi import APIRouter
from pydantic import BaseModel

router = APIRouter(prefix="/metabolism", tags=["metabolism"])

class SMILESSchema(BaseModel):
    smiles: str

@router.post("/predict")
def predict_metabolism(req: SMILESSchema):
    from services.chemistry.metabolism import predict_metabolites
    return {"metabolites": predict_metabolites(req.smiles, max_steps=3)}