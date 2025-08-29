# backend/api/chemistry.py
from fastapi import APIRouter
from .schemas import SMILESSchema

router = APIRouter(prefix="/chemistry", tags=["chemistry"])

@router.post("/metabolism")
def predict_metabolism(req: SMILESSchema):
    from services.chemistry.metabolism import predict_metabolism
    return {"metabolites": predict_metabolism(req.smiles)}

@router.post("/tautomers")
def get_tautomers(req: SMILESSchema):
    from services.chemistry.tautomer import generate_tautomers
    return {"tautomers": generate_tautomers(req.smiles)}