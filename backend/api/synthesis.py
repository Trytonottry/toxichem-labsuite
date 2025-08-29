# backend/api/synthesis.py
from fastapi import APIRouter
from . import SMILESSchema

router = APIRouter(prefix="/synthesis", tags=["synthesis"])

@router.post("/retrosynthesis/brics")
def brics_retro(req: SMILESSchema):
    from services.chemistry.synthesis import retrosynthesis_brics
    return {"routes": retrosynthesis_brics(req.smiles)}

@router.post("/retrosynthesis/recap")
def recap_route(req: SMILESSchema):
    from services.chemistry.synthesis import recap_decomposition
    return {"routes": recap_decomposition(req.smiles)}