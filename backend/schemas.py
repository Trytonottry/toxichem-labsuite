# backend/schemas.py
from pydantic import BaseModel
from typing import Optional

class ToxinBase(BaseModel):
    name: str
    smiles: str
    molecular_weight: float
    logp: float
    ld50: Optional[float] = None
    category: Optional[str] = None
    mechanism: Optional[str] = None

class ToxinCreate(ToxinBase):
    pass

class Toxin(ToxinBase):
    id: int

    class Config:
        from_attributes = True