# backend/models.py
from sqlalchemy import Column, Integer, String, Float, Text
from database import Base

class Toxin(Base):
    __tablename__ = "toxins"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, index=True)
    smiles = Column(String, unique=True, index=True)
    molecular_weight = Column(Float)
    logp = Column(Float)
    ld50 = Column(Float)  # mg/kg
    category = Column(String)  # neurotoxin, hepatotoxin, etc.
    mechanism = Column(Text)