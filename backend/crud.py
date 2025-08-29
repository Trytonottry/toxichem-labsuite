# backend/crud.py
from sqlalchemy.orm import Session
from . import models, schemas
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def create_toxin(db: Session, toxin: schemas.ToxinCreate):
    # Автоматически вычисляем свойства
    mol = Chem.MolFromSmiles(toxin.smiles)
    if not mol:
        raise ValueError("Invalid SMILES")

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)

    db_toxin = models.Toxin(
        name=toxin.name,
        smiles=toxin.smiles,
        molecular_weight=mw,
        logp=logp,
        ld50=toxin.ld50,
        category=toxin.category,
        mechanism=toxin.mechanism
    )
    db.add(db_toxin)
    db.commit()
    db.refresh(db_toxin)
    return db_toxin

def get_toxin(db: Session, toxin_id: int):
    return db.query(models.Toxin).filter(models.Toxin.id == toxin_id).first()

def get_toxins(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Toxin).offset(skip).limit(limit).all()

def search_toxin_by_smiles(db: Session, smiles: str):
    # Используем RDKit Tanimoto similarity
    query = f"""
    SELECT *
    FROM toxins
    WHERE morganbv % morgan_fp('{smiles}')
    ORDER BY morganbv <-> morgan_fp('{smiles}')
    LIMIT 10;
    """
    result = db.execute(query)
    return [dict(row) for row in result]