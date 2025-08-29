# backend/tox_predict.py
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from rdkit.Chem import Descriptors, MolFromSmiles
import joblib
import numpy as np

# Пример: обучение на синтетических данных
def featurize(smiles):
    mol = MolFromSmiles(smiles)
    if mol is None:
        return None
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Lipinski.HBondDonorCount(mol),
        Lipinski.HBondAcceptorCount(mol),
        Descriptors.TPSA(mol)
    ]

# В реальности: загрузка обученной модели
# model = joblib.load("tox_model.pkl")