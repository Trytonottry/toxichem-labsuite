# backend/services/ml_model_loader.py
import joblib
import os
from typing import Dict
from sklearn.ensemble import RandomForestClassifier

class ModelLoader:
    _models: Dict[str, object] = {}

    @classmethod
    def load_model(cls, model_name: str) -> object:
        if model_name in cls._models:
            return cls._models[model_name]

        model_path = f"ml_models/{model_name}.pkl"
        if not os.path.exists(model_path):
            # Создаём заглушку (в реальности — обучение или скачивание)
            model = RandomForestClassifier(n_estimators=10)
            # Здесь можно добавить обучение на синтетических данных
            cls._models[model_name] = model
            return model

        cls._models[model_name] = joblib.load(model_path)
        return cls._models[model_name]

    @classmethod
    def predict(cls, model_name: str, X):
        model = cls.load_model(model_name)
        return model.predict(X), model.predict_proba(X)