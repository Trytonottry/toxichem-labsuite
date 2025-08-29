# backend/api/pbpk.py
import json

def simulate_pbk(smiles: str, dose: float, route: str):
    # В реальности: вызов PK-Sim CLI или SimCyp
    # Здесь — заглушка
    time = np.linspace(0, 24, 100)
    concentration = np.exp(-0.1 * time) * dose  # Простая модель
    return {
        "time_h": time.tolist(),
        "plasma_conc": concentration.tolist(),
        "organ_levels": {"liver": 0.8, "kidney": 0.3, "brain": 0.1}
    }