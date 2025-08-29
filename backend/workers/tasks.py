# backend/workers/tasks.py
from celery import shared_task
import time
import json
import os

@shared_task
def run_docking_task(receptor_pdb: str, ligand_smiles: str):
    """Имитация запуска AutoDock Vina"""
    time.sleep(10)  # имитация вычислений
    return {
        "status": "completed",
        "affinity_kcal": -9.2,
        "output_pdbqt": f"generated_{receptor_pdb}_{ligand_smiles[:6]}.pdbqt"
    }

@shared_task
def simulate_pbk_task(smiles: str, dose: float, duration: int = 24):
    """Асинхронное PBPK-моделирование"""
    import numpy as np
    t = np.linspace(0, duration, 100)
    conc = dose * np.exp(-0.1 * t)
    return {
        "time_h": t.tolist(),
        "plasma_conc_mg_L": conc.tolist(),
        "half_life_h": 6.93,
        "auc": float(np.trapz(conc, t))
    }