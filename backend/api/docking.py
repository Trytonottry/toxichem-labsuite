# backend/api/docking.py
from fastapi import APIRouter, HTTPException
import subprocess
import os
from .schemas import DockingRequest

router = APIRouter()

@router.post("/docking/")
def run_docking(request: DockingRequest):
    # Скачиваем PDB (пример: 1A2C)
    # Генерируем PDBQT
    # Запускаем Vina
    try:
        result = subprocess.run([
            "vina", "--receptor", f"{request.pdb_id}.pdbqt",
            "--ligand", "ligand.pdbqt",
            "--center_x", "0", "--center_y", "0", "--center_z", "0",
            "--size_x", "20", "--size_y", "20", "--size_z", "20",
            "--out", "output.pdbqt"
        ], capture_output=True, text=True)
        return {"success": True, "output": result.stdout}
    except Exception as e:
        raise HTTPException(500, str(e))