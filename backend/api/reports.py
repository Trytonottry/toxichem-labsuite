# backend/api/reports.py
from fastapi import APIRouter, Response
from reports.generator import generate_pdf_report

router = APIRouter(prefix="/report", tags=["reports"])

@router.get("/{toxin_id}.pdf")
def get_toxin_report(toxin_id: int):
    # Здесь — получение данных из БД
    toxin = {"name": "Test Toxin", "smiles": "C1=CC=CC=C1", "molecular_weight": 78.11, "logp": 2.13, "ld50": 100}
    predictions = {"hERG": {"risk": True, "confidence": 0.88}, "hepatotox": {"risk": False, "confidence": 0.76}}
    metabolites = [{"reaction": "Oxidation", "enzyme": "CYP3A4", "smiles": "Oc1ccccc1"}]

    pdf = generate_pdf_report(toxin, predictions, metabolites)
    return Response(content=pdf, media_type="application/pdf", headers={
        "Content-Disposition": f"attachment; filename=toxin_{toxin_id}_report.pdf"
    })