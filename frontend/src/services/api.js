// frontend/src/services/api.js
const API_BASE = 'http://localhost:8000';

export const api = {
  async searchBySmiles(smiles) {
    const res = await fetch(`${API_BASE}/toxins/search/?smiles=${smiles}`);
    return res.json();
  },
  async predictToxicity(smiles) {
    const res = await fetch(`${API_BASE}/predict/toxicity`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles })
    });
    return res.json();
  },
  async runDocking(pdbId, smiles) {
    const res = await fetch(`${API_BASE}/docking/`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ pdb_id: pdbId, ligand_smiles: smiles })
    });
    return res.json();
  },
  async getPbPkSimulation(smiles, dose) {
    const res = await fetch(`${API_BASE}/pbpk/simulate`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles, dose, route: "oral" })
    });
    return res.json();
  }
};