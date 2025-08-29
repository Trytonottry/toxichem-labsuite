// frontend/src/components/ToxPredictForm.jsx
import React, { useState } from 'react';
import { api } from '../services/api';

export default function ToxPredictForm() {
  const [smiles, setSmiles] = useState('');
  const [result, setResult] = useState(null);

  const handleSubmit = async (e) => {
    e.preventDefault();
    const data = await api.predictToxicity(smiles);
    setResult(data);
  };

  return (
    <div>
      <h3>Predict Toxicity</h3>
      <form onSubmit={handleSubmit}>
        <input
          type="text"
          value={smiles}
          onChange={(e) => setSmiles(e.target.value)}
          placeholder="Enter SMILES"
          required
        />
        <button type="submit">Predict</button>
      </form>
      {result && (
        <pre>{JSON.stringify(result, null, 2)}</pre>
      )}
    </div>
  );
}