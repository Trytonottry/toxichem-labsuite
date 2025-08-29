// frontend/src/components/PBPKChart.jsx
import Plot from 'react-plotly.js';
import { useState, useEffect } from 'react';
import { api } from '../services/api';

export default function PBPKChart({ smiles, dose = 10 }) {
  const [data, setData] = useState(null);

  useEffect(() => {
    api.getPbPkSimulation(smiles, dose).then(setData);
  }, [smiles, dose]);

  if (!data) return <p>Simulating PK...</p>;

  return (
    <Plot
      data={[
        {
          x: data.time_h,
          y: data.plasma_conc_mg_L,
          type: 'scatter',
          mode: 'lines+points',
          marker: { color: 'blue' },
          name: 'Plasma Concentration'
        }
      ]}
      layout={{
        title: 'PBPK Simulation',
        xaxis: { title: 'Time (h)' },
        yaxis: { title: 'Concentration (mg/L)' }
      }}
      style={{ width: '100%', height: '400px' }}
    />
  );
}