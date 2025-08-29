// frontend/src/components/MoleculeViewer.jsx
import React, { useEffect, useRef } from 'react';
import { Viewer } from '3dmol';

const MoleculeViewer = ({ smiles }) => {
  const viewerRef = useRef();

  useEffect(() => {
    if (!smiles || !viewerRef.current) return;

    const config = { backgroundColor: 'white' };
    const viewer = new Viewer(viewerRef.current, config);

    fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${smiles}/JSON`)
      .then(res => res.json())
      .then(data => {
        const cid = data.PC_Compounds[0].id.id.cid;
        fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF`)
          .then(r => r.text())
          .then(sdf => {
            viewer.addModel(sdf, "sdf");
            viewer.setStyle({}, { stick: {}, sphere: { scale: 0.3 } });
            viewer.zoomTo();
            viewer.render();
          });
      });

  }, [smiles]);

  return <div ref={viewerRef} style={{ height: '400px', width: '100%' }} />;
};

export default MoleculeViewer;