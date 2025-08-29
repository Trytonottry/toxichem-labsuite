import { use3DMol } from 'react-3dmol';

function ProteinLigandViewer({ pdbId, ligandSdf }) {
  const style = { backgroundColor: 'white' };
  const { ref, viewer } = use3DMol();

  useEffect(() => {
    if (viewer) {
      viewer.addModel(pdbId, "pdb");
      viewer.setStyle({ cartoon: { color: 'spectrum' } });
      viewer.addModel(ligandSdf, "sdf");
      viewer.setStyle({ stick: { colorscheme: 'greenCarbon' } }, { model: 1 });
      viewer.zoomTo();
      viewer.render();
    }
  }, [viewer, pdbId, ligandSdf]);

  return <div ref={ref} style={{ height: '500px' }} />;
}