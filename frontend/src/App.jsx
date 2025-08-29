// frontend/src/App.jsx
import { LangProvider, useLang } from './i18n';
import ToxPredictForm from './components/ToxPredictForm';
import MoleculeViewer from './components/MoleculeViewer';
import PBPKChart from './components/PBPKChart';

export default function App() {
  const { t } = useLang();
  const [smiles, setSmiles] = useState('CN1C=NC2=C1C(=O)N(C(=O)N2C)C'); // кофеин

  return (
    <LangProvider>
      <div style={{ padding: 20 }}>
        <h1>{t('title')}</h1>
        <input value={smiles} onChange={e => setSmiles(e.target.value)} />
        <MoleculeViewer smiles={smiles} />
        <ToxPredictForm />
        <PBPKChart smiles={smiles} dose={100} />
      </div>
    </LangProvider>
  );
}