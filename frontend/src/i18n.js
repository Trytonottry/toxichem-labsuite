// frontend/src/i18n.js
import { useState, createContext, useContext } from 'react';

const translations = {
  en: {
    title: "ToxiChem LabSuite",
    search: "Search toxin by SMILES",
    predict: "Predict Toxicity",
    docking: "Run Docking",
    pathway: "Pathway Analysis"
  },
  ru: {
    title: "ToxiChem LabSuite",
    search: "Поиск токсина по SMILES",
    predict: "Предсказать токсичность",
    docking: "Запустить докинг",
    pathway: "Анализ путей"
  }
};

const LangContext = createContext();

export const useLang = () => useContext(LangContext);

export const LangProvider = ({ children }) => {
  const [lang, setLang] = useState('en');
  const t = (key) => translations[lang][key] || key;
  return (
    <LangContext.Provider value={{ lang, setLang, t }}>
      {children}
    </LangContext.Provider>
  );
};