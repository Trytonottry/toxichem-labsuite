# 🧪 ToxiChem LabSuite — Platform for Organic Toxin Research

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docker Build](https://img.shields.io/docker/cloud/build/toxichem/toxichem-labsuite)](https://hub.docker.com/r/toxichem/toxichem-labsuite)
[![Uptime](https://img.shields.io/badge/Uptime-99.9%25-brightgreen)](https://status.toxichem.org)
[![Release](https://img.shields.io/github/v/release/toxichem/toxichem-labsuite?include_prereleases)](https://github.com/Trytonottry/toxichem-labsuite/releases)
[![Documentation](https://img.shields.io/badge/Docs-GitBook-blue)](https://docs.toxichem.org)
![Code Style: Black](https://img.shields.io/badge/code%20style-black-000000.svg)
[![Type Checked: Mypy](https://img.shields.io/badge/type_checked-mypy-blue)](https://mypy-lang.org)
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://pre-commit.com)

> **ToxiChem LabSuite** — open-source platform for biochemists and toxicologists to study organic toxins, their mechanisms, metabolism, and toxicity. Combines cheminformatics, AI, and lab automation.

---

## 🌟 Features

- ✅ **Chemical Database** — Store, search, and analyze toxins (SMILES, InChI, properties)
- ✅ **3D Visualization** — Interactive molecule & protein-ligand viewer (3DMol.js)
- ✅ **Toxicity Prediction** — ML models for hERG, hepatotoxicity, mutagenicity
- ✅ **Metabolism & Synthesis** — Predict metabolic pathways (CYP450) and retrosynthesis (BRICS/RECAP)
- ✅ **Molecular Docking** — Integrated AutoDock Vina support
- ✅ **PBPK/PD Modeling** — Pharmacokinetic simulation
- ✅ **Pathway Analysis** — KEGG, Reactome integration
- ✅ **Offline Mode** — SQLite + Electron desktop app
- ✅ **PWA & Desktop** — One codebase, multiple platforms
- ✅ **GLP Reports** — PDF/Word report generation with templates
- ✅ **Monitoring** — Prometheus + Grafana metrics
- ✅ **Auto-updates** — Electron app with GitHub Releases

---

## 🚀 Quick Start

### Prerequisites
- Docker & Docker Compose
- Git

### Run with Docker

```bash
git clone https://github.com/Trytonottry/toxichem-labsuite.git
cd toxichem-labsuite
docker-compose up --build
```

Open: 

- 🔬 Web App: http://localhost:3000 
- 📊 API Docs: http://localhost:8000/docs 
- 📈 Metrics: http://localhost:9090  (Prometheus)
- 📊 Grafana: http://localhost:3001 
     
## 💻 Desktop App

```bash
cd desktop
npm install
npm start
```

## 🧪 API Examples

Predict Toxicity
```bash
curl -X POST http://localhost:8000/api/predict/toxicity \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"}'
```

Generate Metabolites
```bash
curl -X POST http://localhost:8000/metabolism/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"}'
```

Get PDF Report
```bash
curl -OJ http://localhost:8000/report/1.pdf
```

## 🛠️ Technologies

Backend - Python, FastAPI, RDKit, PostgreSQL (with RDKit), Celery
Frontend - React, Three.js, Plotly, 3DMol.js
Desktop - Electron, SQLite
ML - Scikit-learn, DeepChem, GNN
DevOps - Docker, Prometheus, Grafana, GitHub Actions
Reporting - Jinja2, WeasyPrint, python-docx

## 📄 License 

MIT © TryToNotTry
