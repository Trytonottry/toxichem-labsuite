# backend/reports/generator.py
from jinja2 import Template
from weasyprint import HTML
import json
from datetime import datetime

REPORT_TEMPLATE = """
<h1>ToxiChem Report: {{ toxin.name }}</h1>
<p><strong>Generated:</strong> {{ now }}</p>
<h2>Properties</h2>
<ul>
  <li><strong>SMILES:</strong> {{ toxin.smiles }}</li>
  <li><strong>MW:</strong> {{ "%.2f"|format(toxin.molecular_weight) }}</li>
  <li><strong>LogP:</strong> {{ "%.2f"|format(toxin.logp) }}</li>
  <li><strong>LD50:</strong> {{ toxin.ld50 }} mg/kg</li>
</ul>

<h2>Toxicity Prediction</h2>
{% for pred, data in predictions.items() %}
  <p>{{ pred }}: Risk={{ data.risk }}, Confidence={{ "%.2f"|format(data.confidence) }}</p>
{% endfor %}

<h2>Metabolites</h2>
<ul>
{% for meta in metabolites %}
  <li>{{ meta.reaction }} ({{ meta.enzyme }}) → {{ meta.smiles }}</li>
{% endfor %}
</ul>
"""

def generate_pdf_report(toxin_data: dict, predictions: dict, metabolites: list) -> bytes:
    template = Template(REPORT_TEMPLATE)
    html_string = template.render(
        toxin=toxin_data,
        predictions=predictions,
        metabolites=metabolites,
        now=datetime.now().isoformat()
    )
    html = HTML(string=html_string)
    return html.write_pdf()