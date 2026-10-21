# DROPLET-RISK: Unified Descriptor Atlas

DROPLET-RISK is a descriptor-based prioritization tool for identifying small molecules that may require leakage, oil-partition, interfacial-transfer, or material-sink controls in droplet microfluidic experiments.

Public site:

https://matteobroketa.github.io/Droplet-Risk-Atlas/

## Dataset scope

The figures and aggregate counts describe an offline analysis of **2,854,815 ChEMBL36 rows**, of which **2,784,347** had all seven descriptors required by the unified score.

The website does not download all 2.8 million rows. It provides:

1. **Instant local atlas:** a compact browser index built from the curated biological, Payne calibration, and ChEMBL priority/clinical records. The current source contains 16,382 records.
2. **Live ChEMBL fallback:** when the local index has no textual match, the visitor may explicitly search ChEMBL. After the visitor chooses the intended ChEMBL molecule or form, the browser retrieves its exact molecular properties and calculates the same score when all required values are present.
3. **Optional archived full export:** split `.jsonl.gz` release assets for reproducibility and bulk analysis. These are not loaded by the website.

A live score is labeled separately because it is based on the current API response and is not necessarily identical to the archived ChEMBL36 record used in the original database-wide analysis.

## Unified score

```text
score = 0.70·ALogP
      − 0.012·PSA
      − 0.18·HBD
      − 0.06·HBA
      − 0.0010·MW
      + 0.035·RTB
      + 0.10·aromatic rings

index = sigmoid(−2.1 + score)
```

Risk bands:

| Band | Score range |
|---|---:|
| Lower | below −0.2563230876636397 |
| Leaky | −0.2563230876636397 to below 0.1453605434241663 |
| High | 0.1453605434241663 to below 2.1 |
| Very high | 2.1 or above |
| Unscored | one or more required descriptors are missing |

The archived Python pipeline rounds the score to 10 decimal places for cutoff comparison. The browser implementation reproduces that behavior.

The score is a transparent ranking and triage proxy. It is not an absolute permeability or leakage measurement and does not replace validation in the exact oil, surfactant, droplet, pH, temperature, and device-material system.

## Website architecture

```text
site/index.html
site/assets/styles.css
site/assets/score.js
site/assets/app.js
site/data/lookup-16k.compact.json   # generated during deployment
site/data/atlas-summary.json        # generated during deployment
```

The compact browser file is generated from:

```text
data/lookup_records_unified.json
```

The builder removes duplicated `searchBlob` values and large `extra` dictionaries, omits standard interpretation strings by default, stores rows under a shared column schema, dictionary-encodes repeated strings, and trims trailing null values. The browser reconstructs only the search text needed for the approximately 16k local records after download.

## Local verification

Requirements:

```text
Python 3.10+
Node.js 18+
```

No pip or npm installation is required.

Run from the repository root:

```bash
python scripts/verify_drop_in.py
python -m http.server 8000 --directory site
```

Open:

```text
http://localhost:8000/
```

The verification also runs offline browser-logic regressions for descriptor formatting, Clear/example interactions, current-query CSV export, and delayed ChEMBL responses. Run these checks alone with `node tests/test_app.mjs`; no external API is contacted.

Optional live API diagnostic:

```bash
python scripts/check_chembl_api.py CHEMBL25
```

Do not open `site/index.html` directly with a `file://` URL because browser module and JSON requests may be blocked.

## GitHub Pages deployment

The workflow in `.github/workflows/pages.yml` builds the compact assets, runs model and site validation, and deploys `site/` on every push to `main`.

In repository settings, select:

```text
Settings → Pages → Build and deployment → Source → GitHub Actions
```

No ChEMBL API key or GitHub secret is required for the live fallback.

## Full archived row-level atlas

The optional exporter accepts the original full ChEMBL property CSV, CSV.GZ, or ZIP:

```bash
python scripts/export_full_atlas.py /path/to/chembl36_molecule_property_table_v2.zip
```

Outputs are written to `release-assets/` by default. Upload the generated `.jsonl.gz` parts and manifest to a GitHub Release. Do not copy them into `site/`.

Use [the local verification script](scripts/verify_drop_in.py) before deployment. The deployment and archived-export commands above describe the supported publishing paths.

## Empirical basis

The score vocabulary is anchored to the 36-compound droplet-transfer dataset reported by Payne et al.:

Payne, E. M. et al. *Analytical Chemistry* 2023, 95, 4662–4670. DOI: 10.1021/acs.analchem.2c04985.

## Author

Matteo Broketa — ORCID 0000-0001-7825-1174
