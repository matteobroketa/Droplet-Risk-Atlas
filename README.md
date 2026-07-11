# DROPLET-RISK: Unified Descriptor Atlas

**DROPLET-RISK** is a specialized tool designed to help researchers identify and prioritize the risk of small-molecule leakage (crosstalk) in droplet microfluidic experiments. By scoring molecules through a unified proxy equation and a bounded index, this atlas provides actionable thresholds to predict when reporters, metabolites, or perturbagens might partition out of droplets.

---

The main tool is located in [DROPLET_RISK.html](https://matteobroketa.github.io/Droplet-Risk-Atlas/)

---

### Key Features
* **Unified Scoring:** Consistent descriptor pipeline applied to over **2.8 million ChEMBL36 records** and curated biological molecules.
* **Empirical Calibration:** Uses the 36-compound Payne ESI-MS dataset as a gold-standard calibration set.
* **Standardized Vocabulary:** Assigned risk thresholds and lookup terms for easy experimental planning.

### Data & Methodology
The tool's underlying measurements and logic are based on the work of **Payne et al.** ([DOI: 10.1021/acs.analchem.2c04985](https://doi.org/10.1021/acs.analchem.2c04985)). 
* **Original Data:** Sourced from ChEMBL; available upon request.
* **Descriptors:** Integrated from PubChem, LINCS, and OPERA/CompTox.

### Author
**Matteo Broketa** ORCID: [0000-0001-7825-1174](https://orcid.org/0000-0001-7825-1174)

---
*Keywords: droplet microfluidics, molecular crosstalk, ChEMBL, chemical descriptors, lookup database.*