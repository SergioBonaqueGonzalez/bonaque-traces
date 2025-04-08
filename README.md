# bonaque-traces
This repository contains the MATLAB scripts used to generate the figures in the paper:

**"Microthermal Residuals from Directed Electromagnetic Emissions as a New Class of Technosignature"**  
*Sergio Bonaque-González, 2025*

---

## 📄 Description

These scripts simulate the thermal signatures left by narrow beams of electromagnetic radiation interacting with various interstellar structures. The figures explore divergence geometries, temperature increases, and detectability thresholds.
For detailed methods and parameter definitions, see the Methods section of the accompanying paper.

### Contents

| Script                         | Description                                                  |
|-------------------------------|--------------------------------------------------------------|
| `generate_figure1_beam_divergence.m`    | Generates divergence plots for collimated and post-focus beams. |
| `BonaqueTrace_simulation_Fig_2_4_5.m.m`| Computes and visualizes temperature increases (ΔT) across ISM structures for different distances, durations, and powers. Generates Figures 2, 4, and 5. 
| `Figure3_SimulatedDetection.m`| Simulates beam detectability over noisy thermal backgrounds with spatial binning. Generates Figure 3. |

---

## 🛠 Requirements

- MATLAB R2016a or later.
- No toolboxes required beyond core MATLAB functions.
- Scripts are designed to be **self-contained** and reproducible, except for random noise in Figure 3.

---

## 🔁 Reproducibility Notes

- **Figure 3** includes random Gaussian noise. The result shown in the paper uses `rng(1)` for reproducibility, but minor visual differences may still appear due to randomness.
- Figures 2–5 include warnings if the parameters (e.g., duration) differ from the canonical ones used in the paper.

---

## 📚 Citation

If you use this code, please cite the corresponding paper:

> Bonaque-González, S. (2025). *Microthermal Residuals from Directed Electromagnetic Emissions as a New Class of Technosignature.* (manuscript under review; preprint available soon)

A Zenodo DOI will be provided once the preprint is available on arXiv.

---
## 📄 License

This repository is licensed under the **Creative Commons Attribution-NonCommercial 4.0 International** license.

You may reuse or adapt the code for academic and non-commercial purposes **as long as explicit attribution is given** to:

> **Sergio Bonaque-González (2025)**  
> *Microthermal Residuals from Directed Electromagnetic Emissions as a New Class of Technosignature*

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

---
For questions, contact: **sbonaque@ull.edu.es**

