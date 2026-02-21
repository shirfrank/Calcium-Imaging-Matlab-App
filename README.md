# OPC Calcium Imaging Analysis Tool
### Final Project in Neuroscience | Research Collaboration

This project was developed as a final graduation project in **Neuroscience**, serving as a specialized analysis pipeline for research conducted in collaboration with **Sophie Shohat's Master's thesis**.

![Main Interface](/images/main_gui.png)

*Figure 1: The landing screen and main interface of the Calcium Imaging Analysis GUI.*

## Overview
This repository contains a MATLAB-based suite designed for the automated analysis of calcium imaging data from Oligodendrocyte Precursor Cells (OPCs). The tool streamlines the transition from raw fluorescence intensity recordings to quantified biological parameters.

## Academic Context
* **Field:** Neuroscience / Neurobiology.
* **Role:** Projectant Helper / Lead Developer of the Analysis GUI.
* **Collaboration:** Developed to support and analyze data for the Master's research of Sophie Shohat.

## Key Features

### 1. Pre-processing Pipeline
![Preprocessing Tab](/images/feature_extraction.png)

*Figure 2: Pre-processing workflow including data loading, trimming, and ΔF/F calculation.*

* **Automated Normalization:** Calculates ΔF/F and applies smoothing (moving average).
* **Drift Correction:** Optional detrending to handle baseline instability.

### 2. Feature Extraction & Classification
![Feature Extraction](/images/main_guipreprocessing_tab.png)

*Figure 3: Automated identification of peaks and signal classification.*

* **Signal Categorization:** Automatically identifies "Low Activity", "Moderate Response", and "Suspected Neuronal Signals".
* **Kinetic Analysis:** Extracts Rise Time, Decay Time, AUC, and Inter-peak Intervals (IPI).

### 3. Visual Sanity Checks
![Sanity Check](/images/sanity_check.png)

*Figure 4: Validation interface to verify automated peak detection against raw traces.*

## File Structure
* `Calcium_Imaging_GUI.mlapp`: The main interactive application interface.
* `analys_oscillations.m`: The core analysis engine containing the peak detection logic.
* `OPC_CaAnalysis_Preprocessing.mlx`: Documentation for data cleaning.
* `OPC_CaAnalysis_FeatureExtraction.mlx`: Batch processing for large-scale datasets.
* `Plot_CalciumTraces.m`: Utility for high-resolution visualization.

## Installation & Usage
1.  **Requirements:** MATLAB (R2021a or later) and the *Signal Processing Toolbox*.
2.  **Run:** Open `Calcium_Imaging_GUI.mlapp` in MATLAB and click "Run".
3.  **Workflow:** Load `.xlsx`/`.mat` files -> Set parameters -> Extract Features -> Export Results.

## Acknowledgments
Special thanks to **Sophie Shohat** for the collaboration and for providing the biological datasets that guided the development of these tools.
