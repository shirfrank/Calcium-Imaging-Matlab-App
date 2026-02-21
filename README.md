# OPC Calcium Imaging Analysis Tool
### Final Project in Neuroscience | Research Collaboration

This project was developed as a final graduation project in **Neuroscience**, serving as a specialized analysis pipeline for research conducted in collaboration with **Sophie Shohat's Master's thesis**.

## Overview
This repository contains a MATLAB-based suite designed for the automated analysis of calcium imaging data from Oligodendrocyte Precursor Cells (OPCs). The tool streamlines the transition from raw fluorescence intensity recordings to quantified biological parameters, supporting the investigation of cellular oscillations and signaling patterns.

## Academic Context
* **Field:** Neuroscience / Neurobiology.
* **Role:** Projectant Helper / Lead Developer of the Analysis GUI.
* **Collaboration:** Developed to support and analyze data for the Master's research of Sophie Shohat.
* **Objective:** To quantify OPC oscillation characteristics, including frequency, amplitude, and kinetics (rise/decay times).

## Key Features
* **Interactive MATLAB GUI:** Built with App Designer to allow researchers to load and analyze data without writing code.
* **Signal Pre-processing:** Automated $\Delta F/F$ calculation, smoothing (moving average), and baseline drift correction.
* **Advanced Peak Detection:** Specialized algorithms to distinguish between OPC oscillations and potential neuronal signal interference based on duration and prominence.
* **Detailed Feature Extraction:** * Peak Prominence & Amplitude.
    * Area Under the Curve (AUC).
    * Onset/Offset Latencies and Rise/Decay times.
    * Inter-peak Intervals (IPI).
* **Visual Sanity Checks:** Integrated plotting tools to verify the accuracy of the automated detection against raw traces.

## File Structure
* `Calcium_Imaging_GUI.mlapp`: The main interactive application interface.
* `analys_oscillations.m`: The core analysis engine containing the peak detection logic and feature calculations.
* `OPC_CaAnalysis_Preprocessing.mlx`: Documentation and workflow for initial data cleaning and normalization.
* `OPC_CaAnalysis_FeatureExtraction.mlx`: Batch processing script for large-scale data analysis across different experimental conditions (e.g., Mutant vs. Control).
* `Plot_CalciumTraces.m`: Utility script for high-resolution visualization of individual ROI signals.
* `model_params.mat`: Pre-configured parameters for the analysis models.

## Installation & Usage
1.  **Requirements:** MATLAB (R2021a or later) and the *Signal Processing Toolbox*.
2.  **Setup:** Clone the repository and add the folder to your MATLAB path.
3.  **Run:** Type `Calcium_Imaging_GUI` in the MATLAB Command Window or open the `.mlapp` file.
4.  **Process:** Load your `.xlsx` or `.mat` data, set the drug administration time (if applicable), and run the pre-processing and feature extraction tabs.

## Acknowledgments
Special thanks to **Sophie Shohat** for the collaboration and for providing the biological datasets and experimental context that guided the development of these tools.
