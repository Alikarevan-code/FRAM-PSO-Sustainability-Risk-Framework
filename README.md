# FRAM-PSO: A Semi-Quantitative Framework for Sustainability and Optimization

This repository contains the data, models, and MATLAB code for the manuscript titled: "A Semi-Quantitative Framework for Integrating Sustainability and Optimization into Systemic Risk Assessment using FRAM-PSO."

## Overview

This project provides a computational framework for integrating the Functional Resonance Analysis Method (FRAM) with Particle Swarm Optimization (PSO) to perform sustainability-driven risk assessment. The MATLAB script reproduces the analysis and figures as presented in the paper. The FRAM models for all three case studies are also provided.

## Repository Structure

- **`/FRAM_Models/`**: This directory contains the `.xfmv` files for the three case studies. These models can be viewed and edited using the FRAM Model Visualizer (FMV).
  - `case_study_1.xfmv`
  - `case_study_2.xfmv`
  - `case_study_3.xfmv`
- **`FRAM_PSO_Analysis.m`**: The main MATLAB script that performs the entire analysis for Case Study 1. It initializes the model, runs the iterative mitigation loop with the PSO, and generates the plots shown in the manuscript.
- **`LICENSE`**: The MIT License for this project.
- **`README.md`**: This file.

## How to Use This Repository

### 1. Viewing the FRAM Models
The `.xfmv` models in the `/FRAM_Models/` directory can be visualized using the free, web-based FRAM Model Visualizer (FMV).

1.  Go to https://functionalresonance.github.io/FMV_Community_Edition/
2.  Download the `.xfmv` file for the case study you wish to view from this repository.
3.  On the FMV website, upload the downloaded `.xfmv` file. You will then have access to the full FRAM model structure.

### 2. Reproducing the Analysis

To reproduce the quantitative analysis and figures for Case Study 1:

#### Prerequisites
-   MATLAB (Version R2022b or newer recommended)
-   No additional toolboxes are required.

#### Running the Script
1.  Clone or download this repository to your local machine.
2.  Open MATLAB.
3.  Navigate to the repository's root directory in MATLAB's "Current Folder" panel.
4.  Open the `FRAM_PSO_Analysis.m` script in the MATLAB editor.
5.  Click the **"Run"** button.

The script will execute all three phases of the analysis and generate the figures from the paper, including the Ranking Metric comparison and the Mitigation Progress plot.

## Contact

For any questions or issues with the code, please open an issue in this repository or contact the  author at: `ali.karevan.1@ens.etsmtl.ca`
