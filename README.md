# GEC_MAR
Github repository for Generalized Entropy Calibration under Missing at Random. This repository contains R scripts and data for the simulation study (Section 7.1) and the real data analysis (Section 7.2) in our paper.

## Contents
* **main.R** – Main script for the simulation study (Section 7.1).
* **nhanes_preprocess.R** – Preprocesses and imputes missing values in the NHANES dataset.
* **nhanes.RData** – Output data from `nhanes_preprocess.R` with imputed values.
* **nhanes.R** – Main script for the real data analysis (Section 7.2).
* **archive/** – Contains old or unused scripts retained for reference.

## Requirements
* **MOSEK** solver and the R package **Rmosek** must be installed to run the scripts.

## Usage
1. Execute `main.R` for the simulation study.
2. Execute `nhanes.R` for the real data analysis.
