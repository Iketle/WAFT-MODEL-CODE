# WAFT Model Code

This repository contains R code to reproduce the results in:

"A Semi-Parametric Weighted Least Squares Approach to Accelerated Failure Time Model: A Suitable Alternative to the Cox Model"

## Description

This project implements the Weighted Accelerated Failure Time (WAFT) model and compares it with the AFT Weibull model using simulation studies and real data analysis.

## Files

* waft_simulation.R
  Performs the simulation study for the WAFT model.

* waft_application.R
  Applies the WAFT model to the veteran dataset.

* aft_weibull_simulation.R
  Simulation study for the AFT Weibull model.

* aft_application.R
  Application of the AFT Weibull model.

## Requirements

* R version 4.4.1
* Required packages:
  survival, quantreg, KernSmooth, SurvRegCensCov, writexl

## How to run

1. Run `waft_simulation.R` to reproduce simulation results
2. Run `waft_application.R` for real data analysis

## Notes

* The dataset used is available in the `survival` package
* All scripts are fully reproducible
