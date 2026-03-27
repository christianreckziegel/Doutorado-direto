# Charm hadronization: radial profile of D0 mesons

## Project Overview
This project focuses on the analysis of a dataset, including signal extraction, background subtraction, data visualization, and efficiency estimation. The main tasks are as follows:

1. **Signal Extraction**: Implement a method to extract the signal from the dataset.
2. **Background Subtraction**: Apply the side-band technique to subtract the background from the dataset.
3. **Plotting and Comparison**: Create visualizations to compare the results of the signal extraction and background subtraction methods.
4. **General Investigations**: Perform additional analyses and investigations on the dataset.
5. **Efficiency Estimation**: Estimate the efficiencies of the various methods used in the project.

## Dataset
The dataset used in this project is denominated 
* Experimental dataset: HF_LHC22o_pass4_minBias_medium_2P3PDstar (HY)
* MC dataset: HF_LHC24d3a_2P3PDstar, HF_LHC24d3b_2P3PDstar
and are processed on Hyperloop generating different output files depending on which O2 table will be produced (which 
dependes on the process function executed in the task).
In order to run the O2 task locally, an `AO2D.root` type file is downloaded with its path specified on file `data.txt`.

## Tasks

### 0. O2 task
The task for analysing ALICE run 3 data called `hffragmenationfunction.cxx` and is located in the O2 open repository in `O2Physics/PWGJE/Tasks/`.
Depending on the process executed, different tables are produced. These tables can be mainly subdivided in collected experimental or simulated data.  

### 1. Signal Treatment
#### 1. Signal Extraction
Implement a method to extract the signal from the dataset.
The signal extraction task is located in the `signal_extraction/` directory. The main script is `signal_extraction.py`, which contains the implementation of the signal extraction method.

#### 2. Background Subtraction
Apply the side-band technique to subtract the background from the dataset.
The background subtraction task is located in the `background_subtraction/` directory. The main script is `background_subtraction.py`, which contains the implementation of the background subtraction method.

### 2. Efficiency Estimation
The efficiency estimation task is located in the `efficiency_estimation/` directory. The main script is `estimate_efficiencies.py`, which calculates the efficiencies of the various methods used in the project.

### 3. Detector Response

### 4. Detector Response

### 5. Systematic Uncertainties

## Usage
To run the different tasks, navigate to the respective directories and execute the corresponding scripts. For example:

