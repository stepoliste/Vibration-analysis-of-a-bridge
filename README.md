# Bridge Dynamics Analysis

This repository contains a MATLAB script for performing structural dynamics and vibration analysis of a bridge model. The script covers tasks such as computing natural frequencies, mode shapes, and frequency response functions (FRFs), as well as analyzing the response to seismic excitations.

## Prerequisites

- MATLAB
- Bridge model data files (`ponte_mkr.mat` and `ponte_dampe_mkr.mat`)
- Seismic displacement data file (`seismic_displ.txt`)

## Usage

1. **Clone the repository:**

   ```bash
   git clone https://github.com/stepoliste/bridge-dynamics-analysis.git
   cd bridge-dynamics-analysis
   ```

2. **Run the script:**

   Open MATLAB, navigate to the repository directory, and run the script:

   ```matlab
   run('bridge_dynamics_analysis.m')
   ```

## Script Overview

The script performs the following main tasks:

1. **Initialization and Data Loading:**
   - Clears previous variables and figures.
   - Loads the bridge model data from the file `ponte_mkr.mat`.
   - Defines the number of degrees of freedom for the free-free and constrained parts of the structure.

2. **Matrices Extraction:**
   - Extracts submatrices from the mass (M), stiffness (K), and damping (R) matrices for free-free and constrained parts.

3. **Natural Frequencies and Modes:**
   - Computes the natural frequencies and mode shapes of the undamped system using eigenvalue analysis.
   - Sorts and organizes the eigenvalues and eigenvectors.

4. **Damped Natural Frequencies:**
   - Calculates the damped natural frequencies using specified damping coefficients (alpha and beta).

5. **Frequency Response Functions (FRF):**
   - Computes the frequency response of the system at specific nodes (node A and node B) for a range of frequencies.
   - Plots the magnitude and phase of the response.

6. **Mode Reduction and Frequency Response:**
   - Reduces the system to its first three modes.
   - Computes the frequency response for the reduced system.

7. **Bending Moment Calculation:**
   - Calculates the frequency response of the bending moment at a specified point (point C).

8. **Seismic Response Analysis:**
   - Loads seismic displacement data from `seismic_displ.txt`.
   - Computes the discrete Fourier transform (DFT) of the seismic data.
   - Computes and plots the frequency response of the bridge due to seismic excitation in terms of displacement and acceleration.

9. **Modified System Analysis:**
   - Loads a modified bridge model with damping from `ponte_dampe_mkr.mat`.
   - Computes and plots the frequency response of the vertical acceleration at point B for the modified system.

## Output

The script generates several plots that illustrate:

- The magnitude and phase of the frequency response at different nodes.
- The vertical displacement and acceleration due to earthquake excitation.
- The bending moment frequency response.

---

Feel free to contribute to this project by forking the repository and submitting pull requests. If you encounter any issues or have suggestions, please open an issue on GitHub.

---

**Note:** Ensure that the data files `ponte_mkr.mat`, `ponte_dampe_mkr.mat`, and `seismic_displ.txt` are in the repository directory before running the script.
