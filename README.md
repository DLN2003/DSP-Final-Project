# Octave Band Filtering: Lab P -14

## Overview
This project focuses on the implementation and analysis of octave band filters using FIR filters with various configurations. Both simple and enhanced filters (using a Hamming window) are explored, including their frequency responses, passband widths, and filter designs for different scenarios.

## Table of Contents
1. **Simple Bandpass Filter**  
   - Filter Design and Frequency Response Analysis  
   - Passband Width Calculation  
   - Effect of Filter Length (L) on Passband Width  
2. **Improved Bandpass Filter with Hamming Window**  
   - Filter Design  
   - Frequency Response Analysis at Specific Frequencies  
   - Passband Width Calculation for Different Filter Lengths  
3. **Application of Octave Band Filters**  
   - Filter Design for Specific Octaves  
   - Visualization of Octave Filter Responses  

---

## Implementation Details

### 1. Simple Bandpass Filter
- **Impulse Response Formula:**  
  $\ h(n) = \frac{2}{L} \cos(\omega_c n), \quad 0 \leq n < L\$  
- **Key Steps:**  
  - Generate and analyze a bandpass filter with $\( \omega_c = 0.4\pi \) and \( L = 40 \)$.  
  - Measure passband width at 50% level.  
  - Observe the effect of changing $\ L \$ (e.g., $\ L = 20 \$ and $\ L = 80 \$) on passband width.  

### 2. Improved Bandpass Filter (Hamming Window)
- **Impulse Response Formula:**  
  $\ h(n) = (0.54 - 0.46 \cos[\frac{2\pi n}{L-1}]) \cos(\omega_c [n - \frac{L-1}{2}]), \quad 0 \leq n < L \$  
- **Key Steps:**  
  - Design filter for $` \omega_c = 0.25\pi `$ and $\ L = 41 $.  
  - Analyze frequency response at specific frequencies $\( \omega = 0, 0.1\pi, 0.25\pi, 0.4\pi, 0.5\pi, 0.75\pi \)$.  
  - Compare passband widths for $\ L = 21 \$, $\ L = 41 \$, and $\ L = 81 \$.  

### 3. Application to Octave Band Filters
- **Goal:** Create bandpass filters for seven-octave bands.  
- **Design:** Determine filter lengths $\ L \$ to match specified bandwidths.  
- **Visualization:** Overlay the magnitude responses of all octave filters on a single plot.

---

## Functions Used
1. **`BPFsimp(wc, L, N)`**  
   - Generates a simple bandpass filter.  
2. **`PBWidth(H, w, threshold)`**  
   - Calculates passband width based on a specified magnitude threshold.  
3. **`BPFbetter(wc, L, N)`**  
   - Generates an improved bandpass filter using a Hamming window.  
4. **`HammingNorm(wc, L, N)`**  
   - Designs normalized Hamming-windowed filters for octave bands.

---

## Results
- **Effect of $\ L \$:**  
  Larger $\ L \$ results in narrower passbands with higher gain.  
- **Improved Filter:**  
  Hamming window reduces side lobes, providing a cleaner passband.  
- **Octave Band Filters:**  
  Successfully designed filters for seven octave bands, visualizing their responses.

---

## Running the Code
1. Ensure MATLAB/Octave is installed.  
2. Place all required files in the working directory, including:
   - `Bandpass_Filters.xlsx`
   - Custom functions: `BPFsimp.m`, `PBWidth.m`, `BPFbetter.m`, `HammingNorm.m`.  
3. Run the script to generate all figures and results.  

---

## Visual Outputs
The script generates:
1. Frequency response (magnitude and phase) plots for all filters.  
2. Passband widths displayed in the MATLAB/Octave command window.  
3. Overlayed responses for octave band filters.

---

## References
- Filter Design Theory  
- Hamming Window Applications  
- MATLAB/Octave Documentation

**Group 8 - Lab P -14**  
Authors:
- Giovanni Gutierrez
- Katie Henn
- Cade Boynton
- David Needens
