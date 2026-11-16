Compressive sensing (here, CS) is a alternative to traditional sampling methods, using mathematical methods
to solve indeterminate linear systems from undersampled data [2]. Compressive sensing is a useful
approach for applications such as Magnetic Resonance Imaging (MRI) reconstruction, image com-
pression, and radar signals. This project will focus on MRI applications and image reconstruction.

FFT

Wavelets

How to Run

## Data

Download the fastMRI single-coil knee dataset from https://fastmri.med.nyu.edu/ and place the `.tar` file in `data/raw/` (gitignored). The function will extract it automatically on first use.

```matlab
addpath('src')
kspace = load_mri_data();           % single-coil, file 1, middle slice
kspace = load_mri_data(2, 10, true) % multi-coil, file 2, slice 10
```
