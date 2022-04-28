# Paleodetectors Unicamp

This repository hosts the code being developed at the University of Campinas (Unicamp) to calculate the recoil rate of nuclei from Coherent Elastic Neutrino-Nucleus Ncattering (CE&nu;NS) in rocks. We expect to use the recoil rate to estimate measurable defects in rock samples.

In this first vesion of the code, we use reactors as our source of antineutrinos. Our idea is to use this accessible source for a test of principle and techniques. For our calculations, we use the Angra-II Brazilian nuclear reactor, as well as the distance of the neutrino laboratory placed near the facility.

## Usage

The code is structured in Python scripts as follows:

* [`xsc.py`](xsc.py): it contains the functions for the CE&nu;NS cross-section;
* [`Reactor_Flux.py`](Reactor_Flux.py): it contains the data and functions to calculate the Angra-II antineutrino flux at the laboratory;
* [`Recoil_Spectrum.py`](Recoil_Spectrum.py): it contains the the functions to calculate the recoil rate of nuclei in the materials being considered.

Some plots and physical descrtiption of our calculations, as well as references, can be found in the iPython notebook [`Paleodetectors.ipynb`](Paleodetectors.ipynb).

## Notes
 
 This code is still a work in progress.
  
 For the neutrino flux from other sources (solar, atmospheric, DSNB) we use the library WIMpy:
 ```
 B. J. Kavanagh and T. D. P. Edwards, WIMpy NREFT v1.0 [Computer Software], doi:10.5281/zenodo.1230503. Available at https://github.com/bradkav/WIMpy_NREFT, (2018)
```
