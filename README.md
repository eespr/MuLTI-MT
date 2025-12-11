# MuLTI-MT
Multimodal Layered Transdimensional Inversion of Magnetotelluric data with Depth Constraints

The MuLTI MT algorithm was developed in Matlab version 2023b. This code can be run on Windows or Linux based platforms. This github repository includes: 1) the main MuLTI_MT.m code and 2) all functions needed to run this code. 

MuLTI_MT is a MATLAB script that runs a transdimensional MCMC inversion for 1-D MT data. The physical model consists of internal layers, each with an associated resistivity (R), sampled as log(R) defined by Voronoi nuclei. The domain is divided into a number of layers, num_layers (which could be one – unconstrained/no depth constraints applied) each with its own prior distribution on R. Each layer has a special nuclei that cannot leave its layer ("confined" nuclei). npt is the number of "floating" nuclei that can change layer. MuLTI MT Inverts the MT data using complex impedances (Z_real, Z_imag) but plots the apparent resistivity/phase. The code inverts either the determinant (Zssq; Rung-Arunwan et al., 2017), ZXY or ZYX.

Loading MT data: 
The function “load_data_edi.m” takes a given EDI filename and loads it into the standard data structure format (SI unit impedance with e^{+iwt} convention). © 2020 Martyn Unsworth Research Group (University of Alberta, Edmonton, Canada)

MuLTI MT core functions: 
“run_multiple_chains_Zinvert.m”, ‘thicknesses_and_priors’ and ‘whichnuclei’ are the core functions of MuLTI MT and based off Killingbeck et al. (2018).

1-D forward model solver:
“MT1D.m” computes the surface impedance Z for a 1D resistivity model. Written by Kerry Key, Scripps Institution of Oceanography.

Citation:
Siobhan Killingbeck. (2025). eespr/MuLTI-MT: MuLTI MT (v1.1). Zenodo. https://doi.org/10.5281/zenodo.17896670

References:
Killingbeck, S. F., Livermore, P. W., Booth, A. D., & West, L. J. (2018). Multimodal layered transdimensional inversion of seismic dispersion curves with depth constraints. Geochemistry, Geophysics, Geosystems, 19(12), 4957-4971.
Rung-Arunwan, T., Siripunvaraporn, W., & Utada, H. (2017). Use of ssq rotational invariant of magnetotelluric impedances for estimating informative properties for galvanic distortion. Earth, Planets and Space, 69(1), 80.
