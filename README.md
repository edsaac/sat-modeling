# sat-modeling

Supplementary code to the research article
**Continuum modeling of bioclogging of soil aquifer treatment systems segregating active and inactive biomass** by Edwin Saavedra Cifuentes, Alex Furman, Ravid Rosenzweig, and Aaron I. Packman (HESS-2024-251)

## Installation

1. Install OpenFOAM v.7.
2. Compile the auxiliary libraries from the `libs` folder using `wmake`.
3. Compile the `unsatRespiration` solver using `wmake`.
4. Check the solver compilation by running `unsatRespiration -help`.