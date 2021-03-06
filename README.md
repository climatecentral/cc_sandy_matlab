# Code Supporting "Economic Damages from Hurricane Sandy Attributable to Sea Level Rise Caused by Anthropogenic Climate Change" by Strauss et al. (2020)
 
This repository contains the code to compute much of the analysis presented in "Economic Damages from Hurricane Sandy Attributable to Sea Level Rise Caused by Anthropogenic Climate Change", which will soon be published at *Nature Communications* (April 2021).

If you have any questions, comments, or feedback on this code, please [contact Daniel Gilford (Climate Central)](mailto:dgilford@climatecentral.org) or open an [Issue](https://github.com/climatecentral/cc_sandy_matlab/issues) in the repository.

## Citation

If you use any or all of this code in your work, please include the citation:


> B. H. Strauss, P. Orton, K. Bittermann, M. K. Buchanan, D. M. Gilford, R. E. Kopp, S. Kulp, C. Massey, H. de Moel, S. Vinogradov, 2021: 
> Economic Damages from Hurricane Sandy Attributable to Sea Level Rise Caused by Anthropogenic Climate Change. 
> Nature Communications (in press, Apr. 2021).

## Data

The data supporting (and in some cases output by) these analyses is archived at Zenodo with the doi:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4543662.svg)](https://doi.org/10.5281/zenodo.4543662)

The code expects this data, which includes simulation water height fields, semi-empirical model analyses, and observations. If you are running these analyses yourself, please download the data and update the appropriate filepaths in the codes you are running/referencing.

## Functionality

This code was developed with MATLAB R2017b and R2018b.

Importantly, this repository includes the code to load and correct the ADCIRC simulations. To do this correction, run the scripts in order:

```
WarpSandy()
SandyCorrectionError()
```

This will load the table "sandyObservationData," and build the "sandyWarpedModelCompleteData" matrix variable, which contains latitude and longitude in the first two columns, and the corrected max water heights for simulations alt01-alt09 (alt06 being historical) in the remaining columns.

Errors can be assessed using the script:
```
EvaluateSandyError( sandyObservationData, sandyWarpedModelCompleteData, 9 )
```

This does not include data nor code to compute exposure of any variable, as these are part of Climate Central's proprietary analysis systems, which due to licensing restrictions, cannot be shared.

## Authors

### Primary Code Developers
* **Scott Kulp** - [GitHub](https://github.com/sckulp)
* **Daniel Gilford** - [GitHub](https://github.com/dgilford)

### Code Contributors
* **Klaus Bittermann**
* **Bob Kopp**
* **Ben Strauss**

### Additional Co-authors
* P. Orton
* M. K. Buchanan
* C. Massey
* H. de Moel
* S. Vinogradov

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

Benjamin Strauss, Scott Kulp, and Maya Buchanan were supported by the Kresge Foundation and the George and Estelle Sands Foundation. Benjamin Strauss, Robert Kopp, Scott Kulp, and Daniel Gilford were supported by NSF grant ICER-1663807 and NASA grant 80NSSC17K0698. Philip Orton was supported by NOAA award NA15OAR4310147. We thank Andrew Cox at Oceanweather Inc. for providing their meteorological reanalysis data for Hurricane Sandy, Chris Schubert, and Tom Suro (USGS) for providing data and interpretive guidance concerning maximum flood elevations. We thank Michael Oppenheimer for thoughtful comments on the manuscript. We acknowledge the World Climate Research Programmeõs Working Group on Coupled Modeling, which is responsible for CMIP, and we thank the climate modeling groups (listed in Supplementary Table 5) for producing and making available their model output. For CMIP, the U.S. DOE's Program for Climate Model Diagnosis and Intercomparison provides coordinating support and led the development of software infrastructure in partnership with the Global Organization for Earth System Science Portals.
