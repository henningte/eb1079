
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Compendium of R code and data for “Prediction of Peat Properties from Transmission Mid-Infrared Spectra”

This repository contains the data and code for our manuscript:

Teickner, H. & Knorr, K.-H. (2025). “Prediction of peat Properties from
transmission mid-infrared spectra”. (unpublished).

### How to cite

Please cite this compendium as:

> Henning Teickner and Klaus-Holger Knorr, (2025). Compendium of R code
> and data for “Prediction of peat Properties from transmission
> mid-infrared spectra”. Accessed 26 Sep 2025.
> <https://github.com/henningte/eb1079>

### How to use

Instructions how to set up the Docker containers to reproduce the
computations are available from the Dockerfile. The Dockerfile also
provides instructions to run the
[`targets`](https://github.com/ropensci/targets) workflow to reproduce
all computations and eventually the manuscript and supporting
information.

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code :** [GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html)

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse. See the sources section for licenses for
data derived from external sources and how to give credit to the
original author(s) and the source.

### Sources

Data in folder `data/raw_data` are derived from different sources. To
use these data and give credit to data authors, please follow the
following information:

- [:file_folder: d9.rds](data/raw_data/d9.rds): This file contains data
  values (porosity and bulk density of organic shales) extracted from
  Fig. 4 in Wang et al. (2015).  
- [:file_folder:
  caldat-wittington2024](data/raw_data/caldat-wittington2024): This
  folder contains data values (bulk density, porosity, saturated
  hydraulic conductivity) extracted from Fig. 1 and 3 from Whittington
  and Koiter (2024).  
- [:file_folder: CHNOSZ](data/raw_data/CHNOSZ): This folder contains
  additional data for P$_4$O$_10$ to be included in the ‘OBIGT’ database
  from the R package ‘CHNOSZ’ (Dick 2019). These additional data values
  are derived from
  <https://webbook.nist.gov/cgi/cbook.cgi?ID=C16752606&Mask=6F>
  (Linstrom 1997).  
- [:file_folder:
  dry_thermal_conductivity](data/raw_data/dry_thermal_conductivity):
  This folder contains corrected versions of the supporting data from
  O’Connor et al. (2020) which I received via email from
  Prof. Cardenas.  
- [:file_folder: pmird](data/raw_data/pmird): This is a template folder
  in which the folder `pmird_prepared_data` from the pmird database
  needs to be stored to reproduce the computations. This folder is
  available from <https://doi.org/10.5281/zenodo.17092587>.
- [:file_folder:
  specific_heat_capacity](data/raw_data/specific_heat_capacity): This
  folder contains the supporting data from Gnatowski et al. (2022).

### Contributions

We welcome contributions from everyone. Please note that the eb1079
project is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

### Funding

This study was funded by the Deutsche Forschungsgemeinschaft (DFG,
German Research Foundation) grant no. KN 929/23-1 to Klaus-Holger Knorr
and grant no. PE 1632/18-1 to Edzer Pebesma.

### References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Dick.2019" class="csl-entry">

Dick, Jeffrey M. 2019. “CHNOSZ: Thermodynamic Calculations and Diagrams
for Geochemistry.” *Frontiers in Earth Science* 7 (July): 180.
<https://doi.org/10.3389/feart.2019.00180>.

</div>

<div id="ref-Gnatowski.2022" class="csl-entry">

Gnatowski, Tomasz, Ewa Ostrowska-Ligęza, Cedric Kechavarzi, Grzegorz
Kurzawski, and Jan Szatyłowicz. 2022. “Heat Capacity of Drained Peat
Soils.” *Applied Sciences* 12 (3): 1579.
<https://doi.org/10.3390/app12031579>.

</div>

<div id="ref-Linstrom.1997" class="csl-entry">

Linstrom, Peter. 1997. “NIST Chemistry WebBook, NIST Standard Reference
Database 69.” National Institute of Standards and Technology.
<https://doi.org/10.18434/T4D303>.

</div>

<div id="ref-OConnor.2020" class="csl-entry">

O’Connor, Michael T., M. Bayani Cardenas, Stephen B. Ferencz, Yue Wu,
Bethany T. Neilson, Jingyi Chen, and George W. Kling. 2020. “Empirical
Models for Predicting Water and Heat Flow Properties of Permafrost
Soils.” *Geophysical Research Letters* 47 (11): e2020GL087646.
<https://doi.org/10.1029/2020GL087646>.

</div>

<div id="ref-Wang.2015b" class="csl-entry">

Wang, Guochang, Yiwen Ju, Zhifeng Yan, and Qingguang Li. 2015. “Pore
Structure Characteristics of Coal-Bearing Shale Using Fluid Invasion
Methods: A Case Study in the Huainan–Huaibei Coalfield in China.”
*Marine and Petroleum Geology* 62 (April): 1–13.
<https://doi.org/10.1016/j.marpetgeo.2015.01.001>.

</div>

<div id="ref-Whittington.2024" class="csl-entry">

Whittington, Pete, and Alex Koiter. 2024. “Evaluation of Hydro-Physical
Properties Along a Northern Boreal Bog Peatland Transect.”
<https://doi.org/10.21203/rs.3.rs-4650224/v1>.

</div>

</div>
