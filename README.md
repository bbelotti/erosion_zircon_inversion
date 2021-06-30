## _Zircon ages from suspended load as tracers for the inversion of subglacial erosion rates._
### Bruno Belotti - Master's thesis - UNIL GSE 2021


Contents of this repository are listed and described below:

* **Code_gorner-nonlin (folder):** Matlab code for the inversion of zircon ages adjusted to my data. The folder contains all files necessary to reproduce my thesis' results. Code from: [De Doncker, F., Herman, F., & Fox, M. (2020). Inversion of provenance data and sediment load into spatially varying  erosion  rates. Earth Surface Processes and Landforms, 45(15), 3879–3901](https://doi.org/10.1002/esp.5008). More info and instructions at: [https://github.com/fdedonck/Erosion_Inversion_ESPL](https://github.com/fdedonck/Erosion_Inversion_ESPL).

* **LA-ICPMS (folder):** Results of LA-ICPMS geochronology on zircon grains. I also made available the self-implemented code I used for the plotting of geochronology results (age signatures and concordia diagrams). Contents of this folder are listed below:
  * **'Plotting_geochronology_results.m'** Self-written code for the visualization of U-Pb geochronology results. Plots age signatures and concordia diagrams;
  * **'data (folder)'** Input for 'Plotting_geochronology_results.m';
  * **'IsoPlot_Bruno_2021.xlsx'** Reduced results of LA-ICPMS analyses on zircon grains (LAMTRACE format);
  * **'SEM_CL_images_DEMO (folder)'** Examples of Scanning Electron Microscope Cathodoluminescence images (.png) of analyzed zircon grains. For the full dataset (all zircon grains for all samples) please contact me (Bruno Belotti) (folder too heavy for GITHUB);
  * **'README_grain_ID_interpretation.txt'** Interpretation key for grain identification in CL images AND in 'IsoPlot_Bruno_2021.xlsx'.

* **Additional_data (folder):** Contains additional data used for the analysis and interpretation of my thesis' results (and other related data):
  * **'Chemical_Final.xlsx'** are the results of Zr concentration (and majors) analyses on bedrock samples (for zircon fertility);
  * **'rockeval_gornergletscher_2020.xlsx'** are Total Organic Carbon (TOC) concentrations analyses;
  * **'Sample_list_and_lithology.xlsx'** contains information on bedrock, sand and water samples used for the thesis (such as weight, number of extracted zircons, interpretation of bedrock geologies etc.;
  * **'SamplingResults_TurbidityCalibration_PumpSampler.xlsx'** are the results of turbidity measurements on the Gornera river. Data: [Prasicek, G., Mettra, F., Lane, S., & Herman, F.(2020). Recent patterns of discharge and sediment output of the gorner glacier, Switzerland. In EGU general assembly conference abstracts (p.9242).](2020EGUGA..22.9242P);
  * **'Discharge.mat'/'Precipitation.mat'/'Temperature.mat'** are for the Gornergletscher area and Gornera river during the study period (summer 2019). Data: [GIN database (Swiss Federal Office for the Environment and MeteoSwiss)](https://www.slf.ch/en/services-and-products/forecasting-and-warning/gin-the-common-natural-hazard-information-platform.html).
