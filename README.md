# Wind-Driven Time-Variable Ocean Circulation 🔎🌊

Welcome to this repositorty! 

Here, you will find data and code related to a study of the wind-driven, time-variable ocean circulation in the Arctic Mediterranean (that's the combined Arctic Ocean and Nordic Seas). Specifically, this repository contains data and code for calculating, analyzing and plotting time variable circulations around depth contours in: 
- Numerical simulations (ROMS-A4).
- A linear idealized model that estimates circulation driven by surface forcing.
  
 A good place to start exploring the code is the [create_figures.ipynb](create_figures.ipynb) notebook. 
 You can have a look at Sjur et al. (submitted) if you are interested in further details on the study.

## Contents
### [code](code)
- [calculate_circulations.py](calculate_circulations.py): Skript for extracting circulation timeseries from ROMS-A4 output.
- [create_contours.ipynb](create_contours.ipynb): Interactive notebook for filtering ROMS bathymetry and constructing depth contours used for calculating circulations in ROMS-A4. Contains plots illustrating the process. 
- [create_contours.py](create_contours.py): A more eficient way to cunstruct contours than the notebook version, but less informative. 
- [create_figures.ipynb](create_figures.ipynb): Notebook for creating figures. This is a good place to start exploring the code. 
- [extract_interpolated_contour.py](extract_interpolated_contour.py): Script for extracting ROMS-A4 model fields interpolated to a contour. 
- [linear_model.ipynb](linear_model.ipynb): Notebook for calculating linear model estimates. This notebook also contain information on the discretization of the linear model. 
- [linear_model.py](linear_model.py): A more efficient way to calculate linear estimates than the notebook version.
- [methods_create_contours.py](methods_create_contours.py): Contains methods used in [create_contours.ipynb](create_contours.ipynb) and [create_contours.py](create_contours.py).
- [methods_linear_model.py](methods_linear_model.py): Contains methoods used in [linear_model.ipynb](linear_model.ipynb) and [linear_model.py](linear_model.py). 
### [data](data)
- **along_contour**: Directory containing files with model fields interpolated to a given contour. Files are created by the script [extract_interpolated_contour.py](extract_interpolated_contour.py).
- **circulations**: Directory containing files with circulation time series extracted from ROMS-A4. Files are created by the script [calculate_circulations.py](calculate_circulations.py).
- **contours**: Directory containing files with contour data used to extract circulaiton time series from ROMS-A4. Files are created by either [create_contours.ipynb](create_contours.ipynb) or [create_contours.py](create_contours.py).
- **estimates**: Directory containing linear model estimates. Files are created by either [linear_model.ipynb](linear_model.ipynb) or [linear_model.py](linear_model.py).
- **full_domain**: Directory containing files with data covering the full ROMS-A4 domain, including the bathymetry and the time-mean velocity field. 


Additionally, the repository includes an `environment.yml` file for creating a conda environment with the necessary Python packages to run the code. To set up the environment, navigate to this directory in your terminal and run the following command:
```
conda env create -f environment.yml
```

Enjoy!