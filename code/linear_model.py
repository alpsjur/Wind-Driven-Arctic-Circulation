"""
See linear_model.ipynb for a description of the linear model.
"""
import matplotlib.pyplot as plt
import glob
import xarray as xr
import scipy as sp
import os

from methods_linear_model import *

Rs = np.arange(2,151,2)*1e-5
rho = 1025
dt = 60*60*24

# path where contour files are stored. 
circulation_filepath = 'data/circulations/'
estimates_path = 'data/estimates/'


circulation_files = glob.glob(circulation_filepath+"*_circ.nc")


# calculate linear model estimates for each R value
for cfile in circulation_files:
    H = int(cfile.split("/")[-1].split("_")[1])
    
    ds = xr.open_dataset(cfile)
    F = ds.circ_tau/(rho*H)
    U0 = ds.circ_ub.isel(time=0).values

    estimates = run_linear_model(F, U0, Rs, H, dt)

    ds["estimates"] = estimates 

    ub = ds.circ_ub.fillna(0).values
    ubar = ds.circ_u.fillna(0).values

    # Calculate Pearson correlation coefficients for each R value
    rs_ub = []
    rs_ubar = []
    for R in ds.R:
        # chose y for correlation calculation. Could also be wind-driven or vorticity-driven components
        y = ds.sel(R=R).estimates.fillna(0).values
        r_ubar = np.corrcoef(ubar, y)[0,1]
        rs_ubar.append(r_ubar)
        
        r_ub = np.corrcoef(ub, y)[0,1]
        rs_ub.append(r_ub)

    # Store correlation coefficients in the dataset
    rs_ubar = np.array(rs_ubar)
    ds["r_ubar"] = ("R", rs_ubar)

    rs_ub = np.array(rs_ub)
    ds["r"] = ("R", rs_ub)

    outname = estimates_path+cfile[len(circulation_filepath):-7]+"estimates.nc"

    if os.path.exists(outname):
        os.remove(outname)

    ds.to_netcdf(outname)