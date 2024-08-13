import matplotlib.pyplot as plt
from skimage import measure
import cmocean

from methods_create_contours import *

# Read time mean fields for the full domain using xarray
filepath = "data/full_domain/"
filename = "A4mean.nc"

outpath = "data/contours/"

ds = xr.open_dataset(filepath+filename).squeeze()
pn = ds.pn 
pm = ds.pm

# select bathymetry
bath = ds.h.where(ds.mask_rho==1)


# create unfiltered contours for norwegian basin
print("Starting: Norwegian unfiltered contours")
depths = [2500, 2750, 3000, 3500]
ns = [0, 0, 0, 1]

for depth, n in zip(depths, ns):
    outname = f'norwegian_{depth}_0.nc'
    
    contours = measure.find_contours(bath.values, depth)
    
    #remove small contours
    contours = [contour for contour in contours if len(contour) > 100]

    contour = contours[n]
    save_contour(contour, bath, pn, pm, outpath, outname)
    

# create contours with various filter sizes for the 2750 m depth contour
# is used for checking filter sensitivity
print("Starting: Contours for filter sensitivity")
basins = ["norwegian", "greenland", "eurasian", "canadian"]
depth = 2750
filter_scales = [0, 20000, 40000, 60000, 80000, 120000, 140000, 160000, 180000, 200000]
ns = [ # filtering scale 
    [ 0,  0, 0, 0, 0, 0, 0, 0, 0, 0],   # norwegian
    [ 3,  3, 2, 2, 2, 3, 3, 3, 2, 2],   # greenland
    [ 6,  5, 4, 3, 3, 4, 4, 4, 3, 3],   # eurasian
    [11, 10, 7, 7, 5, 5, 5, 5, 4, 4]    # canadian
]

for j, filter_scale in enumerate(filter_scales):
    if filter_scale != 0:
        var = filter_bathymetry(ds, filter_scale).values
    else:
        var = bath.values
    
    contours = measure.find_contours(var, depth)
    
    #remove small contours
    contours = [contour for contour in contours if len(contour) > 100]
    
    for i, basin in enumerate(basins):
        n = ns[i][j]
        if np.isfinite(n):
            contour = contours[n]
            
            outname = f'{basin}_{depth}_{filter_scale}.nc'
            save_contour(contour, bath, pn, pm, outpath, outname)
            



# create contours for various depths
print("Starting: Contours for various depths")

basins = ["norwegian", "greenland", "eurasian", "canadian"]
filter_scale = 100000
depths = [2500, 2750, 3000, 3250, 3500, 3750]
ns = [                        # depth
    [     0, 0, 0, np.nan, 1, np.nan],   # norwegian
    [np.nan, 2, 2,      2, 2, np.nan],   # greenland
    [np.nan, 3, 3,      3, 4,      1],   # eurasian
    [np.nan, 5, 6,      7, 8,      4]    # canadian
]

var = filter_bathymetry(ds, filter_scale).values
for j, depth in enumerate(depths):
    contours = measure.find_contours(var, depth)
    
    #remove small contours
    contours = [contour for contour in contours if len(contour) > 100]
    
    for i, basin in enumerate(basins):
        n = ns[i][j]
        if np.isfinite(n):
            contour = contours[n]
            
            outname = f'{basin}_{depth}_{filter_scale}.nc'
            save_contour(contour, bath, pn, pm, outpath, outname)