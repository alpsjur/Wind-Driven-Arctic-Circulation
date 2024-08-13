import xarray as xr   # type: ignore
import numpy as np    # type: ignore
import xgcm           # type: ignore

import warnings
warnings.filterwarnings("ignore")

#contour = "norwegian_2750_100000_interpolate.nc"
contour = "canadian_2750_100000_interpolate.nc"

datapath = "data/"

outpath = datapath+"along_contour/"


# ROMS files
roms_files = sys.argv[0]


def rotate_to_contour(u, v, grid, C):
    # Convert to rho points and interpolate
    u_rho, v_rho = grid.interp(u, "X"), grid.interp(v, "Y")

    uc, vc = u_rho.interp(xi_rho=C.xi_rho, eta_rho=C.eta_rho), v_rho.interp(xi_rho=C.xi_rho, eta_rho=C.eta_rho)

    # calculate along-contour component
    Uc = (uc * C.tx + vc * C.ty)
    
    # calculate cross-component
    Vc = (uc * C.ty - vc * C.tx)

    return Uc, Vc

xchunk = -1
ychunk = -1
tchunk = 1

chunks = {"ocean_time":tchunk, 
          "xi_rho":xchunk, "xi_u":xchunk, 
          "eta_rho":ychunk, "eta_v":ychunk, 
          }

C = xr.open_dataset(datapath+f"contours/{contour}")


datasets = []
for i, file in enumerate(roms_files):
    ds = xr.open_dataset(file)

    # instead of xroms
    ds = ds.rename({'eta_u': 'eta_rho', 'xi_v': 'xi_rho', 'xi_psi': 'xi_u', 'eta_psi': 'eta_v'})
    ds = ds.chunk(chunks)

    coords={'X':{'center':'xi_rho', 'inner':'xi_u'}, 
        'Y':{'center':'eta_rho', 'inner':'eta_v'}, 
        'Z':{'center':'s_rho', 'outer':'s_w'}}

    grid = xgcm.Grid(ds, coords=coords, periodic=[])

    Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
    z_rho = Zo_rho* (ds.zeta + ds.h) + ds.zeta
    Zo_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
    z_w = Zo_w* (ds.zeta + ds.h) + ds.zeta

    ds.coords['z_w'] = z_w.where(ds.mask_rho, 0).transpose('ocean_time', 's_w', 'eta_rho', 'xi_rho')
    ds.coords['z_rho'] = z_rho.where(ds.mask_rho, 0).transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho')


    ds['pm_v'] = grid.interp(ds.pm, 'Y')
    ds['pn_u'] = grid.interp(ds.pn, 'X')
    ds['pm_u'] = grid.interp(ds.pm, 'X')
    ds['pn_v'] = grid.interp(ds.pn, 'Y')
    ds['pm_psi'] = grid.interp(grid.interp(ds.pm, 'Y'),  'X') # at psi points (eta_v, xi_u) 
    ds['pn_psi'] = grid.interp(grid.interp(ds.pn, 'X'),  'Y') # at psi points (eta_v, xi_u)

    ds['dx'] = 1/ds.pm
    ds['dx_u'] = 1/ds.pm_u
    ds['dx_v'] = 1/ds.pm_v
    ds['dx_psi'] = 1/ds.pm_psi

    ds['dy'] = 1/ds.pn
    ds['dy_u'] = 1/ds.pn_u
    ds['dy_v'] = 1/ds.pn_v
    ds['dy_psi'] = 1/ds.pn_psi

    ds['dz'] = grid.diff(ds.z_w, 'Z', boundary='fill')
    ds['dz_w'] = grid.diff(ds.z_rho, 'Z', boundary='fill')
    ds['dz_u'] = grid.interp(ds.dz, 'X')
    ds['dz_w_u'] = grid.interp(ds.dz_w, 'X')
    ds['dz_v'] = grid.interp(ds.dz, 'Y')
    ds['dz_w_v'] = grid.interp(ds.dz_w, 'Y')

    ds['dA'] = ds.dx * ds.dy

    metrics = {
        ('X',): ['dx', 'dx_u', 'dx_v', 'dx_psi'], # X distances
        ('Y',): ['dy', 'dy_u', 'dy_v', 'dy_psi'], # Y distances
        ('Z',): ['dz', 'dz_u', 'dz_v', 'dz_w', 'dz_w_u', 'dz_w_v'], # Z distances
        ('X', 'Y'): ['dA'] # Areas
    }


    grid = xgcm.Grid(ds, coords=coords, metrics=metrics, periodic=[])


    u, v = ds.u, ds.v
    ubar, vbar = ds.ubar, ds.vbar
    ub, vb = ds.u.isel(s_rho = 0), ds.v.isel(s_rho = 0)
    taux, tauy = ds.sustr, ds.svstr


    U, V = rotate_to_contour(u, v, grid, C)
    Ubar, Vbar = rotate_to_contour(ubar, vbar, grid, C)
    Ub, Vb = rotate_to_contour(ub, vb, grid, C)
    tauu, tauv = rotate_to_contour(taux, tauy, grid, C)


    # create dataset with extracted variables
    ds_out = xr.Dataset(
        data_vars = dict(
            Ubar = Ubar.astype("float32"),
            Vbar = Vbar.astype("float32"),
            Ub = Ub.astype("float32"),
            Vb = Vb.astype("float32"),
            U = U.astype("float32"),
            V = V.astype("float32"),
            taux = tauu.astype("float32"),
            tauy = tauv.astype("float32"),
            #distance = C.distance,
        ),
        coords = dict(
            ocean_time = ds.ocean_time.astype("float32"),
        )
    )

    datasets.append(ds_out)
    
# combine all time steps
combined = xr.concat(datasets, dim='ocean_time')

# save only contour-mean u(z) to save disc space
combined["U"] = combined["U"].mean("point")
combined["V"] = combined["V"].mean("point")


combined["distance"] = C.distance

# save as netCDF-file
combined.to_netcdf(outpath+f"U_at_{contour}")