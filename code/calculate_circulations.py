import xarray as xr   # type: ignore
import numpy as np    # type: ignore
import os
import sys

# path where contour files are stored. 
contour_filepath = '../data/contours/'
contour_files = os.listdir(contour_filepath)


# roms files
roms_files = sys.argv[0]

# where calculated circulations will be stored

out_path = '../data/circulations/'


#initialize dictionary to store results in
results = {} 
for cfile in contour_files:
    # open contour file 
    dsc = xr.open_dataset(contour_filepath+cfile)
    
    results[cfile] = dict(
    	time = [],
        circ_u = [],             # circulation of ubar 
        circ_ub = [],            # circulation of bottom velocity
        circ_us = [],            # circulation of surface velocity 
        circ_tau = [],           # circulation of surface stress
        )   
    
# loop over roms files, one for each time step
for rfile in roms_files:
    print('Starting ', rfile)
    
    # open roms file
    dsr = xr.open_dataset(rfile).squeeze()  

    # select depth-averaged velocities
    if 'ubar' and 'vbar' in dsr.keys():
        u = dsr.ubar
        v = dsr.vbar
    else:
    # calculate depth averaged velocities if not available
        hu = 0.5*(dsr.h[:,:-1]+dsr.h[:,1:])
        Su = ((dsr.hc * dsr.s_w + dsr.Cs_w * hu) / (dsr.hc + hu))
        
        hv = 0.5*(dsr.h[:-1,:]+dsr.h[1:,:])
        Sv = ((dsr.hc * dsr.s_w + dsr.Cs_w * hv) / (dsr.hc + hv))
        
        dSu = Su.values[1:]-Su.values[:-1]
        dSv = Sv.values[1:]-Sv.values[:-1]
        
        ubar = np.sum(dsr.u.values*dSu, axis=0)
        vbar = np.sum(dsr.v.values*dSv, axis=0)
        
        dsr['ubar'] = (('eta_u','xi_u'), ubar)
        dsr['vbar'] = (('eta_v','xi_v'), vbar)
        
        u = dsr.ubar
        v = dsr.vbar
    
    # bottom velocities
    ub = dsr.u[0]
    vb = dsr.v[0]
    
    # surface velocities
    ut = dsr.u[-1]
    vt = dsr.v[-1]
    
    # and surface stresses
    if 'sustr' and 'svstr' in dsr.keys():
        tx = dsr.sustr
        ty = dsr.svstr
    else: 
        # set variable to nan
    
        tx = np.full_like(u.values, np.nan, dtype=np.double)
        ty = np.full_like(v.values, np.nan, dtype=np.double)
        
        dsr['sustr'] = (('eta_u','xi_u'), tx)
        dsr['svstr'] = (('eta_v','xi_v'), ty)
        
        tx = dsr.sustr
        ty = dsr.svstr
    
    # loop over contours
    for cfile in contour_files:
        # open contour file with info about contour
        dsc = xr.open_dataset(contour_filepath+cfile)
        
        # read indecies (x, y) and lengths of contour segment
        # sign of dx and dy depend on direction of line segment
        x = dsc.index_x
        y = dsc.index_y
        
        dx_u = dsc.dlx_u
        dy_v = dsc.dly_v
        
        dx_v = dsc.dlx_v
        dy_u = dsc.dly_u
        
        # select variables at the contour
        us = u.isel(eta_u=y, xi_u=x-1)
        vs = v.isel(eta_v=y-1, xi_v=x)
        
        ubs = ub.isel(eta_u=y, xi_u=x-1)
        vbs = vb.isel(eta_v=y-1, xi_v=x)
        
        uts = ut.isel(eta_u=y, xi_u=x-1)
        vts = vt.isel(eta_v=y-1, xi_v=x)
        
        
        txs = tx.isel(eta_u=y, xi_u=x-1)
        tys = ty.isel(eta_v=y-1, xi_v=x)
    
      	
        # circulation 
        uc = np.sum((us*dx_u).values + (vs*dy_v).values)
        tauc = np.sum((txs*dx_u).values + (tys*dy_v).values)
        ubc = np.sum((ubs*dx_u).values + (vbs*dy_v).values) 
        utc = np.sum((uts*dx_u).values + (vts*dy_v).values)


        # find length of contour
        L = np.sum(np.abs(dx_u.values))+np.sum(np.abs(dy_v.values))

        # normalize circulation and fluxes with length of contour
        uc /= L
        tauc /= L
        ubc /= L
        utc /= L
   

        # store result
        results[cfile]['circ_u'].append(uc)
        results[cfile]['circ_ub'].append(ubc)
        results[cfile]['circ_us'].append(utc)
        results[cfile]['circ_tau'].append(tauc)

        results[cfile]['time'].append(np.datetime64(dsr.ocean_time.values))
        
    
    
# loop over results and save as a netCDF file
for cfile, values in results.items():   
    circ_u = np.array(values['circ_u'])
    circ_us = np.array(values['circ_us'])
    circ_ub = np.array(values['circ_ub'])
    circ_tau = np.array(values['circ_tau'])

    time = np.array(values['time'])
    
    ds_out = xr.Dataset(
            data_vars = dict(
                circ_u = (['time'], circ_u),
                circ_ub = (['time'], circ_ub),
                circ_us = (['time'], circ_us),
                circ_tau = (['time'], circ_tau),
            ),
            coords = dict(time = time
            ),
            attrs = dict(
                contour_file = cfile,
            )
        )
    
    outname = cfile[:-3]+'_circ.nc'
    ds_out.to_netcdf(out_path+outname)
