import numpy as np      # type: ignore
import xarray as xr     # type: ignore
import gcm_filters      # type: ignore


def find_closest_indices(i, j, ds):
    """
    Find the closest indices in a dataset based on given coordinates.

    Parameters:
    i (int or float): The target x-coordinate.
    j (int or float): The target y-coordinate.
    ds (xarray.Dataset): The dataset containing 'xi_rho' and 'eta_rho' variables.

    Returns:
    tuple: The indices (ii, jj) in the dataset closest to the given coordinates.
    """
    # Calculate the absolute difference between the target coordinates and the dataset coordinates
    absi = np.abs(ds.xi_rho - i)
    absj = np.abs(ds.eta_rho - j)

    # Compute the maximum of the absolute differences
    c = np.maximum(absi, absj)

    # Find the indices where the computed maximum difference is minimal
    (ii, jj) = np.where(c == np.min(c))

    return ii[0], jj[0]


def save_contour(contour, bath, pn, pm, outpath, outname):
    """
    Create a contour by finding the closest indices and save it as a netCDF file.

    Parameters:
    contour (numpy.ndarray): Array of contour points.
    bath (xarray.Dataset): Bathymetric dataset.
    pn (xarray.DataArray): Inverse metric 'pn' of the x-direction.
    pm (xarray.DataArray): Inverse metric 'pm' of the y-direction.
    outpath (str): Path where the netCDF file will be saved.
    outname (str): Name of the output netCDF file.

    Returns:
    None
    """
    # Extract lists of indices from the contour
    ilist = contour[:, 1]
    jlist = contour[:, 0]

    
    ii = []
    jj = []
    for i, j in zip(ilist, jlist):
        it, jt = find_closest_indices(i, j, bath)

        ii.append(it)
        jj.append(jt)

    ii.append(ii[0])
    jj.append(jj[0])

    
    ii = np.array(ii)
    jj = np.array(jj)
    
    
    if (np.sum(ii[:-1]-ii[1:]>1)>0) or (np.sum(jj[:-1]-jj[1:]>1)>0):
        print(np.sum(ii[:-1]-ii[1:]>1)>0)
        print(np.sum(jj[:-1]-jj[1:]>1)>0)
        print('Spacing more than 1 between indexes')
    
    indx = []
    jndx = []


    # dx and dy for u and v grid
    dx_u = []
    dy_v = []

    dx_v = []
    dy_u = []

    for c in range(len(ii)):
        i = ii[c]
        j = jj[c]

        ip = ii[c-1]
        jp = jj[c-1]

        if i == ip and j == jp:
            pass

        elif i != ip and j != jp:
            # go first x-direction
            indx.append(i)
            jndx.append(jp)

            dxu = 2/(pn.isel(xi_rho=i, eta_rho=jp)+pn.isel(xi_rho=i-1, eta_rho=jp)).values
            dxv = 2/(pn.isel(xi_rho=i, eta_rho=jp)+pn.isel(xi_rho=i, eta_rho=jp-1)).values

            # handle circulation direction, so that u dot dl gives right sign
            if i > ip:
                dx_u.append(dxu)

                dx_v.append(dxv)
            else:
                dx_u.append(-dxu)

                dx_v.append(-dxv)
            dy_v.append(0)
            dy_u.append(0)

            # ... and then y-direction
            indx.append(i)
            jndx.append(j)

            dyu = 2/(pm.isel(xi_rho=i, eta_rho=j)+pn.isel(xi_rho=i-1, eta_rho=j)).values
            dyv = 2/(pm.isel(xi_rho=i, eta_rho=j)+pn.isel(xi_rho=i, eta_rho=j-1)).values

            dx_u.append(0)
            dx_v.append(0)
            # handle circulation direction, so that u dot dl gives right sign
            if j > jp:
                dy_v.append(dyv)

                dy_u.append(dyu)
            else:
                dy_v.append(-dyv)

                dy_u.append(-dyu)
        else:
            indx.append(i)
            jndx.append(j)

            if i!=ip:

                dxu = 2/(pn.isel(xi_rho=i, eta_rho=j)+pn.isel(xi_rho=i-1, eta_rho=j)).values
                dxv = 2/(pn.isel(xi_rho=i, eta_rho=j)+pn.isel(xi_rho=i, eta_rho=j-1)).values

                if i > ip:
                    dx_u.append(dxu)

                    dx_v.append(dxv)
                else:
                    dx_u.append(-dxu)

                    dx_v.append(-dxv)
                dy_v.append(0)
                dy_u.append(0)
            else:

                dyu = 2/(pm.isel(xi_rho=i, eta_rho=j)+pn.isel(xi_rho=i-1, eta_rho=j)).values
                dyv = 2/(pm.isel(xi_rho=i, eta_rho=j)+pn.isel(xi_rho=i, eta_rho=j-1)).values   

                dx_u.append(0)
                dx_v.append(0)

                if j > jp:
                    dy_v.append(dyv)

                    dy_u.append(dyu)
                else:
                    dy_v.append(-dyv)

                    dy_u.append(-dyu)

    indx = np.array(indx)
    jndx = np.array(jndx)

    dx_u = np.array(dx_u)
    dy_v = np.array(dy_v)

    dx_v = np.array(dx_v)
    dy_u = np.array(dy_u)
    
    # save contour as netCDF-file
    step = np.arange(len(indx))

    ds_out = xr.Dataset(
            data_vars = dict(
                index_x = (['step'], indx),
                index_y = (['step'], jndx),
                dlx_u = (['step'], dx_u),
                dly_v = (['step'], dy_v),
                dlx_v = (['step'], dx_v),
                dly_u = (['step'], dy_u),

            ),
            coords = dict(step=step
            ),
        )
    
    ds_out.to_netcdf(outpath + outname)


def filter_bathymetry(ds, filter_scale):
    """
    Apply a Gaussian filter to the bathymetry data in the dataset.

    Parameters:
    ds (xarray.Dataset): The input dataset containing bathymetric data.
    filter_scale (float): The scale of the Gaussian filter.

    Returns:
    xarray.DataArray: The filtered bathymetry.
    """
    # Extract grid info centered at T-points
    wet_mask_t = ds.mask_rho
    dxT = 1 / ds.pn
    dyT = 1 / ds.pm
    area = dxT * dyT

    # Use T-point values for U-points and V-points
    dxCu = dxT
    dyCu = dyT
    dxCv = dxT
    dyCv = dyT

    # Compute minimum and maximum grid spacing
    dx_min = min(dxT.where(wet_mask_t).min(), dyT.where(wet_mask_t).min()).values
    dx_max = max(dxT.max(), dyT.max(), dxCu.max(), dyCu.max(), dxCv.max(), dyCv.max()).values

    # Initialize diffusivities
    kappa_iso = xr.ones_like(dxT)
    kappa_aniso = xr.zeros_like(dyT)
    kappa_w = xr.ones_like(dxCu)
    kappa_s = xr.ones_like(dxCu)

    # Prepare grid variables for the filter
    grid_vars_diff = {
        'wet_mask': wet_mask_t, 
        'dxw': dxCu, 'dyw': dyCu,
        'dxs': dxCv, 'dys': dyCv,
        'area': area, 'kappa_w': kappa_w, 'kappa_s': kappa_s
    }
    
    # Create the filter object
    filter_diff = gcm_filters.Filter(
        filter_scale=filter_scale,
        dx_min=dx_min,
        filter_shape=gcm_filters.FilterShape.GAUSSIAN,
        grid_type=gcm_filters.GridType.IRREGULAR_WITH_LAND,
        grid_vars=grid_vars_diff
    )

    # Apply the filter to the bathymetry data
    h_filtered = filter_diff.apply(ds.h.where(ds.mask_rho == 1), dims=['eta_rho', 'xi_rho'])
    
    return h_filtered

