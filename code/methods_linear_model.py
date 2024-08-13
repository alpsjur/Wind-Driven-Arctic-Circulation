"""
See linear_model.ipynb for a description of the linear model.
"""

import xarray as xr     # type: ignore
import numpy as np      # type: ignore

def create_Ab(F, U0, R, H, dt):
    """
    Create matrix A and vector b for a given forcing time series.

    Parameters:
    F (numpy.ndarray): Forcing time series data.
    U0 (float): Initial value for the calculation.
    R (float): Friction parameter.
    H (float): Mean depth.
    dt (float): Time step.

    Returns:
    tuple: Tuple containing:
        - A (numpy.ndarray): The matrix A used for the calculation.
        - b (numpy.ndarray): The vector b used for the calculation.
    """
    

    n = len(F)
    dts = np.ones_like(F)*dt

    # Construct b 
    b = (np.cumprod(np.ones(n) * np.exp(-R * dts / H)) / np.exp(-R * dts / H)) * U0

    A = np.ones((n, n)) * np.exp(-R * dts / H)
    A = np.triu(A)
    np.fill_diagonal(A, 1)

    for i in range(n):
        A[i, i:] = np.cumprod(A[i, i:]) * dts[i]

        
    A[0,0] *= 0 
    
    return A, b


def run_linear_model(F, U0, Rs, H, dt):
    """
    Run a linear model for different friction parameters and return the estimates as an xarray.DataArray.

    Parameters:
    F (xarray.DataArray): Time series data with possible NaNs.
    U0 (float): Initial value for the circulation.
    Rs (list or numpy.ndarray): Array of friction parameters to evaluate.
    H (float): Mean depth.
    dt (float): Time step.

    Returns:
    xarray.DataArray: The estimated circulation for each friction parameter.
    """
    estimates = []

    for R in Rs:
        # Generate A and b matrices for solving matrix equation
        A, b = create_Ab(F.values, U0, R, H, dt)
     
        # Calculate the estimated values
        Uest = np.matmul(F.values, A) + b

        # Append the estimated values to the estimates list
        estimates.append(Uest)

    # Convert the estimates list to a numpy array
    estimates = np.array(estimates)

    # Create an xarray.DataArray from the estimates
    da = xr.DataArray(
        data=estimates,
        dims=["R", "time"],
        coords=dict(
            time=F.time,
            R=Rs,
        ),
    )

    return da

