from astropy.io import fits
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline, RegularGridInterpolator


"""
the path_to_grid variable must be set to the filepath for the directory in which
the grids are stored such that the exposure time functions can find these files

ex: path_to_grid = '/home/anaconda3/lib/python3.6/site-packages/neid_calculator/'
"""

path_to_grid='/home/afgupta/neid_flask/neid_calculator/'

def NEID_RV_prec(teff, vmag, exptime, use_order=False, order=0):
    """
    calculate expected RV precision in cm/s for the given inputs
    teff:          Effective Temperature (K)
    vmag:          V-band magnitude
    exptime:       Exposure time (s)
    set use_order=True to calculate RV precision for a specific order
    """

    exptime_grid = fits.open(path_to_grid+'photon_grid_exptime.fits')[0].data
    teff_grid = fits.open(path_to_grid+'photon_grid_teff.fits')[0].data
    vmag_grid = fits.open(path_to_grid+'photon_grid_vmag.fits')[0].data
    logexp=np.log10(exptime_grid)
    
    bound_test=True
    if teff<np.min(teff_grid) or teff>np.max(teff_grid):
        print("Temperature out of bounds. The allowed range is %d K to %d K." % (np.min(teff_grid), np.max(teff_grid)))
        bound_test=False
    if vmag<np.min(vmag_grid) or vmag>np.max(vmag_grid):
        print("Magnitude out of bounds. The allowed range is V = %d to V = %d." % (np.min(vmag_grid), np.max(vmag_grid)))
        bound_test=False
    if exptime<np.min(exptime_grid) or exptime>np.max(exptime_grid):
        print("Exposure time out of bounds. The allowed range is %d s to %d s." % (np.min(exptime_grid), np.max(exptime_grid)))
        bound_test=False
    if bound_test==False:
        return np.nan

    
    if use_order==True:
        order_grid = fits.open(path_to_grid+'order_wvl_centers.fits')[0].data[0]
        order_loc=np.where(order_grid==order)[0][0]
        rvprec_grid_order = fits.open(path_to_grid+'dv_uncertainty_master_order.fits')[0].data
        rvprec_grid=rvprec_grid_order[order_loc]
    else:
        rvprec_grid = fits.open(path_to_grid+'dv_uncertainty_master.fits')[0].data
    
    
    teff_index=InterpolatedUnivariateSpline(teff_grid, 
                            np.arange(len(teff_grid), dtype=np.double))(teff)

    vmag_index=InterpolatedUnivariateSpline(vmag_grid, 
                            np.arange(len(vmag_grid), dtype=np.double))(vmag)
    exptime_index=InterpolatedUnivariateSpline(logexp, np.arange(len(exptime_grid),
                                                dtype=np.double))(np.log10(exptime))
    rvprec_interpolator=RegularGridInterpolator((np.arange(len(exptime_grid)),
                          np.arange(len(vmag_grid)),
                          np.arange(len(teff_grid))),
                                rvprec_grid)
    inputs=[exptime_index, vmag_index, teff_index]
    rv_precision=rvprec_interpolator(inputs)[0]

    return rv_precision

def NEID_exptime_RV(teff, vmag, rv_precision,use_order=False, order=0):
    """
    calculate exposure time required to achieve specified precision for given inputs
    teff:          Effective Temperature (K)
    vmag:          V-band magnitude
    rv_precision:  Desired Radial Velocity Precision (m/s)
    set use_order=True to calculate exposure time for a specific order
    """
    
    
    exptime_grid = fits.open(path_to_grid+'photon_grid_exptime.fits')[0].data
    teff_grid = fits.open(path_to_grid+'photon_grid_teff.fits')[0].data
    vmag_grid = fits.open(path_to_grid+'photon_grid_vmag.fits')[0].data
    logexp=np.log10(exptime_grid)
    
    bound_test=True
    if teff<np.min(teff_grid) or teff>np.max(teff_grid):
        print("Temperature out of bounds. The allowed range is %d K to %d K." % (np.min(teff_grid), np.max(teff_grid)))
        bound_test=False
    if vmag<np.min(vmag_grid) or vmag>np.max(vmag_grid):
        print("Magnitude out of bounds. The allowed range is V = %d to V = %d." % (np.min(vmag_grid), np.max(vmag_grid)))
        bound_test=False
    if bound_test==False:
        return np.nan
    
    if use_order==True:
        order_grid = fits.open(path_to_grid+'order_wvl_centers.fits')[0].data[0]
        order_loc=np.where(order_grid==order)[0][0]
        rvprec_grid_order = fits.open(path_to_grid+'dv_uncertainty_master_order.fits')[0].data
        rvprec_grid=rvprec_grid_order[order_loc]
    else:
        rvprec_grid = fits.open(path_to_grid+'dv_uncertainty_master.fits')[0].data

    teff_index=InterpolatedUnivariateSpline(teff_grid, 
                            np.arange(len(teff_grid), dtype=np.double))(teff)

    vmag_index=InterpolatedUnivariateSpline(vmag_grid, 
                            np.arange(len(vmag_grid), dtype=np.double))(vmag)        
    j=0
    eta=1e10
    while eta>rv_precision:
        exptime=2*(j+6)
        if exptime>np.max(exptime_grid):
            print("\nMaximum Exposure Time Exceeded (t>3600s).\n")
            return np.nan
        exptime_index=InterpolatedUnivariateSpline(logexp, np.arange(len(exptime_grid),
                                                dtype=np.double))(np.log10(exptime))
        rvprec_interpolator=RegularGridInterpolator((np.arange(len(exptime_grid)),
                              np.arange(len(vmag_grid)),
                              np.arange(len(teff_grid))),
                                    rvprec_grid)
        inputs=[exptime_index, vmag_index, teff_index]
        eta=rvprec_interpolator(inputs)[0]
        j+=1
    
    return exptime


def NEID_SNR(teff, vmag, exptime, wavelength):
    """
    calculate expected SNR for the given inputs
    teff:          Effective Temperature (K)
    vmag:          V-band magnitude
    exptime:       Exposure time (s)
    wavelength:    wavelength (nm) at which SNR should be calculated
    """
    
    exptime_grid = fits.open(path_to_grid+'photon_grid_exptime.fits')[0].data
    teff_grid = fits.open(path_to_grid+'photon_grid_teff.fits')[0].data
    vmag_grid = fits.open(path_to_grid+'photon_grid_vmag.fits')[0].data
    logexp=np.log10(exptime_grid)
    
    bound_test=True
    if teff<np.min(teff_grid) or teff>np.max(teff_grid):
        print("Temperature out of bounds. The allowed range is %d K to %d K." % (np.min(teff_grid), np.max(teff_grid)))
        bound_test=False
    if vmag<np.min(vmag_grid) or vmag>np.max(vmag_grid):
        print("Magnitude out of bounds. The allowed range is V = %d to V = %d." % (np.min(vmag_grid), np.max(vmag_grid)))
        bound_test=False
    if exptime<np.min(exptime_grid) or exptime>np.max(exptime_grid):
        print("Exposure time out of bounds. The allowed range is %d s to %d s." % (np.min(exptime_grid), np.max(exptime_grid)))
        bound_test=False
    if bound_test==False:
        return np.nan
    
    wavelength_grid = fits.open(path_to_grid+'order_wvl_centers.fits')[0].data[1]
    order_loc=np.where(np.abs(wavelength_grid-wavelength)<0.1)[0][0]
    snr_grid_order=fits.open(path_to_grid+'snr_master_order.fits')[0].data
    snr_grid=snr_grid_order[order_loc]

    
    teff_index=InterpolatedUnivariateSpline(teff_grid, 
                            np.arange(len(teff_grid), dtype=np.double))(teff)

    vmag_index=InterpolatedUnivariateSpline(vmag_grid, 
                            np.arange(len(vmag_grid), dtype=np.double))(vmag)
    exptime_index=InterpolatedUnivariateSpline(logexp, np.arange(len(exptime_grid),
                                                dtype=np.double))(np.log10(exptime))
    snr_interpolator=RegularGridInterpolator((np.arange(len(exptime_grid)),
                          np.arange(len(vmag_grid)),
                          np.arange(len(teff_grid))),
                                snr_grid)
    inputs=[exptime_index, vmag_index, teff_index]
    snr=snr_interpolator(inputs)[0]
    
    return snr
    
def NEID_exptime_SNR(teff, vmag, snr, wavelength):
    """
    calculate exposure time required to achieve specified SNR for given inputs
    teff:          Effective Temperature (K)
    vmag:          V-band magnitude
    snr:           Desired SNR
    wavelength:    wavelength (nm) at which exposure time should be calculated
    """
    
    exptime_grid = fits.open(path_to_grid+'photon_grid_exptime.fits')[0].data
    teff_grid = fits.open(path_to_grid+'photon_grid_teff.fits')[0].data
    vmag_grid = fits.open(path_to_grid+'photon_grid_vmag.fits')[0].data
    logexp=np.log10(exptime_grid)
    
    bound_test=True
    if teff<np.min(teff_grid) or teff>np.max(teff_grid):
        print("Temperature out of bounds. The allowed range is %d K to %d K." % (np.min(teff_grid), np.max(teff_grid)))
        bound_test=False
    if vmag<np.min(vmag_grid) or vmag>np.max(vmag_grid):
        print("Magnitude out of bounds. The allowed range is V = %d to V = %d." % (np.min(vmag_grid), np.max(vmag_grid)))
        bound_test=False
    if bound_test==False:
        return np.nan
    
    wavelength_grid = fits.open(path_to_grid+'order_wvl_centers.fits')[0].data[1]
    order_loc=np.where(np.abs(wavelength_grid-wavelength)<0.1)[0][0]
    snr_grid_order=fits.open(path_to_grid+'snr_master_order.fits')[0].data
    snr_grid=snr_grid_order[order_loc]
        
    
    teff_index=InterpolatedUnivariateSpline(teff_grid, 
                            np.arange(len(teff_grid), dtype=np.double))(teff)

    vmag_index=InterpolatedUnivariateSpline(vmag_grid, 
                            np.arange(len(vmag_grid), dtype=np.double))(vmag)  
    j=0
    eta=0
    while eta<snr:
        exptime=2*(j+6)
        if exptime>np.max(exptime_grid):
            print("\nMaximum Exposure Time Exceeded (t>3600s).\n")
            return np.nan
        exptime_index=InterpolatedUnivariateSpline(logexp, np.arange(len(exptime_grid),
                                                dtype=np.double))(np.log10(exptime))
        snr_interpolator=RegularGridInterpolator((np.arange(len(exptime_grid)),
                              np.arange(len(vmag_grid)),
                              np.arange(len(teff_grid))),
                                    snr_grid)
        inputs=[exptime_index, vmag_index, teff_index]
        eta=snr_interpolator(inputs)[0]
        j+=1
    return exptime

def NEID_max_exptime(teff, vmag, exptime=60.):
    """
    calculate max recommended exposure time (60% full well) for a given target
    
    teff:          Effective Temperature (K)
    vmag:          V-band magnitude
    exptime:       Exposure time (s)
    """
    
    exptime_grid = fits.open(path_to_grid+'photon_grid_exptime.fits')[0].data
    teff_grid = fits.open(path_to_grid+'photon_grid_teff.fits')[0].data
    vmag_grid = fits.open(path_to_grid+'photon_grid_vmag.fits')[0].data
    logexp=np.log10(exptime_grid)
    
    bound_test=True
    if teff<np.min(teff_grid) or teff>np.max(teff_grid):
        print("Temperature out of bounds. The allowed range is %d K to %d K." % (np.min(teff_grid), np.max(teff_grid)))
        bound_test=False
    if vmag<np.min(vmag_grid) or vmag>np.max(vmag_grid):
        print("Magnitude out of bounds. The allowed range is V = %d to V = %d." % (np.min(vmag_grid), np.max(vmag_grid)))
        bound_test=False
    if exptime<np.min(exptime_grid) or exptime>np.max(exptime_grid):
        print("Exposure time out of bounds. The allowed range is %d s to %d s." % (np.min(exptime_grid), np.max(exptime_grid)))
        bound_test=False
    if bound_test==False:
        return np.nan
    
    wavelength_grid = fits.open(path_to_grid+'order_wvl_centers.fits')[0].data[1]
    snr_grid_order=fits.open(path_to_grid+'snr_master_order.fits')[0].data
    
    
    teff_index=InterpolatedUnivariateSpline(teff_grid, 
                            np.arange(len(teff_grid), dtype=np.double))(teff)

    vmag_index=InterpolatedUnivariateSpline(vmag_grid, 
                            np.arange(len(vmag_grid), dtype=np.double))(vmag)
    exptime_index=InterpolatedUnivariateSpline(logexp, np.arange(len(exptime_grid),
                                                dtype=np.double))(np.log10(exptime))
    snr=np.zeros(len(snr_grid_order))
    for o in range(len(snr_grid_order)):
        snr_grid=snr_grid_order[o]
        snr_interpolator=RegularGridInterpolator((np.arange(len(exptime_grid)),
                              np.arange(len(vmag_grid)),
                              np.arange(len(teff_grid))),
                                    snr_grid)
        inputs=[exptime_index, vmag_index, teff_index]
        snr[o]=snr_interpolator(inputs)[0]
    
    softlimit=np.array([497.82248, 500.46307, 499.61078, 495.75638, 495.59677, 497.53268,
       494.92285, 497.87125, 495.24042, 499.50885, 499.4121 , 495.15826,
       498.5816 , 498.41797, 494.5738 , 498.28055, 495.77396, 495.66986,
       497.76923, 495.65466, 495.5136 , 494.9543 , 497.7236 , 496.294  ,
       494.54578, 492.45773, 493.8629 , 495.3245 , 496.50537, 497.20462,
       496.3666 , 493.56155, 495.66016, 495.95648, 495.29364, 493.72418,
       495.75903, 496.19632, 494.97125, 492.46533, 496.4472 , 494.52365,
       493.13345, 494.53046, 495.96414, 494.08124, 491.92264, 494.0127 ,
       495.47772, 495.2882 , 495.43652, 494.86865, 495.158  , 495.0169 ,
       494.81766, 495.0668 , 495.3031 , 493.20224, 494.71045, 493.0483 ,
       495.01123, 494.97443, 492.57593, 493.84567, 490.0411 , 493.2728 ,
       491.09906, 494.50082, 493.13843, 494.31946, 494.4462 , 491.75977,
       493.00977, 491.4752 , 491.64948, 492.39023, 492.10922, 489.91223,
       494.3898 , 494.33267, 490.79178, 493.4196 , 491.94476, 494.1685 ,
       492.94263, 494.19852, 489.60358, 492.3478 , 490.21204, 493.80936,
       489.8438 , 490.6124 , 490.24118, 493.70288, 491.47467])
    
    peak_arg=np.argmax(snr/softlimit)
    wvl=wavelength_grid[peak_arg]
    snr_threshold=softlimit[peak_arg]
    
    max_exp=NEID_exptime_SNR(teff, vmag, snr_threshold, wvl)
    if np.isnan(max_exp):
        max_exp=3600.
    
    return max_exp
