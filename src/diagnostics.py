"""Intermediate data processing such as averaging and other
diagnostics.
"""

import numpy as np
npna = np.newaxis  # alias


# ============================================================ #
# GRID PROPERTIES
# ============================================================ #

def estimate_lon_lat_bnds_1D(lon, lat):
    """"""
    
    lon_bnds = np.zeros((len(lon), 2), dtype=np.float64)
    lat_bnds = np.zeros((len(lat), 2), dtype=np.float64)
    
    for x, x_bnds in zip([lon, lat], [lon_bnds, lat_bnds]):
        
        x_bnds[1:-1,:] = 0.5*np.array([(x[1:-1] + x[:-2]),
                                      (x[2:] + x[1:-1])]).T
    
        x_bnds[0,1] = 0.5*(x[1] + x[0])
        x_bnds[0,0] = 0.5*(3*x[0] - x[1])
        
        x_bnds[-1,0] = 0.5*(x[-1] + x[-2])
        x_bnds[-1,1] = 0.5*(3*x[-1] - x[-2])
    
    return lon_bnds, lat_bnds


# ============================================================ #
# TIME AVERAGING
# ============================================================ #

def year_mean_1D(field, nsteps_year=12, keep_nan=True):
    """Calculate yearly mean of a 1-spatial dimension field
    (e.g., hfbasin) from monthly data by default.
    """
    nt, ny = np.shape(field)
    field_rs = np.reshape(field,
        (nt//nsteps_year, nsteps_year, ny)
    )
    
    if keep_nan:
        return np.mean(field_rs, axis=1)
    else:
        return np.nanmean(field_rs, axis=1)



def year_mean_time_series_multi_member(field, nsteps_year=12,
         keep_nan=True):
    """Calculate yearly mean of a field (nt, n_ens)."""
    return year_mean_1D(field, nsteps_year, keep_nan)



def year_mean_2D(field, nsteps_year=12, keep_nan=True):
    """Calculate yearly mean of a 2D spatial field, from monthly
    data by default.
    """
    nt, ny, nx = np.shape(field)
    field_rs = np.reshape(field,
        (nt//nsteps_year, nsteps_year, ny, nx)
    )
    
    if keep_nan:
        return np.mean(field_rs, axis=1)
    else:
        return np.nanmean(field_rs, axis=1)



# ============================================================ #
# SPATIAL AVERAGING/INTEGRATION
# ============================================================ #

def integrate_horizontally_multi_members(field, areacell, lat,
        ref_lats=[65.0], hemi="n", land_mask=None, verbose=True
    ):
    """Integrate a 2D field [4D array (time, member, y, x)] over
    both spatial dimensions (y, x) with specified cell areas,
    everywhere poleward of a reference latitude. Returns both 
    the integral (no normalisation) and mean (with).
    """
    
    if hemi.lower() not in ["n", "s"]:
        raise Exception("Parameter \'hemi\' must be \'n\' or \'"
                        + "s\' (got \'" + hemi + "\')")
    hemi = hemi.lower()
    
    if land_mask is None:
        # Determine it from the data itself, assuming NaN
        # corresponds to land:
        land_mask = np.where(
            np.all(np.isnan(field), axis=(0, 1)), 0.0, 1.0)
    
    nt, n_ens, ny, nx = np.shape(field)
    n_ref_lats = len(ref_lats)
    
    field_hint = np.zeros((nt, n_ens, n_ref_lats))
    field_mean = np.zeros((nt, n_ens, n_ref_lats))
    
    if verbose:
        print("Calculating area integrals (0%)", end="\r")
    
    for j in range(n_ref_lats):
        
        if hemi == "n":
            lat_msk = np.where(lat >= ref_lats[j], 1., 0.)
        else:
            lat_msk = np.where(lat <= ref_lats[j], 1., 0.)
        
        field_hint[:,:,j] = np.nansum(
            field*areacell[npna,npna,:,:]
                *lat_msk[npna,npna,:,:],
            axis=(2,3)
        )
        
        # Assume land mask is any grid cell that is missing
        # (nan, assumed) for all time and all ensemble members:
        field_mean[:,:,j] = field_hint[:,:,j] / np.nansum(
            areacell*lat_msk*land_mask
        )
        
        if verbose:
            print("Calculating area integrals ("
                + f"{100*(j+1)/n_ref_lats:.0f}%)",
                end="\r"
            )
    
    if verbose:
        print("")
    
    return field_hint, field_mean



def integrate_horizontally_exact(field, lat_bnds, areacell,
        hemi="n", normalise=True, verbose=True
    ):
    """Calculate horizontal integrals poleward of grid cell
    latitude bounds on regular, fixed grids (i.e., independent
    longitude and latitude axes). Assumes cells are
    contiguous.
    
    
    Parameters
    ----------
    field : array (nt, n_ens, n_lat, n_lon)
        Field to integrate for each time (axis=0) and
        ensemble member (axis=1).
    
    lat_bnds : array (n_lat, 2)
        Cell latitude bounds such that lat_bnds[i,0]
        is the southward latitude and lat_bnds[i,1]
        is the northward latitude of cell i.
    
    areacell : array (n_lat, n_lon)
        Cell areas in m2.
    
    """
    
    nt, n_ens, n_lat, _ = np.shape(field)
    
    field_hint = np.zeros((nt, n_ens, n_lat))
    
    if verbose:
        print("Calculating area integrals (0%)", end="\r")
    
    if hemi.lower() == "n":
        
        # All original latitudes from cell edges
        # (assumes cells are contiguous):
        lat_org = lat_bnds[:,0]
        lat_org_bnds = np.array([
            [lat_org[k], lat_bnds[-1,1]]
            for k in range(len(lat_org))
        ])
        
        # Integrate northward of latitudes
        # 
        # Note that lat_org[j] = lat_bnds[j,0], the southward
        # edge of cell j, so that we want to sum all from j
        # onwards (northwards)
        
        for j in range(n_lat):
            
            field_hint[:,:,j] = np.nansum(
                field[:,:,j:,:]*areacell[npna,npna,j:,:],
                axis=(2,3)
            )
            
            if normalise:
                field_hint[:,:,j] /= np.nansum(areacell[j:,:])
                
            if verbose:
                print("Calculating area integrals ("
                    + f"{100*(j+1)/n_lat:.0f}%)",
                    end="\r"
                )
    
    else:
        
        # All original latitudes from cell edges
        # (assumes cells are contiguous):
        lat_org = lat_bnds[:,1]
        lat_org_bnds = np.array([
            [lat_bnds[0,0], lat_org[k]]
            for k in range(len(lat_org))
        ])
        
        # Integrate southward of latitudes
        # 
        # Note that lat_org[j] = lat_bnds[j,1], the northward
        # edge of cell j, so that we want to sum all cells up
        # to and including j
        
        for j in range(n_lat):
            field_hint[:,:,j] = np.nansum(
                field[:,:,:j+1]*areacell[npna,npna,:j+1,:],
                axis=(2,3)
            )
            
            if normalise:
                field_hint[:,:,j] /=np.nansum(areacell[:j+1,:])
                
            if verbose:
                print("Calculating area integrals ("
                    + f"{100*(j+1)/(n_lat):.0f}%)",
                    end="\r"
                )
    
    if verbose:
        print("Calculating area integrals (100%)")
    
    return lat_org, lat_org_bnds, field_hint



def interpolate_to_ref_latitudes(true_lats, true_data,
        ref_lats):
    """Interpolate data to specified reference latitudes
    (usually called after integrate_horizontally_exact).
    
    Uses linear interpolation of the original points
    k, k+1 surrounding each specified reference latitude.
    
    
    Parameters
    ----------
    true_lats : array (n_true_lats,)
        Latitudes of original data.
    
    field : array (nt, n_ens, n_true_lats)
        Field to interpolate to for each time (axis=0)
        and ensemble member (axis=1).
    
    ref_lats : array (n_ref_lats,)
        Latitudes at which to interpolate field to.
    
    
    Returns
    -------
    interp_data : array (nt, n_ens, n_ref_lats)
        Interpolated data.
    
    """
    
    nt, n_ens, n_true_lats = np.shape(true_data)
    n_ref_lats  = len(ref_lats)
    
    interp_data = np.zeros((nt, n_ens, n_ref_lats))
    
    for j in range(n_ref_lats):
        
        weights = np.zeros(n_true_lats)
        
        if all(ref_lats[j] < true_lats):
            # Determine weights assuming persistence:
            weights[0] = 1
        elif all(ref_lats[j] > true_lats):
            # Determine weights assuming persistence:
            weights[-1] = 1
        else:
            # Find the first index k of true_lats
            # such that true_lats[k] >= ref_lats[j].
            # This index and the previous (k-1) bound
            # ref_lats[j]
            k = np.nonzero(true_lats >= ref_lats[j])[0][0]
            
            # Determine the weights assuming a
            # linear interpolation scheme:
            weights[k] = (ref_lats[j] - true_lats[k-1]
                         )/(true_lats[k] - true_lats[k-1])
            
            if k > 0:
                weights[k-1] = 1 - weights[k]
        
        interp_data[:,:,j] = np.nansum(
            true_data*weights[npna,npna,:],
            axis=2
        )
    
    return interp_data
