"""Generic mathematical functions used by plotting scripts and
EBM fitting.
"""

import numpy as np
from scipy import odr
from scipy.stats import linregress, t as students_t


def ordinary_least_squares(x, y):
    """Ordinary least-squares regression of y on x:
    
    y = b0 + b1*x + e
    
    where b0 and b1 are the fit parameters and e is a residual
    term (not returned) here. This is a wrapper function for
    scipy.stats.linregress.
    
    
    Parameters
    ----------
    x, y : 1D arrays
        Input data.
    
    
    Returns
    -------
    b0, b0_se : float
        First fit parameter (intercept) and standard error.
    
    b1, b1_se : float
        Second fit parameter (slope) and standard error.
    
    Note that standard error computation assumes the residuals
    are normally distributed and that there is no error noise
    in the x data.
    """
    
    ols_fit = linregress(x, y)
    b0      = ols_fit.intercept
    b0_se   = ols_fit.intercept_stderr
    b1      = ols_fit.slope
    b1_se   = ols_fit.stderr
    
    return b0, b0_se, b1, b1_se



def orthogonal_distance_regression(x, y):
    """Orthogonal distance regression fit for y against x. Uses
    equal weighting per observation and returns the fit
    parameters. This is a wrapper function for scipy.odr.
    
    
    Parameters
    ----------
    x, y : 1D arrays
        Input data.
    
    
    Returns
    -------
    b0, b0_se : float
        Intercept and its standard error.
    
    b1, b1_se : float
        Linear coefficient (slope) and its standard error.
    
    """
    
    data = odr.Data(x, y)
    odr_fit = odr.ODR(data, model=odr.unilinear).run()
    
    b1, b0 = odr_fit.beta
    b1_se, b0_se = odr_fit.sd_beta
    
    return b0, b0_se, b1, b1_se



def correlation_coefficient(x1, x2):
    """Correlation coefficient between two datasets x1 and x2,
    both arrays of the same size, computed using numpy.corrcoef.
    """
    return np.corrcoef(x1.flatten(), x2.flatten())[0,1]



def correlation_critical_value(n, conf_lev=0.95):
    """Critical value of the correlation coefficient between n
    pairs of values at the confidence level conf_lev [0, 1].
    """
    # Critical values of students t statistic with
    # n-2 degrees of freedom:
    t_crit = students_t.ppf(conf_lev, n-2)
    
    # Solve for correlation critical value (+/-) and return:
    return t_crit / np.sqrt(n-2 + t_crit**2)



def cross_correlation(xs, y):
    """Cross correlation function between inputs xs and y, where
    xs is contains multiple datasets.
    """
    return np.array([correlation_coefficient(xs[:,k], y)
                    for k in range(len(xs[0,:]))])
