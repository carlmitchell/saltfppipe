import numpy as np
from scipy.optimize import curve_fit

def nGauss(x,sig):
    """Normalized gaussian centered at zero, evaluated at x, with st.dev sig"""
    return np.exp(-(x**2)/(2*sig**2))/(sig*np.sqrt(2*np.pi))

def gaussfunc(x,cont,inty,x0,sig):
    return cont+inty*nGauss((x-x0),sig)

def make_sum_gaussians(wavelist):
    """This function generates a function which is the sum of N gaussians
    centered at N wavelengths. The resulting function's arguments are the
    X coordinate, the continuum strength, the N intensities, and the N sigmas.
    
    Inputs:
    wavelist -> List of wavelengths for the gaussian centers
    
    Outputs:
    sum_gaussians -> A function with arguments described above
    
    In other words:
    We want the func. Give us the func.
    
    """
    
    def sum_gaussians(x,cont,*lineparams):
        res = cont
        n=len(wavelist)
        for i in range(n):
            res += lineparams[i]*nGauss(x-wavelist[i], lineparams[n+i])
        return res
    
    return sum_gaussians

def GaussFit(x, y):
    """Attempts to fit a Gaussian profile to a set of data points (x,y).
    Returns the fitted profile.
    
    """
    
    mask = np.logical_not(np.logical_or(np.isnan(y),np.isnan(x)))
    x=x[mask]
    y=y[mask]
    
    guess = np.empty(4)
    guess[0] = np.min(y)
    guess[1] = np.sum(((y[1:]+y[:-1])/2-np.min(y))*np.diff(x))
    guess[2] = x[y==np.max(y)][0]
    guess[3] = (np.max(x)-np.min(x))/10
    
    try:
        fit = curve_fit(gaussfunc,x,y,p0=guess)[0]
    except RuntimeError:
        print "Warning: GaussFit did not converge."
        fit = guess
    except TypeError:
        #print "Warning: GaussFit did not converge."
        fit = guess
    
    return fit