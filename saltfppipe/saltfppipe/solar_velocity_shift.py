from pyraf import iraf
from astropy.io.fits import open as openfits
import sys
from cStringIO import StringIO
import numpy as np

def solar_velocity_shift(waveimagelist, restwave):
    """Opens the wavelength images of a data cube and shifts the wavelengths
    so they are shifted to the solar frame. This is accomplished using the
    IRAF task 'rvcorrect'.
    
    It does it in a really really stupid way, which is to temporarily redirect
    IRAF's output to a string, because IRAF doesn't actually save the output
    anywhere (that I know of), it only prints it to the terminal. Then this
    routine searches through that string for the necessary value. It's a bit
    backwards, but it works!
    
    Note that the boost here is a Lorentz boost, not a Galilean boost, though
    values for either will be similar (i.e. the relativistic formula is used)
    
    Inputs:
    waveimagelist -> List of strings, the paths to *wavelength* images to be
                     updated.
    restwave -> The rest wavelength of the line you're interested in.
    
    """
    
    c = 299792.458 #in km/s
    
    for i in range(len(waveimagelist)):
        sys.stdout = iraf_output = StringIO()
        image = openfits(waveimagelist[i],mode="update")
        iraf.rvcorrect(header="N",
                       input="N",
                       imupdate="N",
                       observatory="SAAO",
                       year=float(image[0].header["date-obs"].split("-")[0]),
                       month=float(image[0].header["date-obs"].split("-")[1]),
                       day=float(image[0].header["date-obs"].split("-")[2]),
                       ut=image[0].header["utc-obs"],
                       ra=image[0].header["ra"],
                       dec=image[0].header["dec"],
                       vobs=0,
                       mode="a")
        image[0].header["fpsolar"] = float(iraf_output.getvalue().split()[-6])
        iraf_output.close()
        beta_earth = ((image[0].data/restwave)**2-1) / ((image[0].data/restwave)**2+1)
        beta_shift = image[0].header["fpsolar"]/c
        beta_helio = (beta_earth+beta_shift)/(1+beta_earth*beta_shift)
        image[0].data = restwave*(1+beta_helio)/np.sqrt(1-beta_helio**2)
        image[0].data[np.isnan(image[0].data)] = 0
        image.close()
        
    sys.stdout = sys.__stdout__
    
    