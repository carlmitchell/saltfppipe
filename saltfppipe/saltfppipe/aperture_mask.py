from astropy.io.fits import open as openfits
import numpy as np

def aperture_mask(fnlist, axcen, aycen, arad):
    """Opens a series of images and masks pixels outside a circular aperture.
    The images are overwritten with the new masked image. The values of axcen,
    aycen, and arad are added to the FITS headers as keywords 'fpaxcen',
    'fpaycen', and 'fparad' respectively.
    
    Inputs:
    fnlist -> A list containing strings, each the path to a fits image.
    axcen -> The aperture x center
    aycen -> The aperture y center
    arad -> The aperture radius
    
    """
    
    for i in range(len(fnlist)):
        print "Masking pixels outside aperture for image "+str(i+1)+" of "+str(len(fnlist))+": "+fnlist[i]
        image = openfits(fnlist[i],mode="update")
        xgrid, ygrid = np.meshgrid(np.arange(image[1].data.shape[1]),np.arange(image[1].data.shape[0]))
        image[1].data[np.power((xgrid-axcen),2) + np.power((ygrid-aycen),2) > np.power(arad,2)] = 0
        image[3].data[np.power((xgrid-axcen),2) + np.power((ygrid-aycen),2) > np.power(arad,2)] = 1
        image[0].header["fpaxcen"] = axcen
        image[0].header["fpaycen"] = aycen
        image[0].header["fparad"] = arad
        image.close()