import numpy as np
from astropy.io import fits
from os.path import join
import voigtfit
from saltfppipe.fp_image_class import FPImage

def fit_velmap_ha_n2_mode(fnlist, outdir, clobber=False):
    """
    """
    
    #Parameters
    c = 299792.458 #Speed of light (km/s)
    ha = 6562.81 #H-Alpha wavelength (ang)
    
    #Open images
    images = [FPImage(fn) for fn in fnlist]
    
    #Create blank arrays for filling with data
    wavecube = np.zeros((images[0].ysize,images[0].xsize,len(fnlist)))
    intycube = np.zeros((images[0].ysize,images[0].xsize,len(fnlist)))
    uncertcube = np.zeros((images[0].ysize,images[0].xsize,len(fnlist)))
    badpixcube = np.zeros((images[0].ysize,images[0].xsize,len(fnlist)))
    
    #Fill the cube arrays with data
    for i in range(len(fnlist)):
        wavecube[:,:,i] = images[i].wave
        intycube[:,:,i] = images[i].inty
        uncertcube[:,:,i] = np.sqrt(images[i].vari)
        badpixcube[:,:,i] = images[i].badp

    #Create blank arrays for the velocity+etc fits
    contarray = np.zeros_like(images[0].inty)
    intyarray = np.zeros_like(images[0].inty)
    wavearray = np.zeros_like(images[0].inty)
    dcontarray = np.zeros_like(images[0].inty)
    dintyarray = np.zeros_like(images[0].inty)
    dwavearray = np.zeros_like(images[0].inty)
    
    #Create x and y arrays for the "progress bar"
    numb_updates=20 #Number of progress updates
    xgrid, ygrid = images[0].xyarrays()
    xgrid = xgrid.flatten()
    ygrid = ygrid.flatten()
    progX = []
    progY = []
    progPerc = []
    for i in range(1,numb_updates+1):
        progX.append(xgrid[np.int(len(xgrid)*i/numb_updates)-1])
        progY.append(ygrid[np.int(len(xgrid)*i/numb_updates)-1])
        progPerc.append(i*100/numb_updates)
    progX = np.array(progX)
    progY = np.array(progY)
    progPerc = np.array(progPerc)
    #Loop over pixels
    for y in range(images[0].ysize):
        for x in range(images[0].xsize):
            #Display the "progress bar" if appropriate
            if np.any(np.logical_and(x==progX,y==progY)):
                print ("Velocity map approximately "+
                       repr(progPerc[np.logical_and(progX==x,progY==y)][0])+
                       "% complete...")
            #Get the spectrum to fit
            badpixs = badpixcube[y,x,:]
            waves = wavecube[y,x,:][badpixs!=1]
            intys = intycube[y,x,:][badpixs!=1]
            uncerts = uncertcube[y,x,:][badpixs!=1]
            
            #Attempt the two-lined fit (halpha and NII)
            #CURRENTLY NOT BEING DONE

            #Fit with voigtfit
            if np.sum(badpixs)<0.5*len(waves):
                fit,err,_chi2 = voigtfit.voigtfit(waves,intys,uncerts,
                                                  voigtfit.voigtguess(waves,intys),
                                                  [True,True,True,True,True])
                contarray[y,x] = fit[0]
                intyarray[y,x] = fit[1]
                wavearray[y,x] = fit[2]
                dcontarray[y,x] = err[0]
                dintyarray[y,x] = err[1]
                dwavearray[y,x] = err[2]
                
    #Convert wavelengths to velocities
    velarray = np.zeros_like(wavearray)
    dvelarray = np.zeros_like(wavearray)
    mask = wavearray != 0
    velarray[mask] = c*( ((wavearray[mask]/ha)**2-1) /
                         ((wavearray[mask]/ha)**2+1) )
    dvelarray[mask] = c*dwavearray[mask]/ha
    
    #Save the output images
    fits.writeto(join(outdir,"velocity.fits"), velarray, clobber=clobber)
    fits.writeto(join(outdir,"intensity.fits"), intyarray, clobber=clobber)
    fits.writeto(join(outdir,"continuum.fits"), contarray, clobber=clobber)
    fits.writeto(join(outdir,"dvelocity.fits"), dvelarray, clobber=clobber)
    fits.writeto(join(outdir,"dcontinuum.fits"), dcontarray, clobber=clobber)
    fits.writeto(join(outdir,"dintensity.fits"), dintyarray, clobber=clobber)

    #Close input images
    for image in images: image.close()
    
    return