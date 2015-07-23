import numpy as np
from astropy.io.fits import open as openfits
from astropy.io.fits import writeto as writefits
#from saltfppipe.voigt import fit_ha_and_n2, fit_only_ha
from os.path import join
import voigtfit

def fit_velmap_ha_n2_mode(wavelist, intylist, outdir, uncertlist=None, clobber=False):
    """
    """
    #Create blank lists for filling with images
    wavefilelist = []
    intyfilelist = []
    if not (uncertlist is None):
        uncertfilelist = []
    #Open images and populate the lists
    for i in range(len(wavelist)):
        wavefilelist.append(openfits(wavelist[i]))
        intyfilelist.append(openfits(intylist[i]))
        if not (uncertlist is None):
            uncertfilelist.append(openfits(uncertlist[i]))
    #Create blank arrays for filling with data
    wavecube = np.zeros((intyfilelist[0][0].data.shape[0],intyfilelist[0][0].data.shape[1],len(wavelist)))
    intycube = np.zeros((intyfilelist[0][0].data.shape[0],intyfilelist[0][0].data.shape[1],len(wavelist)))
    if not (uncertlist is None):
        uncertcube = np.zeros((intyfilelist[0][0].data.shape[0],intyfilelist[0][0].data.shape[1],len(wavelist)))
    #Fill the cube arrays with data
    for i in range(len(wavelist)):
        wavecube[:,:,i] = wavefilelist[i][0].data
        intycube[:,:,i] = intyfilelist[i][0].data
        if not (uncertlist is None):
            uncertcube[:,:,i] = uncertfilelist[i][0].data
    #Create blank arrays for the velocity+etc fits
    contarray = np.zeros_like(wavefilelist[0][0].data)
    intyarray = np.zeros_like(wavefilelist[0][0].data)
    intyratioarray = np.zeros_like(wavefilelist[0][0].data)
    wavearray = np.zeros_like(wavefilelist[0][0].data)
    two_line_success_array = np.zeros_like(wavefilelist[0][0].data,dtype="bool")
    one_line_success_array = np.zeros_like(wavefilelist[0][0].data,dtype="bool")
    #Create x and y arrays for the "progress bar"
    numb_updates=20 #Number of progress updates
    xgrid, ygrid = np.meshgrid(np.arange(contarray.shape[1]),np.arange(contarray.shape[0]))
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
    for y in range(wavearray.shape[0]):
        for x in range(wavearray.shape[1]):
            #Display the "progress bar" if appropriate
            if np.any(np.logical_and(x==progX,y==progY)): print "Velocity map approximately "+repr(progPerc[np.logical_and(progX==x,progY==y)][0])+"% complete..."
            #Get the spectrum to fit
            waves = wavecube[y,x,:]
            intys = intycube[y,x,:]
            if not (uncertlist is None): uncerts = uncertcube[y,x,:]
            else: uncerts = None
            
#             #Attempt the two-lined fit (halpha and NII)
#             if np.sum(waves!=0)>0.5*len(waves): fit, success = fit_ha_and_n2(waves[waves!=0], intys[intys!=0], uncerts[uncerts!=0])
#             else: fit, success = np.zeros(8), False
#             #Was the fit successful
#             if success:
#                 contarray[y,x] = fit[0]
#                 intyarray[y,x] = fit[1]
#                 intyratioarray[y,x] = fit[2]/fit[1]
#                 wavearray[y,x] = fit[3]
#                 two_line_success_array[y,x] = True
#             #If not, let's try only fitting the halpha line
#             else:
#                 if np.sum(waves!=0)>0.5*len(waves): fit, success = fit_only_ha(waves[waves!=0], intys[intys!=0], uncerts[uncerts!=0])
#                 else: fit, success = np.zeros(5), False
#                 if success:
#                     contarray[y,x] = fit[0]
#                     intyarray[y,x] = fit[1]
#                     wavearray[y,x] = fit[2]
#                     one_line_success_array[y,x] = True
            #Fit with voigtfit
            
            if np.sum(waves!=0)>0.5*len(waves):
                fit,_err,_chi2 = voigtfit.voigtfit(waves[waves!=0],intys[intys!=0],uncerts[uncerts!=0],voigtfit.voigtguess(waves[waves!=0],intys[intys!=0]),[True,True,True,True,True])
                contarray[y,x] = fit[0]
                intyarray[y,x] = fit[1]
                wavearray[y,x] = fit[2]
    #Convert wavelengths to velocities
    velarray = np.zeros_like(wavearray)
    velarray[wavearray!=0] = 299792.458*((wavearray[wavearray!=0]/6562.81)**2-1)/((wavearray[wavearray!=0]/6562.81)**2+1)
    #Save the output images
    writefits(join(outdir,"velocity.fits"), velarray, clobber=clobber)
    writefits(join(outdir,"intensity.fits"), intyarray, clobber=clobber)
    writefits(join(outdir,"intyratio.fits"), intyratioarray, clobber=clobber)
    writefits(join(outdir,"continuum.fits"), contarray, clobber=clobber)
    #writefits(join(outdir,"2linesuccess.fits"), two_line_success_array, clobber=clobber)
    #writefits(join(outdir,"1linesuccess.fits"), one_line_success_array, clobber=clobber)
    #Close input images
    for i in range(len(wavelist)):
        wavefilelist[i].close()
        intyfilelist[i].close()
        if not (uncertlist is None):
            uncertfilelist[i].close()
    
    return