import numpy as np
from astropy.io.fits import open as openfits
from astropy.io.fits import writeto as writefits
#from saltfppipe.voigt import fit_ha_and_n2, fit_only_ha
from os.path import join
import voigtfit

def fit_velmap_ha_n2_mode(fnlist, outdir, clobber=False):
    """
    """
    #Create blank lists for filling with images
    imagelist = []
    #Open images and populate the lists
    for i in range(len(fnlist)):
        imagelist.append(openfits(fnlist[i]))
    #Create blank arrays for filling with data
    wavecube = np.zeros((imagelist[0][1].data.shape[0],imagelist[0][1].data.shape[1],len(fnlist)))
    intycube = np.zeros((imagelist[0][1].data.shape[0],imagelist[0][1].data.shape[1],len(fnlist)))
    uncertcube = np.zeros((imagelist[0][1].data.shape[0],imagelist[0][1].data.shape[1],len(fnlist)))
    badpixcube = np.zeros((imagelist[0][1].data.shape[0],imagelist[0][1].data.shape[1],len(fnlist)))
    #Fill the cube arrays with data
    for i in range(len(fnlist)):
        wavecube[:,:,i] = imagelist[i][4].data
        intycube[:,:,i] = imagelist[i][1].data
        uncertcube[:,:,i] = imagelist[i][2].data
        badpixcube[:,:,i] = imagelist[i][3].data
    #Create blank arrays for the velocity+etc fits
    contarray = np.zeros_like(imagelist[0][1].data)
    intyarray = np.zeros_like(imagelist[0][1].data)
    intyratioarray = np.zeros_like(imagelist[0][1].data)
    wavearray = np.zeros_like(imagelist[0][1].data)
    two_line_success_array = np.zeros_like(imagelist[0][1].data,dtype="bool")
    one_line_success_array = np.zeros_like(imagelist[0][1].data,dtype="bool")
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
            badpixs = badpixcube[y,x,:]
            waves = wavecube[y,x,:][badpixs!=1]
            intys = intycube[y,x,:][badpixs!=1]
            uncerts = uncertcube[y,x,:][badpixs!=1]
            
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
            if np.sum(badpixs)<0.5*len(waves):
                fit,_err,_chi2 = voigtfit.voigtfit(waves,intys,uncerts,voigtfit.voigtguess(waves,intys),[True,True,True,True,True])
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
    for i in range(len(imagelist)):
        imagelist[i].close()
    
    return