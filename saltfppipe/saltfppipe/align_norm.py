from astropy.io.fits import open as openfits
import numpy as np
from saltfppipe.fit_sky_level import fit_sky_level
from sys import exit as crash
from os import remove
import matplotlib.pyplot as plt
from photutils import daofind, CircularAperture, aperture_photometry
from pyraf import iraf

class PlotFluxRatios:
    def __init__(self,goodfluxrats,fluxratavg,fluxratsig):
        #Initialize the coordinates and key presses
        self.xcoo = None
        self.ycoo = None
        self.key = None
        #Create the plot
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        self.plot = ax1.scatter(np.arange(len(goodfluxrats[:,0])),goodfluxrats[:,0],color="blue",marker="o")
        ax1.plot(np.arange(len(goodfluxrats[:,0])),goodfluxrats[:,0],color="blue")
        for i in range(1,goodfluxrats.shape[1]):
            ax1.plot(np.arange(len(goodfluxrats[:,i])),goodfluxrats[:,i],color="blue",marker="o",linestyle="-")
        ax1.errorbar(np.arange(len(fluxratavg)), fluxratavg, fluxratsig, color="red", marker="s",linestyle="--")
        plt.title("Press 'd' on a point to delete an obviously bad star from the fit.\n"+
                  "Press 'r' to restore all previously deleted stars.\n"+
                  "Press 'a' to accept the fit.")
        plt.xlabel("Image Number")
        plt.ylabel("Normalization Value")
        self.cid = self.plot.figure.canvas.mpl_connect("key_press_event",self.keypress)
        plt.tight_layout()
        plt.show()
        
    def keypress(self,event):
        if event.key in ["d","r","a"] and event.inaxes == self.plot.axes:
            self.key = event.key
            self.xcoo = event.xdata
            self.ycoo = event.ydata
            plt.close()
        return
    
def calc_norm(counts,dcounts):
    """Calculates the normalization values each image must be divided by
    in order to normalize them all to each other, as well as the uncertainty
    in these values.
    
    The routine is interactive, in that it plots the stellar photometry and
    the user can remove troublesome stars from the fit.
    
    It is very likely you don't want to call this routine on its own, but
    rather have align_norm do it for you.
    
    Inputs:
    counts -> Array of stellar fluxes
    dcounts -> Array of uncertainties in fluxes
    
    Outputs:
    norm -> List containing normalization values for each image
    dnorm -> Normalization uncertainty for each image
    
    """

    for i in range(counts.shape[1]):
        counts[:,i] *= 1/np.median(counts[:,i])
    
    #Interactively plot the star flux ratios
    goodcounts = counts.copy()
    while True:
        countavg = np.average(goodcounts,axis=1)
        countsig = np.std(goodcounts,axis=1)
        normplot = PlotFluxRatios(goodcounts,countavg,countsig)
        if normplot.key == "a": break
        if normplot.key == "r": goodcounts = counts.copy()
        if normplot.key == "d": 
            #Figure out which point was "nearest" the keypress and delete it
            xcoo = np.int(np.rint(normplot.xcoo))
            star = np.argmin(np.abs(goodcounts[xcoo,:]-normplot.ycoo))
            goodcounts = np.delete(goodcounts,star,1)
    
    return countavg, countsig
    
def align_norm(fnlist,tolerance=3,thresh=3.5):
    """Aligns a set of images to each other, as well as normalizing the images
    to the same average brightness.
    
    Both the alignment and normalization are accomplished through stellar
    photometry using the IRAF routine 'daophot'. The centroids of a handful
    of stars are found and used to run the IRAF routine 'imalign'. The
    instrumental magnitudes of the stars are used to determine by how much
    each image must be scaled for the photometry to match across images.
    
    The images are simply updated with their rescaled, shifted selves. This
    overwrites the previous images and adds the header keyword 'fpphot' to
    the images.
    
    A handful of temporary files are created during this process, which should
    all be deleted by the routine at the end. But if it is interrupted, they
    might not be.
    
    If the uncertainty images exist, this routine also shifts them by the same
    amounts as the intensity images, as well as updating the uncertainty values
    for both the new normalization and the uncertainties in normalizing the
    images.
    
    Inputs:
    fnlist -> List of strings, each the path to a fits image.
    tolerance -> How close two objects can be and still be considered the same
                 object. Default is 3 pixels.
    thresh -> Optional. Level above sky background variation to look for objs.
              Default is 3.5 (times SkySigma). Decrease if center positions
              aren't being found accurately. Increase for crowded fields to
              decrease computation time.
    
    """
    
    #Fit for the sky background level
    skyavg, skysig = fit_sky_level(fnlist)
    
    #Get image FWHMs
    fwhm = np.empty(len(fnlist))
    firstimage = openfits(fnlist[0])
    toggle = firstimage[0].header.get("fpfwhm")
    axcen = firstimage[0].header.get("fpaxcen")
    aycen = firstimage[0].header.get("fpaycen")
    arad = firstimage[0].header.get("fparad")
    firstimage.close()
    if axcen == None:
        print "Error! Images have not yet been aperture-masked! Do this first!"
        crash()
    if toggle == None:
        print "Warning: FWHMs have not been measured! Assuming 5 pixel FWHM for all images."
        for i in range(len(fnlist)): fwhm[i] = 5
    else:
        for i in range(len(fnlist)):
            image = openfits(fnlist[i])
            fwhm[i] = image[0].header["fpfwhm"]
            image.close()
    
    #Identify the stars in each image
    xlists = []
    ylists = []
    for i in range(len(fnlist)):
        xlists.append([])
        ylists.append([])
        image = openfits(fnlist[i])
        axcen = image[0].header["fpaxcen"]
        aycen = image[0].header["fpaycen"]
        arad = image[0].header["fparad"]
        sources = daofind(image[1].data-skyavg[i],
                          fwhm=fwhm[i],
                          threshold=thresh*skysig[i]).as_array()
        for j in range(len(sources)):
            if np.logical_and((sources[j][1]-axcen)**2 + (sources[j][2]-aycen)**2 > 50**2,(sources[j][1]-axcen)**2 + (sources[j][2]-aycen)**2 < (0.95*arad)**2):
                xlists[i].append(sources[j][1])
                ylists[i].append(sources[j][2])
        image.close()
    
    #Match objects between fields
    print "Matching objects between images..."
    xcoo = []
    ycoo = []
    for i in range(len(xlists[0])):
        #For each object in the first image
        accept = True
        for j in range(1,len(fnlist)):
            #For each other image
            dist2 = (np.array(xlists[j])-xlists[0][i])**2 + (np.array(ylists[j])-ylists[0][i])**2
            if (min(dist2) > tolerance**2):
                accept = False
                break
        if accept:
            #We found an object at that position in every image
            xcoo.append(xlists[0][i])
            ycoo.append(ylists[0][i])
    
    #Create coordinate arrays for the photometry and shifting
    x = np.zeros((len(fnlist),len(xcoo)))
    y = np.zeros_like(x)
    for i in range(len(xcoo)):
        #For every object found in the first image
        for j in range(len(fnlist)):
            #Find that object in every image
            dist2 = (np.array(xlists[j])-xcoo[i])**2 + (np.array(ylists[j])-ycoo[i])**2
            index = np.argmin(dist2)
            x[j,i] = xlists[j][index]
            y[j,i] = ylists[j][index]

    #Do aperture photometry on the matched objects
    counts = np.zeros_like(x)
    dcounts = np.zeros_like(x)
    for i in range(len(fnlist)):
        image = openfits(fnlist[i])
        apertures = CircularAperture((x[i],y[i]), r=2*fwhm[i])
        phot_table = aperture_photometry(image[1].data-skyavg[i],
                                         apertures,error=np.sqrt(image[2].data))
        counts[i] = phot_table["aperture_sum"]
        dcounts[i] = phot_table["aperture_sum_err"]
        image.close()
    
    #Calculate the normalizations
    norm, dnorm = calc_norm(counts,dcounts)
    
    #Normalize the images and adjust the variance plane
    for i in range(len(fnlist)):
        print "Normalizing image "+fnlist[i]
        image = openfits(fnlist[i],mode="update")
        image[1].data /= norm[i]
        image[2].data = image[2].data/norm[i]**2 + image[1].data**2*dnorm[i]**2/norm[i]**2
        image.close()
        
    #Calculate the shifts
    for i in range(x.shape[1]):
        x[:,i] = -(x[:,i] - x[0,i])
        y[:,i] = -(y[:,i] - y[0,i])
    xshifts = np.average(x,axis=1)
    yshifts = np.average(y,axis=1)
    
    #Shift the images and update headers
    iraf.images(_doprint=0)
    iraf.immatch(_doprint=0)
    for i in range(len(fnlist)):
        #This is incredibly stupid, but works. Iraf is annoying.
        print "Shifting image "+fnlist[i]
        iraf.geotran(input=fnlist[i]+"[1]",
                     output="temp1.fits",
                     geometry="linear",
                     database="",
                     transforms="",
                     xin="INDEF",
                     yin="INDEF",
                     xshift=xshifts[i],
                     yshift=yshifts[i],
                     xout="INDEF",
                     yout="INDEF",
                     xmag="INDEF",
                     ymag="INDEF",
                     xrotation="INDEF",
                     yrotation="INDEF",
                     xmin="INDEF",
                     xmax="INDEF",
                     ymin="INDEF",
                     ymax="INDEF",
                     xscale=1.0,
                     yscale=1.0,
                     ncols="INDEF",
                     nlines="INDEF",
                     xsample=1.0,
                     ysample=1.0,
                     interpolant="linear",
                     boundary="nearest",
                     constant=0.0,
                     fluxconserve="yes",
                     nxblock=512,
                     nyblock=512,
                     verbose="no")
        iraf.geotran(input=fnlist[i]+"[2]",
                     output="temp2.fits",
                     geometry="linear",
                     database="",
                     transforms="",
                     xin="INDEF",
                     yin="INDEF",
                     xshift=xshifts[i],
                     yshift=yshifts[i],
                     xout="INDEF",
                     yout="INDEF",
                     xmag="INDEF",
                     ymag="INDEF",
                     xrotation="INDEF",
                     yrotation="INDEF",
                     xmin="INDEF",
                     xmax="INDEF",
                     ymin="INDEF",
                     ymax="INDEF",
                     xscale=1.0,
                     yscale=1.0,
                     ncols="INDEF",
                     nlines="INDEF",
                     xsample=1.0,
                     ysample=1.0,
                     interpolant="linear",
                     boundary="nearest",
                     constant=0.0,
                     fluxconserve="yes",
                     nxblock=512,
                     nyblock=512,
                     verbose="no")
        iraf.geotran(input=fnlist[i]+"[3]",
                     output="temp3.fits",
                     geometry="linear",
                     database="",
                     transforms="",
                     xin="INDEF",
                     yin="INDEF",
                     xshift=xshifts[i],
                     yshift=yshifts[i],
                     xout="INDEF",
                     yout="INDEF",
                     xmag="INDEF",
                     ymag="INDEF",
                     xrotation="INDEF",
                     yrotation="INDEF",
                     xmin="INDEF",
                     xmax="INDEF",
                     ymin="INDEF",
                     ymax="INDEF",
                     xscale=1.0,
                     yscale=1.0,
                     ncols="INDEF",
                     nlines="INDEF",
                     xsample=1.0,
                     ysample=1.0,
                     interpolant="linear",
                     boundary="nearest",
                     constant=0.0,
                     fluxconserve="yes",
                     nxblock=512,
                     nyblock=512,
                     verbose="no")
        #Open the image in question
        image = openfits(fnlist[i],mode="update")
        #Update the data arrays with the appropriate shifted frames
        temp = openfits("temp1.fits")
        image[1].data = temp[0].data
        temp.close()
        temp = openfits("temp2.fits")
        image[2].data = temp[0].data
        temp.close()
        temp = openfits("temp3.fits")
        image[3].data = temp[0].data
        image[3].data[image[3].data>0]=1 #Bad pixel adjustment
        temp.close()
        #Update the image headers
        image[0].header["fpphot"]="True"
        image[0].header["fpxcen"]+=xshifts[i]
        image[0].header["fpycen"]+=yshifts[i]
        image[0].header["fpaxcen"]+=xshifts[i]
        image[0].header["fpaycen"]+=yshifts[i]
        #Close the image and delete temp files
        image.close()
        remove("temp1.fits")
        remove("temp2.fits")
        remove("temp3.fits")

    return