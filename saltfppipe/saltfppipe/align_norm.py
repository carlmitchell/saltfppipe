from astropy.io.fits import open as openfits
import numpy as np
from pyraf import iraf
from saltfppipe.identify_objects import identify_objects, clean_files
from saltfppipe.fit_sky_level import fit_sky_level
from sys import exit as crash
from os.path import isfile
from os import remove
import matplotlib.pyplot as plt

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
    
def match_objects(coolist,coofile="stars.coo",tolerance=3):
    """This routine opens a list of coordinate files of the type produced by
    IRAF's 'daofind' task and looks for objects that appear at roughly the same
    position across all images. Creates a coordinate file, by default it is
    'stars.coo'. It is very likely you don't want to call this routine on its
    own, but rather have align_norm do it for you.
    
    Inputs:
    coolist -> List of strings, each the path to a coordinate file
    coofile -> Path to the resulting coordinate file. Default 'stars.coo'
    tolerance -> Number of pixels that objects can differ by and still be
                 considered the same object. Default 3 pixels.
    
    Outputs:
    coofile -> The path to the resulting coordinate file.
    
    """
    
    print "Matching objects between images..."
    
    #Check to see if the coordinate file already exists
    if isfile(coofile):
        while True:
            yn = raw_input("Coordinate file "+coofile+" already exists. Delete it? (y/n) ")
            if "n" in yn or "N" in yn: return coofile
            elif "y" in yn or "Y" in yn:
                remove(coofile)
                break
    
    #Open each coordinate file and read in the X, Y coordinates
    xarray = []
    yarray = []
    for i in range(len(coolist)):
        xarray.append([])
        yarray.append([])
        infile = open(coolist[i])
        data = infile.readlines()
        for j in range(41,len(data)):
            line = data[j].split()
            xarray[i].append(float(line[0]))
            yarray[i].append(float(line[1]))
        infile.close()
    
    #Search for objects that are found in all images
    xcoo = []
    ycoo = []
    for i in range(len(xarray[0])):
        #For each object in the first image...
        accept = True
        for j in range(1,len(xarray)):
            #For each other image...
            dist2 = (np.array(xarray[j])-xarray[0][i])**2 + (np.array(yarray[j])-yarray[0][i])**2
            if min(dist2) > tolerance**2:
                #There's no object within "tolerance" pixels of that location
                accept = False
                break
        if accept:
            xcoo.append(xarray[0][i])
            ycoo.append(yarray[0][i])
    
    #Write the coordinates list to file
    print repr(len(xcoo))+" stars will be written to coordinate file "+coofile+"."
    outfile = open(coofile, 'w')
    for i in range(0,len(xcoo)):
        outfile.write(str(xcoo[i])+' '+str(ycoo[i])+'\n')
    outfile.close()
    
    return coofile
    
def do_phot(fnlist,coofile,fwhmlist,skysiglist):
    """Runs the IRAF 'phot' routine on a series of images to determine their
    centroid positions and magnitudes. It is very likely you don't want to call
    this routine on its own, but rather have align_norm do it for you.
    
    Inputs:
    fnlist -> List of strings, each the path to a fits image
    coofile -> Path to the file containing the coordinate list of stars
    fwhmlist -> List of the PSF FWHM's for each image
    skysiglist -> List of the sky background sigmas for each image
    
    Outputs:
    photlist -> List of strings, the paths to the photometry outputs
    
    """
        
    iraf.noao(_doprint=0)
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    photlist = []
    for i in range(len(fnlist)):
        print "Performing photometry on the stars in image "+fnlist[i]
        photlist.append(fnlist[i]+".mag")
        iraf.phot(image=fnlist[i],
                  skyfile="",
                  coords=coofile,
                  output=photlist[i],
                  interactive="N",
                  sigma=skysiglist[i],
                  datamin='INDEF',
                  datamax='INDEF',
                  calgorithm='centroid',
                  cbox=2*fwhmlist[i],
                  salgorithm='centroid',
                  annulus=4*fwhmlist[i],
                  dannulus=6*fwhmlist[i],
                  apertures=4*fwhmlist[i],
                  verify='N')
    
    return photlist
    
def read_phot(photlist):
    """This routine reads the photometry files produced by IRAF's 'phot' and
    picks out the (X,Y) coordinates and magnitudes of each star in each image,
    as well as the uncertainties in magnitudes, and returns them. It is very
    likely you don't want to call this routine on its own, but rather have
    align_norm do it for you.
    
    Inputs:
    photlist -> List of paths to photometry files.
    
    Outputs:
    x -> 2D array of x positions of stars
    y -> 2D array of y positions of stars
    mag -> 2D array of stellar magnitudes
    dmag -> 2D array of magnitude uncertainties
    
    """
    
    x = []
    y = []
    mag = []
    dmag = []
    flag = []
    for i in range(len(photlist)):
        #Open and read the photometry file
        infile = open(photlist[i])
        data = infile.readlines()
        for i in range(75,len(data),5):
            str2 = data[i+1].split()
            str5 = data[i+4].split()
            #Append the values for X, Y, Mag, dMag, and Flag to the lists
            x.append(float(str2[0]))
            y.append(float(str2[1]))
            if str5[4]!='INDEF' and str5[5]!='INDEF':
                mag.append(float(str5[4]))
                dmag.append(float(str5[5]))
                flag.append(True)
            else:
                mag.append(0.0)
                dmag.append(0.0)
                flag.append(False)
        infile.close()
    
    #Figure out the shape of the resultant arrays and reshape them
    newshape = (len(photlist),len(x)/len(photlist))
    x = np.reshape(x,newshape)
    y = np.reshape(y,newshape)
    mag = np.reshape(mag,newshape)
    dmag = np.reshape(dmag,newshape)
    flag = np.reshape(flag,newshape)
    
    #Eliminate any star "columns" that have flags of False.
    badstars = []
    for i in range(flag.shape[1]):
        if sum(flag[:,i])!=len(flag[:,i]):
            badstars.append(i)
    for i in badstars[::-1]:
        x = np.delete(x,i,1)
        y = np.delete(y,i,1)
        mag = np.delete(mag,i,1)
        dmag = np.delete(dmag,i,1)

    return x, y, mag, dmag
    
def calc_norm(mag,dmag):
    """Calculates the normalization values each image must be multiplied by
    in order to normalize them all to each other, as well as the uncertainty
    in these values. 
    
    The routine is interactive, in that it plots the stellar photometry and
    the user can remove troublesome stars from the fit.
    
    It is very likely you don't want to call this routine on its own, but
    rather have align_norm do it for you.
    
    Inputs:
    mag -> Array of stellar magnitudes
    dmag -> Array of uncertainties
    
    Outputs:
    norm -> List containing normalization values for each image
    dnorm -> Normalization uncertainty for each image
    
    """
    

    
    #Convert from magnitudes to flux ratios
    fluxrats = 10**(0.4*mag)
    for i in range(fluxrats.shape[1]):
        fluxrats[:,i] *= 1/np.average(fluxrats[:,i])
    
    #Interactively plot the star flux ratios
    goodfluxrats = fluxrats.copy()
    while True:
        fluxratavg = np.average(goodfluxrats,axis=1)
        fluxratsig = np.std(goodfluxrats,axis=1)
        normplot = PlotFluxRatios(goodfluxrats,fluxratavg,fluxratsig)
        if normplot.key == "a": break
        if normplot.key == "r": goodfluxrats = fluxrats.copy()
        if normplot.key == "d": 
            #Figure out which point was "nearest" the keypress and delete it
            xcoo = np.int(np.rint(normplot.xcoo))
            star = np.argmin(np.abs(goodfluxrats[xcoo,:]-normplot.ycoo))
            goodfluxrats = np.delete(goodfluxrats,star,1)
    
    return fluxratavg, fluxratsig
    
def align_norm(fnlist,uncertlist=None):
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
    uncertlist (optional) -> List of paths to uncertainty images.
    
    """
    
    #Fit for the sky background level
    _skyavg, skysig = fit_sky_level(fnlist)
    
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
    
    #Identify objects in the fields
    coolist = identify_objects(fnlist,skysig,fwhm)
    
    #Match objects between fields
    coofile = match_objects(coolist)
    
    #Do aperture photometry on the matched objects
    photlist = do_phot(fnlist,coofile,fwhm,skysig)
    
    #Read the photometry files
    x, y, mag, dmag = read_phot(photlist)
    
    #Calculate the normalizations
    norm, dnorm = calc_norm(mag,dmag)
    
    #Normalize the images (and optionally, the uncertainty images)
    for i in range(len(fnlist)):
        print "Normalizing image "+fnlist[i]
        image = openfits(fnlist[i],mode="update")
        if not (uncertlist is None):
            uncimage = openfits(uncertlist[i],mode="update")
            uncimage[0].data = np.sqrt(norm[i]**2*uncimage[0].data**2 + dnorm[i]**2*image[0].data**2)
            uncimage.close()
        image[0].data *= norm[i]
        image.close()
    
    #Calculate the shifts
    for i in range(x.shape[1]):
        x[:,i] = -(x[:,i] - x[0,i])
        y[:,i] = -(y[:,i] - y[0,i])
    xshifts = np.average(x,axis=1)
    yshifts = np.average(y,axis=1)
    
    #Shift the images (and optionally, the uncertainty images)
    iraf.images(_doprint=0)
    iraf.immatch(_doprint=0)
    for i in range(len(fnlist)):
        print "Shifting image "+fnlist[i]
        iraf.geotran(input=fnlist[i],
                     output=fnlist[i],
                     geometry="linear",
                     xshift=xshifts[i],
                     yshift=yshifts[i],
                     database="",
                     verbose="no")
        if not (uncertlist is None):
            iraf.geotran(input=uncertlist[i],
                         output=uncertlist[i],
                         geometry="linear",
                         xshift=xshifts[i],
                         yshift=yshifts[i],
                         database="",
                         verbose="no")
    
    #Update the image headers
    for i in range(len(fnlist)):
        image = openfits(fnlist[i],mode="update")
        image[0].header["fpphot"]="True"
        image[0].header["fpxcen"]+=xshifts[i]
        image[0].header["fpycen"]+=yshifts[i]
        image[0].header["fpaxcen"]+=xshifts[i]
        image[0].header["fpaycen"]+=yshifts[i]
        image.close()
    
    #Clean up the coordinate file list
    clean_files(fnlist)
    remove(coofile)
    for i in range(len(photlist)):
        remove(photlist[i])
    
    return