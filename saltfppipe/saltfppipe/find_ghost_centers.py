from astropy.io.fits import open as openfits
from sys import exit as crash
import numpy as np
from saltfppipe.fit_sky_level import fit_sky_level
from saltfppipe.identify_objects import identify_objects, clean_files
import matplotlib.pyplot as plt

class Verify_Center_Plot:
    def __init__(self,fn,data,flag,xcen,ycen,objx,objy,ghox,ghoy):
        self.key = None
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        self.plot = ax1.imshow(data,
                               cmap = "Greys",
                               aspect = "equal",
                               vmin = np.percentile(data,1),
                               vmax = np.percentile(data,99),
                               origin = "lower",
                               interpolation = "nearest")
        plt.title("Image "+fn+", Flagged = "+str(flag)+"\n"+
                  "Press 'q' to quit.\n"+
                  "Use Right/Left Arrow Keys to Navigate \n"+
                  "Up/Down Arrow Keys to Flag Bad Center")
        for i in range(len(objx)):
            plt.plot([objx[i],ghox[i]],[objy[i],ghoy[i]],color="blue")
        plt.scatter([xcen],[ycen],color="red",s=30)
        plt.xlim(0,data.shape[1])
        plt.ylim(0,data.shape[0])
        self.cid = self.plot.figure.canvas.mpl_connect("key_press_event",self.keypress)
        plt.tight_layout()
        plt.show()
        
    def keypress(self,event):
        if event.key in ["left","right","up","down","q"]:
            self.key = event.key
            plt.close()

def likely_center(coolist,axcen,aycen,arad,tolerance=3):
    """Finds the most likely ghost center for a series of images. A 'likely'
    center is found for each image individually, then the most common of these
    is chosen as the 'most likely' center.
    
    Objects too close to the image center or too close to the image edge are
    rejected from the analysis because they have been found to be troublesome
    before.
    
    It's probable you don't want to call this function on its own, but rather
    let 'find_ghost_centers' call it, as that function creates the coordinate
    files for you and loads the aperture information for you.
    
    Inputs:
    coolist -> List of strings, each the path to a coordinate file
    axcen -> Aperture X center
    aycen -> Aperture Y center
    arad -> Aperture radius
    tolerance -> Number of pixels away from each other two centers can be and
                be considered the same point. Default is 3 pixels.
    
    Outputs:
    bestxcen -> Best x center
    bestycen -> Best y center
    
    """
    
    print "Finding likely optical centers for each image..."
    
    xcen = np.zeros(len(coolist))
    ycen = np.zeros(len(coolist))
    goodcen = np.zeros(len(coolist))
    
    #Search for the best fitting center in each image
    for i in range(len(coolist)):
        #Read in the coordinate file
        xlist = []
        ylist = []
        infile = open(coolist[i])
        data = infile.readlines()
        for j in range(41,len(data)):
            line = data[j].split()
            x = float(line[0])
            y = float(line[1])
            #Reject any objects too close to the center or too close to the edge
            if np.logical_and((x-axcen)**2 + (y-aycen)**2 > 100**2,(x-axcen)**2 + (y-aycen)**2 < (0.95*arad)**2):
                xlist.append(x)
                ylist.append(y)
        infile.close()
    
        #Create Xcen, Ycen arrays for every possible combination of 2 objects
        N = len(xlist)
        xcenarray = np.zeros(((N*N-N)/2))
        ycenarray = np.zeros(((N*N-N)/2))
        dist2array = np.zeros(((N*N-N)/2))
        sub1array = np.zeros(((N*N-N)/2))
        index = 0
        for j in range(N):
            for k in range(j+1,N):
                xcenarray[index] = xlist[j]+xlist[k]
                ycenarray[index] = ylist[j]+ylist[k]
                index = index+1 #THIS IS TERRIBLE CODE, I'M SO SORRY

        xcenarray, ycenarray = 0.5*xcenarray, 0.5*ycenarray
    
        #Cross check the various possible centers against each other
        for j in range(len(xcenarray)):
            dist2array = (xcenarray-xcenarray[j])**2 + (ycenarray-ycenarray[j])**2
            sub1array[j] = np.sum(dist2array<tolerance**2)
        #Determine the locations of the "best" centers.
        bestcenloc = np.where(sub1array==max(sub1array))[0]
        
        #Now cross check JUST the best ones against each other
        sub1array = np.zeros(len(bestcenloc))
        xcenarray = xcenarray[bestcenloc]
        ycenarray = ycenarray[bestcenloc]
        for j in range(len(bestcenloc)):
            dist2array = (xcenarray-xcenarray[j])**2 + (ycenarray-ycenarray[j])**2
            sub1array[j] = np.sum(dist2array<tolerance**2)
        #Again, determine the locations of the "best" centers.
        bestcenloc = np.where(sub1array==max(sub1array))[0]
        xcen[i] = np.average(xcenarray[bestcenloc])
        ycen[i] = np.average(ycenarray[bestcenloc])
        
        print "Likely center for file "+coolist[i]+": "+str(np.int(xcen[i]))+', '+str(np.int(ycen[i]))
    
    #Cross-check the various image's centers against each other
    for i in range(len(coolist)):
        dist2array = (xcen-xcen[i])**2 + (ycen-ycen[i])**2
        goodcen[i] = np.sum(dist2array<tolerance**2)
    
    #Determine where in the arrays the best fitting centers are
    bestcenloc = np.where(goodcen==max(goodcen))[0]
    return np.average(xcen[bestcenloc]), np.average(ycen[bestcenloc])

def fit_better_centers(xcen, ycen, coolist, axcen, aycen, arad, tolerance=3):
    """This routine improves the center measurements of 'likely_center' for each
    individual image of a track. It's very likely you never want to call this
    routine on its own; rather you should probably let 'find_ghost_centers' call
    it for you.
    
    Inputs:
    xcen -> Initial 'guess' at the x center location (must be very good guess)
    ycen -> Initial 'guess' at the y center location (must be very good guess)
    coolist -> List of strings, each the path to a coordinate file
    tolerance -> Maximum distance two points can be away from each other and be
                 considered the same point. Default is 3 pixels.
    
    Outputs:
    xcenlist -> List of x centers for each image
    ycenlist -> List of y centers for each image
    
    """

    print "Fitting better optical centers for each image..."
    
    #Initialize the arrays to be filled
    xcenlist = np.zeros(len(coolist))
    ycenlist = np.zeros(len(coolist))
    
    #For each image's coordinate file
    for i in range(len(coolist)):
        #Load the coordinate file
        xlist = []
        ylist = []
        maglist = []
        infile = open(coolist[i])
        data = infile.readlines()
        for j in range(41,len(data)):
            line = data[j].split()
            if np.logical_and((float(line[0])-axcen)**2 + (float(line[1])-aycen)**2 > 100**2,(float(line[0])-axcen)**2 + (float(line[1])-aycen)**2 < (0.95*arad)**2):
                xlist.append(float(line[0]))
                ylist.append(float(line[1]))
                maglist.append(float(line[2]))
        infile.close()
        
        #For each object, calculate the position of its reflected partner
        refxlist = 2.*xcen - np.array(xlist)
        refylist = 2.*ycen - np.array(ylist)
        
        objx = []
        objy = []
        ghox = []
        ghoy = []
        
        #Now search for objects near these
        for j in range(len(xlist)):
            dist2list = (xlist - refxlist[j])**2 + (ylist - refylist[j])**2
            matchlist = dist2list < tolerance**2
            if np.sum(matchlist) >= 1:
                #We found a match! Now we need to know where the match(es) is(are)
                matchloc = np.where(matchlist==1)[0][0]
                #Cool, now, is this thing brighter than the ghost?
                if maglist[matchloc] < maglist[j]:
                    #It is! Record it:
                    objx.append(xlist[matchloc])
                    objy.append(ylist[matchloc])
                    ghox.append(xlist[j])
                    ghoy.append(ylist[j])
                    
        #Center position
        if len(objx) == 0:
            xcenlist[i] = 0
            ycenlist[i] = 0
        else:
            xcenlist[i] = 0.5*(np.average(objx)+np.average(ghox))
            ycenlist[i] = 0.5*(np.average(objy)+np.average(ghoy))
            
    xcenlist[xcenlist==0] = np.average(xcenlist[xcenlist!=0])
    ycenlist[ycenlist==0] = np.average(ycenlist[ycenlist!=0])
    
    xcenlist[np.isnan(xcenlist)] = xcen
    ycenlist[np.isnan(ycenlist)] = ycen
    
    return xcenlist, ycenlist

def verify_center(fnlist, coolist, tolerance=3):
    """Displays images and their likely ghost pairs and allows the user to
    call a crash if the center locations aren't good.
    
    Inputs:
    fnlist -> List of strings, each the path to a fits image
    coolist -> List of strings, the paths to coordinate files for each image
    tolerance -> Maximum distance two points can be away from each other and
                 still be considered the same point. Default is 3 pixels.
    
    Outputs:
    goodcenters -> Boolean. True if centers are approved. False if not.
    
    """

    #Make a ton of empty lists
    xlist = []
    ylist = []
    maglist = []
    objx = []
    objy = []
    ghox = []
    ghoy = []    
    imagelist = []
    xcenlist = []
    ycenlist = []
    datalist = []
    axcenlist = []
    aycenlist = []
    aradlist = []
    flagged = np.zeros(len(fnlist),dtype="bool")
    
    #Open each image, make lists of the objects and ghosts in each image
    for i in range(len(fnlist)):
        imagelist.append(openfits(fnlist[i]))
        axcenlist.append(imagelist[i][0].header["fpaxcen"])
        aycenlist.append(imagelist[i][0].header["fpaycen"])
        aradlist.append(imagelist[i][0].header["fparad"])
        datalist.append(imagelist[i][0].data[aycenlist[i]-aradlist[i]:aycenlist[i]+aradlist[i],axcenlist[i]-aradlist[i]:axcenlist[i]+aradlist[i]])
        xlist.append([])
        ylist.append([])
        maglist.append([])
        objx.append([])
        objy.append([])
        ghox.append([])
        ghoy.append([])
        infile = open(coolist[i])
        data = infile.readlines()
        for j in range(41,len(data)):
            line = data[j].split()
            if np.logical_and((float(line[0])-axcenlist[i])**2 + (float(line[1])-aycenlist[i])**2 > 100**2,(float(line[0])-axcenlist[i])**2 + (float(line[1])-aycenlist[i])**2 < (0.95*aradlist[i])**2):
                xlist[i].append(float(line[0]))
                ylist[i].append(float(line[1]))
                maglist[i].append(float(line[2]))
        infile.close()
        xcenlist.append(imagelist[i][0].header["fpxcen"])
        ycenlist.append(imagelist[i][0].header["fpycen"])
        
        #Where reflected partners *should* be
        refxlist = 2.*xcenlist[i] - np.array(xlist[i])
        refylist = 2.*ycenlist[i] - np.array(ylist[i])
        
        #Now search for objects near these
        for j in range(len(xlist[i])):
            dist2list = (np.array(xlist[i]) - refxlist[j])**2 + (np.array(ylist[i]) - refylist[j])**2
            matchlist = dist2list < tolerance**2
            if np.sum(matchlist) >= 1:
                #We found a match! Now we need to know where the match is
                matchloc = np.where(matchlist==1)[0][0]
                #Is this thing brighter than the ghost?
                if maglist[i][matchloc] < maglist[i][j]:
                    #It is! Record it:
                    objx[i].append(xlist[i][matchloc])
                    objy[i].append(ylist[i][matchloc])
                    ghox[i].append(xlist[i][j])
                    ghoy[i].append(ylist[i][j])
        
    #Adjust object and ghost coordinates because of the aperture shift
    for i in range(len(fnlist)):
        objx[i] = np.array(objx[i])-axcenlist[i]+aradlist[i]
        objy[i] = np.array(objy[i])-aycenlist[i]+aradlist[i]
        ghox[i] = np.array(ghox[i])-axcenlist[i]+aradlist[i]
        ghoy[i] = np.array(ghoy[i])-aycenlist[i]+aradlist[i]
        xcenlist[i] = xcenlist[i]-axcenlist[i]+aradlist[i]
        ycenlist[i] = ycenlist[i]-aycenlist[i]+aradlist[i]
            
    #Interactive plotting routine thing
    frame = 0
    while True:
        plot = Verify_Center_Plot(fnlist[frame],datalist[frame],flagged[frame],xcenlist[frame],ycenlist[frame],objx[frame],objy[frame],ghox[frame],ghoy[frame])
        if plot.key == "right": frame+=1
        if plot.key == "left": frame+=-1
        if plot.key == "up": flagged[frame]=True
        if plot.key == "down": flagged[frame]=False
        if frame == -1:
            while True:
                yn = raw_input("Finished reviewing images? (y/n) ")
                if "n" in yn or "N" in yn:
                    frame = 0
                    break
                elif "y" in yn or "Y" in yn:
                    break
        elif frame == len(fnlist):
            while True:
                yn = raw_input("Finished reviewing images? (y/n) ")
                if "n" in yn or "N" in yn:
                    frame = len(fnlist)-1
                    break
                elif "y" in yn or "Y" in yn:
                    break
        if frame == -1 or frame == len(fnlist) or plot.key == "q": break

    #Close all images
    for i in range(len(fnlist)):
        imagelist[i].close()
    
    #If any images are flagged, return False
    if np.sum(flagged)>0:
        print repr(np.int(np.sum(flagged)))+" centers rejected!"
        return False
    else:
        print "All centers approved!"
        return True
    
def find_ghost_centers(fnlist):
    """This routine finds the most likely optical center for a series of images
    by attempting to match ghost reflections of stars with their stellar
    counterparts. It can be rather slow if the fields are very dense with stars,
    but is also much more likely to succeed in locating the correct center in
    dense fields.
    
    Images should have already been aperture masked at this point. If they
    haven't, run 'aperture_mask' on them first.
    
    It is useful if the images have had their seeing FWHMs measured. If they
    have, these values (in pixels) should be in the fits headers as "FPFWHM".
    If not, the FWHM is assumed to be 5 pixels (and this is not ideal).
    
    The routine creates a number of temporary files while running, which are
    deleted at the end of the routine. Interrupting the routine while running
    is probably not a great idea.
    
    The routine returns a list of image centers as well as appending these to
    the fits headers as "fpxcen" and "fpycen"
    
    Inputs:
    fnlist -> List of strings, each one containing the path to a fits image.
              The images should be roughly aligned to one another or the
              routine will not work.
    
    Outputs:
    xcenlist -> List of image center X coordinates
    ycenlist -> List of image center Y coordinates
    
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
    
    #Identify objects in the images
    coolist = identify_objects(fnlist,skysig,fwhm)
    
    #Determine the likely location of the optical center for the whole track
    likelyxcen, likelyycen = likely_center(coolist,axcen,aycen,arad)
    
    #Now use likely center as a guess to find the best center for each image
    xcenlist, ycenlist = fit_better_centers(likelyxcen, likelyycen, coolist, axcen, aycen, arad)
    
    #Append the values to the image headers
    for i in range(len(fnlist)):
        image = openfits(fnlist[i],mode="update")
        image[0].header["fpxcen"] = xcenlist[i]
        image[0].header["fpycen"] = ycenlist[i]
        image.close()
    
    #Manually verify ghost centers
    while True:
        yn = raw_input("Manually verify ghost centers? (Recommended) (y/n) ")
        if "n" in yn or "N" in yn:
            break
        elif "y" in yn or "Y" in yn:
            goodcenters = verify_center(fnlist,coolist)
            if goodcenters: break
            else:
                print "Centers not approved!"
                clean_files(fnlist)
                crash()
    
    #Delete the coordinate files
    clean_files(fnlist)
    
    return xcenlist, ycenlist
