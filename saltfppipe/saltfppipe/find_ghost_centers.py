from astropy.io.fits import open as openfits
from sys import exit as crash
import numpy as np
from saltfppipe.fit_sky_level import fit_sky_level
import matplotlib.pyplot as plt
from photutils import daofind

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
                  "Use 'a'/'d' Keys to Navigate \n"+
                  "Use 'w'/'s' Keys to Flag Bad Center")
        for i in range(len(objx)):
            plt.plot([objx[i],ghox[i]],[objy[i],ghoy[i]],color="blue")
        plt.scatter([xcen],[ycen],color="red",s=30)
        plt.xlim(0,data.shape[1])
        plt.ylim(0,data.shape[0])
        self.cid = self.plot.figure.canvas.mpl_connect("key_press_event",self.keypress)
        plt.tight_layout()
        plt.show()
        
    def keypress(self,event):
        if event.key in ["w","a","s","d","q"]:
            self.key = event.key
            plt.close()

def verify_center(fnlist, objxs, objys, ghoxs, ghoys, xcens, ycens):
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
    
    #Change pyplot shortcut keys
    oldsavekey = plt.rcParams["keymap.save"]
    plt.rcParams["keymap.save"] = ""
    
    #Make empty lists to load things for plotting
    imagelist = []
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
        datalist.append(imagelist[i][1].data[aycenlist[i]-aradlist[i]:aycenlist[i]+aradlist[i],axcenlist[i]-aradlist[i]:axcenlist[i]+aradlist[i]])
        
    #Adjust object and ghost coordinates because of the aperture shift
    for i in range(len(fnlist)):
        objxs[i] = np.array(objxs[i])-axcenlist[i]+aradlist[i]
        objys[i] = np.array(objys[i])-aycenlist[i]+aradlist[i]
        ghoxs[i] = np.array(ghoxs[i])-axcenlist[i]+aradlist[i]
        ghoys[i] = np.array(ghoys[i])-aycenlist[i]+aradlist[i]
        xcens[i] = xcens[i]-axcenlist[i]+aradlist[i]
        ycens[i] = ycens[i]-aycenlist[i]+aradlist[i]
            
    #Interactive plotting routine thing
    i = 0
    while True:
        plot = Verify_Center_Plot(fnlist[i],datalist[i],flagged[i],xcens[i],ycens[i],objxs[i],objys[i],ghoxs[i],ghoys[i])
        if plot.key == "d": i+=1
        if plot.key == "a": i+=-1
        if plot.key == "w": flagged[i]=True
        if plot.key == "s": flagged[i]=False
        if i == -1:
            while True:
                yn = raw_input("Finished reviewing images? (y/n) ")
                if "n" in yn or "N" in yn:
                    i = 0
                    break
                elif "y" in yn or "Y" in yn:
                    break
        elif i == len(fnlist):
            while True:
                yn = raw_input("Finished reviewing images? (y/n) ")
                if "n" in yn or "N" in yn:
                    i = len(fnlist)-1
                    break
                elif "y" in yn or "Y" in yn:
                    break
        if i == -1 or i == len(fnlist) or plot.key == "q": break
    
    #Close all images
    for i in range(len(fnlist)):
        imagelist[i].close()

    #Restore old save key    
    plt.rcParams["keymap.save"] = oldsavekey
    
    #If any images are flagged, return False
    if np.sum(flagged)>0:
        print repr(np.int(np.sum(flagged)))+" centers rejected!"
        return False
    else:
        print "All centers approved!"
        return True
    
def find_ghost_centers(fnlist, tolerance=3, thresh=4.5):
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
    tolerance -> Optional. How close two pixels can be and be "close enough"
                 Default is 3 pixels.
    thresh -> Optional. Level above sky background variation to look for objs.
              Default is 3.5 (times SkySigma). Decrease if center positions
              aren't being found accurately. Increase for crowded fields to
              decrease computation time.
    
    Outputs:
    xcenlist -> List of image center X coordinates
    ycenlist -> List of image center Y coordinates
    
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
    maglists = []
    for i in range(len(fnlist)):
        xlists.append([])
        ylists.append([])
        maglists.append([])
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
                maglists[i].append(sources[j][-1])
        image.close()
    
    #Use the identified stars to come up with a center guess for each image
    xcen = np.zeros(len(fnlist))
    ycen = np.zeros(len(fnlist))
    goodcen = np.zeros(len(fnlist))
    for i in range(len(fnlist)):
        N = len(xlists[i])
        xcenarray = np.zeros(((N*N-N)/2))
        ycenarray = np.zeros(((N*N-N)/2))
        dist2array = np.zeros(((N*N-N)/2))
        sub1array = np.zeros(((N*N-N)/2))
        index = 0
        #All "possible" centers for all possible pairs of stars
        for j in range(N):
            for k in range(j+1,N):
                xcenarray[index] = xlists[i][j]+xlists[i][k]
                ycenarray[index] = ylists[i][j]+ylists[i][k]
                index = index+1
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
         
    #Cross-check the various image's centers against each other
    for i in range(len(fnlist)):
        dist2array = (xcen-xcen[i])**2 + (ycen-ycen[i])**2
        goodcen[i] = np.sum(dist2array<tolerance**2)
    
    #Determine where in the arrays the best fitting centers are
    bestcenloc = np.where(goodcen==max(goodcen))[0]
    bestxcen, bestycen = np.average(xcen[bestcenloc]), np.average(ycen[bestcenloc])
    #bestxcen, bestycen = 791, 503
    
    #Now we want to improve the center for each image using the best guesses
    xcenlist = np.zeros(len(fnlist))
    ycenlist = np.zeros(len(fnlist))
    objxs = []
    objys = []
    ghoxs = []
    ghoys = []
    for i in range(len(fnlist)):
        #Where would the reflected objects be if they exist
        refxlist = 2.*bestxcen - np.array(xlists[i])
        refylist = 2.*bestycen - np.array(ylists[i])
        #Populate lists of objects and ghosts based on this
        objxs.append([])
        objys.append([])
        ghoxs.append([])
        ghoys.append([])
        for j in range(len(xlists[i])):
            dist2list = (xlists[i] - refxlist[j])**2 + (ylists[i] - refylist[j])**2
            matchlist = dist2list < tolerance**2
            if np.sum(matchlist) >= 1:
                #We found a match! Now we need to know where the match(es) is(are)
                matchloc = np.where(matchlist==1)[0][0]
                #Cool, now, is this thing brighter than the ghost?
                if maglists[i][matchloc] < maglists[i][j]:
                    #It is! Record it:
                    objxs[i].append(xlists[i][matchloc])
                    objys[i].append(ylists[i][matchloc])
                    ghoxs[i].append(xlists[i][j])
                    ghoys[i].append(ylists[i][j])
        #Calculate the centers based on the object / ghost coords
        if len(objxs[i]) == 0:
            xcenlist[i] = 0
            ycenlist[i] = 0
        else:
            xcenlist[i] = 0.5*(np.average(objxs[i])+np.average(ghoxs[i]))
            ycenlist[i] = 0.5*(np.average(objys[i])+np.average(ghoys[i]))
    
    #Fill in the blanks with a "best guess"
    xcenlist[xcenlist==0] = np.average(xcenlist[xcenlist!=0])
    ycenlist[ycenlist==0] = np.average(ycenlist[ycenlist!=0])
    xcenlist[np.isnan(xcenlist)] = bestxcen
    ycenlist[np.isnan(ycenlist)] = bestycen
    
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
            goodcenters = verify_center(fnlist,objxs,objys,ghoxs,ghoys,xcenlist,ycenlist)
            if goodcenters: break
            else:
                print "Centers not approved!"
                crash()
        
    return xcenlist, ycenlist
